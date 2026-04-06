#!/usr/bin/env python3
"""
hyposat_wrapper.py - SeisComP ExternalLocator wrapper for Hyposat

Reads SeisComP XML (EventParameters with Origin + Picks) from stdin,
converts to Hyposat input format, runs the hyposat binary, parses
the output, and writes an updated SeisComP XML Origin to stdout.

Configure in SeisComP (e.g. global.cfg or scolv.cfg):
    plugins = locext
    ExternalLocator.profiles = Hyposat:@CONFIGDIR@/scripts/hyposat_wrapper.py

The wrapper accepts these command-line flags from ExternalLocator:
    --fixed-depth=<km>           Fix depth to this value
    --max-dist=<deg>             Ignore arrivals beyond this distance
    --ignore-initial-location    Let Hyposat determine the starting location

Environment variables (override the defaults below):
    HYPOSAT_BIN           Path to the hyposat executable
    HYPOSAT_DATA_DIR      Directory containing velocity models + stations.dat
    HYPOSAT_STATION_FILE  Path to Hyposat stations.dat  (default: DATA_DIR/stations.dat)
    HYPOSAT_MODEL         Velocity model name            (default: ak135_A)

IMPORTANT - Station coordinates:
    Hyposat reads station coordinates from its own stations.dat file.
    This file must contain every station that appears in your SeisComP picks.
    See seiscomp/make_stations_dat.py for a helper that generates stations.dat
    from the SeisComP inventory.
"""

import sys
import os
import re
import tempfile
import subprocess
import shutil
import argparse
import uuid
import xml.etree.ElementTree as ET
from datetime import datetime, timezone

# ---------------------------------------------------------------------------
# Configuration  (override with environment variables)
# ---------------------------------------------------------------------------

HYPOSAT_BIN = os.environ.get(
    'HYPOSAT_BIN',
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                 'bin', 'hyposat')
)

HYPOSAT_DATA_DIR = os.environ.get(
    'HYPOSAT_DATA_DIR',
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                 'data')
)

HYPOSAT_STATION_FILE = os.environ.get(
    'HYPOSAT_STATION_FILE',
    os.path.join(HYPOSAT_DATA_DIR, 'stations.dat')
)

HYPOSAT_MODEL = os.environ.get('HYPOSAT_MODEL', 'ak135_A')

HYPOMOD_BIN = os.environ.get(
    'HYPOMOD_BIN',
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                 'bin', 'hypomod')
)

# ---------------------------------------------------------------------------
# SeisComP XML helpers
# ---------------------------------------------------------------------------

def detect_namespace(xml_bytes):
    """Return the namespace URI from the root element, or empty string."""
    root = ET.fromstring(xml_bytes)
    m = re.match(r'\{(.+?)\}', root.tag)
    return m.group(1) if m else ''


def q(ns_uri, tag):
    """Qualify a tag name with its namespace URI."""
    return f'{{{ns_uri}}}{tag}' if ns_uri else tag


def get_text(el, ns_uri, *tags):
    """Follow a chain of child tags and return the text of the last one."""
    curr = el
    for tag in tags:
        if curr is None:
            return None
        curr = curr.find(q(ns_uri, tag))
    return curr.text if curr is not None else None


def get_float(el, ns_uri, *tags):
    t = get_text(el, ns_uri, *tags)
    try:
        return float(t) if t is not None else None
    except ValueError:
        return None

# ---------------------------------------------------------------------------
# Hyposat input formatting
# ---------------------------------------------------------------------------

def _fmt_sec(sec):
    """Format seconds as f6.3 (6 chars).  e.g. 44.7 -> '44.700'"""
    return f'{sec:6.3f}'


def _fmt_baz(baz):
    """Format backazimuth as f6.2 (6 chars).  -999 sentinel -> '-999. '"""
    if baz < -900.0:
        return '-999. '          # 6 chars, Fortran f6.2 reads as -999.0
    return f'{baz:6.2f}'


def _fmt_slow(slow):
    """Format slowness as f5.2 (5 chars).  -999 sentinel -> '-999.'"""
    if slow < -900.0:
        return '-999.'            # 5 chars, Fortran f5.2 reads as -999.0
    return f'{slow:5.2f}'


def format_hyposat_in(origin_el, picks_dict, ns_uri,
                      max_dist=None, ignore_init_loc=False):
    """
    Build hyposat-in file content.

    Fortran FORMAT used when line length >= 69 (hyposat.f ~line 2341):
        a5, 1x, a8, 1x, i4, 4(1x,i2), 1x, f6.3, 1x, f5.3,
        1x, f6.2, 3(1x,f5.2), 1x, a7

    Column positions (1-based):
        1-5   station (5 chars, left-justified)
        6     space
        7-14  phase   (8 chars, left-justified)
        15    space
        16-19 year
        21-22 month
        24-25 day
        27-28 hour
        30-31 minute
        33-38 second      f6.3
        40-44 time error  f5.3
        46-51 backazimuth f6.2  (-999. if unknown)
        53-57 baz error   f5.2
        59-63 slowness    f5.2  (-999. if unknown)
        65-69 slow error  f5.2
        71-77 flags       a7    (T=time, A=baz, S=slow, D=tt-diff, M=mag)
    """
    orig_time_str = get_text(origin_el, ns_uri, 'time', 'value') or ''
    title = f'SeisComP event {orig_time_str[:19]}'
    lines = [title]

    for arr_el in origin_el.findall(q(ns_uri, 'arrival')):
        pick_id = get_text(arr_el, ns_uri, 'pickID')
        phase   = (get_text(arr_el, ns_uri, 'phase') or 'P')[:8]
        weight  = get_float(arr_el, ns_uri, 'weight')

        if weight is not None and weight == 0.0:
            continue                      # excluded arrival

        pick_el = picks_dict.get(pick_id)
        if pick_el is None:
            continue

        wf = pick_el.find(q(ns_uri, 'waveformID'))
        if wf is None:
            continue
        station = wf.get('stationCode', '')[:5]

        # Optional distance cutoff
        if max_dist is not None:
            dist = get_float(arr_el, ns_uri, 'distance')
            if dist is not None and dist > max_dist:
                continue

        # Pick time
        time_val = get_text(pick_el, ns_uri, 'time', 'value')
        if not time_val:
            continue
        try:
            if '.' in time_val:
                dt = datetime.strptime(time_val, '%Y-%m-%dT%H:%M:%S.%fZ')
            else:
                dt = datetime.strptime(time_val, '%Y-%m-%dT%H:%M:%SZ')
        except ValueError:
            continue

        sec      = dt.second + dt.microsecond / 1e6
        time_err = get_float(pick_el, ns_uri, 'time', 'uncertainty') or 0.5

        # Backazimuth
        baz_used = (get_text(arr_el, ns_uri, 'backazimuthUsed') == 'true')
        baz      = get_float(pick_el, ns_uri, 'backazimuth', 'value')
        baz_err  = get_float(pick_el, ns_uri, 'backazimuth', 'uncertainty')

        if not baz_used or baz is None:
            baz, baz_err = -999.0, 0.0
        else:
            baz_err = baz_err if baz_err is not None else 20.0

        # Horizontal slowness (s/deg)
        slow_used = (get_text(arr_el, ns_uri, 'horizontalSlownessUsed') == 'true')
        slow      = get_float(pick_el, ns_uri, 'horizontalSlowness', 'value')
        slow_err  = get_float(pick_el, ns_uri, 'horizontalSlowness', 'uncertainty')

        if not slow_used or slow is None:
            slow, slow_err = -999.0, 0.0
        else:
            slow_err = slow_err if slow_err is not None else 2.0

        # Usage flags (7 chars)
        time_used_str = get_text(arr_el, ns_uri, 'timeUsed')
        time_used     = (time_used_str != 'false')   # True by default

        flag_t = 'T' if time_used else '_'
        flag_a = 'A' if (baz_used  and baz  > -900.0) else '_'
        flag_s = 'S' if (slow_used and slow > -900.0) else '_'
        flags  = f'{flag_t}{flag_a}{flag_s}____'     # 7 chars

        line = (
            f'{station:<5s} {phase:<8s} '
            f'{dt.year:4d} {dt.month:02d} {dt.day:02d} '
            f'{dt.hour:02d} {dt.minute:02d} {_fmt_sec(sec)} '
            f'{time_err:5.3f} {_fmt_baz(baz)} {baz_err:5.2f} '
            f'{_fmt_slow(slow)} {slow_err:5.2f} {flags}'
        )
        lines.append(line)

    return '\n'.join(lines) + '\n'


def format_hyposat_parameter(init_lat, init_lon, init_depth,
                              station_file, data_dir, model,
                              fixed_depth=None, ignore_init_loc=False):
    """Build hyposat-parameter file content."""

    if ignore_init_loc or init_lat is None or init_lon is None:
        lat_str = '999.'
        lon_str = '999.'
    else:
        lat_str = f'{init_lat:.4f}'
        lon_str = f'{init_lon:.4f}'

    if fixed_depth is not None:
        depth_str  = f'{fixed_depth:.2f}'
        depth_flag = 'f'
    else:
        depth_str  = f'{init_depth:.2f}' if init_depth is not None else '0.'
        depth_flag = 'b'

    return (
        f'GLOBAL MODEL                       : {model}\n'
        f'\n'
        f'CRUST 1.0                          : 4\n'
        f'OUTPUT OF REGIONAL MODEL           : 1\n'
        f'\n'
        f'STATION FILE                       : {station_file}\n'
        f'\n'
        f'P-VELOCITY TO CORRECT ELEVATION    : 4.5\n'
        f'S-VELOCITY TO CORRECT ELEVATION    : 0.\n'
        f'\n'
        f'LG GROUP-VELOCITY                  : 3.5752\n'
        f'RG GROUP-VELOCITY                  : 2.5\n'
        f'LQ GROUP-VELOCITY                  : 4.4\n'
        f'LR GROUP-VELOCITY                  : 2.85\n'
        f'\n'
        f'STARTING SOURCE LATITUDE           : {lat_str}\n'
        f'STARTING LATITUDE UNCERTAINTY      : 50.\n'
        f'\n'
        f'STARTING SOURCE LONGITUDE          : {lon_str}\n'
        f'STARTING LONGITUDE UNCERTAINTY     : 50.\n'
        f'\n'
        f'STARTING SOURCE DEPTH              : {depth_str}\n'
        f'STARTING DEPTH UNCERTAINTY         : 50.\n'
        f'DEPTH FLAG                         : {depth_flag}\n'
        f'\n'
        f'STARTING SOURCE TIME               : 0.\n'
        f'STARTING TIME UNCERTAINTY          : 600.\n'
        f'\n'
        f'MAXIMUM # OF ITERATIONS            : 80\n'
        f'ITERATIONS TO SEARCH OSCILLATIONS  : 6\n'
        f'\n'
        f'LOCATION ACCURACY                  : 1.\n'
        f'CONSTRAIN SOLUTION                 : 1\n'
        f'\n'
        f'CONFIDENCE LEVEL                   : 68.27\n'
        f'EPICENTER UNCERTAINTY ELLIPSE      : 1\n'
        f'\n'
        f'SLOWNESS [S/DEG]                   : 1\n'
        f'\n'
        f'MAXIMUM BAZ RESIDUAL               : 30.\n'
        f'MAXIMUM SLOWNESS RESIDUAL          : 5.\n'
        f'\n'
        f'FLAG USING TRAVEL-TIME DIFFERENCES : 1\n'
        f'\n'
        f'TT DATA UNCERTAINTY OUT            : 1\n'
        f'ISF_i                              : 0.05\n'
        f'ISF_o                              : 0.5\n'
        f'\n'
        f'MAGNITUDE CALCULATION              : 0\n'
        f'\n'
        f'INPUT FILE NAME                    : _\n'
        f'\n'
        f'FLAG EMERGENCE ANGLE OUTPUT        : 0\n'
        f'\n'
        f'OUTPUT SWITCH                      : 1\n'
        f'OUTPUT FILE NAME                   : _\n'
        f'OUTPUT LEVEL                       : 4\n'
    )

# ---------------------------------------------------------------------------
# HYPOMOD parameter file
# ---------------------------------------------------------------------------

def format_hypomod_parameter(loc, station_file, model):
    """Build a hyposat-parameter file for HYPOMOD.

    HYPOMOD uses the same parameter file as HYPOSAT but with the final
    solution pinned as the starting point (no inversion — it just evaluates
    residuals).  The source time is passed as an ISO-like epochal string;
    HYPOMOD accepts the human-readable form YYYY-MM-DD:HH.MM.SS.sss.
    """
    # Convert ISO time string to HYPOSAT human format: YYYY-MM-DD:HH.MM.SS.sss
    # loc['time_iso'] is like "2026-04-05T09:40:08.374Z"
    t = loc['time_iso'].rstrip('Z').replace('T', ' ')
    # YYYY-MM-DD HH:MM:SS.sss  →  YYYY-MM-DD:HH.MM.SS.sss
    date_part, time_part = t.split(' ')
    time_part = time_part.replace(':', '.')
    hypomod_time = f'{date_part}:{time_part}'

    return (
        f'GLOBAL MODEL                       : {model}\n'
        f'\n'
        f'CRUST 1.0                          : 4\n'
        f'OUTPUT OF REGIONAL MODEL           : 1\n'
        f'\n'
        f'STATION FILE                       : {station_file}\n'
        f'\n'
        f'P-VELOCITY TO CORRECT ELEVATION    : 4.5\n'
        f'S-VELOCITY TO CORRECT ELEVATION    : 0.\n'
        f'\n'
        f'LG GROUP-VELOCITY                  : 3.5752\n'
        f'RG GROUP-VELOCITY                  : 2.5\n'
        f'LQ GROUP-VELOCITY                  : 4.4\n'
        f'LR GROUP-VELOCITY                  : 2.85\n'
        f'\n'
        f'STARTING SOURCE LATITUDE           : {loc["lat"]:.4f}\n'
        f'STARTING LATITUDE UNCERTAINTY      : 0.\n'
        f'\n'
        f'STARTING SOURCE LONGITUDE          : {loc["lon"]:.4f}\n'
        f'STARTING LONGITUDE UNCERTAINTY     : 0.\n'
        f'\n'
        f'STARTING SOURCE DEPTH              : {loc["depth"]:.3f}\n'
        f'STARTING DEPTH UNCERTAINTY         : 0.\n'
        f'DEPTH FLAG                         : f\n'
        f'\n'
        f'STARTING SOURCE TIME               : {hypomod_time}\n'
        f'STARTING TIME UNCERTAINTY          : 0.\n'
        f'\n'
        f'MAXIMUM # OF ITERATIONS            : 0\n'
        f'\n'
        f'SLOWNESS [S/DEG]                   : 1\n'
        f'\n'
        f'MAXIMUM BAZ RESIDUAL               : 180.\n'
        f'MAXIMUM SLOWNESS RESIDUAL          : 999.\n'
        f'\n'
        f'FLAG USING TRAVEL-TIME DIFFERENCES : 1\n'
        f'\n'
        f'FLAG FREE PHASE SEARCH             : 1\n'
        f'\n'
        f'TT DATA UNCERTAINTY OUT            : 1\n'
        f'ISF_i                              : 0.05\n'
        f'ISF_o                              : 0.5\n'
        f'\n'
        f'MAGNITUDE CALCULATION              : 0\n'
        f'\n'
        f'INPUT FILE NAME                    : _\n'
        f'\n'
        f'OUTPUT SWITCH                      : 1\n'
        f'OUTPUT FILE NAME                   : _\n'
        f'OUTPUT LEVEL                       : 4\n'
    )


# ---------------------------------------------------------------------------
# Parse hyposat-out
# ---------------------------------------------------------------------------

def parse_hyposat_out(text):
    """
    Extract the final location solution from hyposat-out.

    The machine-readable summary line immediately follows the header:
        T0                         LAT      LON       Z  VPVS  DLAT   DLON    DZ     DT0  DVPVS DEF  RMS
        1996-06-29 00 36 49.153    1.322  126.294   44.78  1.76  0.0152  0.0360   4.43  0.496  0.04  54  0.654

    When depth is fixed, the DZ column contains the word "Fixed":
        1996-06-29 00 36 42.439    1.410  126.311    0.00  1.76  0.0164  0.0345  Fixed  0.089  0.04  45  0.503

    Returns a dict or None on failure.
    """
    lines = text.splitlines()

    # Find the last occurrence of the summary header
    header_idx = None
    for i, line in enumerate(lines):
        if ('T0' in line and 'LAT' in line and 'LON' in line
                and 'DEF' in line and 'RMS' in line):
            header_idx = i

    if header_idx is None:
        return None

    # The data line follows the header (skip any blanks)
    data_line = None
    for i in range(header_idx + 1, min(header_idx + 5, len(lines))):
        line = lines[i].strip()
        if line and re.match(r'\d{4}-\d{2}-\d{2}', line):
            data_line = line
            break

    if not data_line:
        return None

    # Tokenise:
    # date  HH  MM  SS.sss  LAT  LON  Z  VPVS  DLAT  DLON  DZ  DT0  DVPVS  DEF  RMS
    parts = data_line.split()
    if len(parts) < 13:
        return None

    try:
        date_part = parts[0]   # YYYY-MM-DD
        hh_part   = parts[1]
        mm_part   = parts[2]
        ss_part   = parts[3]
        lat_val   = float(parts[4])
        lon_val   = float(parts[5])
        depth_val = float(parts[6])
        # parts[7] = VPVS
        dlat_val  = float(parts[8])
        dlon_val  = float(parts[9])

        if parts[10].lower() == 'fixed':
            dz_val    = 0.0
            dt0_val   = float(parts[11])
            # parts[12] = DVPVS
            def_count = int(parts[13])   if len(parts) > 13 else 0
            rms_val   = float(parts[14]) if len(parts) > 14 else 0.0
        else:
            dz_val    = float(parts[10])
            dt0_val   = float(parts[11])
            # parts[12] = DVPVS
            def_count = int(parts[13])   if len(parts) > 13 else 0
            rms_val   = float(parts[14]) if len(parts) > 14 else 0.0

    except (ValueError, IndexError) as e:
        print(f'ERROR: parsing summary line: {e}', file=sys.stderr)
        return None

    # Parse origin time
    try:
        ss_f  = float(ss_part)
        ss_i  = int(ss_f)
        us    = int(round((ss_f - ss_i) * 1e6))
        dt    = datetime(
            int(date_part[0:4]), int(date_part[5:7]), int(date_part[8:10]),
            int(hh_part), int(mm_part), ss_i, us,
            tzinfo=timezone.utc
        )
        # Truncate microseconds to 6 digits (already done), format ISO
        time_iso = dt.strftime('%Y-%m-%dT%H:%M:%S.%f') + 'Z'
    except (ValueError, IndexError) as e:
        print(f'ERROR: parsing origin time: {e}', file=sys.stderr)
        return None

    result = {
        'lat':              lat_val,
        'lon':              lon_val,
        'depth':            depth_val,
        'time_iso':         time_iso,
        'dlat':             dlat_val,
        'dlon':             dlon_val,
        'dz':               dz_val,
        'dt0':              dt0_val,
        'def_count':        def_count,
        'rms':              rms_val,
        'major_axis':       None,
        'minor_axis':       None,
        'ellipse_azimuth':  None,
    }

    # Confidence ellipse axes (last occurrence)
    for m in re.finditer(
            r'Major half axis:\s*([\d.]+)\s*\[km\].*?'
            r'Minor half axis:\s*([\d.]+)\s*\[km\]', text):
        result['major_axis'] = float(m.group(1))
        result['minor_axis'] = float(m.group(2))

    for m in re.finditer(r'Azimuth:\s*([\d.]+)\s*\[deg\]', text):
        result['ellipse_azimuth'] = float(m.group(1))

    return result

# ---------------------------------------------------------------------------
# Build output SeisComP XML
# ---------------------------------------------------------------------------

def _copy_element(src, ns_uri, dest_parent):
    """Shallow-copy a parsed element (and its children) into dest_parent."""
    el = ET.SubElement(dest_parent, src.tag, src.attrib)
    el.text = src.text
    el.tail = None
    for child in src:
        _copy_element(child, ns_uri, el)
    return el


def build_output_xml(origin_el, picks_dict, loc, ns_uri):
    """
    Construct an EventParameters XML document with the Hyposat-derived origin.
    """
    ET.register_namespace('', ns_uri)       # suppress 'ns0:' prefix

    root = ET.Element(q(ns_uri, 'seiscomp'))

    # Origin must be a DIRECT child of <seiscomp>.
    # The XMLArchive's findTag() only searches direct children of _current
    # (set to the <seiscomp> node), so wrapping in <EventParameters> hides it.

    # New origin (new publicID so SeisComP treats it as an update)
    # Tag must be 'Origin' (capital O) — XMLArchive.isValidTag() compares
    # against T::ClassName() = "Origin" using case-sensitive xmlStrcmp.
    new_id  = f'Hyposat/{uuid.uuid4()}'
    orig_el = ET.SubElement(root, q(ns_uri, 'Origin'))
    orig_el.set('publicID', new_id)

    # --- Time ---
    time_el = ET.SubElement(orig_el, q(ns_uri, 'time'))
    ET.SubElement(time_el, q(ns_uri, 'value')).text = loc['time_iso']
    ET.SubElement(time_el, q(ns_uri, 'uncertainty')).text = f"{loc['dt0']:.4f}"

    # --- Latitude ---
    lat_el = ET.SubElement(orig_el, q(ns_uri, 'latitude'))
    ET.SubElement(lat_el, q(ns_uri, 'value')).text = f"{loc['lat']:.6f}"
    ET.SubElement(lat_el, q(ns_uri, 'uncertainty')).text = f"{loc['dlat']:.6f}"

    # --- Longitude ---
    lon_el = ET.SubElement(orig_el, q(ns_uri, 'longitude'))
    ET.SubElement(lon_el, q(ns_uri, 'value')).text = f"{loc['lon']:.6f}"
    ET.SubElement(lon_el, q(ns_uri, 'uncertainty')).text = f"{loc['dlon']:.6f}"

    # --- Depth (km) ---
    dep_el = ET.SubElement(orig_el, q(ns_uri, 'depth'))
    ET.SubElement(dep_el, q(ns_uri, 'value')).text = f"{loc['depth']:.3f}"
    if loc['dz'] > 0.0:
        ET.SubElement(dep_el, q(ns_uri, 'uncertainty')).text = f"{loc['dz']:.3f}"

    # --- Method / earth model ---
    ET.SubElement(orig_el, q(ns_uri, 'methodID')).text    = 'HYPOSAT'
    ET.SubElement(orig_el, q(ns_uri, 'earthModelID')).text = HYPOSAT_MODEL

    # --- Quality (compute stats from arrivals that will be written) ---
    used_arrivals = [
        arr_el for arr_el in origin_el.findall(q(ns_uri, 'arrival'))
        if not ((get_float(arr_el, ns_uri, 'weight') or 1.0) == 0.0)
    ]
    all_arrivals = origin_el.findall(q(ns_uri, 'arrival'))

    # Azimuthal gap from azimuth values on used arrivals
    azimuths: list[float] = sorted(
        v for a in used_arrivals
        if (v := get_float(a, ns_uri, 'azimuth')) is not None
    )
    if len(azimuths) >= 2:
        gaps = [azimuths[i+1] - azimuths[i] for i in range(len(azimuths)-1)]
        gaps.append(360.0 - azimuths[-1] + azimuths[0])
        az_gap = max(gaps)
    else:
        az_gap = None

    # Distance stats from used arrivals
    distances: list[float] = sorted(
        v for a in used_arrivals
        if (v := get_float(a, ns_uri, 'distance')) is not None
    )

    # Unique station codes used
    used_sta = set()
    for a in used_arrivals:
        pid = get_text(a, ns_uri, 'pickID')
        # station code encoded in arrival azimuth/distance — use pick lookup below
        # (we count arrivals as a proxy; exact dedup done via pick waveformID)
        if pid:
            used_sta.add(pid)

    qual_el = ET.SubElement(orig_el, q(ns_uri, 'quality'))
    ET.SubElement(qual_el, q(ns_uri, 'associatedPhaseCount')).text = str(len(all_arrivals))
    ET.SubElement(qual_el, q(ns_uri, 'usedPhaseCount')).text       = str(loc['def_count'])
    ET.SubElement(qual_el, q(ns_uri, 'standardError')).text        = f"{loc['rms']:.4f}"
    if az_gap is not None:
        ET.SubElement(qual_el, q(ns_uri, 'azimuthalGap')).text     = f"{az_gap:.2f}"
    if distances:
        ET.SubElement(qual_el, q(ns_uri, 'minimumDistance')).text  = f"{distances[0]:.4f}"
        ET.SubElement(qual_el, q(ns_uri, 'maximumDistance')).text  = f"{distances[-1]:.4f}"
        mid = distances[len(distances) // 2]
        ET.SubElement(qual_el, q(ns_uri, 'medianDistance')).text   = f"{mid:.4f}"

    # --- Uncertainty ellipse ---
    if loc['major_axis'] is not None:
        unc_el = ET.SubElement(orig_el, q(ns_uri, 'uncertainty'))
        ET.SubElement(unc_el, q(ns_uri, 'maxHorizontalUncertainty')).text = \
            f"{loc['major_axis']:.3f}"
        ET.SubElement(unc_el, q(ns_uri, 'minHorizontalUncertainty')).text = \
            f"{loc['minor_axis']:.3f}"
        if loc['ellipse_azimuth'] is not None:
            ET.SubElement(unc_el,
                          q(ns_uri, 'azimuthMaxHorizontalUncertainty')).text = \
                f"{loc['ellipse_azimuth']:.1f}"

    # --- Arrivals (copy non-zero-weight arrivals from input) ---
    for arr_el in origin_el.findall(q(ns_uri, 'arrival')):
        w = get_float(arr_el, ns_uri, 'weight')
        if w is not None and w == 0.0:
            continue
        _copy_element(arr_el, ns_uri, orig_el)

    # Pretty-print if Python >= 3.9
    if sys.version_info >= (3, 9):
        ET.indent(root, space='  ')

    xml_body = ET.tostring(root, encoding='unicode')
    return '<?xml version="1.0" encoding="UTF-8"?>\n' + xml_body

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='SeisComP ExternalLocator wrapper for Hyposat')
    parser.add_argument('--fixed-depth', type=float, default=None,
                        metavar='KM',
                        help='Fix hypocenter depth to this value [km]')
    parser.add_argument('--max-dist', type=float, default=None,
                        metavar='DEG',
                        help='Ignore arrivals beyond this distance [deg]')
    parser.add_argument('--ignore-initial-location', action='store_true',
                        help='Let Hyposat compute the starting location')
    parser.add_argument('--hypomod-log', metavar='PATH', default=None,
                        help='Run HYPOMOD after location and append residuals '
                             'to this log file')
    args = parser.parse_args()

    # --- Read XML from stdin ---
    xml_bytes = sys.stdin.buffer.read()

    ns_uri = detect_namespace(xml_bytes)

    root   = ET.fromstring(xml_bytes)
    ep     = root.find(q(ns_uri, 'EventParameters'))
    if ep is None:
        print('ERROR: No EventParameters element in input XML', file=sys.stderr)
        sys.exit(1)

    picks_dict = {}
    for pick_el in ep.findall(q(ns_uri, 'pick')):
        pid = pick_el.get('publicID')
        if pid:
            picks_dict[pid] = pick_el

    origin_el = ep.find(q(ns_uri, 'origin'))
    if origin_el is None:
        print('ERROR: No origin element in input XML', file=sys.stderr)
        sys.exit(1)

    init_lat   = get_float(origin_el, ns_uri, 'latitude',  'value')
    init_lon   = get_float(origin_el, ns_uri, 'longitude', 'value')
    init_depth = get_float(origin_el, ns_uri, 'depth',     'value')

    # --- Build Hyposat input ---
    hyposat_in_text = format_hyposat_in(
        origin_el, picks_dict, ns_uri,
        max_dist=args.max_dist,
        ignore_init_loc=args.ignore_initial_location
    )
    hyposat_param_text = format_hyposat_parameter(
        init_lat, init_lon, init_depth,
        station_file=HYPOSAT_STATION_FILE,
        data_dir=HYPOSAT_DATA_DIR,
        model=HYPOSAT_MODEL,
        fixed_depth=args.fixed_depth,
        ignore_init_loc=args.ignore_initial_location
    )

    # --- Run Hyposat in a temporary directory ---
    loc = None
    tmpdir = tempfile.mkdtemp(prefix='hyposat_sc_')
    try:
        in_file    = os.path.join(tmpdir, 'hyposat-in')
        param_file = os.path.join(tmpdir, 'hyposat-parameter')
        out_file   = os.path.join(tmpdir, 'hyposat-out')

        with open(in_file,    'w') as f:
            f.write(hyposat_in_text)
        with open(param_file, 'w') as f:
            f.write(hyposat_param_text)

        env = os.environ.copy()
        env['HYPOSAT_DATA'] = HYPOSAT_DATA_DIR

        n_arrivals = hyposat_in_text.count('\n') - 1  # minus title line
        print(f'INFO: sending {n_arrivals} arrivals to hyposat', file=sys.stderr)
        print(f'INFO: running {HYPOSAT_BIN} in {tmpdir}', file=sys.stderr)
        result = subprocess.run(
            [HYPOSAT_BIN],
            cwd=tmpdir,
            env=env,
            stdin=subprocess.DEVNULL,   # prevent inheriting scolv's stdin pipe
            capture_output=True,
            text=True,
            timeout=60
        )

        # Forward Hyposat's stdout (contains warnings, station-not-found, etc.)
        for line in result.stdout.splitlines():
            print(f'HYPOSAT: {line}', file=sys.stderr)

        if result.returncode != 0:
            print(f'ERROR: hyposat exited with code {result.returncode}',
                  file=sys.stderr)
            print(result.stderr, file=sys.stderr)
            sys.exit(1)

        if not os.path.exists(out_file):
            print('ERROR: hyposat-out was not created', file=sys.stderr)
            sys.exit(1)

        with open(out_file, 'r') as f:
            hyposat_out = f.read()

        # --- Parse Hyposat output (inside try so tmpdir is still available) ---
        loc = parse_hyposat_out(hyposat_out)
        if loc is None:
            print('ERROR: Could not parse hyposat-out summary line', file=sys.stderr)
            print(hyposat_out[-2000:], file=sys.stderr)
            sys.exit(1)

        print(f'INFO: solution  lat={loc["lat"]:.4f}  lon={loc["lon"]:.4f}  '
              f'z={loc["depth"]:.1f} km  rms={loc["rms"]:.3f} s  '
              f'def={loc["def_count"]}', file=sys.stderr)

        # --- Optionally run HYPOMOD ---
        if args.hypomod_log and os.path.isfile(HYPOMOD_BIN):
            hypomod_param_text = format_hypomod_parameter(
                loc,
                station_file=HYPOSAT_STATION_FILE,
                model=HYPOSAT_MODEL,
            )
            hypomod_param_file = os.path.join(tmpdir, 'hyposat-parameter')
            hypomod_out_file   = os.path.join(tmpdir, 'hypomod-out')

            # HYPOMOD reads from 'hypomod-in' (not 'hyposat-in')
            shutil.copy(in_file, os.path.join(tmpdir, 'hypomod-in'))

            # Overwrite the parameter file with the HYPOMOD version
            with open(hypomod_param_file, 'w') as f:
                f.write(hypomod_param_text)

            print(f'INFO: running hypomod for residual log', file=sys.stderr)
            hm_result = subprocess.run(
                [HYPOMOD_BIN],
                cwd=tmpdir,
                env=env,
                stdin=subprocess.DEVNULL,
                capture_output=True,
                text=True,
                timeout=30,
            )

            if hm_result.returncode == 0 and os.path.exists(hypomod_out_file):
                with open(hypomod_out_file, 'r') as f:
                    hypomod_out = f.read()

                # Append to the log file with an event header
                try:
                    os.makedirs(os.path.dirname(os.path.abspath(args.hypomod_log)),
                                exist_ok=True)
                    with open(args.hypomod_log, 'a') as f:
                        f.write(f'\n{"="*72}\n')
                        f.write(f'Event: {loc["time_iso"]}  '
                                f'lat={loc["lat"]:.4f}  lon={loc["lon"]:.4f}  '
                                f'z={loc["depth"]:.1f} km  '
                                f'rms={loc["rms"]:.3f} s  '
                                f'model={HYPOSAT_MODEL}\n')
                        f.write(f'{"="*72}\n')
                        f.write(hypomod_out)
                    print(f'INFO: hypomod residuals appended to {args.hypomod_log}',
                          file=sys.stderr)
                except OSError as e:
                    print(f'WARNING: could not write hypomod log: {e}',
                          file=sys.stderr)
            else:
                print(f'WARNING: hypomod failed (exit {hm_result.returncode})',
                      file=sys.stderr)
                if hm_result.stdout:
                    print(hm_result.stdout[:500], file=sys.stderr)

    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    if loc is None:
        print('ERROR: Hyposat solution was not obtained', file=sys.stderr)
        sys.exit(1)

    # --- Build and emit output XML ---
    out_xml = build_output_xml(origin_el, picks_dict, loc, ns_uri)
    sys.stdout.write(out_xml)
    sys.stdout.flush()


if __name__ == '__main__':
    main()
