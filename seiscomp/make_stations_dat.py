#!/usr/bin/env python3
"""
make_stations_dat.py - Generate a Hyposat stations.dat from SeisComP inventory

Reads station coordinates from a SeisComP XML inventory file (scxmldump output)
and writes them in the Hyposat stations.dat format:

    SSSSS DDMMSS.sN/S DDDMMSS.sE/W  ELEV Description

Usage:
    # Dump SeisComP inventory to XML
    scxmldump -fI -o inventory.xml

    # Convert to Hyposat stations.dat
    python3 make_stations_dat.py inventory.xml > stations.dat

    # Or append to the existing NORSAR stations.dat (to add your local network)
    python3 make_stations_dat.py inventory.xml >> /path/to/hyposat/data/stations.dat

The generated file can be used directly as the STATION FILE in hyposat-parameter,
or merged with the NORSAR-provided stations.dat that comes with Hyposat.
"""

import sys
import xml.etree.ElementTree as ET
import argparse


def detect_namespace(root):
    """Extract namespace URI from the root element tag.

    Works with any SeisComP schema version (0.9 through 0.14 and beyond)
    by reading the namespace directly from the parsed XML rather than
    matching against a hardcoded list.
    """
    tag = root.tag
    if tag.startswith('{'):
        return tag[1:tag.index('}')]
    return ''


def q(ns, tag):
    return f'{{{ns}}}{tag}' if ns else tag


def decimal_to_dms(deg):
    """Convert decimal degrees to (degrees, minutes, seconds, hemisphere_sign)."""
    sign = 1 if deg >= 0 else -1
    deg  = abs(deg)
    d    = int(deg)
    m    = int((deg - d) * 60)
    s    = (deg - d - m / 60) * 3600
    return d, m, s, sign


def format_lat(lat_deg):
    """Format latitude as DDMMss.sN/S  (e.g.  694431.2N)"""
    d, m, s, sign = decimal_to_dms(lat_deg)
    hemi = 'N' if sign >= 0 else 'S'
    return f'{d:2d}{m:02d}{s:04.1f}{hemi}'


def format_lon(lon_deg):
    """Format longitude as DDDMMss.sE/W  (e.g. 0614918.0E)"""
    d, m, s, sign = decimal_to_dms(lon_deg)
    hemi = 'E' if sign >= 0 else 'W'
    return f'{d:3d}{m:02d}{s:04.1f}{hemi}'


def parse_inventory(xml_file):
    """
    Parse a SeisComP inventory XML file.
    Returns a list of dicts: {code, lat, lon, elev_m, description}
    One entry per unique station code (last epoch wins).
    """
    tree = ET.parse(xml_file)
    root = tree.getroot()
    ns   = detect_namespace(root)

    stations = {}   # code -> dict

    for net_el in root.iter(q(ns, 'network')):
        net_code = net_el.get('code', '')
        net_desc = net_el.get('description', '')

        for sta_el in net_el.findall(q(ns, 'station')):
            code = sta_el.get('code', '')
            if not code:
                continue

            # Coordinates are child elements, not attributes
            def _child_float(parent, tag, default=0.0):
                el = parent.find(q(ns, tag))
                try:
                    return float(el.text) if el is not None and el.text else default
                except ValueError:
                    return default

            lat  = _child_float(sta_el, 'latitude')
            lon  = _child_float(sta_el, 'longitude')
            elev = _child_float(sta_el, 'elevation')   # metres

            if lat == 0.0 and lon == 0.0:
                continue   # skip stations with no coordinates

            # Description is also a child element
            desc_el = sta_el.find(q(ns, 'description'))
            desc = (desc_el.text if desc_el is not None and desc_el.text else '') \
                   or net_desc or net_code
            # Truncate description to 50 chars (keep the file readable)
            desc = desc[:50]

            stations[code] = {
                'code':  code,
                'lat':   lat,
                'lon':   lon,
                'elev':  elev,   # metres (NEIC format uses metres)
                'desc':  desc,
            }

    return list(stations.values())


def write_stations_dat(stations, out_file=None):
    """
    Write stations in Hyposat stations.dat NEIC format.

    Fortran FORMAT (hyposat_geotab.f line 2154):
        format(A5, a1, i2,i2,f4.1,a1, i3,i2,f4.1,a1, f7.1, a48)
        SSSSS<sep>DDMMSS.sN/SDDDMMSS.sE/W<elev_m><desc>

    Example from NORSAR stations.dat:
        WRA   195633.4S1342021.8E  419.0Warramunga Array Beam Reference Point
    """
    out = open(out_file, 'w') if out_file else sys.stdout

    for sta in sorted(stations, key=lambda x: x['code']):
        code  = f"{sta['code']:<5s}"
        lat_s = format_lat(sta['lat'])
        lon_s = format_lon(sta['lon'])
        elev  = sta['elev']   # metres
        desc  = sta['desc']
        # No space between lat_s and lon_s — position 25 must be E/W for format detection
        line = f"{code} {lat_s}{lon_s}{elev:7.1f}{desc}\n"
        out.write(line)

    if out_file:
        out.close()


def main():
    parser = argparse.ArgumentParser(
        description='Convert SeisComP inventory XML to Hyposat stations.dat')
    parser.add_argument('inventory_xml',
                        help='SeisComP inventory XML file (from scxmldump -fI)')
    parser.add_argument('-o', '--output', default=None,
                        help='Output file (default: stdout)')
    args = parser.parse_args()

    try:
        stations = parse_inventory(args.inventory_xml)
    except ET.ParseError as e:
        print(f'ERROR: Could not parse XML: {e}', file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print(f'ERROR: File not found: {args.inventory_xml}', file=sys.stderr)
        sys.exit(1)

    if not stations:
        print('WARNING: No stations found in inventory', file=sys.stderr)
        sys.exit(0)

    print(f'INFO: Writing {len(stations)} stations', file=sys.stderr)
    write_stations_dat(stations, args.output)


if __name__ == '__main__':
    main()
