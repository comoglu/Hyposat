#!/usr/bin/env python3
"""
check_dfx_baz.py — Verify DFX backazimuth measurements from SeisComP.

Accepts two input formats:

  1. SeisComP XML  (scxmldump -fPAO output)
     Station coordinates are read from stations.dat; the theoretical BAZ is
     computed from the great-circle formula.

  2. scolv arrival CSV export  (semicolon-delimited, header starts with
     "# Used;ID;Status;...")
     The 'Az' column already gives the forward azimuth from the epicentre to
     the station, so theoretical BAZ = Az + 180°.  No station coordinates
     needed.

Usage
-----
  # XML input (preferred — full detail including DFX:rectilinearity)
  python3 check_dfx_baz.py event.xml

  # scolv CSV export
  python3 check_dfx_baz.py arrivals.csv

  # Override the minimum-distance threshold (default 20°)
  python3 check_dfx_baz.py event.xml --min-dist 30

  # Use a specific origin rather than the preferred one (XML only)
  python3 check_dfx_baz.py event.xml --origin-id 20260405013508.279513.88656

  # Override the residual threshold for GOOD/BAD classification (default 30°)
  python3 check_dfx_baz.py event.xml --threshold 20

Output
------
  A table, one row per pick that has a DFX backazimuth:

    Station       dist°  obs_baz   ±unc  theo_baz   resid    rect_raw  status
    AU.WRKA        28.3   358.9      ?     356.0     +2.9         8.2  GOOD
    IM.AS31        27.9   346.6      ?     345.0     +1.6     33655.0  GOOD
    AU.QIS         27.1    21.3      ?     331.0    +50.3        33.2  BAD
    AU.LCRK         7.8   311.4      ?     124.2   -172.9       277.7  NEAR

  Status codes:
    GOOD      — residual within --threshold (default ±30°)
    BAD       — residual exceeds threshold
    NEAR      — station closer than --min-dist; DFX unreliable at steep incidence
    NO_COORDS — station not found in stations.dat  (XML mode only)

  Followed by a summary with RMS (stations beyond --min-dist only).

Notes on rectilinearity
-----------------------
  DFX stores 2 × λ₁ (twice the largest covariance eigenvalue) as the
  DFX:rectilinearity Pick comment.  This is NOT the normalised Jurkevics
  (1988) value (0–1).  Use it only as a relative ranking within the same
  event — higher values = stronger signal polarization.
"""

import argparse
import math
import os
import sys
import xml.etree.ElementTree as ET

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SC_NS = "http://geofon.gfz.de/ns/seiscomp-schema/0.14"


def q(tag: str) -> str:
    return f"{{{SC_NS}}}{tag}"


# ---------------------------------------------------------------------------
# Great-circle geometry
# ---------------------------------------------------------------------------

def _great_circle(lat1_deg, lon1_deg, lat2_deg, lon2_deg):
    """Return (distance_deg, forward_azimuth_deg, backazimuth_deg).

    Forward azimuth  = direction from point-1 toward point-2.
    Backazimuth      = direction from point-2 back toward point-1 —
                       what DFX measures at the station pointing at the source.
    """
    lat1, lon1 = math.radians(lat1_deg), math.radians(lon1_deg)
    lat2, lon2 = math.radians(lat2_deg), math.radians(lon2_deg)
    dlon = lon2 - lon1

    sin_dlon = math.sin(dlon)
    cos_dlon = math.cos(dlon)
    sin_lat1 = math.sin(lat1)
    cos_lat1 = math.cos(lat1)
    sin_lat2 = math.sin(lat2)
    cos_lat2 = math.cos(lat2)

    a = math.sqrt((cos_lat2 * sin_dlon) ** 2 +
                  (cos_lat1 * sin_lat2 - sin_lat1 * cos_lat2 * cos_dlon) ** 2)
    b = sin_lat1 * sin_lat2 + cos_lat1 * cos_lat2 * cos_dlon
    dist_deg = math.degrees(math.atan2(a, b))

    fwd_az = math.atan2(sin_dlon * cos_lat2,
                        cos_lat1 * sin_lat2 - sin_lat1 * cos_lat2 * cos_dlon)
    fwd_az_deg = (math.degrees(fwd_az) + 360.0) % 360.0

    baz = math.atan2(-sin_dlon * cos_lat1,
                     cos_lat2 * sin_lat1 - sin_lat2 * cos_lat1 * cos_dlon)
    baz_deg = (math.degrees(baz) + 360.0) % 360.0

    return dist_deg, fwd_az_deg, baz_deg


def _wrap180(angle):
    return ((angle + 180.0) % 360.0) - 180.0


# ---------------------------------------------------------------------------
# stations.dat reader  (NEIC format)
# ---------------------------------------------------------------------------

def load_stations(path: str) -> dict:
    """Return {STATIONCODE: (lat_deg, lon_deg)} from a NEIC stations.dat file.

    The ZZZZ sentinel line separates sections (NORSAR original from
    SeisComP-derived entries); we skip it and continue reading.
    Duplicates: first occurrence wins (same as Hyposat).
    """
    stations = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if len(line) < 25:
                continue
            code = line[0:5].strip()
            if not code or code.startswith("ZZZZ"):
                continue  # sentinel — skip, don't stop
            lat_s = line[6:15]   # DDMMss.sN/S
            lon_s = line[15:25]  # DDDMMss.sE/W
            if len(lat_s) < 9 or len(lon_s) < 10:
                continue
            hemi_lat = lat_s[8]
            hemi_lon = lon_s[9]
            if hemi_lat not in "NS" or hemi_lon not in "EW":
                continue
            try:
                lat = (int(lat_s[0:2]) + int(lat_s[2:4]) / 60.0
                       + float(lat_s[4:8]) / 3600.0)
                lon = (int(lon_s[0:3]) + int(lon_s[3:5]) / 60.0
                       + float(lon_s[5:9]) / 3600.0)
            except ValueError:
                continue
            if hemi_lat == "S":
                lat = -lat
            if hemi_lon == "W":
                lon = -lon
            if code not in stations:
                stations[code] = (lat, lon)
    return stations


# ---------------------------------------------------------------------------
# SeisComP XML helpers
# ---------------------------------------------------------------------------

def get_text(el, *tags):
    node = el
    for tag in tags:
        if node is None:
            return None
        node = node.find(q(tag))
    return node.text.strip() if node is not None and node.text else None


def get_float(el, *tags):
    t = get_text(el, *tags)
    try:
        return float(t) if t is not None else None
    except ValueError:
        return None


# ---------------------------------------------------------------------------
# Format detection
# ---------------------------------------------------------------------------

def _is_csv(path: str) -> bool:
    """Return True if the file looks like a scolv arrival CSV export."""
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                return line.startswith("# Used;") or line.startswith("Used;")
    except OSError:
        pass
    return False


# ---------------------------------------------------------------------------
# CSV input mode
# ---------------------------------------------------------------------------

# scolv CSV column indices (0-based, semicolon-delimited)
#  0  Used
#  1  ID
#  2  Status
#  3  Phase
#  4  Weight
#  5  Method
#  6  Polarity
#  7  Onset
#  8  Takeoff (°)
#  9  Net
# 10  Sta
# 11  Loc/Cha
# 12  Timeres (s)
# 13  Dis (°)
# 14  Az (°)         ← forward azimuth from epicentre to station
# 15  Time (UTC)
# 16  +/- (s)
# 17  Slo (s/°)
# 18  Slores (s/°)
# 19  Baz (°)        ← DFX-measured backazimuth (may be empty)
# 20  Bazres (°)     ← residual computed by LOCSAT (may be empty)
# 21  Created
# 22  Latency (s)

_COL_NET    = 9
_COL_STA    = 10
_COL_LOCHA  = 11
_COL_DIST   = 13
_COL_AZ     = 14
_COL_BAZ    = 19
_COL_BAZRES = 20


def _float_or_none(s: str):
    s = s.strip()
    try:
        return float(s) if s else None
    except ValueError:
        return None


def process_csv(path: str, min_dist: float, threshold: float):
    rows = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line or line.startswith("#") or line.lower().startswith("used"):
                continue
            parts = line.split(";")
            if len(parts) < 20:
                continue

            net  = parts[_COL_NET].strip()
            sta  = parts[_COL_STA].strip()
            dist = _float_or_none(parts[_COL_DIST])
            az   = _float_or_none(parts[_COL_AZ])
            obs  = _float_or_none(parts[_COL_BAZ])
            locsat_res = _float_or_none(parts[_COL_BAZRES])

            if obs is None:
                continue  # no DFX BAZ

            if az is not None:
                theo = (az + 180.0) % 360.0
                res  = _wrap180(obs - theo)
            else:
                theo = None
                res  = None

            if dist is None:
                status = "NO_DIST"
            elif dist < min_dist:
                status = "NEAR"
            elif res is None:
                status = "NO_AZ"
            elif abs(res) <= threshold:
                status = "GOOD"
            else:
                status = "BAD"

            rows.append({
                "net": net, "sta": sta,
                "dist": dist, "az": az,
                "obs": obs, "theo": theo,
                "res": res, "locsat_res": locsat_res,
                "status": status,
            })

    rows.sort(key=lambda r: (r["dist"] is None, r["dist"] or 0))

    hdr = (f"{'Station':<12} {'dist°':>6}  {'Az→sta':>6}  {'obs_baz':>7}  "
           f"{'theo_baz':>8}  {'resid':>6}  {'LOCSAT_res':>10}  status")
    print()
    print(hdr)
    print("-" * len(hdr))

    for r in rows:
        sta_label  = f"{r['net']}.{r['sta']}"
        dist_s     = f"{r['dist']:6.1f}" if r["dist"] is not None else "     ?"
        az_s       = f"{r['az']:6.1f}"   if r["az"]   is not None else "     ?"
        obs_s      = f"{r['obs']:7.1f}"
        theo_s     = f"{r['theo']:8.1f}" if r["theo"] is not None else "       ?"
        res_s      = f"{r['res']:+6.1f}" if r["res"]  is not None else "     ?"
        lcres_s    = f"{r['locsat_res']:+10.1f}" if r["locsat_res"] is not None else "         ?"
        print(f"{sta_label:<12} {dist_s}  {az_s}  {obs_s}  "
              f"{theo_s}  {res_s}  {lcres_s}  {r['status']}")

    _print_summary(rows, min_dist, threshold, mode="csv")


# ---------------------------------------------------------------------------
# XML input mode
# ---------------------------------------------------------------------------

def process_xml(path: str, stations_path: str, origin_id: str,
                min_dist: float, threshold: float):

    if not os.path.isfile(stations_path):
        sys.exit(f"ERROR: stations.dat not found: {stations_path}")
    stations = load_stations(stations_path)
    print(f"Loaded {len(stations)} stations from {stations_path}", file=sys.stderr)

    try:
        tree = ET.parse(path)
    except ET.ParseError as exc:
        sys.exit(f"ERROR: cannot parse XML: {exc}")

    root = tree.getroot()
    ep = root.find(q("EventParameters"))
    if ep is None:
        sys.exit("ERROR: no <EventParameters> in XML")

    picks_by_id = {}
    for pick in ep.iter(q("pick")):
        pid = pick.get("publicID")
        if pid:
            picks_by_id[pid] = pick

    origins_by_id = {}
    for origin in ep.iter(q("origin")):
        oid = origin.get("publicID")
        if oid:
            origins_by_id[oid] = origin

    preferred_origin = None
    if origin_id:
        for oid, orig in origins_by_id.items():
            if origin_id in oid:
                preferred_origin = orig
                print(f"Using origin: {oid}", file=sys.stderr)
                break
        if preferred_origin is None:
            sys.exit(f"ERROR: origin '{origin_id}' not found in XML")
    else:
        event = ep.find(q("event"))
        if event is not None:
            pref_id = get_text(event, "preferredOriginID")
            if pref_id and pref_id in origins_by_id:
                preferred_origin = origins_by_id[pref_id]
                print(f"Using preferred origin: {pref_id}", file=sys.stderr)
        if preferred_origin is None:
            if origins_by_id:
                preferred_origin = next(iter(origins_by_id.values()))
                print("WARNING: no preferredOriginID — using first origin",
                      file=sys.stderr)
            else:
                sys.exit("ERROR: no origins found in XML")

    epi_lat = get_float(preferred_origin, "latitude",  "value")
    epi_lon = get_float(preferred_origin, "longitude", "value")
    epi_dep = get_float(preferred_origin, "depth",     "value")
    if epi_lat is None or epi_lon is None:
        sys.exit("ERROR: could not read epicentre coordinates from origin")

    method = get_text(preferred_origin, "methodID")     or "?"
    model  = get_text(preferred_origin, "earthModelID") or "?"
    print(f"Epicentre: lat={epi_lat:.4f}°  lon={epi_lon:.4f}°  "
          f"depth={epi_dep:.1f} km  [{method}/{model}]", file=sys.stderr)

    # Build a map: pickID → arrival element (from the preferred origin).
    # This lets us check backazimuthUsed and backazimuthResidual per pick.
    arrivals_by_pick = {}
    for arr in preferred_origin.iter(q("arrival")):
        pid = get_text(arr, "pickID")
        if pid:
            arrivals_by_pick[pid] = arr

    rows = []
    for pid, pick in picks_by_id.items():
        baz_obs = get_float(pick, "backazimuth", "value")
        if baz_obs is None:
            continue

        wf = pick.find(q("waveformID"))
        if wf is None:
            continue
        net  = wf.get("networkCode", "?")
        sta  = wf.get("stationCode",  "?")

        rect = None
        for comment in pick.findall(q("comment")):
            if get_text(comment, "id") == "DFX:rectilinearity":
                rect = get_float(comment, "text")
                break

        baz_unc = get_float(pick, "backazimuth", "uncertainty")

        # Check whether the locator actually used this BAZ observation
        arr_el = arrivals_by_pick.get(pid)
        baz_used     = (get_text(arr_el, "backazimuthUsed") == "true") if arr_el is not None else None
        locsat_bazres = get_float(arr_el, "backazimuthResidual") if arr_el is not None else None

        if sta not in stations:
            rows.append({
                "net": net, "sta": sta,
                "dist": None, "obs": baz_obs, "theo": None,
                "res": None, "rect": rect, "unc": baz_unc,
                "baz_used": baz_used, "locsat_bazres": locsat_bazres,
                "status": "NO_COORDS",
            })
            continue

        slat, slon = stations[sta]
        dist_deg, theo_baz, _ = _great_circle(slat, slon, epi_lat, epi_lon)
        residual = _wrap180(baz_obs - theo_baz)

        if dist_deg < min_dist:
            status = "NEAR"
        elif abs(residual) <= threshold:
            status = "GOOD"
        else:
            status = "BAD"

        rows.append({
            "net": net, "sta": sta,
            "dist": dist_deg, "obs": baz_obs, "theo": theo_baz,
            "res": residual, "rect": rect, "unc": baz_unc,
            "baz_used": baz_used, "locsat_bazres": locsat_bazres,
            "status": status,
        })

    rows.sort(key=lambda r: (r["dist"] is None, r["dist"] or 0))

    hdr = (f"{'Station':<12} {'dist°':>6}  {'obs_baz':>7}  {'±unc':>5}  "
           f"{'theo_baz':>8}  {'resid':>6}  {'lc_res':>6}  {'used':>4}  "
           f"{'rect_raw':>12}  status")
    print()
    print(hdr)
    print("-" * len(hdr))

    for r in rows:
        sta_label  = f"{r['net']}.{r['sta']}"
        dist_s     = f"{r['dist']:6.1f}"  if r["dist"] is not None else "     ?"
        obs_s      = f"{r['obs']:7.1f}"
        unc_s      = f"{r['unc']:5.1f}"   if r["unc"]  is not None else "    ?"
        theo_s     = f"{r['theo']:8.1f}"  if r["theo"] is not None else "       ?"
        res_s      = f"{r['res']:+6.1f}"  if r["res"]  is not None else "     ?"
        rect_s     = f"{r['rect']:12.1f}" if r["rect"] is not None else "           ?"
        lcres_s    = (f"{r['locsat_bazres']:+6.1f}"
                      if r.get("locsat_bazres") is not None else "     ?")
        if r.get("baz_used") is True:
            used_s = " YES"
        elif r.get("baz_used") is False:
            used_s = "  no"
        else:
            used_s = "   -"
        print(f"{sta_label:<12} {dist_s}  {obs_s}  {unc_s}  "
              f"{theo_s}  {res_s}  {lcres_s}  {used_s}  {rect_s}  {r['status']}")

    _print_summary(rows, min_dist, threshold, mode="xml")


# ---------------------------------------------------------------------------
# Shared summary
# ---------------------------------------------------------------------------

def _print_summary(rows, min_dist, threshold, mode):
    valid     = [r for r in rows if r.get("dist") is not None
                 and r["dist"] >= min_dist and r["status"] in ("GOOD", "BAD")]
    n_good    = sum(1 for r in valid if r["status"] == "GOOD")
    n_bad     = sum(1 for r in valid if r["status"] == "BAD")
    n_near    = sum(1 for r in rows  if r["status"] == "NEAR")
    n_other   = len(rows) - len(valid) - n_near

    res_vals  = [r["res"] for r in valid if r.get("res") is not None]
    rms       = (math.sqrt(sum(v**2 for v in res_vals) / len(res_vals))
                 if res_vals else None)

    # Count how many BAZ observations the locator actually used
    n_used    = sum(1 for r in rows if r.get("baz_used") is True)
    n_not_used = sum(1 for r in rows if r.get("baz_used") is False)
    locator_info_available = (n_used + n_not_used) > 0

    print()
    print(f"Picks with DFX backazimuth:  {len(rows)}")
    print(f"  Distance ≥ {min_dist}°  (reliable range): "
          f" {n_good} GOOD  /  {n_bad} BAD  /  {len(valid)} total")
    print(f"  Distance <  {min_dist}° (NEAR — unreliable):  {n_near}")
    if mode == "xml":
        print(f"  No station coordinates:  {n_other}")
    if locator_info_available:
        print(f"  Used by locator (backazimuthUsed=true):  {n_used}  "
              f"/ rejected or not passed:  {n_not_used}")
    else:
        print("  Used by locator:  unknown (no arrival info for this origin)")
    if rms is not None:
        print(f"  BAZ residual RMS (reliable range):  {rms:.1f}°")
    else:
        print("  BAZ residual RMS:  n/a")
    print()

    if n_near > 0 and n_near == len(rows):
        print(
            "WARNING: ALL stations are within the near-distance threshold.\n"
            "         This is likely a local/regional event.  DFX BAZ is\n"
            "         unreliable at steep P-wave incidence angles.\n"
            f"         Retest with a teleseismic event (D ≥ {min_dist}°)."
        )
    elif n_bad > n_good and len(valid) > 3:
        print(
            f"WARNING: More BAD ({n_bad}) than GOOD ({n_good}) measurements.\n"
            "         Possible causes:\n"
            "           • Wrong origin — try --origin-id to select another\n"
            "           • Station orientation errors in SeisComP inventory\n"
            "           • ~180° offset → DFX sign convention issue\n"
            "           • Poor signal quality at those stations"
        )
    elif n_good > 0 and rms is not None and rms <= threshold:
        print(f"DFX backazimuth measurements look GOOD "
              f"(RMS {rms:.1f}°, {n_good}/{len(valid)} within ±{threshold:.0f}°).")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_stations = os.path.join(script_dir, "..", "data", "stations.dat")

    parser = argparse.ArgumentParser(
        description=(
            "Verify DFX backazimuth measurements from a SeisComP XML file "
            "or a scolv arrival CSV export."
        )
    )
    parser.add_argument(
        "input_file",
        help="SeisComP XML (scxmldump -fPAO) or scolv arrival CSV export",
    )
    parser.add_argument(
        "--stations", default=default_stations, metavar="PATH",
        help=f"Path to stations.dat [default: {default_stations}]  (XML mode only)",
    )
    parser.add_argument(
        "--origin-id", default=None, metavar="ID",
        help="Partial origin ID to use instead of preferred origin  (XML mode only)",
    )
    parser.add_argument(
        "--min-dist", type=float, default=20.0, metavar="DEG",
        help="Stations closer than this (°) are NEAR / excluded from RMS  [default: 20]",
    )
    parser.add_argument(
        "--threshold", type=float, default=30.0, metavar="DEG",
        help="Residual > this → BAD  [default: 30]",
    )
    args = parser.parse_args()

    if not os.path.isfile(args.input_file):
        sys.exit(f"ERROR: file not found: {args.input_file}")

    if _is_csv(args.input_file):
        print("Input format: scolv arrival CSV", file=sys.stderr)
        print(f"Min reliable distance: {args.min_dist}°  "
              f"BAD threshold: ±{args.threshold}°", file=sys.stderr)
        process_csv(args.input_file, args.min_dist, args.threshold)
    else:
        print("Input format: SeisComP XML", file=sys.stderr)
        print(f"Min reliable distance: {args.min_dist}°  "
              f"BAD threshold: ±{args.threshold}°", file=sys.stderr)
        process_xml(args.input_file, args.stations, args.origin_id,
                    args.min_dist, args.threshold)


if __name__ == "__main__":
    main()
