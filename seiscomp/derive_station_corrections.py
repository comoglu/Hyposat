#!/usr/bin/env python3
"""
derive_station_corrections.py — Derive empirical station corrections from
SeisComP preferred-origin XML files.

The residuals already stored in each <arrival>'s <timeResidual> element are
extracted directly — no travel-time recomputation needed.  Corrections are
computed as:

    correction = -mean(residuals)

so that adding the correction to the observed travel time removes the
systematic bias.

Usage
-----
  # Pass a directory (recommended for 1000+ files — avoids shell ARG_MAX)
  python3 derive_station_corrections.py --dir /data/xml/ \\
      --hyposat-cor stations.cor \\
      --locsat-cor  locsat_corrections.dat \\
      --report      corrections_report.txt

  # Or pass files explicitly
  python3 derive_station_corrections.py event1.xml event2.xml ...

  # Glob (watch out for ARG_MAX with thousands of files)
  python3 derive_station_corrections.py /data/xml/*.xml

  # Tighten quality filters
  python3 derive_station_corrections.py --dir /data/xml/ \\
      --max-rms 0.8 --max-gap 150 --min-phases 10 --min-events 30

Output files
------------
  stations.cor         HYPOSAT format: STA  Vp  Vs  dtp  dts
  locsat_corrections.dat  SeisComP/LOCSAT format: NET STA PHASE correction
  corrections_report.txt  Human-readable summary per station

Notes
-----
  * A station is only written if it has at least --min-events observations
    for that phase (default: 15).
  * The correction std deviation is written to the report.  Stations with
    σ > --max-sigma (default: 1.5 s) are flagged as UNSTABLE and excluded
    from the output files.
  * Near-surface velocities for the elevation term default to 5.8 / 3.36 km/s
    (IASP91 surface Pg/Sg).  Override per-station with --vp-default / --vs-default.
  * The LOCSAT corrections file format expected by SeisComP is one line per
    station/phase: "NET STA PHASE correction_s"
"""

import argparse
import glob
import math
import os
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict

# ---------------------------------------------------------------------------
# Known SeisComP XML namespaces
# ---------------------------------------------------------------------------

_KNOWN_NS = [
    "http://geofon.gfz-potsdam.de/ns/seiscomp3-schema/0.11",
    "http://geofon.gfz-potsdam.de/ns/seiscomp3-schema/0.12",
    "http://geofon.gfz-potsdam.de/ns/seiscomp3-schema/0.13",
    "http://geofon.gfz.de/ns/seiscomp-schema/0.14",
]


def q(ns: str, tag: str) -> str:
    return f"{{{ns}}}{tag}"


def _text(el, ns: str, *tags):
    node = el
    for t in tags:
        if node is None:
            return None
        node = node.find(q(ns, t))
    return node.text.strip() if node is not None and node.text else None


def _float(el, ns: str, *tags):
    t = _text(el, ns, *tags)
    try:
        return float(t) if t is not None else None
    except ValueError:
        return None


# ---------------------------------------------------------------------------
# Per-event extraction
# ---------------------------------------------------------------------------

def _detect_ns(root) -> str:
    """Detect SeisComP namespace from the already-parsed root element."""
    tag = root.tag
    if tag.startswith("{"):
        return tag[1:tag.index("}")]
    return _KNOWN_NS[-1]


def extract_event(path: str,
                  max_rms: float, max_gap: float, min_phases: int,
                  phases_wanted: set) -> list[dict]:
    """Return a list of {net, sta, phase, residual} dicts for one XML file.

    Returns [] if the preferred origin fails the quality filters.
    Detects namespace from the root element (single parse, no double-read).
    """
    try:
        tree = ET.parse(path)
    except ET.ParseError as e:
        print(f"  SKIP {os.path.basename(path)}: parse error: {e}",
              file=sys.stderr)
        return []

    root = tree.getroot()
    ns = _detect_ns(root)
    ep = root.find(q(ns, "EventParameters"))
    if ep is None:
        return []

    # Build pick → (net, sta) map
    pick_sta = {}
    for pick in ep.iter(q(ns, "pick")):
        pid = pick.get("publicID")
        wf  = pick.find(q(ns, "waveformID"))
        if pid and wf is not None:
            pick_sta[pid] = (
                wf.get("networkCode", "?"),
                wf.get("stationCode",  "?"),
            )

    # Build origin map
    origins = {}
    for orig in ep.iter(q(ns, "origin")):
        oid = orig.get("publicID")
        if oid:
            origins[oid] = orig

    # Find preferred origin
    preferred = None
    event = ep.find(q(ns, "event"))
    if event is not None:
        pref_id = _text(event, ns, "preferredOriginID")
        if pref_id and pref_id in origins:
            preferred = origins[pref_id]
    if preferred is None:
        if origins:
            preferred = next(iter(origins.values()))
        else:
            return []

    # Quality filters
    rms = _float(preferred, ns, "quality", "standardError")
    gap = _float(preferred, ns, "quality", "azimuthalGap")
    n   = _float(preferred, ns, "quality", "usedPhaseCount")

    if rms is not None and rms > max_rms:
        return []
    if gap is not None and gap > max_gap:
        return []
    if n is not None and n < min_phases:
        return []

    lat = _float(preferred, ns, "latitude",  "value")
    lon = _float(preferred, ns, "longitude", "value")
    dep = _float(preferred, ns, "depth",     "value")
    if lat is None or lon is None:
        return []

    rows = []
    for arr in preferred.iter(q(ns, "arrival")):
        # Skip zero-weight (excluded) arrivals
        w = _float(arr, ns, "weight")
        if w is not None and w == 0.0:
            continue

        phase   = _text(arr, ns, "phase")
        res     = _float(arr, ns, "timeResidual")
        pick_id = _text(arr, ns, "pickID")

        if phase is None or res is None or pick_id is None:
            continue

        # Normalise phase to P or S family
        phase_upper = phase.upper()
        if phases_wanted and phase_upper not in phases_wanted:
            continue

        net, sta = pick_sta.get(pick_id, ("?", "?"))
        if sta == "?":
            continue

        dist = _float(arr, ns, "distance")

        rows.append({
            "net":   net,
            "sta":   sta,
            "phase": phase_upper,
            "res":   res,
            "dist":  dist,
            "file":  os.path.basename(path),
            "lat":   lat,
            "lon":   lon,
            "dep":   dep,
        })

    return rows


# ---------------------------------------------------------------------------
# Statistics helpers
# ---------------------------------------------------------------------------

def _stats(values: list[float]) -> tuple[float, float, float]:
    """Return (mean, std_dev, median)."""
    n = len(values)
    if n == 0:
        return 0.0, 0.0, 0.0
    mean = sum(values) / n
    std  = math.sqrt(sum((v - mean) ** 2 for v in values) / n) if n > 1 else 0.0
    s    = sorted(values)
    med  = s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2.0
    return mean, std, med


def _sigma_clip(values: list[float], nsigma: float = 2.5) -> list[float]:
    """Remove outliers beyond nsigma * std from the mean (one pass)."""
    if len(values) < 4:
        return values
    mean, std, _ = _stats(values)
    if std == 0.0:
        return values
    return [v for v in values if abs(v - mean) <= nsigma * std]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Derive empirical station corrections from SeisComP preferred-"
            "origin XML files.  Residuals are taken directly from the "
            "<timeResidual> elements — no travel-time recomputation needed."
        )
    )
    parser.add_argument("xml_files", nargs="*", metavar="FILE",
                        help="SeisComP XML files (preferred-origin exports)")
    parser.add_argument("--dir", metavar="DIR", default=None,
                        help="Directory of *.xml files — use instead of listing "
                             "thousands of paths on the command line")
    parser.add_argument("--hyposat-cor", metavar="PATH", default=None,
                        help="Write HYPOSAT stations.cor to this file")
    parser.add_argument("--locsat-cor", metavar="PATH", default=None,
                        help="Write SeisComP/LOCSAT corrections to this file")
    parser.add_argument("--report", metavar="PATH", default=None,
                        help="Write human-readable report to this file "
                             "(default: print to stdout)")
    parser.add_argument("--max-rms", type=float, default=1.0, metavar="S",
                        help="Reject events with RMS > this (s)  [default: 1.0]")
    parser.add_argument("--max-gap", type=float, default=180.0, metavar="DEG",
                        help="Reject events with azimuthal gap > this (°) "
                             "[default: 180]")
    parser.add_argument("--min-phases", type=int, default=8, metavar="N",
                        help="Reject events with fewer than N defining phases "
                             "[default: 8]")
    parser.add_argument("--min-events", type=int, default=15, metavar="N",
                        help="Only write corrections for stations with ≥ N "
                             "observations  [default: 15]")
    parser.add_argument("--max-sigma", type=float, default=1.5, metavar="S",
                        help="Flag UNSTABLE if residual std dev > this (s) "
                             "after sigma-clipping  [default: 1.5]")
    parser.add_argument("--phases", default="P,S,Pn,Pg,Sn,Sg",
                        help="Comma-separated list of phases to include "
                             "[default: P,S,Pn,Pg,Sn,Sg]")
    parser.add_argument("--vp-default", type=float, default=5.8, metavar="KM/S",
                        help="Default near-surface Vp for HYPOSAT elevation "
                             "correction  [default: 5.8]")
    parser.add_argument("--vs-default", type=float, default=3.36, metavar="KM/S",
                        help="Default near-surface Vs for HYPOSAT elevation "
                             "correction  [default: 3.36]")
    parser.add_argument("--no-sigma-clip", action="store_true",
                        help="Disable sigma-clipping of outlier residuals")
    args = parser.parse_args()

    # Build file list from positional args and/or --dir
    xml_files = list(args.xml_files)
    if args.dir:
        xml_files += sorted(glob.glob(os.path.join(args.dir, "*.xml")))
    if not xml_files:
        sys.exit("ERROR: no XML files specified — use positional arguments or --dir")

    phases_wanted = {p.strip().upper() for p in args.phases.split(",")}

    # -------------------------------------------------------------------
    # Collect residuals from all XML files
    # -------------------------------------------------------------------
    # residuals[(net, sta, phase)] = [res1, res2, ...]
    residuals: dict[tuple, list[float]] = defaultdict(list)
    n_events_used = 0
    n_events_skipped = 0
    n_total = len(xml_files)

    for i, path in enumerate(xml_files, 1):
        if i % 100 == 0 or i == n_total:
            print(f"  [{i}/{n_total}] used={n_events_used} skipped={n_events_skipped}",
                  file=sys.stderr)
        if not os.path.isfile(path):
            print(f"WARNING: not found: {path}", file=sys.stderr)
            continue
        rows = extract_event(path,
                             max_rms=args.max_rms,
                             max_gap=args.max_gap,
                             min_phases=args.min_phases,
                             phases_wanted=phases_wanted)
        if rows:
            n_events_used += 1
            for r in rows:
                key = (r["net"], r["sta"], r["phase"])
                residuals[key].append(r["res"])
        else:
            n_events_skipped += 1

    print(f"Events used: {n_events_used}  skipped: {n_events_skipped}",
          file=sys.stderr)
    print(f"Station/phase combinations seen: {len(residuals)}",
          file=sys.stderr)

    if n_events_used == 0:
        sys.exit("ERROR: no events passed the quality filters — "
                 "try relaxing --max-rms or --max-gap")

    # -------------------------------------------------------------------
    # Compute corrections
    # -------------------------------------------------------------------
    # corrections[(net, sta, phase)] = {mean, std, median, n, correction, stable}
    corrections = {}
    for (net, sta, phase), vals in sorted(residuals.items()):
        if len(vals) < args.min_events:
            continue

        clipped = _sigma_clip(vals) if not args.no_sigma_clip else vals
        mean, std, med = _stats(clipped)
        correction = -mean   # add to observed time to remove bias
        stable = (std <= args.max_sigma)

        corrections[(net, sta, phase)] = {
            "mean":       mean,
            "std":        std,
            "median":     med,
            "n_raw":      len(vals),
            "n_clipped":  len(clipped),
            "correction": correction,
            "stable":     stable,
        }

    n_stable   = sum(1 for v in corrections.values() if v["stable"])
    n_unstable = sum(1 for v in corrections.values() if not v["stable"])
    print(f"Corrections derived: {len(corrections)}  "
          f"({n_stable} stable, {n_unstable} unstable/excluded)",
          file=sys.stderr)

    # -------------------------------------------------------------------
    # Build report
    # -------------------------------------------------------------------
    report_lines = []
    report_lines.append(
        f"Station corrections derived from {n_events_used} events\n"
        f"Filters: max_rms={args.max_rms}s  max_gap={args.max_gap}°  "
        f"min_phases={args.min_phases}  min_events={args.min_events}\n"
        f"Sigma-clipping: {'disabled' if args.no_sigma_clip else '2.5σ'}\n"
    )
    report_lines.append(
        f"{'Station':<14} {'Phase':<6} {'N':>5} {'N_clip':>6} "
        f"{'Mean res':>9} {'Std':>7} {'Median':>8} "
        f"{'Correction':>11}  {'Status'}"
    )
    report_lines.append("-" * 80)

    for (net, sta, phase), c in sorted(corrections.items()):
        sta_label = f"{net}.{sta}"
        status = "OK" if c["stable"] else "UNSTABLE (excluded)"
        report_lines.append(
            f"{sta_label:<14} {phase:<6} {c['n_raw']:>5} {c['n_clipped']:>6} "
            f"{c['mean']:>+9.3f} {c['std']:>7.3f} {c['median']:>+8.3f} "
            f"{c['correction']:>+11.3f}  {status}"
        )

    report_text = "\n".join(report_lines) + "\n"

    if args.report:
        with open(args.report, "w") as f:
            f.write(report_text)
        print(f"Report written to {args.report}", file=sys.stderr)
    else:
        print(report_text)

    # -------------------------------------------------------------------
    # HYPOSAT stations.cor
    # -------------------------------------------------------------------
    if args.hyposat_cor:
        # Build per-station P and S corrections
        sta_p = {}   # (net, sta) → correction
        sta_s = {}

        for (net, sta, phase), c in corrections.items():
            if not c["stable"]:
                continue
            key = (net, sta)
            # Map phase families to P or S
            if phase in ("P", "PN", "PG", "PB", "PKP", "PKKP"):
                if key not in sta_p:
                    sta_p[key] = []
                sta_p[key].append(c["correction"])
            elif phase in ("S", "SN", "SG", "SB", "LG"):
                if key not in sta_s:
                    sta_s[key] = []
                sta_s[key].append(c["correction"])

        all_stas = set(sta_p) | set(sta_s)
        lines = [
            "* HYPOSAT station corrections derived by derive_station_corrections.py",
            f"* Events used: {n_events_used}",
            f"* Format: STA  Vp  Vs  dtp  dts",
            "*",
        ]
        for (net, sta) in sorted(all_stas):
            dtp = sum(sta_p.get((net, sta), [0.0])) / max(1, len(sta_p.get((net, sta), [0.0])))
            dts = sum(sta_s.get((net, sta), [0.0])) / max(1, len(sta_s.get((net, sta), [0.0])))
            vp  = args.vp_default
            vs  = args.vs_default
            lines.append(
                f"{sta:<6}  {vp:.2f}  {vs:.2f}  {dtp:+.3f}  {dts:+.3f}"
                f"  * {net}"
            )

        with open(args.hyposat_cor, "w") as f:
            f.write("\n".join(lines) + "\n")
        print(f"HYPOSAT corrections written to {args.hyposat_cor}  "
              f"({len(all_stas)} stations)", file=sys.stderr)

    # -------------------------------------------------------------------
    # LOCSAT / SeisComP corrections  (.stacor format)
    #
    # SeisComP LOCSAT subtracts the delay from the observed travel time:
    #   t_obs - delay = t_predicted
    # So delay = mean(residuals) = -correction.
    # Format: LOCDELAY NET.STA PHASE N_READINGS DELAY_S
    # File goes in $SEISCOMP_ROOT/share/locsat/tables/<profile>.stacor
    # -------------------------------------------------------------------
    if args.locsat_cor:
        lines = [
            "# SeisComP/LOCSAT station corrections derived by "
            "derive_station_corrections.py",
            f"# Events used: {n_events_used}",
            "# Format: LOCDELAY NET.STA PHASE N_READINGS DELAY_S",
            "# Place this file at $SEISCOMP_ROOT/share/locsat/tables/<profile>.stacor",
            "#",
        ]
        for (net, sta, phase), c in sorted(corrections.items()):
            if not c["stable"]:
                continue
            # delay = mean residual = -correction (LOCSAT subtracts this value)
            delay = c["mean"]
            lines.append(
                f"LOCDELAY {net}.{sta} {phase} {c['n_clipped']} {delay:+.3f}"
            )

        with open(args.locsat_cor, "w") as f:
            f.write("\n".join(lines) + "\n")
        print(f"LOCSAT corrections written to {args.locsat_cor}  "
              f"({n_stable} entries)", file=sys.stderr)


if __name__ == "__main__":
    main()
