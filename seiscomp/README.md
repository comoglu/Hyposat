# SeisComP Integration for HYPOSAT

This directory contains scripts to use HYPOSAT as a seismic event locator
within [SeisComP](https://www.seiscomp.de/) via the **ExternalLocator** plugin
(`locext`). Once configured, HYPOSAT appears as a selectable locator in `scolv`
alongside the built-in locators (LOCSAT, NonLinLoc, etc.).

## How it works

SeisComP's ExternalLocator plugin pipes an event's origin and picks as an XML
document to an external script, waits for the script to return an updated origin
as XML, and displays the result in `scolv`. The wrapper script
`hyposat_wrapper.py` bridges the two:

```
scolv → locext plugin → hyposat_wrapper.py → hyposat binary → scolv
```

1. `locext` serialises the current origin + picks to SeisComP XML and writes it
   to the wrapper's stdin.
2. `hyposat_wrapper.py` converts the picks to a `hyposat-in` file, writes a
   `hyposat-parameter` file, and runs the `hyposat` binary.
3. The wrapper parses `hyposat-out` and returns an updated `<Origin>` element
   as SeisComP XML on stdout.
4. `locext` deserialises the origin and passes it back to `scolv`.

## Requirements

- SeisComP ≥ 5 with the `locext` plugin available  
  (`/opt/seiscomp/lib/plugins/locext.so` or similar)
- Python 3.6+  (standard library only — no extra packages needed)
- The compiled `hyposat` binary (see the build instructions in `../README.md`);
  the pre-compiled binary in `../bin/hyposat` works on Rocky 9 and Ubuntu 24.04 LTS

## Installation

### 1. Build HYPOSAT

Follow the build instructions in the top-level `README.md`. After building,
the `hyposat` binary should be at `../bin/hyposat` relative to this directory
(i.e. `<repo>/bin/hyposat`).

### 2. Enable the ExternalLocator plugin in SeisComP

Edit `/opt/seiscomp/etc/global.cfg` (or your SeisComP operator configuration)
and add `locext` to the plugins list and register the HYPOSAT profile:

```ini
plugins = ..., locext

ExternalLocator.profiles = Hyposat:/path/to/Hyposat/seiscomp/hyposat_wrapper.py
```

Replace `/path/to/Hyposat` with the actual path to this repository.

If you already have other ExternalLocator profiles, separate them with commas:

```ini
ExternalLocator.profiles = OtherLocator:/path/to/other.py, Hyposat:/path/to/Hyposat/seiscomp/hyposat_wrapper.py
```

### 3. Make the wrapper executable

```bash
chmod +x /path/to/Hyposat/seiscomp/hyposat_wrapper.py
```

### 4. Verify the wrapper works standalone

```bash
python3 hyposat_wrapper.py --help
python3 hyposat_wrapper.py < ../examples/hyposat-in.tele   # (not exactly, but for testing)
```

Or pipe a real SeisComP XML event:

```bash
scxmldump -fPAMp -E <eventID> -o event.xml
python3 hyposat_wrapper.py < event.xml
```

### 5. Restart scolv

```bash
seiscomp restart scolv
# or just relaunch scolv
```

HYPOSAT will now appear in the **Locator** drop-down in `scolv`.

## Station coordinates

HYPOSAT reads station coordinates from its own `stations.dat` file in NEIC
format. The file shipped with this repository (`../data/stations.dat`) contains
the NORSAR global station list plus stations extracted from the SeisComP
inventory at Geoscience Australia. You will likely need to add your own network.

### Adding stations from your SeisComP inventory

```bash
# Export your SeisComP inventory
scxmldump -fI -o inventory.xml

# Append your stations to stations.dat
python3 make_stations_dat.py inventory.xml >> ../data/stations.dat
```

To replace all SeisComP-derived entries and regenerate from scratch, keep only
the original NORSAR section (up to and including the `ZZZZ` sentinel line) and
then append:

```bash
# Find the NORSAR sentinel line number
grep -n "^ZZZZ" ../data/stations.dat

# Keep only NORSAR entries, then append freshly generated ones
head -n <line_number> ../data/stations.dat > /tmp/stations_norsar.dat
python3 make_stations_dat.py inventory.xml >> /tmp/stations_norsar.dat
cp /tmp/stations_norsar.dat ../data/stations.dat
```

### NEIC station file format

The format that HYPOSAT detects automatically is:

```
SSSSS<sep>DDMMss.sN/SDDDMMss.sE/W<elev_m><description>
```

Fortran FORMAT: `(A5, a1, i2, i2, f4.1, a1, i3, i2, f4.1, a1, f7.1, a48)`

- Columns 1–5: station code (left-justified)
- Column 6: any separator character (typically a space)
- Columns 7–15: latitude — `DDMMss.sN` or `DDMMss.sS`
- Columns 16–25: longitude — `DDDMMss.sE` or `DDDMMss.sW`
  (3-digit degree field; a leading space appears for longitudes < 100°)
- Columns 26–32: elevation in **metres**, format `f7.1`
- Columns 33–80: free-text description (48 chars)

**Key constraint**: HYPOSAT detects the format by checking that position 15
is `N`/`S` and position 25 is `E`/`W`. Any extra space between the latitude
and longitude fields will shift the `E`/`W` to position 26 and break detection.

## Configuration options

The wrapper accepts the following command-line options, which SeisComP's
ExternalLocator can pass via `scolv` profile parameters:

| Option | Description |
|--------|-------------|
| `--fixed-depth=<km>` | Fix the hypocenter depth to this value |
| `--max-dist=<deg>` | Ignore arrivals beyond this epicentral distance |
| `--ignore-initial-location` | Let HYPOSAT determine the starting location |

These can be set in the `scolv` locator profile configuration.

### Environment variables

| Variable | Default | Description |
|----------|---------|-------------|
| `HYPOSAT_BIN` | `../bin/hyposat` | Path to the hyposat executable |
| `HYPOSAT_DATA_DIR` | `../data` | Directory with velocity models and `stations.dat` |
| `HYPOSAT_STATION_FILE` | `$HYPOSAT_DATA_DIR/stations.dat` | Station file path |
| `HYPOSAT_MODEL` | `ak135_A` | Velocity model name (must exist in `HYPOSAT_DATA_DIR`) |

Paths are resolved relative to the location of `hyposat_wrapper.py`.

## Velocity models

The wrapper generates a `hyposat-parameter` file that enables CRUST 1.0 crustal
corrections by default. Set the model with the `HYPOSAT_MODEL` environment
variable (default: `ak135_A`).

### Global 1-D models (available in `../data/`)

| Model name | Description |
|------------|-------------|
| `ak135_A` | AK135 — recommended default for global and teleseismic locations |
| `ek137_A` | EK137 — extended AK135 variant |
| `iasp91_A` | IASP91 — the original ISC standard model |
| `iasp91a_A` | IASP91a — as IASP91 but with a different inner core |
| `prem_A` | PREM (Preliminary Reference Earth Model) |
| `sp6_A` | SP6 regional model |
| `jb_A` | Jeffreys-Bullen — the historical standard |

### Regional models (available in `../data/`)

These are pre-computed tau-spline tables for specific regions:

| Model name | Region |
|------------|--------|
| `barents16_A` | Barents Sea region |
| `barey_A` | Barents Sea / Svalbard (Y-component) |
| `barez_A` | Barents Sea / Svalbard (Z-component) |
| `bergen_A` | Bergen (Norway) regional model |
| `fescan_A` | Fennoscandian Shield |

Regional models are used as the `GLOBAL MODEL` setting and are best suited for
events within the region they were derived for.

### Selecting a model

```bash
# Use IASP91 instead of AK135
HYPOSAT_MODEL=iasp91_A python3 hyposat_wrapper.py < event.xml

# Or set permanently in your environment / SeisComP profile config
export HYPOSAT_MODEL=iasp91_A
```

### Local / custom velocity model

HYPOSAT also supports a custom flat-layer velocity model via the
`LOCAL OR REGIONAL MODEL` parameter in `hyposat-parameter`. This is not yet
exposed in the wrapper but can be added by extending `format_hyposat_parameter()`
in `hyposat_wrapper.py`. See `../examples/loc.dat` for the format.

### CRUST 1.0 corrections

Crustal corrections using the CRUST 1.0 global model are enabled by default
(`CRUST CORR CODE: 1` in the generated parameter file). This improves location
accuracy particularly for regional distances. The required data files
(`crust1.bnds`, `crust1.vp`, `crust1.vs`) are included in `../data/`.

## Troubleshooting

### "no origin in result document" popup in scolv

The `locext` plugin could not deserialise the wrapper's output. Check that:

- The wrapper exits with code 0 (check `scolv` debug log for `ERROR:` lines)
- The output XML has `<Origin>` as a direct child of `<seiscomp>` (capital O,
  matching SeisComP's internal class name)

### "Cannot find station: XXXXX entry skipped"

The station is missing from `stations.dat`. Export your inventory with
`scxmldump -fI` and run `make_stations_dat.py` to add the missing stations
(see above).

### "Cannot read station file (wrongly formatted)"

A line in `stations.dat` does not conform to the NEIC format. Common cause:
an extra space between the latitude hemisphere character and the longitude
degrees shifts the `E`/`W` character from column 25 to column 26. Verify the
format with:

```bash
python3 -c "
with open('../data/stations.dat') as f:
    for i, line in enumerate(f, 1):
        if len(line) >= 25 and line[14] in 'NS' and line[24] not in 'EW ':
            print(f'line {i}: {repr(line.rstrip())}')
" 
```

### scolv hangs when relocating

Earlier versions of the wrapper inherited `scolv`'s stdin pipe. The current
wrapper uses `stdin=subprocess.DEVNULL` when invoking `hyposat`. If you see
hangs, verify you have the latest `hyposat_wrapper.py`.

### Viewing HYPOSAT diagnostic output

HYPOSAT's stdout (iteration details, station warnings) is forwarded to the
wrapper's stderr. When running from `scolv`, look in the `scolv` terminal or
log file for lines prefixed with `HYPOSAT:`.

To test interactively:

```bash
python3 hyposat_wrapper.py < event.xml > /dev/null
```

## Verifying DFX backazimuth measurements

After enabling `fx = DFX` in scautopick, use `check_dfx_baz.py` to verify
that the measured backazimuths are physically reasonable for a given event:

```bash
# Export the event XML
scxmldump -fPAO -E <eventID> -o event.xml

# Run the check (uses preferred origin and ../data/stations.dat by default)
python3 check_dfx_baz.py event.xml

# Override the minimum distance threshold (default 20°)
python3 check_dfx_baz.py event.xml --min-dist 30

# Use a specific origin instead of the preferred one
python3 check_dfx_baz.py event.xml --origin-id 20260405013508.279513.88656
```

The script computes theoretical backazimuths from great-circle geometry and
compares them against the DFX-measured values, reporting per-station residuals
and an RMS.

**Interpreting the output:**

| Status | Meaning |
|--------|---------|
| `GOOD` | Residual within `--threshold` (default ±30°) |
| `BAD` | Residual exceeds threshold — bad measurement or wrong origin |
| `NEAR` | Station is closer than `--min-dist` — DFX unreliable here |
| `NO_COORDS` | Station missing from stations.dat |

DFX is designed for **teleseismic** P waves (distance ≥ 20°). At short
distances the P wave arrives nearly vertically, so horizontal particle motion
is minimal and the backazimuth estimate is unreliable. Always test DFX on
teleseismic events (M5.5+, distance > 30° at most recording stations).

**Note on the `rect_raw` column:** DFX stores `2 × λ₁` (twice the largest
covariance eigenvalue) in the `DFX:rectilinearity` Pick comment. This is
**not** the normalised Jurkevics (1988) value (0–1). Use it only as a
relative quality ranking within the same event — higher values indicate
stronger signal polarization.

## Files in this directory

| File | Description |
|------|-------------|
| `hyposat_wrapper.py` | SeisComP ExternalLocator wrapper (main integration script) |
| `make_stations_dat.py` | Convert SeisComP inventory XML to HYPOSAT `stations.dat` |
| `check_dfx_baz.py` | Verify DFX backazimuth measurements against theoretical values |
