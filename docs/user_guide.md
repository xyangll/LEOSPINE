# SPINE User Guide

This guide covers all features of SPINE and how to use them effectively.

---

## Table of Contents

1. [Getting Started](#getting-started)
2. [Constellation Design](#constellation-design)
3. [Observation Simulation](#observation-simulation)
4. [Positioning Computation](#positioning-computation)
5. [Settings](#settings)
6. [CLI Tools](#cli-tools)
7. [Tips & Troubleshooting](#tips--troubleshooting)

---

## Getting Started

### Launching the Application

```bash
# Default: Auto-detect GUI backend
python main.py

# Force Qt GUI
python main.py --gui qt

# Debug mode (shows detailed errors)
python main.py --debug
```

### Interface Overview

The main window has four tabs:
- **Constellation**: Design and visualize satellite constellations
- **Simulation**: Generate simulated GNSS observations
- **Positioning**: Compute positions from observations
- **Settings**: Configure appearance and preferences

---

## Constellation Design

### Creating a Walker-Delta Constellation

1. **Set Orbit Parameters**
   - Semi-major axis (m): Orbital radius, e.g., `26560000` for MEO
   - Eccentricity: Orbit shape, typically `0.0` - `0.1`
   - Inclination (°): Orbit tilt, e.g., `55°` for GPS-like
   - RAAN (°): Right ascension of ascending node
   - Arg. of Perigee (°): Perigee orientation
   - Mean Anomaly (°): Initial satellite position

2. **Configure Walker-Delta**
   - **P**: Number of orbital planes
   - **S**: Satellites per plane
   - **F**: Phasing factor (0 to P-1)
   - **RAAN spread**: Angular spread of planes (default 360°)
   - **Epoch**: Reference date for TLE

3. **Generate TLE**
   - Click **✨ Generate** to create the constellation
   - Total satellites = P × S

### Example Configurations

| Constellation | P | S | F | Inc | Semi-major (km) |
|--------------|---|---|---|-----|-----------------|
| GPS-like | 6 | 4 | 1 | 55° | 26,560 |
| Iridium-like | 6 | 11 | 2 | 86° | 7,080 |
| LEO Walker | 8 | 10 | 1 | 53° | 6,878 |

### 3D Visualization

1. **Configure Display**
   - Duration: Simulation time (seconds)
   - Step: Position sampling interval
   - Frame: 
     - *Inertial*: Earth rotates, orbit appears closed
     - *Fixed*: Earth fixed, shows ground track
   - Show Beam: Display antenna coverage cones

2. **Refresh View**
   - Click **🌐 Refresh View** to update Cesium display
   - Click **🚀 Browser** to open in external browser

3. **Color Coding**
   - Each satellite automatically gets a unique color
   - Orbit path and beam use matching colors

### Visibility Analysis

1. **Choose Analysis Type**
   - **Instant**: Single-epoch snapshot
   - **Average**: Time-averaged over period

2. **Set Parameters**
   - Elevation mask: Minimum visible elevation
   - Grid resolution: Lat/lon grid spacing (smaller = slower)

3. **For Average Mode**
   - Set start/end datetime
   - Set time step (seconds between samples)

4. **View Results**
   - Global satellite visibility heatmap
   - PDOP (Position Dilution of Precision) map
   - Switch to "📊 Visibility Analysis" tab

---

## Observation Simulation

### Basic Setup

1. **Load TLE File**
   - Click **📂** next to TLE Path
   - Select a `.tle` file

2. **Select Satellites** (Optional)
   - After loading TLE file, select which satellites to simulate
   - Use "Select All" / "Deselect All" buttons
   - Check/uncheck individual satellites
   - Only selected satellites will be included in simulation

3. **Set Time Parameters**
   - Start: Simulation start time (UTC)
   - Duration: Simulation length (seconds)
   - Step: Sampling interval (seconds)

4. **Configure Station**
   - Choose **LLH** (Lat/Lon/Height) or **XYZ** (ECEF)
   - Enter coordinates
   - Set elevation mask (typically 10°)

### Error Models

Enable/disable error sources:

| Error Model | Description | Typical Value |
|-------------|-------------|---------------|
| Ionosphere | Klobuchar model | 5-15 m |
| Troposphere | Saastamoinen model | 2-5 m |
| Multipath | Simplified model | 0.5-2 m |
| Relativistic | Special/General relativity | ~0.01 m/s |
| Hardware Delay | Receiver delays | Configurable |
| PCO | Phase Center Offset | Configurable |

### Receiver Clock

**Model Mode:**
- Bias: Initial clock offset (meters)
- Drift: Clock rate error (m/s)
- Formula: `clock(t) = bias + drift × t`

**File Mode:**
- Load external clock bias file
- Load external clock drift file
- Format: One value per epoch

### Measurement Noise

Add Gaussian white noise:
- Pseudorange σ: Standard deviation (meters)
- Doppler σ: Standard deviation (m/s)

### Running Simulation

1. Set output file path
2. Click **🚀 Run Simulation**
3. View generated plots:
   - Pseudorange observations
   - Doppler observations

---

## Positioning Computation

### Input Configuration

1. **Load Observations**
   - Click **📂** to select observation file
   - Format: CSV file from simulation

2. **Select Satellites** (Optional)
   - After loading observations, select which satellites to use
   - Use "Select All" / "Deselect All" buttons
   - Check/uncheck individual satellites
   - Only selected satellites will be used for positioning

3. **Set True Position** (for error analysis)
   - Enter known position (LLH or XYZ)
   - Used to compute positioning error

### Positioning Strategy

Choose one of four methods:

| Strategy | Observations | Notes |
|----------|--------------|-------|
| Pseudorange Only | PR | Standard SPP |
| Doppler Only | Doppler | Velocity-aided |
| PR + Doppler (Fixed) | Both | Equal weighting |
| PR + Doppler (Adaptive) | Both | VCE adaptive weighting |

### Error Correction

Enable error models to correct observations:
- Ionosphere, Troposphere, Multipath
- Relativistic, Hardware Delay, PCO (Phase Center Offset)

### Running Positioning

1. Click **🎯 Run Positioning**
2. View results:
   - **Skyplot**: Satellite positions (Az/El)
   - **Error Plot**: Position error vs time
   - **Residuals**: PR and Doppler residuals

### Interpreting Results

- **Skyplot**: Shows satellite geometry
  - Good geometry: Satellites spread across sky
  - Poor geometry: Satellites clustered
  
- **Position Error**: 
  - 3D error = sqrt(ΔN² + ΔE² + ΔU²)
  - Lower is better

- **Residuals**:
  - Should be zero-mean
  - Large residuals indicate problems

---

## Settings

### Theme Configuration

- **Dark Mode**: Dark background, light text
- **Light Mode**: Light background, dark text
- **Accent Color**: UI highlight color

### Cesium Settings

- **Ion Token**: For 3D terrain (optional)
- Get token from [Cesium Ion](https://cesium.com/ion/)
- **High-precision terrain**: Enable world terrain and high-res imagery

### About Dialog

- Click **ℹ️ About SPINE** button to view:
  - Application logo
  - Version information
  - License information
  - Project description

---

## CLI Tools

### Observation Simulation

**Basic usage:**
```bash
python tools/sim_obs.py \
  --tle InputData/LEO_sat.tle \
  --start 0 \
  --duration 3600 \
  --step 10 \
  --lat 39.8 \
  --lon 119.29 \
  --alt 60 \
  --mask 10 \
  --out OutputData/observations.csv \
  --gps_week 2377
```

**With error models and satellite selection:**
```bash
python tools/sim_obs.py \
  --tle InputData/LEO_sat.tle \
  --lat 39.8 --lon 119.29 --alt 60 \
  --iono --tropo --multipath --relativistic --hardware \
  --prnoise 5.0 --dopnoise 0.1 \
  --prns "1,3,5,7"  # Select specific satellites
```

**Parameters:**
| Parameter | Description |
|-----------|-------------|
| `--tle` | TLE file path (required) |
| `--start` | Start GPS time (seconds, default: 0) |
| `--duration` | Simulation duration (seconds, default: 3600) |
| `--step` | Time step (seconds, default: 10) |
| `--lat` | Station latitude (degrees, required) |
| `--lon` | Station longitude (degrees, required) |
| `--alt` | Station altitude (meters, default: 0) |
| `--mask` | Elevation mask (degrees, default: 10) |
| `--out` | Output file path (default: observations.csv) |
| `--gps_week` | GPS week number |
| `--prnoise` | Pseudorange noise sigma (meters) |
| `--dopnoise` | Doppler noise sigma (m/s) |
| `--iono` / `--no-iono` | Enable/disable ionosphere (default: enabled) |
| `--tropo` / `--no-tropo` | Enable/disable troposphere (default: enabled) |
| `--multipath` | Enable multipath effect |
| `--relativistic` | Enable relativistic effect |
| `--hardware` | Enable hardware delay |
| `--clk-model` | Use clock bias/drift model (default) |
| `--clk-file` | Clock bias file (disables model) |
| `--clk-drift-file` | Clock drift file |
| `--clk-bias` | Clock bias in meters (default: 1000.0) |
| `--clk-drift` | Clock drift in m/s (default: 0.05) |
| `--prns` | Comma-separated PRN list (e.g., "1,3,5,7") |

### Positioning

**Basic usage:**
```bash
python tools/run_positioning.py \
  --obs OutputData/observations.csv \
  --tle InputData/LEO_sat.tle \
  --use_pr --use_dop \
  --init_x 0 --init_y 0 --init_z 0
```

**Parameters:**
| Parameter | Description |
|-----------|-------------|
| `--obs` | Observation file path (required) |
| `--nav` | Navigation/ephemeris file (optional) |
| `--tle` | TLE file path (optional, for TLE-based positioning) |
| `--use_pr` | Use pseudorange observations |
| `--use_dop` | Use Doppler observations |
| `--init_x` | Initial X position in ECEF (meters, default: 0) |
| `--init_y` | Initial Y position in ECEF (meters, default: 0) |
| `--init_z` | Initial Z position in ECEF (meters, default: 0) |

**Note:** For advanced options (error correction, satellite selection, etc.), use the GUI or call `run_lsq_positioning()` directly from Python.

---

## Tips & Troubleshooting

### Performance Tips

1. **Large Constellations**
   - Increase CZML step to 60-120s
   - Reduce duration if not needed
   - Use coarser visibility grid (10° vs 5°)

2. **Long Simulations**
   - Consider splitting into segments
   - Increase time step

3. **Slow Visibility Analysis**
   - Use 10° grid resolution
   - Limit time range for average mode

### Common Issues

**Cesium view is blank:**
- Check if TLE file is loaded
- Click "Refresh View"
- Try "Open in Browser"

**No satellites visible:**
- Check elevation mask (try lowering to 5°)
- Verify station coordinates
- Check TLE epoch matches simulation time

**Positioning fails:**
- Need at least 4 visible satellites
- Check observation file format
- Verify true position is reasonable

**Large positioning errors:**
- Enable error correction models
- Check if error models match simulation
- Try different positioning strategy

### File Locations

| Type | Default Location |
|------|------------------|
| TLE Files | `InputData/` |
| Observations | `OutputData/` |
| CZML Output | `app/web/data/` |

---

## Keyboard Shortcuts

| Shortcut | Action |
|----------|--------|
| Ctrl+Q | Quit application |
| F5 | Refresh current view |

---

## Getting Help

- **Issues**: Report bugs on GitHub Issues
- **Discussions**: Ask questions on GitHub Discussions
- **Docs**: Check `docs/` folder for more information
