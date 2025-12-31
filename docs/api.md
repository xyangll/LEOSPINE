# SPINE API Reference

This document provides API documentation for SPINE's core modules.

---

## Table of Contents

1. [core.sim_data](#coresim_data)
2. [core.positioning](#corepositioning)
3. [core.constellation](#coreconstellation)
4. [core.utilities](#coreutilities)
5. [app.czml](#appczml)

---

## core.sim_data

Satellite simulation and observation generation module.

### Class: SatelliteSimulation

```python
class SatelliteSimulation:
    def __init__(self, 
                 enable_iono=True, 
                 enable_tropo=True, 
                 enable_multipath=False,
                 enable_relativistic=False,
                 enable_hardware=False,
                 enable_pco=False)
```

**Parameters:**
- `enable_iono` (bool): Enable ionospheric delay model (Klobuchar)
- `enable_tropo` (bool): Enable tropospheric delay model (Saastamoinen)
- `enable_multipath` (bool): Enable multipath effect model
- `enable_relativistic` (bool): Enable relativistic effect correction
- `enable_hardware` (bool): Enable hardware delay correction
- `enable_pco` (bool): Enable phase center offset (PCO) correction

**Attributes:**
- `tle_satellites` (dict): Loaded TLE satellite data
- `gps_week` (int): Current GPS week number
- `clock` (float): Receiver clock bias (meters)
- `clockdrift` (float): Receiver clock drift (m/s)
- `start_time` (float): Simulation start time (GPS seconds)
- `enable_relativistic` (bool): Relativistic effect flag
- `enable_hardware` (bool): Hardware delay flag
- `enable_pco` (bool): Phase center offset flag

---

### Methods

#### load_tle_file

```python
def load_tle_file(self, tle_path: str) -> int
```

Load satellites from a TLE file.

**Parameters:**
- `tle_path`: Path to TLE file

**Returns:**
- Number of satellites loaded

**Example:**
```python
sim = SatelliteSimulation()
count = sim.load_tle_file("InputData/LEO_sat.tle")
print(f"Loaded {count} satellites")
```

---

#### add_tle_satellite

```python
def add_tle_satellite(self, name: str, line1: str, line2: str, prn: int) -> bool
```

Add a single TLE satellite.

**Parameters:**
- `name`: Satellite name
- `line1`: TLE line 1
- `line2`: TLE line 2
- `prn`: PRN number to assign

**Returns:**
- True if successful

---

#### compute_observations

```python
def compute_observations(self, 
                         prn: int,
                         ephemeris: dict,
                         receiver_pos: np.ndarray,
                         receiver_vel: np.ndarray,
                         gps_time: float,
                         prnoise: float = None,
                         dopnoise: float = None) -> dict
```

Compute observations for a satellite. Supports both TLE-based and ephemeris-based satellite position calculation.

**Parameters:**
- `prn`: Satellite PRN number
- `ephemeris`: Ephemeris data dictionary (None for TLE mode, or dict with 'use_tle' key)
- `receiver_pos`: Receiver ECEF position [x, y, z] (meters)
- `receiver_vel`: Receiver ECEF velocity [vx, vy, vz] (m/s)
- `gps_time`: GPS time of week (seconds)
- `prnoise`: Pseudorange noise sigma (meters, Gaussian white noise)
- `dopnoise`: Doppler noise sigma (m/s, Gaussian white noise)

**Returns:**
```python
{
    'pseudorange': float,      # Pseudorange (meters) with all error models applied
    'doppler': float,          # Doppler (m/s) with relativistic correction if enabled
    'carrier_phase': float,    # Carrier phase (cycles)
    'sat_position': np.ndarray,  # Satellite ECEF position [x, y, z] (meters)
    'sat_velocity': np.ndarray,  # Satellite ECEF velocity [vx, vy, vz] (m/s)
    'emission_time': float,    # Signal emission time (GPS seconds)
    'time_delay': float,       # Signal travel time (seconds)
    'snr': float,              # Signal-to-noise ratio
    'elevation': float,        # Elevation angle (radians)
    'azimuth': float           # Azimuth angle (radians)
}
```

**Note:** Error models (ionosphere, troposphere, multipath, relativistic, hardware delay, PCO) are applied based on the initialization flags.

---

#### compute_iono_delay

```python
def compute_iono_delay(self, 
                       azimuth: float, 
                       elevation: float,
                       user_lat: float, 
                       user_lon: float, 
                       gps_time: float) -> float
```

Compute ionospheric delay using Klobuchar model.

**Parameters:**
- `azimuth`: Satellite azimuth (radians)
- `elevation`: Satellite elevation (radians)
- `user_lat`: User latitude (radians)
- `user_lon`: User longitude (radians)
- `gps_time`: GPS time of week (seconds)

**Returns:**
- Ionospheric delay (meters)

---

#### compute_tropo_delay

```python
def compute_tropo_delay(self, elevation: float, height: float) -> float
```

Compute tropospheric delay using Saastamoinen model.

**Parameters:**
- `elevation`: Satellite elevation (radians)
- `height`: User height above sea level (meters)

**Returns:**
- Tropospheric delay (meters)

---

#### compute_relativistic_effect

```python
def compute_relativistic_effect(self,
                                sat_pos: np.ndarray,
                                sat_vel: np.ndarray,
                                rec_pos: np.ndarray) -> float
```

Compute relativistic effect on Doppler observation.

**Parameters:**
- `sat_pos`: Satellite ECEF position (meters)
- `sat_vel`: Satellite ECEF velocity (m/s)
- `rec_pos`: Receiver ECEF position (meters)

**Returns:**
- Relativistic Doppler contribution (m/s)

---

### Function: main_sim

```python
def main_sim(sim: SatelliteSimulation,
             week: int,
             start_time: float,
             end_time: float,
             step_s: float,
             elev_mask_deg: float,
             output_file: str,
             tle_file: str = None,
             clk_data: np.ndarray = None,
             clk_drift_data: np.ndarray = None,
             receiver_pos: np.ndarray = None,
             prnoise: float = None,
             dopnoise: float = None,
             clk_bias: float = 1000.0,
             clk_drift: float = 0.05,
             selected_prns: list = None)
```

Run observation simulation and write results to CSV file.

**Parameters:**
- `sim`: SatelliteSimulation instance (error models configured in constructor)
- `week`: GPS week number
- `start_time`: Start GPS time of week (seconds)
- `end_time`: End GPS time of week (seconds)
- `step_s`: Time step (seconds)
- `elev_mask_deg`: Elevation mask (degrees)
- `output_file`: Output CSV file path
- `tle_file`: TLE file path (required for TLE-based simulation)
- `clk_data`: Clock bias data array from file (None to use model)
- `clk_drift_data`: Clock drift data array from file (None to use model)
- `receiver_pos`: Receiver ECEF position [x, y, z] (meters)
- `prnoise`: Pseudorange noise sigma (meters, Gaussian white noise)
- `dopnoise`: Doppler noise sigma (m/s, Gaussian white noise)
- `clk_bias`: Clock bias model value (meters, used when clk_data is None)
- `clk_drift`: Clock drift model value (m/s, used when clk_drift_data is None)
- `selected_prns`: List of PRN numbers to simulate (None = all satellites)

**Output Format:**
CSV file with columns: `time`, `prn`, `pseudorange`, `doppler`, `snr`, `elevation`, `azimuth`

**Note:** 
- Error models (ionosphere, troposphere, multipath, relativistic, hardware delay, PCO) are controlled by the `SatelliteSimulation` constructor flags.
- Clock bias/drift can be provided from a file or computed from a model.
- If `selected_prns` is specified, only those satellites are simulated.

---

## core.positioning

Positioning computation module.

### Function: run_lsq_positioning

```python
def run_lsq_positioning(sim_data_file: str,
                        ephemeris_file: str = None,
                        tle_file: str = None,
                        init_pos: tuple = None,
                        height_constraint: float = None,
                        true_position: tuple = None,
                        height_constraint_mode: str = 'hard',
                        adaptive_weights_enable: bool = True,
                        threshold: float = 5.0,
                        min_weight: float = 0.1,
                        enable_iono: bool = True,
                        enable_tropo: bool = True,
                        enable_multipath: bool = False,
                        enable_relativistic: bool = False,
                        enable_hardware: bool = False,
                        enable_pco: bool = False,
                        fit_doppler: bool = False,
                        use_pseudorange: bool = False,
                        use_doppler: bool = True,
                        use_hvce: bool = False,
                        use_DBSCAN: bool = False,
                        progress_callback = None,
                        selected_prns: list = None) -> tuple
```

Run least-squares batch positioning with support for pseudorange, Doppler, or combined strategies.

**Parameters:**
- `sim_data_file`: Path to observation CSV file
- `ephemeris_file`: Ephemeris file path (optional, for ephemeris-based positioning)
- `tle_file`: TLE file path (optional, for TLE-based positioning)
- `init_pos`: Initial position guess [x, y, z] in ECEF (meters)
- `height_constraint`: Height constraint (meters above ellipsoid)
- `true_position`: True position (lat, lon, alt) for error calculation
- `height_constraint_mode`: Height constraint mode: `'hard'`, `'soft'`, or `'both'`
- `adaptive_weights_enable`: Enable adaptive weight adjustment
- `threshold`: Residual threshold for outlier detection
- `min_weight`: Minimum weight limit
- `enable_iono`: Apply ionospheric correction
- `enable_tropo`: Apply tropospheric correction
- `enable_multipath`: Apply multipath correction
- `enable_relativistic`: Apply relativistic correction
- `enable_hardware`: Apply hardware delay correction
- `enable_pco`: Apply phase center offset correction
- `fit_doppler`: Apply Fourier fitting to Doppler observations
- `use_pseudorange`: Use pseudorange for positioning
- `use_doppler`: Use Doppler for positioning
- `use_hvce`: Use hierarchical variance component estimation
- `use_DBSCAN`: Use DBSCAN clustering for outlier detection
- `progress_callback`: Progress callback function `(current, total, message) -> bool`
- `selected_prns`: List of PRN numbers to use (None = all satellites)

**Returns:**
```python
(positioning_object, results_dict)
```

Where `results_dict` contains:
```python
{
    'epochs': list,           # Time epochs
    'positions': list,        # Estimated ECEF positions [x, y, z]
    'clocks': list,           # Clock bias estimates (meters)
    'residuals_pr': list,     # Pseudorange residuals
    'residuals_dop': list,    # Doppler residuals
    'sat_info': list,         # Satellite Az/El info per epoch
    'errors': list,           # Position errors (if true_position given)
    'statistics': dict        # Error statistics (mean, std, RMS, etc.)
}
```

**Example:**
```python
from core.positioning import run_lsq_positioning

pos_obj, result = run_lsq_positioning(
    sim_data_file='observations.csv',
    tle_file='InputData/LEO_sat.tle',
    true_position=(39.8, 119.29, 60.0),  # lat, lon, alt
    use_pseudorange=True,
    use_doppler=True,
    adaptive_weights_enable=True,
    enable_iono=True,
    enable_tropo=True,
    selected_prns=[1, 2, 3, 4, 5]  # Use only these satellites
)

print(f"Mean 3D error: {result['statistics']['mean_3d_error']:.2f} m")
```

---

## core.constellation

Constellation design and analysis module.

### Function: generate_walker_delta_tles

```python
def generate_walker_delta_tles(name_prefix: str,
                               a_m: float,
                               e: float,
                               inc_deg: float,
                               raan0_deg: float,
                               argp_deg: float,
                               mean_anomaly0_deg: float,
                               epoch_dt: datetime,
                               planes: int,
                               sats_per_plane: int,
                               phasing: int,
                               bstar: float = 0.0,
                               catalog_start: int = 80000,
                               raan_spread: float = 360.0) -> list
```

Generate TLE strings for Walker-Delta constellation following STK convention.

**Parameters:**
- `name_prefix`: Satellite name prefix (e.g., "SAT")
- `a_m`: Semi-major axis (meters)
- `e`: Eccentricity
- `inc_deg`: Inclination (degrees)
- `raan0_deg`: Initial RAAN (degrees)
- `argp_deg`: Argument of perigee (degrees)
- `mean_anomaly0_deg`: Initial mean anomaly (degrees)
- `epoch_dt`: Epoch datetime object (UTC)
- `planes`: Number of orbital planes (P)
- `sats_per_plane`: Satellites per plane (S)
- `phasing`: Phasing factor (F, 0 to P-1)
- `bstar`: B* drag coefficient (default: 0.0)
- `catalog_start`: Starting catalog number (default: 80000)
- `raan_spread`: RAAN angular spread (degrees, default: 360.0 for full coverage)

**Returns:**
- List of tuples: `[(name, line1, line2), ...]`

**Walker-Delta Formula:**
- Total satellites: T = P × S
- RAAN for plane p: RAAN₀ + raan_spread × p / P
- Mean anomaly for satellite s in plane p: M₀ + 360° × s/S + 360° × p×F / (P×S)

**Example:**
```python
from core.constellation import generate_walker_delta_tles
from datetime import datetime, timezone

tles = generate_walker_delta_tles(
    name_prefix="LEO",
    a_m=26560000,      # GPS altitude
    e=0.0,
    inc_deg=55.0,
    raan0_deg=0.0,
    argp_deg=0.0,
    mean_anomaly0_deg=0.0,
    epoch_dt=datetime(2025, 1, 1, tzinfo=timezone.utc),
    planes=6,
    sats_per_plane=4,
    phasing=1,
    raan_spread=360.0  # Full 360° coverage
)

print(f"Generated {len(tles)} satellites")
```

---

### Function: compute_visibility_analysis

```python
def compute_visibility_analysis(sim_or_tles,
                                epoch_dt: datetime = None,
                                lat_resolution: float = 5.0,
                                lon_resolution: float = 5.0,
                                elevation_mask: float = 10.0,
                                progress_callback = None) -> dict
```

Compute global satellite visibility and PDOP analysis. High-level function that accepts either a `SatelliteSimulation` object or a TLE list.

**Parameters:**
- `sim_or_tles`: Either a `SatelliteSimulation` object or list of TLE tuples `[(name, line1, line2), ...]`
- `epoch_dt`: Analysis epoch (datetime, UTC). If None, uses current time.
- `lat_resolution`: Latitude grid resolution (degrees, default: 5.0)
- `lon_resolution`: Longitude grid resolution (degrees, default: 5.0)
- `elevation_mask`: Minimum elevation angle (degrees, default: 10.0)
- `progress_callback`: Progress callback function `(percent, message) -> None`

**Returns:**
```python
{
    'lat_grid': np.ndarray,    # Latitude values (degrees)
    'lon_grid': np.ndarray,    # Longitude values (degrees)
    'visibility': np.ndarray,  # Visible satellite count (2D array: lat × lon)
    'pdop': np.ndarray,        # PDOP values (2D array: lat × lon)
    'epoch': datetime,         # Analysis epoch
    'elevation_mask': float    # Elevation mask used
}
```

**Example:**
```python
from core.constellation import compute_visibility_analysis
from core.sim_data import SatelliteSimulation
from datetime import datetime, timezone

# Using SatelliteSimulation object
sim = SatelliteSimulation()
sim.load_tle_file('InputData/LEO_sat.tle')
result = compute_visibility_analysis(
    sim,
    epoch_dt=datetime(2025, 1, 1, 12, 0, 0, tzinfo=timezone.utc),
    lat_resolution=5.0,
    lon_resolution=5.0,
    elevation_mask=10.0
)

# Or using TLE list directly
tles = [("SAT1", "1 12345U ...", "2 12345 ..."), ...]
result = compute_visibility_analysis(tles, epoch_dt=datetime.now(timezone.utc))
```

---

## core.utilities

Coordinate transformation and utility functions.

### Function: geodetic_to_ecef

```python
def geodetic_to_ecef(lat: float, lon: float, alt: float) -> tuple
```

Convert geodetic coordinates to ECEF.

**Parameters:**
- `lat`: Latitude (radians)
- `lon`: Longitude (radians)
- `alt`: Altitude above ellipsoid (meters)

**Returns:**
- Tuple (x, y, z) in meters

---

### Function: ecef_to_geodetic

```python
def ecef_to_geodetic(x: float, y: float, z: float) -> tuple
```

Convert ECEF to geodetic coordinates.

**Parameters:**
- `x`, `y`, `z`: ECEF coordinates (meters)

**Returns:**
- Tuple (lat, lon, alt) where lat/lon in radians, alt in meters

---

### Function: datetime_to_gps_time

```python
def datetime_to_gps_time(dt: datetime) -> tuple
```

Convert datetime to GPS time.

**Parameters:**
- `dt`: Python datetime object (UTC)

**Returns:**
- Tuple (gps_week, gps_seconds)

---

### Function: gps_time_to_datetime

```python
def gps_time_to_datetime(week: int, seconds: float) -> datetime
```

Convert GPS time to datetime.

**Parameters:**
- `week`: GPS week number
- `seconds`: Seconds of week

**Returns:**
- Python datetime object (UTC)

---

## app.czml

CZML generation for Cesium visualization.

### Function: write_czml

```python
def write_czml(tle_path: str,
               start_dt: datetime,
               duration_s: float,
               step_s: float,
               out_path: str,
               accent_rgba: list = None,
               show_beam: bool = False,
               beam_half_angle: float = 21.0,
               use_fixed_frame: bool = False)
```

Generate CZML file for Cesium visualization.

**Parameters:**
- `tle_path`: Path to TLE file
- `start_dt`: Start datetime
- `duration_s`: Duration (seconds)
- `step_s`: Position sampling step (seconds)
- `out_path`: Output CZML file path
- `accent_rgba`: Accent color [R, G, B, A] (auto-assigned if None)
- `show_beam`: Show satellite beam cones
- `beam_half_angle`: Beam half-angle (degrees)
- `use_fixed_frame`: Use ECEF frame (True) or Inertial frame (False)

**Example:**
```python
from app.czml import write_czml
from datetime import datetime, timezone

write_czml(
    tle_path='InputData/LEO_sat.tle',
    start_dt=datetime.now(timezone.utc),
    duration_s=10800,
    step_s=60,
    out_path='app/web/data/satellites.czml',
    show_beam=True,
    beam_half_angle=21.0,
    use_fixed_frame=False
)
```

---

## Constants

### Physical Constants (core.sim_data)

```python
SPEED_OF_LIGHT = 299792458.0      # m/s
MU_EARTH = 3.986004418e14        # m³/s² (Earth gravitational parameter)
OMEGA_EARTH = 7.2921151467e-5    # rad/s (Earth rotation rate)
F_L1 = 1575.42e6                 # Hz (GPS L1 frequency)
```

### WGS84 Ellipsoid (core.utilities)

```python
WGS84_A = 6378137.0              # m (semi-major axis)
WGS84_F = 1/298.257223563        # flattening
WGS84_E2 = 0.00669437999014      # first eccentricity squared
```
