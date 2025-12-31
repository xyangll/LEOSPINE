# Changelog

All notable changes to SPINE will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

### Added
- (Future changes will be listed here)

---

## [1.0.0] - 2025-12-31

### Added

#### Open-Source Release
- Initial open-source release of LEO-SPINE
- Comprehensive documentation suite:
  - Full API documentation in `docs/api.md`
  - Detailed architecture documentation with diagrams
  - User guide with step-by-step instructions
  - Setup guide for all platforms
- Application logo support (window icon and About dialog)
- About dialog accessible from Settings tab
- Updated `pyproject.toml` with Python 3.13 support and improved dev dependencies

#### Constellation Design
- Walker-Delta constellation TLE generation with configurable P/S/F parameters
- RAAN spread configuration for partial coverage constellations
- 3D Cesium visualization with orbit paths and satellite beams
- Automatic color assignment for each satellite
- Reference frame selection (Inertial/Fixed ECEF)
- Global satellite visibility analysis (instant and average modes)
- PDOP heatmap visualization
- Collapsible satellite list in 3D view

#### Observation Simulation
- Pseudorange and Doppler observation generation
- Error models:
  - Ionosphere delay (Klobuchar model)
  - Troposphere delay (Saastamoinen model)
  - Multipath effect
  - Relativistic effect (corrected formula)
  - Hardware delay
  - Phase center offset (PCO)
- Receiver clock modeling:
  - Model-based (bias + drift)
  - File-based (external data)
- Gaussian noise injection for measurements
- Configurable station coordinates (LLH or ECEF)
- Satellite selection: Choose which satellites to include in simulation

#### Positioning Computation
- Four positioning strategies:
  - Pseudorange only
  - Doppler only
  - PR + Doppler (fixed weight)
  - PR + Doppler (adaptive weight / VCE)
- Error correction models:
  - Ionosphere, Troposphere, Multipath, Relativistic, Hardware delay, PCO
- Satellite selection: Choose which satellites to use for positioning
- Real-time visualization:
  - Skyplot (satellite Az/El)
  - Positioning error plot
  - Pseudorange residuals
  - Doppler residuals
- Convergence analysis with statistics

#### GUI
- Modern PySide6-based interface
- Dark/Light theme support
- Scrollable configuration panels
- Tab-based layout for different functions
- Application logo and About dialog
- Settings tab with theme and Cesium configuration options

#### CLI Tools
- `tools/sim_obs.py` for observation simulation
- `tools/run_positioning.py` for positioning



---

## Types of Changes

- **Added** for new features
- **Changed** for changes in existing functionality
- **Deprecated** for soon-to-be removed features
- **Removed** for now removed features
