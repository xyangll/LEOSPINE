# LEO-SPINE - Low Earth Orbit Satellite Positioning, Initialization & Navigation Estimation

<p align="center">
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-blue.svg" alt="License"></a>
  <a href="https://www.python.org/"><img src="https://img.shields.io/badge/python-3.9+-blue.svg" alt="Python"></a>
  <a href="https://doc.qt.io/qtforpython-6/"><img src="https://img.shields.io/badge/GUI-PySide6-green.svg" alt="PySide6"></a>
</p>

**LEO-SPINE** is an open-source desktop application for LEO satellite constellation design, observation simulation, and positioning computation. It supports Walker-Delta constellation design, TLE-based SGP4 propagation, pseudorange/Doppler observations with multiple error models, and satellite positioning with convergence analysis. 

A modern PySide6 GUI and CLI tools are provided. Cesium 3D visualization is integrated via CZML for browser or embedded WebView. Automated pipelines generate observations, run positioning, and plot convergence, enabling a one-stop platform for LEO navigation enhancement and algorithm validation.

---

## ✨ Features

### 🛰️ Constellation Design
- **Walker-Delta constellation** generation with configurable parameters (P/S/F, RAAN spread)
- **TLE file** generation, loading, and saving
- **3D Cesium visualization** with orbit paths, satellite beams, and reference frame options (inertial/fixed)
- **Global visibility analysis** (instant and average) with satellite count and PDOP heatmaps
- **Automatic color assignment** for distinguishing satellites and orbits

### 📡 Observation Simulation
- **Pseudorange and Doppler** observation simulation
- **Error models**: Ionosphere (Klobuchar), Troposphere (Saastamoinen), Multipath, Relativistic, Hardware delay
- **Receiver clock**: Model-based (bias + drift) or file-based
- **Gaussian noise** injection for pseudorange and Doppler
- **SGP4 propagation** using TLE data
- **Satellite selection**: Choose which satellites to include in simulation

### 📍 Positioning Computation
- **Positioning strategies**: Pseudorange-only, Doppler-only, Combined (fixed/adaptive weight)
- **Error correction models**: Ionosphere, Troposphere, Multipath, Relativistic, Hardware delay, PCO
- **Satellite selection**: Choose which satellites to use for positioning
- **Real-time visualization**: Skyplot, positioning error, residual plots
- **Convergence analysis** with statistics

### 🎨 User Interface
- **Modern PySide6 GUI** with tab-based layout
- **Dark/Light themes** with customizable accent colors
- **Application logo** displayed in window icon and About dialog
- **About dialog** accessible from Settings tab
- **Scrollable configuration panels** for better organization

---

## 🖥️ Screenshots

| Constellation Design | Observation Simulation | Positioning |
|:-------------------:|:---------------------:|:-----------:|
| 3D orbit visualization with Cesium | Error model configuration | Skyplot and error plots |

---

## 🚀 Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/xyangll/LEOSPINE.git
cd LEOSPINE

# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # Linux/Mac
# or: venv\Scripts\activate  # Windows

# Install dependencies
pip install -r requirements.txt

# Or install as a package (recommended for development)
pip install -e .

# Install with development dependencies
pip install -e ".[dev]"
```

### Run GUI

```bash
python main.py
```

### CLI Tools

#### Observation Simulation

```bash
python tools/sim_obs.py --tle InputData/LEO_sat.tle \
  --start 0 --duration 3600 --step 10 \
  --lat 39.8 --lon 119.29 --alt 60 --mask 10 \
  --out OutputData/observations.csv --gps_week 2377
```

**Advanced options:**
```bash
# With error models and noise
python tools/sim_obs.py --tle InputData/LEO_sat.tle \
  --lat 39.8 --lon 119.29 \
  --iono --tropo --multipath \
  --prnoise 5.0 --dopnoise 0.1 \
  --prns "1,3,5,7"  # Select specific satellites
```

#### Positioning

```bash
python tools/run_positioning.py --obs OutputData/observations.csv \
  --tle InputData/LEO_sat.tle \
  --use_pr --use_dop \
  --init_x 0 --init_y 0 --init_z 0
```

---

## 📁 Project Structure

```
LEOSPINE/
├── main.py                 # Application entry point
├── app/
│   ├── gui_qt.py          # Main PySide6 GUI (Constellation, Simulation, Positioning tabs)
│   ├── czml.py            # CZML generation for Cesium 3D visualization
│   ├── settings.py        # Application settings and theming
│   ├── webserver.py       # Local web server for Cesium
│   └── web/               # Cesium web assets
├── core/
│   ├── sim_data.py        # Satellite simulation and observation generation
│   ├── positioning.py     # LSQ positioning algorithms
│   ├── constellation.py   # Walker-Delta TLE generation, visibility analysis
│   └── utilities.py       # Coordinate transforms, time conversions
├── tools/
│   ├── sim_obs.py         # CLI observation simulation
│   └── run_positioning.py # CLI positioning
├── docs/                  # Documentation (includes logo.png)
├── build/                 # Build scripts and outputs
│   ├── build_exe.py       # PyInstaller build script
│   ├── build_simple.bat    # Windows build script
│   ├── build_simple.sh     # Linux/Mac build script
│   └── GUI/               # Generated executable location
├── InputData/             # Sample input files (TLE, clock data)
├── OutputData/            # Generated output files
└── tests/                 # Unit tests
```

---

## 📖 Documentation

| Document | Description |
|----------|-------------|
| [Architecture](docs/architecture.md) | System architecture and module design |
| [Setup Guide](docs/setup.md) | Detailed installation instructions |
| [User Guide](docs/user_guide.md) | How to use each feature |
| [API Reference](docs/api.md) | Core module API documentation |

---

## 🔧 Configuration

### GUI Settings
- **Theme**: Dark/Light mode
- **Accent color**: Customizable
- **Cesium Ion Token**: For 3D terrain (optional)
- **About dialog**: View application information and logo (accessible from Settings tab)

### Simulation Parameters
- Orbit parameters (semi-major axis, eccentricity, inclination, etc.)
- Walker-Delta configuration (planes, satellites per plane, phasing)
- Error model toggles (Ionosphere, Troposphere, Multipath, Relativistic, Hardware)
- Measurement noise levels (pseudorange σ, Doppler σ)
- Receiver clock model (bias, drift) or external file

### Positioning Parameters
- Initial position (LLH or ECEF)
- Positioning strategy selection
- Error correction model toggles
- Satellite selection for positioning

---

## 🧪 Testing

```bash
# Install dev dependencies
pip install -r requirements-dev.txt

# Run tests
pytest tests/ -v

# Run with coverage
pytest tests/ -v --cov=core --cov=app
```

---

## 📦 Building Executable

Build standalone executable for distribution:

```bash
# Windows
cd build
build_simple.bat

# Linux/Mac
cd build
chmod +x build_simple.sh
./build_simple.sh
```

The executable will be created in the `build/GUI/` directory. See [Quick Build Guide](build/QUICK_BUILD.md) for detailed instructions.

---

## 🤝 Contributing

We welcome contributions! Please read our [Contributing Guide](CONTRIBUTING.md) before submitting.

### Development Workflow
1. Fork the repository
2. Create a feature branch: `git checkout -b feature/your-feature`
3. Make changes and test
4. Commit: `git commit -m "feat: add your feature"`
5. Push and create a Pull Request

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

Third-party dependencies and their licenses are listed in [NOTICE](NOTICE).

## 🔒 Security

For security vulnerabilities, please see our [Security Policy](SECURITY.md) and report them privately via [GitHub Security Advisories](https://github.com/xyangll/LEOSPINE/security/advisories/new).

---

## 🙏 Acknowledgments

- [SGP4](https://pypi.org/project/sgp4/) - Satellite propagation
- [Skyfield](https://rhodesmill.org/skyfield/) - High-precision astronomy
- [CesiumJS](https://cesium.com/cesiumjs/) - 3D globe visualization
- [PySide6](https://doc.qt.io/qtforpython-6/) - GUI framework
- [NumPy](https://numpy.org/) / [SciPy](https://scipy.org/) - Numerical computing
- [Matplotlib](https://matplotlib.org/) - Plotting

---

## 📬 Contact

- **Issues**: [GitHub Issues](https://github.com/xyangll/LEOSPINE/issues)
- **Discussions**: [GitHub Discussions](https://github.com/xyangll/LEOSPINE/discussions)

---

<p align="center">
  Made with ❤️ for the LEO navigation community
</p>
