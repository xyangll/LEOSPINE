# SPINE Documentation

Welcome to the SPINE documentation!

SPINE (Satellite Positioning & INference Environment) is an open-source tool for satellite constellation design, GNSS observation simulation, and positioning computation.

---

## Documentation Index

### Getting Started

| Document | Description |
|----------|-------------|
| [Setup Guide](setup.md) | Installation and configuration |
| [User Guide](user_guide.md) | How to use each feature |

### Reference

| Document | Description |
|----------|-------------|
| [Architecture](architecture.md) | System design and module structure |
| [API Reference](api.md) | Core module API documentation |

### Project

| Document | Description |
|----------|-------------|
| [README](../README.md) | Project overview |
| [CHANGELOG](../CHANGELOG.md) | Version history |
| [CONTRIBUTING](../CONTRIBUTING.md) | Contribution guidelines |
| [LICENSE](../LICENSE) | MIT License |

---

## Quick Links

- **GitHub Repository**: [LEO-SPINE on GitHub](https://github.com/xyangll/LEOSPINE)
- **Issues**: [Report bugs or request features](https://github.com/xyangll/LEOSPINE/issues)
- **Discussions**: [Community Q&A](https://github.com/xyangll/LEOSPINE/discussions)

---

## Module Overview

```
SPINE/
├── app/           # GUI and visualization
│   ├── gui_qt.py     # Main PySide6 GUI
│   ├── czml.py       # Cesium CZML generation
│   └── web/          # Cesium web assets
│
├── core/          # Core algorithms
│   ├── sim_data.py       # Satellite simulation
│   ├── positioning.py    # Positioning algorithms
│   ├── constellation.py  # Constellation design
│   └── utilities.py      # Coordinate/time utilities
│
├── tools/         # CLI tools
│   ├── sim_obs.py        # Observation simulation
│   └── run_positioning.py # Positioning computation
│
└── docs/          # Documentation (you are here)
```

---

## Features at a Glance

### 🛰️ Constellation Design
- Walker-Delta constellation generation
- 3D Cesium visualization
- Global visibility analysis

### 📡 Observation Simulation
- Pseudorange & Doppler simulation
- Error models (Iono, Tropo, Multipath, etc.)
- Configurable receiver clock

### 📍 Positioning
- Multiple strategies (PR, Doppler, Combined)
- Adaptive weighting (VCE)
- Real-time visualization

---

## Version

Current version: **1.0.0** (see [CHANGELOG](../CHANGELOG.md))

---

## License

SPINE is released under the [MIT License](../LICENSE).
