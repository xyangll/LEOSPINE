# SPINE Setup Guide

Complete installation and setup instructions for SPINE.

---

## System Requirements

### Minimum Requirements

| Component | Requirement |
|-----------|-------------|
| OS | Windows 10+, macOS 10.14+, Linux (Ubuntu 18.04+) |
| Python | 3.9 to 3.13 |
| RAM | 4 GB |
| Disk | 500 MB free space |
| Display | 1280×720 resolution |

### Recommended Requirements

| Component | Requirement |
|-----------|-------------|
| Python | 3.10, 3.11, or 3.12 |
| RAM | 8 GB+ |
| Display | 1920×1080 resolution |
| GPU | OpenGL 3.3+ (for Cesium) |

---

## Installation

### 1. Clone Repository

```bash
git clone https://github.com/xyangll/LEOSPINE.git
cd LEOSPINE
```

### 2. Create Virtual Environment (Recommended)

**Windows:**
```bash
python -m venv venv
venv\Scripts\activate
```

**Linux/macOS:**
```bash
python3 -m venv venv
source venv/bin/activate
```

**Conda:**
```bash
conda create -n spine python=3.11
conda activate spine
```

### 3. Install Dependencies

**Option 1: Install from requirements files**
```bash
pip install -r requirements.txt
```

**Option 2: Install as a package (recommended for development)**
```bash
pip install -e .
# Or with development dependencies:
pip install -e ".[dev]"
```

#### Required Packages

| Package | Version | Purpose |
|---------|---------|---------|
| PySide6 | ≥6.6.0 | GUI framework (includes WebEngine) |
| numpy | ≥1.24.0 | Numerical computing |
| scipy | ≥1.10.0 | Scientific computing |
| matplotlib | ≥3.7.0 | Plotting |
| pandas | ≥2.0.0 | Data handling |
| scikit-learn | ≥1.3.0 | Machine learning (DBSCAN) |
| sgp4 | ≥2.22 | Satellite propagation |
| skyfield | ≥1.46 | High-precision astronomy |
| tqdm | ≥4.66.0 | Progress bars |

### 4. Verify Installation

```bash
# Run tests
pytest tests/ -v

# Launch GUI
python main.py
```

---

## Platform-Specific Notes

### Windows

**Qt Platform Plugin Error:**
If you see "qt.qpa.plugin: Could not find the Qt platform plugin":
```bash
pip install --force-reinstall PySide6
```

**WebEngine Issues:**
Ensure you have Visual C++ Redistributable installed.

### macOS

**M1/M2 Macs:**
```bash
# Use Rosetta if needed
arch -x86_64 pip install PySide6
```

**Certificate Errors:**
```bash
pip install --upgrade certifi
```

### Linux

**Missing Qt Libraries:**
```bash
# Ubuntu/Debian
sudo apt install libxcb-xinerama0 libxcb-cursor0

# Fedora
sudo dnf install xcb-util-cursor
```

**Display Issues:**
```bash
export QT_QPA_PLATFORM=xcb
python main.py
```

---

## Development Setup

### Install Development Dependencies

```bash
pip install -r requirements-dev.txt
```

This includes:
- pytest (testing)
- black (formatting)
- flake8 (linting)
- mypy (type checking)

### Running Tests

```bash
# All tests
pytest tests/ -v

# Specific test
pytest tests/test_imports.py -v

# With coverage
pytest tests/ --cov=core --cov-report=html
```

### Code Formatting

```bash
# Format code (line length: 120, configured in pyproject.toml)
black .

# Check style
flake8 core/ app/ tools/
```

---

## Configuration

### Application Settings

Settings are stored in:
- Windows: `%APPDATA%/SPINE/settings.json`
- macOS: `~/Library/Application Support/SPINE/settings.json`
- Linux: `~/.config/SPINE/settings.json`

### Cesium Ion Token (Optional)

For 3D terrain visualization:

1. Create account at [Cesium Ion](https://cesium.com/ion/)
2. Get your access token
3. Enter in Settings → Cesium Ion Token

### Environment Variables

| Variable | Description |
|----------|-------------|
| `SPINE_DATA_DIR` | Override default data directory |
| `QT_SCALE_FACTOR` | HiDPI scaling (e.g., "1.5") |

---

## Troubleshooting

### GUI Won't Start

1. **Check Python version:**
   ```bash
   python --version  # Should be 3.9 to 3.13
   ```

2. **Reinstall Qt:**
   ```bash
   pip uninstall PySide6 PySide6-Essentials PySide6-Addons
   pip install PySide6
   ```

3. **Run with debug:**
   ```bash
   python main.py --debug
   ```

### Import Errors

```bash
# Verify packages
pip list | grep -E "(PySide6|numpy|sgp4)"

# Reinstall all
pip install -r requirements.txt --force-reinstall
```

### Cesium View Issues

1. Check WebEngine is installed:
   ```bash
   python -c "from PySide6.QtWebEngineWidgets import QWebEngineView; print('OK')"
   ```

2. Try external browser:
   - Click "🚀 Open in Browser"

3. Check console for JavaScript errors

### Performance Issues

1. Reduce constellation size
2. Increase CZML step size (60-120s)
3. Use lower visibility grid resolution

---

## Updating

### From Git

```bash
git pull origin main
pip install -r requirements.txt --upgrade
```

### Check for Updates

See [CHANGELOG.md](../CHANGELOG.md) for version history.

---

## Uninstallation

### Remove Virtual Environment

```bash
# Deactivate first
deactivate

# Remove venv folder
rm -rf venv/  # Linux/macOS
rmdir /s venv  # Windows
```

### Remove Settings

- Windows: Delete `%APPDATA%/SPINE/`
- macOS: Delete `~/Library/Application Support/SPINE/`
- Linux: Delete `~/.config/SPINE/`

---

## Getting Help

- **Documentation:** Check `docs/` folder
- **Issues:** [GitHub Issues](https://github.com/xyangll/LEOSPINE/issues)
- **Discussions:** [GitHub Discussions](https://github.com/xyangll/LEOSPINE/discussions)
