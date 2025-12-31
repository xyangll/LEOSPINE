"""Pytest configuration and fixtures for SPINE tests."""
import os
import sys
import pytest

# Ensure project root is in path
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


@pytest.fixture
def sample_tle_path():
    """Return path to sample TLE file."""
    path = os.path.join(ROOT, "InputData", "LEO_sat.tle")
    if os.path.exists(path):
        return path
    return None


@pytest.fixture
def temp_output_dir(tmp_path):
    """Provide a temporary output directory."""
    return str(tmp_path)


@pytest.fixture
def beijing_position():
    """Return Beijing position (lat, lon, alt)."""
    return {
        'lat_deg': 39.9042,
        'lon_deg': 116.4074,
        'alt_m': 50.0
    }
