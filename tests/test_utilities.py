"""Unit tests for core.utilities module."""
import numpy as np
import pytest
from core.utilities import (
    geodetic_to_ecef,
    ecef_to_geodetic,
)


class TestCoordinateConversions:
    """Test coordinate conversion functions."""

    def test_geodetic_to_ecef_equator(self):
        """Test conversion at equator, prime meridian."""
        lat_rad = 0.0
        lon_rad = 0.0
        alt_m = 0.0
        x, y, z = geodetic_to_ecef(lat_rad, lon_rad, alt_m)
        
        # At equator, prime meridian, x should be ~Earth radius, y and z should be ~0
        assert abs(x - 6378137.0) < 1.0  # WGS84 equatorial radius
        assert abs(y) < 1.0
        assert abs(z) < 1.0

    def test_geodetic_to_ecef_north_pole(self):
        """Test conversion at North Pole."""
        lat_rad = np.pi / 2  # 90 degrees
        lon_rad = 0.0
        alt_m = 0.0
        x, y, z = geodetic_to_ecef(lat_rad, lon_rad, alt_m)
        
        # At North Pole, z should be ~polar radius, x and y should be ~0
        assert abs(x) < 1.0
        assert abs(y) < 1.0
        assert z > 6356000  # WGS84 polar radius is ~6356752

    def test_geodetic_ecef_roundtrip(self):
        """Test roundtrip conversion geodetic -> ECEF -> geodetic."""
        lat_deg = 39.9042  # Beijing
        lon_deg = 116.4074
        alt_m = 100.0
        
        lat_rad = np.deg2rad(lat_deg)
        lon_rad = np.deg2rad(lon_deg)
        
        # Convert to ECEF
        x, y, z = geodetic_to_ecef(lat_rad, lon_rad, alt_m)
        
        # Convert back
        lat2, lon2, alt2 = ecef_to_geodetic(x, y, z)
        
        # Should match within tolerance
        assert abs(lat2 - lat_rad) < 1e-10
        assert abs(lon2 - lon_rad) < 1e-10
        assert abs(alt2 - alt_m) < 0.01  # 1cm accuracy

    def test_geodetic_to_ecef_with_altitude(self):
        """Test conversion with non-zero altitude."""
        lat_rad = np.deg2rad(0.0)  # Equator for simpler test
        lon_rad = np.deg2rad(0.0)
        alt_m = 10000.0  # 10km altitude
        
        x, y, z = geodetic_to_ecef(lat_rad, lon_rad, alt_m)
        
        # At equator, x should be approximately Earth radius + altitude
        # WGS84 equatorial radius is 6378137.0m
        expected_x = 6378137.0 + alt_m
        assert abs(x - expected_x) < 1.0  # Allow 1m tolerance


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_south_pole(self):
        """Test South Pole conversion."""
        lat_rad = -np.pi / 2
        lon_rad = 0.0
        alt_m = 0.0
        x, y, z = geodetic_to_ecef(lat_rad, lon_rad, alt_m)
        
        assert abs(x) < 1.0
        assert abs(y) < 1.0
        assert z < -6356000

    def test_date_line(self):
        """Test conversion at 180 degrees longitude."""
        lat_rad = 0.0
        lon_rad = np.pi  # 180 degrees
        alt_m = 0.0
        x, y, z = geodetic_to_ecef(lat_rad, lon_rad, alt_m)
        
        assert x < -6378000  # Should be negative x
        assert abs(y) < 1.0
        assert abs(z) < 1.0
