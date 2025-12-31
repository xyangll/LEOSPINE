"""Unit tests for core.constellation module."""
import os
import tempfile
from datetime import datetime, timezone
import pytest
from core.constellation import (
    generate_walker_delta_tles,
    write_tle_file,
)


class TestWalkerDeltaGeneration:
    """Test Walker-Delta constellation generation."""

    def test_generate_basic_constellation(self):
        """Test basic constellation generation."""
        epoch = datetime(2025, 1, 1, 0, 0, 0, tzinfo=timezone.utc)
        tles = generate_walker_delta_tles(
            name_prefix="TEST",
            a_m=7178137.0,  # ~800km altitude
            e=0.001,
            inc_deg=86.4,
            raan0_deg=0.0,
            argp_deg=0.0,
            mean_anomaly0_deg=0.0,
            epoch_dt=epoch,
            planes=3,
            sats_per_plane=4,
            phasing=1,
        )
        
        # Should generate 12 satellites (3 planes × 4 satellites)
        assert len(tles) == 12
        
        # Each TLE should have 3 elements (name, line1, line2)
        for tle in tles:
            assert len(tle) == 3
            name, line1, line2 = tle
            assert isinstance(name, str)
            assert line1.startswith('1 ')
            assert line2.startswith('2 ')

    def test_generate_single_plane(self):
        """Test single-plane constellation."""
        epoch = datetime(2025, 1, 1, 0, 0, 0, tzinfo=timezone.utc)
        tles = generate_walker_delta_tles(
            name_prefix="SINGLE",
            a_m=6978137.0,  # ~600km altitude
            e=0.0,
            inc_deg=45.0,
            raan0_deg=0.0,
            argp_deg=0.0,
            mean_anomaly0_deg=0.0,
            epoch_dt=epoch,
            planes=1,
            sats_per_plane=6,
            phasing=0,
        )
        
        assert len(tles) == 6

    def test_generate_with_raan_spread(self):
        """Test constellation with custom RAAN spread."""
        epoch = datetime(2025, 1, 1, 0, 0, 0, tzinfo=timezone.utc)
        tles = generate_walker_delta_tles(
            name_prefix="SPREAD",
            a_m=7378137.0,  # ~1000km altitude
            e=0.0,
            inc_deg=55.0,
            raan0_deg=0.0,
            argp_deg=0.0,
            mean_anomaly0_deg=0.0,
            epoch_dt=epoch,
            planes=6,
            sats_per_plane=2,
            phasing=1,
            raan_spread=180.0,  # Half coverage
        )
        
        assert len(tles) == 12

    def test_tle_format_validity(self):
        """Test TLE format is valid."""
        epoch = datetime(2025, 1, 1, 0, 0, 0, tzinfo=timezone.utc)
        tles = generate_walker_delta_tles(
            name_prefix="FORMAT",
            a_m=7158137.0,  # ~780km altitude
            e=0.001,
            inc_deg=86.0,
            raan0_deg=0.0,
            argp_deg=0.0,
            mean_anomaly0_deg=0.0,
            epoch_dt=epoch,
            planes=2,
            sats_per_plane=2,
            phasing=0,
        )
        
        for tle in tles:
            name, line1, line2 = tle
            
            # Line 1 checks
            assert len(line1) == 69
            assert line1[0] == '1'
            
            # Line 2 checks
            assert len(line2) == 69
            assert line2[0] == '2'


class TestTLEFileOperations:
    """Test TLE file I/O operations."""

    def test_write_tle_file(self):
        """Test writing TLE file."""
        epoch = datetime(2025, 1, 1, 0, 0, 0, tzinfo=timezone.utc)
        tles = generate_walker_delta_tles(
            name_prefix="WRITE",
            a_m=7178137.0,
            e=0.001,
            inc_deg=86.4,
            raan0_deg=0.0,
            argp_deg=0.0,
            mean_anomaly0_deg=0.0,
            epoch_dt=epoch,
            planes=2,
            sats_per_plane=3,
            phasing=1,
        )
        
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, 'test.tle')
            write_tle_file(filepath, tles)  # Note: path first, then tles
            
            assert os.path.exists(filepath)
            
            # Read back and verify
            with open(filepath, 'r') as f:
                content = f.read()
            
            # Should have content
            assert len(content) > 0
            # Should have 6 satellites × 3 lines each
            lines = content.strip().split('\n')
            assert len(lines) == 18  # 6 satellites × 3 lines


class TestParameterValidation:
    """Test parameter validation."""

    def test_minimum_satellites(self):
        """Test constellation with minimum satellites."""
        epoch = datetime(2025, 1, 1, 0, 0, 0, tzinfo=timezone.utc)
        tles = generate_walker_delta_tles(
            name_prefix="MIN",
            a_m=6878137.0,  # ~500km altitude
            e=0.0,
            inc_deg=0.0,
            raan0_deg=0.0,
            argp_deg=0.0,
            mean_anomaly0_deg=0.0,
            epoch_dt=epoch,
            planes=1,
            sats_per_plane=1,
            phasing=0,
        )
        
        assert len(tles) == 1

    def test_various_inclinations(self):
        """Test various inclination values."""
        epoch = datetime(2025, 1, 1, 0, 0, 0, tzinfo=timezone.utc)
        for inc in [0, 45, 63.4, 86.4, 90, 98.7]:
            tles = generate_walker_delta_tles(
                name_prefix="INC",
                a_m=6978137.0,
                e=0.0,
                inc_deg=inc,
                raan0_deg=0.0,
                argp_deg=0.0,
                mean_anomaly0_deg=0.0,
                epoch_dt=epoch,
                planes=1,
                sats_per_plane=1,
                phasing=0,
            )
            assert len(tles) == 1
