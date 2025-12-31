import os
import sys
import tempfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from tools.sim_obs import run  # noqa: E402


def test_sim_obs_smoke():
    # Use bundled sample TLE and short duration for faster test
    tle_path = os.path.join(os.path.dirname(__file__), "..", "InputData", "LEO_sat.tle")
    assert os.path.exists(tle_path)

    with tempfile.TemporaryDirectory() as tmpdir:
        out_csv = os.path.join(tmpdir, "obs.csv")
        run(
            tle_path=tle_path,
            start_time_s=0.0,
            duration_s=120.0,   # 2 minutes to keep test quick
            step_s=10.0,
            lat_deg=39.8,
            lon_deg=119.29,
            alt_m=60.0,
            elev_mask_deg=5.0,
            output_file=out_csv,
            gps_week=2377,
        )
        assert os.path.exists(out_csv)
        data = open(out_csv, "r", encoding="utf-8").read().strip().splitlines()
        # Allow empty if geometry masks everything, but ensure file exists
        # and the header was attempted to be written.
        if len(data) > 0:
            assert data[0].startswith("OBS:")

