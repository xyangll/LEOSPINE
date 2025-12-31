import argparse
import csv
import numpy as np
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from core.utilities import geodetic_to_ecef
from core.sim_data import SatelliteSimulation, main_sim, read_rck_data

def run(tle_path, start_time_s, duration_s, step_s, lat_deg, lon_deg, alt_m, elev_mask_deg, output_file, gps_week=None,
        prnoise=None, dopnoise=None,
        enable_iono=True, enable_tropo=True, enable_multipath=False, 
        enable_relativistic=False, enable_hardware=False,
        use_clk_model=True, clk_bias=1000.0, clk_drift=0.05,
        clk_file=None, clk_drift_file=None,
        selected_prns=None):
    """
    Run simulation observation generation.
    
    Args:
        enable_iono: Enable ionosphere delay
        enable_tropo: Enable troposphere delay
        enable_multipath: Enable multipath effect
        enable_relativistic: Enable relativistic effect
        enable_hardware: Enable hardware delay
        use_clk_model: True to use bias/drift values, False to use files
        clk_bias: Clock bias in meters (used when use_clk_model=True)
        clk_drift: Clock drift in m/s (used when use_clk_model=True)
        clk_file: Clock bias file path (used when use_clk_model=False)
        clk_drift_file: Clock drift file path (used when use_clk_model=False)
        selected_prns: List of PRNs to simulate (None = all satellites)
    """
    sim = SatelliteSimulation(enable_iono=enable_iono, enable_tropo=enable_tropo, enable_multipath=enable_multipath, enable_relativistic=enable_relativistic, enable_hardware=enable_hardware)
    
    # Handle clock data
    if not use_clk_model and clk_file is not None:
        # Use clock data from files
        rck_values, rck_drift_values = read_rck_data(clk_file, clk_drift_file)
    else:
        # Use model values
        rck_values = None
        rck_drift_values = None

    end_time_s = start_time_s + duration_s
    receiver_pos = np.array(geodetic_to_ecef(np.deg2rad(lat_deg), np.deg2rad(lon_deg), alt_m))
    main_sim(sim, gps_week, start_time_s, end_time_s, step_s, elev_mask_deg, output_file, 
             tle_file=tle_path, 
             clk_data=rck_values, clk_drift_data=rck_drift_values,
             receiver_pos=receiver_pos, prnoise=prnoise, dopnoise=dopnoise,
             clk_bias=clk_bias, clk_drift=clk_drift,
             selected_prns=selected_prns)

def parse_args():
    p = argparse.ArgumentParser(description='GNSS Observation Simulation Tool')
    
    # Basic parameters
    p.add_argument('--tle', required=True, help='TLE file path')
    p.add_argument('--start', type=float, default=0.0, help='Start time in GPS seconds')
    p.add_argument('--duration', type=float, default=3600.0, help='Simulation duration in seconds')
    p.add_argument('--step', type=float, default=10.0, help='Time step in seconds')
    p.add_argument('--lat', type=float, required=True, help='Receiver latitude in degrees')
    p.add_argument('--lon', type=float, required=True, help='Receiver longitude in degrees')
    p.add_argument('--alt', type=float, default=0.0, help='Receiver altitude in meters')
    p.add_argument('--mask', type=float, default=10.0, help='Elevation mask in degrees')
    p.add_argument('--out', type=str, default='observations.csv', help='Output file path')
    p.add_argument('--gps_week', type=int, default=None, help='GPS week number')
    
    # Noise parameters
    p.add_argument('--prnoise', type=float, default=None, help='Pseudorange noise sigma (meters)')
    p.add_argument('--dopnoise', type=float, default=None, help='Doppler noise sigma (m/s)')
    
    # Error model flags
    p.add_argument('--iono', action='store_true', default=True, help='Enable ionosphere delay (default: True)')
    p.add_argument('--no-iono', dest='iono', action='store_false', help='Disable ionosphere delay')
    p.add_argument('--tropo', action='store_true', default=True, help='Enable troposphere delay (default: True)')
    p.add_argument('--no-tropo', dest='tropo', action='store_false', help='Disable troposphere delay')
    p.add_argument('--multipath', action='store_true', default=False, help='Enable multipath effect')
    p.add_argument('--relativistic', action='store_true', default=False, help='Enable relativistic effect')
    p.add_argument('--hardware', action='store_true', default=False, help='Enable hardware delay')
    
    # Clock model parameters
    p.add_argument('--clk-model', action='store_true', default=True, help='Use clock bias/drift model (default)')
    p.add_argument('--clk-file', type=str, default=None, help='Clock bias file (disables clk-model)')
    p.add_argument('--clk-drift-file', type=str, default=None, help='Clock drift file')
    p.add_argument('--clk-bias', type=float, default=1000.0, help='Clock bias in meters (for model)')
    p.add_argument('--clk-drift', type=float, default=0.05, help='Clock drift in m/s (for model)')
    
    # Satellite selection
    p.add_argument('--prns', type=str, default=None, 
                   help='Comma-separated list of PRNs to simulate (e.g., "1,3,5,7"). If not specified, all satellites are used.')
    
    return p.parse_args()

def main():
    """Main entry point for observation simulation tool."""
    a = parse_args()
    
    # Parse PRN selection
    selected_prns = None
    if a.prns:
        selected_prns = [p.strip() for p in a.prns.split(',')]
    
    # Determine if using clock model or files
    use_clk_model = True
    if a.clk_file is not None:
        use_clk_model = False
    
    run(
        tle_path=a.tle,
        start_time_s=a.start,
        duration_s=a.duration,
        step_s=a.step,
        lat_deg=a.lat,
        lon_deg=a.lon,
        alt_m=a.alt,
        elev_mask_deg=a.mask,
        output_file=a.out,
        gps_week=a.gps_week,
        prnoise=a.prnoise,
        dopnoise=a.dopnoise,
        enable_iono=a.iono,
        enable_tropo=a.tropo,
        enable_multipath=a.multipath,
        enable_relativistic=a.relativistic,
        enable_hardware=a.hardware,
        use_clk_model=use_clk_model,
        clk_bias=a.clk_bias,
        clk_drift=a.clk_drift,
        clk_file=a.clk_file,
        clk_drift_file=a.clk_drift_file,
        selected_prns=selected_prns
    )

if __name__ == '__main__':
    main()