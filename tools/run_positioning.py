import argparse
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from core.positioning import run_lsq_positioning

def parse_args():
    p = argparse.ArgumentParser(description='GNSS Positioning Tool')
    p.add_argument('--obs', required=True, help='Observation file path')
    p.add_argument('--nav', default=None, help='Navigation file path')
    p.add_argument('--tle', default=None, help='TLE file path')
    p.add_argument('--use_pr', action='store_true', help='Use pseudorange observations')
    p.add_argument('--use_dop', action='store_true', help='Use Doppler observations')
    p.add_argument('--init_x', type=float, default=0.0, help='Initial X position (ECEF)')
    p.add_argument('--init_y', type=float, default=0.0, help='Initial Y position (ECEF)')
    p.add_argument('--init_z', type=float, default=0.0, help='Initial Z position (ECEF)')
    return p.parse_args()

def main():
    """Main entry point for positioning tool."""
    a = parse_args()
    init = [a.init_x, a.init_y, a.init_z]
    run_lsq_positioning(
        a.obs, 
        ephemeris_file=a.nav, 
        tle_file=a.tle, 
        init_pos=init, 
        use_pseudorange=a.use_pr, 
        use_doppler=a.use_dop
    )

if __name__ == '__main__':
    main()