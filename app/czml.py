import json
import os
import math
from datetime import datetime, timezone, timedelta
import numpy as np
from core.utilities import datetime_to_gps_time, gps_time_to_datetime
from core.sim_data import SatelliteSimulation, MU_EARTH, datetime_to_jday

# Earth radius in meters
EARTH_RADIUS = 6378137.0

def teme_to_ecef(pos_teme, dt):
    """
    Convert TEME (True Equator Mean Equinox) coordinates to ECEF.
    
    Args:
        pos_teme: Position in TEME [x, y, z] in meters
        dt: datetime object (UTC)
        
    Returns:
        pos_ecef: Position in ECEF [x, y, z] in meters
    """
    # Compute Julian Date
    jd = dt.toordinal() + 1721424.5
    fr = (dt.hour * 3600 + dt.minute * 60 + dt.second + dt.microsecond / 1e6) / 86400.0
    
    # Days since J2000 (2000-01-01 12:00:00 UTC)
    j2000_jd = 2451545.0
    days_since_j2000 = (jd + fr) - j2000_jd
    
    # Greenwich Mean Sidereal Time (GMST) in degrees
    # Using simplified IAU formula
    gmst_deg = 280.46061837 + 360.98564736629 * days_since_j2000
    gmst_rad = math.radians(gmst_deg % 360.0)
    
    # Rotation matrix: TEME -> ECEF (rotate by -GMST around Z axis)
    cos_g = math.cos(gmst_rad)
    sin_g = math.sin(gmst_rad)
    
    x_teme, y_teme, z_teme = pos_teme
    
    # Apply rotation (TEME to ECEF)
    x_ecef = cos_g * x_teme + sin_g * y_teme
    y_ecef = -sin_g * x_teme + cos_g * y_teme
    z_ecef = z_teme
    
    return [x_ecef, y_ecef, z_ecef]

def get_satellite_color(index, total):
    """Generate a distinct color for each satellite using HSV color space."""
    # Use a set of visually distinct colors
    colors = [
        [255, 87, 87, 255],    # Red
        [255, 165, 87, 255],   # Orange
        [255, 230, 87, 255],   # Yellow
        [165, 255, 87, 255],   # Lime
        [87, 255, 87, 255],    # Green
        [87, 255, 165, 255],   # Spring
        [87, 255, 255, 255],   # Cyan
        [87, 165, 255, 255],   # Sky
        [87, 87, 255, 255],    # Blue
        [165, 87, 255, 255],   # Purple
        [255, 87, 255, 255],   # Magenta
        [255, 87, 165, 255],   # Pink
    ]
    return colors[index % len(colors)]

def write_czml(tle_path, start_dt, duration_s, step_s, out_path, accent_rgba=None,
               show_beam=False, beam_half_angle=21.0, use_fixed_frame=False):
    """
    Generate CZML for Cesium visualization.
    
    Performance notes:
    - step_s: position sampling step, recommend 60-120 s
    - Trail length limited to one orbit period or duration_s
    - Labels only visible when near
    
    Args:
        show_beam: If True, display satellite beam cones
        beam_half_angle: Beam half-angle in degrees (default 21° typical for GPS)
        use_fixed_frame: If True, use ECEF (Fixed) frame; if False, use Inertial frame
    """
    sim = SatelliteSimulation()
    sim.load_tle_file(tle_path)
    gps_week, gps_seconds = datetime_to_gps_time(start_dt.replace(tzinfo=None))
    sim.gps_week = gps_week
    prns = list(sim.tle_satellites.keys())
    start_iso = start_dt.astimezone(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
    end_dt = start_dt + timedelta(seconds=duration_s)
    end_iso = end_dt.astimezone(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
    doc = [{"id":"document","version":"1.0","clock":{"interval":f"{start_iso}/{end_iso}","currentTime":start_iso,"multiplier":60}}]
    
    # Performance: increase sampling step (at least 60s)
    effective_step = max(step_s, 60)
    
    total_sats = len(prns)
    
    for idx, prn in enumerate(prns):
        times = []
        t = gps_seconds
        while t <= gps_seconds + duration_s:
            sat = sim.tle_satellites[prn]['satellite']
            dt = gps_time_to_datetime(sim.gps_week, t)
            jd, fr = datetime_to_jday(dt)
            try:
                err, p_km, v_km_s = sat.sgp4(jd, fr)
            except Exception:
                err = 99
                p_km = None
            if err == 0 and p_km is not None:
                # SGP4 outputs TEME coordinates in km, convert to meters
                p_teme_m = [p_km[0]*1000.0, p_km[1]*1000.0, p_km[2]*1000.0]
                if math.isfinite(p_teme_m[0]) and math.isfinite(p_teme_m[1]) and math.isfinite(p_teme_m[2]):
                    if use_fixed_frame:
                        # Convert TEME to ECEF for Earth-fixed display
                        pos = teme_to_ecef(p_teme_m, dt)
                    else:
                        # Use TEME directly for Inertial display (closed orbit)
                        pos = p_teme_m
                    times.append((t - gps_seconds, pos, dt))
            t += effective_step
        if not times:
            # Fallback: use Kepler propagation if SGP4 fails
            e = sim.tle_satellites[prn]
            l1 = e['line1']
            l2 = e['line2']
            try:
                inc_deg = float(l2[8:16])
                raan_deg = float(l2[17:25])
                ecc_str = l2[26:33].strip()
                ecc = float('0.' + ecc_str) if ecc_str else 0.0
                argp_deg = float(l2[34:42])
                M0_deg = float(l2[43:51])
                n_rev_day = float(l2[52:63])
                n_rad_s = n_rev_day * 2.0 * math.pi / 86400.0
                a = (MU_EARTH / (n_rad_s**2)) ** (1.0/3.0)
                inc = math.radians(inc_deg)
                raan = math.radians(raan_deg)
                argp = math.radians(argp_deg)
                M0 = math.radians(M0_deg)
                def kepler_pos(dt_s):
                    M = M0 + n_rad_s * dt_s
                    E = M
                    for _ in range(20):
                        E = M + ecc * math.sin(E)
                    cosE = math.cos(E)
                    sinE = math.sin(E)
                    r = a * (1 - ecc * cosE)
                    x_orb = a * (cosE - ecc)
                    y_orb = a * math.sqrt(1 - ecc**2) * sinE
                    cosw = math.cos(argp)
                    sinw = math.sin(argp)
                    cosO = math.cos(raan)
                    sinO = math.sin(raan)
                    cosi = math.cos(inc)
                    sini = math.sin(inc)
                    x_eci = (cosO*cosw - sinO*sinw*cosi)*x_orb + (-cosO*sinw - sinO*cosw*cosi)*y_orb
                    y_eci = (sinO*cosw + cosO*sinw*cosi)*x_orb + (-sinO*sinw + cosO*cosw*cosi)*y_orb
                    z_eci = (sinw*sini)*x_orb + (cosw*sini)*y_orb
                    return [x_eci, y_eci, z_eci]
                t_offset = 0.0
                while t_offset <= duration_s:
                    dt = start_dt + timedelta(seconds=t_offset)
                    p_eci = kepler_pos(t_offset)
                    if use_fixed_frame:
                        # Convert ECI to ECEF
                        pos = teme_to_ecef(p_eci, dt)
                    else:
                        # Use ECI directly for Inertial display
                        pos = p_eci
                    times.append((t_offset, pos, dt))
                    t_offset += step_s
            except Exception:
                pass
        
        # Build cartesian array for CZML (ECEF coordinates)
        cart = []
        for item in times:
            dt_s = item[0]
            p = item[1]
            cart.extend([dt_s, p[0], p[1], p[2]])
        e = sim.tle_satellites[prn]
        name = e['name']
        
        # Auto-assign different colors to each satellite
        color = get_satellite_color(idx, total_sats)
        
        # Path length: use full simulation duration to show closed orbit
        # LEO orbital period ~90 minutes = 5400 s; recommend duration_s >= 5400
        trail_time = duration_s
        
        ent = {
            "id": f"SAT-{prn}",
            "name": name,
            "availability": f"{start_iso}/{end_iso}",
            # Performance: show labels only when near
            "label": {
                "text": name, 
                "font": "12px sans-serif", 
                "fillColor": {"rgba": [255,255,255,200]},
                "showBackground": False,
                "distanceDisplayCondition": {"distanceDisplayCondition": [0, 50000000]},  # Visible within 50,000 km
                "pixelOffset": {"cartesian2": [0, -10]}
            },
            "position": {"referenceFrame": "FIXED" if use_fixed_frame else "INERTIAL", "epoch": start_iso, "cartesian": cart},
            # Performance: smaller point
            "point": {"pixelSize": 5, "color": {"rgba": color}},
            # Performance: reduce path density
            "path": {
                "show": True, 
                "material": {"solidColor": {"color": {"rgba": [color[0], color[1], color[2], 180]}}},  # Slightly transparent
                "width": 1.5,           # Thin line
                "leadTime": duration_s,   # Show full future path to keep orbit closed
                "trailTime": trail_time,  # Path length (set to duration_s)
                "resolution": 120         # Sample every 120s to reduce vertices
            }
        }
        
        # Add beam cone if enabled
        if show_beam and len(times) > 0:
            # Estimate satellite altitude from first position
            first_pos = times[0][1]
            sat_dist = math.sqrt(first_pos[0]**2 + first_pos[1]**2 + first_pos[2]**2)
            altitude = sat_dist - EARTH_RADIUS
            
            # Beam length extends from satellite toward Earth center
            beam_length = altitude + 100000  # Extend 100km beyond Earth surface for visibility
            
            # Bottom radius based on beam half-angle
            beam_angle_rad = math.radians(beam_half_angle)
            bottom_radius = beam_length * math.tan(beam_angle_rad)
            
            # Beam color (semi-transparent)
            beam_color = [color[0], color[1], color[2], 60]  # Very transparent
            
            ent["cylinder"] = {
                "show": True,
                "length": beam_length,
                "topRadius": 0,  # Cone apex at satellite
                "bottomRadius": bottom_radius,
                "material": {
                    "solidColor": {
                        "color": {"rgba": beam_color}
                    }
                },
                "outline": True,
                "outlineColor": {"rgba": [color[0], color[1], color[2], 120]},
                "outlineWidth": 1,
                "numberOfVerticalLines": 16,
                "slices": 32
            }
        
        doc.append(ent)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w', encoding='utf-8') as f:
        json.dump(doc, f)