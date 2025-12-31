import numpy as np
import math
from datetime import datetime, timedelta

def ecef_to_geodetic(x, y, z, tol=1e-12, max_iter=100):
    """
    Convert ECEF to WGS-84 geodetic (improved).
    
    Args:
        x, y, z: ECEF (m)
        tol: convergence threshold
        max_iter: max iterations
        
    Returns:
        (lat, lon, h) tuple in (radians, radians, meters)
    """
    # WGS-84 ellipsoid
    a = 6378137.0        # semi-major (m)
    f = 1/298.257223563  # flattening
    b = a*(1 - f)        # semi-minor (m)
    e_sq = f*(2 - f)     # first eccentricity squared
    
    # Auxiliary params
    p = math.hypot(x, y)
    lon = math.atan2(y, x)
    
    # Initial guess
    h = 0.0
    lat = math.atan2(z, p*(1 - e_sq))
    
    # Iterate (Bowring for faster convergence)
    for _ in range(max_iter):
        N = a / math.sqrt(1 - e_sq*math.sin(lat)**2)
        h_prev = h
        h = p / math.cos(lat) - N
        
        # Update latitude estimate
        sin_lat = z / (N*(1 - e_sq) + h)
        lat_new = math.atan((z + e_sq*N*sin_lat)/p)
        
        # Convergence check
        if abs(lat_new - lat) < tol and abs(h - h_prev) < tol:
            break
        lat = lat_new
    
    return lat, lon, h

def geodetic_to_ecef(lat, lon, h):
    """
    Convert WGS-84 geodetic to ECEF (improved).
    
    Args:
        lat: latitude (rad)
        lon: longitude (rad)
        h: height (m)
        
    Returns:
        (x, y, z) in meters
    """
    # WGS-84 ellipsoid
    a = 6378137.0        # semi-major (m)
    f = 1/298.257223563  # flattening
    e_sq = f*(2 - f)     # first eccentricity squared
    
    # Prime vertical radius
    N = a / math.sqrt(1 - e_sq*math.sin(lat)**2)
    
    x = (N + h) * math.cos(lat) * math.cos(lon)
    y = (N + h) * math.cos(lat) * math.sin(lon)
    z = (N*(1 - e_sq) + h) * math.sin(lat)
    
    return x, y, z

def gps_time_to_datetime(gps_week, gps_seconds):
    """
    Convert GPS week + seconds to datetime (UTC).
    
    Args:
        gps_week: GPS week number
        gps_seconds: seconds of week
        
    Returns:
        dt: datetime object (UTC)
    """
    # GPS epoch: 1980-01-06 00:00:00 UTC
    gps_epoch = datetime(1980, 1, 6, 0, 0, 0)
    
    # Elapsed since epoch
    dt = gps_epoch + timedelta(weeks=gps_week, seconds=gps_seconds)
    
    return dt

def datetime_to_gps_time(dt):
    """
    Convert datetime to GPS week/seconds.
    
    Args:
        dt: datetime object
        
    Returns:
        gps_week: GPS week number
        gps_seconds: seconds of GPS week
    """
    # GPS epoch: 1980-01-06 00:00:00 UTC
    gps_epoch = datetime(1980, 1, 6, 0, 0, 0)
    
    # Use total_seconds() to keep microseconds
    total_seconds = (dt - gps_epoch).total_seconds()
    seconds_in_week = 604800.0
    gps_week = int(total_seconds // seconds_in_week)
    gps_seconds = total_seconds - gps_week * seconds_in_week
    return gps_week, gps_seconds

def datetimeUTC_to_gps_time(dt, leap_seconds=18):
    """
    Converts a datetime object to GPS Week and Seconds.
    
    Args:
        dt (datetime): The UTC datetime object.
        leap_seconds (int): Offset between GPS and UTC (currently ~18s).
        
    Returns:
        tuple: (gps_week, gps_seconds)
    """
    # GPS Epoch: 1980-01-06 00:00:00 UTC
    gps_epoch = datetime(1980, 1, 6, 0, 0, 0)
    
    # Calculate time difference
    # Note: Using total_seconds() handles microseconds correctly, 
    # which manual calculation might miss.
    time_diff = dt - gps_epoch
    
    # Add leap seconds because GPS time is continuous and ahead of UTC
    total_seconds = time_diff.total_seconds() + leap_seconds
    
    seconds_in_week = 604800  # 7 * 24 * 60 * 60
    
    gps_week = int(total_seconds // seconds_in_week)
    gps_seconds = total_seconds % seconds_in_week
    
    return gps_week, gps_seconds
