import math
from datetime import datetime, timezone

MU_EARTH = 3.986004418e14


def _doy_fraction(dt):
    """Calculate day of year with fractional portion."""
    year = dt.year
    start = datetime(year, 1, 1, tzinfo=timezone.utc)
    delta = dt - start
    doy = delta.days + 1
    frac = (delta.seconds + delta.microseconds / 1e6) / 86400.0
    return doy + frac


def _tle_checksum(line):
    """
    Calculate TLE checksum (modulo 10).
    Sum all digits, count '-' as 1, ignore all other characters.
    """
    total = 0
    for ch in line:
        if ch.isdigit():
            total += int(ch)
        elif ch == '-':
            total += 1
    return total % 10


def _format_tle_exp(value):
    """
    Format a value in TLE exponential notation (8 characters total).
    Format: ±NNNNN±N (e.g., " 00000-0", "-12345-6", " 12345+0")
    
    Standard TLE exponential format:
    - Position 1: sign (' ' for positive, '-' for negative)
    - Positions 2-6: 5-digit mantissa (leading decimal assumed, e.g., .12345)
    - Position 7: exponent sign ('+' or '-')
    - Position 8: exponent magnitude (single digit)
    """
    if value == 0 or value is None:
        return " 00000-0"
    
    # Determine sign
    sign = ' ' if value >= 0 else '-'
    abs_val = abs(value)
    
    # Convert to scientific notation to get exponent
    if abs_val != 0:
        exp = math.floor(math.log10(abs_val))
    else:
        exp = 0
    
    # Normalize mantissa to be between 0.1 and 1.0
    # TLE format assumes leading decimal (e.g., .12345 = 0.12345)
    mantissa = abs_val / (10 ** exp)
    
    # Adjust so mantissa is in form 0.NNNNN (between 0.1 and 1)
    if mantissa >= 1:
        mantissa /= 10
        exp += 1
    
    # Format mantissa as 5 digits (remove leading "0.")
    mantissa_int = int(round(mantissa * 100000))
    mantissa_str = f"{mantissa_int:05d}"
    
    # Format exponent (single digit with sign)
    exp_sign = '+' if exp >= 0 else '-'
    exp_digit = abs(exp) % 10  # Single digit
    
    return f"{sign}{mantissa_str}{exp_sign}{exp_digit}"


def _format_tle_first_deriv(value):
    """
    Format first time derivative of mean motion (ndot/2) for TLE Line 1.
    Format: 10 characters total including sign and decimal point.
    Standard format: ±.NNNNNNNN (e.g., " .00000000", "-.00012345")
    Columns 34-43 in TLE Line 1.
    """
    if value == 0 or value is None:
        return " .00000000"
    
    sign = ' ' if value >= 0 else '-'
    # Format: ±.NNNNNNNN (8 decimal places after the decimal point)
    formatted = f"{abs(value):.8f}"
    # Remove leading zero before decimal
    if formatted.startswith('0'):
        formatted = formatted[1:]
    return f"{sign}{formatted}"


def _format_eccentricity(e):
    """
    Format eccentricity for TLE Line 2 (7 characters).
    Leading decimal point is assumed, so we just output 7 digits.
    E.g., 0.0001234 -> "0001234"
    """
    # Multiply by 10^7 and format as 7 digits
    ecc_int = int(round(e * 10000000))
    return f"{ecc_int:07d}"


def generate_walker_delta_tles(name_prefix, a_m, e, inc_deg, raan0_deg, argp_deg, mean_anomaly0_deg, epoch_dt, planes, sats_per_plane, phasing, bstar=0.0, catalog_start=80000, raan_spread=360.0):
    """
    Generate TLE data for Walker-Delta constellation.
    
    Strictly follows standard TLE format (69 characters per line).
    
    Walker-Delta T/P/F notation:
        T = total satellites = planes * sats_per_plane
        P = number of orbital planes
        F = phasing factor (0 to P-1)
    
    Standard Walker-Delta formulas:
        RAAN for plane p: RAAN0 + raan_spread * p / P
        Mean anomaly for sat s in plane p: M0 + 360° * s/S + 360° * p*F / (P*S)
    
    Args:
        raan_spread: Total RAAN angular spread in degrees (default 360° for full coverage).
    """
    # Calculate mean motion in revolutions per day
    n_rad_s = math.sqrt(MU_EARTH / (a_m ** 3))
    n_rev_day = n_rad_s * 86400.0 / (2.0 * math.pi)
    
    # Epoch formatting
    epoch_year = epoch_dt.year % 100
    epoch_doy = _doy_fraction(epoch_dt)
    
    name_lines = []
    total_sats = planes * sats_per_plane
    
    for p in range(planes):
        # RAAN spacing: distributed across the specified raan_spread angle
        raan = (raan0_deg + raan_spread * p / planes) % 360.0
        
        for s in range(sats_per_plane):
            prn = p * sats_per_plane + s + 1
            cat = catalog_start + prn
            
            # Walker-Delta phasing formula (STK compatible):
            in_plane_phase = 360.0 * s / sats_per_plane
            inter_plane_phase = 360.0 * p * phasing / total_sats
            phase = (in_plane_phase + inter_plane_phase) % 360.0
            mean_anomaly = (mean_anomaly0_deg + phase) % 360.0
            
            # Generate TLE Line 1 and Line 2
            line1 = _build_tle_line1(cat, epoch_year, epoch_doy, bstar)
            line2 = _build_tle_line2(cat, inc_deg, raan, e, argp_deg, mean_anomaly, n_rev_day)
            
            name = f"{name_prefix}-{prn:03d}"
            name_lines.append((name, line1, line2))
    
    return name_lines


def _build_tle_line1(catalog_num, epoch_year, epoch_doy, bstar, 
                     classification='U', intl_desig='00001A  ',
                     ndot=0.0, nddot=0.0, ephemeris_type=0, element_set=999):
    """
    Build TLE Line 1 strictly following the standard format.
    
    TLE Line 1 Format (69 characters):
    Col  Description
    ---  -----------
    01   Line number (1)
    02   Space
    03-07  Catalog number (5 digits)
    08   Classification (U/C/S)
    09   Space
    10-17  International designator (8 chars: YYNNNPPP)
    18   Space
    19-20  Epoch year (2 digits)
    21-32  Epoch day of year (DDD.DDDDDDDD, 12 chars)
    33   Space
    34-43  First derivative of mean motion / 2 (10 chars)
    44   Space
    45-52  Second derivative of mean motion / 6 (8 chars)
    53   Space
    54-61  BSTAR drag term (8 chars)
    62   Space
    63   Ephemeris type (single digit)
    64   Space
    65-68  Element set number (4 digits)
    69   Checksum
    """
    # Format each field according to TLE specification
    line_num = '1'
    cat_str = f"{catalog_num:5d}"
    
    # Ensure international designator is 8 characters
    intl_str = f"{intl_desig:<8s}"[:8]
    
    # Epoch: 2-digit year + 12-character day of year (DDD.DDDDDDDD)
    epoch_str = f"{epoch_year:02d}{epoch_doy:012.8f}"
    
    # First derivative of mean motion (10 chars)
    ndot_str = _format_tle_first_deriv(ndot)
    
    # Second derivative of mean motion (8 chars, exponential format)
    nddot_str = _format_tle_exp(nddot)
    
    # BSTAR drag term (8 chars, exponential format)
    bstar_str = _format_tle_exp(bstar)
    
    # Ephemeris type (single digit)
    eph_str = str(ephemeris_type % 10)
    
    # Element set number (4 digits, right-justified with leading spaces)
    elset_str = f"{element_set:4d}"
    
    # Build line without checksum (68 characters)
    # Format: "1 NNNNNC NNNNNNNN YYEEEEEEEEE.EEEE ±.DDDDDDDD ±NNNNN±N ±NNNNN±N E NNNN"
    line_body = (
        f"{line_num} "           # Col 1-2: "1 "
        f"{cat_str}"             # Col 3-7: catalog number
        f"{classification}"      # Col 8: classification
        f" "                     # Col 9: space
        f"{intl_str}"            # Col 10-17: international designator
        f" "                     # Col 18: space
        f"{epoch_str}"           # Col 19-32: epoch (14 chars)
        f" "                     # Col 33: space
        f"{ndot_str}"            # Col 34-43: ndot/2 (10 chars)
        f" "                     # Col 44: space
        f"{nddot_str}"           # Col 45-52: nddot/6 (8 chars)
        f" "                     # Col 53: space
        f"{bstar_str}"           # Col 54-61: BSTAR (8 chars)
        f" "                     # Col 62: space
        f"{eph_str}"             # Col 63: ephemeris type
        f" "                     # Col 64: space
        f"{elset_str}"           # Col 65-68: element set number (4 chars)
    )
    
    # Ensure exactly 68 characters before checksum
    line_body = line_body[:68].ljust(68)
    
    # Calculate and append checksum
    checksum = _tle_checksum(line_body)
    
    return line_body + str(checksum)


def _build_tle_line2(catalog_num, inc_deg, raan_deg, ecc, argp_deg, ma_deg, 
                     mean_motion, rev_num=0):
    """
    Build TLE Line 2 strictly following the standard format.
    
    TLE Line 2 Format (69 characters):
    Col  Description
    ---  -----------
    01   Line number (2)
    02   Space
    03-07  Catalog number (5 digits)
    08   Space
    09-16  Inclination (degrees, 8 chars: NNN.NNNN)
    17   Space
    18-25  RAAN (degrees, 8 chars: NNN.NNNN)
    26   Space
    27-33  Eccentricity (7 digits, decimal point assumed)
    34   Space
    35-42  Argument of perigee (degrees, 8 chars: NNN.NNNN)
    43   Space
    44-51  Mean anomaly (degrees, 8 chars: NNN.NNNN)
    52   Space
    53-63  Mean motion (rev/day, 11 chars: NN.NNNNNNNN)
    64-68  Revolution number (5 digits)
    69   Checksum
    """
    line_num = '2'
    cat_str = f"{catalog_num:5d}"
    
    # Angular values: 8 characters with 4 decimal places (NNN.NNNN format)
    inc_str = f"{inc_deg:8.4f}"
    raan_str = f"{raan_deg:8.4f}"
    argp_str = f"{argp_deg:8.4f}"
    ma_str = f"{ma_deg:8.4f}"
    
    # Eccentricity: 7 digits (decimal point assumed)
    ecc_str = _format_eccentricity(ecc)
    
    # Mean motion: 11 characters with 8 decimal places (NN.NNNNNNNN format)
    mm_str = f"{mean_motion:11.8f}"
    
    # Revolution number: 5 digits
    rev_str = f"{rev_num:5d}"
    
    # Build line without checksum (68 characters)
    line_body = (
        f"{line_num} "           # Col 1-2: "2 "
        f"{cat_str}"             # Col 3-7: catalog number
        f" "                     # Col 8: space
        f"{inc_str}"             # Col 9-16: inclination (8 chars)
        f" "                     # Col 17: space
        f"{raan_str}"            # Col 18-25: RAAN (8 chars)
        f" "                     # Col 26: space
        f"{ecc_str}"             # Col 27-33: eccentricity (7 chars)
        f" "                     # Col 34: space
        f"{argp_str}"            # Col 35-42: argument of perigee (8 chars)
        f" "                     # Col 43: space
        f"{ma_str}"              # Col 44-51: mean anomaly (8 chars)
        f" "                     # Col 52: space
        f"{mm_str}"              # Col 53-63: mean motion (11 chars)
        f"{rev_str}"             # Col 64-68: revolution number (5 chars)
    )
    
    # Ensure exactly 68 characters before checksum
    line_body = line_body[:68].ljust(68)
    
    # Calculate and append checksum
    checksum = _tle_checksum(line_body)
    
    return line_body + str(checksum)


def write_tle_file(path, name_lines):
    """
    Write TLE data to file in standard 3-line format.
    Each satellite entry consists of:
    - Line 0: Satellite name (up to 24 characters)
    - Line 1: TLE line 1 (69 characters)
    - Line 2: TLE line 2 (69 characters)
    """
    with open(path, "w", encoding="utf-8") as f:
        for name, l1, l2 in name_lines:
            # Name line (Line 0): typically up to 24 characters
            f.write(f"{name}\n")
            # Line 1: must be exactly 69 characters
            f.write(f"{l1}\n")
            # Line 2: must be exactly 69 characters
            f.write(f"{l2}\n")


# ============================================================
# Global Visibility and PDOP Analysis
# ============================================================

import numpy as np
from sgp4.api import Satrec, jday
from datetime import timedelta, datetime, timezone

# WGS84 constants
WGS84_A = 6378137.0  # semi-major axis (m)
WGS84_F = 1 / 298.257223563  # flattening
WGS84_E2 = 2 * WGS84_F - WGS84_F ** 2  # eccentricity squared


def geodetic_to_ecef_simple(lat_rad, lon_rad, alt_m=0.0):
    """Convert geodetic coordinates (lat, lon, alt) to ECEF (x, y, z)."""
    sin_lat = math.sin(lat_rad)
    cos_lat = math.cos(lat_rad)
    sin_lon = math.sin(lon_rad)
    cos_lon = math.cos(lon_rad)
    
    N = WGS84_A / math.sqrt(1 - WGS84_E2 * sin_lat ** 2)
    x = (N + alt_m) * cos_lat * cos_lon
    y = (N + alt_m) * cos_lat * sin_lon
    z = (N * (1 - WGS84_E2) + alt_m) * sin_lat
    return np.array([x, y, z])


def compute_elevation_azimuth(sat_ecef, rec_ecef, rec_lat_rad, rec_lon_rad):
    """Compute satellite elevation and azimuth from receiver position.
    
    Returns:
        elevation_deg: Elevation angle in degrees
        azimuth_deg: Azimuth angle in degrees (0=North, 90=East)
    """
    # Line of sight vector
    los = sat_ecef - rec_ecef
    los_norm = np.linalg.norm(los)
    if los_norm < 1e-6:
        return 0.0, 0.0
    
    # Local ENU basis vectors
    sin_lat = math.sin(rec_lat_rad)
    cos_lat = math.cos(rec_lat_rad)
    sin_lon = math.sin(rec_lon_rad)
    cos_lon = math.cos(rec_lon_rad)
    
    # East, North, Up unit vectors
    e_east = np.array([-sin_lon, cos_lon, 0.0])
    e_north = np.array([-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat])
    e_up = np.array([cos_lat * cos_lon, cos_lat * sin_lon, sin_lat])
    
    # Project LOS to ENU
    los_e = np.dot(los, e_east)
    los_n = np.dot(los, e_north)
    los_u = np.dot(los, e_up)
    
    # Elevation
    elevation_rad = math.asin(los_u / los_norm)
    elevation_deg = math.degrees(elevation_rad)
    
    # Azimuth (from North, clockwise)
    azimuth_rad = math.atan2(los_e, los_n)
    azimuth_deg = math.degrees(azimuth_rad) % 360.0
    
    return elevation_deg, azimuth_deg


def tle_to_satrec(line1, line2):
    """Create SGP4 Satrec object from TLE lines."""
    return Satrec.twoline2rv(line1, line2)


def propagate_satrec(satrec, dt):
    """Propagate satellite to given datetime, return ECEF position in meters.
    
    Args:
        satrec: SGP4 Satrec object
        dt: datetime object (timezone-aware or naive)
        
    Returns:
        ECEF position array [x, y, z] in meters, or None if propagation fails
    """
    try:
        # Convert datetime to UTC if timezone-aware
        if dt.tzinfo is not None:
            dt = dt.astimezone(timezone.utc).replace(tzinfo=None)
        
        # Use SGP4's built-in jday() function for accurate Julian day calculation
        # This is the recommended and most accurate method
        year = dt.year
        month = dt.month
        day = dt.day
        hour = dt.hour
        minute = dt.minute
        second = dt.second + dt.microsecond / 1e6
        
        # Calculate Julian day using SGP4's jday() function
        jd, fr = jday(year, month, day, hour, minute, second)
        
        # Check TLE epoch - SGP4 may return NaN if propagation time is too far from TLE epoch
        # Get TLE epoch from satrec object (if available)
        # Note: SGP4 is typically valid for ~2 weeks from epoch, but we allow up to 30 days
        try:
            tle_epoch_jd = satrec.jdsatepoch + satrec.jdsatepochF
            prop_jd = jd + fr
            days_from_epoch = abs(prop_jd - tle_epoch_jd)
            
            # Allow up to 30 days from epoch (more lenient than typical 2 weeks)
            if days_from_epoch > 30.0:
                # Too far from epoch, likely to return NaN
                return None
        except (AttributeError, TypeError):
            # If satrec doesn't have epoch info, proceed anyway
            # This shouldn't happen with valid Satrec objects, but handle gracefully
            pass
        
        # SGP4 propagation
        # Note: SGP4 returns (error_code, position_tuple, velocity_tuple)
        # position_tuple is in km in TEME frame
        e, r, v = satrec.sgp4(jd, fr)
        
        # Check SGP4 error code first
        if e != 0:
            # SGP4 error codes: 
            # 1 = mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
            # 2 = mean motion < 0.0
            # 3 = pert elements, ecc < 0.0  or  ecc > 1.0
            # 4 = semi-latus rectum < 0.0
            # 5 = epoch elements are sub-orbital
            # For visibility analysis, we can try to recover from some errors
            if e == 1:  # Mean elements error - try with small time offset
                jd_adj = jd + 1e-6
                e2, r2, v2 = satrec.sgp4(jd_adj, fr)
                if e2 == 0:
                    r, v = r2, v2
                else:
                    # Log error for debugging (only first few to avoid spam)
                    return None
            else:
                # Log error for debugging
                return None
        
        # r is in km (TEME frame), convert to numpy array
        # SGP4 returns a tuple of floats
        # Check for NaN/Inf before converting (more efficient)
        try:
            # Check each element for NaN using identity comparison (NaN != NaN)
            if any(x != x or not np.isfinite(x) for x in r):
                return None
            
            r_km = np.array(r, dtype=np.float64)
            
            # Double-check after conversion (defensive programming)
            if np.any(np.isnan(r_km)) or np.any(np.isinf(r_km)):
                return None
        except (TypeError, ValueError) as ex:
            # If r is not a valid sequence or contains non-numeric values
            return None
        
        # Check for reasonable position values (satellite should be near Earth)
        # Note: LEO satellites are typically 200-2000km, MEO 2000-35786km, GEO ~35786km
        # But allow wider range to handle edge cases and different orbit types
        r_norm = np.linalg.norm(r_km)
        if r_norm < 100.0 or r_norm > 500000.0:  # Less than 100km (too low) or more than 500000km (too high) is suspicious
            return None
        
        # Additional check: ensure position is not all zeros
        if np.allclose(r_km, 0.0):
            return None
        
        # Convert TEME to ECEF (simplified - just Earth rotation)
        # Days since J2000
        j2000_jd = 2451545.0
        days_since_j2000 = (jd + fr) - j2000_jd
        
        # GMST in radians
        gmst_deg = 280.46061837 + 360.98564736629 * days_since_j2000
        gmst_rad = math.radians(gmst_deg % 360.0)
        
        # Rotation matrix TEME -> ECEF
        cos_g = math.cos(gmst_rad)
        sin_g = math.sin(gmst_rad)
        
        x_ecef = cos_g * r_km[0] + sin_g * r_km[1]
        y_ecef = -sin_g * r_km[0] + cos_g * r_km[1]
        z_ecef = r_km[2]
        
        ecef_pos = np.array([x_ecef, y_ecef, z_ecef]) * 1000.0  # Convert km to m
        
        # Final validation of ECEF position
        if np.any(np.isnan(ecef_pos)) or np.any(np.isinf(ecef_pos)):
            return None
        
        return ecef_pos
        
    except Exception as ex:
        # Catch any unexpected errors
        # Don't print here to avoid spam, error will be logged in calling function
        return None


def compute_pdop(sat_positions, rec_ecef):
    """Compute PDOP from satellite positions and receiver position.
    
    Args:
        sat_positions: List of satellite ECEF positions (Nx3 array)
        rec_ecef: Receiver ECEF position (3,)
        
    Returns:
        pdop: Position Dilution of Precision (or inf if < 4 satellites)
    """
    n_sats = len(sat_positions)
    if n_sats < 4:
        return float('inf')
    
    # Build geometry matrix H
    H = np.zeros((n_sats, 4))
    for i, sat_pos in enumerate(sat_positions):
        los = sat_pos - rec_ecef
        rng = np.linalg.norm(los)
        if rng < 1e-6:
            return float('inf')
        
        # Unit vector from receiver to satellite
        H[i, 0] = -los[0] / rng
        H[i, 1] = -los[1] / rng
        H[i, 2] = -los[2] / rng
        H[i, 3] = 1.0  # Clock bias
    
    try:
        # Compute (H^T H)^-1
        HtH = H.T @ H
        Q = np.linalg.inv(HtH)
        
        # PDOP = sqrt(Q[0,0] + Q[1,1] + Q[2,2])
        pdop = math.sqrt(Q[0, 0] + Q[1, 1] + Q[2, 2])
        return pdop
    except np.linalg.LinAlgError:
        return float('inf')


def compute_global_visibility_pdop(tle_list, epoch_dt, lat_grid, lon_grid, 
                                    elevation_mask=10.0, progress_callback=None):
    """Compute global satellite visibility and PDOP.
    
    Args:
        tle_list: List of (name, line1, line2) tuples
        epoch_dt: Datetime for computation
        lat_grid: Array of latitudes in degrees
        lon_grid: Array of longitudes in degrees
        elevation_mask: Minimum elevation angle in degrees
        progress_callback: Optional callback(percent, message)
        
    Returns:
        visibility_map: 2D array of visible satellite counts (lat x lon)
        pdop_map: 2D array of PDOP values (lat x lon)
    """
    n_lat = len(lat_grid)
    n_lon = len(lon_grid)
    
    visibility_map = np.zeros((n_lat, n_lon), dtype=int)
    pdop_map = np.full((n_lat, n_lon), np.inf)
    
    # Create Satrec objects for all satellites
    satrecs = []
    for name, l1, l2 in tle_list:
        try:
            sat = tle_to_satrec(l1, l2)
            # Check if Satrec object was created successfully
            if sat.error != 0:
                print(f"[Visibility] TLE parse error for {name}: error code {sat.error}")
                continue
            satrecs.append(sat)
        except Exception as e:
            print(f"[Visibility] Failed to parse TLE for {name}: {e}")
    
    print(f"[Visibility] Created {len(satrecs)} Satrec objects from {len(tle_list)} TLEs")
    
    # Verify first Satrec object (for debugging)
    if len(satrecs) > 0:
        try:
            first_sat = satrecs[0]
            first_name = tle_list[0][0] if len(tle_list) > 0 else "Unknown"
            epoch_jd = first_sat.jdsatepoch + first_sat.jdsatepochF
            print(f"[Visibility] First satellite ({first_name}) TLE epoch: JD {epoch_jd:.2f}, propagation time: {epoch_dt}")
        except Exception as e:
            print(f"[Visibility] Warning: Could not check first Satrec: {e}")
    
    if len(satrecs) == 0:
        print("[Visibility] No valid Satrec objects!")
        return visibility_map, pdop_map
    
    # Propagate all satellites to epoch
    sat_positions = []
    failed_sats = []
    error_details = {}  # Store error details for first few failures
    
    for i, sat in enumerate(satrecs):
        # Try propagation with detailed error reporting for first few failures
        pos = propagate_satrec(sat, epoch_dt)
        if pos is not None:
            sat_positions.append(pos)
        else:
            failed_sats.append(i + 1)
            sat_name = tle_list[i][0] if i < len(tle_list) else f"Satellite {i+1}"
            
            # For first 3 failures, try to get detailed error info
            if len(failed_sats) <= 3:
                try:
                    # Try propagation again to capture detailed diagnostics
                    year = epoch_dt.year
                    month = epoch_dt.month
                    day = epoch_dt.day
                    hour = epoch_dt.hour
                    minute = epoch_dt.minute
                    second = epoch_dt.second + epoch_dt.microsecond / 1e6
                    jd, fr = jday(year, month, day, hour, minute, second)
                    
                    # Check TLE epoch
                    try:
                        tle_epoch_jd = sat.jdsatepoch + sat.jdsatepochF
                        prop_jd = jd + fr
                        days_from_epoch = abs(prop_jd - tle_epoch_jd)
                        epoch_info = f"TLE epoch: JD {tle_epoch_jd:.2f}, days from epoch: {days_from_epoch:.1f}"
                    except (AttributeError, TypeError):
                        epoch_info = "TLE epoch: unknown"
                    
                    e, r, v = sat.sgp4(jd, fr)
                    
                    # Check what went wrong
                    if e != 0:
                        error_details[sat_name] = f"SGP4 error code: {e}, {epoch_info}"
                    else:
                        r_km = np.array(r)
                        r_norm = np.linalg.norm(r_km) if not np.any(np.isnan(r_km)) else float('nan')
                        has_nan = np.any(np.isnan(r_km))
                        has_inf = np.any(np.isinf(r_km))
                        range_ok = not has_nan and not has_inf and 100.0 <= r_norm <= 500000.0
                        error_details[sat_name] = f"SGP4 OK (e={e}), r_norm={r_norm:.2f}km, NaN={has_nan}, Inf={has_inf}, RangeOK={range_ok}, {epoch_info}"
                except Exception as ex:
                    error_details[sat_name] = f"Exception: {ex}"
            
            if len(failed_sats) <= 10:  # Only print first 10 to avoid spam
                print(f"[Visibility] Propagation failed for satellite {i+1} ({sat_name})")
    
    if failed_sats:
        print(f"[Visibility] Warning: {len(failed_sats)} satellites failed propagation")
        if error_details:
            print(f"[Visibility] Error details (first {len(error_details)}): {error_details}")
        if len(failed_sats) > 10:
            print(f"[Visibility] ... and {len(failed_sats) - 10} more satellites failed")
    
    print(f"[Visibility] Successfully propagated {len(sat_positions)} satellites")
    
    if len(sat_positions) == 0:
        print("[Visibility] No valid satellite positions!")
        return visibility_map, pdop_map
    
    sat_positions = np.array(sat_positions)
    n_sats = len(sat_positions)
    
    total_points = n_lat * n_lon
    processed = 0
    
    # Compute visibility and PDOP for each grid point
    for i, lat_deg in enumerate(lat_grid):
        lat_rad = math.radians(lat_deg)
        
        for j, lon_deg in enumerate(lon_grid):
            lon_rad = math.radians(lon_deg)
            
            # Receiver ECEF position (on surface)
            rec_ecef = geodetic_to_ecef_simple(lat_rad, lon_rad, 0.0)
            
            # Find visible satellites
            visible_sats = []
            for sat_pos in sat_positions:
                el, az = compute_elevation_azimuth(sat_pos, rec_ecef, lat_rad, lon_rad)
                if el >= elevation_mask:
                    visible_sats.append(sat_pos)
            
            visibility_map[i, j] = len(visible_sats)
            
            # Compute PDOP if enough satellites
            if len(visible_sats) >= 4:
                pdop_map[i, j] = compute_pdop(visible_sats, rec_ecef)
            
            processed += 1
            
            # Progress callback
            if progress_callback and processed % 100 == 0:
                percent = int(100 * processed / total_points)
                progress_callback(percent, f"Computing grid point {processed}/{total_points}")
    
    return visibility_map, pdop_map


def compute_visibility_analysis(sim_or_tles, epoch_dt=None, 
                                 lat_resolution=5.0, lon_resolution=5.0,
                                 elevation_mask=10.0, progress_callback=None):
    """High-level function to compute visibility analysis.
    
    Args:
        sim_or_tles: Either a SatelliteSimulation object or list of (name, l1, l2)
        epoch_dt: Datetime for analysis (default: now)
        lat_resolution: Latitude grid resolution in degrees
        lon_resolution: Longitude grid resolution in degrees
        elevation_mask: Minimum elevation angle in degrees
        progress_callback: Optional callback(percent, message)
        
    Returns:
        dict with keys: 'lat_grid', 'lon_grid', 'visibility', 'pdop', 'epoch'
    """
    from datetime import datetime, timezone
    
    if epoch_dt is None:
        epoch_dt = datetime.now(timezone.utc)
    
    # Extract TLE list
    if hasattr(sim_or_tles, 'tle_satellites'):
        # It's a SatelliteSimulation object
        tle_list = []
        for prn, sat_info in sim_or_tles.tle_satellites.items():
            if hasattr(sat_info, 'line1') and hasattr(sat_info, 'line2'):
                tle_list.append((sat_info.name, sat_info.line1, sat_info.line2))
            elif isinstance(sat_info, dict) and 'line1' in sat_info and 'line2' in sat_info:
                tle_list.append((sat_info.get('name', f'SAT-{prn}'), sat_info['line1'], sat_info['line2']))
    else:
        tle_list = sim_or_tles
    
    # Create grid
    lat_grid = np.arange(-90, 90 + lat_resolution, lat_resolution)
    lon_grid = np.arange(-180, 180 + lon_resolution, lon_resolution)
    
    # Compute
    visibility, pdop = compute_global_visibility_pdop(
        tle_list, epoch_dt, lat_grid, lon_grid, 
        elevation_mask, progress_callback
    )
    
    return {
        'lat_grid': lat_grid,
        'lon_grid': lon_grid,
        'visibility': visibility,
        'pdop': pdop,
        'epoch': epoch_dt,
        'elevation_mask': elevation_mask
    }