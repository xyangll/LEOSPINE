import numpy as np
import math
import matplotlib.pyplot as plt
from core.utilities import ecef_to_geodetic, geodetic_to_ecef
from matplotlib import cm
import os
from datetime import datetime, timedelta
from sgp4.api import Satrec
from sgp4.earth_gravity import wgs84
from skyfield.api import load, EarthSatellite
from skyfield.sgp4lib import TEME_to_ITRF


# Constants
SPEED_OF_LIGHT = 299792458.0  # Speed of light (m/s)
F_L1 = 1.199169832000000e+09   # L1 frequency (Hz)
LAMBDA_L1 = SPEED_OF_LIGHT / F_L1  # L1 wavelength (m)
MU = 3.986005e14  # Earth GM (m^3/s^2)
OMEGA_EARTH = 7.2921151467e-5  # Earth rotation rate (rad/s)
MU_EARTH = 3.986004418000002e+14 # Earth GM (m^3/s^2)

def read_tle_file(filename):
    """
    Read TLE data (multi-satellite) from file.
    
    Args:
        filename: path to TLE file
        
    Returns:
        satellites: list of tuples (name, line1, line2)
    """
    satellites = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        
    i = 0
    while i < len(lines):
        # Skip empty lines
        if not lines[i].strip():
            i += 1
            continue
        
        # Read three lines (satellite name and two TLE lines)
        if i + 2 < len(lines):
            name = lines[i].strip()
            line1 = lines[i+1].strip()
            line2 = lines[i+2].strip()
            
            # Check if TLE lines are valid
            if line1.startswith('1 ') and line2.startswith('2 '):
                satellites.append((name, line1, line2))
            
            i += 3
        else:
            break
    
    return satellites

def gps_time_to_datetime(gps_week, gps_seconds):
    """
    Convert GPS (week, seconds) to datetime.
    
    Args:
        gps_week: GPS week number
        gps_seconds: seconds of GPS week
        
    Returns:
        datetime object in UTC
    """
    # GPS epoch: 1980-01-06 00:00:00
    gps_epoch = datetime(1980, 1, 6)
    
    # Total seconds since GPS epoch
    total_seconds = gps_week * 604800 + gps_seconds
    
    # Current time
    return gps_epoch + timedelta(seconds=total_seconds)

def datetime_to_jday(dt):
    """
    Convert datetime to Julian day.
    
    Args:
        dt: datetime object
        
    Returns:
        (jd, fr): Julian day and fractional day
    """
    year = dt.year
    month = dt.month
    day = dt.day
    hour = dt.hour
    minute = dt.minute
    second = dt.second + dt.microsecond / 1e6
    
    # Integer part of JD
    if month <= 2:
        year -= 1
        month += 12
    
    A = year // 100
    B = 2 - A + A // 4
    
    jd = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day + B - 1524.5
    
    # Fractional day
    fr = (hour + minute / 60.0 + second / 3600.0) / 24.0
    
    return jd, fr

def calculate_elevation_angle(sat_pos, rec_pos):
    """
    Compute satellite elevation/azimuth relative to receiver.
    Args:
        sat_pos: satellite position (ECEF)
        rec_pos: receiver position (ECEF)
    Returns:
        elevation_angle: elevation (degrees)
        azimuth_angle: azimuth (degrees)
    """
    los_vector = sat_pos - rec_pos
    los_norm = np.linalg.norm(los_vector)
    
    # Receiver geodetic
    rec_lat, rec_lon, rec_alt = ecef_to_geodetic(*rec_pos)
    
    # Local basis vectors
    north_vector = np.array([
        -np.sin(rec_lat) * np.cos(rec_lon),
        -np.sin(rec_lat) * np.sin(rec_lon),
        np.cos(rec_lat)
    ])
    east_vector = np.array([
        -np.sin(rec_lon),
        np.cos(rec_lon),
        0
    ])
    up_vector = np.array([
        np.cos(rec_lat) * np.cos(rec_lon),
        np.cos(rec_lat) * np.sin(rec_lon),
        np.sin(rec_lat)
    ])
    
    # Project to local level frame
    los_local = np.array([
        np.dot(los_vector, east_vector),
        np.dot(los_vector, north_vector),
        np.dot(los_vector, up_vector)
    ])
    
    # Elevation
    elevation_angle = np.arcsin(los_local[2]/los_norm) * (180.0 / np.pi)
    
    # Azimuth (0-360 deg)
    azimuth_angle = np.arctan2(los_local[0], los_local[1]) * (180.0 / np.pi)
    azimuth_angle = azimuth_angle % 360  # normalize to 0-360
    
    return elevation_angle, azimuth_angle

def teme_to_ecef(teme_position, teme_velocity, dt):
    # Days since J2000 epoch
    j2000_epoch = datetime(2000, 1, 1, 12, 0, 0)
    days_since_j2000 = (dt - j2000_epoch).total_seconds() / 86400.0
    
    # GMST in radians (simplified; higher accuracy may be needed)
    gmst = 280.46061837 + 360.98564736629 * days_since_j2000
    gmst = gmst % 360.0  # wrap to 0-360
    gmst_rad = np.radians(gmst)
    
    # TEME -> ECEF rotation (simplified, only Earth rotation)
    # For high precision, add precession/nutation, etc.
    rotation_matrix = np.array([
        [np.cos(gmst_rad), np.sin(gmst_rad), 0],
        [-np.sin(gmst_rad), np.cos(gmst_rad), 0],
        [0, 0, 1]
    ])
    
    # Position
    ecef_position = np.dot(rotation_matrix, teme_position)
    
    # Velocity (consider Earth rotation)
    earth_rotation_rate = 7.2921150e-5  # rad/s
    omega_earth = np.array([0, 0, earth_rotation_rate])
    velocity_rotation = np.dot(rotation_matrix, teme_velocity)
    velocity_earth_rotation = np.cross(omega_earth, ecef_position)
    ecef_velocity = velocity_rotation - velocity_earth_rotation
    
    return ecef_position, ecef_velocity

class SatelliteSimulation:
    def __init__(self, enable_iono=True, enable_tropo=True, enable_multipath=False, enable_relativistic=False, enable_hardware=False, enable_pco=False):
        """
        Initialize satellite simulation system
        
        Args:
            enable_iono: Enable ionospheric delay
            enable_tropo: Enable tropospheric delay
            enable_multipath: Enable multipath effect
        """
        self.enable_iono = enable_iono
        self.enable_tropo = enable_tropo
        self.enable_multipath = enable_multipath
        self.enable_relativistic = enable_relativistic
        self.enable_hardware = enable_hardware
        self.enable_pco = enable_pco
        self.tle_satellites = {}  # Store TLE satellite objects
        self.gps_week = 2200  # Default GPS week, adjust as needed
        self.clock0 = 0
        self.clockdrift=0
        self.clock = 0
        self.start_time = 0


        self.ionpara = [
            3.2596E-08, 1.4901E-08, -1.7881E-07, -5.9605E-08,
            1.3722E+05, 1.6384E+04, -3.2768E+05, 3.2768E+05
        ]
    
    def add_tle_satellite(self, name, line1, line2, prn):
        """
        Add TLE satellite
        
        Args:
            name: Satellite name
            line1: TLE line 1
            line2: TLE line 2
            prn: PRN number assigned to satellite (for simulation)
        """
        try:
            # Create SGP4 satellite object
            satellite = Satrec.twoline2rv(line1, line2)
            self.tle_satellites[prn] = {
                'name': name,
                'line1': line1,
                'line2': line2,
                'satellite': satellite,
                'prn': prn
            }
            return True
        except Exception as e:
            print(f"Failed to add TLE satellite: {e}")
            return False
    
    def load_tle_file(self, filename):
        """
        Load TLE data from file
        
        Args:
            filename: TLE file path
        
        Returns:
            loaded_count: Number of successfully loaded satellites
        """
        satellites = read_tle_file(filename)
        loaded_count = 0
        
        for i, (name, line1, line2) in enumerate(satellites):
            # Assign PRN number to each satellite, starting from 1
            prn = i + 1
            if self.add_tle_satellite(name, line1, line2, prn):
                loaded_count += 1
        
        print(f"Successfully loaded {loaded_count} satellites from TLE data")
        return loaded_count
    
    def calculate_satellite_position_tle(self, prn, gps_time):
        """
        Calculating satellite position and velocity using TLE data
        
        Args:
            prn: Satellite PRN number
            gps_time: GPS time (seconds)
            
        Returns:
            sat_pos: Satellite ECEF position [x,y,z] (m)
            sat_vel: Satellite ECEF velocity [vx,vy,vz] (m/s)
        """
        if prn not in self.tle_satellites:
            raise ValueError(f"Satellite with PRN={prn} does not exist in TLE data")
        
        # Get satellite object
        satellite = self.tle_satellites[prn]['satellite']
        
        # Convert GPS time to datetime
        dt = gps_time_to_datetime(self.gps_week, gps_time)
        
        # Convert to Julian day
        jd, fr = datetime_to_jday(dt)
        
        # Calculating satellite position and velocity using SGP4
        error_code, position, velocity = satellite.sgp4(jd, fr)
        
        if error_code != 0:
            print(f"SGP4 calculation error, error code: {error_code}, using zero vector instead")
            return np.zeros(3), np.zeros(3)
        
        # Convert to ITRF (ECEF) coordinates
        # ecef_position, ecef_velocity = TEME_to_ITRF(jd, np.array(position), np.array(velocity))
        ecef_position, ecef_velocity = teme_to_ecef(position, velocity, dt)

        # Convertunit: km到m
        sat_pos = np.array(ecef_position) * 1000.0
        sat_vel = np.array(ecef_velocity) * 1000.0
        
        return sat_pos, sat_vel
        
    def calculate_satellite_position(self, ephemeris, gps_time):
        """
        according toephemerisCalculatesatelliteposition和velocity
        
        Args:
            ephemeris: satelliteephemeris字典
            gps_time: GPS time(秒)
            
        Returns:
        
            sat_pos: satelliteECEFposition[x,y,z](m)
            sat_vel: satelliteECEFvelocity[vx,vy,vz](m/s)
        """
        # CheckwhetherusingTLE模式
        if 'use_tle' in ephemeris and ephemeris['use_tle']:
            return self.calculate_satellite_position_tle(ephemeris['prn'], gps_time)
        
        # usingephemerisArgsCalculatesatelliteposition
        """Calculatesatelliteposition和velocity"""
        # 常量定义
        F = -4.442807633e-10  # 相对论效应常数
        
        # 从ephemeris中提取Args
        toc = ephemeris['toc']
        toe = ephemeris['toe']
        sqrtA = ephemeris['sqrta']
        e = ephemeris['es']
        i0 = ephemeris['inc0']
        Omega0 = ephemeris['omega0']
        omega = ephemeris['omega']
        M0 = ephemeris['m0']
        delta_n = ephemeris['dotn']
        i_dot = ephemeris['idot']
        Omega_dot = ephemeris['omegadot']
        Cuc = ephemeris['cuc']
        Cus = ephemeris['cus']
        Crc = ephemeris['crc']
        Crs = ephemeris['crs']
        Cic = ephemeris['cic']
        Cis = ephemeris['cis']
        af0 = ephemeris['af0']
        af1 = ephemeris['af1']
        af2 = ephemeris['af2']
        tgd1 = ephemeris['tgd1']
        
        # Calculate时间差
        dt_toc = gps_time - toc
        if dt_toc > 302400.0:
            dt_toc = dt_toc - 604800.0
        elif dt_toc < -302400.0:
            dt_toc = dt_toc + 604800.0
        
        # Calculatesatelliteclock bias
        clock_bias = af0 + af1 * dt_toc + af2 * dt_toc * dt_toc - tgd1
        
        # Calculateconsideringclock bias后真实时间
        t_corrected = gps_time - clock_bias
        
        # Calculate时间差 dt
        dt = t_corrected - toe
        if dt > 302400.0:
            dt = dt - 604800.0
        elif dt < -302400.0:
            dt = dt + 604800.0
        
        # Calculateaverage角velocity
        A = sqrtA * sqrtA
        n0 = math.sqrt(MU_EARTH / (A * A * A))
        n = n0 + delta_n
        
        # Calculate平近点角
        M = M0 + n * dt
        
        # iterationCalculate偏近点角 E (correction: using与MATLAB相同牛顿-拉夫森iteration)
        E = M
        max_iter = 30
        iter_count = 0
        diff = 1.0  # 初始设置一大于阈值值
        
        while abs(diff) > 1.0e-12 and iter_count < max_iter:
            iter_count += 1
            diff = (M - (E - e * math.sin(E))) / (1.0 - e * math.cos(E))
            E = E + diff
        
        # Calculate偏近点角正余弦
        cos_E = math.cos(E)
        sin_E = math.sin(E)
        
        # consideringrelativistic效应clock biascorrection
        clock_bias = clock_bias + F * e * sqrtA * sin_E
        clock_drift = af1 + 2 * af2 * dt_toc
        
        # Calculate真近点角 (correction: usingatan2函数而非两步Calculate)
        ta = math.atan2(math.sqrt(1.0 - e * e) * sin_E, cos_E - e)
        
        # Calculate升交角距
        phi = ta + omega
        
        # 二次谐波correction
        sin_2phi = math.sin(2.0 * phi)
        cos_2phi = math.cos(2.0 * phi)
        
        # correction升交角距
        delta_u = Cus * sin_2phi + Cuc * cos_2phi
        u = phi + delta_u
        
        # correction半径
        delta_r = Crs * sin_2phi + Crc * cos_2phi
        r = A * (1.0 - e * cos_E) + delta_r
        
        # correctionorbit倾角
        delta_i = Cis * sin_2phi + Cic * cos_2phi
        i = i0 + delta_i + i_dot * dt
        
        # Calculate在orbit平面内position
        x_orb = r * math.cos(u)
        y_orb = r * math.sin(u)
        
        # correction升交点赤经 (完全按照MATLAB代码实现)
        Omega = Omega0 + (Omega_dot - OMEGA_EARTH) * dt - OMEGA_EARTH * toe
        
        # CalculateECEF坐标
        x = x_orb * math.cos(Omega) - y_orb * math.cos(i) * math.sin(Omega)
        y = x_orb * math.sin(Omega) + y_orb * math.cos(i) * math.cos(Omega)
        z = y_orb * math.sin(i)
        
        # CalculateEderivative
        E_dot = n / (1.0 - e * cos_E)
        
        # Calculateuderivative (correction: 正确Calculatephi_dot)
        phi_dot = math.sqrt(1.0 - e * e) * E_dot / (1.0 - e * cos_E)
        u_dot = phi_dot * (1.0 + 2.0 * (Cus * cos_2phi - Cuc * sin_2phi))
        
        # Calculaterderivative
        r_dot = A * e * sin_E * E_dot + 2.0 * (Crs * cos_2phi - Crc * sin_2phi) * phi_dot
        
        # Calculateiderivative
        i_dot_corrected = i_dot + 2.0 * (Cis * cos_2phi - Cic * sin_2phi) * phi_dot
        
        # CalculateOmegaderivative
        Omega_dot_corrected = Omega_dot - OMEGA_EARTH
        
        # Calculateorbit平面内velocity
        x_orb_dot = r_dot * math.cos(u) - r * math.sin(u) * u_dot
        y_orb_dot = r_dot * math.sin(u) + r * math.cos(u) * u_dot
        
        # CalculateECEF坐标系中velocity (完全按照MATLAB代码实现)
        vx = x_orb_dot * math.cos(Omega) - y_orb_dot * math.cos(i) * math.sin(Omega) + z * math.sin(Omega) * i_dot_corrected - y * Omega_dot_corrected
        
        vy = x_orb_dot * math.sin(Omega) + y_orb_dot * math.cos(i) * math.cos(Omega) - z * math.cos(Omega) * i_dot_corrected + x * Omega_dot_corrected
        
        vz = y_orb_dot * math.sin(i) + y_orb * math.cos(i) * i_dot_corrected
        
        return np.array([x, y, z]), np.array([vx, vy, vz])
    
    def cal_sat_posvel(self, time, prn, ephemeris, rec_pos):
        """
        Calculatesatelliteposition和velocity
        
        Args:
            time: 观测时间column表
            prn: satellitePRN号column表
            ephemeris: ephemeris字典
            rec_pos: receiverposition(ECEF坐标)
            
        Returns:
            sat_positions: satellitepositioncolumn表
            sat_velocities: satellitevelocitycolumn表
        """
        # CheckwhetherusingTLE模式
        if (isinstance(ephemeris, dict) and 'use_tle' in ephemeris and ephemeris['use_tle']) or prn in self.tle_satellites:
            return self.calculate_satellite_position_tle(prn, time)
        
        # Calculateeach观测时刻satelliteposition和velocity（iteration法considering传播时延）
        sat_positions = []
        sat_velocities = []
        c = SPEED_OF_LIGHT  # 光速(m/s)

        current_prn = prn  
        # Get该PRNsatelliteephemeris
        if ephemeris and current_prn in ephemeris:
            # 找到最接近current时间toc
            prn_ephs = ephemeris[current_prn]
            closest_toc = min(prn_ephs.keys(), key=lambda x: abs(x - time))
            current_eph = prn_ephs[closest_toc]
            
            # 添加必要Args
            current_eph['toe'] = closest_toc
            current_eph['prn'] = current_prn
        else:
            current_eph = ephemeris
            
        emission_time = time
        
        # usingfinaliteration发射时间Calculatesatelliteposition和velocity
        sat_positions, sat_velocities = self.calculate_satellite_position(current_eph, emission_time)

        return sat_positions, sat_velocities

    def calculate_geometric_range(self, sat_pos, rec_pos):
        """
        Calculatesatellite与receiver之间geometric distance
        
        Args:
            sat_pos: satelliteECEFposition[x,y,z](m)
            rec_pos: receiverECEFposition[x,y,z](m)
            
        Returns:
            range: geometric distance(m)
            los: line of sightvector(unitvector)
        """
        diff = rec_pos - sat_pos
        range_value = np.linalg.norm(diff)
        los = diff / range_value
        return range_value, los
    
    def calculate_sagnac_correction(self, sat_pos, rec_pos):
        """
        CalculateSagnac effectcorrection量（修复版）
        
        Args:
            sat_pos: satelliteECEFposition[x,y,z](m)
            rec_pos: receiverECEFposition[x,y,z](m)
            
        Returns:
            correction: 萨格纳克correction量(m)
        """
        correction = (OMEGA_EARTH / SPEED_OF_LIGHT) * (
            sat_pos[0] * rec_pos[1] - 
            sat_pos[1] * rec_pos[0]
        )
        return correction
    
    def compute_atmospheric_effects(self, time, prn, eph, recpos, height, delt):
        """
        Calculateatmospheric effecttotal影响
        """
        # Calculatetropospheredelay
        trop_delay = self.compute_tropospheric_delay_rate(time, prn, eph, recpos, height, delt)
        
        # Calculateionospheredelay
        iono_delay = self.compute_ionospheric_delay_rate(time, prn, eph, recpos, height, 0.1)
        #iono_delay = 0.0
        # Calculatetotal大气delayrate of change
        total_delay_rate = trop_delay + iono_delay
        
        # Convert为Doppler频移
        atmos_doppler = -total_delay_rate * F_L1/SPEED_OF_LIGHT
        
        return atmos_doppler, trop_delay, iono_delay
        
    def compute_tropospheric_delay_rate(self, time, prn, eph, recpos, height, delt):
        """
        Calculatetropospheredelayrate of change
        """
        
        time_former = time - delt
        time_latter = time + delt
        rec_lat, rec_lon, rec_alt = ecef_to_geodetic(*recpos)

        sat_pos_former, sat_vel_former = self.cal_sat_posvel(time_former, prn, eph, recpos)
        sat_pos_latter, sat_vel_latter = self.cal_sat_posvel(time_latter, prn, eph, recpos)

        # Calculateelevation和azimuth
        elevation_former,__ = calculate_elevation_angle(sat_pos_former, recpos)
        delays_former = self.tropmodel([rec_lat, rec_lon, rec_alt],elevation_former* (np.pi / 180.0), 0.1)

        elevation_latter,__ = calculate_elevation_angle(sat_pos_latter, recpos)
        delays_latter = self.tropmodel([rec_lat, rec_lon, rec_alt],elevation_latter* (np.pi / 180.0), 0.1)
        
        # delayrate of change
        delay_rate = (delays_latter - delays_former)/(2.0 * delt)
        
        return delay_rate

    def tropmodel(self, pos, ele, humi):
        """
        Calculatetropospheredelay(Saastamoinen模型)
        Args:
            time: 观测时间
            pos: 测站position[lat(rad), lon(rad), height(m)]
            azel: azimuth和elevation[az(rad), el(rad)]
            humi: 相对湿度(0-1)
        Returns:
            tropospheredelay(m)
        """
        TEMP0 = 15.0  # 海平面温度(℃)
        PI = np.pi
        
        # Check输入Argsvalid性
        if pos[2] < -100.0 or pos[2] > 1E4 or ele <= 0:
            return 0.0
        
        # 标准大气模型
        # 处理负altitude情况
        hgt = max(0.0, pos[2])
        
        # Calculate气压(mbar)
        pres = 1013.25 * pow(1.0 - 2.2557E-5 * hgt, 5.2568)
        
        # Calculate温度(K)
        temp = TEMP0 - 6.5E-3 * hgt + 273.16
        
        # Calculate水汽压(mbar)
        e = 6.108 * humi * np.exp((17.15 * temp - 4684.0)/(temp - 38.45))
        
        # Saastamoinen模型
        z = PI/2.0 - ele  # 天顶角(rad)
        
        # Calculate干delay
        trph = 0.0022768 * pres / (1.0 - 0.00266 * np.cos(2.0 * pos[0]) - 0.00028 * hgt/1E3) / np.cos(z)
        
        # Calculate湿delay
        trpw = 0.002277 * (1255.0/temp + 0.05) * e / np.cos(z)
        
        # Returnstotaldelay
        return trph + trpw

    def compute_ionospheric_delay_rate(self, time, prn, eph, recpos, height, delt):
        """
        Calculateionospheredelayrate of change
        """
        
        time_former = time - delt
        time_latter = time + delt
        rec_lat, rec_lon, rec_alt = ecef_to_geodetic(*recpos)

        sat_pos_former, sat_vel_former = self.cal_sat_posvel(time_former, prn, eph, recpos)
        sat_pos_latter, sat_vel_latter = self.cal_sat_posvel(time_latter, prn, eph, recpos)

        # Calculateelevation和azimuth
        azel_former = calculate_elevation_angle(sat_pos_former, recpos)
        delays_former = self.ionmodel(time_former, self.ionpara, [rec_lat, rec_lon, rec_alt], np.deg2rad(azel_former))

        azel_latter = calculate_elevation_angle(sat_pos_latter, recpos)
        delays_latter = self.ionmodel(time_latter, self.ionpara, [rec_lat, rec_lon, rec_alt], np.deg2rad(azel_latter))
        
        # delayrate of change
        delay_rate = (delays_latter - delays_former)/(2.0 * delt)
        
        return delay_rate

    def ionmodel(self, t, ion, pos, azel):
        """Calculateionospheredelay(Klobuchar模型)
        
        Args:
            t (float): GPS time
            ion (numpy.ndarray): ionosphere模型Args [a0,a1,a2,a3,b0,b1,b2,b3]
            pos (numpy.ndarray): receiverposition [lat,lon,h] (弧度,米)
            azel (numpy.ndarray): azimuth/elevation angle [az,el] (弧度)
        
        Returns:
            float: L1频段ionospheredelay(米)
        """
        
        # defaultionosphereArgs (2004/1/1)
        ion_default = np.array([
            0.1118E-07, -0.7451E-08, -0.5961E-07, 0.1192E-06,
            0.1167E+06, -0.2294E+06, -0.1311E+06, 0.1049E+07
        ])
        
        CLIGHT = 299792458.0  # 光速(m/s)
        
        # Check输入valid性
        if pos[2] < -1E3 or azel[1] <= 0:
            return 0.0
        
        if np.linalg.norm(ion) <= 0.0:
            ion = ion_default
        
        # Calculate地心角(半圆)
        psi = 0.0137 / (azel[1]/math.pi + 0.11) - 0.022
        
        # Calculateionosphere穿刺点latitude/longitude(半圆)
        phi = pos[0]/math.pi + psi * np.cos(azel[0])
        phi = np.clip(phi, -0.416, 0.416)
        lam = pos[1]/math.pi + psi * np.sin(azel[0])/np.cos(phi*math.pi)
        
        # Calculate地磁latitude(半圆)
        phi += 0.064 * np.cos((lam - 1.617) * math.pi)
        
        # Calculate当地时间(秒)
        tt = 43200.0 * lam + t
        tt -= math.floor(tt/86400.0) * 86400.0  # 确保 0 <= tt < 86400
        
        # Calculate倾斜因子
        f = 1.0 + 16.0 * pow(0.53 - azel[1]/math.pi, 3.0)
        
        # Calculateionospheredelay
        amp = ion[0] + phi * (ion[1] + phi * (ion[2] + phi * ion[3]))
        per = ion[4] + phi * (ion[5] + phi * (ion[6] + phi * ion[7]))
        
        # 限制amplitude和periodrange
        amp = max(0.0, amp)
        per = max(72000.0, per)
        
        # Calculatephase
        x = 2.0 * math.pi * (tt - 50400.0) / per
        
        # Returnsionospheredelay
        if abs(x) < 1.57:
            delay = CLIGHT * f * (5E-9 + amp * (1.0 + x*x * (-0.5 + x*x/24.0)))
        else:
            delay = CLIGHT * f * 5E-9
            
        return delay

    def compute_earth_rotation_effect(self, rec_pos, sat_pos, sat_vel):
        """
        CalculateEarth rotation对Doppler观测影响
        """
        # Earth rotation角velocityvector(ECEF)
        omega_earth = np.array([0, 0, OMEGA_EARTH])
        
        # CalculateEarth rotation引起velocity
        rec_vel_rot = np.cross(omega_earth, rec_pos)
        sat_vel_rot = np.cross(omega_earth, sat_pos)
        
        # Calculate相对运动引起Doppler效应
        los = sat_pos - rec_pos
        los_unit = los / np.linalg.norm(los)
        
        rotation_effect = np.dot(sat_vel_rot - rec_vel_rot, los_unit)
        
        return rotation_effect

    def compute_relativistic_effect(self, sat_pos, sat_vel, rec_pos):
        """
        Calculaterelativistic效应对Doppler观测影响 (unit: m/s)
        
        主要包括: 
        1. satellite钟period性relativistic改正时间derivative
        2. Shapirodelay时间derivative（通常很小, 可忽略）
        
        对于GNSS, 最重要是椭圆orbit带来period性relativistic效应。
        satellite钟relativistic改正: Δt_rel = -2(r·v)/c²
        对Doppler贡献需要对时间求导。
        """
        c = SPEED_OF_LIGHT
        gm = MU_EARTH
        
        r_sat = np.linalg.norm(sat_pos)
        v_sat = np.linalg.norm(sat_vel)
        
        # 1. period性relativistic效应对Doppler贡献
        # satellite钟period性改正: Δt_rel = -2(r·v)/c²
        # 其时间derivative: d(Δt_rel)/dt = -2(|v|² + r·a)/c²
        # 由于向心加velocity r·a ≈ -GM/r (对于近圆orbit)
        # 所以: d(Δt_rel)/dt ≈ -2(|v|² - GM/r)/c²
        # # Convert为Dopplervelocity贡献 (m/s):
        # Δv_rel = c × d(Δt_rel)/dt = -2(|v|² - GM/r)/c
        
        rel_doppler = -2.0 * (v_sat**2 - gm/r_sat) / c  # m/s
        
        # 2. Shapirodelay对Doppler贡献（通常很小, 约0.001 m/s量级）
        # Shapirodelay: Δt_sh = (2GM/c³) × ln((r_sat + r_rec + ρ)/(r_sat + r_rec - ρ))
        # 其时间derivative较小, 这里简化处理
        r_rec = np.linalg.norm(rec_pos)
        los = sat_pos - rec_pos
        range_ = np.linalg.norm(los)
        
        # line of sight方向unitvector
        e_los = los / range_
        # satellite径向velocity分量
        v_radial = np.dot(sat_vel, e_los)
        
        # Shapirodelay时间derivative近似
        # d(Δt_sh)/dt ≈ (2GM/c³) × d/dt[ln((r_sat+r_rec+ρ)/(r_sat+r_rec-ρ))]
        # 简化Calculate（量级很小）
        denom1 = r_sat + r_rec + range_
        denom2 = r_sat + r_rec - range_
        if denom2 > 1.0:  # 避免除零
            shapiro_rate = (4.0 * gm / (c**3)) * v_radial / (denom1 * denom2) * (r_sat + r_rec)
            shapiro_doppler = c * shapiro_rate  # 转换为m/s
        else:
            shapiro_doppler = 0.0
        
        total_rel_doppler = rel_doppler + shapiro_doppler
        
        return total_rel_doppler

    
    def calculate_snr(self, elevation):
            """
            according tosatelliteelevationCalculateSNR(SNR)
            
            Args:
                elevation: satelliteelevation angle(rad)
                
            Returns:
                snr: SNR(dB-Hz)
            """
            # based onelevationSNR模型
            # 典型值: 低elevation angle约35dB-Hz, 高elevation angle约50dB-Hz
            min_snr = 35.0  # 最低SNR(dB-Hz)
            max_snr = 70.0  # 最高SNR(dB-Hz)
            
            # 将elevation从弧度转为度
            elev_deg = elevation * 180 / math.pi
            
            # 线性映射: 0度->min_snr, 90度->max_snr
            snr = min_snr + (max_snr - min_snr) * elev_deg / 90.0
            
            # 加入少量随机噪声(±1dB)
            snr += np.random.uniform(-1.0, 1.0)
            
            # 限制在合理range内
            snr = max(min_snr - 5, min(max_snr + 5, snr))
            
            return snr

    def compute_pseudorange(self, times, sat_pos, sat_vel, receiver_pos,receiver_vel,snr_values):
        """
        Calculatepseudorange, consideringall误差源
        """

        los = receiver_pos - sat_pos
        range_ = np.linalg.norm(los)
        los_unit = los / range_
        
        # Calculateelevation和azimuth
        azel = calculate_elevation_angle(sat_pos, receiver_pos)
        rec_lat, rec_lon, rec_alt = ecef_to_geodetic(*receiver_pos)     

        # 1. Calculateatmospheric effect
        if self.enable_tropo:              
            # tropospheredelay
            trop_delay = self.tropmodel([rec_lat, rec_lon, rec_alt], azel[0]* (np.pi / 180.0), 0.1)     
        else:
            trop_delay = 0.0    
        # ionospheredelay
        if self.enable_iono: 
            iono_delay = self.ionmodel(times, self.ionpara, [rec_lat, rec_lon, rec_alt], np.deg2rad(azel))       
        else:
            iono_delay = 0.0
        
        # 3. Sagnac effect(Earth rotation引起信号传播delay)
        omega_e = OMEGA_EARTH  # 地球自转角速度
        rotation_effect = (omega_e * (sat_pos[0] * receiver_pos[1] - sat_pos[1] * receiver_pos[0])) / SPEED_OF_LIGHT

        # 4. multipath效应
        if self.enable_multipath:           
            if snr_values is not None:
                # based onelevation和signal-to-noise ratio简单multipath模型
                elevation_rad = np.deg2rad(azel[0])
                snr = snr_values
                # multipath影响随elevation增加而减小, 随signal-to-noise ratio增加而减小
                multipath_effect = 2.0 * (1.0 - np.sin(elevation_rad)) * np.exp(-0.05 * (snr - 30))
                multipath_effect = np.clip(multipath_effect, 0, 5.0)  # 限制在5米以内
            else:
                multipath_effect = 0.0
        else:
            multipath_effect = 0.0
        
        # 5. relativistic效应(主要部分已在satelliteclock bias中considering)
        if self.enable_relativistic:
            # 附加relativistic效应: 引力势差异导致时间delay
            gm = MU_EARTH
            r_sat = np.linalg.norm(sat_pos)
            r_rec = np.linalg.norm(receiver_pos)
            relativity_extra = 2 * gm * np.log((r_sat + r_rec + range_)/(r_sat + r_rec - range_)) / (SPEED_OF_LIGHT**2)
        else:
            relativity_extra = 0.0
        
        # 6. hardware delay
        if self.enable_hardware:
            hardware_delay = 0.0  # 默认值, 实际应从外部获取或建模
        else:
            hardware_delay = 0.0
        
        # 7. phase center偏移
        if self.enable_pco:
            # 简单模型: elevation相关天线phase center偏移
            if azel[0] < 15:
                pco_effect = 0.05 * (15 - azel[0])  # 低高度角下额外偏移
        else:
            pco_effect = 0.0
        
        # Calculatereceiverclock bias在current时刻值(consideringclock drift)
        clock_bias_at_t = self.clock
        # clock_bias_at_t += self.clockdrift * (times - self.start_time)
        sat_clock_bias = 0
        # Calculate预测pseudorange(includingall误差源)
        pred_pseudorange = (
            range_ +                    # 几何距离
            (clock_bias_at_t - sat_clock_bias) +  # 钟差(接收机-卫星)
            trop_delay +                # 对流层延迟
            iono_delay +                # 电离层延迟
            rotation_effect +           # 地球自转效应
            multipath_effect +          # 多路径效应
            relativity_extra +  # 附加相对论效应
            hardware_delay +            # 硬件延迟
            pco_effect                  # 相位中心偏移
        )
               
        return pred_pseudorange
    
    def compute_multipath_effect(self, elevation, snr=None):
        """
        估计multipath效应影响
        usingelevation和signal-to-noise ratio(if有)进line简单建模
        """
        # based onelevation简单余弦模型
        elevation_rad = np.deg2rad(elevation)
        base_mp = 0.1 * (1 - np.cos(elevation_rad))  # 基础多路径影响
        
        # if有signal-to-noise ratio信息, 进一步调整
        if snr is not None:
            # SNR越高, multipath影响越小
            snr_factor = np.exp(-0.1 * (snr - 50)) if snr > 50 else 1.0
            base_mp *= snr_factor
        
        return -base_mp * F_L1/SPEED_OF_LIGHT
    
    def compute_doppler(self, times, prnlist,eph, sat_pos, sat_vel, receiver_pos,receiver_vel,snr_values):
        """
        CalculateDoppler, consideringall误差源
        """

        los = receiver_pos - sat_pos
        range_ = np.linalg.norm(los)
        los_unit = los / range_
        
        # Calculateelevation和azimuth
        azel = calculate_elevation_angle(sat_pos, receiver_pos)
        rec_lat, rec_lon, rec_alt = ecef_to_geodetic(*receiver_pos)     

        # 1. Calculateatmospheric effect
        # trop
        if self.enable_tropo:
            trop_delay = self.compute_tropospheric_delay_rate(times, prnlist, eph, receiver_pos, rec_alt, 0.005)
        else:
            trop_delay = 0
        # iono
        if self.enable_iono:
            iono_delay = self.compute_ionospheric_delay_rate(times, prnlist, eph, receiver_pos, rec_alt, 0.1)
        else:
            iono_delay = 0.0
        atmos_doppler = - (trop_delay + iono_delay) * F_L1/SPEED_OF_LIGHT

        if self.enable_relativistic:
            relativistic_effect = self.compute_relativistic_effect(sat_pos, sat_vel, receiver_pos)
            relativistic_effect = -relativistic_effect * F_L1/SPEED_OF_LIGHT
        else:
            relativistic_effect = 0.0

        if self.enable_hardware:
            hardware_delay = 0.0  # 默认值, 实际应从外部获取或建模
        else:
            hardware_delay = 0.0

        if self.enable_pco:
            pco_effect = 0.0
        else:
            pco_effect = 0.0

        # CalculateEarth rotation效应
        rotation_effect = self.compute_earth_rotation_effect(receiver_pos, sat_pos, sat_vel)
        rotation_doppler = -rotation_effect * F_L1/SPEED_OF_LIGHT
        
        # Calculatemultipath效应
        if self.enable_multipath:
            multipath_effect = self.compute_multipath_effect(
                azel[0], 
                snr_values if snr_values is not None else None
            )
        else:
            multipath_effect = 0.0
        
        
        # Calculate几何Doppler
        geometric_doppler = F_L1/SPEED_OF_LIGHT * np.dot(sat_vel, los_unit)

        clock_drift_Hz = self.clockdrift/SPEED_OF_LIGHT*F_L1
        # 合成total预测Doppler
        pred_doppler = (geometric_doppler + 
                        atmos_doppler +      # 大气效应
                        rotation_doppler +    # 地球自转
                        multipath_effect +    # 多路径
                        pco_effect + # 相位中心
                        relativistic_effect + # 相对论
                        hardware_delay +     # 硬件延迟
                        clock_drift_Hz)         # 钟漂
   
        return pred_doppler

    def compute_observations(self, prn, ephemeris, receiver_pos, receiver_vel, gps_time, prnoise=None,dopnoise=None):
        """
        Calculatepseudorange、Doppler和carrier phaseobservation
        
        Args:
            prn: satellitePRN号
            ephemeris: satelliteephemeris字典或TLEsatellite对象
            receiver_pos: receiverECEFposition[x,y,z](m)
            receiver_vel: receiverECEFvelocity[vx,vy,vz](m/s)
            gps_time: GPS time(s), 接收时间
            clock_bias: receiverclock bias(s)
            clock_drift: receiverclock drift(s/s)
            
        Returns:
            observation字典, includingpseudorange、Doppler和carrier phase
        """
        
        # usingiteration法Calculate信号发射时间
        emission_time = gps_time
        max_iter = 10
        
        # initial时延估计
        time_delay = 0.01  # 初始猜测为10ms
        
        for _ in range(max_iter):
            # 在current估计发射时间Calculatesatelliteposition
            # CheckwhetherusingTLE模式
            if prn in self.tle_satellites:
                tle_eph = {'use_tle': True, 'prn': prn}
                sat_pos, sat_vel = self.calculate_satellite_position_tle(prn, emission_time)
            else:
                sat_pos, sat_vel = self.calculate_satellite_position(ephemeris, emission_time)
            
            # Calculatesatellite与receivergeometric distance
            los = receiver_pos - sat_pos
            geometric_range = np.linalg.norm(los)
            
            # Update信号传播时延
            new_time_delay = geometric_range / SPEED_OF_LIGHT
            
            # Update发射时间
            new_emission_time = gps_time - new_time_delay
            
            # Check收敛条件
            if abs(new_emission_time - emission_time) < 1e-6:
                break
            
            emission_time = new_emission_time
            time_delay = new_time_delay
        
        
        elevation, azimuth = calculate_elevation_angle(sat_pos, receiver_pos)

        # CalculateSNR(SNR)
        snr = self.calculate_snr(np.radians(elevation))
        
        # Calculatepseudorange
        pseudorange=self.compute_pseudorange(emission_time, sat_pos, sat_vel, receiver_pos,receiver_vel,snr)
        
        pseudorange+=np.random.normal(0, prnoise if prnoise is not None else 10)

        # CalculateDoppler
        doppler=self.compute_doppler(emission_time, prn,ephemeris, sat_pos, sat_vel, receiver_pos,receiver_vel,snr)
        
        doppler += np.random.normal(0, dopnoise if dopnoise is not None else 0.1)
        
        snr_noise_factor = 1.0 / (10**(snr / 20))
        carrier_phase = pseudorange + np.random.normal(0, 0.01 * snr_noise_factor)
    
        return {
            'pseudorange': pseudorange,
            'doppler': doppler,
            'carrier_phase': carrier_phase,
            'satellite_position': sat_pos,
            'satellite_velocity': sat_vel,
            'emission_time': emission_time,
            'time_delay': time_delay,
            'snr': snr,           
            'elevation': elevation,
            'azimuth': azimuth
        }

def main(week,start_time, end_time, output_file, tle_file=None, clk_data=None, clk_drift_data=None,receiver_pos=None):
    # Create仿真系统
    sim = SatelliteSimulation(enable_iono=False, enable_tropo=False, enable_multipath=False)
    
    # SetGPS week（可according to实际情况修改）
    sim.gps_week = week  # 例如2022年某周
    
    # 也可以直接SetECEF坐标
    if receiver_pos is None:
        receiver_pos = np.array([-2173447.0803, 4397916.0083, 4062634.4708])
        receiver_vel = np.array([0.0, 0.0, 0.0])  # 静态接收机
    receiver_vel = np.array([0.0, 0.0, 0.0])  # 静态接收机
    # receiverclock bias和clock drift
    
    clock_bias = 1000 #0  # 钟差(秒)
    clock_drift = 0.05  # 钟漂(秒/秒)

    # sim.clock0 = clock_bias
    # sim.clockdrift = clock_drift
    sim.start_time = start_time
    # satellitedata来源
    satellites = []
    
    # usingTLEfile
    if tle_file and os.path.exists(tle_file):
        # LoadTLEfile
        print(f"加载TLE文件: {tle_file}")
        count = sim.load_tle_file(tle_file)
        if count > 0:
            # usingTLEdata
            for prn, sat_data in sim.tle_satellites.items():
                satellites.append({
                    'prn': prn,
                    'name': sat_data['name'],
                    'use_tle': True
                })
            print(f"使用TLE模式, 已加载{len(satellites)}卫星")
    else:
        # usingdefaultephemeris
        satellites = [
            {
                'prn': 24,  # 卫星PRN号
                'ephemeris': {
                    'toc': 437100,  # 参考时间(GPS周内秒)
                    'toe': 437100,  # 参考时间(GPS周内秒)
                    'sqrta': 2738.2328242865,  # 轨道半长轴平方根(m^0.5)
                    'es': 0.0021833032,  # 偏心率
                    'inc0': 1.5097466027,  # 轨道倾角(rad)
                    'omega0': -0.4934195498,  # 升交点赤经(rad)
                    'omega': 1.3370709787,  # 近地点角距(rad)
                    'm0': 0.9666686747,  # 平近点角(rad)
                    'dotn': -0.0000002690,  # 平均运动差值(rad/s)
                    'idot': -0.0000000076,  # 轨道倾角变化率(rad/s)
                    'omegadot': -0.0000000705,  # 升交点赤经变化率(rad/s)
                    'cuc': 0.0000007841,  # 纬度幅余弦调和修正项(rad)
                    'cus': 0.0000984017,  # 纬度幅正弦调和修正项(rad)
                    'crc': 1466.5468750000,  # 轨道半径余弦调和修正项(m)
                    'crs': 28.1484375000,  # 轨道半径正弦调和修正项(m)
                    'cic': -0.0000007081,  # 轨道倾角余弦调和修正项(rad)
                    'cis': 0.0000011557,  # 轨道倾角正弦调和修正项(rad)
                    'af0': 0.0,  # 卫星钟偏差(秒)
                    'af1': 0.0,  # 卫星钟漂(秒/秒)
                    'af2': 0.0,  # 卫星钟漂变化(秒/秒^2)
                    'tgd1': 0.0,  # 卫星钟偏差(秒)
                }
            }
        ]
        print("使用默认星历模式")
    i=0
    # 打开输出file
    with open(output_file, 'w', encoding='utf-8') as f:
        # 在时间range内进line仿真
        gps_time = start_time
        while gps_time <= end_time:
            # 对每satellite进line仿真
            for sat in satellites:
                try:
                    # Calculateclock biasclock drift
                    if clk_data is not None:
                        sim.clock = clk_data[i]
                        sim.clockdrift = clk_drift_data[i]                       
                    else:
                        sim.clockdrift = clock_drift
                        sim.clock = clock_bias + sim.clockdrift * (gps_time - sim.start_time)
                    # Calculateobservation
                    if 'use_tle' in sat and sat['use_tle']:
                        # usingTLECalculate
                        observations = sim.compute_observations(
                            sat['prn'],
                            None,  # 使用TLE模式不需要传递ephemeris 
                            receiver_pos,
                            receiver_vel,
                            gps_time
                        )
                    else:
                        # usingephemerisCalculate
                        observations = sim.compute_observations(
                            sat['prn'],
                            sat['ephemeris'], 
                            receiver_pos,
                            receiver_vel,
                            gps_time
                        )
                    
                    # 判断satellitewhether可见（elevation angle大于0度）
                    if observations['elevation'] > 0:
                        # 将carrier phase从米Convert为周
                        carrier_phase_cycles = observations['carrier_phase'] / LAMBDA_L1
                        
                        # Writedataline
                        f.write(f"OBS: 1  {observations['emission_time']:.9f}    "
                                f"{gps_time:.9f}      "
                                f"{observations['pseudorange']:.2f}    "
                                f"{observations['doppler']:.5f}    "
                                f"{observations['snr']:.1f}   "
                                f"6   1    1   "
                                f"{observations['elevation']:.2f}  "
                                f"{sat['prn']}\n")
                except Exception as e:
                    print(f"计算卫星{sat['prn']}观测值时出错: {e}")
                    continue
            
            # 时间前进
            gps_time += 1  # 每1秒采样一次, 可根据需要调整
            i = i + 1

    print(f"仿真完成!结果已保存到 {output_file}")

def create_example_tle_file(filename, num_satellites=10):
    """
    CreateExampleTLEfile, for测试
    
    Args:
        filename: 输出file名
        num_satellites: Number of satellites
    """
    # 一些ExampleTLEdata（实际using时应替换为真实TLEdata）
    example_tles = [
        ("GSAT0101 (GALILEO-PFM)", "1 37846U 11060A   22083.18124909 -.00000044  00000+0  00000+0 0  9997", "2 37846  57.5580 216.3818 0002314 230.1501 129.8866  1.70475366 63964"),
        ("GSAT0102 (GALILEO-FM2)", "1 37847U 11060B   22083.38254294 -.00000044  00000+0  00000+0 0  9999", "2 37847  57.5583 216.3798 0002196 254.1642 105.8666  1.70475486 63943"),
        ("COSMOS 2544", "1 44850U 19096A   22082.84188225  .00000080  00000+0  12196-3 0  9997", "2 44850  67.1490 213.7361 0007131 355.7135 101.9622 14.13412967119638"),
        ("COSMOS 2545", "1 44851U 19096B   22082.85561090  .00000071  00000+0  11115-3 0  9997", "2 44851  67.1488 213.6880 0003874 359.8654 158.2839 14.13412972119593"),
        ("STARLINK-3042", "1 48276U 21044BT  22082.60702309  .00011137  00000+0  73906-3 0  9994", "2 48276  53.0507 185.5547 0001532  45.1644 314.9446 15.06421459 37698"),
        ("STARLINK-3043", "1 48277U 21044BU  22082.60702493  .00009633  00000+0  65245-3 0  9992", "2 48277  53.0493 185.5499 0001426  95.2064 264.9018 15.06382234 37692"),
        ("ISS (ZARYA)", "1 25544U 98067A   22083.78183125  .00006396  00000+0  12403-3 0  9993", "2 25544  51.6423 170.7326 0003825 100.1724 336.5860 15.49558623329741"),
        ("NOAA 19", "1 33591U 09005A   22083.51271506  .00000075  00000+0  65529-4 0  9996", "2 33591  99.1839 142.1071 0013621 185.9449 174.1558 14.12501737675838"),
        ("GOES 16", "1 41866U 16071A   22082.96132094 -.00000225  00000+0  00000+0 0  9996", "2 41866   0.0439 264.4393 0000331 177.0131 301.2163  1.00271433 19695"),
        ("ONEWEB-0051", "1 47438U 21006C   22083.63483587  .00000030  00000+0  34953-4 0  9995", "2 47438  87.8995 128.8143 0001257  90.6790  76.3975 13.02336353 59134"),
        ("GPS SVN78", "1 44506U 19056A   22083.46116370 -.00000050  00000+0  00000+0 0  9990", "2 44506  55.1242  16.5363 0010989 185.5485 174.5115  2.00563122 19013"),
        ("GPS SVN74", "1 41328U 16007A   22083.71851852 -.00000046  00000+0  00000+0 0  9993", "2 41328  54.9798 136.1888 0025699 203.2781 156.6452  2.00566189 45179"),
        ("BEIDOU 3 GEO-1", "1 43001U 17069A   22083.20496540 -.00000323  00000+0  00000+0 0  9992", "2 43001   1.6022 238.9626 0003638 240.0896 280.0574  1.00276562 15878"),
        ("GLONASS-M", "1 46837U 20074A   22083.61553359 -.00000062  00000+0  00000+0 0  9993", "2 46837  64.8250  14.8414 0001941 190.0250 169.9964  2.13101922 11204"),
        ("GSAT0219 (GALILEO 22)", "1 43057U 17079B   22082.92508367 -.00000045  00000+0  00000+0 0  9992", "2 43057  57.1139  96.1701 0002643 344.9193  15.1587  1.70474953 29155")
    ]
    
    # 确保Number of satellites不超过Exampledata
    num_satellites = min(num_satellites, len(example_tles))
    
    # Create输出目录（if不存在）
    os.makedirs(os.path.dirname(filename) if os.path.dirname(filename) else '.', exist_ok=True)
    
    # WriteExampleTLEfile
    with open(filename, 'w') as f:
        for i in range(num_satellites):
            name, line1, line2 = example_tles[i]
            f.write(f"{name}\n{line1}\n{line2}\n")
    
    print(f"创建了示例TLE文件: {filename}, 包含{num_satellites}卫星")

def plot_doppler_observations(obs_file, output_file=None):
    """
    Read观测file并绘制Dopplerobservation随时间变化图, 颜色表示elevation或SNR
    
    Args:
        obs_file: 观测datafilepath
        output_file: 输出图像filepath(optional)
    """
    # Readdata, 处理特定format
    times = []
    pseudoranges = []
    dopplers = []
    snrs = []
    elevations = []
    prns = []
    
    with open(obs_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('OBS:'):
                # 按空格分割, 处理field之间多空格情况
                fields = [field for field in line.strip().split() if field]
                
                if len(fields) >= 12:  # 确保有足够字段
                    # 提取各columndata
                    try:
                        # 第4column是时间, 第5column是pseudorange, 第6column是Doppler, 第7column是SNR, 第11column是elevation, 第12column是satellite号
                        time = float(fields[3])  # 第4列: 时间
                        pseudorange = float(fields[4])  # 第5列: 伪距
                        doppler = float(fields[5])  # 第6列: 多普勒
                        snr = float(fields[6])  # 第7列: SNR
                        elevation = float(fields[10])  # 第11列: 高度角
                        prn = int(fields[11])  # 第12列: 卫星号
                        
                        times.append(time)
                        pseudoranges.append(pseudorange)
                        dopplers.append(doppler)
                        snrs.append(snr)
                        elevations.append(elevation)
                        prns.append(prn)
                    except (ValueError, IndexError) as e:
                        print(f"解析行时出错: {line.strip()}")
                        print(f"错误: {e}")
                        continue
    
    # Convert为NumPy数组以便于后续处理
    times = np.array(times)
    pseudoranges = np.array(pseudoranges)
    dopplers = np.array(dopplers)
    snrs = np.array(snrs)
    elevations = np.array(elevations)
    prns = np.array(prns)
    
    print(f"成功读取数据, 共{len(times)}观测")
    if len(times) == 0:
        print("没有读取到数据, 请检查文件格式")
        return
    
    # Set中文字体
    plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'STSong', 'Arial Unicode MS']
    plt.rcParams['axes.unicode_minus'] = False
    
    # Create图形
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
    
    # Get唯一PRN值
    unique_prns = np.unique(prns)
    
    # 为不同PRN分配不同颜色
    cmap = cm.get_cmap('tab10')
    
    # 绘制Doppler随时间变化, 颜色表示elevation
    sc1 = ax1.scatter(times, dopplers, c=elevations, cmap='viridis', 
                     s=20, alpha=0.8, edgecolor='none')
    
    ax1.set_title('Doppler Shift vs Time (Color indicates Elevation)')
    ax1.set_ylabel('Doppler Shift (Hz)')
    ax1.grid(True, linestyle='--', alpha=0.7)
    
    # 添加颜色条
    cbar1 = plt.colorbar(sc1, ax=ax1)
    cbar1.set_label('高度角 (度)')

    # 第二张图   
    # 绘制elevation随时间变化, 颜色表示SNR
    sc2 = ax2.scatter(times, elevations, c=snrs, cmap='viridis', 
                     s=20, alpha=0.8, edgecolor='none')
    
    ax2.set_title('Elevation vs Time (Color indicates SNR)')
    ax2.set_ylabel('Elevation (deg)')
    ax2.grid(True, linestyle='--', alpha=0.7)
    
    # 找出maximumelevation及其对应时间点
    if len(elevations) > 0:
        max_elevation_idx = np.argmax(elevations)
        max_elevation = elevations[max_elevation_idx]
        max_elevation_time = times[max_elevation_idx]
        max_elevation_prn = prns[max_elevation_idx]
        
        # 在图上标出maximumelevation点
        ax2.plot(max_elevation_time, max_elevation, 'ro', markersize=10)
        ax2.annotate(f'最大高度角: {max_elevation:.1f}°\nPRN: {max_elevation_prn}',
                    xy=(max_elevation_time, max_elevation),
                    xytext=(max_elevation_time, max_elevation - 10),  # 文本位置在点下方
                    arrowprops=dict(facecolor='red', shrink=0.05, width=2.0, headwidth=10, headlength=10),
                    fontsize=12, ha='center')
    
    # 添加颜色条
    cbar2 = plt.colorbar(sc2, ax=ax2)
    cbar2.set_label('SNR (dB-Hz)')
    
    # # 给不同PRN分组并绘制折线
    # for i, sat_prn in enumerate(unique_prns):
    #     mask = prns == sat_prn
    # if np.sum(mask) > 0:  # 确保有data
    #         color = cmap(i % 10)
    #         ax2.plot(times[mask], dopplers[mask], 'o-', 
    #                 color=color, label=f'PRN {int(sat_prn)}',
    #                 markersize=4, linewidth=1.5, alpha=0.8)
    
    # ax2.set_title('Doppler Shift vs Time for Each Satellite')
    # ax2.set_xlabel('Time (seconds)')
    # ax2.set_ylabel('Doppler Shift (Hz)')
    # ax2.grid(True, linestyle='--', alpha=0.7)
    # ax2.legend(loc='best', framealpha=0.7)
    
    # 调整布局
    plt.tight_layout()
    
    # if提供了输出file名, Save图像
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"图像已保存至: {output_file}")
    
    # 显示图像
    plt.show()
    
    # # 绘制SNR与elevation关系图
    # plt.figure(figsize=(10, 6))
    # plt.scatter(elevations, snrs, c=dopplers, cmap='coolwarm', alpha=0.7)
    # plt.colorbar(label='Doppler频移 (Hz)')
    # plt.title('SNR(SNR)与elevation关系')
    # plt.xlabel('Elevation (deg)')
    # plt.ylabel('SNR (dB-Hz)')
    # plt.grid(True, linestyle='--', alpha=0.7)
    
    # # 添加趋势线
    # if len(elevations) > 2:
    #     z = np.polyfit(elevations, snrs, 1)
    #     p = np.poly1d(z)
    #     plt.plot(sorted(elevations), p(sorted(elevations)), 'r--', 
    # label=f'拟合线: SNR = {z[0]:.2f}*elevation + {z[1]:.2f}')
    #     plt.legend()
    
    # plt.tight_layout()
    
    # # Save第二张图
    # if output_file:
    #     snr_file = output_file.replace('.', '_snr.')
    #     plt.savefig(snr_file, dpi=300, bbox_inches='tight')
    # print(f"SNRImage saved to: {snr_file}")
    
    # plt.show()
    
    # 统计信息
    print(f"\n观测数据统计:")
    print(f"总观测点数: {len(times)}")
    print(f"卫星数量: {len(unique_prns)}")
    print(f"时间range: {min(times):.1f}s - {max(times):.1f}s")
    print(f"多普勒range: {min(dopplers):.1f}Hz - {max(dopplers):.1f}Hz")
    print(f"高度角range: {min(elevations):.1f}° - {max(elevations):.1f}°")
    if len(elevations) > 0:
        print(f"最大高度角: {max_elevation:.1f}°, 对应卫星PRN: {max_elevation_prn}, 时间: {max_elevation_time:.1f}s")
    print(f"SNRrange: {min(snrs):.1f}dB-Hz - {max(snrs):.1f}dB-Hz")

def read_rck_data(filename,clk_driftfilename):
    """ReadRCKdatafile"""
    times = []
    rck_gps = []
    
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('*'):  # 跳过注释行
                continue
            try:
                # 解析dataline
                fields = line.strip().split()
                if len(fields) >= 7:  # 确保有足够字段
                    # 构建时间
                    year = int(fields[0])
                    month = int(fields[1])
                    day = int(fields[2])
                    hour = int(fields[3])
                    minute = int(fields[4])
                    second = float(fields[5])
                    
                    # Createdatetime对象
                    dt = datetime(year, month, day, hour, minute, int(second))
                    
                    # Convert为相对时间（分钟）
                    start_time = datetime(year, month, day, 0, 0, 0)
                    time_minutes = (dt - start_time).total_seconds() / 60
                    
                    # ReadRCK(GPS)值
                    rck_value = float(fields[6])
                    
                    times.append(time_minutes)
                    rck_gps.append(rck_value)
            except (ValueError, IndexError) as e:
                print(f"解析行时出错: {line.strip()}")
                continue

        # """Readderivative结果file"""
    times = []
    clk_drift = []
    with open(clk_driftfilename, 'r') as f:
        for i, line in enumerate(f):
            if i == 0 or i==1:
                continue  # 跳过表头
            fields = line.strip().split()
            if len(fields) >= 2:
                times.append(float(fields[0]))
                clk_drift.append(float(fields[1]))


    return np.array(rck_gps), np.array(clk_drift)

def calculate_ddop(sat_positions, receiver_pos):
    """
    CalculateDopplerpositioningDDOP值
    
    Args:
        sat_positions: satellitepositioncolumn表, each元素为[x,y,z]数组
        receiver_pos: receiverposition[x,y,z]
        
    Returns:
        ddop: DDOP值
    """
    # 构建观测matrixH
    H = []
    for sat_pos in sat_positions:
        jac = np.zeros((4))
        # Calculateline of sightvector
        los = receiver_pos - sat_pos[0]
        los_norm = np.linalg.norm(los)
        los_unit = los / los_norm

        # Calculatepositionderivative (Doppler对position偏derivative)
        sat_vel = sat_pos[1]

        dlos_dpos = (np.eye(3) - np.outer(los_unit, los_unit)) / los_norm
        dDop_dpos = F_L1/SPEED_OF_LIGHT * np.dot(dlos_dpos, sat_vel)

        # 填充雅可比matrix
        jac[:3] = -dDop_dpos
        jac[3] = -SPEED_OF_LIGHT/ F_L1

        # 添加观测方程系数
        H.append(jac)
    
    H = np.array(H)
    
    # Calculate法matrix
    N = H.T @ H
    
    # Calculate协方差matrix
    try:
        P = np.linalg.inv(N)
        # CalculateDDOP
        ddop = np.sqrt(np.trace(P))
    except np.linalg.LinAlgError:
        # ifmatrix不可逆, Returns一较大值
        ddop = 999.0
        
    return ddop

def plot_elevation_observations(obs_file,week, receiver_pos,output_file=None):
    """
    Read观测file并绘制elevation随时间变化图, 同时CalculateDDOP
    
    Args:
        obs_file: 观测datafilepath
        output_file: 输出图像filepath(optional)
    """
    # Readdata
    times = []
    elevations = []
    prns = []
    snrs = []
    sat_positions = []  

    # Create仿真系统
    sim = SatelliteSimulation(enable_iono=False, enable_tropo=False, enable_multipath=False)
    
    # SetGPS week（可according to实际情况修改）
    sim.gps_week = week  # 例如2022年某周
    
    
    receiver_vel = np.array([0.0, 0.0, 0.0])  # 静态接收机
   
    satellites = []
 
    # LoadTLEfile
    print(f"加载TLE文件: {tle_file}")
    count = sim.load_tle_file(tle_file)
    if count > 0:
        # usingTLEdata
        for prn, sat_data in sim.tle_satellites.items():
            satellites.append({
                'prn': prn,
                'name': sat_data['name'],
                'use_tle': True
            })
        print(f"使用TLE模式, 已加载{len(satellites)}卫星")

    with open(obs_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('OBS:'):
                fields = [field for field in line.strip().split() if field]
                
                if len(fields) >= 12:
                    try:
                        trtime = float(fields[2])
                        time = float(fields[3])
                        elevation = float(fields[10])
                        prn = int(fields[11])
                        snr = float(fields[6])
                        if elevation < 5:
                            continue
                        # Calculatesatelliteposition
                        if prn in sim.tle_satellites:
                            sat_pos = sim.calculate_satellite_position_tle(prn, trtime)                        
                        
                        times.append(time)
                        elevations.append(elevation)
                        prns.append(prn)
                        snrs.append(snr)
                        sat_positions.append(sat_pos)
                    except (ValueError, IndexError) as e:
                        print(f"解析行时出错: {line.strip()}")
                        continue
    
    # Convert为NumPy数组
    times = np.array(times)
    elevations = np.array(elevations)
    prns = np.array(prns)
    snrs = np.array(snrs)
    sat_positions = np.array(sat_positions)
    
    print(f"成功读取数据, 共{len(times)}观测")
    if len(times) == 0:
        print("没有读取到数据, 请检查文件格式")
        return
    
    # CalculateDDOP
    ddop = calculate_ddop(sat_positions, receiver_pos)
    print(f"DDOP值: {ddop:.2f}")
    
    # Set中文字体
    plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'STSong', 'Arial Unicode MS']
    plt.rcParams['axes.unicode_minus'] = False
    
    # Create图形
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
    
    # Get唯一PRN值
    unique_prns = np.unique(prns)
    
    # 为不同PRN分配不同颜色
    cmap = cm.get_cmap('tab10')
    
    # 绘制elevation随时间变化, 颜色表示SNR
    sc1 = ax1.scatter(times, elevations, c=snrs, cmap='viridis', 
                     s=20, alpha=0.8, edgecolor='none')
    
    ax1.set_title('Elevation vs Time (Color indicates SNR)')
    ax1.set_ylabel('Elevation (deg)')
    ax1.grid(True, linestyle='--', alpha=0.7)
    
    # 添加颜色条
    cbar1 = plt.colorbar(sc1, ax=ax1)
    cbar1.set_label('SNR (dB-Hz)')
    
    # 找出maximumelevation及其对应时间点
    max_elevation_idx = np.argmax(elevations)
    max_elevation = elevations[max_elevation_idx]
    max_elevation_time = times[max_elevation_idx]
    max_elevation_prn = prns[max_elevation_idx]
    
    # 在图上标出maximumelevation点
    ax1.plot(max_elevation_time, max_elevation, 'ro', markersize=10)
    ax1.annotate(f'最大高度角: {max_elevation:.1f}°\nPRN: {max_elevation_prn}',
                xy=(max_elevation_time, max_elevation),
                xytext=(max_elevation_time, max_elevation - 10),  # 文本位置在点下方
                arrowprops=dict(facecolor='red', shrink=0.05, width=2.0, headwidth=10, headlength=10),
                fontsize=12, ha='center')
    
    # 第二张图   
    # 绘制elevation随时间变化, 颜色表示SNR
    sc2 = ax2.scatter(times, elevations, c=snrs, cmap='viridis', 
                     s=20, alpha=0.8, edgecolor='none')
    
    ax2.set_title('Elevation vs Time (Color indicates SNR)')
    ax2.set_ylabel('Elevation (deg)')
    ax2.grid(True, linestyle='--', alpha=0.7)
    
    # 添加颜色条
    cbar2 = plt.colorbar(sc2, ax=ax2)
    cbar2.set_label('SNR (dB-Hz)')
    
    # 添加DDOP信息
    ax2.text(0.02, 0.98, f'DDOP: {ddop:.2f}', transform=ax2.transAxes,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # 调整布局
    plt.tight_layout()
    
    # if提供了输出file名, Save图像
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"图像已保存至: {output_file}")
    
    # 显示图像
    plt.show()
    
    # 统计信息
    print(f"\n观测数据统计:")
    print(f"总观测点数: {len(times)}")
    print(f"卫星数量: {len(unique_prns)}")
    print(f"时间range: {min(times):.1f}s - {max(times):.1f}s")
    print(f"高度角range: {min(elevations):.1f}° - {max(elevations):.1f}°")
    print(f"最大高度角: {max_elevation:.1f}°, 对应卫星PRN: {max_elevation_prn}, 时间: {max_elevation_time:.1f}s")
    print(f"SNRrange: {min(snrs):.1f}dB-Hz - {max(snrs):.1f}dB-Hz")
    print(f"DDOP值: {ddop:.2f}")

def plot_global_ddop(week, tle_file, output_file=None):
    """
    绘制全球Doppler单星positioningDDOP图
    
    Args:
        week: GPS week
        tle_file: TLEfilepath
        output_file: 输出图像filepath(optional)
    """
    # Create仿真系统
    sim = SatelliteSimulation(enable_iono=False, enable_tropo=False, enable_multipath=False)
    sim.gps_week = week
    
    # LoadTLEfile
    print(f"加载TLE文件: {tle_file}")
    count = sim.load_tle_file(tle_file)
    if count == 0:
        print("加载TLE文件失败")
        return
    
    # Set采样网格
    lat_grid = np.linspace(-90, 90, 181)  # 纬度网格, 1度间隔
    lon_grid = np.linspace(-180, 180, 361)  # 经度网格, 1度间隔
    alt = 0.0  # 海平面高度
    
    # InitializeDDOPmatrix
    ddop_matrix = np.zeros((len(lat_grid), len(lon_grid)))
    
    # Set观测时间
    obs_time = 436046  # 示例时间, 可以根据需要调整
    
    # 对each网格点CalculateDDOP
    for i, lat in enumerate(lat_grid):
        for j, lon in enumerate(lon_grid):
            # Convert地理坐标到ECEF
            receiver_pos = geodetic_to_ecef(np.radians(lat), np.radians(lon), alt)
            # for obs_time in range(436046, 436046 + 1000, 100):
                # Getall可见satelliteposition和velocity
            sat_positions = []
            for prn in sim.tle_satellites:
                try:
                    sat_pos, sat_vel = sim.calculate_satellite_position_tle(prn, obs_time)
                    # Calculateelevation
                    elevation, _ = calculate_elevation_angle(sat_pos, receiver_pos)
                    if elevation > 5:  # 只考虑高度角大于5度卫星
                        sat_positions.append((sat_pos, sat_vel))
                except Exception as e:
                    continue
            
            # CalculateDDOP
            if len(sat_positions) >= 4:  # 至少需要4卫星
                ddop = calculate_ddop(sat_positions, receiver_pos)
                ddop_matrix[i, j] = ddop
            else:
                ddop_matrix[i, j] = np.nan
    
    # Set中文字体
    plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'STSong', 'Arial Unicode MS']
    plt.rcParams['axes.unicode_minus'] = False
    
    # Create图形
    plt.figure(figsize=(15, 10))
    
    # 绘制DDOP等值线图
    lon_mesh, lat_mesh = np.meshgrid(lon_grid, lat_grid)
    levels = np.linspace(0, 10, 21)  # 设置等值线级别
    
    # usingpcolormesh绘制基础图
    plt.pcolormesh(lon_mesh, lat_mesh, ddop_matrix, 
                   cmap='jet', shading='auto', 
                   vmin=0, vmax=10)
    
    # 添加等值线
    contour = plt.contour(lon_mesh, lat_mesh, ddop_matrix, 
                         levels=levels, colors='k', linewidths=0.5)
    plt.clabel(contour, inline=True, fontsize=8, fmt='%.1f')
    
    # 添加颜色条
    plt.colorbar(label='DDOP值')
    
    # Set标题和轴标签
    plt.title('全球多普勒单星定位DDOP分布图')
    plt.xlabel('Longitude (deg)')
    plt.ylabel('Latitude (deg)')
    
    # Set网格
    plt.grid(True, linestyle='--', alpha=0.3)
    
    # 调整布局
    plt.tight_layout()
    
    # if提供了输出file名, Save图像
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"图像已保存至: {output_file}")
    
    # 显示图像
    plt.show()
    
    # 打印统计信息
    valid_ddop = ddop_matrix[~np.isnan(ddop_matrix)]
    if len(valid_ddop) > 0:
        print(f"\nDDOP统计信息:")
        print(f"最小值: {np.min(valid_ddop):.2f}")
        print(f"最大值: {np.max(valid_ddop):.2f}")
        print(f"平均值: {np.mean(valid_ddop):.2f}")
        print(f"Std Dev: {np.std(valid_ddop):.2f}")


def main_sim(sim, week, start_time, end_time, step_s, elev_mask_deg, output_file, 
             tle_file=None, clk_data=None, clk_drift_data=None, receiver_pos=None, 
             prnoise=None, dopnoise=None,
             clk_bias=1000.0, clk_drift=0.05,
             selected_prns=None):
    """
    Main simulation function.
    
    Args:
        sim: SatelliteSimulation instance (already configured with iono/tropo/multipath flags)
        week: GPS week number
        start_time: Start GPS time in seconds
        end_time: End GPS time in seconds
        step_s: Time step in seconds
        elev_mask_deg: Elevation mask in degrees
        output_file: Output file path
        tle_file: TLE file path
        clk_data: Clock bias data array (from file), None to use model
        clk_drift_data: Clock drift data array (from file), None to use model
        receiver_pos: Receiver ECEF position
        prnoise: Pseudorange noise sigma in meters
        dopnoise: Doppler noise sigma in m/s
        clk_bias: Clock bias in meters (used when clk_data is None)
        clk_drift: Clock drift in m/s (used when clk_data is None)
        selected_prns: List of PRN numbers to simulate (None = all satellites)
    """

    # SetGPS week（可according to实际情况修改）
    sim.gps_week = week  # 例如2022年某周
    
    receiver_vel = np.array([0.0, 0.0, 0.0])  # 静态接收机
    
    # receiverclock bias和clock drift - use passed values (from model or will be overridden by file data)
    clock_bias_model = clk_bias
    clock_drift_model = clk_drift

    sim.start_time = start_time
    # satellitedata来源
    satellites = []
    
    # usingTLEfile
    if tle_file and os.path.exists(tle_file):
        # LoadTLEfile
        print(f"加载TLE文件: {tle_file}")
        count = sim.load_tle_file(tle_file)
        if count > 0:
            # usingTLEdata
            for prn, sat_data in sim.tle_satellites.items():
                # Filter by selected PRNs if specified
                if selected_prns is not None and prn not in selected_prns:
                    continue
                satellites.append({
                    'prn': prn,
                    'name': sat_data['name'],
                    'use_tle': True
                })
            if selected_prns is not None:
                print(f"使用TLE模式, selected{len(satellites)}/{count}卫星进行仿真")
            else:
                print(f"使用TLE模式, 已加载{len(satellites)}卫星")
    else:
        print("ERROR! 请输入星历文件!")
    i=0
    # 打开输出file
    with open(output_file, 'w', encoding='utf-8') as f:
        # 在时间range内进line仿真
        gps_time = start_time
        while gps_time <= end_time:
            # 对每satellite进line仿真
            for sat in satellites:
                try:
                    # Calculateclock biasclock drift
                    if clk_data is not None:
                        # Use clock data from file
                        sim.clock = clk_data[i]
                        sim.clockdrift = clk_drift_data[i] if clk_drift_data is not None else 0.0
                    else:
                        # Use model: clock = bias + drift * time
                        sim.clockdrift = clock_drift_model
                        sim.clock = clock_bias_model + sim.clockdrift * (gps_time - sim.start_time)
                    # Calculateobservation
                    if 'use_tle' in sat and sat['use_tle']:
                        # usingTLECalculate
                        observations = sim.compute_observations(
                            sat['prn'],
                            None,  # 使用TLE模式不需要传递ephemeris 
                            receiver_pos,
                            receiver_vel,
                            gps_time,
                            prnoise,
                            dopnoise
                        )
                    else:
                        # usingephemerisCalculate
                        observations = sim.compute_observations(
                            sat['prn'],
                            sat['ephemeris'], 
                            receiver_pos,
                            receiver_vel,
                            gps_time,
                            prnoise,
                            dopnoise
                        )
                    
                    # 判断satellitewhether可见（elevation angle大于0度）
                    if observations['elevation'] > elev_mask_deg:
                        # 将carrier phase从米Convert为周
                        carrier_phase_cycles = observations['carrier_phase'] / LAMBDA_L1
                        
                        # Writedataline
                        f.write(f"OBS: 1  {observations['emission_time']:.9f}    "
                                f"{gps_time:.9f}      "
                                f"{observations['pseudorange']:.2f}    "
                                f"{observations['doppler']:.5f}    "
                                f"{observations['snr']:.1f}   "
                                f"6   1    1   "
                                f"{observations['elevation']:.2f}  "
                                f"{sat['prn']}\n")
                except Exception as e:
                    print(f"计算卫星{sat['prn']}观测值时出错: {e}")
                    continue
            
            # 时间前进
            gps_time += step_s  # 每1秒采样一次, 可根据需要调整
            i = i + 1

    print(f"仿真完成!结果已保存到 {output_file}")

if __name__ == "__main__":
    # Checkinput_data目录whether存在, 不存在则Create
    os.makedirs("input_data", exist_ok=True)
    
    # Set输出filepath
    sim_file = os.path.join(".", "OutputData", "sim_spt_TLE.txt")
    
    clkfilename = os.path.join(".", "InputData", "rck_2025001.txt")
    clk_driftfilename = os.path.join(".", "InputData", "rckdt_2025001.txt")
    rck_values, rck_drift_values = read_rck_data(clkfilename,clk_driftfilename)

    # 仿真data
    tle_file = os.path.join(".", "InputData", "LEO_sat.tle")
    # main(436046, 436046 +1000, sim_file, tle_file)  # 仿真30分钟data

    init_pos = geodetic_to_ecef(39.8*np.pi/180, 142.51*np.pi/180, 60.0)   # elevation=10.0

    week = 2377
    
    # 运line仿真
    main_sim(week,459200, 459200 +1000, sim_file, tle_file,rck_values,rck_drift_values,init_pos)
    
    # # 绘制Dopplerobservation
    plot_doppler_observations(sim_file)
    
    # # 绘制elevationobservation
    # plot_elevation_observations(sim_file,week,init_pos)
    
    # 绘制全球DDOP图
    # ddop_output_file = os.path.join(".", "input_data", "global_ddop.png")
    # plot_global_ddop(week, tle_file, ddop_output_file)