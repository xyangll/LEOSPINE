from mimetypes import init
import numpy as np
from scipy.stats import median_abs_deviation
import math
import matplotlib.pyplot as plt
from scipy.optimize import least_squares, curve_fit
from core.utilities import ecef_to_geodetic, geodetic_to_ecef
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import os
from sgp4.api import Satrec
from sgp4.earth_gravity import wgs84
from datetime import datetime, timedelta
from tqdm import tqdm
import csv
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler

# Constants
SPEED_OF_LIGHT = 299792458.0  # Speed of light (m/s)
OMEGA_EARTH = 7.2921151467e-5  # Earth rotation rate (rad/s)
MU_EARTH = 3.986004418000002e+14 # Earth GM (m^3/s^2)
F_L1 = 1.199169832000000e+09   # SPT frequency (Hz)
LAMBDA_L1 = SPEED_OF_LIGHT / F_L1  # L1 wavelength (m)

tle_satellites = {}  # Store TLE satellite objects
# gps_week=2377
Tle_sat_pos = False   # Use TLE to compute satellite position/velocity
global_sigma={'pr':[],'dop':[]}

def calculate_satellite_position(ephemeris, t):
    """Compute satellite position and velocity"""
    # Constants
    F = -4.442807633e-10  # Relativistic correction constant
    
    # Extract ephemeris parameters
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
    
    # Time difference
    dt_toc = t - toc
    if dt_toc > 302400.0:
        dt_toc = dt_toc - 604800.0
    elif dt_toc < -302400.0:
        dt_toc = dt_toc + 604800.0
    
    # Satellite clock bias
    clock_bias = af0 + af1 * dt_toc + af2 * dt_toc * dt_toc - tgd1
    
    # Correct time with clock bias
    t_corrected = t - clock_bias
    
    # Compute dt
    dt = t_corrected - toe
    if dt > 302400.0:
        dt = dt - 604800.0
    elif dt < -302400.0:
        dt = dt + 604800.0
    
    # Mean motion
    A = sqrtA * sqrtA
    n0 = math.sqrt(MU_EARTH / (A * A * A))
    n = n0 + delta_n
    
    # Mean anomaly
    M = M0 + n * dt
    
    # Solve eccentric anomaly E via Newton-Raphson (matches MATLAB logic)
    E = M
    max_iter = 30
    iter_count = 0
    diff = 1.0  # Init larger than threshold
    
    while abs(diff) > 1.0e-12 and iter_count < max_iter:
        iter_count += 1
        diff = (M - (E - e * math.sin(E))) / (1.0 - e * math.cos(E))
        E = E + diff
    
    # Cos/Sin of eccentric anomaly
    cos_E = math.cos(E)
    sin_E = math.sin(E)
    
    # Relativistic correction on clock bias
    clock_bias = clock_bias + F * e * sqrtA * sin_E
    clock_drift = af1 + 2 * af2 * dt_toc
    
    # True anomaly (use atan2)
    ta = math.atan2(math.sqrt(1.0 - e * e) * sin_E, cos_E - e)
    
    # Argument of latitude
    phi = ta + omega
    
    # Second harmonic corrections
    sin_2phi = math.sin(2.0 * phi)
    cos_2phi = math.cos(2.0 * phi)
    
    # Corrected argument of latitude
    delta_u = Cus * sin_2phi + Cuc * cos_2phi
    u = phi + delta_u
    
    # Corrected radius
    delta_r = Crs * sin_2phi + Crc * cos_2phi
    r = A * (1.0 - e * cos_E) + delta_r
    
    # Corrected inclination
    delta_i = Cis * sin_2phi + Cic * cos_2phi
    i = i0 + delta_i + i_dot * dt
    
    # Position in orbital plane
    x_orb = r * math.cos(u)
    y_orb = r * math.sin(u)
    
    # Corrected RAAN (per MATLAB implementation)
    Omega = Omega0 + (Omega_dot - OMEGA_EARTH) * dt - OMEGA_EARTH * toe
    
    # ECEF position
    x = x_orb * math.cos(Omega) - y_orb * math.cos(i) * math.sin(Omega)
    y = x_orb * math.sin(Omega) + y_orb * math.cos(i) * math.cos(Omega)
    z = y_orb * math.sin(i)
    
    # E-dot
    E_dot = n / (1.0 - e * cos_E)
    
    # u-dot (phi_dot computed correctly)
    phi_dot = math.sqrt(1.0 - e * e) * E_dot / (1.0 - e * cos_E)
    u_dot = phi_dot * (1.0 + 2.0 * (Cus * cos_2phi - Cuc * sin_2phi))
    
    # r-dot
    r_dot = A * e * sin_E * E_dot + 2.0 * (Crs * cos_2phi - Crc * sin_2phi) * phi_dot
    
    # i-dot
    i_dot_corrected = i_dot + 2.0 * (Cis * cos_2phi - Cic * sin_2phi) * phi_dot
    
    # Omega-dot
    Omega_dot_corrected = Omega_dot - OMEGA_EARTH
    
    # Velocity in orbital plane
    x_orb_dot = r_dot * math.cos(u) - r * math.sin(u) * u_dot
    y_orb_dot = r_dot * math.sin(u) + r * math.cos(u) * u_dot
    
    # ECEF velocity (per MATLAB implementation)
    vx = x_orb_dot * math.cos(Omega) - y_orb_dot * math.cos(i) * math.sin(Omega) + z * math.sin(Omega) * i_dot_corrected - y * Omega_dot_corrected
    
    vy = x_orb_dot * math.sin(Omega) + y_orb_dot * math.cos(i) * math.cos(Omega) - z * math.cos(Omega) * i_dot_corrected + x * Omega_dot_corrected
    
    vz = y_orb_dot * math.sin(i) + y_orb * math.cos(i) * i_dot_corrected
    
    return np.array([x, y, z]), np.array([vx, vy, vz])

def teme_to_ecef(teme_position, teme_velocity, dt):
    # Days since J2000 epoch
    j2000_epoch = datetime(2000, 1, 1, 12, 0, 0)
    days_since_j2000 = (dt - j2000_epoch).total_seconds() / 86400.0
    
    # Compute GMST (rad); simplified formula
    gmst = 280.46061837 + 360.98564736629 * days_since_j2000
    gmst = gmst % 360.0  # Wrap to 0-360 deg
    gmst_rad = np.radians(gmst)
    
    # TEME -> ECEF rotation (simplified, only Earth rotation)
    # For high precision, add precession/nutation, etc.
    rotation_matrix = np.array([
        [np.cos(gmst_rad), np.sin(gmst_rad), 0],
        [-np.sin(gmst_rad), np.cos(gmst_rad), 0],
        [0, 0, 1]
    ])
    
    # Transform position
    ecef_position = np.dot(rotation_matrix, teme_position)
    
    # Transform velocity (consider Earth rotation)
    earth_rotation_rate = 7.2921150e-5  # rad/s
    omega_earth = np.array([0, 0, earth_rotation_rate])
    velocity_rotation = np.dot(rotation_matrix, teme_velocity)
    velocity_earth_rotation = np.cross(omega_earth, ecef_position)
    ecef_velocity = velocity_rotation - velocity_earth_rotation
    
    return ecef_position, ecef_velocity

def dbscan_time_series_clustering(time_series_data, eps=0.5, min_samples=5, window_size=10):
    """
    Cluster a 1-D time series with DBSCAN, keeping window centers aligned.
    
    Args:
        time_series_data: 1-D array/list
        eps: DBSCAN epsilon
        min_samples: DBSCAN min_samples
        window_size: sliding window size
    Returns:
        n_clusters: number of clusters (excluding noise)
        labels: label per window (-1 is noise)
        features_scaled: standardized feature matrix
        window_center_indices: center index of each window in original series
    """
    ts_data = np.array(time_series_data).reshape(-1, 1)
    n_samples = len(ts_data)
    features = []
    
    # Window centers aligned to original indices
    window_center_indices = []

    # Edge padding by repeating border values
    pad_width = window_size // 2
    ts_data_padded = np.pad(ts_data, ((pad_width, pad_width), (0, 0)), mode='edge')
    features = []
    window_center_indices = []
    # Traverse padded data so every original point is a window center
    for i in range(n_samples):
        window = ts_data_padded[i:i+window_size].flatten()
        window_mean = np.mean(window)
        window_std = np.std(window)
        features.append([window_mean, window_std])
        # Center index corresponds to original data index
        center_index = i
        window_center_indices.append(center_index)
    
    features = np.array(features)
    window_center_indices = np.array(window_center_indices)

    # Standardize features
    scaler = StandardScaler()
    features_scaled = scaler.fit_transform(features)
    
    # DBSCAN clustering
    dbscan = DBSCAN(eps=eps, min_samples=min_samples)
    dbscan.fit(features_scaled)
    
    labels = dbscan.labels_

    noise_label = None
    if np.any(labels == -1):
        noise_label = labels.max() + 1
        labels = labels.copy()
        labels[labels == -1] = noise_label

    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    
    return n_clusters, labels, features_scaled, window_center_indices


class SingleSatDopplerPositioning:
    """Single LEO satellite Doppler positioning"""
    
    def __init__(self, init_pos=None, init_vel=None, height_constraint=None, height_constraint_mode='hard'):
        """
        Initialize positioning system.
        
        Args:
            init_pos: initial ECEF position [x,y,z] (m); default Earth surface near Beijing
            init_vel: initial velocity [vx,vy,vz] (m/s); default zeros
            height_constraint: optional height constraint (m)
            height_constraint_mode: 'hard', 'soft', or 'both'
        """
        
        # Initial state and covariance
        if init_pos is None:
            # Default start near Beijing
            init_pos = geodetic_to_ecef(39.9*np.pi/180, 116.3*np.pi/180, 50.0)
        
        if init_vel is None:
            init_vel = [0.0, 0.0, 0.0]
        
        # State vector: [x, y, z, vx, vy, vz]
        # Fix type mismatch
        init_pos_list = list(init_pos) if isinstance(init_pos, tuple) else list(init_pos)
        init_vel_list = list(init_vel) if isinstance(init_vel, tuple) else list(init_vel)
        self.state = np.array(init_pos_list + init_vel_list, dtype=float)
        
        # Initial covariance
        pos_var = 1e6**2  # position variance (m^2)
        vel_var = 100**2  # velocity variance (m^2/s^2)
        self.cov = np.diag([pos_var, pos_var, pos_var, vel_var, vel_var, vel_var])
        
        # Noise parameters
        self.sigma_dop = 1.0  # Doppler measurement noise (Hz)
        self.q_pos = 0.01     # Position process noise (m^2/s)
        self.q_vel = 0.01     # Velocity process noise (m^2/s^3)
        
        # Height constraint
        self.height_constraint = height_constraint
        self.height_constraint_mode = height_constraint_mode  # mode for height constraint
        
        # Store trajectories/results
        self.trajectory = []
        self.gdop_history = []
        self.residuals = []
        self.times = []
        
        self.last_adaptive_weight = 1.0  # store last weight
        self.last_height_error = 0.0
        
        # Adaptive weights
        self.use_adaptive_weights = True
        self.residual_threshold = 5.0     # residual threshold (Hz) to downweight
        self.min_weight = 0.1             # minimum weight limit
        self.weight_for_obs = []
        self.iternum = 0
        # Placeholder for deep-learning initializer
        self.dl_initializer = None
        self.use_dl_initializer = False

        self.start_time = 436046
        self.atm_cor = True
        self.ionpara = [
            3.2596E-08, 1.4901E-08, -1.7881E-07, -5.9605E-08,
            1.3722E+05, 1.6384E+04, -3.2768E+05, 3.2768E+05
        ]
        
        # HVCE options
        self.use_hvce = False  # Helmert variance component estimation
        self.hvce = HelmertVCE(n_groups=2)  # two groups: PR, Doppler
        self.hvce_history = []  # HVCE iteration history
        self.pr_weight_factor = 1  # pseudorange weight factor
        self.dop_weight_factor = 1  # Doppler weight factor

        self.rece_time = []

        self.enable_tropo = True
        self.enable_iono = True
        self.enable_multipath = False
        self.enable_relativistic = False
        self.enable_hardware = False
        self.enable_pco = False

        self.DBSCAN_result = False
        self.use_DBSCAN = True
    
    def batch_least_squares_position_only(self, times, week, prnlist, ephemeris, doppler_obs, subset_snr, init_pos=None):
        """Estimate receiver position via batch least squares (Doppler only)"""
        # Get initial state
        current_state = np.zeros(4)
        
        # Optional: grid search init
        # if init_pos is not None:
        #     print("Grid search initialization...")
        #     grid_init_pos = self.grid_search_initialization(times, prnlist, ephemeris, doppler_obs, init_pos)
        #     print(f"Grid search initial position: {grid_init_pos}")
        #     current_state[:3] = grid_init_pos
        # else:
        #     current_state[:3] = self.state[:3]

        if init_pos is not None:
            current_state[:3] = init_pos
        else:
            current_state[:3] = self.state[:3]

        # Initial height info
        init_lat, init_lon, init_h = ecef_to_geodetic(*current_state[:3])
        if self.height_constraint is not None:
            print(f"Initial height: {init_h:.2f} m, target constraint: {self.height_constraint:.2f} m")
        
        # If 'both', keep a backup of initial state
        compare_results = {}

        if self.height_constraint is not None:
            
            # Soft-constrained positioning
            if self.height_constraint_mode in ['soft', 'both']:
                soft_result = self._batch_lsq_with_soft_constraint(
                    times, week, prnlist, ephemeris, doppler_obs, subset_snr,
                    current_state.copy(), init_h
                )
                
                if self.height_constraint_mode == 'soft':
                    return soft_result
                else:
                    compare_results['soft'] = soft_result
        else:
            # Unconstrained positioning
            none_constraint_result = self._batch_lsq_with_soft_constraint(
                    times, week, prnlist, ephemeris, doppler_obs,  subset_snr,
                    current_state.copy(), init_h
                )
            return none_constraint_result

    
    
    def _batch_lsq_with_soft_constraint(self, times, week, prnlist, ephemeris, doppler_obs, subset_snr, initial_state, init_h):
        """Least squares with soft height constraint"""
        if self.height_constraint is not None:
            print("\nSoft-constrained positioning...")
        
        # Multi-stage LSQ, gradually reduce height weight
        # Define stage weights
        height_weights = [1.0]
        current_state = initial_state.copy()
        
        # Store per-stage results
        stage_heights = []
        self.iternum = 0

        for stage, weight in enumerate(height_weights):
            # Residual function with current stage weight
            def residual_func(state):
                position = state[:3]
                clock_drift = state[3]  # clock drift (Hz)
                
                # Doppler residuals
                doppler_residuals,weight_for_doppler = self.compute_residuals_1(state, times, week, prnlist, ephemeris, doppler_obs, subset_snr)
                # Apply height constraint if present
                if self.height_constraint is not None:
                    # Current height
                    lat, lon, h = ecef_to_geodetic(*position)
                    
                    # Height error with stage weight
                    height_error = h - self.height_constraint
                    
                    # Height residual
                    height_residual = height_error * weight
                    
                    # Store last error/weight for debugging
                    self.last_height_error = height_error
                    self.last_adaptive_weight = weight
                    
                    # Append height residual
                    all_residuals = np.append(doppler_residuals,height_residual)
                else:
                    all_residuals = np.array(doppler_residuals)

                # Robust IGGIII
                if self.use_adaptive_weights and self.iternum > 3:
                    weight_IGGIII = self.robust_weighting(doppler_residuals, k0=1.5, k1=3.0)
                    weights = weight_IGGIII * weight_for_doppler    
                    weights = weights  / np.max(weights)         
                else:
                    weights = weight_for_doppler  / np.max(weight_for_doppler)   

                if self.height_constraint is not None: 
                    weights = np.append(weights ,self.last_adaptive_weight)

                self.weight_for_obs = weights
                return all_residuals,weights
            
            # Residuals at current state
            once_residuals,weight_ = residual_func(current_state)
            
    # # Adaptive weights (legacy, disabled)
    # if self.use_adaptive_weights and stage > 0:  # first stage uses uniform weights
            #     doppler_res = current_residuals[:-1] if self.height_constraint is not None else current_residuals
            #     _, doppler_weights = self.compute_residuals_with_weights(
            #         current_state, times,prnlist,ephemeris, doppler_obs
            #     )
                
    #     # Build full weight vector (including height constraint)
            #     if self.height_constraint is not None:
            #         weights = np.ones(len(current_residuals))
            #         weights[:-1] = doppler_weights
    #         weights[-1] = 1.0  # keep height weight unchanged
            #     else:
            #         weights = doppler_weights
            # else:
    #     # Use uniform weights
            #     weights = np.ones_like(current_residuals)
            
            

            def _height_jacobian(position):
                """Analytic Jacobian of height w.r.t. ECEF position."""
                x, y, z = position
                a = 6378137.0  # WGS84 semi-major
                b = 6356752.3142  # WGS84 semi-minor
                e2 = 1 - (b/a)**2  # first eccentricity squared
                
                # Auxiliary params
                p = np.sqrt(x**2 + y**2)
                theta = np.arctan2(z * a, p * b)
                
                # Latitude via parametric form
                sin_theta = np.sin(theta)
                cos_theta = np.cos(theta)
                
                # Geodetic latitude
                lat = np.arctan2(z + e2 * b * sin_theta**3, p - e2 * a * cos_theta**3)
                sin_lat = np.sin(lat)
                cos_lat = np.cos(lat)
                
                # Prime vertical radius
                N = a / np.sqrt(1 - e2 * sin_lat**2)
                
                # Height partials wrt ECEF
                dh_dx = cos_lat * np.cos(np.arctan2(y, x))
                dh_dy = cos_lat * np.sin(np.arctan2(y, x))
                dh_dz = sin_lat
                
                return np.array([dh_dx, dh_dy, dh_dz])
                
            # Jacobian function
            def jacobian_func(state):
                position = state[:3]
                
                # Jacobian for Doppler observations
                n_obs = len(times)
                n_params = 4  # 3 position params + 1 clock drift
                
                # Rows depend on height constraint
                if self.height_constraint is not None:
                    jac = np.zeros((n_obs + 1, n_params))
                else:
                    jac = np.zeros((n_obs, n_params))
                
                for i in range(n_obs):
                    if prnlist[i]==0:
                        continue
                    sat_pos,sat_vel = cal_sat_posvel(times[i], week[i], prnlist[i],ephemeris,position)
                    
                    # Line-of-sight
                    los = position - sat_pos
                    range_ = np.linalg.norm(los)
                    los_unit = los / range_
                    
                    # Assume receiver stationary
                    rx_vel = np.zeros(3)

                    # Position derivatives
                    dlos_dpos = (np.eye(3) - np.outer(los_unit, los_unit)) / range_
                    dDop_dpos = F_L1/SPEED_OF_LIGHT * np.dot(dlos_dpos, (sat_vel - rx_vel))
                    
                    # Fill Jacobian
                    jac[i, :3] = -dDop_dpos
                    jac[i, 3] = -SPEED_OF_LIGHT/ F_L1 # clock drift partial
                
                # Height constraint row
                if self.height_constraint is not None:
                    lat, lon, h = ecef_to_geodetic(*position)
                    
                    # Height partials
                    dh_dpos = _height_jacobian(position.copy())
                    
                    # Apply stage weight
                    jac[-1, :3] = dh_dpos * weight
                    jac[-1, 3] = 0.0  # height constraint has no clock drift partial
                
                return jac
            
            # Weighted residual/Jacobian
            def weighted_residual_func(state):
                res,weights = residual_func(state)
                return res * np.sqrt(weights)
                #return res
            
            def weighted_jacobian_func(state):
                jac = jacobian_func(state)
                for i in range(len(self.weight_for_obs)):
                    jac[i, :] *= np.sqrt(self.weight_for_obs[i])
                self.iternum = self.iternum + 1
                return jac
            
            # Run optimization (example)
            # result = least_squares(
            #     weighted_residual_func, 
            #     current_state, 
            #     jac=weighted_jacobian_func, 
            #     method='trf',
            #     loss='linear',  # linear loss since custom weights are applied
            #     f_scale=1.0,
            #     verbose=0
            # )
            result = least_squares(
                weighted_residual_func, 
                current_state,  
                jac=weighted_jacobian_func, 
                method='lm',
                loss='linear',
                ftol=1e-6,      # relaxed convergence tolerance
                xtol=1e-6,      # relaxed convergence tolerance
                gtol=1e-6,      # relaxed convergence tolerance
                max_nfev=20,    # limit function evaluations
                verbose=0       # set 1 to inspect iterations
            )
            # result = least_squares(
            #     weighted_residual_func, 
            #     current_state,  
            #     jac=weighted_jacobian_func, 
            #     method='trf',
            #     loss='linear',
            #     ftol=1e-6,      # relaxed tolerance
            #     xtol=1e-6,      # relaxed tolerance
            #     gtol=1e-6,      # relaxed tolerance
            # max_nfev=20,    # 限制函数评估次数
            # tr_solver='exact',  # using精确求解器
            # tr_options={'initial_trust_radius': 1.0},  # Set较小initial信任区域半径
            # verbose=0       # 设为1可以查看iteration过程, 帮助调试
            # )
            # Updatecurrent状态为本阶段优化结果
            current_state = result.x
            
            # 记录currentaltitude
            if self.height_constraint is not None:
                lat, lon, h = ecef_to_geodetic(*current_state[:3])
                stage_heights.append(h)
                print(f"阶段 {stage+1}/{len(height_weights)} (权重={weight:.1f}): 高度 = {h:.2f}m, 偏差 = {h-self.height_constraint:.2f}m")
        
        # 提取final结果
        estimated_state = current_state
        estimated_position = estimated_state[:3]
        estimated_clock_drift = estimated_state[3]
        
        # Get仅Doppler残差(不包括可能altitude约束残差)
        all_residuals,__ = residual_func(estimated_state)
        if self.height_constraint is not None:
            doppler_residuals = all_residuals[:-1]  # 移除高度约束残差
        else:
            doppler_residuals = all_residuals
        
        # Calculate协方差matrix - using最后一次iteration雅可比matrix
        J = jacobian_func(estimated_state)
        if self.height_constraint is not None:
            J_dop = J[:-1, :]  # 仅使用多普勒部分计算协方差
        else:
            J_dop = J

        # using残差平方和除以自由度作为方差估计
        sigma2 = np.sum(doppler_residuals**2) / len(doppler_residuals)
        
        try:
            cov = sigma2 * np.linalg.inv(J_dop.T @ J_dop)
        except np.linalg.LinAlgError:
            print("警告: 协方差计算中出现奇异矩阵, 使用伪逆")
            cov = sigma2 * np.linalg.pinv(J.T @ J)
        
        # Update状态
        self.state[:3] = estimated_position
        
        # 输出finalaltitude信息
        lat, lon, h = ecef_to_geodetic(*estimated_position)
        if self.height_constraint is not None:
            print(f"软约束定位结果 - 高度: {h:.2f}m, 目标: {self.height_constraint:.2f}m, 残差RMS: {np.sqrt(np.mean(np.square(doppler_residuals))):.3f}Hz")
        
        return estimated_position, estimated_clock_drift, cov, doppler_residuals,once_residuals

    def adjust_weights_based_on_snr(self, weights, snr_values):
        """
        according tosignal-to-noise ratio调整权重（optional功能）
        
        Args:
            weights: based on残差initial权重
            snr_values: 对应signal-to-noise ratio值（dB）
        
        Returns:
            调整后权重
        """
        if snr_values is None:
            return weights
            
        # 定义SNR阈值
        snr_min = 45.0  # dB
        snr_max = 70.0  # dB
        
        # 归一化SNR因子 (0.5-1.0)
        snr_factors = np.ones_like(snr_values)
        mask = snr_values < snr_max
        snr_factors[mask] = 0.5 + 0.5 * (snr_values[mask] - snr_min) / (snr_max - snr_min)
        snr_factors = np.clip(snr_factors, 0.5, 1.0)
        
        # 结合两种权重
        combined_weights = weights * snr_factors
        
        return combined_weights
    
    def elevation_snr_weighting(self, elevation_deg, snr_db, a=0.01, b=1, gamma=0.1):
        """elevation与SNR组合定权函数
        :param elevation_deg: satelliteelevation数组(度)
        :param snr_db: signal-to-noise ratio数组(dB-Hz)
        :param a,b: elevation模型Args
        :param gamma: SNR权重系数
        :return: 归一化权重数组
        """
        # elevation权重(正弦模型)
        ele_rad = np.deg2rad(np.clip(elevation_deg, 1, 90))
        sin_e = np.sin(ele_rad)
        w_ele = 1 / (a**2 + b**2 / (sin_e**2 + 1e-6))
        
        # SNR权重(指数衰减模型)
        # w_snr = 1 / (1 + gamma * 10**(-snr_db / 10))
        w_snr = 1 / ( 1+ gamma * 10**(-snr_db / 10))

        if elevation_deg < 7.0:
            w_ele = 0
        else:
            w_ele = 1
        # 组合权重
        combined_weights = w_ele * w_snr
        return w_snr

    def robust_weighting(self, residuals, k0=1.5, k1=3.0):
        """
        抗差定权函数(改进IGGIII模型)
        :param residuals: 观测残差数组
        :param k0: 轻微异常阈值(default1.5σ)
        :param k1: 严重异常阈值(default3.0σ)
        :return: 抗差权重数组[0-1]
        """
        abs_res = np.abs(residuals)
        sigma = median_abs_deviation(residuals, scale='normal')  # 鲁棒Std Dev估计
        scaled_res = abs_res / (sigma + 1e-6)  # 防止除零
        
        weights = np.ones_like(scaled_res)
        mask = (scaled_res > k0) & (scaled_res <= k1)
        weights[mask] = (k0 / scaled_res[mask]) * ((k1 - scaled_res[mask])/(k1 - k0))**2
        weights[scaled_res > k1] = 0
        return weights

    def compute_residuals_with_weights(self, state, times, prnlist, eph, doppler_obs, snr_values=None):
        """
        CalculateDoppler残差并Generate融合SNR自适应权重
        
        Args:
            state: current状态vector [x, y, z, b_dot]
            times: 观测时间column表
            prnlist: satellitePRN号column表
            eph: satelliteephemeris字典
            doppler_obs: Dopplerobservationcolumn表
            snr_values: signal-to-noise ratio值column表（optional）
        
        Returns:
            residuals: Doppler残差
            weights: 对应权重值
        """
        # Calculate标准残差
        residuals = self.compute_residuals(state, times, prnlist, eph, doppler_obs, snr_values)
        
        # if未enable自适应权重, Returns全1权重
        if not self.use_adaptive_weights:
            return residuals, np.ones_like(residuals)
        
        # based on残差Calculateinitial权重
        if len(residuals) > 3:
            # usingMAD作为鲁棒尺度估计
            mad = np.median(np.abs(residuals - np.median(residuals)))
            scale = 1.4826 * mad  # 使MAD成为正态分布Std Dev无偏估计
        else:
            # 样本太少时usingStd Dev
            scale = max(np.std(residuals), 0.01)
        
        # 标准化残差
        normalized_residuals = np.abs(residuals) / (scale * self.residual_threshold)
        
        # usingHuber型权重函数
        weights = np.ones_like(residuals)
        large_res_mask = normalized_residuals > 1.0
        
        if np.any(large_res_mask):
            weights[large_res_mask] = 1.0 / normalized_residuals[large_res_mask]
            
            # 确保权重不低于minimum值
            weights = np.maximum(weights, self.min_weight)
        
        # if提供了SNR值, 进一步调整权重
        if snr_values is not None and len(snr_values) == len(weights):
            weights = self.adjust_weights_based_on_snr(weights, snr_values)
            
            if np.any(large_res_mask):
                print(f"检测到{np.sum(large_res_mask)}异常值, 权重range: {np.min(weights):.3f}-{np.max(weights):.3f}")
        
        return residuals, weights

    def compute_residuals(self, state, times, week, prnlist, eph, doppler_obs, snr_values=None):
        """
        CalculateDoppler残差, consideringall误差源
        """
        position = state[:3]
        clock_drift = state[4]
        lat, lon, h = ecef_to_geodetic(*position)
        
        residuals = []
        weights_for_doppler = []
        for i in range(len(times)):
            # Calculatesatelliteposition和velocity
            if prnlist[i]==0:
                continue
            sat_pos, sat_vel = cal_sat_posvel(times[i], week[i], prnlist[i], eph, position)
            los = position - sat_pos
            range_ = np.linalg.norm(los)
            los_unit = los / range_
            
            elevation,__ = calculate_elevation_angle(sat_pos, position)

            # 1. Calculateatmospheric effect
            if self.enable_tropo:
                trop_delay = self.compute_tropospheric_delay_rate(times[i], week[i], prnlist[i], eph, position, h, 0.005)
            else:
                trop_delay = 0.0

            if self.enable_iono:
                iono_delay = self.compute_ionospheric_delay_rate(times[i], week[i], prnlist[i], eph, position, h, 0.1)
            else:
                iono_delay = 0.0

            # Calculatetotal大气delayrate of change
            total_delay_rate = trop_delay + iono_delay
            
            # Convert为Doppler频移
            atmos_doppler = -total_delay_rate * F_L1/SPEED_OF_LIGHT

            # 2. CalculateEarth rotation效应
            rotation_effect = self.compute_earth_rotation_effect(position, sat_pos, sat_vel)
            rotation_doppler = -rotation_effect * F_L1/SPEED_OF_LIGHT
            
            # 3. Calculatemultipath效应
            if self.enable_multipath:
                multipath_effect = self.compute_multipath_effect(
                    elevation, 
                    snr_values[i] if snr_values is not None else None
                )
            else:
                multipath_effect = 0.0
            
            # 4. Calculatephase center效应
            #phase_center_effect = self.compute_phase_center_effect(elevation, azimuth)
            if self.enable_pco:
                pco_effect = 0.0
            else:
                pco_effect = 0.0

            # 5. Calculate几何Doppler
            geometric_doppler = F_L1/SPEED_OF_LIGHT * np.dot(sat_vel, los_unit)
            
            # 6. Calculaterelativistic效应
            if self.enable_relativistic:
                relativistic_effect = self.compute_relativistic_effect(sat_pos, sat_vel, position)
                relativistic_effect = -relativistic_effect * F_L1/SPEED_OF_LIGHT
            else:
                relativistic_effect = 0.0

            # 7. Calculatehardware delay效应
            if self.enable_hardware:
                hardware_delay = 0.0  # 默认值, 实际应从外部获取或建模
            else:
                hardware_delay = 0.0
            
            # 合成total预测Doppler
            pred_doppler = (geometric_doppler + 
                           atmos_doppler +      # 大气效应
                           rotation_doppler +    # 地球自转
                           multipath_effect +    # 多路径
                           pco_effect + # 相位中心
                           relativistic_effect + # 相对论
                           hardware_delay +     # 硬件延迟
                           clock_drift)         # 钟漂
            
            # Calculate残差
            residuals.append(doppler_obs[i] - pred_doppler)
            # Calculate权重
            weight_temp=self.elevation_snr_weighting(elevation, snr_values[i], a=0.003, b=0.01, gamma=0.01)
            weights_for_doppler.append(weight_temp)

            # 输出各种效应大小(for调试)
            # if i == 0:  # 只打印第一历元结果
            # print(f"\n第一历元各效应大小:")
            # print(f"几何Doppler: {geometric_doppler:.3f} Hz")
            # print(f"atmospheric effect: {atmos_doppler:.3f} Hz")
            # print(f"Earth rotation: {rotation_doppler:.3f} Hz")
            # print(f"multipath效应: {multipath_effect:.3f} Hz")
            # print(f"phase center效应: {phase_center_effect:.3f} Hz")
            # print(f"relativistic效应: {relativistic_effect:.3f} Hz")
            # print(f"hardware delay效应: {hardware_effect:.3f} Hz")
            # print(f"total预测Doppler: {pred_doppler:.3f} Hz")
            # print(f"观测Doppler: {doppler_obs[i]:.3f} Hz")
            # print(f"残差: {residuals[-1]:.3f} Hz")
        
        return np.array(residuals),np.array(weights_for_doppler)
    
    def compute_residuals_1(self, state, times, week, prnlist, eph, doppler_obs, snr_values=None):
        """
        CalculateDoppler残差, consideringall误差源
        """
        position = state[:3]
        clock_drift = state[3]
        lat, lon, h = ecef_to_geodetic(*position)
        
        residuals = []
        weights_for_doppler = []
        for i in range(len(times)):
            # Calculatesatelliteposition和velocity
            if prnlist[i]==0:
                continue
            sat_pos, sat_vel = cal_sat_posvel(times[i], week[i], prnlist[i], eph, position)
            los = position - sat_pos
            range_ = np.linalg.norm(los)
            los_unit = los / range_
            
            elevation,__ = calculate_elevation_angle(sat_pos, position)

            # 1. Calculateatmospheric effect
            if self.enable_tropo:
                trop_delay = self.compute_tropospheric_delay_rate(times[i],week[i], prnlist[i], eph, position, h, 0.005)
            else:
                trop_delay = 0.0

            if self.enable_iono:
                iono_delay = self.compute_ionospheric_delay_rate(times[i],week[i], prnlist[i], eph, position, h, 0.1)
            else:
                iono_delay = 0.0

            # Calculatetotal大气delayrate of change
            total_delay_rate = trop_delay + iono_delay
            
            # Convert为Doppler频移
            atmos_doppler = -total_delay_rate * F_L1/SPEED_OF_LIGHT

            # 2. CalculateEarth rotation效应
            rotation_effect = self.compute_earth_rotation_effect(position, sat_pos, sat_vel)
            rotation_doppler = -rotation_effect * F_L1/SPEED_OF_LIGHT
            
            # 3. Calculatemultipath效应
            if self.enable_multipath:
                multipath_effect = self.compute_multipath_effect(
                    elevation, 
                    snr_values[i] if snr_values is not None else None
                )
            else:
                multipath_effect = 0.0
            
            # 4. Calculatephase center效应
            #phase_center_effect = self.compute_phase_center_effect(elevation, azimuth)
            if self.enable_pco:
                pco_effect = 0.0
            else:
                pco_effect = 0.0

            # 5. Calculate几何Doppler
            geometric_doppler = F_L1/SPEED_OF_LIGHT * np.dot(sat_vel, los_unit)
            
            # 6. Calculaterelativistic效应
            if self.enable_relativistic:
                relativistic_effect = self.compute_relativistic_effect(sat_pos, sat_vel, position)
                relativistic_effect = -relativistic_effect * F_L1/SPEED_OF_LIGHT
            else:
                relativistic_effect = 0.0

            # 7. Calculatehardware delay效应
            if self.enable_hardware:
                hardware_delay = 0.0  # 默认值, 实际应从外部获取或建模
            else:
                hardware_delay = 0.0
            
            # 合成total预测Doppler
            pred_doppler = (geometric_doppler + 
                           atmos_doppler +      # 大气效应
                           rotation_doppler +    # 地球自转
                           multipath_effect +    # 多路径
                           pco_effect + # 相位中心
                           relativistic_effect + # 相对论
                           hardware_delay +     # 硬件延迟
                           clock_drift)         # 钟漂
            
            # Calculate残差
            residuals.append(doppler_obs[i] - pred_doppler)
            # Calculate权重
            weight_temp=self.elevation_snr_weighting(elevation, snr_values[i], a=0.003, b=0.01, gamma=0.01)
            weights_for_doppler.append(weight_temp)

            # 输出各种效应大小(for调试)
            # if i == 0:  # 只打印第一历元结果
            # print(f"\n第一历元各效应大小:")
            # print(f"几何Doppler: {geometric_doppler:.3f} Hz")
            # print(f"atmospheric effect: {atmos_doppler:.3f} Hz")
            # print(f"Earth rotation: {rotation_doppler:.3f} Hz")
            # print(f"multipath效应: {multipath_effect:.3f} Hz")
            # print(f"phase center效应: {phase_center_effect:.3f} Hz")
            # print(f"relativistic效应: {relativistic_effect:.3f} Hz")
            # print(f"hardware delay效应: {hardware_effect:.3f} Hz")
            # print(f"total预测Doppler: {pred_doppler:.3f} Hz")
            # print(f"观测Doppler: {doppler_obs[i]:.3f} Hz")
            # print(f"残差: {residuals[-1]:.3f} Hz")
        
        return np.array(residuals),np.array(weights_for_doppler)

    def compute_tropospheric_delay_rate(self, time, week, prn, eph, recpos, height, delt):
        """
        Calculatetropospheredelayrate of change
        """
        
        time_former = time - delt
        time_latter = time + delt
        rec_lat, rec_lon, rec_alt = ecef_to_geodetic(*recpos)

        sat_pos_former, sat_vel_former = cal_sat_posvel(time_former, week, prn, eph, recpos)
        sat_pos_latter, sat_vel_latter = cal_sat_posvel(time_latter, week, prn, eph, recpos)

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
    
    def compute_ionospheric_delay_rate(self, time, week, prn, eph, recpos, height, delt):
        """
        Calculateionospheredelayrate of change
        """
        
        time_former = time - delt
        time_latter = time + delt
        rec_lat, rec_lon, rec_alt = ecef_to_geodetic(*recpos)

        sat_pos_former, sat_vel_former = cal_sat_posvel(time_former, week, prn, eph, recpos)
        sat_pos_latter, sat_vel_latter = cal_sat_posvel(time_latter, week, prn, eph, recpos)

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

    def configure_adaptive_weights(self, enable=True, threshold=5.0, min_weight=0.1):
        """
        配置自适应权重Args
        
        Args:
            enable: whetherenable自适应权重
            threshold: 残差阈值(Hz), 超过此值观测权重将降低
            min_weight: minimum权重限制(0-1)
        """
        self.use_adaptive_weights = enable
        self.residual_threshold = threshold
        self.min_weight = min_weight
        
        print(f"自适应权重配置更新: 启用={enable}, 阈值={threshold}Hz, 最小权重={min_weight}")
        
        return self

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

    # def compute_phase_center_effect(self, elevation, azimuth):
    #     """
    # Calculate天线phase center变化影响
    #     """
    # # 简化模型, 实际应该using天线校准data
    #     elevation_rad = np.deg2rad(elevation)
    #     azimuth_rad = np.deg2rad(azimuth)
        
    # # PCO (Phase Center Offset) - Example值
    # pco = np.array([0.001, 0.001, 0.010])  # 米
        
    # # PCV (Phase Center Variations) - 简化模型
    #     pcv = 0.002 * np.cos(elevation_rad) * np.cos(azimuth_rad)
        
    # # Calculatetotalphase center效应
    #     total_effect = np.dot(pco, [np.cos(elevation_rad)*np.cos(azimuth_rad),
    #                                np.cos(elevation_rad)*np.sin(azimuth_rad),
    #                                np.sin(elevation_rad)]) + pcv
        
    #     return total_effect * F_L1/SPEED_OF_LIGHT

    # def compute_hardware_delay_effect(self):
    #     """
    # Calculatehardware delay变化影响
    #     """
    # # 基础hardware delay
    # base_delay = 1e-9  # 1纳秒
        
    # # if有温度信息, considering温度影响
    #     if self.temperature is not None:
    # # 假设标称温度为25度
    # temp_effect = 1e-11 * (self.temperature - 25)  # 每度温差带来额外delay
    #         base_delay += temp_effect
        
    #     return base_delay * F_L1

    def compute_atmospheric_effects(self, time, week, prn, eph, recpos, height, delt):
        """
        Calculateatmospheric effecttotal影响
        """
        # Calculatetropospheredelay
        trop_delay = self.compute_tropospheric_delay_rate(time, week, prn, eph, recpos, height, delt)
        
        # Calculateionospheredelay
        iono_delay = self.compute_ionospheric_delay_rate(time, week, prn, eph, recpos, height, 0.1)
        iono_delay = 0.0
        # Calculatetotal大气delayrate of change
        total_delay_rate = trop_delay + iono_delay
        
        # Convert为Doppler频移
        atmos_doppler = -total_delay_rate * F_L1/SPEED_OF_LIGHT
        
        return atmos_doppler, trop_delay, iono_delay


    def compute_residuals_pseudorange(self, state, times, week, prnlist, eph, pseudorange_obs, snr_values=None):
        """
        Calculatepseudorange残差, consideringall误差源
        """
        position = state[:3]
        clock_bias = state[3]  # 钟差 (m)
        clock_drift = state[4] if len(state) > 4 else 0.0  # 钟漂 (m/秒)
        lat, lon, h = ecef_to_geodetic(*position)
        
        pr_residuals = []
        weights_for_pr = []
        for i in range(len(times)):
            # Calculatesatelliteposition和velocity
            if prnlist[i]==0:
                continue
            sat_pos, sat_vel = cal_sat_posvel(times[i], week[i], prnlist[i], eph, position)
            los = position - sat_pos
            range_ = np.linalg.norm(los)
            los_unit = los / range_
            
            # Calculateelevation和azimuth
            azel = calculate_elevation_angle(sat_pos, position)
            rec_lat, rec_lon, rec_alt = ecef_to_geodetic(*position)   

            # 1. Calculateatmospheric effect
            if self.enable_tropo:                 
                # tropospheredelay
                trop_delay = self.tropmodel([rec_lat, rec_lon, rec_alt], azel[0]* (np.pi / 180.0), 0.1)    
            else:
                trop_delay = 0.0   

            if self.enable_iono:          
                # ionospheredelay
                iono_delay = self.ionmodel(times[i], self.ionpara, [rec_lat, rec_lon, rec_alt], np.deg2rad(azel))
            else:
                iono_delay = 0.0
            
            # 3. Sagnac effect(Earth rotation引起信号传播delay)
            omega_e = OMEGA_EARTH  # 地球自转角速度
            rotation_effect = (omega_e * (sat_pos[0] * position[1] - sat_pos[1] * position[0])) / SPEED_OF_LIGHT

            # 4. multipath效应
            if self.enable_multipath:              
                if snr_values is not None:
                    # based onelevation和signal-to-noise ratio简单multipath模型
                    elevation_rad = np.deg2rad(azel[0])
                    snr = snr_values[i]
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
                r_rec = np.linalg.norm(position)
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
            clock_bias_at_t = clock_bias
            if len(state) > 4:  # 如果状态向量包含钟漂
                clock_bias_at_t += clock_drift * (self.rece_time[i] - self.start_time)
                # clock_bias_at_t += clock_drift * (times[i] - self.start_time)
            sat_clock_bias = 0

            # Calculate预测pseudorange(includingall误差源)
            pred_pseudorange = (
                range_ +                    # 几何距离
                (clock_bias_at_t - sat_clock_bias) +  # 钟差(接收机-卫星)
                trop_delay +                # 对流层延迟
                iono_delay +                # 电离层延迟
                rotation_effect +           # 地球自转效应
                multipath_effect +          # 多路径效应
                relativity_extra +          # 附加相对论效应
                hardware_delay +            # 硬件延迟
                pco_effect                  # 相位中心偏移
            )
            
            # Calculate残差
            pr_residuals.append(pseudorange_obs[i] - pred_pseudorange)
            
            # Calculate权重
            weight_temp = self.elevation_snr_weighting(azel[0], snr_values[i] if snr_values is not None else 45.0, 
                                                     a=0.003, b=0.01, gamma=0.01)
            weights_for_pr.append(weight_temp)     
        
        return np.array(pr_residuals), np.array(weights_for_pr)

    def batch_least_squares_with_pseudorange(self, times,week, prnlist, ephemeris, pseudorange_obs, doppler_obs=None, subset_snr=None, init_pos=None, use_doppler=False):
        """usingpseudorange(optional配合Doppler)估计receiverposition、clock bias和clock drift, usingHVCE优化权重"""
        
        # Getinitial状态
        current_state = np.zeros(5)  # [x,y,z,钟差,钟漂]

        # using提供initialposition或current状态
        if init_pos is not None:
            current_state[:3] = init_pos
        else:
            current_state[:3] = self.state[:3]
        
        # using网格搜索找到更好initialposition
        # if init_pos is not None:
        # print("执line网格搜索Initialize...")
        #     grid_init_pos = self.grid_search_initialization(times, prnlist, ephemeris, doppler_obs, init_pos)
        # print(f"网格搜索找到initialposition: {grid_init_pos}")
        #     current_state[:3] = grid_init_pos
        # else:
        #     current_state[:3] = self.state[:3]

        # initialclock bias和clock drift设为零
        current_state[3] = 0.0  # 钟差 (秒)
        current_state[4] = 0.0  # 钟漂 (秒/秒)
        
        # initialaltitude信息
        init_lat, init_lon, init_h = ecef_to_geodetic(*current_state[:3])
        
        print("\n执行伪距定位...")
        
        # whetherusingaltitude约束
        use_height_constraint = (self.height_constraint is not None)
        
        # Initialize可变维度权重matrix
        # 首次iteration: 2维权重matrix [pseudorange权重, Doppler权重]
        if use_doppler and doppler_obs is not None:
            weights_BDSCAN = np.ones(len(pseudorange_obs) + len(doppler_obs))
        else:
            weights_BDSCAN = np.ones(len(pseudorange_obs))
               
        
        # 重置HVCE对象
        if self.use_hvce and use_doppler:
            self.hvce.reset()
            self.hvce_history = []
            self.pr_weight_factor = 1.0
            self.dop_weight_factor = 1.0
            
        # 开始iteration
        max_iterations = 10
        convergence_threshold = 1e-8
        prev_state = current_state.copy()
        
        for iter_num in range(max_iterations):
            # 定义目标函数(残差函数)
            def residual_func(state):
                position = state[:3]
                
                # Calculatepseudorange残差
                pr_residuals, pr_weights_detailed = self.compute_residuals_pseudorange(
                    state, times,week, prnlist, ephemeris, pseudorange_obs, subset_snr)
                         
                # if同时usingDopplerobservation
                if use_doppler and doppler_obs is not None:
                    # CalculateDoppler观测残差
                    doppler_residuals, dop_weights_detailed = self.compute_residuals(
                        state, times,week, prnlist, ephemeris, doppler_obs, subset_snr)

                    # for i in range(len(pr_weights)):
                    # # 应用pseudorange权重因子
                    #     pr_residuals[i] *= np.sqrt(pr_weights[i])
                    # for i in range(len(dop_weights)):
                    # # 应用Doppler权重因子
                    #     doppler_residuals[i] *= np.sqrt(dop_weights[i])                  
                    # 合并pseudorange和Doppler残差, 应用HVCE权重因子
                    

                    all_residuals = np.concatenate([pr_residuals, doppler_residuals])
                else:
                    all_residuals = pr_residuals
              
                for i in range(len(weights_BDSCAN)):
                    all_residuals[i] *= np.sqrt(weights_BDSCAN[i])

                # if有altitude约束
                if use_height_constraint:
                    lat, lon, h = ecef_to_geodetic(*position)
                    height_error = h - self.height_constraint
                    
                    # 添加altitude约束残差
                    all_residuals = np.append(all_residuals, height_error * 1000.0)  # 权重系数10.0
                
                return all_residuals
            

            def _height_jacobian(position):
                """Calculatealtitude对ECEFposition解析derivative"""
                x, y, z = position
                a = 6378137.0  # WGS84椭球长半轴
                e_sq = 6.69437999014e-3  # 第一偏心率平方
                
                p = np.sqrt(x**2 + y**2)
                theta = np.arctan2(z * a, p * (a * (1 - e_sq)))
                
                sin_theta = np.sin(theta)
                cos_theta = np.cos(theta)
                
                dx = cos_theta * x / p
                dy = cos_theta * y / p
                dz = sin_theta * (a / (a * (1 - e_sq)))
                
                return np.array([dx, dy, dz])
            
            # 定义雅可比matrix函数
            def jacobian_func(state):
                position = state[:3]
                
                # CalculateArgs数量
                n_params = 5  # [x,y,z,钟差,钟漂]
                
                # Calculateobservation数量
                n_pr_obs = sum(1 for prn in prnlist if prn != 0)
                n_dop_obs = n_pr_obs if (use_doppler and doppler_obs is not None) else 0
                n_constraints = 1 if use_height_constraint else 0
                
                # Create雅可比matrix
                jac = np.zeros((n_pr_obs + n_dop_obs + n_constraints, n_params))
                
                # 填充pseudorange观测雅可比matrix
                pr_idx = 0
                for i in range(len(times)):
                    if prnlist[i] == 0:
                        continue
                    
                    # Calculatesatelliteposition和velocity
                    sat_pos, sat_vel = cal_sat_posvel(times[i], week[i], prnlist[i], ephemeris, position)
                    
                    # Calculateline of sightvector
                    los = position - sat_pos
                    range_ = np.linalg.norm(los)
                    los_unit = los / range_
                    
                    # positionArgs偏derivative
                    jac[pr_idx, :3] = -los_unit
                    
                    # clock biasArgs偏derivative
                    jac[pr_idx, 3] = -1
                    
                    # clock driftArgs偏derivative (considering时间差)
                    jac[pr_idx, 4] = -(self.rece_time[i] - self.start_time)
                    
                    pr_idx += 1
                
                # ifusingDopplerobservation, 填充Doppler部分雅可比matrix
                if use_doppler and doppler_obs is not None:
                    dop_idx = n_pr_obs
                    for i in range(len(times)):
                        if prnlist[i] == 0:
                            continue
                        
                        sat_pos, sat_vel = cal_sat_posvel(times[i], week[i], prnlist[i], ephemeris, position)
                        
                        # Calculateline of sightvector
                        los = position - sat_pos
                        range_ = np.linalg.norm(los)
                        los_unit = los / range_
                        
                        # Calculatepositionderivative (Doppler对position偏derivative)
                        dlos_dpos = (np.eye(3) - np.outer(los_unit, los_unit)) / range_
                        dDop_dpos = F_L1/SPEED_OF_LIGHT * np.dot(dlos_dpos, sat_vel)

                        # 填充雅可比matrix
                        jac[dop_idx, :3] = -dDop_dpos
                        jac[dop_idx, 3] = 0.0  # 钟差对多普勒没有影响
                        jac[dop_idx, 4] = -SPEED_OF_LIGHT/ F_L1 # 钟漂偏导数
                        
                        dop_idx += 1
                
                # if有altitude约束, Calculatealtitude对ECEFpositionderivative
                if use_height_constraint:
                    constraint_idx = n_pr_obs + n_dop_obs
                    lat, lon, h = ecef_to_geodetic(*position)
                    
                    # Calculatealtitude对ECEFposition偏derivative (using之前定义辅助函数)
                    dh_dpos = _height_jacobian(position.copy())
                    
                    # 填充altitude约束部分雅可比matrix
                    jac[constraint_idx, :3] = dh_dpos * 1000.0  # 与残差中相同权重
                    jac[constraint_idx, 3:] = 0.0  # 高度约束对钟差和钟漂没有影响
                
                return jac
            
            def weighted_jacobian_func(state):
                jac = jacobian_func(state)
                # 应用动态权重因子                
                for i in range(len(weights_BDSCAN)):
                    jac[i, :] *= np.sqrt(weights_BDSCAN[i])

                return jac

            # usingminimum二乘法求解
            # result = least_squares(
            #     residual_func, 
            #     current_state, 
            #     jac=weighted_jacobian_func, 
            #     method='trf',
            #     loss='linear',
            #     verbose=0
            # )
            result = least_squares(
                residual_func, 
                current_state,  
                jac=weighted_jacobian_func, 
                method='trf',
                loss='linear',
                ftol=1e-6,      # 放宽收敛容差
                xtol=1e-6,      # 放宽收敛容差
                gtol=1e-6,      # 放宽收敛容差
                max_nfev=20,    # 限制函数评估次数
                tr_solver='exact',  # 使用精确求解器
                tr_options={'initial_trust_radius': 1.0},  # 设置较小初始信任区域半径
                verbose=0       # 设为1可以查看迭代过程, 帮助调试
            )
            # 提取结果
            estimated_state = result.x
            
            # Calculatepseudorange和Doppler残差
            position = estimated_state[:3]
            pr_residuals, pr_weights_detailed = self.compute_residuals_pseudorange(
                estimated_state, times, week, prnlist, ephemeris, pseudorange_obs, subset_snr)
            
            if use_doppler and doppler_obs is not None:
                dop_residuals, dop_weights_detailed = self.compute_residuals(
                    estimated_state, times, week, prnlist, ephemeris, doppler_obs, subset_snr)

                all_residuals = np.concatenate([pr_residuals, dop_residuals])
            else:
                all_residuals = pr_residuals

            # 将pseudorange与Doppler残差顺序Write同一file（先pseudorange, 后Doppler）
            # try:
            #     output_path = "residuals_output.txt"
            #     with open(output_path, "w", encoding="utf-8") as f:
            # # Writepseudorange残差
            #         for v in np.ravel(pr_residuals):
            #             f.write(f"{v}\n")
            # # 若存在Doppler残差, 则继续追加Write
            #         if use_doppler and doppler_obs is not None:
            #             for v in np.ravel(dop_residuals):
            #                 f.write(f"{v}\n")
            # print(f"残差已Write {output_path}")
            # except Exception as e:
            # print(f"Write残差filefailed: {e}")

            if self.use_DBSCAN and not self.DBSCAN_result:
                # 应用DBSCAN聚类
                n_clusters, labels, features, window_center_indices  = dbscan_time_series_clustering(
                    np.concatenate([pr_residuals, dop_residuals]), 
                    eps=0.2, min_samples=20, window_size=20
                )
                self.DBSCAN_result = True
                self.hvce.n_groups = n_clusters
                self.hvce.reset()

                unique_labels = set(labels)
                indices_DBSCAN = {}
                for k in unique_labels:              
                    class_member_mask = (labels == k)
                    x_coords = window_center_indices[class_member_mask] 
                    # 将current簇索引按line存储到indices_DBSCAN中
                    indices_DBSCAN[k] = x_coords.tolist()  # 转换为列表格式

            # usingHVCE优化权重
            if self.use_hvce and use_doppler:
                # Get雅可比matrix
                J = jacobian_func(estimated_state)

                # using稳定簇键顺序, 将分组Convert为column表传入 HVCE
                cluster_keys = list(indices_DBSCAN.keys())
                residuals_groups_list = [all_residuals[indices_DBSCAN[k]] for k in cluster_keys]
                jacobian_groups_list = [J[indices_DBSCAN[k]] for k in cluster_keys]
                weights_groups_list = [weights_BDSCAN[indices_DBSCAN[k]] for k in cluster_keys]

                # 执line一次HVCEiteration（column表输入）
                new_weights_list, hvce_converged, sigma2_values = self.hvce.estimate(
                    residuals_groups_list, jacobian_groups_list, weights_groups_list)

                # 将 new_weights_list 按 indices_DBSCAN 对应分组写回权重
                try:
                    for i, k in enumerate(cluster_keys):
                        idxs = indices_DBSCAN[k]
                        w_k = np.asarray(new_weights_list[i]).reshape(-1)
                        if w_k.size == 1:
                            w_k = np.full(len(idxs), float(w_k[0]))
                        elif w_k.size != len(idxs):
                            min_len = min(len(idxs), w_k.size)
                            w_k = w_k[:min_len]
                            idxs = idxs[:min_len]
                        # 就地Update该组权重
                        weights_groups_list[i] = w_k
                        # optional: 同步Update扁平权重数组（若为一维数值数组）
                        try:
                            wb = np.asarray(weights_BDSCAN)
                            if wb.ndim == 1 and np.issubdtype(wb.dtype, np.number):
                                idx_arr = np.asarray(idxs, dtype=int)
                                valid_mask = (idx_arr >= 0) & (idx_arr < wb.size)
                                wb[idx_arr[valid_mask]] = w_k[valid_mask]
                                weights_BDSCAN = wb
                        except Exception as _:
                            pass
                except Exception as e:
                    print(f"更新分组权重失败: {e}")             

                # according to分组数输出权重因子
                self.pr_weight_factor = 1.0  # 默认第一分组为基准
                weight_factors = []
                if len(sigma2_values) > 1:
                    for i in range(len(sigma2_values)):
                        if i == 0:
                            weight_factors.append(1.0)
                        else:
                            # 以第一分组为基准, Calculate相对权重因子
                            weight_factors.append(sigma2_values[0] / sigma2_values[i])
                    self.dop_weight_factor = weight_factors[1]  # 兼容原有变量
                else:
                    weight_factors = [1.0]
                    self.dop_weight_factor = 1.0

                # 记录HVCE历史
                self.hvce_history.append(np.sqrt(sigma2_values))

                # 权重因子调试输出
                print(f"HVCE迭代 {iter_num+1}:")
                for i, sigma2 in enumerate(sigma2_values):
                    print(f"  分组{i+1} 方差因子 = {np.sqrt(sigma2):.6f}, 权重因子 = {weight_factors[i]:.3f}")
                # 兼容原有输出
                if len(sigma2_values) >= 2:
                    print(f"  相对权重(分组2/分组1) = {weight_factors[1]:.3f}")
                
                # if权重变化很小, 认为HVCE已收敛
                if hvce_converged:
                    print(f"HVCE在第{iter_num+1}次迭代收敛")
                    prev_state = estimated_state.copy()
                    current_state = estimated_state.copy()
                    break
            
            # Check收敛 - 恢复并改进收敛判断
            state_change = np.linalg.norm(estimated_state - prev_state) / np.linalg.norm(prev_state)
            
            # Check状态收敛或HVCE收敛
            converged = False
            if state_change < convergence_threshold:
                print(f"状态收敛于第{iter_num+1}次迭代, 状态变化: {state_change:.6f}")
                converged = True
            
            # ifusingHVCE并且HVCE已收敛, 也认为整体收敛
            if self.use_hvce and use_doppler and hvce_converged:
                print(f"HVCE收敛于第{iter_num+1}次迭代")
                converged = True
            
            if converged:
                # Update状态
                prev_state = estimated_state.copy()
                current_state = estimated_state.copy()
                break
                
            # Update状态
            prev_state = estimated_state.copy()
            current_state = estimated_state.copy()
            
            if iter_num == max_iterations - 1:
                print(f"定位未完全收敛, 最大迭代次数已达到, 最终状态变化: {state_change:.6f}")
                
            # 输出本次iteration权重信息
            if self.use_hvce and use_doppler:
                print(f"迭代{iter_num+1}: 权重因子更新完成, 继续下一次迭代")
        
        # # 将pseudorange与Doppler残差顺序Write同一file（先pseudorange, 后Doppler）
        # try:
        #     output_path = "residuals_output.txt"
        #     with open(output_path, "w", encoding="utf-8") as f:
        # # Writepseudorange残差
        #         for v in np.ravel(pr_residuals):
        #             f.write(f"{v}\n")
        # # 若存在Doppler残差, 则继续追加Write
        #         if use_doppler and doppler_obs is not None:
        #             for v in np.ravel(dop_residuals):
        #                 f.write(f"{v}\n")
        # print(f"残差已Write {output_path}")
        # except Exception as e:
        # print(f"Write残差filefailed: {e}")


        # 提取final结果
        # 收敛后再usingminimum二乘法求解
        if self.use_hvce and use_doppler:
            # result = least_squares(
            #     residual_func, 
            #     current_state, 
            #     jac=weighted_jacobian_func, 
            #     method='trf',
            #     loss='linear',
            #     verbose=0
            # )

            result = least_squares(
                residual_func, 
                current_state,  
                jac=weighted_jacobian_func, 
                method='lm',
                loss='linear',
                ftol=1e-6,      # 放宽收敛容差
                xtol=1e-6,      # 放宽收敛容差
                gtol=1e-6,      # 放宽收敛容差
                max_nfev=20,    # 限制函数评估次数
                verbose=0       # 设为1可以查看迭代过程, 帮助调试
            )
            # result = least_squares(
            #     residual_func, 
            #     current_state,  
            #     jac=weighted_jacobian_func, 
            #     method='trf',
            #     loss='linear',
            # ftol=1e-6,      # 放宽收敛容差
            # xtol=1e-6,      # 放宽收敛容差
            # gtol=1e-6,      # 放宽收敛容差
            # max_nfev=20,    # 限制函数评估次数
            # tr_solver='exact',  # using精确求解器
            # tr_options={'initial_trust_radius': 1.0},  # Set较小initial信任区域半径
            # verbose=0       # 设为1可以查看iteration过程, 帮助调试
            # )
            estimated_state = result.x

        estimated_position = estimated_state[:3]
        estimated_clock_bias = estimated_state[3]
        estimated_clock_drift = estimated_state[4]
        
        # Calculate残差
        if use_doppler and doppler_obs is not None:
            all_residuals = np.concatenate([pr_residuals, dop_residuals])
        else:
            all_residuals = pr_residuals
            
        # 区分不同类型残差
        n_pr_obs = len(pr_residuals)
        pr_residuals_final = all_residuals[:n_pr_obs]
        
        # Calculate协方差matrix
        J = weighted_jacobian_func(estimated_state)
        
        if self.use_hvce and use_doppler:
            try:
                cov = np.linalg.inv(J.T @ J)
            except np.linalg.LinAlgError:
                print("警告: 协方差计算中出现奇异矩阵, 使用伪逆")
                cov = np.linalg.pinv(J.T @ J)
        else:
            # using残差平方和除以自由度作为方差估计
            sigma2 = np.sum(all_residuals**2) / len(all_residuals)     
            try:
                cov = sigma2 * np.linalg.inv(J.T @ J)
            except np.linalg.LinAlgError:
                print("警告: 协方差计算中出现奇异矩阵, 使用伪逆")
                cov = sigma2 * np.linalg.pinv(J.T @ J)
        
        # Calculateposition精度指标
        pos_std = np.sqrt(np.diag(cov)[:3])
        clock_bias_std = np.sqrt(cov[3, 3]) 
        clock_drift_std = np.sqrt(cov[4,4])
        
        # ConvertECEFposition到大地坐标
        lat, lon, h = ecef_to_geodetic(*estimated_position)
        lat_deg, lon_deg = math.degrees(lat), math.degrees(lon)
        
        # 输出结果
        print("\n伪距定位结果:")
        print(f"位置 (ECEF): {estimated_position[0]:.2f}, {estimated_position[1]:.2f}, {estimated_position[2]:.2f} m")
        print(f"经纬度: {lat_deg:.6f}°, {lon_deg:.6f}°, 高度: {h:.2f} m")
        print(f"钟差: {estimated_clock_bias:.2f} m")
        print(f"钟漂:{estimated_clock_drift:.2f} m/s")
        print(f"伪距残差 RMS: {np.sqrt(np.mean(pr_residuals_final**2)):.3f} m")
        if use_doppler and doppler_obs is not None:
            dop_residuals_final = all_residuals[n_pr_obs:]
            print(f"多普勒残差 RMS: {np.sqrt(np.mean(dop_residuals_final**2)):.3f} Hz")
        print(f"位置Std Dev: {pos_std[0]:.2f}, {pos_std[1]:.2f}, {pos_std[2]:.2f} m")
        print(f"钟差Std Dev: {clock_bias_std:.2f} m")
        print(f"钟漂Std Dev: {clock_drift_std:.2f} m/s")
        
        # ifusing了HVCE, 输出final权重比例
        if self.use_hvce and use_doppler and len(self.hvce_history) > 0:
            final_pr_factor = self.hvce_history[-1][0]
            final_dop_factor = self.hvce_history[-1][1]
            print(f"\nHVCE最终权重因子:")
            print(f"伪距方差因子: {final_pr_factor:.6f}")
            print(f"多普勒方差因子: {final_dop_factor:.6f}")
            print(f"相对权重比(伪距:多普勒): 1:{final_pr_factor/final_dop_factor:.2f}")
        
        # Update状态
        self.state[:3] = estimated_position
        
        return estimated_position, estimated_clock_bias, estimated_clock_drift, cov, pr_residuals_final,self.hvce_history

    def grid_search_initialization(self, times, week, prnlist, ephemeris, doppler_obs, init_pos):
        """using网格搜索找到更好initialposition
        
        Args:
            times: 观测时间column表
            prnlist: satellitePRN号column表
            ephemeris: ephemerisdata
            doppler_obs: Dopplerobservation
            init_pos: initialposition猜测
            
        Returns:
            best_pos: 网格搜索找到最佳position
        """
        # 网格搜索Args
        grid_size = 20000  # 20km网格range
        grid_points = 7    # 每维度7点
        min_elevation = 10  # 最小仰角（度）
        
        best_cost = float('inf')
        best_pos = init_pos
        
        # Create网格点
        x_grid = np.linspace(init_pos[0]-grid_size, init_pos[0]+grid_size, grid_points)
        y_grid = np.linspace(init_pos[1]-grid_size, init_pos[1]+grid_size, grid_points)
        z_grid = np.linspace(init_pos[2]-grid_size, init_pos[2]+grid_size, grid_points)
        
        total_points = len(x_grid) * len(y_grid) * len(z_grid)
        print(f"网格搜索总点数: {total_points}")
        
        # usingtqdm显示进度
        progress = tqdm(total=total_points, desc="网格搜索进度")
        
        for x in x_grid:
            for y in y_grid:
                for z in z_grid:
                    test_pos = np.array([x, y, z])
                    
                    # Calculate该position下残差
                    residuals = []
                    valid_obs = 0
                    
                    for i in range(len(times)):
                        # Calculatesatelliteposition和velocity
                        sat_pos, sat_vel = cal_sat_posvel(times[i], week[i], prnlist[i], ephemeris, test_pos)
                        
                        # Calculateelevation angle
                        elevation = calculate_elevation_angle(sat_pos, test_pos)
                        elevation_deg = np.degrees(elevation)
                        
                        # 只usingelevation angle大于阈值observation
                        if elevation_deg >= min_elevation:
                            valid_obs += 1
                            # Calculateline of sightvector
                            los = test_pos - sat_pos
                            range_ = np.linalg.norm(los)
                            los_unit = los / range_
                            
                            # Calculate预测Doppler频移
                            rel_vel = sat_vel
                            radial_vel = np.dot(rel_vel, los_unit)
                            pred_doppler = -(-F_L1/SPEED_OF_LIGHT * radial_vel)
                            
                            # Calculate残差
                            residuals.append(doppler_obs[i] - pred_doppler)
                    
                    # 只有当validobservation数量足够时才Calculate代价
                    if valid_obs >= 4:  # 至少需要4有效观测值
                        # using加权残差平方和作为代价函数
                        weights = np.ones_like(residuals)
                        cost = np.sum(weights * np.array(residuals)**2)
                        
                        if cost < best_cost:
                            best_cost = cost
                            best_pos = test_pos
                            print(f"\n找到更好初始位置, 代价: {best_cost:.2f}")
                            print(f"位置: {best_pos}")
                    
                    progress.update(1)
        
        progress.close()
        
        # 输出final结果
        print(f"\n网格搜索完成:")
        print(f"最佳初始位置: {best_pos}")
        print(f"最小代价: {best_cost:.2f}")
        
        return best_pos

    def _compute_pseudorange_jacobian(self, state, times,week, prnlist, ephemeris):
        """Calculatepseudorange观测雅可比matrix"""
        position = state[:3]
        n_params = len(state)  # 5参数: x,y,z,钟差,钟漂
        n_obs = sum(1 for prn in prnlist if prn != 0)
        jacobian = np.zeros((n_obs, n_params))
        
        idx = 0
        for i in range(len(times)):
            if prnlist[i] == 0:
                continue
                
            sat_pos, _ = cal_sat_posvel(times[i], week[i], prnlist[i], ephemeris, position)
            
            # Calculateline of sightvector
            los = position - sat_pos
            range_ = np.linalg.norm(los)
            los_unit = los / range_
            
            # positionArgs偏derivative
            jacobian[idx, :3] = los_unit
            
            # clock biasArgs偏derivative
            jacobian[idx, 3] = 1.0
            
            # clock driftArgs偏derivative
            jacobian[idx, 4] = times[i] - self.start_time
            
            idx += 1
        
        return jacobian

    def _compute_doppler_jacobian(self, state, times,week, prnlist, ephemeris):
        """CalculateDoppler观测雅可比matrix"""
        position = state[:3]
        n_params = len(state)  # 5参数: x,y,z,钟差,钟漂
        n_obs = sum(1 for prn in prnlist if prn != 0)
        jacobian = np.zeros((n_obs, n_params))
        
        idx = 0
        for i in range(len(times)):
            if prnlist[i] == 0:
                continue
                
            sat_pos, sat_vel = cal_sat_posvel(times[i], week[i], prnlist[i], ephemeris, position)
            
            # Calculateline of sightvector
            los = position - sat_pos
            range_ = np.linalg.norm(los)
            los_unit = los / range_
            
            # Calculatepositionderivative
            dlos_dpos = (np.eye(3) - np.outer(los_unit, los_unit)) / range_
            dDop_dpos = F_L1/SPEED_OF_LIGHT * np.dot(dlos_dpos, sat_vel)
            
            # positionArgs偏derivative
            jacobian[idx, :3] = -dDop_dpos  # 注意符号
            
            # clock biasArgs偏derivative - Doppler对clock bias不敏感
            jacobian[idx, 3] = 0.0
            
            # clock driftArgs偏derivative
            jacobian[idx, 4] = -1.0
            
            idx += 1
        
        return jacobian

def read_tle_file(filename):
    """
    从file中ReadTLEdata, 支持多satellite
    
    Args:
        filename: TLEfilepath
        
    Returns:
        satellites: includingsatelliteTLE信息column表, each元素为(name, line1, line2)元组
    """
    satellites = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        
    i = 0
    while i < len(lines):
        # 跳过空line
        if not lines[i].strip():
            i += 1
            continue
        
        # Read三line(satellite名称和两lineTLE)
        if i + 2 < len(lines):
            name = lines[i].strip()
            line1 = lines[i+1].strip()
            line2 = lines[i+2].strip()
            
            # Checkwhether为validTLEline
            if line1.startswith('1 ') and line2.startswith('2 '):
                satellites.append((name, line1, line2))
            
            i += 3
        else:
            break
    
    return satellites

def add_tle_satellite(name, line1, line2, prn):
    """
    添加TLEsatellite
    
    Args:
        name: satellite名称
        line1: TLE第一line
        line2: TLE第二line
        prn: 为satellite分配PRN号(for仿真)
    """
    global tle_satellites
    try:
        # CreateSGP4satellite对象
        satellite = Satrec.twoline2rv(line1, line2)
        tle_satellites[prn] = {
            'name': name,
            'line1': line1,
            'line2': line2,
            'satellite': satellite,
            'prn': prn
        }
        return True
    except Exception as e:
        print(f"添加TLE卫星失败: {e}")
        return False

def load_tle_file(filename):
    """
    从fileLoadTLEdata
    
    Args:
        filename: TLEfilepath
    
    Returns:
        loaded_count: successfulLoadNumber of satellites
    """
    global tle_satellites

    satellites = read_tle_file(filename)
    loaded_count = 0
    
    for i, (name, line1, line2) in enumerate(satellites):
        # 为eachsatellite分配PRN号, 从1开始
        prn = i + 1
        if add_tle_satellite(name, line1, line2, prn):
            loaded_count += 1
    
    print(f"成功加载了{loaded_count}卫星TLE数据")
    return loaded_count

def gps_time_to_datetime(gps_week, gps_seconds):
    """
    将GPS time(周和秒)Convert为datetime对象
    
    Args:
        gps_week: GPS week
        gps_seconds: GPS周内秒
        
    Returns:
        datetime: 对应datetime对象
    """
    # GPS起始时间: 1980-01-06 00:00:00
    gps_epoch = datetime(1980, 1, 6)
    
    # Calculate从GPS纪元开始total秒数
    total_seconds = gps_week * 604800 + gps_seconds
    
    # Calculatecurrent时间
    return gps_epoch + timedelta(seconds=total_seconds)

def datetime_to_jday(dt):
    """
    将datetime对象Convert为Julian day
    
    Args:
        dt: datetime对象
        
    Returns:
        (jd, fr): Julian day和日小数部分
    """
    year = dt.year
    month = dt.month
    day = dt.day
    hour = dt.hour
    minute = dt.minute
    second = dt.second + dt.microsecond / 1e6
    
    # CalculateJulian day(整数部分)
    if month <= 2:
        year -= 1
        month += 12
    
    A = year // 100
    B = 2 - A + A // 4
    
    jd = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day + B - 1524.5
    
    # Calculate日小数部分
    fr = (hour + minute / 60.0 + second / 3600.0) / 24.0
    
    return jd, fr

def calculate_satellite_position_tle(prn,gps_week, gps_time):
    """
    usingTLEdataCalculatesatelliteposition和velocity
    
    Args:
        prn: satellitePRN号
        gps_time: GPS time(秒)
        
    Returns:
        sat_pos: satelliteECEFposition[x,y,z](m)
        sat_vel: satelliteECEFvelocity[vx,vy,vz](m/s)
    """
    if prn not in tle_satellites:
        raise ValueError(f"PRN={prn}卫星不存在于TLE数据中")
    
    # Getsatellite对象
    satellite = tle_satellites[prn]['satellite']
    
    # ConvertGPS time为datetime
    dt = gps_time_to_datetime(gps_week, gps_time)
    
    # Convert为Julian day
    jd, fr = datetime_to_jday(dt)
    
    # usingSGP4Calculatesatelliteposition和velocity
    error_code, position, velocity = satellite.sgp4(jd, fr)
    
    if error_code != 0:
        print(f"SGP4计算错误, 错误代码: {error_code}, 使用零向量替代")
        return np.zeros(3), np.zeros(3)
    
    ecef_position, ecef_velocity = teme_to_ecef(position, velocity, dt)

    # Convertunit: km到m
    sat_pos = np.array(ecef_position) * 1000.0
    sat_vel = np.array(ecef_velocity) * 1000.0
    
    return sat_pos, sat_vel

def load_sim_data(filename, nav_file=None, ephemeris=None, tle_file=None):
    """
    Load仿真datafile并Calculatesatelliteposition和velocity
    
    Args:
        filename: 仿真datafile名
        nav_file: 导航ephemerisfile名(nav.txt)
        ephemeris: satelliteephemeris字典(optional, if不提供则从nav_fileRead)
        
    Returns:
        data: data字典including观测data和satellite信息
    """
    # Readfile, 跳过前两line注释
    try:
        # df = pd.read_csv(filename, delimiter=r'\s+', skiprows=2, 
        #                 names=['gps_time', 'pseudorange', 'carrier_phase', 'doppler', 'elevation', 'prn'])
        df = pd.read_csv(filename, delimiter=r'\s+', skiprows=0,
                        usecols=[2,3,4,5,6,7,9],  
                        names=['gps_time','week','rece_time', 'pseudorange', 'doppler', 'snr', 'prn'])
    except Exception as e:
        print(f"无法读取文件: {filename}")
        print(f"错误: {e}")
        return None

    # if提供了导航file, ReadephemerisArgs
    if nav_file and ephemeris is None:
        try:
            ephemeris = {}  # 按PRN和toc存储星历字典
            current_prn = None
            current_toc = None
            current_eph = {}
            
            with open(nav_file, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line:  # 跳过空行
                        continue
                        
                    # 解析Argsline
                    if 'eph.' in line:
                        parts = line.split('=')
                        if len(parts) != 2:
                            continue
                            
                        param_name = parts[0].strip()  # 例如: "eph.prn" 或 "eph.toc"
                        
                        # correction: 去除值末尾可能存在分号(;)和其他非数字字符
                        param_value_str = parts[1].strip()
                        if param_value_str.endswith(';'):
                            param_value_str = param_value_str[:-1]  # 去除末尾分号
                        
                        # Convert为数值
                        try:
                            param_value = float(param_value_str)
                        except ValueError:
                            print(f"警告: 无法将 '{param_value_str}' 转换为数值, 在参数 {param_name}")
                            continue
                        
                        # 提取Args名（去掉"eph."前缀）
                        param_key = param_name.split('.')[1]
                        
                        # 处理PRN和toc特殊情况
                        if param_key == 'prn':
                            current_prn = int(param_value)
                            if current_prn not in ephemeris:
                                ephemeris[current_prn] = {}
                        elif param_key == 'toe':
                            current_toc = param_value
                            if current_prn is not None:
                                current_eph = {}
                                ephemeris[current_prn][current_toc] = current_eph
                        else:
                            # 存储其他Args
                            if current_prn is not None and current_toc is not None:
                                current_eph[param_key] = param_value
            
            print(f"已加载星历数据:")
            for prn in ephemeris:
                toc_count = len(ephemeris[prn])
                print(f"PRN {prn}: {toc_count}时刻星历")
                
        except Exception as e:
            print(f"读取导航文件时出错: {e}")
            import traceback
            traceback.print_exc()
            return None
    elif tle_file and os.path.exists(tle_file):
        # LoadTLEfile
        print(f"加载TLE文件: {tle_file}")
        count = load_tle_file(tle_file)
        if count > 0:
            # usingTLEdata
            print(f"使用TLE模式, 已加载{count}卫星")
    else:
    # usingdefaultephemeris
        ephemeris = {}  # 按PRN和toc存储星历字典
        
 
    times = df['gps_time'].values   
    data = {
        'times': times,
        'week': df['week'].values,
        'doppler': df['doppler'].values,
        'rece_time': df['rece_time'].values,
        'pseudorange': df['pseudorange'].values,
        'snr': df['snr'].values,
        'prn': df['prn'].values,
    }
    
    print("已加载观测数据:")
    print(f"观测时间range: {times[0]:.1f}s - {times[-1]:.1f}s")
    print(f"观测点数: {len(times)}")
    print(f"平均伪距: {np.mean(df['pseudorange'].values):.2f}m")
    
    return data,ephemeris

def cal_sat_posvel(time, week, prn, ephemeris,rec_pos):
    """
    Calculatesatelliteposition和velocity
    
    Args:
        time: 观测时间column表
        week: 观测周数column表
        prn: satellitePRN号column表
        ephemeris: ephemeris字典
        rec_pos: receiverposition(ECEF坐标)
        
    Returns:
        sat_positions: satellitepositioncolumn表
        sat_velocities: satellitevelocitycolumn表
    """
    # Calculateeach观测时刻satelliteposition和velocity（iteration法considering传播时延）
    sat_positions = []
    sat_velocities = []
    c = SPEED_OF_LIGHT  # 光速(m/s)

    if Tle_sat_pos:
        return calculate_satellite_position_tle(prn, week,time)
    else:
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
        
        # initial信号传播时延10ms
        time_delay = 0.01
        
        # initial发射时间估计
        emission_time = time - time_delay

        emission_time = time
        # iteration求解精确发射时间
        max_iter = 10  # 最大迭代次数
        for iter_count in range(max_iter):
            # 在current估计发射时间Calculatesatelliteposition
            sat_pos, sat_vel = calculate_satellite_position(current_eph, emission_time)
            
            # Calculatesatellite到参考positiongeometric distance
            geometric_range = np.linalg.norm(sat_pos - rec_pos)
            
            # Update信号传播时延
            new_time_delay = geometric_range / c
            
            # Update发射时间
            new_emission_time = time - new_time_delay
            new_emission_time = time
            # Check收敛条件
            if abs(new_emission_time - emission_time) < 1e-6:
                break
                
            emission_time = new_emission_time
        
        # usingfinaliteration发射时间Calculatesatelliteposition和velocity
        sat_positions, sat_velocities = calculate_satellite_position(current_eph, emission_time)

    return sat_positions,sat_velocities

def calculate_elevation_angle(sat_pos, rec_pos):
    """
    Calculatesatellite相对于receiverelevation和azimuth
    Args:
        sat_pos: satelliteposition (ECEF坐标)
        rec_pos: receiverposition (ECEF坐标)
    Returns:
        elevation_angle: elevation (degrees)
        azimuth_angle: azimuth (degrees)
    """
    los_vector = sat_pos - rec_pos
    los_norm = np.linalg.norm(los_vector)
    
    # Convertreceiverposition到地理坐标
    rec_lat, rec_lon, rec_alt = ecef_to_geodetic(*rec_pos)
    
    # Calculate本地坐标系基vector
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
    
    # 投影到本地水平面
    los_local = np.array([
        np.dot(los_vector, east_vector),
        np.dot(los_vector, north_vector),
        np.dot(los_vector, up_vector)
    ])
    
    # Calculateelevation
    elevation_angle = np.arcsin(los_local[2]/los_norm) * (180.0 / np.pi)
    
    # Calculateazimuth（0-360度）
    azimuth_angle = np.arctan2(los_local[0], los_local[1]) * (180.0 / np.pi)
    azimuth_angle = azimuth_angle % 360  # 规范化到0-360度range
    
    return elevation_angle, azimuth_angle


def run_lsq_positioning(sim_data_file, ephemeris_file=None, tle_file=None, init_pos=None, 
                         height_constraint=None, true_position=None, height_constraint_mode='hard',
                         adaptive_weights_enable=True, threshold=5.0, min_weight=0.1, 
                         enable_iono=True, enable_tropo=True, enable_multipath=False, enable_relativistic=False, enable_hardware=False, enable_pco=False,fit_doppler=False, 
                         use_pseudorange=False, use_doppler=True,use_hvce=False,use_DBSCAN=False,
                         progress_callback=None, selected_prns=None):
    """
    运lineminimum二乘批处理positioning, 支持Doppler、pseudorange或两者结合
    
    Args:
        sim_data_file: 仿真datafile
        ephemeris_file: ephemerisfile(optional)
        tle_file: TLEfile(optional)
        init_pos: initialposition猜测(optional)
        height_constraint: altitude约束(optional)
        true_position: 真实positionfor评估(optional)
        height_constraint_mode: altitude约束模式, 'hard', 'soft', 或 'both'
        adaptive_weights_enable: whetherenable自适应权重
        threshold: 残差阈值
        min_weight: minimum权重限制
        atmcor: whether应用大气correction
        fit_doppler: whether对Doppler进line傅里叶拟合
        use_pseudorange: whetherusingpseudorange进linepositioning
        use_doppler: whetherusingDoppler进linepositioning
        progress_callback: 进度回调函数 (current, total, message) -> bool
        selected_prns: 选择satellitePRNcolumn表(None表示全部)
        
    Returns:
        positioning: positioning处理对象
        results: positioning结果
    """


    # Load仿真data
    data, ephemeris = load_sim_data(sim_data_file, ephemeris_file, None, tle_file)
    if data is None:
        return None, None
    
    # Read观测data
    times = data['times']
    week = data['week']
    rece_time = data['rece_time']
    prn_list = data['prn']
    pseudorange_obs = data['pseudorange'] if use_pseudorange else None
    doppler_obs = data['doppler'] if use_doppler else None
    snr_values = data.get('snr', None)
    
    # satellite筛选: if指定了selected_prns, 只保留选中satellite观测
    if selected_prns is not None and len(selected_prns) > 0:
        selected_set = set(str(p) for p in selected_prns)
        mask = np.array([str(prn) in selected_set for prn in prn_list])
        
        times = np.array(times)[mask].tolist()
        week = np.array(week)[mask].tolist()
        rece_time = np.array(rece_time)[mask].tolist()
        prn_list = np.array(prn_list)[mask].tolist()
        if pseudorange_obs is not None:
            pseudorange_obs = np.array(pseudorange_obs)[mask].tolist()
        if doppler_obs is not None:
            doppler_obs = np.array(doppler_obs)[mask].tolist()
        if snr_values is not None:
            snr_values = np.array(snr_values)[mask].tolist()
        
        # Updatedata字典以便后续using
        data['times'] = times
        data['week'] = week
        data['rece_time'] = rece_time
        data['prn'] = prn_list
        if use_pseudorange:
            data['pseudorange'] = pseudorange_obs
        if use_doppler:
            data['doppler'] = doppler_obs
        if snr_values is not None:
            data['snr'] = snr_values
        
        print(f"卫星筛选: 使用 {len(selected_prns)} 卫星 (PRN: {selected_prns})")

    positioning = SingleSatDopplerPositioning(
        init_pos, 
        height_constraint=height_constraint,
        height_constraint_mode=height_constraint_mode
    )
    positioning.rece_time = rece_time
    positioning.use_hvce=use_hvce
    if use_hvce:
        use_DBSCAN = True
    positioning.use_DBSCAN = use_DBSCAN

    positioning.enable_iono = enable_iono
    positioning.enable_tropo = enable_tropo
    positioning.enable_multipath = enable_multipath
    positioning.enable_relativistic = enable_relativistic
    positioning.enable_hardware = enable_hardware
    positioning.enable_pco = enable_pco

    if tle_file and os.path.exists(tle_file):
        global Tle_sat_pos
        Tle_sat_pos = True
        ephemeris = tle_satellites

    positioning.configure_adaptive_weights(adaptive_weights_enable, threshold, min_weight)
    
    # 存储收敛过程data
    convergence_data = {
        'obs_count': [],
        'position_error': [],
        'position_std': [],  
        'height_error': [],
        'latitude': [],
        'longitude': [],
        'height': [],
        'clock_bias': [],
        'clock_drift': [],
        'residual_rms': [],
        'DGDOP': [],
        'cov': [],
        'hvce_history':[],
        'position_history': []  # 添加位置历史记录（ECEF坐标）
    }
    
    # datatotal长度
    total_obs = len(data['times'])
    
    # 定义观测点数增长序column
    # 开始using较少点, 然后逐渐增加
    obs_steps = []
    
    # 前100点以10点为步长
    for i in range(50, min(200, total_obs), 20):
        obs_steps.append(i)
    
    # 前100点以10点为步长
    for i in range(200, min(300, total_obs), 20):
        obs_steps.append(i)

    # 后面以更大步长增加
    step_size = 20
    for i in range(300, total_obs, step_size):
        obs_steps.append(i)
    
    # 确保including全部observation
    if total_obs not in obs_steps:
        obs_steps.append(total_obs)
    
    print(f"总观测数: {total_obs}")
    print(f"收敛分析步骤: {obs_steps}")

    init_clock=None
    # Keep the latest subset so GUI can plot details (skyplot/residuals)
    last_subset_times = None
    last_subset_prn = None
    last_subset_snr = None
    last_subset_pr = None
    last_subset_doppler = None
    last_clock_bias = 0.0
    last_clock_drift = 0.0

    # 逐步增加观测点数量进linepositioning
    for step_idx, obs_count in enumerate(obs_steps):
        print(f"\n处理 {obs_count}/{total_obs} 观测值...")
        
        # 调用进度回调
        if progress_callback is not None:
            should_continue = progress_callback(
                step_idx + 1, 
                len(obs_steps), 
                f"正在处理 {obs_count}/{total_obs} 观测值..."
            )
            if should_continue is False:
                print("用户取消操作")
                break
        
        # 截取前obs_countobservation
        subset_times = np.asarray(data['times'][:obs_count])
        subset_week = np.asarray(data['week'][:obs_count])
        subset_prn = np.asarray(data['prn'][:obs_count])
        subset_snr = np.asarray(data['snr'][:obs_count]) if 'snr' in data else None
        
        # according tousing观测类型准备data
        subset_pr = np.asarray(data['pseudorange'][:obs_count]) if use_pseudorange else None
        subset_doppler = np.asarray(data['doppler'][:obs_count]) if use_doppler else None

        # update "last" references for plotting later
        last_subset_times = subset_times
        last_subset_week = subset_week
        last_subset_prn = subset_prn
        last_subset_snr = subset_snr
        last_subset_pr = subset_pr
        last_subset_doppler = subset_doppler
        
        positioning.DBSCAN_result = False
        try:
            # according tousing观测类型执line不同positioning算法
            if use_pseudorange:
                # usingpseudorange进linepositioning
                position, clock_bias, clock_drift, cov, residuals,hvce_history = positioning.batch_least_squares_with_pseudorange(
                    subset_times, subset_week, subset_prn, ephemeris, subset_pr, 
                    subset_doppler if use_doppler else None, 
                    subset_snr, init_pos, use_doppler
                )
            else:
                # 仅usingDoppler进linepositioning
                result = positioning.batch_least_squares_position_only(
                    subset_times, subset_week, subset_prn, ephemeris, subset_doppler, subset_snr, init_pos=init_pos
                )
                
                # CheckReturns值数量
                if height_constraint_mode == 'both' and len(result) == 6:
                    position, clock_drift, cov, residuals, once_residuals, compare_results = result
                    clock_bias = 0.0  # 多普勒模式下无钟差估计
                else:
                    position, clock_drift, cov, residuals, once_residuals = result
                    clock_bias = 0.0  # 多普勒模式下无钟差估计
                
                hvce_history=None
            
            # Remember latest estimated clocks for plotting later
            last_clock_bias = float(clock_bias) if clock_bias is not None else 0.0
            last_clock_drift = float(clock_drift) if clock_drift is not None else 0.0
            
            # Convert结果
            lat, lon, h = ecef_to_geodetic(*position)
            
            # 确保协方差matrix存在
            if cov is not None:
                position_std = np.sqrt(np.diag(cov)[:3])
                clock_bias_std = np.sqrt(cov[3, 3]) if use_pseudorange else 0.0
                clock_drift_std = np.sqrt(cov[4 if use_pseudorange else 3, 4 if use_pseudorange else 3])
                DGDOP = np.sqrt(np.sum(np.diag(cov)))
            else:
                position_std = np.array([float('nan'), float('nan'), float('nan')])
                clock_bias_std = float('nan')
                clock_drift_std = float('nan')
            
            # Calculatepositioning误差(if有真实position)
            position_error = 0
            if true_position is not None:
                position_error = np.linalg.norm(position - true_position)
            
            # Calculatealtitude误差(if有altitude约束)
            height_error = 0
            if height_constraint is not None:
                height_error = abs(h - height_constraint)
            
            # 记录收敛data
            convergence_data['obs_count'].append(obs_count)
            convergence_data['position_error'].append(position_error)
            convergence_data['position_std'].append(np.linalg.norm(position_std))  # 添加位置Std Dev模值
            convergence_data['height_error'].append(height_error)
            convergence_data['latitude'].append(math.degrees(lat))
            convergence_data['longitude'].append(math.degrees(lon))
            convergence_data['height'].append(h)
            convergence_data['clock_bias'].append(clock_bias)
            convergence_data['clock_drift'].append(clock_drift)
            convergence_data['residual_rms'].append(np.sqrt(np.mean(np.square(residuals))))
            convergence_data['DGDOP'].append(DGDOP)
            convergence_data['cov'].append(cov)
            convergence_data['hvce_history'].append(hvce_history)
            convergence_data['position_history'].append(list(position))  # 记录ECEF位置
            
            # 简要输出current结果
            print(f"观测点数: {obs_count}")
            print(f"位置: {math.degrees(lat):.6f}°, {math.degrees(lon):.6f}°, {h:.2f}m")
            if true_position is not None:
                print(f"位置误差: {position_error:.2f}m")
            print(f"钟差: {clock_bias:.2f}m, 钟漂: {clock_drift:.3f}Hz")
            print(f"残差RMS: {np.sqrt(np.mean(np.square(residuals))):.3f} {'m' if use_pseudorange else 'Hz'}")
            
        except Exception as e:
            print(f"处理 {obs_count} 观测值时失败: {e}")
            import traceback
            traceback.print_exc()
    
    # 绘制收敛过程曲线
    # if len(convergence_data['obs_count']) > 1:
        # # Set中文字体
        # plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'STSong', 'Arial Unicode MS']
        # plt.rcParams['axes.unicode_minus'] = False
        
        # # Create2x3子图布局
        # fig, axs = plt.subplots(2, 3, figsize=(18, 10))
        
        # # 1. position误差随观测点数变化
        # if true_position is not None:
        #     axs[0, 0].plot(convergence_data['obs_count'], convergence_data['position_error'], 'b-o')
        # axs[0, 0].set_title('position误差随观测点数变化')
        # axs[0, 0].set_xlabel('观测点数')
        # axs[0, 0].set_ylabel('position误差 (m)')
        # # Sety轴上限为1000米
        #     axs[0, 0].set_ylim(0, min(1000, max(convergence_data['position_error'])*1.1))
        #     axs[0, 0].grid(True)
              
        # plt.tight_layout()
        # plt.savefig('positioning_convergence.png', dpi=300)
        # plt.show()

        # 绘制单收敛曲线
        # plt.figure(figsize=(12, 8))
        
        # #position误差随观测点数变化                   
        # plt.plot(convergence_data['obs_count'], convergence_data['position_error'], 'bo', label='Doppler', linewidth=2, markersize=8)
        # # Set图表属性
        # plt.ylim(0, min(1000, max(convergence_data['position_error'])*1.1))
        # plt.xlabel('Time (s)', fontsize=20)
        # plt.ylabel('Position Error (m)', fontsize=20)
        # plt.grid(True, linestyle='--', alpha=0.7)
        # plt.tick_params(axis='both', which='major', labelsize=20)
        # plt.tight_layout()
        # plt.savefig('positioning_convergence.png', dpi=300)
        # plt.show()
    # using最后一观测结果绘制残差散点图
    # if obs_count == obs_steps[-1]:
    # # usingcolumn表推导式过滤valid观测data
    #     valid_indices = [i for i, prn in enumerate(subset_prn) if prn != 0]
    #     valid_times = [subset_times[i] for i in valid_indices]
        
    # # 绘制残差散点图
    #     if use_pseudorange:
    #         plot_residual_scatter(residuals, valid_times, f"{'pseudorange' if use_pseudorange else 'doppler'}_residuals.png")
    #     else:
    #         plot_residual_scatter(residuals, valid_times)
    #         plot_residual_scatter_first(once_residuals, valid_times)

    # 构建final结果字典(using全部observation结果)
    if len(convergence_data['obs_count']) > 0:
        idx = -1  # 最后一组结果
        results = {
            'final_position': positioning.state[:3],
            'latitude_deg': convergence_data['latitude'][idx],
            'longitude_deg': convergence_data['longitude'][idx],
            'height_m': convergence_data['height'][idx],
            'clock_bias': convergence_data['clock_bias'][idx],
            'clock_drift': convergence_data['clock_drift'][idx],
            'residual_rms': convergence_data['residual_rms'][idx],
            'observations_count': convergence_data['obs_count'][idx],
            'convergence_data': convergence_data,
            'used_pseudorange': bool(use_pseudorange),
            'used_doppler': bool(use_doppler),
        }

        # ---- Extra outputs for GUI plots (skyplot + per-epoch residuals) ----
        try:
            if last_subset_times is not None and last_subset_prn is not None:
                times_arr = np.asarray(last_subset_times, dtype=float)
                week_arr = np.asarray(last_subset_week, dtype=int)
                prn_arr = np.asarray(last_subset_prn)
                valid_mask = prn_arr != 0

                times_valid = times_arr[valid_mask]
                week_valid = week_arr[valid_mask]
                prn_valid = prn_arr[valid_mask].astype(int)

                # Sky plot (az/el) based on final estimated receiver position
                rec_pos = np.asarray(positioning.state[:3], dtype=float)
                sky_el = []
                sky_az = []
                for t, w, prn in zip(times_valid,week_valid, prn_valid):
                    sat_pos, _sat_vel = cal_sat_posvel(float(t), w, int(prn), ephemeris, rec_pos)
                    el_deg, az_deg = calculate_elevation_angle(sat_pos, rec_pos)
                    sky_el.append(float(el_deg))
                    sky_az.append(float(az_deg))

                # Residuals at final solution
                pr_res = None
                dop_res = None

                if use_pseudorange and last_subset_pr is not None:
                    state5 = np.array([rec_pos[0], rec_pos[1], rec_pos[2], float(last_clock_bias), float(last_clock_drift)])
                    pr_res, _ = positioning.compute_residuals_pseudorange(
                        state5, times_arr, week_arr, prn_arr, ephemeris, last_subset_pr, last_subset_snr
                    )

                if use_doppler and last_subset_doppler is not None:
                    if use_pseudorange:
                        state5 = np.array([rec_pos[0], rec_pos[1], rec_pos[2], float(last_clock_bias), float(last_clock_drift)])
                        dop_res, _ = positioning.compute_residuals(
                            state5, times_arr, week_arr, prn_arr, ephemeris, last_subset_doppler, last_subset_snr
                        )
                    else:
                        state4 = np.array([rec_pos[0], rec_pos[1], rec_pos[2], float(last_clock_drift)])
                        dop_res, _ = positioning.compute_residuals_1(
                            state4, times_arr, week_arr, prn_arr, ephemeris, last_subset_doppler, last_subset_snr
                        )

                # Store arrays for GUI consumption (keep as plain lists to cross thread safely)
                results['skyplot'] = {
                    'times_s': (times_valid - times_valid[0]).tolist() if len(times_valid) else [],
                    'prn': prn_valid.tolist(),
                    'elevation_deg': sky_el,
                    'azimuth_deg': sky_az,
                }
                results['residuals'] = {
                    'times_s': (times_valid - times_valid[0]).tolist() if len(times_valid) else [],
                    'pseudorange_m': pr_res.tolist() if pr_res is not None else None,
                    'doppler_hz': dop_res.tolist() if dop_res is not None else None,
                }
        except Exception:
            # Plot extras are best-effort; do not break positioning result
            pass
        
        # 显示final结果
        print("\n最终定位结果(使用全部观测值):")
        print(f"纬度: {results['latitude_deg']:.6f}°")
        print(f"经度: {results['longitude_deg']:.6f}°")
        print(f"高程: {results['height_m']:.2f} m")
        if use_pseudorange:
            print(f"钟差: {results['clock_bias']:.3f} m")
        print(f"钟漂: {results['clock_drift']:.3f} m/s")
        print(f"处理观测数: {results['observations_count']}")
        
        if true_position is not None:
            position_error = np.linalg.norm(positioning.state[:3] - true_position)
            print(f"位置误差: {position_error:.2f} m")
        
        # # Calculateeach历元elevation和azimuth,绘制天空图
        # elevation_angles = []
        # azimuth_angles = []
        # for i in range(total_obs):
        #     sat_pos, _ = cal_sat_posvel(data['times'][i], data['prn'][i], ephemeris, init_pos)
        #     elevation_angle, azimuth_angle = calculate_elevation_angle(sat_pos, init_pos)
        #     elevation_angles.append(elevation_angle)
        #     azimuth_angles.append(azimuth_angle)
        
        # # 绘制天空图
        # plt.figure(figsize=(12, 10))
        
        # # Create极坐标子图
        # ax = plt.subplot(111, projection='polar')
        
        # # Set正北方在最上方
        # ax.set_theta_zero_location('N')
        # ax.set_theta_direction(-1)  # 顺时针方向
        
        # # Convert角度为弧度
        # azimuth_rad = np.array(azimuth_angles) * np.pi / 180.0
        # elevation_deg = np.array(elevation_angles)
        
        # # 绘制satellite轨迹（elevation作为径向distance, azimuth作为角度）
        # # 天空图中, 中心是头顶（90度elevation）, 边缘是地平线（0度elevation）
        # # 所以径向distance = 90 - elevation
        # radial_distance = 90 - elevation_deg
        
        # # 美化天空图

        # # 绘制轨迹线, using更粗线条和渐变色
        # from matplotlib.collections import LineCollection

        # # Generate轨迹分段坐标
        # points = np.array([azimuth_rad, radial_distance]).T.reshape(-1, 1, 2)
        # segments = np.concatenate([points[:-1], points[1:]], axis=1)
        # # 颜色映射: elevation由低到高渐变
        # norm = plt.Normalize(elevation_deg.min(), elevation_deg.max())
        # lc = LineCollection(segments, cmap='viridis', norm=norm, linewidth=3, alpha=0.85)
        # lc.set_array(elevation_deg)
        # ax.add_collection(lc)

        # # 绘制起始点
        # ax.plot(azimuth_rad[0], radial_distance[0], marker='o', markersize=12, 
        #         markerfacecolor='limegreen', markeredgecolor='k', label='Start', zorder=5)
        # # 绘制结束点
        # ax.plot(azimuth_rad[-1], radial_distance[-1], marker='*', markersize=18, 
        #         markerfacecolor='crimson', markeredgecolor='k', label='End', zorder=5)

        # # Set径向网格（elevation）
        # ax.set_ylim(0, 90)
        # ax.set_yticks([0, 15, 30, 45, 60, 75, 90])
        # # 天空图中心为头顶（90°）, 边缘为地平线（0°）
        # ax.set_yticklabels(['90°', '75°', '60°', '45°', '30°', '15°', '0°'], fontsize=13, color='navy')

        # # Set角度网格（azimuth）
        # ax.set_xticks(np.linspace(0, 2*np.pi, 8, endpoint=False))
        # ax.set_xticklabels(['0°', '45°', '90°', '135°', '180°', '225°', '270°', '315°'], fontsize=13, color='navy')

        # # Set标题和标签
        # ax.set_title('Satellite Sky Plot\n(Elevation-Azimuth)', fontsize=20, pad=25, color='darkblue', weight='bold')

        # # 添加图例
        # ax.legend(loc='upper right', bbox_to_anchor=(1.25, 1.05), fontsize=13, frameon=True, facecolor='white', edgecolor='gray')

        # # 添加网格, 增强可读性
        # ax.grid(True, alpha=0.5, linestyle='--', linewidth=1)

        # # 添加颜色条, 表示elevation
        # cbar = plt.colorbar(lc, ax=ax, pad=0.12, orientation='vertical', shrink=0.8)
        # cbar.set_label('Elevation Angle (°)', fontsize=14, color='navy')
        # cbar.ax.tick_params(labelsize=12)
        
        # # 调整布局
        # plt.tight_layout()
        # plt.savefig('satellite_sky_plot.png', dpi=300, bbox_inches='tight')
        # plt.show()
        # # 输出统计信息
        # print(f"\n天空图统计信息:")
        # print(f"Elevation range: {min(elevation_angles):.2f}° - {max(elevation_angles):.2f}°")
        # print(f"azimuthrange: {min(azimuth_angles):.2f}° - {max(azimuth_angles):.2f}°")
        # print(f"averageelevation: {np.mean(elevation_angles):.2f}°")
        # print(f"averageazimuth: {np.mean(azimuth_angles):.2f}°")

        # # Save为CSVfile
        # with open('skyplot_data.csv', 'w', newline='') as f:
        #     writer = csv.writer(f)
        # writer.writerow(['Time', 'Elevation_Angle', 'Azimuth_angles'])  # Writecolumn标题
        #     for t, elev,azl in zip(data['times'], elevation_angles,azimuth_angles):
        #         writer.writerow([t, elev, azl])

        # # 绘制elevation随历元变化图
        # plt.figure(figsize=(8, 6))
        # plt.scatter(data['times']-data['times'][0], elevation_angles, color='blue', marker='.', s=20, edgecolors='b', alpha=0.6)
        # # plt.title('elevation随历元变化')
        # plt.xlabel('Time (s)',fontsize=14)
        # plt.ylabel('Elevation Angle (°)',fontsize=14)
        # plt.grid(True, linestyle='--', alpha=0.5)
        # plt.tick_params(axis='both', which='major', labelsize=14)
        # plt.savefig('elevation_angle_vs_time.png', dpi=300)
        # plt.show()
        
    return convergence_data, results

def plot_residual_scatter(doppler_residuals, times, output_file='residual_scatter.png'):
    """
    绘制Doppler残差散点图
    
    Args:
        doppler_residuals: Doppler残差数组
        times: 观测时间数组
        output_file: 输出file名
    """
    # Set中文字体
    plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'STSong', 'Arial Unicode MS']
    plt.rcParams['axes.unicode_minus'] = False
    
    # Create图形
    plt.figure(figsize=(12, 6))
    
    # 绘制残差散点图
    plt.scatter(times, doppler_residuals, c=doppler_residuals, cmap='coolwarm', 
                alpha=0.7, s=50, edgecolor='k')
    
    # 添加颜色条
    cbar = plt.colorbar()
    cbar.set_label('残差大小 (Hz)')
    
    # 添加零线
    plt.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    
    # Set标题和标签
    plt.title('多普勒残差散点图')
    plt.xlabel('时间 (s)')
    plt.ylabel('残差 (Hz)')
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # 添加统计信息
    mean_res = np.mean(doppler_residuals)
    std_res = np.std(doppler_residuals)
    plt.text(0.02, 0.98, f'平均残差: {mean_res:.3f} Hz\nStd Dev: {std_res:.3f} Hz',
             transform=plt.gca().transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Save图像
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.show()
    
    # 输出统计信息
    print(f"\n残差统计信息:")
    print(f"观测点数: {len(doppler_residuals)}")
    print(f"平均残差: {mean_res:.3f} Hz")
    print(f"残差Std Dev: {std_res:.3f} Hz")
    print(f"最大残差: {np.max(np.abs(doppler_residuals)):.3f} Hz")
    print(f"残差RMS: {np.sqrt(np.mean(np.square(doppler_residuals))):.3f} Hz")

def plot_residual_scatter_first(doppler_residuals, times, output_file='first_residual_scatter.png'):
    """
    绘制Doppler残差散点图
    
    Args:
        doppler_residuals: Doppler残差数组
        times: 观测时间数组
        output_file: 输出file名
    """
    # Set中文字体
    plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'STSong', 'Arial Unicode MS']
    plt.rcParams['axes.unicode_minus'] = False
    
    # Create图形
    plt.figure(figsize=(12, 6))
    
    # 绘制残差散点图
    plt.scatter(times, doppler_residuals, c=doppler_residuals, cmap='coolwarm', 
                alpha=0.7, s=50, edgecolor='k')
    
    # 添加颜色条
    cbar = plt.colorbar()
    cbar.set_label('残差大小 (Hz)')
    
    # 添加零线
    plt.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    
    # Set标题和标签
    plt.title('首次迭代多普勒残差散点图(O-C)')
    plt.xlabel('时间 (s)')
    plt.ylabel('残差 (Hz)')
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # 添加统计信息
    mean_res = np.mean(doppler_residuals)
    std_res = np.std(doppler_residuals)
    plt.text(0.02, 0.98, f'平均残差: {mean_res:.3f} Hz\nStd Dev: {std_res:.3f} Hz',
             transform=plt.gca().transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Save图像
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.show()
    
    # 输出统计信息
    print(f"\n残差统计信息:")
    print(f"观测点数: {len(doppler_residuals)}")
    print(f"平均残差: {mean_res:.3f} Hz")
    print(f"残差Std Dev: {std_res:.3f} Hz")
    print(f"最大残差: {np.max(np.abs(doppler_residuals)):.3f} Hz")
    print(f"残差RMS: {np.sqrt(np.mean(np.square(doppler_residuals))):.3f} Hz")

class HelmertVCE:
    """赫尔默特方差分量估计方法实现类(correction版)"""
    
    def __init__(self, n_groups=2):
        """Initialize方法
        
        Args:
            n_groups: 观测组数量, default2组(pseudorange和Doppler)
        """
        self.n_groups = n_groups
        self.sigma2_0 = np.ones(n_groups)  # 各组初始方差因子
        self.iteration = 0
        self.convergence_threshold = 1e-8
        self.max_iterations = 10
        
    def reset(self):
        """重置方差因子和iteration计数"""
        self.sigma2_0 = np.ones(self.n_groups)
        self.iteration = 0
        
    def estimate(self, residuals_groups, jacobian_groups, weights_groups=None):
        """执line赫尔默特方差分量估计
        
        Args:
            residuals_groups: 各组残差vectorcolumn表
            jacobian_groups: 各组雅可比matrixcolumn表 
            weights_groups: 各组initial权重matrixcolumn表(optional)
            
        Returns:
            updated_weights: Update后权重column表
            converged: whether收敛
            sigma2_values: 各组方差因子
        """
        self.iteration += 1
        
        # Initialize权重, 若未提供则usingunit权阵
        if weights_groups is None:
            weights_groups = [np.ones(len(res)) for res in residuals_groups]
        
        # 构建total雅可比matrix和权重matrix
        # Check residuals_groups, jacobian_groups, weights_groups whether为column表或字典
        if isinstance(residuals_groups, dict):
            all_residuals = np.concatenate([np.asarray(residuals_groups[k]) for k in residuals_groups])
            all_jacobian = np.vstack([np.asarray(jacobian_groups[k]) for k in jacobian_groups])
            all_weights = np.concatenate([np.asarray(weights_groups[k]) for k in weights_groups])
        else:
            all_residuals = np.concatenate([np.asarray(r) for r in residuals_groups])
            all_jacobian = np.vstack([np.asarray(j) for j in jacobian_groups])
            all_weights = np.concatenate([np.asarray(w) for w in weights_groups])
        
        # 构建权重matrix
        W_total = np.diag(all_weights)
        
        # Calculatetotal法方程matrix
        try:
            N_total = all_jacobian.T @ W_total @ all_jacobian
            N_inv = np.linalg.inv(N_total)
        except np.linalg.LinAlgError:
            print("警告: 法方程矩阵奇异, 使用伪逆")
            N_inv = np.linalg.pinv(all_jacobian.T @ W_total @ all_jacobian)
        
        # Calculate各组方差分量
        new_sigma2 = []
        start_idx = 0
        
        for i in range(self.n_groups):
            res_i = residuals_groups[i]
            jac_i = jacobian_groups[i]
            w_i = weights_groups[i]
            n_obs_i = len(res_i)
            
            # 构建该组权重matrix
            W_i = np.diag(w_i)
            
            # Calculate该组冗余度matrix
            # Q_i = W_i^(-1) - W_i^(-1) * J_i * N_inv * J_i^T * W_i^(-1)
            W_i_inv = np.diag(1.0 / (w_i + 1e-12))  # 防止除零
            
            # Calculate冗余度
            try:
                temp = jac_i @ N_inv @ jac_i.T
                Q_i = W_i_inv - W_i_inv @ temp @ W_i_inv
                r_i = np.trace(W_i @ Q_i)  # 第i组冗余度
            except:
                # 简化Calculate
                r_i = max(n_obs_i - jac_i.shape[1], 1.0)

            r_i = max(n_obs_i - jac_i.shape[1], 1.0)
            # 确保冗余度为正
            r_i = max(r_i, 1.0)
            
            # Calculate加权残差平方和
            vtpv_i = res_i.T @ W_i @ res_i
            
            # Calculate方差因子 (标准Helmert公式)
            sigma2_i = vtpv_i / r_i
            new_sigma2.append(sigma2_i)
            
            start_idx += n_obs_i
        
        # Check收敛性
        rel_change = np.max(np.abs(np.array(new_sigma2) - self.sigma2_0) / (self.sigma2_0 + 1e-10))
        converged = (rel_change < self.convergence_threshold) or (self.iteration >= self.max_iterations)
        
        # Update方差因子
        self.sigma2_0 = np.array(new_sigma2)
        
        # Calculate相对权重 (归一化到第一组)
        relative_weights = self.sigma2_0[0] / self.sigma2_0
        
        # according to相对权重Update各组权重
        updated_weights = []
        for i in range(self.n_groups):
            # 权重调整: 相对于第一组权重比例
            updated_w = weights_groups[i] / self.sigma2_0[i]
            updated_weights.append(updated_w)

        return updated_weights, converged, self.sigma2_0

