from PySide6.QtWidgets import (QApplication, QMainWindow, QWidget, QTabWidget, QVBoxLayout, 
    QHBoxLayout, QGridLayout, QLabel, QLineEdit, QPushButton, QFileDialog, QComboBox, QDateTimeEdit,
    QRadioButton, QCheckBox, QGroupBox, QSizePolicy, QFrame, QSpacerItem, QScrollArea,
    QProgressBar, QDialog, QDialogButtonBox)
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtCore import QUrl, Qt, QThread, Signal, QObject
from PySide6.QtGui import QFont, QIcon, QPixmap
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timezone
from core.constellation import generate_walker_delta_tles, write_tle_file
from core.sim_data import SatelliteSimulation
from app.czml import write_czml
from app.settings import load_settings, save_settings, apply_qss, write_web_theme, write_web_config, hex_to_rgba, get_web_data_path
from app.webserver import start as start_web
from tools.sim_obs import run as run_sim_obs
from core.positioning import run_lsq_positioning
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from core.utilities import datetimeUTC_to_gps_time, geodetic_to_ecef

def get_resource_path(relative_path):
    """
    Get absolute path to resource, works for dev and PyInstaller.
    
    Args:
        relative_path: Path relative to app directory (e.g., 'web', 'styles')
    
    Returns:
        Absolute path to the resource
    """
    if getattr(sys, 'frozen', False):
        # Running as compiled executable
        # Resources are in sys._MEIPASS/app/ (not sys._MEIPASS/)
        base_path = os.path.join(sys._MEIPASS, 'app')
    else:
        # Running as script
        base_path = os.path.dirname(os.path.abspath(__file__))
    
    full_path = os.path.join(base_path, relative_path)
    
    # Debug output in frozen mode
    if getattr(sys, 'frozen', False):
        if not os.path.exists(full_path):
            print(f"[WARNING] Resource path not found: {full_path}")
            print(f"  sys._MEIPASS = {sys._MEIPASS}")
            print(f"  base_path = {base_path}")
            print(f"  relative_path = {relative_path}")
        elif relative_path == 'web':
            # Verify web directory structure
            cesium_path = os.path.join(full_path, 'libs', 'cesium', 'Build', 'Cesium')
            if os.path.exists(cesium_path):
                print(f"[INFO] Cesium resources found at: {cesium_path}")
            else:
                print(f"[WARNING] Cesium resources not found at: {cesium_path}")
                # List what's actually in web directory
                try:
                    web_contents = os.listdir(full_path)
                    print(f"  web directory contents: {web_contents[:10]}")
                except:
                    pass
        print(f"  relative_path = {relative_path}")
    
    return full_path

def get_logo_path():
    """
    Get absolute path to logo.png, works for dev and PyInstaller.
    
    Returns:
        Absolute path to logo.png, or None if not found
    """
    if getattr(sys, 'frozen', False):
        # Running as compiled executable
        # Logo should be in sys._MEIPASS/docs/ or sys._MEIPASS/
        possible_paths = [
            os.path.join(sys._MEIPASS, 'docs', 'logo.png'),
            os.path.join(sys._MEIPASS, 'logo.png'),
        ]
    else:
        # Running as script - logo is in docs/ folder at project root
        # Get project root (parent of app directory)
        app_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(app_dir)
        possible_paths = [
            os.path.join(project_root, 'docs', 'logo.png'),
        ]
    
    for path in possible_paths:
        if os.path.exists(path):
            return path
    
    return None

class ConstellationTab(QWidget):
    def __init__(self):
        super().__init__()
        self.sim = SatelliteSimulation()
        self.tle_path = ''
        self._build()

    def _build(self):
        layout = QHBoxLayout(self)
        layout.setSpacing(12)
        layout.setContentsMargins(12, 12, 12, 12)
        
        # Left control panel with scroll area
        left_panel = QWidget()
        left_panel.setMinimumWidth(380)
        left_panel.setMaximumWidth(420)
        
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.NoFrame)
        scroll_content = QWidget()
        left_layout = QVBoxLayout(scroll_content)
        left_layout.setSpacing(6)
        left_layout.setContentsMargins(0, 0, 8, 0)
        
        # Orbit parameters
        orbit_group = QGroupBox('Orbit Parameters')
        orbit_layout = QVBoxLayout(orbit_group)
        orbit_layout.setSpacing(4)
        orbit_layout.setContentsMargins(8, 4, 8, 4)
        
        self.a = QLineEdit('26560000')
        self.e = QLineEdit('0.01')
        self.inc = QLineEdit('55.0')
        self.raan = QLineEdit('0.0')
        self.argp = QLineEdit('0.0')
        self.m0 = QLineEdit('0.0')
        
        for lab, w, unit in [
            ('Semi-major axis', self.a, 'm'),
            ('Eccentricity', self.e, ''),
            ('Inclination', self.inc, '°'),
            ('RAAN', self.raan, '°'),
            ('Arg. of perigee', self.argp, '°'),
            ('Mean anomaly', self.m0, '°')
        ]:
            row = QHBoxLayout()
            label = QLabel(lab)
            label.setMinimumWidth(90)
            row.addWidget(label)
            row.addWidget(w, 1)
            if unit:
                row.addWidget(QLabel(unit))
            orbit_layout.addLayout(row)
        
        left_layout.addWidget(orbit_group)
        
        # Constellation configuration
        walker_group = QGroupBox('Walker-Delta Config')
        walker_layout = QVBoxLayout(walker_group)
        walker_layout.setSpacing(4)
        walker_layout.setContentsMargins(8, 4, 8, 4)
        
        self.P = QLineEdit('3')
        self.S = QLineEdit('4')
        self.F = QLineEdit('1')
        self.raan_spread = QLineEdit('360')
        self.raan_spread.setToolTip('RAAN spread angle (0-360°). 360° for full global coverage.')
        self.epoch = QLineEdit(datetime.utcnow().strftime('%Y-%m-%d'))
        
        # First row: P, S, F
        psf_row = QHBoxLayout()
        psf_row.addWidget(QLabel('P'))
        psf_row.addWidget(self.P)
        psf_row.addWidget(QLabel('S'))
        psf_row.addWidget(self.S)
        psf_row.addWidget(QLabel('F'))
        psf_row.addWidget(self.F)
        walker_layout.addLayout(psf_row)
        
        # RAAN spread
        raan_row = QHBoxLayout()
        raan_row.addWidget(QLabel('RAAN spread'))
        raan_row.addWidget(self.raan_spread, 1)
        raan_row.addWidget(QLabel('°'))
        walker_layout.addLayout(raan_row)
        
        # Epoch
        epoch_row = QHBoxLayout()
        epoch_row.addWidget(QLabel('Epoch'))
        epoch_row.addWidget(self.epoch, 1)
        walker_layout.addLayout(epoch_row)
        
        left_layout.addWidget(walker_group)
        
        # TLE actions
        btn_group = QGroupBox('TLE Actions')
        btn_layout = QHBoxLayout(btn_group)
        btn_layout.setContentsMargins(8, 4, 8, 4)
        
        gen = QPushButton('✨ Generate')
        gen.setObjectName('primaryBtn')
        sav = QPushButton('💾 Save')
        lod = QPushButton('📂 Load')
        
        btn_layout.addWidget(gen)
        btn_layout.addWidget(sav)
        btn_layout.addWidget(lod)
        
        left_layout.addWidget(btn_group)
        
        # Web visualization config
        web_group = QGroupBox('3D Visualization')
        web_layout = QVBoxLayout(web_group)
        web_layout.setSpacing(4)
        web_layout.setContentsMargins(8, 4, 8, 4)
        
        time_row = QHBoxLayout()
        time_row.addWidget(QLabel('Duration'))
        self.web_duration = QLineEdit('10800')
        self.web_duration.setToolTip('>=5400s for full orbit')
        time_row.addWidget(self.web_duration)
        time_row.addWidget(QLabel('s'))
        time_row.addWidget(QLabel('Step'))
        self.web_step = QLineEdit('60')
        self.web_step.setToolTip('30-120s recommended')
        time_row.addWidget(self.web_step)
        time_row.addWidget(QLabel('s'))
        web_layout.addLayout(time_row)
        
        # Reference frame
        ref_frame_row = QHBoxLayout()
        ref_frame_row.addWidget(QLabel('Frame'))
        self.ref_inertial = QRadioButton('Inertial')
        self.ref_inertial.setToolTip('Orbit closes, Earth rotates')
        self.ref_inertial.setStyleSheet("QRadioButton::indicator{width:18px;height:18px;} QRadioButton::indicator:checked{border:1px solid #3b82f6; background:#3b82f6; border-radius:9px;} QRadioButton::indicator:unchecked{border:1px solid #666; background:#111; border-radius:9px;}")
        self.ref_inertial.setChecked(True)
        self.ref_fixed = QRadioButton('Fixed')
        self.ref_fixed.setToolTip('Orbit drifts, Earth fixed')
        self.ref_fixed.setStyleSheet("QRadioButton::indicator{width:18px;height:18px;} QRadioButton::indicator:checked{border:1px solid #3b82f6; background:#3b82f6; border-radius:9px;} QRadioButton::indicator:unchecked{border:1px solid #666; background:#111; border-radius:9px;}")
        ref_frame_row.addWidget(self.ref_inertial)
        ref_frame_row.addWidget(self.ref_fixed)
        ref_frame_row.addStretch()
        web_layout.addLayout(ref_frame_row)

        # Beam options
        beam_row = QHBoxLayout()
        self.chk_beam = QCheckBox('Show Beam')
        self.chk_beam.setToolTip('Display antenna beam cone')
        self.chk_beam.setStyleSheet("QCheckBox::indicator{width:18px;height:18px;} QCheckBox::indicator:checked{border:1px solid #3b82f6; background:#3b82f6; border-radius:9px;} QCheckBox::indicator:unchecked{border:1px solid #666; background:#111; border-radius:9px;}")
        beam_row.addWidget(self.chk_beam)
        beam_row.addWidget(QLabel('Angle'))
        self.beam_angle = QLineEdit('21.0')
        self.beam_angle.setFixedWidth(55)
        beam_row.addWidget(self.beam_angle)
        beam_row.addWidget(QLabel('°'))
        beam_row.addStretch()
        web_layout.addLayout(beam_row)
        
        ref_btn_row = QHBoxLayout()
        ref = QPushButton('🌐 Refresh View')
        ref.setObjectName('primaryBtn')
        open_browser = QPushButton('🚀 Browser')
        open_browser.setToolTip('Open in external browser')
        ref_btn_row.addWidget(ref)
        ref_btn_row.addWidget(open_browser)
        web_layout.addLayout(ref_btn_row)
        
        left_layout.addWidget(web_group)
        
        # Visibility Analysis group
        vis_group = QGroupBox('Visibility Analysis')
        vis_layout = QVBoxLayout(vis_group)
        vis_layout.setSpacing(4)
        vis_layout.setContentsMargins(8, 4, 8, 4)
        
        # Analysis type selection
        vis_type_row = QHBoxLayout()
        self.vis_instant = QRadioButton('Instant')
        self.vis_instant.setToolTip('Single epoch visibility')
        self.vis_instant.setStyleSheet("QRadioButton::indicator{width:18px;height:18px;} QRadioButton::indicator:checked{border:1px solid #3b82f6; background:#3b82f6; border-radius:9px;} QRadioButton::indicator:unchecked{border:1px solid #666; background:#111; border-radius:9px;}")
        self.vis_instant.setChecked(True)
        self.vis_average = QRadioButton('Average')
        self.vis_average.setToolTip('Average over time period')
        self.vis_average.setStyleSheet("QRadioButton::indicator{width:18px;height:18px;} QRadioButton::indicator:checked{border:1px solid #3b82f6; background:#3b82f6; border-radius:9px;} QRadioButton::indicator:unchecked{border:1px solid #666; background:#111; border-radius:9px;}")
        vis_type_row.addWidget(self.vis_instant)
        vis_type_row.addWidget(self.vis_average)
        vis_type_row.addStretch()
        vis_layout.addLayout(vis_type_row)
        
        # Time range for average visibility
        self.vis_time_widget = QWidget()
        vis_time_layout = QVBoxLayout(self.vis_time_widget)
        vis_time_layout.setContentsMargins(0, 0, 0, 0)
        vis_time_layout.setSpacing(2)
        
        from PySide6.QtCore import QDateTime
        start_row = QHBoxLayout()
        start_row.addWidget(QLabel('Start'))
        self.vis_start_dt = QDateTimeEdit()
        self.vis_start_dt.setDisplayFormat('yyyy-MM-dd HH:mm')
        self.vis_start_dt.setDateTime(QDateTime.currentDateTimeUtc())
        start_row.addWidget(self.vis_start_dt, 1)
        vis_time_layout.addLayout(start_row)
        
        end_row = QHBoxLayout()
        end_row.addWidget(QLabel('End'))
        self.vis_end_dt = QDateTimeEdit()
        self.vis_end_dt.setDisplayFormat('yyyy-MM-dd HH:mm')
        self.vis_end_dt.setDateTime(QDateTime.currentDateTimeUtc().addSecs(3600))
        end_row.addWidget(self.vis_end_dt, 1)
        vis_time_layout.addLayout(end_row)
        
        step_row = QHBoxLayout()
        step_row.addWidget(QLabel('Step'))
        self.vis_time_step = QLineEdit('300')
        self.vis_time_step.setToolTip('Time step in seconds')
        step_row.addWidget(self.vis_time_step, 1)
        step_row.addWidget(QLabel('s'))
        vis_time_layout.addLayout(step_row)
        
        self.vis_time_widget.setVisible(False)
        vis_layout.addWidget(self.vis_time_widget)
        
        # Connect radio buttons
        self.vis_instant.toggled.connect(lambda checked: self.vis_time_widget.setVisible(not checked))
        
        # Elevation mask and resolution
        params_row = QHBoxLayout()
        params_row.addWidget(QLabel('Elev'))
        self.vis_elev_mask = QLineEdit('10')
        self.vis_elev_mask.setFixedWidth(40)
        params_row.addWidget(self.vis_elev_mask)
        params_row.addWidget(QLabel('°'))
        params_row.addWidget(QLabel('Grid'))
        self.vis_resolution = QLineEdit('10')
        self.vis_resolution.setFixedWidth(40)
        params_row.addWidget(self.vis_resolution)
        params_row.addWidget(QLabel('°'))
        params_row.addStretch()
        vis_layout.addLayout(params_row)
        
        # Compute button
        self.vis_compute_btn = QPushButton('🌍 Compute Global Coverage')
        self.vis_compute_btn.setObjectName('primaryBtn')
        vis_layout.addWidget(self.vis_compute_btn)
        
        # Progress bar
        self.vis_progress = QProgressBar()
        self.vis_progress.setMinimum(0)
        self.vis_progress.setMaximum(100)
        self.vis_progress.setValue(0)
        self.vis_progress.setTextVisible(True)
        self.vis_progress.setFixedHeight(16)
        vis_layout.addWidget(self.vis_progress)
        
        left_layout.addWidget(vis_group)
        left_layout.addStretch()
        
        scroll.setWidget(scroll_content)
        
        left_container = QVBoxLayout(left_panel)
        left_container.setContentsMargins(0, 0, 0, 0)
        left_container.addWidget(scroll)
        
        # Right-side: Tabbed panel for 3D view and visibility analysis
        right_tabs = QTabWidget()
        right_tabs.setMinimumWidth(600)
        
        # Tab 1: 3D Visualization
        self.view = QWebEngineView()
        right_tabs.addTab(self.view, '🌐 3D Visualization')
        
        # Tab 2: Visibility Analysis
        vis_widget = QWidget()
        vis_tab_layout = QVBoxLayout(vis_widget)
        vis_tab_layout.setContentsMargins(4, 4, 4, 4)
        vis_tab_layout.setSpacing(4)
        
        # Two matplotlib figures for visibility and PDOP
        chart_bg = '#f5f5f7'
        
        self.fig_visibility = Figure(figsize=(8, 4), facecolor=chart_bg)
        self.ax_visibility = self.fig_visibility.add_subplot(111)
        self.canvas_visibility = FigureCanvasQTAgg(self.fig_visibility)
        
        self.fig_pdop = Figure(figsize=(8, 4), facecolor=chart_bg)
        self.ax_pdop = self.fig_pdop.add_subplot(111)
        self.canvas_pdop = FigureCanvasQTAgg(self.fig_pdop)
        
        vis_tab_layout.addWidget(self.canvas_visibility, 1)
        vis_tab_layout.addWidget(self.canvas_pdop, 1)
        
        right_tabs.addTab(vis_widget, '📊 Visibility Analysis')
        
        # Initialize visibility plots
        self._init_visibility_plots()
        
        layout.addWidget(left_panel, 1)
        layout.addWidget(right_tabs, 2)
        self.right_tabs = right_tabs
        
        # Signal connections
        gen.clicked.connect(self.on_generate)
        sav.clicked.connect(self.on_save)
        lod.clicked.connect(self.on_load)
        ref.clicked.connect(self.on_refresh_web)
        open_browser.clicked.connect(self.on_open_browser)
        self.vis_compute_btn.clicked.connect(self.on_compute_visibility)
        self.load_blank()
    
    def _init_visibility_plots(self):
        """Initialize empty visibility plots with proper styling."""
        # Visibility heatmap
        ax = self.ax_visibility
        ax.clear()
        ax.set_facecolor('#ffffff')
        ax.set_title('Global Satellite Visibility', fontsize=11, color='#1f2937', fontweight='bold', pad=10)
        ax.set_xlabel('Longitude (°)', fontsize=9, color='#4b5563')
        ax.set_ylabel('Latitude (°)', fontsize=9, color='#4b5563')
        ax.tick_params(colors='#4b5563', labelsize=8)
        ax.text(0.5, 0.5, 'Generate or load TLE first,\nthen compute visibility', 
                transform=ax.transAxes, ha='center', va='center', fontsize=10, color='#9ca3af')
        self.fig_visibility.tight_layout()
        self.canvas_visibility.draw()
        
        # PDOP heatmap
        ax = self.ax_pdop
        ax.clear()
        ax.set_facecolor('#ffffff')
        ax.set_title('Global Position DOP (PDOP)', fontsize=11, color='#1f2937', fontweight='bold', pad=10)
        ax.set_xlabel('Longitude (°)', fontsize=9, color='#4b5563')
        ax.set_ylabel('Latitude (°)', fontsize=9, color='#4b5563')
        ax.tick_params(colors='#4b5563', labelsize=8)
        ax.text(0.5, 0.5, 'Generate or load TLE first,\nthen compute visibility', 
                transform=ax.transAxes, ha='center', va='center', fontsize=10, color='#9ca3af')
        self.fig_pdop.tight_layout()
        self.canvas_pdop.draw()
    
    def on_compute_visibility(self):
        """Compute and display global visibility analysis."""
        from core.constellation import compute_visibility_analysis
        
        # Check if TLE data is available
        if not hasattr(self, 'tle_data') and not self.tle_path:
            self._init_visibility_plots()
            return
        
        # Get TLE list
        if hasattr(self, 'tle_data') and self.tle_data:
            tle_list = self.tle_data
            print(f"[Visibility] Using generated TLE data: {len(tle_list)} satellites")
        elif self.tle_path:
            from core.sim_data import read_tle_file
            tle_list = read_tle_file(self.tle_path)
            print(f"[Visibility] Loaded TLE from file '{self.tle_path}': {len(tle_list)} satellites")
        else:
            print("[Visibility] No TLE data available!")
            self._init_visibility_plots()
            return
        
        if len(tle_list) == 0:
            print("[Visibility] TLE list is empty!")
            self._init_visibility_plots()
            return
        
        # Get parameters
        try:
            elev_mask = float(self.vis_elev_mask.text())
            resolution = float(self.vis_resolution.text())
        except ValueError:
            elev_mask = 10.0
            resolution = 10.0
        
        # Progress callback
        def progress_cb(percent, msg):
            self.vis_progress.setValue(percent)
            QApplication.processEvents()
        
        self.vis_compute_btn.setEnabled(False)
        self.vis_progress.setValue(0)
        
        try:
            is_average = self.vis_average.isChecked()
            
            if is_average:
                # Average visibility over time period
                from datetime import timedelta
                start_dt = self.vis_start_dt.dateTime().toPython().replace(tzinfo=timezone.utc)
                end_dt = self.vis_end_dt.dateTime().toPython().replace(tzinfo=timezone.utc)
                try:
                    time_step = int(self.vis_time_step.text())
                except ValueError:
                    time_step = 300
                
                print(f"[Visibility] Average mode: {start_dt} to {end_dt}, step={time_step}s")
                
                # Compute for multiple epochs and average
                result = self._compute_average_visibility(
                    tle_list, start_dt, end_dt, time_step,
                    resolution, elev_mask, progress_cb
                )
            else:
                # Instant visibility at epoch
                epoch_dt = datetime.strptime(self.epoch.text(), '%Y-%m-%d').replace(tzinfo=timezone.utc)
                print(f"[Visibility] Instant mode at {epoch_dt}")
                print(f"[Visibility] TLE count: {len(tle_list)}, Resolution: {resolution}°, Elev mask: {elev_mask}°")
                
                result = compute_visibility_analysis(
                    tle_list, epoch_dt,
                    lat_resolution=resolution,
                    lon_resolution=resolution,
                    elevation_mask=elev_mask,
                    progress_callback=progress_cb
                )
            
            # Debug: print results
            vis = result['visibility']
            pdop = result['pdop']
            print(f"[Visibility] Result - Min vis: {np.min(vis)}, Max vis: {np.max(vis)}, Mean: {np.mean(vis):.2f}")
            valid_pdop = pdop[~np.isinf(pdop)]
            if len(valid_pdop) > 0:
                print(f"[Visibility] PDOP - Min: {np.min(valid_pdop):.2f}, Max: {np.max(valid_pdop):.2f}")
            else:
                print("[Visibility] All PDOP values are infinite")
            
            # Plot results
            self._plot_visibility_results(result)
            
            # Switch to visibility tab
            self.right_tabs.setCurrentIndex(1)
            
        except Exception as e:
            print(f"Visibility computation error: {e}")
            import traceback
            traceback.print_exc()
        finally:
            self.vis_compute_btn.setEnabled(True)
            self.vis_progress.setValue(100)
    
    def _compute_average_visibility(self, tle_list, start_dt, end_dt, time_step, resolution, elev_mask, progress_cb):
        """Compute average visibility over a time period."""
        from core.constellation import compute_visibility_analysis
        from datetime import timedelta
        import numpy as np
        
        # Calculate number of epochs
        total_seconds = (end_dt - start_dt).total_seconds()
        n_epochs = max(1, int(total_seconds / time_step) + 1)
        
        print(f"[Visibility] Computing average over {n_epochs} epochs")
        
        # Initialize accumulators
        accumulated_vis = None
        accumulated_pdop = None
        lat_grid = None
        lon_grid = None
        
        for i in range(n_epochs):
            epoch_dt = start_dt + timedelta(seconds=i * time_step)
            
            # Update progress
            if progress_cb:
                progress_cb(int(100 * i / n_epochs), f"Epoch {i+1}/{n_epochs}")
            
            # Compute for this epoch
            result = compute_visibility_analysis(
                tle_list, epoch_dt,
                lat_resolution=resolution,
                lon_resolution=resolution,
                elevation_mask=elev_mask,
                progress_callback=None  # Don't call inner progress
            )
            
            if accumulated_vis is None:
                lat_grid = result['lat_grid']
                lon_grid = result['lon_grid']
                accumulated_vis = result['visibility'].astype(float)
                accumulated_pdop = np.where(np.isinf(result['pdop']), np.nan, result['pdop'])
            else:
                accumulated_vis += result['visibility']
                pdop = np.where(np.isinf(result['pdop']), np.nan, result['pdop'])
                accumulated_pdop = np.nanmean([accumulated_pdop, pdop], axis=0) if i > 0 else pdop
        
        # Compute averages
        avg_vis = accumulated_vis / n_epochs
        avg_pdop = accumulated_pdop
        avg_pdop[np.isnan(avg_pdop)] = np.inf
        
        return {
            'lat_grid': lat_grid,
            'lon_grid': lon_grid,
            'visibility': avg_vis,
            'pdop': avg_pdop,
            'epoch': start_dt,
            'elevation_mask': elev_mask,
            'is_average': True,
            'n_epochs': n_epochs
        }
    
    def _plot_visibility_results(self, result):
        """Plot visibility and PDOP heatmaps."""
        import matplotlib.pyplot as plt
        
        lat_grid = result['lat_grid']
        lon_grid = result['lon_grid']
        visibility = result['visibility']
        pdop = result['pdop']
        elev_mask = result['elevation_mask']
        is_average = result.get('is_average', False)
        n_epochs = result.get('n_epochs', 1)
        
        # Create meshgrid for plotting
        lon_mesh, lat_mesh = np.meshgrid(lon_grid, lat_grid)
        
        # --- Visibility heatmap ---
        ax = self.ax_visibility
        # Remove old colorbar before clearing axes
        if hasattr(self, '_cbar_vis') and self._cbar_vis is not None:
            try:
                self._cbar_vis.remove()
            except Exception:
                pass
            self._cbar_vis = None
        
        ax.clear()
        ax.set_facecolor('#ffffff')
        
        # Use a nice colormap
        cmap = plt.get_cmap('YlGnBu')
        im1 = ax.pcolormesh(lon_mesh, lat_mesh, visibility, cmap=cmap, shading='auto')
        
        # Add colorbar
        self._cbar_vis = self.fig_visibility.colorbar(im1, ax=ax, pad=0.02, shrink=0.9)
        label = 'Avg Visible Sats' if is_average else 'Visible Satellites'
        self._cbar_vis.set_label(label, fontsize=9, color='#4b5563')
        self._cbar_vis.ax.tick_params(labelsize=8, colors='#4b5563')
        
        ax.set_xlim(-180, 180)
        ax.set_ylim(-90, 90)
        title_suffix = f' (Avg over {n_epochs} epochs)' if is_average else ''
        ax.set_title(f'Global Satellite Visibility (Elev ≥ {elev_mask}°){title_suffix}', fontsize=10, color='#1f2937', fontweight='bold', pad=8)
        ax.set_xlabel('Longitude (°)', fontsize=9, color='#4b5563')
        ax.set_ylabel('Latitude (°)', fontsize=9, color='#4b5563')
        ax.tick_params(colors='#4b5563', labelsize=8)
        ax.grid(True, linestyle=':', alpha=0.3, color='#9ca3af')
        
        # Add continent outlines (simplified)
        ax.axhline(0, color='#d1d5db', linewidth=0.5, linestyle='-')
        ax.axvline(0, color='#d1d5db', linewidth=0.5, linestyle='-')
        
        self.fig_visibility.tight_layout()
        self.canvas_visibility.draw()
        
        # --- PDOP heatmap ---
        ax = self.ax_pdop
        
        # Remove old colorbar before clearing axes
        if hasattr(self, '_cbar_pdop') and self._cbar_pdop is not None:
            try:
                self._cbar_pdop.remove()
            except Exception:
                pass
            self._cbar_pdop = None
        
        ax.clear()
        ax.set_facecolor('#ffffff')
        
        # Clip PDOP for better visualization
        pdop_clipped = np.clip(pdop, 0, 10)
        pdop_clipped[np.isinf(pdop)] = 10  # Set inf to max for display
        
        # Use reversed colormap (low PDOP = good = green, high PDOP = bad = red)
        cmap = plt.get_cmap('RdYlGn_r')
        im2 = ax.pcolormesh(lon_mesh, lat_mesh, pdop_clipped, cmap=cmap, shading='auto', vmin=0, vmax=10)
        
        # Add colorbar
        self._cbar_pdop = self.fig_pdop.colorbar(im2, ax=ax, pad=0.02, shrink=0.9)
        self._cbar_pdop.set_label('PDOP (clipped at 10)', fontsize=9, color='#4b5563')
        self._cbar_pdop.ax.tick_params(labelsize=8, colors='#4b5563')
        
        ax.set_xlim(-180, 180)
        ax.set_ylim(-90, 90)
        ax.set_title(f'Global Position DOP (PDOP)', fontsize=11, color='#1f2937', fontweight='bold', pad=10)
        ax.set_xlabel('Longitude (°)', fontsize=9, color='#4b5563')
        ax.set_ylabel('Latitude (°)', fontsize=9, color='#4b5563')
        ax.tick_params(colors='#4b5563', labelsize=8)
        ax.grid(True, linestyle=':', alpha=0.3, color='#9ca3af')
        
        # Add continent outlines (simplified)
        ax.axhline(0, color='#d1d5db', linewidth=0.5, linestyle='-')
        ax.axvline(0, color='#d1d5db', linewidth=0.5, linestyle='-')
        
        self.fig_pdop.tight_layout()
        self.canvas_pdop.draw()

    def on_generate(self):
        a = float(self.a.text())
        e = float(self.e.text())
        inc = float(self.inc.text())
        raan = float(self.raan.text())
        argp = float(self.argp.text())
        m0 = float(self.m0.text())
        P = int(self.P.text())
        S = int(self.S.text())
        F = int(self.F.text())
        raan_spread = float(self.raan_spread.text())
        epoch = datetime.strptime(self.epoch.text(),'%Y-%m-%d').replace(tzinfo=timezone.utc)
        tles = generate_walker_delta_tles('SAT', a, e, inc, raan, argp, m0, epoch, P, S, F, raan_spread=raan_spread)
        self.tle_data = tles
        self.sim = SatelliteSimulation()
        self.sim.tle_satellites = {}
        for i,(name,l1,l2) in enumerate(tles):
            self.sim.add_tle_satellite(name,l1,l2,i+1)

    def on_save(self):
        if not hasattr(self,'tle_data'):
            return
        path, _ = QFileDialog.getSaveFileName(self, 'Save TLE', filter='TLE (*.tle)')
        if path:
            write_tle_file(path, self.tle_data)
            self.tle_path = path

    def on_load(self):
        path, _ = QFileDialog.getOpenFileName(self, 'Load TLE', filter='TLE (*.tle)')
        if path:
            self.sim = SatelliteSimulation()
            self.sim.load_tle_file(path)
            self.tle_path = path
            # Clear generated TLE data so visibility analysis uses loaded file
            if hasattr(self, 'tle_data'):
                delattr(self, 'tle_data')
            self.on_refresh_web()

    def on_refresh_web(self):
        if not self.tle_path and not hasattr(self,'tle_data'):
            return
        tle = self.tle_path
        if not tle:
            tmp = get_web_data_path('temp.tle')
            write_tle_file(tmp, self.tle_data)
            tle = tmp
        start_dt = datetime.utcnow()
        duration_s = int(self.web_duration.text())
        step_s = int(self.web_step.text())
        out_czml = get_web_data_path('data.czml')
        s = load_settings()
        write_web_theme(s.get('accent','#00c8ff'))
        start_iso = start_dt.astimezone(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
        write_web_config(s.get('ion_token',''), s.get('use_world_terrain', False), duration_s, start_iso)
        accent_rgba = hex_to_rgba(s.get('accent','#00c8ff'))
        
        # Beam options
        show_beam = self.chk_beam.isChecked()
        try:
            beam_angle = float(self.beam_angle.text())
        except ValueError:
            beam_angle = 21.0
        
        # Reference frame
        use_fixed = self.ref_fixed.isChecked()
        
        write_czml(tle, start_dt, duration_s, step_s, out_czml, accent_rgba, 
                   show_beam=show_beam, beam_half_angle=beam_angle, use_fixed_frame=use_fixed)
        web_dir = get_resource_path('web')
        # Debug output in frozen mode
        if getattr(sys, 'frozen', False):
            print(f"[GUI] Loading 3D view, web_dir = {web_dir}")
            print(f"[GUI] web_dir exists: {os.path.exists(web_dir)}")
            if os.path.exists(web_dir):
                try:
                    contents = os.listdir(web_dir)
                    print(f"[GUI] web_dir contents: {contents[:10]}")
                    libs_path = os.path.join(web_dir, 'libs')
                    if os.path.exists(libs_path):
                        libs_contents = os.listdir(libs_path)
                        print(f"[GUI] libs/ contents: {libs_contents[:5]}")
                        cesium_path = os.path.join(libs_path, 'cesium', 'Build', 'Cesium')
                        if os.path.exists(cesium_path):
                            print(f"[GUI] Cesium path exists: {cesium_path}")
                        else:
                            print(f"[GUI] ERROR: Cesium path NOT found: {cesium_path}")
                except Exception as e:
                    print(f"[GUI] Error listing web_dir: {e}")
            else:
                print(f"[GUI] ERROR: web_dir does not exist!")
                print(f"[GUI] sys._MEIPASS = {getattr(sys, '_MEIPASS', 'N/A')}")
        # Get data directory (parent of where files are stored)
        if getattr(sys, 'frozen', False):
            import tempfile
            data_dir = os.path.join(tempfile.gettempdir(), 'LEO-SPINE', 'web_data')
        else:
            data_dir = os.path.join(os.path.dirname(__file__), 'web')
        url_base, _srv = start_web(web_dir, data_dir)
        # Pass duration via query param to avoid stale cached config
        self.view.load(QUrl(f"{url_base}/index.html?dur={duration_s}"))

    def load_blank(self):
        web_dir = get_resource_path('web')
        # Get data directory (parent of where files are stored)
        if getattr(sys, 'frozen', False):
            import tempfile
            data_dir = os.path.join(tempfile.gettempdir(), 'LEO-SPINE', 'web_data')
            os.makedirs(data_dir, exist_ok=True)
        else:
            data_dir = os.path.join(os.path.dirname(__file__), 'web')
        
        try:
            url_base, _srv = start_web(web_dir, data_dir)
            self.view.load(QUrl(f"{url_base}/index.html"))
        except Exception as e:
            error_msg = f"Failed to start web server: {e}\nWeb dir: {web_dir}\nData dir: {data_dir}"
            print(f"[GUI] ERROR: {error_msg}")
            if getattr(sys, 'frozen', False):
                html_error = f"""
                <html><body style="font-family: Arial; padding: 20px; background: #1a1a2e; color: #fff;">
                <h2>Web Server Error</h2>
                <p>{error_msg}</p>
                </body></html>
                """
                self.view.setHtml(html_error)
            else:
                raise

    def on_open_browser(self):
        """Open WebGL view in external browser"""
        import webbrowser
        # Read duration consistently to avoid undefined values
        duration_s = int(self.web_duration.text()) if self.web_duration.text() else 5400
        start_dt = datetime.utcnow()
        step_s = int(self.web_step.text()) if self.web_step.text() else 60
        start_iso = start_dt.astimezone(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')

        # Always write config first (timeline extension needs it) even without TLE
        s = load_settings()
        write_web_theme(s.get('accent', '#00c8ff'))
        write_web_config(s.get('ion_token', ''), s.get('use_world_terrain', False), duration_s, start_iso)

        # Refresh data first
        if self.tle_path or hasattr(self, 'tle_data'):
            tle = self.tle_path
            if not tle:
                tmp = get_web_data_path('temp.tle')
                write_tle_file(tmp, self.tle_data)
                tle = tmp
            out_czml = get_web_data_path('data.czml')
            write_web_config(s.get('ion_token', ''), s.get('use_world_terrain', False), duration_s, start_iso)
            accent_rgba = hex_to_rgba(s.get('accent', '#00c8ff'))
            write_czml(tle, start_dt, duration_s, step_s, out_czml, accent_rgba)
        
        # Start web server and open in browser
        web_dir = get_resource_path('web')
        # Get data directory (parent of where files are stored)
        if getattr(sys, 'frozen', False):
            import tempfile
            data_dir = os.path.join(tempfile.gettempdir(), 'LEO-SPINE', 'web_data')
        else:
            data_dir = os.path.join(os.path.dirname(__file__), 'web')
        url_base, _srv = start_web(web_dir, data_dir)
        webbrowser.open(f"{url_base}/index.html?dur={duration_s}")

class SimulationTab(QWidget):
    def __init__(self):
        super().__init__()
        root = QHBoxLayout(self)
        root.setSpacing(12)
        root.setContentsMargins(12, 12, 12, 12)
        
        # Left parameter panel with scroll area
        left_panel = QWidget()
        left_panel.setMinimumWidth(420)
        left_panel.setMaximumWidth(480)
        
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.NoFrame)
        scroll_content = QWidget()
        left = QVBoxLayout(scroll_content)
        left.setSpacing(6)
        left.setContentsMargins(0, 0, 8, 0)
        
        # TLE file selection
        tle_group = QGroupBox('TLE File')
        tle_layout = QHBoxLayout(tle_group)
        tle_layout.setContentsMargins(8, 4, 8, 4)
        self.tle = QLineEdit()
        self.tle.setPlaceholderText('Select a TLE file...')
        btn = QPushButton('...')
        btn.setFixedWidth(36)
        btn.clicked.connect(self.pick_tle)
        tle_layout.addWidget(self.tle, 1)
        tle_layout.addWidget(btn)
        left.addWidget(tle_group)
        
        # Satellite selection
        sat_group = QGroupBox('Satellite Selection')
        sat_layout = QVBoxLayout(sat_group)
        sat_layout.setSpacing(4)
        sat_layout.setContentsMargins(8, 4, 8, 4)
        
        # Select all / deselect all buttons
        sat_btn_row = QHBoxLayout()
        self.btn_select_all = QPushButton('Select All')
        self.btn_select_all.setFixedHeight(24)
        self.btn_deselect_all = QPushButton('Deselect All')
        self.btn_deselect_all.setFixedHeight(24)
        self.sat_count_label = QLabel('0 satellites loaded')
        self.sat_count_label.setStyleSheet('color: #888;')
        sat_btn_row.addWidget(self.btn_select_all)
        sat_btn_row.addWidget(self.btn_deselect_all)
        sat_btn_row.addWidget(self.sat_count_label)
        sat_btn_row.addStretch()
        sat_layout.addLayout(sat_btn_row)
        
        # Satellite list with checkboxes
        from PySide6.QtWidgets import QListWidget, QListWidgetItem
        self.sat_list = QListWidget()
        self.sat_list.setMaximumHeight(120)
        self.sat_list.setSelectionMode(QListWidget.MultiSelection)
        self.sat_list.setStyleSheet('''
            QListWidget { background: #1a1a2e; border: 1px solid #333; border-radius: 4px; }
            QListWidget::item { padding: 2px 4px; }
            QListWidget::item:selected { background: #3b82f6; color: white; }
        ''')
        sat_layout.addWidget(self.sat_list)
        
        # Connect buttons
        self.btn_select_all.clicked.connect(self._select_all_sats)
        self.btn_deselect_all.clicked.connect(self._deselect_all_sats)
        
        left.addWidget(sat_group)
        
        # Time configuration
        time_group = QGroupBox('Time Configuration')
        time_layout = QVBoxLayout(time_group)
        time_layout.setSpacing(4)
        time_layout.setContentsMargins(8, 4, 8, 4)
        
        from PySide6.QtCore import QDateTime
        self.start_dt = QDateTimeEdit()
        self.start_dt.setDisplayFormat('yyyy-MM-dd HH:mm:ss')
        self.start_dt.setTimeSpec(Qt.UTC)
        self.start_dt.setDateTime(QDateTime.currentDateTimeUtc())
        self.start_dt.setCalendarPopup(True)
        
        time_row = QHBoxLayout()
        time_row.addWidget(QLabel('Start (UTC)'))
        time_row.addWidget(self.start_dt, 1)
        btn_now = QPushButton('Now')
        btn_now.setFixedWidth(40)
        def _set_now():
            self.start_dt.setDateTime(QDateTime.currentDateTimeUtc())
        btn_now.clicked.connect(_set_now)
        time_row.addWidget(btn_now)
        time_layout.addLayout(time_row)
        
        self.duration = QLineEdit('1000')
        self.step = QLineEdit('1')
        
        dur_step_row = QHBoxLayout()
        dur_step_row.addWidget(QLabel('Duration'))
        dur_step_row.addWidget(self.duration, 1)
        dur_step_row.addWidget(QLabel('s'))
        dur_step_row.addWidget(QLabel('Step'))
        dur_step_row.addWidget(self.step, 1)
        dur_step_row.addWidget(QLabel('s'))
        time_layout.addLayout(dur_step_row)
        
        left.addWidget(time_group)
        
        # Station configuration
        station_group = QGroupBox('Station Position')
        station_layout = QVBoxLayout(station_group)
        station_layout.setSpacing(4)
        station_layout.setContentsMargins(8, 4, 8, 4)
        
        # Coordinate type selection
        coord_type_row = QHBoxLayout()
        radio_style = "QRadioButton::indicator{width:16px;height:16px;}"
        self.coord_llh_radio = QRadioButton('LLH')
        self.coord_llh_radio.setStyleSheet(radio_style + " QRadioButton::indicator:checked{border:1px solid #3b82f6; background:#3b82f6; border-radius:9px;} QRadioButton::indicator:unchecked{border:1px solid #666; background:#111; border-radius:9px;}")
        self.coord_xyz_radio = QRadioButton('XYZ')
        self.coord_xyz_radio.setStyleSheet(radio_style + " QRadioButton::indicator:checked{border:1px solid #3b82f6; background:#3b82f6; border-radius:9px;} QRadioButton::indicator:unchecked{border:1px solid #666; background:#111; border-radius:9px;}")
        self.coord_llh_radio.setChecked(True)
        coord_type_row.addWidget(self.coord_llh_radio)
        coord_type_row.addWidget(self.coord_xyz_radio)
        coord_type_row.addStretch()
        station_layout.addLayout(coord_type_row)
        
        # Lat/Lon/Alt input
        self.lat = QLineEdit('39.8')
        self.lon = QLineEdit('119.29')
        self.alt = QLineEdit('60')
        
        llh_row = QHBoxLayout()
        llh_row.addWidget(QLabel('Lat'))
        llh_row.addWidget(self.lat, 1)
        llh_row.addWidget(QLabel('°'))
        llh_row.addWidget(QLabel('Lon'))
        llh_row.addWidget(self.lon, 1)
        llh_row.addWidget(QLabel('°'))
        llh_row.addWidget(QLabel('Alt'))
        llh_row.addWidget(self.alt, 1)
        llh_row.addWidget(QLabel('m'))
        station_layout.addLayout(llh_row)
        
        # ECEF XYZ input
        self.x = QLineEdit()
        self.y = QLineEdit()
        self.z = QLineEdit()
        
        xyz_row = QHBoxLayout()
        xyz_row.addWidget(QLabel('X'))
        xyz_row.addWidget(self.x, 1)
        xyz_row.addWidget(QLabel('Y'))
        xyz_row.addWidget(self.y, 1)
        xyz_row.addWidget(QLabel('Z'))
        xyz_row.addWidget(self.z, 1)
        xyz_row.addWidget(QLabel('m'))
        station_layout.addLayout(xyz_row)
        
        # Elevation mask
        self.mask = QLineEdit('0')
        mask_row = QHBoxLayout()
        mask_row.addWidget(QLabel('Elevation mask'))
        mask_row.addWidget(self.mask, 1)
        mask_row.addWidget(QLabel('°'))
        station_layout.addLayout(mask_row)
        
        left.addWidget(station_group)
        
        # Error Models
        error_group = QGroupBox('Error Models')
        error_layout = QVBoxLayout(error_group)
        error_layout.setSpacing(4)
        error_layout.setContentsMargins(8, 4, 8, 4)
        
        # First row of error checkboxes
        err_row1 = QHBoxLayout()
        self.chk_iono = QCheckBox('Ionosphere')
        self.chk_tropo = QCheckBox('Troposphere')
        self.chk_iono.setChecked(True)
        self.chk_tropo.setChecked(True)
        err_row1.addWidget(self.chk_iono)
        err_row1.addWidget(self.chk_tropo)
        err_row1.addStretch()
        error_layout.addLayout(err_row1)
        
        # Second row of error checkboxes
        err_row2 = QHBoxLayout()
        self.chk_multipath = QCheckBox('Multipath')
        self.chk_relativistic = QCheckBox('Relativistic')
        err_row2.addWidget(self.chk_multipath)
        err_row2.addWidget(self.chk_relativistic)
        err_row2.addStretch()
        error_layout.addLayout(err_row2)
        
        # Third row
        err_row3 = QHBoxLayout()
        self.chk_hardware = QCheckBox('Hardware Delay')
        err_row3.addWidget(self.chk_hardware)
        err_row3.addStretch()
        error_layout.addLayout(err_row3)
        
        left.addWidget(error_group)
        
        # Receiver Clock
        clock_group = QGroupBox('Receiver Clock')
        clock_layout = QVBoxLayout(clock_group)
        clock_layout.setSpacing(4)
        clock_layout.setContentsMargins(8, 4, 8, 4)
        
        clk_src_row = QHBoxLayout()
        self.clk_model_radio = QRadioButton('Model')
        self.clk_file_radio = QRadioButton('From File')
        self.clk_model_radio.setStyleSheet("QRadioButton::indicator{width:18px;height:18px;} QRadioButton::indicator:checked{border:1px solid #3b82f6; background:#3b82f6; border-radius:9px;} QRadioButton::indicator:unchecked{border:1px solid #666; background:#111; border-radius:9px;}")
        self.clk_file_radio.setStyleSheet("QRadioButton::indicator{width:18px;height:18px;} QRadioButton::indicator:checked{border:1px solid #3b82f6; background:#3b82f6; border-radius:9px;} QRadioButton::indicator:unchecked{border:1px solid #666; background:#111; border-radius:9px;}")
        self.clk_model_radio.setChecked(True)
        clk_src_row.addWidget(self.clk_model_radio)
        clk_src_row.addWidget(self.clk_file_radio)
        clk_src_row.addStretch()
        clock_layout.addLayout(clk_src_row)
        
        # Clock model parameters
        clk_model_row = QHBoxLayout()
        clk_model_row.addWidget(QLabel('Bias'))
        self.clk_bias = QLineEdit('1000')
        clk_model_row.addWidget(self.clk_bias, 1)
        clk_model_row.addWidget(QLabel('m'))
        clk_model_row.addWidget(QLabel('Drift'))
        self.clk_drift = QLineEdit('0.05')
        clk_model_row.addWidget(self.clk_drift, 1)
        clk_model_row.addWidget(QLabel('m/s'))
        clock_layout.addLayout(clk_model_row)
        
        # Clock file selection
        clk_file_row = QHBoxLayout()
        self.clk_file = QLineEdit()
        self.clk_file.setPlaceholderText('Clock bias file...')
        btn_clk = QPushButton('📂')
        btn_clk.setFixedWidth(36)
        btn_clk.clicked.connect(self.pick_clk_file)
        clk_file_row.addWidget(self.clk_file, 1)
        clk_file_row.addWidget(btn_clk)
        clock_layout.addLayout(clk_file_row)
        
        clk_drift_file_row = QHBoxLayout()
        self.clk_drift_file = QLineEdit()
        self.clk_drift_file.setPlaceholderText('Clock drift file...')
        btn_clk_drift = QPushButton('📂')
        btn_clk_drift.setFixedWidth(36)
        btn_clk_drift.clicked.connect(self.pick_clk_drift_file)
        clk_drift_file_row.addWidget(self.clk_drift_file, 1)
        clk_drift_file_row.addWidget(btn_clk_drift)
        clock_layout.addLayout(clk_drift_file_row)
        
        left.addWidget(clock_group)
        
        # Measurement Noise
        noise_group = QGroupBox('Measurement Noise (Gaussian)')
        noise_layout = QVBoxLayout(noise_group)
        noise_layout.setSpacing(4)
        noise_layout.setContentsMargins(8, 4, 8, 4)
        
        pr_noise_row = QHBoxLayout()
        pr_noise_row.addWidget(QLabel('Pseudorange σ'))
        self.pr_noise = QLineEdit('1.0')
        pr_noise_row.addWidget(self.pr_noise, 1)
        pr_noise_row.addWidget(QLabel('m'))
        noise_layout.addLayout(pr_noise_row)
        
        dop_noise_row = QHBoxLayout()
        dop_noise_row.addWidget(QLabel('Doppler σ'))
        self.dop_noise = QLineEdit('0.01')
        dop_noise_row.addWidget(self.dop_noise, 1)
        dop_noise_row.addWidget(QLabel('m/s'))
        noise_layout.addLayout(dop_noise_row)
        
        left.addWidget(noise_group)
        
        # Output configuration
        output_group = QGroupBox('Output File')
        output_layout = QHBoxLayout(output_group)
        output_layout.setContentsMargins(8, 4, 8, 4)
        self.out = QLineEdit('observations.csv')
        self.out.setPlaceholderText('Select save path...')
        btn_out = QPushButton('💾')
        btn_out.setFixedWidth(36)
        btn_out.clicked.connect(self.pick_output)
        output_layout.addWidget(self.out, 1)
        output_layout.addWidget(btn_out)
        left.addWidget(output_group)
        
        # Run button
        runbtn = QPushButton('🚀 Start Simulation')
        runbtn.setObjectName('primaryBtn')
        runbtn.clicked.connect(self.run_sim)
        left.addWidget(runbtn)
        
        left.addStretch()
        
        # Wire signals for live coordinate conversion
        self.coord_llh_radio.toggled.connect(self._on_coord_type_changed)
        self.coord_xyz_radio.toggled.connect(self._on_coord_type_changed)
        self.lat.textChanged.connect(self._on_llh_changed)
        self.lon.textChanged.connect(self._on_llh_changed)
        self.alt.textChanged.connect(self._on_llh_changed)
        self.x.textChanged.connect(self._on_xyz_changed)
        self.y.textChanged.connect(self._on_xyz_changed)
        self.z.textChanged.connect(self._on_xyz_changed)
        
        # Wire clock source signals
        self.clk_model_radio.toggled.connect(self._on_clk_source_changed)
        self.clk_file_radio.toggled.connect(self._on_clk_source_changed)
        
        # Initialize states
        self._on_coord_type_changed()
        self._on_llh_changed()
        self._on_clk_source_changed()
        
        scroll.setWidget(scroll_content)
        
        left_container = QVBoxLayout(left_panel)
        left_container.setContentsMargins(0, 0, 0, 0)
        left_container.addWidget(scroll)
        
        # Right chart area
        right = QVBoxLayout()
        right.setSpacing(8)
        
        # MacOS-style light gray background
        chart_bg = '#f5f5f7'
        axis_bg = '#ffffff'
        
        self.fig1 = Figure(figsize=(5, 3), facecolor=chart_bg)
        self.ax_pr = self.fig1.add_subplot(111)
        self.ax_pr.set_facecolor(axis_bg)
        self.canvas1 = FigureCanvasQTAgg(self.fig1)
        right.addWidget(self.canvas1)
        
        self.fig2 = Figure(figsize=(5, 3), facecolor=chart_bg)
        self.ax_dp = self.fig2.add_subplot(111)
        self.ax_dp.set_facecolor(axis_bg)
        self.canvas2 = FigureCanvasQTAgg(self.fig2)
        right.addWidget(self.canvas2)
        
        # Initialize chart styles
        self._style_charts()
        
        root.addWidget(left_panel, 1)
        root.addLayout(right, 2)

    def pick_tle(self):
        path, _ = QFileDialog.getOpenFileName(self, 'Select TLE', filter='TLE (*.tle *.txt);;All Files (*.*)')
        if path:
            self.tle.setText(path)
            self._load_satellite_list(path)
    
    def _load_satellite_list(self, tle_path):
        """Parse TLE file and populate satellite selection list."""
        from PySide6.QtWidgets import QListWidgetItem
        self.sat_list.clear()
        
        try:
            satellites = []
            with open(tle_path, 'r', encoding='utf-8') as f:
                lines = [l.strip() for l in f.readlines() if l.strip()]
            
            i = 0
            prn = 1
            while i < len(lines):
                # Try to detect TLE format
                if i + 2 < len(lines) and lines[i+1].startswith('1 ') and lines[i+2].startswith('2 '):
                    # Three-line TLE format
                    name = lines[i]
                    satellites.append((prn, name))
                    prn += 1
                    i += 3
                elif lines[i].startswith('1 ') and i + 1 < len(lines) and lines[i+1].startswith('2 '):
                    # Two-line TLE format (no name line)
                    # Extract satellite number from line 1
                    try:
                        sat_num = lines[i][2:7].strip()
                        name = f'SAT-{sat_num}'
                    except:
                        name = f'SAT-{prn}'
                    satellites.append((prn, name))
                    prn += 1
                    i += 2
                else:
                    i += 1
            
            # Populate list
            for sat_prn, sat_name in satellites:
                item = QListWidgetItem(f'{sat_prn}: {sat_name}')
                item.setData(Qt.UserRole, sat_prn)  # Store PRN
                item.setSelected(True)  # Select by default
                self.sat_list.addItem(item)
            
            # Update count label
            self.sat_count_label.setText(f'{len(satellites)} satellites loaded')
            
        except Exception as e:
            self.sat_count_label.setText(f'Error: {str(e)[:30]}')
    
    def _select_all_sats(self):
        """Select all satellites in the list."""
        for i in range(self.sat_list.count()):
            self.sat_list.item(i).setSelected(True)
    
    def _deselect_all_sats(self):
        """Deselect all satellites in the list."""
        for i in range(self.sat_list.count()):
            self.sat_list.item(i).setSelected(False)
    
    def _get_selected_prns(self):
        """Get list of selected satellite PRNs."""
        selected = []
        for i in range(self.sat_list.count()):
            item = self.sat_list.item(i)
            if item.isSelected():
                prn = item.data(Qt.UserRole)
                selected.append(prn)
        return selected if selected else None  # None means all satellites
    
    def pick_output(self):
        path, _ = QFileDialog.getSaveFileName(self, 'Save observation file', 'observations.csv', 
                                              filter='CSV Files (*.csv);;OBS Files (*.obs);;All Files (*.*)')
        if path:
            self.out.setText(path)
    
    def pick_clk_file(self):
        path, _ = QFileDialog.getOpenFileName(self, 'Select Clock Bias File', filter='Text Files (*.txt);;All Files (*.*)')
        if path:
            self.clk_file.setText(path)
    
    def pick_clk_drift_file(self):
        path, _ = QFileDialog.getOpenFileName(self, 'Select Clock Drift File', filter='Text Files (*.txt);;All Files (*.*)')
        if path:
            self.clk_drift_file.setText(path)
    
    def _on_clk_source_changed(self):
        """Handle clock source switching"""
        is_model = self.clk_model_radio.isChecked()
        self.clk_bias.setEnabled(is_model)
        self.clk_drift.setEnabled(is_model)
        self.clk_file.setEnabled(not is_model)
        self.clk_drift_file.setEnabled(not is_model)
    
    def _on_coord_type_changed(self):
        """Handle coordinate type switching"""
        is_llh = self.coord_llh_radio.isChecked()
        # Enable/disable input fields
        self.lat.setEnabled(is_llh)
        self.lon.setEnabled(is_llh)
        self.alt.setEnabled(is_llh)
        self.x.setEnabled(not is_llh)
        self.y.setEnabled(not is_llh)
        self.z.setEnabled(not is_llh)
        
        # If switching to LLH and XYZ is available, convert
        if is_llh and self.x.text() and self.y.text() and self.z.text():
            self._on_xyz_changed()
        # If switching to XYZ and LLH is available, convert
        elif not is_llh and self.lat.text() and self.lon.text() and self.alt.text():
            self._on_llh_changed()
    
    def _on_llh_changed(self):
        """Update XYZ when lat/lon/alt changes"""
        if not self.coord_llh_radio.isChecked():
            return
        
        try:
            lat = float(self.lat.text())
            lon = float(self.lon.text())
            alt = float(self.alt.text())
            
            # Import conversion helpers
            from core.utilities import geodetic_to_ecef
            import numpy as np
            
            # Convert to radians
            lat_rad = np.deg2rad(lat)
            lon_rad = np.deg2rad(lon)
            
            # Convert to ECEF
            x, y, z = geodetic_to_ecef(lat_rad, lon_rad, alt)
            
            # Update XYZ display (block reverse trigger)
            self.x.blockSignals(True)
            self.y.blockSignals(True)
            self.z.blockSignals(True)
            self.x.setText(f'{x:.2f}')
            self.y.setText(f'{y:.2f}')
            self.z.setText(f'{z:.2f}')
            self.x.blockSignals(False)
            self.y.blockSignals(False)
            self.z.blockSignals(False)
        except Exception:
            # Skip if invalid input
            pass
    
    def _on_xyz_changed(self):
        """Update lat/lon/alt when XYZ changes"""
        if not self.coord_xyz_radio.isChecked():
            return
        
        try:
            x = float(self.x.text())
            y = float(self.y.text())
            z = float(self.z.text())
            
            # Import conversion helpers
            from core.utilities import ecef_to_geodetic
            
            # Convert to geodetic
            lat_rad, lon_rad, alt = ecef_to_geodetic(x, y, z)
            
            # Convert to degrees
            import numpy as np
            lat = np.rad2deg(lat_rad)
            lon = np.rad2deg(lon_rad)
            
            # Update LLH display (block reverse trigger)
            self.lat.blockSignals(True)
            self.lon.blockSignals(True)
            self.alt.blockSignals(True)
            self.lat.setText(f'{lat:.6f}')
            self.lon.setText(f'{lon:.6f}')
            self.alt.setText(f'{alt:.2f}')
            self.lat.blockSignals(False)
            self.lon.blockSignals(False)
            self.alt.blockSignals(False)
        except Exception:
            # Skip if invalid input
            pass

    def _style_charts(self):
        """Apply macOS-style chart theme"""
        for ax in [self.ax_pr, self.ax_dp]:
            ax.set_facecolor('#ffffff')
            ax.tick_params(colors='#333333', labelsize=9)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_color('#cccccc')
            ax.spines['left'].set_color('#cccccc')
            ax.grid(True, linestyle='--', alpha=0.3, color='#999999')
        self.canvas1.draw()
        self.canvas2.draw()

    def run_sim(self):
        try:
            dt = self.start_dt.dateTime().toUTC().toPython()
        except Exception:
            from datetime import datetime, timezone
            dt = datetime.utcnow().replace(tzinfo=timezone.utc)
        gps_week, gps_seconds = datetimeUTC_to_gps_time(dt.replace(tzinfo=None))
        
        # Get station coordinates
        if self.coord_llh_radio.isChecked():
            # Using lat/lon/alt input
            lat = float(self.lat.text())
            lon = float(self.lon.text())
            alt = float(self.alt.text())
        else:
            # Using ECEF XYZ input
            x = float(self.x.text())
            y = float(self.y.text())
            z = float(self.z.text())
            # Convert XYZ to geodetic
            from core.utilities import ecef_to_geodetic
            lat_rad, lon_rad, alt = ecef_to_geodetic(x, y, z)
            import numpy as np
            lat = np.rad2deg(lat_rad)
            lon = np.rad2deg(lon_rad)
        
        # Collect error model settings
        enable_iono = self.chk_iono.isChecked()
        enable_tropo = self.chk_tropo.isChecked()
        enable_multipath = self.chk_multipath.isChecked()
        enable_relativistic = self.chk_relativistic.isChecked()
        enable_hardware = self.chk_hardware.isChecked()
        
        # Collect measurement noise settings
        try:
            pr_noise = float(self.pr_noise.text())
        except ValueError:
            pr_noise = 1.0
        try:
            dop_noise = float(self.dop_noise.text())
        except ValueError:
            dop_noise = 0.01
        
        # Collect receiver clock settings
        use_clk_model = self.clk_model_radio.isChecked()
        if use_clk_model:
            # Use model with bias and drift values
            try:
                clk_bias = float(self.clk_bias.text())
            except ValueError:
                clk_bias = 1000.0
            try:
                clk_drift = float(self.clk_drift.text())
            except ValueError:
                clk_drift = 0.05
            clk_file = None
            clk_drift_file = None
        else:
            # Use files
            clk_bias = 1000.0  # default, won't be used
            clk_drift = 0.05  # default, won't be used
            clk_file = self.clk_file.text() if self.clk_file.text() else None
            clk_drift_file = self.clk_drift_file.text() if self.clk_drift_file.text() else None
        
        # Get selected satellites
        selected_prns = self._get_selected_prns()
        
        # Call simulation with all parameters
        run_sim_obs(
            self.tle.text(), gps_seconds, float(self.duration.text()), float(self.step.text()),
            lat, lon, alt, float(self.mask.text()), self.out.text(), gps_week,
            prnoise=pr_noise, dopnoise=dop_noise,
            enable_iono=enable_iono, enable_tropo=enable_tropo, enable_multipath=enable_multipath,
            enable_relativistic=enable_relativistic, enable_hardware=enable_hardware,
            use_clk_model=use_clk_model, clk_bias=clk_bias, clk_drift=clk_drift,
            clk_file=clk_file, clk_drift_file=clk_drift_file,
            selected_prns=selected_prns
        )
        try:
            import csv
            times = {}
            pr_map = {}
            dp_map = {}
            with open(self.out.text(), 'r', encoding='utf-8') as f:
                first = ''
                for line in f:
                    if line.strip():
                        first = line.strip()
                        break
                f.seek(0)
                if first.startswith('OBS:') or first.startswith('BS:'):
                    for s in f:
                        s = s.strip()
                        if not s or not (s.startswith('OBS:') or s.startswith('BS:')):
                            continue
                        fields = [x for x in s.split() if x]
                        if len(fields) < 12:
                            continue
                        try:
                            te = float(fields[2])
                            pr = float(fields[4])
                            dp = float(fields[5])
                            ele = float(fields[10])
                            prn = fields[11]
                        except Exception:
                            continue
                        times.setdefault(prn, []).append(te)
                        pr_map.setdefault(prn, []).append(pr)
                        dp_map.setdefault(prn, []).append(dp)
                else:
                    r = csv.DictReader(f)
                    for row in r:
                        te = float(row.get('emission_time', row.get('time_s', '0')))
                        pr = float(row.get('pseudorange', row.get('pseudorange_m', '0')))
                        dp = float(row.get('doppler', row.get('doppler_Hz', '0')))
                        prn = row.get('prn')
                        if prn is None:
                            continue
                        times.setdefault(prn, []).append(te)
                        pr_map.setdefault(prn, []).append(pr)
                        dp_map.setdefault(prn, []).append(dp)
            self.ax_pr.clear()
            self.ax_pr.set_facecolor('#ffffff')
            for prn, ts in times.items():
                self.ax_pr.plot(ts, pr_map.get(prn, []), label=f"PRN {prn}", linewidth=1.5)
            self.ax_pr.set_ylabel('Pseudorange(m)', fontsize=10, color='#333333')
            self.ax_pr.set_xlabel('Time(s)', fontsize=10, color='#333333')
            self.ax_pr.legend(loc='upper right', fontsize=8, framealpha=0.9)
            self.ax_pr.tick_params(colors='#333333', labelsize=9)
            self.ax_pr.spines['top'].set_visible(False)
            self.ax_pr.spines['right'].set_visible(False)
            self.ax_pr.spines['bottom'].set_color('#cccccc')
            self.ax_pr.spines['left'].set_color('#cccccc')
            self.ax_pr.grid(True, linestyle='--', alpha=0.3, color='#999999')
            self.canvas1.draw()
            self.ax_dp.clear()
            self.ax_dp.set_facecolor('#ffffff')
            for prn, ts in times.items():
                self.ax_dp.plot(ts, dp_map.get(prn, []), label=f"PRN {prn}", linewidth=1.5)
            self.ax_dp.set_ylabel('Doppler(Hz)', fontsize=10, color='#333333')
            self.ax_dp.set_xlabel('Time(s)', fontsize=10, color='#333333')
            self.ax_dp.legend(loc='upper right', fontsize=8, framealpha=0.9)
            self.ax_dp.tick_params(colors='#333333', labelsize=9)
            self.ax_dp.spines['top'].set_visible(False)
            self.ax_dp.spines['right'].set_visible(False)
            self.ax_dp.spines['bottom'].set_color('#cccccc')
            self.ax_dp.spines['left'].set_color('#cccccc')
            self.ax_dp.grid(True, linestyle='--', alpha=0.3, color='#999999')
            self.canvas2.draw()
        except Exception:
            pass


class PositioningWorker(QObject):
    """Positioning computation worker thread"""
    progress = Signal(int, str)  # progress percent, status text
    finished = Signal(object, object)  # completion signal: convergence_data, results
    error = Signal(str)  # error signal
    
    def __init__(self, obs_file, eph_file, tle_file, init_pos, strategy, iono=True, tropo=True, multipath=False, relativistic=False, hardware=False, use_hvce=False, selected_prns=None):
        super().__init__()
        self.obs_file = obs_file
        self.eph_file = eph_file
        self.tle_file = tle_file
        self.init_pos = init_pos
        self.strategy = strategy  # 'pr', 'dop', 'fixed', 'adaptive'
        self.iono = iono
        self.tropo = tropo
        self.multipath = multipath
        self.relativistic = relativistic
        self.hardware = hardware
        self.pco = False
        self.use_hvce = use_hvce
        self.selected_prns = selected_prns  # List of PRNs to use, None = all
        self._is_cancelled = False
    
    def cancel(self):
        self._is_cancelled = True
    
    def run(self):
        try:
            self.progress.emit(5, "Loading observation data...")
            
            # Determine use_pseudorange and use_doppler based on strategy
            if self.strategy == 'pr':
                use_pr, use_dop, use_hvce = True, False, False
            elif self.strategy == 'dop':
                use_pr, use_dop, use_hvce = False, True, False
            elif self.strategy == 'fixed':
                use_pr, use_dop, use_hvce = True, True, False
            else:  # adaptive
                use_pr, use_dop, use_hvce = True, True, True
            

            # Call positioning with progress callback
            convergence_data, results = run_lsq_positioning(
                self.obs_file, 
                ephemeris_file=self.eph_file, 
                tle_file=self.tle_file, 
                init_pos=self.init_pos, 
                use_pseudorange=use_pr, 
                use_doppler=use_dop,
                use_hvce=use_hvce, 
                enable_iono=self.iono,
                enable_tropo=self.tropo,
                enable_multipath=self.multipath,
                enable_relativistic=self.relativistic,
                enable_hardware=self.hardware,
                enable_pco=self.pco,
                progress_callback=self._on_progress,
                selected_prns=self.selected_prns
            )
            
            if not self._is_cancelled:
                self.progress.emit(100, "Positioning completed!")
                self.finished.emit(convergence_data, results)
        except Exception as e:
            import traceback
            self.error.emit(f"Positioning failed: {str(e)}\n{traceback.format_exc()}")
    
    def _on_progress(self, current, total, message):
        """Progress callback"""
        if self._is_cancelled:
            return False  # Return False to cancel
        percent = int(10 + (current / total) * 85) if total > 0 else 10
        self.progress.emit(percent, message)
        return True

class PositioningTab(QWidget):
    def __init__(self):
        super().__init__()
        self.worker = None
        self.thread = None
        
        main_layout = QHBoxLayout(self)
        main_layout.setSpacing(12)
        main_layout.setContentsMargins(12, 12, 12, 12)
        
        # Left parameter panel
        left_panel = QWidget()
        left_panel.setMinimumWidth(450)
        layout = QVBoxLayout(left_panel)
        layout.setSpacing(8)
        layout.setContentsMargins(0, 0, 0, 0)
        
        # Scroll area
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.NoFrame)
        content = QWidget()
        scroll_layout = QVBoxLayout(content)
        scroll_layout.setSpacing(8)
        
        # Observation file selection
        obs_group = QGroupBox('Observation Data')
        obs_main_layout = QVBoxLayout(obs_group)
        obs_main_layout.setContentsMargins(8, 6, 8, 6)
        obs_main_layout.setSpacing(6)
        
        obs_file_row = QHBoxLayout()
        self.obs = QLineEdit('observations.csv')
        self.obs.setPlaceholderText('Select observation file...')
        btn_obs = QPushButton('Browse')
        btn_obs.clicked.connect(self.pick_obs)
        for w in (self.obs, btn_obs):
            w.setFixedHeight(30)
        obs_file_row.addWidget(self.obs, 1)
        obs_file_row.addWidget(btn_obs)
        obs_main_layout.addLayout(obs_file_row)
        
        # Satellite selection for positioning
        sat_select_row = QHBoxLayout()
        self.pos_btn_select_all = QPushButton('Select All')
        self.pos_btn_select_all.setFixedHeight(24)
        self.pos_btn_deselect_all = QPushButton('Deselect All')
        self.pos_btn_deselect_all.setFixedHeight(24)
        self.pos_sat_count_label = QLabel('Load obs file to see satellites')
        self.pos_sat_count_label.setStyleSheet('color: #888; font-size: 10px;')
        sat_select_row.addWidget(self.pos_btn_select_all)
        sat_select_row.addWidget(self.pos_btn_deselect_all)
        sat_select_row.addWidget(self.pos_sat_count_label)
        sat_select_row.addStretch()
        obs_main_layout.addLayout(sat_select_row)
        
        # Satellite list with checkboxes
        from PySide6.QtWidgets import QListWidget
        self.pos_sat_list = QListWidget()
        self.pos_sat_list.setMaximumHeight(100)
        self.pos_sat_list.setSelectionMode(QListWidget.MultiSelection)
        self.pos_sat_list.setStyleSheet('''
            QListWidget { background: #1a1a2e; border: 1px solid #333; border-radius: 4px; }
            QListWidget::item { padding: 2px 4px; }
            QListWidget::item:selected { background: #3b82f6; color: white; }
        ''')
        obs_main_layout.addWidget(self.pos_sat_list)
        
        # Connect buttons
        self.pos_btn_select_all.clicked.connect(self._pos_select_all_sats)
        self.pos_btn_deselect_all.clicked.connect(self._pos_deselect_all_sats)
        
        obs_group.setStyleSheet(
            "QGroupBox{margin-top:6px;}"
            "QGroupBox::title{left:8px;padding:2px 4px;}"
        )
        scroll_layout.addWidget(obs_group)

        # Ephemeris source selection
        eph_group = QGroupBox('Ephemeris Source')
        eph_layout = QVBoxLayout(eph_group)
        eph_layout.setSpacing(6)
        
        src_row = QHBoxLayout()
        radio_style = "QRadioButton::indicator{width:18px;height:18px;} QRadioButton::indicator:checked{border:1px solid #3b82f6; background:#3b82f6; border-radius:9px;} QRadioButton::indicator:unchecked{border:1px solid #666; background:#111; border-radius:9px;}"
        self.use_nav_radio = QRadioButton('Use navigation file')
        self.use_nav_radio.setStyleSheet(radio_style)
        self.use_tle_radio = QRadioButton('Use TLE file')
        self.use_tle_radio.setStyleSheet(radio_style)
        self.use_tle_radio.setChecked(True)
        src_row.addWidget(self.use_nav_radio)
        src_row.addWidget(self.use_tle_radio)
        src_row.addStretch()
        eph_layout.addLayout(src_row)

        nav_row = QHBoxLayout()
        self.nav = QLineEdit()
        self.nav.setPlaceholderText('Select navigation file...')
        btn_nav = QPushButton('📂 Browse')
        btn_nav.clicked.connect(self.pick_nav)
        nav_row.addWidget(QLabel('Navigation file'))
        nav_row.addWidget(self.nav, 1)
        nav_row.addWidget(btn_nav)
        eph_layout.addLayout(nav_row)

        tle_row = QHBoxLayout()
        self.tle = QLineEdit()
        self.tle.setPlaceholderText('Select TLE file...')
        btn_tle = QPushButton('📂 Browse')
        btn_tle.clicked.connect(self.pick_tle)
        tle_row.addWidget(QLabel('TLE file'))
        tle_row.addWidget(self.tle, 1)
        tle_row.addWidget(btn_tle)
        eph_layout.addLayout(tle_row)
        
        scroll_layout.addWidget(eph_group)

        # Initial position
        pos_group = QGroupBox('Initial Position')
        pos_layout = QVBoxLayout(pos_group)
        pos_layout.setSpacing(6)
        
        mode_row = QHBoxLayout()
        self.mode_llh = QRadioButton('Lat/Lon/Alt')
        self.mode_llh.setStyleSheet("QRadioButton::indicator{width:18px;height:18px;} QRadioButton::indicator:checked{border:1px solid #3b82f6; background:#3b82f6; border-radius:9px;} QRadioButton::indicator:unchecked{border:1px solid #666; background:#111; border-radius:9px;}")
        self.mode_xyz = QRadioButton('ECEF XYZ')
        self.mode_xyz.setStyleSheet("QRadioButton::indicator{width:18px;height:18px;} QRadioButton::indicator:checked{border:1px solid #3b82f6; background:#3b82f6; border-radius:9px;} QRadioButton::indicator:unchecked{border:1px solid #666; background:#111; border-radius:9px;}")
        self.mode_llh.setChecked(True)
        mode_row.addWidget(self.mode_llh)
        mode_row.addWidget(self.mode_xyz)
        mode_row.addStretch()
        pos_layout.addLayout(mode_row)

        llh_row = QHBoxLayout()
        self.lat_deg = QLineEdit('39.80')
        self.lon_deg = QLineEdit('119.3')
        self.alt_m = QLineEdit('100')
        for lab, w, unit in [('Latitude', self.lat_deg, '°'), ('Longitude', self.lon_deg, '°'), ('Altitude', self.alt_m, 'm')]:
            llh_row.addWidget(QLabel(lab))
            llh_row.addWidget(w, 1)
            llh_row.addWidget(QLabel(unit))
        pos_layout.addLayout(llh_row)

        xyz_row = QHBoxLayout()
        self.x = QLineEdit('0')
        self.y = QLineEdit('0')
        self.z = QLineEdit('0')
        for lab, w in [('X', self.x), ('Y', self.y), ('Z', self.z)]:
            xyz_row.addWidget(QLabel(lab))
            xyz_row.addWidget(w, 1)
            xyz_row.addWidget(QLabel('m'))
        pos_layout.addLayout(xyz_row)
        
        scroll_layout.addWidget(pos_group)

        # True position (for error calculation)
        true_pos_group = QGroupBox('True Position (for error calculation)')
        true_pos_layout = QVBoxLayout(true_pos_group)
        true_pos_layout.setSpacing(6)
        
        true_mode_row = QHBoxLayout()
        radio_style = "QRadioButton::indicator{width:18px;height:18px;} QRadioButton::indicator:checked{border:1px solid #3b82f6; background:#3b82f6; border-radius:9px;} QRadioButton::indicator:unchecked{border:1px solid #666; background:#111; border-radius:9px;}"
        self.true_mode_llh = QRadioButton('Lat/Lon/Alt')
        self.true_mode_llh.setStyleSheet(radio_style)
        self.true_mode_xyz = QRadioButton('ECEF XYZ')
        self.true_mode_xyz.setStyleSheet(radio_style)
        self.true_mode_llh.setChecked(True)
        true_mode_row.addWidget(self.true_mode_llh)
        true_mode_row.addWidget(self.true_mode_xyz)
        true_mode_row.addStretch()
        true_pos_layout.addLayout(true_mode_row)

        true_llh_row = QHBoxLayout()
        self.true_lat_deg = QLineEdit('39.80')
        self.true_lon_deg = QLineEdit('119.29')
        self.true_alt_m = QLineEdit('60')
        for lab, w, unit in [('Latitude', self.true_lat_deg, '°'), ('Longitude', self.true_lon_deg, '°'), ('Altitude', self.true_alt_m, 'm')]:
            true_llh_row.addWidget(QLabel(lab))
            true_llh_row.addWidget(w, 1)
            true_llh_row.addWidget(QLabel(unit))
        true_pos_layout.addLayout(true_llh_row)

        true_xyz_row = QHBoxLayout()
        self.true_x = QLineEdit('0')
        self.true_y = QLineEdit('0')
        self.true_z = QLineEdit('0')
        for lab, w in [('X', self.true_x), ('Y', self.true_y), ('Z', self.true_z)]:
            true_xyz_row.addWidget(QLabel(lab))
            true_xyz_row.addWidget(w, 1)
            true_xyz_row.addWidget(QLabel('m'))
        true_pos_layout.addLayout(true_xyz_row)
        
        scroll_layout.addWidget(true_pos_group)

        # Error Models for positioning
        error_group = QGroupBox('Error Correction Models')
        error_layout = QVBoxLayout(error_group)
        error_layout.setSpacing(4)
        error_layout.setContentsMargins(8, 4, 8, 4)
        
        err_row1 = QHBoxLayout()
        self.pos_chk_iono = QCheckBox('Ionosphere')
        self.pos_chk_tropo = QCheckBox('Troposphere')
        self.pos_chk_iono.setChecked(True)
        self.pos_chk_tropo.setChecked(True)
        err_row1.addWidget(self.pos_chk_iono)
        err_row1.addWidget(self.pos_chk_tropo)
        err_row1.addStretch()
        error_layout.addLayout(err_row1)
        
        err_row2 = QHBoxLayout()
        self.pos_chk_multipath = QCheckBox('Multipath')
        self.pos_chk_relativistic = QCheckBox('Relativistic')
        err_row2.addWidget(self.pos_chk_multipath)
        err_row2.addWidget(self.pos_chk_relativistic)
        err_row2.addStretch()
        error_layout.addLayout(err_row2)
        
        err_row3 = QHBoxLayout()
        self.pos_chk_hardware = QCheckBox('Hardware Delay')
        self.pos_chk_clock = QCheckBox('Clock Bias/Drift')
        self.pos_chk_clock.setChecked(True)
        err_row3.addWidget(self.pos_chk_hardware)
        err_row3.addWidget(self.pos_chk_clock)
        err_row3.addStretch()
        error_layout.addLayout(err_row3)
        
        scroll_layout.addWidget(error_group)

        # Positioning Strategy
        strategy_group = QGroupBox('Positioning Strategy')
        strategy_layout = QVBoxLayout(strategy_group)
        strategy_layout.setSpacing(4)
        strategy_layout.setContentsMargins(8, 4, 8, 4)
        
        radio_style = "QRadioButton::indicator{width:16px;height:16px;}"
        
        self.strategy_pr = QRadioButton('Pseudorange Only')
        self.strategy_pr.setStyleSheet(radio_style + " QRadioButton::indicator:checked{border:1px solid #3b82f6; background:#3b82f6; border-radius:9px;} QRadioButton::indicator:unchecked{border:1px solid #666; background:#111; border-radius:9px;}")
        self.strategy_dop = QRadioButton('Doppler Only')
        self.strategy_dop.setStyleSheet(radio_style + " QRadioButton::indicator:checked{border:1px solid #3b82f6; background:#3b82f6; border-radius:9px;} QRadioButton::indicator:unchecked{border:1px solid #666; background:#111; border-radius:9px;}")
        self.strategy_dop.setChecked(True)
        self.strategy_fixed = QRadioButton('PR + Doppler (Fixed Weight)')
        self.strategy_fixed.setStyleSheet(radio_style + " QRadioButton::indicator:checked{border:1px solid #3b82f6; background:#3b82f6; border-radius:9px;} QRadioButton::indicator:unchecked{border:1px solid #666; background:#111; border-radius:9px;}")
        self.strategy_adaptive = QRadioButton('PR + Doppler (Adaptive Weight / VCE)')
        self.strategy_adaptive.setStyleSheet(radio_style + " QRadioButton::indicator:checked{border:1px solid #3b82f6; background:#3b82f6; border-radius:9px;} QRadioButton::indicator:unchecked{border:1px solid #666; background:#111; border-radius:9px;}")
        
        strategy_layout.addWidget(self.strategy_pr)
        strategy_layout.addWidget(self.strategy_dop)
        strategy_layout.addWidget(self.strategy_fixed)
        strategy_layout.addWidget(self.strategy_adaptive)
        
        scroll_layout.addWidget(strategy_group)

        # Run button
        self.runbtn = QPushButton('🎯 Start positioning')
        self.runbtn.setObjectName('primaryBtn')
        self.runbtn.clicked.connect(self.run_pos)
        scroll_layout.addWidget(self.runbtn)
        
        # Progress area
        progress_group = QGroupBox('Progress')
        progress_layout = QVBoxLayout(progress_group)
        progress_layout.setSpacing(6)
        
        self.progress_bar = QProgressBar()
        self.progress_bar.setMinimum(0)
        self.progress_bar.setMaximum(100)
        self.progress_bar.setValue(0)
        self.progress_bar.setTextVisible(True)
        self.progress_bar.setFormat('%p%')
        progress_layout.addWidget(self.progress_bar)
        
        self.status_label = QLabel('Ready')
        self.status_label.setStyleSheet('color: #666666; font-size: 11px;')
        progress_layout.addWidget(self.status_label)
        
        scroll_layout.addWidget(progress_group)
        
        scroll.setWidget(content)
        layout.addWidget(scroll)
        
        # Right chart area (2 rows):
        # Row 1: Sky plot | Residuals
        # Row 2: Positioning error / convergence (spans 2 columns)
        right_container = QWidget()
        right_panel = QGridLayout(right_container)
        right_panel.setSpacing(10)
        right_panel.setContentsMargins(0, 0, 0, 0)
        
        # MacOS-style light gray background
        chart_bg = '#f5f5f7'
        axis_bg = '#ffffff'
        
        # Sky plot: satellite azimuth/elevation relative to station
        self.fig_sky = Figure(figsize=(6, 4), facecolor=chart_bg)
        self.ax_sky = self.fig_sky.add_subplot(111, projection='polar')
        self.canvas_sky = FigureCanvasQTAgg(self.fig_sky)
        self.canvas_sky.setMinimumHeight(280)

        # Residuals: pseudorange + doppler per epoch (single plot with dual Y axes)
        self.fig_res = Figure(figsize=(6, 4), facecolor=chart_bg)
        self.ax_res = self.fig_res.add_subplot(111)
        self.ax_res2 = self.ax_res.twinx()
        self.canvas_res = FigureCanvasQTAgg(self.fig_res)
        self.canvas_res.setMinimumHeight(280)

        # Convergence / positioning error
        self.fig_conv = Figure(figsize=(6, 5), facecolor=chart_bg)
        self.ax_conv = self.fig_conv.add_subplot(111)
        self.ax_conv.set_facecolor(axis_bg)
        self.canvas_conv = FigureCanvasQTAgg(self.fig_conv)
        self.canvas_conv.setMinimumHeight(320)

        # Size policy to fill grid cells nicely
        for c in (self.canvas_sky, self.canvas_res, self.canvas_conv):
            c.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        # Place widgets
        right_panel.addWidget(self.canvas_sky, 0, 0)
        right_panel.addWidget(self.canvas_res, 0, 1)
        right_panel.addWidget(self.canvas_conv, 1, 0, 1, 2)
        right_panel.setColumnStretch(0, 1)
        right_panel.setColumnStretch(1, 1)
        right_panel.setRowStretch(0, 1)
        right_panel.setRowStretch(1, 2)
        
        # Initialize chart style
        self._style_chart()
        self._style_skyplot()
        self._style_residuals()
        
        main_layout.addWidget(left_panel, 1)
        main_layout.addWidget(right_container, 2)

        # State toggle logic
        def _toggle_src():
            self.nav.setEnabled(self.use_nav_radio.isChecked())
            btn_nav.setEnabled(self.use_nav_radio.isChecked())
            self.tle.setEnabled(self.use_tle_radio.isChecked())
            btn_tle.setEnabled(self.use_tle_radio.isChecked())
        def _toggle_mode():
            llh = self.mode_llh.isChecked()
            self.lat_deg.setEnabled(llh)
            self.lon_deg.setEnabled(llh)
            self.alt_m.setEnabled(llh)
            self.x.setEnabled(not llh)
            self.y.setEnabled(not llh)
            self.z.setEnabled(not llh)
        def _toggle_true_mode():
            llh = self.true_mode_llh.isChecked()
            self.true_lat_deg.setEnabled(llh)
            self.true_lon_deg.setEnabled(llh)
            self.true_alt_m.setEnabled(llh)
            self.true_x.setEnabled(not llh)
            self.true_y.setEnabled(not llh)
            self.true_z.setEnabled(not llh)
        self.use_nav_radio.toggled.connect(_toggle_src)
        self.use_tle_radio.toggled.connect(_toggle_src)
        self.mode_llh.toggled.connect(_toggle_mode)
        self.mode_xyz.toggled.connect(_toggle_mode)
        self.true_mode_llh.toggled.connect(_toggle_true_mode)
        self.true_mode_xyz.toggled.connect(_toggle_true_mode)
        
        # Live conversion for true position
        self.true_lat_deg.textChanged.connect(self._on_true_llh_changed)
        self.true_lon_deg.textChanged.connect(self._on_true_llh_changed)
        self.true_alt_m.textChanged.connect(self._on_true_llh_changed)
        self.true_x.textChanged.connect(self._on_true_xyz_changed)
        self.true_y.textChanged.connect(self._on_true_xyz_changed)
        self.true_z.textChanged.connect(self._on_true_xyz_changed)
        
        _toggle_src()
        _toggle_mode()
        _toggle_true_mode()
        self._on_true_llh_changed()  # Initialize XYZ display

    def _on_true_llh_changed(self):
        """Update XYZ when true lat/lon/alt changes"""
        if not self.true_mode_llh.isChecked():
            return
        
        try:
            lat = float(self.true_lat_deg.text())
            lon = float(self.true_lon_deg.text())
            alt = float(self.true_alt_m.text())
            
            # Convert to radians
            lat_rad = np.deg2rad(lat)
            lon_rad = np.deg2rad(lon)
            
            # Convert to ECEF
            x, y, z = geodetic_to_ecef(lat_rad, lon_rad, alt)
            
            # Update XYZ display (block reverse trigger)
            self.true_x.blockSignals(True)
            self.true_y.blockSignals(True)
            self.true_z.blockSignals(True)
            self.true_x.setText(f'{x:.2f}')
            self.true_y.setText(f'{y:.2f}')
            self.true_z.setText(f'{z:.2f}')
            self.true_x.blockSignals(False)
            self.true_y.blockSignals(False)
            self.true_z.blockSignals(False)
        except Exception:
            pass
    
    def _on_true_xyz_changed(self):
        """Update true lat/lon/alt when XYZ changes"""
        if not self.true_mode_xyz.isChecked():
            return
        
        try:
            x = float(self.true_x.text())
            y = float(self.true_y.text())
            z = float(self.true_z.text())
            
            # Import conversion helper
            from core.utilities import ecef_to_geodetic
            
            # Convert to geodetic
            lat_rad, lon_rad, alt = ecef_to_geodetic(x, y, z)
            
            # Convert to degrees
            lat = np.rad2deg(lat_rad)
            lon = np.rad2deg(lon_rad)
            
            # Update LLH display (block reverse trigger)
            self.true_lat_deg.blockSignals(True)
            self.true_lon_deg.blockSignals(True)
            self.true_alt_m.blockSignals(True)
            self.true_lat_deg.setText(f'{lat:.6f}')
            self.true_lon_deg.setText(f'{lon:.6f}')
            self.true_alt_m.setText(f'{alt:.2f}')
            self.true_lat_deg.blockSignals(False)
            self.true_lon_deg.blockSignals(False)
            self.true_alt_m.blockSignals(False)
        except Exception:
            pass
    
    def _get_true_position(self):
        """Get true position ECEF coordinates"""
        try:
            if self.true_mode_llh.isChecked():
                lat = float(self.true_lat_deg.text())
                lon = float(self.true_lon_deg.text())
                alt = float(self.true_alt_m.text())
                return list(geodetic_to_ecef(np.deg2rad(lat), np.deg2rad(lon), alt))
            else:
                return [float(self.true_x.text()), float(self.true_y.text()), float(self.true_z.text())]
        except Exception:
            return None

    def _style_chart(self):
        """Initialize convergence chart style"""
        ax = self.ax_conv
        ax.set_facecolor('#ffffff')
        ax.set_title('Positioning Convergence & Error', fontsize=11, color='#1f2937', fontweight='bold', pad=12)
        ax.set_xlabel('Iteration Step', fontsize=9, color='#4b5563')
        ax.set_ylabel('Position Std (m)', fontsize=9, color='#2563eb')
        ax.tick_params(colors='#4b5563', labelsize=8)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_color('#e5e7eb')
        ax.spines['left'].set_color('#e5e7eb')
        ax.grid(True, linestyle=':', alpha=0.4, color='#9ca3af')
        self.fig_conv.tight_layout()
        self.canvas_conv.draw()

    def _style_skyplot(self):
        """Initialize sky plot style"""
        ax = self.ax_sky
        ax.clear()
        ax.set_facecolor('#ffffff')
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        ax.set_ylim(0, 90)
        ax.set_yticks([0, 30, 60, 90])
        ax.set_yticklabels(['90°', '60°', '30°', '0°'], fontsize=8, color='#6b7280')
        ax.set_xticks(np.deg2rad([0, 45, 90, 135, 180, 225, 270, 315]))
        ax.set_xticklabels(['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'], fontsize=8, color='#6b7280')
        ax.set_title('Satellite Sky Plot (Az/El)', fontsize=11, color='#1f2937', fontweight='bold', pad=12)
        ax.grid(True, linestyle=':', alpha=0.4, color='#9ca3af')
        
        # Remove polar spine
        ax.spines['polar'].set_visible(False)
        
        self.fig_sky.tight_layout()
        self.canvas_sky.draw()

    def _style_residuals(self):
        """Initialize residual chart style"""
        ax = self.ax_res
        ax2 = self.ax_res2

        ax.clear()
        ax2.clear()

        ax.set_facecolor('#ffffff')
        ax.grid(True, linestyle=':', alpha=0.4, color='#9ca3af')
        
        # Left Y axis (Pseudorange) - blue ticks and label
        ax.tick_params(axis='x', colors='#4b5563', labelsize=8)
        ax.tick_params(axis='y', colors='#2563eb', labelsize=8)
        ax.yaxis.set_label_position('left')
        ax.yaxis.tick_left()
        ax.set_ylabel('Pseudorange (m)', fontsize=9, color='#2563eb')
        
        # Show left and bottom spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(True)
        ax.spines['bottom'].set_color('#e5e7eb')
        ax.spines['left'].set_visible(True)
        ax.spines['left'].set_color('#2563eb')

        # Right Y axis (Doppler) - red ticks and label
        ax2.tick_params(axis='y', colors='#dc2626', labelsize=8)
        ax2.yaxis.set_label_position('right')
        ax2.yaxis.tick_right()
        ax2.set_ylabel('Doppler (Hz)', fontsize=9, color='#dc2626')
        
        # Show right spine for ax2
        ax2.spines['top'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax2.spines['bottom'].set_visible(False)
        ax2.spines['right'].set_visible(True)
        ax2.spines['right'].set_color('#dc2626')

        ax.set_title('Residuals per Epoch', fontsize=11, color='#1f2937', fontweight='bold', pad=12)
        ax.set_xlabel('Time (s)', fontsize=9, color='#4b5563')
        
        self.fig_res.tight_layout()
        self.canvas_res.draw()

    def _plot_skyplot(self, results):
        """Render skyplot on the right panel (best-effort)."""
        sky = results.get('skyplot') if isinstance(results, dict) else None
        
        # Reset style
        self._style_skyplot()
        ax = self.ax_sky

        if not sky:
            ax.text(0.5, 0.5, 'No skyplot data', transform=ax.transAxes,
                    ha='center', va='center', fontsize=9, color='#9ca3af')
            self.canvas_sky.draw()
            return

        try:
            t = np.asarray(sky.get('times_s', []), dtype=float)
            az = np.asarray(sky.get('azimuth_deg', []), dtype=float)
            el = np.asarray(sky.get('elevation_deg', []), dtype=float)
            prn = np.asarray(sky.get('prn', []), dtype=int) if sky.get('prn', None) is not None else None
            
            if len(az) == 0 or len(el) == 0:
                raise ValueError("empty skyplot")

            az_rad = np.deg2rad(az)
            r = 90.0 - el  # center=zenith, edge=horizon
            
            # Only keep within plot range
            m = np.isfinite(az_rad) & np.isfinite(r) & (r >= 0) & (r <= 90)
            az_rad = az_rad[m]
            r = r[m]
            t = t[m] if len(t) else t
            if prn is not None and len(prn) == len(m):
                prn = prn[m]

            if len(az_rad) == 0:
                raise ValueError("no valid skyplot points")

            # Prefer grouping by PRN if available; otherwise color by time
            if prn is not None and len(prn) == len(r) and len(prn) > 0:
                uniq = np.unique(prn)
                # Use a pleasant qualitative colormap
                cmap = plt.get_cmap('Set2')
                for i, p in enumerate(uniq[:12]):  # cap legend clutter
                    idx = prn == p
                    # Cycle through colors
                    color = cmap(i % 8)
                    ax.plot(az_rad[idx], r[idx], '-', linewidth=1.5, alpha=0.6, color=color)
                    ax.scatter(az_rad[idx], r[idx], s=16, alpha=0.9, color=color, edgecolors='white', linewidth=0.5, label=f'{int(p):02d}')
                
                if len(uniq) <= 12:
                    leg = ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1), fontsize=8, frameon=True, framealpha=0.9)
                    leg.get_frame().set_linewidth(0.0)
            else:
                c = t if len(t) == len(r) else None
                sc = ax.scatter(az_rad, r, c=c, s=16, cmap='viridis', alpha=0.8, edgecolors='none')
                if c is not None and len(t) > 0:
                    cbar = self.fig_sky.colorbar(sc, ax=ax, pad=0.1, shrink=0.7)
                    cbar.set_label('Time (s)', fontsize=8, color='#4b5563')
                    cbar.ax.tick_params(labelsize=8, colors='#4b5563')
                    cbar.outline.set_visible(False)

            # Start/end markers with modern colors
            ax.plot(az_rad[0], r[0], marker='o', markersize=6, color='#10b981', markeredgecolor='white', markeredgewidth=1.5, zorder=10)
            ax.plot(az_rad[-1], r[-1], marker='s', markersize=6, color='#ef4444', markeredgecolor='white', markeredgewidth=1.5, zorder=10)
            
        except Exception:
            ax.text(0.5, 0.5, 'Failed to render skyplot', transform=ax.transAxes,
                    ha='center', va='center', fontsize=9, color='#ef4444')

        self.fig_sky.tight_layout()
        self.canvas_sky.draw()

    def _plot_residuals(self, results):
        """Render pseudorange/doppler residuals on the right panel (best-effort)."""
        res = results.get('residuals') if isinstance(results, dict) else None
        self._style_residuals()
        ax = self.ax_res
        ax2 = self.ax_res2

        if not res:
            ax.text(0.5, 0.5, 'No residual data', transform=ax.transAxes,
                    ha='center', va='center', fontsize=9, color='#9ca3af')
            self.canvas_res.draw()
            return

        try:
            t = np.asarray(res.get('times_s', []), dtype=float)
            pr = res.get('pseudorange_m', None)
            dop = res.get('doppler_hz', None)

            # Zero lines
            ax.axhline(0, color='#9ca3af', linestyle='--', linewidth=1, alpha=0.5)

            rms_lines = []

            if pr is not None:
                pr = np.asarray(pr, dtype=float)
                if len(t) == len(pr) and len(t) > 0:
                    ax.plot(t, pr, 'o', markersize=3, color='#2563eb', alpha=0.6, markeredgewidth=0, label='Pseudorange')
                    rms = float(np.sqrt(np.nanmean(pr * pr)))
                    rms_lines.append(f"PR RMS: {rms:.3f} m")
                else:
                    ax.text(0.02, 0.90, 'Pseudorange unavailable', transform=ax.transAxes, fontsize=8, color='#9ca3af')
            else:
                ax.text(0.02, 0.90, 'Pseudorange unused', transform=ax.transAxes, fontsize=8, color='#9ca3af')

            if dop is not None:
                dop = np.asarray(dop, dtype=float)
                if len(t) == len(dop) and len(t) > 0:
                    ax2.plot(t, dop, 'o', markersize=3, color='#dc2626', alpha=0.6, markeredgewidth=0, label='Doppler')
                    rms = float(np.sqrt(np.nanmean(dop * dop)))
                    rms_lines.append(f"Dop RMS: {rms:.3f} Hz")
                else:
                    ax.text(0.02, 0.82, 'Doppler unavailable', transform=ax.transAxes, fontsize=8, color='#9ca3af')
            else:
                ax.text(0.02, 0.82, 'Doppler unused', transform=ax.transAxes, fontsize=8, color='#9ca3af')

            # Combined legend
            h1, l1 = ax.get_legend_handles_labels()
            h2, l2 = ax2.get_legend_handles_labels()
            if h1 or h2:
                leg = ax.legend(h1 + h2, l1 + l2, loc='upper right', fontsize=8, frameon=True, framealpha=0.95, facecolor='white')
                leg.get_frame().set_linewidth(0.0)

            # RMS summary box
            if rms_lines:
                ax.text(0.02, 0.02, "\n".join(rms_lines), transform=ax.transAxes,
                        ha='left', va='bottom', fontsize=8, color='#374151',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='#f3f4f6', alpha=0.9, edgecolor='none'))
        except Exception:
            ax.text(0.5, 0.5, 'Failed to render residuals', transform=ax.transAxes,
                    ha='center', va='center', fontsize=9, color='#ef4444')

        self.fig_res.tight_layout()
        self.canvas_res.draw()

    def pick_obs(self):
        fn, _ = QFileDialog.getOpenFileName(self, 'Select observation CSV/OBS', '', 'Observation Files (*.csv *.obs *.txt);;All Files (*.*)')
        if fn:
            self.obs.setText(fn)
            self._load_obs_satellite_list(fn)
    
    def _load_obs_satellite_list(self, obs_path):
        """Parse observation file and populate satellite selection list."""
        from PySide6.QtWidgets import QListWidgetItem
        self.pos_sat_list.clear()
        
        try:
            prns = set()
            with open(obs_path, 'r', encoding='utf-8') as f:
                # Check file format
                first_line = ''
                for line in f:
                    if line.strip():
                        first_line = line.strip()
                        break
                f.seek(0)
                
                if first_line.startswith('OBS:') or first_line.startswith('BS:'):
                    # Custom OBS format
                    for line in f:
                        line = line.strip()
                        if not line or not (line.startswith('OBS:') or line.startswith('BS:')):
                            continue
                        fields = line.split()
                        if len(fields) >= 12:
                            try:
                                prn = fields[11]  # PRN is in field 11
                                prns.add(prn)
                            except:
                                pass
                else:
                    # CSV format
                    import csv
                    f.seek(0)
                    reader = csv.DictReader(f)
                    for row in reader:
                        prn = row.get('prn', row.get('PRN', row.get('sat_prn', '')))
                        if prn:
                            prns.add(str(prn))
            
            # Sort PRNs and populate list
            sorted_prns = sorted(prns, key=lambda x: (not x.isdigit(), int(x) if x.isdigit() else x))
            for prn in sorted_prns:
                item = QListWidgetItem(f'PRN {prn}')
                item.setData(Qt.UserRole, prn)
                item.setSelected(True)  # Select by default
                self.pos_sat_list.addItem(item)
            
            self.pos_sat_count_label.setText(f'{len(sorted_prns)} satellites found')
            
        except Exception as e:
            self.pos_sat_count_label.setText(f'Error: {str(e)[:25]}')
    
    def _pos_select_all_sats(self):
        """Select all satellites in the positioning list."""
        for i in range(self.pos_sat_list.count()):
            self.pos_sat_list.item(i).setSelected(True)
    
    def _pos_deselect_all_sats(self):
        """Deselect all satellites in the positioning list."""
        for i in range(self.pos_sat_list.count()):
            self.pos_sat_list.item(i).setSelected(False)
    
    def _get_pos_selected_prns(self):
        """Get list of selected satellite PRNs for positioning."""
        selected = []
        for i in range(self.pos_sat_list.count()):
            item = self.pos_sat_list.item(i)
            if item.isSelected():
                prn = item.data(Qt.UserRole)
                selected.append(prn)
        return selected if selected else None  # None means all satellites

    def pick_nav(self):
        fn, _ = QFileDialog.getOpenFileName(self, 'Select navigation file', '', 'Navigation Files (*.nav *.rnx *.txt);;All Files (*.*)')
        if fn:
            self.nav.setText(fn)

    def pick_tle(self):
        fn, _ = QFileDialog.getOpenFileName(self, 'Select TLE file', '', 'TLE Files (*.tle *.txt);;All Files (*.*)')
        if fn:
            self.tle.setText(fn)

    def run_pos(self):
        # If already running, cancel
        if self.thread is not None and self.thread.isRunning():
            if self.worker:
                self.worker.cancel()
            self.thread.quit()
            self.thread.wait()
            self._reset_ui()
            return
        
        if self.mode_llh.isChecked():
            try:
                lat = float(self.lat_deg.text())
                lon = float(self.lon_deg.text())
                alt = float(self.alt_m.text())
                import numpy as np
                init_pos = list(geodetic_to_ecef(np.deg2rad(lat), np.deg2rad(lon), alt))
            except Exception:
                init_pos = [0.0,0.0,0.0]
        else:
            init_pos = [float(self.x.text()), float(self.y.text()), float(self.z.text())]
        eph = self.nav.text() if self.use_nav_radio.isChecked() and self.nav.text() else None
        tle = self.tle.text() if self.use_tle_radio.isChecked() and self.tle.text() else None
        
        # Determine positioning strategy
        if self.strategy_pr.isChecked():
            strategy = 'pr'
        elif self.strategy_dop.isChecked():
            strategy = 'dop'
        elif self.strategy_fixed.isChecked():
            strategy = 'fixed'
        else:
            strategy = 'adaptive'
        
        # Get error correction settings
        atmcor = self.pos_chk_iono.isChecked() or self.pos_chk_tropo.isChecked()
        
        # Init progress
        self.progress_bar.setValue(0)
        self.status_label.setText('Initializing...')
        self.runbtn.setText('⏹ Cancel')
        self.runbtn.setEnabled(True)
        
        # Get selected satellites for positioning
        selected_prns = self._get_pos_selected_prns()
        
        enable_iono = self.pos_chk_iono.isChecked()
        enable_tropo = self.pos_chk_tropo.isChecked()
        enable_multipath = self.pos_chk_multipath.isChecked()
        enable_relativistic = self.pos_chk_relativistic.isChecked()
        enable_hardware = self.pos_chk_hardware.isChecked()

        # Create worker thread
        self.thread = QThread()
        self.worker = PositioningWorker(
            self.obs.text(), eph, tle, init_pos, strategy, enable_iono,enable_tropo,enable_multipath,enable_relativistic,enable_hardware,
            selected_prns=selected_prns
        )

        self.worker.moveToThread(self.thread)
        
        # Wire signals
        self.thread.started.connect(self.worker.run)
        self.worker.progress.connect(self._on_progress)
        self.worker.finished.connect(self._on_finished)
        self.worker.error.connect(self._on_error)
        self.worker.finished.connect(self.thread.quit)
        self.worker.error.connect(self.thread.quit)
        self.thread.finished.connect(self._cleanup_thread)
        
        # Start thread
        self.thread.start()
    
    def _on_progress(self, percent, message):
        """Update progress"""
        self.progress_bar.setValue(percent)
        self.status_label.setText(message)
    
    def _on_finished(self, convergence_data, results):
        """Positioning finished"""
        self._reset_ui()
        self.status_label.setText('✅ Positioning completed!')
        self.status_label.setStyleSheet('color: #22c55e; font-size: 11px; font-weight: bold;')
        
        # Get true position
        true_pos = self._get_true_position()
        
        # Plot convergence
        if convergence_data and len(convergence_data.get('obs_count', [])) > 0:
            self._plot_convergence(convergence_data, results, true_pos)

        # Plot skyplot + residuals
        if results:
            self._plot_skyplot(results)
            self._plot_residuals(results)
    
    def _on_error(self, error_msg):
        """Handle error"""
        self._reset_ui()
        self.status_label.setText('❌ Failed')
        self.status_label.setStyleSheet('color: #ef4444; font-size: 11px; font-weight: bold;')
        print(error_msg)
    
    def _reset_ui(self):
        """Reset UI state"""
        self.runbtn.setText('🎯 Start positioning')
        self.runbtn.setEnabled(True)
    
    def _cleanup_thread(self):
        """Cleanup thread"""
        self.worker = None
        self.thread = None
    
    def _plot_convergence(self, convergence_data, results, true_pos=None):
        """Plot positioning convergence"""
        ax = self.ax_conv
        ax.clear()
        ax.set_facecolor('#ffffff')
        
        obs_count = convergence_data.get('obs_count', [])
        position_std = convergence_data.get('position_std', [])
        residual_rms = convergence_data.get('residual_rms', [])
        position_history = convergence_data.get('position_history', [])  # Position history
        
        if len(obs_count) > 0 and len(position_std) > 0:
            # Plot position STD convergence
            color1 = '#3b82f6'
            ax.plot(obs_count, position_std, 'o-', color=color1, 
                    linewidth=2, markersize=4, label='Position Std (m)')
            ax.set_ylim([0, 1000])
            ax.set_xlabel('Observation Count', fontsize=11, color='#333333')
            ax.set_ylabel('Position Std (m)', fontsize=11, color=color1)
            ax.tick_params(axis='y', labelcolor=color1)
            
            # Secondary Y-axis for absolute position error
            ax2 = ax.twinx()
            color2 = '#ef4444'
            
            # Compute absolute position error
            absolute_errors = []
            if true_pos is not None and len(position_history) > 0:
                true_pos_arr = np.array(true_pos)
                for pos in position_history:
                    if pos is not None and len(pos) >= 3:
                        pos_arr = np.array(pos[:3])
                        error = np.linalg.norm(pos_arr - true_pos_arr)
                        absolute_errors.append(error)
                    else:
                        absolute_errors.append(np.nan)
            
            if len(absolute_errors) > 0 and len(absolute_errors) == len(obs_count):
                ax2.plot(obs_count, absolute_errors, 's--', color=color2,
                         linewidth=1.5, markersize=3, alpha=0.7, label='Positioning Error (m)')
                ax2.set_ylabel('Positioning Error (m)', fontsize=11, color=color2)
                
                # Set Y-axis range
                valid_errors = [e for e in absolute_errors if not np.isnan(e)]
                if valid_errors:
                    max_error = max(valid_errors)
                    ax2.set_ylim([0, min(max_error * 1.2, 5000)])
            elif len(residual_rms) > 0:
                # If no true position, show residual RMS
                ax2.plot(obs_count, residual_rms, 's--', color=color2,
                         linewidth=1.5, markersize=3, alpha=0.7, label='Residual RMS')
                ax2.set_ylabel('Residual RMS', fontsize=11, color=color2)
            
            ax2.tick_params(axis='y', labelcolor=color2)
            ax2.spines['right'].set_color(color2)
            
            # Style settings
            ax.set_title('Positioning Convergence', fontsize=12, color='#333333', fontweight='bold', pad=10)
            ax.tick_params(colors='#333333', labelsize=9)
            ax.spines['top'].set_visible(False)
            ax.spines['bottom'].set_color('#cccccc')
            ax.spines['left'].set_color(color1)
            ax.grid(True, linestyle='--', alpha=0.3, color='#999999')
            
            # Add final result text
            if results:
                final_text = f"Final Result:\n"
                final_text += f"Lat: {results.get('latitude_deg', 0):.6f}°\n"
                final_text += f"Lon: {results.get('longitude_deg', 0):.6f}°\n" 
                final_text += f"Alt: {results.get('height_m', 0):.2f} m"
                
                # Add absolute error info
                if len(absolute_errors) > 0:
                    valid_errors = [e for e in absolute_errors if not np.isnan(e)]
                    if valid_errors:
                        final_error = valid_errors[-1]
                        final_text += f"\nAbs Err: {final_error:.2f} m"
                
                ax.text(0.02, 0.98, final_text, transform=ax.transAxes, fontsize=9,
                        verticalalignment='top', color='#333333',
                        bbox=dict(boxstyle='round', facecolor='#f0f4f8', alpha=0.8, edgecolor='#cccccc'))
            
            # Add legend
            lines1, labels1 = ax.get_legend_handles_labels()
            lines2, labels2 = ax2.get_legend_handles_labels()
            ax.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=8)
        
        self.fig_conv.tight_layout()
        self.canvas_conv.draw()

class SettingsTab(QWidget):
    def __init__(self, app):
        super().__init__()
        self.app = app
        self.s = load_settings()
        
        layout = QVBoxLayout(self)
        layout.setSpacing(16)
        layout.setContentsMargins(16, 16, 16, 16)
        
        # Appearance
        appearance_group = QGroupBox('Appearance')
        appearance_layout = QVBoxLayout(appearance_group)
        appearance_layout.setSpacing(12)
        
        theme_row = QHBoxLayout()
        theme_label = QLabel('UI Theme')
        theme_label.setMinimumWidth(120)
        self.theme = QComboBox()
        self.theme.addItems(['light', 'dark'])
        self.theme.setCurrentText(self.s.get('theme', 'light'))
        theme_row.addWidget(theme_label)
        theme_row.addWidget(self.theme, 1)
        appearance_layout.addLayout(theme_row)
        
        accent_row = QHBoxLayout()
        accent_label = QLabel('Accent Color')
        accent_label.setMinimumWidth(120)
        self.accent = QLineEdit(self.s.get('accent', '#3b82f6'))
        self.accent.setPlaceholderText('#3b82f6')
        accent_row.addWidget(accent_label)
        accent_row.addWidget(self.accent, 1)
        appearance_layout.addLayout(accent_row)
        
        layout.addWidget(appearance_group)
        
        # Cesium configuration
        cesium_group = QGroupBox('Cesium Globe Config')
        cesium_layout = QVBoxLayout(cesium_group)
        cesium_layout.setSpacing(12)
        
        token_row = QHBoxLayout()
        token_label = QLabel('Ion Token')
        token_label.setMinimumWidth(120)
        self.ion = QLineEdit(self.s.get('ion_token', ''))
        self.ion.setPlaceholderText('Enter Cesium Ion Token...')
        self.ion.setEchoMode(QLineEdit.Password)
        token_row.addWidget(token_label)
        token_row.addWidget(self.ion, 1)
        cesium_layout.addLayout(token_row)
        
        terrain_row = QHBoxLayout()
        terrain_label = QLabel('Enable high-precision terrain')
        terrain_label.setMinimumWidth(120)
        self.use_terrain = QCheckBox('Use world terrain and high-res imagery')
        self.use_terrain.setChecked(bool(self.s.get('use_world_terrain', False)))
        terrain_row.addWidget(terrain_label)
        terrain_row.addWidget(self.use_terrain)
        terrain_row.addStretch()
        cesium_layout.addLayout(terrain_row)
        
        layout.addWidget(cesium_group)
        
        # Save button
        btn = QPushButton('💾 Save & Apply')
        btn.setObjectName('primaryBtn')
        btn.clicked.connect(self.save_apply)
        layout.addWidget(btn)
        
        # About button
        about_btn = QPushButton('ℹ️ About SPINE')
        about_btn.clicked.connect(self.show_about)
        layout.addWidget(about_btn)
        
        layout.addStretch()

    def save_apply(self):
        self.s['theme'] = self.theme.currentText()
        self.s['accent'] = self.accent.text()
        self.s['ion_token'] = self.ion.text().strip()
        self.s['use_world_terrain'] = bool(self.use_terrain.isChecked())
        save_settings(self.s)
        apply_qss(self.app, self.s['theme'], self.s['accent'])
    
    def show_about(self):
        """Show About dialog."""
        dialog = AboutDialog(self)
        dialog.exec()

class AboutDialog(QDialog):
    """About dialog showing application information and logo."""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('About SPINE')
        self.setMinimumWidth(400)
        self.setMinimumHeight(300)
        
        layout = QVBoxLayout(self)
        layout.setSpacing(16)
        layout.setContentsMargins(24, 24, 24, 24)
        
        # Logo
        logo_path = get_logo_path()
        if logo_path and os.path.exists(logo_path):
            logo_label = QLabel()
            pixmap = QPixmap(logo_path)
            # Scale logo to reasonable size (max 200px width)
            if pixmap.width() > 200:
                pixmap = pixmap.scaled(200, 200, Qt.KeepAspectRatio, Qt.SmoothTransformation)
            logo_label.setPixmap(pixmap)
            logo_label.setAlignment(Qt.AlignCenter)
            layout.addWidget(logo_label)
        
        # Application name and version
        title_label = QLabel('SPINE')
        title_font = QFont()
        title_font.setPointSize(18)
        title_font.setBold(True)
        title_label.setFont(title_font)
        title_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(title_label)
        
        subtitle_label = QLabel('Satellite Positioning, Initialization & Navigation Estimation')
        subtitle_label.setAlignment(Qt.AlignCenter)
        subtitle_label.setStyleSheet('color: #666;')
        layout.addWidget(subtitle_label)
        
        # Version and info
        info_text = """
        <p style="text-align: center;">
        <b>Version:</b> 1.0.0<br>
        <b>License:</b> MIT<br><br>
        LEO-SPINE is an open-source desktop application for LEO satellite 
        constellation design, observation simulation, and positioning computation.
        </p>
        """
        info_label = QLabel(info_text)
        info_label.setWordWrap(True)
        layout.addWidget(info_label)
        
        # Buttons
        buttons = QDialogButtonBox(QDialogButtonBox.Ok)
        buttons.accepted.connect(self.accept)
        layout.addWidget(buttons)

def main():
    app = QApplication(sys.argv)
    s = load_settings()
    apply_qss(app, s.get('theme', 'light'), s.get('accent', '#3b82f6'))
    
    # Set global font
    font = QFont('Microsoft YaHei UI', 10)
    font.setStyleStrategy(QFont.PreferAntialias)
    app.setFont(font)
    
    w = QMainWindow()
    w.setWindowTitle('SPINE - Satellite Constellation Simulation & Positioning')
    
    # Set window icon
    logo_path = get_logo_path()
    if logo_path and os.path.exists(logo_path):
        icon = QIcon(logo_path)
        w.setWindowIcon(icon)
        app.setWindowIcon(icon)
    
    # Main tabs
    tabs = QTabWidget()
    tabs.setDocumentMode(True)
    tabs.addTab(ConstellationTab(), '🛰️ Constellation Design')
    tabs.addTab(SimulationTab(), '📡 Observation Simulation')
    tabs.addTab(PositioningTab(), '📍 Positioning')
    tabs.addTab(SettingsTab(app), '⚙️ Settings')
    
    w.setCentralWidget(tabs)
    w.resize(1400, 900)
    w.show()
    sys.exit(app.exec())
