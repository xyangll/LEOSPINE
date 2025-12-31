import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
from datetime import datetime, timezone
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from core.constellation import generate_walker_delta_tles, write_tle_file
import webbrowser
import os
from datetime import timedelta
from app.czml import write_czml
from app.settings import write_web_theme, write_web_config, load_settings, hex_to_rgba
from pathlib import Path
from app.webserver import start as start_web
from core.sim_data import SatelliteSimulation, calculate_elevation_angle, teme_to_ecef, gps_time_to_datetime, datetime_to_jday
from tools.sim_obs import run as run_sim_obs
from tools.run_positioning import run_lsq_positioning

class ConstellationPage(ttk.Frame):
    def __init__(self, master):
        super().__init__(master)
        self.sim = SatelliteSimulation()
        self.tle_file = None
        self._build_ui()

    def _build_ui(self):
        style = ttk.Style()
        try:
            style.theme_use('clam')
        except Exception:
            pass
        style.configure('Card.TFrame', background='#12161c')
        style.configure('Card.TLabel', background='#12161c', foreground='#e6e6e6')
        style.configure('Card.TButton', padding=6)
        form = ttk.Frame(self, style='Card.TFrame')
        form.grid(row=0, column=0, sticky="nsew")
        self.columnconfigure(0, weight=0)
        self.columnconfigure(1, weight=1)
        labels = ["Semi-major axis (m)","Eccentricity","Inclination (deg)","RAAN (deg)","Argument of perigee (deg)","Mean anomaly (deg)","Planes P","Satellites per plane S","Phasing F","Epoch YYYY-MM-DD"]
        defaults = [26560000,0.01,55.0,0.0,0.0,0.0,3,4,1,datetime.utcnow().strftime("%Y-%m-%d")]
        self.entries = []
        for i,(lab,val) in enumerate(zip(labels,defaults)):
            ttk.Label(form,text=lab,style='Card.TLabel').grid(row=i,column=0,sticky="w",padx=10,pady=6)
            e = ttk.Entry(form)
            e.insert(0,str(val))
            e.grid(row=i,column=1,sticky="ew",padx=10,pady=6)
            self.entries.append(e)
        btns = ttk.Frame(form, style='Card.TFrame')
        btns.grid(row=len(labels),column=0,columnspan=2,sticky="ew",padx=10,pady=8)
        ttk.Button(btns,text="Generate TLE",command=self.on_generate_tle,style='Card.TButton').pack(side=tk.LEFT,padx=4)
        ttk.Button(btns,text="Save TLE",command=self.on_save_tle,style='Card.TButton').pack(side=tk.LEFT,padx=4)
        ttk.Button(btns,text="Load TLE",command=self.on_load_tle,style='Card.TButton').pack(side=tk.LEFT,padx=4)
        ttk.Button(btns,text="Open WebGL",command=self.open_webgl,style='Card.TButton').pack(side=tk.LEFT,padx=4)
        ttk.Button(btns,text="Refresh WebGL",command=self.open_webgl,style='Card.TButton').pack(side=tk.LEFT,padx=4)
        ttk.Label(btns,text="Duration (s)",style='Card.TLabel').pack(side=tk.LEFT,padx=6)
        self.web_duration = tk.DoubleVar(value=7200.0)
        ttk.Entry(btns,textvariable=self.web_duration,width=8).pack(side=tk.LEFT,padx=4)
        ttk.Label(btns,text="Step (s)",style='Card.TLabel').pack(side=tk.LEFT,padx=6)
        self.web_step = tk.DoubleVar(value=30.0)
        ttk.Entry(btns,textvariable=self.web_step,width=6).pack(side=tk.LEFT,padx=4)
        right = ttk.Frame(self)
        right.grid(row=0,column=1,sticky="nsew")
        right.columnconfigure(0, weight=1)
        right.rowconfigure(0, weight=1)
        self.fig = Figure(figsize=(7,5))
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.ax.set_facecolor('#0f1115')
        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.canvas.get_tk_widget().grid(row=0,column=0,sticky="nsew")
        control = ttk.Frame(right, style='Card.TFrame')
        control.grid(row=1,column=0,sticky="ew")
        ttk.Label(control,text="Time (seconds of week)",style='Card.TLabel').pack(side=tk.LEFT,padx=8)
        self.time_var = tk.DoubleVar(value=0.0)
        self.time_scale = ttk.Scale(control,from_=0,to=3600,variable=self.time_var,command=lambda v:self.draw_scene())
        self.time_scale.pack(side=tk.LEFT,fill=tk.X,expand=True,padx=8,pady=6)
        self.rowconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.draw_scene()

    def on_generate_tle(self):
        try:
            a_m = float(self.entries[0].get())
            e = float(self.entries[1].get())
            inc = float(self.entries[2].get())
            raan = float(self.entries[3].get())
            argp = float(self.entries[4].get())
            m0 = float(self.entries[5].get())
            P = int(self.entries[6].get())
            S = int(self.entries[7].get())
            F = int(self.entries[8].get())
            epoch = datetime.strptime(self.entries[9].get(),"%Y-%m-%d").replace(tzinfo=timezone.utc)
            tles = generate_walker_delta_tles("SAT",a_m,e,inc,raan,argp,m0,epoch,P,S,F)
            self.tle_data = tles
            self.sim = SatelliteSimulation()
            self.sim.tle_satellites = {}
            for i,(name,l1,l2) in enumerate(tles):
                self.sim.add_tle_satellite(name,l1,l2,i+1)
            messagebox.showinfo("Success","TLE generated")
            self.draw_scene()
        except Exception as ex:
            messagebox.showerror("Error",str(ex))

    def on_save_tle(self):
        if not hasattr(self,"tle_data"):
            messagebox.showwarning("Info","Please generate or load TLE first")
            return
        path = filedialog.asksaveasfilename(defaultextension=".tle",filetypes=[("TLE","*.tle")])
        if not path:
            return
        write_tle_file(path,self.tle_data)
        self.tle_file = path
        messagebox.showinfo("Success","TLE saved")

    def on_load_tle(self):
        path = filedialog.askopenfilename(filetypes=[("TLE","*.tle")])
        if not path:
            return
        self.sim = SatelliteSimulation()
        self.sim.load_tle_file(path)
        self.tle_file = path
        self.draw_scene()

    def open_webgl(self):
        if not self.tle_file and not hasattr(self,"tle_data"):
            messagebox.showwarning("Info","Please generate or load TLE first")
            return
        tle_path = self.tle_file
        if not tle_path:
            tmp = filedialog.asksaveasfilename(defaultextension=".tle",filetypes=[("TLE","*.tle")])
            if not tmp:
                return
            write_tle_file(tmp, self.tle_data)
            tle_path = tmp
        start_dt = datetime.utcnow()
        duration_s = int(self.web_duration.get())
        step_s = int(self.web_step.get())
        out_dir = os.path.join(os.path.dirname(__file__), 'web')
        out_czml = os.path.join(out_dir, 'data.czml')
        s = load_settings()
        write_web_theme(s.get('accent','#00c8ff'))
        write_web_config(s.get('ion_token',''), s.get('use_world_terrain', False))
        accent_rgba = hex_to_rgba(s.get('accent','#00c8ff'))
        write_czml(tle_path, start_dt, duration_s, step_s, out_czml, accent_rgba)
        index = os.path.join(out_dir, 'index.html')
        url_base, _srv = start_web(out_dir)
        webbrowser.open_new_tab(f"{url_base}/index.html")
        messagebox.showinfo("Opened","WebGL view opened in browser")

    def draw_scene(self):
        self.ax.clear()
        r = 6378137.0
        u = np.linspace(0,2*np.pi,40)
        v = np.linspace(0,np.pi,20)
        x = r*np.outer(np.cos(u),np.sin(v))
        y = r*np.outer(np.sin(u),np.sin(v))
        z = r*np.outer(np.ones_like(u),np.cos(v))
        self.ax.plot_surface(x,y,z,color='#1f7ae0',alpha=0.25,linewidth=0)
        self.ax.grid(False)
        self.ax.set_xticks([]); self.ax.set_yticks([]); self.ax.set_zticks([])
        self.ax.view_init(elev=20, azim=45)
        t = float(self.time_var.get())
        prns = list(self.sim.tle_satellites.keys())
        for prn in prns:
            pos,vel = self.sim.calculate_satellite_position_tle(prn,t)
            if np.linalg.norm(pos)==0:
                continue
            self.ax.scatter(pos[0],pos[1],pos[2],s=12,c='#00c8ff')
        self.ax.set_box_aspect([1,1,1])
        lim = r*1.2
        self.ax.set_xlim(-lim,lim)
        self.ax.set_ylim(-lim,lim)
        self.ax.set_zlim(-lim,lim)
        self.canvas.draw()

class SimulationPage(ttk.Frame):
    def __init__(self, master):
        super().__init__(master)
        self._build_ui()

    def _build_ui(self):
        form = ttk.Frame(self)
        form.pack(fill=tk.X)
        self.tle_path = tk.StringVar()
        self.out_path = tk.StringVar(value="observations.csv")
        self.lat = tk.DoubleVar(value=40.0)
        self.lon = tk.DoubleVar(value=116.0)
        self.alt = tk.DoubleVar(value=50.0)
        self.start = tk.DoubleVar(value=0.0)
        self.duration = tk.DoubleVar(value=3600.0)
        self.step = tk.DoubleVar(value=10.0)
        self.mask = tk.DoubleVar(value=10.0)
        ttk.Button(form,text="Select TLE",command=self.select_tle).grid(row=0,column=0)
        ttk.Entry(form,textvariable=self.tle_path,width=50).grid(row=0,column=1,sticky="ew")
        ttk.Label(form,text="Start time (s of week)").grid(row=1,column=0)
        ttk.Entry(form,textvariable=self.start).grid(row=1,column=1)
        ttk.Label(form,text="Duration (s)").grid(row=2,column=0)
        ttk.Entry(form,textvariable=self.duration).grid(row=2,column=1)
        ttk.Label(form,text="Step (s)").grid(row=3,column=0)
        ttk.Entry(form,textvariable=self.step).grid(row=3,column=1)
        ttk.Label(form,text="Latitude (deg)").grid(row=4,column=0)
        ttk.Entry(form,textvariable=self.lat).grid(row=4,column=1)
        ttk.Label(form,text="Longitude (deg)").grid(row=5,column=0)
        ttk.Entry(form,textvariable=self.lon).grid(row=5,column=1)
        ttk.Label(form,text="Altitude (m)").grid(row=6,column=0)
        ttk.Entry(form,textvariable=self.alt).grid(row=6,column=1)
        ttk.Label(form,text="Elevation mask (deg)").grid(row=7,column=0)
        ttk.Entry(form,textvariable=self.mask).grid(row=7,column=1)
        ttk.Label(form,text="Output CSV").grid(row=8,column=0)
        ttk.Entry(form,textvariable=self.out_path,width=50).grid(row=8,column=1,sticky="ew")
        ttk.Button(form,text="Start simulation",command=self.run_sim).grid(row=9,column=0,columnspan=2,sticky="ew")
        self.fig2 = Figure(figsize=(6,3))
        self.ax2 = self.fig2.add_subplot(111)
        self.canvas2 = FigureCanvasTkAgg(self.fig2, master=self)
        self.canvas2.get_tk_widget().pack(fill=tk.BOTH,expand=True)

    def select_tle(self):
        p = filedialog.askopenfilename(filetypes=[("TLE","*.tle")])
        if p:
            self.tle_path.set(p)

    def run_sim(self):
        try:
            run_sim_obs(self.tle_path.get(), self.start.get(), self.duration.get(), self.step.get(), self.lat.get(), self.lon.get(), self.alt.get(), self.mask.get(), self.out_path.get())
            import pandas as pd
            df = pd.read_csv(self.out_path.get())
            if not df.empty:
                self.ax2.clear()
                self.ax2.plot(df['time_s'], df['doppler_Hz'])
                self.canvas2.draw()
            messagebox.showinfo("Done","Observations generated")
        except Exception as ex:
            messagebox.showerror("Error",str(ex))

class PositioningPage(ttk.Frame):
    def __init__(self, master):
        super().__init__(master)
        self._build_ui()

    def _build_ui(self):
        form = ttk.Frame(self)
        form.pack(fill=tk.X)
        self.obs = tk.StringVar(value="observations.csv")
        self.nav = tk.StringVar()
        self.tle = tk.StringVar()
        self.use_pr = tk.BooleanVar(value=True)
        self.use_dop = tk.BooleanVar(value=False)
        self.init_x = tk.DoubleVar(value=0.0)
        self.init_y = tk.DoubleVar(value=0.0)
        self.init_z = tk.DoubleVar(value=0.0)
        ttk.Button(form,text="Select observations",command=lambda:self._pick(self.obs)).grid(row=0,column=0)
        ttk.Entry(form,textvariable=self.obs,width=50).grid(row=0,column=1,sticky="ew")
        ttk.Button(form,text="Select ephemeris",command=lambda:self._pick(self.nav)).grid(row=1,column=0)
        ttk.Entry(form,textvariable=self.nav,width=50).grid(row=1,column=1,sticky="ew")
        ttk.Button(form,text="Select TLE",command=lambda:self._pick(self.tle)).grid(row=2,column=0)
        ttk.Entry(form,textvariable=self.tle,width=50).grid(row=2,column=1,sticky="ew")
        ttk.Checkbutton(form,text="Use pseudorange",variable=self.use_pr).grid(row=3,column=0)
        ttk.Checkbutton(form,text="Use Doppler",variable=self.use_dop).grid(row=3,column=1)
        ttk.Label(form,text="Initial X").grid(row=4,column=0)
        ttk.Entry(form,textvariable=self.init_x).grid(row=4,column=1)
        ttk.Label(form,text="Initial Y").grid(row=5,column=0)
        ttk.Entry(form,textvariable=self.init_y).grid(row=5,column=1)
        ttk.Label(form,text="Initial Z").grid(row=6,column=0)
        ttk.Entry(form,textvariable=self.init_z).grid(row=6,column=1)
        ttk.Button(form,text="Start positioning",command=self.run_pos).grid(row=7,column=0,columnspan=2,sticky="ew")
        self.result = tk.StringVar()
        ttk.Label(self,textvariable=self.result).pack(fill=tk.X)

    def _pick(self,var):
        p = filedialog.askopenfilename()
        if p:
            var.set(p)

    def run_pos(self):
        try:
            run_lsq_positioning(self.obs.get(), ephemeris_file=self.nav.get() or None, tle_file=self.tle.get() or None, init_pos=[self.init_x.get(), self.init_y.get(), self.init_z.get()], use_pseudorange=self.use_pr.get(), use_doppler=self.use_dop.get())
            self.result.set("Positioning finished; see console output")
            messagebox.showinfo("Done","Positioning completed")
        except Exception as ex:
            messagebox.showerror("Error",str(ex))

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Satellite Constellation Simulation & Positioning")
        self.geometry("1200x800")
        notebook = ttk.Notebook(self)
        notebook.pack(fill=tk.BOTH,expand=True)
        self.page1 = ConstellationPage(notebook)
        self.page2 = SimulationPage(notebook)
        self.page3 = PositioningPage(notebook)
        notebook.add(self.page1,text="Constellation Design")
        notebook.add(self.page2,text="Observation Simulation")
        notebook.add(self.page3,text="Positioning")

def main():
    app = App()
    app.mainloop()