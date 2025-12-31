import json
import os

def settings_path():
    return os.path.join(os.path.dirname(__file__), 'user_settings.json')

def load_settings():
    p = settings_path()
    if os.path.exists(p):
        with open(p, 'r', encoding='utf-8') as f:
            return json.load(f)
    return {"theme":"light","accent":"#3b82f6","ion_token":"","use_world_terrain":False}

def save_settings(s):
    p = settings_path()
    with open(p, 'w', encoding='utf-8') as f:
        json.dump(s, f)

def qss_path(theme):
    base = os.path.join(os.path.dirname(__file__), 'styles')
    t = 'dark.qss' if theme == 'dark' else 'light.qss'
    return os.path.join(base, t)

def accent_override_qss(accent):
    # Generate accent color light/dark variants
    h = accent.lstrip('#')
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    # Darker variant
    darker = '#%02x%02x%02x' % (max(0, r-30), max(0, g-30), max(0, b-30))
    # Lighter variant  
    lighter = '#%02x%02x%02x' % (min(255, r+40), min(255, g+40), min(255, b+40))
    
    return f'''
/* Accent color overrides */
QPushButton#primaryBtn {{
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1, stop:0 {accent}, stop:1 {darker});
}}
QPushButton#primaryBtn:hover {{
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1, stop:0 {lighter}, stop:1 {accent});
}}
QTabBar::tab:selected {{
    border-bottom: 2px solid {accent};
}}
QGroupBox::title {{
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 {accent}, stop:1 {darker});
}}
QCheckBox::indicator:checked {{
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1, stop:0 {accent}, stop:1 {darker});
    border-color: {accent};
}}
QCheckBox::indicator:checked:hover {{
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1, stop:0 {lighter}, stop:1 {accent});
}}
QRadioButton::indicator:checked {{
    border-color: {accent};
}}
QLineEdit:focus {{
    border-color: {accent};
}}
QComboBox:hover {{
    border-color: {accent};
}}
QDateTimeEdit:focus {{
    border-color: {accent};
}}
QProgressBar::chunk {{
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 {accent}, stop:1 {lighter});
}}
QMenu::item:selected {{
    background-color: {accent};
}}
'''

def apply_qss(app, theme, accent):
    p = qss_path(theme)
    with open(p, 'r', encoding='utf-8') as f:
        base = f.read()
    app.setStyleSheet(base + accent_override_qss(accent))

def get_web_data_path(filename):
    """Get path for web data files that need to be writable"""
    import sys
    import tempfile
    if getattr(sys, 'frozen', False):
        # In frozen mode, use temp directory
        temp_dir = tempfile.gettempdir()
        data_dir = os.path.join(temp_dir, 'LEO-SPINE', 'web_data')
        os.makedirs(data_dir, exist_ok=True)
        return os.path.join(data_dir, filename)
    else:
        # In dev mode, use app/web directory
        return os.path.join(os.path.dirname(__file__), 'web', filename)

def write_web_theme(accent_hex):
    out = get_web_data_path('theme.json')
    data = {"accent": accent_hex, "panelBg": "#12161c", "panelBorder": "#1d2026", "fg": "#e6e6e6", "bg": "#0f1115"}
    with open(out, 'w', encoding='utf-8') as f:
        json.dump(data, f)

def write_web_config(ion_token, use_world_terrain, duration_s=None, start_iso=None):
    out = get_web_data_path('config.json')
    data = {
        "ionToken": ion_token or "",
        "useWorldTerrain": bool(use_world_terrain)
    }
    if duration_s is not None:
        data["durationSeconds"] = duration_s
    if start_iso is not None:
        data["startIso"] = start_iso
    with open(out, 'w', encoding='utf-8') as f:
        json.dump(data, f)

def hex_to_rgba(hex_color):
    h = hex_color.lstrip('#')
    r = int(h[0:2], 16)
    g = int(h[2:4], 16)
    b = int(h[4:6], 16)
    return [r, g, b, 255]