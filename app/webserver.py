import threading
from http.server import ThreadingHTTPServer, SimpleHTTPRequestHandler
import socket
import os
import sys

class _Handler(SimpleHTTPRequestHandler):
    def __init__(self, *args, web_dir=None, data_dir=None, **kwargs):
        self.web_dir = web_dir
        self.data_dir = data_dir
        # Don't pass directory to parent, we'll handle it ourselves
        super().__init__(*args, **kwargs)
    
    def translate_path(self, path):
        """Override to support files from both web_dir and data_dir"""
        # Remove query string
        path = path.split('?')[0]
        # Remove leading slash
        path = path.lstrip('/')
        
        # Check data directory first (for CZML, config, theme files)
        if self.data_dir and path in ['data.czml', 'config.json', 'theme.json', 'temp.tle']:
            data_path = os.path.join(self.data_dir, path)
            if os.path.exists(data_path):
                return data_path
        
        # Fall back to web directory
        if self.web_dir:
            web_path = os.path.join(self.web_dir, path)
            if os.path.exists(web_path):
                return web_path
            # Also check for index.html
            if path == '' or path == 'index.html':
                index_path = os.path.join(self.web_dir, 'index.html')
                if os.path.exists(index_path):
                    return index_path
            
            # Debug: log missing files in frozen mode
            if getattr(sys, 'frozen', False) and not os.path.exists(web_path):
                # Only log first few missing files to avoid spam
                if not hasattr(self, '_missing_logged'):
                    self._missing_logged = set()
                if len(self._missing_logged) < 5 and path not in self._missing_logged:
                    self._missing_logged.add(path)
                    print(f"[WebServer] File not found: {path}")
                    print(f"  Expected at: {web_path}")
                    print(f"  web_dir exists: {os.path.exists(self.web_dir)}")
                    if os.path.exists(self.web_dir):
                        # List some files in web_dir to verify it's correct
                        try:
                            files = os.listdir(self.web_dir)[:5]
                            print(f"  web_dir contents (first 5): {files}")
                        except:
                            pass
        
        # If not found, return a path that will cause 404 (better than crashing)
        return os.path.join(self.web_dir or '.', path)
    
    def log_message(self, format, *args):
        """Override to suppress normal request logs, but keep errors"""
        # Only log errors
        if '404' in format % args or '500' in format % args:
            print(f"[WebServer] {format % args}")

def get_data_dir():
    """Get writable directory for web data files"""
    if getattr(sys, 'frozen', False):
        # In frozen mode, use temp directory
        import tempfile
        temp_dir = tempfile.gettempdir()
        data_dir = os.path.join(temp_dir, 'LEO-SPINE', 'web_data')
        os.makedirs(data_dir, exist_ok=True)
        return data_dir
    else:
        # In dev mode, use app/web directory
        return os.path.join(os.path.dirname(__file__), 'web')

def start(root_dir, data_dir=None):
    """
    Start web server.
    
    Args:
        root_dir: Web root directory (for static files)
        data_dir: Data directory (for CZML, config files) - optional
    
    Returns:
        (url, server) tuple
    """
    if data_dir is None:
        data_dir = get_data_dir()
    
    # Verify directories exist
    if not os.path.exists(root_dir):
        raise FileNotFoundError(f"Web directory not found: {root_dir}")
    
    index_path = os.path.join(root_dir, 'index.html')
    if not os.path.exists(index_path):
        raise FileNotFoundError(f"index.html not found in: {root_dir}")
    
    # Verify Cesium resources exist (important for frozen mode)
    cesium_js = os.path.join(root_dir, 'libs', 'cesium', 'Build', 'Cesium', 'Cesium.js')
    if not os.path.exists(cesium_js):
        if getattr(sys, 'frozen', False):
            # In frozen mode, this is a critical error
            print(f"[WebServer] WARNING: Cesium.js not found at: {cesium_js}")
            print(f"  root_dir = {root_dir}")
            print(f"  root_dir exists: {os.path.exists(root_dir)}")
            if os.path.exists(root_dir):
                try:
                    contents = os.listdir(root_dir)
                    print(f"  root_dir contents: {contents[:10]}")
                    libs_path = os.path.join(root_dir, 'libs')
                    if os.path.exists(libs_path):
                        libs_contents = os.listdir(libs_path)
                        print(f"  libs/ contents: {libs_contents[:5]}")
                except Exception as e:
                    print(f"  Error listing directory: {e}")
            # Don't raise error, but log warning - might still work if path is different
    
    # Find available port
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    sock.bind(('127.0.0.1', 0))
    addr, port = sock.getsockname()
    sock.close()
    
    # Create server
    def handler_factory(*args, **kwargs):
        return _Handler(*args, web_dir=root_dir, data_dir=data_dir, **kwargs)
    
    try:
        server = ThreadingHTTPServer(('127.0.0.1', port), handler_factory)
        t = threading.Thread(target=server.serve_forever, daemon=True)
        t.start()
        
        # Verify server started
        import time
        time.sleep(0.1)  # Give server time to start
        
        url = f"http://127.0.0.1:{port}"
        if getattr(sys, 'frozen', False):
            print(f"[WebServer] Started on {url}, web_dir={root_dir}, data_dir={data_dir}")
        
        return url, server
    except Exception as e:
        raise RuntimeError(f"Failed to start web server: {e}") from e