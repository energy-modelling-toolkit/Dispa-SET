#!/usr/bin/env python3
from http.server import HTTPServer, BaseHTTPRequestHandler
import json
import os
import sys
import urllib.parse
import socketserver
import threading
import time

# Try to import yaml, but allow the server to start without it for reading
try:
    import yaml
    HAS_YAML = True
except ImportError:
    print("Warning: 'yaml' (PyYAML) library not found. Saving functionality will be limited.")
    HAS_YAML = False

# Default port for the API server
API_PORT = 8084

class YAMLHandler(BaseHTTPRequestHandler):
    def _set_response_headers(self, content_type='application/json', status_code=200):
        self.send_response(status_code)
        self.send_header('Content-type', content_type)
        self.send_header('Access-Control-Allow-Origin', '*')  # Allow CORS for local development
        self.send_header('Access-Control-Allow-Methods', 'GET, POST, OPTIONS')
        self.send_header('Access-Control-Allow-Headers', 'Content-Type')
        self.end_headers()
    
    def do_OPTIONS(self):
        # Handle preflight requests for CORS
        self._set_response_headers()
        
    def do_GET(self):
        # Determine the directory containing the config files (parent of ConfigEditor)
        # We use an absolute path for robustness
        config_files_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

        # Handle requests to list YAML files in the ConfigFiles directory
        if self.path == '/list-yaml-files':
            try:
                yaml_files = []
                if os.path.exists(config_files_dir):
                    # List only base filenames ending with .yml or .yaml
                    yaml_files = [f for f in os.listdir(config_files_dir) 
                                     if os.path.isfile(os.path.join(config_files_dir, f)) and f.endswith(('.yml', '.yaml'))]
                
                self._set_response_headers()
                self.wfile.write(json.dumps({'files': yaml_files}).encode())
            except Exception as e:
                print(f"Error listing YAML files: {e}")
                self._set_response_headers(status_code=500)
                self.wfile.write(json.dumps({'error': 'Error listing configuration files'}).encode())
                
        # Handle requests to read a specific YAML file
        elif self.path.startswith('/read-file?'):
            try:
                query = urllib.parse.parse_qs(urllib.parse.urlparse(self.path).query)
                filename = query.get('filename', [''])[0]
                
                if not filename:
                    self._set_response_headers(status_code=400)
                    self.wfile.write(json.dumps({'error': 'Filename parameter is required'}).encode())
                    return
                
                # Resolve the file path relative to config_files_dir
                # If filename starts with '../', we remove it as we are already at the parent
                clean_filename = filename
                if filename.startswith('../'):
                    clean_filename = filename[3:]
                
                full_path = os.path.abspath(os.path.join(config_files_dir, clean_filename))
                
                # Security check: ensure the file is within config_files_dir
                if not full_path.startswith(config_files_dir):
                    self._set_response_headers(status_code=403)
                    self.wfile.write(json.dumps({'error': 'Access denied: Path outside config directory'}).encode())
                    return

                if not os.path.exists(full_path):
                    self._set_response_headers(status_code=404)
                    self.wfile.write(json.dumps({'error': f'File {filename} not found'}).encode())
                    return
                
                with open(full_path, 'r') as f:
                    yaml_content = f.read()
                
                self._set_response_headers(content_type='text/yaml')
                self.wfile.write(yaml_content.encode())
                
            except Exception as e:
                self._set_response_headers(status_code=500)
                self.wfile.write(json.dumps({'error': str(e)}).encode())
        else:
            self._set_response_headers(status_code=404)
            self.wfile.write(json.dumps({'error': 'Endpoint not found'}).encode())
    
    def do_POST(self):
        # Determine the directory containing the config files
        # We use the absolute path of the script's directory's parent
        script_dir = os.path.dirname(os.path.abspath(__file__))
        config_files_dir = os.path.abspath(os.path.join(script_dir, '..'))

        # Handle requests to save YAML content to a file
        if self.path == '/save-file':
            try:
                if not HAS_YAML:
                    self._set_response_headers(status_code=500)
                    self.wfile.write(json.dumps({'error': "The 'yaml' library is not installed on the server. Saving is disabled."}).encode())
                    return

                content_length = int(self.headers.get('Content-Length', 0))
                if content_length == 0:
                    self._set_response_headers(status_code=400)
                    self.wfile.write(json.dumps({'error': 'Empty request body'}).encode())
                    return

                post_data = self.rfile.read(content_length).decode('utf-8')
                data = json.loads(post_data)
                
                filename = data.get('filename')
                content = data.get('content')
                
                if not filename or not content:
                    self._set_response_headers(status_code=400)
                    self.wfile.write(json.dumps({'error': 'Both filename and content are required'}).encode())
                    return
                
                # Resolve the path safely
                # Remove potential '../' or './' prefixes
                clean_filename = filename
                while clean_filename.startswith('../') or clean_filename.startswith('./'):
                    clean_filename = clean_filename.split('/', 1)[-1]
                
                # We only allow saving into the ConfigFiles directory
                full_path = os.path.abspath(os.path.join(config_files_dir, clean_filename))
                
                # Security check: ensure the resolved path is inside config_files_dir
                # We use commonpath to be 100% sure
                try:
                    if os.path.commonpath([config_files_dir, full_path]) != config_files_dir:
                        print(f"Path security violation: {full_path} is not in {config_files_dir}")
                        self._set_response_headers(status_code=403)
                        self.wfile.write(json.dumps({'error': f'Access denied: Path {clean_filename} is outside allowed directory'}).encode())
                        return
                except ValueError: # Paths on different drives or invalid
                    self._set_response_headers(status_code=403)
                    self.wfile.write(json.dumps({'error': 'Access denied: Invalid path'}).encode())
                    return
                
                # Add .yml extension if not present and no extension exists
                if not os.path.splitext(full_path)[1]:
                    full_path = f"{full_path}.yml"
                
                print(f"Attempting to save to: {full_path}")

                # Parse and dump the YAML to ensure it's valid and properly formatted
                try:
                    yaml_data = yaml.safe_load(content)
                    
                    # Special handling for date fields - convert arrays to Python tuples
                    if 'StartDate' in yaml_data and isinstance(yaml_data['StartDate'], list):
                        yaml_data['StartDate'] = tuple(yaml_data['StartDate'])
                    
                    if 'StopDate' in yaml_data and isinstance(yaml_data['StopDate'], list):
                        yaml_data['StopDate'] = tuple(yaml_data['StopDate'])
                    
                    # Ensure the directory exists
                    os.makedirs(os.path.dirname(full_path), exist_ok=True)

                    with open(full_path, 'w') as f:
                        yaml.dump(yaml_data, f, sort_keys=False)
                    
                    print(f"Successfully saved {full_path}")
                    self._set_response_headers()
                    self.wfile.write(json.dumps({'success': True, 'message': f'File {os.path.basename(full_path)} saved successfully'}).encode())
                except Exception as e:
                    print(f"YAML dump error: {e}")
                    self._set_response_headers(status_code=400)
                    self.wfile.write(json.dumps({'error': f'Invalid YAML or write error: {str(e)}'}).encode())
                    
            except Exception as e:
                print(f"General POST error: {e}")
                self._set_response_headers(status_code=500)
                self.wfile.write(json.dumps({'error': str(e)}).encode())
        else:
            self._set_response_headers(status_code=404)
            self.wfile.write(json.dumps({'error': 'Endpoint not found'}).encode())

def start_api_server(port=API_PORT, max_attempts=20):
    """Start the API server in a separate thread, trying consecutive ports if needed."""
    httpd = None
    api_thread = None
    used_port = None

    # Use ThreadingHTTPServer for better handling of concurrent requests
    class ThreadingHTTPServer(socketserver.ThreadingMixIn, HTTPServer):
        pass

    for attempt in range(max_attempts):
        current_port = port + attempt
        try:
            server_address = ('', current_port)
            httpd = ThreadingHTTPServer(server_address, YAMLHandler)
            used_port = current_port
            print(f"Successfully bound API server to port {used_port}.")
            break # Exit loop if successful
        except OSError as e:
            if "Address already in use" in str(e):
                print(f"Info: API Port {current_port} is busy, trying next...")
                if attempt == max_attempts - 1: # Last attempt failed
                    print(f"\nError: Could not find an open API port after {max_attempts} attempts starting from {port}.")
                    return None, None, None
            else:
                print(f"\nError starting API server on port {current_port}: {e}")
                return None, None, None # Other OS error

    if httpd is None or used_port is None:
        print(f"\nError: Failed to bind API server to any port between {port} and {port + max_attempts - 1}.")
        return None, None, None

    print(f"Starting API server on http://localhost:{used_port}")
    api_thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    api_thread.start()
    return httpd, api_thread, used_port

if __name__ == "__main__":
    # If running this script directly, start the API server
    httpd, thread, used_port = start_api_server()

    if httpd:
        # Keep the main thread running to prevent the daemon thread from exiting
        try:
            print(f"API server is running on port {used_port}. Press Ctrl+C to stop.")
            while True:
                time.sleep(1)
        except KeyboardInterrupt:
            print("Stopping API server...")
            httpd.shutdown()
            print("API server stopped.")
    else:
        print("API server failed to start.") 