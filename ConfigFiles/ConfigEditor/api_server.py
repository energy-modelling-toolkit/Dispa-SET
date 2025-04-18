#!/usr/bin/env python3
from http.server import HTTPServer, BaseHTTPRequestHandler
import json
import os
import sys
import yaml
import urllib.parse
import socketserver
import threading
import time

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
        # Handle requests to list YAML files in the ConfigFiles directory
        if self.path == '/list-yaml-files':
            try:
                # List YAML files in parent directory (ConfigFiles)
                yaml_files = []
                parent_dir = '..' 
                if os.path.exists(parent_dir):
                    # List only base filenames ending with .yml or .yaml
                    yaml_files = [f for f in os.listdir(parent_dir) 
                                     if os.path.isfile(os.path.join(parent_dir, f)) and f.endswith(('.yml', '.yaml'))]
                else:
                     # Handle case where parent dir might not be accessible (less likely here)
                    print(f"Warning: Could not access parent directory '{parent_dir}'")
                
                self._set_response_headers()
                self.wfile.write(json.dumps({'files': yaml_files}).encode())
            except Exception as e:
                print(f"Error listing YAML files: {e}") # Log error server-side
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
                
                if not os.path.exists(filename):
                    self._set_response_headers(status_code=404)
                    self.wfile.write(json.dumps({'error': f'File {filename} not found'}).encode())
                    return
                
                with open(filename, 'r') as f:
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
        # Handle requests to save YAML content to a file
        if self.path == '/save-file':
            try:
                content_length = int(self.headers['Content-Length'])
                post_data = self.rfile.read(content_length).decode('utf-8')
                data = json.loads(post_data)
                
                filename = data.get('filename')
                content = data.get('content')
                
                if not filename or not content:
                    self._set_response_headers(status_code=400)
                    self.wfile.write(json.dumps({'error': 'Both filename and content are required'}).encode())
                    return
                
                # Ensure we only save YAML files within the current directory or parent ConfigFiles directory for security
                if ('/..' in filename or filename.startswith('/') or '\\' in filename or filename.endswith('/..') or
                    (filename.startswith('..') and not filename.startswith('../'))):
                    self._set_response_headers(status_code=400)
                    self.wfile.write(json.dumps({'error': 'Invalid filename'}).encode())
                    return
                
                # Resolve the path to make sure it's within ConfigFiles directory
                abs_path = os.path.abspath(filename)
                config_files_dir = os.path.abspath(os.path.join(os.getcwd(), '..'))
                if not abs_path.startswith(config_files_dir):
                    self._set_response_headers(status_code=400)
                    self.wfile.write(json.dumps({'error': 'Path traversal not allowed'}).encode())
                    return
                
                # Add .yml extension if not present
                if not filename.endswith(('.yml', '.yaml')):
                    filename = f"{filename}.yml"
                
                # Parse and dump the YAML to ensure it's valid and properly formatted
                try:
                    yaml_data = yaml.safe_load(content)
                    
                    # Special handling for date fields - convert arrays to Python tuples
                    if 'StartDate' in yaml_data and isinstance(yaml_data['StartDate'], list):
                        yaml_data['StartDate'] = tuple(yaml_data['StartDate'])
                    
                    if 'StopDate' in yaml_data and isinstance(yaml_data['StopDate'], list):
                        yaml_data['StopDate'] = tuple(yaml_data['StopDate'])
                    
                    with open(filename, 'w') as f:
                        yaml.dump(yaml_data, f, sort_keys=False)
                    
                    self._set_response_headers()
                    self.wfile.write(json.dumps({'success': True, 'message': f'File {filename} saved successfully'}).encode())
                except Exception as e:
                    self._set_response_headers(status_code=400)
                    self.wfile.write(json.dumps({'error': f'Invalid YAML: {str(e)}'}).encode())
                    
            except Exception as e:
                self._set_response_headers(status_code=500)
                self.wfile.write(json.dumps({'error': str(e)}).encode())
        else:
            self._set_response_headers(status_code=404)
            self.wfile.write(json.dumps({'error': 'Endpoint not found'}).encode())

def start_api_server(port=API_PORT):
    """Start the API server in a separate thread."""
    try:
        # Use ThreadingHTTPServer for better handling of concurrent requests
        class ThreadingHTTPServer(socketserver.ThreadingMixIn, HTTPServer):
            pass
        
        server_address = ('', port)
        httpd = ThreadingHTTPServer(server_address, YAMLHandler)
        
        print(f"Starting API server on http://localhost:{port}")
        api_thread = threading.Thread(target=httpd.serve_forever, daemon=True)
        api_thread.start()
        return httpd, api_thread
    
    except Exception as e:
        print(f"Error starting API server: {e}")
        return None, None

if __name__ == "__main__":
    # If running this script directly, start the API server
    httpd, thread = start_api_server()
    
    # Keep the main thread running to prevent the daemon thread from exiting
    try:
        print("API server is running. Press Ctrl+C to stop.")
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        print("Stopping API server...")
        httpd.shutdown()
        print("API server stopped.") 