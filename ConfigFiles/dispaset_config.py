#!/usr/bin/env python3
import http.server
import socketserver
import webbrowser
import os
import sys
import threading
import time
import importlib
import importlib.util

# --- Configuration ---
PORT = 8083
API_PORT = 8084
HTML_FILE = "dispaset_config_editor.html"
CONFIG_DIR = "ConfigEditor"  # Corrected path: relative to this script

def get_script_directory():
    """Gets the directory where the script is running."""
    try:
        # Works when running as a script
        return os.path.dirname(os.path.abspath(__file__))
    except NameError:
        # Fallback for interactive interpreters or environments where __file__ is not defined
        return os.path.abspath(os.path.dirname(sys.argv[0]))

def start_server(target_dir, port):
    """Starts the HTTP server in a separate thread."""
    os.chdir(target_dir)
    print(f"Changed working directory to: {os.getcwd()}")
    Handler = http.server.SimpleHTTPRequestHandler
    
    # Try to create the server, handle potential "Address already in use" error
    try:
        httpd = socketserver.TCPServer(("", port), Handler)
    except OSError as e:
        if "Address already in use" in str(e):
            print(f"\nError: Port {port} is already in use.")
            print("Please close the other program using the port or choose a different port.")
        else:
            print(f"\nError starting server: {e}")
        return None  # Indicate failure

    print(f"Serving files from '{target_dir}' on http://localhost:{port}")
    print(f"Target URL: http://localhost:{port}/{HTML_FILE}")
    
    # Start the server in a separate thread
    server_thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    server_thread.start()
    return httpd, server_thread

def start_api_server(script_dir, config_dir):
    """Import and start the API server from api_server.py"""
    try:
        # First, ensure we're in the correct directory
        original_dir = os.getcwd()
        os.chdir(script_dir)
        
        # Try to import the module using importlib
        api_module = None
        api_server_path = os.path.join(script_dir, config_dir, 'api_server.py')
        
        if os.path.exists(api_server_path):
            spec = importlib.util.spec_from_file_location("api_server", api_server_path)
            api_module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(api_module)
        else:
            print(f"API server module not found at {api_server_path}")
            return None, None
        
        # Change to the config directory where the API server will serve files from
        os.chdir(os.path.join(script_dir, config_dir))
        
        # Start the API server
        api_server, api_thread = api_module.start_api_server(API_PORT)
        
        # Restore original directory
        os.chdir(original_dir)
        
        return api_server, api_thread
    
    except Exception as e:
        print(f"Error starting API server: {e}")
        return None, None


if __name__ == "__main__":
    # Determine the absolute path to the script directory and ConfigFiles directory
    script_dir = get_script_directory()
    config_path = os.path.join(script_dir, CONFIG_DIR)

    if not os.path.isdir(config_path):
        print(f"Error: The configuration directory '{config_path}' was not found.")
        print(f"Please ensure this script is in the parent directory of '{CONFIG_DIR}'.")
        sys.exit(1)

    # Start the API server first
    api_server, api_thread = start_api_server(script_dir, CONFIG_DIR)
    if not api_server:
        print("Warning: API server failed to start. File saving functionality will not work.")
    else:
        print(f"API server started on http://localhost:{API_PORT}")

    # Start the static file server
    server_instance, server_thread = start_server(config_path, PORT)
    if not server_instance:
        print("Server failed to start. Exiting.")
        sys.exit(1)

    # Construct the URL to open
    target_url = f"http://localhost:{PORT}/{HTML_FILE}"

    # Use threading to attempt opening the browser shortly after starting the server
    browser_opened = False
    def open_browser_after_delay():
        global browser_opened
        time.sleep(1)  # Short delay to give the server a moment to start listening
        print(f"Attempting to open {target_url} in your browser...")
        try:
           opened = webbrowser.open_new_tab(target_url)
           if opened:
               print("Browser tab opened (or attempted).")
               browser_opened = True
           else:
               print("webbrowser.open_new_tab returned False. Could not open browser automatically.")
        except Exception as e:
           print(f"Could not open browser automatically: {e}")

    # Start the browser-opening thread
    browser_thread = threading.Thread(target=open_browser_after_delay, daemon=True)
    browser_thread.start()

    # Keep the main thread running
    print("Press Ctrl+C in this terminal to stop the servers.")
    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        print("\nKeyboard interrupt received, shutting down servers.")
        
        # Shut down servers
        if server_instance:
            server_instance.shutdown()
            print("File server stopped.")
        
        if api_server:
            api_server.shutdown()
            print("API server stopped.")
            
    print("All servers stopped. Script finished.") 