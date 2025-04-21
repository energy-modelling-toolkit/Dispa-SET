import streamlit as st
import dispaset as ds
import pandas as pd
import os
import datetime
# Import the main UI builder function
from dispaset.postprocessing.streamlit_functions import build_streamlit_app

# --- Configuration ---
st.set_page_config(layout="wide")
st.title("Dispa-SET Postprocessing Dashboard")

# --- Data Loading (Keep Cached Functions Here) ---
@st.cache_data
def load_data(sim_path):
    """Loads simulation inputs and results using ds.get_sim_results."""
    try:
        inputs, results = ds.get_sim_results(path=sim_path, cache=False)
        if not inputs or not results:
            st.error(f"Failed to load data from {sim_path}. Check path and simulation files.")
            return None, None
        return inputs, results
    except FileNotFoundError:
        st.error(f"Directory not found: {sim_path}")
        return None, None
    except Exception as e:
        st.error(f"Error loading data from {sim_path}: {e}")
        # import traceback
        # st.error(traceback.format_exc()) # Uncomment for detailed error
        return None, None

@st.cache_data
def get_pp_indicators(_inputs, _results):
    """Gets power plant indicators, cached."""
    if _inputs and _results:
        try:
            # Check if necessary data exists in results
            if 'OutputPower' not in _results:
                 st.warning("Cannot calculate indicators: 'OutputPower' missing from results.")
                 return None
            if 'OutputCommitted' not in _results:
                 st.warning("Cannot calculate indicators: 'OutputCommitted' missing from results.")
                 return None # Or handle differently if some indicators can still be calculated
            indicators = ds.get_indicators_powerplant(_inputs, _results)
            return indicators
        except KeyError as e:
             st.error(f"Error calculating indicators: Missing key {e}. Results might be incomplete.")
             return None
        except Exception as e:
            st.error(f"Error calculating power plant indicators: {e}")
            # import traceback
            # st.error(traceback.format_exc())
            return None
    return None


# --- Sidebar Controls ---
st.sidebar.header("Simulation Setup")
default_sim_path = "Simulations/simulation_test"
if not os.path.isdir(default_sim_path): default_sim_path = ""
simulation_path = st.sidebar.text_input("Simulation Results Path", default_sim_path)

# Initialize state variables
inputs = None
results = None
indicators = None
data_loaded = False # Flag to track loading status

if simulation_path:
    # Automatically try to load data when path changes
    inputs, results = load_data(simulation_path)
    if inputs and results:
        indicators = get_pp_indicators(inputs, results)
        data_loaded = True # Set flag if loading successful
        if indicators is None:
             st.sidebar.warning("Loaded data, but could not calculate indicators.")
        else:
             st.sidebar.success("Data loaded successfully.") # Moved success message here
    else:
         # Error message is handled within load_data
         data_loaded = False

else:
    st.sidebar.warning("Please enter a valid simulation path.")


# --- Main Area ---
if data_loaded:
    # Call the main function from the other module to build the UI
    build_streamlit_app(inputs, results, indicators)
elif simulation_path:
     st.info("Attempting to load data... Check sidebar for status or error messages.")
else:
    st.info("Enter a simulation path in the sidebar to load data.") 