import streamlit as st
import dispaset as ds
import pandas as pd
import matplotlib.pyplot as plt
import os
import datetime # Added for date input
import plotly.express as px  # Import Plotly Express

# Import specific functions needed
from dispaset.postprocessing.postprocessing import filter_by_zone
from dispaset.postprocessing.postprocessing import get_imports
from dispaset.postprocessing.postprocessing import filter_by_tech_list
from dispaset.postprocessing.plot import plot_rug, plot_dispatchX # Added plot_dispatchX
# Import commons for colors
from dispaset.common import commons

# --- Configuration ---
st.set_page_config(layout="wide")
st.title("Dispa-SET Postprocessing Dashboard")

# --- Data Loading ---
@st.cache_data
def load_data(sim_path):
    """Loads simulation inputs and results using ds.get_sim_results."""
    try:
        inputs, results = ds.get_sim_results(path=sim_path, cache=False) # Let Streamlit handle caching
        return inputs, results
    except Exception as e:
        st.error(f"Error loading data from {sim_path}: {e}")
        return None, None

@st.cache_data
def get_pp_indicators(_inputs, _results):
    """Gets power plant indicators, cached."""
    if _inputs and _results:
        try:
            indicators = ds.get_indicators_powerplant(_inputs, _results)
            return indicators
        except Exception as e:
            st.error(f"Error calculating power plant indicators: {e}")
            return None
    return None

# --- Sidebar Controls ---
st.sidebar.header("Simulation Setup")
# Default to a common simulation path if it exists, otherwise leave empty
default_sim_path = "Simulations/simulation_test"
if not os.path.isdir(default_sim_path):
    default_sim_path = ""

simulation_path = st.sidebar.text_input("Simulation Results Path", default_sim_path)

inputs = None
results = None
indicators = None

if simulation_path:
    if st.sidebar.button("Load Data"):
        inputs, results = load_data(simulation_path)
        if inputs and results:
            st.sidebar.success(f"Data loaded successfully from {simulation_path}")
            # Attempt to load indicators immediately after loading data
            indicators = get_pp_indicators(inputs, results)
            if indicators is not None:
                 st.sidebar.info("Power plant indicators calculated.")
            else:
                 st.sidebar.warning("Could not calculate power plant indicators.")
        else:
            st.sidebar.error("Failed to load data.")
    # --- Try loading data from cache if path is set (avoids clicking button always) ---
    # This relies on the fact that text_input value change reruns the script
    # and @st.cache_data will return cached data if args match.
    # We call load_data directly here. If the button was pressed, it will
    # return the same cached data. If not pressed but path exists in cache, it loads.
    _inputs_cached, _results_cached = load_data(simulation_path)
    if _inputs_cached and _results_cached:
        inputs = _inputs_cached
        results = _results_cached
        # Get indicators from cache as well
        indicators = get_pp_indicators(inputs, results)


else:
    st.sidebar.warning("Please enter a valid simulation path.")

# --- Main Area ---
if inputs and results and indicators is not None:
    st.header("Energy Balance per Zone")

    all_zones = sorted(inputs['sets']['n'])
    default_zones = all_zones # Default to all zones

    selected_zones = st.multiselect(
        "Select Zones to Display",
        options=sorted(inputs['sets']['n']),
        default=sorted(inputs['sets']['n']) # Default to all zones
    )

    if selected_zones:
        st.subheader("Energy Balance per Zone (Plotly)")
        try:
            # --- Data Preparation (adapted from ds.plot_energy_zone_fuel) ---
            zones_in_results = indicators.Zone.unique()
            fuels = indicators[[u in commons['Technologies'] for u in indicators['Technology']]].Fuel.unique()
            GenPerZone = pd.DataFrame(0.0, index=zones_in_results, columns=list(fuels) + ['FlowIn', 'Curtailment', 'ShedLoad'])

            # Ensure all common fuels are present
            for f in commons['Fuels'] + ['FlowIn']:
                if f not in GenPerZone.columns:
                    GenPerZone[f] = 0.0

            # Aggregate generation and calculate imports/demand/storage input
            power_consumption = pd.DataFrame(0.0, index=results['OutputPowerConsumption'].index, columns=zones_in_results)
            storage_input = pd.DataFrame(0.0, index=results['OutputPower'].index, columns=zones_in_results)

            for z in zones_in_results:
                if z in selected_zones:  # Only process selected zones
                    for f in fuels:
                        tmp = indicators[(indicators.Fuel == f) & (indicators.Zone == z)]
                        GenPerZone.loc[z, f] = tmp.Generation.sum()
                    net_imports = get_imports(results['OutputFlow'], z)
                    if net_imports > 0:
                        GenPerZone.loc[z, 'FlowIn'] = net_imports

                    # Power Consumption (P2X etc.)
                    p_cons_zone = filter_by_zone(results['OutputPowerConsumption'], inputs, z)
                    if not p_cons_zone.empty:
                        power_consumption[z] = p_cons_zone.sum(axis=1)

                    # Storage Input (Charging)
                    storage_input_techs = filter_by_tech_list(results['OutputStorageInput'], inputs, commons['tech_storage'])
                    s_in_zone = filter_by_zone(storage_input_techs, inputs, z)
                    if not s_in_zone.empty:
                        storage_input[z] = s_in_zone.sum(axis=1)

                    # Curtailment & ShedLoad
                    if z in results['OutputCurtailedPower']:
                         GenPerZone.loc[z, 'Curtailment'] = -results['OutputCurtailedPower'][z].sum() # Store as negative
                    if z in results['OutputShedLoad']:
                         GenPerZone.loc[z, 'ShedLoad'] = results['OutputShedLoad'][z].sum()

            # Convert to TWh
            GenPerZone = GenPerZone / 1E6

            # Filter only selected zones for plotting
            GenPerZone = GenPerZone.loc[selected_zones]

            # Select and order columns based on MeritOrder (plus Curtailment/ShedLoad)
            plot_cols_gen = [col for col in commons['MeritOrder'] if col in GenPerZone.columns and GenPerZone[col].sum() != 0]
            plot_cols_other = ['Curtailment', 'ShedLoad']
            GenPerZone_plot = GenPerZone[plot_cols_gen + plot_cols_other].copy()

            # Calculate Demand components (TWh)
            demand_da = inputs['param_df']['Demand']['DA'].sum(axis=0).reindex(selected_zones).fillna(0) / 1E6
            p_cons_total = power_consumption.sum(axis=0).reindex(selected_zones).fillna(0) / 1E6
            s_in_total = storage_input.sum(axis=0).reindex(selected_zones).fillna(0) / 1E6

            total_demand_line = demand_da + p_cons_total # DA + P2X
            storage_demand_line = total_demand_line + s_in_total # DA + P2X + Storage Input
            final_demand_line = storage_demand_line + GenPerZone_plot['ShedLoad'] # DA + P2X + StorageInput + ShedLoad

            # Reshape for Plotly (Handle positive and negative separately for stacking)
            GenPerZone_plot['TotalGen'] = GenPerZone_plot[plot_cols_gen].sum(axis=1) # Sum positive generation
            GenPerZone_plot = GenPerZone_plot.drop(columns=['ShedLoad']) # Remove ShedLoad, will be plotted differently

            # Reset index and rename the new column to 'Zone'
            GenPerZone_plot_reset = GenPerZone_plot.reset_index()
            GenPerZone_plot_reset = GenPerZone_plot_reset.rename(columns={'index': 'Zone'})

            energy_long = GenPerZone_plot_reset.melt(
                 id_vars=['Zone', 'TotalGen'], # Now 'Zone' exists
                 var_name='Source', # Fuel or Curtailment
                 value_name='Energy_TWh'
             )
            # Separate positive and negative values
            energy_long_pos = energy_long[energy_long['Energy_TWh'] >= 0].copy()
            energy_long_neg = energy_long[energy_long['Energy_TWh'] < 0].copy()

            # --- Re-sanitize colors here in case capacity plot wasn't run ---
            sanitized_colors = {}
            for fuel, color in commons['colors'].items():
                if isinstance(color, str):
                    if color.startswith('#') and len(color) == 9: # Check for #RRGGBBAA
                        sanitized_colors[fuel] = color[:7] # Truncate to #RRGGBB
                    elif color.startswith('#') and len(color) == 7: # Check for #RRGGBB
                        sanitized_colors[fuel] = color
                    # Add checks for other valid formats if needed
            # Ensure curtailment has a color
            if 'curtailment' not in sanitized_colors:
                 sanitized_colors['curtailment'] = 'red'

            # --- Plotting with Plotly --- 
            fig_energy = px.bar(
                energy_long_pos,
                x='Zone',
                y='Energy_TWh',
                color='Source',
                title='Energy Balance per Zone (TWh)',
                labels={'Energy_TWh': 'Energy [TWh]'},
                color_discrete_map=sanitized_colors # Use locally sanitized colors
            )

            # Add negative bars (Curtailment) - might need adjustment if multiple negative sources exist
            if not energy_long_neg.empty:
                 fig_energy.add_bar(
                     x=energy_long_neg['Zone'],
                     y=energy_long_neg['Energy_TWh'],
                     name='Curtailment',
                     marker_color=sanitized_colors.get('Curtailment', 'red') # Use sanitized map, ensure key matches case in data ('Curtailment')
                 )

            # Add Demand Lines
            zone_order_energy = sorted(selected_zones)
            fig_energy.update_xaxes(categoryorder='array', categoryarray=zone_order_energy)

            shapes_energy = []
            for i, zone in enumerate(zone_order_energy):
                # Black line for DA Demand
                shapes_energy.append(dict(type='line', x0=i-0.4, x1=i+0.4, y0=demand_da[zone], y1=demand_da[zone], xref='x', yref='y', line=dict(color='Black', width=2, dash='solid')))
                # Green line for DA + Storage Input Demand
                shapes_energy.append(dict(type='line', x0=i-0.4, x1=i+0.4, y0=total_demand_line[zone], y1=total_demand_line[zone], xref='x', yref='y', line=dict(color='Green', width=2, dash='dash')))
                # Blue line for DA + Storage + P2X Demand
                shapes_energy.append(dict(type='line', x0=i-0.4, x1=i+0.4, y0=storage_demand_line[zone], y1=storage_demand_line[zone], xref='x', yref='y', line=dict(color='Blue', width=2, dash='dot')))
                # Red line for Final Demand (incl ShedLoad)
                # shapes_energy.append(dict(type='line', x0=i-0.4, x1=i+0.4, y0=final_demand_line[zone], y1=final_demand_line[zone], xref='x', yref='y', line=dict(color='Red', width=2, dash='longdash')))
                # Note: ShedLoad is hard to represent with lines on a stacked bar. Consider adding as text annotation or separate plot if needed.

            fig_energy.update_layout(shapes=shapes_energy, barmode='relative') # Use relative barmode for +/- stacking
            fig_energy.update_layout(legend_title_text='Source/Fuel')

            st.plotly_chart(fig_energy, use_container_width=True)

        except Exception as e:
            st.error(f"Error generating energy balance plot: {e}")
            # Optionally print traceback for debugging
            # import traceback
            # st.error(traceback.format_exc())
            plt.close() # Ensure figure is closed on error
    else:
        st.info("Select at least one zone to display the plot.")

    # --- Phase 2: Step 4 - Detailed Zone Dispatch Plot ---
    st.header("Detailed Zone Analysis")
    zone_expander = st.expander("Configure and View Zone Dispatch", expanded=False)
    with zone_expander:
        selected_zone_detail = st.selectbox(
            "Select Zone for Detailed View",
            options=all_zones,
            index=0 # Default to the first zone
        )

        # Get the actual date range from results if available
        # Use OutputPower as a reference, fallback if not present
        if 'OutputPower' in results and not results['OutputPower'].empty:
            min_date = results['OutputPower'].index.min().date()
            max_date = results['OutputPower'].index.max().date()
        else:
            # Fallback dates if results are empty or OutputPower is missing
            min_date = datetime.date(2023, 1, 1)
            max_date = datetime.date(2023, 12, 31)

        # Date Range Selection
        col1, col2 = st.columns(2)
        with col1:
            date_range_start = st.date_input("Start Date", value=min_date, min_value=min_date, max_value=max_date)
        with col2:
            # Ensure end date is not before start date
            valid_end_date = max(date_range_start, min_date) # End date cannot be before start date
            default_end_date = min(date_range_start + datetime.timedelta(days=6), max_date) # Default to 7 days or max_date
            if default_end_date < valid_end_date:
                default_end_date = valid_end_date # Adjust if default range is invalid

            date_range_end = st.date_input("End Date", value=default_end_date, min_value=valid_end_date, max_value=max_date)

        # Combine dates into a DatetimeIndex range (inclusive end needed for pd.date_range)
        # Add time components to span full days
        start_datetime = datetime.datetime.combine(date_range_start, datetime.time.min)
        end_datetime = datetime.datetime.combine(date_range_end, datetime.time.max)

        # Create the range - assumes hourly frequency from Dispa-SET
        try:
            rng = pd.date_range(start=start_datetime, end=end_datetime, freq='h')
            # Ensure the generated range is within the simulation's actual data range
            if 'OutputPower' in results and not results['OutputPower'].empty:
                sim_index = results['OutputPower'].index
                rng = sim_index.intersection(rng)
        except Exception as e:
             st.error(f"Error creating date range: {e}")
             rng = None

        # Plotting Options
        hide_storage = st.checkbox("Hide Storage Level Plot", value=False)
        show_rug = st.checkbox("Show Rug Plot (Generation)", value=True)

        if rng is not None and not rng.empty:
            st.subheader(f"Dispatch Plot for Zone: {selected_zone_detail}")
            try:
                # Create a figure context for plot_zone
                plt.figure(figsize=(13, 7))
                # ds.plot_zone draws on the current figure, it doesn't accept show_plot
                plot_success = ds.plot_zone(
                    inputs, results,
                    z=selected_zone_detail,
                    rng=rng,
                    rug_plot=False, # Control rug plot separately
                    hide_storage_plot=hide_storage
                )
                if plot_success:
                    fig_dispatch = plt.gcf() # Get the current figure created by plot_zone
                    st.pyplot(fig_dispatch)
                    plt.close(fig_dispatch) # Close the figure to free memory
                else:
                    st.warning("Could not generate dispatch plot.")

            except Exception as e:
                st.error(f"Error generating dispatch plot for zone {selected_zone_detail}: {e}")
                # import traceback
                # st.error(traceback.format_exc())
                plt.close() # Ensure figure is closed on error

            # Rug Plot (if enabled and possible)
            if show_rug:
                st.subheader(f"Generation Rug Plot for Zone: {selected_zone_detail}")
                try:
                    zone_generation = filter_by_zone(results['OutputPower'], inputs, selected_zone_detail)
                    zone_generation_ranged = zone_generation.loc[rng]

                    if not zone_generation_ranged.empty and not zone_generation_ranged.isnull().all().all():
                         # Use the plot_rug function from the library
                         plt.figure(figsize=(16, max(1, 0.25 * len(zone_generation_ranged.columns)))) # Adjust figsize
                         plot_rug(zone_generation_ranged, on_off=False, cmap='gist_heat_r', fig_title=selected_zone_detail)
                         fig_rug = plt.gcf()
                         st.pyplot(fig_rug)
                         plt.close(fig_rug) # Close the figure
                    else:
                        st.info(f"No generation data to display in rug plot for zone {selected_zone_detail} in the selected range.")
                except ImportError:
                     st.warning("Could not import `enlopy` for rug plot. Skipping rug plot.") # In case internal plot_rug uses it
                except KeyError as e:
                     st.warning(f"Could not find necessary data for rug plot (KeyError: {e}). Skipping rug plot.")
                except Exception as e:
                    st.error(f"Error generating rug plot: {e}")
                    # import traceback
                    # st.error(traceback.format_exc())
                    plt.close() # Ensure figure is closed on error
        elif rng is not None and rng.empty:
             st.warning("Selected date range has no overlap with simulation data. Please adjust the dates.")
        else:
            st.warning("Please select a valid date range.")

    # --- Phase 2: Step 5 - Boundary Sector Dispatch Plot ---
    st.header("Boundary Sector Analysis")
    boundary_expander = st.expander("View Boundary Sector Dispatch", expanded=False)
    with boundary_expander:
        # Try to identify boundary sector zones from results keys
        boundary_zones = []
        if 'OutputPowerX' in results and not results['OutputPowerX'].empty:
            if isinstance(results['OutputPowerX'].columns, pd.MultiIndex): # Handle multi-index if present
                boundary_zones = sorted(list(results['OutputPowerX'].columns.get_level_values('Zones').unique()))
            else:
                 # Assuming format like 'Z1_th_1 - UnitName', extract the zone part
                 # This might need adjustment based on actual column format
                 try:
                      boundary_zones = sorted(list(set([col.split(' - ')[0] for col in results['OutputPowerX'].columns])))
                 except Exception:
                     # Fallback if splitting fails
                     boundary_zones = sorted(list(results['OutputPowerX'].columns))
        elif 'OutputSectorXStorageLevel' in results and not results['OutputSectorXStorageLevel'].empty:
             boundary_zones = sorted(list(results['OutputSectorXStorageLevel'].columns))
        # Add other potential sources like OutputSectorXFlow if needed

        if not boundary_zones:
            st.info("No boundary sector data found in results (checked OutputPowerX, OutputSectorXStorageLevel).")
        else:
            selected_boundary_zone = st.selectbox(
                "Select Boundary Sector Zone",
                options=boundary_zones,
                index=0
            )

            # Reuse date range controls (or create new ones if different range needed)
            # Using the same date range logic as Detailed Zone Analysis
            if 'rng' in locals() and rng is not None and not rng.empty:
                st.subheader(f"Boundary Dispatch Plot for: {selected_boundary_zone}")
                try:
                    plt.figure(figsize=(13, 7)) # Create figure context
                    plot_success_bx = plot_dispatchX(
                        inputs, results,
                        z=selected_boundary_zone,
                        rng=rng
                        # Add other plot_dispatchX args like alpha if needed
                        # show_plot=False # plot_dispatchX likely doesn't have this either
                    )
                    if plot_success_bx:
                         fig_bx = plt.gcf()
                         st.pyplot(fig_bx)
                         plt.close(fig_bx)
                    else:
                        st.warning(f"Could not generate boundary dispatch plot for {selected_boundary_zone}.")
                except Exception as e:
                    st.error(f"Error generating boundary dispatch plot for {selected_boundary_zone}: {e}")
                    # import traceback
                    # st.error(traceback.format_exc())
                    plt.close() # Ensure figure closed on error
            else:
                st.warning("Please select a valid date range in the 'Detailed Zone Analysis' section first.")

    # --- Phase 2: Step 6 - Zone Capacities Plot (using Plotly) ---
    st.header("Installed Capacity Analysis")
    capacity_expander = st.expander("View Installed Capacity by Zone and Fuel", expanded=False)
    with capacity_expander:
        try:
            # --- Data Preparation (adapted from ds.plot_zone_capacities) ---
            units_df = inputs['units'].copy()
            valid_technologies = list(commons['Technologies'])
            power_units = units_df[units_df['Technology'].isin(valid_technologies)]
            # Add power consumption units if they exist and have capacity defined
            # (e.g., for P2X, though capacity might not be the right metric here - adjust if needed)
            if 'OutputPowerConsumption' in results:
                 cons_units = units_df[units_df.index.isin(results['OutputPowerConsumption'].columns)]
                 power_units = pd.concat([power_units, cons_units[~cons_units.index.isin(power_units.index)]])

            # Handle CHP units specifically if needed (original function adjusts capacity)
            # For simplicity here, we use the base PowerCapacity. Add CHP logic if complex representation is required.

            # Aggregate capacity
            power_units['TotalCapacity'] = power_units['PowerCapacity'] * power_units['Nunits']
            capacity_by_zone_fuel = power_units.groupby(['Zone', 'Fuel'])['TotalCapacity'].sum().unstack(fill_value=0) / 1000 # GW

            # Ensure all common fuels are present as columns, fill missing with 0
            for fuel in commons['MeritOrder']:
                if fuel not in capacity_by_zone_fuel.columns:
                    capacity_by_zone_fuel[fuel] = 0

            # Select and order columns based on MeritOrder
            plot_cols = [col for col in commons['MeritOrder'] if col in capacity_by_zone_fuel.columns]
            capacity_by_zone_fuel = capacity_by_zone_fuel[plot_cols]

            # Filter out zones with zero total capacity
            capacity_by_zone_fuel = capacity_by_zone_fuel.loc[capacity_by_zone_fuel.sum(axis=1) > 0]

            # Reshape for Plotly Express
            capacity_long = capacity_by_zone_fuel.stack().reset_index()
            capacity_long.columns = ['Zone', 'Fuel', 'Capacity_GW']
            capacity_long = capacity_long[capacity_long['Capacity_GW'] > 0] # Remove zero capacity entries for cleaner plot

            # Calculate Peak Demand for reference lines
            peak_demand_da = inputs['param_df']['Demand']['DA'].max(axis=0) / 1000 # GW
            peak_demand_total = peak_demand_da.copy()
            # Add peak power consumption if relevant
            if 'OutputPowerConsumption' in results:
                peak_consumption = results['OutputPowerConsumption'].max(axis=0) / 1000 #GW
                # Align indices and add, filling missing zones with 0
                peak_demand_total = peak_demand_total.add(peak_consumption, fill_value=0)

            # Filter demand data for zones present in capacity plot
            peak_demand_da = peak_demand_da.reindex(capacity_by_zone_fuel.index).fillna(0)
            peak_demand_total = peak_demand_total.reindex(capacity_by_zone_fuel.index).fillna(0)

            # --- Plotting with Plotly Express ---
            # Sanitize colors for Plotly (remove alpha channel if present, keep valid formats)
            sanitized_colors = {}
            for fuel, color in commons['colors'].items():
                if isinstance(color, str):
                    if color.startswith('#') and len(color) == 9: # Check for #RRGGBBAA
                        sanitized_colors[fuel] = color[:7] # Truncate to #RRGGBB
                    elif color.startswith('#') and len(color) == 7: # Check for #RRGGBB
                        sanitized_colors[fuel] = color
                    # Add checks for other valid formats if needed, e.g., CSS names
                    # For now, only handle hex codes explicitly
            # You might want to add a default color for fuels not in sanitized_colors

            if not capacity_long.empty:
                fig_cap = px.bar(
                    capacity_long,
                    x='Zone',\
                    y='Capacity_GW',\
                    color='Fuel',\
                    title='Installed Capacity per Zone (GW)',\
                    labels={'Capacity_GW': 'Installed Capacity [GW]'},\
                    color_discrete_map=sanitized_colors # Use sanitized colors
                )

                # Add Peak Demand Lines (add shapes for horizontal lines)
                # Plotly needs categoryarray for correct shape positioning on categorical axes
                zone_order = capacity_by_zone_fuel.index.tolist()
                fig_cap.update_xaxes(categoryorder='array', categoryarray=zone_order)

                shapes = []
                for i, zone in enumerate(zone_order):
                    # Black line for Peak DA Demand
                    shapes.append(dict(type='line', x0=i-0.4, x1=i+0.4, y0=peak_demand_da[zone], y1=peak_demand_da[zone], xref='x', yref='y', line=dict(color='Black', width=2)))
                    # Blue line for Peak Total Demand (DA + Consumption)
                    shapes.append(dict(type='line', x0=i-0.4, x1=i+0.4, y0=peak_demand_total[zone], y1=peak_demand_total[zone], xref='x', yref='y', line=dict(color='Blue', width=2)))
                fig_cap.update_layout(shapes=shapes)

                # Improve layout
                fig_cap.update_layout(legend_title_text='Fuel Type')

                st.plotly_chart(fig_cap, use_container_width=True)
            else:
                st.info("No capacity data available to plot.")

        except Exception as e:
            st.error(f"Error generating capacity plot: {e}")
            # import traceback
            # st.error(traceback.format_exc())

elif simulation_path:
    st.info("Click 'Load Data' in the sidebar to begin.")
else:
    st.info("Enter a simulation path in the sidebar and click 'Load Data'.") 