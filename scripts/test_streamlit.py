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
from dispaset.postprocessing.postprocessing import get_result_analysis # Added result analysis
from dispaset.postprocessing.plot import plot_rug, plot_dispatchX, plot_power_flow_tracing_matrix # Added power flow
# Added map plot imports
try:
    from dispaset.postprocessing.geoplot import plot_net_flows_map, plot_line_congestion_map
    GEOPLOT_AVAILABLE = True
except ImportError:
    GEOPLOT_AVAILABLE = False
    # Define dummy functions if import fails to avoid NameError later
    def plot_net_flows_map(*args, **kwargs):
        raise ImportError("Cartopy or other dependency missing for geoplot.")
    def plot_line_congestion_map(*args, **kwargs):
        raise ImportError("Cartopy or other dependency missing for geoplot.")

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

# --- Helper Functions ---
@st.cache_data # Cache the conversion to avoid recomputing
def convert_df_to_csv(df):
    """Converts a DataFrame to CSV bytes for downloading."""
    return df.to_csv(index=True).encode('utf-8')

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

    # --- Create Tabs for Organization ---
    tab_overview, tab_zone, tab_network, tab_unit = st.tabs([
        "Overview",
        "Zone Analysis",
        "Network Analysis",
        "Unit Analysis"
    ])

    # --- Tab 1: Overview ---
    with tab_overview:
        st.header("Simulation Overview")

        # Energy Balance Plot (Plotly)
        st.subheader("Energy Balance per Zone")
        # (Keep the multiselect for zones here or make it global? Keep here for now)
        selected_zones_overview = st.multiselect(
            "Select Zones for Energy Balance & Capacity",
            options=sorted(inputs['sets']['n']),
            default=sorted(inputs['sets']['n']),
            key="overview_zone_select" # Unique key
        )

        if selected_zones_overview:
            # (Existing code for Energy Balance Plotly - needs selected_zones_overview)
            # --- Data Preparation (adapted from ds.plot_energy_zone_fuel) ---
            zones_in_results = indicators.Zone.unique()
            fuels = indicators[[u in commons['Technologies'] for u in indicators['Technology']]].Fuel.unique()
            GenPerZone = pd.DataFrame(0.0, index=zones_in_results, columns=list(fuels) + ['FlowIn', 'Curtailment'])

            # Ensure all common fuels are present
            for f in commons['Fuels'] + ['FlowIn']:
                if f not in GenPerZone.columns:
                    GenPerZone[f] = 0.0

            # Aggregate generation and calculate imports/demand/storage input
            power_consumption = pd.DataFrame(0.0, index=results['OutputPowerConsumption'].index, columns=zones_in_results)
            storage_input = pd.DataFrame(0.0, index=results['OutputPower'].index, columns=zones_in_results)

            for z in zones_in_results:
                if z in selected_zones_overview:  # Only process selected zones
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
            GenPerZone = GenPerZone.loc[selected_zones_overview]

            # Select and order columns based on MeritOrder (plus Curtailment/ShedLoad)
            plot_cols_gen = [col for col in commons['MeritOrder'] if col in GenPerZone.columns and GenPerZone[col].abs().sum() != 0]
            plot_cols_other = ['Curtailment'] # Removed ShedLoad plotting for simplicity
            GenPerZone_plot = GenPerZone[[c for c in plot_cols_gen + plot_cols_other if c in GenPerZone.columns]].copy()

            # Calculate Demand components (TWh)
            demand_da = inputs['param_df']['Demand']['DA'].sum(axis=0).reindex(selected_zones_overview).fillna(0) / 1E6
            p_cons_total = power_consumption.sum(axis=0).reindex(selected_zones_overview).fillna(0) / 1E6
            s_in_total = storage_input.sum(axis=0).reindex(selected_zones_overview).fillna(0) / 1E6

            total_demand_line = demand_da + p_cons_total # DA + P2X
            storage_demand_line = total_demand_line + s_in_total # DA + P2X + Storage Input

            # Reshape for Plotly (Handle positive and negative separately for stacking)
            GenPerZone_plot['TotalGen'] = GenPerZone_plot[plot_cols_gen].clip(lower=0).sum(axis=1) # Sum positive generation

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
                color_discrete_map=sanitized_colors, # Use locally sanitized colors
                category_orders={"Zone": sorted(selected_zones_overview)}
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
            zone_order_energy = sorted(selected_zones_overview)
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

        else:
            st.info("Select at least one zone to display the Energy Balance plot.")

        # Capacity Analysis Plot (Plotly)
        st.subheader("Installed Capacity per Zone")
        if selected_zones_overview:
            # (Existing code for Capacity Plotly using selected_zones_overview)
            # --- Data Preparation (adapted from ds.plot_zone_capacities) ---
            units_df = inputs['units'].copy()
            valid_technologies = list(commons['Technologies'])
            power_units = units_df[units_df['Technology'].isin(valid_technologies)]
            # Add power consumption units if they exist and have capacity defined
            if 'OutputPowerConsumption' in results:
                 cons_units = units_df[units_df.index.isin(results['OutputPowerConsumption'].columns)]
                 power_units = pd.concat([power_units, cons_units[~cons_units.index.isin(power_units.index)]])

            # Handle CHP units specifically if needed (original function adjusts capacity)
            # For simplicity here, we use the base PowerCapacity. Add CHP logic if complex representation is required.

            # Aggregate capacity
            power_units['TotalCapacity'] = power_units['PowerCapacity'] * power_units['Nunits']
            capacity_by_zone_fuel = power_units.groupby(['Zone', 'Fuel'])['TotalCapacity'].sum().unstack(fill_value=0) / 1000 # GW

            # Ensure all common fuels are present as columns, fill missing with 0
            # Use fuels present in the data for color mapping consistency
            present_fuels = capacity_by_zone_fuel.columns.tolist()
            for fuel in commons['MeritOrder']:
                if fuel not in capacity_by_zone_fuel.columns:
                    capacity_by_zone_fuel[fuel] = 0

            # Select and order columns based on MeritOrder
            plot_cols = [col for col in commons['MeritOrder'] if col in capacity_by_zone_fuel.columns]
            capacity_by_zone_fuel = capacity_by_zone_fuel[plot_cols]

            # Filter out zones with zero total capacity
            capacity_by_zone_fuel = capacity_by_zone_fuel.loc[capacity_by_zone_fuel.sum(axis=1) > 0]
            # Filter for selected zones
            capacity_by_zone_fuel = capacity_by_zone_fuel.reindex(selected_zones_overview).dropna(how='all')

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
            # Re-sanitize colors here as well
            sanitized_colors = {}
            for fuel, color in commons['colors'].items():
                if isinstance(color, str):
                    if color.startswith('#') and len(color) == 9:
                        sanitized_colors[fuel] = color[:7]
                    elif color.startswith('#') and len(color) == 7:
                        sanitized_colors[fuel] = color

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
                zone_order_cap = capacity_by_zone_fuel.index.tolist()
                fig_cap.update_xaxes(categoryorder='array', categoryarray=zone_order_cap)

                shapes = []
                for i, zone in enumerate(zone_order_cap):
                    # Black line for Peak DA Demand
                    shapes.append(dict(type='line', x0=i-0.4, x1=i+0.4, y0=peak_demand_da[zone], y1=peak_demand_da[zone], xref='x', yref='y', line=dict(color='Black', width=2)))
                    # Blue line for Peak Total Demand (DA + Consumption)
                    shapes.append(dict(type='line', x0=i-0.4, x1=i+0.4, y0=peak_demand_total[zone], y1=peak_demand_total[zone], xref='x', yref='y', line=dict(color='Blue', width=2)))
                fig_cap.update_layout(shapes=shapes)

                # Improve layout
                fig_cap.update_layout(legend_title_text='Fuel Type')

                st.plotly_chart(fig_cap, use_container_width=True)
            else:
                st.info("No capacity data available to plot for selected zones.")

        else:
            st.info("Select zones to display Capacity plot.")

        # Simulation Summary Metrics and Tables
        st.subheader("Simulation Summary")
        analysis_results = get_result_analysis(inputs, results) # Get analysis results (cached)

        if analysis_results:
            st.markdown("**Aggregated Metrics** provide a high-level overview of the entire simulated system's performance, including total energy demand, peak loads, costs, and indicators of grid stress like load shedding (unmet demand) and curtailment (excess renewable generation). Net imports show the reliance on external zones.")
            # Display key scalar metrics using columns and st.metric
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                st.metric("Total Load (TWh)", f"{analysis_results.get('TotalLoad', 0) / 1e6:.2f}")
                st.metric("Total Curtailment (TWh)", f"{analysis_results.get('Curtailment', 0) / 1e6:.2f}")
            with col2:
                st.metric("Peak Load (GW)", f"{analysis_results.get('PeakLoad', 0) / 1e3:.2f}")
                st.metric("Max Curtailment (GW)", f"{analysis_results.get('MaxCurtailment', 0) / 1e3:.2f}") # Assuming MaxCurtailment is in MW
            with col3:
                st.metric("Net Imports (TWh)", f"{analysis_results.get('NetImports', 0) / 1e6:.2f}")
                st.metric("Total Shed Load (TWh)", f"{analysis_results.get('ShedLoad', 0):.2f}") # Assuming ShedLoad is already TWh from function
            with col4:
                st.metric("Avg Cost (€/MWh)", f"{analysis_results.get('Cost_kwh', 'N/A'):.2f}" if isinstance(analysis_results.get('Cost_kwh'), (int, float)) else "N/A")
                st.metric("Max Shed Load (GW)", f"{analysis_results.get('MaxShedLoad', 0) / 1e3:.2f}") # Assuming MaxShedLoad is in MW

            # Add Shifted Load if it exists
            if 'ShiftedLoad' in analysis_results:
                st.metric("Total Shifted Load (TWh)", f"{analysis_results.get('ShiftedLoad', 0):.2f}") # Assuming TWh

            st.markdown("""
            **Aggregated Metrics** provide a high-level overview of the entire simulated system's performance, including total energy demand, peak loads, costs, and indicators of grid stress like load shedding (unmet demand) and curtailment (excess renewable generation). Net imports show the reliance on external zones.
            """)

            # Summary Tables Expander
            with st.expander("View Detailed Summary Tables"):
                st.markdown("**Zone Specific Data** breaks down key metrics by geographical area, allowing comparison of performance and stress levels across different regions within the simulation.")
                if analysis_results.get('ZoneData') is not None:
                    df_zone_data = analysis_results['ZoneData']
                    st.dataframe(df_zone_data)
                    csv_zone_data = convert_df_to_csv(df_zone_data)
                    st.download_button(
                        label="Download Zone Data as CSV",
                        data=csv_zone_data,
                        file_name='zone_data.csv',
                        mime='text/csv',
                    )
                else:
                    st.info("No Zone Specific Data available.")

                st.markdown("**Storage Data** provides insights into the installed energy storage capacity (MWh) and power rating (MW) per zone. It can also include metrics like average cycle depth if calculated, indicating storage utilization.")
                if analysis_results.get('StorageData') is not None:
                    df_storage_data = analysis_results['StorageData'].dropna(axis=1, how='all')
                    st.dataframe(df_storage_data)
                    csv_storage_data = convert_df_to_csv(df_storage_data)
                    st.download_button(
                        label="Download Storage Data as CSV",
                        data=csv_storage_data,
                        file_name='storage_data.csv',
                        mime='text/csv',
                    )
                else:
                    st.info("No Storage Data available.")

                st.markdown("**Congestion Data** reveals bottlenecks in the transmission network by showing the number of hours each interconnection line operated at its maximum power transfer capacity.")
                if analysis_results.get('Congestion'):
                    st.write("Number of hours of congestion per line:")
                    st.json(analysis_results['Congestion'])
                else:
                    st.info("No Congestion Data available.")

                st.markdown("**Water Consumption Data** summarizes the total water withdrawn and consumed by power plants within each zone, crucial for assessing the environmental footprint, especially in water-scarce regions.")
                if analysis_results.get('WaterConsumptionData'):
                    water_data = analysis_results['WaterConsumptionData']
                    if water_data.get('ZoneLevel'):
                        df_water_w = pd.DataFrame(water_data['ZoneLevel']['WaterWithdrawal'].sum(), columns=['TotalWithdrawal_m3'])
                        df_water_c = pd.DataFrame(water_data['ZoneLevel']['WaterConsumption'].sum(), columns=['TotalConsumption_m3'])
                        st.subheader("Zone Level Water Withdrawal (Total m³ over period)")
                        st.dataframe(df_water_w)
                        st.download_button("Download Water Withdrawal Data", convert_df_to_csv(df_water_w), "water_withdrawal.csv", "text/csv")
                        st.subheader("Zone Level Water Consumption (Total m³ over period)")
                        st.dataframe(df_water_c)
                        st.download_button("Download Water Consumption Data", convert_df_to_csv(df_water_c), "water_consumption.csv", "text/csv")
                else:
                    st.info("No Water Consumption Data available.")
        else:
            st.warning("Could not retrieve simulation analysis results. Check logs for errors.")

    # --- Tab 2: Zone Analysis ---
    with tab_zone:
        st.header("Detailed Zone Analysis")
        all_zones_list = sorted(inputs['sets']['n']) # Get zones again for this tab
        selected_zone_detail = st.selectbox(
            "Select Zone for Detailed View",
            options=all_zones_list,
            index=0, # Default to the first zone
            key="zone_detail_select" # Unique key
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

        # Date Range Selection for Zone Analysis
        col1_z, col2_z = st.columns(2)
        with col1_z:
            date_range_start_z = st.date_input("Start Date", value=min_date, min_value=min_date, max_value=max_date, key="zone_date_start")
        with col2_z:
            # Ensure end date is not before start date
            valid_end_date_z = max(date_range_start_z, min_date) # End date cannot be before start date
            default_end_date_z = min(date_range_start_z + datetime.timedelta(days=6), max_date) # Default to 7 days or max_date
            if default_end_date_z < valid_end_date_z:
                default_end_date_z = valid_end_date_z # Adjust if default range is invalid

            date_range_end_z = st.date_input("End Date", value=default_end_date_z, min_value=valid_end_date_z, max_value=max_date, key="zone_date_end")

        # Combine dates into a DatetimeIndex range (inclusive end needed for pd.date_range)
        # Add time components to span full days
        start_datetime_z = datetime.datetime.combine(date_range_start_z, datetime.time.min)
        end_datetime_z = datetime.datetime.combine(date_range_end_z, datetime.time.max)

        # Create the range - assumes hourly frequency from Dispa-SET
        try:
            rng_zone = pd.date_range(start=start_datetime_z, end=end_datetime_z, freq='h')
            # Ensure the generated range is within the simulation's actual data range
            if 'OutputPower' in results and not results['OutputPower'].empty:
                sim_index = results['OutputPower'].index
                rng_zone = sim_index.intersection(rng_zone)
        except Exception as e:
             st.error(f"Error creating zone date range: {e}")
             rng_zone = None

        # Plotting Options
        hide_storage_z = st.checkbox("Hide Storage Level Plot", value=False, key="zone_hide_storage")
        show_rug_z = st.checkbox("Show Rug Plot (Generation)", value=True, key="zone_show_rug")

        if rng_zone is not None and not rng_zone.empty:
            st.subheader(f"Dispatch Plot for Zone: {selected_zone_detail}")
            try:
                # Create a figure context for plot_zone
                plt.figure(figsize=(13, 7))
                # ds.plot_zone draws on the current figure, it doesn't accept show_plot
                plot_success = ds.plot_zone(
                    inputs, results,
                    z=selected_zone_detail,
                    rng=rng_zone,
                    rug_plot=False, # Control rug plot separately
                    hide_storage_plot=hide_storage_z
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
            if show_rug_z:
                st.subheader(f"Generation Rug Plot for Zone: {selected_zone_detail}")
                try:
                    zone_generation = filter_by_zone(results['OutputPower'], inputs, selected_zone_detail)
                    zone_generation_ranged = zone_generation.loc[rng_zone]

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
        elif rng_zone is not None and rng_zone.empty:
             st.warning("Selected date range has no overlap with simulation data. Please adjust the dates.")
        else:
            st.warning("Please select a valid date range.")

    # --- Tab 3: Network Analysis ---
    with tab_network:
        st.header("Network Analysis")

        # Map Plots
        st.subheader("Network Map Visualizations")
        map_expander = st.expander("View Network Flows and Congestion Maps", expanded=False)
        with map_expander:
            if not GEOPLOT_AVAILABLE:
                 st.error("Module 'cartopy' or other dependency for geoplotting not found. Map plots are unavailable.")
            elif rng_zone is None or rng_zone.empty: # Use range from Zone Analysis tab
                 st.warning("Please select a valid date range in the 'Zone Analysis' tab first to view maps.")
            else:
                 # Net Flows Map
                 st.subheader("Net Flows Map")
                 try:
                     plt.figure(figsize=(12, 8))
                     plot_net_flows_map(inputs, results, idx=rng_zone)
                     fig_net_flows = plt.gcf()
                     st.pyplot(fig_net_flows)
                     plt.close(fig_net_flows)
                 except ImportError as e:
                     st.error(f"Could not generate Net Flows Map. Missing dependency: {e}")
                 except Exception as e:
                     st.error(f"Error generating Net Flows Map: {e}")
                     plt.close()

                 # Line Congestion Map
                 st.subheader("Line Congestion Map")
                 try:
                     plt.figure(figsize=(12, 8))
                     plot_line_congestion_map(inputs, results, idx=rng_zone)
                     fig_congestion = plt.gcf()
                     st.pyplot(fig_congestion)
                     plt.close(fig_congestion)
                 except ImportError as e:
                      st.error(f"Could not generate Line Congestion Map. Missing dependency: {e}")
                 except Exception as e:
                     st.error(f"Error generating Line Congestion Map: {e}")
                     plt.close()

        # Power Flow Tracing
        st.subheader("Power Flow Tracing")
        pft_expander = st.expander("View Power Flow Tracing Matrix", expanded=False)
        with pft_expander:
            if rng_zone is None or rng_zone.empty: # Use range from Zone Analysis tab
                 st.warning("Please select a valid date range in the 'Zone Analysis' tab first.")
            else:
                 # This function likely creates its own figure
                 # Call the function, it should plot to the current figure implicitly
                 data_pft, data_pft_prct = plot_power_flow_tracing_matrix(
                     inputs, results,
                     idx=rng_zone,
                     figsize=(10, 7) # Pass figsize or other args if needed
                 )
                 fig_pft = plt.gcf() # Get the figure created by the function
                 st.pyplot(fig_pft)
                 plt.close(fig_pft) # Close the figure

                 # Optionally display the dataframes
                 show_pft_data = st.checkbox("Show Power Flow Tracing DataFrames", value=False)
                 if show_pft_data:
                     st.subheader("Power Flow Tracing Data (TWh over period)")
                     st.dataframe(data_pft / 1e6) # Assuming MW -> TWh conversion needed
                     st.subheader("Power Flow Tracing Data (% Share)")
                     st.dataframe(data_pft_prct)

                 # Add download buttons for PFT data
                 if 'data_pft' in locals() and 'data_pft_prct' in locals():
                      csv_pft = convert_df_to_csv(data_pft / 1e6) # Convert to TWh for download
                      st.download_button("Download PFT Data (TWh)", csv_pft, "pft_twh.csv", "text/csv")
                      csv_pft_prct = convert_df_to_csv(data_pft_prct)
                      st.download_button("Download PFT Data (% Share)", csv_pft_prct, "pft_percent.csv", "text/csv")

    # --- Tab 4: Unit Analysis ---
    with tab_unit:
        st.header("Unit Analysis")

        # --- Step 11: Filtering ---
        unit_data_available = analysis_results is not None and analysis_results.get('UnitData') is not None
        indicators_available = indicators is not None

        if unit_data_available or indicators_available:
            st.subheader("Filter Unit Data")
            # Get filter options (unique values from relevant columns)
            zone_options = sorted(list(inputs['sets']['n']))
            tech_options = sorted(list(units_df['Technology'].unique()))
            fuel_options = sorted(list(units_df['Fuel'].unique()))

            filt_col1, filt_col2, filt_col3 = st.columns(3)
            with filt_col1:
                sel_zones = st.multiselect("Filter by Zone", zone_options, default=zone_options)
            with filt_col2:
                sel_techs = st.multiselect("Filter by Technology", tech_options, default=tech_options)
            with filt_col3:
                sel_fuels = st.multiselect("Filter by Fuel", fuel_options, default=fuel_options)

            # Apply filters
            filtered_indicators = indicators.copy() if indicators_available else pd.DataFrame()
            filtered_unitdata = analysis_results['UnitData'].copy() if unit_data_available else pd.DataFrame()

            if indicators_available:
                filtered_indicators = filtered_indicators[
                    filtered_indicators['Zone'].isin(sel_zones) &
                    filtered_indicators['Technology'].isin(sel_techs) &
                    filtered_indicators['Fuel'].isin(sel_fuels)
                ]

            if unit_data_available:
                 filtered_unitdata = filtered_unitdata[
                     filtered_unitdata['Zone'].isin(sel_zones) &
                     filtered_unitdata['Technology'].isin(sel_techs) &
                     filtered_unitdata['Fuel'].isin(sel_fuels)
                 ]
        else:
            st.info("Load data to enable unit filtering.")
            filtered_indicators = pd.DataFrame()
            filtered_unitdata = pd.DataFrame()

        # Power Plant Indicators Table
        st.subheader("Power Plant / Unit Indicators")
        st.markdown("Detailed operational metrics per unit (e.g., generation, CF, startups). Use filters above.")
        if indicators_available:
            st.dataframe(filtered_indicators)
            # --- Step 14: Download Button ---
            csv_indicators = convert_df_to_csv(filtered_indicators)
            st.download_button(
                label="Download Filtered Indicators as CSV",
                data=csv_indicators,
                file_name='filtered_indicators.csv',
                mime='text/csv',
            )
        else:
            st.warning("Power plant indicators data is not available.")

        # Unit Specific Data Table (from Summary Analysis)
        st.subheader("Unit Specific Data (Costs, Emissions, Water)")
        st.markdown("Detailed cost, emission, and water usage data per unit. Use filters above.")
        if unit_data_available:
             st.dataframe(filtered_unitdata)
             # --- Step 14: Download Button ---
             csv_unitdata = convert_df_to_csv(filtered_unitdata)
             st.download_button(
                 label="Download Filtered Unit Data as CSV",
                 data=csv_unitdata,
                 file_name='filtered_unit_data.csv',
                 mime='text/csv',
             )
        else:
             st.info("Unit Specific Data not available in summary results.")

elif simulation_path:
    st.info("Click 'Load Data' in the sidebar to begin.")
else:
    st.info("Enter a simulation path in the sidebar and click 'Load Data'.") 