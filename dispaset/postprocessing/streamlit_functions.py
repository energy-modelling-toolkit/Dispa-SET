import streamlit as st
import pandas as pd
import plotly.express as px
import matplotlib.pyplot as plt
import datetime
import logging # Added for logging potential issues

# Import necessary dispaset functions (assuming this file is in postprocessing)
from .postprocessing import (
    filter_by_zone, get_imports, filter_by_tech_list,
    get_result_analysis, filter_by_storage,
    filter_sector # Added filter_sector potentially needed
)
from .plot import (
    plot_rug, plot_dispatchX, plot_power_flow_tracing_matrix,
    plot_energy_zone_fuel, plot_zone # Keep for reference/data prep logic if needed
)
try:
    from .geoplot import plot_net_flows_map, plot_line_congestion_map
    GEOPLOT_AVAILABLE = True
except ImportError:
    logging.warning("Could not import geoplot functions. Map plotting will be disabled. Ensure cartopy and dependencies are installed.")
    GEOPLOT_AVAILABLE = False
    # Define dummy functions if import fails
    def plot_net_flows_map(*args, **kwargs): st.error("Map plotting unavailable (missing dependencies).")
    def plot_line_congestion_map(*args, **kwargs): st.error("Map plotting unavailable (missing dependencies).")

# Import commons
from ..common import commons

# --- Helper Functions ---
@st.cache_data # Cache the conversion to avoid recomputing
def convert_df_to_csv(df):
    """Converts a DataFrame to CSV bytes for downloading."""
    if isinstance(df, pd.DataFrame):
        return df.to_csv(index=True).encode('utf-8')
    elif isinstance(df, pd.Series):
        return df.to_csv(index=True).encode('utf-8')
    return b"" # Return empty bytes if not DataFrame or Series

def _sanitize_colors(color_dict):
    """Removes alpha channel from hex codes for Plotly compatibility."""
    sanitized = {}
    default_colors = px.colors.qualitative.Plotly
    color_index = 0
    for key, color in color_dict.items():
        valid_color = None
        if isinstance(color, str):
            if color.startswith('#') and len(color) == 9: # #RRGGBBAA
                valid_color = color[:7] # Truncate to #RRGGBB
            elif color.startswith('#') and len(color) == 7: # #RRGGBB
                 valid_color = color
            # Could add check for valid CSS names here
            # elif color.lower() in ['red', 'blue', ...]: valid_color = color
        if valid_color:
            sanitized[key] = valid_color
        else:
             # Assign a default color if invalid or missing, cycling through Plotly's list
             sanitized[key] = default_colors[color_index % len(default_colors)]
             color_index += 1
             logging.warning(f"Invalid or missing color for '{key}'. Using default '{sanitized[key]}'. Original value: {color}")

    # Ensure essential colors have fallbacks if still missing after loop
    if 'curtailment' not in sanitized: sanitized['curtailment'] = 'red'
    if 'WAT' not in sanitized: sanitized['WAT'] = 'purple'

    return sanitized

# --- Tab Building Functions ---

def build_overview_tab(inputs, results, indicators, analysis_results):
    st.header("Simulation Overview")
    sanitized_colors = _sanitize_colors(commons['colors'])

    # --- Energy Balance Plot (Plotly) ---
    st.subheader("Energy Balance per Zone")
    all_zones = sorted(inputs['sets']['n'])
    if not all_zones:
        st.warning("No zones found in simulation data.")
        # Avoid showing multiselect if no options
    else:
        selected_zones_overview = st.multiselect(
            "Select Zones for Energy Balance & Capacity",
            options=all_zones,
            default=all_zones,
            key="overview_zone_select"
        )

        if selected_zones_overview:
            try:
                # Data Prep adapted from plot_energy_zone_fuel & previous implementation
                zones_in_results = indicators.Zone.unique()
                fuels = indicators[[u in commons['Technologies'] for u in indicators['Technology']]].Fuel.unique()
                GenPerZone = pd.DataFrame(0.0, index=zones_in_results, columns=list(fuels) + ['FlowIn', 'Curtailment'])
                # Ensure all standard fuel columns exist, fill with 0 if missing
                for f in commons['Fuels'] + ['FlowIn']:
                     if f not in GenPerZone.columns:
                         GenPerZone[f] = 0.0

                power_consumption = pd.DataFrame(0.0, index=results['OutputPowerConsumption'].index, columns=zones_in_results)
                storage_input = pd.DataFrame(0.0, index=results['OutputPower'].index, columns=zones_in_results)

                for z in zones_in_results:
                    if z in selected_zones_overview:
                        for f in fuels: GenPerZone.loc[z, f] = indicators[(indicators.Fuel == f) & (indicators.Zone == z)].Generation.sum()
                        net_imports = get_imports(results['OutputFlow'], z)
                        if net_imports > 0: GenPerZone.loc[z, 'FlowIn'] = net_imports
                        p_cons_zone = filter_by_zone(results['OutputPowerConsumption'], inputs, z); power_consumption[z] = p_cons_zone.sum(axis=1) if not p_cons_zone.empty else 0
                        storage_input_techs = filter_by_tech_list(results['OutputStorageInput'], inputs, commons['tech_storage']); s_in_zone = filter_by_zone(storage_input_techs, inputs, z); storage_input[z] = s_in_zone.sum(axis=1) if not s_in_zone.empty else 0
                        if z in results['OutputCurtailedPower']: GenPerZone.loc[z, 'Curtailment'] = -results['OutputCurtailedPower'][z].sum()

                GenPerZone = GenPerZone / 1E6 # TWh
                GenPerZone = GenPerZone.loc[selected_zones_overview]

                plot_cols_gen = [col for col in commons['MeritOrder'] if col in GenPerZone.columns and GenPerZone[col].abs().sum() > 1e-6] # Check for non-zero values more robustly
                GenPerZone_plot = GenPerZone[[c for c in plot_cols_gen + ['Curtailment'] if c in GenPerZone.columns]].copy()

                demand_da = inputs['param_df']['Demand']['DA'].sum(axis=0).reindex(selected_zones_overview).fillna(0) / 1E6
                p_cons_total = power_consumption.sum(axis=0).reindex(selected_zones_overview).fillna(0) / 1E6
                s_in_total = storage_input.sum(axis=0).reindex(selected_zones_overview).fillna(0) / 1E6
                total_demand_line = demand_da + p_cons_total
                storage_demand_line = total_demand_line + s_in_total

                GenPerZone_plot_reset = GenPerZone_plot.reset_index().rename(columns={'index': 'Zone'})
                energy_long = GenPerZone_plot_reset.melt(id_vars=['Zone'], var_name='Source', value_name='Energy_TWh')
                energy_long_pos = energy_long[energy_long['Energy_TWh'] >= 0].copy()
                energy_long_neg = energy_long[energy_long['Energy_TWh'] < 0].copy()

                fig_energy = px.bar(
                    energy_long_pos, x='Zone', y='Energy_TWh', color='Source',
                    title='Energy Balance per Zone (TWh)', labels={'Energy_TWh': 'Energy [TWh]'},
                    color_discrete_map=sanitized_colors,
                    category_orders={"Zone": sorted(selected_zones_overview)}
                )
                if not energy_long_neg.empty:
                     fig_energy.add_bar(x=energy_long_neg['Zone'], y=energy_long_neg['Energy_TWh'], name='Curtailment', marker_color=sanitized_colors.get('Curtailment', 'red'))

                zone_order_energy = sorted(selected_zones_overview)
                shapes_energy = []
                for i, zone in enumerate(zone_order_energy):
                    # Ensure the zone exists in the demand series before plotting shape
                    if zone in demand_da: shapes_energy.append(dict(type='line', x0=i-0.4, x1=i+0.4, y0=demand_da[zone], y1=demand_da[zone], xref='x', yref='y', line=dict(color='Black', width=2, dash='solid')))
                    if zone in total_demand_line: shapes_energy.append(dict(type='line', x0=i-0.4, x1=i+0.4, y0=total_demand_line[zone], y1=total_demand_line[zone], xref='x', yref='y', line=dict(color='Green', width=2, dash='dash')))
                    if zone in storage_demand_line: shapes_energy.append(dict(type='line', x0=i-0.4, x1=i+0.4, y0=storage_demand_line[zone], y1=storage_demand_line[zone], xref='x', yref='y', line=dict(color='Blue', width=2, dash='dot')))
                fig_energy.update_layout(shapes=shapes_energy, barmode='relative', legend_title_text='Source/Fuel')
                st.plotly_chart(fig_energy, use_container_width=True)
            except Exception as e: st.error(f"Error generating energy balance plot: {e}")
        else: st.info("Select zones for Energy Balance plot.")

    # --- Capacity Analysis Plot (Plotly) ---
    st.subheader("Installed Capacity per Zone")
    if selected_zones_overview:
        try:
            # Data Prep adapted from plot_zone_capacities
            units_df = inputs['units'].copy()
            valid_technologies = list(commons['Technologies'])
            power_units = units_df[units_df['Technology'].isin(valid_technologies)]
            if 'OutputPowerConsumption' in results: cons_units = units_df[units_df.index.isin(results['OutputPowerConsumption'].columns)]; power_units = pd.concat([power_units, cons_units[~cons_units.index.isin(power_units.index)]])

            power_units['TotalCapacity'] = power_units['PowerCapacity'] * power_units['Nunits']
            capacity_by_zone_fuel = power_units.groupby(['Zone', 'Fuel'])['TotalCapacity'].sum().unstack(fill_value=0) / 1000 # GW
            # Ensure all standard fuel columns exist for capacity, fill with 0 if missing
            for fuel in commons['MeritOrder']:
                if fuel not in capacity_by_zone_fuel.columns:
                    capacity_by_zone_fuel[fuel] = 0.0
            plot_cols = [col for col in commons['MeritOrder'] if col in capacity_by_zone_fuel.columns]
            capacity_by_zone_fuel = capacity_by_zone_fuel[plot_cols]
            capacity_by_zone_fuel = capacity_by_zone_fuel.loc[capacity_by_zone_fuel.sum(axis=1) > 0].reindex(selected_zones_overview).dropna(how='all')

            if not capacity_by_zone_fuel.empty:
                capacity_long = capacity_by_zone_fuel.stack().reset_index()
                capacity_long.columns = ['Zone', 'Fuel', 'Capacity_GW']
                capacity_long = capacity_long[capacity_long['Capacity_GW'] > 0]

                peak_demand_da = inputs['param_df']['Demand']['DA'].max(axis=0) / 1000
                peak_demand_total = peak_demand_da.copy()
                if 'OutputPowerConsumption' in results: peak_consumption = results['OutputPowerConsumption'].max(axis=0) / 1000; peak_demand_total = peak_demand_total.add(peak_consumption, fill_value=0)
                peak_demand_da = peak_demand_da.reindex(capacity_by_zone_fuel.index).fillna(0)
                peak_demand_total = peak_demand_total.reindex(capacity_by_zone_fuel.index).fillna(0)

                fig_cap = px.bar(
                    capacity_long, x='Zone', y='Capacity_GW', color='Fuel', title='Installed Capacity per Zone (GW)',
                    labels={'Capacity_GW': 'Installed Capacity [GW]'}, color_discrete_map=sanitized_colors
                )
                zone_order_cap = capacity_by_zone_fuel.index.tolist()
                fig_cap.update_xaxes(categoryorder='array', categoryarray=zone_order_cap)
                shapes = []
                for i, zone in enumerate(zone_order_cap):
                    if zone in peak_demand_da: shapes.append(dict(type='line', x0=i-0.4, x1=i+0.4, y0=peak_demand_da[zone], y1=peak_demand_da[zone], xref='x', yref='y', line=dict(color='Black', width=2)))
                    if zone in peak_demand_total: shapes.append(dict(type='line', x0=i-0.4, x1=i+0.4, y0=peak_demand_total[zone], y1=peak_demand_total[zone], xref='x', yref='y', line=dict(color='Blue', width=2)))
                fig_cap.update_layout(shapes=shapes, legend_title_text='Fuel Type')
                st.plotly_chart(fig_cap, use_container_width=True)
            else: st.info("No capacity data available for selected zones.")
        except Exception as e: st.error(f"Error generating capacity plot: {e}")
    else: st.info("Select zones for Capacity plot.")

    # --- Simulation Summary Metrics and Tables ---
    st.subheader("Simulation Summary")
    if analysis_results:
        # Aggregated Metrics
        st.markdown("**Aggregated Metrics** provide a high-level overview...")
        col1, col2, col3, col4 = st.columns(4)
        # Ensure values exist before formatting, provide default
        total_load_val = analysis_results.get('TotalLoad', 0) / 1e6
        curtailment_val = analysis_results.get('Curtailment', 0) / 1e6
        peak_load_val = analysis_results.get('PeakLoad', 0) / 1e3
        max_curtailment_val = analysis_results.get('MaxCurtailment', 0) / 1e3
        net_imports_val = analysis_results.get('NetImports', 0) / 1e6
        shed_load_val = analysis_results.get('ShedLoad', 0) # Assuming TWh
        cost_kwh_val = analysis_results.get('Cost_kwh', None)
        max_shed_load_val = analysis_results.get('MaxShedLoad', 0) / 1e3
        shifted_load_val = analysis_results.get('ShiftedLoad', 0) # Assuming TWh

        with col1: st.metric("Total Load (TWh)", f"{total_load_val:.2f}"); st.metric("Total Curtailment (TWh)", f"{curtailment_val:.2f}")
        with col2: st.metric("Peak Load (GW)", f"{peak_load_val:.2f}"); st.metric("Max Curtailment (GW)", f"{max_curtailment_val:.2f}")
        with col3: st.metric("Net Imports (TWh)", f"{net_imports_val:.2f}"); st.metric("Total Shed Load (TWh)", f"{shed_load_val:.2f}")
        with col4: st.metric("Avg Cost (€/MWh)", f"{cost_kwh_val:.2f}" if isinstance(cost_kwh_val, (int, float)) else "N/A"); st.metric("Max Shed Load (GW)", f"{max_shed_load_val:.2f}")
        if 'ShiftedLoad' in analysis_results: st.metric("Total Shifted Load (TWh)", f"{shifted_load_val:.2f}")

        # Detailed Summary Tables Expander
        with st.expander("View Detailed Summary Tables"):
            st.markdown("**Zone Specific Data** breaks down key metrics...")
            df_zone_data = analysis_results.get('ZoneData')
            if df_zone_data is not None and not df_zone_data.empty:
                st.dataframe(df_zone_data)
                st.download_button("Download Zone Data", convert_df_to_csv(df_zone_data), "zone_data.csv", "text/csv", key="dl_zone")
            else: st.info("No Zone Specific Data.")

            st.markdown("**Storage Data** provides insights...")
            df_storage_data = analysis_results.get('StorageData')
            if df_storage_data is not None and not df_storage_data.empty:
                df_storage_data = df_storage_data.dropna(axis=1, how='all')
                st.dataframe(df_storage_data)
                st.download_button("Download Storage Data", convert_df_to_csv(df_storage_data), "storage_data.csv", "text/csv", key="dl_storage")
            else: st.info("No Storage Data.")

            st.markdown("**Congestion Data** reveals bottlenecks...")
            congestion_data = analysis_results.get('Congestion')
            if congestion_data: st.write("Hours of congestion:"); st.json(congestion_data)
            else: st.info("No Congestion Data.")

            st.markdown("**Water Consumption Data** summarizes...")
            water_data = analysis_results.get('WaterConsumptionData')
            if water_data and water_data.get('ZoneLevel'):
                df_water_w = pd.DataFrame(water_data['ZoneLevel']['WaterWithdrawal'].sum(), columns=['TotalWithdrawal_m3'])
                df_water_c = pd.DataFrame(water_data['ZoneLevel']['WaterConsumption'].sum(), columns=['TotalConsumption_m3'])
                if not df_water_w.empty:
                    st.dataframe(df_water_w); st.download_button("Download Water Withdrawal", convert_df_to_csv(df_water_w), "water_withdrawal.csv", "text/csv", key="dl_water_w")
                else: st.info("No Water Withdrawal data.")
                if not df_water_c.empty:
                    st.dataframe(df_water_c); st.download_button("Download Water Consumption", convert_df_to_csv(df_water_c), "water_consumption.csv", "text/csv", key="dl_water_c")
                else: st.info("No Water Consumption data.")
            else: st.info("No Water Consumption Data.")
    else: st.warning("Could not retrieve simulation analysis results.")


def build_zone_analysis_tab(inputs, results):
    st.header("Detailed Zone Analysis")
    sanitized_colors = _sanitize_colors(commons['colors'])

    all_zones_list = sorted(inputs['sets']['n'])
    if not all_zones_list: st.warning("No zones found."); return

    selected_zone_detail = st.selectbox("Select Zone", options=all_zones_list, index=0, key="zone_detail_select")

    # Date Range Selection
    min_date, max_date = (results['OutputPower'].index.min().date(), results['OutputPower'].index.max().date()) if 'OutputPower' in results and not results['OutputPower'].empty else (datetime.date(2023, 1, 1), datetime.date(2023, 12, 31))
    col1_z, col2_z = st.columns(2)
    with col1_z: date_range_start_z = st.date_input("Start Date", value=min_date, min_value=min_date, max_value=max_date, key="zone_date_start")
    with col2_z:
        valid_end_date_z = max(date_range_start_z, min_date); default_end_date_z = min(date_range_start_z + datetime.timedelta(days=6), max_date)
        if default_end_date_z < valid_end_date_z: default_end_date_z = valid_end_date_z
        date_range_end_z = st.date_input("End Date", value=default_end_date_z, min_value=valid_end_date_z, max_value=max_date, key="zone_date_end")

    rng_zone = None
    try:
        start_dt, end_dt = datetime.datetime.combine(date_range_start_z, datetime.time.min), datetime.datetime.combine(date_range_end_z, datetime.time.max)
        full_range = pd.date_range(start=start_dt, end=end_dt, freq='h')
        if 'OutputPower' in results and not results['OutputPower'].empty: rng_zone = results['OutputPower'].index.intersection(full_range)
        else: rng_zone = full_range # Use if no results to intersect with
    except Exception as e: st.error(f"Error creating date range: {e}")

    hide_storage_z = st.checkbox("Hide Storage Level Plot", value=False, key="zone_hide_storage")
    show_rug_z = st.checkbox("Show Rug Plot (Generation)", value=True, key="zone_show_rug")

    if rng_zone is not None and not rng_zone.empty:
        st.subheader(f"Dispatch Plot for Zone: {selected_zone_detail}")
        # --- Revert to Matplotlib Dispatch Plot --- 
        try:
            plt.figure(figsize=(13, 7)) # Create figure context
            # Call original plot_zone function, which draws on the current figure
            plot_success = plot_zone(
                inputs, results,
                z=selected_zone_detail,
                rng=rng_zone,
                rug_plot=False, # Controlled separately below
                hide_storage_plot=hide_storage_z
                # colors=sanitized_colors # Pass sanitized colors if plot_zone accepts it
            )
            if plot_success:
                fig_dispatch_mpl = plt.gcf()
                st.pyplot(fig_dispatch_mpl)
                plt.close(fig_dispatch_mpl)
            else:
                st.warning("Could not generate dispatch plot.")
        except Exception as e: st.error(f"Error generating dispatch plot: {e}"); plt.close()

        # Rug Plot (Matplotlib)
        if show_rug_z:
            st.subheader(f"Generation Rug Plot for Zone: {selected_zone_detail}")
            try:
                zone_generation = filter_by_zone(results['OutputPower'], inputs, selected_zone_detail)
                if not rng_zone.intersection(zone_generation.index).empty:
                     zone_generation_ranged = zone_generation.loc[rng_zone]
                     if not zone_generation_ranged.empty and not zone_generation_ranged.isnull().all().all():
                          plt.figure(figsize=(16, max(1, 0.25 * len(zone_generation_ranged.columns))))
                          plot_rug(zone_generation_ranged, on_off=False, cmap='gist_heat_r', fig_title=selected_zone_detail)
                          st.pyplot(plt.gcf())
                          plt.close()
                     else: st.info(f"No generation data for rug plot in selected range.")
                else: st.info(f"Selected range has no overlap with rug plot data.")
            except Exception as e: st.error(f"Error generating rug plot: {e}"); plt.close()
    elif rng_zone is not None and rng_zone.empty: st.warning("Selected date range has no overlap with simulation data.")
    else: st.warning("Please select a valid date range.")


def build_network_tab(inputs, results):
    st.header("Network Analysis")
    # Simplified: Assuming Zone Analysis tab sets the date range
    # A more robust app might use session state or pass rng_zone explicitly
    st.info("Date range for Network Analysis is controlled in the 'Zone Analysis' tab.")

    # Placeholder for range - this needs a better mechanism in a real app
    # Attempting to retrieve from locals() is unreliable in Streamlit's execution model.
    # Let's default to the full range if not explicitly passed or set in session state.
    rng_net = results['OutputPower'].index if 'OutputPower' in results else pd.date_range(start='2023-01-01', periods=24*7, freq='h') # Fallback range
    # Consider adding date inputs here if independent range is desired

    # Map Plots
    st.subheader("Network Map Visualizations")
    with st.expander("View Network Maps", expanded=False):
        if not GEOPLOT_AVAILABLE: st.error("Cartopy dependency missing. Maps unavailable.")
        else:
            st.markdown("**Note:** Map generation can be slow...")
            col_map1, col_map2 = st.columns(2)
            with col_map1:
                if st.button("Generate Net Flows Map", key="gen_flow_map"):
                    with st.spinner("Generating Net Flows Map..."):
                        try: plt.figure(figsize=(10, 8)); plot_net_flows_map(inputs, results, idx=rng_net); st.pyplot(plt.gcf()); plt.close()
                        except Exception as e: st.error(f"Error: {e}"); plt.close()
            with col_map2:
                if st.button("Generate Congestion Map", key="gen_cong_map"):
                     with st.spinner("Generating Congestion Map..."):
                        try: plt.figure(figsize=(10, 8)); plot_line_congestion_map(inputs, results, idx=rng_net); st.pyplot(plt.gcf()); plt.close()
                        except Exception as e: st.error(f"Error: {e}"); plt.close()

    # Power Flow Tracing
    st.subheader("Power Flow Tracing")
    with st.expander("View Power Flow Tracing Matrix", expanded=False):
        try:
            with st.spinner("Calculating Power Flow Tracing..."):
                # Ensure figsize is passed correctly if needed by heatmap function internally
                data_pft, data_pft_prct = plot_power_flow_tracing_matrix(inputs, results, idx=rng_net)
            st.pyplot(plt.gcf()) # plot_power_flow_tracing_matrix calls plt.show(), st.pyplot gets the current fig.
            plt.close()

            show_pft_data = st.checkbox("Show PFT DataFrames", value=False, key="pft_data_check")
            if show_pft_data:
                if data_pft is not None: st.dataframe(data_pft / 1e6) # Assuming MW -> TWh
                if data_pft_prct is not None: st.dataframe(data_pft_prct)
            if data_pft is not None: st.download_button("Download PFT (TWh)", convert_df_to_csv(data_pft / 1e6), "pft_twh.csv", "text/csv", key="dl_pft_twh")
            if data_pft_prct is not None: st.download_button("Download PFT (%)", convert_df_to_csv(data_pft_prct), "pft_percent.csv", "text/csv", key="dl_pft_pct")
        except Exception as e: st.error(f"Error generating PFT plot: {e}"); plt.close()


def build_unit_analysis_tab(inputs, results, indicators, analysis_results):
    st.header("Unit Analysis")

    unit_data_available = analysis_results is not None and analysis_results.get('UnitData') is not None
    indicators_available = indicators is not None

    if not unit_data_available and not indicators_available:
        st.info("Load data to view Unit Analysis.")
        return

    # Display Filtered Tables
    st.subheader("Power Plant / Unit Indicators")
    st.markdown("Detailed operational metrics per unit (e.g., generation, CF, startups).")
    if indicators_available:
        st.dataframe(indicators) # Display unfiltered indicators
        if not indicators.empty: st.download_button("Download Indicators", convert_df_to_csv(indicators), "indicators.csv", "text/csv", key="dl_ind")
    else: st.warning("Indicators data not available.")

    st.subheader("Unit Specific Data (Costs, Emissions, Water)")
    st.markdown("Detailed cost, emission, and water usage per unit.")
    if unit_data_available:
        unit_data_df = analysis_results['UnitData'] # Get unfiltered data
        st.dataframe(unit_data_df)
        if not unit_data_df.empty: st.download_button("Download Unit Data", convert_df_to_csv(unit_data_df), "unit_data.csv", "text/csv", key="dl_unit")
    else: st.info("Unit Specific Data not in summary results.")


# --- New Tab: Data Explorer ---
def build_data_explorer_tab(inputs, results):
    st.header("Raw Data Explorer")
    st.markdown("Select a data source and then choose a specific table/item to view.")

    data_source = st.radio(
        "Select Data Source",
        ("Input Parameters (param_df)", "Output Results"),
        key="data_source_select"
    )

    source_dict = None
    if data_source == "Input Parameters (param_df)":
        if inputs and 'param_df' in inputs:
            source_dict = inputs['param_df']
        else:
            st.warning("Input parameters (`inputs['param_df']`) not found.")
            return
    elif data_source == "Output Results":
        if results:
            source_dict = results
        else:
            st.warning("Results data not found.")
            return

    if source_dict:
        # Get keys and sort them
        available_keys = sorted(list(source_dict.keys()))
        if not available_keys:
            st.info(f"No data items found in {data_source}.")
            return

        selected_key = st.selectbox(
            f"Select Item from {data_source}",
            options=available_keys,
            index=None, # Default to no selection
            placeholder="Choose an item...",
            key="data_item_select"
        )

        if selected_key:
            selected_item = source_dict.get(selected_key)

            st.subheader(f"Data for: `{selected_key}`")

            if isinstance(selected_item, pd.DataFrame):
                st.write(f"Shape: {selected_item.shape}")
                st.dataframe(selected_item)
                # Add download button for DataFrames
                csv_data = convert_df_to_csv(selected_item) # Use helper
                st.download_button(
                    label=f"Download {selected_key} as CSV",
                    data=csv_data,
                    file_name=f'{selected_key}.csv',
                    mime='text/csv',
                    key=f"dl_{selected_key}" # Dynamic key
                )
            elif isinstance(selected_item, pd.Series):
                st.write(f"Shape: {selected_item.shape}")
                st.dataframe(selected_item)
                # Add download button for Series
                csv_data = convert_df_to_csv(selected_item) # Use helper
                st.download_button(
                    label=f"Download {selected_key} as CSV",
                    data=csv_data,
                    file_name=f'{selected_key}.csv',
                    mime='text/csv',
                    key=f"dl_{selected_key}" # Dynamic key
                )
            elif isinstance(selected_item, dict):
                 st.write("Type: Dictionary")
                 st.json(selected_item, expanded=False) # Show JSON collapsed by default
            elif isinstance(selected_item, list):
                 st.write("Type: List")
                 st.write(selected_item)
            elif selected_item is not None:
                 st.write(f"Type: {type(selected_item).__name__}")
                 st.write(selected_item)
            else:
                 st.write("Item is None or empty.")
        else:
            st.info("Select an item from the dropdown above to view its contents.")


# --- Main UI Builder Function ---
@st.cache_data # Cache results analysis if inputs/results don't change
def run_and_get_analysis(_inputs, _results):
     return get_result_analysis(_inputs, _results)

def build_streamlit_app(inputs, results, indicators):
    """Builds the main Streamlit application UI using tabs."""

    if not (inputs and results and indicators is not None): # Check indicators explicitly
        st.warning("Data not fully loaded or indicators missing. Please load data using the sidebar.")
        return

    # Run analysis once
    analysis_results = run_and_get_analysis(inputs, results) # Use cached wrapper

    # Define Tabs
    tab_overview, tab_zone, tab_network, tab_unit, tab_explorer = st.tabs([
        "Overview", "Zone Analysis", "Network Analysis", "Unit Analysis",
        "Data Explorer"
    ])

    # Build Tabs
    with tab_overview:
        build_overview_tab(inputs, results, indicators, analysis_results)

    with tab_zone:
        # Note: This tab currently generates its own date range (rng_zone)
        # For network tab to use it, it ideally needs to be passed or stored in session state.
        build_zone_analysis_tab(inputs, results)

    with tab_network:
        # Pass necessary data. How rng_zone gets here needs careful consideration.
        # Simplification: Re-fetch or assume user navigated.
        build_network_tab(inputs, results)

    with tab_unit:
        build_unit_analysis_tab(inputs, results, indicators, analysis_results)

    with tab_explorer:
        build_data_explorer_tab(inputs, results) 