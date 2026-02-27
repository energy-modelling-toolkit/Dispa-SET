# Dispa-SET Postprocessing Streamlit Workflow Plan

This document outlines the plan to progressively replace the current script-based postprocessing of Dispa-SET results with an interactive Streamlit application.

## Phase 1: Initial Setup and Core Functionality

1.  **Basic Streamlit App Structure (`test_streamlit.py`)**
    *   `[ ]` Create a new Python script `scripts/test_streamlit.py`.
    *   `[ ]` Add basic Streamlit imports (`import streamlit as st`, `import dispaset as ds`).
    *   `[ ]` Add a title (`st.title`) for the application.
    *   `[ ]` Add a sidebar (`st.sidebar`) for controls.
    *   `[ ]` Implement a mechanism to select the simulation results path (e.g., `st.sidebar.text_input` or initially hardcoded for testing).

2.  **Data Loading and Caching**
    *   `[ ]` Implement a function within `test_streamlit.py` that uses `ds.get_sim_results(path)` to load inputs and results.
    *   `[ ]` Apply Streamlit's caching decorator (`@st.cache_data`) to the data loading function to avoid reloading on every interaction. This is crucial as `get_sim_results` can be slow.
    *   `[ ]` Add a button (`st.sidebar.button`) to trigger data loading based on the selected path.
    *   `[ ]` Display a confirmation message (`st.success` or `st.info`) once data is loaded.

3.  **Minimalist Plot: Energy per Zone by Fuel (`plot_energy_zone_fuel` equivalent)**
    *   `[ ]` Add a section in the main app area for this plot.
    *   `[ ]` Get the necessary data:
        *   Call `ds.get_indicators_powerplant(inputs, results)` (cache this result too if needed).
        *   Replicate the data processing logic within `plot_energy_zone_fuel` to prepare the `GenPerZone` DataFrame (or refactor `plot_energy_zone_fuel` to return the figure and axes).
    *   `[ ]` Use `st.multiselect` in the sidebar to allow users to select the zones (`ListZones`) to display. Default to all zones.
    *   `[ ]` Generate the stacked bar plot using matplotlib (as the original function does).
    *   `[ ]` Display the plot in Streamlit using `st.pyplot(fig)`. Ensure `show_plot=False` is passed to the original function if it's used directly, or capture the `fig` object.
    *   `[ ]` **Testing:** Run `streamlit run scripts/test_streamlit.py`. Verify that the app loads, data can be loaded, zones can be selected, and the energy balance plot is displayed correctly.

## Phase 2: Migrating Existing Plots

*Strategy: For each existing plot function, we will add controls to the Streamlit sidebar to configure the plot and display the output using `st.pyplot` initially. This minimizes changes to the original plotting logic.*

4.  **Zone Dispatch Plot (`plot_zone` equivalent)**
    *   `[ ]` Add a section/expander (`st.expander`) for detailed zone analysis.
    *   `[ ]` Add a `st.selectbox` in the sidebar to choose the zone (`z`).
    *   `[ ]` Add controls for the date range (`rng`):
        *   Use `st.slider` with `datetime` objects for selecting a range.
        *   Alternatively, use two `st.date_input` widgets for start and end dates, combined with `st.time_input` if hourly precision is needed.
    *   `[ ]` Add checkboxes (`st.checkbox`) to toggle:
        *   Storage level plot visibility (`hide_storage_plot`).
        *   Rug plot visibility (`rug_plot`).
    *   `[ ]` Replicate the data preparation logic within `plot_zone` (calling `get_plot_data`, preparing `level`, `demand`, `curtailment`, `shed_load`, etc.).
    *   `[ ]` Call the plotting logic (similar to `plot_dispatch` within `plot_zone`) to generate the matplotlib figure.
    *   `[ ]` Display the main dispatch plot using `st.pyplot(fig)`.
    *   `[ ]` If the rug plot is enabled, generate it and display it using a separate `st.pyplot(fig_rug)`.
    *   `[ ]` **Testing:** Run the Streamlit app. Select different zones and date ranges. Toggle the storage and rug plots. Verify correct display and absence of errors.

5.  **Boundary Sector Dispatch Plot (`plot_dispatchX` equivalent)**
    *   `[ ]` Add a separate section/expander (`st.expander`) for boundary sector analysis.
    *   `[ ]` Add a `st.selectbox` in the sidebar to choose the boundary sector zone (`z`).
    *   `[ ]` Reuse the date range selection controls from the `plot_zone` step.
    *   `[ ]` Replicate the data preparation logic from `plot_dispatchX`.
    *   `[ ]` Generate the plot using the existing matplotlib logic within `plot_dispatchX`.
    *   `[ ]` Display the plot using `st.pyplot(fig)`.
    *   `[ ]` **Testing:** Run the Streamlit app. Select different boundary sectors and date ranges. Verify correct display.

6.  **Zone Capacities Plot (`plot_zone_capacities` equivalent)**
    *   `[ ]` Add a section/expander for capacity analysis.
    *   `[ ]` Replicate the data processing logic from `plot_zone_capacities` to get the `PowerCapacity` DataFrame.
    *   `[ ]` Generate the stacked bar plot using matplotlib.
    *   `[ ]` Display the plot using `st.pyplot(fig)`.
    *   `[ ]` **Testing:** Run the Streamlit app. Verify the capacity plot renders correctly.

7.  **Map Plots (`plot_net_flows_map`, `plot_line_congestion_map`)**
    *   `[ ]` Add a section/expander for network/map visualizations.
    *   `[ ]` Determine the plotting library used by the original functions (likely matplotlib, potentially geopandas or folium).
    *   `[ ]` **Option 1 (Matplotlib):**
        *   Replicate data preparation.
        *   Generate the map plot using the existing matplotlib logic.
        *   Display using `st.pyplot(fig)`.
    *   `[ ]` **Option 2 (Folium):**
        *   Ensure `streamlit-folium` is installed.
        *   Replicate data preparation to create a Folium map object (`m`).
        *   Display using `st_folium(m)`. This offers interactivity.
    *   `[ ]` Implement equivalent for `plot_net_flows_map`.
    *   `[ ]` Implement equivalent for `plot_line_congestion_map`.
    *   `[ ]` **Testing:** Run the Streamlit app. Verify map plots render correctly.

8.  **Power Flow Tracing Matrix (`plot_power_flow_tracing_matrix` equivalent)**
    *   `[ ]` Add a section/expander for power flow tracing.
    *   `[ ]` Reuse date range selection controls.
    *   `[ ]` Replicate data preparation calling `get_power_flow_tracing`.
    *   `[ ]` Generate the heatmap using the matplotlib logic in `plot_power_flow_tracing_matrix`.
    *   `[ ]` Display using `st.pyplot(fig)`.
    *   `[ ]` Optionally, display the raw data matrices (`data`, `data_prct`) using `st.dataframe`.
    *   `[ ]` **Testing:** Run the Streamlit app. Verify heatmap renders correctly for a selected date range.

## Phase 3: Migrating Analysis Functions

*Strategy: Display the outputs of the analysis functions using Streamlit's text and data display elements.*

9.  **Result Analysis Summary (`get_result_analysis` equivalent)**
    *   `[ ]` Add a main section/expander for "Simulation Summary".
    *   `[ ]` Call `ds.get_result_analysis(inputs, results)`.
    *   `[ ]` Display scalar metrics (Total Load, Peak Load, Costs, Curtailment, Shed Load) using `st.metric`. Arrange them using `st.columns`.
    *   `[ ]` Display tabular data (`ZoneData`, `StorageData`, `UnitData`, `Congestion` dictionary) using `st.dataframe`. Use `st.write` or `st.subheader` for titles.
    *   `[ ]` Display the `FuelData` dictionary structure perhaps using nested expanders or formatted markdown.
    *   `[ ]` Use `st.expander` to group related tables/metrics if the output becomes too long.
    *   `[ ]` **Testing:** Run the Streamlit app. Verify all summary metrics and tables are displayed without errors.

10. **Power Plant Indicators (`get_indicators_powerplant` equivalent)**
    *   `[ ]` Add a section/expander for "Power Plant Indicators".
    *   `[ ]` Call `ds.get_indicators_powerplant(inputs, results)`.
    *   `[ ]` Display the resulting DataFrame using `st.dataframe`. Streamlit's dataframe display allows sorting and searching.
    *   `[ ]` **Testing:** Run the Streamlit app. Verify the power plant indicator table is displayed.

## Phase 4: Enhancements and New Features

*Strategy: Leverage Streamlit's interactive capabilities to go beyond static plots and tables.*

11. **Interactive Filtering and Exploration**
    *   `[ ]` For tables displayed with `st.dataframe` (e.g., UnitData, PPIndicators), users can already search and sort.
    *   `[ ]` Add `st.multiselect` widgets to filter tables based on criteria (e.g., filter UnitData by Zone, Technology, or Fuel). Update the displayed `st.dataframe` based on the selection.

12. **Interactive Charts**
    *   `[ ]` Replace some `st.pyplot` calls with native Streamlit charts or other interactive libraries where appropriate.
        *   *Dispatch Plot:* Could potentially be replotted using `st.area_chart` (for stacked view) or Plotly (`st.plotly_chart`) for zooming/panning/tooltips. This requires restructuring the data preparation significantly. Start with the main dispatch plot.
        *   *Bar Charts:* (`plot_energy_zone_fuel`, `plot_zone_capacities`) could use `st.bar_chart` or Plotly.
        *   *Storage Level:* Could use `st.line_chart` or Plotly.
    *   `[ ]` **Decision:** Prioritize converting the main `plot_dispatch` equivalent to an interactive Plotly chart first, as it benefits most from zooming/tooltips.

13. **Scenario Comparison**
    *   `[ ]` Modify the data loading mechanism to allow selecting and loading results from *multiple* simulation directories.
    *   `[ ]` Store the loaded results in a dictionary keyed by scenario name/path.
    *   `[ ]` Use `st.tabs` to display plots and analysis for different scenarios side-by-side or allow selecting scenarios to overlay on certain plots (e.g., capacity comparison).
    *   `[ ]` Adapt analysis functions (`get_result_analysis`) to show comparisons using `st.columns` and `st.metric`'s `delta` feature.

14. **Data Export**
    *   `[ ]` Provide options to download processed dataframes (e.g., `ZoneData`, `UnitData`, plot data).
    *   `[ ]` Use `st.download_button` to allow users to download data as CSV files.

15. **Customizable Plot Appearance**
    *   `[ ]` Add sidebar controls (`st.color_picker`, `st.selectbox` for line styles) to allow users to customize the appearance of key plots, especially if migrating to interactive libraries like Plotly.

16. **Map Interactivity (if using `st_folium`)**
    *   `[ ]` Explore using map click events to display detailed information about specific lines or zones in a sidebar or popup.

## Testing Strategy

*   **Incremental Testing:** After implementing each major numbered step above (especially plotting steps), run `streamlit run scripts/test_streamlit.py`.
*   **Manual Verification:** Check for:
    *   Absence of script errors.
    *   Correct rendering of plots and data.
    *   Responsiveness of interactive widgets (selectors, sliders, checkboxes).
    *   Correct filtering/updating of displays based on user input.
*   **Consider `AppTest`:** Once the core functionality is stable, consider using Streamlit's native testing framework (`AppTest`) for automated regression testing, particularly for data loading and basic widget interactions, although manual visual checks will remain important for plots.

