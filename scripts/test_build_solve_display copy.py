# -*- coding: utf-8 -*-
"""
Minimalist example showing how to access the Dispa-SET API to:

* read a configuration file,
* create a simulation environment folder,
* run the simulation using GAMS,
* load the results and display a few plots.

By default this script now uses the short (3-day) integration config
``tests/configs/ultimate.yml`` so that it runs in only a few seconds.
This is the same configuration that powers the automated end-to-end
test in ``tests/ultimate/test_ultimate_full_pipeline.py``.

Run this script from the root folder of the Dispa-SET repository:

    python scripts/test_build_solve_display.py

@author: Sylvain Quoilin
"""
import os
import sys

import pandas as pd

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import dispaset as ds


def main(config_path: str = "tests/configs/ultimate.yml") -> None:
    config = ds.load_config(config_path)

    SimData = ds.build_simulation(config)
    ds.solve_GAMS(config["SimulationDirectory"], config["GAMS_folder"])

    inputs, results = ds.get_sim_results(
        path=config["SimulationDirectory"], cache=False
    )

    rng = pd.date_range(
        start=pd.Timestamp(*config["StartDate"][:3]),
        end=pd.Timestamp(*config["StopDate"][:3]),
        freq="h",
    )
    rng = pd.date_range(start='2016-01-01',end='2016-12-31',freq='h')


    ds.plot_zone(inputs, results, rng=rng)

    boundary_zones = inputs.get("sets", {}).get("nx", []) or []
    if "Z1_h2" in boundary_zones:
        ds.plot_dispatchX(inputs, results, z="Z1_h2")

    ds.plot_zone_capacities(inputs, results)
    ds.plot_energy_zone_fuel(
        inputs, results, ds.get_indicators_powerplant(inputs, results)
    )

    ds.get_result_analysis(inputs, results)


if __name__ == "__main__":
    main()
