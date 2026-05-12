# Dispa-SET test-suite roadmap

This document complements the high-level summary in the project
`README.md`. It describes the layout and intent of the new test
suite, the gaps that remain, and a list of structural improvements
to the Dispa-SET codebase that would make future testing easier.

The audience is twofold:

* developers working on Dispa-SET itself;
* automation tools / LLMs iterating on the codebase to fix or add
  tests, who need a clear map of what already exists and what is
  worth adding next.

## 1. Suite layout

```
tests/
├── conftest.py           # global pytest fixtures (Agg backend, GAMS)
├── _helpers.py           # shared helpers (skip_if_no_gams, build_solve, ...)
├── run_all.py            # standalone runner -> tests/output/TEST_REPORT.md
├── configs/              # YAML configs used by the test cases
│   ├── tiny.yml          # 1 zone, 3 days, LP clustered (Units_tiny.csv)
│   ├── tiny_2zones.yml   # 2 zones with NTC
│   ├── tiny_milp.yml     # Integer clustering (unit commitment)
│   ├── tiny_bs.yml       # H2 boundary sector + electrolyser/fuel-cell
│   ├── tiny_mts.yml      # MTS preprocessing (original 2-zone config)
│   ├── tiny_mts_variants.yml  # MTS variant tests (HDAM/HPHS/BATS, cyclic/imposed)
│   ├── tiny_curtailment.yml  # over-sized VRE -> curtailment
│   ├── tiny_dcpf.yml     # 3-zone triangle, DC-Power-Flow transmission
│   ├── tiny_ntc_3zones.yml   # same topology with NTC (comparison partner)
│   └── ultimate.yml      # 3-day "everything" config (Units_testcase.csv)
├── data/                 # all CSVs the configs above read from
│   ├── Units_tiny.csv    # minimal plant set (no complex constraints)
│   ├── Units_dcpf.csv    # single gas unit in Z1 (DC-PF / NTC-3zone tests)
│   ├── Units_testcase.csv# fuller plant set (CHP, BS, hydro)
│   ├── Units_mts.csv     # MTS test system (1500 MW gas + HDAM/HPHS/BATS, STOHours=20)
│   ├── GridData_tiny.csv # 3-line triangle network for DC-PF tests
│   ├── boundary_sector/  # BS_*.csv used by tiny_bs.yml & ultimate.yml
│   ├── AvailabilityFactors/, Load_RealTime/, ...
│   └── reference/        # pre-baked simulation outputs for unit tests
├── output/               # all artifacts written by the test runs
├── unit/                 # fast, GAMS-free tests
├── integration/          # build / solve / read on tiny mock systems
├── failure/              # invalid config / data -> error checks
└── ultimate/             # full pipeline on ultimate.yml
```

Every test file:

* starts with a header docstring (`What this test does`, `How to run`);
* can be executed both via `pytest <path/to/file.py>` and as a script
  with `python <path/to/file.py>` (each file ends with
  `pytest.main([__file__, "-q"])`);
* uses the `Agg` matplotlib backend so plots never block;
* respects the time budgets discussed below.

## 2. Time budget

| Group       | Per-test target | Group target | Notes                              |
|-------------|-----------------|--------------|------------------------------------|
| Unit        | < 0.1 s         | < 5 s        | No GAMS solve, only Python logic.  |
| Integration | < 5 s           | < 60 s       | Tiny configs, 3-day horizons.      |
| Failure     | < 1 s           | < 5 s        | Mostly preprocessing-level checks. |
| Ultimate    | < 30 s          | < 30 s       | One run, 3 days, 3 zones.          |
| **Total**   | -               | **< 5 min**  | Enforced by pytest-timeout.        |

A current run on a recent laptop with `dispaset2` finishes in
**~20 seconds total** for the full suite, well within budget.

The intended test runner is the `dispaset2` conda environment:

```bash
conda run -n dispaset2 python -m pytest tests/ -q
```

or, using the standalone runner:

```bash
conda run -n dispaset2 python tests/run_all.py
```

## 3. Coverage matrix

Below is a map between Dispa-SET features and the tests that
currently exercise them. "F" = first-order coverage (the test would
fail if the feature was completely broken), "S" = secondary
coverage (the feature is exercised but not strictly asserted on).

| Feature                              | Unit | Integration | Failure | Ultimate |
|--------------------------------------|:---:|:----------:|:-------:|:--------:|
| GAMS Python API present              |     |    F       |         |    S     |
| `package_exists` / `get_gams_path`   |     |    F       |         |          |
| GDX round-trip                       |     |    F       |         |          |
| YAML config loader                   |  F  |            |   F     |    S     |
| Excel config loader                  |  F  |            |   F     |          |
| Path resolution (relative -> abs)    |  F  |            |         |    S     |
| `commons` constants consistency      |  F  |            |         |          |
| `str_handler` (GDX safe names)       |  F  |            |         |          |
| `preprocessing.utils` helpers        |  F  |            |         |          |
| `preprocessing.utils` property tests |  F  |            |         |          |
| `data_check.check_units`             |  F  |            |   F     |    S     |
| `data_check.check_chp`               |  F  |            |   F     |    S     |
| `data_check.check_clustering`        |  F  |            |         |    S     |
| `data_handler.NodeBasedTable`        |  F  |            |         |    S     |
| `data_handler.UnitBasedTable`        |  F  |            |         |    S     |
| `boundary_sector` preprocessing      |  F  |    F       |         |    S     |
| `build_simulation` artifacts         |     |    F       |         |    F     |
| LP clustered solve (single zone)     |     |    F       |         |    S     |
| LP clustered solve (multi-zone NTC)  |     |    F       |         |    S     |
| MILP / unit commitment               |     |    F       |         |    F     |
| Boundary sector (H2) solve           |     |    F       |         |    F     |
| CHP via boundary sector              |     |    F       |         |    F     |
| Curtailment of VRE                   |     |    F       |         |          |
| Mid-Term Scheduling (MTS) — regional | F   |    F       |         |          |
| MTS cyclic boundary condition (MTS=1)|     |    F       |         |          |
| MTS cyclic — dispatch annual (final≥initial) |  |  F  |     |          |
| MTS imposed boundary cond. (MTS=2)   | F   |    F       |         |          |
| MTS profile injection (build-only)   | F   |            |         |          |
| MTS all storage types (HDAM/HPHS/BATS)|    |    F       |         |          |
| MTS profiles respected in dispatch   |     |    F       |         |          |
| MTS with boundary sector (H2)        |     |    F       |         |          |
| `get_sim_results`                    |     |    F       |         |    F     |
| `get_indicators_powerplant`          |  F  |    F       |         |    F     |
| `get_result_analysis`                |  F  |    F       |         |    F     |
| `plot_zone` / `plot_zone_capacities` |  F  |    F       |         |    F     |
| `plot_dispatch` / `plot_dispatchX`   |  F  |    F       |         |    F     |
| Off-screen rendering of plots        |  F  |    F       |         |    F     |
| `PTDF_matrix` (DC-Power-Flow)        |  F  |    F       |         |          |
| NTC transmission (multi-zone)        |     |    F       |         |    S     |
| DC-Power-Flow transmission           |  F  |    F       |         |          |
| NTC vs DC-PF loop-flow comparison    |     |    F       |         |          |
| Reserve legacy aliases (Reserve2U/D, FFRLimit, InertiaLimit) | F |  |  |  |
| `data_check.check_FFRDemand`         |  F  |    S       |         |          |
| `data_check.check_FCRDemand`         |  F  |    S       |         |          |
| Reserve participation (nested dict)  |  F  |    F       |         |          |
| Reserve demand pipeline (all types)  |  F  |    F       |         |          |
| Solver KPI helpers (`get_solver_kpis`, `assert_feasible`) | F | F |  |      |

## 4. Step-by-step roadmap

The current state is **213 tests pass** (5 pre-existing MTS failures remain, all in
`test_solve_mts.py` and caused by a pandas `.str` accessor issue in
`build.py:484`) as of the latest run in the `dispaset2` conda environment
(`conda run -n dispaset2 python -m pytest tests/unit/ tests/integration/ tests/failure/ tests/ultimate/`).
The roadmap below lists the next blocks in suggested execution order.

### 4.1 Short-term (low effort, high value) — COMPLETED

1. **[DONE] Tighten the failure tests.** All `data_check.*` paths that
   previously aborted via `sys.exit(1)` now raise a dedicated
   `DispaSETValidationError` (defined in `dispaset/common.py`). All
   `pytest.raises(SystemExit)` / `pytest.raises((SystemExit, Exception))`
   calls in `tests/unit/test_data_check.py` and
   `tests/failure/test_invalid_units_data.py` have been narrowed to
   `pytest.raises(DispaSETValidationError)`. A latent `TypeError` in the
   CHPType validation (int index concatenated into a string) was also fixed.
2. **[DONE] Add property-based checks** with Hypothesis for
   `_mylogspace`, `_split_list`, `_merge_two_dicts`, and `select_units`.
   Twelve property tests live in `tests/unit/test_utils_property.py`.
   The tests document the actual semantics of `_mylogspace` (first element
   equals `-low`, not `low`) and the duplicate-element edge case in
   `_split_list` (detected via `assume`).
3. **[DONE] Add a `--profile` option to `tests/run_all.py`** that records
   the wall-clock time of each test via JUnit XML and writes a
   "slowest first" table to `TEST_REPORT.md`.

### 4.2 Medium-term

5. **[DONE] Pre-compute reference outputs.** `tests/data/reference/`
   contains pre-built `Inputs.gdx`, `Inputs.p` and `Results.gdx`.
   `tests/unit/test_postprocess_reference.py` exercises the full
   postprocessing pipeline (gdx → dataframe → indicators → plots) without
   running a new GAMS solve.
6. **[DONE] Cover all reserve-related options** (`Reserve2U`, `Reserve2D`,
   `FFRLimit`, `PrimaryReserveLimit`, `InertiaLimit`).
   - `tests/unit/test_reserve_config.py` — `TestLegacyDemandAliases` now
     includes `test_ffrlimit_maps_to_ffrdemand` and
     `test_inertialimit_maps_to_systeminertia`; all six legacy-alias paths
     are tested.
   - `tests/unit/test_reserve_demand_pipeline.py` — two new test classes
     `TestCheckFFRDemand` (4 tests) and `TestCheckFCRDemand` (4 tests)
     cover the zone-aggregated validation logic (valid, negative, exceeds
     load, multi-zone).
   - Integration coverage via `tests/integration/test_reserves_generic.py`
     (cumulative FFR + FCR + mFRRU demand scenarios).
   - `FrequencyStability` is a GAMS parameter; it has no dedicated
     `data_check` path and does not need a unit test.
7. **[DONE] Cover all transmission-grid types** (`NTC`, `DC-Power-Flow`).
   Five new tests cover both modes on an identical 3-zone triangle topology:
   - `tests/unit/test_ptdf.py` — 5 GAMS-free unit tests: PTDF matrix shape,
     hand-calculated values (tol=1e-4), loop-flow prediction, and build-time
     PTDF/NaN assertions for both grid types.
   - `tests/integration/test_solve_dcpf.py` — 2 GAMS integration tests:
     `test_solve_dcpf_3zones` (DC-PF builds and solves, PTDF non-zero) and
     `test_ntc_vs_dcpf_loop_flows` (loop flow on ZZ3→Z2 is ≥ 5 MW larger
     under DC-PF than NTC, confirming Kirchhoff constraints are active).
   New data: `tests/data/GridData_tiny.csv` (asymmetric-reactance triangle),
   `tests/data/Units_dcpf.csv` (single gas unit in Z1),
   `tests/configs/tiny_dcpf.yml`, `tests/configs/tiny_ntc_3zones.yml`.
8. **[DONE] Add MTS variant tests** (cyclic, imposed, all storage types,
   profile injection, dispatch adherence, boundary-sector MTS).
   Six new tests live in `tests/integration/test_solve_mts.py`:
   - `test_mts_profile_injection_build_only` — GAMS-free; injects a flat 0.8
     profile for HDAM and verifies `StorageProfile['val']` mean > 0.7.
   - `test_mts_cyclic_boundary_condition` — HydroScheduling=1; runs
     `mid_term_scheduling` and checks the annual HDAM profile is valid (values
     in [0,1], non-empty).
   - `test_mts_imposed_boundary_condition` — HydroScheduling=2; verifies at
     build level that `StorageProfile` is a linspace from ReservoirLevelInitial
     (0.5) to ReservoirLevelFinal (0.3), monotonically non-increasing; also
     checks that `mid_term_scheduling` runs and returns a valid profile.
   - `test_mts_all_storage_types_full_pipeline` — full MTS + dispatch with
     HDAM, HPHS, and long-duration BATS (STOHours=20 each); asserts all three
     units appear in `asu`, `OutputStorageLevel` is present, BATS has a
     non-zero `StorageProfile`.
   - `test_mts_cyclic_dispatch_annual` — full-year rolling dispatch (monthly
     windows, HorizonLength=30) after MTS=1; asserts
     `OutputStorageLevel[-1] >= StorageInitial_normalised - 0.10`.  Completes
     in ≈ 22 s (1 MTS solve + 12 monthly dispatch solves).
   - `test_mts_profiles_respected_in_dispatch` — asserts that for each storage
     type, `OutputStorageLevel[-1] >= StorageProfile[-1] - tol`.
   - `test_mts_bs_regional` — runs `tiny_bs.yml` with HydroScheduling=1 and
     checks that `OutputSectorXStorageLevel` is produced.
   New test data: `tests/data/Units_mts.csv` (4-unit single-zone system with
   1500 MW gas backup, HDAM 100 MW / 50 000 MWh (STOHours=500, seasonal),
   HPHS and BATS with STOHours=20).
   New config: `tests/configs/tiny_mts_variants.yml` (annual MTS, 7-day
   dispatch, single zone Z1, LP clustered).
9. **[DONE] Add `solve_succeeded`-style assertions on solver KPIs**
   (objective value bounded, no spurious lost load on a feasible
   case).
   - `tests/_helpers.py` gains two new helpers:
     - `get_solver_kpis(results)` — returns `total_cost`,
       `total_lostload_mwh` (sum of all `LostLoad_*` variables except
       storage violation) and `storage_violation_mwh`.
     - `assert_feasible(results, max_lostload_mwh=1.0)` — asserts no
       significant load shedding and a finite non-negative objective value.
   - `tests/unit/test_kpi_helper.py` — 15 GAMS-free tests covering both
     helpers: cost aggregation, per-key lost-load summation, storage
     violation separation, negative-cost detection, non-finite cost
     detection.
   - `tests/integration/test_solve_lp.py` and `test_solve_milp.py` call
     `assert_feasible(results)` and check `kpis["total_cost"] >= 0`.
   - `build_solve()` now populates `out["kpis"]` automatically.

### 4.3 Long-term

10. **Snapshot tests for build outputs.** Compare `Inputs.gdx`
    against a reference produced by a known-good commit. This makes
    refactors much safer.
11. **Schema-validate every YAML config** with `pydantic` or
    `voluptuous` so the user gets a clean error before the
    optimization is built.
12. **Document each feature with a "minimum reproducer" config** in
    `Docs/feature_examples/`. Each reproducer becomes both a doc
    snippet and an integration test.

## 5. Suggestions to clean up the Dispa-SET codebase

The test-suite work surfaced a number of small frictions in the
codebase. None are urgent but each would make future testing
significantly easier:

1. **[DONE] Replace `sys.exit(1)` with raising a custom exception.**
   `DispaSETValidationError` is now raised throughout the entire Dispa-SET
   preprocessing and postprocessing pipeline:
   - `dispaset/preprocessing/data_check.py` — original replacement
   - `dispaset/preprocessing/build.py` — date range, time step, SimulationType,
     reservoir levels, CHPType, missing sim path, GAMS file errors
   - `dispaset/preprocessing/preprocessing.py` — MTS result validation
   - `dispaset/postprocessing/data_handler.py` — GAMS status type, index length,
     parameter dimensionality
   - `dispaset/postprocessing/postprocessing.py` — invalid Inputs type
   - `dispaset/postprocessing/plot.py` — mismatched time-series indices
   - `dispaset/misc/str_handler.py` — unsupported argument types
   - `dispaset/misc/gdx_handler.py` — GAMS not found, GDX file missing,
     variable shape mismatches
   Only commented-out `sys.exit` lines and CLI/entry-point uses remain.
   A `SimulationType` validation guard was also added to `build.py` so that
   unknown values (e.g. `"Banana clustered"`) immediately raise
   `DispaSETValidationError` rather than silently defaulting to MILP.
2. **[DONE] Fix the latent string concatenation bug at
   `data_check.py` (CHPType validation).** The unit index `u` was an
   `int` concatenated into a `str`, producing a `TypeError` instead of
   the intended critical log line. Fixed to use `str(unitname)` and
   `str(plants.loc[u, 'CHPType'])` separately.
3. **Make path resolution explicit.** Currently
   `data_handler.load_config_yaml` calculates `basefolder` as the
   directory of the YAML file and resolves relative paths against
   it. The convention is fine but the implementation is undocumented
   and made the test configs trickier than needed (paths in the
   tests must be relative to `tests/`, but fuel prices live in
   `../Database/`). A `paths_basefolder` field in the YAML, or a
   helper API `resolve_path(cfg, key)`, would remove the foot-gun.
4. **Split `dispaset/preprocessing/build.py`.** It mixes
   set-building, parameter-building, derived-quantity computation
   and per-feature glue. Splitting into smaller modules
   (`build_sets`, `build_parameters`, `build_boundary_sector`, ...)
   would let unit tests target specific concerns without running
   the whole pipeline.
5. **Remove the legacy `heat` / `h2` formulations.** They have been
   replaced by the boundary-sector abstraction but still leave
   traces in column names (`CHPMaxHeat`, `STOHours`, the `nx_CC`
   set) and in default values
   (`default.CostHeatSlack`). Deleting the dead code would shrink
   the surface area to test.
6. **Provide a tiny built-in dataset.** `tests/data/Units_tiny.csv`
   is a good starting point - exposing it via
   `dispaset.examples.load_tiny()` would let users prototype without
   ever copying CSVs around, and would give the test suite a
   well-defined "official" toy system.
7. **Add a `dispaset.api.solve(config)` one-shot entry point.**
   The integration tests have to call `build_simulation` and
   `solve_GAMS` separately, then `get_sim_results`. A single
   `solve(config)` returning `(inputs, results)` would simplify
   both the tests and the documentation examples.
8. **Modernize matplotlib calls.** Many `plot_*` helpers call
   `plt.show()` unconditionally. They should accept a
   `show: bool = True` flag (or return the `Figure`) so that
   automation can disable display without juggling backends.
9. **Adopt `pyproject.toml`-only configuration.** A single
   `pyproject.toml` (with `[tool.pytest.ini_options]`,
   `[tool.ruff]`, ...) would replace `setup.py`, `setup.cfg`, and
   any ad-hoc `tests/conftest.py` settings.
10. **Pin numerical dependencies.** A handful of tests trigger
    `DeprecationWarning`s from the pandas / numpy versions installed
    by `environment.yml`. Bumping to recent compatible versions and
    pinning them would silence the noise and prevent surprise
    failures on fresh installs.

## 6. How to extend the suite

Adding a new test should usually follow this three-step pattern:

1. **Pick the right tier**:
   * unit (`tests/unit/`) if the test does not require GAMS;
   * integration (`tests/integration/`) for build / solve / read;
   * failure (`tests/failure/`) for "wrong inputs -> clear error";
   * ultimate (`tests/ultimate/`) for end-to-end checks.
2. **Reuse existing fixtures.** `tests/_helpers.py` exposes
   `skip_if_no_gams`, `load_test_config`, `fresh_simdir`,
   `build_solve`, `solve_succeeded`. Avoid duplicating GAMS path
   detection in new files.
3. **Keep it short.** A new integration test should solve at most
   3 days, on at most 2 zones, with at most a handful of units.
   If you find yourself loading a yearly profile, you are probably
   in the wrong tier.

When the test exposes a real bug, leave a TODO comment with a clear
explanation rather than skipping the test. The test report
generated by `tests/run_all.py` will surface the failure on each
run.
