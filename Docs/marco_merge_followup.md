# Marco Merge Follow-Up

This note tracks what has been integrated from `marco` into `master` and
what remains to make the full data/config pipeline robust.

## Current Status

- Merged `marco` feature set into working tree (no commit yet).
- Resolved merge conflicts in:
  - `dispaset/GAMS/UCM.gms`
  - `dispaset/preprocessing/{build.py,data_check.py,data_handler.py,preprocessing.py}`
  - `dispaset/postprocessing/{data_handler.py,plot.py}`
  - `dispaset/__init__.py`
- Preserved recent `master` validation pattern (`DispaSETValidationError`)
  in conflict resolution where possible.
- Added fast recursive reserve integration test:
  - `tests/integration/test_reserves_generic.py`

## Short-Term Stabilization (Next)

- Validate exogenous reserve input compatibility for:
  - `aFRRUDemand`, `aFRRDDemand`, `FFRDemand`, `FCRDemand`
- Verify backward compatibility of old reserve names/fields in configs and docs.
- Check `OutputReserveProvision` schema stability in postprocessing and plotting.
- Confirm MTS + reserve paths work together in at least one tiny case.

## Data + Config Pipeline To Complete

- Align all config loaders (`yaml`, `excel`, config editor API) on the same
  reserve field names and defaults.
- Add one shared normalization layer in preprocessing for reserve input names
  (single source of truth, no duplicated mappings).
- Add explicit user-facing validation messages for missing reserve data files
  under exogenous reserve mode.
- Add documentation snippets showing minimal reserve-enabled config examples.

## Additional Tests To Add

- Failure test: exogenous reserve mode with malformed reserve CSV schema.
- Integration test: reserve demand with multi-zone load-share distribution.
- Integration test: reserve demand + MTS enabled.
- Unit test: reserve config-key normalization function (once added).

## Known Risks / Open Questions

- Full behavior parity with pre-merge reserve formulations still needs
  benchmarking on larger systems.
- Some `marco` postprocessing outputs may require naming cleanup for long-term
  API stability.
- Performance impact of all reserve constraints on larger cases not yet profiled.
