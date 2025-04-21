# Refactoring Plan for build.py

This document outlines the steps to refactor the `build.py` file in the Dispa-SET model, focusing on breaking down long functions and improving code organization. The primary focus is on interconnection-related code, which will be moved to a dedicated module.

## Motivation

The `build.py` file has grown too large and contains functions that are too long, making it difficult to understand and maintain. By breaking it down into more focused modules, we can:

1. Improve code organization and readability
2. Make the codebase easier to maintain
3. Enable better reuse of code across different parts of the model
4. Facilitate testing of individual components

## Refactoring Steps

### Step 1: Create a dedicated interconnections module ✓

- [-] Create a new file `dispaset/preprocessing/interconnections.py`
- [-] Move interconnection-related functions from `utils.py` to the new module:
  - [-] `interconnections()` function
  - [-] `merge_lines()` function
  - [-] `incidence_matrix()` function (since it's primarily used for interconnection matrices)
  - [-] `adjust_ntc()` function
- [-] Add a new function `create_price_transmission_data()` to encapsulate the price transmission data creation logic
- [-] Update imports in `build.py` to use the new module
- [-] Update the code in `build.py` to use the new `create_price_transmission_data()` function
- [-] Run tests to ensure everything still works

**Status**: Completed and verified ✓  
**Results**: All tests pass. The new module successfully encapsulates interconnection-related functionality.

### Step 2: Create a dedicated boundary sector module

- [-] Create a new file `dispaset/preprocessing/boundary_sector.py`
- [-] Move boundary sector-related functions from `build.py` and `utils.py` to the new module:
  - [-] Create a function to handle boundary sector interconnections
  - [-] Create a function to handle boundary sector spillage
  - [-] Move the zone_to_bs_mapping function and other BS-related functions
- [-] Update imports in `build.py` to use the new module
- [-] Update the code in `build.py` to use the new functions
- [-] Run tests to ensure everything still works

**Status**: Completed and verified ✓  
**Results**: All tests pass. The new module successfully encapsulates boundary sector-related functionality.

### Step 3: Refactor parameter initialization

- [-] Create a new file `dispaset/preprocessing/parameters.py`
- [-] Move parameter-related initialization code from `build.py` to the new module:
  - [-] Create a function to initialize default parameters
  - [-] Create a function to update parameters from plant data
  - [-] Create a function to set up storage parameters
- [-] Update imports in `build.py` to use the new module
- [-] Update the code in `build.py` to use the new functions
- [-] Run tests to ensure everything still works

**Status**: Completed and verified ✓  
**Results**: All tests pass. The new module successfully encapsulates parameter-related functionality.

### Step 4: Refactor sets creation

- [-] Create a new file `dispaset/preprocessing/sets.py`
- [-] Move set-related initialization code from `build.py` to the new module:
  - [-] Create a function to initialize sets
  - [-] Create a function to update sets based on zones and units
- [-] Update imports in `build.py` to use the new module
- [-] Update the code in `build.py` to use the new functions
- [-] Run tests to ensure everything still works

**Status**: Completed and verified ✓  
**Results**: All tests pass. The new module successfully encapsulates sets creation functionality.

### Step 5: Refactor simulation environment setup

- [-] Create a new file `dispaset/preprocessing/simulation.py`
- [-] Move simulation environment setup code from `build.py` to the new module:
  - [-] Create a function to set up the simulation environment
  - [-] Create a function to write GDX files
  - [-] Create a function to write GAMS files
- [-] Update imports in `build.py` to use the new module
- [-] Update the code in `build.py` to use the new functions
- [-] Run tests to ensure everything still works

**Status**: Completed and verified ✓  
**Results**: All tests pass. The new module successfully encapsulates simulation environment setup functionality.

### Step 6: Update tests and documentation

- [-] Ensure all tests pass with the refactored code
- [-] Update documentation to reflect the new structure
- [-] Add docstrings to new functions
- [-] Create examples of how to use the new modules

**Status**: Completed and verified ✓  
**Results**: All tests pass and documentation has been updated to reflect the new module structure.

## Progress Tracking

- [-] Step 1: Create a dedicated interconnections module ✓
- [-] Step 2: Create a dedicated boundary sector module ✓
- [-] Step 3: Refactor parameter initialization ✓
- [-] Step 4: Refactor sets creation ✓
- [-] Step 5: Refactor simulation environment setup ✓
- [-] Step 6: Update tests and documentation ✓

## Conclusion

The refactoring of the `build.py` file has been successfully completed. The code is now more modular, easier to maintain, and better organized. Each module has a clear responsibility and can be tested independently. Future development can now focus on specific modules without impacting the entire codebase. 
