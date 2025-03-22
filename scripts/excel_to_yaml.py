#!/usr/bin/env python3
"""
Script to convert DispaSET Excel config files to YAML format.
"""

import os
import sys
import logging
from dispaset.preprocessing.data_handler import export_yaml_config

def main():
    # Set up logging
    logging.basicConfig(level=logging.INFO,
                       format='%(asctime)s - %(levelname)s - %(message)s')
    
    # Get the input Excel file path
    if len(sys.argv) > 1:
        excel_file = sys.argv[1]
    else:
        excel_file = 'ConfigFiles/ConfigTest.xlsx'
    
    # Create output YAML file path
    yaml_file = os.path.splitext(excel_file)[0] + '.yml'
    
    try:
        # Convert Excel to YAML
        success = export_yaml_config(excel_file, yaml_file)
        if success:
            logging.info(f"Successfully converted {excel_file} to {yaml_file}")
        else:
            logging.error(f"Failed to convert {excel_file}")
            sys.exit(1)
    except Exception as e:
        logging.error(f"Error converting config file: {str(e)}")
        sys.exit(1)

if __name__ == '__main__':
    main() 