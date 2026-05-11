#!/usr/bin/env bash
# Build the Dispa-SET HTML documentation locally using the dispaset2 conda environment.
# Usage: bash scripts/build_docs.sh [--clean]
#
# Requirements:
#   - conda must be installed and on PATH (or CONDA_PREFIX/etc/profile.d/conda.sh must exist)
#   - The 'dispaset2' conda environment must have been created from environment.yml
#
# The HTML output will be in Docs/_build/html/index.html

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
DOCS_DIR="$REPO_ROOT/Docs"
CLEAN=${1:-}

# ---------------------------------------------------------------------------
# Activate conda environment
# ---------------------------------------------------------------------------
if [ -n "${CONDA_EXE:-}" ]; then
    CONDA_BASE="$(dirname "$(dirname "$CONDA_EXE")")"
elif command -v conda &>/dev/null; then
    CONDA_BASE="$(conda info --base 2>/dev/null)"
else
    echo "ERROR: conda not found. Please install Miniconda/Anaconda and re-run."
    exit 1
fi

CONDA_SH="$CONDA_BASE/etc/profile.d/conda.sh"
if [ ! -f "$CONDA_SH" ]; then
    echo "ERROR: Cannot find $CONDA_SH. Is conda correctly installed?"
    exit 1
fi

# shellcheck source=/dev/null
source "$CONDA_SH"
conda activate dispaset2

# ---------------------------------------------------------------------------
# Install/upgrade documentation dependencies
# ---------------------------------------------------------------------------
echo "Installing documentation dependencies..."
pip install --quiet --upgrade \
    "sphinx>=5" \
    sphinx-rtd-theme \
    docutils \
    numpy \
    pandas \
    scipy \
    matplotlib \
    xarray \
    pyyaml \
    networkx

# ---------------------------------------------------------------------------
# Build the docs
# ---------------------------------------------------------------------------
cd "$DOCS_DIR"

if [ "$CLEAN" = "--clean" ]; then
    echo "Cleaning previous build..."
    make clean
fi

echo "Building HTML documentation..."
make html

echo ""
echo "Documentation built successfully."
echo "Open: file://$DOCS_DIR/_build/html/index.html"
