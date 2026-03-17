#!/bin/bash

# Create virtual environment
python3 -m venv scanpy_env

# Activate environment
source scanpy_env/bin/activate

# Upgrade installer tools
pip install --upgrade pip setuptools wheel

# Install dependencies
pip install scanpy scikit-misc

echo "Environment ready. Activate with:"
echo "source scanpy_env/bin/activate"
