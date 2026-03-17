#!/bin/bash

echo "Creating virtual environment..."
python3 -m venv scanpy_env

echo "Activating environment..."
source scanpy_env/bin/activate

echo "Upgrading pip..."
pip install --upgrade pip setuptools wheel

echo "Installing project dependencies..."
pip install -r requirements.txt

echo "Environment setup complete."
echo "Activate later using:"
echo "source scanpy_env/bin/activate"
