#!/usr/bin/env python3
"""Demo: Generate Popov plot data using the C library and visualize.
Requires: Python 3 with matplotlib and the compiled libcircle.a
Usage: python demos/demo_popov_plot.py
"""
import subprocess, os, sys

# Generate data using the example2 binary
print("Generating Popov criterion data...")
result = subprocess.run(["./build/example2_popov"], capture_output=True, text=True)
print(result.stdout)

print("\nTo visualize, run: python demos/plot_popov.py")
print("(This demo shows how to interface with the C library)")
