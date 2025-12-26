#!/bin/bash
# Run all tests and Python scripts

set -e

echo "=== Step 1: Running LOSAT tests ==="
bash run_comparison.sh

echo ""
echo "=== Step 2: Running Python scripts ==="
echo "2.1: plot_comparison.py"
python3 plot_comparison.py

echo "2.2: plot_execution_time.py"
python3 plot_execution_time.py

echo "2.3: plot_overall_trend.py"
python3 plot_overall_trend.py

echo ""
echo "=== All tests and scripts completed! ==="

