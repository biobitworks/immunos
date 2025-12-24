---
source: /Users/byron/projects/rockbeatspaper/notebooks/example_script.py
relative: rockbeatspaper/notebooks/example_script.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
Example Python Script for THRML-HACK Project

This script demonstrates basic data analysis using numpy and pandas.

To run this script:
1. First, set up the Python environment:
   ./scripts/setup-python-env

2. Activate the virtual environment:
   source .venv/bin/activate

3. Run the script:
   python notebooks/example_script.py
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime


def main():
    print("=" * 60)
    print("THRML-HACK Example Python Script")
    print("=" * 60)
    print(f"Execution time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    # Generate sample data
    print("📊 Generating sample data...")
    np.random.seed(42)
    
    # Create sample thermal data
    time_points = np.linspace(0, 10, 100)
    temperature = 20 + 5 * np.sin(time_points) + np.random.normal(0, 0.5, 100)
    
    # Create a DataFrame
    df = pd.DataFrame({
        'time': time_points,
        'temperature': temperature,
        'status': ['normal' if t < 25 else 'elevated' for t in temperature]
    })
    
    print(f"✓ Generated {len(df)} data points")
    print()

    # Display basic statistics
    print("📈 Temperature Statistics:")
    print("-" * 40)
    print(f"Mean temperature: {df['temperature'].mean():.2f}°C")
    print(f"Min temperature:  {df['temperature'].min():.2f}°C")
    print(f"Max temperature:  {df['temperature'].max():.2f}°C")
    print(f"Std deviation:    {df['temperature'].std():.2f}°C")
    print()

    # Count status
    print("🔍 Status Distribution:")
    print("-" * 40)
    status_counts = df['status'].value_counts()
    for status, count in status_counts.items():
        print(f"{status}: {count} ({count/len(df)*100:.1f}%)")
    print()

    # Display first few rows
    print("📋 Sample Data (first 5 rows):")
    print("-" * 40)
    print(df.head().to_string(index=False))
    print()

    print("=" * 60)
    print("✓ Script completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()

```
