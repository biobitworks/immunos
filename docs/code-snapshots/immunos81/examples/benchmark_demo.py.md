---
source: /Users/byron/projects/immunos81/examples/benchmark_demo.py
relative: immunos81/examples/benchmark_demo.py
generated_at: 2025-12-23 10:28
---

```python
"""
Performance Benchmark Demo

Demonstrates parallel vs sequential performance comparison,
reproducing the ~40% speedup from the Immunos 99 paper.
"""

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from immunos81.data.cleveland_loader import ClevelandDataLoader
from immunos81.utils.benchmark import Benchmark


def main():
    print("="*60)
    print("IMMUNOS-81: Performance Benchmark")
    print("Comparing Sequential vs Parallel Implementations")
    print("="*60)
    print()

    # Load dataset
    print("Loading Cleveland dataset...")
    loader = ClevelandDataLoader()
    antigens, library = loader.load_data()
    print(f"Loaded {len(antigens)} records\n")

    # Split data
    n_train = int(len(antigens) * 0.7)
    train_data = antigens[:n_train]
    test_data = antigens[n_train:]

    print(f"Training set: {len(train_data)} instances")
    print(f"Test set: {len(test_data)} instances")
    print()

    # Benchmark options
    print("Select benchmark type:")
    print("  1. Training performance (3 runs)")
    print("  2. Recognition performance (3 runs)")
    print("  3. End-to-end pipeline (single run)")
    print("  4. All benchmarks")
    print()

    choice = input("Enter choice (1-4) [default: 3]: ").strip() or "3"

    n_workers = 4
    print(f"\nUsing {n_workers} worker threads for parallel version\n")
    print("="*60)
    print()

    if choice in ["1", "4"]:
        print("\n### BENCHMARK 1: Training Performance ###\n")
        Benchmark.benchmark_training(train_data, n_workers=n_workers, n_runs=3)

    if choice in ["2", "4"]:
        print("\n### BENCHMARK 2: Recognition Performance ###\n")
        Benchmark.benchmark_recognition(train_data, test_data, n_workers=n_workers, n_runs=3)

    if choice in ["3", "4"]:
        print("\n### BENCHMARK 3: End-to-End Pipeline ###\n")
        Benchmark.benchmark_end_to_end(train_data, test_data, n_workers=n_workers)

    print("\n" + "="*60)
    print("Benchmark Complete!")
    print("="*60)
    print()
    print("📊 Expected speedup from Immunos 99 paper: ~40% (1.67x)")
    print("💡 Actual speedup depends on CPU cores and dataset size")


if __name__ == "__main__":
    main()

```
