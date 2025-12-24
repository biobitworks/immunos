---
source: /Users/byron/projects/immunos81/examples/quick_viz_test.py
relative: immunos81/examples/quick_viz_test.py
generated_at: 2025-12-23 10:28
---

```python
"""
Quick Visualization Test

Tests visualization functions with a small dataset, saving plots to files.
"""

import sys
import os
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from immunos81.data.cleveland_loader import ClevelandDataLoader
from immunos81.engine.trainer import Trainer
from immunos81.engine.recognizer import Recognizer, RecognitionStrategy
from immunos81.utils.visualizer import Immunos81Visualizer


def main():
    print("Testing Immunos-81 Visualizations...")

    # Load dataset
    loader = ClevelandDataLoader()
    antigens, library = loader.load_data()
    print(f"Loaded {len(antigens)} records")

    # Quick train/test split
    train_data = antigens[:200]
    test_data = antigens[200:250]  # Small test set for speed

    # Train
    print("Training...")
    trainer = Trainer(library)
    tcells = trainer.train(train_data)

    # Test
    print("Testing...")
    recognizer = Recognizer(tcells, strategy=RecognitionStrategy.RHA)
    results = recognizer.batch_recognize(test_data)

    # Create output directory
    output_dir = "/Users/byron/projects/immunos81/examples/plots"
    os.makedirs(output_dir, exist_ok=True)

    # Initialize visualizer
    viz = Immunos81Visualizer()

    # Test each visualization
    print("\nGenerating visualizations...")

    print("  1. Avidity distribution...")
    viz.plot_avidity_distribution(tcells, test_data[:30],
                                  save_path=f"{output_dir}/avidity_dist.png")

    print("  2. Confusion matrix...")
    viz.plot_confusion_matrix(test_data, results,
                             save_path=f"{output_dir}/confusion_matrix.png")

    print("  3. Decision confidence...")
    viz.plot_decision_confidence(results, test_data,
                                save_path=f"{output_dir}/confidence.png")

    print("  4. Variable importance...")
    viz.plot_variable_importance(tcells, test_data,
                                top_n=10,
                                save_path=f"{output_dir}/var_importance.png")

    print(f"\n✓ All visualizations saved to: {output_dir}/")
    print("\nGenerated files:")
    for file in os.listdir(output_dir):
        if file.endswith('.png'):
            print(f"  - {file}")


if __name__ == "__main__":
    main()

```
