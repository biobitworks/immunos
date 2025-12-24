---
source: /Users/byron/projects/immunos81/examples/cleveland_demo.py
relative: immunos81/examples/cleveland_demo.py
generated_at: 2025-12-23 10:28
---

```python
"""
Immunos-81 Cleveland Heart Disease Demo

This demo reproduces the results from the Immunos-81 paper using the
Cleveland coronary artery disease dataset.

Expected results:
- SHA: 83.2% accuracy
- RHA: 85.5% accuracy (excluding uncertain cases)
- RHA: 11.2% uncertain cases
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from immunos81.data.cleveland_loader import ClevelandDataLoader
from immunos81.utils.cross_validator import CrossValidator
from immunos81.engine.trainer import Trainer
from immunos81.engine.recognizer import Recognizer, RecognitionStrategy


def main():
    """Run the Immunos-81 Cleveland demo"""

    print("="*60)
    print("IMMUNOS-81: Cleveland Heart Disease Classification")
    print("="*60)
    print()

    # Load dataset
    print("Loading Cleveland Heart Disease dataset...")
    loader = ClevelandDataLoader()
    antigens, library = loader.load_data()

    print(f"Loaded {len(antigens)} patient records")
    print(f"Variables: {len(library)} features")
    print()

    # Show class distribution
    class_counts = {}
    for antigen in antigens:
        label = antigen.class_label
        class_counts[label] = class_counts.get(label, 0) + 1

    print("Class distribution:")
    for label, count in sorted(class_counts.items()):
        print(f"  {label}: {count} ({count/len(antigens)*100:.1f}%)")
    print()

    # Option 1: Quick demo with train/test split
    print("Option 1: Quick Demo (70/30 split)")
    print("-" * 60)
    quick_demo(antigens)

    # Option 2: Full 10-fold cross-validation (reproduces paper)
    print("\n\nOption 2: 10-Fold Cross-Validation (reproduces paper results)")
    print("-" * 60)
    cross_validation_demo(antigens)


def quick_demo(antigens):
    """Quick demonstration with a single train/test split"""

    # Split data 70/30
    n_train = int(len(antigens) * 0.7)
    train_data = antigens[:n_train]
    test_data = antigens[n_train:]

    print(f"Training on {len(train_data)} samples, testing on {len(test_data)} samples\n")

    # Train
    print("Training Immunos-81...")
    trainer = Trainer(train_data[0].library)
    tcells = trainer.train(train_data)
    print()

    # Test with SHA
    print("Testing with SHA (Simple Highest Avidity)...")
    recognizer_sha = Recognizer(tcells, strategy=RecognitionStrategy.SHA)
    metrics_sha = recognizer_sha.evaluate(test_data)

    print(f"  Accuracy: {metrics_sha['accuracy_overall']:.1%}")
    print(f"  Correct: {metrics_sha['correct']}/{metrics_sha['total']}")
    print()

    # Test with RHA
    print("Testing with RHA (Relative Highest Avidity, 5% threshold)...")
    recognizer_rha = Recognizer(tcells, strategy=RecognitionStrategy.RHA)
    metrics_rha = recognizer_rha.evaluate(test_data)

    print(f"  Accuracy (all): {metrics_rha['accuracy_overall']:.1%}")
    print(f"  Accuracy (certain only): {metrics_rha['accuracy_certain']:.1%}")
    print(f"  Uncertain: {metrics_rha['uncertain_attempts']}/{metrics_rha['total']}")
    print()

    # Show sample predictions
    print("Sample predictions:")
    results = recognizer_rha.batch_recognize(test_data[:5])
    for antigen, result in zip(test_data[:5], results):
        certain = "✓" if result.is_certain else "?"
        correct = "✓" if result.predicted_class == antigen.class_label else "✗"
        print(f"  {correct} True: {antigen.class_label}, "
              f"Pred: {result.predicted_class} ({result.confidence:.4f}) {certain}")


def cross_validation_demo(antigens):
    """Full 10-fold cross-validation to reproduce paper results"""

    validator = CrossValidator(n_folds=10, random_seed=42)
    results = validator.cross_validate(antigens)

    # Print detailed comparison
    validator.compare_strategies(antigens)


if __name__ == "__main__":
    main()

```
