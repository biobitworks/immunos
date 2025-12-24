---
source: /Users/byron/projects/immunos81/examples/visualization_demo.py
relative: immunos81/examples/visualization_demo.py
generated_at: 2025-12-23 10:28
---

```python
"""
Immunos-81 Visualization Demo

Demonstrates all visualization capabilities of the Immunos-81 system.
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from immunos81.data.cleveland_loader import ClevelandDataLoader
from immunos81.engine.trainer import Trainer
from immunos81.engine.recognizer import Recognizer, RecognitionStrategy
from immunos81.utils.cross_validator import CrossValidator
from immunos81.utils.visualizer import Immunos81Visualizer


def main():
    """Run comprehensive visualization demo"""

    print("="*60)
    print("IMMUNOS-81: Visualization Demo")
    print("="*60)
    print()

    # Load dataset
    print("Loading Cleveland Heart Disease dataset...")
    loader = ClevelandDataLoader()
    antigens, library = loader.load_data()
    print(f"Loaded {len(antigens)} patient records\n")

    # Initialize visualizer
    viz = Immunos81Visualizer()

    # Demo 1: Quick train/test split with visualizations
    print("Demo 1: Train/Test Split with Visualizations")
    print("-" * 60)
    quick_visualization_demo(antigens, viz)

    # Demo 2: Cross-validation with comprehensive plots
    print("\n\nDemo 2: Cross-Validation Results Visualization")
    print("-" * 60)
    cross_validation_visualization_demo(antigens, viz)


def quick_visualization_demo(antigens, viz):
    """
    Quick demo with train/test split showing all visualization types.

    Args:
        antigens: List of all antigens
        viz: Visualizer instance
    """
    # Split data 70/30
    n_train = int(len(antigens) * 0.7)
    train_data = antigens[:n_train]
    test_data = antigens[n_train:]

    print(f"Training on {len(train_data)} samples...")

    # Train
    trainer = Trainer(train_data[0].library)
    tcells = trainer.train(train_data)
    print()

    # Test with both strategies
    print("Testing with SHA and RHA strategies...\n")

    recognizer_sha = Recognizer(tcells, strategy=RecognitionStrategy.SHA)
    results_sha = recognizer_sha.batch_recognize(test_data)
    metrics_sha = recognizer_sha.evaluate(test_data)

    recognizer_rha = Recognizer(tcells, strategy=RecognitionStrategy.RHA)
    results_rha = recognizer_rha.batch_recognize(test_data)
    metrics_rha = recognizer_rha.evaluate(test_data)

    # Visualization 1: Avidity Distribution
    print("Generating Avidity Distribution plot...")
    viz.plot_avidity_distribution(tcells, test_data)

    # Visualization 2: Confusion Matrix (SHA)
    print("Generating Confusion Matrix (SHA)...")
    viz.plot_confusion_matrix(test_data, results_sha)

    # Visualization 3: Confusion Matrix (RHA)
    print("Generating Confusion Matrix (RHA)...")
    viz.plot_confusion_matrix(test_data, results_rha)

    # Visualization 4: Decision Confidence (SHA)
    print("Generating Decision Confidence plot (SHA)...")
    viz.plot_decision_confidence(results_sha, test_data)

    # Visualization 5: Decision Confidence (RHA)
    print("Generating Decision Confidence plot (RHA)...")
    viz.plot_decision_confidence(results_rha, test_data)

    # Visualization 6: Variable Importance
    print("Generating Variable Importance plot...")
    viz.plot_variable_importance(tcells, test_data, top_n=10)

    print("\nAll visualizations generated successfully!")


def cross_validation_visualization_demo(antigens, viz):
    """
    Cross-validation demo with result visualization.

    Args:
        antigens: List of all antigens
        viz: Visualizer instance
    """
    print("Running 10-fold cross-validation...")
    print("(This will take a few minutes)\n")

    # Run cross-validation
    validator = CrossValidator(n_folds=10, random_seed=42)
    results = validator.cross_validate(antigens)

    print("\nGenerating Cross-Validation Results visualization...")
    viz.plot_cross_validation_results(results)

    print("\nCross-validation visualization complete!")

    # Print summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}\n")

    sha_results = results["strategies"]["sha"]
    rha_results = results["strategies"]["rha"]

    print(f"SHA Strategy:")
    print(f"  Accuracy: {sha_results['overall_accuracy']:.1%}")
    print(f"  Correct: {sha_results['total_correct']}/{sha_results['total_samples']}\n")

    print(f"RHA Strategy:")
    print(f"  Accuracy (all): {rha_results['overall_accuracy']:.1%}")
    print(f"  Accuracy (certain): {rha_results['accuracy_certain']:.1%}")
    print(f"  Uncertain: {rha_results['total_uncertain']}/{rha_results['total_samples']}")
    print()


if __name__ == "__main__":
    main()

```
