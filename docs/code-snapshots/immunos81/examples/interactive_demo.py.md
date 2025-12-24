---
source: /Users/byron/projects/immunos81/examples/interactive_demo.py
relative: immunos81/examples/interactive_demo.py
generated_at: 2025-12-23 10:28
---

```python
"""
Interactive Visualization Demo

Demonstrates interactive Plotly visualizations with HTML output.
"""

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from immunos81.data.cleveland_loader import ClevelandDataLoader
from immunos81.engine.parallel_trainer import ParallelTrainer
from immunos81.engine.parallel_recognizer import ParallelRecognizer, RecognitionStrategy
from immunos81.utils.interactive_visualizer import InteractiveVisualizer


def main():
    print("="*60)
    print("IMMUNOS-81: Interactive Visualization Demo")
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

    # Train with parallel trainer
    print("Training with parallel engine...")
    trainer = ParallelTrainer(library, n_workers=4)
    tcells = trainer.train(train_data)
    print()

    # Test with both strategies
    print("Testing with SHA strategy...")
    recognizer_sha = ParallelRecognizer(tcells, strategy=RecognitionStrategy.SHA, n_workers=4)
    results_sha = recognizer_sha.batch_recognize(test_data)
    metrics_sha = recognizer_sha.evaluate(test_data)
    print(f"  SHA Accuracy: {metrics_sha['accuracy_overall']:.1%}\n")

    print("Testing with RHA strategy...")
    recognizer_rha = ParallelRecognizer(tcells, strategy=RecognitionStrategy.RHA, n_workers=4)
    results_rha = recognizer_rha.batch_recognize(test_data)
    metrics_rha = recognizer_rha.evaluate(test_data)
    print(f"  RHA Accuracy (all): {metrics_rha['accuracy_overall']:.1%}")
    print(f"  RHA Accuracy (certain): {metrics_rha['accuracy_certain']:.1%}")
    print(f"  RHA Uncertain: {metrics_rha['uncertain_attempts']}/{metrics_rha['total']}\n")

    # Create output directory
    output_dir = "/Users/byron/projects/immunos81/examples/interactive_plots"
    os.makedirs(output_dir, exist_ok=True)

    # Initialize visualizer
    viz = InteractiveVisualizer()

    # Generate all interactive visualizations
    print("Generating interactive HTML visualizations...")
    print()

    print("  1. Interactive Confusion Matrix (SHA)...")
    fig1 = viz.plot_interactive_confusion_matrix(
        test_data, results_sha,
        save_path=f"{output_dir}/confusion_matrix_sha.html"
    )
    print(f"     Saved to: {output_dir}/confusion_matrix_sha.html")

    print("  2. Interactive Confusion Matrix (RHA)...")
    fig2 = viz.plot_interactive_confusion_matrix(
        test_data, results_rha,
        save_path=f"{output_dir}/confusion_matrix_rha.html"
    )
    print(f"     Saved to: {output_dir}/confusion_matrix_rha.html")

    print("  3. Interactive Decision Confidence (SHA)...")
    fig3 = viz.plot_interactive_confidence(
        results_sha, test_data,
        save_path=f"{output_dir}/confidence_sha.html"
    )
    print(f"     Saved to: {output_dir}/confidence_sha.html")

    print("  4. Interactive Decision Confidence (RHA)...")
    fig4 = viz.plot_interactive_confidence(
        results_rha, test_data,
        save_path=f"{output_dir}/confidence_rha.html"
    )
    print(f"     Saved to: {output_dir}/confidence_rha.html")

    print("  5. Interactive Variable Importance...")
    fig5 = viz.plot_interactive_variable_importance(
        tcells, test_data, top_n=13,  # All 13 features
        save_path=f"{output_dir}/variable_importance.html"
    )
    print(f"     Saved to: {output_dir}/variable_importance.html")

    print()
    print("="*60)
    print("✓ All interactive visualizations generated!")
    print("="*60)
    print()
    print(f"📁 Output directory: {output_dir}/")
    print()
    print("📊 Generated files:")
    for file in sorted(os.listdir(output_dir)):
        if file.endswith('.html'):
            filepath = os.path.join(output_dir, file)
            size_kb = os.path.getsize(filepath) / 1024
            print(f"  - {file} ({size_kb:.0f} KB)")
    print()
    print("💡 Open any HTML file in your browser to interact with the visualizations!")
    print("   Features: hover tooltips, zoom, pan, export to PNG")


if __name__ == "__main__":
    main()

```
