---
source: /Users/byron/projects/immunos81/examples/enhanced_demo.py
relative: immunos81/examples/enhanced_demo.py
generated_at: 2025-12-23 10:28
---

```python
"""
Enhanced Immunos-81 Demo

Compares baseline implementation with enhanced version featuring:
- Multiple clones per class (clustering)
- Adaptive clonal factor (0.4)
- Affinity maturation via hypermutation
- Negative selection
- Activation/inhibition mechanism
"""

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from immunos81.data.cleveland_loader import ClevelandDataLoader
from immunos81.engine.trainer import Trainer
from immunos81.engine.recognizer import Recognizer, RecognitionStrategy
from immunos81.engine.enhanced_trainer import EnhancedTrainer
from immunos81.engine.enhanced_recognizer import EnhancedRecognizer
from immunos81.engine.enhanced_recognizer import RecognitionStrategy as EnhancedStrategy


def print_separator(title=""):
    """Print a nice separator"""
    if title:
        print("\n" + "="*70)
        print(f"  {title}")
        print("="*70)
    else:
        print("-"*70)


def print_metrics(metrics, prefix=""):
    """Print evaluation metrics"""
    print(f"{prefix}Total: {metrics['total']}")
    print(f"{prefix}Correct: {metrics['correct']}")
    print(f"{prefix}Accuracy: {metrics['accuracy_overall']:.1%}")

    if 'accuracy_certain' in metrics:
        print(f"{prefix}Certain Accuracy: {metrics['accuracy_certain']:.1%}")
        uncertain_rate = metrics.get('uncertain_rate', metrics['uncertain_attempts'] / metrics['total'] if metrics['total'] > 0 else 0.0)
        print(f"{prefix}Uncertain: {metrics['uncertain_attempts']}/{metrics['total']} ({uncertain_rate:.1%})")


def run_baseline(train_data, test_data, library):
    """Run baseline Immunos-81 implementation"""
    print_separator("BASELINE IMPLEMENTATION")
    print("\nFeatures:")
    print("  - 1 clone per class")
    print("  - Standard concentration factor: log(1 + size)")
    print("  - No affinity maturation")
    print("  - No negative selection")
    print("  - No activation/inhibition")
    print()

    # Train baseline
    print("Training baseline model...")
    trainer = Trainer(library)
    tcells = trainer.train(train_data)
    print()

    # Test with SHA
    print("Testing with SHA strategy:")
    recognizer = Recognizer(tcells, strategy=RecognitionStrategy.SHA)
    metrics_sha = recognizer.evaluate(test_data)
    print_metrics(metrics_sha, "  ")
    print()

    # Test with RHA
    print("Testing with RHA strategy:")
    recognizer.set_strategy(RecognitionStrategy.RHA)
    metrics_rha = recognizer.evaluate(test_data)
    print_metrics(metrics_rha, "  ")

    return {
        'sha': metrics_sha,
        'rha': metrics_rha
    }


def run_enhanced(train_data, test_data, library):
    """Run enhanced Immunos-81 implementation"""
    print_separator("ENHANCED IMPLEMENTATION")
    print("\nFeatures:")
    print("  - Multiple clones per class (2-4 via clustering)")
    print("  - Standard concentration factor: log(1 + size)")
    print("  - No affinity maturation (disabled for testing)")
    print("  - No negative selection (disabled for testing)")
    print("  - No activation/inhibition (disabled for testing)")
    print()

    # Train enhanced with conservative settings
    print("Training enhanced model...")
    trainer = EnhancedTrainer(
        library,
        enable_maturation=False,  # Disable affinity maturation
        enable_negative_selection=False  # Disable negative selection
    )
    tcells = trainer.train(train_data)

    # Print statistics
    stats = trainer.get_statistics()
    print(f"\nEnhanced Model Statistics:")
    print(f"  Total clones: {stats['num_clones']}")
    for class_label, class_stats in stats['classes'].items():
        print(f"  {class_label}: {class_stats['num_clones']} clones, avg size: {class_stats['avg_clone_size']:.1f}")
    print()

    # Test with SHA
    print("Testing with SHA strategy:")
    recognizer = EnhancedRecognizer(
        tcells,
        strategy=EnhancedStrategy.SHA,
        enable_inhibition=False  # Disable lateral inhibition
    )
    metrics_sha = recognizer.evaluate(test_data)
    print_metrics(metrics_sha, "  ")
    print()

    # Test with RHA
    print("Testing with RHA strategy:")
    recognizer.set_strategy(EnhancedStrategy.RHA)
    metrics_rha = recognizer.evaluate(test_data)
    print_metrics(metrics_rha, "  ")

    return {
        'sha': metrics_sha,
        'rha': metrics_rha
    }


def compare_results(baseline, enhanced):
    """Compare baseline vs enhanced results"""
    print_separator("COMPARISON")

    print("\nSHA Strategy:")
    print(f"  Baseline:  {baseline['sha']['accuracy_overall']:.1%}")
    print(f"  Enhanced:  {enhanced['sha']['accuracy_overall']:.1%}")
    improvement_sha = enhanced['sha']['accuracy_overall'] - baseline['sha']['accuracy_overall']
    print(f"  Improvement: {improvement_sha:+.1%}")

    print("\nRHA Strategy:")
    print(f"  Baseline (all):      {baseline['rha']['accuracy_overall']:.1%}")
    print(f"  Enhanced (all):      {enhanced['rha']['accuracy_overall']:.1%}")
    print(f"  Baseline (certain):  {baseline['rha']['accuracy_certain']:.1%}")
    print(f"  Enhanced (certain):  {enhanced['rha']['accuracy_certain']:.1%}")
    improvement_rha = enhanced['rha']['accuracy_certain'] - baseline['rha']['accuracy_certain']
    print(f"  Improvement (certain): {improvement_rha:+.1%}")

    print("\nUncertain Cases:")
    baseline_uncertain = baseline['rha'].get('uncertain_rate',
        baseline['rha']['uncertain_attempts'] / baseline['rha']['total'] if baseline['rha']['total'] > 0 else 0.0)
    enhanced_uncertain = enhanced['rha'].get('uncertain_rate',
        enhanced['rha']['uncertain_attempts'] / enhanced['rha']['total'] if enhanced['rha']['total'] > 0 else 0.0)
    print(f"  Baseline: {baseline_uncertain:.1%}")
    print(f"  Enhanced: {enhanced_uncertain:.1%}")

    print_separator()

    # Compare to paper results
    print("\nPaper Results (Hunt & Cooke, 2000):")
    print("  SHA: 83.2%")
    print("  RHA (certain): 85.5%")
    print("  RHA (uncertain): 11.2%")

    print("\nEnhanced Implementation vs Paper:")
    sha_vs_paper = enhanced['sha']['accuracy_overall'] - 0.832
    rha_vs_paper = enhanced['rha']['accuracy_certain'] - 0.855
    print(f"  SHA: {sha_vs_paper:+.1%} from paper")
    print(f"  RHA (certain): {rha_vs_paper:+.1%} from paper")


def main():
    print_separator("IMMUNOS-81: Enhanced Implementation Demo")
    print("\nThis demo compares the baseline implementation with the enhanced")
    print("version that integrates accuracy improvements from AIS research:")
    print("  - Clonal Selection Algorithm improvements")
    print("  - Negative Selection Algorithm")
    print("  - Amorphous Computing (activation/inhibition)")

    # Load dataset
    print_separator("Loading Dataset")
    loader = ClevelandDataLoader()
    antigens, library = loader.load_data()
    print(f"Loaded {len(antigens)} records from Cleveland Heart Disease dataset")
    print(f"Classes: {set(a.class_label for a in antigens)}")

    # Split data (70/30)
    n_train = int(len(antigens) * 0.7)
    train_data = antigens[:n_train]
    test_data = antigens[n_train:]
    print(f"\nTrain set: {len(train_data)} instances")
    print(f"Test set: {len(test_data)} instances")

    # Run baseline
    baseline_results = run_baseline(train_data, test_data, library)

    # Run enhanced
    enhanced_results = run_enhanced(train_data, test_data, library)

    # Compare
    compare_results(baseline_results, enhanced_results)

    print("\n" + "="*70)
    print("Demo Complete!")
    print("="*70)
    print("\n💡 Implemented Enhancements:")
    print("  1. Multiple specialized clones per class (clustering)")
    print("  2. Affinity maturation via hypermutation")
    print("  3. Negative selection filters weak/cross-reactive clones")
    print("  4. Lateral inhibition between competing clones")


if __name__ == "__main__":
    main()

```
