---
source: /Users/byron/projects/immunos81/examples/weka_comparison.py
relative: immunos81/examples/weka_comparison.py
generated_at: 2025-12-23 10:28
---

```python
"""
Weka Comparison Benchmark

Compares Immunos-81 against standard ML algorithms using Weka:
- Random Forest
- SVM (SMO)
- Naive Bayes
- Logistic Regression
- Decision Tree (J48)
"""

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from immunos81.data.cleveland_loader import ClevelandDataLoader
from immunos81.engine.trainer import Trainer
from immunos81.engine.recognizer import Recognizer, RecognitionStrategy
from immunos81.utils.weka_exporter import export_train_test_split
import weka.core.jvm as jvm
from weka.core.converters import Loader
from weka.classifiers import Classifier, Evaluation
from weka.core.classes import Random


def print_separator(title=""):
    """Print a nice separator"""
    if title:
        print("\n" + "="*80)
        print(f"  {title}")
        print("="*80)
    else:
        print("-"*80)


def run_immunos81(train_data, test_data, library):
    """Run Immunos-81 baseline"""
    print_separator("IMMUNOS-81 (Baseline)")

    # Train
    print("Training Immunos-81...")
    trainer = Trainer(library)
    tcells = trainer.train(train_data)

    # Test with SHA
    print("\nTesting with SHA strategy:")
    recognizer = Recognizer(tcells, strategy=RecognitionStrategy.SHA)
    metrics_sha = recognizer.evaluate(test_data)

    print(f"  Accuracy: {metrics_sha['accuracy_overall']:.1%}")
    print(f"  Correct: {metrics_sha['correct']}/{metrics_sha['total']}")

    # Test with RHA
    print("\nTesting with RHA strategy:")
    recognizer.set_strategy(RecognitionStrategy.RHA)
    metrics_rha = recognizer.evaluate(test_data)

    print(f"  Accuracy (overall): {metrics_rha['accuracy_overall']:.1%}")
    print(f"  Accuracy (certain): {metrics_rha['accuracy_certain']:.1%}")
    print(f"  Uncertain: {metrics_rha['uncertain_attempts']}/{metrics_rha['total']}")

    return {
        'name': 'Immunos-81 (SHA)',
        'accuracy': metrics_sha['accuracy_overall'],
        'correct': metrics_sha['correct'],
        'total': metrics_sha['total']
    }, {
        'name': 'Immunos-81 (RHA)',
        'accuracy': metrics_rha['accuracy_certain'],
        'correct': metrics_rha['correct_certain'],
        'total': metrics_rha['certain_attempts']
    }


def run_weka_classifier(train_path, test_path, classifier_name, classifier_options=None):
    """
    Run a Weka classifier on the exported ARFF data.

    Args:
        train_path: Path to training ARFF file
        test_path: Path to test ARFF file
        classifier_name: Weka classifier class name
        classifier_options: Optional list of classifier options

    Returns:
        Dictionary with results
    """
    try:
        # Load data
        loader = Loader(classname="weka.core.converters.ArffLoader")
        train_data = loader.load_file(train_path)
        train_data.class_is_last()

        test_data = loader.load_file(test_path)
        test_data.class_is_last()

        # Build classifier
        classifier = Classifier(classname=classifier_name, options=classifier_options)
        classifier.build_classifier(train_data)

        # Evaluate
        evaluation = Evaluation(train_data)
        evaluation.test_model(classifier, test_data)

        return {
            'name': classifier_name.split('.')[-1],
            'accuracy': evaluation.percent_correct / 100.0,
            'correct': evaluation.correct,
            'total': test_data.num_instances,
            'error_rate': evaluation.error_rate,
            'confusion_matrix': evaluation.confusion_matrix
        }

    except Exception as e:
        print(f"Error running {classifier_name}: {e}")
        return None


def run_weka_benchmarks(train_path, test_path):
    """Run all Weka classifiers"""
    print_separator("WEKA ALGORITHMS")

    algorithms = [
        {
            'name': 'Random Forest',
            'classname': 'weka.classifiers.trees.RandomForest',
            'options': ['-I', '100', '-K', '0', '-S', '1']
        },
        {
            'name': 'SVM (SMO)',
            'classname': 'weka.classifiers.functions.SMO',
            'options': ['-C', '1.0', '-L', '0.001']
        },
        {
            'name': 'Naive Bayes',
            'classname': 'weka.classifiers.bayes.NaiveBayes',
            'options': []
        },
        {
            'name': 'Logistic Regression',
            'classname': 'weka.classifiers.functions.Logistic',
            'options': []
        },
        {
            'name': 'Decision Tree (J48)',
            'classname': 'weka.classifiers.trees.J48',
            'options': ['-C', '0.25', '-M', '2']
        }
    ]

    results = []

    for algo in algorithms:
        print(f"\n{algo['name']}:")
        result = run_weka_classifier(
            train_path,
            test_path,
            algo['classname'],
            algo['options']
        )

        if result:
            print(f"  Accuracy: {result['accuracy']:.1%}")
            print(f"  Correct: {int(result['correct'])}/{int(result['total'])}")
            results.append(result)
        else:
            print(f"  Failed to run")

    return results


def print_comparison_table(immunos_results, weka_results):
    """Print comparison table"""
    print_separator("COMPARISON RESULTS")

    print("\n{:<30} {:>12} {:>15}".format("Algorithm", "Accuracy", "Correct/Total"))
    print("-"*60)

    # Immunos-81 results
    for result in immunos_results:
        print("{:<30} {:>11.1%} {:>8}/{:<6}".format(
            result['name'],
            result['accuracy'],
            result['correct'],
            result['total']
        ))

    print("-"*60)

    # Weka results
    for result in weka_results:
        print("{:<30} {:>11.1%} {:>8}/{:<6}".format(
            result['name'],
            result['accuracy'],
            int(result['correct']),
            int(result['total'])
        ))

    print("-"*60)

    # Best performer
    all_results = immunos_results + weka_results
    best = max(all_results, key=lambda x: x['accuracy'])
    print(f"\n🏆 Best: {best['name']} at {best['accuracy']:.1%}")

    # Gap analysis
    immunos_sha = immunos_results[0]['accuracy']
    best_weka = max(weka_results, key=lambda x: x['accuracy'])
    gap = best_weka['accuracy'] - immunos_sha

    print(f"\n📊 Performance Gap:")
    print(f"   Immunos-81 (SHA): {immunos_sha:.1%}")
    print(f"   Best Weka ({best_weka['name']}): {best_weka['accuracy']:.1%}")
    print(f"   Gap: {gap:.1%}")


def main():
    print_separator("WEKA COMPARISON BENCHMARK")
    print("\nComparing Immunos-81 against standard ML algorithms on")
    print("the Cleveland Heart Disease dataset using Weka.")

    # Load data
    print_separator("Loading Dataset")
    loader = ClevelandDataLoader()
    antigens, library = loader.load_data()
    print(f"Loaded {len(antigens)} records")

    # Split data (same 70/30 split as other demos)
    n_train = int(len(antigens) * 0.7)
    train_data = antigens[:n_train]
    test_data = antigens[n_train:]

    print(f"Train set: {len(train_data)} instances")
    print(f"Test set: {len(test_data)} instances")

    # Run Immunos-81
    immunos_sha, immunos_rha = run_immunos81(train_data, test_data, library)
    immunos_results = [immunos_sha, immunos_rha]

    # Export to ARFF
    print_separator("Exporting to ARFF")
    output_dir = "/Users/byron/projects/immunos81/examples/weka_data"
    os.makedirs(output_dir, exist_ok=True)

    train_path = os.path.join(output_dir, "cleveland_train.arff")
    test_path = os.path.join(output_dir, "cleveland_test.arff")

    export_train_test_split(train_data, test_data, library, train_path, test_path)

    # Start JVM for Weka
    print_separator("Starting Weka")
    print("Initializing Java VM for Weka...")
    jvm.start()

    try:
        # Run Weka benchmarks
        weka_results = run_weka_benchmarks(train_path, test_path)

        # Print comparison
        print_comparison_table(immunos_results, weka_results)

        print_separator()
        print("\n✓ Benchmark complete!")
        print(f"\n📁 ARFF files saved to: {output_dir}/")
        print("   - cleveland_train.arff")
        print("   - cleveland_test.arff")

    finally:
        # Stop JVM
        jvm.stop()


if __name__ == "__main__":
    main()

```
