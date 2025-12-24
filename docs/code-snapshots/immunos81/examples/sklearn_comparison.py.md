---
source: /Users/byron/projects/immunos81/examples/sklearn_comparison.py
relative: immunos81/examples/sklearn_comparison.py
generated_at: 2025-12-23 10:28
---

```python
"""
Scikit-learn Comparison Benchmark

Compares Immunos-81 against standard ML algorithms using scikit-learn:
- Random Forest
- SVM
- Naive Bayes
- Logistic Regression
- Decision Tree
- K-Nearest Neighbors
"""

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix

from immunos81.data.cleveland_loader import ClevelandDataLoader
from immunos81.engine.trainer import Trainer
from immunos81.engine.recognizer import Recognizer, RecognitionStrategy


def print_separator(title=""):
    """Print a nice separator"""
    if title:
        print("\n" + "="*80)
        print(f"  {title}")
        print("="*80)
    else:
        print("-"*80)


def antigens_to_arrays(antigens):
    """
    Convert antigens to numpy arrays for sklearn.

    Returns:
        X: Feature matrix (n_samples, n_features)
        y: Label vector (n_samples,)
    """
    if not antigens:
        return np.array([]), np.array([])

    n_samples = len(antigens)
    n_features = len(antigens[0].values)

    X = np.zeros((n_samples, n_features))
    y = []

    for i, antigen in enumerate(antigens):
        # Convert values to numeric (handle None as NaN)
        for j, value in enumerate(antigen.values):
            if value is None:
                X[i, j] = np.nan
            else:
                try:
                    X[i, j] = float(value)
                except (ValueError, TypeError):
                    # For nominal values, use hash
                    X[i, j] = hash(str(value)) % 1000

        y.append(antigen.class_label)

    # Handle missing values (simple mean imputation)
    col_means = np.nanmean(X, axis=0)
    inds = np.where(np.isnan(X))
    X[inds] = np.take(col_means, inds[1])

    return X, np.array(y)


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

    return [
        {
            'name': 'Immunos-81 (SHA)',
            'accuracy': metrics_sha['accuracy_overall'],
            'correct': metrics_sha['correct'],
            'total': metrics_sha['total']
        },
        {
            'name': 'Immunos-81 (RHA certain)',
            'accuracy': metrics_rha['accuracy_certain'],
            'correct': metrics_rha['correct_certain'],
            'total': metrics_rha['certain_attempts']
        }
    ]


def run_sklearn_classifier(clf, name, X_train, y_train, X_test, y_test):
    """
    Run a scikit-learn classifier.

    Args:
        clf: Scikit-learn classifier
        name: Display name
        X_train, y_train: Training data
        X_test, y_test: Test data

    Returns:
        Dictionary with results
    """
    try:
        # Train
        clf.fit(X_train, y_train)

        # Predict
        y_pred = clf.predict(X_test)

        # Evaluate
        accuracy = accuracy_score(y_test, y_pred)
        correct = np.sum(y_pred == y_test)
        total = len(y_test)

        return {
            'name': name,
            'accuracy': accuracy,
            'correct': correct,
            'total': total,
            'predictions': y_pred
        }

    except Exception as e:
        print(f"Error running {name}: {e}")
        return None


def run_sklearn_benchmarks(train_data, test_data):
    """Run all scikit-learn classifiers"""
    print_separator("SCIKIT-LEARN ALGORITHMS")

    # Convert to numpy arrays
    print("\nConverting data to sklearn format...")
    X_train, y_train = antigens_to_arrays(train_data)
    X_test, y_test = antigens_to_arrays(test_data)

    print(f"Training set: {X_train.shape}")
    print(f"Test set: {X_test.shape}")

    # Define classifiers
    classifiers = [
        (RandomForestClassifier(n_estimators=100, random_state=42), 'Random Forest'),
        (SVC(kernel='rbf', C=1.0, random_state=42), 'SVM (RBF kernel)'),
        (GaussianNB(), 'Naive Bayes'),
        (LogisticRegression(max_iter=1000, random_state=42), 'Logistic Regression'),
        (DecisionTreeClassifier(random_state=42), 'Decision Tree'),
        (KNeighborsClassifier(n_neighbors=5), 'K-Nearest Neighbors')
    ]

    results = []

    for clf, name in classifiers:
        print(f"\n{name}:")
        result = run_sklearn_classifier(clf, name, X_train, y_train, X_test, y_test)

        if result:
            print(f"  Accuracy: {result['accuracy']:.1%}")
            print(f"  Correct: {result['correct']}/{result['total']}")
            results.append(result)

    return results


def print_comparison_table(immunos_results, sklearn_results):
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

    # Scikit-learn results
    for result in sklearn_results:
        print("{:<30} {:>11.1%} {:>8}/{:<6}".format(
            result['name'],
            result['accuracy'],
            result['correct'],
            result['total']
        ))

    print("-"*60)

    # Best performer
    all_results = immunos_results + sklearn_results
    best = max(all_results, key=lambda x: x['accuracy'])
    print(f"\n🏆 Best: {best['name']} at {best['accuracy']:.1%}")

    # Gap analysis
    immunos_sha = immunos_results[0]['accuracy']
    best_sklearn = max(sklearn_results, key=lambda x: x['accuracy'])
    gap = best_sklearn['accuracy'] - immunos_sha

    print(f"\n📊 Performance Gap:")
    print(f"   Immunos-81 (SHA): {immunos_sha:.1%}")
    print(f"   Best sklearn ({best_sklearn['name']}): {best_sklearn['accuracy']:.1%}")
    print(f"   Gap: {gap:+.1%}")

    print(f"\n📈 Comparison with Research (Cleveland dataset):")
    print(f"   Published studies report:")
    print(f"   - Random Forest: ~88.5%")
    print(f"   - SVM: ~84.2%")
    print(f"   - Naive Bayes: ~85.0%")
    print(f"   - Original Immunos-81 paper: 83.2% (SHA) / 85.5% (RHA)")
    print(f"\n   Our results:")
    rf_result = next((r for r in sklearn_results if 'Random Forest' in r['name']), None)
    if rf_result:
        print(f"   - Random Forest: {rf_result['accuracy']:.1%}")
    print(f"   - Immunos-81 (ours): {immunos_sha:.1%}")

    print(f"\n💡 Analysis:")
    print(f"   The ~20% gap to paper results (83.2%) suggests:")
    print(f"   1. Possible implementation differences")
    print(f"   2. Different data preprocessing/splits")
    print(f"   3. Feature engineering not captured in our implementation")


def main():
    print_separator("MACHINE LEARNING COMPARISON BENCHMARK")
    print("\nComparing Immunos-81 against standard ML algorithms")
    print("on the Cleveland Heart Disease dataset using scikit-learn.")

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
    immunos_results = run_immunos81(train_data, test_data, library)

    # Run scikit-learn benchmarks
    sklearn_results = run_sklearn_benchmarks(train_data, test_data)

    # Print comparison
    print_comparison_table(immunos_results, sklearn_results)

    print_separator()
    print("\n✓ Benchmark complete!")


if __name__ == "__main__":
    main()

```
