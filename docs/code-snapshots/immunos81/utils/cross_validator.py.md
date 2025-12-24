---
source: /Users/byron/projects/immunos81/utils/cross_validator.py
relative: immunos81/utils/cross_validator.py
generated_at: 2025-12-23 10:28
---

```python
"""
Cross-Validation Utilities for Immunos-81

Implements k-fold cross-validation to evaluate and compare SHA and RHA strategies.
"""

from typing import List, Dict, Tuple
import numpy as np
from ..core.antigen import Antigen
from ..engine.trainer import Trainer
from ..engine.recognizer import Recognizer, RecognitionStrategy


class CrossValidator:
    """
    Performs k-fold cross-validation on Immunos-81.

    This class implements stratified k-fold cross-validation to evaluate
    the performance of different recognition strategies.
    """

    def __init__(self, n_folds: int = 10, random_seed: int = 42):
        """
        Initialize cross-validator.

        Args:
            n_folds: Number of folds for cross-validation
            random_seed: Random seed for reproducibility
        """
        self.n_folds = n_folds
        self.random_seed = random_seed

    def split_data(self, antigens: List[Antigen]) -> List[Tuple[List[Antigen], List[Antigen]]]:
        """
        Split data into k folds (stratified by class).

        Args:
            antigens: List of all antigens

        Returns:
            List of (train, test) tuples for each fold
        """
        # Set random seed for reproducibility
        np.random.seed(self.random_seed)

        # Group antigens by class
        class_groups = {}
        for antigen in antigens:
            if antigen.class_label not in class_groups:
                class_groups[antigen.class_label] = []
            class_groups[antigen.class_label].append(antigen)

        # Shuffle each class group
        for class_label in class_groups:
            np.random.shuffle(class_groups[class_label])

        # Create folds
        folds = []
        for fold_idx in range(self.n_folds):
            test_fold = []
            train_fold = []

            for class_label, class_antigens in class_groups.items():
                n_samples = len(class_antigens)
                fold_size = n_samples // self.n_folds

                start_idx = fold_idx * fold_size
                if fold_idx == self.n_folds - 1:
                    # Last fold gets remainder
                    end_idx = n_samples
                else:
                    end_idx = start_idx + fold_size

                # Test samples for this fold
                test_samples = class_antigens[start_idx:end_idx]
                test_fold.extend(test_samples)

                # Train samples (all others)
                train_samples = class_antigens[:start_idx] + class_antigens[end_idx:]
                train_fold.extend(train_samples)

            folds.append((train_fold, test_fold))

        return folds

    def cross_validate(self, antigens: List[Antigen],
                      strategies: List[RecognitionStrategy] = None) -> Dict[str, any]:
        """
        Perform k-fold cross-validation.

        Args:
            antigens: List of all antigens
            strategies: List of recognition strategies to evaluate

        Returns:
            Dictionary containing cross-validation results
        """
        if strategies is None:
            strategies = [RecognitionStrategy.SHA, RecognitionStrategy.RHA]

        folds = self.split_data(antigens)

        results = {
            "n_folds": self.n_folds,
            "n_samples": len(antigens),
            "strategies": {}
        }

        for strategy in strategies:
            print(f"\n{'='*60}")
            print(f"Evaluating strategy: {strategy.value.upper()}")
            print(f"{'='*60}\n")

            strategy_results = self._evaluate_strategy(folds, strategy)
            results["strategies"][strategy.value] = strategy_results

            # Print summary
            print(f"\n{strategy.value.upper()} Summary:")
            print(f"  Overall Accuracy: {strategy_results['overall_accuracy']:.1%}")
            if strategy == RecognitionStrategy.RHA:
                print(f"  Accuracy (certain only): {strategy_results['accuracy_certain']:.1%}")
                print(f"  Uncertain cases: {strategy_results['total_uncertain']}/{strategy_results['total_samples']}")

        return results

    def _evaluate_strategy(self, folds: List[Tuple[List[Antigen], List[Antigen]]],
                          strategy: RecognitionStrategy) -> Dict[str, any]:
        """Evaluate a single strategy across all folds"""

        fold_results = []
        total_correct = 0
        total_attempted = 0
        total_samples = 0
        total_certain = 0
        total_correct_certain = 0

        for fold_idx, (train_data, test_data) in enumerate(folds):
            print(f"Fold {fold_idx + 1}/{len(folds)}:", end=" ")

            # Train on this fold
            trainer = Trainer(train_data[0].library)
            tcells = trainer.train(train_data)

            # Recognize test data
            recognizer = Recognizer(tcells, strategy=strategy)
            metrics = recognizer.evaluate(test_data)

            fold_results.append(metrics)

            # Accumulate totals
            total_samples += metrics["total"]
            total_attempted += metrics["attempted"]
            total_correct += metrics["correct"]

            print(f"{metrics['correct']}/{metrics['attempted']} correct", end="")

            if strategy == RecognitionStrategy.RHA:
                total_certain += metrics["certain_attempts"]
                total_correct_certain += metrics["correct_certain"]
                print(f" ({metrics['certain_attempts']} certain)", end="")

            print()

        # Calculate overall metrics
        overall_results = {
            "fold_results": fold_results,
            "total_samples": total_samples,
            "total_attempted": total_attempted,
            "total_correct": total_correct,
            "overall_accuracy": total_correct / total_samples if total_samples > 0 else 0.0,
            "attempted_accuracy": total_correct / total_attempted if total_attempted > 0 else 0.0
        }

        if strategy == RecognitionStrategy.RHA:
            overall_results.update({
                "total_certain": total_certain,
                "total_correct_certain": total_correct_certain,
                "total_uncertain": total_samples - total_certain,
                "accuracy_certain": total_correct_certain / total_certain if total_certain > 0 else 0.0
            })

        return overall_results

    def compare_strategies(self, antigens: List[Antigen]) -> None:
        """
        Run cross-validation and print comparison between SHA and RHA.

        Args:
            antigens: List of all antigens
        """
        results = self.cross_validate(antigens)

        print(f"\n{'='*60}")
        print("STRATEGY COMPARISON")
        print(f"{'='*60}\n")

        sha_results = results["strategies"]["sha"]
        rha_results = results["strategies"]["rha"]

        print(f"Dataset: {results['n_samples']} samples, {results['n_folds']}-fold CV\n")

        print(f"SHA (Simple Highest Avidity):")
        print(f"  Correct: {sha_results['total_correct']}/{sha_results['total_samples']}")
        print(f"  Accuracy: {sha_results['overall_accuracy']:.1%}\n")

        print(f"RHA (Relative Highest Avidity, 5% threshold):")
        print(f"  Correct (all): {rha_results['total_correct']}/{rha_results['total_samples']}")
        print(f"  Accuracy (all): {rha_results['overall_accuracy']:.1%}")
        print(f"  Correct (certain): {rha_results['total_correct_certain']}/{rha_results['total_certain']}")
        print(f"  Accuracy (certain): {rha_results['accuracy_certain']:.1%}")
        print(f"  Uncertain cases: {rha_results['total_uncertain']}\n")

        # Compare to paper results
        print(f"Expected results from paper (Cleveland dataset):")
        print(f"  SHA: 83.2% accuracy")
        print(f"  RHA: 85.5% accuracy (excluding uncertain cases)")
        print(f"  RHA: 11.2% uncertain cases\n")

```
