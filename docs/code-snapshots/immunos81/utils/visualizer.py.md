---
source: /Users/byron/projects/immunos81/utils/visualizer.py
relative: immunos81/utils/visualizer.py
generated_at: 2025-12-23 10:28
---

```python
"""
Visualization Tools for Immunos-81

Provides comprehensive visualization of Immunos-81 classification results,
including avidity distributions, confusion matrices, and performance metrics.
"""

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from typing import List, Dict, Tuple, Optional
from collections import defaultdict
from ..core.antigen import Antigen
from ..core.tcell import TCell
from ..core.clone import Clone
from ..engine.recognizer import Recognizer, RecognitionResult, RecognitionStrategy


class Immunos81Visualizer:
    """
    Visualization toolkit for Immunos-81 results.

    Provides methods to visualize:
    - Avidity score distributions
    - Confusion matrices
    - Decision confidence plots
    - Cross-validation results
    - Variable importance
    """

    def __init__(self, style='seaborn-v0_8-darkgrid'):
        """
        Initialize visualizer with plotting style.

        Args:
            style: Matplotlib style to use
        """
        try:
            plt.style.use(style)
        except:
            # Fallback if style not available
            sns.set_theme()

        self.colors = sns.color_palette("husl", 10)

    def plot_avidity_distribution(self, tcells: Dict[str, TCell],
                                  test_antigens: List[Antigen],
                                  save_path: Optional[str] = None):
        """
        Plot avidity score distributions for each class.

        Shows how well clones from each class recognize test antigens.

        Args:
            tcells: Dictionary of trained T cells
            test_antigens: List of test antigens
            save_path: Optional path to save figure
        """
        fig, axes = plt.subplots(1, len(tcells), figsize=(6*len(tcells), 5))
        if len(tcells) == 1:
            axes = [axes]

        for idx, (class_label, tcell) in enumerate(tcells.items()):
            avidities = []

            # Calculate avidity for each test antigen
            for antigen in test_antigens:
                if tcell.matches_structure(antigen):
                    for clone in tcell.clones:
                        avidity = clone.calculate_avidity(antigen)
                        avidities.append(avidity)

            # Plot distribution
            axes[idx].hist(avidities, bins=30, alpha=0.7,
                          color=self.colors[idx], edgecolor='black')
            axes[idx].set_title(f'Avidity Distribution: {class_label}',
                               fontsize=12, fontweight='bold')
            axes[idx].set_xlabel('Avidity Score')
            axes[idx].set_ylabel('Frequency')
            axes[idx].grid(True, alpha=0.3)

            # Add statistics
            if avidities:
                mean_avidity = np.mean(avidities)
                median_avidity = np.median(avidities)
                axes[idx].axvline(mean_avidity, color='red', linestyle='--',
                                 label=f'Mean: {mean_avidity:.2f}')
                axes[idx].axvline(median_avidity, color='blue', linestyle='--',
                                 label=f'Median: {median_avidity:.2f}')
                axes[idx].legend()

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()

    def plot_confusion_matrix(self, test_antigens: List[Antigen],
                             results: List[RecognitionResult],
                             save_path: Optional[str] = None):
        """
        Plot confusion matrix for classification results.

        Args:
            test_antigens: List of test antigens with true labels
            results: List of recognition results
            save_path: Optional path to save figure
        """
        # Get unique classes
        true_labels = [ag.class_label for ag in test_antigens]
        pred_labels = [r.predicted_class for r in results]
        classes = sorted(set(true_labels))

        # Build confusion matrix
        n_classes = len(classes)
        cm = np.zeros((n_classes, n_classes), dtype=int)

        for true_label, pred_label in zip(true_labels, pred_labels):
            if pred_label is not None:
                true_idx = classes.index(true_label)
                pred_idx = classes.index(pred_label)
                cm[true_idx, pred_idx] += 1

        # Plot
        fig, ax = plt.subplots(figsize=(8, 7))
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues',
                   xticklabels=classes, yticklabels=classes,
                   cbar_kws={'label': 'Count'}, ax=ax)

        ax.set_title('Confusion Matrix', fontsize=14, fontweight='bold', pad=20)
        ax.set_xlabel('Predicted Class', fontsize=12)
        ax.set_ylabel('True Class', fontsize=12)

        # Add accuracy annotation
        accuracy = np.trace(cm) / np.sum(cm) if np.sum(cm) > 0 else 0
        plt.text(0.5, -0.15, f'Overall Accuracy: {accuracy:.1%}',
                ha='center', transform=ax.transAxes, fontsize=12,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()

    def plot_decision_confidence(self, results: List[RecognitionResult],
                                test_antigens: List[Antigen],
                                save_path: Optional[str] = None):
        """
        Plot decision confidence distribution.

        Shows the distribution of confidence scores for correct vs incorrect
        predictions, and certain vs uncertain predictions (for RHA).

        Args:
            results: List of recognition results
            test_antigens: List of test antigens with true labels
            save_path: Optional path to save figure
        """
        # Categorize results
        correct_confidences = []
        incorrect_confidences = []
        certain_correct = []
        certain_incorrect = []
        uncertain_confidences = []

        for result, antigen in zip(results, test_antigens):
            if result.predicted_class is None:
                continue

            confidence = result.confidence
            is_correct = result.predicted_class == antigen.class_label
            is_certain = result.is_certain

            if is_correct:
                correct_confidences.append(confidence)
                if is_certain:
                    certain_correct.append(confidence)
            else:
                incorrect_confidences.append(confidence)
                if is_certain:
                    certain_incorrect.append(confidence)

            if not is_certain:
                uncertain_confidences.append(confidence)

        # Create subplots
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # Plot 1: Correct vs Incorrect
        axes[0].hist([correct_confidences, incorrect_confidences],
                    bins=20, label=['Correct', 'Incorrect'],
                    color=['green', 'red'], alpha=0.6, edgecolor='black')
        axes[0].set_title('Confidence Distribution: Correct vs Incorrect',
                         fontsize=12, fontweight='bold')
        axes[0].set_xlabel('Confidence Score')
        axes[0].set_ylabel('Frequency')
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)

        # Plot 2: Certain vs Uncertain (for RHA)
        if uncertain_confidences:
            axes[1].hist([certain_correct, certain_incorrect, uncertain_confidences],
                        bins=20, label=['Certain Correct', 'Certain Incorrect', 'Uncertain'],
                        color=['darkgreen', 'darkred', 'orange'], alpha=0.6,
                        edgecolor='black')
            axes[1].set_title('RHA Confidence Distribution',
                             fontsize=12, fontweight='bold')
            axes[1].set_xlabel('Confidence Score')
            axes[1].set_ylabel('Frequency')
            axes[1].legend()
            axes[1].grid(True, alpha=0.3)
        else:
            axes[1].text(0.5, 0.5, 'No uncertain predictions\n(SHA strategy used)',
                        ha='center', va='center', transform=axes[1].transAxes,
                        fontsize=12)
            axes[1].set_title('RHA Confidence Distribution (N/A)',
                             fontsize=12, fontweight='bold')

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()

    def plot_cross_validation_results(self, cv_results: Dict,
                                      save_path: Optional[str] = None):
        """
        Plot cross-validation results comparing SHA and RHA strategies.

        Args:
            cv_results: Results from CrossValidator.cross_validate()
            save_path: Optional path to save figure
        """
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        # Extract results for each strategy
        sha_results = cv_results['strategies'].get('sha', {})
        rha_results = cv_results['strategies'].get('rha', {})

        # Plot 1: Fold-by-fold accuracy (SHA)
        if 'fold_results' in sha_results:
            fold_accuracies = [f['accuracy_overall'] * 100
                             for f in sha_results['fold_results']]
            folds = list(range(1, len(fold_accuracies) + 1))

            axes[0, 0].plot(folds, fold_accuracies, marker='o',
                           color=self.colors[0], linewidth=2, markersize=8)
            axes[0, 0].axhline(sha_results.get('overall_accuracy', 0) * 100,
                              color='red', linestyle='--', label='Mean')
            axes[0, 0].set_title('SHA: Fold-by-Fold Accuracy',
                                fontsize=12, fontweight='bold')
            axes[0, 0].set_xlabel('Fold')
            axes[0, 0].set_ylabel('Accuracy (%)')
            axes[0, 0].grid(True, alpha=0.3)
            axes[0, 0].legend()
            axes[0, 0].set_ylim([0, 100])

        # Plot 2: Fold-by-fold accuracy (RHA)
        if 'fold_results' in rha_results:
            fold_accuracies_all = [f['accuracy_overall'] * 100
                                  for f in rha_results['fold_results']]
            fold_accuracies_certain = [f.get('accuracy_certain', 0) * 100
                                      for f in rha_results['fold_results']]
            folds = list(range(1, len(fold_accuracies_all) + 1))

            axes[0, 1].plot(folds, fold_accuracies_all, marker='o',
                           label='All Predictions', color=self.colors[1],
                           linewidth=2, markersize=8)
            axes[0, 1].plot(folds, fold_accuracies_certain, marker='s',
                           label='Certain Only', color=self.colors[2],
                           linewidth=2, markersize=8)
            axes[0, 1].axhline(rha_results.get('overall_accuracy', 0) * 100,
                              color='red', linestyle='--', label='Mean (All)')
            axes[0, 1].set_title('RHA: Fold-by-Fold Accuracy',
                                fontsize=12, fontweight='bold')
            axes[0, 1].set_xlabel('Fold')
            axes[0, 1].set_ylabel('Accuracy (%)')
            axes[0, 1].grid(True, alpha=0.3)
            axes[0, 1].legend()
            axes[0, 1].set_ylim([0, 100])

        # Plot 3: Strategy comparison
        strategies = []
        accuracies = []
        if sha_results:
            strategies.append('SHA\n(All)')
            accuracies.append(sha_results.get('overall_accuracy', 0) * 100)
        if rha_results:
            strategies.extend(['RHA\n(All)', 'RHA\n(Certain)'])
            accuracies.extend([
                rha_results.get('overall_accuracy', 0) * 100,
                rha_results.get('accuracy_certain', 0) * 100
            ])

        bars = axes[1, 0].bar(strategies, accuracies,
                             color=self.colors[:len(strategies)],
                             edgecolor='black', linewidth=1.5)
        axes[1, 0].set_title('Strategy Comparison',
                            fontsize=12, fontweight='bold')
        axes[1, 0].set_ylabel('Accuracy (%)')
        axes[1, 0].set_ylim([0, 100])
        axes[1, 0].grid(True, alpha=0.3, axis='y')

        # Add value labels on bars
        for bar, acc in zip(bars, accuracies):
            height = bar.get_height()
            axes[1, 0].text(bar.get_x() + bar.get_width()/2., height,
                           f'{acc:.1f}%', ha='center', va='bottom',
                           fontweight='bold')

        # Plot 4: Prediction breakdown
        if rha_results and 'total_certain' in rha_results:
            labels = ['Certain\nCorrect', 'Certain\nIncorrect', 'Uncertain']
            sizes = [
                rha_results.get('total_correct_certain', 0),
                rha_results.get('total_certain', 0) - rha_results.get('total_correct_certain', 0),
                rha_results.get('total_uncertain', 0)
            ]
            colors_pie = ['darkgreen', 'darkred', 'orange']

            axes[1, 1].pie(sizes, labels=labels, autopct='%1.1f%%',
                          colors=colors_pie, startangle=90,
                          textprops={'fontsize': 10, 'fontweight': 'bold'})
            axes[1, 1].set_title('RHA: Prediction Breakdown',
                                fontsize=12, fontweight='bold')
        else:
            axes[1, 1].text(0.5, 0.5, 'RHA results not available',
                           ha='center', va='center',
                           transform=axes[1, 1].transAxes, fontsize=12)

        plt.suptitle(f'Immunos-81 Cross-Validation Results ({cv_results["n_folds"]}-Fold)',
                    fontsize=14, fontweight='bold', y=1.00)
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()

    def plot_variable_importance(self, tcells: Dict[str, TCell],
                                test_antigens: List[Antigen],
                                top_n: int = 10,
                                save_path: Optional[str] = None):
        """
        Analyze and plot variable importance based on affinity contributions.

        Args:
            tcells: Dictionary of trained T cells
            test_antigens: List of test antigens
            top_n: Number of top variables to show
            save_path: Optional path to save figure
        """
        # Collect affinity contributions for each variable
        variable_contributions = defaultdict(list)

        for tcell in tcells.values():
            for antigen in test_antigens[:100]:  # Sample for efficiency
                if tcell.matches_structure(antigen):
                    for clone in tcell.clones:
                        affinity_vector = clone.calculate_affinity_vector(antigen)
                        for var_num, affinity in affinity_vector.items():
                            variable_contributions[var_num].append(affinity)

        # Calculate mean contribution for each variable
        variable_importance = {}
        for var_num, affinities in variable_contributions.items():
            if affinities:
                variable_importance[var_num] = np.mean(affinities)

        # Get library from first antigen
        library = test_antigens[0].library

        # Sort by importance and get top N
        sorted_vars = sorted(variable_importance.items(),
                           key=lambda x: x[1], reverse=True)[:top_n]

        # Plot
        fig, ax = plt.subplots(figsize=(10, 6))

        var_names = [library.get_name(var_num) for var_num, _ in sorted_vars]
        importances = [imp for _, imp in sorted_vars]

        bars = ax.barh(var_names, importances, color=self.colors[0],
                      edgecolor='black', linewidth=1.5)
        ax.set_xlabel('Average Affinity Contribution', fontsize=12)
        ax.set_title(f'Top {top_n} Most Important Variables',
                    fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='x')

        # Add value labels
        for bar, imp in zip(bars, importances):
            width = bar.get_width()
            ax.text(width, bar.get_y() + bar.get_height()/2.,
                   f'{imp:.3f}', ha='left', va='center',
                   fontweight='bold', fontsize=9)

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()

```
