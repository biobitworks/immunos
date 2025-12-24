---
source: /Users/byron/projects/immunos81/utils/interactive_visualizer.py
relative: immunos81/utils/interactive_visualizer.py
generated_at: 2025-12-23 10:28
---

```python
"""
Interactive Visualization Tools for Immunos-81 using Plotly

Provides interactive HTML-based visualizations with hover effects,
zooming, and drill-down capabilities.
"""

import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import numpy as np
from typing import List, Dict, Optional
from collections import defaultdict
from ..core.antigen import Antigen
from ..core.tcell import TCell
from ..engine.recognizer import RecognitionResult


class InteractiveVisualizer:
    """
    Interactive visualization toolkit using Plotly.

    Creates HTML-based interactive plots with:
    - Hover tooltips
    - Zooming and panning
    - Data export capabilities
    - Responsive design
    """

    def __init__(self):
        """Initialize interactive visualizer"""
        self.template = "plotly_white"

    def plot_interactive_confusion_matrix(self,
                                         test_antigens: List[Antigen],
                                         results: List[RecognitionResult],
                                         save_path: Optional[str] = None) -> go.Figure:
        """
        Create interactive confusion matrix with hover details.

        Args:
            test_antigens: List of test antigens with true labels
            results: List of recognition results
            save_path: Optional path to save HTML

        Returns:
            Plotly Figure object
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

        # Create heatmap
        fig = go.Figure(data=go.Heatmap(
            z=cm,
            x=classes,
            y=classes,
            colorscale='Blues',
            text=cm,
            texttemplate='%{text}',
            textfont={"size": 16},
            hovertemplate='True: %{y}<br>Predicted: %{x}<br>Count: %{z}<extra></extra>',
            colorbar=dict(title="Count")
        ))

        # Calculate accuracy
        accuracy = np.trace(cm) / np.sum(cm) if np.sum(cm) > 0 else 0

        fig.update_layout(
            title=dict(
                text=f'Interactive Confusion Matrix<br><sub>Overall Accuracy: {accuracy:.1%}</sub>',
                x=0.5,
                xanchor='center'
            ),
            xaxis_title='Predicted Class',
            yaxis_title='True Class',
            template=self.template,
            height=600,
            width=700
        )

        if save_path:
            fig.write_html(save_path)

        return fig

    def plot_interactive_confidence(self,
                                    results: List[RecognitionResult],
                                    test_antigens: List[Antigen],
                                    save_path: Optional[str] = None) -> go.Figure:
        """
        Create interactive confidence distribution plot.

        Args:
            results: List of recognition results
            test_antigens: List of test antigens with true labels
            save_path: Optional path to save HTML

        Returns:
            Plotly Figure object
        """
        # Categorize results
        correct_data = []
        incorrect_data = []
        certain_data = []
        uncertain_data = []

        for result, antigen in zip(results, test_antigens):
            if result.predicted_class is None:
                continue

            confidence = result.confidence
            is_correct = result.predicted_class == antigen.class_label
            is_certain = result.is_certain

            record = {
                'confidence': confidence,
                'correct': 'Correct' if is_correct else 'Incorrect',
                'certainty': 'Certain' if is_certain else 'Uncertain',
                'true_class': antigen.class_label,
                'pred_class': result.predicted_class
            }

            if is_correct:
                correct_data.append(record)
            else:
                incorrect_data.append(record)

            if is_certain:
                certain_data.append(record)
            else:
                uncertain_data.append(record)

        # Create subplots
        fig = make_subplots(
            rows=1, cols=2,
            subplot_titles=('Correct vs Incorrect', 'Certain vs Uncertain'),
            horizontal_spacing=0.15
        )

        # Plot 1: Correct vs Incorrect
        for data, color, name in [
            ([d['confidence'] for d in correct_data], 'green', 'Correct'),
            ([d['confidence'] for d in incorrect_data], 'red', 'Incorrect')
        ]:
            fig.add_trace(
                go.Histogram(
                    x=data,
                    name=name,
                    opacity=0.7,
                    marker_color=color,
                    hovertemplate='Confidence: %{x}<br>Count: %{y}<extra></extra>'
                ),
                row=1, col=1
            )

        # Plot 2: Certain vs Uncertain
        for data, color, name in [
            ([d['confidence'] for d in certain_data], 'blue', 'Certain'),
            ([d['confidence'] for d in uncertain_data], 'orange', 'Uncertain')
        ]:
            fig.add_trace(
                go.Histogram(
                    x=data,
                    name=name,
                    opacity=0.7,
                    marker_color=color,
                    hovertemplate='Confidence: %{x}<br>Count: %{y}<extra></extra>'
                ),
                row=1, col=2
            )

        fig.update_xaxes(title_text="Confidence Score", row=1, col=1)
        fig.update_xaxes(title_text="Confidence Score", row=1, col=2)
        fig.update_yaxes(title_text="Frequency", row=1, col=1)
        fig.update_yaxes(title_text="Frequency", row=1, col=2)

        fig.update_layout(
            title_text='Interactive Decision Confidence Analysis',
            template=self.template,
            height=500,
            width=1200,
            showlegend=True
        )

        if save_path:
            fig.write_html(save_path)

        return fig

    def plot_interactive_cv_results(self,
                                    cv_results: Dict,
                                    save_path: Optional[str] = None) -> go.Figure:
        """
        Create interactive cross-validation results dashboard.

        Args:
            cv_results: Results from CrossValidator.cross_validate()
            save_path: Optional path to save HTML

        Returns:
            Plotly Figure object
        """
        sha_results = cv_results['strategies'].get('sha', {})
        rha_results = cv_results['strategies'].get('rha', {})

        # Create subplots
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=(
                'SHA: Fold-by-Fold Accuracy',
                'RHA: All vs Certain Accuracy',
                'Strategy Comparison',
                'RHA: Prediction Breakdown'
            ),
            specs=[[{'type': 'scatter'}, {'type': 'scatter'}],
                   [{'type': 'bar'}, {'type': 'pie'}]],
            vertical_spacing=0.12,
            horizontal_spacing=0.12
        )

        # Plot 1: SHA fold-by-fold
        if 'fold_results' in sha_results:
            fold_accuracies = [f['accuracy_overall'] * 100
                             for f in sha_results['fold_results']]
            folds = list(range(1, len(fold_accuracies) + 1))

            fig.add_trace(
                go.Scatter(
                    x=folds,
                    y=fold_accuracies,
                    mode='lines+markers',
                    name='SHA',
                    marker=dict(size=10),
                    line=dict(width=3),
                    hovertemplate='Fold: %{x}<br>Accuracy: %{y:.1f}%<extra></extra>'
                ),
                row=1, col=1
            )

            mean_sha = sha_results.get('overall_accuracy', 0) * 100
            fig.add_hline(
                y=mean_sha,
                line_dash="dash",
                line_color="red",
                annotation_text=f"Mean: {mean_sha:.1f}%",
                row=1, col=1
            )

        # Plot 2: RHA all vs certain
        if 'fold_results' in rha_results:
            fold_accuracies_all = [f['accuracy_overall'] * 100
                                  for f in rha_results['fold_results']]
            fold_accuracies_certain = [f.get('accuracy_certain', 0) * 100
                                      for f in rha_results['fold_results']]
            folds = list(range(1, len(fold_accuracies_all) + 1))

            fig.add_trace(
                go.Scatter(
                    x=folds,
                    y=fold_accuracies_all,
                    mode='lines+markers',
                    name='RHA (All)',
                    marker=dict(size=10),
                    hovertemplate='Fold: %{x}<br>Accuracy: %{y:.1f}%<extra></extra>'
                ),
                row=1, col=2
            )

            fig.add_trace(
                go.Scatter(
                    x=folds,
                    y=fold_accuracies_certain,
                    mode='lines+markers',
                    name='RHA (Certain)',
                    marker=dict(size=10, symbol='square'),
                    hovertemplate='Fold: %{x}<br>Accuracy: %{y:.1f}%<extra></extra>'
                ),
                row=1, col=2
            )

        # Plot 3: Strategy comparison
        strategies = []
        accuracies = []
        colors = []

        if sha_results:
            strategies.append('SHA<br>(All)')
            accuracies.append(sha_results.get('overall_accuracy', 0) * 100)
            colors.append('rgb(99, 110, 250)')

        if rha_results:
            strategies.extend(['RHA<br>(All)', 'RHA<br>(Certain)'])
            accuracies.extend([
                rha_results.get('overall_accuracy', 0) * 100,
                rha_results.get('accuracy_certain', 0) * 100
            ])
            colors.extend(['rgb(239, 85, 59)', 'rgb(0, 204, 150)'])

        fig.add_trace(
            go.Bar(
                x=strategies,
                y=accuracies,
                marker_color=colors,
                text=[f'{acc:.1f}%' for acc in accuracies],
                textposition='outside',
                hovertemplate='%{x}<br>Accuracy: %{y:.1f}%<extra></extra>'
            ),
            row=2, col=1
        )

        # Plot 4: RHA breakdown
        if rha_results and 'total_certain' in rha_results:
            labels = ['Certain<br>Correct', 'Certain<br>Incorrect', 'Uncertain']
            values = [
                rha_results.get('total_correct_certain', 0),
                rha_results.get('total_certain', 0) - rha_results.get('total_correct_certain', 0),
                rha_results.get('total_uncertain', 0)
            ]

            fig.add_trace(
                go.Pie(
                    labels=labels,
                    values=values,
                    marker=dict(colors=['darkgreen', 'darkred', 'orange']),
                    hovertemplate='%{label}<br>Count: %{value}<br>%{percent}<extra></extra>'
                ),
                row=2, col=2
            )

        # Update axes
        fig.update_xaxes(title_text="Fold", row=1, col=1)
        fig.update_xaxes(title_text="Fold", row=1, col=2)
        fig.update_xaxes(title_text="Strategy", row=2, col=1)

        fig.update_yaxes(title_text="Accuracy (%)", row=1, col=1)
        fig.update_yaxes(title_text="Accuracy (%)", row=1, col=2)
        fig.update_yaxes(title_text="Accuracy (%)", row=2, col=1)

        fig.update_layout(
            title_text=f'Interactive Cross-Validation Dashboard ({cv_results["n_folds"]}-Fold)',
            template=self.template,
            height=900,
            width=1400,
            showlegend=True
        )

        if save_path:
            fig.write_html(save_path)

        return fig

    def plot_interactive_variable_importance(self,
                                            tcells: Dict[str, TCell],
                                            test_antigens: List[Antigen],
                                            top_n: int = 10,
                                            save_path: Optional[str] = None) -> go.Figure:
        """
        Create interactive variable importance plot.

        Args:
            tcells: Dictionary of trained T cells
            test_antigens: List of test antigens
            top_n: Number of top variables to show
            save_path: Optional path to save HTML

        Returns:
            Plotly Figure object
        """
        # Collect affinity contributions
        variable_contributions = defaultdict(list)

        for tcell in tcells.values():
            for antigen in test_antigens[:100]:  # Sample
                if tcell.matches_structure(antigen):
                    for clone in tcell.clones:
                        affinity_vector = clone.calculate_affinity_vector(antigen)
                        for var_num, affinity in affinity_vector.items():
                            variable_contributions[var_num].append(affinity)

        # Calculate statistics
        variable_stats = {}
        for var_num, affinities in variable_contributions.items():
            if affinities:
                variable_stats[var_num] = {
                    'mean': np.mean(affinities),
                    'std': np.std(affinities),
                    'min': np.min(affinities),
                    'max': np.max(affinities),
                    'count': len(affinities)
                }

        # Get library and sort
        library = test_antigens[0].library
        sorted_vars = sorted(variable_stats.items(),
                           key=lambda x: x[1]['mean'],
                           reverse=True)[:top_n]

        var_names = [library.get_name(var_num) for var_num, _ in sorted_vars]
        means = [stats['mean'] for _, stats in sorted_vars]
        stds = [stats['std'] for _, stats in sorted_vars]

        # Create horizontal bar chart with error bars
        fig = go.Figure()

        fig.add_trace(go.Bar(
            y=var_names,
            x=means,
            orientation='h',
            error_x=dict(type='data', array=stds),
            marker_color='rgb(55, 83, 109)',
            hovertemplate='%{y}<br>Mean Affinity: %{x:.3f}<br>Std Dev: %{error_x.array:.3f}<extra></extra>'
        ))

        fig.update_layout(
            title_text=f'Interactive Variable Importance (Top {top_n})',
            xaxis_title='Average Affinity Contribution',
            yaxis_title='Variable',
            template=self.template,
            height=max(400, top_n * 40),
            width=900
        )

        if save_path:
            fig.write_html(save_path)

        return fig

```
