---
source: /Users/byron/projects/immunos81/dashboard/app.py
relative: immunos81/dashboard/app.py
generated_at: 2025-12-23 10:28
---

```python
"""
Immunos-81 Interactive Web Dashboard

A Dash web application providing real-time visualization and analysis
of Immunos-81 classification results.
"""

import dash
from dash import dcc, html, Input, Output, State
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np

# Import Immunos-81 components
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from immunos81.data.cleveland_loader import ClevelandDataLoader
from immunos81.engine.parallel_trainer import ParallelTrainer
from immunos81.engine.parallel_recognizer import ParallelRecognizer, RecognitionStrategy
from immunos81.utils.interactive_visualizer import InteractiveVisualizer
from immunos81.utils.benchmark import Benchmark


# Initialize Dash app
app = dash.Dash(__name__, suppress_callback_exceptions=True)
app.title = "Immunos-81 Dashboard"

# Global storage for trained model and data
trained_model = None
test_results = None
dataset = None


# App layout
app.layout = html.Div([
    # Header
    html.Div([
        html.H1("Immunos-81: Interactive Analysis Dashboard",
                style={'textAlign': 'center', 'color': '#2c3e50', 'marginBottom': 30}),
        html.H3("Artificial Immune System for Pattern Recognition",
                style={'textAlign': 'center', 'color': '#7f8c8d', 'marginBottom': 50})
    ]),

    # Control Panel
    html.Div([
        html.H3("Control Panel", style={'color': '#34495e'}),
        html.Hr(),

        html.Div([
            html.Label("Dataset:", style={'fontWeight': 'bold'}),
            dcc.Dropdown(
                id='dataset-selector',
                options=[
                    {'label': 'Cleveland Heart Disease (303 patients)', 'value': 'cleveland'}
                ],
                value='cleveland',
                style={'width': '100%'}
            )
        ], style={'marginBottom': 20}),

        html.Div([
            html.Label("Recognition Strategy:", style={'fontWeight': 'bold'}),
            dcc.RadioItems(
                id='strategy-selector',
                options=[
                    {'label': ' SHA (Simple Highest Avidity)', 'value': 'sha'},
                    {'label': ' RHA (Relative Highest Avidity, 5% threshold)', 'value': 'rha'}
                ],
                value='sha',
                style={'marginTop': 10}
            )
        ], style={'marginBottom': 20}),

        html.Div([
            html.Label("Number of Workers:", style={'fontWeight': 'bold'}),
            dcc.Slider(
                id='workers-slider',
                min=1,
                max=8,
                step=1,
                value=4,
                marks={i: str(i) for i in range(1, 9)},
                tooltip={"placement": "bottom", "always_visible": True}
            )
        ], style={'marginBottom': 30}),

        html.Button('Train Model', id='train-button', n_clicks=0,
                   style={'width': '100%', 'padding': 15, 'fontSize': 16,
                          'backgroundColor': '#3498db', 'color': 'white',
                          'border': 'none', 'borderRadius': 5, 'cursor': 'pointer'}),

        html.Div(id='training-status', style={'marginTop': 20, 'textAlign': 'center',
                                              'fontSize': 14, 'color': '#27ae60'})

    ], style={'padding': 30, 'backgroundColor': '#ecf0f1', 'borderRadius': 10,
              'marginBottom': 30}),

    # Results Section
    html.Div([
        html.H3("Classification Results", style={'color': '#34495e'}),
        html.Hr(),

        # Performance Metrics
        html.Div(id='metrics-display', style={'marginBottom': 30}),

        # Tabs for different visualizations
        dcc.Tabs(id='viz-tabs', value='confusion', children=[
            dcc.Tab(label='Confusion Matrix', value='confusion'),
            dcc.Tab(label='Decision Confidence', value='confidence'),
            dcc.Tab(label='Variable Importance', value='importance'),
            dcc.Tab(label='Performance Benchmark', value='benchmark')
        ]),

        html.Div(id='viz-content', style={'marginTop': 30})

    ], style={'padding': 30, 'backgroundColor': 'white', 'borderRadius': 10,
              'boxShadow': '0 2px 4px rgba(0,0,0,0.1)'})

], style={'maxWidth': '1400px', 'margin': '0 auto', 'padding': '40px 20px',
          'fontFamily': 'Arial, sans-serif', 'backgroundColor': '#f5f5f5'})


@app.callback(
    Output('training-status', 'children'),
    Output('metrics-display', 'children'),
    Output('viz-content', 'children'),
    Input('train-button', 'n_clicks'),
    Input('viz-tabs', 'value'),
    State('dataset-selector', 'value'),
    State('strategy-selector', 'value'),
    State('workers-slider', 'value'),
    prevent_initial_call=True
)
def update_dashboard(n_clicks, active_tab, dataset_name, strategy, n_workers):
    """Main callback to handle training and visualization updates"""
    global trained_model, test_results, dataset

    ctx = dash.callback_context
    if not ctx.triggered:
        return "", "", ""

    trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]

    # Training triggered
    if trigger_id == 'train-button' and n_clicks > 0:
        # Load dataset
        loader = ClevelandDataLoader()
        antigens, library = loader.load_data()
        dataset = antigens

        # Split data
        n_train = int(len(antigens) * 0.7)
        train_data = antigens[:n_train]
        test_data = antigens[n_train:]

        # Train model
        trainer = ParallelTrainer(library, n_workers=n_workers)
        tcells = trainer.train(train_data)

        # Test model
        strategy_enum = RecognitionStrategy.SHA if strategy == 'sha' else RecognitionStrategy.RHA
        recognizer = ParallelRecognizer(tcells, strategy=strategy_enum, n_workers=n_workers)
        results = recognizer.batch_recognize(test_data)
        metrics = recognizer.evaluate(test_data)

        # Store globally
        trained_model = {'tcells': tcells, 'test_data': test_data,
                        'results': results, 'metrics': metrics,
                        'strategy': strategy_enum}

        status = f"✓ Model trained successfully! Using {n_workers} workers."

        # Create metrics display
        metrics_div = create_metrics_display(metrics, strategy)

        # Create initial visualization
        viz_div = create_visualization(active_tab, trained_model)

        return status, metrics_div, viz_div

    # Tab changed
    elif trigger_id == 'viz-tabs' and trained_model is not None:
        viz_div = create_visualization(active_tab, trained_model)
        return dash.no_update, dash.no_update, viz_div

    return "", "", ""


def create_metrics_display(metrics, strategy):
    """Create metrics display cards"""
    cards = []

    # Overall accuracy
    cards.append(html.Div([
        html.H4(f"{metrics['accuracy_overall']:.1%}", style={'fontSize': 36, 'margin': 0, 'color': '#27ae60'}),
        html.P("Overall Accuracy", style={'margin': 5, 'color': '#7f8c8d'})
    ], style={'flex': 1, 'padding': 20, 'backgroundColor': '#ecf0f1',
              'borderRadius': 8, 'textAlign': 'center'}))

    # Correct predictions
    cards.append(html.Div([
        html.H4(f"{metrics['correct']}/{metrics['total']}", style={'fontSize': 36, 'margin': 0, 'color': '#3498db'}),
        html.P("Correct Predictions", style={'margin': 5, 'color': '#7f8c8d'})
    ], style={'flex': 1, 'padding': 20, 'backgroundColor': '#ecf0f1',
              'borderRadius': 8, 'textAlign': 'center', 'marginLeft': 15}))

    # RHA-specific metrics
    if strategy == 'rha' and 'accuracy_certain' in metrics:
        cards.append(html.Div([
            html.H4(f"{metrics['accuracy_certain']:.1%}", style={'fontSize': 36, 'margin': 0, 'color': '#e74c3c'}),
            html.P("Certain Accuracy", style={'margin': 5, 'color': '#7f8c8d'})
        ], style={'flex': 1, 'padding': 20, 'backgroundColor': '#ecf0f1',
                  'borderRadius': 8, 'textAlign': 'center', 'marginLeft': 15}))

        cards.append(html.Div([
            html.H4(f"{metrics['uncertain_attempts']}", style={'fontSize': 36, 'margin': 0, 'color': '#f39c12'}),
            html.P("Uncertain Cases", style={'margin': 5, 'color': '#7f8c8d'})
        ], style={'flex': 1, 'padding': 20, 'backgroundColor': '#ecf0f1',
                  'borderRadius': 8, 'textAlign': 'center', 'marginLeft': 15}))

    return html.Div(cards, style={'display': 'flex', 'marginBottom': 30})


def create_visualization(tab, model_data):
    """Create visualization based on selected tab"""
    if model_data is None:
        return html.Div("Please train a model first.", style={'textAlign': 'center', 'padding': 50, 'color': '#7f8c8d'})

    viz = InteractiveVisualizer()

    if tab == 'confusion':
        fig = viz.plot_interactive_confusion_matrix(
            model_data['test_data'],
            model_data['results']
        )
        return dcc.Graph(figure=fig, style={'height': '600px'})

    elif tab == 'confidence':
        fig = viz.plot_interactive_confidence(
            model_data['results'],
            model_data['test_data']
        )
        return dcc.Graph(figure=fig, style={'height': '500px'})

    elif tab == 'importance':
        fig = viz.plot_interactive_variable_importance(
            model_data['tcells'],
            model_data['test_data'],
            top_n=10
        )
        return dcc.Graph(figure=fig, style={'height': '600px'})

    elif tab == 'benchmark':
        return create_benchmark_view()

    return html.Div("Visualization not available")


def create_benchmark_view():
    """Create benchmark comparison view"""
    return html.Div([
        html.H4("Performance Benchmark", style={'color': '#34495e'}),
        html.P("Comparing sequential vs parallel implementations...",
               style={'color': '#7f8c8d', 'marginBottom': 20}),
        html.Button('Run Benchmark', id='run-benchmark-btn',
                   style={'padding': '10px 30px', 'fontSize': 14,
                          'backgroundColor': '#e74c3c', 'color': 'white',
                          'border': 'none', 'borderRadius': 5, 'cursor': 'pointer'}),
        html.Div(id='benchmark-results', style={'marginTop': 30})
    ])


if __name__ == '__main__':
    print("Starting Immunos-81 Dashboard...")
    print("Open your browser to: http://127.0.0.1:8050")
    app.run_server(debug=True, host='127.0.0.1', port=8050)

```
