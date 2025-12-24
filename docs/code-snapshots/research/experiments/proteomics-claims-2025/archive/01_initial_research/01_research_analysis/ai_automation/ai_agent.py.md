---
source: /Users/byron/projects/research/experiments/proteomics-claims-2025/archive/01_initial_research/01_research_analysis/ai_automation/ai_agent.py
relative: research/experiments/proteomics-claims-2025/archive/01_initial_research/01_research_analysis/ai_automation/ai_agent.py
generated_at: 2025-12-23 10:28
---

```python
"""
AI Agent for Automated Bioinformatics Analysis
Handles proteomic data analysis for finding group evaluations
"""

import logging
from typing import Dict, List, Tuple, Optional, Any
import numpy as np
import pandas as pd
from dataclasses import dataclass
from enum import Enum

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class EvaluationResult(Enum):
    SUPPORTED = "SUPPORTED"
    REFUTED = "REFUTED"
    UNSURE = "UNSURE"

@dataclass
class AnalysisResult:
    statement_id: str
    evaluation: EvaluationResult
    evidence: Dict[str, Any]
    explanation: str
    code_used: str
    confidence: float

class BioinformaticsAgent:
    """
    Main AI agent for automated proteomic analysis
    """

    def __init__(self, data_path: str):
        self.data_path = data_path
        self.adata = None
        self.results = {}
        self.metadata = {}
        logger.info(f"Initializing BioinformaticsAgent with data: {data_path}")

    def load_data(self):
        """Load and validate the H5AD dataset"""
        try:
            import scanpy as sc
            self.adata = sc.read_h5ad(self.data_path)
            self._validate_data()
            self._extract_metadata()
            logger.info(f"Successfully loaded data: {self.adata.shape}")
            return True
        except Exception as e:
            logger.error(f"Failed to load data: {str(e)}")
            return False

    def _validate_data(self):
        """Validate required columns and data integrity"""
        required_obs = ['tau_status', 'MC1', 'pseudotime']
        missing = [col for col in required_obs if col not in self.adata.obs.columns]
        if missing:
            raise ValueError(f"Missing required columns: {missing}")

        logger.info("Data validation passed")

    def _extract_metadata(self):
        """Extract and store metadata information"""
        self.metadata = {
            'n_cells': self.adata.n_obs,
            'n_proteins': self.adata.n_vars,
            'tau_positive': sum(self.adata.obs['tau_status'] == 'positive'),
            'tau_negative': sum(self.adata.obs['tau_status'] == 'negative'),
            'mc1_range': (self.adata.obs['MC1'].min(), self.adata.obs['MC1'].max()),
            'pseudotime_range': (self.adata.obs['pseudotime'].min(),
                                self.adata.obs['pseudotime'].max())
        }
        logger.info(f"Metadata extracted: {self.metadata}")

    def analyze_statement(self, statement_type: str, statement_params: Dict) -> AnalysisResult:
        """
        Analyze a specific statement based on its type
        """
        analysis_methods = {
            'differential_expression': self._analyze_de,
            'correlation': self._analyze_correlation,
            'segmented_regression': self._analyze_segmented,
            'biphasic': self._analyze_biphasic,
            'protein_upregulation': self._analyze_upregulation,
            'sliding_window': self._analyze_sliding_window
        }

        if statement_type not in analysis_methods:
            logger.error(f"Unknown statement type: {statement_type}")
            return self._create_unsure_result(statement_params.get('statement_id', 'unknown'))

        try:
            return analysis_methods[statement_type](statement_params)
        except Exception as e:
            logger.error(f"Analysis failed for {statement_type}: {str(e)}")
            return self._create_unsure_result(statement_params.get('statement_id', 'unknown'))

    def _analyze_de(self, params: Dict) -> AnalysisResult:
        """Perform differential expression analysis"""
        from scipy.stats import ttest_ind
        from statsmodels.stats.multitest import multipletests

        tau_pos = self.adata[self.adata.obs['tau_status'] == 'positive']
        tau_neg = self.adata[self.adata.obs['tau_status'] == 'negative']

        p_values = []
        fold_changes = []

        for gene in self.adata.var_names[:100]:  # Sample for demonstration
            pos_expr = tau_pos[:, gene].X.flatten()
            neg_expr = tau_neg[:, gene].X.flatten()

            _, p_val = ttest_ind(pos_expr, neg_expr)
            p_values.append(p_val)

            fc = np.mean(pos_expr) - np.mean(neg_expr)  # log2 scale
            fold_changes.append(fc)

        # FDR correction
        rejected, p_adjusted, _, _ = multipletests(p_values, method='fdr_bh')
        n_significant = sum(rejected)

        evidence = {
            'n_tested': len(p_values),
            'n_significant': n_significant,
            'percentage': (n_significant / len(p_values)) * 100
        }

        return AnalysisResult(
            statement_id=params['statement_id'],
            evaluation=EvaluationResult.SUPPORTED if evidence['percentage'] > 30 else EvaluationResult.REFUTED,
            evidence=evidence,
            explanation=f"Found {n_significant}/{len(p_values)} significantly altered proteins",
            code_used="differential_expression_analysis",
            confidence=0.8
        )

    def _analyze_correlation(self, params: Dict) -> AnalysisResult:
        """Analyze correlation between proteins"""
        from scipy.stats import pearsonr

        protein1 = params.get('protein1')
        protein2 = params.get('protein2')

        if protein1 in self.adata.var_names and protein2 in self.adata.var_names:
            expr1 = self.adata[:, protein1].X.flatten()
            expr2 = self.adata[:, protein2].X.flatten()

            corr, p_val = pearsonr(expr1, expr2)

            evidence = {
                'correlation': corr,
                'p_value': p_val,
                'n_samples': len(expr1)
            }

            expected_corr = params.get('expected_correlation', 0)
            tolerance = params.get('tolerance', 0.1)

            is_supported = abs(corr - expected_corr) < tolerance

            return AnalysisResult(
                statement_id=params['statement_id'],
                evaluation=EvaluationResult.SUPPORTED if is_supported else EvaluationResult.REFUTED,
                evidence=evidence,
                explanation=f"Correlation: r={corr:.3f}, p={p_val:.3e}",
                code_used="correlation_analysis",
                confidence=0.9
            )

        return self._create_unsure_result(params['statement_id'])

    def _analyze_segmented(self, params: Dict) -> AnalysisResult:
        """Perform segmented regression analysis"""
        from scipy.optimize import curve_fit

        def piecewise_linear(x, x0, y0, b1, b2):
            return np.piecewise(x, [x < x0],
                              [lambda x: y0 + b1*(x-x0),
                               lambda x: y0 + b2*(x-x0)])

        x_var = params.get('x_variable')
        y_var = params.get('y_variable')

        if x_var == 'MC1':
            x_data = self.adata.obs['MC1'].values
        else:
            x_data = self.adata.obs[x_var].values

        if y_var in self.adata.var_names:
            y_data = self.adata[:, y_var].X.flatten()
        else:
            # Calculate composite score if needed
            y_data = self._calculate_composite_score(y_var)

        try:
            popt, _ = curve_fit(piecewise_linear, x_data, y_data,
                              p0=[2.8, np.mean(y_data), -0.003, -0.4])

            breakpoint = popt[0]
            slope1, slope2 = popt[2], popt[3]

            evidence = {
                'breakpoint': breakpoint,
                'slope_before': slope1,
                'slope_after': slope2
            }

            expected_bp = params.get('expected_breakpoint', 2.8)
            is_supported = abs(breakpoint - expected_bp) < 0.2

            return AnalysisResult(
                statement_id=params['statement_id'],
                evaluation=EvaluationResult.SUPPORTED if is_supported else EvaluationResult.REFUTED,
                evidence=evidence,
                explanation=f"Breakpoint at {breakpoint:.3f}, slopes: {slope1:.3f}, {slope2:.3f}",
                code_used="segmented_regression",
                confidence=0.75
            )
        except:
            return self._create_unsure_result(params['statement_id'])

    def _analyze_biphasic(self, params: Dict) -> AnalysisResult:
        """Analyze biphasic behavior along pseudotime"""
        system = params.get('system')  # 'vatpase' or 'proteasome'

        score = self._calculate_composite_score(system)
        pseudotime = self.adata.obs['pseudotime'].values

        # Simple breakpoint detection
        from sklearn.tree import DecisionTreeRegressor

        dt = DecisionTreeRegressor(max_depth=1)
        dt.fit(pseudotime.reshape(-1, 1), score)

        breakpoint = dt.tree_.threshold[0]

        evidence = {
            'breakpoint': breakpoint,
            'system': system
        }

        expected_bp = params.get('expected_breakpoint', 0.65)
        is_supported = abs(breakpoint - expected_bp) < 0.1

        return AnalysisResult(
            statement_id=params['statement_id'],
            evaluation=EvaluationResult.SUPPORTED if is_supported else EvaluationResult.REFUTED,
            evidence=evidence,
            explanation=f"{system} breakpoint at pseudotime {breakpoint:.3f}",
            code_used="biphasic_analysis",
            confidence=0.7
        )

    def _analyze_upregulation(self, params: Dict) -> AnalysisResult:
        """Analyze protein upregulation"""
        protein = params.get('protein')

        if protein not in self.adata.var_names:
            return self._create_unsure_result(params['statement_id'])

        tau_pos = self.adata[self.adata.obs['tau_status'] == 'positive']
        tau_neg = self.adata[self.adata.obs['tau_status'] == 'negative']

        pos_expr = tau_pos[:, protein].X.flatten().mean()
        neg_expr = tau_neg[:, protein].X.flatten().mean()

        log2fc = pos_expr - neg_expr  # Already in log2 scale

        evidence = {
            'protein': protein,
            'log2FC': log2fc,
            'tau_pos_mean': pos_expr,
            'tau_neg_mean': neg_expr
        }

        expected_fc = params.get('expected_log2fc', 3.4)
        tolerance = params.get('tolerance', 0.5)

        is_supported = abs(log2fc - expected_fc) < tolerance

        return AnalysisResult(
            statement_id=params['statement_id'],
            evaluation=EvaluationResult.SUPPORTED if is_supported else EvaluationResult.REFUTED,
            evidence=evidence,
            explanation=f"{protein} log2FC = {log2fc:.3f}",
            code_used="upregulation_analysis",
            confidence=0.85
        )

    def _analyze_sliding_window(self, params: Dict) -> AnalysisResult:
        """Perform sliding window correlation analysis"""
        from scipy.stats import pearsonr

        protein1 = params.get('protein1')
        protein2 = params.get('protein2')
        window_size = params.get('window_size', 20)

        if protein1 not in self.adata.var_names or protein2 not in self.adata.var_names:
            return self._create_unsure_result(params['statement_id'])

        # Sort by pseudotime
        sorted_idx = np.argsort(self.adata.obs['pseudotime'])
        sorted_data = self.adata[sorted_idx]

        correlations = []
        positions = []

        for i in range(len(sorted_data) - window_size + 1):
            window = sorted_data[i:i+window_size]
            expr1 = window[:, protein1].X.flatten()
            expr2 = window[:, protein2].X.flatten()

            corr, _ = pearsonr(expr1, expr2)
            correlations.append(corr)
            positions.append(window.obs['pseudotime'].mean())

        positions = np.array(positions)
        correlations = np.array(correlations)

        # Calculate early and late correlations
        early_corr = correlations[positions < 0.33].mean()
        late_corr = correlations[positions > 0.67].mean()

        # Trend analysis
        trend_corr, trend_p = pearsonr(positions, correlations)

        evidence = {
            'early_correlation': early_corr,
            'late_correlation': late_corr,
            'trend_correlation': trend_corr,
            'trend_p_value': trend_p
        }

        # Check if pattern matches expected
        expected_early = params.get('expected_early', -0.4)
        expected_late = params.get('expected_late', 0.5)

        early_match = abs(early_corr - expected_early) < 0.15
        late_match = abs(late_corr - expected_late) < 0.15

        is_supported = early_match and late_match

        return AnalysisResult(
            statement_id=params['statement_id'],
            evaluation=EvaluationResult.SUPPORTED if is_supported else EvaluationResult.REFUTED,
            evidence=evidence,
            explanation=f"Sliding window: early r={early_corr:.3f}, late r={late_corr:.3f}",
            code_used="sliding_window_correlation",
            confidence=0.8
        )

    def _calculate_composite_score(self, system: str) -> np.ndarray:
        """Calculate composite scores for protein systems"""
        if system == 'vatpase':
            genes = [g for g in self.adata.var_names if g.startswith('ATP6V')]
        elif system == 'proteasome':
            genes = [g for g in self.adata.var_names if 'PSM' in g]
        else:
            genes = []

        if not genes:
            return np.zeros(self.adata.n_obs)

        return np.mean(self.adata[:, genes].X, axis=1)

    def _create_unsure_result(self, statement_id: str) -> AnalysisResult:
        """Create an UNSURE result when analysis cannot be performed"""
        return AnalysisResult(
            statement_id=statement_id,
            evaluation=EvaluationResult.UNSURE,
            evidence={},
            explanation="Unable to perform analysis - missing data or proteins",
            code_used="none",
            confidence=0.0
        )

    def generate_report(self, output_path: str):
        """Generate a comprehensive report of all analyses"""
        report = []
        report.append("# Automated Analysis Report\n")
        report.append(f"## Dataset: {self.data_path}\n")
        report.append(f"## Metadata\n")
        for key, value in self.metadata.items():
            report.append(f"- {key}: {value}\n")

        report.append("\n## Analysis Results\n")
        for statement_id, result in self.results.items():
            report.append(f"\n### {statement_id}\n")
            report.append(f"**Evaluation**: {result.evaluation.value}\n")
            report.append(f"**Confidence**: {result.confidence:.2f}\n")
            report.append(f"**Explanation**: {result.explanation}\n")
            report.append(f"**Evidence**:\n")
            for key, value in result.evidence.items():
                report.append(f"- {key}: {value}\n")

        with open(output_path, 'w') as f:
            f.writelines(report)

        logger.info(f"Report generated: {output_path}")

    def run_finding_group_analysis(self, group_config: Dict):
        """Run analysis for an entire finding group"""
        logger.info(f"Starting analysis for finding group: {group_config['name']}")

        for statement in group_config['statements']:
            result = self.analyze_statement(
                statement['type'],
                statement['params']
            )
            self.results[statement['params']['statement_id']] = result
            logger.info(f"Completed: {statement['params']['statement_id']} - {result.evaluation.value}")

        return self.results
```
