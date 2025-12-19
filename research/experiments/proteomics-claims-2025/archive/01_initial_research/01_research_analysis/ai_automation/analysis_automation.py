"""
Analysis Automation Module
Provides automated analysis functions for proteomic data
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional
from scipy import stats
from scipy.stats import pearsonr, spearmanr, ttest_ind
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
from statsmodels.formula.api import ols
import logging

logger = logging.getLogger(__name__)

class AnalysisAutomation:
    """Automated analysis functions for finding group evaluations"""

    def __init__(self, adata):
        self.adata = adata
        self.results_cache = {}

    def covariate_controlled_de(self, covariates: List[str] = ['age', 'PMI', 'PatientID']) -> pd.DataFrame:
        """
        Perform differential expression with covariate control
        Returns DataFrame with log2FC, p-values, and FDR
        """
        logger.info(f"Starting covariate-controlled DE with covariates: {covariates}")

        results = []

        # Check if covariates exist
        available_covariates = [c for c in covariates if c in self.adata.obs.columns]
        if len(available_covariates) < len(covariates):
            logger.warning(f"Missing covariates: {set(covariates) - set(available_covariates)}")

        for gene_idx, gene in enumerate(self.adata.var_names):
            if gene_idx % 500 == 0:
                logger.info(f"Processing gene {gene_idx}/{self.adata.n_vars}")

            try:
                # Prepare data
                df = pd.DataFrame({
                    'expression': self.adata[:, gene].X.flatten(),
                    'tau_status': self.adata.obs['tau_status'].astype(str)
                })

                # Add available covariates
                for cov in available_covariates:
                    df[cov] = self.adata.obs[cov]

                # Build formula
                formula = 'expression ~ tau_status'
                for cov in available_covariates:
                    if cov == 'PatientID':
                        formula += f' + C({cov})'  # Categorical
                    else:
                        formula += f' + {cov}'  # Continuous

                # Fit model
                model = ols(formula, data=df)
                fit = model.fit()

                # Extract tau coefficient
                if 'tau_status[T.positive]' in fit.params:
                    coef = fit.params['tau_status[T.positive]']
                    p_val = fit.pvalues['tau_status[T.positive]']
                else:
                    # Alternative parameterization
                    coef = -fit.params.get('tau_status[T.negative]', 0)
                    p_val = fit.pvalues.get('tau_status[T.negative]', 1)

                # Calculate simple fold change for comparison
                tau_pos = self.adata[self.adata.obs['tau_status'] == 'positive'][:, gene].X.mean()
                tau_neg = self.adata[self.adata.obs['tau_status'] == 'negative'][:, gene].X.mean()
                simple_fc = tau_pos - tau_neg  # log2 scale

                results.append({
                    'gene': gene,
                    'log2FC_adjusted': coef,
                    'log2FC_simple': simple_fc,
                    'p_value': p_val
                })

            except Exception as e:
                logger.warning(f"Failed to analyze {gene}: {str(e)}")
                results.append({
                    'gene': gene,
                    'log2FC_adjusted': np.nan,
                    'log2FC_simple': np.nan,
                    'p_value': 1.0
                })

        # Convert to DataFrame
        results_df = pd.DataFrame(results)

        # FDR correction
        valid_p = ~np.isnan(results_df['p_value'])
        results_df.loc[valid_p, 'FDR'] = multipletests(
            results_df.loc[valid_p, 'p_value'],
            method='fdr_bh'
        )[1]

        # Sort by adjusted p-value
        results_df = results_df.sort_values('FDR')

        # Cache results
        self.results_cache['de_results'] = results_df

        logger.info(f"DE analysis complete. Significant genes (FDR<0.05): {sum(results_df['FDR'] < 0.05)}")

        return results_df

    def calculate_protein_correlations(self, proteins: List[str],
                                     target: str = 'SQSTM1',
                                     method: str = 'pearson') -> pd.DataFrame:
        """
        Calculate correlations between proteins and a target
        """
        results = []

        for protein in proteins:
            if protein not in self.adata.var_names:
                logger.warning(f"Protein {protein} not found in dataset")
                continue

            expr = self.adata[:, protein].X.flatten()

            if target == 'pseudotime':
                target_data = self.adata.obs['pseudotime'].values
            elif target == 'MC1':
                target_data = self.adata.obs['MC1'].values
            elif target in self.adata.var_names:
                target_data = self.adata[:, target].X.flatten()
            else:
                logger.warning(f"Target {target} not found")
                continue

            if method == 'pearson':
                corr, p_val = pearsonr(expr, target_data)
            else:
                corr, p_val = spearmanr(expr, target_data)

            results.append({
                'protein': protein,
                'target': target,
                'correlation': corr,
                'p_value': p_val,
                'method': method
            })

        results_df = pd.DataFrame(results)

        # Bonferroni correction
        if len(results_df) > 0:
            bonferroni_alpha = 0.05 / len(results_df)
            results_df['bonferroni_significant'] = results_df['p_value'] < bonferroni_alpha

        return results_df

    def segmented_regression(self, x_var: str, y_var: str,
                            initial_breakpoint: float = None) -> Dict:
        """
        Perform segmented regression to find breakpoint
        """
        from scipy.optimize import curve_fit

        # Get x data
        if x_var in self.adata.obs.columns:
            x_data = self.adata.obs[x_var].values
        else:
            logger.error(f"X variable {x_var} not found")
            return {}

        # Get y data
        if y_var in self.adata.var_names:
            y_data = self.adata[:, y_var].X.flatten()
        elif y_var == 'vatpase_score':
            y_data = self._calculate_vatpase_score()
        elif y_var == 'proteasome_score':
            y_data = self._calculate_proteasome_score()
        else:
            logger.error(f"Y variable {y_var} not found")
            return {}

        # Define piecewise function
        def piecewise_linear(x, x0, y0, b1, b2):
            return np.piecewise(x, [x < x0],
                              [lambda x: y0 + b1*(x-x0),
                               lambda x: y0 + b2*(x-x0)])

        try:
            # Initial parameters
            if initial_breakpoint is None:
                initial_breakpoint = np.median(x_data)

            p0 = [initial_breakpoint, np.mean(y_data), 0, 0]

            # Fit
            popt, pcov = curve_fit(piecewise_linear, x_data, y_data, p0=p0)

            # Calculate R-squared
            y_pred = piecewise_linear(x_data, *popt)
            ss_res = np.sum((y_data - y_pred) ** 2)
            ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
            r2 = 1 - (ss_res / ss_tot)

            result = {
                'breakpoint': popt[0],
                'intercept': popt[1],
                'slope_before': popt[2],
                'slope_after': popt[3],
                'r_squared': r2,
                'parameters': popt,
                'covariance': pcov
            }

            logger.info(f"Segmented regression: breakpoint={popt[0]:.3f}, R²={r2:.3f}")

            return result

        except Exception as e:
            logger.error(f"Segmented regression failed: {str(e)}")
            return {}

    def sliding_window_correlation(self, protein1: str, protein2: str,
                                  window_size: int = 20,
                                  step_size: int = 1) -> pd.DataFrame:
        """
        Calculate sliding window correlations along pseudotime
        """
        if protein1 not in self.adata.var_names or protein2 not in self.adata.var_names:
            logger.error(f"Proteins not found: {protein1}, {protein2}")
            return pd.DataFrame()

        # Sort by pseudotime
        sorted_idx = np.argsort(self.adata.obs['pseudotime'])
        sorted_data = self.adata[sorted_idx]

        results = []

        for i in range(0, len(sorted_data) - window_size + 1, step_size):
            window = sorted_data[i:i+window_size]

            expr1 = window[:, protein1].X.flatten()
            expr2 = window[:, protein2].X.flatten()

            corr, p_val = pearsonr(expr1, expr2)

            results.append({
                'window_start': i,
                'window_end': i + window_size,
                'pseudotime_mean': window.obs['pseudotime'].mean(),
                'pseudotime_min': window.obs['pseudotime'].min(),
                'pseudotime_max': window.obs['pseudotime'].max(),
                'correlation': corr,
                'p_value': p_val
            })

        results_df = pd.DataFrame(results)

        # Calculate phase statistics
        if len(results_df) > 0:
            early_phase = results_df[results_df['pseudotime_mean'] < 0.33]
            late_phase = results_df[results_df['pseudotime_mean'] > 0.67]

            summary = {
                'early_mean_correlation': early_phase['correlation'].mean() if len(early_phase) > 0 else np.nan,
                'late_mean_correlation': late_phase['correlation'].mean() if len(late_phase) > 0 else np.nan,
                'total_windows': len(results_df)
            }

            # Test for trend
            if len(results_df) > 2:
                trend_corr, trend_p = pearsonr(
                    results_df['pseudotime_mean'],
                    results_df['correlation']
                )
                summary['trend_correlation'] = trend_corr
                summary['trend_p_value'] = trend_p

            logger.info(f"Sliding window summary: {summary}")
            results_df.attrs['summary'] = summary

        return results_df

    def analyze_biphasic_behavior(self, system: str = 'vatpase') -> Dict:
        """
        Analyze biphasic behavior of protein systems
        """
        # Calculate system score
        if system == 'vatpase':
            score = self._calculate_vatpase_score()
        elif system == 'proteasome':
            score = self._calculate_proteasome_score()
        else:
            logger.error(f"Unknown system: {system}")
            return {}

        pseudotime = self.adata.obs['pseudotime'].values

        # Find breakpoint using segmented regression
        seg_result = self.segmented_regression('pseudotime', f'{system}_score')

        if not seg_result:
            return {}

        # Calculate phase statistics
        breakpoint = seg_result['breakpoint']

        early_mask = pseudotime < breakpoint
        late_mask = pseudotime >= breakpoint

        early_score = score[early_mask]
        late_score = score[late_mask]

        # Test for significant difference
        t_stat, p_val = ttest_ind(early_score, late_score)

        result = {
            'system': system,
            'breakpoint': breakpoint,
            'slope_early': seg_result['slope_before'],
            'slope_late': seg_result['slope_after'],
            'mean_early': np.mean(early_score),
            'mean_late': np.mean(late_score),
            't_statistic': t_stat,
            'p_value': p_val,
            'n_early': len(early_score),
            'n_late': len(late_score)
        }

        logger.info(f"Biphasic analysis for {system}: breakpoint={breakpoint:.3f}")

        return result

    def _calculate_vatpase_score(self) -> np.ndarray:
        """Calculate V-ATPase composite score"""
        vatpase_genes = [g for g in self.adata.var_names if g.startswith('ATP6V')]

        if not vatpase_genes:
            logger.warning("No V-ATPase genes found")
            return np.zeros(self.adata.n_obs)

        logger.info(f"Calculating V-ATPase score from {len(vatpase_genes)} genes")
        return np.mean(self.adata[:, vatpase_genes].X, axis=1)

    def _calculate_proteasome_score(self) -> np.ndarray:
        """Calculate proteasome composite score"""
        # Try multiple patterns for proteasome subunits
        proteasome_patterns = ['PSM', 'PROT']
        proteasome_genes = []

        for pattern in proteasome_patterns:
            proteasome_genes.extend([g for g in self.adata.var_names if pattern in g])

        # Remove duplicates
        proteasome_genes = list(set(proteasome_genes))

        if not proteasome_genes:
            logger.warning("No proteasome genes found")
            return np.zeros(self.adata.n_obs)

        logger.info(f"Calculating proteasome score from {len(proteasome_genes)} genes")
        return np.mean(self.adata[:, proteasome_genes].X, axis=1)

    def cohen_d(self, group1: np.ndarray, group2: np.ndarray) -> float:
        """Calculate Cohen's d effect size"""
        n1, n2 = len(group1), len(group2)
        var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
        pooled_std = np.sqrt(((n1-1)*var1 + (n2-1)*var2) / (n1+n2-2))

        if pooled_std == 0:
            return 0

        return (np.mean(group1) - np.mean(group2)) / pooled_std

    def analyze_protein_by_mc1_groups(self, protein: str,
                                     mc1_threshold_low: float = 2.5,
                                     mc1_threshold_high: float = 3.0) -> Dict:
        """
        Analyze protein expression by MC1 groups
        """
        if protein not in self.adata.var_names:
            logger.error(f"Protein {protein} not found")
            return {}

        expr = self.adata[:, protein].X.flatten()
        mc1 = self.adata.obs['MC1'].values

        # Define groups
        low_mc1_mask = mc1 < mc1_threshold_low
        high_mc1_mask = mc1 >= mc1_threshold_high

        low_mc1_expr = expr[low_mc1_mask]
        high_mc1_expr = expr[high_mc1_mask]

        # Statistics
        t_stat, p_val = ttest_ind(low_mc1_expr, high_mc1_expr)
        cohen_d = self.cohen_d(low_mc1_expr, high_mc1_expr)

        result = {
            'protein': protein,
            'low_mc1_mean': np.mean(low_mc1_expr),
            'low_mc1_std': np.std(low_mc1_expr),
            'low_mc1_n': len(low_mc1_expr),
            'high_mc1_mean': np.mean(high_mc1_expr),
            'high_mc1_std': np.std(high_mc1_expr),
            'high_mc1_n': len(high_mc1_expr),
            't_statistic': t_stat,
            'p_value': p_val,
            'cohens_d': cohen_d
        }

        logger.info(f"{protein} by MC1 groups: t={t_stat:.3f}, p={p_val:.3e}, d={cohen_d:.3f}")

        return result