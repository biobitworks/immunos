---
source: /Users/byron/projects/research/experiments/proteomics-claims-2025/archive/01_initial_research/test_calculations.py
relative: research/experiments/proteomics-claims-2025/archive/01_initial_research/test_calculations.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
Unit Tests for Key Calculations
================================
Tests critical calculations to ensure accuracy and reproducibility.

Author: Bioinformatics Finding Group Evaluation Framework
Date: 2024
"""

import sys
import os
import unittest
import numpy as np
from scipy import stats
from scipy.stats import mannwhitneyu, ttest_ind
from statsmodels.stats.multitest import multipletests

sys.path.append('/Users/byron/project_plan')

from config import load_data, get_tau_groups, DATA_SPECS
from master_analysis_with_fdr import apply_fdr_correction


class TestFoldChangeCalculations(unittest.TestCase):
    """Test fold change and log2 fold change calculations."""

    def test_fold_change_basic(self):
        """Test basic fold change calculation."""
        # Test data
        group1 = np.array([10, 12, 11, 13, 9])
        group2 = np.array([5, 6, 4, 5, 5])

        # Calculate fold change
        mean1 = np.mean(group1)
        mean2 = np.mean(group2)
        fold_change = mean1 / mean2

        # Expected: mean1=11, mean2=5, FC=2.2
        self.assertAlmostEqual(mean1, 11.0, places=5)
        self.assertAlmostEqual(mean2, 5.0, places=5)
        self.assertAlmostEqual(fold_change, 2.2, places=5)

    def test_log2_fold_change(self):
        """Test log2 fold change calculation."""
        # Test with known values
        fold_changes = [2, 4, 8, 0.5, 0.25, 1]
        expected_log2 = [1, 2, 3, -1, -2, 0]

        for fc, expected in zip(fold_changes, expected_log2):
            log2_fc = np.log2(fc)
            self.assertAlmostEqual(log2_fc, expected, places=5)

    def test_fold_change_with_zeros(self):
        """Test fold change when denominator is zero."""
        group1 = np.array([10, 12, 11])
        group2 = np.array([0, 0, 0])

        mean1 = np.mean(group1)
        mean2 = np.mean(group2)

        # Should handle division by zero
        if mean2 == 0:
            fold_change = np.inf if mean1 > 0 else 1.0

        self.assertEqual(mean2, 0)
        self.assertEqual(fold_change, np.inf)

    def test_sqstm1_calculation(self):
        """Test SQSTM1 fold change calculation with actual values."""
        # Values from our analysis
        tau_pos_mean = 14.159
        tau_neg_mean = 10.746

        fold_change = tau_pos_mean / tau_neg_mean
        log2_fc = np.log2(fold_change)

        # Should be ~1.318
        self.assertAlmostEqual(fold_change, 1.318, places=2)
        self.assertAlmostEqual(log2_fc, 0.398, places=2)

        # Compare to claimed value
        claimed = 10.7
        discrepancy = claimed / fold_change
        self.assertAlmostEqual(discrepancy, 8.12, places=1)


class TestStatisticalTests(unittest.TestCase):
    """Test statistical test implementations."""

    def test_mann_whitney_u(self):
        """Test Mann-Whitney U test."""
        # Create two groups with known difference
        np.random.seed(42)
        group1 = np.random.normal(10, 2, 30)
        group2 = np.random.normal(8, 2, 30)

        stat, pval = mannwhitneyu(group1, group2)

        # Should detect significant difference
        self.assertLess(pval, 0.05)
        self.assertGreater(stat, 0)

    def test_t_test_comparison(self):
        """Compare t-test and Welch's t-test."""
        # Equal variance case
        np.random.seed(42)
        group1 = np.random.normal(10, 2, 30)
        group2 = np.random.normal(8, 2, 30)

        t_stat, t_pval = ttest_ind(group1, group2)
        welch_stat, welch_pval = ttest_ind(group1, group2, equal_var=False)

        # Results should be similar for equal variance
        self.assertAlmostEqual(t_pval, welch_pval, places=2)

        # Unequal variance case
        group3 = np.random.normal(8, 5, 30)  # Higher variance

        t_stat2, t_pval2 = ttest_ind(group1, group3)
        welch_stat2, welch_pval2 = ttest_ind(group1, group3, equal_var=False)

        # Welch's should be more conservative
        self.assertLessEqual(welch_stat2, t_stat2)

    def test_cohens_d(self):
        """Test Cohen's d effect size calculation."""
        # Large effect size
        group1 = np.array([10, 12, 11, 13, 14])
        group2 = np.array([5, 6, 4, 3, 2])

        # Calculate Cohen's d
        mean1 = np.mean(group1)
        mean2 = np.mean(group2)
        pooled_std = np.sqrt(((len(group1)-1)*np.std(group1, ddof=1)**2 +
                              (len(group2)-1)*np.std(group2, ddof=1)**2) /
                             (len(group1) + len(group2) - 2))
        cohens_d = (mean1 - mean2) / pooled_std

        # Should show large effect size (>0.8)
        self.assertGreater(cohens_d, 0.8)


class TestFDRCorrection(unittest.TestCase):
    """Test FDR correction implementation."""

    def test_fdr_correction_basic(self):
        """Test basic FDR correction."""
        # Test p-values
        pvals = [0.001, 0.01, 0.03, 0.05, 0.1, 0.5, 0.9]

        result = apply_fdr_correction(pvals, alpha=0.05)

        # Check structure
        self.assertIn('pvals_corrected', result)
        self.assertIn('significant', result)
        self.assertIn('n_significant', result)

        # First few should be significant
        self.assertTrue(result['significant'][0])
        self.assertTrue(result['significant'][1])

        # Last few should not be significant
        self.assertFalse(result['significant'][-1])
        self.assertFalse(result['significant'][-2])

    def test_fdr_with_nan(self):
        """Test FDR correction with NaN values."""
        pvals = [0.001, np.nan, 0.03, 0.05, np.nan, 0.5]

        result = apply_fdr_correction(pvals, alpha=0.05)

        # Should handle NaN values
        self.assertEqual(len(result['pvals_corrected']), len(pvals))
        self.assertTrue(np.isnan(result['pvals_corrected'][1]))
        self.assertTrue(np.isnan(result['pvals_corrected'][4]))

    def test_fdr_vs_raw(self):
        """Test that FDR is more conservative than raw p-values."""
        np.random.seed(42)
        # Generate many p-values
        pvals = np.random.uniform(0, 1, 100)

        # Count significant before correction
        n_sig_raw = np.sum(pvals < 0.05)

        # Apply FDR
        result = apply_fdr_correction(pvals, alpha=0.05)
        n_sig_fdr = result['n_significant']

        # FDR should be more conservative
        self.assertLessEqual(n_sig_fdr, n_sig_raw)

    def test_benjamini_hochberg(self):
        """Test Benjamini-Hochberg implementation."""
        # Known example
        pvals = [0.01, 0.02, 0.03, 0.04, 0.05]

        rejected, corrected, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')

        # At least first few should be rejected
        self.assertTrue(rejected[0])
        self.assertTrue(rejected[1])

        # Corrected p-values should be higher
        for orig, corr in zip(pvals, corrected):
            self.assertGreaterEqual(corr, orig)


class TestDataIntegrity(unittest.TestCase):
    """Test data loading and integrity."""

    def test_data_loading(self):
        """Test that data loads correctly."""
        adata = load_data()

        # Check dimensions
        self.assertEqual(adata.n_obs, DATA_SPECS['n_samples'])
        self.assertEqual(adata.n_vars, DATA_SPECS['n_proteins'])

    def test_tau_groups(self):
        """Test tau group splitting."""
        adata = load_data()
        tau_pos, tau_neg = get_tau_groups(adata)

        # Check counts
        self.assertEqual(sum(tau_pos), DATA_SPECS['n_tau_positive'])
        self.assertEqual(sum(tau_neg), DATA_SPECS['n_tau_negative'])

        # Should be mutually exclusive
        self.assertEqual(sum(tau_pos & tau_neg), 0)

    def test_column_names(self):
        """Test that correct column names are used."""
        adata = load_data()

        # Check required columns exist
        self.assertIn(DATA_SPECS['tau_column'], adata.obs.columns)
        self.assertIn(DATA_SPECS['mc1_column'], adata.obs.columns)
        self.assertIn(DATA_SPECS['pseudotime_column'], adata.obs.columns)

    def test_sqstm1_exists(self):
        """Test that SQSTM1 can be found."""
        adata = load_data()

        sqstm1_mask = adata.var['GeneName'].str.contains('SQSTM1', case=False, na=False)
        self.assertTrue(sqstm1_mask.any(), "SQSTM1 not found in dataset")

        # Get index
        sqstm1_idx = np.where(sqstm1_mask)[0][0]
        self.assertIsNotNone(sqstm1_idx)


class TestProteinGroups(unittest.TestCase):
    """Test protein group analysis."""

    def test_proteasome_proteins(self):
        """Test that proteasome proteins are found."""
        adata = load_data()

        proteasome_proteins = ['PSMA1', 'PSMB1', 'PSMC1', 'PSMD1']
        found = 0

        for protein in proteasome_proteins:
            mask = adata.var['GeneName'].str.contains(protein, case=False, na=False)
            if mask.any():
                found += 1

        # Should find most proteasome proteins
        self.assertGreater(found, len(proteasome_proteins) // 2)

    def test_autophagy_proteins(self):
        """Test autophagy protein detection."""
        adata = load_data()

        autophagy_proteins = ['SQSTM1', 'NBR1', 'MAP1LC3B', 'BECN1']
        found = 0

        for protein in autophagy_proteins:
            mask = adata.var['GeneName'].str.contains(protein, case=False, na=False)
            if mask.any():
                found += 1

        # Should find at least some autophagy proteins
        self.assertGreater(found, 0)


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and error handling."""

    def test_empty_array_handling(self):
        """Test handling of empty arrays."""
        empty = np.array([])

        # Should not crash
        mean = np.mean(empty) if len(empty) > 0 else 0
        self.assertEqual(mean, 0)

    def test_single_value_groups(self):
        """Test statistics with single values."""
        group1 = np.array([10])
        group2 = np.array([5])

        # Mann-Whitney should handle this
        try:
            stat, pval = mannwhitneyu(group1, group2)
            # Should complete without error
            self.assertIsNotNone(pval)
        except ValueError:
            # Some versions may raise error for n=1
            pass

    def test_identical_groups(self):
        """Test statistics with identical groups."""
        group1 = np.array([5, 5, 5, 5, 5])
        group2 = np.array([5, 5, 5, 5, 5])

        stat, pval = mannwhitneyu(group1, group2)

        # Should show no difference
        self.assertAlmostEqual(pval, 1.0, places=1)


# ============================================================================
# PERFORMANCE TESTS
# ============================================================================

class TestPerformance(unittest.TestCase):
    """Test performance and efficiency."""

    def test_fdr_performance(self):
        """Test FDR correction performance with large arrays."""
        import time

        # Generate large p-value array
        np.random.seed(42)
        pvals = np.random.uniform(0, 1, 10000)

        start = time.time()
        result = apply_fdr_correction(pvals)
        duration = time.time() - start

        # Should complete in reasonable time
        self.assertLess(duration, 1.0, "FDR correction too slow")

        # Check result validity
        self.assertEqual(len(result['pvals_corrected']), len(pvals))


# ============================================================================
# MAIN TEST RUNNER
# ============================================================================

def run_tests():
    """Run all tests and generate report."""
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    # Add all test classes
    suite.addTests(loader.loadTestsFromTestCase(TestFoldChangeCalculations))
    suite.addTests(loader.loadTestsFromTestCase(TestStatisticalTests))
    suite.addTests(loader.loadTestsFromTestCase(TestFDRCorrection))
    suite.addTests(loader.loadTestsFromTestCase(TestDataIntegrity))
    suite.addTests(loader.loadTestsFromTestCase(TestProteinGroups))
    suite.addTests(loader.loadTestsFromTestCase(TestEdgeCases))
    suite.addTests(loader.loadTestsFromTestCase(TestPerformance))

    # Run tests with verbose output
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)

    # Generate summary
    print("\n" + "="*60)
    print("TEST SUMMARY")
    print("="*60)
    print(f"Tests run: {result.testsRun}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print(f"Success rate: {(result.testsRun - len(result.failures) - len(result.errors))/result.testsRun*100:.1f}%")

    if result.failures:
        print("\nFailed tests:")
        for test, traceback in result.failures:
            print(f"  - {test}")

    if result.errors:
        print("\nTests with errors:")
        for test, traceback in result.errors:
            print(f"  - {test}")

    return result


if __name__ == "__main__":
    print("="*60)
    print("RUNNING UNIT TESTS FOR KEY CALCULATIONS")
    print("="*60)
    result = run_tests()
```
