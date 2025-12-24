#!/usr/bin/env python3
"""
IMMUNOS Statistical Data Logger
================================

Captures ALL interactions for rigor assessment with p-values and confidence intervals.

Tracks:
- Every task, model used, tokens, response time
- Validation results, confidence scores
- Statistical significance (p-values)
- Rigor metrics for scientific credibility
"""
import sqlite3
import numpy as np
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Any, Optional, Tuple
from dataclasses import dataclass

try:
    from scipy import stats
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    print("⚠️  scipy not installed - p-value calculations will be limited")


@dataclass
class RigorMetrics:
    """Statistical rigor metrics"""
    mean_accuracy: float
    std_accuracy: float
    confidence_interval_95: Tuple[float, float]
    p_value: float
    sample_size: int
    validation_rate: float


class ImmunosDataLogger:
    """Comprehensive data logger with statistical analysis"""

    def __init__(self, db_path: str = None):
        self.db_path = db_path or str(Path.home() / ".immunos" / "db" / "immunos.db")
        self._init_db()

    def _init_db(self):
        """Initialize statistics tables"""
        conn = sqlite3.connect(self.db_path)
        conn.execute('''
            CREATE TABLE IF NOT EXISTS interaction_log (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                timestamp TEXT NOT NULL,
                task_type TEXT NOT NULL,
                task_classification TEXT,
                model_used TEXT NOT NULL,
                provider TEXT NOT NULL,
                tokens_prompt INTEGER,
                tokens_response INTEGER,
                tokens_total INTEGER,
                response_time_ms INTEGER,
                confidence_score REAL,
                validation_result INTEGER,
                validation_method TEXT,
                error_detected INTEGER DEFAULT 0,
                hallucination_detected INTEGER DEFAULT 0,
                routing_reason TEXT,
                context TEXT
            )
        ''')

        conn.execute('''
            CREATE TABLE IF NOT EXISTS rigor_statistics (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                calculated_at TEXT NOT NULL,
                time_window_hours INTEGER,
                sample_size INTEGER,
                mean_accuracy REAL,
                std_accuracy REAL,
                ci_lower REAL,
                ci_upper REAL,
                p_value REAL,
                validation_rate REAL,
                mean_confidence REAL,
                mean_response_time_ms REAL
            )
        ''')

        conn.commit()
        conn.close()

    def log_interaction(self,
                       task_type: str,
                       model_used: str,
                       provider: str,
                       tokens_prompt: int = 0,
                       tokens_response: int = 0,
                       response_time_ms: int = 0,
                       confidence_score: float = None,
                       validation_result: bool = None,
                       validation_method: str = None,
                       task_classification: str = None,
                       error_detected: bool = False,
                       hallucination_detected: bool = False,
                       routing_reason: str = None,
                       context: str = None):
        """
        Log every interaction for statistical analysis

        Args:
            task_type: Type of task (code_review, reasoning, etc.)
            model_used: Which model was used
            provider: 'claude' or 'ollama'
            tokens_prompt: Input tokens
            tokens_response: Output tokens
            response_time_ms: Response time in milliseconds
            confidence_score: Confidence in result (0-1)
            validation_result: True if validated successfully
            validation_method: How it was validated
            task_classification: 'routine', 'creative', 'critical'
            error_detected: Whether an error was caught
            hallucination_detected: Whether hallucination was caught
            routing_reason: Why this model was chosen
            context: Additional context
        """
        conn = sqlite3.connect(self.db_path)
        conn.execute('''
            INSERT INTO interaction_log (
                timestamp, task_type, task_classification, model_used, provider,
                tokens_prompt, tokens_response, tokens_total, response_time_ms,
                confidence_score, validation_result, validation_method,
                error_detected, hallucination_detected, routing_reason, context
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (
            datetime.now().isoformat(),
            task_type,
            task_classification,
            model_used,
            provider,
            tokens_prompt,
            tokens_response,
            tokens_prompt + tokens_response,
            response_time_ms,
            confidence_score,
            1 if validation_result else 0 if validation_result is not None else None,
            validation_method,
            1 if error_detected else 0,
            1 if hallucination_detected else 0,
            routing_reason,
            context
        ))
        conn.commit()
        conn.close()

    def calculate_rigor_metrics(self, hours: int = 24) -> RigorMetrics:
        """
        Calculate statistical rigor metrics

        Returns:
            RigorMetrics with p-values and confidence intervals
        """
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row

        # Get data from last N hours
        cursor = conn.execute('''
            SELECT confidence_score, validation_result
            FROM interaction_log
            WHERE datetime(timestamp) >= datetime('now', '-' || ? || ' hours')
            AND confidence_score IS NOT NULL
        ''', (hours,))

        rows = cursor.fetchall()
        conn.close()

        if not rows or len(rows) < 2:
            return RigorMetrics(
                mean_accuracy=0.0,
                std_accuracy=0.0,
                confidence_interval_95=(0.0, 0.0),
                p_value=1.0,
                sample_size=len(rows) if rows else 0,
                validation_rate=0.0
            )

        # Extract data
        confidence_scores = [r['confidence_score'] for r in rows if r['confidence_score'] is not None]
        validation_results = [r['validation_result'] for r in rows if r['validation_result'] is not None]

        # Calculate statistics
        mean_confidence = np.mean(confidence_scores)
        std_confidence = np.std(confidence_scores, ddof=1)

        # Calculate 95% CI
        if HAS_SCIPY:
            ci = stats.t.interval(0.95, len(confidence_scores)-1,
                                 loc=mean_confidence,
                                 scale=stats.sem(confidence_scores))

            # One-sample t-test: H0: mean = 0.5 (random), H1: mean > 0.5 (better than random)
            t_statistic, p_value = stats.ttest_1samp(confidence_scores, 0.5)
            p_value = p_value / 2  # One-tailed test
        else:
            # Simple CI without scipy
            sem = std_confidence / np.sqrt(len(confidence_scores))
            ci = (mean_confidence - 1.96 * sem, mean_confidence + 1.96 * sem)

            # Simple t-test approximation
            t_statistic = (mean_confidence - 0.5) / (std_confidence / np.sqrt(len(confidence_scores)))
            p_value = 0.001 if abs(t_statistic) > 3 else 0.05

        validation_rate = sum(validation_results) / len(validation_results) if validation_results else 0.0

        metrics = RigorMetrics(
            mean_accuracy=mean_confidence,
            std_accuracy=std_confidence,
            confidence_interval_95=ci,
            p_value=max(0.0, p_value),  # Ensure non-negative
            sample_size=len(confidence_scores),
            validation_rate=validation_rate
        )

        # Save to database
        self._save_rigor_stats(metrics, hours)

        return metrics

    def _save_rigor_stats(self, metrics: RigorMetrics, hours: int):
        """Save calculated rigor statistics"""
        conn = sqlite3.connect(self.db_path)
        conn.execute('''
            INSERT INTO rigor_statistics (
                calculated_at, time_window_hours, sample_size,
                mean_accuracy, std_accuracy, ci_lower, ci_upper,
                p_value, validation_rate, mean_confidence, mean_response_time_ms
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (
            datetime.now().isoformat(),
            hours,
            metrics.sample_size,
            metrics.mean_accuracy,
            metrics.std_accuracy,
            metrics.confidence_interval_95[0],
            metrics.confidence_interval_95[1],
            metrics.p_value,
            metrics.validation_rate,
            metrics.mean_accuracy,
            0  # placeholder for mean response time
        ))
        conn.commit()
        conn.close()

    def get_rigor_summary(self) -> dict:
        """Get human-readable rigor summary"""
        metrics = self.calculate_rigor_metrics(hours=24)

        # Interpret p-value
        if metrics.p_value < 0.001:
            significance = "Highly significant (p < 0.001) - IMMUNOS is performing far better than random"
        elif metrics.p_value < 0.01:
            significance = "Very significant (p < 0.01) - Strong evidence of rigor"
        elif metrics.p_value < 0.05:
            significance = "Significant (p < 0.05) - Statistically rigorous"
        else:
            significance = f"Not significant (p = {metrics.p_value:.3f}) - Need more data"

        return {
            'summary': {
                'mean_accuracy': f"{metrics.mean_accuracy:.1%}",
                'confidence_interval': f"[{metrics.confidence_interval_95[0]:.1%}, {metrics.confidence_interval_95[1]:.1%}]",
                'p_value': f"{metrics.p_value:.4f}",
                'significance': significance,
                'sample_size': metrics.sample_size,
                'validation_rate': f"{metrics.validation_rate:.1%}"
            },
            'raw': metrics
        }


if __name__ == '__main__':
    # CLI interface for checking rigor
    logger = ImmunosDataLogger()

    print("\n" + "="*70)
    print("IMMUNOS Rigor Assessment")
    print("="*70 + "\n")

    summary = logger.get_rigor_summary()

    print(f"Sample Size: {summary['summary']['sample_size']}")
    print(f"Mean Accuracy: {summary['summary']['mean_accuracy']}")
    print(f"95% CI: {summary['summary']['confidence_interval']}")
    print(f"P-Value: {summary['summary']['p_value']}")
    print(f"\n{summary['summary']['significance']}")
    print(f"\nValidation Rate: {summary['summary']['validation_rate']}")
    print()
