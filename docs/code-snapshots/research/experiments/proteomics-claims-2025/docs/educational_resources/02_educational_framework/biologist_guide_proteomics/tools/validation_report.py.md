---
source: /Users/byron/projects/research/experiments/proteomics-claims-2025/docs/educational_resources/02_educational_framework/biologist_guide_proteomics/tools/validation_report.py
relative: research/experiments/proteomics-claims-2025/docs/educational_resources/02_educational_framework/biologist_guide_proteomics/tools/validation_report.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
Protein Validation Report Generator

Creates comprehensive validation reports with visualizations and detailed analysis
of protein annotation quality, UniProt consistency, and UPS interactions.

Author: Proteomics Analysis Framework
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Optional
import json
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')


class ValidationReportGenerator:
    """
    Generate comprehensive validation reports with visualizations
    """

    def __init__(self):
        """Initialize report generator with styling"""
        plt.style.use('seaborn-v0_8-darkgrid')
        sns.set_palette("husl")

    def create_validation_report(self, validation_df: pd.DataFrame,
                               output_dir: str = './',
                               create_plots: bool = True) -> Dict:
        """
        Create comprehensive validation report with visualizations

        Parameters:
        -----------
        validation_df : pd.DataFrame
            Results from ProteinValidator
        output_dir : str
            Directory to save reports and plots
        create_plots : bool
            Whether to generate visualization plots

        Returns:
        --------
        dict
            Comprehensive report summary
        """
        print("Generating comprehensive validation report...")

        # Calculate comprehensive statistics
        report = self._calculate_comprehensive_stats(validation_df)

        # Generate visualizations if requested
        if create_plots:
            self._create_validation_plots(validation_df, output_dir, report)

        # Generate detailed HTML report
        self._generate_html_report(validation_df, report, output_dir)

        # Generate recommendations
        report['recommendations'] = self._generate_recommendations(validation_df, report)

        # Save JSON summary
        report_file = f"{output_dir}/validation_summary.json"
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2, default=str)

        print(f"Report generated successfully in {output_dir}/")
        return report

    def _calculate_comprehensive_stats(self, validation_df: pd.DataFrame) -> Dict:
        """Calculate comprehensive validation statistics"""

        total_proteins = len(validation_df)

        # Basic validation stats
        valid_count = len(validation_df[validation_df['validation_status'] == 'valid'])
        warning_count = len(validation_df[validation_df['validation_status'] == 'warning'])
        invalid_count = len(validation_df[validation_df['validation_status'] == 'invalid'])

        # UniProt stats
        existing_count = len(validation_df[validation_df['uniprot_exists']])
        reviewed_count = len(validation_df[validation_df['uniprot_reviewed']])
        multiple_ids_count = len(validation_df[validation_df['has_multiple_ids']])

        # Gene name matching
        gene_match_count = len(validation_df[validation_df['gene_name_match']])

        # UPS analysis
        ups_components = len(validation_df[validation_df['is_ups_component']])
        core_ups = len(validation_df[validation_df['ups_category'] == 'core_ups'])
        ups_associated = len(validation_df[validation_df['ups_category'] == 'ups_associated'])

        # Confidence distribution
        high_conf_ups = len(validation_df[validation_df['ups_confidence'] == 'high'])
        medium_conf_ups = len(validation_df[validation_df['ups_confidence'] == 'medium'])
        low_conf_ups = len(validation_df[validation_df['ups_confidence'] == 'low'])

        report = {
            'metadata': {
                'generation_time': datetime.now().isoformat(),
                'total_proteins_analyzed': total_proteins
            },
            'validation_summary': {
                'valid_proteins': valid_count,
                'warning_proteins': warning_count,
                'invalid_proteins': invalid_count,
                'validation_rate': (valid_count + warning_count) / total_proteins * 100
            },
            'uniprot_analysis': {
                'existing_in_uniprot': existing_count,
                'reviewed_entries': reviewed_count,
                'multiple_id_entries': multiple_ids_count,
                'existence_rate': existing_count / total_proteins * 100,
                'review_rate': reviewed_count / total_proteins * 100,
                'multiple_id_rate': multiple_ids_count / total_proteins * 100
            },
            'annotation_consistency': {
                'gene_name_matches': gene_match_count,
                'gene_match_rate': gene_match_count / total_proteins * 100
            },
            'ups_analysis': {
                'total_ups_components': ups_components,
                'core_ups_components': core_ups,
                'ups_associated_proteins': ups_associated,
                'ups_coverage_rate': ups_components / total_proteins * 100,
                'confidence_distribution': {
                    'high': high_conf_ups,
                    'medium': medium_conf_ups,
                    'low': low_conf_ups
                }
            },
            'quality_metrics': {
                'overall_quality_score': self._calculate_quality_score(validation_df),
                'data_completeness': self._calculate_completeness(validation_df),
                'annotation_reliability': self._calculate_reliability(validation_df)
            }
        }

        return report

    def _calculate_quality_score(self, validation_df: pd.DataFrame) -> float:
        """Calculate overall data quality score (0-100)"""
        weights = {
            'valid': 1.0,
            'warning': 0.7,
            'invalid': 0.0
        }

        total_score = 0
        for status, weight in weights.items():
            count = len(validation_df[validation_df['validation_status'] == status])
            total_score += count * weight

        return (total_score / len(validation_df)) * 100

    def _calculate_completeness(self, validation_df: pd.DataFrame) -> float:
        """Calculate data completeness score"""
        # Check for non-null values in key columns
        key_columns = ['uniprot_exists', 'gene_name_match', 'uniprot_reviewed']
        completeness_scores = []

        for col in key_columns:
            if col in validation_df.columns:
                non_null_count = validation_df[col].notna().sum()
                completeness_scores.append(non_null_count / len(validation_df))

        return np.mean(completeness_scores) * 100 if completeness_scores else 0

    def _calculate_reliability(self, validation_df: pd.DataFrame) -> float:
        """Calculate annotation reliability score"""
        # Base score on reviewed entries and gene name matches
        reviewed_score = validation_df['uniprot_reviewed'].sum() / len(validation_df)
        gene_match_score = validation_df['gene_name_match'].sum() / len(validation_df)

        return ((reviewed_score + gene_match_score) / 2) * 100

    def _create_validation_plots(self, validation_df: pd.DataFrame,
                               output_dir: str, report: Dict):
        """Create comprehensive validation visualizations"""

        # Create figure with subplots
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle('Protein Validation Analysis Report', fontsize=16, fontweight='bold')

        # 1. Validation Status Distribution
        ax = axes[0, 0]
        status_counts = validation_df['validation_status'].value_counts()
        colors = {'valid': 'green', 'warning': 'orange', 'invalid': 'red'}
        status_colors = [colors.get(status, 'gray') for status in status_counts.index]

        wedges, texts, autotexts = ax.pie(status_counts.values, labels=status_counts.index,
                                         autopct='%1.1f%%', colors=status_colors)
        ax.set_title('Validation Status Distribution')

        # 2. UniProt Analysis
        ax = axes[0, 1]
        uniprot_data = {
            'Exists in UniProt': validation_df['uniprot_exists'].sum(),
            'Reviewed Entries': validation_df['uniprot_reviewed'].sum(),
            'Multiple IDs': validation_df['has_multiple_ids'].sum()
        }

        bars = ax.bar(uniprot_data.keys(), uniprot_data.values(),
                     color=['skyblue', 'lightgreen', 'salmon'])
        ax.set_title('UniProt Database Analysis')
        ax.set_ylabel('Number of Proteins')
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right')

        # Add value labels on bars
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                   f'{int(height)}', ha='center', va='bottom')

        # 3. UPS Component Analysis
        ax = axes[0, 2]
        ups_data = validation_df['ups_category'].value_counts()
        if not ups_data.empty:
            ups_colors = ['lightcoral', 'lightsalmon', 'lightgray']
            ax.pie(ups_data.values, labels=ups_data.index, autopct='%1.1f%%',
                  colors=ups_colors[:len(ups_data)])
        ax.set_title('UPS Component Classification')

        # 4. Gene Name Matching
        ax = axes[1, 0]
        match_data = {
            'Matches': validation_df['gene_name_match'].sum(),
            'Mismatches': len(validation_df) - validation_df['gene_name_match'].sum()
        }
        ax.bar(match_data.keys(), match_data.values(), color=['lightgreen', 'lightcoral'])
        ax.set_title('Gene Name Consistency')
        ax.set_ylabel('Number of Proteins')

        # 5. Quality Metrics Overview
        ax = axes[1, 1]
        quality_metrics = {
            'Overall Quality': report['quality_metrics']['overall_quality_score'],
            'Data Completeness': report['quality_metrics']['data_completeness'],
            'Annotation Reliability': report['quality_metrics']['annotation_reliability']
        }

        bars = ax.barh(list(quality_metrics.keys()), list(quality_metrics.values()),
                      color=['gold', 'lightblue', 'lightgreen'])
        ax.set_title('Quality Metrics (0-100)')
        ax.set_xlabel('Score')

        # Add value labels
        for i, bar in enumerate(bars):
            width = bar.get_width()
            ax.text(width + 1, bar.get_y() + bar.get_height()/2,
                   f'{width:.1f}', ha='left', va='center')

        # 6. UPS Confidence Distribution
        ax = axes[1, 2]
        ups_subset = validation_df[validation_df['is_ups_component']]
        if not ups_subset.empty:
            conf_counts = ups_subset['ups_confidence'].value_counts()
            conf_colors = {'high': 'darkgreen', 'medium': 'orange', 'low': 'lightcoral'}
            conf_plot_colors = [conf_colors.get(conf, 'gray') for conf in conf_counts.index]
            ax.bar(conf_counts.index, conf_counts.values, color=conf_plot_colors)
            ax.set_title('UPS Classification Confidence')
            ax.set_ylabel('Number of Proteins')
        else:
            ax.text(0.5, 0.5, 'No UPS components found',
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title('UPS Classification Confidence')

        plt.tight_layout()
        plt.savefig(f'{output_dir}/validation_analysis.png', dpi=300, bbox_inches='tight')
        plt.close()

        # Create detailed UPS analysis plot if UPS components exist
        ups_proteins = validation_df[validation_df['is_ups_component']]
        if not ups_proteins.empty:
            self._create_ups_analysis_plot(ups_proteins, output_dir)

    def _create_ups_analysis_plot(self, ups_df: pd.DataFrame, output_dir: str):
        """Create detailed UPS analysis visualization"""

        fig, axes = plt.subplots(2, 2, figsize=(15, 10))
        fig.suptitle('Detailed UPS Component Analysis', fontsize=14, fontweight='bold')

        # 1. UPS Category Distribution
        ax = axes[0, 0]
        category_counts = ups_df['ups_category'].value_counts()
        ax.pie(category_counts.values, labels=category_counts.index, autopct='%1.1f%%')
        ax.set_title('UPS Category Distribution')

        # 2. Confidence vs Category
        ax = axes[0, 1]
        conf_category = pd.crosstab(ups_df['ups_confidence'], ups_df['ups_category'])
        conf_category.plot(kind='bar', ax=ax, stacked=True)
        ax.set_title('Confidence by UPS Category')
        ax.set_xlabel('Confidence Level')
        ax.set_ylabel('Number of Proteins')
        ax.legend(title='UPS Category')
        plt.setp(ax.get_xticklabels(), rotation=0)

        # 3. Validation Status of UPS Components
        ax = axes[1, 0]
        ups_validation = ups_df['validation_status'].value_counts()
        colors = {'valid': 'green', 'warning': 'orange', 'invalid': 'red'}
        ups_colors = [colors.get(status, 'gray') for status in ups_validation.index]
        ax.bar(ups_validation.index, ups_validation.values, color=ups_colors)
        ax.set_title('Validation Status of UPS Components')
        ax.set_ylabel('Number of Proteins')

        # 4. Top UPS Proteins by Gene Symbol
        ax = axes[1, 1]
        if len(ups_df) <= 15:  # Show all if reasonable number
            gene_list = ups_df['gene_symbol'].tolist()
            confidence_list = ups_df['ups_confidence'].tolist()

            # Color by confidence
            color_map = {'high': 'darkgreen', 'medium': 'orange', 'low': 'lightcoral'}
            colors = [color_map.get(conf, 'gray') for conf in confidence_list]

            ax.barh(gene_list, [1] * len(gene_list), color=colors)
            ax.set_title('UPS Components by Confidence')
            ax.set_xlabel('Relative Importance')
        else:
            # Show top genes by confidence
            high_conf = ups_df[ups_df['ups_confidence'] == 'high']['gene_symbol'].head(10)
            if not high_conf.empty:
                ax.barh(high_conf, [1] * len(high_conf), color='darkgreen')
                ax.set_title('Top High-Confidence UPS Components')
                ax.set_xlabel('Relative Importance')

        plt.tight_layout()
        plt.savefig(f'{output_dir}/ups_detailed_analysis.png', dpi=300, bbox_inches='tight')
        plt.close()

    def _generate_html_report(self, validation_df: pd.DataFrame,
                            report: Dict, output_dir: str):
        """Generate detailed HTML report"""

        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>Protein Validation Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                .header {{ background-color: #f0f0f0; padding: 20px; border-radius: 5px; }}
                .section {{ margin: 20px 0; }}
                .metric {{ display: inline-block; margin: 10px; padding: 15px;
                         background-color: #e8f4f8; border-radius: 5px; min-width: 150px; }}
                .metric-value {{ font-size: 24px; font-weight: bold; color: #2c3e50; }}
                .metric-label {{ font-size: 14px; color: #7f8c8d; }}
                .warning {{ color: #e67e22; }}
                .error {{ color: #e74c3c; }}
                .success {{ color: #27ae60; }}
                table {{ border-collapse: collapse; width: 100%; margin: 10px 0; }}
                th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
                th {{ background-color: #f2f2f2; }}
                .issue-list {{ background-color: #ffeaa7; padding: 10px; border-radius: 5px; }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>Protein Validation Analysis Report</h1>
                <p>Generated on: {report['metadata']['generation_time']}</p>
                <p>Total proteins analyzed: {report['metadata']['total_proteins_analyzed']}</p>
            </div>

            <div class="section">
                <h2>Executive Summary</h2>
                <div class="metric">
                    <div class="metric-value success">{report['validation_summary']['valid_proteins']}</div>
                    <div class="metric-label">Valid Proteins</div>
                </div>
                <div class="metric">
                    <div class="metric-value warning">{report['validation_summary']['warning_proteins']}</div>
                    <div class="metric-label">Warning Proteins</div>
                </div>
                <div class="metric">
                    <div class="metric-value error">{report['validation_summary']['invalid_proteins']}</div>
                    <div class="metric-label">Invalid Proteins</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{report['validation_summary']['validation_rate']:.1f}%</div>
                    <div class="metric-label">Validation Rate</div>
                </div>
            </div>

            <div class="section">
                <h2>UniProt Database Analysis</h2>
                <div class="metric">
                    <div class="metric-value">{report['uniprot_analysis']['existence_rate']:.1f}%</div>
                    <div class="metric-label">Existence Rate</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{report['uniprot_analysis']['review_rate']:.1f}%</div>
                    <div class="metric-label">Review Rate</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{report['uniprot_analysis']['multiple_id_rate']:.1f}%</div>
                    <div class="metric-label">Multiple ID Rate</div>
                </div>
            </div>

            <div class="section">
                <h2>UPS (Ubiquitin-Proteasome System) Analysis</h2>
                <div class="metric">
                    <div class="metric-value">{report['ups_analysis']['total_ups_components']}</div>
                    <div class="metric-label">Total UPS Components</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{report['ups_analysis']['core_ups_components']}</div>
                    <div class="metric-label">Core UPS Components</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{report['ups_analysis']['ups_coverage_rate']:.1f}%</div>
                    <div class="metric-label">UPS Coverage Rate</div>
                </div>
            </div>

            <div class="section">
                <h2>Quality Metrics</h2>
                <div class="metric">
                    <div class="metric-value">{report['quality_metrics']['overall_quality_score']:.1f}</div>
                    <div class="metric-label">Overall Quality Score</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{report['quality_metrics']['data_completeness']:.1f}</div>
                    <div class="metric-label">Data Completeness</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{report['quality_metrics']['annotation_reliability']:.1f}</div>
                    <div class="metric-label">Annotation Reliability</div>
                </div>
            </div>
        """

        # Add issues section if there are invalid proteins
        invalid_proteins = validation_df[validation_df['validation_status'] == 'invalid']
        if not invalid_proteins.empty:
            html_content += """
            <div class="section">
                <h2>Critical Issues Requiring Attention</h2>
                <div class="issue-list">
            """

            for _, row in invalid_proteins.head(10).iterrows():  # Show top 10 issues
                html_content += f"""
                    <p><strong>{row['gene_symbol']}</strong> ({row['original_uniprot_id']}):
                    {', '.join(row['validation_issues']) if row['validation_issues'] else 'Unknown issue'}</p>
                """

            if len(invalid_proteins) > 10:
                html_content += f"<p><em>... and {len(invalid_proteins) - 10} more issues</em></p>"

            html_content += "</div></div>"

        # Add UPS components table
        ups_proteins = validation_df[validation_df['is_ups_component']]
        if not ups_proteins.empty:
            html_content += """
            <div class="section">
                <h2>Identified UPS Components</h2>
                <table>
                    <tr>
                        <th>Gene Symbol</th>
                        <th>UniProt ID</th>
                        <th>UPS Category</th>
                        <th>Confidence</th>
                        <th>Validation Status</th>
                    </tr>
            """

            for _, row in ups_proteins.iterrows():
                status_class = {'valid': 'success', 'warning': 'warning', 'invalid': 'error'}.get(
                    row['validation_status'], ''
                )
                html_content += f"""
                    <tr>
                        <td>{row['gene_symbol']}</td>
                        <td>{row['primary_uniprot_id']}</td>
                        <td>{row['ups_category']}</td>
                        <td>{row['ups_confidence']}</td>
                        <td class="{status_class}">{row['validation_status']}</td>
                    </tr>
                """

            html_content += "</table></div>"

        html_content += """
        </body>
        </html>
        """

        with open(f'{output_dir}/validation_report.html', 'w') as f:
            f.write(html_content)

    def _generate_recommendations(self, validation_df: pd.DataFrame,
                                report: Dict) -> List[str]:
        """Generate actionable recommendations based on validation results"""

        recommendations = []

        # Critical issues
        if report['validation_summary']['invalid_proteins'] > 0:
            recommendations.append(
                f"🔴 CRITICAL: Review and correct {report['validation_summary']['invalid_proteins']} "
                f"invalid protein entries before proceeding with analysis"
            )

        # Quality thresholds
        if report['quality_metrics']['overall_quality_score'] < 70:
            recommendations.append(
                "🟡 WARNING: Overall data quality score is below 70%. Consider data cleanup"
            )

        if report['uniprot_analysis']['review_rate'] < 80:
            recommendations.append(
                "🟡 Consider prioritizing Swiss-Prot reviewed entries for higher annotation quality"
            )

        # Multiple IDs handling
        if report['uniprot_analysis']['multiple_id_rate'] > 5:
            recommendations.append(
                f"🔵 INFO: {report['uniprot_analysis']['multiple_id_entries']} proteins have multiple "
                f"UniProt IDs. Consider establishing primary ID selection criteria"
            )

        # UPS analysis recommendations
        if report['ups_analysis']['ups_coverage_rate'] > 20:
            recommendations.append(
                f"✅ EXCELLENT: {report['ups_analysis']['ups_coverage_rate']:.1f}% UPS coverage "
                f"indicates strong representation of protein quality control systems"
            )
        elif report['ups_analysis']['ups_coverage_rate'] > 10:
            recommendations.append(
                f"✅ GOOD: {report['ups_analysis']['ups_coverage_rate']:.1f}% UPS coverage "
                f"provides adequate representation for quality control analysis"
            )
        else:
            recommendations.append(
                f"🟡 LOW UPS coverage ({report['ups_analysis']['ups_coverage_rate']:.1f}%). "
                f"Consider expanding protein list to include more UPS components"
            )

        # Analysis impact recommendations
        gene_mismatch_rate = 100 - report['annotation_consistency']['gene_match_rate']
        if gene_mismatch_rate > 10:
            recommendations.append(
                f"🟡 {gene_mismatch_rate:.1f}% of proteins have gene name mismatches. "
                f"This could impact downstream analysis accuracy"
            )

        return recommendations


def main():
    """
    Example usage of validation report generator
    """
    # This would typically use results from ProteinValidator
    # For demonstration, create sample validation results

    print("Validation Report Generator - Demo")
    print("Note: Run uniprot_analysis.py first to generate validation_df")

    # Example of how to use after running validation:
    # validation_df = pd.read_csv('protein_validation_results.csv')
    # reporter = ValidationReportGenerator()
    # report = reporter.create_validation_report(validation_df, output_dir='validation_reports')


if __name__ == "__main__":
    main()
```
