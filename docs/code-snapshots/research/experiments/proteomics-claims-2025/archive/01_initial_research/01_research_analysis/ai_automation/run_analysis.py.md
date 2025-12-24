---
source: /Users/byron/projects/research/experiments/proteomics-claims-2025/archive/01_initial_research/01_research_analysis/ai_automation/run_analysis.py
relative: research/experiments/proteomics-claims-2025/archive/01_initial_research/01_research_analysis/ai_automation/run_analysis.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
Main orchestrator script for automated bioinformatics analysis
Runs the AI agent on finding groups
"""

import os
import sys
import json
import argparse
import logging
from datetime import datetime
from ai_agent import BioinformaticsAgent, EvaluationResult
from analysis_automation import AnalysisAutomation

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f'analysis_{datetime.now():%Y%m%d_%H%M%S}.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Finding Group 1 Configuration
FINDING_GROUP_1 = {
    'name': 'Late-Stage Mitochondrial Dysregulation',
    'statements': [
        {
            'type': 'differential_expression',
            'params': {
                'statement_id': 'G1_S1',
                'description': 'UPS protein alterations',
                'proteins': ['UPS1', 'UPS2'],  # Will need actual UPS protein names
                'expected_result': 'no_significant_change'
            }
        },
        {
            'type': 'protein_upregulation',
            'params': {
                'statement_id': 'G1_S2',
                'description': 'SQSTM1 upregulation',
                'protein': 'SQSTM1',
                'expected_log2fc': 3.413,
                'tolerance': 0.5
            }
        },
        {
            'type': 'correlation',
            'params': {
                'statement_id': 'G1_S5',
                'description': 'SQSTM1-VDAC1 global correlation',
                'protein1': 'SQSTM1',
                'protein2': 'VDAC1',
                'expected_correlation': 0.0536,
                'tolerance': 0.1
            }
        },
        {
            'type': 'sliding_window',
            'params': {
                'statement_id': 'G1_S6',
                'description': 'SQSTM1-VDAC1 running correlation',
                'protein1': 'SQSTM1',
                'protein2': 'VDAC1',
                'window_size': 20,
                'expected_early': -0.417,
                'expected_late': 0.478
            }
        }
    ]
}

# Finding Group 2 Configuration
FINDING_GROUP_2 = {
    'name': 'Sequential Failure of Proteostasis',
    'statements': [
        {
            'type': 'differential_expression',
            'params': {
                'statement_id': 'G2_S1',
                'description': 'Covariate-controlled DE analysis',
                'covariates': ['age', 'PMI', 'PatientID'],
                'expected_significant_percentage': 36.14
            }
        },
        {
            'type': 'protein_upregulation',
            'params': {
                'statement_id': 'G2_S2',
                'description': 'SQSTM1 top upregulated',
                'protein': 'SQSTM1',
                'expected_log2fc': 3.41,
                'check_if_top': True
            }
        },
        {
            'type': 'segmented_regression',
            'params': {
                'statement_id': 'G2_S5',
                'description': 'V-ATPase score segmented analysis',
                'x_variable': 'MC1',
                'y_variable': 'vatpase_score',
                'expected_breakpoint': 2.831
            }
        },
        {
            'type': 'biphasic',
            'params': {
                'statement_id': 'G2_S6',
                'description': 'V-ATPase biphasic behavior',
                'system': 'vatpase',
                'expected_breakpoint': 0.654
            }
        }
    ]
}

def run_finding_group(agent: BioinformaticsAgent, group_config: dict, output_dir: str):
    """
    Run analysis for a finding group
    """
    logger.info(f"="*60)
    logger.info(f"Starting analysis for: {group_config['name']}")
    logger.info(f"="*60)

    # Create output directory for this group
    group_dir = os.path.join(output_dir, group_config['name'].replace(' ', '_'))
    os.makedirs(group_dir, exist_ok=True)

    # Run analysis
    results = agent.run_finding_group_analysis(group_config)

    # Generate detailed report
    report_path = os.path.join(group_dir, 'analysis_report.md')
    agent.generate_report(report_path)

    # Generate evaluation summary
    summary_path = os.path.join(group_dir, 'evaluation_summary.json')
    summary = {
        'group_name': group_config['name'],
        'timestamp': datetime.now().isoformat(),
        'total_statements': len(group_config['statements']),
        'evaluations': {}
    }

    for statement_id, result in results.items():
        summary['evaluations'][statement_id] = {
            'evaluation': result.evaluation.value,
            'confidence': result.confidence,
            'explanation': result.explanation
        }

    # Count evaluations
    evaluation_counts = {
        'SUPPORTED': sum(1 for r in results.values() if r.evaluation == EvaluationResult.SUPPORTED),
        'REFUTED': sum(1 for r in results.values() if r.evaluation == EvaluationResult.REFUTED),
        'UNSURE': sum(1 for r in results.values() if r.evaluation == EvaluationResult.UNSURE)
    }
    summary['evaluation_counts'] = evaluation_counts

    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)

    logger.info(f"Analysis complete for {group_config['name']}")
    logger.info(f"Results: {evaluation_counts}")
    logger.info(f"Report saved to: {report_path}")
    logger.info(f"Summary saved to: {summary_path}")

    return results

def run_automated_analysis(data_path: str, automation_module: AnalysisAutomation, output_dir: str):
    """
    Run automated analysis using the analysis automation module
    """
    logger.info("Running automated analysis module")

    # Create output directory
    auto_dir = os.path.join(output_dir, 'automated_analysis')
    os.makedirs(auto_dir, exist_ok=True)

    # Run differential expression analysis
    logger.info("Performing covariate-controlled differential expression...")
    de_results = automation_module.covariate_controlled_de()
    de_results.to_csv(os.path.join(auto_dir, 'differential_expression.csv'))

    # Analyze correlations
    logger.info("Analyzing protein correlations...")
    proteins = ['TAX1BP1', 'CAT', 'VDAC1', 'CYCS', 'ATP5F1A',
                'UQCRC2', 'COX4I1', 'PRDX1', 'KEAP1', 'TFRC']

    sqstm1_corr = automation_module.calculate_protein_correlations(proteins, 'SQSTM1')
    sqstm1_corr.to_csv(os.path.join(auto_dir, 'sqstm1_correlations.csv'))

    pseudotime_corr = automation_module.calculate_protein_correlations(proteins, 'pseudotime')
    pseudotime_corr.to_csv(os.path.join(auto_dir, 'pseudotime_correlations.csv'))

    # Sliding window analysis
    logger.info("Performing sliding window correlation...")
    sliding = automation_module.sliding_window_correlation('SQSTM1', 'VDAC1')
    if not sliding.empty:
        sliding.to_csv(os.path.join(auto_dir, 'sliding_window_correlation.csv'))

    # Biphasic analysis
    logger.info("Analyzing biphasic behavior...")
    vatpase_biphasic = automation_module.analyze_biphasic_behavior('vatpase')
    proteasome_biphasic = automation_module.analyze_biphasic_behavior('proteasome')

    biphasic_results = {
        'vatpase': vatpase_biphasic,
        'proteasome': proteasome_biphasic
    }
    with open(os.path.join(auto_dir, 'biphasic_analysis.json'), 'w') as f:
        json.dump(biphasic_results, f, indent=2, default=str)

    logger.info(f"Automated analysis results saved to: {auto_dir}")

def main():
    parser = argparse.ArgumentParser(description='Run automated bioinformatics analysis')
    parser.add_argument('--data', type=str, default='pool_processed_v2.h5ad',
                       help='Path to H5AD data file')
    parser.add_argument('--group', type=int, choices=[1, 2], default=None,
                       help='Finding group to analyze (1 or 2). If not specified, both will be analyzed.')
    parser.add_argument('--output', type=str, default='analysis_output',
                       help='Output directory for results')
    parser.add_argument('--automation', action='store_true',
                       help='Run additional automated analysis')

    args = parser.parse_args()

    # Check if data file exists
    if not os.path.exists(args.data):
        logger.error(f"Data file not found: {args.data}")
        sys.exit(1)

    # Create output directory
    os.makedirs(args.output, exist_ok=True)

    # Initialize agent
    logger.info(f"Initializing agent with data: {args.data}")
    agent = BioinformaticsAgent(args.data)

    # Load data
    if not agent.load_data():
        logger.error("Failed to load data")
        sys.exit(1)

    # Run finding groups
    if args.group is None or args.group == 1:
        run_finding_group(agent, FINDING_GROUP_1, args.output)

    if args.group is None or args.group == 2:
        run_finding_group(agent, FINDING_GROUP_2, args.output)

    # Run automated analysis if requested
    if args.automation:
        automation = AnalysisAutomation(agent.adata)
        run_automated_analysis(args.data, automation, args.output)

    logger.info("="*60)
    logger.info("All analyses complete!")
    logger.info(f"Results saved to: {args.output}")
    logger.info("="*60)

if __name__ == "__main__":
    main()
```
