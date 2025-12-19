#!/usr/bin/env python3
"""
Update all claim result files with proper visualization links
and create comprehensive visualization index
"""

import os
import re

# Define visualization mappings
VISUALIZATIONS = {
    # Sequential failure claims
    'claim1_v_atpase': {
        'original': 'claim_1_plot.png',
        'paper_figures': ['figure3_volcano_histogram.png', 'figure4_vatpase_mc1.png'],
        'description': 'V-ATPase differential expression and subunit patterns'
    },
    'claim2_atp6v0a1': {
        'original': 'claim_2_plot.png',
        'paper_figures': ['figure5_vatpase_pseudotime.png', 'figure6_11_vatpase_mc1_segmented.png'],
        'description': 'ATP6V0A1 upregulation and disease stage analysis'
    },
    'claim3_organellar': {
        'original': None,
        'paper_figures': ['figure3_volcano_histogram.png'],
        'description': 'Organellar marker perturbation patterns'
    },
    'claim4_retromer': {
        'original': None,
        'paper_figures': ['figure3_volcano_histogram.png'],
        'description': 'Retromer complex dysfunction analysis'
    },
    'claim5_sos': {
        'original': None,
        'paper_figures': ['figure3_volcano_histogram.png'],
        'description': 'Stress response activation patterns'
    },
    'claim6_breakpoints': {
        'original': None,
        'paper_figures': ['figure5_vatpase_pseudotime.png', 'figure6_11_vatpase_mc1_segmented.png'],
        'description': 'Segmented progression breakpoint detection'
    },
    'claim7_temporal': {
        'original': None,
        'paper_figures': ['figure5_vatpase_pseudotime.png'],
        'description': 'Temporal ordering of system failures'
    },
    'claim8_collapse': {
        'original': None,
        'paper_figures': ['figure10_coordinated_decline.png'],
        'description': 'Network-wide proteostasis collapse'
    },
    # Mitochondrial claims
    'claim1_ups_proteins': {
        'original': 'mito_claim_1_plot.png',
        'paper_figures': ['figure7_autophagy_dysregulation.png'],
        'description': 'UPS protein differential expression'
    },
    'claim2_sqstm1': {
        'original': 'mito_claim_2_plot.png',
        'paper_figures': ['figure7_autophagy_dysregulation.png', 'figure8_sqstm1_vdac1_correlation.png'],
        'description': 'SQSTM1/p62 upregulation validation'
    },
    'claim3_correlation': {
        'original': None,
        'paper_figures': ['figure8_sqstm1_vdac1_correlation.png'],
        'description': 'BECN1-SQSTM1 correlation analysis'
    },
    'claim4_becn1': {
        'original': None,
        'paper_figures': ['figure7_autophagy_dysregulation.png'],
        'description': 'BECN1 reduction in late stages'
    },
    'claim5_mitophagy': {
        'original': None,
        'paper_figures': ['figure8_sqstm1_vdac1_correlation.png'],
        'description': 'Mitophagy pathway impairment'
    },
    'claim6_sliding': {
        'original': None,
        'paper_figures': ['figure8_sqstm1_vdac1_correlation.png'],
        'description': 'Sliding window temporal patterns'
    },
    'claim7_biphasic': {
        'original': None,
        'paper_figures': ['figure9_cycs_expression.png'],
        'description': 'Biphasic expression in VDAC1/CYCS'
    },
    'claim8_tau': {
        'original': None,
        'paper_figures': ['figure9_cycs_expression.png', 'figure10_coordinated_decline.png'],
        'description': 'Tau-mitochondrial correlations'
    }
}

def update_claim_file(filepath, claim_key):
    """Update a claim result file with enhanced visualization section"""

    with open(filepath, 'r') as f:
        content = f.read()

    viz_info = VISUALIZATIONS.get(claim_key, {})

    # Create enhanced visualization section
    viz_section = "## Visualizations\n\n"

    # Add original plot if exists
    if viz_info.get('original'):
        viz_section += f"### Analysis Plot\n"
        viz_section += f"![{viz_info['description']}](../figures/{viz_info['original']})\n\n"

    # Add paper replication figures
    if viz_info.get('paper_figures'):
        viz_section += "### Related Paper Figures\n"
        for fig in viz_info['paper_figures']:
            fig_name = fig.replace('.png', '').replace('_', ' ').title()
            viz_section += f"- [{fig_name}](../paper_replications/{fig})\n"
        viz_section += "\n"

    # Add links for Obsidian navigation
    viz_section += "### Quick Navigation\n"
    viz_section += "- [[../figures/master_analysis_dashboard|View Master Dashboard]]\n"
    viz_section += "- [[../paper_replications/summary_dashboard|View Summary Dashboard]]\n"
    viz_section += "- [[../../INDEX|Back to Main Index]]\n"

    # Replace the visualization section
    pattern = r'## Visualization.*?(?=##|\n---|\Z)'
    replacement = viz_section + "\n"

    content = re.sub(pattern, replacement, content, flags=re.DOTALL)

    # If no visualization section found, add it before Technical Notes
    if "## Visualizations" not in content:
        pattern = r'(## Technical Notes)'
        replacement = viz_section + r'\1'
        content = re.sub(pattern, replacement, content)

    return content

# Update sequential failure claims
print("Updating Sequential Failure claim files...")
seq_dir = '/Users/byron/project_plan/01_research_analysis/results/sequential_failure'

claim_files = [
    ('claim1_v_atpase_results.md', 'claim1_v_atpase'),
    ('claim2_atp6v0a1_results.md', 'claim2_atp6v0a1'),
    ('claim3_organellar_results.md', 'claim3_organellar'),
    ('claim4_retromer_results.md', 'claim4_retromer'),
    ('claim5_results.md', 'claim5_sos'),
    ('claim6_results.md', 'claim6_breakpoints'),
    ('claim7_results.md', 'claim7_temporal'),
    ('claim8_results.md', 'claim8_collapse'),
]

for filename, claim_key in claim_files:
    filepath = os.path.join(seq_dir, filename)
    if os.path.exists(filepath):
        updated_content = update_claim_file(filepath, claim_key)
        with open(filepath, 'w') as f:
            f.write(updated_content)
        print(f"  ✓ Updated {filename}")

# Update mitochondrial claims
print("\nUpdating Mitochondrial Dysregulation claim files...")
mito_dir = '/Users/byron/project_plan/01_research_analysis/results/mitochondrial_dysregulation'

mito_files = [
    ('claim1_ups_proteins_results.md', 'claim1_ups_proteins'),
    ('claim2_sqstm1_results.md', 'claim2_sqstm1'),
    ('claim3_results.md', 'claim3_correlation'),
    ('claim4_results.md', 'claim4_becn1'),
    ('claim5_results.md', 'claim5_mitophagy'),
    ('claim6_results.md', 'claim6_sliding'),
    ('claim7_results.md', 'claim7_biphasic'),
    ('claim8_results.md', 'claim8_tau'),
]

for filename, claim_key in mito_files:
    filepath = os.path.join(mito_dir, filename)
    if os.path.exists(filepath):
        updated_content = update_claim_file(filepath, claim_key)
        with open(filepath, 'w') as f:
            f.write(updated_content)
        print(f"  ✓ Updated {filename}")

# Create comprehensive visualization index
print("\nCreating visualization index...")

index_content = """# Visualization Index

## Overview
Complete index of all visualizations generated for the proteomics analysis project.

## Quick Access by Category

### 📊 Summary Dashboards
- [[figures/master_analysis_dashboard|Master Analysis Dashboard]] - 4-panel overview
- [[paper_replications/summary_dashboard|Paper Replication Dashboard]] - 9-panel comprehensive view
- [[figures/evaluation_summary|Evaluation Summary]] - Pie charts of claim outcomes

### 🔬 Sequential Failure Visualizations

#### Core Analysis Plots
- [[figures/claim_1_plot|V-ATPase Differential Expression]] - Volcano plot and heatmap
- [[figures/claim_2_plot|ATP6V0A1 Stage Analysis]] - Disease progression patterns

#### Paper Replications
- [[paper_replications/figure3_volcano_histogram|Figure 3: Proteome Remodeling]] - 2,115 DE proteins
- [[paper_replications/figure4_vatpase_mc1|Figure 4: V-ATPase Biphasic Pattern]] - MC1 analysis
- [[paper_replications/figure5_vatpase_pseudotime|Figure 5: Pseudotime Breakpoint]] - Sequential failure
- [[paper_replications/figure6_11_vatpase_mc1_segmented|Figure 6: Critical MC1 Threshold]] - MC1 = 2.831

### 🧬 Mitochondrial Dysregulation Visualizations

#### Core Analysis Plots
- [[figures/mito_claim_1_plot|UPS Protein Analysis]] - Differential expression
- [[figures/mito_claim_2_plot|SQSTM1 Upregulation]] - Validation plots

#### Paper Replications
- [[paper_replications/figure7_autophagy_dysregulation|Figure 7: Autophagy Dysregulation]] - Pathway analysis
- [[paper_replications/figure8_sqstm1_vdac1_correlation|Figure 8: SQSTM1-VDAC1 Coupling]] - Dynamic correlation
- [[paper_replications/figure9_cycs_expression|Figure 9: Cytochrome C Decline]] - Biphasic pattern
- [[paper_replications/figure10_coordinated_decline|Figure 10: Coordinated Failure]] - Mito-lyso decline

## Visualization Statistics

### File Counts
- **Core Analysis Figures**: 6
- **Paper Replications**: 9
- **Total Visualizations**: 15

### Key Findings Visualized

#### Critical Thresholds
| System | Breakpoint | Visualization |
|--------|------------|---------------|
| Proteasome | Pseudotime 0.372 | [[paper_replications/figure5_vatpase_pseudotime|Figure 5]] |
| V-ATPase (time) | Pseudotime 0.654 | [[paper_replications/figure5_vatpase_pseudotime|Figure 5]] |
| V-ATPase (MC1) | MC1 2.831 | [[paper_replications/figure6_11_vatpase_mc1_segmented|Figure 6]] |

#### Key Protein Changes
| Protein | Change | Visualization |
|---------|--------|---------------|
| SQSTM1 | +1.32 fold | [[figures/mito_claim_2_plot|SQSTM1 Plot]] |
| V-ATPase | Biphasic | [[figures/claim_1_plot|V-ATPase Plot]] |
| CYCS | -2.58 Cohen's d | [[paper_replications/figure9_cycs_expression|Figure 9]] |

## Navigation

### By Analysis Type
- **Differential Expression**: [[paper_replications/figure3_volcano_histogram|Volcano Plots]]
- **Temporal Analysis**: [[paper_replications/figure5_vatpase_pseudotime|Pseudotime Plots]]
- **Correlation Analysis**: [[paper_replications/figure8_sqstm1_vdac1_correlation|Running Correlations]]
- **Segmented Regression**: [[paper_replications/figure6_11_vatpase_mc1_segmented|Breakpoint Analysis]]

### By Biological System
- **Proteostasis**: Claims 1-8 (Sequential Failure)
- **Mitochondria**: Claims 1-8 (Mitochondrial Dysregulation)
- **Autophagy**: [[paper_replications/figure7_autophagy_dysregulation|Figure 7]]
- **Lysosomes**: [[paper_replications/figure4_vatpase_mc1|V-ATPase Figures]]

## File Locations

```
01_research_analysis/results/
├── figures/                    # Core analysis plots
│   ├── claim_1_plot.png
│   ├── claim_2_plot.png
│   ├── evaluation_summary.png
│   ├── master_analysis_dashboard.png
│   ├── mito_claim_1_plot.png
│   └── mito_claim_2_plot.png
└── paper_replications/         # Paper figure replications
    ├── figure3_volcano_histogram.png
    ├── figure4_vatpase_mc1.png
    ├── figure5_vatpase_pseudotime.png
    ├── figure6_11_vatpase_mc1_segmented.png
    ├── figure7_autophagy_dysregulation.png
    ├── figure8_sqstm1_vdac1_correlation.png
    ├── figure9_cycs_expression.png
    ├── figure10_coordinated_decline.png
    └── summary_dashboard.png
```

## Quick Links

- [[../README|Results README]]
- [[../../INDEX|Project Index]]
- [[../PAPER_FIGURE_REPLICATION_SUMMARY|Replication Summary]]
- [[../sequential_failure/README|Sequential Failure Results]]
- [[../mitochondrial_dysregulation/README|Mitochondrial Results]]

---

*Last updated: 2024-09-28*
*Total visualizations: 15 high-quality figures*
"""

index_path = '/Users/byron/project_plan/01_research_analysis/results/VISUALIZATION_INDEX.md'
with open(index_path, 'w') as f:
    f.write(index_content)

print(f"✓ Created visualization index at {index_path}")

print("\n✅ All visualization updates complete!")