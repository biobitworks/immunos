# Research Projects

This directory contains documentation for active and completed research projects related to aging biology and longevity science.

## Purpose

Track and document:
- Research questions and hypotheses
- Project goals and milestones
- Experimental designs
- Data analysis pipelines
- Results and findings
- Publications and presentations
- Collaborations

## Project Types

### Data Analysis Projects
- Meta-analysis of aging datasets
- Gene network analysis
- Pathway enrichment studies
- Cross-species comparisons

### Computational Projects
- Machine learning models for aging prediction
- Protein interaction network analysis
- Drug target identification
- Biomarker discovery

### Integration Projects
- Database integration and harmonization
- Tool development and benchmarking
- Reproducibility studies
- Resource compilation

### Experimental Designs
- Planned wet-lab experiments
- Clinical study designs
- Collaborative research protocols

## File Naming Convention

Use lowercase with hyphens: `project-name.md`

Examples:
- `genage-cellage-overlap.md`
- `immunosenescence-markers.md`
- `ai-aging-biomarkers.md`
- `senolytic-drug-screen.md`

## Project Template

Use the standard template at `/templates/project-tracker.md` which includes:
- YAML frontmatter with project metadata
- Executive summary
- Research questions
- Background and motivation
- Methods and approach
- Timeline and milestones
- Resources needed (data, tools, collaborators)
- Results and findings
- Publications and outputs
- Links to code, data, and documentation

## Project Lifecycle

### 1. Planning Phase
- Define research question
- Literature review
- Identify datasets and tools needed
- Plan methodology

### 2. Active Development
- Data collection and preprocessing
- Analysis implementation
- Regular progress updates
- Troubleshooting and iteration

### 3. Analysis and Interpretation
- Results generation
- Statistical validation
- Biological interpretation
- Visualization

### 4. Documentation and Publication
- Write up findings
- Code documentation
- Dataset publication
- Paper preparation

### 5. Archival
- Long-term data storage
- Code repository finalization
- Knowledge transfer

## Integration with Vault

Projects should link to:
- **Datasets**: `[[datasets/dataset-name|Dataset]]`
- **Tools**: `[[tools/tool-name|Tool]]`
- **Literature**: `[[literature/paper-citation|Paper]]`
- **Researchers**: `[[researchers/lastname-firstname|Collaborator]]`
- **Organizations**: `[[organizations/org-name|Partner Org]]`

## Example Queries

### Dataview: List active projects
```dataview
TABLE status, start_date, lead_researcher
FROM "research/projects"
WHERE status = "active"
SORT start_date DESC
```

### Dataview: Find projects using specific datasets
```dataview
LIST status, research_question
FROM "research/projects"
WHERE contains(datasets, "GenAge")
```

## Current Projects

- **IMMUNOS-MCP** - Artificial immune system for code security scanning (MVP complete)
- **GenAge-CellAge Integration** - Overlap analysis between aging and senescence genes (in progress)
- **Aging Biology Knowledge Graph** - Comprehensive research infrastructure (just started)
- *(more to be added)*

---

**Directory**: `/research/projects/`
**Template**: `/templates/project-tracker.md`
**Script**: *N/A (manual creation)*
**Last Updated**: 2025-12-03
