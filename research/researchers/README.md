# Researcher Profiles

## Summary
Contains 0 subdirectories and 1 files.


This directory contains detailed profiles of researchers working in aging biology, longevity science, cellular senescence, immunosenescence, and related fields.

## Purpose

Track and document:
- Key researchers in the field
- Their research focus and contributions
- Publications and citations
- Institutional affiliations
- Collaborative networks
- Contact information and social presence

## File Naming Convention

Use lowercase with hyphens: `lastname-firstname.md`

Examples:
- `lidsky-peter.md`
- `de-magalhaes-joao-pedro.md`
- `lopez-otin-carlos.md`

## Profile Template

Use the standard template at `/templates/researcher-profile.md` which includes:
- YAML frontmatter with metadata
- Biography and background
- Research focus areas
- Key publications
- Current projects
- Institutional affiliations
- Links to external profiles (Google Scholar, ORCID, Twitter, etc.)
- Connection to our research projects

## Integration with Projects

Researcher profiles should link to:
- **Organizations**: `[[organizations/institution-name|Institution]]`
- **Publications**: `[[literature/paper-citation|Paper Title]]`
- **Projects**: `[[projects/project-name|Project Name]]`
- **Datasets**: `[[datasets/dataset-name|Dataset]]`

## Search and Discovery

Use Obsidian's graph view and Dataview plugin to:
- Find researchers by research area (using tags)
- Map collaboration networks
- Track citation patterns
- Identify emerging researchers

## Example Queries

### Dataview: List all researchers in cellular senescence
```dataview
TABLE research_focus, institution, h_index
FROM "research/researchers"
WHERE contains(research_focus, "senescence")
SORT h_index DESC
```

### Dataview: Find researchers by institution
```dataview
LIST
FROM "research/researchers"
WHERE institution = "City University of Hong Kong"
```

## Adding New Researchers

### Manual Method
1. Copy template: `cp ../templates/researcher-profile.md lastname-firstname.md`
2. Fill in YAML frontmatter and content
3. Link to relevant papers, organizations, datasets

### Automated Method
```bash
# Use automation script
bash /Users/byron/projects/scripts/create-researcher-profile.sh "FirstName LastName" "Institution" "email@domain.com"
```

## Current Researchers

- **Peter Lidsky** - Virology, RNA biology, aging connections (City University of Hong Kong)
- *(more to be added)*

---

**Directory**: `/research/researchers/`
**Template**: `/templates/researcher-profile.md`
**Script**: `/scripts/create-researcher-profile.sh`
**Last Updated**: 2025-12-03

## Directory Map
```
researchers/
└── lidsky-peter.md
```
