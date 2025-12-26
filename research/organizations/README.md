# Organization Profiles

## Summary
Contains 0 subdirectories and 4 files.


This directory contains profiles of universities, research institutes, companies, and labs working in aging biology and longevity research.

## Purpose

Track and document:
- Academic institutions with aging research programs
- Biotech companies focused on longevity
- Research labs and centers
- Funding organizations and foundations
- Industry partnerships
- Collaborative networks

## Organization Types

### Academic Institutions
- Universities with aging biology departments
- Research institutes (e.g., Buck Institute, Max Planck)
- Medical schools and hospitals

### Companies
- Biotech startups (e.g., Calico, Unity Biotechnology)
- Pharmaceutical companies with aging programs
- Diagnostic and data companies

### Labs and Centers
- Individual research labs within institutions
- Dedicated aging research centers
- Multi-institutional collaborations

### Funding Organizations
- Government agencies (NIH, AFAR)
- Private foundations (Glenn Foundation, Longevity Fund)
- Venture capital focused on longevity

## File Naming Convention

Use lowercase with hyphens: `organization-name.md`

Examples:
- `city-university-hong-kong.md`
- `buck-institute-aging.md`
- `calico-labs.md`
- `longevity-fund.md`

## Profile Template

Use the standard template at `/templates/organization-profile.md` which includes:
- YAML frontmatter with metadata
- Overview and mission
- Key research areas
- Notable researchers and labs
- Major publications and discoveries
- Funding sources
- Partnerships and collaborations
- Contact information

## Integration with Projects

Organization profiles should link to:
- **Researchers**: `[[researchers/lastname-firstname|Researcher Name]]`
- **Publications**: `[[literature/paper-citation|Paper Title]]`
- **Projects**: `[[projects/project-name|Project Name]]`
- **Other Organizations**: `[[organizations/org-name|Partner Organization]]`

## Example Queries

### Dataview: List all biotech companies
```dataview
TABLE type, location, founded, focus_areas
FROM "research/organizations"
WHERE type = "biotech"
SORT founded DESC
```

### Dataview: Find organizations by research area
```dataview
LIST focus_areas
FROM "research/organizations"
WHERE contains(focus_areas, "cellular senescence")
```

## Current Organizations

- **City University of Hong Kong** - Academic institution (Peter Lidsky)
- **Future House** - AI for science research (LAB-Bench creators)
- **longevity-db (HuggingFace)** - Open aging biology datasets
- *(more to be added)*

---

**Directory**: `/research/organizations/`
**Template**: `/templates/organization-profile.md`
**Script**: `/scripts/create-organization-profile.sh`
**Last Updated**: 2025-12-03

## Directory Map
```
organizations/
├── city-university-hong-kong.md
├── future-house.md
├── longevity-db-huggingface.md
└── monarch-initiative.md
```
