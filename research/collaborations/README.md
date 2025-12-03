# Collaborations and Partnerships

This directory documents collaborative relationships, partnerships, and networking activities in aging biology research.

## Purpose

Track and document:
- Ongoing collaborations with researchers and labs
- Partnership agreements with organizations
- Data sharing arrangements
- Co-authorship and joint projects
- Conference connections and networking
- Community involvement

## Collaboration Types

### Research Collaborations
- Joint projects with other researchers
- Multi-institutional studies
- Cross-disciplinary partnerships
- Student mentorship

### Data Partnerships
- Data sharing agreements
- Database integration projects
- Standardization initiatives
- Open science collaborations

### Industry Partnerships
- Biotech company collaborations
- Sponsored research
- Technology transfer
- Consulting arrangements

### Community Engagement
- Open source contributions
- Conference presentations
- Workshop organization
- Online community participation

## File Naming Convention

Use lowercase with hyphens: `collaboration-description.md`

Examples:
- `cityu-lidsky-virology-aging.md`
- `huggingface-longevity-datasets.md`
- `future-house-lab-bench.md`
- `string-db-integration.md`

## Collaboration Template

Document includes:
- YAML frontmatter with metadata
- Collaboration overview
- Partners involved
- Goals and objectives
- Timeline and milestones
- Resources shared
- Outputs and deliverables
- Contact points
- Status updates

## Collaboration Stages

### 1. Initial Contact
- How connection was made
- Initial discussion topics
- Shared interests identified

### 2. Planning
- Define scope and objectives
- Agree on roles and responsibilities
- Establish timelines
- Set communication protocols

### 3. Active Collaboration
- Regular progress meetings
- Data/resource sharing
- Joint analysis work
- Problem solving

### 4. Deliverables
- Publications and papers
- Datasets released
- Software tools
- Presentations

### 5. Ongoing Relationship
- Future collaboration opportunities
- Continued data sharing
- Networking and referrals

## Integration with Vault

Collaboration documents should link to:
- **Researchers**: `[[researchers/lastname-firstname|Collaborator]]`
- **Organizations**: `[[organizations/org-name|Partner Institution]]`
- **Projects**: `[[projects/project-name|Joint Project]]`
- **Datasets**: `[[datasets/dataset-name|Shared Data]]`

## Example Queries

### Dataview: List active collaborations
```dataview
TABLE partners, status, start_date, deliverables
FROM "research/collaborations"
WHERE status = "active"
SORT start_date DESC
```

### Dataview: Find collaborations by partner
```dataview
LIST status, focus_area
FROM "research/collaborations"
WHERE contains(partners, "Peter Lidsky")
```

## Current Collaborations

- **HuggingFace longevity-db** - Open aging datasets (planned)
- **Future House LAB-Bench** - AI evaluation for biology (data usage)
- **HAGR Project** - GenAge and CellAge data integration (data user)
- *(more to be added)*

---

**Directory**: `/research/collaborations/`
**Template**: *N/A (use flexible format)*
**Script**: *N/A (manual creation)*
**Last Updated**: 2025-12-03
