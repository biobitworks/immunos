---
source: /Users/byron/projects/research/projects/bioviztech/lab-resources/gapmap_search.py
relative: research/projects/bioviztech/lab-resources/gapmap_search.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
GAP-MAP Database Search Tool for BioViz Tech

Searches the Convergent Research GAP-MAP database for:
- R&D gaps related to imaging, aging, ECM, biophysics
- Foundational capabilities addressing these gaps
- Resources and publications

Citation: Convergent Research (2025). Gap Map Database.
Retrieved from https://gap-map.org/

Author: Byron P. Lee
Project: BioViz Tech - AI-Enhanced Picosecond Ultrasonic Imaging
Created: 2025-09-30
"""

import json
import sys
from pathlib import Path
from typing import List, Dict, Any

# BioViz Tech research keywords
BIOVIZ_KEYWORDS = {
    'imaging': ['imaging', 'microscopy', 'ultrasonics', 'acoustics', 'phonon',
                'nanoscale', 'resolution', 'non-destructive', 'label-free'],
    'aging': ['aging', 'longevity', 'senescence', 'lifespan', 'regeneration'],
    'ecm': ['extracellular', 'matrix', 'ecm', 'tissue', 'mechanobiology',
            'mechanical properties'],
    'cell_bio': ['cell', 'cellular', 'live cell', 'organoid', 'tissue'],
    'biophysics': ['biophysics', 'biomechanics', 'physical', 'mechanical']
}

def load_json(filepath: Path) -> Any:
    """Load JSON data from file."""
    with open(filepath, 'r') as f:
        return json.load(f)

def search_text(text: str, keywords: List[str]) -> List[str]:
    """Search text for keywords, return matches."""
    text_lower = text.lower()
    return [kw for kw in keywords if kw in text_lower]

def search_gaps(gaps: List[Dict], keywords: Dict[str, List[str]]) -> List[Dict]:
    """Search gaps for BioViz Tech relevant topics."""
    relevant = []

    for gap in gaps:
        name = gap.get('name', '')
        description = gap.get('description', '')
        tags = ' '.join(gap.get('tags', []))
        combined = f"{name} {description} {tags}"

        matches = {}
        for category, kws in keywords.items():
            found = search_text(combined, kws)
            if found:
                matches[category] = found

        if matches:
            relevant.append({
                'id': gap['id'],
                'name': gap['name'],
                'field': gap.get('field', {}).get('name', 'N/A') if gap.get('field') else 'N/A',
                'description': gap.get('description', ''),
                'tags': gap.get('tags', []),
                'capabilities_count': len(gap.get('foundationalCapabilities', [])),
                'matches': matches
            })

    return relevant

def search_capabilities(capabilities: List[Dict], keywords: Dict[str, List[str]]) -> List[Dict]:
    """Search capabilities for BioViz Tech relevant topics."""
    relevant = []

    for cap in capabilities:
        if 'description' not in cap:
            continue

        name = cap.get('name', '')
        description = cap.get('description', '')
        combined = f"{name} {description}"

        matches = {}
        for category, kws in keywords.items():
            found = search_text(combined, kws)
            if found:
                matches[category] = found

        if matches:
            relevant.append({
                'id': cap['id'],
                'name': cap['name'],
                'description': cap.get('description', ''),
                'gaps_count': len(cap.get('gaps', [])),
                'resources_count': len(cap.get('resources', [])),
                'matches': matches
            })

    return relevant

def print_gap_report(gaps: List[Dict], show_limit: int = 20):
    """Print formatted gap analysis."""
    print(f"\n{'='*80}")
    print(f"R&D GAPS RELEVANT TO BIOVIZ TECH: {len(gaps)} found")
    print(f"{'='*80}\n")

    for i, gap in enumerate(gaps[:show_limit], 1):
        print(f"{i}. {gap['name']}")
        print(f"   Field: {gap['field']}")
        print(f"   Capabilities addressing this: {gap['capabilities_count']}")

        # Show matched categories
        categories = list(gap['matches'].keys())
        print(f"   Relevant to: {', '.join(categories)}")

        # Show description (truncated)
        desc = gap['description']
        if len(desc) > 250:
            desc = desc[:250] + '...'
        print(f"   {desc}")
        print()

def print_capability_report(capabilities: List[Dict], show_limit: int = 20):
    """Print formatted capability analysis."""
    print(f"\n{'='*80}")
    print(f"CAPABILITIES RELEVANT TO BIOVIZ TECH: {len(capabilities)} found")
    print(f"{'='*80}\n")

    for i, cap in enumerate(capabilities[:show_limit], 1):
        print(f"{i}. {cap['name']}")
        print(f"   Addresses {cap['gaps_count']} gaps | {cap['resources_count']} resources")

        # Show matched categories
        categories = list(cap['matches'].keys())
        print(f"   Relevant to: {', '.join(categories)}")

        # Show description (truncated)
        desc = cap['description']
        if len(desc) > 250:
            desc = desc[:250] + '...'
        print(f"   {desc}")
        print()

def find_gap_by_name(gaps: List[Dict], search_term: str) -> List[Dict]:
    """Find specific gap by name search."""
    matches = []
    search_lower = search_term.lower()
    for gap in gaps:
        if search_lower in gap['name'].lower():
            matches.append(gap)
    return matches

def main():
    # Check if data directory exists
    data_dir = Path('/Users/byron/bioviztech/03_Lab_Resources/GAPMAP_Data')
    if not data_dir.exists():
        print(f"Error: GAP-MAP data not found at {data_dir}")
        print("Please extract gapmap-data.zip first")
        sys.exit(1)

    # Load data
    print("Loading GAP-MAP database...")
    gaps_file = data_dir / 'gaps.json'
    capabilities_file = data_dir / 'capabilities.json'
    metadata_file = data_dir / 'metadata.json'

    gaps = load_json(gaps_file)
    capabilities = load_json(capabilities_file)
    metadata = load_json(metadata_file)

    print(f"Database version: {metadata['version']}")
    print(f"Export date: {metadata['exportDate']}")
    print(f"Total gaps: {metadata['counts']['gaps']}")
    print(f"Total capabilities: {metadata['counts']['capabilities']}")
    print(f"Total resources: {metadata['counts']['resources']}")

    # Search for relevant gaps
    print("\nSearching for BioViz Tech relevant gaps...")
    relevant_gaps = search_gaps(gaps, BIOVIZ_KEYWORDS)

    # Search for relevant capabilities
    print("Searching for BioViz Tech relevant capabilities...")
    relevant_capabilities = search_capabilities(capabilities, BIOVIZ_KEYWORDS)

    # Print reports
    print_gap_report(relevant_gaps, show_limit=15)
    print_capability_report(relevant_capabilities, show_limit=15)

    # Highlight the key gap
    print(f"\n{'='*80}")
    print("KEY GAP FOR BIOVIZ TECH")
    print(f"{'='*80}\n")

    live_cell_gaps = find_gap_by_name(relevant_gaps, "Live Cell Imaging")
    if live_cell_gaps:
        gap = live_cell_gaps[0]
        print(f"GAP: {gap['name']}")
        print(f"Field: {gap['field']}")
        print(f"\nDescription:\n{gap['description']}")
        print(f"\nThis gap has {gap['capabilities_count']} proposed capabilities")
        print("\n→ BioViz Tech directly addresses this gap with picosecond ultrasonics:")
        print("  • Non-destructive (no photobleaching)")
        print("  • Label-free (no fluorophores)")
        print("  • Nanoscale resolution (<100nm)")
        print("  • Longitudinal tracking (24-48 hours)")
        print("  • Mechanical property mapping (ECM stiffness)")

    # Summary
    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")
    print(f"Found {len(relevant_gaps)} relevant R&D gaps")
    print(f"Found {len(relevant_capabilities)} relevant foundational capabilities")
    print(f"\nCitation: Convergent Research (2025). Gap Map Database.")
    print(f"Retrieved from https://gap-map.org/")

if __name__ == '__main__':
    main()

```
