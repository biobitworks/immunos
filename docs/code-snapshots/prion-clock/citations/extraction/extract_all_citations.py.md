---
source: /Users/byron/projects/prion-clock/citations/extraction/extract_all_citations.py
relative: prion-clock/citations/extraction/extract_all_citations.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
Extract ALL citations from cytoplasmic_inheritance_timer_v4.md
Outputs: full-citations.json with all 53 references
"""
import re
import json
from pathlib import Path

def extract_citations(md_path):
    """Extract all references from markdown file"""
    with open(md_path) as f:
        content = f.read()

    # Find references section
    refs_match = re.search(r'## References\n\n(.+)', content, re.DOTALL)
    if not refs_match:
        raise ValueError("No References section found")

    refs_text = refs_match.group(1)

    # Extract individual citations
    # Pattern: number. Authors. Title. Journal. Year;...
    pattern = r'(\d+)\.\s+(.+?)(?=\n\d+\.|$)'
    citations = []

    for match in re.finditer(pattern, refs_text, re.DOTALL):
        num = int(match.group(1))
        citation_text = match.group(2).strip()

        # Extract DOI/PMID if present
        doi_match = re.search(r'doi:(\S+)', citation_text)
        pmid_match = re.search(r'PMID:\s*(\d+)', citation_text)

        citations.append({
            'number': num,
            'full_text': citation_text,
            'doi': doi_match.group(1) if doi_match else None,
            'pmid': pmid_match.group(1) if pmid_match else None,
            'verified': False,
            'abstract_downloaded': False
        })

    return citations

def main():
    md_path = Path('prion-clock/zenodo-submission/cytoplasmic_inheritance_timer_v4.md')
    output_path = Path('prion-clock/citations/extraction/full-citations.json')

    citations = extract_citations(md_path)

    print(f"Extracted {len(citations)} citations")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(citations, f, indent=2)

    print(f"Saved to {output_path}")

if __name__ == '__main__':
    main()

```
