---
source: /Users/byron/projects/prion-clock/citations/create_v5.py
relative: prion-clock/citations/create_v5.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
Create v5.0 of cytoplasmic inheritance timer with:
1. Corrected reference 53 (eLife paper instead of Dev Cell)
2. Abstracts added to all 19 citations where available
"""
import json
import re
from pathlib import Path

# Correct info for ref 53
REF_53_CORRECT = """53. Bomba-Warczak EK, Velez KM, Zhou LT, Guillermier C, Edassery S, Steinhauser ML, Savas JN, Duncan FE.
    Exceptional longevity of mammalian ovarian and oocyte macromolecules throughout the reproductive lifespan.
    eLife. 2024;13:RP93172.
    PMID: 39480006. PMCID: PMC11527430. DOI: 10.7554/eLife.93172

    Abstract: The mechanisms contributing to age-related deterioration of the female reproductive system
    are complex, however aberrant protein homeostasis is a major contributor. We elucidated exceptionally
    stable proteins, structures, and macromolecules that persist in mammalian ovaries and gametes across
    the reproductive lifespan. Ovaries exhibit localized structural and cell-type-specific enrichment of
    stable macromolecules in both the follicular and extrafollicular environments. Moreover, ovaries and
    oocytes both harbor a panel of exceptionally long-lived proteins, including cytoskeletal, mitochondrial,
    and oocyte-derived proteins. The exceptional persistence of these long-lived molecules suggest a critical
    role in lifelong maintenance and age-dependent deterioration of reproductive tissues.

    Available at: https://pubmed.ncbi.nlm.nih.gov/39480006/ | https://doi.org/10.7554/eLife.93172"""

def main():
    # Load v4.0
    v4_path = Path('prion-clock/zenodo-submission/cytoplasmic_inheritance_timer_v4.md')
    with open(v4_path) as f:
        v4_content = f.read()

    # Load abstracts (use complete version)
    abstracts_path = Path('prion-clock/citations/abstracts/all-abstracts-verified.json')
    with open(abstracts_path) as f:
        abstracts = json.load(f)

    # Load citations (use complete version)
    citations_path = Path('prion-clock/citations/extraction/full-citations-verified.json')
    with open(citations_path) as f:
        citations = json.load(f)

    # Split into sections
    refs_start = v4_content.find('## References\n\n')
    before_refs = v4_content[:refs_start + len('## References\n\n')]
    refs_section = v4_content[refs_start + len('## References\n\n'):]

    # Build new references section
    new_refs = []

    for citation in citations:
        num = citation['number']
        full_text = citation['full_text']

        # Special handling for ref 53 - use corrected version
        if num == 53:
            new_refs.append(REF_53_CORRECT)
            continue

        # For other citations, add abstract if available
        if str(num) in abstracts:
            metadata = abstracts[str(num)]
            abstract = metadata.get('abstract', '')
            # Use PMID and DOI from metadata (more reliable than citation object)
            pmid = metadata.get('pmid', citation.get('pmid', 'N/A'))
            doi = metadata.get('doi', citation.get('doi', 'N/A'))

            # Add abstract and links
            enhanced_ref = f"{num}. {full_text}\n"

            if abstract:
                enhanced_ref += f"\n    Abstract: {abstract}\n"

            if pmid != 'N/A' and doi != 'N/A':
                enhanced_ref += f"\n    Available at: https://pubmed.ncbi.nlm.nih.gov/{pmid}/ | https://doi.org/{doi}"
            elif pmid != 'N/A':
                enhanced_ref += f"\n    Available at: https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            elif doi != 'N/A':
                enhanced_ref += f"\n    Available at: https://doi.org/{doi}"

            new_refs.append(enhanced_ref)
        else:
            # No abstract - just use original
            new_refs.append(f"{num}. {full_text}")

    # Combine
    v5_content = before_refs + '\n\n'.join(new_refs) + '\n'

    # Write v5.0
    v5_path = Path('prion-clock/zenodo-submission/cytoplasmic_inheritance_timer_v5.md')
    with open(v5_path, 'w') as f:
        f.write(v5_content)

    print(f"✓ Created {v5_path}")
    print(f"  - Fixed reference 53 (eLife paper)")
    print(f"  - Added abstracts for {len(abstracts)} citations")
    print(f"  - Total references: {len(citations)}")

if __name__ == '__main__':
    main()

```
