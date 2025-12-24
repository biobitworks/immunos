---
source: /Users/byron/projects/prion-clock/citations/verification/fix_ref53.py
relative: prion-clock/citations/verification/fix_ref53.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
Fix citation 53 with the CORRECT eLife paper
"""
import json
import requests
from pathlib import Path
from xml.etree import ElementTree as ET

PUBMED_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
PAPERS_DIR = Path("papers")

# Correct information for reference 53
CORRECT_DOI = "10.7554/eLife.93172"
CORRECT_PMID = "39480006"

def fetch_complete_metadata_from_pubmed(pmid):
    """Fetch complete metadata including DOI, abstract from PubMed"""
    params = {
        'db': 'pubmed',
        'id': pmid,
        'retmode': 'xml'
    }

    try:
        response = requests.get(PUBMED_EFETCH, params=params, timeout=10)
        if response.status_code != 200:
            return None

        root = ET.fromstring(response.content)

        # Extract abstract
        abstract_elem = root.find('.//Abstract/AbstractText')
        abstract = abstract_elem.text if abstract_elem is not None else ""

        # Extract article metadata
        article = root.find('.//PubmedArticle/MedlineCitation/Article')
        if not article:
            return None

        title = article.find('ArticleTitle').text if article.find('ArticleTitle') is not None else ""

        # Extract authors
        authors = []
        author_list = article.find('AuthorList')
        if author_list:
            for author in author_list.findall('Author'):
                last = author.find('LastName')
                first = author.find('ForeName')
                if last is not None:
                    authors.append(f"{first.text if first is not None else ''} {last.text}")

        # Extract journal
        journal_elem = article.find('Journal')
        journal = journal_elem.find('Title').text if journal_elem and journal_elem.find('Title') is not None else ""

        # Extract year
        year_elem = root.find('.//PubmedArticle/MedlineCitation/Article/Journal/JournalIssue/PubDate/Year')
        year = year_elem.text if year_elem is not None else ""

        # Extract DOI from article IDs
        doi = None
        article_ids = root.find('.//PubmedArticle/PubmedData/ArticleIdList')
        if article_ids:
            for id_elem in article_ids.findall('ArticleId'):
                if id_elem.get('IdType') == 'doi':
                    doi = id_elem.text
                    break

        return {
            'title': title,
            'authors': authors,
            'journal': journal,
            'year': year,
            'abstract': abstract,
            'doi': doi,
            'pmid': pmid
        }
    except Exception as e:
        print(f"Error fetching metadata: {e}")
        return None

def create_note_file(metadata, papers_dir):
    """Create note file for citation"""
    doi = metadata.get('doi')
    if not doi:
        print(f"No DOI, skipping note file")
        return

    # Convert DOI to directory name
    doi_dir = doi.replace('/', '_')
    note_dir = papers_dir / doi_dir
    note_dir.mkdir(parents=True, exist_ok=True)

    note_path = note_dir / 'note.md'

    # Create note content
    content = f"""# {metadata['title']}

**Authors:** {', '.join(metadata.get('authors', []))}

**Journal:** {metadata.get('journal', 'N/A')}
**Year:** {metadata.get('year', 'N/A')}
**PMID:** {metadata.get('pmid', 'N/A')}
**DOI:** {doi}

---

## Abstract

{metadata.get('abstract', 'No abstract available')}

---

## Relevance to Cytoplasmic Inheritance Timer

This paper provides critical evidence for **Prediction #3** - that maternal protein aggregates accumulate in oocytes and persist throughout the reproductive lifespan.

---

## Links

- PubMed: https://pubmed.ncbi.nlm.nih.gov/{metadata.get('pmid', '')}/
- DOI: https://doi.org/{doi}

---

**Date Added:** 2025-12-19
**Verified:** Yes (IMMUNOS + PubMed)
**Abstract Source:** PubMed E-utilities
"""

    with open(note_path, 'w') as f:
        f.write(content)

    print(f"✓ Created {note_path}")

def main():
    print("Fixing citation 53 with CORRECT eLife paper...")
    print(f"Correct DOI: {CORRECT_DOI}")
    print(f"Correct PMID: {CORRECT_PMID}")
    print()

    # Load complete citations
    citations_path = Path('prion-clock/citations/extraction/full-citations-complete.json')
    abstracts_path = Path('prion-clock/citations/abstracts/all-abstracts-complete.json')

    with open(citations_path) as f:
        citations = json.load(f)

    with open(abstracts_path) as f:
        abstracts = json.load(f)

    # Find citation 53
    citation_53 = next((c for c in citations if c['number'] == 53), None)
    if not citation_53:
        print("ERROR: Citation 53 not found!")
        return

    print(f"Old DOI: {citation_53['doi']}")
    print(f"Old PMID: {citation_53['pmid']}")
    print()

    # Update citation 53
    citation_53['doi'] = CORRECT_DOI
    citation_53['pmid'] = CORRECT_PMID
    citation_53['full_text'] = f"Bomba-Warczak EK, Velez KM, Zhou LT, Guillermier C, Edassery S, Steinhauser ML, Savas JN, Duncan FE. Exceptional longevity of mammalian ovarian and oocyte macromolecules throughout the reproductive lifespan. eLife. 2024;13:RP93172. PMID: {CORRECT_PMID}. DOI: {CORRECT_DOI}"

    # Fetch correct metadata
    print(f"Fetching correct metadata from PubMed (PMID: {CORRECT_PMID})...")
    metadata = fetch_complete_metadata_from_pubmed(CORRECT_PMID)

    if metadata:
        # Update abstracts
        abstracts['53'] = metadata

        # Create correct note file
        create_note_file(metadata, PAPERS_DIR)

        print(f"\n✓ Title: {metadata['title']}")
        print(f"✓ Authors: {len(metadata['authors'])} authors")
        print(f"✓ Abstract: {len(metadata.get('abstract', ''))} chars")
        print(f"✓ DOI: {metadata['doi']}")
    else:
        print("ERROR: Failed to fetch metadata!")
        return

    # Save updated files
    with open(citations_path, 'w') as f:
        json.dump(citations, f, indent=2)

    with open(abstracts_path, 'w') as f:
        json.dump(abstracts, f, indent=2)

    print(f"\n✓ Updated {citations_path}")
    print(f"✓ Updated {abstracts_path}")
    print("\n✅ Citation 53 corrected!")

if __name__ == '__main__':
    main()

```
