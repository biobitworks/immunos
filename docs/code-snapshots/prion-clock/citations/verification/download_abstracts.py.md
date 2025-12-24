---
source: /Users/byron/projects/prion-clock/citations/verification/download_abstracts.py
relative: prion-clock/citations/verification/download_abstracts.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
Download abstracts from PubMed for all 53 citations
Create note files in /papers/[DOI]/note.md
"""
import json
import requests
import time
from pathlib import Path
from xml.etree import ElementTree as ET

PUBMED_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
PUBMED_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PAPERS_DIR = Path("papers")

def lookup_pmid_from_doi(doi):
    """Look up PMID from DOI using PubMed E-search"""
    params = {
        'db': 'pubmed',
        'term': doi,
        'retmode': 'json'
    }

    try:
        response = requests.get(PUBMED_ESEARCH, params=params, timeout=10)
        data = response.json()
        id_list = data.get('esearchresult', {}).get('idlist', [])
        if id_list:
            return id_list[0]
    except Exception as e:
        print(f"    Error looking up PMID: {e}")

    return None

def fetch_abstract_from_pubmed(pmid):
    """Fetch abstract from PubMed using E-utilities"""
    params = {
        'db': 'pubmed',
        'id': pmid,
        'retmode': 'xml'
    }

    response = requests.get(PUBMED_EFETCH, params=params)
    if response.status_code != 200:
        return None

    root = ET.fromstring(response.content)

    # Extract abstract
    abstract_elem = root.find('.//Abstract/AbstractText')
    abstract = abstract_elem.text if abstract_elem is not None else ""

    # Extract metadata
    article = root.find('.//PubmedArticle/MedlineCitation/Article')

    title = article.find('ArticleTitle').text if article else ""

    authors = []
    author_list = article.find('AuthorList') if article else None
    if author_list:
        for author in author_list.findall('Author'):
            last = author.find('LastName')
            first = author.find('ForeName')
            if last is not None:
                authors.append(f"{first.text if first is not None else ''} {last.text}")

    journal_elem = article.find('Journal')
    journal = journal_elem.find('Title').text if journal_elem else ""

    year_elem = root.find('.//PubmedArticle/MedlineCitation/Article/Journal/JournalIssue/PubDate/Year')
    year = year_elem.text if year_elem is not None else ""

    return {
        'title': title,
        'authors': authors,
        'journal': journal,
        'year': year,
        'abstract': abstract
    }

def create_note_file(citation, metadata, papers_dir):
    """Create note file for citation"""
    doi = citation.get('doi')
    if not doi:
        print(f"  ⚠️  No DOI for citation {citation['number']}, skipping note file")
        return

    # Convert DOI to directory name (replace / with _)
    doi_dir = doi.replace('/', '_')
    note_dir = papers_dir / doi_dir
    note_dir.mkdir(parents=True, exist_ok=True)

    note_path = note_dir / 'note.md'

    # Create note content
    content = f"""# {metadata['title']}

**Authors:** {', '.join(metadata['authors'])}

**Journal:** {metadata['journal']}
**Year:** {metadata['year']}
**PMID:** {citation.get('pmid', 'N/A')}
**DOI:** {doi}

---

## Abstract

{metadata['abstract']}

---

## Links

- PubMed: https://pubmed.ncbi.nlm.nih.gov/{citation.get('pmid', '')}/
- DOI: https://doi.org/{doi}

---

**Date Added:** 2025-12-19
**Verified:** Yes (IMMUNOS DeepSeek)
**Abstract Source:** PubMed E-utilities
"""

    with open(note_path, 'w') as f:
        f.write(content)

    print(f"  ✓ Created {note_path}")

def main():
    citations_path = Path('prion-clock/citations/extraction/full-citations.json')
    output_path = Path('prion-clock/citations/abstracts/all-abstracts.json')

    with open(citations_path) as f:
        citations = json.load(f)

    all_metadata = {}

    for citation in citations:
        print(f"Downloading abstract for citation {citation['number']}...")

        pmid = citation.get('pmid')

        # If no PMID, try to look it up from DOI
        if not pmid:
            doi = citation.get('doi')
            if doi:
                print(f"  Looking up PMID from DOI: {doi}")
                pmid = lookup_pmid_from_doi(doi)
                if pmid:
                    print(f"  Found PMID: {pmid}")
                    citation['pmid'] = pmid  # Update citation
                else:
                    print(f"  ⚠️  Could not find PMID for DOI, skipping")
                    continue
            else:
                print(f"  ⚠️  No PMID or DOI, skipping")
                continue

        # Fetch from PubMed
        metadata = fetch_abstract_from_pubmed(pmid)

        if metadata:
            all_metadata[citation['number']] = metadata

            # Create note file
            create_note_file(citation, metadata, PAPERS_DIR)

            # Rate limit: PubMed allows 3 requests/second
            time.sleep(0.4)
        else:
            print(f"  ⚠️  Failed to fetch abstract")

    # Save all abstracts
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(all_metadata, f, indent=2)

    print(f"\nDownloaded {len(all_metadata)} abstracts")
    print(f"Saved to {output_path}")

if __name__ == '__main__':
    main()

```
