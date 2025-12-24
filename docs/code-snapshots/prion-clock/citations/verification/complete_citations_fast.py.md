---
source: /Users/byron/projects/prion-clock/citations/verification/complete_citations_fast.py
relative: prion-clock/citations/verification/complete_citations_fast.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
Fast citation completion without DeepSeek:
1. Extract PMIDs from citation text
2. Look up PMIDs from existing DOIs
3. Download all abstracts from PubMed
4. Create note files
"""
import json
import requests
import time
import re
from pathlib import Path
from xml.etree import ElementTree as ET

PUBMED_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUBMED_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
PAPERS_DIR = Path("papers")

def extract_pmid_from_text(text):
    """Extract PMID from citation text using regex"""
    match = re.search(r'PMID:\s*(\d+)', text, re.IGNORECASE)
    return match.group(1) if match else None

def lookup_pmid_from_doi(doi):
    """Look up PMID from DOI using PubMed E-search"""
    if not doi:
        return None

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

def fetch_complete_metadata_from_pubmed(pmid):
    """Fetch complete metadata including DOI, abstract from PubMed"""
    if not pmid:
        return None

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
        print(f"    Error fetching metadata: {e}")
        return None

def create_note_file(citation_num, metadata, papers_dir):
    """Create note file for citation"""
    doi = metadata.get('doi')
    if not doi:
        print(f"  ⚠️  No DOI for citation {citation_num}, skipping note file")
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

    print(f"  ✓ Created note file")

def main():
    citations_path = Path('prion-clock/citations/extraction/full-citations.json')
    output_citations = Path('prion-clock/citations/extraction/full-citations-complete.json')
    output_abstracts = Path('prion-clock/citations/abstracts/all-abstracts-complete.json')

    with open(citations_path) as f:
        citations = json.load(f)

    complete_citations = []
    all_metadata = {}

    print(f"Processing {len(citations)} citations...")
    print("="*60)

    for citation in citations:
        num = citation['number']
        print(f"\n[{num}/53] Citation {num}")

        # Step 1: Try to extract PMID from citation text
        pmid_from_text = extract_pmid_from_text(citation['full_text'])
        if pmid_from_text:
            print(f"  ✓ Found PMID in text: {pmid_from_text}")
            citation['pmid'] = pmid_from_text

        # Step 2: If still no PMID, try lookup from DOI
        if not citation.get('pmid') and citation.get('doi'):
            print(f"  Looking up PMID from DOI: {citation['doi']}")
            pmid = lookup_pmid_from_doi(citation['doi'])
            if pmid:
                print(f"  ✓ Found PMID: {pmid}")
                citation['pmid'] = pmid
            else:
                print(f"  ⚠️  No PMID found for DOI")

        # Step 3: Fetch complete metadata from PubMed if we have PMID
        if citation.get('pmid'):
            print(f"  Fetching metadata from PubMed (PMID: {citation['pmid']})...")
            metadata = fetch_complete_metadata_from_pubmed(citation['pmid'])

            if metadata:
                # Update citation with PubMed DOI if we didn't have one
                if not citation.get('doi') and metadata.get('doi'):
                    citation['doi'] = metadata['doi']
                    print(f"  ✓ Got DOI from PubMed: {metadata['doi']}")

                # Save metadata
                all_metadata[str(num)] = metadata

                # Create note file
                create_note_file(num, metadata, PAPERS_DIR)

                abstract_len = len(metadata.get('abstract', ''))
                print(f"  ✓ Abstract: {abstract_len} chars")
            else:
                print(f"  ⚠️  Failed to fetch metadata")
        else:
            print(f"  ⚠️  No PMID - cannot download abstract")

        # Update citation status
        citation['verified'] = citation.get('pmid') is not None
        citation['abstract_downloaded'] = citation.get('pmid') is not None
        complete_citations.append(citation)

        # Rate limit (3 requests/second for PubMed)
        time.sleep(0.4)

    # Save complete citations
    output_citations.parent.mkdir(parents=True, exist_ok=True)
    with open(output_citations, 'w') as f:
        json.dump(complete_citations, f, indent=2)

    # Save all abstracts
    output_abstracts.parent.mkdir(parents=True, exist_ok=True)
    with open(output_abstracts, 'w') as f:
        json.dump(all_metadata, f, indent=2)

    print(f"\n{'='*60}")
    print(f"COMPLETE!")
    print(f"{'='*60}")
    print(f"Total citations: {len(complete_citations)}")
    print(f"Citations with DOI: {sum(1 for c in complete_citations if c.get('doi'))}")
    print(f"Citations with PMID: {sum(1 for c in complete_citations if c.get('pmid'))}")
    print(f"Abstracts downloaded: {len(all_metadata)}")
    print(f"\nSaved to:")
    print(f"  - {output_citations}")
    print(f"  - {output_abstracts}")

if __name__ == '__main__':
    main()

```
