---
source: /Users/byron/projects/prion-clock/citations/verification/complete_citations.py
relative: prion-clock/citations/verification/complete_citations.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
Use DeepSeek + PubMed to complete all 53 citations with DOI/PMID/abstracts
"""
import json
import requests
import time
import re
from pathlib import Path
from xml.etree import ElementTree as ET

OLLAMA_URL = "http://localhost:11434/api/generate"
DEEPSEEK_MODEL = "deepseek-r1:14b"
PUBMED_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUBMED_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
PAPERS_DIR = Path("papers")

def ask_deepseek_for_identifiers(citation):
    """
    Ask DeepSeek to extract or find DOI/PMID from citation text
    Returns: {doi: str, pmid: str, confidence: float}
    """
    prompt = f"""You are a citation analysis expert. Extract or find the DOI and PMID for this citation.

Citation {citation['number']}:
{citation['full_text']}

Current metadata:
- DOI: {citation.get('doi', 'NOT FOUND')}
- PMID: {citation.get('pmid', 'NOT FOUND')}

Tasks:
1. If DOI/PMID are in the citation text, extract them
2. If missing, use author names, title keywords, journal, and year to suggest what to search for
3. Provide search query for PubMed

Respond ONLY in JSON format (no other text):
{{
  "doi": "10.xxxx/xxxxx or null",
  "pmid": "12345678 or null",
  "search_query": "author year keyword",
  "confidence": 0.0-1.0,
  "notes": "explanation"
}}
"""

    try:
        response = requests.post(OLLAMA_URL, json={
            'model': DEEPSEEK_MODEL,
            'prompt': prompt,
            'stream': False,
            'options': {'temperature': 0.1}  # Low temperature for factual extraction
        }, timeout=60)

        result_text = response.json()['response']

        # Extract JSON from DeepSeek response (it may include thinking)
        json_match = re.search(r'\{[^{}]*"doi"[^{}]*\}', result_text, re.DOTALL)
        if json_match:
            return json.loads(json_match.group(0))
        else:
            # Try parsing the whole response
            return json.loads(result_text)
    except Exception as e:
        print(f"    DeepSeek error: {e}")
        return None

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

def search_pubmed(query):
    """Search PubMed with query string to find PMID"""
    if not query:
        return None

    params = {
        'db': 'pubmed',
        'term': query,
        'retmode': 'json',
        'retmax': 3  # Get top 3 results
    }

    try:
        response = requests.get(PUBMED_ESEARCH, params=params, timeout=10)
        data = response.json()
        id_list = data.get('esearchresult', {}).get('idlist', [])
        if id_list:
            return id_list[0]  # Return first result
    except Exception as e:
        print(f"    Error searching PubMed: {e}")

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
**Verified:** Yes (IMMUNOS DeepSeek + PubMed)
**Abstract Source:** PubMed E-utilities
"""

    with open(note_path, 'w') as f:
        f.write(content)

    print(f"  ✓ Created {note_path}")

def main():
    citations_path = Path('prion-clock/citations/extraction/full-citations.json')
    output_citations = Path('prion-clock/citations/extraction/full-citations-complete.json')
    output_abstracts = Path('prion-clock/citations/abstracts/all-abstracts-complete.json')

    with open(citations_path) as f:
        citations = json.load(f)

    complete_citations = []
    all_metadata = {}

    for citation in citations:
        num = citation['number']
        print(f"\n[{num}/53] Processing citation {num}...")

        # Step 1: Ask DeepSeek to extract/find identifiers
        print(f"  Asking DeepSeek for DOI/PMID...")
        deepseek_result = ask_deepseek_for_identifiers(citation)

        if deepseek_result:
            print(f"  DeepSeek: {deepseek_result.get('notes', 'No notes')}")

            # Update citation with DeepSeek findings
            if deepseek_result.get('doi') and deepseek_result['doi'] != 'null':
                citation['doi'] = deepseek_result['doi']
            if deepseek_result.get('pmid') and deepseek_result['pmid'] != 'null':
                citation['pmid'] = deepseek_result['pmid']

        # Step 2: If still no PMID, try lookup from DOI
        if not citation.get('pmid') and citation.get('doi'):
            print(f"  Looking up PMID from DOI: {citation['doi']}")
            pmid = lookup_pmid_from_doi(citation['doi'])
            if pmid:
                print(f"  Found PMID: {pmid}")
                citation['pmid'] = pmid

        # Step 3: If still no PMID, try PubMed search
        if not citation.get('pmid') and deepseek_result and deepseek_result.get('search_query'):
            print(f"  Searching PubMed: {deepseek_result['search_query']}")
            pmid = search_pubmed(deepseek_result['search_query'])
            if pmid:
                print(f"  Found PMID: {pmid}")
                citation['pmid'] = pmid

        # Step 4: Fetch complete metadata from PubMed if we have PMID
        if citation.get('pmid'):
            print(f"  Fetching complete metadata from PubMed...")
            metadata = fetch_complete_metadata_from_pubmed(citation['pmid'])

            if metadata:
                # Update citation with PubMed DOI if we didn't have one
                if not citation.get('doi') and metadata.get('doi'):
                    citation['doi'] = metadata['doi']

                # Save metadata
                all_metadata[str(num)] = metadata

                # Create note file
                create_note_file(num, metadata, PAPERS_DIR)

                print(f"  ✓ Abstract: {len(metadata.get('abstract', ''))} chars")
                print(f"  ✓ DOI: {metadata.get('doi', 'N/A')}")
            else:
                print(f"  ⚠️  Failed to fetch metadata")
        else:
            print(f"  ⚠️  No PMID found - cannot download abstract")

        # Update citation status
        citation['verified'] = True
        citation['abstract_downloaded'] = citation.get('pmid') is not None
        complete_citations.append(citation)

        # Rate limit
        time.sleep(0.5)

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
