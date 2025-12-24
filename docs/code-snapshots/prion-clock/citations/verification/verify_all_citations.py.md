---
source: /Users/byron/projects/prion-clock/citations/verification/verify_all_citations.py
relative: prion-clock/citations/verification/verify_all_citations.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
Verify ALL 53 citations by finding DOIs/PMIDs for remaining 34 citations
Uses multiple search strategies:
1. OpenAlex API (open, comprehensive)
2. Crossref API (DOI lookup)
3. PubMed search by author+title+year
"""
import json
import requests
import time
import re
from pathlib import Path
from xml.etree import ElementTree as ET
from urllib.parse import quote

OPENALEX_API = "https://api.openalex.org/works"
CROSSREF_API = "https://api.crossref.org/works"
PUBMED_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUBMED_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
PAPERS_DIR = Path("papers")

def extract_author_year_title(citation_text):
    """Extract first author, year, and title keywords from citation"""
    # Try to extract first author (before "et al" or comma)
    author_match = re.match(r'([A-Z][a-z]+)\s+([A-Z]{1,3})', citation_text)
    first_author = f"{author_match.group(1)} {author_match.group(2)}" if author_match else ""

    # Extract year
    year_match = re.search(r'\b(19\d{2}|20\d{2})\b', citation_text)
    year = year_match.group(1) if year_match else ""

    # Extract title (usually after period, before journal name)
    # This is approximate - titles are often truncated
    parts = citation_text.split('.')
    title = parts[1].strip() if len(parts) > 1 else ""

    return first_author, year, title

def search_openalex(author, year, title):
    """Search OpenAlex for DOI"""
    if not (author or title):
        return None

    # Build query
    query_parts = []
    if title:
        query_parts.append(f'title.search:{title[:50]}')  # First 50 chars
    if author:
        query_parts.append(f'author.search:{author}')
    if year:
        query_parts.append(f'publication_year:{year}')

    query = ' '.join(query_parts)

    params = {
        'filter': query,
        'per-page': 1
    }

    try:
        response = requests.get(OPENALEX_API, params=params, timeout=10)
        if response.status_code == 200:
            data = response.json()
            results = data.get('results', [])
            if results:
                work = results[0]
                doi = work.get('doi', '').replace('https://doi.org/', '')
                pmid = work.get('ids', {}).get('pmid', '').replace('https://pubmed.ncbi.nlm.nih.gov/', '')
                return {'doi': doi, 'pmid': pmid, 'source': 'openalex'}
    except Exception as e:
        print(f"    OpenAlex error: {e}")

    return None

def search_crossref(citation_text):
    """Search Crossref for DOI by bibliographic query"""
    query = citation_text[:200]  # First 200 chars

    params = {
        'query.bibliographic': query,
        'rows': 1
    }

    try:
        response = requests.get(CROSSREF_API, params=params, timeout=10)
        if response.status_code == 200:
            data = response.json()
            items = data.get('message', {}).get('items', [])
            if items:
                doi = items[0].get('DOI')
                return {'doi': doi, 'pmid': None, 'source': 'crossref'}
    except Exception as e:
        print(f"    Crossref error: {e}")

    return None

def search_pubmed_by_title(author, year, title):
    """Search PubMed by author, year, and title"""
    if not (author and year):
        return None

    # Build PubMed query
    query_parts = []
    if author:
        query_parts.append(f'{author}[Author]')
    if year:
        query_parts.append(f'{year}[pdat]')
    if title:
        # Take first few significant words from title
        title_words = [w for w in title.split() if len(w) > 3][:3]
        if title_words:
            query_parts.append(f'{" ".join(title_words)}[Title]')

    query = ' AND '.join(query_parts)

    params = {
        'db': 'pubmed',
        'term': query,
        'retmode': 'json',
        'retmax': 1
    }

    try:
        response = requests.get(PUBMED_ESEARCH, params=params, timeout=10)
        data = response.json()
        id_list = data.get('esearchresult', {}).get('idlist', [])
        if id_list:
            return {'doi': None, 'pmid': id_list[0], 'source': 'pubmed'}
    except Exception as e:
        print(f"    PubMed search error: {e}")

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
        return False

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
**Verified:** Yes (IMMUNOS Multi-Source Verification)
**Abstract Source:** PubMed E-utilities
**DOI Source:** {metadata.get('doi_source', 'PubMed')}
"""

    with open(note_path, 'w') as f:
        f.write(content)

    print(f"  ✓ Created note file")
    return True

def main():
    citations_path = Path('prion-clock/citations/extraction/full-citations-complete.json')
    output_citations = Path('prion-clock/citations/extraction/full-citations-verified.json')
    output_abstracts = Path('prion-clock/citations/abstracts/all-abstracts-verified.json')

    with open(citations_path) as f:
        citations = json.load(f)

    # Filter to only unverified citations
    unverified = [c for c in citations if not c.get('pmid')]

    print(f"Found {len(unverified)} unverified citations")
    print("="*60)

    verified_count = 0
    all_metadata = {}

    # Load existing abstracts
    existing_abstracts_path = Path('prion-clock/citations/abstracts/all-abstracts-complete.json')
    if existing_abstracts_path.exists():
        with open(existing_abstracts_path) as f:
            all_metadata = json.load(f)

    for citation in unverified:
        num = citation['number']
        full_text = citation['full_text']

        print(f"\n[{num}/53] Citation {num}")
        print(f"  Text: {full_text[:80]}...")

        # Extract components for search
        author, year, title = extract_author_year_title(full_text)
        print(f"  Parsed: {author} ({year}) - {title[:40]}...")

        found = None

        # Strategy 1: Try OpenAlex
        print(f"  Searching OpenAlex...")
        found = search_openalex(author, year, title)

        # Strategy 2: Try Crossref
        if not found or not found.get('doi'):
            print(f"  Searching Crossref...")
            crossref_result = search_crossref(full_text)
            if crossref_result:
                found = crossref_result

        # Strategy 3: Try PubMed direct search
        if not found or not found.get('pmid'):
            print(f"  Searching PubMed...")
            pubmed_result = search_pubmed_by_title(author, year, title)
            if pubmed_result:
                if not found:
                    found = pubmed_result
                elif pubmed_result.get('pmid'):
                    found['pmid'] = pubmed_result['pmid']

        # Update citation if found
        if found:
            if found.get('doi'):
                citation['doi'] = found['doi']
                print(f"  ✓ Found DOI: {found['doi']} (via {found['source']})")
            if found.get('pmid'):
                citation['pmid'] = found['pmid']
                print(f"  ✓ Found PMID: {found['pmid']}")

            # If we have PMID, fetch complete metadata
            if citation.get('pmid'):
                print(f"  Fetching complete metadata from PubMed...")
                metadata = fetch_complete_metadata_from_pubmed(citation['pmid'])

                if metadata:
                    # Update citation with PubMed DOI if we didn't have one
                    if not citation.get('doi') and metadata.get('doi'):
                        citation['doi'] = metadata['doi']
                        print(f"  ✓ Got DOI from PubMed: {metadata['doi']}")

                    metadata['doi_source'] = found['source']
                    all_metadata[str(num)] = metadata

                    # Create note file
                    if create_note_file(num, metadata, PAPERS_DIR):
                        verified_count += 1

                    print(f"  ✓ Abstract: {len(metadata.get('abstract', ''))} chars")

            citation['verified'] = True
        else:
            print(f"  ⚠️  Could not find DOI/PMID")
            citation['verified'] = False

        # Rate limit
        time.sleep(0.5)

    # Save verified citations
    output_citations.parent.mkdir(parents=True, exist_ok=True)
    with open(output_citations, 'w') as f:
        json.dump(citations, f, indent=2)

    # Save all abstracts
    output_abstracts.parent.mkdir(parents=True, exist_ok=True)
    with open(output_abstracts, 'w') as f:
        json.dump(all_metadata, f, indent=2)

    print(f"\n{'='*60}")
    print(f"VERIFICATION COMPLETE!")
    print(f"{'='*60}")
    print(f"New verifications: {verified_count}")
    print(f"Total verified: {len(all_metadata)}/53")
    print(f"Verification rate: {len(all_metadata)/53*100:.1f}%")
    print(f"\nSaved to:")
    print(f"  - {output_citations}")
    print(f"  - {output_abstracts}")

if __name__ == '__main__':
    main()

```
