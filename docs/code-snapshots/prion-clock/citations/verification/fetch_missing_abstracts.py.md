---
source: /Users/byron/projects/prion-clock/citations/verification/fetch_missing_abstracts.py
relative: prion-clock/citations/verification/fetch_missing_abstracts.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
Fetch ALL missing abstracts - including PMCIDs
"""
import json
import requests
import time
from pathlib import Path
from xml.etree import ElementTree as ET

PUBMED_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUBMED_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Manual PMID additions from user
MANUAL_PMIDS = {
    48: "27304507",  # Barzilai - Metformin
    5: None,  # Williams - 1957 (pre-PubMed)
    6: None,  # Kirkwood - 1977 (check)
}

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
    """Fetch complete metadata including PMCID, DOI, abstract from PubMed"""
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

        # Extract DOI and PMCID from article IDs
        doi = None
        pmcid = None
        article_ids = root.find('.//PubmedArticle/PubmedData/ArticleIdList')
        if article_ids:
            for id_elem in article_ids.findall('ArticleId'):
                id_type = id_elem.get('IdType')
                if id_type == 'doi':
                    doi = id_elem.text
                elif id_type == 'pmc':
                    pmcid = id_elem.text

        return {
            'title': title,
            'authors': authors,
            'journal': journal,
            'year': year,
            'abstract': abstract,
            'doi': doi,
            'pmid': pmid,
            'pmcid': pmcid
        }
    except Exception as e:
        print(f"    Error fetching metadata: {e}")
        return None

def main():
    citations_path = Path('prion-clock/citations/extraction/full-citations-verified.json')
    abstracts_path = Path('prion-clock/citations/abstracts/all-abstracts-verified.json')

    with open(citations_path) as f:
        citations = json.load(f)

    # Load existing abstracts
    with open(abstracts_path) as f:
        all_metadata = json.load(f)

    # Find missing abstracts
    missing = []
    for c in citations:
        num = str(c['number'])
        if num not in all_metadata:
            missing.append(c)
        elif not all_metadata[num].get('abstract'):
            missing.append(c)

    print(f"Found {len(missing)} citations with missing abstracts")
    print("="*60)

    fetched = 0

    for citation in missing:
        num = citation['number']
        print(f"\n[{num}/53] Citation {num}")

        # Check manual PMIDs first
        if num in MANUAL_PMIDS and MANUAL_PMIDS[num]:
            citation['pmid'] = MANUAL_PMIDS[num]
            print(f"  Using manual PMID: {MANUAL_PMIDS[num]}")

        # Try to lookup PMID from DOI if we don't have one
        if not citation.get('pmid') and citation.get('doi'):
            print(f"  Looking up PMID from DOI: {citation['doi']}")
            pmid = lookup_pmid_from_doi(citation['doi'])
            if pmid:
                citation['pmid'] = pmid
                print(f"  ✓ Found PMID: {pmid}")

        # Fetch metadata if we have PMID
        if citation.get('pmid'):
            print(f"  Fetching from PubMed (PMID: {citation['pmid']})...")
            metadata = fetch_complete_metadata_from_pubmed(citation['pmid'])

            if metadata:
                all_metadata[str(num)] = metadata
                fetched += 1

                print(f"  ✓ Title: {metadata['title'][:60]}...")
                print(f"  ✓ Abstract: {len(metadata.get('abstract', ''))} chars")
                print(f"  ✓ PMCID: {metadata.get('pmcid', 'N/A')}")
            else:
                print(f"  ⚠️  Failed to fetch metadata")
        else:
            print(f"  ⚠️  No PMID - cannot fetch from PubMed")
            print(f"     (DOI: {citation.get('doi', 'N/A')})")

        # Rate limit
        time.sleep(0.4)

    # Save updated citations
    with open(citations_path, 'w') as f:
        json.dump(citations, f, indent=2)

    # Save updated abstracts
    with open(abstracts_path, 'w') as f:
        json.dump(all_metadata, f, indent=2)

    print(f"\n{'='*60}")
    print(f"COMPLETE!")
    print(f"{'='*60}")
    print(f"New abstracts fetched: {fetched}")
    print(f"Total abstracts: {len(all_metadata)}/53")
    print(f"Remaining without abstracts: {53 - len(all_metadata)}")

if __name__ == '__main__':
    main()

```
