---
source: /Users/byron/projects/scripts/immunos_cite_verify.py
relative: scripts/immunos_cite_verify.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
IMMUNOS Citation Verification System
Automated citation verification using CrossRef, PubMed, and DOI.org APIs

Usage:
    python immunos_cite_verify.py verify --project prion-clock --priority high
    python immunos_cite_verify.py download-abstracts --project prion-clock
    python immunos_cite_verify.py generate-bibtex --project prion-clock
    python immunos_cite_verify.py status --project prion-clock

IMMUNOS Integration:
    - Memory caching (30-day TTL)
    - Snapshot system for verification sessions
    - Journal entries for daily reports
"""

import argparse
import json
import os
import sys
import time
from dataclasses import dataclass, asdict
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import hashlib
import subprocess

try:
    import requests
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry
except ImportError:
    print("ERROR: requests library required. Install with: pip3 install requests")
    sys.exit(1)


# IMMUNOS Integration
IMMUNOS_MEMORY_DIR = Path.home() / "projects" / ".immunos" / "memory" / "citations"
IMMUNOS_JOURNAL_DIR = Path.home() / "projects" / ".immunos" / "journal"
CACHE_TTL_DAYS = 30


@dataclass
class VerificationResult:
    """Citation verification result"""
    status: str  # "verified", "uncertain", "failed"
    confidence: float  # 0.0-1.0
    metadata: Dict
    issues: List[str]
    sources: List[str]  # APIs used
    verified_date: str
    cache_hit: bool = False


@dataclass
class WebVerificationResult:
    """Web source verification result"""
    status: str
    html_snapshot: Optional[str]
    pdf_snapshot: Optional[str]
    wayback_url: Optional[str]
    credibility_score: float
    issues: List[str]


class ImmunosMemory:
    """T-cell adaptive memory for citation verification results"""

    def __init__(self):
        self.memory_dir = IMMUNOS_MEMORY_DIR
        self.memory_dir.mkdir(parents=True, exist_ok=True)

    def _get_cache_key(self, citation: Dict) -> str:
        """Generate cache key from citation identifiers"""
        if citation.get('doi'):
            return f"doi_{citation['doi'].replace('/', '_')}"
        elif citation.get('pmid'):
            return f"pmid_{citation['pmid']}"
        elif citation.get('url'):
            url_hash = hashlib.md5(citation['url'].encode()).hexdigest()[:12]
            return f"url_{url_hash}"
        else:
            # Fallback: hash of title
            title_hash = hashlib.md5(citation.get('title', '').encode()).hexdigest()[:12]
            return f"title_{title_hash}"

    def get(self, citation: Dict) -> Optional[VerificationResult]:
        """Retrieve cached verification result"""
        cache_key = self._get_cache_key(citation)
        cache_file = self.memory_dir / f"{cache_key}.json"

        if not cache_file.exists():
            return None

        try:
            with open(cache_file, 'r') as f:
                data = json.load(f)

            # Check TTL
            cached_date = datetime.fromisoformat(data['verified_date'])
            if datetime.now() - cached_date > timedelta(days=CACHE_TTL_DAYS):
                cache_file.unlink()  # Expired
                return None

            result = VerificationResult(**data)
            result.cache_hit = True
            return result

        except Exception as e:
            print(f"Warning: Cache read error for {cache_key}: {e}")
            return None

    def store(self, citation: Dict, result: VerificationResult):
        """Store verification result in cache"""
        cache_key = self._get_cache_key(citation)
        cache_file = self.memory_dir / f"{cache_key}.json"

        try:
            with open(cache_file, 'w') as f:
                json.dump(asdict(result), f, indent=2)
        except Exception as e:
            print(f"Warning: Cache write error for {cache_key}: {e}")


class CitationVerifier:
    """Automated citation verification system"""

    def __init__(self, project_name: str):
        self.project_name = project_name
        self.project_dir = Path.home() / "projects" / project_name
        self.citations_dir = self.project_dir / "citations"
        self.memory = ImmunosMemory()

        # Setup session with retry logic
        self.session = requests.Session()
        retry_strategy = Retry(
            total=3,
            backoff_factor=1,
            status_forcelist=[429, 500, 502, 503, 504]
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        self.session.mount("http://", adapter)
        self.session.mount("https://", adapter)

        # Set user agent for polite API usage
        self.session.headers.update({
            'User-Agent': 'IMMUNOS-CitationVerifier/1.0 (byron@immunos.dev)'
        })

    def classify_citation(self, citation: Dict) -> str:
        """Classify citation risk level"""
        risk_level = citation.get('risk_level', '').upper()
        if risk_level in ['HIGH', 'MEDIUM', 'LOW']:
            return risk_level

        # Auto-classify if not provided
        if citation.get('type') == 'web':
            return 'HIGH'

        year = citation.get('year', 9999)
        if year < 2000:
            return 'HIGH'
        elif year < 2010:
            return 'MEDIUM'
        else:
            return 'LOW'

    def verify_via_doi(self, doi: str) -> Tuple[bool, Dict]:
        """Verify citation via DOI.org and CrossRef"""
        try:
            # Try DOI.org content negotiation for metadata
            doi_url = f"https://doi.org/{doi}"
            headers = {'Accept': 'application/vnd.citationstyles.csl+json'}

            response = self.session.get(doi_url, headers=headers, timeout=10)

            if response.status_code == 200:
                metadata = response.json()
                return True, metadata

            # Fallback to CrossRef API
            crossref_url = f"https://api.crossref.org/works/{doi}"
            response = self.session.get(crossref_url, timeout=10)

            if response.status_code == 200:
                data = response.json()
                metadata = data.get('message', {})
                return True, metadata

            return False, {}

        except Exception as e:
            print(f"  DOI verification error: {e}")
            return False, {}

    def verify_via_pubmed(self, pmid: str) -> Tuple[bool, Dict]:
        """Verify citation via PubMed E-utilities"""
        try:
            # PubMed ESummary API
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            params = {
                'db': 'pubmed',
                'id': pmid,
                'retmode': 'json'
            }

            response = self.session.get(url, params=params, timeout=10)

            if response.status_code == 200:
                data = response.json()
                result = data.get('result', {}).get(pmid, {})

                if result and 'error' not in result:
                    return True, result

            return False, {}

        except Exception as e:
            print(f"  PubMed verification error: {e}")
            return False, {}

    def download_abstract_pubmed(self, pmid: str) -> Optional[str]:
        """Download abstract from PubMed"""
        try:
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            params = {
                'db': 'pubmed',
                'id': pmid,
                'rettype': 'abstract',
                'retmode': 'text'
            }

            response = self.session.get(url, params=params, timeout=10)

            if response.status_code == 200:
                abstract_text = response.text.strip()
                if abstract_text and len(abstract_text) > 50:
                    return abstract_text

            return None

        except Exception as e:
            print(f"  PubMed abstract download error: {e}")
            return None

    def verify_citation(self, citation: Dict) -> VerificationResult:
        """Multi-source citation verification"""

        # Check cache first
        cached_result = self.memory.get(citation)
        if cached_result:
            return cached_result

        citation_key = citation.get('citation_key', 'unknown')
        print(f"\nVerifying: {citation_key}")

        issues = []
        sources = []
        metadata = {}
        confidence = 0.0

        # Strategy 1: DOI verification
        if citation.get('doi'):
            print(f"  Checking DOI: {citation['doi']}")
            success, doi_metadata = self.verify_via_doi(citation['doi'])

            if success:
                sources.append('DOI.org/CrossRef')
                metadata.update(doi_metadata)
                confidence += 0.5
            else:
                issues.append(f"DOI {citation['doi']} not found")

        # Strategy 2: PubMed verification
        if citation.get('pmid'):
            print(f"  Checking PMID: {citation['pmid']}")
            time.sleep(0.34)  # Rate limit: 3 requests/sec

            success, pubmed_metadata = self.verify_via_pubmed(citation['pmid'])

            if success:
                sources.append('PubMed')
                metadata['pubmed'] = pubmed_metadata
                confidence += 0.4
            else:
                issues.append(f"PMID {citation['pmid']} not found")

        # Determine status
        if confidence >= 0.8:
            status = "verified"
        elif confidence >= 0.4:
            status = "uncertain"
        else:
            status = "failed"
            if not issues:
                issues.append("Insufficient verification sources")

        result = VerificationResult(
            status=status,
            confidence=confidence,
            metadata=metadata,
            issues=issues,
            sources=sources,
            verified_date=datetime.now().isoformat()
        )

        # Cache result
        self.memory.store(citation, result)

        print(f"  Status: {status.upper()} (confidence: {confidence:.2f})")

        return result

    def download_abstract(self, citation: Dict) -> Optional[str]:
        """Download abstract from available sources"""

        citation_key = citation.get('citation_key', 'unknown')

        # Try PubMed first
        if citation.get('pmid'):
            print(f"  Downloading abstract from PubMed...")
            time.sleep(0.34)  # Rate limit

            abstract = self.download_abstract_pubmed(citation['pmid'])
            if abstract:
                # Save to file
                abstract_file = self.citations_dir / "abstracts" / f"{citation_key}.txt"
                abstract_file.parent.mkdir(parents=True, exist_ok=True)

                with open(abstract_file, 'w') as f:
                    f.write(f"Citation: {citation.get('text', citation_key)}\n")
                    f.write(f"PMID: {citation['pmid']}\n")
                    if citation.get('doi'):
                        f.write(f"DOI: {citation['doi']}\n")
                    f.write(f"\n{abstract}\n")

                print(f"  ✓ Abstract saved to: abstracts/{citation_key}.txt")
                return abstract

        # TODO: Add CrossRef abstract extraction
        # TODO: Add publisher API support

        return None

    def verify_web_source(self, citation: Dict) -> WebVerificationResult:
        """Verify and archive web source"""

        url = citation.get('url', '')
        citation_key = citation.get('citation_key', 'unknown')

        print(f"\nVerifying web source: {citation_key}")
        print(f"  URL: {url}")

        issues = []
        credibility_score = 0.5  # Default

        # Check domain credibility
        if '.gov' in url:
            credibility_score = 0.95
        elif '.edu' in url:
            credibility_score = 0.90
        elif any(org in url for org in ['nih.', 'cdc.', 'who.int']):
            credibility_score = 0.95
        elif any(org in url for org in ['.ac.uk', '.org']):
            credibility_score = 0.70
        else:
            credibility_score = 0.50
            issues.append("Domain credibility uncertain")

        # Try to fetch webpage
        html_snapshot = None
        try:
            response = self.session.get(url, timeout=15)
            if response.status_code == 200:
                # Save HTML snapshot
                web_dir = self.citations_dir / "verification" / "web-snapshots"
                web_dir.mkdir(parents=True, exist_ok=True)

                html_file = web_dir / f"{citation_key}.html"
                with open(html_file, 'w', encoding='utf-8') as f:
                    f.write(response.text)

                html_snapshot = str(html_file)
                print(f"  ✓ HTML snapshot saved")
            else:
                issues.append(f"HTTP {response.status_code}")

        except Exception as e:
            issues.append(f"Fetch error: {e}")

        # Check Wayback Machine
        wayback_url = None
        try:
            wayback_api = f"http://archive.org/wayback/available?url={url}"
            response = self.session.get(wayback_api, timeout=10)

            if response.status_code == 200:
                data = response.json()
                archived = data.get('archived_snapshots', {}).get('closest', {})
                if archived:
                    wayback_url = archived.get('url')
                    print(f"  ✓ Wayback Machine: {wayback_url}")

        except Exception as e:
            print(f"  Wayback check failed: {e}")

        status = "verified" if credibility_score >= 0.7 else "uncertain"

        return WebVerificationResult(
            status=status,
            html_snapshot=html_snapshot,
            pdf_snapshot=None,
            wayback_url=wayback_url,
            credibility_score=credibility_score,
            issues=issues
        )

    def generate_bibtex_entry(self, citation: Dict, verified_metadata: Dict) -> str:
        """Generate BibTeX entry from verified metadata"""

        citation_key = citation.get('citation_key', 'unknown')
        entry_type = citation.get('type', 'academic')

        if entry_type == 'web':
            # @misc for web sources
            authors = citation.get('authors', ['Unknown'])
            author_str = '{' + ' '.join(authors) + '}'

            bibtex = f"@misc{{{citation_key},\n"
            bibtex += f"  title = {{{citation.get('title', 'Untitled')}}},\n"
            bibtex += f"  author = {author_str},\n"
            bibtex += f"  year = {{{citation.get('year', 'n.d.')}}},\n"
            bibtex += f"  url = {{{citation.get('url', '')}}},\n"
            bibtex += f"  note = {{Accessed {datetime.now().strftime('%Y-%m-%d')}}}\n"
            bibtex += "}\n"

        else:
            # @article for academic papers
            authors = citation.get('authors', [])
            if len(authors) >= 2:
                author_str = f"{authors[0]}, {authors[1]}"
                if len(authors) > 2:
                    author_str += " and others"
            else:
                author_str = authors[0] if authors else "Unknown"

            bibtex = f"@article{{{citation_key},\n"
            bibtex += f"  author = {{{author_str}}},\n"
            bibtex += f"  title = {{{citation.get('title', 'Untitled')}}},\n"
            bibtex += f"  journal = {{{citation.get('journal', 'Unknown')}}},\n"
            bibtex += f"  year = {{{citation.get('year', 'n.d.')}}},\n"

            if citation.get('doi'):
                bibtex += f"  doi = {{{citation['doi']}}},\n"
            if citation.get('pmid'):
                bibtex += f"  pmid = {{{citation['pmid']}}},\n"
            if citation.get('pmcid'):
                bibtex += f"  pmcid = {{{citation['pmcid']}}},\n"

            bibtex += f"  verified = {{{datetime.now().strftime('%Y-%m-%d')}}}\n"
            bibtex += "}\n"

        return bibtex

    def verify_all(self, priority: Optional[str] = None) -> Dict:
        """Verify all citations with optional priority filter"""

        # Load citations
        citations_file = self.citations_dir / "extraction" / "inline-citations.json"

        if not citations_file.exists():
            print(f"ERROR: Citations file not found: {citations_file}")
            return {}

        with open(citations_file, 'r') as f:
            data = json.load(f)

        citations = data.get('citations', [])

        # Filter by priority
        if priority:
            priority_upper = priority.upper()
            citations = [c for c in citations if self.classify_citation(c) == priority_upper]
            print(f"Filtering to {priority_upper} priority: {len(citations)} citations")

        results = {
            'verified': [],
            'uncertain': [],
            'failed': [],
            'total': len(citations)
        }

        for citation in citations:
            citation_type = citation.get('type', 'academic')

            if citation_type == 'web':
                web_result = self.verify_web_source(citation)

                result_dict = {
                    'citation': citation,
                    'status': web_result.status,
                    'credibility': web_result.credibility_score,
                    'issues': web_result.issues,
                    'wayback_url': web_result.wayback_url
                }

                results[web_result.status].append(result_dict)

            else:
                verify_result = self.verify_citation(citation)

                result_dict = {
                    'citation': citation,
                    'status': verify_result.status,
                    'confidence': verify_result.confidence,
                    'issues': verify_result.issues,
                    'sources': verify_result.sources,
                    'cache_hit': verify_result.cache_hit
                }

                results[verify_result.status].append(result_dict)

        # Save verification results
        verified_file = self.citations_dir / "verification" / "verified.json"
        verified_file.parent.mkdir(parents=True, exist_ok=True)

        with open(verified_file, 'w') as f:
            json.dump(results, f, indent=2)

        print(f"\n{'='*60}")
        print("VERIFICATION SUMMARY")
        print(f"{'='*60}")
        print(f"Total:      {results['total']}")
        print(f"✓ Verified: {len(results['verified'])} ({len(results['verified'])/results['total']*100:.1f}%)")
        print(f"⚠ Uncertain: {len(results['uncertain'])}")
        print(f"✗ Failed:   {len(results['failed'])}")

        return results

    def download_all_abstracts(self) -> Dict:
        """Download abstracts for all verified citations"""

        # Load citations
        citations_file = self.citations_dir / "extraction" / "inline-citations.json"

        with open(citations_file, 'r') as f:
            data = json.load(f)

        citations = data.get('citations', [])

        # Filter to academic papers only
        academic_citations = [c for c in citations if c.get('type') == 'academic']

        print(f"Downloading abstracts for {len(academic_citations)} academic papers...")

        results = {
            'downloaded': [],
            'failed': [],
            'total': len(academic_citations)
        }

        for citation in academic_citations:
            citation_key = citation.get('citation_key', 'unknown')

            abstract = self.download_abstract(citation)

            if abstract:
                results['downloaded'].append(citation_key)
            else:
                results['failed'].append(citation_key)
                print(f"  ✗ Failed to download abstract for {citation_key}")

        print(f"\n{'='*60}")
        print("ABSTRACT DOWNLOAD SUMMARY")
        print(f"{'='*60}")
        print(f"Total:      {results['total']}")
        print(f"✓ Downloaded: {len(results['downloaded'])} ({len(results['downloaded'])/results['total']*100:.1f}%)")
        print(f"✗ Failed:   {len(results['failed'])}")

        return results

    def generate_bibtex(self) -> str:
        """Generate BibTeX file from verified citations"""

        # Load verification results
        verified_file = self.citations_dir / "verification" / "verified.json"

        if not verified_file.exists():
            print("ERROR: No verification results found. Run 'verify' first.")
            return ""

        with open(verified_file, 'r') as f:
            results = json.load(f)

        # Generate BibTeX entries
        bibtex_entries = []

        for verified in results.get('verified', []):
            citation = verified['citation']
            metadata = verified.get('metadata', {})

            entry = self.generate_bibtex_entry(citation, metadata)
            bibtex_entries.append(entry)

        # Combine into BibTeX file
        bibtex_content = f"""% BibTeX bibliography for {self.project_name}
% Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
% Verified citations: {len(bibtex_entries)}

"""
        bibtex_content += "\n".join(bibtex_entries)

        # Save to project
        bibtex_file = self.citations_dir / f"{self.project_name}-references.bib"

        with open(bibtex_file, 'w') as f:
            f.write(bibtex_content)

        print(f"✓ BibTeX file generated: {bibtex_file}")
        print(f"  Entries: {len(bibtex_entries)}")

        return str(bibtex_file)

    def show_status(self):
        """Show verification status dashboard"""

        # Load verification results if available
        verified_file = self.citations_dir / "verification" / "verified.json"

        if not verified_file.exists():
            print("No verification results found. Run 'verify' first.")
            return

        with open(verified_file, 'r') as f:
            results = json.load(f)

        print(f"\n{'='*60}")
        print(f"{self.project_name.upper()} - CITATION STATUS")
        print(f"{'='*60}\n")

        print(f"Total Citations: {results['total']}")
        print(f"✓ Verified:      {len(results['verified'])} ({len(results['verified'])/results['total']*100:.1f}%)")
        print(f"⚠ Uncertain:     {len(results['uncertain'])}")
        print(f"✗ Failed:        {len(results['failed'])}")

        if results['uncertain']:
            print("\nUNCERTAIN CITATIONS:")
            for item in results['uncertain']:
                citation = item['citation']
                print(f"  - {citation.get('citation_key')}: {', '.join(item['issues'])}")

        if results['failed']:
            print("\nFAILED CITATIONS:")
            for item in results['failed']:
                citation = item['citation']
                print(f"  - {citation.get('citation_key')}: {', '.join(item['issues'])}")


def main():
    parser = argparse.ArgumentParser(
        description="IMMUNOS Citation Verification System",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    subparsers = parser.add_subparsers(dest='command', help='Command to execute')

    # Verify command
    verify_parser = subparsers.add_parser('verify', help='Verify citations')
    verify_parser.add_argument('--project', required=True, help='Project name')
    verify_parser.add_argument('--priority', choices=['high', 'medium', 'low'], help='Filter by priority')

    # Download abstracts command
    download_parser = subparsers.add_parser('download-abstracts', help='Download abstracts')
    download_parser.add_argument('--project', required=True, help='Project name')

    # Generate BibTeX command
    bibtex_parser = subparsers.add_parser('generate-bibtex', help='Generate BibTeX file')
    bibtex_parser.add_argument('--project', required=True, help='Project name')

    # Status command
    status_parser = subparsers.add_parser('status', help='Show verification status')
    status_parser.add_argument('--project', required=True, help='Project name')

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return

    verifier = CitationVerifier(args.project)

    if args.command == 'verify':
        verifier.verify_all(priority=args.priority)

    elif args.command == 'download-abstracts':
        verifier.download_all_abstracts()

    elif args.command == 'generate-bibtex':
        verifier.generate_bibtex()

    elif args.command == 'status':
        verifier.show_status()


if __name__ == '__main__':
    main()

```
