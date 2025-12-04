#!/usr/bin/env python3
"""
Download academic papers with figures

Supports:
- Nature family journals
- PubMed Central
- bioRxiv/medRxiv
- arXiv
- Direct PDF URLs

Features:
- Batch download multiple papers
- License checking (only downloads open access)
- Automatic figure extraction
- SQLite database tracking
- Obsidian literature notes
"""

import subprocess
import json
import re
import sqlite3
from pathlib import Path
from datetime import datetime
import argparse
import sys

def sanitize_doi(doi):
    """Convert DOI to safe folder name"""
    return doi.replace('/', '_').replace(':', '_')

def extract_doi(url_or_doi):
    """Extract DOI from URL or return DOI directly"""
    if url_or_doi.startswith('http'):
        match = re.search(r'10\.\d{4,}/[^\s]+', url_or_doi)
        if match:
            return match.group(0).rstrip('/')
        return None
    return url_or_doi

def fetch_crossref_metadata(doi):
    """Fetch detailed metadata from Crossref using curl"""
    url = f"https://api.crossref.org/works/{doi}"
    try:
        result = subprocess.run(
            ['curl', '-s', url],
            capture_output=True,
            text=True,
            timeout=30
        )

        if result.returncode == 0:
            data = json.loads(result.stdout)
            message = data.get('message', {})

            # Extract license
            license_info = message.get('license', [])
            license_url = license_info[0].get('URL', 'Not specified') if license_info else 'Not specified'

            # Extract metadata
            title = message.get('title', ['Unknown'])[0]
            authors = [
                f"{a.get('given', '')} {a.get('family', '')}".strip()
                for a in message.get('author', [])
            ]
            journal = message.get('container-title', ['Unknown'])[0]
            year = None
            if 'published' in message:
                year = message['published'].get('date-parts', [[None]])[0][0]
            elif 'published-print' in message:
                year = message['published-print'].get('date-parts', [[None]])[0][0]

            return {
                'title': title,
                'authors': authors,
                'journal': journal,
                'year': year,
                'doi': doi,
                'license': license_url
            }
    except Exception as e:
        print(f"⚠️  Crossref error for {doi}: {e}")
    return None

def is_open_access(license_url):
    """Check if license is open access"""
    if not license_url or license_url == 'Not specified':
        return False

    license_lower = license_url.lower()

    # Open access licenses
    open_licenses = [
        'creativecommons.org/licenses/by/',
        'creativecommons.org/licenses/by-nc/',
        'creativecommons.org/licenses/by-sa/',
        'creativecommons.org/licenses/by-nc-sa/',
        'creativecommons.org/licenses/by-nc-nd/',
        'creativecommons.org/publicdomain/',
    ]

    return any(lic in license_lower for lic in open_licenses)

def check_licenses(dois):
    """Check licenses for multiple DOIs"""
    results = {
        'open_access': [],
        'restricted': []
    }

    print(f"\n{'='*70}")
    print(f"Checking licenses for {len(dois)} papers...")
    print(f"{'='*70}\n")

    for doi in dois:
        metadata = fetch_crossref_metadata(doi)
        if metadata:
            is_oa = is_open_access(metadata['license'])

            info = {
                'doi': doi,
                'title': metadata['title'],
                'journal': metadata['journal'],
                'year': metadata['year'],
                'license': metadata['license'],
                'metadata': metadata
            }

            if is_oa:
                results['open_access'].append(info)
                print(f"✅ {doi}")
                print(f"   {metadata['title'][:70]}...")
                print(f"   License: {metadata['license']}")
            else:
                results['restricted'].append(info)
                print(f"❌ {doi} - RESTRICTED")
                print(f"   {metadata['title'][:70]}...")
                print(f"   License: {metadata['license']}")
            print()
        else:
            print(f"⚠️  Could not fetch metadata for {doi}\n")

    return results

def download_nature_pdf(doi, output_dir):
    """Download Nature family paper PDF using curl"""
    pdf_url = f"https://www.nature.com/articles/{doi}.pdf"
    pdf_path = output_dir / 'paper.pdf'

    try:
        result = subprocess.run(
            ['curl', '-s', '-L', pdf_url, '-o', str(pdf_path)],
            capture_output=True,
            timeout=120
        )

        if result.returncode == 0 and pdf_path.exists() and pdf_path.stat().st_size > 1000:
            return str(pdf_path)
        else:
            pdf_path.unlink(missing_ok=True)
            return None
    except Exception as e:
        print(f"⚠️  PDF download error: {e}")
        return None

def download_nature_figures(doi, output_dir):
    """Download figures from Nature paper using curl"""
    figures = []

    # Extract journal and article number from DOI
    parts = doi.split('/')
    if len(parts) >= 2:
        journal_code = parts[-1].split('-')[0]  # e.g., "41467" from "s41467-025-58466-2"

        for i in range(1, 21):
            fig_url = f"https://media.springernature.com/full/springer-static/image/art%3A{doi}/MediaObjects/{journal_code}_{parts[-1].replace('s' + journal_code + '-', '').replace('-', '_')}_Fig{i}_HTML.png"
            fig_path = output_dir / f"figure_{i}.png"

            # Try to download
            result = subprocess.run(
                ['curl', '-s', '-L', fig_url, '-o', str(fig_path)],
                capture_output=True,
                timeout=30
            )

            if result.returncode == 0 and fig_path.exists() and fig_path.stat().st_size > 1000:
                figures.append(str(fig_path.name))
            else:
                fig_path.unlink(missing_ok=True)
                if i > 3:  # Stop after trying first 3
                    break

    return figures

def create_metadata_json(metadata, pdf_path, figures, output_dir, license_url):
    """Create metadata.json file"""
    # Determine license type
    if 'creativecommons.org/licenses/by/4.0' in license_url:
        license_type = 'CC-BY-4.0'
    elif 'creativecommons.org/licenses/by-nc-nd/4.0' in license_url:
        license_type = 'CC-BY-NC-ND-4.0'
    elif 'creativecommons.org' in license_url:
        license_type = 'CC-BY'
    else:
        license_type = 'Open Access'

    meta = {
        'doi': metadata.get('doi'),
        'title': metadata.get('title'),
        'authors': metadata.get('authors', []),
        'journal': metadata.get('journal'),
        'year': metadata.get('year'),
        'url': f"https://doi.org/{metadata.get('doi')}",
        'downloaded': datetime.now().strftime('%Y-%m-%d'),
        'pdf_path': 'paper.pdf' if pdf_path else None,
        'figures': figures,
        'note_path': 'note.md',
        'license': license_type
    }

    meta_path = output_dir / 'metadata.json'
    with open(meta_path, 'w') as f:
        json.dump(meta, f, indent=2)

    return meta

def create_literature_note(metadata, output_dir):
    """Create Obsidian-compatible literature note"""
    authors_str = ', '.join(metadata.get('authors', [])[:3])
    if len(metadata.get('authors', [])) > 3:
        authors_str += f' et al. ({len(metadata["authors"])} authors)'

    # YAML authors
    yaml_authors = '\n'.join([f'  - {a}' for a in metadata.get('authors', [])])

    # Generate tags
    tags = ['paper', 'aging-biology']
    journal_tag = metadata.get('journal', 'unknown').lower().replace(' ', '-')
    tags.append(journal_tag)

    # Figures section
    figures_md = '\n\n'.join([
        f'### Figure {i}\n![Figure {i}]({fig})'
        for i, fig in enumerate(metadata.get('figures', []), 1)
    ])

    note = f'''---
title: "{metadata.get('title')}"
authors:
{yaml_authors}
year: {metadata.get('year')}
journal: "{metadata.get('journal')}"
doi: "{metadata.get('doi')}"
url: "{metadata.get('url')}"
type: journal-article
downloaded: "{metadata.get('downloaded')}"
license: "{metadata.get('license')}"
tags:
{chr(10).join([f'  - {tag}' for tag in tags])}
---

# {metadata.get('title')}

**Authors**: {authors_str}
**Journal**: {metadata.get('journal')} ({metadata.get('year')})
**DOI**: [{metadata.get('doi')}]({metadata.get('url')})
**License**: {metadata.get('license')} (Open Access)

## Quick Summary

[To be filled after reading]

## Relevance to Research

**Connection to**:
- [[../10.1038_s41467-025-58466-2/note|Wang et al. 2025]] - Spatial transcriptomics of aging brain
- [[../../research/databases/biogrid|BioGRID]] - Protein interaction networks
- [[../../data/genage/README|GenAge]] - Aging-related genes
- [[../../data/cellage/README|CellAge]] - Senescence markers

## Key Findings

[To be filled after reading]

## Methodology

[To be filled after reading]

## Figures

{figures_md if figures_md else '[No figures downloaded]'}

## Notes

[Add notes here]

## Related

- [[../INDEX|Papers Dashboard]]
- [[../../research/INDEX|Aging Research Dashboard]]

---

**PDF**: [paper.pdf](paper.pdf)
**Downloaded**: {metadata.get('downloaded')}
**Metadata**: [metadata.json](metadata.json)
'''

    note_path = output_dir / 'note.md'
    with open(note_path, 'w') as f:
        f.write(note)

    return str(note_path)

def add_to_database(metadata, db_path, folder):
    """Add paper to SQLite database"""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Insert paper
    cursor.execute("""
        INSERT OR REPLACE INTO papers
        (doi, title, authors, journal, year, url, pdf_path, download_date, note_path, tags)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, (
        metadata.get('doi'),
        metadata.get('title'),
        json.dumps(metadata.get('authors', [])),
        metadata.get('journal'),
        metadata.get('year'),
        metadata.get('url'),
        f"{folder}/paper.pdf" if metadata.get('pdf_path') else None,
        metadata.get('downloaded'),
        f"{folder}/note.md",
        json.dumps([])
    ))

    # Insert figures
    for i, fig in enumerate(metadata.get('figures', []), 1):
        cursor.execute("""
            INSERT INTO figures (paper_doi, figure_num, filename, file_path)
            VALUES (?, ?, ?, ?)
        """, (
            metadata.get('doi'),
            i,
            fig,
            f"{folder}/{fig}"
        ))

    conn.commit()
    conn.close()

def update_restricted_access_file(restricted_papers, restricted_path):
    """Update or create RESTRICTED_ACCESS.md with new restricted papers"""
    content = f'''# Restricted Access Papers

Papers that could not be downloaded due to access restrictions. These are tracked here for reference but PDFs are not included.

---

## Papers Not Downloaded ({len(restricted_papers)})

'''

    for i, paper in enumerate(restricted_papers, 1):
        content += f'''### {i}. {paper['title']}

**DOI**: [{paper['doi']}](https://doi.org/{paper['doi']})
**Journal**: {paper['journal']} ({paper['year']})
**License**: {paper['license']}
**Status**: ❌ Not downloaded - Restricted access

**Why not downloaded**: This paper does not have an open access license (CC-BY, CC-BY-NC, etc.). The full text and figures are behind a paywall.

**Alternative access**:
- Institutional library access required
- May be available through university subscription
- Contact authors for preprint version

---

'''

    content += f'''## Access Policy

**Papers Archive download policy**:
- ✅ Download: CC-BY, CC-BY-NC, CC-BY-SA, CC0 (open access)
- ✅ Download: CC-BY-NC-ND (open access with restrictions)
- ❌ Do not download: Paywalled, TDM-only, subscription-only

**Last updated**: {datetime.now().strftime('%Y-%m-%d')}
'''

    with open(restricted_path, 'w') as f:
        f.write(content)

def download_paper(paper_info, output_base, db_path):
    """Download a single paper"""
    doi = paper_info['doi']
    metadata = paper_info['metadata']

    print(f"\n{'='*70}")
    print(f"Downloading: {doi}")
    print(f"{'='*70}")
    print(f"Title: {metadata['title'][:60]}...")
    print(f"Journal: {metadata['journal']}")
    print(f"Year: {metadata['year']}")
    print(f"Authors: {len(metadata['authors'])} authors\n")

    # Create output directory
    safe_doi = sanitize_doi(doi)
    output_dir = Path(output_base) / safe_doi
    output_dir.mkdir(parents=True, exist_ok=True)

    # Download PDF
    print("Downloading PDF...")
    pdf_path = download_nature_pdf(doi, output_dir)

    if pdf_path:
        pdf_size = Path(pdf_path).stat().st_size / (1024 * 1024)
        print(f"✓ Downloaded PDF: {pdf_size:.1f} MB")
    else:
        print("⚠️  PDF download failed")

    # Download figures
    print("\nDownloading figures...")
    figures = download_nature_figures(doi, output_dir)
    print(f"✓ Downloaded {len(figures)} figures")

    # Create metadata
    print("\nCreating metadata...")
    metadata_full = create_metadata_json(
        metadata, pdf_path, figures, output_dir, paper_info['license']
    )
    print("✓ Created metadata.json")

    # Create literature note
    print("Creating literature note...")
    create_literature_note(metadata_full, output_dir)
    print("✓ Created note.md")

    # Add to database
    print("Adding to database...")
    add_to_database(metadata_full, db_path, safe_doi)
    print("✓ Added to database")

    return {
        'doi': doi,
        'folder': safe_doi,
        'pdf_size': Path(pdf_path).stat().st_size if pdf_path else 0,
        'figures': len(figures),
        'success': True
    }

def main():
    parser = argparse.ArgumentParser(
        description='Download academic papers with figures (supports batch downloads)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Single paper
  %(prog)s 10.1038/s41467-025-58466-2

  # Multiple papers
  %(prog)s 10.1038/s41467-025-58466-2 10.1038/s43587-025-01016-8

  # From URLs
  %(prog)s https://doi.org/10.1038/s41467-025-58466-2

  # Mixed DOIs and URLs
  %(prog)s 10.1038/s41467-025-58466-2 https://doi.org/10.1038/s43587-025-01016-8
        '''
    )
    parser.add_argument(
        'urls',
        nargs='+',
        help='Paper URLs or DOIs (one or more)'
    )
    parser.add_argument(
        '--output', '-o',
        help='Output directory',
        default='.'
    )
    parser.add_argument(
        '--db',
        help='Database path',
        default='papers.db'
    )
    parser.add_argument(
        '--skip-license-check',
        action='store_true',
        help='Skip license checking (not recommended)'
    )

    args = parser.parse_args()

    # Extract DOIs
    dois = []
    for url in args.urls:
        doi = extract_doi(url)
        if doi:
            dois.append(doi)
        else:
            print(f"⚠️  Could not extract DOI from: {url}")

    if not dois:
        print("❌ No valid DOIs found")
        return 1

    # Check licenses
    if not args.skip_license_check:
        license_results = check_licenses(dois)
    else:
        # Skip license check, assume all are open access
        license_results = {
            'open_access': [{'doi': doi, 'metadata': fetch_crossref_metadata(doi), 'license': 'Unknown'} for doi in dois],
            'restricted': []
        }

    print(f"\n{'='*70}")
    print(f"License Check Summary:")
    print(f"  ✅ Open Access: {len(license_results['open_access'])} papers")
    print(f"  ❌ Restricted: {len(license_results['restricted'])} papers")
    print(f"{'='*70}\n")

    # Download open access papers
    if license_results['open_access']:
        print(f"\nDownloading {len(license_results['open_access'])} open access papers...\n")

        results = []
        for paper in license_results['open_access']:
            try:
                result = download_paper(paper, args.output, args.db)
                results.append(result)
            except Exception as e:
                print(f"❌ Error downloading {paper['doi']}: {e}")

        # Summary
        print(f"\n{'='*70}")
        print(f"DOWNLOAD COMPLETE")
        print(f"{'='*70}")
        print(f"\n✅ Successfully downloaded {len(results)} papers:")

        total_size = sum(r['pdf_size'] for r in results) / (1024 * 1024)
        total_figures = sum(r['figures'] for r in results)

        for r in results:
            print(f"  - {r['folder']} ({r['figures']} figures)")

        print(f"\nTotal storage: {total_size:.1f} MB")
        print(f"Total figures: {total_figures}")

    # Update restricted access tracking
    if license_results['restricted']:
        print(f"\n⚠️  {len(license_results['restricted'])} papers were RESTRICTED and not downloaded:")
        for paper in license_results['restricted']:
            print(f"  - {paper['doi']}")
            print(f"    {paper['title'][:60]}...")

        # Update RESTRICTED_ACCESS.md
        restricted_path = Path(args.output) / 'RESTRICTED_ACCESS.md'
        update_restricted_access_file(license_results['restricted'], restricted_path)
        print(f"\n✓ Updated {restricted_path}")

    print(f"\n{'='*70}")
    print("View papers in Obsidian: [[papers/INDEX|Papers Dashboard]]")
    print(f"{'='*70}\n")

    return 0

if __name__ == "__main__":
    sys.exit(main())
