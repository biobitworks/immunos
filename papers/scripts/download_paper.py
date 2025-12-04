#!/usr/bin/env python3
"""
Download academic papers with figures

Supports:
- Nature family journals
- PubMed Central
- bioRxiv/medRxiv
- arXiv
- Direct PDF URLs
"""

import requests
import json
import re
import sqlite3
from pathlib import Path
from datetime import datetime
from urllib.parse import urlparse
import argparse

def sanitize_doi(doi):
    """Convert DOI to safe folder name"""
    return doi.replace('/', '_').replace(':', '_')

def fetch_doi_metadata(doi):
    """Fetch metadata from DOI.org API"""
    url = f"https://doi.org/api/handles/{doi}"
    try:
        response = requests.get(url, headers={'Accept': 'application/json'})
        if response.status_code == 200:
            return response.json()
    except:
        pass
    return None

def fetch_crossref_metadata(doi):
    """Fetch detailed metadata from Crossref"""
    url = f"https://api.crossref.org/works/{doi}"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            message = data.get('message', {})

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
                'doi': doi
            }
    except Exception as e:
        print(f"⚠️  Crossref error: {e}")
    return None

def download_nature_paper(doi, output_dir):
    """Download Nature family paper PDF"""
    # Nature provides PDF at predictable URLs
    doi_clean = doi.replace('/', '_')
    pdf_url = f"https://www.nature.com/articles/{doi}.pdf"

    try:
        response = requests.get(pdf_url, headers={
            'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36'
        })

        if response.status_code == 200:
            pdf_path = output_dir / 'paper.pdf'
            with open(pdf_path, 'wb') as f:
                f.write(response.content)
            print(f"✓ Downloaded PDF: {pdf_path}")
            return str(pdf_path)
        else:
            print(f"⚠️  PDF download failed: HTTP {response.status_code}")
            return None
    except Exception as e:
        print(f"❌ Error downloading PDF: {e}")
        return None

def download_nature_figures(doi, output_dir):
    """Download figures from Nature paper"""
    # Nature figure URLs pattern
    base_url = f"https://www.nature.com/articles/{doi}/figures/"

    # Try common figure numbers (1-20)
    figures = []
    for i in range(1, 21):
        fig_url = f"https://media.springernature.com/full/springer-static/image/art%3A{doi}/MediaObjects/{doi.split('/')[-1]}_Fig{i}_HTML.png"

        try:
            response = requests.head(fig_url)
            if response.status_code == 200:
                # Download figure
                fig_response = requests.get(fig_url)
                fig_path = output_dir / f"figure_{i}.png"
                with open(fig_path, 'wb') as f:
                    f.write(fig_response.content)
                figures.append(str(fig_path))
                print(f"✓ Downloaded: figure_{i}.png")
        except:
            # Figure doesn't exist, stop trying
            if i > 3:  # At least try first 3
                break

    return figures

def create_metadata_json(metadata, pdf_path, figures, output_dir):
    """Create metadata.json file"""
    meta = {
        'doi': metadata.get('doi'),
        'title': metadata.get('title'),
        'authors': metadata.get('authors', []),
        'journal': metadata.get('journal'),
        'year': metadata.get('year'),
        'url': f"https://doi.org/{metadata.get('doi')}",
        'downloaded': datetime.now().strftime('%Y-%m-%d'),
        'pdf_path': pdf_path,
        'figures': [Path(f).name for f in figures],
        'note_path': str(output_dir / 'note.md')
    }

    meta_path = output_dir / 'metadata.json'
    with open(meta_path, 'w') as f:
        json.dump(meta, f, indent=2)

    print(f"✓ Created: metadata.json")
    return meta

def create_literature_note(metadata, output_dir):
    """Create Obsidian-compatible literature note"""

    authors_str = ', '.join(metadata.get('authors', [])[:3])
    if len(metadata.get('authors', [])) > 3:
        authors_str += ' et al.'

    note = f"""---
title: "{metadata.get('title')}"
authors: {json.dumps(metadata.get('authors', []))}
year: {metadata.get('year')}
journal: "{metadata.get('journal')}"
doi: "{metadata.get('doi')}"
url: "https://doi.org/{metadata.get('doi')}"
type: journal-article
downloaded: "{metadata.get('downloaded')}"
tags:
  - paper
  - {metadata.get('journal', 'unknown').lower().replace(' ', '-')}
---

# {metadata.get('title')}

**Authors**: {authors_str}
**Journal**: {metadata.get('journal')} ({metadata.get('year')})
**DOI**: [{metadata.get('doi')}](https://doi.org/{metadata.get('doi')})

## Abstract

[Add abstract here]

## Key Findings

[Summarize key findings]

## Figures

"""

    # Add figure references
    for i, fig in enumerate(metadata.get('figures', []), 1):
        note += f"### Figure {i}\n![Figure {i}]({fig})\n\n"

    note += """## Notes

[Your notes here]

## Related

- [[../INDEX|Papers Dashboard]]
- [[../../research/INDEX|Aging Research Dashboard]]

---

**PDF**: [paper.pdf](paper.pdf)
**Downloaded**: {downloaded}
**Metadata**: [metadata.json](metadata.json)
""".format(downloaded=metadata.get('downloaded'))

    note_path = output_dir / 'note.md'
    with open(note_path, 'w') as f:
        f.write(note)

    print(f"✓ Created: note.md")
    return str(note_path)

def add_to_database(metadata, db_path="../papers.db"):
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
        metadata.get('pdf_path'),
        metadata.get('downloaded'),
        metadata.get('note_path'),
        json.dumps([])  # tags placeholder
    ))

    # Insert figures
    for i, fig in enumerate(metadata.get('figures', []), 1):
        cursor.execute("""
            INSERT INTO figures (paper_doi, figure_num, filename, file_path)
            VALUES (?, ?, ?, ?)
        """, (
            metadata.get('doi'),
            i,
            Path(fig).name,
            fig
        ))

    conn.commit()
    conn.close()
    print(f"✓ Added to database: papers.db")

def download_paper(url_or_doi, output_base="../", tags=None):
    """Main download function"""

    # Extract DOI from URL or use directly
    if url_or_doi.startswith('http'):
        # Extract DOI from URL
        match = re.search(r'10\.\d{4,}/[^\s]+', url_or_doi)
        if match:
            doi = match.group(0).rstrip('/')
        else:
            print("❌ Could not extract DOI from URL")
            return False
    else:
        doi = url_or_doi

    print(f"\n{'='*70}")
    print(f"Downloading: {doi}")
    print(f"{'='*70}\n")

    # Fetch metadata
    print("Fetching metadata from Crossref...")
    metadata = fetch_crossref_metadata(doi)

    if not metadata:
        print("❌ Could not fetch metadata")
        return False

    print(f"✓ Title: {metadata.get('title')}")
    print(f"✓ Journal: {metadata.get('journal')}")
    print(f"✓ Year: {metadata.get('year')}")
    print(f"✓ Authors: {len(metadata.get('authors', []))} authors\n")

    # Create output directory
    safe_doi = sanitize_doi(doi)
    output_dir = Path(output_base) / safe_doi
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"✓ Created directory: {output_dir}\n")

    # Download PDF
    print("Downloading PDF...")
    pdf_path = download_nature_paper(doi, output_dir)

    if not pdf_path:
        print("⚠️  PDF download failed, continuing with metadata only")

    # Download figures
    print("\nDownloading figures...")
    figures = download_nature_figures(doi, output_dir)
    print(f"✓ Downloaded {len(figures)} figures\n")

    # Create metadata
    print("Creating metadata...")
    metadata_full = create_metadata_json(
        metadata, pdf_path, figures, output_dir
    )

    # Create literature note
    print("Creating literature note...")
    note_path = create_literature_note(metadata_full, output_dir)

    # Add to database
    print("Adding to database...")
    add_to_database(metadata_full)

    print(f"\n{'='*70}")
    print(f"✓ COMPLETE: {safe_doi}/")
    print(f"{'='*70}\n")
    print(f"Files created:")
    print(f"  - {output_dir}/paper.pdf")
    print(f"  - {output_dir}/metadata.json")
    print(f"  - {output_dir}/note.md")
    print(f"  - {len(figures)} figure(s)")
    print(f"\nView in Obsidian: [[papers/{safe_doi}/note|{metadata.get('title')[:50]}...]]")

    return True

def main():
    parser = argparse.ArgumentParser(description='Download academic paper with figures')
    parser.add_argument('url', help='Paper URL or DOI')
    parser.add_argument('--tags', help='Comma-separated tags', default='')
    parser.add_argument('--output', '-o', help='Output directory', default='../')

    args = parser.parse_args()

    tags = [t.strip() for t in args.tags.split(',') if t.strip()]
    download_paper(args.url, args.output, tags)

if __name__ == "__main__":
    main()
