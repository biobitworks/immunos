---
source: /Users/byron/projects/papers/scripts/import_from_zotero.py
relative: papers/scripts/import_from_zotero.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
Import papers from Zotero backup into papers directory

Extracts metadata from Zotero SQLite database and creates:
- DOI-based directory structure
- metadata.json files
- Obsidian literature notes (note.md)
- Copies PDFs from Zotero storage
- Updates papers.db

Usage:
    python import_from_zotero.py \\
      --zotero-db PATH_TO_zotero.sqlite \\
      --bibtex-db PATH_TO_better-bibtex.sqlite \\
      --storage PATH_TO_storage/ \\
      --papers-db PATH_TO_papers.db \\
      --output PATH_TO_papers/ \\
      [--skip-duplicates] \\
      [--dry-run]
"""

import sqlite3
import json
import re
import shutil
import argparse
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Tuple


class ZoteroParser:
    """Parse Zotero SQLite database to extract papers"""

    def __init__(self, zotero_db_path: str, bibtex_db_path: str, storage_path: str):
        self.zotero_db = zotero_db_path
        self.bibtex_db = bibtex_db_path
        self.storage_path = Path(storage_path)

    def extract_all_papers(self) -> List[Dict]:
        """Extract all journal articles, books, book sections, theses"""
        conn = sqlite3.connect(self.zotero_db)
        cursor = conn.cursor()

        # Get all importable items
        cursor.execute("""
            SELECT i.itemID, it.typeName, i.key, i.dateAdded
            FROM items i
            JOIN itemTypes it ON i.itemTypeID = it.itemTypeID
            WHERE it.typeName IN ('journalArticle', 'book', 'bookSection', 'thesis')
            ORDER BY i.itemID
        """)

        items = cursor.fetchall()
        papers = []

        for item_id, type_name, key, date_added in items:
            paper = self.get_paper_metadata(item_id, type_name, key, date_added)
            if paper:
                papers.append(paper)

        conn.close()
        return papers

    def get_paper_metadata(self, item_id: int, type_name: str, key: str, date_added: str) -> Optional[Dict]:
        """Extract complete metadata for a single item"""
        conn = sqlite3.connect(self.zotero_db)
        cursor = conn.cursor()

        # Get metadata fields
        cursor.execute("""
            SELECT f.fieldName, idv.value
            FROM itemData id
            JOIN fields f ON id.fieldID = f.fieldID
            JOIN itemDataValues idv ON id.valueID = idv.valueID
            WHERE id.itemID = ?
        """, (item_id,))

        metadata_raw = dict(cursor.fetchall())

        # Get authors
        authors = self.get_authors(item_id)

        # Get PDF path
        pdf_info = self.get_pdf_path(item_id)

        # Get citation key from Better BibTeX
        citation_key = self.get_citation_key(item_id)

        conn.close()

        # Build paper object
        paper = {
            'item_id': item_id,
            'type': type_name,
            'key': key,
            'doi': metadata_raw.get('DOI'),
            'title': metadata_raw.get('title'),
            'authors': authors,
            'abstract': metadata_raw.get('abstractNote'),
            'journal': metadata_raw.get('publicationTitle'),
            'year': self.extract_year(metadata_raw.get('date', '')),
            'volume': metadata_raw.get('volume'),
            'issue': metadata_raw.get('issue'),
            'pages': metadata_raw.get('pages'),
            'url': metadata_raw.get('url'),
            'issn': metadata_raw.get('ISSN'),
            'isbn': metadata_raw.get('ISBN'),
            'pdf_storage_key': pdf_info['key'] if pdf_info else None,
            'pdf_filename': pdf_info['filename'] if pdf_info else None,
            'citation_key': citation_key,
            'date_added': date_added
        }

        return paper

    def get_authors(self, item_id: int) -> List[str]:
        """Extract authors in order"""
        conn = sqlite3.connect(self.zotero_db)
        cursor = conn.cursor()

        cursor.execute("""
            SELECT c.firstName, c.lastName, c.fieldMode
            FROM itemCreators ic
            JOIN creators c ON ic.creatorID = c.creatorID
            WHERE ic.itemID = ? AND ic.creatorTypeID = 8
            ORDER BY ic.orderIndex
        """, (item_id,))

        authors = []
        for first, last, field_mode in cursor.fetchall():
            if field_mode == 1:
                # Single field name
                authors.append(last)
            else:
                # Standard first + last
                name = f"{first} {last}".strip()
                if name:
                    authors.append(name)

        conn.close()
        return authors

    def get_pdf_path(self, item_id: int) -> Optional[Dict]:
        """Find PDF attachment in storage"""
        conn = sqlite3.connect(self.zotero_db)
        cursor = conn.cursor()

        cursor.execute("""
            SELECT ia.itemID, i.key, ia.path
            FROM itemAttachments ia
            JOIN items i ON ia.itemID = i.itemID
            WHERE ia.parentItemID = ?
            AND ia.contentType = 'application/pdf'
            LIMIT 1
        """, (item_id,))

        result = cursor.fetchone()
        conn.close()

        if result:
            attachment_id, storage_key, path = result
            # path format: "storage:filename.pdf"
            if path and path.startswith('storage:'):
                filename = path.replace('storage:', '')
                return {'key': storage_key, 'filename': filename}

        return None

    def get_citation_key(self, item_id: int) -> Optional[str]:
        """Get Better BibTeX citation key"""
        try:
            conn = sqlite3.connect(self.bibtex_db)
            cursor = conn.cursor()

            cursor.execute("""
                SELECT citationKey
                FROM citationkey
                WHERE itemID = ?
            """, (item_id,))

            result = cursor.fetchone()
            conn.close()

            return result[0] if result else None
        except:
            return None

    def extract_year(self, date_string: str) -> Optional[int]:
        """Extract 4-digit year from date string"""
        if not date_string:
            return None
        match = re.search(r'(\d{4})', date_string)
        return int(match.group(1)) if match else None


def sanitize_doi(doi: str) -> str:
    """Convert DOI to safe folder name"""
    return doi.replace('/', '_').replace(':', '_')


def generate_folder_name(paper: Dict) -> str:
    """Generate folder name for papers without DOI"""
    if paper.get('citation_key'):
        return paper['citation_key']

    # Fallback: author_year_title-slug
    first_author = paper['authors'][0].split()[-1].lower() if paper.get('authors') else 'unknown'
    year = paper.get('year') or 'unknown'
    title_slug = re.sub(r'[^\w\s-]', '', paper.get('title', 'untitled')[:30].lower())
    title_slug = re.sub(r'[-\s]+', '_', title_slug)

    return f"{first_author}_{year}_{title_slug}"


def check_duplicate(paper: Dict, papers_db_path: str) -> Tuple[bool, Optional[str]]:
    """
    Check if paper already exists in database

    Returns:
        (is_duplicate, match_type)
    """
    if not paper.get('doi'):
        return False, None

    conn = sqlite3.connect(papers_db_path)
    cursor = conn.cursor()

    # Check DOI match
    cursor.execute("SELECT doi FROM papers WHERE doi = ?", (paper['doi'],))
    result = cursor.fetchone()

    conn.close()

    if result:
        return True, 'doi_match'

    return False, None


def copy_pdf_from_storage(paper: Dict, output_dir: Path, storage_base: Path) -> Optional[str]:
    """Copy PDF from Zotero storage to paper directory"""
    if not paper.get('pdf_storage_key'):
        return None

    # Zotero storage path format: storage/{KEY}/{filename}.pdf
    storage_dir = storage_base / paper['pdf_storage_key']

    if not storage_dir.exists():
        print(f"⚠️  Storage directory not found: {storage_dir}")
        return None

    # Find PDF in storage directory
    pdf_files = list(storage_dir.glob('*.pdf'))

    if not pdf_files:
        print(f"⚠️  No PDF found in {storage_dir}")
        return None

    source_pdf = pdf_files[0]
    dest_pdf = output_dir / 'paper.pdf'

    try:
        shutil.copy2(source_pdf, dest_pdf)
        print(f"    ✓ Copied PDF: {source_pdf.name}")
        return str(dest_pdf)
    except Exception as e:
        print(f"    ⚠️  PDF copy failed: {e}")
        return None


def create_metadata_json(paper: Dict, pdf_path: Optional[str], output_dir: Path) -> Dict:
    """Create metadata.json following existing format"""
    metadata = {
        'doi': paper.get('doi'),
        'title': paper['title'],
        'authors': paper['authors'],
        'journal': paper.get('journal'),
        'year': paper.get('year'),
        'url': paper.get('url') or (f"https://doi.org/{paper['doi']}" if paper.get('doi') else None),
        'downloaded': datetime.now().strftime('%Y-%m-%d'),
        'pdf_path': 'paper.pdf' if pdf_path else None,
        'figures': [],
        'note_path': 'note.md'
    }

    # Add optional fields
    if paper.get('abstract'):
        metadata['abstract'] = paper['abstract']
    if paper.get('volume'):
        metadata['volume'] = paper['volume']
    if paper.get('issue'):
        metadata['issue'] = paper['issue']
    if paper.get('pages'):
        metadata['pages'] = paper['pages']
    if paper.get('citation_key'):
        metadata['citation_key'] = paper['citation_key']
    if paper.get('issn'):
        metadata['issn'] = paper['issn']
    if paper.get('isbn'):
        metadata['isbn'] = paper['isbn']

    # Save to file
    meta_path = output_dir / 'metadata.json'
    with open(meta_path, 'w') as f:
        json.dump(metadata, f, indent=2)

    return metadata


def create_literature_note(metadata: Dict, output_dir: Path):
    """Create Obsidian-compatible literature note"""
    # Format authors
    authors_str = ', '.join(metadata['authors'][:3])
    if len(metadata['authors']) > 3:
        authors_str += f' et al. ({len(metadata["authors"])} authors)'

    # YAML authors
    yaml_authors = '\n'.join([f'  - {a}' for a in metadata['authors']])

    # Generate tags
    tags = ['paper', 'zotero-import']
    if metadata.get('journal'):
        journal_tag = metadata['journal'].lower().replace(' ', '-').replace('&', 'and')
        tags.append(journal_tag)

    # Build note content
    note = f'''---
title: "{metadata['title']}"
authors:
{yaml_authors}
year: {metadata.get('year')}
journal: "{metadata.get('journal', 'Unknown')}"
doi: "{metadata.get('doi', 'N/A')}"
url: "{metadata.get('url', '')}"
type: journal-article
downloaded: "{metadata['downloaded']}"
citation_key: "{metadata.get('citation_key', '')}"
tags:
{chr(10).join([f'  - {tag}' for tag in tags])}
---

# {metadata['title']}

**Authors**: {authors_str}
**Journal**: {metadata.get('journal', 'Unknown')} ({metadata.get('year', 'N/A')})
**DOI**: [{metadata.get('doi', 'N/A')}]({metadata.get('url', '')})

## Abstract

{metadata.get('abstract', '[No abstract available]')}

## Quick Summary

[To be filled after reading]

## Relevance to Research

**Connection to**:
- [[../INDEX|Papers Dashboard]]
- [[../../research/INDEX|Aging Research Dashboard]]

## Key Findings

[To be filled after reading]

## Methodology

[To be filled after reading]

## Notes

[Add notes here]

## Related

- [[../INDEX|Papers Dashboard]]

---

**PDF**: {f'[paper.pdf](paper.pdf)' if metadata.get('pdf_path') else 'No PDF'}
**Downloaded**: {metadata['downloaded']}
**Metadata**: [metadata.json](metadata.json)
**Citation**: `@{metadata.get('citation_key', '')}`
'''

    note_path = output_dir / 'note.md'
    with open(note_path, 'w') as f:
        f.write(note)

    return str(note_path)


def add_to_database(metadata: Dict, db_path: str, folder_name: str):
    """Add paper to SQLite database"""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute("""
        INSERT OR REPLACE INTO papers
        (doi, title, authors, journal, year, url, pdf_path, download_date, note_path, tags)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, (
        metadata.get('doi'),
        metadata['title'],
        json.dumps(metadata['authors']),
        metadata.get('journal'),
        metadata.get('year'),
        metadata.get('url'),
        f"{folder_name}/paper.pdf" if metadata.get('pdf_path') else None,
        metadata['downloaded'],
        f"{folder_name}/note.md",
        json.dumps([])
    ))

    conn.commit()
    conn.close()


def import_paper(paper: Dict, output_base: Path, papers_db_path: str, storage_base: Path,
                 skip_duplicates: bool = True, dry_run: bool = False) -> Dict:
    """
    Import single paper from Zotero

    Returns status dict with: status, reason, folder_name
    """
    # Check duplicate
    is_duplicate, match_type = check_duplicate(paper, papers_db_path)
    if is_duplicate and skip_duplicates:
        return {
            'status': 'skipped',
            'reason': match_type,
            'doi': paper.get('doi', 'NO_DOI'),
            'title': paper['title'][:60]
        }

    # Generate folder name
    folder_name = sanitize_doi(paper['doi']) if paper.get('doi') else generate_folder_name(paper)
    output_dir = output_base / folder_name

    if dry_run:
        return {
            'status': 'would_import',
            'folder_name': folder_name,
            'doi': paper.get('doi', 'NO_DOI'),
            'title': paper['title'][:60],
            'has_pdf': bool(paper.get('pdf_storage_key'))
        }

    # Create directory
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n📄 {folder_name}")
    print(f"   {paper['title'][:70]}...")

    # Copy PDF
    pdf_path = copy_pdf_from_storage(paper, output_dir, storage_base)

    # Generate metadata.json
    metadata = create_metadata_json(paper, pdf_path, output_dir)
    print(f"    ✓ Created metadata.json")

    # Generate note.md
    create_literature_note(metadata, output_dir)
    print(f"    ✓ Created note.md")

    # Update database
    add_to_database(metadata, papers_db_path, folder_name)
    print(f"    ✓ Updated papers.db")

    return {
        'status': 'imported',
        'folder_name': folder_name,
        'doi': paper.get('doi', 'NO_DOI'),
        'title': paper['title'][:60]
    }


def main():
    parser = argparse.ArgumentParser(description='Import papers from Zotero backup')
    parser.add_argument('--zotero-db', required=True, help='Path to zotero.sqlite')
    parser.add_argument('--bibtex-db', required=True, help='Path to better-bibtex.sqlite')
    parser.add_argument('--storage', required=True, help='Path to Zotero storage directory')
    parser.add_argument('--papers-db', required=True, help='Path to papers.db')
    parser.add_argument('--output', required=True, help='Output directory for papers')
    parser.add_argument('--skip-duplicates', action='store_true', default=True,
                        help='Skip papers already in database (default: True)')
    parser.add_argument('--dry-run', action='store_true', help='Preview import without making changes')

    args = parser.parse_args()

    # Initialize parser
    parser_obj = ZoteroParser(args.zotero_db, args.bibtex_db, args.storage)

    print(f"\n{'='*70}")
    print(f"Zotero Import {'(DRY RUN)' if args.dry_run else ''}")
    print(f"{'='*70}\n")

    # Extract papers
    print("Extracting papers from Zotero database...")
    papers = parser_obj.extract_all_papers()
    print(f"Found {len(papers)} items\n")

    # Process each paper
    output_base = Path(args.output)
    storage_base = Path(args.storage)

    results = {
        'imported': [],
        'skipped': [],
        'would_import': [],
        'errors': []
    }

    for paper in papers:
        try:
            result = import_paper(
                paper,
                output_base,
                args.papers_db,
                storage_base,
                skip_duplicates=args.skip_duplicates,
                dry_run=args.dry_run
            )

            status = result['status']
            results[status].append(result)

        except Exception as e:
            print(f"\n❌ Error importing {paper.get('doi', 'unknown')}: {e}")
            results['errors'].append({
                'doi': paper.get('doi', 'NO_DOI'),
                'title': paper.get('title', 'Unknown')[:60],
                'error': str(e)
            })

    # Print summary
    print(f"\n{'='*70}")
    print(f"{'DRY RUN ' if args.dry_run else ''}Summary")
    print(f"{'='*70}\n")

    print(f"Total papers in Zotero: {len(papers)}")

    if args.dry_run:
        print(f"Would import: {len(results['would_import'])}")
    else:
        print(f"Imported: {len(results['imported'])}")

    print(f"Skipped (duplicates): {len(results['skipped'])}")
    print(f"Errors: {len(results['errors'])}")

    if results['skipped']:
        print(f"\n{'='*70}")
        print("Skipped (already in database):")
        print(f"{'='*70}")
        for item in results['skipped'][:10]:
            print(f"  ⏭️  {item['doi']}: {item['title']}...")
        if len(results['skipped']) > 10:
            print(f"  ... and {len(results['skipped']) - 10} more")

    if results['errors']:
        print(f"\n{'='*70}")
        print("Errors:")
        print(f"{'='*70}")
        for item in results['errors']:
            print(f"  ❌ {item['doi']}: {item['error']}")

    if args.dry_run and results['would_import']:
        print(f"\n{'='*70}")
        print("Would import:")
        print(f"{'='*70}")
        for item in results['would_import'][:20]:
            pdf_status = "with PDF" if item['has_pdf'] else "metadata only"
            print(f"  ✓ {item['doi']}: {item['title']}... ({pdf_status})")
        if len(results['would_import']) > 20:
            print(f"  ... and {len(results['would_import']) - 20} more")

    print()


if __name__ == '__main__':
    main()

```
