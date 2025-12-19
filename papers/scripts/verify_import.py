#!/usr/bin/env python3
"""
Verify import completeness and integrity

Checks:
- Database count matches directory count
- All directories have required files (metadata.json, note.md)
- PDFs exist where specified
- No broken references

Usage:
    python verify_import.py --papers-dir ~/projects/papers --papers-db ~/projects/papers/papers.db
"""

import sqlite3
import json
import argparse
from pathlib import Path
from typing import List, Dict


def verify_import(papers_dir: Path, papers_db: Path) -> Dict:
    """
    Verify import completeness

    Returns dict with verification results
    """
    results = {
        'total_dirs': 0,
        'total_db_entries': 0,
        'missing_metadata': [],
        'missing_note': [],
        'missing_pdf': [],
        'broken_metadata': [],
        'valid_papers': 0
    }

    # Count database entries
    conn = sqlite3.connect(papers_db)
    cursor = conn.cursor()
    cursor.execute("SELECT COUNT(*) FROM papers")
    results['total_db_entries'] = cursor.fetchone()[0]

    # Get all papers from database
    cursor.execute("SELECT doi, title, pdf_path FROM papers")
    db_papers = {row[0]: {'title': row[1], 'pdf_path': row[2]} for row in cursor.fetchall()}
    conn.close()

    # Check directories
    dirs = [d for d in papers_dir.iterdir()
            if d.is_dir() and not d.name.startswith('.')
            and d.name not in ['scripts', 'references', 'templates']]

    results['total_dirs'] = len(dirs)

    for paper_dir in sorted(dirs):
        has_metadata = (paper_dir / 'metadata.json').exists()
        has_note = (paper_dir / 'note.md').exists()

        # Check metadata.json
        if not has_metadata:
            results['missing_metadata'].append(paper_dir.name)
        else:
            # Verify metadata is valid JSON
            try:
                with open(paper_dir / 'metadata.json') as f:
                    metadata = json.load(f)

                # Check PDF existence if specified
                if metadata.get('pdf_path'):
                    pdf_file = paper_dir / metadata['pdf_path']
                    if not pdf_file.exists():
                        results['missing_pdf'].append({
                            'folder': paper_dir.name,
                            'expected': metadata['pdf_path']
                        })

            except (json.JSONDecodeError, Exception) as e:
                results['broken_metadata'].append({
                    'folder': paper_dir.name,
                    'error': str(e)
                })

        # Check note.md
        if not has_note:
            results['missing_note'].append(paper_dir.name)

        # Count valid papers
        if has_metadata and has_note:
            results['valid_papers'] += 1

    return results


def print_results(results: Dict):
    """Print verification results"""
    print(f"\n{'='*70}")
    print("Verification Results")
    print(f"{'='*70}\n")

    print(f"Database entries: {results['total_db_entries']}")
    print(f"Paper directories: {results['total_dirs']}")
    print(f"Valid papers: {results['valid_papers']}")

    if results['total_db_entries'] != results['total_dirs']:
        print(f"\n⚠️  Mismatch: Database has {results['total_db_entries']} entries but "
              f"found {results['total_dirs']} directories")

    # Issues
    issues_found = (
        len(results['missing_metadata']) +
        len(results['missing_note']) +
        len(results['missing_pdf']) +
        len(results['broken_metadata'])
    )

    if issues_found == 0:
        print(f"\n✅ All papers verified successfully - no issues found!")
    else:
        print(f"\n⚠️  Found {issues_found} issue(s):\n")

        if results['missing_metadata']:
            print(f"Missing metadata.json ({len(results['missing_metadata'])}):")
            for folder in results['missing_metadata'][:10]:
                print(f"  - {folder}")
            if len(results['missing_metadata']) > 10:
                print(f"  ... and {len(results['missing_metadata']) - 10} more")

        if results['missing_note']:
            print(f"\nMissing note.md ({len(results['missing_note'])}):")
            for folder in results['missing_note'][:10]:
                print(f"  - {folder}")
            if len(results['missing_note']) > 10:
                print(f"  ... and {len(results['missing_note']) - 10} more")

        if results['missing_pdf']:
            print(f"\nMissing PDFs ({len(results['missing_pdf'])}):")
            for item in results['missing_pdf'][:10]:
                print(f"  - {item['folder']}: expected {item['expected']}")
            if len(results['missing_pdf']) > 10:
                print(f"  ... and {len(results['missing_pdf']) - 10} more")

        if results['broken_metadata']:
            print(f"\nBroken metadata.json ({len(results['broken_metadata'])}):")
            for item in results['broken_metadata'][:10]:
                print(f"  - {item['folder']}: {item['error']}")
            if len(results['broken_metadata']) > 10:
                print(f"  ... and {len(results['broken_metadata']) - 10} more")

    print()


def main():
    parser = argparse.ArgumentParser(description='Verify paper import integrity')
    parser.add_argument('--papers-dir', required=True, help='Path to papers directory')
    parser.add_argument('--papers-db', required=True, help='Path to papers.db')

    args = parser.parse_args()

    papers_dir = Path(args.papers_dir)
    papers_db = Path(args.papers_db)

    if not papers_dir.exists():
        print(f"Error: Papers directory not found: {papers_dir}")
        return 1

    if not papers_db.exists():
        print(f"Error: Database not found: {papers_db}")
        return 1

    results = verify_import(papers_dir, papers_db)
    print_results(results)

    # Exit code: 0 if no issues, 1 if issues found
    issues = (
        len(results['missing_metadata']) +
        len(results['missing_note']) +
        len(results['missing_pdf']) +
        len(results['broken_metadata'])
    )

    return 1 if issues > 0 else 0


if __name__ == '__main__':
    exit(main())
