#!/usr/bin/env python3
"""
IMMUNOS Inbox Intake System

Processes files from the inbox directory and routes them to appropriate destinations.

Supports:
- PDFs → /papers/{doi-or-name}/
- Images (png/jpg/svg) → appropriate assets folders
- Datasets (csv/json/parquet) → /data/ or /research/datasets/
- Citations (bib/ris) → /research/citations/

Usage:
    # Preview without moving files
    python3 immunos_inbox.py --dry-run

    # Process inbox with verbose output
    python3 immunos_inbox.py --verbose

    # Process specific directory
    python3 immunos_inbox.py --inbox /path/to/inbox

    # Skip journal logging
    python3 immunos_inbox.py --no-journal

Examples:
    # Daily intake workflow
    python3 immunos_inbox.py --verbose

    # Preview what would happen
    python3 immunos_inbox.py --dry-run --verbose

    # Process and skip daily journal
    python3 immunos_inbox.py --no-journal

File Type Detection:
    - PDFs: Extracts DOI from filename pattern or content
    - Images: Uses filename hints to determine destination
    - Datasets: Routes based on content and filename
    - Citations: Moves to research/citations/

DOI Extraction:
    - Filename patterns: name_DOI_10.1234-5678.pdf
    - PDF content: Basic text search for DOI patterns
    - Fallback: Uses filename without extension

Logging:
    - Daily journal: /daily/YYYY-MM-DD.md
    - Publications log: /docs/reference/publications-log.md
    - Console: Optional verbose output
"""

import os
import re
import sys
import json
import argparse
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Tuple


class InboxProcessor:
    """Process files from inbox and route to appropriate destinations"""

    # File type mappings
    PDF_EXTENSIONS = {'.pdf'}
    IMAGE_EXTENSIONS = {'.png', '.jpg', '.jpeg', '.svg', '.gif', '.webp'}
    DATASET_EXTENSIONS = {'.csv', '.json', '.jsonl', '.parquet', '.tsv', '.xlsx'}
    CITATION_EXTENSIONS = {'.bib', '.ris', '.enw', '.nbib'}

    # Files to ignore
    IGNORE_FILES = {'.DS_Store', 'README.md', '.gitkeep', 'Thumbs.db'}

    # DOI regex patterns
    # Standard DOI format: 10.1234/suffix
    DOI_PATTERN = re.compile(r'10\.\d{4,9}/[-._;()/:A-Za-z0-9]+')
    # Filename-friendly DOI: 10.1234-suffix (slash replaced with hyphen/underscore)
    DOI_FILENAME_PATTERN = re.compile(r'10\.\d{4,9}[-_][-._;()A-Za-z0-9]+')

    def __init__(self,
                 inbox_path: Path,
                 projects_root: Path,
                 dry_run: bool = False,
                 verbose: bool = False,
                 no_journal: bool = False):
        self.inbox_path = inbox_path
        self.projects_root = projects_root
        self.dry_run = dry_run
        self.verbose = verbose
        self.no_journal = no_journal

        # Define destination paths
        self.papers_dir = projects_root / "papers"
        self.data_dir = projects_root / "data"
        self.research_dir = projects_root / "research"
        self.assets_dir = projects_root / "assets"
        self.daily_dir = projects_root / "daily"
        self.docs_dir = projects_root / "docs" / "reference"

        # Processing stats
        self.stats = {
            'processed': 0,
            'skipped': 0,
            'errors': 0,
            'pdfs': 0,
            'images': 0,
            'datasets': 0,
            'citations': 0,
            'other': 0
        }

        # Processed items for logging
        self.processed_items: List[Dict] = []

    def log(self, message: str, level: str = 'info'):
        """Log message if verbose mode enabled"""
        if self.verbose:
            prefix = {
                'info': '[INFO]',
                'warn': '[WARN]',
                'error': '[ERROR]',
                'success': '[OK]'
            }.get(level, '[INFO]')
            print(f"{prefix} {message}")

    def extract_doi_from_filename(self, filename: str) -> Optional[str]:
        """
        Extract DOI from filename patterns.

        Handles both standard DOIs (10.1234/suffix) and filename-friendly
        versions where slash is replaced with hyphen or underscore.
        """
        # Remove extension for cleaner matching
        name_no_ext = filename.rsplit('.', 1)[0]

        # Try standard DOI pattern first (with slash)
        doi_match = self.DOI_PATTERN.search(name_no_ext)
        if doi_match:
            doi = doi_match.group(0)
            doi = doi.rstrip('_').rstrip('-')
            return doi

        # Try filename-friendly pattern (slash replaced with hyphen/underscore)
        doi_match = self.DOI_FILENAME_PATTERN.search(name_no_ext)
        if doi_match:
            doi = doi_match.group(0)
            doi = doi.rstrip('_').rstrip('-')
            # Convert back to standard DOI format (first hyphen/underscore to slash)
            # Find the first hyphen or underscore after the registrant code
            parts = doi.split('-', 1) if '-' in doi else doi.split('_', 1)
            if len(parts) == 2:
                doi = f"{parts[0]}/{parts[1]}"
            return doi

        return None

    def extract_doi_from_pdf(self, pdf_path: Path) -> Optional[str]:
        """
        Extract DOI from PDF content using basic text extraction.

        Note: Uses simple binary search for DOI patterns without heavy dependencies.
        For production use, consider PyPDF2 or pdfplumber for better extraction.
        """
        try:
            # Read first 50KB of PDF (usually contains metadata/first page)
            with open(pdf_path, 'rb') as f:
                content = f.read(50000)

            # Decode as latin-1 (works for most PDFs)
            try:
                text = content.decode('latin-1', errors='ignore')
            except Exception:
                return None

            # Search for DOI pattern
            doi_match = self.DOI_PATTERN.search(text)
            if doi_match:
                return doi_match.group(0)

        except Exception as e:
            self.log(f"Failed to extract DOI from {pdf_path.name}: {e}", 'warn')

        return None

    def sanitize_doi_for_path(self, doi: str) -> str:
        """Convert DOI to valid directory name"""
        # Replace invalid path characters
        sanitized = doi.replace('/', '_')
        sanitized = re.sub(r'[<>:"|?*]', '_', sanitized)
        return sanitized

    def get_pdf_destination(self, pdf_path: Path) -> Tuple[Path, Optional[str]]:
        """
        Determine destination for PDF file.

        Returns:
            (destination_dir, doi_or_none)
        """
        # Try to extract DOI
        doi = self.extract_doi_from_filename(pdf_path.name)

        if not doi:
            doi = self.extract_doi_from_pdf(pdf_path)

        if doi:
            # Create directory using DOI
            dir_name = self.sanitize_doi_for_path(doi)
            dest_dir = self.papers_dir / dir_name
        else:
            # Use filename without extension
            dir_name = pdf_path.stem
            dest_dir = self.papers_dir / dir_name

        return dest_dir, doi

    def get_image_destination(self, image_path: Path) -> Path:
        """
        Determine destination for image file based on filename hints.

        Patterns:
            - *figure*, *fig* → /assets/figures/
            - *diagram*, *chart* → /assets/diagrams/
            - *screenshot*, *screen* → /assets/screenshots/
            - *logo*, *icon* → /assets/logos/
            - default → /assets/images/
        """
        name_lower = image_path.name.lower()

        if 'figure' in name_lower or 'fig' in name_lower:
            return self.assets_dir / "figures"
        elif 'diagram' in name_lower or 'chart' in name_lower:
            return self.assets_dir / "diagrams"
        elif 'screenshot' in name_lower or 'screen' in name_lower:
            return self.assets_dir / "screenshots"
        elif 'logo' in name_lower or 'icon' in name_lower:
            return self.assets_dir / "logos"
        else:
            return self.assets_dir / "images"

    def get_dataset_destination(self, dataset_path: Path) -> Path:
        """
        Determine destination for dataset file.

        Patterns:
            - *research*, *experiment* → /research/datasets/
            - *raw*, *processed* → /data/
            - default → /data/
        """
        name_lower = dataset_path.name.lower()

        if 'research' in name_lower or 'experiment' in name_lower:
            return self.research_dir / "datasets"
        else:
            return self.data_dir

    def process_pdf(self, pdf_path: Path) -> bool:
        """Process PDF file and move to papers directory"""
        try:
            dest_dir, doi = self.get_pdf_destination(pdf_path)
            dest_file = dest_dir / pdf_path.name

            if not self.dry_run:
                dest_dir.mkdir(parents=True, exist_ok=True)
                pdf_path.rename(dest_file)

            self.stats['pdfs'] += 1
            self.stats['processed'] += 1

            self.processed_items.append({
                'type': 'pdf',
                'filename': pdf_path.name,
                'source': str(pdf_path),
                'destination': str(dest_file),
                'doi': doi,
                'timestamp': datetime.now().isoformat()
            })

            action = "Would move" if self.dry_run else "Moved"
            self.log(f"{action} PDF: {pdf_path.name} → {dest_file}", 'success')
            if doi:
                self.log(f"  DOI: {doi}", 'info')

            return True

        except Exception as e:
            self.stats['errors'] += 1
            self.log(f"Failed to process PDF {pdf_path.name}: {e}", 'error')
            return False

    def process_image(self, image_path: Path) -> bool:
        """Process image file and move to appropriate assets folder"""
        try:
            dest_dir = self.get_image_destination(image_path)
            dest_file = dest_dir / image_path.name

            if not self.dry_run:
                dest_dir.mkdir(parents=True, exist_ok=True)
                image_path.rename(dest_file)

            self.stats['images'] += 1
            self.stats['processed'] += 1

            self.processed_items.append({
                'type': 'image',
                'filename': image_path.name,
                'source': str(image_path),
                'destination': str(dest_file),
                'timestamp': datetime.now().isoformat()
            })

            action = "Would move" if self.dry_run else "Moved"
            self.log(f"{action} Image: {image_path.name} → {dest_file}", 'success')

            return True

        except Exception as e:
            self.stats['errors'] += 1
            self.log(f"Failed to process image {image_path.name}: {e}", 'error')
            return False

    def process_dataset(self, dataset_path: Path) -> bool:
        """Process dataset file and move to data directory"""
        try:
            dest_dir = self.get_dataset_destination(dataset_path)
            dest_file = dest_dir / dataset_path.name

            if not self.dry_run:
                dest_dir.mkdir(parents=True, exist_ok=True)
                dataset_path.rename(dest_file)

            self.stats['datasets'] += 1
            self.stats['processed'] += 1

            self.processed_items.append({
                'type': 'dataset',
                'filename': dataset_path.name,
                'source': str(dataset_path),
                'destination': str(dest_file),
                'timestamp': datetime.now().isoformat()
            })

            action = "Would move" if self.dry_run else "Moved"
            self.log(f"{action} Dataset: {dataset_path.name} → {dest_file}", 'success')

            return True

        except Exception as e:
            self.stats['errors'] += 1
            self.log(f"Failed to process dataset {dataset_path.name}: {e}", 'error')
            return False

    def process_citation(self, citation_path: Path) -> bool:
        """Process citation file and move to research/citations"""
        try:
            dest_dir = self.research_dir / "citations"
            dest_file = dest_dir / citation_path.name

            if not self.dry_run:
                dest_dir.mkdir(parents=True, exist_ok=True)
                citation_path.rename(dest_file)

            self.stats['citations'] += 1
            self.stats['processed'] += 1

            self.processed_items.append({
                'type': 'citation',
                'filename': citation_path.name,
                'source': str(citation_path),
                'destination': str(dest_file),
                'timestamp': datetime.now().isoformat()
            })

            action = "Would move" if self.dry_run else "Moved"
            self.log(f"{action} Citation: {citation_path.name} → {dest_file}", 'success')

            return True

        except Exception as e:
            self.stats['errors'] += 1
            self.log(f"Failed to process citation {citation_path.name}: {e}", 'error')
            return False

    def process_file(self, file_path: Path) -> bool:
        """Route file to appropriate processor based on extension"""
        ext = file_path.suffix.lower()

        if ext in self.PDF_EXTENSIONS:
            return self.process_pdf(file_path)
        elif ext in self.IMAGE_EXTENSIONS:
            return self.process_image(file_path)
        elif ext in self.DATASET_EXTENSIONS:
            return self.process_dataset(file_path)
        elif ext in self.CITATION_EXTENSIONS:
            return self.process_citation(file_path)
        else:
            self.stats['other'] += 1
            self.stats['skipped'] += 1
            self.log(f"Skipped unknown file type: {file_path.name} ({ext})", 'warn')
            return False

    def log_to_daily_journal(self):
        """Append processed items to today's daily journal"""
        if self.no_journal or self.dry_run:
            return

        if not self.processed_items:
            return

        try:
            today = datetime.now().strftime('%Y-%m-%d')
            journal_file = self.daily_dir / f"{today}.md"

            # Create journal file if it doesn't exist
            if not journal_file.exists():
                journal_file.parent.mkdir(parents=True, exist_ok=True)
                journal_file.write_text(f"# {today}\n\n")

            # Append inbox intake entries
            with open(journal_file, 'a') as f:
                f.write(f"\n## Inbox Intake ({datetime.now().strftime('%H:%M')})\n\n")

                for item in self.processed_items:
                    if item['type'] == 'pdf':
                        doi_info = f" (DOI: {item['doi']})" if item['doi'] else ""
                        f.write(f"- PDF: `{item['filename']}`{doi_info} → `{item['destination']}`\n")
                    else:
                        f.write(f"- {item['type'].upper()}: `{item['filename']}` → `{item['destination']}`\n")

            self.log(f"Updated daily journal: {journal_file}", 'success')

        except Exception as e:
            self.log(f"Failed to update daily journal: {e}", 'error')

    def log_to_publications(self):
        """Append PDFs with DOIs to publications log"""
        if self.dry_run:
            return

        pdf_items = [item for item in self.processed_items
                     if item['type'] == 'pdf' and item['doi']]

        if not pdf_items:
            return

        try:
            pub_log = self.docs_dir / "publications-log.md"

            if not pub_log.exists():
                self.log(f"Publications log not found: {pub_log}", 'warn')
                return

            # Read existing content
            content = pub_log.read_text()

            # Find or create today's section
            today = datetime.now().strftime('%Y-%m-%d')
            section_marker = f"## {today}"

            if section_marker not in content:
                # Find the first ## section and insert before it
                lines = content.split('\n')
                insert_idx = 0
                for i, line in enumerate(lines):
                    if line.startswith('## '):
                        insert_idx = i
                        break

                # Insert new section
                new_section = f"\n{section_marker}\n"
                lines.insert(insert_idx, new_section)
                content = '\n'.join(lines)

            # Append entries to today's section
            with open(pub_log, 'w') as f:
                f.write(content)
                if not content.endswith('\n'):
                    f.write('\n')

                for item in pdf_items:
                    entry = (
                        f"- {item['filename']}\n"
                        f"  DOI: https://doi.org/{item['doi']}\n"
                        f"  Path: {item['destination']}\n"
                        f"  Intake: {item['timestamp']}\n"
                    )
                    f.write(entry)

            self.log(f"Updated publications log: {pub_log}", 'success')

        except Exception as e:
            self.log(f"Failed to update publications log: {e}", 'error')

    def process_inbox(self) -> Dict:
        """Process all files in inbox directory"""
        if not self.inbox_path.exists():
            self.log(f"Inbox directory not found: {self.inbox_path}", 'error')
            return self.stats

        if not self.inbox_path.is_dir():
            self.log(f"Inbox path is not a directory: {self.inbox_path}", 'error')
            return self.stats

        self.log(f"Processing inbox: {self.inbox_path}", 'info')
        if self.dry_run:
            self.log("DRY RUN MODE - No files will be moved", 'warn')

        # Get all files in inbox
        files = [f for f in self.inbox_path.iterdir()
                if f.is_file() and f.name not in self.IGNORE_FILES]

        self.log(f"Found {len(files)} files to process", 'info')

        # Process each file
        for file_path in files:
            self.process_file(file_path)

        # Update logs
        if not self.dry_run:
            self.log_to_daily_journal()
            self.log_to_publications()

        return self.stats

    def print_summary(self):
        """Print processing summary"""
        print("\n" + "="*60)
        print("INBOX PROCESSING SUMMARY")
        print("="*60)
        print(f"Total processed: {self.stats['processed']}")
        print(f"  PDFs:          {self.stats['pdfs']}")
        print(f"  Images:        {self.stats['images']}")
        print(f"  Datasets:      {self.stats['datasets']}")
        print(f"  Citations:     {self.stats['citations']}")
        print(f"  Other:         {self.stats['other']}")
        print(f"Skipped:         {self.stats['skipped']}")
        print(f"Errors:          {self.stats['errors']}")
        print("="*60)

        if self.dry_run:
            print("\nDRY RUN COMPLETE - No files were moved")
        else:
            print(f"\nProcessed {self.stats['processed']} files successfully")


def main():
    parser = argparse.ArgumentParser(
        description="IMMUNOS Inbox Intake System - Process and route inbox files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --dry-run --verbose     Preview what would happen
  %(prog)s                         Process inbox files
  %(prog)s --no-journal            Skip daily journal logging

File Types:
  PDFs:      → /papers/{doi-or-name}/
  Images:    → /assets/{figures,diagrams,screenshots,etc}/
  Datasets:  → /data/ or /research/datasets/
  Citations: → /research/citations/
        """
    )

    parser.add_argument(
        '--inbox',
        type=Path,
        default=Path('/Users/byron/projects/inbox'),
        help='Inbox directory path (default: /Users/byron/projects/inbox)'
    )

    parser.add_argument(
        '--projects',
        type=Path,
        default=Path('/Users/byron/projects'),
        help='Projects root directory (default: /Users/byron/projects)'
    )

    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Preview what would happen without moving files'
    )

    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose output'
    )

    parser.add_argument(
        '--no-journal',
        action='store_true',
        help='Skip daily journal logging'
    )

    args = parser.parse_args()

    # Create processor
    processor = InboxProcessor(
        inbox_path=args.inbox,
        projects_root=args.projects,
        dry_run=args.dry_run,
        verbose=args.verbose,
        no_journal=args.no_journal
    )

    # Process inbox
    stats = processor.process_inbox()

    # Print summary
    processor.print_summary()

    # Exit with error code if there were errors
    sys.exit(1 if stats['errors'] > 0 else 0)


if __name__ == '__main__':
    main()
