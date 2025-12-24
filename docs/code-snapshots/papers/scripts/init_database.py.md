---
source: /Users/byron/projects/papers/scripts/init_database.py
relative: papers/scripts/init_database.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
Initialize papers.db SQLite database

Creates tables for tracking downloaded papers and figures.
"""

import sqlite3
from pathlib import Path

def init_database(db_path="papers.db"):
    """Initialize SQLite database with schema"""

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Create papers table
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS papers (
            doi TEXT PRIMARY KEY,
            title TEXT NOT NULL,
            authors TEXT,
            journal TEXT,
            year INTEGER,
            url TEXT,
            pdf_path TEXT,
            download_date TEXT,
            note_path TEXT,
            tags TEXT
        )
    """)

    # Create figures table
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS figures (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            paper_doi TEXT,
            figure_num INTEGER,
            filename TEXT,
            file_path TEXT,
            FOREIGN KEY (paper_doi) REFERENCES papers(doi)
        )
    """)

    # Create index for faster lookups
    cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_papers_year ON papers(year)
    """)

    cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_papers_journal ON papers(journal)
    """)

    conn.commit()
    conn.close()

    print(f"✓ Database initialized: {db_path}")

if __name__ == "__main__":
    import sys
    db_path = sys.argv[1] if len(sys.argv) > 1 else "../papers.db"
    init_database(db_path)

```
