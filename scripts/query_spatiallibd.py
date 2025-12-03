#!/usr/bin/env python3
"""
Query spatialLIBD for Aging Biology Questions

Python interface to query the spatialLIBD brain spatial transcriptomics
dataset for aging-related genes WITHOUT downloading the full 2GB dataset.

Usage:
    python scripts/query_spatiallibd.py
    or import: from scripts.query_spatiallibd import *
"""

import webbrowser
import pandas as pd
from typing import List, Dict, Optional
from pathlib import Path

# Project paths
PROJECT_ROOT = Path(__file__).parent.parent
GENAGE_PATH = PROJECT_ROOT / "data" / "genage" / "human" / "genage_human.csv"
CELLAGE_PATH = PROJECT_ROOT / "data" / "cellage" / "cellage3.tsv"


class SpatialLIBDQuery:
    """Query system for spatialLIBD spatial transcriptomics data"""

    BASE_URL = "http://spatial.libd.org/spatialLIBD/"

    def __init__(self):
        """Initialize query system"""
        self.aging_genes = self._load_aging_genes()
        self.senescence_genes = self._load_senescence_genes()

    def _load_aging_genes(self) -> List[str]:
        """Load aging genes from GenAge"""
        try:
            genage = pd.read_csv(GENAGE_PATH)
            return genage['symbol'].tolist()
        except FileNotFoundError:
            print("⚠️  GenAge file not found, using example genes")
            return ["TP53", "FOXO3", "IGF1R", "SIRT1", "APOE", "TERT"]

    def _load_senescence_genes(self) -> List[str]:
        """Load senescence genes from CellAge"""
        try:
            cellage = pd.read_csv(CELLAGE_PATH, sep='\t')
            return cellage['Gene symbol'].tolist()
        except FileNotFoundError:
            print("⚠️  CellAge file not found, using example genes")
            return ["CDKN2A", "TP53", "RB1", "CDKN1A", "BCL6"]

    def query_genes(self, genes: List[str], open_browser: bool = True) -> Dict:
        """
        Query genes in spatialLIBD web interface

        Args:
            genes: List of gene symbols to query
            open_browser: If True, open URLs in web browser

        Returns:
            Dictionary with query results and URLs
        """
        print("\n" + "="*70)
        print("spatialLIBD Aging Biology Query")
        print("="*70 + "\n")

        print(f"Querying {len(genes)} genes in brain spatial transcriptomics data\n")
        print("Method: Web interface (no download required) ✓\n")

        urls = []
        for gene in genes:
            url = f"{self.BASE_URL}?gene={gene}"
            urls.append((gene, url))
            print(f"  • {gene:10s} → {url}")

        print("\n" + "-"*70)
        print("What you can explore in web interface:")
        print("  - Spatial expression patterns across tissue")
        print("  - Layer-specific enrichment (L1-L6, WM)")
        print("  - Expression in different subjects")
        print("  - Export high-resolution figures")
        print("-"*70 + "\n")

        if open_browser:
            print("Opening first gene in browser...")
            webbrowser.open(urls[0][1])
            print("  (Open other URLs manually to avoid overwhelming browser)\n")

        return {
            'genes': genes,
            'urls': urls,
            'base_url': self.BASE_URL
        }

    def ask_question(self, question_id: str, open_browser: bool = True) -> Dict:
        """
        Ask predefined aging biology question

        Args:
            question_id: Question ID (q1, q2, q3, q4, q5)
            open_browser: Open results in browser

        Returns:
            Query results dictionary
        """
        questions = self.get_aging_questions()

        if question_id not in questions:
            print(f"❌ Invalid question_id: {question_id}")
            print(f"Available: {', '.join(questions.keys())}\n")
            self.list_questions()
            return {}

        q = questions[question_id]

        print("\n" + "="*70)
        print(f"QUESTION {question_id}")
        print("="*70 + "\n")
        print(f"{q['question']}\n")
        print(f"Description: {q['description']}")
        print(f"Genes: {', '.join(q['genes'])}")
        print("-"*70 + "\n")

        return self.query_genes(q['genes'], open_browser=open_browser)

    def get_aging_questions(self) -> Dict:
        """Get predefined aging biology questions"""
        return {
            'q1': {
                'question': 'Are aging genes layer-specific in the brain?',
                'genes': ['TP53', 'FOXO3', 'IGF1R', 'SIRT1'],
                'description': 'Check if GenAge aging-related genes show preferential expression in specific cortical layers'
            },
            'q2': {
                'question': 'Where are cellular senescence genes expressed in DLPFC?',
                'genes': ['CDKN2A', 'TP53', 'RB1', 'CDKN1A', 'BCL6'],
                'description': 'Map CellAge senescence markers to spatial locations in prefrontal cortex'
            },
            'q3': {
                'question': 'Do longevity genes show distinct spatial patterns?',
                'genes': ['FOXO3', 'APOE', 'SIRT1', 'TERT'],
                'description': 'Examine spatial distribution of pro-longevity genes'
            },
            'q4': {
                'question': 'Are immune/inflammation genes enriched in specific layers?',
                'genes': ['IL6', 'TNF', 'NFKB1', 'CD68'],
                'description': 'Spatial pattern of aging-associated inflammation markers'
            },
            'q5': {
                'question': 'How is BCL6 (senescence inhibitor) distributed in brain?',
                'genes': ['BCL6'],
                'description': 'Following up on Shvarts 2002 paper about BCL6 and senescence'
            }
        }

    def list_questions(self):
        """List all available questions"""
        print("\n" + "="*70)
        print("Available Aging Biology Questions for spatialLIBD")
        print("="*70 + "\n")

        for qid, q in self.get_aging_questions().items():
            print(f"[{qid}] {q['question']}")
            print(f"     Genes: {', '.join(q['genes'])}")
            print(f"     {q['description']}\n")

        print("Usage:")
        print('  query.ask_question("q1")  # Layer-specific aging genes')
        print('  query.ask_question("q5")  # BCL6 senescence inhibitor')
        print('  query.query_genes(["TP53", "CDKN2A"])  # Custom query\n')

    def query_genage_top(self, n: int = 10, open_browser: bool = False) -> Dict:
        """Query top N genes from GenAge"""
        top_genes = self.aging_genes[:n]
        print(f"\nQuerying top {n} GenAge aging genes in brain spatial data\n")
        return self.query_genes(top_genes, open_browser=open_browser)

    def query_cellage_sample(self, n: int = 10, open_browser: bool = False) -> Dict:
        """Query sample of CellAge senescence genes"""
        sample_genes = self.senescence_genes[:n]
        print(f"\nQuerying {n} CellAge senescence genes in brain spatial data\n")
        return self.query_genes(sample_genes, open_browser=open_browser)


def main():
    """Main execution when run as script"""
    print("\n")
    print("╔" + "═"*68 + "╗")
    print("║" + " "*20 + "spatialLIBD Aging Biology Query System" + " "*9 + "║")
    print("║" + " "*10 + "Query brain spatial transcriptomics for aging genes" + " "*7 + "║")
    print("╚" + "═"*68 + "╝")
    print("\n")

    # Initialize query system
    query = SpatialLIBDQuery()

    print(f"✓ Loaded {len(query.aging_genes)} aging genes from GenAge")
    print(f"✓ Loaded {len(query.senescence_genes)} senescence genes from CellAge\n")

    # Show available questions
    query.list_questions()

    print("\nExamples:")
    print("  from scripts.query_spatiallibd import SpatialLIBDQuery")
    print("  query = SpatialLIBDQuery()")
    print('  query.ask_question("q1")  # Layer-specific aging genes')
    print('  query.ask_question("q5")  # BCL6 senescence inhibitor')
    print('  query.query_genes(["TP53", "CDKN2A"])  # Custom query\n')

    print("💡 TIP: All queries use web interface (no download required)")
    print("   No 2GB dataset download needed!\n")


if __name__ == "__main__":
    main()
