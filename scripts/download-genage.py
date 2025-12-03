#!/usr/bin/env python3
"""
Download GenAge Datasets

Downloads human genes, model organism genes, and gene expression signatures
from the Human Ageing Genomic Resources (HAGR) database.

Usage:
    python scripts/download-genage.py
"""

import json
import urllib.request
import zipfile
from datetime import datetime
from pathlib import Path


def download_file(url: str, dest_path: Path) -> bool:
    """Download a file from URL to destination path."""
    print(f"Downloading: {url}")
    print(f"        to: {dest_path}")

    try:
        with urllib.request.urlopen(url) as response:
            content = response.read()
            dest_path.parent.mkdir(parents=True, exist_ok=True)
            with open(dest_path, 'wb') as f:
                f.write(content)
        print(f"✓ Downloaded: {len(content):,} bytes")
        return True
    except Exception as e:
        print(f"✗ Error: {e}")
        return False


def extract_zip(zip_path: Path, extract_dir: Path) -> list[Path]:
    """Extract zip file and return list of extracted files."""
    print(f"Extracting: {zip_path.name}")

    extracted_files = []
    try:
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(extract_dir)
            extracted_files = [extract_dir / name for name in zip_ref.namelist()]
        print(f"✓ Extracted {len(extracted_files)} file(s)")
        for f in extracted_files:
            if f.is_file():
                print(f"  - {f.name} ({f.stat().st_size:,} bytes)")
        return extracted_files
    except Exception as e:
        print(f"✗ Error extracting: {e}")
        return []


def create_metadata(dataset_name: str, url: str, files: list[Path]) -> dict:
    """Create metadata dictionary for dataset."""
    return {
        "dataset": dataset_name,
        "source": "Human Ageing Genomic Resources (HAGR)",
        "source_url": "https://genomics.senescence.info/",
        "download_url": url,
        "downloaded_date": datetime.now().isoformat(),
        "build": "21",
        "build_date": "2023-08-28",
        "license": "Creative Commons Attribution 3.0 Unported License",
        "files": [
            {
                "name": f.name,
                "path": str(f.relative_to(Path.cwd())),
                "size_bytes": f.stat().st_size if f.is_file() else 0,
                "format": f.suffix[1:] if f.suffix else "unknown"
            }
            for f in files
        ]
    }


def save_metadata(metadata: dict, dest_path: Path):
    """Save metadata as JSON file."""
    with open(dest_path, 'w') as f:
        json.dump(metadata, f, indent=2)
    print(f"✓ Metadata saved: {dest_path.name}")


def main():
    """Download all GenAge datasets."""
    print("=" * 70)
    print("GenAge Dataset Download")
    print("=" * 70)
    print()

    # Base directory
    data_dir = Path("/Users/byron/projects/data/genage")

    # Dataset URLs (these are the expected URLs - may need verification)
    datasets = {
        "human": {
            "name": "GenAge Human Genes",
            "urls": [
                "https://genomics.senescence.info/genes/human_genes.zip",
                "https://genomics.senescence.info/genes/genage_human.zip",
                "http://genomics.senescence.info/genes/human_genes.zip",
            ],
            "dest_dir": data_dir / "human"
        },
        "models": {
            "name": "GenAge Model Organisms",
            "urls": [
                "https://genomics.senescence.info/genes/models_genes.zip",
                "https://genomics.senescence.info/genes/genage_models.zip",
                "http://genomics.senescence.info/genes/models_genes.zip",
            ],
            "dest_dir": data_dir / "models"
        },
        "expression": {
            "name": "Gene Expression Signatures",
            "urls": [
                "https://genomics.senescence.info/gene_expression/signatures_supplement.zip",
                "http://genomics.senescence.info/gene_expression/signatures_supplement.zip",
            ],
            "dest_dir": data_dir / "expression"
        }
    }

    success_count = 0
    total_count = len(datasets)

    for key, info in datasets.items():
        print(f"\n{'='*70}")
        print(f"Dataset: {info['name']}")
        print(f"{'='*70}\n")

        dest_dir = info['dest_dir']
        dest_dir.mkdir(parents=True, exist_ok=True)

        # Try each URL until one works
        downloaded = False
        zip_path = None
        successful_url = None

        for url in info['urls']:
            zip_path = dest_dir / Path(url).name
            if download_file(url, zip_path):
                downloaded = True
                successful_url = url
                break

        if not downloaded:
            print(f"✗ Failed to download {info['name']} from any URL")
            print(f"  Please manually download from:")
            print(f"  https://genomics.senescence.info/download.html")
            print(f"  and place in: {dest_dir}/")
            continue

        # Extract zip file
        extracted_files = extract_zip(zip_path, dest_dir)
        if not extracted_files:
            print(f"✗ Failed to extract {info['name']}")
            continue

        # Create metadata
        all_files = [zip_path] + extracted_files
        metadata = create_metadata(info['name'], successful_url, all_files)
        metadata_path = dest_dir / "metadata.json"
        save_metadata(metadata, metadata_path)

        success_count += 1
        print(f"✓ {info['name']} complete")

    # Summary
    print(f"\n{'='*70}")
    print("Download Summary")
    print(f"{'='*70}")
    print(f"Successful: {success_count}/{total_count} datasets")
    print(f"Location: {data_dir}")
    print(f"\nNext steps:")
    print(f"1. Review files in {data_dir}")
    print(f"2. Run analysis: python scripts/analyze-genage.py")
    print(f"3. Check Obsidian documentation in /research/datasets/")
    print()


if __name__ == "__main__":
    main()
