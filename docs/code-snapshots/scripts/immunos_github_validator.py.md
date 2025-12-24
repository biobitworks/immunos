---
source: /Users/byron/projects/scripts/immunos_github_validator.py
relative: scripts/immunos_github_validator.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
IMMUNOS GitHub URL Validator
=============================

Validates GitHub repository URLs in documents to prevent broken links
and incorrect username references.

Features:
- Verifies GitHub repository existence via API
- Checks username accuracy
- Validates repository accessibility
- Integrates with IMMUNOS citation verification system

Usage:
    python3 immunos_github_validator.py document.md
    python3 immunos_github_validator.py --check-url github.com/user/repo
"""

import re
import sys
import requests
from typing import Dict, List, Tuple, Optional
from pathlib import Path
from dataclasses import dataclass


@dataclass
class GitHubValidation:
    """Result of GitHub URL validation"""
    url: str
    username: str
    repo: str
    exists: bool
    accessible: bool
    stars: Optional[int]
    description: Optional[str]
    error: Optional[str]


class GitHubValidator:
    """Validates GitHub repository URLs"""

    # Regex patterns for GitHub URLs
    GITHUB_PATTERNS = [
        r'github\.com/([^/\s]+)/([^/\s\)]+)',  # Standard format
        r'https?://github\.com/([^/\s]+)/([^/\s\)]+)',  # With protocol
        r'www\.github\.com/([^/\s]+)/([^/\s\)]+)',  # With www
    ]

    def __init__(self, github_token: Optional[str] = None):
        """
        Initialize validator

        Args:
            github_token: Optional GitHub API token for higher rate limits
        """
        self.github_token = github_token
        self.api_base = "https://api.github.com"

    def extract_github_urls(self, text: str) -> List[Tuple[str, str, str]]:
        """
        Extract all GitHub URLs from text

        Args:
            text: Text to search

        Returns:
            List of (full_url, username, repo) tuples
        """
        urls = []
        for pattern in self.GITHUB_PATTERNS:
            matches = re.finditer(pattern, text)
            for match in matches:
                username = match.group(1)
                repo = match.group(2)

                # Clean up repo name (remove trailing punctuation, fragments)
                repo = re.sub(r'[)\].,;:!?]+$', '', repo)
                repo = re.sub(r'#.*$', '', repo)  # Remove URL fragments

                full_url = f"github.com/{username}/{repo}"
                urls.append((full_url, username, repo))

        return list(set(urls))  # Remove duplicates

    def validate_repository(self, username: str, repo: str) -> GitHubValidation:
        """
        Validate a GitHub repository exists and is accessible

        Args:
            username: GitHub username
            repo: Repository name

        Returns:
            GitHubValidation object with results
        """
        url = f"github.com/{username}/{repo}"
        api_url = f"{self.api_base}/repos/{username}/{repo}"

        headers = {}
        if self.github_token:
            headers['Authorization'] = f'token {self.github_token}'

        try:
            response = requests.get(api_url, headers=headers, timeout=10)

            if response.status_code == 200:
                data = response.json()
                return GitHubValidation(
                    url=url,
                    username=username,
                    repo=repo,
                    exists=True,
                    accessible=True,
                    stars=data.get('stargazers_count'),
                    description=data.get('description'),
                    error=None
                )
            elif response.status_code == 404:
                return GitHubValidation(
                    url=url,
                    username=username,
                    repo=repo,
                    exists=False,
                    accessible=False,
                    stars=None,
                    description=None,
                    error="Repository not found (404)"
                )
            elif response.status_code == 403:
                return GitHubValidation(
                    url=url,
                    username=username,
                    repo=repo,
                    exists=None,  # Unknown
                    accessible=False,
                    stars=None,
                    description=None,
                    error=f"Rate limit exceeded or access denied (403)"
                )
            else:
                return GitHubValidation(
                    url=url,
                    username=username,
                    repo=repo,
                    exists=None,
                    accessible=False,
                    stars=None,
                    description=None,
                    error=f"HTTP {response.status_code}"
                )

        except requests.exceptions.RequestException as e:
            return GitHubValidation(
                url=url,
                username=username,
                repo=repo,
                exists=None,
                accessible=False,
                stars=None,
                description=None,
                error=f"Network error: {str(e)}"
            )

    def validate_document(self, filepath: Path) -> Dict[str, GitHubValidation]:
        """
        Validate all GitHub URLs in a document

        Args:
            filepath: Path to document

        Returns:
            Dict mapping URLs to validation results
        """
        with open(filepath, 'r', encoding='utf-8') as f:
            text = f.read()

        urls = self.extract_github_urls(text)
        results = {}

        for full_url, username, repo in urls:
            validation = self.validate_repository(username, repo)
            results[full_url] = validation

        return results

    def print_validation_report(self, results: Dict[str, GitHubValidation],
                               verbose: bool = False):
        """
        Print validation results

        Args:
            results: Dict of validation results
            verbose: Show all details including successful validations
        """
        print(f"\n{'='*70}")
        print(f"IMMUNOS GitHub URL Validation Report")
        print(f"{'='*70}\n")

        total = len(results)
        valid = sum(1 for v in results.values() if v.exists and v.accessible)
        invalid = sum(1 for v in results.values() if v.exists == False)
        unknown = sum(1 for v in results.values() if v.exists is None)

        print(f"Total URLs found: {total}")
        print(f"✅ Valid: {valid}")
        print(f"❌ Invalid: {invalid}")
        print(f"⚠️  Unknown: {unknown}\n")

        # Show problematic URLs
        if invalid > 0 or unknown > 0:
            print("Issues Found:")
            print("-" * 70)
            for url, validation in results.items():
                if not validation.exists or not validation.accessible:
                    status = "❌ NOT FOUND" if validation.exists == False else "⚠️  ERROR"
                    print(f"\n{status}: {url}")
                    print(f"  Username: {validation.username}")
                    print(f"  Repository: {validation.repo}")
                    if validation.error:
                        print(f"  Error: {validation.error}")

        # Show valid URLs if verbose
        if verbose and valid > 0:
            print("\n\nValid Repositories:")
            print("-" * 70)
            for url, validation in results.items():
                if validation.exists and validation.accessible:
                    print(f"\n✅ {url}")
                    if validation.description:
                        print(f"  Description: {validation.description}")
                    if validation.stars is not None:
                        print(f"  Stars: {validation.stars}")


def main():
    """CLI entry point"""
    import argparse

    parser = argparse.ArgumentParser(
        description='IMMUNOS GitHub URL Validator - Prevent broken repository links'
    )
    parser.add_argument('document', nargs='?',
                       help='Document to validate')
    parser.add_argument('--check-url',
                       help='Validate a single GitHub URL (format: user/repo)')
    parser.add_argument('--token',
                       help='GitHub API token for higher rate limits')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Show all validation details')

    args = parser.parse_args()

    validator = GitHubValidator(github_token=args.token)

    # Single URL check
    if args.check_url:
        # Parse user/repo format
        if '/' in args.check_url:
            parts = args.check_url.replace('github.com/', '').split('/')
            username, repo = parts[0], parts[1]
        else:
            print("❌ Error: URL must be in format 'username/repo' or 'github.com/username/repo'")
            sys.exit(1)

        print(f"\n🔍 Validating: github.com/{username}/{repo}")
        validation = validator.validate_repository(username, repo)

        if validation.exists and validation.accessible:
            print(f"\n✅ Repository exists and is accessible")
            if validation.description:
                print(f"Description: {validation.description}")
            if validation.stars is not None:
                print(f"Stars: {validation.stars}")
        elif validation.exists == False:
            print(f"\n❌ Repository not found")
        else:
            print(f"\n⚠️  Could not verify: {validation.error}")

        sys.exit(0)

    # Document validation
    if args.document:
        filepath = Path(args.document)
        if not filepath.exists():
            print(f"❌ Error: File not found: {filepath}")
            sys.exit(1)

        print(f"\n🔍 Scanning document: {filepath}")
        results = validator.validate_document(filepath)

        if not results:
            print("\n✅ No GitHub URLs found in document")
            sys.exit(0)

        validator.print_validation_report(results, verbose=args.verbose)

        # Exit code based on validation
        invalid_count = sum(1 for v in results.values() if v.exists == False)
        sys.exit(1 if invalid_count > 0 else 0)

    # No arguments
    parser.print_help()
    sys.exit(1)


if __name__ == '__main__':
    main()

```
