---
source: /Users/byron/projects/prion-clock/citations/verification/deepseek_verify.py
relative: prion-clock/citations/verification/deepseek_verify.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
Use DeepSeek-R1 to verify citations and detect hallucinations
Negative reinforcement: flag anything that seems suspicious
"""
import json
import requests
from pathlib import Path

OLLAMA_URL = "http://localhost:11434/api/generate"
MODEL = "deepseek-r1:14b"

def verify_with_deepseek(citation):
    """
    Ask DeepSeek to verify citation and flag issues
    Returns: {suspicious: bool, issues: [], confidence: float}
    """
    prompt = f"""You are a citation verification system using negative reinforcement.
Your job is to be SKEPTICAL and flag ANY issues with this citation.

Citation {citation['number']}:
{citation['full_text']}

Extracted metadata:
- DOI: {citation.get('doi', 'NOT FOUND')}
- PMID: {citation.get('pmid', 'NOT FOUND')}

Check for:
1. Multiple papers by same first author (could be wrong paper)
2. Missing or malformed DOI/PMID
3. Journal name inconsistencies
4. Year mismatches
5. Author list truncation issues
6. Any signs this could be the WRONG paper

Respond in JSON format ONLY, with no thinking or explanation before or after:
{{
  "suspicious": true/false,
  "issues": ["list of specific issues"],
  "confidence": 0.0-1.0,
  "recommendation": "verify/accept/investigate"
}}
"""

    try:
        response = requests.post(OLLAMA_URL, json={
            'model': MODEL,
            'prompt': prompt,
            'stream': False
        }, timeout=60)

        result = response.json()
        response_text = result.get('response', '{}')

        # Extract JSON from response (DeepSeek might include thinking)
        # Look for JSON block
        import re
        json_match = re.search(r'\{[\s\S]*\}', response_text)
        if json_match:
            return json.loads(json_match.group())
        else:
            # Default to flagging if we can't parse
            return {
                "suspicious": True,
                "issues": ["Could not parse DeepSeek response"],
                "confidence": 0.0,
                "recommendation": "verify"
            }
    except Exception as e:
        print(f"  ⚠️  Error calling DeepSeek: {e}")
        return {
            "suspicious": True,
            "issues": [f"DeepSeek error: {str(e)}"],
            "confidence": 0.0,
            "recommendation": "verify"
        }

def main():
    citations_path = Path('prion-clock/citations/extraction/full-citations.json')
    output_path = Path('prion-clock/citations/verification/deepseek-flags.json')

    with open(citations_path) as f:
        citations = json.load(f)

    flagged = []

    for citation in citations:
        print(f"Verifying citation {citation['number']}...")
        result = verify_with_deepseek(citation)

        if result.get('suspicious', True):
            flagged.append({
                'citation': citation,
                'verification': result
            })
            issues_str = ', '.join(result.get('issues', []))
            print(f"  ⚠️  FLAGGED: {issues_str}")
        else:
            print(f"  ✓ OK")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(flagged, f, indent=2)

    print(f"\nFlagged {len(flagged)} suspicious citations")
    print(f"Saved to {output_path}")

if __name__ == '__main__':
    main()

```
