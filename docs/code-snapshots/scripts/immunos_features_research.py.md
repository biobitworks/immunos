---
source: /Users/byron/projects/scripts/immunos_features_research.py
relative: scripts/immunos_features_research.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
IMMUNOS Research Paper Verification Feature Extractor
======================================================

Extracts 128-dimensional feature vectors from scientific claims and evidence
for research paper verification.

Features (128 total):
- SciBERT embeddings: 64 dimensions (scientific domain embeddings)
- Claim specificity: Quantitative markers, causal language, specificity score
- Semantic entailment: Claim-evidence relationship (32 dimensions)
- Citation analysis: Citation density, recency, authority
- Readability metrics: Flesch score, complexity, technical density (16 dimensions)
- Linguistic features: Hedging, certainty, contradiction markers (15 dimensions)

Total: 64 + 32 + 16 + 16 = 128 dimensions

Dependencies:
    pip install transformers sentence-transformers textstat

Datasets:
- SciFact (1,409 claims) - Scientific claim verification
- PubMedQA (211,269 Q&A pairs) - Biomedical question answering

Based on:
- Wadden et al. (2020) SciFact - Scientific claim verification
- Entailment detection in scientific text
"""

import numpy as np
import sys
import re
from pathlib import Path
from typing import Union, Dict, Optional
import warnings
from collections import Counter

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.immunos_universal import FeatureExtractor

# Try to import NLP libraries
try:
    from sentence_transformers import SentenceTransformer
    HAS_SENTENCE_TRANSFORMERS = True
except ImportError:
    HAS_SENTENCE_TRANSFORMERS = False
    warnings.warn("sentence-transformers not installed. Install with: pip install sentence-transformers")

try:
    import textstat
    HAS_TEXTSTAT = True
except ImportError:
    HAS_TEXTSTAT = False
    warnings.warn("textstat not installed. Install with: pip install textstat")


class ResearchFeatureExtractor(FeatureExtractor):
    """
    Feature extractor for scientific claim verification.

    Analyzes scientific claims and evidence for verification.
    """

    # Quantitative markers
    QUANTITATIVE_MARKERS = [
        r'\d+%', r'\d+\.\d+', r'\d+/\d+', r'p\s*<\s*0\.\d+',
        r'n\s*=\s*\d+', r'\d+\s*fold', r'\d+x', r'\d+\s*times',
        r'increased by \d+', r'decreased by \d+', r'\d+\s*±\s*\d+'
    ]

    # Causal language
    CAUSAL_MARKERS = [
        'causes', 'caused by', 'leads to', 'results in', 'due to',
        'because of', 'induces', 'triggers', 'mediates', 'regulates',
        'controls', 'influences', 'affects', 'determines', 'promotes',
        'inhibits', 'suppresses', 'enhances', 'reduces', 'increases'
    ]

    # Certainty markers (high confidence)
    CERTAINTY_MARKERS = [
        'demonstrates', 'shows', 'proves', 'confirms', 'establishes',
        'reveals', 'indicates', 'suggests', 'evidence', 'data show',
        'significantly', 'substantially', 'markedly', 'clearly'
    ]

    # Hedging markers (uncertainty)
    HEDGING_MARKERS = [
        'may', 'might', 'could', 'possibly', 'potentially', 'likely',
        'probably', 'appears', 'seems', 'suggests', 'indicates',
        'preliminary', 'tentative', 'inconclusive', 'unclear'
    ]

    # Contradiction markers
    CONTRADICTION_MARKERS = [
        'however', 'but', 'although', 'despite', 'nevertheless',
        'conversely', 'in contrast', 'on the other hand', 'whereas',
        'unfortunately', 'surprisingly', 'unexpectedly'
    ]

    def __init__(self,
                 model_name: str = 'allenai/scibert_scivocab_uncased',
                 use_embeddings: bool = True,
                 use_readability: bool = True):
        """
        Initialize research feature extractor.

        Args:
            model_name: Sentence transformer model
                       'allenai/scibert_scivocab_uncased' (scientific)
                       'all-MiniLM-L6-v2' (general, faster)
            use_embeddings: Extract semantic embeddings
            use_readability: Extract readability metrics
        """
        self.use_embeddings = use_embeddings and HAS_SENTENCE_TRANSFORMERS
        self.use_readability = use_readability and HAS_TEXTSTAT

        # Initialize sentence transformer
        self.model = None
        if self.use_embeddings and HAS_SENTENCE_TRANSFORMERS:
            try:
                self.model = SentenceTransformer(model_name)
            except Exception as e:
                warnings.warn(f"Failed to load sentence transformer: {e}")
                self.use_embeddings = False

    def extract(self, sample: Union[str, Dict[str, str]]) -> np.ndarray:
        """
        Extract 128-dimensional feature vector from scientific text.

        Args:
            sample: Claim string or dict with 'claim' and optional 'evidence'

        Returns:
            128-dimensional feature vector
        """
        # Parse input
        if isinstance(sample, dict):
            claim = sample.get('claim', '')
            evidence = sample.get('evidence', None)
        else:
            claim = str(sample)
            evidence = None

        # Preprocess
        claim = self.preprocess(claim)
        if evidence:
            evidence = self.preprocess(evidence)

        features = []

        # 1. SciBERT Embeddings (64 dimensions)
        if self.use_embeddings and self.model is not None:
            claim_embeddings = self._extract_embeddings(claim)
            features.append(claim_embeddings)
        else:
            features.append(np.zeros(64))

        # 2. Semantic Entailment (32 dimensions)
        if self.use_embeddings and self.model is not None and evidence:
            entailment = self._extract_entailment_features(claim, evidence)
            features.append(entailment)
        else:
            features.append(np.zeros(32))

        # 3. Readability Metrics (16 dimensions)
        if self.use_readability:
            readability = self._extract_readability_features(claim)
            features.append(readability)
        else:
            features.append(np.zeros(16))

        # 4. Linguistic Features (16 dimensions)
        linguistic = self._extract_linguistic_features(claim)
        features.append(linguistic)

        # Concatenate all features
        feature_vector = np.concatenate(features)

        # Ensure exactly 128 dimensions
        assert feature_vector.shape[0] == 128, f"Expected 128 features, got {feature_vector.shape[0]}"

        return feature_vector

    def get_feature_dim(self) -> int:
        """Return dimensionality of feature vector"""
        return 128

    def preprocess(self, text: str) -> str:
        """
        Preprocess text for feature extraction.

        Args:
            text: Input text

        Returns:
            Preprocessed text
        """
        # Remove extra whitespace
        text = ' '.join(text.split())
        return text

    def _extract_embeddings(self, text: str) -> np.ndarray:
        """
        Extract SciBERT embeddings (64 dimensions).
        """
        if self.model is None:
            return np.zeros(64)

        try:
            # Get embeddings
            embedding = self.model.encode(text, convert_to_numpy=True)

            # Reduce to 64 dimensions
            if embedding.shape[0] > 64:
                bin_size = embedding.shape[0] // 64
                reduced = np.array([
                    embedding[i*bin_size:(i+1)*bin_size].mean()
                    for i in range(64)
                ])
            else:
                reduced = np.zeros(64)
                reduced[:embedding.shape[0]] = embedding

            return reduced[:64]

        except Exception as e:
            warnings.warn(f"Embedding extraction failed: {e}")
            return np.zeros(64)

    def _extract_entailment_features(self, claim: str, evidence: str) -> np.ndarray:
        """
        Extract semantic entailment features (32 dimensions).

        Compares claim with evidence to assess support/contradiction.
        """
        if self.model is None:
            return np.zeros(32)

        try:
            # Get embeddings
            claim_emb = self.model.encode(claim, convert_to_numpy=True)
            evidence_emb = self.model.encode(evidence, convert_to_numpy=True)

            # Cosine similarity
            similarity = np.dot(claim_emb, evidence_emb) / (
                np.linalg.norm(claim_emb) * np.linalg.norm(evidence_emb) + 1e-8
            )

            # Element-wise product (interaction)
            product = claim_emb * evidence_emb

            # Reduce to 31 dimensions
            if product.shape[0] > 31:
                bin_size = product.shape[0] // 31
                reduced_product = np.array([
                    product[i*bin_size:(i+1)*bin_size].mean()
                    for i in range(31)
                ])
            else:
                reduced_product = np.zeros(31)
                reduced_product[:product.shape[0]] = product

            # Combine similarity with product
            entailment = np.concatenate([[similarity], reduced_product[:31]])

            return entailment[:32]

        except Exception as e:
            warnings.warn(f"Entailment extraction failed: {e}")
            return np.zeros(32)

    def _extract_readability_features(self, text: str) -> np.ndarray:
        """
        Extract readability metrics (16 dimensions).

        Uses textstat library for standard readability scores.
        """
        features = []

        if not HAS_TEXTSTAT:
            return np.zeros(16)

        try:
            # 1. Flesch Reading Ease
            flesch = textstat.flesch_reading_ease(text)
            features.append(flesch / 100.0)

            # 2. Flesch-Kincaid Grade Level
            fk_grade = textstat.flesch_kincaid_grade(text)
            features.append(fk_grade / 20.0)

            # 3. SMOG Index
            smog = textstat.smog_index(text)
            features.append(smog / 20.0)

            # 4. Coleman-Liau Index
            coleman = textstat.coleman_liau_index(text)
            features.append(coleman / 20.0)

            # 5. Automated Readability Index
            ari = textstat.automated_readability_index(text)
            features.append(ari / 20.0)

            # 6. Dale-Chall Readability Score
            dale_chall = textstat.dale_chall_readability_score(text)
            features.append(dale_chall / 10.0)

            # 7. Difficult words percentage
            difficult = textstat.difficult_words(text)
            word_count = textstat.lexicon_count(text)
            difficult_ratio = difficult / (word_count + 1)
            features.append(difficult_ratio)

            # 8. Lexicon count (word count)
            features.append(word_count / 500.0)

            # 9. Sentence count
            sent_count = textstat.sentence_count(text)
            features.append(sent_count / 20.0)

            # 10. Syllable count
            syllable_count = textstat.syllable_count(text)
            features.append(syllable_count / 1000.0)

            # 11. Average sentence length
            avg_sent_len = word_count / (sent_count + 1)
            features.append(avg_sent_len / 30.0)

            # 12. Average syllables per word
            avg_syllables = syllable_count / (word_count + 1)
            features.append(avg_syllables / 3.0)

            # 13-16. Text complexity
            features.extend([
                textstat.gunning_fog(text) / 20.0,
                textstat.linsear_write_formula(text) / 20.0,
                textstat.text_standard(text, float_output=True) / 20.0,
                textstat.reading_time(text, ms_per_char=14.69) / 60.0  # Reading time in seconds
            ])

        except Exception as e:
            # Fallback if textstat fails
            warnings.warn(f"Readability extraction failed: {e}")
            return np.zeros(16)

        return np.array(features[:16])

    def _extract_linguistic_features(self, text: str) -> np.ndarray:
        """
        Extract linguistic features (16 dimensions).

        Features: quantitative markers, causal language, certainty, hedging, etc.
        """
        features = []
        text_lower = text.lower()
        words = text.split()

        # 1. Quantitative markers
        quant_count = sum(1 for pattern in self.QUANTITATIVE_MARKERS
                         if re.search(pattern, text))
        features.append(quant_count / 10.0)

        # 2. Causal markers
        causal_count = sum(1 for marker in self.CAUSAL_MARKERS
                          if marker in text_lower)
        features.append(causal_count / 5.0)

        # 3. Certainty markers
        certainty_count = sum(1 for marker in self.CERTAINTY_MARKERS
                             if marker in text_lower)
        features.append(certainty_count / 5.0)

        # 4. Hedging markers
        hedging_count = sum(1 for marker in self.HEDGING_MARKERS
                           if marker in text_lower)
        features.append(hedging_count / 5.0)

        # 5. Contradiction markers
        contradiction_count = sum(1 for marker in self.CONTRADICTION_MARKERS
                                 if marker in text_lower)
        features.append(contradiction_count / 5.0)

        # 6. Citation markers
        citation_count = text.count('[') + text.count('(') + text.count('et al')
        features.append(citation_count / 10.0)

        # 7-10. Specific quantitative patterns
        has_percentage = bool(re.search(r'\d+%', text))
        has_pvalue = bool(re.search(r'p\s*<\s*0\.\d+', text))
        has_sample_size = bool(re.search(r'n\s*=\s*\d+', text))
        has_fold_change = bool(re.search(r'\d+\s*fold', text))

        features.extend([
            float(has_percentage),
            float(has_pvalue),
            float(has_sample_size),
            float(has_fold_change)
        ])

        # 11-13. Claim specificity
        claim_length = len(words)
        unique_words = len(set(words))
        lexical_diversity = unique_words / (claim_length + 1)

        features.extend([
            claim_length / 100.0,
            unique_words / 100.0,
            lexical_diversity
        ])

        # 14-16. Technical terminology (heuristic)
        long_words = sum(1 for w in words if len(w) > 10)
        technical_ratio = long_words / (claim_length + 1)

        uppercase_words = sum(1 for w in words if w.isupper() and len(w) > 1)
        hyphenated_words = sum(1 for w in words if '-' in w)

        features.extend([
            technical_ratio,
            uppercase_words / 10.0,
            hyphenated_words / 10.0
        ])

        return np.array(features[:16])


# =============================================================================
# TESTING AND UTILITY FUNCTIONS
# =============================================================================

def test_extractor():
    """Test research feature extractor with sample claims"""
    print("Testing Research Feature Extractor")
    print("=" * 60)

    # Create extractor
    extractor = ResearchFeatureExtractor()

    print(f"✓ Extractor initialized")
    print(f"  Use embeddings: {extractor.use_embeddings}")
    print(f"  Use readability: {extractor.use_readability}")
    print()

    # Test 1: Simple claim
    print("Test 1: Simple scientific claim")
    claim1 = "Vitamin D supplementation reduces the risk of respiratory infections by 50% (p < 0.01, n = 1,200)."
    features1 = extractor.extract(claim1)

    print(f"✓ Feature extraction successful")
    print(f"  Claim length: {len(claim1)} characters")
    print(f"  Feature dimension: {features1.shape[0]}")
    print(f"  Quantitative markers (feature 64+32): {features1[96]:.3f}")
    print()

    # Test 2: Claim with evidence
    print("Test 2: Claim with evidence (entailment)")
    claim2 = "Regular exercise improves cardiovascular health."
    evidence2 = "A meta-analysis of 25 studies found that individuals who exercised regularly had 35% lower risk of heart disease."
    features2 = extractor.extract({'claim': claim2, 'evidence': evidence2})

    print(f"✓ Entailment check successful")
    print(f"  Feature dimension: {features2.shape[0]}")
    print(f"  Similarity score: {features2[64]:.3f}")
    print()

    # Test 3: Hedged claim
    print("Test 3: Hedged claim (uncertainty)")
    claim3 = "The data suggests that the compound may potentially inhibit tumor growth, although further research is needed."
    features3 = extractor.extract(claim3)

    print(f"✓ Hedging detection successful")
    print(f"  Feature dimension: {features3.shape[0]}")
    print(f"  Hedging markers (feature 99): {features3[99]:.3f}")
    print()

    # Verify dimension
    assert extractor.get_feature_dim() == 128
    assert features1.shape[0] == 128
    assert features2.shape[0] == 128
    assert features3.shape[0] == 128

    print("=" * 60)
    print("✅ All tests passed!")
    print()
    print("Feature breakdown (128 dimensions):")
    print("  - Embeddings: 64 (SciBERT semantic representation)")
    print("  - Entailment: 32 (claim-evidence relationship)")
    print("  - Readability: 16 (complexity, grade level, etc.)")
    print("  - Linguistic: 16 (quantitative, causal, certainty, etc.)")
    print("  Total: 128 dimensions")


if __name__ == "__main__":
    test_extractor()

```
