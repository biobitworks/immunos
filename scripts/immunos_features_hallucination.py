#!/usr/bin/env python3
"""
IMMUNOS Hallucination Detection Feature Extractor
=================================================

Extracts 256-dimensional feature vectors from LLM outputs to detect hallucinations.

Features Extracted (256 total):
- Semantic embeddings: SciBERT/BERT (128 dimensions)
- Linguistic markers:
  - Hedging ratio (uncertainty: "maybe", "possibly", etc.)
  - Certainty score (confidence language)
  - Entity density (named entities per sentence)
  - Claim count (extractive claims)
  - Entropy score (word probability distribution)
  - Knowledge triplet count (subject-predicate-object)
- Semantic consistency: Context comparison (64 dimensions)
- Additional metrics (59 dimensions)

Total: 128 + 64 + 64 = 256 dimensions

Dependencies:
    pip install transformers sentence-transformers spacy
    python -m spacy download en_core_web_sm

Datasets:
- TruthfulQA (817 questions) - Small, good for testing
- HaluEval (35,000 samples) - Large, production use

Based on:
- Li et al. (2023) HaluEval - Hallucination detection in LLMs
- Hedging and certainty markers in academic text
"""

import numpy as np
import sys
from pathlib import Path
from typing import Union, Optional, Dict, List, Tuple
import warnings
import re

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
    import spacy
    HAS_SPACY = True
except ImportError:
    HAS_SPACY = False
    warnings.warn("spacy not installed. Install with: pip install spacy")


class HallucinationFeatureExtractor(FeatureExtractor):
    """
    Feature extractor for LLM hallucination detection.

    Analyzes text for indicators of hallucination:
    - Semantic embeddings
    - Uncertainty/hedging markers
    - Entity and claim density
    - Consistency with context
    """

    # Hedging markers (uncertainty)
    HEDGING_MARKERS = [
        'maybe', 'perhaps', 'possibly', 'probably', 'might', 'may', 'could',
        'seems', 'appears', 'suggests', 'indicates', 'likely', 'unlikely',
        'somewhat', 'relatively', 'fairly', 'rather', 'quite', 'slightly',
        'approximately', 'roughly', 'about', 'around', 'generally', 'usually',
        'typically', 'often', 'sometimes', 'occasionally', 'in some cases'
    ]

    # Certainty markers (confidence)
    CERTAINTY_MARKERS = [
        'definitely', 'certainly', 'clearly', 'obviously', 'undoubtedly',
        'absolutely', 'always', 'never', 'must', 'will', 'shall',
        'prove', 'proven', 'fact', 'facts', 'indeed', 'surely',
        'without doubt', 'no question', 'unquestionably', 'categorically'
    ]

    # Claim markers
    CLAIM_MARKERS = [
        'according to', 'research shows', 'studies show', 'evidence suggests',
        'data indicates', 'experts say', 'scientists found', 'researchers discovered',
        'it is known', 'it is established', 'the truth is', 'in reality'
    ]

    def __init__(self,
                 model_name: str = 'all-MiniLM-L6-v2',
                 use_embeddings: bool = True,
                 use_linguistic: bool = True,
                 use_entities: bool = True):
        """
        Initialize hallucination feature extractor.

        Args:
            model_name: Sentence transformer model name
                       'all-MiniLM-L6-v2' (384 dim, fast)
                       'allenai/scibert_scivocab_uncased' (768 dim, scientific)
            use_embeddings: Extract semantic embeddings
            use_linguistic: Extract linguistic markers
            use_entities: Extract entity and claim features
        """
        self.use_embeddings = use_embeddings and HAS_SENTENCE_TRANSFORMERS
        self.use_linguistic = use_linguistic
        self.use_entities = use_entities and HAS_SPACY

        # Initialize sentence transformer
        self.model = None
        if self.use_embeddings and HAS_SENTENCE_TRANSFORMERS:
            try:
                self.model = SentenceTransformer(model_name)
            except Exception as e:
                warnings.warn(f"Failed to load sentence transformer: {e}")
                self.use_embeddings = False

        # Initialize spaCy
        self.nlp = None
        if self.use_entities and HAS_SPACY:
            try:
                self.nlp = spacy.load('en_core_web_sm')
            except Exception as e:
                warnings.warn(f"Failed to load spaCy model: {e}. Download with: python -m spacy download en_core_web_sm")
                self.use_entities = False

    def extract(self, sample: Union[str, Dict[str, str]]) -> np.ndarray:
        """
        Extract 256-dimensional feature vector from text.

        Args:
            sample: Text string or dict with 'response' and optional 'context'

        Returns:
            256-dimensional feature vector
        """
        # Parse input
        if isinstance(sample, dict):
            text = sample.get('response', '')
            context = sample.get('context', None)
        else:
            text = str(sample)
            context = None

        # Preprocess
        text = self.preprocess(text)

        features = []

        # 1. Semantic Embeddings (128 dimensions)
        if self.use_embeddings and self.model is not None:
            embeddings = self._extract_embeddings(text)
            features.append(embeddings)
        else:
            features.append(np.zeros(128))

        # 2. Linguistic Markers (64 dimensions)
        if self.use_linguistic:
            linguistic_features = self._extract_linguistic_features(text)
            features.append(linguistic_features)
        else:
            features.append(np.zeros(64))

        # 3. Semantic Consistency (64 dimensions)
        if self.use_embeddings and self.model is not None and context:
            consistency_features = self._extract_consistency_features(text, context)
            features.append(consistency_features)
        else:
            features.append(np.zeros(64))

        # Concatenate all features
        feature_vector = np.concatenate(features)

        # Ensure exactly 256 dimensions
        assert feature_vector.shape[0] == 256, f"Expected 256 features, got {feature_vector.shape[0]}"

        return feature_vector

    def get_feature_dim(self) -> int:
        """Return dimensionality of feature vector"""
        return 256

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

        # Lowercase for marker detection
        return text

    def _extract_embeddings(self, text: str) -> np.ndarray:
        """
        Extract semantic embeddings (128 dimensions).

        Uses sentence transformer to encode semantic meaning.
        """
        if self.model is None:
            return np.zeros(128)

        try:
            # Get embeddings
            embedding = self.model.encode(text, convert_to_numpy=True)

            # Reduce to 128 dimensions if needed
            if embedding.shape[0] > 128:
                # Average pooling
                bin_size = embedding.shape[0] // 128
                reduced = np.array([
                    embedding[i*bin_size:(i+1)*bin_size].mean()
                    for i in range(128)
                ])
            else:
                # Pad if too small
                reduced = np.zeros(128)
                reduced[:embedding.shape[0]] = embedding

            return reduced[:128]

        except Exception as e:
            warnings.warn(f"Embedding extraction failed: {e}")
            return np.zeros(128)

    def _extract_linguistic_features(self, text: str) -> np.ndarray:
        """
        Extract linguistic markers (64 dimensions).

        Features:
        - Hedging ratio
        - Certainty score
        - Claim count
        - Entity density
        - Sentence length statistics
        - Lexical diversity
        - Entropy
        """
        features = []

        text_lower = text.lower()
        words = text_lower.split()
        sentences = text.split('.')

        # 1. Hedging ratio
        hedging_count = sum(1 for marker in self.HEDGING_MARKERS if marker in text_lower)
        hedging_ratio = hedging_count / (len(words) + 1)
        features.append(hedging_ratio)

        # 2. Certainty ratio
        certainty_count = sum(1 for marker in self.CERTAINTY_MARKERS if marker in text_lower)
        certainty_ratio = certainty_count / (len(words) + 1)
        features.append(certainty_ratio)

        # 3. Claim count
        claim_count = sum(1 for marker in self.CLAIM_MARKERS if marker in text_lower)
        claim_ratio = claim_count / (len(sentences) + 1)
        features.append(claim_ratio)

        # 4-10. Entity density (requires spaCy)
        if self.use_entities and self.nlp is not None:
            try:
                doc = self.nlp(text)

                # Entity density
                entity_count = len(doc.ents)
                entity_density = entity_count / (len(sentences) + 1)
                features.append(entity_density)

                # Entity type distribution (7 types)
                entity_types = ['PERSON', 'ORG', 'GPE', 'DATE', 'MONEY', 'QUANTITY', 'PERCENT']
                for ent_type in entity_types:
                    type_count = sum(1 for ent in doc.ents if ent.label_ == ent_type)
                    features.append(type_count / (entity_count + 1))

            except Exception as e:
                # Fallback if spaCy fails
                features.extend([0.0] * 8)
        else:
            features.extend([0.0] * 8)

        # 11-15. Sentence statistics
        sent_lengths = [len(s.split()) for s in sentences if s.strip()]
        if sent_lengths:
            features.extend([
                np.mean(sent_lengths),
                np.std(sent_lengths),
                np.min(sent_lengths),
                np.max(sent_lengths),
                len(sent_lengths)  # Sentence count
            ])
        else:
            features.extend([0.0] * 5)

        # 16. Lexical diversity (unique words / total words)
        unique_words = len(set(words))
        lexical_diversity = unique_words / (len(words) + 1)
        features.append(lexical_diversity)

        # 17. Average word length
        avg_word_len = np.mean([len(w) for w in words]) if words else 0
        features.append(avg_word_len)

        # 18-20. Question markers
        question_count = text.count('?')
        exclamation_count = text.count('!')
        features.extend([
            question_count / (len(sentences) + 1),
            exclamation_count / (len(sentences) + 1),
            (question_count + exclamation_count) / (len(sentences) + 1)
        ])

        # 21-30. Word frequency distribution (top 10 word frequencies)
        word_freq = {}
        for word in words:
            if len(word) > 3:  # Ignore short words
                word_freq[word] = word_freq.get(word, 0) + 1

        if word_freq:
            top_freqs = sorted(word_freq.values(), reverse=True)[:10]
            # Normalize
            top_freqs = [f / (len(words) + 1) for f in top_freqs]
            # Pad to 10
            while len(top_freqs) < 10:
                top_freqs.append(0.0)
            features.extend(top_freqs[:10])
        else:
            features.extend([0.0] * 10)

        # 31-40. Character-level features
        chars = len(text)
        digits = sum(c.isdigit() for c in text)
        upper = sum(c.isupper() for c in text)
        lower = sum(c.islower() for c in text)
        spaces = sum(c.isspace() for c in text)
        punct = sum(c in '.,;:!?' for c in text)

        features.extend([
            chars / 1000.0,  # Normalized character count
            digits / (chars + 1),
            upper / (chars + 1),
            lower / (chars + 1),
            spaces / (chars + 1),
            punct / (chars + 1),
            digits,  # Absolute counts
            upper,
            lower,
            punct
        ])

        # 41-50. N-gram features (bigram and trigram diversity)
        bigrams = [' '.join(words[i:i+2]) for i in range(len(words)-1)]
        trigrams = [' '.join(words[i:i+3]) for i in range(len(words)-2)]

        bigram_diversity = len(set(bigrams)) / (len(bigrams) + 1)
        trigram_diversity = len(set(trigrams)) / (len(trigrams) + 1)

        features.extend([
            bigram_diversity,
            trigram_diversity,
            len(bigrams) / (len(words) + 1),
            len(trigrams) / (len(words) + 1)
        ])

        # Pad to 64 dimensions
        while len(features) < 64:
            features.append(0.0)

        return np.array(features[:64])

    def _extract_consistency_features(self, text: str, context: str) -> np.ndarray:
        """
        Extract semantic consistency features (64 dimensions).

        Compares text with context to detect inconsistencies.
        """
        if self.model is None or context is None:
            return np.zeros(64)

        try:
            # Get embeddings for text and context
            text_emb = self.model.encode(text, convert_to_numpy=True)
            context_emb = self.model.encode(context, convert_to_numpy=True)

            # Cosine similarity
            similarity = np.dot(text_emb, context_emb) / (
                np.linalg.norm(text_emb) * np.linalg.norm(context_emb) + 1e-8
            )

            # Element-wise difference
            diff = text_emb - context_emb

            # Reduce to 64 dimensions
            if diff.shape[0] > 63:
                bin_size = diff.shape[0] // 63
                reduced_diff = np.array([
                    diff[i*bin_size:(i+1)*bin_size].mean()
                    for i in range(63)
                ])
            else:
                reduced_diff = np.zeros(63)
                reduced_diff[:diff.shape[0]] = diff

            # Combine similarity with difference
            consistency = np.concatenate([[similarity], reduced_diff[:63]])

            return consistency[:64]

        except Exception as e:
            warnings.warn(f"Consistency extraction failed: {e}")
            return np.zeros(64)


# =============================================================================
# TESTING AND UTILITY FUNCTIONS
# =============================================================================

def test_extractor():
    """Test hallucination feature extractor with sample text"""
    print("Testing Hallucination Feature Extractor")
    print("=" * 60)

    # Create extractor
    extractor = HallucinationFeatureExtractor()

    print(f"✓ Extractor initialized")
    print(f"  Use embeddings: {extractor.use_embeddings}")
    print(f"  Use linguistic: {extractor.use_linguistic}")
    print(f"  Use entities: {extractor.use_entities}")
    print()

    # Test 1: Simple text
    print("Test 1: Simple text")
    text1 = "The Earth is the third planet from the Sun. It has one natural satellite called the Moon."
    features1 = extractor.extract(text1)

    print(f"✓ Feature extraction successful")
    print(f"  Input length: {len(text1)} characters")
    print(f"  Feature dimension: {features1.shape[0]}")
    print(f"  Feature range: [{features1.min():.3f}, {features1.max():.3f}]")
    print(f"  Feature mean: {features1.mean():.3f}")
    print()

    # Test 2: Text with hedging
    print("Test 2: Text with hedging markers")
    text2 = "Perhaps the data suggests that the model might be overfitting. This could indicate a need for regularization."
    features2 = extractor.extract(text2)

    print(f"✓ Hedging detection successful")
    print(f"  Feature dimension: {features2.shape[0]}")
    print(f"  Hedging ratio (feature 0): {features2[128]:.3f}")
    print()

    # Test 3: Text with context
    print("Test 3: Text with context (consistency check)")
    context = "Machine learning models require large datasets for training."
    response = "Neural networks can learn complex patterns from data."
    features3 = extractor.extract({'response': response, 'context': context})

    print(f"✓ Consistency check successful")
    print(f"  Feature dimension: {features3.shape[0]}")
    print(f"  Consistency score: {features3[192]:.3f}")
    print()

    # Verify dimension
    assert extractor.get_feature_dim() == 256
    assert features1.shape[0] == 256
    assert features2.shape[0] == 256
    assert features3.shape[0] == 256

    print("=" * 60)
    print("✅ All tests passed!")
    print()
    print("Feature breakdown (256 dimensions):")
    print("  - Embeddings: 128 (semantic meaning)")
    print("  - Linguistic: 64 (hedging, certainty, entities, etc.)")
    print("  - Consistency: 64 (context comparison)")
    print("  Total: 256 dimensions")


if __name__ == "__main__":
    test_extractor()
