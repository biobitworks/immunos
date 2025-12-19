#!/usr/bin/env python3
"""
IMMUNOS Negative Selection Module
==================================

Implements the NegSl-AIS algorithm for code security anomaly detection.
Based on: Umair et al. (2025) - Negative selection-based artificial immune system

Core Concept:
- Train detectors on "safe" code patterns (self)
- Detectors flag anything that doesn't match (non-self/threats)
- Biological metaphor: T-cell maturation via thymic selection

Usage:
    from immunos_negsel import NegativeSelectionClassifier, CodeFeatureExtractor

    # Extract features from safe code
    extractor = CodeFeatureExtractor()
    safe_features = extractor.extract_from_files(safe_code_files)

    # Train detector
    classifier = NegativeSelectionClassifier(num_detectors=20, r_self=0.85)
    classifier.fit(safe_features)

    # Classify new code
    test_features = extractor.extract_from_file('test.py')
    classification, confidence = classifier.predict_single(test_features)
"""

import numpy as np
import sqlite3
import pickle
import hashlib
import ast
from pathlib import Path
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass
from datetime import datetime
from enum import Enum


# =============================================================================
# DATA STRUCTURES
# =============================================================================

class ClassLabel(Enum):
    """Classification labels for code security"""
    SAFE = "SAFE"
    MALICIOUS = "MALICIOUS"
    SUSPICIOUS = "SUSPICIOUS"
    UNKNOWN = "UNKNOWN"


@dataclass
class Detector:
    """
    Represents a valid detector that passed negative selection.

    Biological: Immunocompetent T-cell
    Security: Verified threat signature
    """
    center: np.ndarray  # Center point in feature space
    radius: float       # Detection radius (r^j = R_q - R_self)
    class_label: str    # Which class this detector identifies

    def __post_init__(self):
        if self.radius < 0:
            raise ValueError("Detector radius must be non-negative")


# =============================================================================
# NEGATIVE SELECTION CLASSIFIER
# =============================================================================

class NegativeSelectionClassifier:
    """
    Negative Selection Classifier for code security.

    Algorithm:
    1. Training: Generate random detectors, discard those matching safe code
    2. Testing: Flag code that matches remaining detectors as potentially malicious

    Parameters from NegSl-AIS paper:
    - num_detectors: 15-25 (optimal range)
    - r_self: 0.87-1.34 (threshold for self-radius)
    """

    def __init__(self, num_detectors: int = 20, r_self: float = 0.85,
                 class_label: str = "SAFE", db_path: str = None):
        """
        Initialize classifier.

        Args:
            num_detectors: Number of detectors to generate (15-25 optimal)
            r_self: Self-radius threshold (0.85-1.0 recommended for IMMUNOS)
            class_label: Target class to learn
            db_path: Path to IMMUNOS database for persistence
        """
        self.num_detectors = num_detectors
        self.r_self = r_self
        self.class_label = class_label
        self.db_path = db_path or "/Users/byron/projects/.immunos/db/immunos.db"

        self.valid_detectors: List[Detector] = []
        self.self_samples: Optional[np.ndarray] = None
        self.feature_dim: int = 0

    def _euclidean_distance(self, point_a: np.ndarray, point_b: np.ndarray) -> float:
        """
        Calculate Euclidean distance.

        Formula: R_q = sqrt(sum((a_k - b_k)^2))
        """
        return np.sqrt(np.sum((point_a - point_b) ** 2))

    def _get_nearest_self_distance(self, detector_center: np.ndarray) -> float:
        """
        Calculate minimum distance from detector to nearest self sample.

        Formula: R_q = min(||Self_i, d_j||) for i = 1 to N
        """
        min_distance = float('inf')

        for self_sample in self.self_samples:
            distance = self._euclidean_distance(self_sample, detector_center)
            min_distance = min(min_distance, distance)

        return min_distance

    def _is_valid_detector(self, detector_center: np.ndarray) -> Tuple[bool, float]:
        """
        Check if candidate detector is valid.

        Validation Rule (Equation 20 from paper):
            Detector = {
                Valid,   if R_q > R_self
                Invalid, if R_q < R_self
            }

        Returns:
            (is_valid, R_q)
        """
        r_q = self._get_nearest_self_distance(detector_center)
        is_valid = r_q > self.r_self
        return is_valid, r_q

    def _generate_random_detector(self) -> np.ndarray:
        """Generate random candidate detector in [0, 1] range."""
        return np.random.uniform(0, 1, self.feature_dim)

    def fit(self, self_samples: np.ndarray, max_attempts: int = 1000) -> 'NegativeSelectionClassifier':
        """
        Train classifier by generating valid detectors.

        Training Phase:
        1. Generate random candidate detectors
        2. Test if they bind to self samples
        3. Keep only detectors that DON'T bind (R_q > R_self)

        Args:
            self_samples: Training samples (m x n matrix of safe code features)
            max_attempts: Maximum attempts to generate required detectors

        Returns:
            self (fitted classifier)
        """
        self.self_samples = np.array(self_samples)
        self.feature_dim = self_samples.shape[1]
        self.valid_detectors = []

        attempts = 0

        while len(self.valid_detectors) < self.num_detectors and attempts < max_attempts:
            candidate_center = self._generate_random_detector()
            is_valid, r_q = self._is_valid_detector(candidate_center)

            if is_valid:
                # Calculate detector radius: r^j = R^q - R^self
                detector_radius = r_q - self.r_self

                detector = Detector(
                    center=candidate_center,
                    radius=detector_radius,
                    class_label=self.class_label
                )

                # Check for duplicates
                is_duplicate = any(
                    np.allclose(d.center, candidate_center)
                    for d in self.valid_detectors
                )

                if not is_duplicate:
                    self.valid_detectors.append(detector)

            attempts += 1

        print(f"Generated {len(self.valid_detectors)} valid detectors for class {self.class_label}")
        print(f"  Attempts: {attempts}, Success rate: {len(self.valid_detectors)/attempts*100:.1f}%")

        return self

    def predict_single(self, sample: np.ndarray) -> Tuple[str, float]:
        """
        Classify a single code sample.

        Testing Phase:
        1. Calculate distance to each detector
        2. If within R_self of any detector → classify as "self" (safe)
        3. Otherwise → classify as "non-self" (potentially malicious)

        Args:
            sample: Feature vector to classify

        Returns:
            (classification, confidence)
        """
        if not self.valid_detectors:
            raise ValueError("Classifier not fitted. Call fit() first.")

        min_distance = float('inf')
        binding_detector = None

        for detector in self.valid_detectors:
            distance = self._euclidean_distance(sample, detector.center)

            if distance < min_distance:
                min_distance = distance
                binding_detector = detector

        # Classification logic
        if min_distance < self.r_self:
            # Sample is within self-radius → Safe
            confidence = 1.0 - (min_distance / self.r_self)
            return "self", confidence
        else:
            # Sample is outside self-radius → Non-self (threat)
            confidence = min(1.0, (min_distance - self.r_self) / self.r_self)
            return "non-self", confidence

    def predict(self, samples: np.ndarray) -> List[Tuple[str, float]]:
        """Classify multiple samples."""
        return [self.predict_single(sample) for sample in samples]

    def save_to_db(self, accuracy: float = None, false_positive_rate: float = None):
        """
        Persist trained detectors to database.

        Args:
            accuracy: Validation accuracy
            false_positive_rate: False positive rate on validation set
        """
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        for detector in self.valid_detectors:
            # Serialize numpy array
            center_blob = pickle.dumps(detector.center)

            cursor.execute("""
                INSERT INTO detectors
                (class_label, center_vector, radius, r_self, feature_dim,
                 training_date, accuracy, false_positive_rate)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                self.class_label,
                center_blob,
                detector.radius,
                self.r_self,
                self.feature_dim,
                datetime.now().isoformat(),
                accuracy,
                false_positive_rate
            ))

        conn.commit()
        conn.close()
        print(f"✓ Saved {len(self.valid_detectors)} detectors to database")

    @classmethod
    def load_from_db(cls, class_label: str, db_path: str = None) -> 'NegativeSelectionClassifier':
        """
        Load trained detectors from database.

        Args:
            class_label: Which class to load detectors for
            db_path: Path to IMMUNOS database

        Returns:
            Fitted classifier with loaded detectors
        """
        db_path = db_path or "/Users/byron/projects/.immunos/db/immunos.db"
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        cursor.execute("""
            SELECT center_vector, radius, r_self, feature_dim
            FROM detectors
            WHERE class_label = ?
            ORDER BY training_date DESC
        """, (class_label,))

        rows = cursor.fetchall()
        conn.close()

        if not rows:
            raise ValueError(f"No detectors found for class '{class_label}'")

        # Reconstruct classifier
        r_self = rows[0][2]
        feature_dim = rows[0][3]

        classifier = cls(
            num_detectors=len(rows),
            r_self=r_self,
            class_label=class_label,
            db_path=db_path
        )

        classifier.feature_dim = feature_dim
        classifier.valid_detectors = []

        for row in rows:
            center = pickle.loads(row[0])
            radius = row[1]

            detector = Detector(
                center=center,
                radius=radius,
                class_label=class_label
            )
            classifier.valid_detectors.append(detector)

        print(f"✓ Loaded {len(classifier.valid_detectors)} detectors from database")
        return classifier


# =============================================================================
# MODALITY BIASING (Scanner Fusion)
# =============================================================================

class ScannerFusion:
    """
    Implements modality biasing for combining multiple security scanners.

    Weights scanners by their individual accuracy (must sum to 1.0).

    IMMUNOS scanners:
    - Bandit: Python security issues
    - Semgrep: Pattern matching
    - CodeQL: Semantic analysis
    - Custom: AI-based detection
    """

    # Default weights based on expected scanner accuracy
    DEFAULT_WEIGHTS = {
        'bandit': 0.25,    # Python-specific security
        'semgrep': 0.30,   # Pattern matching (highest weight)
        'codeql': 0.25,    # Semantic analysis
        'custom': 0.20     # AI/LLM detection
    }

    def __init__(self, weights: Optional[Dict[str, float]] = None):
        """
        Initialize scanner fusion.

        Args:
            weights: Dictionary mapping scanner names to weights.
                    If None, uses default weights.
        """
        self.weights = weights or self.DEFAULT_WEIGHTS.copy()
        self._validate_weights()

    def _validate_weights(self):
        """Ensure weights sum to 1.0"""
        total = sum(self.weights.values())
        if not np.isclose(total, 1.0):
            # Normalize weights
            for key in self.weights:
                self.weights[key] /= total

    def fuse_detections(self, scanner_results: Dict[str, float]) -> float:
        """
        Fuse detection scores from multiple scanners.

        Args:
            scanner_results: Dict mapping scanner name to detection score (0-1)

        Returns:
            Weighted combined score (0-1)
        """
        weighted_sum = 0.0
        total_weight = 0.0

        for scanner, score in scanner_results.items():
            if scanner in self.weights:
                weight = self.weights[scanner]
                weighted_sum += score * weight
                total_weight += weight

        if total_weight == 0:
            return 0.0

        return weighted_sum / total_weight

    @classmethod
    def from_accuracies(cls, accuracies: Dict[str, float]) -> 'ScannerFusion':
        """
        Create fusion weights from scanner accuracies.

        Args:
            accuracies: Dict mapping scanner names to accuracy percentages

        Returns:
            ScannerFusion instance with normalized weights
        """
        total = sum(accuracies.values())
        weights = {k: v / total for k, v in accuracies.items()}
        return cls(weights)


# =============================================================================
# CODE FEATURE EXTRACTION
# =============================================================================

class CodeFeatureExtractor:
    """
    Extract statistical and structural features from code files.

    Features extracted:
    - Structural: Lines, complexity, nesting depth, function count
    - Statistical: Character entropy, keyword density, comment ratio
    - AST-based: Node counts, import patterns, call patterns
    """

    def extract_from_file(self, file_path: str) -> np.ndarray:
        """
        Extract feature vector from a single code file.

        Args:
            file_path: Path to code file

        Returns:
            Feature vector (numpy array)
        """
        try:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
        except Exception:
            # Return zero vector for unreadable files
            return np.zeros(30)

        features = []

        # Structural features
        lines = content.splitlines()
        features.append(len(lines))  # Total lines
        features.append(len([l for l in lines if l.strip()]))  # Non-empty lines
        features.append(len([l for l in lines if l.strip().startswith('#')]))  # Comment lines

        # Character-level statistics
        features.append(len(content))  # Total characters
        features.append(content.count(' ') / max(len(content), 1))  # Space density
        features.append(content.count('\t') / max(len(content), 1))  # Tab density
        features.append(self._calculate_entropy(content))  # Shannon entropy

        # Python keyword density (if Python file)
        if file_path.endswith('.py'):
            keywords = ['def', 'class', 'import', 'if', 'for', 'while', 'try', 'except']
            for keyword in keywords:
                features.append(content.count(keyword) / max(len(lines), 1))
        else:
            features.extend([0.0] * 8)  # Padding for non-Python files

        # AST features (for Python)
        if file_path.endswith('.py'):
            ast_features = self._extract_ast_features(content)
            features.extend(ast_features)
        else:
            features.extend([0.0] * 15)  # Padding

        return np.array(features, dtype=np.float32)

    def extract_from_files(self, file_paths: List[str]) -> np.ndarray:
        """
        Extract features from multiple files.

        Args:
            file_paths: List of file paths

        Returns:
            Feature matrix (n_files x n_features)
        """
        features = [self.extract_from_file(fp) for fp in file_paths]
        return np.array(features)

    def _calculate_entropy(self, text: str) -> float:
        """Calculate Shannon entropy of text."""
        if not text:
            return 0.0

        # Count character frequencies
        freq = {}
        for char in text:
            freq[char] = freq.get(char, 0) + 1

        # Calculate probabilities and entropy
        total = len(text)
        entropy = 0.0
        for count in freq.values():
            p = count / total
            if p > 0:
                entropy -= p * np.log2(p)

        return entropy

    def _extract_ast_features(self, code: str) -> List[float]:
        """Extract features from Python AST."""
        features = [0.0] * 15  # Default zero features

        try:
            tree = ast.parse(code)

            # Count node types
            node_counts = {
                'FunctionDef': 0,
                'ClassDef': 0,
                'Import': 0,
                'ImportFrom': 0,
                'If': 0,
                'For': 0,
                'While': 0,
                'Try': 0,
                'With': 0,
                'Call': 0,
                'BinOp': 0,
                'Compare': 0,
                'Assign': 0,
                'Return': 0,
                'Raise': 0
            }

            for node in ast.walk(tree):
                node_type = type(node).__name__
                if node_type in node_counts:
                    node_counts[node_type] += 1

            features = list(node_counts.values())

        except SyntaxError:
            # Invalid Python syntax - return zeros
            pass

        return features


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def train_safe_code_detector(safe_code_dir: str, num_detectors: int = 20,
                            r_self: float = 0.85) -> NegativeSelectionClassifier:
    """
    Train a detector on a directory of safe code.

    Args:
        safe_code_dir: Path to directory containing safe code samples
        num_detectors: Number of detectors to generate
        r_self: Self-radius threshold

    Returns:
        Fitted classifier
    """
    # Find all Python files
    safe_files = list(Path(safe_code_dir).rglob('*.py'))
    print(f"Found {len(safe_files)} safe code files")

    # Extract features
    extractor = CodeFeatureExtractor()
    safe_features = extractor.extract_from_files([str(f) for f in safe_files])
    print(f"Extracted features: {safe_features.shape}")

    # Train classifier
    classifier = NegativeSelectionClassifier(
        num_detectors=num_detectors,
        r_self=r_self,
        class_label="SAFE"
    )
    classifier.fit(safe_features)

    # Save to database
    classifier.save_to_db()

    return classifier


if __name__ == "__main__":
    import sys
    import argparse

    parser = argparse.ArgumentParser(description='IMMUNOS Negative Selection Trainer')
    parser.add_argument('command', choices=['train', 'info'], help='Command to execute')
    parser.add_argument('safe_code_dir', nargs='?', help='Directory containing safe code samples')
    parser.add_argument('--num-detectors', type=int, default=20,
                       help='Number of detectors to generate (default: 20)')
    parser.add_argument('--r-self', type=float, default=0.85,
                       help='Self-radius threshold (default: 0.85)')

    args = parser.parse_args()

    if args.command == 'info':
        print("IMMUNOS Negative Selection Module")
        print("=" * 50)
        print("\nThis module provides the core NegSl-AIS algorithm for IMMUNOS.")
        print("\nUsage:")
        print("  Train:  python immunos_negsel.py train /path/to/safe/code")
        print("  Info:   python immunos_negsel.py info")
        print("\nExample:")
        print("  python scripts/immunos_negsel.py train scripts/ --num-detectors 20")
        print("\nIntegration:")
        print("  python scripts/immunos_nk_scan.py --mode negsel")
        print("  python scripts/immunos_nk_scan.py --mode hybrid")

    elif args.command == 'train':
        if not args.safe_code_dir:
            print("Error: Please provide a directory containing safe code samples")
            print("Usage: python immunos_negsel.py train <safe_code_dir>")
            sys.exit(1)

        print("IMMUNOS Detector Training")
        print("=" * 50)
        print(f"Safe code directory: {args.safe_code_dir}")
        print(f"Number of detectors: {args.num_detectors}")
        print(f"R_self threshold: {args.r_self}")
        print()

        classifier = train_safe_code_detector(
            safe_code_dir=args.safe_code_dir,
            num_detectors=args.num_detectors,
            r_self=args.r_self
        )

        print("\n✓ Training complete!")
        print(f"  Detectors generated: {len(classifier.valid_detectors)}")
        print(f"  Feature dimensionality: {classifier.feature_dim}")
        print("\nNext steps:")
        print("  1. Test the detector: python scripts/immunos_nk_scan.py --mode negsel")
        print("  2. Use hybrid mode: python scripts/immunos_nk_scan.py --mode hybrid")
