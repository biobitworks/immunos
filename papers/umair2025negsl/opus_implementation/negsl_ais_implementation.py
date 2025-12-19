"""
NegSl-AIS Core Algorithms for IMMUNOS Implementation
=====================================================

This module implements the core Negative Selection Artificial Immune System
algorithms based on the research paper by Umair et al. (2025).

These algorithms can be adapted for AI-assisted development security scanning.
"""

import numpy as np
from typing import List, Tuple, Dict, Set, Optional
from dataclasses import dataclass
from enum import Enum
import math

# =============================================================================
# DATA STRUCTURES
# =============================================================================

class ClassLabel(Enum):
    """Classification labels based on Circumplex Model"""
    LOW_AROUSAL = "LA"
    HIGH_AROUSAL = "HA"
    LOW_VALENCE = "LV"
    HIGH_VALENCE = "HV"
    
    # IMMUNOS equivalents
    SAFE = "SAFE"
    MALICIOUS = "MALICIOUS"
    SUSPICIOUS = "SUSPICIOUS"
    UNKNOWN = "UNKNOWN"


@dataclass
class Detector:
    """
    Represents a valid detector in the Negative Selection algorithm.
    
    In biological terms: An immunocompetent T-cell that has passed
    negative selection in the thymus.
    
    In security terms: A verified pattern signature that identifies
    non-self (potentially malicious) code.
    """
    center: np.ndarray  # d^j - Center point in feature space
    radius: float       # r^j - Detection radius
    class_label: str    # Which class this detector belongs to
    
    def __post_init__(self):
        """Validate detector parameters"""
        if self.radius < 0:
            raise ValueError("Detector radius must be non-negative")


@dataclass
class DetectorConfig:
    """
    Configuration parameters for detector generation.
    Based on optimal values from the paper.
    """
    num_detectors: int      # Number of detectors to generate
    r_self: float           # Self-radius threshold
    class_label: str        # Target class
    
    # Optimal configurations from paper (Table 6)
    @classmethod
    def get_optimal_config(cls, class_label: str) -> 'DetectorConfig':
        """Return optimal configuration for each class"""
        configs = {
            'LA': cls(num_detectors=15, r_self=0.89, class_label='LA'),
            'HA': cls(num_detectors=15, r_self=0.91, class_label='HA'),
            'LV': cls(num_detectors=25, r_self=1.31, class_label='LV'),
            'HV': cls(num_detectors=20, r_self=1.33, class_label='HV'),
        }
        return configs.get(class_label, cls(num_detectors=20, r_self=1.0, class_label=class_label))


# =============================================================================
# CORE ALGORITHMS
# =============================================================================

class NegativeSelectionClassifier:
    """
    Negative Selection-based Artificial Immune System Classifier
    
    Implements the NegSl-AIS algorithm for binary classification using
    the self/non-self discrimination principle from immunology.
    
    Key Concepts:
    - Self samples: Training samples from the target class
    - Non-self samples: Samples from other classes (to be detected)
    - Detectors: Random patterns that don't match self samples
    - R_self: Threshold radius defining the "self" region
    
    Algorithm Steps:
    1. Generate random candidate detectors
    2. Remove detectors that match (bind to) self samples
    3. Remaining detectors can identify non-self samples
    """
    
    def __init__(self, config: DetectorConfig):
        """
        Initialize the classifier with configuration.
        
        Args:
            config: DetectorConfig with num_detectors, r_self, class_label
        """
        self.config = config
        self.valid_detectors: List[Detector] = []
        self.self_samples: Optional[np.ndarray] = None
        self.feature_dim: int = 0
        
    def _euclidean_distance(self, point_a: np.ndarray, point_b: np.ndarray) -> float:
        """
        Calculate Euclidean distance between two points.
        
        Formula: R_q = sqrt(sum((Ω_k^i - d_k^j)^2))
        
        Args:
            point_a: First point (e.g., self sample center Ω^i)
            point_b: Second point (e.g., detector center d^j)
            
        Returns:
            Euclidean distance
        """
        return np.sqrt(np.sum((point_a - point_b) ** 2))
    
    def _get_nearest_self_distance(self, detector_center: np.ndarray) -> float:
        """
        Calculate minimum distance from detector to nearest self sample.
        
        Formula: R_q = min(||Self_i, d_j||) for i = 1 to N
        
        Args:
            detector_center: Center point of candidate detector (d^j)
            
        Returns:
            R_q - distance to nearest self sample
        """
        min_distance = float('inf')
        
        for self_sample in self.self_samples:
            distance = self._euclidean_distance(self_sample, detector_center)
            min_distance = min(min_distance, distance)
            
        return min_distance
    
    def _is_valid_detector(self, detector_center: np.ndarray) -> Tuple[bool, float]:
        """
        Check if a candidate detector is valid (doesn't bind to self).
        
        Validation Rule (from paper):
            Detector = {
                Valid,   if R_q > R_self
                Invalid, if R_q < R_self
            }
        
        Args:
            detector_center: Center point of candidate detector
            
        Returns:
            Tuple of (is_valid, R_q)
        """
        r_q = self._get_nearest_self_distance(detector_center)
        is_valid = r_q > self.config.r_self
        return is_valid, r_q
    
    def _generate_random_detector(self) -> np.ndarray:
        """
        Generate a random candidate detector in the feature space.
        
        Returns:
            Random vector of shape (feature_dim,) with values in [0, 1]
        """
        return np.random.uniform(0, 1, self.feature_dim)
    
    def fit(self, self_samples: np.ndarray, max_attempts: int = 1000) -> 'NegativeSelectionClassifier':
        """
        Train the classifier by generating valid detectors.
        
        This implements the TRAINING phase from the paper:
        1. Initialize empty detector set D
        2. Generate random candidate detectors
        3. For each candidate, check if it binds to self samples
        4. Keep only detectors that DON'T bind to self (R_q > R_self)
        
        Args:
            self_samples: Training samples for target class (m x n matrix)
            max_attempts: Maximum attempts to generate required detectors
            
        Returns:
            self (fitted classifier)
        """
        self.self_samples = np.array(self_samples)
        self.feature_dim = self_samples.shape[1]
        self.valid_detectors = []
        
        attempts = 0
        
        while len(self.valid_detectors) < self.config.num_detectors and attempts < max_attempts:
            # Generate random candidate detector
            candidate_center = self._generate_random_detector()
            
            # Check if valid (doesn't bind to self)
            is_valid, r_q = self._is_valid_detector(candidate_center)
            
            if is_valid:
                # Calculate detector radius: r^j = R^q - R^self
                detector_radius = r_q - self.config.r_self
                
                # Create valid detector
                detector = Detector(
                    center=candidate_center,
                    radius=detector_radius,
                    class_label=self.config.class_label
                )
                
                # Check for duplicates
                is_duplicate = any(
                    np.allclose(d.center, candidate_center) 
                    for d in self.valid_detectors
                )
                
                if not is_duplicate:
                    self.valid_detectors.append(detector)
            
            attempts += 1
        
        print(f"Generated {len(self.valid_detectors)} valid detectors for class {self.config.class_label}")
        return self
    
    def predict_single(self, sample: np.ndarray) -> Tuple[str, float]:
        """
        Classify a single sample using valid detectors.
        
        This implements the TESTING phase from the paper:
        1. For each valid detector, calculate distance to sample
        2. If any detector is within R_self of sample, classify as "self"
        3. Otherwise, classify as "non-self"
        
        Args:
            sample: Feature vector to classify
            
        Returns:
            Tuple of (classification, confidence)
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
        
        # If sample is within R_self of any detector, it's classified as "self"
        if min_distance < self.config.r_self:
            confidence = 1.0 - (min_distance / self.config.r_self)
            return "self", confidence
        else:
            confidence = min(1.0, (min_distance - self.config.r_self) / self.config.r_self)
            return "non-self", confidence
    
    def predict(self, samples: np.ndarray) -> List[Tuple[str, float]]:
        """
        Classify multiple samples.
        
        Args:
            samples: Array of feature vectors (m x n matrix)
            
        Returns:
            List of (classification, confidence) tuples
        """
        return [self.predict_single(sample) for sample in samples]


# =============================================================================
# MODALITY BIASING
# =============================================================================

class ModalityBiasing:
    """
    Implements modality biasing for feature fusion.
    
    Based on unimodal classification accuracies, each modality
    is assigned a weight that sums to 1.0.
    
    Weights from paper (Table 4):
    - EEG: 0.28 (58.25% accuracy)
    - ECG: 0.26 (53.85% accuracy)
    - Respiration: 0.25 (46.38% accuracy)
    - GSR: 0.14 (28.89% accuracy)
    - Body Temp: 0.07 (27.99% accuracy)
    """
    
    # Default weights from the paper
    DEFAULT_WEIGHTS = {
        'EEG': 0.28,
        'ECG': 0.26,
        'RESPIRATION': 0.25,
        'GSR': 0.14,
        'BODY_TEMP': 0.07
    }
    
    def __init__(self, weights: Optional[Dict[str, float]] = None):
        """
        Initialize with modality weights.
        
        Args:
            weights: Dictionary mapping modality names to weights.
                    If None, uses default weights from paper.
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
    
    def apply_bias(self, features: Dict[str, np.ndarray]) -> np.ndarray:
        """
        Apply modality biasing to features from different sources.
        
        Args:
            features: Dictionary mapping modality names to feature arrays
            
        Returns:
            Weighted and concatenated feature vector
        """
        weighted_features = []
        
        for modality, feature_array in features.items():
            if modality in self.weights:
                weight = self.weights[modality]
                weighted = feature_array * weight
                weighted_features.append(weighted)
            else:
                print(f"Warning: No weight defined for modality '{modality}'")
                weighted_features.append(feature_array)
        
        return np.concatenate(weighted_features)
    
    @classmethod
    def from_accuracies(cls, accuracies: Dict[str, float]) -> 'ModalityBiasing':
        """
        Create biasing weights from unimodal classification accuracies.
        
        Args:
            accuracies: Dictionary mapping modality names to accuracy percentages
            
        Returns:
            ModalityBiasing instance with normalized weights
        """
        total = sum(accuracies.values())
        weights = {k: v / total for k, v in accuracies.items()}
        return cls(weights)


# =============================================================================
# FEATURE EXTRACTION
# =============================================================================

class FeatureExtractor:
    """
    Feature extraction utilities based on the paper.
    
    Implements statistical and signal processing features for
    multimodal data analysis.
    """
    
    @staticmethod
    def standard_deviation(x: np.ndarray) -> float:
        """
        Calculate standard deviation.
        
        Formula: SD = sqrt(1/(n-1) * sum((x_i - x')^2))
        """
        return np.std(x, ddof=1)
    
    @staticmethod
    def skewness(x: np.ndarray) -> float:
        """
        Calculate skewness (asymmetry of distribution).
        
        Formula: Skewness = (n/((n-1)(n-2))) * sum(((x_i - x')/SD)^3)
        """
        n = len(x)
        if n < 3:
            return 0.0
        
        mean = np.mean(x)
        std = np.std(x, ddof=1)
        
        if std == 0:
            return 0.0
        
        return (n / ((n - 1) * (n - 2))) * np.sum(((x - mean) / std) ** 3)
    
    @staticmethod
    def zero_crossing_rate(x: np.ndarray) -> float:
        """
        Calculate zero crossing rate.
        
        Formula: ZCR = (1/n) * sum(|Sign(x_i) - Sign(x_{i+1})|)
        """
        signs = np.sign(x)
        # Handle zeros - treat as positive
        signs[signs == 0] = 1
        
        crossings = np.abs(np.diff(signs))
        return np.sum(crossings > 0) / len(x)
    
    @staticmethod
    def kurtosis(x: np.ndarray) -> float:
        """
        Calculate kurtosis (tailedness of distribution).
        
        Formula: Kurtosis = (1/n) * sum(((x_i - μ)/SD)^4) - 3
        """
        n = len(x)
        mean = np.mean(x)
        std = np.std(x, ddof=0)
        
        if std == 0:
            return 0.0
        
        return (1 / n) * np.sum(((x - mean) / std) ** 4) - 3
    
    @staticmethod
    def waveform_length(x: np.ndarray) -> float:
        """
        Calculate waveform length.
        
        Formula: WL = sum(|x_{i+1} - x_i|)
        """
        return np.sum(np.abs(np.diff(x)))
    
    @staticmethod
    def energy(x: np.ndarray) -> float:
        """
        Calculate signal energy.
        
        Formula: Energy = sum(|x_i|^2)
        """
        return np.sum(np.abs(x) ** 2)
    
    @staticmethod
    def entropy(x: np.ndarray, bins: int = 10) -> float:
        """
        Calculate Shannon entropy.
        
        Formula: Entropy = -sum(P(x_i) * log2(P(x_i)))
        """
        hist, _ = np.histogram(x, bins=bins, density=True)
        hist = hist[hist > 0]  # Remove zero bins
        
        # Normalize to get probabilities
        probs = hist / np.sum(hist)
        
        return -np.sum(probs * np.log2(probs))
    
    @staticmethod
    def power_spectral_density(x: np.ndarray) -> np.ndarray:
        """
        Calculate Power Spectral Density using FFT.
        
        Formula: PSD = |FFT(x)|^2
        """
        fft_result = np.fft.fft(x)
        return np.abs(fft_result) ** 2
    
    def extract_all_features(self, signal: np.ndarray) -> Dict[str, float]:
        """
        Extract all statistical features from a signal.
        
        Args:
            signal: 1D numpy array of signal values
            
        Returns:
            Dictionary of feature names to values
        """
        return {
            'mean': np.mean(signal),
            'std': self.standard_deviation(signal),
            'max': np.max(signal),
            'min': np.min(signal),
            'skewness': self.skewness(signal),
            'kurtosis': self.kurtosis(signal),
            'zcr': self.zero_crossing_rate(signal),
            'waveform_length': self.waveform_length(signal),
            'energy': self.energy(signal),
            'entropy': self.entropy(signal)
        }


# =============================================================================
# GENERALIZATION ERROR
# =============================================================================

def calculate_generalization_error(
    training_predictions: List[Tuple[str, float]],
    training_labels: List[str],
    test_predictions: List[Tuple[str, float]],
    test_labels: List[str]
) -> float:
    """
    Calculate generalization error.
    
    Formula: E_gen = |training_error - test_error|
    
    Args:
        training_predictions: List of (prediction, confidence) from training
        training_labels: True labels for training samples
        test_predictions: List of (prediction, confidence) from testing
        test_labels: True labels for test samples
        
    Returns:
        Generalization error (absolute difference between train and test error)
    """
    # Calculate training error
    train_correct = sum(
        1 for (pred, _), true in zip(training_predictions, training_labels)
        if pred == true
    )
    train_error = 1.0 - (train_correct / len(training_labels))
    
    # Calculate test error
    test_correct = sum(
        1 for (pred, _), true in zip(test_predictions, test_labels)
        if pred == true
    )
    test_error = 1.0 - (test_correct / len(test_labels))
    
    return abs(train_error - test_error)


# =============================================================================
# EXAMPLE USAGE FOR IMMUNOS
# =============================================================================

def immunos_security_example():
    """
    Example of how to adapt NegSl-AIS for IMMUNOS security scanning.
    
    In this context:
    - "Self" = Known-safe code patterns
    - "Non-self" = Potentially malicious patterns
    - Detectors = Security signatures
    """
    
    # Simulated safe code feature vectors (training data)
    # In practice, these would be extracted from known-good code
    np.random.seed(42)
    safe_code_features = np.random.rand(100, 50)  # 100 samples, 50 features
    
    # Configure detector for "safe" class
    config = DetectorConfig(
        num_detectors=20,
        r_self=0.8,
        class_label="SAFE"
    )
    
    # Train the classifier
    classifier = NegativeSelectionClassifier(config)
    classifier.fit(safe_code_features)
    
    # Test on new code samples
    test_safe = np.random.rand(10, 50) * 0.5  # Similar to safe patterns
    test_malicious = np.random.rand(10, 50) + 0.5  # Different patterns
    
    # Classify
    safe_results = classifier.predict(test_safe)
    malicious_results = classifier.predict(test_malicious)
    
    print("\n=== IMMUNOS Security Classification Results ===")
    print(f"\nSafe code samples:")
    for i, (label, confidence) in enumerate(safe_results):
        print(f"  Sample {i+1}: {label} (confidence: {confidence:.3f})")
    
    print(f"\nSuspicious code samples:")
    for i, (label, confidence) in enumerate(malicious_results):
        print(f"  Sample {i+1}: {label} (confidence: {confidence:.3f})")


if __name__ == "__main__":
    immunos_security_example()
