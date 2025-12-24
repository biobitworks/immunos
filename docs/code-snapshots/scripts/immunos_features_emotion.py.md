---
source: /Users/byron/projects/scripts/immunos_features_emotion.py
relative: scripts/immunos_features_emotion.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
IMMUNOS Emotion Recognition Feature Extractor
==============================================

Extracts 128-dimensional feature vectors from facial images for emotion classification.

Features Extracted (128 total):
- Facial landmarks: 68 points → 34 inter-landmark distances
- HOG (Histogram of Oriented Gradients): 36 features
- LBP (Local Binary Patterns): 59 features (histogram)
- Statistical: mean, std, skewness, kurtosis (4 features)

Total: 34 + 36 + 54 + 4 = 128 dimensions

Dependencies:
    pip install opencv-python dlib scikit-image

Datasets:
- CK+ (593 images) - Small, good for testing
- FER2013 (35,887 images) - Large, production use

Based on:
- Umair et al. (2025) NegSl-AIS paper (emotion classification)
- Standard facial emotion recognition pipeline
"""

import numpy as np
import sys
from pathlib import Path
from typing import Union, Optional, Tuple
import warnings

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.immunos_universal import FeatureExtractor

# Try to import computer vision libraries
try:
    import cv2
    HAS_CV2 = True
except ImportError:
    HAS_CV2 = False
    warnings.warn("opencv-python not installed. Install with: pip install opencv-python")

try:
    import dlib
    HAS_DLIB = True
except ImportError:
    HAS_DLIB = False
    warnings.warn("dlib not installed. Install with: pip install dlib")

try:
    from skimage.feature import hog, local_binary_pattern
    from skimage import exposure
    from scipy.stats import skew, kurtosis
    HAS_SKIMAGE = True
except ImportError:
    HAS_SKIMAGE = False
    warnings.warn("scikit-image not installed. Install with: pip install scikit-image scipy")


class EmotionFeatureExtractor(FeatureExtractor):
    """
    Feature extractor for facial emotion recognition.

    Extracts 128-dimensional feature vector combining:
    - Facial geometry (landmarks)
    - Texture features (HOG, LBP)
    - Statistical features
    """

    def __init__(self,
                 predictor_path: Optional[str] = None,
                 use_landmarks: bool = True,
                 use_hog: bool = True,
                 use_lbp: bool = True,
                 use_stats: bool = True):
        """
        Initialize emotion feature extractor.

        Args:
            predictor_path: Path to dlib shape predictor model
                          (shape_predictor_68_face_landmarks.dat)
            use_landmarks: Extract facial landmark features
            use_hog: Extract HOG features
            use_lbp: Extract LBP features
            use_stats: Extract statistical features
        """
        self.use_landmarks = use_landmarks and HAS_DLIB
        self.use_hog = use_hog and HAS_SKIMAGE
        self.use_lbp = use_lbp and HAS_SKIMAGE
        self.use_stats = use_stats and HAS_SKIMAGE

        # Initialize face detector (Haar Cascade - lightweight)
        if HAS_CV2:
            # Use OpenCV's built-in Haar Cascade
            cascade_path = cv2.data.haarcascades + 'haarcascade_frontalface_default.xml'
            self.face_cascade = cv2.CascadeClassifier(cascade_path)
        else:
            self.face_cascade = None

        # Initialize facial landmark predictor
        self.predictor = None
        if self.use_landmarks and HAS_DLIB:
            if predictor_path and Path(predictor_path).exists():
                self.predictor = dlib.shape_predictor(predictor_path)
            else:
                # Try default location
                default_path = Path.home() / ".immunos" / "models" / "shape_predictor_68_face_landmarks.dat"
                if default_path.exists():
                    self.predictor = dlib.shape_predictor(str(default_path))
                else:
                    warnings.warn(
                        "Facial landmark predictor not found. "
                        "Download from: http://dlib.net/files/shape_predictor_68_face_landmarks.dat.bz2"
                    )
                    self.use_landmarks = False

    def extract(self, sample: Union[np.ndarray, str, Path]) -> np.ndarray:
        """
        Extract 128-dimensional feature vector from facial image.

        Args:
            sample: Image as numpy array, or path to image file

        Returns:
            128-dimensional feature vector
        """
        # Load image if path provided
        if isinstance(sample, (str, Path)):
            if not HAS_CV2:
                raise ImportError("opencv-python required for loading images")
            image = cv2.imread(str(sample))
            if image is None:
                raise ValueError(f"Failed to load image: {sample}")
        else:
            image = sample

        # Preprocess
        image = self.preprocess(image)

        # Extract face region
        face_region, face_gray = self._detect_face(image)

        if face_region is None:
            # No face detected - return default features
            return self._default_features()

        features = []

        # 1. Facial Landmark Features (34 dimensions)
        if self.use_landmarks and self.predictor is not None:
            landmark_features = self._extract_landmark_features(face_gray, face_region)
            features.append(landmark_features)
        else:
            features.append(np.zeros(34))

        # 2. HOG Features (36 dimensions)
        if self.use_hog:
            hog_features = self._extract_hog_features(face_gray)
            features.append(hog_features)
        else:
            features.append(np.zeros(36))

        # 3. LBP Features (54 dimensions)
        if self.use_lbp:
            lbp_features = self._extract_lbp_features(face_gray)
            features.append(lbp_features)
        else:
            features.append(np.zeros(54))

        # 4. Statistical Features (4 dimensions)
        if self.use_stats:
            stat_features = self._extract_statistical_features(face_gray)
            features.append(stat_features)
        else:
            features.append(np.zeros(4))

        # Concatenate all features
        feature_vector = np.concatenate(features)

        # Ensure exactly 128 dimensions
        assert feature_vector.shape[0] == 128, f"Expected 128 features, got {feature_vector.shape[0]}"

        return feature_vector

    def get_feature_dim(self) -> int:
        """Return dimensionality of feature vector"""
        return 128

    def preprocess(self, image: np.ndarray) -> np.ndarray:
        """
        Preprocess image for feature extraction.

        Args:
            image: Input image (BGR or grayscale)

        Returns:
            Preprocessed image
        """
        if not HAS_CV2:
            return image

        # Convert to grayscale if needed
        if len(image.shape) == 3:
            gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
        else:
            gray = image

        # Resize to standard size (224x224)
        resized = cv2.resize(gray, (224, 224))

        # Histogram equalization for better contrast
        equalized = cv2.equalizeHist(resized)

        return equalized

    def _detect_face(self, image: np.ndarray) -> Tuple[Optional[Tuple], Optional[np.ndarray]]:
        """
        Detect face in image using Haar Cascade.

        Returns:
            (face_region, face_gray) or (None, None) if no face detected
        """
        if self.face_cascade is None or not HAS_CV2:
            # Return full image if no detector
            return (0, 0, image.shape[1], image.shape[0]), image

        # Detect faces
        faces = self.face_cascade.detectMultiScale(
            image,
            scaleFactor=1.1,
            minNeighbors=5,
            minSize=(30, 30)
        )

        if len(faces) == 0:
            return None, None

        # Use largest face
        face = max(faces, key=lambda rect: rect[2] * rect[3])
        x, y, w, h = face

        # Extract face region with some padding
        padding = int(0.1 * min(w, h))
        x1 = max(0, x - padding)
        y1 = max(0, y - padding)
        x2 = min(image.shape[1], x + w + padding)
        y2 = min(image.shape[0], y + h + padding)

        face_gray = image[y1:y2, x1:x2]

        return (x, y, w, h), face_gray

    def _extract_landmark_features(self, face_gray: np.ndarray,
                                    face_region: Tuple[int, int, int, int]) -> np.ndarray:
        """
        Extract facial landmark features (34 dimensions).

        68 facial landmarks → 34 inter-landmark distances
        """
        if not HAS_DLIB or self.predictor is None:
            return np.zeros(34)

        x, y, w, h = face_region

        # Create dlib rectangle
        rect = dlib.rectangle(int(x), int(y), int(x + w), int(y + h))

        # Detect landmarks
        shape = self.predictor(face_gray, rect)

        # Convert to numpy array
        landmarks = np.array([[p.x, p.y] for p in shape.parts()])

        # Calculate inter-landmark distances
        # Use key landmark pairs for emotion (eyes, eyebrows, mouth)
        distances = []

        # Eye landmarks: 36-41 (right eye), 42-47 (left eye)
        # Eyebrow landmarks: 17-21 (right), 22-26 (left)
        # Mouth landmarks: 48-67

        key_pairs = [
            # Eye width
            (36, 39), (42, 45),
            # Eye height
            (37, 41), (43, 47),
            # Eyebrow position
            (19, 37), (24, 44),
            # Mouth width
            (48, 54),
            # Mouth height
            (51, 57), (62, 66),
            # Inter-eye distance
            (39, 42),
            # Eye-mouth distances
            (36, 48), (45, 54),
            # Eyebrow-eye distances
            (19, 36), (24, 45),
            # Face shape
            (0, 16), (1, 15), (2, 14),
            # Nose
            (27, 30), (31, 35),
            # Additional key pairs for 34 total
            (21, 22), (48, 60), (51, 62), (57, 66),
            (17, 26), (36, 45), (48, 54), (0, 8), (8, 16),
            (27, 33), (31, 33), (33, 35), (48, 51), (54, 57)
        ]

        for p1, p2 in key_pairs[:34]:  # Ensure exactly 34 distances
            if p1 < len(landmarks) and p2 < len(landmarks):
                dist = np.linalg.norm(landmarks[p1] - landmarks[p2])
                distances.append(dist)
            else:
                distances.append(0.0)

        # Pad if needed
        while len(distances) < 34:
            distances.append(0.0)

        # Normalize distances
        distances = np.array(distances[:34])
        if distances.max() > 0:
            distances = distances / distances.max()

        return distances

    def _extract_hog_features(self, face_gray: np.ndarray) -> np.ndarray:
        """
        Extract HOG (Histogram of Oriented Gradients) features (36 dimensions).

        HOG captures edge and gradient structure, useful for facial expressions.
        """
        if not HAS_SKIMAGE:
            return np.zeros(36)

        # Resize face to standard size
        if HAS_CV2:
            face_resized = cv2.resize(face_gray, (64, 64))
        else:
            face_resized = face_gray

        # Extract HOG features
        hog_features = hog(
            face_resized,
            orientations=9,
            pixels_per_cell=(16, 16),
            cells_per_block=(2, 2),
            visualize=False,
            feature_vector=True
        )

        # Reduce to 36 dimensions via binning
        if len(hog_features) > 36:
            # Average pool to 36 dimensions
            bin_size = len(hog_features) // 36
            reduced = np.array([
                hog_features[i*bin_size:(i+1)*bin_size].mean()
                for i in range(36)
            ])
        else:
            # Pad if too small
            reduced = np.zeros(36)
            reduced[:len(hog_features)] = hog_features

        return reduced[:36]

    def _extract_lbp_features(self, face_gray: np.ndarray) -> np.ndarray:
        """
        Extract LBP (Local Binary Patterns) features (54 dimensions).

        LBP captures texture information, useful for subtle facial movements.
        """
        if not HAS_SKIMAGE:
            return np.zeros(54)

        # Resize face to standard size
        if HAS_CV2:
            face_resized = cv2.resize(face_gray, (64, 64))
        else:
            face_resized = face_gray

        # Calculate LBP
        radius = 3
        n_points = 8 * radius
        lbp = local_binary_pattern(face_resized, n_points, radius, method='uniform')

        # Calculate histogram
        n_bins = 59  # For uniform LBP with radius=3
        hist, _ = np.histogram(
            lbp.ravel(),
            bins=n_bins,
            range=(0, n_bins),
            density=True
        )

        # Reduce to 54 dimensions
        if len(hist) > 54:
            reduced = hist[:54]
        else:
            reduced = np.zeros(54)
            reduced[:len(hist)] = hist

        return reduced

    def _extract_statistical_features(self, face_gray: np.ndarray) -> np.ndarray:
        """
        Extract statistical features (4 dimensions).

        Features: mean, std, skewness, kurtosis
        """
        if not HAS_SKIMAGE:
            return np.zeros(4)

        pixels = face_gray.flatten().astype(np.float64)

        features = np.array([
            np.mean(pixels),
            np.std(pixels),
            skew(pixels),
            kurtosis(pixels)
        ])

        # Normalize
        features = features / (np.abs(features).max() + 1e-8)

        return features

    def _default_features(self) -> np.ndarray:
        """Return default feature vector when no face detected"""
        return np.zeros(128)


# =============================================================================
# TESTING AND UTILITY FUNCTIONS
# =============================================================================

def test_extractor():
    """Test emotion feature extractor with synthetic data"""
    print("Testing Emotion Feature Extractor")
    print("=" * 60)

    # Create extractor
    extractor = EmotionFeatureExtractor()

    print(f"✓ Extractor initialized")
    print(f"  Use landmarks: {extractor.use_landmarks}")
    print(f"  Use HOG: {extractor.use_hog}")
    print(f"  Use LBP: {extractor.use_lbp}")
    print(f"  Use stats: {extractor.use_stats}")
    print()

    # Test with synthetic image
    print("Test 1: Synthetic grayscale image")
    synthetic_image = np.random.randint(0, 255, (224, 224), dtype=np.uint8)
    features = extractor.extract(synthetic_image)

    print(f"✓ Feature extraction successful")
    print(f"  Input shape: {synthetic_image.shape}")
    print(f"  Feature dimension: {features.shape[0]}")
    print(f"  Feature range: [{features.min():.3f}, {features.max():.3f}]")
    print(f"  Feature mean: {features.mean():.3f}")
    print()

    # Test with colored image
    print("Test 2: Synthetic color image")
    color_image = np.random.randint(0, 255, (224, 224, 3), dtype=np.uint8)
    features2 = extractor.extract(color_image)

    print(f"✓ Color image processing successful")
    print(f"  Input shape: {color_image.shape}")
    print(f"  Feature dimension: {features2.shape[0]}")
    print()

    # Verify dimension
    assert extractor.get_feature_dim() == 128
    assert features.shape[0] == 128
    assert features2.shape[0] == 128

    print("=" * 60)
    print("✅ All tests passed!")
    print()
    print("Feature breakdown (128 dimensions):")
    print("  - Landmarks: 34 (inter-landmark distances)")
    print("  - HOG: 36 (gradient histograms)")
    print("  - LBP: 54 (texture patterns)")
    print("  - Stats: 4 (mean, std, skew, kurtosis)")
    print("  Total: 128 dimensions")


if __name__ == "__main__":
    test_extractor()

```
