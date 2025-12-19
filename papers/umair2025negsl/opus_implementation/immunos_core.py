"""
IMMUNOS: Universal Artificial Immune System Framework
=====================================================

A biologically-inspired self vs. non-self detection system for:
- Emotion Recognition
- LLM Hallucination Detection
- Network Intrusion Detection
- Code Vulnerability Detection
- Research Paper Verification
- Logic/Reasoning Validation

Optimized for Apple M1/M2/M3 with MLX and TensorFlow Metal
"""

from __future__ import annotations
import numpy as np
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Tuple, Union
from enum import Enum
import json
from pathlib import Path
import sqlite3
from datetime import datetime

# Try to import Apple Silicon optimized libraries
try:
    import mlx.core as mx
    HAS_MLX = True
except ImportError:
    HAS_MLX = False
    
try:
    import tensorflow as tf
    HAS_TF = True
except ImportError:
    HAS_TF = False


# =============================================================================
# ENUMS AND DATA CLASSES
# =============================================================================

class CellType(Enum):
    """Types of immune cells in the system"""
    T_CELL = "t_cell"           # Negative selection detectors
    B_CELL = "b_cell"           # Clonal selection detectors
    NK_CELL = "nk_cell"         # Natural killer (immediate response)
    DENDRITIC = "dendritic"     # Context/signal integration
    MEMORY = "memory"           # Long-term memory cells


class DetectionResult(Enum):
    """Detection outcome categories"""
    SELF = "self"               # Known-good pattern
    NON_SELF = "non_self"       # Unknown/anomalous pattern
    DANGER = "danger"           # Confirmed threat
    UNCERTAIN = "uncertain"     # Insufficient evidence


@dataclass
class Detector:
    """Represents an immune system detector"""
    id: str
    cell_type: CellType
    center: np.ndarray
    radius: float
    r_self: float
    domain: str
    created_at: datetime = field(default_factory=datetime.now)
    match_count: int = 0
    confidence: float = 1.0
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def matches(self, sample: np.ndarray) -> Tuple[bool, float]:
        """Check if sample matches this detector"""
        distance = np.linalg.norm(sample - self.center)
        if distance < self.r_self:
            confidence = 1.0 - (distance / self.r_self)
            return True, confidence
        return False, 0.0
    
    def to_dict(self) -> dict:
        return {
            'id': self.id,
            'cell_type': self.cell_type.value,
            'center': self.center.tolist(),
            'radius': self.radius,
            'r_self': self.r_self,
            'domain': self.domain,
            'created_at': self.created_at.isoformat(),
            'match_count': self.match_count,
            'confidence': self.confidence,
            'metadata': self.metadata
        }


@dataclass
class DetectionResponse:
    """Response from the immune system"""
    result: DetectionResult
    confidence: float
    matched_detectors: List[str]
    danger_signal: float
    context: Dict[str, Any]
    timestamp: datetime = field(default_factory=datetime.now)
    
    def to_dict(self) -> dict:
        return {
            'result': self.result.value,
            'confidence': self.confidence,
            'matched_detectors': self.matched_detectors,
            'danger_signal': self.danger_signal,
            'context': self.context,
            'timestamp': self.timestamp.isoformat()
        }


@dataclass
class DomainConfig:
    """Configuration for a specific domain"""
    name: str
    num_detectors: int
    r_self: float
    feature_dim: int
    algorithms: List[str]
    thresholds: Dict[str, float]
    
    @classmethod
    def from_dict(cls, data: dict) -> 'DomainConfig':
        return cls(**data)


# =============================================================================
# ABSTRACT BASE CLASSES
# =============================================================================

class FeatureExtractor(ABC):
    """Abstract base class for domain-specific feature extraction"""
    
    @abstractmethod
    def extract(self, sample: Any) -> np.ndarray:
        """Extract feature vector from sample"""
        pass
    
    @abstractmethod
    def get_feature_dim(self) -> int:
        """Return dimensionality of feature vector"""
        pass
    
    @abstractmethod
    def preprocess(self, sample: Any) -> Any:
        """Preprocess sample before feature extraction"""
        pass


class ImmuneAlgorithm(ABC):
    """Abstract base class for AIS algorithms"""
    
    @abstractmethod
    def train(self, self_samples: np.ndarray, config: DomainConfig) -> List[Detector]:
        """Train detectors on self samples"""
        pass
    
    @abstractmethod
    def detect(self, sample: np.ndarray, detectors: List[Detector]) -> Tuple[bool, float]:
        """Detect if sample is non-self"""
        pass


# =============================================================================
# NEGATIVE SELECTION ALGORITHM
# =============================================================================

class NegativeSelection(ImmuneAlgorithm):
    """
    Negative Selection Algorithm (NSA)
    
    Biological analogy: T-cell maturation in thymus
    - Generate random detector candidates
    - Eliminate those that match self (negative selection)
    - Survivors detect non-self patterns
    """
    
    def __init__(self, random_seed: Optional[int] = None):
        self.rng = np.random.default_rng(random_seed)
    
    def train(self, self_samples: np.ndarray, config: DomainConfig) -> List[Detector]:
        """Train negative selection detectors"""
        detectors = []
        feature_dim = self_samples.shape[1]
        max_attempts = config.num_detectors * 100
        
        # Normalize self samples
        self_normalized = self._normalize(self_samples)
        
        attempts = 0
        while len(detectors) < config.num_detectors and attempts < max_attempts:
            # Generate random candidate detector
            candidate = self.rng.uniform(0, 1, size=feature_dim)
            
            # Check if candidate passes negative selection
            if self._passes_selection(candidate, self_normalized, config.r_self):
                # Calculate optimal radius for this detector
                radius = self._calculate_radius(candidate, self_normalized, config.r_self)
                
                detector = Detector(
                    id=f"nsa_{len(detectors):04d}",
                    cell_type=CellType.T_CELL,
                    center=candidate,
                    radius=radius,
                    r_self=config.r_self,
                    domain=config.name,
                    metadata={'algorithm': 'negative_selection'}
                )
                detectors.append(detector)
            
            attempts += 1
        
        return detectors
    
    def detect(self, sample: np.ndarray, detectors: List[Detector]) -> Tuple[bool, float]:
        """Detect if sample is non-self using trained detectors"""
        sample_norm = self._normalize(sample.reshape(1, -1))[0]
        
        max_confidence = 0.0
        is_non_self = False
        
        for detector in detectors:
            matches, confidence = detector.matches(sample_norm)
            if matches:
                is_non_self = True
                max_confidence = max(max_confidence, confidence)
        
        return is_non_self, max_confidence
    
    def _passes_selection(self, candidate: np.ndarray, 
                          self_samples: np.ndarray, r_self: float) -> bool:
        """Check if candidate detector doesn't match any self sample"""
        for self_sample in self_samples:
            distance = np.linalg.norm(candidate - self_sample)
            if distance < r_self:
                return False  # Matches self - reject
        return True  # Doesn't match self - valid detector
    
    def _calculate_radius(self, candidate: np.ndarray,
                          self_samples: np.ndarray, r_self: float) -> float:
        """Calculate optimal radius for detector"""
        min_distance = float('inf')
        for self_sample in self_samples:
            distance = np.linalg.norm(candidate - self_sample)
            min_distance = min(min_distance, distance)
        return min_distance - r_self
    
    def _normalize(self, samples: np.ndarray) -> np.ndarray:
        """Normalize samples to [0, 1] range"""
        min_vals = samples.min(axis=0)
        max_vals = samples.max(axis=0)
        range_vals = max_vals - min_vals
        range_vals[range_vals == 0] = 1  # Avoid division by zero
        return (samples - min_vals) / range_vals


# =============================================================================
# CLONAL SELECTION ALGORITHM
# =============================================================================

class ClonalSelection(ImmuneAlgorithm):
    """
    Clonal Selection Algorithm (CSA)
    
    Biological analogy: B-cell clonal expansion and affinity maturation
    - High-affinity antibodies proliferate
    - Mutations create diversity
    - Selection pressure improves detection
    """
    
    def __init__(self, 
                 clone_rate: float = 0.1,
                 mutation_rate: float = 0.05,
                 generations: int = 50,
                 random_seed: Optional[int] = None):
        self.clone_rate = clone_rate
        self.mutation_rate = mutation_rate
        self.generations = generations
        self.rng = np.random.default_rng(random_seed)
    
    def train(self, self_samples: np.ndarray, config: DomainConfig) -> List[Detector]:
        """Train using clonal selection with evolution"""
        feature_dim = self_samples.shape[1]
        
        # Initialize antibody population
        population = [self.rng.uniform(0, 1, size=feature_dim) 
                     for _ in range(config.num_detectors * 2)]
        
        for gen in range(self.generations):
            # Calculate affinity for each antibody (inverse of average distance to self)
            affinities = []
            for antibody in population:
                avg_dist = np.mean([np.linalg.norm(antibody - s) for s in self_samples])
                affinities.append(avg_dist)  # Higher distance = better detector
            
            # Select top performers
            sorted_indices = np.argsort(affinities)[::-1]
            num_select = int(len(population) * 0.5)
            selected = [population[i] for i in sorted_indices[:num_select]]
            
            # Clone and mutate
            clones = []
            for i, antibody in enumerate(selected):
                # More clones for better antibodies
                n_clones = max(1, int(self.clone_rate * len(population) * (1 - i/len(selected))))
                for _ in range(n_clones):
                    # Mutation rate inversely proportional to affinity
                    mut_rate = self.mutation_rate * (1 + i/len(selected))
                    clone = self._mutate(antibody, mut_rate)
                    clones.append(clone)
            
            # New population
            population = selected + clones
            population = population[:config.num_detectors * 2]
        
        # Create detectors from final population
        detectors = []
        for i, antibody in enumerate(population[:config.num_detectors]):
            detector = Detector(
                id=f"csa_{i:04d}",
                cell_type=CellType.B_CELL,
                center=antibody,
                radius=config.r_self,
                r_self=config.r_self,
                domain=config.name,
                confidence=affinities[sorted_indices[i]] if i < len(sorted_indices) else 0.5,
                metadata={'algorithm': 'clonal_selection', 'generation': self.generations}
            )
            detectors.append(detector)
        
        return detectors
    
    def detect(self, sample: np.ndarray, detectors: List[Detector]) -> Tuple[bool, float]:
        """Detect using trained B-cell detectors"""
        max_confidence = 0.0
        is_non_self = False
        
        for detector in detectors:
            matches, confidence = detector.matches(sample)
            if matches:
                is_non_self = True
                max_confidence = max(max_confidence, confidence)
        
        return is_non_self, max_confidence
    
    def _mutate(self, antibody: np.ndarray, mutation_rate: float) -> np.ndarray:
        """Apply mutation to antibody"""
        mutated = antibody.copy()
        for i in range(len(mutated)):
            if self.rng.random() < mutation_rate:
                mutated[i] += self.rng.normal(0, 0.1)
                mutated[i] = np.clip(mutated[i], 0, 1)
        return mutated


# =============================================================================
# DANGER THEORY
# =============================================================================

class DangerTheory:
    """
    Danger Theory Implementation
    
    Biological analogy: Danger signals from stressed/dying cells
    - Not just non-self, but context matters
    - Danger signals indicate actual threat
    - Prevents autoimmune-like false positives
    """
    
    def __init__(self, danger_threshold: float = 0.7):
        self.danger_threshold = danger_threshold
        self.danger_zones: List[np.ndarray] = []
    
    def calculate_danger_signal(self, sample: np.ndarray, 
                                 context: Dict[str, Any]) -> float:
        """Calculate danger signal from context"""
        signals = []
        
        # Entropy-based danger
        if 'entropy' in context:
            max_entropy = context.get('max_entropy', 1.0)
            signals.append(context['entropy'] / max_entropy)
        
        # Change rate danger
        if 'change_rate' in context:
            signals.append(min(context['change_rate'] / 10.0, 1.0))
        
        # Anomaly neighborhood
        if 'anomaly_neighbors' in context:
            total = context.get('total_neighbors', 1)
            signals.append(context['anomaly_neighbors'] / total)
        
        # Confidence drop
        if 'confidence_drop' in context:
            signals.append(context['confidence_drop'])
        
        # Pattern violation
        if 'pattern_violation' in context:
            signals.append(context['pattern_violation'])
        
        return np.mean(signals) if signals else 0.0
    
    def should_respond(self, is_non_self: bool, danger_signal: float) -> DetectionResult:
        """Determine response based on non-self detection and danger signal"""
        if not is_non_self:
            return DetectionResult.SELF
        
        if danger_signal > self.danger_threshold:
            return DetectionResult.DANGER
        elif danger_signal > self.danger_threshold * 0.5:
            return DetectionResult.NON_SELF
        else:
            return DetectionResult.UNCERTAIN
    
    def update_danger_zones(self, confirmed_threat: np.ndarray):
        """Add new danger zone around confirmed threat"""
        self.danger_zones.append(confirmed_threat)
    
    def in_danger_zone(self, sample: np.ndarray, radius: float = 0.5) -> bool:
        """Check if sample is in a danger zone"""
        for zone in self.danger_zones:
            if np.linalg.norm(sample - zone) < radius:
                return True
        return False


# =============================================================================
# DENDRITIC CELL ALGORITHM
# =============================================================================

class DendriticCell:
    """
    Dendritic Cell Algorithm (DCA)
    
    Biological analogy: Antigen-presenting cells integrating signals
    - Collect PAMP (pathogen), safe, and danger signals
    - Present antigens based on maturation context
    """
    
    def __init__(self, maturation_threshold: float = 100.0):
        self.maturation_threshold = maturation_threshold
        self.pamp: float = 0.0      # Pathogen-associated signals
        self.safe: float = 0.0      # Safe signals
        self.danger: float = 0.0    # Danger signals
        self.antigens: List[Any] = []
    
    def reset(self):
        """Reset cell state"""
        self.pamp = 0.0
        self.safe = 0.0
        self.danger = 0.0
        self.antigens = []
    
    def collect_signals(self, pamp: float, safe: float, 
                        danger: float, antigen: Any):
        """Collect signals from environment"""
        self.pamp += pamp
        self.safe += safe
        self.danger += danger
        self.antigens.append(antigen)
    
    def get_maturation_context(self) -> float:
        """
        Calculate maturation context
        Returns: -1 (fully safe) to +1 (fully dangerous)
        """
        total = self.pamp + self.safe + self.danger
        if total == 0:
            return 0.0
        return (self.pamp + self.danger - self.safe) / total
    
    def is_mature(self) -> bool:
        """Check if cell has reached maturation threshold"""
        return (self.pamp + self.safe + self.danger) > self.maturation_threshold
    
    def get_presented_antigens(self) -> List[Tuple[Any, float]]:
        """Get antigens with their context scores"""
        context = self.get_maturation_context()
        return [(ag, context) for ag in self.antigens]


# =============================================================================
# IMMUNE NETWORK
# =============================================================================

class ImmuneNetwork:
    """
    Immune Network Algorithm
    
    Biological analogy: Idiotypic network of antibodies
    - Antibodies recognize each other
    - Self-organizing network structure
    - Suppression and stimulation
    """
    
    def __init__(self, suppression_threshold: float = 0.5):
        self.suppression_threshold = suppression_threshold
        self.antibodies: List[np.ndarray] = []
        self.connections: Dict[int, Dict[int, float]] = {}
    
    def add_antibody(self, antibody: np.ndarray):
        """Add antibody and compute connections"""
        self.antibodies.append(antibody)
        idx = len(self.antibodies) - 1
        self.connections[idx] = {}
        
        # Compute connections to existing antibodies
        for i, existing in enumerate(self.antibodies[:-1]):
            affinity = self._compute_affinity(antibody, existing)
            if affinity > self.suppression_threshold:
                self.connections[idx][i] = affinity
                if i not in self.connections:
                    self.connections[i] = {}
                self.connections[i][idx] = affinity
    
    def get_network_response(self, antigen: np.ndarray) -> float:
        """Get combined network response to antigen"""
        if not self.antibodies:
            return 0.0
        
        responses = []
        for i, ab in enumerate(self.antibodies):
            # Direct affinity to antigen
            direct = self._compute_affinity(ab, antigen)
            
            # Suppression from connected antibodies
            suppression = 0.0
            if i in self.connections:
                suppression = sum(self.connections[i].values()) / len(self.antibodies)
            
            responses.append(max(0, direct - suppression))
        
        return np.mean(responses)
    
    def _compute_affinity(self, a: np.ndarray, b: np.ndarray) -> float:
        """Compute affinity between two vectors"""
        # Use inverse of normalized distance
        distance = np.linalg.norm(a - b)
        max_distance = np.sqrt(len(a))  # Maximum possible distance in unit hypercube
        return 1.0 - (distance / max_distance)


# =============================================================================
# IMMUNE MEMORY
# =============================================================================

class ImmuneMemory:
    """
    Immune Memory System
    
    Biological analogy: Memory T and B cells for rapid secondary response
    - Store successful detection patterns
    - Enable faster future detection
    """
    
    def __init__(self, db_path: str = "immunos_memory.db"):
        self.db_path = db_path
        self._init_db()
    
    def _init_db(self):
        """Initialize SQLite database"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS memory_cells (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                domain TEXT NOT NULL,
                pattern_hash TEXT NOT NULL,
                pattern_vector BLOB NOT NULL,
                detection_type TEXT NOT NULL,
                confidence REAL NOT NULL,
                hit_count INTEGER DEFAULT 1,
                first_seen TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                last_seen TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                metadata TEXT
            )
        ''')
        
        cursor.execute('''
            CREATE INDEX IF NOT EXISTS idx_domain_hash 
            ON memory_cells(domain, pattern_hash)
        ''')
        
        conn.commit()
        conn.close()
    
    def store(self, domain: str, pattern: np.ndarray, 
              detection_type: str, confidence: float,
              metadata: Optional[Dict] = None):
        """Store pattern in immune memory"""
        pattern_hash = self._hash_pattern(pattern)
        
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        # Check if pattern already exists
        cursor.execute('''
            SELECT id, hit_count FROM memory_cells 
            WHERE domain = ? AND pattern_hash = ?
        ''', (domain, pattern_hash))
        
        result = cursor.fetchone()
        if result:
            # Update existing
            cursor.execute('''
                UPDATE memory_cells 
                SET hit_count = hit_count + 1, last_seen = CURRENT_TIMESTAMP
                WHERE id = ?
            ''', (result[0],))
        else:
            # Insert new
            cursor.execute('''
                INSERT INTO memory_cells 
                (domain, pattern_hash, pattern_vector, detection_type, confidence, metadata)
                VALUES (?, ?, ?, ?, ?, ?)
            ''', (domain, pattern_hash, pattern.tobytes(), detection_type, 
                  confidence, json.dumps(metadata or {})))
        
        conn.commit()
        conn.close()
    
    def recall(self, domain: str, pattern: np.ndarray, 
               threshold: float = 0.9) -> Optional[Dict]:
        """Recall similar pattern from memory"""
        pattern_hash = self._hash_pattern(pattern)
        
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        # Exact match first
        cursor.execute('''
            SELECT detection_type, confidence, hit_count, metadata
            FROM memory_cells 
            WHERE domain = ? AND pattern_hash = ?
        ''', (domain, pattern_hash))
        
        result = cursor.fetchone()
        conn.close()
        
        if result:
            return {
                'detection_type': result[0],
                'confidence': result[1],
                'hit_count': result[2],
                'metadata': json.loads(result[3])
            }
        return None
    
    def _hash_pattern(self, pattern: np.ndarray) -> str:
        """Create hash of pattern for quick lookup"""
        import hashlib
        return hashlib.sha256(pattern.tobytes()).hexdigest()[:32]


# =============================================================================
# UNIFIED IMMUNE SYSTEM
# =============================================================================

class ImmuneSystem:
    """
    Main IMMUNOS System
    
    Combines all AIS algorithms into unified framework
    """
    
    def __init__(self, config_path: Optional[str] = None):
        self.domains: Dict[str, Dict] = {}
        self.algorithms = {
            'negative_selection': NegativeSelection(),
            'clonal_selection': ClonalSelection()
        }
        self.danger_theory = DangerTheory()
        self.immune_network = ImmuneNetwork()
        self.memory = ImmuneMemory()
        self.dendritic_cells: List[DendriticCell] = []
    
    def add_domain(self, name: str, config: DomainConfig, 
                   feature_extractor: FeatureExtractor):
        """Add a detection domain"""
        self.domains[name] = {
            'config': config,
            'feature_extractor': feature_extractor,
            'detectors': []
        }
    
    def train(self, domain: str, self_samples: List[Any],
              algorithms: Optional[List[str]] = None) -> int:
        """Train detectors for a domain"""
        if domain not in self.domains:
            raise ValueError(f"Unknown domain: {domain}")
        
        domain_data = self.domains[domain]
        config = domain_data['config']
        extractor = domain_data['feature_extractor']
        
        # Extract features from self samples
        features = np.array([extractor.extract(s) for s in self_samples])
        
        # Use specified algorithms or all from config
        algo_names = algorithms or config.algorithms
        
        all_detectors = []
        for algo_name in algo_names:
            if algo_name in self.algorithms:
                detectors = self.algorithms[algo_name].train(features, config)
                all_detectors.extend(detectors)
        
        domain_data['detectors'] = all_detectors
        return len(all_detectors)
    
    def detect(self, domain: str, sample: Any, 
               context: Optional[Dict[str, Any]] = None) -> DetectionResponse:
        """Detect if sample is self or non-self"""
        if domain not in self.domains:
            raise ValueError(f"Unknown domain: {domain}")
        
        domain_data = self.domains[domain]
        extractor = domain_data['feature_extractor']
        detectors = domain_data['detectors']
        
        # Extract features
        features = extractor.extract(sample)
        
        # Check immune memory first
        memory_result = self.memory.recall(domain, features)
        if memory_result:
            return DetectionResponse(
                result=DetectionResult(memory_result['detection_type']),
                confidence=memory_result['confidence'],
                matched_detectors=['memory'],
                danger_signal=0.0,
                context={'source': 'immune_memory', 'hit_count': memory_result['hit_count']}
            )
        
        # Run detection algorithms
        is_non_self = False
        max_confidence = 0.0
        matched = []
        
        for detector in detectors:
            matches, confidence = detector.matches(features)
            if matches:
                is_non_self = True
                matched.append(detector.id)
                max_confidence = max(max_confidence, confidence)
        
        # Calculate danger signal
        context = context or {}
        danger_signal = self.danger_theory.calculate_danger_signal(features, context)
        
        # Determine final result
        if is_non_self:
            result = self.danger_theory.should_respond(True, danger_signal)
        else:
            result = DetectionResult.SELF
        
        # Store in memory if non-self
        if is_non_self:
            self.memory.store(domain, features, result.value, max_confidence)
        
        return DetectionResponse(
            result=result,
            confidence=max_confidence,
            matched_detectors=matched,
            danger_signal=danger_signal,
            context=context
        )
    
    def get_statistics(self, domain: str) -> Dict[str, Any]:
        """Get statistics for a domain"""
        if domain not in self.domains:
            raise ValueError(f"Unknown domain: {domain}")
        
        domain_data = self.domains[domain]
        detectors = domain_data['detectors']
        
        return {
            'domain': domain,
            'num_detectors': len(detectors),
            'detector_types': {ct.value: sum(1 for d in detectors if d.cell_type == ct) 
                             for ct in CellType},
            'avg_confidence': np.mean([d.confidence for d in detectors]) if detectors else 0,
            'total_matches': sum(d.match_count for d in detectors)
        }
    
    def save_detectors(self, domain: str, path: str):
        """Save trained detectors to file"""
        if domain not in self.domains:
            raise ValueError(f"Unknown domain: {domain}")
        
        detectors = self.domains[domain]['detectors']
        data = [d.to_dict() for d in detectors]
        
        with open(path, 'w') as f:
            json.dump(data, f, indent=2)
    
    def load_detectors(self, domain: str, path: str):
        """Load detectors from file"""
        if domain not in self.domains:
            raise ValueError(f"Unknown domain: {domain}")
        
        with open(path, 'r') as f:
            data = json.load(f)
        
        detectors = []
        for d in data:
            detector = Detector(
                id=d['id'],
                cell_type=CellType(d['cell_type']),
                center=np.array(d['center']),
                radius=d['radius'],
                r_self=d['r_self'],
                domain=d['domain'],
                created_at=datetime.fromisoformat(d['created_at']),
                match_count=d['match_count'],
                confidence=d['confidence'],
                metadata=d['metadata']
            )
            detectors.append(detector)
        
        self.domains[domain]['detectors'] = detectors


# =============================================================================
# DEFAULT DOMAIN CONFIGURATIONS
# =============================================================================

DEFAULT_CONFIGS = {
    'emotion': DomainConfig(
        name='emotion',
        num_detectors=25,
        r_self=0.9,
        feature_dim=128,
        algorithms=['negative_selection', 'clonal_selection'],
        thresholds={'danger': 0.7, 'uncertainty': 0.3}
    ),
    'hallucination': DomainConfig(
        name='hallucination',
        num_detectors=30,
        r_self=0.85,
        feature_dim=256,
        algorithms=['negative_selection', 'clonal_selection'],
        thresholds={'danger': 0.75, 'uncertainty': 0.35}
    ),
    'network': DomainConfig(
        name='network',
        num_detectors=20,
        r_self=0.95,
        feature_dim=41,
        algorithms=['negative_selection'],
        thresholds={'danger': 0.8, 'uncertainty': 0.2}
    ),
    'code': DomainConfig(
        name='code',
        num_detectors=35,
        r_self=0.80,
        feature_dim=64,
        algorithms=['negative_selection', 'clonal_selection'],
        thresholds={'danger': 0.65, 'uncertainty': 0.4}
    ),
    'research': DomainConfig(
        name='research',
        num_detectors=20,
        r_self=0.90,
        feature_dim=128,
        algorithms=['negative_selection'],
        thresholds={'danger': 0.7, 'uncertainty': 0.3}
    )
}


# =============================================================================
# EXAMPLE USAGE
# =============================================================================

if __name__ == "__main__":
    # Example: Simple test with synthetic data
    print("IMMUNOS Framework Test")
    print("=" * 50)
    
    # Create system
    immune = ImmuneSystem()
    
    # Create simple feature extractor for testing
    class TestFeatureExtractor(FeatureExtractor):
        def extract(self, sample):
            return np.array(sample)
        def get_feature_dim(self):
            return 10
        def preprocess(self, sample):
            return sample
    
    # Add test domain
    config = DomainConfig(
        name='test',
        num_detectors=10,
        r_self=0.5,
        feature_dim=10,
        algorithms=['negative_selection'],
        thresholds={'danger': 0.7}
    )
    immune.add_domain('test', config, TestFeatureExtractor())
    
    # Generate synthetic self samples (normal patterns)
    np.random.seed(42)
    self_samples = [np.random.uniform(0.3, 0.7, 10) for _ in range(100)]
    
    # Train
    num_detectors = immune.train('test', self_samples)
    print(f"Trained {num_detectors} detectors")
    
    # Test with self sample (should be detected as SELF)
    test_self = np.random.uniform(0.3, 0.7, 10)
    result = immune.detect('test', test_self)
    print(f"Self sample: {result.result.value} (confidence: {result.confidence:.3f})")
    
    # Test with non-self sample (should be detected as NON_SELF)
    test_non_self = np.random.uniform(0.8, 1.0, 10)  # Outside self range
    result = immune.detect('test', test_non_self, context={'entropy': 0.9, 'max_entropy': 1.0})
    print(f"Non-self sample: {result.result.value} (confidence: {result.confidence:.3f})")
    
    # Get statistics
    stats = immune.get_statistics('test')
    print(f"Statistics: {stats}")
