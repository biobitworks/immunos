# IMMUNOS: Universal Artificial Immune System Framework

## Claude Code System Prompt

You are developing IMMUNOS, a universal Artificial Immune System (AIS) framework implementing the biological self vs. non-self paradigm. This system must be:

1. **Domain-Agnostic**: Work across multiple use cases (emotion recognition, LLM hallucination detection, network security, code vulnerability detection, research paper verification, logic validation)
2. **Portable**: Run from USB drive, VM, or local installation like a biological immune system
3. **Trainable on M1 MacBook Pro**: Optimized for 16GB RAM, 1TB storage, with optional cloud GPU support
4. **Biologically-Inspired**: Implement multiple AIS algorithms that mirror human immune function

---

## Core Architecture: The Immune System Metaphor

```
┌─────────────────────────────────────────────────────────────────────────┐
│                          IMMUNOS FRAMEWORK                               │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────┐               │
│  │   THYMUS     │    │  BONE MARROW │    │   LYMPH      │               │
│  │  (Training)  │    │  (Detectors) │    │   NODES      │               │
│  │              │    │              │    │  (Analysis)  │               │
│  │ - Self       │    │ - T-Cells    │    │              │               │
│  │   Samples    │    │ - B-Cells    │    │ - Pattern    │               │
│  │ - Negative   │    │ - NK Cells   │    │   Matching   │               │
│  │   Selection  │    │ - Memory     │    │ - Affinity   │               │
│  │ - Tolerance  │    │   Cells      │    │   Scoring    │               │
│  └──────────────┘    └──────────────┘    └──────────────┘               │
│         │                   │                   │                        │
│         └───────────────────┼───────────────────┘                        │
│                             ▼                                            │
│  ┌─────────────────────────────────────────────────────────────────┐    │
│  │                    BLOOD STREAM (Runtime)                        │    │
│  │                                                                  │    │
│  │   Input → Feature Extraction → Detector Matching → Response     │    │
│  └─────────────────────────────────────────────────────────────────┘    │
│                             │                                            │
│         ┌───────────────────┼───────────────────┐                        │
│         ▼                   ▼                   ▼                        │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────┐               │
│  │   DANGER     │    │   DENDRITIC  │    │   IMMUNE     │               │
│  │   SIGNALS    │    │    CELLS     │    │   MEMORY     │               │
│  │              │    │              │    │              │               │
│  │ - Anomaly    │    │ - Context    │    │ - Known      │               │
│  │   Thresholds │    │   Signals    │    │   Patterns   │               │
│  │ - Stress     │    │ - Antigen    │    │ - Historical │               │
│  │   Indicators │    │   Presenting │    │   Responses  │               │
│  └──────────────┘    └──────────────┘    └──────────────┘               │
│                                                                          │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## Algorithm Implementations

### 1. Negative Selection Algorithm (NSA) - Primary Detection

```python
"""
Negative Selection: Train on SELF, detect NON-SELF
- Generate random detector candidates
- Eliminate any that match self-samples (negative selection)
- Remaining detectors identify anything NOT in training set
"""

class NegativeSelectionDetector:
    def __init__(self, r_self: float = 0.9, num_detectors: int = 20):
        self.r_self = r_self  # Detection radius threshold
        self.num_detectors = num_detectors
        self.detectors = []
        self.self_samples = []
    
    def train(self, self_samples: np.ndarray):
        """Train on known-good (self) samples"""
        self.self_samples = self_samples
        candidates = self._generate_candidates(n=self.num_detectors * 10)
        
        for candidate in candidates:
            if self._passes_negative_selection(candidate):
                self.detectors.append(candidate)
                if len(self.detectors) >= self.num_detectors:
                    break
    
    def _passes_negative_selection(self, candidate) -> bool:
        """Detector is valid if it doesn't match any self-sample"""
        for self_sample in self.self_samples:
            distance = np.linalg.norm(candidate - self_sample)
            if distance < self.r_self:
                return False  # Too close to self - reject
        return True  # Doesn't match self - valid detector
    
    def detect(self, sample) -> tuple[bool, float]:
        """Returns (is_anomaly, confidence)"""
        for detector in self.detectors:
            distance = np.linalg.norm(sample - detector)
            if distance < self.r_self:
                return True, 1.0 - (distance / self.r_self)
        return False, 0.0
```

### 2. Clonal Selection Algorithm (CSA) - Adaptive Learning

```python
"""
Clonal Selection: Detectors that succeed get cloned and mutated
- High-affinity detectors proliferate
- Mutations create diversity
- Affinity maturation improves detection over time
"""

class ClonalSelectionOptimizer:
    def __init__(self, population_size: int = 100, clone_rate: float = 0.1):
        self.population_size = population_size
        self.clone_rate = clone_rate
        self.antibodies = []
    
    def evolve(self, antigens: np.ndarray, generations: int = 50):
        """Evolve antibody population to match antigens"""
        self._initialize_population()
        
        for gen in range(generations):
            # Calculate affinity for each antibody
            affinities = [self._affinity(ab, antigens) for ab in self.antibodies]
            
            # Select top performers
            sorted_indices = np.argsort(affinities)[::-1]
            selected = [self.antibodies[i] for i in sorted_indices[:int(self.population_size * 0.5)]]
            
            # Clone and mutate based on affinity
            clones = []
            for i, ab in enumerate(selected):
                n_clones = int(self.clone_rate * self.population_size * (1 - i/len(selected)))
                for _ in range(n_clones):
                    mutation_rate = 1.0 / (affinities[sorted_indices[i]] + 1)
                    clones.append(self._mutate(ab, mutation_rate))
            
            self.antibodies = selected + clones
```

### 3. Danger Theory Algorithm (DT) - Context-Aware Detection

```python
"""
Danger Theory: Response triggered by danger signals, not just non-self
- Danger signals indicate tissue damage/stress
- Context matters: same antigen may be safe or dangerous
- Prevents autoimmune-like false positives
"""

class DangerTheoryDetector:
    def __init__(self, danger_threshold: float = 0.7):
        self.danger_threshold = danger_threshold
        self.danger_zones = []
    
    def calculate_danger_signal(self, sample, context: dict) -> float:
        """
        Calculate danger signal based on context
        - High entropy = potential danger
        - Rapid changes = potential danger
        - Unusual patterns = potential danger
        """
        signals = []
        
        if 'entropy' in context:
            signals.append(context['entropy'] / context.get('max_entropy', 1.0))
        
        if 'change_rate' in context:
            signals.append(min(context['change_rate'] / 10.0, 1.0))
        
        if 'anomaly_neighbors' in context:
            signals.append(context['anomaly_neighbors'] / context.get('total_neighbors', 1))
        
        return np.mean(signals) if signals else 0.0
    
    def should_respond(self, sample, context: dict) -> bool:
        """Only respond if both non-self AND danger signal present"""
        danger = self.calculate_danger_signal(sample, context)
        return danger > self.danger_threshold
```

### 4. Dendritic Cell Algorithm (DCA) - Signal Integration

```python
"""
Dendritic Cell Algorithm: Multi-signal integration for decision making
- PAMP signals: Pathogen-associated patterns (definite threat indicators)
- Safe signals: Normal operation indicators
- Danger signals: Stress/damage indicators
"""

class DendriticCell:
    def __init__(self):
        self.pamp = 0.0      # Pathogen signals
        self.safe = 0.0      # Safe signals
        self.danger = 0.0    # Danger signals
        self.antigens = []   # Collected antigens
        self.maturation_threshold = 100.0
    
    def collect_signals(self, pamp: float, safe: float, danger: float, antigen):
        """Collect signals from environment"""
        self.pamp += pamp
        self.safe += safe
        self.danger += danger
        self.antigens.append(antigen)
    
    def get_context(self) -> float:
        """Calculate maturation context: -1 (safe) to +1 (dangerous)"""
        total = self.pamp + self.safe + self.danger
        if total == 0:
            return 0.0
        return (self.pamp + self.danger - self.safe) / total
    
    def is_mature(self) -> bool:
        return (self.pamp + self.safe + self.danger) > self.maturation_threshold
```

### 5. Immune Network Algorithm (INA) - Pattern Relationships

```python
"""
Immune Network: Antibodies recognize each other, forming a network
- Idiotypic network theory
- Clustering and pattern relationships
- Self-organizing detection
"""

class ImmuneNetwork:
    def __init__(self, suppression_threshold: float = 0.5):
        self.antibodies = []
        self.connections = {}  # Adjacency matrix
        self.suppression_threshold = suppression_threshold
    
    def add_antibody(self, antibody):
        """Add antibody and compute connections to existing network"""
        self.antibodies.append(antibody)
        idx = len(self.antibodies) - 1
        self.connections[idx] = {}
        
        for i, existing in enumerate(self.antibodies[:-1]):
            affinity = self._compute_affinity(antibody, existing)
            if affinity > self.suppression_threshold:
                self.connections[idx][i] = affinity
                self.connections[i][idx] = affinity
    
    def get_network_response(self, antigen) -> float:
        """Get combined network response to antigen"""
        responses = []
        for i, ab in enumerate(self.antibodies):
            direct = self._compute_affinity(ab, antigen)
            # Suppression from connected antibodies
            suppression = sum(self.connections.get(i, {}).values()) / len(self.antibodies)
            responses.append(max(0, direct - suppression))
        return np.mean(responses)
```

---

## Domain-Specific Feature Extractors

### Emotion Recognition Features
```python
class EmotionFeatureExtractor:
    """Extract features from facial images for emotion recognition"""
    
    def extract(self, image: np.ndarray) -> np.ndarray:
        features = []
        
        # Facial landmarks (68 points)
        landmarks = self._detect_landmarks(image)
        features.extend(self._landmark_distances(landmarks))
        
        # Histogram of Oriented Gradients
        hog_features = self._compute_hog(image)
        features.extend(hog_features)
        
        # Local Binary Patterns
        lbp_features = self._compute_lbp(image)
        features.extend(lbp_features)
        
        # Statistical features
        features.extend([
            np.mean(image), np.std(image),
            skew(image.flatten()), kurtosis(image.flatten())
        ])
        
        return np.array(features)
```

### LLM Hallucination Features
```python
class HallucinationFeatureExtractor:
    """Extract features to detect LLM hallucinations"""
    
    def extract(self, response: str, context: str = None) -> np.ndarray:
        features = []
        
        # Semantic consistency
        if context:
            features.append(self._semantic_similarity(response, context))
        
        # Confidence indicators
        features.append(self._hedging_ratio(response))
        features.append(self._certainty_score(response))
        
        # Factual density
        features.append(self._entity_density(response))
        features.append(self._claim_count(response))
        
        # Self-consistency (if multiple samples)
        features.append(self._entropy_score(response))
        
        # Knowledge triplet verification
        triplets = self._extract_triplets(response)
        features.append(len(triplets))
        
        return np.array(features)
```

### Network Traffic Features
```python
class NetworkFeatureExtractor:
    """Extract features from network flows for intrusion detection"""
    
    def extract(self, flow: dict) -> np.ndarray:
        features = [
            flow.get('duration', 0),
            flow.get('protocol_type', 0),
            flow.get('service', 0),
            flow.get('src_bytes', 0),
            flow.get('dst_bytes', 0),
            flow.get('flag', 0),
            flow.get('land', 0),
            flow.get('wrong_fragment', 0),
            flow.get('urgent', 0),
            flow.get('count', 0),
            flow.get('srv_count', 0),
            flow.get('serror_rate', 0),
            flow.get('srv_serror_rate', 0),
            flow.get('rerror_rate', 0),
            # ... 41 features from NSL-KDD
        ]
        return np.array(features)
```

### Code Vulnerability Features
```python
class CodeFeatureExtractor:
    """Extract features from source code for vulnerability detection"""
    
    def extract(self, code: str) -> np.ndarray:
        features = []
        
        # Structural features
        ast = self._parse_ast(code)
        features.append(self._cyclomatic_complexity(ast))
        features.append(self._ast_depth(ast))
        features.append(self._function_count(ast))
        
        # Security patterns
        features.append(self._dangerous_function_count(code))
        features.append(self._buffer_operation_count(code))
        features.append(self._user_input_handling(code))
        
        # Data flow features
        features.append(self._taint_propagation_score(ast))
        features.append(self._sanitization_ratio(code))
        
        # Statistical features
        features.append(len(code))
        features.append(self._entropy(code))
        
        return np.array(features)
```

### Research Paper Features
```python
class ResearchPaperFeatureExtractor:
    """Extract features for scientific claim verification"""
    
    def extract(self, claim: str, evidence: str = None) -> np.ndarray:
        features = []
        
        # Claim analysis
        features.append(self._claim_specificity(claim))
        features.append(self._quantitative_claims(claim))
        features.append(self._causal_language(claim))
        
        # Evidence matching (if provided)
        if evidence:
            features.append(self._semantic_entailment(claim, evidence))
            features.append(self._contradiction_score(claim, evidence))
            features.append(self._evidence_strength(evidence))
        
        # Citation patterns
        features.append(self._citation_density(claim))
        
        return np.array(features)
```

---

## Optimal Parameters by Domain

Based on research benchmarks:

| Domain | Detectors | R_self | Expected Accuracy |
|--------|-----------|--------|-------------------|
| Emotion (FER2013) | 25 | 0.85-0.95 | 73-87% |
| Hallucination (HaluEval) | 30 | 0.80-0.90 | 80-90% |
| Network (NSL-KDD) | 20 | 0.90-1.00 | 93-99% |
| Code (DiverseVul) | 35 | 0.75-0.85 | 70-89% |
| Research (SciFact) | 20 | 0.85-0.95 | 80-88% |

---

## Recommended Datasets (M1 MacBook Compatible)

### Small (Quick Testing - Minutes)
| Domain | Dataset | Size | Download |
|--------|---------|------|----------|
| Emotion | CK+ | 593 images | Direct |
| Hallucination | TruthfulQA | 817 questions | HuggingFace |
| Network | NSL-KDD | 125K records | UNB |
| Code | Devign | 27K functions | GitHub |
| Research | SciFact | 1.4K claims | AllenAI |

### Medium (Full Training - Hours)
| Domain | Dataset | Size | Download |
|--------|---------|------|----------|
| Emotion | FER2013 | 35K images | Kaggle |
| Hallucination | HaluEval | 35K samples | GitHub |
| Network | CICIDS2017 | 2.8M flows | UNB |
| Code | BigVul | 188K functions | GitHub |
| Research | SciFact-Open | 500K abstracts | AllenAI |

### Large (Cloud GPU Recommended)
| Domain | Dataset | Size | Notes |
|--------|---------|------|-------|
| Emotion | AffectNet | 1M images | Academic license |
| Network | UNSW-NB15 | 2.5M flows | Complex features |
| Code | MegaVul | 17K vulnerabilities | Multi-representation |

---

## Installation Options

### 1. Portable USB (Tails-like)
```bash
# Create bootable IMMUNOS USB
# Base: Ubuntu 24.04 LTS Minimal

# Install core dependencies
apt install python3.11 python3-pip sqlite3 git

# Install ML frameworks (Apple Silicon optimized)
pip install mlx tensorflow-metal torch torchvision

# Install IMMUNOS
git clone https://github.com/your-org/immunos.git
cd immunos && pip install -e .

# Create persistent encrypted storage
cryptsetup luksFormat /dev/sdb2
cryptsetup open /dev/sdb2 immunos-data
mkfs.ext4 /dev/mapper/immunos-data
```

### 2. Virtual Machine
```bash
# Recommended: Ubuntu 22.04/24.04 on VirtualBox/UTM
# Allocate: 8GB RAM, 50GB disk, shared clipboard

# For M1: Use UTM with ARM Ubuntu
# For Intel: VirtualBox with Nested VT-x
```

### 3. Local Installation (macOS)
```bash
# Install Homebrew
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install Python and dependencies
brew install python@3.11 miniforge

# Create conda environment
conda create -n immunos python=3.11
conda activate immunos

# Install Apple-optimized ML
pip install mlx tensorflow-macos tensorflow-metal
pip install torch torchvision torchaudio

# Install IMMUNOS
pip install immunos
```

---

## Project Structure

```
immunos/
├── core/
│   ├── __init__.py
│   ├── algorithms/
│   │   ├── negative_selection.py    # NSA implementation
│   │   ├── clonal_selection.py      # CSA implementation
│   │   ├── danger_theory.py         # DT implementation
│   │   ├── dendritic_cell.py        # DCA implementation
│   │   └── immune_network.py        # INA implementation
│   ├── detectors/
│   │   ├── base.py                  # Abstract detector class
│   │   ├── t_cell.py                # T-cell detectors
│   │   ├── b_cell.py                # B-cell detectors
│   │   └── memory_cell.py           # Memory cells
│   └── features/
│       ├── base.py                  # Abstract feature extractor
│       ├── emotion.py               # Emotion recognition
│       ├── hallucination.py         # LLM hallucination
│       ├── network.py               # Network traffic
│       ├── code.py                  # Code vulnerability
│       └── research.py              # Research verification
├── domains/
│   ├── emotion_recognition/
│   │   ├── config.yaml
│   │   ├── train.py
│   │   └── infer.py
│   ├── hallucination_detection/
│   │   ├── config.yaml
│   │   ├── train.py
│   │   └── infer.py
│   ├── network_security/
│   │   ├── config.yaml
│   │   ├── train.py
│   │   └── infer.py
│   ├── code_security/
│   │   ├── config.yaml
│   │   ├── train.py
│   │   └── infer.py
│   └── research_verification/
│       ├── config.yaml
│       ├── train.py
│       └── infer.py
├── dashboard/
│   ├── app.py                       # Flask application
│   ├── websocket.py                 # Real-time updates
│   ├── templates/
│   └── static/
├── data/
│   ├── self_samples/                # Known-good samples
│   ├── detectors/                   # Trained detectors
│   └── memory/                      # Immune memory
├── config/
│   ├── default.yaml
│   └── domains/
├── scripts/
│   ├── create_usb.sh               # Create bootable USB
│   ├── setup_vm.sh                 # Setup VM
│   └── train_all.sh                # Train all domains
├── tests/
└── docs/
```

---

## Training Commands

```bash
# Train emotion recognition
immunos train emotion --dataset fer2013 --detectors 25 --r-self 0.9

# Train hallucination detection
immunos train hallucination --dataset halueval --detectors 30 --r-self 0.85

# Train network security
immunos train network --dataset nsl-kdd --detectors 20 --r-self 0.95

# Train code security
immunos train code --dataset devign --detectors 35 --r-self 0.80

# Train research verification
immunos train research --dataset scifact --detectors 20 --r-self 0.90

# Train all domains
immunos train all --quick  # Use small datasets
```

---

## API Usage

```python
from immunos import ImmuneSystem
from immunos.domains import EmotionRecognition, HallucinationDetection

# Initialize immune system
immune = ImmuneSystem()

# Add domain modules
immune.add_domain('emotion', EmotionRecognition())
immune.add_domain('hallucination', HallucinationDetection())

# Train on self samples
immune.train('emotion', self_samples=emotion_dataset)

# Detect anomalies
result = immune.detect('emotion', test_image)
print(f"Anomaly: {result.is_anomaly}, Confidence: {result.confidence}")
print(f"Matched Detector: {result.detector_id}")

# Get immune memory
memory = immune.get_memory('emotion')
```

---

## Dashboard Integration

```python
# Run IMMUNOS dashboard
immunos dashboard --port 5000

# Features:
# - Real-time detection monitoring
# - Detector visualization
# - Training progress
# - Immune memory explorer
# - Multi-domain switching
```

---

## Cloud GPU Training (Optional)

For larger datasets, use cloud GPU:

```bash
# Google Colab
!pip install immunos[cloud]
from immunos.cloud import CloudTrainer

trainer = CloudTrainer(provider='colab')
trainer.train('emotion', dataset='affectnet', epochs=50)

# Download trained detectors
trainer.download_detectors('emotion', path='./detectors/')
```

---

## Key Implementation Notes

1. **Self vs Non-Self Philosophy**: Always train on KNOWN-GOOD data. The system learns what's normal and flags everything else.

2. **Multi-Algorithm Fusion**: Combine multiple AIS algorithms for better accuracy:
   - NSA for initial detection
   - CSA for adaptive improvement
   - DT for context-aware filtering
   - DCA for signal integration

3. **Feature Engineering**: Domain-specific feature extraction is critical. Generic features won't work.

4. **Threshold Tuning**: R_self threshold controls sensitivity vs. specificity. Start high (0.95) and tune down.

5. **Memory Management**: For M1 with 16GB, batch processing is essential for larger datasets.

---

## References

- Umair et al. (2025): NegSl-AIS for emotion classification
- Forrest et al. (1994): Original NSA algorithm
- de Castro & Von Zuben (2002): CLONALG algorithm
- Matzinger (2002): Danger Theory
- Greensmith et al. (2005): Dendritic Cell Algorithm
- Wadden et al. (2020): SciFact dataset
- Li et al. (2023): HaluEval benchmark
- Sharafaldin et al. (2018): CICIDS2017 dataset
