# IMMUNOS Architecture Based on NegSl-AIS

## Biological Immune System → AI Security Mapping

This document maps the NegSl-AIS biological immune system concepts to your IMMUNOS security dashboard for AI-assisted development.

---

## 1. CORE METAPHOR MAPPING

### Your Original IMMUNOS Metaphor (from blog)

| Immune Component | Security Function |
|------------------|------------------|
| NK Cells | Static Scanners |
| T Cells | Pattern Memory |
| B Cells | Verification/Antibodies |
| Dendritic Cells | Reporting/Alerting |

### Enhanced with NegSl-AIS Concepts

| Biological Component | NegSl-AIS Role | IMMUNOS Implementation |
|---------------------|----------------|----------------------|
| **Thymus** | Training Environment | Signature Training Pipeline |
| **Self Antigens** | Known-good samples | Safe code pattern database |
| **Non-Self Antigens** | Threat samples | Malicious/vulnerable patterns |
| **Naïve T-Cells** | Candidate detectors | Untrained pattern signatures |
| **Negative Selection** | Detector validation | Signature validation process |
| **Immunocompetent T-Cells** | Valid detectors | Verified threat detectors |
| **Blood Stream** | Production environment | Live code scanning |
| **R_self Threshold** | Binding affinity | Detection sensitivity |

---

## 2. NEGATIVE SELECTION FOR SECURITY

### The Key Insight

Traditional security: "Learn what's bad, flag the bad"
NegSl-AIS approach: "Learn what's good, flag everything else"

This is powerful because:
1. Malware constantly mutates (like pathogens)
2. Good code patterns are more stable
3. Novel attacks are automatically flagged as "non-self"

### Implementation Strategy

```
IMMUNOS Training Phase:
1. Collect known-safe code patterns (self-samples)
2. Extract features from safe code
3. Generate random detector patterns
4. DISCARD detectors that match safe code
5. Keep detectors that don't match safe code
6. These detectors will flag unknown/suspicious code

IMMUNOS Detection Phase:
1. New code submitted for scanning
2. Extract features from new code
3. Compare against valid detectors
4. If ANY detector matches → Flag as suspicious
5. If NO detectors match → Code resembles "self" (safe)
```

---

## 3. FEATURE EXTRACTION FOR CODE ANALYSIS

### Adapted from NegSl-AIS Features

| Original Feature | Code Analysis Equivalent |
|------------------|-------------------------|
| Waveform Length | Code complexity metric |
| Energy | Instruction density |
| Entropy | Code randomness/obfuscation |
| PSD | Frequency of operations |
| Zero Crossing Rate | Control flow changes |
| Skewness | Distribution of instructions |
| Kurtosis | Outlier patterns |

### Suggested Code Features

```python
CODE_FEATURES = {
    # Structural Features
    'cyclomatic_complexity': "Number of independent paths",
    'function_count': "Number of functions",
    'import_count': "Number of imports",
    'line_count': "Lines of code",
    
    # Pattern Features
    'api_call_frequency': "Frequency of API calls",
    'string_entropy': "Entropy of string literals",
    'control_flow_complexity': "Nested loops/conditions",
    
    # Security-Specific Features
    'dangerous_function_count': "eval, exec, os.system, etc.",
    'network_operation_count': "HTTP, socket operations",
    'file_operation_count': "File I/O operations",
    'injection_pattern_score': "SQL/command injection indicators",
    
    # AST Features
    'ast_depth': "Maximum AST depth",
    'node_type_distribution': "Distribution of AST node types",
    'identifier_entropy': "Entropy of variable names"
}
```

---

## 4. SCANNER FUSION (Modality Biasing)

### Applying the Paper's Approach

Just as NegSl-AIS weights different physiological modalities, IMMUNOS can weight different scanners:

```python
SCANNER_WEIGHTS = {
    # Based on scanner reliability/accuracy
    'bandit': 0.25,        # High accuracy for Python security
    'semgrep': 0.30,       # Comprehensive pattern matching
    'codeql': 0.25,        # Deep semantic analysis
    'custom_ai': 0.20,     # AI-based detection
    
    # Sum = 1.0
}

def fuse_scanner_results(results: Dict[str, float]) -> float:
    """
    Fuse results from multiple scanners using modality biasing.
    
    Args:
        results: Dict mapping scanner name to threat score (0-1)
        
    Returns:
        Fused threat score
    """
    fused_score = 0.0
    for scanner, score in results.items():
        weight = SCANNER_WEIGHTS.get(scanner, 0.1)
        fused_score += weight * score
    return fused_score
```

---

## 5. DETECTOR DATABASE SCHEMA

```sql
-- Valid Detectors (Immunocompetent T-Cells)
CREATE TABLE detectors (
    id INTEGER PRIMARY KEY,
    class_label TEXT NOT NULL,          -- 'SAFE', 'MALICIOUS', etc.
    center_vector BLOB NOT NULL,        -- Feature vector (d^j)
    radius REAL NOT NULL,               -- r^j = R^q - R^self
    r_self REAL NOT NULL,               -- Threshold used
    created_at TIMESTAMP DEFAULT NOW,
    last_match_at TIMESTAMP,
    match_count INTEGER DEFAULT 0
);

-- Self Samples (Training Data)
CREATE TABLE self_samples (
    id INTEGER PRIMARY KEY,
    class_label TEXT NOT NULL,
    feature_vector BLOB NOT NULL,
    source_file TEXT,                   -- Where this pattern came from
    created_at TIMESTAMP DEFAULT NOW
);

-- Scan Results
CREATE TABLE scan_results (
    id INTEGER PRIMARY KEY,
    file_path TEXT NOT NULL,
    feature_vector BLOB,
    classification TEXT,                -- 'self' or 'non-self'
    confidence REAL,
    matching_detector_id INTEGER,
    scanned_at TIMESTAMP DEFAULT NOW,
    FOREIGN KEY (matching_detector_id) REFERENCES detectors(id)
);
```

---

## 6. API DESIGN

### Training Endpoint

```python
@app.route('/api/train', methods=['POST'])
def train_detectors():
    """
    Train new detectors from safe code samples.
    
    Request Body:
    {
        "class_label": "SAFE",
        "samples": [
            {"code": "...", "source": "verified_lib.py"},
            ...
        ],
        "config": {
            "num_detectors": 20,
            "r_self": 0.85
        }
    }
    """
    # Extract features from code samples
    # Apply negative selection
    # Store valid detectors
    pass
```

### Scanning Endpoint

```python
@app.route('/api/scan', methods=['POST'])
def scan_code():
    """
    Scan code against trained detectors.
    
    Request Body:
    {
        "code": "import os\nos.system('rm -rf /')",
        "file_path": "malicious.py"
    }
    
    Response:
    {
        "classification": "non-self",
        "confidence": 0.92,
        "threat_level": "HIGH",
        "matching_detectors": [...],
        "recommendations": [...]
    }
    """
    pass
```

---

## 7. DASHBOARD COMPONENTS

### Real-time Metrics (from paper's evaluation)

```javascript
const DASHBOARD_METRICS = {
    // Overall Performance
    overall_accuracy: 0.94,
    cohens_kappa: 0.919,
    mcc: 0.920,
    
    // Per-Scanner Metrics
    scanners: {
        'bandit': { precision: 0.92, recall: 0.88, f1: 0.90 },
        'semgrep': { precision: 0.95, recall: 0.91, f1: 0.93 },
        'codeql': { precision: 0.89, recall: 0.94, f1: 0.91 }
    },
    
    // Detector Stats
    detectors: {
        total: 75,
        by_class: {
            'safe_patterns': 15,
            'injection_patterns': 20,
            'xss_patterns': 15,
            'rce_patterns': 25
        }
    },
    
    // Detection Performance
    detection: {
        true_positives: 1069,
        false_positives: 25,
        true_negatives: 965,
        false_negatives: 79
    }
};
```

### Visualization: Generalization Error Over Time

```javascript
// Chart showing generalization error convergence
// (Similar to Figure 11 in paper)
const generalizationChart = {
    title: "Model Generalization Error",
    xAxis: "Training Epochs",
    yAxis: "Error Rate",
    series: [
        { name: "Training Error", color: "blue" },
        { name: "Validation Error", color: "orange" },
        { name: "Generalization Error", color: "gray" }
    ],
    annotation: "Minimum Error: 0.00051"
};
```

---

## 8. IMPLEMENTATION CHECKLIST

### Phase 1: Core Algorithm
- [ ] Implement `NegativeSelectionClassifier` class
- [ ] Implement `Detector` data structure
- [ ] Implement Euclidean distance calculation
- [ ] Implement detector validation (R_q > R_self)
- [ ] Unit tests for core algorithm

### Phase 2: Feature Extraction
- [ ] Code complexity features
- [ ] AST-based features
- [ ] Security-specific features
- [ ] Feature normalization

### Phase 3: Training Pipeline
- [ ] Safe code sample collection
- [ ] Feature extraction pipeline
- [ ] Detector generation
- [ ] Detector validation and storage

### Phase 4: Scanning Pipeline
- [ ] Code submission endpoint
- [ ] Feature extraction from new code
- [ ] Detector matching
- [ ] Result aggregation (scanner fusion)

### Phase 5: Dashboard
- [ ] Real-time metrics display
- [ ] Detector visualization
- [ ] Generalization error tracking
- [ ] Scan history

---

## 9. KEY PARAMETERS TO TUNE

| Parameter | Paper Value | Suggested for IMMUNOS | Notes |
|-----------|-------------|----------------------|-------|
| num_detectors | 15-25 per class | 20 per threat type | More for common threats |
| r_self | 0.89-1.33 | 0.85-1.0 | Lower = more sensitive |
| feature_dim | 7690 | 50-100 | Depends on feature set |
| window_size | 1536 samples | N/A | For streaming code analysis |

---

## 10. REFERENCES

- Umair et al., "NegSl-AIS", Results in Engineering 27 (2025)
- Original Negative Selection: Forrest et al. (1994)
- MAHNOB-HCI Dataset documentation
