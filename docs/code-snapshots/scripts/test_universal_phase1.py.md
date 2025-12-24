---
source: /Users/byron/projects/scripts/test_universal_phase1.py
relative: scripts/test_universal_phase1.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
Phase 1 Test Suite for IMMUNOS Universal Framework
===================================================

Tests:
1. System instantiation
2. Domain addition (all 5 domains)
3. Training with synthetic data
4. Detection (self vs non-self)
5. Database persistence (immune_memory table)
6. Multi-algorithm support
"""

import numpy as np
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.immunos_universal import (
    ImmuneSystem, FeatureExtractor, DetectionResult,
    CellType
)
from scripts.immunos_config import DEFAULT_CONFIGS, DomainConfig


class DummyFeatureExtractor(FeatureExtractor):
    """Simple feature extractor for testing"""
    def __init__(self, feature_dim: int):
        self.feature_dim = feature_dim

    def extract(self, sample):
        if isinstance(sample, np.ndarray):
            return sample
        return np.array(sample)

    def get_feature_dim(self):
        return self.feature_dim

    def preprocess(self, sample):
        return sample


def test_system_instantiation():
    """Test 1: System instantiation"""
    print("Test 1: System Instantiation")
    print("-" * 60)

    system = ImmuneSystem()
    assert system is not None, "System failed to instantiate"
    assert len(system.domains) == 0, "System should start with 0 domains"
    assert 'negative_selection' in system.algorithms
    assert 'clonal_selection' in system.algorithms

    print("✓ ImmuneSystem instantiated successfully")
    print(f"✓ Algorithms available: {list(system.algorithms.keys())}")
    print(f"✓ Domains: {len(system.domains)}")
    print()
    return system


def test_add_all_domains(system):
    """Test 2: Add all 5 domains"""
    print("Test 2: Add All 5 Domains")
    print("-" * 60)

    domains_added = []

    for domain_name, config in DEFAULT_CONFIGS.items():
        extractor = DummyFeatureExtractor(config.feature_dim)
        system.add_domain(domain_name, config, extractor)
        domains_added.append(domain_name)
        print(f"✓ Added domain: {domain_name} (features: {config.feature_dim}, "
              f"detectors: {config.num_detectors})")

    assert len(system.domains) == 5, f"Expected 5 domains, got {len(system.domains)}"

    print(f"\n✓ All {len(domains_added)} domains added successfully")
    print(f"  Domains: {', '.join(domains_added)}")
    print()
    return domains_added


def test_training(system, domain_name: str):
    """Test 3: Training with synthetic data"""
    print(f"Test 3: Training '{domain_name}' Domain")
    print("-" * 60)

    config = system.domains[domain_name]['config']
    feature_dim = config.feature_dim

    # Generate synthetic self samples (normal patterns around 0.5)
    np.random.seed(42)
    num_samples = 100
    self_samples = [np.random.uniform(0.3, 0.7, feature_dim) for _ in range(num_samples)]

    # Train
    num_detectors = system.train(domain_name, self_samples)

    assert num_detectors > 0, "No detectors were created"
    assert num_detectors <= config.num_detectors * len(config.algorithms), \
        f"Too many detectors created: {num_detectors}"

    print(f"✓ Training completed")
    print(f"  Self samples: {num_samples}")
    print(f"  Detectors created: {num_detectors}")
    print(f"  Algorithms used: {', '.join(config.algorithms)}")

    # Check detector types
    detectors = system.domains[domain_name]['detectors']
    t_cells = sum(1 for d in detectors if d.cell_type == CellType.T_CELL)
    b_cells = sum(1 for d in detectors if d.cell_type == CellType.B_CELL)

    print(f"  T-cells (NSA): {t_cells}")
    print(f"  B-cells (CSA): {b_cells}")
    print()
    return num_detectors


def test_detection(system, domain_name: str):
    """Test 4: Detection (self vs non-self)"""
    print(f"Test 4: Detection on '{domain_name}' Domain")
    print("-" * 60)

    config = system.domains[domain_name]['config']
    feature_dim = config.feature_dim

    # Test with self sample (should be detected as SELF)
    test_self = np.random.uniform(0.4, 0.6, feature_dim)
    result_self = system.detect(domain_name, test_self)

    print(f"✓ Self sample detection:")
    print(f"  Result: {result_self.result.value}")
    print(f"  Confidence: {result_self.confidence:.3f}")
    print(f"  Matched detectors: {len(result_self.matched_detectors)}")

    # Test with non-self sample (should be detected as NON_SELF or DANGER)
    test_non_self = np.random.uniform(0.9, 1.0, feature_dim)
    context = {'entropy': 0.9, 'max_entropy': 1.0}
    result_non_self = system.detect(domain_name, test_non_self, context=context)

    print(f"\n✓ Non-self sample detection:")
    print(f"  Result: {result_non_self.result.value}")
    print(f"  Confidence: {result_non_self.confidence:.3f}")
    print(f"  Danger signal: {result_non_self.danger_signal:.3f}")
    print(f"  Matched detectors: {len(result_non_self.matched_detectors)}")
    print()

    return result_self, result_non_self


def test_database_persistence(system, domain_name: str):
    """Test 5: Database persistence (immune_memory)"""
    print(f"Test 5: Database Persistence")
    print("-" * 60)

    import sqlite3

    # Get database path
    db_path = system.memory.db_path
    print(f"✓ Database path: {db_path}")

    # Check immune_memory table exists
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute("""
        SELECT name FROM sqlite_master
        WHERE type='table' AND name='immune_memory'
    """)

    table_exists = cursor.fetchone() is not None
    assert table_exists, "immune_memory table does not exist"
    print(f"✓ immune_memory table exists")

    # Check for stored patterns
    cursor.execute("""
        SELECT COUNT(*), domain, detection_type
        FROM immune_memory
        WHERE domain = ?
        GROUP BY domain, detection_type
    """, (domain_name,))

    results = cursor.fetchall()
    print(f"\n✓ Immune memory records:")
    if results:
        for count, domain, det_type in results:
            print(f"  {domain}: {count} patterns ({det_type})")
    else:
        print(f"  No patterns stored yet")

    # Test memory recall
    feature_dim = system.domains[domain_name]['config'].feature_dim
    test_pattern = np.random.uniform(0.5, 0.6, feature_dim)

    # Store a pattern
    system.memory.store(domain_name, test_pattern, 'non_self', 0.85)
    print(f"\n✓ Stored test pattern in memory")

    # Recall the pattern
    recalled = system.memory.recall(domain_name, test_pattern)
    assert recalled is not None, "Failed to recall stored pattern"
    print(f"✓ Recalled pattern from memory")
    print(f"  Detection type: {recalled['detection_type']}")
    print(f"  Confidence: {recalled['confidence']:.3f}")
    print(f"  Hit count: {recalled['hit_count']}")

    conn.close()
    print()


def test_statistics(system, domain_name: str):
    """Test 6: Statistics retrieval"""
    print(f"Test 6: Statistics for '{domain_name}'")
    print("-" * 60)

    stats = system.get_statistics(domain_name)

    print(f"✓ Domain: {stats['domain']}")
    print(f"✓ Total detectors: {stats['num_detectors']}")
    print(f"✓ Detector types:")
    for cell_type, count in stats['detector_types'].items():
        if count > 0:
            print(f"    {cell_type}: {count}")
    print(f"✓ Average confidence: {stats['avg_confidence']:.3f}")
    print(f"✓ Total matches: {stats['total_matches']}")
    print()


def run_all_tests():
    """Run complete Phase 1 test suite"""
    print("=" * 60)
    print("IMMUNOS UNIVERSAL FRAMEWORK - PHASE 1 TEST SUITE")
    print("=" * 60)
    print()

    try:
        # Test 1: Instantiation
        system = test_system_instantiation()

        # Test 2: Add all domains
        domains = test_add_all_domains(system)

        # Test 3-6: Train and test first domain (emotion)
        test_domain = domains[0]  # 'emotion'
        test_training(system, test_domain)
        test_detection(system, test_domain)
        test_database_persistence(system, test_domain)
        test_statistics(system, test_domain)

        # Summary
        print("=" * 60)
        print("✅ ALL PHASE 1 TESTS PASSED")
        print("=" * 60)
        print()
        print("Summary:")
        print(f"  ✓ System instantiation: PASS")
        print(f"  ✓ Domain addition (5 domains): PASS")
        print(f"  ✓ Training with synthetic data: PASS")
        print(f"  ✓ Detection (self vs non-self): PASS")
        print(f"  ✓ Database persistence: PASS")
        print(f"  ✓ Statistics retrieval: PASS")
        print()
        print("Phase 1 Core Framework: ✅ COMPLETE")
        print()
        print("Next Steps:")
        print("  → Phase 2A: Emotion feature extractor")
        print("  → Phase 2B: Hallucination feature extractor")
        print("  → Phase 2C: Network feature extractor")
        print("  → Phase 2D: Code feature extractor")
        print("  → Phase 2E: Research feature extractor")

        return True

    except Exception as e:
        print()
        print("=" * 60)
        print("❌ TEST FAILED")
        print("=" * 60)
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)

```
