#!/usr/bin/env python3
"""
IMMUNOS Domain Configuration
=============================

Default configurations for all 5 domains in the universal framework.

Domains:
- Emotion Recognition
- LLM Hallucination Detection
- Network Intrusion Detection
- Code Vulnerability Detection
- Research Paper Verification
"""

from dataclasses import dataclass
from typing import List, Dict


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
# PARAMETER REFERENCE FROM PAPER
# =============================================================================

PAPER_OPTIMAL_PARAMS = {
    'emotion': {
        'paper': 'Umair et al. (2025) NegSl-AIS',
        'dataset': 'FER2013 / CK+',
        'num_detectors_range': (15, 25),
        'r_self_range': (0.85, 0.95),
        'expected_accuracy': (0.73, 0.87),
        'notes': 'Facial emotion classification from images'
    },

    'hallucination': {
        'paper': 'Li et al. (2023) HaluEval',
        'dataset': 'HaluEval / TruthfulQA',
        'num_detectors_range': (25, 35),
        'r_self_range': (0.80, 0.90),
        'expected_accuracy': (0.80, 0.90),
        'notes': 'LLM output verification and hallucination detection'
    },

    'network': {
        'paper': 'Tavallaee et al. (2009) NSL-KDD',
        'dataset': 'NSL-KDD / CICIDS2017',
        'num_detectors_range': (15, 25),
        'r_self_range': (0.90, 1.00),
        'expected_accuracy': (0.93, 0.99),
        'notes': 'Network intrusion detection, high accuracy required'
    },

    'code': {
        'paper': 'Zhou et al. (2019) Devign',
        'dataset': 'Devign / DiverseVul',
        'num_detectors_range': (30, 40),
        'r_self_range': (0.75, 0.85),
        'expected_accuracy': (0.70, 0.89),
        'notes': 'Code vulnerability detection in C/C++ code'
    },

    'research': {
        'paper': 'Wadden et al. (2020) SciFact',
        'dataset': 'SciFact / PubMedQA',
        'num_detectors_range': (15, 25),
        'r_self_range': (0.85, 0.95),
        'expected_accuracy': (0.80, 0.88),
        'notes': 'Scientific claim verification'
    }
}


# =============================================================================
# DATASET CONFIGURATIONS
# =============================================================================

DATASET_CONFIGS = {
    # Emotion Recognition
    'ckplus': {
        'domain': 'emotion',
        'size': 593,
        'format': 'images',
        'download_url': 'http://www.jeffcohn.net/Resources/',
        'license': 'Academic (manual download required)',
        'm1_training_time': '5 minutes',
        'storage': '~50MB'
    },

    'fer2013': {
        'domain': 'emotion',
        'size': 35887,
        'format': 'csv',
        'download_url': 'https://www.kaggle.com/datasets/msambare/fer2013',
        'license': 'Kaggle account required',
        'm1_training_time': '30-60 minutes',
        'storage': '~100MB'
    },

    # Hallucination Detection
    'truthfulqa': {
        'domain': 'hallucination',
        'size': 817,
        'format': 'json',
        'download_url': 'https://github.com/sylinrl/TruthfulQA',
        'license': 'MIT',
        'm1_training_time': '10 minutes',
        'storage': '~5MB'
    },

    'halueval': {
        'domain': 'hallucination',
        'size': 35000,
        'format': 'json',
        'download_url': 'https://github.com/RUCAIBox/HaluEval',
        'license': 'MIT',
        'm1_training_time': '1-2 hours',
        'storage': '~500MB'
    },

    # Network Intrusion
    'nsl-kdd': {
        'domain': 'network',
        'size': 148517,
        'format': 'csv',
        'download_url': 'https://www.unb.ca/cic/datasets/nsl.html',
        'license': 'Academic',
        'm1_training_time': '15-30 minutes',
        'storage': '~25MB'
    },

    'cicids2017': {
        'domain': 'network',
        'size': 2800000,
        'format': 'csv',
        'download_url': 'https://www.unb.ca/cic/datasets/ids-2017.html',
        'license': 'Academic',
        'm1_training_time': '2-4 hours',
        'storage': '~8GB'
    },

    # Code Vulnerability
    'devign': {
        'domain': 'code',
        'size': 27318,
        'format': 'json',
        'download_url': 'https://sites.google.com/view/devign',
        'license': 'Academic',
        'm1_training_time': '30-60 minutes',
        'storage': '~500MB'
    },

    'diversevul': {
        'domain': 'code',
        'size': 18945,
        'format': 'json',
        'download_url': 'https://github.com/wagner-group/diversevul',
        'license': 'MIT',
        'm1_training_time': '1-2 hours',
        'storage': '~1GB'
    },

    # Research Verification
    'scifact': {
        'domain': 'research',
        'size': 1400,
        'format': 'json',
        'download_url': 'https://github.com/allenai/scifact',
        'license': 'Apache 2.0',
        'm1_training_time': '20-30 minutes',
        'storage': '~50MB'
    },

    'pubmedqa': {
        'domain': 'research',
        'size': 211269,
        'format': 'json',
        'download_url': 'https://pubmedqa.github.io/',
        'license': 'MIT',
        'm1_training_time': '1-2 hours',
        'storage': '~500MB'
    }
}


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def get_config(domain: str) -> DomainConfig:
    """Get configuration for a domain"""
    if domain not in DEFAULT_CONFIGS:
        raise ValueError(f"Unknown domain: {domain}. Available: {list(DEFAULT_CONFIGS.keys())}")
    return DEFAULT_CONFIGS[domain]


def get_dataset_config(dataset_name: str) -> dict:
    """Get dataset configuration"""
    if dataset_name not in DATASET_CONFIGS:
        raise ValueError(f"Unknown dataset: {dataset_name}. Available: {list(DATASET_CONFIGS.keys())}")
    return DATASET_CONFIGS[dataset_name]


def list_domains() -> List[str]:
    """List all available domains"""
    return list(DEFAULT_CONFIGS.keys())


def list_datasets(domain: str = None) -> List[str]:
    """List datasets for a domain (or all if domain=None)"""
    if domain is None:
        return list(DATASET_CONFIGS.keys())

    return [name for name, config in DATASET_CONFIGS.items()
            if config['domain'] == domain]


if __name__ == "__main__":
    print("IMMUNOS Domain Configurations")
    print("=" * 60)

    for domain_name in list_domains():
        config = get_config(domain_name)
        print(f"\n{domain_name.upper()}")
        print(f"  Detectors: {config.num_detectors}")
        print(f"  R_self: {config.r_self}")
        print(f"  Feature dim: {config.feature_dim}")
        print(f"  Algorithms: {', '.join(config.algorithms)}")
        print(f"  Datasets: {', '.join(list_datasets(domain_name))}")
