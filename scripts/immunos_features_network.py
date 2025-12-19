#!/usr/bin/env python3
"""
IMMUNOS Network Intrusion Detection Feature Extractor
======================================================

Extracts 41-dimensional feature vectors from network traffic for intrusion detection.

Features: NSL-KDD standard 41 features
- Basic features (9): duration, protocol_type, service, flag, src_bytes, dst_bytes, land, wrong_fragment, urgent
- Content features (13): hot, num_failed_logins, logged_in, num_compromised, root_shell, su_attempted, num_root, num_file_creations, num_shells, num_access_files, num_outbound_cmds, is_host_login, is_guest_login
- Traffic features (9): count, srv_count, serror_rate, srv_serror_rate, rerror_rate, srv_rerror_rate, same_srv_rate, diff_srv_rate, srv_diff_host_rate
- Host-based features (10): dst_host_count, dst_host_srv_count, dst_host_same_srv_rate, dst_host_diff_srv_rate, dst_host_same_src_port_rate, dst_host_srv_diff_host_rate, dst_host_serror_rate, dst_host_srv_serror_rate, dst_host_rerror_rate, dst_host_srv_rerror_rate

Dependencies:
    pip install pandas scikit-learn

Datasets:
- NSL-KDD (148,517 records) - Standard benchmark
- CICIDS2017 (2.8M flows) - Larger, more realistic

Based on:
- Tavallaee et al. (2009) NSL-KDD dataset
- Standard intrusion detection feature set
"""

import numpy as np
import sys
from pathlib import Path
from typing import Union, Dict, List
import warnings

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.immunos_universal import FeatureExtractor

# Try to import data processing libraries
try:
    import pandas as pd
    from sklearn.preprocessing import LabelEncoder, StandardScaler
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    warnings.warn("pandas/scikit-learn not installed. Install with: pip install pandas scikit-learn")


class NetworkFeatureExtractor(FeatureExtractor):
    """
    Feature extractor for network intrusion detection.

    Extracts 41 NSL-KDD standard features from network traffic.
    """

    # NSL-KDD feature names (41 total)
    FEATURE_NAMES = [
        # Basic features (9)
        'duration', 'protocol_type', 'service', 'flag',
        'src_bytes', 'dst_bytes', 'land', 'wrong_fragment', 'urgent',

        # Content features (13)
        'hot', 'num_failed_logins', 'logged_in', 'num_compromised',
        'root_shell', 'su_attempted', 'num_root', 'num_file_creations',
        'num_shells', 'num_access_files', 'num_outbound_cmds',
        'is_host_login', 'is_guest_login',

        # Traffic features (9)
        'count', 'srv_count', 'serror_rate', 'srv_serror_rate',
        'rerror_rate', 'srv_rerror_rate', 'same_srv_rate',
        'diff_srv_rate', 'srv_diff_host_rate',

        # Host-based features (10)
        'dst_host_count', 'dst_host_srv_count', 'dst_host_same_srv_rate',
        'dst_host_diff_srv_rate', 'dst_host_same_src_port_rate',
        'dst_host_srv_diff_host_rate', 'dst_host_serror_rate',
        'dst_host_srv_serror_rate', 'dst_host_rerror_rate',
        'dst_host_srv_rerror_rate'
    ]

    # Categorical features that need encoding
    CATEGORICAL_FEATURES = ['protocol_type', 'service', 'flag']

    def __init__(self):
        """Initialize network feature extractor."""
        self.encoders = {}
        self.scaler = None

        if HAS_SKLEARN:
            # Initialize label encoders for categorical features
            for feature in self.CATEGORICAL_FEATURES:
                self.encoders[feature] = LabelEncoder()

            # Initialize scaler
            self.scaler = StandardScaler()

            # Fit encoders with common values
            self._initialize_encoders()

    def _initialize_encoders(self):
        """Initialize encoders with common NSL-KDD values"""
        # Common protocol types
        self.encoders['protocol_type'].fit(['tcp', 'udp', 'icmp'])

        # Common services (subset of NSL-KDD services)
        services = [
            'http', 'smtp', 'finger', 'domain_u', 'auth', 'telnet',
            'ftp', 'eco_i', 'ntp_u', 'ecr_i', 'other', 'private',
            'pop_3', 'ftp_data', 'rje', 'time', 'mtp', 'link',
            'remote_job', 'gopher', 'ssh', 'name', 'whois', 'domain',
            'login', 'imap4', 'daytime', 'ctf', 'nntp', 'shell',
            'IRC', 'nnsp', 'http_443', 'exec', 'printer', 'efs',
            'courier', 'uucp', 'klogin', 'kshell', 'echo', 'discard',
            'systat', 'supdup', 'iso_tsap', 'hostnames', 'csnet_ns',
            'pop_2', 'sunrpc', 'uucp_path', 'netbios_ns', 'netbios_ssn',
            'netbios_dgm', 'sql_net', 'vmnet', 'bgp', 'Z39_50', 'ldap',
            'netstat', 'urh_i', 'X11', 'urp_i', 'pm_dump', 'tftp_u',
            'tim_i', 'red_i'
        ]
        self.encoders['service'].fit(services)

        # Common flags
        self.encoders['flag'].fit([
            'SF', 'S0', 'REJ', 'RSTR', 'RSTO', 'SH', 'S1', 'S2', 'RSTOS0', 'S3', 'OTH'
        ])

    def extract(self, sample: Union[Dict, List, np.ndarray]) -> np.ndarray:
        """
        Extract 41-dimensional feature vector from network flow.

        Args:
            sample: Network flow as dict, list, or array with 41 features

        Returns:
            41-dimensional feature vector
        """
        # Convert to dict if needed
        if isinstance(sample, (list, np.ndarray)):
            if len(sample) != 41:
                raise ValueError(f"Expected 41 features, got {len(sample)}")
            # Create dict from list/array
            flow = dict(zip(self.FEATURE_NAMES, sample))
        elif isinstance(sample, dict):
            flow = sample
        else:
            raise ValueError(f"Unsupported sample type: {type(sample)}")

        # Preprocess
        flow = self.preprocess(flow)

        # Extract features
        features = []

        for feature_name in self.FEATURE_NAMES:
            value = flow.get(feature_name, 0)

            # Encode categorical features
            if feature_name in self.CATEGORICAL_FEATURES and HAS_SKLEARN:
                try:
                    encoded = self.encoders[feature_name].transform([value])[0]
                    features.append(float(encoded))
                except:
                    # Unknown category - use default
                    features.append(0.0)
            else:
                # Numerical feature
                try:
                    features.append(float(value))
                except:
                    features.append(0.0)

        feature_vector = np.array(features)

        # Normalize (optional)
        if self.scaler is not None and HAS_SKLEARN:
            try:
                # Reshape for scaler
                feature_vector = feature_vector.reshape(1, -1)
                feature_vector = self.scaler.transform(feature_vector)
                feature_vector = feature_vector.flatten()
            except:
                # Scaler not fitted - use raw features
                pass

        # Ensure exactly 41 dimensions
        assert feature_vector.shape[0] == 41, f"Expected 41 features, got {feature_vector.shape[0]}"

        return feature_vector

    def get_feature_dim(self) -> int:
        """Return dimensionality of feature vector"""
        return 41

    def preprocess(self, flow: Dict) -> Dict:
        """
        Preprocess network flow.

        Args:
            flow: Network flow dictionary

        Returns:
            Preprocessed flow
        """
        # Ensure all features present with defaults
        processed = {}

        for feature in self.FEATURE_NAMES:
            if feature in flow:
                processed[feature] = flow[feature]
            else:
                # Default values
                if feature in self.CATEGORICAL_FEATURES:
                    processed[feature] = 'other'
                else:
                    processed[feature] = 0

        return processed


# =============================================================================
# TESTING AND UTILITY FUNCTIONS
# =============================================================================

def test_extractor():
    """Test network feature extractor with sample data"""
    print("Testing Network Feature Extractor")
    print("=" * 60)

    # Create extractor
    extractor = NetworkFeatureExtractor()

    print(f"✓ Extractor initialized")
    print(f"  Encoders: {len(extractor.encoders)}")
    print()

    # Test 1: Sample network flow (dict)
    print("Test 1: Network flow as dictionary")
    flow1 = {
        'duration': 0,
        'protocol_type': 'tcp',
        'service': 'http',
        'flag': 'SF',
        'src_bytes': 181,
        'dst_bytes': 5450,
        'land': 0,
        'wrong_fragment': 0,
        'urgent': 0,
        'hot': 0,
        'num_failed_logins': 0,
        'logged_in': 1,
        'num_compromised': 0,
        'root_shell': 0,
        'su_attempted': 0,
        'num_root': 0,
        'num_file_creations': 0,
        'num_shells': 0,
        'num_access_files': 0,
        'num_outbound_cmds': 0,
        'is_host_login': 0,
        'is_guest_login': 0,
        'count': 8,
        'srv_count': 8,
        'serror_rate': 0.0,
        'srv_serror_rate': 0.0,
        'rerror_rate': 0.0,
        'srv_rerror_rate': 0.0,
        'same_srv_rate': 1.0,
        'diff_srv_rate': 0.0,
        'srv_diff_host_rate': 0.0,
        'dst_host_count': 9,
        'dst_host_srv_count': 9,
        'dst_host_same_srv_rate': 1.0,
        'dst_host_diff_srv_rate': 0.0,
        'dst_host_same_src_port_rate': 0.11,
        'dst_host_srv_diff_host_rate': 0.0,
        'dst_host_serror_rate': 0.0,
        'dst_host_srv_serror_rate': 0.0,
        'dst_host_rerror_rate': 0.0,
        'dst_host_srv_rerror_rate': 0.0
    }

    features1 = extractor.extract(flow1)

    print(f"✓ Feature extraction successful")
    print(f"  Feature dimension: {features1.shape[0]}")
    print(f"  Feature range: [{features1.min():.3f}, {features1.max():.3f}]")
    print(f"  Protocol encoded: {features1[1]:.1f}")
    print(f"  Service encoded: {features1[2]:.1f}")
    print()

    # Test 2: Sample network flow (list)
    print("Test 2: Network flow as list")
    flow2_list = [0, 'udp', 'domain_u', 'SF', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1, 1, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    features2 = extractor.extract(flow2_list)

    print(f"✓ List extraction successful")
    print(f"  Feature dimension: {features2.shape[0]}")
    print()

    # Test 3: Minimal flow (missing features)
    print("Test 3: Minimal flow (missing features)")
    flow3 = {
        'protocol_type': 'icmp',
        'service': 'eco_i',
        'flag': 'SF'
    }
    features3 = extractor.extract(flow3)

    print(f"✓ Minimal flow handled (defaults applied)")
    print(f"  Feature dimension: {features3.shape[0]}")
    print()

    # Verify dimension
    assert extractor.get_feature_dim() == 41
    assert features1.shape[0] == 41
    assert features2.shape[0] == 41
    assert features3.shape[0] == 41

    print("=" * 60)
    print("✅ All tests passed!")
    print()
    print("Feature breakdown (41 dimensions):")
    print("  - Basic: 9 (duration, protocol, service, etc.)")
    print("  - Content: 13 (login attempts, compromises, etc.)")
    print("  - Traffic: 9 (connection counts, error rates, etc.)")
    print("  - Host-based: 10 (destination host statistics)")
    print("  Total: 41 dimensions (NSL-KDD standard)")


if __name__ == "__main__":
    test_extractor()
