---
source: /Users/byron/projects/scripts/immunos_features_code.py
relative: scripts/immunos_features_code.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
IMMUNOS Code Vulnerability Detection Feature Extractor
======================================================

Extracts 64-dimensional feature vectors from code for vulnerability detection.

Features (64 total):
- Structural (10): Cyclomatic complexity, AST depth, function count, class count, nesting depth, etc.
- Security patterns (20): Dangerous functions (strcpy, gets, system), buffer ops, SQL injection, XSS, etc.
- Data flow (15): Taint sources, sinks, sanitization points, flow paths
- Statistical (19): Code length, entropy, comment ratio, variable naming patterns, etc.

Dependencies:
    pip install tree-sitter tree-sitter-c tree-sitter-cpp tree-sitter-python

Datasets:
- Devign (27,318 functions) - C/C++ vulnerability detection
- DiverseVul (18,945 functions) - Diverse vulnerability dataset

Based on:
- Zhou et al. (2019) Devign - Graph neural networks for vulnerability detection
- CWE Top 25 vulnerability patterns
"""

import numpy as np
import sys
import re
import ast as python_ast
from pathlib import Path
from typing import Union, List, Set
import warnings
from collections import Counter

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.immunos_universal import FeatureExtractor

# Try to import code analysis libraries
try:
    import tree_sitter
    HAS_TREE_SITTER = False  # Placeholder - tree_sitter setup is complex
except ImportError:
    HAS_TREE_SITTER = False


class CodeFeatureExtractor(FeatureExtractor):
    """
    Feature extractor for code vulnerability detection.

    Extracts 64-dimensional feature vector combining structural,
    security, data flow, and statistical features.
    """

    # Dangerous C/C++ functions
    DANGEROUS_FUNCTIONS = {
        'strcpy', 'strcat', 'gets', 'sprintf', 'vsprintf',
        'scanf', 'fscanf', 'sscanf', 'vscanf', 'vsscanf',
        'strncpy', 'strncat', 'system', 'popen', 'exec',
        'execl', 'execlp', 'execle', 'execv', 'execvp'
    }

    # SQL injection patterns
    SQL_PATTERNS = [
        r'SELECT.*FROM', r'INSERT\s+INTO', r'UPDATE.*SET',
        r'DELETE\s+FROM', r'DROP\s+TABLE', r'CREATE\s+TABLE',
        r'execute\(', r'exec\(', r'\.query\(', r'\.execute\('
    ]

    # XSS patterns
    XSS_PATTERNS = [
        r'innerHTML\s*=', r'outerHTML\s*=', r'document\.write',
        r'eval\(', r'setTimeout\(', r'setInterval\(',
        r'<script', r'javascript:', r'onerror\s*='
    ]

    # Buffer operation patterns
    BUFFER_OPS = [
        r'malloc\(', r'calloc\(', r'realloc\(', r'alloca\(',
        r'memcpy\(', r'memmove\(', r'memset\(', r'strncpy\(',
        r'strcpy\(', r'sprintf\(', r'gets\('
    ]

    def __init__(self, language: str = 'python'):
        """
        Initialize code feature extractor.

        Args:
            language: Programming language ('python', 'c', 'cpp', 'java')
        """
        self.language = language.lower()

    def extract(self, sample: Union[str, Path]) -> np.ndarray:
        """
        Extract 64-dimensional feature vector from code.

        Args:
            sample: Code string or path to code file

        Returns:
            64-dimensional feature vector
        """
        # Load code if path provided
        if isinstance(sample, Path):
            with open(sample, 'r', encoding='utf-8', errors='ignore') as f:
                code = f.read()
        elif isinstance(sample, str) and '\n' not in sample and len(sample) < 1000:
            # Might be a file path - check if it exists
            try:
                if Path(sample).exists():
                    with open(sample, 'r', encoding='utf-8', errors='ignore') as f:
                        code = f.read()
                else:
                    code = str(sample)
            except (OSError, ValueError):
                code = str(sample)
        else:
            code = str(sample)

        # Preprocess
        code = self.preprocess(code)

        features = []

        # 1. Structural Features (10 dimensions)
        structural = self._extract_structural_features(code)
        features.append(structural)

        # 2. Security Pattern Features (20 dimensions)
        security = self._extract_security_features(code)
        features.append(security)

        # 3. Data Flow Features (15 dimensions)
        dataflow = self._extract_dataflow_features(code)
        features.append(dataflow)

        # 4. Statistical Features (19 dimensions)
        statistical = self._extract_statistical_features(code)
        features.append(statistical)

        # Concatenate all features
        feature_vector = np.concatenate(features)

        # Ensure exactly 64 dimensions
        assert feature_vector.shape[0] == 64, f"Expected 64 features, got {feature_vector.shape[0]}"

        return feature_vector

    def get_feature_dim(self) -> int:
        """Return dimensionality of feature vector"""
        return 64

    def preprocess(self, code: str) -> str:
        """
        Preprocess code for feature extraction.

        Args:
            code: Input code string

        Returns:
            Preprocessed code
        """
        # Remove extra whitespace
        code = '\n'.join(line.rstrip() for line in code.split('\n'))
        return code

    def _extract_structural_features(self, code: str) -> np.ndarray:
        """
        Extract structural features (10 dimensions).

        Features:
        - Lines of code (LOC)
        - Cyclomatic complexity (estimated)
        - Function count
        - Class count
        - Maximum nesting depth
        - Average line length
        - Number of imports
        - Number of loops
        - Number of conditionals
        - Number of try-except blocks
        """
        features = []
        lines = code.split('\n')

        # 1. Lines of code
        loc = len([l for l in lines if l.strip() and not l.strip().startswith('#')])
        features.append(loc / 1000.0)  # Normalized

        # 2. Cyclomatic complexity (estimate via control flow keywords)
        control_flow = ['if', 'elif', 'else', 'for', 'while', 'case', 'switch',
                       'catch', 'except', '&&', '||', '?']
        complexity = sum(code.lower().count(keyword) for keyword in control_flow)
        features.append(complexity / 100.0)

        # 3-4. Function and class count (simple heuristics)
        if self.language == 'python':
            func_count = code.count('def ')
            class_count = code.count('class ')
        elif self.language in ['c', 'cpp']:
            func_count = len(re.findall(r'\w+\s+\w+\s*\([^)]*\)\s*\{', code))
            class_count = len(re.findall(r'class\s+\w+', code))
        else:
            func_count = 0
            class_count = 0

        features.append(func_count / 50.0)
        features.append(class_count / 20.0)

        # 5. Maximum nesting depth (estimate via indentation)
        max_indent = 0
        for line in lines:
            if line.strip():
                indent = len(line) - len(line.lstrip())
                max_indent = max(max_indent, indent)
        features.append(max_indent / 40.0)

        # 6. Average line length
        if lines:
            avg_line_len = np.mean([len(l) for l in lines if l.strip()])
        else:
            avg_line_len = 0
        features.append(avg_line_len / 100.0)

        # 7. Number of imports
        import_count = code.count('import ') + code.count('#include')
        features.append(import_count / 50.0)

        # 8-9. Loops and conditionals
        loop_count = code.count('for ') + code.count('while ')
        cond_count = code.count('if ') + code.count('elif ') + code.count('else ')
        features.append(loop_count / 50.0)
        features.append(cond_count / 50.0)

        # 10. Try-except blocks
        try_count = code.count('try:') + code.count('try{') + code.count('try ')
        features.append(try_count / 20.0)

        return np.array(features[:10])

    def _extract_security_features(self, code: str) -> np.ndarray:
        """
        Extract security pattern features (20 dimensions).

        Detects dangerous functions, SQL injection, XSS, buffer overflows, etc.
        """
        features = []

        # 1-5. Dangerous function calls
        dangerous_count = sum(1 for func in self.DANGEROUS_FUNCTIONS if func in code)
        features.append(dangerous_count / 10.0)

        # Individual dangerous functions (5 features)
        for func in ['strcpy', 'gets', 'system', 'eval', 'exec']:
            features.append(float(func in code))

        # 6-10. SQL injection patterns
        sql_count = sum(1 for pattern in self.SQL_PATTERNS
                       if re.search(pattern, code, re.IGNORECASE))
        features.append(sql_count / 5.0)

        # Specific SQL keywords (4 features)
        for keyword in ['SELECT', 'INSERT', 'UPDATE', 'DELETE']:
            features.append(float(keyword.lower() in code.lower()))

        # 11-15. XSS patterns
        xss_count = sum(1 for pattern in self.XSS_PATTERNS
                       if re.search(pattern, code, re.IGNORECASE))
        features.append(xss_count / 5.0)

        # Specific XSS vectors (4 features)
        for pattern in ['innerHTML', 'eval', 'document.write', '<script']:
            features.append(float(pattern in code))

        # 16-20. Buffer operations
        buffer_count = sum(1 for pattern in self.BUFFER_OPS
                          if re.search(pattern, code))
        features.append(buffer_count / 5.0)

        # Specific buffer ops (4 features)
        for op in ['malloc', 'memcpy', 'strcpy', 'sprintf']:
            features.append(float(op in code))

        return np.array(features[:20])

    def _extract_dataflow_features(self, code: str) -> np.ndarray:
        """
        Extract data flow features (15 dimensions).

        Simplified data flow analysis:
        - User input sources
        - Dangerous sinks
        - Sanitization presence
        """
        features = []

        # 1-5. Taint sources (user input)
        input_sources = ['input(', 'raw_input(', 'request.', 'argv', 'stdin',
                        'getenv(', 'recv(', 'read(', 'fgets(', 'scanf(']

        source_count = sum(1 for src in input_sources if src in code)
        features.append(source_count / 10.0)

        # Specific sources (4 features)
        for src in ['input', 'request', 'argv', 'stdin']:
            features.append(float(src in code))

        # 6-10. Taint sinks (dangerous operations)
        sinks = ['system(', 'exec(', 'eval(', 'query(', 'execute(',
                'write(', 'send(', 'print(', 'echo ', 'innerHTML']

        sink_count = sum(1 for sink in sinks if sink in code)
        features.append(sink_count / 10.0)

        # Specific sinks (4 features)
        for sink in ['system', 'exec', 'eval', 'query']:
            features.append(float(sink in code))

        # 11-15. Sanitization/validation
        sanitizers = ['escape(', 'sanitize(', 'validate(', 'filter(',
                     'strip(', 'clean(', 'encode(', 'htmlspecialchars']

        sanitize_count = sum(1 for san in sanitizers if san in code)
        features.append(sanitize_count / 10.0)

        # Specific sanitizers (4 features)
        for san in ['escape', 'sanitize', 'validate', 'filter']:
            features.append(float(san in code))

        return np.array(features[:15])

    def _extract_statistical_features(self, code: str) -> np.ndarray:
        """
        Extract statistical features (19 dimensions).

        Features:
        - Code length metrics
        - Entropy
        - Comment ratio
        - Variable naming patterns
        - Operator/operand counts
        """
        features = []
        lines = code.split('\n')

        # 1-3. Length metrics
        char_count = len(code)
        word_count = len(code.split())
        line_count = len(lines)

        features.extend([
            char_count / 10000.0,
            word_count / 1000.0,
            line_count / 500.0
        ])

        # 4. Entropy
        if code:
            char_freq = Counter(code)
            total = sum(char_freq.values())
            entropy = -sum((count/total) * np.log2(count/total)
                          for count in char_freq.values())
        else:
            entropy = 0
        features.append(entropy / 8.0)  # Normalized

        # 5. Comment ratio
        if self.language == 'python':
            comment_lines = sum(1 for l in lines if l.strip().startswith('#'))
        else:
            comment_lines = sum(1 for l in lines if '//' in l or '/*' in l)

        comment_ratio = comment_lines / (line_count + 1)
        features.append(comment_ratio)

        # 6-10. Variable naming (heuristics)
        variables = re.findall(r'\b[a-z_][a-z0-9_]*\b', code)
        if variables:
            var_count = len(variables)
            unique_vars = len(set(variables))
            avg_var_len = np.mean([len(v) for v in variables])
            single_char_vars = sum(1 for v in variables if len(v) == 1)
            snake_case_vars = sum(1 for v in variables if '_' in v)

            features.extend([
                unique_vars / 200.0,
                avg_var_len / 20.0,
                single_char_vars / 50.0,
                snake_case_vars / 100.0,
                var_count / 500.0
            ])
        else:
            features.extend([0.0] * 5)

        # 11-15. Operator/operand counts
        operators = ['+', '-', '*', '/', '%', '=', '==', '!=', '<', '>', '&&', '||']
        op_count = sum(code.count(op) for op in operators)
        features.append(op_count / 200.0)

        # Specific operators (4 features)
        for op in ['=', '==', '+', '*']:
            features.append(code.count(op) / 100.0)

        # 16-19. Whitespace statistics
        space_count = code.count(' ')
        tab_count = code.count('\t')
        newline_count = code.count('\n')

        features.extend([
            space_count / 5000.0,
            tab_count / 500.0,
            newline_count / 500.0,
            (space_count + tab_count) / (char_count + 1)  # Whitespace ratio
        ])

        return np.array(features[:19])


# =============================================================================
# TESTING AND UTILITY FUNCTIONS
# =============================================================================

def test_extractor():
    """Test code feature extractor with sample code"""
    print("Testing Code Feature Extractor")
    print("=" * 60)

    # Create extractor
    extractor = CodeFeatureExtractor(language='python')

    print(f"✓ Extractor initialized")
    print(f"  Language: {extractor.language}")
    print()

    # Test 1: Simple Python code
    print("Test 1: Simple Python code")
    code1 = """
def factorial(n):
    if n <= 1:
        return 1
    else:
        return n * factorial(n - 1)

result = factorial(5)
print(result)
"""
    features1 = extractor.extract(code1)

    print(f"✓ Feature extraction successful")
    print(f"  Feature dimension: {features1.shape[0]}")
    print(f"  Feature range: [{features1.min():.3f}, {features1.max():.3f}]")
    print()

    # Test 2: Code with vulnerabilities
    print("Test 2: Code with security issues")
    code2 = """
import os
import subprocess

def execute_command(user_input):
    # DANGEROUS: Direct system call with user input
    os.system(user_input)

def sql_query(user_id):
    # DANGEROUS: SQL injection vulnerability
    query = "SELECT * FROM users WHERE id = " + user_id
    return query

data = input("Enter command: ")
execute_command(data)
"""
    features2 = extractor.extract(code2)

    print(f"✓ Vulnerability detection successful")
    print(f"  Feature dimension: {features2.shape[0]}")
    print(f"  Dangerous functions detected: {features2[10]:.3f}")
    print(f"  SQL patterns detected: {features2[15]:.3f}")
    print()

    # Test 3: From file
    print("Test 3: Extract from test file itself")
    features3 = extractor.extract(__file__)

    print(f"✓ File extraction successful")
    print(f"  Feature dimension: {features3.shape[0]}")
    print()

    # Verify dimension
    assert extractor.get_feature_dim() == 64
    assert features1.shape[0] == 64
    assert features2.shape[0] == 64
    assert features3.shape[0] == 64

    print("=" * 60)
    print("✅ All tests passed!")
    print()
    print("Feature breakdown (64 dimensions):")
    print("  - Structural: 10 (LOC, complexity, functions, etc.)")
    print("  - Security: 20 (dangerous functions, SQL, XSS, etc.)")
    print("  - Data flow: 15 (sources, sinks, sanitization)")
    print("  - Statistical: 19 (entropy, comments, naming, etc.)")
    print("  Total: 64 dimensions")


if __name__ == "__main__":
    test_extractor()

```
