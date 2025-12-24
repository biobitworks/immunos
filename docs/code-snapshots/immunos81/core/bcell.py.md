---
source: /Users/byron/projects/immunos81/core/bcell.py
relative: immunos81/core/bcell.py
generated_at: 2025-12-23 10:28
---

```python
"""
B Cell representation for Immunos-81

B cells are entities that represent instances of a particular class during learning.
After learning, B cells are grouped into clones (populations) for pattern matching.
"""

from typing import Dict, Any
from dataclasses import dataclass
from .antigen import Antigen


@dataclass
class BCell:
    """
    Represents a B cell that recognizes a specific antigen instance.

    B cells are created during the learning phase for each training instance.
    Each B cell stores affinity values that represent how well it recognizes
    specific variable values.

    Attributes:
        class_label: The class this B cell represents
        antigen: The training antigen this B cell was created from
        affinities: Dictionary mapping variable numbers to affinity values
        identifier: Optional unique identifier
    """
    class_label: str
    antigen: Antigen
    affinities: Dict[int, float]
    identifier: str = None

    @classmethod
    def from_antigen(cls, antigen: Antigen, identifier: str = None) -> "BCell":
        """
        Create a B cell from a training antigen.

        Args:
            antigen: Training antigen
            identifier: Optional identifier

        Returns:
            BCell instance

        Note:
            Affinity calculation depends on the data type:
            - Numeric: Based on value similarity/distance
            - Nominal: 1.0 for exact match, 0.0 otherwise
            - Ordinal: Based on ordinal distance
        """
        if antigen.class_label is None:
            raise ValueError("Cannot create B cell from antigen without class label")

        # Initialize affinities for each variable
        # For now, we use a simple scheme:
        # - Affinity = 1.0 for all variables initially
        # - This will be refined during clone creation
        affinities = {var_num: 1.0 for var_num in antigen.variable_sequence}

        return cls(
            class_label=antigen.class_label,
            antigen=antigen,
            affinities=affinities,
            identifier=identifier or f"bcell_{id(antigen)}"
        )

    def calculate_affinity(self, test_antigen: Antigen, variable_number: int) -> float:
        """
        Calculate affinity between this B cell and a test antigen for a specific variable.

        Args:
            test_antigen: Antigen to test against
            variable_number: Variable number to calculate affinity for

        Returns:
            Affinity value (0.0 to 1.0)

        The affinity calculation depends on the data type:
        - Numeric: 1.0 - normalized distance
        - Nominal: 1.0 if match, 0.0 otherwise
        - Ordinal: Based on ordinal distance
        """
        # Get values from both antigens
        train_value = self.antigen.get_value(variable_number)
        test_value = test_antigen.get_value(variable_number)

        # If either doesn't have the variable, affinity is 0
        if train_value is None or test_value is None:
            return 0.0

        # Get variable type
        data_type = self.antigen.library.get_type(variable_number)

        if data_type.value == "numeric":
            # For numeric: use inverse distance (closer = higher affinity)
            # This is a simple scheme; could be improved with normalization
            try:
                train_val = float(train_value)
                test_val = float(test_value)
                # Simple affinity: 1.0 for exact match, decreases with distance
                # Using exponential decay: exp(-|diff|)
                import math
                distance = abs(train_val - test_val)
                affinity = math.exp(-distance / (1.0 + abs(train_val)))  # Normalize by scale
                return affinity
            except (ValueError, TypeError):
                return 0.0

        elif data_type.value == "nominal":
            # For nominal: exact match = 1.0, no match = 0.0
            return 1.0 if train_value == test_value else 0.0

        elif data_type.value == "ordinal":
            # For ordinal: based on distance in ordering
            # Assuming values are comparable
            try:
                distance = abs(train_value - test_value)
                max_distance = max(abs(train_value), abs(test_value), 1)
                affinity = 1.0 - (distance / max_distance)
                return max(0.0, affinity)
            except (TypeError, ValueError):
                # Fall back to nominal matching
                return 1.0 if train_value == test_value else 0.0

        return 0.0

    def __repr__(self) -> str:
        return f"BCell(class={self.class_label}, id={self.identifier})"

```
