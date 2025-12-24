---
source: /Users/byron/projects/immunos81/core/clone.py
relative: immunos81/core/clone.py
generated_at: 2025-12-23 10:28
---

```python
"""
Clone representation for Immunos-81

A clone is a mathematical representation of a population of B cells that recognize
the same pattern. Clones are the basic unit of recognition in Immunos-81.
"""

from typing import List, Dict, Tuple
from dataclasses import dataclass, field
import math
from .bcell import BCell
from .antigen import Antigen


@dataclass
class Clone:
    """
    Represents a population of B cells (clone) for pattern recognition.

    A clone groups B cells that recognize similar patterns and calculates
    an overall "avidity" (recognition strength) for test antigens.

    Attributes:
        class_label: The class this clone represents
        bcells: List of B cells in this clone
        identifier: Optional unique identifier
        concentration: Clone size factor (affects avidity calculation)
    """
    class_label: str
    bcells: List[BCell] = field(default_factory=list)
    identifier: str = None
    concentration: float = 1.0

    def __post_init__(self):
        """Initialize clone properties"""
        if self.identifier is None:
            self.identifier = f"clone_{id(self)}"

    def add_bcell(self, bcell: BCell):
        """
        Add a B cell to this clone.

        Args:
            bcell: B cell to add

        Raises:
            ValueError: If B cell has different class label
        """
        if bcell.class_label != self.class_label:
            raise ValueError(
                f"Cannot add B cell with class {bcell.class_label} "
                f"to clone with class {self.class_label}"
            )
        self.bcells.append(bcell)

    def size(self) -> int:
        """Return the number of B cells in this clone"""
        return len(self.bcells)

    def calculate_avidity(self, test_antigen: Antigen) -> float:
        """
        Calculate the total avidity (recognition strength) for a test antigen.

        Avidity is calculated as the sum of affinities across all paratopic sites
        (variables), adjusted by clone size and concentration factors.

        Args:
            test_antigen: Antigen to calculate avidity for

        Returns:
            Total avidity score (higher = stronger recognition)

        Formula:
            avidity = (sum of affinities for all variables) * concentration_factor

        The concentration factor represents the population size effect:
            concentration_factor = log(1 + clone_size) * concentration
        """
        if not self.bcells:
            return 0.0

        # Calculate affinity for each variable across all B cells
        total_affinity = 0.0
        variable_count = {}  # Track how many B cells recognize each variable

        for bcell in self.bcells:
            for var_num in test_antigen.variable_sequence:
                affinity = bcell.calculate_affinity(test_antigen, var_num)
                total_affinity += affinity

                if var_num not in variable_count:
                    variable_count[var_num] = 0
                variable_count[var_num] += 1

        # Apply concentration factor (clone size effect)
        # Larger clones have more recognition power
        concentration_factor = math.log(1 + self.size()) * self.concentration

        # Total avidity = sum of affinities * concentration factor
        avidity = total_affinity * concentration_factor

        return avidity

    def calculate_affinity_vector(self, test_antigen: Antigen) -> Dict[int, float]:
        """
        Calculate affinity for each variable in the test antigen.

        This provides a detailed breakdown of which variables contribute
        most to the recognition.

        Args:
            test_antigen: Antigen to calculate affinities for

        Returns:
            Dictionary mapping variable numbers to average affinity values
        """
        if not self.bcells:
            return {}

        # Accumulate affinities for each variable
        affinity_sums = {}
        affinity_counts = {}

        for bcell in self.bcells:
            for var_num in test_antigen.variable_sequence:
                affinity = bcell.calculate_affinity(test_antigen, var_num)

                if var_num not in affinity_sums:
                    affinity_sums[var_num] = 0.0
                    affinity_counts[var_num] = 0

                affinity_sums[var_num] += affinity
                affinity_counts[var_num] += 1

        # Calculate average affinities
        affinity_vector = {}
        for var_num in affinity_sums:
            affinity_vector[var_num] = affinity_sums[var_num] / affinity_counts[var_num]

        return affinity_vector

    def get_representative_antigen(self) -> Antigen:
        """
        Get a representative antigen for this clone.

        Returns the antigen from the first B cell as a representative.
        """
        if not self.bcells:
            raise ValueError("Clone has no B cells")
        return self.bcells[0].antigen

    def merge(self, other: "Clone"):
        """
        Merge another clone into this one.

        Args:
            other: Clone to merge

        Raises:
            ValueError: If clones have different class labels
        """
        if other.class_label != self.class_label:
            raise ValueError(
                f"Cannot merge clone with class {other.class_label} "
                f"into clone with class {self.class_label}"
            )

        self.bcells.extend(other.bcells)

    def __repr__(self) -> str:
        return f"Clone(class={self.class_label}, size={self.size()}, id={self.identifier})"

```
