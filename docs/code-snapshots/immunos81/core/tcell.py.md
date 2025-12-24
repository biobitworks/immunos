---
source: /Users/byron/projects/immunos81/core/tcell.py
relative: immunos81/core/tcell.py
generated_at: 2025-12-23 10:28
---

```python
"""
T Cell representation for Immunos-81

T cells are control agents that represent a particular class and determine:
1. The sequence and types of variables within an antigen
2. Whether an antigen has been previously encountered
3. Which B-cell clones should be consulted for recognition
"""

from typing import List, Tuple, Optional
from dataclasses import dataclass, field
from .amino_acid_library import AminoAcidLibrary, DataType
from .antigen import Antigen
from .clone import Clone


@dataclass
class TCell:
    """
    Represents a T cell that controls recognition for a specific class.

    T cells learn the "primary structure" of antigens belonging to their class,
    which is the ordered sequence of variable types. They validate that test
    antigens match this structure before allowing clones to perform recognition.

    Attributes:
        class_label: The class this T cell represents
        primary_structure: Ordered list of (variable_number, data_type) tuples
        library: Reference to the Amino Acid Library
        clones: List of clones associated with this T cell
    """
    class_label: str
    primary_structure: List[Tuple[int, DataType]] = field(default_factory=list)
    library: AminoAcidLibrary = None
    clones: List[Clone] = field(default_factory=list)

    def learn_structure(self, antigen: Antigen):
        """
        Learn the primary structure from a training antigen.

        This should be called during the training phase to establish the
        expected variable sequence and types for this class.

        Args:
            antigen: Training antigen to learn from

        Raises:
            ValueError: If antigen has different class label
        """
        if antigen.class_label != self.class_label:
            raise ValueError(
                f"Cannot learn from antigen with class {antigen.class_label} "
                f"for T cell with class {self.class_label}"
            )

        # Set library reference if not already set
        if self.library is None:
            self.library = antigen.library

        # If no structure learned yet, use this antigen's structure
        if not self.primary_structure:
            self.primary_structure = antigen.get_primary_structure()
        else:
            # Verify structure matches
            antigen_structure = antigen.get_primary_structure()
            if antigen_structure != self.primary_structure:
                raise ValueError(
                    f"Antigen structure {antigen_structure} does not match "
                    f"learned structure {self.primary_structure}"
                )

    def matches_structure(self, antigen: Antigen) -> bool:
        """
        Check if an antigen matches this T cell's learned primary structure.

        Args:
            antigen: Antigen to check

        Returns:
            True if the antigen matches the primary structure, False otherwise
        """
        if not self.primary_structure:
            return False

        antigen_structure = antigen.get_primary_structure()
        return antigen_structure == self.primary_structure

    def add_clone(self, clone: Clone):
        """
        Add a clone to this T cell.

        Args:
            clone: Clone to add

        Raises:
            ValueError: If clone has different class label
        """
        if clone.class_label != self.class_label:
            raise ValueError(
                f"Cannot add clone with class {clone.class_label} "
                f"to T cell with class {self.class_label}"
            )
        self.clones.append(clone)

    def recognize(self, antigen: Antigen) -> Optional[Tuple[Clone, float]]:
        """
        Attempt to recognize an antigen using this T cell's clones.

        This method:
        1. Checks if the antigen matches the primary structure
        2. If yes, calculates avidity for all clones
        3. Returns the clone with highest avidity

        Args:
            antigen: Antigen to recognize

        Returns:
            Tuple of (winning_clone, avidity) if structure matches, None otherwise
        """
        # First, check if antigen matches primary structure
        if not self.matches_structure(antigen):
            return None

        # Calculate avidity for all clones
        if not self.clones:
            return None

        best_clone = None
        best_avidity = -1.0

        for clone in self.clones:
            avidity = clone.calculate_avidity(antigen)
            if avidity > best_avidity:
                best_avidity = avidity
                best_clone = clone

        if best_clone is None:
            return None

        return (best_clone, best_avidity)

    def get_all_avidities(self, antigen: Antigen) -> List[Tuple[Clone, float]]:
        """
        Get avidity scores for all clones.

        Args:
            antigen: Antigen to calculate avidities for

        Returns:
            List of (clone, avidity) tuples, sorted by avidity (descending)
        """
        if not self.matches_structure(antigen):
            return []

        avidities = []
        for clone in self.clones:
            avidity = clone.calculate_avidity(antigen)
            avidities.append((clone, avidity))

        # Sort by avidity (descending)
        avidities.sort(key=lambda x: x[1], reverse=True)
        return avidities

    def get_variable_sequence(self) -> List[int]:
        """Get the sequence of variable numbers for this T cell's structure"""
        return [var_num for var_num, _ in self.primary_structure]

    def get_structure_info(self) -> str:
        """Get a human-readable description of the primary structure"""
        if not self.library or not self.primary_structure:
            return "No structure learned"

        parts = []
        for var_num, data_type in self.primary_structure:
            var_name = self.library.get_name(var_num)
            parts.append(f"{var_name}({data_type.value})")

        return " -> ".join(parts)

    def __repr__(self) -> str:
        return f"TCell(class={self.class_label}, clones={len(self.clones)}, structure_len={len(self.primary_structure)})"

```
