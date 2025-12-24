---
source: /Users/byron/projects/immunos81/core/antigen.py
relative: immunos81/core/antigen.py
generated_at: 2025-12-23 10:28
---

```python
"""
Antigen representation for Immunos-81

An antigen is a data grouping from the real world consisting of multiple variables
of any type, presented in a specific order.
"""

from typing import List, Any, Optional, Dict
from dataclasses import dataclass
from .amino_acid_library import AminoAcidLibrary, DataType


@dataclass
class Antigen:
    """
    Represents a data record (antigen) to be classified.

    An antigen consists of an ordered sequence of variable values and
    an optional class label (for training data).

    Attributes:
        class_label: The class this antigen belongs to (e.g., "CAD+", "CAD-")
        values: Ordered list of variable values
        variable_sequence: Ordered list of amino acid numbers (variable IDs)
        library: Reference to the Amino Acid Library
        identifier: Optional unique identifier for this antigen
    """
    class_label: Optional[str]
    values: List[Any]
    variable_sequence: List[int]
    library: AminoAcidLibrary
    identifier: Optional[str] = None

    def __post_init__(self):
        """Validate antigen structure"""
        if len(self.values) != len(self.variable_sequence):
            raise ValueError(
                f"Number of values ({len(self.values)}) must match "
                f"number of variables ({len(self.variable_sequence)})"
            )

    @classmethod
    def from_dict(cls, data: Dict[str, Any], library: AminoAcidLibrary,
                  class_label: Optional[str] = None,
                  identifier: Optional[str] = None) -> "Antigen":
        """
        Create an antigen from a dictionary of variable names and values.

        Args:
            data: Dictionary mapping variable names to values
            library: Amino Acid Library
            class_label: Optional class label
            identifier: Optional identifier

        Returns:
            Antigen instance

        Example:
            >>> library = AminoAcidLibrary()
            >>> library.add_variable("Age", "numeric")
            >>> library.add_variable("Sex", "nominal")
            >>> antigen = Antigen.from_dict(
            ...     {"Age": 45, "Sex": "M"},
            ...     library,
            ...     class_label="CAD+"
            ... )
        """
        # Get variable numbers in order
        variable_names = list(data.keys())
        variable_sequence = [library.get_number(name) for name in variable_names]

        # Check for unknown variables
        if None in variable_sequence:
            unknown = [name for name, num in zip(variable_names, variable_sequence) if num is None]
            raise ValueError(f"Unknown variables: {unknown}")

        values = list(data.values())

        return cls(
            class_label=class_label,
            values=values,
            variable_sequence=variable_sequence,
            library=library,
            identifier=identifier
        )

    @classmethod
    def from_sequence(cls, values: List[Any], variable_sequence: List[int],
                     library: AminoAcidLibrary,
                     class_label: Optional[str] = None,
                     identifier: Optional[str] = None) -> "Antigen":
        """
        Create an antigen from ordered values and variable sequence.

        Args:
            values: List of variable values
            variable_sequence: List of amino acid numbers
            library: Amino Acid Library
            class_label: Optional class label
            identifier: Optional identifier

        Returns:
            Antigen instance
        """
        return cls(
            class_label=class_label,
            values=values,
            variable_sequence=variable_sequence,
            library=library,
            identifier=identifier
        )

    def get_value(self, variable_number: int) -> Optional[Any]:
        """Get the value for a specific variable by its amino acid number"""
        try:
            idx = self.variable_sequence.index(variable_number)
            return self.values[idx]
        except ValueError:
            return None

    def get_value_by_name(self, variable_name: str) -> Optional[Any]:
        """Get the value for a specific variable by its name"""
        variable_number = self.library.get_number(variable_name)
        if variable_number is None:
            return None
        return self.get_value(variable_number)

    def has_variable(self, variable_number: int) -> bool:
        """Check if this antigen contains a specific variable"""
        return variable_number in self.variable_sequence

    def get_primary_structure(self) -> List[tuple]:
        """
        Get the primary structure (variable types and sequence) of this antigen.

        Returns:
            List of (variable_number, data_type) tuples
        """
        structure = []
        for var_num in self.variable_sequence:
            data_type = self.library.get_type(var_num)
            structure.append((var_num, data_type))
        return structure

    def to_dict(self) -> Dict[str, Any]:
        """Convert antigen to dictionary of variable names and values"""
        result = {}
        for var_num, value in zip(self.variable_sequence, self.values):
            var_name = self.library.get_name(var_num)
            if var_name:
                result[var_name] = value
        return result

    def __len__(self) -> int:
        """Return the number of variables in this antigen"""
        return len(self.values)

    def __repr__(self) -> str:
        data_dict = self.to_dict()
        return f"Antigen(class={self.class_label}, data={data_dict}, id={self.identifier})"

```
