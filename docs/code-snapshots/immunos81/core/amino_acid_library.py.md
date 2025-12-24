---
source: /Users/byron/projects/immunos81/core/amino_acid_library.py
relative: immunos81/core/amino_acid_library.py
generated_at: 2025-12-23 10:28
---

```python
"""
Amino Acid Library for Immunos-81

The Amino Acid Library (AALib) is a registry of all variables encountered by the system.
Each variable is assigned a unique identifier, name, and data type to ensure consistent
definitions throughout the learning and recognition process.
"""

from enum import Enum
from typing import Dict, List, Optional, Union
from dataclasses import dataclass


class DataType(Enum):
    """Data types supported by Immunos-81"""
    NUMERIC = "numeric"      # Continuous numerical values
    NOMINAL = "nominal"      # Categorical without order (e.g., color, sex)
    ORDINAL = "ordinal"      # Categorical with order (e.g., severity levels)


@dataclass
class AminoAcid:
    """
    Represents a single variable (amino acid) in the library.

    Attributes:
        number: Unique identifier for this amino acid
        name: Human-readable variable name
        data_type: Type of data (numeric, nominal, or ordinal)
    """
    number: int
    name: str
    data_type: DataType

    def __repr__(self) -> str:
        return f"AminoAcid({self.number}, '{self.name}', {self.data_type.value})"


class AminoAcidLibrary:
    """
    Registry of all variables in the Immunos-81 system.

    The library maintains a consistent mapping between variable names and their
    properties, ensuring that all antigens use the same variable definitions.
    """

    def __init__(self):
        self._amino_acids: Dict[int, AminoAcid] = {}
        self._name_to_number: Dict[str, int] = {}
        self._next_number = 1

    def add_variable(self, name: str, data_type: Union[DataType, str]) -> int:
        """
        Add a new variable to the library or return existing ID.

        Args:
            name: Variable name
            data_type: Type of data (DataType enum or string)

        Returns:
            The amino acid number (ID) for this variable

        Example:
            >>> lib = AminoAcidLibrary()
            >>> age_id = lib.add_variable("Age", "numeric")
            >>> sex_id = lib.add_variable("Sex", "nominal")
        """
        # Check if variable already exists
        if name in self._name_to_number:
            return self._name_to_number[name]

        # Convert string to DataType enum if needed
        if isinstance(data_type, str):
            data_type = DataType(data_type)

        # Create new amino acid
        number = self._next_number
        amino_acid = AminoAcid(number, name, data_type)

        self._amino_acids[number] = amino_acid
        self._name_to_number[name] = number
        self._next_number += 1

        return number

    def get_by_number(self, number: int) -> Optional[AminoAcid]:
        """Get amino acid by its number (ID)"""
        return self._amino_acids.get(number)

    def get_by_name(self, name: str) -> Optional[AminoAcid]:
        """Get amino acid by its name"""
        number = self._name_to_number.get(name)
        if number is not None:
            return self._amino_acids[number]
        return None

    def get_number(self, name: str) -> Optional[int]:
        """Get the number (ID) for a variable name"""
        return self._name_to_number.get(name)

    def get_name(self, number: int) -> Optional[str]:
        """Get the name for a variable number (ID)"""
        amino_acid = self._amino_acids.get(number)
        return amino_acid.name if amino_acid else None

    def get_type(self, number: int) -> Optional[DataType]:
        """Get the data type for a variable number (ID)"""
        amino_acid = self._amino_acids.get(number)
        return amino_acid.data_type if amino_acid else None

    def get_all_variables(self) -> List[AminoAcid]:
        """Get all variables in the library, ordered by number"""
        return [self._amino_acids[num] for num in sorted(self._amino_acids.keys())]

    def __len__(self) -> int:
        """Return the number of variables in the library"""
        return len(self._amino_acids)

    def __repr__(self) -> str:
        lines = ["AminoAcidLibrary:"]
        for aa in self.get_all_variables():
            lines.append(f"  {aa}")
        return "\n".join(lines)

```
