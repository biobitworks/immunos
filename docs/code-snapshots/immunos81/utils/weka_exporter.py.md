---
source: /Users/byron/projects/immunos81/utils/weka_exporter.py
relative: immunos81/utils/weka_exporter.py
generated_at: 2025-12-23 10:28
---

```python
"""
ARFF Exporter for Weka Integration

Exports Immunos-81 antigens to Weka's ARFF (Attribute-Relation File Format).
"""

from typing import List, TextIO
from ..core.antigen import Antigen
from ..core.amino_acid_library import AminoAcidLibrary, DataType


class ARFFExporter:
    """
    Exports Immunos-81 data to Weka ARFF format.

    ARFF format structure:
    - @relation <name>
    - @attribute declarations
    - @data
    - data rows
    """

    def __init__(self, library: AminoAcidLibrary):
        """
        Initialize exporter with amino acid library.

        Args:
            library: Amino acid library defining variables
        """
        self.library = library

    def export(self, antigens: List[Antigen], output_path: str, relation_name: str = "cleveland"):
        """
        Export antigens to ARFF file.

        Args:
            antigens: List of antigens to export
            output_path: Path to output ARFF file
            relation_name: Name for the relation (dataset name)
        """
        with open(output_path, 'w') as f:
            self._write_header(f, antigens, relation_name)
            self._write_data(f, antigens)

        print(f"Exported {len(antigens)} instances to {output_path}")

    def _write_header(self, f: TextIO, antigens: List[Antigen], relation_name: str):
        """Write ARFF header with relation and attribute declarations"""
        # Relation
        f.write(f"@relation {relation_name}\n\n")

        # Get variable sequence from first antigen
        if not antigens:
            raise ValueError("Cannot export empty antigen list")

        first_antigen = antigens[0]

        # Attribute declarations
        for var_num in first_antigen.variable_sequence:
            var_def = self.library.get_by_number(var_num)

            if var_def.data_type == DataType.NUMERIC:
                f.write(f"@attribute {var_def.name} numeric\n")

            elif var_def.data_type == DataType.NOMINAL:
                # Get all possible values from the data
                values = self._get_nominal_values(antigens, var_num)
                values_str = ','.join(str(v) for v in sorted(values))
                f.write(f"@attribute {var_def.name} {{{values_str}}}\n")

            elif var_def.data_type == DataType.ORDINAL:
                # Treat ordinal as nominal in ARFF
                if var_def.valid_values:
                    values_str = ','.join(str(v) for v in var_def.valid_values)
                else:
                    # Get from data
                    values = self._get_nominal_values(antigens, var_num)
                    values_str = ','.join(str(v) for v in sorted(values))
                f.write(f"@attribute {var_def.name} {{{values_str}}}\n")

        # Class attribute
        class_labels = sorted(set(a.class_label for a in antigens if a.class_label))
        class_str = ','.join(class_labels)
        f.write(f"@attribute class {{{class_str}}}\n")

        f.write("\n@data\n")

    def _write_data(self, f: TextIO, antigens: List[Antigen]):
        """Write data rows"""
        for antigen in antigens:
            row = []

            # Write attribute values
            for i, var_num in enumerate(antigen.variable_sequence):
                value = antigen.values[i]

                if value is None:
                    row.append('?')
                else:
                    var_def = self.library.get_by_number(var_num)

                    if var_def.data_type == DataType.NUMERIC:
                        row.append(str(value))
                    else:
                        # Nominal/ordinal - quote if contains special characters
                        value_str = str(value)
                        if ',' in value_str or ' ' in value_str or "'" in value_str:
                            value_str = f"'{value_str}'"
                        row.append(value_str)

            # Write class label
            row.append(antigen.class_label if antigen.class_label else '?')

            f.write(','.join(row) + '\n')

    def _get_nominal_values(self, antigens: List[Antigen], var_num: int) -> set:
        """Get all distinct non-null values for a nominal/ordinal variable"""
        values = set()

        for antigen in antigens:
            # Find position of this variable
            try:
                idx = antigen.variable_sequence.index(var_num)
                value = antigen.values[idx]
                if value is not None:
                    values.add(value)
            except ValueError:
                continue

        return values


def export_train_test_split(train_data: List[Antigen],
                           test_data: List[Antigen],
                           library: AminoAcidLibrary,
                           train_path: str,
                           test_path: str):
    """
    Export train and test sets to separate ARFF files.

    Args:
        train_data: Training antigens
        test_data: Test antigens
        library: Amino acid library
        train_path: Path for training ARFF file
        test_path: Path for test ARFF file
    """
    exporter = ARFFExporter(library)

    print("Exporting training set...")
    exporter.export(train_data, train_path, relation_name="cleveland_train")

    print("Exporting test set...")
    exporter.export(test_data, test_path, relation_name="cleveland_test")

    print(f"\nExported split:")
    print(f"  Train: {len(train_data)} instances -> {train_path}")
    print(f"  Test:  {len(test_data)} instances -> {test_path}")

```
