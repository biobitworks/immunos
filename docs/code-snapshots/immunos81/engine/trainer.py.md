---
source: /Users/byron/projects/immunos81/engine/trainer.py
relative: immunos81/engine/trainer.py
generated_at: 2025-12-23 10:28
---

```python
"""
Training Engine for Immunos-81

The Trainer handles the learning phase where T cells, B cells, and clones are
created from training data.
"""

from typing import List, Dict, Set
from collections import defaultdict
from ..core.amino_acid_library import AminoAcidLibrary
from ..core.antigen import Antigen
from ..core.tcell import TCell
from ..core.bcell import BCell
from ..core.clone import Clone


class Trainer:
    """
    Handles the learning phase of Immunos-81.

    The training process:
    1. Creates T cells for each class
    2. Generates B cells for each training instance
    3. Groups B cells into clones
    4. Associates clones with their T cells
    """

    def __init__(self, library: AminoAcidLibrary):
        """
        Initialize the trainer.

        Args:
            library: Amino Acid Library containing variable definitions
        """
        self.library = library
        self.tcells: Dict[str, TCell] = {}
        self.bcells: List[BCell] = []
        self.clones: List[Clone] = []

    def train(self, training_data: List[Antigen]) -> Dict[str, TCell]:
        """
        Train the Immunos-81 system on a set of training antigens.

        Args:
            training_data: List of training antigens with class labels

        Returns:
            Dictionary mapping class labels to T cells

        Raises:
            ValueError: If training data is empty or contains antigens without labels
        """
        if not training_data:
            raise ValueError("Training data cannot be empty")

        # Verify all antigens have class labels
        for i, antigen in enumerate(training_data):
            if antigen.class_label is None:
                raise ValueError(f"Training antigen at index {i} has no class label")

        # Step 1: Create T cells for each class
        print(f"Training on {len(training_data)} instances...")
        self._create_tcells(training_data)
        print(f"Created {len(self.tcells)} T cells for classes: {list(self.tcells.keys())}")

        # Step 2: Create B cells for each training instance
        self._create_bcells(training_data)
        print(f"Created {len(self.bcells)} B cells")

        # Step 3: Group B cells into clones
        self._create_clones()
        print(f"Created {len(self.clones)} clones")

        # Step 4: Associate clones with T cells
        self._associate_clones()
        print("Associated clones with T cells")

        # Print summary
        for class_label, tcell in self.tcells.items():
            print(f"  {class_label}: {len(tcell.clones)} clones")

        return self.tcells

    def _create_tcells(self, training_data: List[Antigen]):
        """Create a T cell for each class in the training data"""
        # Get unique class labels
        class_labels = set(antigen.class_label for antigen in training_data)

        # Create T cell for each class
        for class_label in class_labels:
            tcell = TCell(class_label=class_label, library=self.library)

            # Learn primary structure from first antigen of this class
            class_antigens = [a for a in training_data if a.class_label == class_label]
            if class_antigens:
                tcell.learn_structure(class_antigens[0])

            self.tcells[class_label] = tcell

    def _create_bcells(self, training_data: List[Antigen]):
        """Create a B cell for each training instance"""
        for i, antigen in enumerate(training_data):
            bcell = BCell.from_antigen(antigen, identifier=f"bcell_{i}")
            self.bcells.append(bcell)

    def _create_clones(self):
        """
        Group B cells into clones.

        For now, we use a simple grouping strategy: create one clone per class.
        In a more sophisticated implementation, we could cluster B cells based on
        similarity of their antigen patterns.
        """
        # Group B cells by class
        bcells_by_class = defaultdict(list)
        for bcell in self.bcells:
            bcells_by_class[bcell.class_label].append(bcell)

        # Create one clone per class
        for class_label, bcells_list in bcells_by_class.items():
            clone = Clone(
                class_label=class_label,
                bcells=bcells_list,
                identifier=f"clone_{class_label}"
            )
            self.clones.append(clone)

    def _associate_clones(self):
        """Associate clones with their corresponding T cells"""
        for clone in self.clones:
            tcell = self.tcells.get(clone.class_label)
            if tcell:
                tcell.add_clone(clone)

    def add_instance(self, antigen: Antigen):
        """
        Add a new training instance online (without full retraining).

        This enables incremental learning - a key feature of Immunos-81.

        Args:
            antigen: New training antigen to add

        Raises:
            ValueError: If antigen has no class label
        """
        if antigen.class_label is None:
            raise ValueError("Cannot add antigen without class label")

        # Create or get T cell for this class
        if antigen.class_label not in self.tcells:
            tcell = TCell(class_label=antigen.class_label, library=self.library)
            tcell.learn_structure(antigen)
            self.tcells[antigen.class_label] = tcell
        else:
            # Verify structure matches
            tcell = self.tcells[antigen.class_label]
            if not tcell.matches_structure(antigen):
                raise ValueError(
                    f"New antigen structure does not match existing structure "
                    f"for class {antigen.class_label}"
                )

        # Create new B cell
        bcell = BCell.from_antigen(antigen, identifier=f"bcell_{len(self.bcells)}")
        self.bcells.append(bcell)

        # Add to appropriate clone (or create new clone)
        # For simplicity, add to existing clone for this class
        class_clones = [c for c in self.clones if c.class_label == antigen.class_label]
        if class_clones:
            class_clones[0].add_bcell(bcell)
        else:
            # Create new clone
            clone = Clone(
                class_label=antigen.class_label,
                bcells=[bcell],
                identifier=f"clone_{antigen.class_label}_{len(self.clones)}"
            )
            self.clones.append(clone)
            tcell.add_clone(clone)

    def get_tcells(self) -> Dict[str, TCell]:
        """Get all trained T cells"""
        return self.tcells

    def get_statistics(self) -> Dict[str, any]:
        """Get training statistics"""
        stats = {
            "num_classes": len(self.tcells),
            "num_bcells": len(self.bcells),
            "num_clones": len(self.clones),
            "classes": {}
        }

        for class_label, tcell in self.tcells.items():
            class_bcells = [b for b in self.bcells if b.class_label == class_label]
            stats["classes"][class_label] = {
                "num_bcells": len(class_bcells),
                "num_clones": len(tcell.clones),
                "structure": tcell.get_structure_info()
            }

        return stats

```
