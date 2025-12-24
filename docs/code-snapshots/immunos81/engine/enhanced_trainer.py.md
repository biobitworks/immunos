---
source: /Users/byron/projects/immunos81/engine/enhanced_trainer.py
relative: immunos81/engine/enhanced_trainer.py
generated_at: 2025-12-23 10:28
---

```python
"""
Enhanced Training Engine for Immunos-81

Implements accuracy improvements from AIS and amorphous computing research:
- Multiple clones per class (clustering-based)
- Adaptive clonal factor (optimal 0.4)
- Affinity maturation via hypermutation
- Negative selection for outlier detection
- Activation/inhibition between competing clones
"""

from typing import List, Dict, Set, Tuple
from collections import defaultdict
import random
import math
from ..core.amino_acid_library import AminoAcidLibrary, DataType
from ..core.antigen import Antigen
from ..core.tcell import TCell
from ..core.bcell import BCell
from ..core.clone import Clone


class EnhancedTrainer:
    """
    Enhanced trainer with accuracy improvements from AIS research.

    Key improvements over basic Trainer:
    1. Multiple clones per class (3-5 clones via clustering)
    2. Adaptive clonal factor (0.4 optimal value)
    3. Affinity maturation via hypermutation
    4. Negative selection to filter outliers
    5. Activation/inhibition between competing clones
    """

    # Conservative parameters - only apply proven improvements
    MIN_CLONES_PER_CLASS = 2
    MAX_CLONES_PER_CLASS = 4
    HYPERMUTATION_RATE = 0.05  # Reduced mutation rate
    MATURATION_ITERATIONS = 2  # Fewer iterations
    NEGATIVE_SELECTION_THRESHOLD = 0.01
    MIN_CLONES_PER_CLASS_AFTER_SELECTION = 2  # Keep at least 2 clones per class

    def __init__(self, library: AminoAcidLibrary,
                 clones_per_class: int = None,
                 enable_maturation: bool = True,
                 enable_negative_selection: bool = True):
        """
        Initialize the enhanced trainer.

        Args:
            library: Amino Acid Library containing variable definitions
            clones_per_class: Number of clones per class (default: adaptive 3-5)
            enable_maturation: Enable affinity maturation via hypermutation
            enable_negative_selection: Enable negative selection filtering
        """
        self.library = library
        self.tcells: Dict[str, TCell] = {}
        self.bcells: List[BCell] = []
        self.clones: List[Clone] = []
        self.clones_per_class = clones_per_class
        self.enable_maturation = enable_maturation
        self.enable_negative_selection = enable_negative_selection

    def train(self, training_data: List[Antigen]) -> Dict[str, TCell]:
        """
        Train the Immunos-81 system with enhanced algorithms.

        Args:
            training_data: List of training antigens with class labels

        Returns:
            Dictionary mapping class labels to T cells
        """
        if not training_data:
            raise ValueError("Training data cannot be empty")

        # Verify all antigens have class labels
        for i, antigen in enumerate(training_data):
            if antigen.class_label is None:
                raise ValueError(f"Training antigen at index {i} has no class label")

        print(f"Enhanced Training on {len(training_data)} instances...")

        # Step 1: Create T cells for each class
        self._create_tcells(training_data)
        print(f"Created {len(self.tcells)} T cells for classes: {list(self.tcells.keys())}")

        # Step 2: Create B cells for each training instance
        self._create_bcells(training_data)
        print(f"Created {len(self.bcells)} B cells")

        # Step 3: Apply affinity maturation (hypermutation)
        if self.enable_maturation:
            self._apply_affinity_maturation(training_data)
            print(f"Applied affinity maturation ({self.MATURATION_ITERATIONS} iterations)")

        # Step 4: Create multiple clones per class via clustering
        self._create_clones_clustered()
        print(f"Created {len(self.clones)} clones via clustering")

        # Step 5: Negative selection (filter outliers)
        if self.enable_negative_selection:
            initial_clones = len(self.clones)
            self._apply_negative_selection(training_data)
            filtered = initial_clones - len(self.clones)
            print(f"Negative selection: filtered {filtered}/{initial_clones} clones")

        # Step 6: Associate clones with T cells
        self._associate_clones()
        print("Associated clones with T cells")

        # Print summary
        for class_label, tcell in self.tcells.items():
            print(f"  {class_label}: {len(tcell.clones)} clones")

        return self.tcells

    def _create_tcells(self, training_data: List[Antigen]):
        """Create a T cell for each class in the training data"""
        class_labels = set(antigen.class_label for antigen in training_data)

        for class_label in class_labels:
            tcell = TCell(class_label=class_label, library=self.library)
            class_antigens = [a for a in training_data if a.class_label == class_label]
            if class_antigens:
                tcell.learn_structure(class_antigens[0])
            self.tcells[class_label] = tcell

    def _create_bcells(self, training_data: List[Antigen]):
        """Create a B cell for each training instance"""
        for i, antigen in enumerate(training_data):
            bcell = BCell.from_antigen(antigen, identifier=f"bcell_{i}")
            self.bcells.append(bcell)

    def _apply_affinity_maturation(self, training_data: List[Antigen]):
        """
        Apply affinity maturation via hypermutation.

        Creates mutated variants of B cells and selects those with better
        recognition capabilities. Based on clonal selection algorithm improvements.
        """
        for iteration in range(self.MATURATION_ITERATIONS):
            new_bcells = []

            for bcell in self.bcells:
                # Create mutated variant
                mutated = self._hypermutate_bcell(bcell)

                # Evaluate fitness: average affinity to same-class antigens
                same_class = [a for a in training_data if a.class_label == bcell.class_label]

                if same_class:
                    # Calculate average affinity for original and mutated
                    orig_fitness = self._calculate_fitness(bcell, same_class)
                    mut_fitness = self._calculate_fitness(mutated, same_class)

                    # Select the better one
                    if mut_fitness > orig_fitness:
                        new_bcells.append(mutated)
                    else:
                        new_bcells.append(bcell)
                else:
                    new_bcells.append(bcell)

            self.bcells = new_bcells

    def _hypermutate_bcell(self, bcell: BCell) -> BCell:
        """
        Create a hypermutated variant of a B cell.

        Applies small random perturbations to numeric values.
        For nominal/ordinal, occasionally swaps to nearby values.
        """
        mutated_values = bcell.antigen.values.copy()

        for i, (var_num, value) in enumerate(zip(bcell.antigen.variable_sequence, mutated_values)):
            # Skip None values (missing data)
            if value is None:
                continue

            if random.random() < self.HYPERMUTATION_RATE:
                var_def = self.library.get_by_number(var_num)

                if var_def.data_type == DataType.NUMERIC:
                    # Add small Gaussian noise
                    noise = random.gauss(0, 0.1 * abs(value) if value != 0 else 0.1)
                    mutated_values[i] = value + noise

                elif var_def.data_type == DataType.ORDINAL:
                    # Shift to adjacent ordinal value
                    if var_def.valid_values and value in var_def.valid_values:
                        idx = var_def.valid_values.index(value)
                        if random.random() < 0.5 and idx > 0:
                            mutated_values[i] = var_def.valid_values[idx - 1]
                        elif idx < len(var_def.valid_values) - 1:
                            mutated_values[i] = var_def.valid_values[idx + 1]

        # Create new antigen with mutated values
        mutated_antigen = Antigen(
            values=mutated_values,
            variable_sequence=bcell.antigen.variable_sequence.copy(),
            library=self.library,
            class_label=bcell.class_label,
            identifier=f"{bcell.identifier}_mut"
        )

        # Create new B cell
        return BCell.from_antigen(mutated_antigen, identifier=f"{bcell.identifier}_mut")

    def _calculate_fitness(self, bcell: BCell, antigens: List[Antigen]) -> float:
        """Calculate average affinity of a B cell to a set of antigens"""
        if not antigens:
            return 0.0

        total_affinity = 0.0
        count = 0

        for antigen in antigens:
            for var_num in antigen.variable_sequence:
                total_affinity += bcell.calculate_affinity(antigen, var_num)
                count += 1

        return total_affinity / count if count > 0 else 0.0

    def _create_clones_clustered(self):
        """
        Create multiple clones per class via clustering.

        Uses a simple distance-based clustering to group similar B cells
        into different clones. This creates specialized clones that recognize
        different patterns within each class.
        """
        # Group B cells by class
        bcells_by_class = defaultdict(list)
        for bcell in self.bcells:
            bcells_by_class[bcell.class_label].append(bcell)

        # For each class, cluster B cells
        for class_label, bcells_list in bcells_by_class.items():
            # Determine number of clones for this class
            n_clones = self.clones_per_class
            if n_clones is None:
                # Adaptive: 3-5 clones based on class size
                if len(bcells_list) < 10:
                    n_clones = self.MIN_CLONES_PER_CLASS
                else:
                    n_clones = min(self.MAX_CLONES_PER_CLASS,
                                  max(self.MIN_CLONES_PER_CLASS,
                                      len(bcells_list) // 20))

            # Simple k-means-like clustering
            clusters = self._cluster_bcells(bcells_list, n_clones)

            # Create clone for each cluster
            for i, cluster in enumerate(clusters):
                if cluster:  # Only create non-empty clones
                    clone = Clone(
                        class_label=class_label,
                        bcells=cluster,
                        identifier=f"clone_{class_label}_{i}"
                    )
                    self.clones.append(clone)

    def _cluster_bcells(self, bcells: List[BCell], n_clusters: int) -> List[List[BCell]]:
        """
        Cluster B cells using a simple distance-based approach.

        This is a lightweight k-means variant that groups similar B cells together.
        """
        if len(bcells) <= n_clusters:
            # If we have fewer B cells than clusters, each gets its own
            return [[bcell] for bcell in bcells]

        # Initialize clusters with random seeds
        random.shuffle(bcells)
        centroids = bcells[:n_clusters]

        # Run a few iterations of clustering
        for iteration in range(5):
            # Assign B cells to nearest centroid
            clusters = [[] for _ in range(n_clusters)]

            for bcell in bcells:
                min_dist = float('inf')
                best_cluster = 0

                for i, centroid in enumerate(centroids):
                    dist = self._bcell_distance(bcell, centroid)
                    if dist < min_dist:
                        min_dist = dist
                        best_cluster = i

                clusters[best_cluster].append(bcell)

            # Update centroids (pick medoid from each cluster)
            for i, cluster in enumerate(clusters):
                if cluster:
                    centroids[i] = self._find_medoid(cluster)

        return clusters

    def _bcell_distance(self, bcell1: BCell, bcell2: BCell) -> float:
        """
        Calculate distance between two B cells based on their antigen values.
        """
        dist = 0.0
        count = 0

        # Compare common variables (assume both have same variable sequence)
        for i, var_num in enumerate(bcell1.antigen.variable_sequence):
            if i >= len(bcell2.antigen.values):
                break

            val1 = bcell1.antigen.values[i]
            val2 = bcell2.antigen.values[i]

            # Skip if either value is None
            if val1 is None or val2 is None:
                continue

            var_def = self.library.get_by_number(var_num)

            if var_def.data_type == DataType.NUMERIC:
                # Normalized distance for numeric values
                if val1 != val2:
                    dist += abs(val1 - val2)
            elif var_def.data_type == DataType.NOMINAL:
                # Binary distance for nominal
                dist += 0.0 if val1 == val2 else 1.0
            elif var_def.data_type == DataType.ORDINAL:
                # Ordinal distance
                if var_def.valid_values:
                    try:
                        idx1 = var_def.valid_values.index(val1)
                        idx2 = var_def.valid_values.index(val2)
                        dist += abs(idx1 - idx2) / len(var_def.valid_values)
                    except ValueError:
                        dist += 1.0

            count += 1

        return dist / count if count > 0 else 0.0

    def _find_medoid(self, cluster: List[BCell]) -> BCell:
        """Find the medoid (most central B cell) of a cluster"""
        if len(cluster) == 1:
            return cluster[0]

        min_total_dist = float('inf')
        medoid = cluster[0]

        for candidate in cluster:
            total_dist = sum(self._bcell_distance(candidate, other)
                           for other in cluster)
            if total_dist < min_total_dist:
                min_total_dist = total_dist
                medoid = candidate

        return medoid

    def _apply_negative_selection(self, training_data: List[Antigen]):
        """
        Apply negative selection to filter out overly broad or weak clones.

        Removes clones that:
        1. Have very low affinity to their own class (weak)
        2. Are overly cross-reactive (recognize wrong class too strongly)

        Ensures at least MIN_CLONES_PER_CLASS_AFTER_SELECTION clones remain per class.
        """
        # Group clones by class
        clones_by_class = defaultdict(list)
        for clone in self.clones:
            clones_by_class[clone.class_label].append(clone)

        filtered_clones = []

        for class_label, class_clones in clones_by_class.items():
            # Get antigens from same and different classes
            same_class = [a for a in training_data if a.class_label == class_label]
            other_class = [a for a in training_data if a.class_label != class_label]

            # Score each clone
            clone_scores = []
            for clone in class_clones:
                # Calculate average affinity to own class
                own_affinity = 0.0
                if same_class:
                    own_affinity = sum(clone.calculate_avidity(a) for a in same_class[:10]) / min(10, len(same_class))

                # Calculate average affinity to other class
                cross_affinity = 0.0
                if other_class:
                    cross_affinity = sum(clone.calculate_avidity(a) for a in other_class[:10]) / min(10, len(other_class))

                # Score based on specificity (own affinity - cross affinity)
                specificity = own_affinity - cross_affinity
                clone_scores.append((clone, specificity, own_affinity, cross_affinity))

            # Sort by specificity (higher is better)
            clone_scores.sort(key=lambda x: x[1], reverse=True)

            # Keep the best clones, ensuring minimum per class
            kept = 0
            for clone, specificity, own_affinity, cross_affinity in clone_scores:
                # Always keep minimum clones per class
                if kept < self.MIN_CLONES_PER_CLASS_AFTER_SELECTION:
                    filtered_clones.append(clone)
                    kept += 1
                # Otherwise only keep if meets quality thresholds
                elif own_affinity >= self.NEGATIVE_SELECTION_THRESHOLD and cross_affinity < own_affinity * 0.9:
                    filtered_clones.append(clone)
                    kept += 1

        self.clones = filtered_clones

    def _associate_clones(self):
        """Associate clones with their corresponding T cells"""
        for clone in self.clones:
            tcell = self.tcells.get(clone.class_label)
            if tcell:
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
            "maturation_enabled": self.enable_maturation,
            "negative_selection_enabled": self.enable_negative_selection,
            "classes": {}
        }

        for class_label, tcell in self.tcells.items():
            class_bcells = [b for b in self.bcells if b.class_label == class_label]
            stats["classes"][class_label] = {
                "num_bcells": len(class_bcells),
                "num_clones": len(tcell.clones),
                "avg_clone_size": sum(c.size() for c in tcell.clones) / len(tcell.clones) if tcell.clones else 0,
                "structure": tcell.get_structure_info()
            }

        return stats

```
