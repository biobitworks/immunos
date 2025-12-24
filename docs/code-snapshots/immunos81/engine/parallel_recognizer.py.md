---
source: /Users/byron/projects/immunos81/engine/parallel_recognizer.py
relative: immunos81/engine/parallel_recognizer.py
generated_at: 2025-12-23 10:28
---

```python
"""
Parallel Recognition Engine for Immunos-81

Implements parallel affinity calculations for faster recognition,
inspired by Immunos 99 optimizations.
"""

from typing import Dict, List
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing as mp
from ..core.tcell import TCell
from ..core.antigen import Antigen
from ..core.clone import Clone
from .recognizer import Recognizer, RecognitionResult, RecognitionStrategy


class ParallelRecognizer(Recognizer):
    """
    Parallelized version of Immunos-81 recognizer.

    Key improvement: Affinity calculations are distributed across a thread pool,
    significantly reducing recognition time for large datasets.
    """

    def __init__(self, tcells: Dict[str, TCell],
                 strategy: RecognitionStrategy = RecognitionStrategy.SHA,
                 n_workers: int = None):
        """
        Initialize parallel recognizer.

        Args:
            tcells: Dictionary of trained T cells
            strategy: Recognition strategy (SHA or RHA)
            n_workers: Number of worker threads (defaults to CPU count)
        """
        super().__init__(tcells, strategy)
        self.n_workers = n_workers or mp.cpu_count()

    def recognize(self, antigen: Antigen) -> RecognitionResult:
        """
        Recognize an antigen with parallel affinity calculations.

        Args:
            antigen: Antigen to classify

        Returns:
            RecognitionResult
        """
        # Collect avidity scores using parallel affinity calculations
        class_avidities = {}
        best_clone_per_class = {}

        # Process each T cell's clones in parallel
        with ThreadPoolExecutor(max_workers=self.n_workers) as executor:
            # Submit affinity calculation tasks for all T cells
            futures = {}
            for class_label, tcell in self.tcells.items():
                if tcell.matches_structure(antigen):
                    future = executor.submit(self._calculate_best_avidity, tcell, antigen)
                    futures[future] = class_label

            # Collect results
            for future in as_completed(futures):
                class_label = futures[future]
                result = future.result()
                if result is not None:
                    clone, avidity = result
                    class_avidities[class_label] = avidity
                    best_clone_per_class[class_label] = clone

        # No matches
        if not class_avidities:
            return RecognitionResult(
                predicted_class=None,
                confidence=0.0,
                winner_clone=None,
                all_avidities={},
                is_certain=False,
                strategy=self.strategy
            )

        # Apply recognition strategy (SHA or RHA)
        if self.strategy == RecognitionStrategy.SHA:
            return self._apply_sha(class_avidities, best_clone_per_class)
        else:
            return self._apply_rha(class_avidities, best_clone_per_class)

    def _calculate_best_avidity(self, tcell: TCell, antigen: Antigen):
        """
        Calculate best avidity for a T cell's clones.

        This is executed in parallel for each T cell.

        Args:
            tcell: T cell to evaluate
            antigen: Antigen to evaluate against

        Returns:
            Tuple of (best_clone, best_avidity) or None
        """
        if not tcell.clones:
            return None

        best_clone = None
        best_avidity = -1.0

        for clone in tcell.clones:
            avidity = clone.calculate_avidity(antigen)
            if avidity > best_avidity:
                best_avidity = avidity
                best_clone = clone

        if best_clone is None:
            return None

        return (best_clone, best_avidity)

    def batch_recognize(self, antigens: List[Antigen]) -> List[RecognitionResult]:
        """
        Recognize a batch of antigens in parallel.

        This provides maximum speedup by parallelizing across both
        antigens and affinity calculations.

        Args:
            antigens: List of antigens to classify

        Returns:
            List of RecognitionResults
        """
        results = [None] * len(antigens)

        with ThreadPoolExecutor(max_workers=self.n_workers) as executor:
            # Submit recognition tasks for all antigens
            futures = {
                executor.submit(self.recognize, antigen): idx
                for idx, antigen in enumerate(antigens)
            }

            # Collect results
            for future in as_completed(futures):
                idx = futures[future]
                results[idx] = future.result()

        return results

```
