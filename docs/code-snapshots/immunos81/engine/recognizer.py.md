---
source: /Users/byron/projects/immunos81/engine/recognizer.py
relative: immunos81/engine/recognizer.py
generated_at: 2025-12-23 10:28
---

```python
"""
Recognition Engine for Immunos-81

The Recognizer handles classification of unknown antigens using trained T cells
and clones. Supports both SHA (Simple Highest Avidity) and RHA (Relative Highest
Avidity) recognition strategies.
"""

from typing import Dict, List, Tuple, Optional
from enum import Enum
from dataclasses import dataclass
from ..core.tcell import TCell
from ..core.antigen import Antigen
from ..core.clone import Clone


class RecognitionStrategy(Enum):
    """Recognition strategy options"""
    SHA = "sha"  # Simple Highest Avidity
    RHA = "rha"  # Relative Highest Avidity (5% threshold)


@dataclass
class RecognitionResult:
    """
    Result of antigen recognition.

    Attributes:
        predicted_class: The predicted class label
        confidence: Confidence score (avidity of winning clone)
        winner_clone: The clone that won
        all_avidities: Dictionary mapping class labels to their best avidity scores
        is_certain: Whether the prediction is certain (for RHA strategy)
        strategy: The recognition strategy used
    """
    predicted_class: Optional[str]
    confidence: float
    winner_clone: Optional[Clone]
    all_avidities: Dict[str, float]
    is_certain: bool = True
    strategy: RecognitionStrategy = RecognitionStrategy.SHA

    def __repr__(self) -> str:
        certain_str = "certain" if self.is_certain else "uncertain"
        return (f"RecognitionResult(class={self.predicted_class}, "
                f"confidence={self.confidence:.4f}, {certain_str})")


class Recognizer:
    """
    Handles classification of unknown antigens.

    The recognition process:
    1. Present antigen to all T cells
    2. T cells check for structure match
    3. Matching T cells' clones calculate avidity
    4. Winner determined by strategy (SHA or RHA)
    """

    def __init__(self, tcells: Dict[str, TCell], strategy: RecognitionStrategy = RecognitionStrategy.SHA):
        """
        Initialize the recognizer.

        Args:
            tcells: Dictionary mapping class labels to T cells
            strategy: Recognition strategy to use (SHA or RHA)
        """
        self.tcells = tcells
        self.strategy = strategy
        self.rha_threshold = 0.05  # 5% threshold for RHA

    def recognize(self, antigen: Antigen) -> RecognitionResult:
        """
        Recognize (classify) an unknown antigen.

        Args:
            antigen: Antigen to classify

        Returns:
            RecognitionResult with prediction and confidence

        The recognition process:
        - SHA: Winner is clone with highest avidity
        - RHA: Winner must have avidity > 1.05 * second_best, else "uncertain"
        """
        # Collect avidity scores from all matching T cells
        class_avidities = {}
        best_clone_per_class = {}

        for class_label, tcell in self.tcells.items():
            result = tcell.recognize(antigen)
            if result is not None:
                clone, avidity = result
                class_avidities[class_label] = avidity
                best_clone_per_class[class_label] = clone

        # If no matches, return uncertain result
        if not class_avidities:
            return RecognitionResult(
                predicted_class=None,
                confidence=0.0,
                winner_clone=None,
                all_avidities={},
                is_certain=False,
                strategy=self.strategy
            )

        # Apply recognition strategy
        if self.strategy == RecognitionStrategy.SHA:
            return self._apply_sha(class_avidities, best_clone_per_class)
        else:  # RHA
            return self._apply_rha(class_avidities, best_clone_per_class)

    def _apply_sha(self, class_avidities: Dict[str, float],
                   best_clone_per_class: Dict[str, Clone]) -> RecognitionResult:
        """
        Apply Simple Highest Avidity (SHA) strategy.

        Winner is the class with the highest avidity score.

        Args:
            class_avidities: Dictionary mapping class labels to avidity scores
            best_clone_per_class: Dictionary mapping class labels to best clones

        Returns:
            RecognitionResult
        """
        # Find class with highest avidity
        winner_class = max(class_avidities, key=class_avidities.get)
        winner_avidity = class_avidities[winner_class]
        winner_clone = best_clone_per_class[winner_class]

        return RecognitionResult(
            predicted_class=winner_class,
            confidence=winner_avidity,
            winner_clone=winner_clone,
            all_avidities=class_avidities,
            is_certain=True,
            strategy=RecognitionStrategy.SHA
        )

    def _apply_rha(self, class_avidities: Dict[str, float],
                   best_clone_per_class: Dict[str, Clone]) -> RecognitionResult:
        """
        Apply Relative Highest Avidity (RHA) strategy.

        Winner must have avidity > (1 + threshold) * second_best.
        If not, the result is "uncertain".

        Args:
            class_avidities: Dictionary mapping class labels to avidity scores
            best_clone_per_class: Dictionary mapping class labels to best clones

        Returns:
            RecognitionResult
        """
        # Sort classes by avidity (descending)
        sorted_classes = sorted(class_avidities.items(), key=lambda x: x[1], reverse=True)

        if len(sorted_classes) == 1:
            # Only one class matches - it wins
            winner_class, winner_avidity = sorted_classes[0]
            return RecognitionResult(
                predicted_class=winner_class,
                confidence=winner_avidity,
                winner_clone=best_clone_per_class[winner_class],
                all_avidities=class_avidities,
                is_certain=True,
                strategy=RecognitionStrategy.RHA
            )

        # Get top two classes
        first_class, first_avidity = sorted_classes[0]
        second_class, second_avidity = sorted_classes[1]

        # Check if winner is clearly better (5% threshold)
        threshold_avidity = second_avidity * (1.0 + self.rha_threshold)

        if first_avidity > threshold_avidity:
            # Clear winner
            return RecognitionResult(
                predicted_class=first_class,
                confidence=first_avidity,
                winner_clone=best_clone_per_class[first_class],
                all_avidities=class_avidities,
                is_certain=True,
                strategy=RecognitionStrategy.RHA
            )
        else:
            # Too close to call - uncertain
            return RecognitionResult(
                predicted_class=first_class,  # Still report best guess
                confidence=first_avidity,
                winner_clone=best_clone_per_class[first_class],
                all_avidities=class_avidities,
                is_certain=False,
                strategy=RecognitionStrategy.RHA
            )

    def batch_recognize(self, antigens: List[Antigen]) -> List[RecognitionResult]:
        """
        Recognize a batch of antigens.

        Args:
            antigens: List of antigens to classify

        Returns:
            List of RecognitionResults
        """
        return [self.recognize(antigen) for antigen in antigens]

    def evaluate(self, test_data: List[Antigen]) -> Dict[str, any]:
        """
        Evaluate recognition performance on test data.

        Args:
            test_data: List of antigens with known class labels

        Returns:
            Dictionary containing evaluation metrics
        """
        results = self.batch_recognize(test_data)

        # Calculate metrics
        total = len(test_data)
        correct = 0
        attempted = 0
        correct_certain = 0
        certain_attempts = 0

        for antigen, result in zip(test_data, results):
            if result.predicted_class is not None:
                attempted += 1
                if result.predicted_class == antigen.class_label:
                    correct += 1

            if result.is_certain and result.predicted_class is not None:
                certain_attempts += 1
                if result.predicted_class == antigen.class_label:
                    correct_certain += 1

        metrics = {
            "strategy": self.strategy.value,
            "total": total,
            "attempted": attempted,
            "not_attempted": total - attempted,
            "correct": correct,
            "incorrect": attempted - correct,
            "accuracy_overall": correct / total if total > 0 else 0.0,
            "accuracy_attempted": correct / attempted if attempted > 0 else 0.0,
        }

        # RHA-specific metrics
        if self.strategy == RecognitionStrategy.RHA:
            metrics.update({
                "certain_attempts": certain_attempts,
                "uncertain_attempts": attempted - certain_attempts,
                "correct_certain": correct_certain,
                "accuracy_certain": correct_certain / certain_attempts if certain_attempts > 0 else 0.0
            })

        return metrics

    def set_strategy(self, strategy: RecognitionStrategy):
        """Change the recognition strategy"""
        self.strategy = strategy

    def set_rha_threshold(self, threshold: float):
        """
        Set the RHA threshold.

        Args:
            threshold: Threshold value (default 0.05 = 5%)
        """
        self.rha_threshold = threshold

```
