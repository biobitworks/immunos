---
source: /Users/byron/projects/immunos81/engine/enhanced_recognizer.py
relative: immunos81/engine/enhanced_recognizer.py
generated_at: 2025-12-23 10:28
---

```python
"""
Enhanced Recognition Engine for Immunos-81

Implements activation/inhibition mechanism inspired by amorphous computing
to improve classification accuracy.
"""

from typing import Dict, List, Tuple, Optional
from enum import Enum
from dataclasses import dataclass
import math
from ..core.tcell import TCell
from ..core.antigen import Antigen


class RecognitionStrategy(Enum):
    """Recognition strategies for classification"""
    SHA = "sha"  # Simple Highest Avidity
    RHA = "rha"  # Relative Highest Avidity


@dataclass
class RecognitionResult:
    """Result of classifying an antigen"""
    predicted_class: Optional[str]
    actual_class: Optional[str]
    confidence: float
    avidity_scores: Dict[str, float]
    is_uncertain: bool = False
    winner_tcell: Optional[str] = None
    second_tcell: Optional[str] = None


class EnhancedRecognizer:
    """
    Enhanced recognizer with activation/inhibition mechanism.

    Key improvements:
    1. Lateral inhibition between competing clones
    2. Activation boost for strongly matching clones
    3. Confidence-based uncertainty detection
    """

    # Inhibition parameters (from amorphous computing)
    INHIBITION_STRENGTH = 0.15  # 15% suppression from competitors
    ACTIVATION_THRESHOLD = 0.7   # Threshold for activation boost
    RHA_THRESHOLD = 0.05         # 5% threshold for RHA strategy

    def __init__(self, tcells: Dict[str, TCell],
                 strategy: RecognitionStrategy = RecognitionStrategy.SHA,
                 enable_inhibition: bool = True):
        """
        Initialize the enhanced recognizer.

        Args:
            tcells: Dictionary of trained T cells
            strategy: Recognition strategy (SHA or RHA)
            enable_inhibition: Enable lateral inhibition mechanism
        """
        self.tcells = tcells
        self.strategy = strategy
        self.enable_inhibition = enable_inhibition

    def recognize(self, test_antigen: Antigen) -> RecognitionResult:
        """
        Recognize/classify a test antigen using enhanced mechanisms.

        Args:
            test_antigen: Antigen to classify

        Returns:
            RecognitionResult with prediction and confidence
        """
        # Step 1: Calculate raw avidity scores for all T cells
        raw_avidity_scores = {}
        for class_label, tcell in self.tcells.items():
            # Check structural compatibility
            if not tcell.matches_structure(test_antigen):
                raw_avidity_scores[class_label] = 0.0
                continue

            # Get best avidity from all clones
            best_avidity = 0.0
            for clone in tcell.clones:
                avidity = clone.calculate_avidity(test_antigen)
                if avidity > best_avidity:
                    best_avidity = avidity

            raw_avidity_scores[class_label] = best_avidity

        # Step 2: Apply activation/inhibition mechanism
        if self.enable_inhibition:
            adjusted_scores = self._apply_activation_inhibition(raw_avidity_scores)
        else:
            adjusted_scores = raw_avidity_scores

        # Step 3: Apply recognition strategy
        if self.strategy == RecognitionStrategy.SHA:
            result = self._apply_sha(test_antigen, adjusted_scores)
        else:
            result = self._apply_rha(test_antigen, adjusted_scores)

        return result

    def _apply_activation_inhibition(self, avidity_scores: Dict[str, float]) -> Dict[str, float]:
        """
        Apply lateral inhibition and activation boost.

        Inspired by amorphous computing and neural lateral inhibition:
        - Strong competitors inhibit weaker ones
        - Very strong responses get activation boost
        - Creates sharper decision boundaries
        """
        if not avidity_scores:
            return avidity_scores

        # Find max and total avidity
        max_avidity = max(avidity_scores.values()) if avidity_scores else 0.0
        total_avidity = sum(avidity_scores.values())

        if total_avidity == 0:
            return avidity_scores

        adjusted_scores = {}

        for class_label, avidity in avidity_scores.items():
            # Normalize avidity
            normalized_avidity = avidity / total_avidity if total_avidity > 0 else 0.0

            # Calculate inhibition from competitors
            competitor_strength = (total_avidity - avidity) / total_avidity
            inhibition = competitor_strength * self.INHIBITION_STRENGTH

            # Calculate activation boost if this is a strong response
            activation = 0.0
            if normalized_avidity > self.ACTIVATION_THRESHOLD:
                activation = (normalized_avidity - self.ACTIVATION_THRESHOLD) * 0.2

            # Apply adjustments
            adjusted = avidity * (1 - inhibition + activation)
            adjusted_scores[class_label] = max(0.0, adjusted)

        return adjusted_scores

    def _apply_sha(self, test_antigen: Antigen,
                   avidity_scores: Dict[str, float]) -> RecognitionResult:
        """
        Apply Simple Highest Avidity (SHA) strategy.

        Winner is T cell with highest avidity.
        """
        if not avidity_scores:
            return RecognitionResult(
                predicted_class=None,
                actual_class=test_antigen.class_label,
                confidence=0.0,
                avidity_scores={},
                is_uncertain=True
            )

        # Find winner (highest avidity)
        winner_class = max(avidity_scores.items(), key=lambda x: x[1])[0]
        winner_avidity = avidity_scores[winner_class]

        # Calculate confidence (winner's proportion of total avidity)
        total_avidity = sum(avidity_scores.values())
        confidence = winner_avidity / total_avidity if total_avidity > 0 else 0.0

        # Find second place for analysis
        sorted_classes = sorted(avidity_scores.items(), key=lambda x: x[1], reverse=True)
        second_class = sorted_classes[1][0] if len(sorted_classes) > 1 else None

        return RecognitionResult(
            predicted_class=winner_class,
            actual_class=test_antigen.class_label,
            confidence=confidence,
            avidity_scores=avidity_scores,
            is_uncertain=False,
            winner_tcell=winner_class,
            second_tcell=second_class
        )

    def _apply_rha(self, test_antigen: Antigen,
                   avidity_scores: Dict[str, float]) -> RecognitionResult:
        """
        Apply Relative Highest Avidity (RHA) strategy.

        Winner must exceed second place by at least 5% threshold.
        Otherwise marked as uncertain.
        """
        if not avidity_scores:
            return RecognitionResult(
                predicted_class=None,
                actual_class=test_antigen.class_label,
                confidence=0.0,
                avidity_scores={},
                is_uncertain=True
            )

        # Sort by avidity
        sorted_classes = sorted(avidity_scores.items(), key=lambda x: x[1], reverse=True)

        if len(sorted_classes) < 2:
            # Only one class, use it
            winner_class = sorted_classes[0][0]
            return RecognitionResult(
                predicted_class=winner_class,
                actual_class=test_antigen.class_label,
                confidence=1.0,
                avidity_scores=avidity_scores,
                is_uncertain=False,
                winner_tcell=winner_class
            )

        # Get top two
        first_class, first_avidity = sorted_classes[0]
        second_class, second_avidity = sorted_classes[1]

        # Calculate relative difference
        total = first_avidity + second_avidity
        if total == 0:
            relative_diff = 0
        else:
            relative_diff = (first_avidity - second_avidity) / total

        # Apply RHA threshold
        if relative_diff >= self.RHA_THRESHOLD:
            # Clear winner
            confidence = relative_diff
            is_uncertain = False
            predicted_class = first_class
        else:
            # Too close, mark as uncertain
            confidence = relative_diff
            is_uncertain = True
            predicted_class = None

        return RecognitionResult(
            predicted_class=predicted_class,
            actual_class=test_antigen.class_label,
            confidence=confidence,
            avidity_scores=avidity_scores,
            is_uncertain=is_uncertain,
            winner_tcell=first_class,
            second_tcell=second_class
        )

    def batch_recognize(self, test_data: List[Antigen]) -> List[RecognitionResult]:
        """
        Classify a batch of test antigens.

        Args:
            test_data: List of antigens to classify

        Returns:
            List of recognition results
        """
        return [self.recognize(antigen) for antigen in test_data]

    def evaluate(self, test_data: List[Antigen]) -> Dict[str, any]:
        """
        Evaluate classification performance on test data.

        Args:
            test_data: List of test antigens with known labels

        Returns:
            Dictionary of evaluation metrics
        """
        results = self.batch_recognize(test_data)

        # Overall metrics
        total = len(results)
        correct = sum(1 for r in results if r.predicted_class == r.actual_class)
        uncertain = sum(1 for r in results if r.is_uncertain)

        # Certain-only metrics (for RHA)
        certain_results = [r for r in results if not r.is_uncertain]
        certain_total = len(certain_results)
        certain_correct = sum(1 for r in certain_results
                            if r.predicted_class == r.actual_class)

        metrics = {
            "total": total,
            "correct": correct,
            "incorrect": total - correct,
            "uncertain_attempts": uncertain,
            "accuracy_overall": correct / total if total > 0 else 0.0,
            "accuracy_certain": certain_correct / certain_total if certain_total > 0 else 0.0,
            "certain_rate": certain_total / total if total > 0 else 0.0,
            "uncertain_rate": uncertain / total if total > 0 else 0.0,
            "strategy": self.strategy.value,
            "inhibition_enabled": self.enable_inhibition
        }

        # Per-class accuracy
        per_class = {}
        for class_label in self.tcells.keys():
            class_results = [r for r in results if r.actual_class == class_label]
            if class_results:
                class_correct = sum(1 for r in class_results
                                  if r.predicted_class == r.actual_class)
                per_class[class_label] = {
                    "total": len(class_results),
                    "correct": class_correct,
                    "accuracy": class_correct / len(class_results)
                }

        metrics["per_class"] = per_class

        return metrics

    def set_strategy(self, strategy: RecognitionStrategy):
        """Change the recognition strategy"""
        self.strategy = strategy

    def enable_lateral_inhibition(self, enable: bool):
        """Enable or disable lateral inhibition"""
        self.enable_inhibition = enable

```
