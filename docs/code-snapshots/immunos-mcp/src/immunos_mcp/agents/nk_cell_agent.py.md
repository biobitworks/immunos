---
source: /Users/byron/projects/immunos-mcp/src/immunos_mcp/agents/nk_cell_agent.py
relative: immunos-mcp/src/immunos_mcp/agents/nk_cell_agent.py
generated_at: 2025-12-23 10:28
---

```python
"""
NK Cell Agent for IMMUNOS-MCP

Anomaly detection agent using Negative Selection Algorithm.
Inspired by Natural Killer cells and T-cell negative selection in the thymus.
"""

from typing import Any, Dict, List, Optional, Tuple
import time
import numpy as np
import random
import json
from pathlib import Path

from ..core.antigen import Antigen
from ..core.affinity import AffinityCalculator, DistanceMetric
from ..core.protocols import Detector, Pattern, AnomalyResult, AgentResponse


class NKCellAgent:
    """
    NK Cell Agent - Anomaly Detector via Negative Selection

    Mimics NK cell and thymic selection:
    - Trains on "self" (normal) data
    - Generates detectors that DON'T match self
    - Detects "non-self" (anomalies) without explicit training
    - Provides zero-shot threat detection
    """

    def __init__(self, agent_name: str = "nk_cell_001",
                 detection_threshold: float = 0.5,
                 num_detectors: int = 100):
        """
        Initialize NK Cell agent.

        Args:
            agent_name: Unique identifier
            detection_threshold: Threshold for anomaly detection (0-1)
            num_detectors: Number of negative selection detectors to generate
        """
        self.agent_name = agent_name
        self.detection_threshold = detection_threshold
        self.num_detectors = num_detectors
        self.self_patterns: List[Pattern] = []  # Normal (self) patterns
        self.detectors: List[Detector] = []  # Negative selection detectors
        self.affinity_calculator = AffinityCalculator(method="hybrid")
        self.embedding_cache: Dict[str, np.ndarray] = {}

    def train_on_self(self, normal_data: List[Antigen],
                      embeddings: Optional[List[np.ndarray]] = None):
        """
        Train on normal ("self") data.

        The negative selection algorithm learns what is normal,
        then detects deviations from normal as anomalies.

        Args:
            normal_data: List of normal (non-anomalous) antigens
            embeddings: Optional embeddings for normal data
        """
        print(f"[{self.agent_name}] Training on {len(normal_data)} normal (self) patterns...")

        # Store self patterns
        for i, antigen in enumerate(normal_data):
            pattern = Pattern(
                pattern_id=f"{self.agent_name}_self_{i}",
                class_label="self",  # All normal data is "self"
                example_data=antigen.data,
                embedding=embeddings[i].tolist() if embeddings and i < len(embeddings) else None,
                features=antigen.features,
                creation_time=time.time()
            )
            self.self_patterns.append(pattern)

            if embeddings and i < len(embeddings):
                self.embedding_cache[pattern.pattern_id] = embeddings[i]

        # Generate negative selection detectors
        self._generate_detectors()

        print(f"[{self.agent_name}] Generated {len(self.detectors)} detectors")

    def _generate_detectors(self):
        """
        Generate detectors using negative selection.

        Detectors are patterns that DON'T match self.
        They are generated randomly and filtered to ensure they don't match normal data.
        """
        if not self.self_patterns:
            print(f"[{self.agent_name}] Warning: No self patterns to generate detectors from")
            return

        detectors_generated = 0
        attempts = 0
        max_attempts = self.num_detectors * 10  # Prevent infinite loop

        while detectors_generated < self.num_detectors and attempts < max_attempts:
            attempts += 1

            # Generate random detector (in embedding space if available)
            if self.embedding_cache:
                # Random point in embedding space
                sample_embedding = random.choice(list(self.embedding_cache.values()))
                detector_embedding = sample_embedding + np.random.normal(0, 0.5, sample_embedding.shape)
                detector_embedding = detector_embedding / np.linalg.norm(detector_embedding)  # Normalize

                # Check if it matches any self pattern
                matches_self = False
                for pattern_id, self_embedding in self.embedding_cache.items():
                    similarity = 1 - DistanceMetric.cosine_distance(detector_embedding, self_embedding)
                    if similarity > self.detection_threshold:
                        matches_self = True
                        break

                if not matches_self:
                    # This detector doesn't match self - keep it
                    detector = Detector(
                        detector_id=f"{self.agent_name}_detector_{detectors_generated}",
                        self_patterns=self.self_patterns,
                        threshold=self.detection_threshold,
                        creation_time=time.time()
                    )
                    self.detectors.append(detector)
                    # Cache the detector embedding
                    self.embedding_cache[detector.detector_id] = detector_embedding
                    detectors_generated += 1

        if detectors_generated < self.num_detectors:
            print(f"[{self.agent_name}] Warning: Only generated {detectors_generated}/{self.num_detectors} detectors")

    def detect_novelty(self, antigen: Antigen,
                      antigen_embedding: Optional[np.ndarray] = None) -> AnomalyResult:
        """
        Detect if antigen is novel/anomalous using negative selection.

        An antigen is anomalous if it:
        1. Doesn't match self patterns (non-self)
        2. Matches one or more negative selection detectors

        Args:
            antigen: Antigen to test
            antigen_embedding: Optional embedding

        Returns:
            AnomalyResult with detection outcome
        """
        start_time = time.time()

        # Step 1: Check similarity to self patterns
        max_self_similarity = 0.0
        best_self_pattern = None

        for pattern in self.self_patterns:
            pattern_embedding = self.embedding_cache.get(pattern.pattern_id)
            affinity = self.affinity_calculator.calculate(
                antigen.data,
                pattern.example_data,
                embeddings1=antigen_embedding,
                embeddings2=pattern_embedding
            )
            if affinity.score > max_self_similarity:
                max_self_similarity = affinity.score
                best_self_pattern = pattern

        # Step 2: Check if matches detectors (non-self)
        detector_matches = []
        max_detector_similarity = 0.0

        if antigen_embedding is not None:
            for detector in self.detectors:
                detector_embedding = self.embedding_cache.get(detector.detector_id)
                if detector_embedding is not None:
                    similarity = 1 - DistanceMetric.cosine_distance(antigen_embedding, detector_embedding)
                    if similarity > self.detection_threshold:
                        detector_matches.append(detector.detector_id)
                        max_detector_similarity = max(max_detector_similarity, similarity)

        # Decision logic: Anomaly if low self-similarity OR high detector-similarity
        is_anomaly = (
            max_self_similarity < self.detection_threshold or
            len(detector_matches) > 0
        )

        # Calculate anomaly score
        anomaly_score = max(
            1.0 - max_self_similarity,  # Distance from self
            max_detector_similarity  # Similarity to non-self detectors
        )

        # Confidence based on how clear the decision is
        confidence = abs(max_self_similarity - max_detector_similarity)
        confidence = min(1.0, max(0.1, confidence))  # Clamp between 0.1 and 1.0

        explanation = self._generate_explanation(
            is_anomaly, max_self_similarity, max_detector_similarity,
            best_self_pattern, detector_matches
        )

        execution_time = time.time() - start_time

        return AnomalyResult(
            is_anomaly=is_anomaly,
            anomaly_score=anomaly_score,
            confidence=confidence,
            detected_patterns=detector_matches,
            explanation=explanation,
            detector_id=self.agent_name,
            metadata={
                "execution_time": execution_time,
                "num_self_patterns": len(self.self_patterns),
                "num_detectors": len(self.detectors),
                "max_self_similarity": max_self_similarity,
                "max_detector_similarity": max_detector_similarity,
                "detector_matches": len(detector_matches)
            }
        )

    def _generate_explanation(self, is_anomaly: bool,
                            max_self_sim: float,
                            max_detector_sim: float,
                            best_pattern: Optional[Pattern],
                            detector_matches: List[str]) -> str:
        """Generate human-readable explanation"""
        if is_anomaly:
            reasons = []
            if max_self_sim < self.detection_threshold:
                reasons.append(f"low similarity to normal patterns ({max_self_sim:.3f})")
            if detector_matches:
                reasons.append(f"matched {len(detector_matches)} non-self detectors ({max_detector_sim:.3f})")
            return f"Anomaly detected: {', '.join(reasons)}"
        else:
            return f"Normal pattern: high similarity to self ({max_self_sim:.3f})"

    def assess_danger(self, antigen: Antigen,
                     context_signals: Optional[Dict[str, float]] = None,
                     antigen_embedding: Optional[np.ndarray] = None) -> AnomalyResult:
        """
        Assess danger level using context signals (danger theory).

        Combines negative selection with contextual danger signals.

        Args:
            antigen: Antigen to assess
            context_signals: Dict of signal_type -> strength (0-1)
            antigen_embedding: Optional embedding

        Returns:
            AnomalyResult with danger assessment
        """
        # Base detection using negative selection
        base_result = self.detect_novelty(antigen, antigen_embedding)

        # Adjust based on context signals (danger theory)
        if context_signals:
            danger_adjustment = 0.0

            # PAMP signals (known bad patterns)
            if "pamp" in context_signals:
                danger_adjustment += 0.3 * context_signals["pamp"]

            # Danger signals (contextual warnings)
            if "danger" in context_signals:
                danger_adjustment += 0.4 * context_signals["danger"]

            # Safe signals (reduce anomaly score)
            if "safe" in context_signals:
                danger_adjustment -= 0.3 * context_signals["safe"]

            # Adjust anomaly score
            adjusted_score = max(0.0, min(1.0, base_result.anomaly_score + danger_adjustment))

            return AnomalyResult(
                is_anomaly=adjusted_score > 0.5,
                anomaly_score=adjusted_score,
                confidence=base_result.confidence,
                detected_patterns=base_result.detected_patterns,
                explanation=f"{base_result.explanation} (danger-adjusted: {danger_adjustment:+.3f})",
                detector_id=self.agent_name,
                metadata={
                    **base_result.metadata,
                    "danger_signals": context_signals
                }
            )

        return base_result

    def zero_shot_detect(self, antigen: Antigen, description: str) -> AnomalyResult:
        """
        Zero-shot anomaly detection using LLM reasoning.

        For cases where embeddings aren't available or additional reasoning is needed.

        Args:
            antigen: Antigen to assess
            description: Text description for LLM analysis

        Returns:
            AnomalyResult from LLM reasoning
        """
        # Placeholder for LLM integration
        # Would call LLM with prompt like:
        # "Given this normal behavior: [...], is this input anomalous: [...]?"

        # For now, return simple heuristic
        is_anomaly = "suspicious" in description.lower() or "malicious" in description.lower()

        return AnomalyResult(
            is_anomaly=is_anomaly,
            anomaly_score=0.8 if is_anomaly else 0.2,
            confidence=0.6,  # Lower confidence for heuristic
            detected_patterns=[],
            explanation=f"Zero-shot detection based on description",
            detector_id=f"{self.agent_name}_zero_shot"
        )

    def get_statistics(self) -> Dict[str, Any]:
        """Get agent statistics"""
        return {
            "agent_name": self.agent_name,
            "num_self_patterns": len(self.self_patterns),
            "num_detectors": len(self.detectors),
            "detection_threshold": self.detection_threshold
        }

    def to_agent_response(self, result: AnomalyResult, success: bool = True,
                         execution_time: float = 0.0) -> AgentResponse:
        """Convert anomaly result to agent response"""
        return AgentResponse(
            agent_name=self.agent_name,
            agent_type="nk_cell",
            result=result.to_dict(),
            success=success,
            execution_time=execution_time,
            metadata=result.metadata
        )

    def save_state(self, path: str) -> None:
        """Persist self patterns, detectors, and embeddings to disk."""
        state = {
            "agent_name": self.agent_name,
            "detection_threshold": self.detection_threshold,
            "num_detectors": self.num_detectors,
            "self_patterns": [pattern.to_dict() for pattern in self.self_patterns],
            "detectors": [
                {
                    "detector_id": detector.detector_id,
                    "threshold": detector.threshold,
                    "creation_time": detector.creation_time,
                }
                for detector in self.detectors
            ],
            "embedding_cache": {
                key: value.tolist() for key, value in self.embedding_cache.items()
            },
        }
        out_path = Path(path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(state, indent=2, ensure_ascii=True), encoding="utf-8")

    @classmethod
    def load_state(cls, path: str) -> "NKCellAgent":
        """Load self patterns, detectors, and embeddings from disk."""
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        agent = cls(
            agent_name=data.get("agent_name", "nk_cell_001"),
            detection_threshold=data.get("detection_threshold", 0.5),
            num_detectors=data.get("num_detectors", 100),
        )
        agent.self_patterns = [Pattern(**pattern) for pattern in data.get("self_patterns", [])]
        agent.detectors = [
            Detector(
                detector_id=detector.get("detector_id", f"{agent.agent_name}_detector"),
                self_patterns=agent.self_patterns,
                threshold=detector.get("threshold", agent.detection_threshold),
                creation_time=detector.get("creation_time"),
            )
            for detector in data.get("detectors", [])
        ]
        agent.embedding_cache = {
            key: np.array(value) for key, value in data.get("embedding_cache", {}).items()
        }
        return agent

```
