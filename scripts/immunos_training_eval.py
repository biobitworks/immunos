#!/usr/bin/env python3
"""
Evaluate IMMUNOS training accuracy on available datasets.

Runs lightweight, offline evaluations for:
- Hallucination (TruthfulQA)
- Research verification (SciFact)
- Network intrusion (NSL-KDD)
"""

from __future__ import annotations

import argparse
import json
import math
import random
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple
from collections import Counter
import statistics

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
MCP_SRC = REPO_ROOT / "immunos-mcp" / "src"
sys.path.insert(0, str(MCP_SRC))

from immunos_mcp.agents.bcell_agent import BCellAgent
from immunos_mcp.agents.nk_cell_agent import NKCellAgent
from immunos_mcp.agents.nk_cell_enhanced import EnhancedNKCellAgent
from immunos_mcp.core.antigen import Antigen
from immunos_mcp.embeddings.simple_text_embedder import SimpleTextEmbedder
from immunos_mcp.embeddings.ollama_embedder import OllamaEmbedder
from immunos_mcp.training.hallucination_trainer import load_truthfulqa
from immunos_mcp.training.research_trainer import load_scifact_claims, apply_retrieved_evidence


@dataclass
class BinaryMetrics:
    accuracy: float
    precision: float
    recall: float
    f1: float
    tp: int
    fp: int
    tn: int
    fn: int


@dataclass
class MultiClassMetrics:
    accuracy: float
    macro_precision: float
    macro_recall: float
    macro_f1: float
    total: int
    correct: int
    per_class_accuracy: Dict[str, float]
    per_class_precision: Dict[str, float]
    per_class_recall: Dict[str, float]
    per_class_f1: Dict[str, float]


def _shuffle(items: List, rng: random.Random) -> List:
    shuffled = items[:]
    rng.shuffle(shuffled)
    return shuffled


def split_list(items: List, train_ratio: float, rng: random.Random) -> Tuple[List, List]:
    if not items:
        return [], []
    shuffled = _shuffle(items, rng)
    split_idx = max(1, int(len(shuffled) * train_ratio))
    return shuffled[:split_idx], shuffled[split_idx:]


def stratified_split(
    items: List,
    labels: List[str],
    train_ratio: float,
    rng: random.Random,
) -> Tuple[List, List, List[str], List[str]]:
    if not items:
        return [], [], [], []
    by_label: Dict[str, List[int]] = {}
    for idx, label in enumerate(labels):
        by_label.setdefault(label, []).append(idx)

    train_indices: List[int] = []
    test_indices: List[int] = []
    for label, indices in by_label.items():
        indices = _shuffle(indices, rng)
        split_idx = max(1, int(len(indices) * train_ratio))
        train_indices.extend(indices[:split_idx])
        test_indices.extend(indices[split_idx:])

    train_items = [items[i] for i in train_indices]
    test_items = [items[i] for i in test_indices]
    train_labels = [labels[i] for i in train_indices]
    test_labels = [labels[i] for i in test_indices]
    return train_items, test_items, train_labels, test_labels


def stratified_kfold_indices(
    labels: List[str],
    num_folds: int,
    rng: random.Random,
) -> List[List[int]]:
    by_label: Dict[str, List[int]] = {}
    for idx, label in enumerate(labels):
        by_label.setdefault(label, []).append(idx)

    folds: List[List[int]] = [[] for _ in range(num_folds)]
    for indices in by_label.values():
        rng.shuffle(indices)
        for i, idx in enumerate(indices):
            folds[i % num_folds].append(idx)

    return folds


def binary_metrics(y_true: Iterable[bool], y_pred: Iterable[bool]) -> BinaryMetrics:
    tp = fp = tn = fn = 0
    for truth, pred in zip(y_true, y_pred):
        if truth and pred:
            tp += 1
        elif truth and not pred:
            fn += 1
        elif not truth and pred:
            fp += 1
        else:
            tn += 1
    total = tp + tn + fp + fn
    accuracy = (tp + tn) / total if total else 0.0
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) else 0.0
    return BinaryMetrics(accuracy, precision, recall, f1, tp, fp, tn, fn)


def multiclass_metrics(y_true: List[str], y_pred: List[str]) -> MultiClassMetrics:
    total = len(y_true)
    correct = sum(1 for truth, pred in zip(y_true, y_pred) if truth == pred)
    accuracy = correct / total if total else 0.0

    labels = sorted(set(y_true) | set(y_pred))
    per_class_accuracy: Dict[str, float] = {}
    per_class_precision: Dict[str, float] = {}
    per_class_recall: Dict[str, float] = {}
    per_class_f1: Dict[str, float] = {}

    for label in labels:
        tp = sum(1 for truth, pred in zip(y_true, y_pred) if truth == label and pred == label)
        fp = sum(1 for truth, pred in zip(y_true, y_pred) if truth != label and pred == label)
        fn = sum(1 for truth, pred in zip(y_true, y_pred) if truth == label and pred != label)
        support = sum(1 for truth in y_true if truth == label)
        per_class_accuracy[label] = tp / support if support else 0.0
        precision = tp / (tp + fp) if (tp + fp) else 0.0
        recall = tp / (tp + fn) if (tp + fn) else 0.0
        f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) else 0.0
        per_class_precision[label] = precision
        per_class_recall[label] = recall
        per_class_f1[label] = f1

    macro_precision = sum(per_class_precision.values()) / len(per_class_precision) if per_class_precision else 0.0
    macro_recall = sum(per_class_recall.values()) / len(per_class_recall) if per_class_recall else 0.0
    macro_f1 = sum(per_class_f1.values()) / len(per_class_f1) if per_class_f1 else 0.0

    return MultiClassMetrics(
        accuracy=accuracy,
        macro_precision=macro_precision,
        macro_recall=macro_recall,
        macro_f1=macro_f1,
        total=total,
        correct=correct,
        per_class_accuracy=per_class_accuracy,
        per_class_precision=per_class_precision,
        per_class_recall=per_class_recall,
        per_class_f1=per_class_f1,
    )


def evaluate_hallucination(
    data_path: Path,
    train_ratio: float,
    max_samples: int,
    nk_threshold: float,
    nk_detectors: int,
    seed: int,
) -> Dict[str, object]:
    truthful, hallucinated = load_truthfulqa(
        path=data_path,
        include_question=True,
        max_samples=max_samples,
    )

    rng = random.Random(seed)
    truthful_train, truthful_test = split_list(truthful, train_ratio, rng)
    halluc_train, halluc_test = split_list(hallucinated, train_ratio, rng)

    embedder = SimpleTextEmbedder()
    bcell = BCellAgent(agent_name="hallucination_bcell_eval")
    nk_cell = NKCellAgent(
        agent_name="hallucination_nk_eval",
        detection_threshold=nk_threshold,
        num_detectors=nk_detectors,
    )

    train_antigens = truthful_train + halluc_train
    train_embeddings = [embedder.embed(a.get_text_content()) for a in train_antigens]
    bcell.train(train_antigens, embeddings=train_embeddings)

    truthful_embeddings = [embedder.embed(a.get_text_content()) for a in truthful_train]
    nk_cell.train_on_self(truthful_train, embeddings=truthful_embeddings)

    test_antigens = truthful_test + halluc_test
    test_labels = ["truthful"] * len(truthful_test) + ["hallucinated"] * len(halluc_test)
    test_embeddings = [embedder.embed(a.get_text_content()) for a in test_antigens]

    bcell_preds: List[str] = []
    for antigen, embedding in zip(test_antigens, test_embeddings):
        result = bcell.recognize(antigen, antigen_embedding=embedding, strategy="sha")
        bcell_preds.append(result.predicted_class or "unknown")

    bcell_metrics = multiclass_metrics(test_labels, bcell_preds)

    nk_truth = [label == "hallucinated" for label in test_labels]
    nk_preds = []
    for antigen, embedding in zip(test_antigens, test_embeddings):
        result = nk_cell.detect_novelty(antigen, antigen_embedding=embedding)
        nk_preds.append(result.is_anomaly)
    nk_metrics = binary_metrics(nk_truth, nk_preds)

    return {
        "dataset": str(data_path),
        "train_samples": len(train_antigens),
        "test_samples": len(test_antigens),
        "bcell_metrics": bcell_metrics,
        "nk_metrics": nk_metrics,
    }


def balance_pairs(
    antigens: List[Antigen],
    labels: List[str],
    rng: random.Random,
) -> Tuple[List[Antigen], List[str]]:
    by_label: Dict[str, List[int]] = {}
    for idx, label in enumerate(labels):
        by_label.setdefault(label, []).append(idx)
    min_count = min(len(items) for items in by_label.values()) if by_label else 0
    balanced_indices: List[int] = []
    for label, indices in by_label.items():
        rng.shuffle(indices)
        balanced_indices.extend(indices[:min_count])
    rng.shuffle(balanced_indices)
    return [antigens[i] for i in balanced_indices], [labels[i] for i in balanced_indices]


def build_embedder(name: str, ollama_model: Optional[str], ollama_base_url: Optional[str]):
    if name == "ollama":
        return OllamaEmbedder(
            model=ollama_model or "nomic-embed-text",
            base_url=ollama_base_url
        )
    return SimpleTextEmbedder()


def evaluate_research_split(
    train_antigens: List[Antigen],
    train_labels: List[str],
    test_antigens: List[Antigen],
    test_labels: List[str],
    embedder: object,
    nk_threshold: float,
    nk_detectors: int,
    nk_mode: str,
    nk_threshold_method: str,
    nk_detectors_per_class: int,
    bcell_affinity: str,
    bcell_embedding_weight: float,
    bcell_strategy: str,
    split_label: str,
    include_unknown: bool,
    balance_classes: bool,
    seed: int,
    embedder_name: str,
    ollama_model: Optional[str],
    include_evidence: bool,
    evidence_max: int,
) -> Dict[str, object]:
    rng = random.Random(seed)

    if not include_unknown:
        filtered = [(a, l) for a, l in zip(train_antigens, train_labels) if l != "unknown"]
        train_antigens = [item[0] for item in filtered]
        train_labels = [item[1] for item in filtered]

    if balance_classes and train_antigens:
        train_antigens, train_labels = balance_pairs(train_antigens, train_labels, rng)

    bcell = BCellAgent(
        agent_name="research_bcell_eval",
        affinity_method=bcell_affinity,
        embedding_weight=bcell_embedding_weight,
    )
    if nk_mode == "enhanced":
        nk_cell = EnhancedNKCellAgent(
            agent_name="research_nk_eval_enhanced",
            threshold_method=nk_threshold_method,
            detectors_per_class=nk_detectors_per_class,
        )
    else:
        nk_cell = NKCellAgent(
            agent_name="research_nk_eval",
            detection_threshold=nk_threshold,
            num_detectors=nk_detectors,
        )

    train_embeddings = [embedder.embed(a.get_text_content()) for a in train_antigens]
    bcell.train(train_antigens, embeddings=train_embeddings)

    support_pairs = [(a, l) for a, l in zip(train_antigens, train_labels) if l == "support"]
    support_antigens = [item[0] for item in support_pairs]
    support_embeddings = [embedder.embed(a.get_text_content()) for a in support_antigens]
    nk_cell.train_on_self(support_antigens, embeddings=support_embeddings)

    test_embeddings = [embedder.embed(a.get_text_content()) for a in test_antigens]

    bcell_preds: List[str] = []
    for antigen, embedding in zip(test_antigens, test_embeddings):
        result = bcell.recognize(antigen, antigen_embedding=embedding, strategy=bcell_strategy)
        bcell_preds.append(result.predicted_class or "unknown")
    bcell_metrics = multiclass_metrics(test_labels, bcell_preds)

    nk_truth = [label != "support" for label in test_labels]
    nk_preds = []
    for antigen, embedding in zip(test_antigens, test_embeddings):
        result = nk_cell.detect_novelty(antigen, antigen_embedding=embedding)
        nk_preds.append(result.is_anomaly)
    nk_metrics = binary_metrics(nk_truth, nk_preds)

    return {
        "split": split_label,
        "balanced": balance_classes,
        "bcell_affinity": bcell_affinity,
        "bcell_embedding_weight": bcell_embedding_weight,
        "bcell_strategy": bcell_strategy,
        "embedder": embedder_name,
        "embedder_model": ollama_model if embedder_name == "ollama" else "simple_text",
        "evidence": include_evidence,
        "evidence_max": evidence_max if include_evidence else 0,
        "nk_mode": nk_mode,
        "nk_threshold_method": nk_threshold_method if nk_mode == "enhanced" else None,
        "nk_detectors_per_class": nk_detectors_per_class if nk_mode == "enhanced" else None,
        "train_label_counts": dict(Counter(train_labels)),
        "test_label_counts": dict(Counter(test_labels)),
        "train_samples": len(train_antigens),
        "test_samples": len(test_antigens),
        "bcell_metrics": bcell_metrics,
        "nk_metrics": nk_metrics,
    }


def evaluate_research(
    train_path: Optional[Path],
    dev_path: Optional[Path],
    train_ratio: float,
    max_samples: int,
    nk_threshold: float,
    nk_detectors: int,
    nk_mode: str,
    nk_threshold_method: str,
    nk_detectors_per_class: int,
    bcell_affinity: str,
    bcell_embedding_weight: float,
    bcell_strategy: str,
    split_label_override: Optional[str],
    seed: int,
    include_unknown: bool,
    balance_classes: bool,
    force_random_split: bool,
    embedder_name: str,
    ollama_model: Optional[str],
    ollama_base_url: Optional[str],
    include_evidence: bool,
    evidence_mode: str,
    corpus_path: Optional[Path],
    evidence_max: int,
    evidence_candidates_per_doc: int,
) -> Dict[str, object]:
    rng = random.Random(seed)

    if train_path and dev_path and train_path.exists() and dev_path.exists() and not force_random_split:
        train_antigens, train_labels = load_scifact_claims(
            train_path,
            max_samples=0,
            include_evidence=include_evidence,
            corpus_path=corpus_path,
            max_evidence_sentences=evidence_max,
        )
        test_antigens, test_labels = load_scifact_claims(
            dev_path,
            max_samples=0,
            include_evidence=include_evidence,
            corpus_path=corpus_path,
            max_evidence_sentences=evidence_max,
        )
        split_label = split_label_override or "train/dev"
        if max_samples:
            train_antigens = train_antigens[:max_samples]
            train_labels = train_labels[:max_samples]
    else:
        antigens, labels = load_scifact_claims(
            dev_path or train_path,
            max_samples=max_samples,
            include_evidence=include_evidence,
            corpus_path=corpus_path,
            max_evidence_sentences=evidence_max,
        )
        train_antigens, test_antigens, train_labels, test_labels = stratified_split(
            antigens, labels, train_ratio, rng
        )
        split_label = split_label_override or "random"

    embedder = build_embedder(embedder_name, ollama_model, ollama_base_url)
    if evidence_mode == "retrieval":
        train_antigens = apply_retrieved_evidence(
            train_antigens,
            embedder,
            corpus_path=corpus_path,
            max_evidence_sentences=evidence_max,
            max_candidates_per_doc=evidence_candidates_per_doc,
        )
        test_antigens = apply_retrieved_evidence(
            test_antigens,
            embedder,
            corpus_path=corpus_path,
            max_evidence_sentences=evidence_max,
            max_candidates_per_doc=evidence_candidates_per_doc,
        )
    result = evaluate_research_split(
        train_antigens=train_antigens,
        train_labels=train_labels,
        test_antigens=test_antigens,
        test_labels=test_labels,
        embedder=embedder,
        nk_threshold=nk_threshold,
        nk_detectors=nk_detectors,
        nk_mode=nk_mode,
        nk_threshold_method=nk_threshold_method,
        nk_detectors_per_class=nk_detectors_per_class,
        bcell_affinity=bcell_affinity,
        bcell_embedding_weight=bcell_embedding_weight,
        bcell_strategy=bcell_strategy,
        split_label=split_label,
        include_unknown=include_unknown,
        balance_classes=balance_classes,
        seed=seed,
        embedder_name=embedder_name,
        ollama_model=ollama_model,
        include_evidence=include_evidence,
        evidence_max=evidence_max,
    )
    result["dataset"] = str(train_path or dev_path)
    result["evidence_mode"] = evidence_mode
    if evidence_mode != "none":
        result["evidence_max"] = evidence_max
    return result


def summarize_fold_metrics(values: List[float]) -> Dict[str, float]:
    if not values:
        return {"mean": 0.0, "std": 0.0}
    if len(values) == 1:
        return {"mean": values[0], "std": 0.0}
    return {
        "mean": statistics.mean(values),
        "std": statistics.pstdev(values),
    }


def evaluate_research_kfold(
    train_path: Path,
    dev_path: Path,
    num_folds: int,
    max_samples: int,
    nk_threshold: float,
    nk_detectors: int,
    nk_mode: str,
    nk_threshold_method: str,
    nk_detectors_per_class: int,
    bcell_affinity: str,
    bcell_embedding_weight: float,
    bcell_strategy: str,
    seed: int,
    include_unknown: bool,
    balance_classes: bool,
    embedder_name: str,
    ollama_model: Optional[str],
    ollama_base_url: Optional[str],
    include_evidence: bool,
    evidence_mode: str,
    corpus_path: Optional[Path],
    evidence_max: int,
    evidence_candidates_per_doc: int,
) -> Dict[str, object]:
    rng = random.Random(seed)

    train_antigens, train_labels = load_scifact_claims(
        train_path,
        max_samples=0,
        include_evidence=include_evidence,
        corpus_path=corpus_path,
        max_evidence_sentences=evidence_max,
    )
    dev_antigens, dev_labels = load_scifact_claims(
        dev_path,
        max_samples=0,
        include_evidence=include_evidence,
        corpus_path=corpus_path,
        max_evidence_sentences=evidence_max,
    )

    antigens = train_antigens + dev_antigens
    labels = train_labels + dev_labels
    if max_samples:
        antigens = antigens[:max_samples]
        labels = labels[:max_samples]

    folds = stratified_kfold_indices(labels, num_folds, rng)
    embedder = build_embedder(embedder_name, ollama_model, ollama_base_url)
    if evidence_mode == "retrieval":
        antigens = apply_retrieved_evidence(
            antigens,
            embedder,
            corpus_path=corpus_path,
            max_evidence_sentences=evidence_max,
            max_candidates_per_doc=evidence_candidates_per_doc,
        )

    fold_results = []
    for fold_idx, test_indices in enumerate(folds):
        test_set = set(test_indices)
        train_indices = [idx for idx in range(len(antigens)) if idx not in test_set]
        fold_train_antigens = [antigens[idx] for idx in train_indices]
        fold_train_labels = [labels[idx] for idx in train_indices]
        fold_test_antigens = [antigens[idx] for idx in test_indices]
        fold_test_labels = [labels[idx] for idx in test_indices]

        fold_result = evaluate_research_split(
            train_antigens=fold_train_antigens,
            train_labels=fold_train_labels,
            test_antigens=fold_test_antigens,
            test_labels=fold_test_labels,
            embedder=embedder,
            nk_threshold=nk_threshold,
            nk_detectors=nk_detectors,
            nk_mode=nk_mode,
            nk_threshold_method=nk_threshold_method,
            nk_detectors_per_class=nk_detectors_per_class,
            bcell_affinity=bcell_affinity,
            bcell_embedding_weight=bcell_embedding_weight,
            bcell_strategy=bcell_strategy,
            split_label=f"kfold-{fold_idx + 1}/{num_folds}",
            include_unknown=include_unknown,
            balance_classes=balance_classes,
            seed=seed,
            embedder_name=embedder_name,
            ollama_model=ollama_model,
            include_evidence=include_evidence,
            evidence_max=evidence_max,
        )
        fold_results.append(fold_result)

    bcell_acc = [fold["bcell_metrics"].accuracy for fold in fold_results]
    bcell_f1 = [fold["bcell_metrics"].macro_f1 for fold in fold_results]
    bcell_prec = [fold["bcell_metrics"].macro_precision for fold in fold_results]
    bcell_rec = [fold["bcell_metrics"].macro_recall for fold in fold_results]
    nk_acc = [fold["nk_metrics"].accuracy for fold in fold_results]
    nk_prec = [fold["nk_metrics"].precision for fold in fold_results]
    nk_rec = [fold["nk_metrics"].recall for fold in fold_results]
    nk_f1 = [fold["nk_metrics"].f1 for fold in fold_results]

    return {
        "split": f"kfold-{num_folds}",
        "balanced": balance_classes,
        "bcell_affinity": bcell_affinity,
        "bcell_embedding_weight": bcell_embedding_weight,
        "bcell_strategy": bcell_strategy,
        "embedder": embedder_name,
        "embedder_model": ollama_model if embedder_name == "ollama" else "simple_text",
        "evidence": include_evidence,
        "evidence_max": evidence_max if include_evidence else 0,
        "evidence_mode": evidence_mode,
        "nk_mode": nk_mode,
        "nk_threshold_method": nk_threshold_method if nk_mode == "enhanced" else None,
        "nk_detectors_per_class": nk_detectors_per_class if nk_mode == "enhanced" else None,
        "folds": fold_results,
        "bcell_accuracy": summarize_fold_metrics(bcell_acc),
        "bcell_macro_f1": summarize_fold_metrics(bcell_f1),
        "bcell_macro_precision": summarize_fold_metrics(bcell_prec),
        "bcell_macro_recall": summarize_fold_metrics(bcell_rec),
        "nk_accuracy": summarize_fold_metrics(nk_acc),
        "nk_precision": summarize_fold_metrics(nk_prec),
        "nk_recall": summarize_fold_metrics(nk_rec),
        "nk_f1": summarize_fold_metrics(nk_f1),
    }


KDD_CATEGORICAL_INDEXES = {1, 2, 3}


def _encode_category(value: str, mapping: Dict[str, int]) -> int:
    if value not in mapping:
        mapping[value] = len(mapping) + 1
    return mapping[value]


def load_kdd(
    path: Path,
    max_samples: int,
    category_maps: Optional[Dict[int, Dict[str, int]]] = None,
) -> Tuple[List[List[float]], List[str], Dict[int, Dict[str, int]]]:
    if category_maps is None:
        category_maps = {idx: {} for idx in KDD_CATEGORICAL_INDEXES}

    features: List[List[float]] = []
    labels: List[str] = []

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            raw = line.strip()
            if not raw:
                continue
            parts = raw.split(",")
            if len(parts) < 42:
                continue
            label = parts[-2].strip()

            row: List[float] = []
            for i, value in enumerate(parts[:-2]):
                if i in KDD_CATEGORICAL_INDEXES:
                    encoded = _encode_category(value, category_maps[i])
                    row.append(float(encoded))
                else:
                    try:
                        row.append(float(value))
                    except ValueError:
                        row.append(0.0)

            features.append(row)
            labels.append(label)

            if max_samples and len(features) >= max_samples:
                break

    return features, labels, category_maps


def normalize_with_stats(data: np.ndarray, mins: np.ndarray, maxs: np.ndarray) -> np.ndarray:
    ranges = np.where(maxs - mins == 0, 1.0, maxs - mins)
    return (data - mins) / ranges


def evaluate_network(
    train_path: Path,
    test_path: Path,
    train_max_samples: int,
    test_max_samples: int,
    seed: int,
) -> Dict[str, object]:
    rng = random.Random(seed)

    train_features_all, train_labels_all, category_maps = load_kdd(
        train_path,
        max_samples=0,
        category_maps=None,
    )

    normal_indices = [i for i, label in enumerate(train_labels_all) if label == "normal"]
    rng.shuffle(normal_indices)
    if train_max_samples:
        normal_indices = normal_indices[:train_max_samples]

    train_features = [train_features_all[i] for i in normal_indices]
    train_antigens = [Antigen.from_dict({"features": row}, class_label="normal") for row in train_features]

    train_array = np.array(train_features, dtype=np.float32)
    mins = train_array.min(axis=0)
    maxs = train_array.max(axis=0)
    train_embeddings = normalize_with_stats(train_array, mins, maxs)

    nk_cell = EnhancedNKCellAgent(agent_name="network_nk_eval")
    nk_cell.train_on_self(train_antigens, embeddings=list(train_embeddings))

    test_features, test_labels, _ = load_kdd(
        test_path,
        max_samples=test_max_samples,
        category_maps=category_maps,
    )
    test_array = np.array(test_features, dtype=np.float32)
    test_embeddings = normalize_with_stats(test_array, mins, maxs)
    test_antigens = [Antigen.from_dict({"features": row}, class_label="unknown")
                     for row in test_features]

    nk_truth = [label != "normal" for label in test_labels]
    nk_preds = []
    for antigen, embedding in zip(test_antigens, test_embeddings):
        result = nk_cell.detect_novelty(antigen, antigen_embedding=embedding)
        nk_preds.append(result.is_anomaly)
    nk_metrics = binary_metrics(nk_truth, nk_preds)

    return {
        "train_samples": len(train_antigens),
        "test_samples": len(test_antigens),
        "nk_metrics": nk_metrics,
    }


def metrics_to_dict(metrics: object) -> Dict[str, object]:
    if isinstance(metrics, BinaryMetrics):
        return {
            "accuracy": metrics.accuracy,
            "precision": metrics.precision,
            "recall": metrics.recall,
            "f1": metrics.f1,
            "tp": metrics.tp,
            "fp": metrics.fp,
            "tn": metrics.tn,
            "fn": metrics.fn,
        }
    if isinstance(metrics, MultiClassMetrics):
        return {
            "accuracy": metrics.accuracy,
            "macro_precision": metrics.macro_precision,
            "macro_recall": metrics.macro_recall,
            "macro_f1": metrics.macro_f1,
            "total": metrics.total,
            "correct": metrics.correct,
            "per_class_accuracy": metrics.per_class_accuracy,
            "per_class_precision": metrics.per_class_precision,
            "per_class_recall": metrics.per_class_recall,
            "per_class_f1": metrics.per_class_f1,
        }
    return {}


def format_percent(value: float) -> str:
    return f"{value * 100:.1f}%"


def build_report(results: Dict[str, object]) -> str:
    lines = []
    lines.append(f"# IMMUNOS Training Evaluation ({datetime.now().strftime('%Y-%m-%d %H:%M')})")
    lines.append("")
    lines.append("## Summary")
    lines.append("")

    if "hallucination" in results:
        h = results["hallucination"]
        bcell = h["bcell_metrics"]
        nk = h["nk_metrics"]
        lines.append("- Hallucination (TruthfulQA)")
        lines.append(f"  - Train/Test: {h['train_samples']} / {h['test_samples']}")
        lines.append(f"  - B Cell accuracy: {format_percent(bcell.accuracy)}")
        lines.append(f"  - NK accuracy: {format_percent(nk.accuracy)} (precision {format_percent(nk.precision)}, recall {format_percent(nk.recall)})")
        lines.append("")

    for key in ("research_dev", "research_test"):
        if key in results:
            r = results[key]
            bcell = r["bcell_metrics"]
            nk = r["nk_metrics"]
            label = "Research (SciFact)"
            split_name = r.get("split", "unknown")
            if key == "research_test":
                label = "Research (SciFact Test)"
            lines.append(f"- {label}")
            lines.append(f"  - Split: {split_name} (balanced={r.get('balanced', False)})")
            lines.append(f"  - B Cell: affinity={r.get('bcell_affinity')}, "
                         f"strategy={r.get('bcell_strategy')}, "
                         f"embed_weight={r.get('bcell_embedding_weight')}")
            lines.append(f"  - Embedder: {r.get('embedder', 'simple_text')} ({r.get('embedder_model', 'simple_text')})")
            evidence_mode = r.get("evidence_mode") or ("gold" if r.get("evidence") else "none")
            lines.append(f"  - Evidence: {evidence_mode} (max={r.get('evidence_max', 0)})")
            if r.get("nk_mode"):
                if r.get("nk_mode") == "enhanced":
                    lines.append("  - NK mode: enhanced "
                                 f"(threshold={r.get('nk_threshold_method')}, "
                                 f"detectors/class={r.get('nk_detectors_per_class')})")
                else:
                    lines.append("  - NK mode: classic")
            if r.get("train_label_counts") is not None:
                lines.append(f"  - Train labels: {r.get('train_label_counts')}")
            if r.get("test_label_counts") is not None:
                lines.append(f"  - Test labels: {r.get('test_label_counts')}")
                if len(r.get("test_label_counts", {})) <= 1:
                    lines.append("  - Warning: test split has a single label; "
                                 "test metrics are not publication-grade.")
            lines.append(f"  - Train/Test: {r['train_samples']} / {r['test_samples']}")
            lines.append(f"  - B Cell accuracy: {format_percent(bcell.accuracy)} "
                         f"(macro-F1 {format_percent(bcell.macro_f1)})")
            lines.append(f"  - NK accuracy: {format_percent(nk.accuracy)} "
                         f"(precision {format_percent(nk.precision)}, recall {format_percent(nk.recall)})")
            lines.append("")

    if "research_kfold" in results:
        r = results["research_kfold"]
        lines.append("- Research (SciFact K-Fold)")
        lines.append(f"  - Split: {r.get('split')} (balanced={r.get('balanced', False)})")
        lines.append(f"  - B Cell: affinity={r.get('bcell_affinity')}, "
                     f"strategy={r.get('bcell_strategy')}, "
                     f"embed_weight={r.get('bcell_embedding_weight')}")
        lines.append(f"  - Embedder: {r.get('embedder', 'simple_text')} ({r.get('embedder_model', 'simple_text')})")
        evidence_mode = r.get("evidence_mode") or ("gold" if r.get("evidence") else "none")
        lines.append(f"  - Evidence: {evidence_mode} (max={r.get('evidence_max', 0)})")
        if r.get("nk_mode"):
            if r.get("nk_mode") == "enhanced":
                lines.append("  - NK mode: enhanced "
                             f"(threshold={r.get('nk_threshold_method')}, "
                             f"detectors/class={r.get('nk_detectors_per_class')})")
            else:
                lines.append("  - NK mode: classic")
        bcell_acc = r["bcell_accuracy"]
        bcell_f1 = r["bcell_macro_f1"]
        nk_acc = r["nk_accuracy"]
        nk_f1 = r["nk_f1"]
        lines.append(f"  - B Cell accuracy: {format_percent(bcell_acc['mean'])} ± {format_percent(bcell_acc['std'])}")
        lines.append(f"  - B Cell macro-F1: {format_percent(bcell_f1['mean'])} ± {format_percent(bcell_f1['std'])}")
        lines.append(f"  - NK accuracy: {format_percent(nk_acc['mean'])} ± {format_percent(nk_acc['std'])}")
        lines.append(f"  - NK F1: {format_percent(nk_f1['mean'])} ± {format_percent(nk_f1['std'])}")
        lines.append("")

    if "network" in results:
        n = results["network"]
        nk = n["nk_metrics"]
        lines.append("- Network (NSL-KDD)")
        lines.append(f"  - Train/Test: {n['train_samples']} / {n['test_samples']}")
        lines.append(f"  - NK accuracy: {format_percent(nk.accuracy)} (precision {format_percent(nk.precision)}, recall {format_percent(nk.recall)})")
        lines.append("")

    lines.append("## Detailed Metrics")
    lines.append("")
    for key, entry in results.items():
        lines.append(f"### {key.title()}")
        for metric_key in ("bcell_metrics", "nk_metrics"):
            if metric_key in entry:
                metrics = entry[metric_key]
                lines.append(f"- {metric_key}: `{json.dumps(metrics_to_dict(metrics), indent=2)}`")
        lines.append("")

    return "\n".join(lines).strip() + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description="Evaluate IMMUNOS training accuracy")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--train-ratio", type=float, default=0.8, help="Train/test split ratio for text tasks")
    parser.add_argument("--hallucination-max", type=int, default=500, help="Max TruthfulQA samples (0 = all)")
    parser.add_argument("--research-max", type=int, default=1000, help="Max SciFact claims (0 = all)")
    parser.add_argument("--network-train-max", type=int, default=5000, help="Max normal samples for KDD training (0 = all)")
    parser.add_argument("--network-test-max", type=int, default=2000, help="Max KDD test rows (0 = all)")
    parser.add_argument("--include-unknown", action="store_true", help="Include unknown claims in research B cell training")
    parser.add_argument("--research-balance", action="store_true", help="Balance research labels for training")
    parser.add_argument("--research-random-split", action="store_true", help="Use random split instead of train/dev")
    parser.add_argument("--research-include-evidence", action="store_true",
                        help="Append SciFact evidence sentences to claim text")
    parser.add_argument("--research-evidence-mode", type=str, default=None,
                        choices=["gold", "retrieval", "none"],
                        help="Evidence augmentation mode (default: gold if --research-include-evidence)")
    parser.add_argument("--research-evidence-max", type=int, default=3,
                        help="Max evidence sentences appended per claim")
    parser.add_argument("--research-corpus", type=str, default=None,
                        help="Path to SciFact corpus.jsonl (defaults to data root)")
    parser.add_argument("--research-evidence-candidates-per-doc", type=int, default=8,
                        help="Max sentences per cited document to consider for retrieval")
    parser.add_argument("--research-embedder", type=str, default="simple_text",
                        choices=["simple_text", "ollama"],
                        help="Embedding backend for research evaluation")
    parser.add_argument("--research-bcell-affinity", type=str, default="hybrid",
                        choices=["traditional", "embedding", "hybrid"],
                        help="B Cell affinity method for research evaluation")
    parser.add_argument("--research-bcell-embedding-weight", type=float, default=0.7,
                        help="Embedding weight for hybrid affinity (0-1)")
    parser.add_argument("--research-bcell-strategy", type=str, default="sha",
                        choices=["sha", "rha"],
                        help="B Cell recognition strategy")
    parser.add_argument("--ollama-model", type=str, default="nomic-embed-text",
                        help="Ollama embedding model name")
    parser.add_argument("--ollama-base-url", type=str, default=None,
                        help="Override Ollama base URL (default: http://localhost:11434)")
    parser.add_argument("--research-nk-mode", type=str, default="classic",
                        choices=["classic", "enhanced"],
                        help="NK detector implementation for research evaluation")
    parser.add_argument("--research-nk-threshold", type=float, default=0.85,
                        help="NK detection threshold for research evaluation")
    parser.add_argument("--research-nk-detectors", type=int, default=100,
                        help="Number of NK detectors for research evaluation")
    parser.add_argument("--research-nk-threshold-method", type=str, default="min_distance",
                        choices=["min_distance", "mean_distance", "percentile"],
                        help="Enhanced NK threshold method")
    parser.add_argument("--research-nk-detectors-per-class", type=int, default=25,
                        help="Enhanced NK detectors per class")
    parser.add_argument("--research-eval-test", action="store_true",
                        help="Evaluate research agents on the SciFact test split")
    parser.add_argument("--research-kfolds", type=int, default=0,
                        help="Run stratified k-fold evaluation on train+dev (0 = off)")
    parser.add_argument("--skip-hallucination", action="store_true", help="Skip hallucination evaluation")
    parser.add_argument("--skip-network", action="store_true", help="Skip network evaluation")
    parser.add_argument("--output", type=str, default=None, help="Optional markdown output path")
    args = parser.parse_args()

    results: Dict[str, object] = {}

    hallu_path = REPO_ROOT / "data" / "immunos_data" / "hallucination" / "truthfulqa" / "TruthfulQA.csv"
    if not args.skip_hallucination:
        if hallu_path.exists():
            results["hallucination"] = evaluate_hallucination(
                hallu_path,
                train_ratio=args.train_ratio,
                max_samples=args.hallucination_max,
                nk_threshold=0.8,
                nk_detectors=100,
                seed=args.seed,
            )
        else:
            print(f"[WARN] Hallucination dataset not found at {hallu_path}")

    research_train_path = REPO_ROOT / "data" / "immunos_data" / "research" / "scifact" / "data" / "claims_train.jsonl"
    research_dev_path = REPO_ROOT / "data" / "immunos_data" / "research" / "scifact" / "data" / "claims_dev.jsonl"
    research_test_path = REPO_ROOT / "data" / "immunos_data" / "research" / "scifact" / "data" / "claims_test.jsonl"
    research_corpus_path = Path(args.research_corpus) if args.research_corpus else (
        REPO_ROOT / "data" / "immunos_data" / "research" / "scifact" / "data" / "corpus.jsonl"
    )
    evidence_mode = args.research_evidence_mode or ("gold" if args.research_include_evidence else "none")
    include_gold_evidence = evidence_mode == "gold"
    if (research_train_path.exists() or research_dev_path.exists()):
        results["research_dev"] = evaluate_research(
            research_train_path if research_train_path.exists() else None,
            research_dev_path if research_dev_path.exists() else None,
            train_ratio=args.train_ratio,
            max_samples=args.research_max,
            nk_threshold=args.research_nk_threshold,
            nk_detectors=args.research_nk_detectors,
            nk_mode=args.research_nk_mode,
            nk_threshold_method=args.research_nk_threshold_method,
            nk_detectors_per_class=args.research_nk_detectors_per_class,
            bcell_affinity=args.research_bcell_affinity,
            bcell_embedding_weight=args.research_bcell_embedding_weight,
            bcell_strategy=args.research_bcell_strategy,
            split_label_override="train/dev",
            seed=args.seed,
            include_unknown=args.include_unknown,
            balance_classes=args.research_balance,
            force_random_split=args.research_random_split,
            embedder_name=args.research_embedder,
            ollama_model=args.ollama_model,
            ollama_base_url=args.ollama_base_url,
            include_evidence=include_gold_evidence,
            evidence_mode=evidence_mode,
            corpus_path=research_corpus_path,
            evidence_max=args.research_evidence_max,
            evidence_candidates_per_doc=args.research_evidence_candidates_per_doc,
        )
        if args.research_kfolds and args.research_kfolds >= 2:
            if research_train_path.exists() and research_dev_path.exists():
                results["research_kfold"] = evaluate_research_kfold(
                    research_train_path,
                    research_dev_path,
                    num_folds=args.research_kfolds,
                    max_samples=args.research_max,
                    nk_threshold=args.research_nk_threshold,
                    nk_detectors=args.research_nk_detectors,
                    nk_mode=args.research_nk_mode,
                    nk_threshold_method=args.research_nk_threshold_method,
                    nk_detectors_per_class=args.research_nk_detectors_per_class,
                    bcell_affinity=args.research_bcell_affinity,
                    bcell_embedding_weight=args.research_bcell_embedding_weight,
                    bcell_strategy=args.research_bcell_strategy,
                    seed=args.seed,
                    include_unknown=args.include_unknown,
                    balance_classes=args.research_balance,
                    embedder_name=args.research_embedder,
                    ollama_model=args.ollama_model,
                    ollama_base_url=args.ollama_base_url,
                    include_evidence=include_gold_evidence,
                    evidence_mode=evidence_mode,
                    corpus_path=research_corpus_path,
                    evidence_max=args.research_evidence_max,
                    evidence_candidates_per_doc=args.research_evidence_candidates_per_doc,
                )
            else:
                print("[WARN] K-fold evaluation requires both train and dev splits.")
        if args.research_eval_test and research_test_path.exists():
            results["research_test"] = evaluate_research(
                research_train_path if research_train_path.exists() else None,
                research_test_path,
                train_ratio=args.train_ratio,
                max_samples=args.research_max,
                nk_threshold=args.research_nk_threshold,
                nk_detectors=args.research_nk_detectors,
                nk_mode=args.research_nk_mode,
                nk_threshold_method=args.research_nk_threshold_method,
                nk_detectors_per_class=args.research_nk_detectors_per_class,
                bcell_affinity=args.research_bcell_affinity,
                bcell_embedding_weight=args.research_bcell_embedding_weight,
                bcell_strategy=args.research_bcell_strategy,
                split_label_override="train/test",
                seed=args.seed,
                include_unknown=args.include_unknown,
                balance_classes=args.research_balance,
                force_random_split=args.research_random_split,
                embedder_name=args.research_embedder,
                ollama_model=args.ollama_model,
                ollama_base_url=args.ollama_base_url,
                include_evidence=include_gold_evidence,
                evidence_mode=evidence_mode,
                corpus_path=research_corpus_path,
                evidence_max=args.research_evidence_max,
                evidence_candidates_per_doc=args.research_evidence_candidates_per_doc,
            )
    else:
        print(f"[WARN] Research dataset not found at {research_path}")

    train_path = REPO_ROOT / "data" / "immunos_data" / "network" / "KDDTrain+.txt"
    test_path = REPO_ROOT / "data" / "immunos_data" / "network" / "KDDTest+.txt"
    if not args.skip_network:
        if train_path.exists() and test_path.exists():
            results["network"] = evaluate_network(
                train_path,
                test_path,
                train_max_samples=args.network_train_max,
                test_max_samples=args.network_test_max,
                seed=args.seed,
            )
        else:
            print("[WARN] Network dataset not found for evaluation")

    if not results:
        print("No datasets available for evaluation.")
        return

    report = build_report(results)
    if args.output:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(report, encoding="utf-8")
        print(f"Saved report to {output_path}")
    else:
        print(report)


if __name__ == "__main__":
    main()
