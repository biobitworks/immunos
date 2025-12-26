---
title: "immunOS Replication v1: SciFact Baseline and Negative-Selection AIS Framing"
author:
  - name: Byron
    affiliation: BiobitWorks
    email: byron@biobitworks.com
ai_contributors:
  - name: Cursor Auto Agent
    role: Codebase exploration, MCP integration
  - name: Claude Opus 4.5
    role: Documentation, analysis, coordination
software:
  - name: Claude Code
    provider: Anthropic
    role: AI-assisted development environment
  - name: ChatGPT Codex
    provider: OpenAI
    role: Code generation assistance
  - name: VS Code
    provider: Microsoft
    role: Primary IDE
  - name: Obsidian
    provider: Obsidian.md
    role: Knowledge management and notes
  - name: Continue
    provider: Continue.dev
    role: AI-powered code assistant for VS Code
date: 2025-12-26
version: v1.1
type: preprint
status: draft
license: CC-BY-4.0
repository: https://github.com/biobitworks/immunos
---

# immunOS Replication v1: SciFact Baseline and Negative-Selection AIS Framing

**Author**: Byron (BiobitWorks)
**AI Contributors**: Cursor Auto Agent, Claude Opus 4.5
**Software**: Claude Code, ChatGPT Codex, VS Code, Obsidian, Continue
**Date**: December 2025
**Version**: 1.1

## Abstract

We present immunOS, a local-first research verification system inspired by artificial immune system (AIS) principles, specifically the negative selection algorithm for anomaly detection. This preprint documents replication progress across two tracks: (1) a SciFact baseline verification pipeline achieving sentence selection F1 of 0.477 and abstract-level label F1 of 0.510 on open-retrieval dev sets, and (2) alignment with the NegSl-AIS multimodal classification framework. The SciFact baseline provides a foundation for claim verification using retrieval-augmented pipelines, while the negative selection framing motivates detector training against trusted references for downstream hallucination and inconsistency detection. We report oracle and open-retrieval metrics, document reproducibility artifacts, and identify dataset access constraints (MAHNOB-HCI) as blocking dependencies for full NegSl-AIS replication. All code and logs are available locally with paths documented for reproducibility.

## Introduction
immunOS is a local-first research verification system framed as an artificial immune system (AIS) that detects non-self (hallucinations, inconsistencies, and unsupported claims) against a learned self model. The negative selection paradigm motivates detector training against trusted references and subsequent anomaly detection across new inputs. This preprint documents replication progress for the SciFact baseline and the ongoing NegSl-AIS alignment work.

## Related Work (PubMed citations)
- Umair, M. et al. (2025). Negative selection-based artificial immune system (NegSl-AIS)--A hybrid multimodal emotional effect classification model. Results in Engineering 27, 106601. DOI: https://doi.org/10.1016/j.rineng.2025.106601. PubMed: https://pubmed.ncbi.nlm.nih.gov/?term=10.1016%2Fj.rineng.2025.106601
- Ismail, A. et al. (2024). Retrieval augmented scientific claim verification. JAMIA Open. DOI: https://doi.org/10.1093/jamiaopen/ooae021. PubMed: https://pubmed.ncbi.nlm.nih.gov/38455840/
- Pasumarthi, R. et al. (2022). Aggregating pairwise semantic differences for few-shot claim verification. PeerJ Computer Science. DOI: https://doi.org/10.7717/peerj-cs.1137. PubMed: https://pubmed.ncbi.nlm.nih.gov/36426249/
- Li, Y. et al. (2024). Cross-modal credibility modelling for EEG-based multimodal emotion recognition. Journal of Neural Engineering. DOI: https://doi.org/10.1088/1741-2552/ad3987. PubMed: https://pubmed.ncbi.nlm.nih.gov/38565099/
- Ahmad, J. et al. (2018). Artificial Immune System-Negative Selection Classification Algorithm (NSCA) for Four Class Electroencephalogram (EEG) Signals. Frontiers in Human Neuroscience. DOI: https://doi.org/10.3389/fnhum.2018.00439. PubMed: https://pubmed.ncbi.nlm.nih.gov/30524257/

These works ground the claim verification baseline in retrieval-augmented pipelines evaluated on SciFact and few-shot settings, while the AIS literature demonstrates negative selection applied to classification tasks. Recent multimodal emotion recognition studies on MAHNOB-HCI provide the empirical context for NegSl-AIS, motivating why dataset access is a blocking dependency for full replication.

## Methods
- SciFact baseline pipeline: sentence selection followed by sentence labeling, with abstract-level labels computed with and without rationale constraints.
- Environment: local execution with data and logs stored under `/Users/byron/projects/data/immunos_data/research/scifact`.
- Reference notes: `/Users/byron/projects/docs/reference/scifact-baseline-replication.md` and `/Users/byron/projects/docs/reference/negsl-ais-2025-baseline.md`.

### Experimental Setup
- Conda env (x86_64): `/Users/byron/miniconda3-x86/envs/scifact-py38-legacy`
- Key pins: numpy 1.18.2; scipy 1.4.1; scikit-learn 0.22.2.post1; torch 1.5.0; transformers 2.7.0; scispacy 0.2.5; spacy 2.3.9; sentencepiece 0.1.85
- Command refs: pipeline steps in `/Users/byron/projects/docs/reference/scifact-baseline-replication.md`

Figure 1: SciFact pipeline overview (placeholder)

## Results
Dev metrics from `/Users/byron/projects/data/immunos_data/research/scifact/logs/metrics_dev.json`:

| Dev Split | Sentence selection F1 | Sentence label F1 | Abstract label (only) F1 | Abstract label (rationalized) F1 |
| --- | --- | --- | --- | --- |
| Oracle dev | 0.677 | 0.608 | 0.775 | 0.725 |
| Open dev | 0.477 | 0.426 | 0.510 | 0.485 |

Open dev metrics from `/Users/byron/projects/data/immunos_data/research/scifact/logs/metrics_dev_open.json`.

Oracle dev (precision/recall/F1):

| Task | Precision | Recall | F1 |
| --- | --- | --- | --- |
| Sentence selection | 0.794 | 0.590 | 0.677 |
| Sentence label | 0.713 | 0.530 | 0.608 |
| Abstract label (only) | 0.910 | 0.675 | 0.775 |
| Abstract label (rationalized) | 0.852 | 0.632 | 0.725 |

Open dev (precision/recall/F1):

| Task | Precision | Recall | F1 |
| --- | --- | --- | --- |
| Sentence selection | 0.525 | 0.437 | 0.477 |
| Sentence label | 0.469 | 0.391 | 0.426 |
| Abstract label (only) | 0.553 | 0.474 | 0.510 |
| Abstract label (rationalized) | 0.525 | 0.450 | 0.485 |

## Limitations
- MAHNOB-HCI data for NegSl-AIS replication is currently unavailable locally.
- Full NegSl-AIS reproduction (including multimodal feature extraction) is pending.
- PubMed coverage is limited to a targeted subset of related-work queries logged in the publications log.

## Reproducibility
- Metrics log: `/Users/byron/projects/data/immunos_data/research/scifact/logs/metrics_dev.json`
- Dataset root: `/Users/byron/projects/data/immunos_data/research/scifact`
- Replication notes: `/Users/byron/projects/docs/reference/scifact-baseline-replication.md`
- Preprint project: `/Users/byron/projects/immunos-replication-preprint/`

### Reproducibility Checklist
- `metrics_dev.json` and `eval_dev.log` in `/Users/byron/projects/data/immunos_data/research/scifact/logs/`
- Replication log in `/Users/byron/projects/docs/reference/scifact-baseline-replication.md`
- Dataset root in `/Users/byron/projects/data/immunos_data/research/scifact`

## Data/Code Availability
- Data paths and logs are local under `/Users/byron/projects/data/immunos_data`.
- Code and documentation live under `/Users/byron/projects` with readable snapshots under `/Users/byron/projects/docs/code-snapshots/`.

## Citation Relevance Summary

| Citation | Role | Justification |
| --- | --- | --- |
| Umair (2025) NegSl-AIS | core | Foundational negative selection AIS for multimodal classification framing |
| Ismail (2024) CliVER | core | Direct SciFact baseline with PubMed-backed claim verification |
| Pasumarthi (2022) SEED | core | SciFact-focused few-shot claim verification method |
| Zhu (2025) Valsci | supporting | Recent large-batch SciFact verification with LLM workflow |
| Qiu (2024) Multimodal EmoRec | supporting | MAHNOB-HCI benchmark context for NegSl-AIS scope |
| Li (2024) Cross-modal credibility | supporting | MAHNOB-HCI multimodal fusion baselines |
| Tadas (2018) EEG+Face fusion | supporting | Early MAHNOB-HCI multimodal baseline |
| Ahmad (2018) NSCA-EEG | supporting | Negative selection in EEG classification |
| Du (2023) 3D CNN-CRNN | supporting | Recent MAHNOB-HCI results for target performance |
| Zhao (2021) BiLSTM attention | supporting | MAHNOB-HCI multimodal baseline for replication context |
| Ji (2007) NSA review | supporting | Foundational NSA survey for methodological framing |
| Torkamani-Azar (2016) RAIRS2 | adjacent | AIS classifier in biomedical diagnosis (different domain) |
| Cerulo (2013) Negative selection heuristic | adjacent | Negative selection in bioinformatics (different task) |
| Li (2025) Anomaly detection NSA | adjacent | NSA detector generation in intrusion detection |
| Karakose (2013) RL AIS | adjacent | AIS classifier variant for broader methodological context |

---

## References
Sources: `/Users/byron/projects/docs/reference/publications-log.md`. BibTeX: `/Users/byron/projects/research/citations/immunos-references.bib`. Zotero: local library (Zotero app).

- Umair, M. et al. (2025). Negative selection-based artificial immune system (NegSl-AIS)--A hybrid multimodal emotional effect classification model. DOI: https://doi.org/10.1016/j.rineng.2025.106601. PubMed: https://pubmed.ncbi.nlm.nih.gov/?term=10.1016%2Fj.rineng.2025.106601. Summary: Proposes a negative selection AIS for multimodal emotion classification, mapping self/non-self separation to detector training and tuning thresholds for generalization.
- Ismail, A. et al. (2024). Retrieval augmented scientific claim verification. DOI: https://doi.org/10.1093/jamiaopen/ooae021. PubMed: https://pubmed.ncbi.nlm.nih.gov/38455840/. PMCID: TBD. Summary: Presents CliVER, a retrieval-augmented claim verification system using PubMed abstracts and PICO-based rationale extraction, evaluated on CoVERt and SciFact.
- Zhu, J. et al. (2025). Valsci: an open-source, self-hostable literature review utility for automated large-batch scientific claim verification using large language models. DOI: https://doi.org/10.1186/s12859-025-06159-4. PubMed: https://pubmed.ncbi.nlm.nih.gov/40437377/. PMCID: TBD. Summary: Describes a self-hosted, retrieval-augmented LLM workflow for claim verification with bibliometric scoring and lower citation hallucination rates on SciFact.
- Pasumarthi, R. et al. (2022). Aggregating pairwise semantic differences for few-shot claim verification. DOI: https://doi.org/10.7717/peerj-cs.1137. PubMed: https://pubmed.ncbi.nlm.nih.gov/36426249/. PMCID: TBD. Summary: Introduces SEED, aggregating semantic difference vectors to improve few-shot claim verification on FEVER and SciFact.
- Qiu, X. et al. (2024). Multimodal Emotion Recognition Based on Facial Expressions, Speech, and EEG. DOI: https://doi.org/10.1109/OJEMB.2023.3240280. PubMed: https://pubmed.ncbi.nlm.nih.gov/38899017/. PMCID: TBD. Summary: Proposes a three-branch multimodal system with decision-level fusion across facial, speech, and EEG features, including MAHNOB-HCI evaluations.
- Li, Y. et al. (2024). Cross-modal credibility modelling for EEG-based multimodal emotion recognition. DOI: https://doi.org/10.1088/1741-2552/ad3987. PubMed: https://pubmed.ncbi.nlm.nih.gov/38565099/. Summary: Uses intra/inter-modal attention and sequential consistency to improve fusion credibility, reporting gains on DEAP and MAHNOB-HCI.
- Tadas, C. et al. (2018). Emotion Recognition from EEG and Facial Expressions: a Multimodal Approach. DOI: https://doi.org/10.1109/EMBC.2018.8512407. PubMed: https://pubmed.ncbi.nlm.nih.gov/30440451/. Summary: Applies feature-level fusion of EEG and facial microexpressions on MAHNOB-HCI and reports accuracy gains over unimodal baselines.
- Ahmad, J. et al. (2018). Artificial Immune System-Negative Selection Classification Algorithm (NSCA) for Four Class Electroencephalogram (EEG) Signals. DOI: https://doi.org/10.3389/fnhum.2018.00439. PubMed: https://pubmed.ncbi.nlm.nih.gov/30524257/. PMCID: TBD. Summary: Uses MFCC features, stacked auto-encoders, and GA-optimized negative selection detectors to classify motor imagery EEG.
- Torkamani-Azar, R. et al. (2016). RAIRS2 a new expert system for diagnosing tuberculosis with real-world tournament selection mechanism inside artificial immune recognition system. DOI: https://doi.org/10.1007/s11517-015-1323-6. PubMed: https://pubmed.ncbi.nlm.nih.gov/26081904/. Summary: Introduces a hybrid AIRS classifier with tournament selection for TB diagnosis and reports high accuracy on a clinical dataset.
- Cerulo, L. et al. (2013). A negative selection heuristic to predict new transcriptional targets. DOI: https://doi.org/10.1186/1471-2105-14-S1-S3. PubMed: https://pubmed.ncbi.nlm.nih.gov/23368951/. PMCID: TBD. Summary: Selects reliable negative examples from unlabeled data to improve supervised transcriptional target prediction.
- Du, Y. et al. (2023). Attention-based 3D convolutional recurrent neural network model for multimodal emotion recognition. DOI: https://doi.org/10.3389/fnins.2023.1330077. PubMed: https://pubmed.ncbi.nlm.nih.gov/38268710/. PMCID: TBD. Summary: Proposes a 3D CNN-CRNN with attention for EEG and visual fusion, reporting strong accuracy on DEAP and MAHNOB-HCI.
- Zhao, Y. et al. (2021). Expression EEG Multimodal Emotion Recognition Method Based on the Bidirectional LSTM and Attention Mechanism. DOI: https://doi.org/10.1155/2021/9967592. PubMed: https://pubmed.ncbi.nlm.nih.gov/34055043/. PMCID: TBD. Summary: Uses BCN feature extraction and BiLSTM attention fusion for expression+EEG signals, evaluated on MAHNOB-HCI and DEAP.
- Ji, Z. et al. (2007). Revisiting negative selection algorithms. DOI: https://doi.org/10.1162/evco.2007.15.2.223. PubMed: https://pubmed.ncbi.nlm.nih.gov/17535140/. Summary: Reviews negative selection algorithms, detailing representation, matching rules, and hybrid variants in AIS anomaly detection.
- Li, Z. et al. (2025). Generating detectors from anomaly samples via negative selection for network intrusion detection. DOI: https://doi.org/10.1038/s41598-025-20516-6. PubMed: https://pubmed.ncbi.nlm.nih.gov/41107399/. PMCID: TBD. Summary: Generates detectors from anomaly samples to improve coverage in low-dimensional subspaces, outperforming baselines on NSL-KDD and UNSW-NB15.
- Karakose, M. et al. (2013). Reinforcement learning based artificial immune classifier. DOI: https://doi.org/10.1155/2013/581846. PubMed: https://pubmed.ncbi.nlm.nih.gov/23935424/. PMCID: TBD. Summary: Introduces a reinforcement learning AIS classifier and compares performance against negative selection and other AIS variants.

## Zenodo Metadata

| Field | Value |
|-------|-------|
| Title | immunOS Replication v1: SciFact Baseline and Negative-Selection AIS Framing |
| Upload type | Publication (preprint) |
| Publication date | 2025-12-24 |
| Authors | Byron (BiobitWorks) |
| Version | v1.0 |
| Keywords | immunOS, SciFact, claim verification, artificial immune system, negative selection, research integrity |
| License | CC-BY-4.0 |
| Language | en |
| Access rights | Open |
| Related identifiers | GitHub: https://github.com/biobitworks/immunos |
| DOI | Pending Zenodo reservation |

## Ethics/Compliance
- Work is designed for local execution with a local LLM and an air-gapped intent to avoid external data exposure.
- Only public metadata and abstracts are used for literature references.
