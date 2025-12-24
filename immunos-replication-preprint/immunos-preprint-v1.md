# immunOS Replication v1: SciFact Baseline and Negative-Selection AIS Framing

## Abstract
We report initial replication progress for immunOS, focusing on a SciFact baseline verification pipeline and an AIS-inspired negative selection framing. The SciFact replication currently provides dev-set metrics for sentence selection, sentence labeling, and abstract-level labeling, while the NegSl-AIS replication remains in progress due to data availability constraints. This draft captures the current experimental status, logs, and reproducibility hooks suitable for a bioRxiv/Zenodo preprint update.

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

## Results
Dev metrics from `/Users/byron/projects/data/immunos_data/research/scifact/logs/metrics_dev.json`:

| Task | Precision | Recall | F1 |
| --- | --- | --- | --- |
| Sentence selection | 0.794 | 0.590 | 0.677 |
| Sentence label | 0.713 | 0.530 | 0.608 |
| Abstract label (only) | 0.910 | 0.675 | 0.775 |
| Abstract label (rationalized) | 0.852 | 0.632 | 0.725 |

## Limitations
- MAHNOB-HCI data for NegSl-AIS replication is currently unavailable locally.
- Full NegSl-AIS reproduction (including multimodal feature extraction) is pending.
- PubMed coverage is limited to a targeted subset of related-work queries logged in the publications log.

## Reproducibility
- Metrics log: `/Users/byron/projects/data/immunos_data/research/scifact/logs/metrics_dev.json`
- Dataset root: `/Users/byron/projects/data/immunos_data/research/scifact`
- Replication notes: `/Users/byron/projects/docs/reference/scifact-baseline-replication.md`
- Preprint project: `/Users/byron/projects/immunos-replication-preprint/`

## Data/Code Availability
- Data paths and logs are local under `/Users/byron/projects/data/immunos_data`.
- Code and documentation live under `/Users/byron/projects` with readable snapshots under `/Users/byron/projects/docs/code-snapshots/`.

## Ethics/Compliance
- Work is designed for local execution with a local LLM and an air-gapped intent to avoid external data exposure.
- Only public metadata and abstracts are used for literature references.
