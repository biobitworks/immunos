---
project: immunos
type: reference-log
title: Publications log (non-BibTeX)
last_updated: 2025-12-23
---

# Publications Log (Non-BibTeX)

This is a running list of publications and references outside BibTeX.
When the user requests publication additions, log them here and in the daily journal.

## Transfer

- Planned: convert this log to BibTeX via a helper script (TBD).
- Current BibTeX location for IMMUNOS: `/Users/byron/projects/research/citations/immunos-references.bib`

## 2025-12-23
- Umair, M. et al. (2025). Negative selection-based artificial immune system (NegSl-AIS)--A hybrid multimodal emotional effect classification model. Results in Engineering 27, 106601.
  DOI: https://doi.org/10.1016/j.rineng.2025.106601
  PubMed: https://pubmed.ncbi.nlm.nih.gov/?term=10.1016%2Fj.rineng.2025.106601
  Summary: Proposes a negative selection AIS for multimodal emotion classification, mapping self/non-self separation to detector training and tuning thresholds for generalization.
  Relevance: core — foundational NegSl-AIS study anchoring the replication framing and dataset context.
- Ismail, A. et al. (2024). Retrieval augmented scientific claim verification. JAMIA Open.
  DOI: https://doi.org/10.1093/jamiaopen/ooae021
  PubMed: https://pubmed.ncbi.nlm.nih.gov/38455840/
  Summary: Presents CliVER, a retrieval-augmented claim verification system using PubMed abstracts, PICO-based rationale extraction, and an ensemble labeler; evaluated on CoVERt and SciFact with clinician comparisons.
  Relevance: core — direct SciFact/PubMed claim verification baseline aligned with the replication pipeline.
- Zhu, J. et al. (2025). Valsci: an open-source, self-hostable literature review utility for automated large-batch scientific claim verification using large language models. BMC Bioinformatics.
  DOI: https://doi.org/10.1186/s12859-025-06159-4
  PubMed: https://pubmed.ncbi.nlm.nih.gov/40437377/
  Summary: Describes a self-hosted, retrieval-augmented LLM workflow for claim verification with bibliometric scoring; SciFact benchmarks show lower citation hallucination rates versus base LLM outputs.
  Relevance: supporting — recent large-batch verification system evaluated on SciFact, useful for comparative context.
- Pasumarthi, R. et al. (2022). Aggregating pairwise semantic differences for few-shot claim verification. PeerJ Computer Science.
  DOI: https://doi.org/10.7717/peerj-cs.1137
  PubMed: https://pubmed.ncbi.nlm.nih.gov/36426249/
  Summary: Introduces SEED, which aggregates semantic difference vectors for few-shot claim verification and improves performance on FEVER and SciFact.
  Relevance: core — SciFact-focused claim verification method relevant to baseline comparisons.
- Qiu, X. et al. (2024). Multimodal Emotion Recognition Based on Facial Expressions, Speech, and EEG. IEEE Open Journal of Engineering in Medicine and Biology.
  DOI: https://doi.org/10.1109/OJEMB.2023.3240280
  PubMed: https://pubmed.ncbi.nlm.nih.gov/38899017/
  Summary: Proposes a three-branch multimodal system with decision-level fusion across facial, speech, and EEG features; evaluated on CK+, EMO-DB, and MAHNOB-HCI.
  Relevance: supporting — MAHNOB-HCI benchmark context for NegSl-AIS replication scope.
- Li, Y. et al. (2024). Cross-modal credibility modelling for EEG-based multimodal emotion recognition. Journal of Neural Engineering.
  DOI: https://doi.org/10.1088/1741-2552/ad3987
  PubMed: https://pubmed.ncbi.nlm.nih.gov/38565099/
  Summary: Uses intra/inter-modal attention and sequential consistency to improve fusion credibility, reporting gains on DEAP and MAHNOB-HCI.
  Relevance: supporting — establishes MAHNOB-HCI multimodal baselines relevant to NegSl-AIS framing.
- Tadas, C. et al. (2018). Emotion Recognition from EEG and Facial Expressions: a Multimodal Approach. EMBC 2018.
  DOI: https://doi.org/10.1109/EMBC.2018.8512407
  PubMed: https://pubmed.ncbi.nlm.nih.gov/30440451/
  Summary: Applies feature-level fusion of EEG and facial microexpressions with NN and random forest classifiers on MAHNOB-HCI, improving accuracy over unimodal baselines.
  Relevance: supporting — early MAHNOB-HCI fusion baseline for replication comparisons.
- Ahmad, J. et al. (2018). Artificial Immune System-Negative Selection Classification Algorithm (NSCA) for Four Class Electroencephalogram (EEG) Signals. Frontiers in Human Neuroscience.
  DOI: https://doi.org/10.3389/fnhum.2018.00439
  PubMed: https://pubmed.ncbi.nlm.nih.gov/30524257/
  Summary: Uses MFCC features, stacked auto-encoders, and GA-optimized negative selection detectors to classify four-class motor imagery EEG with reported mean accuracy of 86.39%.
  Relevance: supporting — demonstrates negative selection in EEG classification, informing AIS design parallels.
- Torkamani-Azar, R. et al. (2016). RAIRS2 a new expert system for diagnosing tuberculosis with real-world tournament selection mechanism inside artificial immune recognition system. Medical & Biological Engineering & Computing.
  DOI: https://doi.org/10.1007/s11517-015-1323-6
  PubMed: https://pubmed.ncbi.nlm.nih.gov/26081904/
  Summary: Introduces a hybrid AIRS classifier with tournament selection for TB diagnosis, reporting high accuracy on a 175-record clinical dataset.
  Relevance: adjacent — AIS classifier application in biomedical diagnosis outside the current replication scope.
- Cerulo, L. et al. (2013). A negative selection heuristic to predict new transcriptional targets. BMC Bioinformatics.
  DOI: https://doi.org/10.1186/1471-2105-14-S1-S3
  PubMed: https://pubmed.ncbi.nlm.nih.gov/23368951/
  Summary: Improves supervised transcriptional target prediction by selecting reliable negative examples from unlabeled data, boosting classification performance in regulatory network inference.
  Relevance: adjacent — negative selection heuristic in bioinformatics, conceptually related but different task domain.
- Du, Y. et al. (2023). Attention-based 3D convolutional recurrent neural network model for multimodal emotion recognition. Frontiers in Neuroscience.
  DOI: https://doi.org/10.3389/fnins.2023.1330077
  PubMed: https://pubmed.ncbi.nlm.nih.gov/38268710/
  Summary: Proposes a 3D CNN-CRNN with attention for EEG and visual fusion, reporting strong accuracy on DEAP and MAHNOB-HCI.
  Relevance: supporting — recent MAHNOB-HCI results that frame NegSl-AIS target performance.
- Zhao, Y. et al. (2021). Expression EEG Multimodal Emotion Recognition Method Based on the Bidirectional LSTM and Attention Mechanism. Computational and Mathematical Methods in Medicine.
  DOI: https://doi.org/10.1155/2021/9967592
  PubMed: https://pubmed.ncbi.nlm.nih.gov/34055043/
  Summary: Uses BCN feature extraction and BiLSTM attention fusion for expression+EEG signals, evaluated on MAHNOB-HCI and DEAP.
  Relevance: supporting — additional MAHNOB-HCI multimodal baseline for replication context.
- Ji, Z. et al. (2007). Revisiting negative selection algorithms. Evolutionary Computation.
  DOI: https://doi.org/10.1162/evco.2007.15.2.223
  PubMed: https://pubmed.ncbi.nlm.nih.gov/17535140/
  Summary: Reviews negative selection algorithms, detailing representation, matching rules, and hybrid variants in AIS anomaly detection.
  Relevance: supporting — foundational NSA survey supporting methodological framing.
- Li, Z. et al. (2025). Generating detectors from anomaly samples via negative selection for network intrusion detection. Scientific Reports.
  DOI: https://doi.org/10.1038/s41598-025-20516-6
  PubMed: https://pubmed.ncbi.nlm.nih.gov/41107399/
  Summary: Generates detectors from anomaly samples to better cover low-dimensional subspaces, outperforming baselines on NSL-KDD and UNSW-NB15.
  Relevance: adjacent — NSA detector generation in intrusion detection, relevant method but different domain.
- Karakose, M. et al. (2013). Reinforcement learning based artificial immune classifier. TheScientificWorldJournal.
  DOI: https://doi.org/10.1155/2013/581846
  PubMed: https://pubmed.ncbi.nlm.nih.gov/23935424/
  Summary: Introduces a reinforcement learning AIS classifier and compares performance against negative selection and other AIS variants.
  Relevance: adjacent — AIS classifier variant for broader context rather than direct replication.

## 2025-12-22
- V-Detector disease diagnosis (book chapter).
  DOI: https://doi.org/10.1016/B978-0-12-803468-2.00011-4
  Summary: User-provided abstract; V-Detector NSA variant for disease detection with high reported accuracy.
- Applied Soft Computing AIS review.
  DOI: https://doi.org/10.1016/j.asoc.2021.108129
  PubMed: https://pubmed.ncbi.nlm.nih.gov/?term=10.1016%2Fj.asoc.2021.108129
  Summary: User-provided reference for AIS survey/modeling.
- Immunological computation (NCBI Bookshelf, Chapter 23).
  Link: https://www.ncbi.nlm.nih.gov/books/NBK459484/
  Summary: Immunology background for AIS modeling.
- Negative Selection Algorithms review.
  DOI: https://doi.org/10.1162/106365601750190975
  PubMed: https://pubmed.ncbi.nlm.nih.gov/?term=10.1162%2F106365601750190975
  Summary: Core NSA survey referenced in project docs.
- SciFact dataset paper.
  Link: https://arxiv.org/abs/2004.14974
  Summary: Scientific claim verification dataset baseline; test labels not public.
- Predicting emotion with biosignals (valence/arousal).
  DOI: https://doi.org/10.3390/s23031598
  PubMed: https://pubmed.ncbi.nlm.nih.gov/?term=10.3390%2Fs23031598
  Summary: Confirmed DOI from user; cited for emotion recognition baselines.
