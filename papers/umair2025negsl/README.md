# Negative selection-based artificial immune system (NegSl-AIS)

## Summary
Contains 1 subdirectories and 32 files. Key subfolders: opus_implementation/.


**FOUNDATIONAL PAPER FOR IMMUNOS**

## Citation

Umair, M., Rashid, N., Khan, U.S., Hamza, A., Zeb, A., Nawaz, T.H., & Ansari, A.R. (2025).
Negative selection-based artificial immune system (NegSl-AIS)--A hybrid multimodal emotional effect classification model.
*Results in Engineering*, 27, 106601.
https://doi.org/10.1016/j.rineng.2025.106601

## Significance to IMMUNOS

This paper provides the **theoretical and algorithmic foundation** for the IMMUNOS artificial immune system.
The negative selection algorithm, originally designed for multimodal emotion classification, translates directly
to code security anomaly detection.

**Key algorithm**: Negative selection mimics biological T-cell maturation to distinguish "self" (safe) from "non-self" (threat).

## Files

- `1-s2.0-S2590123025026702-main.pdf` - Main paper
- `supplementary.pdf` - Supplementary materials
- `1-s2.0-S2590123025026702-gr*.jpg` - Figures 1-15

## Key Concepts

- **Negative selection**: Train detectors on "safe" samples, detect anything that doesn't match
- **Multimodal fusion**: Combine multiple data sources with weighted biasing
- **Feature extraction**: Statistical, frequency-domain, and LSTM-based features
- **Detector optimization**: Empirical tuning of detector count and thresholds
- **Generalization error**: Track overfitting with train/test error differences

## Direct Mapping to IMMUNOS

| Paper (Emotion) | IMMUNOS (Security) |
|-----------------|-------------------|
| Self samples = Specific emotion | Self samples = Safe code |
| Non-self = Other emotions | Non-self = Threats/anomalies |
| EEG, ECG, GSR modalities | Code, structure, metadata modalities |
| Rself threshold = 0.87-1.34 | Severity thresholds (HIGH/MEDIUM/LOW) |
| 96.48% accuracy | Target: >95% threat detection |

## Full Analysis

See: `/Users/byron/projects/.immunos/journal/IMMUNOS_FOUNDATION_PAPER.md`

**Date added**: 2025-12-18

## Directory Map
```
umair2025negsl/
├── opus_implementation/
├── 1-s2.0-S2590123025026702-gr1.jpg
├── 1-s2.0-S2590123025026702-gr10.jpg
├── 1-s2.0-S2590123025026702-gr10_lrg.jpg
├── 1-s2.0-S2590123025026702-gr11.jpg
├── 1-s2.0-S2590123025026702-gr11_lrg.jpg
├── 1-s2.0-S2590123025026702-gr12.jpg
├── 1-s2.0-S2590123025026702-gr12_lrg.jpg
├── 1-s2.0-S2590123025026702-gr13.jpg
├── 1-s2.0-S2590123025026702-gr13_lrg.jpg
├── 1-s2.0-S2590123025026702-gr14.jpg
├── 1-s2.0-S2590123025026702-gr14_lrg.jpg
├── 1-s2.0-S2590123025026702-gr15.jpg
├── 1-s2.0-S2590123025026702-gr15_lrg.jpg
├── 1-s2.0-S2590123025026702-gr1_lrg.jpg
├── 1-s2.0-S2590123025026702-gr2.jpg
├── 1-s2.0-S2590123025026702-gr2_lrg.jpg
├── 1-s2.0-S2590123025026702-gr3.jpg
├── 1-s2.0-S2590123025026702-gr3_lrg.jpg
├── 1-s2.0-S2590123025026702-gr4.jpg
├── 1-s2.0-S2590123025026702-gr4_lrg.jpg
├── 1-s2.0-S2590123025026702-gr5.jpg
├── 1-s2.0-S2590123025026702-gr5_lrg.jpg
├── 1-s2.0-S2590123025026702-gr6.jpg
├── 1-s2.0-S2590123025026702-gr6_lrg.jpg
├── 1-s2.0-S2590123025026702-gr7.jpg
├── 1-s2.0-S2590123025026702-gr7_lrg.jpg
├── 1-s2.0-S2590123025026702-gr8.jpg
├── 1-s2.0-S2590123025026702-gr8_lrg.jpg
├── 1-s2.0-S2590123025026702-gr9.jpg
├── 1-s2.0-S2590123025026702-gr9_lrg.jpg
├── 1-s2.0-S2590123025026702-main.pdf
└── supplementary.pdf
```
