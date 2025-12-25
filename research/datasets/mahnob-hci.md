---
title: MAHNOB-HCI
tags: [dataset, emotion, multimodal, biosignals]
type: dataset
---

# MAHNOB-HCI

## Overview
Multimodal emotion recognition dataset with synchronized physiological signals and visual stimuli responses.
Commonly used for valence/arousal classification and multimodal fusion baselines.

## Source
- Portal root: https://mahnob-db.eu/ (root page loads)
- HCI-Tagging page: https://mahnob-db.eu/hci-tagging/ (currently returns 404)
- iBUG MAHNOB-Implicit-Tagging page: https://ibug.doc.ic.ac.uk/research/implicit-tagging-database/ (active)
- Access: manual request + license acceptance

## Access Note (404 status)
As of 2025-12-24, the root MAHNOB page loads but the HCI-Tagging subpage returns 404.
If the page remains down, use the Internet Archive for the request details or contact
the dataset maintainers listed on the root page or the iBUG research page.

## Access (Stage 1: Request)
1. Visit the MAHNOB-HCI portal and register with an academic email.
2. Complete the dataset access request form (project summary + intended use).
3. Accept the license terms (non-commercial academic research).
4. Await approval and download credentials.

### Request Email (if required)
Subject: MAHNOB-HCI Dataset Access Request - immunOS Replication

Hello MAHNOB-HCI team,

I am requesting access to the MAHNOB-HCI dataset for academic research.
Project: immunOS replication of NegSl-AIS multimodal emotion classification results.
Institution: [Institution Name]
PI/Advisor: [Name, if applicable]
Intended use: Reproducing published baselines and validating feature-fusion methods.

Please let me know if any additional documentation is needed.

Thank you,
[Full Name]
[Role]
[Institution]
[Contact Email]

## Storage Plan (Stage 2)
- Local root: `/Users/byron/projects/data/immunos_data/emotion/mahnob-hci/`
- Do not commit data to Git.
- Record download date, size, and checksum in this file after retrieval.

## Verification (Stage 3)
- Confirm all modalities present (EEG, ECG, GSR, respiration, temperature).
- Verify sample counts and subject list match documentation.
- Capture a minimal manifest (file list + sizes).

## Preprocessing Plan (Stage 4)
- Follow NegSl-AIS paper: segmentation by status channel, ICA artifact removal,
  modality-specific filtering, 6s windows with overlap, per-modality normalization.
- Record any deviations from the paper with rationale.

## Citation
- Dataset citation as provided on the MAHNOB-HCI portal.
- Link to the NegSl-AIS baseline paper in the replication preprint.

## References
- MAHNOB-HCI-Tagging Database Manual (PDF): https://mahnob-db.eu/hci-tagging/media/uploads/manual.pdf
- Soleymani, M., Lichtenauer, J., Pun, T., Pantic, M. A Multimodal Database for Affect Recognition and Implicit Tagging. IEEE Transactions on Affective Computing, 3(1), 42-55 (2012). DOI: 10.1109/T-AFFC.2011.25.
- TorchEEG MAHNOBDataset docs: https://torcheeg.readthedocs.io/en/v1.1.0/generated/torcheeg.datasets.MAHNOBDataset.html
- iBUG MAHNOB-Implicit-Tagging page: https://ibug.doc.ic.ac.uk/research/implicit-tagging-database/

## Local Assets
- Manual: `/Users/byron/projects/research/datasets/mahnob-hci/assets/mahnob-hci-tagging-manual.pdf`
- Figures (iBUG page):
  - `/Users/byron/projects/research/datasets/mahnob-hci/assets/ImplicitTagging5_admin_big.jpg`
  - `/Users/byron/projects/research/datasets/mahnob-hci/assets/ImplicitTagging6_admin_big.jpg`
  - `/Users/byron/projects/research/datasets/mahnob-hci/assets/ImplicitTagging7_admin_big.jpg`
  - `/Users/byron/projects/research/datasets/mahnob-hci/assets/ImplicitTagging8_admin_big.png`

## TorchEEG Notes (dataset loader)
- Access still requires MAHNOB-HCI Sessions download; TorchEEG expects the unzipped `Sessions/` folder.
- Typical root path: `./Sessions` (contains `Part_*_Trial*_emotion.bdf` + `session.xml`).
- Signals: EEG 32 channels (512 Hz), peripheral physio (ECG, GSR, Temp, Resp at 256 Hz), eye gaze (60 Hz).
- Labels: arousal/valence/control/predictability on 1-9 scale; TorchEEG examples binarize at 5.0.
- Common transforms: band differential entropy + grid mapping, 2D tensor conversion, or GNN graph transform.

## Replication Summary (paper + manual)
- Participants: 27-30 participants; mixed gender, ages ~19-40; academic volunteers.
- Modalities: face video (6 cameras), audio (room + head-worn mics), eye gaze (Tobii), physiological (EEG 32ch, ECG 3ch, GSR, respiration, temperature).
- Experiment 1 (emotion elicitation): 20 emotional videos + neutral clips; self-report on emotion label, arousal, valence, dominance, predictability.
- Experiment 2 (implicit tagging): images/videos shown without tag, then with correct/incorrect tags; participant agreement recorded.
- File formats: physiological data in BDF; eye gaze in TSV; session metadata in XML; status channel encodes stimulus timing pulses.
- Key timing: status channel rising edge from 0 to 16 marks start/stop of stimuli; tagging has additional pulse for tag timing.
- Access: web-based system with EULA; SSL-protected downloads; storage must be firewall-protected (per manual).

## Status
- Access: pending (request not yet submitted)
- Local copy: not present
