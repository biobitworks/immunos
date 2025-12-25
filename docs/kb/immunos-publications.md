# IMMUNOS Publications & References

## Citation Format

Citations follow IEEE style for software/systems publications:
```
[#] A. Author, B. Author, "Title," Journal/Conference, vol. X, no. Y, pp. ZZ-ZZ, Year. DOI: url
```

---

## Core References

### Negative Selection AIS

[1] M. Umair et al., "Negative selection-based artificial immune system (NegSl-AIS)—A hybrid multimodal emotional effect classification model," *Results in Engineering*, vol. 27, art. 106601, 2025.
DOI: https://doi.org/10.1016/j.rineng.2025.106601

[2] Z. Ji et al., "Revisiting negative selection algorithms," *Evolutionary Computation*, vol. 15, no. 2, pp. 223-251, 2007.
DOI: https://doi.org/10.1162/evco.2007.15.2.223
PubMed: https://pubmed.ncbi.nlm.nih.gov/17535140/

[3] S. Forrest et al., "Self-nonself discrimination in a computer," *IEEE Symposium on Security and Privacy*, pp. 202-212, 1994.
DOI: https://doi.org/10.1109/RISP.1994.296580

### Claim Verification

[4] A. Ismail et al., "Retrieval augmented scientific claim verification," *JAMIA Open*, vol. 7, no. 1, art. ooae021, 2024.
DOI: https://doi.org/10.1093/jamiaopen/ooae021
PubMed: https://pubmed.ncbi.nlm.nih.gov/38455840/

[5] J. Zhu et al., "Valsci: an open-source, self-hostable literature review utility for automated large-batch scientific claim verification using large language models," *BMC Bioinformatics*, vol. 26, art. 159, 2025.
DOI: https://doi.org/10.1186/s12859-025-06159-4
PubMed: https://pubmed.ncbi.nlm.nih.gov/40437377/

[6] R. Pasumarthi et al., "Aggregating pairwise semantic differences for few-shot claim verification," *PeerJ Computer Science*, vol. 8, art. e1137, 2022.
DOI: https://doi.org/10.7717/peerj-cs.1137
PubMed: https://pubmed.ncbi.nlm.nih.gov/36426249/

[7] D. Wadden et al., "Fact or Fiction: Verifying Scientific Claims," *arXiv preprint*, arXiv:2004.14974, 2020.
arXiv: https://arxiv.org/abs/2004.14974

### Multimodal Emotion Recognition

[8] X. Qiu et al., "Multimodal Emotion Recognition Based on Facial Expressions, Speech, and EEG," *IEEE Open Journal of Engineering in Medicine and Biology*, vol. 5, pp. 396-403, 2024.
DOI: https://doi.org/10.1109/OJEMB.2023.3240280
PubMed: https://pubmed.ncbi.nlm.nih.gov/38899017/

[9] Y. Li et al., "Cross-modal credibility modelling for EEG-based multimodal emotion recognition," *Journal of Neural Engineering*, vol. 21, no. 2, art. 026040, 2024.
DOI: https://doi.org/10.1088/1741-2552/ad3987
PubMed: https://pubmed.ncbi.nlm.nih.gov/38565099/

[10] Y. Du et al., "Attention-based 3D convolutional recurrent neural network model for multimodal emotion recognition," *Frontiers in Neuroscience*, vol. 17, art. 1330077, 2023.
DOI: https://doi.org/10.3389/fnins.2023.1330077
PubMed: https://pubmed.ncbi.nlm.nih.gov/38268710/

### AIS in Biomedical Applications

[11] J. Ahmad et al., "Artificial Immune System-Negative Selection Classification Algorithm (NSCA) for Four Class Electroencephalogram (EEG) Signals," *Frontiers in Human Neuroscience*, vol. 12, art. 439, 2018.
DOI: https://doi.org/10.3389/fnhum.2018.00439
PubMed: https://pubmed.ncbi.nlm.nih.gov/30524257/

[12] R. Torkamani-Azar et al., "RAIRS2 a new expert system for diagnosing tuberculosis with real-world tournament selection mechanism inside artificial immune recognition system," *Medical & Biological Engineering & Computing*, vol. 54, pp. 385-399, 2016.
DOI: https://doi.org/10.1007/s11517-015-1323-6
PubMed: https://pubmed.ncbi.nlm.nih.gov/26081904/

---

## Datasets

### Public Datasets Used

| Dataset | Domain | Source | License |
|---------|--------|--------|---------|
| **SciFact** | Claim Verification | [AllenAI](https://github.com/allenai/scifact) | Apache 2.0 |
| **TruthfulQA** | Hallucination Detection | [HuggingFace](https://huggingface.co/datasets/truthful_qa) | Apache 2.0 |
| **NSL-KDD** | Network Intrusion | [Canadian Institute for Cybersecurity](https://www.unb.ca/cic/datasets/nsl.html) | Research |
| **MAHNOB-HCI** | Emotion Recognition | [Mahnob Database](https://mahnob-db.eu/hci-tagging/) | Academic |
| **DEAP** | Emotion Recognition | [QMUL](http://www.eecs.qmul.ac.uk/mmv/datasets/deap/) | Academic |

### Dataset Locations

```
data/immunos_data/
├── hallucination/     # TruthfulQA samples
├── network/           # NSL-KDD samples
├── research/          # SciFact samples
│   └── scifact/       # Full SciFact dataset
└── code/              # Code samples
```

### SciFact Dataset Details

- **Claims**: 1,409 scientific claims
- **Corpus**: 5,183 abstracts from S2ORC
- **Labels**: SUPPORTS, REFUTES, NOT_ENOUGH_INFO
- **Note**: Test labels not publicly released

Download:
```bash
git clone https://github.com/allenai/scifact.git data/immunos_data/research/scifact
```

---

## Software & Tools

### IMMUNOS Components

| Component | Version | Purpose |
|-----------|---------|---------|
| immunos_dashboard.py | v20251225 | Web dashboard |
| immunos_recover.py | v1.0 | Context recovery |
| immunos_snapshot.py | v1.0 | Session snapshots |
| immunos_memory.py | v1.0 | Memory management |
| immunos_nk_scan.py | v1.0 | NK Cell scanner |

### External Dependencies

| Package | Version | License | URL |
|---------|---------|---------|-----|
| Flask | 3.0+ | BSD-3 | https://flask.palletsprojects.com/ |
| Flask-SocketIO | 5.0+ | MIT | https://flask-socketio.readthedocs.io/ |
| Ollama | 0.3+ | MIT | https://ollama.ai/ |
| TailwindCSS | 3.0+ | MIT | https://tailwindcss.com/ |
| Chart.js | 4.0+ | MIT | https://www.chartjs.org/ |

### Model Sources

| Model | Source | License | Size |
|-------|--------|---------|------|
| qwen2.5:1.5b | [Ollama](https://ollama.ai/library/qwen2.5) | Apache 2.0 | 986 MB |
| qwen2.5:7b | [Ollama](https://ollama.ai/library/qwen2.5) | Apache 2.0 | 4.7 GB |
| qwen2.5-coder:7b | [Ollama](https://ollama.ai/library/qwen2.5-coder) | Apache 2.0 | 4.7 GB |
| deepseek-r1:14b | [Ollama](https://ollama.ai/library/deepseek-r1) | DeepSeek | 9.0 GB |
| bge-m3 | [Ollama](https://ollama.ai/library/bge-m3) | MIT | 1.2 GB |

---

## Related Projects

### Prion Clock

Aging biology hypothesis exploring protein-based cellular timekeeping.

- Location: `/Users/byron/projects/prion-clock/`
- Status: Ready for bioRxiv submission
- Citations: 30 verified

### Papers Archive

Literature collection for aging biology research.

- Location: `/Users/byron/projects/papers/`
- Integration: Zotero, automated figure download

---

## BibTeX Export

For LaTeX documents, citations available at:
```
/Users/byron/projects/research/citations/immunos-references.bib
```

Generate BibTeX from DOI:
```bash
curl -LH "Accept: application/x-bibtex" https://doi.org/10.1016/j.rineng.2025.106601
```

---

## Adding New References

1. Add to publications log:
   ```bash
   vim docs/reference/publications-log.md
   ```

2. Format as IEEE citation

3. Verify DOI resolves

4. Check PubMed if applicable:
   ```
   https://pubmed.ncbi.nlm.nih.gov/?term=DOI
   ```

5. Update this KB file

---

**Last Updated**: 2025-12-25
**Format**: IEEE Style
**Total References**: 12 core, multiple supporting
