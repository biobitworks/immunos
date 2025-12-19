# IMMUNOS Dataset Reference Guide

## M1 MacBook Pro Compatibility Summary
- **RAM**: 16GB unified memory
- **Storage**: 1TB SSD
- **GPU**: 8-core Apple GPU (MLX/Metal optimized)
- **Recommended batch sizes**: 32-128 depending on model size

---

## 1. EMOTION RECOGNITION DATASETS

### Small: CK+ (Extended Cohn-Kanade)
| Property | Value |
|----------|-------|
| **Size** | 593 image sequences |
| **Format** | PNG images |
| **Download** | http://www.jeffcohn.net/Resources/ |
| **Classes** | 7 emotions (anger, contempt, disgust, fear, happy, sadness, surprise) |
| **M1 Training** | ~5 minutes |
| **Notes** | Lab-controlled, high quality, good for initial testing |

### Medium: FER2013
| Property | Value |
|----------|-------|
| **Size** | 35,887 images (48x48 grayscale) |
| **Format** | CSV (pixel values) |
| **Download** | https://www.kaggle.com/datasets/msambare/fer2013 |
| **Classes** | 7 emotions |
| **M1 Training** | ~30-60 minutes |
| **Storage** | ~100MB |
| **Notes** | Standard benchmark, some label noise |

### Large: AffectNet
| Property | Value |
|----------|-------|
| **Size** | 1M+ images |
| **Format** | JPEG images |
| **Download** | http://mohammadmahoor.com/affectnet/ (academic) |
| **Classes** | 8 emotions + valence/arousal |
| **M1 Training** | Recommend cloud GPU or subset |
| **Storage** | ~120GB |
| **Notes** | Academic license required, most comprehensive |

---

## 2. LLM HALLUCINATION DETECTION DATASETS

### Small: TruthfulQA
| Property | Value |
|----------|-------|
| **Size** | 817 questions |
| **Format** | JSON |
| **Download** | https://github.com/sylinrl/TruthfulQA |
| **Classes** | True/False/Informative |
| **M1 Training** | ~10 minutes |
| **Notes** | Tests if LLMs mimic human falsehoods |

### Medium: HaluEval
| Property | Value |
|----------|-------|
| **Size** | 35,000 samples |
| **Format** | JSON |
| **Download** | https://github.com/RUCAIBox/HaluEval |
| **Tasks** | QA, Dialogue, Summarization |
| **M1 Training** | ~1-2 hours |
| **Storage** | ~500MB |
| **Notes** | Comprehensive hallucination benchmark |

### Additional: FEVER
| Property | Value |
|----------|-------|
| **Size** | 185,445 claims |
| **Format** | JSON |
| **Download** | https://fever.ai/dataset/fever.html |
| **Classes** | Supported, Refuted, NEI |
| **M1 Training** | ~2-3 hours |
| **Notes** | Wikipedia-based fact verification |

### Medical: MedHallu
| Property | Value |
|----------|-------|
| **Size** | 10,000 QA pairs |
| **Format** | JSON |
| **Download** | Derived from PubMedQA |
| **Notes** | Medical domain hallucinations |

---

## 3. NETWORK INTRUSION DETECTION DATASETS

### Small: NSL-KDD
| Property | Value |
|----------|-------|
| **Size** | 148,517 records |
| **Format** | CSV/ARFF |
| **Download** | https://www.unb.ca/cic/datasets/nsl.html |
| **Features** | 41 features |
| **Classes** | Normal + 4 attack types (DoS, Probe, R2L, U2R) |
| **M1 Training** | ~15-30 minutes |
| **Storage** | ~25MB |
| **Notes** | Cleaned version of KDD Cup 99, still widely used |

### Medium: CICIDS2017
| Property | Value |
|----------|-------|
| **Size** | 2.8M flows |
| **Format** | CSV (per-day files) |
| **Download** | https://www.unb.ca/cic/datasets/ids-2017.html |
| **Features** | 80+ features |
| **Attacks** | Brute Force, DoS, DDoS, Web Attacks, Botnet, etc. |
| **M1 Training** | ~2-4 hours |
| **Storage** | ~8GB |
| **Notes** | Modern attacks, realistic traffic |

### Large: UNSW-NB15
| Property | Value |
|----------|-------|
| **Size** | 2.5M records |
| **Format** | CSV |
| **Download** | https://research.unsw.edu.au/projects/unsw-nb15-dataset |
| **Features** | 49 features |
| **Classes** | 9 attack categories |
| **M1 Training** | ~3-5 hours |
| **Storage** | ~2GB |
| **Notes** | Most modern patterns, includes hybrid attacks |

---

## 4. CODE VULNERABILITY DETECTION DATASETS

### Small: Devign
| Property | Value |
|----------|-------|
| **Size** | 27,318 functions |
| **Format** | JSON |
| **Download** | https://sites.google.com/view/devign |
| **Language** | C/C++ |
| **Projects** | FFmpeg, QEMU |
| **M1 Training** | ~30-60 minutes |
| **Notes** | High-quality manual labels, limited scope |

### Medium: BigVul
| Property | Value |
|----------|-------|
| **Size** | 188,636 functions |
| **Format** | CSV/JSON |
| **Download** | https://github.com/ZeoVan/MSR_20_Code_vulnerability_CSV_Dataset |
| **Language** | C/C++ |
| **CWEs** | 91 types |
| **M1 Training** | ~2-4 hours |
| **Storage** | ~1GB |
| **Notes** | Large scale but has label noise (~25% accuracy) |

### High-Quality: DiverseVul
| Property | Value |
|----------|-------|
| **Size** | 18,945 vulnerable functions |
| **Format** | JSON |
| **Download** | https://github.com/wagner-group/diversevul |
| **Language** | C/C++ |
| **CWEs** | 150 types |
| **M1 Training** | ~1-2 hours |
| **Notes** | Best label quality, recommended for production |

### Synthetic: SARD (Juliet Test Suite)
| Property | Value |
|----------|-------|
| **Size** | 64,000+ test cases |
| **Format** | Source files |
| **Download** | https://samate.nist.gov/SARD/ |
| **Language** | C/C++, Java |
| **Notes** | Synthetic but comprehensive CWE coverage |

---

## 5. RESEARCH PAPER VERIFICATION DATASETS

### Small: SciFact
| Property | Value |
|----------|-------|
| **Size** | 1,400 claims + 5,183 abstracts |
| **Format** | JSON |
| **Download** | https://github.com/allenai/scifact |
| **Classes** | Supports, Refutes, NEI |
| **M1 Training** | ~20-30 minutes |
| **Notes** | Expert-written scientific claims |

### Medium: SciFact-Open
| Property | Value |
|----------|-------|
| **Size** | 500K abstracts |
| **Format** | JSON |
| **Download** | https://github.com/allenai/scifact |
| **Notes** | Open-domain evaluation |

### Medical: PubMedQA
| Property | Value |
|----------|-------|
| **Size** | 211,269 QA pairs |
| **Format** | JSON |
| **Download** | https://pubmedqa.github.io/ |
| **Classes** | Yes, No, Maybe |
| **M1 Training** | ~1-2 hours |
| **Notes** | Biomedical research questions |

### Multi-modal: SCITAB
| Property | Value |
|----------|-------|
| **Size** | 1,200 claims with tables |
| **Format** | JSON + Tables |
| **Download** | Research paper |
| **Notes** | Compositional reasoning required |

---

## 6. LOGIC & REASONING DATASETS

### Small: LogiQA
| Property | Value |
|----------|-------|
| **Size** | 8,678 QA pairs |
| **Format** | JSON |
| **Download** | https://github.com/lgw863/LogiQA-dataset |
| **Notes** | Logical reasoning in natural language |

### Medium: ReClor
| Property | Value |
|----------|-------|
| **Size** | 6,138 questions |
| **Format** | JSON |
| **Download** | https://whyu.me/reclor/ |
| **Notes** | Reading comprehension + logical reasoning |

### Mathematical: GSM8K
| Property | Value |
|----------|-------|
| **Size** | 8,500 math word problems |
| **Format** | JSON |
| **Download** | https://github.com/openai/grade-school-math |
| **Notes** | Grade school math, tests reasoning chains |

---

## Quick Start Download Script

```bash
#!/bin/bash
# download_datasets.sh - Download starter datasets for IMMUNOS

mkdir -p ~/immunos_data/{emotion,hallucination,network,code,research}

# FER2013 (requires Kaggle API)
echo "Note: FER2013 requires Kaggle account"
# kaggle datasets download -d msambare/fer2013 -p ~/immunos_data/emotion/

# TruthfulQA
git clone https://github.com/sylinrl/TruthfulQA ~/immunos_data/hallucination/truthfulqa

# HaluEval
git clone https://github.com/RUCAIBox/HaluEval ~/immunos_data/hallucination/halueval

# NSL-KDD
cd ~/immunos_data/network
wget https://www.unb.ca/cic/datasets/nsl.html -O nsl_info.html
echo "Download NSL-KDD from: https://www.unb.ca/cic/datasets/nsl.html"

# SciFact
git clone https://github.com/allenai/scifact ~/immunos_data/research/scifact

# Devign
echo "Download Devign from: https://sites.google.com/view/devign"

# DiverseVul
git clone https://github.com/wagner-group/diversevul ~/immunos_data/code/diversevul

echo "Datasets downloaded to ~/immunos_data/"
ls -la ~/immunos_data/*/
```

---

## M1 MacBook Training Time Estimates

| Dataset | Size | Est. Training Time | RAM Usage |
|---------|------|-------------------|-----------|
| CK+ | 593 | 5 min | 2GB |
| FER2013 | 35K | 30-60 min | 4GB |
| TruthfulQA | 817 | 10 min | 2GB |
| HaluEval | 35K | 1-2 hr | 6GB |
| NSL-KDD | 148K | 15-30 min | 4GB |
| CICIDS2017 | 2.8M | 2-4 hr | 8GB |
| Devign | 27K | 30-60 min | 4GB |
| DiverseVul | 18K | 1-2 hr | 6GB |
| SciFact | 1.4K | 20-30 min | 3GB |
| PubMedQA | 211K | 1-2 hr | 6GB |

---

## Cloud GPU Training (Optional)

For larger datasets or faster training:

### Google Colab (Free Tier)
- GPU: T4 (16GB)
- Suitable for: Medium datasets
- Time limit: ~12 hours continuous

### RunPod
- A100 40GB: ~$1.89/hr
- Suitable for: Large datasets
- No time limits

### Lambda Labs
- A100 80GB: ~$1.10/hr
- Best for: Very large models

### Training Command (Cloud)
```bash
# Install IMMUNOS
pip install immunos[cloud]

# Train with GPU
immunos train all \
  --datasets-dir /content/datasets \
  --output-dir /content/detectors \
  --device cuda \
  --batch-size 256 \
  --epochs 100
```

---

## Dataset Preprocessing Notes

### Emotion Recognition
- Resize images to 48x48 or 224x224
- Normalize pixel values to [0, 1]
- Apply data augmentation (rotation, flip, brightness)

### Network Traffic
- Encode categorical features (protocol, service)
- Normalize numerical features (min-max)
- Handle missing values

### Code Vulnerability
- Parse AST using tree-sitter or Joern
- Extract code property graphs
- Tokenize with CodeBERT tokenizer

### Research Papers
- Use SciBERT for embeddings
- Extract claim-evidence pairs
- Handle citation context

---

## Recommended Learning Path

1. **Start Small**: CK+, TruthfulQA, NSL-KDD, Devign, SciFact
2. **Scale Up**: FER2013, HaluEval, CICIDS2017, DiverseVul
3. **Production**: AffectNet subset, full UNSW-NB15, MegaVul, SciFact-Open
