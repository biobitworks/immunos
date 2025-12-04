# 🧬 Aging Biology Research Dashboard

**Project Goal**: Build comprehensive knowledge base for aging biology research, integrate longevity databases, track researchers and organizations, enable systems-level analysis.

---

## 👥 Researchers

### Active Tracking

**[[researchers/lidsky-peter|Peter V. Lidsky]]** ⭐ Priority
- Professor, City University of Hong Kong
- **Novel theory**: Aging as evolved pathogen control mechanism
- **Research**: Cellular senescence as antiviral defense ("immune militia")
- **Relevance**: High - reframes CellAge senescence genes as immune system
- **Status**: Potential collaboration
- **Connection**: [[literature/shvarts-2002-bcl6-senescence|BCL6 senescence paper]]

### Index
📋 [[researchers/README|All Researchers]] - Complete profiles and tracking

---

## 🏛️ Organizations

### Universities & Institutes

**[[organizations/city-university-hong-kong|City University of Hong Kong]]**
- Peter Lidsky's institution
- Biomedical Sciences department
- Novel aging biology perspectives

**[[organizations/monarch-initiative|Monarch Initiative]]**
- 33-database knowledge graph consortium
- 2.9M nodes, 160M edges
- Cross-species phenotype-genotype integration
- **Use**: Enrich GenAge with phenotypes

### Data Repositories

**[[organizations/longevity-db-huggingface|longevity-db (HuggingFace)]]**
- Aging datasets on HuggingFace Hub
- **Action**: Upload GenAge and CellAge here
- Community engagement platform

**[[organizations/future-house|Future House]]**
- AI for science research
- Created LAB-Bench (2,457 biology questions)
- **Use**: Evaluate AI analysis tools

### Index
📋 [[organizations/README|All Organizations]] - Complete list and profiles

---

## 💾 Databases

### Data Hosting & Distribution

**[[databases/huggingface-hub|HuggingFace Hub]]** 🎯 Next Action
- ML dataset hosting platform
- **Plan**: Upload GenAge (307 genes) and CellAge (950 genes)
- **Benefit**: Community access, version control, citations
- **Status**: Documentation complete, ready for upload

### Spatial Transcriptomics

**[[databases/spatiallibd|spatialLIBD]]** 🔍 Query System Available
- Brain spatial transcriptomics (DLPFC)
- **Size**: ~2 GB (33,538 genes, 47,681 spots)
- **Decision**: Use web interface + query scripts (NO download)
- **Scripts**:
  - [[../scripts/query_spatiallibd.py|Python query script]]
  - [[../scripts/query_spatiallibd.R|R query script]]
- **Questions**: 5 predefined aging biology questions

### Protein Interactions

**[[databases/biogrid|BioGRID]]** 🔗 Network Analysis
- 2.9M+ protein/genetic interactions
- **Size**: ~500 MB (50 MB human subset)
- **Use**: Map aging gene networks, identify hubs
- **API**: Free with key registration
- **Priority**: High for systems analysis

### Index
📋 [[databases/README|All Databases]] - Complete database documentation

---

## 📊 Datasets

### Core Aging Databases

**GenAge Human** (307 genes)
- Location: `/data/genage/human/genage_human.csv`
- Manually curated aging-related genes
- Build 21 (2024)
- [[../data/genage/human/README|Full documentation]]
- 📤 **Ready for HuggingFace upload**

**CellAge** (950 genes)
- Location: `/data/cellage/cellage3.tsv`
- Cellular senescence genes
- Build 3 (2023)
- [[../data/cellage/README|Full documentation]]
- 📤 **Ready for HuggingFace upload**

**CellAge Expression Signatures** (1,259 genes)
- Location: `/data/cellage/signatures1.csv`
- Meta-analysis across datasets
- Senescence expression patterns

### Index
📋 [[../data/README|All Datasets]] - Data directory overview

---

## 📚 Literature

### Key Papers

**[[literature/shvarts-2002-bcl6-senescence|BCL6 as Senescence Inhibitor]]** (2002)
- Shvarts et al., Genes & Development
- BCL6 bypasses p19ARF-p53 senescence pathway
- **Relevance**: CellAge gene, Lidsky's research interest
- **Finding**: BCL6 immortalizes cells via cyclin D1 induction

### Index
📋 [[literature/README|All Literature]] - Complete literature notes
📋 [[../research/citations/README|Citation Management]] - BibTeX files

---

## 🔧 Tools & Scripts

### Query Systems

**spatialLIBD Query** [[../scripts/README_SPATIALLIBD_QUERIES|Documentation]]
- Query brain spatial data without 2GB download
- 5 predefined aging biology questions
- Python & R interfaces
- Integrates GenAge/CellAge genes

### Planned Tools

- [ ] BioGRID network analysis script
- [ ] GenAge-CellAge overlap analyzer
- [ ] HuggingFace upload automation
- [ ] Monarch Initiative integration

---

## 🎯 Current Projects

### Active

**1. Aging Biology Infrastructure** 📍 In Progress
- **Status**: Phase 2 complete
- **Done**:
  - ✅ Directory structure (6 dirs)
  - ✅ Templates (4 files)
  - ✅ Researcher profiles (1 complete)
  - ✅ Organization profiles (4 complete)
  - ✅ Database docs (3 complete)
  - ✅ spatialLIBD query system
- **Next**:
  - HuggingFace authentication & upload
  - More database documentation
  - LAB-Bench dataset docs

**2. GenAge-CellAge Integration** 📍 Planned
- Cross-reference aging and senescence genes
- Identify overlap and unique genes
- Network analysis with BioGRID
- Pathway enrichment

**3. HuggingFace Data Publication** 🎯 Next
- Upload GenAge (307 genes)
- Upload CellAge (950 genes + 1,259 signatures)
- Create dataset cards with citations
- Enable community access

### Index
📋 [[projects/README|All Projects]] - Project tracking

---

## 🤝 Collaborations

### Potential

**Peter Lidsky** (CityU Hong Kong)
- Pathogen control theory validation with our data
- CellAge analysis under immune defense framework
- Evolutionary signatures in senescence genes

### Index
📋 [[collaborations/README|All Collaborations]] - Partnership tracking

---

## 📈 Metrics & Progress

**Infrastructure**:
- 6 directories created
- 4 templates
- 21 documentation files
- 5,865 lines documented
- 3 git commits

**Databases Tracked**:
- 3 documented (HuggingFace, spatialLIBD, BioGRID)
- 3 pending (Monarch KG, Dipper, RegulomeDB)

**Data Assets**:
- 307 human aging genes
- 2,205 model organism longevity genes
- 950 cellular senescence genes
- 1,259 senescence expression signatures
- Total: 4,721 genes tracked

---

## 🔍 Quick Search

**By Tag**:
- #researcher - Researcher profiles
- #organization - Institutions and companies
- #database - External databases
- #aging-biology - Aging-related content
- #cellular-senescence - Senescence research

**By Status**:
- status:active - Currently tracking
- status:pending - To be added
- priority:high - Important resources

---

## 📝 Templates

- [[../templates/researcher-profile|Researcher Profile Template]]
- [[../templates/organization-profile|Organization Profile Template]]
- [[../templates/database-tool-profile|Database/Tool Profile Template]]
- [[../templates/project-tracker|Project Tracker Template]]

---

## 🔗 External Links

- [Human Ageing Genomic Resources (HAGR)](https://genomics.senescence.info/)
- [CellAge Database](https://genomics.senescence.info/cells/)
- [BioGRID](https://thebiogrid.org/)
- [HuggingFace Hub](https://huggingface.co/)
- [Monarch Initiative](https://monarchinitiative.org/)

---

**Last Updated**: 2025-12-03
**Total Files**: 27 markdown files
**Network Connections**: 150+ wiki-links
**Next Milestone**: HuggingFace data publication

[[HOME|← Back to Projects Hub]]
