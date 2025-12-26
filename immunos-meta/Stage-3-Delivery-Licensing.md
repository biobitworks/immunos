# Stage 3 - Delivery + Licensing (Draft)

## A) immunOS Replication Experiment + Preprint Submission

### Required Artifacts
- Reproducible pipeline scripts (SciFact baseline: oracle + open)
- Run manifests (env + model hashes + dataset versions)
- Metrics JSON (oracle + open)
- Figures/tables in preprint (metrics table, pipeline diagram)
- Public repository or archive for replication bundle
- Reproducibility notes: OS/CPU, Python/conda versions, GPU (if any)

### Minimal Submission Checklist
1) Preprint text complete and proofread
2) Methods section contains full dataset and environment details
3) Exact command list for each pipeline stage
4) Results table with oracle and open metrics
5) Limitations section (dataset constraints, test labels unavailable)
6) Reproduction appendix (run order, time, compute requirements)

### Suggested Packaging
- `immunos-preprint/` as the canonical draft location
- `docs/reference/` as the authoritative log of experiments
- `data/immunos_data/research/scifact/` as data root
- `docs/reference/scifact-baseline-replication.md` as the single metrics source of truth

### Reproducibility Enhancements (Industry Standard)
- Create a single `run.sh` script for the pipeline
- Add `environment.yml` or `requirements.txt` snapshot for the legacy env
- Store model checksums and dataset download hashes
- Store a copy of the exact command transcript used in the final run

## B) immunos-mcp Productization (VM/USB Builder + Website)

### Product Concept
- A prebuilt, local-first AIS workbench that runs offline by default.
- Optionally includes MCP server integration for external tools.
- Distributed as VM image + USB builder to lower setup friction.

### Distribution Targets
- VM images (VirtualBox/VMware/QEMU)
- USB installer (persistent storage + local data)
- Optional hosted builder website to assemble custom images
- Optional container image (for CI and dev)

### Core Release Components
- CLI + orchestrator
- Default agents + minimal examples
- Local docs and onboarding
- Safe defaults and explicit opt-in for cloud/LLM providers
- Local-first telemetry settings (off by default)

### Optional Components (paid or add-on)
- Enterprise connectors (GitHub/GitLab, Jira, SIEM)
- Hosted syncing or telemetry (opt-in)
- Advanced datasets + domain packs (biomed, security, research)
- SLA support + compliance packs

## C) Licensing Strategy for ROI

### Options Summary
1) **Open-core + Commercial Add-ons**
   - Core under Apache-2.0 or MIT
   - Paid add-ons under commercial license
   - Pros: strong adoption, easy contributions
   - Cons: weaker moat unless add-ons are significant

2) **Dual License (AGPL + Commercial)**
   - Open source for community use
   - Commercial license for businesses needing proprietary integration
   - Pros: clear ROI path for enterprise
   - Cons: AGPL can deter some adopters

3) **Source-Available (BUSL / Elastic-style)**
   - Time-delayed open source or restricted use
   - Pros: strong ROI control
   - Cons: reduced open-source goodwill

### Recommendation (default)
- **Chosen**: Dual license for immunos-mcp (AGPL + commercial),
  with optional open-core add-ons if adoption needs increase.
- Preprint and replication materials remain permissive/open.
- Consider trademark policy to protect branding for hosted offerings.

## Open Questions
- Which add-ons should be paid vs open?
- Is a hosted builder website part of v1 or v2?
- Preferred distribution channel (direct download, marketplace)?
- Will binaries/VMs be distributed under the same license as code?

## D) VM/USB Builder MVP (v1 Scope)
- Goal: one-click local environment for offline use.
- Inputs: minimal config (models, datasets, update channel).
- Outputs: bootable USB + VM image + manifest.
- Constraints: avoid bundling restricted datasets.
- Security: checksum verification + signed manifests.
- MVP includes:
  - Base OS image + install script
  - Preloaded core app + default configs
  - Offline docs + quickstart
  - Verified build manifest (hashes, build date)
- Out of scope (v1):
  - User accounts
  - Telemetry dashboards
  - Paid add-on marketplace

## E) MVP Dependencies (v1)
- Build system for VM/USB images (Packer or custom scripts)
- Signing infrastructure for artifacts
- Storage for releases (GitHub Releases + checksum file)
- Minimal website landing page with downloads

## F) Website Architecture (v1)
- Landing page + docs + builder UI.
- Storage for build manifests (no user data by default).
- Download center for signed images.
- Optional account system only if telemetry or sync is enabled.

## G) Pricing Model (Draft)
- **Community**: free local-first core.
- **Pro**: paid add-ons (domain packs, advanced analytics, local compliance pack).
- **Enterprise**: commercial license + support SLA + connectors.
- Optional: per-seat vs per-site pricing; evaluate based on target org size.

## H) Open Source Path (Alternative)
- Fully open-source core under Apache-2.0 or MIT.
- Monetize via support, training, and hosted services.
- Keep optional proprietary packs separate (connectors, compliance tooling).

## I) Partner Channels
- GitHub Marketplace (tooling extensions)
- Research communities (preprint, OSF/Zenodo bundles)
- Security/AI tooling marketplaces (if enterprise connectors land)

## J) Compliance Targets (for enterprise)
- Data residency: local-first by default
- Logging/retention policy for audit trails
- SBOM generation and signed release artifacts
- Optional SOC2 path (later stage)

## K) Launch Timeline (Draft)
- Phase 0 (Week 0-1): finalize license + artifact list + build tooling choice.
- Phase 1 (Week 2-3): working VM image + offline docs + checksum/signing.
- Phase 2 (Week 4-5): USB installer + installer validation + release checklist.
- Phase 3 (Week 6): website landing page + download center.
- Phase 4 (Week 7+): add-ons + marketplace exploration.

## L) Future Features (Candidate Roadmap)

### Core Platform
- Policy engine for risk scoring + thresholds
- Pluggable feature extractors (code, text, biomedical)
- Local model registry + version pinning
- Agent test harness with golden datasets

### Productization
- VM/USB builder wizard (web)
- Update channel for offline deployments
- Marketplace for domain packs
- Optional hosted sync for teams

### Research Track
- Extended SciFact experiments (open retrieval variants)
- NegSl-AIS replication on additional datasets
- Benchmark suite across multiple AIS algorithms
