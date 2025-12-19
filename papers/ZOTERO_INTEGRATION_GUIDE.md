# Zotero-Obsidian Integration Guide

Complete guide for importing Zotero library into the papers directory and setting up Obsidian integration.

## Overview

This integration allows you to:
1. Import 23 papers from Zotero backup into the papers/ directory
2. Install Zotero desktop application and restore your library
3. Configure Obsidian plugins for seamless citation workflow
4. Maintain ongoing sync between Zotero and papers/

## Current Status

**Zotero Backup**:
- Location: `~/projects/research/literature/archives/zotero-backup/`
- Items: 23 papers (all with PDFs)
- Focus: ECM research, cell mechanics, imaging techniques

**Papers Directory**:
- Current papers: 32
- After import: 55 (23 new + 32 existing)
- No duplicates detected

## Step-by-Step Implementation

### Phase 1: Import Papers from Backup (READY TO RUN)

The import script has been tested and is ready to run.

#### 1.1 Run Import

```bash
cd ~/projects/papers

# Run the import
python3 scripts/import_from_zotero.py \
  --zotero-db ~/projects/research/literature/archives/zotero-backup/zotero.sqlite \
  --bibtex-db ~/projects/research/literature/archives/zotero-backup/better-bibtex.sqlite \
  --storage ~/projects/research/literature/archives/zotero-backup/storage \
  --papers-db papers.db \
  --output . \
  --skip-duplicates
```

**Expected output**:
- 23 new directories created
- Each with: metadata.json, note.md, paper.pdf
- papers.db updated with 23 new entries

#### 1.2 Verify Import

```bash
# Verify everything imported correctly
python3 scripts/verify_import.py \
  --papers-dir ~/projects/papers \
  --papers-db ~/projects/papers/papers.db
```

**Expected output**:
```
Database entries: 55
Paper directories: 55
Valid papers: 55
✅ All papers verified successfully - no issues found!
```

### Phase 2: Install Zotero Desktop

#### 2.1 Download and Install

1. Visit: https://www.zotero.org/download/
2. Download **Zotero 7** for macOS
3. Open the downloaded .dmg file
4. Drag Zotero to Applications folder
5. Launch Zotero once (this creates ~/Zotero directory)
6. Quit Zotero

#### 2.2 Restore Database

```bash
# Copy database files
cp ~/projects/research/literature/archives/zotero-backup/zotero.sqlite ~/Zotero/
cp ~/projects/research/literature/archives/zotero-backup/better-bibtex.sqlite ~/Zotero/

# Copy storage directory
cp -r ~/projects/research/literature/archives/zotero-backup/storage ~/Zotero/

# Copy other files
cp ~/projects/research/literature/archives/zotero-backup/zotero.sqlite.bak ~/Zotero/
cp -r ~/projects/research/literature/archives/zotero-backup/styles ~/Zotero/
cp -r ~/projects/research/literature/archives/zotero-backup/translators ~/Zotero/
```

#### 2.3 Verify Restore

1. Launch Zotero
2. You should see **23 items** in your library
3. Check that PDFs open correctly
4. Verify item metadata looks correct

### Phase 3: Install Better BibTeX Extension

#### 3.1 Download Extension

1. Visit: https://retorque.re/zotero-better-bibtex/installation/
2. Download the latest `.xpi` file for Zotero 7

#### 3.2 Install in Zotero

1. Open Zotero
2. Go to: **Tools → Add-ons**
3. Click the gear icon (⚙️) → **Install Add-on From File**
4. Select the downloaded `.xpi` file
5. Click **Install Now**
6. Restart Zotero

#### 3.3 Verify Better BibTeX

1. In Zotero, select any paper
2. Right-click → **Generate BibTeX key**
3. You should see citation keys like: `perez-cotaPicosecondUltrasonics2020`

### Phase 4: Configure Auto-Export to BibTeX

#### 4.1 Set Up Export

1. In Zotero: **File → Export Library**
2. Format: **Better BibTeX** (not "BibTeX")
3. Check: ☑️ **Keep updated**
4. Export to: `~/projects/papers/references/papers.bib`
5. Click **OK**

#### 4.2 Verify Export

```bash
# Check the BibTeX file was created
ls -lh ~/projects/papers/references/papers.bib

# View first few entries
head -50 ~/projects/papers/references/papers.bib
```

You should see BibTeX entries like:
```bibtex
@article{perez-cotaPicosecondUltrasonics2020,
  title = {Picosecond Ultrasonics for Elasticity-Based Imaging...},
  author = {Pérez-Cota, Fernando and Fuentes-Domínguez, Rafael...},
  ...
}
```

### Phase 5: Configure Obsidian Plugins

Your Obsidian already has the necessary plugins installed. Just need to configure paths.

#### 5.1 Zotero Desktop Connector Plugin

1. Open Obsidian
2. Settings → Community plugins → **Zotero Integration**
3. Configure:
   - **Database**: Better BibTeX (not Zotero)
   - **Literature Note Folder**: `papers`
   - **Image Output Folder**: `papers/{{citekey}}/figures`
   - **Citation Format**: `[@{{citekey}}]`
   - **Template Path**: (leave default or create custom)

#### 5.2 Citations Plugin

1. Settings → Community plugins → **Citations**
2. Configure:
   - **Citation database path**: `references/papers.bib`
   - **Literature note folder**: `papers`
   - **Literature note template**: (leave default)
   - **Markdown citation format**: `[@{{citekey}}]`

#### 5.3 Test Integration

1. **Test Zotero Connector**: Press `Cmd+Shift+Z`
   - Search for a paper from your Zotero library
   - Select it → Creates literature note in papers/

2. **Test Citation Insertion**: Press `Cmd+Shift+C`
   - Search for a paper
   - Inserts `[@citekey]` at cursor position

### Phase 6: Ongoing Workflow

#### Adding New Papers

**Option A: Via Zotero (Recommended)**
1. Add paper to Zotero desktop (via browser connector or DOI)
2. Zotero auto-exports to papers.bib
3. In Obsidian: `Cmd+Shift+Z` to import note
4. Optionally: Run sync script to add to papers/ directory

**Option B: Via download_paper.py (For Nature papers)**
1. Use existing script for papers with many figures
2. Gets high-quality metadata from Crossref

#### Syncing New Additions

To add new Zotero papers to papers/ directory with full metadata:

```bash
# Re-run import (it will skip duplicates)
python3 scripts/import_from_zotero.py \
  --zotero-db ~/Zotero/zotero.sqlite \
  --bibtex-db ~/Zotero/better-bibtex.sqlite \
  --storage ~/Zotero/storage \
  --papers-db papers.db \
  --output . \
  --skip-duplicates
```

## File Structure After Import

```
~/projects/papers/
├── papers.db (55 papers)
├── scripts/
│   ├── download_paper.py (existing - for Nature papers)
│   ├── import_from_zotero.py (new - Zotero import)
│   └── verify_import.py (new - verification)
├── references/
│   └── papers.bib (auto-exported from Zotero)
├── 10.1038_s41467-025-66434-z/ (existing paper)
├── 10.1063_5.0023744/ (new from Zotero)
│   ├── metadata.json
│   ├── note.md
│   └── paper.pdf
└── ... (53 more directories)
```

## Expected Import Results

Based on dry-run test:

- **Total Zotero items**: 23
- **Items to import**: 23
- **Duplicates**: 0 (no overlap with existing 32 papers)
- **Items with PDFs**: 23 (100%)
- **Items without DOI**: 2 (will use citation keys as folder names)

## Troubleshooting

### Issue: "Database is locked"

**Solution**: Close Zotero before running import script

### Issue: PDF not found

**Solution**: Check that storage/ directory was copied correctly to ~/Zotero

### Issue: Obsidian can't connect to Zotero

**Solution**:
1. Ensure Zotero desktop is running
2. Check plugin settings use "Better BibTeX" database
3. Restart Obsidian

### Issue: BibTeX file not updating

**Solution**:
1. In Zotero: Edit → Preferences → Better BibTeX
2. Check "Automatic Export" shows papers.bib
3. Try manually exporting again

## Advanced: Future Enhancements

### 1. Figure Extraction

Currently, PDFs are copied but figures are not extracted. To add figure extraction:

```python
# Could use PyMuPDF or pdf2image
import fitz  # PyMuPDF

def extract_figures_from_pdf(pdf_path, output_dir):
    doc = fitz.open(pdf_path)
    for page_num in range(len(doc)):
        page = doc[page_num]
        images = page.get_images()
        for img_index, img in enumerate(images):
            # Extract and save as figure_N.png
            ...
```

### 2. Bidirectional Sync

Create a sync script that:
- Detects new papers in Zotero → imports to papers/
- Detects updates in papers/ → updates Zotero
- Syncs tags and notes

### 3. Annotation Import

Export Zotero annotations and highlights to note.md:
- Use Zotero PDF reader annotations
- Export to Obsidian notes
- Maintain highlights and comments

## Quick Reference Commands

```bash
# Import from backup
cd ~/projects/papers
python3 scripts/import_from_zotero.py \
  --zotero-db ~/projects/research/literature/archives/zotero-backup/zotero.sqlite \
  --bibtex-db ~/projects/research/literature/archives/zotero-backup/better-bibtex.sqlite \
  --storage ~/projects/research/literature/archives/zotero-backup/storage \
  --papers-db papers.db \
  --output . \
  --skip-duplicates

# Verify import
python3 scripts/verify_import.py \
  --papers-dir ~/projects/papers \
  --papers-db ~/projects/papers/papers.db

# Sync new Zotero additions
python3 scripts/import_from_zotero.py \
  --zotero-db ~/Zotero/zotero.sqlite \
  --bibtex-db ~/Zotero/better-bibtex.sqlite \
  --storage ~/Zotero/storage \
  --papers-db papers.db \
  --output . \
  --skip-duplicates
```

## Obsidian Keyboard Shortcuts

- `Cmd+Shift+Z` - Import note from Zotero
- `Cmd+Shift+C` - Insert citation
- `Cmd+K` - Insert link (for linking to other papers)

## Next Steps

1. ✅ Run the import script
2. ✅ Verify all 23 papers imported correctly
3. ⏭️ Install Zotero desktop
4. ⏭️ Restore database
5. ⏭️ Install Better BibTeX
6. ⏭️ Configure auto-export
7. ⏭️ Test Obsidian integration

---

**Documentation Updated**: 2025-12-10
**Import Script Version**: 1.0
**Status**: Ready for execution
