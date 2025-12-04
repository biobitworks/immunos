# Obsidian Workspace Configuration Guide

This guide explains how to set up separate workspaces for each project in the Byron Projects vault.

## What Are Workspaces?

Obsidian Workspaces let you save different layouts (open files, sidebar configurations, graph filters) and switch between them instantly. Perfect for separating different projects in one vault!

## Setup Instructions

### 1. Enable Workspaces Core Plugin

1. Open Obsidian Settings (`Ctrl/Cmd + ,`)
2. Go to **Core Plugins**
3. Enable **Workspaces** plugin
4. Close settings

### 2. Create Project Workspaces

#### Workspace 1: Aging Research

**Setup**:
1. Open [[HOME|HOME.md]]
2. Click [[research/INDEX|Aging Research Dashboard]]
3. Open graph view (`Ctrl/Cmd + G`)
4. In graph view, click filters icon
5. Add filter: `path:research/`
6. This shows only research/* files in graph
7. Arrange your preferred layout (left: file explorer, right: graph)
8. Open command palette (`Ctrl/Cmd + P`)
9. Type "Workspace: Save current workspace"
10. Name it: **"Aging Research"**

**Result**: Workspace focused on aging biology files only.

#### Workspace 2: IMMUNOS-MCP

**Setup**:
1. Open [[immunos-mcp/INDEX|IMMUNOS-MCP Dashboard]]
2. Open graph view
3. Add filter: `path:immunos-mcp/`
4. Arrange layout
5. Save workspace as: **"IMMUNOS-MCP"**

#### Workspace 3: IMMUNOS81

**Setup**:
1. Open [[immunos81/INDEX|IMMUNOS81 Dashboard]]
2. Open graph view
3. Add filter: `path:immunos81/`
4. Arrange layout
5. Save workspace as: **"IMMUNOS81"**

#### Workspace 4: All Projects

**Setup**:
1. Open [[HOME|HOME.md]]
2. Open graph view with NO filters
3. Shows entire vault network
4. Save workspace as: **"All Projects"**

### 3. Switch Between Workspaces

**Method 1**: Command Palette
- Press `Ctrl/Cmd + P`
- Type "Workspace: Load"
- Select workspace from list

**Method 2**: Workspace Switcher
- Press `Ctrl/Cmd + L` (Windows/Mac)
- Select workspace

**Method 3**: Ribbon Icon
- Click workspace icon in left ribbon
- Select workspace

## Graph Filters Reference

Use these filters in graph view to isolate projects:

```
path:research/          # Aging biology research only
path:immunos-mcp/       # IMMUNOS-MCP only
path:immunos81/         # IMMUNOS81 only
path:research/researchers/  # Just researcher profiles
path:research/databases/    # Just database docs
tag:#aging-biology      # By tag
```

## Recommended Workspace Layouts

### Aging Research Workspace

**Left Sidebar**:
- File Explorer (filtered to research/)
- Search
- Tags

**Main Area**:
- research/INDEX.md (pinned)
- Other working notes

**Right Sidebar**:
- Graph view (filtered to research/)
- Backlinks
- Outline

### IMMUNOS Workspaces

**Left Sidebar**:
- File Explorer
- Search

**Main Area**:
- INDEX.md (pinned)
- Code documentation
- API docs

**Right Sidebar**:
- Graph view (filtered to project)
- Outline

## Advanced Tips

### 1. Auto-Load Workspace

Settings → Workspaces → "Startup workspace"
- Choose which workspace opens when you start Obsidian
- Recommendation: "All Projects" or "Aging Research"

### 2. Workspace-Specific Settings

Each workspace remembers:
- Open files and their positions
- Sidebar visibility
- Graph filters and layout
- Pane sizes

### 3. Folder Colors (Community Plugin)

Install "Folder Note" or "Colorful Folders" plugin to color-code projects:
- 🧬 Green: research/
- 🛡️ Blue: immunos-mcp/
- 🔬 Purple: immunos81/

### 4. Hotkeys

Set custom hotkeys for workspaces:
1. Settings → Hotkeys
2. Search "Workspace"
3. Assign keys like:
   - `Ctrl/Cmd + 1`: Aging Research
   - `Ctrl/Cmd + 2`: IMMUNOS-MCP
   - `Ctrl/Cmd + 3`: IMMUNOS81

## Network Graph Visualization

### What Shows in Graph

With filter `path:research/`, the network will show:

**Nodes** (files):
- Researcher profiles (Peter Lidsky)
- Organization profiles (CityUHK, Monarch, etc.)
- Database docs (HuggingFace, spatialLIBD, BioGRID)
- Literature notes (BCL6 paper)
- Templates
- READMEs

**Edges** (links):
- Lidsky ↔ CityUHK (affiliation)
- Lidsky ↔ BCL6 paper (research interest)
- spatialLIBD ↔ GenAge (integration)
- BioGRID ↔ CellAge (network analysis)
- Templates ↔ Profiles (documentation structure)

**Visual**: Each project forms a cluster, connections between projects appear when linked.

### Graph View Settings

Recommended settings for clean graph:

**Filters**:
- ✅ Existing files only
- ❌ Orphans (unless debugging)
- ✅ Tags

**Display**:
- Node size: 5-8
- Link thickness: 1-2
- Repel force: 200-500
- Link distance: 150-250
- Center force: 0.3-0.5

**Groups**:
- Group by folder or tags
- Color-code by project

## Troubleshooting

### Graph Shows Too Many Files

**Solution**: Apply path filter
- In graph view, click filter icon
- Add: `path:research/` for aging research only

### Links Not Working

**Issue**: Wiki-links need proper format
**Fix**: Use `[[folder/file|Display Text]]` format

### Workspace Not Saving

**Issue**: Plugin not enabled
**Fix**: Settings → Core Plugins → Enable "Workspaces"

### Cross-Project Links

**Intentional**: Some files link between projects
- Research ↔ IMMUNOS (algorithm papers)
- This is fine! Use "All Projects" workspace to see full network

## Maintenance

### Weekly

- Review orphaned notes (files with no links)
- Update INDEX dashboards with new content
- Check that workspaces load correctly

### Monthly

- Clean up unused tags
- Consolidate duplicate notes
- Update graph filters if project structure changes

## Export and Backup

Workspaces are saved in:
```
/Users/byron/projects/.obsidian/workspaces.json
```

Backup this file to preserve your workspace configurations.

---

**Quick Reference**:
- Switch workspace: `Ctrl/Cmd + L`
- Open graph: `Ctrl/Cmd + G`
- Command palette: `Ctrl/Cmd + P`
- Quick search: `Ctrl/Cmd + O`

**Dashboard Links**:
- [[HOME|Projects Hub]]
- [[research/INDEX|Aging Research]]
- [[immunos-mcp/INDEX|IMMUNOS-MCP]]
- [[immunos81/INDEX|IMMUNOS81]]

---

**Created**: 2025-12-03
**Last Updated**: 2025-12-03
**Vault**: `/Users/byron/projects/`
