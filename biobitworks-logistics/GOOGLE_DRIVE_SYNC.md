# Google Drive Protocol Sync

Guide for syncing protocols from Google Drive to the local Biobitworks Logistics project.

## Google Drive Setup

**Protocol Folder URL**: [Add your Google Drive folder URL here]

**Access**: Ensure the folder is shared with appropriate permissions

## Sync Strategy

### Core Protocols
- **Action**: Download and maintain locally
- **Location**: `protocols/core/`
- **Reason**: Critical protocols should be available offline
- **Update Frequency**: Manual review for changes

### Shared Protocols
- **Action**: Selective download based on active collaborations
- **Location**: `protocols/shared/`
- **Reason**: Only download what's actively being used
- **Update Frequency**: As needed

### Lab-Specific Protocols
- **Action**: Link to Drive, download on request
- **Location**: `protocols/lab-specific/`
- **Reason**: Each lab may have unique variations
- **Update Frequency**: As requested by labs

## Manual Download Process

### Option 1: Using Google Drive Desktop

1. Install Google Drive for Desktop if not already installed
2. Add the shared folder to "My Drive"
3. Make the folder "Available offline"
4. Copy files from the synced folder to project:
   ```bash
   cp ~/GoogleDrive/BiobitworksProtocols/* ~/projects/biobitworks-logistics/protocols/
   ```

### Option 2: Direct Download

1. Open the Google Drive folder in browser
2. For each protocol:
   - Download the file
   - Move to appropriate category folder
   - Update `protocols.json` with local path

### Option 3: Using `gdown` (Python)

Install gdown:
```bash
pip install gdown
```

Download individual file:
```bash
# Get the file ID from the Google Drive URL
# URL format: https://drive.google.com/file/d/FILE_ID/view
gdown https://drive.google.com/uc?id=FILE_ID -O protocols/core/filename.md
```

Download entire folder:
```bash
# Get folder ID from URL
# URL format: https://drive.google.com/drive/folders/FOLDER_ID
gdown --folder https://drive.google.com/drive/folders/FOLDER_ID -O protocols/downloads/
```

## After Download: Update Protocol Registry

After downloading protocols, update `protocols/protocols.json`:

```json
{
  "id": "protocol-XXX",
  "title": "Protocol Title",
  "category": "core|shared|lab_specific",
  "status": "complete",
  "source": "google_drive",
  "google_drive_url": "https://drive.google.com/...",
  "local_path": "protocols/core/filename.md",
  "applicable_labs": ["all"],
  "version": "1.0",
  "last_updated": "2025-12-10",
  "tags": ["keyword1", "keyword2"]
}
```

## Sync Script (Future Enhancement)

Create an automated sync script:

```python
#!/usr/bin/env python3
"""
Sync protocols from Google Drive

Usage:
    python sync_protocols.py --folder-id FOLDER_ID
"""

import gdown
import json
from pathlib import Path
from datetime import datetime

def sync_protocols(folder_id, output_dir):
    # Download from Google Drive
    gdown.download_folder(
        f"https://drive.google.com/drive/folders/{folder_id}",
        output=str(output_dir),
        quiet=False
    )

    # Update protocols.json
    # ... implementation ...
```

## Tracking Protocol Updates

### Manual Tracking

Keep a log of when protocols were last synced:

```bash
# Update sync log
echo "$(date): Synced protocols from Google Drive" >> protocols/sync_log.txt
```

### Using protocols.json

Update the `google_drive` section:

```json
{
  "google_drive": {
    "folder_url": "https://drive.google.com/drive/folders/YOUR_FOLDER_ID",
    "sync_status": "synced",
    "last_sync": "2025-12-10T19:00:00Z"
  }
}
```

## Best Practices

1. **Version Control**: Use git to track changes to downloaded protocols
   ```bash
   cd ~/projects/biobitworks-logistics
   git add protocols/
   git commit -m "Synced protocols from Google Drive"
   ```

2. **Change Detection**: Before syncing, check Google Drive folder's "Last modified" date

3. **Backup**: Keep a backup of downloaded protocols before syncing new versions

4. **Documentation**: Update protocol version numbers when changes are made

5. **Notification**: Notify lab network when critical protocols are updated

## Folder Structure on Google Drive

Recommended organization in Google Drive:

```
Biobitworks Protocols/
├── Core/
│   ├── Safety_BSL1.pdf
│   ├── Agar_Plates.pdf
│   └── Equipment_Maintenance.pdf
├── Shared/
│   ├── PCR_Protocol.pdf
│   └── Gel_Electrophoresis.pdf
└── Lab_Specific/
    ├── CounterCultureLabs/
    ├── BiopunkLab/
    └── BioCurious/
```

## Troubleshooting

### Issue: Permission denied

**Solution**: Ensure you have view/download access to the Google Drive folder

### Issue: Files not downloading

**Solution**:
- Check internet connection
- Verify folder ID is correct
- Check file size limits

### Issue: Duplicate files

**Solution**: Check `protocols.json` before downloading to avoid duplicates

## Next Steps

1. Add Google Drive folder URL to `protocols/protocols.json`
2. Download core protocols from Drive
3. Set up git repository for version control
4. Create sync schedule (monthly review)
5. Notify lab network of protocol updates

---

**Last Updated**: 2025-12-10
**Maintained by**: Byron
