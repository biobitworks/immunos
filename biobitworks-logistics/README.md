# Biobitworks Logistics

## Summary
Contains 2 subdirectories and 1 files. Key subfolders: labs/, protocols/.


Community lab network coordination and logistics management for the San Francisco Bay Area.

## Overview

Biobitworks Logistics is a project to create and coordinate a network of community biology labs across the Bay Area, fostering collaboration, resource sharing, and accessibility to biotechnology education and research.

## Network Labs

### Active Labs
- **Counter Culture Labs** (Oakland) - Existing community lab
- **Biopunk Lab** (San Francisco) - Existing community lab
- **BioCurious** (Santa Clara) - Established community lab

### In Development
- **Oakland BioLab** - New lab in planning stages
- **Berkeley BioLab** - Potential future location

## Project Structure

```
biobitworks-logistics/
├── labs/                    # Lab database and profiles
│   ├── labs.json           # Structured lab data
│   └── profiles/           # Detailed lab profiles (future)
├── protocols/              # Shared protocols
│   ├── protocols.json      # Protocol registry
│   ├── core/              # Essential cross-lab protocols
│   ├── shared/            # Shared between specific labs
│   └── lab-specific/      # Lab-unique protocols
├── networks/              # Network coordination
│   ├── bay-area/         # Regional network info
│   └── partnerships/     # Collaborations and partnerships
├── resources/            # Shared resources
│   ├── equipment/        # Equipment databases and sharing
│   ├── funding/          # Grants and funding opportunities
│   └── legal/            # Legal docs, insurance, etc.
├── docs/                # Documentation
└── data/                # Data and analytics
```

## Key Features

### Lab Database
- Structured JSON database of all network labs
- Contact information, equipment, specializations
- Membership and access details
- Status tracking (active, planned, potential)

### Protocol Library
- Core protocols used across all labs
- Shared protocols for specific collaborations
- Lab-specific protocols
- Hybrid sync with Google Drive

### Network Coordination
- Inter-lab communication and collaboration
- Resource sharing and equipment access
- Joint events and workshops
- Collective purchasing and funding

## Getting Started

### View Lab Network

```bash
# View all labs
cat labs/labs.json | jq '.labs[] | {name, city: .location.city, status}'

# View active labs only
cat labs/labs.json | jq '.labs[] | select(.status == "active")'
```

### Protocol Management

```bash
# List all protocols
cat protocols/protocols.json | jq '.protocols[] | {id, title, category}'

# View protocols to download from Google Drive
cat protocols/protocols.json | jq '.protocols[] | select(.status == "to_download")'
```

## Next Steps

1. **Populate Lab Data**: Add detailed information for each lab
   - Contact details
   - Equipment lists
   - Specializations and focus areas
   - Membership information

2. **Protocol Sync**: Download protocols from Google Drive
   - Identify core protocols to download
   - Set up sync mechanism
   - Organize by category

3. **Network Planning**: Develop coordination strategy
   - Communication channels
   - Resource sharing policies
   - Event coordination
   - Collective initiatives

4. **Oakland Lab Development**: Plan new lab
   - Location scouting
   - Equipment needs
   - Funding strategy
   - Community engagement

## Google Drive Integration

**Protocol Folder**: [Add Google Drive folder URL]

**Sync Strategy**:
- Core protocols: Download and version locally
- Shared protocols: Link to Drive, selective download
- Periodic sync to check for updates

## Resources

- [Counter Culture Labs Website](https://counterculturelabs.org)
- [BioCurious Website](https://biocurious.org)
- [DIYbio.org](https://diybio.org)
- [Biosafety in Microbiological and Biomedical Laboratories (CDC)](https://www.cdc.gov/labs/pdf/SF__19_308133-A_BMBL6_00-BOOK-WEB-final-3.pdf)

## Related Projects

- **BioViz Research**: [[../research/projects/bioviztech/README|BioViz Project]]
- **Research Network**: [[../research/INDEX|Research Dashboard]]

## Contact

**Project Lead**: Byron
**Updated**: 2025-12-10

---

**Status**: Initial setup complete
**Next Action**: Add Google Drive folder URL and begin protocol sync

## Directory Map
```
biobitworks-logistics/
├── labs/
├── protocols/
└── GOOGLE_DRIVE_SYNC.md
```
