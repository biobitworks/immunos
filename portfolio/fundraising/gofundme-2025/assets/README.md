# Assets & Analysis Folder

This directory stores all media files, documents, and data for Claude to analyze and incorporate into the GoFundMe campaign.

## Folder Structure

```
assets/
├── photos/              # Campaign photos, event images, community pictures
├── videos/              # Raw video footage, video scripts, final edits
├── documents/           # PDFs, financial docs, contracts, invoices
├── data/                # Spreadsheets, donor lists, network analysis
└── analysis/            # Claude's analysis reports of uploaded files
```

---

## 📸 Photos Folder (`photos/`)

**Upload photos here for Claude to analyze and select the best campaign image**

### What to Include:
- **Community event photos:**
  - You cooking at Viva Frontier Tower
  - Longevity meetup photos at 1900 Broadway
  - Vitalist Bay popup photos
  - Group shots with community members

- **Personal photos:**
  - Professional headshot (face visible, approachable)
  - Working on HealthClock or bioviz.tech
  - At your current housing location

- **Venue photos:**
  - 1900 Broadway event spaces
  - Community kitchen/dining areas
  - Lab or research workspace

### What Claude Will Analyze:
- ✅ Image quality and lighting
- ✅ Emotional connection (authentic, warm, engaging)
- ✅ GoFundMe best practices (faces perform better than objects)
- ✅ Story alignment (does it match your narrative?)
- ✅ Technical requirements (resolution, aspect ratio)

### Recommended Photo Format:
- **Resolution:** At least 1200x800px
- **Format:** JPG, PNG
- **Aspect ratio:** 16:9 or 4:3
- **File size:** Under 10MB

---

## 🎥 Videos Folder (`videos/`)

**Upload video files or scripts here for Claude to review**

### What to Include:
- **Raw video footage:**
  - 90-second pitch recording (multiple takes)
  - B-roll footage (community events, workspace, etc.)
  - Testimonial clips from community members (if available)

- **Video scripts:**
  - Draft scripts for different video lengths (60s, 90s, 2min)
  - Talking points outline
  - Call-to-action variations

### What Claude Will Analyze:
- ✅ Script effectiveness (hook, story, urgency, CTA)
- ✅ Video quality assessment (lighting, audio, framing)
- ✅ Emotional authenticity
- ✅ Length optimization (GoFundMe sweet spot is 2-3 minutes, but 90s is better for social)
- ✅ Platform-specific optimization (vertical for Instagram, horizontal for LinkedIn)

### Recommended Video Format:
- **Length:** 90 seconds (optimal for attention span)
- **Resolution:** 1080p minimum
- **Format:** MP4, MOV
- **Audio:** Clear voice audio, minimal background noise
- **Captions:** Add text overlay or upload SRT file

---

## 📄 Documents Folder (`documents/`)

**Upload financial docs, contracts, invoices, medical records for Claude to analyze**

### What to Include:
- **Housing documents:**
  - Lease agreement (current and new)
  - Eviction notice or deposit demand letter
  - Landlord communication (emails, texts)
  - Rental payment history

- **Financial records:**
  - SSDI award letter (redact SSN if needed)
  - Contract work invoices (outstanding and paid)
  - Bank statements (last 3 months, redact account numbers)
  - Student loan statements (servicer name, balance, status)
  - Credit card statements (balances, minimum payments)

- **Medical documentation:**
  - Hospital discharge summary (July cellulitis)
  - Medical bills from hospitalization
  - Insurance EOB (Explanation of Benefits)
  - Doctor's notes on recovery timeline

- **Grant applications:**
  - Foresight Fast Grant application (if shareable)
  - FRO proposal excerpts
  - Any pending funding applications

### What Claude Will Analyze:
- ✅ **Financial forensics:** Exact income, expenses, debt breakdown
- ✅ **Legal opportunities:** SSDI TPD discharge eligibility, tenant rights, tax optimization
- ✅ **Credibility verification:** Document authenticity for campaign transparency
- ✅ **Emergency resource matching:** Which Oakland programs you qualify for based on income
- ✅ **Wealth strategy:** Income stabilization opportunities, debt prioritization

### Privacy & Security:
- **Redact sensitive info:** SSN, full account numbers, personal health details
- Claude can analyze documents with partial redaction
- Files stay local on your machine - not shared with external services

---

## 📊 Data Folder (`data/`)

**Upload spreadsheets, contact lists, network analysis for Claude to optimize outreach**

### What to Include:
- **Network mapping:**
  - Viva Frontier Tower contact list (names, donation likelihood, relationship strength)
  - Vitalist Bay attendee list
  - 1900 Broadway longevity meetup attendees
  - LinkedIn connections (export CSV)
  - Telegram group member counts

- **Donor tracking:**
  - Past fundraising results (if any)
  - Pledge commitments (pre-campaign asks)
  - Contact response tracking (who you've DM'd, responses)

- **Financial modeling:**
  - Monthly budget breakdown (income vs expenses)
  - Debt payment schedule
  - Contract work hourly rate analysis
  - Grant funding timeline projection

- **Campaign performance:**
  - GoFundMe analytics (once live)
  - Social media engagement metrics
  - Referral source tracking

### What Claude Will Analyze:
- ✅ **Network size assessment:** Realistic fundraising potential based on reach
- ✅ **Contact prioritization:** Who to DM first (relationship strength + donation capacity)
- ✅ **Financial optimization:** Budget cuts, income maximization, debt payoff strategy
- ✅ **Campaign performance:** Which channels are converting, where to double down
- ✅ **Predictive modeling:** Likelihood of hitting $2,400 vs $4,200 goal

### Recommended Data Format:
- **Spreadsheets:** CSV, XLSX, Google Sheets (export as CSV)
- **Contact lists:** Name, relationship, contact method, estimated donation capacity, last contact date
- **Financial data:** Monthly breakdown, line-item detail

---

## 🔍 Analysis Folder (`analysis/`)

**Claude's generated reports will be saved here automatically**

### What Gets Created:
- **Photo analysis reports:**
  - `photo-selection-report.md` - Best campaign images ranked by effectiveness
  - `image-optimization-guide.md` - Cropping, lighting, and editing recommendations

- **Video analysis reports:**
  - `video-script-analysis.md` - Script effectiveness breakdown
  - `video-editing-suggestions.md` - Cut recommendations, pacing, B-roll needs

- **Document analysis reports:**
  - `financial-forensics-report.md` - Complete income/expense/debt breakdown
  - `legal-opportunities-report.md` - SSDI discharge, tenant rights, tax optimization
  - `emergency-resources-match.md` - Oakland programs you qualify for with application steps

- **Network analysis reports:**
  - `contact-prioritization-matrix.md` - Who to DM first, expected donation amounts
  - `fundraising-projection-model.md` - Realistic goal assessment based on network size
  - `campaign-performance-analysis.md` - Real-time tracking of what's working

---

## How to Use This Folder

### 1. Upload Files
Drop your files into the appropriate subfolder:
```bash
# Example: Upload a photo
cp ~/Pictures/viva-cooking-event.jpg ~/gofundme-campaign-2025/assets/photos/

# Example: Upload financial docs
cp ~/Documents/ssdi-award-letter.pdf ~/gofundme-campaign-2025/assets/documents/

# Example: Upload contact list
cp ~/Downloads/linkedin-connections.csv ~/gofundme-campaign-2025/assets/data/
```

### 2. Ask Claude to Analyze
Tell Claude which files to analyze:

**Examples:**
- "Analyze all photos in assets/photos/ and recommend the best 3 for the GoFundMe campaign image"
- "Review my SSDI award letter in assets/documents/ and tell me if I qualify for TPD student loan discharge"
- "Look at my network contact list in assets/data/ and create a prioritized outreach plan"
- "Analyze my bank statements and create a monthly budget optimization plan"

### 3. Review Claude's Analysis
Claude will create detailed reports in `assets/analysis/` with:
- ✅ Specific recommendations
- ✅ Action items with priority levels
- ✅ Data-driven insights (e.g., "Based on 150 LinkedIn connections, you can realistically raise $1,200-1,800")
- ✅ Step-by-step implementation guides

### 4. Implement Recommendations
Use Claude's analysis to:
- Select the best campaign photo
- Refine video script
- Optimize outreach strategy
- Apply for emergency assistance programs
- Negotiate with landlord using financial data
- Set realistic fundraising goal

---

## Privacy & Confidentiality

**Your files stay local:**
- All analysis happens on your machine
- Claude doesn't upload files to external servers
- You control what gets shared in the final campaign

**Redaction Guidelines:**
- **Always redact:** Social Security Numbers, full account numbers
- **Optional redact:** Full address, phone numbers, personal health details
- **Keep visible:** Amounts, dates, names of institutions, timelines

**What Claude Needs to See:**
- ✅ Dollar amounts (income, expenses, debts)
- ✅ Dates (timelines, deadlines)
- ✅ Institutions (SSDI, loan servicers, landlord)
- ✅ Relationships (who you know, how well)

---

## Quick Start Examples

### Example 1: Photo Selection
```bash
# Upload 5-10 photos
cp ~/Pictures/community-events/* ~/gofundme-campaign-2025/assets/photos/

# Ask Claude:
"Analyze all photos in assets/photos/ and rank them for GoFundMe effectiveness.
Consider: emotional connection, story alignment, image quality, faces visible."
```

### Example 2: Financial Forensics
```bash
# Upload SSDI letter and bank statements
cp ~/Documents/ssdi-award.pdf ~/gofundme-campaign-2025/assets/documents/
cp ~/Documents/bank-statement-sept.pdf ~/gofundme-campaign-2025/assets/documents/

# Ask Claude:
"Review my SSDI award letter and bank statements. Calculate:
1. Exact monthly income (SSDI + contract work average)
2. Essential expenses breakdown
3. How much I can realistically save per month
4. Whether I qualify for TPD student loan discharge
5. Which Oakland emergency programs I'm eligible for"
```

### Example 3: Network Analysis
```bash
# Export LinkedIn connections to CSV
# Upload contact list
cp ~/Downloads/linkedin-connections.csv ~/gofundme-campaign-2025/assets/data/

# Ask Claude:
"Analyze my LinkedIn network (assets/data/linkedin-connections.csv).
Create a prioritized outreach plan:
1. Segment by relationship strength (close friends, colleagues, acquaintances)
2. Estimate realistic donation amounts per segment
3. Calculate total fundraising potential
4. Generate personalized DM templates for top 20 contacts
5. Recommend optimal contact timing and sequence"
```

### Example 4: Video Script Review
```bash
# Upload your draft video script
cp ~/Documents/gofundme-video-script.txt ~/gofundme-campaign-2025/assets/videos/

# Ask Claude:
"Review my video script (assets/videos/gofundme-video-script.txt).
Analyze for:
1. Hook effectiveness (first 5 seconds)
2. Emotional authenticity
3. Story clarity (problem → solution → ask)
4. Call-to-action strength
5. Length optimization (is it under 90 seconds?)
Provide specific rewrite suggestions."
```

---

## File Naming Conventions

**Use descriptive, dated filenames:**

**Photos:**
- `2025-09-15-viva-cooking-dinner.jpg`
- `2025-08-10-1900-broadway-meetup.jpg`
- `headshot-professional-2025.jpg`

**Documents:**
- `ssdi-award-letter-2023.pdf`
- `bank-statement-sept-2025.pdf`
- `lease-agreement-current.pdf`
- `hospital-bill-july-cellulitis.pdf`

**Data:**
- `linkedin-connections-export-2025-10-02.csv`
- `viva-contact-list-donation-likelihood.xlsx`
- `monthly-budget-sept-2025.csv`

**Videos:**
- `gofundme-pitch-take1-90sec.mp4`
- `gofundme-pitch-take2-revised.mp4`
- `community-broll-viva-events.mov`

---

## Next Steps

1. **Upload your first batch of files** (start with photos and any financial docs you have)
2. **Ask Claude to analyze** them using the example prompts above
3. **Review Claude's reports** in `assets/analysis/`
4. **Implement recommendations** to optimize your campaign
5. **Iterate:** Upload more files as you gather them, refine based on Claude's feedback

**Ready to start?** Upload files and tell Claude what you need analyzed!
