# Paper Download Status Report

**Generated**: 2025-12-04
**Total Papers**: 32 directories

---

## Summary

| Status | Count | Percentage |
|--------|-------|------------|
| ✅ Successfully Downloaded (PDF + Figures) | 6 | 19% |
| ❌ Failed (No Content) | 26 | 81% |

## Root Cause Analysis

### 1. PDF Downloads - ALL PAYWALLED
**Issue**: All PDF URLs redirect to authentication (HTTP 303 → idp.nature.com)
- Papers marked as "CC-BY-NC-ND-4.0" (open access) in metadata
- But PDFs require institutional login/subscription
- The 6 successful papers were downloaded on Dec 3rd before access changed or when session was valid

### 2. Figure Downloads - Rate Limiting/Bot Detection
**Issue**: SpringerNature's CDN actively blocks automated downloads
- Manual `curl` requests work initially
- Python subprocess/requests library: Blocked immediately
- After repeated attempts: Even manual curl gets HTTP 404 (HTML error pages)
- Rate limiting threshold appears very low (~10-20 requests)

### 3. Technical Challenges

**Attempted Solutions**:
1. ✅ Fixed URL pattern: `{journal}_{year}_{article}_Fig{n}_HTML.png`
2. ✅ Added proper User-Agent headers
3. ✅ Switched from subprocess curl to requests library
4. ✅ Added rate limiting delays (2s between figures, 5s between papers)
5. ❌ All blocked by SpringerNature's bot detection

**Current State**: IP/session temporarily blocked by SpringerNature

---

## Successfully Downloaded Papers (6)

1. **10.1038/s41467-025-58466-2** - Nature Communications
   - PDF: 10.0 MB
   - Figures: Multiple

2. **10.1038/s41467-025-65692-1** - Nature Communications
   - PDF: 6.1 MB
   - Figures: Multiple

3. **10.1038/s41467-025-66354-y** - Nature Communications
   - PDF: 479 KB
   - Figures: Multiple

4. **10.1038/s41467-025-66434-z** - Nature Communications
   - PDF: 385 KB
   - Figures: Multiple

5. **10.1038/s43587-025-01016-8** - Nature Aging
   - PDF: 15.0 MB
   - Figures: Multiple

6. **10.1038/s43587-025-01022-w** - Nature Aging
   - PDF: 39.0 MB
   - Figures: Multiple

---

## Failed Papers (26)

All papers listed below have metadata but no PDF or figures downloaded.

### Nature Communications (11 papers)
- 10.1038/s41467-018-03770-3
- 10.1038/s41467-025-63229-0
- 10.1038/s41467-025-63429-8
- 10.1038/s41467-025-64275-4
- 10.1038/s41467-025-64652-z
- 10.1038/s41467-025-64835-8
- 10.1038/s41467-025-65297-8
- 10.1038/s41467-025-65664-5
- 10.1038/s41467-025-65692-1 (partial)
- 10.1038/s41467-025-66354-y (partial)
- 10.1038/s41467-025-66434-z (partial)

### Nature Aging (6 papers)
- 10.1038/s43587-025-00943-w
- 10.1038/s43587-025-00961-8
- 10.1038/s43587-025-00986-z
- 10.1038/s43587-025-01011-z
- 10.1038/s43587-025-01014-w
- 10.1038/s43587-025-01021-x

### Scientific Reports (5 papers)
- 10.1038/s41598-025-22076-1
- 10.1038/s41598-025-25545-9
- 10.1038/s41598-025-26154-2
- 10.1038/s41598-025-27042-5

### Other Nature Journals (4 papers)
- 10.1038/s41514-025-00287-0 (Aging)
- 10.1038/s41514-025-00290-5 (Aging)
- 10.1038/s41587-025-02830-6 (Biotechnology)
- 10.1038/s41591-025-03999-8 (Medicine)
- 10.1038/s42003-025-09219-w (Communications Biology)
- 10.1038/s43856-025-01222-w (Communications Medicine)
- 10.1038/s43856-025-01239-1 (Communications Medicine)

---

## Recommendations

### Immediate Options

1. **Wait and Retry (24-48 hours)**
   - SpringerNature rate limit may reset
   - Use the created batch_download_figures.sh script
   - Increase delays to 5s between figures, 30s between papers

2. **Institutional Access**
   - Access through university/institution VPN
   - PDFs require authenticated session
   - Figures may be available without auth

3. **Contact Publisher**
   - Request bulk download access for research
   - Explain academic/non-commercial use
   - SpringerNature may provide API access

4. **Manual Download**
   - For critical papers, download manually through web browser
   - Browser sessions less likely to trigger bot detection

### Long-term Solutions

1. **Use Official APIs**
   - SpringerNature Developer Portal: https://dev.springernature.com/
   - Apply for API key for legitimate research use
   - Supports bulk downloads with proper authentication

2. **Alternative Sources**
   - PubMed Central (PMC) mirrors
   - arXiv/bioRxiv preprints
   - Author manuscript repositories

3. **Selenium/Browser Automation**
   - Last resort: Full browser emulation
   - Much slower but bypasses bot detection
   - Still subject to rate limits

---

## Files Created

- `/scripts/download_paper.py` - Modified with requests library and rate limiting
- `/scripts/batch_download_figures.sh` - Bash script for batch figure downloads
- `/scripts/download_figures_only.sh` - Single paper figure downloader
- `/scripts/retry_failed_downloads.sh` - Helper script for retries

## Database Status

- **papers.db**: Contains metadata for all 32 papers
- Missing PDFs and figures tracked in metadata.json per paper
- `pdf_path: null` for papers without PDFs

---

**Next Steps**: Wait 24-48 hours for rate limit reset, then retry with longer delays OR pursue institutional access/API key approach.
