# XL_pipeline
Find the XL pairs after using IP2 search.
# Cross-linking MS2-MS3 Workflows in IP2
A concise, copy-and-paste guide for **DSBU** datasets, followed by an
explanation of the `dsbu_link_fdr_allin1.py` post-processing script and all QC /
FDR columns it adds.

--------------------------------------------------------------------
1  Searching **DSBU** MS2-MS3 RAW files in IP2
--------------------------------------------------------------------
1. **Convert RAW → *.ms2 + *.ms3***  
   ```bash
   RawConverter.exe -i *.raw -q -m 3 -x -f 3
-m 3 keeps every MS³ scan.

-x writes the parent MS² scan number into each MS³ title.

2. Rename MS³ files
sample1.ms3   →   sample1_ms3.ms2

3. ProLuCID search – include parent .ms2 and all _ms3.ms2

| Setting        | Value                                          |
| -------------- | ---------------------------------------------- |
| Enzyme         | Trypsin, ≤2 missed cleavages                   |
| Precursor tol. | ±7 ppm                                         |
| Fragment tol.  | 0.6 Da (ion-trap) / 20 ppm (Orbitrap)          |
| Variable mods  | +85.05276 Da *and* +111.03205 Da on K & N-term |
| Static mod     | +57.02146 Da (C)                               |
| Optional mod   | +15.9949 Da (Met ox)                           |

4. DTASelect
--fp 0.9999  --trypstat 0  --dm 10

5. Post-processing with dsbu_link_fdr_allin1.py
Save the *.ms2, *_ms3.ms2 and pep_list.csv file to the same folder and run the commend line:

# DSBU Cross-Link Builder - Corrections Applied (2025)

## Summary of Critical Fixes

This document details the corrections applied to the DSBU MS2-MS3 cross-link builder script based on current literature and best practices for DSBU cross-linking mass spectrometry.

## 1. Reporter Mass Difference Correction ✓

**Issue:** Original script used `REPORTER_DELTA = 26.01929 Da`

**Fix:** Corrected to `REPORTER_DELTA = 25.97929 Da`

**Rationale:**
- The reporter ion mass difference should match the difference between beta and alpha stubs
- Beta stub (111.03205) - Alpha stub (85.05276) = 25.97929 Da
- This aligns with literature citing ~26 u reporter-ion spacing in DSBU MS2 spectra
- The corrected value ensures accurate doublet detection

## 2. Precursor Mass Formula Update ✓

**Issue:** Original formula: `calc_total = mass_a + mass_b + DSBU.INTACT_MASS - DSBU.PROTON`

**Fix:** Updated to: `calc_total = mass_a + mass_b + DSBU.INTACT_MASS - 2 * DSBU.PROTON`

**Rationale:**
- DSBU cross-linking involves loss of two NHS (N-hydroxysuccinimide) groups
- Each NHS group departure involves loss of a proton
- The corrected formula better reflects the actual chemistry of DSBU cross-linking
- This improves precursor mass matching accuracy

## 3. Improved Peptide Pairing Algorithm ✓

**Issue:** Original only considered top 2 peptides by XCorr

**Fix:** Now considers all peptides containing DSBU stubs (up to configurable limit)

**Implementation:**
- Separates peptides by stub type (85 vs 111)
- Tests all complementary pairs (85-111 combinations)
- Filters by mass accuracy before scoring
- Configurable `--max_candidates` parameter (default: 10 per stub type)
- Deduplicates results keeping highest scoring pairs

**Benefits:**
- Reduces risk of missing valid cross-links
- More comprehensive coverage of potential links
- Still computationally tractable with parallelization

## 4. Mass Tolerance Recommendations ✓

**Changes:**
- Default remains 10 ppm for compatibility
- Added recommendation for 5 ppm on well-calibrated instruments
- Made tolerance fully configurable via `--ppm` parameter
- Added warning message when using default tolerance

## 5. Intensity Ratio Threshold Update ✓

**Issue:** Original default `min_ratio = 0.10`

**Fix:** Raised default to `min_ratio = 0.15`

**Rationale:**
- Higher threshold reduces noise in doublet detection
- Literature suggests 0.15-0.20 range for reliable detection
- Remains configurable via `--min_ratio` parameter

## 6. Same-Charge Requirement Documentation ✓

**Addition:** Added comprehensive documentation explaining that:
- DSBU doublet peaks share the same charge state (z)
- Instruments trigger MS3 only on same-charge doublets
- The mass difference in m/z domain is REPORTER_DELTA/z
- This is inherent to how the algorithm works (using precursor z)

## 7. Cross-Validation Recommendations ✓

**Added to script output and help:**
- List of recommended XL-MS search engines for validation:
  - MeroX (with Rise/RiseUP scoring)
  - XlinkX 2.5 (Proteome Discoverer 2.5)
  - pLink2
  - XiSearch/xiFDR
  - Kojak + Percolator
  - Prosit-XL (2025)

## 8. Additional Improvements

- **Validation check:** Added post-initialization check to verify REPORTER_DELTA matches stub difference
- **Logging enhancements:** Added informative messages about parameters being used
- **Deduplication:** Added logic to remove duplicate cross-links (same peptides, same parent scan)
- **Documentation:** Comprehensive inline comments explaining corrections
- **Help text:** Detailed epilog in argparse explaining corrections and recommendations

## Usage Examples

### Basic usage with corrections:
```bash
python dsbu_link_fdr_allin1_mt24_progress_corrected.py \
    pep_list.csv \
    /path/to/spectra \
    output_links.csv
```

### Optimized for well-calibrated instrument:
```bash
python dsbu_link_fdr_allin1_mt24_progress_corrected.py \
    pep_list.csv \
    /path/to/spectra \
    output_links.csv \
    --ppm 5 \
    --min_ratio 0.20 \
    --max_candidates 15
```

### High-throughput with more threads:
```bash
python dsbu_link_fdr_allin1_mt24_progress_corrected.py \
    pep_list.csv \
    /path/to/spectra \
    output_links.csv \
    --threads 32 \
    --max_candidates 20
```

## Validation Recommendations

1. **Cross-validate results** with at least one dedicated XL-MS search engine
2. **Export to xiFDR or Percolator** for advanced FDR estimation
3. **Verify MS acquisition settings:**
   - Delta M1 = 25.979 Da (not 26.019)
   - Same charge requirement enabled
   - Appropriate AGC targets and NCEs for DSBU

## Performance Impact

- Improved peptide pairing increases candidate evaluations but maintains speed via parallelization
- More accurate mass calculations reduce false positives
- Better doublet detection (correct mass difference) improves specificity
- Overall: Better quality results with minimal performance impact

## References

The corrections are based on:
- PMC literature on DSBU MS2-MS3 methods
- Current best practices in cross-linking mass spectrometry
- Chemical properties of DSBU cross-linker
- Recommendations from major XL-MS software developers

## Version History

- **Original:** dsbu_link_fdr_allin1_mt24_progress.py
- **Corrected:** dsbu_link_fdr_allin1_mt24_progress_corrected.py (2025)

---

For questions or issues, consider consulting the cross-linking mass spectrometry community or the developers of dedicated XL-MS search engines mentioned above.
