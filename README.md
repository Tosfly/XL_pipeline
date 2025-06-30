# XL_pipeline
Find the XL pairs using IP2
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
python dsbu_link_fdr_allin1.py  pep_list.csv  /path/to/spectra  DSBU_links_FDR.csv

6. Output file details
| Column                  | Calculation                                                                                             | Meaning                                                                                                        |
| ----------------------- | ------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------- |
| **doublet\_present**    | Parent MS² has two peaks Δm/z = 26.016 Da/charge within ±10 ppm and ≥10 % base-peak.                    | Verifies DSBU diagnostic doublet.                                                                              |
| **doublet\_int\_ratio** | (Weaker reporter ÷ base-peak) intensity (0 – 1).                                                        | Strength of the doublet.                                                                                       |
| **decoy\_A / decoy\_B** | `True` if peptide’s protein locus contains “decoy / rev\_ / reverse”.                                   | Target–decoy bookkeeping.                                                                                      |
| **decoy\_link**         | `decoy_A OR decoy_B`.                                                                                   | A link is false if **either** peptide is decoy.                                                                |
| **link\_score**         | `XCorr_A + XCorr_B`.                                                                                    | Unified score for ranking.                                                                                     |
| **rank**                | Row number after sorting by *descending* link\_score (1 = best).                                        |                                                                                                                |
| **cum\_links**          | Running total (= rank).                                                                                 | Denominator for FDR.                                                                                           |
| **cum\_decoys**         | Running total of decoy links.                                                                           | Numerator for FDR (before ×2).                                                                                 |
| **FDR**                 | `2 × cum_decoys ÷ cum_links` (unit-less, 0–1+).                                                         | Estimated false-link fraction at or above the threshold. Values >1 appear only far down where decoys dominate. |
| **Note**                | Stub pattern: `true XL` (+85/+111), `Rare, may not be true` (+85/+85 or +111/+111), `Mono`, `Not sure`. | Quick triage before modelling.                                                                                 |



