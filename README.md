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

python dsbu_link_fdr_allin1.py  pep_list.csv  /path/to/spectra  DSBU_links_FDR.csv

Now it takes ~2 hours using 24 CPUs to finish the analysis. All tempt using GPU failed :-( too much memory used...

6. Output file details
| Column                     | How it is computed                                                                                                                                                                                                                                                                                                                                                                                                                  | Why it matters                                                                                                                                                                                                                            |            |                                                                                         |
| -------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------- | --------------------------------------------------------------------------------------- |
| **doublet\_present**       | *Boolean* returned by the `has_doublet()` function.<br>1. Read the **parent MS² spectrum** (centroid peaks).<br>2. Take the parent charge *z* from the “Z …” line.<br>3. Theoretical m/z gap for a DSBU α/β pair = 26 .016 Da / *z*.<br>4. Scan the peak list for **any two peaks** separated by that gap within ± 10 ppm.<br>5. Require both peaks to be ≥ 10 % of the spectrum’s base peak.<br>→ If found ⇢ `True`, else `False`. | Confirms the instrument actually saw the DSBU diagnostic pair that should trigger MS³; a quick sanity check on each cross-link.                                                                                                           |            |                                                                                         |
| **doublet\_int\_ratio**    | In step 4 above, take the *lower* of the two reporter intensities and divide by the **base-peak** intensity; round to three decimals.<br>E.g. if the reporters are 7 e4 and 5 e4, base peak is 2 e5 → ratio = 0.25.                                                                                                                                                                                                                 | A rough quality score for the reporter doublet; low values can mean weak cleavage or interference.                                                                                                                                        |            |                                                                                         |
| **decoy\_A**, **decoy\_B** | For each peptide sequence, look it up in `seq_info[sequence]['is_decoy']`, which is `True` if **any** of its protein IDs (column *Proteins*) or descriptions contains “decoy”, “rev\_”, “reverse”, etc. (regex \`(decoy                                                                                                                                                                                                             | reverse                                                                                                                                                                                                                                   | rev\_)\`). | Marks whether each peptide matched a target or a decoy sequence in the database search. |
| **decoy\_link**            | Logical OR: `decoy_A or decoy_B`                                                                                                                                                                                                                                                                                                                                                                                                    | A cross-link is counted as *decoy* if **either side** hit a decoy protein.                                                                                                                                                                |            |                                                                                         |
| **link\_score**            | `XCorr_A + XCorr_B` (sum of the ProLuCID cross-correlation scores of the two MS³ PSMs).                                                                                                                                                                                                                                                                                                                                             | A single ranking value that combines evidence from both peptides; higher ⇢ more likely correct.                                                                                                                                           |            |                                                                                         |
| **rank**                   | Links are sorted by **descending** `link_score`; `rank = row_index + 1` so the best-scoring link is rank 1.                                                                                                                                                                                                                                                                                                                         | Gives a stable order for cumulative FDR statistics.                                                                                                                                                                                       |            |                                                                                         |
| **cum\_links**             | Identical to `rank` (running total of links encountered).                                                                                                                                                                                                                                                                                                                                                                           | Needed for cumulative fractions.                                                                                                                                                                                                          |            |                                                                                         |
| **cum\_decoys**            | Running sum of `decoy_link == True` from the top of the table down to the current row.                                                                                                                                                                                                                                                                                                                                              | Tracks how many false (decoy) links are accumulating as you lower the score threshold.                                                                                                                                                    |            |                                                                                         |
| **FDR**                    | Target–decoy estimate, **unit-less probability**:<br>$FDR = 2 × (cum_decoys / cum_links)$<br>The factor 2 projects the observed decoy frequency onto the target population because the database contains a 1:1 mix of target and decoy sequences.                                                                                                                                                                                   | At any row it answers: “If I keep every link **at this score or higher**, what fraction do I expect to be wrong?” 0.01 = 1 %. Values can exceed 1 when decoys outnumber targets far down the list; those rows are never kept in practice. |            |                                                                                         |


