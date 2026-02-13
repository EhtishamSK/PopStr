# Population Structure Analysis Pipeline

**Purpose**  
This R script processes output files from the classic **STRUCTURE** software (used in population genetics to infer genetic clusters / admixture). It helps determine the **optimal number of clusters (K)** and creates clean, publication-quality bar plots showing population structure.

**Main Workflow**

1. **Package setup**  
   - Installs (if needed) and loads:  
     - `devtools`, `ggplot2`  
     - `label.switching` (for dealing with label switching in STRUCTURE runs)  
     - `pophelper` (via devtools/github – main tool for reading + plotting STRUCTURE results)

2. **Data input**  
   - Sets working directory  
   - Reads multiple STRUCTURE output files (usually `_f` files for different K values)

3. **Summarization & Quality check**  
   - Uses `readQ()` and `summariseQ()` from pophelper to combine and summarize the results across replicates and K values

4. **Selecting the best K**  
   - Applies the widely-used **Evanno method** (ΔK statistic) to identify the most supported number of genetic clusters

5. **Visualization**  
   - Creates a **STRUCTURE-style bar plot** with `plotQ()`  
   - Customizations include:  
     - Sorted individuals (by highest ancestry proportion)  
     - Clear colors, labels, and grouping by populations  
     - Publication-ready aesthetics (clean theme, good resolution)

**Key Features / Strengths**
- Automates post-processing of STRUCTURE runs
- Handles label switching problem
- Uses Evanno's ΔK method to objectively choose K
- Produces nice-looking plots with minimal extra code

**Typical Use Case**  
Run STRUCTURE for K=2 to K=10 (with replicates). Point this script at the output folder, get optimal K + beautiful admixture plots for papers/thesis.
Quick to adapt for other datasets, just update file paths!! 

**Final Plot**  
<img width="5400" height="2100" alt="STRUCTURE_K3_Sorted" src="https://github.com/user-attachments/assets/772d158e-d664-4570-bf5f-da3ccd3b0b4a" />
