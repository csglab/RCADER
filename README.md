# RCADER
### Recognition Code-Assisted Discovery of regulatory Elements by Regression

RCADER uses ChIP-seq data and the C2H2-ZF recognition code to identify motifs that represent the binding preferences of C2H2-ZF proteins. It uses an iterative regression-based approach to optimize the motifs that are predicted by the recognition code, so that the correlation of motif scores and the ChIP-seq peak heights are maximized.

## Requirements 
- Unix-compatible OS
- [R](http://www.r-project.org/) version 3.0.1 or later
- R [randomForest](https://cran.r-project.org/web/packages/randomForest/index.html) library
- [GNU Scientific Library](https://www.gnu.org/software/gsl/) version 1.15 or later

### Installation

From the command line, go to the root folder of RCADER, and run `make`. This will compile the codes and will create the RCADER executables in `./bin/`.

### Running RCADER


#### Input files

For running RCADER, you need the following files:

* **C2H2-ZFP sequence file:** The protein sequence of the C2H2-ZF protein that is targeted by ChIP-seq (FASTA format).
* **Peak sequences:** The DNA sequences of the peak regions (FASTA format). We recommend the use of &plusmn;250 bp region around the peak summits.
* **Peak scores:** The scores (e.g. fold-enrichment) associated with peaks (tabular format). The first column represents peak names (must match the names in the peak sequence file), and the second column represents the scores (must be >0).


#### Usage

To run RCADER, go to the main folder of the package, where the file `RACDER.sh` is located. Then simply use the following command:

```bash
bash RCADER.sh <jobID> <C2H2_ZFP.fasta> <ChIP_seq.fasta> <ChIP_seq_scores.txt>
```

The `<jobID>` parameter is a name that will be used to create the output folders.

#### Output

RCADER creates the following main output files in the `./out/<jobID>` folder:

* `results.report.txt`: A table with a summary of optimization results for each of the recognition code seed motifs.
* `results.ps`: A PostScript file containing a graphical representation of each seed and optimized motif.
* `results.PFM.txt`: A text file containing the optimized motifs as position frequency matrices (PFMs).
* `results.scores.txt`: A table with the original peak scores, the nucleotide and dinucleotide frequences of each peak, and the score of each optimized motif for each peak.

#### Example datasets

Example datasets are provided at `./examples`, including a subsampled set of peaks for CTCF ChIP-seq. Use the following command to run RCADER on this example dataset:
```bash
bash RCADER.sh example_CTCF ./examples/CTCF/CTCF.fasta ./examples/CTCF/summits.fa ./examples/CTCF/macs.subsampled.uniform.txt
```
This will create the `./out/example_CTCF` folder, which will contain the RCADER output files.
