# AML-RoaMeR

This repository contains preset instructions to reproduce relative mutation rate (RMR) statistics and figures reported in "[Germline loss of *MBD4* predisposes to leukaemia due to a mutagenic cascade driven by 5mC](https://www.biorxiv.org/content/early/2017/11/01/180588)" (Sanders & Chew et al., 2017). Little to no adjustments will be made over the course of time as to keep the results in line with the manuscript. The current development snapshots of [RoaMeR](https://github.com/MathijsSanders/RoaMeR) and [superFreq](https://github.com/ChristofferFlensburg/superFreq) can be found at their respective repositories. 

## How do I run it?

The statistics and figures reported in the manuscript can be generated by two separate means dependent on the wishes of the user.

### The recommended way

The following command downloads precalculated genome-wide methylation tables based on CD34+ WGBS data from [Spencer et al., 2017, Cell](http://www.cell.com/cell/fulltext/S0092-8674(17)30106-X). In addition superFreq is installed and all statistics and figures are automatically produced.

```bash
./reproduceResults.sh -i input_directory -o output_directory
```

-i: Input directory containing RoaMeR.jar and the directories "Auxiliary" and "Intervals"
-o: Output directory where all results are going to be stored

*Dependencies*

- Java JRE 1.7+
- R 3.3+
- Bioconductor

### The entire way (Requires large memory)

The following command compiles RoaMeR through Maven, constructs a CG genome content table and in conjunction with the CD34+ WGBS data from [Spencer et al., 2017, Cell](http://www.cell.com/cell/fulltext/S0092-8674(17)30106-X) constructs a genome-wide 5mCG table. This step requires downloading of CD34+ WGBS bigWig files from [here](https://wustl.app.box.com/s/6fvifntz6golh4vu52epkfqxeyd0vqmk/folder/11145725870) (both methylation-specific as well as coverage bigWig files). Finally, RMR statistics and figures are produced following the construction of all necessary tables. 

```bash
./reproduceResults_complete.sh -i input_directory -m wgbs_directory -o output_directory [-g reference_genome]
```

-i: Input directory containing RoaMeR.jar and the directories "Auxiliary" and "Intervals"
-o: Output directory where all constructed tables and results are going to be stored
-m: Directory containing the CD34+ bigWig files (coverage as well as methylation)
-g: (Optional) FASTA file containing the reference sequence. When undefined it will automatically start downloading hg19

*Dependencies*
- Maven version 3+
- Java JDK 1.7+
- R 3.3+
- Bioconductor
