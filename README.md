# genetic-drift-evolution
Data and code for the manuscript, "The evolution of genetic drift over 50,000 generations" by Ascensao et al.

# Description of repository
 
### code/

This folder contains all of the scripts used to analyze the barcode sequencing data. Use `run_pipeline.sh` to run the entire analysis pipeline--from raw .fastq sequencing files to sampling from posteriors of the genetic drift parameters.

- `run_MultiCodes.sh`: run MultiCodes.pl to parse barcodes from raw .fastq sequencing files
- `error_correction.sh`: perform error correction to correct for sequencing errors in barcodes
- `count_bc.sh`: count how many times we see each barcode for each sample
- `kt_errors.sh`: compute variance parameters for outlier detection
- `outliers.sh`: detect outlier barcodes--those with trajectories that deviate significantly from neutrality
- `filter_bcs.sh`: filter out outlier barcodes, merge low-frequency barcodes together
- `HMM_stan.sh`: MCMC inference of genetic drift parameters

### experiments.csv

Table of experiments--metadata that gives details of each strain.

### E1/, E23/, E4/, E5/, E6/

Folders that contain data from each batch of experiments. Data from experiments 2 and 3 are combined in the same folder--E23.
- `{.}_meta.csv`: Metadata for the experiment showing the map between Illumina dual index and sample
- `{.}_exps.csv`: All strains used in the batch of experiments
- `{.}_cfus.csv`: Colony forming unit counts
- `data/MultiCode_tables/`: Output of `MultiCodes.pl`--parsed barcodes from raw .fastq files
- `data/errorcorrectedcodes/`: Output of error correction, merge barcodes into "super-barcodes"
- `data/bc_counts/`: Counts and frequencies of super-barcodes
- `data/kt_errors/`: Variance parameters for outlier detection
- `data/outliers/`: Inferred fitness of super-barcodes for outlier detection
- `data/drift/`: Samples from posteriors of genetic drift parameters
