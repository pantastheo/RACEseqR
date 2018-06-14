## RACEseqR

A Next Generation Sequencing data analysis pipeline packed as an R package. The pipeline uses popular command line tools for the downstream analysis of a 5-RACE sequencig experiment.

## Dependencies

The following programs and packages must be installed on your computer.

Package            | Version
--------------     | --------------
Cutadapt (optional)| >=1.12
Bowtie             | 1.0.0
Samtools           | >=1.3.1
Bedtools           | >=2.17.0
Tmap (optional)    | >=3.4.1

## Installation

The best way to install RACEseqR is to use the install_github command from the devtools package.

```
library("devtools")
install_github("pantastheo/RACEseqR")
library("RACEseqR")
```


## Usage

Either use RACEseqR to output graphs dirrectly to your workinf directory.

```
RACEseq(input_data, replicon_ref, mismatch = 0, RACE_adapter = NULL,
  str = NULL, end = NULL, filename, tmap = NULL)
```
Or use RACEseqR inside Rstudio to manually generate tables and graphs.

```
dataframe <- RACEseq(input_data, replicon_ref, mismatch = 0, RACE_adapter = NULL,
  str = NULL, end = NULL, filename, tmap = NULL)
```


## Arguments

Arguments       | Description
--------------  | --------------
input_data      | The input .fastq files to be used for the alignment.
replicon_ref    | The reference sequence to be used in the alignment.
mismatch        | The number of mismatches to be accepted during alignment.
RACE_adapter    | The RACE adapter to be trimmied using cutadapt .
str             | The start position of the binding region.
end             | The end position of the binding region.
filename        | The filename that will be used for the output pdf file.
tmap            | The option to use the Tmap aligner. Defaults NULL

