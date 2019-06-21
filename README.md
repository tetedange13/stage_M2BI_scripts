# Jan 2019: 6-months M2 (BIB) internship


<p align="center"><img src="toblerone_team.png" width="70%"></p>  

## Description of this respository
This Github repository contains scripts developped within the scope of


## Installation
- As this project has been coded to be used on a closed environment 
(cluster) containing all generated fastq files, the different scripts cannot
be run on an external machine

- The Python packages required are:
  ```
  python       3.6.0 or higher
  numpy        1.15.2 or higher
  biopython    1.72 or higher
  docopt       0.6.2 or higher
  pandas       0.22.0 or higher
  sklearn      0.17 or higher
  pysam        0.15.2 or higher
  ```

- The R librairies required are:
```
  ggplot2
  dplyr
  reshape2
  fmsb
```


## Usage
The program is composed of 4 python scripts (`0-solve_SILVA.py`, `1-pipeline.py` 
, `2-prim_analysis.py`, `3-second_analysis.py`), and 2 R scripts 
(`make_stackbar.R`, `make_radar_plot.R`). All in the base directory.  
All modules are stored in the `src/` directory. <br>



- Show help: <br>
  ```
    Mapping statistics computation

    Usage:
      stats_tax.py (-i <inFile>) (-l <taxoCut>)
      
    Options:
      -h --help                  help
      --version                  version of the script
      -i --inFile=input_file     input file
      -l --taxoCut=taxo_cutoff   cutoff for the taxonomic level 
  ```


## Pre-processing
- The `0-solve_SILVA.py` script was dedicated to solve all issues linked to
taxonomy. <br>
It contains several functions, that have to be called by modifying directly
the 'main'. <br> 
- For example to run the function producing 


## Results of Benchmark
- Example of stdout produce with  plot that can be produced:

- Example of R plots that can be produced using R scripts:


## Features
- General pipeline to trim adapters ([Porechop](https://github.com/rrwick/Porechop)), 
detect chimeras ([yacrd](https://github.com/natir/yacrd)), proceed to taxonomic
determination (either Minimap2 or Centrifuge), against different possible
databases (ZYMO, RRN, SILVA or p_compressed, with Centrifuge only)


## Main contributor
- [Félix Vandermeeren](https://github.com/tetedange13)

## Troubleshootings
