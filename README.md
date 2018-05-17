# Precision-Recall Assessment of Protein Function Prediction

## Introduction
Critical Assessment of Function Annotation (CAFA), is a community-wide challenge designed to provide a large-scale assessment of computational methods dedicated to predicting protein function.

More information can be found at http://biofunctionprediction.org/cafa/ as well as the CAFA2 paper (Jiang et al, 2016)

This toolset provides an assessment for CAFA submissions based on precision and recall. 

For bug reports, comments or questions, please email nzhou[AT]iastate.edu.

## Dependencies
 - Python 2.7 or Python 3
 - Python packages can be downloaded from their sites or installed from repositories:
    1. [Biopython](http://biopython.org)
    2. [yaml](http://www.yaml.org/download.html)
    3. [matplotlib](https://matplotlib.org/) 
    4. [seaborn](https://seaborn.pydata.org/)
    
`$ sudo apt install python-biopython python-yaml python-matplotlib python-seaborn`

## Main Functions
 We provide two main functions to assist in the evaluation of GO-term prediction within the scope of CAFA, the main assessment function and the plot function.
 - `assessment.py` 
	- Only input needed is the configuration file `config.yaml`, where the following parameters are specified in the first section `assessment`.
	- First parameter `file`           : prediction file formatted according to [CAFA3 formats](https://www.synapse.org/#!Synapse:syn5840147/wiki/402192)
	- Second parameter `ic_path`       : path of the calculated infromation content map. This must be made prior to running `assessment.py`. This is achieved by running `IC.py` in the ICTool directory, instructions below in auxillary functions.
	- Third parameter  `obo_path`      : path of the gene ontology obo file. The latest version can be downloaded [here](http://purl.obolibrary.org/obo/go.obo). Note that the obo file used here should not be older than the one used in the prediction.
	- Fourth parameter `benchmark_path`: directory of the benchmark folder. Specific formats are required for the benchmark folder, including two sub-directories: groundtruth and lists. Please refer to auxiliary function `benchmark_folder.py` for the creation of this folder, as well as the genral creation of benchmarks. An example benchmark folder is given in this repository `./assessment/benchmark`
	- Fifth parameter  `results_path`  : Folder where results are saved. A `rawdata` folder will be created within the results folder.
	- Seventh parameter `gaf_path`     : Location of the GAF File with respect to the ICTool directory
	- Note that only the first section `assessment` of the configuration file is used here, the rest of the configuration file can be ignored for this function	
 - `plot.py`
	- Only input needed is the configuration file `config.yaml`, where the following parameters are specified in the second section `plot`.
	- First parameter `results`: the results from the `assessment.py` function.
	- Second parameter `title`: title of the plot. Optional.
	- Third parameter `smooth`: whether the precision-recall curves should be smoothed. Input 'Y' or 'N'. 
	- Fourth parameter(s) `fileN`: name of the result file to be plotted. Can add up to 12 files. These results will be drawn on the same plot.
	- Example: if the prediction file is `ZZZ_1_9606.txt`, the result file in the results folder will be `ZZZ_1_9606_results.txt`. Only input `ZZZ_1_9606` in the above parameter for plotting. 



## Auxiliary Functions 
 CAFA3 released its [protein targets](https://www.synapse.org/#!Synapse:syn6172284) in September 2016. Each protein target has a unique CAFA3 ID. To run the above assessment function, each protein should be represented by its CAFA3 ID. However, the benchmark proteins generated by the [benchmark creation tool](https://github.com/nguyenngochuy91/CAFA_benchmark) are identified by UniProt Accession IDs. 
Therefore, we here provide functions to convert between UniProt IDs and CAFA3 IDs. We also provide a function that converts benchmark files generated by the [benchmark creation tool](https://github.com/nguyenngochuy91/CAFA_benchmark) to a benchmark folder that can feed into this program.
 - `benchmark_folder.py`
	- Refer to `python benchmark_folder.py -h` for syntax of using this function by itself.
	- If using our [benchmark creation tool](https://github.com/nguyenngochuy91/CAFA_benchmark), then the `benchmark_pipeline.sh` file is a good example of how to generate a benchmark folder for `assessment.py` from the raw benchmarks.
	- Input your own folder names and different gaf file names in the blanks left in `benchmark_pipeline.sh`.
 	
 - `./ID_conversion/ID_conversion.py` 
	- Two functions are written in this python script, one converts UniProt Accessions to CAFA3 IDs, the other function converts the other way around.
	- First function `uniprotac_to_cafaid(taxon, uniprotacs)`.
	- Second function `cafaid_to_uniprot(taxon, cafaids)`.
	- Refer to comments in the script `./ID_conversion/ID_conversion.py` and third example below for usage.
 - `./ICTool/IC-GAF.py` 
	- Calculates the Information content of a given OBO
	- How to use: From the ICTool directory `IC.py ../config.yaml`.
	-  It uses the assessement config file, so once you have that setup, this will run on the same file set.

 - `./ICTool/IC-LIST.py` 
	- Calculates the Information content of a given a list of Proteins, Terms and Ontologies.
	- How to use: From the ICTool directory `IC.py ../config.yaml`.
	-  It uses the assessement config file, so once you have that setup, this will run on the same file set.

## Examples	
 - `./assessment.py config.yaml`
 - `./plot.py config.yaml`
 - `./ID_conversion/ID_conversion.py ./ID_conversion/example_uniprot_accession_8355.txt 8355 ./ID_conversion/example_output.txt`
 - `./benchmark_pipeline.sh`
	 


## References
Jiang, Yuxiang, et al. "An expanded evaluation of protein function prediction methods shows an improvement in accuracy." Genome biology 17.1 (2016): 184.

http://biofunctionprediction.org/cafa/
