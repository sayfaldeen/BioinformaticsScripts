# nanosim-bac-wrapper.py

- Script that wraps around [nanosim](https://github.com/bcgsc/NanoSim) to automate the creation of synthetic metagenomes composed of purely bacterial species

## Installation
- `wget https://raw.githubusercontent.com/sayfaldeen/BioinformaticsScripts/main/LongReads/nanosim-bac-wrapper.py` to download the script
- `wget https://raw.githubusercontent.com/sayfaldeen/BioinformaticsScripts/main/LongReads/sayf-nanosim.yml` to download the environment
- `conda env create -f sayf-nanosim.yml`

## Quick usage
- `nanosim-bac-wrapper.py --n-genomes 50 --n-samples 200 --read-count 650000 --variable`
	- Creates a metagenomic community of 200 samples with 50 bacterial genomes in each sample. Each sample will also have a total of 650,000 total reads and the relative abundances of each bacterial genome will be randomly variable.
		- If you wsh to keep the relative abundances uniform for all bacterial genomes within all samples, do not include the `--variable` option
