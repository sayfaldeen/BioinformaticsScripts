#!/bin/bash

#$1: FQ file
#$2: Threads
#$3: Output file


# Set up the option parser
bname=$(basename $0) # Pull the filename

# Define the help function
function help {

	echo " "
	echo "------------------------------------"
	echo "| ${bname} 'help' |"
	echo "------------------------------------"

	echo "${bname} -f <input file.fa/fq> -t <number of threads>"
	echo ""

	echo "-f: Input file name (fa|fq)"
	echo "-t: Number of threads to use (default: 2)"
	echo " "

}

# Create the option parser
clear # Clear the output
	# Put in the default arguments
THREADS=2

optstring=":hf:t:"
while getopts ${optstring} opt
do
	case ${opt} in
		h)
			help
			exit 0;;
		f)
			F="${OPTARG}"
			echo "Input file: ${F}";;
		t)
			THREADS="${OPTARG}"
			echo "Number of threads to use: ${THREADS}";;

		?)
			echo ""
			echo "Invalid input: -${OPTARG}"
			help
			exit 1;;
		:)
			echo ""
			echo "Must provide input"
			help
			exit 1;;
	esac
done
echo ""


source /home/sayf/miniconda3/etc/profile.d/conda.sh
conda activate kraken2

rep=$(echo ${F} | awk -F "." '{print $1".kreport"}')
out=$(echo ${F} | awk -F "." '{print $1".kraken.out"}')
kraken2 --use-names --db /home/sayf/utils/Databases/Kraken/ ${F} --threads ${THREADS} --report ${rep} > $out


title=$(echo $F | awk -F ".f" '{print $1}')
./KrakenPiePlots.py -f ${rep} --title ${title} --out "${title}.png"
