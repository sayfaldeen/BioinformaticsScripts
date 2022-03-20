#!/bin/bash

esearch -db biosample -query ${1} | elink -target nuccore | efetch -format fasta | gzip -f > ${1}.fa.gz

if [[ $? == 0 ]]
then
	echo "${1}.fa.gz has been downloaded"
else
	echo "Biosample ${1} had a non-zero exit status"
fi
