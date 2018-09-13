#!/bin/bash

# $1 is the data directory 
# $2 is the metadata file (inside aforementioned data directory) - must contain headers "RNAfilename" and "DNAfilename"

# output: generates filenames.txt containing a single filename per line, for all files listed in the metadata
#         ...this is used as input to obtain all the pathogen report files downstream.

echo $1
echo $2

file="$1$2"; shift

awk -v cols="RNAfilename DNAfilename" '
BEGIN{
	split(cols,C)
	OFS=FS=","
	getline
	split($0,H)
	for(c in C){
		for(h in H){
			if(C[c]==H[h])F[i++]=h
		}
	}
}
{ l="";for(f in F){l=l $F[f] OFS}print l }

' "$file" > filenames


(cut -d',' -f1 filenames; cut -d',' -f2 filenames) | cat > filenames.txt

rm filenames
