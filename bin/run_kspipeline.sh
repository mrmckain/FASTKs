#!/bin/bash

cd `pwd`
CURRDIR=`pwd`
for f in *.cdna
do
        echo $f
        NAMESP=$(echo $f | cut -d "." -f 1)
        mkdir $NAMESP
        cd $NAMESP
        cp ../$NAMESP.pep .
        cp ../$NAMESP.cdna .
        time makeblastdb -in $f -dbtype nucl
        time blastn -num_threads 2 -query $f -db $f -outfmt 6 -evalue 1e-40 > $f.blastn
        perl ~/bin/FAST_kspipe.pl --query $f --subj $f --blast $f.blastn --proc 4 --pairs 100
        cd $CURRDIR
done
