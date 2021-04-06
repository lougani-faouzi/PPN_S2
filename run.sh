#!/bin/bash

#Nettoyage de l'environnement
make clean

#Compilation
make

#Execution du split
./bin/split2 sequences.fasta

#Pour chaque fichier sequence on retire les retours a la ligne, \n
for f in `ls 6* 7* L* MN* MT* MW* N*`
do
   tr -d "\n" < $f > tmp
   mv tmp $f
done

#Execution des analyses
./bin/genysis LC528232.1 6XEZ_T
