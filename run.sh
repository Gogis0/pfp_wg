#!/bin/bash

set -euxo pipefail

# in: 	yeast.fasta
# out: 	yeast.fasta.dicz  yeast.fasta.dicz.len  yeast.fasta.parse  yeast.fasta.parse_old
../bigrepair/ctph/newscanNT.x 	yeast.fasta -w 4 -p 11 -c

# in: 	yeast.fasta.dicz  yeast.fasta.dicz.len
# out: 	yeast.fasta.dicz.int
../bigrepair/procdic 		yeast.fasta.dicz

# in:	yeast.fasta.dicz.int
# out:	yeast.fasta.dicz.int.C  yeast.fasta.dicz.int.R
../bigrepair/largeb_repair/irepair yeast.fasta.dicz.int 15005

# in:	yeast.fasta.parse
# out:	yeast.fasta.parse.C  yeast.fasta.parse.R
../bigrepair/largeb_repair/irepair yeast.fasta.parse 15005

# in:	yeast.fasta.dicz.int.C  yeast.fasta.dicz.int.R  yeast.fasta.parse.C  yeast.fasta.parse.R
# out:	yeast.fasta.C  yeast.fasta.R
../bigrepair/postproc 		yeast.fasta

# in:	yeast.fasta.C  yeast.fasta.R
# out:	yeast.fasta.out
../bigrepair/repair/despair 	yeast.fasta

if [ $(diff yeast.fasta yeast.fasta.out | wc -l) ]
then
	echo "test successful"
else
	echo "test failed"
fi
