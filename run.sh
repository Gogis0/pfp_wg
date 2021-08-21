#!/bin/bash

set -euxo pipefail

bigrepair/ctph/newscanNT.x data/yeast.fasta -w 4 -p 11 -c
bigrepair/procdic data/yeast.fasta.dicz
bigrepair/largeb_repair/irepair data/yeast.fasta.dicz.int 7482
bigrepair/largeb_repair/irepair data/yeast.fasta.parse 7482

# bigrepair/repair/despair yeast.fasta
