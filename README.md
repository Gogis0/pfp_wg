Original work:                              \
https://gitlab.com/manzai/bigrepair         \
https://github.com/simongog/sdsl-lite       \
https://github.com/waYne1337/BWT-Tunneling  \

This repository is not used \
https://gitlab.com/manzai/Big-BWT

# How to install and run test
```
git clone https://github.com/simongog/sdsl-lite.git 
cd sdsl-lite
./install.sh ..
cd ..

cd bigrepair
make
cd ..

cmake .
make

conda create -n bigrepair psutil
conda activate bigrepair

./bigrepair/bigrepair -w 4 -p 11 data/yeast.fasta
./tfm_index_construct           \
    -i                          \
    -sa DIVSUFSORT              \
    data/yeast.fasta.parse      \
    data/yeast.fasta.tunnel
```