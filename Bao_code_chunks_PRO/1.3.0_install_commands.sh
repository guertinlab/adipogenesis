#!/bin/bash

#Retrieve packages and put on the $PATH

#cutadapt
python3 -m pip install --user --upgrade cutadapt
cp /home/bhn9by/.local/bin/cutadapt /home/bhn9by/bin

#fastx-toolkit
mkdir fastx_bin
cd fastx_bin
wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
tar -xjf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
cp ./bin/* /home/bhn9by/bin

#fqdedup
