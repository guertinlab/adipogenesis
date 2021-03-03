#!/bin/bash

#on Rivanna
#Install packages and put on $PATH

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
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
#restart console so cargo is added to $PATH automatically
git clone https://github.com/guertinlab/fqdedup.git
cd fqdedup
cargo build --release 
cp /target/release/fqdedup /home/bhn9by/bin

