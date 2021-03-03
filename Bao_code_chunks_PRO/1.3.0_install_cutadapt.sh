#!/bin/bash

#Install cutadapt binary for trimming adapter

#Retrieve cutadapt binary and put on the $PATH
python3 -m pip install --user --upgrade cutadapt
cp /home/bhn9by/.local/bin/cutadapt /home/bhn9by/bin
