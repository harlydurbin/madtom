#!/bin/bash

while IFS='' read -r LINE || [ -n "${LINE}" ]; do
    zcat ${LINE} | head -n 1 >> ../output/readnames.txt
done < ../output/fastq_filenames.txt
