#!/bin/bash

#########
This script loops through the row from the list file downloded from The 1000 Genomes Project and takes the ftp address from first column to wget the vcf file and download it to your defined directory.
#########

### Define the directory where you want to save the downloaded files
download_dir="download/"
output_file="stdout.out"

### Path to your TSV file
tsv_file="ftp_samples_list.tsv"

### Check if the download directory exists, if not, create it
mkdir -p "$download_dir"

### Loop through each row in the TSV file and take the first column
while IFS=$'\t' read -r ftp_address
do
    # Only use the first column (ftp_address) from the TSV
    ftp_url="ftp://$ftp_address/$filename"
    #echo "Downloading $ftp_address to $download_dir" | tee -a "$output_file"
    echo $ftp_address
    # Use wget to download the file
    wget $ftp_address -P $download_dir

done < "$tsv_file"
