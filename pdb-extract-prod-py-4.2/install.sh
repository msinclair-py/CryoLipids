#!/usr/bin/env bash

## For Bourne shell users:

function highlight_text {
    echo -ne "\e[93m$1\e[0m"
}

function download_file {
    local url=$1
    local filename="${url##*/}"

    if [[ ! -f "$filename" ]]; then
        curl $url -O
    else
        echo "[*] file $(highlight_text $filename) already exists, skipping download"
    fi
}

if [ ! -d packages ]
then
echo "[*] Make packages folder"
mkdir packages
else
echo "[*] Packages folder exists"
fi
echo "[*] Change directory to Maxit installation folder"
cd packages

maxit_name="maxit-v11.100-prod-src"
maxit_filename=${maxit_name}.tar.gz
url_maxit=https://sw-tools.rcsb.org/apps/MAXIT/$maxit_filename
echo "[*] Downloading Maxit suite from $(highlight_text $url_maxit)"
download_file $url_maxit
echo "[*] File $(highlight_text $maxit_filename) downloaded"
echo "[*] Decompressing $(highlight_text $maxit_filename) starts"
tar -zxvf $maxit_filename
echo "[*] Decompressing $(highlight_text $maxit_filename) finished"
echo "[*] Maxit installation starts"
cd $maxit_name
make binary

cd ../..
echo "[*] Maxit installation finished"
