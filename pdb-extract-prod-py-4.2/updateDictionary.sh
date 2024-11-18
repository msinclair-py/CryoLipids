#!/usr/bin/env bash

## For Bourne shell users:
## Must run "source setup.sh" before executing this script to set up enviroment variables.

function highlight_text {
    echo -ne "\e[93m$1\e[0m"
}

function download_file {
    local url=$1
    local filename="${url##*/}"

    if [[ ! -f "$filename" ]]; then
        curl $url -O
    else
        echo "[*] file $(highlight_text $filename) already exists, removing old file"
		rm $filename
		wget $url
    fi
}

if [ -z "$PDB_EXTRACT_PY" ] 
then
    echo "[*] Please run source setup.sh before executing this script"
	exit 1
else 
	echo "[*] Change directory to dictionary folder of"
	cd $PDB_EXTRACT_PY/data/dictionary
	pwd
fi

dict_name="mmcif_pdbx_v5_next.dic.gz"
url_dict=https://mmcif.wwpdb.org/dictionaries/ascii/$dict_name
echo "[*] Downloading latest mmCIF dictionary from $(highlight_text $url_dict)"
download_file $url_dict
echo "[*] Decompressing $(highlight_text $dict_name)"
gunzip -f $dict_name
echo "[*] mmCIF dictionary file has been updated"

