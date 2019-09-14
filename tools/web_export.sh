#!/bin/bash
#
# $Id: web_export.sh 95 2015-01-11 21:14:27Z jansieber $
#
set -e 
#
# Check arguments
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 username folder"
    echo " rsync folder contents to sourceforge web host htdocs folder"
    exit
fi
username=$1
folder=$2
rsync -vturp -e "ssh -l $username" $folder"/" web.sourceforge.net:/home/project-web/ddebiftool/htdocs/
