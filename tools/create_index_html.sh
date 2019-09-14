#!/bin/bash
#
# $Id: create_index_html.sh 91 2015-01-11 17:18:05Z jansieber $
#
# Modified from 
# http://stackoverflow.com/questions/21395159/shell-script-to-create-a-static-html-directory-listing
set -e 
if [[ $# -lt 3 ]]; then
    echo "Usage: $0 template version folder"
    echo "  Creates index.html in {folder} containing a bullet list of file names"
    echo "  by inserting the files names into the html frame template"
    exit
fi
template=$1
version=$2
ROOT=$3
HTTP="/"
OUTPUT=$3"/index.html" 
startpat='<!-- start_insert -->'
sed "s/|version|/$version/" $template > $OUTPUT
fullpath=`readlink -f "$ROOT"`
path=`basename "$fullpath"`
sed -i "s/$startpat/<h2>$path<\/h2>\n<ul>\n$startpat/" $OUTPUT
i=0
for i in `find "$ROOT" -maxdepth 1 -mindepth 1 -type d| sort`; do
    file=`basename "$i"`
    sed -i  "s/$startpat/<li><a href=\"$file\/index.html\">$file<\/a><\/li>\n$startpat/1"  $OUTPUT
done
sed -i "s/$startpat/<hr>\n$startpat/1"  $OUTPUT
for i in `find "$ROOT" -maxdepth 1 -mindepth 1 -type f| sort`; do
    file=`basename "$i"`
    sed -i  "s/$startpat/<li><a href=\"$file\">$file<\/a><\/li>\n$startpat/1"  $OUTPUT
done
sed -i "s/$startpat/<\/ul>\n$startpat/1" $OUTPUT
