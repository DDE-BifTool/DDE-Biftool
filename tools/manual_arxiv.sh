#!/bin/bash
#
# Create ddebiftool_manual.zip in current folder for uploading to arxiv
#
# $Id: manual_arxiv.sh 98 2015-01-11 21:44:46Z jansieber $
#
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 version"
    echo " creates ddebiftool_manual.zip in current folder,"
    echo " inserting version for all mentions of \version in tex file"
    exit
fi
set -e 
version=$1
curdir=`pwd`
tooldir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
msrc="neuron_sys_deri.m  sd_dtau.m  sd_rhs.m  sd_tau.m"
tex="manual.tex manual.bbl"
figs="fig/*.pdf"
license=$tooldir"/license.txt"
arxdir="ddebiftool_manual"
zip=$curdir"/ddebiftool_manual.zip"
cd $curdir
rm -rf $arxdir $zip
mkdir -p $arxdir $arxdir"/fig"
cd $tooldir"/../manual"
cp -p $msrc  $tex $license "$curdir/$arxdir"
cp -p $figs "$curdir/$arxdir/fig"
cd "$curdir/$arxdir"
# insert date into license
year=`date +%G`
sed -i -- "s/|year|/$year/g" license.txt
echo '\newcommand{\version}{'$version'}' >version.tex
sed -i -- 's/%\\pdfoutput=1/\\pdfoutput=1/g' manual.tex
cd ..
zip -r $zip $arxdir
