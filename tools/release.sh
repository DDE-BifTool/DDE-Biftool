#!/bin/bash
#
# $Id: release.sh 346 2019-05-13 05:41:50Z jansieber $
#
set -e 
#
# Check arguments
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 basefolder codefolder releasefolder version [testprog]"
    echo " * calls svn export putting export of basefolder/codefolder into exportdir=basefolder/tags/dde_biftool_r{revision}"
    echo " * copies export to maindir=releasefolder/dde_biftool_v{version}"
    echo " * compiles manual, cover and other docs and copies them into newly created folder doc,"
    echo " * replaces Id lines and updates (c) lines with {version}(commit)"
    echo " * insert {version} and {revision} into Readme.html and creates Readme.txt from Readme.html"
    echo " * if testprog is given applies testprog to all demos (in temporary folder test)"
    echo "    possible choices for testprog: {matlab octave}"
    echo " * removes folders not intended for distribution"
    echo " * zips maindir into dde_biftool_v{version}.zip"
    echo " "
    echo " used programs:"
    echo " bash, svn, pdflatex, bibtex, python (tested with 2.6), unix2dos, html2text"
    echo " for testing {matlab, octave}"
    exit
fi
curdir=`pwd`
base=`cd $1;pwd`
tag=`cd $base/tags;pwd`
codefolder=$2
cd $curdir
releasebase=`cd $3;pwd`
version=$4
# Redirect stdout ( > ) into a named pipe ( >() ) running "tee"
exec > >(tee "$releasebase/logfile_v$version.txt")
if [[ $# -lt 5 ]]; then
    dotest=0
else
    dotest=1
    cmd=$5
fi
# obtain current revision number
cd $base
revision=`svn -R info | awk -- '/Revision/{print $2}' | sort -n | tail -1`
# generate names for folders and files
name="dde_biftool_v"$version
exportdir=$tag"/"$name"_r"$revision
destdir="$name"
destdir=$releasebase/$destdir
docdir=$destdir/doc
license=$destdir/tools/license.txt
zip=$name".zip"
files="*/*.*"
if [[ ! ( -e $exportdir ) ]]; then
    svn export "^/"$codefolder  $exportdir --native-eol CRLF
fi
rm -rf $destdir
cp -urp $exportdir $destdir
#
# set year in license.txt
cd $destdir
year=`date +%G`
sed -i -- "s/|year|/$year/g" $license
#
# compile manuals
mkdir -p $docdir
cd $destdir/manual
cp -p $license ./license.txt
tex="manual Changes-v3"
echo '\newcommand{\version}{'$version'}' >version.tex

for x in $tex; do
    pdflatex $x && pdflatex $x && pdflatex $x && \
    bibtex $x && pdflatex $x && pdflatex $x && mv $x".pdf" $docdir
done
mv Addendum_Manual_DDE-BIFTOOL_2_03.pdf $docdir

nmfm="nmfm_extension_manual"
cd $destdir/manual/$nmfm
tex="nmfm_extension_description"
pdflatex $tex && pdflatex $tex && pdflatex $tex && \
    bibtex $tex && pdflatex $tex && pdflatex $tex && mv $tex".pdf" $docdir

cd $destdir

#
# set (c) line in all m and html files that have (c) or Id
python $destdir/tools/c_insert.py $destdir $version
#
# insert license and set version in readme
nr=`awk -- '/\|license\|/{print NR}' Readme.html`
awk -- "NR<$nr"'{print $0}' Readme.html >tmp.txt
awk -- 'BEGIN{RS="\r\n"}{if(NF==0)print "<br><br>";else print $0}' $license >>tmp.txt
awk -- "NR>$nr"'{print $0}' Readme.html >>tmp.txt
sed  -- "s/|version|/$version/g" tmp.txt > Readme.html
rm -f tmp.txt
html2text Readme.html >Readme.txt
unix2dos Readme.txt

#
# if testing required perform tests of demos
if [[ $dotest -eq 1 ]]; then
    mkdir $destdir/test
    python $destdir/tools/test_demos.py $destdir $cmd
fi
#
# remove folders not intended for distribution (except tools)
rm -rf  manual FilesChangedAndAdded_V203 system test development
#
# create index.html in folders which don't have them
i=0
for i in `find "$destdir" -mindepth 1 -type d`; do
    index=$i"/index.html"
    if [[ ! ( -f $index ) ]]; then
	echo creating $index
	bash $destdir/tools/create_index_html.sh $destdir/tools/index_template.html $version $i
    fi
done
rm -rf tools etc

# zip
cd $releasebase
rm -rf $zip
zip -r $zip $name
