## Create directional derivatives with maple
#
# If command line maple is in the path this file can be called with
# >maple maple_dirderi.mpl
#
# $Id$
#
restart:
# if the file is called inside a worksheet, uncomment the interface line
curdir:=currentdir();
#curdir:=interface(worksheetdir);
# load library for creating directional derivatives of DDE-Biftool
# functionals
read(FileTools[JoinPath]([curdir,"..","..","external_tools","ddebiftool_deri.mpl"]));
# enter desired file name for matlab function file
filename:="humphriesetal_dirderi.m";
file:=FileTools[JoinPath]([curdir,filename]);
# enter functional (for sd-ddes in full functional form, for
# constant-delay dde x[j,k] is the jth component with delay number k-1)
f0:=Vector([
-p[5]*x[1](0)-p[1]*x[1](-p[3]-p[6]*x[1](0))-p[2]*x[1](-p[4]-p[6]*x[1](0))]);
# call of derivative genrating routine
# fcnstr is text of matlab file, formula is maple expression with derivatives
fcnstr,formula:=ddebiftool_deri[sd_dirderi](f0,5):
formula;
# save string to file
fd := fopen(file, WRITE):
fprintf(fd,fcnstr):
fclose(fd):
quit

