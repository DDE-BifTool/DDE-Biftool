#!/usr/bin/python
#
# $Id: test_demos.py 176 2017-03-13 00:25:33Z jansieber $
#
usage='''
Usage: {0:s} basefolder command [demo_names]
This tests demos in basefolder/demos/ by copying each subfolder to
basefolder/test/ changing into the subfolder and calling
if command is "octave"
octave --eval [mainscript]
if command is "matlab"
matlab -nosplash -nodesktop -r [mainscript]

If demo_names are ommitted all demos are test one after another.
Valid demo_names are 
neuron
sd_demo
minimal_demo
nmfm_demo
Mackey-Glass
nested
rotsym_demo
humphriesetal
'''
import sys, os, subprocess, shutil
demos=[['neuron', 'rundemo'], 
       ['sd_demo', 'rundemo'], 
       ['minimal_demo', 'rundemo'], 
       ['nmfm_demo', 'nmfm_demo'], 
       ['hom_demo', 'hom_demo'], 
       ['Mackey-Glass', 'MackeyGlass_demo'], 
       ['nested', 'nested_demo'], 
       ['rotsym_demo', 'rotsym_demo'], 
       ['humphriesetal', 'rundemo'],
       ['phase_oscillator','phase_oscillator'],
       ['Holling-Tanner','HollingTanner_demo'],
       ['cusp','cusp_demo']]
sep='\n=====================\n'
frame=lambda x: 'try;{0:s};exit(0);catch;exit(-1);end'.format(x)
commands={'octave':['octave', '--eval'],
    'matlab': ['matlab', '-nodesktop','-nosplash', '-r']}
# main function performing all steps
def main():
    nargs=len(sys.argv)
    if nargs<3:
        print usage.format(sys.argv[0])
        print 'nargs:', nargs
        print 'argv:', sys.argv
        return -1
    rootdir=sys.argv[1]
    cmd=commands[sys.argv[2]]
    tdir=rootdir+'/test/'
    ddir=rootdir+'/demos/'
    bdir=os.getcwd()
    try:
        shutil.rmtree(tdir)
    except:
        None
    os.mkdir(tdir)
    retval={}
    demosel=sys.argv[3:]
    if len(demosel)>0:
        demolst=[]
        for ind in range(len(demosel)):
            demolst.append([k for k in range(len(demos)) if demos[k][0]==demosel[ind]][0])
    else:
        demolst=range(len(demos))
    for q in demolst:
        p=demos[q]
        shutil.copytree(ddir+p[0], tdir+p[0])
        try:
            os.chdir(tdir+p[0])
        except:
            print tdir+p[0]+' not found'
            continue
        print sep+'testing '+p[0]
        prog=frame(p[1])
        call=list(cmd)
        call.append(prog)
        retval[p[0]]=subprocess.call(call)
        if not retval[p[0]]:
            print p[0]+' ok'+sep
            retval[p[0]]='ok'
        else:
            print p[0]+' failed'+sep
            retval[p[0]]='failed'
        os.chdir(bdir)
        shutil.rmtree(tdir+p[0])
    print sep
    print retval

# scipt body
if __name__ == '__main__':
    sys.exit(main())
