#!/usr/bin/python
#
# $Id$
#
usage='''
Usage: {0:s} basefolder command [test_names]
This runs parametrics in basefolder/tests/ by copying each subfolder to
basefolder/test_run/ changing into the subfolder and calling
if command is "octave"
octave --eval [mainscript]
if command is "matlab"
matlab -nosplash -nodesktop -r [mainscript]

If test_names are omitted all tests are run one after another.
Valid test_names are 
gen_sym_minimal_demo
minimal_demo_2d
gen_sym_RoseHindmarsh
rose_hindmarsh
'''
import sys, os, subprocess, shutil

sep='\n=====================\n'


# main function performing all steps
def main():
    nargs=len(sys.argv)
    if nargs<3:
        print usage.format(sys.argv[0])
        print 'nargs:', nargs
        print 'argv:', sys.argv
        return -1
    rootdir=sys.argv[1]
    tdir=rootdir+'/test_run/'
    ddir=rootdir+'/tests/'
    sys.path.append(os.path.abspath(ddir))
    testdescriptions=__import__('descriptions')
    tests=testdescriptions.tests
    cmd=testdescriptions.commands[sys.argv[2]]
    frame=testdescriptions.frame
    environment=testdescriptions.environment
    bdir=os.getcwd()
    try:
        shutil.rmtree(tdir)
    except:
        None
    os.mkdir(tdir)
    retval={}
    demosel=sys.argv[3:]
    # find out which tests to perform over which ranges
    testnames=tests.keys()
    demos=[testnames[k] for k in range(len(tests))]
    ranges=[tests[k]['args'] for k in tests.keys()]
    if len(demosel)>0:
        demolst=[]
        rangelist=[]
        for ind in range(len(demosel)):
            if demosel[ind] in demos:
                demolst.append(demosel[ind])
                rangelist.append([])
            else:
                rangelist[-1].append(int(demosel[ind]))
        for ind in range(len(rangelist)):
            if len(rangelist[ind])==0:
                rangelist[ind]=tests[demolst[ind]]['args']
    else:
        # do all demos
        demolst=demos
        rangelist=ranges
    # copy tools folder over
    shutil.copytree(ddir+'tools', tdir+'tools')
    # create results folder
    results=tdir+'results/'
    try:
        shutil.rmtree(results)
    except:
        None
    os.mkdir(results)
    # run each test
    for q in range(len(demolst)):
        testname=demolst[q]
        test=tests[testname]
        folder=test['folder']
        shutil.copytree(ddir+folder, tdir+folder)
        try:
            os.chdir(tdir+folder)
        except:
            print tdir+folder+' not found'
            continue
        print sep+'testing {0:s} with value {1:s}'.format(testname,
                           ','.join([str(k) for k in rangelist[q]]))
        retval[q]=list()
        for k in range(len(rangelist[q])):
            prog=frame(test['init'],test['testnum'],rangelist[q][k],test['file'])
            call=list(cmd)
            call.append(prog)
            fname=testname+'-{0:02d}.out'.format(rangelist[q][k])
            fname=os.path.join('..','results',fname)
            f=open(fname,'w',1)
            print >>f, call
            ret=subprocess.call(call,env=environment,stdout=f,stderr=f)
            out=testname+' test {0:d} with value {1:d}'.format(k,rangelist[q][k])
            if not ret:
                out=out+' ok'+sep
            else:
                out=out+' failed'+sep
            print >>f, out
            retval[q].append(out)
            f.close()
        os.chdir(bdir)
    print sep
    for t in retval.keys():
        print t
        for res in retval[t]:
            print res
        print sep
# scipt body
if __name__ == '__main__':
    sys.exit(main())
