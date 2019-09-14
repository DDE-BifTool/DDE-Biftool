from os import environ
tests={
    'gen_sym_minimal_demo': {
        'folder':'minimal_demo',
        'file':'gen_sym_minimal_demo',
        'testnum':'{0:d};',
        'args':range(1,2),
        'init':'pkg load symbolic'},
    'minimal_demo': {
        'folder':'minimal_demo',
        'file': 'minimal_demo_2d',
        'testnum':'indfuncs={0:d};',
        'args':range(1,7),
        'init':''},
    'gen_sym_RoseHindmarsh': {
        'folder':'rose_hindmarsh',
        'file': 'gen_sym_RoseHindmarsh',
        'testnum':'%{0:d};',
        'args': range(1,2),
        'init':'pkg load symbolic'},
    'rose_hindmarsh': {
        'folder': 'rose_hindmarsh',
        'file': 'rose_hindmarsh',
        'testnum': '{0:d};',
        'args': range(1,7),
        'init': ''}
    }
commands={'octave':['octave-cli', '--eval'],
    'octave4':['/usr/local/octave/4.0.0/bin/octave-cli', '--eval'],
    'matlab': ['matlab', '-nodesktop','-nosplash', '-r']}

matlabframe='try;{0:s};end;more off;try;{1:s};{2:s};exit(0);catch;exit(-1);end'

frame=lambda prep,val,rgi,prog: matlabframe.format(prep,val.format(rgi),prog)

environment=environ.copy()
environment['LIBGL_ALWAYS_SOFTWARE']='1'
