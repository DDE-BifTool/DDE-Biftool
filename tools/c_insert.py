#!/usr/bin/python
#
# $Id: c_insert.py 91 2015-01-11 17:18:05Z jansieber $
#
import os, sys
import re, tempfile
usage='''\
Usage: python {0:s} rootdir version [date]
Includes version and (c) info line in all files. 
If (c) line is present but no Id line then (c) line is left untouched.
If (c) line is absent but Id line is present, 
    then Id line is replaced with (c) line using argument version
If (c) line is present and Id line is present, then version and date of (c)
    line are overwritten with version argument(commit number from Id line) 
    and date from Id line. Id line is removed.
'''

#matching patterns
filetypes=['.*[.]m', '.*[.]html']
cbif_re=r'[(]c[)]\W+DDE-BIFTOOL'
id_re=r'[$]Id[0-9.a-zA-Z_()]*:\W*[0-9.a-zA-Z_()]*'

# convert matched regexp to sortable date integer
getdate=lambda x, order:int(x[order[2]])*10000+int(x[order[1]])*100+int(x[order[0]])

# main function performing all steps
def main():
    nargs=len(sys.argv)
    if nargs<3:
        print usage.format(sys.argv[0])
        print 'nargs:', nargs
        print 'argv:', sys.argv
        return -1
    version=sys.argv[2]
    rootdir=sys.argv[1]
    files=filelist(rootdir, filetypes)
    for k in range(len(files)):
        f=open(files[k], 'r')
        text=f.readlines()
        f.close()
        c_info=c_from_file(text)
        repl=replace(text, c_info, version)
        text=repl['text']
        f=open(files[k], 'w')
        f.writelines(text)
        f.close()
        if repl['line']:
            line=text[repl['line']]
        else:
            line=''
        print files[k], line
    return 0

# extract version, commit and date from line
def v_info(line):
    delim=r'[^<> ,/-]+'
    c_info={}
    if re.search(cbif_re, line):
        c_info['kind']='(c)'
        s_list=re.findall(delim,line)
        # find 'v.'
        pos=[i for i in range(len(s_list))if re.match('v', s_list[i])][0]
        version=s_list[pos+1]
        if re.search(r'[(]', version):
            version=re.findall('[^()]+', version)
            c_info['ver']=version[0]
            c_info['commit']=int(version[1])
        else:
            c_info['ver']=version
            c_info['commit']=None
        c_info['date']=getdate(s_list[pos+2:pos+5], [0, 1, 2])
        return c_info
    if re.search(id_re, line):
        c_info['kind']='Id'
        s_list=re.findall(delim,line)
        # find '$Id:'
        pos=[i for i in range(len(s_list))if re.match('[$]Id', s_list[i])][0]
        c_info['commit']=int(s_list[pos+2])
        if len(s_list)>pos+6:
            c_info['date']=getdate(s_list[pos+3:pos+6], [2, 1, 0])
        else:
            c_info['date']=None
        return c_info
    return None

# extract (c) and Id lines of file
def c_from_file(text):
    c_list=[v_info(k) for k in text]
    c_info={'(c)':None, 'Id':None}
    c_ind=filter(lambda x:c_list[x] is not None and c_list[x]['kind']=='(c)', range(len(c_list)))
    id_ind=filter(lambda x:c_list[x] is not None and c_list[x]['kind']=='Id', range(len(c_list)))
    if c_ind:
        c_info['(c)']=c_list[c_ind[0]]
        c_info['(c)']['line']=c_ind[0]
    if id_ind:
        c_info['Id']=c_list[id_ind[0]]
        c_info['Id']['line']=id_ind[0]
    return c_info

# replace everything between patterns[0] and patterns[1] with repl
def replace_c_id(line, patterns, repl=None):
    idstart=re.search(patterns[0], line)
    if not repl:
        repl=''
    if idstart:
        p0=idstart.start(0)
        plen=re.search(patterns[1], line[p0+1:]).end(0)
        line=line[:p0]+repl+line[p0+plen+1:]
    return line

# replace file version with new version if Id present
# and remove Id line
def replace(text, c_info, version, date=None):
    if not c_info['Id'] or c_info['Id']['commit']<0: # no id line -> return
        return {'text': text, 'line':None}
    id=c_info['Id']
    idpat=[id_re, r'[$]']
    cpat=[cbif_re, r'[/][0-9][0-9][0-9][0-9]']
    if not c_info['(c)']: # no (c) line (yet) -> create one in place of Id line
        repl=c_create(version, id['date'], id['commit'])
        text[id['line']]=replace_c_id(text[id['line']], idpat, repl)
        return {'text': text, 'line': id['line']}
    else:
        # otherwise replace (c) with info of Id line and version argument
        cbif=c_info['(c)']
        repl=c_create(version, id['date'], id['commit'])
        text[cbif['line']]=replace_c_id(text[cbif['line']], cpat, repl)
        text[id['line']]=replace_c_id(text[id['line']], idpat)    
        return {'text': text, 'line': cbif['line']}
    
# recursive list of files with names matching fpatterns in path
def filelist(path, fpatterns):
    '''
    returns a list of all files in path (recursively) with filename pattern fpattern and 
    containing a line matching regex
    list elements have format [filename_with_path, number_of_first_matching_line]
    '''
    fp_re=[re.compile(k) for k in fpatterns]
    res = []
    for root, dirs, fnames in os.walk(path):
        for fname in fnames:
            for fp in fp_re:
                if fp.match(fname) is None:
                    continue
                res.append(os.path.join(root,fname));
    return res

# create (c) string with version (as string), commit number (if not None) and date
def c_create(version,  date, commit=None):
    year=date/10000
    month=(date-year*10000)/100
    day=date-year*10000-month*100
    if commit:
        s='(c) DDE-BIFTOOL v. {0:s}({1:d}), {2:02d}/{3:02d}/{4:04d}'.format(version, commit, day, month, year)
    else:
        s='(c) DDE-BIFTOOL v. {0:s}, {2:02d}/{3:02d}/{4:04d}'.format(version, commit, day, month, year)
    return s

# scipt body
if __name__ == '__main__':
    sys.exit(main())
