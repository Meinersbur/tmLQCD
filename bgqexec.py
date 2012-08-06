 #! /usr/bin/env python
 
import sys
import operator
import numbthy
import getopt
import os
import argparse
import shutil
import fnmatch


# globals
jobdir = ''
jobname = ''
jobscript = ''

workfs = os.path.realpath(os.path.expandvars('${WORK}'))



def isInDir(dir, ancestor):
    return ancestor.startswith(dir)

def genJobname(i):
    if i==0:
        return jobname
    else:
        return jobname + '_' + i
    
    
    
def instTemplate(file):
    fr=open(file, 'r')
    template=fr.read()
    fr.close()
    
    # do replacements
    
    sext = os.path.splitext(file);
    newfile = sext[0]
    fw.open(newfile, 'w')
    fw.write(template)
    fw.close()
    
def prepareJob():
    # replace the variable in all .template files
    for root, dirs, files in os.walk(jobdir):
        for file in fnmatch.filter(files, '*.template'):
            absfile = os.path.join(root,file)
            instTemplate(absfile)

    # Execute the pre-submit script; it should also submit the job
    os.system('job.sh')


def submitJob():
    pass



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--jobscript', help="A script that prepares and then submits the job")
    args = parser.parse_args()
    
    global jobscript
    global jobdir
    if args.jobscript:
        jobscript = os.path.realpath(args.jobscript)
    else:
        jobscript = os.path.join(os.path.realpath(os.curdir), 'job.sh')
    jobdir = os.path.split(jobscript)[0]
    
    if not isInDir(workfs, jobdir):
        # copy everything into the workfs before continuing
        
        # find a unique name
        i = 0
        while True:
            newjobname = genJobname(i)
            newjobdir = os.path.join(workfs, newjobname)
            if not os.path.exists(newjobdir):
                break
            i += 1
            
        # Copy everything to that new dir
        shutil.copytree(jobdir, newjobdir,  symlinks=True)
        
        # start the script in the new directory
        jobdir = newjobdir
        jobname = newjobname
        jobscript = os.path.relpath(jobscript, start)
        
    os.chdir(jobdir)
    prepareJob()
    startJob()
        

if __name__ == "__main__":
	main()
