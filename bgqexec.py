#! /usr/bin/env python
 
import sys
import operator
import getopt
import os
import optparse
import shutil
import fnmatch


# globals
jobdir = ''
jobname = ''
jobscript = ''

workfs = os.path.realpath(os.path.expandvars('${WORK}'))



def isInDir(dir, ancestor):
    return os.path.realpath(ancestor).startswith(os.path.realpath(dir))

def genJobname(i):
    return jobname + '_' + ("%03d"%i)
    
    
    
def instTemplate(file):
    sext = os.path.splitext(file);
    newfile = sext[0]
    
    print "Configure " + file + " to " + newfile + " ..."
    fr = open(file, 'r')
    template = fr.read()
    fr.close()
    
    # do replacements
    
    fw = open(newfile, 'w')
    fw.write(template)
    fw.close()
    
    shutil.copystat(file, newfile)
    
    
def prepareJob():
    # replace the variable in all .template files
    print "Configuration..."
    for root, dirs, files in os.walk(jobdir):
        for file in fnmatch.filter(files, '*.template'):
            absfile = os.path.join(root, file)
            instTemplate(absfile)

    # Execute the pre-submit script; it should also submit the job
    print "Run the jobscript"
    os.system(jobscript)



   



def main():
    parser = optparse.OptionParser()
    parser.add_option('--jobscript', help="A script that prepares and then submits the job")
    (options, args) = parser.parse_args()
    
    global jobscript
    global jobdir
    global jobname
    if options.jobscript:
        jobscript = os.path.realpath(args.jobscript)
    else:
        jobscript = os.path.join(os.path.realpath(os.curdir), 'job.sh')
    jobdir = os.path.split(jobscript)[0]
    jobname = os.path.split(jobdir)[1]
    
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
        print "Copying job files to " + newjobdir + "..."
        shutil.copytree(jobdir, newjobdir, symlinks=True)
        
        # start the script in the new directory
        jobscript = os.path.join(newjobdir, os.path.relpath(jobscript, jobdir))
        jobdir = newjobdir
        jobname = newjobname
        
    print "Jobdir: " + jobdir
    print "Jobname: " + jobname
    print "Jobscript: " + jobscript
        
    os.chdir(jobdir)
    prepareJob()
    
        

if __name__ == "__main__":
	main()
