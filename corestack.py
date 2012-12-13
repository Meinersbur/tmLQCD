#! /usr/bin/env python
import os
import sys
import glob
import os.path
import re
import operator
import shutil

prog = re.compile(r'^Program\s*:\s*(?P<program>[\w/.]+)$');
reason = re.compile(r'^\*\*\*FAULT');
failaddr = re.compile(r'^While executing instruction at[^0x]*(?P<addr>\w+)$');

stackstart = re.compile(r'^\+\+\+STACK$');
stackframe = re.compile(r'^\w+\s*(?P<pc>\w+)$');



exe = None

def addr2line(addr):
	#print 'addr2line -e ' + exe + ' ' + addr
	sys.stdout.flush()
	os.system('addr2line -e ' + exe + ' ' + addr)
	sys.stdout.flush()


def process_file(filename):
	global exe
	f = open(filename, 'r')
	folder = os.path.split(filename)[0]
	for line in f:
		m = prog.match(line)
		if m:
			exe = m.group('program')
			exe = os.path.join(folder, exe)
			exe = os.path.realpath(exe)
			print 'Program is', exe
			
		m = reason.match(line);
		if m:
			print line
			
		m = failaddr.match(line);
		if m:
			addr = m.group('addr')
			sys.stdout.write('Fault at ')
			sys.stdout.flush()
			addr2line(addr)
			print ''

		m = stackstart.match(line)
		if m:
			print ''
			print 'Thread Stack: '
		
		m = stackframe.match(line)
		if m:
			addr = m.group('pc')
			if addr != '0000000000000000':
				addr2line('0x' + addr)

	
def main(argv):
	process_file(argv[1])


if __name__ == "__main__":
    main(sys.argv)
