from dateutil.parser import *
from datetime import *
import commands
import sys
import os

now = parse(commands.getoutput("date"))
today = now.date()
this_minute = now.time()

date_and_time = str(today) + '_' + str(this_minute)
computer_name = 'aravali2'
size = 8
Hs = [3000]
for H in Hs:
	cmd = 'scp -r current-GLL belliappa@' + computer_name + ':; ssh belliappa@' + computer_name + ' "cd current-GLL; rm -rf *.dat; python just-automate.py 1 ' + str(size) + ' ' + str(H) + '"; mkdir ' + date_and_time + '_' + str(size) + '_' + str(H) + '; mkdir ' + date_and_time + str(size) + '_' + str(H) + '/data; scp -r belliappa@' + computer_name + ':current-GLL/*.dat ./' + date_and_time + str(size) + '_' + str(H) + '/data'
	print cmd
	os.system(cmd)
	cmd = 'cp current-GLL/just-automate.py ' + date_and_time + str(size) + '_' + str(H) + '/data; cp current-GLL/anal-hist.c ' + date_and_time + str(size) + '_' + str(H) + '/data; cd ' + date_and_time + str(size) + '_' + str(H) + '/data; python just-automate.py 2 ' + str(size) + ' ' + str(H) + '; cd ..; mkdir plots; mv data/*.ps plots'
	print cmd
	os.system(cmd)

