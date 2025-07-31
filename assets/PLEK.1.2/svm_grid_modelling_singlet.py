#!/usr/bin/env python
# updated on June 26, 2014
# Related to: PLEK_howto_generate_scripts.doc, PLEK_generate_scripts.R, 
# Called by: PLEKModelling.py

__all__ = ['find_parameters']

import os, sys, traceback, getpass, time, re, datetime,subprocess


class GridOption:
	def __init__(self, dataset_pathname, options):
		dirname = os.path.dirname(__file__)
		if sys.platform != 'win32':
			self.svmtrain_pathname = os.path.join(dirname, './svm-train')
			self.gnuplot_pathname = '/usr/bin/gnuplot'
		else:
			# example for windows
			self.svmtrain_pathname = os.path.join(dirname, r'..\windows\svm-train.exe')
			# svmtrain_pathname = r'c:\Program Files\libsvm\windows\svm-train.exe'
			self.gnuplot_pathname = r'c:\tmp\gnuplot\binary\pgnuplot.exe'
		self.fold = 5 # cross validation
		self.c_begin, self.c_end, self.c_step = -3,  7,  2 # origin: -5,  15,  2, by aimin
		self.g_begin, self.g_end, self.g_step =  1, -7, -2 # origin:  3, -15, -2, by aimin
		self.grid_with_c, self.grid_with_g = True, True
		self.dataset_pathname = dataset_pathname #data, input
		self.dataset_title = os.path.split(dataset_pathname)[1]
		self.out_pathname = '{0}.grid.out'.format(self.dataset_title) #output ,  grid.out by aimin
		self.png_pathname = '{0}.grid.png'.format(self.dataset_title) #output ,  grid.png by aimin		
		self.pass_through_string = ' '
		self.thread=10
		self.resume_pathname = None
		self.logfile=""
		self.parse_options(options)
		
		self.err_pathname = '{0}.grid.err'.format(self.dataset_title) #output ,  grid.png by aimin
		err_file = open(self.err_pathname, 'w')		
		#nowTime = datetime.datetime.now()
		#strStartTime = time.strftime('%Y-%m-%d %H:%M:%S', nowTime)
		#err_file.write(strStartTime)
		err_file.flush()
		err_file.close()
		

	def parse_options(self, options):
		if type(options) == str:
			options = options.split()
		i = 0
		pass_through_options = []
		
		while i < len(options):
			if options[i] == '-log2c': # parameter c
				i = i + 1
				if options[i] == 'null':
					self.grid_with_c = False
				else:
					self.c_begin, self.c_end, self.c_step = map(float,options[i].split(','))
			elif options[i] == '-log2g': # parameter g
				i = i + 1
				if options[i] == 'null':
					self.grid_with_g = False
				else:
					self.g_begin, self.g_end, self.g_step = map(float,options[i].split(','))
			elif options[i] == '-v': # parameter v, n-fold cross validation
				i = i + 1
				self.fold = options[i]
			elif options[i] in ('-c','-g'):
				raise ValueError('Use -log2c and -log2g.')
			elif options[i] == '-svmtrain': # parameter svmtrain
				i = i + 1
				self.svmtrain_pathname = options[i]
			elif options[i] == '-gnuplot': # parameter gnuplot
				i = i + 1
				if options[i] == 'null':
					self.gnuplot_pathname = None
				else:	
					self.gnuplot_pathname = options[i]
			elif options[i] == '-out': # parameter out
				i = i + 1
				if options[i] == 'null':
					self.out_pathname = None
				else:
					self.out_pathname = options[i]
			elif options[i] == '-png': # parameter png
				i = i + 1
				self.png_pathname = options[i]
			elif options[i] == '-logfile': # logfile
				i = i + 1
				self.logfile = options[i]
			elif options[i] == '-thread': # parameter thread
				i = i + 1
				self.thread = int(options[i])
				nr_local_worker= int(options[i])
			elif options[i] == '-resume': # parameter resume
				if i == (len(options)-1) or options[i+1].startswith('-'):
					self.resume_pathname = self.dataset_title + '.grid.out' # by aimin
				else:
					i = i + 1
					self.resume_pathname = options[i]
			else:
				pass_through_options.append(options[i])
			i = i + 1

		self.pass_through_string = ' '.join(pass_through_options)
		if not os.path.exists(self.svmtrain_pathname):
			raise IOError('svm-train executable not found')
		if not os.path.exists(self.dataset_pathname):
			raise IOError('dataset not found')
		if self.resume_pathname and not os.path.exists(self.resume_pathname):
			raise IOError('file for resumption not found')
		if not self.grid_with_c and not self.grid_with_g:
			raise ValueError('-log2c and -log2g should not be null simultaneously')
		if self.gnuplot_pathname and not os.path.exists(self.gnuplot_pathname):
			sys.stderr.write('gnuplot executable not found\n')
			self.gnuplot_pathname = None
		if os.path.isfile(self.out_pathname):
			os.remove(self.out_pathname)


def find_parameters(dataset_pathname, options=''):
	options = GridOption(dataset_pathname, options);

	cmdline = options.svmtrain_pathname
	if options.grid_with_c: 
		cmdline += ' -c {0} '.format(2.0**options.c_begin)
	if options.grid_with_g: 
		cmdline += ' -g {0} '.format(2.0**options.g_begin)
	cmdline += ' -v {0} {1} {2} '.format\
		(options.fold, options.pass_through_string, options.dataset_pathname)
		
	p = subprocess.Popen(cmdline, stdout=subprocess.PIPE, shell=True)
	(output, err) = p.communicate()
	p_status = p.wait()
	#print "\t\t", cmdline
	if p_status!=0: # error
		print "\t\t", err
		print "Command output : ", output
	#print "Command output : ", output
	#print "Command exit status/return code : ", p_status

	# svm-train -c 0.25 -g -0.125 -v 2 predicted_temp_2_scaled
	# Cross Validation Accuracy = 88.2%
	rate=0.0
	for line in output.split(os.linesep):
		if str(line).find('Cross') != -1:
			rate=float(line.split()[-1][0:-1])
	
	#
	if p_status!=0:
		options.c_begin=-1
		options.g_begin=-1
		rate=-1
	
	#	2.0 -3.0 91.5 (best c=4.0, g=0.125, rate=91.5)
	stdout_str = '	{0} {1} {2} (best '.format(options.c_begin, options.g_begin, rate) # Aimin, output  
	stdout_str += 'c={0}, '.format(2.0**options.c_begin)
	stdout_str += 'g={0}, '.format(2.0**options.g_begin)
	stdout_str += 'rate={0})'.format(rate)
		
	print(stdout_str + ' [{0}]'.format(time.strftime('%Y-%m-%d %H:%M:%S'))) #output to stdout
	
	if options.out_pathname:
		result_file = open(options.out_pathname, 'w')
		result_file.write(stdout_str + "\n")
		result_file.flush()
		result_file.close()
		
	
if __name__ == '__main__':

	def exit_with_help():
		print("""\
Usage: grid.py [grid_options] [svm_options] dataset

grid_options :
-log2c {begin,end,step | "null"} : set the range of c (default -5,15,2)
    begin,end,step -- c_range = 2^{begin,...,begin+k*step,...,end}
    "null"         -- do not grid with c
-log2g {begin,end,step | "null"} : set the range of g (default 3,-15,-2)
    begin,end,step -- g_range = 2^{begin,...,begin+k*step,...,end}
    "null"         -- do not grid with g
-v n : n-fold cross validation (default 5)
-svmtrain pathname : set svm executable path and name
-gnuplot {pathname | "null"} :
    pathname -- set gnuplot executable path and name
    "null"   -- do not plot 
-out {pathname | "null"} : (default dataset.out)
    pathname -- set output file path and name
    "null"   -- do not output file
-png pathname : set graphic output file path and name (default dataset.png)
-resume [pathname] : resume the grid task using an existing output file (default pathname is dataset.out)
    This is experimental. Try this option only if some parameters have been checked for the SAME data.

svm_options : additional options for svm-train""")
		sys.exit(1)
	
	if len(sys.argv) < 2:
		exit_with_help()
	dataset_pathname = sys.argv[-1]
	options = sys.argv[1:-1]
	try:        
		find_parameters(dataset_pathname, options)
	except (IOError,ValueError) as e:
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write('Try "grid.py" for more information.\n')
		sys.exit(1)
