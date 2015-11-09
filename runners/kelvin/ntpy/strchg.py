## ## ## strchg.py v.0.1
## ## ## This program simplifies and executes
## ## ## common text manipulation operations.
## ## ## Created: 07/29/2012 - KDP
## ## ## Last Edited: 02/01/2013 - KDP
import os
import numpy as np

def sed(strings, orig, new):
	"""
	Simplifies and executes the sed terminal command.

	ntpy.strchg.sed(strings, orig, new)
	Parameters
	----------
		strings : dict of type str
			A dictionary holding key:value pairs. Each key
			is original string to be replaced and each value
			is the new string to be subsituted.
		orig : str
			A string containing the original file name. If
			the file is not included in the pathway then it
			can be the absolute or relative pathway to the
			file.
		new : str
			A string containing the new file name. If the
			file is not included in the pathway then it can
			be the absolute or relative pathway to the file.
			An empty string activates the --in-place (-i)
			switch, which edits the original file without
			creating a new one.
	"""
	## Create concatenated string of command expressions
	commands = ''
	for key in strings:
		commands = commands +'-e \'s/' + str(key) + '/' + \
				   str(strings[key]).replace(r'/', r'\/') + '/g\' '
	
	## Execute the sed command for new file or --in-place
	if new == '':
		os.system('sed -i ' + commands + str(orig))
	else:
		os.system('sed ' + commands + str(orig) + ' > ' + str(new))
#-- END Sed --#


def grep(pattern, lines, orig, new, flag='-A'):
	"""
	Simplifies and executes the grep terminal command.
	Supports -A and -v flags.

	ntpy.strchg.sed(pattern, lines, orig, new, flag='-A')
	Parameters
	----------
		pattern : str
			A string of the pattern that grep should search
			for.
		lines : int
			An integer of how many lines to pull after
			pattern match.
		orig : str
			A string containing the original file name. If
			the file is not included in the pathway then it
			can be the absolute or relative pathway to the
			file.
		new : str
			A string containing the new file name. If the
			file is not included in the pathway then it can
			be the absolute or relative pathway to the file.
		flag : {'-A','-v'}, optional
			Specify grep flag. If flag is '-A' (default), then
			grep will extract lines after pattern. If flag is
			'-v', then grep will reverse match all instances 
			of pattern, and the lines parameter is arbitrary.
	"""	
	if flag == '-A':
		os.system('grep -A ' + str(int(lines)) + ' \"' + str(pattern) + '\" ' + 
					str(orig) +	' > ' + str(new))
	elif flag == '-v':
		os.system('grep -v ' + '\"' + str(pattern) + '\" ' + str(orig) +	
				  ' > ' + str(new))
#-- END Grep --#

