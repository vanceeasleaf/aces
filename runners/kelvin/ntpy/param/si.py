## ## ## si.py v.0.1
## ## ## This program returns lj parameters
## ## ## Created: 04/22/2013 - KDP

import numpy as np

def value(string):
	"""
	Returns common parameters for silicon.

	ntpy.param.si.value(string)
	Parameters
	----------
		string : str
		A string that corresponds to a lennard jones argon parameter.
	"""
	
	siparams = dict({
					'lat' : 5.43,	# Lat constant in ang
					'latamor' : 5.43*1.0438, # amorphous lat constant in ang
					'sound' : 8433, # Speed of sound
					'mass' : 28.0855	# Mass constant in grams
					})	

	try:
		return siparams[string]
	except KeyError, e:
		print "KeyError: %s is not a valid key for ntpy.param.si()." % e
		raise	
##### END VALUE

