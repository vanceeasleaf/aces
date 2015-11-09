## ## ## lj.py v.0.1
## ## ## This program returns lj parameters
## ## ## Created: 08/07/2012 - KDP

import numpy as np

def value(string):
	"""
	Returns common parameters for lj argon.

	ntpy.param.lj.value(string)
	Parameters
	----------
		string : str
		A string that corresponds to a lennard jones argon parameter.
	"""
	
	ljparams = dict({
					'lat0' : 5.269,		# Lat constant at 0 K in ang
					'lat20' : 5.315,	# Lat constant at 20 K in ang
					'lat35' : 5.355,	# Lat constant at 35 K in ang
					'lat50' : 5.401,	# Lat constant at 50 K in ang
					'lat65' : 5.455,	# Lat constant at 65 K in ang
					'lat80' : 5.527,	# Lat constant at 80 K in ang
					'epsilon' : 1.67e-21,	# Epsilon constant in joules
					'sigma' : 3.40e-10,	# Sigma constant in meters
					'mass' : 6.63e-26,	# Mass constant in kilograms
					'tau' : 2.14e-12	# Tau constant in seconds
					})	

	try:
		return ljparams[string]
	except KeyError, e:
		print "KeyError: %s is not a valid key for ntpy.param.lj()." % e
		raise	
##### END DICT

def latfit(string):
	"""
	Returns lattice constants for lj argon based off of a polyfit

ntpy.param.lj.latfit(string)
	Parameters
	----------
		string : str
			A string that corresponds to a temperature between 0 and 80K
	"""
	usertemp = float(string)

	assert usertemp <= 80 and usertemp >= 0, \
		"Error: %s is not between 0 and 80 exclusive" % usertemp

	temp = [0, 20, 35, 50, 65, 80]
	constants = [5.269, 5.315, 5.355, 5.401, 5.455, 5.527]

	fitparams = np.polyfit(temp, constants, 4)
	fit = np.poly1d(fitparams)
	return fit(usertemp)
##### END LATFIT

def const(string):
	"""
	Returns common constant parameters.

	ntpy.param.const(string)
	Parameters
	----------
		string : str
			A string that corresponds to a constant parameter.
	"""

	constparams = dict({
					'kb' : 1.3806e-23,	# Boltzmann's constant
					'hbar' : 1.054e-34,	# Planck's constant
					'topeta' : 1e15,	# To peta-
					'totera' : 1e12,	# To tera-
					'togiga' : 1e9,	# To giga-
					'tomega' : 1e6,	# To mega-
					'tokilo' : 1e3,	# To kilo-
					'tocenti' : 1e-2,	# To centi-
					'tomilli' : 1e-3,	# To milli-
					'tomicro' : 1e-6,	# To micro-
					'tonano' : 1e-9,	# To nano-
					'topico' : 1e-12,	# To pico-
					'tofemto' : 1e-15,	# To femto-
					})

	try:
		return constparams[string]
	except KeyError, e:
		print "KeyError: %s is not a valid key for ntpy.param.const()." % e
		raise
	
##### END LJ

