## ## ## const.py v.0.1
## ## ## This program returns common parameters
## ## ## Created: 08/07/2012 - KDP

def value(string):
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
					'avog' : 6.0221414e23, # Avogadro's number
					'c' : 299792485.0, # Speed of light m/s

					'J/eV' : 1.60217657e-19, # Joules per eV
					'eV/J' : 6.24150934e18, # eV per Joule

					'topeta' : 1e-15,	# To peta-
					'totera' : 1e-12,	# To tera-
					'togiga' : 1e-9,	# To giga-
					'tomega' : 1e-6,	# To mega-
					'tokilo' : 1e-3,	# To kilo-
					'tocenti' : 1e2,	# To centi-
					'tomilli' : 1e3,	# To milli-
					'tomicro' : 1e6,	# To micro-
					'tonano' : 1e9,		# To nano-
					'toang' : 1e10,		# To anstrom
					'topico' : 1e12,	# To pico-
					'tofemto' : 1e15,	# To femto-

					'peta' : 1e15,	# peta-
					'tera' : 1e12,	# tera-
					'giga' : 1e9,	# giga-
					'mega' : 1e6,	# mega-
					'kilo' : 1e3,	# kilo-
					'centi' : 1e-2,	# centi-
					'milli' : 1e-3,	# milli-
					'micro' : 1e-6,	# micro-
					'nano' : 1e-9,	# nano-
					'ang' : 1e-10,	# angstrom
					'pico' : 1e-12,	# pico-
					'femto' : 1e-15,	# femto-
					})

	try:
		return constparams[string]
	except KeyError, e:
		print "KeyError: %s is not a valid key for ntpy.param.const()." % e
		raise
	
##### END VALUE

