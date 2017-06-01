if __name__=='__main__':
	from aces.tools import passthru
	from aces import config
	passthru(config.python+"setup.py build --build-lib=. --build-platlib=.")