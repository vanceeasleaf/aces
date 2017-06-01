from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy
include_dirs=[numpy.get_include()]
extensions=[]
extensions.append(Extension("gett", ["gett.pyx"],include_dirs=include_dirs))
setup(
  name = 'gett',
  ext_modules=cythonize(extensions)
)