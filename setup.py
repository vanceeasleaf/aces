from setuptools import setup, find_packages
setup(
      name="mytest",
      version="0.10",
      description="My test module",
      author="Robin Hood",
      url="http://www.csdn.net",
      license="LGPL",
      packages= find_packages(),
      scripts=["script/ae"],
      )