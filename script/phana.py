from aces import config
from aces.tools import *
def compile():
	clapack='CLAPACK-3.2.1'
	spglib='spglib-0.7.1'
	tricubic='tricubic-1.0'
	cd(config.phanalib)
	p=passthru
	if not exists(clapack):
		p('tar xvf clapack.tgz')
	if not exists(spglib):
		p('tar xvf spglib-0.7.1.tar.gz')
	if not exists(tricubic):
		p('tar xvf tricubic-1.0.tgz')
	if not exists('include'):mkdir('include')
	if not exists('lib'):mkdir('lib')
	cd(clapack)
	p("cp make.inc.example make.inc")
	p("make -j")
	cp('lapack_LINUX.a','../lib/liblapack.a')
	cp('blas_LINUX.a','../lib/libblas.a')
	cp('F2CLIBS/libf2c.a','../lib/')
	p("cp INCLUDE/*.h ../include")
	cd('..')

	cd(spglib)
	p("./configure --prefix=`pwd`")
	p("make install")
	p("cp lib/* ../lib")
	p("cp src/*.h ../include")
	cd('..')
	cd(tricubic)

	p("./configure --prefix=`pwd`")
	p("make install")
	p("cp lib/* ../lib")
	p("cp include/*.h ../include")	
	cd('..')
	phana()
def phana():
	dir=dirname(config.phana)
	p=passthru
	cd(dir)
	if not exists('oldMakefile'):
		cp("Makefile",'oldMakefile')
	s=read('oldMakefile')
	s=s.replace('-Wno-unused-result','#-Wno-unused-result')
	lib=config.phanalib+"/lib"
	include=config.phanalib+'/include'
	s=s.replace('$(SPGINC)','$(SPGINC) -I%s'%include)
	s=s.replace('$(SPGLIB)','$(SPGLIB) -L%s'%lib)
	s=s.replace('-lclapack','-llapack')
	write(s,'Makefile')
	p('make clean')
	p('make -j')
