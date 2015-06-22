#encoding:utf8
import json
from aces.tools import shell_exec,write
import os,sys
from aces.tools import exit
import aces.config as config
import time
from aces.env import *
def genPbs(path,disp,queue,nodes,procs):
	content="cd $PBS_O_WORKDIR\n"
	content+=config.python+" %s/App.py\n"%SRCHOME	
	s="""#!/bin/bash -x
#PBS -l nodes=%s:ppn=%s
#PBS -l walltime=240:00:00
#PBS -j oe
#PBS -q %s
#PBS -N %s
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/mkl/10.0.013/lib/em64t:/opt/intel/mpi/openmpi/1.6.3/icc.ifort/lib
%s
"""%(nodes,procs,queue,disp,content)

	write(s,path+"/lammps.pbs")
def genSh(path,disp,procs):
	SRCHOME=os.path.dirname(__file__);
	qloop=' "%s/qloop.php" "%s/species.php" >input '%(path,path)
	s=["#!/bin/bash -x"
	,"#%s"%disp
	,"cd %s/minimize"%path
	,config.php+SRCHOME+"/minimize/input.php"+qloop
	,"mpirun   -np %s %s <input > %s/minimize/log.out 2>/dev/null"%(procs,config.lammps,path)
	,"cd %s"%path
	,"%s %s/input.php"%(config.php,SRCHOME)+qloop
	,"mpirun -np %s %s <input >%s/log.out 2>/dev/null"%(procs,config.lammps,path)
	,"exit 0"]
	pbs=open("%s/run.sh"%path,"w");
	pbs.write('\n'.join(s))

def genPbss(path,disp,queue,nodes,procs,start,ucores):
	nodes=int(nodes);procs=int(procs);ucores=int(ucores);
	len=int(nodes*procs/float(ucores));
	real_cores=len*ucores;
	end=start+len-1;
	if(not os.path.exist(path+"/pbs/minimize")):shell_exe("mkdir -p "+path+"/pbs/minimize");
	pbsfile=path+"/pbs/lammps%s-%s.pbs"%(start,end);
	pbs=open(pbsfile,"w");
	if(not os.path.exist(pbsfile)):
		exit("create file  pbsfile failed!\n");
	mfile=path+"/pbs/minimize/in.%s-%s.pbs"%(start,end);
	input_m=open(mfile,"w");
	ifile=path+"/pbs/in.%s-%s.pbs"%(start,end);
	input=open(ifile,"w");
	b=len;
	for i in range(start,end+1):
		p=i-start+1;
		print >>input_m,"partition yes p log %s/%s/minimize/log.lammps"%(path,i)
		print >>input_m,"partition yes p jump %s/%s/minimize/input"%(path,i)
		print >>input,"partition yes p  log %s/%s/log.lammps"%(path,i)
		print >>input,"partition yes p 	jump %s/%s/input"%(path,i)



	print >>pbs,"#!/bin/bash -x"
	print >>pbs,"#PBS -l nodes=%s:ppn=%s"%(nodes,procs)
	print >>pbs,"#PBS -l walltime=240:00:00"
	print >>pbs,"#PBS -j oe"
	print >>pbs,"#PBS -q %s"%queue
	print >>pbs,"#PBS -N %s%s-%s"%(disp,start,end)
	print >>pbs,"# Setup the OpenMPI topology:"
	print >>pbs,"n_proc=$(cat $PBS_NODEFILE | wc -l)"
	print >>pbs,"contexts=`~/bin/get_psm_sharedcontexts_max.sh`"
	print >>pbs," if [ '?' = '0' ] ; then"
	print >>pbs,"  export PSM_SHAREDCONTEXTS_MAX=contexts"
	print >>pbs," fi"

	for i in range(start,end+1):
		print >>pbs,"cd  %s/%s/minimize"%(path,i)
		print >>pbs,"%s %s/minimize/input.php \"%s/%s/qloop.php\" \"%s/%s/species.php\">input  "%(config.php,SRCHOME,path,i,path,i)
		print >>pbs,config.mpirun+real_cores+config.lammps+" -plog "+path+"/pbs/minimize/log.%s-%s "%(start,end)+"-partition %sx%s "(len,ucores)+"-pscreen "+path+"/pbs/minimize/screen.%s-%s -i %s"%(start,end,mfile)
	for i in range(start,end+1):
		print >>pbs,"cd %s/%s"%(path,i)
		print >>pbs,"%s %s/input.php \"%s/%s/qloop.php\" \"%s/%s/species.php\">input  "%(config.php,SRCHOME,path,i,path,i)
		print >>pbs,config.mpirun+real_cores+config.lammps+" -plog "+path+"/pbs/log.%s-%s "%(start,end)+"-partition %sx%s "(len,ucores)+"-pscreen "+path+"/pbs/screen.%s-%s -i %s"%(start,end,mfile)
	print >>pbs,"exit 0"


	
def makeLoopFile(idx,app):
	dir="%s/%s"%(PROJHOME,idx)
	pro="zy_%s_%s"%(PROJNAME,idx)
	universe,single,procs,nodes,unodes,uprocs,uqueue,queue=app['universe'],app['single'],app['procs'],app['nodes'],app['unodes'],app['uprocs'],app['uqueue'],app['queue']
	if(universe):
		print("prepared:"+dir );
		cores=procs*nodes;
		ucores=unodes*uprocs;
		lenn=int(ucores/(cores+0.0));
		if(idx%lenn==0):
			genPbss(PROJHOME,"zy_%s_"%PROJNAME,uqueue,unodes,uprocs,idx,cores);


	shell_exec("mkdir -p "+dir);
	if(single):
		genSh(dir,pro,procs);
	
	if(not universe and not single):genPbs(dir,pro,queue,nodes,procs);
	cores=procs*nodes
	app['cores']=cores
	write(json.dumps(app),dir+"/app.json");


def setSubProject(index,single):
	if(single):pid=''#exec.background("sh %s/%s/run.sh"%(PROJHOME,index));
	else:
		pid=shell_exec("cd %s/%s;qsub lammps.pbs;"%(PROJHOME,index))[0:5];
		print "submit: %s\t%s/%s"%(pid,PROJHOME,index);
	#sleep(1);
	return pid;
def toolsub(idx,app):
	makeLoopFile(idx,app)
	if(not app['universe']):pid=setSubProject(idx,app['single']);
	
	json_obj={
		"id":idx,
		"pid":pid,
		"time":time.strftime('%Y-%m-%d %H:%M:%S'),
	}
	json_obj=dict(json_obj,**app)


	write(json.dumps(json_obj),'qloops.txt','a+','\n')
