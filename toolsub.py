#encoding:utf8
import json
from aces.tools import shell_exec,write
import os,sys
from aces.input import exit
import aces.config as config
import time
def genPbs(path,disp,queue,nodes,procs,bte):
	#define MPI PATH
	home=os.path.dirname(__file__);	
	global usepy
	usepy=1

	#if not 'usepy' in vars():usepy=0
	if usepy:
		input='exe.py'
		exe=config.python
	else:
		input='input.php  "%s/qloop.php" "%s/species.php" '%(path,path)
		exe=config.php
	qloop=input+' >input'
	if bte==True:		
		content="cd %s/minimize\n"%path
		content+="%s %s/minimize/%s\n"%(exe,home,qloop)
		content+=config.mpirun+" $n_proc "+config.lammps+" <input &>"+path+"/minimize/log.out"
		content+="\ncd "+path+"\n"
		content+=exe+"%s/%s\n"%(home,qloop)
		content+=config.mpirun+"$n_proc "+config.phonts+"  &>"+path+"/log.out"
	elif bte=="correlation":
		content="cd %s/minimize\n"%path
		content+=exe+"%s/minimize/%s\n"%(home,qloop)
		content+=config.mpirun+"$n_proc "+config.lammps+" <input &>"+path+"/minimize/log.out"
		content+="\ncd "+path+"\n"
		content+=exe+"%s/%s\n"%(home,qloop)
	else:
		content="cd %s/minimize\n"%path
		content+=exe+"%s/minimize/%s\n"%(home,qloop)
		content+=config.mpirun+"$n_proc "+config.lammps+" <input &>"+path+"/minimize/log.out"
		content+="\ncd "+path+"\n"
		content+=exe+"%s/%s\n"%(home,qloop)
		content+=config.mpirun+"$n_proc "+config.lammps+" <input &>"+path+"/log.out"
		
	s="""#!/bin/bash -x
#PBS -l nodes=%s:ppn=%s
#PBS -l walltime=240:00:00
#PBS -j oe
#PBS -q %s
#PBS -N %s
# Setup the OpenMPI topology
n_proc=$(cat $PBS_NODEFILE | wc -l)
contexts=`~/bin/get_psm_sharedcontexts_max.sh`
if [ '?' = '0' ] ; then
  export PSM_SHAREDCONTEXTS_MAX=contexts
 fi
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/mkl/10.0.013/lib/em64t:/opt/intel/mpi/openmpi/1.6.3/icc.ifort/lib
%s
exit 0
"""%(nodes,procs,queue,disp,content)
	pbs=open(path+"/lammps.pbs","w");
	pbs.write(s)
	
def genSh(path,disp,procs):
	home=os.path.dirname(__file__);
	qloop=' "%s/qloop.php" "%s/species.php" >input '%(path,path)
	s=["#!/bin/bash -x"
	,"#%s"%disp
	,"cd %s/minimize"%path
	,config.php+home+"/minimize/input.php"+qloop
	,"mpirun   -np %s %s <input > %s/minimize/log.out 2>/dev/null"%(procs,config.lammps,path)
	,"cd %s"%path
	,"%s %s/input.php"%(config.php,home)+qloop
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


	home=os.path.dirname(__file__);
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
		print >>pbs,"%s %s/minimize/input.php \"%s/%s/qloop.php\" \"%s/%s/species.php\">input  "%(config.php,home,path,i,path,i)
		print >>pbs,config.mpirun+real_cores+config.lammps+" -plog "+path+"/pbs/minimize/log.%s-%s "%(start,end)+"-partition %sx%s "(len,ucores)+"-pscreen "+path+"/pbs/minimize/screen.%s-%s -i %s"%(start,end,mfile)
	for i in range(start,end+1):
		print >>pbs,"cd %s/%s"%(path,i)
		print >>pbs,"%s %s/input.php \"%s/%s/qloop.php\" \"%s/%s/species.php\">input  "%(config.php,home,path,i,path,i)
		print >>pbs,config.mpirun+real_cores+config.lammps+" -plog "+path+"/pbs/log.%s-%s "%(start,end)+"-partition %sx%s "(len,ucores)+"-pscreen "+path+"/pbs/screen.%s-%s -i %s"%(start,end,mfile)
	print >>pbs,"exit 0"


	
def makeLoopFile(cmd,idx,projHome,projName,species,units,method,queue ,nodes ,procs,universe,uqueue,single,unodes,uprocs,jj,bte):
	dir="%s/%s"%(projHome,idx)
	pro="zy_%s_%s"%(projName,idx)
	if(universe):
		print("prepared:"+dir );
		cores=procs*nodes;
		if(unodes==''):unodes=20;
		if(uprocs==''):uprocs=1;
		ucores=unodes*uprocs;
		lenn=int(ucores/(cores+0.0));
		if(uqueue==''):uqueue="q3.4";
		if(idx%lenn==0):
			genPbss(projHome,"zy_%s_"%projName,uqueue,unodes,uprocs,idx,cores);


	shell_exec("mkdir -p "+dir+";cd "+dir+";mkdir -p minimize");
	if(single):
		genSh(dir,pro,procs);
	
	if(not universe and not single):genPbs(dir,pro,queue,nodes,procs,bte);
	write('<?php\n%s;\n$projHome="%s/%s";\n?>'%(cmd,projHome,idx),dir+"/qloop.php");
	write('<?php\n$species="%s";\n$units="%s";\n$method="%s";\n?>'%(species,units,method),dir+"/species.php");
	eobj=json.loads(jj)
	if len(eobj)>1:
		opt=eobj[1].copy()
		opt['projHome']='%s/%s'%(projHome,idx)
		opt['species']=species
		opt['units']=units
		opt['method']=method
		opt['bte']=bte
		write(json.dumps(opt),dir+"/app.json");


def setSubProject(index,projHome,single):
	if(single):pid=''#exec.background("sh %s/%s/run.sh"%(projHome,index));
	else:
		pid=shell_exec("cd %s/%s;qsub lammps.pbs;"%(projHome,index))[0:5];
		print "submit: %s\t%s/%s"%(pid,projHome,index);
	#sleep(1);
	return pid;
def toolsub(cmd,idx,projHome,projName,species,units,method,queue ,nodes ,procs ,runTime,jj,universe,uqueue,single,unodes,uprocs,bte):
	makeLoopFile(cmd,idx,projHome,projName,species,units,method,queue ,nodes ,procs,universe,uqueue,single,unodes,uprocs,jj,bte)
	if(universe==''):pid=setSubProject(idx,projHome,single);

	json_obj={
		"id":idx,
		"pid":pid,
		"time":time.strftime('%Y-%m-%d %H:%M:%S'),
		"cmd":cmd,
		"nodes":nodes,
		"procs":procs,
		"species":species,
		"method":method,
		"units":units,
		"runTime":runTime
	}
	eobj=json.loads(jj)
	if len(eobj)>1:
		json_obj=dict(json_obj.items()+eobj[1].items())
	loops=open('qloops.txt','a')
	print >>loops,json.dumps(json_obj)
if __name__=='__main__':

	projHome ,projName ,cmd ,idx ,nodes ,procs ,species ,method ,units ,universe ,queue ,runTime,uqueue,single,unodes,uprocs,jj=sys.argv[1:]
	toolsub(cmd,idx,projHome,projName,species,units,method,queue ,nodes ,procs ,runTime,jj,universe,uqueue,single,unodes,uprocs,bte)