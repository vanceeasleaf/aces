#encoding:utf8
import json
from query import shell_exec
import os,sys
from input import exit
import config
def genPbs(path,disp,queue,nodes,procs):
	#define MPI PATH
	home=os.path.dirname(__file__);	
	qloop= " \""+path+"/qloop.php\""+" \""+path+"/species.php\" >input"
	s=["#!/bin/bash -x"
	,"#PBS -l nodes=%s:ppn=%s"%(nodes,procs)
	,"#PBS -l walltime=240:00:00"
	,"#PBS -j oe"
	,"#PBS -q %s"%queue
	,"#PBS -N %s"%disp
	,"# Setup the OpenMPI topology"
	,"n_proc=$(cat $PBS_NODEFILE | wc -l)"
	,"contexts=`~/bin/get_psm_sharedcontexts_max.sh`"
	," if [ '?' = '0' ] ; then"
	,"  export PSM_SHAREDCONTEXTS_MAX=contexts"
	," fi"
	,"cd %s/minimize"%path
	,config.php+"%s/minimize/input.php"%home+qloop
	,config.OMPI_HOME+"/bin/mpirun  -machinefile $PBS_NODEFILE -np $n_proc "+config.APP_PATH+" <input &>"+path+"/minimize/log.out"
	,"cd "+path
	,config.php+'%s/input.php'%home+qloop
	,config.OMPI_HOME+"/bin/mpirun -machinefile $PBS_NODEFILE -np $n_proc "+config.APP_PATH+" <input &>"+path+"/log.out"
	,"exit 0"]
	pbs=open(path+"/lammps.pbs","w");
	pbs.write('\n'.join(s))
	
def genSh(path,disp,procs):
	home=os.path.dirname(__file__);
	qloop= " \""+path+"/qloop.php\""+" \""+path+"/species.php\" >input"
	s=["#!/bin/bash -x"
	,"#%s"%disp
	,"cd %s/minimize"%path
	,config.php+home+"/minimize/input.php"+qloop
	,"mpirun   -np %s %s <input > %s/minimize/log.out 2>/dev/null"%(procs,config.APP_PATH,path)
	,"cd %s"%path
	,"%s %s/input.php"%(php,home)+qloop
	,"mpirun -np %s %s <input >%s/log.out 2>/dev/null"%(procs,config.APP_PATH,path)
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
		print >>pbs,config.OMPI_HOME+"/bin/mpirun  -machinefile $PBS_NODEFILE -np "+real_cores+config.APP_PATH+" -plog "+path+"/pbs/minimize/log.%s-%s "%(start,end)+"-partition %sx%s "(len,ucores)+"-pscreen "+path+"/pbs/minimize/screen.%s-%s -i %s"%(start,end,mfile)
	for i in range(start,end+1):
		print >>pbs,"cd %s/%s"%(path,i)
		print >>pbs,"%s %s/input.php \"%s/%s/qloop.php\" \"%s/%s/species.php\">input  "%(config.php,home,path,i,path,i)
		print >>pbs,config.OMPI_HOME+"/bin/mpirun  -machinefile $PBS_NODEFILE -np "+real_cores+config.APP_PATH+" -plog "+path+"/pbs/log.%s-%s "%(start,end)+"-partition %sx%s "(len,ucores)+"-pscreen "+path+"/pbs/screen.%s-%s -i %s"%(start,end,mfile)
	print >>pbs,"exit 0"


	
def makeLoopFile(cmd,idx,projHome,projName,species,units,method,queue ,nodes ,procs,universe,uqueue,single,unodes,uprocs):
	dir="%s/%s"%(projHome,idx)
	pro="zy_%s_%s"%(projName,idx)
	if(universe):
		print("prepared:"+dir );
		cores=procs*nodes;
		if(unodes==''):unodes=20;
		if(uprocs==''):uprocs=1;
		ucores=unodes*uprocs;
		len=int(ucores/(cores+0.0));
		if(uqueue==''):uqueue="q3.4";
		if(idx%len==0):
			genPbss("projHome","zy_%s_"%projName,uqueue,unodes,uprocs,idx,cores);


	shell_exec("mkdir -p "+dir+";cd "+dir+";mkdir -p minimize");
	if(single):
		genSh(dir,pro,procs);
	
	if(not universe or not single):genPbs(dir,pro,queue,nodes,procs);
	write("<?php\n%s;\n$projHome=\"%s/%s\";\n?>"%(cmd,projHome,idx),dir+"/qloop.php");
	write("<?php\n$species=\"%s\";\n$units=\"%s\";\n$method=\"%s\";\n?>"%(species,units,method),dir+"/species.php");
def write(cmd,fileName):
	file=open(fileName,"w");
	file.write(cmd+"\n");
	file.close();

def setSubProject(index,projHome,single):
	if(single):pid=''#exec.background("sh %s/%s/run.sh"%(projHome,index));
	else:
		pid=shell_exec("cd %s/%s;qsub lammps.pbs;"%(projHome,index));
		print "submit: %s\t%s/%s"%(pid,projHome,index);
	#sleep(1);
	return pid;
def toolsub(cmd,idx,projHome,projName,species,units,method,queue ,nodes ,procs,universe,uqueue,single,unodes,uprocs):
	makeLoopFile(cmd,idx,projHome,projName,species,units,method,queue ,nodes ,procs,universe,uqueue,single,unodes,uprocs)
	if(universe==''):pid=setSubProject(idx,projHome,single);
	import time
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
		json_obj.extent(eobj[1:])
	loops=open('qloops.txt','a')
	print >>loops,json.dumps(json_obj)
if __name__=='__main__':

	projHome ,projName ,cmd ,idx ,nodes ,procs ,species ,method ,units ,universe ,queue ,runTime,uqueue,single,unodes,uprocs,jj=sys.argv[1:]
	toolsub(cmd,idx,projHome,projName,species,units,method,queue ,nodes ,procs,universe,uqueue,single,unodes,uprocs)