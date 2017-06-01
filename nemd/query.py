#encoding:utf8
import json
import os,sys
import aces.config as config
from aces.inequality import inequality

from aces.tools import *
from ase.io import read
from aces.env import *
from aces.App import App
def getObjs():
	obj=[json.loads(json_string) for  json_string in open("qloops.txt")]
	return obj;
	
def getRatio(path):
	if(not os.path.exists(path)):return 0.0;
	fp=open(path,"r");
	fp.next();
	natom=int(fp.next().split()[0])
	ntype=int(fp.next().split()[0])
	if ntype==1:return 0.0;
	n=0;
	label=""
	while(label!="Atoms" and n<20):
		label=fp.next().strip();
		n+=1
	
	fp.next()
	a=[0.0]*ntype
	for line in fp:
		type=int(line.split()[1])
		a[type-1]+=1
	return float(a[1])/natom;
	
###
# effective parameters 
# @author zhouy
##
def getParas(obj):
	x=reduce(lambda x,y:set(x)|set(y),obj)
	a=["id","pid","time","cmd","project","nodes","procs","species","units","method"]
	return list(set(x).difference(set(a))) 
     

	
def getQueryInfo(workPath,pid,runTime,ob):
	pid=filter(str.isdigit, str(pid))
	lastline=shell_exec("tail -3 %s/log.out"%workPath);
	qstat=shell_exec("qstat %s 2>&1|tail -1 "%pid )

	step=lastline.split()[0]
	if step.isdigit():
		percent="%.1f%%"%(float(step)/runTime*100)
	else:percent="0"
	if(qstat.find("Unknown Job Id")>=0):#该任务已从任务队列中去除*/
		time="complete";
		if(lastline.find("builds")>=0):
			status="C";
			percent="100%";
		else:#异常退出*/
			status="E";
		
		queue=""#ob["queue"];
		nodes=""#ob["nodes"];
		procs=""#ob["procs"];
	else:#正在运行或等待R&Q&C*/
		time,status,queue=qstat.split()[3:6]
		info=shell_exec("qstat -f %s 2>&1|grep nodes"%pid)
		info=info.split()[2]
		nnn=info.split(":ppn=");
		nodes=nnn[0];
		procs=nnn[1];
	
	return (percent,status,queue,nodes,procs);
	

def kappa(result):
	kappaline=shell_exec("tail -1 result.txt 2>err;");
	kappa=kappaline.split('=');
	if len(kappa)>1:
		kappa=kappa[1]
		pwrite(result,"%s"%kappa);	
		
def tEnerty(result):
	# 总能量*/
	totalEline=shell_exec("cd minimize;tail -22 log.out| head -1;");
	totalE=totalEline.split()[1]
	pwrite(result,"\t%s"%totalE);	
	
def nAtom(result):
	# 原子数和平均能量*/
	Natomline=shell_exec("cd minimize;grep atoms log.out ;");
	Natom=Natomline.split()[0]
	if(Natom.isdigit() and Natom>0):
		pwrite(result,"\t%s"%Natom);
		#epn=float(totalE)/float(Natom);        	          
		#pwrite(result,"\t%f"%epn);	
def tDisorder(result):
	
	# 无序度*/
	now=pwd()
	cd('minimize')
	mkdir('disorder');cd('disorder')
	disorderLine=shell_exec("cp %s"%SRCHOME+"/in.disorder .;"+config.lammps+" <in.disorder 2>err 1>log;tail -1 disorder.txt  2>err;");
	k=disorderLine.split()[1:3]				
	if len(k)==1:
		k.append("")
	disorder,rd=k
	cd(now)
	pwrite(result,"\t%s\t%s"%(disorder,rd));
	
def drawStructure():
	atoms=read('minimize/range',format='lammps')
	if len(atoms)<200:
		atoms.write('minimize.png')		
	
def ineq(ob,result):
	now=pwd()
	species=ob["species"];
	if(not (species in ["CN-small"])):return
	cd('minimize')
	mkdir('nonequ')
	cd('nonequ')

	ie=inequality()
	nonequ5= ie.run()
	cd(now)
	pwrite(result,"\t%s"%nonequ5);	
	
def item(ob,result,paras):
	id=ob["id"];
	curPath="%s/%s"%(PROJHOME,id)
	pid=ob["pid"];
	runTime=ob["runTime"];

	# work的固有属性*/
	percent,status,queue,nodes,procs=getQueryInfo(curPath,pid,runTime,ob)
	print '\t'.join([str(id),percent,status,queue]),
	result.write("%d\t"%id);

	# work覆盖过的参数*/
	for key in paras:
		if ob[key]=="":ob[key]="def"
		print ob[key],
		if(float(percent.replace('%',''))>0.5):
			result.write("%s\t"%ob[key])
		 
	# work的计算结果*/
	if(float(percent.replace('%',''))>0):
		cd(curPath)
		App().result()
		kappa(result)
		#tEnerty(result)
		nAtom(result)
		#tDisorder(result)
		drawStructure()
		ineq(ob,result)
		pwrite(result,"\n");	
		
def query(universe):


	result=open("result.txt","w");

	  #*
	  # qloops.txt的处理分为以下几种情况:
	  # 正常情况：通过php sub.php投任务并生成qloops.txt以后
	  # 工程被复制以后
	  # 工程被移动以后
	  # 该文件不存在时
	  #/
	obj=getObjs();

	# work的固有属性*/
	print "id\tpercent\tstatus\tqueue\tprocs";
	result.write("id");

	# work覆盖过的参数*/
	paras=getParas(obj);
	pwrite(result,"\t%s"%('\t'.join(paras)));
	# work的计算结果*/
	pwrite(result,"\tkappa\ttotalE\tNatom\tE/N\tdisorder\trd\tdisorderC\tratio\trdfs\n");

	checkUniverse(universe,obj)

	# 遍历projet中的所有work*/
	for ob in obj:
		item(ob,result,paras)
			
def checkUniverse(universe,obj):
	if(universe==''):return;
	ff=open("%s/pbs/info"%PROJHOME);
	uid=[];upid=[];ulog=[];uscreen=[];
	for line in ff:
		a=line.split()
		uid.append(int(a[0]));upid.append(a[1]);
		ulog.append(a[2]);uscreen.append(a[3]);
	i=0;
	for ob in obj:
		id=ob[id];ob[pid]=int(upid[i])
		if(not os.path.exists("%s/pbs/%s"%(PROJHOME,ulog[i]))):continue;
		if(not os.path.exists("%s/pbs/%s"%(PROJHOME,uscreen[i]))):continue;
		shell_exec("cp %s/pbs/%s %s/%s/log.lammps"%(PROJHOME,ulog[i],PROJHOME,id));
		shell_exec("cp %s/pbs/%s %s/%s/log.out"%(PROJHOME,uscreen[i],PROJHOME,id));
		shell_exec("cp %s/pbs/minimize/%s %s/%s/minimize/log.lammps"%(PROJHOME,ulog[i],PROJHOME,id));
		shell_exec("cp %s/pbs/minimize/%s %s/%s/minimize/log.out"%(PROJHOME,uscreen[i],PROJHOME,id));
		i+=1


def clean(single):
	if PROJHOME=='':print "error PROJHOME"
	#/*删除原始代码以外的文件*/	
	files=shell_exec("cd %s;ls "%PROJHOME);
	deleteFiles=[]
	files=files.split('\n')
	for ls in files:
		if(ls in ["sub.php","post.php","data","","sub.py"] ):continue;
		#print "deleting:%s"%ls;
		deleteFiles.append(ls)
	if len(deleteFiles)>0:
		print "All the files in "+str(deleteFiles)+" will be deleted."
		print "Comfirm ?[y/n]",	
		sys.stdout.flush()
		s=raw_input();
		if(s!="y"):exit("exit with no change.")
	comfirmStop(single);
	for ls in deleteFiles:
		shell_exec("cd %s;rm -r %s"%(PROJHOME,ls));
def stop(single):
	print "Comfirm to stop all the simulation in this project?[y/n]",	
	sys.stdout.flush()
	s=raw_input(); 
	if(s!="y"):exit("exit with no change.")
	comfirmStop(single);
			

def comfirmStop(single):
	
	#/* 容易kill掉同名工程程序*/
	if(single):
		obj=getObjs(PROJHOME+"/qloops.txt");
		for pa in obj:
			pid=pa["pid"];
         		print "kill:%s"%pid;
         	 	#exec::kill(pid);
		return;
	
	tarname="zy_%s_"%PROJNAME;
	if(len(tarname)<12):
		works=shell_exec("qstat|grep %s"%tarname).split('\n');
		for work in works:
			if work.strip()=="":continue;
			if ' C ' in work:continue
			pid=work.split()[0]
			print "qdel:%s"%pid
			sys.stdout.flush()
			shell_exec("qdel %s"%pid)
	else:
		works=shell_exec("qstat|grep xggong").split('\n');
		for work in works:
			if work.strip()=="":continue;
			if ' C ' in work:continue
			pid=work.split()[0]
			jobnameString=shell_exec("qstat -f %s |grep Job_Name"%pid);
			jobname=jobnameString.strip().split()[2]
			if(jobname.find(tarname)>=0):
				print "qdel:%s"%pid
				sys.stdout.flush()
				shell_exec("qdel %s "%pid)
				
	
		
	
