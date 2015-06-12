#encoding:utf8
import json
import os,sys
import aces.config as config
from aces.inequality import inequality
from ase.io import read
from aces.tools import shell_exec,mkdir,cd,passthru,pwrite
from aces.profile import proc
def getObjs():
	
	#执行程序之前会清理上次的，读取qloops.txt，如果没有上次的文件就不清理。*/
	fileName="qloops.txt"
	qloop=open(fileName,"r");
	obj=[]
	for  json_string in qloop:
        	obj.append(json.loads(json_string));
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
# 获取有效参数，排除那些已经考虑过的
# @author zhouy
##
def getParas(obj):
	paras=[]
	for pa in obj:
		for key in pa:
			if key in ["id","pid","time","cmd","project","nodes","procs","species","units","method"]:continue;
			if not(key in paras):paras.append(key)
	return paras;
     

	
def getQueryInfo(workPath,pid,runTime,ob):
	lastline=shell_exec("tail -1 %s/log.out"%workPath);
	qstat=shell_exec("qstat %s 2>&1|tail -1 "%pid);
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
	
def postProc(curPath,srcHome):
	""""php=config.php
	# 准备参数列表并调用后处理*/
	dir=curPath;
	dir=dir.replace("//","\\\/");
	sed=" sed 's@projHome=.\+@projHome=\""+dir+"\";@g' qloop.php > qloop.php1";
	if(os.path.exists("post.php")):postfile= "../post.php";
	else: postfile="";
	cmd="cd %s;"%curPath+sed+";cat qloop.php1 "+postfile+" > qloop.php2;"+php+" %s/profile.php \"%s/qloop.php2\" \"%s/species.php\";  "%(srcHome,curPath,curPath)
	#print cmd
	shell_exec(cmd)"""
	cd(curPath)
	proc()
	
def kappa(curPath,result):
	# 取出后处理结果，热导率*/
	#print curPath
	kappaline=shell_exec("cd %s;tail -1 result.txt 2>err;"%curPath);
	kappa=kappaline.split('=');
	if len(kappa)>1:
		kappa=kappa[1]
		pwrite(result,"%s"%kappa);	
		
def tEnerty(curPath,result):
	# 总能量*/
	totalEline=shell_exec("cd %s/minimize;tail -22 log.out| head -1;"%curPath);
	totalE=totalEline.split()[1]
	pwrite(result,"\t%s"%totalE);	
	
def nAtom(curPath,result):
	# 原子数和平均能量*/
	Natomline=shell_exec("cd %s/minimize;grep atoms log.out ;"%curPath);
	Natom=Natomline.split()[0]
	if(Natom.isdigit() and Natom>0):
		pwrite(result,"\t%s"%Natom);
		#epn=float(totalE)/float(Natom);        	          
		#pwrite(result,"\t%f"%epn);	
def tDisorder(curPath,srcHome,result):
	
	# 无序度*/
	cd('%s/minimize'%curPath)
	mkdir('disorder');cd('disorder')
	disorderLine=shell_exec("cp %s"%srcHome+"/in.disorder .;"+config.lammps+" <in.disorder 2>err 1>log;tail -1 disorder.txt  2>err;");
	k=disorderLine.split()[1:3]				
	if len(k)==1:
		k.append("")
	disorder,rd=k
	cd(curPath)
	pwrite(result,"\t%s\t%s"%(disorder,rd));
	
def drawStructure(curPath):
	cd('%s/minimize'%curPath)
	atoms=read('range',format='lammps')
	atoms.write('../structure.png')		
	cd(curPath)
	
def ineq(ob,curPath,result):

	species=ob["species"];
	if(not (species in ["CN-small"])):return
	cd('%s/minimize'%curPath)
	mkdir('nonequ')
	cd('nonequ')

	ie=inequality()
	nonequ5= ie.run()
	#nonequ5=shell_exec("cd nonequ;python %s/inequality.py;"%srcHome);
	cd(curPath)
	pwrite(result,"\t%s"%nonequ5);	
	
def item(projHome,ob,result,paras,srcHome):
	php=config.php
	id=ob["id"];
	curPath="%s/%s"%(projHome,id)
	pid=ob["pid"];
	runTime=ob["runTime"];
	if(runTime==""):runTime=10000000;

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
		postProc(curPath,srcHome)
		kappa(curPath,result)
		tEnerty(curPath,result)
		nAtom(curPath,result)
		#tDisorder(curPath,srcHome,result)
		drawStructure(curPath)
		ineq(ob,curPath,result)
		pwrite(result,"\n");	
		
def query(projHome,srcHome,universe):


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
	for para in paras:
		pwrite(result,"\t%s"%para);

	# work的计算结果*/
	pwrite(result,"\tkappa\ttotalE\tNatom\tE/N\tdisorder\trd\tdisorderC\tratio\trdfs\n");

	checkUniverse(projHome,universe,obj)

	# 遍历projet中的所有work*/
	for ob in obj:
		item(projHome,ob,result,paras,srcHome)
			
def checkUniverse(projHome,universe,obj):
	if(universe==''):return;
	ff=open("%s/pbs/info"%projHome);
	uid=[];upid=[];ulog=[];uscreen=[];
	for line in ff:
		a=line.split()
		uid.append(int(a[0]));upid.append(a[1]);
		ulog.append(a[2]);uscreen.append(a[3]);
	i=0;
	for ob in obj:
		id=ob[id];ob[pid]=int(upid[i])
		if(not os.path.exists("%s/pbs/%s"%(projHome,ulog[i]))):continue;
		if(not os.path.exists("%s/pbs/%s"%(projHome,uscreen[i]))):continue;
		shell_exec("cp %s/pbs/%s %s/%s/log.lammps"%(projHome,ulog[i],projHome,id));
		shell_exec("cp %s/pbs/%s %s/%s/log.out"%(projHome,uscreen[i],projHome,id));
		shell_exec("cp %s/pbs/minimize/%s %s/%s/minimize/log.lammps"%(projHome,ulog[i],projHome,id));
		shell_exec("cp %s/pbs/minimize/%s %s/%s/minimize/log.out"%(projHome,uscreen[i],projHome,id));
		i+=1


def clean(projHome,projName,single):
	if projHome=='':print "error projHome"
	#/*删除原始代码以外的文件*/	
	files=shell_exec("cd %s;ls "%projHome);
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
	comfirmStop(projHome,projName,single);
	for ls in deleteFiles:
		shell_exec("cd %s;rm -r %s"%(projHome,ls));
def stop(projHome,projName,single):
	print "Comfirm to stop all the simulation in this project?[y/n]",	
	sys.stdout.flush()
	s=raw_input(); 
	if(s!="y"):exit("exit with no change.")
	comfirmStop(projHome,projName,single);
			

def comfirmStop(projHome,projName,single):
	
	#/* 容易kill掉同名工程程序*/
	if(single):
		obj=getObjs(projHome+"/qloops.txt");
		for pa in obj:
			pid=pa["pid"];
         		print "kill:%s"%pid;
         	 	#exec::kill(pid);
		return;
	
	tarname="zy_%s_"%projName;
	if(len(tarname)<10):
		works=shell_exec("qstat|grep %s"%tarname).split('\n');
		for work in works:
			if work.strip()=="":continue;
			pid=work.split()[0]
			print "qdel:%s"%pid
			sys.stdout.flush()
			shell_exec("qdel %s &2>/dev/null"%pid)
	else:
		works=shell_exec("qstat|grep xggong").split('\n');
		for work in works:
			if work.strip()=="":continue;
			pid=work.split()[0]
			jobnameString=shell_exec("qstat -f %s |grep Job_Name"%pid);
			jobname=jobnameString.strip().split()[2]
			if(jobname.find(tarname)>=0):
				print "qdel:%s"%pid
				sys.stdout.flush()
				shell_exec("qdel %s "%pid)
				
	
		
	
		
if __name__=='__main__':
	cmd,projHome,srcHome,universe,projName,single=sys.argv[1:]
	if cmd=='q':
		query(projHome,srcHome,universe)
	elif cmd=='clean':
		clean(projHome,projName,single)
	elif cmd=='stop':
		stop(projHome,projName,single)