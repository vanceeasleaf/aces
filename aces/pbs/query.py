#encoding:utf8
import json
import os,sys

from aces.tools import *
from aces.env import *
def getObjs():
	obj=[json.loads(json_string) for  json_string in open("qloops.txt")]
	return obj;

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
				
	
		
	
