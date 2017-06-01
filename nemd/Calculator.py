	
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
def tDisorder(curPath,result):
	
	# 无序度*/
	cd('%s/minimize'%curPath)
	mkdir('disorder');cd('disorder')
	disorderLine=shell_exec("cp %s"%SRCHOME+"/in.disorder .;"+config.lammps+" <in.disorder 2>err 1>log;tail -1 disorder.txt  2>err;");
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
	#nonequ5=shell_exec("cd nonequ;python %s/inequality.py;"%SRCHOME);
	cd(curPath)
	pwrite(result,"\t%s"%nonequ5);	
		