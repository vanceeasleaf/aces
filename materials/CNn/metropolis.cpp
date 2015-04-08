#include<stdio.h>
#include<math.h>
#include <string.h>
#include<stdlib.h>
#include<time.h>
void createTopo(int* nnei,int neighbor[][20],int natom,double * x,double* y,double *z
,double lx,double ly,double lz,double cut);
bool within(double x,double y,double z,double cut);
char* getString(int m,char* temp);
void drawTopo(char* comment,int *nnei,int  neighbor[][20],int natom,double* x,double* y,double* z,int * type);
double ave(int n,double arr);
double  metropolis(int* nnei,int neighbor[][20],int* fnnei,int fneighbor[][20],int natom,double * x,double *y ,double *z ,int * type);
double random01();
int main(){
	srand((unsigned)time(NULL));
	FILE* fp=fopen("structure","r");
	char s[100];
	char comment[200];
	fgets(comment,200,fp);//comment
	int natom,ntype;
	double xlo,xhi,ylo,yhi,zlo,zhi;
	fscanf(fp,"%d%s",&natom,s);
	fscanf(fp,"%d%s",&ntype,s);
	fgets(s,100,fp);

	fscanf(fp,"%lf%lf%s%s",&xlo,&xhi,s,s);
	fscanf(fp,"%lf%lf%s%s",&ylo,&yhi,s,s);
	fscanf(fp,"%lf%lf%s%s",&zlo,&zhi,s,s);
	double lx=xhi-xlo;
	double ly=yhi-ylo;
	double lz=zhi-zlo;
	for(int i=0;i<9;i++)fgets(s,100,fp);
	int id;
		int type[natom];
			double x[natom],y[natom],z[natom];
	for(int i=0;i<natom;i++){
		fscanf(fp,"%d%d%lf%lf%lf",&id,&type[i],&x[i],&y[i],&z[i]);
		//printf("%d\t%d\t%f\t%f\t%f\n",id,type[i],x[i],y[i],z[i]);
	}
	fclose(fp);
	int neighbor[natom][20];
	int nnei[natom];
	double cut=2;
	createTopo(nnei,neighbor,natom,x,y,z,lx,ly,lz,cut);
	printf("the average neighbor of topo is %f.\n",ave(natom,nnei));
	//use to calculate the energy.
	double fcut=4;
	int fneighbor[natom][20];
	int fnnei[natom];
	createTopo(fnnei,fneighbor,natom,x,y,z,lx,ly,lz,fcut);
	printf("the average neighbor of force is %f.\n",ave(natom,fnnei));
//drawTopo(comment,nnei,neighbor,natom,x,y,z,type);
printf("step\tenergy\tdenergy\tsigma\n");
double energy=totalEnergy(fnnei,fneighbor,x,y,z,type);
printf("0\t%f\t0.0\t0.0\n",energy);
int step=0;
for(int i=0;i<1e6;i++){
	step++;
	double datae=metropolis(nnei,neighbor,fnnei,fneighbor,natom,x,y,z,type);
	if(!step%1000){
	printf("%d\t%f\t%f\t%f\n",step,energy,datae,fabs(datae/energy));
	}
	energy+=datae;
}
	return 0;
}
double ave(int n,double arr){
		double sum=0;
	for(int i=0;i<n;i++){
		sum+=arr[i];
	}
	return sum/n;
}
double  metropolis(int* nnei,int neighbor[][20],int* fnnei,int fneighbor[][20],int natom,double * x,double *y ,double *z ,int * type){
	int selsite=random01()*natom;
	int selnei_index=random01()*nnei[selsite];
	int selnei=neighbor[selsite][selnei_index];
	double eorigin=energy(selsite,selnei,fnnei,fneighbor,x,y,z,type);
	swap(type[selsite],type[selnei]);
	double efinal=energy(selnei,selsite,fnnei,fneighbor,x,y,z,type);
	double denergy=efinal-eorigin;
	if(denergy>0){
		double r=exp(-denergy/(300*8.617343E-5));
		if(random01()>r){
			swap(type[selsite],type[selnei]);
		return 0;
	}
	}
	return denergy;
}
double energy(int selsite,int selnei,double* fnnei,double fneighbor[][20],double* x,double* y,double* z,int *type){
	
}
double random01(){
	return rand()/(RAND_MAX+1);
}
void drawTopo(char* comment,int *nnei,int  neighbor[][20],int natom,double* x,double* y,double* z,int * type){
		FILE* fp=fopen("topo.pdb","w");
		if(fp==NULL){
			printf("fail to open topo.pdb!\n");exit(1);
		}
fprintf(fp,"COMPND    %s",comment);	
fprintf(fp,"COMPND   Created by metropolis\n");

		for(int i=0;i<natom;i++){
			char name[2]="";//至少为2,否则fp有问题，why?	
			if(type[i]==1)sprintf(name,"C");
			if(type[i]==2)sprintf(name,"N");
						
			char temp[20]="";

			sprintf(temp,"%d",i+1);
	char* a_id=getString(5,temp);
	char* a_name=getString(4,name);		
char* a_A=getString(14,"");
sprintf(temp,"%.3f",x[i]);
char* a_x=getString(8,temp);
sprintf(temp,"%.3f",y[i]);
char* a_y=getString(8,temp);
sprintf(temp,"%.3f",z[i]);
char* a_z=getString(8,temp);
char* a_occ=getString(6,"1.00");
char* a_tem=getString(6,"1.00");
char* a_ft=getString(3,"");
char* a_seg=getString(4,"");
char* a_ele=getString(2,name);
char* a_cha=getString(2,"");
char ss[100];
fprintf(fp,"HETATM%s %s%s%s%s%s%s %s%s  %s%s%s\n",a_id,a_name,a_A,a_x,a_y,a_z,a_occ,a_tem,a_ft,a_seg,a_ele,a_cha);
		}
		for(int i=0;i<natom;i++){
						char temp[20];
						sprintf(temp,"%d",i+1);
							char* a_id=getString(5,temp);
fprintf(fp,"CONECT%s",a_id);
for(int j=0;j<nnei[i];j++){
	char temp[20];
	sprintf(temp,"%d",neighbor[i][j]+1);
							char* a_id=getString(5,temp);
	fprintf(fp,"%s",a_id);
}
fprintf(fp,"\n");
		}
		fprintf(fp,"END\n");
		fclose(fp);
	}
char* getString(int m,char* temp){
	char* a_id=new char[m];
				int n=strlen(temp);
for(int k=0;k<m-n;k++){a_id[k]=' ';}
for(int k=m-n;k<m;k++){a_id[k]=temp[k-(m-n)];}
return a_id;
}
void createTopo(int* nnei,int neighbor[][20] ,int natom,double * x,double* y,double *z
,double lx,double ly,double lz,double cut){
	for(int i=0;i<natom;i++){
		nnei[i]=0;
		for(int j=0;j<natom;j++){
			if(i==j)continue;
			if(
				within(x[i]-x[j],y[i]-y[j],z[i]-z[j],cut)
				||within(x[i]-x[j]+lx,y[i]-y[j],z[i]-z[j],cut)
				||within(x[i]-x[j]-lx,y[i]-y[j],z[i]-z[j],cut)
				||within(x[i]-x[j],y[i]-y[j]+ly,z[i]-z[j],cut)
				||within(x[i]-x[j],y[i]-y[j]-ly,z[i]-z[j],cut)
				||within(x[i]-x[j],y[i]-y[j],z[i]-z[j]+lz,cut)
				||within(x[i]-x[j],y[i]-y[j],z[i]-z[j]-lz,cut)
			||within(x[i]-x[j]+lx,y[i]-y[j]+ly,z[i]-z[j],cut)
			||within(x[i]-x[j]-lx,y[i]-y[j]+ly,z[i]-z[j],cut)
			||within(x[i]-x[j]+lx,y[i]-y[j]-ly,z[i]-z[j],cut)
			||within(x[i]-x[j]-lx,y[i]-y[j]-ly,z[i]-z[j],cut)	){
				neighbor[i][nnei[i]++]=j;
			}
		}
	}
}
bool within(double x,double y,double z,double cut){
	return x*x+y*y+z*z<cut*cut;
}
