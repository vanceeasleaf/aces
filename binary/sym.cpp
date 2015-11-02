#include<stdio.h>
int main(){
char filename[20]="mesh.yaml";
	char s[100];
	FILE* fp=fopen(filename,"r");
	long* lineNumber=new long[1000000000];
	long nLine=0;
	lineNumber[nLine++]=0;
	fseek(fp,0,SEEK_SET);
	printf("Preparing Lines...\n");
	while(fgets(s,100,fp)){
		lineNumber[nLine++]=ftell(fp);
	}
	printf("Done\n");
	printf("Total Lines=%d\n",nLine);

	fseek(fp,lineNumber[0],SEEK_SET);
	int mesh[3];
	fscanf(fp,"mesh: [%d,%d,%d ]\n",&mesh[0],&mesh[1],&mesh[2]);
	printf("mesh=[%d,%d,%d]\n",mesh[0],mesh[1],mesh[2]);
	int nqpoint;
	fscanf(fp,"nqpoint:%d\n",&nqpoint);
	printf("nqpoint=%d\n",nqpoint);
	int natom;

	fscanf(fp,"natom: %d\n",&natom);
	printf("natom=%d\n",natom);
		int nbranch=3*natom;
	long npr=nqpoint*nbranch;
	double pr[npr];
	double fr[npr];
	char c;
	int rec=fscanf(fp,"reciprocal_lattice%c\n",&c); 
	for(long ipr=0;ipr<npr;ipr++){
				pr[ipr]=0.0;
			}

			for(int iqp=0;iqp<nqpoint;iqp++){
				long qline=4+iqp*(3+nbranch*(3+4*natom)+1);
				if (rec>0){qline+=4;}
				for(int ibranch=0;ibranch<nbranch;ibranch++){
					long branchline=qline+3+ibranch*(3+4*natom);
					fseek(fp,lineNumber[branchline+1],SEEK_SET);
					double freq;
					fscanf(fp,"    frequency:%lf\n",&freq);
					//printf("%f\n",freq);
				//	fseek(fp,lineNumber[branchline+1],SEEK_SET);
				//	fgets(s,100,fp);
				//	printf("%d\t%s",branchline+2,s);
						for(int iatom=0;iatom<natom;iatom++){

					long atomline=branchline+3+iatom*4;
					fseek(fp,lineNumber[atomline+1],SEEK_SET);
					double e2=0;
					for(int i=0;i<3;i++){
						double a,b;
						fscanf(fp,"      - [ %lf,%lf ]\n",&a,&b);
						e2+=a*a+b*b;//printf("%d\t%f\n",ifreq,ldos[ifreq]);					
					}
					fr[ibranch+nbranch*iqp]=freq;
					pr[ibranch+nbranch*iqp]+=e2*e2;
					//fgets(s,100,fp);
					//printf("%s",s);
				}
			}
				


	}
	FILE* fpr=fopen("pr.txt","w");
		for(int ipr=0;ipr<npr;ipr++){
				pr[ipr]*=natom;
				pr[ipr]=1/pr[ipr];
				fprintf(fpr,"%f\t%f\n",fr[ipr],pr[ipr]);
			}
	//double a,b;
	//fscanf(fp,"      - [ %lf,%lf]\n",&a,&b);
	//printf("%f\t%f\n",a,b);
	return 0;
}
