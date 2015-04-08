#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include <fftw3.h>

	void  getCorr(int N,double *autoCorrAtom,double *v){
		    fftw_complex* in, *out;
		    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
         out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
             		fftw_plan p;	
			p=fftw_plan_dft_1d(N,in,out,FFTW_FORWARD,FFTW_MEASURE);
			    fftw_plan p1;
    p1=fftw_plan_dft_1d(N,out,in,FFTW_BACKWARD,FFTW_MEASURE);
			for(int i=0;i<N;i++) {
				in[i][0]=v[i];
				in[i][1]=0.0;
			}
   

    fftw_execute(p);
    for(int i=0;i<N;i++){
        out[i][0]=out[i][0]*out[i][0]+out[i][1]*out[i][1];
		out[i][1]=0;
		
    }
    


    fftw_execute(p1);
    for(int i=0;i<N;i++){
    	autoCorrAtom[i]=in[i][0]/N;
    }
        fftw_destroy_plan(p);
    fftw_destroy_plan(p1);
    fftw_free(in);
    fftw_free(out);
	}

int main(int argc,char** argv){
	double factor=atof(argv[1]);
	char s[100];
	FILE* fp=fopen("jin.txt","r");
	for(int i=0;i<2;i++)fgets(s,100,fp);
	int totalStep=10000006;
	double *autoCorr=new double[totalStep];
	double *autox=new double[totalStep];
	double *autoy=new double[totalStep];
	double *autoz=new double[totalStep];
	double *vx=new double[totalStep];
	double *vy=new double[totalStep];
	double *vz=new double[totalStep];
				double *autoxAtom=new double[totalStep];
		double *autoyAtom=new double[totalStep];
		double *autozAtom=new double[totalStep];
	for(int i=0;i<totalStep;i++){
		autox[i]=0.0;
		autoy[i]=0.0;
		autoz[i]=0.0;
	}
		int step=0;
		int id;
		while(fscanf(fp,"%d%lf%lf%lf\n",&id,&vx[step],&vy[step],&vz[step])>0){
				//	if(step%100000==0)printf("%f\n",vx[step]);
			step++;
		}
		totalStep=step+1;
		//printf("%d\n",step);
		//printf("read done!\n");
		getCorr(totalStep,autoxAtom,vx);
		//getCorr(totalStep,autoyAtom,vy);
		//getCorr(totalStep,autozAtom,vz);
		for(int j=0;j<step;j++){
			autox[j]+=autoxAtom[j];
		//	autoy[j]+=autoyAtom[j];
		//	autoz[j]+=autozAtom[j];
		}
	
	FILE* f1=fopen("jcor.txt","w");
	int maxStep=int(totalStep/2+1);
		FILE* f2=fopen("jcortc.txt","w");
		double tc=0;
	for(int i=0;i<=maxStep;i++){
	//	autoCorr[i]=autox[i]+autoy[i]+autoz[i];
		tc+=autox[i]*factor;
		fprintf(f1,"%f\n",autox[i]);
		fprintf(f2,"%f\n",tc);
	}
	printf("%f\n",tc);

return 0;
}
