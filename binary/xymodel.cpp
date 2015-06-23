#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace std;
void xymodel(int l,int q);
int START=50000;
int END=100000;
int COUNT=END-START+1;
#define MAXL 129;
double *e1;
void calculate(int l,int q,double t,double *sumh,double *sumh2,double *summ,double *summ2);
double coss[72];
double sins[72];
int main(int argc, char ** argv)
{
	for(int i=0;i<72;i++){
	coss[i]=cos(i*5*3.1415/180);
	sins[i]=sin(i*5*3.1415/180);
	}
	int l=10;
	if(argc>1)l=atoi(argv[1]);
    xymodel(l,72);
    return 0;
}

void xymodel(int l,int q)
{
    double t=.01,tstart=.1,tend=1.50,deltat=.01,sumh=0,sumh2=0,c=0,summ=0,summ2=0,x=0;
    printf("#condition:(l=%d,q=%d),start simulating,please wait.\n",l,q);
    printf("temperature\tenergy\theat_capacity\tmagnetization\tsusceptibility\n");
    char s[1000]="";
   	sprintf(s,"xy_L(%d).txt",l,q);
    FILE *fp=fopen(s,"w");
    //fprintf(fp,"#resultof(L=%d,Q=%d):\n",l,q);
    fprintf(fp,"temperature\tenergy\theat_capacity\tmagnetization\tsusceptibility\n");
    e1=new double[END];
    for(t=tstart; t<=tend; t+=deltat)
        {
            deltat=(t<.6)?.05:((t<1.0)?0.005:0.05);
            calculate(l,q,t,&sumh,&sumh2,&summ,&summ2);
            c=(sumh2-sumh*sumh)/(4*t*t);
            double s=0;
            double s2=0;
            for(int i=START;i<=END;i++){
            s+=e1[i];

            }
            s/=COUNT;
            for(int i=START;i<=END;i++){
            s2+=(e1[i]-s)*(e1[i]-s);
            }
            s2/=COUNT*(4*t*t);
            x=(summ2-summ*summ)/t;
            fprintf(fp,"%f\t%f\t%f\t%f\t%f\n",t,sumh/2,c,summ,x);
            printf("%f\t%f\t%f\t%f\t%f\n",t,sumh/2,c,s2/2,x);
        }
    fclose(fp);
}
int isMK(int x,int y)
{
    return(x==0&&y==0)?1:0;
}

double e_spin(int s1,int s2){
	return -coss[abs(s1-s2)];
}
double e_ij(int at,int x,int y,int l,int a[][129]){
 	int b=y%2?1:-1;
	double e=0;
	e+=e_spin(at,a[x][(y+b+l)%l]);
	e+=e_spin(at,a[x][(y-b+l)%l]);
	e+=e_spin(at,a[(x+b+l)%l][(y-b+l)%l]);
	return e;
}
double random1()
{
    int i;
    double x=0;
    for(i=1; i<6; i++)x=((double)(rand()%10)/10+x/10);
    return x;
}
void calculate(int l,int q,double t,double *sumh,double *sumh2,double *summ,double *summ2)
{
    int a[129][129];
    int b[10];
    double Mc;
    double Ms;
    int i,j,x,y,at,count=0;
    double p,dh,dm;
    for(i=0; i<l; i++)
        for(j=0; j<l; j++)a[i][j]=0;//rand()%q;
    while(1)
        {
            x=rand()%l;
            y=rand()%l;
            at=rand()%q;
            if(at!=a[x][y])
                {
                    dh=e_ij(at,x,y,l,a)-e_ij(a[x][y],x,y,l,a);
                    p=random1();
                    if(p<exp(-dh/t)) a[x][y]=at;
                }
            else if(count<START)
                {
                    count+=isMK(x,y);
                    if(count==START)
                        {
                            (*sumh)=0;
                            (*sumh2)=0;
                            (*summ)=0;
                            (*summ2)=0;
                        }
                    continue;
                }
            count+=isMK(x,y);
            if(count>=START&&1==isMK(x,y))
                {
                    if(count==START)
                        {
                            (*sumh)=0;
                            (*sumh2)=0;
                            (*summ)=0;
                            (*summ2)=0;
                        }
                    dh=0;
                    dm=0;
                    Mc=0;
                    Ms=0;
                    //for(i=0; i<q; i++)b[i]=0;
                    //for(x=rand()%4; x<l; x+=4)
                    for(x=0; x<l; x++)
                        {
                            //for(y=rand()%4; y<l; y+=4)
                            for(y=0; y<l; y++)
                                {
                                    dh+=e_ij(a[x][y],x,y,l,a);
                                    //b[a[x][y]]++;
                                    Mc+=coss[a[x][y]];
                                    Ms+=sins[a[x][y]];
                                }
                        }
                    //for(i=1; i<q; i++)b[0]=(b[i]>b[0])?b[i]:b[0];
                    //dm=((double)(b[0]*q)/(double)(l*l)-1)/(double)(q-1);
                    dm=sqrt(Mc*Mc+Ms*Ms)/(double)(l*l);
                    double e=(dh/(double)(l*l));
                    e1[count]=e;
                    (*sumh)+=(dh/(double)(l*l));
                    (*sumh2)+=(dh*dh)/(double)(l*l*l*l);
                    (*summ)+=dm;
                    (*summ2)+=dm*dm;
                }
            if(count>=END)
                {
                    (*sumh)/=COUNT;
                    (*sumh2)/=COUNT;
                    (*summ)/=COUNT;
                    (*summ2)/=COUNT;
                    break;
                }
        }
}
