#include"mathepar.h"
#include"const.h" 
#include"mcsim.h"
#include"paver.h"
#include<math.h>
#include<assert.h>
#include<stdlib.h>
#include<fstream.h>
#include<string>
#include <string.h>
#include <iostream>
#include <time.h>
int PARTYP[3];
int ngpxa;//using in INDTD and INDS   
double iold[3],ioldd[3][3];

inline int sign(double x){int i=fabs(x)/x;return i;}

#define MBIG 0x0FFFFFFF
#define MSEED 0x03C04FED
#define MZ 0
#define FAC (1.0/MBIG) 

float RanFun(long *seed)
{
	static int inext,inextp,iff=0;
	static long ma[56];
	long mj,mk;
	int i,j,k;

	if(*seed<0||iff==0)
	{ 
		iff=1;
		mj=MSEED-(*seed<0 ? -*seed : *seed);
		mj=mj%MBIG;
		ma[55]=mj;
		for(mk=1,i=1;i<=54;i++){
			j=(21*i)%55;
			ma[j]=mk;
			mk=mj-mk;
			if(mk<MZ)mk+=MBIG;
			mj=ma[j];
		}
		for(k=1;k<=4;k++)
			for(i=1;i<=55;i++){
				ma[i]-=ma[1+(i+30)%55];
				if(ma[i]<MZ)ma[i]+=MBIG;
			}
		inext=0;
		inextp=31;
		*seed=1;
	}
	if(++inext==56)inext=1;
	if(++inextp==56)inextp=1;
	mj=ma[inext]-ma[inextp];
	if(mj<MZ)mj+=MBIG;
	ma[inext]=mj;
	return mj*FAC;
}
inline double RANDNR(void)
{	
	static long reed=1;
	return RanFun(&reed);
}
inline double MAX(double p1,double p2)
{
	double a;
	a=p1>p2?p1:p2;
	return a;
}
inline double MAX(double p1,double p2,double p3)
{
	double a;
	a=p1>p2?p1:p2;
	a=a>p3?a:p3;
	return a;
}
inline double MAX(double p1,double p2,double p3,double p4)
{
	double a;
    a=p1>p2?p1:p2;
	a=a>p3?a:p3;
	a=a>p4?a:p4;
	return a;
}
inline double MAX(double p1,double p2,double p3,double p4,double p5)
{
	double a;
	a=p1>p2?p1:p2;
    a=a>p3?a:p3;
	a=a>p4?a:p4;
	a=a>p5?a:p5;
	return a;
}
inline double MIN(double p1,double p2)
{
	double a;
	a=p1<p2?p1:p2;
	return a;
}
inline double MIN(double p1,double p2,double p3)
{
	double a;
	a=p1<p2?p1:p2;
	a=a<p3?a:p3;
	return a;
}
inline double MIN(double p1,double p2,double p3,double p4)
{
	double a;
    a=p1<p2?p1:p2;
	a=a<p3?a:p3;
	a=a<p4?a:p4;
	return a;
}
inline double MIN(double p1,double p2,double p3,double p4,double p5)
{
	double a;
	a=p1<p2?p1:p2;
    a=a<p3?a:p3;
	a=a<p4?a:p4;
	a=a<p5?a:p5;
	return a;
}
inline double MIN(double p1,double p2,double p3,double p4,double p5,double p6)
{
	double a;
	a=p1<p2?p1:p2;
    a=a<p3?a:p3;
	a=a<p4?a:p4;
	a=a<p5?a:p5;
	a=a<p6?a:p6;
	return a;
}
inline int INDS(int i,int j)
{
	int ij;
	ij=ngpxa*j+i;
	return ij;
}
inline void INDTD(int ij,int &i,int &j)
{
	i=ij%ngpxa;
	j=(ij-i)/ngpxa;
	return;
}

#include"band.cpp"
#include"partical.cpp"
#include"device.cpp"
	Band bd;
	Partical par;
	DevSimulator dev;

void input(char *file)
{
ifstream ftp;
char name[6];
double vol;
int contnumber;
int i,j,k;
char tempchar[255];

ftp.open(file);

ftp>>tempchar;
ftp>>contnumber;
for(i=0;i<contnumber;i++)
{
	ftp>>name;
	ftp>>vol;
	dev.BIAST(name,vol,0,0,1E-8);
}
ftp>>tempchar;
ftp>>ndt;
ftp>>tempchar;
ftp>>outmod;
ftp>>tempchar;
ftp>>stridestat;
ftp>>tempchar;
ftp>>stridemr;
ftp>>tempchar;
ftp>>stridejb;
ftp>>tempchar;
ftp>>dt;
dt=dt/time0;
ftp>>tempchar;
ftp>>starttime;
ftp>>tempchar;
ftp>>statoutmod;
ftp>>tempchar;
ftp>>fhznumber;

for(j=0;j<MNCONT;j++)
{
	for(i=0;i<MNTIMES;i++)
	{
		currentarray[j][i]=0;
		for(k=0;k<3;k++)	currentarrayy[k][j][i]=0;
	}
	averagecurrent[i]=0;
	for(k=0;k<3;k++)	averagecurrentt[k][i]=0;
}
for(i=0;i<MNTIMES;i++)
{
	idtime[i]=0;
	didtime[i]=0;
	idtimeguiyi[i]=0;
	CI[i]=0;
	SI[i]=0;
	fhz[i]=0;
}
for(i=0;i<MNREGION;i++)
{
	evaverage[i]=0;
	hvaverage[i]=0;
	cdaverage[i]=0;
}
dev.bhscatter=0;
dev.phscatter=0;
dev.sscatter=0;
dev.oldsscatter=0;
dev.srscatter=0;
dev.spscatter=0;
dev.siscatter=0;
return;
}

char * getusedtime(void)
{
	// define 'now'. time_t is probably a typedef 
	time_t now; 
	// Get the system time and put it into 'now' as 'calender time' 
	now=time((time_t *)NULL);
	return ctime(&now);
}

void main()
{
	ofstream ftptime;
	ftptime.open("output/simulationtime.txt");
	time_t a,b;
	//start time
	ftptime<<"The simulation starts at "<<getusedtime()<<endl;
	cout<<"The simulation starts at "<<getusedtime()<<endl;
	a=time(NULL);

	masterfl=true;
	dev.READINFILE("in.txt");
	dev.INIST(par);
	bd.IELEC();
	dev.STRUC(1,"run.str");
	dev.CONFI(1,"run.con");
	ngpxa=dev.ngpx;
	dev.OUTPUTMOTRULES();
	dev.GETCUTPOSITION();
	input("input.txt");
	dev.RUN(par,bd);
//	dev.statnoise();
	dev.CONFI(2,"run_new.con");

	//end time
	ftptime<<"The simulation ends at  "<<getusedtime()<<endl;
	cout<<endl<<"The simulation ends at   "<<getusedtime()<<endl;
	b=time(NULL);
	ftptime<<"The simulation costs "<<(b-a)<<" seconds! "<<endl;
	ftptime.close();
	cout<<"The simulation costs "<<(b-a)<<" seconds! "<<endl<<endl;
	cout<<"All simulations are completed! "<<endl<<endl;

	return;
}