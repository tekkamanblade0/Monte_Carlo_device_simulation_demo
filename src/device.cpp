#include<iomanip.h>

void DevSimulator::READINFILE(char *file)
{
	ifstream ftp;
	int i,j;
	ftp.open(file);

	ftp>>tempchar;
	ftp>>material;

	if(material!=0&&material!=1)
	{
		cout<<"This kind of material can't be simulated by this program!!"<<endl;
		exit(0);
	}
	if(material==0)	cout<<"------******   Silicon	"<<endl;
	else if(material==1)	cout<<"------******   Germanium   "<<endl;
	//read in direction
	ftp>>tempchar;
	ftp>>direction;
	if(direction!=100&&direction!=110)
	{
		cout<<"This kind of direction can't be simulated by this program!!"<<endl;
	}
	if(direction==100)	cout<<"------******   (100)   "<<endl;
	else if(direction==110)	cout<<"------******   (110)   "<<endl;
	//read in surface_scattering
	ftp>>tempchar;
	ftp>>tempchar;
	ftp>>tempdifprelectron;
	ftp>>tempchar;
	ftp>>tempdifprhole;
	//read in quantum effect
	//broaden effect
	ftp>>tempchar;
	ftp>>tempchar;
	ftp>>quantumeffectcbchar;
	ftp>>quantumeffectqpchar;
	if(quantumeffectcbchar!=1&&quantumeffectcbchar!=0)
	{
		cout<<"Quantum effect of broaden effect is wrong!!"<<endl;
		exit(0);
	}
	if(quantumeffectcbchar==1)	quantumeffectcbfl=true;
	else if(quantumeffectcbchar==0)	quantumeffectcbfl=false;
	if(quantumeffectcbfl)	cout<<"------******   Quantum effect of cb             YES   "<<endl;
	else	cout<<"------******   Quantum effect of cb             NO    "<<endl;
	cout<<endl<<endl;
	//quantum potential
	if(quantumeffectqpchar!=1&&quantumeffectqpchar!=0&&quantumeffectqpchar!=2&&quantumeffectqpchar!=3&&quantumeffectqpchar!=4)
	{
		cout<<"Quantum effect of quantum potential is wrong!!"<<endl;
		exit(0);
	}
	if(quantumeffectqpchar==1)	quantumeffectqpfl=true;
	else if(quantumeffectqpchar==0)	quantumeffectqpfl=false;
	if(quantumeffectqpchar==2)	quantumeffectbmfl=true;
	else if(quantumeffectqpchar==0)	quantumeffectbmfl=false;
	if(quantumeffectqpchar==3)	quantumeffectqnfl=true;
	else if(quantumeffectqpchar==0)	quantumeffectqnfl=false;
	if(quantumeffectqpchar==4)	quantumeffectfmfl=true;
	else						quantumeffectfmfl=false;
	if(quantumeffectqpfl)	cout<<"------******   Quantum effect of Potential                       YES   "<<endl;
	else	cout<<"------******   Quantum effect of Potential                       NO    "<<endl;
	if(quantumeffectbmfl)	cout<<"------******   Quantum effect of Bohm Potential                  YES   "<<endl;
	else	cout<<"------******   Quantum effect of Bohm Potential                  NO    "<<endl;
	if(quantumeffectqnfl)	cout<<"------******   Quantum effect of Ndens                           YES   "<<endl;
	else	cout<<"------******   Quantum effect of Ndens                           NO    "<<endl;
	if(quantumeffectfmfl)	cout<<"------******   Quantum effect of Effective Potential             YES   "<<endl;
	else	cout<<"------******   Quantum effect of Effective Potential             NO    "<<endl;
	cout<<endl<<endl;
	cout<<"------difpr-electron                  "<<tempdifprelectron<<endl<<endl;
	cout<<"------difpr-hole                      "<<tempdifprhole<<endl<<endl;

	//read in cutregion
	ftp>>tempchar;//cutregion
	ftp>>tempchar;
	ftp>>cutchar;
	if(cutchar==1)	cutfl=true;
	else if(cutchar==0)	cutfl=false;
	else cout<<"Cut flag for getting Vxy and Energy is error!!"<<endl;
	if(cutfl)	cout<<"------******   Cut                         	YES   "<<endl<<endl;
	else	cout<<"------******   Cut                         	NO    "<<endl<<endl;
	ftp>>tempchar;//xcut
	ftp>>xcutnumber;
	for(i=0;i<xcutnumber;i++)	ftp>>positionx[i];
	for(i=0;i<xcutnumber;i++)	ftp>>positionx[i];
	ftp>>tempchar;//ycut
	ftp>>ycutnumber;
	for(i=0;i<ycutnumber;i++)	ftp>>positiony[i];
	for(i=0;i<ycutnumber;i++)	ftp>>positiony[i];

	//read in mindopforcoulomb
	ftp>>tempchar;
	ftp>>mindopforcoulomb;
	cout<<endl<<"------******min dope for coulomb               "<<mindopforcoulomb<<endl;
	//read in minelectronenergy
	ftp>>tempchar;
	ftp>>minelectronenergy;
	cout<<endl<<"------******min energy for coulomb             "<<minelectronenergy<<endl;

	//read in noise stat region
	ftp>>tempchar;
	ftp>>tempchar;
	ftp>>statnoiseregionnumber;
	ftp>>tempchar;
	ftp>>tempchar;
	ftp>>tempchar;
	ftp>>tempchar;
	for(i=0;i<statnoiseregionnumber;i++)
	{
		ftp>>tempchar;
		ftp>>xminnoise[i];
		ftp>>xmaxnoise[i];
		ftp>>yminnoise[i];
		ftp>>ymaxnoise[i];
	}

	//read in surface scatter region
	ftp>>tempchar;
	ftp>>tempchar;
	ftp>>tempsssfl;
	if(tempsssfl==1)	sssfl=true;
	else if(tempsssfl==0)	sssfl=false;
	cout<<endl<<endl;
	if(sssfl)
	cout<<"------******surface scattering                 YES"<<endl;
	else
	cout<<"------******surface scattering                 NO"<<endl;
	ftp>>tempchar;//xmin1
	ftp>>xminss1;
	ftp>>tempchar;//xmax1
	ftp>>xmaxss1;
	ftp>>tempchar;//ymin1
	ftp>>yminss1;
	ftp>>tempchar;//ymax1
	ftp>>ymaxss1;
	ftp>>tempchar;//xmin2
	ftp>>xminss2;
	ftp>>tempchar;//xmax2
	ftp>>xmaxss2;
	ftp>>tempchar;//ymin2
	ftp>>yminss2;
	ftp>>tempchar;//ymax2
	ftp>>ymaxss2;
	ftp>>tempchar;//surface number
	ftp>>surfacenumber;
	if(surfacenumber!=1&&surfacenumber!=2)
	{
		cout<<"------******wrong surface number!!"<<endl;
		exit(0);
	}
	ftp>>tempchar;//surface position 1
	ftp>>surfaceposition1;
	ftp>>tempchar;//surface position 2
	ftp>>surfaceposition2;

	ftp>>tempchar;
	ftp>>ail;
	ftp>>tempchar;
	ftp>>delta;

	cout<<endl<<"------ail                                   "<<ail<<endl;
	cout<<endl<<"------delta                                 "<<delta<<endl;

	ftp>>tempchar;
	ftp>>tempsrfl;
	srfl=false;
	if(tempsrfl==1)	srfl=true;
	else if(tempsrfl==0)	srfl=false;

	if(srfl)
	cout<<endl<<"------******Surface Roughness                  YES"<<endl;
	else
	cout<<endl<<"------******Surface Roughness                  NO"<<endl;

	ftp>>tempchar;
	ftp>>tempspfl;
	spfl=false;
	if(tempspfl==1)	spfl=true;
	else if(tempspfl==0)	spfl=false;

	if(spfl)
	cout<<endl<<"------******Surface Phonon                     YES"<<endl;
	else
	cout<<endl<<"------******Surface Phonon                     NO"<<endl;

	ftp>>tempchar;
	ftp>>tempscfl;
	scfl=false;
	if(tempscfl==1)	scfl=true;
	else if(tempscfl==0)	scfl=false;

	if(scfl)
	cout<<endl<<"------******Surface Coulomb                    YES"<<endl;
	else
	cout<<endl<<"------******Surface Coulomb                    NO"<<endl;

	ftp>>tempchar;
	ftp>>tempballisticfl;
	ballisticfl=false;
	if(tempballisticfl==1)	ballisticfl=true;
	else if(tempballisticfl==0)	ballisticfl=false;

	if(ballisticfl)
	cout<<endl<<"------******Ballistic                          YES"<<endl;
	else
	cout<<endl<<"------******Ballistic                          NO"<<endl;

	ftp>>tempchar>>tempselfconsistentchar;
	if(tempselfconsistentchar==1)	selfconsistentfl=true;
	else if(tempselfconsistentchar==0)	selfconsistentfl=false;

	if(selfconsistentfl)
	cout<<endl<<"------******Self consistent simulation         YES"<<endl;
	else
	cout<<endl<<"------******Self consistent simulation         NO"<<endl;

	ftp>>tempchar>>tempbtbtchar;
	if(tempbtbtchar==1)	btbtfl=true;
	else if(tempbtbtchar==0)	btbtfl=false;

	if(btbtfl)
	cout<<endl<<"------******Band to band tunneling             YES"<<endl;
	else
	cout<<endl<<"------******Band to band tunneling             NO"<<endl;

	ftp>>tempchar>>tempmethodforgetisedatachar;
	if(tempmethodforgetisedatachar==1)	methodforgetisedatafl=true;
	else if(tempmethodforgetisedatachar==0)	methodforgetisedatafl=false;

	if(methodforgetisedatafl)
	cout<<endl<<"------******Method for get data from ise       Rectangle"<<endl;
	else
	cout<<endl<<"------******Method for get data from ise       Four points"<<endl;

	ftp>>tempchar;
	ftp>>rr;

	cout<<endl<<"------******Recombination_constant_parameter_r "<<rr<<endl;

	ftp>>tempchar;
	ftp>>tempchar>>tempchar>>tempchar>>tempchar;

	ftp>>ifminx>>ifmaxx>>ifminy>>ifmaxy;
	cout<<endl<<"------Interface region of X min             "<<ifminx<<endl;
	cout<<endl<<"------Interface region of X max             "<<ifmaxx<<endl;
	cout<<endl<<"------Interface region of Y min             "<<ifminy<<endl;
	cout<<endl<<"------Interface region of Y max             "<<ifmaxy<<endl;

	ftp.close();
	cout<<endl<<"------OK!"<<endl<<endl<<endl;

	ftp.open("fermi.txt");
	nfpoint=1000;
	ftp>>nfpoint;
	for(i=0;i<nfpoint;i++)
	{
		ftp>>fermi[0][i]>>fermi[1][i];
	}
	ftp.close();

	for(i=0;i<MNREGION;i++)
	{
		for(j=0;j<MNTIMES;j++)
		{
			ev[i][j]=0;
			evsum[i][j]=0;
			evnumber[i][j]=0;
			hv[i][j]=0;
			hvsum[i][j]=0;
			hvnumber[i][j]=0;
			cd[i][j]=0;
			SEV[i][j]=0;
			SHV[i][j]=0;
			SCD[i][j]=0;
			CEV[i][j]=0;
			CHV[i][j]=0;
			CCD[i][j]=0;
		}
	}

	return;
}

void DevSimulator::GETISEDATA(void)
{
//  m为isegridnumber,n为输入的pot点number，实际后面用到的主要是isegridnumber也就是m;
//	double isegridx1[MNGP],isegridy1[MNGP];//isepot1[MNGP],isendens1[MNGP],iseef1[MNGP];
	cout<<"------Reading ISE data into this program!!"<<endl;

	int	order[100];

	int datanumber1[100],datanumber2[100];
	
	int isegridnumber1;
	int i,j,k,ii,jj,ij,m,n;
	int regionnb;
	int tempcharlength;
	char singlechar;
	int	charnumber;
	char gridchar[255],datachar[255];

	bool readfl;
	readfl=true;

	ifstream ftp;
	ftp.open("isegrid.txt");
	for(i=0;i<1000;i++)	
	{
		ftp>>tempchar;
		if(i==3)	ftp>>gridchar;
		if(strcmp(tempchar,"nb_regions")==0)
		{
			ftp>>tempchar;
			ftp>>regionnb;
			break;
		}
	}
	for(i=0;i<36+2*regionnb;i++)	ftp>>tempchar;
	ftp>>isegridnumber;
	ftp>>tempchar;
	ftp>>tempchar;
	m=isegridnumber;
//读取ise的格点坐标
	for(ij=0;ij<m;ij++)
	{
		ftp>>isegridy[ij];
		isegridy[ij]=isegridy[ij]*1e-6/spr0;
		ftp>>isegridx[ij];
		isegridx[ij]=isegridx[ij]*1e-6/spr0;
	}
	ftp.close();

	ifstream ftp1;
	ftp1.open("isedata.txt");
	ftp1>>tempchar;
	ftp1>>tempchar;
	ftp1>>singlechar;
	ftp1>>tempchar;
	if(strcmp(tempchar,gridchar)!=0)
	{
		cout<<"The isegrid and isedata file is not the same structure!!";
		exit(0);
	}
	for(i=0;i<4;i++)	ftp1>>tempchar;
	i=0;
	order[0]=0;
	order[1]=1;
	while(strcmp(tempchar,"TEXT")!=0)
	{
		ftp1>>tempchar;
		tempcharlength=strlen(tempchar);
		for(j=1;j<tempcharlength-1;j++)	chartemp[j-1]=tempchar[j];
		if		(strcmp(chartemp,"X")==0)							order[0]=i;//Xnb=i;
		else if	(strcmp(chartemp,"Y")==0)							order[1]=i;//Ynb=i;
		else if	(strcmp(chartemp,"ElectrostaticPotential")==0)		order[2]=i;//ElectrostaticPotential=i;
		else if	(strcmp(chartemp,"eDensity")==0)					order[3]=i;//eDensity=i;
		else if	(strcmp(chartemp,"hDensity")==0)					order[4]=i;//hDensity=i;
		else if	(strcmp(chartemp,"ElectricField-X")==0)				order[5]=i;//ElectricFieldX=i;
		else if	(strcmp(chartemp,"ElectricField-Y")==0)				order[6]=i;//ElectricFieldY=i;
		else if	(strcmp(chartemp,"Abs(ElectricField)")==0)			order[7]=i;//ElectricField=i;
		else if	(strcmp(chartemp,"DopingConcentration")==0)			order[8]=i;//DopingConcentration=i;
		else if	(strcmp(chartemp,"DonorConcentration")==0)			order[9]=i;//DonorConcentration=i;
		else if	(strcmp(chartemp,"AcceptorConcentration")==0)		order[10]=i;//AcceptorConcentration=i;
		else if	(strcmp(chartemp,"eVelocity-X")==0)					order[11]=i;//eVelocityX=i;
		else if	(strcmp(chartemp,"eVelocity-Y")==0)					order[12]=i;//eVelocityY=i;
		else if	(strcmp(chartemp,"Abs(eVelocity)")==0)				order[13]=i;//eVelocity=i;
		else if	(strcmp(chartemp,"hVelocity-X")==0)					order[14]=i;//hVelocityX=i;
		else if	(strcmp(chartemp,"hVelocity-Y")==0)					order[15]=i;//hVelocityY=i;
		else if	(strcmp(chartemp,"Abs(hVelocity)")==0)				order[16]=i;//hVelocity=i;
		else if	(strcmp(chartemp,"eQuasiFermiPotential")==0)		order[17]=i;//hVelocity=i;
		else if	(strcmp(chartemp,"hQuasiFermiPotential")==0)		order[18]=i;//hVelocity=i;
		else if	(strcmp(chartemp,"ConductionBandEnergy")==0)		order[19]=i;//hVelocity=i;
		else if	(strcmp(chartemp,"ValenceBandEnergy")==0)			order[20]=i;//hVelocity=i;		
		else	i=i;

		for(j=0;j<255;j++)	chartemp[j]=NULL;
		i++;
		charnumber=i-1;
	}
	for(i=0;i<18;i++)	ftp1>>tempchar;
	i=0;
	int timesnumber;
	timesnumber=0;
	while(readfl)
	{
		for(j=0;j<2;j++)	ftp1>>tempchar;
		for(j=0;j<2;j++)	ftp1>>singlechar;
		ftp1>>datanumber1[i];
		for(j=0;j<3;j++)	ftp1>>singlechar;
		ftp1>>datanumber2[i];
		ftp1>>tempchar;
		ftp1>>tempchar;
		if(strcmp(tempchar,"ET=Quadrilateral")==0)	readfl=true;
		else if(strcmp(tempchar,"ET=Triangle")==0)	readfl=false;
		for(j=0;j<charnumber+1;j++)	ftp1>>tempchar;
		for(k=timesnumber;k<timesnumber+datanumber1[i];k++)
		{
			for(ii=0;ii<charnumber;ii++)
			{
				if		(order[0]==ii)	{ftp1>>isegridy1[k];
										 isegridy1[k]=isegridy1[k]*1e-6/spr0;}
				else if	(order[1]==ii)	{ftp1>>isegridx1[k];
										 isegridx1[k]=isegridx1[k]*1e-6/spr0;}
				else if	(order[2]==ii)	{ftp1>>isepot1[k];
										 isepot1[k]=isepot1[k]/pot0;}
				else if	(order[3]==ii)	ftp1>>iseedens1[k];
				else if	(order[4]==ii)	ftp1>>isehdens1[k];
				else if	(order[5]==ii)	{ftp1>>iseefy1[k];
										 iseefy1[k]=iseefy1[k]*100/field0;}
				else if	(order[6]==ii)	{ftp1>>iseefx1[k];
										 iseefx1[k]=iseefx1[k]*100/field0;}
				else if	(order[7]==ii)	ftp1>>iseef1[k];
				else if	(order[8]==ii)	ftp1>>isedopeconcentration1[k];
				else if	(order[9]==ii)	ftp1>>isedonorconcentration1[k];
				else if	(order[10]==ii)	ftp1>>iseaccepconcentration1[k];
				else if	(order[11]==ii)	ftp1>>iseevy1[k];
				else if	(order[12]==ii)	ftp1>>iseevx1[k];
				else if	(order[13]==ii)	ftp1>>iseev1[k];
				else if	(order[14]==ii)	ftp1>>isehvy1[k];
				else if	(order[15]==ii)	ftp1>>isehvx1[k];
				else if	(order[16]==ii)	ftp1>>isehv1[k];
				else if (order[17]==ii)	{ftp1>>iseEfn1[k];
										 iseEfn1[k]*=-1.0;}
				else if (order[18]==ii)	{ftp1>>iseEfp1[k];
										 iseEfp1[k]*=-1.0;}
				else if (order[19]==ii)	ftp1>>iseEc1[k];
				else if (order[20]==ii)	ftp1>>iseEv1[k];
				else					ftp1>>tempdouble;
			}
		}
		for(k=0;k<datanumber2[i];k++)
		{
			for(jj=0;jj<4;jj++)	ftp1>>tempint;
		}
		timesnumber+=datanumber1[i];
		i++;
	}
	isegridnumber1=0;
	for(k=0;k<i-1;k++)	isegridnumber1+=datanumber1[k];
	n=isegridnumber1;
//读取ise的网格格点坐标以及其电势值
//通过对比坐标，赋予各个格点相应的电势。
	int eqnumber;
	eqnumber=0;
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			if((fabs(isegridx[i]-isegridx1[j])<=1e-10)&&(fabs(isegridy[i]-isegridy1[j])<=1e-10))	
			{
				eqnumber++;
				isepot[i]=isepot1[j];
				iseedens[i]=iseedens1[j];
				isehdens[i]=isehdens1[j];
				iseefx[i]=iseefx1[j];
				iseefy[i]=iseefy1[j];
				iseevx[i]=iseevx1[j];
				iseevy[i]=iseevy1[j];
				iseev[i]=iseev1[j];
				isehvx[i]=isehvx1[j];
				isehvy[i]=isehvy1[j];
				isehv[i]=isehv1[j];
				iseEfn[i]=iseEfn1[j];
				iseEfp[i]=iseEfp1[j];
				iseEc[i]=iseEc1[j];
				iseEv[i]=iseEv1[j];
				break;
			}
		}
	}
//输出电势赋予后的各个格点的电势
/*
	ofstream ftp2;
	ftp2.open("output/isedataoutput.txt");
	for(i=0;i<m;i++)
	{
		ftp2<<isegridx[i]<<" "<<isegridy[i]<<" "<<isepot[i]<<" "<<iseedens[i]
		    <<" "<<isehdens[i]<<" "<<iseefx[i]<<" "<<iseefy[i]<<" "<<iseev[i]<<" "<<isehv[i]
			<<" "<<iseEfn[i]<<" "<<iseEfp[i]<<" "<<iseEc[i]<<" "<<iseEv[i]<<endl;
	}
	ftp2.close();
*/
	cout<<"------OK!"<<endl<<endl;
}

void DevSimulator::DATAFROMISETOISEGOOD(void)
{
	cout<<"------Transfering normal ISE data to good!!"<<endl;
	int i,j,ij,k,m;
	// translate isepot to well-regulated isepot;				
//	double x[4],y[4];
	double tempgridx,tempgridy;
	double temppot,tempedens,temphdens,tempefx,tempefy,tempEfn,tempEfp,tempEc,tempEv;
		m=isegridnumber;
		if(m>MNGP)
		{
			cout<<"there is no enough grid number in x and y !!!!"<<endl;
			exit(0);
		}
//从左上角开始排列各点
//		1	5	9	13...
//		2	6	10	14...
//		3	7	11	15...
//		4	8	12	16...

		for(ij=0;ij<m;ij++)
		{
			mingridx1[ij]=isegridx[ij];
			mingridy1[ij]=isegridy[ij];
			minpot1[ij]=isepot[ij];
			minedens1[ij]=iseedens[ij];
			minhdens1[ij]=isehdens[ij];
			minefx1[ij]=iseefx[ij];
			minefy1[ij]=iseefy[ij];
			minEfn1[ij]=iseEfn[ij];
			minEfp1[ij]=iseEfp[ij];
			minEc1[ij]=iseEc[ij];
			minEv1[ij]=iseEv[ij];
		}

		for(i=0;i<m-1;i++)
		{
			k=i;
			for(j=i+1;j<m;j++)
			{
				if(mingridy1[j]<mingridy1[k])
				{
						k=j;
				}
			}
			tempgridx=mingridx1[k];mingridx1[k]=mingridx1[i];mingridx1[i]=tempgridx;
			tempgridy=mingridy1[k];mingridy1[k]=mingridy1[i];mingridy1[i]=tempgridy;
			temppot=minpot1[k];minpot1[k]=minpot1[i];minpot1[i]=temppot;
			tempedens=minedens1[k];minedens1[k]=minedens1[i];minedens1[i]=tempedens;
			temphdens=minhdens1[k];minhdens1[k]=minhdens1[i];minhdens1[i]=temphdens;
			tempefx=minefx1[k];minefx1[k]=minefx1[i];minefx1[i]=tempefx;
			tempefy=minefy1[k];minefy1[k]=minefy1[i];minefy1[i]=tempefy;
			tempEfn=minEfn1[k];minEfn1[k]=minEfn1[i];minEfn1[i]=tempEfn;
			tempEfp=minEfp1[k];minEfp1[k]=minEfp1[i];minEfp1[i]=tempEfp;
			tempEc=minEc1[k];minEc1[k]=minEc1[i];minEc1[i]=tempEc;
			tempEv=minEv1[k];minEv1[k]=minEv1[i];minEv1[i]=tempEv;
		}

		for(i=0;i<m-1;i++)
		{
			k=i;
			for(j=i+1;j<m;j++)
			{
				if(mingridy1[j]==mingridy1[k])
				{
					if(mingridx1[j]<mingridx1[k])
					{
						k=j;
					}
				}
			}
			tempgridx=mingridx1[k];mingridx1[k]=mingridx1[i];mingridx1[i]=tempgridx;
			tempgridy=mingridy1[k];mingridy1[k]=mingridy1[i];mingridy1[i]=tempgridy;
			temppot=minpot1[k];minpot1[k]=minpot1[i];minpot1[i]=temppot;
			tempedens=minedens1[k];minedens1[k]=minedens1[i];minedens1[i]=tempedens;
			temphdens=minhdens1[k];minhdens1[k]=minhdens1[i];minhdens1[i]=temphdens;
			tempefx=minefx1[k];minefx1[k]=minefx1[i];minefx1[i]=tempefx;
			tempefy=minefy1[k];minefy1[k]=minefy1[i];minefy1[i]=tempefy;
			tempEfn=minEfn1[k];minEfn1[k]=minEfn1[i];minEfn1[i]=tempEfn;
			tempEfp=minEfp1[k];minEfp1[k]=minEfp1[i];minEfp1[i]=tempEfp;
			tempEc=minEc1[k];minEc1[k]=minEc1[i];minEc1[i]=tempEc;
			tempEv=minEv1[k];minEv1[k]=minEv1[i];minEv1[i]=tempEv;
		}
/*
		ofstream ftp1;
		ftp1.open("output/data1out.txt");
		for(i=0;i<m;i++)
		{
			ftp1<<mingridx1[i]<<" "<<mingridy1[i]<<" "<<minpot1[i]<<endl;
		}
		ftp1.close();
*/
//从右下角开始排列各点
//		16	12	8	4...
//		15	11	7	3...
//		14	10	6	2...
//		13	9	5	1...

		for(ij=0;ij<m;ij++)
		{
			mingridx2[ij]=mingridx1[m-ij-1];
			mingridy2[ij]=mingridy1[m-ij-1];
			minpot2[ij]=minpot1[m-ij-1];
			minedens2[ij]=minedens1[m-ij-1];
			minhdens2[ij]=minhdens1[m-ij-1];
			minefx2[ij]=minefx1[m-ij-1];
			minefy2[ij]=minefy2[m-ij-1];
			minEfn2[ij]=minEfn1[m-ij-1];
			minEfp2[ij]=minEfp1[m-ij-1];
			minEc2[ij]=minEc1[m-ij-1];
			minEv2[ij]=minEv1[m-ij-1];
		}
/*
		ofstream ftp2;
		ftp2.open("output/data2out.txt");
		for(i=0;i<m;i++)
		{
			ftp2<<mingridx2[i]<<" "<<mingridy2[i]<<" "<<minpot2[i]<<endl;
		}
		ftp2.close();
*/

//从左下角开始排列各点
//		13	14	15	16...
//		9	10	11	12...
//		5	6	7	8...
//		1	2	3	4...
		for(ij=0;ij<m;ij++)
		{
			mingridx3[ij]=isegridx[ij];
			mingridy3[ij]=isegridy[ij];
			minpot3[ij]=isepot[ij];
			minedens3[ij]=iseedens[ij];
			minhdens3[ij]=isehdens[ij];
			minefx3[ij]=iseefx[ij];
			minefy3[ij]=iseefy[ij];
			minEfn3[ij]=iseEfn[ij];
			minEfp3[ij]=iseEfp[ij];
			minEc3[ij]=iseEc[ij];
			minEv3[ij]=iseEv[ij];
		}

		for(i=0;i<m-1;i++)
		{
			k=i;
			for(j=i+1;j<m;j++)
			{
				if(mingridx3[j]>mingridx3[k])
				{
						k=j;
				}
			}
			tempgridx=mingridx3[k];mingridx3[k]=mingridx3[i];mingridx3[i]=tempgridx;
			tempgridy=mingridy3[k];mingridy3[k]=mingridy3[i];mingridy3[i]=tempgridy;
			temppot=minpot3[k];minpot3[k]=minpot3[i];minpot3[i]=temppot;
			tempedens=minedens3[k];minedens3[k]=minedens3[i];minedens3[i]=tempedens;
			temphdens=minhdens3[k];minhdens3[k]=minhdens3[i];minhdens3[i]=temphdens;
			tempefx=minefx3[k];minefx3[k]=minefx3[i];minefx3[i]=tempefx;
			tempefy=minefy3[k];minefy3[k]=minefy3[i];minefy3[i]=tempefy;
			tempEfn=minEfn3[k];minEfn3[k]=minEfn3[i];minEfn3[i]=tempEfn;
			tempEfp=minEfp3[k];minEfp3[k]=minEfp3[i];minEfp3[i]=tempEfp;
			tempEc=minEc3[k];minEc3[k]=minEc3[i];minEc3[i]=tempEc;
			tempEv=minEv3[k];minEv3[k]=minEv3[i];minEv3[i]=tempEv;
		}

		for(i=0;i<m-1;i++)
		{
			k=i;
			for(j=i+1;j<m;j++)
			{
				if(mingridx3[j]==mingridx3[k])
				{
					if(mingridy3[j]<mingridy3[k]) 
					{
						k=j;
					}
				}
			}
			tempgridx=mingridx3[k];mingridx3[k]=mingridx3[i];mingridx3[i]=tempgridx;
			tempgridy=mingridy3[k];mingridy3[k]=mingridy3[i];mingridy3[i]=tempgridy;
			temppot=minpot3[k];minpot3[k]=minpot3[i];minpot3[i]=temppot;
			tempedens=minedens3[k];minedens3[k]=minedens3[i];minedens3[i]=tempedens;
			temphdens=minhdens3[k];minhdens3[k]=minhdens3[i];minhdens3[i]=temphdens;
			tempefx=minefx3[k];minefx3[k]=minefx3[i];minefx3[i]=tempefx;
			tempefy=minefy3[k];minefy3[k]=minefy3[i];minefy3[i]=tempefy;
			tempEfn=minEfn3[k];minEfn3[k]=minEfn3[i];minEfn3[i]=tempEfn;
			tempEfp=minEfp3[k];minEfp3[k]=minEfp3[i];minEfp3[i]=tempEfp;
			tempEc=minEc3[k];minEc3[k]=minEc3[i];minEc3[i]=tempEc;
			tempEv=minEv3[k];minEv3[k]=minEv3[i];minEv3[i]=tempEv;
		}
/*
		ofstream ftp3;
		ftp3.open("output/data3out.txt");
		for(i=0;i<m;i++)
		{
			ftp3<<mingridx3[i]<<" "<<mingridy3[i]<<" "<<minpot3[i]<<endl;
		}
		ftp3.close();
*/
//从右上角开始排列各点
//		4	3	2	1...
//		8	7	6	5...
//		12	11	10	9...
//		16	15	14	13...
		for(ij=0;ij<m;ij++)
		{
			mingridx4[ij]=mingridx3[m-ij-1];
			mingridy4[ij]=mingridy3[m-ij-1];
			minpot4[ij]=minpot3[m-ij-1];
			minedens4[ij]=minedens3[m-ij-1];
			minhdens4[ij]=minhdens3[m-ij-1];
			minefx4[ij]=minefx3[m-ij-1];
			minefy4[ij]=minefy3[m-ij-1];
			minEfn4[ij]=minEfn3[m-ij-1];
			minEfp4[ij]=minEfp3[m-ij-1];
			minEc4[ij]=minEc3[m-ij-1];
			minEv4[ij]=minEv3[m-ij-1];
		}
/* 
		ofstream ftp4;
		ftp4.open("output/data4out.txt");
		for(i=0;i<m;i++)
		{
			ftp4<<mingridx4[i]<<" "<<mingridy4[i]<<" "<<minpot4[i]<<endl;
		}
		ftp4.close();
*/

//得到各种排列下的排数，以及各排对应的端点

		//从左上角开始的排列
		minn1=0;
		bpx1up[minn1]=mingridx1[minn1];
		for(i=0;i<m;i++)
		{
			if(mingridy1[i+1]!=mingridy1[i])
			{
				bpx1down[minn1]=mingridx1[i];
				bpx1up[minn1+1]=mingridx1[i+1];
				bpy1[minn1]=mingridy1[i];
				minn1++;
			}

		}
		bpx1down[minn1]=mingridx1[m];
/*
		ofstream ftp5;
		ftp5.open("output/1.txt");
		ftp5<<minn1<<endl;
		for(i=0;i<minn1;i++)
		{
			ftp5<<bpx1up[i]<<" "<<bpx1down[i]<<" "<<bpy1[i]<<endl;
		}
		ftp5.close();
*/
		//从右下角开始的排列
		minn2=0;
		bpx2down[minn2]=mingridx2[minn2];
		for(i=0;i<m;i++)
		{
			if(mingridy2[i+1]!=mingridy2[i])
			{
				bpx2up[minn2]=mingridx2[i];
				bpx2down[minn2+1]=mingridx2[i+1];
				bpy2[minn2]=mingridy2[i];
				minn2++;
			}

		}
		bpx2up[minn2]=mingridx2[m];
/*
		ofstream ftp6;
		ftp6.open("output/2.txt");
		ftp6<<minn2<<endl;
		for(i=0;i<minn2;i++)
		{
			ftp6<<bpx2down[i]<<" "<<bpx2up[i]<<" "<<bpy2[i]<<endl;
		}
		ftp6.close();
*/
		//从左下角开始的排列
		minn3=0;
		bpy3left[minn3]=mingridy3[minn3];
		for(i=0;i<m;i++)
		{
			if(mingridx3[i+1]!=mingridx3[i])
			{
				bpy3right[minn3]=mingridy3[i];
				bpy3left[minn3+1]=mingridy3[i+1];
				bpx3[minn3]=mingridx3[i];
				minn3++;
			}

		}
		bpy3right[minn3]=mingridy3[m];
/*
		ofstream ftp7;
		ftp7.open("output/3.txt");
		ftp7<<minn3<<endl;
		for(i=0;i<minn3;i++)
		{
			ftp7<<bpy3left[i]<<" "<<bpy3right[i]<<" "<<bpx3[i]<<endl;
		}
		ftp7.close();
*/
		//从右上角开始的排列
		minn4=0;
		bpy4right[minn4]=mingridy4[minn4];
		for(i=0;i<m;i++)
		{
			if(mingridx4[i+1]!=mingridx4[i])
			{
				bpy4left[minn4]=mingridy4[i];
				bpy4right[minn4+1]=mingridy4[i+1];
				bpx4[minn4]=mingridx4[i];
				minn4++;
			}

		}
		bpy4left[minn4]=mingridy4[m];
/*
		ofstream ftp8;
		ftp8.open("output/4.txt");
		ftp8<<minn4<<endl;
		for(i=0;i<minn3;i++)
		{
			ftp8<<bpy4right[i]<<" "<<bpy4left[i]<<" "<<bpx4[i]<<endl;
		}
		ftp8.close();
*/
/*
		ofstream ftp111;
		ftp111.open("output/isenumber.txt");
		ftp111<<"y minn1  "<<"y minn2  "<<"x minn3  "<<"x minn4"<<endl;
		ftp111<<minn1<<"  "<<minn2<<"  "<<minn3<<"  "<<minn4<<endl;
		ftp111.close();
*/

		cout<<"------OK!"<<endl<<endl;
return;
}

void DevSimulator::DATAFROMISETOMC()
{
	cout<<"------Transfering good ISE data to Monte Carlo data!!"<<endl;
	int i,j,ii,jj;
	int k,l,m,n,p;
// Method fro getting data from ise is rectangle!
if(methodforgetisedatafl)
{
	int singlenumber,bordernumber,rectanglenumber;
	int nn1,nn2,nn3;
	nn1=0;
	nn2=0;
	nn3=0;

	singlenumber=0;
	bordernumber=0;
	rectanglenumber=0;
	bool singlepoint,borderpoint;
	bool borderpointupdown,borderpointleftright;
	bool voidpoint;
	double xx[4],yy[4];
	double temppot[4];
	double tempedens[4];
	double temphdens[4];
	double tempefx[4];
	double tempefy[4];
	double tempEfn[4];
	double tempEfp[4];
	double tempEc[4];
	double tempEv[4];
	double potup,potdown,potleft,potright;
	double edensup,edensdown,edensleft,edensright;
	double hdensup,hdensdown,hdensleft,hdensright;
	double efxup,efxdown,efxleft,efxright;
	double efyup,efydown,efyleft,efyright;
	double Efnup,Efndown,Efnleft,Efnright;
	double Efpup,Efpdown,Efpleft,Efpright;
	double Ecup,Ecdown,Ecleft,Ecright;
	double Evup,Evdown,Evleft,Evright;

	for(ii=0;ii<ngpx-1;ii++)
	{
		gridxx[ii]=(gridx[ii]+gridx[ii+1])/2.0;
	}
	for(jj=0;jj<ngpy-1;jj++)
	{
		gridyy[jj]=(gridy[jj]+gridy[jj+1])/2.0;
	}

	m=isegridnumber;

	for(jj=0;jj<ngpy;jj++)
	for(ii=0;ii<ngpx;ii++)
	{
		minflag1=false;
		minflag2=false;
		minflag3=false;
		minflag4=false;
		singlepoint=false;
		borderpoint=false;
		borderpointupdown=false;
		borderpointleftright=false;

		//验证是否与某个点重合;
		//因为本身选取的是Monte Carlo网格内中心点的坐标，而且排除了最下边和最右边
		//的两排网格，所以现在剩下的网格内中心点的坐标肯定在世纪的器件内部，所以
		//如果和某个ISE的网格格点重合的话，也不会在最外面，即不会在世纪器件的边界
		//点上，只需要判断是在那条边上即可。
		if(!singlepoint)
		{
			for(i=0;i<m;i++)
			{
				//先判断是否与minn1中的某个点重合(第i个点)
				if(((gridxx[ii]-mingridx1[i])*(gridxx[ii]-mingridx1[i])+(gridyy[jj]-mingridy1[i])*(gridyy[jj]-mingridy1[i]))<1e-20)
				{
					n=i;
					singlenumber++;
					singlepoint=true;
					break;
				}//end if
			}//end for
		}
		if(singlepoint)
		{
			pot[jj*ngpx+ii]=minpot1[n];
			edens[jj*ngpx+ii]=minedens1[n];
			hdens[jj*ngpx+ii]=minhdens1[n];
			xfield[jj*ngpx+ii]=minefx1[n];
			yfield[jj*ngpx+ii]=minefy1[n];
			Efn[jj*ngpx+ii]=minEfn1[n];
			Efp[jj*ngpx+ii]=minEfp1[n];
			Ec[jj*ngpx+ii]=minEc1[n];
			Ev[jj*ngpx+ii]=minEv1[n];
			Ei[jj*ngpx+ii]=(Ec[jj*ngpx+ii]+Ev[jj*ngpx+ii])/2.0;
		}
		//验证是否在线上
		else if(!singlepoint)
		{
			//验证该点是否在updown线上
			for(i=0;i<minn1;i++)
			{
				if((fabs(bpy1[i]-gridyy[jj])<1e-10)&&(bpx1up[i]<gridxx[ii])&&(bpx1down[i]>gridxx[ii]))
				{
					bordernumber++;
					borderpoint=true;
					borderpointupdown=true;
					break;
				}//end if
			}//end for

			if(borderpoint&&borderpointupdown)
			{
				//找出是在这条线上哪两个点围住它
				for(j=0;j<m;j++)
				{
					if((fabs(mingridy1[j]-gridyy[jj])<1e-10)&&(fabs(mingridy1[j-1]-gridyy[jj])<1e-10)&&(mingridx1[j]>gridxx[ii])&&(mingridx1[j-1]<gridxx[ii]))
					{
						p=j;//p指明了围住这个点的这两个点在从左上角排序方法中的位置
						break;
					}//end if
				}//end for

				pot[jj*ngpx+ii]   =minpot1[p-1]+(minpot1[p]-minpot1[p-1])
								  *(gridxx[ii]-mingridx1[p-1])/(mingridx1[p]-mingridx1[p-1]);
				edens[jj*ngpx+ii] =minedens1[p-1]+(minedens1[p]-minedens1[p-1])
								  *(gridxx[ii]-mingridx1[p-1])/(mingridx1[p]-mingridx1[p-1]);
				hdens[jj*ngpx+ii] =minhdens1[p-1]+(minhdens1[p]-minhdens1[p-1])
								  *(gridxx[ii]-mingridx1[p-1])/(mingridx1[p]-mingridx1[p-1]);
				xfield[jj*ngpx+ii]=minefx1[p-1]+(minefx1[p]-minefx1[p-1])
								  *(gridxx[ii]-mingridx1[p-1])/(mingridx1[p]-mingridx1[p-1]);
				yfield[jj*ngpx+ii]=minefy1[p-1]+(minefy1[p]-minefy1[p-1])
								  *(gridxx[ii]-mingridx1[p-1])/(mingridx1[p]-mingridx1[p-1]);
				Efn[jj*ngpx+ii]   =minEfn1[p-1]+(minEfn1[p]-minEfn1[p-1])
								  *(gridxx[ii]-mingridx1[p-1])/(mingridx1[p]-mingridx1[p-1]);
				Efp[jj*ngpx+ii]   =minEfp1[p-1]+(minEfp1[p]-minEfp1[p-1])
								  *(gridxx[ii]-mingridx1[p-1])/(mingridx1[p]-mingridx1[p-1]);
				Ec[jj*ngpx+ii]    =minEc1[p-1]+(minEc1[p]-minEc1[p-1])
								  *(gridxx[ii]-mingridx1[p-1])/(mingridx1[p]-mingridx1[p-1]);
				Ev[jj*ngpx+ii]    =minEv1[p-1]+(minEv1[p]-minEv1[p-1])
								  *(gridxx[ii]-mingridx1[p-1])/(mingridx1[p]-mingridx1[p-1]);
				Ei[jj*ngpx+ii]    =(Ec[jj*ngpx+ii]+Ev[jj*ngpx+ii])/2.0;
			}//end if

			if(!borderpointupdown)
			{
			//验证是否在leftright线上
				for(i=0;i<minn3;i++)
				{
					if((fabs(bpx3[i]-gridxx[ii])<1e-10)&&(bpy3left[i]<gridyy[jj])&&(bpy3right[i]>gridyy[jj]))
					{
						bordernumber++;
						borderpoint=true;
						borderpointleftright=true;
						break;
					}//end if
				}//end for
			}

			if(borderpoint&&borderpointleftright)
			{
				//找出是在这条线上哪两个点围住它
				for(j=0;j<m;j++)
				{
					if((fabs(mingridx3[j]-gridxx[ii])<1e-10)&&(fabs(mingridx3[j-1]-gridxx[ii])<1e-10)&&(mingridy3[j]>gridyy[jj])&&(mingridy3[j-1]<gridyy[jj]))
					{
						p=j;//p指明了围住这个点的这两个点在从左下角排序方法中的位置
						break;
					}//end if
				}//end for
				pot[jj*ngpx+ii]   =minpot3[p-1]+(minpot3[p]-minpot3[p-1])
								  *(gridyy[jj]-mingridy3[p-1])/(mingridy3[p]-mingridy3[p-1]);
				edens[jj*ngpx+ii] =minedens3[p-1]+(minedens3[p]-minedens3[p-1])
								  *(gridyy[jj]-mingridy3[p-1])/(mingridy3[p]-mingridy3[p-1]);
				hdens[jj*ngpx+ii] =minhdens3[p-1]+(minhdens3[p]-minhdens3[p-1])
								  *(gridyy[jj]-mingridy3[p-1])/(mingridy3[p]-mingridy3[p-1]);
				xfield[jj*ngpx+ii]=minefx3[p-1]+(minefx3[p]-minefx3[p-1])
								  *(gridyy[jj]-mingridy3[p-1])/(mingridy3[p]-mingridy3[p-1]);
				yfield[jj*ngpx+ii]=minefy3[p-1]+(minefy3[p]-minefy3[p-1])
								  *(gridyy[jj]-mingridy3[p-1])/(mingridy3[p]-mingridy3[p-1]);
				Efn[jj*ngpx+ii]   =minEfn3[p-1]+(minEfn3[p]-minEfn3[p-1])
								  *(gridyy[jj]-mingridy3[p-1])/(mingridy3[p]-mingridy3[p-1]);
				Efp[jj*ngpx+ii]   =minEfp3[p-1]+(minEfp3[p]-minEfp3[p-1])
								  *(gridyy[jj]-mingridy3[p-1])/(mingridy3[p]-mingridy3[p-1]);
				Ec[jj*ngpx+ii]   =minEc3[p-1]+(minEc3[p]-minEc3[p-1])
								  *(gridyy[jj]-mingridy3[p-1])/(mingridy3[p]-mingridy3[p-1]);
				Ev[jj*ngpx+ii]   =minEv3[p-1]+(minEv3[p]-minEv3[p-1])
								  *(gridyy[jj]-mingridy3[p-1])/(mingridy3[p]-mingridy3[p-1]);
				Ei[jj*ngpx+ii]   =(Ec[jj*ngpx+ii]+Ev[jj*ngpx+ii])/2.0;
			}//end if

		}//end if

		if(!singlepoint&&!borderpoint)
		{
			rectanglenumber++;
			//下面的四个大循环目的是为了得到包围该点的最小矩形的四个端点坐标;
			for(i=0;i<minn1-1;i++)
			{
				if((bpy1[i+1]>gridyy[jj])&&(bpy1[i]<gridyy[jj]))
				{
					minflag1=true;
					n=i;
					for(j=n;j>=0;j--)
					{
						if((bpx1up[j]<gridxx[ii])&&(bpx1down[j]>gridxx[ii]))
						{
							yy[0]=bpy1[j];
							yy[1]=bpy1[j];
							break;
						}
					}
					break;
				}
			}//end for

			for(i=0;i<minn2-1;i++)
			{
				if((bpy2[i+1]<gridyy[jj])&&(bpy2[i]>gridyy[jj]))
				{
					minflag2=true;
					n=i;
					for(j=n;j>=0;j--)
					{
						if((bpx2down[j]>gridxx[ii])&&(bpx2up[j]<gridxx[ii]))
						{
							yy[2]=bpy2[j];
							yy[3]=bpy2[j];
							break;
						}
					}
					break;
				}
			}//end for

			for(i=0;i<minn3-1;i++)
			{
				if((bpx3[i+1]<gridxx[ii])&&(bpx3[i]>gridxx[ii]))
				{
					minflag3=true;
					n=i;
					for(j=n;j>=0;j--)
					{
						if((bpy3left[j]<gridyy[jj])&&(bpy3right[j]>gridyy[jj]))
						{
							xx[1]=bpx3[j];
							xx[2]=bpx3[j];
							break;
						}
					}
					break;
				}
			}//end for

			for(i=0;i<minn4-1;i++)
			{
				if((bpx4[i+1]>gridxx[ii])&&(bpx4[i]<gridxx[ii]))
				{
					minflag4=true;
					n=i;
					for(j=n;j>=0;j--)
					{
						if((bpy4right[j]>gridyy[jj])&&(bpy4left[j]<gridyy[jj]))
						{
							xx[3]=bpx4[j];
							xx[0]=bpx4[j];
							break;
						}
					}
					break;
				}
			}//end for


			//得到四个点处的电势值
			if(minflag1&&minflag2&&minflag3&&minflag4)
			{
			for(i=0;i<m;i++)
			{
				if((fabs(isegridx[i]-xx[0])<1e-10)&&(fabs(isegridy[i]-yy[0])<1e-10))
				{
					temppot[0]=isepot[i];
					tempedens[0]=iseedens[i];
					temphdens[0]=isehdens[i];
					tempefx[0]=iseefx[i];
					tempefy[0]=iseefy[i];
					tempEfn[0]=iseEfn[i];
					tempEfp[0]=iseEfp[i];
					tempEc[0]=iseEc[i];
					tempEv[0]=iseEv[i];
					continue;
				}
				if((fabs(isegridx[i]-xx[1])<1e-10)&&(fabs(isegridy[i]-yy[1])<1e-10))
				{
					temppot[1]=isepot[i];
					tempedens[1]=iseedens[i];
					temphdens[1]=isehdens[i];
					tempefx[1]=iseefx[i];
					tempefy[1]=iseefy[i];
					tempEfn[1]=iseEfn[i];
					tempEfp[1]=iseEfp[i];
					tempEc[1]=iseEc[i];
					tempEv[1]=iseEv[i];
					continue;
				}
				if((fabs(isegridx[i]-xx[2])<1e-10)&&(fabs(isegridy[i]-yy[2])<1e-10))
				{
					temppot[2]=isepot[i];
					tempedens[2]=iseedens[i];
					temphdens[2]=isehdens[i];
					tempefx[2]=iseefx[i];
					tempefy[2]=iseefy[i];
					tempEfn[2]=iseEfn[i];
					tempEfp[2]=iseEfp[i];
					tempEc[2]=iseEc[i];
					tempEv[2]=iseEv[i];
					continue;
				}
				if((fabs(isegridx[i]-xx[3])<1e-10)&&(fabs(isegridy[i]-yy[3])<1e-10))
				{
					temppot[3]=isepot[i];
					tempedens[3]=iseedens[i];
					temphdens[3]=isehdens[i];
					tempefx[3]=iseefx[i];
					tempefy[3]=iseefy[i];
					tempEfn[3]=iseEfn[i];
					tempEfp[3]=iseEfp[i];
					tempEc[3]=iseEc[i];
					tempEv[3]=iseEv[i];
					continue;
				}
			}//end for
			potup=temppot[0]+(temppot[3]-temppot[0])*(gridyy[jj]-yy[0])/(yy[3]-yy[0]);
			potdown=temppot[1]+(temppot[2]-temppot[1])*(gridyy[jj]-yy[1])/(yy[2]-yy[1]);
			potleft=temppot[0]+(temppot[1]-temppot[0])*(gridxx[ii]-xx[0])/(xx[1]-xx[0]);
			potright=temppot[3]+(temppot[2]-temppot[3])*(gridxx[ii]-xx[3])/(xx[2]-xx[3]);
			pot[jj*ngpx+ii]=(potleft+(potright-potleft)*(gridyy[jj]-yy[0])/(yy[3]-yy[0])
							+potup+(potdown-potup)*(gridxx[ii]-xx[0])/(xx[1]-xx[0]))/2.0;

			edensup=tempedens[0]+(tempedens[3]-tempedens[0])*(gridyy[jj]-yy[0])/(yy[3]-yy[0]);
			edensdown=tempedens[1]+(tempedens[2]-tempedens[1])*(gridyy[jj]-yy[1])/(yy[2]-yy[1]);
			edensleft=tempedens[0]+(tempedens[1]-tempedens[0])*(gridxx[ii]-xx[0])/(xx[1]-xx[0]);
			edensright=tempedens[3]+(tempedens[2]-tempedens[3])*(gridxx[ii]-xx[3])/(xx[2]-xx[3]);
			edens[jj*ngpx+ii]=(edensleft+(edensright-edensleft)*(gridyy[jj]-yy[0])/(yy[3]-yy[0])
							  +edensup+(edensdown-edensup)*(gridxx[ii]-xx[0])/(xx[1]-xx[0]))/2.0;

			hdensup=temphdens[0]+(temphdens[3]-temphdens[0])*(gridyy[jj]-yy[0])/(yy[3]-yy[0]);
			hdensdown=temphdens[1]+(temphdens[2]-temphdens[1])*(gridyy[jj]-yy[1])/(yy[2]-yy[1]);
			hdensleft=temphdens[0]+(temphdens[1]-temphdens[0])*(gridxx[ii]-xx[0])/(xx[1]-xx[0]);
			hdensright=temphdens[3]+(temphdens[2]-temphdens[3])*(gridxx[ii]-xx[3])/(xx[2]-xx[3]);
			hdens[jj*ngpx+ii]=(hdensleft+(hdensright-hdensleft)*(gridyy[jj]-yy[0])/(yy[3]-yy[0])
							  +hdensup+(hdensdown-hdensup)*(gridxx[ii]-xx[0])/(xx[1]-xx[0]))/2.0;

			efxup=tempefx[0]+(tempefx[3]-tempefx[0])*(gridyy[jj]-yy[0])/(yy[3]-yy[0]);
			efxdown=tempefx[1]+(tempefx[2]-tempefx[1])*(gridyy[jj]-yy[1])/(yy[2]-yy[1]);
			efxleft=tempefx[0]+(tempefx[1]-tempefx[0])*(gridxx[ii]-xx[0])/(xx[1]-xx[0]);
			efxright=tempefx[3]+(tempefx[2]-tempefx[3])*(gridxx[ii]-xx[3])/(xx[2]-xx[3]);
			xfield[jj*ngpx+ii]=(efxleft+(efxright-efxleft)*(gridyy[jj]-yy[0])/(yy[3]-yy[0])
							+efxup+(efxdown-efxup)*(gridxx[ii]-xx[0])/(xx[1]-xx[0]))/2.0;

			efyup=tempefy[0]+(tempefy[3]-tempefy[0])*(gridyy[jj]-yy[0])/(yy[3]-yy[0]);
			efydown=tempefy[1]+(tempefy[2]-tempefy[1])*(gridyy[jj]-yy[1])/(yy[2]-yy[1]);
			efyleft=tempefy[0]+(tempefy[1]-tempefy[0])*(gridxx[ii]-xx[0])/(xx[1]-xx[0]);
			efyright=tempefy[3]+(tempefy[2]-tempefy[3])*(gridxx[ii]-xx[3])/(xx[2]-xx[3]);
			yfield[jj*ngpx+ii]=(efyleft+(efyright-efyleft)*(gridyy[jj]-yy[0])/(yy[3]-yy[0])
							+efyup+(efydown-efyup)*(gridxx[ii]-xx[0])/(xx[1]-xx[0]))/2.0;

			Efnup=tempEfn[0]+(tempEfn[3]-tempEfn[0])*(gridyy[jj]-yy[0])/(yy[3]-yy[0]);
			Efndown=tempEfn[1]+(tempEfn[2]-tempEfn[1])*(gridyy[jj]-yy[1])/(yy[2]-yy[1]);
			Efnleft=tempEfn[0]+(tempEfn[1]-tempEfn[0])*(gridxx[ii]-xx[0])/(xx[1]-xx[0]);
			Efnright=tempEfn[3]+(tempEfn[2]-tempEfn[3])*(gridxx[ii]-xx[3])/(xx[2]-xx[3]);
			Efn[jj*ngpx+ii]=(Efnleft+(Efnright-Efnleft)*(gridyy[jj]-yy[0])/(yy[3]-yy[0])
							+Efnup+(Efndown-Efnup)*(gridxx[ii]-xx[0])/(xx[1]-xx[0]))/2.0;

			Efpup=tempEfp[0]+(tempEfp[3]-tempEfp[0])*(gridyy[jj]-yy[0])/(yy[3]-yy[0]);
			Efpdown=tempEfp[1]+(tempEfp[2]-tempEfp[1])*(gridyy[jj]-yy[1])/(yy[2]-yy[1]);
			Efpleft=tempEfp[0]+(tempEfp[1]-tempEfp[0])*(gridxx[ii]-xx[0])/(xx[1]-xx[0]);
			Efpright=tempEfp[3]+(tempEfp[2]-tempEfp[3])*(gridxx[ii]-xx[3])/(xx[2]-xx[3]);
			Efp[jj*ngpx+ii]=(Efpleft+(Efpright-Efpleft)*(gridyy[jj]-yy[0])/(yy[3]-yy[0])
							+Efpup+(Efpdown-Efpup)*(gridxx[ii]-xx[0])/(xx[1]-xx[0]))/2.0;

			Ecup=tempEc[0]+(tempEc[3]-tempEc[0])*(gridyy[jj]-yy[0])/(yy[3]-yy[0]);
			Ecdown=tempEc[1]+(tempEc[2]-tempEc[1])*(gridyy[jj]-yy[1])/(yy[2]-yy[1]);
			Ecleft=tempEc[0]+(tempEc[1]-tempEc[0])*(gridxx[ii]-xx[0])/(xx[1]-xx[0]);
			Ecright=tempEc[3]+(tempEc[2]-tempEc[3])*(gridxx[ii]-xx[3])/(xx[2]-xx[3]);
			Ec[jj*ngpx+ii]=(Ecleft+(Ecright-Ecleft)*(gridyy[jj]-yy[0])/(yy[3]-yy[0])
							+Ecup+(Ecdown-Ecup)*(gridxx[ii]-xx[0])/(xx[1]-xx[0]))/2.0;

			Evup=tempEv[0]+(tempEv[3]-tempEv[0])*(gridyy[jj]-yy[0])/(yy[3]-yy[0]);
			Evdown=tempEv[1]+(tempEv[2]-tempEv[1])*(gridyy[jj]-yy[1])/(yy[2]-yy[1]);
			Evleft=tempEv[0]+(tempEv[1]-tempEv[0])*(gridxx[ii]-xx[0])/(xx[1]-xx[0]);
			Evright=tempEv[3]+(tempEv[2]-tempEv[3])*(gridxx[ii]-xx[3])/(xx[2]-xx[3]);
			Ev[jj*ngpx+ii]=(Evleft+(Evright-Evleft)*(gridyy[jj]-yy[0])/(yy[3]-yy[0])
							+Evup+(Evdown-Evup)*(gridxx[ii]-xx[0])/(xx[1]-xx[0]))/2.0;

			Ei[jj*ngpx+ii]=(Ec[jj*ngpx+ii]+Ev[jj*ngpx+ii])/2.0;

			}
			else if(!minflag1||!minflag2||!minflag3||!minflag4)
			{
				xfield[jj*ngpx+ii]=0;
				yfield[jj*ngpx+ii]=0;
			}

			if(yfield[jj*ngpx+ii]*field0<-1e10)
				nn3++;

		}//end if
		
		if(mat[ii+jj*ngpx]==VOID)
		{
			pot[jj*ngpx+ii]=0;
			edens[jj*ngpx+ii]=0;
			hdens[jj*ngpx+ii]=0;
			xfield[jj*ngpx+ii]=0;
			yfield[jj*ngpx+ii]=0;
			Efn[jj*ngpx+ii]=0;
			Efp[jj*ngpx+ii]=0;
			Ec[jj*ngpx+ii]=0;
			Ev[jj*ngpx+ii]=0;
			Ei[jj*ngpx+ii]=0;
		}
		if(mat[ii+jj*ngpx]==OXIDE)
		{
			Ec[jj*ngpx+ii]=0;
			Ev[jj*ngpx+ii]=0;
		}
	}//end for for
}
// Method fro getting data from ise is four points!
else
{
	int number0;
	double temppot[MNGP],tempedens[MNGP],temphdens[MNGP],tempxfield[MNGP],tempyfield[MNGP];
	double tempEfn[MNGP],tempEfp[MNGP],tempEc[MNGP],tempEv[MNGP];
	double tempisegridx[MNGP],tempisegridy[MNGP];
	double tempdistance[MNGP];
	double tdistance;
	double tpot,tedens,thdens,txfield,tyfield,tEfn,tEfp,tEc,tEv,tisegridx,tisegridy;
	double ttdistance;

	for(i=0;i<MNGP;i++)
	{
		temppot[i]=0;
		tempedens[i]=0;
		temphdens[i]=0;
		tempxfield[i]=0;
		tempyfield[i]=0;
		tempEfn[i]=0;
		tempEfp[i]=0;
		tempEc[i]=0;
		tempEv[i]=0;
		tempdistance[i]=0;
	}
	tdistance=0;
	ttdistance=0;

	for(ii=0;ii<ngpx-1;ii++)
	{
		gridxx[ii]=(gridx[ii]+gridx[ii+1])/2.0;
	}
	for(jj=0;jj<ngpy-1;jj++)
	{
		gridyy[jj]=(gridy[jj]+gridy[jj+1])/2.0;
	}
	m=isegridnumber;

for(j=0;j<ngpy-1;j++)
for(i=0;i<ngpx-1;i++)
{
	cout<<i<<"      "<<j<<"      ngpx= "<<ngpx<<"      ngpy= "<<ngpy<<endl;
	if(mat[i+j*ngpx]==VOID)
	{
		pot[i+j*ngpx]=0;
		edens[i+j*ngpx]=0;
		hdens[i+j*ngpx]=0;
		xfield[i+j*ngpx]=0;
		yfield[i+j*ngpx]=0;
		Efn[i+j*ngpx]=0;
		Efp[i+j*ngpx]=0;
		Ei[i+j*ngpx]=0;
		Ec[i+j*ngpx]=0;
		Ev[i+j*ngpx]=0;
	}
	else
	{
		for(k=0;k<m;k++)
		{
			tempdistance[k]=(gridxx[i]-isegridx[k])*(gridxx[i]-isegridx[k])+
						    (gridyy[j]-isegridy[k])*(gridyy[j]-isegridy[k]);
			tempisegridx[k]=isegridx[k];
			tempisegridy[k]=isegridy[k];
			temppot[k]=isepot[k];
			tempedens[k]=iseedens[k];
			temphdens[k]=isehdens[k];
			tempxfield[k]=iseefx[k];
			tempyfield[k]=iseefy[k];
			tempEfn[k]=iseEfn[k];
			tempEfp[k]=iseEfp[k];
			tempEc[k]=iseEc[k];
			tempEv[k]=iseEv[k];
		}
		for(ii=0;ii<m;ii++)
		{
			k=ii;
			for(jj=ii+1;jj<m;jj++)
			{
				if(tempdistance[jj]<tempdistance[k])	
				{
					k=jj;
				}
			}
			if(k!=ii)
			{
				tdistance=tempdistance[ii];tempdistance[ii]=tempdistance[k];tempdistance[k]=tdistance;
				tisegridx=tempisegridx[ii];tempisegridx[ii]=tempisegridx[k];tempisegridx[k]=tisegridx;
				tisegridy=tempisegridy[ii];tempisegridy[ii]=tempisegridy[k];tempisegridy[k]=tisegridy;
				tpot=temppot[ii];temppot[ii]=temppot[k];temppot[k]=tpot;
				tedens=tempedens[ii];tempedens[ii]=tempedens[k];tempedens[k]=tedens;
				thdens=temphdens[ii];temphdens[ii]=temphdens[k];temphdens[k]=thdens;
				txfield=tempxfield[ii];tempxfield[ii]=tempxfield[k];tempxfield[k]=txfield;
				tyfield=tempyfield[ii];tempyfield[ii]=tempyfield[k];tempyfield[k]=tyfield;
				tEfn=tempEfn[ii];tempEfn[ii]=tempEfn[k];tempEfn[k]=tEfn;
				tEfp=tempEfp[ii];tempEfp[ii]=tempEfp[k];tempEfp[k]=tEfp;
				tEc=tempEc[ii];tempEc[ii]=tempEc[k];tempEc[k]=tEc;
				tEv=tempEv[ii];tempEv[ii]=tempEv[k];tempEv[k]=tEv;
			}
		}
		number0=0;
		ttdistance=0;
		tpot=0;
		tedens=0;
		thdens=0;
		txfield=0;
		tyfield=0;
		tEfn=0;
		tEfp=0;
		tEc=0;
		tEv=0;
		for(p=0;p<4;p++)
		{
			if(tempdistance[p]==0.0)
			{
				number0++;
			}
		}
		if(number0!=0)
		{
			for(p=0;p<number0;p++)
			{
				tpot+=temppot[p];
				tedens+=tempedens[p];
				thdens+=temphdens[p];
				txfield+=tempxfield[p];
				tyfield+=tempyfield[p];
				tEfn+=tempEfn[p];
				tEfp+=tempEfp[p];
				tEc+=tempEc[p];
				tEv+=tempEv[p];
			}
			pot[i+j*ngpx]=tpot/number0;
			edens[i+j*ngpx]=tedens/number0;
			hdens[i+j*ngpx]=thdens/number0;
			xfield[i+j*ngpx]=txfield/number0;
			yfield[i+j*ngpx]=tyfield/number0;
			Efn[i+j*ngpx]=tEfn/number0;
			Efp[i+j*ngpx]=tEfp/number0;
			Ec[i+j*ngpx]=tEc/number0;
			Ev[i+j*ngpx]=tEv/number0;
			Ei[i+j*ngpx]=(Ec[i+j*ngpx]+Ev[i+j*ngpx])/2.0;
		}
		else if(number0==0)
		{
			//for(p=0;p<4;p++)
			p=0;
			l=0;
			while(p<4)
			{
				if(tempisegridx[p+l]>gridx[ibreg[numreg[i+j*ngpx]]]&&tempisegridx[p+l]<gridx[iereg[numreg[i+j*ngpx]]]&&
				   tempisegridy[p+l]>gridy[jbreg[numreg[i+j*ngpx]]]&&tempisegridy[p+l]<gridy[jereg[numreg[i+j*ngpx]]])
				{
					tpot+=temppot[p+l]/tempdistance[p+l];
					tedens+=tempedens[p+l]/tempdistance[p+l];
					thdens+=temphdens[p+l]/tempdistance[p+l];
					txfield+=tempxfield[p+l]/tempdistance[p+l];
					tyfield+=tempyfield[p+l]/tempdistance[p+l];
					tEfn+=tempEfn[p+l]/tempdistance[p+l];
					tEfp+=tempEfp[p+l]/tempdistance[p+l];
					tEc+=tempEc[p+l]/tempdistance[p+l];
					tEv+=tempEv[p+l]/tempdistance[p+l];
					ttdistance+=1.0/tempdistance[p+l];
					p++;
				}
				else
				{
					l++;
				}
			}
			pot[i+j*ngpx]=tpot/ttdistance;
			edens[i+j*ngpx]=tedens/ttdistance;
			hdens[i+j*ngpx]=thdens/ttdistance;
			xfield[i+j*ngpx]=txfield/ttdistance;
			yfield[i+j*ngpx]=tyfield/ttdistance;
			Efn[i+j*ngpx]=tEfn/ttdistance;
			Efp[i+j*ngpx]=tEfp/ttdistance;
			Ec[i+j*ngpx]=tEc/ttdistance;
			Ev[i+j*ngpx]=tEv/ttdistance;
			Ei[i+j*ngpx]=(Ec[i+j*ngpx]+Ev[i+j*ngpx])/2.0;
		}
	}
}
}
	cout<<"------OK!"<<endl<<endl;
/*
	//Calc Nc Nv
	for(j=0;j<ngpy-1;j++)
	for(i=0;i<ngpx-1;i++)
	{
		dc[i+j*ngpx]=(Efn[i+j*ngpx]-Ec[j*ngpx+i])/0.0258;
		dv[i+j*ngpx]=(Ev[i+j*ngpx]-Efp[j*ngpx+i])/0.0258;

		if(dc[i+j*ngpx]<fermi[0][0])	
		{
			dc[i+j*ngpx]=fermi[0][0];
			Nc[i+j*ngpx]=fermi[1][0];
		}
		else if(dc[i+j*ngpx]>fermi[0][nfpoint-1])
		{
			dc[i+j*ngpx]=fermi[0][nfpoint-1];
			Nc[i+j*ngpx]=fermi[1][nfpoint-1];
		}
		else
		{
			for(p=0;p<nfpoint-1;p++)
			{
				if(dc[i+j*ngpx]>=fermi[0][p]&&dc[i+j*ngpx]<=fermi[0][p+1])
				{
					m=p;
					Nc[i+j*ngpx]=exp(log(fermi[1][m])+(log(fermi[1][m+1])-log(fermi[1][m]))
						         *(dc[i+j*ngpx]-fermi[0][m])/(fermi[0][m+1]-fermi[0][m]));
					Nc[i+j*ngpx]=edens[i+j*ngpx]/Nc[i+j*ngpx];
					break;
				}
			}
		}

		if(dv[i+j*ngpx]<fermi[0][0])	
		{
			dv[i+j*ngpx]=fermi[0][0];
			Nv[i+j*ngpx]=fermi[1][0];
		}
		else if(dv[i+j*ngpx]>fermi[0][nfpoint-1])
		{
			dv[i+j*ngpx]=fermi[0][nfpoint-1];
			Nv[i+j*ngpx]=fermi[1][nfpoint-1];
		}
		else
		{
			for(p=0;p<nfpoint-1;p++)
			{
				if(dv[i+j*ngpx]>=fermi[0][p]&&dv[i+j*ngpx]<=fermi[0][p+1])
				{
					n=p;
					Nv[i+j*ngpx]=exp(log(fermi[1][n])+(log(fermi[1][n+1])-log(fermi[1][n]))
							     *(dv[i+j*ngpx]-fermi[0][n])/(fermi[0][n+1]-fermi[0][n]));
					Nv[i+j*ngpx]=hdens[i+j*ngpx]/Nv[i+j*ngpx];
					break;
				}
			}
		}

	}
	//Calc Nc[MNGP],Nv[MNGP] end
*/
	ofstream ftpout;
	ftpout.open("output/isedataout.txt");
	for(j=0;j<ngpy-1;j++)
	for(i=0;i<ngpx-1;i++)
	{
		ftpout<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9<<"  "<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9<<"  "<<pot[j*ngpx+i]*pot0
			  <<"  "<<edens[j*ngpx+i]<<"  "<<hdens[j*ngpx+i]<<"  "<<xfield[j*ngpx+i]*field0/100.0
			  <<"  "<<yfield[j*ngpx+i]*field0/100.0<<"  "<<Efn[j*ngpx+i]<<"  "<<Efp[j*ngpx+i]
			  <<"  "<<Ec[j*ngpx+i]<<"  "<<Ev[j*ngpx+i]<<"  "<<Ei[j*ngpx+i]
			  <<"  "<<Nc[j*ngpx+i]<<"  "<<Nv[j*ngpx+i]<<endl;
	}
	ftpout.close();

	ofstream pxftp("output/ise_potcutx.txt");//		V		1
	ofstream exftp("output/ise_edenscutx.txt");//	cm-3	2
	ofstream hxftp("output/ise_hdenscutx.txt");//	cm-3	3
	ofstream exxftp("output/ise_xfieldcutx.txt");//	V/cm	4
	ofstream eyxftp("output/ise_yfieldcutx.txt");//	V/cm	5
	ofstream Efnxftp("output/ise_Efncutx.txt");//	V/cm	6
	ofstream Efpxftp("output/ise_Efpcutx.txt");//	V/cm	7
	ofstream Ecxftp("output/ise_Eccutx.txt");//	V/cm		8
	ofstream Evxftp("output/ise_Evcutx.txt");//	V/cm		9
	ofstream Eixftp("output/ise_Eicutx.txt");//	V/cm		10

// output xcut
	for(j=0;j<ngpy-1;j++)
	{
		pxftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;	//1
		exftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;	//2
		hxftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;	//3
		exxftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;	//4
		eyxftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;	//5
		Efnxftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;	//6
		Efpxftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;	//7
		Ecxftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;	//8
		Evxftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;	//9
		Eixftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;	//10
		for(i=0;i<xcutnumber;i++)
		{
			pxftp<<setw(15)<<pot[j*ngpx+ipcutx[i]]*pot0;
			exftp<<setw(15)<<edens[j*ngpx+ipcutx[i]];
			hxftp<<setw(15)<<hdens[j*ngpx+ipcutx[i]];
			exxftp<<setw(15)<<xfield[j*ngpx+ipcutx[i]]*field0/100.0;
			eyxftp<<setw(15)<<yfield[j*ngpx+ipcutx[i]]*field0/100.0;
			Efnxftp<<setw(15)<<Efn[j*ngpx+ipcutx[i]];
			Efpxftp<<setw(15)<<Efp[j*ngpx+ipcutx[i]];
			Ecxftp<<setw(15)<<Ec[j*ngpx+ipcutx[i]];
			Evxftp<<setw(15)<<Ev[j*ngpx+ipcutx[i]];
			Eixftp<<setw(15)<<Ei[j*ngpx+ipcutx[i]];
		}
		pxftp<<endl;
		exftp<<endl;
		hxftp<<endl;
		exxftp<<endl;
		eyxftp<<endl;
		Efnxftp<<endl;
		Efpxftp<<endl;
		Ecxftp<<endl;
		Evxftp<<endl;
		Eixftp<<endl;
	}	

// output ycut

	ofstream pyftp("output/ise_potcuty.txt");//		V		1
	ofstream eyftp("output/ise_edenscuty.txt");//	cm-3	2
	ofstream hyftp("output/ise_hdenscuty.txt");//	cm-3	3
	ofstream exyftp("output/ise_xfieldcuty.txt");//	V/cm	4
	ofstream eyyftp("output/ise_yfieldcuty.txt");//	V/cm	5
	ofstream Efnyftp("output/ise_Efncuty.txt");//	V/cm	6
	ofstream Efpyftp("output/ise_Efpcuty.txt");//	V/cm	7
	ofstream Ecyftp("output/ise_Eccuty.txt");//	V/cm		8
	ofstream Evyftp("output/ise_Evcuty.txt");//	V/cm		9
	ofstream Eiyftp("output/ise_Eicuty.txt");//	V/cm		10

	for(i=0;i<ngpx-1;i++)
	{
		pyftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;	//1
		eyftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;	//2
		hyftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;	//3
		exyftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;	//4
		eyyftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;	//5
		Efnyftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;	//6
		Efpyftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;	//7
		Ecyftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;	//8
		Evyftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;	//9
		Eiyftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;	//10
		for(j=0;j<ycutnumber;j++)
		{
			pyftp<<setw(15)<<pot[ipcuty[j]*ngpx+i]*pot0;
			eyftp<<setw(15)<<edens[ipcuty[j]*ngpx+i];
			hyftp<<setw(15)<<hdens[ipcuty[j]*ngpx+i];
			exyftp<<setw(15)<<xfield[ipcuty[j]*ngpx+i]*field0/100.0;
			eyyftp<<setw(15)<<yfield[ipcuty[j]*ngpx+i]*field0/100.0;
			Efnyftp<<setw(15)<<Efn[ipcuty[j]*ngpx+i];
			Efpyftp<<setw(15)<<Efp[ipcuty[j]*ngpx+i];
			Ecyftp<<setw(15)<<Ec[ipcuty[j]*ngpx+i];
			Evyftp<<setw(15)<<Ev[ipcuty[j]*ngpx+i];
			Eiyftp<<setw(15)<<Ei[ipcuty[j]*ngpx+i];
		}
		pyftp<<endl;
		eyftp<<endl;
		hyftp<<endl;
		exyftp<<endl;
		eyyftp<<endl;
		Efnyftp<<endl;
		Efpyftp<<endl;
		Ecyftp<<endl;
		Evyftp<<endl;
		Eiyftp<<endl;
	}
	return;
}

void DevSimulator::CALCFERMI(void)
{
	int i,j,m,n,p;
	double eedens[MNGP],hhdens[MNGP];
	for(i=0;i<MNGP;i++)	eedens[i]=0;
	for(i=0;i<MNGP;i++)	hhdens[i]=0;

	for(j=0;j<ngpy-1;j++)
	for(i=0;i<ngpx-1;i++)
	{
		if(j==80)
		{
			j=j;
		}
		eedens[i+j*ngpx]=edens[i+j*ngpx]+quadndenss[1][i+j*ngpx]*conc0/1e6+quadndenss[2][i+j*ngpx]*conc0/1e6;
		if(Nc[i+j*ngpx]!=0)	dc[i+j*ngpx]=eedens[i+j*ngpx]/Nc[i+j*ngpx];
		else dc[i+j*ngpx]=0;

		if(dc[i+j*ngpx]<fermi[1][0])	
		{
			dc[i+j*ngpx]=fermi[1][0];
			Efn[i+j*ngpx]=fermi[0][0]+Ec[i+j*ngpx];
		}
		else if(dc[i+j*ngpx]>fermi[1][nfpoint-1])
		{
			dc[i+j*ngpx]=fermi[1][nfpoint-1];
			Efn[i+j*ngpx]=fermi[0][nfpoint-1]+Ec[i+j*ngpx];
		}
		else
		{
			for(p=0;p<nfpoint-1;p++)
			{
				if(dc[i+j*ngpx]>=fermi[1][p]&&dc[i+j*ngpx]<=fermi[1][p+1])
				{
					m=p;
					Efn[i+j*ngpx]=Ec[i+j*ngpx]+0.0258*(fermi[0][m]+(fermi[0][m+1]-fermi[0][m])
							      *(log(dc[i+j*ngpx])-log(fermi[1][m]))/(log(fermi[1][m+1])-log(fermi[1][m])));
					break;
				}
			}
		}

		hhdens[i+j*ngpx]=hdens[i+j*ngpx]+quadpdenss[1][i+j*ngpx]*conc0/1e6+quadpdenss[2][i+j*ngpx]*conc0/1e6;
		if(Nv[i+j*ngpx]!=0)	dv[i+j*ngpx]=hhdens[i+j*ngpx]/Nv[i+j*ngpx];
		else dv[i+j*ngpx]=0;

		if(dv[i+j*ngpx]<fermi[1][0])	
		{
			dv[i+j*ngpx]=fermi[1][0];
			Efp[i+j*ngpx]=Ev[i+j*ngpx]-fermi[0][0];
		}
		else if(dv[i+j*ngpx]>fermi[1][nfpoint-1])
		{
			dv[i+j*ngpx]=fermi[1][nfpoint-1];
			Efp[i+j*ngpx]=Ev[i+j*ngpx]-fermi[0][nfpoint-1];
		}
		else
		{
			for(p=0;p<nfpoint-1;p++)
			{
				if(dv[i+j*ngpx]>=fermi[1][p]&&dv[i+j*ngpx]<=fermi[1][p+1])
				{
					n=p;
					Efp[i+j*ngpx]=Ev[i+j*ngpx]-0.0258*(fermi[0][n]+(fermi[0][n+1]-fermi[0][n])
							      *(log(dv[i+j*ngpx])-log(fermi[1][n]))/(log(fermi[1][n+1])-log(fermi[1][n])));
					break;
				}
			}
		}
	}

	ofstream ftpout;
	ftpout.open("output/self_consistent.txt");
	for(j=0;j<ngpy-1;j++)
	for(i=0;i<ngpx-1;i++)
	{
		ftpout<<(gridx[i]+gridx[i+1])/2.0*spr0*1e6<<"  "<<(gridy[j]+gridy[j+1])/2.0*spr0*1e6
			  <<"  "<<eedens[j*ngpx+i]<<"  "<<hhdens[j*ngpx+i]<<"  "<<Efn[j*ngpx+i]<<"  "<<Efp[j*ngpx+i]
			  <<"  "<<Ec[j*ngpx+i]<<"  "<<Ev[j*ngpx+i]<<"  "<<Ei[j*ngpx+i]
			  <<"  "<<Nc[j*ngpx+i]<<"  "<<Nv[j*ngpx+i]<<endl;
	}
	ftpout.close();
}

void DevSimulator::CALCRDENS(void)
{
	int i,j;
	for(i=0;i<ngpx-1;i++)
	for(j=0;j<ngpy-1;j++)
	{
		RDENS[i+j*ngpx]=rr*(quadndens[i+j*ngpx]*quadpdens[i+j*ngpx]*conc0*conc0*1e-12-1.5*1.5*1e20);
		RDENS[i+j*ngpx]=RDENS[i+j*ngpx]*dt*time0*1e6*(gridx[i+1]-gridx[i])*(gridy[j+1]-gridy[j])*spr0*spr0*spr0;
		NRDENS[0][i+j*ngpx]=RDENS[i+j*ngpx]*quadndenss[0][i+j*ngpx]/quadndens[i+j*ngpx];
		NRDENS[1][i+j*ngpx]=RDENS[i+j*ngpx]*quadndenss[1][i+j*ngpx]/quadndens[i+j*ngpx];
		NRDENS[2][i+j*ngpx]=RDENS[i+j*ngpx]*quadndenss[2][i+j*ngpx]/quadndens[i+j*ngpx];
		PRDENS[0][i+j*ngpx]=RDENS[i+j*ngpx]*quadpdenss[0][i+j*ngpx]/quadpdens[i+j*ngpx];
		PRDENS[1][i+j*ngpx]=RDENS[i+j*ngpx]*quadpdenss[1][i+j*ngpx]/quadpdens[i+j*ngpx];
		PRDENS[2][i+j*ngpx]=RDENS[i+j*ngpx]*quadpdenss[2][i+j*ngpx]/quadpdens[i+j*ngpx];
	}
	return;
}

void DevSimulator::recombination(void)
{
	int i,j,ij;
	double pc,pcold;
	int carriertype;
	int ipt;
	for(i=0;i<npar0;i++)
	{
		ij=ifield[2][i];
		pc=dfield[5][i];
		pcold=pc;
		carriertype=ifield[IVPP][i];
		if(pc>0)	ipt=0;
		else if(pc<0)	ipt=1;
		if(quadndenss[carriertype][ij]==0)	pc=pc;
		if(ipt==0)	
		{
			if(quadndenss[carriertype][ij]==0)	pc=pc;
			else	pc=pc-NRDENS[carriertype][ij]*pc/(quadndenss[carriertype][ij]/quadarea[ij]);
		}
		else if(ipt==1)
		{	
			if(quadpdenss[carriertype][ij]==0)	pc=pc;
			else	pc=pc-PRDENS[carriertype][ij]*pc/(quadpdenss[carriertype][ij]/quadarea[ij]);
		}

		if(pcold*pc<0)	
		{
			cout<<"Recombination is error!!";
			exit(0);
		}
	}
	return;
}

void DevSimulator::CALCBTBTDENS(void)
{
	int i,j;
	double nfermi,pfermi;
	double xxfield,yyfield;
	double meff,noftrap,Et;
	double mm,trapratio1,trapratio2;
	double aa,bb;
	double alpha;
	alpha=2.0;
	aa=2.9e20;
	bb=2.56e7;
	meff=0.4;
	noftrap=1e15;
	Et=0;
	ifstream ftp;
	ftp.open("m.txt");
	ftp>>tempchar>>meff>>tempchar>>alpha>>tempchar>>noftrap>>tempchar>>Et>>tempchar>>aa>>tempchar>>bb;
	mm=sqrt(meff/0.4);
	trapratio1=noftrap/(5.18e18);
	trapratio2=(1.12-Et)/1.12*sqrt((1.12-Et)/1.12);
	for(i=0;i<ngpx-1;i++)
	for(j=0;j<ngpy-1;j++)
	{
		if(mat[i+j*ngpx]==SILICON)
		{
			nfermi=1/(1+(exp((Ei[i+j*ngpx]-Efn[i+j*ngpx])/0.0258)));
			pfermi=1/(1+(exp((Ei[i+j*ngpx]-Efp[i+j*ngpx])/0.0258)));
			DE[i+j*ngpx]=fabs(nfermi-pfermi);
			yyfield=fabs(yfield[i+j*ngpx]*field0/100.0);
			xxfield=fabs(xfield[i+j*ngpx]*field0/100.0);
			RBBL0[i+j*ngpx]=aa*mm*pow(yyfield,alpha)*exp(-bb*mm/yyfield)*DE[i+j*ngpx];
			RBBT0[i+j*ngpx]=aa*mm*pow(xxfield,alpha)*exp(-bb*mm/xxfield)*DE[i+j*ngpx];
			if((gridx[i]+gridx[i+1])/2.0*spr0*1e6>=ifminx&&(gridx[i]+gridx[i+1])/2.0*spr0*1e6<=ifmaxx
			 &&(gridy[j]+gridy[j+1])/2.0*spr0*1e6>=ifminy&&(gridy[j]+gridy[j+1])/2.0*spr0*1e6<=ifmaxy)
			{
				RBBL1[i+j*ngpx]=aa*mm*pow(yyfield,alpha)*exp(-bb*trapratio2*mm/yyfield)*DE[i+j*ngpx]*trapratio1;
				RBBT1[i+j*ngpx]=aa*mm*pow(xxfield,alpha)*exp(-bb*trapratio2*mm/xxfield)*DE[i+j*ngpx]*trapratio1;
			}
			else
			{
				RBBL1[i+j*ngpx]=0;
				RBBT1[i+j*ngpx]=0;
			}
			RBBL[i+j*ngpx]=RBBL0[i+j*ngpx]+RBBL1[i+j*ngpx];
			RBBT[i+j*ngpx]=RBBT0[i+j*ngpx]+RBBT1[i+j*ngpx];
		}
		else
		{
			RBBL0[i+j*ngpx]=0;
			RBBL1[i+j*ngpx]=0;
			RBBL[i+j*ngpx]=RBBL0[i+j*ngpx]+RBBL1[i+j*ngpx];
			RBBT0[i+j*ngpx]=0;
			RBBT1[i+j*ngpx]=0;
			RBBT[i+j*ngpx]=RBBT0[i+j*ngpx]+RBBT1[i+j*ngpx];
		}
		RBBL0DENS[i+j*ngpx]=RBBL0[i+j*ngpx]*dt*time0*1e6*(gridx[i+1]-gridx[i])*(gridy[j+1]-gridy[j])*spr0*spr0*spr0;
		RBBT0DENS[i+j*ngpx]=RBBT0[i+j*ngpx]*dt*time0*1e6*(gridx[i+1]-gridx[i])*(gridy[j+1]-gridy[j])*spr0*spr0*spr0;
		RBBL1DENS[i+j*ngpx]=RBBL1[i+j*ngpx]*dt*time0*1e6*(gridx[i+1]-gridx[i])*(gridy[j+1]-gridy[j])*spr0*spr0*spr0;
		RBBT1DENS[i+j*ngpx]=RBBT1[i+j*ngpx]*dt*time0*1e6*(gridx[i+1]-gridx[i])*(gridy[j+1]-gridy[j])*spr0*spr0*spr0;
		RBBLDENS[i+j*ngpx]=RBBL[i+j*ngpx]*dt*time0*1e6*(gridx[i+1]-gridx[i])*(gridy[j+1]-gridy[j])*spr0*spr0*spr0;
		RBBTDENS[i+j*ngpx]=RBBT[i+j*ngpx]*dt*time0*1e6*(gridx[i+1]-gridx[i])*(gridy[j+1]-gridy[j])*spr0*spr0*spr0;
	}

	ofstream oftp;
	oftp.open("output/RBBLout.txt");
	
	oftp<<setw(15)<<"X(nm)"<<setw(15)<<"Y(nm)"<<setw(15)<<"DE"<<setw(15)<<"RBBLnoT"<<setw(15)<<"RBBLwithT"<<setw(15)<<"RBBLAll"<<setw(15)<<"RBBLPC"<<endl;
	for(i=0;i<ngpx-1;i++)
	{
		for(j=0;j<ngpy-1;j++)
		{
			oftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
			oftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
			oftp<<setw(15)<<fabs(DE[i+j*ngpx]);
			oftp<<setw(15)<<fabs(RBBL0[i+j*ngpx]);
			oftp<<setw(15)<<fabs(RBBL1[i+j*ngpx]);
			oftp<<setw(15)<<fabs(RBBL[i+j*ngpx]);
			oftp<<setw(15)<<fabs(RBBLDENS[i+j*ngpx])<<endl;
		}
	}
	oftp.close();

	oftp.open("output/RBBTout.txt");
	
	oftp<<setw(15)<<"X(nm)"<<setw(15)<<"Y(nm)"<<setw(15)<<"DE"<<setw(15)<<"RBBTnoT"<<setw(15)<<"RBBTwithT"<<setw(15)<<"RBBTAll"<<setw(15)<<"RBBTPC"<<endl;
	for(i=0;i<ngpx-1;i++)
	{
		for(j=0;j<ngpy-1;j++)
		{
			oftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
			oftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
			oftp<<setw(15)<<fabs(DE[i+j*ngpx]);
			oftp<<setw(15)<<fabs(RBBT0[i+j*ngpx]);
			oftp<<setw(15)<<fabs(RBBT1[i+j*ngpx]);
			oftp<<setw(15)<<fabs(RBBT[i+j*ngpx]);
			oftp<<setw(15)<<fabs(RBBTDENS[i+j*ngpx])<<endl;
		}
	}
	oftp.close();
}

void DevSimulator::statnoise(void)
{
	int i,j,m;
	double fmin,fmax;
	int tempint;
	int tsn;
	int statnumber;

	tsn=ndt-starttime;
	if(tsn>MNTIMES)
	{
		cout<<" tsn has greater than max time steps MNTIMES! "<<endl;
		exit(0);
	}
	fmin=1.0/(1.0*tsn*dt*time0);
	fmax=1.0/(2.0*dt*time0);
	ifstream ftp;

	//read in every time's current
	ftp.open("output/currentarray.txt");
	for(i=0;i<tsn;i++)
	{
		ftp>>tempint;
		ftp>>idtime[i];
	}
	ftp>>idba;
	ftp.close();

	ftp.open("output/electronvelocityarray.txt");
	for(j=0;j<tsn;j++)
	{
		ftp>>tempint;
		for(i=0;i<statnoiseregionnumber;i++)
		{
			ftp>>ev[i][j];
			ev[i][j]=ev[i][j]/velo0;
		}
	}
	for(i=0;i<statnoiseregionnumber;i++)
	{
		ftp>>evaverage[i];
		evaverage[i]=evaverage[i]/velo0;
	}
	ftp.close();

	ftp.open("output/holevelocityarray.txt");
	for(j=0;j<tsn;j++)
	{
		ftp>>tempint;
		for(i=0;i<statnoiseregionnumber;i++)
		{
			ftp>>hv[i][j];
			hv[i][j]=hv[i][j]/velo0;
		}
	}
	for(i=0;i<statnoiseregionnumber;i++)
	{
		ftp>>hvaverage[i];
		hvaverage[i]=hvaverage[i]/velo0;
	}
	ftp.close();

	ftp.open("output/carrierdensityarray.txt");
	for(j=0;j<tsn;j++)
	{
		ftp>>tempint;
		for(i=0;i<statnoiseregionnumber;i++)
		{
			ftp>>cd[i][j];
			cd[i][j]=cd[i][j]/conc0;
		}
	}
	for(i=0;i<statnoiseregionnumber;i++)
	{
		ftp>>cdaverage[i];
		cdaverage[i]=cdaverage[i]/conc0;
	}
	ftp.close();

	//calculate for id(t)-idaverage
	for(i=0;i<tsn;i++)
	{
		didtime[i]=idtime[i]-idba;
	}
	for(i=0;i<statnoiseregionnumber;i++)
	{
		for(j=0;j<tsn;j++)
		{
			ev[i][j]=ev[i][j]-evaverage[i];
			hv[i][j]=hv[i][j]-hvaverage[i];
			cd[i][j]=cd[i][j]-cdaverage[i];
		}
	}
	//calculate for normalized Id noise 
	for(i=0;i<tsn;i++)
	{
		idtimeguiyi[i]=idtime[i]/idba;
	}
	//output normalized Id noise in output/idnoise.txt
	ofstream ff;
	ff.open("output/idnoise.txt");
	for(i=0;i<tsn;i++)
	{
		if((fmod(i+1,statoutmod)==0)||(i==tsn-1))
		{
			ff<<i*dt*time0<<" "<<idtimeguiyi[i]<<endl;
		}
	}
	ff.close();
	//calculate CI(t)
	for(m=0;m<tsn/2;m++)//m*dt=t in CI(t=m*dt)
	{
		for(i=0;i<tsn/2;i++)//i is the stat position
		{
			CI[m]=CI[m]+didtime[i+m]*didtime[i];
		}
		CI[m]=2.0*CI[m]/tsn;
	}
	for(j=0;j<statnoiseregionnumber;j++)
	{
		for(m=0;m<tsn/2;m++)//m*dt=t in CI(t=m*dt)
		{
			for(i=0;i<tsn/2;i++)//i is the stat position
			{
				CEV[j][m]=CEV[j][m]+ev[j][i+m]*ev[j][i];
				CHV[j][m]=CHV[j][m]+hv[j][i+m]*hv[j][i];
				CCD[j][m]=CCD[j][m]+cd[j][i+m]*cd[j][i];
			}
			CEV[j][m]=2.0*CEV[j][m]/tsn;
			CHV[j][m]=2.0*CHV[j][m]/tsn;
			CCD[j][m]=2.0*CCD[j][m]/tsn;
		}
	}
	//output CI(t)
	ofstream Ciout;
	Ciout.open("output/CI(t).txt");
	for(i=0;i<tsn/2;i++)
	{
		Ciout<<i*dt*time0<<" "<<CI[i]/CI[0]<<" "<<CI[i]<<" "<<i<<endl;
	}
	Ciout.close();

	ofstream Cevout;
	ofstream Chvout;
	ofstream Ccdout;
	Cevout.open("output/CEV(t).txt");
	Chvout.open("output/CHV(t).txt");
	Ccdout.open("output/CCD(t).txt");
	
	for(i=0;i<tsn/2;i++)
	{
		Cevout<<i*dt*time0<<" ";
		Chvout<<i*dt*time0<<" ";
		Ccdout<<i*dt*time0<<" ";
		for(j=0;j<statnoiseregionnumber;j++)
		{
			Cevout<<CEV[j][i]/CEV[j][0]<<" "<<CEV[j][i]<<" ";
			Chvout<<CHV[j][i]/CHV[j][0]<<" "<<CHV[j][i]<<" ";
			Ccdout<<CCD[j][i]/CCD[j][0]<<" "<<CCD[j][i]<<" ";
		}
		Cevout<<endl;
		Chvout<<endl;
		Ccdout<<endl;
	}
	Cevout.close();
	Chvout.close();
	Ccdout.close();

	//calculate SI(f)
	for(i=0;i<fhznumber;i++)
	{
		fhz[i]=fmin+i*(fmax-fmin)/fhznumber;//liner
//		fhz[i]=fmin*pow(10,i*(log10(fmax)-log10(fmin))/fhznumber);//log
		for(m=0;m<tsn/2;m++)
		{
			SI[i]=SI[i]+4.0*CI[m]*cos(2*PI*fhz[i]*m*dt*time0);
		}
		SI[i]=SI[i]*2.0/(1.0*tsn);
	}

	for(i=0;i<fhznumber;i++)
	{
		fhz[i]=fmin+i*(fmax-fmin)/fhznumber;//liner
//		fhz[i]=fmin*pow(10,i*(log10(fmax)-log10(fmin))/fhznumber);//log
		for(m=0;m<tsn/2;m++)
		{
			for(j=0;j<statnoiseregionnumber;j++)
			{
				SEV[j][i]=SEV[j][i]+4.0*CEV[j][m]*cos(2*PI*fhz[i]*m*dt*time0);
				SHV[j][i]=SHV[j][i]+4.0*CHV[j][m]*cos(2*PI*fhz[i]*m*dt*time0);
				SCD[j][i]=SCD[j][i]+4.0*CCD[j][m]*cos(2*PI*fhz[i]*m*dt*time0);
			}
		}
		SEV[j][i]=SEV[j][i]*2.0/(1.0*tsn)*velo0*velo0;
		SHV[j][i]=SHV[j][i]*2.0/(1.0*tsn)*velo0*velo0;
		SCD[j][i]=SCD[j][i]*2.0/(1.0*tsn)*conc0*conc0;
	}
	//output SI(f)
	ofstream Sout;
	Sout.open("output/SI(f).txt");
	Sout<<"fhznumber is "<<fhznumber<<" fmin is "<<fmin<<" fmax is "<<fmax<<endl;
	for(i=0;i<fhznumber;i++)
	{
		Sout<<fhz[i]<<" "<<SI[i]<<" "<<i<<endl;
	}
	Sout.close();

	ofstream Sevout;
	ofstream Shvout;
	ofstream Scdout;
	Sevout.open("output/SEV(f).txt");
	Shvout.open("output/SHV(f).txt");
	Scdout.open("output/SCD(f).txt");
	Sevout<<"fhznumber is "<<fhznumber<<" fmin is "<<fmin<<" fmax is "<<fmax<<endl;
	Shvout<<"fhznumber is "<<fhznumber<<" fmin is "<<fmin<<" fmax is "<<fmax<<endl;
	Scdout<<"fhznumber is "<<fhznumber<<" fmin is "<<fmin<<" fmax is "<<fmax<<endl;
	for(i=0;i<fhznumber;i++)
	{
		Sevout<<fhz[i]<<" ";
		Shvout<<fhz[i]<<" ";
		Scdout<<fhz[i]<<" ";
		for(j=0;j<statnoiseregionnumber;j++)
		{
			Sevout<<SEV[j][i]<<" ";
			Shvout<<SHV[j][i]<<" ";
			Scdout<<SCD[j][i]<<" ";
		}
		Sevout<<endl;
		Shvout<<endl;
		Scdout<<endl;
	}
	Sevout.close();
	Shvout.close();
	Scdout.close();

	return;
}

void DevSimulator::GETCUTPOSITION(void)
{
		int i,j;
// intput the cut position
	if(cutfl)
	{
		for(i=0;i<xcutnumber;i++)
		{
			positionx[i]=positionx[i]*1e-6/spr0;
			if(positionx[i]<gridx[0]||positionx[i]>gridx[ngpx-1])
			{
				cout<<"error region input of xcut borden!!!!"<<endl;
				exit(0);
			}
		}
		for(i=0;i<ycutnumber;i++)
		{
			positiony[i]=positiony[i]*1e-6/spr0;
			if(positiony[i]<gridy[0]||positiony[i]>gridy[ngpy-1])
			{
				cout<<"error region input of ycut borden!!!!"<<endl;
				exit(0);
			}
		}

// get the position of cut x
		for(j=0;j<xcutnumber;j++)
		{
			for(i=0;i<ngpx;i++)
			{
				if(j<xcutnumber-1)
				{
					if(gridx[i+1]>positionx[j]&&positionx[j]>=gridx[i])
					{
						ipcutx[j]=i;
						break;
					}
				}
				else
				{
					if(gridx[i+1]>=positionx[j]&&positionx[j]>gridx[i])
					{
						ipcutx[j]=i;
						break;
					}
				}
			}
		}
// get the position of cut y
		for(j=0;j<ycutnumber;j++)
		{
			for(i=0;i<ngpy;i++)
			{
				if(j<ycutnumber-1)
				{
					if(gridy[i+1]>positiony[j]&&positiony[j]>=gridy[i])
					{
						ipcuty[j]=i;
						break;
					}
				}
				else
				{
					if(gridy[i+1]>=positiony[j]&&positiony[j]>gridy[i])
					{
						ipcuty[j]=i;
						break;
					}
				}
			}
		}
	}

// get noise region position
		for(i=0;i<statnoiseregionnumber;i++)
		{
			xminnoise[i]=xminnoise[i]*1e-6/spr0;
			xmaxnoise[i]=xmaxnoise[i]*1e-6/spr0;
			yminnoise[i]=yminnoise[i]*1e-6/spr0;
			ymaxnoise[i]=ymaxnoise[i]*1e-6/spr0;
			if(xminnoise[i]<gridx[0]||xminnoise[i]>gridx[ngpx-1])
			{
				cout<<"error region input of noise region borden!!!!"<<endl;
				exit(0);
			}
			if(xmaxnoise[i]<gridx[0]||xmaxnoise[i]>gridx[ngpx-1])
			{
				cout<<"error region input of noise region borden!!!!"<<endl;
				exit(0);
			}
			if(xminnoise[i]>=xmaxnoise[i])
			{
				cout<<"xminnoise can't be greater than xmaxnoise!";
				exit(0);
			}

			if(yminnoise[i]<gridy[0]||yminnoise[i]>gridy[ngpy-1])
			{
				cout<<"error region input of noise region borden!!!!"<<endl;
				exit(0);
			}
			if(ymaxnoise[i]<gridy[0]||ymaxnoise[i]>gridy[ngpy-1])
			{
				cout<<"error region input of noise region borden!!!!"<<endl;
				exit(0);
			}
			if(yminnoise[i]>=ymaxnoise[i])
			{
				cout<<"yminnoise can't be greater than ymaxnoise!";
				exit(0);
			}
		}

		for(j=0;j<statnoiseregionnumber;j++)
		{
			for(i=0;i<ngpx-1;i++)
			{
				if(xminnoise[j]==gridx[0])
				{
					xminnoisegrid[j]=0;
					break;
				}
				else if(gridx[i+1]>=xminnoise[j]&&xminnoise[j]>gridx[i])
				{
					xminnoisegrid[j]=i;
					break;
				}
			}
			for(i=0;i<ngpx-1;i++)
			{
				if(gridx[i+1]>xmaxnoise[j]&&xmaxnoise[j]>=gridx[i])
				{
					xmaxnoisegrid[j]=i;
					break;
				}
				else if(xmaxnoise[j]==gridx[ngpx-2])
				{
					xmaxnoisegrid[j]=ngpx-2;
					break;
				}
			}
			for(i=0;i<ngpy-1;i++)
			{
				if(yminnoise[j]==gridy[0])
				{
					yminnoisegrid[j]=0;
					break;
				}
				else if(gridy[i+1]>=yminnoise[j]&&yminnoise[j]>gridy[i])
				{
					yminnoisegrid[j]=i;
					break;
				}
			}
			for(i=0;i<ngpy-1;i++)
			{
				if(gridy[i+1]>ymaxnoise[j]&&ymaxnoise[j]>=gridy[i])
				{
					ymaxnoisegrid[j]=i;
					break;
				}
				else if(ymaxnoise[j]==gridy[ngpx-2])
				{
					ymaxnoisegrid[j]=ngpy-2;
					break;
				}
			}
		}
	return;
}

void DevSimulator::OUTPUTMOTRULES(void)
{
		int i,j;
	ofstream ftp0;
	ftp0.open("motrules0.txt");
	for(i=0;i<ngpx;i++)
	{
		for(j=0;j<ngpy;j++)
		{
			ftp0<<motrules[0][i+j*ngpx]<<" ";
		}
		ftp0<<endl;
	}
	ftp0.close();

	ofstream ftp1;
	ftp1.open("motrules1.txt");
	for(i=0;i<ngpx;i++)
	{
		for(j=0;j<ngpy;j++)
		{
			ftp1<<motrules[1][i+j*ngpx]<<" ";
		}
		ftp1<<endl;
	}
	ftp1.close();

	ofstream ftp2;
	ftp2.open("motrules2.txt");
	for(i=0;i<ngpx;i++)
	{
		for(j=0;j<ngpy;j++)
		{
			ftp2<<motrules[2][i+j*ngpx]<<" ";
		}
		ftp2<<endl;
	}
	ftp2.close();

	ofstream ftp3;
	ftp3.open("motrules3.txt");
	for(i=0;i<ngpx;i++)
	{
		for(j=0;j<ngpy;j++)
		{
			ftp3<<motrules[3][i+j*ngpx]<<" ";
		}
		ftp3<<endl;
	}
	ftp3.close();
	return;
}

void DevSimulator::INIST(Partical &par)
{
//     process the initialize command!

//     initialize certain variables and simulation conditions

//____local variables
      int i,j,k,ncal,nhyp,idirec,ier,ij,iptype;
      int  ntit,ic;
      double  dd,rfield;

	  isymflag=false;

	  for(i=0;i<MNGP;i++)
	  {
		  for(j=0;j<4;j++)
		  {
			  ebandsum[i][j]=0;
			  ebandratio[i][j]=0;
		  }
	  }
	  for(i=0;i<MNGP;i++)
	  {
		  for(j=0;j<3;j++)
		  {
			  hbandsum[i][j]=0;
			  hbandratio[i][j]=0;
		  }
	  }

	  btbtTn=0;
	  btbtTp=0;
	  btbtLn=0;
	  btbtLp=0;

	  for(i=0;i<ngp;i++)
	  {
		  evelox[i]=0;
		  eveloy[i]=0;
		  eveloxsum[i]=0;
		  eveloysum[i]=0;
		  enumber[i]=0;

		  hvelox[i]=0;
		  hveloy[i]=0;
		  hveloxsum[i]=0;
		  hveloysum[i]=0;
		  hnumber[i]=0;

		  for(j=0;j<3;j++)
		  {
			  eveloxx[j][i]=0;
			  eveloyy[j][i]=0;
			  eveloxsumm[j][i]=0;
			  eveloysumm[j][i]=0;
			  enumberr[j][i]=0;

			  hveloxx[j][i]=0;
			  hveloyy[j][i]=0;
			  hveloxsumm[j][i]=0;
			  hveloysumm[j][i]=0;
			  hnumberr[j][i]=0;

			  eenergyy[j][i]=0;
			  eenergysumm[j][i]=0;
			  henergyy[j][i]=0;
			  henergysumm[j][i]=0;
		  }

		  mobility[i]=0;
		  ndensaver[i]=0;
		  pdensaver[i]=0;
		  ndenssum[i]=0;
		  pdenssum[i]=0;
		  numberndens[i]=0;
		  numberpdens[i]=0;

		  for(j=0;j<3;j++)
		  {
			  ndensaverr[j][i]=0;
			  pdensaverr[j][i]=0;
			  ndenssumm[j][i]=0;
			  pdenssumm[j][i]=0;
			  numberndenss[j][i]=0;
			  numberpdenss[j][i]=0;
			  NRDENS[j][i]=0;
			  PRDENS[j][i]=0;
		  }
		  eenergy[i]=0;
		  eenergysum[i]=0;
		  henergy[i]=0;
		  henergysum[i]=0;

	      RBBT0[i]=0;
		  RBBT1[i]=0;
		  RBBT[i]=0;
	      RBBL0[i]=0;
		  RBBL1[i]=0;
		  RBBL[i]=0;
		  RBBTDENS[i]=0;
		  RBBLDENS[i]=0;
		  RBBT0DENS[i]=0;
		  RBBT1DENS[i]=0;
		  RBBL0DENS[i]=0;
		  RBBL1DENS[i]=0;
		  RBBDENS[i]=0;
		  RDENS[i]=0;
		  iseEfn[i]=0;
		  iseEfp[i]=0;
		  iseEi[i]=0;
		  iseEc[i]=0;
		  iseEv[i]=0;
		  iseefx[i]=0;
		  iseefy[i]=0;
		  isepot[i]=0;
		  pot[i]=0;
		  edens[i]=0;
		  hdens[i]=0;
		  xfield[i]=0;
		  yfield[i]=0;
		  Efn[i]=0;
		  Efp[i]=0;
		  Ei[i]=0;
		  Ec[i]=0;
		  Ev[i]=0;
		  DE[i]=0;
		  Nc[i]=0;
		  Nv[i]=0;
		  dc[i]=0;
		  dv[i]=0;
	  }

//____card parameters
      testfl    =false;// lval(1,linum)
      odxfl     =false;// lval(2,linum)
      bulkfl    =false;// lval(3,linum)
      pvmfl     =false;// lval(4,linum)
      resumefl  =false;// lval(5,linum)
      wlbfl     =true;// lval(6,linum)
      ningfl    =false;// lval(7,linum)
      qyfl      =false;// lval(8,linum)
      sifl      =true;// lval(9,linum)
      gaasfl    =false;// lval(10,linum)
      rnlongfl  =false;// lval(11,linum)

//____check semiconductor type 
      ic=0;
      if(sifl)ic=ic+1;
      if(gaasfl)ic=ic+1;
      if(ic!=1)
	  {
		  cout<<"error INIST: specify semiconductor type";exit(0);
      }//endif

      if(bulkfl&&odxfl)
	  {
		  cout<<"error INIST: Bulk and 1D simulation together not possible";exit(0);
      }//endif


      temp =300;// dval(1,linum)

      if(gaasfl)
	  {
//_______gaas bulk parameters
//       lattice constant
         sia0 =5.64e-10;
//       mass density
         sirho=5.36e3;
//       sound velocity
         siul =5.24e3;
         siut =2.47e3;

      }
	  else if(sifl)
	  {//then
//_______si bulk parameters
//       lattice constant
		  if(material==0)
		  {
//          lattice constant
			sia0 =5.43e-10;
//          mass density
			sirho=2.33e3;
//          sound velocity
			siul =9.05e3;
			siut =9.05e3;
		  }
		  else if(material==1)
		  {
//          sia0 =5.43e-10;
			sia0 =5.658e-10;
//          mass density
//          sirho=2.33e3;
			sirho=5.32e3;
//          sound velocity
			siul =5.4e3;
			siut =3.2e3;
		  }
      }//endif

      cfilter=0.250;// dval(3,linum)

      idirec=111;// ival(1,linum)

      ncal=0;//ival(2,linum)
      nhyp=0;//ival(3,linum)
      ntit=0;//ival(4,linum)

//____only one initialize card
      if(inifl)
	  {
		  cout<<"error INIST: Too many initialize cards";exit(0);
      }
	  else
	  {
         inifl=true;
      }//endif

//____check that the temperature is not negativ or zero
      if(temp<=0.0)
	  {
		  cout<<"error INIST: invalid temperature";exit(0);
      }//endif

//____Ning experiment is 1D
      if(ningfl&&!odxfl)
	  {
		  cout<<"error INIST: ning requires one dimensional structure";exit(0);
      }//endif

//    get the material coefficients and normalize them

//    temperatur [K]
      T0 =temp;

//    energy [eV]/electon rest mass [kg]/Planck's constant [eVs] /
//    electron charge [As]
      eV0   =BOLTZ*T0;
      em0   =EM;
      hq0   =PLANCK;
      ec0   =EC;

//    momentum [eVs/m]/r-space [m]/k-space [1/m]/time [s] /
//    velocity [m/s]
      rmom0 =sqrt(em0/ec0*eV0);
      spr0  =hq0/rmom0;
      spk0  =1.0/spr0;
      time0 =hq0/eV0;
      velo0 =spr0/time0;
      cvr   =CLIGHT/velo0;

//    el. potential [V]/el. field [V/m]/concentration [1/m**3]
      pot0  =eV0;
      field0=pot0/spr0;
      conc0 =spk0*spk0*spk0;

//    mass density [kg/m**3]
      dens0 =em0*conc0;

//    deformation potential constant [eV/m]/scattering rate [1/s]
      dpc0  =field0;
      scrt0 =1.0/time0;

//____current [A/m]
      curr0 =ec0/(time0*spr0);

//    Si bulk parameters
      sia0=sia0/spr0;
      sirho=sirho/dens0;
      siul=siul/velo0;
      siut=siut/velo0;

//  temperature dependent band gap
	if(material==0)
	{
      if(T0<190.0)
	  {//then
         sieg=(1.170+1.059e-5*T0-6.05e-7*T0*T0)/eV0;
      }
	  else if(T0<250.0)
	  {//then
         sieg=(1.17850-9.025e-5*T0-3.05e-7*T0*T0)/eV0;
      }
	  else
	  {
         sieg=(1.2060-2.730e-4*T0)/eV0;
      }//endif
	}
	else if(material==1)
	{
	  sieg=0.70e0/eV0;
	}

//    Defines a0pi
      a0pi=TWOPI/sia0;

//    get relative dielectric constant and normalize (eps_vacuum <> 1 within the program !!!)
      eps[VACUUM]= 1.00/(4*PI*cvr*FSC);
      eps[OXIDE]= 3.90/(4*PI*cvr*FSC);
      if(gaasfl)
	  {//then
         eps[SILICON]=12.90/(4*PI*cvr*FSC);
      }
	  else if (sifl)
	  {//then
		  if(material==0)
		  {
			  eps[SILICON]=11.70/(4*PI*cvr*FSC);
		  }
		  else if(material==1)
		  {
			  eps[SILICON]=16.0/(4*PI*cvr*FSC);
		  }
      }//endif

//____bulk electric field
      xfieldbulk=0.0;
      yfieldbulk=0.0;
      zfieldbulk=0.0;

//____in case of one dimensional simulation it is possible to apply 
//     a homogeneous field in y-direction (perpendicular to the 1D-structure)
      if(odxfl) yfieldbulk=0;//dval(2,linum)/field0
         
//____Quantum yield simulation
      if(qyfl)
	  {//then
         if(!bulkfl)
		 {
			 cout<<"error INIST: Quantum yield only in bulk sim";exit(0);
         }//endif
         if(qyenergy<0.0)
		 {
			 cout<<"error INIST: Qunatum yield energy < 0";exit(0);
         }//endif
         qyenergy=3./eV0;//dval(4,linum)/eV0
      }//endif

//    initialize ierp,used for scalar list for Multiple Refresh
      for(iptype=0;iptype<NPARTYP;iptype++)
	  {//	  DO iptype=1,NPARTYP
         for(ij=0;ij<=MNGP;ij++)
		 {//DO ij=0,MNGP
            for(ier=0;ier<MNER;ier++)
			{//DO ier=1,MNER
               ierp[ier][ij][iptype]=-1;
			}//ENDDO
         }//ENDDO
      }//ENDDO

//    initialize fitb
	  par.fitb[0]=1.0;
      par.fitb[1]=0.20;
      par.fitb[2]=0.040;
      par.fitb[3]=0.080;
      par.fitb[4]=0.150;
      
      par.dopmin=1e22/conc0;
      par.dopmax=1e26/conc0;

      par.frickel=1.0;	

//    end of INIST
      return;
}

void DevSimulator::QUANTARRAY(void)
{
	int i,j,ij,ityp;
	double dx,dy,dx1,dy1,qdtemp;
	for(ij=0;ij<ngp;ij++)
	{
		INDTD(ij,i,j);
		qtarray[ij][0]=0.0e0;
		qtarray[ij][1]=0.0e0;
		qtarray[ij][2]=0.0e0;
		qtarray[ij][3]=0.0e0;
		qtarray[ij][4]=0.0e0;
		dx=0.0e0;
		dx1=0.0e0;
		dy=0.0e0;
		dy1=0.0e0;
		if(true)
		{
			if(true)//mat[ij]==SILICON)
			{
				if(i>0)			dx =(gridx[i]-gridx[i-1])*spr0;
				if(i<ngpx-1)	dx1=(gridx[i+1]-gridx[i])*spr0;
				if(j>0)			dy =(gridy[j]-gridy[j-1])*spr0;
				if(j<ngpy-1)	dy1=(gridy[j+1]-gridy[j])*spr0;

				if(i>0)
				{
					if(true)//mat[ij]==SILICON)
					{
						qdtemp=2.0e0/(dx*(dx+dx1));
						qtarray[ij][0]=qtarray[ij][0]-qdtemp;
						qtarray[ij][1]=qdtemp;
					}
				}
				if(i<ngpx-1)
				{
					if(true)//mat[ij]==SILICON)
					{
						qdtemp=2.0e0/(dx1*(dx+dx1));
						qtarray[ij][0]=qtarray[ij][0]-qdtemp;
						qtarray[ij][2]=qdtemp;
					}
				}
				if(j>0)
				{
					if(true)//mat[ij]==SILICON)
					{
						qdtemp=2.0e0/(dy*(dy+dy1));
						qtarray[ij][0]=qtarray[ij][0]-qdtemp;
						qtarray[ij][3]=qdtemp;
					}
				}
				if(j<ngpy-1)
				{
					if(true)//mat[ij]==SILICON)
					{
						qdtemp=2.0e0/(dy1*(dy+dy1));
						qtarray[ij][0]=qtarray[ij][0]-qdtemp;
						qtarray[ij][4]=qdtemp;
					}
				}
			}
		}
		qdtemp=qtarray[ij][0]+qtarray[ij][1]+qtarray[ij][2]+qtarray[ij][3]+qtarray[ij][4];
	}//end do
	return;
}

//quantum effect
void DevSimulator::QUANTPOT(int instat)
{
	int i,j,ij,k,l;
	double templn[5],ptemp,lndens[MNGP];
	double a;
	a=hq0*ec0/sqrt(12.0*0.26*em0*0.026*ec0);
	if(quantumeffectqpfl)	ptemp=26.0*hq0*hq0*ec0/(12.0*0.26*em0*0.1);
	else if(quantumeffectbmfl)	ptemp=hq0*hq0*ec0/(2.0*0.26*em0);
	else if(quantumeffectqnfl)	ptemp=hq0*hq0*ec0/(12.0*0.26*em0);

	//search Si region
	siregion=0;
	for(i=0;i<nreg;i++)
	{
		if(regtype[i]==SILICON)
		{
			siregion=i;
		}
	}

if(quantumeffectqpfl||quantumeffectbmfl||quantumeffectqnfl)
{

	templn[0]=0.0e0;
	templn[1]=0.0e0;
	templn[2]=0.0e0;
	templn[3]=0.0e0;
	templn[4]=0.0e0;

	double oxndens,tempndens;
	oxndens=1e15;
	tempndens=log(1e15);

	for(i=ibreg[siregion];i<iereg[siregion]+1;i++)
	{
		for(j=jbreg[siregion];j<jereg[siregion]+1;j++)
		{
			if(quantumeffectqpfl)		lndens[i+j*ngpx]=statpot[i+j*ngpx]*pot0/(cdt+1);
			else if(quantumeffectbmfl)	lndens[i+j*ngpx]=sqrt(fabs(spchar[i+j*ngpx])*conc0);
			else if(quantumeffectqnfl)	lndens[i+j*ngpx]=log(fabs(spchar[i+j*ngpx])*conc0);
		}
	}

	for(i=ibreg[siregion];i<iereg[siregion]+1;i++)
	{
		for(j=jbreg[siregion];j<jereg[siregion]+1;j++)
		{
			templn[0]=lndens[i+j*ngpx];
			templn[1]=lndens[i+j*ngpx-1];
			templn[2]=lndens[i+j*ngpx+1];
			templn[3]=lndens[i+(j-1)*ngpx];
			templn[4]=lndens[i+(j+1)*ngpx];

			if(i==ibreg[siregion])
			{
				templn[0]=lndens[i+j*ngpx];
				templn[1]=tempndens;
				templn[2]=lndens[i+j*ngpx+1];
				templn[3]=lndens[i+(j-1)*ngpx];
				templn[4]=lndens[i+(j+1)*ngpx];
			}
			if(i==iereg[siregion])
			{
				templn[0]=lndens[i+j*ngpx];
				templn[1]=lndens[i+j*ngpx-1];
				templn[2]=tempndens;
				templn[3]=lndens[i+(j-1)*ngpx];
				templn[4]=lndens[i+(j+1)*ngpx];
			}

			qpot[i+j*ngpx][PELEC]=ptemp*(templn[0]*qtarray[i+j*ngpx][0]
							      +templn[1]*qtarray[i+j*ngpx][1]
								  +templn[2]*qtarray[i+j*ngpx][2]
								  +templn[3]*qtarray[i+j*ngpx][3]
								  +templn[4]*qtarray[i+j*ngpx][4]);
			if(quantumeffectbmfl)
			{
				if(fabs(spchar[i+j*ngpx]*conc0)!=0)
				{
					qpot[i+j*ngpx][PELEC]=qpot[i+j*ngpx][PELEC]/lndens[i+j*ngpx];
				}
			}
			if(fabs(qpot[i+j*ngpx][PELEC])>1.0e-20)	continue;
			else	qpot[i+j*ngpx][PELEC]=0.0e0;
		}
	}
}
else if(quantumeffectfmfl)
{
	double rdy;
	double rdx;
	double rd;
	double rx,ry;
	for(ij=0;ij<ngp;ij++)
	{
		qpot[ij][PELEC]=0.0e0;
		qpot[ij][PHOLE]=0.0e0;
		lndens[ij]=statpot[ij]*pot0/(cdt+1-starttime);
//		lndens[ij]=pot[ij]*pot0;
	}
	for(j=0;j<ngpy;j++)
	{
		for(i=0;i<ngpx;i++)
		{
			if(true)//mat[i+j*ngpx]==SILICON)
			{
				for(k=0;k<ngpy;k++)
				{
					for(l=0;l<ngpx;l++)
					{
						if(true)//mat[l+k*ngpx]==SILICON)
						{
							rdy=fabs(gridy[k+1]-gridy[k])*spr0;
							rdx=fabs(gridx[l+1]-gridx[l])*spr0;
							ry=fabs(gridy[k+1]-gridx[j])*spr0;
							rx=fabs(gridx[l+1]-gridx[i])*spr0;
							rd=sqrt(rx*rx+ry*ry);
							qpot[i+j*ngpx][PELEC]+=1/(2*PI*a*a)*rdx*rdy*lndens[l+k*ngpx]*exp(-rd*rd/(2*a*a));
						}
					}
				}
				qpot[i+j*ngpx][PELEC]=qpot[i+j*ngpx][PELEC];//-lndens[i+j*ngpx];
			}
		}
	}
}

		ofstream ftpqp;
		ftpqp.open("output/quantumpotentialout.txt");

		for (j=0;j<ngpy;j++)
			for (i=0;i<ngpx;i++)
			{
				ftpqp<<gridx[i]*spr0*1e9<<"  ";
				ftpqp<<gridy[j]*spr0*1e9<<"  ";
				ftpqp<<qpot[j*ngpx+i][PELEC]<<"  ";
				ftpqp<<endl;
			}
		ftpqp.close();

	return;
}

void DevSimulator::QUANTFIELD(void)
{
	int i,j,ij,ityp;
	int gp;
	double dx,dy;

	for(ityp=PELEC;ityp<PHOLE+1;ityp++)
	{
		for(ij=0;ij<ngp;ij++)
		{
			xqf[ij][ityp]=0.0e0;
            yqf[ij][ityp]=0.0e0;
		}
	}

	for(ityp=PELEC;ityp<PHOLE+1;ityp++)
	{
		for(i=ibreg[siregion];i<iereg[siregion];i++)
		for(j=jbcont[2];j<jecont[2];j++)//notice
		{
			dx=gridx[i+1]-gridx[i];
			dy=gridy[j+1]-gridy[j];
			xqf[i+j*ngpx][ityp]=(qpot[i+j*ngpx][ityp]+qpot[i+j*ngpx+ngpx][ityp]
				          -qpot[i+j*ngpx+ngpx+1][ityp]-qpot[i+j*ngpx+1][ityp])
						 /(2.0e0*dx);
			yqf[i+j*ngpx][ityp]=(qpot[i+j*ngpx][ityp]-qpot[i+j*ngpx+ngpx][ityp]
						  -qpot[i+j*ngpx+ngpx+1][ityp]+qpot[i+j*ngpx+1][ityp])
						 /(2.0e0*dy);
		}
	}

	  	int itypp;
		ofstream ftpqf;
		ftpqf.open("output/quantumelectricfieldout.txt");

		for (j=0;j<ngpy-1;j++)
			for (i=0;i<ngpx-1;i++)
			{
				ftpqf<<gridx[i]*spr0*1e9<<"  ";
				ftpqf<<gridy[j]*spr0*1e9<<"  ";
//				for(itypp=PELEC;itypp<PHOLE+1;itypp++)
				for(itypp=PELEC;itypp<PELEC+1;itypp++)
				{
					ftpqf<<xqf[j*ngpx+i][itypp]*field0/100.0<<"  ";
					ftpqf<<yqf[j*ngpx+i][itypp]*field0/100.0<<"  ";
				}
				ftpqf<<endl;
			}
		ftpqf.close();
	return;
}
//-----------------------------------------------------------------------
void DevSimulator::BIAST(char *name,double voltage,double delvol,double stime, double etime)
{ 
//     process bias command
//____local variables
      //CHARACTER*10 name
      int icont, ic;
      //double voltage, delvol, stime, etime
//____voltages must be applied after structure card
      if (! strfl)
	  {
		  cout<<"erro !sttrfl;biast";
		  cout<<"error BIAST: apply voltage after structure card";exit(0);
      }//ENDIF
      voltage=voltage/pot0;
      delvol=delvol/pot0;
      stime=stime/time0;
      etime=etime/time0;
      if (stime>=etime)
	  {
		  cout<<"stime>=etime;biast";
		  cout<<"error BIAST: start time is larger than end time";exit(0);
      }//ENDIF
      
//____apply voltage at contact and add work function difference
      ic=0;
	  for(icont=0;icont<ncont;icont++)
	  {// DO icont=1, ncont
         if (strstr(contname[icont],name)!=NULL)
		 {//THEN
            ic=ic+1;
            contpot[icont]=voltage+phims[icont];
//        ... calculate difference of voltage
            st[icont]=stime;
            ft[icont]=etime;
            diffvol[icont]=0.0;
            diffvol[icont]=delvol/(ft[icont]-st[icont]);
         }//ENDIF
      }//ENDDO
      if (ic==0) 
	  {
		  cout<<"no contact";
		  cout<<"error BIAST: named contact not found";exit(0);
      }//ENDIF
//____end of BIAST
      return;
}     
         
//====
void DevSimulator::GETCONTPOT(void)
{
//    get timedependent contact voltage
      int idt, icont;
      double nowtime;
//____initialize contact potential
      for(icont=0;icont<ncont;icont++)
	  {
		for(idt=0;idt<ndt;idt++)
		{
			ctpotnow[icont][idt]=0.0;
		}
      }
//____assign timedependent contact potential
      for(icont=0;icont<ncont;icont++)
	  {
		for(idt=0;idt<ndt;idt++)
		{
		    nowtime=(idt+1)*dt;
		 if ((nowtime>st[icont])&&
             (nowtime<ft[icont])) 
		 {
            ctpotnow[icont][idt]=contpot[icont] 
					+diffvol[icont]*(nowtime-st[icont]);
		 }
         else if(nowtime>=ft[icont])
		 {
            ctpotnow[icont][idt]=contpot[icont]
				+diffvol[icont]*(ft[icont]-st[icont]);
		 }
         else
            ctpotnow[icont][idt]=contpot[icont];
		}
      }
//____end of GETCONTVOL
      return;
}

void DevSimulator::CONFI(int rdwr,char *file)
{

//      int rdwr,lu
//      char * file

//     Purpose: reads/stores the configuration space from/to file
//     --------
//     parameter:
//       rdwr   : 1 to read from,2 to store to file
//       lu     : logical unit number
//       file   : name of file

      int i,ipar;

  //____read from the file
      if(rdwr==1)
	  {
		 ifstream ftp;
         if (!strfl)
		 {
			 cout<<"error CONFI: Load confi after strucutre";exit(0);
         }

         if (masterfl)
		 {
//__________open the file
			ftp.open(file);
//__________number of particles
			ftp>>npar0;
			for(int ii=0;ii<NPARTYP;ii++)ftp>>npar[ii];

//          load all particles into the memory
			ifstream ftpp;
			ftpp.open(file);
			for(i=0;i<4;i++)	ftpp>>tempint;
			ftpp>>tempint;
			for(i=0;i<9;i++)	ftpp>>tempdouble;
			ftpp>>tempdouble;
			if(tempint!=0&&tempint!=1&&tempint!=2&&fmod(tempdouble,1)>0)	firstrcfl=true;
			else	firstrcfl=false;
			ftpp.close();

            for(ipar= 0;ipar<npar0;ipar++)
			{
				if(firstrcfl)
				{
					for(i=0;i<IVPP;i++)ftp>>ifield[i][ipar];
					ifield[IVPP][ipar]=0;
				}
				else if(!firstrcfl)
				{
					ftp>>ifield[IVPP][ipar];
					for(i=0;i<IVPP;i++)ftp>>ifield[i][ipar];
				}
				else
				{

				}
				ifield[2][ipar]-=1;
				ifield[1][ipar]-=1;
				ifield[0][ipar]-=1;
				for(i=0;i<DVPP;i++)		ftp>>dfield[i][ipar];
			}

//__________state of random generator (in the case of parallel mode only
//          the master random generator is initialized with old variables)
//          don not read of random generator

            ftp>>iseedl;
			ftp>>iseedl;

//__________close the file
			ftp.close();
         }
      }

//____write to the file
	  else if(rdwr==2)
	  {//THEN
			ofstream ftp;
//__________open the file		
			ftp.open(file);
//__________number of particles
			ftp<<npar0<<" ";
			for(int ii=0;ii<NPARTYP;ii++)ftp<<npar[ii]<<" ";
			ftp<<endl;
            for(ipar=0;ipar<npar0;ipar++)
			{
				ftp<<ifield[IVPP][ipar]<<" ";
				for(i=0;i<IVPP;i++)ftp<<(ifield[i][ipar]+1)<<" ";
				for(i=0;i<DVPP;i++)ftp<<dfield[i][ipar]<<" ";
				ftp<<endl;
			}
//__________integers for state of random generator
			ftp<<iseedl<<" ";
			ftp<<iseedl<<" ";
//__________close the file
			ftp.close();

      }//ENDIF
	  if(rdwr==3)
	  {
		 ifstream ftp;
         if (!strfl)
		 {
			 cout<<"error CONFI: Load confi after strucutre";exit(0);
         }

         if (masterfl)
		 {
//__________open the file
			ftp.open(file);
//__________number of particles
			ftp>>npar0;
			for(int ii=0;ii<NPARTYP;ii++)ftp>>npar[ii];

//             ..load all particles into the memory
               for(ipar= 0;ipar<npar0;ipar++)
			   {
				   for(i=0;i<IVPP;i++)ftp>>ifield[i][ipar];
				   for(i=0;i<DVPP;i++)ftp>>dfield[i][ipar];
               }



//__________state of random generator (in the case of parallel mode only
//           the master random generator is initialized with old variables)
//           don not read of random generator
            ftp>>iseedl;
			ftp>>iseedl;

//__________close the file
			ftp.close();
         }
      }
//____end of CONFI
      return;
}
//=====

void DevSimulator::STRUC(int rdwr,char *file)
{
 //     int rdwr,lu
 //     char * file
//     Purpose: reads/stores the structure from/to file
//     --------
//     parameter:
//       rdwr   : 1 to read from,2 to store to file
//       lu     : logical unit number
//       file   : name of file

      int ij,i,j,iptype,ier,ipsc;
//____read from the file
      if (rdwr==1)
	  {
		ifstream ftp;
//_______structure from file or build
         if(strfl)
		 {
			 cout<<"error STRUC: structure already defined";exit(0);
         }
         if (grifl)
		 {
			 cout<<"error STRUC: grid already defined";exit(0);
         }

//_______initialize certain flags
         strfl=true;
         grifl=true;
		 int inttemp;
		 char inttempchar;

//__________open the file
			ftp.open(file);
//__________number of grid points
            ftp>>ngp;
			ftp>>ngpx;
			ftp>>ngpy;
			ftp>>nbwp;
			if (ngp>MNGP)
			{
				cout<<"error STRUC: too many grid points";exit(0);
            }
            if (ngpx>MNGPX)
			{
				cout<<"error STRUC: too many grid points in x-dir";exit(0);
            }
            if (ngpy>MNGPY)
			{
				cout<<"error STRUC: too many grid points in y-dir";exit(0); 
            }

			for(i=0;i<ngpx;i++)ftp>>gridx[i];
            for(j=0;j<ngpy;j++)ftp>>gridy[j];
            for(ij=0;ij<ngp;ij++)
				ftp>>quadarea[ij];
            for(i=0;i<NPARTYP;i++)
				for(ij=0;ij<ngp;ij++)	ftp>>boxarea[ij][i];
            for(ij=0;ij<ngp;ij++)ftp>>mat[ij];
            for(ij=0;ij<ngp;ij++){ftp>>inttemp;numreg[ij]=inttemp-1;}
            ftp>>nreg;
            for(i=0;i<nreg;i++)
			{ftp>>inttemp;ibreg[i]=inttemp-1;}
            for(i=0;i<nreg;i++)
			{ftp>>inttemp;iereg[i]=inttemp-1;}
            for(i=0;i<nreg;i++)
			{ftp>>inttemp;jbreg[i]=inttemp-1;}
            for(i=0;i<nreg;i++)
			{ftp>>inttemp;jereg[i]=inttemp-1;}
            for(i=0;i<nreg;i++)ftp>>regtype[i];
            for(i=0;i<nreg;i++)ftp>>regname[i];

//          for(i=0;i<nreg;i++){ftp>>inttemp;selfforcefl[i]=((inttemp!='F')&&(inttemp!=0));}
//          ftp>>inttemp;selffl=((inttemp!='F')&&(inttemp!=0));
//			start	transfer F to 0 and T to 1
            for(i=0;i<nreg;i++){ftp>>inttempchar;selfforcefl[i]=((inttempchar!='F')&&(inttempchar!=0));}
            ftp>>inttempchar;selffl=((inttempchar!='F')&&(inttempchar!=0));
//			end	transfer F to 0 and T to 1

            ftp>>ncont;
            for(i=0;i<ncont;i++)
			{ftp>>inttemp;ibcont[i]=inttemp-1;}
            for(i=0;i<ncont;i++)
			{ftp>>inttemp;iecont[i]=inttemp-1;}
            for(i=0;i<ncont;i++)
			{ftp>>inttemp;jbcont[i]=inttemp-1;}
            for(i=0;i<ncont;i++)
			{ftp>>inttemp;jecont[i]=inttemp-1;}
            for(i=0;i<ncont;i++)ftp>>conttype[i];
            for(i=0;i<ncont;i++)ftp>>contname[i];
            for(i=0;i<ncont;i++)ftp>>contpos[i];
            for(i=0;i<ncont;i++)ftp>>contpot[i];
            for(i=0;i<ncont;i++)ftp>>phims[i];
            for(ij=0;ij<ngp;ij++){ftp>>inttemp;gridcont[ij]=inttemp-1;}

//          for(i=0;i<ncont;i++)
//			{ftp>>inttemp;equisicfl[i]=((inttemp!='F')&&(inttemp!=0));}
//			start	transfer F to 0 and T to 1
			for(i=0;i<ncont;i++)
			{ftp>>inttempchar;equisicfl[i]=((inttempchar!='F')&&(inttempchar!=0));}//READ (lu) (equisicfl(i),i=1,ncont)
//			end	transfer F to 0 and T to 1

            ftp>>nqss;
            for(i=0;i<nqss;i++)
			{ftp>>inttemp;ibqss[i]=inttemp-1;}
            for(i=0;i<nqss;i++)
			{ftp>>inttemp;ieqss[i]=inttemp-1;}
            for(i=0;i<nqss;i++)
			{ftp>>inttemp;jbqss[i]=inttemp-1;}
            for(i=0;i<nqss;i++)
			{ftp>>inttemp;jeqss[i]=inttemp-1;}
            for(i=0;i<nqss;i++)ftp>>qssdens[i];
            for(ij=0;ij<ngp;ij++)
			{
               for(i=0;i<4;i++)ftp>>cont[i][ij];
               for(i=0;i<4;i++)ftp>>motrules[i][ij];
               for(i=MNBWP-nbwp;i<=MNBWP+nbwp;i++)ftp>>dpe[i][ij];
            }
            for(ij=0;ij<ngp;ij++)
				for(i=0;i<4;i++)	ftp>>donor[i][ij];
            for(ij=0;ij<ngp;ij++)
				for(i=0;i<4;i++)	ftp>>accep[i][ij];
            for(ij=0;ij<ngp;ij++)ftp>>doprhspe[ij];

			for(ij=0;ij<ngp;ij++)
			{
				if(mat[ij]!=SILICON)
				{
					for(i=0;i<4;i++)	donor[i][ij]=0;;
					for(i=0;i<4;i++)	accep[i][ij]=0;
					doprhspe[ij]=0;
				}
			}

            ftp>>eni;

//			ftp>>inttemp;mulfl0=((inttemp!='F')&&(inttemp!=0));
//          for(i=0;i<NPARTYP;i++){ftp>>inttemp;mulfl[i]=((inttemp!='F')&&(inttemp!=0));}
//			for(i=0;i<NPARTYP;i++){ftp>>inttemp;enermulfl[i]=((inttemp!='F')&&(inttemp!=0));}
//			start	transfer F to 0 and T to 1
			ftp>>inttempchar;mulfl0=((inttempchar!='F')&&(inttempchar!=0));
            for(i=0;i<NPARTYP;i++){ftp>>inttempchar;mulfl[i]=((inttempchar!='F')&&(inttempchar!=0));}
			for(i=0;i<NPARTYP;i++){ftp>>inttempchar;enermulfl[i]=((inttempchar!='F')&&(inttempchar!=0));}
//			end	transfer F to 0 and T to 1

			ftp>>rpar;
			ftp>>rquad;
			ftp>>npsc;
            ftp>>pckill;
			ftp>>pcaver;
            for(iptype=0;iptype<NPARTYP;iptype++)
			{
               ftp>>nerrf[MNGP][iptype];
			   for(ij=0;ij<ngp;ij++)ftp>>nerrf[ij][iptype];
               ftp>>erfmax[MNGP][iptype];
			   for(ij=0;ij<ngp;ij++)ftp>>erfmax[ij][iptype];
               ij=MNGP;
			   for(ier=0;ier<=nerrf[ij][iptype];ier++)
			   {
				   ftp>>inttemp;
				   ierp[ier][ij][iptype]=inttemp-1;
			   }
			   for(ij=0;ij<ngp;ij++)
			   {
                  for(ier=0;ier<=nerrf[ij][iptype];ier++)
				  {
					  ftp>>inttemp;
					  ierp[ier][ij][iptype]=inttemp-1;
				  }
               }
            }
            for(ipsc=0;ipsc<npsc;ipsc++) ftp>>ndpper[ipsc];

//          ftp>>inttemp;bulkfl=((inttemp!='F')&&(inttemp!=0));
//			ftp>>inttemp;odxfl=((inttemp!='F')&&(inttemp!=0));
//			ftp>>inttemp;ningfl=((inttemp!='F')&&(inttemp!=0));
//			start	transfer F to 0 and T to 1
            ftp>>inttempchar;bulkfl=((inttempchar!='F')&&(inttempchar!=0));
			ftp>>inttempchar;odxfl=((inttempchar!='F')&&(inttempchar!=0));
			ftp>>inttempchar;ningfl=((inttempchar!='F')&&(inttempchar!=0));
//			end	transfer F to 0 and T to 1

            ftp>>inttemp;ibinj=inttemp-1;
			ftp>>inttemp;ieinj=inttemp-1;
			ftp>>inttemp;jbinj=inttemp-1;
			ftp>>inttemp;jeinj=inttemp-1;
			ftp>>injdir;
			ftp>>nparinj;
			ftp>>injcur;
			ftp>>itypinj;

//			ftp>>inttemp;injfl=((inttemp!='F')&&(inttemp!=0));
//			ftp>>inttemp;cbhfl=((inttemp!='F')&&(inttemp!=0));
//			start	transfer F to 0 and T to 1
			ftp>>inttemp;injfl=((inttemp!='F')&&(inttemp!=0));
			ftp>>inttemp;cbhfl=((inttemp!='F')&&(inttemp!=0));
//			end	transfer F to 0 and T to 1

//__________close the file
			ftp.close();
      }//end if

	  else if(rdwr==2)
	  {
		 ofstream ftp;
//       only master can save structure
         if (masterfl) 
		 {
			ftp.open(file);
//__________number of grid points
            ftp<<ngp;
			ftp<<ngpx;
			ftp<<ngpy;
			ftp<<nbwp;
			
            for(i=0;i<ngpx;i++)ftp<<gridx[i];
            for(j=0;j<ngpy;j++)ftp<<gridy[j];
            for(ij=0;ij<ngp;ij++)ftp<<quadarea[ij];
            for(i=0;i<NPARTYP;i++)
				for(ij=0;ij<ngp;ij++)	ftp<<boxarea[ij][i];
            for(ij=0;ij<ngp;ij++)ftp<<mat[ij];
            for(ij=0;ij<ngp;ij++)ftp<<numreg[ij];
            ftp<<nreg;
            for(i=0;i<nreg;i++)ftp<<ibreg[i];
            for(i=0;i<nreg;i++)ftp<<iereg[i];
            for(i=0;i<nreg;i++)ftp<<jbreg[i];
            for(i=0;i<nreg;i++)ftp<<jereg[i];
            for(i=0;i<nreg;i++)ftp<<regtype[i];
            for(i=0;i<nreg;i++)ftp<<regname[i];
            for(i=0;i<nreg;i++)ftp<<selfforcefl[i];
            ftp<<selffl;
            ftp<<ncont;
            for(i=0;i<ncont;i++)ftp<<ibcont[i];
            for(i=0;i<ncont;i++)ftp<<iecont[i];
            for(i=0;i<ncont;i++)ftp<<jbcont[i];
            for(i=0;i<ncont;i++)ftp<<jecont[i];
            for(i=0;i<ncont;i++)ftp<<conttype[i];
            for(i=0;i<ncont;i++)ftp<<contname[i];
            for(i=0;i<ncont;i++)ftp<<contpos[i];
            for(i=0;i<ncont;i++)ftp<<contpot[i];
            for(i=0;i<ncont;i++)ftp<<phims[i];
            for(ij=0;ij<ngp;ij++)ftp<<gridcont[ij];
            for(i=0;i<ncont;i++)ftp<<equisicfl[i];
            ftp<<nqss;
            for(i=0;i<nqss;i++)ftp<<ibqss[i];
            for(i=0;i<nqss;i++)ftp<<ieqss[i];
            for(i=0;i<nqss;i++)ftp<<jbqss[i];
            for(i=0;i<nqss;i++)ftp<<jeqss[i];
            for(i=0;i<nqss;i++)ftp<<qssdens[i];
            for(ij=0;ij<ngp;ij++)
			{
               for(i=0;i<4;i++)ftp<<cont[i][ij];
               for(i=0;i<4;i++)ftp<<motrules[i][ij];
               for(i=MNBWP-nbwp;i<=MNBWP+nbwp;i++)ftp<<dpe[i][ij];
            }
            for(ij=0;ij<ngp;ij++)
				for(i=0;i<4;i++)	ftp<<donor[i][ij];
            for(ij=0;ij<ngp;ij++)
				for(i=0;i<4;i++)	ftp<<accep[i][ij];
            for(ij=0;ij<ngp;ij++)ftp<<doprhspe[ij];
            ftp<<eni;
			ftp<<mulfl0;
            for(i=0;i<NPARTYP;i++)ftp<<mulfl[i];
			for(i=0;i<NPARTYP;i++)ftp<<enermulfl[i];
			ftp<<rpar;
			ftp<<rquad;
			ftp<<npsc;
            ftp<<pckill;
			ftp<<pcaver;
            for(iptype=0;iptype<NPARTYP;iptype++)
			{
               for(ij=0;ij<ngp;ij++)ftp<<nerrf[ij][iptype];
               for(ij=0;ij<ngp;ij++)ftp<<erfmax[ij][iptype];
               for(ij=0;ij<ngp;ij++)
			   {
                  for(ier=0;ier<=nerrf[ij][iptype];ier++)
					  ftp<<ierp[ier][ij][iptype];
               }
            }
            for(ipsc=0;ipsc<npsc;ipsc++)ftp<<ndpper[ipsc];
            ftp<<bulkfl;
			ftp<<odxfl;
			ftp<<ningfl;
            ftp<<ibinj;
			ftp<<ieinj;
			ftp<<jbinj;
			ftp<<jeinj;
			ftp<<injdir;
			ftp<<nparinj;
			ftp<<injcur;
			ftp<<itypinj;
			ftp<<injfl;
			ftp<<cbhfl;

//__________close the file
			ftp.close();
         }//ENDIF masterfl
      }//ENDIF rwdr

//____end of STRUC
      return;
}
//-----------------------------------------------------------------------
void DevSimulator::CALCPOT(bool wrfl)
{
//____local variables
      int icont,i,j,ij,ireal,icpu,ierr;
      double vdum[MNGP],sum,tcpu,drhspe[MNGP];
	  static double contpotl[MNGP],cc,dopeb[MNGP];
	  static bool flag=false;


//____calculate right hand side of Poisson equation
      for(ij=0;ij<ngp;ij++)
	  {
         drhspe[ij]=doprhspe[ij]+spchar[ij];
      }

	  if(!flag)
	  {
		  flag=true;
//____cal. contact voltages (Silicon contacts)
//    get doping on grid points
      for(ij=0;ij<ngp;ij++)
	  {
         dopeb[ij]=0.0;
      }

      for(i=0;i<ngpx-1;i++)
	  {
        for(j=0;j<ngpy-1;j++)
		{
         ij=INDS(i,j);
         if (mat[ij]==SILICON)
		 {
            dopeb[ij]=dopeb[ij]
                  +0.25*quadarea[ij]*(donor[UL][ij]-accep[UL][ij]);
            dopeb[ij+1]=dopeb[ij+1]
                  +0.25*quadarea[ij]*(donor[LL][ij]-accep[LL][ij]);
            dopeb[ij+ngpx]=dopeb[ij+ngpx]
                  +0.25*quadarea[ij]*(donor[UR][ij]-accep[UR][ij]);
            dopeb[ij+ngpx+1]=dopeb[ij+ngpx+1]
                  +0.25*quadarea[ij]*(donor[LR][ij]-accep[LR][ij]);
         }
		}
      }
      for(ij=0;ij<ngp;ij++)
	  {
         if (boxarea[ij][PELEC]!=0.0)
		 {
            dopeb[ij]=dopeb[ij]/boxarea[ij][PELEC];
         }
      }
	  for(ij=0;ij<ngp;ij++)
	  {
         contpotl[ij]=0.0;
      }
      for(icont=0;icont<ncont;icont++)
	  {
         for(i=ibcont[icont];i<=iecont[icont];i++)
		 {
          for(j=jbcont[icont];j<=jecont[icont];j++)
		  {
            ij=INDS(i,j);
            contpotl[ij]=0.0;//ctpotnow[icont][cdt];
//          cal. contact voltage offset for silicon
//          doping should BE CONSTANT within and around silicon contacts !!!!
            if (conttype[icont]==SICONT)
			{
               cc=dopeb[ij]/(eni*2.0);
               if (cc>0.0)
			   {
                  contpotl[ij]=log(cc+sqrt(1.0+cc*cc));
				  //ctpotnow[icont][cdt]+
               }
			   else
			   {
                  contpotl[ij]=-log(sqrt(1.0+cc*cc)-cc);
				  //ctpotnow[icont][cdt] -
               }//ENDIF
            }//ENDIF
          }//ENDDO
         }//ENDDO
      }//ENDDO
	}//dopeb  and Sicontact just need cal once
	  

//____set Dirichlet points (contacts)
      for(icont=0;icont<ncont;icont++)
	  {
         for(i=ibcont[icont];i<=iecont[icont];i++)
		 {
           for(j=jbcont[icont];j<=jecont[icont];j++)
		   {
            ij=INDS(i,j);
            drhspe[ij]=ctpotnow[icont][cdt];
			  //contpotl[ij];
//            cal. contact voltage offset for silicon
//            doping should BE CONSTANT within and around silicon contacts !!!!
			if (conttype[icont]==SICONT)
			drhspe[ij]=ctpotnow[icont][cdt]+contpotl[ij];
		   }//ENDDO
         }//ENDDO
      }//ENDDO
    
//____solve: L-matrix times vdum-vector=drhspe-vector
//     remember L has 1's on the main diagonal
//     discrete Poisson equation
      for(ij=0;ij<ngp;ij++)
	  {
         vdum[ij]=drhspe[ij];
         for(i=1;i<=MIN(ij,nbwp);i++)
		 {
            vdum[ij]=vdum[ij]-dpe[MNBWP-i][ij]*vdum[ij-i];
         }//ENDDO
      }//ENDDO
      
//____perform back substitution (solve : R-matrix times potential-vector
//                                       =vdum-vector)

      for(ij=ngp;ij>=1;ij--)
	  {
         sum=vdum[ij-1];
         for(i=1;i<=MIN(nbwp,ngp-ij);i++)
		 {
            sum=sum-dpe[i+MNBWP][ij-1]*pot[ij-1+i];
         }//ENDDO
         pot[ij-1]=sum/dpe[0+MNBWP][ij-1];
      }//ENDDO

//____statistic for potential
      for(ij=0;ij<ngp;ij++)
	  {
		  if(cdt>=starttime)
		  {
			statpot[ij]=statpot[ij]+pot[ij];
		  }
      }//ENDDO
//____end of CALCPOT
return;
}

void DevSimulator::CALCRHO(bool wrfl)
{
//     Calculate space charge and concentration of ionized dopants
//____local variables
      int ij,i,j,k,ipar,ireal,icpu,ierr;
      double tcpu,dxrr,dyrr,pc;
	  int carriertype;

//____initialize space charge (electron and holes only)
      for(ij=0;ij<ngp;ij++)
	  {
         spchar[ij]=0.0;
         quaddens[ij]=0.0;
		 quadndens[ij]=0.0;
		 quadpdens[ij]=0.0;
		 for(i=0;i<3;i++)
		 {
			 quaddenss[i][ij]=0;
			 quadndenss[i][ij]=0;
			 quadpdenss[i][ij]=0;
		 }
      }//ENDDO

//____get particle space charge
      for(ipar=0;ipar<npar0;ipar++)
	  {
         ij=ifield[2][ipar];
         INDTD(ij,i,j);
         pc=dfield[5][ipar];
         dxrr=(dfield[0][ipar]-gridx[i])/(gridx[i+1]-gridx[i]);
         dyrr=(dfield[1][ipar]-gridy[j])/(gridy[j+1]-gridy[j]);
         spchar[ij]			+=pc*(1.0-dxrr)*(1.0-dyrr);
         spchar[ij+1]		+=pc*(       dxrr)*(1.0-dyrr);
         spchar[ij+ngpx]	+=pc*(1.0-dxrr)*(       dyrr);
         spchar[ij+ngpx+1]	+=pc*(       dxrr)*(       dyrr);
         quaddens[ij]+=fabs(pc);

		 if(pc<0)
		 {
			 quadndens[ij]+=fabs(pc);
		 }
		 else if(pc>0)
		 {
			 quadpdens[ij]+=fabs(pc);
		 }

		 carriertype=ifield[IVPP][ipar];
		 if(carriertype==0||carriertype==1||carriertype==2)
		 {
			quaddenss[carriertype][ij]+=fabs(pc);
			if(pc<0)
			{
				quadndenss[carriertype][ij]+=fabs(pc);
			}
			else if(pc>0)
			{
				quadpdenss[carriertype][ij]+=fabs(pc);
			}
		 }
      }//ENDDO
//____normalize particle density within a quadrant
     
      for(i=0;i<ngpx-1;i++)
	  {
       for(j=0;j<ngpy-1;j++)
	   {
         ij=INDS(i,j);
         quaddens[i+j*ngpx]=quaddens[i+j*ngpx]/quadarea[i+j*ngpx];
		 quadndens[i+j*ngpx]=quadndens[i+j*ngpx]/quadarea[i+j*ngpx];
		 quadpdens[i+j*ngpx]=quadpdens[i+j*ngpx]/quadarea[i+j*ngpx];
		 for(k=0;k<3;k++)
		 {
			 quaddenss[k][i+j*ngpx]=quaddenss[k][i+j*ngpx]/quadarea[i+j*ngpx];
			 quadndenss[k][i+j*ngpx]=quadndenss[k][i+j*ngpx]/quadarea[i+j*ngpx];
			 quadpdenss[k][i+j*ngpx]=quadpdenss[k][i+j*ngpx]/quadarea[i+j*ngpx];
		 }
       }//ENDDO
      }//ENDDO
//____end of CALCRHO
      return;
}

void DevSimulator::CALCFIELD(void)
{ 
//____local variables
      int ij,i,j;
      double dx,dy,dxcen,dycen,xfcl,xfcr,yfcu,yfcl;
      double epsup,epslow,epsright,epsleft,epscen;
      double xful,xfur,xflr,xfll;
      double yful,yfur,yflr,yfll;
      double potl[MNGP];

//____get potential
      if(potgalfl)
	  {
         for(ij=0;ij<ngp;ij++)
		 {
            potl[ij]=galpot[ij];
         }//ENDDO
      }
	  else
	  {
         for(ij=0;ij<ngp;ij++)
		 {
            potl[ij]=pot[ij];
         }//ENDDO
      }//ENDIF

//____loop over all quadrants
      for(ij=0;ij<ngp;ij++)
	  {
         INDTD(ij,i,j);
         if (mat[ij]!=VOID)
		 {
            dx=gridx[i+1]-gridx[i];
            dy=gridy[j+1]-gridy[j];
            xfield[ij]=( potl[ij]+potl[ij+ngpx]
                          -potl[ij+ngpx+1]-potl[ij+1]
                      +imagepot[ij]+imagepot[ij+ngpx]
                          -imagepot[ij+ngpx+1]-imagepot[ij+1])
                      /(2.0*dx);
            yfield[ij]=( potl[ij]-potl[ij+ngpx]
                          -potl[ij+ngpx+1]+potl[ij+1]
                      +imagepot[ij]-imagepot[ij+ngpx]
                          -imagepot[ij+ngpx+1]+imagepot[ij+1])
                      /(2.0*dy);
         }
		 else
		 {
            xfield[ij]=0.0;
            yfield[ij]=0.0;
         }//ENDIF
      }//ENDDO

//____field coefficients in y-direction
      if (selffl)
	  {
         for(ij=0;ij<ngp;ij++)
		 {
            INDTD(ij,i,j);
            if (mat[ij]!=VOID)
			{
               dxcen=gridx[i+1]-gridx[i];
               dycen=gridy[j+1]-gridy[j];
               epscen=eps[mat[ij]];
               xfcl=(potl[ij]-potl[ij+1])/dxcen;
               xfcr=(potl[ij+ngpx]-potl[ij+ngpx+1])/dxcen;
               yfcu=(potl[ij]-potl[ij+ngpx])/dycen;
               yfcl=(potl[ij+1]-potl[ij+ngpx+1])/dycen;

//          ..upper two corners
               //if ((i==1)||(mat[ij-1]==VOID)){// THEN
               if ((i==0)||(mat[ij-1]==VOID))
			   {
                  xful=-xfcl;
                  xfur=-xfcr;
                  epsup=epscen;
               }
			   else
			   {
                  xful=(potl[ij-1]-potl[ij]) 
                   /(gridx[i]-gridx[i-1]);
                  xfur=(potl[ij-1+ngpx]-potl[ij+ngpx]) 
                   /(gridx[i]-gridx[i-1]);
                  epsup=eps[mat[ij-1]];
               }//ENDIF
               xfself[UL][ij]=0.50*(xfcl+xful*epsup/epscen);
               xfself[UR][ij]=0.50*(xfcr+xfur*epsup/epscen);

//          ..lower two corners
               //if ((i==ngpx-1)||(mat[ij+1]==VOID))
			   if ((i==ngpx-2)||(mat[ij+1]==VOID))
			   {
                  xfll=-xfcl;
                  xflr=-xfcr;
                  epslow=epscen;
               }
			   else
			   {
                  xfll=(potl[ij+1]-potl[ij+2]) 
                   /(gridx[i+2]-gridx[i+1]);
                  xflr=(potl[ij+1+ngpx]-potl[ij+2+ngpx]) 
                   /(gridx[i+2]-gridx[i+1]);
                  epslow=eps[mat[ij+1]];
               }//ENDIF
               xfself[LL][ij]=0.50*(xfcl+xfll*epslow/epscen);
               xfself[LR][ij]=0.50*(xfcr+xflr*epslow/epscen);

//          ..left two corners
               //if ((j==1)||(mat[ij-ngpx]==VOID)){// THEN
               if ((j==0)||(mat[ij-ngpx]==VOID))
			   {
                  yful=-yfcu;
                  yfll=-yfcl;
                  epsleft=epscen;
               }
			   else
			   {
                  yful=(potl[ij-ngpx]-potl[ij]) 
                   /(gridy[j]-gridy[j-1]);
                  yfll=(potl[ij+1-ngpx]-potl[ij+1]) 
                   /(gridy[j]-gridy[j-1]);
                  epsleft=eps[mat[ij-ngpx]];
               }//ENDIF
               yfself[UL][ij]=0.50*(yfcu+yful*epsleft/epscen);
               yfself[LL][ij]=0.50*(yfcl+yfll*epsleft/epscen);

//          ..right two corners
               //if ((j==ngpy-1)||(mat[ij+ngpx]==VOID)){// THEN
               if ((j==ngpy-2)||(mat[ij+ngpx]==VOID))
			   {
                  yfur=- yfcu;
                  yflr=- yfcl;
                  epsright=epscen;
               }
			   else
			   {
                  yfur=(potl[ij+ngpx]-potl[ij+2*ngpx]) 
                   /(gridy[j+2]-gridy[j+1]);
                  yflr=(potl[ij+1+ngpx]-potl[ij+1+2*ngpx]) 
                   /(gridy[j+2]-gridy[j+1]);
                  epsright=eps[mat[ij+ngpx]];
               }//ENDIF
               yfself[UR][ij]=0.50*(yfcu+yfur*epsright/epscen);
               yfself[LR][ij]=0.50*(yfcl+yflr*epsright/epscen);
   
            }
			else
			{
               xfself[UL][ij]=0.0;
               xfself[UR][ij]=0.0;
               xfself[LR][ij]=0.0;
               xfself[LL][ij]=0.0;
               yfself[UL][ij]=0.0;
               yfself[UR][ij]=0.0;
               yfself[LR][ij]=0.0;
               yfself[LL][ij]=0.0;
            }//ENDIF
         }//ENDDO
      }//ENDIF

//____end of CALCFIELD
      return;
}
void DevSimulator::CALCRS(void)
{
//____local variables
      int icont,i,j,ij,ireal,icpu,ierr;;
      double vdum[MNGP],sum,tcpu,dx,dy;
      double rhsrs[MNGP];
      double rs[MNGP];

//____loop over all contacts
      for(icont=0;icont<ncont;icont++)
	  {
//_______calculate right hand side of Laplace equation
         for(ij=0;ij<ngp;ij++)
		 {
            rhsrs[ij]=0.0;
            xgradrs[ij][icont]=0.0;
            ygradrs[ij][icont]=0.0;
         }//ENDDO

//_______set Dirichlet points (contacts)
         for(i=ibcont[icont];i<=iecont[icont];i++)
		 {
           for(j=jbcont[icont];j<=jecont[icont];j++)
		   {
            ij=INDS(i,j);
            rhsrs[ij]=1.0;
		   }//ENDDO
         }//ENDDO
    
//_______solve: L-matrix times vdum-vector=drhspe-vector
//        remember L has 1's on the main diagonal
//        discrete Poisson equation
         for(ij=0;ij<ngp;ij++)
		 {
            vdum[ij]=rhsrs[ij];
            for(i=1;i<=MIN(ij,nbwp);i++)//////be carful
			{//DO i=1,MIN(ij-1,nbwp)
               vdum[ij]=vdum[ij]-dpe[MNBWP-i][ij]*vdum[ij-i];
            }//ENDDO
         }//ENDDO
      
//_______perform back substitution (solve : R-matrix times vector
//                                       =vdum-vector)
         for(ij=ngp;ij>=1;ij--)
		 {//1DO ij=ngp,1,-1
            sum=vdum[ij-1];
            for(i=1;i<=MIN(nbwp,ngp-ij);i++)///be careful
			{//DO i=1,MIN(nbwp,ngp-ij)
               sum=sum-dpe[i+MNBWP][ij-1]*rs[ij-1+i];
            }//ENDDO
            rs[ij-1]=sum/dpe[0+MNBWP][ij-1];
         }//ENDDO

//_______loop over all quadrants
         for(ij=0;ij<ngp;ij++)
		 {
            INDTD(ij,i,j);
            if (mat[ij]!=VOID)
			{
               dx=gridx[i+1]-gridx[i];
               dy=gridy[j+1]-gridy[j];
               xgradrs[ij][icont]=(rs[ij]+rs[ij+ngpx]
                        -rs[ij+ngpx+1]-rs[ij+1])/(2.0*dx);
               ygradrs[ij][icont]=(rs[ij]-rs[ij+ngpx]
                        -rs[ij+ngpx+1]+rs[ij+1])/(2.0*dy);
            }
			else
			{
               xgradrs[ij][icont]=0.0;
               ygradrs[ij][icont]=0.0;
            }//ENDIF
         }//ENDDO

      }//ENDDO
//____end of CALCRS
      return;
}
//=====

void DevSimulator::CALCCAP(void)
{
//    calculate a capacitance between two contacts (for R-S)
//____local variables
      int icont,jcont,ij,i,j;//,INDS

//____loop over all contacts
      for(jcont=0;jcont<ncont;jcont++)
	  {
		for(icont=0;icont<ncont;icont++)
		{
			cap[icont][jcont]=0.0;
			for(i=0;i<ngpx-1;i++)
			{
				for(j=0;j<ngpy-1;j++)
				{
					ij=INDS(i,j);
					if (mat[ij]!=VOID)
					{
						cap[icont][jcont]=cap[icont][jcont]  
					        +eps[mat[ij]]*quadarea[ij]
							*(xgradrs[ij][icont]*xgradrs[ij][jcont]
							+ygradrs[ij][icont]*ygradrs[ij][jcont]);
					}//ENDIF
				}//ENDDO
			}//ENDDO
		}//ENDDO
      }//ENDDO
//____end of CALCCAP
      return;
}
//-----------------------------------------------------------------------
void DevSimulator::RUN(Partical &par,Band &bd)
{
//____local variables 
      bool wrfl,opalfl,averfl;
      int idt,ireal,icpu,ierr,instat,injb;
      //int stridestat,stridemr,stridejb;
      double tcpu,rcpu,rtime;
	  
	  int nparold0,nparolde,nparoldh;
	  int ii;
	  int i,j,k;
	  int ipt;
//____local variables
      bool avfl;
//	  int stridestat,stridemr,stridejb;

//____only one run card
      if(runfl)
	  {
		  cout<<"error RUNST: Too many run cards !!";exit(0);
      }
	  else
	  {
         runfl=true;
      }

      if(ningfl&&!mulfl0)
	  {
		  cout<<"error RUNST: ning only possible with refresh";exit(0);
      }

//____card parameter
//	  ndt=10000;//ival(1,linum)
//	  outmod=10;//ival(2,linum)
//	  stridestat=10;//ival(3,linum)
//	  stridemr=20;//ival(4,linum)
//	  stridejb=10;//ival(5,linum)
//	  dt=2.5e-16/time0;//dval(1,linum)/time0
      edfemax[PELEC]=5.0/eV0;//dval(2,linum)/eV0
      edfemax[PHOLE]=5.0/eV0;//dval(3,linum)/eV0
      edfemax[POXEL]=5.0/eV0;//dval(4,linum)/eV0
      facstopii=0.1;//dval(5,linum)
      avfl    =true;//lval(1,linum)
      potgalfl=false;//lval(2,linum)
      jbscfl  =false;//lval(3,linum)
      elecisubfl=false;//lval(4,linum)
      holeisubfl=false;//lval(5,linum)
      ghdmfl=false;//lval(6,linum)
      elecghdmfl=false;//lval(7,linum)
      holeghdmfl=false;//lval(8,linum)
      jbedffl=false;//lval(9,linum)
      quantfl=false;//lval(10,linum)
      gatedepfl=false;//lval(11,linum)

//____check input

//____cal. image potentials in the oxide (only x-direction)
      GETIMAGEPOT(bd);

	  ofstream ftp;
	  ftp.open("output/cur.txt");

	  ofstream currentout;
	  currentout.open("output/currentarray.txt");

	  ofstream evout;
	  evout.open("output/electronvelocityarray.txt");

	  ofstream hvout;
	  hvout.open("output/holevelocityarray.txt");

	  ofstream cdout;
	  cdout.open("output/carrierdensityarray.txt");

//____initials
      catime=0.0;
      instat=0;
      injb  =0;
      rtime=0.0;
	  eni=bd.eni;



//____get Ramo-Shockley testfunctions
      CALCRS();
//____get contact voltage
      GETCONTPOT();
//____get capacitance
      CALCCAP();
      
//____get particle charge
      CALCRHO(true);

	  GETISEDATA();
	  if(btbtfl)
	  {
		if(methodforgetisedatafl)
		{
			DATAFROMISETOISEGOOD();
		}
		DATAFROMISETOMC();
//		CALCFERMI();
		CALCBTBTDENS();
	  }

//____Calculate the quantum potential
	  QUANTARRAY();


//____calculate potential for the first time
      cdt=0;

//xiazl      if(!potgalfl)CALCPOT(true);
	  if(selfconsistentfl&&!potgalfl&&cdt>=starttime)	
	  {
		  CALCPOT(true);
		  CALCFIELD();
	  }

//____calculate electric field
//xiazl      CALCFIELD();

//____initialize particle simulation
      EINIT(par);

//____Initialize averages for electrons and holes
      if(avfl)AVERINIT();

//=====for each observation time
      idt=0;
      while(idt<ndt)
	  {
         cdt=idt;
         if (masterfl)
		 {
            if(CONVER(injb)) ndt=idt;
         }

         if((fmod(idt+1,outmod)==0)||(idt==ndt))
		 {
			 wrfl=true;
         }
		 else
		 {
            wrfl=false;
         }

//_______initialize just before scattering statistics
         if(fmod(idt,stridejb)==0)ZEROJB();

//_______get number of statistic evaluations
         if(avfl&&(fmod(idt+1,stridestat)==0))
		 {
            instat=instat+1;
            averfl=true;
         }
		 else
		 {
            averfl=false;
         }

//_______set propagation time
		 CALCRDENS();
		 recombination();
		 par.GENBTBTDENS(*this,bd);
         SETDT();

//_______Propagate electrons during this time step
         ELEC(idt,wrfl,averfl,par,bd);

		 if(cdt>=starttime)
		 {
			 for(i=0;i<statnoiseregionnumber;i++)
			 {
				 if(evnumber[i][cdt-starttime]!=0)
				 {
					 ev[i][cdt-starttime]=evsum[i][cdt-starttime]/evnumber[i][cdt-starttime];
					 hv[i][cdt-starttime]=hvsum[i][cdt-starttime]/hvnumber[i][cdt-starttime];
				 }
				 else
				 {
					 ev[i][cdt-starttime]=0;
					 hv[i][cdt-starttime]=0;
				 }
			 }
			 for(k=0;k<statnoiseregionnumber;k++)
			 {
				 for(i=xminnoisegrid[k];i<xmaxnoisegrid[k]+1;i++)
				 {
					 for(j=yminnoisegrid[k];j<ymaxnoisegrid[k]+1;j++)
					 {
						 cd[k][cdt-starttime]+=quaddens[i+j*ngpx];
					 }
				 }
				 cd[k][cdt-starttime]=cd[k][cdt-starttime]
					 /((xmaxnoisegrid[k]-xminnoisegrid[k]+1)
					  *(ymaxnoisegrid[k]-yminnoisegrid[k]+1));
			 }
		 }

		 if(wrfl)
		 {

		 }

//_______Add just before scattering statistics (ADDJB contains mesage passing.
//       Thus do not include in real time calculation)
         if(fmod(idt+1,stridejb)==0)
		 {
            injb=injb + 1;
            ADDJB(stridejb,bd);
         }

//_______get particle charge
//xiazl	 CALCRHO(wrfl);
		 CALCRHO(true);
		 if(btbtfl)
		 {
//			CALCFERMI();
//			CALCBTBTDENS();
		 }
		 if(cdt>=starttime)
		 {
			 for(i=0;i<ngp;i++)
			 {
				 ndenssum[i]+=quadndens[i];
				 numberndens[i]++;
				 pdenssum[i]+=quadpdens[i];
				 numberpdens[i]++;
				 for(k=0;k<3;k++)
				 {
					 ndenssumm[k][i]+=quadndenss[k][i];
					 numberndenss[k][i]++;
					 pdenssumm[k][i]+=quadpdenss[k][i];
					 numberpdenss[k][i]++;
				 }
			 }
		 }

//_______Multiple Refresh
		if(mulfl0&&(fmod(idt+1,stridemr)==0)) //xiazl
		{
//			REFRESH(wrfl,bd,par);//xiazl
			for(i=0;i<3;i++)	nparbufall[i]=0;
			for(i=0;i<3;i++)	nparbufe[i]=0;
			for(i=0;i<3;i++)	nparbufh[i]=0;

			for(ii=0;ii<npar0;ii++)
			{
				ipt=bd.partyp[bd.ibt[ifield[0][ii]]];
				if(ifield[IVPP][ii]==0)
				{
					for(i=0;i<IVPP+1;i++)	ifieldbuf[0][i][nparbufall[0]]=ifield[i][ii];
					for(i=0;i<DVPP;i++)		dfieldbuf[0][i][nparbufall[0]]=dfield[i][ii];
					if(ipt==PELEC)	nparbufe[0]++;
					else if(ipt==PHOLE)	nparbufh[0]++;
					nparbufall[0]++;
					if(nparbufe[0]+nparbufh[0]!=nparbufall[0])	cout<<endl<<"ERROR!"<<endl;
				}
				else if(ifield[IVPP][ii]==1)
				{
					for(i=0;i<IVPP+1;i++)	ifieldbuf[1][i][nparbufall[1]]=ifield[i][ii];
					for(i=0;i<DVPP;i++)		dfieldbuf[1][i][nparbufall[1]]=dfield[i][ii];
					if(ipt==PELEC)	nparbufe[1]++;
					else if(ipt==PHOLE)	nparbufh[1]++;
					nparbufall[1]++;
					if(nparbufe[1]+nparbufh[1]!=nparbufall[1])	cout<<endl<<"ERROR!"<<endl;
				}
				else if(ifield[IVPP][ii]==2)
				{
					for(i=0;i<IVPP+1;i++)	ifieldbuf[2][i][nparbufall[2]]=ifield[i][ii];
					for(i=0;i<DVPP;i++)		dfieldbuf[2][i][nparbufall[2]]=dfield[i][ii];
					if(ipt==PELEC)	nparbufe[2]++;
					else if(ipt==PHOLE)	nparbufh[2]++;
					nparbufall[2]++;
					if(nparbufe[2]+nparbufh[2]!=nparbufall[2])	cout<<endl<<"ERROR!"<<endl;
				}
			}

			for(j=0;j<3;j++)
			{
				npar0=nparbufall[j];
				npar[PELEC]=nparbufe[j];
				npar[PHOLE]=nparbufh[j];
				for(ii=0;ii<npar0;ii++)
				{
					for(i=0;i<IVPP+1;i++)	ifield[i][ii]=ifieldbuf[j][i][ii];
					for(i=0;i<DVPP;i++)		dfield[i][ii]=dfieldbuf[j][i][ii];
				}
				REFRESH(wrfl,bd,par);//xiazl

				for(ii=0;ii<npar0;ii++)
				{
					for(i=0;i<IVPP+1;i++)	ifieldbuf[j][i][ii]=ifield[i][ii];
					for(i=0;i<DVPP;i++)		dfieldbuf[j][i][ii]=dfield[i][ii];
				}
				nparbufall[j]=npar0;
				nparbufe[j]=npar[PELEC];
				nparbufh[j]=npar[PHOLE];
			}
			npar0=0;npar[PELEC]=0;npar[PHOLE]=0;
			nparold0=0;nparolde=0;nparoldh=0;
			for(j=0;j<3;j++)
			{
				npar0+=nparbufall[j];
				npar[PELEC]+=nparbufe[j];
				npar[PHOLE]+=nparbufh[j];
				for(ii=nparold0;ii<npar0;ii++)
				{
					for(i=0;i<IVPP+1;i++)	ifield[i][ii]=ifieldbuf[j][i][ii-nparold0];
					for(i=0;i<DVPP;i++)		dfield[i][ii]=dfieldbuf[j][i][ii-nparold0];
				}
				nparold0+=nparbufall[j];
				nparolde+=nparbufe[j];
				nparoldh+=nparbufh[j];
			}
		}

//_______calculate potential 
//xiazl  if(!potgalfl)CALCPOT(wrfl);

		 if((quantumeffectqpfl||quantumeffectfmfl||quantumeffectbmfl||quantumeffectqnfl)&&cdt>=starttime)//averfl)
		 {
			 CALCRHO(true);
			 QUANTPOT(instat);
			 QUANTFIELD();
		 }

//_______calculate electric field
//xiazl  if(!potgalfl)CALCFIELD();
		 if(selfconsistentfl&&!potgalfl&&cdt>=starttime)	
		 {
			 CALCPOT(true);
			 CALCFIELD();
		 }

//_______simulation time duration
         catime=catime+dt;

		 if(cdt>=starttime)
		 {
			for(i=0;i<ncont;i++)
			{
				if(strstr(contname[i],"DRAIN")!=NULL)
				{
					currentout<<cdt-starttime<<" "<<currentarray[i][idt]*curr0*1e-2
							  <<" "<<currentarrayy[0][i][idt]*curr0*1e-2
							  <<" "<<currentarrayy[1][i][idt]*curr0*1e-2
							  <<" "<<currentarrayy[2][i][idt]*curr0*1e-2<<endl;
				}
			}
			evout<<cdt-starttime<<" ";
			hvout<<cdt-starttime<<" ";
			cdout<<cdt-starttime<<" ";

			for(i=0;i<statnoiseregionnumber;i++)
			{
				evout<<ev[i][cdt-starttime]*velo0<<" ";
				hvout<<hv[i][cdt-starttime]*velo0<<" ";
				cdout<<cd[i][cdt-starttime]*conc0/1e6<<" ";

				evaverage[i]+=ev[i][cdt-starttime];
				hvaverage[i]+=hv[i][cdt-starttime];
				cdaverage[i]+=cd[i][cdt-starttime];
			}
			evout<<endl;
			hvout<<endl;
			cdout<<endl;
		 }

         idt=idt+1;
		 btbtLn=0;
		 btbtTn=0;
		 btbtLp=0;
		 btbtTp=0;
		 for(i=0;i<npar0;i++)
		 {
			 ipt=bd.partyp[bd.ibt[ifield[0][i]]];
			 if(ipt==PELEC)
			 {
				if(ifield[IVPP][i]==1)	btbtLn++;
				else if(ifield[IVPP][i]==2)	btbtTn++;
			 }
			 else if(ipt==PHOLE)
			 {
				if(ifield[IVPP][i]==1)	btbtLp++;
				else if(ifield[IVPP][i]==2)	btbtTp++;
			 }
		 }
		 cout<<"------"<<setw(10)<<idt<<"      Step numbers=    "<<setw(10)<<ndt<<endl<<endl
			 <<"                 "<<setw(10)<<npar[0]-btbtLn-btbtTn<<setw(10)<<npar[1]-btbtLp-btbtTp
			 <<setw(10)<<btbtLn<<setw(10)<<btbtLp<<setw(10)<<btbtTn<<setw(10)<<btbtTp<<endl<<endl;
		 ftp<<idt<<endl;
		 if(wrfl) 
		 {
			 cout<<"idt= "<<setw(10)<<idt<<endl<<endl;
			 ftp<<"idt= "<<setw(10)<<idt<<endl<<endl;

			for(i=0;i<ncont;i++)
			{
			
				ftp<<setw(8)<<contname[i]<<" "<<"RM current= "<<endl<<endl
					<<"             "
					<<setw(15)<<curcont[i]*curr0*1e-2
					<<setw(15)<<curcontt[0][i]*curr0*1e-2
					<<setw(15)<<curcontt[1][i]*curr0*1e-2
					<<setw(15)<<curcontt[2][i]*curr0*1e-2
					<<" "<<" A/cm"<<endl<<endl;
				cout<<setw(8)<<contname[i]<<" "<<"RM current= "<<endl<<endl
					<<"             "
					<<setw(15)<<curcont[i]*curr0*1e-2
					<<setw(15)<<curcontt[0][i]*curr0*1e-2
					<<setw(15)<<curcontt[1][i]*curr0*1e-2
					<<setw(15)<<curcontt[2][i]*curr0*1e-2
					<<" "<<" A/cm"<<endl<<endl;

				ftp<<setw(8)<<contname[i]<<" "<<"catch current= "<<endl<<endl
				   <<"             "
				   <<setw(15)<<(dgen[i]-dcatch[i])/dt*curr0*1e-2
				   <<setw(15)<<(dgenn[0][i]-dcatchh[0][i])/dt*curr0*1e-2
			       <<setw(15)<<(dgenn[1][i]-dcatchh[1][i])/dt*curr0*1e-2
				   <<setw(15)<<(dgenn[2][i]-dcatchh[2][i])/dt*curr0*1e-2
				   <<" "<<" A/cm"<<endl<<endl;
				cout<<setw(8)<<contname[i]<<" "<<"catch current= "<<endl<<endl
					<<"             "
					<<setw(15)<<(dgen[i]-dcatch[i])/dt*curr0*1e-2
				    <<setw(15)<<(dgenn[0][i]-dcatchh[0][i])/dt*curr0*1e-2
			    	<<setw(15)<<(dgenn[1][i]-dcatchh[1][i])/dt*curr0*1e-2
				    <<setw(15)<<(dgenn[2][i]-dcatchh[2][i])/dt*curr0*1e-2
			        <<" "<<" A/cm"<<endl<<endl;
			}
		 }

      }//ENDDO (Big While)

	  if(avfl)
	  {
		  for(int ij=0;ij<ngp;ij++)
		  {
            statpot[ij]=statpot[ij]/(double)(ndt-starttime);
		  }
		  if (instat> 1)
		  {
				ftp<<endl<<"The average current is "<<endl<<endl;
				cout<<endl<<"The average current is "<<endl<<endl;
			 for(int icont=0;icont<ncont;icont++)
			 {
				meancurcont[icont]=meancurcont[icont]/(double)(ndt-starttime);
				dcurrent[icont]=dcurrent[icont]/((double)(ndt-starttime)*dt);
				for(j=0;j<3;j++)	meancurcontt[j][icont]=meancurcontt[j][icont]/(double)(ndt-starttime);
				for(j=0;j<3;j++)	dcurrentt[j][icont]=dcurrentt[j][icont]/((double)(ndt-starttime)*dt);
				squcurcont[icont]=squcurcont[icont]/(double)(ndt-starttime)/(double)(ndt-starttime)
					-meancurcont[icont]*meancurcont[icont];
				if (meancurcont[icont]!=0.0)
				{
					ftp<<setw(10)<<contname[icont]<<setw(20)<<"   meancur= "<<endl<<endl
					   <<"             "
					   <<setw(15)<<meancurcont[icont]*curr0*1e-2
					   <<setw(15)<<meancurcontt[0][icont]*curr0*1e-2
					   <<setw(15)<<meancurcontt[1][icont]*curr0*1e-2
					   <<setw(15)<<meancurcontt[2][icont]*curr0*1e-2
					   <<" "<<" A/cm"<<endl<<endl;
					cout<<setw(10)<<contname[icont]<<setw(20)<<"   meancur= "<<endl<<endl
						<<"             "
						<<setw(15)<<meancurcont[icont]*curr0*1e-2
						<<setw(15)<<meancurcontt[0][icont]*curr0*1e-2
					    <<setw(15)<<meancurcontt[1][icont]*curr0*1e-2
					    <<setw(15)<<meancurcontt[2][icont]*curr0*1e-2
						<<" "<<" A/cm"<<endl<<endl;
					ftp<<setw(10)<<contname[icont]<<setw(20)<<"catch current= "<<endl<<endl
					   <<"             "
					   <<setw(15)<<dcurrent[icont]*curr0*1e-2
					   <<setw(15)<<dcurrentt[0][icont]*curr0*1e-2
				       <<setw(15)<<dcurrentt[1][icont]*curr0*1e-2
					   <<setw(15)<<dcurrentt[2][icont]*curr0*1e-2
					   <<" "<<" A/cm"<<endl<<endl;
					cout<<setw(10)<<contname[icont]<<setw(20)<<"catch current= "<<endl<<endl
						<<"             "
						<<setw(15)<<dcurrent[icont]*curr0*1e-2
						<<setw(15)<<dcurrentt[0][icont]*curr0*1e-2
						<<setw(15)<<dcurrentt[1][icont]*curr0*1e-2
						<<setw(15)<<dcurrentt[2][icont]*curr0*1e-2
						<<" "<<" A/cm"<<endl<<endl;
		
		        }//ENDIF
			 }//ENDDO
		  }//ENDIF
	  }//ENDIF

	  for(i=0;i<statnoiseregionnumber;i++)
	  {
		  evaverage[i]=evaverage[i]/(ndt-starttime);
		  hvaverage[i]=hvaverage[i]/(ndt-starttime);
		  cdaverage[i]=cdaverage[i]/(ndt-starttime);
		  evout<<evaverage[i]*velo0<<" ";
		  hvout<<hvaverage[i]*velo0<<" ";
		  cdout<<cdaverage[i]*conc0/1e6<<" ";
	  }
	  for(i=0;i<ngp;i++)
	  {
		  if(mat[i]==SILICON)
		  {
			 evelox[i]=eveloxsum[i]/enumber[i];
			 eveloy[i]=eveloysum[i]/enumber[i];
			 hvelox[i]=hveloxsum[i]/hnumber[i];
			 hveloy[i]=hveloysum[i]/hnumber[i];

			 for(j=0;j<3;j++)
			 {
				if(enumberr[j][i]==0)
				{
					eveloxx[j][i]=0;
					eveloyy[j][i]=0;
				}
				else
				{
					eveloxx[j][i]=eveloxsumm[j][i]/enumberr[j][i];
					eveloyy[j][i]=eveloysumm[j][i]/enumberr[j][i];
				}
				if(hnumberr[j][i]==0)
				{
					hveloxx[j][i]=0;
					hveloyy[j][i]=0;
				}
				else
				{
					hveloxx[j][i]=hveloxsumm[j][i]/hnumberr[j][i];
					hveloyy[j][i]=hveloysumm[j][i]/hnumberr[j][i];
				}
			 }

			 ndensaver[i]=ndenssum[i]/numberndens[i];
			 pdensaver[i]=pdenssum[i]/numberpdens[i];

			 for(j=0;j<3;j++)
			 {
				 if(numberndenss[j][i]==0)	ndensaverr[j][i]=0;
				 else	ndensaverr[j][i]=ndenssumm[j][i]/numberndenss[j][i];
				 if(numberpdenss[j][i]==0)	pdensaverr[j][i]=0;
				 else	pdensaverr[j][i]=pdenssumm[j][i]/numberpdenss[j][i];
			 }

			 eenergy[i]=eenergysum[i]/enumber[i];
			 henergy[i]=henergysum[i]/hnumber[i];

			 for(j=0;j<3;j++)
			 {
				if(enumberr[j][i]==0)
				{
					eenergyy[j][i]=0;
				}
				else
				{
					eenergyy[j][i]=eenergysumm[j][i]/enumberr[j][i];
				}
				if(hnumberr[j][i]==0)
				{
					henergyy[j][i]=0;
				}
				else
				{
					henergyy[j][i]=henergysumm[j][i]/hnumberr[j][i];
				}
			 }
		  }
	  }
	  double ebandallsum=0,hbandallsum=0;
	  for(i=0;i<ngp;i++)
	  {
		  if(mat[i]==SILICON)
		  {
			  ebandallsum=ebandsum[i][0]+ebandsum[i][1]+ebandsum[i][2]+ebandsum[i][3];
			  ebandratio[i][0]=ebandsum[i][0]/ebandallsum*100.0;
			  ebandratio[i][1]=ebandsum[i][1]/ebandallsum*100.0;
			  ebandratio[i][2]=ebandsum[i][2]/ebandallsum*100.0;
			  ebandratio[i][3]=ebandsum[i][3]/ebandallsum*100.0;
			  hbandallsum=hbandsum[i][4]+hbandsum[i][5]+hbandsum[i][6];
			  hbandratio[i][4]=hbandsum[i][4]/hbandallsum*100.0;
			  hbandratio[i][5]=hbandsum[i][5]/hbandallsum*100.0;
			  hbandratio[i][6]=hbandsum[i][6]/hbandallsum*100.0;
		  }
	  }
      //if(avfl)AVEREVAL(instat,injb);
		output();
		outputcut();
		OUTPUTSCATTER();

//____end of RUN
		for(i=0;i<ncont;i++)
		{
			if(strstr(contname[i],"DRAIN")!=NULL)
			{
				for(k=starttime;k<ndt;k++)
				{
					averagecurrent[i]=averagecurrent[i]+currentarray[i][k];
					for(j=0;j<3;j++)	averagecurrentt[j][i]=averagecurrentt[j][i]+currentarrayy[j][i][k];
				}
				averagecurrent[i]=averagecurrent[i]/(ndt-starttime);
				for(j=0;j<3;j++)	averagecurrentt[j][i]=averagecurrentt[j][i]/(ndt-starttime);
				currentout<<endl
						  <<"                   "
						  <<averagecurrent[i]*curr0*1e-2<<" "
						  <<averagecurrentt[0][i]*curr0*1e-2<<" "
						  <<averagecurrentt[1][i]*curr0*1e-2<<" "
						  <<averagecurrentt[2][i]*curr0*1e-2<<" "<<endl;
			}
		}
	  currentout.close();
	  ftp.close();
      return;
//end of RUN
}

void DevSimulator::GETIMAGEPOT(Band &bd)
{
//    Purpose: cal. image potential in oxide (only x-dir)
//____local variables
      int ij,i,j,ireg;

      for(ij=0;ij<ngp;ij++)
	  {
         imagepot[ij]=0.0;
      }

      for(ireg=0;ireg<nreg;ireg++)
	  {
         if(regtype[ireg]==OXIDE)
		 {
            if((iereg[ireg]-ibreg[ireg])<2)
			{
				cout<<"error GETIMAGEPOT: Oxide must be resolved with more than 2 gridpoints";exit(0);
            }
            for(j=jbreg[ireg];j<=jereg[ireg];j++)
			{
              for(i=ibreg[ireg]+1;i<=iereg[ireg]-1;i++)
			  {
               ij=INDS(i,j);
               imagepot[ij]=- MAX(0.0,bd.sioxbgo-0.250*bd.beta*bd.beta*
                                (1.0/(gridx[i]-gridx[ibreg[ireg]]) +
                                1.0/(-gridx[i]+gridx[iereg[ireg]])));
			   }//ENDDO
            }//ENDDO
         }//ENDIF
      }//ENDDO

//____end of GETIMAGEPOT
      return;
}
//=====

void DevSimulator::SETDT(void)
{
//    Purpose: set propagation time for each particle
//____local variables
      int ipar;

//____set time of time slice length for every particle
      for(ipar=0;ipar<npar0;ipar++)
	  {
         dfield[6][ipar]=dt;//defield(7,0)
      }

//____end of SETDT
      return;
}
//=====
bool DevSimulator::CONVER(int idt)
{
//    int idt
//    Purpose: Stopps simulation when results are converged within
//    +/- stopii with a probabillity of 95%.
//____local variables
      bool stopfl;
      double errest,mvexp,svexp;

      stopfl=false;

      if(idt>10)
	  {
         if(holeghdmfl)
		 {
            mvexp=meanghexp[GHV][PHOLE]/(double)(idt);
            svexp=squghexp[GHV][PHOLE]/(double)(idt);
            if(mvexp!=0.0)
			{
                errest=sqrt((svexp-mvexp*mvexp)/(double)(idt-1))/fabs(mvexp)*2.0;
               if(errest<facstopii)stopfl=true;
            }//ENDIF
         }//ENDIF

         if(elecghdmfl)
		 {
            mvexp=meanghexp[GHV][PELEC]/(double)(idt);
            svexp=squghexp[GHV][PELEC]/(double)(idt);
            if(mvexp!=0.0)
			{
                errest=sqrt((svexp-mvexp*mvexp)/(double)(idt-1))/fabs(mvexp)*2.0;
               if(errest < facstopii)stopfl=true;
            }//ENDIF
         }//ENDIF

         if(elecisubfl)
		 {
            if(meaniicurjb[PELEC]!=0.0)
			{
                errest=sqrt((squiicurjb[PELEC]-meaniicurjb[PELEC]*meaniicurjb[PELEC]/(double)(idt))
					 /((double)idt*(double)(idt-1))) 
                  /fabs(meaniicurjb[PELEC]/(double)idt/2.0); 
               if(errest < facstopii) stopfl=true;
            }//ENDIF
         }//ENDIF

         if(holeisubfl)
		 {
            if(meaniicurjb[PHOLE]!=0.0)
			{
                errest=sqrt((squiicurjb[PHOLE]-meaniicurjb[PHOLE]*meaniicurjb[PHOLE]/(double)(idt)) 
                    /((double)(idt)*(double)(idt-1))) 
                  /fabs(meaniicurjb[PHOLE]/(double)(idt)/2.0); 
               if(errest < facstopii) stopfl=true;
            }//ENDIF
         }//ENDIF

      }//ENDIF

//____end of CONVER
      return stopfl;
}
//------------------------------------------------------------------
void DevSimulator::ELEC(int idt,bool wrfl,bool averfl,Partical &par,Band &bd)
{
//    Purpose: - loops over all particles for current observation time
//    intervall
//____local variables
      bool catchfl,phfl,quadfl,tetfl,impfl;
      int ipar,i,j,jj,ireal,icpu,ierr,iptype,icont,jcont;
      int nparold[NPARTYP],nparold0;
      double rtotal[NPARTYP],tcpu;
      double tetruef,quadtf,phtf,imptf,tf;
      double xfs,yfs,phrnl,imprnl;
      double gamimpmax;
      double discur[MNCONT];
	  int nndt;
	  
	  ofstream ftp("output/sca.txt",ios::app);
//____save number of particles
      for(iptype=0;iptype<NPARTYP;iptype++)
	  {
         nparold[iptype]=npar[iptype];
      }
	  nparold0=npar0;
//    contact current
      for(icont=0;icont<ncont;icont++)
	  {
         curcont[icont]=0.0;
		 for(jj=0;jj<3;jj++)	curcontt[jj][icont]=0;
		 dgen[icont]=0.0;
		 dcatch[icont]=0.0;
		 for(jj=0;jj<3;jj++)	dgenn[jj][icont]=0;
		 for(jj=0;jj<3;jj++)	dcatchh[jj][icont]=0;
      }
	  
	  ipar=0;
	  while(ipar<npar0)
	  {
		  par.singleele(ipar,catchfl,averfl,*this,bd);
		  if(catchfl)
		     catchfl=false;
		  else
		     ipar=ipar+1;
	  }

//____calculate band-population
	  if (wrfl)
	  {
         for(iptype=0;iptype<NPARTYP;iptype++)
		 {
            rtotal[iptype]=0.0;
            for(i=bd.bandof[iptype];i<bd.nband[iptype]+bd.bandof[iptype];i++)
			{
               rtotal[iptype]=rtotal[iptype]+par.dtotp[i];
            }
            if (rtotal[iptype]<=0.0) rtotal[iptype]=1.0;
         }
      }
      
//____calculate displacement current
      for(icont=0;icont<ncont;icont++)
	  {
         discur[icont]=0.0;
         for(jcont=0;jcont<ncont;jcont++)
		 {
            if (cdt>=2)
			{
               discur[icont]=discur[icont]+cap[icont][jcont]
               *(ctpotnow[jcont][cdt]-ctpotnow[jcont][cdt-1])/dt;
            }
         }
      }

//____statistic for R-S (every time steps)
      for(icont=0;icont<ncont;icont++)
	  {
		 currentarray[icont][idt]=curcont[icont];
		 for(jj=0;jj<3;jj++)	currentarrayy[jj][icont][idt]=curcontt[jj][icont];
		 if(cdt>=starttime)
		 {
			meancurcont[icont]=meancurcont[icont]
				               +curcont[icont] - discur[icont];
			dcurrent[icont]=dcurrent[icont]+dgen[icont]-dcatch[icont];
			for(jj=0;jj<3;jj++)
			{
				meancurcontt[jj][icont]=meancurcontt[jj][icont]
								    +curcontt[jj][icont] - discur[icont];
				dcurrentt[jj][icont]=dcurrentt[jj][icont]+dgenn[jj][icont]-dcatchh[jj][icont];
			}
			squcurcont[icont]=squcurcont[icont] 
				              +(curcont[icont] - discur[icont])*(curcont[icont] - discur[icont]);
		 }
      }

//____write scattering event table for electrons
      if (wrfl)
	  {
		 nndt=idt+1;
		 ftp<<nndt<<endl;
		 ftp<<endl;
         if (npar[PELEC]>0)
		 {
            if (sifl)
			{
				ftp.width(10);
				ftp<<"Number of electrons (before/after propagation) :"<<nparold[PELEC]<<"  "<<npar[PELEC]<<endl;
				ftp<<" S c a t t e r i n g   E v e n t   T a b l e (Electrons)"<<endl;
				ftp<<endl;
				ftp<<"             band 1    band 2    band 3    band 4"<<endl;
				ftp<<endl;
				ftp<<"   tet    :  ";
				for(i=0;i<bd.nband[PELEC];i++)ftp<<par.ntet[i+bd.bandof[PELEC]]<<"  ";
				ftp<<endl;
				ftp<<"   bz     :  ";
				for(i=0;i<bd.nband[PELEC];i++)ftp<<par.nbz[i+bd.bandof[PELEC]]<<"  ";
				ftp<<endl;
				ftp<<"   quad   :  ";
				for(i=0;i<bd.nband[PELEC];i++)ftp<<par.nquad[i+bd.bandof[PELEC]]<<"  ";
				ftp<<endl;
				ftp<<"   tot    :  ";
				for(i=0;i<bd.nband[PELEC];i++)ftp<<par.ntotp[i+bd.bandof[PELEC]]<<"  ";
				ftp<<endl;
				ftp<<"   self   :  ";
				for(i=0;i<bd.nband[PELEC];i++)ftp<<par.nslfp[i+bd.bandof[PELEC]]<<"  ";
				ftp<<endl;
				ftp<<"   real ep:  ";
				for(i=0;i<bd.nband[PELEC];i++)ftp<<par.nreap[i+bd.bandof[PELEC]]<<"  ";
				ftp<<endl;
				ftp<<"   totbh  :  ";
				for(i=0;i<bd.nband[PELEC];i++)ftp<<par.ntotbh[i+bd.bandof[PELEC]]<<"  ";
				ftp<<endl;
				ftp<<"   self BH:  ";
				for(i=0;i<bd.nband[PELEC];i++)ftp<<par.nslfbh[i+bd.bandof[PELEC]]<<"  ";
				ftp<<endl;
				ftp<<"   real BH:  ";
				for(i=0;i<bd.nband[PELEC];i++)ftp<<par.nreabh[i+bd.bandof[PELEC]]<<"  ";
				ftp<<endl;
				ftp<<"   real II:  ";
				for(i=0;i<bd.nband[PELEC];i++)ftp<<par.nreaii[i+bd.bandof[PELEC]]<<"  ";
				ftp<<endl;
				ftp<<endl;
				ftp<<endl;
				for(i=0;i<bd.nband[PELEC];i++)
				{
					if(bd.jacophfl)
					{
						ftp<<"   ac     :";
						for(j=0;j<bd.nband[PELEC];j++)ftp<<par.nsctype[i*bd.scpre+0][j]<<"  ";
						ftp<<endl;
						ftp<<"   g-TA-ab:";
						for(j=0;j<bd.nband[PELEC];j++)ftp<<par.nsctype[i*bd.scpre+1][j]<<"  ";
						ftp<<endl;
						ftp<<"   g-TA-em:";
						for(j=0;j<bd.nband[PELEC];j++)ftp<<par.nsctype[i*bd.scpre+2][j]<<"  ";
						ftp<<endl;
						ftp<<"   g-LA-ab:";
						for(j=0;j<bd.nband[PELEC];j++)ftp<<par.nsctype[i*bd.scpre+3][j]<<"  ";
						ftp<<endl;
						ftp<<"   g-LA-em:";
						for(j=0;j<bd.nband[PELEC];j++)ftp<<par.nsctype[i*bd.scpre+4][j]<<"  ";
						ftp<<endl;
						ftp<<"   g-LO-ab:";
						for(j=0;j<bd.nband[PELEC];j++)ftp<<par.nsctype[i*bd.scpre+5][j]<<"  ";
						ftp<<endl;
						ftp<<"   g-LO-em:";
						for(j=0;j<bd.nband[PELEC];j++)ftp<<par.nsctype[i*bd.scpre+6][j]<<"  ";
						ftp<<endl;
						ftp<<"   f-TA-ab:";
						for(j=0;j<bd.nband[PELEC];j++)ftp<<par.nsctype[i*bd.scpre+7][j]<<"  ";
						ftp<<endl;
						ftp<<"   f-TA-em:";
						for(j=0;j<bd.nband[PELEC];j++)ftp<<par.nsctype[i*bd.scpre+8][j]<<"  ";
						ftp<<endl;
						ftp<<"   f-LA-ab:";
						for(j=0;j<bd.nband[PELEC];j++)ftp<<par.nsctype[i*bd.scpre+9][j]<<"  ";
						ftp<<endl;
						ftp<<"   f-LA-em:";
						for(j=0;j<bd.nband[PELEC];j++)ftp<<par.nsctype[i*bd.scpre+10][j]<<"  ";
						ftp<<endl;
						ftp<<"   f-TO-ab:";
						for(j=0;j<bd.nband[PELEC];j++)ftp<<par.nsctype[i*bd.scpre+11][j]<<"  ";
						ftp<<endl;
						ftp<<"   f-TO-em:";
						for(j=0;j<bd.nband[PELEC];j++)ftp<<par.nsctype[i*bd.scpre+12][j]<<"  ";
						ftp<<endl;
						ftp<<"   II     :";
						for(j=0;j<bd.nband[PELEC];j++)ftp<<par.nsctype[i*bd.scpre+13][j]<<"  ";
						ftp<<endl;
						ftp<<endl;
					}
					else if(bd.fiscphfl)
					{
						;
					}
				}
			}
		 }
		 
         if (npar[PHOLE]>0)
		 {
            if (sifl)
			{
				ftp.width(10);
				ftp<<"Number of holes (before/after propagation) :"<<nparold[PHOLE]<<"  "<<npar[PHOLE]<<endl;
				ftp<<" S c a t t e r i n g   E v e n t   T a b l e (Holes)"<<endl;
				ftp<<endl;
				ftp<<"             band 1    band 2    band 3"<<endl;
				ftp<<endl;
				ftp<<"   tet    :  ";
				for(i=0;i<bd.nband[PHOLE];i++)ftp<<par.ntet[i+bd.bandof[PHOLE]]<<"  ";
				ftp<<endl;
				ftp<<"   bz     :  ";
				for(i=0;i<bd.nband[PHOLE];i++)ftp<<par.nbz[i+bd.bandof[PHOLE]]<<"  ";
				ftp<<endl;
				ftp<<"   quad   :  ";
				for(i=0;i<bd.nband[PHOLE];i++)ftp<<par.nquad[i+bd.bandof[PHOLE]]<<"  ";
				ftp<<endl;
				ftp<<"   tot    :  ";
				for(i=0;i<bd.nband[PHOLE];i++)ftp<<par.ntotp[i+bd.bandof[PHOLE]]<<"  ";
				ftp<<endl;
				ftp<<"   self   :  ";
				for(i=0;i<bd.nband[PHOLE];i++)ftp<<par.nslfp[i+bd.bandof[PHOLE]]<<"  ";
				ftp<<endl;
				ftp<<"   real ep:  ";
				for(i=0;i<bd.nband[PHOLE];i++)ftp<<par.nreap[i+bd.bandof[PHOLE]]<<"  ";
				ftp<<endl;
				ftp<<"   totbh  :  ";
				for(i=0;i<bd.nband[PHOLE];i++)ftp<<par.ntotbh[i+bd.bandof[PHOLE]]<<"  ";
				ftp<<endl;
				ftp<<"   self BH:  ";
				for(i=0;i<bd.nband[PHOLE];i++)ftp<<par.nslfbh[i+bd.bandof[PHOLE]]<<"  ";
				ftp<<endl;
				ftp<<"   real BH:  ";
				for(i=0;i<bd.nband[PHOLE];i++)ftp<<par.nreabh[i+bd.bandof[PHOLE]]<<"  ";
				ftp<<endl;
				ftp<<"   real II:  ";
				for(i=0;i<bd.nband[PHOLE];i++)ftp<<par.nreaii[i+bd.bandof[PHOLE]]<<"  ";
				ftp<<endl;
				ftp<<endl;
				ftp<<endl;
				for(i=0;i<bd.nband[PHOLE];i++)
				{
					if(bd.jacophfl)
					{
						ftp<<"   ac     :";
						for(j=0;j<NBH;j++)ftp<<par.nsctyph[i*bd.scprh+0][j]<<"  ";
						ftp<<endl;
						ftp<<"   opt  ab:";
						for(j=0;j<NBH;j++)ftp<<par.nsctyph[i*bd.scprh+1][j]<<"  ";
						ftp<<endl;
						ftp<<"   opt  ac:";
						for(j=0;j<NBH;j++)ftp<<par.nsctyph[i*bd.scprh+2][j]<<"  ";
						ftp<<endl;
						ftp<<endl;
					}
					else if(bd.fiscphfl)
					{
						;
					}
				}
			}
		 }
	}
//____end of ELEC
      return;
}
//=====

void DevSimulator::EINIT(Partical &par)
{
//    Purpose: - initialization of electron related values for each run
//____local variables
      int i,j,icont;
	  par.EINIT();
//____set cumulative variables to zero

//____contact statistics
      for(icont=0;icont<ncont;icont++)
	  {
         ngen[icont]=0;
         dgen[icont]=0.0;
         ncatch[icont]=0;
         dcatch[icont]=0.0;
		 dcurrent[icont]=0.0;
		 for(i=0;i<3;i++)	dgenn[i][icont]=0.0;	
		 for(i=0;i<3;i++)	dcatchh[i][icont]=0.0;
		 for(i=0;i<3;i++)	dcurrentt[i][icont]=0.0;
      }
//____end of EINIT
      return;
}

void DevSimulator::ZEROJB(void)
{
//    Purpose: initalize just before scattering statistics
//____local variables
      int iptype,itab,ighexp;

//____initialize just before values
      for(iptype=0;iptype<NPARTYP;iptype++)
	  {
//       impact ionization
         iicurjb[iptype]=0.0;
         iicurts[iptype]=0.0;
//       velocity,energy,...
         for(ighexp=0;ighexp<NGHEXP;ighexp++)
		 {
            ghexp[ighexp][iptype]=0.0;
         }
//       homogeneous distribution functions
         for(itab=0;itab<MTAB;itab++)
		 {
            edfjb[itab][iptype]=0.0;
            edishom[itab][iptype]=0.0;
         }
      }
//____end of ZEROJB
      return;
}

void DevSimulator::ADDJB(int stridejb,Band &bd)
{
//    Purpose: add just before scattering statistics
//____local variables
      int iptype,itab,ighexp;
//____update mean and square of the sampled data
      for(iptype=0;iptype<NPARTYP;iptype++)
	  {
         iicurjb[iptype]=iicurjb[iptype]/dt/(double)(stridejb);
         meaniicurjb[iptype]=meaniicurjb[iptype]+iicurjb[iptype];
         squiicurjb[iptype]=squiicurjb[iptype]+iicurjb[iptype]*iicurjb[iptype];
         meaniicurts[iptype]=meaniicurts[iptype]+iicurts[iptype];
         squiicurts[iptype]=squiicurts[iptype]+iicurts[iptype]*iicurts[iptype];

         for(ighexp=0;ighexp<NGHEXP;ighexp++)
		 {
            ghexp[ighexp][iptype]=ghexp[ighexp][iptype] 
                                /dt/(double)(stridejb);
            meanghexp[ighexp][iptype]=meanghexp[ighexp][iptype] 
                                        +ghexp[ighexp][iptype];
            squghexp[ighexp][iptype]=squghexp[ighexp][iptype] 
                                      +ghexp[ighexp][iptype]*ghexp[ighexp][iptype];
         }

         for(itab=0;itab<MTAB;itab++)
		 {
            edfjb[itab][iptype]=edfjb[itab][iptype]/dt/bd.dtable
                                                   /(double)(stridejb);
            meanedfjb[itab][iptype]=meanedfjb[itab][iptype] 
                                    +edfjb[itab][iptype] ;
            squedfjb[itab][iptype]=squedfjb[itab][iptype] 
                                   +edfjb[itab][iptype]*edfjb[itab][iptype];
            edishom[itab][iptype]=edishom[itab][iptype]/dt/bd.dtable
                                                   /(double)(stridejb);
            meanedishom[itab][iptype]=meanedishom[itab][iptype] 
                                    +edishom[itab][iptype];
            squedishom[itab][iptype]=squedishom[itab][iptype] 
                                   +edishom[itab][iptype]*edishom[itab][iptype];
         }
      }
//____end of ADDJB
	  return;
}
void DevSimulator::REFRESH(bool wrfl,Band &bd,Partical &par)
{

//    Purpose: perform refresh
//____local variables 
      int ipar, iptype, ireal, icpu, ierr;
	  static int i, j, ij, ier, nempty[NPARTYP], nrefresh[NPARTYP];
	  static int nparold[NPARTYP],ifree,nsum,ipsc,nparold0;
      double tcpu,pc;
	  int ipt,ijqp;
	  int carriertype;

//____clean particle memory (eliminate all particles with a very small charge)
      for(iptype=0;iptype<NPARTYP;iptype++)	nparold[iptype]=npar[iptype];
	  nparold0=npar0;
      ipar=0;
      while(ipar<npar0)
	  {
         if(fabs(dfield[5][ipar])<pckill*pcaver)
		 {
            while((fabs(dfield[5][npar0-1])<pckill*pcaver)&&(npar0>0))
			{
               iptype = bd.partyp[bd.ibt[ifield[0][npar0-1]]];
               npar[iptype]=npar[iptype]-1;
               npar0=npar0-1;
            }
            if(ipar<(npar0-1))
			{
			   iptype=bd.partyp[bd.ibt[ifield[0][ipar]]];
               npar[iptype]=npar[iptype] - 1;
			   ifield[IVPP][ipar]=ifield[IVPP][npar0-1];
               ifield[1][ipar] = ifield[1][npar0-1];
               ifield[2][ipar] = ifield[2][npar0-1];
               ifield[0][ipar] = ifield[0][npar0-1];
               dfield[1][ipar] = dfield[1][npar0-1];
               dfield[2][ipar] = dfield[2][npar0-1];
               dfield[3][ipar] = dfield[3][npar0-1];
               dfield[4][ipar] = dfield[4][npar0-1];
               dfield[5][ipar] = dfield[5][npar0-1];
               dfield[0][ipar] = dfield[0][npar0-1];
               dfield[7][ipar] = dfield[7][npar0-1];
               npar0 = npar0 - 1;
            }
         }
         ipar=ipar+1;
      }
//____initialize list for first particle of chain of particles of quadrant ij
//    and other local lists and counters
      for(ipsc=0;ipsc<npsc;ipsc++)
	  {
         ifirst[ipsc]=-1;
         npper[ipsc]=0;
         sumpc[ipsc]=0.0;
         sumquad[ipsc]=0.0;
	  }
      for(iptype=0;iptype<NPARTYP;iptype++)
	  {
         nempty[iptype]=0;
         nrefresh[iptype]=0;
      }
      ifree=-1;

//____loop over all particles
      for(ipar=0;ipar<npar0;ipar++)
	  {
         ipt=bd.partyp[bd.ibt[ifield[0][ipar]]];
         if(mulfl[ipt]||enermulfl[ipt])
		 {
//          regions also in real space
            if(mulfl[ipt])ijqp=ifield[2][ipar];
//          real space treated as one region
            if(enermulfl[ipt])ijqp=MNGP;/////////////////be carful!!!!!!!!!!!!!
//          region index in energy
            ier=(int)(dfield[7][ipar]/erfmax[ijqp][ipt]*nerrf[ijqp][ipt]);      
            if(ier>=nerrf[ijqp][ipt])ier=nerrf[ijqp][ipt];
            ipsc=ierp[ier][ijqp][ipt];
            if(ningfl&&(ier==0)&&(ifield[2][ipar]==ijning))
			{
//             remove particles in the lowest region of Ning experiment
               dfield[5][ipar] = 0.0;
               npar[ipt]=npar[ipt]-1;
			}
            else if(ipsc>=0)
			{
//             make list of particles in region ipsc
               inext[ipar]=ifirst[ipsc];
               ifirst[ipsc]=ipar;
               npper[ipsc]= npper[ipsc]+1;
               pc= dfield[5][ipar];
               sumpc[ipsc]=sumpc[ipsc]+fabs(pc);
               sumquad[ipsc]=sumquad[ipsc]+pc*pc;
			}
         }
      }//end for

//____perform refresh
//____loop over all particle types
      for(iptype=0;iptype<NPARTYP;iptype++)
	  {
//_______energy refresh
         if(enermulfl[iptype])
		 {
            ij=MNGP;
            REPLACE(ij,iptype,ifree,nempty,nrefresh,par);
         }
//_______phase space refresh
         if(mulfl[iptype])
		 {
//_______loop over all quadrants
            for(i=0;i<ngpx-1;i++)
			for(j=0;j<ngpy-1;j++)
			{
               ij=INDS(i,j);
               REPLACE(ij,iptype,ifree,nempty,nrefresh,par);
			}
		 }
	  }

//____clean particle memory (eliminate all particles with zero charge)
      ipar=0;
      while(ipar<npar0)
	  {
         if(dfield[5][ipar]==0.0)
		 {
            while((dfield[5][npar0-1]==0.0)&&(npar0>0))
				npar0=npar0-1;
            if(ipar<npar0-1)
			{
			   ifield[IVPP][ipar]=ifield[IVPP][npar0-1];
               ifield[1][ipar] = ifield[1][npar0-1];
               ifield[2][ipar] = ifield[2][npar0-1];
               ifield[0][ipar] = ifield[0][npar0-1];
               dfield[1][ipar] = dfield[1][npar0-1];
               dfield[2][ipar] = dfield[2][npar0-1];
               dfield[3][ipar] = dfield[3][npar0-1];
               dfield[4][ipar] = dfield[4][npar0-1];
               dfield[5][ipar] = dfield[5][npar0-1];
               dfield[0][ipar] = dfield[0][npar0-1];
               dfield[7][ipar] = dfield[7][npar0-1];
               npar0=npar0-1;
			}
         }
         ipar=ipar+1;
      }
//____calculate particle number
      nsum=0;
      for(iptype=0;iptype<NPARTYP;iptype++)
         nsum=nsum+npar[iptype];
	  if(nsum!=npar0)
	  {
		  cout<<"error!"<<endl;
		  exit(0);
	  }
//____check refresh for loss of particle charge
	  for(ipsc=0;ipsc<npsc;ipsc++)
	  {
		  sumquad[ipsc]=0.0;
	  }
            
//____loop over all particles
         for(ipar=0;ipar<npar0;ipar++)
		 {
            ipt=bd.partyp[bd.ibt[ifield[0][ipar]]];
            if(mulfl[ipt]||enermulfl[ipt])
			{
			  ijqp = ifield[2][ipar];
               if(enermulfl[ipt]) ijqp =MNGP;
               ier = (int)(dfield[7][ipar] / erfmax[ijqp][ipt]
                                        * double(nerrf[ijqp][ipt]));// + 1
               if(ier>nerrf[ijqp][ipt])ier=nerrf[ijqp][ipt];// + 1
               ipsc= ierp[ier][ijqp][ipt];
               if(ipsc>=0)
			   {
				   pc=fabs(dfield[5][ipar]);
				   sumquad[ipsc]+=pc;//fabs(dfield[5][ipar]);
			   }       
            }
         }
         for(ipsc=0;ipsc<npsc;ipsc++)
		 {
            if(fabs(sumpc[ipsc]-sumquad[ipsc])>1e-10*sumpc[ipsc])
			{
				cout<<"REFRESH: Loss of particle charge"<<endl;
				cout<<sumpc[ipsc]<<endl;
				cout<<sumpc[ipsc]<<endl;
				cout<<ipsc<<endl;
            }
         }        
//____end of REFRESH
    return;
}
//===

void DevSimulator::REPLACE(int ij,int iptype,int &ifree,
					 int *nempty,int *nrefresh,Partical &par)
{
//    Purpose: perform actual refresh
//____local variables
      const int MNRP=10000;
      int ipar;
	  int il, ir, ier;
	  static int ilist[MNRP],ibuf[IVPP+1][MNRP];
	  static int iflist[MNPAR],ipsc;
	  static double pcbuf[MNRP], dbuf[DVPP][MNRP];
      double rpc,pcsum;
	  double pc;
	  int carriertype;

	for(ier=0;ier<=nerrf[ij][iptype];ier++)
	{
      ipsc=ierp[ier][ij][iptype];
	  if(ipsc!=-2)
	  {
	     if (ndpper[ipsc]>0)
		 {
//          refresh only possible, if atleast one particle is within quadrant
			if(npper[ipsc]==0)
			{
				nempty[iptype]=nempty[iptype]+1;
			}
//          check if refresh is necessary
//          particle cirterium
			else if ((npper[ipsc]>ndpper[ipsc]*rpar)
					||(npper[ipsc]<ndpper[ipsc]/rpar)
//                  square of particle weight cirterium
					||(sumquad[ipsc]*npper[ipsc]>sumpc[ipsc]*sumpc[ipsc]*rquad))
			{

				nrefresh[iptype]=nrefresh[iptype]+1;
				nrfper[ipsc]=nrfper[ipsc]+1;
//__________make list of particles in chain and add them to list of free
//          particle cells
				il=-1;
				pcsum=0.0;
				ipar=ifirst[ipsc];
				while(ipar!=-1)
				{
				   il=il+1;
	               ilist[il]=ipar;
				   pc=fabs(dfield[5][ipar]);
		           pcsum=pcsum+pc;//+fabs(dfield[5][ipar]);
			       pcbuf[il]=pcsum;
	
		           iflist[ipar]= ifree;
			       ifree=ipar;
				   dfield[5][ipar]=0.0;

	               npar[iptype]=npar[iptype]-1;
		           ipar=inext[ipar];
				}
              
//__________check number of particles in actual cell against value 
//          determined above

//__________chose new particle states from the old ones
	            for(ir=0;ir<ndpper[ipsc];ir++)
				{
			       rpc=pcsum * RANDNR();
				   il=0;
	               while(rpc>pcbuf[il])
				   {
			          il=il+1;
				   }
	               ipar=ilist[il];
				   ibuf[IVPP][ir]=ifield[IVPP][ipar];
		           ibuf[1][ir] = ifield[1][ipar];
			       ibuf[2][ir] = ifield[2][ipar];
				   ibuf[0][ir] = ifield[0][ipar];
	               dbuf[1][ir] = dfield[1][ipar];
		           dbuf[2][ir] = dfield[2][ipar];
			       dbuf[3][ir] = dfield[3][ipar];
				   dbuf[4][ir] = dfield[4][ipar];
	               dbuf[0][ir] = dfield[0][ipar];
		           dbuf[7][ir] = dfield[7][ipar];
			    }

//__________calculate new particle weight
		double pc=par.charsign[iptype]*pcsum/ndpper[ipsc];

//__________copy new particles from buffer into particle memory
	            for(ir=0;ir<ndpper[ipsc];ir++)
				{
//                 push particles on free cell heap
	               if(!(ifree==-1))
				   {
			          ipar=ifree;
				      ifree=iflist[ipar];
					  ifield[IVPP][ipar]=ibuf[IVPP][ir];
				      ifield[1][ipar] = ibuf[1][ir];
					  ifield[2][ipar] = ibuf[2][ir];
	                  ifield[0][ipar] = ibuf[0][ir];
		              dfield[1][ipar] = dbuf[1][ir];
			          dfield[2][ipar] = dbuf[2][ir];
				      dfield[3][ipar] = dbuf[3][ir];
					  dfield[4][ipar] = dbuf[4][ir];
	                  dfield[0][ipar] = dbuf[0][ir];
		              dfield[5][ipar] = pc;
			          dfield[7][ipar] = dbuf[7][ir];
				      npar[iptype] = npar[iptype]+1;
				   }
//                 push particles in particle memory behind last particle
		           else
				   {
					  npar[iptype] = npar[iptype]+1;
					  npar0=npar0+1;
					  ifield[IVPP][npar0-1]=ibuf[IVPP][ir];
	                  ifield[1][npar0-1] = ibuf[1][ir];
		              ifield[2][npar0-1] = ibuf[2][ir];
			          ifield[0][npar0-1] = ibuf[0][ir];
				      dfield[1][npar0-1] = dbuf[1][ir];
					  dfield[2][npar0-1] = dbuf[2][ir];
	                  dfield[3][npar0-1] = dbuf[3][ir];
		              dfield[4][npar0-1] = dbuf[4][ir];
			          dfield[0][npar0-1] = dbuf[0][ir];
				      dfield[5][npar0-1] = pc;
					  dfield[7][npar0-1] = dbuf[7][ir];
	               }
	            }//end for
		     }       
		}       
      }       
    }//end for

//____end of REPLACE
	return;
}
void DevSimulator::AVERINIT(void)
  {
//____local variables
      int istat,ij,iptype,ireal,icpu,ierr,ier,i,j,jj,ipsc;
      int itab,icont,ighexp;
      double tcpu;
//____initialize statistics
      for(iptype=0;iptype<NPARTYP;iptype++)
	  {
         for(ij=0;ij<ngp;ij++)
		 {
            for(istat=0;istat<NSTAT;istat++)
			{
               statis[istat][ij][iptype]=0.0;
            }
            for(ier=0;ier<NEREDF;ier++)
			{
               edf[ier][ij][iptype]=0.0;
            }
            curx[ij][iptype]=0.0;
            cury[ij][iptype]=0.0;
         }
      }
      for(ij=0;ij<ngp;ij++)
	  {
         statpot[ij]=0.0;
      }
      for(iptype=0;iptype<NPARTYP;iptype++)
	  {
         for(i=0;i<ngpx;i++)
		 {
            for(ier=0;ier<NEREDF;ier++)
			{
               bedfx[ier][i][iptype]=0.0;
            }
         }
         for( j=0;j<ngpy;j++)
		 {
            for(ier=0;ier<NEREDF;ier++)
			{
               bedfy[ier][j][iptype]=0.0;
            }
         }
      }
      if (mulfl0)
	  {
         for(ipsc=0;ipsc<MNPSC;ipsc++)
		 {
            nrfper[ipsc]=0;
         }
      }

//____initialize just before values for impact ionization current
      for(icont=0;icont<ncont;icont++)
	  {
         meancurcont[icont]=0.0;
         squcurcont[icont] =0.0;
		 for(jj=0;jj<3;jj++)	meancurcontt[jj][icont]=0;
      }
	  return;
  }
void DevSimulator::output(void)
{
	int i,j;
	ofstream ftp("output/gridf1.txt");
	ofstream dftp("output/dens.txt");//	cm-3
	ofstream daftp("output/densaver.txt");//	cm-3
	ofstream dpftp("output/doping.txt");//  cm-3
	ofstream exftp("output/xfield.txt");//	V/cm
	ofstream eyftp("output/yfield.txt");//	V/cm
	ofstream spftp("output/statpot.txt");//	V
	ofstream pftp("output/pot.txt");//		V
	ofstream vftp("output/vxyout.txt");//	m/s
	ofstream vvftp("output/vxyout_btbt.txt");//	m/s
	ofstream eeftp("output/energyout.txt");//   eV
	ofstream eeeftp("output/energyout_btbt.txt");//   eV
	ofstream ebandftp("output/ebandratio.txt");//  percent
	ofstream hbandftp("output/hbandratio.txt");//  percent

	ftp<<ngpx<<endl;
	for(i=0;i<ngpx;i++)
	{
		ftp<<gridx[i]*spr0*1e9<<endl;
	}
	ftp<<ngpy<<endl;
	for(j=0;j<ngpy;j++)	
	{
		ftp<<gridy[j]*spr0*1e9<<endl;
	}

	for(i=0;i<ngpx-1;i++)
	{
		for(j=0;j<ngpy-1;j++)
		{
			dftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9
				<<setw(15)<<quadndens[i+j*ngpx]*conc0/1e6<<setw(15)<<quadpdens[i+j*ngpx]*conc0/1e6
				<<setw(15)<<quadndenss[0][i+j*ngpx]*conc0/1e6<<setw(15)<<quadpdenss[0][i+j*ngpx]*conc0/1e6
				<<setw(15)<<quadndenss[1][i+j*ngpx]*conc0/1e6<<setw(15)<<quadpdenss[1][i+j*ngpx]*conc0/1e6
				<<setw(15)<<quadndenss[2][i+j*ngpx]*conc0/1e6<<setw(15)<<quadpdenss[2][i+j*ngpx]*conc0/1e6<<endl;
			daftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9
				 <<setw(15)<<ndensaver[i+j*ngpx]*conc0/1e6<<setw(15)<<pdensaver[i+j*ngpx]*conc0/1e6
				 <<setw(15)<<ndensaverr[0][i+j*ngpx]*conc0/1e6<<setw(15)<<pdensaverr[0][i+j*ngpx]*conc0/1e6
				 <<setw(15)<<ndensaverr[1][i+j*ngpx]*conc0/1e6<<setw(15)<<pdensaverr[1][i+j*ngpx]*conc0/1e6
				 <<setw(15)<<ndensaverr[2][i+j*ngpx]*conc0/1e6<<setw(15)<<pdensaverr[2][i+j*ngpx]*conc0/1e6<<endl;
			dpftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9<<setw(15)
				 <<donor[0][i+j*ngpx]*conc0/1e6<<setw(15)
				 <<accep[0][i+j*ngpx]*conc0/1e6<<endl;
			spftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9<<setw(15)
				<<statpot[i+j*ngpx]*pot0<<endl;
			pftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9<<setw(15)
				<<pot[i+j*ngpx]*pot0<<endl;
			exftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9<<setw(15)
				<<xfield[i+j*ngpx]*field0/100.0<<endl;
			eyftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9<<setw(15)
				<<yfield[i+j*ngpx]*field0/100.0<<endl;
			vftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9<<setw(15)
				<<evelox[i+j*ngpx]*velo0*100.0<<setw(15)<<eveloy[i+j*ngpx]*velo0<<setw(15)
				<<hvelox[i+j*ngpx]*velo0*100.0<<setw(15)<<hveloy[i+j*ngpx]*velo0<<endl;
			vvftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9<<setw(15)
				<<eveloxx[0][i+j*ngpx]*velo0*100.0<<setw(15)<<eveloyy[0][i+j*ngpx]*velo0<<setw(15)
				<<hveloxx[0][i+j*ngpx]*velo0*100.0<<setw(15)<<hveloyy[0][i+j*ngpx]*velo0<<setw(15)
				<<eveloxx[1][i+j*ngpx]*velo0*100.0<<setw(15)<<eveloyy[1][i+j*ngpx]*velo0<<setw(15)
				<<hveloxx[1][i+j*ngpx]*velo0*100.0<<setw(15)<<hveloyy[1][i+j*ngpx]*velo0<<setw(15)
				<<eveloxx[2][i+j*ngpx]*velo0*100.0<<setw(15)<<eveloyy[2][i+j*ngpx]*velo0<<setw(15)
				<<hveloxx[2][i+j*ngpx]*velo0*100.0<<setw(15)<<hveloyy[2][i+j*ngpx]*velo0<<endl;
			eeftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9<<setw(15)
				 <<eenergy[i+j*ngpx]<<setw(15)<<henergy[i+j*ngpx]<<endl;
			eeeftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9<<setw(15)
				 <<eenergyy[0][i+j*ngpx]<<setw(15)<<henergyy[0][i+j*ngpx]<<setw(15)
				 <<eenergyy[1][i+j*ngpx]<<setw(15)<<henergyy[1][i+j*ngpx]<<setw(15)
				 <<eenergyy[2][i+j*ngpx]<<setw(15)<<henergyy[2][i+j*ngpx]<<endl;
			ebandftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9<<setw(15)
				    <<ebandratio[i+j*ngpx][0]<<setw(15)<<ebandratio[i+j*ngpx][1]<<setw(15)
					<<ebandratio[i+j*ngpx][2]<<setw(15)<<ebandratio[i+j*ngpx][3]<<endl;
			hbandftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9<<setw(15)
				    <<hbandratio[i+j*ngpx][4]<<setw(15)<<hbandratio[i+j*ngpx][5]<<setw(15)
					<<hbandratio[i+j*ngpx][6]<<endl;
		}
	}
	dftp.close();
	daftp.close();
	dpftp.close();
	spftp.close();
	pftp.close();
	exftp.close();
	eyftp.close();
	vftp.close();
	vvftp.close();
	eeftp.close();
	eeeftp.close();
	ebandftp.close();
	hbandftp.close();
	return;
}

void DevSimulator::outputcut(void)
{
	int i,j;
	ofstream dxftp("output/dens_cutx.txt");//	cm-3
	ofstream ndxftp("output/dens_n_cutx.txt");//	cm-3
	ofstream pdxftp("output/dens_p_cutx.txt");//	cm-3
	ofstream ndxaftp("output/densaver_n_cutx.txt");//	cm-3
	ofstream pdxaftp("output/densaver_p_cutx.txt");//	cm-3

	ofstream ndxa0ftp("output/densaver_0_n_cutx.txt");//	cm-3
	ofstream pdxa0ftp("output/densaver_0_p_cutx.txt");//	cm-3
	ofstream ndxa1ftp("output/densaver_1_n_cutx.txt");//	cm-3
	ofstream pdxa1ftp("output/densaver_1_p_cutx.txt");//	cm-3
	ofstream ndxa2ftp("output/densaver_2_n_cutx.txt");//	cm-3
	ofstream pdxa2ftp("output/densaver_2_p_cutx.txt");//	cm-3

	ofstream dpxftp("output/doping_cutx.txt");//	cm-3
	ofstream exxftp("output/field_x_cutx.txt");//	V/cm
	ofstream eyxftp("output/field_y_cutx.txt");//	V/cm
	ofstream spxftp("output/potaver_cutx.txt");//V
	ofstream pxftp("output/pot_cutx.txt");//		V
	ofstream evxxftp("output/vout_x_electron_cutx.txt");//	cm/s
	ofstream evyxftp("output/vout_y_electron_cutx.txt");//	cm/s
	ofstream hvxxftp("output/vout_x_hole_cutx.txt");//	cm/s
	ofstream hvyxftp("output/vout_y_hole_cutx.txt");//	cm/s

	ofstream evxx0ftp("output/vout_x_0_electron_cutx.txt");//	cm/s
	ofstream evxx1ftp("output/vout_x_1_electron_cutx.txt");//	cm/s
	ofstream evxx2ftp("output/vout_x_2_electron_cutx.txt");//	cm/s
	ofstream evyx0ftp("output/vout_y_0_electron_cutx.txt");//	cm/s
	ofstream evyx1ftp("output/vout_y_1_electron_cutx.txt");//	cm/s
	ofstream evyx2ftp("output/vout_y_2_electron_cutx.txt");//	cm/s
	ofstream hvxx0ftp("output/vout_x_0_hole_cutx.txt");//	cm/s
	ofstream hvxx1ftp("output/vout_x_1_hole_cutx.txt");//	cm/s
	ofstream hvxx2ftp("output/vout_x_2_hole_cutx.txt");//	cm/s
	ofstream hvyx0ftp("output/vout_y_0_hole_cutx.txt");//	cm/s
	ofstream hvyx1ftp("output/vout_y_1_hole_cutx.txt");//	cm/s
	ofstream hvyx2ftp("output/vout_y_2_hole_cutx.txt");//	cm/s

	ofstream eexftp("output/energy_electrn_cutx.txt");// eV
	ofstream hexftp("output/energy_hole_cutx.txt");// eV

	ofstream eex0ftp("output/energy_0_electrn_cutx.txt");// eV
	ofstream eex1ftp("output/energy_1_electrn_cutx.txt");// eV
	ofstream eex2ftp("output/energy_2_electrn_cutx.txt");// eV
	ofstream hex0ftp("output/energy_0_hole_cutx.txt");// eV
	ofstream hex1ftp("output/energy_1_hole_cutx.txt");// eV
	ofstream hex2ftp("output/energy_2_hole_cutx.txt");// eV

	ofstream eb0xftp("output/electron_band_0_cutx.txt");
	ofstream eb1xftp("output/electron_band_1_cutx.txt");
	ofstream eb2xftp("output/electron_band_2_cutx.txt");
	ofstream eb3xftp("output/electron_band_3_cutx.txt");
	ofstream hb0xftp("output/hole_band_0_cutx.txt");
	ofstream hb1xftp("output/hole_band_1_cutx.txt");
	ofstream hb2xftp("output/hole_band_2_cutx.txt");


// output xcut
	for(j=0;j<ngpy-1;j++)
	{
		dxftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		ndxftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		pdxftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		ndxaftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		pdxaftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;

		ndxa0ftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		pdxa0ftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		ndxa1ftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		pdxa1ftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		ndxa2ftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		pdxa2ftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;

		dpxftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		exxftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		eyxftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		spxftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		pxftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		evxxftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		evyxftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		hvxxftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		hvyxftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;

		evxx0ftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		evxx1ftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		evxx2ftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		evyx0ftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		evyx1ftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		evyx2ftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		hvxx0ftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		hvxx1ftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		hvxx2ftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		hvyx0ftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		hvyx1ftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		hvyx2ftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;

		eexftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		hexftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;

		eex0ftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		eex1ftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		eex2ftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		hex0ftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		hex1ftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		hex2ftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;

		eb0xftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		eb1xftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		eb2xftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		eb3xftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		hb0xftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		hb1xftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		hb2xftp<<setw(15)<<(gridy[j]+gridy[j+1])/2.0*spr0*1e9;
		for(i=0;i<xcutnumber;i++)
		{
			dxftp<<setw(15)<<quaddens[j*ngpx+ipcutx[i]]*conc0/1e6;
			ndxftp<<setw(15)<<quadndens[j*ngpx+ipcutx[i]]*conc0/1e6;
			pdxftp<<setw(15)<<quadpdens[j*ngpx+ipcutx[i]]*conc0/1e6;
			ndxaftp<<setw(15)<<ndensaver[j*ngpx+ipcutx[i]]*conc0/1e6;
			pdxaftp<<setw(15)<<pdensaver[j*ngpx+ipcutx[i]]*conc0/1e6;

			ndxa0ftp<<setw(15)<<ndensaverr[0][j*ngpx+ipcutx[i]]*conc0/1e6;
			pdxa0ftp<<setw(15)<<pdensaverr[0][j*ngpx+ipcutx[i]]*conc0/1e6;
			ndxa1ftp<<setw(15)<<ndensaverr[1][j*ngpx+ipcutx[i]]*conc0/1e6;
			pdxa1ftp<<setw(15)<<pdensaverr[1][j*ngpx+ipcutx[i]]*conc0/1e6;
			ndxa2ftp<<setw(15)<<ndensaverr[2][j*ngpx+ipcutx[i]]*conc0/1e6;
			pdxa2ftp<<setw(15)<<pdensaverr[2][j*ngpx+ipcutx[i]]*conc0/1e6;

			dpxftp<<setw(15)<<donor[0][j*ngpx+ipcutx[i]]*conc0/1e6<<setw(15)
				  <<accep[0][j*ngpx+ipcutx[i]]*conc0/1e6;
			exxftp<<setw(15)<<xfield[j*ngpx+ipcutx[i]]*field0/100.0;
			eyxftp<<setw(15)<<yfield[j*ngpx+ipcutx[i]]*field0/100.0;
			spxftp<<setw(15)<<statpot[j*ngpx+ipcutx[i]]*pot0;
			pxftp<<setw(15)<<pot[j*ngpx+ipcutx[i]]*pot0;
			evxxftp<<setw(15)<<evelox[j*ngpx+ipcutx[i]]*velo0;
			evyxftp<<setw(15)<<eveloy[j*ngpx+ipcutx[i]]*velo0;
			hvxxftp<<setw(15)<<hvelox[j*ngpx+ipcutx[i]]*velo0;
			hvyxftp<<setw(15)<<hveloy[j*ngpx+ipcutx[i]]*velo0;

			evxx0ftp<<setw(15)<<eveloxx[0][j*ngpx+ipcutx[i]]*velo0;
			evxx1ftp<<setw(15)<<eveloxx[1][j*ngpx+ipcutx[i]]*velo0;
			evxx2ftp<<setw(15)<<eveloxx[2][j*ngpx+ipcutx[i]]*velo0;
			evyx0ftp<<setw(15)<<eveloyy[0][j*ngpx+ipcutx[i]]*velo0;
			evyx1ftp<<setw(15)<<eveloyy[1][j*ngpx+ipcutx[i]]*velo0;
			evyx2ftp<<setw(15)<<eveloyy[2][j*ngpx+ipcutx[i]]*velo0;
			hvxx0ftp<<setw(15)<<hveloxx[0][j*ngpx+ipcutx[i]]*velo0;
			hvxx1ftp<<setw(15)<<hveloxx[1][j*ngpx+ipcutx[i]]*velo0;
			hvxx2ftp<<setw(15)<<hveloxx[2][j*ngpx+ipcutx[i]]*velo0;
			hvyx0ftp<<setw(15)<<hveloyy[0][j*ngpx+ipcutx[i]]*velo0;
			hvyx1ftp<<setw(15)<<hveloyy[1][j*ngpx+ipcutx[i]]*velo0;
			hvyx2ftp<<setw(15)<<hveloyy[2][j*ngpx+ipcutx[i]]*velo0;

			eexftp<<setw(15)<<eenergy[j*ngpx+ipcutx[i]];
			hexftp<<setw(15)<<henergy[j*ngpx+ipcutx[i]];

			eex0ftp<<setw(15)<<eenergyy[0][j*ngpx+ipcutx[i]];
			eex1ftp<<setw(15)<<eenergyy[1][j*ngpx+ipcutx[i]];
			eex2ftp<<setw(15)<<eenergyy[2][j*ngpx+ipcutx[i]];
			hex0ftp<<setw(15)<<henergyy[0][j*ngpx+ipcutx[i]];
			hex1ftp<<setw(15)<<henergyy[1][j*ngpx+ipcutx[i]];
			hex2ftp<<setw(15)<<henergyy[2][j*ngpx+ipcutx[i]];

			eb0xftp<<setw(15)<<ebandratio[j*ngpx+ipcutx[i]][0];
			eb1xftp<<setw(15)<<ebandratio[j*ngpx+ipcutx[i]][1];
			eb2xftp<<setw(15)<<ebandratio[j*ngpx+ipcutx[i]][2];
			eb3xftp<<setw(15)<<ebandratio[j*ngpx+ipcutx[i]][3];
			hb0xftp<<setw(15)<<hbandratio[j*ngpx+ipcutx[i]][4];
			hb1xftp<<setw(15)<<hbandratio[j*ngpx+ipcutx[i]][5];
			hb2xftp<<setw(15)<<hbandratio[j*ngpx+ipcutx[i]][6];
		}
		dxftp<<endl;
		ndxftp<<endl;
		pdxftp<<endl;
		ndxaftp<<endl;
		pdxaftp<<endl;

		ndxa0ftp<<endl;
		pdxa0ftp<<endl;
		ndxa1ftp<<endl;
		pdxa1ftp<<endl;
		ndxa2ftp<<endl;
		pdxa2ftp<<endl;

		dpxftp<<endl;
		exxftp<<endl;
		eyxftp<<endl;
		spxftp<<endl;
		pxftp<<endl;
		evxxftp<<endl;
		evyxftp<<endl;
		hvxxftp<<endl;
		hvyxftp<<endl;

		evxx0ftp<<endl;
		evxx1ftp<<endl;
		evxx2ftp<<endl;
		evyx0ftp<<endl;
		evyx1ftp<<endl;
		evyx2ftp<<endl;
		hvxx0ftp<<endl;
		hvxx1ftp<<endl;
		hvxx2ftp<<endl;
		hvyx0ftp<<endl;
		hvyx1ftp<<endl;
		hvyx2ftp<<endl;

		eexftp<<endl;
		hexftp<<endl;

		eex0ftp<<endl;
		eex1ftp<<endl;
		eex2ftp<<endl;
		hex0ftp<<endl;
		hex1ftp<<endl;
		hex2ftp<<endl;

		eb0xftp<<endl;
		eb1xftp<<endl;
		eb2xftp<<endl;
		eb3xftp<<endl;
		hb0xftp<<endl;
		hb1xftp<<endl;
		hb2xftp<<endl;
	}

		dxftp.close();
		ndxftp.close();
		pdxftp.close();
		ndxaftp.close();
		pdxaftp.close();

		ndxa0ftp.close();
		pdxa0ftp.close();
		ndxa1ftp.close();
		pdxa1ftp.close();
		ndxa2ftp.close();
		pdxa2ftp.close();

		dpxftp.close();
		exxftp.close();
		eyxftp.close();
		spxftp.close();
		pxftp.close();
		evxxftp.close();
		evyxftp.close();
		hvxxftp.close();
		hvyxftp.close();

		evxx0ftp.close();
		evxx1ftp.close();
		evxx2ftp.close();
		evyx0ftp.close();
		evyx1ftp.close();
		evyx2ftp.close();
		hvxx0ftp.close();
		hvxx1ftp.close();
		hvxx2ftp.close();
		hvyx0ftp.close();
		hvyx1ftp.close();
		hvyx2ftp.close();

		eexftp.close();
		hexftp.close();

		eex0ftp.close();
		eex1ftp.close();
		eex2ftp.close();
		hex0ftp.close();
		hex1ftp.close();
		hex2ftp.close();

		eb0xftp.close();
		eb1xftp.close();
		eb2xftp.close();
		eb3xftp.close();
		hb0xftp.close();
		hb1xftp.close();
		hb2xftp.close();

// output ycut
	ofstream dyftp("output/dens_cuty.txt");//	cm-3
	ofstream ndyftp("output/dens_n_cuty.txt");//	cm-3
	ofstream pdyftp("output/dens_p_cuty.txt");//	cm-3
	ofstream ndyaftp("output/densaver_n_cuty.txt");//	cm-3
	ofstream pdyaftp("output/densaver_p_cuty.txt");//	cm-3

	ofstream ndya0ftp("output/densaver_0_n_cuty.txt");//	cm-3
	ofstream pdya0ftp("output/densaver_0_p_cuty.txt");//	cm-3
	ofstream ndya1ftp("output/densaver_1_n_cuty.txt");//	cm-3
	ofstream pdya1ftp("output/densaver_1_p_cuty.txt");//	cm-3
	ofstream ndya2ftp("output/densaver_2_n_cuty.txt");//	cm-3
	ofstream pdya2ftp("output/densaver_2_p_cuty.txt");//	cm-3

	ofstream dpyftp("output/doping_cuty.txt");//	cm-3
	ofstream exyftp("output/field_x_cuty.txt");//	V/cm
	ofstream eyyftp("output/field_y_cuty.txt");//	V/cm
	ofstream spyftp("output/potaver_cuty.txt");//V
	ofstream pyftp("output/pot_cuty.txt");//		V
	ofstream evxyftp("output/vout_x_electron_cuty.txt");//	cm/s
	ofstream evyyftp("output/vout_y_electron_cuty.txt");//	cm/s
	ofstream hvxyftp("output/vout_x_hole_cuty.txt");//	cm/s
	ofstream hvyyftp("output/vout_y_hole_cuty.txt");//	cm/s

	ofstream evxy0ftp("output/vout_x_0_electron_cuty.txt");//	cm/s
	ofstream evxy1ftp("output/vout_x_1_electron_cuty.txt");//	cm/s
	ofstream evxy2ftp("output/vout_x_2_electron_cuty.txt");//	cm/s
	ofstream evyy0ftp("output/vout_y_0_electron_cuty.txt");//	cm/s
	ofstream evyy1ftp("output/vout_y_1_electron_cuty.txt");//	cm/s
	ofstream evyy2ftp("output/vout_y_2_electron_cuty.txt");//	cm/s
	ofstream hvxy0ftp("output/vout_x_0_hole_cuty.txt");//	cm/s
	ofstream hvxy1ftp("output/vout_x_1_hole_cuty.txt");//	cm/s
	ofstream hvxy2ftp("output/vout_x_2_hole_cuty.txt");//	cm/s
	ofstream hvyy0ftp("output/vout_y_0_hole_cuty.txt");//	cm/s
	ofstream hvyy1ftp("output/vout_y_1_hole_cuty.txt");//	cm/s
	ofstream hvyy2ftp("output/vout_y_2_hole_cuty.txt");//	cm/s

	ofstream eeyftp("output/energy_electrn_cuty.txt");// eV
	ofstream heyftp("output/energy_hole_cuty.txt");// eV

	ofstream eey0ftp("output/energy_0_electrn_cuty.txt");// eV
	ofstream eey1ftp("output/energy_1_electrn_cuty.txt");// eV
	ofstream eey2ftp("output/energy_2_electrn_cuty.txt");// eV
	ofstream hey0ftp("output/energy_0_hole_cuty.txt");// eV
	ofstream hey1ftp("output/energy_1_hole_cuty.txt");// eV
	ofstream hey2ftp("output/energy_2_hole_cuty.txt");// eV

	ofstream eb0yftp("output/electron_band_0_cuty.txt");
	ofstream eb1yftp("output/electron_band_1_cuty.txt");
	ofstream eb2yftp("output/electron_band_2_cuty.txt");
	ofstream eb3yftp("output/electron_band_3_cuty.txt");
	ofstream hb0yftp("output/hole_band_0_cuty.txt");
	ofstream hb1yftp("output/hole_band_1_cuty.txt");
	ofstream hb2yftp("output/hole_band_2_cuty.txt");

	for(i=0;i<ngpx-1;i++)
	{
		dyftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		ndyftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		pdyftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		ndyaftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		pdyaftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;

		ndya0ftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		pdya0ftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		ndya1ftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		pdya1ftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		ndya2ftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		pdya2ftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;

		dpyftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		exyftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		eyyftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		spyftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		pyftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		evxyftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		evyyftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		hvxyftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		hvyyftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;

		evxy0ftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		evxy1ftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		evxy2ftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		evyy0ftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		evyy1ftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		evyy2ftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		hvxy0ftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		hvxy1ftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		hvxy2ftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		hvyy0ftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		hvyy1ftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		hvyy2ftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;

		eeyftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		heyftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;

		eey0ftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		eey1ftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		eey2ftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		hey0ftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		hey1ftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		hey2ftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;

		eb0yftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		eb1yftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		eb2yftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		eb3yftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		hb0yftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		hb1yftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		hb2yftp<<setw(15)<<(gridx[i]+gridx[i+1])/2.0*spr0*1e9;
		for(j=0;j<ycutnumber;j++)
		{
			dyftp<<setw(15)<<quaddens[ipcuty[j]*ngpx+i]*conc0/1e6;
			ndyftp<<setw(15)<<quadndens[ipcuty[j]*ngpx+i]*conc0/1e6;
			pdyftp<<setw(15)<<quadpdens[ipcuty[j]*ngpx+i]*conc0/1e6;
			ndyaftp<<setw(15)<<ndensaver[ipcuty[j]*ngpx+i]*conc0/1e6;
			pdyaftp<<setw(15)<<pdensaver[ipcuty[j]*ngpx+i]*conc0/1e6;

			ndya0ftp<<setw(15)<<ndensaverr[0][ipcuty[j]*ngpx+i]*conc0/1e6;
			pdya0ftp<<setw(15)<<pdensaverr[0][ipcuty[j]*ngpx+i]*conc0/1e6;
			ndya1ftp<<setw(15)<<ndensaverr[1][ipcuty[j]*ngpx+i]*conc0/1e6;
			pdya1ftp<<setw(15)<<pdensaverr[1][ipcuty[j]*ngpx+i]*conc0/1e6;
			ndya2ftp<<setw(15)<<ndensaverr[2][ipcuty[j]*ngpx+i]*conc0/1e6;
			pdya2ftp<<setw(15)<<pdensaverr[2][ipcuty[j]*ngpx+i]*conc0/1e6;

			dpyftp<<setw(15)<<donor[0][ipcuty[j]*ngpx+i]*conc0/1e6<<setw(15)
				  <<accep[0][ipcuty[j]*ngpx+i]*conc0/1e6;
			exyftp<<setw(15)<<xfield[ipcuty[j]*ngpx+i]*field0/100.0;
			eyyftp<<setw(15)<<yfield[ipcuty[j]*ngpx+i]*field0/100.0;
			spyftp<<setw(15)<<statpot[ipcuty[j]*ngpx+i]*pot0;
			pyftp<<setw(15)<<pot[ipcuty[j]*ngpx+i]*pot0;
			evxyftp<<setw(15)<<evelox[ipcuty[j]*ngpx+i]*velo0;
			evyyftp<<setw(15)<<eveloy[ipcuty[j]*ngpx+i]*velo0;
			hvxyftp<<setw(15)<<hvelox[ipcuty[j]*ngpx+i]*velo0;
			hvyyftp<<setw(15)<<hveloy[ipcuty[j]*ngpx+i]*velo0;

			evxy0ftp<<setw(15)<<eveloxx[0][ipcuty[j]*ngpx+i]*velo0;
			evxy1ftp<<setw(15)<<eveloxx[1][ipcuty[j]*ngpx+i]*velo0;
			evxy2ftp<<setw(15)<<eveloxx[2][ipcuty[j]*ngpx+i]*velo0;
			evyy0ftp<<setw(15)<<eveloyy[0][ipcuty[j]*ngpx+i]*velo0;
			evyy1ftp<<setw(15)<<eveloyy[1][ipcuty[j]*ngpx+i]*velo0;
			evyy2ftp<<setw(15)<<eveloyy[2][ipcuty[j]*ngpx+i]*velo0;
			hvxy0ftp<<setw(15)<<hveloxx[0][ipcuty[j]*ngpx+i]*velo0;
			hvxy1ftp<<setw(15)<<hveloxx[1][ipcuty[j]*ngpx+i]*velo0;
			hvxy2ftp<<setw(15)<<hveloxx[2][ipcuty[j]*ngpx+i]*velo0;
			hvyy0ftp<<setw(15)<<hveloyy[0][ipcuty[j]*ngpx+i]*velo0;
			hvyy1ftp<<setw(15)<<hveloyy[1][ipcuty[j]*ngpx+i]*velo0;
			hvyy2ftp<<setw(15)<<hveloyy[2][ipcuty[j]*ngpx+i]*velo0;

			eeyftp<<setw(15)<<eenergy[ipcuty[j]*ngpx+i];
			heyftp<<setw(15)<<henergy[ipcuty[j]*ngpx+i];

			eey0ftp<<setw(15)<<eenergyy[0][ipcuty[j]*ngpx+i];
			eey1ftp<<setw(15)<<eenergyy[1][ipcuty[j]*ngpx+i];
			eey2ftp<<setw(15)<<eenergyy[2][ipcuty[j]*ngpx+i];
			hey0ftp<<setw(15)<<henergyy[0][ipcuty[j]*ngpx+i];
			hey1ftp<<setw(15)<<henergyy[1][ipcuty[j]*ngpx+i];
			hey2ftp<<setw(15)<<henergyy[2][ipcuty[j]*ngpx+i];

			eb0yftp<<setw(15)<<ebandratio[ipcuty[j]*ngpx+i][0];
			eb1yftp<<setw(15)<<ebandratio[ipcuty[j]*ngpx+i][1];
			eb2yftp<<setw(15)<<ebandratio[ipcuty[j]*ngpx+i][2];
			eb3yftp<<setw(15)<<ebandratio[ipcuty[j]*ngpx+i][3];
			hb0yftp<<setw(15)<<hbandratio[ipcuty[j]*ngpx+i][4];
			hb1yftp<<setw(15)<<hbandratio[ipcuty[j]*ngpx+i][5];
			hb2yftp<<setw(15)<<hbandratio[ipcuty[j]*ngpx+i][6];
		}
		dyftp<<endl;
		ndyftp<<endl;
		pdyftp<<endl;
		ndyaftp<<endl;
		pdyaftp<<endl;

		ndya0ftp<<endl;
		pdya0ftp<<endl;
		ndya1ftp<<endl;
		pdya1ftp<<endl;
		ndya2ftp<<endl;
		pdya2ftp<<endl;

		dpyftp<<endl;
		exyftp<<endl;
		eyyftp<<endl;
		spyftp<<endl;
		pyftp<<endl;
		evxyftp<<endl;
		evyyftp<<endl;
		hvxyftp<<endl;
		hvyyftp<<endl;

		evxy0ftp<<endl;
		evxy1ftp<<endl;
		evxy2ftp<<endl;
		evyy0ftp<<endl;
		evyy1ftp<<endl;
		evyy2ftp<<endl;
		hvxy0ftp<<endl;
		hvxy1ftp<<endl;
		hvxy2ftp<<endl;
		hvyy0ftp<<endl;
		hvyy1ftp<<endl;
		hvyy2ftp<<endl;

		eeyftp<<endl;
		heyftp<<endl;

		eey0ftp<<endl;
		eey1ftp<<endl;
		eey2ftp<<endl;
		hey0ftp<<endl;
		hey1ftp<<endl;
		hey2ftp<<endl;

		eb0yftp<<endl;
		eb1yftp<<endl;
		eb2yftp<<endl;
		eb3yftp<<endl;
		hb0yftp<<endl;
		hb1yftp<<endl;
		hb2yftp<<endl;
	}
		dyftp.close();
		ndyftp.close();
		pdyftp.close();
		ndyaftp.close();
		pdyaftp.close();

		ndya0ftp.close();
		pdya0ftp.close();
		ndya1ftp.close();
		pdya1ftp.close();
		ndya2ftp.close();
		pdya2ftp.close();

		dpyftp.close();
		exyftp.close();
		eyyftp.close();
		spyftp.close();
		pyftp.close();
		evxyftp.close();
		evyyftp.close();
		hvxyftp.close();
		hvyyftp.close();

		evxy0ftp.close();
		evxy1ftp.close();
		evxy2ftp.close();
		evyy0ftp.close();
		evyy1ftp.close();
		evyy2ftp.close();
		hvxy0ftp.close();
		hvxy1ftp.close();
		hvxy2ftp.close();
		hvyy0ftp.close();
		hvyy1ftp.close();
		hvyy2ftp.close();

		eeyftp.close();
		heyftp.close();

		eey0ftp.close();
		eey1ftp.close();
		eey2ftp.close();
		hey0ftp.close();
		hey1ftp.close();
		hey2ftp.close();

		eb0yftp.close();
		eb1yftp.close();
		eb2yftp.close();
		eb3yftp.close();
		hb0yftp.close();
		hb1yftp.close();
		hb2yftp.close();

return;
}

void DevSimulator::OUTPUTSCATTER(void)
{
//	int i,j;
//	int ijk;
		ofstream ftp;
		ftp.open("output/All_Scattering.txt");
		cout<<"------Impurity scattering =                             "<<bhscatter<<endl;
		ftp<<"------Impurity scattering =                            "<<bhscatter<<endl;
		cout<<"------Phonon scattering =                               "<<phscatter<<endl<<endl;
		ftp<<"------Phonon scattering =                              "<<phscatter<<endl<<endl;

		cout<<"------Old Surface scattering =                          "<<oldsscatter<<endl<<endl;
		ftp<<"------Old Surface scattering =                           "<<oldsscatter<<endl<<endl;

		cout<<"------Surface scattering =                              "<<sscatter<<endl<<endl;
		ftp<<"------Surface scattering =                             "<<sscatter<<endl<<endl;

		cout<<"------Surface roughness scattering =                    "<<srscatter<<endl;
		ftp<<"------Surface roughness scattering =                   "<<srscatter<<endl;
		cout<<"------Surface phonon scattering =                       "<<spscatter<<endl;
		ftp<<"------Surface phonon scattering =                      "<<spscatter<<endl;
		cout<<"------Surface impurity scattering =                     "<<siscatter<<endl;
		ftp<<"------Surface impurity scattering =                    "<<siscatter<<endl;
		cout<<endl<<endl;

		cout<<"------Simulation is completed!!"<<endl<<endl;
		ftp.close();
/*
		ofstream stfm1ftp;
		stfm1ftp.open("output/stfm1output.txt");
		for(i=0;i<ngp;i++)
		{
			for(j=0;j<MNEDF;j++)
			{
				stfm1ftp<<stfm1[0][j][i]<<"  ";
				stfm1ftp<<stfm1[1][j][i]<<"  ";
				stfm1ftp<<stfm1[2][j][i]<<"  ";
				stfm1ftp<<stfm1[3][j][i]<<"  ";
				stfm1ftp<<stfm1[4][j][i]<<"  ";
				stfm1ftp<<stfm1[5][j][i]<<endl;
			}
		}
		stfm1ftp.close();

		ofstream stfm2ftp;
		stfm2ftp.open("output/stfm2output.txt");
		for(i=0;i<ngp;i++)
		{
			for(j=0;j<MNEDF;j++)
			{
				stfm2ftp<<stfm2[0][j][i]<<"  ";
				stfm2ftp<<stfm2[1][j][i]<<"  ";
				stfm2ftp<<stfm2[2][j][i]<<"  ";
				stfm2ftp<<stfm2[3][j][i]<<"  ";
				stfm2ftp<<stfm2[4][j][i]<<"  ";
				stfm2ftp<<stfm2[5][j][i]<<endl;
			}
		}
		stfm2ftp.close();
*/
}