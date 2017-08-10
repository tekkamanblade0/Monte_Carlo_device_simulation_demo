class Partical
{
	 double oldxf,oldyf,oldzf;
	 double oldxv,oldyv,oldzv;
	 double oldxk,oldyk,oldzk;
	 double oldisym;
//   energy band index of the current particle
     int iband;
//   carrier type
	 int ctype;
//   tetrahedron index of the current particle
     int itet;
//   symmetry index of the current particle
     int isym;
//   quadrant index of the current particle
     int ijqp;
//   quadrant index of the current particle in x-direction
     int iqp;
//   quadrant index of the current particle in y-direction
     int jqp;
//   particle type index of the current particle
     int ipt;
//   direction of the boundary hit next in real space by the current particle
     int idir;
//   direction of the boundary hit next in k-space by the current particle
     int itetdir;

//   real space position of the current particle
     double xr,yr;
//   velocity of the current particle
     double xv,yv,zv;
//   k-vector of the current particle
     double xk,yk,zk;
//   electric field at the position of the current particle
     double xf,yf,zf;
//   particle charge of the current particle
     double pc;
//   energy of the current particle
     double ee;
//   remaining simulation time of the current particle
     double dts;
//   charge sign of the different particle types
     double charsign[NPARTYP];
//    energy of the base of the tetrahedron
     double ebzp;
//   groupvelocity within the tetrahedron (irreducible wedge)
     double vgx,vgy,vgz;
//   minimum doping concentration of the fitb table
     double dopmin;
//   maximum doping concentration of the fitb table
     double dopmax;
//   doping density in the current quadrant
     double dope;
//   sum of particle densities in the current quadrant
     double dens;
//   empirical correction factor for BH impurity scattering for the current
//   doping concentration
     double frickel;
//	 doping dependent empirical correction factor for BH impurity scattering
     double fitb[NFITB+1];////////////be careful !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//   energy band dependent before scattering value of charge
     double dtotp[MNB];

//   flag for generation of secondary particles
//   bool bd.seciifl;

//   flag for tunneling of electrons in the oxide (only negative x-direction)
//   bool bd.tunxfl;

//   empirical correction factor for BH impurity scattering
//   bool bd.frickfl;

//   flag for particle transition from silicon into oxide
//   bool bd.injoxfl;

//   self scattering event
     bool selfscfl;
//   flag for Brooks-Herring impurity scattering
//   bool bd.bhfl;
//   flag for use of the inverse microscopic relaxation time instead 
//   of scattering rate
//   bool bd.bhmrtfl;///taf.f scatt.f

//     number of bd.tet changes
       int ntet[MNB];

//     number of BZ changes
       int nbz[MNB];

//     number of quadrant boundary hit events
       int nquad[MNB];

//     number of phonon real and ficious scat. events
       int ntotp[MNB];

//     number of phonon real scattering events
       int  nreap[MNB];

//     number of phonon self scattering events
       int  nslfp[MNB];

//     number of total BH scattering events
       int ntotbh[MNB];

//     number of real BH scattering events
       int  nreabh[MNB];

//     number of selfscattering BH scattering events
       int  nslfbh[MNB];

//     number of II real scattering events
       int nreaii[MNB];

//     number of different phonon scattering events (electrons)
       int nsctype[MSCPRE*NBE][NBE];

//     number of different phonon scattering events (holes)
       int nsctyph[MSCPRH*NBH][NBH];

//     number of different phonon scattering events (oxide electrons)
       int nsctypoe[MSCPROE*NBOE][NBOE];

//     TOTAL Number of OXide hits above barrier energy
       int noxtotal[NPARTYP];

//     Number of transfers between silicon and oxide
       int noxreal[NPARTYP];

//     Number of reflected particles
       int noxrefl[NPARTYP];

//     charge tranfered between silicon and oxide
       double doxreal[NPARTYP];
//     accumulative tunnel charge in the oxide
       double pctun;

	void EENER(Band &bd);
	void TETTIME(double &tf,Band &bd);
	void QUADTIME(double &tf,DevSimulator &dev);
	void INWEDGE(double &xout,double &yout,double &zout,
		double xin,double yin,double zin,int matsym[6][48],int isym);
	
	void GETSYM(Band &bd);
	void GETGVE(Band &bd);
	void TETHIT(bool tunposfl,DevSimulator &dev,Band &bd);
	void GETSYMBZ(Band &bd);
	void GETSTAT(DevSimulator &dev,Band &bd);
	void QUADHIT(bool &catchfl,DevSimulator &dev,Band &bd);
	void CATCHBC(DevSimulator &dev);
	void CATCHGATEBC(DevSimulator &dev);
	void INJECT(int imot,bool &hitfl,DevSimulator &dev,Band &bd);
	void GETSYMINJECT(Band &bd);
	void GETCURR(DevSimulator &dev);
	void REFLECTBC(Band &bd);
	void PASSBC(DevSimulator &dev);
	void PERIODBC(DevSimulator &dev);
	void OETUNNEL(int ctypeold,int isymold,int itetold,DevSimulator &dev,Band &bd);
	void FORCEINJECT(bool &transfl,DevSimulator &dev,Band &bd);
	void DIFFUSEBC(Band &bd);
	void GENERATIONBC(int &icont,DevSimulator &dev);
	void GETINTER(DevSimulator &dev);
	void SELFF(double&xfs,double&yfs,DevSimulator &dev);
	void GETDD(DevSimulator &dev,Band &bd);
	void GETFIELD(double xfs,double yfs,DevSimulator &dev);
	void GETFRICKEL(Band &bd);
	void TALIND(int ityp);
	void TALINDL(int ityp,Band &bd);
	inline double EOLINT(double qq);
	void FINALSTATE(Band &bd);
	void FINALK(Band &bd);

	//  calculate the total surface scattering rate
	void CALSCATTSURFACE(DevSimulator &dev,Band &bd,int itype);
	void ESSCRT(DevSimulator &dev,Band &bd);
	void HSSCRT(DevSimulator &dev,Band &bd);
	void OESSCRT(DevSimulator &dev,Band &bd);


	double EBHMAXSCRT(DevSimulator &dev,Band &bd);
	double EBHSCRT(double eel,DevSimulator &dev,Band &bd);
	void EBHSCTR(double gamimpmax,DevSimulator &dev,Band &bd);
	void EBHSCAT(DevSimulator &dev,Band &bd);
	void EFPSCAT(int iscat,Band &bd);
	void EPSCAT(int iscat,Band &bd);
	void ESCATII(int iscat,DevSimulator &dev,Band &bd);
	void ESCTR(DevSimulator &dev,Band &bd);
	void FPFSAB(Band &bd);
	void FPFSEM(Band &bd);
	void HPSCAT(int iscat,Band &bd);
	void HSCATII(int iscat,DevSimulator &dev,Band &bd);
	void HSCTR(DevSimulator &dev,Band &bd);
	void OEPSCAT(int iscat,Band &bd);
	void OESCTR(Band &bd);
	void STRATTAU(double &tau,double &dtau,double eed,
		             int iinital,DevSimulator &dev,Band &bd);
	  
	
public:

	void FROM110TO100(double &xout,double &yout,double &zout,
				  double xin,double yin,double zin); 
	void FROM100TO110(double &xout,double &yout,double &zout,
				  double xin,double yin,double zin);
	void GENBTBTDENS(DevSimulator &dev,Band &bd);
	double FEKLOEM(double eet,Band &bd);
	double FEKLOAB(double eet,Band &bd);
	double FEKLAAB(double eet,Band &bd);
	double FEKLAEM(double eet,Band &bd);
	double FEKTOAB(double eet,Band &bd);
	double FEKTOEM(double eet,Band &bd);
	void GSEST(int rdwr,int ipar,bool checkfl,DevSimulator &dev,Band &bd);
	void singleele(int ipar,bool &pcatch,bool averfl,DevSimulator &dev,Band &bd);
	void EINIT(void);
	friend class DevSimulator;
};

void Partical::FROM110TO100(double &xout,double &yout,double &zout,
						double xin,double yin,double zin)
{
	double localx,localy,localz;
	localx=xin;
	localy=yin;
	localz=zin;
	xout=localx*SQRT2/2.0-localy*SQRT2/2.0;
	yout=localx*SQRT2/2.0+localy*SQRT2/2.0;
	zout=localz;
	return;
}

void Partical::FROM100TO110(double &xout,double &yout,double &zout,
						double xin,double yin,double zin)
{
	double localx,localy,localz;
	localx=xin;
	localy=yin;
	localz=zin;
	xout=localx*SQRT2/2.0+localy*SQRT2/2.0;
	yout=localy*SQRT2/2.0-localx*SQRT2/2.0;
	zout=localz;
	return;
}

void Partical::GSEST(int rdwr,int ipar,bool checkfl,DevSimulator &dev,Band &bd)
	{
//  purpose  : get or save particle state
//  ----------
//  parameter:
//  ----------
//    rdwr   : 1 to read , 2 to store
//    ipar   : particle number

//____local variables
      bool flag;
	  double eeold;

//____get state of current particle
      if(rdwr==1)
	  {
		 ctype=dev.ifield[IVPP][ipar];
         itet=dev.ifield[0][ipar];
         isym=dev.ifield[1][ipar];
         ijqp=dev.ifield[2][ipar];
		 xr=dev.dfield[0][ipar];
         yr=dev.dfield[1][ipar];
         xk=dev.dfield[2][ipar];
         yk=dev.dfield[3][ipar];
         zk=dev.dfield[4][ipar];
         pc=dev.dfield[5][ipar];
         dts=dev.dfield[6][ipar];
         ee=dev.dfield[7][ipar];

//_______get band index
         iband=bd.ibt[itet];

//_______determine particle type
         ipt=bd.partyp[iband];

//_______calculate two dimensional quadrant index
         jqp=ijqp/dev.ngpx;
         iqp=ijqp-dev.ngpx*jqp;


//_______get group velocity and energy at origin of cube
         GETGVE(bd);
//____save state of current particle
      }
	  else
	  {
		 dev.ifield[IVPP][ipar]=ctype;
         dev.ifield[0][ipar]=itet;
         dev.ifield[1][ipar]=isym;
         dev.ifield[2][ipar]=ijqp;

         dev.dfield[0][ipar]=xr;
         dev.dfield[1][ipar]=yr;
         dev.dfield[2][ipar]=xk;
         dev.dfield[3][ipar]=yk;
         dev.dfield[4][ipar]=zk;
         dev.dfield[5][ipar]=pc;
         dev.dfield[6][ipar]=dts;
         dev.dfield[7][ipar]=ee;
	  }

//____end of GSEST
      return;
	}


void Partical::EENER(Band &bd)
{
//      purpose:	update electron energy
	    double dx,dy,dz;
	
		
//_____transform k-vector into irreducible wedge
 		INWEDGE(dx,dy,dz,xk,yk,zk,bd.matsym,isym);///////memeber of bandry!!!!!!!!!!
//_____momentum relative to origin of tetrahedron
		dx=dx-bd.xkk[bd.tet[0][itet]];///////memeber of bandry maybe bd.tet[1]?!!!!!!!!!!
		dy=dy-bd.ykk[bd.tet[0][itet]];///////memeber of bandry!!!!!!!!!!
		dz=dz-bd.zkk[bd.tet[0][itet]];///////memeber of bandry!!!!!!!!!!
//_____new energy
		ee=ebzp+vgx*dx+vgy*dy+vgz*dz;
		return;
}

void Partical::QUADTIME(double &tf,DevSimulator &dev)
{

//     purpose: calculate time when the particle reaches the boarder
//     -------- of quadrant

//____local variables
      int idirx,idiry;
      double tfx,tfy;

//____time to reach border in x-direction
	//for 100 direction
	if(direction==100)
	{
      if(xv>0)
	  {
		  idirx=LOW;
		  tfx=(dev.gridx[iqp+1]-xr)/xv;
	  }
      else if(xv<0)
	  {
		  idirx=UP;
          tfx=-(xr-dev.gridx[iqp])/xv;
	  }
	  else
	  {
		  idirx=UP;
          tfx=1e99;
      }

//____time to reach border in y-direction
      if(yv>0)
	  {
		  idiry=RIGHT;
		  tfy=(dev.gridy[jqp+1]-yr)/yv;
	  }
      else if(yv<0)
	  {
		  idiry=LEFT;
          tfy=-(yr-dev.gridy[jqp])/yv;
	  }
	  else
	  {
		  idiry=RIGHT;
          tfy=1e99;
      }
//____chose the smaller time
	  if(tfx>=tfy)
	  {
		  idir=idiry;
		  tf=tfy;
	  }
	  else
	  {
		  idir=idirx;
		  tf=tfx;
	  }
	}

	//for 110 direction
	else if(direction==110)
	{
      if(oldxv>0)
	  {
		  idirx=LOW;
		  tfx=(dev.gridx[iqp+1]-xr)/oldxv;
	  }
      else if(oldxv<0)
	  {
		  idirx=UP;
          tfx=-(xr-dev.gridx[iqp])/oldxv;
	  }
	  else
	  {
		  idirx=UP;
          tfx=1e99;
      }

//____time to reach border in y-direction
      if(oldyv>0)
	  {
		  idiry=RIGHT;
		  tfy=(dev.gridy[jqp+1]-yr)/oldyv;
	  }
      else if(oldyv<0)
	  {
		  idiry=LEFT;
          tfy=-(yr-dev.gridy[jqp])/oldyv;
	  }
	  else
	  {
		  idiry=RIGHT;
          tfy=1e99;
      }
//____chose the smaller time
	  if(tfx>=tfy)
	  {
		  idir=idiry;
		  tf=tfy;
	  }
	  else
	  {
		  idir=idirx;
		  tf=tfx;
	  }
	}
	
	if(tf<=0)
	{
		cout<<"why  QUADTIME<0 !"<<endl;
	}
	    
//____end of QUADTIME
	return;
}

void Partical::TETTIME(double &tf,Band &bd)
{
//  purpose:calculate time when the particle reaches the boarder of tetrahedron   
//____local variables
     double  xfl,yfl,zfl,xkl,ykl,zkl;
     double tfp1,tfp2,tfp3,tfp4,nf1,nf2,nf3,nf4;
     int ibase;


// ..transform k-vector into irreduzible wedge
     INWEDGE(xkl,ykl,zkl,xk,yk,zk,bd.matsym,isym);

// ..transform electric field into irreduzible wedge
     INWEDGE(xfl,yfl,zfl,xf,yf,zf,bd.matsym,isym);

//____cal. time when the surface of the tetrahedron is reached
	ibase=itet;
	nf1=bd.datantlin[0][0][ibase]*xfl+bd.datantlin[1][0][ibase]*yfl
		+bd.datantlin[2][0][ibase]*zfl;
    nf2=bd.datantlin[0][1][ibase]*xfl+bd.datantlin[1][1][ibase]*yfl
		+bd.datantlin[2][1][ibase]*zfl;
    nf3=bd.datantlin[0][2][ibase]*xfl+bd.datantlin[1][2][ibase]*yfl
		+bd.datantlin[2][2][ibase]*zfl;
    nf4=bd.datantlin[0][3][ibase]*xfl+bd.datantlin[1][3][ibase]*yfl
      +bd.datantlin[2][3][ibase]*zfl;
    if(nf1<0)
		tfp1=(bd.datantlin[3][0][ibase]-(bd.datantlin[0][0][ibase]*xkl
			+bd.datantlin[1][0][ibase]*ykl+bd.datantlin[2][0][ibase]*zkl))/nf1;
    else
        tfp1=1e99;
 

    if(nf2<0)
        tfp2=(bd.datantlin[3][1][ibase]-(bd.datantlin[0][1][ibase]*xkl
            +bd.datantlin[1][1][ibase]*ykl+bd.datantlin[2][1][ibase]*zkl))/nf2;
    else
		tfp2=1e99;
    if(nf3<0)
        tfp3=(bd.datantlin[3][2][ibase]-(bd.datantlin[0][2][ibase]*xkl
            +bd.datantlin[1][2][ibase]*ykl+bd.datantlin[2][2][ibase]*zkl))/nf3;
    else
        tfp3=1e99;
    if(nf4<0)
		tfp4=(bd.datantlin[3][3][ibase]-(bd.datantlin[0][3][ibase]*xkl
            +bd.datantlin[1][3][ibase]*ykl+bd.datantlin[2][3][ibase]*zkl))/nf4;
    else
		tfp4=1e99;
    tf=tfp1;
    itetdir=0;//1;
	if(tfp2<tf)
	{
		itetdir=1;//2;
        tf=tfp2;
	}
	if(tfp3<tf)
	{
		itetdir=2;//3;
        tf=tfp3;
	}
	if(tfp4<tf)
	{
		itetdir=3;//4;
        tf=tfp4;
	}
	if(tf<=0)
	{
		cout<<"why  TETIME<0 !"<<endl;
	}

//____end of TETTIME
    return;
}

void Partical::INWEDGE(double &xout,double &yout,
		double &zout,double xin,double yin,double zin,int matsym[6][48],int isym)
{
// purpose map k-vector into the irreducible wedge
		//local variables
		double vsort[3];////vsort(3)
		int ibase;

		ibase=isym;
		vsort[matsym[0][ibase]-1]=xin*matsym[3][ibase];
		vsort[matsym[1][ibase]-1]=yin*matsym[4][ibase];
		vsort[matsym[2][ibase]-1]=zin*matsym[5][ibase];

		xout=vsort[0];
		yout=vsort[1];
		zout=vsort[2];

//_____End of INWEDGE
		return;
}


void Partical::GETSYM(Band &bd)
{
//  purpose:   update symmetry operation and bd.tet index
//____local variables
      int i,j,pos[4],poshilf[3],it,itl,ie,ibase;
		  //pos(0:3) so be careful it needn't chage index
      int isymold,icx,icy,icz,ic;
      double ksort[4],d1,d2,d3,d4;
		//ksort(0:3) so be careful it needn't chage index
      double xkl,ykl,zkl;
      bool foundfl;
   
      isymold=isym;
	  while((fabs(xk)+fabs(yk)+fabs(zk))>=1.50*a0pi)
	  {
		  xk*=0.999;
		  yk*=0.999;
		  zk*=0.999;
	  }

//____transform k-vector into first quadrant (xk>0, yk>0, zk>0)
      ksort[1]=fabs(xk);
      ksort[2]=fabs(yk);
      ksort[3]=fabs(zk);

//____position variable
      pos[1]=1;
      pos[2]=2;
      pos[3]=3;

//____Sort koordinates in wedge
      for(i=1;i<=2;i++)
		  for(j=2;j<=3;j++)
		  {
			if(ksort[i]<ksort[j])
			{
              ksort[0]=ksort[i];
              ksort[i]=ksort[j];
              ksort[j]=ksort[0];

              pos[0]=pos[i];
              pos[i]=pos[j];
              pos[j]=pos[0];
			}
		  };

//____determine symmetrie transformation
      poshilf[pos[1]-1]=1-1;///indmax[3][3][3][3][3][3]
      poshilf[pos[2]-1]=2-1;
      poshilf[pos[3]-1]=3-1;

	  isym=bd.indmat[sign(xk)+1][sign(yk)+1][sign(zk)+1][poshilf[0]][poshilf[1]][poshilf[2]];

//____find tetrahedron index
      INWEDGE(xkl,ykl,zkl,xk,yk,zk,bd.matsym,isym);
      foundfl=false;
//____check original tetrahedron
      if(iband==bd.ibt[itet])
	  {
         ibase=itet;
         d1=bd.datantlin[0][0][ibase]*xkl +bd.datantlin[1][0][ibase]*ykl 
            +bd.datantlin[2][0][ibase]*zkl -bd.datantlin[3][0][ibase];
         d2=bd.datantlin[0][1][ibase]*xkl +bd.datantlin[1][1][ibase]*ykl 
            +bd.datantlin[2][1][ibase]*zkl -bd.datantlin[3][1][ibase];
		 d3=bd.datantlin[0][2][ibase]*xkl +bd.datantlin[1][2][ibase]*ykl 
            +bd.datantlin[2][2][ibase]*zkl -bd.datantlin[3][2][ibase];
		 d4=bd.datantlin[0][3][ibase]*xkl +bd.datantlin[1][3][ibase]*ykl 
            +bd.datantlin[2][3][ibase]*zkl -bd.datantlin[3][3][ibase];
         if((d1>=-1e-10)&&(d2>=-1e-10)&&(d3>=-1e-10)&&(d4>=-1e-10))foundfl=true;
      }
//____search only for tetraheda contained in the given cube
//  search fine list
      if((!foundfl)&&(iband==bd.bandof[PELEC]))///iband==bd.bandof[PELEC]+1????
	  {
		  icx=(int)((xkl/a0pi-0.7)/0.3*MCX);
		  icy=(int)(ykl/(a0pi*0.075)*MCY);
		  icz=(int)(zkl/(a0pi*0.075)*MCZ);
	      if (((icx>=0)&&(icx<MCX))&&((icy>=0)&&(icy<MCY))
			&&((icz>=0)&&(icz<MCZ)))
		  {
			ic=icz+ MCZ*(icy+MCY*icx);//1+icz+ MCZ*(icy+MCY*icx);maybe error!!!!! 
			for(itl=bd.pclist[ic];itl<=bd.pclist[ic]+bd.nclist[ic]-1;itl++)
			{
				it=bd.clist[itl];
				ibase=it;
				d1=bd.datantlin[0][0][ibase]*xkl +bd.datantlin[1][0][ibase]*ykl 
			        +bd.datantlin[2][0][ibase]*zkl -bd.datantlin[3][0][ibase];
				d2=bd.datantlin[0][1][ibase]*xkl +bd.datantlin[1][1][ibase]*ykl 
			        +bd.datantlin[2][1][ibase]*zkl -bd.datantlin[3][1][ibase];
				d3=bd.datantlin[0][2][ibase]*xkl +bd.datantlin[1][2][ibase]*ykl 
			        +bd.datantlin[2][2][ibase]*zkl -bd.datantlin[3][2][ibase];
				d4=bd.datantlin[0][3][ibase]*xkl +bd.datantlin[1][3][ibase]*ykl 
					+bd.datantlin[2][3][ibase]*zkl -bd.datantlin[3][3][ibase];
				if((d1>=-1e-10)&&(d2>=-1e-10)&&(d3>=-1e-10)&&(d4>=-1e-10))
				{
					foundfl=true;
					itet=it;
					break;
				}
			}
		  }
	  }
//____search only for tetraheda containing the given energy ee
//  search fine list
		if(!foundfl)
      //ie=INT( ee * bd.dlists ) + 1
		{
			ie=(int)(ee*bd.dlists);/////???maybe (int)(ee*dalists)+1
			if(ie<0)ie=0;
			if(ie<MWLES)
			for(itl=bd.ptlists[ie][iband];
				itl<=bd.ptlists[ie][iband]+bd.ntlists[ie][iband]-1;itl++)
			{
				it=bd.tlists[itl];
				if(iband==bd.ibt[it])
				{
					ibase=it;
					d1=bd.datantlin[0][0][ibase]*xkl +bd.datantlin[1][0][ibase]*ykl 
				        +bd.datantlin[2][0][ibase]*zkl -bd.datantlin[3][0][ibase];
					d2=bd.datantlin[0][1][ibase]*xkl +bd.datantlin[1][1][ibase]*ykl 
				        +bd.datantlin[2][1][ibase]*zkl -bd.datantlin[3][1][ibase];
					d3=bd.datantlin[0][2][ibase]*xkl +bd.datantlin[1][2][ibase]*ykl 
				        +bd.datantlin[2][2][ibase]*zkl -bd.datantlin[3][2][ibase];
					d4=bd.datantlin[0][3][ibase]*xkl +bd.datantlin[1][3][ibase]*ykl 
						+bd.datantlin[2][3][ibase]*zkl -bd.datantlin[3][3][ibase];
					if((d1>=-1e-10)&&(d2>=-1e-10)&&(d3>=-1e-10)&&(d4>=-1e-10))
					{
						foundfl=true;
						itet=it;
						break;
					}
				}
			}
		}
		
//____search only for tetraheda containing the given energy ee
   		if(!foundfl)
		{
			ie=(int)(ee*bd.dlist);/////???maybe (int)(ee*dalists)+1
			if(ie<0)ie=0;
			if(ie>=MWLE)ie=MWLE-1;/////maybe ie=MWLE;??????
			for(itl=bd.ptlist[ie][iband];
				itl<=bd.ptlist[ie][iband]+bd.ntlist[ie][iband]-1;itl++)
			{
	  			it=bd.tlist[itl];
				if(iband==bd.ibt[it])
				{
					ibase=it;
					d1=bd.datantlin[0][0][ibase]*xkl +bd.datantlin[1][0][ibase]*ykl 
				        +bd.datantlin[2][0][ibase]*zkl -bd.datantlin[3][0][ibase];
					d2=bd.datantlin[0][1][ibase]*xkl +bd.datantlin[1][1][ibase]*ykl 
				        +bd.datantlin[2][1][ibase]*zkl -bd.datantlin[3][1][ibase];
					d3=bd.datantlin[0][2][ibase]*xkl +bd.datantlin[1][2][ibase]*ykl 
				        +bd.datantlin[2][2][ibase]*zkl -bd.datantlin[3][2][ibase];
					d4=bd.datantlin[0][3][ibase]*xkl +bd.datantlin[1][3][ibase]*ykl 
						+bd.datantlin[2][3][ibase]*zkl -bd.datantlin[3][3][ibase];
					if((d1>=-1e-10)&&(d2>=-1e-10)&&(d3>=-1e-10)&&(d4>=-1e-10))
					{
						foundfl=true;
						itet=it;
						break;
					}
				}
			}
		}  
//____search only for tetraheda containing the given energy ee + dee
		ie=ie+1;
		if((!foundfl)&&(ie<MWLE))
			for(itl=bd.ptlist[ie][iband];
				itl<=bd.ptlist[ie][iband]+bd.ntlist[ie][iband]-1;itl++)
			{
	  			it=bd.tlist[itl];
				if(iband==bd.ibt[it])
				{
					ibase=it;
					d1=bd.datantlin[0][0][ibase]*xkl +bd.datantlin[1][0][ibase]*ykl 
				        +bd.datantlin[2][0][ibase]*zkl -bd.datantlin[3][0][ibase];
					d2=bd.datantlin[0][1][ibase]*xkl +bd.datantlin[1][1][ibase]*ykl 
				        +bd.datantlin[2][1][ibase]*zkl -bd.datantlin[3][1][ibase];
					d3=bd.datantlin[0][2][ibase]*xkl +bd.datantlin[1][2][ibase]*ykl 
				        +bd.datantlin[2][2][ibase]*zkl -bd.datantlin[3][2][ibase];
					d4=bd.datantlin[0][3][ibase]*xkl +bd.datantlin[1][3][ibase]*ykl 
						+bd.datantlin[2][3][ibase]*zkl -bd.datantlin[3][3][ibase];
					if((d1>=-1e-10)&&(d2>=-1e-10)&&(d3>=-1e-10)&&(d4>=-1e-10))
					{
						foundfl=true;
						itet=it;
						break;
					}
				}
			}

//____search only for tetraheda containing the given energy ee - dee
        ie=ie-2;
		if((!foundfl)&&(ie>=0))
			for(itl=bd.ptlist[ie][iband];
				itl<=bd.ptlist[ie][iband]+bd.ntlist[ie][iband]-1;itl++)
			{
	  			it=bd.tlist[itl];
				if(iband==bd.ibt[it])
				{
					ibase=it;
					d1=bd.datantlin[0][0][ibase]*xkl +bd.datantlin[1][0][ibase]*ykl 
				        +bd.datantlin[2][0][ibase]*zkl -bd.datantlin[3][0][ibase];
					d2=bd.datantlin[0][1][ibase]*xkl +bd.datantlin[1][1][ibase]*ykl 
				        +bd.datantlin[2][1][ibase]*zkl -bd.datantlin[3][1][ibase];
					d3=bd.datantlin[0][2][ibase]*xkl +bd.datantlin[1][2][ibase]*ykl 
				        +bd.datantlin[2][2][ibase]*zkl -bd.datantlin[3][2][ibase];
					d4=bd.datantlin[0][3][ibase]*xkl +bd.datantlin[1][3][ibase]*ykl 
						+bd.datantlin[2][3][ibase]*zkl -bd.datantlin[3][3][ibase];
					if((d1>=-1e-10)&&(d2>=-1e-10)&&(d3>=-1e-10)&&(d4>=-1e-10))
					{
						foundfl=true;
						itet=it;
						break;
					}
				}
			}

//____search all tetraheda, if not found above
		if(!foundfl)
		for(it=0;it<bd.nt;it++)
		{
			if(iband==bd.ibt[it])
			{
				ibase=it;
				d1=bd.datantlin[0][0][ibase]*xkl +bd.datantlin[1][0][ibase]*ykl 
				    +bd.datantlin[2][0][ibase]*zkl -bd.datantlin[3][0][ibase];
				d2=bd.datantlin[0][1][ibase]*xkl +bd.datantlin[1][1][ibase]*ykl 
				    +bd.datantlin[2][1][ibase]*zkl -bd.datantlin[3][1][ibase];
				d3=bd.datantlin[0][2][ibase]*xkl +bd.datantlin[1][2][ibase]*ykl 
				    +bd.datantlin[2][2][ibase]*zkl -bd.datantlin[3][2][ibase];
				d4=bd.datantlin[0][3][ibase]*xkl +bd.datantlin[1][3][ibase]*ykl 
					+bd.datantlin[2][3][ibase]*zkl -bd.datantlin[3][3][ibase];
				if((d1>=-1e-10)&&(d2>=-1e-10)&&(d3>=-1e-10)&&(d4>=-1e-10))
				{
					foundfl=true;
					itet=it;
					break;
				}
			}
		}

  	  if(!foundfl)
	  {
		  cout<<"GETSYM: no bd.tet. found"<<endl;
		  cout<<"xk="<<xk<<endl;
		  cout<<"xk="<<xk<<endl;
		  cout<<"yk="<<yk<<endl;
		  cout<<"zk="<<zk<<endl;
		  cout<<"xk-yk+zk="<<xk-yk+zk<<endl;
		  cout<<"1.5a0pi"<<1.50*a0pi<<endl;
		  cout<<"isym"<<isym<<endl;
		  cout<<"isymold"<<isymold<<endl;
		  ofstream ftp;
		  ftp.open("output/error.txt");
		  ftp<<"GETSYM: no bd.tet. found"<<endl;
		  ftp<<"xk="<<xk<<endl;
		  ftp<<"xk="<<xk<<endl;
		  ftp<<"yk="<<yk<<endl;
		  ftp<<"zk="<<zk<<endl;
		  ftp<<"xk-yk+zk="<<xk-yk+zk<<endl;
		  ftp<<"1.5a0pi"<<1.50*a0pi<<endl;
		  ftp<<"isym"<<isym<<endl;
		  ftp<<"isymold"<<isymold<<endl;
		  exit(0);
	  }

//____End of GETSYM
    return;
}

void Partical::GETGVE(Band &bd)
{
//  purpose:  get values of group velocity and energy for actual bd.tet.
//____local variables
      int ibase;
//____get energy band parameters for actual bd.tet. (within irreduzible wedge)
      ibase=4*itet;
      vgx =bd.vt[ibase];
      vgy =bd.vt[ibase+1];
      vgz =bd.vt[ibase+2];
      ebzp=bd.vt[ibase+3];

//____new electron velocity
      bd.OUTWEDGE(xv,yv,zv,vgx,vgy,vgz,bd.matsym,isym);

	  if(direction==110)	FROM100TO110(oldxv,oldyv,oldzv,xv,yv,zv);
//____End of GETGVE
      return;
}

void Partical::TETHIT(bool tunposfl,DevSimulator &dev,Band &bd)
{
//  : Perform change of tetrahedron

//____local variables
		int isymold,itetold,ctypeold;

		ctypeold=ctype;
		isymold=isym;
		itetold=itet;
		ntet[iband]=ntet[iband]+1;

//____check for tunneling
		if(bd.tunxfl&&tunposfl&&(ipt==POXEL)&&(fabs(xk)<1e-15*a0pi)&&(xf>0)&&(xv<=0))
			OETUNNEL(ctypeold,isymold,itetold,dev,bd);
      
		if(bd.nlt[itetdir][itet]>=0)///be careful nlt!!!!!!!!!!!!!!!!
// .. only change of tetrahedron
			itet=bd.nlt[itetdir][itet];
		//else if((bd.nlt[itetdir][itet]==-4)||(bd.nlt[itetdir][itet]==-5))
		else if((bd.nlt[itetdir][itet]==-5)||(bd.nlt[itetdir][itet]==-6))
		{
//       .. change of tetrahedron, wedge and BZ
            int ii=-bd.nlt[itetdir][itet]-2;
			/////ii=-bd.nlt[itetdir][itet]-1; *nlt<0 need special care
			nbz[iband]=nbz[iband]+1;
            //xk=xk+xbz[-bd.nlt[itetdir][itet]-1][isym];
            //yk=yk+ybz[-bd.nlt[itetdir][itet]-1][isym];
            //zk=zk+zbz[-bd.nlt[itetdir][itet]-1][isym];
            xk=xk+bd.xbz[ii][isym];
            yk=yk+bd.ybz[ii][isym];
            zk=zk+bd.zbz[ii][isym];
			isym=bd.newisym[ii][isym];
            GETSYMBZ(bd);
            GETGVE(bd);
//       ..linear interpolation destroys symmetry and energy must be 
//         calculated anew
            EENER(bd);
		}
        else
//       .. change of tetrahedron and wedge
			//isym=bd.newisym[-bd.nlt[itetdir][itet]-1][isym]; nlt be careful!!!!!!!!!
            isym=bd.newisym[-bd.nlt[itetdir][itet]-2][isym];
         
      
//____update particle state
		GETGVE(bd);

//____End of TETHIT
		return;
	}
	void Partical::GETSYMBZ(Band &bd)
	{
//     purpose:   find bd.tet index after BZ scattering
//_____local variables
		int it,itl,ie,ibase;
		double d1,d2,d3,d4;
		double xkl,ykl,zkl;
		bool foundfl;
   
//_____transform into the irreduzibile wedge
		INWEDGE(xkl,ykl,zkl,xk,yk,zk,bd.matsym,isym);
		foundfl=false;
//_____check original tetrahedron
		if (iband ==bd.ibt[itet]) 
		{  
			ibase=itet;
			d1=bd.datantlin[0][0][ibase]*xkl +bd.datantlin[1][0][ibase]*ykl 
			    +bd.datantlin[2][0][ibase]*zkl -bd.datantlin[3][0][ibase];
			d2=bd.datantlin[0][1][ibase]*xkl +bd.datantlin[1][1][ibase]*ykl 
			    +bd.datantlin[2][1][ibase]*zkl -bd.datantlin[3][1][ibase];
			d3=bd.datantlin[0][2][ibase]*xkl +bd.datantlin[1][2][ibase]*ykl 
			    +bd.datantlin[2][2][ibase]*zkl -bd.datantlin[3][2][ibase];
			d4=bd.datantlin[0][3][ibase]*xkl +bd.datantlin[1][3][ibase]*ykl 
				+bd.datantlin[2][3][ibase]*zkl -bd.datantlin[3][3][ibase];
			if((d1>=-1e-10)&&(d2>=-1e-10)&&(d3>=-1e-10)&&(d4>=-1e-10))	foundfl=true;
		}
//_____first: search only for tetraheda containing the given energy ee and
//     which are part of the surface of the BZ
//_____energy index of tetraheda list
      if(!foundfl)
		{
			ie=(int)(ee*bd.dlist);/////???maybe (int)(ee*dalists)+1
			if(ie<0)ie=0;
			if(ie>=MWLE)ie=MWLE-1;/////maybe ie=MWLE;??????
			for(itl=bd.ptlist[ie][iband];
				itl<=bd.ptlist[ie][iband]+bd.ntlist[ie][iband]-1;itl++)
			{
	  			it=bd.tlist[itl];
				if(iband==bd.ibt[it])
				{
					ibase=it;
					d1=bd.datantlin[0][0][ibase]*xkl +bd.datantlin[1][0][ibase]*ykl 
				        +bd.datantlin[2][0][ibase]*zkl -bd.datantlin[3][0][ibase];
					d2=bd.datantlin[0][1][ibase]*xkl +bd.datantlin[1][1][ibase]*ykl 
				        +bd.datantlin[2][1][ibase]*zkl -bd.datantlin[3][1][ibase];
					d3=bd.datantlin[0][2][ibase]*xkl +bd.datantlin[1][2][ibase]*ykl 
				        +bd.datantlin[2][2][ibase]*zkl -bd.datantlin[3][2][ibase];
					d4=bd.datantlin[0][3][ibase]*xkl +bd.datantlin[1][3][ibase]*ykl 
						+bd.datantlin[2][3][ibase]*zkl -bd.datantlin[3][3][ibase];
					if((d1>=-1e-10)&&(d2>=-1e-10)&&(d3>=-1e-10)&&(d4>=-1e-10))
					{
						foundfl=true;
						itet=it;
						break;
					}
				}
			}
		}

//_____this should not happen, but who knows...
		if(!foundfl)GETSYM(bd);
		
//_____End of GETSYMBZ
      return;
}

void Partical::QUADHIT(bool &catchfl,DevSimulator &dev,Band &bd)
{
//____Particle has arrived at a boundary. Perform necessary action
//____local variables
		bool transfl;
		int imot,icont;
		double dd,ddmax;

		nquad[iband]=nquad[iband]+1;

//____get rules of motion for boundary
		imot=dev.motrules[idir][ijqp];
//      pass border
		if(imot==PASS)//1
		{
			GETCURR(dev);
			PASSBC(dev);
		}
// ..reflecte at border
	    else if(imot==REFLECT)//2
			REFLECTBC(bd);
// ..reflecte or scatter at border or inject into oxide
	    else if(imot==SCATTOX)//3
		{
			GETINTER(dev);
	        transfl=false;
		    if(bd.injoxfl)
			{
//    ..if possible, try transistion
				if(ipt==POXEL)
				{
					INJECT(imot,transfl,dev,bd);
//       ..force transition for particles with energy lower than 1eV
					if ((ee<1.0/eV0)&&!transfl)
					{
						noxtotal[ipt]=noxtotal[ipt]-1;
						FORCEINJECT(transfl,dev,bd);
					}
				}
				else if (ipt==PELEC)
				{
					if (ee>1.0/eV0)  INJECT(imot,transfl,dev,bd);
				}
			}

			if(sssfl)
			{
				REFLECTBC(bd);
			}
			else
			{
				if(!transfl)
				{
					double tempran=RANDNR();
					if(tempran>bd.difpr[ipt])
//					..reflect
					REFLECTBC(bd);
					else
//				  ..diffusive back scattering
					{
						dev.oldsscatter++;
						DIFFUSEBC(bd);
					}
				}
			}
			GETGVE(bd);
		 }

// ..periodic boundary condition in quadrant
		else if(imot==PERIOD)//4
			PERIODBC(dev);
// ..periodic boundary condition in quadrant plus generation of a new particle
		else if(imot==GENERATE)//5
		{
			GETCURR(dev);
	        GENERATIONBC(icont,dev);
		    PERIODBC(dev);
		}
// ..periodic boundary condition in quadrant plus generation of a new particle
		else if(imot==GENREF)//8
		{
			GETCURR(dev);
			GENERATIONBC(icont,dev);
	        REFLECTBC(bd);
		}

// ..catch particle at contact
		else if(imot==CATCH)//6
		{
			GETCURR(dev);
	        catchfl=true;
		    CATCHBC(dev);
		}

// ..catch particle at a gate contact
		else if(imot==CATCHGATE)//7
		{
			GETCURR(dev);
			catchfl=true;
			CATCHGATEBC(dev);
		}
		else
		{
			cout<<"error QUADHIT: No such rule of motion"<<endl;
			exit(0);
		}

//____End of QUADHIT
      return;
}

void Partical::CATCHBC(DevSimulator &dev)
{
//  catch boundary condition
//____local variables
      int icont;

      if(idir==UP)icont=dev.gridcont[ijqp-1];////////????or [ijqp]?????????
      else if(idir==RIGHT)icont=dev.gridcont[ijqp+2*dev.ngpx];
      else if(idir==LOW)icont=dev.gridcont[ijqp+2];
      else if(idir==LEFT)icont=dev.gridcont[ijqp-dev.ngpx];
	  int carriertype;
	  carriertype=ctype;
	  dev.ncatch[icont]=dev.ncatch[icont]+1;
	  dev.dcatch[icont]=dev.dcatch[icont]+pc;

	  dev.dcatchh[carriertype][icont]=dev.dcatchh[carriertype][icont]+pc;

//____End of CATCHBC
      return;
}

void Partical::CATCHGATEBC(DevSimulator &dev)
{
//  catchgate boundary condition

//____local variables
      int icont;

	  if(idir==UP)icont=dev.gridcont[ijqp];////////????or [ijqp]?????????
      else if(idir==RIGHT)icont=dev.gridcont[ijqp+dev.ngpx];
      else if(idir==LOW)icont=dev.gridcont[ijqp+1];
      else if(idir==LEFT)icont=dev.gridcont[ijqp];

	  int carriertype;
	  carriertype=ctype;
	  dev.ncatch[icont]=dev.ncatch[icont]+1;
	  dev.dcatch[icont]=dev.dcatch[icont]+pc;

	  dev.dcatchh[carriertype][icont]=dev.dcatchh[carriertype][icont]+pc;

//____End of CATCHGATEBC
      return;
	}
void Partical::INJECT(int imot,bool &hitfl,DevSimulator &dev,Band &bd)
{

//____local variables
      bool finishfl;
      int ibandold,itetold,isymold,iptold,ib;
      double eenew,xkold,ykold,zkold;
      double xfold,yfold,zfold,tettf,tf,ts;
      double tt,vf,tfmax;

      noxtotal[ipt]=noxtotal[ipt]+1;

//____save inital particle state, electric field
      ibandold=iband;
      itetold=itet;
      isymold=isym;
      iptold =ipt;
      xkold=xk;
      ykold=yk;
      zkold=zk;
      xfold=xf;
      yfold=yf;
      zfold=zf;

//____new particle energy
      eenew=ee;

//____new particle type
      if((imot==SCATTOX)&&(dev.mat[ijqp]==SILICON))
	  {
         if(ipt==PELEC)ipt=POXEL;
	  }
      else if((imot== SCATTOX)&&(dev.mat[ijqp]==OXIDE))
	  {
         if(ipt==POXEL)ipt=PELEC;
	  }
      
//____search in direction perpindicular to the interface
      if((idir==UP)||(idir==LOW))
	  {
         xf=1;
         yf=0;
         zf=0;
	  }
      else if((idir==RIGHT)||(idir==LEFT))
	  {
		 xf=0;
         yf=1;
         zf=0;
	  }
//____maximum flight time through BZ (electric field=1 in MC units)
      tfmax=2*a0pi;

//____traveltime of particle in k-space (a complete pass of BZ takes tfmax)
      hitfl=false;
	  
	  for(iband=bd.bandof[ipt];iband<bd.bandof[ipt]+bd.nband[ipt];iband++)/////////be careful!!!!
	  {
		 ts=0;
         finishfl=false;
         xk=xkold;
         yk=ykold;
         zk=zkold;
         if ((idir==UP)||(idir==LOW))xk=0;
         else if((idir==RIGHT)||(idir==LEFT))yk=0;
         GETSYMINJECT(bd);

		 while(!finishfl)
		 {
            GETGVE(bd);
            EENER(bd);
            TETTIME(tettf,bd);
            tf=tettf;
            vf=xv*xf+yv*yf+zv*zf;
            if (!(vf==0))
               tt=(eenew-ee)/vf;
            else
               tt=tfmax;
            if(tt<0) tt=tfmax;
            if((tt<tf)&&!finishfl) 
			{
               ib=iband;
               xk=xk+xf*tt;
               yk=yk+yf*tt;
               zk=zk+zf*tt;
               hitfl=true;
			   break;
       		}
            if (!finishfl)
			{
			   xk=xk+xf*tf;
               yk=yk+yf*tf;
               zk=zk+zf*tf;
               ts=ts+tf;
               TETHIT(false,dev,bd);
               if(ts>=tfmax)finishfl=true;
          	}
         }
      }

//____transfer between materials possible (new state found)
	  if(hitfl)
	  {
         noxreal[iptold]=noxreal[iptold]+1;
         doxreal[iptold]=doxreal[iptold]+pc;
//    ..restore old field
         xf=xfold;
         yf=yfold;
         zf=zfold;
         iband=ib;
//    ..get new energy
         GETGVE(bd);
         EENER(bd);
         if(fabs(ee-eenew)>1e-10) 
		 {
			 cout<<"error INJECT: New state has wrong energy";exit(0);
		 }
//    ..change quadrant
         if(idir==UP)
		 {
            xr=dev.gridx[iqp];
            iqp=iqp-1;
            ijqp=ijqp-1;
            if(xv>0.0) 
			{
               xk=-xk;     
               isym=bd.indmat[-bd.matsym[3][isym]+1][bd.matsym[4][isym]+1]
				   [bd.matsym[5][isym]+1][bd.matsym[0][isym]-1]
				   [bd.matsym[1][isym]-1][bd.matsym[2][isym]-1];
            }
		 }
         else if(idir==RIGHT) 
		 {
			yr=dev.gridy[jqp+1];
            jqp=jqp+1;
            ijqp=ijqp+dev.ngpx;
            if(yv<0.0) 
			{
               yk=-yk;
			   /*
               isym=bd.indmat( bd.matsym(4,isym),-bd.matsym(5,isym), 
     >                        bd.matsym(6,isym),      
     >                        bd.matsym(1,isym), bd.matsym(2,isym), 
     >                        bd.matsym(3,isym))
            ENDIF
			*/
			   isym=bd.indmat[bd.matsym[3][isym]+1][-bd.matsym[4][isym]+1]
				   [bd.matsym[5][isym]+1][bd.matsym[0][isym]-1]
				   [bd.matsym[1][isym]-1][bd.matsym[2][isym]-1];
			}
		 }
         else if(idir==LOW)
		 {
            xr=dev.gridx[iqp+1];
            iqp=iqp+1;
            ijqp=ijqp+1;
            if(xv<0.0)
			{
               xk=-xk;
			   /*
               isym=bd.indmat(-bd.matsym(4,isym), bd.matsym(5,isym), 
     >                        bd.matsym(6,isym),      
     >                        bd.matsym(1,isym), bd.matsym(2,isym), 
     >                        bd.matsym(3,isym))
            ENDIF
			*/
			   isym=bd.indmat[-bd.matsym[3][isym]+1][bd.matsym[4][isym]+1]
				   [bd.matsym[5][isym]+1][bd.matsym[0][isym]-1]
				   [bd.matsym[1][isym]-1][bd.matsym[2][isym]-1];
			}
		 }
         else if(idir==LEFT) 
		 {
            yr=dev.gridy[jqp];
            jqp=jqp-1;
            ijqp=ijqp-dev.ngpx;
            if(yv>0.0)
			{
               yk=-yk;
			   /*
               isym=bd.indmat( bd.matsym(4,isym),-bd.matsym(5,isym), 
     >                        bd.matsym(6,isym),      
     >                        bd.matsym(1,isym), bd.matsym(2,isym), 
     >                        bd.matsym(3,isym))
            */
			   isym=bd.indmat[bd.matsym[3][isym]+1][-bd.matsym[4][isym]+1]
				   [bd.matsym[5][isym]+1][bd.matsym[0][isym]-1]
				   [bd.matsym[1][isym]-1][bd.matsym[2][isym]-1];
			}
		 }
         if(!(dev.mat[ijqp]==bd.typemat[ipt]))
		 {
			 cout<<"error INJECT: Particle not in ";exit(0);
         }
	  }
//____no new state found, restore old values
      else
	  {
         noxrefl[iptold]=noxrefl[iptold]+1;
         iband=ibandold;
         itet=itetold;
         isym=isymold;
         ipt=iptold;
         xk=xkold;
         yk=ykold;
         zk=zkold;
         xf=xfold;
         yf=yfold;
         zf=zfold;
         GETGVE(bd);
         EENER(bd);
       }

//____adjust particle number for particle types
      dev.npar[iptold]=dev.npar[iptold]-1;
      dev.npar[ipt]=dev.npar[ipt]+1;

//____End of INJECT
	  return;
}
	
void Partical::GETSYMINJECT(Band &bd)
{
//  purpose:   update symmetry operation and bd.tet index
//  --------   after injection into silicon or oxide

//____local variables
      int i,j,pos[4],poshilf[3],it,itl,ibase;
      double  ksort[4],d1,d2,d3,d4;
      double  xkl,ykl,zkl;
      bool foundfl;
   
//____transform k-vector into first quadrant (xk>0, yk>0, zk>0)
      ksort[1]=fabs(xk);
      ksort[2]=fabs(yk);
      ksort[3]=fabs(zk);

//____position variable
      pos[1]=1;
      pos[2]=2;
      pos[3]=3;

//____Sort koordinates in wedge
      for(i=1;i<=2;i++)
	  {
		  for(j=2;j<=3;j++)
		  {
			if(ksort[i]<ksort[j])
			{
              ksort[0]=ksort[i];
              ksort[i]=ksort[j];
              ksort[j]=ksort[0];

              pos[0]=pos[i];
              pos[i]=pos[j];
              pos[j]=pos[0];
			}
		  }
	  }

//____determine symmetrie transformation
      poshilf[pos[1]-1]=1-1;///indmax[3][3][3][3][3][3]
      poshilf[pos[2]-1]=2-1;
      poshilf[pos[3]-1]=3-1;
      /*
	  isym=bd.indmat( NINT(SIGN(1.d0,xk)), NINT(SIGN(1.d0,yk)),
     >               NINT(SIGN(1.d0,zk)), 
     >               poshilf((1)), poshilf((2)), poshilf((3)) )//??????????
		*/
	  isym=bd.indmat[sign(xk)+1][sign(yk)+1][sign(zk)+1][poshilf[0]][poshilf[1]][poshilf[2]];

//____find tetrahedron index
      INWEDGE(xkl,ykl,zkl,xk,yk,zk,bd.matsym,isym);
      foundfl=false;
	  
//____check original tetrahedron
      if(iband==bd.ibt[itet])
	  {
         ibase=itet;
         d1=bd.datantlin[0][0][ibase]*xkl +bd.datantlin[1][0][ibase]*ykl 
            +bd.datantlin[2][0][ibase]*zkl -bd.datantlin[3][0][ibase];
         d2=bd.datantlin[0][1][ibase]*xkl +bd.datantlin[1][1][ibase]*ykl 
            +bd.datantlin[2][1][ibase]*zkl -bd.datantlin[3][1][ibase];
		 d3=bd.datantlin[0][2][ibase]*xkl +bd.datantlin[1][2][ibase]*ykl 
            +bd.datantlin[2][2][ibase]*zkl -bd.datantlin[3][2][ibase];
		 d4=bd.datantlin[0][3][ibase]*xkl +bd.datantlin[1][3][ibase]*ykl 
            +bd.datantlin[2][3][ibase]*zkl -bd.datantlin[3][3][ibase];
         if((d1>=-1e-10)&&(d2>=-1e-10)&&(d3>=-1e-10)&&(d4>=-1e-10))	foundfl=true;
       
      }

//____search kz=0 plane
      if(!foundfl)
      for(itl=bd.pzlist[iband];itl<=bd.pzlist[iband]+bd.nzlist[iband]-1;itl++)
	  {  
		 it=bd.zlist[itl];
         ibase=it;

		 d1=bd.datantlin[0][0][ibase]*xkl +bd.datantlin[1][0][ibase]*ykl 
            +bd.datantlin[2][0][ibase]*zkl -bd.datantlin[3][0][ibase];
         d2=bd.datantlin[0][1][ibase]*xkl +bd.datantlin[1][1][ibase]*ykl 
            +bd.datantlin[2][1][ibase]*zkl -bd.datantlin[3][1][ibase];
		 d3=bd.datantlin[0][2][ibase]*xkl +bd.datantlin[1][2][ibase]*ykl 
            +bd.datantlin[2][2][ibase]*zkl -bd.datantlin[3][2][ibase];
		 d4=bd.datantlin[0][3][ibase]*xkl +bd.datantlin[1][3][ibase]*ykl 
            +bd.datantlin[2][3][ibase]*zkl -bd.datantlin[3][3][ibase];
         if((d1>=-1e-10)&&(d2>=-1e-10)&&(d3>=-1e-10)&&(d4>=-1e-10))
		 {
			 if(foundfl) 
			 {

			 }
			 foundfl=true;
			 itet=it;
         }
      }
      
//____this should not happen, but who knows
      if(!foundfl)GETSYM(bd);

//____End of GETSYMINJECT
      return;
}

void Partical::GETINTER(DevSimulator &dev)
{      
//  Purpose: get Si/SiO2 interface statistics
//  -------  

//____local variables
      int ier;
      
//    energy index of energy distribution function
//    ier=(int)(ee/dev.edfemax[ipt]*NEREDF)+1
//    IF (ier .GT. NEREDF) ier=NEREDF  
	  ier=(int)(ee/dev.edfemax[ipt]*NEREDF);//////////////////be careful!!!!!!!!!
      if(ier>NEREDF-1) ier=NEREDF-1;      
      if((idir==RIGHT)||(idir==LEFT))
	  {
		  if((xr-dev.gridx[iqp])>(dev.gridx[iqp+1] - xr))
			  dev.bedfx[ier][iqp+1][ipt]=dev.bedfx[ier][iqp+1][ipt]+fabs(pc);
          else
			  dev.bedfx[ier][iqp][ipt]=dev.bedfx[ier][iqp][ipt]+fabs(pc);
      }
      if((idir==UP)||(idir==LOW))
	  {
		  if((yr-dev.gridy[jqp])>(dev.gridy[jqp+1]-yr))
		  	  dev.bedfy[ier][jqp+1][ipt]=dev.bedfy[ier][jqp+1][ipt]+fabs(pc);
          else
			  dev.bedfy[ier][jqp][ipt]=dev.bedfy[ier][jqp][ipt]+fabs(pc);
       }

//____end of GETINTER
      return;
}

void Partical::GETCURR(DevSimulator &dev)
{
//    Purpose: get current statistics

      if(idir==RIGHT)
	  {
         if((xr-dev.gridx[iqp])>(dev.gridx[iqp+1]-xr))
            dev.cury[ijqp+dev.ngpx+1][ipt]+=pc;
         else
            dev.cury[ijqp+dev.ngpx][ipt]+=pc;
	  }
      else if(idir==LEFT)
	  {
		 if((xr-dev.gridx[iqp])>(dev.gridx[iqp+1] - xr))
            dev.cury[ijqp+1][ipt]-=pc;
         else
            dev.cury[ijqp][ipt]-=pc;
      }
      else if(idir==LOW)
	  {
         if((yr-dev.gridy[jqp])>(dev.gridy[jqp+1]-yr))
            dev.curx[ijqp+dev.ngpx+1][ipt]+=pc;
         else
            dev.curx[ijqp+1][ipt]+=pc;
      }
      else if(idir==UP)
	  {
         if((yr-dev.gridy[jqp])>(dev.gridy[jqp+1]-yr))
            dev.curx[ijqp+dev.ngpx][ipt]-=pc;
         else
            dev.curx[ijqp][ipt]-=pc;
      }
      
//____end of GETCURR
      return;
}

void Partical::OETUNNEL(int ctypeold,int isymold,int itetold,DevSimulator &dev,Band &bd)
{
//____tunneling in oxide (only in negative x-direction)

//____local variables
      int ij,i,j,iq,nq,idl[MNGPX+MNGPY];
      double potw[MNGPX+MNGPY+1],dl[MNGPX+MNGPY];

//    potw[MNGPX+MNGPY] is a temp cell use in c ; in fortan it is potw[0];
	  
      double prob,xrp,bef,dup,dlow,potwl,potwu;

//____tunneling only in negative x-direction
      ij=ijqp;
      INDTD(ij,i,j);
// .. potential at the particle location
      if(potgalfl)
	  {
		potw[MNGPX+MNGPY]=(dev.galpot[ij]+dev.imagepot[ij]+dev.galpot[ij+dev.ngpx]+dev.imagepot[ij+dev.ngpx])*0.50
               * (dev.gridx[i+1]-xr)/(dev.gridx[i+1]-dev.gridx[i])
               +(dev.galpot[ij+1]+dev.imagepot[ij+1]+dev.galpot[ij+dev.ngpx+1]+dev.imagepot[ij+dev.ngpx+1])*0.50
               *(xr-dev.gridx[i])/(dev.gridx[i+1]-dev.gridx[i]);
	  } 
	  else
	  {
		potw[MNGPX+MNGPY]=(dev.pot[ij]+dev.imagepot[ij]+dev.pot[ij+dev.ngpx]+dev.imagepot[ij+dev.ngpx])*0.50
               * (dev.gridx[i+1]-xr)/(dev.gridx[i+1]-dev.gridx[i])
               +(dev.pot[ij+1]+dev.imagepot[ij+1]+dev.pot[ij+dev.ngpx+1]+dev.imagepot[ij+dev.ngpx+1])*0.50
               *(xr-dev.gridx[i])/(dev.gridx[i+1]-dev.gridx[i]);
	  }

      iq=0;
	  for(i=0;i<=iqp-1;i++)
	  {
         ij=INDS(i,j);
         if(potgalfl) 
		 {
            potwu=(dev.galpot[ij]+dev.imagepot[ij])*0.50 
                    + (dev.galpot[ij+dev.ngpx]+dev.imagepot[ij+dev.ngpx])*0.50;
            potwl=(dev.galpot[ij+1]+dev.imagepot[ij+1])*0.50 
                    + (dev.galpot[ij+dev.ngpx+1]+dev.imagepot[ij+dev.ngpx+1])*0.50;
         }
		 else
		 {
            potwu=(dev.pot[ij]+dev.imagepot[ij])*0.50 
                 + (dev.pot[ij+dev.ngpx]+dev.imagepot[ij+dev.ngpx])*0.50;
            potwl=(dev.pot[ij+1]+dev.imagepot[ij+1])*0.50 
                 + (dev.pot[ij+dev.ngpx+1]+dev.imagepot[ij+dev.ngpx+1])*0.5;
         }
         if((dev.mat[ij]==OXIDE)&&((potwu<=potw[MNGPX+MNGPY])||(potwl<= potw[MNGPX+MNGPY])))
		 {
            iq=iq+1;
            idl[iq]=i;
            if(potwu>potw[MNGPX+MNGPY])
			{
               xrp=(potw[MNGPX+MNGPY]-potwl)/(potwl-potwu)*(dev.gridx[i+1]-dev.gridx[i]); 
               dl[iq]=-xrp;
               potw[iq]=potw[MNGPX+MNGPY];
               potw[iq+1]=potwl;
			}
			else if(potwl>potw[MNGPX+MNGPY])
			{
               xrp=(potw[MNGPX+MNGPY]-potwu)/(potwl-potwu)*(dev.gridx[i+1]-dev.gridx[i]);
               dl[iq]=xrp;
               potw[iq]=potwu;
               potw[iq+1]=potw[MNGPX+MNGPY];
			}
            else
			{
               potw[iq]=potwu;
               potw[iq+1]=potwl;
               dl[iq]=dev.gridx[i+1]-dev.gridx[i];
            }
         }
      }

// ..cal. tunneling probability with WKB method for a linearly interpolated 
//   potential
      if(iq>1)
	  {
         nq=iq;
//    ..cal. position, where the energy is again positive
         //xrp=dev.gridx(idl(2)) - dl(1)??????????????????????????
	     xrp=dev.gridx[idl[1]]-dl[0];///////////???
//    ..tunneling prob.
         prob=0;
         //DO iq=1, nq
		 for(iq=0;iq<=nq-1;iq++)
		 {
            bef=(potw[iq+1]-potw[iq])/dl[iq];
            //dup= SQRT(-(potw(iq)-potw[0]))**3
            dup=pow((-(potw[iq]-potw[MNGPX+MNGPY])),3.0/2.0);
			//dlow=SQRT(-(potw(iq+1)-potw(0)))**3
			dlow=pow((-(potw[iq+1]-potw[MNGPX+MNGPY])),3.0/2.0);
            prob=prob+(dup-dlow)/bef;
         }
         prob=exp(-4.0/3.0*sqrt(2.0*bd.dmox)*prob);
//    ..split the primary particle in one which tunnels and one wich does nit
//      tunnel and assign appropriate particle charges
         if(prob>1e-30)
		 {
			 pctun=pctun+pc*prob;
             if((dev.npar0+1)>=MNPAR)
			 {

			 }
			 else
			 {
				dev.npar0=dev.npar0+1;
                dev.npar[ipt]=dev.npar[ipt]+1;

				dev.ifield[IVPP][dev.npar0]=ctypeold;
                dev.ifield[0][dev.npar0]=itetold;
                dev.ifield[1][dev.npar0]=isymold;
                dev.ifield[2][dev.npar0]=INDS(idl[0],j);/////////////maybe idl[1]????
   
				dev.dfield[0][dev.npar0]=xrp;
				dev.dfield[1][dev.npar0]=yr;
				dev.dfield[2][dev.npar0]=xk;
				dev.dfield[3][dev.npar0]=yk;
				dev.dfield[4][dev.npar0]=zk;
				dev.dfield[5][dev.npar0]=pc*prob;
				dev.dfield[6][dev.npar0]=dts;
				dev.dfield[7][dev.npar0]=ee;
				pc=pc*(1.0-prob);
			}
		 }
	  }

//____End of OETUNNEL
     return;
}
 
void Partical::FORCEINJECT(bool &transfl,DevSimulator &dev,Band &bd)
{
//-----: Force injection of an oxide electron into silicon
//____local variables
      int ibandold,isymold,itetold,iptold;
      double xkold,ykold,zkold;

      noxtotal[ipt]=noxtotal[ipt]+1;

//____search finalstate on equi energy surface within the conduction bands
// ..save initial particle state
      ibandold=iband;
      iptold=ipt;
      itetold=itet;
      isymold=isym;
      xkold=xk;
      ykold=yk;
      zkold=zk;
// ..chose random direction on equi energy surface
      transfl=true;
      ipt=PELEC;
      iband=bd.bandof[ipt];

	  while((!transfl)&&(iband<(bd.nband[ipt]+bd.bandof[ipt])))
	  {
      if(material==0)	TALIND(1);
	  else if(material==1)	TALINDL(1,bd);
      FINALSTATE(bd);
// ..if FINALSTATE fails restore inital particle state
      if(selfscfl)
	  {
         iband=ibandold;
         itet=itetold;
         isym=isymold;
         xk=xkold;
         yk=ykold;
         zk=zkold;
         transfl=false;
      }
      GETGVE(bd);
// ..chose particle velocity to point into actual device region
      if(idir==UP)
	  {
            if(xv==0.0)
			{
				iband=ibandold;
				itet=itetold;
	            isym=isymold;
		        xk=xkold;
			    yk=ykold;
				zk=zkold;
	            transfl=false;
		        GETGVE(bd);
			}
			if(xv<0.0) 
			{
               xk=-xk;     
               isym=bd.indmat[-bd.matsym[3][isym]+1][bd.matsym[4][isym]+1]
				   [bd.matsym[5][isym]+1][bd.matsym[0][isym]-1]
				   [bd.matsym[1][isym]-1][bd.matsym[2][isym]-1];
            }
		 }
         else if(idir==RIGHT) 
		 {
			if(yv==0.0)
			{
				iband=ibandold;
				itet=itetold;
	            isym=isymold;
		        xk=xkold;
			    yk=ykold;
				zk=zkold;
	            transfl=false;
		        GETGVE(bd);
			}
			if(yv>0.0) 
			{
               yk=-yk;
			   isym=bd.indmat[bd.matsym[3][isym]+1][-bd.matsym[4][isym]+1]
				   [bd.matsym[5][isym]+1][bd.matsym[0][isym]-1]
				   [bd.matsym[1][isym]-1][bd.matsym[2][isym]-1];
			}
		 }
         else if(idir==LOW)
		 {
            if(xv==0.0)
			{
				iband=ibandold;
				itet=itetold;
	            isym=isymold;
		        xk=xkold;
			    yk=ykold;
				zk=zkold;
	            transfl=false;
		        GETGVE(bd);
			}
			if(xv>0.0)
			{
               xk=-xk;
			   isym=bd.indmat[-bd.matsym[3][isym]+1][bd.matsym[4][isym]+1]
				   [bd.matsym[5][isym]+1][bd.matsym[0][isym]-1]
				   [bd.matsym[1][isym]-1][bd.matsym[2][isym]-1];
			}
		 }
         else if(idir==LEFT) 
		 {
            if(yv==0.0)
			{
				iband=ibandold;
				itet=itetold;
	            isym=isymold;
		        xk=xkold;
			    yk=ykold;
				zk=zkold;
	            transfl=false;
		        GETGVE(bd);
			}
			if(yv<0.0)
			{
               yk=-yk;
			   isym=bd.indmat[bd.matsym[3][isym]+1][-bd.matsym[4][isym]+1]
				   [bd.matsym[5][isym]+1][bd.matsym[0][isym]-1]
				   [bd.matsym[1][isym]-1][bd.matsym[2][isym]-1];
			}
		 }
      GETGVE(bd);
	  iband+=1;
	  }///end while

//____pass particle into silicon
      if(transfl)
	  {
         noxreal[iptold]=noxreal[iptold]+1;
         doxreal[iptold]=doxreal[iptold]+pc;
         if(idir==UP)
		  {
		     xr=dev.gridx[iqp];
			 iqp=iqp-1;
	         ijqp=ijqp-1;
		  }
	      else if(idir==RIGHT)
		  {
	         yr=dev.gridy[jqp+1];
		     jqp=jqp+1;
			 ijqp=ijqp+dev.ngpx;
		  }
		  else if(idir==LOW)
		  {
		     xr=dev.gridx[iqp+1];
			 iqp=iqp+1;
	         ijqp=ijqp+1;
		  }
	      else if(idir==LEFT)
		  {
			 yr=dev.gridy[jqp];
			jqp=jqp-1;
	         ijqp=ijqp-dev.ngpx;
		  }
	  }
      else
	  {
         ipt=iptold;
	  }
      
//____adjust particle number for particle types
      dev.npar[iptold]-=1;
      dev.npar[ipt]+=1;

//____End of FORCEINJECT
      return;
}

void Partical::REFLECTBC(Band &bd)
{
//  Refelct the particle at the boundary
	  double ooxk,ooyk,oozk;
	  if(direction==110)
	  {
		  ooxk=xk;
		  ooyk=yk;
		  oozk=zk;
	  }	
	  if(idir==UP)
	  {
		    if(direction==100)
			{
				xk=-xk;
				isym=bd.indmat[-bd.matsym[3][isym]+1][bd.matsym[4][isym]+1]
					[bd.matsym[5][isym]+1][bd.matsym[0][isym]-1]
					[bd.matsym[1][isym]-1][bd.matsym[2][isym]-1];
			}
			else if(direction==110)
			{
				xk=-ooyk;
				yk=-ooxk;
				isym=bd.indmat[bd.matsym[4][isym]+1][bd.matsym[3][isym]+1]
					[bd.matsym[5][isym]+1][bd.matsym[1][isym]-1]
					[bd.matsym[0][isym]-1][bd.matsym[2][isym]-1];
				isym=bd.indmat[-bd.matsym[3][isym]+1][bd.matsym[4][isym]+1]
					[bd.matsym[5][isym]+1][bd.matsym[0][isym]-1]
					[bd.matsym[1][isym]-1][bd.matsym[2][isym]-1];
				isym=bd.indmat[bd.matsym[3][isym]+1][-bd.matsym[4][isym]+1]
					[bd.matsym[5][isym]+1][bd.matsym[0][isym]-1]
					[bd.matsym[1][isym]-1][bd.matsym[2][isym]-1];
			}
        
	  }
      else if(idir==RIGHT) 
	  {
		  if(direction==100)
		  {
			  yk=-yk;
			  isym=bd.indmat[bd.matsym[3][isym]+1][-bd.matsym[4][isym]+1]
				  [bd.matsym[5][isym]+1][bd.matsym[0][isym]-1]
				  [bd.matsym[1][isym]-1][bd.matsym[2][isym]-1];
		  }
		  else if(direction==110)
		  {
			  xk=ooyk;
			  yk=ooxk;
			  isym=bd.indmat[bd.matsym[4][isym]+1][bd.matsym[3][isym]+1]
				  [bd.matsym[5][isym]+1][bd.matsym[1][isym]-1]
				  [bd.matsym[0][isym]-1][bd.matsym[2][isym]-1];
		  }
	  }
      else if(idir==LOW)
	  {
		    if(direction==100)
			{
				xk=-xk;
				isym=bd.indmat[-bd.matsym[3][isym]+1][bd.matsym[4][isym]+1]
					[bd.matsym[5][isym]+1][bd.matsym[0][isym]-1]
					[bd.matsym[1][isym]-1][bd.matsym[2][isym]-1];
			}
			else if(direction==110)
			{
				xk=-ooyk;
				yk=-ooxk;
				isym=bd.indmat[bd.matsym[4][isym]+1][bd.matsym[3][isym]+1]
					[bd.matsym[5][isym]+1][bd.matsym[1][isym]-1]
					[bd.matsym[0][isym]-1][bd.matsym[2][isym]-1];
				isym=bd.indmat[-bd.matsym[3][isym]+1][bd.matsym[4][isym]+1]
					[bd.matsym[5][isym]+1][bd.matsym[0][isym]-1]
					[bd.matsym[1][isym]-1][bd.matsym[2][isym]-1];
				isym=bd.indmat[bd.matsym[3][isym]+1][-bd.matsym[4][isym]+1]
					[bd.matsym[5][isym]+1][bd.matsym[0][isym]-1]
					[bd.matsym[1][isym]-1][bd.matsym[2][isym]-1];
			}
	  }
      else if(idir==LEFT) 
	  {
		  if(direction==100)
		  {
			  yk=-yk;
			  isym=bd.indmat[bd.matsym[3][isym]+1][-bd.matsym[4][isym]+1]
				  [bd.matsym[5][isym]+1][bd.matsym[0][isym]-1]
				  [bd.matsym[1][isym]-1][bd.matsym[2][isym]-1];
		  }
		  else if(direction==110)
		  {
			  xk=ooyk;
			  yk=ooxk;
			  isym=bd.indmat[bd.matsym[4][isym]+1][bd.matsym[3][isym]+1]
				  [bd.matsym[5][isym]+1][bd.matsym[1][isym]-1]
				  [bd.matsym[0][isym]-1][bd.matsym[2][isym]-1];
		  }
	  }
      GETGVE(bd);

//____End of REFLECTBC
      return;
}
    
void Partical::PERIODBC(DevSimulator &dev)
{
//  periodic boundary condition
      if(idir==UP)
         xr=dev.gridx[iqp+1];
      else if(idir==RIGHT)
         yr=dev.gridy[jqp];
      else if(idir==LOW)
         xr=dev.gridx[iqp];
      else if(idir==LEFT)
         yr=dev.gridy[jqp+1];
//____End of PERIODBC
	  return;
}
//====

void Partical::PASSBC(DevSimulator &dev)
{
//  pass boundary condition
      if(idir==UP)
	  {
         xr=dev.gridx[iqp];
         iqp=iqp-1;
         ijqp=ijqp-1;
	  }
      else if(idir==RIGHT)
	  {
         yr=dev.gridy[jqp+1];
         jqp=jqp+1;
         ijqp=ijqp+dev.ngpx;
	  }
      else if(idir==LOW)
	  {
         xr=dev.gridx[iqp+1];
         iqp=iqp+1;
         ijqp=ijqp+1;
	  }
      else if(idir==LEFT)
	  {
         yr=dev.gridy[jqp];
         jqp=jqp-1;
         ijqp=ijqp-dev.ngpx;
	  }
//____End of PASSBC
      return;
}
//====

void Partical::DIFFUSEBC(Band &bd)
{
//  diffusive boundary condition
//____local variables
      int ibandold,itetold,isymold;
      double xkold,ykold,zkold;

	  double ooxk,ooyk,oozk;

//    save initial particle state
      ibandold=iband;
      itetold=itet;
      isymold=isym;
      xkold=xk;
      ykold=yk;
      zkold=zk;
//    chose random direction on equi energy surface
	  if(material==0)	TALIND(1);
	  else if(material==1)	TALINDL(1,bd);
      FINALSTATE(bd);
	  if(direction==100)
	  {

	  }
	  else if(direction==110)
	  {
		  ooxk=xk;
		  ooyk=yk;
		  oozk=zk;
	  }
//    if FINALSTATE fails restore inital particle state
      if(selfscfl)
	  {
         iband=ibandold;
         itet=itetold;
         isym=isymold;
         xk=xkold;
         yk=ykold;
         zk=zkold;
         //transfl=false;
      }
      GETGVE(bd);
// ..chose particle velocity to point into actual device region
	if(direction==100)
	{
      if(idir==UP)
	  {
            if(xv==0.0)
			{
				iband=ibandold;
				itet=itetold;
	            isym=isymold;
		        xk=xkold;
			    yk=ykold;
				zk=zkold;
	       //     transfl=false;
		        GETGVE(bd);
			}
			if(xv<0.0) 
			{
               xk=-xk;     
               isym=bd.indmat[-bd.matsym[3][isym]+1][bd.matsym[4][isym]+1]
				   [bd.matsym[5][isym]+1][bd.matsym[0][isym]-1]
				   [bd.matsym[1][isym]-1][bd.matsym[2][isym]-1];
            }
		 }
         else if(idir==RIGHT) 
		 {
			if(yv==0.0)
			{
				iband=ibandold;
				itet=itetold;
	            isym=isymold;
		        xk=xkold;
			    yk=ykold;
				zk=zkold;
	         //   transfl=false;
		        GETGVE(bd);
			}
			if(yv>0.0) 
			{
               yk=-yk;
			   isym=bd.indmat[bd.matsym[3][isym]+1][-bd.matsym[4][isym]+1]
				   [bd.matsym[5][isym]+1][bd.matsym[0][isym]-1]
				   [bd.matsym[1][isym]-1][bd.matsym[2][isym]-1];
			}
		 }
         else if(idir==LOW)
		 {
            if(xv==0.0)
			{
				iband=ibandold;
				itet=itetold;
	            isym=isymold;
		        xk=xkold;
			    yk=ykold;
				zk=zkold;
//              transfl=false;
		        GETGVE(bd);
			}
			if(xv>0.0)
			{
               xk=-xk;
			   isym=bd.indmat[-bd.matsym[3][isym]+1][bd.matsym[4][isym]+1]
				   [bd.matsym[5][isym]+1][bd.matsym[0][isym]-1]
				   [bd.matsym[1][isym]-1][bd.matsym[2][isym]-1];
			}
		 }
         else if(idir==LEFT) 
		 {
            if(yv==0.0)
			{
				iband=ibandold;
				itet=itetold;
	            isym=isymold;
		        xk=xkold;
			    yk=ykold;
				zk=zkold;
//              transfl=false;
		        GETGVE(bd);
			}
			if(yv<0.0)
			{
               yk=-yk;
			   isym=bd.indmat[bd.matsym[3][isym]+1][-bd.matsym[4][isym]+1]
				   [bd.matsym[5][isym]+1][bd.matsym[0][isym]-1]
				   [bd.matsym[1][isym]-1][bd.matsym[2][isym]-1];
			}
		 }
	}
	else if(direction==110)
	{
      if(idir==UP)
	  {
            if(oldxv==0.0)
			{
				iband=ibandold;
				itet=itetold;
	            isym=isymold;
		        xk=xkold;
			    yk=ykold;
				zk=zkold;
//              transfl=false;
		        GETGVE(bd);
			}
			if(oldxv<0.0) 
			{
//				xk=-xk;     
				xk=-ooyk;
				yk=-ooxk;
				isym=bd.indmat[bd.matsym[4][isym]+1][bd.matsym[3][isym]+1]
					[bd.matsym[5][isym]+1][bd.matsym[1][isym]-1]
					[bd.matsym[0][isym]-1][bd.matsym[2][isym]-1];
				isym=bd.indmat[-bd.matsym[3][isym]+1][bd.matsym[4][isym]+1]
					[bd.matsym[5][isym]+1][bd.matsym[0][isym]-1]
					[bd.matsym[1][isym]-1][bd.matsym[2][isym]-1];
				isym=bd.indmat[bd.matsym[3][isym]+1][-bd.matsym[4][isym]+1]
					[bd.matsym[5][isym]+1][bd.matsym[0][isym]-1]
					[bd.matsym[1][isym]-1][bd.matsym[2][isym]-1];
            }
		 }
         else if(idir==RIGHT) 
		 {
			if(oldyv==0.0)
			{
				iband=ibandold;
				itet=itetold;
	            isym=isymold;
		        xk=xkold;
			    yk=ykold;
				zk=zkold;
//              transfl=false;
		        GETGVE(bd);
			}
			if(oldyv>0.0) 
			{
//				yk=-yk;
				xk=ooyk;
				yk=ooxk;
				isym=bd.indmat[bd.matsym[4][isym]+1][bd.matsym[3][isym]+1]
					[bd.matsym[5][isym]+1][bd.matsym[1][isym]-1]
					[bd.matsym[0][isym]-1][bd.matsym[2][isym]-1];
			}
		 }
         else if(idir==LOW)
		 {
            if(oldxv==0.0)
			{
				iband=ibandold;
				itet=itetold;
	            isym=isymold;
		        xk=xkold;
			    yk=ykold;
				zk=zkold;
//              transfl=false;
		        GETGVE(bd);
			}
			if(oldxv>0.0)
			{
//				xk=-xk;     
				xk=-ooyk;
				yk=-ooxk;
				isym=bd.indmat[bd.matsym[4][isym]+1][bd.matsym[3][isym]+1]
					[bd.matsym[5][isym]+1][bd.matsym[1][isym]-1]
					[bd.matsym[0][isym]-1][bd.matsym[2][isym]-1];
				isym=bd.indmat[-bd.matsym[3][isym]+1][bd.matsym[4][isym]+1]
					[bd.matsym[5][isym]+1][bd.matsym[0][isym]-1]
					[bd.matsym[1][isym]-1][bd.matsym[2][isym]-1];
				isym=bd.indmat[bd.matsym[3][isym]+1][-bd.matsym[4][isym]+1]
					[bd.matsym[5][isym]+1][bd.matsym[0][isym]-1]
					[bd.matsym[1][isym]-1][bd.matsym[2][isym]-1];
			}
		 }
         else if(idir==LEFT) 
		 {
            if(oldyv==0.0)
			{
				iband=ibandold;
				itet=itetold;
	            isym=isymold;
		        xk=xkold;
			    yk=ykold;
				zk=zkold;
//              transfl=false;
		        GETGVE(bd);
			}
			if(oldyv<0.0)
			{
//				yk=-yk;
				xk=ooyk;
				yk=ooxk;
				isym=bd.indmat[bd.matsym[4][isym]+1][bd.matsym[3][isym]+1]
					[bd.matsym[5][isym]+1][bd.matsym[1][isym]-1]
					[bd.matsym[0][isym]-1][bd.matsym[2][isym]-1];
			}
		 }
	}
	GETGVE(bd);

//____End of DIFFUSEBC
      return;
}

void Partical::GENERATIONBC(int &icont,DevSimulator &dev)
{
//  generation boundary condition
      dev.npar0=dev.npar0+1;
      dev.npar[ipt]=dev.npar[ipt]+1;
	  int carriertype;
	  carriertype=ctype;
      if(dev.npar0>MNPAR)
	  {
	  }
      if(idir==UP)
	  {
         icont=dev.gridcont[ijqp+1];
         dev.ngen[icont]=dev.ngen[icont]+1;
         dev.dgen[icont]=dev.dgen[icont]+pc;
		 dev.dgenn[carriertype][icont]=dev.dgenn[carriertype][icont]+pc;
         dev.ifield[2][dev.npar0-1]=ijqp-1;
//       dev.ifield[3][dev.npar0]=ijqp;// no need
         dev.dfield[0][dev.npar0-1]=dev.gridx[iqp];
         dev.dfield[1][dev.npar0-1]=yr;
	  }
      else if(idir==RIGHT)
	  {
         icont=dev.gridcont[ijqp];
         dev.ngen[icont]=dev.ngen[icont]+1;
         dev.dgen[icont]=dev.dgen[icont]+pc;
		 dev.dgenn[carriertype][icont]=dev.dgenn[carriertype][icont]+pc;
         dev.ifield[2][dev.npar0-1]=ijqp+dev.ngpx;
         dev.dfield[0][dev.npar0-1]=xr;
         dev.dfield[1][dev.npar0-1]=dev.gridy[jqp+1];
	  }
      else if(idir==LOW)
	  {
         icont=dev.gridcont[ijqp];
         dev.ngen[icont]=dev.ngen[icont]+1;
         dev.dgen[icont]=dev.dgen[icont]+pc;
		 dev.dgenn[carriertype][icont]=dev.dgenn[carriertype][icont]+pc;
         dev.ifield[2][dev.npar0-1]=ijqp+1;
         dev.dfield[0][dev.npar0-1]=dev.gridx[iqp+1];
         dev.dfield[1][dev.npar0-1]=yr;
	  }
      else if(idir==LEFT)
	  {
         icont=dev.gridcont[ijqp+dev.ngpx];
         dev.ngen[icont]=dev.ngen[icont]+1;
         dev.dgen[icont]=dev.dgen[icont]+pc;
		 dev.dgenn[carriertype][icont]=dev.dgenn[carriertype][icont]+pc;
         dev.ifield[2][dev.npar0-1]=ijqp-dev.ngpx;
         dev.dfield[0][dev.npar0-1]=xr;
         dev.dfield[1][dev.npar0-1]=dev.gridy[jqp];
      }
	  dev.ifield[IVPP][dev.npar0-1]=ctype;
      dev.ifield[0][dev.npar0-1]=itet;
      dev.ifield[1][dev.npar0-1]=isym;
      dev.dfield[2][dev.npar0-1]=xk;
      dev.dfield[3][dev.npar0-1]=yk;
      dev.dfield[4][dev.npar0-1]=zk;
      dev.dfield[5][dev.npar0-1]=pc;
      dev.dfield[6][dev.npar0-1]=dts;
      dev.dfield[7][dev.npar0-1]=ee;

//____End of GENERATIONBC
      return;
}

void Partical::singleele(int ipar,bool &pcatch,bool averfl,DevSimulator &dev,Band &bd)
{
	bool catchfl,phfl,quadfl,tetfl,impfl,ssfl;
	int i, j, ireal, icpu, ierr, iptype, icont, jcont;
	double tettf, quadtf, phtf, imptf, sstf, tf, tempdts;
    double xfs, yfs, phrnl, imprnl,ssnl;
    double gamimpmax;
    
    GSEST(1,ipar,true,dev,bd);

//____calculate field without self force (approx.)
	if(selffl)SELFF(xfs,yfs,dev);

//____init flags
    catchfl=false;
    tetfl=true;
    quadfl=true;
    phfl=true;
    impfl=true;
	ssfl=true;

	ssregionfl=false;
	if((xr<=(xmaxss1*1e-6/spr0)&&xr>=(xminss1*1e-6/spr0)&&yr<=(ymaxss1*1e-6/spr0)&&yr>=(yminss1*1e-6/spr0))
	 ||(xr<=(xmaxss2*1e-6/spr0)&&xr>=(xminss2*1e-6/spr0)&&yr<=(ymaxss2*1e-6/spr0)&&yr>=(yminss2*1e-6/spr0)))
	{
		ssregionfl=true;
	}
	else
	{
		ssregionfl=false;
	}
	channelfl=false;
	if(ballisticfl&&yr>=dev.gridy[dev.jbcont[2]-1]&&yr<=dev.gridy[dev.jecont[2]+1])
	{
		channelfl=true;
	}
	else
	{
		channelfl=false;
	}

	GETFIELD(xfs,yfs,dev);
//  get doping density and particle density
    GETDD(dev,bd);

//__no change of quadrant

//__time until next tetrahedron change
	if(tetfl)
	{
		TETTIME(tettf,bd);
        tetfl=false;
	}

//____time until next quadrant change
    if(quadfl)
	{
		if(bulkfl)
			quadtf=dt;
        else
            QUADTIME(quadtf,dev);
        quadfl=false;
	}

//__time until next phonon scattering process (variable Gamma scheme)
    if(phfl)
	{
		phrnl=-log(RANDNR());
        phfl=false;
	}
	phtf=phrnl/bd.gamtet[itet];
	if(channelfl)
	{
		phtf=2*dt;
	}

//__time until next impurity scattering process (variable Gamma scheme)
//  impurity scattering only in the lowest conduction band
    if(bd.bhfl&&(iband==bd.bandof[PELEC]))///////////??bd.bandof[PELEC]+1?????????????
									///iband= the first band of electric
	{
		if(impfl)
		{
            imprnl=-log(RANDNR());
            impfl=false;
        }
        gamimpmax=EBHMAXSCRT(dev,bd);
        imptf=imprnl/gamimpmax;
	}
    else
	{
		imptf=2*dt;
        gamimpmax=1.0/scrt0;
    }
	if(channelfl)
	{
		imptf=2*dt;
	}

//__time until next surface scattering
	if(sssfl&&ssregionfl)//total surface scattering flag
	{
		if(ssfl)//surface scattering flag
		{
			ssnl=-log(RANDNR());
			ssfl=false;
		}
		CALSCATTSURFACE(dev,bd,ipt);
		surfacemax=ssrv;
		if(ssrv==0)	surfacemax=dt;
		sstf=ssnl/surfacemax;
	}
	else
	{
		sstf=2*dt;
	}

	if(channelfl)
	{
		sstf=2*dt;
	}

//__free flight time
    tf=MIN(tettf,quadtf,phtf,imptf,sstf,dts);

	if(tf<0)
	{
		catchfl=true;
		dts=0.0;
		pcatch=catchfl;
		cout<<"why tf<0"<<endl;
		tf=0.0;
	}

//__adjust free flight times
    tettf=tettf-tf;
	quadtf=quadtf-tf;
    phtf=phtf-tf;
    imptf=imptf-tf;
    dts=dts-tf;
	sstf=sstf-tf;

//__variable gamma schemes
    phrnl=phrnl-bd.gamtet[itet]*tf;
    imprnl=imprnl-gamimpmax*tf;
	ssnl=ssnl-surfacemax*tf;

//__update k-vektor
    xk=xk+xf*tf;
    yk=yk+yf*tf;
    zk=zk+zf*tf;

//__update energy
	if(direction==100)	ee=ee+(xv*xf+yv*yf+zv*zf)*tf;
	else if(direction==110)	ee=ee+(oldxv*oldxf+oldyv*oldyf+oldzv*oldzf)*tf;

	if(ee<0)	
		cout<<"error ee<0"<<endl;
		

//__update r-vektor
	if(direction==100)
	{
		if(!bulkfl)xr=xr+xv*tf;
		if((!odxfl)&&(!bulkfl))yr=yr+yv*tf;
	}
	else if(direction==110)
	{
		if(!bulkfl)xr=xr+oldxv*tf;
		if((!odxfl)&&(!bulkfl))yr=yr+oldyv*tf;
	}
//  big while
	while(dts>1e-10)
	{
		if((xr<=(xmaxss1*1e-6/spr0)&&xr>=(xminss1*1e-6/spr0)&&yr<=(ymaxss1*1e-6/spr0)&&yr>=(yminss1*1e-6/spr0))
		||(xr<=(xmaxss2*1e-6/spr0)&&xr>=(xminss2*1e-6/spr0)&&yr<=(ymaxss2*1e-6/spr0)&&yr>=(yminss2*1e-6/spr0)))
		{
			ssregionfl=true;
		}
		else
		{
			ssregionfl=false;
		}

		/*
		if(ballisticfl&&yr>=dev.gridx[dev.jbcont[2]]&&yr<=dev.gridx[dev.jecont[2]])
		{
			channelfl=true;
		}
		else
		{
			channelfl=false;
		}
		*/

//_______Tetrahedron change
         if(tettf<=1e-10)
		 {
            TETHIT(true,dev,bd);
            tetfl=true;
            quadfl=true;
          }
         
//_______Quadrant  change
         else if(quadtf<=1e-10)
		 {
			QUADHIT(catchfl,dev,bd);
//          particle catched by contact ?
            if(catchfl)
			{
				dts=0.0;
				pcatch=catchfl;
				break;
			}
            tetfl=true;
            quadfl=true;

            GETFIELD(xfs,yfs,dev);
//    ..get doping density and particle density
			GETDD(dev,bd);
         
		}
//_______phonon scattering?
        else if(phtf<=1e-10)
		{
			dev.phscatter++;
//       ..just before value for impact ionization current
            if (ipt==PELEC)ESCTR(dev,bd);
            if (ipt==PHOLE)HSCTR(dev,bd);
            if (ipt==POXEL)OESCTR(bd);
            if (!selfscfl)
			{
               tetfl=true;
               quadfl=true;
            }
            phfl=true;
		}

//_______impurity scattering?
		else if(imptf<=1e-10)
		{   
			dev.bhscatter++;
			selfscfl=true;
            if(ipt==PELEC)	
			{
				EBHSCTR(gamimpmax,dev,bd);
			}
            if(!selfscfl)
			{
               tetfl=true;
               quadfl=true;
            }
            impfl=true;
		}

//_______sourface scattering
		else if(sstf<=1e-10)
		{
			if(sssfl)
			{
				dev.sscatter++;
				if (ipt==PELEC)	ESSCRT(dev,bd);
				if (ipt==PHOLE)	HSSCRT(dev,bd);
//				if (ipt==POXEL) OESSCRT(dev,bd);
				if (!selfscfl)
				{
					tetfl=true;
					quadfl=true;
				}
				ssfl=true;
			}
			else
			{
			}
		}

//____time until next tetrahedron change
		if(tetfl)
		{
			TETTIME(tettf,bd);
	        tetfl=false;
		}

//____time until next quadrant change
	    if(quadfl)
		{
			if(bulkfl)
				quadtf=dt;
	        else
		        QUADTIME(quadtf,dev);
			quadfl=false;
		}
//____time until next phonon scattering process (variable Gamma scheme)
	    if(phfl)
		{
			phrnl=-log(RANDNR());
	        phfl=false;
		}
		phtf=phrnl/bd.gamtet[itet];

		if(channelfl)
		{
			phtf=2*dt;
		}

//____time until next impurity scattering process (variable Gamma scheme)
//    ..impurity scattering only in the lowest conduction band
		if(bd.bhfl&&(iband==bd.bandof[PELEC]))///////////??bd.bandof[PELEC]+1?????????????
									///iband= the first band of electric
		{
			if(impfl)
			{
				imprnl=-log(RANDNR());
	            impfl=false;
		    }
			gamimpmax=EBHMAXSCRT(dev,bd);
	        imptf=imprnl/gamimpmax;
		}
	    else
		{
			imptf=2*dt;
	        gamimpmax=1.0/scrt0;
		}

		if(channelfl)
		{
			imptf=2*dt;
		}

//____time until next surface scattering
		if(sssfl&&ssregionfl)//total surface scattering flag
		{
			if(ssfl)//surface scattering flag
			{
				ssnl=-log(RANDNR());
				ssfl=false;
			}
			CALSCATTSURFACE(dev,bd,ipt);
			surfacemax=ssrv;
			if(ssrv==0)	surfacemax=dt;
			sstf=ssnl/surfacemax;
		}
		else
		{
			sstf=2*dt;
		}

		if(channelfl)
		{
			sstf=2*dt;
		}

//____free flight time
		tf=MIN(tettf,quadtf,phtf,imptf,sstf,dts);
	
		if(tf<0)
		{
			catchfl=true;
			dts=0.0;
			pcatch=catchfl;
			cout<<"why tf<0"<<endl;
			tf=0.0;
		}	

//____adjust free flight times
	    tettf=tettf-tf;
		quadtf=quadtf-tf;
	    phtf=phtf-tf;
		imptf=imptf-tf;
		sstf=sstf-tf;
	    dts=dts-tf;

//____variable gamma schemes
		phrnl=phrnl-bd.gamtet[itet]*tf;
		imprnl=imprnl-gamimpmax*tf;
		ssnl=ssnl-surfacemax*tf;

//____update k-vektor
	    xk=xk+xf*tf;
		yk=yk+yf*tf;
	    zk=zk+zf*tf;

//____update energy
		if(direction==100)	ee=ee+(xv*xf+yv*yf+zv*zf)*tf;
		else if(direction==110)	ee=ee+(oldxv*oldxf+oldyv*oldyf+oldzv*oldzf)*tf;

		if(ee<0)	cout<<"error ee<0"<<endl;

//____update r-vektor
		if(direction==100)
		{
			if(!bulkfl)xr=xr+xv*tf;
			if((!odxfl)&&(!bulkfl))yr=yr+yv*tf;
		}
		else if(direction==110)
		{
			if(!bulkfl)xr=xr+oldxv*tf;
			if((!odxfl)&&(!bulkfl))yr=yr+oldyv*tf;
		}
	}//end big while

	if(catchfl)
	{
//     last particle?
       if(ipar==(dev.npar0-1))
	   {
           dev.npar0=dev.npar0-1;
           dev.npar[ipt]=dev.npar[ipt]-1;
	    }
        else
		{
//         get last particle and procced with it
           dev.npar[ipt]=dev.npar[ipt]-1;
           dev.npar0=dev.npar0-1;
		   GSEST(1,dev.npar0,true,dev,bd);
           GSEST(2,ipar,false,dev,bd);
		}
	}
	else
	{
//____get particle statistic
		if(averfl)
		{
			GETSTAT(dev,bd);
		}
//____save particle state
		GSEST(2,ipar,false,dev,bd);
	
		//from starttime sum all the velocity in x direction and y direction
		if(cdt>=starttime)
		{
			if(ipt==PELEC)
			{
				ebandsum[ijqp][iband]=ebandsum[ijqp][iband]+fabs(pc);
				if(direction==100)
				{
					eveloxsum[ijqp]+=xv*pc;
					eveloysum[ijqp]+=yv*pc;

					eveloxsumm[ctype][ijqp]+=xv*pc;
					eveloysumm[ctype][ijqp]+=yv*pc;
				}
				else if(direction==110)
				{
					eveloxsum[ijqp]+=oldxv*pc;
					eveloysum[ijqp]+=oldyv*pc;

					eveloxsumm[ctype][ijqp]+=oldxv*pc;
					eveloysumm[ctype][ijqp]+=oldyv*pc;
				}

				for(i=0;i<statnoiseregionnumber;i++)
				{
					if(xr>=xminnoise[i]&&xr<=xmaxnoise[i]&&yr>=yminnoise[i]&&yr<=ymaxnoise[i])
					{
						evsum[i][cdt-starttime]+=sqrt(xv*xv+yv*yv)*pc;
						evnumber[i][cdt-starttime]+=pc;
					}
				}
				eenergysum[ijqp]+=ee*pc;
				eenergysumm[ctype][ijqp]+=ee*pc;
				enumber[ijqp]+=pc;
				enumberr[ctype][ijqp]+=pc;//sum times
			}
			else if(ipt==PHOLE)
			{
				hbandsum[ijqp][iband]=hbandsum[ijqp][iband]+fabs(pc);
				if(direction==100)
				{
					hveloxsum[ijqp]+=xv*pc;
					hveloysum[ijqp]+=yv*pc;

					hveloxsumm[ctype][ijqp]+=xv*pc;
					hveloysumm[ctype][ijqp]+=yv*pc;
				}
				else if(direction==110)
				{
					hveloxsum[ijqp]+=oldxv*pc;
					hveloysum[ijqp]+=oldyv*pc;

					hveloxsumm[ctype][ijqp]+=oldxv*pc;
					hveloysumm[ctype][ijqp]+=oldyv*pc;
				}

				for(i=0;i<statnoiseregionnumber;i++)
				{
					if(xr>=xminnoise[i]&&xr<=xmaxnoise[i]&&yr>=yminnoise[i]&&yr<=ymaxnoise[i])
					{
						hvsum[i][cdt-starttime]+=sqrt(xv*xv+yv*yv)*pc;
						hvnumber[i][cdt-starttime]+=pc;
					}
				}
				henergysum[ijqp]+=ee*pc;
				henergysumm[ctype][ijqp]+=ee*pc;
				hnumber[ijqp]+=pc;
				hnumberr[ctype][ijqp]+=pc;//sum times
			}
		}
	}
	return;
}

void Partical::SELFF(double&xfs,double&yfs,DevSimulator &dev)
{
//   Purpose: cal. selfforce free electric field
	double rx,ry;

    rx=(xr-dev.gridx[iqp])/(dev.gridx[iqp+1]-dev.gridx[iqp]);
    ry=(yr-dev.gridy[jqp])/(dev.gridy[jqp+1]-dev.gridy[jqp]);
    xfs=dev.xfself[UL][ijqp]*(1.0-rx)*(1.0-ry)
        +dev.xfself[UR][ijqp]*(1.0-rx)*ry
        +dev.xfself[LR][ijqp]*rx*ry
        +dev.xfself[LL][ijqp]*rx*(1.0-ry);
    yfs=dev.yfself[UL][ijqp]*(1.0-rx)*(1.0-ry)
        +dev.yfself[UR][ijqp]*(1.0-rx)*ry
        +dev.yfself[LR][ijqp]*rx*ry
        +dev.yfself[LL][ijqp]*rx*(1.0-ry);

//_____end of SELFF
    return;
}

void Partical::GETDD(DevSimulator &dev,Band &bd)
{
//     Purpose: get doping density and particle density

    dope=0.25*(dev.donor[UL][ijqp]+dev.donor[LL][ijqp]
		+dev.donor[UR][ijqp]+dev.donor[LR][ijqp]
        +dev.accep[UL][ijqp]+dev.accep[LL][ijqp]
        +dev.accep[UR][ijqp]+dev.accep[LR][ijqp]);

    if(potgalfl)
	{
		dens=MAX(0.25*(dev.galecon[ijqp]+dev.galecon[ijqp+1]
					+dev.galecon[ijqp+dev.ngpx]+dev.galecon[ijqp+dev.ngpx+1]
					+dev.galhcon[ijqp]+dev.galhcon[ijqp+1]
			        +dev.galhcon[ijqp+dev.ngpx]+dev.galhcon[ijqp+dev.ngpx+1]),
		            0.50*dope);
	}
    else
		dens = MAX(dev.quaddens[ijqp],0.5*dope);
   GETFRICKEL(bd);

//_____end of GETDD
	return;
}

void Partical::GETFIELD(double xfs,double yfs,DevSimulator &dev)
{ 
//     Purpose: cal. electric field
//     -------  
	if(bulkfl)
	{
        xf=charsign[ipt]*xfieldbulk;
        yf=charsign[ipt]*yfieldbulk;
        zf=charsign[ipt]*zfieldbulk;
	}
	else
	{
		if(dev.selfforcefl[dev.numreg[ijqp]])
			xf=charsign[ipt]*xfs;
		else
		{
            xf=charsign[ipt]*dev.xfield[ijqp];
            //if(quantfl)xf=xf+charsign[ipt]*xqf[ijqp][ipt];
			if(quantumeffectqpfl||quantumeffectfmfl||quantumeffectbmfl||quantumeffectqnfl)xf=xf+charsign[ipt]*dev.xqf[ijqp][ipt];
		}
        if(odxfl)
            yf=charsign[ipt]*yfieldbulk;
        else
            if(dev.selfforcefl[dev.numreg[ijqp]])
               yf=charsign[ipt]*yfs;
            else
			{
               yf=charsign[ipt]*dev.yfield[ijqp];
               //if(quantfl)yf=yf+charsign[ipt]*yqf[ijqp][ipt];
			   if(quantumeffectqpfl||quantumeffectfmfl||quantumeffectbmfl||quantumeffectqnfl)yf=yf+charsign[ipt]*dev.yqf[ijqp][ipt];
            }
        zf=0.0;
	}
	// if direction 110, need transfer electric field from 110 to 100
	if(direction==110)
	{
		oldxf=xf;
		oldyf=yf;
		oldzf=zf;
		FROM110TO100(xf,yf,zf,oldxf,oldyf,oldzf);
	}
//_____end of GETFIELD
	return;
}

void Partical::GETFRICKEL(Band &bd)
{    
//    Purpose: get the frickel parameter
//____local variables
	double xx;
    int itab;

    if(bd.frickfl)
	{
		if(dope<=dopmin)
            frickel=fitb[0];
        else if(dope>=dopmax)
            frickel=fitb[NFITB];
        else
		{
            xx=log(dope/dopmin)/log(dopmax/dopmin)*NFITB;
            itab=(int)xx;
            xx=xx-itab;
            frickel=fitb[itab]*(1.0-xx)+fitb[itab+1]*xx;
		}
	}
    else
		frickel=1.0;
   
//_____end of GETFRICKEL
    return;
}

void Partical::TALIND(int ityp)
{
	  //int ityp;

	  //Purpose: Calculation of the symetrie-operation
	  //Parameter :  ityp = typ off scattering-prozess
	  //	     1 = intra-valley
	  //         2 = inter-valley g-Phonon
	  //         3 = inter-valley f-Phonon
	  //         isym = indize of symetrie-operation

	  //local variables
      int talend[5],ital;
      //double RANDNR;

	  //number of valley 
      if((fabs(xk)>=fabs(yk))&&(fabs(xk)>=fabs(zk)))
	  {	
         if(xk>=0.0)
		 {	 
			//location valley  1
            if(ityp==1)
			{
				ital=1;
			}
            if(ityp==2)
			{
				ital=2;
			}
            if(ityp==3)
			{
                talend[1]=3;
                talend[2]=4;
                talend[3]=5;
                talend[4]=6;
                ital=talend[(int)(RANDNR()*4.0+1.0)];
            }
         }
         if(xk<0.0)
		 {	
			//location valley  2
            if(ityp==1)
			{
				ital=2;
			}
            if(ityp==2)
			{
				ital=1;
			}
            if(ityp==3)
			{
                talend[1]=3;
                talend[2]=4;
                talend[3]=5;
                talend[4]=6;
                ital=talend[(int)(RANDNR()*4.0+1.0)];
			}
		 }
	  }	
      if((fabs(yk)>=fabs(zk))&&(fabs(yk)>=fabs(xk)))
      {
		  if(yk>=0.0)
		  {
			  //location valley  3
              if(ityp==1)
			  {
				  ital=3;
			  }
              if(ityp==2)
			  {
				  ital=4;
			  }
              if(ityp==3)
			  {
				  talend[1]=1;
                  talend[2]=2;
                  talend[3]=5;
                  talend[4]=6;
                  ital=talend[(int)(RANDNR()*4.0+1.0)];
			  }
		  }
         if(yk<0.0)
		 {
			 //location valley  4
            if(ityp==1)
			{
				ital=4;
			}
            if(ityp==2)
			{
				ital=3;
			}
            if(ityp==3)
			{
               talend[1]=1;
               talend[2]=2;
               talend[3]=5;
               talend[4]=6;
               ital=talend[(int)(RANDNR()*4.0+1.0)];
			}
		 }
	  }	

      if((fabs(zk)>fabs(yk))&&(fabs(zk)>fabs(xk)))
	  {	
         if(zk>=0.0)
		 {
			 //location valley  5
            if(ityp==1)
			{
				ital=5;
			}
            if(ityp==2)
			{
				ital=6;
			}
            if(ityp==3)
			{
               talend[1]=1;
               talend[2]=2;
               talend[3]=3;
               talend[4]=4;
               ital=talend[(int)(RANDNR()*4.0+1.0)];
            }
		 }	
         if(zk<0.0)
		 {	
			//location valley  6
            if(ityp==1)
			{
				ital=6;
			}
            if(ityp==2)
			{
				ital=5;
			}
            if(ityp==3)
			{
               talend[1]=1;
               talend[2]=2;
               talend[3]=3;
               talend[4]=4;
               ital=talend[(int)(RANDNR()*4.0+1.0)];
			}
		 }	
	  }	
      //isym=(int)(RANDNR()*8.0+1.0)+8*(ital-1);
	  isym=(int)(RANDNR()*8.0)+8*(ital-1);//maxim=47 minim=0;????????????????

//____end of TALIND
	  return;
}

void Partical::TALINDL(int ityp,Band &bd)
{
//      SUBROUTINE TALINDL (ityp)
//      INTEGER ityp

//     Purpose: Calculation of the symetrie-operation of L point
//     =======
//     Parameter :  ityp = typ off scattering-prozess
//                     1 = intra-valley
//                     2 = inter-valley
//                  isym = indize of symetrie-operation

//      INCLUDE "mcgii.par" 
//      INCLUDE "mcsim.par" 
//      INCLUDE "fullband.com"
//      INCLUDE "ele.com" 

//______local variables
//      INTEGER talend(4), ital,i
		int talend[5],ital,i;
//      DOUBLE PRECISION RANDNR
      
//_____number of valley ??
		if(xk>=0)
		{
			if(yk>=0)
			{
				if(zk>=0)
				{
					if(ityp==1) isym=bd.mapsym[0][(int)(RANDNR()*6e0+43e0)-1]-1;
					if(ityp==2)	isym=bd.mapsym[0][(int)(RANDNR()*42e0+1e0)-1]-1;
				}
				else
				{
					if(ityp==1) isym=bd.mapsym[1][(int)(RANDNR()*6e0+43e0)-1]-1;
					if(ityp==2)	isym=bd.mapsym[1][(int)(RANDNR()*42e0+1e0)-1]-1;

				}
			}
			else
			{
				if(zk>=0)
				{
					if(ityp==1) isym=bd.mapsym[3][(int)(RANDNR()*6e0+43e0)-1]-1;
					if(ityp==2)	isym=bd.mapsym[3][(int)(RANDNR()*42e0+1e0)-1]-1;
				}
				else
				{
					if(ityp==1) isym=bd.mapsym[2][(int)(RANDNR()*6e0+43e0)-1]-1;
					if(ityp==2)	isym=bd.mapsym[2][(int)(RANDNR()*42e0+1e0)-1]-1;
				}
			}
		}
		else
		{
			if(yk>=0)
			{
				if(zk>=0)
				{
					if(ityp==1) isym=bd.mapsym[4][(int)(RANDNR()*6e0+43e0)-1]-1;
					if(ityp==2)	isym=bd.mapsym[4][(int)(RANDNR()*42e0+1e0)-1]-1;
				}
				else
				{
					if(ityp==1) isym=bd.mapsym[5][(int)(RANDNR()*6e0+43e0)-1]-1;
					if(ityp==2)	isym=bd.mapsym[5][(int)(RANDNR()*42e0+1e0)-1]-1;
				}
			}
			else
			{
				if(zk>=0)
				{
					if(ityp==1) isym=bd.mapsym[7][(int)(RANDNR()*6e0+43e0)-1]-1;
					if(ityp==2)	isym=bd.mapsym[7][(int)(RANDNR()*42e0+1e0)-1]-1;
				}
				else
				{
					if(ityp==1) isym=bd.mapsym[6][(int)(RANDNR()*6e0+43e0)-1]-1;
					if(ityp==2)	isym=bd.mapsym[6][(int)(RANDNR()*42e0+1e0)-1]-1;
				}
			}
		}
//_____end of TALINDL
	return;
}

void Partical::FINALSTATE(Band &bd)
{
//Purpose: - calculate the final state for ee and iband on an aequienergy surface
      bool smallfl;
      int  itab,nted,count2;
      double dostet;
      double maxdostet;

//    all bands have energies larger than 0eV
//    energies smaller than 0eV are sometimes chosen due 
//    to discretization errors in the scattering rates

      if(ee<=0.0)
	  {	
//    no final state possible for this energy
         selfscfl =true;
	  }

//    get maximum dos of tetraheda for the given energy
	  else
	  {
		  maxdostet=bd.CALDOSTETMAX(ee,iband);

//        energy index of tetraheder list
//        fine list below 120meV
	      itab=(int)(ee*bd.dlists);//+1;///////minim=0;
	      if(itab<MWLES)//maxim=MWLES-1
		  {	
			 smallfl=true;
		  }
		  else
		  {
//           otherwise coarse list
		     itab=(int)(ee*bd.dlist);//+1;//minim=0;
			 if(itab>=MWLE)
			 {
				 itab=MWLE-1;//maxim=MWLE-1
			 }
	         smallfl=false;
		  }	
//        Number of tetraheder in list
	      if(smallfl)
		  {	
			 nted=bd.ntlists[itab][iband];
		  }	
		  else
		  {	
		     nted=bd.ntlist[itab][iband];
		  }	
		  count2=0;
		
		  if(nted==0) //nted could be 0;
		  {
//________no final state possible for this band and energy
			 selfscfl=true;
	      }

//________chose random tetraheder from list
		  else
		  {	
			  bool nfindfl=true;
			  while(nfindfl)
			  {
				  if(smallfl)
				  {
						itet=bd.tlists[(int)(nted*RANDNR())+bd.ptlists[itab][iband]];
				  }	
				  else
				  {
						itet=bd.tlist[(int)(nted*RANDNR())+bd.ptlist[itab][iband]];
				  }	
				  count2=count2+1;

//              tetrahedra which do not contain the energy are rejected)
				if((count2<10000)&&!((bd.eek[bd.tet[0][itet]]<ee)&&(bd.eek[bd.tet[3][itet]]>ee)))
					continue;
				dostet=bd.FASTSURF(ee,itet)/bd.vgt[itet];
				if(((RANDNR()*maxdostet)>dostet)&&(count2<10000))continue;
				nfindfl=false;
//              acception/rejection due to dos in tetrahedron
			  }
//             check band index
			  if(iband!=bd.ibt[itet])
			  {
			  }	
//            no final state found
		      if(count2>=10000)
			  {	 
			     selfscfl=true;
			  }
			  else
			  {	  
//            at this point a tetrahedron has been selected
//            now we are looking for a point on the surface within the tetrahedron
				  selfscfl=false;
			      FINALK(bd);
			  }
		  }
	  }
//end of FINALSTATE
	  return;
}

double Partical::EBHMAXSCRT(DevSimulator &dev,Band &bd)
{
	  //local variables
      double eemax,gammax,rv;

      if(bd.bhmrtfl)
	  {	  
         eemax=frickel*dens/(8.0*bd.meld*dev.eps[SILICON])*0.56875;
	  }	
      else
	  {	
         eemax=frickel*dens/(8.0*bd.melt*dev.eps[SILICON]);
	  }	
      if(eemax<bd.eek[bd.tet[0][itet]])
	  {	
         gammax=EBHSCRT(bd.eek[bd.tet[0][itet]],dev,bd);
	  }
	  else
	  {
		  if(eemax>bd.eek[bd.tet[3][itet]])
		  {
			  gammax=EBHSCRT(bd.eek[bd.tet[3][itet]],dev,bd);
		  }	
		  else
		  {	
			  gammax=EBHSCRT(eemax,dev,bd);
		  }
	  }
      if(bd.bhmrtfl)
	  {
		  gammax=1.01*gammax;
	  }	
      //EBHMAXSCRT=MAX(gammax,1.0/scrt0);

	  //end of EBHMAXSCRT
      rv=MAX(gammax,1.0/scrt0);
	  return rv;
}

void Partical::ESCTR(DevSimulator &dev,Band &bd)
{
//Purpose:   performs electron scattering

	  //local variables
      int iscat,ibold;
      double drand,dscat;//,EIIRATE;
      double sca[MSCPRE*NBE],scasum;
	  int isigh;

	  //only electrons are allowed in this routine
      if(ipt!=PELEC)
	  {	
      }	
	  //save inital band
      ibold=iband;

	  //scattering statistics
      ntotp[iband]+=1;
      dtotp[iband]+=fabs(pc)/bd.gamtet[itet];

	  //just before value for impact ionization current
      if(jbscfl)
	  {
		  dev.iicurjb[PELEC]=dev.iicurjb[PELEC]-pc*bd.EIIRATE(ee)/bd.gamtet[itet];
	  }	
	  //first check for fictious scattering
      scasum=bd.CALSCATTSUM(ee,iband);
      if(scasum>bd.gamtet[itet])
	  {
         drand=RANDNR()*scasum;
	  }	
      else
      {
		  drand=RANDNR()*bd.gamtet[itet];
	  }	
      if(scasum<drand)
	  {
         selfscfl=true;
	  }	
      else
	  {	
         selfscfl=false;

		 //calculate scattering-rate for all prozesses at energy ee
         bd.CALSCATTE(ee,sca,iband);

		 //select the scattering prozess with AcRejection
         iscat=-1;
         dscat=0.0;
		 //90      
		 while(drand>dscat)	
         {
			 iscat=iscat+1;
			 if(iscat>=(bd.scpre*bd.nband[PELEC]))
			{
				selfscfl=true;
				break;
			}
			dscat=dscat+sca[iscat];
		 }
		 if(!selfscfl)
		 {
		 //iscat is number of chosen scatt.pro.
		 //II
			if((fmod(iscat,bd.scpre)+1)==bd.scpre)
			{	
				ESCATII(iscat,dev,bd);
			} 
			//Phonon
			else 
			{
			 if(bd.jacophfl)
             {
				 EPSCAT(iscat,bd);
			 }
			 else
			 {	
				 if(bd.fiscphfl)
				 {
					  EFPSCAT(iscat,bd);
				 }	
				 else
				 {
					 cout<<"error Warning: ESCTR: Specify phonon type";exit(0);
				 }
			 }	
			}
		 }
      }
	  //statistic of scattering-process
      if(selfscfl)
	  {	
         nslfp[ibold]=nslfp[ibold]+1;
	  }	
      else
      {
		  if((fmod(iscat,bd.scpre)+1)==bd.scpre)
          {
			  nreaii[ibold]+=1;
		  }
          else
		  {	
              nreap[ibold]= nreap[ibold]+1;
          }	
         nsctype[iscat][ibold-bd.bandof[ipt]]=nsctype[iscat][ibold-bd.bandof[ipt]]+1;
      }	
//    end of ESCTR
	  return;
}

void Partical::HSCTR(DevSimulator &dev,Band &bd)
{
//Purpose:   performs hole scattering

	  //local variables
      int iscat,ibold;
      double drand,dscat; 
      double sca[MSCPRH*NBH], scasum;

	  //only holes are allowed in this routine
      if(ipt!=PHOLE)
      {
		  cout<<"error HSCTR: Wrong particle type";exit(0);
      }	
	  //save intial band
      ibold=iband;

	  //scattering statistics
      ntotp[iband]=ntotp[iband]+1;
      dtotp[iband]=dtotp[iband]+fabs(pc)/bd.gamtet[itet];

	  //just before value for impact ionization current
      if(jbscfl) 
	  {	  
		  dev.iicurjb[PHOLE]=dev.iicurjb[PHOLE]+pc*bd.HIIRATE(ee)/bd.gamtet[itet];
	  }	

	  //first check for fictious scattering
      scasum=bd.CALSCATTSUM(ee,iband);
      if(scasum>bd.gamtet[itet])
	  {
         drand=RANDNR()*scasum;
	  }	
      else
	  {	
         drand=RANDNR()*bd.gamtet[itet];
      }
      if(scasum<drand)
      {
		 selfscfl =true;
	  }	
      else
      {
         selfscfl =false;

	     //calculate scattering-rate for all prozesses at energy ee
         bd.CALSCATTH(ee,sca,iband);

	     //select the scattering prozess with AcRejection
         iscat=-1;
         dscat=0.0;
		 while(drand>dscat)
		 {   
			 iscat=iscat+1;
			 if(iscat>=(bd.scprh*bd.nband[PHOLE]))
			 {
				 selfscfl =true;
				 break;
			 }
			 dscat=dscat+sca[iscat];
		 }

		 //iscat is number of chosen scatt.pro.
         if(!selfscfl)
		 {
			if((fmod(iscat,bd.scprh)+1)==bd.scprh)
			{
				 HSCATII(iscat,dev,bd);
			 }	
		     else
			 {
				 HPSCAT(iscat,bd);
			 }
		 }
	  }	
	
	  //statistic of scattering-process
      if(selfscfl)
	  {
		   nslfp[ibold]=nslfp[ibold]+1;
	  }
      else
	  {
		   if((fmod(iscat,bd.scprh)+1)==bd.scprh)
		   {
				nreaii[ibold]=nreaii[ibold]+1;
		   }	
           else
		   {
				nreap[ibold]=nreap[ibold]+1;
		   }	
	    
		   nsctyph[iscat][ibold-bd.bandof[ipt]]+=1;
      }	
//    end of HSCTR
	  return;
}

void Partical::OESCTR(Band &bd) 
{
//     Purpose:   performs oxide electron scattering
	  //local variables
      int iscat,ibold;
      double drand,dscat; 
      double sca[MSCPROE*NBOE],scasum;

      //only oxide electrons are allowed in this routine
	  if(ipt!=POXEL)
	  {
		  cout<<"error OESCTR: Wrong particle type";exit(0);
	  }	
      
      //scattering statistics
	  ntotp[iband]=ntotp[iband]+1;	
      //dtotp(iband) = dtotp(iband)+ ABS(pc)/bd.gamtet[itet]
	  dtotp[iband]=dtotp[iband]+fabs(pc)/bd.gamtet[itet];
      //save band index
      ibold=iband;

      //first check for fictious scattering
      scasum=bd.CALSCATTSUM(ee,iband);
	  if(scasum>bd.gamtet[itet])
      {
          drand = RANDNR()*scasum;
	  }
      else
	  {	
          drand = RANDNR()*bd.gamtet[itet];
	  }	
	  if(scasum<drand)
      {
		  selfscfl=true;
      }
      else
      {
         selfscfl=false;
	  
	  //calculate scattering-rate for all prozesses at energy ee
	      bd.CALSCATTOE(ee,sca,iband);
	
		  //select the scattering prozess with AcRejection
	      iscat=-1;
		  dscat=0.0;
		  while(drand>dscat)
		  {	
			  iscat++;
		 	  if(iscat>=bd.scproe*bd.nband[POXEL])
			  {
				selfscfl=true;
				break;
			 }
			else
			dscat=dscat+sca[iscat];
		}	
	  //iscat is number of chosen scatt.pro.
		if(!selfscfl)OEPSCAT(iscat,bd);
      }

      //statistic of scattering-process
      if(selfscfl)
      {
		  nslfp[ibold]=nslfp[ibold]+1;
      }
	  else
	  {	
          nreap[ibold]=nreap[ibold]+1;
          nsctypoe[iscat][ibold-bd.bandof[ipt]]+=1;
      }
	  return;
}

void Partical::CALSCATTSURFACE(DevSimulator &dev,Band &bd,int itype)
{
	double zzz,eebeta;
	double f1,f2,f3,f4,f5;
	double effmass=0;
//	xminss=xminss*1e-6/spr0;
//	xmaxss=xmaxss*1e-6/spr0;
//	yminss=yminss*1e-6/spr0;
//	ymaxss=ymaxss*1e-6/spr0;
	surfaceposition1=surfaceposition1*1e-6/spr0;
	surfaceposition2=surfaceposition2*1e-6/spr0;
	if(material==0)
	{
		if(itype==PELEC)
		{
			effmass=0.286;
		}
		else if(itype==PHOLE)
		{
			effmass=0.8;
		}
		else
		{
			effmass=0.286;
		}
	}
	if(material==1)
	{
		if(itype==PELEC)
		{
			effmass=0.12;
		}
		else if(itype==PHOLE)
		{
			effmass=0.34;
		}
		else
		{
			effmass=0.12;
		}
	}

	if(fabs(xf)<1e-20)	xf=1e-20;
//  surface roughness scattering rate

	srrv=pow(dev.eps[SILICON]/(dev.eps[SILICON]+dev.eps[OXIDE]),2)*2*PI*effmass*EM*EC*EC
		/pow(hq0*EC,3)*pow(delta*ail*1e-18,2)*pow(fabs(xf*field0),2);

	srrv=srrv/scrt0;
	if(!srfl)	srrv=0;

//	surface phonon scattering rate
	f1=pow((T0/300.0),2)*(8.6173468e-5)*T0/(fabs(xf*field0)*0.01);
	f2=pow(hq0,2.0/3.0)/pow((12.0*effmass*(9.1095e-35)/EC),1.0/3.0)
		*pow(fabs(xf*field0)*0.01,-1.0/3.0);
	f3=0.7*T0/300.0+0.3*pow(dens*conc0*1e-6/1e17,-0.6)*300.0/T0;
	f4=(1.5*f1+2.857703*f2)/f3;
	sprv=effmass*(9.109534e-35)/(1.6021892e-19)*(8.6173468e-5)*T0*144
		/(pow(6.5821731e-16,3)*(2330e-10)/(EC)*pow(843302,2)*f4)
		+1.0/(effmass*(9.109534e-35)/(1.6021892e-19)*3600*pow(T0/300.0,2.24));
	sprv=sprv/scrt0;
	if(!spfl)	sprv=0;

//	surface impurity scattering rate
	zzz=fabs(xr-surfaceposition1)*spr0/1e-6*1e3;
	f1=0.5+10/(1+exp((1e18-dope*conc0/1e6)/(0.08*1e18)));
	f2=-f1*pow(zzz,6)/((pow(zzz,4)+pow(1,4))*pow(1,4));
	f3=1-exp(f2);
	if(fabs(f3)<1e-20)
	{
		f3=1e-20;
	}
    if((ee*eV0<0.120)&&(dope*conc0>=1e20))
	{
		eebeta=frickel*dens/(2.0*bd.meld*dev.eps[SILICON]);
		sirv=dope*bd.meld*sqrt(2.0*bd.meld*ee)/(PI*pow(frickel*dens/f3,2.0)
				*(1.0+4.0*ee*bd.melt/(eebeta*bd.meld)));
	}
	else
	{
		sirv=0;
	}
	if(!scfl)	sirv=0;

	ssrv=srrv+sprv+sirv;

	return;
}

void Partical::ESSCRT(DevSimulator &dev,Band &bd)
{
      int iscat;
      double drand,dscat;
	  double scasum;
	  double sscat[3];

	  scasum=ssrv;
	  sscat[0]=srrv;
	  sscat[1]=sprv;
	  sscat[2]=sirv;

	  if(scasum>surfacemax)
	  {
		  drand=RANDNR()*scasum;
	  }
	  else
	  {
		  drand=RANDNR()*surfacemax;
	  }
	  if(scasum<drand)
	  {
         selfscfl=true;
	  }
      else
	  {	
         selfscfl=false;
		 iscat=-1;
         dscat=0.0;
		 while(drand>dscat)
         {
			 iscat=iscat+1;
			 dscat=dscat+sscat[iscat];
		 }
		 if(!selfscfl)
		 {
			 if(iscat==0)
			 {
				 dev.srscatter++;
				 DIFFUSEBC(bd);
			 }
			 else if(iscat==1)
			 {
				 dev.spscatter++;
				 ESCTR(dev,bd);
			 }
			 else if(iscat==2)
			 {
				 dev.siscatter++;
				 EBHSCAT(dev,bd);
			 }
			 else
			 {
				 cout<<"Wrong Electron surface scattering type!"<<endl;
				 exit(0);
			 }
		 }
	  }
	  return;
}

void Partical::HSSCRT(DevSimulator &dev,Band &bd)
{
      int iscat;
      double drand,dscat;
	  double scasum;
	  double sscat[3];

	  scasum=ssrv;
	  sscat[0]=srrv;
	  sscat[1]=sprv;
	  sscat[2]=sirv;

	  if(scasum>surfacemax)
	  {
		  drand=RANDNR()*scasum;
	  }
	  else
	  {
		  drand=RANDNR()*surfacemax;
	  }
	  if(scasum<drand)
	  {
         selfscfl=true;
	  }
      else
	  {	
         selfscfl=false;
		 iscat=-1;
         dscat=0.0;
		 while(drand>dscat)
         {
			 iscat=iscat+1;
			 dscat=dscat+sscat[iscat];
		 }
		 if(!selfscfl)
		 {
			 if(iscat==0)
			 {
				 dev.srscatter++;
				 DIFFUSEBC(bd);
			 }
			 else if(iscat==1)
			 {
				 dev.spscatter++;
				 HSCTR(dev,bd);
			 }
			 else if(iscat==2)
			 {
//				 dev.siscatter++;
//				 EBHSCAT(dev,bd);
			 }
			 else
			 {
				 cout<<"Wrong Electron surface scattering type!"<<endl;
				 exit(0);
			 }
		 }
	  }
	  return;
}

void Partical::OESSCRT(DevSimulator &dev,Band &bd)
{

}

double Partical::EBHSCRT(double eel,DevSimulator &dev,Band &bd)
{
	  //local variables
      double eebeta,rvscrt;
	  if(eel<0)	eel=1e-10;

      if((eel*eV0<minelectronenergy)&&(dope*conc0>=mindopforcoulomb*1e6))//1e22
//	  if((eel*eV0<0.120)&&(dope*conc0>=1e22))
	  {
         eebeta=frickel*dens/(2.0*bd.meld*dev.eps[SILICON]);

         if(bd.bhmrtfl)
		 {	 
            eel=MAX(eel,(1e-5));
			/*
            EBHSCRT=dope/(16.0*PI*(dev.eps[SILICON])**2 
                           * SQRT(2.0*bd.meld*eel**3))
                   * (LOG((eebeta+4.0*eel)/eebeta)
                   - 4.0*eel/(eebeta+4.0*eel));
			*/
			rvscrt=dope/(16.0*PI*pow(dev.eps[SILICON],2.0) 
                   *sqrt(2.0*bd.meld*pow(eel,3.0)))
                   * (log((eebeta+4.0*eel)/eebeta)
                   - 4.0*eel/(eebeta+4.0*eel));
		 }
         else
		 {
			/* 
			EBHSCRT=dope*bd.meld*sqrt(2.0*bd.meld*eel)/(PI*(frickel*dens)**2 
				*(1.0+4.0*eel*bd.melt/(eebeta*bd.meld)));
			*/
			rvscrt=dope*bd.meld*sqrt(2.0*bd.meld*eel)/(PI*pow(frickel*dens,2.0) 
				*(1.0+4.0*eel*bd.melt/(eebeta*bd.meld)));
		 }
	  }
      else
	  {	
         //EBHSCRT=0.0;
		  rvscrt=0.0;
      }
//    end of EBHSCRT
	  return rvscrt;
}

void Partical::EBHSCAT(DevSimulator &dev,Band &bd)
{
	  //Purpose:   performs electron impurity scattering
	  //local variables
      double xkl,ykl,zkl,ca,cb;
      double xkll,ykll,zkll;
      double xklll,yklll,zklll;
      double cosphi,sinphi,costheta,sintheta;
      double cospr,sinpr,costr,sintr;
      double eebeta,pr,alfa,betaq;
      double dn,dnmax;

	  //Transform k-vector into irreducible wedge
      INWEDGE (xkl,ykl,zkl,xk,yk,zk,bd.matsym,isym);

	  //calculate k-vector relative to band minimum
	  if(material==0)
	  {
		  xkl=xkl-0.85*a0pi;
	  }
	  else if(material==1)
	  {
		xkl=xkl-0.5*a0pi;
		ykl=ykl-0.5*a0pi;
		zkl=zkl-0.5*a0pi;
	  }

	  //apply BH-transformation
      xkl=xkl*sqrt(bd.meld/bd.mell);
      ykl=ykl*sqrt(bd.meld/bd.melt);
      zkl=zkl*sqrt(bd.meld/bd.melt);

	  //screening energy
      betaq=frickel*dens/dev.eps[SILICON];
      eebeta=betaq/(2.0*bd.meld);

	  //cal. angle between old and new k-vector
      alfa=eebeta*bd.meld/(2.0*ee*bd.melt);
      pr=RANDNR()/(1.0+0.5*alfa);
      costr=1.0-alfa*pr/(1.0-pr);
	  //???
      sintr=sqrt(1.0-costr*costr);

	  //cal. angle of the new k-vector in the plane perpendicular to 
	  //the old k-vector
      pr=2.0*PI*RANDNR();
      cospr=cos(pr);
      sinpr=sin(pr);
      
	  //cal. new k-vector in a coordinate frame, where the old k-vector points 
	  //into the positve z-direction
      ca=sqrt(xkl*xkl+ykl*ykl);
      cb=sqrt(ca*ca+zkl*zkl);
      xkll=cb*sintr*cospr;
	  ykll=cb*sintr*sinpr;
      zkll=cb*costr;

	  //transform new k-vector into the BH-transformed coordinate frame
      cosphi=xkl/ca;
      sinphi=ykl/ca;

      costheta=zkl/cb;
      sintheta=ca/cb;

      xklll=cosphi*costheta*xkll-sinphi*ykll+cosphi*sintheta*zkll;           
      yklll=sinphi*costheta*xkll+cosphi*ykll+sinphi*sintheta*zkll;
      zklll=-sintheta*xkll+costheta*zkll;
      
	  //Acception Rejection
      //dn=(betaq+bd.mell/bd.meld*(xkl-xklll)**2+bd.melt/bd.meld*
	  //  ((ykl-yklll)**2+(zkl-zklll)**2))**(-2);
	  dn=pow(betaq+bd.mell/bd.meld*(xkl-xklll)*(xkl-xklll)+bd.melt/bd.meld*
	    ((ykl-yklll)*(ykl-yklll)+(zkl-zklll)*(zkl-zklll)),-2);	
      //dnmax=(betaq+bd.melt/bd.meld* (xkl-xklll)**2
		//  +bd.melt/bd.meld*((ykl-yklll)**2+(zkl-zklll)**2))**(-2);
	  dnmax=pow(betaq+bd.melt/bd.meld* (xkl-xklll)*(xkl-xklll)
		  +bd.melt/bd.meld*((ykl-yklll)*(ykl-yklll)+(zkl-zklll)*(zkl-zklll)),-2);
      if(dn>dnmax)
	  {
		  cout<<"error EBHSCAT: AcRe - Error";exit(0);

	  }	
      if(dn>RANDNR()*dnmax)
	  {	  

		 //inverse BH-transformation
         xkl=xklll*sqrt(bd.mell/bd.meld);
         ykl=yklll*sqrt(bd.melt/bd.meld);
         zkl=zklll*sqrt(bd.melt/bd.meld);

		 //move origin of the coordinate frame back into the Gamma point
		 if(material==0)
		 {
			 xkl=xkl+0.85*a0pi;
		 }
		 else if(material==1)
		 {
			 xkl=xkl+0.5*a0pi;
			 ykl=ykl+0.5*a0pi;
			 zkl=zkl+0.5*a0pi;
		 }

		 //transform vector back into the BZ
         bd.OUTWEDGE(xk,yk,zk,xkl,ykl,zkl,bd.matsym,isym);

		 //find the new symmetry operation and tetrahedron index
         bd.CONFBZ(xk,yk,zk);
         GETSYM(bd);

		 //cal. new energy (since the assumed elliptical parabolic band strucutre
		 //only approximates the real band strucutre)
         GETGVE(bd);
         EENER(bd);

		 //real scattering event
         selfscfl=false;

		 //statistic of scattering-process
         nreabh[iband]=nreabh[iband]+1;
	  }	
      else
	  {
	     //self scattering event
         nslfbh[iband]=nslfbh[iband]+1;
         selfscfl=true;

	  }	
//    end of EBHSCAT
	  return;
}

void Partical::EBHSCTR(double gamimpmax,DevSimulator &dev,Band &bd)
{
	  //Purpose:   performs electron impurity scattering
	  //local variables
      int itetold,isymold;
      double gamimp;

	  //only electrons are allowed in this routine
      if(ipt!=PELEC)
	  {
	  }
	  //scattering statistics
      ntotbh[iband]+=1;

	  //check maximum scattering rate
      gamimp=EBHSCRT(ee,dev,bd);
      if(gamimp>(gamimpmax+1e-12))
	  {
		  cout<<"error EBHSCTR: Variable Gamma wrong";exit(0);
	  }	
	  //first check for fictious scattering
      if(gamimp<RANDNR()*gamimpmax)
	  {	  
         nslfbh[iband]+=1;
         selfscfl=true;
	  }	
      else
	  {
		 //Chose final state
         if(bd.bhmrtfl)
		 {		 
            isymold=isym;
            itetold=itet;
            TALIND(1);
            FINALSTATE(bd);
            if(selfscfl)
			{		
               isym=isymold;
               itet=itetold;
               nslfbh[iband]+=1;
			}
            else
			{
               GETGVE(bd);
               nreabh[iband]+=1;
			}
		 }
         else
		 {	
            EBHSCAT(dev,bd); 
		 }
	 }
	 //end of EBHSCTR
     return;
}

double Partical::FEKLOEM(double eet,Band &bd)
{
//	local variables
  int i,j,k,it,ei,ej;
  double ie,edeta,ran,ksum,efinal;
  if(material==0)
  {
	ie=eet/bd.templog;
	if(ie<1.0e0)
	{
		efinal=1.0e-10;
		return efinal;
	}
	if(ie<3.05e0)
	{
		it=int((ie-1.0e0)*10e0+0.5e0);
		if(it<1)	it=1;
		ei=it;
		edeta=1e0;
	}
	else if(ie<5.625e0)
	{
		it=int((ie-3.0e0)*4e0+0.5e0);
		if(it<1)	it=1;
		ei=20+it;
		edeta=2e0;
	}
	else if(ie<13.25e0)
	{
		it=int((ie-5.5e0)*2e0+0.5e0);
		if(it<1)	it=1;
		ei=30+it;
		edeta=4e0;
	}
	else if(ie<28.5e0)
	{
		it=int(ie-13e0+0.5e0);
		if(it<1)	it=1;
		ei=45+it;
		edeta=8e0;
	}
	else
	{
		it=int((ie-28.0e0)+0.5e0);
		if(it<1)	it=1;
		if(it>60)	it=60;
		ei=60+it;
		edeta=12e0;
	}
	ran=RANDNR();
	i=0;
	ksum=0.0e0;
	while(ksum<ran)
	{
		i++;
		if(i==200)
		{
			ksum=1.1e0;
			i=1;
		}
		ksum+=bd.kloem[i-1][ei-1];
	}
	j=fmod(i,2);
	k=int(i/2);
	if(j==1)	efinal=eet-bd.templog*(1.0e0-0.01*edeta*k);
	else	efinal=eet-bd.templog*(1.0e0+0.01*edeta*k);
  }
  else if(material==1)
  {
	ie=eet/bd.templog;
	if(ie<1.0e0)
	{
		efinal=1.0e-10;
		return efinal;
	}
	it=int((ie-1e0)*10e0+0.5e0);
	if(it<1)	it=1;
	ei=it;
	edeta=1.0e0;
	ran=RANDNR();
	i=0;
	ksum=0.0e0;
	while(ksum<ran)
	{
		i++;
		if(i==100)
		{
			ksum=1.1e0;
			i=1;
		}
		ksum+=bd.kloem[i-1][ei-1];
	}
	j=fmod(i,2);
	k=int(i/2);
	if(j==1)	efinal=eet-bd.templog*(1.0e0-0.02*edeta*k);
	else	efinal=eet-bd.templog*(1.0e0+0.02*edeta*k);
  }
	return efinal;
}

double Partical::FEKLOAB(double eet,Band &bd)
{
//	local variables
	int i,j,k,it,ei,ej;
	double ie,edeta,ran,ksum,efinal;
	it=0;
  if(material==0)
  {
	ie=eet/bd.templog;
	if(ie<1.05e0)
	{
		it=int(it*10e0+0.5e0);
		if(it<1)	it=1;
		ei=it;
		edeta=1e0;
	}
	else if(ie<3.625e0)
	{
		it=int((ie-1.0e0)*4e0+0.5e0);
		if(it<1)	it=1;
		ei=10+it;
		edeta=2e0;
	}
	else if(ie<11.25e0)
	{
		it=int((ie-3.5e0)*2e0+0.5e0);
		if(it<1)	it=1;
		ei=20+it;
		edeta=4e0;
	}
	else if(ie<26.5e0)
	{
		it=int((ie-11e0)+0.5e0);
		if(it<1)	it=1;
		ei=35+it;
		edeta=8e0;
	}
	else
	{
		it=int((ie-26.0e0)+0.5e0);
		if(it<1)	it=1;
		if(it>60)	it=60;
		ei=50+it;
		edeta=12e0;
	}
	ran=RANDNR();
	i=0;
	ksum=0e0;
	while(ksum<ran)
	{
		i++;
		ksum+=bd.kloab[i-1][ei-1];
		if(i==200)
		{
			ksum=1.1e0;
			i=1;
		}
	}
	j=fmod(i,2);
	k=int(i/2);
	if(j==1)	efinal=eet+bd.templog*(1.0e0+0.01*edeta*k);
	else	efinal=eet+bd.templog*(1.0e0-0.01*edeta*k);
  }
  else if(material==1)
  {
	ie=eet/bd.templaf;
	if(ie<1.05e0)
	{
		it=int(it*10e0+0.5e0);
		if(it<1)	it=1;
		ei=it;
		edeta=0.25e0;
	}
	else if(ie<3.625e0)
	{
		it=int((ie-1.0e0)*4e0+0.5e0);
		if(it<1)	it=1;
		ei=10+it;
		edeta=0.25e0;
	}
	else if(ie<11.25e0)
	{
		it=int((ie-3.5e0)*2e0+0.5e0);
		if(it<1)	it=1;
		ei=20+it;
		edeta=0.5e0;
	}
	else
	{
		it=int((ie-11.0e0)+0.5e0);
		if(it<1)	it=1;
		if(it>60)	it=115;
		ei=35+it;
		edeta=1e0;
	}
	ran=RANDNR();
	i=0;
	ksum=0e0;
	while(ksum<ran)
	{
		i++;
		ksum+=bd.kloab[i-1][ei-1];
		if(i==100)
		{
			ksum=1.1e0;
			i=1;
		}
	}
	j=fmod(i,2);
	k=int(i/2);
	if(j==1)	efinal=eet+bd.templaf*(1.0e0+0.01*edeta*k);
	else	efinal=eet+bd.templaf*(1.0e0-0.01*edeta*k);
  }
	return efinal;
}

double Partical::FEKLAAB(double eet,Band &bd)
{
//	local variables
	int i,j,k,it,ei,ej;
	double ie,edeta,ran,ksum,efinal;
	it=0;
  if(material==0)
  {
	ie=eet/bd.templaf;
	if(ie<1.05e0)
	{
		it=int(it*10e0+0.5e0);
		if(it<1)	it=1;
		ei=it;
		edeta=0.25e0;
	}
	else if(ie<3.625e0)
	{
		it=int((ie-1.0e0)*4e0+0.5e0);
		if(it<1)	it=1;
		ei=10+it;
		edeta=0.25e0;
	}
	else if(ie<11.25e0)
	{
		it=int((ie-3.5e0)*2e0+0.5e0);
		if(it<1)	it=1;
		ei=20+it;
		edeta=0.5e0;
	}
	else
	{
		it=int((ie-11.0e0)+0.5e0);
		if(it<1)	it=1;
		if(it>60)	it=115;
		ei=35+it;
		edeta=1.0e0;
	}
	ran=RANDNR();
	i=0;
	ksum=0e0;
	while(ksum<ran)
	{
		i++;
		ksum+=bd.klaab[i-1][ei-1];
		if(i==50)
		{
			ksum=1.1e0;
			i=1;
		}
	}
	j=fmod(i,2);
	k=int(i/2);
	if(j==1)	efinal=eet+bd.templaf*(1.0e0+0.04*edeta*k);
	else	efinal=eet+bd.templaf*(1.0e0-0.04*edeta*k);
	if(efinal<=0e0)	efinal=(1.0e-5)*bd.templaf;
  }
  else if(material==1)
  {
	ie=eet/bd.templaf;
	if(ie<1.05e0)
	{
		it=int(it*10e0+0.5e0);
		if(it<1)	it=1;
		ei=it;
		edeta=0.25e0;
	}
	else if(ie<3.625e0)
	{
		it=int((ie-1.0e0)*4e0+0.5e0);
		if(it<1)	it=1;
		ei=10+it;
		edeta=0.25e0;
	}
	else if(ie<11.25e0)
	{
		it=int((ie-3.5e0)*2e0+0.5e0);
		if(it<1)	it=1;
		ei=20+it;
		edeta=0.5e0;
	}
	else
	{
		it=int((ie-11.0e0)+0.5e0);
		if(it<1)	it=1;
		if(it>60)	it=115;
		ei=35+it;
		edeta=1.0e0;
	}
	ran=RANDNR();
	i=0;
	ksum=0e0;
	while(ksum<ran)
	{
		i++;
		ksum+=bd.klaab[i-1][ei-1];
		if(i==100)
		{
			ksum=1.1e0;
			i=1;
		}
	}
	j=fmod(i,2);
	k=int(i/2);
	if(j==1)	efinal=eet+bd.templaf*(1.0e0+0.02*edeta*k);
	else	efinal=eet+bd.templaf*(1.0e0-0.02*edeta*k);
	if(efinal<=0e0)	efinal=(1.0e-5)*bd.templaf;
  }
	return efinal;
}

double Partical::FEKLAEM(double eet,Band &bd)
{
//	local variables
	int i,j,k,it,ei,ej;
	double ie,edeta,ran,ksum,efinal;
	it=0;
  if(material==0)
  {
	ie=eet/bd.templaf;
	if(ie<1e0)
	{
		efinal=-1e-10;
		return efinal;
	}
	if(ie<3.05e0)
	{
		it=int((it-1)*10e0+0.5e0);
		if(it<1)	it=1;
		ei=it;
		edeta=0.25e0;
	}
	else if(ie<5.625e0)
	{
		it=int((ie-3.0e0)*4e0+0.5e0);
		if(it<1)	it=1;
		ei=20+it;
		edeta=0.5e0;
	}
	else if(ie<13.25e0)
	{
		it=int((ie-5.5e0)*2e0+0.5e0);
		if(it<1)	it=1;
		ei=30+it;
		edeta=0.5e0;
	}
	else
	{
		it=int((ie-13.0e0)+0.5e0);
		if(it<1)	it=1;
		if(it>115)	it=115;
		ei=45+it;
		edeta=1.0e0;
	}
	ran=RANDNR();
	i=0;
	ksum=0e0;
	while(ksum<ran)
	{
		i++;
		ksum+=bd.klaab[i-1][ei-1];
		if(i==50)
		{
			ksum=1.1e0;
			i=1;
		}
	}
	j=fmod(i,2);
	k=int(i/2);
	if(j==1)	efinal=eet+bd.templaf*(1.0e0+0.04*edeta*k);
	else	efinal=eet+bd.templaf*(1.0e0-0.04*edeta*k);
	if(efinal<=0e0)	efinal=(1.0e-5)*bd.templaf;
  }
  if(material==1)
  {
	ie=eet/bd.templaf;
	if(ie<0e0)
	{
		efinal=-1e-10;
		return efinal;
	}
	it=int((ie-1e0)*10e0+0.5e0);
	if(it<1)	it=1;
	ei=it;
	edeta=1e0;
	ran=RANDNR();
	i=0;
	ksum=0e0;
	while(ksum<ran)
	{
		i++;
		if(i==100)
		{
			ksum=1.1e0;
			i=1;
		}
		ksum+=bd.klaem[i-1][ei-1];
	}
	j=fmod(i,2);
	k=int(i/2);
	if(j==1)	efinal=eet-bd.templaf*(1.0e0-0.02*edeta*k);
	else	efinal=eet-bd.templaf*(1.0e0+0.02*edeta*k);
	if(efinal<=0e0)	efinal=(1.0e-5)*bd.templaf;
  }
	return efinal;
}

double Partical::FEKTOAB(double eet,Band &bd)
{
//	local variables
  int i,j,k,it,ei,ej;
  double ie,edeta,ran,ksum,efinal;
  it=0;
  if(material==0)
  {
	ie=eet/bd.temptof;
	if(ie<1.05e0)
	{
		it=int(it*10e0+0.5e0);
		if(it<1)	it=1;
		ei=it;
		edeta=0.25e0;
	}
	else if(ie<3.625e0)
	{
		it=int((ie-1.0e0)*4e0+0.5e0);
		if(it<1)	it=1;
		ei=10+it;
		edeta=0.25e0;
	}
	else if(ie<11.25e0)
	{
		it=int((ie-3.5e0)*2e0+0.5e0);
		if(it<1)	it=1;
		ei=20+it;
		edeta=0.5e0;
	}
	else
	{
		it=int((ie-11.0e0)+0.5e0);
		if(it<1)	it=1;
		if(it>60)	it=115;
		ei=35+it;
		edeta=1.0e0;
	}
	ran=RANDNR();
	i=0;
	ksum=0e0;
	while(ksum<ran)
	{
		i++;
		ksum+=bd.ktoab[i-1][ei-1];
		if(i==50)
		{
			ksum=1.1e0;
			i=1;
		}
	}
	j=fmod(i,2);
	k=int(i/2);
	if(j==1)	efinal=eet+bd.temptof*(1.0e0+0.04*edeta*k);
	else	efinal=eet+bd.temptof*(1.0e0-0.04*edeta*k);
  }
  else if(material==1)
  {
	ie=eet/bd.temptof;
	if(ie<1.05e0)
	{
		it=int(it*10e0+0.5e0);
		if(it<1)	it=1;
		ei=it;
		edeta=0.25e0;
	}
	else if(ie<3.625e0)
	{
		it=int((ie-1.0e0)*4e0+0.5e0);
		if(it<1)	it=1;
		ei=10+it;
		edeta=0.25e0;
	}
	else if(ie<11.25e0)
	{
		it=int((ie-3.5e0)*2e0+0.5e0);
		if(it<1)	it=1;
		ei=20+it;
		edeta=0.5e0;
	}
	else
	{
		it=int((ie-11.0e0)+0.5e0);
		if(it<1)	it=1;
		if(it>60)	it=115;
		ei=35+it;
		edeta=1.0e0;
	}
	ran=RANDNR();
	i=0;
	ksum=0e0;
	while(ksum<ran)
	{
		i++;
		ksum+=bd.ktoab[i-1][ei-1];
		if(i==50)
		{
			ksum=1.1e0;
			i=1;
		}
	}
	j=fmod(i,2);
	k=int(i/2);
	if(j==1)	efinal=eet+bd.templaf*(1.0e0+0.04*edeta*k);
	else	efinal=eet+bd.templaf*(1.0e0-0.04*edeta*k);
  }
	return efinal;
}

double Partical::FEKTOEM(double eet,Band &bd)
{
//	local variables
	int i,j,k,it,ei,ej;
	double ie,edeta,ran,ksum,efinal;
  if(material==0)
  {
	ie=eet/bd.temptof;
	if(ie<1.0e0)
	{
		efinal=1.0e-10;
		return efinal;
	}
	else if(ie<3.05e0)
	{
		it=int((ie-1.0e0)*10e0+0.5e0);
		if(it<1)	it=1;
		ei=10+it;
		edeta=0.25e0;
	}
	else if(ie<5.625e0)
	{
		it=int((ie-3.0e0)*4e0+0.5e0);
		if(it<1)	it=1;
		ei=20+it;
		edeta=0.5e0;
	}
	else if(ie<13.25e0)
	{
		it=int((ie-5.5e0)*2e0+0.5e0);
		if(it<1)	it=1;
		ei=30+it;
		edeta=0.5e0;
	}
	else
	{
		it=int((ie-13.0e0)+0.5e0);
		if(it<1)	it=1;
		if(it>115)	it=115;
		ei=45+it;
		edeta=1.0e0;
	}
	ran=RANDNR();
	i=0;
	ksum=0e0;
	while(ksum<ran)
	{
		i++;
		ksum+=bd.ktoab[i-1][ei-1];
		if(i==50)
		{
			ksum=1.1e0;
			i=1;
		}
	}
	j=fmod(i,2);
	k=int(i/2);
	if(j==1)	efinal=eet-bd.temptof*(1.0e0-0.04*edeta*k);
	else	efinal=eet-bd.temptof*(1.0e0+0.04*edeta*k);
  }
  else if(material==1)
  {
	ie=eet/bd.temptof;
	if(ie<1.0e0)
	{
		efinal=1.0e-10;
		return efinal;
	}
	else if(ie<3.05e0)
	{
		it=int((ie-1.0e0)*10e0+0.5e0);
		if(it<1)	it=1;
		ei=10+it;
		edeta=0.25e0;
	}
	else if(ie<5.625e0)
	{
		it=int((ie-3.0e0)*4e0+0.5e0);
		if(it<1)	it=1;
		ei=20+it;
		edeta=0.5e0;
	}
	else if(ie<13.25e0)
	{
		it=int((ie-5.5e0)*2e0+0.5e0);
		if(it<1)	it=1;
		ei=30+it;
		edeta=0.5e0;
	}
	else
	{
		it=int((ie-13.0e0)+0.5e0);
		if(it<1)	it=1;
		if(it>115)	it=115;
		ei=45+it;
		edeta=1.0e0;
	}
	ran=RANDNR();
	i=0;
	ksum=0e0;
	while(ksum<ran)
	{
		i++;
		ksum+=bd.ktoab[i-1][ei-1];
		if(i==51)
		{
			ksum=1.1e0;
			i=1;
		}
	}
	j=fmod(i,2);
	k=int(i/2);
	if(j==1)	efinal=eet-bd.temptof*(1.0e0-0.04*edeta*k);
	else	efinal=eet-bd.temptof*(1.0e0+0.04*edeta*k);
  }
	return efinal;
}

void Partical::EPSCAT(int iscat,Band &bd)
{
	  //int iscat;

	  //Purpose: Calculation of the electron state just after scattering
      //Parameter :  iscat = number of actual scattering-prozess
	  //local variables
      int ibold,itetold,isymold;
      double eeold,xkold,ykold,zkold;

	  //save particle state
      eeold=ee;
      xkold=xk;
      ykold=yk;
      zkold=zk;
      ibold=iband;
      itetold=itet;
      isymold=isym;

	  //calculate final band and type of scattering
      iband=iscat/bd.scpre+bd.bandof[ipt];
      //GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14),MOD(iscat-1,bd.scpre)+1
	  int ntemp=fmod(iscat,bd.scpre);

	if(!quantumeffectcbfl)
	{
		if(material==0)
		{
			switch(ntemp)
			{
				case 0:ee=ee;TALIND(1);break;//elastic intravalley
				case 1:ee=ee+bd.temptag;TALIND(2);break;//transversal acoustic g-Phonon
				case 2:ee=ee-bd.temptag;TALIND(2);break;//transversal acoustic g-Phonon
				case 3:ee=ee+bd.templag;TALIND(2);break;//longitudinal acoustic g-Phonon
				case 4:ee=ee-bd.templag;TALIND(2);break;//longitudinal acoustic g-Phonon
				case 5:ee=ee+bd.templog;TALIND(2);break;//longitudinal optical g-Phonon
				case 6:ee=ee-bd.templog;TALIND(2);break;//longitudinal optical g-Phonon
				case 7:ee=ee+bd.temptaf;TALIND(3);break;//transversal acoustic f-Phonon
				case 8:ee=ee-bd.temptaf;TALIND(3);break;//transversal acoustic f-Phonon
				case 9:ee=ee+bd.templaf;TALIND(3);break;//longitudinal acoustic f-Phonon
				case 10:ee=ee-bd.templaf;TALIND(3);break;//longitudinal acoustic f-Phonon
				case 11:ee=ee+bd.temptof;TALIND(3);break;//transversal optical f-Phonon
				case 12:ee=ee-bd.temptof;TALIND(3);break;//transversal optical f-Phonon
				case 13:exit(0);break;//Impact Ionization
				default:break;
			}
		}
		else if(material==1)
		{
			switch(ntemp){
				case 0:ee=ee;TALINDL(1,bd);break;//elastic intravalley
				case 1:ee=ee+bd.temptag;TALINDL(2,bd);break;//transversal acoustic g-Phonon
				case 2:ee=ee-bd.temptag;TALINDL(2,bd);break;//transversal acoustic g-Phonon
				case 3:ee=ee+bd.templag;TALINDL(2,bd);break;//longitudinal acoustic g-Phonon
				case 4:ee=ee-bd.templag;TALINDL(2,bd);break;//longitudinal acoustic g-Phonon
				case 5:ee=ee+bd.templog;TALINDL(2,bd);break;//longitudinal optical g-Phonon
				case 6:ee=ee-bd.templog;TALINDL(2,bd);break;//longitudinal optical g-Phonon
				case 7:exit(0);break;//Impact Ionization
				default:break;
			}
		}
	}
	else if(quantumeffectcbfl)
	{
		if(material==0)
		{
			switch(ntemp)
			{
				case 0:ee=ee;TALIND(1);break;//elastic intravalley
				case 1:ee=ee+bd.temptag;TALIND(2);break;//transversal acoustic g-Phonon
				case 2:ee=ee-bd.temptag;TALIND(2);break;//transversal acoustic g-Phonon
				case 3:ee=ee+bd.templag;TALIND(2);break;//longitudinal acoustic g-Phonon
				case 4:ee=ee-bd.templag;TALIND(2);break;//longitudinal acoustic g-Phonon
				case 5:ee=FEKLOAB(ee,bd);TALIND(2);break;//longitudinal optical g-Phonon
				case 6:ee=FEKLOEM(ee,bd);TALIND(2);break;//longitudinal optical g-Phonon
				case 7:ee=ee+bd.temptaf;TALIND(3);break;//transversal acoustic f-Phonon
				case 8:ee=ee-bd.temptaf;TALIND(3);break;//transversal acoustic f-Phonon
				case 9:ee=FEKLAAB(ee,bd);TALIND(3);break;//longitudinal acoustic f-Phonon
				case 10:ee=FEKLAEM(ee,bd);TALIND(3);break;//longitudinal acoustic f-Phonon
				case 11:ee=FEKTOAB(ee,bd);TALIND(3);break;//transversal optical f-Phonon
				case 12:ee=FEKTOEM(ee,bd);TALIND(3);break;//transversal optical f-Phonon
				case 13:exit(0);break;//Impact Ionization
				default:break;
			}
		}
		else if(material==1)
		{
			switch(ntemp)
			{
				case 0:ee=ee;TALINDL(1,bd);break;//elastic intravalley
				case 1:ee=ee+bd.temptag;TALINDL(2,bd);break;//transversal acoustic g-Phonon
				case 2:ee=ee-bd.temptag;TALINDL(2,bd);break;//transversal acoustic g-Phonon
				case 3:ee=FEKLAAB(ee,bd);TALINDL(2,bd);break;//longitudinal acoustic g-Phonon
				case 4:ee=FEKLAEM(ee,bd);TALINDL(2,bd);break;//longitudinal acoustic g-Phonon
				case 5:ee=FEKLOAB(ee,bd);TALINDL(2,bd);break;//longitudinal optical g-Phonon
				case 6:ee=FEKLOEM(ee,bd);TALINDL(2,bd);break;//longitudinal optical g-Phonon
				case 7:exit(0);break;//Impact Ionization
				default:break;
			}
		}
	}
	  
//    calculate the final state for energy ee and bandindex iband on an aeuquienergysurface
	  FINALSTATE(bd);
      if(selfscfl)
      {
		 ee=eeold;
         xk=xkold;
         yk=ykold;
         zk=zkold;
         iband=ibold;
         itet=itetold;
         isym=isymold;
      }
      GETGVE(bd);

	  //end of EPSCAT
      return;
}

void Partical::EFPSCAT(int iscat,Band &bd)
{
      //int iscat

	  //Purpose: Calculation of the electron state just after scattering
	  //	  Fischetti's phonon system
	  //Parameter :  iscat = number of actual scattering-prozess

	  //local variables
      bool aefl;
      int ibold,itetold,isymold;
      double eeold,xkold,ykold,zkold;
      double qq,qact,eep,xkctl,ykctl,zkctl;
      double koef,koefmax,aovtet;
      double xqct,yqct,zqct,eepmax;

	  //save particle state
      eeold=ee;
      xkold=xk;
      ykold=yk;
      zkold=zk;
      ibold=iband;
      itetold=itet;
      isymold=isym;

      selfscfl=false;

      iband=iscat/bd.scpre+1+bd.bandof[ipt];
      //GOTO (1,2,3,4,5,6,7),MOD(iscat-1,bd.scpre)+1
	  int cntemp=fmod(iscat,bd.scpre);
	  switch(cntemp){
	  case 0:ee=ee+bd.efopee;
			 isym=fmod(int(48.0*RANDNR()),48);
			 FINALSTATE(bd);
			 if(!selfscfl)
			 {		  
				  xqct=xk-xkold;
				  yqct=yk-ykold;
				  zqct=zk-zkold;
		 //confine the phonon vector to the first BZ
				  bd.CONFBZ(xqct,yqct,zqct);
				  qq=sqrt(xqct*xqct+yqct*yqct+zqct*zqct);
				  if(RANDNR()>EOLINT(qq))
				  {
					  selfscfl=true;
				  }
			 }
		     break;
	  case 1:ee=ee-bd.efopee;
		     isym=fmod(int(48.0*RANDNR()),48);//+1;mini=0
			 FINALSTATE(bd);
		     if(!selfscfl)
			 {		  
				xqct=xk-xkold;
				yqct=yk-ykold;
				zqct=zk-zkold;
		 //confine the phonon vector to the first BZ
				bd.CONFBZ(xqct,yqct,zqct);
				qq=sqrt(xqct*xqct+yqct*yqct+zqct*zqct);
				if(RANDNR()>EOLINT(qq))
				{
				  selfscfl=true;
				}
			 }
			 break;
	  case 2:case 3:
	  case 4:case 5:
					switch(cntemp){
					case 2:FPFSAB(bd);
						   //select maximum phonon energy
						   eepmax=bd.efapeet;
					       //select absorption
						   aefl=true;
						   break;
					case 3:FPFSEM(bd);
						   eepmax=bd.efapeet;
						   aefl=false;
						   break;
					case 4:FPFSAB(bd);
						   eepmax=bd.efapeel;
						   aefl=true;
						   break;
					case 5:FPFSEM(bd);
						   eepmax=bd.efapeel;
						   aefl=false;
						   break;
					default:break;
					}
					if(!selfscfl)
					{
					  //select wedge
				      isym=fmod(int(48.0*RANDNR()),48);//+1;minim=0;
					  //transform center of chosen tetrahedron into wedge
				      bd.OUTWEDGE(xkctl,ykctl,zkctl,bd.xkct[itet],bd.ykct[itet],bd.zkct[itet],bd.matsym,isym);
					  //calculata momentum transfer (approximation)
				      xqct=xkctl-xkold;
					  yqct=ykctl-ykold;
				      zqct=zkctl-zkold;
					  //confine the phonon vector to the first BZ
				      bd.CONFBZ(xqct,yqct,zqct);
					  //calculate phonon wave vector
				      qact=sqrt(xqct*xqct+yqct*yqct+zqct*zqct)*sia0;
					  //calculate phonon energy
				      if(qact>TWOPI)
					  {	  
				           eep=eepmax;
					  }	
				      else
					  {
						   if(qact<0.1)
						   {
							  eep=eepmax*0.176776695*qact;
						   }
						   else
						   {	
							  eep=eepmax*sqrt(1.0-cos(0.25*qact));
						   }
				      }
					  //calculate final particle energy (absorption/emisson)
				      if(aefl)
					  {
				         ee=ee+eep;
					  }
				      else
					  {	
				         ee=ee-eep;
					  }	
					  //tetrahedra which do not contain the energy are rejected
				      if((bd.eek[bd.tet[0][itet]]<ee)||(bd.eek[bd.tet[3][itet]]>ee))
					  {	  
						 //calculate the DOS within the selected tetrahedron (area 
						 //over velocity) for the final energy
				         aovtet=bd.FASTSURF(ee,itet)/bd.vgt[itet];
						 //check the upper bound of the DOS
				         if(aovtet>bd.maxaovtet[itet])
						 {	/* 
				            WRITE (LUTTYO,*) 'EFPSCAT: Max of aovtet too high',
                             aovtet/bd.maxaovtet[itet]
				            WRITE (LUOUT ,*) 'EFPSCAT: Max of aovtet too high',
                             aovtet/bd.maxaovtet[itet]
							 */
						 }
						 //selfscattering due to the upper bound of DOS
				         if(aovtet>RANDNR()*bd.maxaovtet[itet])
						 {	 
							//calculate the coupling factor and the upper bound
				            if(aefl)
							{
				               koefmax=35.0/(eepmax*eepmax);  
				               koef=qact*qact/eep/(exp(eep)-1.0);
							}
				            else
							{
				               koefmax=6.0*PI*PI/eepmax*(1.0/(exp(eepmax)-1.0)+1.0);
				               koef=qact*qact/eep*(1.0/(exp(eep)-1.0)+1.0);
							}
							//check the upper bound
				            if(koef>koefmax)
							{	
								/*
				               WRITE (LUTTYO,*) 'EFPSCAT: Max of Scat. rate too high',
				                               koef/koefmax
				               WRITE (LUOUT ,*) 'EFPSCAT: Max of Scat. rate too high',
				                               koef/koefmax
											   */
							}
							//selfscattering to account for the upper bound of coupling factor
				            if(RANDNR()*koefmax>koef)
							{	
				               selfscfl = true;
							}
				            else
							{
							//calculate final k-vector in the tetrahedron
								FINALK(bd);
								//selscattering to account for the overlap integral
								if(RANDNR()>EOLINT(qact/sia0))
								{
									selfscfl = true;
								}
							}
						 }
				         else
						 {
				            selfscfl = true;
						 }
					  }
				      else
					  {	
				         selfscfl = true;
					  }
					}
					break;
	  case 6:/*
		  WRITE (LUTTYO,*) 'EFPSCAT: No II in EFPSCAT'
		  WRITE (LUOUT ,*) ''
		  STOP
		  */
		  cout<<"error EFPSCAT: No II in EFPSCAT";exit(0);
		  break;
	  default:break;
	  }
	  //If selscattering restore inital particle state
	  if(selfscfl)
	  {	 
			ee=eeold;
			xk=xkold;
			yk=ykold;
			zk=zkold;
			iband=ibold;
			itet=itetold;
			isym=isymold;
	  }
	  GETGVE(bd);
	  //end of EFPSCAT
      return;
}

void Partical::GENBTBTDENS(DevSimulator &dev,Band &bd)
{
	int i,j,k;
	int m;
	int ei,hi,ej,hj;
	int eiinumber,hiinumber,ejjnumber,hjjnumber;
	int eii[10],hii[10],ejj[10],hjj[10];
	bool cfindfl,vfindfl;
	int changefl;//0_变小，1_变大

	//沿着沟道方向BTBT
	for(i=0;i<dev.ngpx-1;i++)
	{
		for(j=0;j<dev.ngpy-1;j++)
		{
//无界面态情况
			for(k=0;k<10;k++)
			{
				ejj[k]=0;
				hjj[k]=0;
			}
			ej=0;hj=0;
			ejjnumber=0;hjjnumber=0;
			cfindfl=false;
			vfindfl=false;
			pc=fabs(RBBL0DENS[i+j*dev.ngpx]);
			if(fabs(pc)>1e-30)
			{
				for(m=0;m<dev.ngpy-2;m++)
				{
					if(dev.mat[i+m*dev.ngpx]==SILICON&&dev.mat[i+(m+1)*dev.ngpx]==SILICON
					 &&(Ei[i+j*dev.ngpx]>=Ec[i+m*dev.ngpx]&&Ei[i+j*dev.ngpx]<Ec[i+(m+1)*dev.ngpx])
					 ||(Ei[i+j*dev.ngpx]<Ec[i+m*dev.ngpx]&&Ei[i+j*dev.ngpx]>=Ec[i+(m+1)*dev.ngpx]))
					{
						ejjnumber++;
						ejj[ejjnumber]=m;
						if(ejjnumber==1)
						{
							if(fabs(dev.gridy[ejj[ejjnumber]]-dev.gridy[j])<fabs(100000000-dev.gridy[j]))
							{
								cfindfl=true;
								ej=m;
							}
						}
						else if(ejjnumber>1)
						{
							if(fabs(dev.gridy[ejj[ejjnumber]]-dev.gridy[j])<fabs(dev.gridy[ejj[ejjnumber-1]]-dev.gridy[j]))
							{
								cfindfl=true;
								ej=m;
							}
						}
					}
				}
				for(m=0;m<dev.ngpy-2;m++)
				{
					if(dev.mat[i+m*dev.ngpx]==SILICON&&dev.mat[i+(m+1)*dev.ngpx]==SILICON
					 &&(Ei[i+j*dev.ngpx]>=Ev[i+m*dev.ngpx]&&Ei[i+j*dev.ngpx]<Ev[i+(m+1)*dev.ngpx])
					 ||(Ei[i+j*dev.ngpx]<Ev[i+m*dev.ngpx]&&Ei[i+j*dev.ngpx]>=Ev[i+(m+1)*dev.ngpx]))
					{
						hjjnumber++;
						hjj[hjjnumber]=m;
						if(hjjnumber==1)
						{
							if(fabs(dev.gridy[hjj[hjjnumber]]-dev.gridy[j])<fabs(100000000-dev.gridy[j]))
							{
								cfindfl=true;
								hj=m;
							}
						}
						else if(hjjnumber>1)
						{
							if(fabs(dev.gridy[hjj[hjjnumber]]-dev.gridy[j])<fabs(dev.gridy[hjj[hjjnumber-1]]-dev.gridy[j]))
							{
								vfindfl=true;
								hj=m;
							}
						}
					}
				}
				if(dev.mat[i+j*dev.ngpx]==SILICON)//cfindfl&&vfindfl)
				{
					pc=-1.0*fabs(pc);
					ipt=PELEC;
					ctype=1;
					xr=(dev.gridx[i]+dev.gridx[i+1])/2.0;
					yr=(dev.gridy[j]+dev.gridy[j+1])/2.0;
					ijqp=i+j*dev.ngpx;
					ee=0.1;
					isym=fmod(int(RANDNR()*48.0),48);
					iband=bd.bandof[ipt];//iband=bd.ibt[itet];
					FINALSTATE(bd);
					dev.npar0+=1;
					dev.npar[ipt]+=1;
					GSEST(2,dev.npar0-1,false,dev,bd);

					pc=fabs(pc);
					ipt=PHOLE;
					ctype=1;
					xr=(dev.gridx[i]+dev.gridx[i+1])/2.0;
					yr=(dev.gridy[j]+dev.gridy[j+1])/2.0;
					ijqp=i+j*dev.ngpx;
					ee=0.1;
					isym=fmod(int(RANDNR()*48.0),48);
					iband=bd.bandof[ipt]+1;//iband=bd.ibt[itet];
					FINALSTATE(bd);
					dev.npar0+=1;
					dev.npar[ipt]+=1;
					GSEST(2,dev.npar0-1,false,dev,bd);
				}
			}
//有界面态情况
			
			for(k=0;k<10;k++)
			{
				ejj[k]=0;
				hjj[k]=0;
			}
			ej=0;hj=0;
			ejjnumber=0;hjjnumber=0;
			cfindfl=false;
			vfindfl=false;
			pc=fabs(RBBL1DENS[i+j*dev.ngpx]);
			if(fabs(pc)>1e-30)
			{
				for(m=0;m<dev.ngpy-2;m++)
				{
					if(dev.mat[i+m*dev.ngpx]==SILICON&&dev.mat[i+(m+1)*dev.ngpx]==SILICON
					 &&(Ei[i+j*dev.ngpx]>=Ec[i+m*dev.ngpx]&&Ei[i+j*dev.ngpx]<Ec[i+(m+1)*dev.ngpx])
					 ||(Ei[i+j*dev.ngpx]<Ec[i+m*dev.ngpx]&&Ei[i+j*dev.ngpx]>=Ec[i+(m+1)*dev.ngpx]))
					{
						ejjnumber++;
						ejj[ejjnumber]=m;
						if(ejjnumber==1)
						{
							if(fabs(dev.gridy[ejj[ejjnumber]]-dev.gridy[j])<fabs(100000000-dev.gridy[j]))
							{
								cfindfl=true;
								ej=m;
							}
						}
						else if(ejjnumber>1)
						{
							if(fabs(dev.gridy[ejj[ejjnumber]]-dev.gridy[j])<fabs(dev.gridy[ejj[ejjnumber-1]]-dev.gridy[j]))
							{
								cfindfl=true;
								ej=m;
							}
						}
					}
				}
				for(m=0;m<dev.ngpy-2;m++)
				{
					if(dev.mat[i+m*dev.ngpx]==SILICON&&dev.mat[i+(m+1)*dev.ngpx]==SILICON
					 &&(Ei[i+j*dev.ngpx]>=Ev[i+m*dev.ngpx]&&Ei[i+j*dev.ngpx]<Ev[i+(m+1)*dev.ngpx])
					 ||(Ei[i+j*dev.ngpx]<Ev[i+m*dev.ngpx]&&Ei[i+j*dev.ngpx]>=Ev[i+(m+1)*dev.ngpx]))
					{
						hjjnumber++;
						hjj[hjjnumber]=m;
						if(hjjnumber==1)
						{
							if(fabs(dev.gridy[hjj[hjjnumber]]-dev.gridy[j])<fabs(100000000-dev.gridy[j]))
							{
								cfindfl=true;
								hj=m;
							}
						}
						else if(hjjnumber>1)
						{
							if(fabs(dev.gridy[hjj[hjjnumber]]-dev.gridy[j])<fabs(dev.gridy[hjj[hjjnumber-1]]-dev.gridy[j]))
							{
								vfindfl=true;
								hj=m;
							}
						}
					}
				}
				if(dev.mat[i+j*dev.ngpx]==SILICON)//cfindfl)
				{
					pc=-1.0*fabs(pc);
					ipt=PELEC;
					ctype=1;
					xr=(dev.gridx[i]+dev.gridx[i+1])/2.0;
					yr=(dev.gridy[j]+dev.gridy[j+1])/2.0;
					ijqp=i+j*dev.ngpx;
					ee=0.1;
					isym=fmod(int(RANDNR()*48.0),48);
					iband=bd.bandof[ipt];//iband=bd.ibt[itet];
					FINALSTATE(bd);
					dev.npar0+=1;
					dev.npar[ipt]+=1;
					GSEST(2,dev.npar0-1,false,dev,bd);
				}
				if(dev.mat[i+j*dev.ngpx]==SILICON)//vfindfl)
				{
					pc=fabs(pc);
					ipt=PHOLE;
					ctype=1;
					xr=(dev.gridx[i]+dev.gridx[i+1])/2.0;
					yr=(dev.gridy[j]+dev.gridy[j+1])/2.0;
					ijqp=i+j*dev.ngpx;
					ee=0.1;
					isym=fmod(int(RANDNR()*48.0),48);
					iband=bd.bandof[ipt]+1;//iband=bd.ibt[itet];
					FINALSTATE(bd);
					dev.npar0+=1;
					dev.npar[ipt]+=1;
					GSEST(2,dev.npar0-1,false,dev,bd);
				}
			}
		}
	}

	//垂直沟道方向BTBT
	for(j=0;j<dev.ngpy-2;j++)
	{
		for(i=0;i<dev.ngpx-2;i++)
		{			
			for(k=0;k<10;k++)
			{
				eii[k]=0;
				hii[k]=0;
			}
			ei=0;hi=0;
			eiinumber=0;hiinumber=0;
			cfindfl=false;
			vfindfl=false;

			if(dev.mat[i+j*dev.ngpx]==SILICON)
			{
			ipt=PELEC;
			pc=-1.0*fabs(RBBT1DENS[i+j*dev.ngpx]);
			if(fabs(pc)>1e-30)
			{
				for(m=0;m<dev.ngpx-2;m++)
				{
					if(dev.mat[m+j*dev.ngpx]==SILICON&&dev.mat[m+1+j*dev.ngpx]==SILICON
					 &&(Ei[i+j*dev.ngpx]>=Ec[m+j*dev.ngpx]&&Ei[i+j*dev.ngpx]<Ec[m+1+j*dev.ngpx])
					 ||(Ei[i+j*dev.ngpx]<Ec[m+j*dev.ngpx]&&Ei[i+j*dev.ngpx]>=Ec[m+1+j*dev.ngpx]))
					{
						eiinumber++;
						eii[eiinumber]=m;
						if(eiinumber==1)
						{
							if(fabs(dev.gridx[eii[eiinumber]]-dev.gridx[i])<fabs(100000000-dev.gridx[i]))
							{
								cfindfl=true;
								ei=m;
							}
						}
						else if(eiinumber>1)
						{
							if(fabs(dev.gridx[eii[eiinumber]]-dev.gridx[i])<fabs(dev.gridx[eii[eiinumber-1]]-dev.gridx[i]))
							{
								cfindfl=true;
								ei=m;
							}
						}
					}
				}
				for(m=0;m<dev.ngpx-2;m++)
				{
					if(dev.mat[m+j*dev.ngpx]==SILICON&&dev.mat[m+1+j*dev.ngpx]==SILICON
					 &&(Ei[i+j*dev.ngpx]>=Ev[m+j*dev.ngpx]&&Ei[i+j*dev.ngpx]<Ev[m+1+j*dev.ngpx])
					 ||(Ei[i+j*dev.ngpx]<Ev[m+j*dev.ngpx]&&Ei[i+j*dev.ngpx]>=Ev[m+1+j*dev.ngpx]))
					{
						hiinumber++;
						hii[hiinumber]=m;
						if(hiinumber==1)
						{
							if(fabs(dev.gridx[hii[hiinumber]]-dev.gridx[i])<fabs(100000000-dev.gridx[i]))
							{
								vfindfl=true;
								hi=m;
							}
						}
						else if(hiinumber>1)
						{
							if(fabs(dev.gridx[hii[hiinumber]]-dev.gridx[i])<fabs(dev.gridx[hii[hiinumber-1]]-dev.gridx[i]))
							{
								vfindfl=true;
								hi=m;
							}
						}
					}
				}
///*				
				if(dev.mat[i+j*dev.ngpx]==SILICON)//cfindfl)
				{
					ctype=2;
					xr=(dev.gridx[i]+dev.gridx[i+1])/2.0;
					yr=(dev.gridy[j]+dev.gridy[j+1])/2.0;
					ijqp=i+j*dev.ngpx;
					ee=0.1;
					isym=fmod(int(RANDNR()*48.0),48);
					iband=bd.bandof[ipt];//iband=bd.ibt[itet];
					FINALSTATE(bd);
					dev.npar0+=1;
					dev.npar[ipt]+=1;
					GSEST(2,dev.npar0-1,false,dev,bd);
				}
//*/
			}
/*			
			ipt=PHOLE;
			pc=fabs(RBBT1DENS[i+j*dev.ngpx]);
			if(fabs(pc)>1e-30)
			{
				if(dev.mat[i+j*dev.ngpx]==SILICON)//vfindfl)
				{
					ctype=2;
					xr=(dev.gridx[i]+dev.gridx[i+1])/2.0;
					yr=(dev.gridy[j]+dev.gridy[j+1])/2.0;
					ijqp=i+j*dev.ngpx;
					ee=0.1;
					isym=fmod(int(RANDNR()*48.0),48);
					iband=bd.bandof[ipt]+1;//iband=bd.ibt[itet];
					FINALSTATE(bd);		
					dev.npar0+=1;
					dev.npar[ipt]+=1;
					GSEST(2,dev.npar0-1,false,dev,bd);
				}
			}
*/
		}
		}
	}
	
	return;
}
void Partical::ESCATII(int iscat,DevSimulator &dev,Band &bd)
{
      //int iscat

	  //Purpose: Calculation of the electron state just after scattering
	  //=======  for II

	  //local variables
      int ibold,itetold,isymold;
      double eeold,xkold,ykold,zkold;

      nsctype[iscat][iband-bd.bandof[ipt]]+=1;

	  //save particle state
      eeold=ee;
      xkold=xk;
      ykold=yk;
      zkold=zk;
      ibold=iband;
      itetold=itet;
      isymold=isym;

      iband=iscat/bd.scpre+bd.bandof[ipt];

//    Impact Ionization
//    calculate the final state for energy ee and band index iband on an a euquienergysurface
      ee=(ee-sieg)/3.0;
      isym=fmod(int(RANDNR()*48.0),48);
      FINALSTATE(bd);

      if(selfscfl)
      {   
		 ee=eeold;
         xk=xkold;
         yk=ykold;
         zk=zkold;
         iband=ibold;
         itet=itetold;
         isym=isymold;
	  }	
      else
	  {	 
		 if(bd.seciifl)
		 {
		 if((dev.npar0+2)>MNPAR)
		 {
		 }	
         //generate secondary electron
         dev.npar0=dev.npar0+1;
         dev.npar[ipt]+=1;
		 GSEST(2,dev.npar0-1,false,dev,bd);

		 //generate secondary hole
         ipt=PHOLE;
         iband=bd.bandof[ipt]-1;
         selfscfl=true;
         while(selfscfl)
         {
			iband=iband+1;
            if((iband-bd.bandof[ipt])>=bd.nband[ipt])
            {   
				/*
				WRITE (LUTTYO,*)  'ESCATII: No final state for secondary hole'     
                WRITE (LUOUT ,*)  'ESCATII: No final state for secondary hole'     
                STOP
				*/
            }
            isym=fmod(int(RANDNR()*48.0),48);
            FINALSTATE(bd);
         }
         pc=-pc;
         dev.npar0=dev.npar0+1;
         dev.npar[ipt]+=1;
         GSEST(2,dev.npar0-1,false,dev,bd);

		 //get state of primary electron
         GSEST(1,dev.npar0-2,false,dev,bd);
		 //primary electron has negative k-vector of secondary electron
         xk=-xk;     
         yk=-yk;    
         zk=-zk;
		 isym=bd.indmat[-bd.matsym[3][isym]+1][-bd.matsym[4][isym]+1][-bd.matsym[5][isym]+1]      
             [bd.matsym[0][isym]-1][bd.matsym[1][isym]-1][bd.matsym[2][isym]-1];
		 /*
         isym=bd.indmat(-bd.matsym(4,isym),-bd.matsym(5,isym),-bd.matsym(6,isym),      
		    bd.matsym(1,isym),bd.matsym(2,isym),bd.matsym(3,isym))*/
		 }
	  }
      GETGVE(bd);
      return;
      //END
}

void Partical::FINALK(Band &bd)
{

//Purpose: - calculate the final k-vector in the chosen tetrahedron an equienergy surface

      double rv1,rv2,rfac,xkl,ykl,zkl;
      double ax,ay,az,bx,by,bz,cx,cy,cz;
      double dx,dy,dz,a1,a2,a3,a4,ar;

	  //get pair of random numbers, which is homogeneously distributed in 
	  //a triangle
      rv1=1.0-sqrt(RANDNR());
      rv2=RANDNR()*(1.0-rv1);

	  //triangular intersection with energy lower e2
      if(ee<=bd.eek[bd.tet[1][itet]])
	  {		  
		 //cal. vectors (relative to origin of tetrahedron) pointing to the 
		 //three nodes of the triangle
         rfac=(ee-bd.eek[bd.tet[0][itet]])/(bd.eek[bd.tet[1][itet]]-bd.eek[bd.tet[0][itet]]);
         ax=rfac*(bd.xkk[bd.tet[1][itet]]-bd.xkk[bd.tet[0][itet]]);
         ay=rfac*(bd.ykk[bd.tet[1][itet]]-bd.ykk[bd.tet[0][itet]]);
         az=rfac*(bd.zkk[bd.tet[1][itet]]-bd.zkk[bd.tet[0][itet]]);
         if(rfac<0.0||rfac>1.0)
		 {
		 }
         rfac=(ee-bd.eek[bd.tet[0][itet]])/(bd.eek[bd.tet[2][itet]]-bd.eek[bd.tet[0][itet]]);
         bx=rfac*(bd.xkk[bd.tet[2][itet]]-bd.xkk[bd.tet[0][itet]]);
         by=rfac*(bd.ykk[bd.tet[2][itet]]-bd.ykk[bd.tet[0][itet]]);
         bz=rfac*(bd.zkk[bd.tet[2][itet]]-bd.zkk[bd.tet[0][itet]]);
         if(rfac<0.0||rfac>1.0)
		 {
		 }
   
         rfac=(ee-bd.eek[bd.tet[0][itet]])/(bd.eek[bd.tet[3][itet]]-bd.eek[bd.tet[0][itet]]);
         cx=rfac*(bd.xkk[bd.tet[3][itet]]-bd.xkk[bd.tet[0][itet]]);
         cy=rfac*(bd.ykk[bd.tet[3][itet]]-bd.ykk[bd.tet[0][itet]]);
         cz=rfac*(bd.zkk[bd.tet[3][itet]]-bd.zkk[bd.tet[0][itet]]);
         if(rfac<0.0||rfac>1.0)
		 {
		 }	
         xkl=rv1*(bx-ax)+rv2*(cx-ax)+ax+bd.xkk[bd.tet[0][itet]];
         ykl=rv1*(by-ay)+rv2*(cy-ay)+ay+bd.ykk[bd.tet[0][itet]];
         zkl=rv1*(bz-az)+rv2*(cz-az)+az+bd.zkk[bd.tet[0][itet]];
	}
	//quadrangular intersection with energy between e2 and e3
	else if(ee<bd.eek[bd.tet[2][itet]])
	{
			//cal. vectors (relative to origin of k-space) pointing to the 
			//four nodes of the quadrangle
			rfac=(ee-bd.eek[bd.tet[0][itet]])/(bd.eek[bd.tet[2][itet]]-bd.eek[bd.tet[0][itet]]);
			ax=rfac*(bd.xkk[bd.tet[2][itet]]-bd.xkk[bd.tet[0][itet]])+bd.xkk[bd.tet[0][itet]];
			ay=rfac*(bd.ykk[bd.tet[2][itet]]-bd.ykk[bd.tet[0][itet]])+bd.ykk[bd.tet[0][itet]];
			az=rfac*(bd.zkk[bd.tet[2][itet]]-bd.zkk[bd.tet[0][itet]])+bd.zkk[bd.tet[0][itet]];
			if(rfac<0.0||rfac>1.0)
			{
			}
			rfac=(ee-bd.eek[bd.tet[1][itet]])/(bd.eek[bd.tet[2][itet]]-bd.eek[bd.tet[1][itet]]);
			bx=rfac*(bd.xkk[bd.tet[2][itet]]-bd.xkk[bd.tet[1][itet]])+bd.xkk[bd.tet[1][itet]];
			by=rfac*(bd.ykk[bd.tet[2][itet]]-bd.ykk[bd.tet[1][itet]])+bd.ykk[bd.tet[1][itet]];
			bz=rfac*(bd.zkk[bd.tet[2][itet]]-bd.zkk[bd.tet[1][itet]])+bd.zkk[bd.tet[1][itet]];
			if(rfac<0.0||rfac>1.0)
            {
			}
			rfac=(ee-bd.eek[bd.tet[0][itet]])/(bd.eek[bd.tet[3][itet]]-bd.eek[bd.tet[0][itet]]);
			cx=rfac*(bd.xkk[bd.tet[3][itet]]-bd.xkk[bd.tet[0][itet]])+bd.xkk[bd.tet[0][itet]];
			cy=rfac*(bd.ykk[bd.tet[3][itet]]-bd.ykk[bd.tet[0][itet]])+bd.ykk[bd.tet[0][itet]];
			cz=rfac*(bd.zkk[bd.tet[3][itet]]-bd.zkk[bd.tet[0][itet]])+bd.zkk[bd.tet[0][itet]];
			if(rfac<0.0||rfac>1.0)
			{
			}
			rfac=(ee-bd.eek[bd.tet[1][itet]])/(bd.eek[bd.tet[3][itet]]-bd.eek[bd.tet[1][itet]]);
			dx=rfac*(bd.xkk[bd.tet[3][itet]]-bd.xkk[bd.tet[1][itet]])+bd.xkk[bd.tet[1][itet]];
			dy=rfac*(bd.ykk[bd.tet[3][itet]]-bd.ykk[bd.tet[1][itet]])+bd.ykk[bd.tet[1][itet]];
			dz=rfac*(bd.zkk[bd.tet[3][itet]]-bd.zkk[bd.tet[1][itet]])+bd.zkk[bd.tet[1][itet]];
			if(rfac<0.0||rfac>1.0)
			{
			}
			//cal. area of the four possible triangles made from the four nodes
			//of the quadrangle
			a1=bd.FII(ax,ay,az,bx,by,bz,cx,cy,cz);
			a2=bd.FII(dx,dy,dz,bx,by,bz,cx,cy,cz);
			a3=bd.FII(ax,ay,az,dx,dy,dz,cx,cy,cz);
			a4=bd.FII(ax,ay,az,bx,by,bz,dx,dy,dz);
			ar=RANDNR()*(a1+a2+a3+a4);     

			//select one of the four triangles and calculate k-vector
			if(ar<a1)
			{
				 xkl=rv1*(bx-ax)+rv2*(cx-ax)+ax;
				 ykl=rv1*(by-ay)+rv2*(cy-ay)+ay;
	             zkl=rv1*(bz-az)+rv2*(cz-az)+az;
			 }
			else if(ar<(a1+a2))
			{
				xkl=rv1*(bx-dx)+rv2*(cx-dx)+dx;
				ykl=rv1*(by-dy)+rv2*(cy-dy)+dy;
				zkl=rv1*(bz-dz)+rv2*(cz-dz)+dz;
			}
			else if(ar<(a1+a2+a3))
			{
				xkl=rv1*(dx-ax)+rv2*(cx-ax)+ax;
				ykl=rv1*(dy-ay)+rv2*(cy-ay)+ay;
				zkl=rv1*(dz-az)+rv2*(cz-az)+az;
			} 
			else
			{
				xkl=rv1*(bx-ax)+rv2*(dx-ax)+ax;
				ykl=rv1*(by-ay)+rv2*(dy-ay)+ay;
				zkl=rv1*(bz-az)+rv2*(dz-az)+az;
			}
			//triangular intersection with energy higher e3
		}
		else
		{
			//cal. vectors (relative to origin of tetrahedron) pointing to the 
			//three nodes of the triangle
			rfac=(ee-bd.eek[bd.tet[3][itet]])/(bd.eek[bd.tet[0][itet]]-bd.eek[bd.tet[3][itet]]);
			ax=rfac*(bd.xkk[bd.tet[0][itet]]-bd.xkk[bd.tet[3][itet]]);
			ay=rfac*(bd.ykk[bd.tet[0][itet]]-bd.ykk[bd.tet[3][itet]]);
			az=rfac*(bd.zkk[bd.tet[0][itet]]-bd.zkk[bd.tet[3][itet]]);
			if(rfac<0.0&&rfac>1.0)
			{
			}
			rfac=(ee-bd.eek[bd.tet[3][itet]])/(bd.eek[bd.tet[1][itet]]-bd.eek[bd.tet[3][itet]]);
			bx=rfac*(bd.xkk[bd.tet[1][itet]]-bd.xkk[bd.tet[3][itet]]);
			by=rfac*(bd.ykk[bd.tet[1][itet]]-bd.ykk[bd.tet[3][itet]]);
			bz=rfac*(bd.zkk[bd.tet[1][itet]]-bd.zkk[bd.tet[3][itet]]);
			if(rfac<0.0&&rfac>1.0)
			{
			}
	        rfac=(ee-bd.eek[bd.tet[3][itet]])/(bd.eek[bd.tet[2][itet]]-bd.eek[bd.tet[3][itet]]);
			cx=rfac*(bd.xkk[bd.tet[2][itet]]-bd.xkk[bd.tet[3][itet]]);
			cy=rfac*(bd.ykk[bd.tet[2][itet]]-bd.ykk[bd.tet[3][itet]]);
			cz=rfac*(bd.zkk[bd.tet[2][itet]]-bd.zkk[bd.tet[3][itet]]);
			if(rfac<0.0&&rfac>1.0)
			{
			}	
			xkl=rv1*(bx-ax)+rv2*(cx-ax)+ax+bd.xkk[bd.tet[3][itet]];
			ykl=rv1*(by-ay)+rv2*(cy-ay)+ay+bd.ykk[bd.tet[3][itet]];
			zkl=rv1*(bz-az)+rv2*(cz-az)+az+bd.zkk[bd.tet[3][itet]];
		}
	  //end of finding the final state
	  //apply the symmetry transformation
      bd.OUTWEDGE (xk,yk,zk,xkl,ykl,zkl,bd.matsym,isym);
	  //End of FINALK
      return;
}

void Partical::FPFSAB(Band &bd)
{
//Purpose: - calculate the final state for ee and iband on
//           an equienergy surface for Fischetti 
//           acoustic phonon scattering (absorbtion)

	  //local variables
      int itab,nted,count2;
      double maxdostet;

	  //energy index of tetraheder list
      itab=int(MAX(0.0,ee)*bd.dlistfp);//+1//minim=0;
      if(itab>=MWLE)
	  {
		  itab=MWLEFP-1;//maxim=MWLEFP-1
	  }	

	  //get maximum dos of tetraheda
	  //for the given energy
      maxdostet=bd.maxaovfpab[iband][itab];

	  //Number of tetraheder in list
      nted=bd.ntlistfpab[itab][iband];
      count2=0;
      if(nted==0)
	  {
	  //no final state possible for this band and energy
         selfscfl=true;
	  }
	  else
	  {
		  //chose random tetraheder from list
		  itet=bd.tlistfpab[int(nted*RANDNR())+bd.ptlistfpab[itab][iband]];
		  while(((RANDNR()*maxdostet)>bd.maxaovtet[itet])&&(count2<10000))
		  {
			itet=bd.tlistfpab[int(nted*RANDNR())+bd.ptlistfpab[itab][iband]];
		    count2=count2+1;
		  //acception/rejection due to dos in tetrahedron
		  }

		  //check band index
	      if(iband!=bd.ibt[itet])
		  {
			  cout<<"error FPFSAB: iband <> bd.ibt";exit(0);
 		  }
		  //no final state found
		  if(count2>=10000)
		  {
			  cout<<"error FPFSAB: count2 > 10000 iband = ";exit(0);
	      }
	  }
	  //end of FPFSAB
      return;

}

void Partical::FPFSEM(Band &bd)
{
//Purpose: - calculate the final state for ee and iband on
	  //    an aequienergy surface for Fischetti 
	  //    acoustic phonon scattering (emisson)

	  //local variables
      int itab,nted,count2;
      double maxdostet;

	  //energy index of tetraheder list
      itab=int(MAX(0.0,ee)*bd.dlistfp);//+1;///minim=0
      if(itab>=MWLE) 
	  {
		  itab=MWLEFP-1;//maxim=MWLEFP-1
	  }	
	  //get maximum dos of tetraheda
	  //for the given energy
      maxdostet=bd.maxaovfpem[iband][itab];

	  //Number of tetraheder in list
      nted=bd.ntlistfpem[itab][iband];
      count2=0;
      if(nted==0)
	  {
		 //no final state possible for this band and energy
         selfscfl=true;
      }
	  else
	  {
	  //chose random tetraheder from list
		  itet= bd.tlistfpem[int(nted*RANDNR())+bd.ptlistfpem[itab][iband]];
		  while(((RANDNR()*maxdostet)>bd.maxaovtet[itet])&&(count2<10000))
		  {
			itet= bd.tlistfpem[int(nted*RANDNR())+bd.ptlistfpem[itab][iband]];
			count2=count2+1;
		  }	
	  //acception/rejection due to dos in tetrahedron
 
	  //check band index
		  if(iband!=bd.ibt[itet])
		  {
			  cout<<"error FPFSEM: iband <> bd.ibt";exit(0);
		  }	
	  //no final state found
	      if(count2>=10000)
		  {
			  cout<<"error FPFSEM: count2 > 10000 iband = ";exit(0);
		  }
	  }
	  //end of FPFSEM
      return;
}

void Partical::HPSCAT(int iscat,Band &bd)
{
	  //int iscat;

	  //Purpose: Calculation of the hole state just after scattering
	  //Parameter :  iscat = number of actual scattering-prozess

	  //local variables
      int ibold,itetold,isymold;
      double eeold,xkold,ykold,zkold;

	  //save inital particle state
      eeold=ee;
      xkold=xk;
      ykold=yk;
      zkold=zk;
      ibold=iband;
      itetold=itet;
      isymold=isym;

      //calculate final band and type of scattering
      iband=iscat/bd.scprh+bd.bandof[ipt];
      //GOTO (1,2,3,4),MOD(iscat-1,bd.scprh)+1
	  int ntemp=fmod(iscat,bd.scprh);
	  switch(ntemp)
	  {
		case 0:ee=ee;break;//elastic intravalley
		case 1:ee=ee+bd.temphop;break;//absorption
		case 2:ee=ee-bd.temphop;break;//emission
		default:break;//Impact Ionization
	  }
	      
//    calculate the final state for energy ee and bandindex iband on an aeuquienergy surface

	  isym=fmod(int(RANDNR()*48.0),48);
      FINALSTATE(bd);

      if(selfscfl)
	  {
          ee=eeold;
          xk=xkold;
          yk=ykold;
          zk=zkold;
          iband=ibold;
          itet=itetold;
          isym=isymold;
      }
      GETGVE(bd);
	  return;
}

void Partical::HSCATII(int iscat,DevSimulator &dev,Band &bd)
{
	  //int iscat

	  //Purpose: Calculation of the hole state just after scattering
	  // for Impact Ionization scatter
	  //Parameter :  iscat = number of actual scattering-prozess

	  //local variables
      int ibold,itetold,isymold;
      double eeold,xkold,ykold,zkold;

	  //save inital particle state
      eeold=ee;
      xkold=xk;
      ykold=yk;
      zkold=zk;
      ibold=iband;
      itetold=itet;
      isymold=isym;

	  //calculate final band and type of scattering
      iband=iscat/bd.scprh+bd.bandof[ipt];

	  //Impact Ionization
	  //calculate the final state for energy
	  //ee and bandindex iband on an aeuquienergy surface
      ee=(ee-sieg)/3.0;
	
      isym=fmod(int(RANDNR()*48.0),48);
      FINALSTATE(bd);

      if(selfscfl)
      {
		 ee=eeold;
         xk=xkold;
         yk=ykold;
         zk=zkold;
         iband=ibold;
         itet=itetold;
         isym=isymold;
	  }	
      else if(bd.seciifl)
	  {
 		 if((dev.npar0+2)>MNPAR)
		 {
		 }
			 //generate secondary hole
		 dev.npar0=dev.npar0+1;
		 dev.npar[ipt]=dev.npar[ipt]+1;
		 GSEST(2,dev.npar0-1,false,dev,bd);

			 //generate secondary electron
		 ipt=PELEC;
		 iband=bd.bandof[ipt]-1;
		 selfscfl=true;
		 while(selfscfl)
		 {
			 iband=iband+1;
			 if((iband-bd.bandof[ipt])>=bd.nband[ipt])
			 {
			 }
			 isym=fmod(int(RANDNR()*48.0),48);
			 FINALSTATE(bd);
		 }
		 pc=-pc;
         dev.npar0=dev.npar0+1;
         dev.npar[ipt]=dev.npar[ipt]+1;
         GSEST(2,dev.npar0-1,false,dev,bd);

		 //get state of primary hole
         GSEST(1,dev.npar0-2,false,dev,bd);
		 //primary hole has negative k-vector of secondary hole
         xk=-xk;     
         yk=-yk;     
         zk=-zk;    
         isym=bd.indmat[-bd.matsym[3][isym]+1][-bd.matsym[4][isym]+1][-bd.matsym[5][isym]+1]      
             [bd.matsym[0][isym]-1][bd.matsym[1][isym]-1][bd.matsym[2][isym]-1];
      
	  }
      GETGVE(bd);

      return;

}

void Partical::OEPSCAT(int iscat,Band &bd)
{
      //int iscat

      //Calculation of the oxide electron state just after scattering
      //Parameter :  iscat = number of actual scattering-prozess

      //local variables
      //INTEGER ibold, itetold, isymold
	  int ibold, itetold, isymold;
      //DOUBLE PRECISION RANDNR, eeold, xkold, ykold, zkold
	  double eeold, xkold, ykold, zkold;
	  //save old particle state
      eeold=ee;
      xkold=xk;
      ykold=yk;
      zkold=zk;
      ibold=iband;
      itetold=itet;
      isymold=isym;

      //calculate new band and type of scattering
      iband=iscat/bd.scproe+ bd.bandof[ipt];
      //GOTO (1,2,3),MOD(iscat-1,bd.scproe)+1
	  //elastic intravalley
	  //1
	  int ntemp=fmod(iscat,bd.scproe);
	  switch(ntemp)
	  {
		case 0:ee=ee;break;//elastic intravalley
		case 1:ee=ee+bd.tempoeop;break;//absorption
		case 2:ee=ee-bd.tempoeop;break;//emission
		default:break;
	  }

//    calculate the final state for energy ee and bandindex iband on an aeuquienergy surface
	  isym=fmod(int(RANDNR()*48.0),48);
      //CALL FINALSTATE
	  FINALSTATE(bd);	
      if(selfscfl)
      {
		 ee=eeold;
         xk=xkold;
         yk=ykold;
         zk=zkold;
         iband=ibold;
         itet=itetold;
         isym=isymold;
	  }	
      GETGVE(bd);
      return;
}

inline double Partical::EOLINT(double qq)
{

	  //Purpose: Calculation of the electron overlap integral
	  //local variables
      double qa,rtnv;
      
	  qa=qq*sia0*0.3907963210;
      if(qa>0.01)
	  {
         //EOLINT=(3.0*(sin(qa)-qa*cos(qa))/pow(qa,3.0))**2;
		 rtnv=pow((3.0*(sin(qa)-qa*cos(qa))/pow(qa,3.0)),2.0); 
      }
      else
	  {
		  rtnv=1.0; 
	  }	

//    End of EOLINT
	  return rtnv;
}

void Partical::STRATTAU(double &tau,double &dtau,
						double eed,int iinital,DevSimulator &dev,Band &bd)
{
 //   double tau,dtau,eed
 //   int iinital
//____Purpose : 
//    Parameter : eed=energy of actual electron
//           : ip=particle type
      int itab;
      double intp;

//____smallest allowed energy
      if(eed<=0.0)
	  {
         tau=bd.sumscatt[0][iinital];
         dtau=(bd.sumscatt[1][iinital]-bd.sumscatt[0][iinital])
              /(bd.energy[1]-bd.energy[0]);
//____largest allowed energy
      }
	  else if(eed>=bd.emax)
	  {
         tau=bd.sumscatt[MTAB][iinital];
         dtau=(bd.sumscatt[MTAB][iinital]-bd.sumscatt[MTAB-1][iinital])
			 /(bd.energy[MTAB]-bd.energy[MTAB-1]);
      }
	  else
	  {
         itab=(int)(eed/bd.dtable);
         intp=(eed-bd.energy[itab])/(bd.energy[itab+1]-bd.energy[itab]);
         tau=bd.sumscatt[itab][iinital]+intp 
          *(bd.sumscatt[itab+1][iinital]
             -bd.sumscatt[itab][iinital]);
         dtau=(bd.sumscatt[itab+1][iinital]-bd.sumscatt[itab][iinital])
              /(bd.energy[itab+1]-bd.energy[itab]);
      }

	  if(bd.bhfl&&bd.bhmrtfl&&iinital==bd.bandof[PELEC]){//THEN be careful!!!!!!!!!
         tau=tau+EBHSCRT(eed,dev,bd);
         dtau=dtau 
            +1e3*(EBHSCRT(eed+1e-3*bd.dtable,dev,bd)-EBHSCRT(eed,dev,bd))/bd.dtable;
      }

      if(tau>0.0)
	  {
         if(dtau!=0.0)dtau=-dtau/(tau*tau);
         tau=1e0/tau;
      }
//____End of STRATTAU
      return;
}

void Partical::EINIT(void)
{
//     Purpose: - initialization of electron related values for each run
//     -------
//____local variables
      int i,j,icont;
      charsign[PELEC]=-1.0;
      charsign[PHOLE]=+1.0;
      charsign[POXEL]=-1.0;


//____set cumulative variables to zero
      for(i=0;i<MNB;i++)
	  {
//_______integer
         ntet[i]=0;
         nbz[i]=0;
         nquad[i]=0;

         ntotp[i]=0;
         nreap[i]=0;
         nslfp[i]=0;

         nreaii[i]=0;

//_______double
         dtotp[i]=0.0;

      }

      for(i=0;i<NBE;i++)
	  {
		  for(j=0;j<MSCPRE*NBE;j++)
		  {
            nsctype[j][i]=0;
          }
      }

       for(i=0;i<NBH;i++)
	  {
          for(j=0;j<MSCPRH*NBH;j++)
		  {
            nsctyph[j][i]=0;
          }
      }

      for(i=0;i<NBOE;i++)
	  {
          for(j=0;j<MSCPROE*NBOE;j++)
		  {
            nsctypoe[j][i]=0;
          }
      }
	  //____Si/SiO2 injection quantities
      for(i=0;i<NPARTYP;i++)
	  {
         noxtotal[i]=0;
         noxreal[i]=0;
         noxrefl[i]=0;
         doxreal[i]=0;
      }

//____tunneling charge
      pctun=0.0;	
	  return;
}

void Partical::GETSTAT(DevSimulator &dev,Band &bd)
{
//____local variables
      int ier,ibase,itab,icont;
      double dxrr,dyrr,vv,pii;
      double fac1,fac2,fac3,fac4;//,HIIRATE
      double imxx,imxy,imxz,imyy,imyz,imzz;
      double msort[3][3];
	  int carriertype;
//    .. square of velocity
      vv=bd.vgt[itet]*bd.vgt[itet];
//    .. avalanche generation
      if (ipt==PELEC)
	  {
         pii=bd.EIIRATE(ee);
      }
	  else if (ipt==PHOLE)
	  {
         pii=bd.HIIRATE(ee);
      }
	  else
	  {
         pii=0.0;
      }
//    ..enesemble average of II
      dev.iicurts[ipt]=dev.iicurts[ipt]+pii*fabs(pc);
//    ..energy index of energy distribution function
      ier=(int)(ee/dev.edfemax[ipt]*(double)(NEREDF));//+1
      if (ier>=NEREDF)ier=NEREDF-1;//maxim=NEREDF-1
//    ..homogeneous energy distribution
      itab=(int)(ee/bd.dtable);//+1;////////////becareful itab
      if (itab<MTAB) 
         dev.edishom[itab][ipt]=dev.edishom[itab][ipt]+fabs(pc);
      
//    ..relative position of particle in quadrant
      dxrr=(xr-dev.gridx[iqp])/(dev.gridx[iqp+1]-dev.gridx[iqp]);
      dyrr=(yr-dev.gridy[jqp])/(dev.gridy[jqp+1]-dev.gridy[jqp]);
      fac1=(1.0-dxrr)*(1.0-dyrr)*fabs(pc);
      fac2=(     dxrr)*(1.0-dyrr)*fabs(pc);
      fac3=(1.0-dxrr)*(     dyrr)*fabs(pc);
      fac4=(     dxrr)*(     dyrr)*fabs(pc);
	  carriertype=ctype;
      for(icont=0;icont<dev.ncont;icont++)
	  {
		  if(direction==100)
		  {
			  dev.curcont[icont]+=pc*(xv*dev.xgradrs[ijqp][icont]
								  +yv*dev.ygradrs[ijqp][icont]);
			  dev.curcontt[carriertype][icont]+=pc*(xv*dev.xgradrs[ijqp][icont]
								                +yv*dev.ygradrs[ijqp][icont]);
		  }
		  else if(direction==110)
		  {
			  dev.curcont[icont]+=pc*(oldxv*dev.xgradrs[ijqp][icont]
								  +oldyv*dev.ygradrs[ijqp][icont]);
			  dev.curcontt[carriertype][icont]+=pc*(oldxv*dev.xgradrs[ijqp][icont]
											    +oldyv*dev.ygradrs[ijqp][icont]);
		  }
      }
//____bilinear mapping of quantities to grid points
//    ..left upper grid point
      //ibase=NSTAT*(ijqp-1+MNGP*(ipt-1)) 
	  
      dev.statis[ARHO][ijqp][ipt]=dev.statis[ARHO][ijqp][ipt] 
                               +fac1;
      dev.statis[AVX ][ijqp][ipt]=dev.statis[AVX ][ijqp][ipt] 
                               +fac1*xv;
      dev.statis[AVY ][ijqp][ipt]=dev.statis[AVY ][ijqp][ipt] 
                               +fac1*yv;
      dev.statis[AVZ ][ijqp][ipt]=dev.statis[AVZ ][ijqp][ipt] 
                               +fac1*zv;
      dev.statis[AVV ][ijqp][ipt]=dev.statis[AVV ][ijqp][ipt] 
                               +fac1*vv;
      dev.statis[AE  ][ijqp][ipt]=dev.statis[AE  ][ijqp][ipt] 
                               +fac1*ee;
      dev.statis[AII ][ijqp][ipt]=dev.statis[AII ][ijqp][ipt] 
                               +fac1*pii;
      dev.statis[AVXX ][ijqp][ipt]=dev.statis[AVXX ][ijqp][ipt] 
                           +fac1*xv*xv;
      dev.statis[AVXY ][ijqp][ipt]=dev.statis[AVXY ][ijqp][ipt] 
                           +fac1*xv*yv;
      dev.statis[AVXZ ][ijqp][ipt]=dev.statis[AVXZ ][ijqp][ipt] 
                           +fac1*xv*zv;
      dev.statis[AVYY ][ijqp][ipt]=dev.statis[AVYY ][ijqp][ipt] 
                           +fac1*yv*yv;
      dev.statis[AVYZ ][ijqp][ipt]=dev.statis[AVYZ ][ijqp][ipt] 
                           +fac1*yv*zv;
      dev.statis[AVZZ ][ijqp][ipt]=dev.statis[AVZZ ][ijqp][ipt] 
                           +fac1*zv*zv;
      dev.statis[AEVX ][ijqp][ipt]=dev.statis[AEVX ][ijqp][ipt] 
                           +fac1*ee*xv;
      dev.statis[AEVY ][ijqp][ipt]=dev.statis[AEVY ][ijqp][ipt] 
                           +fac1*ee*yv;
//    ..left lower grid point
      //ibase=NSTAT*(ijqp+1-1+MNGP*(ipt-1))  
      
      dev.statis[ARHO][ijqp+1][ipt]=dev.statis[ARHO][ijqp+1][ipt] 
                               +fac2;
      dev.statis[AVX ][ijqp+1][ipt]=dev.statis[AVX ][ijqp+1][ipt] 
                               +fac2*xv;
      dev.statis[AVY ][ijqp+1][ipt]=dev.statis[AVY ][ijqp+1][ipt] 
                               +fac2*yv;
      dev.statis[AVZ ][ijqp+1][ipt]=dev.statis[AVZ ][ijqp+1][ipt] 
                               +fac2*zv;
      dev.statis[AVV ][ijqp+1][ipt]=dev.statis[AVV ][ijqp+1][ipt] 
                               +fac2*vv;
      dev.statis[AE  ][ijqp+1][ipt]=dev.statis[AE  ][ijqp+1][ipt] 
                               +fac2*ee;
      dev.statis[AII ][ijqp+1][ipt]=dev.statis[AII ][ijqp+1][ipt] 
                               +fac2*pii;
      dev.statis[AVXX ][ijqp+1][ipt]=dev.statis[AVXX ][ijqp+1][ipt] 
                           +fac2*xv*xv;
      dev.statis[AVXY ][ijqp+1][ipt]=dev.statis[AVXY ][ijqp+1][ipt] 
                           +fac2*xv*yv;
      dev.statis[AVXZ ][ijqp+1][ipt]=dev.statis[AVXZ ][ijqp+1][ipt] 
                           +fac2*xv*zv;
      dev.statis[AVYY ][ijqp+1][ipt]=dev.statis[AVYY ][ijqp+1][ipt] 
                           +fac2*yv*yv;
      dev.statis[AVYZ ][ijqp+1][ipt]=dev.statis[AVYZ ][ijqp+1][ipt] 
                           +fac2*yv*zv;
      dev.statis[AVZZ ][ijqp+1][ipt]=dev.statis[AVZZ ][ijqp+1][ipt] 
                           +fac2*zv*zv;
      dev.statis[AEVX ][ijqp+1][ipt]=dev.statis[AEVX ][ijqp+1][ipt] 
                           +fac2*ee*xv;
      dev.statis[AEVY ][ijqp+1][ipt]=dev.statis[AEVY ][ijqp+1][ipt] 
                           +fac2*ee*yv;
//    ..right upper grid point
      //ibase=NSTAT*(ijqp+ngpx-1+MNGP*(ipt-1)) 
	  
      dev.statis[ARHO][ijqp+dev.ngpx][ipt]=dev.statis[ARHO][ijqp+dev.ngpx][ipt] 
                               +fac3;
      dev.statis[AVX ][ijqp+dev.ngpx][ipt]=dev.statis[AVX ][ijqp+dev.ngpx][ipt] 
                               +fac3*xv;
      dev.statis[AVY ][ijqp+dev.ngpx][ipt]=dev.statis[AVY ][ijqp+dev.ngpx][ipt] 
                               +fac3*yv;
      dev.statis[AVZ ][ijqp+dev.ngpx][ipt]=dev.statis[AVZ ][ijqp+dev.ngpx][ipt] 
                               +fac3*zv;
      dev.statis[AVV ][ijqp+dev.ngpx][ipt]=dev.statis[AVV ][ijqp+dev.ngpx][ipt] 
                               +fac3*vv;
      dev.statis[AE  ][ijqp+dev.ngpx][ipt]=dev.statis[AE  ][ijqp+dev.ngpx][ipt] 
                               +fac3*ee;
      dev.statis[AII ][ijqp+dev.ngpx][ipt]=dev.statis[AII ][ijqp+dev.ngpx][ipt] 
                               +fac3*pii;
      dev.statis[AVXX ][ijqp+dev.ngpx][ipt]=dev.statis[AVXX ][ijqp+dev.ngpx][ipt] 
                           +fac3*xv*xv;
      dev.statis[AVXY ][ijqp+dev.ngpx][ipt]=dev.statis[AVXY ][ijqp+dev.ngpx][ipt] 
                           +fac3*xv*yv;
      dev.statis[AVXZ ][ijqp+dev.ngpx][ipt]=dev.statis[AVXZ ][ijqp+dev.ngpx][ipt] 
                           +fac3*xv*zv;
      dev.statis[AVYY ][ijqp+dev.ngpx][ipt]=dev.statis[AVYY ][ijqp+dev.ngpx][ipt] 
                           +fac3*yv*yv;
      dev.statis[AVYZ ][ijqp+dev.ngpx][ipt]=dev.statis[AVYZ ][ijqp+dev.ngpx][ipt] 
                           +fac3*yv*zv;
      dev.statis[AVZZ ][ijqp+dev.ngpx][ipt]=dev.statis[AVZZ ][ijqp+dev.ngpx][ipt] 
                           +fac3*zv*zv;
      dev.statis[AEVX ][ijqp+dev.ngpx][ipt]=dev.statis[AEVX ][ijqp+dev.ngpx][ipt] 
                           +fac3*ee*xv;
      dev.statis[AEVY ][ijqp+dev.ngpx][ipt]=dev.statis[AEVY ][ijqp+dev.ngpx][ipt] 
                           +fac3*ee*yv;
//    ..right lower grid point
      //ibase=NSTAT*(ijqp+ngpx+1-1+MNGP*(ipt-1)) 
	  
      dev.statis[ARHO][ijqp+dev.ngpx+1][ipt]=dev.statis[ARHO][ijqp+dev.ngpx+1][ipt] 
                               +fac4;
      dev.statis[AVX ][ijqp+dev.ngpx+1][ipt]=dev.statis[AVX ][ijqp+dev.ngpx+1][ipt] 
                               +fac4*xv;
      dev.statis[AVY ][ijqp+dev.ngpx+1][ipt]=dev.statis[AVY ][ijqp+dev.ngpx+1][ipt] 
                               +fac4*yv;
      dev.statis[AVZ ][ijqp+dev.ngpx+1][ipt]=dev.statis[AVZ ][ijqp+dev.ngpx+1][ipt] 
                               +fac4*zv;
      dev.statis[AVV ][ijqp+dev.ngpx+1][ipt]=dev.statis[AVV ][ijqp+dev.ngpx+1][ipt] 
                               +fac4*vv;
      dev.statis[AE  ][ijqp+dev.ngpx+1][ipt]=dev.statis[AE  ][ijqp+dev.ngpx+1][ipt] 
                               +fac4*ee;
      dev.statis[AII ][ijqp+dev.ngpx+1][ipt]=dev.statis[AII ][ijqp+dev.ngpx+1][ipt] 
                               +fac4*pii;
      dev.statis[AVXX ][ijqp+dev.ngpx+1][ipt]=dev.statis[AVXX ][ijqp+dev.ngpx+1][ipt] 
                           +fac4*xv*xv;
      dev.statis[AVXY ][ijqp+dev.ngpx+1][ipt]=dev.statis[AVXY ][ijqp+dev.ngpx+1][ipt] 
                           +fac4*xv*yv;
      dev.statis[AVXZ ][ijqp+dev.ngpx+1][ipt]=dev.statis[AVXZ ][ijqp+dev.ngpx+1][ipt] 
                           +fac4*xv*zv;
      dev.statis[AVYY ][ijqp+dev.ngpx+1][ipt]=dev.statis[AVYY ][ijqp+dev.ngpx+1][ipt] 
                           +fac4*yv*yv;
      dev.statis[AVYZ ][ijqp+dev.ngpx+1][ipt]=dev.statis[AVYZ ][ijqp+dev.ngpx+1][ipt] 
                           +fac4*yv*zv;
      dev.statis[AVZZ ][ijqp+dev.ngpx+1][ipt]=dev.statis[AVZZ ][ijqp+dev.ngpx+1][ipt] 
                           +fac4*zv*zv;
      dev.statis[AEVX ][ijqp+dev.ngpx+1][ipt]=dev.statis[AEVX ][ijqp+dev.ngpx+1][ipt] 
                           +fac4*ee*xv;
      dev.statis[AEVY ][ijqp+dev.ngpx+1][ipt]=dev.statis[AEVY ][ijqp+dev.ngpx+1][ipt] 
                           +fac4*ee*yv;
//____energy distribution funktion
      //ibase=NEREDF*(ijqp-1+MNGP*(ipt-1))        
      dev.edf[ier][ijqp][ipt]=dev.edf[ier][ijqp][ipt]+fac1;
      //ibase=NEREDF*(ijqp+1-1+MNGP*(ipt-1))        
      dev.edf[ier][ijqp+1][ipt]=dev.edf[ier][ijqp+1][ipt]+fac2;
      //ibase=NEREDF*(ijqp+ngpx-1+MNGP*(ipt-1))        
      dev.edf[ier][ijqp+dev.ngpx][ipt]=dev.edf[ier][ijqp+dev.ngpx][ipt]+fac3;
      //ibase=NEREDF*(ijqp+ngpx+1-1+MNGP*(ipt-1))        
      dev.edf[ier][ijqp+dev.ngpx+1][ipt]=dev.edf[ier][ijqp+dev.ngpx+1][ipt]+fac4;
//____end of GETSTAT
      return;
}