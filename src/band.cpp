class Band{
//	C     colision broadening kloem
      double kloem[201][201];
//C     colision broadening kloab
      double kloab[201][201];
//C     colision broadening klaem
      double klaem[201][201];
//C     colision broadening klaab
      double klaab[201][201];
//C     colision broadening ktoem
      double ktoem[201][201];
//C     colision broadening ktoab
      double ktoab[201][201];
		 int mapsym[8][48];
	//    contains Name for each particle TYPE
      //CHARACTER*8 Typename[NPARTYP]
	  char Typename[NPARTYP][9];
//     number of scattering processes (electrons)
       int scpre;

//     number of scattering processes (holes)
       int scprh;

//     number of scattering processes (oxide electrons)
       int scproe;


//    field with 48 Symetrie-operations
       int matsym[6][48];

//    inverse of matsym
      // int indmat[-1:1,-1:1,-1:1,3,3,3]
	   int indmat[3][3][3][3][3][3];

//    matsym with a scalar argument
      //int matsymlin[6*48]

//    number tetraeder in in the energy-fields MWLE
       int ntlist[MWLE][MNB];

//    pointer to tlist
       int ptlist[MWLE][MNB];
//     particle type of the energy bands
      int partyp[MNB];


//    sequentiell list of tedraeder
       int tlist[MWLI];

//    number of surface tetraeder in in the energy-fields MWLE
       int nslist[MWLE][MNB];

//    pointer to slist
       int pslist[MWLE][MNB];

//    sequentiell list of surface tedraeder
       int slist[MSLI];

//    number of surface tetraeder in the energy-fields MWLE
       int nzlist[MNB];

//    pointer to slist
       int pzlist[MNB];

//    sequentiell list of surface tedraeder
       int zlist[MZLI];

//    number of tetraeder in the energy-fields MWLES
//    [finer energy spacing than tlist]
       int ntlists[MWLES][MNB];

//    pointer to tlists
       int ptlists[MWLES][MNB];

//    sequentiell list of tedraedra
       int tlists[MWLIS];

//    number of tetraeder in each cube
       int nclist[MCLE];

//    pointer to clist
       int pclist[MCLE];

//    sequentiell list of tetrahedra
       int clist[MCLI];

//    number tetraeder in the energy-fields MWLEFP [Fischetti, absorbtion]
       int ntlistfpab[MWLEFP][MNB];

//    pointer to tlistfpab [Fischetti, absorbtion]
       int ptlistfpab[MWLEFP][MNB];

//    sequentiell list of tedraeder [Fischetti, absorbtion]
       int tlistfpab[MWLIFP];

//    number tetraeder in the energy-fields MWLEFP [Fischetti, emisson]
       int ntlistfpem[MWLEFP][MNB];

//    pointer to tlistfpab [Fischetti, emisson]
       int ptlistfpem[MWLEFP][MNB];

//    sequentiell list of tedraeder [Fischetti, emisson]
       int tlistfpem[MWLIFP];

//    contains total number of bands
       int nbt;

//    contains number of bands for each particle type
       int nband[NPARTYP];

//    OFfset of BAND index
       int bandof[NPARTYP];

//    contains material type for each particle TYPE
       int typemat[NPARTYP];

//    indices of the for nodes of each tetrahedron
       int tet[4][MNT];

//    number of k-space grid points
       int nk;

//    number of tetraheda
       int nt;

//    energy band index of each tetrahdron
       int ibt[MNT];

//    list of neighbour tetraheda
       int nlt[4][MNT];

//    in the case of a change of the wedge newisym gives the new symmetry trans.
       int newisym[5][48];
//     flag for writing the scattering files to a rgraph file
       bool opcalc;

//     flag for Brooks-Herring impurity scattering
      bool bhfl;
//     flag for Sano formula for electron impact ionization rate
       bool sanoiifl;
//     flag for use of the inverse microscopic relaxation time instead 
//     of scattering rate
       bool bhmrtfl;///taf.f scatt.f

//     flag for Kane formula for electron impact ionization rate
      bool kaneiifl;

//     flag for Thoma formula for electron impact ionization rate
      bool thomiifl;

//     flag for Fischetti formula for electron impact ionization rate
      bool fisciifl;

//     flag for impact ionization scattering
       bool iifl;

//    evaluate the inverse mass tensor
       bool massfl;

//    self scattering event
	  // bool selfscfl;
//     flag for generation of secondary particles
      bool seciifl;

//     flag for tunneling of electrons in the oxide (only negative x-direction)
      bool tunxfl;
//     flag for Jacoboni phonon system
       bool jacophfl;

//     falg for fischetti phonon system
       bool fiscphfl;

//     empirical correction factor for BH impurity scattering
      bool frickfl;

//     flag for particle transition from silicon into oxide
       bool injoxfl;
//    energy of the base of the tetrahedron
       double ebzp;

//    groupvelocity within the tetrahedron [irreducible wedge]
       double vgx,vgy,vgz;

//    maximum of density of states for each particle type
       double DOSMAX[NPARTYP];

//    table of density of states for each band
       double dos[MNB][MTAB+1];

//    table of the sum of DOS over all bands for each particle type
       double sumdos[MTAB+1][NPARTYP];

//    DOS of all tetrahedra in the list, each tetrahedron
//    contributing with its maximum DOS [Fischetti, absorbtion]
       double dosfpab[MNB][MWLEFP];

//    DOS of all tetrahedra in the list, each tetrahedron [Fischetti, emisson]
//    contributing with its maximum DOS
       double dosfpem[MNB][MWLEFP];

//    maximum AOV for list [Fischetti, aborbtion]
       double maxaovfpab[MNB][MWLEFP];

//    maximum AOV for list [Fischetti, emisson]
       double maxaovfpem[MNB][MWLEFP];

//    energy for tables [scatt and DOS]
       double energy[MTAB+1];

//    energy spacing of tables
       double dtable;

//    minimum energy in table
       double emin;

//    maximum energy in table
       double emax;

//    energy spacing for treaheder list
       double dlist;

//    energy spacing for treaheder list
       double dlists;

//    energy spacing for treaheder list [Fischetti]
       double dlistfp;

//    array of grid points in k-space
       double xkk[MNK],ykk[MNK],zkk[MNK];

//    array of energy of grid points in k-space
       double eek[MNK];

//    array of comp. of velocity in tetrahedron
       double vt[4*MNT];

//    absolute value of the group velocity
       double vgt[MNT];

//    center coordinates of each tetrahedron
       double xkct[MNT],ykct[MNT],zkct[MNT];

//    maximum dos in tet. for each band and energy
       double dostetmax[MNB][MTAB+1];

//    maximum area over velocity for each tet.
       double maxaovtet[MNT];
//    inverse mass [1:xx, 2:xy, 3:xz, 4:yy, 5:yz, 6:zz]
//    centered in the middle of the tetrahedron
       double massinv[6][MNT];
//    array of normal vectors of tetrahedra surfaces
//    first index = tet. surface
//    datantlin[1,*,*] = x-comp
//    datantlin[2,*,*] = y-comp
//    datantlin[3,*,*] = z-comp
//    datantlin[4,*,*] = distances of tet. surfaces from origin
      //double datantlin[4][4][MNT];!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
//    same as datant with a scalar argument
       double datantlin[4][4][MNT]; 
//    factor for FASTSURF for energy between e1 and e2
       double faclow[MNT];
//    factor for FASTSURF for energy between e2 and e3
       double facmedium[MNT];
//    factor for FASTSURF for energy between e3 and e4
       double fachigh[MNT];

//     prefactor for hole II
       double iifachole;
//     hole II threshold energy
       double hiithresh;
//     hole II energy exponent
       double hiiexp;
//     enhancement factor for e-P scattering except for intraband scattering 
//     in the first conduction band
       double ephb;
//     oxide electron mass
       double dmox;
//     expectation value of the trace of the inverse mass tensor for equilibrium
       double dmc[NPARTYP];
//     longitudinal mass in the minimum of the first conduction band
       double mell;
//    transversal mass in the minimum of the first conduction band
       double melt;
//    DOS mass in the minimum of the first conduction band
       double meld;
//    temperature of the transversal acoustic g-phonon (electrons)
       double temptag;
//    temperature of the longitudinal acoustic g-phonon (electrons)
       double templag;
//    temperature of the longitudinal optical g-phonon (electrons)
       double templog;
//    temperature of the transversal acoustic f-phonon (electrons)
       double temptaf;
//    temperature of the longitudinal acoustic f-phonon (electrons)
       double templaf;
//    temperature of the transversal optical f-phonon (electrons)
       double temptof;
//    temperature of the optical phonon (holes)
       double temphop;
//    temperature of the optical phonon (oxide electrons)
       double tempoeop;
//    deformation potential constant of the transversal acoustic g-phonon 
//    (electrons)
       double dftag;
//    deformation potential constant of the longitudinal acoustic g-phonon 
//    (electrons)
       double dflag;
//    deformation potential constant of the longitudinal optical g-phonon 
//    (electrons)
       double dflog;
//    deformation potential constant of the transversal acoustic f-phonon 
//    (electrons)
       double dftaf;
//    deformation potential constant of the longitudinal acoustic f-phonon 
//    (electrons)
       double dflaf;
//    deformation potential constant of the transversal optical f-phonon 
//    (electrons)
       double dftof;
//    deformation potential for elastic acoustic phonon scattering (electrons)
       double dfelast;
//    deformation potential constant for optical phonon scattering (holes)
       double dfhop;
//    deformation potential for elastic acoustic phonon scattering (holes)
       double dfhelast;
//    deformation potential constant for optical phonon scattering 
//    (oxide electrons)
       double dfoeop;
//    deformation potential for elastic acoustic phonon scattering 
//    (oxide electrons)
       double dfoeelast;
//    prefactor for electron impact ionization
       double iifacelec;

//    coupling constant for optical phonons low band
       double efoplow;

//    coupling constant for optical phonons high band
       double efophigh;

//    energy for optical phonons 
       double efopee;

//    coupling constant for optical phonons low band
       double efaplow;

//    coupling constant for optical phonons high band
       double efaphigh;

//    energy for transversal acoustic phonons 
       double efapeet;

//    energy for longitudinal acoustic phonons 
       double efapeel;
//     intrinsic carrier concentration
       double eni;
//    lattice vectors for change of BZ
       double xbz[5][48],ybz[5][48],zbz[5][48];
//     maximaum scattering rate for particle scattering
       double gamma[NPARTYP];
//     maximum scattering rate in a tetrahedron
       double gamtet[MNT];

//     fraction of diffusive scattering at Si/SiO2 interface
       double difpr[NPARTYP];
//     phonon scattering rate prefactor table for electrons
       double scatte[MSCPRE][NBE][NBE];

//     scattering rate table for electron impact ionization
       double scattiie[NBE][NBE][MTAB+1]; 

//     phonon scattering rate prefactor table for holes
       double scatth[MSCPRH][NBH][NBH];

//     scattering rate table for hole impact ionization
       double scattiih[NBH][NBH][MTAB+1];

//     phonon scattering rate prefactor table for oxide electrons
       double scattoe[MSCPROE][NBOE][NBOE];

//     DOS for the different phonon scattering processes final energies 
//     (electrons)
       double dose[MSCPRE][NBE][MTAB+1];

//     DOS for the different phonon scattering processes final energies (holes)
       double dosh[MSCPRH][NBH][MTAB+1];

//     DOS for the different phonon scattering processes final energies 
//     (oxide electrons)
       double dosoe[MSCPROE][NBOE][MTAB+1];

//     sum of scattering rates for each band (II + Phonon)
       double sumscatt[MTAB+1][MNB];

	  void BUILDPHSCATT(void);
	  double CALDOSTETMAX(double ee,int iband);
	  void CALSCATTE(double ee,double*sca,int iband);
	  void CALSCATTH(double ee,double*sca, int iband);
	  void CALSCATTOE(double ee,double*sca,int iband);
	  double CALSCATTSUM(double ee,int iband);
	  void CALCFAC(void);
	  double CALDOS(double eed,int ib);
	  double CALDOSSUM(double eed,int ip);
	  void CONFBZ(double &xkl,double &ykl,double &zkl);
	  double EIIRATE(double ee);
	  double FASTSURF(double ee,int itet);
	  double FII(double x1,double y1,double z1,
				  double x2,double y2,double z2,
				  double x3,double y3,double z3);
	  void GETMAXAOVTET(void);
	  void GETINDENS(void);
	  double HIIRATE(double ee);
	  void HESSTET(void);
	  void SETNEIBSYM(void);
	  void SETMATSYM(void);
	  double SURF(double eps,int it);
	  void VELOTET(void);
	  void ZUDI(void);
	  void OUTWEDGE(double &xout,double &yout,double &zout,
		  double xin,double yin,double zin,int matsym[6][48],int isym);
	  void OVERLAP(int it,int&nct,int*ict);
	  void  READBS(void);
	  void GETDOS(bool calcdos);
	  void BUILDLISTS(void);
public :
	friend class Partical;
	friend class DevSimulator;
	//friend void DevSimulator::GETIMAGEPOT(Band &bd);
//    silicon to silicon oxide conduction band discontinuity
      double sioxbgo;
//     image potential factor for Schottky barrier lowering
      double beta;	
	void IELEC(void);
};
void Band::VELOTET(void)
{
//____Purpose : Calculate velocities within each tetrahedron
//              
//____local variables
      int it,ibase;
      double m[4][4],ar[4],vv,mi[4][4],det;
	  double ee2,ee3,ee4;
//____loop over all tetrahedra
      for(it=0;it<nt;it++)
	  {
	  //DO it=1,nt
//      ..build matrix
         m[1][1]=xkk[tet[1][it]]-xkk[tet[0][it]];         
         m[1][2]=ykk[tet[1][it]]-ykk[tet[0][it]];        
         m[1][3]=zkk[tet[1][it]]-zkk[tet[0][it]];         
         m[2][1]=xkk[tet[2][it]]-xkk[tet[0][it]];         
         m[2][2]=ykk[tet[2][it]]-ykk[tet[0][it]];         
         m[2][3]=zkk[tet[2][it]]-zkk[tet[0][it]];         
         m[3][1]=xkk[tet[3][it]]-xkk[tet[0][it]];         
         m[3][2]=ykk[tet[3][it]]-ykk[tet[0][it]];         
         m[3][3]=zkk[tet[3][it]]-zkk[tet[0][it]];         
//      ..determinat
         det=m[1][1]*(m[2][2]*m[3][3]-m[3][2]*m[2][3]) 
           -m[1][2]*(m[2][1]*m[3][3]-m[3][1]*m[2][3])
           +m[1][3]*(m[2][1]*m[3][2]-m[3][1]*m[2][2]);
//      ..inverse matrix
         mi[1][1]=(m[2][2]*m[3][3]-m[3][2]*m[2][3])/det;
         mi[2][1] =-(m[2][1]*m[3][3]-m[3][1]*m[2][3])/det;
         mi[3][1]=(m[2][1]*m[3][2]-m[3][1]*m[2][2])/det;
         mi[1][2] =-(m[1][2]*m[3][3]-m[3][2]*m[1][3])/det;
         mi[2][2]=(m[1][1]*m[3][3]-m[3][1]*m[1][3])/det;
         mi[3][2] =-(m[1][1]*m[3][2]-m[3][1]*m[1][2])/det;
         mi[1][3]=(m[1][2]*m[2][3]-m[2][2]*m[1][3])/det;
         mi[2][3] =-(m[1][1]*m[2][3]-m[2][1]*m[1][3])/det;
         mi[3][3]=(m[1][1]*m[2][2]-m[2][1]*m[1][2])/det;
         ar[1]=eek[tet[1][it]]-eek[tet[0][it]];
         ar[2]=eek[tet[2][it]]-eek[tet[0][it]];
         ar[3]=eek[tet[3][it]]-eek[tet[0][it]];
//      ..group velocity within tetrahedron
         ibase=4*it;// 4*(it-1)  it changed
         vt[ibase+0]=mi[1][1]*ar[1]+mi[1][2]*ar[2]+mi[1][3]*ar[3];
         vt[ibase+1]=mi[2][1]*ar[1]+mi[2][2]*ar[2]+mi[2][3]*ar[3];
         vt[ibase+2]=mi[3][1]*ar[1]+mi[3][2]*ar[2]+mi[3][3]*ar[3];
         vt[ibase+3]=eek[tet[0][it]];
//      ..absolute value of group velocity
         vv=sqrt(vt[ibase+0]*vt[ibase+0]+vt[ibase+1]*vt[ibase+1]+vt[ibase+2]*vt[ibase+2]);
         vgt[it]=vv;
	  } 
//____End of VELOTET
      return;
}
void Band::ZUDI(void)
{
//________ Purpose   : calculate the DOS-table and maximum DOS in tet.
//____local variables
      int itab,it,ib,ie,itl,is;
      double eps,area,epsu,epsm,epsl;
      double  areal,areau;
//____set all dos-values
	  for(itab=0;itab<=MTAB;itab++)
	  {
      //DO itab=0,MTAB
		  for(ib=0;ib<MNB;ib++)
		  {
         //DO ib=1,MNB
            dos[ib][itab]=0.0;
            dostetmax[ib][itab]=0.0;
		  }
	  }
//____first: claculate DOS
//____loop over all energy-tab-steps
      for(itab=0;itab<=MTAB;itab++)
	  {
	  //DO itab=0,MTAB
         eps=energy[itab];
//_______energy index of tetraheda list
         ie=eps*dlist;//+1//minim=0;
         if(ie<0)ie=0;//1;
         if(ie>=MWLE)ie=MWLE-1;//maxim=
//_______cal. DOS for each tetrahedron
//      ..loop over all tetraheda,which may contain the energy eps
         for(ib=0;ib<nbt;ib++)
		 {
		 //DO ib=1,nbt
			for(itl=ptlist[ie][ib];itl<=ptlist[ie][ib]+ntlist[ie][ib]-1;itl++)
			{
        // DO itl=
				it=tlist[itl];
//         ..check if the tetrahedron contains the energy eps
				if(eps>eek[tet[0][it]]&&eps<eek[tet[3][it]])
				{
//             ..calculate area of intersection
					area=FASTSURF(eps,it);
//             ..cal. density of states
					dos[ibt[it]][itab]+=area/vgt[it]*48.0/(8.0*PI*PI*PI);//**3
				}
            }
         }
//____end of energy-loop
      }

//____cal. maximum dos in tetrahedron
//____loop over all energy-tab-steps

	  for(itab=1;itab<=MTAB;itab++)///itab=0 is set to 0;
	  {
	  //DO itab=1,MTAB
//_______cal. maximum area for each tetrahedron between the energies
//       epsl and epsu
//      .. loop over all tetraheda
         for(it=0;it<nt;it++)
		 {
		 //DO it=1,nt
            epsl=energy[itab-1];
            epsu=energy[itab];
//         ..check,if the tetrahedron overlapps with the given energy range
            if((epsu>eek[tet[0][it]])&&(epsl<eek[tet[3][it]]))
			{
               if(epsu<=eek[tet[1][it]])
//            ..the area grows quadratically between the energy of
//              node 1 and 2
                  area=FASTSURF(epsu,it);
               else if(epsl>=eek[tet[2][it]])
//            ..the area decreases quadratically between the energy of
//              node 3 and 4
                  area=FASTSURF(epsl,it);
               else
			   {
//               ..between the energy of node 2 and 3 the area has a maximum
//                 the maximum in the interval between epsl and epsu is found
//                 with bisection (it has been assumed,that 64 iterations are 
//                 sufficient).
//               ..limit the energy interval to the range between the energy
//                 of node 2 and 3
                  epsl=MAX(epsl,eek[tet[1][it]]);
                  epsu=MIN(epsu,eek[tet[2][it]]);
                  areal=FASTSURF(epsl,it);
                  areau=FASTSURF(epsu,it);
                  for(is=0;is<64;is++)
				  {
				  //DO is=1,64
                     epsm=0.50*(epsl+epsu);
                     if(areau>areal)
					 {
                        epsl=epsm;
                        areal=FASTSURF(epsl,it);
					 }
                     else
					 {
                        epsu=epsm;
                        areau=FASTSURF(epsu,it);
					 }
                  }
                  area=MAX(areal,areau);
               }
//            ..take maximum
               dostetmax[ibt[it]][itab]=MAX(area/vgt[it],
                                          dostetmax[ibt[it]][itab]);
            }
         }
//____end of energy-loop
      }
//____set maximum area for zero energy
      for(ib=0;ib<MNB;ib++)
	  // DO ib=1,MNB
         dostetmax[ib][0]=dostetmax[ib][1];
//____End of ZUDI
      return;
}
       
//====
double Band::SURF(double eps,int it)
{
 //     double eps
  //    int it
//    : Calculate the area of the intersection of an equi energy plane with
//      a tetrahedron
//
//    : between the energies of the node 1 and 2 the intersection is always
//      of triangular shape
//      between node 2 and 3 always a quadrangle and between 3 and 4 again
//      a triangle
//____local variables
      double rfac,ax,ay,az,bx,by,bz,cx,cy,cz;
      double dx,dy,dz,a1,a2,a3,a4;//,FII;
	  double rtnv;
      if(eps>eek[tet[0][it]]&&eps<eek[tet[3][it]])
	  {
//_______triangular intersection with energy lower e2
         if (eps<=eek[tet[1][it]])
		 {
//         ..cal. vectors (relative to origin of tetrahedron) pointing to the 
//           three nodes of the triangle
            rfac= (eps-eek[tet[0][it]])/(eek[tet[1][it]]-eek[tet[0][it]]);
            ax=rfac*(xkk[tet[1][it]]-xkk[tet[0][it]]);
            ay=rfac*(ykk[tet[1][it]]-ykk[tet[0][it]]);
            az=rfac*(zkk[tet[1][it]]-zkk[tet[0][it]]);
            if (rfac<0.0||rfac>1.0)
			{
				cout<<"error SURF: rfac< 0 or>1";exit(0);
			}
			rfac=(eps-eek[tet[0][it]])/(eek[tet[2][it]]-eek[tet[0][it]]);
            bx=rfac*(xkk[tet[2][it]]-xkk[tet[0][it]]);
            by=rfac*(ykk[tet[2][it]]-ykk[tet[0][it]]);
            bz=rfac*(zkk[tet[2][it]]-zkk[tet[0][it]]);
            if (rfac<0.0||rfac>1.0)
			{
				cout<<"error SURF: rfac< 0 or>1";exit(0);
			}
            rfac=(eps-eek[tet[0][it]])/(eek[tet[3][it]]-eek[tet[0][it]]);
            cx=rfac*(xkk[tet[3][it]]-xkk[tet[0][it]]);
            cy=rfac*(ykk[tet[3][it]]-ykk[tet[0][it]]);
            cz=rfac*(zkk[tet[3][it]]-zkk[tet[0][it]]);
            if (rfac<0.0||rfac>1.0)
			{
				cout<<"error SURF: rfac< 0 or>1";exit(0);
			}
//         ..cal. area
            rtnv=FII(ax,ay,az,bx,by,bz,cx,cy,cz);
		 }
      
//_______quadrangular intersection with energy between e2 and e3
         else if(eps>eek[tet[1][it]]&&eps<eek[tet[2][it]])
		 {
//         ..cal. vectors (relative to origin of k-space) pointing to the 
//           four nodes of the quadrangle
            rfac=(eps-eek[tet[0][it]])/(eek[tet[2][it]]-eek[tet[0][it]]);
            ax=rfac*(xkk[tet[2][it]]-xkk[tet[0][it]])+xkk[tet[0][it]];
            ay=rfac*(ykk[tet[2][it]]-ykk[tet[0][it]])+ykk[tet[0][it]];
            az=rfac*(zkk[tet[2][it]]-zkk[tet[0][it]])+zkk[tet[0][it]];
            if (rfac< 0.0||rfac> 1.0)
			{
			   cout<<"error SURF: rfac< 0 or>1";exit(0);
			}
            rfac=(eps-eek[tet[1][it]])
               /(eek[tet[2][it]]-eek[tet[1][it]]);
            bx=rfac*(xkk[tet[2][it]]-xkk[tet[1][it]])+xkk[tet[1][it]];
            by=rfac*(ykk[tet[2][it]]-ykk[tet[1][it]])+ykk[tet[1][it]];
            bz=rfac*(zkk[tet[2][it]]-zkk[tet[1][it]])+zkk[tet[1][it]];
            if (rfac< 0.0||rfac> 1.0)
			{
				cout<<"error SURF: rfac< 0 or>1";exit(0);
			}
            rfac= (eps-eek[tet[0][it]])
               /(eek[tet[3][it]]-eek[tet[0][it]]);
            cx=rfac*(xkk[tet[3][it]]-xkk[tet[0][it]])+xkk[tet[0][it]];
            cy=rfac*(ykk[tet[3][it]]-ykk[tet[0][it]])+ykk[tet[0][it]];
            cz=rfac*(zkk[tet[3][it]]-zkk[tet[0][it]])+zkk[tet[0][it]];
            if (rfac< 0.0||rfac> 1.0)
			{
				cout<<"error SURF: rfac< 0 or>1";exit(0);
			}
			rfac= (eps-eek[tet[1][it]])
               /(eek[tet[3][it]]-eek[tet[1][it]]);
            dx=rfac*(xkk[tet[3][it]]-xkk[tet[1][it]])+xkk[tet[1][it]];
            dy=rfac*(ykk[tet[3][it]]-ykk[tet[1][it]])+ykk[tet[1][it]];
            dz=rfac*(zkk[tet[3][it]]-zkk[tet[1][it]])+zkk[tet[1][it]];
            if (rfac< 0.0||rfac> 1.0)
			{
				cout<<"error SURF: rfac< 0 or>1";exit(0);
			}
//      ..cal. area of the four possible triangles made from the four nodes
//        of the quadrangle
            a1=FII(ax,ay,az,bx,by,bz,cx,cy,cz);
            a2=FII(dx,dy,dz,bx,by,bz,cx,cy,cz);
            a3=FII(ax,ay,az,dx,dy,dz,cx,cy,cz);
            a4=FII(ax,ay,az,bx,by,bz,dx,dy,dz);
            rtnv=0.50*(a1+a2+a3+a4);
            
		 }      
//_______triangular intersection with energy higher e3
         else if(eps>eek[tet[2][it]]&&eps<eek[tet[3][it]])
		 {
//         ..cal. vectors (relative to origin of tetrahedron) pointing to the 
//           three nodes of the triangle
            rfac= (eps-eek[tet[3][it]])/(eek[tet[0][it]]-eek[tet[3][it]]);
            ax=rfac*(xkk[tet[0][it]]-xkk[tet[3][it]]);
            ay=rfac*(ykk[tet[0][it]]-ykk[tet[3][it]]);
            az=rfac*(zkk[tet[0][it]]-zkk[tet[3][it]]);
            if (rfac< 0.0||rfac> 1.0)
			{
				cout<<"error SURF: rfac< 0 or>1";exit(0);
			}
            rfac= (eps-eek[tet[3][it]])/(eek[tet[1][it]]-eek[tet[3][it]]);
            bx=rfac*(xkk[tet[1][it]]-xkk[tet[3][it]]);
            by=rfac*(ykk[tet[1][it]]-ykk[tet[3][it]]);
            bz=rfac*(zkk[tet[1][it]]-zkk[tet[3][it]]);
            if (rfac< 0.0||rfac> 1.0)
			{
				cout<<"error SURF: rfac< 0 or>1";exit(0);
			}
            rfac= (eps-eek[tet[3][it]])/(eek[tet[2][it]]-eek[tet[3][it]]);
            cx=rfac*(xkk[tet[2][it]]-xkk[tet[3][it]]);
            cy=rfac*(ykk[tet[2][it]]-ykk[tet[3][it]]);
            cz=rfac*(zkk[tet[2][it]]-zkk[tet[3][it]]);
            if (rfac< 0.0||rfac> 1.0)
			{
				cout<<"error SURF: rfac< 0 or>1";exit(0);
			}
            rtnv=FII(ax,ay,az,bx,by,bz,cx,cy,cz);
		 }
         else
            rtnv=0.0;
      }
      else
         rtnv=0.0;
//____End of SURF
      return rtnv;
}
////====
double Band:: FII(double x1,double y1,double z1,
				  double x2,double y2,double z2,
				  double x3,double y3,double z3)
{
      //double x1,y1,z1,x2,y2,z2,x3,y3,z3,
	  double sbx,sby,sbz,rtnv;
//____Purpose : calculate the square in a triangle
//    Parameter : values of the three vectors pointing to the nodes
//    : The area is given by the absolute vale of the cross product
//      of the two vectors along the edges between the nodes 2/1 and 3/1
     
      sbx=(y2-y1)*(z3-z1)-(y3-y1)*(z2-z1);
      sby=(x3-x1)*(z2-z1)-(x2-x1)*(z3-z1);
      sbz=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
      rtnv=0.50*sqrt(sbx*sbx+sby*sby+sbz*sbz);
//____End of FII
      return rtnv;
}
//====
void Band::HESSTET(void)
{
//____Purpose : cal. Hessian form for the four surfaces of the tetraheda,
//              the normal vectors are not normalized,
//              they allways point into the tetrahedron.
//____local variables
      double ax,ay,az,bx,by,bz,cx,cy,cz;
      double dx,dy,dz,ex,ey,ez,fx,fy,fz;
      int it;
     
//____loop over all tetraheda
      for(it=0;it<nt;it++)
	  {
	  //DO it=1,nt
//      ..cal. the six vectors along the edges of the tetrahedron
         ax=xkk[tet[1][it]]-xkk[tet[0][it]];
         ay=ykk[tet[1][it]]-ykk[tet[0][it]];
         az=zkk[tet[1][it]]-zkk[tet[0][it]];
         bx=xkk[tet[2][it]]-xkk[tet[1][it]];
         by=ykk[tet[2][it]]-ykk[tet[1][it]];
         bz=zkk[tet[2][it]]-zkk[tet[1][it]];
         cx=xkk[tet[0][it]]-xkk[tet[2][it]];
         cy=ykk[tet[0][it]]-ykk[tet[2][it]];
         cz=zkk[tet[0][it]]-zkk[tet[2][it]];
         dx=xkk[tet[3][it]]-xkk[tet[0][it]];
         dy=ykk[tet[3][it]]-ykk[tet[0][it]];
         dz=zkk[tet[3][it]]-zkk[tet[0][it]];
         ex=xkk[tet[3][it]]-xkk[tet[1][it]];
         ey=ykk[tet[3][it]]-ykk[tet[1][it]];
         ez=zkk[tet[3][it]]-zkk[tet[1][it]];
         fx=xkk[tet[3][it]]-xkk[tet[2][it]];
         fy=ykk[tet[3][it]]-ykk[tet[2][it]];
         fz=zkk[tet[3][it]]-zkk[tet[2][it]];
//      .. surface 1
//      .. cal. vector perpendicular to the surface (cross product)
         datantlin[0][0][it]=ey*fz-ez*fy;
         datantlin[1][0][it]=ez*fx-ex*fz;
         datantlin[2][0][it]=ex*fy-ey*fx;
//      .. cal. distance of surface to origin of tetrahedron
         datantlin[3][0][it]=datantlin[0][0][it]*xkk[tet[1][it]] 
                 +datantlin[1][0][it]*ykk[tet[1][it]] 
                 +datantlin[2][0][it]*zkk[tet[1][it]]; 
//      .. force direction of perpendicular vector to point into the tetrahedron
         if ( datantlin[0][0][it]*xkk[tet[0][it]] 
          +datantlin[1][0][it]*ykk[tet[0][it]] 
          +datantlin[2][0][it]*zkk[tet[0][it]]<datantlin[3][0][it])
		 {
            datantlin[0][0][it] =-datantlin[0][0][it];
            datantlin[1][0][it] =-datantlin[1][0][it];
            datantlin[2][0][it] =-datantlin[2][0][it];
            datantlin[3][0][it] =-datantlin[3][0][it];
         }
            
//      .. surface 2
         datantlin[0][1][it]=dy*fz-dz*fy;
         datantlin[1][1][it]=dz*fx-dx*fz;
         datantlin[2][1][it]=dx*fy-dy*fx;
         datantlin[3][1][it]=datantlin[0][1][it]*xkk[tet[2][it]] 
                 +datantlin[1][1][it]*ykk[tet[2][it]] 
                 +datantlin[2][1][it]*zkk[tet[2][it]]; 
         if ( datantlin[0][1][it]*xkk[tet[1][it]] 
          +datantlin[1][1][it]*ykk[tet[1][it]] 
          +datantlin[2][1][it]*zkk[tet[1][it]]<datantlin[3][1][it])
		 {
            datantlin[0][1][it] =-datantlin[0][1][it];
            datantlin[1][1][it] =-datantlin[1][1][it];
            datantlin[2][1][it] =-datantlin[2][1][it];
            datantlin[3][1][it] =-datantlin[3][1][it];
         }
//      .. surface 3
         datantlin[0][2][it]=ey*dz-ez*dy;
         datantlin[1][2][it]=ez*dx-ex*dz;
         datantlin[2][2][it]=ex*dy-ey*dx;
         datantlin[3][2][it]=datantlin[0][2][it]*xkk[tet[3][it]] 
                 +datantlin[1][2][it]*ykk[tet[3][it]] 
                 +datantlin[2][2][it]*zkk[tet[3][it]]; 
         if ( datantlin[0][2][it]*xkk[tet[2][it]] 
          +datantlin[1][2][it]*ykk[tet[2][it]] 
          +datantlin[2][2][it]*zkk[tet[2][it]]<datantlin[3][2][it])
		 {
            datantlin[0][2][it] =-datantlin[0][2][it];
            datantlin[1][2][it] =-datantlin[1][2][it];
            datantlin[2][2][it] =-datantlin[2][2][it];
            datantlin[3][2][it] =-datantlin[3][2][it];
         }
//      .. surface 4
         datantlin[0][3][it]=ay*bz-az*by;
         datantlin[1][3][it]=az*bx-ax*bz;
         datantlin[2][3][it]=ax*by-ay*bx;
         datantlin[3][3][it]=datantlin[0][3][it]*xkk[tet[0][it]] 
                  +datantlin[1][3][it]*ykk[tet[0][it]] 
                  +datantlin[2][3][it]*zkk[tet[0][it]] ;
         if ( datantlin[0][3][it]*xkk[tet[3][it]] 
           +datantlin[1][3][it]*ykk[tet[3][it]] 
           +datantlin[2][3][it]*zkk[tet[3][it]]<datantlin[3][3][it])
		 {
            datantlin[0][3][it] =-datantlin[0][3][it];
            datantlin[1][3][it] =-datantlin[1][3][it];
            datantlin[2][3][it] =-datantlin[2][3][it];
            datantlin[3][3][it] =-datantlin[3][3][it];
		 }
	  }
         
//____End of HESSTET
      return ;
}

//=====

void Band::SETNEIBSYM(void)
{
//____Purpose : cal. the five neighbouring wedges for each of the 48 wedges
//              
//____local variables
	  int jsym,iv,isym;
      int i,j,pos[4],poshilf[4];//(3)
      double ksort[4],ax,ay,az,xkl,ykl,zkl;
      double xtrans[5],ytrans[5],ztrans[5];
	  double xk,yk,zk;
//____five vectors which are outside of the irreducible wedge
//    and within the five neighbour wedges
//   .. kz=0
      xtrans[0]= 0.5*a0pi;
      ytrans[0]= 0.3*a0pi;
      ztrans[0]=-0.01*a0pi;
//   .. kx-ky=0
      xtrans[1]= 0.49*a0pi;
      ytrans[1]= 0.51*a0pi;
      ztrans[1]= 0.1*a0pi;
//   .. ky-kz=0
      xtrans[2]= 0.5*a0pi;
      ytrans[2]= 0.24*a0pi;
      ztrans[2]= 0.26*a0pi;
//   .. kx+ky+kz=1.5
      xtrans[3]= 0.71*a0pi;
      ytrans[3]= 0.51*a0pi;
      ztrans[3]= 0.31*a0pi;
//   .. kx=1
      xtrans[4]= 1.01*a0pi;
      ytrans[4]= 0.1*a0pi;
      ztrans[4]= 0.05*a0pi;
//____loop over all 48 wedges of the first BZ
      for(isym=0;isym<48;isym++)
	  {
	  //DO isym=1,48
//      ..loop over all five neighbours
         for(iv=0;iv<=4;iv++)/////////becareful
		 {
		  //DO iv=1,5
            ax=xtrans[iv];
            ay=ytrans[iv];
            az=ztrans[iv];
            OUTWEDGE(xk,yk,zk,ax,ay,az,matsym,isym);
            xkl=xk;
            ykl=yk;
            zkl=zk;
//         ..in the case that the vector is outside of the first BZ
//           reduce it to the first BZ
            CONFBZ(xk,yk,zk);
//         .. save k-vector displacement
            xbz[iv][isym]=xk-xkl;
            ybz[iv][isym]=yk-ykl;
            zbz[iv][isym]=zk-zkl;
            if(fabs(xbz[iv][isym])<1e-10) xbz[iv][isym]=0.0;
            if(fabs(ybz[iv][isym])<1e-10) ybz[iv][isym]=0.0;
            if(fabs(zbz[iv][isym])<1e-10) zbz[iv][isym]=0.0;
//         .. cal. index of the new wedge
            ksort[1]=fabs(xk);
            ksort[2]=fabs(yk);
            ksort[3]=fabs(zk);
            pos[1]=1;
            pos[2]=2;
            pos[3]=3;
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
		  /*
            DO i=1,2
               DO j=2,3
                  if (ksort(i)<ksort(j)) THEN
      
                    ksort(0)=ksort(i)
                    ksort(i)=ksort(j)
                    ksort(j)=ksort(0)
      
                    pos(0)=pos(i)
                    pos(i)=pos(j)
                    pos(j)=pos(0)
      
                  ENDif
               ENDDO
            ENDDO
			*/
			poshilf[pos[1]]=1-1;///indmax[3][3][3][3][3][3]
			poshilf[pos[2]]=2-1;
			poshilf[pos[3]]=3-1;
            //poshilf(pos(1))=1
            //poshilf(pos(2))=2
            //poshilf(pos(3))=3
            //jsym=indmat[sign(xk)),NINT(SIGN(1.d0,yk)),
             //              NINT(SIGN(1.d0,zk)),
             //              poshilf((1)),poshilf((2)),poshilf((3)) )
		    jsym=indmat[sign(xk)+1][sign(yk)+1][sign(zk)+1][poshilf[1]][poshilf[2]][poshilf[3]];
            newisym[iv][isym]=jsym;
         }
	  }
          
         
//____End of SETNEIBSYM
      return;
}
void Band::OUTWEDGE(double &xout,double &yout,
		double &zout,double xin,double yin,double zin,int matsym[6][48],int isym)
{

//transform k-vector from the irreducible wedge into the BZ

//_____local variables
      double  vsort[3];
	  int ibase;

      vsort[0]=xin;
      vsort[1]=yin;
      vsort[2]=zin;

      ibase=isym;
      xout=vsort[matsym[0][ibase]-1]*matsym[3][ibase];
      yout=vsort[matsym[1][ibase]-1]*matsym[4][ibase];
      zout=vsort[matsym[2][ibase]-1]*matsym[5][ibase];

//_____End of OUTWEDGE
      return;
}
void Band::CONFBZ (double &xkl,double &ykl,double &zkl)
 {
	 //     Purpose: confine k-vector to first BZ
//     -------
//____local variables
    bool flag,ndfl;

//____find the equivalent k-vektor in the 
//     first brillouin zone
	ndfl=false;
    while(!ndfl)
	{
	  flag=false;
      if(fabs(xkl)>a0pi) 
               flag=true;
      else if (fabs(ykl)>a0pi)
               flag=true;
      else if (fabs(zkl)>a0pi)
               flag=true;
      else if (fabs(xkl)+fabs(ykl)+fabs(zkl)>1.50*a0pi)
               flag=true;
     
      if(flag)
	  {
         if (xkl>0.0)
            xkl=xkl-a0pi;
         else
            xkl=xkl+a0pi;
         if (ykl>0.0)
            ykl=ykl-a0pi;
         else
            ykl=ykl+a0pi;
         if (zkl>0.0)
            zkl=zkl-a0pi;
         else
            zkl=zkl+a0pi;
      }
      else
	     ndfl=true;   
    }//GOTO 10

//____end of CONFBZ
      return;
 }
void Band::CALCFAC(void)
{
//____Purpose : Calculate factors for area calculation within each tetrahedron
//              used in FASTSURF
//              
//    : In the energy intervals between the energy of the nodes the
//      dependence of the area on energy is given by a quadrati//relation
//____local variables
      int it;
      double eps;
//____loop over all tetrahedra
      for(it=0;it<nt;it++)
	  {
	  //DO it=1,nt
//      .. between e1 and e2 the intersection is triangular
         if (eek[tet[1][it]]>eek[tet[0][it]])
		 {
            eps=0.50*(eek[tet[0][it]]+eek[tet[1][it]]);
            faclow[it]=SURF(eps,it)/((eps-eek[tet[0][it]])*(eps-eek[tet[0][it]]));
		 }
         else
            faclow[it]=0.0;
//      .. between e3 and e4 the intersection is triangular
         if (eek[tet[3][it]]>eek[tet[2][it]])
		 {
            eps=0.50*(eek[tet[2][it]]+eek[tet[3][it]]);
            fachigh[it]=SURF(eps,it)/((eps-eek[tet[3][it]])*(eps-eek[tet[3][it]]));
		 }
         else
            fachigh[it]=0.0;
//      .. between e2 and e3 the intersection is a quadrangle
         if (eek[tet[2][it]]>eek[tet[1][it]])
		 {
            eps=0.5*(eek[tet[2][it]]+eek[tet[1][it]]);
            facmedium[it]=(faclow[it]*(eps-eek[tet[0][it]])*(eps-eek[tet[0][it]])
                         -SURF(eps,it))/((eps-eek[tet[1][it]])*(eps-eek[tet[1][it]]));
         }
		 else
            facmedium[it]=0.0;
	  }
//____End of CALCFAC
      return ;
}
//====
	double Band::FASTSURF(double eps,int it)
	{
      double rtnv;
      //int it;
//    : Calculate the area of the intersection of a equi energy plane with
//      a tetrahedron
      if (eps>eek[tet[0][it]]&&eps<eek[tet[3][it]])
	  {
//_______triangular intersection with energy lower than e2
         if (eps<=eek[tet[1][it]])
            rtnv=faclow[it]*(eps-eek[tet[0][it]])*(eps-eek[tet[0][it]]);
//_______quadrangular intersection with energy between e2 and e3
         else if (eps<eek[tet[2][it]])
		 {
			
            if (eek[tet[1][it]]-eek[tet[0][it]]>1e-10)
			   rtnv=faclow[it]*(eps-eek[tet[0][it]])*(eps-eek[tet[0][it]])
                       -facmedium[it]*(eps-eek[tet[1][it]])*(eps-eek[tet[1][it]]);
            else
               rtnv=SURF(eps,it);
		 }	      
//_______triangular intersection with energy higher than e3
         else
            rtnv=fachigh[it]*(eps-eek[tet[3][it]])*(eps-eek[tet[3][it]]);
	  }
      else
         rtnv=0.0;
      if (rtnv<0.0)
	  {/*
		 WRITE (LUTTYO,*) 'FASTSURF: Area<0'
         WRITE (LUOUT ,*) ''
         STOP
		*/
		  cout<<"error FASTSURF: Area<0";exit(0);
	  }
//____End of FASTSURF
      return rtnv;
	}
//====
void Band::GETMAXAOVTET(void)
{
//________ Purpose   : calculate the maximum AOV in tet.
//                     this quantity is area over velocity
//____local variables
      int it,is;
      double  area,epsu,epsm,epsl;
      double  areal,areau;
//____loop over all tetrahedra
      for(it=0;it<nt;it++)
	  {
	  //DO it=1,nt
//_______cal. maximum area for each tetrahedron between the energies e2 and e3
//      ..between the energy of node 2 and 3 the area has a maximum
//        the maximum in the interval between epsl and epsu is found
//        with bisection (it has been assumed,that 20 iterations are 
//        sufficient).
         epsl=eek[tet[1][it]];
         epsu=eek[tet[2][it]];
         areal=FASTSURF(epsl,it);
         areau=FASTSURF(epsu,it);
         for(is=0;is<20;is++)
		 {
		 //DO is=1,20
            epsm=0.50*(epsl+epsu);
            if (areau>areal)
			{
               epsl=epsm;
               areal=FASTSURF(epsl,it);
			}
            else
			{
               epsu=epsm;
               areau=FASTSURF(epsu,it);
			}
		 }
		 area=MAX(areal,areau);
//   ..take maximum
         maxaovtet[it]=area/vgt[it];
	  }
//____End of GETMAXAOVTET
      return;
}

void  Band::READBS(void)
{
	  //READ band structure data
      int  it,ik,iptype;
      int  ntl;
	  //FILE*lutmp;
      ifstream ftp;
      if(masterfl)
	  {
	  //read tetrahedra data (not normalized), the nodes of the tetrahedra
	  //must be sorted in such a way, that their energy is in ascending order
	  //with the node number (e1 < e2 < e3 < e4)!
         if(sifl)
		 {
//xiazl		 ftp.open("bs_si.asc");
			 if(material==0)	ftp.open("c:/montecarlo/input/bs.si.asc");
			 else if(material==1)	ftp.open("c:/montecarlo/input/bs.ge.trs");
			 assert(ftp);
		 }
         else if(gaasfl)
		 { 	
			 ftp.open("bs_gaas.asc");
		 } 	
         //ENDIF
		 
	     //Number of energy bands
         //READ(lutmp) nband
		 for(int i=0;i<NPARTYP;i++)ftp>>nband[i];
         if(nband[PELEC]>NBE)
		 {
			 cout<<"error ETABF: nbe > MNBE";exit(0);
		 }
         if(nband[PHOLE]>NBH)
		 {
			 cout<<"error ETABF: nbh > MNBH";exit(0);
		 }
         if(nband[POXEL]>NBOE)
		 {
			 cout<<"error ETABF: nboe > MNBOE";exit(0);
		 }
         nbt=nband[PELEC]+nband[PHOLE]+nband[POXEL];

		 //set band pointers
         bandof[PELEC]=0;
         bandof[PHOLE]=nband[PELEC];
         bandof[POXEL]=nband[PELEC]+nband[PHOLE];
		 for(iptype=0;iptype<NPARTYP;iptype++)
		 {
         //DO iptype=1, NPARTYP
			for(int iband=bandof[iptype];iband<bandof[iptype]+nband[iptype];iband++)
			{
            //DO iband=bandof[iptype]+1, bandof[iptype]+nband[iptype]
               partyp[iband]=iptype;
            //ENDDO
			}
         //ENDDO
		 }

		 //Number of k-space points
         //READ(lutmp) nk
		 ftp>>nk;
         if(nk>MNK)
		 {

		 }
		 //read k-vector and energy
		 for(ik=0;ik<nk;ik++)ftp>>xkk[ik];
		 for(ik=0;ik<nk;ik++)ftp>>ykk[ik];
		 for(ik=0;ik<nk;ik++)ftp>>zkk[ik];
		 for(ik=0;ik<nk;ik++)ftp>>eek[ik];
		 for(ik=0;ik<nk;ik++)
		 {
         //DO ik=1, nk
            if(eek[ik]<0.0)
			{ 
				cout<<"error ETABF: eek < 0";exit(0);
			}
		 }
         //ENDDO

		 //read number of tetrahedra
         //READ(lutmp) nt
		 ftp>>nt;

         if(nt>MNT)
		 {  
			 cout<<"error ETABF: nt > MNT";exit(0);
		 }
		 //read grid point indices for each of the four corners of a tetrahedron
		 int iitemp;
		 for(it=0;it<nt;it++){ftp>>iitemp;tet[0][it]=iitemp-1;}
		 for(it=0;it<nt;it++){ftp>>iitemp;tet[1][it]=iitemp-1;}
		 for(it=0;it<nt;it++){ftp>>iitemp;tet[2][it]=iitemp-1;}
		 for(it=0;it<nt;it++){ftp>>iitemp;tet[3][it]=iitemp-1;}
		 
		 //read list of neighbouring tetrahedra
		 for(it=0;it<nt;it++){ftp>>iitemp;nlt[0][it]=iitemp-1;}
		 for(it=0;it<nt;it++){ftp>>iitemp;nlt[1][it]=iitemp-1;}
		 for(it=0;it<nt;it++){ftp>>iitemp;nlt[2][it]=iitemp-1;}
		 for(it=0;it<nt;it++){ftp>>iitemp;nlt[3][it]=iitemp-1;}
		 //read band index of tetrahedron
         //READ(lutmp) (ibt[it],it=1,nt)
		 for(it=0;it<nt;it++){ftp>>iitemp;ibt[it]=iitemp-1;}
         ftp.close();
		 //normalize k-space grid point data
		 for(ik=0;ik<nk;ik++)
		 {
            xkk[ik]=xkk[ik]*a0pi;
            ykk[ik]=ykk[ik]*a0pi;
            zkk[ik]=zkk[ik]*a0pi;
            eek[ik]=eek[ik]/eV0;
		 }
		 //check order of the energy of at the corners of each tetrahedron 
		 for(it=0;it<nt;it++)
		 {
            if(!((eek[tet[0][it]]<=eek[tet[1][it]])||
                      (eek[tet[1][it]]<=eek[tet[2][it]])||
                      (eek[tet[2][it]]<=eek[tet[3][it]])))
			{  
				cout<<"error ETABF: energy order wrong";exit(0);
			}
		 }

      }//endif(masterfl)
		 

	   //read inverse mass tensor
      if(massfl)
	  {	
         if(masterfl)
		 {
            if(sifl)
			{
				ftp.open("imt_si.asc");
				assert(ftp);
			}
            else
			{
				if(sifl)
				{
					ftp.open("imt_gaas.asc");
					assert(ftp);
				}
			}
			//read number of tetrahedra
			ftp>>ntl;
			//check if number of tetrahedra is equivalent to the previous 
			//loaded one
            if(nt!=ntl)
			{
				cout<<"error ETABF: nt <> ntl";exit(0);
			}
			for(it=0;it<nt;it++)
			{
				for(int i=0;i<6;i++)ftp>>massinv[i][it];
			}
            ftp.close();
		 }//end if(masterfl)
	  }
         


	  //cal. center of tetrahedron
      for(it=0;it<nt;it++)
	  {	  
	  //DO it=1, nt
         xkct[it]=0.250*(xkk[tet[0][it]]+xkk[tet[1][it]]
			+xkk[tet[2][it]]+xkk[tet[3][it]]);
		 ykct[it]=0.250*(ykk[tet[0][it]]+ykk[tet[1][it]]+ykk[tet[2][it]]
			+ykk[tet[3][it]]);
         zkct[it]=0.250*(zkk[tet[0][it]]+zkk[tet[1][it]]+zkk[tet[2][it]]
			+zkk[tet[3][it]]);
	  //ENDDO
	  }	
	  //cal. velocities within each tetrahedron
      VELOTET();

	  //get vectors normal to the four surfaces of each tetrahedron
	  //(Hessian form)
      HESSTET();

	  //matrix operations (the 48 symmetry transformations of the BZ)
      SETMATSYM();

	  //cal. the five neighbouring wedges for each of the 48 wedges
      SETNEIBSYM();

	  //cal. factors for the determination of area within a tetrahedron
	  //used in FASTSURF
      CALCFAC();

	  //calculate maximum area over velocity for each tetrahedron
      GETMAXAOVTET();

	  //check tetrahedra

//C     collision broadening
//	  int it,ik;
//c	   kloem
	if(material==0)
	{
	  ifstream ftpkloem;
	  ftpkloem.open("c:/montecarlo/input/kloem.txt");
	  for(it=0;it<120;it++)
		  for(ik=0;ik<201;ik++)
			  ftpkloem>>kloem[ik][it];
	  ftpkloem.close();
//c	   kloab
	  ifstream ftpkloab;
	  ftpkloab.open("c:/montecarlo/input/kloab.txt");
	  for(it=0;it<110;it++)
		  for(ik=0;ik<201;ik++)
			  ftpkloab>>kloab[ik][it];
	  ftpkloab.close();
//c	   ktoem
	  ifstream ftpktoem;
	  ftpktoem.open("c:/montecarlo/input/ktoem.txt");
	  for(it=0;it<160;it++)
		  for(ik=0;ik<51;ik++)
			  ftpktoem>>ktoem[ik][it];
	  ftpktoem.close();
//c	   ktoab
	  ifstream ftpktoab;
	  ftpktoab.open("c:/montecarlo/input/ktoab.txt");
	  for(it=0;it<150;it++)
		  for(ik=0;ik<51;ik++)
			  ftpktoab>>ktoab[ik][it];
	  ftpktoab.close();
//c	   klaem
	  ifstream ftpklaem;
	  ftpklaem.open("c:/montecarlo/input/klaem.txt");
	  for(it=0;it<160;it++)
		  for(ik=0;ik<51;ik++)
			  ftpklaem>>klaem[ik][it];
	  ftpklaem.close();
//c	   klaab
	  ifstream ftpklaab;
	  ftpklaab.open("c:/montecarlo/input/klaab.txt");
	  for(it=0;it<150;it++)
		  for(ik=0;ik<51;ik++)
			  ftpkloem>>klaab[ik][it];
	  ftpklaab.close();
	}
	else if(material==1)
	{
	  ifstream ftpkloem;
	  ftpkloem.open("c:/montecarlo/input/kloemfi.txt");
	  for(it=0;it<150;it++)
		  for(ik=0;ik<51;ik++)
			  ftpkloem>>kloem[ik][it];
	  ftpkloem.close();
//c	   kloab
	  ifstream ftpkloab;
	  ftpkloab.open("c:/montecarlo/input/kloabfi.txt");
	  for(it=0;it<150;it++)
		  for(ik=0;ik<51;ik++)
			  ftpkloab>>kloab[ik][it];
	  ftpkloab.close();
//c	   ktoab
//c	   klaem
	  ifstream ftpklaem;
	  ftpklaem.open("c:/montecarlo/input/klaemfi.txt");
	  for(it=0;it<160;it++)
		  for(ik=0;ik<51;ik++)
			  ftpklaem>>klaem[ik][it];
	  ftpklaem.close();
//c	   klaab
	  ifstream ftpklaab;
	  ftpklaab.open("c:/montecarlo/input/klaabfi.txt");
	  for(it=0;it<150;it++)
		  for(ik=0;ik<51;ik++)
			  ftpkloem>>klaab[ik][it];
	  ftpklaab.close();
	  }
	  //End of READBS 
      return;
      //END
}
void Band::SETMATSYM(void)
{
	  //Purpose:   initialize the 48 symmetry transformations of silicon
	  //--------
	  //local variables
      int isym;
	  //matrix operations (the 48 symmetry transformations of the BZ)
      matsym[3][0]=1;
      matsym[4][0]=1;
      matsym[5][0]=1;
      matsym[0][0]=1;
      matsym[1][0]=2;
      matsym[2][0]=3;
      matsym[3][1]=1;
      matsym[4][1]=1;
      matsym[5][1]=-1;
      matsym[0][1]=1;
      matsym[1][1]=2;
      matsym[2][1]=3;
      matsym[3][2]=1;
      matsym[4][2]=-1;
      matsym[5][2]=1;
      matsym[0][2]=1;
      matsym[1][2]=2;
      matsym[2][2]=3;
      matsym[3][3]=1;
      matsym[4][3]=-1;
      matsym[5][3]=-1;
      matsym[0][3]=1;
      matsym[1][3]=2;
      matsym[2][3]=3;
      matsym[3][4]=1;
      matsym[4][4]=1;
      matsym[5][4]=1;
      matsym[0][4]=1;
      matsym[1][4]=3;
      matsym[2][4]=2;
      matsym[3][5]=1;
      matsym[4][5]=1;
      matsym[5][5]=-1;
      matsym[0][5]=1;
      matsym[1][5]=3;
      matsym[2][5]=2;
      matsym[3][6]=1;
      matsym[4][6]=-1;
      matsym[5][6]=1;
      matsym[0][6]=1;
      matsym[1][6]=3;
      matsym[2][6]=2;
      matsym[3][7]=1;
      matsym[4][7]=-1;
      matsym[5][7]=-1;
      matsym[0][7]=1;
      matsym[1][7]=3;
      matsym[2][7]=2;
      matsym[3][8]=-1;
      matsym[4][8]=1;
      matsym[5][8]=1;
      matsym[0][8]=1;
      matsym[1][8]=2;
      matsym[2][8]=3;
      matsym[3][9]=-1;
      matsym[4][9]=1;
      matsym[5][9]=-1;
      matsym[0][9]=1;
      matsym[1][9]=2;
      matsym[2][9]=3;
      matsym[3][10]=-1;
      matsym[4][10]=-1;
      matsym[5][10]=1;
      matsym[0][10]=1;
      matsym[1][10]=2;
      matsym[2][10]=3;
      matsym[3][11]=-1;
      matsym[4][11]=-1;
      matsym[5][11]=-1;
      matsym[0][11]=1;
      matsym[1][11]=2;
      matsym[2][11]=3;
      matsym[3][12]=-1;
      matsym[4][12]=1;
      matsym[5][12]=1;
      matsym[0][12]=1;
      matsym[1][12]=3;
      matsym[2][12]=2;
      matsym[3][13]=-1;
      matsym[4][13]=1;
      matsym[5][13]=-1;
      matsym[0][13]=1;
      matsym[1][13]=3;
      matsym[2][13]=2;
      matsym[3][14]=-1;
      matsym[4][14]=-1;
      matsym[5][14]=1;
      matsym[0][14]=1;
      matsym[1][14]=3;
      matsym[2][14]=2;
      matsym[3][15]=-1;
      matsym[4][15]=-1;
      matsym[5][15]=-1;
      matsym[0][15]= 1;
      matsym[1][15]= 3;
      matsym[2][15]= 2;
      matsym[3][16]=1;
      matsym[4][16]=1;
      matsym[5][16]=1;
      matsym[0][16]=2;
      matsym[1][16]=1;
      matsym[2][16]=3;
      matsym[3][17]=1;
      matsym[4][17]=1;
      matsym[5][17]=-1;
      matsym[0][17]=2;
      matsym[1][17]=1;
      matsym[2][17]=3;
      matsym[3][18]=-1;
      matsym[4][18]= 1;
      matsym[5][18]=-1;
      matsym[0][18]= 2;
      matsym[1][18]= 1;
      matsym[2][18]= 3;
      matsym[3][19]=-1;
      matsym[4][19]= 1;
      matsym[5][19]= 1;
      matsym[0][19]= 2;
      matsym[1][19]= 1;
      matsym[2][19]= 3;
      matsym[3][20]=1;
      matsym[4][20]=1;
      matsym[5][20]=1;
      matsym[0][20]=3;
      matsym[1][20]=1;
      matsym[2][20]=2;
      matsym[3][21]=1;
      matsym[4][21]=1;
      matsym[5][21]=-1;
      matsym[0][21]=3;
      matsym[1][21]=1;
      matsym[2][21]=2;
      matsym[3][22]=-1;
      matsym[4][22]= 1;
      matsym[5][22]=-1;
      matsym[0][22]= 3;
      matsym[1][22]= 1;
      matsym[2][22]= 2;
      matsym[3][23]=-1;
      matsym[4][23]= 1;
      matsym[5][23]= 1;
      matsym[0][23]= 3;
      matsym[1][23]= 1;
      matsym[2][23]= 2;
      matsym[3][24]=1;
      matsym[4][24]=-1;
      matsym[5][24]=1;
      matsym[0][24]=2;
      matsym[1][24]=1;
      matsym[2][24]=3;
      matsym[3][25]=1;
      matsym[4][25]=-1;
      matsym[5][25]=-1;
      matsym[0][25]=2;
      matsym[1][25]=1;
      matsym[2][25]=3;
      matsym[3][26]=-1;
      matsym[4][26]=-1;
      matsym[5][26]=-1;
      matsym[0][26]= 2;
      matsym[1][26]= 1;
      matsym[2][26]= 3;
      matsym[3][27]=-1;
      matsym[4][27]=-1;
      matsym[5][27]= 1;
      matsym[0][27]= 2;
      matsym[1][27]= 1;
      matsym[2][27]= 3;
      matsym[3][28]=1;
      matsym[4][28]=-1;
      matsym[5][28]=1;
      matsym[0][28]=3;
      matsym[1][28]=1;
      matsym[2][28]=2;
      matsym[3][29]=1;
      matsym[4][29]=-1;
      matsym[5][29]=-1;
      matsym[0][29]=3;
      matsym[1][29]=1;
      matsym[2][29]=2;
      matsym[3][30]=-1;
      matsym[4][30]=-1;
      matsym[5][30]=-1;
      matsym[0][30]= 3;
      matsym[1][30]= 1;
      matsym[2][30]= 2;
      matsym[3][31]=-1;
      matsym[4][31]=-1;
      matsym[5][31]= 1;
      matsym[0][31]= 3;
      matsym[1][31]= 1;
      matsym[2][31]= 2;
      matsym[3][32]=1;
      matsym[4][32]=1;
      matsym[5][32]=1;
      matsym[0][32]=2;
      matsym[1][32]=3;
      matsym[2][32]=1;
      matsym[3][33]=1;
      matsym[4][33]=-1;
      matsym[5][33]=1;
      matsym[0][33]=2;
      matsym[1][33]=3;
      matsym[2][33]=1;
      matsym[3][34]=-1;
      matsym[4][34]=-1;
      matsym[5][34]= 1;
      matsym[0][34]= 2;
      matsym[1][34]= 3;
      matsym[2][34]= 1;
      matsym[3][35]=-1;
      matsym[4][35]= 1;
      matsym[5][35]= 1;
      matsym[0][35]= 2;
      matsym[1][35]= 3;
      matsym[2][35]= 1;
      matsym[3][36]=1;
      matsym[4][36]=1;
      matsym[5][36]=1;
      matsym[0][36]=3;
      matsym[1][36]=2;
      matsym[2][36]=1;
      matsym[3][37]=1;
      matsym[4][37]=-1;
      matsym[5][37]=1;
      matsym[0][37]=3;
      matsym[1][37]=2;
      matsym[2][37]=1;
      matsym[3][38]=-1;
      matsym[4][38]=-1;
      matsym[5][38]= 1;
      matsym[0][38]= 3;
      matsym[1][38]= 2;
      matsym[2][38]= 1;
      matsym[3][39]=-1;
      matsym[4][39]= 1;
      matsym[5][39]= 1;
      matsym[0][39]= 3;
      matsym[1][39]= 2;
      matsym[2][39]= 1;
      matsym[3][40]=1;
      matsym[4][40]=1;
      matsym[5][40]=-1;
      matsym[0][40]=2;
      matsym[1][40]=3;
      matsym[2][40]=1;
      matsym[3][41]=1;
      matsym[4][41]=-1;
      matsym[5][41]=-1;
      matsym[0][41]=2;
      matsym[1][41]=3;
      matsym[2][41]=1;
      matsym[3][42]=-1;
      matsym[4][42]=-1;
      matsym[5][42]=-1;
      matsym[0][42]= 2;
      matsym[1][42]= 3;
      matsym[2][42]= 1;
      matsym[3][43]=-1;
      matsym[4][43]= 1;
      matsym[5][43]=-1;
      matsym[0][43]= 2;
      matsym[1][43]= 3;
      matsym[2][43]= 1;
      matsym[3][44]=1;
      matsym[4][44]=1;
      matsym[5][44]=-1;
      matsym[0][44]=3;
      matsym[1][44]=2;
      matsym[2][44]=1;
      matsym[3][45]=1;
      matsym[4][45]=-1;
      matsym[5][45]=-1;
      matsym[0][45]=3;
      matsym[1][45]=2;
      matsym[2][45]=1;
      matsym[3][46]=-1;
      matsym[4][46]=-1;
      matsym[5][46]=-1;
      matsym[0][46]= 3;
      matsym[1][46]= 2;
      matsym[2][46]= 1;
      matsym[3][47]=-1;
      matsym[4][47]= 1;
      matsym[5][47]=-1;
      matsym[0][47]= 3;
      matsym[1][47]= 2;
      matsym[2][47]= 1;

	  for(isym=0;isym<48;isym++)
	  {	
      //DO isym=1,48 
         indmat[matsym[3][isym]+1][matsym[4][isym]+1][matsym[5][isym]+1]
			     [matsym[0][isym]-1][matsym[1][isym]-1][matsym[2][isym]-1] 
              =isym;
		 mapsym[0][isym]=isym+1;
		 mapsym[1][isym]=isym+1;
		 mapsym[2][isym]=isym+1;
		 mapsym[3][isym]=isym+1;
		 mapsym[4][isym]=isym+1;
		 mapsym[5][isym]=isym+1;
		 mapsym[6][isym]=isym+1;
		 mapsym[7][isym]=isym+1;
	  }	

	  mapsym[0][0]=43;
      mapsym[0][42]=1;
	  mapsym[0][4]=44;
      mapsym[0][43]=5;
	  mapsym[0][16]=45;
      mapsym[0][44]=17;
	  mapsym[0][20]=46;
      mapsym[0][45]=21;
	  mapsym[0][32]=47;
      mapsym[0][46]=33;
	  mapsym[0][36]=48;
      mapsym[0][47]=37;

	  mapsym[1][1]=43;
      mapsym[1][42]=2;
	  mapsym[1][5]=44;
      mapsym[1][43]=6;
	  mapsym[1][17]=45;
      mapsym[1][44]=18;
	  mapsym[1][21]=46;
      mapsym[1][45]=22;
	  mapsym[1][40]=47;
      mapsym[1][46]=41;
      mapsym[1][47]=45;

	  mapsym[2][3]=43;
      mapsym[2][42]=4;
	  mapsym[2][7]=44;
      mapsym[2][43]=8;
	  mapsym[2][25]=45;
      mapsym[2][44]=26;
	  mapsym[2][29]=46;
      mapsym[2][45]=30;
	  mapsym[2][41]=47;
      mapsym[2][46]=42;
      mapsym[2][47]=46;

	  mapsym[3][2]=43;
      mapsym[3][42]=3;
	  mapsym[3][6]=44;
      mapsym[3][43]=7;
	  mapsym[3][24]=45;
      mapsym[3][44]=25;
	  mapsym[3][28]=46;
      mapsym[3][45]=29;
	  mapsym[3][33]=47;
      mapsym[3][46]=34;
	  mapsym[3][37]=48;
      mapsym[3][47]=38;

//c Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	  mapsym[4][8]=43;
      mapsym[4][42]=9;
	  mapsym[4][12]=44;
      mapsym[4][43]=13;
	  mapsym[4][19]=45;
      mapsym[4][44]=20;
	  mapsym[4][23]=46;
      mapsym[4][45]=24;
	  mapsym[4][35]=47;
      mapsym[4][46]=36;
	  mapsym[4][39]=48;
      mapsym[4][47]=40;

	  mapsym[5][9]=43;
      mapsym[5][42]=10;
	  mapsym[5][13]=44;
      mapsym[5][43]=14;
	  mapsym[5][18]=45;
      mapsym[5][44]=19;
	  mapsym[5][22]=46;
      mapsym[5][45]=23;
      mapsym[5][46]=44;
	  mapsym[5][47]=48;

	  mapsym[6][11]=43;
      mapsym[6][42]=12;
	  mapsym[6][15]=44;
      mapsym[6][43]=16;
	  mapsym[6][26]=45;
      mapsym[6][44]=27;
	  mapsym[6][30]=46;
      mapsym[6][45]=31;
	  mapsym[6][46]=43;
	  mapsym[6][47]=47;

	  mapsym[7][10]=43;
      mapsym[7][42]=11;
	  mapsym[7][14]=44;
      mapsym[7][43]=15;
	  mapsym[7][27]=45;
      mapsym[7][44]=28;
	  mapsym[7][31]=46;
      mapsym[7][45]=32;
	  mapsym[7][34]=47;
	  mapsym[7][46]=35;
	  mapsym[7][38]=48;
	  mapsym[7][47]=39;
	  //end of SETMATSYM
      return;
      //END
}

//rewrite by dugang DEC.14 
void Band::BUILDPHSCATT(void)//(cscat)
{
//____calculate phonon scattering rates

     //CHARACTER*(CPCVL) cscat
	  int  itab,iscat,it;
      int  itabl,itabu,itabfp;
      int  iinital,ifinal;
      double k1em,k2em,k3em,k4em,k5em,k6em ;
      double k1ab,k2ab,k3ab,k4ab,k5ab,k6ab, 
                       kinth,kemh,kabh,besh,
                       kintoe,kemoe,kaboe,besoe,
                      bes[6],dosd1,dosd2,kint,
                      eeafter,sii,
                      kopemlow,kopablow,kopemhigh,kopabhigh,
                      kaptemlow,kaptablow,kaptemhigh,kaptabhigh,
                      kaplemlow,kaplablow,kaplemhigh,kaplabhigh;
     //double CALDOS,CALDOSSUM,HIIRATE,
	  double dosfac;
	  double ee;
	  int ipt;

      if(jacophfl)
	  {
//_______Bose-Einstein for optical phonons
		if(material==0)
		{
//_______Bose-Einstein for optical phonons
         bes[0] =1.0/(exp(temptag)-1.0); 
         bes[1] =1.0/(exp(templag)-1.0);
         bes[2] =1.0/(exp(templog)-1.0); 
         bes[3] =1.0/(exp(temptaf)-1.0); 
         bes[4] =1.0/(exp(templaf)-1.0); 
         bes[5] =1.0/(exp(temptof)-1.0); 


//_______optical phonon constants for emission und absorbtion(electrons)
//     (phonons+II)
         scpre=14; 
         kint=2.0*PI*dfelast*dfelast/sirho/6e0/(siul*siul);

         k1em=(dftag)*(dftag)*PI/temptag/sirho/6e0;
         k2em=(dflag)*(dflag)*PI/templag/sirho/6e0;
         k3em=(dflog)*(dflog)*PI/templog/sirho/6e0;
         k4em=(dftaf)*(dftaf)*PI/temptaf/sirho*0.666e0;
         k5em=(dflaf)*(dflaf)*PI/templaf/sirho*0.666e0;
         k6em=(dftof)*(dftof)*PI/temptof/sirho*0.666e0;

         k1ab=bes[0]*k1em;
         k2ab=bes[1]*k2em;
         k3ab=bes[2]*k3em;
         k4ab=bes[3]*k4em;
         k5ab=bes[4]*k5em;
         k6ab=bes[5]*k6em;

         k1em=(bes[0]+1.0)*k1em;
         k2em=(bes[1]+1.0)*k2em;
         k3em=(bes[2]+1.0)*k3em;
         k4em=(bes[3]+1.0)*k4em;
         k5em=(bes[4]+1.0)*k5em;
         k6em=(bes[5]+1.0)*k6em;

//_______calculate scattering-rates(electrons)
         ipt=PELEC;
		 for(ifinal=0;ifinal<nband[ipt];ifinal++)
		 {
			 //DO ifinal  =1,nband[ipt]
		  for(itab=0;itab<=MTAB;itab++)
		  {//DO itab=0,MTAB
            ee=energy[itab];

//__________intra-valley elastic Phonon
            dose[0][ifinal][itab]=CALDOS(ee,ifinal+bandof[ipt]);
   
//__________transversal acoustic g-Phonon(Absorption)
            dose[1][ifinal][itab]=CALDOS(ee+temptag,ifinal+bandof[ipt]);

//__________transversal acoustic g-Phonon(Emission)
            dose[2][ifinal][itab]=CALDOS(ee-temptag,ifinal+bandof[ipt]);

//__________longitudinal acoustic g-Phonon(Absorption)
            dose[3][ifinal][itab]=CALDOS(ee+templag,ifinal+bandof[ipt]);

//__________longitudinal acoustic g-Phonon(Emission)
            dose[4][ifinal][itab]=CALDOS(ee-templag,ifinal+bandof[ipt]);

//__________longitudinal optical g-Phonon(Absorption)
            dose[5][ifinal][itab]=CALDOS(ee+templog,ifinal+bandof[ipt]);

//__________longitudinal optical g-Phonon(Emission)
            dose[6][ifinal][itab]=CALDOS(ee-templog,ifinal+bandof[ipt]);

//__________transversal acoustic f-Phonon(Absorption)
            dose[7][ifinal][itab]=CALDOS(ee+temptaf,ifinal+bandof[ipt]);

//__________transversal acoustic f-Phonon(Emission)
            dose[8][ifinal][itab]=CALDOS(ee-temptaf,ifinal+bandof[ipt]);

//__________longitudinal acoustic f-Phonon(Absorption)
            dose[9][ifinal][itab]=CALDOS(ee+templaf,ifinal+bandof[ipt]);

//__________longitudinal acoustic f-Phonon(Emission)
            dose[10][ifinal][itab]=CALDOS(ee-templaf,ifinal+bandof[ipt]);

//__________transversal optical f-Phonon(Absorption)
            dose[11][ifinal][itab]=CALDOS(ee+temptof,ifinal+bandof[ipt]);

//__________transversal optical f-Phonon(Emission)
            dose[12][ifinal][itab]=CALDOS(ee-temptof,ifinal+bandof[ipt]);

          }//ENDDO
         }//ENDDO

         for(iinital=0;iinital<nband[ipt];iinital++)
		 {//DO iinital =1,nband[ipt]
           for(ifinal=0;ifinal<nband[ipt];ifinal++)
		   {//DO ifinal  =1,nband[ipt]

            if(iinital==0&&ifinal==0)///be caareful not 1 but 0
               dosfac=1.0;
            else 
               dosfac=ephb;
            
//__________intra-valley elastic Phonon
            scatte[0][ifinal][iinital]=kint*dosfac;

//__________transversal acoustic g-Phonon(Absorption)
            scatte[1][ifinal][iinital]=k1ab*dosfac;

//__________transversal acoustic g-Phonon(Emission)
            scatte[2][ifinal][iinital]=k1em*dosfac;

//__________longitudinal acoustic g-Phonon(Absorption)
            scatte[3][ifinal][iinital]=k2ab*dosfac;

//__________longitudinal acoustic g-Phonon(Emission)
            scatte[4][ifinal][iinital]=k2em*dosfac;

//__________longitudinal optical g-Phonon(Absorption)
            scatte[5][ifinal][iinital]=k3ab*dosfac;

//__________longitudinal optical g-Phonon(Emission)
            scatte[6][ifinal][iinital]=k3em*dosfac;

//__________transversal acoustic f-Phonon(Absorption)
            scatte[7][ifinal][iinital]=k4ab*dosfac;

//__________transversal acoustic f-Phonon(Emission)
            scatte[8][ifinal][iinital]=k4em*dosfac;

//__________longitudinal acoustic f-Phonon(Absorption)
            scatte[9][ifinal][iinital]=k5ab*dosfac;

//__________longitudinal acoustic f-Phonon(Emission)
            scatte[10][ifinal][iinital]=k5em*dosfac;

//__________transversal optical f-Phonon(Absorption)
            scatte[11][ifinal][iinital]=k6ab*dosfac;

//__________transversal optical f-Phonon(Emission)
            scatte[12][ifinal][iinital]=k6em*dosfac;

		   }//ENDDO
         }//ENDDO
		}
		else if(material==1)
		{
         bes[0] =1.0/(exp(temptag)-1.0); 
         bes[1] =1.0/(exp(templag)-1.0);
         bes[2] =1.0/(exp(templog)-1.0); 
//       bes[3] =1.0/(exp(temptaf)-1.0); 
//         bes[4] =1.0/(exp(templaf)-1.0); 
//         bes[5] =1.0/(exp(temptof)-1.0); 


//_______optical phonon constants for emission und absorbtion(electrons)
//     (phonons+II) 
		 //200405
         scpre=8; 
//         kint=2.0*PI*dfelast*dfelast/sirho/6e0/(siul*siul);
		 kint=2.0*PI*dfelast*dfelast/sirho*0.125/(siul*siul);
//         k1em=(dftag)*(dftag)*PI/temptag/sirho/6e0;
		 k1em=(dftag)*(dftag)*PI/temptag/sirho*0.875e0;
//         k2em=(dflag)*(dflag)*PI/templag/sirho/6e0;
		 k2em=(dflag)*(dflag)*PI/templag/sirho*0.875e0;
//         k3em=(dflog)*(dflog)*PI/templog/sirho/6e0;
		 k3em=(dflog)*(dflog)*PI/templog/sirho*0.875e0;
//         k4em=(dftaf)*(dftaf)*PI/temptaf/sirho*0.666e0;
//         k5em=(dflaf)*(dflaf)*PI/templaf/sirho*0.666e0;
//         k6em=(dftof)*(dftof)*PI/temptof/sirho*0.666e0;

         k1ab=bes[0]*k1em;
         k2ab=bes[1]*k2em;
         k3ab=bes[2]*k3em;
//         k4ab=bes[3]*k4em;
//         k5ab=bes[4]*k5em;
//         k6ab=bes[5]*k6em;

         k1em=(bes[0]+1.0)*k1em;
         k2em=(bes[1]+1.0)*k2em;
         k3em=(bes[2]+1.0)*k3em;
//         k4em=(bes[3]+1.0)*k4em;
//         k5em=(bes[4]+1.0)*k5em;
//         k6em=(bes[5]+1.0)*k6em;

//_______calculate scattering-rates(electrons)
         ipt=PELEC;
		 for(ifinal=0;ifinal<nband[ipt];ifinal++)
		 {
			 //DO ifinal  =1,nband[ipt]
		  for(itab=0;itab<=MTAB;itab++)
		  {//DO itab=0,MTAB
            ee=energy[itab];

//__________intra-valley elastic Phonon
            dose[0][ifinal][itab]=CALDOS(ee,ifinal+bandof[ipt]);
   
//__________transversal acoustic g-Phonon(Absorption)
            dose[1][ifinal][itab]=CALDOS(ee+temptag,ifinal+bandof[ipt]);

//__________transversal acoustic g-Phonon(Emission)
            dose[2][ifinal][itab]=CALDOS(ee-temptag,ifinal+bandof[ipt]);

//__________longitudinal acoustic g-Phonon(Absorption)
            dose[3][ifinal][itab]=CALDOS(ee+templag,ifinal+bandof[ipt]);

//__________longitudinal acoustic g-Phonon(Emission)
            dose[4][ifinal][itab]=CALDOS(ee-templag,ifinal+bandof[ipt]);

//__________longitudinal optical g-Phonon(Absorption)
            dose[5][ifinal][itab]=CALDOS(ee+templog,ifinal+bandof[ipt]);

//__________longitudinal optical g-Phonon(Emission)
            dose[6][ifinal][itab]=CALDOS(ee-templog,ifinal+bandof[ipt]);

//__________transversal acoustic f-Phonon(Absorption)
//            dose[7][ifinal][itab]=CALDOS(ee+temptaf,ifinal+bandof[ipt]);

//__________transversal acoustic f-Phonon(Emission)
//           dose[8][ifinal][itab]=CALDOS(ee-temptaf,ifinal+bandof[ipt]);

//__________longitudinal acoustic f-Phonon(Absorption)
//            dose[9][ifinal][itab]=CALDOS(ee+templaf,ifinal+bandof[ipt]);

//__________longitudinal acoustic f-Phonon(Emission)
//            dose[10][ifinal][itab]=CALDOS(ee-templaf,ifinal+bandof[ipt]);

//__________transversal optical f-Phonon(Absorption)
//            dose[11][ifinal][itab]=CALDOS(ee+temptof,ifinal+bandof[ipt]);

//__________transversal optical f-Phonon(Emission)
//            dose[12][ifinal][itab]=CALDOS(ee-temptof,ifinal+bandof[ipt]);
          }//ENDDO
         }//ENDDO

         for(iinital=0;iinital<nband[ipt];iinital++)
		 {//DO iinital =1,nband[ipt]
           for(ifinal=0;ifinal<nband[ipt];ifinal++)
		   {//DO ifinal  =1,nband[ipt]

            if(iinital==0&&ifinal==0)///be caareful not 1 but 0
               dosfac=1.0;
            else 
               dosfac=ephb;           
//__________intra-valley elastic Phonon
            scatte[0][ifinal][iinital]=kint*dosfac;

//__________transversal acoustic g-Phonon(Absorption)
            scatte[1][ifinal][iinital]=k1ab*dosfac;

//__________transversal acoustic g-Phonon(Emission)
            scatte[2][ifinal][iinital]=k1em*dosfac;

//__________longitudinal acoustic g-Phonon(Absorption)
            scatte[3][ifinal][iinital]=k2ab*dosfac;

//__________longitudinal acoustic g-Phonon(Emission)
            scatte[4][ifinal][iinital]=k2em*dosfac;

//__________longitudinal optical g-Phonon(Absorption)
            scatte[5][ifinal][iinital]=k3ab*dosfac;

//__________longitudinal optical g-Phonon(Emission)
            scatte[6][ifinal][iinital]=k3em*dosfac;

//__________transversal acoustic f-Phonon(Absorption)
//            scatte[7][ifinal][iinital]=k4ab*dosfac;

//__________transversal acoustic f-Phonon(Emission)
//            scatte[8][ifinal][iinital]=k4em*dosfac;

//__________longitudinal acoustic f-Phonon(Absorption)
//            scatte[9][ifinal][iinital]=k5ab*dosfac;

//__________longitudinal acoustic f-Phonon(Emission)
//            scatte[10][ifinal][iinital]=k5em*dosfac;

//__________transversal optical f-Phonon(Absorption)
//            scatte[11][ifinal][iinital]=k6ab*dosfac;

//__________transversal optical f-Phonon(Emission)
//            scatte[12][ifinal][iinital]=k6em*dosfac;
		   }//ENDDO
         }//ENDDO
		 }
//_______hole constants(phonons+II)
         scprh=4; 
         besh =1.0/(exp(temphop)-1.0); 
         kinth=2e0*PI*dfhelast*dfhelast/sirho/2.0/(siul*siul);
         kemh=dfhop*dfhop*PI/temphop/sirho;
         kabh=besh*kemh;
         kemh=(besh+1.0)*kemh;

//_______calculate scattering-rates(holes)
         ipt=PHOLE;
         for(ifinal=0;ifinal<nband[ipt];ifinal++)
		 {//DO ifinal  =1,nband[ipt]
           for(itab=0;itab<=MTAB;itab++)
		   {//DO itab=0,MTAB
            ee=energy[itab];

//__________intra-valley elastic Phonon
            dosh[0][ifinal][itab]=CALDOS(ee,ifinal+bandof[ipt]);

//__________(Absorption)
            dosh[1][ifinal][itab]=CALDOS(ee+temphop,ifinal+bandof[ipt]);

//__________(Emission)
            dosh[2][ifinal][itab]=CALDOS(ee-temphop,ifinal+bandof[ipt]);

		   }//ENDDO
         }//ENDDO

         for(iinital=0;iinital<nband[ipt];iinital++)
		 {//DO iinital =1,nband[ipt]
           for(ifinal=0;ifinal<nband[ipt];ifinal++)
		   {//DO ifinal  =1,nband[ipt]

//__________intra-valley elastic Phonon
            scatth[0][ifinal][iinital]=kinth;

//__________(Absorption)
            scatth[1][ifinal][iinital]=kabh;

//__________(Emission)
            scatth[2][ifinal][iinital]=kemh;

		   }//ENDDO
         }//ENDDO
      }//if(jcophfl)
      else if(fiscphfl)
	  {
//_______Bose-Einstein for optical phonons
         bes[0] =1.0/(exp(efopee)-1.0); 

//_______seven processes(phonons+II)
         scpre=7; 

//_______optical phonon constants for emission und absorbtion(electrons)
         kopemlow=efoplow*efoplow*PI/efopee/sirho; 
         kopemhigh=efophigh*efophigh*PI/efopee/sirho; 

         kopablow =bes[0]*kopemlow;
         kopabhigh=bes[0]*kopemhigh;
         kopemlow =(bes[0]+1.0)*kopemlow;
         kopemhigh=(bes[0]+1.0)*kopemhigh;

//_______acoustic phonon constants for emission und absorbtion(electrons)
//      including 1 TA and 2 LA modes(LA and TA modes are exchanged)

         kaptablow =35e0*efaplow*efaplow*PI/(efapeet*efapeet*sirho*sia0*sia0);
         kaptabhigh=35e0*(efaphigh*efaphigh)*PI/(efapeet*efapeet)/sirho/(sia0*sia0);

         kaptemlow =(efaplow*efaplow)*(PI*PI*PI)/efapeet/sirho/(sia0*sia0)
                 *(1.0/(exp(efapeet)-1.0)+1.0)*6e0;
         kaptemhigh=(efaphigh*efaphigh)*(PI*PI*PI)/efapeet/sirho/(sia0*sia0)
                 *(1.0/(exp(efapeet)-1.0)+1.0)*6e0;

         kaplablow =35e0*(efaplow*efaplow)*PI/(efapeel*efapeel)/sirho/(sia0*sia0)
                  *2e0;
         kaplabhigh=35e0*(efaphigh*efaphigh)*PI/(efapeel*efapeel)/sirho/(sia0*sia0)
                  *2e0;

         kaplemlow =(efaplow*efaplow)*(PI*PI*PI)/efapeel/sirho/(sia0*sia0)
                *(1.0/(exp(efapeel)-1.0)+1.0)*6e0
                 *2e0;
         kaplemhigh=(efaphigh*efaphigh)*(PI*PI*PI)/efapeel/sirho/(sia0*sia0)
                 *(1.0/(exp(efapeel)-1.0)+1.0)*6e0
                  *2e0;

//_______calculate scattering-rates(electrons)
         ipt=PELEC;
         for(ifinal=0;ifinal<nband[ipt];ifinal++)
		 {//DO ifinal  =1,nband[ipt]
           for(itab=0;itab<=MTAB;itab++)
		   {//DO itab=0,MTAB
            ee=energy[itab];
            itabfp=(int)(ee*dlistfp);//+ 1 be careful!!!!!
            if(itabfp>MWLEFP) itabfp=MWLEFP;

//__________optical Phonon(Absorption)
            dose[0][ifinal][itab]=CALDOS(ee+efopee,ifinal+bandof[ipt]);

//__________optical Phonon(Emission)
            dose[1][ifinal][itab]=CALDOS(ee-efopee,ifinal+bandof[ipt]);

//__________transversal acoustic  Phonon(Absorption)
            dose[2][ifinal][itab]=dosfpab[ifinal+bandof[ipt]][itabfp];

//__________transversal acoustic  Phonon(Emission)
            dose[3][ifinal][itab]=dosfpem[ifinal+bandof[ipt]][itabfp];

//__________longitudinal acoustic  Phonon(Absorption)
            dose[4][ifinal][itab]=dosfpab[ifinal+bandof[ipt]][itabfp];

//__________longitudinal acoustic  Phonon(Emission)
            dose[5][ifinal][itab]=dosfpem[ifinal+bandof[ipt]][itabfp];

		   }//ENDDO
         }//ENDDO

         for(iinital=0;iinital<nband[ipt];iinital++)
		 {//DO iinital =1,nband[ipt]
           for(ifinal=0;ifinal<nband[ipt];ifinal++)
		   {//DO ifinal  =1,nband[ipt]

            if(iinital==0&&ifinal==0)///be caredul
			{
//_____________optical Phonon(Absorption)
               scatte[0][ifinal][iinital]=kopablow;

//_____________optical Phonon(Emission)
               scatte[1][ifinal][iinital]=kopemlow;

//_____________transversal acoustic Phonon(Absorption)
               scatte[2][ifinal][iinital]=kaptablow;

//_____________transversal acoustic Phonon(Emission)
               scatte[3][ifinal][iinital]=kaptemlow;

//_____________longitudinal acoustic Phonon(Absorption)
               scatte[4][ifinal][iinital]=kaplablow;

//_____________longitudinal acoustic Phonon(Emission)
               scatte[5][ifinal][iinital]=kaplemlow;

            }
			else
			{ 

//_____________transversal optical f-Phonon(Absorption)
               scatte[0][ifinal][iinital]=kopabhigh;

//_____________transversal optical f-Phonon(Emission)
               scatte[1][ifinal][iinital]=kopemhigh;

//_____________transversal acoustic Phonon(Absorption)
               scatte[2][ifinal][iinital]=kaptabhigh;

//_____________transversal acoustic Phonon(Emission)
               scatte[3][ifinal][iinital]=kaptemhigh;

//_____________longitudinal acoustic Phonon(Absorption)
               scatte[4][ifinal][iinital]=kaplabhigh;

//_____________longitudinal acoustic Phonon(Emission)
               scatte[5][ifinal][iinital]=kaplemhigh;
            }//endif
		   }//ENDDO
         }//ENDDO

//_______hole constants(phonons+II)
         scprh=4; 
         besh =1.0/(exp(temphop)-1.0); 
         kinth=2e0*PI*(dfhelast*dfhelast)/sirho/2e0/(siul*siul);
         kemh=(dfhop*dfhop)*PI/temphop/sirho;
         kabh=besh*kemh;
         kemh=(besh+1.0)*kemh;
//_______calculate scattering-rates(holes)
         ipt=PHOLE;
         for(ifinal=0;ifinal<nband[ipt];ifinal++)
		 {//DO ifinal  =1,nband[ipt]
           for(itab=0;itab<=MTAB;itab++)
		   {//DO itab=0,MTAB
            ee=energy[itab];

//__________intra-valley elastic Phonon
            dosh[0][ifinal][itab]=CALDOS(ee,ifinal+bandof[ipt]);

//__________(Absorption)
            dosh[1][ifinal][itab]=CALDOS(ee+temphop,ifinal+bandof[ipt]);

//__________(Emission)
            dosh[2][ifinal][itab]=CALDOS(ee-temphop,ifinal+bandof[ipt]);

		   }//ENDDO
         }//ENDDO

         for(iinital=0;iinital<nband[ipt];iinital++)
		 {//DO iinital =1,nband[ipt]
           for(ifinal=0;ifinal<nband[ipt];ifinal++)
		   {//DO ifinal  =1,nband[ipt]

//__________intra-valley elastic Phonon
            scatth[0][ifinal][iinital]=kinth;

//__________(Absorption)
            scatth[1][ifinal][iinital]=kabh;

//__________(Emission)
            scatth[2][ifinal][iinital]=kemh;
		   }//ENDDO
         }//ENDDO
      }
	  else
	  {
		  cout<<"error BUILDPHSCAT: Specify phonon system";exit(0);
      }//endif

//____II-Rate
//____calculate scattering-rates(electrons)
      ipt=PELEC;
      for(iinital=0;iinital<nband[ipt];iinital++)
	  {//DO iinital =1,nband[ipt]
        for(ifinal=0;ifinal<nband[ipt];ifinal++)
		{//DO ifinal  =1,nband[ipt]

         for(itab=0;itab<=MTAB;itab++)
		 {//DO itab=0,MTAB
            ee=energy[itab];
            if(iifl)
			   sii=EIIRATE(ee);
            else 
               sii=0.0;
//     ..calculate final energy after scattering
//     .. formula by Taniguchi
            eeafter=(ee-sieg)/3.0;
            if(eeafter>0.0)
			{
               dosd1=CALDOS(eeafter,ifinal+bandof[ipt]);
               dosd2=CALDOSSUM(eeafter,ipt);
               if(dosd2==0.0)
			   {
                  dosd1=0.0;
                  dosd2=1.0;
               }//endif
            }
			else
			{ 
               dosd1=0.0;
               dosd2=1.0;
            }//endif
            scattiie[ifinal][iinital][itab]=dosd1/dosd2*sii;
         }//ENDDO
		}//ENDDO
      }//ENDDO

//____II-Rate
//____calculate scattering-rates(holes)
      ipt=PHOLE;

      for(iinital=0;iinital<nband[ipt];iinital++)
	  {//DO iinital =1,nband[ipt]
        for(ifinal=0;ifinal<nband[ipt];ifinal++)
		{//DO ifinal  =1,nband[ipt]
         for(itab=0;itab<=MTAB;itab++)
		 {//DO itab=0,MTAB
            ee=energy[itab];
            if(iifl)
               sii=HIIRATE(ee);
            else 
               sii=0.0;
//     .. formula by Taniguchi
            eeafter=(ee-sieg)/3.0;
            if(eeafter>0.0)
			{
               dosd1=CALDOS(eeafter,ifinal+bandof[ipt]);
               dosd2=CALDOSSUM(eeafter,ipt);
               if(dosd2==0e0)
			   {
                  dosd1=0.0;
                  dosd2=1.0;
               }//endif
            }
			else
			{ 
               dosd1=0.0;
               dosd2=1.0;
            }//endif
            scattiih[ifinal][iinital][itab]=dosd1/dosd2*sii;
         }//ENDDO

      }//ENDDO
      }//ENDDO

//____oxide electron constants
      scproe=3;
      besoe =1.0/(exp(tempoeop)-1.0);
      kintoe=2.0*PI*(dfoeelast*dfoeelast )/sirho/2e0/(siul*siul);
      kemoe=(dfoeop*dfoeop)*1.5e0*PI/tempoeop/sirho;
      kaboe=besoe*kemoe;
      kemoe=(besoe+1.0)*kemoe;

//____calculate scattering-rates(oxide electrons)
      ipt=POXEL;
      for(ifinal=0;ifinal<nband[ipt];ifinal++)
	  {//DO ifinal  =1,nband[ipt]
        for(itab=0;itab<=MTAB;itab++)
		{//DO itab=0,MTAB
         ee=energy[itab];

//_______intra-valley elastic Phonon
         dosoe[0][ifinal][itab]=CALDOS(ee,ifinal+bandof[ipt]);

//_______(Absorption)
         dosoe[1][ifinal][itab]=CALDOS(ee+tempoeop,ifinal+bandof[ipt]);

//_______(Emission)
         dosoe[2][ifinal][itab]=CALDOS(ee-tempoeop,ifinal+bandof[ipt]);

		}//ENDDO
      }//ENDDO

      for(iinital=0;iinital<nband[ipt];iinital++)
	  {//DO iinital =1,nband[ipt]
        for(ifinal=0;ifinal<nband[ipt];ifinal++)
		{//DO ifinal  =1,nband[ipt]

//_______intra-valley elastic Phonon
         scattoe[0][ifinal][iinital]=kintoe;

//_______(Absorption)
         scattoe[1][ifinal][iinital]=kaboe;

//_______(Emission)
         scattoe[2][ifinal][iinital]=kemoe;
		}//ENDDO
      }//ENDDO

//____Sum of scattering rates
      for(iinital=0;iinital<MNB;iinital++)
	  {//DO iinital =1,MNB
        for(itab=0;itab<=MTAB;itab++)
		{//DO itab=0,MTAB
         sumscatt[itab][iinital]=0.0;
		}//ENDDO
      }//ENDDO

//  ..electrons
      ipt=PELEC;
      for(itab=0;itab<=MTAB;itab++)
	  {//DO itab=0,MTAB
       for(iinital=0;iinital<nband[ipt];iinital++)
	   {//DO iinital =1,nband[ipt]
        for(ifinal=0;ifinal<nband[ipt];ifinal++)
		{//DO ifinal  =1,nband[ipt]
         for( iscat=0;iscat<scpre-1;iscat++)
		 {//DO iscat=1,scpre-1
		   sumscatt[itab][iinital+bandof[ipt]]=
				 sumscatt[itab][iinital+bandof[ipt]] 
           +scatte[iscat][ifinal][iinital] 
            *dose[iscat][ifinal][itab];
		 }//ENDDO
         sumscatt[itab][iinital+bandof[ipt]]=
            sumscatt[itab][iinital+bandof[ipt]] 
           +scattiie[ifinal][iinital][itab]; 
		}//ENDDO
       }//ENDDO
      }//ENDDO

//  ..holes
      ipt=PHOLE;
      for(itab=0;itab<=MTAB;itab++)
	  {//DO itab=0,MTAB
       for(iinital=0;iinital<nband[ipt];iinital++)
	   {//DO iinital =1,nband[ipt] 
        for(ifinal=0;ifinal<nband[ipt];ifinal++)
		{//DO ifinal  =1,nband[ipt]
         for(iscat=0;iscat<scprh-1;iscat++)
		 {//DO iscat=1,scprh-1
           sumscatt[itab][iinital+bandof[ipt]]=
				sumscatt[itab][iinital+bandof[ipt]] 
	           +scatth[iscat][ifinal][iinital]
		        *dosh[iscat][ifinal][itab];
		 }//ENDDO
         sumscatt[itab][iinital+bandof[ipt]]=
            sumscatt[itab][iinital+bandof[ipt]] 
           +scattiih[ifinal][iinital][itab];
		}//ENDDO
       }//ENDDO
      }//ENDDO

//  ..oxide electrons
      ipt=POXEL;
      for(itab=0;itab<=MTAB;itab++)
	  {//DO itab=0,MTAB
       for(iinital=0;iinital<nband[ipt];iinital++)
	   {//DO iinital =1,nband[ipt]
        for(ifinal=0;ifinal<nband[ipt];ifinal++)
		{//DO ifinal  =1,nband[ipt]
         for(iscat=0;iscat<scproe-1;iscat++)
		 {//DO iscat=1,scproe
           sumscatt[itab][iinital+bandof[ipt]]=
				 sumscatt[itab][iinital+bandof[ipt]] 
                 +scattoe[iscat][ifinal][iinital]
                *dosoe[iscat][ifinal][itab];
		 }//ENDDO
		}//ENDDO
       }//ENDDO
      }//ENDDO

//____cal. maximum scattering rate in each tetrahedron
      for(it=0;it<nt;it++)
	  {//DO it=1,nt
         ipt=partyp[ibt[it]];
         itabl=MAX(0,(int)(eek[tet[0][it]]/dtable));
         itabu=MIN(MTAB,(int)(eek[tet[3][it]]/dtable)+1);
         gamtet[it]=0.0;
         for(itab=itabl;itab<=itabu;itab++)
		 {//DO itab=itabl,itabu
            gamtet[it]=MAX(gamtet[it],sumscatt[itab][ibt[it]]);
         }//ENDDO
         if(gamtet[it]==0.0) gamtet[it]=sumscatt[MTAB][ibt[it]];
      }//ENDDO
//____End of BUILDPHSCATT
      return;
}

double Band::EIIRATE(double ee)
{
//____Purpose : cal. electron II-rate
      double sii,eel,rtnv;
	if(material==0)
	{
//____account for temperature dependent band gap
      eel=ee-sieg+1.1241e0/eV0;
      sii=0.0;

      if(kaneiifl)
	  {
//__Kane
		if(eel*eV0>1.1241e0)
		{
		  sii=7.216e11/scrt0*pow((eel*eV0-1.12410),3.540);
		}//endif
      }//endif

      if(thomiifl)
//__Thoma
		if((eel*eV0>1.128e0)&&(eel*eV0<1.75e0))
		{
         sii=1.25e12/scrt0*pow((eel*eV0-1.128e0),3);
		}
		else if(eel*eV0>1.75e0)
		{
         sii=9.49e12/scrt0*(eel*eV0-1.572e0)*(eel*eV0-1.572e0);
		}//endif
      
      if(sanoiifl)
//__sano
		if(eel*eV0>1.242e0)
		{
		  sii=4.25e11/scrt0*pow((eel*eV0-1.242e0),4.188);
		}//endif
      
      if(fisciifl)
	  {
//___Fischetti
	      if(eel*eV0>1.2e0) 
		  {
		     sii=6.25e10/scrt0*pow(((eel*eV0-1.2e0)/1.20),2);
		  }//endif
          if(eel*eV0>1.8e0)
		  {
			 sii=3e12/scrt0*pow(((eel*eV0-1.80)/1.80),2)+sii;
		  }//endif
          if(eel*eV0>3.45e0)
		  {
			 sii=6.8e14/scrt0*pow(((eel*eV0-3.45e0)/3.45e0),2)+sii;
		  }//endif
      }//endif
	}
	else if(material==1)
	{
//____account for temperature dependent band gap
//      eel=ee-sieg+1.1241e0/eV0;Si
	  eel=ee-sieg+0.7605e0/eV0;
//	  eel=ee;
      sii=0.0;

      if(kaneiifl)
	  {
//__Kane
		if(eel*eV0>1.1241e0)
		{
		  sii=7.216e11/scrt0*pow((eel*eV0-1.12410),3.540);
		}//endif
      }//endif

      if(thomiifl)
//__Thoma
//		if((eel*eV0>1.128e0)&&(eel*eV0<1.75e0)) si
		if((eel*eV0>0.7605e0)&&(eel*eV0<1.3e0))
		{
//         sii=1.25e12/scrt0*pow((eel*eV0-1.128e0),3);
			sii=1.8e11/scrt0*pow((eel*eV0-0.7605e0),3);
		}
//		else if(eel*eV0>1.75e0)
		else if(eel*eV0>1.3e0)
		{
//         sii=9.49e12/scrt0*(eel*eV0-1.572e0)*(eel*eV0-1.572e0);
			sii=1.6*9.49e12/scrt0*(eel*eV0-0.96e0)*(eel*eV0-0.96e0);
		}//endif
      
      if(sanoiifl)
//__sano
		if(eel*eV0>1.242e0)
		{
		  sii=4.25e11/scrt0*pow((eel*eV0-1.242e0),4.188);
		}//endif
      
      if(fisciifl)
	  {
//___Fischetti
	      if(eel*eV0>1.2e0) 
		  {
		     sii=6.25e10/scrt0*pow(((eel*eV0-1.2e0)/1.20),2);
		  }//endif
          if(eel*eV0>1.8e0)
		  {
			 sii=3e12/scrt0*pow(((eel*eV0-1.80)/1.80),2)+sii;
		  }//endif
          if(eel*eV0>3.45e0)
		  {
			 sii=6.8e14/scrt0*pow(((eel*eV0-3.45e0)/3.45e0),2)+sii;
		  }//endif
      }//endif
	}
      rtnv=sii*iifacelec;

//____End of EIIRATE
      return rtnv;
}

double Band::HIIRATE(double ee)
{

//____Purpose : cal. hole II-rate
//____local variable
      double eel,rtnv;

	if(material==0)
	{
//____account for temperature dependent band gap
      eel=ee-sieg+1.1241e0/eV0;

      if(eel>hiithresh)
	  {
         rtnv=iifachole*pow((eel-hiithresh),hiiexp);
      }
	  else
	  { 
         rtnv=0.0;
      }//endif
	}

	else if(material==1)
	{
//____account for temperature dependent band gap
	  eel=ee-sieg+0.7605e0/eV0;

      if(eel>hiithresh)
	  {
         rtnv=iifachole*pow((eel-hiithresh),hiiexp);
      }
	  else
	  { 
         rtnv=0.0;
      }//endif
	}

//____End of HIIRATE
      return rtnv;
}

double Band::CALSCATTSUM(double eed,int iinital)
{

//____Purpose : interpolate the dos-rate between 
//             two tab-values
//   Parameter : eed=energy of actual electron
//             : ip=particle type

      int itab;
      double intp,rtnv;

//____smallest allowed energy
      if(eed<=0.0)
	  {
         rtnv=sumscatt[0][iinital];
//____largest allowed energy
      }
	  else if(eed>=emax) 
	  {
         rtnv=sumscatt[MTAB][iinital];
      }
	  else
	  { 
         itab=(int)(eed/dtable);  
         intp=(eed-energy[itab])/(energy[itab+1]-energy[itab]);
         rtnv=sumscatt[itab][iinital]+intp 
                  *(sumscatt[itab+1][iinital] 
                     -sumscatt[itab][iinital]);
      }//endif

//____End of CALSCATTSUM
      return rtnv;
}
//====

void Band::CALSCATTOE(double eed,double*sca,int ib)
{

//____Purpose : interpolate the scattering rates between 
//             two tab-values for oxide electrons
//   Parameter : eed=energy of actual hole
//             : sca=array of scattering rates for eed

      int itab,is,iinital,ifinal;
      double intp;//,sca(SCPROE*NBOE)

//   
      iinital=ib-bandof[POXEL];

//____smallest allowed energy
      if(eed<=0e0) 
	  {
         for(ifinal=0;ifinal<nband[POXEL];ifinal++)
		 {////DO ifinal=1,nband[POXEL]
            for( is=0;is<scproe;is++)
			{//DO is=1,scproe
				sca[is+scproe*(ifinal)]=scattoe[is][ifinal][iinital]
                                    *dosoe[is][ifinal][0];
			}//ENDDO
         }//ENDDO
//____largest allowed energy
      }
	  else if(eed>=emax)
	  {
         for(ifinal=0;ifinal<nband[POXEL];ifinal++)
		 {////DO ifinal=1,nband[POXEL]
			for( is=0;is<scproe;is++)
			{//DO is=1,scproe
				sca[is+scproe*(ifinal)]=scattoe[is][ifinal][iinital]
                                    *dosoe[is][ifinal][MTAB];
			}//ENDDO
         }//ENDDO
      }
	  else
	  { 
         itab=(int)(eed/dtable);  
         intp=(eed-energy[itab])/(energy[itab+1]-energy[itab]);
         
         for(ifinal=0;ifinal<nband[POXEL];ifinal++)
		 {////DO ifinal=1,nband[POXEL]
			for( is=0;is<scproe;is++)
			{//DO is=1,scproe
				sca[is+scproe*(ifinal-1)]=scattoe[is][ifinal][iinital]
                *(dosoe[is][ifinal][itab]+intp 
                *(dosoe[is][ifinal][itab+1] 
                   -dosoe[is][ifinal][itab]));
			}//ENDDO
         }//ENDDO
      }//endif
//____End of CALSCATTOE
      return;
}

void Band::CALSCATTH(double eed,double*sca,int ib)
{
//____Purpose : interpolate the scattering rates between 
//             two tab-values for holes
//   Parameter : eed=energy of actual hole
//             : sca=array of scattering rates for eed

      int itab,is,iinital,ifinal;
      double intp;//,sca(SCPRH*NBH)

//   
      iinital=ib-bandof[PHOLE];

//____smallest allowed energy
      if(eed<=0.0) 
	  {
         for(ifinal=0;ifinal<nband[PHOLE];ifinal++)
		 {//         DO ifinal=1,nband[PHOLE]
			for(is=0;is<scprh-1;is++)
			{
			   sca[is+scprh*(ifinal)]=scatth[is][ifinal][iinital]
                                 *dosh[is][ifinal][0];
			}//ENDDO
            sca[scprh-1+scprh*(ifinal)]=scattiih[ifinal][iinital][0];
         }//ENDDO
//____largest allowed energy
      }
	  else if(eed>=emax)
	  {
         for(ifinal=0;ifinal<nband[PHOLE];ifinal++)
		 {//         DO ifinal=1,nband[PHOLE]
			for(is=0;is<scprh-1;is++)
			{//DO is=1,scprh-1
            sca[is+scprh*(ifinal)]=scatth[is][ifinal][iinital]
                                 *dosh[is][ifinal][MTAB];
			}//ENDDO
            sca[scprh-1+scprh*(ifinal)]=scattiih[ifinal][iinital][MTAB];
         }//ENDDO
      }
	  else
	  { 
         itab=(int)(eed/dtable);  
         intp=(eed-energy[itab])/(energy[itab+1]-energy[itab]);
         
		 for(ifinal=0;ifinal<nband[PHOLE];ifinal++)
		 {//         DO ifinal=1,nband[PHOLE]
			for(is=0;is<scprh-1;is++)
			{//DO is=1,scprh-1
            sca[is+scprh*(ifinal)]=scatth[is][ifinal][iinital]
                *(dosh[is][ifinal][itab]+intp 
                *(dosh[is][ifinal][itab+1] 
                   -dosh[is][ifinal][itab]));
			}//ENDDO
            sca[scprh-1+scprh*(ifinal)]=
                      scattiih[ifinal][iinital][itab]+intp 
                *(scattiih[ifinal][iinital][itab+1] 
                   -scattiih[ifinal][iinital][itab]);
         }//ENDDO
      }//endif

//____End of CALSCATTH
      return;
}

void Band::CALSCATTE(double eed,double*sca,int ib)
{

//____Purpose : interpolate the scattering rates between 
//             two tab-values for electrons
//   Parameter : eed=energy of actual electron
//             : sca=array of scattering rates for eed
     int itab,is,iinital,ifinal;
      double intp;//,sca(MSCPRE*NBE)

//   
      iinital=ib-bandof[PELEC];

//____smallest allowed energy
      if(eed<=0e0)
	  {
         for(ifinal=0;ifinal<nband[PELEC];ifinal++)
		 {//DO ifinal=1,nband[PELEC]
			for(is=0;is<scpre-1;is++)
			{//DO is=1,scpre-1
				sca[is+scpre*(ifinal)]=scatte[is][ifinal][iinital]
                                 *dose[is][ifinal][0];
			}//ENDDO
            sca[scpre-1+scpre*(ifinal)]=scattiie[ifinal][iinital][0];
         }//ENDDO
//____largest allowed energy
      }
	  else if(eed>=emax)
	  {
         for(ifinal=0;ifinal<nband[PELEC];ifinal++)
		 {//DO ifinal=1,nband[PELEC]
			for(is=0;is<scpre-1;is++)
			{//DO is=1,scpre-1
				sca[is+scpre*(ifinal)]=scatte[is][ifinal][iinital]
                                 *dose[is][ifinal][MTAB];
			}//ENDDO
            sca[scpre-1+scpre*(ifinal)]=scattiie[ifinal][iinital][MTAB];
         }//ENDDO
      }
	  else
	  { 
         itab=(int)(eed/dtable);  
         intp=(eed-energy[itab])/(energy[itab+1]-energy[itab]);
         
         for(ifinal=0;ifinal<nband[PELEC];ifinal++)
		 {//DO ifinal=1,nband[PELEC]
			for(is=0;is<scpre-1;is++)
			{//DO is=1,scpre-1
				sca[is+scpre*(ifinal)]=scatte[is][ifinal][iinital]
                *(dose[is][ifinal][itab]+intp 
                *(dose[is][ifinal][itab+1] 
                   -dose[is][ifinal][itab]));
			}//ENDDO
            sca[scpre-1+scpre*(ifinal)]=
                      scattiie[ifinal][iinital][itab]+intp 
                *(scattiie[ifinal][iinital][itab+1] 
                   -scattiie[ifinal][iinital][itab]);
         }//ENDDO
      }//endif

//____End of CALSCATTE
      return;
}
//====

double Band::CALDOS(double eed,int ib)
{
//____Purpose : interpolate the dos between 
//             two tab-values
//   Parameter : eed=energy of actual electron
//             : ib=bandindex
      int itab;
      double intp,rtnv;

//____smallest allowed energy
      if(eed<=0e0)
	  {
         rtnv=dos[ib][0];
//____largest allowed energy
      }
	  else if(eed>=emax)
	  {
         rtnv=dos[ib][MTAB];
      }
	  else
	  { 
         itab=(int)(eed/dtable);  
         intp=(eed-energy[itab])/(energy[itab+1]-energy[itab]);
         rtnv=dos[ib][itab]+intp*(dos[ib][itab+1]-dos[ib][itab]);
      }//endif

//____End of CALDOS
      return rtnv;
}
//====

double Band::CALDOSSUM(double eed,int ip)
{

//____Purpose : interpolate the total dos between 
//             two tab-values
//   Parameter : eed=energy of actual electron
//             : ip=particle type

      int itab;
      double rtnv,intp;
//____smallest allowed energy
      if(eed<=0e0)
		  rtnv=sumdos[0][ip];
//____largest allowed energy
      else if(eed>=emax)
          rtnv=sumdos[MTAB][ip];
      else
	  { 
          itab=(int)(eed/dtable);  
          intp=(eed-energy[itab])/(energy[itab+1]-energy[itab]);
          rtnv=sumdos[itab][ip]+intp 
                *(sumdos[itab+1][ip]-sumdos[itab][ip]);
      }//endif

//____End of CALDOSSUM
      return rtnv;
}
//====

double Band::CALDOSTETMAX(double eed,int ib)
{
//____Purpose : the maximum area between two tab-values
//   Parameter : eed=energy of actual electron
//             : ib=bandindex
      int itab;
      double rtnv;
//____smallest allowed energy
      if(eed<=0.0)
	     rtnv=dostetmax[ib][0];
//____largest allowed energy
      else if(eed>=emax)
         rtnv=dostetmax[ib][MTAB];
      else
	  { 
         itab=(int)(eed/dtable);  
         rtnv=MAX(dostetmax[ib][itab],dostetmax[ib][itab+1]);
      }//endif

//____End of CALDOSTETMAX
      return rtnv;
}
//====

void Band::OVERLAP(int it,int&nct,int*ict)
{
// Purpose:   Find cubes which contain the given tetrahedron.
//--------  For each cube a list is built,which gives the tetrahedrons
//            contained in the cube. The cube grid extends from
//            0.7<=x<1.0,0<=y<0.075 and 0<=z<0.075 and is made only
//            for the first conduction band.
//____local variables
      bool foundfl;
      //int it,nct,ict(MCLE),
	  int itc,ibase;
      int ic,icx,icy,icz;
      int icxmin,icymin,iczmin;
      int icxmax,icymax,iczmax;
      double xkl,xkh,ykl,ykh,zkl,zkh;
      double d1,d2,d3,d4;
      double ex,ey,ez,el,ll;
      double xs,ys,zs;
      double xk1,xk2,xk3,xk4;
      double yk1,yk2,yk3,yk4;
      double zk1,zk2,zk3,zk4;
      double xkmin,ykmin,zkmin;
      double xkmax,ykmax,zkmax;
      double dnt1,dnt2,dnt3,dnt4;
      double dnt5,dnt6,dnt7,dnt8;
      double dnt9,dnt10,dnt11,dnt12;
      double dnt13,dnt14,dnt15,dnt16;

//____initializes numbers of cubes within the tetrahedron  
      nct=0;

//____check if tetrahedron is outside of all cubes
   foundfl=false;
   for(itc=0;itc<4;itc++)
   {//DO itc=1,4
        if(ykk[tet[itc][it]]<0.075*a0pi)foundfl=true;
   }//ENDDO
   if(!foundfl)//!foundfl 1
   {//RETURN
      foundfl=false;
      for(itc=0;itc<4;itc++){//DO itc=1,4
         if(xkk[tet[itc][it]]>0.7*a0pi)foundfl=true;
      }//ENDDO
   }//!foundfl 1
   if(!foundfl)///!foundfl 2
   {//RETURN

//____
      xk1=xkk[tet[0][it]];
      xk2=xkk[tet[1][it]];
      xk3=xkk[tet[2][it]];
      xk4=xkk[tet[3][it]];

      yk1=ykk[tet[0][it]];
      yk2=ykk[tet[1][it]];
      yk3=ykk[tet[2][it]];
      yk4=ykk[tet[3][it]];

      zk1=zkk[tet[0][it]];
      zk2=zkk[tet[1][it]];
      zk3=zkk[tet[2][it]];
      zk4=zkk[tet[3][it]];

      xkmin=MIN(xk1,xk2,xk3,xk4);
      xkmax=MAX(xk1,xk2,xk3,xk4);

      ykmin=MIN(yk1,yk2,yk3,yk4);
      ykmax=MAX(yk1,yk2,yk3,yk4);

      zkmin=MIN(zk1,zk2,zk3,zk4);
      zkmax=MAX(zk1,zk2,zk3,zk4);

      icxmin=MAX(0,(int)((xkmin/a0pi-0.7e0)/0.3e0*(MCX)));
      icymin=MAX(0,(int)(ykmin/(a0pi*0.075e0)*(MCY)));
      iczmin=MAX(0,(int)(zkmin/(a0pi*0.075e0)*(MCZ)));

      icxmax=MIN(MCX-1,(int)((xkmax/a0pi-0.7e0)/0.3e0*(MCX)));
      icymax=MIN(MCY-1,(int)(ykmax/(a0pi*0.075e0)*(MCY)));
      iczmax=MIN(MCZ-1,(int)(zkmax/(a0pi*0.075e0)*(MCZ)));

      //ibase=16*(it-1)
	  ibase=it;
      dnt1=datantlin[0][0][it];//(1+ibase)
      dnt2=datantlin[1][0][it];//(2+ibase)
      dnt3=datantlin[2][0][it];//(3+ibase)
      dnt4=datantlin[3][0][it];//(4+ibase)
      dnt5=datantlin[0][1][it];//(5+ibase)
      dnt6=datantlin[1][1][it];//(6+ibase)
      dnt7=datantlin[2][1][it];//(7+ibase)
      dnt8=datantlin[3][1][it];//(8+ibase)
      dnt9=datantlin[0][2][it];//(9+ibase)
      dnt10=datantlin[1][2][it];//(10+ibase)
      dnt11=datantlin[2][2][it];//(11+ibase)
      dnt12=datantlin[3][2][it];//(12+ibase)
      dnt13=datantlin[0][3][it];//(13+ibase)
      dnt14=datantlin[1][3][it];//(14+ibase)
      dnt15=datantlin[2][3][it];//(15+ibase)
      dnt16=datantlin[3][3][it];//(16+ibase)

//____check every cube
      for(icx=icxmin;icx<=icxmax;icx++)
	  {//DO icx=icxmin,icxmax
       for(icy=icymin;icy<=icymax;icy++)
	   {//DO icy=icymin,icymax
        for(icz=iczmin;icz<=iczmax;icz++)
		{//DO icz=iczmin,iczmax
         ic=icz+MCZ*(icy+MCY*icx);//+1 minim=0!!!!!!!!!!!ict[]=ic;ict[] is used as index???????
//_______check if the cube under investigation contains the tetrahedron
         foundfl=false;
         xkl=(0.3e0/(MCX)*(icx)+0.7e0)*a0pi;
         if(xk1>xkl)foundfl=true;
         if(xk2>xkl)foundfl=true;
         if(xk3>xkl)foundfl=true;
         if(xk4>xkl)foundfl=true;
         if(!foundfl)continue;

         foundfl=false;
         xkh=(0.3e0/(MCX)*(icx+1)+0.7e0)*a0pi;
         if(xk1<xkh)foundfl=true;
         if(xk2<xkh)foundfl=true;
         if(xk3<xkh)foundfl=true;
         if(xk4<xkh)foundfl=true;
         if(!foundfl)continue;

         foundfl=false;
         ykl=(0.075e0/(MCY)*(icy))*a0pi;
         if(yk1>ykl)foundfl=true;
         if(yk2>ykl)foundfl=true;
         if(yk3>ykl)foundfl=true;
         if(yk4>ykl)foundfl=true;
         if(!foundfl)continue;

         foundfl=false;
         ykh=(0.075e0/(MCY)*(icy+1))*a0pi;
         if(yk1<ykh)foundfl=true;
         if(yk2<ykh)foundfl=true;
         if(yk3<ykh)foundfl=true;
         if(yk4<ykh)foundfl=true;
         if(!foundfl)continue;

         foundfl=false;
         zkl=(0.075e0/(MCZ)*(icz))*a0pi;
         if(zk1>zkl)foundfl=true;
         if(zk2>zkl)foundfl=true;
         if(zk3>zkl)foundfl=true;
         if(zk4>zkl)foundfl=true;
         if(!foundfl)continue;

         foundfl=false;
         zkh=(0.075e0/(MCZ)*(icz+1))*a0pi;
         if(zk1<zkh)foundfl=true;
         if(zk2<zkh)foundfl=true;
         if(zk3<zkh)foundfl=true;
         if(zk4<zkh)foundfl=true;
         if(!foundfl)continue;
        
//_______check if one of the nodes of the tetrahedron is within the cube
         foundfl=false;
         if((xk1>=xkl)&& 
            (xk1<=xkh)&&
            (yk1>=ykl)&& 
            (yk1<=ykh)&&
            (zk1>=zkl)&& 
            (zk1<=zkh))foundfl=true;
         if((xk2>=xkl)&& 
            (xk2<=xkh)&&
            (yk2>=ykl)&& 
            (yk2<=ykh)&&
            (zk2>=zkl)&& 
            (zk2<=zkh))foundfl=true;
         if((xk3>=xkl)&& 
            (xk3<=xkh)&&
            (yk3>=ykl)&& 
            (yk3<=ykh)&&
            (zk3>=zkl)&& 
            (zk3<=zkh))foundfl=true;
         if((xk4>=xkl)&& 
            (xk4<=xkh)&&
            (yk4>=ykl)&& 
            (yk4<=ykh)&&
            (zk4>=zkl)&& 
            (zk4<=zkh))foundfl=true;
         if(foundfl)
		 {//THEN
            nct=nct+1;
            if(nct>MCLE)
			{
				cout<<"error OVERLAP: nct >MCLE";exit(0);
            }//ENDIF
            ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
            continue;
         }//ENDIF

//_______check if the cube is sorounded by a tetrahedron

//   ..xkl,ykl,zkl
         d1=dnt1*xkl+dnt2*ykl 
           +dnt3*zkl-dnt4;
         d2=dnt5*xkl+dnt6*ykl 
            +dnt7*zkl-dnt8;
         d3=dnt9*xkl+dnt10*ykl 
            +dnt11*zkl-dnt12;
         d4=dnt13*xkl+dnt14*ykl 
            +dnt15*zkl-dnt16;
         if(d1>0.0&&d2>0.0&&
             d3>0.0&&d4>0.0)
		 {//THEN
            nct=nct+1;
            if(nct>MCLE)
			{//THEN
				/*
               WRITE(LUTTYO,*)'OVERLAP: nct >MCLE'
               WRITE(LUOUT ,*)''
               STOP
			 */
				cout<<"error OVERLAP: nct >MCLE";exit(0);
            }//ENDIF
            ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
            continue;
         }//ENDIF

//   ..xkl,ykl,zkh
         d1=dnt1*xkl+dnt2*ykl 
            +dnt3*zkh-dnt4;
         d2=dnt5*xkl+dnt6*ykl 
            +dnt7*zkh-dnt8;
         d3=dnt9*xkl+dnt10*ykl 
            +dnt11*zkh-dnt12;
         d4=dnt13*xkl+dnt14*ykl 
            +dnt15*zkh-dnt16;
         if(d1>0.0&&d2>0.0&&
             d3>0.0&&d4>0.0)
		 {//THEN
            nct=nct+1;
            if(nct>MCLE)
			{
				cout<<"error OVERLAP: nct >MCLE";exit(0);
            }//ENDIF
            ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
            continue;
         }//ENDIF

//   ..xkl,ykh,zkl
         d1=dnt1*xkl+dnt2*ykh 
            +dnt3*zkl-dnt4;
         d2=dnt5*xkl+dnt6*ykh 
            +dnt7*zkl-dnt8;
         d3=dnt9*xkl+dnt10*ykh 
            +dnt11*zkl-dnt12;
         d4=dnt13*xkl+dnt14*ykh 
            +dnt15*zkl-dnt16;
         if(d1>0.0&&d2>0.0&&
             d3>0.0&&d4>0.0)
		 {//THEN
            nct=nct+1;
            if(nct>MCLE)
			{//THEN
				/*
               WRITE(LUTTYO,*)'OVERLAP: nct >MCLE'
               WRITE(LUOUT ,*)'OVERLAP: nct >MCLE'
               STOP
			 */
				cout<<"error OVERLAP: nct >MCLE";exit(0);
            }//ENDIF
            ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
            continue;
         }//ENDIF

//   ..xkl,ykh,zkh
         d1=dnt1*xkl+dnt2*ykh 
            +dnt3*zkh-dnt4;
         d2=dnt5*xkl+dnt6*ykh 
            +dnt7*zkh-dnt8;
         d3=dnt9*xkl+dnt10*ykh 
            +dnt11*zkh-dnt12;
         d4=dnt13*xkl+dnt14*ykh 
            +dnt15*zkh-dnt16;
         if(d1>0.0&&d2>0.0&&
             d3>0.0&&d4>0.0)
		 {//THEN
            nct=nct+1;
            if(nct>MCLE)
			{
				cout<<"error OVERLAP: nct >MCLE";exit(0);
            }//ENDIF
            ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
            continue;
         }//ENDIF

//   ..xkh,ykl,zkl
         d1=dnt1*xkh+dnt2*ykl 
            +dnt3*zkl-dnt4;
         d2=dnt5*xkh+dnt6*ykl 
            +dnt7*zkl-dnt8;
         d3=dnt9*xkh+dnt10*ykl 
            +dnt11*zkl-dnt12;
         d4=dnt13*xkh+dnt14*ykl 
            +dnt15*zkl-dnt16;
         if(d1>0.0&&d2>0.0&&
             d3>0.0&&d4>0.0)
		 {//THEN
            nct=nct+1;
            if(nct>MCLE)
			{
				cout<<"error OVERLAP: nct >MCLE";exit(0);
            }//ENDIF
            ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
            continue;
         }//ENDIF

//   ..xkh,ykl,zkh
         d1=dnt1*xkh+dnt2*ykl 
            +dnt3*zkh-dnt4;
         d2=dnt5*xkh+dnt6*ykl 
            +dnt7*zkh-dnt8;
         d3=dnt9*xkh+dnt10*ykl 
            +dnt11*zkh-dnt12;
         d4=dnt13*xkh+dnt14*ykl
            +dnt15*zkh-dnt16;
         if(d1>0.0&&d2>0.0&&
             d3>0.0&&d4>0.0)
		 {//THEN
            nct=nct+1;
            if(nct>MCLE)
			{
				cout<<"error OVERLAP: nct >MCLE";exit(0);
            }//ENDIF
            ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
            continue;
         }//ENDIF

//   ..xkh,ykh,zkl
         d1=dnt1*xkh+dnt2*ykh 
            +dnt3*zkl-dnt4;
         d2=dnt5*xkh+dnt6*ykh 
            +dnt7*zkl-dnt8;
         d3=dnt9*xkh+dnt10*ykh 
            +dnt11*zkl-dnt12;
         d4=dnt13*xkh+dnt14*ykh 
            +dnt15*zkl-dnt16;
         if(d1>0.0&&d2>0.0&&
             d3>0.0&&d4>0.0)
		 {//THEN
            nct=nct+1;
            if(nct>MCLE)
			{
				cout<<"error OVERLAP: nct >MCLE";exit(0);
            }//ENDIF
            ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
            continue;
         }//ENDIF

//   ..xkh,ykh,zkh
         d1=dnt1*xkh+dnt2*ykh 
            +dnt3*zkh-dnt4;
         d2=dnt5*xkh+dnt6*ykh 
            +dnt7*zkh-dnt8;
         d3=dnt9*xkh+dnt10*ykh 
            +dnt11*zkh-dnt12;
         d4=dnt13*xkl+dnt14*ykh
            +dnt15*zkh-dnt16;
         if(d1>0.0&&d2>0.0&&
             d3>0.0&&d4>0.0)
		 {//THEN
            nct=nct+1;
            if(nct>MCLE)
			{
				cout<<"error OVERLAP: nct >MCLE";exit(0);
            }//ENDIF
            ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
            continue;
         }//ENDIF

//_______intersection between edge 1 of tetrahedron and surface of cube
         ex=xk2-xk1;
         ey=yk2-yk1;
         ez=zk2-zk1;
         el=sqrt(ex*ex+ey*ey+ez*ez);
         ex=ex/el;
         ey=ey/el;
         ez=ez/el;

         if(fabs(ex)>0.0)
		 {//THEN
            ll=(xkl-xk1)/ex;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               ys=ll*ey+yk1;
               zs=ll*ez+zk1;
               if((ys>ykl)&&(ys<ykh)&&
                  (zs>zkl)&&(zs<zkh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
         if(fabs(ex)>0.0)
		 {//THEN
            ll=(xkh-xk1)/ex;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               ys=ll*ey+yk1;
               zs=ll*ez+zk1;
               if((ys>ykl)&&(ys<ykh)&&
                  (zs>zkl)&&(zs<zkh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
         if(fabs(ey)>0.0)
		 {//THEN
            ll=(ykl-yk1)/ey;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               xs=ll*ex+xk1;
               zs=ll*ez+zk1;
               if((xs>xkl)&&(xs<xkh)&&
                  (zs>zkl)&&(zs<zkh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
         if(fabs(ey)>0.0)
		 {//THEN
            ll=(ykh-yk1)/ey;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               xs=ll*ex+xk1;
               zs=ll*ez+zk1;
               if((xs>xkl)&&(xs<xkh)&&
                  (zs>zkl)&&(zs<zkh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
         if(fabs(ez)>0.0)
		 {//THEN
            ll=(zkl-zk1)/ez;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               xs=ll*ex+xk1;
               ys=ll*ey+yk1;
               if((xs>xkl)&&(xs<xkh)&&
                  (ys>ykl)&&(ys<ykh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
         if(fabs(ez)>0.0)
		 {//THEN
            ll=(zkh-zk1)/ez;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               xs=ll*ex+xk1;
               ys=ll*ey+yk1;
               if((xs>xkl)&&(xs<xkh)&&
                  (ys>ykl)&&(ys<ykh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
//_______intersection between edge 2 of tetrahedron and surface of cube
         ex=xk3-xk1;
         ey=yk3-yk1;
         ez=zk3-zk1;
         el=sqrt(ex*ex+ey*ey+ez*ez);
         ex=ex/el;
         ey=ey/el;
         ez=ez/el;

         if(fabs(ex)>0.0)
		 {//THEN
            ll=(xkl-xk1)/ex;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               ys=ll*ey+yk1;
               zs=ll*ez+zk1;
               if((ys>ykl)&&(ys<ykh)&&
                  (zs>zkl)&&(zs<zkh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
         if(fabs(ex)>0.0)
		 {//THEN
            ll=(xkh-xk1)/ex;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               ys=ll*ey+yk1;
               zs=ll*ez+zk1;
               if((ys>ykl)&&(ys<ykh)&&
                  (zs>zkl)&&(zs<zkh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
         if(fabs(ey)>0.0)
		 {//THEN
            ll=(ykl-yk1)/ey;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               xs=ll*ex+xk1;
               zs=ll*ez+zk1;
               if((xs>xkl)&&(xs<xkh)&&
                  (zs>zkl)&&(zs<zkh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
         if(fabs(ey)>0.0)
		 {//THEN
            ll=(ykh-yk1)/ey;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               xs=ll*ex+xk1;
               zs=ll*ez+zk1;
               if((xs>xkl)&&(xs<xkh)&&
                  (zs>zkl)&&(zs<zkh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
         if(fabs(ez)>0.0)
		 {//THEN
            ll=(zkl-zk1)/ez;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               xs=ll*ex+xk1;
               ys=ll*ey+yk1;
               if((xs>xkl)&&(xs<xkh)&&
                  (ys>ykl)&&(ys<ykh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
         if(fabs(ez)>0.0)
		 {//THEN
            ll=(zkh-zk1)/ez;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               xs=ll*ex+xk1;
               ys=ll*ey+yk1;
               if((xs>xkl)&&(xs<xkh)&&
                  (ys>ykl)&&(ys<ykh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
//_______intersection between edge 3 of tetrahedron and surface of cube
         ex=xk4-xk1;
         ey=yk4-yk1;
         ez=zk4-zk1;
         el=sqrt(ex*ex+ey*ey+ez*ez);
         ex=ex/el;
         ey=ey/el;
         ez=ez/el;

         if(fabs(ex)>0.0)
		 {//THEN
            ll=(xkl-xk1)/ex;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               ys=ll*ey+yk1;
               zs=ll*ez+zk1;
               if((ys>ykl)&&(ys<ykh)&&
                  (zs>zkl)&&(zs<zkh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
         if(fabs(ex)>0.0)
		 {//THEN
            ll=(xkh-xk1)/ex;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               ys=ll*ey+yk1;
               zs=ll*ez+zk1;
               if((ys>ykl)&&(ys<ykh)&&
                  (zs>zkl)&&(zs<zkh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
         if(fabs(ey)>0.0)
		 {//THEN
            ll=(ykl-yk1)/ey;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               xs=ll*ex+xk1;
               zs=ll*ez+zk1;
               if((xs>xkl)&&(xs<xkh)&&
                  (zs>zkl)&&(zs<zkh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
         if(fabs(ey)>0.0)
		 {//THEN
            ll=(ykh-yk1)/ey;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               xs=ll*ex+xk1;
               zs=ll*ez+zk1;
               if((xs>xkl)&&(xs<xkh)&&
                  (zs>zkl)&&(zs<zkh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
         if(fabs(ez)>0.0)
		 {//THEN
            ll=(zkl-zk1)/ez;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               xs=ll*ex+xk1;
               ys=ll*ey+yk1;
               if((xs>xkl)&&(xs<xkh)&&
                  (ys>ykl)&&(ys<ykh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
         if(fabs(ez)>0.0)
		 {//THEN
            ll=(zkh-zk1)/ez;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               xs=ll*ex+xk1;
               ys=ll*ey+yk1;
               if((xs>xkl)&&(xs<xkh)&&
                  (ys>ykl)&&(ys<ykh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
//_______intersection between edge 4 of tetrahedron and surface of cube
         ex=xk3-xk2;
         ey=yk3-yk2;
         ez=zk3-zk2;
         el=sqrt(ex*ex+ey*ey+ez*ez);
         ex=ex/el;
         ey=ey/el;
         ez=ez/el;

         if(fabs(ex)>0.0)
		 {//THEN
            ll=(xkl-xk2)/ex;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               ys=ll*ey+yk2;
               zs=ll*ez+zk2;
               if((ys>ykl)&&(ys<ykh)&&
                  (zs>zkl)&&(zs<zkh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
         if(fabs(ex)>0.0)
		 {//THEN
            ll=(xkh-xk2)/ex;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               ys=ll*ey+yk2;
               zs=ll*ez+zk2;
               if((ys>ykl)&&(ys<ykh)&&
                  (zs>zkl)&&(zs<zkh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
         if(fabs(ey)>0.0)
		 {//THEN
            ll=(ykl-yk1)/ey;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               xs=ll*ex+xk2;
               zs=ll*ez+zk2;
               if((xs>xkl)&&(xs<xkh)&&
                  (zs>zkl)&&(zs<zkh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
         if(fabs(ey)>0.0)
		 {//THEN
            ll=(ykh-yk1)/ey;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               xs=ll*ex+xk2;
               zs=ll*ez+zk2;
               if((xs>xkl)&&(xs<xkh)&&
                  (zs>zkl)&&(zs<zkh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
         if(fabs(ez)>0.0)
		 {//THEN
            ll=(zkl-zk2)/ez;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               xs=ll*ex+xk2;
               ys=ll*ey+yk2;
               if((xs>xkl)&&(xs<xkh)&&
                  (ys>ykl)&&(ys<ykh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
         if(fabs(ez)>0.0)
		 {//THEN
            ll=(zkh-zk2)/ez;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               xs=ll*ex+xk2;
               ys=ll*ey+yk2;
               if((xs>xkl)&&(xs<xkh)&&
                  (ys>ykl)&&(ys<ykh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
//_______intersection between edge 5 of tetrahedron and surface of cube
         ex=xk4-xk2;
         ey=yk4-yk2;
         ez=zk4-zk2;
         el=sqrt(ex*ex+ey*ey+ez*ez);
         ex=ex/el;
         ey=ey/el;
         ez=ez/el;

         if(fabs(ex)>0.0)
		 {//THEN
            ll=(xkl-xk2)/ex;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               ys=ll*ey+yk2;
               zs=ll*ez+zk2;
               if((ys>ykl)&&(ys<ykh)&&
                  (zs>zkl)&&(zs<zkh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
         if(fabs(ex)>0.0)
		 {//THEN
            ll=(xkh-xk2)/ex;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               ys=ll*ey+yk2;
               zs=ll*ez+zk2;
               if((ys>ykl)&&(ys<ykh)&&
                  (zs>zkl)&&(zs<zkh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
         if(fabs(ey)>0.0)
		 {//THEN
            ll=(ykl-yk2)/ey;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               xs=ll*ex+xk2;
               zs=ll*ez+zk2;
               if((xs>xkl)&&(xs<xkh)&&
                  (zs>zkl)&&(zs<zkh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
         if(fabs(ey)>0.0)
		 {//THEN
            ll=(ykh-yk2)/ey;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               xs=ll*ex+xk2;
               zs=ll*ez+zk2;
               if((xs>xkl)&&(xs<xkh)&&
                  (zs>zkl)&&(zs<zkh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
         if(fabs(ez)>0.0)
		 {//THEN
            ll=(zkl-zk2)/ez;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               xs=ll*ex+xk2;
               ys=ll*ey+yk2;
               if((xs>xkl)&&(xs<xkh)&&
                  (ys>ykl)&&(ys<ykh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
         if(fabs(ez)>0.0)
		 {//THEN
            ll=(zkh-zk2)/ez;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               xs=ll*ex+xk2;
               ys=ll*ey+yk2;
               if((xs>xkl)&&(xs<xkh)&&
                  (ys>ykl)&&(ys<ykh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
//_______intersection between edge 6 of tetrahedron and surface of cube
         ex=xk4-xk3;
         ey=yk4-yk3;
         ez=zk4-zk3;
         el=sqrt(ex*ex+ey*ey+ez*ez);
         ex=ex/el;
         ey=ey/el;
         ez=ez/el;

         if(fabs(ex)>0.0)
		 {//THEN
            ll=(xkl-xk3)/ex;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               ys=ll*ey+yk3;
               zs=ll*ez+zk3;
               if((ys>ykl)&&(ys<ykh)&&
                  (zs>zkl)&&(zs<zkh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
         if(fabs(ex)>0.0)
		 {//THEN
            ll=(xkh-xk3)/ex;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               ys=ll*ey+yk3;
               zs=ll*ez+zk3;
               if((ys>ykl)&&(ys<ykh)&&
                  (zs>zkl)&&(zs<zkh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
         if(fabs(ey)>0.0)
		 {//THEN
            ll=(ykl-yk3)/ey;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               xs=ll*ex+xk3;
               zs=ll*ez+zk3;
               if((xs>xkl)&&(xs<xkh)&&
                  (zs>zkl)&&(zs<zkh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
         if(fabs(ey)>0.0)
		 {//THEN
            ll=(ykh-yk3)/ey;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               xs=ll*ex+xk3;
               zs=ll*ez+zk3;
               if((xs>xkl)&&(xs<xkh)&&
                  (zs>zkl)&&(zs<zkh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
         if(fabs(ez)>0.0)
		 {//THEN
            ll=(zkl-zk3)/ez;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               xs=ll*ex+xk3;
               ys=ll*ey+yk3;
               if((xs>xkl)&&(xs<xkh)&&
                  (ys>ykl)&&(ys<ykh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
             
         if(fabs(ez)>0.0)
		 {//THEN
            ll=(zkh-zk3)/ez;
            if((ll>0.0)&&(ll<1e0))
			{//THEN
               xs=ll*ex+xk3;
               ys=ll*ey+yk3;
               if((xs>xkl)&&(xs<xkh)&&
                  (ys>ykl)&&(ys<ykh))
			   {//THEN
                  nct=nct+1;
                  if(nct>MCLE)
				  {
					 cout<<"error OVERLAP: nct >MCLE";exit(0);
                  }//ENDIF
                  ict[nct-1]=ic;//nct is number not inedx so use nct-1; 
                  continue;
               }//ENDIF
            }//ENDIF
         }//ENDIF
//10      CONTINUE
        }//ENDDO
       }//ENDDO
      }//ENDDO
  }//!foundfl 2
//____end of OVERLAP
  return;
}

//====
//////////have a rest
void Band::GETINDENS(void)
{
//____Calculate intrinsic carrier density(DOS is only for one spin
// direction)

//____local variables
      int itab;
      double en0,ep0,eel,eelmo;//,CALDOSSUM
      double dosl,doslmo;

      en0=0.0;
      ep0=0.0;

//____calculate intrinsic carrier densities using Boltzmann statistics
      eelmo=energy[0];
      for(itab=1;itab<=MTAB;itab++)
	  {//DO itab=1,MTAB
         eel=energy[itab];
         dosl=CALDOSSUM(eel,PELEC);
         doslmo=CALDOSSUM(eelmo,PELEC);
         en0=en0+((doslmo*eel-dosl*eelmo)
                  *(exp(-eelmo)-exp(-eel))
            +(dosl-doslmo)*((eelmo+1e0)*exp(-eelmo)
                            -(eel +1e0)*exp(-eel)))
             /(eel-eelmo);
         dosl  =CALDOSSUM(eel  ,PHOLE);
         doslmo=CALDOSSUM(eelmo,PHOLE);
         ep0=ep0+((doslmo*eel-dosl*eelmo)
                  *(exp(-eelmo)-exp(-eel))
            +(dosl-doslmo)*((eelmo+1e0)*exp(-eelmo)
                            -(eel +1e0)*exp(-eel  )))
             /(eel-eelmo);

         eelmo=eel;
      }//ENDDO
      en0=2e0*en0;
      ep0=2e0*ep0;
      eni=sqrt(en0*ep0*exp(-sieg));
//____end of GETINDENS
      return;
}

void Band::GETDOS( bool calcdos)
{
	  //get the density of states

      //LOGICAL calcdos
      int  itab,iptype,ip,iband;
      double dum, sum, maxaov;
	  ifstream ftp;

	  //calculate the DOS table
      if(calcdos)
      {
		  ZUDI();
      //ENDIF
	  }

	  //read the DOS table from disk
      if(masterfl)
	  {	
         if(sifl)
		 {
//			ftp.open("zd_si.asc");
			 if(material==0)	ftp.open("c:/montecarlo/input/zd.si.asc");	
			 else if(material==1)	ftp.open("c:/montecarlo/input/zd.ge.trs");
			assert(ftp);
		 }
         else
		 {
			if(gaasfl)
			{ 
				ftp.open("zd_gaas.asc");
				assert(ftp);
			}
         //ENDIF
		 }
         //READ (lutmp) itab
         ftp>>itab;//data is MTAB+1
		 if(itab!=MTAB)
		 {
         //ENDIF
//			 cout<<"error ETABF: wrong number of table entries in file";exit(0);
		 }
		 //read dos
		 for(itab=0;itab<=MTAB;itab++)
		 {
         //DO itab=0,MTAB
            for(iband=0;iband<MNB;iband++)
			{
			//DO iband=1, MNB
               dos[iband][itab]=0.0;
			//ENDDO
			}
            //READ (lutmp) dum, (dos[iband][itab],iband=1, nbt)
			ftp>>dum;
			for(iband=0;iband<nbt;iband++)ftp>>dos[iband][itab];
			//calculate the sum of dos
			for(iptype=0;iptype<NPARTYP;iptype++)
			{
            //DO iptype=1, NPARTYP
               sumdos[itab][iptype]=0.0;
			   for(iband=bandof[iptype];iband<nband[iptype]+bandof[iptype];iband++)
			   {	
               //DO iband=1+bandof[iptype],nband[iptype]+bandof[iptype]
                  sumdos[itab][iptype]=sumdos[itab][iptype]+dos[iband][itab];
				  //ENDDO
			   }
			 //ENDDO
			}
		 //ENDDO
		 }
		 //read maximum area
		 ftp>>itab;
		 for(itab=0;itab<=MTAB;itab++)
		 {
         //DO itab=0,MTAB
			for(iband=0;iband<MNB;iband++)
			{
            //DO iband=1, MNB
               dostetmax[iband][itab]=0.0;
			   //ENDDO
			}
            //READ (lutmp) dum, (dostetmax[iband][itab], iband=1, nbt)
            ftp>>dum;
			for(iband=0;iband<nbt;iband++)ftp>>dostetmax[iband][itab];
		 //ENDDO
		 }
		ftp.close();
      }//ENDIF(masterfl) 

	  //normalize density of states and maximum areas of intersection
	  for(iptype=0;iptype<NPARTYP;iptype++)
	  {	  
      //DO iptype=1, NPARTYP
         DOSMAX[iptype]=0.0;
	  //ENDDO
	  }
	  
	  for(itab=0;itab<=MTAB;itab++)
	  {	
      //DO itab=0,MTAB
		 for(int iband=0;iband<nbt;iband++)
		 {
         //DO iband=1,nbt
            dos[iband][itab]=dos[iband][itab]*eV0*pow(spr0,3.0);
            dostetmax[iband][itab]=dostetmax[iband][itab]*velo0/(spk0*spk0);
		 //ENDDO
		 }
		 for(iptype=0;iptype<NPARTYP;iptype++)
		 {
         //DO iptype=1, NPARTYP
            sumdos[itab][iptype]=sumdos[itab][iptype]*eV0*pow(spr0,3.0);
            DOSMAX[iptype]=MAX(DOSMAX[iptype],sumdos[itab][iptype]);
	     //ENDDO
		 }
	  //ENDDO
	  }

	  //calculate maximum DOS in listfp for Fischetti phonon scattering
	  for(iband=0;iband<nbt;iband++)
	  {	
      //DO iband=1, nbt
		 for(itab=0;itab<MWLEFP;itab++)
		 {
         //DO itab=1, MWLEFP
            sum=0.0;
            maxaov=0.0;
			for(ip=ptlistfpab[itab][iband];ip<=ptlistfpab[itab][iband]+ntlistfpab[itab][iband]-1;ip++)
			{
            //DO ip=ptlistfpab[itab][iband], ptlistfpab[itab][iband] + ntlistfpab[itab][iband] - 1                  
               sum   =sum+maxaovtet[tlistfpab[ip]];
               maxaov=MAX(maxaov,maxaovtet[tlistfpab[ip]]);
               if(ibt[tlistfpab[ip]]!=iband)
			   {	
				   cout<<"error GETDOS: Wrong band";exit(0);
			   }
            //ENDDO
			}
            dosfpab[iband][itab]=48.0/pow((2.0*PI),3.0)*sum;
            maxaovfpab[iband][itab]=maxaov;
            sum=0.0;
            maxaov=0.0;
			for(ip=ptlistfpem[itab][iband];ip<=ptlistfpem[itab][iband]+ntlistfpem[itab][iband]-1;ip++)
			{
            //DO ip=ptlistfpem[itab][iband], ptlistfpem[itab][iband] + ntlistfpem[itab][iband] - 1                  
               sum   =sum +      maxaovtet[tlistfpem[ip]];
               maxaov=MAX(maxaov,maxaovtet[tlistfpem[ip]]);
               if(ibt[tlistfpem[ip]]!=iband)
			   {
				   cout<<"error GETDOS: Wrong band";exit(0);
               //ENDIF
			   }
			//ENDDO
			}
            dosfpem[iband][itab]=48.0/pow((2.0*PI),3.0)*sum;
            maxaovfpem[iband][itab]=maxaov;
		 //ENDDO
		 }
      //ENDDO
	  }
	  //End of GETDOS
      return;
      //END
}

void Band::BUILDLISTS(void)
{
      bool chsfl;
      int itab,ibeg,iend,ipoint,it;
      int ipoins,ipoinz,ipoinss,ipointfpab,ipointfpem;
      int nct,ict[MCLE],ipoinc;

	  int iband,idir;
	  //list of tetraheder (table spacing)
      dlist=eV0/0.040;
      dlists=eV0/0.0020;
/*
	  Build list of tetraheda within a certain energy range
	  ntlist : Number of tetraheda within the energy range
	  ptlist : pointer to the first tetrahedron in a list
	   tlist : list of tetraheda
	   slist : the same list containing only tetraheda from
	  the surface of the BZ
	   zlist : list containing only tetrahedra from
	  the surface kz = 0 (used for injection)
	   clist : list containing tetrahedra which overlap
	  with a given cube
*/
      for(iband=0;iband<MNB;iband++)//minim=0;maxim=MNB-1
	  {
	  	 for(itab=0;itab<MWLE;itab++)//minm=0;maxim=MWLE-1
		 {
            ntlist[itab][iband] = 0;
            nslist[itab][iband] = 0;
		 }
		 nzlist[iband]=0;
		 itab=1;
		 for(itab=0;itab<MWLES;itab++)//minm=0;maxim=MWLES-1
		 {
            ntlists[itab][iband] = 0;
		 }
	  }
	  for(itab=0;itab<MCLE;itab++)//minim=0 maxim=MCLE-1
	     nclist[itab]=0;
	  for(it=0;it<nt;it++)//minim=0;maxim=nt-1
	  {
      //DO it= 1, nt
		 //Count tetraheda per energy range
		 //get minimum and maximum energy index within the tetrahedron
         ibeg=eek[tet[0][it]]*dlist;// + 1;//minim=0;
         iend=eek[tet[3][it]]*dlist;// + 1;//minim=0;
         if(iend>=MWLE)
		 {	
			 cout<<"error ETABF: iend >=MWLE";exit(0);
		 }
         if(ibeg>=0)
		 {
            for(itab=ibeg;itab<=iend;itab++)
			{
			//DO itab = ibeg,iend
			//add tetrahedron to each list within the energy interval
               ntlist[itab][ibt[it]] = ntlist[itab][ibt[it]] + 1;
               chsfl=false;
			//(-4) and (-5) are the surfaces of the wedge which belong to the
			//surface of the BZ.
			   //maybe (-5) and (-6)
               for(idir=0;idir<4;idir++)
			   {
			      if((nlt[idir][it]==-5)||(nlt[idir][it]==-6))
				  {
					  chsfl =true;
				  }
			   }
               if(chsfl)
			   {
				   nslist[itab][ibt[it]] = nslist[itab][ibt[it]] + 1;
			   }
			//ENDDO
			}
         //ENDIF
		 }
         chsfl =false;
		 //(-1) is the surface with kz=0
		 //maybe -2
		 for(idir=0;idir<4;idir++)
		 {	
         //DO idir = 1, 4
            if(nlt[idir][it]==(-2))
			{
				chsfl =true;
			}
		 }
         if(chsfl)
		 {
			 nzlist[ibt[it]] = nzlist[ibt[it]] + 1;
		 }	
		 //list for low energies with finer spacing than tlist
         ibeg =eek[tet[0][it]]*dlists;// ) + 1;
         iend =eek[tet[3][it]]*dlists;// ) + 1;
         if(iend>=MWLES)
		 {
			 iend=MWLES-1;//maxim=
		 } 
         if(ibeg>=0)
		 {
			 for(itab=ibeg;itab<=iend;itab++)
			 {
			 //DO itab = ibeg,iend
               ntlists[itab][ibt[it]] = ntlists[itab][ibt[it]] + 1;
			 }
         //ENDIF
		 }
		 
		 //list of tetrahedra in regular cubes (to find states in the X-Minimum)
         if(ibt[it]==bandof[PELEC])////if(ibt[it]==bandof[PELEC]+1)
		 { 
            OVERLAP(it,nct,ict);
			for(itab=0;itab<nct;itab++)
			{
            //DO itab = 1, nct
               nclist[ict[itab]]+=1;
			}
         }
	  }
	  //set pointers	  
	  ipoint = 0;
      ipoins = 0;
      ipoinz = 0;
      ipoinss = 0;
      ipoinc = 0;
	  
	  //set pointer to point to the element behind the last element of the list
      
	  for(iband=0;iband<nbt;iband++)
	  {	
	  //DO iband=1,nbt
		 for(itab=0;itab<MWLE;itab++)
		 {
         //DO itab=1,MWLE
            ipoint=ipoint + ntlist[itab][iband];
            ptlist[itab][iband]=ipoint;// + 1;
            ipoins=ipoins+nslist[itab][iband];
            pslist[itab][iband]=ipoins;// + 1;
            if(ipoint>=MWLI)
			{	
				cout<<"error ETABF: ipoint > MWLI";exit(0);
			}
            if(ipoins>=MSLI)
			{
				cout<<"error ETABF: ipoins > MSLI";exit(0);
			}
		 }
         ipoinz = ipoinz + nzlist[iband];
         pzlist[iband] = ipoinz;// + 1;
         if(ipoinz>=MZLI)
		 {
			 cout<<"error ETABF: ipoinz > MZLI";exit(0);
		 }
		 for(itab=0;itab<MWLES;itab++)
		 {
         //DO itab=1,MWLES
            ipoinss = ipoinss + ntlists[itab][iband];
            ptlists[itab][iband] = ipoinss;// + 1;
            if(ipoinss>=MWLIS)
			{
				cout<<"error ETABF: ipoinss > MWLIS";exit(0);
			}
		 }
	  }
	  for(itab=0;itab<MCLE;itab++)
	  {	
      //DO itab=1,MCLE
         ipoinc = ipoinc + nclist[itab];
         pclist[itab] = ipoinc;// + 1;
         if(ipoinc>=MCLI)
		 {
			 cout<<"error ETABF: ipoinc > MCLI";exit(0);
		 }
	  }
      
	  for(it=0;it<nt;it++)
	  {
	  //DO it= 1, nt
	     //make list
         ibeg=eek[tet[0][it]]*dlist;// + 1;minim=0;
         iend=eek[tet[3][it]]*dlist;// + 1;minim=0;
	     //reduce pointer by one and assign the elements
		 //after this loop the pointer point to the first element of each list
         if(ibeg>=0)
		 {
            for(itab=ibeg;itab<=iend;itab++)
			{
			//DO itab = ibeg,iend
               ipoint = ptlist[itab][ibt[it]] - 1;
               ptlist[itab][ibt[it]] = ipoint;
               tlist[ipoint] = it;
               chsfl =false;
               for(idir=0;idir<4;idir++)
			   {	
			   //DO idir = 1, 4
                  //if(nlt[idir][it]==(-4)||nlt[idir][it]==(-5))
				  if(nlt[idir][it]==(-5)||nlt[idir][it]==(-6))//be careful!!!!!!!!
				  {
					  chsfl =true;
				  }	
			   }
               if(chsfl)
			   {	
                  ipoins = pslist[itab][ibt[it]] - 1;
                  pslist[itab][ibt[it]] = ipoins;
                  slist[ipoins]=it;
               //ENDIF
			   }
			}
         //ENDIF
		 }
         chsfl =false;
		 for(idir=0;idir<4;idir++)
		 {
         //DO idir = 1, 4
            //if(nlt[idir][it]==-1)
			if(nlt[idir][it]==-2)//nlt[][] be careful!!!!!!!!!!!!!
			{
				chsfl =true;
			} 
         }
         if(chsfl)
		 {
            ipoinz = pzlist[ibt[it]] - 1;
            pzlist[ibt[it]] = ipoinz;
            zlist[ipoinz] = it;
         }
         ibeg=eek[tet[0][it]]*dlists;// ) + 1;
         iend=eek[tet[3][it]]*dlists;// ) + 1;
         if(iend>=MWLES)
		 {
			 iend=MWLES-1;//maxim=MWLES-1;
		 }
         if(ibeg>=0)
		 {
            for(itab=ibeg;itab<=iend;itab++)
			{
			//DO itab = ibeg,iend
               ipoinss = ptlists[itab][ibt[it]] - 1;
               ptlists[itab][ibt[it]] = ipoinss;
               tlists[ipoinss] = it;
			   itab++;
            //ENDDO
			}
         //ENDIF
		 }
         if(ibt[it]==bandof[PELEC])/////////////be careful 3.29
		 {
            OVERLAP(it,nct,ict);
			for(itab=0;itab<nct;itab++)
			{
            //DO itab = 1, nct
               ipoinc = pclist[ict[itab]] - 1;
               pclist[ict[itab]] = ipoinc;
               clist[ipoinc] = it;
			}
         //ENDIF
		 }
	  }
	  //list of tetraheder
      dlistfp = 1.0 / efapeet;

//     Build list of tetraheda within a certain energy range
//     listfpab : This list contains all tetrahedrons which are in an
//     energy interval of the width two times the transversal phonon
//     energy. The spacing of the list is the energy of the transversal phonon.
//     In the case of absorbtion all tetrahedrons are included which contain
//     the inital energy plus the phonon energy. In the case of emisson
//     it is the inital energy minus the phonon energy.

	  for(iband=0;iband<MNB;iband++)
	  {	
      //DO iband=1,MNB
		 for(itab=0;itab<MWLEFP;itab++)
		 {
         //DO itab=1,MWLEFP
            ntlistfpab[itab][iband]=0;
            ntlistfpem[itab][iband]=0;
		 }
	  }	
	  //Count tetraheda per energy range
      for(it=0;it<nt;it++)
	  {
	  //DO it= 1, nt
	     //get minimum and maximum energy index within the tetrahedron
         ibeg =eek[tet[0][it]]*dlistfp;// ) + 1;
         iend =eek[tet[3][it]]*dlistfp;// ) + 1;
         if(iend>=MWLEFP)
		 {
			 cout<<"error ETABF: iend >=MWLEFP";exit(0);
		 }
         if(ibeg>=0)
		 {
			for(itab=ibeg;itab<=iend;itab++)
			{
            //DO itab = ibeg,iend
			   //add tetrahedron to each list within the energy interval
               ntlistfpab[itab][ibt[it]] = ntlistfpab[itab][ibt[it]] + 1;
               ntlistfpem[itab][ibt[it]] = ntlistfpem[itab][ibt[it]] + 1;
			}
            if(ibeg-1>=0)//be careful!!!!!!!!!!!!!!!!!! 
			{
			   //add tetrahedron to each list within the energy interval
               ntlistfpab[ibeg-1][ibt[it]] = ntlistfpab[ibeg-1][ibt[it]]+1;
            //ENDIF
			}
            if(iend+1<MWLEFP)//be careful!!!!!!!!!!!
			{
			   //add tetrahedron to each list within the energy interval
               ntlistfpem[iend+1][ibt[it]] = ntlistfpem[iend+1][ibt[it]]+1;
            //ENDIF
			}
         //ENDIF
		 }
	  }	
	  //set pointers
      ipointfpab = 0;
      ipointfpem = 0;

	  //set pointer to point to the element behind the last element of the list
	  for(iband=0;iband<nbt;iband++)
	  {	
      //DO iband=1,nbt 
		 
		 for(itab=0;itab<MWLEFP;itab++)
		 {
         //DO itab=1,MWLEFP
            ipointfpab = ipointfpab + ntlistfpab[itab][iband];
            ptlistfpab[itab][iband]=ipointfpab;// + 1;
            if(ipointfpab>=MWLIFP)
			{
				cout<<"error ETABF: ipointfpab > MWLIFP";exit(0);
            //ENDIF
			}
            ipointfpem = ipointfpem + ntlistfpem[itab][iband];
            ptlistfpem[itab][iband] = ipointfpem;// + 1;
            if(ipointfpem>=MWLIFP)
			{
				cout<<"error ETABF: ipointfpem > MWLIFP";exit(0);
			}
		 }
	  }	

	  //make list
	  for(it=0;it<nt;it++)
	  {	
      //DO it= 1, nt
         ibeg=eek[tet[0][it]]*dlistfp;// ) + 1;
         iend=eek[tet[3][it]]*dlistfp;// ) + 1;
		 //reduce pointer by one and assign the elements
		 //after this loop the pointer point to the first element of each list
         if(ibeg>=0)
		 {
			for(itab=ibeg;itab<=iend;itab++)
			{
            //DO itab = ibeg,iend
               ipointfpab = ptlistfpab[itab][ibt[it]] - 1;
               ptlistfpab[itab][ibt[it]] = ipointfpab;
               tlistfpab[ipointfpab] = it;
               ipointfpem = ptlistfpem[itab][ibt[it]] - 1;
               ptlistfpem[itab][ibt[it]] = ipointfpem;
               tlistfpem[ipointfpem] = it;
			}
            if(ibeg-1>=0)//be careful !!!!!!!!!!!
			{
               ipointfpab = ptlistfpab[ibeg-1][ibt[it]] - 1;
               ptlistfpab[ibeg-1][ibt[it]] = ipointfpab;
               tlistfpab[ipointfpab] = it;
            //ENDIF
			}
            if(iend+1<MWLEFP)//be careful !!!!!!!!!!
			{ 
               ipointfpem = ptlistfpem[iend+1][ibt[it]] - 1;
               ptlistfpem[iend+1][ibt[it]] = ipointfpem;
               tlistfpem[ipointfpem] = it;
            //ENDIF
			}
         //ENDIF
		 }
	  }	
	  //End of BUILDLISTS
      return;
}
void Band::IELEC(void)
{
//     process elec.scatt command
//____local variables
      bool calcdos,calcband;
      int i,itab,iptype,ireal,ierr,icpu;
      double sck,tcpu;
     // char cscat[CPCVL];

//____set parameters for bandstructure and particle type
      //Typename[PELEC]="electron";
      //Typename[PHOLE]="    hole";
      //Typename[POXEL]="oxidelec";
      typemat[PELEC]=SILICON;
      typemat[PHOLE]=SILICON;
      typemat[POXEL]=OXIDE;

//____set values for energy discretization [0:MTAB] -> [emin,emax]
      dtable=0.0010/eV0;
      emin=0.0;
      emax=(double)(MTAB)*dtable;
      for(itab=0;itab<=MTAB;itab++)
	  {//DO itab=0,MTAB
         energy[itab]=(double)(itab)*dtable+emin;
      }//ENDDO

//____read card parameters
      opcalc   =false;//lval(1,linum)
      calcband =false;//lval(2,linum)
      calcdos  =false;//lval(3,linum)
      kaneiifl =false;//lval(4,linum)
      thomiifl =true;//lval(5,linum)
      fisciifl =false;//lval(6,linum)
      sanoiifl =false;//lval(7,linum)
      iifl    =true;//lval(8,linum)
      seciifl =false;//true;//lval(9,linum)
      tunxfl  =false;//lval(10,linum)
      massfl  =false;//lval(11,linum)
      bhfl    =true;//true;//lval(12,linum)
      frickfl =false;//true;//lval(13,linum)
      bhmrtfl =false;//lval(14,linum)
      injoxfl =false;//lval(15,linum)
      jacophfl=true;//lval(16,linum)
      fiscphfl=false;//lval(17,linum)

//____normalize phonon temperature and deformation potentials
	if(material==0)
	{
      temptag=140.0/T0;//dval( 1,linum)/T0
      templag=215.0/T0;//dval( 2,linum)/T0
      templog=720.0/T0;//dval( 3,linum)/T0
      temptaf=220.0/T0;//dval( 4,linum)/T0
      templaf=550.0/T0;//dval( 5,linum)/T0
      temptof=685.0/T0;//dval( 6,linum)/T0
      dftag  =4.65e9/dpc0;//dval( 7,linum)/dpc0
      dflag  =7.44e9/dpc0;//dval( 8,linum)/dpc0
      dflog  =1.023e11/dpc0;//dval( 9,linum)/dpc0
      dftaf  =2.79e9/dpc0;//dval(10,linum)/dpc0
      dflaf  =1.86e10/dpc0;//dval(11,linum)/dpc0
      dftof  =1.86e10/dpc0;//dval(12,linum)/dpc0
      dfelast=8.7/eV0;//dval(13,linum)/eV0

      temphop=680.0/T0;//dval(14,linum)/T0
      dfhop  =4e10/dpc0;//dval(15,linum)/dpc0
      dfhelast=9.3/eV0;// dval(16,linum)/eV0

      tempoeop=735.0/T0;//dval(17,linum)/T0
      dfoeop  =2.2e10/dpc0;//dval(18,linum)/dpc0
      dfoeelast=25.0/eV0;// dval(19,linum)/eV0

      iifacelec   =0.16;//dval(20,linum)

	  difpr[PELEC]=tempdifprelectron;
	  difpr[PHOLE]=tempdifprhole;
//    difpr[PELEC]=0.16;//dval(21,linum)/200405 interface scattering
//    difpr[PHOLE]=0.35;//dval(22,linum)
	  difpr[POXEL]=0.16;//dval(23,linum)

      sioxbgo     =3.2/eV0;//dval(24,linum)/eV0
      beta        =2.15e-5/eV0*sqrt(field0);//dval(25,linum)/eV0*SQRT(field0)

      iifachole   =1.14e12;//dval(26,linum)
      hiithresh=1.49/eV0;//dval(27,linum)/eV0
      hiiexp   =3.4;//dval(28,linum)
      iifachole  =iifachole*pow(eV0,hiiexp)*time0;

      ephb=1.6;//dval(29,linum)

      dmox  =0.5;//dval(30,linum)

      mell=0.9116;//dval(31,linum)
      melt=0.1946;//dval(32,linum)
      meld=pow((mell*melt*melt),(1.0/3.0));

      efoplow =2.47e10/dpc0;//dval(33,linum)/dpc0
      efophigh=2.97e10/dpc0;//dval(34,linum)/dpc0
      efopee  =6.2e-2/eV0;//dval(35,linum)/eV0

      efaplow =1.7/eV0;//dval(36,linum)/eV0
      efaphigh=2.4/eV0;//dval(37,linum)/eV0
      efapeet =4.43e-2/eV0;//dval(38,linum)/eV0
      efapeel =2.21e-2/eV0;//dval(39,linum)/eV0

      dmc[PELEC]=0.289;//dval(40,linum)
      dmc[PHOLE]=0.349;//dval(41,linum)
      dmc[POXEL]=0.5;//dval(42,linum)
	}
	else if(material==1)
	{
      temptag=120.0/T0;//dval( 1,linum)/T0	//200405
      templag=320.0/T0;//dval( 2,linum)/T0
      templog=430.0/T0;//dval( 3,linum)/T0
      temptaf=120.0/T0;//dval( 4,linum)/T0
      templaf=320.0/T0;//dval( 5,linum)/T0
      temptof=320.0/T0;//dval( 6,linum)/T0
      dftag  =2.0e9/dpc0;//dval( 7,linum)/dpc0
      dflag  =4.0e10/dpc0;//dval( 8,linum)/dpc0
      dflog  =2.6e10/dpc0;//dval( 9,linum)/dpc0
      dftaf  =2.e9/dpc0;//dval(10,linum)/dpc0
      dflaf  =0.0e0/dpc0;//dval(11,linum)/dpc0
      dftof  =0.0e0/dpc0;//dval(12,linum)/dpc0
      dfelast=13.0/eV0;//dval(13,linum)/eV0

      temphop=430.0/T0;//dval(14,linum)/T0
      dfhop  =4.1e10/dpc0;//dval(15,linum)/dpc0
      dfhelast=3.5/eV0;// dval(16,linum)/eV0

      tempoeop=735.0/T0;//dval(17,linum)/T0
      dfoeop  =22e10/dpc0;//dval(18,linum)/dpc0
      dfoeelast=25.0/eV0;// dval(19,linum)/eV0

	  iifacelec   =0.16;//dval(20,linum)

	  difpr[PELEC]=tempdifprelectron;
	  difpr[PHOLE]=tempdifprhole;
//    difpr[PELEC]=0.16;//dval(21,linum)/200405 interface scattering
//    difpr[PHOLE]=0.35;//dval(22,linum)
      difpr[POXEL]=0.16;//dval(23,linum)

      sioxbgo     =3.2/eV0;//dval(24,linum)/eV0
      beta        =2.15e-5/eV0*sqrt(field0);//dval(25,linum)/eV0*SQRT(field0)

      iifachole   =1.4e12;//dval(26,linum)
      hiithresh=0.86/eV0;//dval(27,linum)/eV0
      hiiexp   =3.5;//dval(28,linum)
      iifachole  =iifachole*pow(eV0,hiiexp)*time0;

      ephb=2.6;//dval(29,linum)

      dmox  =0.5;//dval(30,linum)

      mell=1.59;//dval(31,linum)
      melt=0.08;//dval(32,linum)
      meld=pow((mell*melt*melt),(1.0/3.0));

      efoplow =2.47e10/dpc0;//dval(33,linum)/dpc0
      efophigh=2.97e10/dpc0;//dval(34,linum)/dpc0
      efopee  =6.2e-2/eV0;//dval(35,linum)/eV0

      efaplow =1.7/eV0;//dval(36,linum)/eV0
      efaphigh=2.4/eV0;//dval(37,linum)/eV0
      efapeet =4.43e-2/eV0;//dval(38,linum)/eV0
      efapeel =2.21e-2/eV0;//dval(39,linum)/eV0

      dmc[PELEC]=0.289;//dval(40,linum)
      dmc[PHOLE]=0.349;//dval(41,linum)
      dmc[POXEL]=0.5;//dval(42,linum)
	}

      //cscat ="scatt.dat";//cval(1,linum)

//____only one impact ionization model can be specified
      i=0; 
      if (kaneiifl) i=i+1;
      if (thomiifl) i=i+1; 
      if (sanoiifl) i=i+1;
      if (fisciifl) i=i+1;
      if (i>1) 
	  {//THEN 
		  /*
         WRITE (LUTTYO,*) 'IELEC: Only one II-Model'
         WRITE (LUOUT ,*) 'IELEC: Only one II-Model'
         STOP
		 */
      }//ENDIF

//____calculate certain bandstructure values
      if (calcband) 
	  {

	  }//ENDIF

//____read bandstruc-values from files
      READBS(); 

//____build lists for finding states in k-space 
      BUILDLISTS();

//____get density of states
      GETDOS(calcdos);

//____build lists for phonon scattering
      BUILDPHSCATT();

//____intrinsic carrier density
      GETINDENS();

//____find the maximum of scattering rate
      for(iptype=0;iptype<NPARTYP;iptype++)
	  {//DO iptype=1,NPARTYP
         if (nband[iptype]!=0) 
		 {//THEN
            gamma[iptype]=0.0;
            for(itab=0;itab<=MTAB;itab++)
			{//DO itab=0,MTAB
               sck=CALSCATTSUM(energy[itab],bandof[iptype]);
				   //sck=CALSCATTSUM(energy[itab],1+bandof[iptype])
               gamma[iptype]=MAX(gamma[iptype],sck);
            }//ENDDO
		 }//ENDIF
      }//ENDDO

//____end of IELEC
      return;
}

class Partical;
class DevSimulator
{
//----------------------------------------------------------------------
//    COMMON for device parameters
//
//             ij UL    UP    ij+ngpx UR
//                  +-------+
//                  |ij     |
//              LEFT| quad  |RIGHT
//                  | rant  |
//                  +-------+ 
//             ij+1 LL LOW    ij+ngpx+1 LR
//
//A quadrant is labeled with the upper][ left grid point
//----------------------------------------------------------------------
 
//____Variables
//    Grid coordinates in X-directions
       double gridx[MNGPX];
//    Grid coordinates in Y-directions
       double gridy[MNGPY];
//    QUADrant AREA
       double quadarea[MNGP];
//    BOX AREA for each particle type
       double boxarea[MNGP][NPARTYP];
//    dielectri//constant of material
       double eps[MNREGTYP];
//    work function difference between metal and semiconductor
       double phims[MNCONT];
//    charge generated by contact
       double dgen[MNCONT],dgenn[3][MNCONT];
//    charge catched by contact
       double dcatch[MNCONT],dcatchh[3][MNCONT];
//	  catch current
	   double dcurrent[MNCONT],dcurrentt[3][MNCONT];
//    CONTact POTential
       double contpot[MNCONT];
//    Discrete Poisson Equation
       double dpe[2*MNBWP+1][MNGP];
//    electro stati//POTential
       double pot[MNGP];
//    electri//field in x-direction
       double xfield[MNGP];
//    electri//field in y-direction
       double yfield[MNGP];

	   double etemp[MNGP];
//     isegrid的格点数;
	   int isegridnumber;
//    grid point X from ISE
	   double isegridx[MNGP],isegridx1[MNGP];
//    grid point Y from ISE
       double isegridy[MNGP],isegridy1[MNGP];
//    EP point XY from ISE
	   double isepot[MNGP],isepot1[MNGP];
//    EF point XY from ISE
	   double iseefx[MNGP],iseefx1[MNGP];
	   double iseefy[MNGP],iseefy1[MNGP];
	   double iseef[MNGP],iseef1[MNGP];
//    velocity from ISE
	   double iseevx[MNGP],iseevx1[MNGP];
	   double iseevy[MNGP],iseevy1[MNGP];
	   double iseev[MNGP],iseev1[MNGP];
	   double isehvx[MNGP],isehvx1[MNGP];
	   double isehvy[MNGP],isehvy1[MNGP];
	   double isehv[MNGP],isehv1[MNGP];
//	  Ndens point XY from ISE
	   double iseedens[MNGP],iseedens1[MNGP];
	   double isehdens[MNGP],isehdens1[MNGP];
//    Doping from ISE
	   double isedopeconcentration[MNGP],isedopeconcentration1[MNGP];
	   double isedonorconcentration[MNGP],isedonorconcentration1[MNGP];
	   double iseaccepconcentration[MNGP],iseaccepconcentration1[MNGP];
//    Fermi level from ISE
	   double iseEfn[MNGP],iseEfn1[MNGP];
	   double iseEfp[MNGP],iseEfp1[MNGP];
//    band energy from ISE
	   double iseEi[MNGP],iseEi1[MNGP];
	   double iseEc[MNGP],iseEc1[MNGP];
	   double iseEv[MNGP],iseEv1[MNGP];

//    Grid coordinates in X-directions
       double gridxx[MNGPX];
//    Grid coordinates in Y-directions
       double gridyy[MNGPY];

//    imagepot in oxide
       double imagepot[MNGP];
//    electri//field in x-direction [eleminated self forces];
       double xfself[4][MNGP];
//    electri//field in y-direction [eleminated self forces];
       double yfself[4][MNGP];
//    DONOR concentration
       double donor[4][MNGP];
//    ACCEPtor concentration
       double accep[4][MNGP];
//    DOPing set up for right hand side of Poisson Equation
       double doprhspe[MNGP];
//    contains density of interface charge
       double qssdens[MNQSS];
//    DONOR from Galene
       double galdonor[MNGP];
//    ACCEPtor from Galene
       double galaccep[MNGP];
//    intrinsi//carrier concentration
       double eni;
//    SPace CHARge 
       double spchar[MNGP];
//    sum of electron and hole particle density within a quadrant
       double quaddens[MNGP],quadndens[MNGP],quadpdens[MNGP];
	   double quaddenss[3][MNGP],quadndenss[3][MNGP],quadpdenss[3][MNGP];
//    GALene POTential
       double galpot[MNGP];
//    GALene Electron CONentration
       double galecon[MNGP];
//    GALene Hole CONentration
       double galhcon[MNGP];
//    GALene Electron TEMPerature
       double galetemp[MNGP];
//    GALene Hole TEMPerature
       double galhtemp[MNGP];
//    GALene Electron AVALanche generation rate   
       double galeaval[MNGP];
//    GALene Hole AVALanche generation rate   
       double galhaval[MNGP];
//    GALene Electron CURrent in X-direcction
       double galecurx[MNGP];
//    GALene Hole CURrent in X-direcction
       double galhcurx[MNGP];
//    GALene Electron CURrentin Y-direction
       double galecury[MNGP];
//    GALene Hole CURrent in Y-direction
       double galhcury[MNGP];
//    Number of particles inject in each time step by an injection contact
       double nparinj;
//    Current inject by an injection contact
       double injcur;
//    Gradient of the Ramo-Shockley testfunctions
       double xgradrs[MNGP][MNCONT];
       double ygradrs[MNGP][MNCONT];
//    difference of voltage
       double diffvol[MNCONT];
//    Time dependent contact potential
       double ctpotnow[MNCONT][MNDT];
//    the time to start changing voltage
       double st[MNCONT];
//    the time to end changing voltage
       double ft[MNCONT];
//    capacitance between two contacts[for RS-theory];
       double cap[MNCONT][MNCONT];
//    Number of GridPoints in X-direction
      //int ngpx;
//    Number of GridPoints in Y-direction
       int ngpy;
//    Number of GridPoints
       int ngp;
//    Number of Band With Points
       int nbwp;
//    contains material for each quadrant
       int mat[MNGP];
//    Number of the region to which the grid points belong
       int numreg[MNGP];
//    contains types of contacts for each quadrant
       int cont[4][MNGP];
//    contains type of region [material];
       int regtype[MNREG];
//    conttype : contains TYPE of CONTact
       int conttype[MNCONT];
//    contains CONTACT POSition with respect to quadrant
       int contpos[MNCONT];
//    index of CONTact attached to grid point
       int gridcont[MNGP];
//    number of particles generated by contact
       int ngen[MNCONT];
//    number of particles catched by contact
       int ncatch[MNCONT];
//    Number of REGions
       int nreg;
//    Index of Begin of REGion in x-direction
       int ibreg[MNREG];
//    Index of End   of REGion in x-direction
       int iereg[MNREG];
//    index of Begin of REGion in y-direction
       int jbreg[MNREG];
//    index of End   of REGion in y-direction
       int jereg[MNREG];
//    Number of CONTacts
       int ncont;
//    Index of Begin of CONTact in x-direction
       int ibcont[MNCONT];
//    Index of End   of CONTact in x-direction
       int iecont[MNCONT];
//    index of Begin of CONTact in y-direction
       int jbcont[MNCONT];
//    index of End   of CONTact in y-direction
       int jecont[MNCONT];
//    Number of interface charges
       int nqss;
//    Index of Begin of interface charge in x-direction
       int ibqss[MNQSS];
//    Index of End   of interface charge in x-direction
       int ieqss[MNQSS];
//    index of Begin of interface charge in y-direction
       int jbqss[MNQSS];
//    index of End   of interface charge in y-direction
       int jeqss[MNQSS];
//    Index of Begin of injection contact in x-direction
       int ibinj;
//    Index of End of injection contact in x-direction
       int ieinj;
//    Index of Begin of injection contact in y-direction
       int jbinj;
//    Index of End of injection contact in y-direction
       int jeinj;
//    Direction of particle injection for injection contact
       int injdir;
//    particle type to be injected by injection contact
       int itypinj;
//    RULES of MOTion for quadrant borders
       int motrules[4][MNGP];
//    flag for selfforce correction in each region
       bool selfforcefl[MNREG];
//    give current by hand for injection contact
       bool cbhfl;
//    SI-CONT in equilibrium
       bool equisicfl[MNCONT];
//    REGion NAME
      char regname[MNREG][10];
//    CONTact NAME
      char contname[MNCONT][10];
	  //miltrefreh
//     Number of Desired Particles Per Energy Range
       int ndpper[MNPSC];

//     Number of ReFreshs Per Energy Range [statistics];
       int nrfper[MNPSC];

//     scalar Index Pointer for each Energy Range
       int ierp[MNER][MNGP+1][NPARTYP];

//     Number of Energy ranges [ReFresh];
       int nerrf[MNGP+1][NPARTYP];

//     Number of Phase Space Cells
       int npsc;

//     pointer to the first particle in chain for the given energy range
       int ifirst[MNPSC];

//     pointers to the next particle
       int inext[MNPAR];

//     current number or particles per energy range
       int npper[MNPSC];

//     current sum of particle weights
       double sumpc[MNPSC];

//     current sum of the squares of particle weight
       double sumquad[MNPSC];

//     Ratio for PARticle criterium
       double rpar;

//     Ratio for sQUAre statistical weight criterium
       double rquad;

//     maximum energy for k-space partitioning
       double erfmax[MNGP+1][NPARTYP];

//     all particles with a charge less than pckill times pcaver are eliminated
       double pckill;

//     average particle charge
       double pcaver;
//     variables for the statistics of substrate current evaluated at the
//     end of the time steps
	   double iicurts[NPARTYP];
       double squiicurts[NPARTYP], meaniicurts[NPARTYP];

//     variables for the statistics of substrate current evaluated by
//     just before scattering statistics
       double iicurjb[NPARTYP];
       double meaniicurjb[NPARTYP], squiicurjb[NPARTYP];
//     statistics for potetnial
       double statpot[MNGP];

//     statistics over the grid
       double statis[NSTAT][MNGP][NPARTYP];

//     scalar version of array statis
      // double statlin(NSTAT*MNGP*NPARTYP)

//     energy distribution function
       double edf[NEREDF][MNGP][NPARTYP];

//     maximum energy for EDF sampling
       double edfemax[NPARTYP];

//     scalar version of array edf
      // double edflin(NEREDF*MNGP*NPARTYP)

//     Si/SiO2 boundary statistics in x and y direction
       double bedfx[NEREDF][MNGPX][NPARTYP];
       double bedfy[NEREDF][MNGPY][NPARTYP];

//     particle current densities (sampled whenever a particle changes quadrant)
       double curx[MNGP][NPARTYP],cury[MNGP][NPARTYP];


//     variables for the statistics of homogeneous EDF evaluated at the end
//     of the time steps
       double edishom[MTAB][NPARTYP];
       double meanedishom[MTAB][NPARTYP];
       double squedishom[MTAB][NPARTYP];

//     variables for the statistics of homogeneous EDF evaluated by
//     just before scattering statistics
       double edfjb[MTAB][NPARTYP];
       double meanedfjb[MTAB][NPARTYP];
       double squedfjb[MTAB][NPARTYP];

//     variables for the statistics of GHDM quantities evaluated by
//     just before scattering statistics
//     (nexp is a normalization quantity)
       double nexp[NPARTYP];
       double ghexp[NGHEXP][NPARTYP];
       double meanghexp[NGHEXP][NPARTYP];
       double squghexp[NGHEXP][NPARTYP];

//     variables for the statistics of contact currents evaluated 
       double curcont[MNCONT],meancurcont[MNCONT];
	   double curcontt[3][MNCONT],meancurcontt[3][MNCONT];
       double squcurcont[MNCONT];
//	int ngpx,ngpy;
//    double gridx[100],gridy[100];
	  
//   $Source: /export/home3/junge/monte/SRC/RCS/pvars.com,v $
//   $Author: junge $
//   $Revision: 1.6 $     $Date: 1995/10/24 06:54:01 $

//-----file: P V A R S . // O M-------------------------------------------
//
//  All Rights Reserved, Copyright (//) FUJITSU Ltd. 1995
//
//-----------------------------------------------------------------------
//     Common for particle variables
//        IVPP    : Integer Variables Per Particle
//          1     : band number
//          2     : cube number
//          3     : tetraeder number
//          4     : number of symetrie-operation
//
//        DVPP    : Double Variables Per Particle
//          1,2,3 : x,y,z coordinates of k-space
//          4     : particle charge/weight
//          5     : left propagation time (if secondary particle)
//
//        MNPAR   : maximum number of particles
//
//        ifield  : integer field for particle variables
//        dfield  : double field for particle variables
//
//        npar    : number of particles
//-----------------------------------------------------------------------
       double dfield[DVPP][MNPAR];
       int ifield[IVPP+1][MNPAR],npar[NPARTYP],npar0;

	  void GETCONTPOT(void);
	  void CALCFIELD(void);
	  void CALCRHO(bool wrfl);
	  void CALCPOT(bool wrfl);
	  void CALCCAP(void);
	  void CALCRS(void);
	 
	  bool CONVER(int idt);
	  void SETDT(void);
	  void GETIMAGEPOT(Band &bd);
	  void ZEROJB(void);
	  void EINIT(Partical &par);
	  void ELEC(int idt,bool wrfl,bool averfl,Partical &par,Band &bd);
	  void ADDJB(int stridejb,Band &bd);
	  void REPLACE(int ij,int iptype,int &ifree,int *nempty,int *nrefresh,Partical &par);
	  void REFRESH(bool wrfl,Band &bd,Partical &par);
	  void AVERINIT(void);
	  void output(void);
	  void outputcut(void);
	  void OUTPUTSCATTER(void);
	  
public:

	// cut position
	int uml,umr,umu,umd;

	int bhscatter,phscatter;
	int sscatter,srscatter,spscatter,siscatter;
	int oldsscatter;

	double qpot[MNGP][NPARTYP];
	double qtarray[MNGP][5];
	double xqf[MNGP][NPARTYP],yqf[MNGP][NPARTYP];


	   double mingridx1[MNGP],mingridy1[MNGP],minpot1[MNGP],minedens1[MNGP],minhdens1[MNGP],
		   minefx1[MNGP],minefy1[MNGP],minevx1[MNGP],minevy1[MNGP],minEfn1[MNGP],minEfp1[MNGP],
		   minEc1[MNGP],minEv1[MNGP];
	   double mingridx2[MNGP],mingridy2[MNGP],minpot2[MNGP],minedens2[MNGP],minhdens2[MNGP],
		   minefx2[MNGP],minefy2[MNGP],minevx2[MNGP],minevy2[MNGP],minEfn2[MNGP],minEfp2[MNGP],
		   minEc2[MNGP],minEv2[MNGP];
	   double mingridx3[MNGP],mingridy3[MNGP],minpot3[MNGP],minedens3[MNGP],minhdens3[MNGP],
		   minefx3[MNGP],minefy3[MNGP],minevx3[MNGP],minevy3[MNGP],minEfn3[MNGP],minEfp3[MNGP],
		   minEc3[MNGP],minEv3[MNGP];
	   double mingridx4[MNGP],mingridy4[MNGP],minpot4[MNGP],minedens4[MNGP],minhdens4[MNGP],
		   minefx4[MNGP],minefy4[MNGP],minevx4[MNGP],minevy4[MNGP],minEfn4[MNGP],minEfp4[MNGP],
		   minEc4[MNGP],minEv4[MNGP];
//用来判断某一点在ISE网格内部的位置，即找出其所在的矩形。
	   bool minflag1,minflag2,minflag3,minflag4;
	   double xx[4],yy[4];
	   int minn1,minn2,minn3,minn4;
	   double bpx1up[500],bpx1down[500],bpy1[500];
	   double bpx2down[500],bpx2up[500],bpy2[500];
	   double bpy3left[500],bpy3right[500],bpx3[500];
	   double bpy4right[500],bpy4left[500],bpx4[250];

	   double edens[MNGP],hdens[MNGP];
	   int ixmin1,ixmax1,jymin1,jymax1;
	   int ixmin2,ixmax2,jymin2,jymax2;
	   int leftcatch,rightcatch,leftgen,rightgen,subcatch;
	   int upcatch,downcatch,upgen,downgen;

	//    Number of GridPoints in X-direction
    int ngpx;
	friend class Partical;
	void BIAST(char *name,double voltage,double delvol,double stime, double etime);
	void INIST(Partical &par);
	void OUTPUTMOTRULES(void);

	void statnoise(void);

	//quantum effect
	void QUANTPOT(int instat);
	void QUANTARRAY(void);
	void QUANTFIELD(void);
	//quantum effect end

	void READINFILE(char *file);
	void GETCUTPOSITION(void);

	void GETISEDATA(void);
	void DATAFROMISETOISEGOOD(void);
	void DATAFROMISETOMC(void);

	void CALCFERMI(void);
	void CALCBTBTDENS(void);
	void CALCRDENS(void);
	void recombination(void);

	void STRUC(int rdwr,char *file);
	void CONFI(int rdwr,char *file);
	//void RUNST(Partical &par,Band &bd);
	 void RUN(Partical &par,Band &bd);
};