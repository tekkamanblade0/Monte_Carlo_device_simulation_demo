// ..key number definition
      const int INIKEY=1, GRIKEY=2, REGKEY=3,
				CONKEY=4, DOPKEY=5, QSSKEY=6,
		        STRKEY=7, ELEKEY=8,
			    BIAKEY=9, 
		        RUNKEY=11, FILKEY=12, MULKEY=10,
			    ENDKEY=18, WRTKEY=13, OSVKEY=14,
		        EDFKEY=15, BDFKEY=16, INJKEY=17;

// ..Integer Variables Per Particle
      const int IVPP=3;
  //    PARAMETER (IVPP=3)

// ..Double Variables Per Particle
      const int DVPP=8;
      //PARAMETER (DVPP=8)

// ..Maximum Number of PARticles
      const int MNPAR=400000;
      //PARAMETER (MNPAR=400000)
//    PARAMETER (MNPAR=300000)

// ..const int Variables of All Particles
      const int IVAP=IVPP*MNPAR;
      //PARAMETER (IVAP=IVPP*MNPAR)

// ..Double Variables of All Particles
      const int DVAP=DVPP*MNPAR;

//   ..Parameters for particle type
     const int NPARTYP=3,PELEC=0,PHOLE=1,POXEL=2;

//  Maximum number of scattering processes (electrons)
     const int MSCPRE=14;
//  Maximum number of scattering processes (holes)
     const int MSCPRH=4;

//  Maximum number of scattering processes (oxide electrons)
     const int MSCPROE=3;

//  Number of table elements for the doping dependent frickel factor
      const int NFITB=4;

//  maximum number of conduction bands
      const int NBE=4;
 
//  maximum number of valence bands
      const int NBH=3;
 
//  maximum number of conduction bands in the oxide
      const int NBOE=1;
 
//  maximum number energy bands
      const int MNB=NBE+NBH+NBOE;

//  maximum number of tetraeder in the onedimensional list
      const int MWLI=1470000;//470000;
 
//  maximum number of energy-fields for tetrdaederlist
      const int MWLE=325;
 
//  maximum number of tetraeder in the onedimensional list
//  of tetrahedrons at the surface of the wedge
      const int MSLI=36000;//16000;
 
      const int MZLI=5500;//4500;
//  number of cubes in x-direction
      const int MCX=25;
 
//  number of cubes in y-direction
      const int MCY=14;
 
//  number of cubes in z-direction
      const int MCZ=14;

//  maxnumber of energy dependent dos-table and scattering-table
      const int MTAB=15000;//7500;
//	  const int MTAB1=2000;
 
//  the same as above for Fischetti type phonon scattering
      const int MWLIFP=1080000;//480000;
 
//  the same as above for Fischetti type phonon scattering
      const int MWLEFP=325;

//  the same as above for a finer list
      const long int MWLIS=301000;//101000;
 
//  the same as above for a finer list
      const int MWLES=60;

//  the same as above for the cube oriented list
      const int MCLI=50000;
 
//  the same as above for the cube oriented list
      const int MCLE=MCX*MCY*MCZ;

//  Maximum number of k-space points
      const int MNK=36000;//13000;
 
//  Maximum number of tetraheda
      const int MNT=175000;//56000;
//____Parameters
//    Maximum Number of GridPoints in X-direction
      const int MNGPX=65;

//    Maximum Number of GridPoints in Y-direction
      const int MNGPY=120;

//    Maximum Number of GridPoints
      const int MNGP=MNGPX*MNGPY;

//    Maximum Number of Band Width Points
      const int MNBWP=MNGPX;
//    Maximum Number of X Y cut
	  const int MAXXCUTNUMBER=100;
	  const int MAXYCUTNUMBER=100;
//    Token for a material that represents nothing (Neumann boundary conditions
//    on the surface for every thing
      const int VOID=0;

//    Token for vacuum
      const int VACUUM=1;

//    Token for oxide
      const int OXIDE=2;

//    Token for silicon
      const int SILICON=3;

//    Token for surface at Si02-Si
      const int SILINE=4;

//    Token for surface at Polygate-SiO2
      const int GATELINE=5;

//    Token for no cont
      const int NOCONT=0;

//    Token for silicon contact
      const int SICONT=1;

//    Token for a gate contact
      const int GATE=2;

//    Token for a catch contact
      const int CCONT=3;

//    Maximum Number of REGions
      const int MNREG=16;

//    Maximum Number of REGions TYPs
      const int MNREGTYP=6;

//    Maximum Number of CONTacts
      const int MNCONT=16;

//    Maximum Number of CONTact TYPs
      const int MNCONTTYP=2;

//    Maximum Number of interface charges
      const int MNQSS=16;

//    Token for direction up
      const int UP=0;

//    Token for direction right
      const int RIGHT=1;

//    Token for direction low
      const int LOW=2;

//    Token for direction left
      const int LEFT=3;

//    Token for upper left corner
      const int UL=1-1;

//    Token for upper right corner
      const int UR=3-1;

//    Token for lower right corner
      const int LR=4-1;

//    Token for lower left corner
      const int LL=2-1;

//    Token for motion rule error
      const int MOTERR=0;

//    Token for motion rule pass (particle passes from one quadrant to another)
      const int PASS=1;

//    Token for motion rule reflect (particle is reflected at the boundary)
      const int REFLECT=2;

//    Token for motion rule scattox (particle is reflected or scattered at the
//    SiO2 boundary)
      const int SCATTOX=3;

//    Token for motion rule period (particle is moved to the opposite boundary 
//    of the quadrant)
      const int PERIOD=4;

//    Token for motion rule generate (new particle is generated at the other 
//    side of the boundary, original particle experiences period boundary 
//    condition)
      const int GENERATE=5;
      const int GENREF=8;

//    Token for motion rule catch (particle is annihilated)
      const int CATCH=6;

//    Token for motion rule catchgate (particle is annihilated by a gate)
      const int CATCHGATE=7;

//    Maximum number of time slice
      const int MNDT=1000000;
//
//-----------------------------------------------------------------------
//     Common for silicon parameters
//
//   ..silicon bulk parameters:
//     sia0=5.43d-10  ! m        ! lattice parameter
//     sirho=2.33d3   ! kg/m**3  ! crystal density
//     siul=9.0d3     ! m/s      ! longitudinal sound velocity
//     siut=5.3d3     ! m/s      ! transverse sound velocity
//     a0pi           ! 1/m      ! 2*PI/sia0
//     sieg           ! eV       ! silicon band gap
//
//------------------------------------------------------------------------
double sia0,sirho,siul,siut,a0pi,sieg;

//-----------------------------------------------------------------------
//     Common for normalization
//
//        T0     = temperature . . . . . . . . . [K]
//        eV0    = energy. . . . . . . . . . . . [eV]
//        em0    = electon rest mass . . . . . . [kg]
//        hq0    = Planck constant . . . . . . . [eVs]
//        ec0    = electron charge . . . . . . . [As]
//        rmom0  = momentum. . . . . . . . . . . [eVs/m]
//        spr0   = r-space . . . . . . . . . . . [m]
//        spk0   = k-space . . . . . . . . . . . [1/m]
//        time0  = time. . . . . . . . . . . . . [s]
//        velo0  = velocity. . . . . . . . . . . [m/s]
//        cvr    = speed of light / velo0
//        pot0   = el. potential . . . . . . . . [V]
//        field0 = el. field . . . . . . . . . . [V/m]
//        conc0  = concentration . . . . . . . . [1/m**3]
//        dens   = mass density. . . . . . . . . [kg/m**3]
//        dpc0   = optical deformationpotential. [eV/m]
//        scrt0  = scattering rate . . . . . . . [1/s]
//        curr0  = terminal current (2D) . . . . [A/m]
//-----------------------------------------------------------------------

      double T0,eV0,em0,hq0,ec0,rmom0,spr0,spk0,time0,velo0; 
      double cvr,pot0,field0,conc0,dens0,dpc0,scrt0,curr0;

//-----------------------------------------------------------------------
//     physical constants:
//        ec     = electron charge         [ As ]
//        em     = electron rest mass      [ kg ]
//        planck = Planck s constant       [ eVs ]
//        clight = speed of ligth          [ m/s ]
//        avoga  = Avogadro sconstant      [ 1/mol ]
//        boltz  = Boltzmann sconstant     [ eV/K ]
//        eps0   = dielectric constant     [ As/(Vm; ]
//        fsc    = fine structure constant
//-----------------------------------------------------------------------

      const double EC     = 1.602192e-19;
      const double EM     = 9.109558e-31;
      const double PLANCK = 6.582183e-16;
      const double CLIGHT = 2.997925e8;
      const double AVOGA  = 6.022169e23;
      const double BOLTZ  = 8.617084e-5;
      const double EPS0   = 8.854185e-12;
      const double FSC    = EC/(PLANCK*CLIGHT*4.0*PI*EPS0);

//-----------------------------------------------------------------------
//     Maximum Number of Energy Ranges
      const int MNER=70;
//     Maximu Number of Phase Space Cells
      const int MNPSC=MNGP*10;
//     Maximu Number of points
	  const int MNREGION=10;
//     Maximu Number of times
	  const int MNTIMES=150000;
//      PARAMETER (MNER = 70, MNPSC=MNGP*10)
	  int iseedl;
