
//   $Source: /export/home3/junge/monte/SRC/RCS/paver.com,v $
//   $Author: junge $
//   $Revision: 1.41 $     $Date: 1997/01/14 00:33:26 $

//-----file: P A V E R . // O M-------------------------------------------
//
//  All Rights Reserved, Copyright (//) FUJITSU Ltd. 1995
//
//-----------------------------------------------------------------------
//
//     COMMON for particle averages
//-----------------------------------------------------------------------

//     Number of STATistic quantities
      const int NSTAT=21;
 
//     Number of Energy Ranges of Energy Distribution Function
      const int NEREDF=50;

//     number of GHDM statistics quantities 
      const int NGHEXP=21;

//     particle density
      const int ARHO=1;

//     Velocity in X-direction
      const int AVX=2;

//     Velocity in Y-direction
      const int AVY=3;

//     Velocity in Z-direction
      const int AVZ=4;

//     Velocity square
      const int AVV=5;

//     Energy
      const int AE=6;

//     Impact Ionization
      const int AII=7;

//     X-Velocity times X-Velocity
      const int AVXX=8;

//     X-Velocity times Y-Velocity
      const int AVXY=9;

//     X-Velocity times Z-Velocity
      const int AVXZ=10;

//     Y-Velocity times Y-Velocity
      const int AVYY=11;

//     Y-Velocity times Z-Velocity
      const int AVYZ=12;

//     Z-Velocity times Z-Velocity
      const int AVZZ=13;

//     Energy times velocity in x-direction 
      const int AEVX=14;

//     Energy times velocity in y-direction 
      const int AEVY=15;


//     invers mass tensor
      const int AIMXX=16;
      const int AIMXY=17;
      const int AIMXZ=18;
      const int AIMYY=19;
      const int AIMYZ=20;
      const int AIMZZ=0;//21;

//     velocity in the direction of the electric field for GHDM quantities
      const int GHV=1;

//     x-velocity for GHDM quantities
      const int GHVX=2;

//     y-velocity for GHDM quantities
      const int GHVY=3;

//     z-velocity for GHDM quantities
      const int GHVZ=4;

//     velocity times velocity for GHDM quantities
      const int GHVV=5;

//     energy for GHDM quantities
      const int GHE=6;

//     energy times velocity for GHDM quantities
      const int GHEV=7;

//     trace of inverse mass tensor
      const int GHSIM=8;

//     trace of inverse mass tensor times energy plus velocity dyadic velocity
      const int GHSIMEVV=9;

//     Stratton mobility tensor
      const int GHMXX=10;
      const int GHMXY=11;
      const int GHMXZ=12;
      const int GHMYY=13;
      const int GHMYZ=14;
      const int GHMZZ=15;

//     Stratton diffusion tensor
      const int GHDXX=16;
      const int GHDXY=17;
      const int GHDXZ=18;
      const int GHDYY=19;
      const int GHDYZ=20;
      const int GHDZZ=0;//21;
/*
      PARAMETER (NSTAT=21)
      PARAMETER (NEREDF = 50)
      PARAMETER (NGHEXP = 21)

      PARAMETER (ARHO=1)
      PARAMETER (AVX =2)
      PARAMETER (AVY =3)
      PARAMETER (AVZ =4)
      PARAMETER (AVV =5)
      PARAMETER (AE  =6)
      PARAMETER (AII =7)
      PARAMETER (AVXX=8)
      PARAMETER (AVXY=9)
      PARAMETER (AVXZ=10)
      PARAMETER (AVYY=11)
      PARAMETER (AVYZ=12)
      PARAMETER (AVZZ=13)
      PARAMETER (AEVX=14)
      PARAMETER (AEVY=15)
      PARAMETER (AIMXX=16)
      PARAMETER (AIMXY=17)
      PARAMETER (AIMXZ=18)
      PARAMETER (AIMYY=19)
      PARAMETER (AIMYZ=20)
      PARAMETER (AIMZZ=21)
      PARAMETER (GHV=1, GHVX=2, GHVY=3, GHVZ=4, GHVV=5, 
     >           GHE=6, GHEV=7, GHSIM=8, GHSIMEVV=9)
      PARAMETER (GHMXX=10, GHMXY=11, GHMXZ=12, 
     >           GHMYY=13, GHMYZ=14, GHMZZ=15,
     >           GHDXX=16, GHDXY=17, GHDXZ=18, 
     >           GHDYY=19, GHDYZ=20, GHDZZ=21)
*/
//     variables for the statistics of substrate current evaluated at the
//     end of the time steps
     // double iicurts[NPARTYP];
     // double squiicurts[NPARTYP], meaniicurts[NPARTYP];

//     variables for the statistics of substrate current evaluated by
//     just before scattering statistics
     // double iicurjb[NPARTYP];
     // double meaniicurjb[NPARTYP], squiicurjb[NPARTYP];
/*
//     statistics for potetnial
      double statpot[MNGP];

//     statistics over the grid
      double statis[NSTAT][MNGP][NPARTYP];

//     scalar version of array statis
      //double statlin(NSTAT*MNGP*NPARTYP)

//     energy distribution function
      double edf[NEREDF][MNGP][NPARTYP];

//     maximum energy for EDF sampling
      double edfemax[NPARTYP];

//     scalar version of array edf
      //double edflin(NEREDF*MNGP*NPARTYP)

//     Si/SiO2 boundary statistics in x and y direction
      double bedfx[NEREDF][MNGPX][NPARTYP];
      double bedfy[NEREDF][MNGPY]NPARTYP];

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
      double squcurcont[MNCONT];
/*
      COMMON /paver/ iicurjb, meaniicurjb, squiicurjb,
     >               iicurts, meaniicurts, squiicurts,
     >               statpot, statis, edf, edfemax, bedfx, bedfy,
     >               curx, cury, 
     >               edishom, meanedishom, squedishom,
     >               nexp,
     >               edfjb, meanedfjb, squedfjb,
     >               curcont, meancurcont, squcurcont,
     >               ghexp, meanghexp, squghexp

      EQUIVALENCE (statis, statlin)
      EQUIVALENCE (edf, edflin)
	  */
