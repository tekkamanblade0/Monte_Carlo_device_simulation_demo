	  double evelox[MNGP],eveloy[MNGP];
	  double hvelox[MNGP],hveloy[MNGP];
	  double eveloxsum[MNGP],eveloysum[MNGP];
	  double hveloxsum[MNGP],hveloysum[MNGP];
	  double enumber[MNGP],hnumber[MNGP];

	  double eveloxx[3][MNGP],eveloyy[3][MNGP];
	  double hveloxx[3][MNGP],hveloyy[3][MNGP];
	  double eveloxsumm[3][MNGP],eveloysumm[3][MNGP];
	  double hveloxsumm[3][MNGP],hveloysumm[3][MNGP];
	  double enumberr[3][MNGP],hnumberr[3][MNGP];


	  double mobility[MNGP];

	  double vvxx=0,vvyy=0;
	  double xxf=0,yyf=0;

	  double ndensaver[MNGP],pdensaver[MNGP];
	  double ndenssum[MNGP],pdenssum[MNGP];
	  int numberndens[MNGP],numberpdens[MNGP];
	  double ndensaverr[3][MNGP],pdensaverr[3][MNGP];
	  double ndenssumm[3][MNGP],pdenssumm[3][MNGP];
	  int numberndenss[3][MNGP],numberpdenss[3][MNGP];

	  double eenergy[MNGP],henergy[MNGP];
	  double eenergysum[MNGP],henergysum[MNGP];

	  double eenergyy[3][MNGP],henergyy[3][MNGP];
	  double eenergysumm[3][MNGP],henergysumm[3][MNGP];

//     COMMON for statement parameters
//-----------------------------------------------------------------------
	  char		tempchar[255];
	  char		chartemp[255];
	  int		tempint;
	  double	tempdouble;
	  //material
	  int material;
	  //direction
	  int direction;
	  //surface_scattering
	  double tempdifprelectron,tempdifprhole;
	  //flag for quantum effect
	  int quantumeffectcbchar;
	  bool quantumeffectcbfl;
	  int quantumeffectqpchar;
	  bool quantumeffectqpfl;

	  bool quantumeffectfmfl;
	  bool quantumeffectbmfl;
	  bool quantumeffectqnfl;
	  int siregion;//for the calculation of quantum potential correction

	  //cutregion
	  int cutchar;
	  bool cutfl;
	  int xcutnumber,ycutnumber;
	  double positionx[MAXXCUTNUMBER],positiony[MAXYCUTNUMBER];
	  int ipcutx[MAXXCUTNUMBER],ipcuty[MAXYCUTNUMBER];

	  //the Minimum Doping concentration which will not take into acount Coulomb Scattering
	  double mindopforcoulomb;
	  //the Minimum electron energy which will allow taking into acount Coulomb Scattering
	  double minelectronenergy;



	  //surface scattering
	  //______surface scattering flag 
	  //______including 
	  //______surface roughness scattering
	  //______surface phonon scattering
	  //______surface impurity scattering
	  bool srsumfl;

	  //______xmin xmax ymin ymax of surface scattering region
	  double xminss1,xmaxss1,yminss1,ymaxss1;
	  double xminss2,xmaxss2,yminss2,ymaxss2;
	  //______max surface scattering rate
	  double surfacemax;
	  int tempsssfl;
	  bool sssfl;
	  bool ssregionfl;
	  int surfacenumber;
	  double surfaceposition1,surfaceposition2;
	  double srrv,sprv,sirv,ssrv;
	  //flag for the surface-roughness flag
	  int tempsrfl;
	  bool srfl;
	  //flag for the surface-phonon flag
	  int tempspfl;
	  bool spfl;
	  //flag for the surface-coulomb flag
	  int tempscfl;
	  bool scfl;
	  double ail,delta;

	  int tempballisticfl;
	  bool ballisticfl,channelfl;

	  int tempselfconsistentchar;
	  bool selfconsistentfl;

	  int tempbtbtchar;
	  bool btbtfl;

	  int tempmethodforgetisedatachar;
	  bool methodforgetisedatafl;

	  double RBBT[MNGP],RBBT0[MNGP],RBBT1[MNGP];
	  double RBBL[MNGP],RBBL0[MNGP],RBBL1[MNGP];
	  double Efn[MNGP],Efp[MNGP];
	  double Ei[MNGP];
	  double Ec[MNGP],Ev[MNGP];
	  double DE[MNGP];
	  double RBBDENS[MNGP],RBBTDENS[MNGP],RBBLDENS[MNGP];
	  double RBBT0DENS[MNGP],RBBT1DENS[MNGP];
	  double RBBL0DENS[MNGP],RBBL1DENS[MNGP];

	  double ifminx,ifmaxx,ifminy,ifmaxy;

	  double rr,RDENS[MNGP],NRDENS[3][MNGP],PRDENS[3][MNGP];

	  double Nc[MNGP],Nv[MNGP];
	  double dc[MNGP],dv[MNGP];
	  int nfpoint;
	  double fermi[2][2000];

	  int btbtTn,btbtTp;
	  int btbtLn,btbtLp;

	  double ebandsum[MNGP][4],hbandsum[MNGP][3];
	  double ebandratio[MNGP][4],hbandratio[MNGP][3];

	  bool firstrcfl;

	  int ifieldbuf[3][IVPP+1][MNPAR];
	  double dfieldbuf[3][DVPP][MNPAR];
	  int nparbufall[3],nparbufe[3],nparbufh[3];

	  bool	isymflag;

//    Maximum number of slaves under PVM
      const int MNSLAV=16;

//    number of the command in the input file currentl;y processed
      int linum;

//    Data path for band structure
      char cdirbs[100];

//    flags whether the initialize statement has been processed or not
      bool inifl;

//     enables CPU time consuming testing of MC simulation (for debugging)
      bool testfl;

//    enables 1D real space simulation (device is in x-direction)
      bool odxfl;

//    enables bulk simulation
      bool bulkfl;

//    enables parallel execution under PVM
      bool pvmfl;

//    true for the master process under PVM (false: slave)
      bool masterfl;

//     resume simulation after recoverable error
      bool resumefl;

//    workload balancing und PVM
      bool wlbfl;

//    simulate ning experiment (low energetic particles beneath the oxide are removed)
      bool ningfl;

//    quantum yield experiment (particles are distributed monoenergetically)
      bool qyfl;

//    load GaAs band strucutre
      bool gaasfl;

//    load silicon band structure
      bool sifl;

//    use SUN 48Bit random number generator
      bool rnlongfl;

//    lattice temperature
      double temp;

//    electric field for bulk simulation (yfield can be used in 1D simulation)
      double xfieldbulk, yfieldbulk, zfieldbulk;

//    inital energy of particles for quantum yield experiment
      double qyenergy;

//    flags whether the grid statement has been processed or not
      bool grifl;

//    flags whether the structure statement has been processed or not
      bool strfl;

//    position of the quadrant where the particles are removed in the case
//    of Ning's experiment
      int ijning;

//    flags whether the run statement has been processed or not
      bool runfl;

//    use potential from galene instead of the selfconsistent one
      bool potgalfl;

//    flags whether selfforce correction 8is used or not
      bool selffl;
//    enables injection into the oxide
      bool injfl;

//    enables just before scattering statisitcs
      bool jbscfl;

//    stop if electron or hole substrate current is converged
      bool elecisubfl, holeisubfl;

//    sample GHDM quantities
      bool ghdmfl;

//    sample energy distribution function with just before scattering statistics
//    (only homogeneous energy distribution functiuons)
      bool jbedffl;

//    quantum effects in the channel area
      bool quantfl;

//    gate depletion effect in the gatedep area
      bool gatedepfl;

//    stop if electron or hole velocity is converged
      bool elecghdmfl, holeghdmfl;

//    current time steps
      int cdt;

//    number of time steps
      int ndt;

//    current array of every dt
	  double currentarray[MNCONT][MNTIMES];
	  double averagecurrent[MNCONT];

	  double currentarrayy[3][MNCONT][MNTIMES];
	  double averagecurrentt[3][MNCONT];

	  double idtime[MNTIMES];
	  double didtime[MNTIMES];
	  double idtimeguiyi[MNTIMES];
	  double fhz[MNTIMES];
	  double idba;

	  int fhznumber;

//    for every outmod-th time step output data is printed onto the screen
      int outmod;

//    for stat time 
	  int stridestat;

//    for refresh time
	  int stridemr;

//	  for ADDJB time
	  int stridejb;

//    length of time step
      double dt;

//    start time for current stat
	  int starttime;

//    stat outmod
	  int statoutmod;

//    autocorrelation function CI(T)
	  double CI[MNTIMES];

//    Current Noise Spectral density
	  double SI[MNTIMES];

//    Electron Velocity Noise Spectral density
	  double SEV[MNREGION][MNTIMES];
	  double CEV[MNREGION][MNTIMES];
	  double ev[MNREGION][MNTIMES],evsum[MNREGION][MNTIMES],evnumber[MNREGION][MNTIMES];
	  double evaverage[MNREGION];

//    Hole Velocity Noise Spectral density
	  double SHV[MNREGION][MNTIMES];
	  double CHV[MNREGION][MNTIMES];
	  double hv[MNREGION][MNTIMES],hvsum[MNREGION][MNTIMES],hvnumber[MNREGION][MNTIMES];
	  double hvaverage[MNREGION];

//	  Carrier Density Noise Spectral density
	  double SCD[MNREGION][MNTIMES];
	  double CCD[MNREGION][MNTIMES];
	  double cd[MNREGION][MNTIMES];
	  double cdaverage[MNREGION];

	  int statnoiseregionnumber;
	  double xminnoise[MNREGION];
	  double xmaxnoise[MNREGION];
	  double yminnoise[MNREGION];
	  double ymaxnoise[MNREGION];

	  double xminnoisegrid[MNREGION];
	  double xmaxnoisegrid[MNREGION];
	  double yminnoisegrid[MNREGION];
	  double ymaxnoisegrid[MNREGION];

//    total simulated time
      double catime;

//    relative half width of the confidence interval
      double facstopii;

//    enables Multiple Refresh 
      bool mulfl[NPARTYP],mulfl0;// int mulfl2[NPARTYP+1];

//    enables MR in k-space only
      bool  enermulfl[NPARTYP];//int enermulfl2[NPARTYP];

//    doping profile has been loaded from Galene
      bool galdopfl;

//    electron (hole) density has been loaded from Galene
      bool galndfl, galpdfl;

//    task identifier for master (0) and slaves (>=1) under PVM
      int mtid[MNSLAV+1];

//    number of slaves under PVM
      int nslav;

//    TID of current process under PVM
      int mynum;

//    boundaries of the different processes in the real space grid under PVM
	  int ijbpvm[MNSLAV+1],ijepvm[MNSLAV+1];

//    first order low pass filter variable
      double rcdsv[MNSLAV+1];

//     filter constant of low pass for damping of workload balancing
      double cfilter;

//     >               xfieldbulk, yfieldbulk, zfieldbulk,
//     >               dt, ctime, rcdsv, cfilter,
//     >               qyenergy, facstopii,
//     >               inifl, testfl, runfl, potgalfl, selffl,
//     >               odxfl, bulkfl, pvmfl, masterfl, gaasfl, sifl,
//     >               grifl, strfl, resumefl, wlbfl, ningfl,
//     >               qyfl, injfl, jbscfl,
//     >               elecisubfl, holeisubfl, ghdmfl, jbedffl, 
//     >               quantfl, gatedepfl,
//     >               elecghdmfl, holeghdmfl, 
//     >               cdt, ndt, outmod, linum,
//     >               mulfl, enermulfl,
//     >               galdopfl, galndfl, galpdfl, rnlongfl,
//     >               mtid, nslav, mynum, ijbpvm, ijepvm, ijning,
//     >               cdirbs
