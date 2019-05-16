C***********************************************************************
C
      PROGRAM SHAW
C
C
C     Last Modified:  October 10, 2006
C
C     THIS IS THE DRIVER FOR THE SHAW MODEL (VERSION 2.4 Beta)
C
C     (PROVISIONS FOR SATURATED CONDITIONS)
C
C     DEVELOPED BY: GERALD N. FLERCHINGER
C                   USDA - ARS
C                   800 PARK BLVD., SUITE 105
C                   BOISE, IDAHO 83712
C                   TELEPHONE: 208-422-0716
C                   E-MAIL: gerald.flerchinger@ars.usda.gov
C
C***********************************************************************
C
C
C     LF    latent heat of fusion
C     LV    latent heat of vaporization
C     LS    latent heat of sublimation
C     G     acceleration of gravity (m/s^2)
C     UGAS  universal gas constant (J/mole/K)
C     RHOL  density of liquid water  (kg/m^3)
C     RHOA  density of air  (kg/m^3)
C     RHOI     density of ice   (kg/m^3)
C     RHOB(I)  bulk density of soil  (kg/m^3)
C     RHOM     density of the soil minerals (particle density)  (kg/m^3)
C     RHOOM    density of organic matter fraction of soil (kg/m^3)
C     RHOR(I)  density of residue in matted layer  (kg/m^3)
C     CL       specific heat capacity of liquid water  (J/kg/K)
C     CI       specific heat capacity of ice  (J/kg/K)
C     CM       specific heat capacity of mineral fraction of soil (J/kg/K)
C     COM      specific heat capacity of organic matter (J/kg/K)
C     CA       specific heat capacity of air  (J/kg/K)
C     CV       specific heat of water vapor  (J/kg/K)
C     CR       specific heat capacity of residue  (J/kg/K)
C     VONKRM   von Karman constant
C     VDIFF    vapor diffusivity in still air (m^2/s)
C     PRESUR   average atmospheric pressure at site  (N/m^2)
C     P0       standard atmospheric pressure at sea level  (101,300 Pa)
C     TKL      thermal conductivity of liquid (W/m/K)
C     TKI      thermal conductivity of ice (W/m/K)
C     TKA      thermal conductivity of air (W/m/K)
C     TKR      thermal conductivity of residue (W/m/K)
C     TKSI     thermal conductivity of silt fraction of soil (W/m/K)
C     TKCL     thermal conductivity of clay fraction of soil (W/m/K)
C     TKSA     thermal conductivity of sand  fraction of soil (W/m/K)
C     TKOM     thermal conductivity of organic matter fraction of soil (W/m/K)
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      REAL LF,LV,LS
C
C     SUNHOR(I)   radiation from sun upon a horizontal surface for hour i (W/m^2)
C     TMPDAY(I)   air temperature at hour i (C)
C     WINDAY(I)   average windspeed over hour i (m/s)
C     HUMDAY(I)   humidity at hour i
C     PRECIP(I)   total precipitation occurring at end of hour i (m)
C     SNODEN(I)   density of newly fallen snow for hour i of current day (g/cm^3)
C     SOITMP(I)   temperature at lower boundary for hour i (C)
C     VLCDAY(I)   volumetric moisture content (ice+water) at hour i (m/m)
C     SOILXT(I,J) soil water extracted from soil layer i for hour j based on user input (m)
      REAL SUNHOR(24),TMPDAY(24),WINDAY(24),HUMDAY(24),PRECIP(24),
     >     SNODEN(24),SOITMP(24),VLCDAY(24),SOILXT(50,24)
C
C     PLTWGT(J)   dry biomass of plant species j (kg/m^2)
C     PLTHGT(J)   plant height for plant species j (m)
C     PLTLAI(J)   total leaf area index of plant species j
C     ROOTDP(J)   rooting depth of plant species J (m)
C     DCHAR(J,I)  characteristic dimension (leaf width) of plant species j in canopy
C                 layer i (m)
C     TCCRIT(J)   minimum temperature required for plant j to transpire (C)
C     RSTOM0(J)   stomatal resistance of plant j with no water stress (s/m)
C     RSTEXP(J)   exponent for calculating stomatal resistance
C     PLEAF0(J)   critical leaf potential for plant j (m)
C     RLEAF0(J)   overall leaf conductivity for plant j (kg/m^3/s)
C     RROOT0(J)   overall root conductivity for plant j (kg/m^3/s)
C     PCANDT(J)   intercepted precipitation on plant species j at time t+dt (m)
C     CANALB(J)   albedo of plant j within the canopy (decimal)
      REAL PLTHGT(8),PLTWGT(8),PLTLAI(8),ROOTDP(8),DCHAR(8),
     >     TCCRIT(8),RSTOM0(8),RSTEXP(8),PLEAF0(8),
     >     RLEAF0(8),RROOT0(8),PCANDT(8),CANALB(8)

C
C     ZC(I)       depth of node i in the canopy from the top of the canopy (m)
C     TCDT(I)     temperature of canopy layer i at time t + dt (C)
C     VAPCDT(I)   vapor density at i-th canopy node at time t+dt  (kg/m^3)
C     WCANDT(I)   water content of dead plant material in canopy at time t+dt (kg/kg)
C     ZSP(I)      depth of node i in the snow from the top of the snowpack (m)
C     DZSP(I)     thickness of i-th snowpack node
C     RHOSP(I)    density of ice fraction of snowpack  (kg/m^3)
C     TSPDT(I)    temperature of snowpack layer at time t + dt (C)
C     DLWDT(I)    depth of liquid water for i-th snow pack node at end of time step
C     WLAG(I)     depth of water for i-th lag of snowpack outflow (m)
C     ZR(I)       depth of node i in the residue from the top of the residue (m)
C     RHOR(I)     density of residue in matted layer  (kg/m^3)
C     TRDT(I)     temperature of residue layer at time t + dt (C)
C     VAPRDT(I)   vapor density at i-th residue node at time t+dt (kg/m^3)
C     GMCDT(I)    gravimetric moisture content of i-th residue node at t+dt

      REAL ZC(11),TCDT(11),VAPCDT(11),WCANDT(10)
      REAL ZSP(100),DZSP(100),RHOSP(100),TSPDT(100),DLWDT(100),WLAG(11)
      REAL ZR(10),RHOR(10),TRDT(10),VAPRDT(10),GMCDT(10)
C
C     ZS(I)       depth of node i in the soil from the soil surface (m)
C     TSDT(I)     temperature of soil layer at time t + dt (C)
C     MATDT(I)    matric potential at t+dt (m)
C     CONCDT(J,I) concentration of j-th solute in soil solution at time t+dt (eq/kg)
C     VLCDT(I)    volumetric liquid content at t+dt (m/m)
C     VICDTS(I)   volumetric ice content at t+dt (m/m)
C     SALTDT(J,I) total amount of j-th solute in i-th soil node at time t+dt (eq/kg)
C     SLTDIF(J)   diffusion rate of the j-th solute in water (m^2/s)
C     ASALT(I)    coefficient relating solute diffusion in water to that in soil
C     DISPER(I)   hydrodynamic dispersion coefficient for solutes in i-th soil node
      REAL ZS(50),TSDT(50),MATDT(50),CONCDT(10,50),
     >     VLCDT(50),VICDT(50),SALTDT(10,50)
      REAL DGRADE(10),SLTDIF(10),ASALT(50),DISPER(50)
C

C
C     YRSTAR      year in which simulation is to start
C     HRSTAR      hour at which simulation is to start
C     YREND       year in which simulation is to end
C     HOUR        current hour
C     YEAR        current year
C     LEVEL(I)    array for desired level of output and when to switch levels
C     LVLOUT(I)   array containing flags for level of output for output files
C     LANGLE(J)   code for leaf orientation of plant j (0=verticle, 1=random)
C     ITYPE(J)    code for plant type for plant j (0=dead plant; 1=transpiring plant)
C     ICESDT(I)   1 if soil layer contains both ice and water at time t+dt
C     ICESPT(I)   1 if snow layer contains both ice and water at time t+dt
C     INITAL      flag for initializing the profile on first time into GOSHAW
      INTEGER YRSTAR,HRSTAR,YREND,HOUR,YEAR
      INTEGER LEVEL(6),LVLOUT(15)
      INTEGER LANGLE(8),ITYPE(8),ICESDT(50),ICESPT(100)
C
C     SET FLAG INDICATING FIRST TIME THROUGH SUBROUTINES FOR 
C     INITIALIZATION
      INITAL=0
C
C
C**** INPUT GENERAL INFORMATION, PROPERTIES, AND INITIAL CONDITIONS
C
C     MTSTEP     flag indicating time increment and format for input weather data
C     MZCINP     flag indicating user input of node spacing within the canopy
C     MPLTGRO    flag indicating whether plant growth curves are input
C     INPH2O     flag indicating whether soil water input is content or potential
C     MWATRXT    flag indicating if a sink term for water extraction by soil layer
C                (SOILXT) is supplied by the user
C     IVLCBC     flag indicating lower boundary condition for volumetric water
C                content (0 = specified input water content; 1 = unit gradient in
C                water potential)
C     ITMPBC     flag indicating the lower boundary condition for temperature
C     WCMAX      maximum allowable water content of dead plant material in canopy  (kg/kg)
C     WCAN(I)    water content of dead plant material in canopy at time t (kg/kg)
C     WCANDT(I)  water content of dead plant material in canopy at time t+dt (kg/kg)
C     GMC (I)    gravimetric moisture content of i-th residue node at time t (kg/kg)
C     GMCMAX (I) maximum allowable gravimetric moisture content of i-th residue node (kg/kg)
C     GMCDT (I)  gravimetric moisture content of i-th residue node at t+dt
C     COVER      fraction of ground surface covered by residue
C     RESCOF(I)  resistance to evaporation from residue elements for i-th node (s/m)
C     SLTDIF(J)  diffusion rate of the j-th solute in water (m^2/s)
C     ASALTS(I)  coefficient relating solute diffusion in water to that in soil
C     DISPER(I)  hydrodynamic dispersion coefficient for solutes in i-th soil node


C     ZERO       plane of zero displacement in wind profile calculations (m)
C     ZERSRF     plane of zero displacement in wind profile calculations above surface
C                residue or soil (m)
C     ZH         surface roughness parameter for heat transfer (m)
C     ZHSP       surface roughness parameter for heat transfer above snowpack (m)
C     ZHSRF      surface roughness parameter for heat transfer above surface residue or
C                soil snowcover (m)
C     ZM         surface roughness parameter for momentum transfer (m)
C     ZMSP       surface roughness parameter for momentum transfer above snowpack (m)
C     ZMSRF      surface roughness parameter for momentum transfer above surface
C                residue or soil (m)
C     HEIGHT     height of windspeed and temperature measurements (m)
C     NHRPDT     desired minimum number of hours per time step
C                Number of hours per time step
C     WDT        weighting factor for end of time step values
C     ALATUD       latitude of the simulation site (radians)
      CALL INPUT (NC,NSP,NR,NS,TOLER,LEVEL,MTSTEP,MZCINP,MPLTGRO,
     >  INPH2O,MWATRXT,LVLOUT,IVLCBC,ITMPBC,NPLANT,PLTHGT,PLTWGT,PLTLAI,
     >  ROOTDP,DCHAR,TCCRIT,RSTOM0,RSTEXP,PLEAF0,RLEAF0,RROOT0,
     >  CANALB,CANMA,CANMB,WCMAX,LANGLE,ITYPE,ZC,WCANDT,
     >  ZSP,DZSP,RHOSP,TSPDT,DLWDT,ICESPT,SNOTMP,
     >  GMCDT,GMCMAX,RLOAD,ZRTHIK,COVER,ALBRES,RESCOF,
     >  ZS,TSDT,VLCDT,VICDT,MATDT,CONCDT,ICESDT,SALTDT,ALBDRY,ALBEXP,
     >  DGRADE,SLTDIF,ASALT,DISPER,
     >  ZMSRF,ZHSRF,ZERSRF,ZMSP,ZHSP,HEIGHT,PONDMX,
     >  JSTART,YRSTAR,HRSTAR,JEND,YREND,NHRPDT,WDT,
     >  ALATUD,SLOPE,ASPECT,HRNOON)
C
C     JULIAN       current julian day
C     HOUR         current hour
C     YEAR         current year
C     MAXJUL       number of julian days in current year
      JULIAN=JSTART
      HOUR=HRSTAR
      YEAR=YRSTAR
      MAXJUL=365
      IF (MOD(YEAR,4) .EQ. 0)  MAXJUL=366
C
C**** DEFINE INITIAL BOUNDARY CONDITIONS FOR THE DAY
      IF (INPH2O.NE.1) THEN
C        INPUT SOIL MOISTURE IS WATER CONTENT
         VLCDAY(HOUR)=VLCDT(NS) + VICDT(NS)*RHOI/RHOL
        ELSE
C        INPUT SOIL MOISTURE IS MATRIC POTENTIAL
         VLCDAY(HOUR)=MATDT(NS)
      END IF
      SOITMP(HOUR)=TSDT(NS)
C
C     MWATRXT  flag indicating if a sink term for water extraction by soil layer
C             (SOILXT) is supplied by the user
      CALL DAYINP (JULIAN,YEAR,MAXJUL,HOUR,NHRPDT,MTSTEP,INITAL,
     >    MPLTGRO,MWATRXT,NS,ALATUD,HRNOON,SUNHOR,TMPDAY,WINDAY,HUMDAY,
     >    PRECIP,SNODEN,SOITMP,VLCDAY,SOILXT,NPLANT,PLTHGT,DCHAR,PLTWGT,
     >    PLTLAI,ROOTDP)
      CALL CLOUDY (CLOUDS,ALATUD,DECLIN,HAFDAY,SUNHOR,JULIAN,NHRPDT)
C
C**** BEGIN SIMULATION
      IF (LVLOUT(13) .NE. 0) WRITE (6,525)
C
C
C**** START OF A NEW HOUR - UPDATE VALUES ASSUMED CONSTANT OVER THE HOUR
  100 DT=NHRPDT*3600.
      HOUR=HOUR+NHRPDT
      IF (HOUR .GT. 24) THEN
C        START A NEW DAY
         HOUR=HOUR-24
         JULIAN=JULIAN + 1
         IF (JULIAN .GT. MAXJUL) THEN
            JULIAN=JULIAN-MAXJUL
            YEAR=YEAR + 1
            MAXJUL=365
            IF (MOD(YEAR,4) .EQ. 0)  MAXJUL=366
         END IF
         IF (YEAR .EQ. YREND  .AND.  JULIAN .GT. JEND)   STOP
         IF (YEAR .GT. YREND) STOP
         CALL DAYINP (JULIAN,YEAR,MAXJUL,0,NHRPDT,MTSTEP,1,
     >    MPLTGRO,MWATRXT,NS,ALATUD,HRNOON,SUNHOR,TMPDAY,WINDAY,HUMDAY,
     >    PRECIP,SNODEN,SOITMP,VLCDAY,SOILXT,NPLANT,PLTHGT,DCHAR,PLTWGT,
     >    PLTLAI,ROOTDP)
         CALL CLOUDY (CLOUDS,ALATUD,DECLIN,HAFDAY,SUNHOR,JULIAN,NHRPDT)
      END IF
C
      IF (JULIAN .EQ. LEVEL(2)  .AND.  HOUR .EQ. LEVEL(3)) THEN
C        CHANGE LEVEL OF OUTPUT
         LVL = LEVEL(1)
         LEVEL(1) = LEVEL(4)
         LEVEL(2) = LEVEL(5)
         LEVEL(3) = LEVEL(6)
         LEVEL(4) = LVL
         LEVEL(5) = JULIAN
         LEVEL(6) = HOUR
      END IF
C
C 
      CALL GOSHAW (JULIAN,HOUR,YEAR,NHRPDT,WDT,DT,INITAL,
     >  NC,NSP,NR,NS,TOLER,LEVEL,MZCINP,
     >  INPH2O,MWATRXT,LVLOUT,IVLCBC,ITMPBC,NPLANT,PLTHGT,PLTWGT,PLTLAI,
     >  ROOTDP,DCHAR,TCCRIT,RSTOM0,RSTEXP,PLEAF0,RLEAF0,RROOT0,PCANDT,
     >  CANALB,CANMA,CANMB,WCMAX,LANGLE,ITYPE,ZC,TCDT,VAPCDT,WCANDT,
     >  ZSP,DZSP,RHOSP,TSPDT,DLWDT,ICESPT,WLAG,STORE,SNOWEX,SNOTMP,
     >  ZR,RHOR,TRDT,VAPRDT,GMCDT,GMCMAX,RLOAD,ZRTHIK,COVER,ALBRES,
     >  RESCOF,DIRRES,
     >  ZS,TSDT,VLCDT,VICDT,MATDT,CONCDT,ICESDT,SALTDT,ALBDRY,ALBEXP,
     >  DGRADE,SLTDIF,ASALT,DISPER,        
     >  ZMSRF,ZHSRF,ZERSRF,ZMSP,ZHSP,HEIGHT,POND,PONDMX,
     >  ALATUD,SLOPE,ASPECT,HRNOON,CLOUDS,DECLIN,HAFDAY,
     >  SUNHOR(HOUR),TMPDAY(HOUR),WINDAY(HOUR),HUMDAY(HOUR),
     >  PRECIP(HOUR),SNODEN(HOUR),SOITMP(HOUR),VLCDAY(HOUR),
     >  SOILXT(1,HOUR))
C
      INITAL=1
C
      GO TO 100
C
  525 FORMAT(/20X,'Day   Hour   Year   Min Steps  Max Steps'/)
      END
C***********************************************************************
C
      SUBROUTINE GOSHAW (JULIAN,HOUR,YEAR,NHRPDT,WWDT,DTIME,INITAL,
     >  NC,NSP,NR,NS,TOLER,LEVEL,MZCINP,
     >  INPH2O,MWATRXT,LVLOUT,IVLCBC,ITMPBC,NPLANT,PLTHGT,PLTWGT,PLTLAI,
     >  ROOTDP,DCHAR,TCCRIT,RSTOM0,RSTEXP,PLEAF0,RLEAF0,RROOT0,PCANDT,
     >  CANALB,CANMA,CANMB,WCMAX,LANGLE,ITYPE,ZC,TCDT,VAPCDT,WCANDT,
     >  ZSP,DZSP,RHOSP,TSPDT,DLWDT,ICESPT,WLAG,STORE,SNOWEX,SNOTMP,
     >  ZR,RHOR,TRDT,VAPRDT,GMCDT,GMCMAX,RLOAD,ZRTHIK,COVER,ALBRES,
     >  RESCOF,DIRRES,
     >  ZS,TSDT,VLCDT,VICDT,MATDT,CONCDT,ICESDT,SALTDT,ALBDRY,ALBEXP,
     >  DGRADE,SLTDIF,ASALT,DISPER,        
     >  ZMSRF,ZHSRF,ZERSRF,ZMSP,ZHSP,HEIGHT,POND,PONDMX,
     >  ALATUD,SLOPE,ASPECT,HRNOON,CLOUDS,DECLIN,HAFDAY,
     >  SUNHOR,TMPDAY,WINDAY,HUMDAY,PRECIP,SNODEN,SOITMP,VLCDAY,SOILXT)
C
C***********************************************************************
      COMMON /TIMEWT/ WT,WDT,DT
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
C     SATK(I)     saturated conductivity of soil  (m/s)
C     NSALT       number of different types of solutes to be simulated
C     SALTKQ(J,I) adsorption characteristics of solute j for soil node i  (kg/kg)
C     VAPCOF(I)   coefficient for calculating vapor diffusivity in soil layer
C     VAPEXP(I)   exponent for calculating vapor diffusivity in soil layer
C     SAND(20)    fraction of sand in soil layer I

      COMMON /SLPARM/ B(50),ENTRY(50),RHOB(50),SAT(50),SATK(50),
     >                NSALT,SALTKQ(10,50),VAPCOF(50),VAPEXP(50),
     >                SAND(50),SILT(50),CLAY(50),OM(50)
      COMMON /MATRIX/ A1(150),B1(150),C1(150),D1(150),
     >                A2(150),B2(150),C2(150),D2(150)
      REAL LF,LV,LS
C
      REAL SOILXT(50)
C
      REAL DCHAR(8),CANALB(8),TCCRIT(8),RSTOM0(8),RSTEXP(8),PLEAF0(8),
     >     RLEAF0(8),RROOT0(8),ZC(11)
      REAL DGRADE(10),SLTDIF(10),ASALT(50),DISPER(50)
C
C     TLWCAN      net long-wave radiation for the canopy for time step (kJ)
C     LWSNOW(1)   net long-wave radiation absorbed by snow surface node (W/m^2)
C     LWSNOW(2)   net long-wave radiation absorbed by bottom snow node (W/m^2)
C     LWRES(I)    net long-wave radiation absorbed by i-th residue node (W/m^2)
C     LWSOIL      net long-wave radiation absorbed by soil surface (W/m^2)
C     SWCAN(I)    net short-wave radiation absorbed by i-th canopy node (W/m^2)
C     SWSNOW(I)   net short-wave radiation absorbed by i-th snow node (W/m^2)
C     SWRES(I)    net short-wave radiation absorbed by i-th residue node (W/m^2)
      REAL LWCAN(9,10),LWSNOW(2),LWRES(10),LWSOIL,SWCAN(9,10),
     >     SWSNOW(100),SWRES(10)

C
C     ZS(I)       depth of node i in the soil from the soil surface (m)
C     TS(I)       temperature of soil layer at time t (C)
C     TSDT(I)     temperature of soil layer at time t + dt (C)
C     MAT(I)      matric potential  (m)
C     MATDT(I)    matric potential at t+dt (m)
C     CONC(J,I)   concentration of j-th solute in soil solution at time t  (eq/kq)
C     CONCDT(J,I) concentration of j-th solute in soil solution at time t+dt (eq/kg)
C     VLC(I)      volumetric liquid content (m/m)
C     VICDT(I)    volumetric ice content at t+dt (m/m)
C     VIC(I)      volumetric ice content (m/m)
C     QSL(I)      average liquid flux for time step between node i and i+1 of the soil (m/s)
C     QSV(I)      vapor flux between node i and i+1 of the soil (kg/s/m^2)
C     XTRACT(I)   water extracted by roots for soil layer i (m/s)
C     US(I)       source term in for soil water balance (m/s)
C     SS(I)       heat souce for soil layer (W/m^2)
C     TK(I)       average thermal conductivity over the time step (W/m/K)
C     CS(I)       average heat capacity of the soil over time step  (J/kg/K)
C     SALT(J,I)   total amount of j-th solute in i-th soil node at time t (eq/kg)
C     SINK(J,I)   sink term for j-th solute in i-th soil node (eq/m^2/s)
C     TOTFLO(I)   total liquid water flow between soil layers i and i+1 for time step (m)
C     BTS(I)      value of TS(I) at beginning of user-specified time step (C)
C     BMAT(I)     value of MAT (I) at beginning of user-specified time step (m)
C     BCONC(J,I)  value of CONC (J,I) at beginning of user-specified time step (eq/kg)
      REAL ZS(50),TS(50),TSDT(50),MAT(50),MATDT(50),CONC(10,50),
     >     CONCDT(10,50),VLC(50),VLCDT(50),VIC(50),VICDT(50),QSL(50),
     >     QSV(50),XTRACT(50),US(50),SS(50),TK(50),CS(50),
     >     SALT(10,50),SALTDT(10,50),SINK(10,50),ROOTXT(50),TOTFLO(50)
      REAL BTS(50),BMAT(50),BCONC(10,50),BVLC(50),BVIC(50),BSALT(10,50)
      REAL RHOR(10),ZR(10)
C
C     VAPR(I)     vapor density at i-th residue node at time t  (kg/m^3)
C     GMC(I)      gravimetric moisture content of i-th residue node at time t (kg/kg)
C     UR(I)       source term in water balance for i-th residue node (kg/kg)
C     SR(I)       source term in energy balance for i-th residue node (W/m^2)
C     QVR(I)      vapor flux between nodes i and i+1 of the residue (kg/s/m^2)
      REAL TR(10),TRDT(10),VAPR(10),VAPRDT(10),GMC(10),GMCDT(10),
     >     UR(10),SR(10),QVR(10)
      REAL BTR(10),BVAPR(10),BGMC(10)

C
C     VAPC(I)     vapor density at i-th canopy node at time t  (kg/m^3)
C     WCAN(I)     water content of dead plant material in canopy at time t (kg/kg)
C     PCAN(J)     intercepted precipitation on plant species j at time t (m)
C     PCANDT(J)   intercepted precipitation on plant species j at time t+dt (m)
C     TRNSP(J)    total water transpired by plant j for the time step (kg)
C     UC
C     SC
C     QVC(I)      vapor flux between node i and i+1 of the canopy (kg/s/m^2)
C     BPCAN(I)    value of PCAN (I) at beginning of user-specified time step (m)
      REAL TC(11),TCDT(11),VAPC(11),VAPCDT(11),WCAN(10),WCANDT(10),
     >     PCAN(8),PCANDT(8),TRNSP(9),UC(10),SC(10),
     >     PLTHGT(8),PLTWGT(8),PLTLAI(8),ROOTDP(8),QVC(10)
      REAL BTC(11),BVAPC(11),BWCAN(10),BPCAN(8)
C
C     DZSP(I)     thickness of i-th snowpack node
C     DLW(I)      depth of liquid water for i-th snow pack node
C     QVSP(I)     vapor flux between nodes i and i+1 of the snowpack (kg/s/m^2)
C     TQVSP(I)    summation of vapor flux in snow pack over the time step (kg/m^2)
C     WLAG(I)     depth of water for i-th lag of snowpack outflow (m)

      REAL ZSP(100),DZSP(100),RHOSP(100)
      REAL TSP(100),TSPDT(100),DLW(100),DLWDT(100),USP(100),SSP(100),
     >     QVSP(100),TQVSP(100),WLAG(11)
      REAL BTSP(100),BDLW(100)
C
C     DELTA(I)    matrix solution for change in parameter for i-th node
      REAL DELTA(150),DELNRG(150),DELWTR(150)
C
C     HOUR         current hour
C     HRSTAR       hour at which simulation is to start
C     YEAR        current year
      INTEGER HOUR,HRSTRT,YEAR
      INTEGER LEVEL(6),LVLOUT(15)
      INTEGER LANGLE(8),ITYPE(8)
C
C     IBICES(I)   value of ICES(I) at beginning of user-specified time step
C     ICES (I)     1 if soil layer contains both ice and water at time t
      INTEGER IBICES(50),ICES(50),ICESDT(50),ICESP(100),ICESPT(100)
C
      SAVE NPRINT,MAXSTP,MINSTP
      DATA NPRINT/ 0/ MAXSTP/ 0/ MINSTP/ 0/
C
C
C--------------------初始化
      IF (INITAL .EQ. 0) THEN
C        FIRST TIME INTO SUBROUTINE FOR CURRENT PROFILE -- INITIALIZE
C        POORLY DEFINED STATE VARIABLES 
C
C        DEFINE MATRIC POTENTIAL OR WATER CONTENT OF SOIL AND DETERMINE
C        WHETHER SOIL IS FROZEN
         DO 10 I=1,NS
            IF (INPH2O.NE.1) THEN
C              INPUT SOIL MOISTURE IS WATER CONTENT
               CALL MATVL1 (I,MATDT(I),VLCDT(I),DUMMY)
              ELSE
C              INPUT SOIL MOISTURE IS MATRIC POTENTIAL
               CALL MATVL2 (I,MATDT(I),VLCDT(I),DUMMY)
               VICDT(I)=0.0
               ICESDT(I)=0
            END IF
            IF (TSDT(I) .LE. 0.0) CALL FROZEN (I,VLCDT,VICDT,MATDT,
     >                                 CONCDT,TSDT,SALTDT,ICESDT)
   10    CONTINUE
C
C        INITIALIZE PONDING AND SNOWPACK LAG AND STORAGE 
C     STORE        excess water in snowpack currently being attenuated (m)
C     POND         current depth of ponding (m)
         POND=0.0
         STORE=0.0
         SNOWEX=0.0
         DO 15 I=1,11
            WLAG(I)=0.0
   15    CONTINUE
         IF (NR .GT. 0) THEN
C           SET INITIAL VAPOR DENSITY AND TEMPERATURE OF RESIDUE EQUAL 
C           TO SOIL SURFACE NODE.
C
C
C
C     CONCDT (J,I) concentration of j-th solute in soil solution at time t+dt (eq/kg)
C     VAPDT        vapor density of surface node at time t+dt   (kg/m^3)
            TLCONC=0.0
            DO 20 J=1,NSALT
               TLCONC=TLCONC+CONCDT(J,1)
   20       CONTINUE
            TOTPOT=MATDT(1)-TLCONC*UGAS*(TSDT(1)+273.16)/G
            CALL VSLOPE (DUMMY,SATVDT,TSDT(1))
            VAPDT=SATVDT*EXP(.018*G/UGAS/(TSDT(1)+273.16)*TOTPOT)
            HUM=VAPDT/SATVDT
C
C           INITIALIZE CONDITIONS AND PARAMETERS FOR EACH LAYER
C
C     DIRRES       coefficient for direct radiation transmissivity through residue
            RESDEN=RLOAD/ZRTHIK
            IF (COVER .LT. 0.999) THEN
               DIRRES=-ALOG(1.-COVER)/RLOAD/10.
              ELSE
C              AVOID NUMERICAL OVERLOAD
               DIRRES = 6.9/RLOAD/10.
            END IF
C           INITIALIZE RESIDUE WATER CONTENT DEPENDING ON INPUT
            IF(GMCDT(1).LE.0.0)
     >         CALL RESHUM(2,HUM,DUMMY,GMCDT(1),TSDT(1))
            DO 25 I=1,NR
C     ZRTHIK       thickness of residue layer (m).
               ZR(I)=ZRTHIK*I/(NR+1)
               GMCDT(I)=GMCDT(1)
               RHOR(I)=RESDEN
               VAPRDT(I)=VAPDT
               TRDT(I)= TSDT(1)
   25       CONTINUE
            ZR(1)=0.0
            ZR(NR+1)=ZRTHIK
         END IF
C
         IF (NPLANT.NE.0) THEN
C           IF NODE SPACING NOT USER SPECIFIED, INITIALIZE NC
            IF (MZCINP.EQ.0) NC=1
C           INITIALIZE TEMPERATURE AND VAPOR DENSITY OF TOP LAYER
C           (REMAINDER WILL BE SET IN SUBROUTINE CANLAY)
            CALL VSLOPE (DUMMY,SATV,TMPDAY)
            TCDT(1)=TMPDAY
            VAPCDT(1)=SATV*HUMDAY
            IF (WCANDT(1) .LE. 0.0) THEN
C             INITIALIZE WATER CONTENT BASED ON ATMOSPHERIC HUMIDITY
              CALL CANHUM (2,HUMDAY,DUMMY,WCANDT(1),TCDT(1),CANMA,CANMB)
            END IF
C           INITIALIZE PRECIP INTERCEPTION FOR CANOPY
            DO 30 I=1,NPLANT
               PCANDT(I)=0.0
   30       CONTINUE
C
C           DEFINE LAYERING OF CANOPY AND ROOT DISTRIBUTION
C           CANOPY对灌层进行分层，并初始化温度和水分状态
            CALL CANOPY (NPLANT,NC,NSP,0,MZCINP,ITYPE,INITAL,
     >        PLTHGT,PLTWGT,PLTLAI,RLEAF0,ZSP,ZC,TC,TCDT,VAPC,VAPCDT,
     >        WCAN,WCANDT,WCMAX,PCANDT,CANMA,CANMB,TMPDAY,HUMDAY)
C           DETERMINE ROOT DISTRIBUTION (IF NOT DEFINED BY USER)
            IF (MZCINP.NE.1 .AND. NC.GT.0) 
     >         CALL RTDIST (NPLANT,ITYPE,NS,ZS,ROOTDP,RROOT0)
           ELSE
            NC=0
         END IF
C
         HRSTRT=HOUR-NHRPDT!//
C        INITIALIZE THE WATER BALANCE SUMMARY
C
C
C     TOTFLO (I)   total liquid water flow between soil layers i and i+1 for time step (m)
C     RUNOFF       runoff calculated for the time step (m)

         IF (LVLOUT(7) .NE. 0)
     >   CALL WBALNC(NPLANT,NC,NSP,NR,NS,LVLOUT(7),JULIAN,HRSTRT,YEAR,
     >   ITYPE,INITAL,ZC,WCAN,WCANDT,PCAN,PCANDT,VAPC,VAPCDT,RHOSP,DZSP,
     >   DLWDT,WLAG,STORE,ZR,GMC,GMCDT,VAPR,VAPRDT,RHOR,ZS,VLC,VLCDT,
     >   VIC,VICDT,TOTFLO,PRECIP,RUNOFF,POND,EVAP1,ETSUM)
C
C        PRINT OUT INITIAL CONDITIONS
         CALL OUTPUT (NPLANT,NC,NSP,NR,NS,LVLOUT,0,INPH2O,
     >    JULIAN,HRSTRT,YEAR,INITAL,ZC,TCDT,VAPCDT,WCANDT,ROOTXT,RHOSP,
     >    ZSP,TSPDT,DLWDT,ZR,TRDT,VAPRDT,GMCDT,ZS,TSDT,VLCDT,VICDT,
     >    MATDT,TOTFLO,CONCDT,SALTDT)
C
C        PRINT OUT INITIAL FROST AND SNOW DEPTH
         IF (LVLOUT(10) .NE. 0) CALL FROST (NSP,NS,LVLOUT(10),JULIAN,
     >      HRSTRT,YEAR,INITAL,ZSP,RHOSP,DZSP,DLWDT,WLAG,STORE,
     >      ZS,VLCDT,VICDT,ICESDT)
C
      ELSE 
        IF (NPLANT .GT. 0) THEN
C        DEFINE LAYERING OF CANOPY AND ROOT DISTRIBUTION
C     CANOPY对灌层进行分层，并初始化温度和水分状态
         CALL CANOPY (NPLANT,NC,NSP,0,MZCINP,ITYPE,1,
     >        PLTHGT,PLTWGT,PLTLAI,RLEAF0,ZSP,ZC,TC,TCDT,VAPC,VAPCDT,
     >        WCAN,WCANDT,WCMAX,PCANDT,CANMA,CANMB,TMPDAY,HUMDAY)
C        DETERMINE ROOT DISTRIBUTION (IF NOT DEFINED BY USER)
         IF (MZCINP.NE.1 .AND. NC.GT.0) 
     >        CALL RTDIST (NPLANT,ITYPE,NS,ZS,ROOTDP,RROOT0)
        ELSE
         NC=0
        END IF
C-----------------------------------------------------------------------
        IF (LEVEL(1) .GE. 1) THEN
C        PRINT CONDITIONS AT BEGINNING OF TIME STEP FOR DEBUGGING
         CALL OUTPUT (NPLANT,NC,NSP,NR,NS,LVLOUT,1,INPH2O,JULIAN,
     >   HOUR-NHRPDT,YEAR,1,ZC,TCDT,VAPCDT,WCANDT,ROOTXT,RHOSP,
     >   ZSP,TSPDT,DLWDT,ZR,TRDT,VAPRDT,GMCDT,ZS,TSDT,VLCDT,VICDT,
     >   MATDT,TOTFLO,CONCDT,SALTDT)
        END IF
C-----------------------------------------------------------------------
      END IF
C
C
C     WWDT   desired weighting factor (WDT) for end of time step values
      WDT=WWDT
      WT=1.-WDT
      DT=DTIME
      NTIMES=1
      NDT=1
      MAXNDT=NHRPDT*80
      MAXDBL=MAXNDT
C
C
C
C**** UPDATE BOUNDARY CONDITIONS FOR NEXT TIME STEP
      TA = TMPDAY
      TADT=TMPDAY
      HUM = HUMDAY
      HUMDT=HUMDAY
      WIND=WINDAY
C
C
C**** SET LOWER BOUNDARY CONDITIONS
      IF (ITMPBC .GT. 0) THEN
C....    ESTIMATE TEMPERATURE AT LOWER BOUNDARY
         CALL SOILTK (NS,TK,VLCDT,VICDT)
         CALL SOILHT (NS,CS,VLCDT,VICDT,TSDT,MATDT,CONCDT)
C        CALCULATED DAMPING DEPTH
         DAMPNG=SQRT(2.*TK(NS)/CS(NS)/1.99238E-07)
C        CALCULATE WEIGHTING COEFFICIENT BASED ON DAILY TIME STEP
         AA=(-0.00082+0.00983957*DAMPNG/(ZS(NS)-ZS(NS-1)))
     >       *(ZS(NS)/DAMPNG)**(-0.381266)
         IF (AA .LT. 0.0) AA=0.0
C        ADJUST FOR ACTUAL TIME STEP
         AA=AA*NHRPDT/24
         TSDT(NS)=(1.-AA)*TSDT(NS) + AA*TSDT(NS-1)
        ELSE
C        LOWER TEMPERATURE BOUNDARY SPECIFIED FROM INPUT FILE
         TSDT(NS)=SOITMP
      END IF
C
      IF (IVLCBC .GT. 0) THEN
C....    UNIT GRADIENT IS ASSUMED FOR WATER FLUX AT BOTTOM BOUNDARY;
C        SET MATRIC POTENTIAL FOR LAST NODE EQUAL TO SECOND TO LAST NODE
C        IF FREEZING FRONT IS IN NS-1, DO NOT CHANGE POTENTIAL AT BOTTOM 
         IF(VICDT(NS-1).LE.0. .OR. VICDT(NS).GT.0.)MATDT(NS)=MATDT(NS-1)
         VLC1=VLCDT(NS)+VICDT(NS)*RHOI/RHOL
         VICDT(NS)=0.0
         CALL MATVL2 (NS,MATDT(NS),VLCDT(NS),DUMMY)
         IF (TSDT(NS) .LT. 0.0) VICDT(NS)=(VLC1-VLCDT(NS))*RHOL/RHOI
         IF (VICDT(NS) .LT. 0.0) VICDT(NS)=0.0
        ELSE
C        INPUT WATER CONTENT SPECIFIED FOR WATER FLUX AT BOTTOM BOUNDARY
         IF (INPH2O.NE.1) THEN
C           INPUT SOIL MOISTURE IS WATER CONTENT
            VLCDT(NS)=VLCDAY
            CALL MATVL1 (NS,MATDT(NS),VLCDT(NS),DUMMY)
          ELSE
C           INPUT SOIL MOISTURE IS MATRIC POTENTIAL
            MATDT(NS)=VLCDAY
            CALL MATVL2 (NS,MATDT(NS),VLCDT(NS),DUMMY)
         END IF
         VICDT(NS)=0.0
      END IF
      ICESDT(NS)=0
      DO 105 J=1,NSALT
        CONCDT(J,NS)=SALTDT(J,NS)/(SALTKQ(J,NS)+VLCDT(NS)*RHOL/RHOB(NS))
  105 CONTINUE
      IF (TSDT(NS) .LE. 0.0) CALL FROZEN (NS,VLCDT,VICDT,MATDT,
     >                                   CONCDT,TSDT,SALTDT,ICESDT)
C
C
      IF (MWATRXT.EQ.1) THEN
C****    SOIL SINK TERM INPUT FOR EACH LAYER
C        CHECK IF SOIL LAYERS CAN SATIFY SINK TERM
         AVAILA=VLCDT(1)*(ZS(2)-ZS(1))/2
         IF (SOILXT(1)*DT.GT.AVAILA) SOILXT(1)=0.0
         AVAILA=VLCDT(NS)*(ZS(NS)-ZS(NS-1))/2
         IF (SOILXT(NS)*DT.GT.AVAILA) SOILXT(NS)=0.0
         DO 1105 I=2,NS-1
            AVAILA=VLCDT(I)*(ZS(I+1)-ZS(I-1))/2
            IF (SOILXT(I)*DT.GT.AVAILA) SOILXT(I)=0.0
 1105    CONTINUE
        ELSE
C        SET SOIL SINK TERM TO ZERO
         DO 1106 I=1,NS
            SOILXT(I)=0.0
 1106    CONTINUE
      END IF
C
C
C**** DEFINE THE ORDER IN WHICH THE MATERIALS ARE LAYERED
C     定义系统的ZM及ZH,当有植被时，根据植被就算；当物质时，根据雪面计算。
      CALL SYSTM (NC,NPLANT,NSP,ZC,ZSP,ZMSRF,ZHSRF,ZERSRF,ZMSP,
     >            ZHSP,HEIGHT,SUNHOR,TMPDAY,TCCRIT)
C
C
C**** CALCULATE THE SHORT-WAVE ENERGY BALANCE
C     确定冠层、腐殖质层、雪层和土壤层的短波辐射量，分别保存于swcan,swsnow,swres和swsoil
      CALL SWRBAL (NPLANT,NC,NSP,NR,LANGLE,SWCAN,SWSNOW,SWRES,
     > SWSOIL,SUNHOR,CANALB,ZSP,DZSP,RHOSP,ZR,RHOR,ALBRES,DIRRES,
     > ALBDRY,ALBEXP,VLCDT(1),ALATUD,SLOPE,ASPECT,HRNOON,HAFDAY,DECLIN,
     >  HOUR,NHRPDT)
C
C
C**** SAVE VALUES AT THE BEGINNING OF THE TIME STEP
      CALL UPDATE (NS,NR,NSP,NC,NPLANT,NSALT,IBICES,ICESDT,BTS,TSDT,
     > BMAT,MATDT,BCONC,CONCDT,BVLC,VLCDT,BVIC,VICDT,BSALT,SALTDT,
     > BTR,TRDT,BVAPR,VAPRDT,BGMC,GMCDT,BTC,TCDT,BVAPC,VAPCDT,
     > BWCAN,WCANDT,BPCAN,PCANDT,BTSP,TSPDT,BDLW,DLWDT,ICESP,ICESPT)
C
C**** UPDATE PARAMETERS FOR NEXT TIME STEP
  110 CALL UPDATE (NS,NR,NSP,NC,NPLANT,NSALT,ICES,ICESDT,TS,TSDT,
     > MAT,MATDT,CONC,CONCDT,VLC,VLCDT,VIC,VICDT,SALT,SALTDT,
     > TR,TRDT,VAPR,VAPRDT,GMC,GMCDT,TC,TCDT,VAPC,VAPCDT,WCAN,WCANDT,
     > PCAN,PCANDT,TSP,TSPDT,DLW,DLWDT,ICESP,ICESPT)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C**** START ITERATIVE PROCEDURE TO SOLVE ENERGY AND MOISTURE BALANCE
C
C     ITER         current number of iterations for energy and water balance

  120 ITER = 0
      ITRSLT = 0
  200 ITER = ITER + 1
      IF (WDT.LE.1.0) THEN
C        IF SATURATED OR EXTREMELY DRY CONDITIONS EXIST, 
C        SPECIFY FULLY IMPLICIT SOLUTION
         DO 100 I=NS-1,1,-1
            IF (MATDT(I).GT.ENTRY(I) .OR. VLCDT(I).LT.0.001) THEN
               WDT=1.0
               WT=0.0
               GO TO 101
            END IF
  100    CONTINUE          
  101    CONTINUE
      END IF
C
C
C     TS表示模拟前的土壤温度
C**** DETERMINE LONG-WAVE RADIATION BALANCE FOR EACH NODE
      CALL LWRBAL (NC,NSP,NR,NPLANT,TA,TADT,TC,TCDT,TSP,TSPDT,
     >         TR,TRDT,TS,TSDT,CLOUDS,LWCAN,LWSNOW,LWRES,LWSOIL)
C
C**** CALCULATE LONG-WAVE RAD. CONTRIBUTION TO ENERGY BALANCE MATRIX
      CALL LWRMAT (NC,NSP,NR,TC,TSP,TR,TS,ICESPT)
C
C**** SUM THE SOURCE-SINK TERMS FOR EACH NODE
      CALL SOURCE (NC,NSP,NR,NS,NSALT,NPLANT,UC,SC,USP,SSP,UR,SR,US,SS,
     >  SINK,SWCAN,SWSNOW,SWRES,SWSOIL,LWCAN,LWSNOW,LWRES,LWSOIL,
     >  SOILXT)
C
      N = 1
C
C**** DETERMINE THE BOUNDARY CONDITION FOR THE SURFACE MATERIAL
C
C     HTOLER       desired tolerance for heat flux calculation at surface
C     MATERL       number of materials above the material currently considering

      HTOLER=0.001
      MATERL=2

      IF (NC .GT. 0) THEN
C        CANOPY IS THE SURFACE MATERIAL
         CALL ATSTAB (NPLANT,NC,NSP,NR,TA,TADT,TC(1),TCDT(1),HUM,HUMDT,
     >     VAPC(1),VAPCDT(1),WIND,HEIGHT,HTOLER,HFLUX,VFLUX,ZC,ZR,RHOR,
     >     0,ITER)
C        GO THROUGH CALCULATION OF THE VAPOR DENSITY AT THE SOIL SURFACE
         GO TO 205
      END IF
C
      IF (NSP .GT. 0) THEN
C        SNOW IS THE SURFACE MATERIAL
         CALL VSLOPE (DUMMY,SATV,TSP(1))
         CALL VSLOPE (DUMMY,SATVDT,TSPDT(1))
C        CALCULATE THE SATURATED VAPOR DENSITY OVER ICE
         TMP = TSP(1) + 273.16
         TMPDT = TSPDT(1) + 273.16
         VAP = SATV*EXP(0.018/(UGAS*TMP)*LF*TSP(1)/TMP)
         VAPDT = SATVDT*EXP(0.018/(UGAS*TMPDT)*LF*TSPDT(1)/TMPDT)
         CALL ATSTAB(NPLANT,NC,NSP,NR,TA,TADT,TSP(1),TSPDT(1),HUM,HUMDT,
     >     VAP,VAPDT,WIND,HEIGHT,HTOLER,HFLUX,VFLUX,ZC,ZR,RHOR,
     >     ICESPT(1),ITER)
C        GO THROUGH CALCULATION OF THE VAPOR DENSITY AT THE SOIL SURFACE
         GO TO 205
      END IF
C
      IF (NR .GT. 0) THEN
C        RESIDUE IS THE SURFACE MATERIAL
         CALL ATSTAB (NPLANT,NC,NSP,NR,TA,TADT,TR(1),TRDT(1),HUM,HUMDT,
     >     VAPR(1),VAPRDT(1),WIND,HEIGHT,HTOLER,HFLUX,VFLUX,ZC,ZR,RHOR,
     >     0,ITER)
      END IF
C
C     SOIL IS THE SURFACE MATERIAL
C     ---- FIRST CALCULATE THE TOTAL WATER POTENTIAL, THEN THE HUMIDITY
C          MAY BE DETERMINED.  CALL VSLOPE TO OBTAIN THE SATURATED
C          VAPOR DENSITY.  USING HUMIDITY AND SATURATED VAPOR DENSITY,
C          THE VAPOR PRESSURE AT THE SOIL SURFACE MAY BE CALCULATED.
  205 TLCONC=0.0
      TLCNDT=0.0
      DO 210 J=1,NSALT
         TLCONC=TLCONC+CONC(J,1)
         TLCNDT=TLCNDT+CONCDT(J,1)
  210 CONTINUE
      TOTPOT=MAT(1)-TLCONC*UGAS*(TS(1)+273.16)/G
      TOTPDT=MATDT(1)-TLCNDT*UGAS*(TSDT(1)+273.16)/G
      CALL VSLOPE (DUMMY,SATV,TS(1))
      CALL VSLOPE (DUMMY,SATVDT,TSDT(1))
      VAP=SATV*EXP(.018*G/UGAS/(TS(1)+273.16)*TOTPOT)
      VAPDT=SATVDT*EXP(.018*G/UGAS/(TSDT(1)+273.16)*TOTPDT)
      IF (NC.EQ.0 .AND. NR.EQ.0 .AND. NSP.EQ.0)
     >   CALL ATSTAB (NPLANT,NC,NSP,NR,TA,TADT,TS(1),TSDT(1),HUM,HUMDT,
     >     VAP,VAPDT,WIND,HEIGHT,HTOLER,HFLUX,VFLUX,ZC,ZR,RHOR,
     >     ICESDT(1),ITER)
C
C
C**** DETERMINE ENERGY BALANCE FOR THE CANOPY LAYERS
      IF (NC .GT. 0) THEN
C        DEFINE VAPOR DENSITY FOR THE LOWER BOUNDARY OF CANOPY
         IF (NSP .GT. 0) THEN
C           CANOPY IS OVERLYING SNOWPACK
            CALL VSLOPE (DUMMY,SATV,TSP(1))
            CALL VSLOPE (DUMMY,SATVDT,TSPDT(1))
            TMP = TSP(1) + 273.16
            TMPDT = TSPDT(1) + 273.16
            VAPC(NC+1) = SATV*EXP(0.018/(UGAS*TMP)*LF*TSP(1)/TMP)
           VAPCDT(NC+1)=SATVDT*EXP(0.018/(UGAS*TMPDT)*LF*TSPDT(1)/TMPDT)
            TC(NC+1)=TSP(1)
            TCDT(NC+1)=TSPDT(1)
          ELSE
            IF (NR .GT. 0) THEN
C              CANOPY IS OVERLYING RESIDUE
               VAPC(NC+1) = VAPR(1)
               VAPCDT(NC+1) = VAPRDT(1)
               TC(NC+1) = TR(1)
               TCDT(NC+1) = TRDT(1)
             ELSE
C              CANOPY IS OVERLYING BARE SOIL
               VAPC(NC+1) = VAP
               VAPCDT(NC+1) = VAPDT
               TC(NC+1) = TS(1)
               TCDT(NC+1) = TSDT(1)
            END IF
         END IF
         CALL EBCAN (N,NPLANT,NC,NSP,NR,NS,ZC,TC,TCDT,VAPC,VAPCDT,
     >        WCAN,WCANDT,PCAN,PCANDT,MAT,MATDT,SWCAN,LWCAN,CANMA,CANMB,
     >        DCHAR,RSTOM0,RSTEXP,PLEAF0,RLEAF0,ITYPE,ITER)
         MATERL = MATERL + 1
      END IF
C
C**** DETERMINE ENERGY BALANCE FOR THE SNOW LAYERS
      IF (NSP .GT. 0) THEN
C        DEFINE VAPOR DENSITY AND TEMPERATURE FOR LOWER BOUNDARY OF 
C        SNOWPACK
         IF (NR .GT. 0) THEN
C           SNOW IS OVERLYING RESIDUE
            VAPSP = VAPR(1)
            VAPSPT = VAPRDT(1)
            TSP(NSP+1) = TR(1)
            TSPDT(NSP+1) = TRDT(1)
           ELSE
C           SNOW IS OVERLYING BARE SOIL
            VAPSP = VAP
            VAPSPT = VAPDT
            TSP(NSP+1) = TS(1)            
            TSPDT(NSP+1) = TSDT(1)
         END IF
         CALL EBSNOW (N,NSP,NR,ICESPT,TSP,TSPDT,DLW,DLWDT,RHOSP,ZSP,
     >                DZSP,QVSP,VAPSP,VAPSPT,SSP,ITER)
         MATERL = MATERL + 1
      END IF
C
C**** DETERMINE ENERGY BALANCE FOR THE RESIDUE LAYERS
      IF (NR .GT. 0) THEN
         VAPR(NR+1) = VAP
         VAPRDT(NR+1)=VAPDT
         TR(NR+1)=TS(1)
         TRDT(NR+1)=TSDT(1)
         CALL EBRES (N,NR,NSP,ZR,TR,TRDT,VAPR,VAPRDT,GMC,GMCDT,
     >               GMCMAX,RHOR,RESCOF,SR,UR,QVR,RHOSP,ITER)
         MATERL = MATERL + 1
      END IF
C
C**** SOLVE FOR ENERGY BALANCE OF SOIL
      CALL EBSOIL (N,NS,ZS,TS,TSDT,MAT,MATDT,CONC,CONCDT,
     >             VLC,VLCDT,VIC,VICDT,ICES,ICESDT,QSL,QSV,SS,ITER)
C
C-----------------------------------------------------------------------
      IF (LEVEL(1) .GE. 2) THEN
         WRITE (21,*) 'HFLUX (W/M2), VFLUX (KG/S): ', HFLUX,VFLUX
         WRITE (21,*) ' VAPOR FLUXES (KG/S)'
         WRITE (21,*) (QSV(K), K=1,NS-1)
         WRITE (21,*) ' LIQUID FLUXES (M/S)'
         WRITE (21,*) (QSL(K), K=1,NS-1)
         WRITE (21,*) ' JACOBIAN MATRIX'
         DO 215 K=1,N
            WRITE (21,*) A1(K),B1(K),C1(K),D1(K)
  215    CONTINUE
         WRITE (21,*)
         WRITE (21,*) ' VALUES FOR ENERGY BALANCE ITERATION ',ITER
      END IF
C-----------------------------------------------------------------------
C
C**** SOLVE THE ENERGY BALANCE MATRIX
      CALL TDMA (N,A1,B1,C1,D1,DELTA)
C
C**** SORT OUT THE SOLUTION OF THE MATRIX INTO THE PROPER MATERIALS
      MATERL=2
      N=1
      IEFLAG=0
C
      IF (NC .GT. 0) THEN
C**** CANOPY LAYERS
      DO 220 I=1,NC
         IF (ABS(DELTA(N)) .GT. 25.  .AND.  NDT .LT. MAXNDT) THEN
C           TOO LARGE OF A CHANGE -- CONVERGENCE WILL NOT LIKELY BE
C           BE MET WITH CURRENT TIME STEP.  CUT TIME STEP IN HALF
            IEFLAG=1
            ITER=11
            GO TO 350
         END IF
         IF (ABS(DELTA(N)) .GT. TOLER) IEFLAG = IEFLAG+1
         IF (ITER .GT. 3 .AND. DELNRG(N)*DELTA(N) .LT. 0) THEN
C           DELTA IS JUMPING BETWEEN NEG. AND POS. -- CUT DELTA IN HALF
C           TO SUPPRESS TENDENCY TO JUMP BACK AND FORTH AROUND SOLUTION
            DELTA(N)=DELTA(N)/2.
         END IF
         DELNRG(N)=DELTA(N)
         TCDT(I)=TCDT(I)-DELTA(N)
         IF (LEVEL(1) .GE. 2) WRITE (21,*) DELTA(N),TCDT(I),VAPCDT(I),
     >                                    WCANDT(I)
         N=N+1
  220 CONTINUE
      MATERL=MATERL+1
      END IF
C
      IF (NSP .GT. 0) THEN
C**** SNOW PACK LAYERS
      DO 230 I=1,NSP
         IF (ABS(DELTA(N)) .GT. 25.  .AND.  NDT .LT. MAXNDT) THEN
C           TOO LARGE OF A CHANGE -- CONVERGENCE WILL NOT LIKELY BE
C           BE MET WITH CURRENT TIME STEP.  CUT TIME STEP IN HALF
            IEFLAG=1
            ITER=11
            GO TO 350
         END IF
         IF (ICESPT(I) .EQ. 0) THEN
C           NO LIQUID WATER IN CURRENT LAYER AT END OF TIME STEP
            IF (ABS(DELTA(N)) .GT. TOLER) IEFLAG = IEFLAG + 1
            IF (ITER .GT. 3 .AND. DELNRG(N)*DELTA(N) .LT. 0) THEN
C              DELTA IS JUMPING BETWEEN NEG. AND POS. -- CUT DELTA IN HALF
C              TO SUPPRESS TENDENCY TO JUMP BACK AND FORTH AROUND SOLUTION
               DELTA(N)=DELTA(N)/2.
            END IF
            DELNRG(N)=DELTA(N)
            TSPDT(I) = TSPDT(I) - DELTA(N)
C           CHECK IF LAYER HAS BEGUN MELTING
            IF (TSPDT(I) .GT. 0) THEN
C              LAYER HAS GONE ABOVE 0 C - ADJUST WATER CONTENT IN LAYER
               ICESPT(I) = 1
CCCC           DLWDT(I) = RHOI*CI*DZSP(I)*TSPDT(I)/(RHOL*LF)
               TSPDT(I) = 0.0
            END IF
C
          ELSE
C           LAYER CONTAINS LIQUID WATER
C           CONVERT TOLERANCE FOR TEMPERATURE TO LIQUID EQUIVALENT
C           (DELTA LIQUID FRACTION) => 0.001*(DELTA TEMP.)
            IF (ABS(DELTA(N))/DZSP(I) .GT. 0.001*TOLER) THEN
C              IF DELTA < E-07*DLWDT, PRECISION OF COMPUTER IS EXCEEDED,
C              AND  ADDING DELTA TO DLWDT WILL NOT CHANGE DLWDT
               IF (ABS(DELTA(N)) .GT. DLWDT(I)*10E-07) IEFLAG=IEFLAG+1
            END IF
            IF (ITER .GT. 3 .AND. DELNRG(N)*DELTA(N) .LT. 0) THEN
C              DELTA IS JUMPING BETWEEN NEG. AND POS. -- CUT DELTA IN HALF
C              TO SUPPRESS TENDENCY TO JUMP BACK AND FORTH AROUND SOLUTION
               DELTA(N)=DELTA(N)/2.
            END IF
            DELNRG(N)=DELTA(N)
            DLWDT(I) = DLWDT(I) - DELTA(N)
C           CHECK IF ALL THE LIQUID HAS FROZEN
            IF (DLWDT(I) .LT. 0.0) THEN
C              LAYER HAS FROZEN COMPLETELY - ADJUST TEMPERATURE
               ICESPT(I) = 0
CCCC           TSPDT(I) = RHOL*LF*DLWDT(I)/(RHOI*CI)
               DLWDT(I) = 0.0
            END IF
         END IF
         IF (LEVEL(1) .GE. 2) WRITE (21,*) DELTA(N),TSPDT(I),DLWDT(I),
     >                                  DZSP(I)
         N=N+1
  230 CONTINUE
      MATERL=MATERL+1
      END IF
C
      IF (NR .GT. 0) THEN
C**** RESIDUE LAYERS
      DO 240 I=1,NR
         IF (ABS(DELTA(N)) .GT. 25.  .AND.  NDT .LT. MAXNDT) THEN
C           TOO LARGE OF A CHANGE -- CONVERGENCE WILL NOT LIKELY BE
C           BE MET WITH CURRENT TIME STEP.  CUT TIME STEP IN HALF
            IEFLAG=1
            ITER=11
            GO TO 350
         END IF
         IF (ABS(DELTA(N)) .GT. TOLER) IEFLAG = IEFLAG+1
         IF (ITER .GT. 3 .AND. DELNRG(N)*DELTA(N) .LT. 0) THEN
C           DELTA IS JUMPING BETWEEN NEG. AND POS. -- CUT DELTA IN HALF
C           TO SUPPRESS TENDENCY TO JUMP BACK AND FORTH AROUND SOLUTION
            DELTA(N)=DELTA(N)/2.
         END IF
         DELNRG(N)=DELTA(N)
         TRDT(I)=TRDT(I)-DELTA(N)
         IF (LEVEL(1) .GE. 2) WRITE (21,*) DELTA(N),TRDT(I),VAPRDT(I),
     >                                    GMCDT(I)
         N=N+1
  240 CONTINUE
      MATERL=MATERL+1
      END IF
C
C**** SOIL LAYERS
      DO 250 I=1,NS-1
         IF (ABS(DELTA(N)) .GT. 25.  .AND.  NDT .LT. MAXNDT) THEN
C           TOO LARGE OF A CHANGE -- CONVERGENCE WILL NOT LIKELY BE
C           BE MET WITH CURRENT TIME STEP.  CUT TIME STEP IN HALF
            IEFLAG=1
            ITER=11
            GO TO 350
         END IF
         IF (ABS(DELTA(N)) .GT. TOLER) IEFLAG = IEFLAG+1
         IF (ITER .GT. 3 .AND. DELNRG(N)*DELTA(N) .LT. 0) THEN
C           DELTA IS JUMPING BETWEEN NEG. AND POS. -- CUT DELTA IN HALF
C           TO SUPPRESS TENDENCY TO JUMP BACK AND FORTH AROUND SOLUTION
            DELTA(N)=DELTA(N)/2.
         END IF
         DELNRG(N)=DELTA(N)
         TSDT(I)=TSDT(I)-DELTA(N)
C
C        CHECK IF LAYER IS BELOW 0 C
         IF (TSDT(I) .LE. 0.0) THEN
C           ICE MAY BE PRESENT - CALL FROZEN TO DETERMINE IF ICE PRESENT
            ICE=ICESDT(I)
            CALL FROZEN (I,VLCDT,VICDT,MATDT,CONCDT,TSDT,SALTDT,ICESDT)
C           CHECK IF LAYER HAS CROSSED THE FREEZING POINT AND ADJUST THE
C           TEMPERATURE FOR LATENT HEAT IF SO
            IF (ICE .EQ. 0  .AND.  ICESDT(I) .EQ. 1)
     >      CALL ADJUST (I,VICDT,VLCDT,MATDT,CONCDT,TSDT,SALTDT,ICESDT)
C
           ELSE
C           NO ICE IS PRESENT
            IF (ICESDT(I) .EQ. 1) THEN
C              CONVERT ANY REMAINING ICE TO WATER
               VLCDT(I) = VLCDT(I) + VICDT(I)*RHOI/RHOL
               IF (VLCDT(I) .GT. SAT(I)) VLCDT(I)=SAT(I)
               VICDT(I) = 0.0
               ICESDT(I) = 0
               CALL MATVL1 (I,MATDT(I),VLCDT(I),DUMMY)
            END IF
         END IF
C
         IF (LEVEL(1) .GE. 2) WRITE (21,*) DELTA(N),TSDT(I),VLCDT(I),
     >                           VICDT(I),MATDT(I),CONCDT(1,I),ICESDT(I)
         N=N+1
  250 CONTINUE
C
      IF (NC .GT. 0) THEN
C        DEFINE TEMPERATURE OF BOTTOM BOUNDARY OF CANOPY
         IF (NSP .GT. 0) THEN
C           SNOW IS UNDERLYING CANOPY
            TCDT(NC+1) = TSPDT(1)
           ELSE
            IF (NR .GT. 0) THEN
C              RESIDUE IS UNDERLYING CANOPY
               TCDT(NC+1) = TRDT(1)
              ELSE
C              SOIL IS UNDERLYING CANOPY
               TCDT(NC+1) = TSDT(1)
            END IF
         END IF
      END IF
C
      IF (NSP .GT. 0) THEN
C        DEFINE TEMPERATURE OF BOTTOM BOUNDARY OF SNOWPACK
         IF (NR .GT. 0) THEN
C           SNOW IS OVERLYING RESIDUE
            TSPDT(NSP+1) = TRDT(1)
           ELSE
C           SNOW IS OVERLYING BARE SOIL
            TSPDT(NSP+1) = TSDT(1)
         END IF
      END IF
C
C     DEFINE TEMPERATURE OF BOTTOM BOUNDARY OF RESIDUE
      IF (NR .GT. 0) TRDT(NR+1) = TSDT(1)
C
C-----------------------------------------------------------------------
      IF (LEVEL(1) .GE. 1) THEN
         WRITE (21,*) ' ENERGY BALANCE AT ITER = ',ITER,JULIAN,HOUR,
     >                NTIMES,NDT
         CALL OUTPUT (NPLANT,NC,NSP,NR,NS,LVLOUT,1,INPH2O,
     >    JULIAN,HOUR,YEAR,1,ZC,TCDT,VAPCDT,WCANDT,ROOTXT,RHOSP,
     >    ZSP,TSPDT,DLWDT,ZR,TRDT,VAPRDT,GMCDT,ZS,TSDT,VLCDT,VICDT,
     >    MATDT,TOTFLO,CONCDT,SALTDT)
      END IF
C-----------------------------------------------------------------------
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C**** BEGIN CALCULATIONS FOR THE MOISTURE BALANCE OF THE SYSTEM
      N = 1
C
C**** SOLVE WATER BALANCE FOR CANOPY LAYERS
      IF (NC .GT. 0) THEN
C     QVC (I)      vapor flux between node i and i+1 of the canopy (kg/s/m^2)
C     TRNSP (J)    total water transpired by plant j for the time step (kg)
         CALL WBCAN(N,NPLANT,NC,NSP,NR,NS,ZC,TC,TCDT,VAPC,VAPCDT,WCAN,
     >     WCANDT,PCAN,PCANDT,QVC,TRNSP,MAT,MATDT,XTRACT,SWCAN,LWCAN,
     >     CANMA,CANMB,DCHAR,RSTOM0,RSTEXP,PLEAF0,RLEAF0,ITYPE,
     >     ICESDT(1),ITER)
       ELSE
C        SET ROOT EXTRACTION AND PLANT TRANSPIRATION TO ZERO
         DO 300 I=1,NS
            XTRACT(I)=0.0
  300    CONTINUE
C        TRNSP(NPLANT+1) IS TOTAL TRANSPIRATION FOR ALL PLANTS
         DO 302 I=1,NPLANT+1
            TRNSP(I)=0.0
  302    CONTINUE
      END IF
C
C**** SET BOUNDARY CONDITIONS IF THERE IS A SNOWPACK PRESENT
C     (SNOWPACK IS NOT PART OF THE WATER BALANCE MATRIX - THE WATER
C     BALANCE FOR THE SNOWPACK IS DONE AT THE END OF THE HOUR)
      IF (NSP .GT. 0.0) THEN
         A2(N) = 0.0
         CALL SNOWBC (N,NSP,NR,ZSP,QVSP,VAPSPT,TSDT,TSPDT,ICESDT(1))
      END IF
C
C**** SOLVE FOR WATER BALANCE OF THE RESIDUE LAYERS
      IF (NR .GT. 0) 
     >   CALL WBRES (N,NR,NSP,ZR,TR,TRDT,VAPR,VAPRDT,GMC,GMCDT,RHOR,
     >               QVR,RHOSP,ICESDT(1),ITER)
C
C**** SOLVE FOR THE WATER BALANCE OF SOIL LAYERS
      CALL WBSOIL (N,NS,ZS,TS,TSDT,MAT,MATDT,VLC,VLCDT,VIC,VICDT,
     >                   CONC,CONCDT,ICESDT,QSL,QSV,XTRACT,SEEP,US,ITER)
C
C-----------------------------------------------------------------------
      IF (LEVEL(1) .GE. 2) THEN
         WRITE (21,*) 'HFLUX (W/M2), VFLUX (KG/S): ', HFLUX,VFLUX
         WRITE (21,*) ' VAPOR FLUXES (KG/S)'
         WRITE (21,*) (QSV(K), K=1,NS-1)
         WRITE (21,*) ' LIQUID FLUXES (M/S)'
         WRITE (21,*) (QSL(K), K=1,NS-1)
         WRITE (21,*) ' JACOBIAN MATRIX'
         DO 305 I=1,N
            WRITE (21,*) A2(I),B2(I),C2(I),D2(I)
  305    CONTINUE
         WRITE (21,*)
         WRITE (21,*) ' VALUES AT WATER BALANCE ITERATION ', ITER
      END IF
C-----------------------------------------------------------------------
C
C**** SOLVE THE WATER BALANCE MATRIX
      CALL TDMA (N,A2,B2,C2,D2,DELTA)
C
C**** SORT OUT THE SOLUTION OF THE MATRIX INTO THE PROPER MATERIALS
      MATERL=2
      N=1
      IWFLAG=0
C
      IF (NC .GT. 0) THEN
C**** CANOPY LAYERS
      DO 310 I=1,NC
         IF (ABS(DELTA(N)/VAPCDT(I)) .GT. TOLER) IWFLAG = IWFLAG+1
         IF (ITER .GT. 3 .AND. DELWTR(N)*DELTA(N) .LT. 0) THEN
C           DELTA IS JUMPING BETWEEN NEG. AND POS. -- CUT DELTA IN HALF
C           TO SUPPRESS TENDENCY TO JUMP BACK AND FORTH AROUND SOLUTION
            DELTA(N)=DELTA(N)/2.
         END IF
         DELWTR(N)=DELTA(N)
         VAPCDT(I)=VAPCDT(I)-DELTA(N)
         IF (LEVEL(1) .GE. 2) WRITE (21,*) DELTA(N),VAPCDT(I),TCDT(I),
     >                                    WCANDT(I)
         N=N+1
  310 CONTINUE
      MATERL=MATERL+1
      END IF
C
      IF (NSP .GT. 0) THEN
C**** SNOW PACK LAYERS
C     SNOW IS NOT PART OF WATER BALANCE SOLUTION
      MATERL=MATERL+1
      END IF
C
      IF (NR .GT. 0) THEN
C**** RESIDUE LAYERS
      DO 320 I=1,NR
         IF (ABS(DELTA(N)/VAPRDT(I)) .GT. TOLER) IWFLAG = IWFLAG+1
         IF (ITER .GT. 3 .AND. DELWTR(N)*DELTA(N) .LT. 0) THEN
C           DELTA IS JUMPING BETWEEN NEG. AND POS. -- CUT DELTA IN HALF
C           TO SUPPRESS TENDENCY TO JUMP BACK AND FORTH AROUND SOLUTION
            DELTA(N)=DELTA(N)/2.
         END IF
         DELWTR(N)=DELTA(N)
         VAPRDT(I)=VAPRDT(I)-DELTA(N)
         IF (LEVEL(1) .GE. 2) WRITE (21,*) DELTA(N),VAPRDT(I),TRDT(I),
     >                                    GMCDT(I)
         N=N+1
  320 CONTINUE
      MATERL=MATERL+1
      END IF
C
C**** SOIL LAYERS
      DO 330 I=1,NS-1
         IF (ICESDT(I) .EQ. 1) THEN
            IF (ABS(DELTA(N)) .GT. 1.0  .AND.  NDT .LT. MAXNDT) THEN
C              TOO LARGE OF A CHANGE -- CONVERGENCE WILL NOT LIKELY BE
C              BE MET WITH CURRENT TIME STEP.  CUT TIME STEP IN HALF
               IWFLAG=1
               ITER=11
               GO TO 350
            END IF
            IF (ABS(DELTA(N)) .GT. TOLER) IWFLAG = IWFLAG+1
            IF (ITER .GT. 3 .AND. DELWTR(N)*DELTA(N) .LT. 0) THEN
C              DELTA IS JUMPING BETWEEN NEG. AND POS.--CUT DELTA IN HALF
C              TO SUPPRESS TENDENCY TO JUMP AROUND SOLUTION
               DELTA(N)=DELTA(N)/2.
            END IF
            DELWTR(N)=DELTA(N)
            VICDT(I)=VICDT(I)-DELTA(N)
            IF (VICDT(I) .LT. 0.0) THEN
C              CHANGE IN ICE IS GREATER THAN ICE CONTENT -- ADJUST WATER
C              CONTENT FOR THE DIFFERENCE, IF NOT GREATER THAN VLCDT(I)
               IF ((VLCDT(I)+VICDT(I)*RHOI/RHOL) .GT. VLCDT(I)/2.) THEN
                  VLCDT(I) = VLCDT(I) + VICDT(I)*RHOI/RHOL
                  CALL MATVL1 (I,MATDT(I),VLCDT(I),DUMMY)
                 ELSE
                  VLCDT(I) = VLCDT(I)/2.
                  CALL MATVL1 (I,MATDT(I),VLCDT(I),DUMMY)
               END IF
               VICDT(I)=0.0
               ICESDT(I)=0
            END IF
          ELSE
C
C           SAVE VALUE TO COMPARE RELATIVE CHANGE IN MATRIC POTENTIAL;
C           IF MATRIC POTENTIAL IS NEAR ZERO, SET CHECK VALUE TO 1.0
            CHKMAT=MATDT(I)
            IF (ABS(CHKMAT).LT.1.0) CHKMAT=1.0
C
            IF (ABS(DELTA(N)/CHKMAT).GT.100. .AND. NDT.LT.MAXNDT) THEN
C              TOO LARGE OF A CHANGE -- CONVERGENCE WILL NOT LIKELY BE
C              BE MET WITH CURRENT TIME STEP.  CUT TIME STEP IN HALF
               IWFLAG=1
               ITER=11
               GO TO 350
            END IF
            IF (ABS(DELTA(N)/CHKMAT) .GT. TOLER) IWFLAG = IWFLAG+1
            IF (ITER .GT. 3 .AND. DELWTR(N)*DELTA(N) .LT. 0) THEN
C              DELTA IS JUMPING BETWEEN NEG. AND POS.--CUT DELTA IN HALF
C              TO SUPPRESS TENDENCY TO JUMP AROUND SOLUTION
               DELTA(N)=DELTA(N)/2.
            END IF
            DELWTR(N)=DELTA(N)
            MATDT(I)=MATDT(I)-DELTA(N)
         END IF
         CALL MATVL2 (I,MATDT(I),VLCDT(I),DLDM)
         IF (LEVEL(1) .GE. 2) WRITE (21,*) DELTA(N),MATDT(I),VLCDT(I),
     >                                    VICDT(I)
         N=N+1
  330 CONTINUE
C     LIMIT WATER POTENTIAL AT SURFACE TO FREE WATER CONDITION 
C     TO ALLOW SEEPAGE TO OCCUR
      IF (MATDT(1).GT.0.0) THEN
C       CHANGE ANY SUBSEQUENT SATURATED NODES BY SAME AMOUNT
        DO 335 I=2,NS
      	  IF (MATDT(I).GT.0.0) THEN
      	    MATDT(I) = MATDT(I) - MATDT(1)
            CALL MATVL2 (I,MATDT(I),VLCDT(I),DUMMY)
      	   ELSE
C           DONE WITH SATURATED NODES
      	    GO TO 340
      	  END IF
  335   CONTINUE
  340   MATDT(1)=0.0
        CALL MATVL2 (1,MATDT(1),VLCDT(1),DUMMY)
      END IF
C
C-----------------------------------------------------------------------
      IF (LEVEL(1) .GE. 1) THEN
         WRITE (21,*) ' WATER BALANCE AT ITERWB = ',ITER,JULIAN,HOUR,
     >                NTIMES,NDT
         CALL OUTPUT (NPLANT,NC,NSP,NR,NS,LVLOUT,1,INPH2O,
     >    JULIAN,HOUR,YEAR,1,ZC,TCDT,VAPCDT,WCANDT,ROOTXT,RHOSP,
     >    ZSP,TSPDT,DLWDT,ZR,TRDT,VAPRDT,GMCDT,ZS,TSDT,VLCDT,VICDT,
     >    MATDT,TOTFLO,CONCDT,SALTDT)
         WRITE (21,*) ' ITERATION =',ITER,'ENERGY FLAGS =',IEFLAG,
     >               'WATER FLAGS =',IWFLAG
      END IF
C-----------------------------------------------------------------------
C
  350 IF (IWFLAG .GT. 0 .OR. IEFLAG .GT. 0) THEN
C        CONVERGENCE HAS NOT BEEN MET - IF ITERATIONS ARE UNDER 10, GO
C        BACK AND START NEXT ITERATION
         IF (ITER .LE. 10) GO TO 200
C
C        HAVING PROBLEMS REACHING CONVERGENCE WITH CURRENT TIME STEP
         IF (NDT .LT. MAXNDT) THEN
C           CUT TIME STEP IN HALF AND TRY AGAIN
            NDT=NDT*2
            NTIMES=NTIMES*2 - 1
            DT=DTIME/NDT
C           REDEFINE END-OF-TIME-STEP VALUES
            CALL BACKUP (NS,NR,NSP,NC,NPLANT,NSALT,ICES,ICESDT,TS,TSDT,
     >           MAT,MATDT,CONC,CONCDT,VLC,VLCDT,VIC,VICDT,SALT,SALTDT,
     >           TR,TRDT,VAPR,VAPRDT,GMC,GMCDT,TC,TCDT,VAPC,VAPCDT,WCAN,
     >           WCANDT,PCAN,PCANDT,TSP,TSPDT,DLW,DLWDT,ICESP,ICESPT)
C           GO BACK AND CALCULATE BOUNDARY CONDITIONS FOR HALF TIME-STEP
            GO TO 120
         END IF
C
C        NDT > MAXNDT  --> INDICATE CONVERGENCE PROBLEMS AND CONTINUE
         IF (LVLOUT(13).NE.0) WRITE (6,535) JULIAN,HOUR,YEAR,NTIMES,NDT
         WRITE (21,*) ' CONVERGENCE PROBLEMS AT : ',
     >            JULIAN,HOUR,YEAR,NTIMES
         IF (IEFLAG.GT.0) WRITE (21,*)
     >		' ENERGY BALANCE WILL NOT CONVERGE'
         IF (IWFLAG.GT.0) WRITE (21,*)
     >		' WATER BALANCE WILL NOT CONVERGE'
         ITER = 0
         IF (LEVEL(1) .EQ. 1 .OR. LEVEL(1) .EQ. 2) STOP
      END IF
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C**** CALCULATE SOLUTE CONCENTRATIONS AT THE END OF THE TIME STEP
      IF (NSALT .LE. 0) GO TO 400
      CALL SOLUTE (NS,ZS,CONC,CONCDT,SALT,SALTDT,VLC,VLCDT,TS,
     >             TSDT,QSL,XTRACT,SINK,DGRADE,SLTDIF,ASALT,DISPER)
C
C-----------------------------------------------------------------------
      IF (LEVEL(1) .GE. 1) THEN
         WRITE(21,522) (HOUR,ZS(I),TSDT(I),VLCDT(I),VICDT(I),
     >                MATDT(I),CONCDT(1,I),TA,TADT,SALTDT(1,I), I=1,NS)
      END IF
C-----------------------------------------------------------------------
C
      IF (ITER .GT. 1) THEN
C        IF IT TOOK MORE THAN ONE ITERATION FOR THE ENERGY AND WATER
C        BALANCE, END-OF-TIME-STEP VALUES HAVE CHANGED ENOUGH THAT THE
C        SOLUTE CONDITIONS ASSUMED WERE NOT ACCURATE.  GO BACK THROUGH
C        ENERGY AND WATER BALANCE AGAIN
         ITRSLT=ITRSLT+1
         IF (ITRSLT .LE. 4) THEN
            ITER = 0
            GO TO 200
         END IF
C
C        CHANGES IN SOLUTE BALANCE ARE CAUSING SUCCESSIVE SOLUTIONS TO
C        THE WATER-ENERGY BALANCE FOR THIS TIME STEP TO DIFFER TOO MUCH
         IF (NDT .LT. MAXNDT) THEN
C           CUT TIME STEP IN HALF AND TRY AGAIN
            NDT=NDT*2
            NTIMES=NTIMES*2 - 1
            DT=DTIME/NDT
C           REDEFINE END-OF-TIME-STEP VALUES
            CALL BACKUP (NS,NR,NSP,NC,NPLANT,NSALT,ICES,ICESDT,TS,TSDT,
     >           MAT,MATDT,CONC,CONCDT,VLC,VLCDT,VIC,VICDT,SALT,SALTDT,
     >           TR,TRDT,VAPR,VAPRDT,GMC,GMCDT,TC,TCDT,VAPC,VAPCDT,WCAN,
     >           WCANDT,PCAN,PCANDT,TSP,TSPDT,DLW,DLWDT,ICESP,ICESPT)
C           GO BACK AND CALCULATE BOUNDARY CONDITIONS FOR HALF TIME-STEP
            GO TO 120
         END IF
C
C        NDT > MAXNDT  ---> INDICATE CONVERGENCE PROBLEMS AND CONTINUE
         IF (LVLOUT(13).NE.0) WRITE (6,535) JULIAN,HOUR,YEAR,NTIMES,NDT
         WRITE (21,*) ' CONVERGENCE PROBLEMS AT : ',
     >                JULIAN,HOUR,YEAR,NTIMES
         WRITE (21,*)' SOLUTES CAUSING CONVERGENCE PROBLEMS'
         IF (LEVEL(1) .EQ. 1 .OR. LEVEL(1) .EQ. 2) STOP
      END IF
      ITRSLT = 0
  400 CONTINUE
C
C**** END OF ITERATION FOR TIME STEP ***********************************
C
      IF (MINSTP.EQ.0) THEN
         MAXSTP=NDT
         MINSTP=NDT
      END IF
      IF (MAXSTP .LT. NDT) MAXSTP=NDT
      IF (MINSTP .GT. NDT) MINSTP=NDT
C
C     SUM THE NECESSARY FLUXES OCCURRING OVER THE TIME STEP
      CALL SUMDT (NC,NPLANT,NSP,NR,NS,HFLUX,VFLUX,LWCAN,LWSNOW,
     >     LWRES,LWSOIL,SWCAN,SWSNOW,SWRES,SWSOIL,QVC,QVR,QVSP,QSL,QSV,
     >     TRNSP,XTRACT,SEEP,TSWSNO,TLWSNO,TSWCAN,TLWCAN,TSWRES,TLWRES,
     >     TSWSOI,TLWSOI,THFLUX,EVAP1,ETSUM,TSEEP,ROOTXT,TOTFLO,NTIMES,
     >     TOPSNO,TQVSP)

C
      IF (NTIMES .NE. NDT) THEN
C        NOT REACHED END OF TIME STEP -- GO BACK THROUGH SUB-TIME-STEP
         IF (MOD(NTIMES,2).EQ.0 .AND. (IEFLAG+IWFLAG).EQ.0) THEN
C           CHECK IF IT'S WORTH IT TO TRY TO DOUBLE THE TIME STEP
            IF (NDT .LE. MAXDBL) THEN
               IF (NDT .LT. MAXDBL) MAXTRY=0
               IF (NDT .LE. 8) MAXTRY=MAXTRY+2
               MAXTRY=MAXTRY+1
               MAXDBL=NDT
            END IF
            IF (NDT.GT.MAXDBL .OR. MAXTRY.LE.3 .OR. NDT.GE.MAXNDT) THEN
C              ATTEMPT TO INCREASE LENGTH OF SUB-TIME-STEP
               NDT=NDT/2
               NTIMES=NTIMES/2
               DT=DTIME/NDT
            END IF
         END IF
         NTIMES=NTIMES + 1
         GO TO 110
        ELSE
C        END OF TIME STEPS - INFILTRATE RAIN, DETERMINE OUTFLOW FROM 
C        SNOW,AND CALCULATE WATER BALANCE.  SET TIME STEP TO ONE HOUR 
C        FOR PRECIP CALULATIONS
         DT = DTIME
         RAIN=PRECIP
         TADT=TMPDAY
         HUMDT=HUMDAY
         TA=TMPDAY
         HUM=HUMDAY
         NSPLST=NSP
         CALL PRECP(NPLANT,NC,NSP,NR,NS,TA,TADT,HUM,HUMDT,
     >   ZC,LANGLE,ITYPE,WCANDT,WCMAX,PCANDT,
     >   ZSP,DZSP,TSPDT,BDLW,DLWDT,RHOSP,TQVSP,TOPSNO,WLAG,STORE,SNOWEX,
     >   ZR,TRDT,GMCDT,GMCMAX,RHOR,
     >   ZS,TSDT,VLCDT,VICDT,MATDT,TOTFLO,SALTDT,CONCDT,ICESPT,ICESDT,
     >   RAIN,POND,RUNOFF,EVAP1,PONDMX,SNOTMP,SNODEN,DIRRES,
     >   SLOPE)
C
C        ADD SEEPAGE TO RUNOFF
         RUNOFF=RUNOFF+TSEEP
C
C        ADJUST DEPTHS OF CANOPY LAYERS FOR ANY CHANGE IN SNOWPACK
         IF (NSP.GT.0 .OR. NSPLST.GT.0) THEN
            IF (NPLANT.NE.0) CALL CANOPY (NPLANT,NC,NSP,NSPLST,MZCINP,
     >       ITYPE,1,PLTHGT,PLTWGT,PLTLAI,RLEAF0,ZSP,ZC,TC,TCDT,BVAPC,
     >       VAPCDT,BWCAN,WCANDT,WCMAX,PCANDT,CANMA,CANMB,TMPDAY,HUMDAY)
         END IF
C
C        PRINT OUT WATER BALANCE FOR THIS HOUR
         IF (LVLOUT(7) .NE. 0)
     >    CALL WBALNC (NPLANT,NC,NSP,NR,NS,LVLOUT(7),JULIAN,HOUR,YEAR,
     >     ITYPE,1,ZC,BWCAN,WCANDT,BPCAN,PCANDT,BVAPC,VAPCDT,RHOSP,DZSP,
     >     DLWDT,WLAG,STORE,ZR,BGMC,GMCDT,BVAPR,VAPRDT,RHOR,ZS,BVLC,
     >     VLCDT,BVIC,VICDT,TOTFLO,PRECIP,RUNOFF,POND,EVAP1,ETSUM)
C
C        PRINT OUT ENERGY BALANCE FOR THIS HOUR
         IF (LVLOUT(6) .NE. 0)
     >     CALL ENERGY (NSPLST,LVLOUT(6),JULIAN,HOUR,YEAR,INITAL,
     >      TSWSNO,TLWSNO,TSWCAN,TLWCAN,TSWRES,TLWRES,TSWSOI,TLWSOI,
     >      THFLUX,EVAP1)
C
C        PRINT OUT FROST AND SNOW DEPTH
         IF (LVLOUT(10) .NE. 0) CALL FROST (NSP,NS,LVLOUT(10),JULIAN,
     >      HOUR,YEAR,1,ZSP,RHOSP,DZSP,DLWDT,WLAG,STORE,
     >      ZS,VLCDT,VICDT,ICESDT)
C
C        PRINT OUTPUT FOR THIS HOUR
         CALL OUTPUT (NPLANT,NC,NSP,NR,NS,LVLOUT,0,INPH2O,
     >    JULIAN,HOUR,YEAR,1,ZC,TCDT,VAPCDT,WCANDT,ROOTXT,RHOSP,
     >    ZSP,TSPDT,DLWDT,ZR,TRDT,VAPRDT,GMCDT,ZS,TSDT,VLCDT,VICDT,
     >    MATDT,TOTFLO,CONCDT,SALTDT)
C
C        PRINT OUTPUT TO SCREEN
         IF (LVLOUT(13) .NE. 0) THEN
            NPRINT=NPRINT+1
            IF (MOD(NPRINT,LVLOUT(13)).EQ.0) THEN
               WRITE (6,530) JULIAN,HOUR,YEAR,MINSTP,MAXSTP
               NPRINT=0
               MAXSTP=0
               MINSTP=0
            END IF
         END IF
C
         RETURN
C
      END IF
C
C
  522 FORMAT (I5,F5.2,F8.3,F6.3,F6.3,F8.2,E14.5,F8.3,F8.3,E14.5)
  530 FORMAT('+ Completed :      ',I4,3X,I4,3X,I4,5X,I4,6X,I4)
  535 FORMAT('+ Converg. Prob. : ',I4,3X,I4,3X,I4,5X,I4,6X,I4)
      END
C
C***********************************************************************
C
      BLOCK DATA
C
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /SWRCOE/ SOLCON,DIFATM,DIFRES,SNOCOF,SNOEXP
      COMMON /LWRCOF/ STEFAN,EMATM1,EMATM2,EMITC,EMITR,EMITSP,EMITS
C
      COMMON /SLPARM/ B(50),ENTRY(50),RHOB(50),SAT(50),SATK(50),
     >                NSALT,SALTKQ(10,50),VAPCOF(50),VAPEXP(50),
     >                SAND(50),SILT(50),CLAY(50),OM(50)
      COMMON /RSPARM/ RESMA,RESMB,RESMC,RESTKA,RESTKB
      COMMON /SPPROP/ G1,G2,G3,EXTSP,TKSPA,TKSPB,TKSPEX,VDIFSP,VAPSPX
      COMMON /SPWATR/ CLAG1,CLAG2,CLAG3,CLAG4,PLWMAX,PLWDEN,PLWHC,
     >                THICK,CTHICK
      COMMON /METASP/ CMET1,CMET2,CMET3,CMET4,CMET5,SNOMAX
C
      REAL LF,LS,LV
C
C**** DEFINE PHYSICAL CONSTANTS
C
      DATA  LF/335000./ LV/2500000./ LS/2835000./ G/9.81/ UGAS/8.3143/
     >      VONKRM/0.4/ VDIFF/.0000212/ P0/101300./
      DATA  RHOM/2650./ RHOL/1000./ RHOI/920./ RHOA/1.25/ RHOOM/1300./
      DATA  TKL/0.57/ TKI/2.2/ TKA/0.025/ TKR/ 0.05/ TKSA/8.8/
     >      TKSI/2.92/ TKCL/2.92/ TKOM/0.25/
      DATA  CM/900./ COM/1920./ CL/4200./ CI/2100./ CA/1006./ CV/1860./
     >      CR/1900./
      DATA  SOLCON/1360./ STEFAN/5.6697E-08/ 
C
C
C**** DEFINE HARD-CODED PARAMETERS
C
      DATA  EMATM1/0.261/ EMATM2/0.000777/ 
     >      EMITC/1.0/ EMITR/1.0/ EMITSP/1.0/ EMITS/1.0/
C
      DATA  DIFATM/0.76/ DIFRES/0.667/ SNOCOF/0.0/ SNOEXP/1.0/
C
C     DEFINE VAPOR COEFFICIENTS
      DATA VAPCOF/50*0.66/ VAPEXP/50*1.0/
C
      DATA RESMA/-53.72/ RESMB/1.32/ RESMC/0.0/ RESTKA/0.007/ 
     >     RESTKB/4.0/
C
C     EXTSP=1.77 FROM BARRY ET AL. 1990 WRR 26:1079-1092
      DATA  G1/0.16/ G2/0.0/ G3/110.0/ EXTSP/1.77/
     >      TKSPA/0.021/ TKSPB/2.51/ TKSPEX/2.0/ VDIFSP/0.00009/
     >      VAPSPX/14./
      DATA CLAG1/10.0/ CLAG2/1.0/ CLAG3/5.0/ CLAG4/450./ PLWMAX/0.1/
     >     PLWDEN/200./ PLWHC/0.03/ THICK/0.025/ CTHICK/0.05/ 
     >     CMET1/0.01/
      DATA CMET2/21./ CMET3/0.01/ CMET4/0.04/ CMET5/2.0/ SNOMAX/150./
C
      END
C***********************************************************************
C
      SUBROUTINE IOFILES (MTSTEP,INPH2O,MWATRXT,LVLOUT)
C
C     THIS SUBROUTINE OPENS ALL INPUT FILES AND OUTPUT FILES.
C     INPUT/OUTPUT FILES MUST BE LISTED IN A FILE (WHICH IS INPUT FROM
C     THE SCREEN) IN THE FOLLOWING ORDER:
C
C               FILE                                         UNIT
C         INPUT FROM SCREEN  --------------------------------   5
C         OUTPUT TO SCREEN  ---------------------------------   6
C         LIST OF INPUT/OUTPUT FILES  -----------------------  10
C         INPUT SITE DESCRIPTION DATA  ----------------------  11
C         INPUT WEATHER DATA  -------------------------------  12
C         INPUT MOISTURE PROFILES  --------------------------  13
C         INPUT TEMPERATURE PROFILES  -----------------------  14
C         INPUT WATER EXTRACTION  ---------------------------  15
C         INPUT CANOPY GROWTH  ------------------------------  41-48
C         OUTPUT GENERAL INFORMATION  -----------------------  21
C         OUTPUT MEASURED AND PREDICTED PROFILES  -----------  22
C         OUTPUT PREDICTED TEMPERATURE PROFILES  ------------  23
C         OUTPUT PREDICTED WATER CONTENT PROFILES -----------  24
C         OUTPUT PREDICTED WATER POTENTIAL PROFILES ---------  31
C         OUTPUT SURFACE ENERGY FLUX SUMMARY  ---------------  25
C         OUTPUT WATER BALANCE SUMMARY   --------------------  26
C         OUTPUT WATER FLUX BETWEEN SOIL NODES --------------  32
C         OUTPUT ROOT EXTRACTION ----------------------------  30
C         OUTPUT PREDICTED FROST DEPTH AND ICE CONTENT  -----  27
C         OUTPUT TOTAL SALT CONCENTRATION PROFILES   --------  28
C         OUTPUT SOLUTE CONCENTRATION PROFILES   ------------  29
C
C***********************************************************************
C
      INTEGER      IOS,LVLOUT(15)
      CHARACTER*1  ANSWER
      CHARACTER*80 IFILE,OFILE
C
C     -- Check for file created by ModShell.  If it exists, suppress
C     -- printing the screen message.
      OPEN (UNIT=99, FILE='SHAWMOD.FOF', STATUS='OLD', ERR=900,
     >               IOSTAT=IOS)
C
  900 IF (IOS.EQ.0) THEN
         CLOSE (99)
         OPEN (UNIT=10,FILE='SHAWMOD.FOF',STATUS='OLD')
        ELSE         
         WRITE (6,105)
         READ (5,100) IFILE
         WRITE (6,*)
         OPEN (UNIT=10,FILE=IFILE,STATUS='OLD',ERR=3,IOSTAT=INFIL)
    3    IF (INFIL.NE.0) THEN
C           CANNOT OPEN FILE FOR LIST OF INPUT/OUTPUT FILES
            WRITE (6,110) IFILE
            STOP
         END IF
      END IF
C
C     SPECIFY DATA FORMAT OF WEATHER DATA (HOURLY OR DAILY), WHETHER
C     THERE WILL BE ANY INPUT FOR A SINK TERM IN THE SOIL
      READ (10,*) MTSTEP,INPH2O,MWATRXT
C
C     OPEN INPUT FILES
      NINFIL=0
      DO 10 I=11,14
         READ (10,100) IFILE
         OPEN (UNIT=I,FILE=IFILE,STATUS='OLD',ERR=5,IOSTAT=INFIL)
    5    IF (INFIL.NE.0) THEN
C           CANNOT OPEN INPUT FILE
            WRITE (6,115) IFILE
            NINFIL=1
         END IF
   10 CONTINUE
C
C     OPEN WATER EXTRACTION FILE IF SPECIFIED
      IF (MWATRXT.GT.0) THEN
         READ (10,100) IFILE
         OPEN (UNIT=I,FILE=IFILE,STATUS='OLD',ERR=6,IOSTAT=INFIL)
    6    IF (INFIL.NE.0) THEN
C           CANNOT OPEN INPUT FILE
            WRITE (6,115) IFILE
            NINFIL=1
         END IF
      END IF
C
      IF (NINFIL.NE.0) THEN
C        SOME INPUT FILES COULD NOT BE OPENED
         WRITE (6,116) 
         STOP
      END IF
C
C
C**** ALLOW USER TO SPECIFY WHICH OUTPUT FILES ARE DESIRED
      READ (10,*) (LVLOUT(I), I=1,13)
      NEWFIL=0
C
C     OPEN GENERAL OUTPUT FILE
      READ (10,100) OFILE
      OPEN (UNIT=21,FILE=OFILE,STATUS='NEW',ERR=11,IOSTAT=NEW)
   11 IF (NEW .NE. 0) THEN
         OPEN (UNIT=21,FILE=OFILE)
         NEWFIL=1
         WRITE (6,120) 
         WRITE (6,130) OFILE
      END IF
C
C     OPEN PROFILES FILE
      READ (10,100) OFILE
      IF (LVLOUT(2) .GT. 0) THEN 
         OPEN (UNIT=22,FILE=OFILE,STATUS='NEW',ERR=12,IOSTAT=NEW)
   12    IF (NEW .NE. 0) THEN
            OPEN (UNIT=22,FILE=OFILE)
            IF (NEWFIL.EQ.0) WRITE (6,120) 
            NEWFIL=1
            WRITE (6,130) OFILE
         END IF
      END IF
C
C     OPEN SOIL TEMPERATURE FILE
      READ (10,100) OFILE
      IF (LVLOUT(3) .GT. 0) THEN
         OPEN (UNIT=23,FILE=OFILE,STATUS='NEW',ERR=13,IOSTAT=NEW)
   13    IF (NEW .NE. 0) THEN
            OPEN (UNIT=23,FILE=OFILE)
            IF (NEWFIL.EQ.0) WRITE (6,120) 
            NEWFIL=1
            WRITE (6,130) OFILE
         END IF
      END IF
C
C     OPEN SOIL WATER CONTENT FILE FILE
      READ (10,100) OFILE
      IF (LVLOUT(4) .GT. 0) THEN
         OPEN (UNIT=24,FILE=OFILE,STATUS='NEW',ERR=14,IOSTAT=NEW)
   14    IF (NEW .NE. 0) THEN
            OPEN (UNIT=24,FILE=OFILE)
            IF (NEWFIL.EQ.0) WRITE (6,120) 
            NEWFIL=1
            WRITE (6,130) OFILE
         END IF
      END IF
C
C     OPEN SOIL MATRIC POTENTIAL FILE
      READ (10,100) OFILE
      IF (LVLOUT(5) .GT. 0) THEN
         OPEN (UNIT=31,FILE=OFILE,STATUS='NEW',ERR=15,IOSTAT=NEW)
   15    IF (NEW .NE. 0) THEN
            OPEN (UNIT=31,FILE=OFILE)
            IF (NEWFIL.EQ.0) WRITE (6,120) 
            NEWFIL=1
            WRITE (6,130) OFILE
         END IF
      END IF
C
C     OPEN SURFACE ENERGY BALANCE FILE
      READ (10,100) OFILE
      IF (LVLOUT(6) .GT. 0) THEN
         OPEN (UNIT=25,FILE=OFILE,STATUS='NEW',ERR=16,IOSTAT=NEW)
   16    IF (NEW .NE. 0) THEN
            OPEN (UNIT=25,FILE=OFILE)
            IF (NEWFIL.EQ.0) WRITE (6,120) 
            NEWFIL=1
            WRITE (6,130) OFILE
         END IF
      END IF
C
C     OPEN WATER BALANCE SUMMARY FILE
      READ (10,100) OFILE
      IF (LVLOUT(7) .GT. 0) THEN
         OPEN (UNIT=26,FILE=OFILE,STATUS='NEW',ERR=17,IOSTAT=NEW)
   17    IF (NEW .NE. 0) THEN
            OPEN (UNIT=26,FILE=OFILE)
            IF (NEWFIL.EQ.0) WRITE (6,120) 
            NEWFIL=1
            WRITE (6,130) OFILE
         END IF
      END IF
C
C     OPEN FILE FOR WATER FLUX BETWEEN SOIL NODES
      READ (10,100) OFILE
      IF (LVLOUT(8) .GT. 0) THEN
         OPEN (UNIT=32,FILE=OFILE,STATUS='NEW',ERR=18,IOSTAT=NEW)
   18    IF (NEW .NE. 0) THEN
            OPEN (UNIT=32,FILE=OFILE)
            IF (NEWFIL.EQ.0) WRITE (6,120) 
            NEWFIL=1
            WRITE (6,130) OFILE
         END IF
      END IF
C
C     OPEN ROOT EXTRACTION FILE
      READ (10,100) OFILE
      IF (LVLOUT(9) .GT. 0) THEN
         OPEN (UNIT=30,FILE=OFILE,STATUS='NEW',ERR=19,IOSTAT=NEW)
   19    IF (NEW .NE. 0) THEN
            OPEN (UNIT=30,FILE=OFILE)
            IF (NEWFIL.EQ.0) WRITE (6,120) 
            NEWFIL=1
            WRITE (6,130) OFILE
         END IF
      END IF
C
C     OPEN FROST DEPTH AND ICE CONTENT FILE
      READ (10,100) OFILE
      IF (LVLOUT(10) .GT. 0) THEN
         OPEN (UNIT=27,FILE=OFILE,STATUS='NEW',ERR=20,IOSTAT=NEW)
   20    IF (NEW .NE. 0) THEN
            OPEN (UNIT=27,FILE=OFILE)
            IF (NEWFIL.EQ.0) WRITE (6,120) 
            NEWFIL=1
            WRITE (6,130) OFILE
         END IF
      END IF
C
C     OPEN SALT CONCENTRATION FILE
      READ (10,100) OFILE
      IF (LVLOUT(11) .GT. 0) THEN
         OPEN (UNIT=28,FILE=OFILE,STATUS='NEW',ERR=21,IOSTAT=NEW)
   21    IF (NEW .NE. 0) THEN
            OPEN (UNIT=28,FILE=OFILE)
            IF (NEWFIL.EQ.0) WRITE (6,120) 
            NEWFIL=1
            WRITE (6,130) OFILE
         END IF
      END IF
C
C     OPEN SOLUTE CONCENTRATION FILE
      READ (10,100) OFILE
      IF (LVLOUT(12) .GT. 0) THEN
         OPEN (UNIT=29,FILE=OFILE,STATUS='NEW',ERR=22,IOSTAT=NEW)
   22    IF (NEW .NE. 0) THEN
            OPEN (UNIT=29,FILE=OFILE)
            IF (NEWFIL.EQ.0) WRITE (6,120) 
            NEWFIL=1
            WRITE (6,130) OFILE
         END IF
      END IF
C
      IF (NEWFIL .EQ. 1) THEN
C        OUTPUT FILES AREADY EXIST -- CHECK IF OK TO WRITE OVER THEM
         WRITE (6,140)
         READ (5,'(A)') ANSWER
         IF (ANSWER.NE.'Y' .AND. ANSWER.NE.'y') THEN
            WRITE (6,145) 
            STOP
         END IF
      END IF
C
      IF (IOS .NE. 0) WRITE (6,150)
      CLOSE (UNIT=10)
C
  100 FORMAT(A80)
  105 FORMAT(//,
     >' >>>>>>> Simultaneous Heat And Water (SHAW) Model <<<<<<<',/,
     >'                        Version 2.4',//,
     >' Enter the file containing the list of input/output files: ')
  110 FORMAT(/,' THE FOLLOWING INPUT FILE CANNOT BE OPENED: ',A80,
     >//,' Check that the filename and path for the file containing'
     > /,' the list of input/output files are correct.',/)
  115 FORMAT(/,' THE FOLLOWING INPUT FILE CANNOT BE OPENED: ',A80)
  116 FORMAT(/,' Check that filenames and paths listed in the file',
     >/' containing the list of input/output files are correct.',/)
  120 FORMAT(//,' The following output files already exist:',/)
  130 FORMAT(5X,A80)
  140 FORMAT(/,' DO YOU WISH TO WRITE OVER THESE FILES? (Y/N): ')
  145 FORMAT(/,' Rename above files or change output filenames in ',/,
     >       ' the file containing the list of input/output files.',/)
  150 FORMAT(/,' Simulation in progress . . .')
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE INPUT (NC,NSP,NR,NS,TOLER,LEVEL,MTSTEP,MZCINP,MPLTGRO,
     >  INPH2O,MWATRXT,LVLOUT,IVLCBC,ITMPBC,NPLANT,PLTHGT,PLTWGT,PLTLAI,
     >  ROOTDP,DCHAR,TCCRIT,RSTOM0,RSTEXP,PLEAF0,RLEAF0,RROOT0,
     >  CANALB,CANMA,CANMB,WCMAX,LANGLE,ITYPE,ZC,WCANDT,
     >  ZSP,DZSP,RHOSP,TSPDT,DLWDT,ICESPT,SNOTMP,
     >  GMCDT,GMCMAX,RLOAD,ZRTHIK,COVER,ALBRES,RESCOF,
     >  ZS,TSDT,VLCDT,VICDT,MATDT,CONCDT,ICESDT,SALTDT,ALBDRY,ALBEXP,
     >  DGRADE,SLTDIF,ASALT,DISPER,
     >  ZMSRF,ZHSRF,ZERSRF,ZMSP,ZHSP,HEIGHT,PONDMX,
     >  JSTART,YRSTAR,HRSTAR,JEND,YREND,NHRPDT,WDT,
     >  ALATUD,SLOPE,ASPECT,HRNOON)
C
C     THIS SUBROUTINE IS USED TO INPUT ALL GENERAL INFORMATION AND
C     INITIAL CONDITIONS.
C
C***********************************************************************
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /SWRCOE/ SOLCON,DIFATM,DIFRES,SNOCOF,SNOEXP
      COMMON /LWRCOF/ STEFAN,EMATM1,EMATM2,EMITC,EMITR,EMITSP,EMITS
      COMMON /CLAYRS/ DRYCAN(8,10),CANLAI(8,10),TOTLAI(8),IEVAP(8),
     >                RLEAF(8,10),RROOT(8,50),ROOTDN(8,50),TOTROT(8)
      COMMON /SPPROP/ G1,G2,G3,EXTSP,TKSPA,TKSPB,TKSPEX,VDIFSP,VAPSPX
      COMMON /SPWATR/ CLAG1,CLAG2,CLAG3,CLAG4,PLWMAX,PLWDEN,PLWHC,
     >                THICK,CTHICK
      COMMON /METASP/ CMET1,CMET2,CMET3,CMET4,CMET5,SNOMAX
      COMMON /RSPARM/ RESMA,RESMB,RESMC,RESTKA,RESTKB
      COMMON /SLPARM/ B(50),ENTRY(50),RHOB(50),SAT(50),SATK(50),
     >                NSALT,SALTKQ(10,50),VAPCOF(50),VAPEXP(50),
     >                SAND(50),SILT(50),CLAY(50),OM(50)
      COMMON /MEASUR/ MEASDY,MEASHR,VLCMES(50),TSMEAS(50)
      REAL LF,LS,LV
C
      REAL DCHAR(8),CANALB(8),TCCRIT(8),RSTOM0(8),RSTEXP(8),PLEAF0(8),
     >     RLEAF0(8),RROOT0(8),ZC(11)
      REAL ZSP(100),DZSP(100),RHOSP(100)
      REAL DGRADE(10),SLTDIF(10),ASALT(50),DISPER(50)
C
      INTEGER LANGLE(8),ITYPE(8)
C
C
      REAL ZS(50),TSDT(50),MATDT(50),CONCDT(10,50),VLCDT(50),VICDT(50),
     >     SALTDT(10,50)
      REAL GMCDT(10)
      REAL TSPDT(100),DLWDT(100)
      REAL WCANDT(10)
      REAL PLTHGT(8),PLTWGT(8),PLTLAI(8),ROOTDP(8)
      INTEGER YRSTAR,HRSTAR,YREND
      INTEGER ICESDT(50),ICESPT(100),LEVEL(6),LVLOUT(15)
      CHARACTER*80 TITLE,IFILE
C
C     OPEN INPUT AND OUTPUT FILES
      CALL IOFILES (MTSTEP,INPH2O,MWATRXT,LVLOUT)
C
C**** INPUT GENERAL INFORMATION FOR THE RUN
      READ (11,100) TITLE
      READ (11,*) JSTART,HRSTAR,YRSTAR,JEND,YREND
      READ (11,*) ALTDEG,ALTMIN,SLP,ASPEC,HRNOON,ELEV
      READ (11,*) NPLANT,NSP,NR,NS,NSALT,TOLER,NHRPDT,
     >           (LEVEL(I), I=1,6)
      IF (HRSTAR .EQ. 0) THEN
         HRSTAR=24
         JSTART=JSTART-1
      END IF
      PRESUR=101300.*EXP(-ELEV/8278.)
C
      WRITE (21,102) TITLE
      WRITE (21,105) JSTART,HRSTAR,YRSTAR,JEND,YREND
      WRITE (21,110) ALTDEG,ALTMIN,SLP,ASPEC,HRNOON,ELEV,PRESUR/1000.
      IF (ALTDEG.LT.0.0.AND.ALTMIN.GT.0.0) ALTMIN=-ALTMIN
      ALATUD=(ALTDEG + ALTMIN/60.)*3.14159/180.
      SLOPE=ATAN(SLP/100.)
      ASPECT=ASPEC*3.14159/180.
      IF (NHRPDT.EQ.1) THEN
         WDT=0.6
        ELSE
         WDT=1.0
      END IF
C
C     CHECK IF NHRPDT IS VALID
      IF (MOD(24,NHRPDT) .NE. 0) THEN
        WRITE (6,*)
        WRITE (6,*) ' ************************************************'
        WRITE (6,*) ' PARAMETER FOR THE NUMBER OF HOURS PER TIME STEP'
        WRITE (6,*) ' (NHRPDT; LINE D OF INPUT FILE) IS NOT VALID. '
        WRITE (6,*) ' NHRPDT MUST BE EVENLY DIVISIBLE INTO 24 HOURS.'
        WRITE (6,*) ' ************************************************'
        WRITE(21,*)
        WRITE(21,*)
        WRITE(21,*) ' ************************************************'
        WRITE(21,*) ' THE PARAMETER FOR THE NUMBER OF HOURS PER TIME'
        WRITE(21,*) ' (NHRPDT; LINE D OF INPUT FILE) IS NOT VALID. '
        WRITE(21,*) ' NHRPDT MUST BE EVENLY DIVISIBLE INTO 24 HOURS.'
        WRITE(21,*) ' ************************************************'
        STOP
      END IF 
C
C     CHECK IF HRSTAR IS COMPATILBE WITH NHRPDT
      IF (MOD(HRSTAR,NHRPDT) .NE. 0) THEN
        WRITE (6,*)
        WRITE (6,*)' **************************************************'
        WRITE (6,*)' THE PARAMETER FOR THE NUMBER OF HOURS PER TIME'
        WRITE (6,*)' STEP (NHRPDT; LINE D OF INPUT FILE) IS NOT' 
        WRITE (6,*)' COMPATIBLE WITH THE BEGINNING HOUR FOR THE'
        WRITE (6,*)' SIMULATION (HRSTAR; LINE B OF INPUT FILE). ' 
        WRITE (6,*)' HRSTAR MUST BE ZERO OR EVENLY DIVISIBLE BY NHRPDT.'
        WRITE (6,*)' **************************************************'
        WRITE(21,*)
        WRITE(21,*)
        WRITE(21,*)' **************************************************'
        WRITE(21,*)' THE PARAMETER FOR THE NUMBER OF HOURS PER TIME'
        WRITE(21,*)' STEP (NHRPDT; LINE D OF INPUT FILE) IS NOT' 
        WRITE(21,*)' COMPATIBLE WITH THE BEGINNING HOUR FOR THE'
        WRITE(21,*)' SIMULATION (HRSTAR; LINE B OF INPUT FILE).'
        WRITE(21,*)' HRSTAR MUST BE ZERO OR EVENLY DIVISIBLE BY NHRPDT.'
        WRITE(21,*)' **************************************************'
        STOP
      END IF
C
C**** INPUT AERODYNAMIC ATMOSPHERIC AND SURFACE PROPERTIES
      READ (11,*) ZMCM,HEIGHT,PONDCM
      ZMSRF=ZMCM/100.
      ZHSRF=0.2*ZMSRF
      ZHCM=ZHSRF*100.
      PONDMX=PONDCM/100.
      ZERSRF=0.0
      WRITE (21,115) ZMCM,ZHCM,ZERSRF,HEIGHT,PONDCM,DIFATM,EMATM1,EMATM2
C
C**** INPUT PROPERTIES OF CANOPY
      IF (NPLANT .GT. 0) THEN
         READ (11,*) MPLTGRO,CANMA,CANMB,WCANDT(1)
C
C        DETERMINE MAXIMUM WATER CONTENT OF STANDING DEAD (AT 99.9% RH)
         CALL CANHUM (2,0.999,DUMMY,WCMAX,0.0,CANMA,CANMB)
         IF (WCANDT(1) .GT. WCMAX) WCANDT(1)=WCMAX
C
         IF (MPLTGRO .EQ. 2) THEN
C           SET FLAG FOR USER INPUT OF NODE SPACING
            MPLTGRO=0
            MZCINP=1
           ELSE
C           MODEL WILL GENERATE SPACING OF NODES WITHIN THE CANOPY
            MZCINP=0
         END IF
         DO 1010 J=1,NPLANT
            READ (11,*) ITYPE(J),LANGLE(J),CANALB(J),
     >      TCCRIT(J),RSTOM0(J),RSTEXP(J),PLEAF0(J),RLEAF0(J),RROOT0(J)
 1010    CONTINUE
         WRITE (21,170) EMITC,CANMA,CANMB,WCANDT(1),WCMAX,
     >                  ('  PLANT #',J, J=1,NPLANT)
         WRITE (21,171) (ITYPE(J), J=1,NPLANT)
         WRITE (21,172) (LANGLE(J), J=1,NPLANT)
         WRITE (21,173) (CANALB(J), J=1,NPLANT)
         WRITE (21,174) (TCCRIT(J), J=1,NPLANT)
         WRITE (21,175) (RSTOM0(J), J=1,NPLANT)
         WRITE (21,176) (RSTEXP(J), J=1,NPLANT)
         WRITE (21,177) (PLEAF0(J), J=1,NPLANT)
         WRITE (21,178) (RLEAF0(J), J=1,NPLANT)
         WRITE (21,179) (RROOT0(J), J=1,NPLANT)
C
C        CONVERT RLEAF0 AND RROOT0 TO FROM RESISTANCE TO CONDUCTANCE
         DO 1012 J=1,NPLANT
            RLEAF0(J)=1./RLEAF0(J)
            RROOT0(J)=1./RROOT0(J)
 1012    CONTINUE
C
         IF (MZCINP.NE.0)THEN
C           USER INPUT OF ROOT DISTRIBUTION AND NODE SPACING WITHIN THE
C           CANOPY (ASSUMED TO BE CONSTANT FOR ENTIRE SIMULATION)
            MPLTGRO=0
            WRITE (21,180) (I, I=1,NS)
            DO 1015 J=1,NPLANT
               READ (11,*) (ROOTDN(J,I), I=1,NS)
               WRITE (21,181) J,(ROOTDN(J,I), I=1,NS)
 1015       CONTINUE
C
            WRITE (21,185) ('        PLANT SPECIES #',J, J=1,NPLANT)
            WRITE (21,186) '     DEPTH  DEAD MATL ',
     >                     ('  CHAR.DIM. BIOMASS  LAI', J=1,NPLANT)
            WRITE (21,186) '      (M)    (KG/KG)  ',
     >                     ('     (CM)   (KG/M2)     ', J=1,NPLANT)
            READ (11,*) NC
            DO 1020 I=1,NC
               READ (11,*) ZC(I),
     >              (DCHAR(J),DRYCAN(J,I),CANLAI(J,I), J=1,NPLANT)
               WCANDT(I)=WCANDT(1)
               WRITE (21,187) ZC(I),WCANDT(I),
     >              (DCHAR(J),DRYCAN(J,I),CANLAI(J,I), J=1,NPLANT)
 1020       CONTINUE
            READ (11,*) ZC(NC+1)
            WRITE (21,188) ZC(NC+1)
C
C           INITIALIZE PLANT HEIGHT, WEIGHT, LAI, AND ROOTS
            DO 1060 J=1,NPLANT
               PLTHGT(J)=0.0
               PLTWGT(J)=0.0
               PLTLAI(J)=0.0
               TOTROT(J)=0.0
               DCHAR(J)=DCHAR(J)/100.
               DO 1025 I=NC,1,-1
                  IF (CANLAI(J,I).GT.0.0) PLTHGT(J)=ZC(NC+1)-ZC(I)
                  PLTWGT(J)=PLTWGT(J)+DRYCAN(J,I)
                  PLTLAI(J)=PLTLAI(J)+CANLAI(J,I)
C                 CONVERT LEAF DIMENSION FROM CM TO METERS
 1025          CONTINUE
               TOTLAI(J)=PLTLAI(J)
               IF (ITYPE(J) .NE. 0) THEN
C                 TRANSPIRING PLANT -- CALCULATE TOTAL ROOT DENSITY
C                 (FRACTION OF ROOTS IN SOIL COLUMN)
                  DO 1030 I=1,NS
                    TOTROT(J)=TOTROT(J) + ROOTDN(J,I)
 1030             CONTINUE
C                 CALCULATE EFFECTIVE LEAF CONDUCT. FOR EACH CANOPY LAYR
                  DO 1040 I=1,NC
                     RLEAF(J,I)=RLEAF0(J)*CANLAI(J,I)/TOTLAI(J)
 1040             CONTINUE
C                 CALCULATE EFFECTIVE ROOT CONDUC. FOR EACH SOIL LAYER
                  DO 1050 I=1,NS
                     RROOT(J,I)=RROOT0(J)*ROOTDN(J,I)/TOTROT(J)
 1050             CONTINUE
               END IF
 1060       CONTINUE
           ELSE
C
            IF (MPLTGRO.GT.0) THEN
C              OPEN FILES FOR PLANT GROWTH
               NINFIL = 0
               DO 1070 I=1,NPLANT
                  READ (11,100) IFILE
                  OPEN (UNIT=40+I,FILE=IFILE,STATUS='OLD',
     >                  ERR=1065,IOSTAT=INFIL)
 1065             IF (INFIL.NE.0) THEN
C                    CANNOT OPEN INPUT FILE
                     WRITE (6,190) IFILE
                     NINFIL=1
                  END IF
 1070          CONTINUE
               IF (NINFIL.NE.0) THEN
C                 SOME INPUT FILES COULD NOT BE OPENED
                  WRITE (6,192) 
                  STOP
               END IF
C
             ELSE
C              READ IN PLANT PARAMETERS (ASSUMED TO BE CONSTANT FOR
C              ENTIRE SIMULATION)
               DO 1080 J=1,NPLANT
                  READ (11,*) PLTHGT(J),DCHAR(J),PLTWGT(J),PLTLAI(J),
     >                        ROOTDP(J)
C                 CONVERT LEAF DIMENSION FROM CM TO METERS
                  DCHAR(J)=DCHAR(J)/100.
 1080          CONTINUE
            END IF
         END IF
      END IF
C
C
C**** INPUT PROPERTIES FOR SNOWPACK
      IF (NSP .GE. 0) THEN
C        NSP = 0  ---> SNOW IS NOT PRESENT, BUT IT MAY SNOW
         READ (11,*) SNOTMP,ZMSPCM
         ZMSP=ZMSPCM/100.
         ZHSP=0.2*ZMSP
         ZHSPCM=ZHSP*100.
         WRITE (21,120) SNOTMP,ZMSPCM,ZHSPCM,G1,G2,G3,EXTSP,EMITSP,
     >               SNOCOF,SNOEXP,TKSPA,TKSPB,TKSPEX,VDIFSP,VAPSPX,
     >               CLAG1,CLAG2,CLAG3,CLAG4,PLWMAX,PLWHC,PLWDEN,
     >               CMET1,CMET2,CMET3,CMET4,CMET5,SNOMAX
         IF (NSP .GT. 0) THEN
C           SNOWPACK IS PRESENT AT THE START OF THE SIMULATION
            DO 5 I=1,NSP
               READ (11,*) DZSP(I),TSPDT(I),DLWDT(I),RHOSP(I)
               ICESPT(I) = 0
               IF (DLWDT(I) .GT. 0.0) ICESPT(I) = 1
    5       CONTINUE
            ZSP(1)=0.0
            TDEPTH = DZSP(1)
            DO 10 I=2,NSP
               ZSP(I) = TDEPTH + DZSP(I)/2.
               TDEPTH = TDEPTH + DZSP(I)
   10       CONTINUE
            ZSP(NSP+1)=TDEPTH
         END IF
C
       ELSE
C        NSP < 0  ---> SNOW WILL NOT BE CONSIDERED FOR THE SIMULATION
         SNOTMP = -200.
         NSP = 0
      END IF
C
C
C**** INPUT PROPERTIES FOR RESIDUE LAYERS
      IF (NR .GT. 0) THEN
         READ (11,*) COVER,ALBRES,RLOAD,ZRTHIK,
     >               GMCDT(1),RESCOF
C        DETERMINE MAXIMUM WATER CONTENT OF RESIDUE (AT 99.9% RH)
         CALL RESHUM (2,0.999,DUMMY,GMCMAX,0.0)
         IF (GMCDT(1) .GT. GMCMAX) GMCDT(1)=GMCMAX
C
         WRITE (21,125) COVER,ALBRES,EMITR,RLOAD,ZRTHIK,
     >                  RESMA,RESMB,RESTKA,RESTKB,RESCOF,GMCMAX
C        CONVERT RLOAD FROM KG/HA TO KG/M2
         RLOAD=RLOAD/10000.
C        CONVERT ZRTHIK FROM CM TO M
         ZRTHIK=ZRTHIK/100.
      END IF
C
C**** INPUT PROPERTIES AND INITIAL CONDITIONS FOR EACH TYPE OF SOLUTE
      IF (NSALT .GT. 0) WRITE (21,135)
      DO 30 I=1,NSALT
         READ (11,*) SLTDIF(I),HAFLIF
         READ (11,*) (SALTKQ(I,J), J=1,NS)
         READ (11,*) (SALTDT(I,J), J=1,NS)
         WRITE (21,140) I,SLTDIF(I),HAFLIF
         WRITE (21,145) (SALTKQ(I,J), J=1,NS)
C        CALCULATE EXPONENT FOR DEGRADATION BASED ON HALF LIFE
         IF (HAFLIF .EQ. 0.0) THEN
C           INDICATES NO DEGRADATION OF SOLUTE SPECIES
            DGRADE(I)=0.0
           ELSE
            DGRADE(I)=0.693147/HAFLIF
         END IF
   30 CONTINUE
C
C**** INPUT INITIAL CONDITIONS OF TEMP AND MOISTURE FOR SOIL LAYERS
   40 READ (14,*,END=90) JDAY,JHR,JYR,(TSDT(I), I=1,NS)
      IF (JDAY.NE.JSTART .OR. JHR.NE.HRSTAR .OR. JYR.NE.YRSTAR) THEN
         IF (HRSTAR.EQ.24 .AND. JHR.EQ.0) THEN 
            IF (JSTART.NE.JDAY-1 .OR. YRSTAR.NE.JYR) GOTO 40
           ELSE 
            GO TO 40
         END IF
      END IF
   50 READ (13,*,END=95) JDAY,JHR,JYR,(VLCDT(I), I=1,NS)
      IF (JDAY.NE.JSTART .OR. JHR.NE.HRSTAR .OR. JYR.NE.YRSTAR) THEN
         IF (HRSTAR.EQ.24 .AND. JHR.EQ.0) THEN 
            IF (JSTART.NE.JDAY-1 .OR. YRSTAR.NE.JYR) GOTO 50
           ELSE
            GO TO 50
         END IF
      END IF
C     SAVE DAY AND HOUR OF LAST TEMPERATURE AND MOISTURE MEASUREMENTS
      MEASDY=JSTART
      MEASHR=HRSTAR
C
C**** INPUT PROPERTIES FOR SOIL LAYERS
      READ (11,*) IVLCBC,ITMPBC,ALBDRY,ALBEXP
      WRITE (21,150)
      IF (IVLCBC.LE.0) THEN
         WRITE (21,151)
        ELSE
         WRITE (21,152)
      END IF
      IF (ITMPBC.LE.0) THEN
         WRITE (21,153)
        ELSE
         WRITE (21,154)
      END IF
      WRITE (21,155) ALBDRY,ALBEXP,EMITS
      IF (NSALT.EQ.0) THEN
         WRITE (21,156)
        ELSE
         WRITE (21,157)
      END IF
      DO 60 I=1,NS
         IF (NSALT .EQ. 0) THEN
            READ (11,*) ZS(I),B(I),ENTRY(I),SATCON,RHOB(I),SAT(I),
     >                SAND(I),SILT(I),CLAY(I),OM(I)
           ELSE
            READ (11,*) ZS(I),B(I),ENTRY(I),SATCON,RHOB(I),SAT(I),
     >                SAND(I),SILT(I),CLAY(I),OM(I),ASALT(I),DISPER(I)
         END IF
C        ASSURE THAT TEXTURE SUMS TO 100% AND CONVERT TO FRACTION
         TEXTOT=SAND(I)+SILT(I)+CLAY(I)
         CF=1./TEXTOT
         SAND(I)=SAND(I)*CF
         SILT(I)=SILT(I)*CF
         CLAY(I)=CLAY(I)*CF
         OM(I)=OM(I)/100.
C        CHECK IF SATURATION VALUE IS NOT TOO HIGH FOR DENSITY
C        (ADJUST SPECIFIC DENSITY FOR ORGANIC MATTER)
         IF (SAT(I) .GT. (1.- RHOB(I)*((1.-OM(I))/RHOM+OM(I)/RHOOM)))
     >             SAT(I)=1.- RHOB(I)*((1.-OM(I))/RHOM+OM(I)/RHOOM)
C        CONVERT CONDUCTIVITY FROM CM/HR TO M/SEC
         SATK(I) = SATCON/360000.
         IF (NSALT .EQ. 0) THEN
            WRITE (21,158)ZS(I),B(I),ENTRY(I),SATCON,RHOB(I),
     >      SAND(I)*100.,SILT(I)*100.,CLAY(I)*100.,OM(I)*100.
          ELSE
            WRITE (21,158)ZS(I),B(I),ENTRY(I),SATCON,RHOB(I),
     >      SAND(I)*100.,SILT(I)*100.,CLAY(I)*100.,OM(I)*100.,
     >      ASALT(I),DISPER(I)
         END IF
C        DEFINE MATRIC POTENTIAL OR WATER CONTENT OF SOIL AND DETERMINE
C        WHETHER SOIL IS FROZEN
         IF (INPH2O.NE.1) THEN
C           INPUT SOIL MOISTURE IS WATER CONTENT
            IF (VLCDT(I) .GT. SAT(I)) VLCDT(I)=SAT(I)
            CALL MATVL1 (I,MATDT(I),VLCDT(I),DUMMY)
C           SAVE INITIAL CONDITIONS FOR SUBROUTINE OUTPUT
            TSMEAS(I)=TSDT(I)
            VLCMES(I)=VLCDT(I)
           ELSE
C           INPUT SOIL MOISTURE IS MATRIC POTENTIAL
            MATDT(I)=VLCDT(I)
            CALL MATVL2 (I,MATDT(I),VLCDT(I),DUMMY)
C           SAVE INITIAL CONDITIONS FOR SUBROUTINE OUTPUT
            TSMEAS(I)=TSDT(I)
            VLCMES(I)=MATDT(I)
         END IF 
         VICDT(I)=0.0
         ICESDT(I)=0
C===     (do not adjust if frozen - this is done on entry into GOSHAW 
C===     and problems result if it is done twice when INPH2O=1
C===     IF (TSDT(I) .LE. 0.0) CALL FROZEN (I,VLCDT,VICDT,MATDT,
C=== >                                 CONCDT,TSDT,SALTDT,ICESDT)
C        DEFINE SOLUTE CONCENTRATION
         DO 55 J=1,NSALT
            CONCDT(J,I)=SALTDT(J,I)/(SALTKQ(J,I)+VLCDT(I)*RHOL/RHOB(I))
   55    CONTINUE
   60 CONTINUE
C
      WRITE (21,160)
      RETURN
C
   90 WRITE (6,*) 'CANNOT FIND SOIL TEMPERATURE DATA FOR INITIAL HOUR'
      WRITE (21,*) 'CANNOT FIND SOIL TEMPERATURE DATA FOR INITIAL HOUR'
      STOP
   95 WRITE (6,*) 'CANNOT FIND SOIL MOISTURE DATA FOR INITIAL HOUR'
      WRITE (21,*) 'CANNOT FIND SOIL MOISTURE DATA FOR INITIAL HOUR'
      STOP
C
  100 FORMAT (A80)
  102 FORMAT (' ',A80)
  105 FORMAT (/5X,'SIMULATION BEGINS ON DAY',I4,', HOUR',I3,', OF ',I4,
     >        /5X,'SIMULATION ENDS   ON DAY',I4,', HOUR 24  OF ',I4)
  110 FORMAT (//,' GENERAL SITE DESCRIPTION',//5X,'LATITUDE  :',F5.1,
     >' DEGREES ',F5.2,' MINUTES',/5X,'SLOPE     : ',F5.2,' %',
     >/5X,'ASPECT    : ',F5.1,' DEGREES FROM NORTH',/5X,'SOLAR NOON: ',
     >F5.2,' HOURS',/5X,'ELEVATION : ',F5.0,' M'/5X,'PRESSURE  : ',
     >F5.1,' KPA')
  115 FORMAT (//,' WIND PROFILE AND SURFACE PARAMETERS',//5X,
     >'ZM =',F5.2,' CM (WHEN SNOW AND CANOPY ARE NOT PRESENT)',/5X,
     >'ZH =',F5.2,' CM (WHEN SNOW AND CANOPY ARE NOT PRESENT)',/5X,
     >'ZERO PLANE DISPLACEMENT =',F5.2,' M',/5X,
     >'HEIGHT OF INSTRUMENTATION =',F5.1,' M',/5X,
     >'MAXIMUM DEPTH OF PONDING =',F4.1,' CM',///1X,
     >'ATMOSPHERIC CONDITIONS',//5X,'MAXIMUM CLEAR-SKY TRANSMISSIVITY ='
     >,F5.2,/5X,'CLEAR-SKY LONG-WAVE EMISSIVITY PARAMETERS :',
     >F6.3,E15.3)
  170 FORMAT (//,' CANOPY PARAMETERS',//5X,
     >'EMISSIVITY OF PLANT MATERIAL :',F5.2,/5X,
     >'MOISTURE PARAMETERS FOR ANY DEAD PLANT MATERIAL :',F7.2,F5.2,
     >/5X,'INITIAL MOISTURE CONTENT FOR DEAD PLANT MATERIAL :',F5.2,
     >' KG/KG',
     >/5X,'MAXIMUM MOISTURE CONTENT FOR DEAD PLANT MATERIAL :',F5.2,
     >' KG/KG',//40X,8(A,I1))
  171 FORMAT (5X,'PLANT TYPE (0=DEAD,1=TRANSPIRING) :',8I10)
  172 FORMAT (5X,'LEAF ANGLE (0=VERTICLE,1=RANDOM)  :',8I10)
  173 FORMAT (5X,'ALBEDO                            :',8F10.2)
  174 FORMAT (5X,'MIN. TRANSPIRATION TEMPERATURE (C):',8F10.1)
  175 FORMAT (5X,'MINIMUM STOMATAL RESISTANCE (S/M) :',8F10.0)
  176 FORMAT (5X,'STOMATAL RESISTANCE EXPONENT      :',8F10.1)
  177 FORMAT (5X,'CRITICAL LEAF WATER POTENTIAL (M) :',8F10.0)
  178 FORMAT (5X,'LEAF RESISTANCE (KG/M2-S)         :',8E10.3)
  179 FORMAT (5X,'ROOT RESISTANCE (KG/M2-S)         :',8E10.3)
  180 FORMAT (/5X,'ROOT DENSITY PROFILE:',
     >/5X,'SOIL LAYER            :',50I5)
  181 FORMAT ( 5X,'PLANT #',I1,' ROOT FRACTION: ',50F5.2)
  185 FORMAT (/5X,'        H2O IN',8(A,I1))
  186 FORMAT (A,8A)
  187 FORMAT (5X,F5.2,F9.2,4X,8(F8.2,F9.2,F7.2))
  188 FORMAT (5X,F5.2,'  ==>  BOTTOM DEPTH OF CANOPY')
  190 FORMAT(/,' THE FOLLOWING INPUT FILE CANNOT BE OPENED: ',A80)
  192 FORMAT(/,' Check that filenames and paths for the plant growth',
     >/' files listed in the site characteristics file are correct.',/)
  120 FORMAT (//,' SNOW PARAMETERS',//5X,
     >'MAXIMUM TEMPERATURE FOR SNOWFALL =',F5.2,' C',/5X,
     >'WIND-PROFILE PARAMETERS FOR SNOW: ZM =',F5.2,' CM    ZH =',F5.2,
     >' CM',/5X,
     >'GRAIN-SIZE DIAMETER PARAMETERS :',F5.2,F6.2,F8.1,/5X,
     >'SOLAR RADIATION EXTINCTION COEFFICIENT PARAMETER =',F6.2,/5X,
     >'EMISSIVITY OF SNOWPACK =',F5.2,/5X,
     >'COEFFICIENT AND EXPONENT FOR ALBEDO :',F5.2,F7.2,/5X,
     >'COEFFICIENTS AND EXPONENT FOR THEMAL COND.:',F6.3,F7.2,F7.1,/5X,
     >'VAPOR DIFFUSIVITY AND TEMP-DEPENDENCE EXPONENT :',F8.5,' M**2',
     >F7.1,/5X,'LAG COEFFICIENTS FOR SNOWCOVER OUTFLOW :',F5.1,F6.1,
     >F6.1,F8.1,/5X,'MAX AND MIN WATER HOLDING CAPACITY AND WHC DENSITY
     >=',F5.2,' M/M',F7.3,' M/M',F8.1,' KG/M**3',/5X,
     >'COMPACTION AND SETTLING PARAMETERS :',F5.2,F7.1,F7.2,F7.2,F7.1,
     >F8.1)
  125 FORMAT (//,' RESIDUE PARAMETERS',//5X,
     >'FRACTION OF GROUND COVERED BY RESIDUE:',F5.2,/5X,
     >'ALBEDO AND EMISSIVITY OF RESIDUE :',F5.2,F7.2,/5X,
     >'RESIDUE LOADING =',F7.0,' KG/HA',/5X,
     >'THICKNESS OF RESIDUE LAYER =',F5.2,' CM',/5X,
     >'MOISTURE PARAMETERS FOR RESIDUE :',F7.2,F5.2,/5X,
     >'THERMAL CONVECTION PARAMETERS :',F6.3,F6.2,/5X,
     >'VAPOR TRANSFER RESISTANCE FROM RESIDUE:',F7.0,/5X,
     >'MAXIMUM MOISTURE CONTENT FOR RESIDUE :',F5.2,' KG/KG')
  135 FORMAT (//,' SOLUTE PROPERTIES')
  140 FORMAT (/5X,'SALT SPECIES #',I1,':',
     >/5X,'DIFFUSIVITY  =',E10.3,' M**2/S'
     >,/5X,'HALF-LIFE =',F6.1,' DAYS (ZERO INDICATES NO DEGRADATION)',
     >/5X,'SOIL MATRIX-SOIL SOLUTION PARTITIONING COEEF (Kd) ',
     >'FOR EACH SOIL NODE:')
  145 FORMAT (5X,15F7.2,/5X,5F7.2)
  150 FORMAT (//,' SOIL PROPERTIES')
  151 FORMAT (/5X,
     >'INPUT WATER CONTENT SPECIFIED FOR LOWER BOUNDARY OF WATER FLUX')
  152 FORMAT (/5X,
     >'UNIT GRADIENT SPECIFIED FOR LOWER BOUNDARY OF WATER FLUX')
  153 FORMAT (5X,'INPUT SOIL TEMPERATURE SPECIFIED FOR LOWER BOUNDARY')
  154 FORMAT (5X,
     >'SOIL TEMPERATURE AT LOWER BOUNDARY ESTIMATED BY MODEL')
  155 FORMAT (/5X,
     >'ALBEDO OF DRY SOIL AND MOISTURE-DEPENDENCE EXPONENT :',F5.2,F7.2,
     >/5X,'EMISSIVITY OF SOIL =',F5.2)
  156 FORMAT (/5X,
     >'DEPTH  B-VALUE  AIR ENTRY  K-SAT  DENSITY  SAND SILT CLAY  OM',
     >/5X,
     >' (M)               (M)    (CM/HR) (KG/M3)   (%)  (%)  (%)  (%)')
  157 FORMAT (/5X,
     >'DEPTH  B-VALUE  AIR ENTRY  K-SAT  DENSITY  SAND SILT CLAY  OM  ',
     >' SALT PARAMS',/5X,
     >' (M)               (M)    (CM/HR) (KG/M3)   (%)  (%)  (%)  (%) ',
     >'-------------')
  158 FORMAT (5X,F5.2,F8.2,F10.2,F9.3,F8.0,2X,4F5.1,F6.1,F8.3)
  160 FORMAT (//,' SIMULATION BEGINS',/)
      END
C***********************************************************************
C
      SUBROUTINE DAYINP (JULIAN,YEAR,MAXJUL,INHOUR,NHRPDT,MTSTEP,INITAL,
     >    MPLTGRO,MWATRXT,NS,ALATUD,HRNOON,SUNHOR,TMPDAY,WINDAY,HUMDAY,
     >    PRECIP,SNODEN,SOITMP,VLCDAY,SOILXT,NPLANT,PLTHGT,DCHAR,PLTWGT,
     >    PLTLAI,ROOTDP)
C
C     THIS SUBROUTINE READS THE TIME VARYING DATA REQUIRED TO RUN THE
C     PROGRAM, AND STORES IT IN ARRAYS WITH VALUES FOR THE END OF EACH
C     HOUR.  DATA SUCH AS SOIL TEMPERATURE AND MOISTURE CONTENT AT
C     DEPTH, WHICH ARE SELDOM AVAILABLE AT EVERY HOUR, ARE
C     INTERPOLATED TO GET HOURLY VALUES.
C     IN:   JULIAN,YEAR,MAXJUL,INHOUR,NHRPDT,MTSTEP,INITAL,
C           MPLTGRO,MWATRXT,NS,ALATUD,HRNOON,NPLANT
C     INOUT:  SUNHOR,TMPDAY,WINDAY,HUMDAY,PRECIP,SNODEN,SOITMP,VLCDAY,SOILXT,PLTHGT,DCHAR,PLTWGT,
C           PLTLAI,ROOTDP
C***********************************************************************
      COMMON /MEASUR/ MEASDY,MEASHR,VLCMES(50),TSMEAS(50)
C
C     SUNHOR(I)   radiation from sun upon a horizontal surface for hour i (W/m^2)
C     TMPDAY(I)   air temperature at hour i (C)
C     WINDAY(I)   average windspeed over hour i (m/s)
C     PRECIP(I)   total precipitation occurring at end of hour i (m)
C     HUMDAY(I)   humidity at hour i
C     SNODEN(I)   density of newly fallen snow for hour i of current day (g/cm^3)
C     SOITMP(I)   temperature at lower boundary for hour i (C)
C     VLCDAY(I)   volumetric moisture content (ice+water) at hour i (m/m)
C     SOILXT(I,J) soil water extracted from soil layer i for hour j based on user input (m)

      REAL SUNHOR(24),TMPDAY(24),WINDAY(24),HUMDAY(24),PRECIP(24),
     >     SNODEN(24),SOITMP(24),VLCDAY(24),SOILXT(50,24)
      REAL VLC1(50),TMP1(50),FLUX(50)
      REAL PLTHGT(8),DCHAR(8),PLTWGT(8),PLTLAI(8),ROOTDP(8),
     >     ZCLST(8),DCHLST(8),WLST(8),TLALST(8),RDPLST(8),
     >     DAYLST(8),DAYSC(8)
      INTEGER IFLAGC(8),LCANDY(8),LCANYR(8)
      INTEGER YEAR
      SAVE NSTEP,NSTART,
     >     LTMPDY,LTMPHR,LTMPYR,TMPLST,TMP1,TMP,
     >     LVLCDY,LVLCHR,LVLCYR,VLCLST,VLC1,VLC,
     >     LWTRDY,LWTRHR,LWTRYR,LWTRMX,FLUX,
     >     LCANDY,LCANYR,ZCLST,DCHLST,WLST,TLALST,RDPLST,IFLAGC
C
      DATA  THRESH / 0.5/
C     THRESHHOLD WINDSPEED IS ASSUMED TO BE 0.5 M/S
C
C*********** READ WEATHER DATA**************
C
C     CHECK FLAG (MTSTEP) FOR TYPE OF WEATHER DATA:  0 = HOURLY;
C     1 = DAILY;   2 ==> DATA MATCHES NHRPDT
      IF (MTSTEP.EQ.1) THEN
         CALL DAY2HR (JULIAN,YEAR,INITAL,SUNHOR,TMPDAY,WINDAY,
     >                HUMDAY,PRECIP,SNODEN,ALATUD,HRNOON)
         NSTART=1
         NSTEP=1
         GO TO 25
      END IF

      IF (INITAL .EQ. 0) THEN
C        FIRST TIME INTO SUBROUTINE - FIND CORRECT PLACE IN DATA SET
         NSTEP=1
         IF (MTSTEP.EQ.2) NSTEP=NHRPDT
   10    READ(12,*,END=100) JD,JH,JYR,TMPDAY(NSTEP),WINDAY(NSTEP),
     >           HUMDAY(NSTEP),PRECIP(NSTEP),SNODEN(NSTEP),SUNHOR(NSTEP)
         IF (JH .EQ. 24) JD = JD + 1
         IF (JD.NE.JULIAN .OR. JYR.NE.YEAR .OR. JH.NE.NSTEP) GO TO 10
         NSTART=2*NSTEP
        ELSE
         NSTART=NSTEP
      END IF
C
      DO 20 I=NSTART,24,NSTEP
         READ (12,*) JD,JH,JYR,TMPDAY(I),WINDAY(I),HUMDAY(I),PRECIP(I),
     >              SNODEN(I),SUNHOR(I)
   20 CONTINUE
C
C     MAKE SURE IT IS AT THE RIGHT POINT IN THE DATA FILE
      IF (JD .NE. (JULIAN+1)  .OR.  JH .NE. 0) THEN
C        CHECK FOR EXCEPTIONS -> END OF YEAR OR CALLING HOUR 0, HOUR 24
         IF (JD .EQ. JULIAN  .AND.  JH .EQ. 24) GO TO 25
         IF (JULIAN.EQ.MAXJUL .AND. JD.EQ.1 .AND. JH.EQ.0) GO TO 25
C        NOT AT THE CORRECT POINT IN THE DATA FILE
         WRITE(6,*) ' ENCOUNTERED PROBLEMS READING HOURLY WEATHER DATA'
         WRITE (6,*) ' FOR JULIAN DAY ',JULIAN,' IN SUBROUTINE DAYINP'
         WRITE(21,*) ' ENCOUNTERED PROBLEMS READING HOURLY WEATHER DATA'
         WRITE (21,*) ' FOR JULIAN DAY ',JULIAN,' IN SUBROUTINE DAYINP'
         STOP
      END IF
C
   25 DO 27 I=NSTART,24,NSTEP
C        CONVERT HUMIDITY FROM % TO DECIMAL; PRECIP FROM INCHES TO
C        METERS; AND WIND FROM MILES PER HOUR TO M/S
         HUMDAY(I)=HUMDAY(I)/100.
         PRECIP(I)=PRECIP(I)*.0254
         WINDAY(I)=WINDAY(I)*0.447
C        SET WIND TO MINIMUM THRESHHOLD VALUE TO AVOID ZERO WINDSPEED
         IF (WINDAY(I) .LT. THRESH) WINDAY(I) = THRESH
   27 CONTINUE
C
C     AVERAGE HOURLY WEATHER DATA IF TIME STEPS ARE LARGER THAN 1 HOUR
      IF (NHRPDT.GT.1 .AND. NSTEP.EQ.1) THEN
         DO 35 I=NHRPDT,24,NHRPDT
            ISNOW=0
            TMPDAY(I)=TMPDAY(I)/NHRPDT
            WINDAY(I)=WINDAY(I)/NHRPDT
            HUMDAY(I)=HUMDAY(I)/NHRPDT
            IF (SNODEN(I)*PRECIP(I) .GT.0) THEN
               ISNOW=1
              ELSE
               SNODEN(I)=1.0
            END IF
            SNODEN(I)=SNODEN(I)*PRECIP(I)
            SUNHOR(I)=SUNHOR(I)/NHRPDT
C           SUM UP ALL HOURS WITHIN TIME STEP
            DO 30 J=I-NHRPDT+1,I-1
               TMPDAY(I)=TMPDAY(I)+TMPDAY(J)/NHRPDT
               WINDAY(I)=WINDAY(I)+WINDAY(J)/NHRPDT
               HUMDAY(I)=HUMDAY(I)+HUMDAY(J)/NHRPDT
               PRECIP(I)=PRECIP(I)+PRECIP(J)
               IF (SNODEN(J)*PRECIP(J) .GT.0) THEN
                  ISNOW=1
                 ELSE
                  SNODEN(J)=1.0
               END IF
               SNODEN(I)=SNODEN(I)+SNODEN(J)*PRECIP(J)
               SUNHOR(I)=SUNHOR(I)+SUNHOR(J)/NHRPDT
   30       CONTINUE
            IF (ISNOW.EQ.0 .OR. PRECIP(I).LT.0.0001) THEN
               SNODEN(I)=0.0
              ELSE
               SNODEN(I)=SNODEN(I)/PRECIP(I)
            END IF
   35    CONTINUE
      END IF
C
C**** READ TEMPERATURE DATA FOR LOWER BOUNDARY CONDITION
C
      IF (INITAL .EQ. 0) THEN
         LTMPDY=JULIAN
         LTMPHR=INHOUR
         LTMPYR=YEAR
         LVLCDY=JULIAN
         LVLCHR=INHOUR
         LVLCYR=YEAR
         TMPLST=SOITMP(INHOUR)
        ELSE
         TMPLST=SOITMP(24)
      END IF
      INHR=INHOUR
      NSTEP1=1
   40 IF (LTMPDY.EQ.JULIAN.AND.LTMPHR.EQ.INHR.AND.LTMPYR.EQ.YEAR) THEN
C        DAY AND HOUR ARE SAME AS THE LAST DAY AND HOUR READ FOR
C        TEMPERATURE DATA - READ TEMPERATUE DATA FOR NEXT SAMPLING DATE;
C        BUT FIRST SAVE THE LAST SAMPLING DATE SO IT CAN PRINTED
C        OUT ALONG WITH THE MOISTURE DATA IN SUBROUTINE 'OUTPUT'
         IF (LTMPDY .EQ. LVLCDY  .AND.  LTMPHR .EQ. LVLCHR  .AND. 
     >      INITAL .NE. 0) THEN
C           SAVE TEMPERATURE DATA
            DO 45 I=1,NS
               TSMEAS(I) = TMP1(I)
   45       CONTINUE
         END IF
         READ (14,*) LTMPDY,LTMPHR,LTMPYR,(TMP1(I), I=1,NS)
         IF (LTMPHR .EQ. 24) THEN
C           ALGORITHM ASSUMES HOUR 24 IS HOUR 0 -- ADJUST DAY AND HOUR
            LTMPHR=0
            LTMPDY=LTMPDY+1
            IF (LTMPDY .GT. MAXJUL) THEN
               LTMPDY=LTMPDY-MAXJUL
               LTMPYR=LTMPYR+1
            END IF
         END IF
         TMP = TMP1(NS)
      END IF
      DAYS=(LTMPYR-YEAR)*MAXJUL + LTMPDY - JULIAN
      LASTHR=LTMPHR
      IF (DAYS .GT. 0) THEN
         LASTHR=24
C        ALLOW FOR STEPPING MORE THAN 1 HOUR IF NO DATA ON THIS DAY
         IF (INHR.EQ.INHOUR) NSTEP1=NHRPDT
        ELSE
         IF (DAYS .LT. 0) THEN
            WRITE (6,110) LTMPDY,LTMPHR,LTMPYR
            WRITE (21,110) LTMPDY,LTMPHR,LTMPYR
            STOP
         END IF
      END IF
C     INTERPOLATE FOR EACH HOUR OF DAY BETWEEN THE CURRENT HOUR AND
C     THE LAST SAMPLING DATE READ.
      J=LASTHR
      DO 50 I=INHR+NSTEP1,LASTHR,NSTEP1
         J=I
         SOITMP(I) = TMPLST
     >             +(I-INHR)*(TMP-TMPLST)/(DAYS*24+LTMPHR-INHR)
   50 CONTINUE
      INHR=J
      TMPLST=SOITMP(INHR)
      IF (LASTHR .NE. 24) THEN
C        HAVE NOT REACH END OF CURRENT DAY - GO BACK AND READ DATA AGAIN
         GO TO 40
      END IF
C     AVERAGE HOURLY VALUES IF TIME STEPS ARE LARGER THAN 1 HOUR AND
C     DATA WERE GIVEN DURING MIDDLE OF CURRENT DAY
      IF (NHRPDT.GT.1 .AND. NSTEP1.EQ.1) THEN
         DO 55 I=INHOUR+NHRPDT,24,NHRPDT
            SOITMP(I)=SOITMP(I)/NHRPDT
            DO 54 J=I-NHRPDT+1,I-1
               SOITMP(I)=SOITMP(I)+SOITMP(J)/NHRPDT
   54       CONTINUE
   55    CONTINUE
      END IF
C
C**** READ MOISTURE CONTENT DATA FOR THE LOWER BOUNDARY CONDITION
C
      IF (INITAL .EQ. 0) THEN
         VLCLST=VLCDAY(INHOUR)
        ELSE
         VLCLST=VLCDAY(24)
      END IF
      INHR=INHOUR
      NSTEP2=1
   60 IF (LVLCDY.EQ.JULIAN.AND.LVLCHR.EQ.INHR.AND.LVLCYR.EQ.YEAR) THEN
C        DAY AND HOUR ARE SAME AS THE LAST DAY AND HOUR READ FOR
C        MOISTURE DATA - READ MOISTURE DATA FOR NEXT SAMPLING DATE;
C        BUT FIRST SAVE THE DATA FOR THE LAST SAMPLING DATE READ SO
C        SUBROUTINE 'OUTPUT' CAN PRINT IT OUT WITH PREDICTED VALUES
         IF (INITAL .NE. 0) THEN
C           IF INITAL =0, VLCMES ALREADY SET IN 'INPUT'
            MEASHR = LVLCHR
            MEASDY = LVLCDY
            DO 65 I=1,NS
               VLCMES(I) = VLC1(I)
   65       CONTINUE
         END IF
         READ (13,*) LVLCDY,LVLCHR,LVLCYR,(VLC1(I), I=1,NS)
         IF (LVLCHR .EQ. 24) THEN
C           ALGORITHM ASSUMES HOUR 24 IS HOUR 0 -- ADJUST DAY AND HOUR
            LVLCHR=0
            LVLCDY=LVLCDY+1
            IF (LVLCDY .GT. MAXJUL) THEN
               LVLCDY=LVLCDY-MAXJUL
               LVLCYR=LVLCYR+1
            END IF
         END IF
         VLC = VLC1(NS)
      END IF
      DAYS=(LVLCYR-YEAR)*MAXJUL + LVLCDY - JULIAN
      LASTHR=LVLCHR
      IF (DAYS .GT. 0) THEN
         LASTHR=24
C        ALLOW FOR STEPPING MORE THAN 1 HOUR IF NO DATA ON THIS DAY
         IF (INHR.EQ.INHOUR) NSTEP2=NHRPDT
        ELSE
         IF (DAYS .LT. 0) THEN
            WRITE (6,120) LVLCDY,LVLCHR,LVLCYR
            WRITE (21,120) LVLCDY,LVLCHR,LVLCYR
            STOP
         END IF
      END IF
C     INTERPOLATE FOR EACH HOUR OF DAY BETWEEN THE CURRENT HOUR AND
C     THE LAST SAMPLING DATE READ.
      J=LASTHR
      DO 70 I=INHR+NSTEP2,LASTHR,NSTEP2
         J=I
         VLCDAY(I) = VLCLST
     >             +(I-INHR)*(VLC-VLCLST)/(DAYS*24+LVLCHR-INHR)
   70 CONTINUE
      INHR=J
      VLCLST=VLCDAY(INHR)
      IF (LASTHR .NE. 24) THEN
C        HAVE NOT REACH END OF CURRENT DAY - GO BACK AND READ DATA AGAIN
         GO TO 60
      END IF
C     AVERAGE HOURLY VALUES IF TIME STEPS ARE LARGER THAN 1 HOUR AND
C     DATA WERE GIVEN DURING MIDDLE OF CURRENT DAY
      IF (NHRPDT.GT.1 .AND. NSTEP2.EQ.1) THEN
         DO 75 I=INHOUR+NHRPDT,24,NHRPDT
            VLCDAY(I)=VLCDAY(I)/NHRPDT
            DO 74 J=I-NHRPDT+1,I-1
               VLCDAY(I)=VLCDAY(I)+VLCDAY(J)/NHRPDT
   74       CONTINUE
   75    CONTINUE
      END IF
C
C
C     READ SOIL SINK DATA
      IF (MWATRXT .EQ. 0) GO TO 200
      IF (INITAL .EQ. 0) THEN
C        FIRST TIME INTO SUBROUTINE - FIND CORRECT PLACE IN DATA SET
         LWTRDY=JULIAN
         LWTRHR=INHOUR
         LWTRYR=YEAR
         LWTRMX=MAXJUL
         JDY1=JULIAN
         JHR1=INHOUR
         JYR1=YEAR
   80    READ (15,*,END=102) JDY,JHR,JYR,(FLUX(J),J=1,NS)
         DAYS=(JYR-LWTRYR)*LWTRMX + JDY - LWTRDY + (JHR-LWTRHR)/24.
         IF (DAYS .LE. 0) THEN
            JYR1=JYR
            JDY1=JDY
            JHR1=JHR
            GO TO 80
         END IF
         DAYS=(JYR-JYR1)*LWTRMX + JDY - JDY1 + (JHR-JHR1)/24.
         DO 82 J=1,NS
            FLUX(J)=FLUX(J)/DAYS/24./3600.
   82    CONTINUE
         LWTRDY=JDY
         LWTRHR=JHR
         LWTRYR=JYR
      END IF
      DO 90 I=INHOUR,24,NHRPDT
   84    DAYS=(YEAR-LWTRYR)*LWTRMX + JULIAN - LWTRDY + (I-LWTRHR)/24.
         IF (DAYS .GT. 0.0) THEN
C           TIME IS NOW PAST LAST LINE OF DATA READ -- READ DATA AGAIN
            READ (15,*) JDY,JHR,JYR,(FLUX(J),J=1,NS)
            DAYS=(JYR-LWTRYR)*LWTRMX + JDY - LWTRDY + (JHR-LWTRHR)/24.
            DO 85 J=1,NS
               FLUX(J)=FLUX(J)/DAYS/24./3600.
   85       CONTINUE
            LWTRDY=JDY
            LWTRHR=JHR
            LWTRYR=JYR
            LWTRMX=MAXJUL
            GO TO 84
         END IF
         DO 88 J=1,NS
            SOILXT(J,I)=FLUX(J)
   88    CONTINUE
   90 CONTINUE
C
C
C***  READ CANOPY INFORMATION
C
  200 IF (NPLANT .EQ. 0 .OR. MPLTGRO .EQ. 0) GO TO 250
      IF (INITAL .EQ. 0) THEN
C        FIRST TIME INTO SUBROUTINE - FIND CORRECT PLACE IN DATA SET
         DO 212 J=1,NPLANT
            IFLAGC(J)=-1
            DAYLST(J)=0.0
  210       READ (40+J,*,END=104) LCANDY(J),LCANYR(J),
     >            ZCLST(J),DCHLST(J),WLST(J),TLALST(J),RDPLST(J)
C           CONVERT LEAF DIMENSION FROM CM TO METERS
            DCHLST(J)=DCHLST(J)/100.
            DAYSC(J)=(LCANYR(J)-YEAR)*MAXJUL + LCANDY(J) - JULIAN
            IF (DAYSC(J) .LE. 0) THEN
C              SAVE LAST VALUES READ IN
               DAYLST(J)=-DAYSC(J)
               PLTHGT(J) = ZCLST(J)
               DCHAR(J) = DCHLST(J)
               PLTWGT(J) = WLST(J)
               PLTLAI(J) = TLALST(J)
               ROOTDP(J) = RDPLST(J)
               IFLAGC(J)=0
               GO TO 210
             ELSE
C              INTERPOLATE VEGETATION DATA FOR STARTING DAY
C              (CHECK IF DATA START ON OR BEFORE STARTING DAY
               IF (IFLAGC(J).LT.0) GO TO 104
               IFLAGC(J)=0
               RATIO=DAYLST(J)/(DAYSC(J)+DAYLST(J))
               PLTHGT(J)=PLTHGT(J)+(ZCLST(J)-PLTHGT(J))*RATIO
               DCHAR(J)=DCHAR(J)+(DCHLST(J)-DCHAR(J))*RATIO
               PLTWGT(J)=PLTWGT(J)+(WLST(J)-PLTWGT(J))*RATIO
               PLTLAI(J)=PLTLAI(J)+(TLALST(J)-PLTLAI(J))*RATIO
               ROOTDP(J)=ROOTDP(J)+(RDPLST(J)-ROOTDP(J))*RATIO
            END IF
  212    CONTINUE
        ELSE
C
         DO 215 J=1,NPLANT
C           CHECK IF WE HAVE REACHED END OF VEGETATION DATA FILE
            IF (IFLAGC(J) .EQ. 1) GO TO 215
C
C           INTERPOLATE BETWEEN CURRENT DAY AND NEXT DAY W/ PLANT DATA -
C           CALCULATE # OF DAYS TO NEXT DAY WITH PLANT DATA
            DAYSC(J)=(LCANYR(J)-YEAR)*MAXJUL + LCANDY(J) - JULIAN
C
C           CHECK TIME COMPARED TO DATA -- INCLUDE COMPARISON OF YEAR 
C           FOR SITUATION WHEN INPUT DATA INCLUDE DAYS 366 AND DAY 1 AT
C           THE END OF A LEAP YEAR.  (MAXJUL PERTAINS TO YEAR RATHER 
C           THAN LCANYR)
            IF (DAYSC(J).LT.0.0 .OR. LCANYR(J).LT.YEAR) THEN
C              TIME IS NOW PAST LAST LINE OF DATA READ - READ DATA AGAIN
               READ (40+J,*,END=214) LCANDY(J),LCANYR(J),
     >              ZCLST(J),DCHLST(J),WLST(J),TLALST(J),RDPLST(J)
C              CONVERT LEAF DIMENSION FROM CM TO METERS
               DCHLST(J)=DCHLST(J)/100.
               DAYSC(J)=(LCANYR(J)-YEAR)*MAXJUL + LCANDY(J) - JULIAN
            END IF
            PLTHGT(J) = PLTHGT(J) + (ZCLST(J)-PLTHGT(J))/(DAYSC(J)+1)
            DCHAR(J)= DCHAR(J) + (DCHLST(J)-DCHAR(J))/(DAYSC(J)+1)
            PLTWGT(J) = PLTWGT(J) + (WLST(J)-PLTWGT(J))/(DAYSC(J)+1)
            PLTLAI(J) = PLTLAI(J) + (TLALST(J)-PLTLAI(J))/(DAYSC(J)+1)
            ROOTDP(J) = ROOTDP(J) + (RDPLST(J)-ROOTDP(J))/(DAYSC(J)+1)
            GO TO 215
C
C           END OF VEGETATION DATA FILE -- ISSUE WARNING AND SET FLAG
  214       WRITE (21,140) J,JULIAN,YEAR
            WRITE (6,140) J,JULIAN,YEAR
            IFLAGC(J)=1
C
  215    CONTINUE
      END IF
C
  250 RETURN
C
  100 WRITE (6,*) ' CANNOT FIND WEATHER DATA FOR BEGINNING DAY'
      WRITE (21,*) ' CANNOT FIND WEATHER DATA FOR BEGINNING DAY'
      STOP
  102 WRITE (6,*) ' CANNOT FIND SOIL SINK DATA FOR BEGINNING DAY'
      WRITE (21,*) ' CANNOT FIND SOIL SINK DATA FOR BEGINNING DAY'
      STOP
  104 WRITE (6,105) J
      WRITE (21,105) J
  105 FORMAT (' INITIAL CONDITIONS FOR PLANT #',I2,
     >' CANNOT BE INTERPOLATED FROM PLANT GROWTH FILE')
      STOP
  110 FORMAT(' PROBLEM ENCOUNTERED IN SUBROUTINE DAYINP:',/,
     >' SOIL TEMPERATURE DATA IS NOT CHRONOLOGICAL',/,
     >' AT JULIAN DAY ',I4,', HOUR ',I4,', YEAR',I6)
  120 FORMAT(' PROBLEM ENCOUNTERED IN SUBROUTINE DAYINP:',/,
     >' SOIL WATER CONTENT DATA IS NOT CHRONOLOGICAL',/,
     >' AT JULIAN DAY ',I4,', HOUR ',I4,', YEAR',I6)
  140 FORMAT(' *** LAST DATA FOR GROWTH OF PLANT #',I1,'WAS ON DAY',
     >I4,' OF YEAR ',I4,' ***'/)
      END
C***********************************************************************
C
      SUBROUTINE DAY2HR (JULIAN,YEAR,INITAL,SUNHOR,TMPDAY,WINDAY,
     >                   HUMDAY,PRECIP,SNODEN,ALATUD,HRNOON)
C
C     THIS SUBROUTINE SEPARATES DAILY VALUES FOR MAXIMUM-MINIMUM AIR
C     TEMPERATURE (C), DEW-POINT TEMPERATURE (C), WIND RUN (MILES),
C     PRECIPITATION (INCHES) AND SOLAR RADIATION (W/M2) INTO HOURLY
C     VALUES OF TEMPERATURE, WINDSPEED, RELATIVE HUMIDITY, PRECIPITATION
C     AND SOLAR RADIATION.
C
C     THE INPUT FILE IS READ WITH OPEN FORMAT, AND DATA FOR EACH DAY IN
C     THE INPUT FILE MUST BE IN THE FOLLOWING ORDER:
C
C     DAY, YEAR, MAX TEMP, MIN TEMP, DEW-POINT, WIND, PRECIP, SOLAR
C
C     For example:
C
C      65  87   12.3    4.6    2.6   109.0   0.09    38.9
C      66  87   13.2    1.6    0.6    77.6   0.00   171.1
C      67  87   15.8    1.4    1.4   125.0   0.06   168.8
C      68  87   12.6    5.0    3.0   119.2   0.08   110.4
C***********************************************************************
C
      COMMON /SWRCOE/ SOLCON,DIFATM,DIFRES,SNOCOF,SNOEXP
C
      REAL SUNHOR(24),TMPDAY(24),HUMDAY(24),WINDAY(24),
     >     PRECIP(24),DYPREC,SNODEN(24)
      INTEGER HOUR,YEAR,YEAR3
      SAVE TMAX1,TMAX,TMIN,TDEW,WIND,DYPREC,SOLRAD,JDAY,JYR
      DATA PI/3.14159/ DTMIN/0.5/
C
      IF (INITAL.EQ.0) THEN
         TMAX=-999.
C        SAVE MAXIMUM TEMPERATURE FROM PREVIOUS DAY TO ALLOW FOR
C        INTERPOLATION OF MORNING TEMPERATURES
   10    TMAX1=TMAX
         READ (12,*,END=50) JDAY,JYR,TMAX,TMIN,TDEW,WIND,DYPREC,SOLRAD
         IF (JDAY.NE.JULIAN .OR. JYR.NE.YEAR) GO TO 10
C        CHECK IF TMAX1 IS SET TO PREVIOUS DAY'S MAX TEMP (IF START TIME
C        IS FIRST LINE OF DATA IN FILE, SET TMAX1 TO CURRENT DAY'S TEMP)
         IF (TMAX1.LT.-998.) TMAX1=TMAX
      END IF
C
      READ (12,*) JULAN3,YEAR3,TMAX3,TMIN3,TDEW3,WIND3,PREC3,SOLAR3
C
C     MAKE SURE IT IS AT THE RIGHT POINT IN THE DATA FILE
      IF (JDAY.NE.JULIAN .OR. JYR.NE.YEAR) THEN
C        NOT AT THE CORRECT POINT IN THE DATA FILE
         WRITE(6,*) ' ENCOUNTERED PROBLEMS READING DAILY WEATHER DATA'
         WRITE (6,*) ' AT JULIAN DAY ',JULIAN,' IN SUBROUTINE DAY2HR'
         WRITE(21,*) ' ENCOUNTERED PROBLEMS READING DAILY WEATHER DATA'
         WRITE (21,*) ' AT JULIAN DAY ',JULIAN,' IN SUBROUTINE DAY2HR'
         STOP
      END IF
C
C**** CALCULATE SUN'S DECLINATION ANGLE (DECLIN), HALF-DAY ANGLE
C     (HAFDAY), MAXIMUM SOLAR RADIATION AT 100% TRANMISSIVITY (SUNMAX),
C     TOTAL TRANSMISSIVITY TO RADIATION (TTOTAL), AND HOUR OF SUNRISE
C     AND SUNSET
      DECLIN = 0.4102*SIN(2*PI*(JULIAN-80)/365.)
      COSHAF=-TAN(ALATUD)*TAN(DECLIN)
      IF (COSHAF .GE. 1.0) THEN
         HAFDAY=PI
       ELSE
         HAFDAY=ACOS(COSHAF)
      END IF
      SUNMAX=SOLCON*(HAFDAY*SIN(ALATUD)*SIN(DECLIN)
     >                  +COS(ALATUD)*COS(DECLIN)*SIN(HAFDAY))/3.14159
      TTOTAL=SOLRAD/SUNMAX
      SUNRIS=HRNOON - HAFDAY/0.261799
      SUNSET=HRNOON + HAFDAY/0.261799
C
      TIMMIN = SUNRIS - DTMIN
      TIMMAX = (SUNSET + HRNOON)/2.
      DAY = TIMMAX - TIMMIN
      ANIGHT = 24. - DAY
      CALL VSLOPE (DUMMY,VAPOR,TDEW)

      DO 30 HOUR=1,24
C**** SET HOUR TO MID-POINT OF THE HOUR, OR IF THE SUN RISES OR SETS IN
C     THIS HOUR, SET HOUR TO MID-POINT OF THE HOUR AND TIME AT WHICH IT
C     RISES OR SETS.
      HALFHR=HOUR-0.5
      IF (HOUR .GT. SUNRIS  .AND.  (HOUR-1.0) .LT. SUNRIS)
     >   HALFHR=(SUNRIS + HOUR)/2.
      IF (HOUR .GT. SUNSET  .AND.  (HOUR-1.0) .LT. SUNSET)
     >   HALFHR=SUNSET - (SUNSET-(HOUR-1))/2.
C
C**** DETERMINE THE GEOMETRY OF THE SUN ANGLE AT THE HOUR MID=POINT
      HRANGL=0.261799*(HALFHR-HRNOON)
      ALTITU=ASIN( SIN(ALATUD)*SIN(DECLIN)
     >           + COS(ALATUD)*COS(DECLIN)*COS(HRANGL))
      SUNMAX=SOLCON*SIN(ALTITU)
      SUNHOR(HOUR) = TTOTAL*SUNMAX
      IF (SUNHOR(HOUR) .LT. 0) SUNHOR(HOUR) = 0.0
C
C
C**** INTERPOLATE TEMPERATURE BETWEEN MAXIMUM AND MINIMUM TEMPERATURE
C     DEPENDING ON TIME OF DAY
         IF (HOUR .LT. TIMMIN) THEN
            TMPDAY(HOUR) = 0.5*((TMAX1+TMIN)
     >                   + (TMAX1-TMIN)*COS(PI/ANIGHT*(24+HOUR-TIMMAX)))
            GO TO 25
         END IF
         IF (HOUR .GT. TIMMAX) THEN
            TMPDAY(HOUR) = 0.5*((TMAX+TMIN3)
     >                   + (TMAX-TMIN3)*COS(PI/ANIGHT*(HOUR-TIMMAX)))
            GO TO 25
         END IF
         TMPDAY(HOUR) = 0.5*((TMAX+TMIN)
     >                + (TMAX-TMIN)*COS(PI/DAY*(HOUR-TIMMAX)))
C
C****    CALCULATE HUMIDITY
   25    CALL VSLOPE (DUMMY,SATV,TMPDAY(HOUR))
         HUMDAY(HOUR) = 100.*VAPOR/SATV
         IF (HUMDAY(HOUR) .GT. 100.) HUMDAY(HOUR) = 100.
C
C****    CALCULATE HOURLY WINDSPEED
         WINDAY(HOUR) = WIND/24.0
C
C****    ASSIGN DAILY PRECIPITATION TO THE FIRST HOUR
         IF (HOUR.EQ.1) THEN
            PRECIP(HOUR) = DYPREC
            SNODEN(HOUR) = 0.0
           ELSE
            PRECIP(HOUR) = 0.0
            SNODEN(HOUR) = 0.0
         END IF
C
   30 CONTINUE
C
      TMAX1 = TMAX
      TMAX = TMAX3
      TMIN = TMIN3
      TDEW = TDEW3
      WIND = WIND3
      DYPREC = PREC3
      SOLRAD = SOLAR3
      JDAY = JULAN3
      JYR = YEAR3
C
      RETURN
C
C**** CANNOT FIND STARTING POINT IN DATA FILE
   50 WRITE (6,120)
      WRITE (21,120)
  120 FORMAT (' ERROR READING DAILY WEATHER DATA IN SUBROUTINE DAY2HR:',
     >/,' YEAR AND/OR DAY DO NOT MATCH')
      STOP
      END
C***********************************************************************
C
      SUBROUTINE CLOUDY (CLOUDS,ALATUD,DECLIN,HAFDAY,SUNHOR,JULIAN,
     >  NHRPDT)
C
C     THIS SUBROUTINE TOTALS THE SOLAR RADIATION FOR THE DAY.  IT THEN
C     ESTIMATES THE CLOUD COVER, WHICH WILL BE USED IN THE CALCULATION
C     FOR THE EMMISSIVITY OF THE ATMOSPHERE.  (DECLINATION AND HALF-DAY
C     LENGTH FOR THE CURRENT JULIAN DAY ARE ALSO CALCULATED.)
C
C***********************************************************************
      COMMON /SWRCOE/ SOLCON,DIFATM,DIFRES,SNOCOF,SNOEXP
C
      REAL SUNHOR(24)
C
      TOTSUN=0.0
      DO 10, I=NHRPDT,24
         IF (SUNHOR(I) .LT. 0.0) SUNHOR(I)=0.0
         TOTSUN=TOTSUN + SUNHOR(I)
   10 CONTINUE
      DECLIN=0.4102*SIN(2*3.14159*(JULIAN-80)/365.)
      COSHAF=-TAN(ALATUD)*TAN(DECLIN)
      IF (ABS(COSHAF) .GE. 1.0) THEN
         IF (COSHAF .GE. 1.0) THEN
C           SUN DOES NOT COME UP ON THIS DAY (WINTER IN ARCTIC CIRCLE)
            HAFDAY=0.0
          ELSE
C           SUN DOES NOT SET ON THIS DAY (SUMMER IN THE ARCTIC CIRCLE)
            HAFDAY=2*3.14159
         END IF
       ELSE
         HAFDAY=ACOS(COSHAF)
      END IF
      SUNMAX=24.*SOLCON*(HAFDAY*SIN(ALATUD)*SIN(DECLIN)
     >                  +COS(ALATUD)*COS(DECLIN)*SIN(HAFDAY))/3.14159
      TTOTAL=TOTSUN/SUNMAX
      CLOUDS=2.4 - 4.*TTOTAL
      IF (CLOUDS .GT. 1.0) CLOUDS=1.0
      IF (CLOUDS .LT. 0.0) CLOUDS=0.0
C
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE CANOPY (NPLANT,NC,NSP,NSPLST,MZCINP,ITYPE,INITAL,
     >   PLTHGT,PLTWGT,PLTLAI,RLEAF0,ZSP,ZC,TC,TCDT,VAPC,VAPCDT,
     >   WCAN,WCANDT,WCMAX,PCANDT,CANMA,CANMB,TA,HUMID)
C
C     THIS SUBROUTINE DETERMINES WHICH LAYERS OF THE CANOPY ARE COVERED
C     WITH SNOW
C
C***********************************************************************
      COMMON /CLAYRS/ DRYCAN(8,10),CANLAI(8,10),TOTLAI(8),IEVAP(8),
     >                RLEAF(8,10),RROOT(8,50),ROOTDN(8,50),TOTROT(8)
      REAL PLTHGT(8),PLTWGT(8),PLTLAI(8),ZSP(100),RLEAF0(8),ZC(11),
     >     TC(11),TCDT(11),VAPC(11),VAPCDT(11),WCAN(10),WCANDT(10),
     >     PCANDT(8)
      REAL ZZC(11),PLANTZ(8),PLANTW(8)
      INTEGER ITYPE(8)
      SAVE NNC,ZZC
C
      IF (INITAL .EQ. 0) THEN
         IF (MZCINP .EQ. 1) THEN
C           LAYERING OF CANOPY NODES SPECIFIED BY USER -- SAVE SPACING
            NNC=NC
            ZZC(1)=0.0
            DO 5 I=2,NC+1
               ZZC(I)=ZC(I)
               TCDT(I)=TCDT(I-1)
               VAPCDT(I)=VAPCDT(I-1)
    5       CONTINUE
C           IF NO NEED TO ADJUST FOR SNOWPACK, RETURN TO MAIN PROGRAM
            IF (NSP .EQ. 0) GO TO 80
         END IF
      END IF
C
C     CHECK IF CANOPY HAS EMERGED
      NCCHK=0
      DO 8 J=1,NPLANT
         IF (PLTLAI(J).GT.0) NCCHK=1
    8 CONTINUE
      IF (NCCHK.EQ.0) THEN
         NC=0
         GO TO 80
      END IF
C
C     CANOPY IS PRESENT - CHECK IF COVERED BY SNOW
      NCLAST=NC
C
C     NSPLST       number of snow pack nodes at end of last time step
      IF (NCLAST .EQ. 0) THEN
C        NO CANOPY LAST TIME STEP
C        DEFINE STATE VARIABLES FOR TOP OF CANOPY
         TC(1)=TA
         TCDT(1)=TA
         CALL VSLOPE (DUMMY,SATV,TA)
         VAPC(1)=HUMID*SATV
         VAPCDT(1)=HUMID*SATV
         IF (NSPLST.GT.0) THEN
C           SNOW MELTED AND EXPOSED TOP OF CANOPY--SET MAX WATER CONTENT
            WCAN(1)=WCMAX
            WCANDT(1)=WCMAX
           ELSE
C           SET WATER CONTENT TO EQUILIBRIUM WITH HUMIDITY
            CALL CANHUM (2,HUMID,DUMMY,WCANDT(1),TCDT(1),CANMA,CANMB)
            WCAN(1)=WCANDT(1)
         END IF
         NC=1
      END IF
C
      IF (NSP .GT. 0) THEN
C        DETERMINE IF CANOPY IS COVERED WITH SNOW
         DO 10 J=1,NPLANT
            IF (PLTHGT(J).GT.ZSP(NSP+1) .AND. PLTLAI(J).GT.0.0) GO TO 15
   10    CONTINUE
C        SNOW COVERS CANOPY
         NC=0
         GO TO 80
C
   15    IF (MZCINP .EQ. 1) GO TO 25
C        ALLOW MODEL TO REDEFINE LAYERING AND NODES WITHIN CANOPY
         DO 20 J=1,NPLANT
            PLANTZ(J)=PLTHGT(J)-ZSP(NSP+1)
            IF (PLANTZ(J).LT.0.0) THEN
               PLANTZ(J)=0.0
               PLANTW(J)=0.0
               TOTLAI(J)=0.0
              ELSE
               PLANTW(J)=PLTWGT(J)*PLANTZ(J)/PLTHGT(J)
               TOTLAI(J)=PLTLAI(J)*PLANTZ(J)/PLTHGT(J)
            END IF
   20    CONTINUE
C
C        DETERMINE LAYERING OF CANOPY ACCOUNTING FOR SNOW
         CALL CANLAY (NC,NPLANT,ITYPE,ZC,WCANDT,TCDT,VAPCDT,PLANTZ,
     >                PLANTW,TOTLAI,RLEAF0)
         GO TO 55
C
C        CANOPY LAYERING WAS INPUT - DO NOT ALLOW MOVING OF NODES
   25    DO 30 I=NNC,1,-1
            IF ((ZZC(NNC+1)-ZZC(I)) .GT. ZSP(NSP+1)) THEN
               NC=I
               ZC(NC+1) = ZZC(NNC+1)-ZSP(NSP+1)
               IF (NC .EQ. 1) THEN
                  ZC(NC)=0.0
                  GO TO 35
               END IF
               ZMID1=(ZZC(NC)+ZZC(NC-1))/2
               ZMID2=(ZZC(NC+1)+ZZC(NC))/2
               IF (ZC(NC+1) .GT. ZMID2) THEN
                  ZC(NC)=ZZC(NC)
                 ELSE
                  ZC(NC)=ZZC(NC) - (ZZC(NC)-ZMID1)*(ZMID2-ZC(NC+1))
     >                             /(ZMID2-ZZC(NC))
               END IF
               GO TO 35
            END IF
   30    CONTINUE
        ELSE
C
C        NO SNOW ON GROUND
         IF (MZCINP .EQ. 1) THEN
C           RESET NUMBER OF NODES AND POSITION OF BOTTOM TWO NODES
C           (POSITION OF OTHER NODES RESET BELOW)
            NC=NNC
            ZC(NC)=ZZC(NC)
            ZC(NC+1)=ZZC(NC+1)
           ELSE
C           DETERMINE LAYERING OF CANOPY WITH NO SNOWPACK (ADJUSTING FOR
C           ANY CHANGES IN CANOPY CHARACTERISTICS)
            CALL CANLAY (NC,NPLANT,ITYPE,ZC,WCANDT,TCDT,VAPCDT,PLTHGT,
     >                   PLTWGT,PLTLAI,RLEAF0)
         END IF
      END IF
C
   35 IF (MZCINP.EQ.1) THEN
C        CALCULATE TOTAL LEAF AREA INDEX ABOVE SNOW FOR EACH PLANT TYPE
         IF (NC .NE. NCLAST) THEN
            DO 45 J=1,NPLANT
               TOTLAI(J)=0.0
               DO 40 I=1,NC
                  TOTLAI(J)=TOTLAI(J)+CANLAI(J,I)
   40          CONTINUE
   45       CONTINUE
C
            IF (NC .GT. NCLAST) THEN
C              SNOW HAS MELTED AND EXPOSED ADDITIONAL CANOPY LAYERS --
C              DEFINE STATE VARIABLES OF NEWLY EXPOSED LAYERS
               NC1=NCLAST
               IF (NCLAST .EQ. 0) NC1=1
               DO 50 I=NC1+1,NC
                  TCDT(I)=TCDT(I-1)
                  VAPCDT(I)=VAPCDT(I-1)
                  ZC(I-1)=ZZC(I-1)
   50          CONTINUE
            END IF
C
         END IF
      END IF
   55 IF (NC .GT. NCLAST) THEN
C        SET WATER CONTENT, TEMPERATURE AND VAPOR OF NEW LAYERS
C        FOR BEGINNING OF TIME STEP
         NC1=NCLAST
         IF (NCLAST .EQ. 0) NC1=1
         DO 60 I=NC1+1,NC
            IF (NSPLST .GT. 0) THEN
C              SNOW HAS MELTED EXPOSING LAYERS - WATER CONTENT = MAXIMUM
               WCANDT(I)=WCMAX
              ELSE
               WCANDT(I)=WCANDT(I-1)
            END IF
            WCAN(I)=WCANDT(I)
            TC(I)=TCDT(I)
            VAPC(I)=VAPCDT(I)
   60    CONTINUE
      END IF
C
   80 DO 85 J=1,NPLANT
         IF (PLTLAI(J).LE. 0.0 .OR. TOTLAI(J) .LE. 0.0) PCANDT(J)=0.0
   85 CONTINUE
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE CANLAY (NC,NPLANT,ITYPE,ZC,WCANDT,TCDT,VAPCDT,PLANTZ,
     >                   PLANTW,PLTLAI,RLEAF0)
C
C      This subroutine splits the canopy into layers.  Each layer
C      will have the same total leaf area index, but the actual layer
C      dimension will vary.  The different plant variables are
C      dimensioned by plant type and layer, and their values are
C      apportioned according to their representation in the layer.
C
C***********************************************************************
      COMMON /CLAYRS/ DRYCAN(8,10),CANLAI(8,10),TOTLAI(8),IEVAP(8),
     >                RLEAF(8,10),RROOT(8,50),ROOTDN(8,50),TOTROT(8)
C
      REAL ZC(11),WCANDT(10),TCDT(11),VAPCDT(11),
     >     PLANTZ(8),PLANTW(8),PLTLAI(8),RLEAF0(8)
      REAL LAISUM,LAILYR,LYRTOP,LYRTOT,DZMAX,TMP(10),A,DZ(10),
     >     FACTOR,TMPLAI
C
      INTEGER NC,HTINDX,HTORDR(8),J,LAYER,OLDNC,ITYPE(8)
C
C     -- Local variables:
C     A         Leaf area index per layer.
C     DZMAX     Maximum theoretical layer thickness for current values.
C     HTINDX    Index for ordered plant height - in increasing order.
C     HTORDR    Array of plant numbers ordered by height (1 = smallest)
C     LAILYR    Current accumulated LAI for current layer.
C     LAISUM    The delta LAI/unit distance for current calculations.
C     LAYER     Current layer number, used for loop control.
C     LYRTOP    Height of top of current layer.
C     LYRTOT    Current accumulated layer thickness for current layer.
C     TMPLAI    Sum of total leaf area index for all plants.

C**   --Initialization
      LAILYR = 0.0
      LYRTOT = 0.0
      LAISUM = 0.0
      TMPLAI = 0.0
      NZEROS = 0
      LYRTOP = 0.0
      OLDNC = NC

C     -- Initialize temporary variable used for sorting,
C        sum of total leaf area index for all plants,
C        and sum of delta LAI per unit distance.
      DO 10 J=1, NPLANT
         IF (PLANTZ(J)*PLTLAI(J) .NE. 0.0) THEN
            TMP(J) = PLANTZ(J)
            TOTLAI(J) = PLTLAI(J)
            TMPLAI = TMPLAI + TOTLAI(J)
            LAISUM = LAISUM + TOTLAI(J) / PLANTZ(J)
           ELSE
            TMP(J) = 9999.
            NZEROS=NZEROS+1
            HTORDR(NZEROS)=J
            TOTLAI(J) = 0.0
         END IF
   10 CONTINUE

C**   Arrange the plants by height - in increasing order
      DO 40 HTINDX = NZEROS+1, NPLANT
         HTORDR(HTINDX) = 1
         DO 20 J = 1, NPLANT
            IF (TMP(HTORDR(HTINDX)).GT.TMP(J)) HTORDR(HTINDX) = J
   20    CONTINUE
         TMP(HTORDR(HTINDX)) = 9999.0
   40 CONTINUE
C

      NC = NINT (TMPLAI / 0.75)
C     -- Set upper limit for number of canopy layers.
      IF (NC.EQ.0) THEN
         NC=1
        ELSE
         IF (NC.GT.10) NC = 10
      END IF
C
      HTINDX = NZEROS+1
C
      A  = TMPLAI / NC
      LAYER  = NC
      ZC(NC+1) = PLANTZ(HTORDR(NPLANT))


C**   --Main loop   ***

C     --Calculate theoretical maximum layer based on current LAISUM.
  100 IF (LAYER .NE. 1) THEN
         DZMAX = (A - LAILYR) / LAISUM
        ELSE
         DZMAX = PLANTZ(HTORDR(NPLANT)) - LYRTOP
      END IF
C
      IF(LYRTOP+DZMAX.LE.PLANTZ(HTORDR(HTINDX)).OR.LAYER.EQ.1)THEN
C**   -- Top of layer is below or equal to next tallest plant.

         LYRTOP = LYRTOP + DZMAX
         DZ(LAYER) = LYRTOT + DZMAX

C**      --Apportioning routine start.
         DO 110 J = 1, NPLANT
           IF (PLANTZ(J).GE.LYRTOP) THEN
C               -- Plant fully contained in layer, full apportioning.
                FACTOR = DZ(LAYER) / PLANTZ(J)
           ELSE IF (PLANTZ(J).GT.LYRTOP-DZ(LAYER)) THEN
C          -- Plant partially contained in layer, partial apportioning.
                FACTOR = (PLANTZ(J)-(LYRTOP-DZ(LAYER))) / PLANTZ(J)
           ELSE
C          -- Plant is not contained in layer, set values to zero.
                FACTOR = 0.0
           END IF

           CANLAI(J,LAYER) = TOTLAI(J) * FACTOR
           DRYCAN(J,LAYER) = PLANTW(J) * FACTOR
           IF (ITYPE(J).NE.0) THEN
C          -- If plant is not dead.
              IF (FACTOR*TOTLAI(J).NE.0) THEN
C                -- Plant is within layer.
                 RLEAF(J,LAYER) = RLEAF0(J)*(CANLAI(J,LAYER)/TOTLAI(J))
               ELSE
C                -- Plant is outside of layer.
                 RLEAF(J,LAYER) = 0.0
              END IF
            ELSE
C             -- Plant is dead.
              RLEAF (J,LAYER) = 0.0
           END IF

C**        -- End apportioning routine.
  110    CONTINUE

C          --Update variables if more layers were added.
           IF (OLDNC.LT.LAYER) THEN
                TCDT(LAYER) = TCDT(OLDNC)
                VAPCDT(LAYER) = VAPCDT(OLDNC)
                WCANDT(LAYER) = WCANDT(OLDNC)
           END IF

C          --Calc midpoint of the layer
           IF (LAYER.EQ.1) THEN
               ZC(1) = 0.0
           ELSE
               ZC(LAYER) = PLANTZ(HTORDR(NPLANT))-(LYRTOP-DZ(LAYER)/2.)
           END IF

C          Check if done.  If so, drop out of loop
           IF (LAYER.EQ.1) THEN
                IF (ABS(PLANTZ(HTORDR(NPLANT))-LYRTOP).GT.0.0001) THEN
                     WRITE (6,*)'ERROR MATCHING TOP OF CANOPY'
                     WRITE (6,*)'  LYRTOP =  ', LYRTOP
                END IF
                GO TO 190
           END IF

C**        -- If top of layer equaled top of the next tallest plant
           IF (LYRTOP .GE. PLANTZ(HTORDR(HTINDX))-.0001) THEN
                LAISUM = LAISUM -
     >             TOTLAI(HTORDR(HTINDX))/PLANTZ(HTORDR(HTINDX))
                HTINDX = HTINDX + 1
           END IF

           LYRTOT = 0.0
           LAILYR = 0.0
           LAYER = LAYER - 1

      ELSE
C**   -- Top of layer extends beyond top of plant being considered.
C          -- Update calcs based on limitation of top of next plant.
           LYRTOT = LYRTOT + PLANTZ(HTORDR(HTINDX))-LYRTOP
           LAILYR = LAILYR + (A-LAILYR) *
     >              ((PLANTZ(HTORDR(HTINDX))-LYRTOP)/DZMAX)
           LAISUM = LAISUM -
     >             TOTLAI(HTORDR(HTINDX))/PLANTZ(HTORDR(HTINDX))
           LYRTOP = PLANTZ(HTORDR(HTINDX))
           HTINDX = HTINDX + 1
      END IF
      GO TO 100
C     -- End of loop

  190 CONTINUE
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE RTDIST (NPLANT,ITYPE,NS,ZS,ROOTDP,RROOT0)
C
C     THIS PROGRAM CALCULATES THE ROOT RESISTANCES AND FRACTION OF ROOTS
C     WITHIN EACH SOIL LAYER GIVEN THE MAXIMUM ROOTING DEPTH.  A
C     TRIANGULAR ROOTING DENSITY IS ASSUMED HAVING A MAXIMUM DENSITY AT
C     A DEPTH OF ZMXDEN AND ZERO DENSITY AT THE SURFACE AND THE MAXIMUM
C     ROOTING DEPTH.  ZMXDEN IS ASSUMED TO BE A CONSTANT 
C     FRACTION (RMXDEN) OF THE MAXIMUM ROOTING DEPTH.
C     (CURRENTLY RMXDEN IS SET IS DATA STATEMENT)
C***********************************************************************
      COMMON /CLAYRS/ DRYCAN(8,10),CANLAI(8,10),TOTLAI(8),IEVAP(8),
     >                RLEAF(8,10),RROOT(8,50),ROOTDN(8,50),TOTROT(8)
      REAL ZS(50),ROOTDP(8),RROOT0(8)
      INTEGER ITYPE(8)
C
      DATA RMXDEN/ 0.1/
C
      DO 30 J=1,NPLANT
        TOTROT(J)=0.0
        IF (ITYPE(J).NE.0 .AND. TOTLAI(J).NE.0.0) THEN
C         TRANSPIRING PLANT -- CALCULATE FRACTION OF ROOTS IN EACH
C         SOIL LAYER; START BY COMPUTING DEPTH OF MAXIMUM ROOT DENSITY
          ZMXDEN=RMXDEN*ROOTDP(J)
C
          ZMID1=0.0
          DO 10 I=1,NS
C           CALCULATE MID-POINT BETWEEN THIS AND NEXT NODE, I.E. THE 
C           LOWER BOUNDARY OF THIS LAYER
            IF (I.LT.NS) THEN
               ZMID2=(ZS(I+1)+ZS(I))/2.
              ELSE
               ZMID2=ZS(NS)
            END IF            
C
            IF (ZMID2.LT.ZMXDEN) THEN
C             BOTTOM OF LAYER IS LESS THAN DEPTH OF MAXIMUM DENSITY
              ROOTDN(J,I)=(ZMID2-ZMID1)*(ZMID2+ZMID1)/ZMXDEN
     >                      /ROOTDP(J)
             ELSE
C             BOTTOM OF LAYER IS BEYOND DEPTH OF MAXIMUM DENSITY
              IF (ZMID2.LT.ROOTDP(J)) THEN
C               BOTTOM OF LAYER IS WITHIN ROOTING DEPTH
                IF (ZMID1.LT.ZMXDEN) THEN
C                 LAYER STRATTLES DEPTH OF MAXIMUM DENSITY
                  AREA1=(ZMXDEN-ZMID1)*(ZMXDEN+ZMID1)/ZMXDEN
                  AREA2=(ZMID2-ZMXDEN)
     >                  *(2.-(ZMID2-ZMXDEN)/(ROOTDP(J)-ZMXDEN))
                  ROOTDN(J,I)=(AREA1+AREA2)/ROOTDP(J)
                 ELSE
C                 LAYER IS BEYOND DEPTH OF MAXIMUM ROOTING DENSITY BUT 
C                 IS FULLY WITHIN THE ROOTING DEPTH
                  ROOTDN(J,I)=(ZMID2-ZMID1)*(2.*ROOTDP(J)-ZMID2-ZMID1)         
     >                  /(ROOTDP(J)-ZMXDEN)/ROOTDP(J)
                END IF
               ELSE

C               BOTTOM OF LAYER IS BEYOND ROOTING DEPTH
                IF (ZMID1.LT.ROOTDP(J)) THEN
C                 TOP OF LAYER IS STILL WITHIN ROOTING DEPTH; THE 
C                 REMAINING FRACTION OF ROOTS ARE WITHIN THIS LAYER
                  ROOTDN(J,I)=1.0-TOTROT(J)
                 ELSE
C                 LAYER IS BEYOND THE ROOTING DEPTH -- NO ROOTS
                  ROOTDN(J,I)=0.0
                END IF
              END IF
            END IF
C           SUM THE TOTAL FRACTION OF ROOTS
            TOTROT(J)=TOTROT(J) + ROOTDN(J,I)
            ZMID1=ZMID2
   10     CONTINUE
C         CALCULATE EFFECTIVE ROOT CONDUCTANCE FOR EACH SOIL LAYER
          DO 20 I=1,NS
             RROOT(J,I)=RROOT0(J)*ROOTDN(J,I)/TOTROT(J)
   20     CONTINUE
        END IF
   30 CONTINUE
C
      RETURN
      END

C***********************************************************************
C
      SUBROUTINE SYSTM (NC,NPLANT,NSP,ZC,ZSP,ZMSRF,ZHSRF,ZERSRF,ZMSP,
     >                  ZHSP,HEIGHT,SOLR,TA,TCCRIT)
C
C     THIS SUBROUTINE DEFINES THE ROUGHNESS PARAMETERS FOR THE SYSTEM, 
C     AND DETERMINES WHETHER THE CANOPY WILL TRANSPIRE
C
C***********************************************************************
      COMMON /WINDV/ ZH,ZM,ZERO,USTAR,STABLE,WINDC(11),WINDR(10)
      COMMON /CLAYRS/ DRYCAN(8,10),CANLAI(8,10),TOTLAI(8),IEVAP(8),
     >                RLEAF(8,10),RROOT(8,50),ROOTDN(8,50),TOTROT(8)
C
      REAL ZC(11),ZSP(100),TCCRIT(8)
C
C     ZERO         plane of zero displacement in wind profile calculations (m)
      ZM=ZMSRF
      ZH=ZHSRF
      ZERO=ZERSRF
C
      IF (NC .GT. 0) THEN
         ZERO = 0.77*ZC(NC+1)
         ZM = 0.13*ZC(NC+1)
         ZH = 0.2*ZM
C        NO TRANPIRATION AT NIGHT
         IF (SOLR .LT. 10.0) THEN
            DO 10 J=1,NPLANT
               IEVAP(J)=0
   10       CONTINUE
         ELSE
            DO 20 J=1,NPLANT
C     IEVAP (J)    code for if plant j will transpire this time step (0=no; 1=yes)
               IEVAP(J)=1
C              CHECK IF TEMPERATURE IS TOO COLD FOR TRANSPIRATION
               IF (TA.LE.TCCRIT(J)) IEVAP(J)=0
   20       CONTINUE
         END IF
      ELSE
C
         IF (NSP .GT. 0) THEN
C           SNOW IS SURFACE MATERIAL -- DEFINE ROUGHNESS ELEMENTS AND 
C           ZERO PLANE OF DISPLACEMENT
            ZM=ZMSP
            ZH=ZHSP
            ZERO=ZSP(NSP+1)
C           DO NOT LET DISPLACEMENT PLANE BE ABOVE HEIGHT OF INSTRUMENTS
            IF (ZERO+ZM .GT. HEIGHT) ZERO=HEIGHT/2.
         END IF
      END IF
C
      RETURN
      END
C***********************************************************************        
C                                                                               
      SUBROUTINE SWRBAL (NPLANT,NC,NSP,NR,LANGLE,SWCAN,SWSNOW,SWRES,
     >  SWSOIL,SUNHOR,CANALB,ZSP,DZSP,RHOSP,ZR,RHOR,ALBRES,DIRRES,
     >  ALBDRY,ALBEXP,VLC,ALATUD,SLOPE,ASPECT,HRNOON,HAFDAY,DECLIN,
     >  HOUR,NHRPDT)
C                                                                               
C     THIS SUBROUTINE DETERMINES THE SHORT-WAVE RADIATION BALANCE FOR
C     FOR EACH NODE ABOVE THE GROUND SURFACE.
C
C     HAFDAY，DECLIN were calculated in SUBROUTINE DAY2HR
C***********************************************************************        
      COMMON /RADCAN/ TDIRCC(10),TDIFFC(10),DIRKL(9,10),DIFKL(9,10)
      COMMON /RADRES/ TDIREC(10),TDIFFU(10)
C
C     TDIRCC (I)   transmittance coefficient for direct radiation through canopy
C     TDIFFC (I)   transmittance coefficient for diffuse radiation through canopy
C     DIRKL (J,I)  product of direct radiation attenuation and leaf area index for
C                  plant j of canopy i. (NPLANT+1,I) is sum of all plants in layer i.
C     DIFKL (J,I)  product of diffuse radiation attenuation and leaf area index for
C                  plant j of canopy i. (NPLANT+1,I) is sum of all plants in layer i.

      REAL SWCAN(9,10),SWSNOW(100),SWRES(10),CANALB(8),
     >  ZSP(100),DZSP(100),RHOSP(100),ZR(10),RHOR(10)
      INTEGER HOUR,LANGLE(8)                                                   
C                                                                               
C**** SEPARATE THE TOTAL SOLAR RADIATION INTO DIRECT AND DIFFUSE                
      CALL SOLAR (DIRECT,DIFFUS,SUNSLP,ALTITU,SUNHOR,ALATUD,SLOPE,
     >                  ASPECT,HRNOON,HAFDAY,DECLIN,HOUR,NHRPDT)
      IF ((DIRECT+DIFFUS) .LE. 0.0) GO TO 60                                    
C
C**** DETERMINE THE SOIL AND SNOW ALBEDO
      ALBSOI=ALBDRY*EXP(-ALBEXP*VLC)
      IF (NSP.GT.0) THEN
         IF (NR .GT. 0) THEN
            ALBNXT=ALBRES
           ELSE
            ALBNXT=ALBSOI
         END IF
         CALL SNOALB (NSP,ALBSNO,ZSP,RHOSP,SUNSLP,ALBNXT)
      END IF
C
C
C**** DETERMINE THE SOLAR RADIATION BALANCE FOR EACH MATERIAL, STARTING         
C     FROM THE SURFACE MATERIAL AND WORKING DOWNWARD.                           
C                                                                               
C**** SOLAR RADIATION BALANCE OF THE CANOPY                                     
      IF (NC .GT. 0) THEN
         IF (NSP.EQ.0 .AND. NR.GT.0) THEN
C           RESIDUE LIES BENEATH CANOPY - DETERMINE SOLAR RADIATION 
C           BALANCE THROUGH CANOPY AND RESIDUE TOGETHER
            CALL SWRCAN (NPLANT,NC,NSP,NR,LANGLE,SWCAN,SWRES,DIRECT,
     >        DIFFUS,SUNSLP,ALTITU,CANALB,ZR,RHOR,ALBRES,DIRRES,ALBSOI)
         ELSE
C           CALCULATE RADIATION THROUGH CANOPY DOWN TO SNOW OR SOIL
            IF (NSP .GT. 0) THEN
               ALBNXT=ALBSNO
              ELSE
               ALBNXT=ALBSOI
            END IF
            CALL SWRCAN (NPLANT,NC,NSP,0,LANGLE,SWCAN,SWRES,DIRECT,
     >        DIFFUS,SUNSLP,ALTITU,CANALB,ZR,RHOR,ALBRES,DIRRES,ALBNXT)
C
         END IF
      END IF
C                                                                               
C**** SOLAR RADIATION BALANCE OF THE SNOWPACK
      IF (NSP.GT.0) CALL SWRSNO (NSP,SWSNOW,DZSP,RHOSP,DIRECT,DIFFUS,
     >                           SUNSLP,ALBSNO)
C
C
C**** SOLAR RADIATION BALANCE OF THE RESIDUE IF NO CANOPY PRESENT OR IF
C     RESIDUE IS COVERED BY SNOW 
      IF (NR.GT.0) THEN
         IF (NSP.GT.0 .OR. NC.EQ.0) THEN
            CALL SWRCAN (NPLANT,0,NSP,NR,LANGLE,SWCAN,SWRES,DIRECT,
     >        DIFFUS,SUNSLP,ALTITU,CANALB,ZR,RHOR,ALBRES,DIRRES,ALBSOI)
         END IF
      END IF
C                                                                               
C**** SOLAR RADIATION BALANCE FOR SOIL SURFACE                                  
      SWSOIL=(DIRECT+DIFFUS)*(1.-ALBSOI)                                        
C                                                                               
      RETURN                                                                    
C                                                                               
C**** NO SOLAR RADIATION -- SET ABSORBED SOLAR RADIATION TO ZERO                
C                                                                               
C     CANOPY NODES                                                              
   60 DO 70 I=1,NC                                                              
         DO 65 J=1,NPLANT+1
            SWCAN(J,I)=0.0
   65    CONTINUE
   70 CONTINUE                                                                  
C     DEFINE THE DIFFUSE RADIATION TRANSMISSION FOR LONG-WAVE BALANCE           
      IF (NC .GT. 0) CALL TRANSC (NPLANT,NC,LANGLE,TDIRCC,TDIFFC,DIFKL,
     >                            DIRKL,0.0,0.0)
C                                                                               
C     SNOWPACK NODES                                                            
      DO 80 I=1,NSP                                                             
         SWSNOW(I)=0.0                                                          
   80 CONTINUE                                                                  
C                                                                               
C     RESIDUE NODES                                                             
      DO 90 I=1,NR                                                              
         SWRES(I)=0.0                                                           
   90 CONTINUE                                                                  
C     DEFINE THE DIFFUSE RADIATION TRANSMISSION FOR LONG-WAVE BALANCE           
      IF (NR .GT. 0) CALL TRANSR (NR,TDIREC,TDIFFU,ZR,RHOR,0.785,DIRRES)                       
C                                                                               
C     SOIL SURFACE                                                              
      SWSOIL=0.0                                                                
C                                                                               
      RETURN                                                                    
      END                                                                       
C***********************************************************************
C
      SUBROUTINE SOLAR (DIRECT,DIFFUS,SUNSLP,ALTITU,SUNHOR,ALATUD,SLOPE,
     >                  ASPECT,HRNOON,HAFDAY,DECLIN,HOUR,NHRPDT)
C
C     THIS SUBROUTINE SEPARATES THE TOTAL RADIATION MEASURED ON THE
C     HORIZONTAL INTO THE DIRECT AND DIFFUSE ON THE LOCAL SLOPE.
C

C     HAFDAY，DECLIN were calculated in SUBROUTINE DAY2HR
C***********************************************************************
      COMMON /TIMEWT/ WT,WDT,DT
      COMMON /SWRCOE/ SOLCON,DIFATM,DIFRES,SNOCOF,SNOEXP
      INTEGER HOUR
C
C     SOLCON     solar constant (solar radiation at outer surface of atmosphere, W/m^2)
C     DIFATM       coefficient for atmospheric transmissivity of diffuse radiation
C     DIFRES       coefficient for diffuse radiation transmissivity through residue
C     SNOCOF       minimum snow depth for complete snowcover used in albedo calculations (m)
C     SNOEXP       exponent for the albedo of the snowpack
C     DECLIN       solar declination angle (radians)
C**** CHECK IF SUN HAS RISEN YET (OR IF IT HAS ALREADY SET)
      IF (SUNHOR .LE. 0.0) THEN
         DIRECT=0.0
         DIFFUS=0.0
         RETURN
      END IF
      SUNRIS=HRNOON - HAFDAY/0.261799
      SUNSET=HRNOON + HAFDAY/0.261799
C

C**** CALCULATE HOUR ANGLE AT WHICH THE SUN WILL BE DUE EAST/WEST IN
C     ORDER TO ADJUST AZIMUTH ANGLE FOR SOUTHERN AZIMUTHS
C     -- SIN(AZIMUTH) TELLS YOU ONLY THE EAST/WEST DIRECTION - NOT
C     WHETHER THE SUN IS NORTH/SOUTH.
      IF (ABS(DECLIN).GE.ABS(ALATUD)) THEN
C        LATITUDE IS WITHIN THE TROPICS (EQUATION WON'T WORK)
         HRWEST=3.14159
        ELSE
         HRWEST=ACOS(TAN(DECLIN)/TAN(ALATUD))
      END IF
C
C**** SUM UP VALUES AND FIND AVERAGE SUN POSITION FOR TIME STEP
      SINAZM=0.0
      COSAZM=0.0
      SUMALT=0.0
      COSALT=0.0
      SUNMAX=0.0
      DO 10 IHR=HOUR-NHRPDT,HOUR
         THOUR=IHR
C
C****    DETERMINE THE GEOMETRY OF THE SUN'S RAYS AT CURRENT TIME
         HRANGL=0.261799*(IHR-HRNOON)
         IF (THOUR .GT. SUNRIS  .AND.  THOUR .LT. SUNSET) THEN
C           SUN IS ABOVE HORIZON -- CALCULATE ITS ALTITUDE ABOVE THE
C           HORIZON (ALTITU) AND ANGLE FROM DUE NORTH (AZMUTH)
            SINALT=SIN(ALATUD)*SIN(DECLIN)
     >             + COS(ALATUD)*COS(DECLIN)*COS(HRANGL)
            ALTITU=ASIN(SINALT)
            AZM = ASIN(-COS(DECLIN)*SIN(HRANGL)/COS(ALTITU))
C           CORRECT AZIMUTH FOR SOUTHERN ANGLES
            IF (ALATUD-DECLIN .GT. 0.0) THEN
C              NORTHERN LATITUDES   (HRANGL=0.0 AT NOON)
               IF (ABS(HRANGL).LT.HRWEST) AZM=3.14159-AZM
              ELSE
C              SOUTHERN LATITUDES
               IF (ABS(HRANGL).GE.HRWEST) AZM=3.14159-AZM
            END IF
C           SUM CONDITIONS TO GET AVERAGE ALTITUDE AND AZMUTH
C           (OBTAIN AVERAGE BY SUMMING VECTOR COMPONENTS)
            SUN=SOLCON*SINALT
            SUMALT=SUMALT+SUN*SINALT
            COSALT=COSALT+SUN*COS(ALTITU)
            SINAZM=SINAZM+SUN*SIN(AZM)
            COSAZM=COSAZM+SUN*COS(AZM)
            SUNMAX=SUNMAX+SUN
         END IF
C
   10 CONTINUE
C
C**** DETERMINE AVERAGE SOLAR RADIATION, AVERAGE ALTITUDE AND AZIMUTH OF
C     THE SUN AND ANGLE ON LOCAL SLOPE
      IF (SUNMAX .EQ. 0) THEN
         ALTITU=0.0
         SUNSLP=0.0
        ELSE
         ALTITU=ATAN(SUMALT/COSALT)
         AZMUTH=ATAN2(SINAZM,COSAZM)
         SUNMAX=SUNMAX/(NHRPDT+1)
         SUNSLP=ASIN( SIN(ALTITU)*COS(SLOPE)
     >           + COS(ALTITU)*SIN(SLOPE)*COS(AZMUTH-ASPECT))
      END IF
C
C**** SEPARATE THE SOLAR RADIATION INTO DIRECT AND DIFFUSE COMPONENTS
      IF (ALTITU .LE. 0.0) THEN
C     SUN IS BELOW THE HORIZON - ALL RADIATION MUST BE DIFFUSE
         DIFFUS=SUNHOR
         DIRECT=0.0
         RETURN
      END IF
      TTOTAL=SUNHOR/SUNMAX
C     LIMIT TOTAL TRANSMISSIVITY TO MAXIMUM (DIFATM) WHICH WILL
C     CAUSE TDIFFU TO BE 0.0
      IF (TTOTAL .GT. DIFATM) TTOTAL = DIFATM
      TDIFFU=TTOTAL*(1. - EXP(0.6*(1.-DIFATM/TTOTAL)/(DIFATM-0.4)))
      DIFFUS=TDIFFU*SUNMAX
      DIRHOR=SUNHOR-DIFFUS
C
C**** NOW CALCULATE THE DIRECT SOLAR RADIATION ON THE LOCAL SLOPE
      IF (SUNSLP .LE. 0.0) THEN
C        SUN HAS NOT RISEN ON THE LOCAL SLOPE -- NO DIRECT RADIATION
         DIRECT=0.0
       ELSE
         DIRECT=DIRHOR*SIN(SUNSLP)/SIN(ALTITU)
C        IF THE SUN'S ALTITUDE IS NEAR ZERO, THE CALCULATED DIRECT
C        RADIATION ON THE SLOPING SURFACE MAY BE UNREALISTICALLY LARGE --
C        LIMIT DIRECT TO 5*DIRHOR (THIS IS APPROXIMATELY THE CASE WHEN
C        THE SUN IS 10 DEGREES ABOVE THE HORIZON AND THE SLOPING SURFACE
C        IS PERPENDICULAR TO THE SUN`S RAYS
         IF (DIRECT .GT. 5.*DIRHOR) DIRECT=5.*DIRHOR
      END IF
C
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE SNOALB (NSP,ALBSNO,ZSP,RHOSP,SUNSLP,ALBNXT)
C
C     THIS SUBROUTINE COMPUTES THE ALBEDO, THE EXTINCTION COEFFICIENT
C     AND THE AMOUNT OF SOLAR RADIATION ABSORBED BY EACH LAYER
C
C     OUT:ALBSNO
C***********************************************************************
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /SWRCOE/ SOLCON,DIFATM,DIFRES,SNOCOF,SNOEXP
      COMMON /SPPROP/ G1,G2,G3,EXTSP,TKSPA,TKSPB,TKSPEX,VDIFSP,VAPSPX
      REAL LF,LS,LV
C
      REAL RHOSP(100),ZSP(100)
C
C     DETERMINE THE ALBEDO OF THE SNOW
      SPGRAV = RHOSP(1)/RHOL
      GRAIN = G1 + G2*SPGRAV*SPGRAV + G3*SPGRAV**4
      ALBSNO = 1. - 0.206*EXTSP*SQRT(GRAIN)
      IF (ALBSNO .LT. 0.35) ALBSNO = 0.35
C
      IF (ZSP(NSP+1) .LE. 0.04) THEN
C        SNOWPACK IS LESS THAN 4.0 CM -- ALBEDO IS AFFECTED BY
C        UNDERLYING MATERIAL
         ABS = 1. - ALBNXT
         W = 2.*(1. - ALBSNO)/(1. + ALBSNO)
C        EXC = EXTINCTION COEFFICIENT --> CONVERT FROM 1/CM TO 1/M
         EXC = EXTSP*SPGRAV*SQRT(1.0/GRAIN)*100.
         Y = EXP(-EXC*ZSP(NSP+1))*(ABS + W*(ABS/2.-1.))
     >   / (W*(ABS/2.-1.)*COSH(EXC*ZSP(NSP+1))-ABS*SINH(EXC*ZSP(NSP+1)))
         ALBSNO = (1. - W*(1.-Y)/2)/(1. + W*(1.-Y)/2)
      END IF
C
      IF (ZSP(NSP+1) .LE. SNOCOF) THEN
C        SNOCOF IS MINIMUM DEPTH OF SNOW FOR COMPLETE GROUND COVER --
C        SNOW IS NOT COMPLETELY COVERED WITH SNOW
         ALBSNO = ALBNXT + (ALBSNO-ALBNXT)*(ZSP(NSP+1)/SNOCOF)**SNOEXP
      END IF
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE SWRSNO (NSP,SWSNOW,DZSP,RHOSP,DIRECT,DIFFUS,
     >                   SUNSLP,ALBSNO)
C
C     THIS SUBROUTINE COMPUTES THE AMOUNT OF OF SOLAR RADIATION
C     ABSORBED BY EACH LAYER WITHIN THE SNOWPACK
C
C***********************************************************************
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /SPPROP/ G1,G2,G3,EXTSP,TKSPA,TKSPB,TKSPEX,VDIFSP,VAPSPX
      REAL LF,LS,LV
C
      REAL SWSNOW(100),RHOSP(100),DZSP(100)
C
C     COMPUTE THE TOTAL ABSORBED FOR THE ENTIRE SNOWPACK
      TOTAL = (DIRECT + DIFFUS)*(1. - ALBSNO)
C
C     CALCULATE EXTINCTION COEFF. AND RADIATION ABSORBED BY EACH LAYER
      DO 10 I= 1,NSP
         SPGRAV = RHOSP(I)/RHOL
         GRAIN = G1 + G2*SPGRAV*SPGRAV + G3*SPGRAV**4
C        EXC = EXTINCTION COEFFICIENT --> CONVERT FROM 1/CM TO 1/M
         EXC = EXTSP*SPGRAV*SQRT(1.0/GRAIN)*100.
         PERCAB = 1.0 - EXP(-EXC*DZSP(I))
         SWSNOW(I) = TOTAL*PERCAB
         TOTAL = TOTAL - SWSNOW(I)
   10 CONTINUE
C
C     DEFINE THE DIFFUSE RADIATION OUT THE BOTTOM OF THE SNOWPACK
C     (THIS WILL BE ABSORBED BY EITHER THE RESIDUE OR THE SOIL.)
      DIFFUS = TOTAL
      DIRECT = 0.0
      RETURN
      END
C***********************************************************************        
C                                                                               
      SUBROUTINE SWRCAN (NPLANT,NC,NSP,NR,LANGLE,SWCAN,SWRES,DIRECT,
     >  DIFFUS,SUNSLP,ALTITU,CANALB,ZR,RHOR,ALBRES,DIRRES,ALBNXT)
C                                                                               
C     THIS SUBROUTINE CALCULATES THE AMOUNT OF SOLAR RADIATION ABSORBED         
C     BY EACH CANOPY AND RESIDUE NODE
C                                                                               
C***********************************************************************        
      COMMON /RADCAN/ TDIRCC(10),TDIFFC(10),DIRKL(9,10),DIFKL(9,10)
      COMMON /RADRES/ TDIREC(10),TDIFFU(10)
C
      REAL SWCAN(9,10),SWRES(10),CANALB(8),DIRALB(10),DIFALB(10),
     >     ZR(10),RHOR(10),DIR(21),DOWN(21),UP(21)
      INTEGER LANGLE(8)
C
C**** OBTAIN TRANSMISSIVITY TO DIRECT AND DIFFUSE RADIATION THROUGH 
C     CANOPY AND/OR RESIDUE DEPENDING ON LAYERING OF MATERIALS
C
C     SUNSLP       angle that the sun makes with the local slope (radians)

      IF (NC .GT. 0) THEN
         CALL TRANSC (NPLANT,NC,LANGLE,TDIRCC,TDIFFC,DIFKL,DIRKL,
     >                SUNSLP,ALTITU)
C        CALCULATE WEIGHTED ALBEDO FOR DIRECT AND DIFFUSE RADIATION 
C        FOR EACH OF THE CANOPY LAYERS
         DO 15 I=1,NC
            TOTDIR=0.0
            TOTDIF=0.0
            DO 10 J=1,NPLANT
               TOTDIR = TOTDIR + CANALB(J)*DIRKL(J,I)
               TOTDIF = TOTDIF + CANALB(J)*DIFKL(J,I)
   10       CONTINUE
            DIRALB(I)=TOTDIR/DIRKL(NPLANT+1,I)
            DIFALB(I)=TOTDIF/DIFKL(NPLANT+1,I)
   15    CONTINUE
      END IF
C
      IF (NR .GT. 0) 
     >   CALL TRANSR (NR,TDIREC,TDIFFU,ZR,RHOR,SUNSLP,DIRRES)
      NR1=NR
      IF (NR.EQ.1) THEN
C        RADIATION EXCHANGE IS NOT ACCURATE WITH ONLY ON RESIDUE LAYER.
C        SPLIT RESIDUE LAYER IN TWO AND REDEFINE TRANSMISSION COEFF.
C        FOR TWO RESIDUE LAYERS.
         NR1=2
         TDIREC(1)=SQRT(TDIREC(1))
         TDIFFU(1)=SQRT(TDIFFU(1))
         TDIREC(2)=TDIREC(1)
         TDIFFU(2)=TDIFFU(1)
      END IF
C                                                                               
C**** RADIATION EXCHANGE WITHIN THE CANOPY AND RESIDUE
C                                                                               
C     CALCULATE THE DIRECT RADIATION PASSING THROUGH EACH LAYER, AND            
C     INITIALIZE THE UPWELLING DIFFUSE COMPONENT                                
      DIR(1)=DIRECT                                                             
      DOWN(1)=DIFFUS                                                         
      UP(1)=0.0                                                                 
C     RADIATION THROUGH CANOPY LAYERS      
      DO 25 I=2,NC+1
         DIR(I)=DIR(I-1)*TDIRCC(I-1)
         UP(I)=0.0
   25 CONTINUE
C     RADIATION THROUGH RESIDUE LAYERS      
      DO 30 I=2,NR1+1
         K=I+NC
         DIR(K)=DIR(K-1)*TDIREC(I-1)                                            
         UP(K)=0.0                                                              
   30 CONTINUE                                                                  
C                                                                               
C**** START ITERATIVE LOOP TO BALANCE THE RADIATION AT EACH LAYER - 
C     SUFFICIENT CONVERGENCE SHOULD BE MET AFTER GOING THROUGH LOOP             
C     TWICE TIMES THE NUMBER OF LAYERS                                          
C                                                                               
      DO 60 ITER=1,2*(NC+NR1)
C                                                                               
C****    CALCULATE THE DOWNWELLING DIFFUSE RADIATION 
C        CANOPY LAYERS 
         DO 35 I=2,NC+1                                                         
            DOWN(I)=DOWN(I-1)*TDIFFC(I-1) 
     >             +UP(I)*(1.-TDIFFC(I-1))*DIFALB(I-1)
   35    CONTINUE                                                               
C        RESIDUE LAYERS 
         DO 40 I=2,NR1+1                                                         
            K=I+NC
            DOWN(K)=DOWN(K-1)*TDIFFU(I-1) +UP(K)*(1.-TDIFFU(I-1))*ALBRES
   40    CONTINUE                                                               
C                                                                               
C****    CALCULATE THE UPWELLING DIFFUSE RADIATION STARTING WITH THAT 
C        REFLECTED FROM THE SNOW OR SOIL BENEATH
         UP(NC+NR1+1)=(DOWN(NC+NR1+1) + DIR(NC+NR1+1))*ALBNXT
         DO 45 I=NR1,1,-1                                                        
            K=I+NC
            UP(K)= UP(K+1)*TDIFFU(I) + DOWN(K)*(1-TDIFFU(I))*ALBRES             
     >           + DIR(K)*(1.-TDIREC(I))*ALBRES                                 
   45    CONTINUE                                                               
         DO 50 I=NC,1,-1                                                        
            UP(I)= UP(I+1)*TDIFFC(I) + DOWN(I)*(1-TDIFFC(I))*DIFALB(I)
     >           + DIR(I)*(1.-TDIRCC(I))*DIRALB(I)
   50    CONTINUE                                                               
C                                                                               
   60 CONTINUE                                                                  
C                                                                               
C**** DETERMINE THE AMOUNT OF RADIATION ABSORBED AT EACH CANOPY NODE BY 
C     WEIGHTING ACCORDING TO ALBEDO, LEAF AREA, AND ATTENUATION FACTOR
      DO 65 I=1,NC
         TOTDIR=0.0
         TOTDIF=0.0
         DO 63 J=1,NPLANT
            TOTDIR = TOTDIR + (1.-CANALB(J))*DIRKL(J,I)
            TOTDIF = TOTDIF + (1.-CANALB(J))*DIFKL(J,I)
   63    CONTINUE
         SWCAN(NPLANT+1,I)=0.0
         DO 64 J=1,NPLANT
            FRACT1=(1.-CANALB(J))*DIRKL(J,I)/TOTDIR
            FRACT2=(1.-CANALB(J))*DIFKL(J,I)/TOTDIF
            SWCAN(J,I)=FRACT1*(1-DIRALB(I))*(1-TDIRCC(I))*DIR(I)
     >            + FRACT2*(1-DIFALB(I))*(1-TDIFFC(I))*(DOWN(I)+UP(I+1))
            SWCAN(NPLANT+1,I)=SWCAN(NPLANT+1,I)+SWCAN(J,I)
   64    CONTINUE
   65 CONTINUE
C**** DETERMINE THE AMOUNT OF RADIATION ABSORBED AT EACH RESIDUE NODE
      ABSORB=0.0                                                                
      DO 70 I=1,NR1                                                              
         K=I+NC
         SWRES(I)=DIR(K)+DOWN(K)+UP(K+1)-DIR(K+1)-DOWN(K+1)-UP(K)               
         ABSORB=ABSORB + SWRES(I)                                               
   70 CONTINUE                                                                  
      IF (NR .EQ. 1) THEN
         SWRES(1)=SWRES(1)+SWRES(2)
         SWRES(2)=0.0
      END IF
C                                                                               
C**** IF SNOW LAYER LIES ABOVE THE RESIDUE, ALL RADIATION LEAVING THE           
C     THE RESIDUE IS ASSUMED TO BE REFLECTED BACK AND ABSORBED.                 
C     SPLIT UP THIS RADIATION PROPORTIONATELY BETWEEN THE LAYERS                
      IF (NR.GT.0 .AND. NSP.GT.0) THEN 
         IF (ABSORB .GT. 0.0) THEN
         DO 80 I=1,NR
            SWRES(I)=SWRES(I) + UP(1)*SWRES(I)/ABSORB
   80    CONTINUE
         END IF
      END IF
C                                                                               
C**** DEFINE THE AMOUNT OF RADIATION AVAILABLE AT THE SOIL SURFACE              
      DIRECT=DIR(NC+NR1+1)
      DIFFUS=DOWN(NC+NR1+1)
C                                                                               
      RETURN                                                                    
      END                                                                       
C***********************************************************************
C
      SUBROUTINE TRANSC (NPLANT,NC,LANGLE,TDIRCC,TDIFFC,DIFKL,DIRKL,
     >                   SUNSLP,ALTITU)
C
C     THIS SUBROUTINE CALCULATES THE DIRECT AND DIFFUSE RADIATION
C     TRANSMISSION COEFFICIENTS AND THE SUM OF THE ATTENUATION FACTOR
C     TIMES LEAF AREA INDEX FOR EACH CANOPY LAYER. (THIS SUM IS USED IN
C     CALCULATING AN EFFECTIVE ALBEDO FOR THE CANOPY LAYER AND IN
C     PROPORTIONING ABSORBED RADIATION INTO EACH OF THE PLANT SPECIES.
C
C     INOUT:TDIRCC,TDIFFC,DIFKL,DIRKL
C***********************************************************************
      COMMON /CLAYRS/ DRYCAN(8,10),CANLAI(8,10),TOTLAI(8),IEVAP(8),
     >                RLEAF(8,10),RROOT(8,50),ROOTDN(8,50),TOTROT(8)
C
      REAL TDIRCC(10),TDIFFC(10),DIRKL(9,10),DIFKL(9,10)
      INTEGER LANGLE(8)
C
      REAL DIFCAN(2)
      DATA DIFCAN(1)/0.6807 /DIFCAN(2)/0.7815/
C
C**** DETERMINE THE TRANSMISSIVITY TO DIRECT AND DIFFUSE RADIATION
C     BETWEEN CANOPY NODES
      DO 20 I=1,NC
         DIRKL(NPLANT+1,I)=0.0
         DIFKL(NPLANT+1,I)=0.0
         IF (SUNSLP .LE. 0.0) THEN
C           SUN IS ON OR BELOW HORIZON OF LOCAL SLOPE - SET DIRKL(J,I)
C           ALL EQUAL AND SET DIRECT TRANSMISSIVITY TO ZERO.
            DO 10 J=1,NPLANT
               DIRKL(J,I)=1.
               DIFKL(J,I)=CANLAI(J,I)*DIFCAN(LANGLE(J)+1)
               DIRKL(NPLANT+1,I)=DIRKL(NPLANT+1,I)+DIRKL(J,I)
               DIFKL(NPLANT+1,I)=DIFKL(NPLANT+1,I)+DIFKL(J,I)
   10       CONTINUE
            TDIRCC(I)=0.0
            TDIFFC(I)=EXP(-DIFKL(NPLANT+1,I))
C
           ELSE
C           CALCULATE TRANSMISSIVITY BASED ON SUN AND LEAF ANGLE
            DO 15 J=1,NPLANT
               IF (LANGLE(J) .EQ. 1) THEN
C                 LEAF ANGLE INCLINATION IS RANDOM
                  DIRKL(J,I)=CANLAI(J,I)/(2*SIN(SUNSLP))
                ELSE
C                 LEAF ANGLE INCLINATION IS VERTICLE
                  IF (ALTITU .GE. 1.57) THEN
                     DIRKL(J,I) = 0.0
                    ELSE
                     DIRKL(J,I)=CANLAI(J,I)/(1.571*TAN(ALTITU))
                  END IF
               END IF
               DIFKL(J,I)=CANLAI(J,I)*DIFCAN(LANGLE(J)+1)
               DIRKL(NPLANT+1,I)=DIRKL(NPLANT+1,I)+DIRKL(J,I)
               DIFKL(NPLANT+1,I)=DIFKL(NPLANT+1,I)+DIFKL(J,I)
   15       CONTINUE
C           TRANSMISSIVITY OF LAYER IS EQUAL TO EXPONENT OF THE SUM OF
C           ATTENUATION FACTOR TIMES LEAF AREA INDEX
            TDIRCC(I)=EXP(-DIRKL(NPLANT+1,I))
            TDIFFC(I)=EXP(-DIFKL(NPLANT+1,I))
         END IF
   20 CONTINUE
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE TRANSR (NR,TDIREC,TDIFFU,ZR,RHOR,ANGLE,DIRRES)
C
C     THIS SUBROUTINE CALCULATES THE DIRECT AND DIFFUSE TRANSMISSION
C     COEFFICIENTS FOR THE RESIDUE LAYERS, GIVEN THE ANGLE BETWEEN THE
C     SURFACE AND THE LINE OF DIRECT APPROACH (SUN-SLOPE ANGLE)
C
C     TRANSMISSION COEFFICIENTS ARE BASED ON THE TONNES/HA OF RESIDUE
C***********************************************************************
      COMMON /SWRCOE/ SOLCON,DIFATM,DIFRES,SNOCOF,SNOEXP
C
      REAL TDIREC(10),TDIFFU(10),ZR(10),RHOR(10)

C**** DETERMINE THE TRANSMISSIVITY TO DIRECT AND DIFFUSE RADIATION
C     BETWEEN RESIDUE NODES BASED ON FRACTION OF SURFACE COVERED
C     CALCULATED FROM THE WEIGHT OF RESIDUE (TON/HA)
      DO 10 I=1,NR
C        DETERMINE THE WEIGHT OF THE RESIDUE AT EACH NODE
         IF (I .EQ. 1) THEN
            W=RHOR(I)*(ZR(2)-ZR(1))*10.
            IF (NR .NE. 1) W=W/2.
           ELSE
            IF (I .EQ. NR) THEN
               W=RHOR(I)*(ZR(NR+1)-ZR(NR)+(ZR(NR)-ZR(NR-1))/2)*5.
              ELSE
               W=RHOR(I)*(ZR(I+1)-ZR(I-1))*5.
            END IF
         END IF
C        NOW DETERMINE THE DIRECT AND DIFFUSE TRANMISSIVITIES
         IF (ANGLE .LE. 0.0) THEN
C           SUN IS NOT ABOVE HORIZON OF LOCAL SLOPE -- SET DIRECT
C           TRANSMISSIVITY TO ZERO
            TDIREC(I)=0.0
           ELSE
            TDIREC(I)=EXP(-DIRRES*W)*SIN(ANGLE)
         END IF
         TDIFFU(I)=DIFRES*EXP(-DIRRES*W)
   10 CONTINUE
C
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE UPDATE (NS,NR,NSP,NC,NPLANT,NSALT,ICES,ICESDT,TS,TSDT,
     > MAT,MATDT,CONC,CONCDT,VLC,VLCDT,VIC,VICDT,SALT,SALTDT,
     > TR,TRDT,VAPR,VAPRDT,GMC,GMCDT,TC,TCDT,VAPC,VAPCDT,WCAN,WCANDT,
     > PCAN,PCANDT,TSP,TSPDT,DLW,DLWDT,ICESP,ICESPT)
C
C     THIS SUBROUTINE UPDATES THE BEGINNING OF TIME STEP VALUES FOR THE
C     NEW TIME STEP.
C
C***********************************************************************
      REAL TS(50),TSDT(50),MAT(50),MATDT(50),CONC(10,50),CONCDT(10,50),
     >     VLC(50),VLCDT(50),VIC(50),VICDT(50),SALT(10,50),SALTDT(10,50)
      REAL TR(10),TRDT(10),VAPR(10),VAPRDT(10),GMC(10),GMCDT(10)
      REAL TC(11),TCDT(11),VAPC(11),VAPCDT(11),WCAN(10),WCANDT(10),
     >     PCAN(8),PCANDT(8)
      REAL TSP(100),TSPDT(100),DLW(100),DLWDT(100)
      INTEGER ICES(50),ICESDT(50),ICESP(100),ICESPT(100)
C
C     SOIL PROPERTIES
      DO 15 I=1,NS
         TS(I)=TSDT(I)
         MAT(I)=MATDT(I)
         VLC(I)=VLCDT(I)
         VIC(I)=VICDT(I)
         ICES(I)=ICESDT(I)
         DO 10 J=1,NSALT
            SALT(J,I)=SALTDT(J,I)
            CONC(J,I)=CONCDT(J,I)
   10    CONTINUE
   15 CONTINUE
C
C     RESIDUE PROPERTIES
      DO 20 I=1,NR
         TR(I)=TRDT(I)
         VAPR(I)=VAPRDT(I)
         GMC(I)=GMCDT(I)
   20 CONTINUE
C
C     SNOWPACK PROPERTIES
      DO 30 I=1,NSP
         ICESP(I) = ICESPT(I)
         TSP(I) = TSPDT(I)
         DLW(I) = DLWDT(I)
   30 CONTINUE
C
C     CANOPY PROPERTIES
      DO 40 I=1,NC
         TC(I)=TCDT(I)
         VAPC(I)=VAPCDT(I)
         WCAN(I)=WCANDT(I)
   40 CONTINUE
      DO 45 J=1,NPLANT
         PCAN(J)=PCANDT(J)
   45 CONTINUE
C
      RETURN
      END
C***********************************************************************        
C                                                                               
      SUBROUTINE LWRBAL (NC,NSP,NR,NPLANT,TA,TADT,TC,TCDT,TSP,TSPDT,
     >         TR,TRDT,TS,TSDT,CLOUDS,LWCAN,LWSNOW,LWRES,LWSOIL)
C                                                                               
C     THIS SUBROUTINE CALLS SUBROUTINES TO SET UP THE LONG-WAVE 
C     RADIATION BALANCE FOR SYSTEM.
C                                                                               
C***********************************************************************        
      COMMON /TIMEWT/ WT,WDT,DT                                                 
      COMMON /LWRCOF/ STEFAN,EMATM1,EMATM2,EMITC,EMITR,EMITSP,EMITS             
C
      REAL TC(11),TCDT(11),TSP(100),TSPDT(100),TR(10),TRDT(10),                 
     >     TS(50),TSDT(50),LWCAN(9,10),LWSNOW(2),LWRES(10),LWSOIL
C                                                                               
C**** DETERMINE THE LONGWAVE RADIATION FLUX ABOVE EACH NODE.  START WITH        
C     ATMOSPHERE AND WORK DOWNWARDS.                                            
C                                                                               
      CALL LWRATM (ABOVE,CLOUDS,TA,TADT)                                        
C                                                                               
      IF (NC .GT. 0) THEN
C        RADIATION BALANCE FOR CANOPY NODES                                        
         DO 25 I=1,NC                                                              
           DO 20 J=1,NPLANT+1
              LWCAN(J,I)=0.0
   20      CONTINUE
   25    CONTINUE                                                                  
         CALL LWRCAN (LWCAN,NPLANT,ABOVE,1,NC,1,TC,TCDT)
      END IF
C                                                                               
      IF (NSP .GT. 0) THEN
C        RADIATION BALANCE FOR SNOW PACK                                           
         LWSNOW(1) = 0.0
         LWSNOW(2) = 0.0
         CALL LWRSNO (LWSNOW(1),LWSNOW(2),ABOVE,1,TSP,TSPDT)
C        NO LWR BENEATH THE SNOW                                                   
         LWSNOW(1) = LWSNOW(1) + LWSNOW(2)
         LWSNOW(2) = 0.0
         LWSOIL = 0.0                                                              
         DO 30 I=1,NR                                                              
            LWRES(I) = 0.0                                                         
   30    CONTINUE                                                                  
         BELOW=ABOVE
         GO TO 70
      END IF
C                                                                               
      IF (NR .GT. 0) THEN
C        RADIATION BALANCE FOR THE RESIDUE LAYERS                                  
         DO 40 I=1,NR                                                              
            LWRES(I)=0.0                                                           
   40    CONTINUE                                                                  
         CALL LWRRES (LWRES,ABOVE,1,NR,1,TR,TRDT)                                  
      END IF
C                                                                               
C     RADIATION BALANCE FOR THE SOIL SURFACE                                    
      TSK=TS(1) + 273.16                                                        
      BELOW=EMITS*STEFAN*((TSK**4) + 4.*WDT*(TSK**3)*(TSDT(1)-TS(1)))           
      LWSOIL=EMITS*ABOVE - BELOW                                                
C                                                                               
C**** DETERMINE THE LONGWAVE RADIATION FLUX BELOW EACH NODE.  START AT          
C     THE SOIL SURFACE AND WORK UPWARDS.                                        
C                                                                               
      IF (NR .GT. 0) THEN
C        RADIATION BALANCE FOR THE RESIDUE LAYERS                                  
         CALL LWRRES (LWRES,BELOW,NR,1,-1,TR,TRDT)                                 
      END IF
C
   70 IF (NC .GT. 0) THEN
C        RADIATION BALANCE FOR CANOPY NODES                                        
         CALL LWRCAN (LWCAN,NPLANT,BELOW,NC,1,-1,TC,TCDT)
      END IF
C                                                                               
      RETURN                                                                    
      END                                                                       
C***********************************************************************
C
      SUBROUTINE LWRCAN (LWCAN,NPLANT,IN,NSTART,NSTOP,NCOUNT,TC,TCDT)
C
C     THIS SUBROUTINE SUMS THE NET LONG-WAVE RADIATION EXCHANGE
C     FOR A GIVEN SIDE OF EACH CANOPY LAYER (ABOVE OR BELOW).  THE
C     VALUE OF "IN" COMING INTO THE SUBROUTINE IS LONG-WAVE FLUX
C     ENTERING A GIVEN SIDE OF THE CANOPY, AND THE VALUE OF "IN"
C     LEAVING THE SUBROUTINE IS THE LONG-WAVE FLUX LEAVING THE CANOPY
C     ON THE OPPOSITE SIDE.
C
C***********************************************************************
      COMMON /TIMEWT/ WT,WDT,DT
      COMMON /LWRCOF/ STEFAN,EMATM1,EMATM2,EMITC,EMITR,EMITSP,EMITS
      COMMON /RADCAN/ TDIRCC(10),TDIFFC(10),DIRKL(9,10),DIFKL(9,10)
      REAL IN,LWCAN(9,10),TC(11),TCDT(11)
C
      DO 10 I=NSTART,NSTOP,NCOUNT
         TRK=TC(I)+273.16
         OUT=(1-TDIFFC(I))*EMITC*STEFAN
     >        *(TRK**4 + 4.*WDT*(TRK**3)*(TCDT(I)-TC(I)))
         LWCAN(NPLANT+1,I)=LWCAN(NPLANT+1,I)+(1-TDIFFC(I))*EMITC*IN-OUT
         IN = TDIFFC(I)*IN + OUT
         DO 5 J=1,NPLANT
            LWCAN(J,I)=LWCAN(NPLANT+1,I)*DIFKL(J,I)/DIFKL(NPLANT+1,I)
    5    CONTINUE
   10 CONTINUE
C
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE LWRATM (ABOVE,CLOUDS,TA,TADT)
C
C     THIS SUBROUTINE DETERMINES THE LONG-WAVE RADIATIVE FLUX FROM THE
C     ATMOSPHERE.
C
C***********************************************************************
      COMMON /TIMEWT/ WT,WDT,DT
      COMMON /LWRCOF/ STEFAN,EMATM1,EMATM2,EMITC,EMITR,EMITSP,EMITS
C
C**** DETERMINE THE EMISSIVITY OF CLEAR-SKY ATMOSPHERE
      AVGTMP = WT*TA + WDT*TADT
      EMITAT=1. - EMATM1*EXP(-EMATM2*AVGTMP*AVGTMP)
C
C**** ADJUST CLEAR-SKY EMISSIVITY FOR FRACION OF CLOUD COVER
CC    EMITAT=EMITAT + CLOUDS*(1-EMITAT-8./(273.16+AVGTMP))
CC    USE EQUATION GIVEN CAMPBELL, 1985 (SOIL PHYSICS WITH BASIC)
      EMITAT=(1.-0.84*CLOUDS)*EMITAT + 0.84*CLOUDS
C
C**** NOW CALCULATE THE LONG-WAVE RADIATIVE FLUX
      ABOVE = EMITAT*STEFAN*(AVGTMP+273.16)**4
C
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE LWRSNO (LWRTOP,LWRBOT,IN,IBOTM,TSP,TSPDT)
C
C     THIS SUBROUTINE SUMS THE NET LONG-WAVE RADIATION EXCHANGE
C     FOR A GIVEN SIDE OF THE SNOWPACK (TOP OR BOTTOM).  THE
C     VALUE OF "IN" COMING INTO THE SUBROUTINE IS LONG-WAVE FLUX
C     ENTERING A GIVEN SIDE OF THE SNOWPACK, AND THE VALUE OF "IN"
C     LEAVING THE SUBROUTINE IS THE LONG-WAVE FLUX EMITTED BY THE
C     THE SIDE OF THE SNOWPACK DESIGNATED BY IBOTM.
C     (IBOTM = 1 --> TOP OF SNOW;  IBOTM = NSP --> BOTTOM OF SNOWPACK)
C
C***********************************************************************
      COMMON /TIMEWT/ WT,WDT,DT
      COMMON /LWRCOF/ STEFAN,EMATM1,EMATM2,EMITC,EMITR,EMITSP,EMITS
      REAL LWRTOP,LWRBOT,IN,TSP(100),TSPDT(100)
C
      TSPK = TSP(IBOTM)+273.16
      LWRTOP = LWRTOP + EMITSP*IN
      OUT = EMITSP*STEFAN
     >      *(TSPK**4 + 4.*WDT*(TSPK**3)*(TSPDT(IBOTM)-TSP(IBOTM)))
      LWRBOT = LWRBOT - OUT
      IN = OUT
C
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE LWRRES (LWRES,IN,NSTART,NSTOP,NCOUNT,TR,TRDT)
C
C     THIS SUBROUTINE SUMS THE NET LONG-WAVE RADIATION EXCHANGE
C     FOR A GIVEN SIDE OF EACH RESIDUE LAYER (ABOVE OR BELOW).  THE
C     VALUE OF "IN" COMING INTO THE SUBROUTINE IS LONG-WAVE FLUX
C     ENTERING A GIVEN SIDE OF THE RESIDUE, AND THE VALUE OF "IN"
C     LEAVING THE SUBROUTINE IS THE LONG-WAVE FLUX LEAVING THE RESIDUE
C     ON THE OPPOSITE SIDE.
C
C***********************************************************************
      COMMON /TIMEWT/ WT,WDT,DT
      COMMON /LWRCOF/ STEFAN,EMATM1,EMATM2,EMITC,EMITR,EMITSP,EMITS
      COMMON /RADRES/ TDIREC(10),TDIFFU(10)
C
      REAL IN,LWRES(10),TR(10),TRDT(10)
C
      DO 10 I=NSTART,NSTOP,NCOUNT
         TRK=TR(I)+273.16
         OUT=(1-TDIFFU(I))*EMITR*STEFAN
     >        *(TRK**4 + 4.*WDT*(TRK**3)*(TRDT(I)-TR(I)))
         LWRES(I)=LWRES(I) + (1-TDIFFU(I))*EMITR*IN - OUT
         IN = TDIFFU(I)*IN + OUT
   10 CONTINUE
C
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE LWRMAT (NC,NSP,NR,TC,TSP,TR,TS,ICESPT)
C
C     THIS SUBROUTINE CALCULATES THE CONTRIBUTION TO THE ENERGY BALANCE
C     MATRIX BY LONG-WAVE RADIATION.
C
C***********************************************************************
      COMMON /TIMEWT/ WT,WDT,DT
      COMMON /MATRIX/ A1(150),B1(150),C1(150),D1(150),
     >                A2(150),B2(150),C2(150),D2(150)
      COMMON /LWRCOF/ STEFAN,EMATM1,EMATM2,EMITC,EMITR,EMITSP,EMITS
      COMMON /RADCAN/ TDIRCC(10),TDIFFC(10),DIRKL(9,10),DIFKL(9,10)
      COMMON /RADRES/ TDIREC(10),TDIFFU(10)
      REAL TC(11),TSP(100),TR(10),TS(50)
      INTEGER ICESPT(100)
C
      N = 1
      MATERL=2
      B1(N) = 0.0
C
C     RADIATION BALANCE MATRIX COEFFICIENTS FOR CANOPY NODES
C
C
C     EMITC        emissivity of the canopy elements
      IF (NC .GT. 0) THEN
         EMITNC=(1.-TDIFFC(1))*EMITC
         TC3=STEFAN*4.*WDT*(TC(1)+273.16)**3
         B1(N) =-2.*EMITNC*TC3
         C1(N) = EMITNC
         A1(N+1) = EMITNC*TC3

         N = N + 1
C
         DO 10 I=2,NC
            EMITNC=(1.-TDIFFC(I))*EMITC
            TC3=STEFAN*4.*WDT*(TC(I)+273.16)**3
            C1(N-1)=C1(N-1)*EMITNC*TC3
            A1(N) = A1(N)*EMITNC
            B1(N) = -2.*EMITNC*TC3
            C1(N)= EMITNC
            A1(N+1) = EMITNC*TC3
            N = N + 1
   10    CONTINUE
         MATERL=MATERL + 1
      END IF
C
C     RADIATION BALANCE MATRIX COEFFICIENTS FOR SNOW PACK
      IF (NSP .GT. 0) THEN
      	IF (ICESPT(1) .GT. 0) THEN
C           TEMP IS KNOWN - ENERGY BALANCE BASED ON LIQUID CONTENT
            IF (MATERL .GT. 2) THEN
               A1(N) = A1(N)*EMITSP
               C1(N-1)=0.0
            END IF
            B1(N) = 0.0
      	ELSE
C           ENERGY BALANCE BASED ON TEMPERATURE
            TSP3=STEFAN*4.*WDT*(TSP(1)+273.16)**3
            IF (MATERL .GT. 2) THEN
C              ADJUST COEFFICIENTS FOR MATERIAL ABOVE SNOWPACK
               A1(N) = A1(N)*EMITSP
               C1(N-1)=C1(N-1)*EMITSP*TSP3
      	END IF
            B1(N) = -EMITSP*TSP3
      END IF
         N = N + NSP
         MATERL=MATERL + 1
C        INITIALIZE COEFFICIENTS FOR BOTTOM NODE AND UNDERLYING MATERIAL
         IF (NSP .GT. 1) B1(N-1) = 0.0
         C1(N-1) = 0.0
         DO 5 I=1,NR+1
C           SNOW OVERLYING RESIDUE AND SOIL - INITIALIZE MATRIX COEFF.
            A1(N) = 0.0
            B1(N) = 0.0
            C1(N) = 0.0
            N = N + 1
    5    CONTINUE
         GO TO 50
      END IF
C
C     RADIATION BALANCE MATRIX COEFFICIENTS FOR RESIDUE NODES
C
      IF (NR .GT. 0) THEN
         EMITNC=(1-TDIFFU(1))*EMITR
         TR3=STEFAN*4.*WDT*(TR(1)+273.16)**3
         IF (MATERL .GT. 2) THEN
C           ADJUST COEFFICIENTS FOR MATERIAL ABOVE RESIDUE
            A1(N) = A1(N)*EMITNC
            C1(N-1)=C1(N-1)*EMITNC*TR3
         END IF
         B1(N) =-2.*EMITNC*TR3
         C1(N) = EMITNC
         A1(N+1) = EMITNC*TR3
         N = N + 1
C
         DO 20 I=2,NR
            EMITNC=(1-TDIFFU(I))*EMITR
            TR3=STEFAN*4.*WDT*(TR(I)+273.16)**3
            C1(N-1)=C1(N-1)*EMITNC*TR3
            A1(N) = A1(N)*EMITNC
            B1(N) = -2.*EMITNC*TR3
            C1(N) = EMITNC
            A1(N+1) =EMITNC*TR3
            N = N + 1
   20    CONTINUE
         MATERL = MATERL + 1
      END IF
C
C     RADIATION BALANCE MATRIX COEFFICIENTS FOR SOIL SURFACE
      TS3=STEFAN*4.*WDT*(TS(1)+273.16)**3
      IF (MATERL .GT. 2) THEN
C        ADJUST COEFFICIENTS FOR MATERIAL ABOVE SOIL
         A1(N) = A1(N)*EMITS
         C1(N-1)=C1(N-1)*EMITS*TS3
      END IF
      B1(N)=-EMITS*TS3
   50 RETURN
      END
C***********************************************************************
C
      SUBROUTINE SOURCE (NC,NSP,NR,NS,NSALT,NPLANT,UC,SC,USP,SSP,UR,SR,
     >   US,SS,SINK,SWCAN,SWSNOW,SWRES,SWSOIL,LWCAN,LWSNOW,LWRES,LWSOIL,
     >   SOILXT)
C
C     THIS SUBROUTINE SUMS UP THE SOURCE-SINK TERMS FOR EACH NODE.
C
C***********************************************************************
      REAL UC(10),SC(10),USP(100),SSP(100),UR(10),SR(10),US(50),
     >     SS(50),SINK(10,50),SWCAN(9,10),SWSNOW(100),SWRES(10),SWSOIL,
     >     LWCAN(9,10),LWSNOW(2),LWRES(10),LWSOIL,SOILXT(50)
C
C**** CANOPY NODES
      DO 10 I=1,NC
         UC(I)=0.0
         SC(I)=SWCAN(NPLANT+1,I)+LWCAN(NPLANT+1,I)
   10 CONTINUE
C
C**** SNOWPACK NODES
      IF (NSP .GT. 0) THEN
         DO 20 I=1,NSP
            USP(I)=0.0
            SSP(I)=SWSNOW(I)
   20    CONTINUE
         SSP(1)=SSP(1)+LWSNOW(1)
         SSP(NSP)=SSP(NSP)+LWSNOW(2)
      END IF
C
C**** RESIDUE NODES
      DO 30 I=1,NR
         UR(I)=0.0
         SR(I)=SWRES(I)+LWRES(I)
   30 CONTINUE
C**** SOIL NODES
      DO 50 I=1,NS
         US(I)=-SOILXT(I)
         SS(I)=0.0
         DO 40 J=1,NSALT
            SINK(J,I)=0.0
   40    CONTINUE
   50 CONTINUE
      SS(1)=SS(1)+SWSOIL+LWSOIL
C
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE ATSTAB (NPLANT,NC,NSP,NR,TA,TADT,T,TDT,HUM,HUMDT,
     >           VAP,VAPDT,WIND,HEIGHT,HTOLER,HFLUX,VFLUX,ZC,ZR,RHOR,
     >           ICE,ITER)
C
C     THIS SUBROUTINE CALCULATES THE ATMOSPHERIC STABILITY, THE TRANSFER
C     COEFFICIENTS FOR THE BOTH HEAT AND VAPOR FROM THE SURFACE,
C     DEFINES THE BOUNDARY CONDITIONS FOR THE MATRIX SOLUTIONS, AND
C     CALCULATES THE WINDSPEED PROFILES IN THE RESIDUE AND CANOPY
C
C***********************************************************************
      COMMON /TIMEWT/ WT,WDT,DT
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /MATRIX/ A1(150),B1(150),C1(150),D1(150),
     >                A2(150),B2(150),C2(150),D2(150)
      COMMON /CLAYRS/ DRYCAN(8,10),CANLAI(8,10),TOTLAI(8),IEVAP(8),
     >                RLEAF(8,10),RROOT(8,50),ROOTDN(8,50),TOTROT(8)
      COMMON /WINDV/ ZH,ZM,ZERO,USTAR,STABLE,WINDC(11),WINDR(10)
      REAL LF,LS,LV
C
      REAL ZC(11),ZR(10),RHOR(10)
      SAVE TMPAIR,VAPAIR,ZMLOG,ZHLOG,HFLUX1,PSIM,PSIH
C
      IF (ITER .LE. 2) THEN
         IF (ITER .LE. 1) THEN
C****       CALCULATE CONSTANTS FOR TIME STEP
            CALL VSLOPE (DUMMY,SATV,TA)
            CALL VSLOPE (DUMMY,SATVDT,TADT)
            TMPAIR=WT*TA + WDT*TADT
            VAPAIR=WT*HUM*SATV + WDT*HUMDT*SATVDT
            ZMLOG=ALOG((HEIGHT+ZM-ZERO)/ZM)
            ZHLOG=ALOG((HEIGHT+ZH-ZERO)/ZH)
         END IF
C
C****    DEFINE INITIAL ASSUMPTIONS FOR CURRENT TIME STEP
C     HFLUX        turbulent heat flux at surface node (W/m^2)
C     USTAR        friction velocity (m/s)
         HFLUX1=0.0
         PSIM=0.0
         PSIH=0.0
C     Refer to eq.(19) in the document.
         USTAR=WIND*VONKRM/(ZMLOG + PSIM)
      END IF
C
      TMPSFC=WT*T + WDT*TDT
      VAPSFC=WT*VAP + WDT*VAPDT
C     VFLUX        turbulent vapor transfer at surface node (kg/m^2/s)
      HFLUX=0.0
      VFLUX=0.0
      IF (WIND .LE. 0.0) GO TO 20
C
      ITER1=0
C
C**** START ITERATIVE PROCESS TO OBTAIN HEAT FLUX
   10 ITER1=ITER1+1
C     Refer to Eq.(21) In the Document
      STABLE=VONKRM*(HEIGHT-ZERO)*G*HFLUX1
     >       /(RHOA*CA*(TMPAIR+273.16)*USTAR**3)
C
      IF (STABLE .GE. 0.0) THEN
C****    ATMOSPHERE IS STABLE
C        IF STABILITY IS GREATER THAN ZMLOG/9.4, COMPUTED HFLUX WILL
C        ACTUALLY DECREASE WITH INCREASING TEMPERATURE DIFFERENTIAL,
C        CAUSING NUMERICAL INSTABILITY AND CONVERGENCE PROBLEMS -- AVOID
C        THIS SITUATION BY LIMITING STABLE TO ZMLOG/9.4
         IF (STABLE .GT. ZMLOG/9.4) STABLE=ZMLOG/9.4
         PSIH=4.7*STABLE
         PSIM=PSIH
      ELSE
C****    ATMOSPHERE IS UNSTABLE
         PSIH=-2.*ALOG((1 + (1-16.*STABLE)**0.5)/2.0)
         PSIM=0.6*PSIH
C        IF ATMOSPHERE IS SUFFICIENTLY UNSTABLE, EQUATIONS MAY RESULT
C        IN FRICTION VELOCITY LESS THEN ZERO !!??!!?? (LIMIT PSIM)
         IF (PSIM/ZMLOG .LT. -0.50) THEN
            PSIM=-0.50*ZMLOG
            PSIH=PSIM/0.6
            STABLE = -((2.*EXP(-PSIH/2.) - 1)**2 -1)/16.
         END IF
      END IF
      USTAR=WIND*VONKRM/(ZMLOG + PSIM)
      RH=(ZHLOG + PSIH)/(USTAR*VONKRM)
      HFLUX2=RHOA*CA*(TMPAIR-TMPSFC)/RH
      ERROR=ABS(HFLUX2-HFLUX1)
      HFLUX1=HFLUX2
      IF(ABS(HFLUX1).LT.ABS(TKA*(TMPAIR-TMPSFC)/(HEIGHT-ZERO))) GO TO 15
C        WIND IS SUFFICIENTLY SMALL, OR ATMOSPHERE IS SUFFICIENTLY
C        STABLE THAT TURBULENCE DOES NOT INCREASE FLUX -- NO NEED TO
C        ITERATE FURTHER
C
C**** CHECK IF TOLERANCE HAS BEEN MET
      IF (ITER1 .GT. 40) THEN
C        CONVERGENCE NOT MET, BUT PROBABLY CLOSE ENOUGH
         GO TO 15
      END IF
      IF (ERROR .GT. HTOLER) GO TO 10
C
C**** HEAT FLUX IS WITHIN TOLERABLE ERROR -- NOW CALCULATE VAPOR FLUX
   15 HFLUX=HFLUX1
      RV=RH
      VFLUX=(VAPAIR-VAPSFC)/RV
C
C**** COMPARE FLUX WITH MINIMUM FLUX AS CALCULATED FROM THERMAL
C     CONDUCTIVITY AND VAPOR DIFFUSIVITY OF STILL AIR -- HEAT FLUX FIRST
   20 HMIN=TKA*(TMPAIR-TMPSFC)/(HEIGHT-ZERO)
      IF (ABS(HMIN) .LE. ABS(HFLUX)) THEN
C        WIND IS SUFFICIENT THAT HEAT FLUX IS ENHANCED
         CON=RHOA*CA/RH
        ELSE
C        WIND IS SUFFICIENTLY LOW TO BE CONSIDERED AS "STILL AIR"
         HFLUX=HMIN
         CON=TKA/(HEIGHT-ZERO)
      END IF
C
C**** NOW COMPARE VAPOR FLUXES
      DV=VDIFF*(((TMPAIR+273.16)/273.16)**2)*(P0/PRESUR)
      VMIN=DV*(VAPAIR-VAPSFC)/(HEIGHT-ZERO)
      IF (ABS(VMIN) .LE. ABS(VFLUX)) THEN
C        WIND IS SUFFICIENT THAT VAPOR FLUX IS ENHANCED
         CONV=1/RV
       ELSE
C        WIND IS SUFFICIENTLY LOW TO BE CONSIDERED AS "STILL AIR"
         VFLUX=VMIN
         CONV=DV/(HEIGHT-ZERO)
       END IF
C
C**** DEFINE MATRIX COEFFICIENTS FOR SURFACE MATERIAL PRESENT
C
      IF (NC .GT. 0) THEN
C****    SURFACE MATERIAL IS CANOPY
         B1(1)=B1(1) - WDT*CON
         D1(1)=HFLUX
         B2(1)=-WDT*CONV
         D2(1)=VFLUX
         WINDCH=USTAR*ALOG((ZC(NC+1)+ZM-ZERO)/ZM)/VONKRM
C        CALCULATE WINDSPEED AT THE CANOPY NODES ASSUMING AN
C        EXPONENTIAL DECREASE IN WINDSPEED WITH DEPTH
C        AWIND = 0  FOR SPARSE CANOPY: AWIND >= 4 FOR DENSE CANOPY
         SUMLAI=0.0
         DO 40 J=1,NPLANT
            SUMLAI = SUMLAI + TOTLAI(J)
   40    CONTINUE
         AWIND = SUMLAI
         DO 45 I=1,NC+1
            WINDC(I) = WINDCH*EXP(-AWIND*ZC(I)/ZC(NC+1))
   45    CONTINUE
         IF (NR.GT.0 .AND. NSP.EQ.0) THEN
C           RESIDUE LIES BENEATH THE CANOPY -- CALCULATE WINDRH
            WINDRH = WINDC(NC+1)
           ELSE
C           NO RESIDUE LIES BENEATH CANOPY -- PERHAPS IT IS SNOW-COVERED
            WINDRH=0.0
         END IF
      ELSE
C
C
C           ==================================
      IF (NSP .GT. 0) THEN
C****    SURFACE MATERIAL IS SNOW
         CALL VSLOPE (S,DUMMY,TDT)
         IF (ICE .EQ. 0) THEN
C           SNOW IS NOT MELTING - ENERGY BALANCE BASED ON TEMPERATURE
C     LS           latent heat of sublimation
            B1(1) = B1(1) - WDT*CON - WDT*LS*S*CONV
          ELSE
C           SNOW IS MELTING - ENERGY BALANCE BASED ON LIQUID WATER
            B1(1) = B1(1)
         END IF
         D1(1)=HFLUX + LS*VFLUX
         B2(1)=0.0
         D2(1)=0.0
C        IF RESIDUE LIES BENEATH THE SNOW, WINDRH = 0.0
         WINDRH=0.0
      ELSE
C
      IF (NR .GT. 0) THEN
C****    SURFACE MATERIAL IS RESIDUE
         B1(1)=B1(1) - WDT*CON
         D1(1)=HFLUX
         B2(1)=-WDT*CONV
         D2(1)=VFLUX
C        CALCULATE WINDSPEED AT THE TOP OF THE RESIDUE LAYER
         WINDRH=USTAR*ALOG((ZR(NR+1)+ZM-ZERO)/ZM)/VONKRM
      ELSE
C
C**** SURFACE MATERIAL IS BARE SOIL
         CALL VSLOPE (S,DUMMY,TDT)
         D1(1)=HFLUX+LV*VFLUX
         D2(1)=VFLUX/RHOL
C        CHECK IF ICE IS PRESENT -- IF SO, THE WATER BALANCE IS WRITTEN
C        IN TERMS OF ICE CONTENT.  THE ENERGY BALANCE IS WRITTEN WITH
C        LIQUID WATER AS THE REFERENCE STATE, THEREFORE LATENT HEAT OF
C        SUBLIMATION IS NEVER USED -- ONLY HEAT OF VAPORIZATION.
         IF (ICE .EQ. 0) THEN
C           NO ICE IS PRESENT AT SOIL SURFACE
            B1(1) = B1(1) - WDT*CON - WDT*LV*S*CONV
            B2(1) = -WDT*VAPDT*CONV*0.018*G/(UGAS*(TDT+273.16))/RHOL
           ELSE
C           ICE IS PRESENT AT SOIL SURFACE
            B1(1) = B1(1)
            B2(1) = 0.0
         END IF
C
      END IF
      END IF
      END IF
C
      IF (NR .GT. 0) THEN
C        CALCULATE THE WINDSPEED AT THE MID-POINT OF EACH RESIDUE 
C        LAYER ASSUMING AN EXPONENTIAL DECREASE WITH DEPTH
C        AWIND = 0  FOR SPARSE RESIDUE: AWIND >= 4 FOR DENSE RESIDUE 
         AWIND = 0.1*(RHOR(1) - 20.)
         IF (AWIND .LT. 0.0) AWIND=0.0
         WINDR(1)=WINDRH*EXP(-AWIND*(ZR(1)+ZR(2))/4/ZR(NR+1))
         DO 65 I=2,NR+1
            WINDR(I) = WINDRH*EXP(-AWIND*ZR(I)/ZR(NR+1))
   65    CONTINUE
      END IF
C
      RETURN
C
      END
C***********************************************************************
C
      SUBROUTINE EBCAN (N,NPLANT,NC,NSP,NR,NS,ZC,TC,TCDT,VAPC,VAPCDT,
     >  WCAN,WCANDT,PCAN,PCANDT,MAT,MATDT,SWCAN,LWCAN,CANMA,CANMB,
     >  DCHAR,RSTOM0,RSTEXP,PLEAF0,RLEAF0,ITYPE,ITER)
C
C     THIS SUBROUTINE CALCULATES THE JACOBIAN MATRIX COEFFICIENTS FOR
C     THE CANOPY PORTION OF THE NEWTON-RAPHSON SOLUTION OF THE ENERGY
C     BALANCE
C
C***********************************************************************
      COMMON /TIMEWT/ WT,WDT,DT
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /MATRIX/ A1(150),B1(150),C1(150),D1(150),
     >                A2(150),B2(150),C2(150),D2(150)
      REAL LF,LS,LV
C
      REAL ZC(11),TC(11),TCDT(11),VAPC(11),VAPCDT(11),
     >     WCAN(10),WCANDT(10),PCAN(8),PCANDT(8),MAT(50),MATDT(50),
     >     SWCAN(9,10),LWCAN(9,10)
      REAL DCHAR(8),RSTOM0(8),RSTEXP(8),PLEAF0(8),RLEAF0(8)
      REAL CON(10),CONT(10),CONDT(10),
     >     TRNSP(9),XTRACT(50),HEATC(10),ETLYR(10),DETLYR(10),
     >     DHEATC(10),DTLDTC(10)
      INTEGER ITYPE(8)
      SAVE CONT
C
C**** CALCULATE LEAF TEMPERATURE AND HEAT TRANSFER FROM CANOPY
      CALL LEAFT (NPLANT,NC,NS,ITER,ITYPE,TC,TCDT,VAPC,VAPCDT,
     > WCAN,WCANDT,PCAN,PCANDT,MAT,MATDT,TRNSP,XTRACT,SWCAN,LWCAN,
     > HEATC,ETLYR,DETLYR,CANMA,CANMB,DCHAR,RSTOM0,RSTEXP,PLEAF0,RLEAF0,
     > DHEATC,DTLDTC)
C
C**** DETERMINE THE EDDY CONDUCTANCE TERM BETWEEN NODES
      IF (ITER .EQ. 1) CALL CANTK (NC,CONT,TC,ZC)
      CALL CANTK (NC,CONDT,TCDT,ZC)
C
C**** CALCULATE THE AVERAGE CONDUCTANCE TERM OVER THE TIME STEP
      CALL WEIGHT (NC,CON,CONT,CONDT)
C
C**** DETERMINE THE MATRIX COEFFICIENTS FOR THE TOP LAYER
C
      A1(N+1)=A1(N+1) + WDT*CON(1)
      B1(N)=B1(N) - WDT*(CON(1)+DHEATC(1)*(1.-DTLDTC(1)))
     >            - (ZC(2)-ZC(1))/(2.*DT)*RHOA*CA
      C1(N)=C1(N) + WDT*CON(1)
      D1(N)=D1(N) - CON(1)*(WT*(TC(1)-TC(2))+WDT*(TCDT(1)-TCDT(2)))
     >    + HEATC(1) -(ZC(2)-ZC(1))/(2.*DT)*RHOA*CA*(TCDT(1)-TC(1))
C
C**** DETERMINE THE COEFFICIENTS FOR THE REST OF THE LAYERS
      DO 10 I=N+1,N+NC-1
         J=I-N+1
         A1(I+1)=A1(I+1) + WDT*CON(J)
         B1(I)=B1(I) - WDT*(CON(J-1) +CON(J) +DHEATC(J)*(1.-DTLDTC(J)))
     >               - (ZC(J+1)-ZC(J-1))/(2.*DT)*RHOA*CA
         C1(I)=C1(I) + WDT*CON(J)
         D1(I)=CON(J-1)*(WDT*(TCDT(J-1)-TCDT(J)) + WT*(TC(J-1)-TC(J)))
     >        -CON(J)*(WDT*(TCDT(J)-TCDT(J+1)) + WT*(TC(J)-TC(J+1)))
     >        -(ZC(J+1)-ZC(J-1))/(2.*DT)*RHOA*CA*(TCDT(J)-TC(J))
     >        +HEATC(J)
   10 CONTINUE
      N=N+NC
C
C**** DETERMINE THE COEFFICIENTS FOR THE TOP LAYER OF THE NEXT MATERIAL
      IF (NSP.GT.0 .OR. NR.EQ.0) THEN
C        NEXT MATERIAL IS SNOW OR SOIL -- NEED TO CONSIDER LATENT HEAT
C        TRANSFER TO THE NEXT NODE.
C        DETERMINE THE CONVECTIVE VAPOR TRANSPORT FROM EDDY CONDUCTANCE
C     CONV (I)     vapor conductance between node i and i+1
         CONV=CON(NC)/RHOA/CA
C        CALCULATE THE VAPOR TRANSFER BETWEEN NODES
         QVCAN=CONV*(WDT*(VAPCDT(NC)-VAPCDT(NC+1))
     >                 + WT*(VAPC(NC)-VAPC(NC+1)))
C        CALCULATE HUMIDITY AND SLOPE OF SAT. VAPOR CURVE
         CALL VSLOPE (SLPNC,SATV,TC(NC))
         HUMNC=VAPCDT(NC)
         CALL VSLOPE (SLPNC1,SATV,TC(NC+1))
         HUMNC1=VAPCDT(NC+1)/SATV
      END IF
C
      IF (NSP.GT.0) THEN
C****    NEXT MATERIAL IS SNOW
         A1(N) = A1(N) + WDT*CONV*LS*SLPNC*HUMNC
         B1(N) = B1(N) - WDT*(CON(NC) + CONV*LS*SLPNC1*HUMNC1)
         D1(N) =CON(NC)*(WT*(TC(NC)-TC(NC+1))+WDT*(TCDT(NC)-TCDT(NC+1)))
     >          + LS*QVCAN
        ELSE
C
      IF (NR.GT.0) THEN
C****    NEXT MATERIAL IS RESIDUE
         B1(N) = B1(N) - WDT*CON(NC)
         D1(N) =CON(NC)*(WT*(TC(NC)-TC(NC+1))+WDT*(TCDT(NC)-TCDT(NC+1)))
        ELSE
C
C**** NEXT MATERIAL IS SOIL
      A1(N) = A1(N) + WDT*CONV*LV*SLPNC*HUMNC
      B1(N) = B1(N) - WDT*(CON(NC) + CONV*LV*SLPNC1*HUMNC1)
      D1(N) = CON(NC)*(WT*(TC(NC)-TC(NC+1)) + WDT*(TCDT(NC)-TCDT(NC+1)))
     >       + LV*QVCAN
C
      END IF
      END IF
C
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE LEAFT (NPLANT,NC,NS,ITER,ITYPE,TC,TCDT,VAPC,VAPCDT,
     > WCAN,WCANDT,PCAN,PCANDT,MAT,MATDT,TRNSP,XTRACT,SWCAN,LWCAN,
     > HEATC,ETLYR,DETLYR,CANMA,CANMB,DCHAR,RSTOM0,RSTEXP,PLEAF0,RLEAF0,
     > DHEATC,DTLDTC)
C
C     THIS SUBROUTINE COMPUTES LEAF TEMPERATURE OF EACH CANOPY TYPE AND
C     THE TOTAL HEAT AND WATER TRANSFERRED FROM CANOPY TO SURROUNDING
C     AIR SPACE IN EACH CANOPY LAYER
C
C***********************************************************************
C
C     ETLYR (I)    evapotranspiration within canopy layer i
C     DETLYR (I)   average conductance (stomatal plus leaf boundary layer) for transpiration
C             in canopy layer i, which is equal to the derivative (m/s)
C     HEATC (I)    heat flux term from plants leaves within canopy layer i
C     TRNSP (J)    total water transpired by plant j for the time step (kg)
C     RSTEXP (J)   exponent for calculating stomatal resistance
C     DHEATC (I)   conductance term for heat flux term from plants leaves within canopy
C             layer i, which is equal to derivative with respect to temperature
C     DTLDTC (I)   derivative of canopy leaf temperature to canopy air temperature

      COMMON /TIMEWT/ WT,WDT,DT
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /CLAYRS/ DRYCAN(8,10),CANLAI(8,10),TOTLAI(8),IEVAP(8),
     >                RLEAF(8,10),RROOT(8,50),ROOTDN(8,50),TOTROT(8)
      COMMON /WINDV/ ZH,ZM,ZERO,USTAR,STABLE,WINDC(11),WINDR(10)
      REAL LF,LV,LS
C
      REAL TC(11),TCDT(11),VAPC(11),VAPCDT(11),WCAN(10),WCANDT(10),
     >     PCAN(8),PCANDT(8),MAT(50),MATDT(50),TRNSP(9),XTRACT(50),
     >     SWCAN(9,10),LWCAN(9,10),HEATC(10),ETLYR(10),
     >     DETLYR(10),DHEATC(10),DTLDTC(10)
      REAL DCHAR(8),RSTOM0(8),RSTEXP(8),PLEAF0(8),RLEAF0(8)
C     AVGTMP (I)   average temperature of layer for time step (C)

      REAL AVGTMP(10),AVGVAP(10),RSTOM(10),PLEAF(10),PEVAP(10),
     >     PXYLEM(8),RHCAN(8,10),ETCAN(8,10),TLC(8,10),AVGMAT(50),
     >     HUMID(10),VTSLOP(10)
      REAL F1(10),DF1DP(10),DF1DT(10),DF1DX(10),DP(10),AA1(10),CC1(10),
     >     F2(10),DF2DP(10),DF2DT(10),DF2DX(10),DTLC(10),DF3DP(10)
      INTEGER INIT(8),ITYPE(8)
      SAVE INIT,PXYLEM,RHCAN,ETCAN,TLC
C
C
C     INITIALIZE ROOT EXTRACTION AND COMPUTE AVERAGE MATRIC POTENTIAL
C     XTRACT(I)   water extracted by roots for soil layer i (m/s)
C     初始化根系提水
      DO 5 I=1,NS
         XTRACT(I)=0.0
         AVGMAT(I)=WT*MAT(I)+WDT*MATDT(I)
C        LIMIT THE WATER POTENTIALS THAT THE PLANT SEES 
         IF (AVGMAT(I) .GT. 0.0) AVGMAT(I) = 0.0
    5 CONTINUE
C
C     TRNSP (J)    total water transpired by plant j for the time step (kg)
C     INITIALIZE THE TOTAL TRANSP. FROM THE ENTIRE TRANSPIRING CANOPY
C     初始化冠层呼吸蒸腾
      TRNSP(NPLANT+1)=0.0
C
C     初始化灌层的温度、水汽密度、热量和水分通量
C     HEATC (I)    heat flux term from plants leaves within canopy layer i
C     INITIALIZE CANOPY TEMP, VAPOR DENSITY, AND HEAT AND WATER FLUXES.
C     HEATC(I) AND ETLYR ARE HEAT AND WATER FLUXES; DHEATC AND DETLYR
C     ARE DERIVATIVES OF FLUX TERMS.

C     AVGTMP (I)   average temperature of layer for time step (C)
      DO 10 I=1,NC
         AVGTMP(I)=WT*TC(I) + WDT*TCDT(I)     !一个步长中冠层空气的平均温度
         AVGVAP(I)=WT*VAPC(I) + WDT*VAPCDT(I) !一个步长中冠层空气的平均水汽
         HEATC(I)=0.0
         DHEATC(I)=0.0  !第i冠层叶片向空气传输热量的导热系数
         ETLYR(I)=0.0   !第i冠层蒸散发量
         DETLYR(I)=0.0  !叶片向空气输水的传导系数
C     DTLDTC (I)   derivative of canopy leaf temperature to canopy air temperature
         DTLDTC(I)=0.0  !叶片温度对冠层气温的偏导
C        INITIALIZE RESISTANCE TO TRANSPORT FROM CANOPY LEAVES
C        (NOT NECESSARY IF ALREADY CALCULATED THIS TIME STEP)
         IF (ITER .EQ. 1) THEN
            DO 8 J=1,NPLANT
               RHCAN(J,I)=307.*SQRT(DCHAR(J)/WINDC(I))
               INIT(J)=1
    8       CONTINUE
         END IF
   10 CONTINUE
C
C
      DO 60 J=1,NPLANT
      	SUMET=0.0
      	SUMEV = 0.0
C     "FRACTN" IS FRACTION OF TIME STEP PLANTS WILL TRANSPIRE
      	FRACTN = 1.0
C
      	IF (TOTLAI(J) .EQ. 0.0) THEN
C           NO LEAF AREA FOR THIS PLANT(PERHAPS DORMANT OR SNOW-COVERED)
            TRNSP(J)=0.0
            GO TO 60
      	END IF
C
      IF (ITYPE(J) .NE. 0) THEN
C
C********   TRANSPIRING PLANT - CHECK IF CONDITIONS ARE SUCH THAT PLANT
C           WILL TRANSPIRE

C     PCANDT (J)   intercepted precipitation on plant species j at time t+dt (m)
C     IEVAP (J)    code for if plant j will transpire this time step (0=no; 1=yes)
C     PCAN (J)     intercepted precipitation on plant species j at time t (m)

            PCANDT(J)=0.0  !t+dt时刻冠层i雨水截留
      		IF (IEVAP(J) .EQ. 0 .OR. PCAN(J) .GT. 0.0) THEN
            !冠层不呼吸或者冠层有雨水截留
C*****         PLANT ISN'T TRANSPIRING - PERHAPS NO SUNLIGHT, TOO COLD,
C              OR INTERCEPTED PRECIP AVAILABLE ON PLANT LEAVES
C              "FRACTN" IS FRACTION OF TIME STEP PLANTS WILL TRANSPIRE
               FRACTN=0.0
C              CALCULATE LEAF TEMPERATURE
               DO 12 I=1,NC
C                 CHECK IF PLANT HAS ANY LEAF AREA IN LAYER
                  IF (CANLAI(J,I) .LE. 0.0) GO TO 12
                  RSTOM(I)=0.0  !冠层气孔导度（阻抗）
                  ETCAN(J,I)=0.0
                  IF (PCAN(J) .GT. 0.0) THEN
C***                 INTERCEPTED PRECIP AVAILABLE FOR EVAPORATION
C     HUMID        relative humidity (decimal)
                     HUMID(I)=1.0
                     CALL VSLOPE (VTSLOP(I),SATV,AVGTMP(I))
                     VAPDEF = CANLAI(J,I)*(SATV-AVGVAP(I))
                     TLC(J,I)=AVGTMP(I) + RHCAN(J,I)
     >                     *(SWCAN(J,I)+LWCAN(J,I)-LV*VAPDEF/RHCAN(J,I))
     >                     /CANLAI(J,I)/(LV*VTSLOP(I)+RHOA*CA)
C                    DETERMINE AMOUNT OF INTERCEPTED PRECIP THAT
C                    EVAPORATED FROM PLANT SURFACES
C                    PEVAP雨水蒸发
C                    SUMEV总蒸散发
                     PEVAP(I)=(CANLAI(J,I)*VTSLOP(I)
     >                        *(TLC(J,I)-AVGTMP(I))+VAPDEF)/RHCAN(J,I)
                     SUMEV=SUMEV+PEVAP(I)
                    ELSE
C***                 NO WATER AVAILABLE FOR EVAPORATION
                     HUMID(I)=0.0
                     VTSLOP(I)=0.0
                     PEVAP(I)=0.0
                     SUMEV=0.0
                     TLC(J,I)=AVGTMP(I) + RHCAN(J,I)*
     >                     (SWCAN(J,I)+LWCAN(J,I))/(CANLAI(J,I)*RHOA*CA)
                  END IF
   12          CONTINUE

      			IF (PCAN(J) .GT. 0.0) THEN
       !调整冠层雨水截留
C*****            CALCULATE WATER ON PLANTS AT END OF TIME STEP
                  PCANDT(J)=PCAN(J)-DT*SUMEV/RHOL
      				IF (PCANDT(J) .LT. 0.0) THEN
C                    NO WATER REMAINING ON PLANTS - COMPUTE FRACTION OF
C                    TIME STEP PLANTS WILL BE TRANSPIRING AND ADJUST
C                    AMOUNT EVAPORATED
                     PCANDT(J)=0.0
                     FRACTN= 1 - (RHOL*PCAN(J)/DT)/SUMEV
                     SUMEV = (1.-FRACTN)*SUMEV
                     DO 13 I=1,NC
                        PEVAP(I)=PEVAP(I)*(1.-FRACTN)
   13                CONTINUE
      				END IF
      			END IF
      		END IF
C
            IF (PCANDT(J) .LE. 0.0 .AND. IEVAP(J) .NE. 0) THEN
            !t至t+dt内，有冠层呼吸
C********      PLANT IS TRANSPIRING
C
C              FIND THE EXTREMES OF MATRIC POTENTIAL SEEN BY ROOTS
C              TAKING INTO ACCOUNT ROOT RESISTANCE FOR EACH LAYER
               MAX=0
               DO 14 I=1,NS
                  IF (ROOTDN(J,I).GT.0.0) THEN
                    IF (MAX .EQ. 0) THEN
                       MAX=I
                       MIN=I
                     ELSE
                       IF (AVGMAT(I).LT.AVGMAT(MIN)) MIN=I
                       IF (AVGMAT(I)/RROOT(J,I) .GT.
     >                     AVGMAT(MAX)/RROOT(J,MAX)) MAX=I
                    END IF
                  END IF
   14          CONTINUE
C   
C              DETERMINE LEAF POTENTIAL AND TEMP、
C              确定植被叶片水势和温度
               VAPDEF=0.0
               RHAVG=0.0
C              INITIZE VARIABLES FOR CANOPY
               DO 15 I=1,NC
C                 CHECK IF PLANT HAS ANY LEAF AREA IN THIS LAYER
                  IF (CANLAI(J,I) .GT. 0.0) THEN
                     IF (FRACTN .GT. 0.999) PEVAP(I)=0.0
                     TLC(J,I) = AVGTMP(I)
                     HUMID(I)=1.0
C                    CALCULATE AVERAGE CONDITIONS IN CANOPY FOR AN 
C                    INTIIAL APPROXIMATION TO PLEAF AND PXYLEM
                     CALL VSLOPE (VTSLOP(I),SATV,TLC(J,I))
                     VAPDEF = VAPDEF + CANLAI(J,I)*(SATV-AVGVAP(I))
                     RHAVG = RHAVG + CANLAI(J,I)/RHCAN(J,I)
                  END IF
   15          CONTINUE
               VAPDEF=VAPDEF/TOTLAI(J)!植被j平均单位叶面积上的水汽压差
               IF (VAPDEF .LT. 0.0) VAPDEF=0.0
               RHAVG=TOTLAI(J)/RHAVG
C
C*****         BEGIN ITERATION TO FIND INITIAL ESTIMATE FOR PXYLEM
               ITER0=0
C              CALCULATE INITIAL GUESS FOR PXYLEM IF THIS IS FIRST TIME
C              FOR THIS TIME STEP
   18          IF (INIT(J) .EQ. 1) PXYLEM(J)=2.*AVGMAT(MAX)
C***           COMPUTE SUM OF ROOT CONDUCTANCE TIMES MATRIC POTENTIAL
   19          RSOIL=0.0 !植被根区导水系数乘基质势
               SRROOT=0.0!植被根区导水率
               IROOT=0
               DO 20 I=1,NS
C                 DO NOT CONSIDER IF NO PLANT J ROOTS IN SOIL LAYER
                  IF (ROOTDN(J,I).GT. 0.0) THEN
C                    DO NOT INCLUDE SOIL LAYERS DRYER THAN PLANT XYLEM
C     RROOT (J,I)  root conductivity of plant j in canopy layer i (kg/m^3/s)
C     RROOT0 (J)   overall root conductivity for plant j (kg/m^3/s)
                     IF (AVGMAT(I) .GE. PXYLEM(J)) THEN
                        RSOIL=RSOIL + RROOT(J,I)*AVGMAT(I)
                        SRROOT= SRROOT + RROOT(J,I)
                       ELSE
                        IROOT=1
                     END IF
                  END IF
   20          CONTINUE
               IF (INIT(J) .GT. 1) THEN
C                 USE VALUES FOR PXYLEM, ETCAN AND TLC FROM PREVIOUS
C                 CALCULATIONS FOR THIS TIME STEP IF AVAILABLE
C                 (CALCULATE SUMET AND GO DIRECTLY TO ITERATIVE SCHEME)
                  SUMET=RSOIL - PXYLEM(J)*SRROOT
                  IF (SUMET*SRROOT .LE. 0.0) THEN
C                    PREVIOUS VALUE OF PXYLEM WILL NOT WORK FOR UPDATED
C                    END-OF-TIME-STEP CONDITIONS
                     INIT(J) = 1
                     GO TO 18
                  END IF
                  GO TO 24
               END IF
C***           CALC. EFFECTIVE MATRIC POT. AND TOTAL RESISTANCE OF PLANT
C             计算有效基质势以及植被的总阻抗
               IF (SRROOT .EQ. 0.0) THEN
                  RSOIL=RROOT(J,MAX)*AVGMAT(MAX)
                  SRROOT=RROOT(J,MAX)
               END IF
               SOIMAT=RSOIL/SRROOT !土壤有效基质势（针对植被）
               RESIST= 1/SRROOT + 1/RLEAF0(J)!总阻抗
               IF (ITER0 .EQ. 0) THEN
C                 ESTIMATE STARTING POINT FOR ITERATION
C     PLEAF0 (J)   critical leaf potential for plant j (m)
C     RSTOM0 (J)   stomatal resistance of plant j with no water stress (s/m)
                  PLEAF1=PXYLEM(J)
                  IF (PLEAF1/PLEAF0(J) .GT. 40.) THEN
C                    LIKELIHOOD OF ARITHMETIC OVERFLOW -- PROCEED WITH
C                    CAUTION BY TAKING LOGARITHM
                     RSLOG=ALOG10(RSTOM0(J))
     >                     + RSTEXP(J)*ALOG10(PLEAF1/PLEAF0(J))!eq.(50)
                     IF (RSLOG .GT. 20) THEN
C                       LOG OF STOMATAL RESISTANCE EXTREMELY LARGE --
C                       TRANSP IS ESSENTIALLY ZERO - CALC LEAF TEMP
                        SUMET=0.0
                        DO 21 I=1,NC
                           IF (CANLAI(J,I) .LE. 0.0) GO TO 21
                           HUMID(I)=0.0
                           RSTOM(I)=1.0E20
                           ETCAN(J,I)=0.0
                           TLC(J,I)=AVGTMP(I) + (SWCAN(J,I)+LWCAN(J,I))
     >                                 *RHCAN(J,I)/(CANLAI(J,I)*RHOA*CA)
   21                   CONTINUE
                        GO TO 48
                     END IF
                  END IF
                  RSTOM1=RSTOM0(J)*(1. + (PLEAF1/PLEAF0(J))**RSTEXP(J))
                  SUMET=TOTLAI(J)*VAPDEF/(RSTOM1+RHAVG)
                  PLEAF1=SOIMAT-SUMET*RESIST/2.
               END IF
C***           UPDATE STOMATAL RESISTANCE AND TRANSPIRATION
   22          RSTOM1=RSTOM0(J)*(1. + (PLEAF1/PLEAF0(J))**RSTEXP(J))
               SUMET=TOTLAI(J)*VAPDEF/(RSTOM1+RHAVG)
C              CALCULATE ERROR IN ET ESTIMATE, DERIVATIVE WITH RESPECT
C              TO LEAF POTENTIAL, AND NEW APPROX. TO LEAF POTENTIAL
               ERROR = (SOIMAT-PLEAF1)/RESIST - SUMET
               DERIV = -1/RESIST + SUMET*RSTOM0(J)*RSTEXP(J)
     >                            *(PLEAF1/PLEAF0(J))**(RSTEXP(J)-1)
     >                            /PLEAF0(J)/(RSTOM1+RHAVG)
               DELTA=ERROR/DERIV
C***           DEPENDING ON MAGNITUDE OF RSTEXP, A DRASTIC POINT OF
C              INFLECTION OCCURS IN THE ERROR FUNCTION AT PLEAF0. IF
C              UPDATED PLEAF1 CROSSES THIS POINT, CUT DELTA IN HALF
               ANEG=(PLEAF1-PLEAF0(J))*(PLEAF1-DELTA-PLEAF0(J))
               IF (ANEG .LT. 0.0) THEN
                  PLEAF1=PLEAF1-DELTA/2.
                 ELSE
                  PLEAF1=PLEAF1-DELTA
               END IF
C              CALCULATE UPDATED ET AND XYLEM POTENTIAL
               SUMET = (SOIMAT-PLEAF1)/RESIST
               PXYLEM(J)=(RSOIL-SUMET)/SRROOT
               IF (ABS(PLEAF1) .LT. 1.0) THEN
C                 AVOID DIVISION BY ZERO
                  COMPAR = 1.0
                 ELSE
                  COMPAR = PLEAF1
               END IF                 
               IF (ABS(DELTA/COMPAR).GT.0.01 .AND. ITER0.LE.20) THEN
C***              PLEAF AND PXYLEM NOT CLOSE ENOUGH
                  ITER0=ITER0+1
C                 IF PXYLEM > MINIMUM SOIL POTENTIAL, RECALCULATE
C                 RSOIL, SRROOT AND PXYLEM TO EXCLUDE DRY LAYERS
                  IF (PXYLEM(J).GT.AVGMAT(MIN) .OR. IROOT.GT.0) GO TO 19
                  GO TO 22
               END IF
C*****         ESTIMATE TRANSP (ETCAN) WITHIN EACH LAYER TO BEGIN ITER
               DO 23 I=1,NC
                  ETCAN(J,I)=CANLAI(J,I)*SUMET/TOTLAI(J)
   23          CONTINUE
               INIT(J) = 2
C
C*****         BEGIN ITERATION TO FIND LEAF TEMPERATURE, LEAF POTENTIAL
C              AND TRANSPIRATION FROM EACH CANOPY LAYER FOR PLANT
   24          ITER1 = 0
               DO 25 I=1,NC
C                 INTIAL ESTIMATE OF LEAF POTENTIAL IN EACH LAYER
C                 AND DEFINE NEWTON-RAPHSON COEFF. THAT ARE CONSTANT
                  IF (CANLAI(J,I) .GT. 0.0) THEN
                     PLEAF(I)=PXYLEM(J) - ETCAN(J,I)/RLEAF(J,I)
                     DF1DX(I) = -RLEAF(J,I)
                     DF2DT(I) = -CANLAI(J,I)*RHOA*CA/RHCAN(J,I)
                  END IF
   25          CONTINUE
C
   26          IFLAG= 0
               SUMET=0.0
               DF3DX=0.0
C              SET UP COEFFICIENTS FOR NEWTON RAPHSON SOLUTION
               DO 30 I=1,NC
                  IF (CANLAI(J,I) .GT. 0.0) THEN
                  CALL VSLOPE (VTSLOP(I),SATV,TLC(J,I))
                  RSTOM(I)=RSTOM0(J)*(1. +
     >                               (PLEAF(I)/PLEAF0(J))**RSTEXP(J))
                  F1(I) =CANLAI(J,I)*(SATV-AVGVAP(I))
     >                    /(RSTOM(I)+RHCAN(J,I))
                  IF (F1(I) .LT. 0.0) THEN
C***                 NO TRANSPIRATION
                     VTSLOP(I)=0.0
                      ETCAN(J,I)=0.0
                     DF1DP(I) = RLEAF(J,I)
                     DF1DT(I) = 0.0
                     DF2DP(I) = 0.0
                     DF2DX(I) = 0.0
                     DF3DP(I) = 0.0
C                    FORCE LEAF POTENTIAL EQUAL TO PXYLEM POTENTIAL,
C                    I.E. FORCE TRANSPIRATION IN LAYER TO ZERO
                     F1(I) = 0.0
                     PLEAF(I) = PXYLEM(J)
                    ELSE
C***                 CALCULATE TRANPIRATION IN EACH LAYER AND SET UP 
C***                 MATRIX FOR NEWTON-RAPHSON APPROX. OF UPDATED XYLEM 
C                    POTENTIAL, LEAF TEMPERATURE AND LEAF POTENTIAL
                     ETCAN(J,I)=RLEAF(J,I)*(PXYLEM(J)-PLEAF(I))
                     SUMET = SUMET +ETCAN(J,I)
                     DF1DP(I) = RLEAF(J,I) - F1(I)*RSTOM0(J)*RSTEXP(J)
     >                              *(PLEAF(I)/PLEAF0(J))**(RSTEXP(J)-1)
     >                              /PLEAF0(J)/(RSTOM(I)+RHCAN(J,I))
                     DF1DT(I) = CANLAI(J,I)*VTSLOP(I)
     >                          /(RSTOM(I)+RHCAN(J,I))
                     DF2DP(I) = LV*RLEAF(J,I)
                     DF2DX(I) = -DF2DP(I)
                     DF3DP(I) = RLEAF(J,I)
                     DF3DX = DF3DX - RLEAF(J,I)
                     F1(I) = F1(I) - ETCAN(J,I)
                  END IF
                  F2(I) = -LV*ETCAN(J,I) + DF2DT(I)*(TLC(J,I)-AVGTMP(I))
     >                 + SWCAN(J,I) + LWCAN(J,I)
                  END IF
   30          CONTINUE
               F3 = RSOIL - PXYLEM(J)*SRROOT - SUMET
               DF3DX = DF3DX - SRROOT
               CC3=DF3DX
C
               IF (SUMET .LE. 0.0) THEN
C                 SET TRANPIRATION TO ZERO AND CALCULATE LEAF TEMP
                  SUMET=0.0
                  DO 31 I=1,NC
                     IF (CANLAI(J,I) .GT. 0.0) THEN
                     HUMID(I)=0.0
                     RSTOM(I)=1.0E20
                     ETCAN(J,I)=0.0
                     TLC(J,I)=AVGTMP(I) + (SWCAN(J,I)+LWCAN(J,I))
     >                              *RHCAN(J,I)/(CANLAI(J,I)*RHOA*CA)
                     INIT(J)=1
                     END IF
   31             CONTINUE
                  GO TO 48
               END IF
C
C              SOLVE MATRIX FOR CHANGE IN PXYLEM(J) 
               DO 32 I=1,NC  
                  IF (CANLAI(J,I) .GT. 0.0) THEN
                  AA1(I)=DF1DP(I)-(DF1DT(I)/DF2DT(I))*DF2DP(I)
                  CC1(I)=DF1DX(I)-(DF1DT(I)/DF2DT(I))*DF2DX(I)
                  F1(I) =  F1(I) -(DF1DT(I)/DF2DT(I))*F1(I)
                  CC3=CC3-(DF3DP(I)/AA1(I))*CC1(I)
                  F3 = F3-(DF3DP(I)/AA1(I))*F1(I)
                  END IF
   32          CONTINUE
               DX=F3/CC3
               PXYLEM(J)=PXYLEM(J)-DX
C
C              SOLVE MATRIX FOR CHANGE IN PLEAF(I) AND TLC(J,I) 
               DO 33 I=1,NC
                  IF (CANLAI(J,I) .GT. 0.0) THEN
                  DP(I)=(F1(I)-CC1(I)*DX)/AA1(I)
                  DTLC(I)=(F2(I)-DF2DX(I)*DX-DF2DP(I)*DP(I))/DF2DT(I)
C                 ADJUST LEAF TEMP & POTENTIAL
                  PLEAF(I)=PLEAF(I)-DP(I)
                  IF (PLEAF(I) .GT. PXYLEM(J))
     >                PLEAF(I)=(PLEAF(I)+DP(I)+PXYLEM(J))/2.
                  TLC(J,I)=TLC(J,I)-DTLC(I)
C                 CHECK IF TEMPERATURE CHANGE IS WITHIN 0.01 C
                  IF (ABS(DTLC(I)) .GT. 0.01) IFLAG=IFLAG+1
                  END IF
   33          CONTINUE
C
C*****         CHECK IF TOLERANCES HAVE BEEN MET
               IF (IFLAG .GT. 0) THEN
                  ITER1 = ITER1 + 1
                  IF (ITER1 .LT. 20) GO TO 26
               END IF
C
C*****         SOLUTION HAS BEEN FOUND FOR LEAF TEMP AND TRANPIRATION
C              FIND FINAL XYLEM POTENTIAL FOR USE IN ROOT EXTRACTION
               IF (SRROOT*SUMET .NE. 0.0) THEN
                  PXYLEM(J)=(RSOIL-SUMET)/SRROOT
                  IF (PXYLEM(J).GT.AVGMAT(MIN) .OR. IROOT.GT.0) THEN
                     RSOIL=0.0
                     SRROOT=0.0
                     DO 34 I=1,NS
                     IF (ROOTDN(J,I).GT. 0.0) THEN
C                       DON'T INCLUDE SOIL LAYERS DRYER THAN PLANT XYLEM
                        IF (AVGMAT(I) .GE. PXYLEM(J)) THEN
                           RSOIL=RSOIL + RROOT(J,I)*AVGMAT(I)
                           SRROOT= SRROOT + RROOT(J,I)
                        END IF
                     END IF
   34                CONTINUE
C                    RECALCULATE WATER POTENTIAL IN XYLEM
                     PXYLEM(J)=(RSOIL-SUMET)/SRROOT
                  END IF
C
C                 CALCULATE ROOT EXTRACTION FROM EACH SOIL LAYER
                  DO 35 I=1,NS
                     IF (AVGMAT(I) .GT. PXYLEM(J)) XTRACT(I) = XTRACT(I)
     >               + FRACTN*TOTROT(J)*RROOT(J,I)*(AVGMAT(I)-PXYLEM(J))
   35             CONTINUE
               ELSE
C                 SOIL TOO DRY FOR PLANTS TO EXTACT WATER
                  SUMET=0.0
               END IF
C
C              SUM TRANSPIRATION FROM ENTIRE TRANSPIRING CANOPY
               TRNSP(NPLANT+1)=TRNSP(NPLANT+1)+FRACTN*SUMET
            END IF
C
      ELSE !枯死的植被
C********   DEAD PLANT MATERIAL -- COMPUTE WATER CONTENT AND TEMPERATURE
            DO 45 I=1,NC
               IF (CANLAI(J,I) .LE. 0.0) GO TO 45
C
               PEVAP(I)=0.0
               RSTOM(I)=0.0
               B1LV = DRYCAN(J,I)/DT
               A1 = -CANLAI(J,I)*RHOA*CA/RHCAN(J,I)
               B1=LV*B1LV
C              CALCULATE WATER CONTENT BASED ON VAPOR DENSITY OF AIR
               CALL VSLOPE (DUMMY,SATV,AVGTMP(I))
               HUMCAN=AVGVAP(I)/SATV
               CALL CANHUM (2,HUMCAN,DUMMY,WCHUM,TCDT(I),CANMA,CANMB)
C
               IF (INIT(J) .EQ. 1) THEN
C                 INITIALIZE VARIABLES FOR THIS TIME STEP
                  ETCAN(J,I)= -B1LV*(WCANDT(I)-WCAN(I))
                  IF (ETCAN(J,I) .EQ. 0.0) THEN
                     TLC(J,I) = AVGTMP(I)
C                    COMPUTE HUMIDITY IN PLANT MATERIAL
                     CALL CANHUM (1,HUMID(I),DUMMY,WCANDT(I),TCDT(I),
     >                            CANMA,CANMB)
                     CALL VSLOPE (DUMMY,SATV,TLC(J,I))
                     ETCAN(J,I)=CANLAI(J,I)*
     >                          (HUMID(I)*SATV-AVGVAP(I))/RHCAN(J,I)
                     WCANDT(I) = WCAN(I) - ETCAN(J,I)/B1LV
C                    CHECK IF WATER CONTENT IS REASONABLE --
                     IF ((WCANDT(I)-WCHUM)*(WCAN(I)-WCHUM).LT.0.0) THEN
C                       WATER CONTENT WENT BEYOND EQUILIBRUIM WITH AIR
                        ETCAN(J,I)= -B1LV*(WCHUM-WCAN(I))
                        HUMID(I)=(AVGVAP(I) +
     >                       ETCAN(J,I)*RHCAN(J,I)/CANLAI(J,I))/SATV
                        CALL CANHUM (2,HUMID(I),DUMMY,WCANDT(I),TCDT(I),
     >                               CANMA,CANMB)
                        ETCAN(J,I)= -B1LV*(WCANDT(I)-WCAN(I))
                     END IF
                  END IF
                  TLC(J,I) = AVGTMP(I) +
     >                       (LV*ETCAN(J,I)-SWCAN(J,I)-LWCAN(J,I))/A1
               END IF
C
               CALL VSLOPE (VTSLOP(I),SATV,TLC(J,I))
C
C*****         BEGIN ITERATIONS TO FIND LEAF TEMP AND WATER CONTENT
               ITER1 = 0
   40          ETCAN(J,I)= -B1LV*(WCANDT(I)-WCAN(I))
C              COMPUTE HUMIDITY IN PLANT MATERIAL AT END OF TIME STEP
               CALL CANHUM(1,HUMID(I),HWSLOP,WCANDT(I),TCDT(I),
     >                     CANMA,CANMB)
C***           SET UP AND SOLVE 2X2 MATRIC FOR NEWTON-RAPHSON APPROX.
C              FOR LEAF TEMP AND WATER CONTENT
               A2 = CANLAI(J,I)*HUMID(I)*VTSLOP(I)/RHCAN(J,I)
               B2 = CANLAI(J,I)*SATV*HWSLOP/RHCAN(J,I) + B1LV
               FF1 = -LV*ETCAN(J,I) + A1*(TLC(J,I)-AVGTMP(I))
     >              + SWCAN(J,I) + LWCAN(J,I)
               FF2=CANLAI(J,I)*(HUMID(I)*SATV-AVGVAP(I))/RHCAN(J,I)
     >                        - ETCAN(J,I)
               DET = A1*B2 - A2*B1
               D1 = (FF1*B2 - FF2*B1)/DET
               D2 = (A1*FF2 - A2*FF1)/DET
C
C***           UPDATE VALUES
               TLC(J,I)=TLC(J,I) - D1
               WCANDT(I) = WCANDT(I) - D2
C              CHECK IF WATER CONTENT IS REASONABLE --
               CALL VSLOPE (VTSLOP(I),SATV,TLC(J,I))
               HUMCAN=AVGVAP(I)/SATV
               CALL CANHUM (2,HUMCAN,DUMMY,WCHUM,TCDT(I),CANMA,CANMB)
               DELMAX=WCAN(I)-WCHUM
               DELTA=WCANDT(I)-WCHUM
CCCCC          IF(DELTA*DELMAX.LT.0.0 .OR.ABS(DELMAX).LT.ABS(DELTA))THEN
               IF(DELTA*DELMAX.LT.0.0)THEN
C                 WATER CONTENT WENT BEYOND EQUILIBRUIM WITH HUMIDITY
                  ETCAN(J,I)= -B1LV*(WCHUM-WCAN(I))
                  CALL CANHUM (1,HUM,DUMMY,WCAN(I),TCDT(I),CANMA,CANMB)
                  ETMAX=CANLAI(J,I)*(HUM*SATV-AVGVAP(I))/RHCAN(J,I)
                  IF (ABS(ETCAN(J,I)) .GT. ABS(ETMAX)) ETCAN(J,I)=ETMAX
                  HUMID(I)=(AVGVAP(I) +
     >                       ETCAN(J,I)*RHCAN(J,I)/CANLAI(J,I))/SATV
                  CALL CANHUM (2,HUMID(I),DUMMY,WCANDT(I),TCDT(I),
     >                         CANMA,CANMB)
               END IF
               IF (ABS(D1) .GT. 0.01) THEN
                  ITER1 = ITER1 +1
                  IF (ITER1 .LT. 10) GO TO 40
               END IF
C              CALCULATE EVAPORATION FROM THIS LAYER
               ETCAN(J,I)= -B1LV*(WCANDT(I)-WCAN(I))
               SUMET=SUMET+ETCAN(J,I)
   45       CONTINUE
            INIT(J) = 2
         END IF
C
C********STORE EVAPORATION/TRANSPIRATION FROM PLANT SPECIES
   48    TRNSP(J)=FRACTN*SUMET+SUMEV
C
C        SUM HEAT AND WATER TRANSFER IN EACH LAYER FROM ALL CANOPY TYPES
         DO 50 I=1,NC
            IF (CANLAI(J,I) .LE. 0.0) GO TO 50
C     叶片向空气的热通量
            HEATC(I)=HEATC(I)
     >             + CANLAI(J,I)*RHOA*CA*(TLC(J,I)-AVGTMP(I))/RHCAN(J,I)
C     ETLYR (I)    evapotranspiration within canopy layer i
            ETLYR(I)=ETLYR(I) + FRACTN*ETCAN(J,I) + PEVAP(I)
C     DHEATC (I)   conductance term for heat flux term from plants leaves within canopy
            DHEATC(I)=DHEATC(I)
     >             + CANLAI(J,I)*RHOA*CA/RHCAN(J,I)
C     average conductance (stomatal plus leaf boundary layer) for transpiration
            DETLYR(I)=DETLYR(I) + CANLAI(J,I)/(RSTOM(I)+RHCAN(J,I))
            DTLDTC(I)=DTLDTC(I) + LV*HUMID(I)*VTSLOP(I)*
     >                            CANLAI(J,I)/(RSTOM(I)+RHCAN(J,I))
   50    CONTINUE
   60 CONTINUE
C
C     CALCULATE CHANGE OF LEAF TEMP WITH RESPECT TO CANOPY TEMP
      DO 70 I=1,NC
         DTLDTC(I)=DHEATC(I)/(DHEATC(I)+DTLDTC(I))
   70 CONTINUE
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE CANHUM (NCALC,HUM,DHDW,WCAN,TC,CANMA,CANMB)
C
C     THIS SUBROUTINE DEFINES THE RELATION BETWEEN THE HUMIDITY AND THE
C     MOISTURE CONTENT OF THE DEAD CANOPY MATERIAL USING THE CONVERSION
C     BETWEEN WATER POTENTIAL AND HUMIDITY.  DEPENDING ON THE VALUE OF
C     NCALC, THE SUBROUTINE WILL EITHER CALCULATE HUMIDITY AND THE
C     DERIVATIVE OF HUMIDTY WITH RESPECT TO WATER CONTENT FOR A GIVEN
C     WATER CONTENT, OR AN EQUILIBRIUM WATER CONTENT GIVEN THE HUMIDITY.
C     (NCALC = 1: CALCULATE HUMIDITY; OTHERWISE CALCULATE WATER CONTENT)
C***********************************************************************
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      REAL LF,LV,LS
C
C
      IF (NCALC .EQ. 1) THEN
C        CALCULATE HUMIDITY AND SLOPE BASED ON WATER CONTENT
         IF (WCAN .LE. 0.0) THEN
C           EQUATIONS CANNOT HANDLE CASE WHERE WCAN=0.0 -- SET HUM=0.0
            HUM=0.0
            DHDW=0.0
           ELSE
C
            HUM=EXP(0.018*G/(UGAS*(TC+273.16))*CANMA*WCAN**(-CANMB))
            DHDW=-0.018*G*CANMA*CANMB*HUM*WCAN**(-CANMB-1.)
     >           /(UGAS*(TC+273.16))
         END IF
        ELSE
C        CALCULATE WATER CONTENT BASED ON HUMIDITY
         IF (HUM .LT. 0.999) THEN
          IF (HUM .LE. 0.0) THEN
            WCAN = 0.0
           ELSE
            WCAN=(ALOG(HUM)*UGAS*(TC+273.16)/0.018/G/CANMA)**(-1./CANMB)
          END IF
         ELSE
          WCAN=(ALOG(0.999)*UGAS*(TC+273.16)/0.018/G/CANMA)**(-1./CANMB)
         END IF
      END IF
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE CANTK (NC,CON,TC,ZC)
C
C     THIS SUBROUTINE COMPUTES THE EDDY CONDUCTANCE COEFFICIENT THROUGH
C     THE CANOPY
C
C***********************************************************************
C
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /WINDV/ ZH,ZM,ZERO,USTAR,STABLE,WINDC(11),WINDR(10)
      REAL LF,LV,LS
      REAL CON(10),TC(11),ZC(11)
C
      RICHRD=STABLE
      IF (RICHRD .GE. 0) THEN
C        STABLE CONDITIONS IN CANOPY
         DIABAT=1
        ELSE
C        UNSTABLE CONDITIONS IN CANOPY
         IF (RICHRD .LT. -2.) RICHRD=-2
         DIABAT=1/SQRT(1-16*RICHRD)
      END IF
C
      DO 10 I=1,NC
         AVGZ=ZC(NC+1)-(ZC(I+1)+ZC(I))/2
         IF (AVGZ .LT. ZERO) THEN
            TKCAN=VONKRM*RHOA*CA*USTAR*ZH/DIABAT
           ELSE
            TKCAN=VONKRM*RHOA*CA*USTAR*(AVGZ+ZH-ZERO)
         END IF
         IF (TKCAN .LT. TKA) TKCAN=TKA
C     CON (I)      conductance term between nodes I and I+1  K/(Z(I+1)-Z(I))
         CON(I) = TKCAN/(ZC(I+1)-ZC(I))
cccc         WRITE (21,*) 'LAYER ',I
cccc         WRITE (21,*) 'WINDC(I),WINDC(I+1):',WINDC(I),WINDC(I+1)
cccc         WRITE (21,*) 'TC(I),TC(I+1),AVGZ: ',TC(I),TC(I+1),AVGZ
cccc         WRITE (21,*) 'RICHRD,FREE,DIABAT: ',RICHRD,FREE,DIABAT
cccc         WRITE (21,*) 'TKCAN,CON(I),USTAR: ',TKCAN,CON(I),USTAR
   10 CONTINUE
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE EBSNOW (N,NSP,NR,ICESPT,TSP,TSPDT,DLW,DLWDT,RHOSP,ZSP,
     >                   DZSP,QVSP,VAPSP,VAPSPT,SSP,ITER)
C
C     THIS SUBROUTINE CALCULATES THE JACOBIAN MATRIX COEFFICIENTS FOR
C     THE SNOW PORTION OF THE NEWTON-RAPHSON SOLUTION OF THE ENERGY
C     BALANCE
C
C***********************************************************************
      COMMON /TIMEWT/ WT,WDT,DT
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /MATRIX/ A1(150),B1(150),C1(150),D1(150),
     >                A2(150),B2(150),C2(150),D2(150)
      REAL TSP(100),TSPDT(100),DLW(100),DLWDT(100),RHOSP(100),
     >     ZSP(100),DZSP(100),QVSP(100),SSP(100)
      REAL QVSPT(100),QVSPDT(100),TK(100),CON(100),CONV(100),SLOPE(100),
     >     CSP(100),CSPT(100),CSPDT(100)
      INTEGER ICESPT(100)
      REAL LF,LS,LV
      SAVE QVSPT,CON,CSPT
C
C**** DETERMINE THE VAPOR FLUX BETWEEN NODES
C     CONV (I)     vapor conductance between node i and i+1

      IF (ITER .EQ. 1) CALL QVSNOW (NSP,QVSPT,CONV,SLOPE,TSP,ZSP,VAPSP)
      CALL QVSNOW (NSP,QVSPDT,CONV,SLOPE,TSPDT,ZSP,VAPSPT)
C     IF UNDERLYING MATERIAL IS RESIDUE, VAPOR DENSITY IS NOT A FUNCTION
C     OF TEMPERATURE, I.E.  SLOPE(NSP+1) = 0.0
      IF (NR .GT. 0) SLOPE(NSP+1) = 0.0
C
C**** OBTAIN THE AVERAGE VAPOR FLUX OVER THE TIME STEP
C     QVSP (I)     vapor flux between nodes i and i+1 of the snowpack (kg/s/m^2)
C     QVSPT (I)    vapor flux between nodes i and i+1 of the snowpack at time t (kg/s/m^2)

      CALL WEIGHT (NSP,QVSP,QVSPT,QVSPDT)
C
      IF (ITER .EQ. 1) THEN
C        (DENSITY AND THERMAL CONDUCTIVITY ARE CONSTANT OVER TIME STEP)
C****    DETERMINE THE THERMAL CONDUCTIVITY OF EACH NODE
         CALL SNOWTK (NSP,TK,RHOSP)
         TK(NSP+1) = TK(NSP)
C
C****    DETERMINE THE CONDUCTANCE TERM BETWEEN NODES
         NSP1 = NSP+1
         CALL CONDUC (NSP1,ZSP,TK,CON)
      END IF
C
C**** CALCULATE THE SPECIFIC HEAT OF EACH NODE
      IF (ITER .EQ. 1) CALL SNOWHT (NSP,CSPT,TSP,RHOSP)
      CALL SNOWHT (NSP,CSPDT,TSPDT,RHOSP)
C
C**** OBTAIN THE AVERAGE SPECIFIC HEAT OVER THE TIME STEP
      CALL WEIGHT (NSP,CSP,CSPT,CSPDT)
C
C
C**** DETERMINE THE COEFFICIENTS FOR THE SURFACE LAYER
      D1(N)= D1(N) - CON(1)*(WT*(TSP(1)-TSP(2))
     >                      + WDT*(TSPDT(1)-TSPDT(2))) + SSP(1)
     >           - CSP(1)*(TSPDT(1)-TSP(1))*DZSP(1)/DT
     >           - RHOL*LF*(DLWDT(1)-DLW(1))/DT - LS*QVSP(1)
C
      IF (ICESPT(1) .EQ. 0) THEN
C        LAYER IS NOT MELTING - ENERGY BUDGET BASED ON TEMPERATURE
         A1(N+1) = WDT*(CON(1) + CONV(1)*LS*SLOPE(1))
         B1(N) = B1(N) - WDT*(CON(1) + CONV(1)*LS*SLOPE(1))
     >                 - DZSP(1)*CSP(1)/DT
        ELSE
C
C        LAYER IS MELTING - ENERGY BUDGET BASED ON WATER CONTENT
         A1(N+1)= 0.0
         B1(N)= B1(N) - RHOL*LF/DT
C        IF SNOW IS NOT FIRST MATERIAL, SET DERIV. FOR LAST NODE TO 0.0
         IF (N.GT.1) C1(N-1)=0.0
      END IF
C
C**** DETERMINE THE COEFFICIENTS FOR THE REMAINDER OF THE SNOWPACK
      DO 20 I=N+1,N+NSP-1
         J=I-N+1
         D1(I)= CON(J-1)*(WT*(TSP(J-1)-TSP(J))
     >                   + WDT*(TSPDT(J-1)-TSPDT(J)))
     >         - CON(J)*(WT*(TSP(J)-TSP(J+1))
     >                   + WDT*(TSPDT(J)-TSPDT(J+1))) + SSP(J)
     >         - CSP(J)*(TSPDT(J)-TSP(J))*DZSP(J)/DT
     >         - RHOL*LF*(DLWDT(J)-DLW(J))/DT - LS*(QVSP(J)-QVSP(J-1))
C
         IF (ICESPT(J) .EQ. 0) THEN
C           LAYER IS NOT MELTING - ENERGY BUDGET BASED ON TEMPERATURE
            A1(I+1)= WDT*(CON(J) + CONV(J)*LS*SLOPE(J))
            B1(I)= -WDT*(CON(J-1)+CON(J))
     >             - WDT*LS*SLOPE(J)*(CONV(J-1) + CONV(J))
     >             - DZSP(J)*CSP(J)/DT
            C1(I-1)= WDT*(CON(J-1) + CONV(J-1)*LS*SLOPE(J))
           ELSE
C
C           LAYER IS MELTING - ENERGY BUDGET BASED ON WATER CONTENT
            A1(I+1)= 0.0
            B1(I)= -RHOL*LF/DT
            C1(I-1)= 0.0
         END IF
20    CONTINUE
C
C
C**** DETERMINE THE BOUNDARY CONDITIONS FOR TOP LAYER OF NEXT MATERIAL
      N=N+NSP
C
      IF (NR .GT. 0) THEN
C        SNOW OVERLYING RESIDUE
C        CHECK IF LAST SNOW NODE IS MELTING - IF SO, ENERGY BALANCE
C        IS BASED ON WATER CONTENT, NOT TEMPERATURE AND A1(N)=0.0
         IF (ICESPT(NSP) .EQ. 0) A1(N) =  WDT*CON(NSP)
         B1(N) = -WDT*CON(NSP)
         C1(N-1) = C1(N-1) + WDT*CON(NSP)
         D1(N) = CON(NSP)* (WT*(TSP(NSP) - TSP(NSP+1))
     >                     + WDT*(TSPDT(NSP) - TSPDT(NSP+1)))
C
        ELSE
C        SNOW IS LYING ON BARE SOIL - INCLUDE LATENT HEAT TRANSFER
C        AND VAPOR FLUX DEPENDENCE ON TEMPERATURE OF SOIL SURFACE
C        CHECK IF LAST SNOW NODE IS MELTING - IF SO, ENERGY BALANCE
C        IS BASED ON WATER CONTENT, NOT TEMPERATURE AND A1(N)=0.0
         IF (ICESPT(NSP) .GT. 0)
     >      A1(N)=WDT*(CON(NSP)+CONV(NSP)*LV*SLOPE(NSP))
         B1(N) =-WDT*(CON(NSP) + CONV(NSP)*LV*SLOPE(NSP+1))
         C1(N-1) = WDT*(CON(NSP) + CONV(NSP)*LS*SLOPE(NSP+1))
         D1(N) = CON(NSP)* (WT*(TSP(NSP) - TSP(NSP+1))
     >                 + WDT*(TSPDT(NSP) - TSPDT(NSP+1))) + LV*QVSP(NSP)
      END IF
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE QVSNOW (NSP,QVSP,CONV,SLOPE,TSP,ZSP,VAPSP)
C
C     THIS SUBROUTINE CALCULATES THE VAPOR DIFFUSION IN THE SNOWPACK
C
C***********************************************************************
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /SPPROP/ G1,G2,G3,EXTSP,TKSPA,TKSPB,TKSPEX,VDIFSP,VAPSPX
      REAL QVSP(100),TSP(100),ZSP(100),VAPSNO(100),
     >     VAPICE(100),SLOPE(100),CONV(100)
      REAL LF,LS,LV
C
C**** DETERMINE THE VAPOR DIFFUSIVITY, DENSITY AND SLOPE OF EACH NODE
      DO 10 I=1,NSP
         VAPSNO(I) = VDIFSP*(P0/PRESUR)*(1.+TSP(I)/273.16)**VAPSPX
         CALL VSLOPE (SLOPE(I),SATV,TSP(I))
C        CALCULATE THE SATURATED VAPOR DENSITY OVER ICE
         TMP = TSP(I) + 273.16
         HUMID = EXP(0.018/(UGAS*TMP)*LF*TSP(I)/TMP)
         VAPICE(I) = HUMID*SATV
         SLOPE(I) = HUMID*SLOPE(I)
   10 CONTINUE
      VAPICE(NSP+1) = VAPSP
      VAPSNO(NSP+1) = VAPSNO(NSP)
      CALL VSLOPE (SLOPE(NSP+1),SATV,TSP(NSP+1))
      SLOPE(NSP+1) = SLOPE(NSP+1)*VAPSP/SATV
C
C**** DETERMINE THE CONDUCTANCE TERM BETWEEN NODES
      NSP1 = NSP+1
      CALL CONDUC (NSP1,ZSP,VAPSNO,CONV)
C
C**** CALCULATE THE VAPOR FLUX BETWEEN EACH NODE
      DO 20 I=1,NSP
         QVSP(I) = CONV(I)*(VAPICE(I)-VAPICE(I+1))
   20 CONTINUE
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE SNOWTK (NSP,TK,RHOSP)
C
C     THIS SUBROUTINE CALCULATES THE THERMAL CONDUCTIVITY OF THE
C     SNOWPACK LAYERS.  THE EQUATION USED IS OF THE FORM:
C
C           K = A + B*(RHOSP/RHOL)**C
C
C     WHERE:    A = TKSPA      B = TKSPB     C = TKSPEX
C
C***********************************************************************
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /SPPROP/ G1,G2,G3,EXTSP,TKSPA,TKSPB,TKSPEX,VDIFSP,VAPSPX
      REAL TK(100),RHOSP(100)
      REAL LF,LS,LV
C
      DO 10 I=1,NSP
         TK(I) = TKSPA + TKSPB*(RHOSP(I)/RHOL)**TKSPEX
   10 CONTINUE
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE SNOWHT (NSP,CSP,TSP,RHOSP)
C
C     THIS SUBOUTINE CALCULATES THE VOLUMETRIC SPECIFIC HEAT FOR EACH
C     RESIDUE NODE, INCLUDING LATENT HEAT EFFECTS.
C
C***********************************************************************
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      REAL CSP(100),TSP(100),RHOSP(100)
      REAL LF,LS,LV
C
      DO 10 I=1,NSP
C        CALCULATE THE SPECIFIC HEAT OF THE SNOW
         SPHEAT = 92.96 + 7.37*(TSP(I)+273.16)
C        INCLUDE THE LATENT HEAT TERM
         CALL VSLOPE (S,DUMMY,TSP(I))
C        CALCULATE THE VOLUMETRIC SPECIFIC HEAT INCLUDING LATENT HEAT
         CSP(I) = RHOSP(I)*SPHEAT + (1.-RHOSP(I)/RHOI)*LS*S
   10 CONTINUE
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE EBRES (N,NR,NSP,ZR,TR,TRDT,VAPR,VAPRDT,GMC,GMCDT,
     >                  GMCMAX,RHOR,RESCOF,SR,UR,QVR,RHOSP,ITER)
C
C     THIS SUBOUTINE CALCULATES THE NEWTON-RAPHSON COEFFICIENTS FOR THE
C     ENERGY BALANCE OF THE RESIDUE LAYERS.
C
C***********************************************************************
      COMMON /TIMEWT/ WT,WDT,DT
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /MATRIX/ A1(150),B1(150),C1(150),D1(150),
     >                A2(150),B2(150),C2(150),D2(150)
      COMMON /RESIDU/ EVAP(10),EVAPK(10)
      REAL LF,LS,LV
C
      REAL ZR(10),TR(10),TRDT(10),VAPR(10),VAPRDT(10),GMC(10),GMCDT(10),
     >     RHOR(10),SR(10),UR(10),QVR(10),RHOSP(100)
      REAL CRES(10),CREST(10),CRESDT(10),TKRES(10),CONVEC(10),
     >     VAPCON(10),CON(10),CONV(10)
      SAVE CREST
C
C**** DETERMINE THE EVAPORATION FROM RESIDUE ELEMENTS AT EACH NODE
      CALL RESVAP (NR,EVAP,EVAPK,ZR,TR,TRDT,VAPR,VAPRDT,GMC,GMCDT,
     >             GMCMAX,RHOR,RESCOF,UR)
C
C**** DETERMINE THE HEAT TRANSFER COEFFICIENT FOR EACH NODE
      CALL RESTK (NR,NSP,TKRES,CONVEC,TR,TRDT,GMC,GMCDT,RHOR,RHOSP)
      NR1=NR+1
      TKRES(NR1)=TKRES(NR)
C
C**** DETERMINE THE CONDUCTANCE TERM FOR HEAT TRANSPORT BETWEEN NODES
      CALL CONDUC (NR1,ZR,TKRES,CON)
C
C**** DETERMINE THE VAPOR TRANSPORT FROM THE THERMAL CONVECTION
      CALL RESVK (NR,NSP,TR,TRDT,CONVEC,VAPCON)
      VAPCON(NR1)=VAPCON(NR)
C
C**** DETERMINE THE CONDUCTANCE TERM FOR CONVECTIVE VAPOR TRANSPORT
      CALL CONDUC (NR1,ZR,VAPCON,CONV)
C
C**** DETERMINE THE VAPOR FLUX BETWEEN RESIDUE NODES
      CALL QVRES (NR,QVR,CONV,VAPR,VAPRDT)
C
C**** DETERMINE THE VOLUMETRIC HEAT CAPACITY OF EACH NODE
      IF (ITER .EQ. 1) CALL RESHT (NR,CREST,GMC,RHOR)
      CALL RESHT (NR,CRESDT,GMCDT,RHOR)
C
C**** CALCULATE THE AVERAGE HEAT CAPACITY OVER THE TIME STEP
      CALL WEIGHT (NR,CRES,CREST,CRESDT)
C
C
C**** DETERMINE THE MATRIX COEFFICIENTS FOR THE TOP LAYER
C
      A1(N+1)=A1(N+1) + WDT*CON(1)
      C1(N)=C1(N) + WDT*CON(1)
      IF (NR.GT.1) THEN
         B1(N)=B1(N) - WDT*CON(1) - (ZR(2)-ZR(1))/(2*DT)*CRES(1)
         D1(N)=D1(N) - CON(1)*(WT*(TR(1)-TR(2))+WDT*(TRDT(1)-TRDT(2)))
     >            + LV*EVAP(1) + SR(1)
     >            - (ZR(2)-ZR(1))/(2*DT)*CRES(1)*(TRDT(1)-TR(1))
       ELSE
         B1(N)=B1(N) - WDT*CON(1) - (ZR(2)-ZR(1))/DT*CRES(1)
         D1(N)=D1(N) - CON(1)*(WT*(TR(1)-TR(2))+WDT*(TRDT(1)-TRDT(2)))
     >            + LV*EVAP(1) + SR(1)
     >            - (ZR(2)-ZR(1))/DT*CRES(1)*(TRDT(1)-TR(1))
      END IF
C
C**** DETERMINE THE COEFFICIENTS FOR THE REST OF THE LAYERS
      DO 10 I=N+1,N+NR-1
         J=I-N+1
         A1(I+1)=A1(I+1) + WDT*CON(J)
         C1(I)=C1(I) + WDT*CON(J)
         IF (J .NE. NR) THEN
            B1(I)=B1(I) - WDT*(CON(J-1)+CON(J))
     >                  - (ZR(J+1)-ZR(J-1))/(2*DT)*CRES(J)
            D1(I)=CON(J-1)*(WDT*(TRDT(J-1)-TRDT(J))+WT*(TR(J-1)-TR(J)))
     >        -CON(J)*(WDT*(TRDT(J)-TRDT(J+1)) + WT*(TR(J)-TR(J+1)))
     >        +LV*EVAP(J) + SR(J)
     >        -(ZR(J+1)-ZR(J-1))/(2*DT)*CRES(J)*(TRDT(J)-TR(J))
        ELSE
            B1(I)=B1(I) - WDT*(CON(J-1)+CON(J))
     >                  - (ZR(J+1)-ZR(J)+(ZR(J)-ZR(J-1))/2)/DT*CRES(J)
            D1(I)=CON(J-1)*(WDT*(TRDT(J-1)-TRDT(J))+WT*(TR(J-1)-TR(J)))
     >        -CON(J)*(WDT*(TRDT(J)-TRDT(J+1)) + WT*(TR(J)-TR(J+1)))
     >        +LV*EVAP(J) + SR(J)-(ZR(J+1)-ZR(J)+(ZR(J)-ZR(J-1))/2)
     >        /DT*CRES(J)*(TRDT(J)-TR(J))
         END IF
   10 CONTINUE
C
C**** DETERMINE THE COEFFICIENTS FOR THE SOIL SURFACE
C
      N=N+NR
      B1(N)=B1(N) - WDT*CON(NR)
      D1(N)=CON(NR)*(WDT*(TRDT(NR)-TRDT(NR+1))+WT*(TR(NR)-TR(NR+1)))
     >      + LV*QVR(NR) 
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE RESVAP (NR,EVAP,EVAPK,ZR,TR,TRDT,VAPR,VAPRDT,GMC,GMCDT,
     >                   GMCMAX,RHOR,RESCOF,UR)
C
C     THIS SUBROUTINE IS TO DETERMINE THE EVAPORATION WITHIN THE
C     RESIDUE LAYERS ASSUMING THE AVERAGE MOISTURE CONTENT OF THE
C     RESIDUE IS THE MOISTURE CONTENT AT THE BEGINNING OF THE TIME STEP
C
C***********************************************************************
      COMMON /TIMEWT/ WT,WDT,DT
C
      REAL EVAP(10),EVAPK(10),ZR(10),TR(10),TRDT(10),VAPR(10),
     >     VAPRDT(10),GMC(10),GMCDT(10),RHOR(10),UR(10)
C
      DO 20 I=1,NR
C
C****    DETERMINE THE RELATIVE HUMIDITY IN THE RESIDUE
         CALL RESHUM (1,HUMT,DUMMY,GMC(I),TR(I))
         CALL RESHUM (1,HUMDT,DHDW,GMCDT(I),TRDT(I))
C
C****    DETERMINE THE SATURATED VAPOR DENSITY OF THE RESIDUE
         CALL VSLOPE (DUMMY,SATV,TR(I))
         CALL VSLOPE (DUMMY,SATVDT,TRDT(I))
C
         IF (I .EQ. 1) THEN
            DZ=ZR(2)-ZR(1)
            IF (NR .NE. 1) DZ=DZ/2.
           ELSE
            IF (I .EQ. NR) THEN
               DZ=ZR(I+1)-ZR(I)+(ZR(I)-ZR(I-1))/2.
              ELSE
               DZ=(ZR(I+1)-ZR(I-1))/2.
            END IF
         END IF
C
C****    DEFINE EVAPK(I) -- FOR NOW, IT WILL SIMPLY BE SET EQUAL TO
C        THE INVERSE OF RESCOF, ASSUMING THAT RESCOF IS THE VAPOR
C        RESISTANCE.  LATER WORK MAY REQUIRE THAT THE VAPOR TRANSPORT
C        RESISTANCE BE CALCULATED OR OBTAINED FROM A SUBROUTINE.
         EVAPK(I)=1./RESCOF
C
         EVAP(I)=EVAPK(I)*(WDT*(VAPRDT(I)-HUMDT*SATV)
     >                     +WT*(VAPR(I)-HUMT*SATV))
         GMCDT(I)= GMC(I) + EVAP(I)*DT/DZ/RHOR(I) + UR(I)
         IF (GMCDT(I) .LT. 0.0) THEN
            GMCDT(I)=0.0
            EVAP(I)= DZ*RHOR(I)*(GMCDT(I)-GMC(I)-UR(I))/DT
         END IF
         IF (GMCDT(I) .GT. GMCMAX) THEN
            GMCDT(I)=GMCMAX
            EVAP(I)= DZ*RHOR(I)*(GMCDT(I)-GMC(I)-UR(I))/DT
         END IF
C
   20 CONTINUE
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE RESHUM (NCALC,HUM,DHDW,GMC,TR)
C
C     THIS SUBROUTINE DEFINES THE RELATION BETWEEN THE HUMIDITY AND THE
C     MOISTURE CONTENT OF THE RESIDUE USING THE CONVERSION
C     BETWEEN WATER POTENTIAL AND HUMIDITY.  DEPENDING ON THE VALUE OF
C     NCALC, THE SUBROUTINE WILL EITHER CALCULATE HUMIDITY AND THE
C     DERIVATIVE OF HUMIDTY WITH RESPECT TO WATER CONTENT FOR A GIVEN
C     WATER CONTENT, OR AN EQUILIBRIUM WATER CONTENT GIVEN THE HUMIDITY.
C     (NCALC = 1: CALCULATE HUMIDITY; OTHERWISE CALCULATE WATER CONTENT)
C***********************************************************************
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /RSPARM/ RESMA,RESMB,RESMC,RESTKA,RESTKB
      REAL LF,LS,LV
C
C**** THIS SUBROUTINE USES FOLLOWING RELATION:
C                TENSION = A*W**B
C                H=EXP(M*TENSION/(UGAS*TEMP))
C     THEREFORE  H=EXP(M*A*W**B/(UGAS*TEMP))
C
      IF (NCALC .EQ. 1) THEN
C        CALCULATE HUMIDITY AND SLOPE BASED ON WATER CONTENT
         IF (GMC .LE. 0.0) THEN
C           EQUATIONS CANNOT HANDLE CASE WHERE GMC=0.0 -- SET HUM=0.0
            HUM=0.0
            DHDW=0.0
           ELSE
            HUM=EXP(0.018*G/(UGAS*(TR+273.16))*RESMA*GMC**(-RESMB))
            DHDW=-0.018*G*RESMA*RESMB*HUM*GMC**(-RESMB-1.)
     >           /(UGAS*(TR+273.16))
         END IF
        ELSE
C        CALCULATE WATER CONTENT BASED ON HUMIDITY
         IF (HUM .GT. 0.999) THEN
           GMC=(ALOG(0.999)*UGAS*(TR+273.16)/0.018/G/RESMA)**(-1./RESMB)
          ELSE
           GMC=(ALOG(HUM)*UGAS*(TR+273.16)/0.018/G/RESMA)**(-1./RESMB)
         END IF
      END IF
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE RESTK(NR,NSP,TKRES,CONVEC,TR,TRDT,GMC,GMCDT,RHOR,RHOSP)
C
C     THIS SUBROUTINE CALCULATES THE THERMAL CONDUCTANCE TERM BETWEEN
C     RESIDUE NODES USING THE CALCULATED WINDSPEED AT A NODE FOR THE
C     THE CONVECTIVE TRANSPORT TERM, AND THE THERMAL CONDUCTIVITIES OF
C     OF THE RESIDUE AND WATER THE THERMAL CONDUCTIVITY TERM.
C     CONVECTIVE AND CONDUCTIVE TERMS ARE THEN WEIGHTED ACCORDING TO
C     THE VOLUMETRIC FRACTIONS OF AIR, WATER AND RESIDUE.
C
C***********************************************************************
      COMMON /TIMEWT/ WT,WDT,DT
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /RSPARM/ RESMA,RESMB,RESMC,RESTKA,RESTKB
      COMMON /WINDV/ ZH,ZM,ZERO,USTAR,STABLE,WINDC(11),WINDR(10)
      REAL LF,LS,LV
C
      REAL TR(10),TRDT(10),GMC(10),GMCDT(10),RHOR(10),RHOSP(100)
      REAL CONVEC(10),TKRES(10),TKSP(100)
      DATA RESDEN/ 170./ 
C
      IF (NSP .GT. 0) THEN
C        ASSUME THE SNOW FILTERS DOWN INTO THE RESIDUE - THEREFORE
C        REDEFINE THE THERMAL CONVECTION AS A CONDUCTION THRU SNOW
         CALL SNOWTK (NSP,TKSP,RHOSP)
      END IF
C
      DO 10 I=1,NR
C
C        AVERAGE THE TEMP AND WATER CONTENT OVER TIME
         AVGTMP = WDT*TRDT(I) + WT*TR(I)
         AVGGMC = WDT*GMCDT(I) + WT*GMC(I)
C
C        CALCULATE THE VOLUME FRACTION OF EACH MATERIAL
         RESVOL = RHOR(I)/RESDEN
         WATER = AVGGMC*RHOR(I)/RHOL
         AIR = 1. - RESVOL - WATER
C
C        IF SNOW IS PRESENT - DO NOT CALCULATE CONVECTION 
         IF (NSP .GT. 0) THEN
            CONVEC(I) = TKSP(NSP)
           ELSE
C
C           CALCULATE THE CONVECTIVE HEAT TRANSFER COEFFICIENT
            CONVEC(I) = TKA*(1.+RESTKA*AVGTMP)*(1.+RESTKB*WINDR(I))
         END IF
C
C        ADJUST CONVECTIVE TRANSPORT FOR AIR POROSITY
C        (IN THE CASE OF SNOW, AIR POROSITY IS THE FRACTION OF SNOW)
         CONVEC(I) = CONVEC(I)*AIR
C
C        CALCULATE THERMAL CONDUCTIVITY THROUGH THE RESIDUE AND WATER
         TK =  WATER*TKL + RESVOL*TKR
C
C        CALCULATE THE EFFECTIVE HEAT TRANSFER COEFFICIENT
         TKRES(I) = CONVEC(I) + TK
   10 CONTINUE
      RETURN
      END 
C***********************************************************************
C
      SUBROUTINE RESVK (NR,NSP,TR,TRDT,CONVEC,VAPCON)
C
C     THIS SUBROUTINE CALCULATES THE VAPOR CONDUCTANCE TERM (CONV)
C     BETWEEN RESIDUE NODES USING THE CALCULATED THERMAL CONDUCTANCE.
C
C***********************************************************************
      COMMON /TIMEWT/ WT,WDT,DT
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /SPPROP/ G1,G2,G3,EXTSP,TKSPA,TKSPB,TKSPEX,VDIFSP,VAPSPX
      REAL CONVEC(10),VAPCON(10),TR(10),TRDT(10)
      REAL LF,LS,LV
C
      DO 10 I=1,NR
         VAPCON(I)=CONVEC(I)/RHOA/CA
C
         IF (NSP .GT. 0) THEN
C           SNOW FILTERS DOWN THRU RESIDUE - VAPOR DIFFUSIVITY IS THAT
C           THROUGH SNOW
            AVGTMP = WT*TR(I) + WDT*TRDT(I)
            VAPCON(I) = VDIFSP*(P0/PRESUR)*(1.+AVGTMP/273.16)**VAPSPX
         END IF
   10 CONTINUE
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE QVRES (NR,QVR,CONV,VAPR,VAPRDT)
C
C     THIS SUBOUTINE CALCULATES VAPOR FLUX BETWEEN RESIDUE NODES.
C
C***********************************************************************
      COMMON /TIMEWT/ WT,WDT,DT
      REAL QVR(10),CONV(10),VAPR(10),VAPRDT(10)
C
      DO 10 I=1,NR
         QVR(I) = CONV(I)*(WT*(VAPR(I)-VAPR(I+1))
     >                 + WDT*(VAPRDT(I)-VAPRDT(I+1)))
   10 CONTINUE
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE RESHT (NR,CRES,GMC,RHOR)
C
C     THIS SUBOUTINE CALCULATES THE VOLUMETRIC SPECIFIC HEAT FOR EACH
C     RESIDUE NODE.
C
C***********************************************************************
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      REAL LF,LS,LV
C
      REAL CRES(10),GMC(10),RHOR(10)
C
      DO 10 I=1,NR
         CRES(I)=RHOR(I)*(CR + GMC(I)*CL)
   10 CONTINUE
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE QVSOIL (NS,QSV,ZS,TS,MAT,VLC,VIC,CONC,DCONVP)
C
C     THIS SUBROUTINE CALCULATES THE VAPOR FLUX BETWEEN SOIL NODES DUE
C     TO GRADIENTS IN POTENTIAL AND TEMPERATURE GRADIENTS. (FLUX DUE TO
C     TEMPERATURE GRADIENTS ARE MULTIPLIED BY AN ENHANCEMENT FACTOR.
C
C***********************************************************************
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /SLPARM/ B(50),ENTRY(50),RHOB(50),SAT(50),SATK(50),
     >                NSALT,SALTKQ(10,50),VAPCOF(50),VAPEXP(50),
     >                SAND(50),SILT(50),CLAY(50),OM(50)
      REAL QSV(50),ZS(50),TS(50),MAT(50),VLC(50),VIC(50),CONC(10,50)
      REAL HUMID(50),DVT(50),DVP(50),CONVT(50),CONVP(50),DCONVP(50)
      REAL LF,LS,LV
C
      DO 20 I=1,NS
         VAC=1. - RHOB(I)*((1.-OM(I))/RHOM+OM(I)/RHOOM) - VLC(I)- VIC(I)
         IF (MAT(I).LT.ENTRY(I) .AND. VAC.GT.0.0) THEN
            DV=VDIFF*(((TS(I)+273.16)/273.16)**2)*(P0/PRESUR)
            DV=DV*VAPCOF(I)*(VAC**VAPEXP(I))
            CALL ENHANC (I,EN,VLC)
            CALL VSLOPE (S,SATVAP,TS(I))
C
C****       DETERMINE THE HUMIDITY FROM THE TOTAL WATER POTENTIAL
            TLCONC=0.0
            DO 10 J=1,NSALT
               TLCONC= TLCONC+CONC(J,I)
   10       CONTINUE
            TOTPOT=MAT(I) - TLCONC*UGAS*(TS(I)+273.16)/G
            HUMID(I)=EXP(.018*G/UGAS/(TS(I)+273.16)*TOTPOT)
            DVT(I)=DV*EN*HUMID(I)*S
            DVP(I)=DV*SATVAP
          ELSE
C           NO AIR POROSITY OR SATURATED --> NO VAPOR DIFFUSION
            DVT(I)=0.0
            DVP(I)=0.0
            HUMID(I)=1.0
         END IF
C
   20 CONTINUE
C
      CALL CONDUC (NS,ZS,DVT,CONVT)
      CALL CONDUC (NS,ZS,DVP,CONVP)
C
      DO 30 I=1,NS-1
         QSV(I)=CONVT(I)*(TS(I)-TS(I+1))+ CONVP(I)*(HUMID(I)-HUMID(I+1))
         DCONVP(I)=.018*G/UGAS/(TS(I)+273.16)*CONVP(I)
   30 CONTINUE
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE ENHANC (I,EN,VLC)
C
C     THIS SUBROUTINE CALCULATES THE ENHANCEMENT COEFFICIENT FOR VAPOR
C     TRANSPORT DUE TO THERMAL GRADIENTS FOR NODE I
C
C***********************************************************************
      COMMON /SLPARM/ B(50),ENTRY(50),RHOB(50),SAT(50),SATK(50),
     >                NSALT,SALTKQ(10,50),VAPCOF(50),VAPEXP(50),
     >                SAND(50),SILT(50),CLAY(50),OM(50)
C
      REAL VLC(50)
C
      DATA  EN1/9.5/ EN2/3.0/ EN4/1.0/ EN5/4.0/
C
      IF (CLAY(I) .LE. 0.02) THEN
C        EN3 BECOMES INFINITELY LARGE WHEN CLAY IS ZERO
C        (SMALLEST CLAY CONTENT IN DATA FROM CASS ET. AL. WAS 2%)
         EN3=SAT(I)*(1+2.6/SQRT(0.02))
        ELSE
         EN3=SAT(I)*(1+2.6/SQRT(CLAY(I)))
      END IF
      EXPON=-(EN3*VLC(I)/SAT(I))**EN5
C
C     EXP(-50.) APPROXIMATELY = 0.0, BUT UNDERFLOW MAY RESULT IF YOU TRY
C     TO USE LARGE NEGATIVE NUMBERS  --->  CHECK IF  < -50.
      IF (EXPON .LE. -50.) THEN
         EXPON=0.0
       ELSE
         EXPON=EXP(EXPON)
      END IF
C
C**** CALCULATE ENHANCEMENT FACTOR
      EN = EN1 + EN2*VLC(I)/SAT(I) - (EN1-EN4)*EXPON
      CONTINUE
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE EBSOIL (N,NS,ZS,TS,TSDT,MAT,MATDT,CONC,CONCDT,
     >                  VLC,VLCDT,VIC,VICDT,ICES,ICESDT,QSL,QSV,SS,ITER)
C
C     THIS SUBROUTINE CALCULATES THE COEFFICIENTS FOR THE SOIL IN THE
C     NEWTON-RAPHSON PROCEDURE OF THE ENERGY BALANCE
C
C***********************************************************************
      COMMON /TIMEWT/ WT,WDT,DT
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /SLPARM/ B(50),ENTRY(50),RHOB(50),SAT(50),SATK(50),
     >                NSALT,SALTKQ(10,50),VAPCOF(50),VAPEXP(50),
     >                SAND(50),SILT(50),CLAY(50),OM(50)
      COMMON /MATRIX/ A1(150),B1(150),C1(150),D1(150),
     >                A2(150),B2(150),C2(150),D2(150)
      COMMON /SPHEAT/ CS(50)
      REAL ZS(50),TS(50),TSDT(50),MAT(50),MATDT(50),CONC(10,50),
     >     CONCDT(10,50),VLC(50),VLCDT(50),VIC(50),VICDT(50),QSLT(50),
     >     QSLDT(50),QSL(50),QSV(50),SS(50),TK(50),TKT(50),TKDT(50),
     >     CON(50),CST(50),CSDT(50),HKT(50),HKDT(50),CONH(50),CONHT(50),
     >     QSVT(50),QSVDT(50),DCONVP(50),DELTA(50),EPSLON(50)
      INTEGER ICES(50),ICESDT(50)
      REAL LF,LV,LS
      SAVE QSVT,QSLT,TKT,CST
C
C
C**** DETERMINE THE AVERAGE VAPOR FLUX BETWEEN SOIL NODES
      IF (ITER.EQ.1) CALL QVSOIL (NS,QSVT,ZS,TS,MAT,VLC,VIC,CONC,DCONVP)
      CALL QVSOIL (NS,QSVDT,ZS,TSDT,MATDT,VLCDT,VICDT,CONCDT,DCONVP)
      CALL WEIGHT (NS-1,QSV,QSVT,QSVDT)
C
C**** DETERMINE THE LIQUID MOISTURE FLUX BETWEEN NODES
      IF (ITER .EQ. 1) THEN
C****    CALCULATE FLUX AT BEGINNING OF TIME STEP
C****    DETERMINE THE HYDRAULIC CONDUCTIVITY OF EACH SOIL NODE
         CALL SOILHK (NS,HKT,MAT,VLC,VIC)
C****    DETERMINE THE CONDUCTANCE TERM BETWEEN SOIL NODES
         CALL CONDUC (NS,ZS,HKT,CONHT)
         CALL QLSOIL (NS,ZS,QSLT,CONHT,MAT)
      END IF
      CALL SOILHK (NS,HKDT,MATDT,VLCDT,VICDT)
      CALL CONDUC (NS,ZS,HKDT,CONH)
      CALL QLSOIL (NS,ZS,QSLDT,CONH,MATDT)
C
C**** OBTAIN THE AVERAGE MOISTURE FLUX OVER THE TIME STEP
      CALL WEIGHT (NS-1,QSL,QSLT,QSLDT)
C
C     LIMIT NET WATER FLUX INTO NODES IF FLUX FILLS NODE BEYOND 
C     AVAILABLE POROSITY AS A RESULT DOWNWARD FLUX INTO LAYERS WITH
C     LIMITED PERMEABILITY (PERHAPS DUE TO ICE)
C     (START AT BOTTOM OF PROFILE AND WORK UP SO THAT ANY CHANGES IN THE
C     FLUXES CAN WORK THEIR UP THROUGH THE PROFILE)     
      IFLAG=0
      DO 5 I=NS-1,2,-1
         IF (VICDT(I).GT.0.0) THEN
            QMAX=(SAT(I) - VLC(I) - VIC(I))*(ZS(I+1) - ZS(I-1))/2./DT
            IF (QMAX.LT.0.0) QMAX=0.0
            IF (QSL(I-1)-QSL(I) .GT. QMAX) THEN
               IFLAG=1
               IF (QSL(I-1).GT.0.0) THEN
                  CONH(I-1)=0.0
                  QSL(I-1)=QSL(I)+QMAX
                  IF (QSL(I-1).LT.0.0) THEN 
                     QSL(I-1)=0.0
                     QSL(I)=QSL(I-1)-QMAX
                     IF (QSL(I).GT.0.0) QSL(I)=0.0
                  END IF
                 ELSE
                  CONH(I)=0.0
                  QSL(I)=QSL(I-1)-QMAX
                  IF (QSL(I).GT.0.0) QSL(I)=0.0
               END IF
            END IF
         END IF
C        SET UP LIQUID FLUX COEFFICIENTS SO THAT HEAT IS CARRIED IN
C        DIRECTION OF MOISTURE MOVEMENT.
         IF (QSL(I) .GT. 0.0) THEN
            DELTA(I)=1.0
            EPSLON(I)=0.0
          ELSE
            DELTA(I)=0.0
            EPSLON(I)=1.0
         END IF
    5 CONTINUE
C
      IF (QSL(1) .GT. 0.0) THEN
         DELTA(1)=1.0
         EPSLON(1)=0.0
       ELSE
         DELTA(1)=0.0
         EPSLON(1)=-1.0
      END IF
C
      IF (VICDT(1).GT.0.0) THEN
         QMAX=(SAT(I) - VLC(1) - VIC(1))*(ZS(2) - ZS(1))/2./DT
         IF (QMAX.LT.0.0) QMAX=0.0
         IF (-QSL(1) .GT. QMAX) THEN
            IFLAG=1
            CONH(1)=0.0
            QSL(1)=-QMAX
         END IF
      END IF
C
      IF (IFLAG.GT.0) THEN
      DO 10 I=2,NS
         IF (VICDT(I).GT.0.0) THEN
            QMAX=(SAT(I) - VLC(I) - VIC(I))*(ZS(I+1) - ZS(I-1))/2./DT
            IF (QMAX.LT.0.0) QMAX=0.0
            IF (QSL(I-1)-QSL(I) .GT. QMAX) THEN
               IF (QSL(I).LT.0.0) THEN
                  CONH(I)=0.0
                  QSL(I)=QSL(I-1)-QMAX
                  IF (QSL(I).GT.0.0) THEN 
                     QSL(I)=0.0
                     QSL(I-1)=QSL(I)+QMAX
                     IF (QSL(I-1).LT.0.0) QSL(I-1)=0.0
                  END IF
                 ELSE
                  CONH(I-1)=0.0
                  QSL(I-1)=QSL(I)+QMAX
                  IF (QSL(I-1).LT.0.0) QSL(I-1)=0.0
               END IF
            END IF
         END IF
   10 CONTINUE
      END IF
C
C**** CALCULATE THERMAL CONDUCTIVITY
      IF (ITER .EQ. 1) CALL SOILTK (NS,TKT,VLC,VIC)
      CALL SOILTK (NS,TKDT,VLCDT,VICDT)
C
C**** OBTAIN AVERAGE CONDUCTIVITY OVER THE TIME STEP
      CALL WEIGHT (NS,TK,TKT,TKDT)
C
C**** CALCULATE CONDUCTANCE TERMS BETWEEN NODES
      CALL CONDUC (NS,ZS,TK,CON)
C
C**** CALCULATE THE EFFECTIVE SPECIFIC HEAT AT EACH NODE
      IF (ITER .EQ. 1) CALL SOILHT (NS,CST,VLC,VIC,TS,MAT,CONC)
      CALL SOILHT (NS,CSDT,VLCDT,VICDT,TSDT,MATDT,CONCDT)
C
C**** OBTAIN AVERAGE SPECIFIC HEAT OVER THE TIME STEP
      CALL WEIGHT (NS,CS,CST,CSDT)
C
C
C**** DETERMINE THE COEFFICIENTS FOR THE SURFACE LAYER
C
      A1(N+1) = WDT*(CON(1) + DELTA(1)*(RHOL*CL*QSL(1)+CV*QSV(1)))
      B1(N) = B1(N)- WDT*(CON(1) + EPSLON(1)*(RHOL*CL*QSL(1)+CV*QSV(1)))
     >             - (ZS(2)-ZS(1))*CSDT(1)/(2.*DT)
      C1(N) = WDT*(CON(1) + EPSLON(1)*(RHOL*CL*QSL(1)+CV*QSV(1)))
      D1(N) = D1(N)-(CON(1) + EPSLON(1)*(RHOL*CL*QSL(1)+CV*QSV(1)))
     >           *(WT*(TS(1)-TS(2))+ WDT*(TSDT(1)-TSDT(2))) 
     >        - LV*QSV(1) + SS(1) - (ZS(2)-ZS(1))/(2.*DT)
     >           *(CS(1)*(TSDT(1)-TS(1)) - RHOI*LF*(VICDT(1)-VIC(1)))
C
C**** CHECK IF ICE IS PRESENT AT THE END OF THE TIME STEP
      IF (ICESDT(1) .NE. 1) GO TO 20
C
C**** ICE IS PRESENT IN LAYER - - ADJUST COEFFICIENTS FOR LATENT HEAT
C**** TRANSFER AND THE SLOPE OF THE WATER CONTENT-TEMPERATURE CURVE
      IF (ICES(1) .EQ. 1) THEN
C        ICE IS PRESENT FOR THE ENTIRE TIME STEP
         T=273.16+TS(1)
        ELSE
C        ICE IS PRESENT ONLY AT THE END OF THE TIME STEP - DETERMINE
C        THE TEMPERATURE AT WHICH THE SOIL WILL BEGIN TO FREEZE
         TLCONC=0.0
         DO 15 K=1,NSALT
            TLCONC=TLCONC + CONC(K,1)
   15    CONTINUE
         TMPFRZ=273.16*LF/G/(LF/G-MAT(1)+TLCONC*UGAS*(TS(1)+273.17)/G)
         T=TMPFRZ
      END IF
      TDT=273.16+TSDT(1)
C     DETERMINE SLOPE OF LIQUID CONTENT-TEMPERATURE CURVE
      CALL FSLOPE (1,DLDTDT,DUMMY,TDT,MATDT,CONCDT,VLCDT)
      CALL FSLOPE (1,DLDT,DUMMY,T,MAT,CONC,VLC)
C     ENTER MATVLC FOR SLOPE OF LIQUID-MATRIC POTENTIAL CURVE
      CALL MATVL3 (1,MATDT(1),VLCDT(1),DLDM)
      B1(N)=B1(N) - (ZS(2)-ZS(1))/(2.*DT)* 0.5*RHOL*LF*(DLDT+DLDTDT)
     >     - WDT*RHOL*LF*CONH(1)*DLDTDT/DLDM
C
C**** DETERMINE THE MATRIX COEFFICIENTS FOR THE REST OF THE PROFILE
C     DETLYR (I)   average conductance (stomatal plus leaf boundary layer) for transpiration
C             in canopy layer i, which is equal to the derivative (m/s)
   20 DO 30 I=N+1,N+NS-2
         J=I-N+1
         A1(I+1)=WDT*(CON(J) + DELTA(J)*(RHOL*CL*QSL(J)+CV*QSV(J)))
         B1(I)= -WDT*(CON(J-1)+CON(J)
     >         +(DELTA(J-1)*(RHOL*CL*QSL(J-1)+CV*(QSV(J-1)))
     >         + EPSLON(J)*(RHOL*CL*QSL(J)+CV*QSV(J))))
     >         - (ZS(J+1)-ZS(J-1))*CSDT(J)/(2.*DT)
         C1(I)=WDT*(CON(J) + EPSLON(J)*(RHOL*CL*QSL(J)+CV*QSV(J)))
         D1(I)=(CON(J-1) + DELTA(J-1)*(RHOL*CL*QSL(J-1) + CV*QSV(J-1)))
     >            *(WT*(TS(J-1)-TS(J)) + WDT*(TSDT(J-1)-TSDT(J)))
     >        -(CON(J) + EPSLON(J)*(RHOL*CL*QSL(J) + CV*QSV(J)))
     >            *(WT*(TS(J)-TS(J+1)) + WDT*(TSDT(J)-TSDT(J+1)))
     >        -LV*(QSV(J)-QSV(J-1)) + SS(J)
     >        -(ZS(J+1)-ZS(J-1))/(2.*DT)
     >            *(CS(J)*(TSDT(J)-TS(J)) - RHOI*LF*(VICDT(J)-VIC(J)))
C
C****    CHECK IF ICE IS PRESENT AT THE END OF THE TIME STEP
         IF (ICESDT(J) .NE. 1) GO TO 30
C
C****    ICE IS PRESENT IN LAYER - - ADJUST COEFFICIENTS FOR LATENT HEAT
C****    TRANSFER AND THE SLOPE OF THE WATER CONTENT-TEMPERATURE CURVE
         IF (ICES(J) .EQ. 1) THEN
C           ICE IS PRESENT FOR THE ENTIRE TIME STEP
            T=273.16+TS(J)
          ELSE
C           ICE IS PRESENT ONLY AT THE END OF THE TIME STEP - DETERMINE
C           THE TEMPERATURE AT WHICH THE SOIL WILL BEGIN TO FREEZE
            TLCONC=0.0
            DO 25 K=1,NSALT
               TLCONC=TLCONC + CONC(K,J)
   25       CONTINUE
            TMPFRZ=273.16*LF/G
     >             /(LF/G-MAT(J)+TLCONC*UGAS*(TS(J)+273.16)/G)
            T=TMPFRZ
         END IF
         TDT=273.16+TSDT(J)
C        DETERMINE SLOPE OF LIQUID CONTENT-TEMPERATURE CURVE
         CALL FSLOPE (J,DLDTDT,DUMMY,TDT,MATDT,CONCDT,VLCDT)
         CALL FSLOPE (J,DLDT,DUMMY,T,MAT,CONC,VLC)
C        ENTER MATVLC FOR SLOPE OF LIQUID-MATRIC POTENTIAL CURVE
         CALL MATVL3 (J,MATDT(J),VLCDT(J),DLDM)
         B1(I)=B1(I) - (ZS(J+1)-ZS(J-1))/(2.*DT)
     >        *0.5*RHOL*LF*(DLDT+DLDTDT)
     >       - WDT*RHOL*LF*(CONH(J-1)+CONH(J))*DLDTDT/DLDM
   30 CONTINUE
      N=N+NS-2
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE SOILTK (NS,TK,VLC,VIC)
C
C     THIS SUBROUTINE CALCULATES THE THERMAL CONDUCTIVITY OF THE SOIL
C     LAYERS USING DE VRIES'S METHOD.
C
C***********************************************************************
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /SLPARM/ B(50),ENTRY(50),RHOB(50),SAT(50),SATK(50),
     >                NSALT,SALTKQ(10,50),VAPCOF(50),VAPEXP(50),
     >                SAND(50),SILT(50),CLAY(50),OM(50)
      REAL TK(50),VLC(50),VIC(50)
      REAL LF,LV,LS
C
      SAVE IFIRST,WFAIRD,WFSAD,WFSID,WFCLD,WFOMD,WFICED,TKMA,
     >     WFL,WFSA,WFSI,WFCL,WFOM,WFICE
C
C     WEIGHTING FACTOR FOR ICE IS NOW SET TO 1.0 TO CORRECT PROBLEMS 
C     POINTED OUT BY KENNEDY & SHARRATT, SOIL SCI. 163(8):636-645.
C     DATA GAI, GAOM, GASA, GASI, GAC /0.33, 0.5, 0.144, 0.144, 0.125/
      DATA GAOM, GASA, GASI, GAC /0.5, 0.144, 0.144, 0.125/
      DATA IFIRST/ 0/
C
C     DEFINE FUNCTION FOR CALCULATING DEVRIES WEIGHTING FACTORS
      DeVWF(TK0,TK1,GA) = 2.0/3.0/(1.0+(TK1/TK0-1.0)*GA)
     >                   + 1.0/3.0/(1.0+(TK1/TK0-1.0)*(1.0-2.0*GA))
C
      IF (IFIRST .EQ. 0) THEN
C        FIRST TIME INTO SUBROUTINE -- CALCULATE WEIGHTING FACTORS FOR
C        SOIL COMPONENTS WHEN SOIL IS COMPLETELY DRY AND WHEN WET
C     TKA          thermal conductivity of air (W/m/K)
C     TKSA         thermal conductivity of sand  fraction of soil (W/m/K)
         WFAIRD = 1.0
         WFSAD = DeVWF(TKA,TKSA,GASA)
         WFSID = DeVWF(TKA,TKSI,GASI)
         WFCLD = DeVWF(TKA,TKCL,GAC)
         WFOMD = DeVWF(TKA,TKOM,GAOM)
         WFICED = 1.0
C
C        SET THERMAL CONDUCTIVITY OF MOIST AIR EQUAL TO THAT OF DRY AIR
C        (VAPOR TRANSFER IS CALCULATED SEPARATELY AND INCLUDED IN HEAT
C        TRANSFER THROUGH SOIL; THEREFORE IT SHOULD NOT BE COMPENSATED
C        FOR IN COMPUTATION FOR THERMAL CONDUCTIVITY AS IS USUALLY DONE
C        WHEN USING DEVRIES METHOD.)
         TKMA = TKA
         WFL = 1.0
         WFSA = DeVWF(TKL,TKSA,GASA)
         WFSI = DeVWF(TKL,TKSI,GASI)
         WFCL = DeVWF(TKL,TKCL,GAC)
         WFOM = DeVWF(TKL,TKOM,GAOM)
         WFICE = 1.0
         IFIRST=1
      END IF
C
      DO 10 I=1,NS
C       CORRECTION OF TEXTURE FOR ORGANIC MATTER
        XSM = 1.0 - OM(I)
        XSAND=XSM*SAND(I)
        XSILT=XSM*SILT(I)
        XCLAY=XSM*CLAY(I)
C
C       CONVERT WEIGHT FRACTION TO VOLUMETRIC FRACTION
        VSAND = XSAND*RHOB(I)/RHOM
        VSILT = XSILT*RHOB(I)/RHOM
        VCLAY = XCLAY*RHOB(I)/RHOM
        VOM = OM(I)*RHOB(I)/RHOOM
        PORO = 1.-VSAND-VSILT-VCLAY-VOM
        VAC = PORO-VIC(I)-VLC(I)
        IF (VAC .LT. 0.0) VAC=0.0
C
C       DETERMINE LIMIT OF DEVRIES METHOD FROM TEXTURE
C       SAND ==> 0.05;  LOAM ==> 0.10;  CLAY ==> 0.15
        VLMT = 0.10 + 0.2*CLAY(I) - 0.1*SAND(I)
        IF (VLMT .LT. 0.05) VLMT = 0.05
        IF (VLMT .GT. 0.15) VLMT = 0.15
C
        IF (VLC(I) .GT. VLMT) THEN
C         MOIST SOIL -- USE DEVRIES METHOD DIRECTLY
          GAAIR = 0.035+0.298*(VLC(I)-VLMT)/(PORO-VLMT)
          WFAIR= DeVWF(TKL,TKMA,GAAIR)
          TK(I)=(WFSA*VSAND*TKSA + WFL*VLC(I)*TKL + WFICE*VIC(I)*TKI
     >    + WFAIR*VAC*TKMA + WFSI*VSILT*TKSI + WFCL*VCLAY*TKCL
     >    + WFOM*VOM*TKOM)  /  (WFSA*VSAND + WFL*VLC(I) + WFICE*VIC(I)
     >    + WFAIR*VAC + WFSI*VSILT + WFCL*VCLAY + WFOM*VOM)
        ELSE
C         INTERPOLATE THERMAL CONDUCTIVITY BETWEEN WATER CONTENT AT 0
C         AND LIMIT OF DEVRIES METHOD
          TK0=1.25*(WFSAD*VSAND*TKSA + WFICED*VIC(I)*TKI 
     >    + WFAIRD*(VAC+VLC(I))*TKA
     >    + WFSID*VSILT*TKSI + WFCLD*VCLAY*TKCL + WFOMD*VOM*TKOM)
     >    /(WFSAD*VSAND + WFICED*VIC(I) + WFAIRD*(VAC+VLC(I))
     >    + WFSID*VSILT + WFCLD*VCLAY + WFOMD*VOM)
C
          WFAIR= DeVWF(TKL,TKMA,0.035)
          TKLMT =(WFSA*VSAND*TKSA + WFL*VLMT*TKL + WFICE*VIC(I)*TKI
     >    + WFAIR*VAC*TKMA + WFSI*VSILT*TKSI + WFCL*VCLAY*TKCL
     >    + WFOM*VOM*TKOM) / (WFSA*VSAND + WFL*VLMT + WFICE*VIC(I)
     >    + WFAIR*VAC + WFSI*VSILT + WFCL*VCLAY + WFOM*VOM)
C
          TK(I) = TK0 + (TKLMT-TK0)*VLC(I)/VLMT
        ENDIF
   10 CONTINUE
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE SOILHT (NS,CS,VLC,VIC,TS,MAT,CONC)
C     THIS SUBROUTINE CALCULATES THE SPECIFIC HEAT THE SOIL LAYERS
C     -- THE LATENT HEAT OF VAPORIZATION IS INCLUDED.
C
C***********************************************************************
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /SLPARM/ B(50),ENTRY(50),RHOB(50),SAT(50),SATK(50),
     >                NSALT,SALTKQ(10,50),VAPCOF(50),VAPEXP(50),
     >                SAND(50),SILT(50),CLAY(50),OM(50)
      REAL LF,LV,LS
C
      REAL CS(50),VLC(50),VIC(50),TS(50),MAT(50),CONC(10,50)
C
      DO 20 I=1,NS
C        CORRECTION OF MINERAL FRACTION FOR ORGANIC MATTER 
         XSM = 1.0 - OM(I)
C
         VAC = 1. -XSM*RHOB(I)/RHOM -OM(I)*RHOB(I)/RHOOM -VLC(I) -VIC(I)
         IF (VAC .LT. 0.0) VAC=0.0
         CS(I)=XSM*RHOB(I)*CM + OM(I)*RHOB(I)*COM + VLC(I)*RHOL*CL 
     >          + VIC(I)*RHOI*CI + VAC*RHOA*CA
C
C****    INCLUDE THE LATENT HEAT OF VAPORIZATION IN THE HEAT CAPACITY
C        IF LAYER IS NOT SATURATED
         IF (MAT(I).LT.ENTRY(I) .AND. VAC.GT.0.0) THEN
C           DETERMINE HUMIDITY OF LAYER FROM TOTAL WATER POTENTIAL
            TLCONC=0.0
            DO 10 J=1,NSALT
               TLCONC= TLCONC+CONC(J,I)
   10       CONTINUE
            TOTPOT=MAT(I) - TLCONC*UGAS*TS(I)/G
            HUMID=EXP(.018*G/UGAS/(TS(I)+273.16)*TOTPOT)
C           OBTAIN SLOPE OF VAPOR DENSITY CURVE
            CALL VSLOPE (S,RHOV,TS(I))
            CS(I)=CS(I)+VAC*LV*HUMID*S
         END IF
   20 CONTINUE
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE FSLOPE(I,DLDT,PTDLDT,TMP,MAT,CONC,VLC)
C
C     THIS SUBROUTINE DETERMINES THE SLOPE OF THE LIQUID CONTENT -
C     TEMPERATURE CURVE (DLDT) AS WELL AS THE PARTIAL OF T*DLDT.
C     (OSMOTIC EFFECTS ARE INCLUDED.)
C     ---- TMP MUST BE IN DEGREES KELVIN
C
C***********************************************************************
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /SLPARM/ B(50),ENTRY(50),RHOB(50),SAT(50),SATK(50),
     >                NSALT,SALTKQ(10,50),VAPCOF(50),VAPEXP(50),
     >                SAND(50),SILT(50),CLAY(50),OM(50)
      REAL MAT(50),CONC(10,50),VLC(50)
      REAL LF,LS,LV
C
      CNCGAS=0.0
      CRT=0.0
      DO 10 J=1,NSALT
         DKQ=SALTKQ(J,I) + VLC(I)*RHOL/RHOB(I)
         CNCGAS=CNCGAS + CONC(J,I)*UGAS
         CRT=CRT + CONC(J,I)*UGAS*TMP*RHOL/RHOB(I)/DKQ
   10 CONTINUE
      CALL MATVL3 (I,MAT(I),VLC(I),DLDM)
      IF (DLDM .GT. 0.0) THEN
         DLDT=(LF/TMP + CNCGAS)/(G/DLDM + CRT)
        ELSE
C        NO CHANGE IN WATER CONTENT WITH MATRIC POTENTIAL OR TEMPERATURE
C        LAYER IS PROBABLY SATURATED
         DLDT=0.0
      END IF
C
C**** THIS PART OF THE SUBROUTINE CALCULATES THE PARTIAL OF T*DLDT
C
      CRKQ=0.0
      DO 20 J=1,NSALT
         DKQ=SALTKQ(J,I) + VLC(I)*RHOL/RHOB(I)
         CRKQ=CRKQ +
     >          CONC(J,I)*UGAS/DKQ*(1.-2.*TMP*(RHOL/RHOB(I))*DLDT/DKQ)
   20 CONTINUE
CXXXX PTDLDT=-LF/ ((-B(I)*G*MAT(I)/VLC(I) + CRT)**2)
CXXXX>      * (B(I)*G*MAT(I)*(B(I)+1)*DLDT/(VLC(I)**2) + CRKQ)
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE FROZEN (I,VLC,VIC,MAT,CONC,TS,SALT,ICES)
C
C     THIS SUBROUTINE DETERMINES THE LIQUID AND ICE WATER CONTENT OF THE
C     SOIL FOR TMP < 0 C  AND DETERMINES IF ANY ICE MAY BE PRESENT.
C     IT THEN UPDATES THE MATRIC POTENTIAL AND SOLUTE CONCENTRATION.
C
C     IN: TS,I
C     INOUT:VLC,VIC,MAT,CONC,ICES

C***********************************************************************
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /SLPARM/ B(50),ENTRY(50),RHOB(50),SAT(50),SATK(50),
     >                NSALT,SALTKQ(10,50),VAPCOF(50),VAPEXP(50),
     >                SAND(50),SILT(50),CLAY(50),OM(50)
      REAL VLC(50),VIC(50),MAT(50),CONC(10,50),TS(50),SALT(10,50)
      INTEGER ICES(50)
      REAL LF,LS,LV
      REAL MAT1
C
C     ITERATE TO FIND THE MAXIMUM WATER CONTENT FROM THE TOTAL, MATRIC
C     AND OSMOTIC POTENTIALS:
C               TOTAL - MATRIC - OSMOTIC = 0 => F
C     MATRIC AND OSMOTIC POTENTIALS ARE FUNCTIONS OF WATER CONTENT.
C     SOLVE BY NEWTON-RAPHSON ITERATION  =>  D(VLC) = -F/(DF/D(VLC)
C
      ICONV=0
    5 TMP=TS(I)+273.16
      TOTPOT = LF*TS(I)/TMP/G
      CALL MATVL2 (I,TOTPOT,VLC1,DUMMY)
C
C**** BEGIN ITERATIONS
      ITER=0
C     CALCULATE OSMOTIC POTENTIAL AND DERIVATIVE WITH RESPECT TO LIQUID
   10 OSMPOT=0.0
      DERIV=0.0
      DO 20 J=1,NSALT
         TERM = SALT(J,I)*UGAS*TMP/(SALTKQ(J,I)+VLC1*RHOL/RHOB(I))
         OSMPOT= OSMPOT-TERM
         DERIV= DERIV+TERM*RHOL/RHOB(I)/(SALTKQ(J,I)+VLC1*RHOL/RHOB(I))
   20 CONTINUE
      OSMPOT= OSMPOT/G
      DERIV= DERIV/G
C
C     CALCULATE MATRIC POTENTIAL AND DERIVATIVE WITH RESPECT TO LIQUID
      CALL MATVL1 (I,MAT1,VLC1,DUMMY)
      CALL MATVL3 (I,MAT1,VLC1,DLDM)
C
C     DERTERMINE ERROR IN WATER POTENTIAL AND ADJUST WATER CONTENT
      ERROR = TOTPOT - OSMPOT - MAT1
      DELTA = ERROR/(DERIV + 1/DLDM)
      VLC2 = VLC1 + DELTA
      IF (VLC2 .GT. SAT(I)) VLC2=SAT(I)
      IF (VLC2 .LE. VLC1/2) VLC2=VLC1/2.
      DELTA = VLC2 - VLC1
      VLC1=VLC2
      IF (ABS(DELTA) .LT. .00001) GO TO 30
      ITER=ITER+1
      IF (ICONV .GT. 0)
     >WRITE (21,*) ITER,OSMPOT,MAT1,VLC1,DERIV,1/DLDM,ERROR
      IF (ITER .GT. 30) THEN
         ICONV = ICONV+1
         WRITE (21,100)
         WRITE (21,*) I,TS(I),VLC(I),VIC(I),TOTPOT
         IF (ICONV .EQ. 1) GO TO 5
         STOP
      END IF
      GO TO 10
C
C**** IF ACTUAL LIQUID PLUS ICE CONTENT IS LESS THAN THAT CALCULATED
C     FROM ABOVE, THERE IS NO ICE PRESENT IN THE LAYER.  (THE TOTAL
C     POTENTIAL IS SUFFICIENTLY NEGATIVE THAT THE WATER WILL NOT
C     FREEZE.)
   30 IF ((VLC(I)+VIC(I)*RHOI/RHOL).LT.VLC2) THEN
C*       NO ICE IS PRESENT AND LAYER IS UNSATURATED
         ICES(I)=0
         VLC(I)=VLC(I) + VIC(I)*RHOI/RHOL
         VIC(I)=0.0
         CALL MATVL1 (I,MAT(I),VLC(I),DUMMY)
        ELSE
C*       ICE MAY BE PRESENT; CONVERT ANY WATER PRESENT ABOVE THE MAXIMUM
C        LIQUID WATER CONTENT TO ICE
         VIC(I)=VIC(I)+(VLC(I)-VLC2)*RHOL/RHOI
         VLC(I)=VLC2
         IF (VIC(I) .GT. 0.0) THEN
C           ICE IS PRESENT --> VLC AND MAT ARE FUNCTION OF TEMP
            ICES(I)=1
C           COMPUTE MATRIC POTENTIAL FROM LIQUID WATER CONTENT
            CALL MATVL1 (I,MAT(I),VLC(I),DUMMY)
          ELSE
C           NO ICE IS PRESENT AND WATER CONTENT IS AT SATURATION
            VIC(I)=0.0
            ICES(I)=0
         END IF
      END IF
C     REDEFINE SOLUTE CONCENTRATIONS
      DO 40 J=1,NSALT
         CONC(J,I)=SALT(J,I)/(SALTKQ(J,I) + VLC(I)*RHOL/RHOB(I))
   40 CONTINUE
      RETURN
  100 FORMAT (//10X,'*****  PROGRAM WAS STOPPED IN SUBROUTINE FROZEN DUE
     > TO CONVERGENCE PROBLEMS *****')
      END

C***********************************************************************
C
      SUBROUTINE ADJUST(I,VICDT,VLCDT,MATDT,CONCDT,TSDT,SALTDT,ICESDT)
C
C     THIS SUBROUTINE ADJUSTS THE SOIL LAYER TO ACCOUNT FOR THE LATENT
C     HEAT OF FUSION ON THE FIRST ITERATION THAT THE LAYER CROSSES THE
C     FREEZING POINT.
C
C***********************************************************************
      COMMON /TIMEWT/ WT,WDT,DT
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /SLPARM/ B(50),ENTRY(50),RHOB(50),SAT(50),SATK(50),
     >                NSALT,SALTKQ(10,50),VAPCOF(50),VAPEXP(50),
     >                SAND(50),SILT(50),CLAY(50),OM(50)
      COMMON /SPHEAT/ CS(50)
      REAL VICDT(50),VLCDT(50),MATDT(50),CONCDT(10,50),
     >     TSDT(50),SALTDT(10,50)
      INTEGER ICESDT(50)
      REAL LF,LS,LV
C
C**** INITIALLY ASSUME THAT ALL WATER IS LIQUID
      VLCDT(I)=VLCDT(I) + (RHOI/RHOL)*VICDT(I)
      CALL MATVL1 (I,MATDT(I),VLCDT(I),DUMMY)
      DO 5 J=1,NSALT
         CONCDT(J,I)=SALTDT(J,I)/(SALTKQ(J,I) + VLCDT(I)*RHOL/RHOB(I))
    5 CONTINUE
      VICDT(I)=0.0
C
C**** CALCULATE THE FREEZING POINT TEMPERATURE ACCORDING TO TOTAL WATER
C     POTENTIAL AT THE END OF THE TIME STEP
      TLCONC=0.0
      DO 10 J=1,NSALT
         TLCONC= TLCONC+CONCDT(J,I)
   10 CONTINUE
      TOTPOT=MATDT(I) - TLCONC*UGAS*273.16/G
      TMPFRZ=273.16*TOTPOT/(LF/G-TOTPOT)
      TMPF=TMPFRZ+273.16
C
C**** CALCULATE THE ENERGY AVAILABLE TO FREEZE WATER
      ENERGY=CS(I)*(TMPFRZ-TSDT(I))
C
C**** CALCULATE THE EFFECTIVE HEAT CAPACITY INCLUDING LATENT HEAT TERM
      CALL FSLOPE (I,DLDT,DUMMY,TMPF,MATDT,CONCDT,VLCDT)
      EFFCS=CS(I)+RHOL*LF*DLDT
C
C**** CALCULATE THE TEMPERATURE AT THE END OF THE TIME STEP
      TSDT(I)=TMPFRZ - ENERGY/EFFCS
C
C**** CALL SUBROUTINE FROZEN TO DETERMINE THE LIQUID, ICE, AND SOLUTES
      CALL FROZEN (I,VLCDT,VICDT,MATDT,CONCDT,TSDT,SALTDT,ICESDT)
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE WBCAN (N,NPLANT,NC,NSP,NR,NS,ZC,TC,TCDT,VAPC,VAPCDT,   
     >  WCAN,WCANDT,PCAN,PCANDT,QVC,TRNSP,MAT,MATDT,XTRACT,SWCAN,LWCAN,
     >  CANMA,CANMB,DCHAR,RSTOM0,RSTEXP,PLEAF0,RLEAF0,ITYPE,ICESDT,ITER)
C
C     THIS SUBROUTINE CALCULATES THE JACOBIAN MATRIX COEFFICIENTS FOR
C     THE CANOPY PORTION OF THE NEWTON-RAPHSON SOLUTION OF THE ENERGY
C     BALANCE
C
C***********************************************************************
      COMMON /TIMEWT/ WT,WDT,DT
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /MATRIX/ A1(150),B1(150),C1(150),D1(150),
     >                A2(150),B2(150),C2(150),D2(150)
      REAL LF,LS,LV
C
      REAL ZC(11),TC(11),TCDT(11),VAPC(11),VAPCDT(11),
     >     WCAN(10),WCANDT(10),PCAN(8),PCANDT(8),QVC(10),TRNSP(9),
     >     MAT(50),MATDT(50),XTRACT(50),SWCAN(9,10),LWCAN(9,10)
      REAL DCHAR(8),RSTOM0(8),RSTEXP(8),PLEAF0(8),RLEAF0(8)
      REAL CON(10),CONT(10),CONDT(10),
     >     HEATC(10),ETLYR(10),DETLYR(10),DHEATC(10),DTLDTC(10)
      INTEGER ITYPE(8)
      SAVE CONT
C
C**** CALCULATE TRANSPIRATION FROM CANOPY
C     QVC (I)      vapor flux between node i and i+1 of the canopy (kg/s/m^2)
      CALL LEAFT (NPLANT,NC,NS,ITER,ITYPE,TC,TCDT,VAPC,VAPCDT,
     > WCAN,WCANDT,PCAN,PCANDT,MAT,MATDT,TRNSP,XTRACT,SWCAN,LWCAN,
     > HEATC,ETLYR,DETLYR,CANMA,CANMB,DCHAR,RSTOM0,RSTEXP,PLEAF0,RLEAF0,
     > DHEATC,DTLDTC)
C
C**** DETERMINE THE EDDY CONDUCTANCE TERM BETWEEN NODES
      IF (ITER .EQ. 1) CALL CANTK (NC,CONT,TC,ZC)
      CALL CANTK (NC,CONDT,TCDT,ZC)
C
C**** CALCULATE THE AVERAGE CONDUCTANCE TERM OVER THE TIME STEP
      CALL WEIGHT (NC,CON,CONT,CONDT)
C
C**** CONVERT THERMAL EDDY CONDUCTANCE TO VAPOR CONDUCTANCE
C     AND CALCULATE LIQUID FLUX BETWEEN LAYERS
      DO 5 I=1,NC
         CON(I)=CON(I)/RHOA/CA
         QVC(I)=CON(I)*(WDT*(VAPCDT(I)-VAPCDT(I+1))
     >                 + WT*(VAPC(I)-VAPC(I+1)))
    5 CONTINUE
C
C**** DETERMINE THE MATRIX COEFFICIENTS FOR THE TOP LAYER
C     DETLYR (I)   average conductance (stomatal plus leaf boundary layer) for transpiration
C             in canopy layer i, which is equal to the derivative (m/s)
      A2(N+1)=WDT*CON(1)
      B2(N)=B2(N) - WDT*(CON(1) + DETLYR(1)) - (ZC(2)-ZC(1))/2/DT
      C2(N)=WDT*CON(1)
      D2(N)=D2(N)-QVC(1)+ETLYR(1)-(ZC(2)-ZC(1))*(VAPCDT(1)-VAPC(1))/2/DT
C
C**** DETERMINE THE COEFFICIENTS FOR THE REST OF THE LAYERS
      DO 10 I=N+1,N+NC-1
         J=I-N+1
         A2(I+1)=WDT*CON(J)
         B2(I)=-WDT*(CON(J-1) + CON(J) + DETLYR(J))
     >         - (ZC(J+1)-ZC(J-1))/2/DT
         C2(I)=WDT*CON(J)
         D2(I)=QVC(J-1) - QVC(J) + ETLYR(J)
     >        - (ZC(J+1)-ZC(J-1))*(VAPCDT(J)-VAPC(J))/2/DT
   10 CONTINUE
      N=N+NC
C
C**** DETERMINE THE COEFFICIENTS FOR THE TOP LAYER OF THE NEXT MATERIAL
C
      IF (NSP.GT.0) THEN
C****    NEXT MATERIAL IS SNOW -- NO WATER BALANCE MATRIX FOR SNOW
         A2(N) = 0.0
         B2(N) = 0.0
         C2(N-1)=0.0
         D2(N) = 0.0
        ELSE
C
      IF (NR.GT.0) THEN
C****    NEXT MATERIAL IS RESIDUE
         B2(N) = -WDT*CON(NC)
         D2(N) = QVC(NC)
        ELSE
C
C**** NEXT MATERIAL IS SOIL -- WATER BALANCE BASED ON WATER CONTENT AND
C     MATRIC POTENTIAL
      A2(N) = A2(N)/RHOL
      IF (ICESDT .EQ. 0) THEN
C       NO ICE PRESENT AT SOIL SURFACE - WATER BALANCE FOR MATRIC POT.
        B2(N) = -WDT*CON(NC)*VAPCDT(NC+1)*0.018*G
     >           /(UGAS*(TCDT(NC+1)+273.16)*RHOL)
        C2(N-1)=C2(N-1)*VAPCDT(NC+1)*0.018*G/(UGAS*(TCDT(NC+1)+273.16))
       ELSE
C       ICE IS PRESENT AT SOIL SURFACE - WATER BALANCE FOR ICE CONTENT
        B2(N)=0.0
        C2(N-1)=0.0
      END IF
      D2(N) = QVC(NC)/RHOL
C
      END IF
      END IF
C
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE SNOWBC (N,NSP,NR,ZSP,QVSP,VAPSPT,TSDT,TSPDT,ICESDT)
     >                   
C
C     THIS SUBROUTINE SETS UP THE UPPER BOUNDARY CONDITION FOR THE WATER
C     BALANCE OF THE RESIDUE-SOIL SYSTEM WHEN THERE IS A SNOWPACK
C***********************************************************************
      COMMON /TIMEWT/ WT,WDT,DT
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /MATRIX/ A1(150),B1(150),C1(150),D1(150),
     >                A2(150),B2(150),C2(150),D2(150)
      COMMON /SPPROP/ G1,G2,G3,EXTSP,TKSPA,TKSPB,TKSPEX,VDIFSP,VAPSPX
      REAL LF,LS,LV
C
      REAL ZSP(100),QVSP(100),TSPDT(100),TSDT(50)
C
C**** DETERMINE THE VAPOR DIFFUSIVITY OF BOTTOM SNOWPACK NODE
      VDIF= VDIFSP*(P0/PRESUR)*(1.+TSPDT(NSP)/273.16)**VAPSPX
      CON = VDIF/(ZSP(NSP+1)-ZSP(NSP))
C
C**** DEFINE COEFFICIENTS
      IF (NR .GT. 0) THEN
C        CALCULATE BOUNDARY COEFFICIENTS FOR SURFACE RESIDUE NODE
         B2(N) = -WDT*CON
         D2(N) = QVSP(NSP)
       ELSE
C        CALCULATE BOUNDARY COEFFICIENTS FOR SURFACE SOIL NODE
         IF (ICESDT .EQ. 0) THEN
C           NO ICE PRESENT AT SOIL SURFACE - MATRIC POTENTIAL UNKNOWN
            B2(N) = -WDT*CON*VAPSPT*0.018*G/(UGAS*(TSDT(1)+273.16))/RHOL
           ELSE
C           ICE PRESENT AT SOIL SURFACE - ICE CONTENT UNKNOWN
            B2(N) = 0.0
         END IF
         D2(N) = QVSP(NSP)/RHOL
      END IF
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE WBRES (N,NR,NSP,ZR,TR,TRDT,VAPR,VAPRDT,GMC,GMCDT,RHOR,
     >                  QVR,RHOSP,ICESDT,ITER)
C
C     THIS SUBOUTINE CALCULATES THE NEWTON-RAPHSON COEFFICIENTS FOR THE
C     WATER BALANCE OF THE RESIDUE LAYERS.
C
C***********************************************************************
      COMMON /TIMEWT/ WT,WDT,DT
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /MATRIX/ A1(150),B1(150),C1(150),D1(150),
     >                A2(150),B2(150),C2(150),D2(150)
      COMMON /RESIDU/ EVAP(10),EVAPK(10)
      REAL LF,LS,LV
C
      REAL ZR(10),TR(10),TRDT(10),VAPR(10),VAPRDT(10),GMC(10),GMCDT(10),
     >     RHOR(10),QVR(10),RHOSP(100),CONV(10),VAPCON(10),CONVEC(10),
     >     TKRES(10)
C
C**** DETERMINE THE HEAT TRANSFER COEFFICIENT FOR EACH NODE
      CALL RESTK (NR,NSP,TKRES,CONVEC,TR,TRDT,GMC,GMCDT,RHOR,RHOSP)
C
C**** DETERMINE THE VAPOR TRANSPORT FROM THE THERMAL CONVECTION
      CALL RESVK (NR,NSP,TR,TRDT,CONVEC,VAPCON)
      NR1=NR+1
      VAPCON(NR1)=VAPCON(NR)
C
C**** DETERMINE THE CONDUCTANCE TERM FOR CONVECTIVE VAPOR TRANSPORT
      CALL CONDUC (NR1,ZR,VAPCON,CONV)
C
C**** DETERMINE THE VAPOR FLUX BETWEEN RESIDUE NODES
      CALL QVRES (NR,QVR,CONV,VAPR,VAPRDT)
C
C**** DETERMINE THE COEFFICIENTS FOR THE SURFACE LAYER
      C2(N)=WDT*CONV(1)
      IF (NR.GT.1) THEN
         B2(N)=B2(N) - WDT*(CONV(1) + EVAPK(1)) - (ZR(2)-ZR(1))/(2*DT)
         D2(N)= D2(N) - QVR(1)- EVAP(1)
     >                - (ZR(2)-ZR(1))/(2*DT)*(VAPRDT(1)-VAPR(1))
        ELSE
         B2(N)=B2(N) - WDT*(CONV(1) + EVAPK(1)) - (ZR(2)-ZR(1))/DT
         D2(N)= D2(N) - QVR(1)- EVAP(1)
     >                - (ZR(2)-ZR(1))/DT*(VAPRDT(1)-VAPR(1))
      END IF
C
C**** DETERMINE THE COEFFICIENTS FOR THE REST OF THE RESIDUE LAYERS
      DO 10 I=N+1,N+NR-1
         J=I-N+1
         A2(I)=WDT*CONV(J-1)
         C2(I)=WDT*CONV(J)
         IF (J .NE. NR) THEN
            B2(I)=-WDT*(CONV(J-1)+CONV(J)+EVAPK(J))
     >            - (ZR(J+1)-ZR(J-1))/2./DT
            D2(I)= QVR(J-1) - QVR(J) - EVAP(J)
     >                   - (ZR(J+1)-ZR(J-1))/(2*DT)*(VAPRDT(J)-VAPR(J))
          ELSE
            B2(I)=-WDT*(CONV(J-1)+CONV(J)+EVAPK(J))
     >            - (ZR(J+1)-ZR(J)+(ZR(J)-ZR(J-1))/2)/DT
            D2(I)= QVR(J-1) - QVR(J) - EVAP(J)
     >       - (ZR(J+1)-ZR(J)+(ZR(J)-ZR(J-1))/2)/DT*(VAPRDT(J)-VAPR(J))
         END IF 
   10 CONTINUE
      N=N+NR
C
C**** DETERMINE THE COEFFICIENTS SOIL SURFACE
      A2(N)=WDT*CONV(NR)/RHOL
      IF (ICESDT .EQ. 0) THEN
C        NO ICE IS PRESENT AT SOIL SURFACE
         B2(N)=-WDT*CONV(NR)*VAPRDT(NR+1)
     >        *0.018*G/UGAS/(TRDT(NR+1)+273.16)/RHOL
         C2(N-1)=C2(N-1)*VAPRDT(NR+1)*.018*G/UGAS/(TRDT(NR+1)+273.16)
        ELSE
C        ICE IS PRESENT AT SOIL SURFACE - WATER BALANCE FOR ICE CONTENT
         B2(N)=0.0
         C2(N-1)=0.0
      END IF
      D2(N) = QVR(NR)/RHOL
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE WBSOIL (N,NS,ZS,TS,TSDT,MAT,MATDT,VLC,VLCDT,VIC,VICDT,
     >                   CONC,CONCDT,ICESDT,QSL,QSV,XTRACT,SEEP,US,ITER)
C
C     THIS SUBROUTINE CALCULATES THE COEFFICIENTS FOR THE SOIL IN THE
C     NEWTON-RAPHSON PROCEDURE OF THE ENERGY BALANCE
C
C***********************************************************************
      COMMON /TIMEWT/ WT,WDT,DT
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /SLPARM/ B(50),ENTRY(50),RHOB(50),SAT(50),SATK(50),
     >                NSALT,SALTKQ(10,50),VAPCOF(50),VAPEXP(50),
     >                SAND(50),SILT(50),CLAY(50),OM(50)
      COMMON /MATRIX/ A1(150),B1(150),C1(150),D1(150),
     >                A2(150),B2(150),C2(150),D2(150)
      REAL LF,LS,LV
C
      REAL ZS(50),TS(50),TSDT(50),MAT(50),MATDT(50),VLC(50),VLCDT(50),
     >     VIC(50),VICDT(50),CONC(10,50),CONCDT(10,50),QSL(50),QSLT(50),
     >     QSLDT(50),QSV(50),XTRACT(50),US(50),HKT(50),HKDT(50),
     >     CON(50),CONT(50),QSVT(50),QSVDT(50),DCONVP(50)
      INTEGER ICESDT(50),NODE(50)
      SAVE QSLT,QSVT
C
C     INITIALIZE A2(1) IF SOIL IS FIRST MATERIAL
      IF (N .EQ. 1) A2(N)=0.0
C
C**** DETERMINE THE AVERAGE VAPOR FLUX BETWEEN SOIL NODES
      IF (ITER.EQ.1) CALL QVSOIL (NS,QSVT,ZS,TS,MAT,VLC,VIC,CONC,DCONVP)
      CALL QVSOIL (NS,QSVDT,ZS,TSDT,MATDT,VLCDT,VICDT,CONCDT,DCONVP)
      CALL WEIGHT (NS-1,QSV,QSVT,QSVDT)
C
C**** DETERMINE THE LIQUID MOISTURE FLUX BETWEEN NODES
      IF (ITER .EQ. 1) THEN
C****    CALCULATE FLUX AT BEGINNING OF TIME STEP
C****    DETERMINE THE HYDRAULIC CONDUCTIVITY OF EACH SOIL NODE
         CALL SOILHK (NS,HKT,MAT,VLC,VIC)
C****    DETERMINE THE CONDUCTANCE TERM BETWEEN SOIL NODES
         CALL CONDUC (NS,ZS,HKT,CONT)
         CALL QLSOIL (NS,ZS,QSLT,CONT,MAT)
      END IF
      CALL SOILHK (NS,HKDT,MATDT,VLCDT,VICDT)
      CALL CONDUC (NS,ZS,HKDT,CON)
      CALL QLSOIL (NS,ZS,QSLDT,CON,MATDT)
C
C**** OBTAIN THE AVERAGE MOISTURE FLUX OVER THE TIME STEP
      CALL WEIGHT (NS-1,QSL,QSLT,QSLDT)
C
C     LIMIT NET WATER FLUX INTO NODES IF FLUX FILLS NODE BEYOND 
C     AVAILABLE POROSITY AS A RESULT DOWNWARD FLUX INTO LAYERS WITH
C     LIMITED PERMEABILITY (PERHAPS DUE TO ICE)
C     (START AT BOTTOM OF PROFILE AND WORK UP SO THAT ANY CHANGES IN THE
C     FLUXES CAN WORK THEIR UP THROUGH THE PROFILE)     
      IFLAG=0
      DO 5 I=NS-1,2,-1
         IF (VICDT(I).GT.0.0) THEN
            QMAX=(SAT(I) - VLC(I) - VIC(I))*(ZS(I+1) - ZS(I-1))/2./DT
            IF (QMAX.LT.0.0) QMAX=0.0
            IF (QSL(I-1)-QSL(I) .GT. QMAX) THEN
               IFLAG=1
               IF (QSL(I-1).GT.0.0) THEN
                  CON(I-1)=0.0
                  QSL(I-1)=QSL(I)+QMAX
                  IF (QSL(I-1).LT.0.0) THEN 
                     QSL(I-1)=0.0
                     QSL(I)=QSL(I-1)-QMAX
                     IF (QSL(I).GT.0.0) QSL(I)=0.0
                  END IF
                 ELSE
                  CON(I)=0.0
                  QSL(I)=QSL(I-1)-QMAX
                  IF (QSL(I).GT.0.0) QSL(I)=0.0
               END IF
            END IF
         END IF
    5 CONTINUE
      IF (VICDT(1).GT.0.0) THEN
         QMAX=(SAT(I) - VLC(1) - VIC(1))*(ZS(2) - ZS(1))/2./DT
         IF (QMAX.LT.0.0) QMAX=0.0
         IF (-QSL(1) .GT. QMAX) THEN
            IFLAG=1
            CON(1)=0.0
            QSL(1)=-QMAX
         END IF
      END IF
C
      IF (IFLAG.GT.0) THEN
      DO 10 I=2,NS
         IF (VICDT(I).GT.0.0) THEN
            QMAX=(SAT(I) - VLC(I) - VIC(I))*(ZS(I+1) - ZS(I-1))/2./DT
            IF (QMAX.LT.0.0) QMAX=0.0
            IF (QSL(I-1)-QSL(I) .GT. QMAX) THEN
               IF (QSL(I).LT.0.0) THEN
                  CON(I)=0.0
                  QSL(I)=QSL(I-1)-QMAX
                  IF (QSL(I).GT.0.0) THEN 
                     QSL(I)=0.0
                     QSL(I-1)=QSL(I)+QMAX
                     IF (QSL(I-1).LT.0.0) QSL(I-1)=0.0
                  END IF
                 ELSE
                  CON(I-1)=0.0
                  QSL(I-1)=QSL(I)+QMAX
                  IF (QSL(I-1).LT.0.0) QSL(I-1)=0.0
               END IF
            END IF
         END IF
   10 CONTINUE
      END IF
C
C     COUNT NUMBER OF SATURATED NODS
      ISAT=0
C
C**** DETERMINE THE COEFFICIENTS FOR THE SURFACE LAYER
      D2(N)=D2(N) - QSL(1) - QSV(1)/RHOL + US(1) - XTRACT(1)/RHOL -
     > (ZS(2)-ZS(1))*(VLCDT(1)-VLC(1)+RHOI/RHOL*(VICDT(1)-VIC(1)))/2./DT
      IF (ICESDT(1) .EQ. 1) THEN
C        ICE IS PRESENT -- CALCULATE COEFFICIENTS FOR ICE CONTENT
         B2(N)=-(ZS(2)-ZS(1))/(2.*DT)*RHOI/RHOL
         A2(N+1)=0.0
         SEEP=0.0
       ELSE
C        NO ICE IS PRESENT -- CALCULATE COEFFICIENTS FOR LIQUID BALANCE
         IF (MATDT(1).LT.0.0 .OR. D2(N).LT.0.0) THEN
C           NO SEEPAGE
            IF (MATDT(1) .GT. ENTRY(1)) THEN
C              NODE IS SATURATED - CHECK IF BOUNDARY IS POORLY DEFINED 
               IF (-B2(N).LT.CON(1)*1.E-6 .OR. CON(1).LE.0.0) THEN
C                 SOIL WATER FLOW OVERWHELMS SURFACE BOUNDARY CONDITION
                  ISAT=ISAT+1
                  NODE(ISAT)=1
               END IF
            END IF
            CALL MATVL3 (1,MATDT(1),VLCDT(1),DLDM)
            B2(N)=B2(N) - WDT*(CON(1)+DCONVP(1)/RHOL) 
     >            - (ZS(2)-ZS(1))/(2.*DT)*DLDM
            A2(N+1)= WDT*(CON(1)+DCONVP(1)/RHOL) 
            SEEP=0.0
           ELSE 
C           ALLOW FOR SEEPAGE FROM SURFACE NODE - BOUNDARY IS PROPERLY 
C           DEFINED (DO NOT INCLUDE AS A SATURATED NODE)
            B2(N)=1.0
            A2(N+1)= WDT*CON(1)
C           NET FLOW INTO NODE IS EQUAL TO SEEPAGE FROM SURFACE
            SEEP=D2(N)
            D2(N)=0.0
         END IF
      END IF
C
C**** DETERMINE THE COEFFICIENTS FOR THE REST OF THE PROFILE
      DO 20 I=N+1,N+NS-2
         J=I-N+1
         D2(I)= QSL(J-1) - QSL(J) + (QSV(J-1)-QSV(J))/RHOL + US(J)
     >          - XTRACT(J)/RHOL - (ZS(J+1)-ZS(J-1))
     >          /(2.*DT)*(VLCDT(J)-VLC(J) + RHOI/RHOL*(VICDT(J)-VIC(J)))
         IF (ICESDT(J) .EQ. 1) THEN
C           ICE IS PRESENT -- CALCULATE COEFFICIENTS FOR ICE CONTENT
            B2(I)=-(ZS(J+1)-ZS(J-1))/(2.*DT)*RHOI/RHOL
            C2(I-1)=0.0
            A2(I+1)=0.0
          ELSE
C           NO ICE IS PRESENT -- CALCULATE COEFF. FOR LIQUID BALANCE
            IF (MATDT(J) .GT. ENTRY(J)) THEN
C              NODE IS SATURATED
               ISAT=ISAT+1
               NODE(ISAT)=J
            END IF                              
            CALL MATVL3 (J,MATDT(J),VLCDT(J),DLDM)
            B2(I)= -WDT*(CON(J-1)+CON(J)+(DCONVP(J-1)+DCONVP(J))/RHOL)
     >             -(ZS(J+1)-ZS(J-1))/(2.*DT)*DLDM
            C2(I-1)=WDT*(CON(J-1)+DCONVP(J-1)/RHOL)
            A2(I+1)=WDT*(CON(J)+DCONVP(J)/RHOL)
         END IF
   20 CONTINUE
C
C     ADJUST MATRIX IF SEEPAGE OCCURS
      IF (SEEP.GT.0.0) THEN
         A2(N)=0.0
         C2(N)=0.0
      END IF
C
      IF (ISAT .GT. 0) THEN
C        SATURATED NODES ARE PRESENT - CHECK FOR MATRIX SINGULARITY
         JTOP=0
         NODE(ISAT+1)=NS
         DO 30 J=1,ISAT
            I=NODE(J)+N-1
C           CHECK FOR MATRIX DISCONTINUITY AND SAVE TOP NODE
            IF (NODE(J).EQ.1) THEN
               JTOP=NODE(J)
              ELSE
cxxxx          IF (CON(NODE(J)-1).LE.0.0) JTOP=NODE(J)
               IF (CON(NODE(J)-1).LE.-B2(I)*1.E-7) JTOP=NODE(J)
            END IF
            IF (JTOP.NE.0) THEN
cxxxx          IF (CON(NODE(J)).LE.0.0) THEN
               IF (CON(NODE(J)).LE.-B2(I)*1.E-7) THEN
C                SATURATED NODES SURROUNDED BY ZERO CONDUCTIVITY; 
C                SOLUTION IS UNDEFINED - RESET B2(I) OF TOP SATURATED
C                LAYER AS IF IT IS JUST UNDER AIR ENTRY POTENTIAL
                 ITOP=JTOP+N-1
                 CLOSE = 1.0001*ENTRY(JTOP)
                 CALL MATVL3 (JTOP,CLOSE,VLCDT(JTOP),DLDM)
                 IF (JTOP.EQ.1) THEN
                  B2(ITOP)=B2(ITOP)-(ZS(JTOP+1)-ZS(JTOP))/(2.*DT)*DLDM
                 ELSE
                  B2(ITOP)=B2(ITOP)-(ZS(JTOP+1)-ZS(JTOP-1))/(2.*DT)*DLDM
                 END IF
                 JTOP=0
               END IF
            END IF
C           IF NEXT NODE IS NOT SATURATED, RESET JTOP
            IF (NODE(J)+1.NE.NODE(J+1)) JTOP=0
   30    CONTINUE
      END IF
C
      N=N+NS-2
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE SOILHK (NS,HK,MAT,VLC,VIC)
C
C     THIS SUBROUTINE CALCULATES THE HYDRAULIC CONDUCTIVITY OF EACH SOIL
C     NODE
C
C***********************************************************************
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /SLPARM/ B(50),ENTRY(50),RHOB(50),SAT(50),SATK(50),
     >                NSALT,SALTKQ(10,50),VAPCOF(50),VAPEXP(50),
     >                SAND(50),SILT(50),CLAY(50),OM(50)
      REAL HK(50),MAT(50),VLC(50),VIC(50)
      REAL LF,LS,LV
C
      DO 10 I=1,NS
         IF (ENTRY(I) .GT. MAT(I)) THEN
            HK(I)=SATK(I)*(ENTRY(I)/MAT(I))**(2.+3/B(I))
          ELSE
            HK(I)=SATK(I)
         END IF
         IF (VIC(I) .GT. 0.0) THEN
C           LIMIT CONDUCTIVITY IF ICE CONTENT IS TOO LARGE (BASED ON
C           RESULTS FROM BLOOMSBURG AND WANG, 1969
C           (CALCULATE POROSITY ADJUSTING SPEC. DENSITY FOR ORG. MATTER
            POROS = 1. - RHOB(I)*((1.-OM(I))/RHOM+OM(I)/RHOOM)
            IF (POROS-VIC(I) .LT. 0.13) THEN
               HK(I) = 0.0
              ELSE
C              LINEARLY REDUCE CONDUCTIVITY TO ZERO AT
C              AN AVAILABE POROSITY OF 0.13
               FRACT=(POROS-VIC(I)-0.13)/(POROS-0.13)
               HK(I)=FRACT*HK(I)
            END IF
         END IF
   10 CONTINUE
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE QLSOIL (NS,ZS,QSL,CON,MAT)
C
C     THIS SUBROUTINE CALCULATES THE LIQUID MOISTURE FLUX BETWEEN
C     SOIL NODES
C
C***********************************************************************
      REAL QSL(50),ZS(50),CON(50),MAT(50)
C
      DO 10 I=1,NS-1
         QSL(I)=CON(I)*(MAT(I)-MAT(I+1)+ZS(I+1)-ZS(I))
   10 CONTINUE
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE BACKUP (NS,NR,NSP,NC,NPLANT,NSALT,ICES,ICESDT,TS,TSDT,
     > MAT,MATDT,CONC,CONCDT,VLC,VLCDT,VIC,VICDT,SALT,SALTDT,
     > TR,TRDT,VAPR,VAPRDT,GMC,GMCDT,TC,TCDT,VAPC,VAPCDT,WCAN,WCANDT,
     > PCAN,PCANDT,TSP,TSPDT,DLW,DLWDT,ICESP,ICESPT)
C
C     THIS SUBROUTINE SETS END OF TIME STEP VALUES BACK TO THOSE FOR
C     THE BEGINNING OF THE TIME STEP.
C
C***********************************************************************
      REAL TS(50),TSDT(50),MAT(50),MATDT(50),CONC(10,50),CONCDT(10,50),
     >     VLC(50),VLCDT(50),VIC(50),VICDT(50),SALT(10,50),SALTDT(10,50)
      REAL TR(10),TRDT(10),VAPR(10),VAPRDT(10),GMC(10),GMCDT(10)
      REAL TC(11),TCDT(11),VAPC(11),VAPCDT(11),WCAN(10),WCANDT(10),
     >     PCAN(8),PCANDT(8)
      REAL TSP(100),TSPDT(100),DLW(100),DLWDT(100)
      INTEGER ICES(50),ICESDT(50),ICESP(100),ICESPT(100)
C
C     SOIL PROPERTIES
      DO 15 I=1,NS
         TSDT(I)=TS(I)
         MATDT(I)=MAT(I)
         VLCDT(I)=VLC(I)
         VICDT(I)=VIC(I)
         ICESDT(I)=ICES(I)
         DO 10 J=1,NSALT
            SALTDT(J,I)=SALT(J,I)
            CONCDT(J,I)=CONC(J,I)
   10    CONTINUE
   15 CONTINUE
C
C     RESIDUE PROPERTIES
      DO 20 I=1,NR
         TRDT(I)=TR(I)
         VAPRDT(I)=VAPR(I)
         GMCDT(I)=GMC(I)
   20 CONTINUE
C
C     SNOWPACK PROPERTIES
      DO 30 I=1,NSP
         ICESPT(I) = ICESP(I)
         TSPDT(I) = TSP(I)
         DLWDT(I) = DLW(I)
   30 CONTINUE
C
C     CANOPY PROPERTIES
      DO 40 I=1,NC
         TCDT(I)=TC(I)
         VAPCDT(I)=VAPC(I)
         WCANDT(I)=WCAN(I)
   40 CONTINUE
      IF (NC .GT. 0) THEN
         DO 45 J=1,NPLANT
            PCANDT(J)=PCAN(J)
   45    CONTINUE
      END IF
C
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE SOLUTE (NS,ZS,CONC,CONCDT,SALT,SALTDT,VLC,VLCDT,TS,
     >                  TSDT,QSL,XTRACT,SINK,DGRADE,SLTDIF,ASALT,DISPER)
C
C     THIS SUBROUTINE SOLVES THE SOLUTE BALANCE OF THE SOIL FOR EACH
C     TYPE OF SOLUTE.  THE SOLUTION FOR THE SOLUTE BALANCE IS A LINEAR
C     SET OF EQUATIONS, AND CAN THEREFORE BE SOLVED DIRECTLY.
C     (NO NEWTON-RAPHSON ITERATION REQUIRED)
C
C***********************************************************************
      COMMON /TIMEWT/ WT,WDT,DT
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /SLPARM/ B(50),ENTRY(50),RHOB(50),SAT(50),SATK(50),
     >                NSALT,SALTKQ(10,50),VAPCOF(50),VAPEXP(50),
     >                SAND(50),SILT(50),CLAY(50),OM(50)
      REAL LF,LS,LV
C
      REAL DGRADE(10),SLTDIF(10),ASALT(50),DISPER(50)
      REAL ZS(50),CONC(10,50),CONCDT(10,50),SALT(10,50),SALTDT(10,50),
     >     VLC(50),VLCDT(50),TS(50),TSDT(50),QSL(50),XTRACT(50),
     >     SINK(10,50),DIFF(50),DIFFT(50),DIFFDT(50),CON(50),CONCNW(50),
     >     A3(50),B3(50),C3(50),D3(50),DELTA(50),EPSLON(50)
C
C**** SET UP LIQUID FLUX COEFFICIENTS SO THAT SOLUTES ARE CARRIED IN
C**** DIRECTION OF MOISTURE MOVEMENT.
      DO 10 I=1,NS-1
         IF (QSL(I) .GT. 0.0) THEN
            DELTA(I)=1.0
            EPSLON(I)=0.0
          ELSE
            DELTA(I)=0.0
            EPSLON(I)=1.0
         END IF
   10 CONTINUE
C
C**** START LOOP FOR EACH TYPE OF SOLUTE
      DO 50 J=1,NSALT
C        CALCULATE DIFFUSION COEFFICIENT FOR EACH LAYER
         CALL SALTK1 (NS,J,DIFFT,ZS,VLC,TS,QSL,SLTDIF,ASALT,DISPER)
         CALL SALTK1 (NS,J,DIFFDT,ZS,VLCDT,TSDT,QSL,SLTDIF,ASALT,DISPER)
C
C****    DETERMINE AVERAGE DIFFUSIVITY OVER TIME STEP, AND CONDUCTANCE
C****    TERM BETWEEN NODES
         CALL WEIGHT (NS,DIFF,DIFFT,DIFFDT)
         CALL CONDUC (NS,ZS,DIFF,CON)
         CALL SALTK2 (NS,J,CON,ZS,VLC,VLCDT,QSL,SLTDIF,ASALT,DISPER)
C
C****    DETERMINE THE COEFFICIENTS FOR THE SURFACE LAYER
         BB=RHOL*(CON(1)+QSL(1)*DELTA(1)) + XTRACT(1)
         CC=RHOL*(CON(1)-QSL(1)*EPSLON(1))
         B3(1)=-WDT*BB -  RHOB(1)*(ZS(2)-ZS(1))/(2.*DT)
     >             *(SALTKQ(J,1)+VLCDT(1)*RHOL/RHOB(1))
         C3(1)=WDT*CC
         D3(1)=-WT*CC*CONC(J,2)+(WT*BB-RHOB(1)*(ZS(2)-ZS(1))/(2.*DT)
     >           *(SALTKQ(J,1)+VLC(1)*RHOL/RHOB(1)))*CONC(J,1)+SINK(J,1)
C
C****    DETERMINE THE COEFFICIENTS FOR THE REST OF THE PROFILE
         DO 30 I=2,NS-1
           AA=RHOL*(CON(I-1)+QSL(I-1)*DELTA(I-1))
           BB=RHOL*(CON(I-1)-QSL(I-1)*EPSLON(I-1)
     >             +CON(I) + QSL(I)*DELTA(I)) + XTRACT(I)
           CC=RHOL*(CON(I)-QSL(I)*EPSLON(I))
           A3(I)=WDT*AA
           B3(I)=-WDT*BB -  RHOB(I)*(ZS(I+1)-ZS(I-1))/(2.*DT)
     >               *(SALTKQ(J,I)+VLCDT(I)*RHOL/RHOB(I))
           C3(I)=WDT*CC
           D3(I)=-WT*AA*CONC(J,I-1) - WT*CC*CONC(J,I+1) + SINK(J,I)
     >           + (WT*BB-RHOB(I)*(ZS(I+1)-ZS(I-1))/(2.*DT)
     >               *(SALTKQ(J,I)+VLC(I)*RHOL/RHOB(I)))*CONC(J,I)
C
   30    CONTINUE
C
C****    ADJUST D3(NS-1) FOR KNOWN BOUNDARY CONDITION AT NODE NS
C        D3(NS-1)=D3(NS-1) - C3(NS)*CONCDT(J,NS)
         B3(NS-1)=B3(NS-1) + C3(NS-1)
C
C****    SOLVE SET OF EQUATIONS FOR CONCENTRATION AT END OF TIME STEP
         CALL TDMA (NS-1,A3,B3,C3,D3,CONCNW)
         CONCNW(NS)=CONCNW(NS-1)
C
C****    UPDATE SOLUTE CONCENTRATIONS AND TOTAL SALTS FOR EACH NODE
         DO 40 I=1,NS
            CONCDT(J,I)=CONCNW(I)
            SALTDT(J,I)=(SALTKQ(J,I)+VLCDT(I)*RHOL/RHOB(I))*CONCDT(J,I)
C           ADJUST TOTAL SALTS AND SOLUTE CONCENTRATION FOR DEGRADATION
            SALTDT(J,I)=SALTDT(J,I)*EXP(-DGRADE(J)*DT/86400.)
            CONCDT(J,I)=SALTDT(J,I)/(SALTKQ(J,I)+VLCDT(I)*RHOL/RHOB(I))
   40    CONTINUE
   50 CONTINUE
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE SALTK (NS,J,DIFF,ZS,VLC,TS,QSL,SLTDIF,ASALT,DISPER)
C
C     THIS SUBROUTINE IS DIVIDED INTO TWO PARTS.  THE FIRST DETERMINES
C     THE MOLECULAR DIFFUSION OF SOLUTES AT EACH OF THE SOIL NODES.
C     AFTER WEIGHTING THE MOLECULAR DIFFUSION OVER TIME AND CALCULATING
C     THE CONDUCTANCE DUE TO MOLECULAR DIFFUSION BETWEEN NODES IN THE
C     MAIN PORTION OF THE PROGRAM, THE SECOND PORTION OF THE SUBROUTINE
C     MAY BE CALLED TO ADD THE CONDUCTANCE BETWEEN NODES DUE TO
C     HYDRODYNAMIC DISPERSION (WHICH IS A FUNCTION OF LIQUID FLUX
C     BETWEEN NODES).
C
C***********************************************************************
      COMMON /TIMEWT/ WT,WDT,DT
C
      REAL ZS(50),DIFF(50),VLC(50),VLCDT(50),TS(50),QSL(50),SLTDIF(10),
     >     ASALT(50),DISPER(50),CON(50)
C****
C     THIS PORTION OF THE SUBROUTINE CALCULATES THE MOLECULAR DIFFUSION
C     OF THE SOLUTES AT THE SOIL NODES.  (SINCE LIQUID FLUX IS
C     DETERMINED FOR BETWEEN NODES, THE HYDRODYNAMIC DISPERSION CANNOT
C     BE ADDED ON UNTIL THE MOLECULAR DIFFUSION IS DETERMINED FOR
C     BETWEEN NODES IN THE "CONDUC" SUBROUTINE)
C
      ENTRY SALTK1 (NS,J,DIFF,ZS,VLC,TS,QSL,SLTDIF,ASALT,DISPER)
      DO 10 I=1,NS
         DIFF(I)=SLTDIF(J)*ASALT(I)*VLC(I)**3*((TS(I)+273.16)/273.16)
   10 CONTINUE
      RETURN
C****
C     THIS PORTION OF THE SUBROUTINE ADDS THE CONDUCTANCE DUE TO
C     HYDRODYNAMIC DISPERSION TO THE CONDUCTANCE DUE TO MOLECULAR
C     DIFFUSION (BETWEEN NODES).
C
      ENTRY SALTK2 (NS,J,CON,ZS,VLC,VLCDT,QSL,SLTDIF,ASALT,DISPER)
      DO 20 I=1,NS-1
         AVGVLC=WT*(VLC(I)+VLC(I+1))/2. + WDT*(VLCDT(I)+VLCDT(I+1))/2.
         CON(I)=CON(I) + DISPER(I)*ABS(QSL(I))/AVGVLC/(ZS(I+1)-ZS(I))
   20 CONTINUE
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE SUMDT (NC,NPLANT,NSP,NR,NS,HFLUX,VFLUX,LWCAN,LWSNOW,
     >     LWRES,LWSOIL,SWCAN,SWSNOW,SWRES,SWSOIL,QVC,QVR,QVSP,QSL,QSV,
     >     TRNSP,XTRACT,SEEP,TSWSNO,TLWSNO,TSWCAN,TLWCAN,TSWRES,TLWRES,
     >     TSWSOI,TLWSOI,THFLUX,EVAP1,ETSUM,TSEEP,ROOTXT,TOTFLO,NTIMES,
     >     TOPSNO,TQVSP)
C
C     THIS SUBROUTINE SUMS THE NECESSARY FLUXES EACH TIME STEP FOR USE
C     IN THE WATER AND ENERGY BALANCE SUMMARIES AND IN THE SNOWPACK 
C     ADJUSTMENTS DUE FOR VAPOR TRANSFER
C
C***********************************************************************
      COMMON /TIMEWT/ WT,WDT,DT
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
C
      REAL LWCAN(9,10),LWSNOW(2),LWRES(10),LWSOIL,SWCAN(9,10),
     >     SWSNOW(100),SWRES(10)
C
      REAL QVC(10),QVR(10),QVSP(100),QSL(50),QSV(50),TRNSP(9),
     >     XTRACT(50),ROOTXT(50),TOTFLO(50),TQVSP(100),TRANSP(8)
      REAL LF,LS,LV
C
      SAVE TRANSP
C
      IF (NTIMES .EQ. 1) THEN
         TSWCAN=0.0
         TLWCAN=0.0
         TSWSNO=0.0
         TLWSNO=0.0
         TSWRES=0.0
         TLWRES=0.0
         TSWSOI=0.0
         TLWSOI=0.0
         THFLUX=0.0
         EVAP1=0.0
         ETSUM=0.0
         TSEEP=0.0
         TOPSNO=0.0
         DO 6 J=1,NPLANT
            TRANSP(J)=0.0
   6     CONTINUE
         DO 8 I=1,NSP
            TQVSP(I)=0.0
   8    CONTINUE
        DO 10 I=1,NS
            TOTFLO(I)=0.0
            ROOTXT(I)=0.0
  10    CONTINUE
      END IF
C
C     SUM RADIATION ABSORBED BY THE CANOPY (IN KILOJOULES)
      DO 20 I=1,NC
         TSWCAN = TSWCAN + SWCAN(NPLANT+1,I)*DT/1000.
         TLWCAN = TLWCAN + LWCAN(NPLANT+1,I)*DT/1000.
   20 CONTINUE
C
C     SUM RADIATION ABSORBED BY THE SNOWPACK (IN KILOJOULES)
      DO 30 I=1,NSP
         TSWSNO = TSWSNO + SWSNOW(I)*DT/1000.
   30 CONTINUE
      IF (NSP .GT. 0) TLWSNO = TLWSNO + (LWSNOW(1) + LWSNOW(2))*DT/1000.
C
C     SUM RADIATION ABSORBED BY THE RESIDUE (IN KILOJOULES)
      DO 40 I=1,NR
         TSWRES = TSWRES + SWRES(I)*DT/1000.
         TLWRES = TLWRES + LWRES(I)*DT/1000.
   40 CONTINUE
C
C     SUM RADIATION ABSORBED BY THE SOIL (IN KILOJOULES)
      TSWSOI = TSWSOI + SWSOIL*DT/1000.
      TLWSOI = TLWSOI + LWSOIL*DT/1000.
C
C     SUM SENSIBLE HEAT FLUX
      THFLUX = THFLUX + HFLUX*DT/1000.
C
C     SUM THE EVAPORATION, TRANSPIRATION FROM EACH PLANT SPECIES
      EVAP1 = EVAP1 + VFLUX/RHOL*DT
      DO 50 J=1,NPLANT
         TRANSP(J)=TRANSP(J) + TRNSP(J)*DT
   50 CONTINUE
      ETSUM = ETSUM + TRNSP(NPLANT+1)/RHOL*DT
C
C     SUM THE VAPOR FLUX OCCURRING IN, ABOVE AND BELOW THE SNOWPACK
      IF (NSP .GT. 0) THEN
         IF (NC .GT. 0) THEN
C           VAPOR FLUX ABOVE SNOWPACK IS FROM BOTTOM OF CANOPY
            TOPSNO = TOPSNO + QVC(NC)*DT
           ELSE
C           VAPOR FLUX ABOVE SNOWPACK IS FROM ATMOSPHERE
            TOPSNO = TOPSNO + VFLUX*DT
         END IF
         DO 60 I=1,NSP
            TQVSP(I) = TQVSP(I) + QVSP(I)*DT
   60    CONTINUE
C        TQVSP(NSP) IS THE VAPOR LEAVING BOTTOM OF SNOWPACK
      END IF
C
C     SUM THE ROOT EXTRACTION FROM EACH SOIL LAYER
      IF (NPLANT .GT. 0) THEN
         DO 70 I=1,NS
            ROOTXT(I)=ROOTXT(I)+XTRACT(I)*DT
   70    CONTINUE
      END IF
C
C     SUM THE TOTAL WATER FLUX BETWEEN SOIL NODES
      DO 80 I=1,NS-1
         TOTFLO(I)=TOTFLO(I) + (QSL(I) + QSV(I)/RHOL)*DT
   80 CONTINUE
C
C     SUM THE SEEPAGE FROM THE SURFACE NODE
      TSEEP=SEEP*DT
C
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE PRECP (NPLANT,NC,NSP,NR,NS,TA,TADT,HUM,HUMDT,
     > ZC,LANGLE,ITYPE,WCANDT,WCMAX,PCANDT,
     > ZSP,DZSP,TSPDT,DLW,DLWDT,RHOSP,TQVSP,TOPSNO,WLAG,STORE,SNOWEX,
     > ZR,TRDT,GMCDT,GMCMAX,RHOR,
     > ZS,TSDT,VLCDT,VICDT,MATDT,TOTFLO,SALTDT,CONCDT,ICESPT,ICESDT,
     > RAIN,POND,RUNOFF,EVAP1,PONDMX,SNOTMP,SNODEN,DIRRES,SLOPE)
C
C     THIS SUBROUTINE DETERMINES WHERE THE PRECIPITATION AND SNOWMELT
C     SHOULD GO
C
C***********************************************************************
      COMMON /TIMEWT/ WT,WDT,DT
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      REAL LF,LS,LV
C
      REAL ZC(11),WCANDT(10),PCANDT(8),ZSP(100),DZSP(100),TSPDT(100),
     >    DLW(100),DLWDT(100),RHOSP(100),TQVSP(100),WLAG(11),
     >    ZR(10),TRDT(10),GMCDT(10),RHOR(10),ZS(50),TSDT(50),VLCDT(50),
     >    VICDT(50),MATDT(50),TOTFLO(50),SALTDT(10,50),CONCDT(10,50)
      REAL MELT
      INTEGER LANGLE(8),ITYPE(8),ICESPT(100),ICESDT(50)
C
      AVGTA= WT*TA + WDT*TADT
      AVGHUM= WT*HUM + WDT*HUMDT
      TRAIN=AVGTA
C
C     RAINFALL INTERCEPTION BY THE CANOPY
      IF (RAIN.GT.0.0 .AND. NC.GT.0) 
     >   CALL RAINC (NPLANT,NC,LANGLE,ITYPE,RAIN,PCANDT,WCANDT,WCMAX,
     >               SLOPE)
C
      MELT = 0.0
      IF (RAIN .GT. 0.0) THEN
C        RAIN TEMPERATURE IS ASSUMED TO BE EQUAL TO WET BULB TEMPERATURE
         CALL WTBULB (TRAIN,AVGTA,AVGHUM)
         IF (AVGTA .LE. SNOTMP  .OR.  SNODEN .GT. 0.0) THEN
C           PRECIPITATION IS ASSUMED TO BE SNOW
C           (ADD ANY SNOW EXCESS REMAINING FROM ELIMINATING A SHALLOW
C           SNOWCOVER)
            SNOW = RAIN + SNOWEX
            POND = POND - SNOWEX
            SNOWEX = 0.0
            RAIN = 0.0
            IF (TRAIN .GT. 0.0) TRAIN = 0.0
            CALL NWSNOW (NSP,NC,NR,NS,NPLANT,ICESPT,ZSP,DZSP,RHOSP,
     >        TSPDT,DLW,DLWDT,TRAIN,SNOW,SNODEN,TQVSP,MELT,ZC,PCANDT,
     >        ZR,TRDT,GMCDT,RHOR,ZS,VLCDT,VICDT,TSDT,MATDT,CONCDT)
            IF (NSP .LE. 0) THEN
C              SNOW DID NOT STICK (MELT IS NOT INTERCEPTED BY CANOPY)
               TRAIN = 0.0
            END IF
           ELSE
C           PRECIP IS ASSUMED TO BE RAIN - DO NOT ALLOW RAIN TEMP < 0.0
            IF (TRAIN .LT. 0.0) TRAIN = 0.0
         END IF
      END IF
C
C     RAINFALL ABSORPTION BY SNOWPACK AND CALCULATION OF SNOWMELT
      IF (NSP .GT. 0) THEN
         CALL WBSNOW (NSP,ICESPT,ZSP,DZSP,RHOSP,TSPDT,DLW,DLWDT,RAIN,
     >                TRAIN,TOPSNO,EVAP1,TQVSP,WLAG,STORE,SCOUT)
         RAIN = SCOUT
         TRAIN = 0.0
      END IF
C
C     DO NOT INFILTRATE PONDED WATER OR SNOW EXCESS (SNOWEX) IF
C     TEMPERATURE IS BELOW FREEZING
      IF ((RAIN + POND + MELT).GT.0.0 .AND. TRAIN.GE.0.0) THEN
C
C        ADD MELT (FROM SNOW NOT STICKING) BACK INTO RAIN
         RAIN=RAIN+MELT
C
C        RAINFALL INTERCEPTION BY RESIDUE LAYERS
         IF(NR.GT.0) CALL RAINRS (NR,TRAIN,RAIN,ZR,TRDT,GMCDT,RHOR,
     >                            GMCMAX,DIRRES,SLOPE)
C
         IF (RAIN .GT. 0.0) THEN
C           ADD ANY EXCESS FROM ELIMINATING SHALLOW SNOWPACK BACK 
C           INTO RAIN
            RAIN=RAIN+SNOWEX
            POND=POND-SNOWEX 
            SNOWEX=0.0
         END IF
         RAIN = RAIN + POND
         POND = 0.0
C
C        RAINFALL INFILTRATION INTO THE SOIL
         CALL RAINSL (NS,TRAIN,RAIN,ZS,TSDT,VLCDT,VICDT,ICESDT,
     >            MATDT,TOTFLO,SALTDT,CONCDT,POND,PONDMX)
C
         IF (SNOWEX .GT. 0.0) THEN
C           DO NOT ALLOW EXCESS WATER FROM ELIMINATING SHALLOW SNOWPACK
C           TO RUNOFF - PUT THIS WATER BACK INTO PONDING AND SNOW EXCESS
C           (THIS MAY ALLOW PONDING DEPTH TO BECOME GREATER THAN PONDMX)
            POND = POND + RAIN
            RAIN = 0.0
            IF (POND .LT. SNOWEX) SNOWEX=POND
         END IF
      END IF
C
      IF (NSP .GT. 0) THEN
C
C     CALCULATE WATER EQUIVALENT - IF SUFFICIENTLY SMALL, ASSUME ENTIRE
C     SNOWPACK IS GONE
      WE = 0.0
      DO 20 I=1,NSP
         WE = WE + DZSP(I)*RHOSP(I)/RHOL + DLWDT(I)
   20 CONTINUE
      IF (WE .LT. 0.0005) THEN
C        WATER EQUIVALENT IS SUFFICIENTLY SMALL-ASSUME SNOWPACK IS GONE
         SCOUT = WE
         NSP = 0
C        ADD WATER CURRENTLY BEING LAGGED TO SNOW COVER OUTFLOW
         CALL SNOMLT (NSP,ICESPT,DZSP,RHOSP,TSPDT,DLWDT,WLAG,
     >                STORE,SCOUT)
C           
         WRITE (21,*) 
         WRITE (21,*)'SNOWCOVER IS ASSUMED GONE: WATER-EQUIVALENT ',
     >               'OF ICE AND TOTAL WATER =',WE,SCOUT,' METERS'
C
C        SET SNOW WATER EXCESS (SNOWEX) TO THE TOTAL WATER EQUIV.
C        REMAINING IN THE SNOW COVER - THIS WILL BE ADDED TO
C        PRECIPITATION OR ALLOWED TO INFITRATE AT A LATER TIME STEP
         SNOWEX=SNOWEX + SCOUT
C        INCLUDE SNOW EXCESS IN AMOUNT PONDED SO IT CAN BE INCLUDED IN
C        WATER BALANCE
         POND=POND + SCOUT
c
        ELSE
C        SNOWPACK STILL PRESENT -- DEFINE TEMPERATURE FOR BOTTOM BOUNDARY
         IF (NR .GT. 0) THEN
C           BOTTOM TEMPERATURE BOUNDARY IS TOP OF RESIDUE
            TSPDT(NSP+1) = TRDT(1)
          ELSE
C           BOTTOM TEMPERATURE BOUNDARY IS SOIL SURFACE
            TSPDT(NSP+1) = TSDT(1)
         END IF
      END IF
C
      END IF
C
      RUNOFF=RAIN
C
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE WTBULB (TW,TA,HUM)
C
C     THIS SUBROUTINE CALCULATES THE WETBULB TEMPERATURE (TW IN CELCIUS)
C     FROM THE AIR TEMPERATURE AND RELATIVE HUMIDITY
C     ---- PROCEDURE FOR FINDING WETBULB TEMPERATURE WAS ADOPTED FROM
C          THE ANDERSON SNOWMELT MODEL.
C
C***********************************************************************
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      REAL LF,LS,LV
C
C     CONVERT ATMOSPHERIC PRESSURE FROM PASCALS TO MBARS
      PA=0.01*PRESUR
C
C     CALCULATE VAPOR PRESSURE AND SATURATED VAPOR PRESSURE IN MBARS
      EAS= 2.7489E8*EXP(-4278.63/(TA+242.792))
      EA=HUM*EAS
C
      DELTA= (4278.63/((TA+242.792)*(TA+242.792)))*EAS
      DO 10 I=1,3
         TW= DELTA*TA+6.6E-4*PA*TA+7.59E-7*PA*TA*TA+EA-EAS
         TW= TW/(DELTA+6.6E-4*PA+7.59E-7*PA*TA)
         TAV= (TA+TW)*0.5
         EAV= 2.7489E8*EXP(-4278.63/(TAV+242.792))
         DELTA= (4278.63/((TAV+242.792)*(TAV+242.792)))*EAV
   10 CONTINUE
C
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE NWSNOW (NSP,NC,NR,NS,NPLANT,ICESPT,ZSP,DZSP,RHOSP,
     >  TSPDT,DLW,DLWDT,TWBULB,SNOW,SNODEN,TQVSP,MELT,ZC,PCANDT,ZR,TRDT,
     >  GMCDT,RHOR,ZS,VLCDT,VICDT,TSDT,MATDT,CONCDT)
C
C     THIS SUBROUTINES ADDS SNOW LAYERS FOR NEWLY FALLEN SNOW
C
C***********************************************************************
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /SPWATR/ CLAG1,CLAG2,CLAG3,CLAG4,PLWMAX,PLWDEN,PLWHC,
     >                THICK,CTHICK
      REAL LF,LS,LV
C
      REAL ZSP(100),DZSP(100),RHOSP(100),TSPDT(100),DLW(100),DLWDT(100),
     >  TQVSP(100)
      REAL ZC(11),PCANDT(8)
      REAL ZR(10),GMCDT(10),TRDT(10),RHOR(10),CRES(10)
      REAL ZS(50),VLCDT(50),VICDT(50),TSDT(50),MATDT(50),CONCDT(10,50),
     >     CS(50)
      INTEGER ICESPT(100)
      REAL MELT
C
C     PRECIPITATION IS SNOW
      IF (SNODEN .GT. 0.0) THEN
C        DENSITY OF NEW SNOW IS KNOWN -- CONVERT TO KG/M3
         DENSTY = SNODEN*RHOL
       ELSE
C        COMPUTE DENSITY OF THE NEW SNOW BASED ON WETBULB TEMPERATURE
         IF (TWBULB .LE. -15.) THEN
            DENSTY = 50.
          ELSE
            DENSTY = 50. + 1.7*((TWBULB + 15.)**1.5)
         END IF
      END IF
      DEPTH = SNOW*RHOL/DENSTY
C
      IF (NC .GT. 0) THEN
C        ADJUST CANOPY INTERCEPTION FOR DEPTH OF CANOPY COVERED BY SNOW
C        (NO CONSIDERATION ABOUT WHERE LAI IS WITHIN THE CANOPY)
         TOTAL=0.0
         DO 3 J=1,NPLANT
            CHANGE=PCANDT(J)*DEPTH/ZC(NC+1)
            IF (CHANGE .GT. PCANDT(J)) CHANGE=PCANDT(J)
            PCANDT(J)=PCANDT(J)-CHANGE
            TOTAL=TOTAL+CHANGE
    3    CONTINUE
         SNOW=SNOW+TOTAL
         DENSTY=SNOW*RHOL/DEPTH
      END IF
C
      MELT = 0.0
      IF (NSP .EQ. 0) THEN
C        SNOWCOVER DOES NOT EXIST -- MELT SOME SNOW TO GET UNDERLYING
C        MATERIAL TO 0 DEGREES
         IF (NR .GT. 0) THEN
C           SNOW IS FALLING ON RESIDUE
            CALL RESHT (NR,CRES,GMCDT,RHOR)
            IF (TRDT(1) .LE. 0.0) THEN
C              NO MELTING OF SNOW IS NEEDED TO GET RESIDUE TO 0 DEGREES
               MELT=0.0
              ELSE
               MELT = CRES(1)*ZR(2)/2.*TRDT(1)/(RHOL*LF-RHOL*CI*TWBULB)
               SNOW = SNOW - MELT
               DEPTH = DEPTH - MELT*RHOL/DENSTY
               IF (SNOW .LE. 0.0) THEN
C                 ALL SNOW THAT FELL IS MELTED - ADJUST RESIDUE FOR
C                 ENERGY REQUIRED TO MELT SNOW (MELTWATER ENERGY WILL
C                 BE CONSIDERED LATER IN SUBROUTINE RAINRS)
                  MELT = MELT + SNOW
                  TRDT(1)= TRDT(1) - (RHOL*LF-RHOL*CI*TWBULB)*MELT
     >                               /(CRES(1)*ZR(2)/2.)
                  RETURN
               END IF
               TRDT(1) = 0.0
            END IF
C
          ELSE
C           SNOW IS FALLING ON BARE SOIL
            IF (TSDT(1) .LE. 0.0) THEN
C              NO MELTING OF SNOW IS NEEDED TO GET SOIL TO 0 DEGREES
               MELT=0.0
              ELSE
               CALL SOILHT(NS,CS,VLCDT,VICDT,TSDT,MATDT,CONCDT)
               MELT = CS(1)*ZS(2)/2.*TSDT(1)/(RHOL*LF-RHOL*CI*TWBULB)
               SNOW = SNOW - MELT
               DEPTH = DEPTH - MELT*RHOL/DENSTY
               IF (SNOW .LE. 0.0) THEN
C                 ALL SNOW THAT FELL IS MELTED - ADJUST SOIL TEMP FOR
C                 ENERGY REQUIRED TO MELT SNOW (MELTWATER ENERGY WILL
C                 BE CONSIDERED LATER IN SUBROUTINE RAINSL)
                  MELT = MELT + SNOW
                  TSDT(1)= TSDT(1) - (RHOL*LF-RHOL*CI*TWBULB)*MELT
     >                               /(CS(1)*ZS(2)/2.)
                  RETURN
               END IF
               TSDT(1) = 0.0
            END IF
         END IF
C
       ELSE
C        NEW SNOW IS FALLING ON OLD SNOW - CHECK IF NEW SNOW CAUSES THE
C        THE SURFACE LAYER TO EXCEED THE SPECIFIED DEPTH LIMIT.
         IF((DZSP(1) + DEPTH) .LE. (1.55*THICK)) THEN
C           INCORPORATE ALL OF NEW SNOW INTO OLD SURFACE LAYER.
            DNS = DEPTH
            DEPTH = 0.0
            NWLAYR = 0
           ELSE
C           DETERMINE HOW MUCH SHOULD BE COMBINED WITH OLD SURFACE LAYER
            IF (DZSP(1) .GT. THICK) THEN
C              LAYER IS ALREADY GREATER THAN DESIRED THICKNESS - DO NOT
C              ADD NEW SNOW TO OLD SURFACE LAYER
               DNS = 0.0
C              GO SET UP NEW LAYERS FOR THE NEW SNOW
               GO TO 5
            END IF
            DNS = THICK - DZSP(1)
            DEPTH = DEPTH - DNS
         END IF
C
C        INCORPORATE NEW SNOW INTO OLD SURFACE LAYER - DEFINE PROPERTIES
         WEI=RHOSP(1)*DZSP(1)/RHOL
         WENS=DENSTY*DNS/RHOL
         WEL=WEI+WENS
         DZSP(1) = DZSP(1) + DNS
         RHOSP(1)=RHOL*WEL/DZSP(1)
         TSPDT(1)=(WENS*TWBULB + WEI*TSPDT(1))/WEL
         IF (ICESPT(1) .EQ. 1)  THEN
C           LIQUID-WATER AND SUBFREEZING TEMP MAY EXIST IN NEW LAYER
            FREEZE=-(TSPDT(1)*CI*WEL)/(LF*RHOL)
            IF (DLWDT(1) .LE. FREEZE) THEN
               TSPDT(1)=TSPDT(1)+(DLWDT(1)*LF*RHOL)/(CI*WEL)
               WEL=WEL+DLWDT(1)
               RHOSP(1)=RHOL*WEL/DZSP(1)
               DLWDT(1)=0.0
               ICESPT(1)=0
              ELSE
               DLWDT(1)=DLWDT(1)-FREEZE
               WEL=WEL+FREEZE
               RHOSP(1)=RHOL*WEL/DZSP(1)
               TSPDT(1)= 0.0
            END IF
         END IF
      END IF
C
C
C     COMPUTE THE NUMBER OF LAYERS WITH ENTIRELY NEW SNOW
    5 NWLAYR = DEPTH/THICK
      EXTRA = DEPTH - NWLAYR*THICK
      IF (EXTRA .GT. 0.55*THICK) NWLAYR = NWLAYR + 1
      IF (EXTRA .GT. 0.0  .AND.  NWLAYR .EQ. 0) NWLAYR = 1
C
      IF (NSP .GT. 0) THEN
         DOWN = DEPTH + DNS
         IF (NWLAYR .EQ. 0) THEN
C           ADJUST DEPTHS OF THE LAYERS AND RETURN
            ZSP(1) = 0.0
            DO 10 I=2,NSP+1
               ZSP(I) = ZSP(I) + DOWN
   10       CONTINUE
            RETURN
          ELSE
C           MOVE LAYERS DOWN TO MAKE ROOM FOR NEW SNOW LAYERS
            ZSP(NSP+NWLAYR+1) = ZSP(NSP+1) + DOWN
            TSPDT(NSP+NWLAYR+1) = TSPDT(NSP+1)
            DO 15 I=NSP,1,-1
               RHOSP(NWLAYR + I)= RHOSP(I)
               TSPDT(NWLAYR + I)= TSPDT(I)
               DLW(NWLAYR + I)= DLW(I)
               DLWDT(NWLAYR + I)= DLWDT(I)
               ICESPT(NWLAYR + I)= ICESPT(I)
               ZSP(NWLAYR + I) = ZSP(I) + DOWN
               DZSP(NWLAYR + I) = DZSP(I)
               TQVSP(NWLAYR + I) = TQVSP(I)
   15       CONTINUE
C           ADJUST DEPTH OF OLD SURFACE LAYER TO BE IN MIDDLE OF LAYER
            ZSP(NWLAYR+1)=ZSP(NWLAYR+1)+DZSP(NWLAYR+1)/2.
         END IF
      END IF
C
C     LEFT-OVER SNOW GOES IN TOP LAYER
      EXTRA=DEPTH - (NWLAYR-1)*THICK
      ZSP(1)=0.0
      DZSP(1)=EXTRA
      RHOSP(1)=DENSTY
      TSPDT(1)=TWBULB
      DLW(1)=0.0
      DLWDT(1)=0.0
      ICESPT(1)= 0
      TQVSP(1)=0.0
      TDEPTH=DZSP(1)
C
C     FULL SIZE NEW SNOW LAYERS.
      DO 20 J=2,NWLAYR
         DZSP(J)=THICK
         ZSP(J)= TDEPTH + DZSP(J)/2.
         TDEPTH=TDEPTH + DZSP(J)
         RHOSP(J)=DENSTY
         TSPDT(J)=TWBULB
         DLW(J)=0.0
         DLWDT(J)=0.0
         ICESPT(J)= 0
         TQVSP(J)=0.0
   20 CONTINUE
C
      IF (NSP .EQ. 0) THEN
C        DEFINE DEPTH OF BOTTOM BOUNDARY OF SNOWPACK
         ZSP(NWLAYR+1) = TDEPTH
C        NEW SNOW COVER - DEFINE TEMP OF BOTTOM BOUNDARY OF SNOWPACK
C        THIS NODE IS SHARED BY THE MATERIAL BELOW
         IF (NR .GT. 0) THEN
C           UNDERLYING MATERIAL IS RESIDUE
            TSPDT(NWLAYR+1) = TRDT(1)
           ELSE
C           UNDERLYING MATERIAL IS SOIL
            TSPDT(NWLAYR+1) = TSDT(1)
         END IF
      END IF
C
      NSP = NSP + NWLAYR
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE WBSNOW (NSP,ICESPT,ZSP,DZSP,RHOSP,TSPDT,DLW,DLWDT,RAIN,
     >           TRAIN,TOPSNO,EVAP1,TQVSP,WLAG,STORE,SCOUT)
C
C     THIS SUBROUTINE PERFORMS A WATER BALANCE OF THE SNOWPACK BY
C     ADSUSTING THE DENSITY FOR VAPOR FLUX AND ANY MELT WHICH OCCURRED
C     OVER THE TIME STEP.  CHECKS ARE MADE TO SEE IF ANY LAYERS HAVE
C     DISAPPEARED DUE TO MELT, OR IF ANY LAYERS ARE OUTSIDE THE
C     ACCEPTABLE THICKNESS RANGE.
C
C***********************************************************************
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /SPWATR/ CLAG1,CLAG2,CLAG3,CLAG4,PLWMAX,PLWDEN,PLWHC,
     >                THICK,CTHICK
      REAL LF,LS,LV
C
      REAL RHOSP(100),DZSP(100),ZSP(100),TSPDT(100),DLW(100),DLWDT(100),
     >  TQVSP(100),WLAG(11)
      INTEGER ICESPT(100)
C
      SCOUT = 0.0
C
C     ADJUST THE DENSITY FOR VAPOR FLUX
C     (THICKNESS OF LAYER IS ADJUSTED FOR SURFACE LAYER)
C
      CHANGE = (TOPSNO - TQVSP(1))/RHOSP(1)
      DZSP(1) = DZSP(1) + CHANGE
      IF (DZSP(1).LE.0.0 .AND. NSP.EQ.1) THEN
C        SNOWPACK IS LOST TO SUBLIMATION -- ADJUST TOPSNO AND EVAP1 FOR
C        THE WATER BALANCE TO CLOSE AND RETURN TO CALLING SUBROUTINE
         TOPSNO=TOPSNO - DZSP(1)*RHOSP(1)         
         EVAP1=EVAP1 - DZSP(1)*RHOSP(1)/RHOL
         NSP=0
         RETURN
      END IF
      DO 10 I=2,NSP
         CHANGE = (TQVSP(I-1) - TQVSP(I))/DZSP(I)
         RHOSP(I) = RHOSP(I) + CHANGE
   10 CONTINUE
C
C     ADJUST LAYERS FOR MELT AND RAINFALL
C        ABOVE = THE RAIN OR THE MELTWATER FROM ABOVE LAYERS
C        EXCESS= THIS IS THE MELTWATER THAT WOULD HAVE BEEN PRODUCED
C                BY ABOVE LAYERS HAD THERE BEEN ENOUGH ICE IN THE LAYER
C                FOR ALL THE ENERGY ABSORBED.  (FOR THE FIRST LAYER,
C                THIS TERM INCLUDES THE ENERGY TRANSFERRED BY RAIN.)
      NL=NSP
      ABOVE= RAIN
      EXCESS= RAIN*(CL*TRAIN+LF)/LF
C
      DO 20 I=1,NSP
      IF (ICESPT(I) .EQ. 0) THEN
C
C        TEMPERATURE IS UNKNOWN--SOME EXCESS LIQUID-WATER FROM
C        THE ABOVE LAYER MUST BE FROZEN.
C
         IF (EXCESS .NE. 0.0) THEN
C           COMPUTE AMOUNT TO BE FROZE IN ORDER TO RAISE THE TEMPERATURE
C           TO ZERO DEGREES CELSIUS. (USE SPECIFIC HEAT OF ICE, CI)
            WEL=RHOSP(I)*DZSP(I)/RHOL
            FREEZE=-TSPDT(I)*WEL*CI/LF
            IF (EXCESS .GT. FREEZE) THEN
C              EXCESS EXCEEDS REFREEZE.
               TSPDT(I)=0.0
               WEL=WEL+FREEZE
               RHOSP(I)=RHOL*WEL/DZSP(I)
               EXCESS=EXCESS-FREEZE
               ABOVE=ABOVE-FREEZE
               ICESPT(I) = 1
             ELSE
C              EXCESS IS ALL FROZEN IN THIS LAYER.
               TSPDT(I)=TSPDT(I) +(EXCESS*LF*RHOL)/(CI*RHOSP(I)*DZSP(I))
               WEL=WEL+EXCESS
               RHOSP(I)=RHOL*WEL/DZSP(I)
               ABOVE = ABOVE - EXCESS
               EXCESS=0.0
            END IF
         END IF
      END IF
C
      DLWDT(I)= DLWDT(I) + EXCESS
      CHANGE= DLWDT(I) - DLW(I) - ABOVE
C
      IF (ABS(CHANGE).LT.0.0000001) THEN
C        THIS IS BEYOND PRECISION OF COMPUTER -- ASSUME ZERO
         ABOVE=0.0
         EXCESS=0.0
       ELSE
         IF (CHANGE .LE. 0.0) THEN
C           SOME LIQUID-WATER IS FROZEN. ADD TO ICE CONTENT OF LAYER.
            WEL=RHOSP(I)*DZSP(I)/RHOL
            WEL=WEL-CHANGE
            RHOSP(I)=RHOL*WEL/DZSP(I)
            ABOVE=0.0
            EXCESS=0.0
          ELSE
C           MELT HAS OCCURRED,SUBTRACT FROM ICE CONTENT
            WEL=RHOSP(I)*DZSP(I)/RHOL
C           IF MELT EXCEEDS ICE CONTENT OF LAYER--LAYER IS GONE.
            IF (CHANGE .GE. WEL) THEN
C              LAYER IS GONE
C              LIQUID-WATER IS ADDED TO THE LAYER BELOW (EXCESS).  THE 
C              ICE CONTENT PLUS DLW OF THE LAYER CANNOT BE TAKEN FROM THE
C              NEXT LAYER, BUT IS STILL ADDED TO THE LIQUID-WATER 
C              CONTENT OF THE NEXT LAYER (ABOVE).
               NL=NL-1
               ABOVE= ABOVE + WEL + DLW(I)
               EXCESS=DLWDT(I)
               DZSP(I)=0.0
             ELSE
C              LAYER REMAINS
               WEL=WEL-CHANGE
               DZSP(I)=RHOL*WEL/RHOSP(I)
               ABOVE=0.0
               EXCESS=0.0
            END IF
         END IF
      END IF
C
   20 CONTINUE
C
C***********************************************************************
C     CHECK TO SEE IF ENTIRE SNOW COVER IS GONE.
      IF (NL .LE. 0) THEN
         NSP = 0
         WE = 0.0
         TDEPTH = 0.0
         SCOUT = ABOVE
         CALL SNOMLT (NSP,ICESPT,DZSP,RHOSP,TSPDT,DLWDT,WLAG,
     >                STORE,SCOUT)
         RETURN
      END IF
C
C     ELIMINATE LAYERS WHICH ARE GONE.
      IF (NL .NE. NSP) THEN
         DO 30 I=1,NSP
   23       IF (DZSP(I) .GT. 0.0) GO TO 30
C           LAYER GONE. MOVE OTHER LAYERS UP.
            NEXT=I+1
            DO 25 J=NEXT,NSP
               L=J-1
               DZSP(L)=DZSP(J)
               RHOSP(L)=RHOSP(J)
               TSPDT(L)=TSPDT(J)
               DLWDT(L)=DLWDT(J)
               ICESPT(L)=ICESPT(J)
   25       CONTINUE
            NSP = NSP - 1
            IF (NSP .EQ. NL) GO TO 35
            GO TO 23
   30    CONTINUE
   35    DLWDT(NSP)= DLWDT(NSP) + ABOVE
      END IF
C
C***********************************************************************
C     CHECK THICKNESS OF EACH LAYER. IF NOT WITHIN SPECIFIED LIMITS,
C     DIVIDE (IF TOO LARGE) OR ADD TO AN ADJACENT LAYER (IF TOO SMALL).
      TDEPTH=0.0
      I=1
      IF (NSP .EQ. 1) GO TO 80
C
   40 ZZ=TDEPTH+DZSP(I)*0.5
C     COMPUTED DESIRED THICKNESS FOR THE LAYER.
      DZ1=THICK
      IF (ZZ .GT. 0.30)  DZ1=CTHICK*(ZZ-0.30)+THICK
C     CHECK ACTUAL THICKNESS AGAINST DESIRED THICKNESS.
      IF (DZSP(I) .GT. 1.55*DZ1) GO TO 50
      IF (DZSP(I) .LT. 0.55*DZ1 .AND. NSP.GT.1) GO TO 60
C     THICKNESS IS WITHIN SPECIFIED LIMITS.
      TDEPTH=TDEPTH+DZSP(I)
      GO TO 70
C***********************************************************************
C     THICKNESS IS GREATER THAN SPECIFIED LIMTS.
C          SUB-DIVIDE LAYER,THUS CREATING A NEW LAYER.
C          PROPERTIES ARE THE SAME FOR BOTH LAYERS.
C          MOVE OTHER LAYERS DOWN.
   50 NEXT=I+1
      IF (NEXT .LE. NSP) THEN
         DO 55 J=NEXT,NSP
            L=NSP-J+NEXT
            LL=L+1
            DZSP(LL)=DZSP(L)
            RHOSP(LL)=RHOSP(L)
            TSPDT(LL)=TSPDT(L)
            DLWDT(LL)=DLWDT(L)
            ICESPT(LL)=ICESPT(L)
   55    CONTINUE
      END IF
      TDEPTH=TDEPTH+DZSP(I)
      DZSP(NEXT)=DZSP(I)*0.5
      RHOSP(NEXT)=RHOSP(I)
      TSPDT(NEXT)=TSPDT(I)
      DLWDT(NEXT)=DLWDT(I)*0.5
      ICESPT(NEXT)=ICESPT(I)
      DZSP(I)=DZSP(I)*0.5
      DLWDT(I)=DLWDT(I)*0.5
      I=I+1
      NSP=NSP+1
      GO TO 70
C***********************************************************************
C     THICKNESS IS SMALLER THAN SPECIFIED LIMITS.
C          ADD TO THE SMALLEST ADJACENT LAYER, THUS
C          LOSING A LAYER. PROPERIES OF THE NEW LAYER
C          ARE THE WEIGHTED AVERAGE OF THE TWO FORMER
C          LAYERS. MOVE OTHER LAYERS UP.
   60 IF (I .EQ. 1  .OR.  I .EQ. NSP) THEN
         IF (I .EQ. 1) THEN
            NL=2
           ELSE
            NL=NSP-1
         END IF
       ELSE
         NL=I-1
         IF (DZSP(I+1) .LT. DZSP(I-1)) NL=I+1
      END IF
C
C     NL IS THE SMALLEST ADJACENT LAYER - ADD LAYER I TO LAYER NL
      WEI=RHOSP(I)*DZSP(I)/RHOL
      WENL=RHOSP(NL)*DZSP(NL)/RHOL
      WEL=WEI+WENL
      DZSP(NL)=DZSP(NL)+DZSP(I)
      RHOSP(NL)=RHOL*WEL/DZSP(NL)
      DLWDT(NL)= DLWDT(NL) + DLWDT(I)
      TSPDT(NL)=(WENL*TSPDT(NL)+WEI*TSPDT(I))/WEL
C
      IF (ICESPT(I) .NE. ICESPT(NL)) THEN
C        UNKNOWNS ARE DIFFERENT. COMPUTE THE UNKNOWN FOR THE NEW LAYER.
         FREEZE=-(TSPDT(NL)*CI*RHOSP(NL)*DZSP(NL))/(LF*RHOL)
         IF (DLWDT(NL) .LE. FREEZE) THEN
C           TEMPERATURE IS UNKNOWN
            TSPDT(NL)=TSPDT(NL) + (DLWDT(NL)*LF*RHOL)
     >                            /(CI*RHOSP(NL)*DZSP(NL))
            WEL=WEL+DLWDT(NL)
            RHOSP(NL)=RHOL*WEL/DZSP(NL)
            ICESPT(NL)=0
          ELSE
C           LIQUID-WATER IS UNKNOWN.
            DLWDT(NL)=DLWDT(NL)-FREEZE
            WEL=WEL+FREEZE
            RHOSP(NL)=RHOL*WEL/DZSP(NL)
            ICESPT(NL)=1
         END IF
      END IF
C
      IF (ICESPT(NL) .EQ. 1) TSPDT(NL)=0.0
      IF (ICESPT(NL) .EQ. 0) DLWDT(NL)=0.0
C     MOVE OTHER LAYERS UP.
      IF (NL .LT. I) TDEPTH=TDEPTH+DZSP(I)
      NEXT=I+1
C
      IF (NEXT .LE. NSP) THEN
         DO 65 J=NEXT,NSP
            L=J-1
            DZSP(L)=DZSP(J)
            RHOSP(L)=RHOSP(J)
            TSPDT(L)=TSPDT(J)
            DLWDT(L)=DLWDT(J)
            ICESPT(L)=ICESPT(J)
   65    CONTINUE
      END IF
C
      I = I-1
      NSP = NSP-1
C
C***********************************************************************
C     CHECK FOR LAST LAYER.
   70 IF (I .LT. NSP) THEN
C        START NEXT LAYER
         I=I+1
         GO TO 40
      END IF
C
C***********************************************************************
C     ADJUST THE DENSITY FOR METAMORPHISM AND COMPUTE ANY SNOMELT
   80 CALL META (NSP,TSPDT,RHOSP,DZSP,DLWDT)
      CALL SNOMLT (NSP,ICESPT,DZSP,RHOSP,TSPDT,DLWDT,WLAG,STORE,SCOUT)
C
C***********************************************************************
C     CALCULATE THE DEPTH FOR EACH NODE
      ZSP(1) = 0.0
      TDEPTH = DZSP(1)
      WE = DZSP(1)*RHOSP(1)/RHOL + DLWDT(1)
      DO 85 I=2,NSP
         ZSP(I)= TDEPTH + DZSP(I)/2.
         TDEPTH= TDEPTH + DZSP(I)
         WE = WE + DZSP(I)*RHOSP(I)/RHOL + DLWDT(I)
   85 CONTINUE
      ZSP(NSP+1) = TDEPTH
C
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE META (NSP,TSPDT,RHOSP,DZSP,DLWDT)
C
C     COMPUTES THE CHANGE IN DENSITY OF THE SNOW COVER CAUSED BY
C     DESTRUCTIVE (EQUI-TEMPERATURE) METAMORPHISM, COMPACTION, AND THE
C     PRESENCE OF LIQUID-WATER.
C***********************************************************************
      COMMON /TIMEWT/ WT,WDT,DT
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /METASP/ CMET1,CMET2,CMET3,CMET4,CMET5,SNOMAX
      REAL TSPDT(100),RHOSP(100),DZSP(100),DLWDT(100)
      REAL LF,LS,LV
C
      IF (CMET1 .LE. 0.0  .OR.  CMET3 .LE. 0.0) THEN
         IF (CMET5 .LE. 0.0) RETURN
      END IF
C
C     WEIGHT IS THE WATER-EQUIVALENT (CM) ABOVE THE LAYER.
      WEIGHT=0.0
C
      DO 10 I=1,NSP
C        IF DENSITY IS THAT OF ICE, DO NOT INCREASE IT ANY MORE
         IF (RHOSP(I) .GE. RHOI) GO TO 10
         WEL= RHOSP(I)*DZSP(I)/RHOL
C
C****    DESTRUCTIVE METAMORPHISM TERM.
C
         TERM=EXP(CMET4*TSPDT(I))
         IF (RHOSP(I) .LE. SNOMAX) THEN
            T1=TERM*CMET3
           ELSE
            T1=TERM*CMET3*EXP(-46.0*(RHOSP(I)-SNOMAX)/RHOL)
         END IF
C
C****    COMPACTION TERM.
C
         TERM=EXP(0.08*TSPDT(I))
         T2=WEIGHT*CMET1*EXP(-CMET2*RHOSP(I)/RHOL)*TERM
C
C***     LIQUID-WATER TERM.
C
         IF (DLWDT(I) .GT. 0.0) T1=CMET5*T1
C
C****    DENSIFICATION OF THE LAYER.  (WATER-EQUIVALENT STAYS THE SAME.)
C
         RHOSP(I)=RHOSP(I)*(1.0 + DT*(T1+T2)/3600.)
         DZSP(I)=RHOL*WEL/RHOSP(I)
         WEIGHT=WEIGHT + (WEL+DLWDT(I))*100.
   10 CONTINUE
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE SNOMLT (NSP,ICESPT,DZSP,RHOSP,TSPDT,DLWDT,WLAG,
     >                   STORE,SCOUT)
C
C     THIS SUBROUTINE DETERMINES THE SNOW COVER OUTFLOW DURING EACH 
C     PERIOD BASED ON THE CHANGE IN LIQUID-WATER DURING THE PERIOD AND
C     THE PHYSICAL CHARTERISTICS OF THE SNOW COVER.
C***********************************************************************
      COMMON /TIMEWT/ WT,WDT,DT
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /SPWATR/ CLAG1,CLAG2,CLAG3,CLAG4,PLWMAX,PLWDEN,PLWHC,
     >                THICK,CTHICK
      REAL LF,LS,LV
C
      REAL DZSP(100),RHOSP(100),TSPDT(100),DLWDT(100),WLAG(11)
      INTEGER ICESPT(100)
      SAVE IFIRST,NLAG
      DATA IFIRST/0/
C
      IF(IFIRST .EQ. 0) THEN
C        INITIALIZE VARIABLES.
C        1. STORE IS THE AMOUNT OF LIQUID-WATER(THAT HAS ALREADY BEEN
C           LAGGED) THAT IS IN STORAGE IN THE SNOW-COVER.
C        2. WLAG() IS THE AMOUNT OF LIQUID-WATER IN THE PROCESS OF
C           BEING LAGGED.  NLAG IS THE NUMBER OF ARRAY ELEMENTS USED.
         JDT=CLAG1+0.01
         NLAG=JDT+2
         IF (NLAG.GT.11) NLAG=11
         IFIRST=1
      END IF
C
      IF (NSP .LE. 0) THEN
C        SNOW COVER HAS JUST DISAPPEARED
         DO 10 J=1,NLAG
            SCOUT=SCOUT+WLAG(J)
            WLAG(J)=0.0
   10    CONTINUE
         SCOUT=SCOUT+STORE
         STORE=0.0
         RETURN
      END IF
C
C     DETERMINE THE EXCESS LIQUID-WATER (IN EXCESS OF LIQUID-WATER
C     HOLDING CAPACITY) GENERATED DURING THIS TIME PERIOD.
      EXCESS=0.0
C
C     EXCESS WATER IN BOTTOM LAYER IS NOT TO BE ROUTED.
      IF (NSP .EQ. 1) GO TO 12
      BOTTOM=0.0
      IF (ICESPT(NSP) .EQ. 1) THEN
         WEL=RHOSP(NSP)*DZSP(NSP)/RHOL
         PLW=(PLWMAX-PLWHC)*(PLWDEN-RHOSP(NSP))/PLWDEN+PLWHC
         IF (RHOSP(NSP) .GE. PLWDEN) PLW=PLWHC
         WMAX=PLW*WEL
         IF (DLWDT(NSP) .GT. WMAX) BOTTOM= DLWDT(NSP) - WMAX
      END IF
C
   12 WE=0.0
      TDEPTH=0.0
      DO 20 I=1,NSP
         WEL=RHOSP(I)*DZSP(I)/RHOL
         IF(ICESPT(I) .EQ. 0) THEN
C           LAYER IS BELOW ZERO DEGREES CELSIUS.
            IF(EXCESS .EQ. 0.0) GO TO 15
C           FREEZE SOME OF THE LIQUID-WATER.
            FREEZE=-TSPDT(I)*CI*WEL/LF
            IF (EXCESS .LE. FREEZE) THEN
C              EXCESS IS ALL FROZEN IN THIS LAYER.
               TSPDT(I)=TSPDT(I)+(EXCESS*LF*RHOL)/(CI*WEL*RHOL)
               WEL=WEL+EXCESS
               EXCESS=0.0
               RHOSP(I)=RHOL*WEL/DZSP(I)
               GO TO 15
              ELSE
C              EXCESS EXCEEDS REFREEZE.
               TSPDT(I)= 0.0
               WEL=WEL+FREEZE
               RHOSP(I)=RHOL*WEL/DZSP(I)
               EXCESS=EXCESS-FREEZE
               ICESPT(I)=1
            END IF
         END IF
         PLW=(PLWMAX-PLWHC)*(PLWDEN-RHOSP(I))/PLWDEN+PLWHC
         IF (RHOSP(I) .GE. PLWDEN) PLW=PLWHC
         WMAX= PLW*WEL
         W= DLWDT(I) + EXCESS
         IF (W .LE. WMAX) THEN
C           LIQUID-WATER HOLDING CAPACITY IS NOT SATISFIED.
            DLWDT(I)=W
            EXCESS=0.0
           ELSE
C           LIQUID-WATER HOLDING CAPACITY IS EXCEEDED.
            DLWDT(I)=WMAX
            EXCESS=W-WMAX
         END IF
   15    WE=WE+WEL
         TDEPTH=TDEPTH+DZSP(I)
   20 CONTINUE
C
      IF (NSP .EQ. 1) THEN
C        WATER NOT LAGGED IF ONLY ONE NODE IN THE SNOWPACK
         SCOUT = EXCESS
         DO 50 J=1,NLAG
            SCOUT=SCOUT+WLAG(J)
            WLAG(J)=0.0
   50    CONTINUE
         SCOUT=SCOUT+STORE
         STORE=0.0
         RETURN
      END IF
C
      EXCESS=EXCESS-BOTTOM
      IF (ICESPT(NSP) .EQ. 0) THEN
C        DO NOT ALLOW ANY WATER LAGGED OR STORED TO LEAVE SNOWPACK
C        (EXCESS MUST BE EQUAL 0.0)
         OUTHR=0.0
         SCOUT=0.0
         RETURN
      END IF
C
C     ROUTE EXCESS WATER THROUGH THE SNOW COVER.
C     EMPIRICAL LAG AND ATTENUATION EQUATIONS - ONE HOUR TIME STEP USED.
      EXCSHR = EXCESS*3600./DT
      IDT=DT/3600.
      IF (AMOD(DT,3600.).GT.0.00001) IDT=IDT+1
      IF (IDT .EQ. 0) IDT=1
      HOURS=1.0
C
      DENSE = WE/TDEPTH
      SCOUT=0.0
      DO 35 IHR=1,IDT
C        IF THIS IS THE LAST TIME INCREMENT, ACCOUNT FOR ANY REMAINING 
C        FRACTION OF AN HOUR IN THE TIME STEP
         IF (IHR.EQ.IDT) HOURS=DT/3600-(IDT-1)
         EXCESS=EXCSHR*HOURS
C
         OUTHR=0.0
C        LAG-FUNCTION OF DEPTH,DENSITY,AND EXCESS WATER.
         IF (EXCSHR .GE. 0.00001) THEN
            NI= ((EXCSHR*10000.0)**0.3)+0.5
            IF (NI .LT. 1) NI=1
            FN=NI
            FLMAX=CLAG1*(1.0-EXP(-0.25*TDEPTH/DENSE))
            DO 25 J=1,NI
               FJ=J
               FLAG=FLMAX/(CLAG2*EXCSHR*100.*(FJ-0.5)/FN+1.0)
               K= FLAG+1.0
               POR=K-FLAG
               WINC=1.0/FN
               WLAG(K)=WLAG(K)+EXCESS*WINC*POR
               WLAG(K+1)=WLAG(K+1)+EXCESS*WINC*(1.0-POR)
   25       CONTINUE
           ELSE
            WLAG(1)=WLAG(1)+EXCESS
         END IF
C
C        ATTENUATION-FUNCTION OF DENSITY AND PREVIOUS OUTFLOW.
         IF ((STORE+WLAG(1)) .NE. 0.0) THEN
            R=1.0/(CLAG3*EXP(-CLAG4*WLAG(1)*DENSE/TDEPTH)+1.0)
            OUTHR=(STORE+WLAG(1))*R
            STORE=STORE+(WLAG(1)-OUTHR)*HOURS
            SCOUT=SCOUT+OUTHR*HOURS
            IF(STORE .LE. 0.00001) THEN
               OUTHR=OUTHR+STORE
               SCOUT=SCOUT+STORE
               STORE=0.0
            END IF
         END IF
         NI=NLAG-1
         DO 30 J=1,NI
            WLAG(J)=WLAG(J)*(1.-HOURS) + WLAG(J+1)*HOURS
            IF (WLAG(J+1).LE. 0.000001) THEN
C              TIME STEPS OF LESS THAN 1 HOUR RESULT IN CONTINUALLY
C              TAKING FRACTIONS OF THE LAGGED DEPTH; DEPTH IS TOO
C              SMALL -- PASS THE ENTIRE DEPTH ONTO THE NEXT LAG
               WLAG(J)=WLAG(J)+WLAG(J+1) + WLAG(J+1)*(1.-HOURS)
               WLAG(J+1)=0.0
            END IF
   30    CONTINUE
         WLAG(NLAG)=0.0
   35 CONTINUE
C
      SCOUT=SCOUT+BOTTOM
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE RAINC (NPLANT,NC,LANGLE,ITYPE,RAIN,PCANDT,WCANDT,WCMAX,
     >                  SLOPE)
C
C     THIS SUBROUTINE CALCULATES THE RAINFALL DEPTH INTERCEPTED BY THE
C     CANOPY LAYERS AND ADJUST THE WATER CONTENT OF THE DEAD PLANT
C
C***********************************************************************
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /CLAYRS/ DRYCAN(8,10),CANLAI(8,10),TOTLAI(8),IEVAP(8),
     >                RLEAF(8,10),RROOT(8,50),ROOTDN(8,50),TOTROT(8)
C
      REAL PCANDT(8),WCANDT(10)
      REAL TDIRCC(10),TDIFFC(10),DIRKL(9,10),DIFKL(9,10)
      REAL LF,LS,LV
      REAL MAXINT
      INTEGER LANGLE(8),ITYPE(8)
C
C     USE THE TRANSMITTANCE TO DIRECT RADIATION TO CALCULATE THE
C     FRACTION OF RAIN INTERCEPTED BY THE CANOPY
      ANGLE=1.5708-SLOPE
      CALL TRANSC (NPLANT,NC,LANGLE,TDIRCC,TDIFFC,DIFKL,DIRKL,
     >             ANGLE,1.5708)
C
C**** CALCULATE RAINFALL INTERCEPTED BY THE ENTIRE CANOPY
      TRANS=1.0
      DO 5 I=1,NC
         TRANS=TRANS*TDIRCC(I)
    5 CONTINUE
      RINTER=RAIN*(1.-TRANS)
C
      SUMLAI=0.0
      DO 10 J=1,NPLANT
         SUMLAI=SUMLAI+TOTLAI(J)
   10 CONTINUE
C
      TOTAL = 0.0
C**** DETERMINE HOW MUCH RAINFALL INTERCEPTED BY EACH PLANT SPECIES
      DO 40 J=1,NPLANT
         IF (ITYPE(J) .NE. 0) THEN
C****       PRECIP INTERCEPTED BY LIVING PLANTS
            SLAI=0.0
            DO 20 I=1,NC
               SLAI=SLAI+CANLAI(J,I)
   20       CONTINUE
C
C           CALCULATE PRECIP INTERCEPTED BY ALL LAYERS
            RINT=RINTER*SLAI/SUMLAI
C           ADD TO PRECIP REMAINING ON PLANT LEAVES FROM BEFORE
            PCANDT(J)=PCANDT(J)+RINT
C           ESTIMATE HOW MUCH PRECIP CAN BE HELD ON PLANTS
C           -- ASSUME LIVING PLANTS CAN INTERCEPT AND HOLD 1MM OF WATER
C           FOR EACH UNIT OF LEAF AREA INDEX
            MAXINT=0.001*SLAI
C           CHECK IF PRECIP INTERCEPTED EXCEEDS AMOUNT THAT CAN BE HELD
            IF (PCANDT(J) .GT. MAXINT) THEN
               RINT = MAXINT - (PCANDT(J) - RINT)
               PCANDT(J) = MAXINT
            END IF
            TOTAL = TOTAL + RINT
          ELSE
C****       DEAD PLANT MATERIAL
            EXCESS = 0.0
            DO 30 I=1,NC
C           CHECK IF PLANT EXISTS IN LAYER
            IF (CANLAI(J,I).GT.0.0) THEN
C              PRECIP INTERCEPTED BY LAYER
               RINTLY=RINTER*CANLAI(J,I)/SUMLAI + EXCESS
               WCANDT(I) = WCANDT(I) + RINTLY*RHOL/DRYCAN(J,I)
               EXCESS = 0.0
C              CHECK IF THE LAYER CAN HOLD THIS MUCH WATER
               IF (WCANDT(I) .GT. WCMAX) THEN
                  EXCESS = (WCANDT(I)-WCMAX)*DRYCAN(J,I)/RHOL
                  RINTLY = RINTLY-EXCESS
                  WCANDT(I)=WCMAX
               END IF
               TOTAL=TOTAL+RINTLY
            END IF
   30       CONTINUE
         END IF
   40 CONTINUE
C
C     CALCULATE THE AMOUNT OF PRECIP THAT MAKES IT THROUGH CANOPY
      RAIN = RAIN - TOTAL
C
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE RAINRS (NR,TRAIN,RAIN,ZR,TRDT,GMCDT,RHOR,
     >                   GMCMAX,DIRRES,SLOPE)
C
C     THIS SUBROUTINE CALCULATES THE RAINFALL DEPTH INTERCEPTED BY THE
C     RESIDUE LAYERS, AND ADJUSTS THE TEMPERATURE TO ACCOUNT FOR THE
C     HEAT (OR LACK THEREOF) INTRODUCED BY THE RAINFALL
C***********************************************************************
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      REAL LF,LS,LV
C
      REAL ZR(10),TRDT(10),GMCDT(10),RHOR(10),TDIREC(10),TDIFFU(10),
     >  CRES(10)
C
C     USE THE TRANSMITTANCE TO DIRECT RADIATION TO CALCULATE THE
C     FRACTION OF RAIN INTERCEPTED BY THE RESIDUE
      ANGLE=1.5708-SLOPE
      CALL TRANSR (NR,TDIREC,TDIFFU,ZR,RHOR,ANGLE,DIRRES)
C
C     CALCULATE THE SPECIFIC HEAT OF THE RESIDUE BEFORE RAIN IS ADDED
      CALL RESHT (NR,CRES,GMCDT,RHOR)
C
      DO 10 I=1,NR
C
         IF (I .EQ. 1) THEN
            DZ=ZR(2)-ZR(1)
            IF (NR .NE. 1) DZ=DZ/2.
           ELSE
            IF (I .EQ. NR) THEN
               DZ=ZR(I+1)-ZR(I)+(ZR(I)-ZR(I-1))/2.
              ELSE
               DZ=(ZR(I+1)-ZR(I-1))/2.
            END IF
         END IF
C
C
C****    CALCULATE RAINFALL INTERCEPTED BY THE NODE
         RINTER=RAIN*(1.-TDIREC(I))
         GMCDT(I)=GMCDT(I) + RINTER*RHOL/(DZ*RHOR(1))
         RAIN=RAIN-RINTER
C        CHECK IF THE NODE CAN HOLD THIS MUCH WATER
         IF (GMCDT(I) .GT. GMCMAX) THEN
            EXCESS= RHOR(I)*(GMCDT(I)-GMCMAX)*DZ/RHOL
            RINTER= RINTER-EXCESS
            RAIN= RAIN+EXCESS
            GMCDT(I)= GMCMAX
         END IF
C
C****    CALCULATE THE TEMPERATURE OF THE NODE
         RAINHT=RINTER*RHOL*CL
         RESIHT=CRES(I)*DZ
         TRDT(I)=TRDT(I) + RAINHT/(RAINHT+RESIHT)*(TRAIN-TRDT(I))
C
   10 CONTINUE
      END
C***********************************************************************
C
      SUBROUTINE RAINSL (NS,TRAIN,RAIN,ZS,TSDT,VLCDT,VICDT,ICESDT,
     >              MATDT,TOTFLO,SALTDT,CONCDT,POND,PONDMX)
C
C     THIS SUBROUTINE CALCULATES THE INFILTRATION INTO THE SOIL AS WELL
C     AS THE FINAL SOIL TEMPERATURE AND SOLUTE CONCENTRATION RESULTING
C     FROM THE HEAT AND SALTS CARRIED BY THE WATER.
C     ---- CURRENTLY THE SUBROUTINE WILL NOT ALLOW THE WATER CONTENT TO
C          EXCEED 90% OF THE SATURATED VALUE (PSAT=0.9)
C***********************************************************************
      COMMON /TIMEWT/ WT,WDT,DT
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /SLPARM/ B(50),ENTRY(50),RHOB(50),SAT(50),SATK(50),
     >                NSALT,SALTKQ(10,50),VAPCOF(50),VAPEXP(50),
     >                SAND(50),SILT(50),CLAY(50),OM(50)
      REAL LF,LS,LV
C
      REAL ZS(50),TSDT(50),VLCDT(50),VICDT(50),MATDT(50),
     >     TOTFLO(50),SALTDT(10,50),CONCDT(10,50),SATKMX(50),MATRIC(50),
     >     CS(50)
      INTEGER ICESDT(50)
      REAL INFMAX,INFIL,PSAT(50),VLC9(50),MAT9(50),SATK9(50)
C
      DATA PSATK/0.9/ PSAT/50*0.9/
C
C     DETERMINE THE MAXIMUM CONDUCTIVITY FOR EACH NODE -- THIS WILL BE
C     DEFINED BY THE AIR ENTRY POTENTIAL UNLESS THE SOIL CONTAINS ICE
      DO 10 I=1,NS
         IF (ICESDT(I) .GE. 1) THEN
C           SATURATED CONDUCTIVITY IS LIMITED BY THE PRESENCE OF ICE
            VLCMAX= SAT(I) - VICDT(I)
            IF (VLCDT(I) .GT. VLCMAX) VLCMAX=VLCDT(I)
            CALL MATVL1 (I,MATRIC(I),VLCMAX,DUMMY)
            VLC9(I)=PSATK*VLCMAX
           ELSE
C           SET MATRIC POTENTIAL FOR COMPUTING SATURATED CONDUCTIVITY 
C           AND INITIALIZE 90% (OR PSATK) OF MAXIMUM WATER CONENT
            MATRIC(I)=ENTRY(I)
            VLC9(I)=PSATK*SAT(I)
         END IF
         CALL MATVL1 (I,MAT9(I),VLC9(I),DUMMY)
   10 CONTINUE
      CALL SOILHK (NS,SATKMX,MATRIC,VLCDT,VICDT)
C
C     INITIALIZE INFILTRATION PARAMETERS
      DTIME=DT
      INFIL=0.0
      ZSTAR=0.0
      SUMZ=0.0
      SUMZK=0.0
      PSAT(1)=PSATK
C
      I = 1
C     CHECK IF CONDUCTIVITY OF SURFACE IS ZERO (POSSIBLY DUE TO ICE)
      IF (SATKMX(1) .LE. 0) GO TO 40
C     CHECK IF WETTING FRONT HAS REACHED SATURATED CONDITIONS
      IF (MATDT(1) .GE. 0.0) GO TO 40
C
C     CALCULATE THE MAXIMUM INFILTRATION THAT THE FIRST LAYER CAN HOLD
      INFMAX=(PSAT(1)*SAT(1) - VLCDT(1) - VICDT(1))*(ZS(2) - ZS(1))/2.
C
C     IF CURRENT LAYER IS ALREADY SATURATED -- GO TO NEXT LAYER
      IF (INFMAX .LE. 0.0) GO TO 30
C
C     CALCULATE THE POTENTIAL INFILTRATION FOR THE LAYER
   20 IF (ZSTAR .LE. 1.0) THEN
C        ASSUMPTION OF SATURATED FLOW BEHIND WETTING FRONT IS VALID
         DTSTAR= SATKMX(I)*DTIME/(PSAT(I)*SAT(I) - VLCDT(I) - VICDT(I))
     >                          /(-MATDT(I)+SUMZ)
         POTINF=(DTSTAR-2*ZSTAR + SQRT((DTSTAR-2*ZSTAR)**2 +8*DTSTAR))/2
         POTINF=POTINF*(PSAT(I)*SAT(I)-VLCDT(I)-VICDT(I))
     >                *(-MATDT(I)+SUMZ)
       ELSE
C        INFILTRATION FLOW HAS BECOME UNSATURATED -- USE FLOW THROUGH 
C        SATURATED MEDIUM 
         POTINF=SATINF*DTIME
      END IF
C
      IF (INFMAX .GE. (RAIN-INFIL)) THEN
C        LAYER CAN HOLD REMAINDER OF WATER - CHECK IF ALL CAN INFILTRATE
C
         IF (POTINF .GE. (RAIN-INFIL)) THEN
C           LAYER CAN INFILTRATE REMAINDER OF RAIN
            INFIL=RAIN
           ELSE
C           PROFILE CANNOT INFILTRATE ALL OF THE RAIN
            INFIL=INFIL+POTINF
         END IF
C
C        INFILTRATION CALCULATION IS COMPLETE - GO TO ENERGY CALCULATION
         GO TO 40
      END IF
C
C     LAYER CANNOT HOLD REMAINDER OF RAIN - CHECK IF IT CAN HOLD POTINF
C
      IF (POTINF .LT. INFMAX) THEN
C        PROFILE CANNOT INFILTRATE REMAINDER OF THE RAIN
         INFIL=INFIL+POTINF
C        INFILTRATION CALCULATION COMPLETE - GO TO ENERGY CALCULATION
         GO TO 40
      END IF
C
C     POTENTIAL INFILTRATION AND RAIN ARE SUFFICIENT TO FILL LAYER
C      --- CALCULATE TIME REQUIRED TO FILL THE LAYER AND ADJUST TIME
      IF (ZSTAR .LE. 1.0) THEN
C        ASSUMPTION OF SATURATED FLOW BEHIND WETTING FRONT IS VALID     
         DIMINF=INFMAX/(PSAT(I)*SAT(I)-VLCDT(I)-VICDT(I))
     >                /(-MATDT(I)+SUMZ)
         IF (DIMINF .GT. 0.01) THEN
            DTSTAR= (ZSTAR-1)*ALOG(1. + DIMINF) + DIMINF
           ELSE
C           IF DIMINF IS SMALL, USE ALTERNATE METHOD OF CALCULATING LOG
C           TERM BECAUSE SERIOUS ERRORS MAY RESULT WHEN TAKING THE LOG 
C           OF A NUMBER VERY CLOSE TO 1.0, I.E. (1. + 10^-5)
            DTSTAR= (ZSTAR-1)*2*DIMINF/(2+DIMINF) + DIMINF
         END IF
         TIME= DTSTAR*(PSAT(I)*SAT(I)-VLCDT(I)-VICDT(I))
     >               *(-MATDT(I)+SUMZ)/SATKMX(I)
       ELSE
C        INFILTRATION FLOW HAS BECOME UNSATURATED -- USE FLOW THROUGH
C        SATURATED MEDIUM
         TIME=INFMAX/SATINF
      END IF
C
      DTIME= DTIME-TIME
      INFIL=INFIL+INFMAX
      IF (DTIME .LE. 0.0) GO TO 40
C
C     UPDATE PARAMETERS FOR NEXT LAYER
   30 I = I + 1
C     CHECK IF HYDRAULIC CONDUCTIVITY IS ZERO (POSSIBLY DUE TO ICE)
      IF (SATKMX(I) .LE. 0) GO TO 40
C     CHECK IF WETTING FRONT HAS REACHED SATURATED CONDITIONS
      IF (MATDT(I) .GT. SUMZ) GO TO 40
C
      SUMZ2= (ZS(I)+ZS(I-1))/2.
      DZ= SUMZ2-SUMZ
      SUMZ= SUMZ2
      SUMZK= SUMZK + DZ/SATKMX(I-1)
C
      IF (I.LT.NS) THEN
         IF (ZSTAR .LE. 1.0) THEN
C           SATURATED FLOW STILL EXISTS -- SETUP ZSTAR FOR NEXT LAYER
            ZSTAR=SATKMX(I)*SUMZK/(-MATDT(I) + SUMZ)
            PSAT(I)=PSATK
            IF (ZSTAR .GT. 1.0) THEN
C              HYDRAULIC CONDUCTIVITY OF THIS LAYER IS SUFFICIENTLY HIGH
C              THAT SATURATED CONDITIONS NO LONGER EXIST -- CALCULATE 
C              INFILTRATION USING SATURATED RELATIONSHIPS BEHIND
C              CURRENT LAYER ASSUMING GRAVITY FLOW.
               SATINF=(-MATDT(I) + SUMZ)/SUMZK
C              COMPUTE CONDUCTIVITY OF ALL LAYERS AT PSATK (CURRENTLY 
C              SET AT 0.90) OF MAXIMUM WATER CONENT FOR INTERPOLATION
               CALL SOILHK (NS,SATK9,MAT9,VLC9,VICDT)
C              DETERMINE WATER CONTENT OF LAYER SUCH THAT CONDUCTIVITY 
C              MATCHES FLOW INTO LAYER - USE LOGARITHMIC INTERPOLATION
C              BETWEEN KNOWN POINTS OF VLC AND CONDUCTIVITY
               PSAT(I)=10.**(ALOG10(PSATK)*
     >             ALOG10(SATINF/SATKMX(I))/ALOG10(SATK9(I)/SATKMX(I)))
               IF (PSAT(I).GT. PSATK) PSAT(I)=PSATK
C              DEFINE MAXIMUM SATURATED CONDUCTIVITY FOR CURRENT LAYER
               SATKMX(I)=SUMZ/SUMZK
            END IF 
           ELSE
C           SATURATED CONDITIONS NO LONGER EXIST AT WETTING FRONT;
C           COMPARE FLUX TO SATURATED CONDUCTIVITY OF CURRENT LAYER
            SATINF=(-MATDT(I) + SUMZ)/SUMZK
            IF (SATKMX(I) .LT. SATINF) SATINF=SATKMX(I)
C           DETERMINE WATER CONTENT OF LAYER SUCH THAT CONDUCTIVITY 
C           MATCHES FLOW INTO LAYER - USE LOGARITHMIC INTERPOLATION
C           BETWEEN KNOWN POINTS OF VLC AND CONDUCTIVITY
            PSAT(I)=10.**(ALOG10(PSATK)*
     >          ALOG10(SATINF/SATKMX(I))/ALOG10(SATK9(I)/SATKMX(I)))
            IF (PSAT(I).GT. PSATK) PSAT(I)=PSATK
C           DEFINE MAXIMUM SATURATED CONDUCTIVITY FOR CURRENT LAYER
            SATKMX(I)=SUMZ/SUMZK
         END IF
C        CALCULATE THE MAXIMUM INFILTRATION THAT LAYER CAN HOLD
         INFMAX=(PSAT(I)*SAT(I)-VLCDT(I)-VICDT(I))*(ZS(I+1)-ZS(I-1))/2.
C        IF LAYER CANNOT HOLD ADDITIONAL WATER, GO TO NEXT LAYER
         IF (INFMAX .LT. 0.0) GO TO 30
         GO TO 20
        ELSE
C        INFILTRATION HAS REACHED BOTTOM OF PROFILE
         IF (ZSTAR .LE. 1.0) THEN
C           COMPUTE FLOW THROUGH SATURATED PROFILE      
            INFIL=INFIL + DTIME*SUMZ/SUMZK
           ELSE
C           COMPARE FLUX TO SATURATED CONDUCIVITY OF CURRENT LAYER
            IF (SATKMX(I) .LT. SATINF) SATINF=SATKMX(I)
            INFIL=INFIL + DTIME*SATINF
         END IF
         IF (INFIL .GT. RAIN) INFIL=RAIN
      END IF
C
C-----------------------------------------------------------------------
C     INFILTRATION CALCULATION COMPLETE -- CALCULATE WATER TO BE PONDED,
C     THEN DETERMINE TEMPERATURE, MOISTURE CONTENTS AND SOLUTES.
C
   40 POND= RAIN-INFIL
      RAIN= 0.0
      IF (POND .GT. PONDMX) THEN
C        SET RAIN EQUAL TO THE EXCESS RAIN, I.E. THE RUNOFF
         RAIN= POND-PONDMX
         POND= PONDMX
      END IF
C
      IF (INFIL.LE.0.0) THEN
C        UPDATE OF WATER CONTENTS, TEMPERATURES, AND SOLUTES COMPLETE
         INFMAX=0.0
         GO TO 110
      END IF
C
C     CALCULATE THE SPECIFIC HEAT OF THE LAYERS PRIOR TO ADDING WATER
      CALL SOILHT (NS,CS,VLCDT,VICDT,TSDT,MATDT,CONCDT)
C
C     DEFINE CONDITIONS FOR SURFACE NODE
C     SLOSS = SALTS LOST OR LEACHED FROM NODES ABOVE
      I = 1
      SLOSS = 0.0
      DZ=(ZS(2)-ZS(1))/2.
C
C     START ENERGY AND SOLUTE CALCULATIONS FOR EACH NODE
C
C     CALCULATE THE WATER INFILTRATED INTO THE CURRENT NODE
   50 INFMAX=(PSAT(I)*SAT(I) - VLCDT(I) - VICDT(I))*DZ
      IF (INFMAX .LT. 0.0) INFMAX=0.0
      IF (INFMAX .GE. INFIL) INFMAX=INFIL
C
C     ADJUST THE SOLUTE CONCENTRATION IN NODE TO ACCOUNT FOR LEACHING
      DO 60 J=1,NSALT
C        CALCULATE SOLUTE CONCENTRATION OF WATER COMING INTO NODE
C        AND ADJUST SALT CONCENTRATION IN NODE FOR THE WATER
C        REQUIRED TO GET TO MAXIMUM WATER CONTENT
         IF (INFIL .LE. 0.0) THEN
            CLEACH = 0.0
          ELSE
            CLEACH = SLOSS/RHOL/INFIL
         END IF
         SALTDT(J,I) = SALTDT(J,I) + CLEACH*RHOL*INFMAX/RHOB(I)
C        CALCULATE SALT PRESENT IN NODE AFTER LEACHING
         SKQVLC= SALTKQ(J,I) + (VLCDT(I)+INFMAX/DZ)*RHOL/RHOB(I)
         SLTNEW= (CLEACH*RHOL*(INFIL-INFMAX)/DZ +
     >           SALTDT(J,I)*(RHOB(I)-WT*RHOL*(INFIL-INFMAX)/DZ/SKQVLC))
     >                      /(RHOB(I)+WDT*RHOL*(INFIL-INFMAX)/DZ/SKQVLC)
         IF (SLTNEW .LT. 0.0) SLTNEW=0.0
         SLOSS= (SALTDT(J,I) - SLTNEW)*RHOB(I)*DZ
         SALTDT(J,I)= SLTNEW
   60 CONTINUE
C
C     CALCULATE THE HEAT CAPACITY OF THE RAIN AND OF THE SOIL
C     --RAIN HEAT CAPACITY IS BASED ON TOTAL WATER ABSORBED AND PASSING
C     THROUGH CURRENT NODE (TEMP WATER LEAVING NODE IS SOIL TEMP)
      SHEAT= CS(I)*DZ
      RAINHT= INFIL*RHOL*CL
C
      IF (ICESDT(I) .LE. 0) THEN
C        LAYER IS NOT FROZEN -- CALCULATE THE TEMPERATURE DIRECTLY AND
C        ADJUST MOISTURE CONTENT, MATRIC POTENTIAL, AND SOLUTES
         TSDT(I)= TSDT(I) + RAINHT/(RAINHT+SHEAT)*(TRAIN-TSDT(I))
         VLCDT(I)= VLCDT(I) + INFMAX/DZ
C        DO NOT ALLOW ADJUSTMENT OF SATURATED POTENTIALS
         IF (MATDT(I).LT.ENTRY(I)) 
     >      CALL MATVL1 (I,MATDT(I),VLCDT(I),DUMMY)
         DO 70 J=1,NSALT
            CONCDT(J,I)=SALTDT(J,I)/(SALTKQ(J,I)+VLCDT(I)*RHOL/RHOB(I))
   70    CONTINUE
C
       ELSE
C        LAYER CONTAINS ICE -- MUST ACCOUNT FOR LATENT HEAT OF FUSION.
C        CALCULATE TEMPERATURE AT WHICH LAYER WILL BE COMPLETELY THAWED
         OLDVLC= VLCDT(I)
         VLCDT(I)= VLCDT(I) + VICDT(I)*RHOI/RHOL + INFMAX/DZ
         CALL MATVL1 (I,MATDT(I),VLCDT(I),DUMMY)
         TLCONC=0.0
         DO 80 J=1,NSALT
            CONCDT(J,I)= SALTDT(J,I)/(SALTKQ(J,I)+VLCDT(I)*RHOL/RHOB(I))
            TLCONC=TLCONC + CONCDT(J,I)
   80    CONTINUE
         TOTPOT= MATDT(I) - TLCONC*UGAS*273.16/G
         TMPFRZ= 273.16*TOTPOT/(LF/G-TOTPOT)
C
C        CHECK IF HEAT OF RAIN IS SUFFICIENT TO MELT ALL OF THE ICE
         AVHEAT= RAINHT*(TRAIN-TMPFRZ)
         RQHEAT= RHOI*LF*VICDT(I)*DZ + SHEAT*(TMPFRZ-TSDT(I))
         IF (AVHEAT .GE. RQHEAT) THEN
C           AVAILABLE HEAT IS MORE THAN REQUIRED TO MELT THE ICE.
C           CALCULATE TEMPERATURE OF SOIL
            TSDT(I)= TMPFRZ + (AVHEAT-RQHEAT)/(SHEAT+RAINHT)
            VICDT(I)= 0.0
C           DONE WITH THIS LAYER -- GO TO NEXT LAYER
            GO TO 110
         END IF
C
C        LAYER REMAINS FROZEN - ASSUME TEMPERATURE IS AT FREEZING POINT
C        AND CALCULATE THE ICE CONTENT FROM AN ENERGY BALANCE:
C        (LATENT HEAT)*(DELTA ICE) =
C                        RAINHT*(TRAIN-T(NEW)) + SOILHT(T(NEW) - TSDT)
         TSOIL= TSDT(I)
         OLDVIC= VICDT(I)
         VICMAX= VICDT(I) + INFMAX/DZ*(RHOL/RHOI)
         VICDT(I)= OLDVIC - (RAINHT*(TRAIN-TMPFRZ)-SHEAT*(TMPFRZ-TSOIL))
     >                      /(RHOI*LF*DZ)
         IF (VICDT(I) .GT. VICMAX) VICDT(I) = VICMAX
         VLCDT(I)= OLDVLC + (VICMAX-VICDT(I))*RHOI/RHOL
C
C        DETERMINE MATRIC, SOLUTES AND TEMP FOR THIS WATER CONTENT
         CALL MATVL1 (I,MATDT(I),VLCDT(I),DUMMY)
         TLCONC=0.0
         DO 90 J=1,NSALT
            CONCDT(J,I)= SALTDT(J,I)/(SALTKQ(J,I)+VLCDT(I)*RHOL/RHOB(I))
            TLCONC=TLCONC + CONCDT(J,I)
   90    CONTINUE
         TOTPOT= MATDT(I) - TLCONC*UGAS*273.16/G
         TSDT(I)= 273.16*TOTPOT/(LF/G-TOTPOT)
C
C        CALCULATE ICE CONTENT USING UPDATED TEMPERATURE
         VICDT(I)= OLDVIC-(RAINHT*(TRAIN-TSDT(I))-SHEAT*(TSDT(I)-TSOIL))
     >                      /(RHOI*LF*DZ)
         IF (VICDT(I) .GT. VICMAX) VICDT(I) = VICMAX
         VLCDT(I)= OLDVLC + (VICMAX-VICDT(I))*RHOI/RHOL
C
C        DETERMINE MATRIC, SOLUTES AND TEMP FOR THIS WATER CONTENT
         CALL MATVL1 (I,MATDT(I),VLCDT(I),DUMMY)
         TLCONC=0.0
         DO 100 J=1,NSALT
            CONCDT(J,I)= SALTDT(J,I)/(SALTKQ(J,I)+VLCDT(I)*RHOL/RHOB(I))
            TLCONC=TLCONC + CONCDT(J,I)
  100    CONTINUE
         TOTPOT= MATDT(I) - TLCONC*UGAS*273.16/G
         TSDT(I)= 273.16*(LF/G/(LF/G-TOTPOT) - 1.0)
      END IF
C
C     CALCULATE CONDITIONS OF NEXT LAYER
  110 I = I + 1
      INFIL=INFIL-INFMAX
C     ADD INFILTRATION BETWEEN NODES TO TOTAL FLOW BETWEEN NODES
      TOTFLO(I-1)=TOTFLO(I-1)+INFIL
      IF (I .LT. NS) THEN
         DZ=(ZS(I+1) - ZS(I-1))/2.
         TRAIN = TSDT(I-1)
         IF (INFIL .GT. 0.0) GO TO 50
      END IF
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE WBALNC (NPLANT,NC,NSP,NR,NS,LVLOUT,JULIAN,HOUR,YEAR,
     >  ITYPE,INITAL,ZC,WCAN,WCANDT,PCAN,PCANDT,VAPC,VAPCDT,RHOSP,DZSP,
     >  DLWDT,WLAG,STORE,ZR,GMC,GMCDT,VAPR,VAPRDT,RHOR,ZS,VLC,VLCDT,
     >  VIC,VICDT,TOTFLO,PRECIP,RUNOFF,POND,EVAP1,ETSUM)
C
C     THIS SUBROUTINE SUMS THE EVAPORATION AND DEEP PERCOLATION AT THE
C     END OF EACH HOUR, THEN PRINTS THE SUMMARY AT THE DESIRED OUTPUT
C     INTERVAL
C     对指定时间间隔进行水量平衡计算
C***********************************************************************
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /CLAYRS/ DRYCAN(8,10),CANLAI(8,10),TOTLAI(8),IEVAP(8),
     >                RLEAF(8,10),RROOT(8,50),ROOTDN(8,50),TOTROT(8)
      REAL LF,LS,LV
C
      REAL ZC(11),WCAN(10),WCANDT(10),PCAN(8),PCANDT(8),VAPC(11),
     >     VAPCDT(11),RHOSP(100),DZSP(100),DLWDT(100),WLAG(11),ZR(10),
     >     GMC(10),GMCDT(10),VAPR(10),VAPRDT(10),RHOR(10),ZS(50),
     >     VLC(50),VLCDT(50),VIC(50),VICDT(50),TOTFLO(50)
      INTEGER HOUR,YEAR,ITYPE(8)
      SAVE RAIN,DPCAN,DCAN,DSNOW,DRES,DSOIL,TRUNOF,POND2,TPERC,TETSUM,
     >     TEVAP,CUMVAP,SWE
      DATA RAIN,DPCAN,DCAN,DSNOW,DRES,DSOIL,TRUNOF,POND2,TPERC,TETSUM,
     >     TEVAP,CUMVAP/12*0.0/
C
      IF (INITAL .EQ. 0) THEN
         WRITE (26,100)
C        CALCULATE THE SNOW WATER EQUIVALENT FOR THE INITIAL CONDITIONS
         SWE=0.0
         DO 5 I=1,NSP
            SWE=SWE + RHOSP(I)*DZSP(I)/RHOL + DLWDT(I)
    5    CONTINUE
         RETURN
      END IF
C
C
C     END OF THE HOUR -- DETERMINE THE CHANGE IN STORAGE
C
C     CHANGE IN STORAGE OF CANOPY
      DO 15 I=1,NC
         DO 10 J=1,NPLANT
C           CALCULATE CHANGE IN WATER CONTENT OF DEAD PLANT MATERIAL
            IF (ITYPE(J) .EQ. 0) DCAN = DCAN +
     >      (WCANDT(I)-WCAN(I))*DRYCAN(J,I)*1000./RHOL
   10    CONTINUE
C        INCLUDE CHANGE IN VAPOR DENSITY OF AIR SPACE
         IF (I .EQ. 1)  THEN
           DZ = (ZC(2) - ZC(1))/2.
          ELSE
           DZ = (ZC(I+1) - ZC(I-1))/2.
         END IF
         DCAN = DCAN + DZ*(VAPCDT(I)-VAPC(I))*1000./RHOL
   15 CONTINUE
C
C     CHANGE IN SNOWPACK
C     DLWDT (I)    depth of liquid water for i-th snow pack node at end of time step
      SWEDT = 0.0
      IF (NSP .GT. 0) THEN
         DO 20 I=1,NSP
            SWEDT = SWEDT + RHOSP(I)*DZSP(I)/RHOL + DLWDT(I)
   20    CONTINUE
C        WATER IN PROCESS OF BEING LAGGED THROUGH THE SNOWPACK
         SWEDT = SWEDT + STORE
         DO 25 I=1,11
            SWEDT = SWEDT + WLAG(I)
   25    CONTINUE
      END IF
C     CALCULATE CHANGE IN WATER CONTENT OF SNOWPACK IN MILLIMETERS
      DSNOW = DSNOW + (SWEDT - SWE)*1000.
      SWE = SWEDT
C
C     CHANGE IN STORAGE OF RESIDUE
      DO 30 I=1,NR
         IF (I .EQ. 1 .OR.  I .EQ. NR) THEN
            IF (NR .EQ.1) THEN
               DZ = ZR(NR+1)
              ELSE
               IF (I .EQ. 1)  DZ = (ZR(2) - ZR(1))/2.
               IF (I .EQ. NR) DZ = ZR(NR+1)-ZR(NR)+(ZR(NR)-ZR(NR-1))/2.
            END IF
          ELSE
           DZ = (ZR(I+1) - ZR(I-1))/2.
         END IF
         DRES = DRES + DZ*((GMCDT(I)-GMC(I))*RHOR(I)
     >                    +(VAPRDT(I)-VAPR(I)))*1000./RHOL
   30 CONTINUE
C
C     CHANGE IN STORAGE FOR THE SOIL
      DO 40 I=1,NS-1
         IF (I .EQ. 1) THEN
           DZ = (ZS(2) - ZS(1))/2.
          ELSE
           DZ = (ZS(I+1) - ZS(I-1))/2.
         END IF
         DSOIL= DSOIL + DZ*1000*(VLCDT(I) - VLC(I)
     >                        + (VICDT(I) - VIC(I))*RHOI/RHOL)
   40 CONTINUE
C
C     COMPUTE THE CHANGE IN PRECIP INTERCEPTED ON PLANT LEAVES
      TPCAN=0.0
      DO 50 J=1,NPLANT
         IF (ITYPE(J).NE.0) TPCAN=TPCAN + PCANDT(J)
         IF (ITYPE(J).NE.0) DPCAN=DPCAN + PCANDT(J)-PCAN(J)
   50 CONTINUE
C
      RAIN = RAIN + PRECIP
      TRUNOF = TRUNOF + RUNOFF
      TETSUM = TETSUM + ETSUM
      TEVAP = TEVAP + EVAP1
      TPERC = TPERC + TOTFLO(NS-1)
C
C     RETURN IF HOURLY OUTPUT IS NOT REQUIRED AND IT IS NOT END OF DAY
      IF (MOD(HOUR,LVLOUT) .NE. 0) RETURN
C
C     CONVERT TO MILLIMETERS
      TPCAN=TPCAN*1000.
      DPCAN=DPCAN*1000.
      DPOND = POND*1000. - POND2
      POND2 = POND*1000.
      RAIN = RAIN*1000.
      TRUNOF = TRUNOF*1000.
      TETSUM = TETSUM*1000.
      TEVAP = TEVAP*1000.
      TPERC = TPERC*1000.
      CUMVAP = CUMVAP + TEVAP
      ERROR = RAIN+TEVAP-TPERC-TRUNOF-DPOND-DPCAN-DCAN-DSNOW-DRES-DSOIL
C
      WRITE (26,110)JULIAN,HOUR,YEAR,RAIN,TPCAN,-TEVAP,TETSUM,DCAN,
     >              DSNOW,DRES,DSOIL,TPERC,TRUNOF,POND2,-CUMVAP,ERROR
C
      RAIN=0.0
      DPCAN=0.0
      TRUNOF=0.0
      TETSUM=0.0
      TEVAP=0.0
      TPERC=0.0
      DCAN=0.0
      DSNOW=0.0
      DRES=0.0
      DSOIL=0.0
      RETURN
C
  100 FORMAT (/35X,'SUMMARY OF WATER BALANCE',//,54X,'CHANGE IN STORAGE'
     >,/,21X,'  PRECIP',8X,'  PLANT   ',30('-'),2X,'DEEP',20X,' CUM.',/,
     >' DAY HR  YR   PRECIP   INTRCP    ET    TRANSP  CANOPY   SNOW  RES
     >IDUE   SOIL   PERC   RUNOFF   PONDED    ET    ERROR',
     >/,11X,13(6X,'MM'))
  110 FORMAT (1X,I3,I3,I5,11(1X,F7.2),1X,F7.1,1X,F7.2)
      END
C***********************************************************************
C
      SUBROUTINE ENERGY (NSP,LVLOUT,JULIAN,HOUR,YEAR,INITAL,
     >    TSWSNO,TLWSNO,TSWCAN,TLWCAN,TSWRES,TLWRES,TSWSOI,TLWSOI,
     >    THFLUX,EVAP1)
C
C     THIS SUBROUTINE SUMS THE SOLAR AND LONG-WAVE RADIATION AND THE
C     TURBULENT CONVECTION AT THE SURFACE, AND PRINT THE SUMMARY OF THE
C     ENERGY TRANSFER AT THE SURFACE AT THE DESIRED OUTPUT INTERVAL.
C
C***********************************************************************
      COMMON /TIMEWT/ WT,WDT,DT
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      REAL LF,LV,LS
      INTEGER HOUR,YEAR
C
      SAVE SSWSNO,SLWSNO,SSWCAN,SLWCAN,SSWRES,SLWRES,SSWSOI,SLWSOI,
     >     SHFLUX,SLATNT
C
      DATA SSWSNO,SLWSNO,SSWCAN,SLWCAN,SSWRES,SLWRES,SSWSOI,SLWSOI,
     >     SHFLUX,SLATNT / 10*0.0/
C
      IF (INITAL .EQ. 0) THEN
         WRITE (25,100)
      END IF
C
C
      SSWCAN = SSWCAN + TSWCAN
      SLWCAN = SLWCAN + TLWCAN
      SSWSNO = SSWSNO + TSWSNO
      SLWSNO = SLWSNO + TLWSNO
      SSWRES = SSWRES + TSWRES
      SLWRES = SLWRES + TLWRES
      SSWSOI = SSWSOI + TSWSOI
      SLWSOI = SLWSOI + TLWSOI
      SHFLUX = SHFLUX + THFLUX
C     LATENT HEAT
      IF (NSP .GT. 0) THEN
C        USE LATENT HEAT OF SUBLIMATION
         SLATNT = SLATNT + LS*EVAP1*RHOL/1000.
        ELSE
         SLATNT = SLATNT + LV*EVAP1*RHOL/1000.
      END IF

C     RETURN IF HOURLY OUTPUT IS NOT REQUIRED AND IT IS NOT END OF DAY
      IF (MOD(HOUR,LVLOUT) .NE. 0) RETURN
C
      TOTSWR = SSWCAN + SSWSNO + SSWRES + SSWSOI
      TOTLWR = SLWCAN + SLWSNO + SLWRES + SLWSOI
C
      WRITE (25,110) JULIAN,HOUR,YEAR,SSWCAN,SSWSNO,SSWRES,SSWSOI,
     >      TOTSWR,SLWCAN,SLWSNO,SLWRES,SLWSOI,TOTLWR,SHFLUX,SLATNT
C
      SSWCAN = 0.0
      SSWSNO = 0.0
      SSWRES = 0.0
      SSWSOI = 0.0
      SLWCAN = 0.0
      SLWSNO = 0.0
      SLWRES = 0.0
      SLWSOI = 0.0
      SHFLUX = 0.0
      SLATNT = 0.0
      RETURN
C
  100 FORMAT (/33X,' SUMMARY OF ENERGY TRANSFER AT THE SURFACE',///,30X,
     >' SOLAR RADIATION',25X,'LONG-WAVE RADIATION',/,15X,42('-'),4X,41
     >('-'),/,' DAY HR  YR    CANOPY     SNOW   RESIDUE    SOIL   TOTAL'
     >,'     CANOPY    SNOW   RESIDUE    SOIL    TOTAL  SENSIBLE  ',
     >'LATENT',/,12X,10('    KJ/M2'),1X,2('   KJ/M2'))
  110 FORMAT (I4,I3,I5,12(1X,F8.0))
      END
C***********************************************************************
C
      SUBROUTINE OUTPUT (NPLANT,NC,NSP,NR,NS,LVLOUT,NOTALL,INPH2O,
     >   JULIAN,HOUR,YEAR,INITAL,ZC,TCDT,VAPCDT,WCANDT,ROOTXT,RHOSP,
     >   ZSP,TSPDT,DLWDT,ZR,TRDT,VAPRDT,GMCDT,ZS,TSDT,VLCDT,VICDT,
     >   MATDT,TOTFLO,CONCDT,SALTDT)
C
C     THIS SUBROUTINE PRINTS THE TEMPERATURE AND MOISTURE PROFILES AT
C     DESIRED INTERVALS.
C
C***********************************************************************
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      COMMON /CLAYRS/ DRYCAN(8,10),CANLAI(8,10),TOTLAI(8),IEVAP(8),
     >                RLEAF(8,10),RROOT(8,50),ROOTDN(8,50),TOTROT(8)
      COMMON /SLPARM/ B(50),ENTRY(50),RHOB(50),SAT(50),SATK(50),
     >                NSALT,SALTKQ(10,50),VAPCOF(50),VAPEXP(50),
     >                SAND(50),SILT(50),CLAY(50),OM(50)
      COMMON /MEASUR/ MEASDY,MEASHR,VLCMES(50),TSMEAS(50)
      REAL LF,LS,LV
C
      REAL ZC(11),TCDT(11),VAPCDT(11),WCANDT(10),ROOTXT(50),RHOSP(100),
     >     ZSP(100),TSPDT(100),DLWDT(100),ZR(10),TRDT(10),VAPRDT(10),
     >     GMCDT(10),ZS(50),TSDT(50),VLCDT(50),VICDT(50),MATDT(50),
     >     TOTFLO(50),CONCDT(10,50),SALTDT(10,50)
      REAL AVGTMP(50),AVGMAT(50),AVGVLC(50),AVGSLT(10,50),AVGCON(10,50),
     >     SUMFLO(50),SUMRXT(50)
      INTEGER HOUR,YEAR,LVLOUT(15)
      SAVE NOUT,NAVTMP,NAVMAT,NAVVLC,NAVSLT,NAVCON,
     >     AVGTMP,AVGMAT,AVGVLC,SUMFLO,SUMRXT,AVGSLT,AVGCON
C
C     IF THIS IS THE FIRST TIME INTO SUBROUTINE, INITIALIZE SUMMATIONS 
C     FOR AVG TEMPERATURE, MOISTURE & SALTS, WATER FLOW BETWEEN SOIL 
C     LAYERS, AND TOTAL ROOT EXTRACTION
      IF (INITAL .EQ. 0) THEN
         NOUT=1
         DO 10 I=1,NS
            AVGTMP(I)=0.0
            AVGMAT(I)=0.0
            AVGVLC(I)=0.0
            TOTFLO(I)=0.0
            SUMFLO(I)=0.0
            ROOTXT(I)=0.0
            SUMRXT(I)=0.0
            DO 5 J=1,NSALT
               AVGSLT(J,I)=0.0
               AVGCON(J,I)=0.0
    5       CONTINUE
   10    CONTINUE
         NAVTMP=0
         NAVMAT=0
         NAVVLC=0
         NAVSLT=0
         NAVCON=0
C
         IF (LVLOUT(2) .GT. 0) WRITE (22,175)
         IF (LVLOUT(3) .GT. 0) THEN
            WRITE (23,*) 
     >      'AVERAGE PREDICTED TEMPERATURES FOR EACH NODE (C)'
            WRITE (23,180) (ZS(J), J=1,NS)
         END IF
         IF (LVLOUT(4) .GT. 0) THEN
           WRITE(24,*)
     >     'AVERAGE PREDICTED MOISTURE CONTENT FOR EACH NODE (M3/M3)'
           WRITE (24,180) (ZS(J), J=1,NS)
         END IF
         IF (LVLOUT(5) .GT. 0) THEN
           WRITE(31,*)
     >     'AVERAGE PREDICTED MATRIC POTENTIAL FOR EACH NODE (METERS)'
           WRITE (31,182) (ZS(J), J=1,NS)
         END IF
         IF (LVLOUT(8) .GT. 0) THEN
           WRITE (32,*)
     >     'TOTAL WATER FLUX BETWEEN SOIL NODES (MM; DOWNWARD POSITIVE)'
           WRITE (32,181) (ZS(J), J=1,NS)
         END IF
         IF (LVLOUT(9) .GT. 0) THEN
           WRITE (30,*) 'WATER EXTRACTED BY ROOTS FROM EACH LAYER (M)'
           WRITE (30,182) (ZS(J), J=1,NS)
         END IF
         IF (NSALT .GT. 0) THEN
C           WRITE HEADERS FOR SOLUTE OUTPUT FILES
            IF (LVLOUT(11) .GT. 0) THEN
               WRITE (28,*) 'AVERAGE PREDICTED SALT CONCENTRATION ',
     >         'FOR EACH NODE (MOLES/KG OF SOIL)'
               WRITE (28,185) (ZS(J), J=1,NS)
            END IF
            IF (LVLOUT(12) .GT. 0) THEN
               WRITE (29,*) 'AVERAGE PREDICTED SOLUTE CONCENTRATION ', 
     >         'FOR EACH NODE (MOLES/LITER OF SOIL WATER)'
               WRITE (29,185) (ZS(J), J=1,NS)
            END IF
         END IF
      END IF
C
C     IF NOTALL = 1, THIS INDICATES THAT NOT EVERYTHING IS TO PRINTED --
C     ONLY THE FULL PROFILE TO THE GENERAL OUTPUT FILE
      IF (NOTALL .EQ. 1) GO TO 70
C
C
      IF (HOUR.EQ.MEASHR .AND.JULIAN.EQ.MEASDY .AND.LVLOUT(2).GT.0) THEN
C        PRINT OUT COMPARISON OF MEASURED VS PREDICTED PROFILES
         DO 15 I=1,NS
            TOTVLC = VLCDT(I) + VICDT(I)*RHOI/RHOL
            IF (INPH2O .EQ. 1) THEN
C              INPUT SOIL MOISTURE IS MATRIC POTENTIAL - CONVERT
               SOIMAT=VLCMES(I)
               CALL MATVL2 (I,SOIMAT,VLCMES(I),DUMMY)
            END IF
            IF (NSALT .GT. 0) THEN
               WRITE (22,110) NOUT,JULIAN,HOUR,YEAR,ZS(I),VLCMES(I),
     >                 TOTVLC,TSMEAS(I),TSDT(I),(SALTDT(J,I), J=1,NSALT)
             ELSE
               WRITE (22,110) NOUT,JULIAN,HOUR,YEAR,ZS(I),VLCMES(I),
     >                 TOTVLC,TSMEAS(I),TSDT(I)
            END IF
   15    CONTINUE
         NOUT=NOUT+1
         WRITE (22,*)
      END IF
C
      DO 25 I=1,NS
         AVGTMP(I)=AVGTMP(I)+TSDT(I)
         AVGMAT(I)=AVGMAT(I)+MATDT(I)
         AVGVLC(I)=AVGVLC(I)+VLCDT(I)+VICDT(I)*RHOI/RHOL !====
         SUMFLO(I)=SUMFLO(I)+TOTFLO(I)
         SUMRXT(I)=SUMRXT(I)+ROOTXT(I)
         DO 20 J=1,NSALT
            AVGSLT(J,I)=AVGSLT(J,I)+SALTDT(J,I)
            AVGCON(J,I)=AVGCON(J,I)+CONCDT(J,I)
   20    CONTINUE
   25 CONTINUE
      NAVTMP=NAVTMP+1
      NAVMAT=NAVMAT+1
      NAVVLC=NAVVLC+1
      NAVSLT=NAVSLT+1
      NAVCON=NAVCON+1
C
C     PRINT OUT SOIL TEMPERATURES IF AT DESIRED INTERVAL
      IF (LVLOUT(3) .GT. 0) THEN
         IF (MOD(HOUR,LVLOUT(3)).EQ.0 .OR. INITAL.EQ.0) THEN
            WRITE (23,120) JULIAN,HOUR,YEAR,(AVGTMP(I)/NAVTMP, I=1,NS)
            DO 30 I=1,NS
               AVGTMP(I)=0.0
   30       CONTINUE
            NAVTMP=0
         END IF
      END IF
C
C     PRINT OUT SOIL MOISTURE CONTENT IF AT DESIRED INTERVAL
      IF (LVLOUT(4) .GT. 0) THEN
         IF (MOD(HOUR,LVLOUT(4)).EQ.0 .OR. INITAL.EQ.0) THEN
            WRITE (24,130) JULIAN,HOUR,YEAR,(AVGVLC(I)/NAVVLC, I=1,NS)
            DO 35 I=1,NS
               AVGVLC(I)=0.0
   35       CONTINUE
            NAVVLC=0
         END IF
      END IF
C
C     PRINT OUT SOIL MATRIC POTENTIAL IF AT DESIRED INTERVAL
      IF (LVLOUT(5) .GT. 0) THEN
         IF (MOD(HOUR,LVLOUT(5)).EQ.0 .OR. INITAL.EQ.0) THEN
            WRITE (31,135) JULIAN,HOUR,YEAR,(AVGMAT(I)/NAVMAT, I=1,NS)
            DO 40 I=1,NS
               AVGMAT(I)=0.0
   40       CONTINUE
            NAVMAT=0
         END IF
      END IF
C
C     PRINT OUT WATER FLOW BETWEEN SOIL NODES IF AT DESIRED INTERVAL
      IF (LVLOUT(8) .GT. 0) THEN
         IF (MOD(HOUR,LVLOUT(8)).EQ.0 .AND. INITAL.NE.0) THEN
            WRITE (32,125) JULIAN,HOUR,YEAR,(SUMFLO(I)*1000.,I=1,NS-1)
            DO 42 I=1,NS
               SUMFLO(I)=0.0
   42       CONTINUE
         END IF
      END IF
C
C     PRINT OUT EXTRACTED BY ROOTS FROM EACH NODE
      IF (LVLOUT(9) .GT. 0 .AND. NPLANT .GT. 0) THEN
         IF (MOD(HOUR,LVLOUT(9)).EQ.0 .OR. INITAL.EQ.0) THEN
            WRITE (30,137) JULIAN,HOUR,YEAR,(SUMRXT(I)/RHOL, I=1,NS)
            DO 43 I=1,NS
               SUMRXT(I)=0.0
   43       CONTINUE
         END IF
      END IF
C
C     PRINT OUT SALT CONCENTRATION IF AT DESIRED INTERVAL
      IF (LVLOUT(11) .GT. 0) THEN
        IF (MOD(HOUR,LVLOUT(11)).EQ.0 .OR. INITAL.EQ.0) THEN
          DO 50 J=1,NSALT
            WRITE (28,170) JULIAN,HOUR,YEAR,J,
     >                     (AVGSLT(J,I)/NAVSLT, I=1,NS)
            DO 45 I=1,NS
              AVGSLT(J,I)=0.0
   45       CONTINUE
   50     CONTINUE
          NAVSLT=0
        END IF
      END IF
C
C     PRINT OUT SOLUTE CONCENTRATION IF AT DESIRED INTERVAL
      IF (LVLOUT(12) .GT. 0) THEN
        IF (MOD(HOUR,LVLOUT(12)).EQ.0 .OR. INITAL.EQ.0) THEN
          DO 60 J=1,NSALT
            WRITE (29,170) JULIAN,HOUR,YEAR,J,
     >                     (AVGCON(J,I)/NAVCON, I=1,NS)
            DO 55 I=1,NS
              AVGCON(J,I)=0.0
   55       CONTINUE
   60     CONTINUE
          NAVCON=0
        END IF
      END IF
C
C
      IF (INITAL .EQ. 0) GO TO 70
      IF (LVLOUT(1).EQ.0) GO TO 80
      IF (MOD(HOUR,LVLOUT(1)).NE.0) GO TO 80
C
   70 WRITE (21,*)
C
C     CANOPY LAYERS
      IF (NC .GT. 0) THEN
         WRITE (21,*) 'CANOPY LAYERS :'
         WRITE (21,145) (' LAI#',J, J=1,NPLANT)
         WRITE (21,*) '              (M)     (C)  (KG/M3)  (KG/KG)'

         DO 75 I=1,NC
            WRITE (21,150) JULIAN,HOUR,YEAR,ZC(I),TCDT(I),VAPCDT(I),
     >                  WCANDT(I),(CANLAI(J,I), J=1,NPLANT)
   75    CONTINUE
         WRITE (21,150) JULIAN,HOUR,YEAR,ZC(NC+1)
      END IF
C
C     SNOWPACK LAYERS
      IF (NSP .GT. 0) THEN
         WRITE (21,*) 'SNOW LAYERS :'
         WRITE (21,*) ' DAY HR  YR  DEPTH  TEMP   LIQUID  DENSITY'
         WRITE (21,*) '              (M)    (C)     (M)   (KG/M3)'
         WRITE (21,140) (JULIAN,HOUR,YEAR,ZSP(I),TSPDT(I),DLWDT(I),
     >                  RHOSP(I), I=1,NSP)
         WRITE (21,140) JULIAN,HOUR,YEAR,ZSP(NSP+1)
      END IF
C
C      RESIDUE LAYERS
      IF (NR .GT. 0) THEN
         WRITE (21,*) 'RESIDUE LAYERS :'
         WRITE (21,*) ' DAY HR  YR  DEPTH   TEMP   VAPOR  MOISTURE'
         WRITE (21,*) '              (M)     (C)  (KG/M3)  (KG/KG)'
         DO 76 I=1,NR
            WRITE (21,150) JULIAN,HOUR,YEAR,ZR(I),TRDT(I),VAPRDT(I),
     >                  GMCDT(I)
   76    CONTINUE
         WRITE (21,150) JULIAN,HOUR,YEAR,ZR(NR+1)
      END IF
C
C     SOIL LAYER
      WRITE (21,*) 'SOIL LAYERS :'
      IF (NSALT.EQ.0) THEN
         WRITE (21,*) ' DAY HR  YR  DEPTH   TEMP WATER  ICE   MATRIC'
         WRITE (21,*) '              (M)     (C) (M/M) (M/M)    (M)'
         WRITE (21,160) (JULIAN,HOUR,YEAR,ZS(I),TSDT(I),VLCDT(I),
     >   VICDT(I),MATDT(I), I=1,NS)
       ELSE
      WRITE (21,*) ' DAY HR  YR  DEPTH   TEMP WATER  ICE   MATRIC    SOL
     >UTE #1   TOTAL SALT #1'
      WRITE (21,*) '              (M)     (C) (M/M) (M/M)    (M)      (M
     >OLE/L)     (MOLES/KG)'
      WRITE (21,161) (JULIAN,HOUR,YEAR,ZS(I),TSDT(I),VLCDT(I),VICDT(I),
     >              MATDT(I),CONCDT(1,I),SALTDT(1,I), I=1,NS)
      END IF
C
   80 RETURN
C
  110 FORMAT (1X,I2,I4,I3,I5,F6.2,3X,2F6.3,3X,2F6.1,3X,10(E10.3))
  120 FORMAT (1X,2I3,I5,50F6.1)
  125 FORMAT (1X,2I3,I5,50F6.2)
  130 FORMAT (1X,2I3,I5,50F6.3)
  135 FORMAT (I4,I3,I5,50F10.1)
  137 FORMAT (I4,I3,I5,50E10.2)
  140 FORMAT (I4,I3,I5,F6.3,F7.3,F8.4,F8.2)
  145 FORMAT (' DAY HR  YR  DEPTH   TEMP   VAPOR  MOISTURE ',8(A,I1))
  150 FORMAT (I4,I3,I5,F6.3,F7.3,F10.6,F7.3,1X,8F6.2)
  160 FORMAT (I4,I3,I5,F6.3,F7.3,F6.3,F6.3,F8.2)
  161 FORMAT (I4,I3,I5,F6.3,F7.3,F6.3,F6.3,F8.2,E14.5,E14.5)
  170 FORMAT (1X,2I3,I5,I3,50E10.3)
  175 FORMAT (' COMPARISON OF PREDICTED AND MEASURED PROFILES',//,
     >'  N  DY HR  YR   DEPTH     MOISTURE      TEMPERATURE    SOLUTES',
     >/,26X,'MEAS  PRED     MEAS  PRED    PREDICTED')
  180 FORMAT ('  DY HR  YR ',50F6.2)
  181 FORMAT ('  DY HR  YR ',F3.1,50F6.2)
  182 FORMAT ('  DY HR  YR ',50F10.2)
  185 FORMAT ('  DY HR  YR   #',50F10.2)
      END
C***********************************************************************
C
      SUBROUTINE FROST (NSP,NS,LVLOUT,JULIAN,HOUR,YEAR,INITAL,
     >           ZSP,RHOSP,DZSP,DLWDT,WLAG,STORE,ZS,VLCDT,VICDT,ICESDT)
C
C     THIS SUBROUTINE INTERPOLATES BETWEEN NODES TO DETERMINE THE FROST
C     DEPTH, THEN PRINT THE FROST DEPTH AND SNOW DEPTH AT THE DESIRED
C     OUTPUT INTERVAL
C
C***********************************************************************
C
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
C
      REAL ZSP(100),RHOSP(100),DZSP(100),DLWDT(100),WLAG(11),
     >     ZS(50),VLCDT(50),VICDT(50)
      INTEGER HOUR,YEAR,ICESDT(50)
C
      SAVE NPRINT,LAST,LDAY,LHOUR,LYR,ZERO,FDEPTH,TDEPTH
C
      IF (INITAL .EQ. 0) THEN
C        LIMIT PRINT OUT OF ICE CONTENT TO TOP 15 SOIL NODES -- ANY MORE
C        COULD NOT FIT ON ONE LINE.
         NPRINT=NS
         IF (NPRINT .GT. 15) NPRINT=15
C
         WRITE (27,100) (ZS(I), I=1,NPRINT)
C        LAST => FLAG STATING WHETHER THERE WAS FROST LAST TIME STEP
C        LDAY,LHOUR,LYR => LAST TIME THAT FROST CHECKED
         LAST = 0
         LDAY = 0
         LHOUR = 0
         LYR = 0
         ZERO = 0.0
         FDEPTH = 0.0
         TDEPTH = 0.0
         GO TO 5
      END IF
C
      IF (MOD(HOUR,LVLOUT) .NE. 0) RETURN
C
C     FIND LAYERS OF MAXIMUM THAW AND FROST
    5 NTHAW = 0
      NFROST = 0
      DO 10 I=NS,1,-1
         IF (ICESDT(I) .EQ. 1) THEN
            IF (NFROST .EQ. 0) NFROST = I
          ELSE
            IF (NFROST .GT. 0 .AND. NTHAW .EQ. 0) THEN
               NTHAW = I+1
               GO TO 15
            END IF
         END IF
   10 CONTINUE
C
   15 IF (NFROST .EQ. 0  .AND.  LAST .EQ. 0  .AND.  NSP .LE. 0) THEN
C        SAVE TIME ON WHICH FROST WAS LAST CHECKED AND RETURN
         LDAY=JULIAN
         LHOUR=HOUR
         LYR=YEAR
C     changes commented
C         RETURN
C     changes commented
      END IF

C     changes commented
C     IF LAST = 0, SNOW OR FROST WAS NOT PRESENT LAST TIME STEP BUT IS
C     PRESENT THIS TIME STEP --> PRINT OUT ZEROES FOR LAST TIME STEP
C      IF (LAST .EQ. 0  .AND.  LDAY .NE. 0)
C     >    WRITE (27,110) LDAY,LHOUR,LYR,(ZERO,I=1,NPRINT+4)
C      LDAY=JULIAN
C      LHOUR=HOUR
C      LYR=YEAR
C     changes commented

      IF (NFROST .GT. 0) THEN
C        CALCULATE DEPTH OF FROST
         LAST = 1
         NF = NFROST
         FRACTN = VICDT(NF)/(VLCDT(NF) + VICDT(NF))
         IF (NF .EQ. 1  .OR.  NF .EQ. NS) THEN
            IF (NF .EQ. 1) THEN
              FDEPTH=  FRACTN*(ZS(2)-ZS(1))/2.
             ELSE
              FDEPTH= (ZS(NS)+ZS(NS-1))/2. + FRACTN*(ZS(NS)-ZS(NS-1))/2.
            END IF
          ELSE
            FDEPTH= (ZS(NF)+ZS(NF-1))/2. + FRACTN*(ZS(NF+1)-ZS(NF-1))/2.
          END IF
C
         TDEPTH = 0.0
         IF (NTHAW .GT. 0) THEN
C           CALCULATE DEPTH OF THAW
            NT = NTHAW
            FRACTN = VLCDT(NT)/(VLCDT(NT) + VICDT(NT))
            IF (NT .EQ. 1) THEN
              TDEPTH =  FRACTN*(ZS(2)-ZS(1))/2.
             ELSE
              IF (NT .EQ. NFROST) THEN
               TDEPTH=(ZS(NT)+ZS(NT-1))/2. +
     >                FRACTN*(FDEPTH-(ZS(NT)+ZS(NT-1))/2.)
                ELSE
               TDEPTH=(ZS(NT)+ZS(NT-1))/2.+FRACTN*(ZS(NT+1)-ZS(NT-1))/2.
              END IF
            END IF
         END IF
C
       ELSE
C        NO FROST EXISTS IN THE PROFILE - IF FROST EXISTED AT THE LAST
C        OUTPUT, INDICATE THAT FROST HAS LEFT
         IF (LAST .GT. 0) THEN
C           FROST EXISTED IN PROFILE AT LAST OUTPUT -- ENSURE THAT
C           FROST DEPTH AND THAW DEPTH LINES COME TOGETHER
            IF (TDEPTH .EQ. 0.0) THEN
C              GROUND WAS THAWING FROM THE BOTTOM UP
               FDEPTH = 0.0
              ELSE
C              THAW AND FROST DEPTH COME TOGETHER SOMEWHERE BETWEEN THE
C              DEPTHS AT LAST OUTPUT - USUALLY 2/3 BETWEEN
               FDEPTH = (2.*FDEPTH + TDEPTH)/3.
               TDEPTH = FDEPTH
            END IF
          ELSE
C           NO FROST IN PROFILE AT LAST OUTPUT
            FDEPTH = 0.0
            TDEPTH = 0.0
         END IF
         LAST = 0
      END IF
C
      SWE = 0.0
      IF (NSP .LE. 0) THEN
         SNOW = 0.0
        ELSE
C        SAVE SNOW DEPTH AND CALCULATE SNOW WATER EQUIVALENT
         SNOW = ZSP(NSP+1)
         DO 20 I=1,NSP
            SWE = SWE + RHOSP(I)*DZSP(I)/RHOL + DLWDT(I)
   20    CONTINUE
C        WATER IN PROCESS OF BEING LAGGED THROUGH THE SNOWPACK
         SWE = SWE + STORE
         DO 25 I=1,11
            SWE = SWE + WLAG(I)
   25    CONTINUE
         LAST = 1
      END IF
C
C     PRINT FROST, THAW AND SNOW DEPTH IN CENTIMETERS AND SWE IN MM
      WRITE (27,110) JULIAN,HOUR,YEAR,TDEPTH*100,FDEPTH*100,SNOW*100,
     >               SWE*1000.,(VICDT(I), I=1,NPRINT)

C     changes commented
C      IF (LAST .EQ. 0) WRITE (27,*)
C     changes commented
      RETURN
C
  100 FORMAT (/4X,'FROST AND SNOW DEPTH',///,' DAY HR  YR',
     >5X,'THAW   FROST    SNOW     SWE      ICE CONTENT OF EACH NODE ',
     >'(M3/M3)',/,17X,'CM      CM      CM       MM',6X,20F6.2)
  110 FORMAT (I4,I3,I5,4F8.1,6X,15F6.3)
C
      END
C***********************************************************************
C
      SUBROUTINE TDMA (N,A,B,C,D,X)
C
C     THIS SUBROUTINE SOLVES A TRI-DIAGONAL MATRIX OF THE FORM :
C
C              | B1  C1   0   0  .  .   . | |X1|   |D1|
C              | A2  B2  C2   0  .  .   . | |X2|   |D2|
C              |  0  A3  B3  C3  .  .   . | |X3| = |D3|
C              |  .   .   .   .     .   . | | .|   | .|
C              |  0   0   0   0  . AN  BN | |XN|   |DN|
C
C***********************************************************************
      DIMENSION A(N),B(N),C(N),D(N),X(N)
C
      DO 10 I=2,N
         C(I-1)=C(I-1)/B(I-1)
         D(I-1)=D(I-1)/B(I-1)
         B(I)=B(I)-A(I)*C(I-1)
         D(I)=D(I)-A(I)*D(I-1)
   10 CONTINUE
      X(N)=D(N)/B(N)
      DO 20 I=N-1,1,-1
         X(I)=D(I)-C(I)*X(I+1)
   20 CONTINUE
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE WEIGHT (N,AVG,BEGIN,END)
C
C     THIS SUBROUTINE CALCULATES THE WEIGHTED AVERAGE OF VARIABLES AT
C     THE BEGINNING AND END OF THE TIME STEP
C
C        WT = WEIGHTING FOR VALUES AT THE BEGINNING OF THE TIME STEP
C        WDT = WEIGHTING FOR VALUES AT THE END OF THE TIME STEP
C
C***********************************************************************
      COMMON /TIMEWT/ WT,WDT,DT
      DIMENSION BEGIN(N), END(N), AVG(N)
C
      DO 10 I=1,N
         AVG(I) = WT*BEGIN(I) + WDT*END(I)
   10 CONTINUE
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE CONDUC (N,Z,RK,CON)
C
C     THIS SUBROUTINE CALCULATES THE CONDUCTANCE TERM BETWEEN TWO NODES
C     USING THE GEOMETRIC MEAN AND THE SPACE INCREMENT
C
C                                                     K(I)
C        K(I) = ( K(I)*K(I+1) )**1/2 =>  CON(I) = -------------
C                                                 Z(I+1) - Z(I)
C
C***********************************************************************
      DIMENSION Z(N),RK(N),CON(N)
C
      DO 10 I=1,N-1
            CON(I)= SQRT(RK(I)*RK(I+1))/(Z(I+1) - Z(I))
   10 CONTINUE
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE VSLOPE (S,SATV,T)
C
C     THIS SUBROUTINE CALCULATES THE SATURATED VAPOR DENSITY AND THE
C     SLOPE OF THE VAPOR DENSITY CURVE  (SATV IS IN KG/M**3)
C
C***********************************************************************
      COMMON /CONSTN/ LF,LV,LS,G,UGAS,RHOL,RHOI,RHOM,RHOOM,RHOA,CL,CI,
     >                CM,COM,CA,CV,CR,VONKRM,VDIFF,PRESUR,P0,
     >                TKL,TKI,TKA,TKR,TKSA,TKSI,TKCL,TKOM
      REAL LF,LS,LV
C
      TMP=T+273.16
      SATV=EXP(52.57633-6790.4985/TMP-5.02808*ALOG(TMP))
     >     *1000./(UGAS*TMP/.018)
      S=.0000165 + 4944.43*SATV/(TMP**2)
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE MATVLC (I,MAT,VLC,DLDM)
C
C     THIS SUBROUTINE DEFINES THE RELATION BETWEEN THE VOLUMETRIC
C     LIQUID CONTENT AND THE MATRIC POTENTIAL. THE SUBROUTINE IS DIVIDED
C     INTO THREE PART, AND THE OUTPUT DEPENDS ON WHICH PART IS CALLED
C           MATVL1 : THE MATRIC POTENTIAL IS CALCULATED FROM MOISTURE
C           MATVL2 : THE MOISTURE CONTENT IS CALCULATED FROM MATRIC
C           MATVL3 : THE DERIVATIVE OF MOISTURE CONTENT WITH RESPECT
C                    TO MATRIC POTENTIAL IS CALCULTED.
C                I = NODE NUMBER
C
C***********************************************************************
      COMMON /SLPARM/ B(50),ENTRY(50),RHOB(50),SAT(50),SATK(50),
     >                NSALT,SALTKQ(10,50),VAPCOF(50),VAPEXP(50),
     >                SAND(50),SILT(50),CLAY(50),OM(50)
      REAL MAT,VLC,DLDM
C****
      ENTRY MATVL1 (I,MAT,VLC,DLDM)
C
C     DETERMINE THE MATRIC POTENTIAL FROM THE MOISTURE CONTENT
      IF (VLC .LT. SAT(I)) THEN
         MAT=ENTRY(I)*(VLC/SAT(I))**(-B(I))
        ELSE
         MAT=ENTRY(I)
      END IF
      RETURN
C****
      ENTRY MATVL2 (I,MAT,VLC,DLDM)
C
C     DETERMINE THE MOISTURE CONTENT FROM THE MATRIC POTENTIAL
      IF (ENTRY(I) .GT. MAT) THEN
         VLC=SAT(I)*(MAT/ENTRY(I))**(-1./B(I))
       ELSE
         VLC=SAT(I)
      END IF
      RETURN
C****
      ENTRY MATVL3 (I,MAT,VLC,DLDM)
C
C     DETERMINE THE DERIVATIVE OF THE MOISTURE CONTENT WITH RESPECT
C     TO MATRIC POTENTIAL
      IF (MAT .GT. ENTRY(I)) THEN
C        LAYER IS SATURATED
         DLDM=0.0
       ELSE
C        UNSATURATED CONDITIONS
         DLDM=-VLC/B(I)/MAT
      END IF
      RETURN
C****
      END
