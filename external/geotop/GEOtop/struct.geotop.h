
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.0 Public - Version "Montebello" - Update 2 (29 April 2010)
 
 Copyright (c), 2010 - Stefano Endrizzi and Riccardo Rigon
 
 This file is part of GEOtop 1.0 Public
 
 GEOtop 1.0 Public is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.0 Public is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community at the following E-mail address: geotopusers@googlegroups.com to which you can subscribe at  http://groups.google.com/group/geotopusers/
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors
 
 */
    
 #include "turtle.h"
 #include "t_datamanipulation.h"
 #include "t_utilities.h"
 #include "tensor3D.h"
 //#include "turtle2umfpack.h"

 
/*---------------------------------------------------------------------------*/
typedef struct {

/*egy->Dlay = new_doublevector( Nl + par->snowlayer_max + par->glaclayer_max );
egy->wliq = new_doublevector( Nl + par->snowlayer_max + par->glaclayer_max );
egy->wice = new_doublevector( Nl + par->snowlayer_max + par->glaclayer_max );
egy->Temp = new_doublevector( Nl + par->snowlayer_max + par->glaclayer_max );
egy->deltaw = new_doublevector( Nl + par->snowlayer_max + par->glaclayer_max );

egy->SWlayer = new_doublevector( par->snowlayer_max+1 );

egy->soil_transp_layer = new_doublevector(land->root_fraction->nch);
initialize_doublevector(egy->soil_transp_layer, 0.);

egy->dFenergy = new_doublevector( Nl + par->snowlayer_max + par->glaclayer_max );
egy->Kth0=new_doublevector( Nl + par->snowlayer_max + par->glaclayer_max );
egy->Kth1=new_doublevector( Nl + par->snowlayer_max + par->glaclayer_max );
egy->Fenergy=new_doublevector( Nl + par->snowlayer_max + par->glaclayer_max );
egy->Newton_dir=new_doublevector( Nl + par->snowlayer_max + par->glaclayer_max );
egy->T0=new_doublevector( Nl + par->snowlayer_max + par->glaclayer_max );
egy->T1=new_doublevector( Nl + par->snowlayer_max + par->glaclayer_max );
egy->Tstar=new_doublevector(Nl); //soil temperature at which freezing begins
egy->THETA=new_doublevector(Nl);	//water content (updated in the iterations)
*/
    DOUBLEMATRIX *Rn_mean; 
	DOUBLEMATRIX *Rn_max;
	DOUBLEMATRIX *Rn_min;   
	DOUBLEMATRIX *LW_in;
	DOUBLEMATRIX *LW;
	DOUBLEMATRIX *LW_max;
	DOUBLEMATRIX *LW_min;	
	DOUBLEMATRIX *SW;	
	DOUBLEMATRIX *SW_max;			
    DOUBLEMATRIX *ET_mean; 
	DOUBLEMATRIX *ET_max;
	DOUBLEMATRIX *ET_min;
    DOUBLEMATRIX *H_mean; 
	DOUBLEMATRIX *H_max;
	DOUBLEMATRIX *H_min;    
    DOUBLEMATRIX *SEB_mean; 
	DOUBLEMATRIX *G_max;
	DOUBLEMATRIX *G_min;
	DOUBLEMATRIX *G_snowsoil;	
    DOUBLEMATRIX *Ts_mean;  /*averaged surface Temperature(on nDt_output_basin Dt time intervals)*/
	DOUBLEMATRIX *Ts_max;
	DOUBLEMATRIX *Ts_min;
	DOUBLEMATRIX *Rswdown_mean;
	DOUBLEMATRIX *Rswdown_max;
	DOUBLEMATRIX *Ta_mean;
	DOUBLEMATRIX *Ta_max;
	DOUBLEMATRIX *Ta_min;	
	
	DOUBLEMATRIX *Rswbeam;
	LONGMATRIX *nDt_shadow;/*The step numbers in shadow*/
	LONGMATRIX *nDt_sun;/*The step numbers in Sun*/
	
	DOUBLEMATRIX *Hgplot;
	DOUBLEMATRIX *LEgplot;
	DOUBLEMATRIX *Hvplot;
	DOUBLEMATRIX *LEvplot;	
	
	DOUBLEMATRIX *SWinplot;
	DOUBLEMATRIX *SWgplot;
	DOUBLEMATRIX *SWvplot;
	
	DOUBLEMATRIX *LWinplot;
	DOUBLEMATRIX *LWgplot;
	DOUBLEMATRIX *LWvplot;
	
	DOUBLEMATRIX *Tgplot;
	DOUBLEMATRIX *Tvplot;
	DOUBLEMATRIX *Tsplot;	

	DOUBLEMATRIX *SWin;
	DOUBLEMATRIX *LWin;
	
	double hsun;/*solar elevation angle in rad*/
	double dsun;/*solar azimuth angle, from N clockwise, radiants*/

	DOUBLEVECTOR *Dlay;/*Layer Thickness [m] for energy calculation,ns,ng,nl */
	DOUBLEVECTOR *wliq;/*liquid water content for each layer,Total layers ns+ng+nl*/
	DOUBLEVECTOR *wice;/*ice content for each layer, Total layers ns+ng+nl*/
	DOUBLEVECTOR *Temp;/*Temperature for each layer in energy balance, Total layer ns+ng+nl*/
	DOUBLEVECTOR *deltaw;/*mass that melts (if>0) or freezes (if<0), phase change*/
	DOUBLEVECTOR *SWlayer;/* Absorbed short-wave radiation in each snow layer*/
	DOUBLEVECTOR *soil_transp_layer;
	DOUBLEVECTOR *dFenergy;
	DOUBLEVECTOR *Kth0;/*Average thermal conductivity between two layers,Kth0: the first layer*/
	DOUBLEVECTOR *Kth1;/*Average thermal conductivity between two layers,Kth1: the second layer*/
	DOUBLEVECTOR *Fenergy;/*Net energy balance for each layer in each time step*/
	DOUBLEVECTOR *Newton_dir;
	DOUBLEVECTOR *T0;/*Initial Temp for each iterate step*/
	DOUBLEVECTOR *T1;
	DOUBLEVECTOR *Tstar;/*max temperature at which the first particle of ice comes up*/
	DOUBLEVECTOR *THETA;/*water content to be modified at each iteration*/
	DOUBLEVECTOR *soil_evap_layer_bare;/*Soil water evaporates from the first %ld layers*/
	DOUBLEVECTOR *soil_evap_layer_veg;/*Soil water transpires from the first %ld layers*/
	
} ENERGY;


/*---------------------------------------------------------------------------*/
typedef struct {

	LONGMATRIX *type;
	DOUBLETENSOR *pa;/** contains the soil parameters from the __soil.txt, pa[SoilType][Property][Layer] */
	DOUBLETENSOR *P;/**Psi,P[layer][Row][Cols]*/
	DOUBLETENSOR *Ptot;/**Psi*/
	DOUBLETENSOR *T;	/**Temperature*/
	DOUBLETENSOR *thice;/**Volume metric ice fraction*/
	DOUBLETENSOR *th;	/**Volume metric water fraction*/
	DOUBLEMATRIX *Jinf;
	DOUBLEMATRIX *Tv;/**vegetation temperature*/
	SHORTMATRIX *bc;
	DOUBLETENSOR *ET;
	
	DOUBLEMATRIX *th_av;
	DOUBLEMATRIX *thice_av;
	DOUBLEMATRIX *T_av;

} SOIL;


/*---------------------------------------------------------------------------*/
typedef struct {
    DOUBLEMATRIX *Z0;         /*elevetions of each pixel (DEM)*/
	DOUBLEMATRIX *Z1;
	DOUBLEMATRIX *Z0dp;		  /*DEM depitted*/
	DOUBLEMATRIX *Z0ext;      /*DEM extended (to avoid curvature problems)*/
	DOUBLETENSOR *Z;		  /*土壤各层中心（包括地表）的相对高度位置，绝对高度-流域平均海拔*/
    DOUBLEMATRIX *sky;        /*view factor (of the sky) for each pixel*/
    SHORTMATRIX *pixel_type; 
    
	SHORTMATRIX *DD;         /*Drainage Directions for each pixel; ex matr_ev->slopes*/
	LONGMATRIX *DDup;
	LONGVECTOR *DDdown;
	
    DOUBLEMATRIX *i_DD;       /*slope along Drainage Direction for each pixel*/
    DOUBLEMATRIX *top_index; /*topographic index for water content distribution*/
    DOUBLEMATRIX *area;       /*area of a pixel considering the slopeex matr_ev->area*/
    DOUBLEMATRIX *aspect;     /*aspect; ex: matr_ev->azimuth*/
    DOUBLEMATRIX *slopes;     /*slope of the pixels; ex: matr_ev->slopes*/
	double ****horizon_height;
	
	//for micromet
	DOUBLEMATRIX *Zm;/**/
	DOUBLEMATRIX *curv_m;/*curvature*/
	DOUBLEMATRIX *slope_m;/*terrain_slope*/
	DOUBLEMATRIX *slopeaz_m;/* Compute the slope azimuth, making sure that north has zero azimuth.
	Also note that for the Ryan wind rotation, the azimuth values must range from 0 to 360.*/
	
	//Cont for Richards 3D/
	/*Layer 0 is water on the surface. This allows a more robust description of infiltration processes.
	However, surface water flow is described separately in the subroutine supflow
	Layers > 0 , represent subsurface*/
	long ***i_cont;	/**node number for each active node,icont[Nl+1][Nr][Nc] 参考input.c i_lrc_cont()*/
	LONGMATRIX *lrc_cont;	/** lrc_cont[(Nl+1)*par->total_pixel + cnet->r->nh][3]*/
							/** lrc_cont[i_cont][1]=l*/
							/** lrc_cont[i_cont][2]=r*/
							/** lrc_cont[i_cont][3]=c*/
	//Cont for overland flow
	long **j_cont;
	LONGMATRIX *rc_cont;	/**rc_cont->co[cont][1]=r;
							   rc_cont->co[cont][2]=c;*/

		
	LONGVECTOR *Lp;
	LONGVECTOR *Li;
	LONGVECTOR *Up;
	LONGVECTOR *Ui;
			
} TOPO;


/*---------------------------------------------------------------------------*/
typedef struct {
    DOUBLEMATRIX *LC;            //land cover map (wood,lake,town,...) for each pixel*/
	//SHORTMATRIX *LC2;			  //FURTHER LAND CLASSIFICATION 	
    DOUBLEMATRIX *albedo;         //albedo calculated for each pixel*/

    /**Calculated in void shadow_haiden(short point, TOPO *top, double alpha, double direction, SHORTMATRIX *shadow)*/
	SHORTMATRIX *shadow;		  //shadows matrix (1 shadow, 0 sun)*/
	LONGMATRIX *cont;
	DOUBLEMATRIX *ty;/**各种土地利用类型的参数表,行方向:土地利用类型;列方向:23种参数值,ty[lu][1]对应参数文件的第二行*/
					 /**ty[n_landuses][nlandprop]参考__parameters.txt文件*/
	
	/**vegparp,vegpars,vegparv,vegpar used to read vegetation specified parameter files*/
	LONGMATRIX *vegparp;
	double ***vegpars;
	double **vegparv;
	DOUBLEVECTOR *vegpar;/**(r,c)点的植被参数,parameter文件植被块从第4行起*/
	
	DOUBLEMATRIX *root_fraction;/**每一种土地利用类型的各层的根分布,root_fraction[n_landuses][layer]*/
								/**weighting factor assigned to each layer for transpiration,
								 * sum of all root_fraction components gives 1
								 */
	
} LAND;/*all this data are calculated on the basis of land use data and some other par*/


/*---------------------------------------------------------------------------*/
typedef struct {/*nch=number of channel-pixel,ns=number of virtual stretches of channel,L=number of layers,
                  R=number of rows of the basin,C=number of columns in the basin*/
    LONGVECTOR *r;          /*array of rows of the channel-pixels; dimension=nch*/
    LONGVECTOR *c;          /*array of columns of the channel-pixels; dimension=nch*/
	LONGMATRIX *ch;			/*保存在搜索标记河道网格时每个网格被搜索的次序*/

	//input.c Initialization of the matrix with the sub/superficial flows which goes in a channel-cell ("cnet->q_sup"):*/
	DOUBLEVECTOR *Qsup;		/**？？？*/
	DOUBLEVECTOR *Qsub;		/**？？？*/
    //DOUBLEVECTOR *s0;       /*distance of each channel-pixel from outlet; dimension=nch*/
    //DOUBLEMATRIX *fraction_spread;/*fraction of flow for each s(=virtual distance) of all channel-pixel*/
    //DOUBLEVECTOR *Q_sup_s;   /*derived from q_sup[mc/s]; dimension=ns*/
    //DOUBLEVECTOR *Q_sub_s;   /*derived from q_sub[mc/s]; dimension=ns*/
	//DOUBLEVECTOR *Qsup_spread;
    //DOUBLEVECTOR *Qsub_spread;
	
	DOUBLEVECTOR *h_sup;	/**？？？*/
	DOUBLEVECTOR *dh_sup;	/**？？？*/
	
	DOUBLEMATRIX *hsupav;	/**？？？*/
	
	DOUBLEVECTOR *length;/**是河网的每一个网格内河段的长度*/
	
	double Q_out;			/**？？？*/
			
} CHANNEL;




/*---------------------------------------------------------------------------*/
typedef struct { /*nstations=number of all the rain-stations,npixel=number of all the pixels of the basin R*C,
                   R=number of rows,C=number of columns,nt=number of time-step of the whole similation*/
    DOUBLEMATRIX *weights_Kriging; /*dimension=npixel*nstations Kriging weights to model the rainfall distribution*/
    DOUBLEMATRIX *PrecTot;    /*total(snow+rain) precipitation in mm (in a Dt)*/
    DOUBLEMATRIX *Pnet;       /*liquid precipitation which reaches the sl surface in mm in a Dt as input
                              of "punctual_energy" subroutine, rain intensity in mm/s as output of the
                              same subroutine and in "water.balance.c" module*/
    DOUBLEMATRIX *wcan_rain;       /*liquid precipitation intercepted by vegetation in mm*/
    DOUBLEMATRIX *wcan_snow;       /*snow intercepted by vegetation in mm*/

    DOUBLEMATRIX *PrTOT_mean;  /*Total precipitation [mm](on nDt_output_basin Dt time intervals) total precipitation (rain+snow) precipitation*/
    DOUBLEMATRIX *PrSNW_mean; /*output of total precipitation and interception*/

	DOUBLEMATRIX *hsupav;
	DOUBLEMATRIX *dh_sup;
		
	DOUBLEMATRIX *error;
		
	//UMFPACK_REAL_TRIPLET *Jtriplet;
	//UMFPACK_REAL_MATRIX *Jmatrix;	
	
	DOUBLEVECTOR *Lx;
	DOUBLEVECTOR *Ux;
	
	DOUBLEVECTOR *H0;/*potential head for each active node in the basin at the beginning of iteration, vec length is Nl*Nr*Nc */
	DOUBLEVECTOR *H1;/*potential head for each active node in the basin, vec length is Nl*Nr*Nc */
	DOUBLEVECTOR *dH;
	DOUBLEVECTOR *B;
	DOUBLEVECTOR *f;
	DOUBLEVECTOR *df;
	
		
} WATER;


/*---------------------------------------------------------------------------*/
typedef struct {
    short iter;    /*current iteration of a time interval for egy-mass balance*/
    short n_iter;  /*n_iter=number of iterations to do for each time-step*/
    double TH;     /*TH=The numbers of days of simulation*/
    long i_pixel;  /*counter for the output of a pixel*/
    long n_pixel;  /*nDt_output_pixel=number of Dt after which the output of a pixel are printed*/
	long i_basin;  /*counter for the output of a pixel*/
	long n_basin;  /*THE NUMBER OF Dt AFTER WHICH THE OUTPUT FOR THE BASIN ARE PRINTED*/
	long i_discharge;/*counter for the output of a pixel*/
	long n_discharge;/*THE NUMBER OF Dt AFTER WHICH THE DISCHARGE IS PRINTED*/
	long i_plot;
	long n_plot;
	long nt_plot;
	long d_plot;
    double JD;      /*day=current Julian day during the simulation*/
    long day;       /*current day of the month during the simulation*/
    long month;       /*current month of the year during the simulation*/
    long year;     /*current year*/
    long hour;       /*current hour of day during the simulation*/
    long min;       /*current minute of hour during the simulation*/
    double time;    /*time=current time from the begin of simulation [s]*/
    double egy;  /*the time of egy subroutine [s]*/
    double vert_wb; /*the time of water-vertical-distribution subroutine [s]*/
    double horiz_wb;/*the time of water-horizontal-distribution subroutine [s]*/
    double writeout;/*the time of write-output subroutine [s]*/
} TIMES;


/*---------------------------------------------------------------------------*/
typedef struct {
    double Dt;      /*Dt=the integration time interval [s]*/
	double JD0;		/*Decimal julian day at the beginning of simulation (0.0-365.99)*/
	long year0;		/*Year at the beginning of simulation*/
	double ST;		/*Standard time to which all the output data are referred (difference respect UMT, in hour)*/
    short print;         /*1 IF YOU WANT TO PRINT MATRICES WITH INTERMEDIATE RESULTS, 0 OTHERWISE*/
    short monin_obukhov;  /**#2      monin_obuhkov    1 stability and instability considered, 2 stability not considered, 3 instability not considered, 4 always neutrality,
			 5 Oerlemans&Grisogono (2002) on glacier*/

    double gamma_m;   /**5 Exponent of the law of uniform motion on the surface*/
    double T_rain;    /*TEMPERATURE ABOVE WICH ALL PRECIPITAION IS RAIN [C]*/
    double T_snow;    /*TEMPERATURE BELOW WICH ALL PRECIPITAION IS SNOW [C]*/
    double aep;       /*ALBEDO EXTINCTION PARAMETER [m] if snow depth < aep, albedo is interpolated between soil and snow*/
    double avo;       /*NEW SNOW VISIBLE BAND REFLECTANCE*/
    double airo;      /*NEW NEAR INFRARED BAND REFLECTANCE*/
    double Sr;		  /*WATER FRACTION RETAINED BY CAPILLARY FORCES IN SNOW IRREDUCIBLE WATER SATURATION [-] - from Colbeck (0.02 - 0.07)*/
    double rho_ice;     /*Ice density [kg/mc]*/
    long total_pixel;    /*The number of the valid pixel of the whole basin*/
	long total_channel;	/**The total number of channel-pixel in the basin */
 	long snowlayer_max; /*MAXIMUM NUMBER OF SNOW LAYERS*/
	long snowlayer_inf; /*layer of illimitate thickness (beginning from below)*/
	DOUBLEVECTOR *Dmin; /*MINIMUM SNOW LAYER THICKNESS*/
	DOUBLEVECTOR *Dmax; /*MAXIMUM SNOW LAYER THICKNESS*/
	double Sr_glac;
	long glaclayer_max;
	DOUBLEVECTOR *Dmin_glac;
	DOUBLEVECTOR *Dmax_glac;

	short state_turb;

	/** 9 block - BASE AND ADVANCED PARAMETERS
	#1      state_lwrad	 Which formula for incoming longwave radiation:
		1 (Brutsaert, 1975), 2 (Satterlund, 1979), 3 (Idso, 1981), 4(Idso+Hodges),
		5 (Koenig-Langlo & Augstein, 1994), 6 (Andreas & Ackley, 1982), 7 (Konzelmann, 1994),
		8 (Prata, 1996), 9 (Dilley 1998)*/
	short state_lwrad;

	double imp;/**1 Impedence factor for (partially) frozen soil  1.Reduction factor of the hydraulic conductivity in partially frozen soil (K=K_no_ice*10^(impedence*Q), where Q is the ice ratio)
	suggested values: from 0 to 7, best value 4*/
	double f_bound_Richards;

	double epsilon_snow;/**SNOW LONGWAVE EMISSIVITY [-]*/
	
	double output_Txy;
	double output_TETAxy;
	double output_TETAICExy;
	double output_PSIxy;
	double output_snow;
	double output_glac;
	double output_h_sup;
    double output_Rn;
	double output_G;
	double output_H;
	double output_ET;
	double output_Ts;
	double output_P;
	double output_Wr;
	double output_balancesn;
	double output_balancegl;
	double output_Rswdown;
	double output_meteo;
	



	/**
	 * chkpt COORDINATES of the points for which the simulations is printed
	 * rc Rows and Cols of the points for which the simulations is printed
	 */
	DOUBLEMATRIX *chkpt;
	LONGMATRIX *rc;
	/**
	 * state_px_coord	 1 IF ALL COORDINATES ARE IN FORMAT (EAST,NORTH), 0 IF IN FORMAT ROW AND COLUMS (r,c)
	 */
	short state_px_coord;
	
	double integr_scale_rain;
	double variance_rain;
	
	short recover;		/*1 if you want to recover a simulation, 0 otherwise*/

	double Vmin;/**MINIMUM WIND VELOCITY (m/s)*/
		
	double snowcorrfact;/*factor multiplying the snow precipitation*/
	double raincorrfact;/*factor multiplying the rain precipitation*/
	
	double RHmin;/**MINIMUM RELATIVE HUMIDITY VALUE (%)*/
	
	short format_out;	/*OUTPUT MAPS in fluidturtle format (=1), GRASS ASCII (=2), ESRI ASCII (=3)*/
	short nsky;			/*Multiplying factor decreasing the dem resolution for the calculation of the sky view factor*/
	double channel_thres;/*Value of the threshold for definition of pixel channel [1.0-1000.0] (used if the network map is not provided)*/
	
	/**
	 * saving_points saving points expressed as numbers of day after the beginning
	 */
	DOUBLEVECTOR *saving_points;
	
	//double Vis; //visibility in km (>5 km)
	//double Lozone; //thickness of the stratospheric ozone layer (in cm normal conditions)
	
	short point_sim;/* =0 distributed simulation, =1 point simulation (the parameter files are different in the two cases) */
			
	double snow_maxpor;/*MAXIMUM SNOW POROSITY ALLOWED*/
	double snow_density_cutoff;
	double drysnowdef_rate;/*SNOW COMPACTION (% per hour) DUE TO DESTRUCTIVE METAMORPHISM for SNOW DENSITY<snow_density_cutoff and DRY SNOW */
	double wetsnowdef_rate;/*ENHANCEMENT FACTOR IN PRESENCE OF WET SNOW*/
	double snow_viscosity; /*SNOW VISCOSITY COEFFICIENT (kg s m^-2) at T=0 C and snow density=0*/
	
	double latitude;/**latitude in radiance*/
	double longitude;
	
	double z0_snow;/*roughness length over snow (mm)*/
	long n_landuses;/**The Maximum number of land use types*/
	
	/**JD_plots It is possible to select some days for which energy balance and
	 * meteo data are plotted with a very short time step.
	 */
	DOUBLEVECTOR *JD_plots;
	
	short micromet;/**Use Micromet (=1), otherwise (=0)*/
	short blowing_snow;/**PBSM	 Use PBSM (snow horizontal transport) (=1), otherwise (=0)*/
	
	LONGVECTOR *r_points;
	LONGVECTOR *c_points;
	
	double psimin;
	double stmin;
		    	
	/**wat_balance	water balance 	 1 calculate, 0 do not calculate
	 * en_balance	energy balance   1 calculate, 0 do not calculate
	 */
	short wat_balance;
	short en_balance;
	
	short distr_stat;
	
	double Dpsi;
	double dtmin;
	
	double PsiInf;
	double TolPsiInf;
	
	double q1;
	double q2;

	double ***transect;
	double **vtrans;
	
	LONGVECTOR *cont_trans;
	LONGVECTOR *ibeg;
	
	long nLC;	
	
	double snow_fetch;/*MINIMUM HORIZONTAL FETCH THAT ALLOWS FULL SNOW TRANSPORT*/
	
	short ifill;
	short iobsint;
	double dn;
	double curve_len_scale;/**__parameter.txt block 3 16*/
	double slopewt;
	double curvewt;	
	short topoflag;/**if=1 top->Zm->co[r][c] = top->Z0->co[r][c] + DEPTH(r, c, snow->lnum, snow->Dzl);*/
	
	short LRflag;/**LRflag=1,has a lapse rates file;LRflag=0,without a lapse rates file*/
	
	SHORTVECTOR *vegflag;/**par->vegflag->co[i]=1,There is a specific vegetation parameter file for land cover type i;
						    par->vegflag->co[i]=0,There is NOT a specific vegetation parameter file for land cover type i*/
	
	short harm_or_arit_mean;
	long MaxiterTol;/**Max iterations for the integration of Richards' equation*/
	double MaxErrWb;/**MaximalErrorRichards*/
	double TolVWb;/**2  Absolute Tolerance for the integration of Richards' equation*/
	
	double incr_V_factor;/*Factor multiplying the averaged wind speed in order to run the blowing snow model*/

	double alpha_snow;
	double tol_energy;
	long maxiter_energy;
	long maxiter_canopy;
	long maxiter_Ts;
	long maxiter_Loc;
	long maxiter_Businger;
	short stabcorr_incanopy;
	
	short state_pixel;/* if n_pixel>0, state_pixel=1,else state_pixel=0 */
	short state_basin;/* if n_basin>0, state_basin=1,else state_basin=0 */
	short state_discharge;/* if n_discharge>0, state_discharge=1,else state_discharge=0 */
			
	double DtminWb;
	long nredDtWb;
	
	double TolCG;/**Initial forcing term of Newton method */
	long MaxiterCorr;
	short UpdateK;
	
	short channel_network;
	short bedrock;/** indicate whether the bedrock file is input,Yes=1,No=0 */
	
	double thres_hsup;/**Threshold on h_sup [mm] above which C_m is independent from h_sup*/
	double thres_hchannel;
	
	//double cm_hsup_0;
	//double max_Courant_sup_sub;
	
	double Kch_b;/**hydraulic conductivity of the sediments divided by sediment thickness*/
	double w_dx;/**channel width / pixel width*/
	
	double RelTolVWb;/**RelativeErrorRichards*/
	
	double snow_smin;/*MINIMUM SLOPE [degree] TO ADJUST PRECIPITATION REDUCTION*/
	double snow_smax;/*MAXIMUM SLOPE [degree] TO ADJUST PRECIPITATION REDUCTION*/
	double snow_curv;/*SHAPE PARAMETER FOR PRECIPITATION REDUCTION,correction of snow precipitation in case of steep slopes (contribution by Stephan Gruber)*/
	
	double Zboundary;/**0 temperature amplitude depth in the soil [mm]*/
	double Tboundary;/**temperature at the 0 temperature amplitude depth [C]*/
	
	double Ks_channel;/**Cm coefficient for the channel flow (the same gamma for surface flow is used)*/
	double depr_channel;/**Depression of the channel bed with respect to the neighbouring slopes*/
	
} PAR;



/*---------------------------------------------------------------------------*/
typedef struct { 
	SHORTMATRIX  *type;/**/
	LONGMATRIX	 *lnum;/*snow layer number*/
	DOUBLETENSOR *Dzl;/**snow thickness for each layer*/
	DOUBLETENSOR *w_liq;/**liquid water content for each snow layer*/
	DOUBLETENSOR *w_ice;/**ice content for each snow layer*/
	DOUBLETENSOR *T;/*Temperature for each snow layer*/
	DOUBLEMATRIX *nondimens_age;
	DOUBLEMATRIX *dimens_age;	
	DOUBLEMATRIX *max;
	DOUBLEMATRIX *average;	
	DOUBLEMATRIX *MELTED;
	DOUBLEMATRIX *SUBL;
	DOUBLEMATRIX *t_snow;
	
	DOUBLEMATRIX *DDF;
	DOUBLEMATRIX *DDF1;
	DOUBLEMATRIX *DDFvar;	
	LONGMATRIX   *DDFcont;
	DOUBLEMATRIX *DDFTmin;
	DOUBLEMATRIX *DDFmeltTL0;
	DOUBLEMATRIX *DDFmelt;
	
	DOUBLEMATRIX *rho_newsnow;
	DOUBLEMATRIX *Qsub;
	DOUBLEMATRIX *Wtrans;
	DOUBLEMATRIX *Qtrans;
	DOUBLEMATRIX *Qtrans_x;
	DOUBLEMATRIX *Qtrans_y;	
	DOUBLEMATRIX *Wtot;
	DOUBLEMATRIX *Wsubl_cum;
	DOUBLEMATRIX *Wsusp_cum;
	DOUBLEMATRIX *Wtrans_cum;
	DOUBLEMATRIX *Wsubgrid_cum;
	DOUBLEMATRIX *out_bs;
	DOUBLEMATRIX *ListonSWE;
	DOUBLEMATRIX *softSWE;
	DOUBLEMATRIX *softSWE1;	
	DOUBLEMATRIX *Dplot;	
	
	DOUBLEVECTOR *CR1;
	DOUBLEVECTOR *CR2;
	DOUBLEVECTOR *CR3;
	
	LONGVECTOR *change_dir_wind;
	DOUBLEMATRIX *Psnow;/* The total snow fall to the ground,Dripp + through */
	
} SNOW;

typedef struct { 
	LONGMATRIX	 *lnum;
	DOUBLETENSOR *Dzl;
	DOUBLETENSOR *w_liq;
	DOUBLETENSOR *w_ice;
	DOUBLETENSOR *T;	
	DOUBLEMATRIX *MELTED;
	DOUBLEMATRIX *SUBL;	
} GLACIER;

typedef struct{
	DOUBLEVECTOR *E;		/**East(m)*/
	DOUBLEVECTOR *N;		/**North(m)*/
	DOUBLEVECTOR *lat;		/**Latitude(degree)*/
	DOUBLEVECTOR *lon;		/**Longitude(degree)*/
	DOUBLEVECTOR *Z;		/**Elevation (m a.s.l.)*/
	DOUBLEVECTOR *sky;		/**Sky view factor (-)*/
	DOUBLEVECTOR *ST;		/**Standard time minus UTM for the data*/
	DOUBLEVECTOR *Vheight;	/**Wind velocity measurement height (m a.g.l.)*/
	DOUBLEVECTOR *Theight;	/**Temperature measurement height (m a.g.l.)*/
	DOUBLEVECTOR *JD0;		/**Decimal Julian Day of the first data*/
	LONGVECTOR *Y0;			/**Year of the first data*/
	DOUBLEVECTOR *Dt;		/**Dt data (hour)*/
	LONGVECTOR *offset;		/**offset column*/
} METEO_STATIONS;


typedef struct {
	METEO_STATIONS *st;

	/** data[stationNum][time][properties]*/
	double ***data;/**Data read from the __meteo0001.txt file,keeping the column order*/
	long **column;/** Row:the number of stations,Cols:the columns of the meteo input file
										 start from 0 用来保存每一个站点的驱动变量在输入文件中的列数*/
	double ***horizon;
	double **var;/*meteo variables for the current instant 13 var[stationNum][properties]*/
				/**start from 0 var中保存了模型模拟的某一时刻时所有气象占的驱动数据 行方向为气象行占，列方向为气象变量*/
	double **LRs;
	double *LRv;
	LONGVECTOR *LRp;/**LAPSE RATES*/
	
	DOUBLEMATRIX *Tgrid;/**模型运行到某一时刻的驱动数据的空间分布，可有micromet插值得到，也可直接输入模型*/
	DOUBLEMATRIX *Pgrid;
	DOUBLEMATRIX *Vgrid;
	DOUBLEMATRIX *Vdir;
	DOUBLEMATRIX *RHgrid;
	
	DOUBLEMATRIX *Vspdmean;
	DOUBLEMATRIX *Vdirmean;
	DOUBLEMATRIX *RHmean;
	
	DOUBLEMATRIX *Taplot;
	DOUBLEMATRIX *Vspdplot;
	DOUBLEMATRIX *Vdirplot;
	DOUBLEMATRIX *RHplot;
	
	double V;
	double RH;	
		
	DOUBLEMATRIX *Tday;	
	DOUBLEMATRIX *Tvar;
		
	long nstsrad;/**nstsrad找到的具有段波辐射的第一个气象占的序号 1~气象站数目*/
	long nstlrad;/**nstsrad找到的具有长波辐射的第一个气象占的序号 1~气象站数目*/
	long nstcloud;/**nstsrad找到的具有云数据的第一个气象占的序号 1~气象站数目*/
	
} METEO;


typedef struct {
	SOIL *S;
	WATER *W;
	LAND *L;
	PAR *P;
	TOPO *T;
	CHANNEL *C;
	ENERGY *E;
	SNOW *N;
	GLACIER *G;
	METEO *M;
	TIMES *I;	
}ALLDATA;
