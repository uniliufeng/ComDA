/******************************************************************************
// $Id: LAKE.h,v 5.9.2.12 2009/10/08 02:03:06 vicadmin Exp $
  Modifications:
  2006-Nov-07 Removed LAKE_MODEL option.				TJB
  2007-Aug-16 Made return value of initialize_prcp an int.		JCA
  2007-Aug-21 Added features for EXCESS_ICE option.			JCA
  2007-Oct-24 Changed get_sarea, get_volume, and get_depth to return exit
	      status so that errors can be trapped and communicated up the
	      chain of function calls.					KAC via TJB
  2007-Oct-24 Changed lakeice() to return exit status.			KAC via TJB
  2007-Nov-06 New lake physics parameters.  Modified argument lists for
	      various functions.  Moved get_dist() to vicNl.h.		LCB via TJB
  2008-Apr-21 Added argument to alblake.				LCB via TJB
  2008-Sep-10 Updated values of CONDS and lamwlw to match Laura
	      Bowling's lake work.					LCB via TJB
  2009-Jul-31 Removed lakemain() and wetland_energy(); initialize_lake
	      no longer takes a snow structure as input.		TJB
  2009-Sep-28 Removed initialize_prcp and update_prcp.  Modified
	      argument list of initialize_lake.				TJB
  2009-Sep-30 Miscellaneous fixes for lake model.			TJB
  2009-Oct-05 Added functions for updating/rescaling lake and wetland
	      fluxes and storages when lake area changes.		TJB
******************************************************************************/

//#ifndef LAKE_SET

#define LAKE_SET
#define TMELT 0.0
#define EMICE 0.97      /* Ice emissivity */
#define EMH2O .98
#define RHOSNOW   250.  /* densities of water and snow */
#define RHOICE   917.   /* ice density*/
#define rhosurf 1.275   /* surface air density */
#define MAX_SURFACE_LAKE   .6  /* max. surface layer thickness for E-B (m) */
#define BETA 0.001       /* Curve shape parameter for lake profile. */
#define FRACMIN  0.10   /* min ice thickness in meters */
#define FRACLIM   0.02  /* lower limit on fractional ice cover */
#define CPW_ICE   4200. /* specific heat of ice */
#define DM 1.38889E-07  /* molecular diffusivity of water */
#define SNOWCRIT   0.05  /* for albedo, in m */
//#define G 9.80616
#define ZWATER 0.0045    // 0.004 - original value
#define ZSNOW 0.005
#define CONDI 2.3        /* thermal conductivity of ice */
#define CONDS 0.7       /* thermal conductivity of snow */ 

// attenuation of short and longwave radiation through ice (1/m)
#define lamisw 1.5 // 1.5 in Patterson & Hamblin
#define lamilw 20  // 20.0 in Patterson & Hamblin
// attenuation of short and longwave radiation through snow (1/m)
#define lamssw 6.0 // 6.0 in Patterson & Hamblin
#define lamslw 20  // 20.0 in Patterson & Hamblin
// attenuation of short and longwave radiation through water (1/m)
#define lamwsw .3  // San Fran Bay data: 0.31 - 29.9 1/m (visible)
#define lamwlw 1.4 // Hostetler and Bartlein assume 0.85 1/m (total)
#define  a1 0.7        /* Percent of radiation in visible band. */
#define  a2 0.3        /* Percent of radiation in infrared band. */
#define QWTAU 86400./2.   /* D. Pollard sub-ice time constant. */
#define RADIUS 6371.228 /* Earth radius in km. */

//#endif // LAKE_SET

