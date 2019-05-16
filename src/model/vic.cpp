#include "vic.h"
using namespace ldas;
using namespace std;

Vic::Vic()
{
    //ctor
    init();
}

Vic::~Vic()
{
    //dtor
}

Vic::Vic(const Vic& other)
{
    //copy ctor
    init();
}

Vic& Vic::operator=(const Vic& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}

void Vic::init()
{
//null now
    version = "4.1.1";

    optstring = "g:vo";

#if QUICK_FS
    temps[0]=-1.e-5;
    temps[1]=-0.075;
    temps[2]=-0.20;
    temps[3]=-0.50;
    temps[4]=-1.00;
    temps[5]=-2.50;
    temps[6]=-5;
    temps[7]=-10;
#endif


    /**************************************************************************
      Define some reference landcover types that always exist regardless
      of the contents of the library (mainly for potential evap calculations):
      Non-natural:
        satsoil = saturated bare soil
        h2osurf = open water surface (deep enough to have albedo of 0.08)
        short   = short reference crop (grass)
        tall    = tall reference crop (alfalfa)
      Natural:
        natveg  = current vegetation
        vegnocr = current vegetation with canopy resistance set to 0
      NOTE: these are external variables, declared in vicNl_def.h.
      NOTE2: bare soil roughness and displacement will be overwritten by the
             values found in the soil parameter file; bare soil wind_h will
       be overwritten by the value specified in the global param file.
    **************************************************************************/

    /* One element for each non-natural PET type */
    for(int i=0; i<4; i++)
    {
        ref_veg_over[i]=0;
        ref_veg_rad_atten[i]=0;
        ref_veg_wind_atten[i]=0;
        ref_veg_trunk_ratio[i]=0;
    }
    ref_veg_rarc[0]=0;
    ref_veg_rarc[1]=0;
    ref_veg_rarc[2]=25;
    ref_veg_rarc[3]=25;

    ref_veg_rmin[0]=0;
    ref_veg_rmin[1]=0;
    ref_veg_rmin[2]=100;
    ref_veg_rmin[3]=100;

    ref_veg_lai[0]=1;
    ref_veg_lai[1]=1;
    ref_veg_lai[2]=2.88;
    ref_veg_lai[3]=4.45;

    ref_veg_albedo[0]=BARE_SOIL_ALBEDO;
    ref_veg_albedo[1]=H2O_SURF_ALBEDO;
    ref_veg_albedo[2]=0.23;
    ref_veg_albedo[3]=0.23;

    ref_veg_rough[0]=0.001;
    ref_veg_rough[1]=0.001;
    ref_veg_rough[2]=0.0148;
    ref_veg_rough[3]=0.0615;

    ref_veg_displ[0]=0.0054;
    ref_veg_displ[1]=0.0054;
    ref_veg_displ[2]=0.08;
    ref_veg_displ[3]=0.3333;

    ref_veg_wind_h[0]=10;
    ref_veg_wind_h[1]=10;
    ref_veg_wind_h[2]=10;
    ref_veg_wind_h[3]=10;

    ref_veg_RGL[0]=0;
    ref_veg_RGL[1]=0;
    ref_veg_RGL[2]=100;
    ref_veg_RGL[3]=100;
    /* One element for each PET type (non-natural or natural) */
    ref_veg_ref_crop[0]=FALSE;
    ref_veg_ref_crop[1]=FALSE;
    ref_veg_ref_crop[2]=TRUE;
    ref_veg_ref_crop[3]=TRUE;
    ref_veg_ref_crop[4]=FALSE;
    ref_veg_ref_crop[5]=FALSE;
    initialize_global();
}

void Vic::config(const char* fn)
{
    strcpy(filenames.global, fn);
}

void Vic::run()
{
    /** Read Model Options **/
    //filenames = cmd_proc(argc, argv);

#if VERBOSE
    display_current_settings(DISP_VERSION,(filenames_struct*)NULL,(global_param_struct*)NULL);
#endif
    /** Read Global Control File **/
    filep.globalparam = open_file(filenames.global,"r");
    global_param = get_global_param(&filenames, filep.globalparam);
    /** Set up output data structures **/
    out_data = create_output_list();
    out_data_files = set_output_defaults(out_data);
    fclose(filep.globalparam);
    filep.globalparam = open_file(filenames.global,"r");
    parse_output_info(&filenames, filep.globalparam, &out_data_files, out_data);

    /** Check and Open Files **/
    check_files(&filep, &filenames);

#if !OUTPUT_FORCE

    /** Check and Open Debugging Files **/
#if LINK_DEBUG
    open_debug();
#endif

    /** Read Vegetation Library File **/
    veg_lib = read_veglib(filep.veglib,&Nveg_type);

#endif // !OUTPUT_FORCE

    /** Initialize Parameters **/
    if(options.DIST_PRCP) Ndist = 2;
    else Ndist = 1;
    cellnum = -1;

    /** Make Date Data Structure **/
    dmy      = make_dmy(&global_param);

    /** allocate memory for the atmos_data_struct **/
    alloc_atmos(global_param.nrecs, &atmos);

    /** Initial state **/
    startrec = 0;
#if !OUTPUT_FORCE
    if ( options.INIT_STATE )
        filep.init_state = check_state_file(filenames.init_state, dmy,
                                            &global_param, options.Nlayer,
                                            options.Nnode, &startrec);

    /** open state file if model state is to be saved **/
    if ( options.SAVE_STATE && strcmp( filenames.statefile, "NONE" ) != 0 )
        filep.statefile = open_state_file(&global_param, filenames, options.Nlayer,
                                          options.Nnode);
    else filep.statefile = NULL;

#endif // !OUTPUT_FORCE

    /************************************
      Run Model for all Active Grid Cells
      ************************************/
    MODEL_DONE = FALSE;
    cell_cnt=0;
    int flag;
    while(!MODEL_DONE)
    {
        if(!options.ARC_SOIL)
        {
            if((fscanf(filep.soilparam, "%d", &flag))!=EOF)
            {
                if(flag) RUN_MODEL=TRUE;
                else     RUN_MODEL=FALSE;
            }
            else
            {
                MODEL_DONE = TRUE;
                RUN_MODEL = FALSE;
            }
            if(!MODEL_DONE) soil_con = read_soilparam(filep.soilparam, RUN_MODEL);

        }
        else
        {
            soil_con = read_soilparam_arc(filep.soilparam,
                                          filenames.soil_dir, &Ncells,
                                          &RUN_MODEL, cell_cnt);
            cell_cnt++;
            if(cell_cnt==Ncells) MODEL_DONE = TRUE;
        }
        if(RUN_MODEL)
        {
#if LINK_DEBUG
            if(debug.PRT_SOIL) write_soilparam(&soil_con);
#endif

#if QUICK_FS
            /** Allocate Unfrozen Water Content Table **/
            if(options.FROZEN_SOIL)
            {
                for(i=0; i<MAX_LAYERS; i++)
                {
                    soil_con.ufwc_table_layer[i] = (double **)malloc((QUICK_FS_TEMPS+1)*sizeof(double *));
                    for(j=0; j<QUICK_FS_TEMPS+1; j++)
                        soil_con.ufwc_table_layer[i][j] = (double *)malloc(2*sizeof(double));
                }
                for(i=0; i<MAX_NODES; i++)
                {
                    soil_con.ufwc_table_node[i] = (double **)malloc((QUICK_FS_TEMPS+1)*sizeof(double *));

                    for(j=0; j<QUICK_FS_TEMPS+1; j++)
                        soil_con.ufwc_table_node[i][j] = (double *)malloc(2*sizeof(double));
                }
            }
#endif /* QUICK_FS */

            NEWCELL=TRUE;
            cellnum++;

#if !OUTPUT_FORCE

            /** Read Grid Cell Vegetation Parameters **/
            veg_con = read_vegparam(filep.vegparam, soil_con.gridcel,
                                    Nveg_type);
            calc_root_fractions(veg_con, &soil_con);
#if LINK_DEBUG
            if(debug.PRT_VEGE) write_vegparam(veg_con);
#endif /* LINK_DEBUG*/

            if ( options.LAKES )
                lake_con = read_lakeparam(filep.lakeparam, soil_con, veg_con);

#endif // !OUTPUT_FORCE

            /** Build Gridded Filenames, and Open **/
            make_in_and_outfiles(&filep, &filenames, &soil_con, out_data_files);

            if (options.PRT_HEADER)
            {
                /** Write output file headers **/
                write_header(out_data_files, out_data, dmy, global_param);
            }

#if !OUTPUT_FORCE

            /** Read Elevation Band Data if Used **/
            read_snowband(filep.snowband, &soil_con);

            /** Make Precipitation Distribution Control Structure **/
            prcp     = make_dist_prcp(veg_con[0].vegetat_type_num);

#endif // !OUTPUT_FORCE

            /**************************************************
               Initialize Meteological Forcing Values That
               Have not Been Specifically Set
             **************************************************/

#if VERBOSE
            fprintf(stderr,"Initializing Forcing Data\n");
#endif /* VERBOSE */

            initialize_atmos(atmos, dmy, filep.forcing,
                             (double)soil_con.time_zone_lng, (double)soil_con.lng,
                             (double)soil_con.lat, soil_con.elevation,
                             soil_con.annual_prec, global_param.wind_h,
                             soil_con.rough, soil_con.avgJulyAirTemp,
                             soil_con.Tfactor,
#if OUTPUT_FORCE
                             soil_con.AboveTreeLine, out_data_files, out_data);
#else /* OUTPUT_FORCE */
                             soil_con.AboveTreeLine);
#endif /* OUTPUT_FORCE */
            printf("after Initialize_atoms");
#if !OUTPUT_FORCE
#if LINK_DEBUG
            if(debug.PRT_ATMOS) write_atmosdata(atmos, global_param.nrecs);
#endif

            /**************************************************
              Initialize Energy Balance and Snow Variables
            **************************************************/

#if VERBOSE
            fprintf(stderr,"Model State Initialization\n");
#endif /* VERBOSE */
            ErrorFlag = initialize_model_state(&prcp, dmy[0], &global_param, filep,
                                               soil_con.gridcel, veg_con[0].vegetat_type_num,
                                               options.Nnode, Ndist,
                                               atmos[0].air_temp[NR],
                                               &soil_con, veg_con, lake_con,
                                               &init_STILL_STORM, &init_DRY_TIME, &save_data);
            if ( ErrorFlag == ERROR )
            {
                if ( options.CONTINUEONERROR == TRUE )
                {
                    // Handle grid cell solution error
                    fprintf(stderr, "ERROR: Grid cell %i failed in record %i so the simulation has not finished.  An incomplete output file has been generated, check your inputs before rerunning the simulation.\n", soil_con.gridcel, rec);
                    break;
                }
                else
                {
                    // Else exit program on cell solution error as in previous versions
                    sprintf(ErrStr, "ERROR: Grid cell %i failed in record %i so the simulation has ended. Check your inputs before rerunning the simulation.\n", soil_con.gridcel, rec);
                    vicerror(ErrStr);
                }
            }

#if VERBOSE
            fprintf(stderr,"Running Model\n");
#endif /* VERBOSE */

            /** Update Error Handling Structure **/
            Error.filep = filep;
            Error.out_data_files = out_data_files;

            /** Initialize the storage terms in the water and energy balances **/
            /** Sending a negative record number (-global_param.nrecs) to dist_prec() will accomplish this **/
            ErrorFlag = dist_prec(&atmos[0], &prcp, &soil_con, veg_con,
                                  &lake_con, dmy, &global_param, &filep, out_data_files,
                                  out_data, &save_data, -global_param.nrecs, cellnum,
                                  NEWCELL, LASTREC, init_STILL_STORM, init_DRY_TIME);

            /******************************************
            Run Model in Grid Cell for all Time Steps
            ******************************************/

            for ( rec = startrec ;
                    rec < global_param.nrecs;
                    rec++ )
            {

                if ( rec == global_param.nrecs - 1 ) LASTREC = TRUE;
                else LASTREC = FALSE;

                ErrorFlag = dist_prec(&atmos[rec], &prcp, &soil_con, veg_con,
                                      &lake_con, dmy, &global_param, &filep,
                                      out_data_files, out_data, &save_data, rec, cellnum,
                                      NEWCELL, LASTREC, init_STILL_STORM, init_DRY_TIME);

                if ( ErrorFlag == ERROR )
                {
                    if ( options.CONTINUEONERROR == TRUE )
                    {
                        // Handle grid cell solution error
                        fprintf(stderr, "ERROR: Grid cell %i failed in record %i so the simulation has not finished.  An incomplete output file has been generated, check your inputs before rerunning the simulation.\n", soil_con.gridcel, rec);
                        break;
                    }
                    else
                    {
                        // Else exit program on cell solution error as in previous versions
                        sprintf(ErrStr, "ERROR: Grid cell %i failed in record %i so the simulation has ended. Check your inputs before rerunning the simulation.\n", soil_con.gridcel, rec);
                        vicerror(ErrStr);
                    }
                }

                NEWCELL=FALSE;
                for ( veg = 0; veg <= veg_con[0].vegetat_type_num; veg++ )
                    init_DRY_TIME[veg] = -999;

            }	/* End Rec Loop */

#endif /* !OUTPUT_FORCE */

            close_files(&filep,out_data_files,&filenames);

#if !OUTPUT_FORCE

#if QUICK_FS
            if(options.FROZEN_SOIL)
            {
                for(i=0; i<MAX_LAYERS; i++)
                {
                    for(j=0; j<6; j++)
                        free((char *)soil_con.ufwc_table_layer[i][j]);
                    free((char *)soil_con.ufwc_table_layer[i]);
                }
                for(i=0; i<MAX_NODES; i++)
                {
                    for(j=0; j<6; j++)
                        free((char *)soil_con.ufwc_table_node[i][j]);
                    free((char *)soil_con.ufwc_table_node[i]);
                }
            }
#endif /* QUICK_FS */
            free_dist_prcp(&prcp,veg_con[0].vegetat_type_num);
            free_vegcon(&veg_con);
            free((char *)soil_con.AreaFract);
            free((char *)soil_con.BandElev);
            free((char *)soil_con.Tfactor);
            free((char *)soil_con.Pfactor);
            free((char *)soil_con.AboveTreeLine);
            free((char*)init_STILL_STORM);
            free((char*)init_DRY_TIME);
#endif /* !OUTPUT_FORCE */
        }	/* End Run Model Condition */
    } 	/* End Grid Loop */

    /** cleanup **/
    free_atmos(global_param.nrecs, &atmos);
    free_dmy(&dmy);
    free_out_data_files(&out_data_files);
    free_out_data(&out_data);
#if !OUTPUT_FORCE
    free_veglib(&veg_lib);
#endif /* !OUTPUT_FORCE */
    fclose(filep.soilparam);
#if !OUTPUT_FORCE
    fclose(filep.vegparam);
    fclose(filep.veglib);
    if (options.SNOW_BAND>1)
        fclose(filep.snowband);
    if (options.LAKES)
        fclose(filep.lakeparam);
#endif /* !OUTPUT_FORCE */
}
