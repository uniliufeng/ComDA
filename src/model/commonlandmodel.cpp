#include "commonlandmodel.h"

namespace ldas
{
CoLMConfig::CoLMConfig():ModelConfig(),row_num(0),col_num(0),sub_row_num(0),sub_col_num(0),sub_start_row_num(0),sub_start_col_num(0),ensemble_size(0)
{
}
void CoLMConfig::open(const char* fn)
{
    char buf[MAXBUFFER];
    ifstream in(fn);
    if(!in)
    {
        string error_msg="can't open config data file: ";
        error_msg+=fn;
        throw Exception("CoLMConfig","open(fn)",error_msg.c_str());
    }

    //skip the first 3 lines, for comments
    in.getline(buf,MAXBUFFER);
    in.getline(buf,MAXBUFFER);
    in.getline(buf,MAXBUFFER);

    //vegetation data path
    getline(in,vegetation_path,' ');
    cout<<vegetation_path<<endl;
    in.getline(buf,MAXBUFFER);

    //SSMI data path, for ldas
    getline(in,ssmi_path,' ');
    cout<<ssmi_path<<endl;
    in.getline(buf,MAXBUFFER);

    //AMSR data path
    getline(in,amsr_path,' ');
    cout<<amsr_path<<endl;
    in.getline(buf,MAXBUFFER);

    //meteorological data path
    getline(in,forcing_path,' ');
    cout<<forcing_path<<endl;
    in.getline(buf,MAXBUFFER);

    //output data path
    getline(in,output_path,' ');
    cout<<output_path<<endl;
    in.getline(buf,MAXBUFFER);

    //restart data path
    getline(in,restart_path,' ');
    cout<<restart_path<<endl;
    in.getline(buf,MAXBUFFER);

    in.getline(buf,MAXBUFFER);
    in.getline(buf,MAXBUFFER);
    in.getline(buf,MAXBUFFER);

    //file name of time-invariant file
    getline(in,constant_path,' ');
    cout<<constant_path<<endl;
    in.getline(buf,MAXBUFFER);

    //file name of time-variant file
    getline(in,variable_path,' ');
    cout<<variable_path<<endl;
    in.getline(buf,MAXBUFFER);

    //time step (senconds)
    in.getline(buf,MAXBUFFER,' ');
    std::istringstream iss;
    iss.str(buf);
    iss>>time_step;
    cout<<time_step<<endl;
    in.getline(buf,MAXBUFFER);

    //model step for simulation [-]
    in.getline(buf,MAXBUFFER,' ');
    iss.clear();
    iss.str(buf);
    iss>>total_steps;
    cout<<total_steps<<endl;
    in.getline(buf,MAXBUFFER);

    //output variable file span
    in.getline(buf,MAXBUFFER,' ');
    iss.clear();
    iss.str(buf);
    iss>>output_span;
    cout<<output_span<<endl;
    in.getline(buf,MAXBUFFER);

    //output restart file span (second)
    in.getline(buf,MAXBUFFER,' ');
    iss.clear();
    iss.str(buf);
    iss>>restart_span;
    cout<<restart_span<<endl;
    in.getline(buf,MAXBUFFER);

    //model grid row number
    in.getline(buf,MAXBUFFER,' ');
    iss.clear();
    iss.str(buf);
    iss>>row_num;
    cout<<row_num<<endl;
    in.getline(buf,MAXBUFFER);

    //model grid column number
    in.getline(buf,MAXBUFFER,' ');
    iss.clear();
    iss.str(buf);
    iss>>col_num;
    cout<<col_num<<endl;
    in.getline(buf,MAXBUFFER);

    //subregion grid row number
    in.getline(buf,MAXBUFFER,' ');
    iss.clear();
    iss.str(buf);
    iss>>sub_row_num;
    cout<<sub_row_num<<endl;
    in.getline(buf,MAXBUFFER);

    //subregion grid column number
    in.getline(buf,MAXBUFFER,' ');
    iss.clear();
    iss.str(buf);
    iss>>sub_col_num;
    cout<<sub_col_num<<endl;
    in.getline(buf,MAXBUFFER);

    //start row of subregion grid
    in.getline(buf,MAXBUFFER,' ');
    iss.clear();
    iss.str(buf);
    iss>>sub_start_row_num;
    cout<<sub_start_row_num<<endl;
    in.getline(buf,MAXBUFFER);

    //start col of subregion grid
    in.getline(buf,MAXBUFFER,' ');
    iss.clear();
    iss.str(buf);
    iss>>sub_start_col_num;
    cout<<sub_start_col_num<<endl;
    in.getline(buf,MAXBUFFER);

    //patch number of model grid
    in.getline(buf,MAXBUFFER,' ');
    iss.clear();
    iss.str(buf);
    iss>>total_patches;
    cout<<total_patches<<endl;
    in.getline(buf,MAXBUFFER);

    //ensemble size
    in.getline(buf,MAXBUFFER,' ');
    iss.clear();
    iss.str(buf);
    iss>>ensemble_size;
    cout<<ensemble_size<<endl;
    in.getline(buf,MAXBUFFER);

    in.getline(buf,MAXBUFFER);
    in.getline(buf,MAXBUFFER);
    in.getline(buf,MAXBUFFER);

    //clay percent
    getline(in,clay_path,' ');
    cout<<clay_path<<endl;
    in.getline(buf,MAXBUFFER);
    //sand percent
    getline(in,sand_path,' ');
    cout<<sand_path<<endl;
    //mask path
    in.getline(buf,MAXBUFFER);
    in.getline(buf,MAXBUFFER);
    in.getline(buf,MAXBUFFER);

    in.getline(buf,MAXBUFFER);
    getline(in,mask_path,' ');
    cout<<mask_path<<endl;
    in.getline(buf,MAXBUFFER);
    iss.str(buf);
    iss>>lat_span>>lon_span;
    in.getline(buf,MAXBUFFER);
    iss.str(buf);
    iss>>lat_corner>>lon_corner;
    in.close();
}

CoLMConstant::CoLMConstant():ModelParameter()
{
}
CoLMConstant::CoLMConstant(const int patches):ModelParameter(patches)
{
    fcon=new double*[numpatch];
    for(int i=0; i<numpatch; i++)
        fcon[i]=new double [NFCON];
}
CoLMConstant::~CoLMConstant()
{
    if (numpatch>0)
	{
        for(int i=0; i<numpatch; i++)
            delete []fcon[i];
		delete []fcon;
	}
}
void CoLMConstant::patches(const unsigned int p)
{
    ModelParameter::patches(p);
    fcon=new double*[numpatch];
    for(int i=0; i<numpatch; i++)
        fcon[i]=new double [NFCON];
}

void CoLMConstant::open(const char* fn)
{
    ifstream in;
    in.open(fn,ios_base::binary);

    if(!in)
    {
        string error_msg="can't open time invariant constant data file: ";
        error_msg+=fn;
        throw Exception("CoLMConstant","open(fn)",error_msg.c_str());
    }

    //*************************************************************************
    //  1    read patch number
    //*************************************************************************
    in.read((char *) &numpatch,sizeof(numpatch));
    if (is_bigendian()) SWAP_32(numpatch);

    /*  2     read longitude index, latitude index and
       this should be skiped in next version.
    */
    int lon_index[numpatch],lat_index[numpatch];
    double fpatch_patch[numpatch];
    in.read((char *) lon_index,numpatch*sizeof(int));
    in.read((char *) lat_index,numpatch*sizeof(int));
    in.read((char *) fpatch_patch,numpatch*sizeof(double));
    if (is_bigendian())  for(int i=0; i<numpatch; i++)
        {
            SWAP_32(lon_index[i]);
            SWAP_32(lat_index[i]);
            SWAP_FLOAT(fpatch_patch[i]);
        }

    //*************************************************************************
    //  3     read tunalbe variable
    //*************************************************************************
    in.read((char *) ftune,NFTUNE*sizeof(double));   //ftune[NFTUNE]
    if (is_bigendian()) for(int i=0; i<NFTUNE; i++) SWAP_FLOAT(ftune[i]);
    cout<<ftune[0]<<" "<<ftune[1]<<endl;

    //*************************************************************************
    //  4     read constant variable
    //*************************************************************************
    for(int i = 0; i<numpatch; i++)
    {
        in.read((char *) fcon[i],sizeof(double)*NFCON);
        if (is_bigendian()) for (int j=0; j<NFCON; j++) SWAP_FLOAT(fcon[i][j]);
    }
    in.close();
}
void CoLMConstant::save(const char* fn)
{
    ofstream out;
    out.open(fn,ios_base::binary);
    if(!out)
    {
        string error_msg="can't save time invariant constant data file: ";
        error_msg+=fn;
        throw Exception("CoLMConstant","save(fn)",error_msg.c_str());
    }

    //*************************************************************************
    //  1     write patch number
    //*************************************************************************
    out.write((char *) &numpatch,sizeof(numpatch));

    /*  2     write longitude index, latitude index and
    TimeConst_file.write((char *) lon_index,numpatch*sizeof(int));
    TimeConst_file.write((char *) lat_index,numpatch*sizeof(int));
    TimeConst_file.write((char *) fpatch_patch,numpatch*sizeof(float));
    */

    //*************************************************************************
    //  3     write tunalbe variable: ftune[NFTUNE]
    //*************************************************************************
    out.write((char *) ftune,NFTUNE*sizeof(double));


    //*************************************************************************
    //  4     write constant variable
    //*************************************************************************
    for(int i = 0; i<numpatch; i++)
    {
        out.write((char *) &fcon[i],NFCON*sizeof(double));
    }
    out.close();
}

CoLMParameter::CoLMParameter():ModelParameter(0),m_use_lai(false)
{
}
CoLMParameter::CoLMParameter(const int patches):ModelParameter(patches),m_use_lai(false)
{
    fvar=new double* [numpatch];
    for(int i=0; i<numpatch; i++)
        fvar[i]=new double[NFVAR];
}
CoLMParameter::~CoLMParameter()
{
    if (numpatch>0)
	{
		for(int i=0; i<numpatch; i++)
            delete []fvar[i];
		delete []fvar;
	}
    if (m_use_lai)
        delete []lai;
}
void CoLMParameter::patches(const unsigned int p)
{
    ModelParameter::patches(p);
    fvar=new double*[numpatch];
    for(int i=0; i<numpatch; i++)
        fvar[i]=new double [NFVAR];
}
void CoLMParameter::open(const char* fn)
{
    int year, jday, second;
    ifstream in;
    in.open(fn,ios_base::binary);
    if(!in)
    {
        string error_msg="can't open time invariant constant data file: ";
        error_msg+=fn;
        throw Exception("CoLMParameter","open(fn)",error_msg.c_str());
    }

    //read time variant variable
    in.read((char *) &year,sizeof(int));
    in.read((char *) &jday,sizeof(int));
    in.read((char *) &second,sizeof(int));
    if (is_bigendian())
    {
        SWAP_32(year);
        SWAP_32(jday);
        SWAP_32(second);
    }
    //time=Time();
    time.julian(year,jday,second);
    /*
    idate[0] = year;
    idate[1] = jday;
    idate[2] = second;*/

    for(int i = 0; i<numpatch; i++)
    {
        in.read((char *) fvar[i],NFVAR*sizeof(double));
        for(int j=0; j<NFVAR; j++)
        {
            if (is_bigendian()) SWAP_FLOAT(fvar[i][j]);
        }
    }
    in.close();
}
void CoLMParameter::open(const char* fn, const DATA_TYPE dt, const int patches)
{
    ifstream in;
    switch(dt)
    {
    case MODIS_LAI:
        in.open(fn,ios_base::binary);
        if(!in)
        {
            string error_msg="can't open data file: ";
            error_msg+=fn;
            throw Exception("CoLMParameter","open(fn,DATA_TYPE)",error_msg.c_str());
        }
        int p=(patches>0?patches:numpatch);
        lai=new float [p];
        m_use_lai=true;
        in.read((char *) lai,  p*sizeof(float));
        if (is_bigendian())
            for(int i=0; i<p; i++) SWAP_FLOAT(lai[i]);
        in.close();
        break;
    }
}
void CoLMParameter::save(const char* fn)
{
    ofstream out;
    out.open(fn,ios_base::binary);
    if(!out)
    {
        string error_msg="can't save time variant parameter data file: ";
        error_msg+=fn;
        throw Exception("CoLMParameter","save(fn)",error_msg.c_str());
    }

    int year,jday,second;
    year=time.getYear();
    jday=time.getDays();
    second=time.getSeconds();
    /*
    year = idate[0];
    jday = idate[1];
    second = idate[2];*/

    /* is this required?
    if (is_bigendian())
    {
        SWAP_32(year);
        SWAP_32(jday);
        SWAP_32(second);
    }*/

    //write date
    out.write((char *) &year,sizeof(year));
    out.write((char *) &jday,sizeof(jday));
    out.write((char *) &second,sizeof(second));

    //write restart variable
    for(int i=0; i<numpatch; i++)
    {
        // is this required? please test this in the big endian machine.
        //for(int j=0; j<NFVAR; j++)   if (is_bigendian()) SWAP_FLOAT(fvar_patch[i][j]);
        out.write((char *) fvar[i],sizeof(double)*NFVAR);
    }
    out.close();
}

CoLMForcing::CoLMForcing():numpatch(0)
{
}
CoLMForcing::CoLMForcing(const int sum):numpatch(sum)
{
    forc=new float* [NFORC_IN];
    for(int i=0; i<NFORC_IN; i++)
        forc[i]=new float[numpatch];
}
CoLMForcing::~CoLMForcing()
{
    if (numpatch>0)
    {
        for(int i=0; i<NFORC_IN; i++)
            delete []forc[i];
        delete []forc;
    }
}
void CoLMForcing::patches(const unsigned int p)
{
    numpatch=p;
    forc=new float* [NFORC_IN];
    for(int i=0; i<NFORC_IN; i++)
        forc[i]=new float[numpatch];
}
void CoLMForcing::open(const char* fn)
{
    // 0  u10   10 m u wind (m/sec)
    // 1   v10   10 m v wind (m/sec)
    // 2   swd   downward short wave rad (W/m2)
    // 3   lwd   downward long wave rad (W/m2)
    // 4   rc   accum conv pcn (mm/s)
    // 5   rn   accum non-c pcn (mm/s)
    // 6   t2m  2 m temperature (K)
    // 7   q2m  2 m mixing ratio (kg/kg)
    // 8   psfc Surface pressure(hPa)
    ifstream in;
    in.open(fn,ios_base::binary);
    if(!in)
    {
        string error_msg="can't open forcing data file: ";
        error_msg+=fn;
        throw Exception("CoLMForcing","open(fn)",error_msg.c_str());
    }
    for(int i=0; i<NFORC_IN; i++)
    {
        in.read((char *) forc[i],sizeof(float)*numpatch);
        if (is_bigendian())for (int j=0; j<numpatch; j++)SWAP_FLOAT(forc[i][j]);
    }
    in.close();
}

CommonLandModel::CommonLandModel()
{
}

CommonLandModel::~CommonLandModel()
{
    if (m_config.total_patches>0)
    {
        for(int i=0; i<m_config.total_patches; i++)
        {
            delete []flux_out[i];
            delete []fvar_out[i];
        }
        delete []flux_out;
        delete []fvar_out;
    }
}

CommonLandModel::CommonLandModel(const CommonLandModel& other)
{
    //copy ctor
    //todo: deep copy, copy all local variable
}

CommonLandModel& CommonLandModel::operator=(const CommonLandModel& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    //todo: copy variable
    return *this;
}
CommonLandModel::CommonLandModel(const std::string config_filename)
{
    m_config.open(config_filename.c_str());
}

void CommonLandModel::init()
{
    //m_config.open("input.txt");
    grid.rows(m_config.row_num);
    grid.cols(m_config.col_num);
    //this should be configured in the model config class
    grid.lat_span(m_config.lat_span);
    grid.lon_span(m_config.lon_span);
    grid.corner(LandPoint(m_config.lat_corner,m_config.lon_corner));
    grid.mask(m_config.mask_path.c_str());
    parameter.patches(grid.patches());
    parameter.open(m_config.variable_path.c_str());
    constant.patches(grid.patches());
    constant.open(m_config.constant_path.c_str());
    //the forcing data is a full grid now.
    forcing.patches(m_config.row_num*m_config.col_num);
    //dynamic array define
    flux_out=new double*[m_config.total_patches];
    fvar_out=new double*[m_config.total_patches];
    for(int i=0; i<m_config.total_patches; i++)
    {
        flux_out[i]=new double[OUTFLUX];
        fvar_out[i]=new double[VARNUM];
    }
}
void CommonLandModel::run()
{
    int month=0;//to judge with lai month data
    //*************************************************************************
    // begin time stepping loop
    // calculate land state one by one
    //*************************************************************************
    for (int istep=1; istep<=m_config.total_steps; istep++)
    {
        Time t_time=parameter.time;
        cout<<"***************************************************************"<<endl;
        cout<<"Current Date (Year, Day, Hour):  "<<parameter.time.getYear()<<"  "<<parameter.time.getDays()<<"   "<<parameter.time.getHour()<<endl;

        //*************************************************************************
        //  step 1     read  meteorological forcing data and LAI data
        //*************************************************************************
        cout<<"  Read meteorolgical forcing data and LAI data..."<<endl;
        int ym = parameter.time.getYear()*100+parameter.time.getMonth();
        int ydh=parameter.time.getYear()*100000+parameter.time.getDays()*100+parameter.time.getHour();
        std::stringstream sst;
        sst<<m_config.forcing_path<<ym<<"/"<<ydh<<".forc";
        forcing.open(sst.str().c_str());
        if (month!=parameter.time.getMonth())
        {
            sst.str("");
            sst<<m_config.vegetation_path<<"LAI_"<<ym<<".900s";
            parameter.open(sst.str().c_str(),MODIS_LAI,m_config.row_num*m_config.col_num);
            month=parameter.time.getMonth();
        }
        // Calendar for NEXT time step
        parameter.time.addSeconds(m_config.time_step);

        //********************************************************
        //step 3 simulated patch by patch
        //********************************************************
#pragma omp parallel for
        for (int ipatch=0; ipatch<m_config.total_patches; ipatch++) run(ipatch);

        //*************************************************************************
        // step 4    output simulated variable
        //*************************************************************************
        cout<<" OutPut results..."<<endl;
        if ((istep%m_config.output_span==0)||(istep==m_config.total_steps))
        {
            // output temperature
            for (int ipatch=0; ipatch<m_config.total_patches; ipatch++)
            {

                for (int i=0; i<15; i++)
                    fvar_out[ipatch][i] = parameter.fvar[ipatch][30+i];
                // output water content
                for (int i=0; i<15; i++)
                    fvar_out[ipatch][i+15] = parameter.fvar[ipatch][45+i];
                // output ice content
                for (int i=0; i<15; i++)
                    fvar_out[ipatch][i+30] = parameter.fvar[ipatch][60+i];

                // ground surface temperature
                fvar_out[ipatch][45] = parameter.fvar[ipatch][75];
                // snow cover
                fvar_out[ipatch][46] = parameter.fvar[ipatch][80];
                // snow depth
                fvar_out[ipatch][47] = parameter.fvar[ipatch][81];
            }
            sst.str("");
            sst<<m_config.output_path<<ym<<"/"<<ydh<<".grid";
            output(sst.str().c_str(),VARNUM,fvar_out);
            // output flux
            sst.str("");
            sst<<m_config.output_path<<ym<<"/"<<ydh<<".flux";
            output(sst.str().c_str(),FLUXNUM,flux_out);
        }

        //*************************************************************************
        // step 5    write restart file
        //*************************************************************************
        if ((istep%m_config.restart_span==0)||(istep==m_config.total_steps))
        {
            std::stringstream ss;
            ss<<m_config.restart_path<<ydh<<".patch";
            parameter.save(ss.str().c_str());
        }
        //*************************************************************************
        cout<<"caculation is finished in current time"<<endl;
        //*************************************************************************
    }
}
void CommonLandModel::run(const int ipatch)
{
    //*************************************************************************
    // 3.1   allocate atmospheric data to patch
    //*************************************************************************
    //int location = row*ColNum + col;
    int location=grid.location(ipatch);
    double swdown = forcing.forc[2][location];
    double lat = constant.fcon[ipatch][0]; //弧度单位，非度数单位
    double lon = constant.fcon[ipatch][1];
    double sols   = 0.0;
    double soll   = 0.0;
    double solsd  = 0.0;
    double solld  = 0.0;
    //todo: is this method required in Focing class?
    radiation(parameter.time.getDays(),parameter.time.getSeconds(),lon,lat, swdown,sols,soll,solsd, solld);

    double forc_data[NFORC];
    forc_data[0] = forcing.forc[8][location]*355.e-06;
    forc_data[1] = forcing.forc[8][location]*0.209;
    forc_data[2] = forcing.forc[0][location];
    forc_data[3] = forcing.forc[1][location];
    forc_data[4] = forcing.forc[6][location];
    forc_data[5] = forcing.forc[7][location];
    forc_data[6] = forcing.forc[4][location]/10.0;
    forc_data[7] = forcing.forc[5][location]/10.0;
    forc_data[8] = forcing.forc[8][location];
    forc_data[9] = forcing.forc[8][location];
    forc_data[10] = sols;
    forc_data[11] = soll;
    forc_data[12] = solsd;
    forc_data[13] = solld;
    forc_data[14] = forcing.forc[3][location];
    forc_data[15] = 2;
    forc_data[16] = 2;
    forc_data[17] = 2;

    //10  Savanna
    //11  Deciduous Broadleaf Forest
    //12  Deciduous Needleleaf Forest
    //13  Evergreen Broadleaf Forest
    //14  Evergreen Needleleaf Forest
    //15  Mixed Forest
    if ((constant.fcon[ipatch][3]>=10)&&(constant.fcon[ipatch][3]<=15))
    {
        forc_data[15] = 15;
        forc_data[16] = 15;
        forc_data[17] = 15;
    }
    //18  Wooded Wetland
    //19  Barren or Sparsely Vegetated
    //20  Herbaceous Tundra
    //21  Wooded Tundra
    if ((constant.fcon[ipatch][3]>=18)&&(constant.fcon[ipatch][3]<=24))
    {
        forc_data[15] = 15;
        forc_data[16] = 15;
        forc_data[17] = 15;
    }

    //*************************************************************************
    //3.2 add lai from MODIS LAI production 0.25*0.25deg
    //*************************************************************************
    parameter.fvar[ipatch][86] = parameter.lai[location];
    //********************************************************************
    // step 3.3  calculate forecast variable of patch
    //********************************************************************
    run(constant.ftune,constant.fcon[ipatch],forc_data,
        parameter.fvar[ipatch],flux_out[ipatch]);
}
int CommonLandModel::snowlayer() const
{
	return snl;
}
void CommonLandModel::run (double ftune[NFTUNE], double fcon[NFCON], double forc[NFORC], double fvar[NFVAR],double flux[OUTFLUX])
{
    // ----------------------------------------------------------------
    // I. Time invariant model variables
    // ----------------------------------------------------------------
    double   dlat;      // latitude in radians
    double   dlon;      // longitude in radians
    int      itypwat ;  // land water type
    int      ivt;       // land cover type of classification of USGS etc.
    snl=0; //snow layer number
    double oro=1; //ocean(0)/seaice(2)/ flag, it is all land patches here.
    // true if time for time-varying vegetation paramter
    bool dolai=true;
    // true if time for surface albedo calculation
    bool doalb=true;
    // true if time for update sst/ice/snow
    bool dosst=true;

    // Soil physical parameters
    double   albsol;          // soil albedo for different coloured soils [-]
    double   csol  [MAXSOILL]; // heat capacity of soil solids [J/(m3 K)]
    double   porsl [MAXSOILL]; // fraction of soil that is voids [-]
    double   phi0  [MAXSOILL]; // minimum soil suction [mm]
    double   bsw   [MAXSOILL]; // clapp and hornbereger "b" parameter [-]
    double   dkmg  [MAXSOILL]; // thermal conductivity of soil minerals [W/m-K]
    double   dksatu[MAXSOILL]; // thermal conductivity of saturated soil [W/m-K]
    double   dkdry [MAXSOILL]; // thermal conductivity for dry soil  [W/(m-K)]
    double   hksati[MAXSOILL]; // hydraulic conductivity at saturation [mm h2o/s]

    // Vegetation static parameters
    double   z0m      ; // aerodynamic roughness length [m]
    double   displa   ; // displacement height [m]
    double   sqrtdi   ; // inverse sqrt of leaf dimension [m**-0.5]
    double   effcon   ; // quantum efficiency of RuBP regeneration (molCO2/molquanta)
    double   vmax25   ; // maximum carboxylation rate at 25 C at canopy top
    double   slti     ; // s3: slope of low temperature inhibition function
    double   hlti     ; // s4: 1/2 point of low temperature inhibition function
    double   shti     ; // s1: slope of high temperature inhibition function
    double   hhti     ; // s2: 1/2 point of high temperature inhibition function
    double   trda     ; // s5: temperature coefficient in gs-a model
    double   trdm     ; // s6: temperature coefficient in gs-a model
    double   trop     ; // temperature coefficient in gs-a model
    double   gradm    ; // conductance-photosynthesis slope parameter
    double   binter   ; // conductance-photosynthesis intercep
    double   extkn    ; // coefficient of leaf nitrogen allocation
    double   chil     ; // leaf angle distribution factor
    double   ref[2][2]; // leaf reflectance (iw=iband, il=life and dead)
    double   tran[2][2]; // leaf transmittance (iw=iband, il=life and dead)
    double   rootfr[MAXSOILL];  // fraction of roots in each soil layer

    // CLM time step and TUNABLE constants
    double   dtime=m_config.time_step  ;            // CLM time step [seconds]
    double   zlnd   ;            // roughness length for soil [m]
    double   zsno   ;            // roughness length for snow [m]
    double   csoilc ;            // drag coefficient for soil under canopy [-]
    double   dewmx  ;            // maximum dew
    double   wtfact ;            // fraction of model area with high water table
    double   capr   ;            // tuning factor to turn first layer T into surface T
    double   cnfac  ;            // Crank Nicholson factor between 0 and 1
    double   ssi    ;            // irreducible water saturation of snow
    double   wimp   ;            // water impremeable if porosity less than wimp
    double   pondmx ;            // ponding depth (mm)
    double   smpmax ;            // wilting point potential in mm
    double   smpmin ;            // restriction for min of soil poten. (mm)
    double   trsmx0 ;            // max transpiration for moist soil+100% veg.  [mm/s]
    double   tcrit  ;            // critical temp. to determine rain or snow

    // -----------------------------------------------------------------
    // II. Time-varying state variables which reaquired by restart run
    // -----------------------------------------------------------------
    int year=parameter.time.getYear();                 // current year of model run
    int jday=parameter.time.getDays();                 // current julian day of model run
    int msec=parameter.time.getSeconds();              // current seconds of model run (0 - 86400)

    // Main land surface variables
    double   z   [MAXSNL+MAXSOILL]; // node depth [m]
    double   dz  [MAXSNL+MAXSOILL]; // interface depth [m]
    double   tss [MAXSNL+MAXSOILL]; // soil temperature [K]
    double   wliq[MAXSNL+MAXSOILL]; // liquid water in layers [kg/m2]
    double   wice[MAXSNL+MAXSOILL]; // ice lens in layers [kg/m2]
    double   tg       ; // ground surface temperature [K]
    double   tlsun    ; // sunlit leaf temperature [K]
    double   tlsha    ; // shaded leaf temperature [K]
    double   ldew     ; // depth of water on foliage [mm]
    double   sag      ; // non dimensional snow age [-]
    double   scv      ; // snow cover, water equivalent [mm]
    double   snowdp   ; // snow depth [meter]

    // Vegetation dynamic parameters
    double   fveg     ; // fraction of vegetation cover
    double   fsno     ; // fraction of snow cover on ground
    double   sigf     ; // fraction of veg cover, excluding snow-covered veg [-]
    double   green    ; // leaf greenness
    double   lai      ; // leaf area index
    double   sai      ; // stem area index

    // Radiation  related (albedoes)
    double   coszen   ; // cosine of solar zenith angle
    double   albg[2][2]; // albedo, ground [-]
    double   albv[2][2];// albedo, vegetation [-]
    double   alb[2][2];// averaged albedo [-]
    double   ssun[2][2]; // sunlit canopy absorption for solar radiation (0-1)
    double   ssha[2][2]; // shaded canopy absorption for solar radiation (0-1)
    double   thermk   ; // canopy gap fraction for tir radiation
    double   extkb    ; // (k, g(mu)/mu) direct solar extinction coefficient
    double   extkd    ; // diffuse and scattered diffuse PAR extinction coefficient

    // Additional variables required by reginal model (WRF & RSM)
    double   trad     ; // radiative temperature of surface [K]
    double   tref     ; // 2 m height air temperature [kelvin]
    double   qref     ; // 2 m height air specific humidity
    double   rst      ; // canopy stomatal resistance (s/m)

    double   emis     ; // averaged bulk surface emissivity
    double   z0ma     ; // effective roughness [m]
    double   zol      ; // dimensionless height (z/L) used in Monin-Obukhov theory
    double   rib      ; // bulk Richardson number in surface layer
    double   ustar    ; // u* in similarity theory [m/s]
    double   qstar    ; // q* in similarity theory [kg/kg]
    double   tstar    ; // t* in similarity theory [K]
    double   u10m     ; // 10m u-velocity
    double   v10m     ; // 10m v-velocity
    double   f10m     ; // integral of profile function for momentum at 10m
    double   fm       ; // integral of profile function for momentum
    double   fh       ; // integral of profile function for heat
    double   fq       ; // integral of profile function for moisture

    // -----------------------------------------------------------------
    // III. Forcing
    // -----------------------------------------------------------------
    // Forcing
    double   pco2m    ; // CO2 concentration in atmos. (35 pa)
    double   po2m     ; // O2 concentration in atmos. (20900 pa)
    double   us       ; // wind in eastward direction [m/s]
    double   vs       ; // wind in northward direction [m/s]
    double   tm       ; // temperature at reference height [kelvin]
    double   qm       ; // specific humidity at reference height [kg/kg]
    double   prc      ; // convective precipitation [mm/s]
    double   prl      ; // large scale precipitation [mm/s]
    double   psrf     ; // atmospheric pressure at the surface [pa]
    double   pbot     ; // atm bottom level pressure (or reference height) (pa)
    double   sols     ; // atm vis direct beam solar rad onto srf [W/m2]
    double   soll     ; // atm nir direct beam solar rad onto srf [W/m2]
    double   solsd    ; // atm vis diffuse solar rad onto srf [W/m2]
    double   solld    ; // atm nir diffuse solar rad onto srf [W/m2]
    double   frl      ; // atmospheric infrared (longwave) radiation [W/m2]
    double   hu       ; // observational height of wind [m]
    double   ht       ; // observational height of temperature [m]
    double   hq       ; // observational height of humidity [m]
    double   rhoair   ; // air density [kg/m3]

    // -----------------------------------------------------------------
    // IV. Fluxes
    // -----------------------------------------------------------------
    double   taux     ; // wind stress: E-W [kg/m/s2]
    double   tauy     ; // wind stress: N-S [kg/m/s2]
    double   fsena    ; // sensible heat from canopy height to atmosphere [W/m2]
    double   lfevpa   ; // latent heat flux from canopy height to atmosphere [W/m2]
    double   fevpa    ; // evapotranspiration from canopy to atmosphere [mm/s]
    double   fsenl    ; // sensible heat from leaves [W/m2]
    double   fevpl    ; // evaporation+transpiration from leaves [mm/s]
    double   etr      ; // transpiration rate [mm/s]
    double   fseng    ; // sensible heat flux from ground [W/m2]
    double   fevpg    ; // evaporation heat flux from ground [mm/s]
    double   fgrnd    ; // ground heat flux [W/m2]
    double   sabvsun  ; // solar absorbed by sunlit vegetation [W/m2]
    double   sabvsha  ; // solar absorbed by shaded vegetation [W/m2]
    double   sabg     ; // solar absorbed by ground  [W/m2]
    double   olrg     ; // outgoing long-wave radiation from ground+canopy [W/m2]
//	double   rnet     ; // net radiation by surface [W/m2]
    double   xerr=0     ; // the error of water banace [mm/s]
    double   zerr     ; // the error of energy balance [W/m2]

    double   rsur     ; // surface runoff (mm h2o/s)
    double   rnof=0     ; // total runoff (mm h2o/s)
    double   assim    ; // canopy assimilation rate (mol m-2 s-1)
    double   respc    ; // canopy respiration (mol m-2 s-1)


    // -----------------------------------------------------------------
    // V. Local declaration
    // -----------------------------------------------------------------
    int j;         // loop/array indices
    int nl_soil,nl_snow;
    nl_soil = MAXSOILL;
    nl_snow = MAXSNL;

    double   parsun   ; // PAR by sunlit leaves [W/m2]
    double   parsha   ; // PAR by shaded leaves [W/m2]
    double   sabvg    ; // solar absorbed by ground + vegetation [W/m2]

    // ======================================================================
    //  [1] Transfer the time invariant and time-varying variables
    // ======================================================================

    // Time invariant model variables
    dlat   = fcon[0];                 //0
    dlon   = fcon[1];                //1
    itypwat= (int)(fcon[2]);         //2
    ivt    = (int)(fcon[3]);         //3
    albsol = fcon[4];                  //4

    for (j=0; j<MAXSOILL; j++)               //(5-14)
        csol[j] = fcon[5+j];

    for (j=0; j<MAXSOILL; j++)               //(15-24)
        porsl[j] = fcon[15+j];

    for (j=0; j<MAXSOILL; j++)               //(25-34)
        phi0[j] = fcon[25+j] ;

    for (j=0; j<MAXSOILL; j++)               //(35-44)
        bsw[j] = fcon[35+j] ;

    for (j=0; j<MAXSOILL; j++)               //(45-54)
        dkmg[j] = fcon[45+j] ;

    for (j=0; j<MAXSOILL; j++)               //(55-64)
        dksatu[j] = fcon[55+j] ;

    for (j=0; j<MAXSOILL; j++)               //(65-74)
        dkdry[j] = fcon[65+j] ;

    for (j=0; j<MAXSOILL; j++)               //(75-84)
        hksati[j] = fcon[75+j] ;

    //ub = ub + 1;
    z0m       = fcon[85];             //85
    displa    = fcon[86];             //86
    sqrtdi    = fcon[87];             //87
    effcon    = fcon[88];             //88
    vmax25    = fcon[89];             //89
    slti      = fcon[90];             //90
    hlti      = fcon[91];             //91
    shti      = fcon[92];             //92
    hhti      = fcon[93];             //93
    trda      = fcon[94];             //94
    trdm      = fcon[95];             //95
    trop      = fcon[96];             //96
    gradm     = fcon[97];             //97
    binter    = fcon[98];             //98
    extkn     = fcon[99];             //99
    chil      = fcon[100];             //100
    ref [0][0] = fcon[101];             //101
    ref [0][1] = fcon[102];             //102
    ref [1][0] = fcon[103];             //103
    ref [1][1] = fcon[104];             //104
    tran[0][0] = fcon[105];             //105
    tran[0][1] = fcon[106];             //106
    tran[1][0] = fcon[107];             //107
    tran[1][1] = fcon[108];          //108

    for (j=0; j<MAXSOILL; j++)               //(109-118)
        rootfr[j] = fcon[109+j];

    // CLM time step and TUNABLE constants
    // dtime  = deltim;
    zlnd   = ftune[0];
    zsno   = ftune[1];
    csoilc = ftune[2];
    dewmx  = ftune[3];
    wtfact = ftune[4];
    capr   = ftune[5];
    cnfac  = ftune[6];
    ssi    = ftune[7];
    wimp   = ftune[8];
    pondmx = ftune[9];
    smpmax = ftune[10];
    smpmin = ftune[11];
    trsmx0 = ftune[12];
    tcrit  = ftune[13];

    // Time-varying variables
    for (j=0; j<15; j++)                       //(0-14)
        z[j] = fvar[0+j];

    for (j=0; j<15; j++)                       //(15-34)
        dz[j] = fvar[15+j];

    for (j=0; j<15; j++)                       //(30-44)
        tss[j] = fvar[30+j] ;

    for (j=0; j<15; j++)                       //(45-59)
        wliq[j] = fvar[45+j] ;


    for (j=0; j<15; j++)                       //(60-74)
        wice[j] = fvar[60+j] ;


    tg      = fvar[75];                    //75
    tlsun   = fvar[76];                    //76
    tlsha   = fvar[77];                    //77
    ldew    = fvar[78];                    //78
    sag     = fvar[79];                    //79
    scv     = fvar[80];                    //80
    snowdp  = fvar[81];                    //81
    fveg    = fvar[82];                    //82
    fsno    = fvar[83];                    //83
    sigf    = fvar[84];                    //84
    green   = fvar[85];                    //85
    lai     = fvar[86];                    //86
    sai     = fvar[87];                    //87
    coszen  = fvar[88];                    //88
    albg[0][0] = fvar[89];                  //89
    albg[0][1] = fvar[90];                  //90
    albg[1][0] = fvar[91];                  //91
    albg[1][1] = fvar[92];                  //92
    albv[0][0] = fvar[93];                  //93
    albv[0][1] = fvar[94];                  //94
    albv[1][0] = fvar[95];                  //95
    albv[1][1] = fvar[96];                  //96
    alb[0][0] = fvar[97];                  //97
    alb[0][1] = fvar[98];                  //98
    alb[1][0] = fvar[99];                  //99
    alb[1][1] = fvar[100];                  //100
    ssun[0][0] = fvar[101];                  //101
    ssun[0][1] = fvar[102];                  //102
    ssun[1][0] = fvar[103];                  //103
    ssun[1][1] = fvar[104];                  //104
    ssha[0][0] = fvar[105];                  //105
    ssha[0][1] = fvar[106];                  //106
    ssha[1][0] = fvar[107];                  //107
    ssha[1][1] = fvar[108];                  //108
    thermk  = fvar[109];                    //109
    extkb   = fvar[110];                    //110
    extkd   = fvar[111];                    //111

    trad      = fvar[112];                  //112
    tref      = fvar[113];               //113
    qref      = fvar[114] ;               //114
    rst       = fvar[115] ;               //115
    emis      = fvar[116] ;               //116
    z0ma      = fvar[117] ;               //117
    zol       = fvar[118] ;               //118
    rib       = fvar[119]  ;               //119
    ustar     = fvar[120] ;               //120
    qstar     = fvar[121]  ;               //121
    tstar     = fvar[122]  ;               //122
    fm        = fvar[123]  ;               //123
    fh        = fvar[124]  ;               //124
    fq        = fvar[125]  ;               //125

    // ======================================================================
    //  [2] atmospheric fields to force colm
    // ======================================================================
    pco2m   = forc[0]   ;   //CO2 concentration in atmos. (35 pa)
    po2m    = forc[1]   ;   //O2 concentration in atmos. (20900 pa)
    us      = forc[2]   ;   //wind in eastward direction [m/s]
    vs      = forc[3]   ;   //wind in northward direction [m/s]
    tm      = forc[4]   ;   //temperature at reference height [kelvin]
    qm      = forc[5]   ;   //specific humidity at reference height [kg/kg]
    prc     = forc[6]   ;   //convective precipitation [mm/s]
    prl     = forc[7]   ;   //large scale precipitation [mm/s]
    pbot    = forc[8]   ;   //atm bottom level pressure (or reference height) (pa)
    psrf    = forc[9]   ;   //atmospheric pressure at the surface [pa]
    sols    = forc[10]  ;   //atm vis direct beam solar rad onto srf [W/m2]
    soll    = forc[11]  ;   //atm nir direct beam solar rad onto srf [W/m2]
    solsd   = forc[12]  ;   //atm vis diffuse solar rad onto srf [W/m2]
    solld   = forc[13]  ;   //atm nir diffuse solar rad onto srf [W/m2]
    frl     = forc[14]  ;   //atmospheric infrared (longwave) radiation [W/m2]
    hu      = forc[15]  ;   //observational height of wind [m]
    ht      = forc[16]  ;   //observational height of temperature [m]
    hq      = forc[17]  ;   //observational height of humidity [m]
    //cout<<"rhoair  "<<qm<<" "<<pbot<<" "<<tm<<endl;
    rhoair= (pbot-0.378*qm*pbot/(0.622+0.378*qm))/(rgas*tm);

    // ======================================================================
    // [2] Main driver for CLM
    // ======================================================================
    CLMMAIN (dtime, doalb, dolai, dosst, nl_soil,
             nl_snow, dlon , dlat , itypwat , ivt , oro ,

             // soil information
             albsol , csol , porsl , phi0 ,
             bsw , dkmg , dksatu , dkdry , hksati ,

             // vegetation information
             z0m , displa , sqrtdi ,
             effcon, vmax25 , slti , hlti , shti , hhti ,
             trda , trdm , trop , gradm , binter , extkn ,
             chil , ref ,tran , rootfr ,

             // atmospheric forcing
             frl , sols , soll , solsd , solld ,
             pco2m , po2m , us , vs , tm , qm ,
             prc , prl , psrf , rhoair ,
             hu , ht , hq ,

             // model variables needed by restart run
             year, jday, msec,
             z , dz , tss ,
             wliq , wice ,
             tg , tlsun , tlsha , ldew ,
             sag , scv , snowdp ,
             fveg , fsno , sigf , green , lai , sai ,
             coszen , albg , albv , alb ,
             ssun , ssha , thermk , extkb , extkd ,

             // fluxes
             taux , tauy ,
             fsena , fevpa , lfevpa , fsenl , fevpl , etr ,
             fseng , fevpg , olrg , fgrnd , trad , tref , qref ,
             rsur , rnof , rst , assim , respc ,
             parsun ,parsha ,sabvsun ,sabvsha ,sabg ,sabvg ,
             xerr , zerr ,

             // TUNABLE modle constants
             zlnd, zsno, csoilc, dewmx, wtfact,
             capr, cnfac, ssi, wimp, pondmx,
             smpmax, smpmin, trsmx0, tcrit,

             // additional variables required by coupling with regional model (WRF & RSM)
             emis , z0ma , zol , rib , ustar , qstar , tstar ,
             u10m , v10m , f10m , fm , fh , fq,
             snl
            );

    // ======================================================================
    // [3] the model variables for restart run
    // ======================================================================
    for (j=0; j<MAXSNL+MAXSOILL; j++)                       //(0-14)
        fvar[0+j]  = z[j] ;

    for (j=0; j<MAXSNL+MAXSOILL; j++)                       //(15-29)
        fvar[15+j] = dz[j] ;

    for (j=0; j<MAXSNL+MAXSOILL; j++)                       //(30-44)
        fvar[30+j] = tss[j] ;

    for (j=0; j<MAXSNL+MAXSOILL; j++)                       //(45-59)
        fvar[45+j] = wliq[j] ;

    for (j=0; j<MAXSNL+MAXSOILL; j++)                       //(60-74)
        fvar[60+j] = wice[j] ;

    fvar[75] 	=	tg 	;              //75
    fvar[76] 	=	tlsun 	;          //76
    fvar[77] 	=	tlsha 	;          //77
    fvar[78] 	=	ldew 	;          //78
    fvar[79] 	=	sag 	;          //79
    fvar[80] 	=	scv 	;          //80
    fvar[81] 	=	snowdp 	;          //81
    fvar[82] 	=	fveg 	;              //82
    fvar[83] 	=	fsno 	;              //83
    fvar[84] 	=	sigf 	;              //84
    fvar[85] 	=	green 	;          //85
    fvar[86] 	=	lai 	;              //86
    fvar[87] 	=	sai 	;              //87
    fvar[88] 	=	coszen 	;          //88
    fvar[89] 	=	albg[0][0] 	;      //89
    fvar[90] 	=	albg[0][1] 	;      //90
    fvar[91] 	=	albg[1][0] 	;      //91
    fvar[92] 	=	albg[1][1] 	;      //92
    fvar[93] 	=	albv[0][0] 	;      //93
    fvar[94] 	=	albv[0][1] 	;      //94
    fvar[95] 	=	albv[1][0] 	;      //95
    fvar[96] 	=	albv[1][1] 	;      //96
    fvar[97] 	=	alb[0][0] 	;      //97
    fvar[98] 	=	alb[0][1] 	;      //98
    fvar[99] 	=	alb[1][0] 	;      //99
    fvar[100] 	=	alb[1][1] 	;      //100
    fvar[101] 	=	ssun[0][0] 	;      //101
    fvar[102] 	=	ssun[0][1] 	;      //102
    fvar[103] 	=	ssun[1][0] 	;      //103
    fvar[104] 	=	ssun[1][1] 	;      //104
    fvar[105] 	=	ssha[0][0] 	;      //105
    fvar[106] 	=	ssha[0][1] 	;      //106
    fvar[107] 	=	ssha[1][0] 	;      //107
    fvar[108] 	=	ssha[1][1] 	;      //108
    fvar[109] 	=	thermk 	;          //109
    fvar[110] 	=	extkb 	;          //110
    fvar[111] 	=	extkd 	;          //111

    // Additional variables required by reginal model (WRF & RSM)
    fvar[112]    = trad     ;               //112
    fvar[113]    = tref     ;               //113
    fvar[114]    = qref     ;               //114
    fvar[115]    = rst      ;               //115
    fvar[116]    = emis     ;               //116
    fvar[117]    = z0ma     ;               //117
    fvar[118]    = zol      ;               //118
    fvar[119]    = rib      ;               //119
    fvar[120]    = ustar    ;               //120
    fvar[121]    = qstar    ;               //121
    fvar[122]    = tstar    ;               //122
    fvar[123]    = fm       ;               //123
    fvar[124]    = fh       ;               //124
    fvar[125]    = fq       ;               //125

    // ======================================================================
    // [4] output variable and flux for assimilation
    //  add by Huang C L   Nov 9, 2005
    // ======================================================================
    flux[0]  = taux      ;   //wind stress: E-W [kg/m/s2]
    flux[1]  = tauy      ;   //wind stress: N-S [kg/m/s2]
    flux[2]  = fsena     ;   //sensible heat from canopy height to atmosphere [W/m2]
    flux[3]  = lfevpa    ;   //latent heat flux from canopy height to atmosphere [W/m2]
    flux[4]  = fevpa     ;   //evapotranspiration from canopy to atmosphere [mm/s]
    flux[5]  = fsenl     ;   //sensible heat from leaves [W/m2]
    flux[6]  = fevpl     ;   //evaporation+transpiration from leaves [mm/s]
    flux[7]  = etr       ;   //transpiration rate [mm/s]
    flux[8]  = fseng     ;   //sensible heat flux from ground [W/m2]
    flux[9]  = fevpg     ;   //evaporation heat flux from ground [mm/s]
    flux[10] = fgrnd     ;   //ground heat flux [W/m2]
    flux[11] = sabvsun   ;   //solar absorbed by sunlit canopy [W/m2]
    flux[12] = sabvsha   ;   //solar absorbed by shaded [W/m2]
    flux[13] = sabg      ;   //solar absorbed by ground  [W/m2]
    flux[14] = olrg      ;   //outgoing long-wave radiation from ground+canopy [W/m2]
    flux[15] = sabvg    + frl  - olrg  ;   //net radiation [W/m2]
    flux[16] = xerr      ;   //the error of water banace [mm/s]
    flux[17] = zerr      ;   //the error of energy balance [W/m2]
    flux[18] = rsur      ;   //surface runoff [mm/s]
    flux[19] = rnof      ;   //total runoff [mm/s]
    flux[20] = assim     ;   //canopy assimilation rate [mol m-2 s-1]
    flux[21] = respc     ;   //respiration (plant+soil) [mol m-2 s-1]
}
void CommonLandModel::output(const char* fn,const int num, double** outgrid)
{
    double *tmp_grid  = new double [m_config.row_num*m_config.col_num];
    for(int i=0; i<m_config.row_num*m_config.col_num; i++)
        tmp_grid[i]=NULL_VALUE;
    ofstream out;
    out.open(fn,ios_base::binary);
    if(!out)
    {
        string error_msg="can't save output data file: ";
        error_msg+=fn;
        throw Exception("CommonLandModel","output()",error_msg.c_str());
    }

    for(int i=0; i<num; i++)
    {
        for(int j=0; j<m_config.total_patches; j++)
            tmp_grid[grid.location(j)] = outgrid[j][i];
        out.write((char *) tmp_grid,m_config.row_num*m_config.col_num*sizeof(double));
    }

    out.close();
    delete [] tmp_grid;
}
void CommonLandModel::config(const char* fn)
{
    m_config.open(fn);
    init();
}
void CommonLandModel::albland(int    itypwat,   double albsol,   double chil,      double ref[2][2], double tran[2][2],
                              double fveg,      double green,    double lai,       double sai,       double coszen,
                              double wt,        double fsno,     double scv,       double sag,       double ssw,
                              double tg,        double alb[2][2],double albg[2][2],double albv[2][2],double ssun[2][2],
                              double ssha[2][2],double &thermk,   double &extkb,     double &extkd)
{

    //=======================================================================
    // Calculates fragmented albedos (direct and diffuse) in
    // wavelength regions split at 0.7um.
    //
    // (1) soil albedos: as in BATS formulations, which are the function of
    //     soil color and moisture in the surface soil layer
    // (2) snow albedos: as in BATS formulations, which are inferred from
    //     the calculations of Wiscombe and Warren (1980) and the snow model
    //     and data of Anderson(1976), and the function of snow age, grain size,
    //     solar zenith angle, pollution, the amount of the fresh snow
    // (3) canopy albedo: has been developed to capture the essential features of
    //     a two-stream approximation model while forgoing the complexity of the full
    //     treatment. It combines soil and canopy albedo by simple rules that are
    //     formulated to reduce to correct asymptotic limits for thick and thin canopies
    //     and provide reasonable results for intermediate values of leaf area index
    // (4) glacier albedos: as in BATS, which are set to constants (0.8 for visible beam,
    //     0.55 for near-infrared)
    // (5) lake and wetland albedos: as in BATS, which depend on cosine solar zenith angle,
    //     based on data in Henderson-Sellers (1986). The frozen lake and wetland albedos
    //     are set to constants (0.6 for visible beam, 0.4 for near-infrared)
    // (6) over the snow covered tile, the surface albedo is estimated by a linear
    //     combination of albedos for snow, canopy and bare soil (or lake, wetland, glacier).
    //
    // Original author : Yongjiu Dai, 09/15/1999; 08/30/2002
    //=======================================================================

    //use precision
    // use phycon_module, only : tfrz
    // implicit none

    //------------------------- Dummy Arguments -----------------------------
    // ground cover index
    //integer, INTENT(in) :: &
    //     itypwat     // land water type (0=soil, 1=urban or built-up, 2=wetland,
    //                 // 3=land ice, 4=deep lake, 5=shallow lake)

    //                 // parameters
    //real(r8), INTENT(in) :: &
    //     albsol,    &// soil albedo for different coloured soils [-]
    //     chil,      &// leaf angle distribution factor
    //     ref(2,2),  &// leaf reflectance (iw=iband, il=life and dead)
    //     tran(2,2), &// leaf transmittance (iw=iband, il=life and dead)
    //     fveg,      &// fractional vegetation cover [-]
    //     green,     &// green leaf fraction
    //     lai,       &// leaf area index (LAI+SAI) [m2/m2]
    //     sai,       &// stem area index (LAI+SAI) [m2/m2]

    //                 // variables
    //     coszen,    &// cosine of solar zenith angle [-]
    //     wt,        &// fraction of vegetation covered by snow [-]
    //     fsno,      &// fraction of soil covered by snow [-]
    //     ssw,       &// water volumetric content of soil surface layer [m3/m3]
    //     scv,       &// snow cover, water equivalent [mm]
    //     sag,       &// non dimensional snow age [-]
    //     tg          // ground surface temperature [K]

    //real(r8), INTENT(out) :: &
    //     alb(2,2),  &// averaged albedo [-]
    //     albg(2,2), &// albedo, ground
    //     albv(2,2), &// albedo, vegetation [-]
    //     ssun(2,2), &// sunlit canopy absorption for solar radiation
    //     ssha(2,2), &// shaded canopy absorption for solar radiation,
    //                 // normalized by the incident flux
    //     thermk,    &// canopy gap fraction for tir radiation
    //     extkb,     &// (k, g(mu)/mu) direct solar extinction coefficient
    //     extkd       // diffuse and scattered diffuse PAR extinction coefficient

    //-------------------------- Local variables ----------------------------
//	int             //
//		iw,        // wavelength (1=visible, 2=near-infrared)
//		id,        // 1=direct, 2=diffuse
//		k;           // looping indx
    int i,j;

    double 	age;       // factor to reduce visible snow alb due to snow age [-]
    double 		albg0;     // temporary varaiable [-]
    double 		albsno[2][2];// snow albedo [-]
//		albv0[2];  // vegetation albedo [-]
    double 		alwet;     // decrease in soil albedo due to wetness [-]
    double 		beta0;     // upscattering parameter for direct beam [-]
    double 		cff;       // snow alb correction factor for zenith angle > 60 [-]
    double 		conn;      // constant (=0.5) for visible snow alb calculation [-]
    double 		cons;      // constant (=0.2) for nir snow albedo calculation [-]
    double 		czen;      // cosine of solar zenith angle > 0 [-]
    double 		czf;       // solar zenith correction for new snow albedo [-]
    double 		dfalbl;    // snow albedo for diffuse nir radiation [-]
    double 		dfalbs;    // snow albedo for diffuse visible solar radiation [-]
    double 		dralbl;    // snow albedo for visible radiation [-]
    double 		dralbs;    // snow albedo for near infrared radiation [-]
    double 		fsol1;     // solar flux fraction for wavelength < 0.7 micron [-]
    double 		fsol2;     // solar flux fraction for wavelength > 0.7 micron [-]
    double 		lsai;      // leaf and stem area index (LAI+SAI) [m2/m2]
    double 		scat[2];   // single scattering albedo for vir/nir beam [-]
    double 		sl;        // factor that helps control alb zenith dependence [-]
    double 		snal0;     // alb for visible;incident on new snow (zen ang<60) [-]
    double 		snal1;     // alb for NIR; incident on new snow (zen angle<60) [-]
//		tdiffs;    // difference of air temperature and freezing temp [K]
//		tff;       // std::exp(-LSAI)
//		tffd;      // std::exp(-0.5*LSAI/czen)
//		ti;        // correction due to scattering
//		upscat;    // upward scattered fraction for direct beam [-]
    double 		tranc[2][2];// canopy transmittances for solar radiation
//		zkat[2],   // temporary
//		zkatd[2];    // temporary

    // ----------------------------------------------------------------------
    // 1. Initial set
    // ----------------------------------------------------------------------
    // division of solar flux for wavelength less or greater than 0.7 micron
    fsol1 = 0.5;      // shortwave
    fsol2 = 0.5;      // longwave

    // short and long wave albedo for new snow
    snal0 = 0.85;     // shortwave
    snal1 = 0.65;     // long wave

    // set initial leaf scattering reflectance. Note: "scat" may use different
    // value for different vegetation latter
    beta0 = 0.5;
    scat[0] = 0.15;
    scat[1] = 0.85;

    // ----------------------------------------------------------------------
    // set default soil and vegetation albedos and solar absorption
    for (i=0; i<2; i++)
        for (j=0; j<2; j++)
        {
            alb[i][j] = 0.0; // averaged
            albg[i][j] = 0.0; // ground
            albv[i][j] = 0.0; // vegetation
            albsno[i][j]=0.0; //set initial snow albedo
            ssun[i][j] = 0.0;
            ssha[i][j] = 0.0;
            tranc[i][j] = 0.0;
        }
    thermk = 1.e-3;
    extkb = 1.e-6;
    extkd = 0.718;

    lsai=lai+sai;
    if (coszen<=0.0)
        return;  //only do albedo when coszen > 0

    czen=max(coszen,0.001);

    // ----------------------------------------------------------------------
    // 2. albedo for snow cover.
    //    snow albedo depends on snow-age, zenith angle, and thickness
    //    of snow age gives reduction of visible radiation
    // ----------------------------------------------------------------------
    if (scv>0.)
    {
        cons = 0.2;
        conn = 0.5;
        sl  = 2.0 ;          //sl helps control albedo zenith dependence

        // correction for snow age
        age = 1.-1./(1.+sag) ;//correction for snow age
        dfalbs = snal0*(1.-cons*age);

        // czf corrects albedo of new snow for solar zenith
        cff    = ((1.+1./sl)/(1.+czen*2.*sl )- 1./sl);
        cff    = max(cff,0.);
        czf    = 0.4*cff*(1.-dfalbs);
        dralbs = dfalbs+czf;
        dfalbl = snal1*(1.-conn*age);
        czf    = 0.4*cff*(1.-dfalbl);
        dralbl = dfalbl+czf;

        albsno[0][0] = dralbs;
        albsno[1][0] = dralbl;
        albsno[0][1] = dfalbs;
        albsno[1][1] = dfalbl;
    }

    // ----------------------------------------------------------------------
    // 3. get albedo over land
    // ----------------------------------------------------------------------
    // 3.1 bare soil albedos, depends on moisture
    if (itypwat<=1)   // not wetland, permanent ice and water
    {
        alwet = max((11.-40.0*ssw),0.)*0.01;
        alwet = min(alwet,albsol);
        albg0 = albsol+alwet;
        albg[0][0] = albg0;
        albg[1][0] = 2.*albg0;
        albg[0][1] = albg[0][0];
        albg[1][1] = albg[1][0];//diffused albedos for bare soil
    }

    // 3.2 albedos for permanent ice sheet.
    else if (itypwat==3)        //permanent ice sheet
    {
        albg[0][0] = 0.8;
        albg[0][1] = 0.8;
        albg[1][0] = 0.55;
        albg[1][1] = 0.55;
    }

    // 3.3 albedo for wetland (swamps, rice paddies etc) and inland water
    else if ((itypwat==2)||(itypwat>=4))
    {
        albg0 = 0.05/(czen+0.15);
        albg[0][0] = albg0;
        albg[0][1] = albg0;
        albg[1][0] = albg0;
        albg[1][1] = albg0;
        if (tg<tfrz)             //frozen lake and wetland
        {
            albg[0][0] = 0.6;
            albg[0][1] = 0.6;
            albg[1][0] = 0.4;
            albg[1][1] = 0.4;
        }
    }

    // 3.4 correction due to snow cover
    for (i=0; i<2; i++)
        for (j=0; j<2; j++)
        {
            albg[i][j] = (1.-fsno)*albg[i][j] + fsno*albsno[i][j];
            alb[i][j] = albg[i][j];
        }

    // ----------------------------------------------------------------------
    // 4. canopy albedos : two stream approximation
    // ----------------------------------------------------------------------
    if (fveg>0.001)
    {
        twostream (chil,  ref,  tran,  green, lai,
                   sai   ,czen, albg,  albv,  tranc,
                   thermk,extkb,extkd, ssun,  ssha) ;
        for (i=0; i<2; i++)
            for (j=0; j<2; j++)
            {
                albv[i][j] = (1.-wt)*albv[i][j]+wt*albsno[i][j];
                alb[i][j] = (1.-fveg)*albg[i][j]+fveg*albv[i][j];
            }
    }
    //-----------------------------------------------------------------------
}
void CommonLandModel::albocean (double oro,double scv,double  coszrs,double alb[2][2])
{

    //-----------------------------------------------------------------------
    //
    // Compute surface albedos
    //
    // Computes surface albedos for direct/diffuse incident radiation for
    // two spectral intervals:
    //   s = 0.2-0.7 micro-meters
    //   l = 0.7-5.0 micro-meters
    //
    // Albedos specified as follows:
    //
    // Ocean           Uses solar zenith angle to compute albedo for direct
    //                 radiation; diffuse radiation values constant; albedo
    //                 independent of spectral interval and other physical
    //                 factors such as ocean surface wind speed.
    //
    // Ocean with      Surface albs specified; combined with overlying snow
    //   sea ice
    //
    // For more details , see Briegleb, Bruce P., 1992: Delta-Eddington
    // Approximation for Solar Radiation in the NCAR Community Climate Model,
    // Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
    //
    // yongjiu dai and xin-zhong liang (08/01/2001)
    //-----------------------------------------------------------------------

    // use precision
    //implicit none

    //------------------------------Arguments--------------------------------

    // real(r8), INTENT(in) :: oro       // /ocean(0)/seaice(2) flag
    // real(r8), INTENT(in) :: scv       // snow water equivalent) [mm]
    // real(r8), INTENT(in) :: coszrs    // Cosine solar zenith angle

    // real(r8), INTENT(out) :: alb(2,2) // srf alb for direct (diffuse) rad 0.2-0.7 micro-ms
    // Srf alb for direct (diffuse) rad 0.7-5.0 micro-ms

    //---------------------------Local variables-----------------------------

    double  frsnow;       // horizontal fraction of snow cover
    double  snwhgt;       // physical snow height
    double  rghsnw;       // roughness for horizontal snow cover fractn

    double  sasdir;       // snow alb for direct rad  0.2-0.7 micro-ms
    double  saldir;       // snow alb for direct rad  0.7-5.0 micro-ms
    double  sasdif;       // snow alb for diffuse rad  0.2-0.7 micro-ms
    double  saldif;       // snow alb for diffuse rad  0.7-5.0 micro-ms

    const double   asices = 0.70; // sea ice albedo for 0.2-0.7 micro-meters [-]
    const double   asicel = 0.50; // sea ice albedo for 0.7-5.0 micro-meters [-]
    const double   asnows = 0.95; // snow    albedo for 0.2-0.7 micro-meters [-]
    const double   asnowl = 0.70; // snow    albedo for 0.7-5.0 micro-meters

    //-----------------------------------------------------------------------
    // initialize all ocean/sea ice surface albedos to zero

    alb[0][0] = 0.0;
    alb[0][1] = 0.0;
    alb[1][0] = 0.0;
    alb[1][1] = 0.0;
    if (coszrs<=0.0)
        return;

    if ((int)(oro)==2)
    {
        alb[0][0] = asices;
        alb[1][0] = asicel;
        alb[0][1] = alb[0][0];
        alb[1][1] = alb[1][0];
        sasdif = asnows;
        saldif = asnowl;

        if (scv>0.0)
        {
            if (coszrs<0.5)
            {
                // zenith angle regime 1 ( coszrs < 0.5 ).
                // set direct snow albedos (limit to 0.98 max)
                sasdir = min(0.98,sasdif+(1.-sasdif)*0.5*(3./(1.+4.*coszrs)-1.));
                saldir = min(0.98,saldif+(1.-saldif)*0.5*(3./(1.+4.*coszrs)-1.));
            }
            else
                // zenith angle regime 2 ( coszrs >= 0.5 )
            {
                sasdir = asnows;
                saldir = asnowl;
            }

            // compute both diffuse and direct total albedos
            snwhgt = 20.*scv / 1000.;
            rghsnw = 0.25;
            frsnow = snwhgt/(rghsnw+snwhgt);
            alb[0][0] = alb[0][0]*(1.-frsnow) + sasdir*frsnow;
            alb[1][0] = alb[1][0]*(1.-frsnow) + saldir*frsnow;
            alb[0][1] = alb[0][1]*(1.-frsnow) + sasdif*frsnow;
            alb[1][1] = alb[1][1]*(1.-frsnow) + saldif*frsnow;
        }
    }

    // ice-free ocean albedos function of solar zenith angle only, and
    // independent of spectral interval:

    if ((int)oro==0)
    {
        alb[1][0] = 0.026/(pow(coszrs,1.7)+0.065)
                    + 0.15*(coszrs-0.1)*(coszrs-0.5)*(coszrs-1.) ;
        alb[0][0] = alb[1][0];
        alb[0][1] = 0.06;
        alb[1][1] = 0.06;
    }
}
double CommonLandModel::sign(double a, double b)
{
    if (b >= 0) return fabs(a);
    else return -fabs(a);
}
void CommonLandModel::combo (double &dz,   double &wliq, double &wice,double &t,double dz2,
                             double wliq2,double wice2,double t2 )
{

//=======================================================================
// Original author: Yongjiu Dai, September 15, 1999
//
// combines two elements and returns the following combined
// variabless: dz, t, wliq, wice.
// the combined temperature is based on the equation:
// the sum of the enthalpies of the two elements = that of the combined element.
//
//=======================================================================

    //use precision
    //use phycon_module, only : cpice, cpliq, hfus, tfrz
    //implicit none

//-------------------------- Dummy argument -----------------------------

    //real(r8), INTENT(in) :: dz2     // nodal thickness of 2 elements being combined [m]
    //real(r8), INTENT(in) :: wliq2   // liquid water of element 2 [kg/m2]
    //real(r8), INTENT(in) :: wice2   // ice of element 2 [kg/m2]
    //real(r8), INTENT(in) :: t2      // nodal temperature of element 2 [K]

    //real(r8), INTENT(inout) :: dz   // nodal thickness of 1 elements being combined [m]
    //real(r8), INTENT(inout) :: wliq // liquid water of element 1
    //real(r8), INTENT(inout) :: wice // ice of element 1 [kg/m2]
    //real(r8), INTENT(inout) :: t    // nodel temperature of elment 1 [K]

//----------------------- Local variables ------------------------------

    double   dzc  ;  // Total thickness of nodes 1 and 2 (dzc=dz+dz2).
    double   wliqc;  // Combined liquid water [kg/m2]
    double   wicec;  // Combined ice [kg/m2]
    double   tc   ;  // Combined node temperature [K]
    double   h    ;  // enthalpy of element 1 [J/m2]
    double   h2   ;  // enthalpy of element 2 [J/m2]
    double   hc   ;  // temporary

//-----------------------------------------------------------------------

    dzc = dz+dz2;
    wicec = (wice+wice2);
    wliqc = (wliq+wliq2);
    h   = (cpice*wice+cpliq*wliq)*(t-tfrz)+hfus*wliq;
    h2  = (cpice*wice2+cpliq*wliq2)*(t2-tfrz)+hfus*wliq2;

    hc = h + h2;
    if (hc < 0.0)
        tc = tfrz + hc/(cpice*wicec+cpliq*wliqc);
    else if (hc<=(hfus*wliqc))
        tc = tfrz;
    else
        tc = tfrz + (hc - hfus*wliqc)/(cpice*wicec+cpliq*wliqc);

    dz = dzc;
    wice = wicec;
    wliq = wliqc;
    t = tc;
}
void CommonLandModel::dewfraction (double sigf,double lai,double sai,double dewmx,double ldew,
                                   double &fwet,double &fdry)
{

//=======================================================================
// Original author: Yongjiu Dai, September 15, 1999
//
// determine fraction of foliage covered by water and
// fraction of foliage that is dry and transpiring
//
//=======================================================================

//use precision
//implicit none

    //real(r8), INTENT(in) :: sigf   // fraction of veg cover, excluding snow-covered veg [-]
    //real(r8), INTENT(in) :: lai    // leaf area index  [-]
    //real(r8), INTENT(in) :: sai    // stem area index  [-]
    //real(r8), INTENT(in) :: dewmx  // maximum allowed dew [0.1 mm]
    //real(r8), INTENT(in) :: ldew   // depth of water on foliage [kg/m2/s]

    //real(r8), INTENT(out) :: fwet  // fraction of foliage covered by water [-]
    //real(r8), INTENT(out) :: fdry  // fraction of foliage that is green and dry [-]

    //real(r8) lsai                  // lai + sai
    //real(r8) dewmxi                // inverse of maximum allowed dew [1/mm]
    //real(r8) vegt                  // sigf*lsai
//
//-----------------------------------------------------------------------
// Fwet is the fraction of all vegetation surfaces which are wet
// including stem area which contribute to evaporation
    double lsai ;                 // lai + sai
    double dewmxi;                // inverse of maximum allowed dew [1/mm]
    double vegt  ;                // sigf*lsai

    lsai = lai + sai;
    dewmxi  = 1.0/dewmx;
    vegt   =  sigf*lsai;

    fwet = 0;
    if (ldew > 0.0)
    {
        fwet = pow(((dewmxi/vegt)*ldew),0.666666666666);

// Check for maximum limit of fwet
        fwet = min(fwet,1.0);

    }

// fdry is the fraction of lai which is dry because only leaves can
// transpire. Adjusted for stem area which does not transpire
    fdry = (1.-fwet)*lai/lsai;
//
//-----------------------------------------------------------------------
}

void CommonLandModel::eroot(int lb,int    nl_soil,           double trsmx0,             double porsl[MAXSOILL],     double bsw[MAXSOILL],        double phi0[MAXSOILL],
                            double rootfr[MAXSOILL],  double dz[MAXSNL+MAXSOILL],double tss[MAXSNL+MAXSOILL],double wliq[MAXSNL+MAXSOILL],double rootr[MAXSOILL],
                            double &etrc,              double &rstfac)
{

//=======================================================================
// effective root fraction and maximum possible transpiration rate
// Original author : Yongjiu Dai, 08/30/2002
//=======================================================================

//  use precision
//  use phycon_module, only : tfrz
//  implicit none

//-----------------------Argument-----------------------------------------

    //integer, INTENT(in) :: nl_soil            // upper bound of array

    //real(r8), INTENT(in) :: trsmx0            // max transpiration for moist soil+100% veg.[mm/s]
    //real(r8), INTENT(in) :: porsl(1:nl_soil)  // soil porosity [-]
    //real(r8), INTENT(in) :: bsw(1:nl_soil)    // Clapp-Hornberger "B"
    //real(r8), INTENT(in) :: phi0(1:nl_soil)   // saturated soil suction (mm)
    //real(r8), INTENT(in) :: rootfr(1:nl_soil) // fraction of roots in a layer,
    //real(r8), INTENT(in) :: dz(1:nl_soil)     // layer thickness (m)
    //real(r8), INTENT(in) :: tss(1:nl_soil)    // soil/snow skin temperature (K)
    //real(r8), INTENT(in) :: wliq(1:nl_soil)   // liquid water (kg/m2)

    //real(r8), INTENT(out) :: rootr(1:nl_soil) // root resistance of a layer, all layers add to 1
    //real(r8), INTENT(out) :: etrc             // maximum possible transpiration rate (mm h2o/s)
    //real(r8), INTENT(out) :: rstfac           // factor of soil water stress for photosynthesis

//-----------------------Local Variables------------------------------

    double  roota ;            // accumulates root resistance factors
    double  rresis[MAXSOILL];            // soil water contribution to root resistance
    double  s_node ;           // vol_liq/porosity
    double  smpmax;            // wilting point potential in mm
    double  smp_node;          // matrix potential

    int i;                    // loop counter

//-----------------------End Variables list---------------------------

    // transpiration potential(etrc) and root resistance factors (rstfac)

    roota = 1.e-10 ;        // must be non-zero to begin
    // do i = 1, nl_soil
    for (i = 0; i<nl_soil; i++)
    {
        if ((tss[lb+i]>tfrz)&&(porsl[i]>=1.e-6))
        {
            smpmax = -1.5e+5;
            s_node = max(wliq[lb+i]/(1000.*dz[lb+i]*porsl[i]),0.001);
            s_node = min(1., s_node);
            smp_node = max(smpmax, -phi0[i]*pow(s_node,(-bsw[i])));
            rresis[i]   =(1.-smp_node/smpmax)/(1.+phi0[i]/smpmax) ;
            rootr[i] = rootfr[i]*rresis[i];
            roota = roota + rootr[i];
        }
        else
            rootr[i] = 0.0;
    }

    // normalize root resistances to get layer contribution to ET
    for (i = 0; i<nl_soil; i++)
        rootr[i] = rootr[i]/roota;

    // determine maximum possible transpiration rate
    etrc = trsmx0*roota;
    rstfac = roota;
}
void CommonLandModel::groundfluxes (double zlnd,  double zsno,  double hu,      double ht,     double hq,
                                    double us,    double vs,    double tm,      double qm,     double rhoair,
                                    double psrf,  double ur,    double thm,     double th,     double thv,
                                    double tg,    double qg,    double dqgdT,   double htvp,   double fsno,
                                    double sigf,  double &cgrnd,double &cgrndl, double &cgrnds,double &taux,
                                    double &tauy, double &fsena,double &fevpa,  double &fseng, double &fevpg,
                                    double &tref, double &qref, double &z0ma,   double &zol,   double &rib,
                                    double &ustar,double &qstar,double &tstar,  double &f10m,  double &fm,
                                    double &fh,   double &fq)
{

    ////=======================================================================
    //// this is the main subroutine to execute the calculation of thermal processes
    //// and surface fluxes
    ////
    //// Original author : Yongjiu Dai, 09/15/1999; 08/30/2002
    ////=======================================================================
    //
    //  use precision
    //  use phycon_module, only : cpair,vonkar,grav
    //  implicit none
    //
    ////----------------------- Dummy argument --------------------------------
    //  double, INTENT(in) :: &
    //        zlnd,     // roughness length for soil [m]
    //        zsno,     // roughness length for snow [m]
    //
    //        // atmospherical variables and observational height
    //        hu,       // observational height of wind [m]
    //        ht,       // observational height of temperature [m]
    //        hq,       // observational height of humidity [m]
    //        us,       // wind component in eastward direction [m/s]
    //        vs,       // wind component in northward direction [m/s]
    //        tm,       // temperature at agcm reference height [kelvin] [not used]
    //        qm,       // specific humidity at agcm reference height [kg/kg]
    //        rhoair,   // density air [kg/m3]
    //        psrf,     // atmosphere pressure at the surface [pa] [not used]
    //
    //        fsno,     // fraction of ground covered by snow
    //        sigf,     // fraction of veg cover, excluding snow-covered veg [-]
    //
    //        ur,       // wind speed at reference height [m/s]
    //        thm,      // intermediate variable (tm+0.0098*ht)
    //        th,       // potential temperature (kelvin)
    //        thv,      // virtual potential temperature (kelvin)
    //
    //        tg,       // ground surface temperature [K]
    //        qg,       // ground specific humidity [kg/kg]
    //        dqgdT,    // d(qg)/dT
    //        htvp       // latent heat of vapor of water (or sublimation) [j/kg]
    //
    //  double, INTENT(out) :: &
    //        taux,     // wind stress: E-W [kg/m/s**2]
    //        tauy,     // wind stress: N-S [kg/m/s**2]
    //        fsena,    // sensible heat from canopy height to atmosphere [W/m2]
    //        fevpa,    // evapotranspiration from canopy height to atmosphere [mm/s]
    //        fseng,    // sensible heat flux from ground [W/m2]
    //        fevpg,    // evaporation heat flux from ground [mm/s]
    //        cgrnd,    // deriv. of soil energy flux wrt to soil temp [w/m2/k]
    //        cgrndl,   // deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
    //        cgrnds,   // deriv of soil latent heat flux wrt soil temp [w/m**2/k]
    //        tref,     // 2 m height air temperature [kelvin]
    //        qref,     // 2 m height air humidity
    //
    //        z0ma,     // effective roughness [m]
    //        zol,      // dimensionless height (z/L) used in Monin-Obukhov theory
    //        rib,      // bulk Richardson number in surface layer
    //        ustar,    // friction velocity [m/s]
    //        tstar,    // temperature scaling parameter
    //        qstar,    // moisture scaling parameter
    //        f10m,     // integral of profile function for momentum at 10m
    //        fm,       // integral of profile function for momentum
    //        fh,       // integral of profile function for heat
    //        fq         // integral of profile function for moisture

    //------------------------ LOCAL VARIABLES ------------------------------
    //  long i;
    long niters, // maximum number of iterations for surface temperature
         iter,      // iteration index
         nmozsgn;     // number of times moz changes sign

    double
    beta,      // coefficient of conective velocity [-]
    displax,   // zero-displacement height [m]
    dth,       // diff of virtual temp. between ref. height and surface
    dqh,       // diff of humidity between ref. height and surface
    dthv,      // diff of vir. poten. temp. between ref. height and surface
    obu,       // monin-obukhov length (m)
    obuold,    // monin-obukhov length from previous iteration
    ram,       // aerodynamical resistance [s/m]
    rah,       // thermal resistance [s/m]
    raw,       // moisture resistance [s/m]
    raih,      // temporary variable [kg/m2/s]
    raiw,      // temporary variable [kg/m2/s]
    temp1,     // relation for potential temperature profile
    temp2,     // relation for specific humidity profile
    temp12m,   // relation for temperature at 2m
    temp22m,   // relation for specific humidity at 2m
    thvstar,   // virtual potential temperature scaling parameter
    um,        // wind speed including the stablity effect [m/s]
    wc,        // convective velocity [m/s]
    wc2,       // wc**2
    zeta,      // dimensionless height used in Monin-Obukhov theory
    zii,       // convective boundary height [m]
    zldis,     // reference height "minus" zero displacement heght [m]
    z0mg,      // roughness length over ground, momentum [m]
    z0hg,      // roughness length over ground, sensible heat [m]
    z0qg ;     // roughness length over ground, latent heat [m]

    //----------------------- Dummy argument --------------------------------
    // initial roughness length
    if (fsno > 0.0)
    {
        z0mg = zsno;
        z0hg = z0mg;
        z0qg = z0mg;
    }
    else
    {
        z0mg = zlnd;
        z0hg = z0mg;
        z0qg = z0mg;
    }
    // potential temperatur at the reference height
    beta = 1.0;      // -  (in computing W_*)
    zii = 1000.0;    // m  (pbl height)
    z0ma = z0mg;

    //-----------------------------------------------------------------------
    //     Compute sensible and latent fluxes and their derivatives with respect
    //     to ground temperature using ground temperatures from previous time step.
    //-----------------------------------------------------------------------
    // Initialization variables
    nmozsgn = 0;
    obuold = 0.0;

    dth   = thm-tg;
    dqh   = qm-qg;
    dthv  = dth*(1.+0.61*qm)+0.61*th*dqh;
    zldis = hu-0.0;
    moninobukini(ur,th,thm,thv,dth,dqh,dthv,zldis,z0mg,um,obu);

    // Evaluated stability-dependent variables using moz from prior iteration
    niters=6;

    //----------------------------------------------------------------
    // do iter = 1, niters         // begin stability iteration
    //----------------------------------------------------------------
    for (iter=1; iter<=niters; iter++)
    {
        displax = 0.0;
        moninobuk(hu,ht,hq,displax,z0mg,
                  z0hg,z0qg,obu,um,ustar,
                  temp1,temp2,temp12m,temp22m,f10m,
                  fm,fh,fq);
        tstar = temp1*dth;
        qstar = temp2*dqh;

        z0hg = z0mg/std::exp(0.13*pow((ustar*z0mg/1.5e-5),0.45));
        z0qg = z0hg;

        thvstar=tstar+0.61*th*qstar;
        zeta=zldis*vonkar*grav*thvstar/(ustar*ustar*thv);
        if (zeta >= 0.0)     //stable
            zeta = min(2.0,max(zeta,1.e-6));
        else                    //unstable
            zeta = max(-100.,min(zeta,-1.e-6));

        obu = zldis/zeta;

        if (zeta >= 0.)
            um = max(ur,0.1);
        else
        {
            wc = pow((-grav*ustar*thvstar*zii/thv),(1.0/3.0));
            wc2 = beta*beta*(wc*wc);
            um = std::sqrt(ur*ur+wc2);
        }
        if ((obuold*obu) < 0.0)
            nmozsgn = nmozsgn+1;
        if (nmozsgn >= 4)
            break;

        obuold = obu;

        //----------------------------------------------------------------
        //enddo                       // end stability iteration
        //----------------------------------------------------------------
    }
    // Get derivative of fluxes with repect to ground temperature
    ram    = 1./(ustar*ustar/um);
    rah    = 1./(temp1*ustar);
    raw    = 1./(temp2*ustar) ;

    raih   = (1.-sigf)*rhoair*cpair/rah;
    raiw   = (1.-sigf)*rhoair/raw;
    cgrnds = raih;
    cgrndl = raiw*dqgdT;
    cgrnd  = cgrnds + htvp*cgrndl;
    zol = zeta;
//	rib = min(5.,zeta*pow(vonkar,3)*ustar*ustar/(temp1*um*um));

    // correct 27, December, 2006
    rib = min(5.,zol*ustar*ustar/(vonkar*temp1*um*um));

    // surface fluxes of momentum, sensible and latent
    // using ground temperatures from previous time step
    taux   = -(1.-sigf)*rhoair*us/ram ;
    tauy   = -(1.-sigf)*rhoair*vs/ram;

    fseng  = -raih*dth;
    fevpg  = -raiw*dqh ;
    //fseng = max(fseng,0.0);
    //fevpg = max(fevpg,0.0);

    fsena  = fseng;
    fevpa  = fevpg;

    // 2 m height air temperature
    tref   = (1.-sigf)*(thm + temp1*dth * (1./temp12m - 1./temp1));
    qref   = (1.-sigf)*( qm + temp2*dqh * (1./temp22m - 1./temp2));
}
void CommonLandModel::groundtem(int    itypwat,              int     lb,                    int     nl_soil,              double  dtime,                  double  capr,
                                double  cnfac,               double  csol[MAXSOILL],        double  porsl[MAXSOILL],      double  dkmg[MAXSOILL],         double  dkdry[MAXSOILL],
                                double  dksatu[MAXSOILL],    double  sigf,                  double  dz[MAXSNL+MAXSOILL],  double  z[MAXSNL+MAXSOILL],     double  zi[MAXSNL+MAXSOILL+1],
                                double  tss[MAXSNL+MAXSOILL],double  wice[MAXSNL+MAXSOILL], double  wliq[MAXSNL+MAXSOILL], double  &scv,                    double  &snowdp,
                                double  frl,                 double  dlrad,                 double  sabg,                 double  fseng,                  double  fevpg,
                                double  cgrnd,               double  htvp,                  double  emg,                  int     imelt[MAXSNL+MAXSOILL], double & sm,
                                double  &xmf,                 double  fact[MAXSNL+MAXSOILL])
{

//=======================================================================
// Snow and soil temperatures
// o The volumetric heat capacity is calculated as a linear combination
//   in terms of the volumetric fraction of the constituent phases.
// o The thermal conductivity of soil is computed from
//   the algorithm of Johansen (as reported by Farouki 1981), and of snow is from
//   the formulation used in SNTHERM (Jordan 1991).
// o Boundary conditions:
//   F = Rnet - Hg - LEg (top),  F= 0 (base of the soil column).
// o Soil / snow temperature is predicted from heat conduction
//   in 10 soil layers and up to 5 snow layers.
//   The thermal conductivities at the interfaces between two neighbor layers
//   (j, j+1) are derived from an assumption that the flux across the interface
//   is equal to that from the node j to the interface and the flux from the
//   interface to the node j+1. The equation is solved using the Crank-Nicholson
//   method and resulted in a tridiagonal system equation.
//
// Phase change (see meltf.F90)
//
// Original author : Yongjiu Dai, 09/15/1999; 08/30/2002
//=======================================================================

    //use precision
    //use phycon_module, only : stefnc
    //implicit none

// integer, INTENT(in) :: lb           //lower bound of array
// integer, INTENT(in) :: nl_soil      //upper bound of array
// integer, INTENT(in) :: itypwat      //land water type (0=soil,1=urban or built-up,2=wetland,
//                                     //3=land ice, 4=deep lake, 5=shallow lake)
//    real(r8), INTENT(in) :: dtime    //model time step [second]
//    real(r8), INTENT(in) :: capr     //tuning factor to turn first layer T into surface T
//    real(r8), INTENT(in) :: cnfac    //Crank Nicholson factor between 0 and 1

//    real(r8), INTENT(in) :: csol(1:nl_soil)  //heat capacity of soil solids [J/(m3 K)]
//    real(r8), INTENT(in) :: porsl(1:nl_soil) //soil porosity [-]
//    real(r8), INTENT(in) :: dkmg(1:nl_soil)  //thermal conductivity of soil minerals [W/m-K]
//    real(r8), INTENT(in) :: dkdry(1:nl_soil) //thermal conductivity of dry soil [W/m-K]
//    real(r8), INTENT(in) :: dksatu(1:nl_soil)//thermal conductivity of saturated soil [W/m-K]

//    real(r8), INTENT(in) :: sigf     //fraction of veg cover, excluding snow-covered veg [-]
//    real(r8), INTENT(in) :: dz[MAXSNL+MAXSOILL]   //layer thickiness [m]
//    real(r8), INTENT(in) :: z [MAXSNL+MAXSOILL]   //node depth [m]
//    real(r8), INTENT(in) :: zi(lb-1:nl_soil) //interface depth [m]

//    real(r8), INTENT(in) :: sabg     //solar radiation absorbed by ground [W/m2]
//    real(r8), INTENT(in) :: frl      //atmospheric infrared (intwave) radiation [W/m2]
//    real(r8), INTENT(in) :: dlrad    //downward intwave radiation blow the canopy [W/m2]
//    real(r8), INTENT(in) :: fseng    //sensible heat flux from ground [W/m2]
//    real(r8), INTENT(in) :: fevpg    //evaporation heat flux from ground [mm/s]
//    real(r8), INTENT(in) :: cgrnd    //deriv. of soil energy flux wrt to soil temp [w/m2/k]
//    real(r8), INTENT(in) :: htvp     //latent heat of vapor of water (or sublimation) [j/kg]
//    real(r8), INTENT(in) :: emg      //ground emissivity (0.97 for snow,

// real(r8), INTENT(inout) :: tss [MAXSNL+MAXSOILL] //soil temperature [K]
// real(r8), INTENT(inout) :: wice[MAXSNL+MAXSOILL] //ice lens [kg/m2]
// real(r8), INTENT(inout) :: wliq[MAXSNL+MAXSOILL] //liqui water [kg/m2]
// real(r8), INTENT(inout) :: scv      //snow cover, water equivalent [mm, kg/m2]
// real(r8), INTENT(inout) :: snowdp   //snow depth [m]

//   real(r8), INTENT(out) :: sm       //rate of snowmelt [kg/(m2 s)]
//   real(r8), INTENT(out) :: xmf      //total latent heat of phase change of ground water
//   real(r8), INTENT(out) :: fact[MAXSNL+MAXSOILL] //used in computing tridiagonal matrix
//integer, INTENT(out) :: imelt[MAXSNL+MAXSOILL]    //flag for melting or freezing [-]


//------------------------ local variables ------------------------------
    double  cv[MAXSNL+MAXSOILL];     // heat capacity [J/(m2 K)]
    double  tk[MAXSNL+MAXSOILL];     // thermal conductivity [W/(m K)]

    double  fn[MAXSNL+MAXSOILL];   // heat diffusion through the layer interface [W/m2]
    double  fn1[MAXSNL+MAXSOILL];   // heat diffusion through the layer interface [W/m2]
    double  dzm;                // used in computing tridiagonal matrix
    double  dzp;                // used in computing tridiagonal matrix

    double  tssbef[MAXSNL+MAXSOILL]; // soil/snow temperature before update
    double  hs;                 // net energy flux into the surface (w/m2)
    double  dhsdT;              // d(hs)/dT
    double  brr[MAXSNL+MAXSOILL];    // temporay set

    //double  at[MAXSNL+MAXSOILL];     //"a" vector for tridiagonal matrix
    //double  bt[MAXSNL+MAXSOILL];     //"b" vector for tridiagonal matrix
    //double  ct[MAXSNL+MAXSOILL];     //"c" vector for tridiagonal matrix
    //double  rt[MAXSNL+MAXSOILL];     //"r" vector for tridiagonal solution

    double *at    = new double[nl_soil+MAXSNL-lb];
    double *bt    = new double[nl_soil+MAXSNL-lb];
    double *ct    = new double[nl_soil+MAXSNL-lb];
    double *rt    = new double[nl_soil+MAXSNL-lb];
    double *tss_1 = new double[nl_soil+MAXSNL-lb];
    int i,j;
    int t_flag = 1;

//=======================================================================
    //cout<<fseng<<" "<<fevpg<<" "<<cgrnd<<endl;
// heat capacity
    hCapacity (itypwat,lb,nl_soil,csol,porsl,
               wice,wliq,scv,dz,cv);
// thermal conductivity
    hConductivity (itypwat,lb,nl_soil,dkmg,dkdry,
                   dksatu,porsl,dz,z,zi,
                   tss,wice,wliq,tk);
    //for(int i= lb;i<MAXSNL+MAXSOILL;i++)
    //	  cout<<i<<"  "<<cv[i]<<"     "<<tk[i]<<endl;
// net ground heat flux into the surface and its temperature derivative
    //cout<<emg*stefnc*pow(tss[lb],4)<<" "<<fevpg*htvp<<endl;
    hs = sabg + dlrad
         + (1.-sigf)*emg*frl - emg*stefnc*pow(tss[lb],4)
         - (fseng+fevpg*htvp) ;
    dhsdT = - cgrnd - 4.*emg * stefnc * pow(tss[lb],3);
    for (i= lb; i<MAXSNL+nl_soil; i++)
    {
        tssbef[i] = tss[i];
    }

    j       = lb;
    fact[j] = dtime / cv[j]
              * dz[j] / (0.5*(z[j]-zi[j]+capr*(z[j+1]-zi[j])));

    //cout<<fact[lb]<<endl;
    for (j = lb + 1; j<MAXSNL+nl_soil; j++)
        fact[j] = dtime/cv[j];
    for (j = lb; j<MAXSNL+nl_soil - 1; j++)
    {
        //cout<<j<<"  "<<tk[j]<<" "<<tss[j+1]<<" "<<tss[j]<<" "<<z[j+1]<<"  "<<z[j]<<endl;

        fn[j] = tk[j]*(tss[j+1]-tss[j])/(z[j+1]-z[j]);
        //cout<<j<<" "<<fn[j]<<endl;
    }
    fn[MAXSNL+nl_soil-1] = 0.0;



// set up vector r and vectors a, b, c that define tridiagonal matrix
    j     = lb;
    dzp   = z[j+1]-z[j];
    at[j-lb] = 0.;
    bt[j-lb] = 1+(1.-cnfac)*fact[j]*tk[j]/dzp-fact[j]*dhsdT;
    ct[j-lb] =  -(1.-cnfac)*fact[j]*tk[j]/dzp;
    rt[j-lb] = tss[j] + fact[j]*( hs - dhsdT*tss[j] + cnfac*fn[j] );
    for (j = lb+1; j<MAXSNL+nl_soil - 1; j++)
    {
        dzm   = (z[j]-z[j-1]);
        dzp   = (z[j+1]-z[j]);
        at[j-lb] =   - (1.-cnfac)*fact[j]* tk[j-1]/dzm;
        bt[j-lb] = 1.+ (1.-cnfac)*fact[j]*(tk[j]/dzp + tk[j-1]/dzm);
        ct[j-lb] =   - (1.-cnfac)*fact[j]* tk[j]/dzp;
        rt[j-lb] = tss[j] + cnfac*fact[j]*( fn[j] - fn[j-1] );
    }
    j     =  MAXSNL+nl_soil-1;
    dzm   = (z[j]-z[j-1]);
    at[j-lb] =   - (1.-cnfac)*fact[j]*tk[j-1]/dzm;
    bt[j-lb] = 1.+ (1.-cnfac)*fact[j]*tk[j-1]/dzm;
    ct[j-lb] = 0.;
    rt[j-lb] = tss[j] - cnfac*fact[j]*fn[j-1];

    for (j = lb; j<MAXSNL+nl_soil; j++)   tss_1[j-lb] = tss[j];
    //for(j = 0; j<10;j++)
    //cout<<j<<" "<<at[j]<<" "<<bt[j]<<" "<<ct[j]<<" "<<rt[j]<<endl;
// solve for tss
    tridia (nl_soil+(MAXSNL-lb),at ,bt ,ct ,rt ,tss_1);
    for (j = lb; j<MAXSNL+nl_soil; j++)   tss[j] = tss_1[j-lb];
    /* for (j = lb; j<MAXSNL+nl_soil;j++)
     {
         if (fabs(tss_1[j-lb]-tss[j])<15)
         {
             tss[j] = tss_1[j-lb];
         }

     }*/

    delete []at;
    delete []bt;
    delete []ct;
    delete []rt;
    delete []tss_1;
//=======================================================================
// melting or freezing
//=======================================================================

    for (j = lb; j<MAXSNL+MAXSOILL - 1; j++)
        fn1[j] = tk[j]*(tss[j+1]-tss[j])/(z[j+1]-z[j]);
    fn1[MAXSNL+MAXSOILL-1] = 0.;

    j = lb;
    brr[j] = cnfac*fn[j] + (1.-cnfac)*fn1[j];

    for (j = lb+1; j<MAXSNL+MAXSOILL; j++)
        brr[j] = cnfac*(fn[j]-fn[j-1]) + (1.-cnfac)*(fn1[j]-fn1[j-1]);
    meltf ( lb, nl_soil, dtime, fact, brr,
            hs, dhsdT,  tssbef, tss, wliq,
            wice, imelt,scv, snowdp, sm, xmf );
//-----------------------------------------------------------------------
}
void CommonLandModel::hCapacity (int    itypwat,        int lb,                       int   nl_soil1,                 double csol[MAXSOILL],
                                 double porsl[MAXSOILL],double wice[MAXSNL+MAXSOILL], double wliq[MAXSNL+MAXSOILL],   double scv,
                                 double dz[MAXSNL+MAXSOILL],           double cv[MAXSNL+MAXSOILL])
{


//-----------------------------------------------------------------------
// Original author : Yongjiu Dai, September 15, 1999
//
// calculation of heat capacities of snow / soil layers
// the volumetric heat capacity is calculated as a linear combination
// in terms of the volumetric fraction of the constituent phases.
//-----------------------------------------------------------------------

    //use precision
// use phycon_module, only : cpice,cpliq
//implicit none

    //integer, INTENT(in) :: lb           // lower bound of array
    //integer, INTENT(in) :: nl_soil1      // upper bound of array
    //integer, INTENT(in) :: itypwat      // land water type (0=soil, 1=urban, 2=wetland,
    //real(r8), INTENT(in) :: csol(1:nl_soil1)  // heat capacity of soil soilds [J/(m3 K)]
    //real(r8), INTENT(in) :: porsl(1:nl_soil1) // soil porosity
    //real(r8), INTENT(in) :: wice(lb:nl_soil1) // ice lens [kg/m2]
    //real(r8), INTENT(in) :: wliq(lb:nl_soil1) // liqui water [kg/m2]
    //real(r8), INTENT(in) :: dz(lb:nl_soil1)   // layer thickiness [m]
    //real(r8), INTENT(in) :: scv              // snow water equivalent [mm]
    //real(r8), INTENT(out) :: cv(lb:nl_soil1)  // heat capacity [J/(m2 K)]

//-----------------------------------------------------------------------
    int i;
// Snow heat capacity


    // Soil heat capacity, which from de Vires (1963)
    if (itypwat<=1) // soil ground
    {
        for (i = 0; i<nl_soil1; i++)
        {
            //cout<<csol[i]<<" "<<porsl[i]<<" "<<dz[MAXSNL+i]<<" "<<wice[MAXSNL+i]<<" "<<wliq[MAXSNL+i]<<endl;
            cv[MAXSNL+i] = csol[i]*(1.-porsl[i])*dz[MAXSNL+i] + wice[MAXSNL+i]*cpice + wliq[MAXSNL+i]*cpliq;
            //cout<<cv[MAXSNL+i]<<" "<<csol[i]<<" "<<porsl[i]<<" "<<dz[MAXSNL+i]<<" "<<wice[MAXSNL+i]<<" "<<wliq[MAXSNL+i]<<endl;
        }
        // getchar();

    }
    else               // wet land or glacier
    {
        for (i = 0; i<nl_soil1; i++)
        {
            cv[MAXSNL+i] = wice[MAXSNL+i]*cpice + wliq[MAXSNL+i]*cpliq;
            //cout<<cv[MAXSNL+i]<<" "<<wice[MAXSNL+i]<<" "<<wliq[MAXSNL+i]<<endl;
        }
    }
    if ((lb==MAXSNL)&&(scv>0.0))
        cv[MAXSNL] = cv[MAXSNL] + cpice*scv;

    if (lb<MAXSNL)
    {
        for (i =lb; i<MAXSNL; i++)
        {
            cv[i] = cpliq*wliq[i] + cpice*wice[i];
            //cout<<i<<" "<<cv[i]<<" "<<wliq[i]<<" "<<wice[i]<<endl;
        }
    }
    //for(int i=0;i<10;i++)
    //  cout<<cv[i+5]<<" "<<csol[i]<<" "<<porsl[i]<<" "<<dz[i+5]<<" "<<wice[i+5]<<" "<<wliq[i+5]<<endl;
}
void CommonLandModel::hConductivity (int   itypwat,                int    lb,                     int    nl_soil1,              double dkmg[MAXSOILL],       double dkdry[MAXSOILL],
                                     double dksatu[MAXSOILL],      double porsl[MAXSOILL],        double dz[MAXSNL+MAXSOILL],   double z[MAXSNL+MAXSOILL],   double zi[MAXSNL+MAXSOILL+1],
                                     double tss[MAXSNL+MAXSOILL],  double wice[MAXSNL+MAXSOILL],  double wliq[MAXSNL+MAXSOILL], double tk[MAXSNL+MAXSOILL])
{

    //-----------------------------------------------------------------------
    // Original author : Yongjiu Dai, September 15, 1999
    //
    // calculation of thermal conductivities of snow / soil layers
    // The thermal conductivity of soil is computed from
    // the algorithm of Johansen (as reported by Farouki 1981), and of snow is from
    // the formulation used in SNTHERM (Jordan 1991).
    //
    // The thermal conductivities at the interfaces between two neighbor layers
    // (j, j+1) are derived from an assumption that the flux across the interface
    // is equal to that from the node j to the interface and the flux from the
    // interface to the node j+1.
    //-----------------------------------------------------------------------

    //use precision
    //use phycon_module, only : denh2o,denice,tfrz,tkwat,tkice,tkair
    //implicit none

    //integer, INTENT(in) :: lb            // lower bound of array
    //integer, INTENT(in) :: nl_soil1       // upper bound of array
    //integer, INTENT(in) :: itypwat       // land water type (0=soil, 1=urban, 2=wetland,
    //                                     // 3=land ice, 4=deep lake, 5=shallow lake)
    //real(r8), INTENT(in) ::   dkmg(1:nl_soil1)  // dkm**dmvol, where dkm is the mineral cond.
    //real(r8), INTENT(in) ::  dkdry(1:nl_soil1)  // thermal conductivity for dry soil [W/m-K]
    //real(r8), INTENT(in) :: dksatu(1:nl_soil1)  // Thermal conductivity of saturated soil [W/m-K]
    //real(r8), INTENT(in) ::  porsl(1:nl_soil1)  // fractional volume between soil grains=1.-dmvol
    //real(r8), INTENT(in) ::   dz(lb:nl_soil1)   // layer thickiness [m]
    //real(r8), INTENT(in) ::    z(lb:nl_soil1)   // node depth [m]
    //real(r8), INTENT(in) ::   zi(lb-1:nl_soil1) // interface depth [m]
    //real(r8), INTENT(in) ::  tss(lb:nl_soil1)   // Nodal temperature [K]
    //real(r8), INTENT(in) :: wice(lb:nl_soil1)   // ice lens [kg/m2]
    //real(r8), INTENT(in) :: wliq(lb:nl_soil1)   // liqui water [kg/m2]

    //real(r8), INTENT(out) :: tk(lb:nl_soil1)    // thermal conductivity [W/(m K)]

    // local
    double  rhosnow;      // partitial density of water (ice + liquid)
    double  dksat;        // thermal conductivity for saturated soil (j/(k s m))
    double  dke ;         // kersten number
    double  fl ;          // fraction of liquid or unfrozen water to total water
    double  satw ;        // relative total water content of soil.
    double  thk[MAXSNL+MAXSOILL];  // thermal conductivity of layer

    int i;

    dke = 0.0;
    dksat = 0.0;
    //-----------------------------------------------------------------------
    // Thermal conductivity of soil from Farouki (1981),
    for (i=0; i<nl_soil1; i++)
    {
        if (itypwat<=1)        //soil ground
        {
            thk[MAXSNL+i] = dkdry[i];       //rock or dry soil
            //cout<<i<<" thk  "<<thk[MAXSNL+i]<<endl;
            if (porsl[i]>1.e-05)
            {
                satw = (wliq[MAXSNL+i]/denh2o+wice[MAXSNL+i]/denice)/(dz[MAXSNL+i]*porsl[i]);
                satw = min(1., satw);
                //cout<<wliq[MAXSNL+i]<<" "<<wice[MAXSNL+i]<<" "<<dz[MAXSNL+i]<<" "<<porsl[i]<<endl;
                if (satw>1e-6)
                {
                    fl = wliq[MAXSNL+i]/(wice[MAXSNL+i]+wliq[MAXSNL+i]);
                    if (tss[MAXSNL+i]>= tfrz)      // Unfrozen soil
                    {
                        dke = std::log10(satw) + 1.0;
                        dke = max(dke, 0.);
                        //dke = min(dke, 0.);
                        dksat = dksatu[i];
                    }
                    else                          // Frozen soil
                    {
                        dke = satw;
                        dksat = dkmg[i]*pow(0.249,(fl*porsl[i]))*pow(2.29,porsl[i]);
                    }
                    thk[MAXSNL+i] = dke*dksat + (1.-dke)*dkdry[i];
                }
            }

        }
        else           // wetland or glacier
        {
            thk[MAXSNL+i] = tkwat;
            if (tss[MAXSNL+i]<tfrz)
                thk[MAXSNL+i] = tkice;
        }
    }

    // Thermal conductivity of snow, which from Jordan (1991) pp. 18
    if (lb<MAXSNL)
    {
        for (i = lb; i<MAXSNL; i++)
        {

            rhosnow = (wice[i]+wliq[i])/dz[i];
            //cout<<wice[i]<<" "<<wliq[i]<<" "<<dz[i]<<" "<<rhosnow<<endl;
            //cout<<tkair<<" "<<rhosnow<<" "<<tkice<<endl;
            thk[i] = tkair+(7.75e-5*rhosnow+1.105e-6*rhosnow*rhosnow)*(tkice-tkair);
            //cout<<i<<" "<<wice[i]<<" "<<wliq[i]<<" "<<dz[i]<<" "<<thk[i]<<endl;
        }
    }
    // Thermal conductivity at the layer interface
    for (i = lb; i<MAXSNL+MAXSOILL-1; i++)
        // the following consideration is try to avoid the snow conductivity
        // to be dominant in the thermal conductivity of the interface.
        // Because when the distance of bottom snow node to the interfacee
        // is larger than that of interface to top soil node,
        // the snow thermal conductivity will be dominant, and the result is that
        // lees heat tranfer between snow and soil
    {
        //cout<<i<<"  "<<z[i+1]<<"  "<<zi[i+1]<<"  "<<z[i]<<endl;
        if ((i==MAXSNL-1)&&((z[i+1]-zi[i+1])<(zi[i+1]-z[i])))
        {
            tk[i] = 2.*thk[i]*thk[i+1]/(thk[i]+thk[i+1]);
            tk[i] = max(0.5*thk[i+1],tk[i]);
        }
        else
            tk[i] = thk[i]*thk[i+1]*(z[i+1]-z[i])
                    /(thk[i]*(z[i+1]-zi[i+1])+thk[i+1]*(zi[i+1]-z[i]));
        //cout<<i<<"  "<<thk[i]<<" "<<thk[i+1]<<" "<<z[i+1]<<" "<<z[i]<<" "<<zi[i+1]<<endl;
        //cout<<i<<" "<<tk[i]<<endl;
    }
    tk[MAXSNL+MAXSOILL-1] = 0.0;
}
void CommonLandModel::lai_empirical(int    ivt, int    nl_soil,double rootfr[MAXSOILL],double t[MAXSOILL+MAXSNL],double &lai,
                                    double &sai, double &fveg,   double &green)

//-----------------------------------------------------------------------
// provides leaf and stem area parameters
// Original author : Yongjiu Dai, 08/31/2002
//-----------------------------------------------------------------------

//use precision
//implicit none

//integer, intent(in)  :: ivt      //land cover type
//integer, intent(in)  :: nl_soil  //number of soil layers

// real(r8), intent(in)  :: rootfr(1:nl_soil)  //root fraction
// real(r8), intent(in)  :: t(1:nl_soil)  //soil temperature
// real(r8), intent(out) :: lai     //leaf area index
// real(r8), intent(out) :: sai     //Stem area index
// real(r8), intent(out) :: fveg    //fractional cover of vegetation
// real(r8), intent(out) :: green   //greenness
{

//local variable
    double  f;      //
    double  roota;  //accumulates root fraction
    int     jrt;     //number of soil layers with 90% root fraction
    int     j;       //number of soil layers

//-----------------------------------------------------------------------
// Maximum fractional cover of vegetation [-]
    double  vegc[24]= {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
                       1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,
                       1.0,1.0,0.0,1.0,1.0,1.0,0.0,0.0
                      };

// Maximum leaf area index, the numbers are based on the data of
// "worldwide histrorical estimates of leaf area index, 1932-2000" :
// http://www.daac.ornl.gov/global_vegetation/HistoricalLai/data"
    double  xla[24]= {1.50, 3.29, 4.18, 3.50, 2.50, 3.60, 2.02, 1.53,
                      2.00, 0.85, 4.43, 4.42, 4.56, 3.95, 4.50, 0.00,
                      4.00, 3.63, 0.00, 0.64, 1.60, 1.00, 0.00, 0.00
                     };

// Minimum leaf area index
//  real(r8), dimension(24), parameter :: &
// xla0=(/1.00, 0.50, 0.50, 0.50, 1.00, 0.50, 0.50, 0.50, &
//        0.50, 0.50, 0.50, 0.50, 4.00, 4.00, 4.00, 0.00, &
//        3.00, 3.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /)
    double xla0[24]= {1.00, 0.50, 0.50, 0.50, 1.00, 0.50, 0.50, 0.50,
                      0.50, 0.30, 0.50, 0.50, 4.00, 4.00, 4.00, 0.00,
                      3.00, 3.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00
                     };

// Stem area index [-]
//  real(r8), dimension(24), parameter :: &
// sai0=(/0.50, 0.50, 0.50, 0.50, 1.00, 1.00, 2.00, 1.00, &
//        1.00, 1.00, 2.00, 2.00, 2.00, 2.00, 2.00, 0.00, &
//        2.00, 2.00, 0.00, 0.50, 0.50, 0.50, 0.00, 0.00 /)
    double sai0[24]= {0.20, 0.20, 0.30, 0.30, 0.50, 0.50, 1.00, 0.50,
                      1.00, 0.50, 2.00, 2.00, 2.00, 2.00, 2.00, 0.00,
                      2.00, 2.00, 0.00, 0.10, 0.10, 0.10, 0.00, 0.00
                     };

//-----------------------------------------------------------------------
    roota = 0.0;
    jrt = 0;
    for (j=0; j<nl_soil; j++)
    {
        roota = roota + rootfr[j];
        if (roota>0.9)
        {
            jrt = j;
            break;
        }
    }

// Adjust leaf area index for seasonal variation

    f = max(0.0,1.-0.0016*max(298.-t[MAXSNL+jrt],0.0)*max(298.-t[MAXSNL+jrt],0.0));

    //利用MODIS LAI 产品，这里不再由经验公式计算LAI
    lai = xla[ivt-1] + (xla0[ivt-1]-xla[ivt-1])*(1.0-f);

// Sum leaf area index and stem area index
    sai = sai0[ivt-1];

// Fractional vegetation cover
    fveg = vegc[ivt-1];
    green = 0.0;
    if (fveg > 0.0)
        green = 1.0;
}
void CommonLandModel::lake (int    nl_lake ,              int    itypwat ,                double dlat,                 double dtime                 , double zlak[MAXSNL+MAXSOILL],
                            double dzlak[MAXSNL+MAXSOILL],double zilak[MAXSNL+MAXSOILL+1],double hu  ,                   double ht,                     double hq ,
                            double us   ,                 double vs   ,                   double tm  ,                 double qm ,                    double prc,
                            double prl ,                  double rhoair  ,                double psrf,                 double sabg  ,                 double frl,
                            double &tg ,                   double tlak[MAXSNL+MAXSOILL],   double wliq[MAXSNL+MAXSOILL],  double wice[MAXSNL+MAXSOILL]  ,double &scv,
                            double &snowdp,                double &trad,                    double &tref,                 double &qref ,                  double &taux,
                            double &tauy ,                 double &fsena,                   double &fevpa,                double &lfevpa,                 double &fseng,
                            double &fevpg,                 double &olrg ,                   double &fgrnd,                double &tcrit ,                 double &emis ,
                            double &z0ma ,                 double &zol ,                    double &rib ,                 double &ustar ,                 double &qstar,
                            double &tstar,                 double &u10m ,                   double &v10m,                 double &f10m  ,                 double &fm ,
                            double &fh   ,                 double &fq)
{
// ------------------------ code history ---------------------------
// purpose:           lake temperature and snow on frozen lake
// initial author:    Gordon Bonan
// revised:           Yongjiu Dai, September 15, 1999
//
// ------------------------ notes ----------------------------------
// calculate lake temperatures from one-dimensional thermal
// stratification model based on eddy diffusion concepts to
// represent vertical mixing of heat
//
// d ts    d            d ts     1 ds
// ---- = -- [(km + ke) ----] + -- --
//  dt    dz             dz     cw dz

// where: ts = temperature (kelvin)
//         t = time (s)
//         z = depth (m)
//        km = molecular diffusion coefficient (m**2/s)
//        ke = eddy diffusion coefficient (m**2/s)
//        cw = heat capacity (j/m**3/kelvin)
//         s = heat source term (w/m**2)

// there are two types of lakes:
//    deep lakes are 50 m. shallow lakes are 10 m deep.
//    for unfrozen deep lakes:    ke > 0 and    convective mixing
//    for unfrozen shallow lakes: ke = 0 and no convective mixing

// use crank-nicholson method to set up tridiagonal system of equations to
// solve for ts at time n+1, where the temperature equation for layer i is
// r_i = a_i [ts_i-1] n+1 + b_i [ts_i] n+1 + c_i [ts_i+1] n+1

// the solution conserves energy as

// cw*([ts(  1)] n+1 - [ts(  1)] n)*dz(  1)/dt + ... +
// cw*([ts(nl_lake)] n+1 - [ts(nl_lake)] n)*dz(nl_lake)/dt = fin

// where
// [ts] n   = old temperature (kelvin)
// [ts] n+1 = new temperature (kelvin)
// fin      = heat flux into lake (w/m**2)
//          = beta*sabg+frl-olrg-fsena-lfevpa-hm + phi(1) + ... + phi(nl_lake)

// -----------------------------------------------------------------

    // use precision
    // use phycon_module, only : tfrz,hvap,hfus,hsub,tkwat,tkice,stefnc,&
    //vonkar,grav,cpliq,cpair,denh2o,rgas
    // implicit none

// ------------------------ input/output variables -----------------

// integer, INTENT(in) :: &
//       nl_lake,   &//number of soil layers
//       itypwat     //land water type (4=deep lake, 5=shallow lake)

// real(r8), INTENT(in) :: &
//       dlat,      &//latitude (radians)
//       dtime,     &//time step (s)
//       tcrit,     &//critical temp. to determine rain or snow
//       dzlak(nl_lake),&//soil layer thickness (m)
//       zlak(nl_lake), &//depth (m)
//       zilak(0:nl_lake), &//
//       hu,        &//observational height of wind [m]
//       ht,        &//observational height of temperature [m]
//       hq,        &//observational height of humidity [m]
//       us,        &//wind component in eastward direction [m/s]
//       vs,        &//wind component in northward direction [m/s]
//       tm,        &//temperature at agcm reference height [kelvin]
//       qm,        &//specific humidity at agcm reference height [kg/kg]
//       prc,       &//convective precipitation [mm/s]
//       prl,       &//large scale precipitation [mm/s]
//       rhoair,    &//density air [kg/m3]
//       psrf,      &//atmosphere pressure at the surface [pa]
//       sabg,      &//solar radiation absorbed by ground [W/m2]
//       frl         //atmospheric infrared (longwave) radiation [W/m2]

// real(r8), INTENT(inout) :: &
//       tg,        &//surface temperature (kelvin)
//       tlak(nl_lake), &//lake temperature (kelvin)
//       wliq(nl_lake), &//
//       wice(nl_lake), &//
//       scv,       &//snow water equivalent [mm]
//       snowdp      //snow depth [mm]

// real(r8), INTENT(out) :: &
//       taux,      &//wind stress: E-W [kg/m/s**2]
//       tauy,      &//wind stress: N-S [kg/m/s**2]
//       fsena,     &//sensible heat from canopy height to atmosphere [W/m2]
//       fevpa,     &//evapotranspiration from canopy height to atmosphere [mm/s]
//       lfevpa,    &//latent heat flux from canopy height to atmosphere [W/m2]
//       fseng,     &//sensible heat flux from ground [W/m2]
//       fevpg,     &//evaporation heat flux from ground [mm/s]
//       olrg,      &//outgoing long-wave radiation from ground+canopy
//       fgrnd,     &//ground heat flux [W/m2]

//       tref,      &//2 m height air temperature [kelvin]
//       qref,      &//2 m height air specific humidity
//       trad,      &//radiative temperature [K]

    //emis,      &//averaged bulk surface emissivity
//       z0ma,      &//effective roughness [m]
//       zol,       &//dimensionless height (z/L) used in Monin-Obukhov theory
//       rib,       &//bulk Richardson number in surface layer
//       ustar,     &//u* in similarity theory [m/s]
//       qstar,     &//q* in similarity theory [kg/kg]
//       tstar,     &//t* in similarity theory [K]
//       u10m,      &//10m u-velocity
//       v10m,      &//10m v-velocity
//       f10m,      &//integral of profile function for momentum at 10m
//       fm,        &//integral of profile function for momentum
//       fh,        &//integral of profile function for heat
//       fq          //integral of profile function for moisture

// ------------------------ local variables ------------------------
    int
    niters,    //maximum number of iterations for surface temperature
    iter,      //iteration index
    nmozsgn,   //number of times moz changes sign
    idlak;       //index of lake, 1 = deep lake, 2 = shallow lake

    double  ax,    //
         bx,        //
         beta1,     //coefficient of conective velocity [-]
         degdT,     //d(eg)/dT
         displax,   //zero- displacement height [m]
         dqh,       //diff of humidity between ref. height and surface
         dth,       //diff of virtual temp. between ref. height and surface
         dthv,      //diff of vir. poten. temp. between ref. height and surface
         dzsur,     //
         eg,        //water vapor pressure at temperature T [pa]
         emg,       //ground emissivity (0.97 for snow,
         errore,    //lake temperature energy conservation error (w/m**2)
         hm,        //energy residual [W/m2]
         htvp,      //latent heat of vapor of water (or sublimation) [j/kg]
         obu,       //monin-obukhov length (m)
         obuold,    //monin-obukhov length of previous iteration
         qsatg,     //saturated humidity [kg/kg]
         qsatgdT,   //d(qsatg)/dT
         qseva,     //ground surface evaporation rate (mm h2o/s)
         qsdew,     //ground surface dew formation (mm h2o /s) [+]
         qsubl,     //sublimation rate from snow pack (mm h2o /s) [+]
         qfros,     //surface dew added to snow pack (mm h2o /s) [+]
         qmelt,     //snow melt [mm/s]
         ram,       //aerodynamical resistance [s/m]
         rah,       //thermal resistance [s/m]
         raw,       //moisture resistance [s/m]
         snowrate,  //rate of snowfall [mm/s]
         stftg3,    //
         temp1,     //relation for potential temperature profile
         temp2,     //relation for specific humidity profile
         temp12m,   // relation for temperature at 2m
         temp22m,   // relation for specific humidity at 2m
         tgbef,     //
         thm,       //intermediate variable (tm+0.0098*ht)
         th,        //potential temperature (kelvin)
         thv,       //virtual potential temperature (kelvin)
         thvstar,   //virtual potential temperature scaling parameter
         tksur,     //thermal conductivity of snow/soil (w/m/kelvin)

         um,        //wind speed including the stablity effect [m/s]
         ur,        //wind speed at reference height [m/s]
         visa,      // kinematic viscosity of dry air [m2/s]
         wc,        //convective velocity [m/s]
         wc2,       //wc*wc
         xt,        //
         xq,        //
         zeta,      //dimensionless height used in Monin-Obukhov theory
         zii,       //convective boundary height [m]
         zldis,     //reference height "minus" zero displacement heght [m]
         z0mg,      //roughness length over ground, momentum [m]
         z0hg,      //roughness length over ground, sensible heat [m]
         z0qg,      //roughness length over ground, latent heat [m]

         beta[2],   //fraction solar rad absorbed at surface: depends on lake type
         za[2],     //base of surface absorption layer (m): depends on lake type
         eta[2],    //light extinction coefficient (/m): depends on lake type
         p0,        //neutral value of turbulent prandtl number

         //a(nl_lake),    //"a" vector for tridiagonal matrix
         //b(nl_lake),    //"b" vector for tridiagonal matrix
         //c(nl_lake),    //"c" vector for tridiagonal matrix
         //r(nl_lake),    //"r" vector for tridiagonal solution
         //rhow(nl_lake), //density of water (kg/m**3)
         //phi(nl_lake),  //solar radiation absorbed by layer (w/m**2)
         //kme(nl_lake),  //molecular + eddy diffusion coefficient (m**2/s)

         cwat,      //specific heat capacity of water (j/m**3/kelvin)
         ws,        //surface friction velocity (m/s)
         ks,        //coefficient
         in,        //relative flux of solar radiation into layer
         out,       //relative flux of solar radiation out of layer
         ri,        //richardson number
         fin,       //heat flux into lake - flux out of lake (w/m**2)
         ocvts,     //(cwat*(tlak[n  ])*dzlak
         ncvts,     //(cwat*(tlak[n+1])*dzlak

         m1,        //intermediate variable for calculating r, a, b, c
         m2,        //intermediate variable for calculating r, a, b, c
         m3,        //intermediate variable for calculating r, a, b, c
         ke,        //eddy diffusion coefficient (m**2/s)
         km,        //molecular diffusion coefficient (m**2/s)
         zin,       //depth at top of layer (m)
         zout,      //depth at bottom of layer (m)
         drhodz,    //d [rhow] /dz (kg/m**4)
         n2,        //brunt-vaisala frequency (/s**2)
         num,       //used in calculating ri
         den,       //used in calculating ri
         tav,       //used in aver temp for convectively mixed layers
         nav,       //used in aver temp for convectively mixed layers
         phidum,    //temporary value of phi
         u2m;         //2 m wind speed (m/s)

    int i,j;       //do loop or array index

    //a(nl_lake),    //"a" vector for tridiagonal matrix
    //b(nl_lake),    //"b" vector for tridiagonal matrix
    //c(nl_lake),    //"c" vector for tridiagonal matrix
    //r(nl_lake),    //"r" vector for tridiagonal solution
    //rhow(nl_lake), //density of water (kg/m**3)
    //phi(nl_lake),  //solar radiation absorbed by layer (w/m**2)
    //kme(nl_lake),  //molecular + eddy diffusion coefficient (m**2/s)
    double *a    = new double [MAXSNL+MAXSOILL];
    double *b    = new double [MAXSNL+MAXSOILL];
    double *c    = new double [MAXSNL+MAXSOILL];
    double *r    = new double [MAXSNL+MAXSOILL];
    double *rhow = new double [MAXSNL+MAXSOILL];
    double *phi  = new double [MAXSNL+MAXSOILL];
    double *kme  = new double [MAXSNL+MAXSOILL];



// -----------------------------------------------------------------
//*[1] constants and model parameters
// -----------------------------------------------------------------
// constants for lake temperature model
    beta[0] = 0.4;
    beta[1] = 0.4;                              // (deep lake, shallow lake)
    za[0]   = 0.6;
    za[1]   = 0.5;
    eta[0]  = 0.1;
    eta[1]  = 0.5;
    p0      = 1.0;

// latent heat
    if (tm > tfrz)
        htvp = hvap;
    else
        htvp = hsub;
// surface emissivity
    emg = 0.97;

// deep lake or shallow
    idlak = 1;
    if (itypwat==5)
        idlak = 2;

    snowrate = 0;
    if (((prc+prl)>0.0)&&(tm<=(tfrz+tcrit)))
        snowrate = prc+prl;

// ----------------------------------------------------------------------
//*[2] surface temperature and fluxes
// ----------------------------------------------------------------------

    dzsur = dzlak[MAXSNL] + snowdp;


    //cout<<"tg1 :"<<tg<<endl;
    qsadv(tg,psrf,eg,degdT,qsatg,qsatgdT);
    //cout<<"tg2 :"<<tg<<endl;
    //cout<<"qsadv:"<<tg<<" "<<psrf<<" "<<eg<<" "<<degdT<<" "<<qsatg<<" "<<qsatgdT<<endl;

// potential temperatur at the reference height

    beta1=1.0;       // -  (in computing W_*)
    zii = 1000.0 ;   // m  (pbl height)

    thm = tm + 0.0098*ht;              // intermediate variable equivalent to
    // tm*(pgcm/psrf)**(rgas/cpair)
    th = tm*pow((100000./psrf),(rgas/cpair)); // potential T
    thv = th*(1.+0.61*qm) ;            // virtual potential T
    ur = max(0.1,std::sqrt(us*us+vs*vs));   // limit set to 1

// Initialization variables

    nmozsgn = 0;
    obuold = 0.0;

    dth   = thm-tg;
    dqh   = qm-qsatg;
    dthv  = dth*(1.+0.61*qm)+0.61*th*dqh;
    zldis = hu-0.0;

// aerodynamical roughness
// Kinematic viscosity of dry air (m2/s)- Andreas (1989) CRREL Rep. 89-11

    visa=1.326e-5*(1.+6.542e-3*tm + 8.301e-6*tm*tm - 4.84e-9*tm*tm*tm);
    ustar=0.06;
    wc=0.5;

    if (tg>=tfrz)         // unfrozen lake
        // loop to obtain initial and good ustar and zo
    {
        if (dthv>=0.0)
            um=max(ur,0.1);
        else
            um=std::sqrt(ur*ur+wc*wc);
        for (i=1; i<=5; i++)
        {
            z0mg=0.013*ustar*ustar/grav+0.11*visa/ustar;
            ustar=vonkar*um/std::log(zldis/z0mg);
        }
    }
    else                        // frozen lake
    {
        z0mg = 0.04;
    }
    z0qg = z0mg;
    z0hg = z0mg;

    //cout<<ur<<" "<<th<<" "<<thm<<" "<<thv<<" "<<dth<<endl;
    //cout<<dqh<<" "<<dthv<<" "<<zldis<<" "<<z0mg<<endl;
    moninobukini(ur,th,thm,thv,dth,dqh,dthv,zldis,z0mg,um,obu);

// ----------------------------------------------------------------------
    niters = 3;

    for (iter = 1; iter<= niters; iter++)      // begin stability iteration
    {
        tgbef = tg;
        if (tg>=tfrz)
            tksur = tkwat;
        else
            tksur = tkice;

        if (tg>=tfrz)      // unfrozen lake
        {
            z0mg=0.013*ustar*ustar/grav + 0.11*visa/ustar;
            xq=2.67*pow((ustar*z0mg/visa),0.25) - 2.57;
            xt= xq;
            z0qg=z0mg/std::exp(xq);
            z0hg=z0mg/std::exp(xt);
        }

// Evaluated stability-dependent variables using moz from prior iteration
        displax = 0.0;
        moninobuk(hu,   ht,   hq,   displax,z0mg,
                  z0hg, z0qg, obu,  um,     ustar,
                  temp1,temp2,temp12m,temp22m,f10m,
                  fm,   fh,   fq);

        obuold = obu;
//
// Get derivative of fluxes with repect to ground temperature

        ram    = 1.0/(ustar*ustar/um);
        rah    = 1.0/(temp1*ustar);
        raw    = 1.0/(temp2*ustar);

        stftg3 = emg*stefnc*tgbef*tgbef*tgbef;

        ax  = sabg + emg*frl + 3.*stftg3*tgbef
              + rhoair*cpair/rah*thm
              - htvp*rhoair/raw*(qsatg-qsatgdT*tgbef - qm)
              + tksur*tlak[MAXSNL]/dzsur;

        bx  = 4.*stftg3 + rhoair*cpair/rah
              + htvp*rhoair/raw*qsatgdT + tksur/dzsur;

        tg = ax/bx;
        if (fabs(tg-273)>50)
            tg = 280;

// surface fluxes of momentum, sensible and latent
// using ground temperatures from previous time step

        fseng = rhoair*cpair*(tg-thm)/rah;
        fevpg = rhoair*(qsatg+qsatgdT*(tg-tgbef)-qm)/raw;
        qsadv(tg,psrf,eg,degdT,qsatg,qsatgdT);
        dth=thm-tg;
        dqh=qm-qsatg;

        tstar = temp1*dth;
        qstar = temp2*dqh;

        thvstar=tstar+0.61*th*qstar;
        zeta=zldis*vonkar*grav*thvstar/(ustar*ustar*thv);
        if (zeta >= 0.0)     //stable
            zeta = min(2.,max(zeta,1.e-6));
        else                    //unstable
            zeta = max(-100.,min(zeta,-1.e-6));

        obu = zldis/zeta;

        if (zeta >= 0.)
            um = max(ur,0.1);
        else
        {
            wc = pow((-grav*ustar*thvstar*zii/thv),(1.0/3.0));
            wc2 = beta1*beta1*(wc*wc);
            um = std::sqrt(ur*ur+wc2);
        }

        if ((obuold*obu)< 0.0)
            nmozsgn = nmozsgn+1;
        if (nmozsgn >= 4)
            break;
    }
// ----------------------------------------------------------------------

// if snow on ground and tg > tfrz: reset tg = tfrz. reevaluate ground fluxes.
// energy inbalance used to melt snow. scv > 0.5 prevents spurious fluxes

    if ((scv > 0.5)&&(tg > tfrz))
    {
        tg = tfrz;
        fseng = rhoair*cpair*(tg-thm)/rah;
        fevpg = rhoair*(qsatg+qsatgdT*(tg-tgbef)-qm)/raw; //*qsatg and qsatgdT
        //*should be f(tgbef)
    }

// net longwave from ground to atmosphere
    olrg = (1.-emg)*frl + stftg3*(-3.*tgbef+4.*tg);

// additional variables for WRF model
    emis = emg;
    z0ma = z0mg;
    zol  = zeta;
    //rib  = min(5.,zol*pow(vonkar,3)*ustar*ustar/(temp1*um*um));

    // correct 27, December, 2006
    rib = min(5.,zol*ustar*ustar/(vonkar*temp1*um*um));
// radiative temperature
    trad = pow((olrg/stefnc),0.25);

// ground heat flux
    fgrnd = sabg + frl - olrg - fseng - htvp*fevpg;

    taux   = -rhoair*us/ram;
    tauy   = -rhoair*vs/ram;

    fsena  = fseng;
    fevpa  = fevpg;
    lfevpa = htvp*fevpg;

// 2 m height air temperature
    tref   = thm + temp1*dth * (1./temp12m - 1./temp1);
    qref   =  qm + temp2*dqh * (1./temp22m - 1./temp2);

// 10 m wind
    u10m = us/max(0.1,ur) * ustar/vonkar * f10m;
    v10m = vs/max(0.1,ur) * ustar/vonkar * f10m;

// energy residual for snow melting
    if ((scv > 0.0)&&(tg >= tfrz))
        hm = min( scv*hfus/dtime, max(fgrnd,0.) );
    else
        hm = 0.0;
    qmelt = hm/hfus;             // snow melt (mm/s)

// ----------------------------------------------------------------------
//*[3] lake layer temperature
// ----------------------------------------------------------------------

// lake density

    for (j=0; j<nl_lake; j++)
        rhow[MAXSNL+j] = 1000.*( 1.0 - 1.9549e-05*pow((fabs(tlak[MAXSNL+j]-277.)),1.68) );

// eddy diffusion +  molecular diffusion coefficient:
// eddy diffusion coefficient used for unfrozen deep lakes only

    cwat = cpliq*denh2o;
    km = tkwat/cwat;

    fin = beta[idlak-1]*sabg + frl - (olrg+fsena+lfevpa+hm);
    u2m = max(1.0,ustar/vonkar*std::log(2./z0mg));

    ws = 1.2e-03 * u2m;
    ks = 6.6 * std::sqrt( fabs(std::sin(dlat)) ) * (pow(u2m,(-1.84)));

    for (j = 0; j<nl_lake-1; j++)
    {
        drhodz = (rhow[MAXSNL+j+1]-rhow[MAXSNL+j]) / (zlak[MAXSNL+j+1]-zlak[MAXSNL+j]);
        n2 = -grav / rhow[MAXSNL+j] * drhodz;
        num = 40. * n2 * pow((vonkar*zlak[MAXSNL+j]),2);
        den = max( (ws*ws) * std::exp(-2.*ks*zlak[MAXSNL+j]), 1.e-10 );
        ri = ( -1. + std::sqrt( max(1.+num/den, 0.) ) ) / 20.0;
        if ((idlak == 1)&&(tg>tfrz))
            ke = vonkar*ws*zlak[MAXSNL+j]/p0 * std::exp(-ks*zlak[MAXSNL+j]) / (1.+37.*ri*ri);
        else
            ke = 0.0;
        kme[MAXSNL+j] = km + ke;
    }

    kme[nl_lake-1] = kme[nl_lake-2];

// heat source term: unfrozen lakes only

    for (j=0; j<nl_lake; j++)
    {
        zin  = zlak[MAXSNL+j] - 0.5*dzlak[MAXSNL+j];
        zout = zlak[MAXSNL+j] + 0.5*dzlak[MAXSNL+j];
        in  = std::exp( -eta[idlak-1]*max( zin-za[idlak-1],0.0 ) );
        out = std::exp( -eta[idlak-1]*max( zout-za[idlak-1],0.0 ) );
        //assumed solar absorption is only in the considered depth
        if (j == nl_lake-1)
            out = 0.0;
        if (tg > tfrz)
            phidum = (in-out) * sabg * (1.-beta[idlak-1]);
        else if (j ==0)
            phidum= sabg * (1.-beta[idlak-1]);
        else
            phidum = 0.0;
        phi[MAXSNL+j] = phidum;
    }

// sum cwat*tlak*dzlak for energy check

    ocvts = 0.0;
    for (j=0; j<nl_lake; j++)
        ocvts = ocvts + cwat*tlak[MAXSNL+j]*dzlak[MAXSNL+j];


// set up vector r and vectors a, b, c that define tridiagonal matrix

    j = 0;
    m2 = dzlak[MAXSNL+j]/kme[MAXSNL+j] + dzlak[MAXSNL+j+1]/kme[MAXSNL+j+1];
    m3 = dtime/dzlak[MAXSNL+j];
    r[MAXSNL+j] = tlak[MAXSNL+j] + (fin+phi[MAXSNL+j])*m3/cwat - (tlak[MAXSNL+j]-tlak[MAXSNL+j+1])*m3/m2;
    a[MAXSNL+j] = 0.0;
    b[MAXSNL+j] = 1. + m3/m2;
    c[MAXSNL+j] = -m3/m2;

    j = nl_lake-1;
    m1 = dzlak[MAXSNL+j-1]/kme[MAXSNL+j-1] + dzlak[MAXSNL+j]/kme[MAXSNL+j];
    m3 = dtime/dzlak[MAXSNL+j];
    r[MAXSNL+j] = tlak[MAXSNL+j] + phi[MAXSNL+j]*m3/cwat + (tlak[MAXSNL+j-1]-tlak[MAXSNL+j])*m3/m1;
    a[MAXSNL+j] = -m3/m1;
    b[MAXSNL+j] = 1. + m3/m1;
    c[MAXSNL+j] = 0.0;

    for (j=1; j<nl_lake-1; j++)
    {
        m1 = dzlak[MAXSNL+j-1]/kme[MAXSNL+j-1] + dzlak[MAXSNL+j]/kme[MAXSNL+j];
        m2 = dzlak[MAXSNL+j]/kme[MAXSNL+j] + dzlak[MAXSNL+j+1]/kme[MAXSNL+j+1];
        m3 = dtime/dzlak[MAXSNL+j];
        r[MAXSNL+j] = tlak[MAXSNL+j] + phi[MAXSNL+j]*m3/cwat
                      +(tlak[MAXSNL+j-1]-tlak[MAXSNL+j])*m3/m1 - (tlak[MAXSNL+j]-tlak[MAXSNL+j+1])*m3/m2;

        a[MAXSNL+j] = -m3/m1;
        b[MAXSNL+j] = 1. + m3/m1 + m3/m2;
        c[MAXSNL+j] = -m3/m2;
    }

// solve for tlak: a, b, c, r, u go from 1 to nsoi. tlak = 1 to npt

    double *a1 = new double [nl_lake];
    double *b1 = new double [nl_lake];
    double *c1 = new double [nl_lake];
    double *r1 = new double [nl_lake];
    double *tlak1 = new double [nl_lake];
    for (j=0; j<nl_lake; j++)
    {
        a1[j] = a[MAXSNL+j];
        b1[j] = b[MAXSNL+j];
        c1[j] = c[MAXSNL+j];
        r1[j] = r[MAXSNL+j];
    }


    tridia(nl_lake,a1 ,b1 ,c1 ,r1 ,tlak1) ;
    for (j=0; j<nl_lake; j++)
    {
        if ((tlak1[j]>220)&&(tlak1[j]<320))
            tlak[MAXSNL+j]  = tlak1[j];
    }

    delete []a1;
    delete []b1;
    delete []c1;
    delete []r1;
    delete []tlak1;

// convective mixing: make sure cwat*dzlak*ts is conserved. mixing
// is only allowed for unfrozen deep lakes. mix every 3 time steps

    if ((idlak == 1)&&(tg > tfrz))
    {
        for (j=1; j<=nl_lake-1; j++)
        {
            if (rhow[MAXSNL+j-1] > rhow[MAXSNL+j+1-1])
            {
                tav = 0.0;
                nav = 0.0;
                for (i = 1; i<=j+1; i++)
                {
                    tav = tav + tlak[MAXSNL+i-1]*dzlak[MAXSNL+i-1];
                    nav = nav + dzlak[MAXSNL+i-1];
                }
                tav = tav/nav;

                for (i=1; i<=j+1; i++)
                {
                    tlak[MAXSNL+i-1] = tav;
                    rhow[MAXSNL+i-1] = 1000.*( 1.0
                                               - 1.9549e-05*pow((fabs(tlak[MAXSNL+i-1]-277.)),1.68) );
                }
            }
        }
    }
// sum cwat*tlak*dzlak and total energy into lake for energy check

    ncvts = 0.0;
    for (j=0; j<nl_lake; j++)
    {
        ncvts = ncvts + cwat*tlak[MAXSNL+j]*dzlak[MAXSNL+j];
        fin = fin + phi[MAXSNL+j];
    }
    errore = (ncvts-ocvts) / dtime - fin;

// ----------------------------------------------------------------------
// [4] snow on the lake ice
// ----------------------------------------------------------------------

    qseva = 0.0;
    qsubl = 0.0;
    qfros = 0.0;
    qsdew = 0.0;

    if (fevpg >= 0.0)
// sublimation. do not allow for more sublimation than there is snow
// after melt. remaining surface evaporation used for infiltration
    {
        qsubl = min( fevpg, scv/dtime-qmelt );
        qseva = fevpg - qsubl;
    }
    else

        if (tg < (tfrz-0.1))
            qfros = fabs(fevpg);
        else
            qsdew = fabs(fevpg);


    // update snow pack
    scv = scv + (snowrate-qmelt-qsubl+qfros)*dtime;
    scv = max( scv, 0.0 );
    //to avoid the micro value
    scv = scv < 1e-7? 0 :scv;

    scv = min(scv,250.0);

// no snow if lake unfrozen
    if (tg > tfrz) scv = 0.0;

// snow height and fractional coverage
    snowdp = scv/250.0  ;     //assumed a constant snow bulk density = 250.

// null water mass
    for (j=0; j<nl_lake; j++)
    {
        wice[MAXSNL+j] = 0.0;
        wliq[MAXSNL+j] = 0.0;
    }
    delete [] a    ;
    delete [] b    ;
    delete [] c    ;
    delete [] r    ;
    delete [] rhow ;
    delete [] phi  ;
    delete [] kme  ;
}
void CommonLandModel::leafinterception (double dtime,double dewmx,double chil, double prc,
                                        double prl,  double tm,   double scv,  double sigf,
                                        double lai,  double sai,
                                        double &ldew,double &pg)
{

    //=======================================================================
    //
    // calculation of  interception and drainage of precipitation
    // the treatment are based on Sellers et al. (1996)
    //
    // modified by Yongjiu Dai, 08/31/2002
    //----------------------------------------------------------------------

    /*

    //-----------------------Arguments---------------------------------------

    real(r8), INTENT(in) :: dtime   // time step [second]
    real(r8), INTENT(in) :: dewmx   // maximum dew [mm]
    real(r8), INTENT(in) :: chil    // leaf angle distribution factor
    real(r8), INTENT(in) :: prc     // convective precipitation rate [mm/s]
    real(r8), INTENT(in) :: prl     // large-scale precipitation rate [mm/s]
    real(r8), INTENT(in) :: tm      // air temperature at reference height [K]
    real(r8), INTENT(in) :: scv     // snow mass (mm)
    real(r8), INTENT(in) :: sigf    // fraction of veg cover, excluding snow-covered veg [-]
    real(r8), INTENT(in) :: lai     // leaf area index [-]
    real(r8), INTENT(in) :: sai     // stem area index [-]

    real(r8), INTENT(inout) :: ldew // depth of water on foliage [mm]
    real(r8), INTENT(out) :: pg     // water onto ground including canopy runoff [kg/(m2 s)]

    */

    //-----------------------Local Variables---------------------------------

    double satcap;   // maximum allowed water on canopy [mm]
    double lsai;     // sum of leaf area index and stem area index [-]
    double chiv;     // leaf angle distribution factor
    double ppc;      // convective precipitation in time-step [mm]
    double ppl;      // large-scale precipitation in time-step [mm]
    double p0;       // precipitation in time-step [mm]
    double fpi;      // coefficient of interception
    double pinf;     // interception of precipitation in time step [mm]
    double tti;      // direct throughfall in time step [mm]
    double tex ;     // canopy drainage in time step [mm]
    double vegt ;    // sigf*lsai
    double xs;       // proportion of the grid area where the intercepted rainfall
    // plus the preexisting canopy water storage

    double ap, cp, bp, aa, bb, exrain, arg, thru, xsc, w;
    double pcoefs[2][2];


    //-----------------------End Variable List-------------------------------

    if (sigf>=0.001)
    {
        pcoefs[0][0] = 20.0;
        pcoefs[0][1] = 0.206e-8;
        pcoefs[1][0] = 0.0001;
        pcoefs[1][1] = 0.9999;
        bp = 20.0;
        lsai = lai + sai;
        vegt = sigf*lsai;
        satcap = dewmx*vegt;
        ppc = prc*dtime;
        ppl = prl*dtime;
        p0 =  ppc + ppl;
        w = ldew+p0;

        if ((scv>0.0)||(tm<tfrz))
            ppc = 0.0;

        ppl = p0 - ppc;
        xsc = max(0., (ldew-satcap));
        ldew = ldew - xsc;
        ap = pcoefs[1][0];
        cp = pcoefs[1][1];

        if (p0>1.e-8)
        {
            ap = ppc/p0 * pcoefs[0][0] + ppl/p0 * pcoefs[1][0];
            cp = ppc/p0 * pcoefs[0][1] + ppl/p0 * pcoefs[1][1];

            //----------------------------------------------------------------------
            //      proportional saturated area (xs) and leaf drainage(tex)
            //-----------------------------------------------------------------------

            chiv = chil;
            if ( fabs(chiv)<= 0.01 )
                chiv = 0.01;

            aa = 0.5 - 0.633 * chiv - 0.33 * chiv * chiv;
            bb = 0.877 * ( 1. - 2. * aa );
            exrain = aa + bb;
            fpi = ( 1.-std::exp(-exrain*lsai) ) * sigf;
            tti = p0 * ( 1.-fpi );
            xs = 1.;

            if ((p0*fpi)>1.e-9)
            {
                arg = (satcap-ldew)/(p0*fpi*ap) - cp/ap;
                if (arg>1.e-9)
                {
                    xs = -1./bp * std::log( arg );
                    xs = min( xs, 1. );
                    xs = max( xs, 0. );
                }
            }
            tex = p0 * fpi * ( ap/bp*(1.-std::exp(-bp*xs)) + cp*xs ) - ( satcap - ldew ) * xs;
            tex = max( tex, 0. );
            //cout<<p0<<" "<<satcap<<"  "<<ldew<<" "<<xs<<endl;
            //getchar();
            if ((tex+tti) > p0)
            {
                cout<<tex<<" "<<tti<<" "<<p0<<" "<<ldew<<endl;
                cout<<"tex + tti > p0 in interception code :"<<endl;

                //tti = p0-tex;
                //getchar();
                //exit(0);
            }

        }
        else
        {
            tti = 0.;
            tex = 0.;
        }

        //----------------------------------------------------------------------
        //   total throughfall (thru) and store augmentation
        //----------------------------------------------------------------------

        thru = tti + tex;
        pinf = p0 - thru;
        // do not know the reason here. (wlx,20090907)
        //if (pinf<1.0e-9) pinf = 0.0;
        ldew = ldew + pinf;

        pg = (xsc + thru) / dtime;

        w = w - ldew - pg*dtime;
        /*	if(pg>(ppl+ppc))
        	{
        		cout<<p0<<" "<<thru<<" "<<xsc<<endl;
        		cout<<pg<<"  "<<ppl<<" "<<ppc<<endl;
        		cout<<"error pg!"<<endl;
        		getchar();
        	}*/
        if (fabs(w)>1.e-6)
        {
            cout<<w<<" "<<ldew<<" "<<pg*dtime<<" "<<satcap<<endl;
            cout<<"something wrong in interception code :"<<endl;
            //ldew = ldew + w;
            //getchar();
            //exit(0);
        }
    }

    else
    {
        ldew=0.;
        pg = prc + prl;
    }
    //cout<<ldew<<"  "<<pg<<endl;
    //getchar();
}
void CommonLandModel::leaftemone(double dtime , double csoilc ,double dewmx  ,double htvp   ,double lai    ,
                                 double sai   , double displa ,double sqrtdi ,double z0m    ,double effcon ,
                                 double vmax25 ,double slti   ,double hlti   ,double shti   ,double hhti   ,
                                 double trda   ,double trdm   ,double trop   ,double gradm  ,double binter ,
                                 double extkn  ,double extkb  ,double extkd  ,double hu     ,double ht 	,
                                 double hq     ,double us     ,double vs     ,double thm    ,double th      ,
                                 double thv    ,double qm     ,double psrf   ,double rhoair ,double par    ,
                                 double sabv   ,double frl    ,double thermk ,double rstfac ,double po2m   ,
                                 double pco2m  ,double sigf   ,double etrc   ,double tg     ,double qg     ,
                                 double dqgdT  ,double emg    ,double & tl     ,double & ldew   ,double & taux   ,
                                 double & tauy   ,double & fseng  ,double & fevpg  ,double & cgrnd  ,double & cgrndl ,
                                 double & cgrnds ,double & tref   ,double & qref   ,double & rst    ,double & assim   ,
                                 double & respc  ,double & fsenl  ,double & fevpl  ,double & etr    ,double & dlrad  ,
                                 double &  ulrad  ,double & z0ma   ,double & zol    ,double & rib    ,double & ustar  ,
                                 double & qstar  ,double & tstar  ,double & f10m   ,double & fm     ,double & fh     ,double & fq)
{

//=======================================================================
// Original author : Yongjiu Dai, August 15, 2001
//
// Foliage energy conservation is given by foliage energy budget equation
//                      Rnet - Hf - LEf = 0
// The equation is solved by Newton-Raphson iteration, in which this iteration
// includes the calculation of the photosynthesis and stomatal resistance, and the
// integration of turbulent flux profiles. The sensible and latent heat
// transfer between foliage and atmosphere and ground is linked by the equations:
//                      Ha = Hf + Hg and Ea = Ef + Eg
//
//=======================================================================

//  use precision
//  use phycon_module, only : vonkar, grav, hvap, cpair, stefnc
//  implicit none
//
////-----------------------Arguments---------------------------------------
//
//  double , INTENT(in) :: &
//        dtime,    // time step [second]
//        csoilc,   // drag coefficient for soil under canopy [-]
//        dewmx,    // maximum dew
//        htvp         // latent heat of evaporation (/sublimation) [J/kg]
//
//// vegetation parameters
//  double , INTENT(in) :: &
//        lai,      // adjusted leaf area index for seasonal variation [-]
//        sai,      // stem area index  [-]
//        displa,   // displacement height [m]
//        sqrtdi,   // inverse sqrt of leaf dimension [m**-0.5]
//        z0m,      // roughness length, momentum [m]
//
//        effcon,   // quantum efficiency of RuBP regeneration (mol CO2 / mol quanta)
//        vmax25,   // maximum carboxylation rate at 25 C at canopy top
//                     // the range : 30.e-6 <-> 100.e-6 (mol co2 m-2 s-1)
//        shti,     // slope of high temperature inhibition function     (s1)
//        hhti,     // 1/2 point of high temperature inhibition function (s2)
//        slti,     // slope of low temperature inhibition function      (s3)
//        hlti,     // 1/2 point of low temperature inhibition function  (s4)
//        trda,     // temperature coefficient in gs-a model             (s5)
//        trdm,     // temperature coefficient in gs-a model             (s6)
//        trop,     // temperature coefficient in gs-a model         (273+25)
//        gradm,    // conductance-photosynthesis slope parameter
//        binter,   // conductance-photosynthesis intercept
//        extkn        // coefficient of leaf nitrogen allocation
//
//// input variables
//  double , INTENT(in) :: &
//        hu,       // observational height of wind [m]
//        ht,       // observational height of temperature [m]
//        hq,       // observational height of humidity [m]
//        us,       // wind component in eastward direction [m/s]
//        vs,       // wind component in northward direction [m/s]
//        thm,      // intermediate variable (tm+0.0098*ht)
//        th,       // potential temperature (kelvin)
//        thv,      // virtual potential temperature (kelvin)
//        qm,       // specific humidity at reference height [kg/kg]
//        psrf,     // pressure at reference height [pa]
//        rhoair,   // density air [kg/m**3]
//
//        par,      // par absorbed per unit lai [w/m**2]
//        sabv,     // solar radiation absorbed by vegetation [W/m2]
//        frl,      // atmospheric infrared (intwave) radiation [W/m2]
//
//        extkb,    // (k, g(mu)/mu) direct solar extinction coefficient
//        extkd,    // diffuse and scattered diffuse PAR extinction coefficient
//        thermk,   // canopy gap fraction for tir radiation
//        rstfac,   // factor of soil water stress to plant physiologocal processes
//
//        po2m,     // atmospheric partial pressure  o2 (pa)
//        pco2m,    // atmospheric partial pressure co2 (pa)
//
//        sigf,     // fraction of veg cover, excluding snow-covered veg [-]
//        etrc,     // maximum possible transpiration rate (mm/s)
//        tg,       // ground surface temperature [K]
//        qg,       // specific humidity at ground surface [kg/kg]
//        dqgdT,    // temperature derivative of "qg"
//        emg          // vegetation emissivity
//
//  double , INTENT(inout) :: &
//        tl,       // leaf temperature [K]
//        ldew,     // depth of water on foliage [mm]
//        taux,     // wind stress: E-W [kg/m/s**2]
//        tauy,     // wind stress: N-S [kg/m/s**2]
//        fseng,    // sensible heat flux from ground [W/m2]
//        fevpg,    // evaporation heat flux from ground [mm/s]
//        cgrnd,    // deriv. of soil energy flux wrt to soil temp [w/m2/k]
//        cgrndl,   // deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
//        cgrnds,   // deriv of soil latent heat flux wrt soil temp [w/m**2/k]
//        tref,     // 2 m height air temperature (kelvin)
//        qref         // 2 m height air specific humidity
//
//  double , INTENT(out) :: &
//        rst,      // stomatal resistance
//        assim,    // rate of assimilation
//        respc,    // rate of respiration
//        fsenl,    // sensible heat from leaves [W/m2]
//        fevpl,    // evaporation+transpiration from leaves [mm/s]
//        etr,      // transpiration rate [mm/s]
//        dlrad,    // downward intwave radiation blow the canopy [W/m2]
//        ulrad,    // upward intwave radiation above the canopy [W/m2]
//
//        z0ma,     // effective roughness [m]
//        zol,      // dimensionless height (z/L) used in Monin-Obukhov theory
//        rib,      // bulk Richardson number in surface layer
//        ustar,    // friction velocity [m/s]
//        tstar,    // temperature scaling parameter
//        qstar,    // moisture scaling parameter
//        f10m,     // integral of profile function for momentum at 10m
//        fm,       // integral of profile function for momentum
//        fh,       // integral of profile function for heat
//        fq           // integral of profile function for moisture

//-----------------------Local Variables---------------------------------

    double
    zldis,    // reference height "minus" zero displacement heght [m]
    zii,      // convective boundary layer height [m]
    z0mv,     // roughness length, momentum [m]
    z0hv,     // roughness length, sensible heat [m]
    z0qv,     // roughness length, latent heat [m]
    zeta,     // dimensionless height used in Monin-Obukhov theory
    beta,     // coefficient of conective velocity [-]
    wc,       // convective velocity [m/s]
    wc2,      // wc**2
    dth,      // diff of virtual temp. between ref. height and surface
    dthv,     // diff of vir. poten. temp. between ref. height and surface
    dqh,      // diff of humidity between ref. height and surface
    obu,      // monin-obukhov length (m)
    um,       // wind speed including the stablity effect [m/s]
    ur,       // wind speed at reference height [m/s]
    uaf,      // velocity of air within foliage [m/s]
    temp1,    // relation for potential temperature profile
    temp2,    // relation for specific humidity profile
    temp12m,  // relation for temperature at 2m
    temp22m,  // relation for specific humidity at 2m
    thvstar,  // virtual potential temperature scaling parameter
    taf,      // air temperature within canopy space [K]
    qaf,      // humidity of canopy air [kg/kg]
    eah,      // canopy air vapor pressure (pa)
//        pco2g,    // co2 pressure (pa) at ground surface (pa)
    pco2a,    // canopy air co2 pressure (pa)

    fdry,     // fraction of foliage that is green and dry [-]
    fwet,     // fraction of foliage covered by water [-]
    cf,       // heat transfer coefficient from leaves [-]
    rb,       // leaf boundary layer resistance [s/m]
    rbone,    // canopy bulk boundary layer resistance
    rd,       // aerodynamical resistance between ground and canopy air
    ram,      // aerodynamical resistance [s/m]
    rah,      // thermal resistance [s/m]
    raw,      // moisture resistance [s/m]
    clai,     // canopy heat capacity [Jm-2K-1]
    cah,      // heat conduactance for air [m/s]
    cgh,      // heat conduactance for ground [m/s]
    cfh,      // heat conduactance for leaf [m/s]
    caw,      // latent heat conduactance for air [m/s]
    cgw,      // latent heat conduactance for ground [m/s]
    cfw,      // latent heat conduactance for leaf [m/s]
    wtshi,    // sensible heat resistance for air, grd and leaf [-]
    wtsqi,    // latent heat resistance for air, grd and leaf [-]
    wta0,     // normalized heat conduactance for air [-]
    wtg0,     // normalized heat conduactance for ground [-]
    wtl0,     // normalized heat conductance for air and leaf [-]
    wtaq0,    // normalized latent heat conduactance for air [-]
    wtgq0,    // normalized heat conduactance for ground [-]
    wtlq0,    // normalized latent heat cond. for air and leaf [-]

    ei,       // vapor pressure on leaf surface [pa]
    deiDT,    // derivative of "ei" on "tl" [pa/K]
    qsatl,    // leaf specific humidity [kg/kg]
    qsatlDT,  // derivative of "qsatl" on "tlef"

    del,      // absolute change in leaf temp in current iteration [K]
    del2,     // change in leaf temperature in previous iteration [K]
    dele,     // change in heat fluxes from leaf [K]
    dele2,    // change in heat fluxes from leaf [K]
    det,      // maximum leaf temp. change in two consecutive iter [K]
    //del,      // maximum heat fluxes change in two consecutive iter [K]
    delmax,   // maximum change in  leaf temperature [K]
    dtmin,    // max limit for temperature convergence [K]
    dlemin,   // max limit for energy flux convergence [w/m2]

    obuold,   // monin-obukhov length from previous iteration
    tlbef,    // leaf temperature from previous iteration [K]
    //ecidif,   // excess energies [W/m2]
    err,      // balance error

    rs,       // leaf stomatal resistance [s/m]
    rsoil,    // soil respiration
    gah2o,    // conductance between canopy and atmosphere
    gdh2o,    // conductance between canopy and ground
    tprcor;   //  // tf*psur*100./1.013e5

    int it, itmax, itmin, nmozsgn;

    double  delta, fac;
    double  evplwet, evplwet_dtl, etr_dtl, elwmax, elwdif;
    double  irab, dirab_dtl, fsenl_dtl, fevpl_dtl;
    double  w, csoilcn, z0mg, cint[3];

    if (abs(tl-280)>40)
        tl = 280;

    // double , dimension(:), allocatable :: dtl // difference of tl between two iterative step

//-----------------------End Variable List-------------------------------

// initialization of errors and  iteration parameters
    it     = 1 ;   // counter for leaf temperature iteration
    del    = 0.0;  // change in leaf temperature from previous iteration
    dele   = 0.0;  // latent head flux from leaf for previous iteration
    um=0;
    rs=0;
// assign iteration parameters
    itmax  = 40 ;  // maximum number of iteration
    itmin  = 6 ;   // minimum number of iteration
    delmax = 3.0 ; // maximum change in  leaf temperature
    dtmin  = 0.01; // max limit for temperature convergence
    dlemin = 0.1;  // max limit for energy flux convergence

    double *dtl = new double [itmax+2];
    dtl[0] = 0.0;

//-----------------------------------------------------------------------
// scaling-up coefficients from leaf to canopy
//-----------------------------------------------------------------------

    cint[0] = (1.-std::exp(-extkn*lai))/extkn;
    cint[1] = (1.-std::exp(-extkd*lai))/extkd;
    cint[2] = lai;

//-----------------------------------------------------------------------
// get fraction of wet and dry canopy surface (fwet & fdry)
// initial saturated vapor pressure and humidity and their derivation
//-----------------------------------------------------------------------

    //clai = 4.2 * 1000. * 0.2
    clai = 0.0;

    dewfraction(sigf,lai,sai,dewmx,ldew,fwet,fdry);

    qsadv(tl,psrf,ei,deiDT,qsatl,qsatlDT);

//-----------------------------------------------------------------------
// initial for fluxes profile
//-----------------------------------------------------------------------

    nmozsgn = 0;    // number of times moz changes sign
    obuold = 0. ;   // monin-obukhov length from previous iteration
    zii = 1000. ;   // m  (pbl height)
    beta = 1.;      // -  (in computing W_*)
    det=0;

    z0mv = z0m;
    z0hv = z0m;
    z0qv = z0m;

    taf = 0.5 * (tg + thm);
    qaf = 0.5 * (qm + qg);

    pco2a = pco2m;
    tprcor = 44.6*273.16*psrf/1.013e5;
    rsoil = 0.;  //respiration (mol m-2 s-1)
//      rsoil = 1.22e-6*std::exp(308.56*(1./56.02-1./(tg-227.13)))
//      rsoil = rstfac * 0.23 * 15. * 2.**((tg-273.16-10.)/10.) * 1.e-6
//      rsoil = 5.22 * 1.e-6
    rsoil = 0.22 * 1.e-6;

    ur = max(0.1, std::sqrt(us*us+vs*vs));    // limit set to 0.1
    dth = thm - taf;
    dqh = qm - qaf;
    dthv = dth*(1.+0.61*qm) + 0.61*th*dqh;
    zldis = hu - displa;
    moninobukini(ur,th,thm,thv,dth,dqh,dthv,zldis,z0mv,um,obu);

// ======================================================================
//     BEGIN stability iteration
// ======================================================================
    while (it<=itmax)
    {
        tlbef = tl;
        del2 = del;
        dele2 = dele;

//-----------------------------------------------------------------------
// Aerodynamical resistances
//-----------------------------------------------------------------------
// Evaluate stability-dependent variables using moz from prior iteration
        moninobuk(hu,ht,hq,displa,z0mv,z0hv,
                  z0qv,obu,um,ustar,temp1,
                  temp2,temp12m,temp22m,f10m,
                  fm,fh,fq);
// Aerodynamic resistance
        ram = 1./(ustar*ustar/um);
        rah = 1./(temp1*ustar);
        raw = 1./(temp2*ustar);

// Bulk boundary layer resistance of leaves
        uaf = um*std::sqrt( 1./(ram*um) );
        cf = 0.01*sqrtdi/std::sqrt(uaf);
        rb = 0.5/(cf*uaf);

        // rd = 1./(csoilc*uaf)                 // BATS legacy
        // w = std::exp(-0.5*(lai+sai))              // Dickinson's modification :
        // csoilc = ( 1.-w + w*um/uaf)/rah      // "rah" here is the resistance over
        // rd = 1./(csoilc*uaf)                 // bare ground fraction

// modified by Xubin Zeng's suggestion at 08-07-2002
        z0mg = 0.01;
        w = std::exp(-(lai+sai));
        csoilcn = (vonkar/(0.13*pow((z0mg*uaf/1.5e-5),0.45)))*w + csoilc*(1.-w);
        rd = 1./(csoilcn*uaf);
//-----------------------------------------------------------------------
// stomatal resistances
//-----------------------------------------------------------------------

        if (lai>0.001)
        {
            rbone = rb / lai;
            eah = qaf * psrf / ( 0.622 + 0.378 * qaf );    // pa
            stomata(vmax25 ,effcon ,slti  ,hlti   ,shti ,
                    hhti  ,trda ,trdm   ,trop   ,gradm ,
                    binter ,thm  ,psrf  ,po2m ,pco2m  ,
                    pco2a  ,eah   ,ei     ,tl   ,par  ,
                    rbone ,raw  ,rstfac ,cint   ,assim ,
                    respc  ,rs );
        }
        else
        {
            rs = 2.e4;
            assim = 0.;
            respc = 0.;
        }
// above stomatal resistances are for the canopy, the stomatal rsistances
// and the "rb" in the following calculations are the average for single leaf. thus,
        rs = rs * lai;

//-----------------------------------------------------------------------
// dimensional and non-dimensional sensible and latent heat conductances
// for canopy and soil flux calculations.
//-----------------------------------------------------------------------

        delta = 0.0;
        if ( (fwet<0.99)&&((sigf*etrc)>1.e-12))
            if ((qsatl-qaf)>0.0)
                delta = 1.0;

        cah = sigf / rah;
        cgh = sigf / rd;
        cfh = sigf * (lai + sai) / rb;

        caw = sigf / raw;
        cgw = sigf / rd;
        cfw = sigf * ( (1.-delta*(1.-fwet))*(lai+sai)/rb + (1.-fwet)*delta*lai/(rb+rs) );
        wtshi = 1. / ( cah + cgh + cfh );
        wtsqi = 1. / ( caw + cgw + cfw );
        wta0 = cah * wtshi;
        wtg0 = cgh * wtshi;
        wtl0 = cfh * wtshi;
        wtaq0 = caw * wtsqi;
        wtgq0 = cgw * wtsqi;
        wtlq0 = cfw * wtsqi;

//-----------------------------------------------------------------------
// IR radiation, sensible and latent heat fluxes and their derivatives
//-----------------------------------------------------------------------
// the partial derivatives of areodynamical resistance are ignored
// which cannot be determined analtically
        fac = sigf * (1. - thermk);

// intwave absorption and their derivatives
        irab = (frl - 2. * stefnc * pow(tl,4) + emg*stefnc*pow(tg,4) ) * fac ;
        dirab_dtl = - 8. * stefnc * pow(tl,3)* fac;

// sensible heat fluxes and their derivatives
        fsenl = rhoair * cpair * cfh * ( (wta0 + wtg0)*tl - wta0*thm - wtg0*tg );
        fsenl_dtl = rhoair * cpair * cfh * (wta0 + wtg0);

// latent heat fluxes and their derivatives
        etr =   sigf * rhoair * (1.-fwet) * delta * lai / (rb + rs)*( (wtaq0 + wtgq0)*qsatl - wtaq0*qm - wtgq0*qg );
        etr_dtl =  sigf * rhoair * (1.-fwet) * delta * lai / (rb + rs)*(wtaq0 + wtgq0)*qsatlDT ;

        if (etr>=(sigf*etrc))
        {
            etr = sigf*etrc;
            etr_dtl = 0.;
        }
        evplwet = sigf * rhoair * (1.-delta*(1.-fwet)) * (lai+sai) / rb
                  * ( (wtaq0 + wtgq0)*qsatl - wtaq0*qm - wtgq0*qg );
        evplwet_dtl = sigf * rhoair * (1.-delta*(1.-fwet)) * (lai+sai) / rb
                      * (wtaq0 + wtgq0)*qsatlDT ;
        if (evplwet>=(ldew/dtime))
        {
            evplwet = ldew/dtime;
            evplwet_dtl = 0.;
        }
        fevpl = etr + evplwet;
        fevpl_dtl = etr_dtl + evplwet_dtl;

//-----------------------------------------------------------------------
// difference of temperatures by quasi-newton-raphson method for the non-linear system equations
//-----------------------------------------------------------------------

        dtl[it] = (sabv + irab - fsenl - hvap*fevpl)/((lai+sai)*clai/dtime - dirab_dtl + fsenl_dtl + hvap*fevpl_dtl);
        // check magnitude of change in leaf temperature limit to maximum allowed value
        if (it<itmax)
        {
            // put brakes on large temperature excursions
            if(abs(dtl[it])>delmax)
                dtl[it] = delmax*int(dtl[it]/abs(dtl[it]));
            //if (dtl[it]>delmax)
            //     dtl[it] = delmax;
            //if (dtl[it]<-delmax)
            //    dtl[it] = -delmax;
            if ((it>=2)&&((dtl[it-1]*dtl[it]<0.0)))
                dtl[it] = 0.5*(dtl[it-1] + dtl[it]);
        }
        tl = tlbef + dtl[it];

//-----------------------------------------------------------------------
// square roots differences of temperatures and fluxes for use as the condition of convergences
//-----------------------------------------------------------------------

        del  = abs(dtl[it]);
        dele = abs(dtl[it]) * std::sqrt( dirab_dtl*dirab_dtl + fsenl_dtl*fsenl_dtl + hvap*fevpl_dtl*fevpl_dtl );

//-----------------------------------------------------------------------
//  saturated vapor pressures and canopy air temperature, canopy air humidity
//-----------------------------------------------------------------------
// Recalculate leaf saturated vapor pressure (ei_)for updated leaf temperature
// and adjust specific humidity (qsatl_) proportionately
        qsadv(tl,psrf,ei,deiDT,qsatl,qsatlDT);

// update vegetation/ground surface temperature, canopy air temperature,
// canopy air humidity
        taf = wta0*thm + wtg0*tg + wtl0*tl;
        qaf = wtaq0*qm + wtgq0*qg + wtlq0*qsatl;

// update co2 partial pressure within canopy air
        gah2o = 1.0/raw * tprcor/thm ;                    // mol m-2 s-1
        gdh2o = 1.0/rd  * tprcor/thm ;                    // mol m-2 s-1
        pco2a = pco2m - 1.37*psrf/max(0.446,gah2o) * (assim - respc - rsoil);

//-----------------------------------------------------------------------
// Update monin-obukhov length and wind speed including the stability effect
//-----------------------------------------------------------------------

        dth = thm - taf;
        dqh = qm - qaf;
        tstar = temp1*dth;
        qstar = temp2*dqh;
        thvstar = tstar + 0.61*th*qstar;
        zeta = zldis*vonkar*grav*thvstar / (ustar*ustar*thv);

        if (zeta>=0.0)                            //stable
            zeta = min(2.,max(zeta,1.e-6));
        else                                             //unstable
            zeta = max(-100.,min(zeta,-1.e-6));
        obu = zldis/zeta;
        if (zeta>0.0)
            um = max(ur,0.1);
        else
        {
            wc = pow((-grav*ustar*thvstar*zii/thv),(1./3.));
            wc2 = beta*beta*(wc*wc);
            um = std::sqrt(ur*ur+wc2);
        }
        if ((obuold*obu)<0.0)
            nmozsgn = nmozsgn+1;
        if (nmozsgn>=4)
            obu = zldis/(-0.01);
        obuold = obu;

//-----------------------------------------------------------------------
// Test for convergence
//-----------------------------------------------------------------------
        it = it+1;
        if (it>itmin)
        {
            det = max(del,del2);
            del = max(dele,dele2);
            if ((det<dtmin)&&(del<dlemin))
            {
                break;
            }
        }

    }
// ======================================================================
//     END stability iteration
// ======================================================================

    z0ma = z0mv;
    zol = zeta;
    //rib = min(5.,zeta*pow(vonkar,3)*ustar*ustar/(temp1*um*um));
    // correct 27, December, 2006
    rib = min(5.,zol*ustar*ustar/(vonkar*temp1*um*um));


// canopy fluxes and total assimilation amd respiration

    if (lai>0.001)
        rst = rs / lai;
    else
    {
        rst = 2.0e+4;
        assim = 0.0;
        respc = 0.0;
    }
    respc = respc + rsoil;
// canopy fluxes and total assimilation amd respiration
    fsenl = fsenl + fsenl_dtl*dtl[it-1];

    etr     = etr     +     etr_dtl*dtl[it-1];
    evplwet = evplwet + evplwet_dtl*dtl[it-1];
    fevpl   = fevpl   +   fevpl_dtl*dtl[it-1];

    elwmax = ldew/dtime;
    elwdif = max(0., evplwet-elwmax);
    evplwet = min(evplwet, elwmax);

    fevpl = fevpl - elwdif;
    fsenl = fsenl + hvap*elwdif;

    // wind stresses
    taux = taux - sigf*rhoair*us/ram;
    tauy = tauy - sigf*rhoair*vs/ram;

//-----------------------------------------------------------------------
// fluxes from ground to canopy space
//-----------------------------------------------------------------------
    fseng = fseng + cpair*rhoair*cgh*(tg-taf);
    fevpg = fevpg +    rhoair*cgw*(qg-qaf);


    //cout<<cpair<<" "<<rhoair<<" "<<cgh<<" "<<cgw<<" "<<tg<<" "<<taf<<" "<<qg<<" "<<qaf<<endl;

    //fseng = max(fseng,0.0);
    //fevpg = max(fevpg,0.0);

//-----------------------------------------------------------------------
// downward (upward) intwave radiation below (above) the canopy
//-----------------------------------------------------------------------

    dlrad = sigf * thermk * frl
            + stefnc * fac * pow(tlbef,3) * (tlbef + 4.*dtl[it-1]);
    ulrad = stefnc * ( fac * pow(tlbef,3) * (tlbef + 4.*dtl[it-1])
                       + sigf*thermk*emg*pow(tg,4 ));

//-----------------------------------------------------------------------
// Derivative of soil energy flux with respect to soil temperature (cgrnd)
//-----------------------------------------------------------------------
    cgrnds = cgrnds + cpair*rhoair*cgh*(1.-wtg0);
    cgrndl = cgrndl + rhoair*cgw*(1.-wtgq0)*dqgdT;
    cgrnd  = cgrnds + cgrndl*htvp;
//-----------------------------------------------------------------------
// balance check
// (the computational error was created by the assumed 'dtl' in line 406-408)
//-----------------------------------------------------------------------

    err = sabv + irab + dirab_dtl*dtl[it-1] - fsenl - hvap*fevpl;

    //if(fabs(err)>0.2)
    //  printf("energy imbalance in leaftemone.F90 %ld, %f, %f,%f,%f,%f\n",it-1,err,sabv,irab,fsenl,hvap*fevpl);

//-----------------------------------------------------------------------
// Update dew accumulation (kg/m2)
//-----------------------------------------------------------------------

    //if (evplwet>0.0)
    ldew = max(0., ldew-evplwet*dtime);


//-----------------------------------------------------------------------
// 2 m height air temperature
//-----------------------------------------------------------------------

    tref = tref + sigf*(thm + temp1*dth * (1./temp12m - 1./temp1)) ;
    qref = qref + sigf*( qm + temp2*dqh * (1./temp22m - 1./temp2));

    delete [] dtl;
}
void CommonLandModel::leaftemtwo (double dtime ,    double csoilc  ,    double dewmx  ,    double htvp   ,    double lai    ,
                                  double sai   ,    double displa   ,   double sqrtdi  ,   double z0m    ,    double effcon ,
                                  double vmax25 ,   double slti    ,    double hlti     ,  double shti    ,   double hhti   ,
                                  double trda   ,   double trdm   ,     double trop    ,   double gradm    ,  double binter  ,
                                  double extkn  ,   double extkb  ,     double extkd  ,    double hu      ,   double ht      ,
                                  double hq      ,  double us     ,     double vs     ,    double thm    ,    double th      ,
                                  double thv      , double qm      ,    double psrf   ,    double rhoair ,    double parsun ,
                                  double parsha  ,  double sabvsun  ,   double sabvsha ,   double frl    ,    double fsun   ,
                                  double thermk ,   double rstfac  ,    double po2m     ,  double pco2m   ,   double sigf   ,
                                  double etrc   ,   double tg     ,     double qg      ,   double dqgdT    ,  double emg     ,
                                  double & tlsun  ,   double & tlsha  ,     double & ldew   ,    double & taux    ,   double & tauy     ,
                                  double & fseng   ,  double & fevpg  ,     double & cgrnd  ,    double & cgrndl ,    double & cgrnds  ,
                                  double & tref     , double & qref    ,    double & rst    ,    double & assim  ,    double & respc  ,
                                  double & fsenl   ,  double & fevpl    ,   double & etr     ,   double & dlrad  ,    double & ulrad  ,
                                  double & z0ma   ,   double & zol     ,    double & rib      ,  double & ustar   ,   double & qstar  ,
                                  double & tstar  ,   double & f10m   ,     double & fm      ,   double & fh       ,  double & fq  )
{

    //=======================================================================
    // Original author : Yongjiu Dai, August 15, 2001
    //
    // Foliage energy conservation is given by foliage energy budget equation
    //                   sunlit:  [Rnet - Hf - LEf] = 0
    //                   shaded:  [Rnet - Hf - LEf] = 0
    // The equation is solved by Newton-Raphson iteration, in which this iteration
    // includes the calculation of the photosynthesis and stomatal resistance, and the
    // integration of turbulent flux profiles. The sensible and latent heat
    // transfer between foliage and atmosphere and ground is linked by the equations:
    //                   Ha = [Hf]sunlit + [Hf]shaded + Hg
    //                   Ea = [Ef]sunlit + [Ef]shaded + Eg
    //
    //=======================================================================


    //-----------------------Arguments---------------------------------------

    //  real(r8), INTENT(in) :: &
    //        dtime,      &// time step [second]
    //        csoilc,     &// drag coefficient for soil under canopy [-]
    //        dewmx,      &// maximum dew
    //        htvp         // latent heat of evaporation (/sublimation) [J/kg]
    //
    //// vegetation parameters
    //  real(r8), INTENT(in) :: &
    //        lai,        &// adjusted leaf area index for seasonal variation [-]
    //        sai,        &// stem area index  [-]
    //        displa,     &// displacement height [m]
    //        sqrtdi,     &// inverse sqrt of leaf dimension [m**-0.5]
    //        z0m,        &// roughness length, momentum [m]
    //
    //        effcon,     &// quantum efficiency of RuBP regeneration (mol CO2 / mol quanta)
    //        vmax25,     &// maximum carboxylation rate at 25 C at canopy top
    //                     // the range : 30.e-6 <-> 100.e-6 (mol co2 m-2 s-1)
    //        shti,       &// slope of high temperature inhibition function     (s1)
    //        hhti,       &// 1/2 point of high temperature inhibition function (s2)
    //        slti,       &// slope of low temperature inhibition function      (s3)
    //        hlti,       &// 1/2 point of low temperature inhibition function  (s4)
    //        trda,       &// temperature coefficient in gs-a model             (s5)
    //        trdm,       &// temperature coefficient in gs-a model             (s6)
    //        trop,       &// temperature coefficient in gs-a model             (273+25)
    //        gradm,      &// conductance-photosynthesis slope parameter
    //        binter,     &// conductance-photosynthesis intercept
    //        extkn        // coefficient of leaf nitrogen allocation
    //
    //// input variables
    //  real(r8), INTENT(in) :: &
    //        hu,         &// observational height of wind [m]
    //        ht,         &// observational height of temperature [m]
    //        hq,         &// observational height of humidity [m]
    //        us,         &// wind component in eastward direction [m/s]
    //        vs,         &// wind component in northward direction [m/s]
    //        thm,        &// intermediate variable (tm+0.0098*ht)
    //        th,         &// potential temperature (kelvin)
    //        thv,        &// virtual potential temperature (kelvin)
    //        qm,         &// specific humidity at reference height [kg/kg]
    //        psrf,       &// pressure at reference height [pa]
    //        rhoair,     &// density air [kg/m**3]
    //
    //        parsun,     &// par absorbed per unit sunlit lai [w/m**2]
    //        parsha,     &// par absorbed per unit shaded lai [w/m**2]
    //        sabvsun,    &// solar radiation absorbed by vegetation [W/m2]
    //        sabvsha,    &// solar radiation absorbed by vegetation [W/m2]
    //        frl,        &// atmospheric infrared (intwave) radiation [W/m2]
    //
    //        fsun,       &// sunlit fraction of canopy
    //        extkb,      &// (k, g(mu)/mu) direct solar extinction coefficient
    //        extkd,      &// diffuse and scattered diffuse PAR extinction coefficient
    //        thermk,     &// canopy gap fraction for tir radiation
    //        rstfac,     &// factor of soil water stress to transpiration
    //
    //        po2m,       &// atmospheric partial pressure  o2 (pa)
    //        pco2m,      &// atmospheric partial pressure co2 (pa)
    //
    //        sigf,       &// fraction of veg cover, excluding snow-covered veg [-]
    //        etrc,       &// maximum possible transpiration rate (mm/s)
    //        tg,         &// ground surface temperature [K]
    //        qg,         &// specific humidity at ground surface [kg/kg]
    //        dqgdT,      &// temperature derivative of "qg"
    //        emg          // vegetation emissivity
    //
    //  real(r8), INTENT(inout) :: &
    //        tlsun,      &// sunlit leaf temperature [K]
    //        tlsha,      &// shaded leaf temperature [K]
    //        ldew,       &// depth of water on foliage [mm]
    //        taux,       &// wind stress: E-W [kg/m/s**2]
    //        tauy,       &// wind stress: N-S [kg/m/s**2]
    //        fseng,      &// sensible heat flux from ground [W/m2]
    //        fevpg,      &// evaporation heat flux from ground [mm/s]
    //        cgrnd,      &// deriv. of soil energy flux wrt to soil temp [w/m2/k]
    //        cgrndl,     &// deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
    //        cgrnds,     &// deriv of soil latent heat flux wrt soil temp [w/m**2/k]
    //        tref,       &// 2 m height air temperature (kelvin)
    //        qref         // 2 m height air specific humidity
    //
    //  real(r8), INTENT(out) :: &
    //        rst,        &// total stomatal resistance (s m-1)
    //        assim,      &// total assimilation (mol m-2 s-1)
    //        respc,      &// total respiration (mol m-2 s-1)
    //        fsenl,      &// sensible heat from leaves [W/m2]
    //        fevpl,      &// evaporation+transpiration from leaves [mm/s]
    //        etr,        &// transpiration rate [mm/s]
    //        dlrad,      &// downward intwave radiation blow the canopy [W/m2]
    //        ulrad,      &// upward intwave radiation above the canopy [W/m2]
    //
    //        z0ma,       &// effective roughness [m]
    //        zol,        &// dimensionless height (z/L) used in Monin-Obukhov theory
    //        rib,        &// bulk Richardson number in surface layer
    //        ustar,      &// friction velocity [m/s]
    //        tstar,      &// temperature scaling parameter
    //        qstar,      &// moisture scaling parameter
    //        f10m,       &// integral of profile function for momentum at 10m
    //        fm,         &// integral of profile function for momentum
    //        fh,         &// integral of profile function for heat
    //        fq           // integral of profile function for moisture

    //-----------------------Local Variables---------------------------------

    double
    zldis,      // reference height "minus" zero displacement heght [m]
    zii,        // convective boundary layer height [m]
    z0mv,       // roughness length, momentum [m]
    z0hv,       // roughness length, sensible heat [m]
    z0qv,       // roughness length, latent heat [m]
    zeta,       // dimensionless height used in Monin-Obukhov theory
    beta,       // coefficient of conective velocity [-]
    wc,         // convective velocity [m/s]
    wc2,        // wc**2
    dth,        // diff of virtual temp. between ref. height and surface
    dthv,       // diff of vir. poten. temp. between ref. height and surface
    dqh,        // diff of humidity between ref. height and surface
    obu,        // monin-obukhov length (m)
    um,         // wind speed including the stablity effect [m/s]
    ur,         // wind speed at reference height [m/s]
    uaf,        // velocity of air within foliage [m/s]
    temp1,      // relation for potential temperature profile
    temp2,      // relation for specific humidity profile
    temp12m,    // relation for temperature at 2m
    temp22m,    // relation for specific humidity at 2m
    thvstar,    // virtual potential temperature scaling parameter
    taf,        // air temperature within canopy space [K]
    qaf,        // humidity of canopy air [kg/kg]
    eah,        // canopy air vapor pressure (pa)
    pco2a,      // canopy air co2 pressure (pa)

    fdry,       // fraction of foliage that is green and dry [-]
    fwet,       // fraction of foliage covered by water [-]
    cf,         // heat transfer coefficient from leaves [-]
    rb,         // leaf boundary layer resistance [s/m]
    rbsun,      // bulk boundary layer resistance of sunlit fraction of canopy
    rbsha,      // bulk boundary layer resistance of shaded fraction of canopy
    rd,         // aerodynamical resistance between ground and canopy air
    ram,        // aerodynamical resistance [s/m]
    rah,        // thermal resistance [s/m]
    raw,        // moisture resistance [s/m]
    cah,        // heat conduactance for air [m/s]
    cgh,        // heat conduactance for ground [m/s]
    cfsunh,     // heat conduactance for sunlit leaf [m/s]
    cfshah,     // heat conduactance for shaded leaf [m/s]
    caw,        // latent heat conduactance for air [m/s]
    cgw,        // latent heat conduactance for ground [m/s]
    cfsunw,     // latent heat conduactance for sunlit leaf [m/s]
    cfshaw,     // latent heat conduactance for shaded leaf [m/s]
    wtshi,      // sensible heat resistance for air, grd and leaf [-]
    wtsqi,      // latent heat resistance for air, grd and leaf [-]
    wta0,       // normalized heat conduactance for air [-]
    wtg0,       // normalized heat conduactance for ground [-]
    wtlsun0,    // normalized heat conductance for air and sunlit leaf [-]
    wtlsha0,    // normalized heat conductance for air and shaded leaf [-]
    wtaq0,      // normalized latent heat conduactance for air [-]
    wtgq0,      // normalized heat conduactance for ground [-]
    wtlsunq0,   // normalized latent heat cond. for air and sunlit leaf [-]
    wtlshaq0,   // normalized latent heat cond. for air and shaded leaf [-]

    del,        // absolute change in leaf temp in current iteration [K]
    del2,       // change in leaf temperature in previous iteration [K]
    dele,       // change in heat fluxes from leaf [K]
    dele2,      // change in heat fluxes from leaf [K]
    det,        // maximum leaf temp. change in two consecutive iter [K]
    //        del,        // maximum heat fluxes change in two consecutive iter [K]
    delmax,     // maximum change in  leaf temperature
    dtmin,      // max limit for temperature convergence [K]
    dlemin,     // max limit for energy flux convergence [w/m2]

    obuold,     // monin-obukhov length from previous iteration
    tlsunbef,   // sunlit leaf temperature from previous iteration [K]
    tlshabef,   // shaded leaf temperature from previous iteration [K]
//		ecidif,     // excess energies [W/m2]
    err,        // balance error

    fsha,       // shaded fraction of canopy
    laisun,     // sunlit leaf area index, one-sided
    laisha,     // shaded leaf area index, one-sided
    rssun,      // sunlit leaf stomatal resistance [s/m]
    rssha,      // shaded leaf stomatal resistance [s/m]
    assimsun,   // sunlit leaf assimilation rate [umol co2 /m**2/ s] [+]
    assimsha,   // shaded leaf assimilation rate [umol co2 /m**2/ s] [+]
    respcsun,   // sunlit leaf respiration rate [umol co2 /m**2/ s] [+]
    respcsha,   // shaded leaf respiration rate [umol co2 /m**2/ s] [+]
    rsoil,      // soil respiration
    gah2o,      //
    tprcor ;      //

    int  it, itmax, itmin, nmozsgn;
    double   delta1, delta2, fac, fac1, fac2 ;
    double   cintsun[3], cintsha[3];
    double   etrsun, etrsha, evplwetsun, evplwetsha, evplwet;
    double   etrsun_dtlsun, etrsha_dtlsun, evplwetsun_dtlsun, evplwetsha_dtlsun;
    double   etrsun_dtlsha, etrsha_dtlsha, evplwetsun_dtlsha, evplwetsha_dtlsha;

    double   irabsun, senlsun, evplsun, ftsun;
    double   irabsha, senlsha, evplsha, ftsha;
    double   dirabsun_dtlsun, senlsun_dtlsun, dirabsun_dtlsha, senlsun_dtlsha ;
    double   dirabsha_dtlsun, senlsha_dtlsun, dirabsha_dtlsha, senlsha_dtlsha;
    double   evplsun_dtlsun, evplsun_dtlsha, dftsunDTlsun, dftsunDTlsha;
    double   evplsha_dtlsun, evplsha_dtlsha, dftshaDTlsun, dftshaDTlsha ;

    double   eisun, eisha, deisunDT, deishaDT;
    double   qsatlsun, qsatlsha, qsatlsunDT, qsatlshaDT;
    double   z0mg, w, csoilcn;
    double   clai, elwmax, elwdif;

    // real(r8), dimension(:), allocatable :: &
    //      dtlsun,     &// difference of tlsun between two iterative step
    //      dtlsha       // difference of tlsha between two iterative step

    //-----------------------------------------------------------------------
    // initialization of errors and  iteration parameters
    //-----------------------------------------------------------------------

    it     = 1 ;   // counter for leaf temperature iteration
    del    = 0.0 ; // change in leaf temperature from previous iteration
    dele   = 0.0;  // latent head flux from leaf for previous iteration
    det=0;

    // assign iteration parameters
    itmax  = 40 ;  // maximum number of iteration
    itmin  = 6  ;  // minimum number of iteration
    delmax = 1.0;  // maximum change in  leaf temperature
    dtmin  = 0.01; // max limit for temperature convergence
    dlemin = 0.1 ; // max limit for energy flux convergence

    double *dtlsun = new double [itmax+2];
    double *dtlsha = new double [itmax+2];

    dtlsun[0] = 0.0;
    dtlsha[0] = 0.0;

    //-----------------------------------------------------------------------
    // leaf area index
    //-----------------------------------------------------------------------
    // partion visible canopy absorption to sunlit and shaded fractions
    // to get average absorbed par for sunlit and shaded leaves
    fsha = 1. - fsun;

    laisun = lai*fsun;
    laisha = lai*fsha;

    // scaling-up coefficients from leaf to canopy
    cintsun[0] = (1.-std::exp(-(extkn+extkb)*lai))/(extkn+extkb);
    cintsun[1] = (1.-std::exp(-(extkb+extkd)*lai))/(extkb+extkd);
    cintsun[2] = (1.-std::exp(-extkb*lai))/extkb;

    cintsha[0] = (1.-std::exp(-extkn*lai))/extkn - cintsun[0];
    cintsha[1] = (1.-std::exp(-extkd*lai))/extkd - cintsun[1];
    cintsha[2] = lai - cintsun[2];

    //-----------------------------------------------------------------------
    // get fraction of wet and dry canopy surface (fwet & fdry)
    // initial saturated vapor pressure and humidity and their derivation
    //-----------------------------------------------------------------------
    //暂时取消边界控制(wlx,20090912)
    //if (abs(tlsun-280)>40)  tlsun = 280;
    //if (abs(tlsha-280)>40)  tlsha = 280;

    //clai = 4.2 * 1000. * 0.2
    clai = 0.0;
    dewfraction (sigf,lai,sai,dewmx,ldew,fwet,fdry);
    qsadv(tlsun,psrf,eisun,deisunDT,qsatlsun,qsatlsunDT);
    qsadv(tlsha,psrf,eisha,deishaDT,qsatlsha,qsatlshaDT);
    //cout<<eisun<<" "<<deisunDT<<" "<<qsatlsun<<" "<<qsatlsunDT<<endl;
    //cout<<eisha<<" "<<deishaDT<<" "<<qsatlsha<<" "<<qsatlshaDT<<endl;

    //-----------------------------------------------------------------------
    // initial for fluxes profile
    //-----------------------------------------------------------------------
    nmozsgn = 0;    // number of times moz changes sign
    obuold = 0.  ;  // monin-obukhov length from previous iteration
    zii = 1000. ;   // m  (pbl height)
    beta = 1.  ;    // -  (in computing W_*)

    z0mv = z0m;
    z0hv = z0m;
    z0qv = z0m;

    taf = 0.5 * (tg + thm);
    qaf = 0.5 * (qm + qg);

    pco2a = pco2m;
    tprcor = 44.6*273.16*psrf/1.013e+5;
    rsoil = 0.; //soil respiration (mol m-2 s-1)
    //      rsoil = 1.22e-6*std::exp(308.56*(1./56.02-1./(tg-227.13)))
    //      rsoil = rstfac * 0.23 * 15. * 2.**((tg-273.16-10.)/10.) * 1.e-6
    //      rsoil = 5.22 * 1.e-6
    rsoil = 0.22 * 1.e-6;

    ur = max(0.1, std::sqrt(us*us+vs*vs));    // limit set to 0.1
    dth = thm - taf;
    dqh = qm - qaf;
    dthv = dth*(1.+0.61*qm) + 0.61*th*dqh;
    zldis = hu - displa;

    moninobukini(ur,th,thm,thv,dth,dqh,dthv,zldis,z0mv,um,obu);
    //cout<<ur<<" "<<th<<" "<<thm<<" "<<thv<<" "<<dth<<" "<<dqh<<endl;
    //cout<<dthv<<" "<<zldis<<" "<<z0mv<<" "<<um<<" "<<obu<<endl;

    // ======================================================================
    //     BEGIN stability iteration
    // ======================================================================
    while (it<=itmax)
    {
        tlsunbef = tlsun;
        tlshabef = tlsha;

        del2 = del;
        dele2 = dele;
        //-----------------------------------------------------------------------
        // Aerodynamical resistances
        //-----------------------------------------------------------------------
        // Evaluate stability-dependent variables using moz from prior iteration
        moninobuk(hu,ht,hq,displa,z0mv,
                  z0hv,z0qv,obu,um,ustar,
                  temp1,temp2,temp12m,temp22m,f10m ,
                  fm,fh,fq);
        //cout<<hu<<" "<<ht<<" "<<hq<<" "<<displa<<" "<<z0mv<<endl;
        //cout<<z0hv<<" "<<z0qv<<" "<<obu<<" "<<um<<" "<<ustar<<endl;
        //cout<<temp1<<" "<<temp2<<" "<<temp12m<<" "<<temp22m<<endl;
        //cout<<f10m<<" "<<fm<<" "<<fh<<" "<<fq<<endl;

        // Aerodynamic resistance
        ram = 1./(ustar*ustar/um) ;
        rah = 1./(temp1*ustar);
        raw = 1./(temp2*ustar);

        // Bulk boundary layer resistance of leaves
        uaf = um*std::sqrt( 1./(ram*um) );
        cf = 0.01*sqrtdi/std::sqrt(uaf);
        rb = 0.5/(cf*uaf);

        //cout<<um<<"  "<<ustar<<"  "<<ram<<endl;
        //cout<<uaf<<"  "<<cf<<"  "<<rb<<endl;
        //<->    rd = 1./(csoilc*uaf)          // legacy from BATS
        // modified by Xubin Zeng's suggestion at 08-07-2002
        z0mg = 0.01;
        w = std::exp(-(lai+sai));
        csoilcn = (vonkar/(0.13*pow((z0mg*uaf/1.5e-5),0.45)))*w + csoilc*(1.-w);
        rd = 1./(csoilcn*uaf);
        //cout<<w<<"  "<<csoilcn<<"  "<<uaf<<"  "<<rd<<endl;
        //-----------------------------------------------------------------------
        // stomatal resistances for sunlit and shaded fractions of canopy.
        // should do each iteration to account for differences in eah, leaf temperatures.
        //-----------------------------------------------------------------------
        if (lai>0.001)
        {
            rbsun = rb / laisun;
            rbsha = rb / laisha;

            eah = qaf * psrf / ( 0.622 + 0.378 * qaf );    // pa

            // Sunlit leaves
            stomata(vmax25 ,effcon  ,slti     ,hlti     ,shti   ,
                    hhti  , trda ,   trdm   ,  trop    , gradm  ,
                    binter ,thm  ,   psrf  ,   po2m ,    pco2m  ,
                    pco2a   ,eah    ,eisun    ,tlsun  ,   parsun ,
                    rbsun ,  raw  ,  rstfac ,cintsun ,    assimsun ,
                    respcsun ,rssun);
            // Shaded leaves
            stomata( vmax25 ,effcon  ,slti     ,hlti     ,shti ,
                     hhti  , trda ,   trdm   ,  trop    , gradm,
                     binter ,thm  ,   psrf  ,   po2m ,    pco2m,
                     pco2a   ,eah    ,eisha   , tlsha  ,  parsha ,
                     rbsha ,  raw  ,  rstfac ,  cintsha ,  assimsha ,
                     respcsha ,rssha );

        }

        else
        {
            rssun = 2.0e+4 ;
            rssha = 2.0e+4;
            assimsun = 0. ;
            assimsha = 0.;
            respcsun = 0. ;
            respcsha = 0.;
        }

        // above stomatal resistances are for the canopy, the stomatal resistance
        // the "rb" in the following calculations are the averages for single leaf. thus,
        rssun = rssun * laisun;
        rssha = rssha * laisha;
        //-----------------------------------------------------------------------
        // dimensional and non-dimensional sensible and latent heat conductances
        // for canopy and soil flux calculations.
        //-----------------------------------------------------------------------
        delta1 = 0.0;
        delta2 = 0.0;
        if ( (fwet<0.99)&&((sigf*etrc)>1.e-12))
        {
            if ((qsatlsun-qaf)>0.0)
                delta1 = 1.0;
            if ((qsatlsha-qaf)> 0.0)
                delta2 = 1.0;
        }

        cah = sigf / rah;
        cgh = sigf / rd;
        cfsunh = sigf * laisun / rb;
        cfshah = sigf * (laisha + sai) / rb;

        caw = sigf / raw;
        cgw = sigf / rd;
        cfsunw = sigf * ( (1.-delta1*(1.-fwet)) * laisun / rb
                          + (1. - fwet) * delta1 * laisun / (rb + rssun) );
        cfshaw = sigf * ( (1.-delta2*(1.-fwet)) * (laisha + sai) / rb
                          + (1. - fwet) * delta2 * laisha / (rb + rssha) );

        wtshi = 1. / ( cah + cgh + cfsunh + cfshah );
        wtsqi = 1. / ( caw + cgw + cfsunw + cfshaw );

        wta0 = cah * wtshi;
        wtg0 = cgh * wtshi;
        wtlsun0 = cfsunh * wtshi;
        wtlsha0 = cfshah * wtshi;

        wtaq0 = caw * wtsqi;
        wtgq0 = cgw * wtsqi;
        wtlsunq0 = cfsunw * wtsqi;
        wtlshaq0 = cfshaw * wtsqi;
        //-----------------------------------------------------------------------
        // IR radiation, sensible and latent heat fluxes and their derivatives
        //-----------------------------------------------------------------------
        // the partial derivatives of areodynamical resistance are ignored
        // which cannot be determined analtically
        fac = sigf * (1. - thermk);
        fac1 = fac * fsun;
        fac2 = fac * fsha;

        // intwave absorption and their derivatives
        irabsun = (frl - 2. * stefnc * pow(tlsun,4) + emg*stefnc*pow(tg,4) ) * fac1 ;
        dirabsun_dtlsun = - 8.* stefnc * pow(tlsun,3) * fac1;
        dirabsun_dtlsha = 0.;

        irabsha = (frl - 2. * stefnc * pow(tlsha,4) + emg*stefnc*pow(tg,4) ) * fac2;
        dirabsha_dtlsha = - 8.* stefnc * pow(tlsha,3)                    * fac2 ;
        dirabsha_dtlsun = 0.;

        // sensible heat fluxes and their derivatives
        senlsun = rhoair * cpair * cfsunh
                  * ( (wta0 + wtg0 + wtlsha0)*tlsun
                      - wta0*thm - wtg0*tg - wtlsha0*tlsha );
        senlsun_dtlsun = rhoair * cpair * cfsunh * (wta0 + wtg0 + wtlsha0);
        senlsun_dtlsha = rhoair * cpair * cfsunh * ( - wtlsha0 );


        senlsha = rhoair * cpair * cfshah
                  * ( (wta0 + wtg0 + wtlsun0)*tlsha
                      - wta0*thm - wtg0*tg - wtlsun0*tlsun );
        senlsha_dtlsun = rhoair * cpair * cfshah * ( - wtlsun0 );
        senlsha_dtlsha = rhoair * cpair * cfshah * ( wta0 + wtg0 + wtlsun0 );

        //cout<<it<<"  "<<senlsun<<"  "<<senlsun_dtlsun<<"  "<<senlsun_dtlsha<<endl;
        //cout<<it<<"  "<<senlsha<<"  "<<senlsha_dtlsun<<"  "<<senlsha_dtlsha<<endl;

        // latent heat fluxes and their derivatives
        etrsun = sigf * rhoair * (1.-fwet) * delta1 * laisun / (rb + rssun)
                 * ( (wtaq0 + wtgq0 + wtlshaq0)*qsatlsun
                     - wtaq0*qm - wtgq0*qg - wtlshaq0*qsatlsha );
        etrsun_dtlsun = sigf * rhoair * (1.-fwet) * delta1 * laisun / (rb + rssun)
                        * (wtaq0 + wtgq0 + wtlshaq0)*qsatlsunDT ;
        etrsun_dtlsha = sigf * rhoair * (1.-fwet) * delta1 * laisun / (rb + rssun)
                        * ( - wtlshaq0*qsatlshaDT );

        if (etrsun>=(sigf*etrc*laisun/lai))
        {
            etrsun =  sigf*etrc*laisun/lai;
            etrsun_dtlsun = 0.;
            etrsun_dtlsha = 0.;
        }

        evplwetsun = sigf * rhoair * (1.-delta1*(1.-fwet)) * laisun / rb
                     * ( (wtaq0 + wtgq0 + wtlshaq0)*qsatlsun
                         - wtaq0*qm - wtgq0*qg - wtlshaq0*qsatlsha );
        evplwetsun_dtlsun = sigf * rhoair * (1.-delta1*(1.-fwet)) * laisun / rb
                            * (wtaq0 + wtgq0 + wtlshaq0)*qsatlsunDT ;
        evplwetsun_dtlsha = sigf * rhoair * (1.-delta1*(1.-fwet)) * laisun / rb
                            * ( - wtlshaq0*qsatlshaDT );

        if (evplwetsun>=(ldew/dtime*laisun/lai))
        {
            evplwetsun = ldew/dtime*laisun/lai;
            evplwetsun_dtlsun = 0.;
            evplwetsun_dtlsha = 0.;
        }

        evplsun = etrsun + evplwetsun;
        evplsun_dtlsun = etrsun_dtlsun + evplwetsun_dtlsun;
        evplsun_dtlsha = etrsun_dtlsha + evplwetsun_dtlsha;

        // ---------------

        etrsha = sigf * rhoair * (1.-fwet) * delta2 * laisha / (rb + rssha)
                 * ( (wtaq0 + wtgq0 + wtlsunq0)*qsatlsha
                     - wtaq0*qm - wtgq0*qg - wtlsunq0*qsatlsun );
        etrsha_dtlsun = sigf * rhoair * (1.-fwet) * delta2 * laisha / (rb + rssha)
                        * ( - wtlsunq0*qsatlsunDT );
        etrsha_dtlsha = sigf * rhoair * (1.-fwet) * delta2 * laisha / (rb + rssha)
                        * ( (wtaq0 + wtgq0 + wtlsunq0)*qsatlshaDT );

        if (etrsha>=(sigf*etrc*laisha/lai))
        {
            etrsha = sigf*etrc*laisha/lai;
            etrsha_dtlsun = 0.;
            etrsha_dtlsha = 0.;
        }
        evplwetsha = sigf * rhoair * (1.-delta2*(1.-fwet)) * (laisha+sai) / (rb)
                     * ( (wtaq0 + wtgq0 + wtlsunq0)*qsatlsha
                         - wtaq0*qm - wtgq0*qg - wtlsunq0*qsatlsun );
        evplwetsha_dtlsun = sigf * rhoair * (1.-delta2*(1.-fwet)) * (laisha+sai) / rb
                            * ( - wtlsunq0*qsatlsunDT );
        evplwetsha_dtlsha = sigf * rhoair * (1.-delta2*(1.-fwet)) * (laisha+sai) / rb
                            * ( (wtaq0 + wtgq0 + wtlsunq0)*qsatlshaDT );

        if (evplwetsha>=(ldew/dtime*(laisha+sai)/(lai+sai)))
        {
            evplwetsha = ldew/dtime*(laisha+sai)/(lai+sai);
            evplwetsha_dtlsun = 0.;
            evplwetsha_dtlsha = 0;
        }
        evplsha = etrsha + evplwetsha;

        evplsha_dtlsun = etrsha_dtlsun + evplwetsha_dtlsun;
        evplsha_dtlsha = etrsha_dtlsha + evplwetsha_dtlsha;

        // functions and their derivatives with respect to temperatures
        ftsun = sabvsun + irabsun - senlsun - hvap*evplsun ;
        ftsha = sabvsha + irabsha - senlsha - hvap*evplsha;
        //cout<<it<<" "<<sabvsun<<"  "<<irabsun<<"  "<<senlsun<<"  "<<hvap<<"  "<<evplsun<<endl;
        //cout<<it<<" "<<sabvsha<<"  "<<irabsha<<"  "<<senlsha<<"  "<<hvap<<"  "<<evplsha<<endl;

        dftsunDTlsun = dirabsun_dtlsun - senlsun_dtlsun - hvap*evplsun_dtlsun;
        dftsunDTlsha = dirabsun_dtlsha - senlsun_dtlsha - hvap*evplsun_dtlsha;
        dftshaDTlsun = dirabsha_dtlsun - senlsha_dtlsun - hvap*evplsha_dtlsun;
        dftshaDTlsha = dirabsha_dtlsha - senlsha_dtlsha - hvap*evplsha_dtlsha;

        //cout<<it<<" "<<ftsun<<"  "<<dftsunDTlsun<<"  "<<dftshaDTlsun<<endl;
        //cout<<it<<" "<<ftsha<<"  "<<dftsunDTlsha<<"  "<<dftshaDTlsha<<endl;

        //-----------------------------------------------------------------------
        // difference of temperatures by quasi-newton-raphson method for the non-linear system equations
        //-----------------------------------------------------------------------
        dtlsun[it] = - (ftsun * (dftshaDTlsha-(laisha+sai)*clai/dtime) - ftsha * dftsunDTlsha)
                     / ((dftsunDTlsun-laisun*clai/dtime) * (dftshaDTlsha-(laisha+sai)*clai/dtime)
                        - dftsunDTlsha * dftshaDTlsun);
        dtlsha[it] = - (ftsun * dftshaDTlsun - ftsha * (dftsunDTlsun-laisun*clai/dtime))
                     / (dftsunDTlsha * dftshaDTlsun
                        - (dftsunDTlsun-laisun*clai/dtime) * (dftshaDTlsha-(laisha+sai)*clai/dtime));
        //cout<<it<<" "<<tlsun<<"  "<<tlsunbef<<"  "<<dtlsun[it]<<endl;
        //cout<<it<<" "<<tlsha<<"  "<<tlshabef<<"  "<<dtlsha[it]<<endl;
        if (it<itmax)
        {

            // put brakes on large temperature excursions
            if (abs(dtlsun[it])>delmax)
                dtlsun[it] = delmax*dtlsun[it]/abs(dtlsun[it]);

            if (abs(dtlsha[it])>delmax)
                dtlsha[it] = delmax*dtlsha[it]/abs(dtlsha[it]);


            if (it>1)
            {
                if ((dtlsun[it-1]*dtlsun[it])<0.)
                    dtlsun[it] = 0.5*(dtlsun[it-1] + dtlsun[it]);
                if ((dtlsha[it-1]*dtlsha[it])<0.)
                    dtlsha[it] = 0.5*(dtlsha[it-1] + dtlsha[it]);
            }
        }

        tlsun = tlsunbef + dtlsun[it];
        tlsha = tlshabef + dtlsha[it];
        //cout<<it<<" "<<tlsun<<"  "<<tlsunbef<<"  "<<dtlsun[it]<<endl;
        //cout<<it<<" "<<tlsha<<"  "<<tlshabef<<"  "<<dtlsha[it]<<endl;


        //-----------------------------------------------------------------------
        // square roots differences of temperatures and fluxes for use as the condition of convergences
        //-----------------------------------------------------------------------
        del  = std::sqrt( dtlsun[it]*dtlsun[it] + dtlsha[it]*dtlsha[it] );
        dele = pow((    dirabsun_dtlsun * dtlsun[it] +     dirabsun_dtlsha * dtlsha[it]),2)
               + pow((    dirabsha_dtlsun * dtlsun[it] +     dirabsha_dtlsha * dtlsha[it]),2)
               + pow((     senlsun_dtlsun * dtlsun[it] +      senlsun_dtlsha * dtlsha[it]),2)
               + pow((     senlsha_dtlsun * dtlsun[it] +      senlsha_dtlsha * dtlsha[it]),2)
               + pow( ( hvap*evplsun_dtlsun * dtlsun[it] +  hvap*evplsun_dtlsha * dtlsha[it]),2)
               + pow( ( hvap*evplsha_dtlsun * dtlsun[it] +  hvap*evplsha_dtlsha * dtlsha[it]),2)  ;
        dele = std::sqrt(dele);

        //-----------------------------------------------------------------------
        //  saturated vapor pressures and canopy air temperature, canopy air humidity
        //-----------------------------------------------------------------------
        // Recalculate leaf saturated vapor pressure (ei_)for updated leaf temperature
        // and adjust specific humidity (qsatl_) proportionately
        qsadv(tlsun,psrf,eisun,deisunDT,qsatlsun,qsatlsunDT);

        qsadv(tlsha,psrf,eisha,deishaDT,qsatlsha,qsatlshaDT);

        // update vegetation/ground surface temperature, canopy air temperature,
        // canopy air humidity

        taf = wta0*thm + wtg0*tg + wtlsun0*tlsun +  wtlsha0*tlsha;

        //cout<<thm<<"  "<<tg<<"  "<<tlsun<<"  "<<tlsha<<endl;

        qaf = wtaq0*qm + wtgq0*qg + wtlsunq0*qsatlsun +  wtlshaq0*qsatlsha;

        // update co2 partial pressure within canopy air
        gah2o  = 1.0/raw * tprcor/thm   ;                  // mol m-2 s-1
        pco2a = pco2m - 1.37*psrf/max(0.446,gah2o)
                * (assimsun + assimsha - respcsun - respcsha - rsoil);

        //-----------------------------------------------------------------------
        // Update monin-obukhov length and wind speed including the stability effect
        //-----------------------------------------------------------------------
        dth = thm - taf;
        dqh = qm - qaf;

        tstar = temp1*dth;
        qstar = temp2*dqh;

        thvstar = tstar + 0.61*th*qstar;
        zeta = zldis*vonkar*grav*thvstar / (ustar*ustar*thv);

        if (zeta>=0.)                            //stable
            zeta = min(2.,max(zeta,1.e-6));
        else                                             //unstable
            zeta = max(-100.,min(zeta,-1.e-6));

        obu = zldis/zeta;

        if (zeta>=0.)
            um = max(ur,0.1);
        else
        {
            wc = pow((-grav*ustar*thvstar*zii/thv),(1./3.));
            wc2 = beta*beta*wc*wc;
            um = std::sqrt(ur*ur+wc2);
        }

        if (obuold*obu <0.)
            nmozsgn = nmozsgn+1;
        if (nmozsgn>=4)
            obu = zldis/(-0.01);
        obuold = obu;

        //-----------------------------------------------------------------------
        // Test for convergence
        //-----------------------------------------------------------------------
        it++;

        if (it >itmin)
        {
            det = max(del,del2);
            del = max(dele,dele2);
            //cout<<it<<"  "<<dtmin<<"  "<<dlemin<<endl;
            if ((det<dtmin)&&(del<dlemin))
                break;
        }
    }
    // ITERATION

    // ======================================================================
    //     END stability iteration
    // ======================================================================

    z0ma = z0mv;
    zol = zeta;
    //rib = min(5.,zeta*pow(vonkar,3)*ustar*ustar/(temp1*um*um));

    // correct 27, December, 2006
    rib = min(5.,zol*ustar*ustar/(vonkar*temp1*um*um));


    // canopy fluxes and total assimilation amd respiration

    if (lai> 0.001)
        rst = 1./(laisun/rssun + laisha/rssha);
    else
    {
        rssun = 2.0e+4 ;
        rssha = 2.0e+4;
        assimsun = 0. ;
        assimsha = 0.;
        respcsun = 0. ;
        respcsha = 0.;
        rst = 2.0e+4;
    }
    assim = assimsun + assimsha;
    respc = respcsun + respcsha + rsoil;

    fsenl = senlsun + senlsha + senlsun_dtlsun*dtlsun[it-1] + senlsun_dtlsha*dtlsha[it-1]
            + senlsha_dtlsun*dtlsun[it-1] + senlsha_dtlsha*dtlsha[it-1];
    etr = etrsun + etrsha     +  etrsun_dtlsun*dtlsun[it-1] +  etrsun_dtlsha*dtlsha[it-1]
          +  etrsha_dtlsun*dtlsun[it-1] +  etrsha_dtlsha*dtlsha[it-1];
    evplwet = evplwetsun + evplwetsun_dtlsun*dtlsun[it-1] + evplwetsun_dtlsha*dtlsha[it-1]
              + evplwetsha + evplwetsha_dtlsun*dtlsun[it-1] + evplwetsha_dtlsha*dtlsha[it-1];
    fevpl = evplsun + evplsha + evplsun_dtlsun*dtlsun[it-1] + evplsun_dtlsha*dtlsha[it-1]
            + evplsha_dtlsun*dtlsun[it-1] + evplsha_dtlsha*dtlsha[it-1];

    elwmax = ldew/dtime;
    elwdif = max(0., evplwet-elwmax);
    evplwet = min (evplwet, elwmax);

    fevpl = fevpl - elwdif;
    fsenl = fsenl + hvap*elwdif;

    // wind stresses
    taux  = taux - sigf*rhoair*us/ram;
    tauy  = tauy - sigf*rhoair*vs/ram;

    //-----------------------------------------------------------------------
    // fluxes from ground to canopy space
    //-----------------------------------------------------------------------

    fseng = fseng + cpair*rhoair*cgh*(tg-taf);
    fevpg = fevpg +    rhoair*cgw*(qg-qaf);

    //status1 = _isnan(fseng);
    //status2 = _isnan(fevpg);
    //if((status1!=0)||(status2!=0))
    //{
    //	cout<<"error"<<endl;
    //	getchar();
    //}

    //fseng = max(fseng,0.0);
    //fevpg = max(fevpg,0.0);
    //-----------------------------------------------------------------------
    // downward (upward) intwave radiation below (above) the canopy
    //-----------------------------------------------------------------------
    dlrad = sigf * thermk * frl
            + stefnc * ( fac1*pow(tlsunbef,3)*(tlsunbef+4.*dtlsun[it-1])
                         + fac2*pow(tlshabef,3)*(tlshabef+4.*dtlsha[it-1]) );
    ulrad = stefnc * ( fac1*pow(tlsunbef,3)*(tlsunbef+4.*dtlsun[it-1])
                       + fac2*pow(tlshabef,3)*(tlshabef+4.*dtlsha[it-1])
                       + sigf*emg*thermk*pow(tg,4) ) ;

    //-----------------------------------------------------------------------
    // Derivative of soil energy flux with respect to soil temperature (cgrnd)
    //-----------------------------------------------------------------------
    cgrnds = cgrnds + cpair*rhoair*cgh*(1.-wtg0);
    cgrndl = cgrndl + rhoair*cgw*(1.-wtgq0)*dqgdT;
    cgrnd  = cgrnds + cgrndl*htvp;

    //-----------------------------------------------------------------------
    // balance check
    //-----------------------------------------------------------------------
    err = sabvsun+sabvsha
          + irabsun+dirabsun_dtlsun*dtlsun[it-1]+dirabsun_dtlsha*dtlsha[it-1]
          + irabsha+dirabsha_dtlsun*dtlsun[it-1]+dirabsha_dtlsha*dtlsha[it-1]
          - fsenl-hvap*fevpl;

    //if(fabs(err)>0.3)
    //	printf("leaftemtwo.F90 : energy imbalance %f, %d\n",err, it-1);

    //-----------------------------------------------------------------------
    // Update dew accumulation (kg/m2)
    //-----------------------------------------------------------------------

    //if (evplwet>0.0)
    ldew = max(0.,ldew-evplwet*dtime);

    //-----------------------------------------------------------------------
    // 2 m height air temperature
    //-----------------------------------------------------------------------
    tref = tref + sigf*(thm + temp1*dth * (1./temp12m - 1./temp1));
    qref = qref + sigf*( qm + temp2*dqh * (1./temp22m - 1./temp2));
    delete []dtlsun;
    delete [] dtlsha;
}
void CommonLandModel::CLMMAIN (double dtime, bool doalb, bool dolai, bool dosst,
                               int nl_soil, int maxsnl, double dlon, double dlat,
                               int itypwat, int ivt, double &oro,

                               // soil information
                               double albsol, double csol[MAXSOILL], double porsl[MAXSOILL], double phi0[MAXSOILL], double bsw[MAXSOILL],
                               double dkmg[MAXSOILL], double dksatu[MAXSOILL], double dkdry[MAXSOILL], double hksati[MAXSOILL],

                               // vegetation information
                               double z0m, double displa, double sqrtdi,
                               double effcon, double vmax25, double slti,
                               double hlti, double shti, double hhti,
                               double trda, double trdm, double trop,
                               double gradm, double binter, double extkn,
                               double chil, double ref[2][2], double tran[2][2], double rootfr[MAXSOILL],

                               // atmospheric forcing
                               double frl, double sols, double soll, double solsd, double solld,
                               double pco2m, double po2m, double us, double vs, double tm, double qm,
                               double prc, double prl, double psrf, double rhoair,
                               double hu, double ht, double hq,

                               // land surface variables required for restart
                               int year, int jday, int msec,
                               double z[MAXSNL+MAXSOILL], double dz[MAXSNL+MAXSOILL],
                               double tss[MAXSNL+MAXSOILL], double wliq[MAXSNL+MAXSOILL], double wice[MAXSNL+MAXSOILL],
                               double &tg, double &tlsun, double &tlsha,
                               double &ldew, double &sag, double &scv, double &snowdp,
                               double &fveg, double &fsno, double &sigf,
                               double &green, double &lai, double &sai,
                               double &coszen, double albg[2][2], double albv[2][2],
                               double alb[2][2], double ssun[2][2], double ssha[2][2], double &thermk,
                               double &extkb, double &extkd,

                               // fluxes
                               double &taux,  double &tauy,
                               double &fsena, double &fevpa, double &lfevpa,
                               double &fsenl, double &fevpl, double &etr,
                               double &fseng, double &fevpg, double &olrg,
                               double &fgrnd, double &trad, double &tref, double &qref,
                               double &rsur, double &rnof, double &rst,
                               double &assim, double &respc,
                               double &parsun, double &parsha, double &sabvsun,
                               double &sabvsha, double &sabg, double &sabvg, double &xerr, double &zerr,

                               // TUNABLE modle constants
                               double &zlnd,   double &zsno,   double &csoilc, double &dewmx,  double &wtfact,
                               double &capr,   double &cnfac,  double &ssi,    double &wimp,   double &pondmx,
                               double &smpmax, double &smpmin, double &trsmx0, double &tcrit,

#ifndef EcoDynamics
                               // time-varying vegetation parameters from read-in file
                               double &lai_r, double &sai_r, double &green_r, double &fveg_r,
#endif

                               // additional variables required by coupling with WRF model
                               double &emis, double &z0ma, double &zol,
                               double &rib, double &ustar, double &qstar, double &tstar,
                               double &u10m, double &v10m, double &f10m, double &fm, double &fh, double &fq,
                               int    &snl)
{
    // ----------------------- Local  Variables -----------------------------
    double calday  ; // Julian cal day (1.xx to 365.xx)
    double endwb      ; // water mass at the end of time step
    double errore     ; // energy balnce errore (Wm-2)
    double errorw     ; // water balnce errore (mm)
    double fiold[MAXSNL+MAXSOILL]; // fraction of ice relative to the total water
    //	double obs_coszen ; // cosine of the solar zenith angle
    double pg         ; // water onto ground including canopy runoff [kg/(m2 s)]
    double pg_rain    ; // liquid water onto ground [kg/(m2 s)]
    double pg_snow    ; // ice onto ground [kg/(m2 s)]
    double qseva      ; // ground surface evaporation rate (mm h2o/s)
    double qsdew      ; // ground surface dew formation (mm h2o /s) [+]
    double qsubl      ; // sublimation rate from snow pack (mm h2o /s) [+]
    double qfros      ; // surface dew added to snow pack (mm h2o /s) [+]
    double rootr[MAXSOILL]      ; // root resistance of a layer, all layers add to 1.0
    double scvold     ; // snow cover for previous time step [mm]
    double sm         ; // rate of snowmelt [kg/(m2 s)]
    double ssw        ; // water volumetric content of soil surface layer [m3/m3]
    double tssub[7]   ; // surface/sub-surface temperatures [K]
    double tssea      ; // sea surface temperature [K]
    double totwb      ; // water mass at the begining of time step
    double wt         ; // fraction of vegetation buried (covered) by snow [-]
    double zi[MAXSOILL+MAXSNL+1]; // interface level below a "z" level (m)

    //int snl      ;      // number of snow layers
    int imelt[MAXSOILL+MAXSNL];	  // flag for: melting=1, freezing=2, Nothing happended=0
    int lb ;           // lower bound of arrays
    int j;             // do looping index


    double snowdp_old;
    double my_temp = 0.0;

    /* 暂时取消控制，以进行调式。20090831
    if ((tg-tm)>30.0)
        tm = tg - 15.0;
        */
    //======================================================================
    //  [1] Solar absorbed by vegetation and ground
    //======================================================================

//	if(jday==4&&msec == 10*3600)
//		getchar();

    //cout<<"[1] Solar absorbed by vegetation and ground"<<endl;
    netsolar (itypwat,sigf,albg,albv,alb,ssun,ssha,
              sols,soll,solsd,solld,
              parsun,parsha,sabvsun,sabvsha,sabg,sabvg);
    //cout<<"solar absorbed"<<endl;
    //cout<<parsun<<"  "<<parsha<<" "<<sabvsun<<"  "<<sabvsha<<"  "<<sabg<<" "<<sabvg<<endl;
    //getchar();
    //cout<<itypwat<<endl;
    //getchar();
    //======================================================================
    if (itypwat<=3)    // <=== not lake and ocean (itypwat = 0, 1, 2, 3)
    {
        //======================================================================
        //initial set

        scvold = scv;        //snow mass at previous time step
        snowdp_old = snowdp;
        snl = MAXSNL-1;

        //calculate snow layer number
        for (j=0; j<MAXSNL; j++)
        {
            if ((wliq[j]+wice[j])>0.0)
                snl=snl-1;
        }

        zi[MAXSNL] = 0;

        if (snl<MAXSNL)
        {
            for (j = MAXSNL-1; j>=snl; j--)
            {
                zi[j]=zi[j+1]-dz[j];
                //cout<<j<<"  "<<zi[j]<<endl;
            }

        }

        for (int j=MAXSNL+1; j<=MAXSNL+nl_soil; j++)
        {
            zi[j]=zi[(j-1)]+dz[j-1];
            //cout<<j<<"  "<<zi[j]<<endl;
        }
        double temp=0.0;
        for (int i=0; i<nl_soil; i++)
        {
            temp += wice[MAXSNL+i]+wliq[MAXSNL+i];
        }
        totwb = ldew + scv + temp;
        if (snl<MAXSNL-1)
        {
            for (int i=snl+1; i<MAXSNL; i++)
            {
                fiold[i] = wice[i]/(wliq[i]+wice[i]);
            }
        }
        //----------------------------------------------------------------------
        // [2] Canopy interception and precipitation onto ground surface
        //----------------------------------------------------------------------
        //cout<<"[2] Canopy interception and precipitation onto ground surface"<<endl;
        leafinterception (dtime,dewmx,chil,prc,prl,tm,scv,sigf,lai,sai,ldew,pg);
        //cout<<"prc  prl  pg :   "<<prc<<" "<<prl<<" "<<pg<<endl;
        //cout<<dewmx<<" "<<chil<<" "<<prc<<" "<<prl<<endl;
        //cout<<tm<<"  "<<scv<<" "<<sigf<<" "<<lai<<endl;
        // cout<<sai<<"  "<<ldew<<"  "<<pg<<endl;

        //----------------------------------------------------------------------
        // [3] Initilize new snow nodes for snowfall / sleet
        //----------------------------------------------------------------------
        //cout<<"[3] Initilize new snow nodes for snowfall // slee"<<endl;
        newsnow (itypwat,maxsnl,dtime,tm,tg,pg,tcrit,
                 zi,z,dz,tss,wliq,wice,fiold,
                 snl,sag,scv,snowdp,pg_rain,pg_snow);
        // added by huangcl
        /*
        if(snowdp>MAXSNOWDP)
        {
        	dz[maxsnl-snl]   = dz[maxsnl-snl] - (snowdp - MAXSNOWDP);
        	zi[maxsnl-snl]   = -MAXSNOWDP;
        	z[maxsnl-snl]    = zi[maxsnl-snl] + 0.5*dz[maxsnl-snl];
        	wice[maxsnl-snl] = wice[maxsnl-snl] - (snowdp - MAXSNOWDP)*1000.0;
        	snowdp  = MAXSNOWDP;
        	scv              = snowdp*250.0;
        }

          */
        //----------------------------------------------------------------------
        // [4] Energy AND Water balance
        //----------------------------------------------------------------------
        //cout<<"[4] Energy AND Water balance"<<endl;
        lb  = snl+1;	     // lower bound of array
        THERMAL(itypwat    ,lb       ,nl_soil  ,dtime   ,trsmx0   ,
                zlnd       ,zsno     ,csoilc   ,dewmx   ,capr     ,
                cnfac      ,csol     ,porsl    ,phi0    ,bsw      ,
                dkmg       ,dkdry    ,dksatu   ,lai     ,sai      ,
                z0m        ,displa   ,sqrtdi   ,rootfr  ,effcon   ,
                vmax25     ,slti     ,hlti     ,shti    ,hhti     ,
                trda       ,trdm     ,trop     ,gradm   ,binter   ,
                extkn      ,hu       ,ht       ,hq      ,us       ,
                vs         ,tm       ,qm       ,rhoair  ,psrf     ,
                pco2m      ,po2m     ,coszen   ,parsun  ,parsha   ,
                sabvsun    ,sabvsha  ,sabg     ,frl     ,extkb    ,
                extkd      ,thermk   ,fsno     ,sigf    ,dz  ,
                z          ,zi       ,tlsun    ,tlsha   ,tss ,
                wice       ,wliq     ,ldew     ,scv     ,snowdp   ,
                imelt      ,taux     ,tauy     ,fsena   ,fevpa    ,
                lfevpa     ,fsenl    ,fevpl    ,etr     ,fseng    ,
                fevpg      ,olrg     ,fgrnd    ,rootr   ,qseva    ,
                qsdew      ,qsubl    ,qfros    ,sm      ,tref     ,
                qref       ,trad     ,rst      ,assim   ,respc    ,
                errore     ,emis     ,z0ma     ,zol     ,rib      ,
                ustar      ,qstar    ,tstar    ,u10m    ,v10m     ,
                f10m       ,fm       ,fh       ,fq  );
        //	if(tss[MAXSNL]<0)
        //	{
        //		cout<<tss[MAXSNL]<<endl;
        //		getchar();
        //	}
        WATER ( itypwat     ,lb       ,nl_soil ,dtime    ,
                z           ,dz       ,zi      ,bsw     ,porsl    ,
                phi0       ,hksati   ,rootr    ,tss     ,wliq     ,
                wice       ,pg_rain  ,sm       ,etr     ,qseva    ,
                qsdew      ,qsubl    ,qfros    ,rsur    ,rnof     ,
                wtfact     ,pondmx   ,ssi      ,wimp    ,smpmin   );
        //stop here huang //4.26
        if (snl<MAXSNL-1)
        {
            // Compaction rate for snow
            // Natural compaction and metamorphosis. The compaction rate
            // is recalculated for every new timestep
            lb  = snl+1;  			// lower bound of array
            snowcompaction (lb,dtime,imelt,fiold,tss, wliq,wice,dz);


            // Combine thin snow elements
            lb = 0;
            //cout<<"combine  in: "<<snowdp<<endl;
            snowlayerscombine(lb,snl,z,dz,zi,
                              wliq,wice,tss,scv,snowdp);

            //cout<<"combine  out: "<<snowdp<<endl;
            // Divide thick snow elements
            if (snl<MAXSNL-1)
                snowlayersdivide (lb,snl,z,dz,zi,
                                  wliq,wice,tss);
        }
        // Set zero to the empty node
        if (snl==MAXSNL-1) sag=0.0;
        for (int kk=0; kk<=snl; kk++)
        {
            wice[kk]=0.0;
            wliq[kk]=0.0;
            tss [kk]=0.0;
            z[kk]   =0.0;
            dz[kk]  =0.0;
        }

        lb=snl+1;
        tg = tss[lb];//MAXSNL-snl
        ssw = min(1.,1.e-3*wliq[MAXSNL]/dz[MAXSNL]);
        // ----------------------------------------
        // Update the snow age
        // ----------------------------------------
        snowage (dtime, tg, scv, scvold, sag);
        // ----------------------------------------
        // energy balance
        // ----------------------------------------
        zerr=errore;
        //if(fabs(errore)>0.2)
        //	cout<<"Warning: energy balance violation "<<errore<<"   <<ivt<<endl;

        // ----------------------------------------
        // water balance
        // ----------------------------------------

        /*        double wliq_max,wice_max,wliq_error,wice_error;

                for (int i=nl_soil-1;i>=0;i--)
                {
                    wliq_max = dz[MAXSNL+i]*porsl[i]*denh2o;
                    wice_max = dz[MAXSNL+i]*porsl[i]*denice;
                    //cout<<porsl[i]<<" "<<dz[MAXSNL+i]<<"  "<<wliq_max<<"  "<<wice_max<<endl;
                    if (i>0)
                    {
                        if (wliq[MAXSNL+i]>wliq_max)
                        {
                            wliq_error  = wliq[MAXSNL+i] - wliq_max;
                            wliq[MAXSNL+i] = wliq_max;
                            wliq[MAXSNL+i-1] = wliq[MAXSNL+i-1] + wliq_error;
                        }

                        if (wice[MAXSNL+i]>wice_max)
                        {
                            wice_error  = wice[MAXSNL+i] - wice_max;
                            wice[MAXSNL+i] = wice_max;
                            wice[MAXSNL+i-1] = wice[MAXSNL+i-1] + wice_error;
                        }
                    }
                    else   // i==0
                    {
                        if (wliq[MAXSNL]>wliq_max)
                        {
                            wliq_error  = wliq[MAXSNL] - wliq_max;
                            wliq[MAXSNL] = wliq_max;
                        }

                        if (wice[MAXSNL]>wice_max)
                        {
                            wice_error  = wice[MAXSNL] - wice_max;
                            wice[MAXSNL+i] = wice_max;
                        }
                        rnof       = rnof + wice_error + wliq_error;
                    }
                }
        */
        endwb = 0.0;
        for (int j=0; j<nl_soil; j++)
        {
            endwb+= wice[MAXSNL+j]+wliq[MAXSNL+j];
        }
        endwb = endwb + ldew+scv;
        errorw=(endwb-totwb)-(prc+prl-fevpa-rnof)*dtime;
        if (itypwat>1) errorw=0.0;      //wetland, glacier
        xerr=errorw/dtime;
        // if(abs(errorw)>1.e-3) write(6,*) 'Warning: water balance violation', errorw,ivt

        // ----------------------------------------
        // temperature check
        // ----------------------------------------
        /*        for (int i=MAXSNL-snl;i<MAXSNL;i++)
                {
                    if (tss[i]>273.16)
                        tss[i] = 273.16;

                    if ((tss[i]-273)<-30)
                        tss[i] = 250;
                }

                for (int i=MAXSNL;i<MAXSNL+nl_soil;i++)
                {
                    if ((tss[i]-280)>40)
                        tss[i] = 320;

                    if ((tss[i]-273)<-30)
                        tss[i] = 250;
                }

                if ((tg-280)>40)
                    tg = 320;

                if ((tg-273)<-30)
                    tg = 250;

                if (fabs(tlsun-280)>40)
                    tlsun = 280;

                if (fabs(tlsha-280)>40)
                    tlsha = 280;
        */
        // ----------------------------------------
        // snow check
        // ----------------------------------------
        /*
        		if(snowdp>1.0)
        		{
        			snowdp = 1.0;
        			scv    = 800.0;

        			z[MAXSNL-1]  = -0.5;
        			dz[MAXSNL-1] = 1.0;
        		}
        */
        if (snowdp<0.0)
        {
            snowdp = 0.0;
            scv    = 0.0;

            z[MAXSNL-1]  = 0.0;
            dz[MAXSNL-1] = 0.0;
        }
    }
    //======================================================================

    else if (itypwat <= 5) // <=== is lake (itypwat = 4 or 5)
        //======================================================================
    {
        scvold = scv ;         // snow mass at previous time step
        zi[MAXSNL] =0.0;
        for (j = 1; j<=nl_soil; j++)
            zi[MAXSNL+j] =zi[MAXSNL+j-1]+dz[MAXSNL-1+j];

        lake(nl_soil,itypwat,dlat,dtime,z,
             dz,     zi,     hu,  ht,   hq,
             us,     vs,     tm,  qm,   prc,
             prl,    rhoair, psrf,sabg, frl,
             tg,     tss,    wliq,wice, scv,
             snowdp, trad,   tref,qref, taux,
             tauy,   fsena,  fevpa,lfevpa,fseng,
             fevpg,  olrg,   fgrnd,tcrit, emis,
             z0ma,   zol,    rib,   ustar,qstar,
             tstar,  u10m,   v10m,  f10m,  fm,
             fh,    fq);

        //cout<<"tg1:  "<<tg<<endl;
        // null data for lake component
        snl = 0;
        for (j=0; j<MAXSNL; j++)
        {
            z[j]  = 0.0;
            dz[j] = 0.0;
            tss[j]= 0.0;
        }
        for (j=0; j<MAXSNL+MAXSOILL; j++)
        {
            wliq[j] = 0.0;
            wice[j] = 0.0;
        }
        tlsun = tm;
        tlsha = tm;
        ldew  = 0.0;

        fsenl = 0.0;
        fevpl = 0.0;
        etr   = 0.0;
        rsur  = 0.0;
        rnof  = 0.0;
        rst   = -9999.;
        assim = 0.0;
        respc = 0.0;
        zerr=0.;
        xerr=0.;

        // Update the snow age
        snowage (dtime, tg, scv, scvold, sag);
        ssw   = 0.0;
    }
    //======================================================================
    else                     // <=== is ocean (itypwat >= 99)
        //======================================================================
        // simple ocean-sea ice model
    {
        tssea = tg;
        for (j=0; j<7; j++)
            tssub[j] = tss[MAXSNL+j];
        SOCEAN(dosst,  dtime,  oro,   hu,    ht,
               hq,     us,     vs,    tm,    qm,
               rhoair, psrf,   sabg,  frl,   tssea,
               tssub,  scv,    taux,  tauy,  fsena,
               fevpa,  lfevpa, fseng, fevpg, tref,
               qref,   z0ma,   zol,   rib,   ustar,
               qstar,  tstar,  u10m,  v10m,  f10m,
               fm,     fh,     fq,    emis,  olrg);

        // null data for sea component
        for (j=0; j<MAXSNL+MAXSOILL; j++)
        {
            z[j]    = 0.0;
            dz[j] = 0.0;
            tss[j] = 0.0;
            wliq[j] = 0.0;
            wice[j] = 0.0;
        }

        for (j=0; j<7; j++)
            tss[MAXSNL+j] = tssub[j] ;
        tg      = tssea;
        tlsun   = tm;
        tlsha   = tm;
        ldew    = 0.0;
        sag     = 0.0;
        snowdp  = scv/1000.*20.;

        trad    = tssea;
        fsenl   = 0.0;
        fevpl   = 0.0 ;
        etr     = 0.0;
        fgrnd   = 0.0;
        rsur    = 0.0;
        rnof    = 0.0;
        rst     = -9999.;
        assim   = 0.0;
        respc   = 0.0;

        xerr    = 0.0;
        zerr    = 0.0;
    }

    //======================================================================
    //======================================================================
    // Preparation for the next time step
    // 1) time-varying parameters for vegatation
    // 2) fraction of snow cover
    // 3) solar zenith angle and
    // 4) albedos
    //======================================================================

    // cosine of solar zenith angle
    calday = (double)jday + (double)(msec/86400.0);
    coszen = orb_coszen(calday,dlon,dlat);
    if (itypwat <= 5)
        // land grid
    {
        // need to update lai and sai, fveg, green, they are done once in a day only
#if (defined EcoDynamics)
        if (dolai)
            // call EcoModel(ivt)
            lai_empirical(ivt,nl_soil,rootfr,tss,lai,sai,fveg,green);
#else
        {
            lai=lai_r;
            sai=sai_r;
            green=green_r;
            fveg=fveg_r;
        }
#endif
        // fraction of snow cover.
        snowfraction (fveg,z0m,snowdp,wt,sigf,fsno);

        // albedos
        // we supposed call it every time-step, just because
        // other vegeation related parameters are needed to create
        //*if(doalb)then
        albland (itypwat,albsol,  chil,  ref,   tran,
                 fveg,    green,  lai,   sai,   coszen,
                 wt,      fsno,   scv,   sag,   ssw,
                 tg,      alb,    albg,  albv,  ssun,
                 ssha,    thermk, extkb, extkd);

    }
    //*endif
    else                   // ocean grid
    {

        //*if(doalb)then
        albocean (oro,scv,coszen,alb);
        //*endif

        // null data for sea component
        lai = 0.0;
        sai = 0.0;
        green = 0.0;
        fveg = 0.0;
        sigf = 0.0;
        fsno = 0.0;

        for (j=0; j<2; j++)
            for (int k=0; k<2; k++)
            {
                albg[j][k] = alb[j][k];
                albv[j][k] = 0.0;
                ssun[j][k] = 0.0;
                ssha[j][k] = 0.0;
            }
        thermk = 0.0;
        extkb = 0.0;
        extkd = 0.0;
    }
    //----------------------------------------------------------------------
}
void CommonLandModel::meltf (int   lb,                    int    nl_soil,               double dtime,                  double fact[MAXSNL+MAXSOILL], double brr[MAXSNL+MAXSOILL] ,
                             double hs,                   double dhsdT,                  double tssbef[MAXSNL+MAXSOILL],double tss[MAXSNL+MAXSOILL],  double wliq[MAXSNL+MAXSOILL],
                             double wice[MAXSNL+MAXSOILL],int    imelt[MAXSNL+MAXSOILL], double &scv,                    double &snowdp, double &sm, double &xmf)
{

    //-----------------------------------------------------------------------
    // Original author : Yongjiu Dai, September 15, 1999
    //
    // calculation of the phase change within snow and soil layers:
    //
    // (1) check the conditions which the phase change may take place,
    //     i.e., the layer temperature is great than the freezing point
    //     and the ice mass is not equal to zero (i.e., melting),
    //     or layer temperature is less than the freezing point
    //     and the liquid water mass is not equal to zero (i.e., freezing);
    // (2) assess the rate of phase change from the energy excess (or deficit)
    //     after setting the layer temperature to freezing point;
    // (3) re-adjust the ice and liquid mass, and the layer temperature
    //
    //-----------------------------------------------------------------------


    //-----------------------------------------------------------------------

    // integer, INTENT(in) :: nl_soil             // upper bound of array (i.e., soil layers)
    // integer, INTENT(in) :: lb                  // lower bound of array (i.e., snl +1)
    //real(r8), INTENT(in) :: dtime               // time step [second]
    //real(r8), INTENT(in) :: tssbef(lb:nl_soil)  // temperature at previous time step [K]
    //real(r8), INTENT(in) :: brr   (lb:nl_soil)  //
    //real(r8), INTENT(in) :: fact  (lb:nl_soil)  // temporary variables
    //real(r8), INTENT(in) :: hs                  // net ground heat flux into the surface
    //real(r8), INTENT(in) :: dhsdT               // temperature derivative of "hs"

    //real(r8), INTENT(inout) :: tss (lb:nl_soil) // temperature at current time step [K]
    //real(r8), INTENT(inout) :: wice(lb:nl_soil) // ice lens [kg/m2]
    //real(r8), INTENT(inout) :: wliq(lb:nl_soil) // liquid water [kg/m2]
    //real(r8), INTENT(inout) :: scv              // snow mass [kg/m2]
    //real(r8), INTENT(inout) :: snowdp           // snow depth [m]

    //real(r8), INTENT(out) :: sm                 // rate of snowmelt [mm/s, kg/(m2 s)]
    //real(r8), INTENT(out) :: xmf                // total latent heat of phase change
    // integer, INTENT(out) :: imelt(lb:nl_soil)  // flag for melting or freezing [-]

    // Local
    double  hm[MAXSNL+MAXSOILL];                  // energy residual [W/m2]
    double  xm[MAXSNL+MAXSOILL];                  // metling or freezing within a time step [kg/m2]
    double  heatr;                                // energy residual or loss after melting or freezing
    double  temp1;                                // temporary variables [kg/m2]
    //	double  temp2 ;                               // temporary variables [kg/m2]

    double wmass0[MAXSNL+MAXSOILL];
    double wice0[MAXSNL+MAXSOILL];
    double wliq0[MAXSNL+MAXSOILL];
    double  propor,tinc, we, scvold;
    int    j;

    //-----------------------------------------------------------------------

    sm = 0.;
    xmf = 0.;
    for (j= lb; j<MAXSNL+nl_soil; j++)
    {
        imelt[j] = 0;
        hm[j] = 0.;
        xm[j] = 0.;
        wice0[j] = wice[j];
        wliq0[j] = wliq[j];
        wmass0[j] = wice[j] + wliq[j];
    }
    scvold=scv;
    we=0.;
    if (lb<MAXSNL)
    {
        for (j=lb; j<MAXSNL; j++)
            we = we + wice[j]+wliq[j];
    }


    for (j = lb; j<MAXSNL+nl_soil; j++)

        // Melting identification
        // if ice exists above melt point, melt some to liquid.
    {
        if ((wice[j] > 0)&&(tss[j] > tfrz))
        {
            imelt[j] = 1;
            tss[j] = tfrz;
        }
        // Freezing identification
        // if liquid exists below melt point, freeze some to ice.

        //×××××××××××××××××××××××××××××××××××××××××××××××××××××××

        if ((wliq[j] > 0.0)&&(tss[j] < tfrz))
        {
            imelt[j] = 2;
            tss[j] = tfrz;
        }

        //×××××××××××××××××××××××××××××××××××××××××××××××××××××××
    }

    // If snow exists, but its thickness less than the critical value (0.01 m)

    if ((lb == MAXSNL)&&(scv > 0.0))
    {
        if (tss[MAXSNL]> tfrz)
        {
            imelt[MAXSNL] = 1;
            tss[MAXSNL] = tfrz;
        }
    }

    // Calculate the energy surplus and loss for melting and freezing
    for (j = lb; j<MAXSNL+nl_soil; j++)
    {
        if (imelt[j]> 0)
        {
            tinc = tss[j]-tssbef[j];
            if (j > lb)
                hm[j] = brr[j] - tinc/fact[j] ;
            else
                hm[j] = hs + dhsdT*tinc + brr[j] - tinc/fact[j] ;
        }
    }
    for (j = lb; j<MAXSNL+nl_soil; j++)
    {
        if ((imelt[j] == 1)&&(hm[j] < 0.))
        {
            hm[j] = 0.;
            imelt[j] = 0;
        }

        // this error was checked carefully, it results from the the computed error
        // of "Tridiagonal-Matrix" in subroutine "thermal".

        if ((imelt[j] == 2)&&(hm[j] > 0.))
        {
            hm[j] = 0.;
            imelt[j] = 0;
        }
    }

    // The rate of melting and freezing

    for (j = lb; j<MAXSNL+nl_soil; j++)
    {
        if ((imelt[j] > 0)&&(abs(hm[j]) > 0.0))
        {
            xm[j] = hm[j]*dtime/hfus;                        // kg/m2
            // if snow exists, but its thickness less than the critical value (1 cm)
            // Note: more work is need on how to tune the snow depth at this case
            if (j == MAXSNL)
            {
                if ((lb == MAXSNL)&&(scv > 0.)&&(xm[j] > 0.))
                {
                    temp1 = scv ;                                // kg/m2
                    scv = max(0.,temp1-xm[j]);
                    propor = scv/temp1;
                    if (propor<1.0)
                        snowdp = propor * snowdp;
                    heatr = hm[j] - hfus*(temp1-scv)/dtime ;      // W/m2
                    if (heatr > 0.)
                    {
                        xm[j] = heatr*dtime/hfus;                  // kg/m2
                        hm[j] = heatr;                            // W/m2
                    }
                    else
                    {
                        xm[j] = 0.;
                        hm[j] = 0.;
                    }
                    sm = max(0.,(temp1-scv))/dtime;              // kg/(m2 s)
                    xmf = hfus*sm;
                }
            }
            heatr = 0.;

            if (xm[j] > 0.)
            {
                wice[j] = max(0., wice0[j]-xm[j]);
                heatr = hm[j] - hfus*(wice0[j]-wice[j])/dtime;
            }
            else if (xm[j] < 0.)
            {
                wice[j] = min(wmass0[j], wice0[j]-xm[j]);
                wice[j] = max(0.0, wice[j]);
                heatr = hm[j] - hfus*(wice0[j]-wice[j])/dtime  ;
            }
            wliq[j] = max(0.,wmass0[j]-wice[j]);
            //if (abs(heatr) <500.)//为什么？wlx20090902
            if(abs(heatr) >0.)
            {
                //cout<<j<<" "<<fact[j]<<" "<<heatr<<" "<<wliq[j]<<" "<<wice[j]<<endl;
                if (j > lb)
                    tss[j] = tss[j] + fact[j]*heatr;
                else
                    tss[j] = tss[j] + fact[j]*heatr/(1.-fact[j]*dhsdT);
                if ((wliq[j]*wice[j]) > 0.)
                    tss[j] = tfrz;
            }
            xmf = xmf + hfus * (wice0[j]-wice[j])/dtime;
            if ((imelt[j] == 1)&&(j <MAXSNL) )
                sm = sm + max(0.,(wice0[j]-wice[j]))/dtime;
        }
    }

    //scvold=scv
    if (lb<=MAXSNL-1)
    {
        double temp= 0.0;
        for (j=lb; j<MAXSNL; j++)
            temp = temp + wice[j]+wliq[j];
        we = temp - we;
        if (abs(we)>1.e-6)
        {
            std::cout<<"meltf err   "<<we<<std::endl;
            //exit(0);
        }
    }
}
void CommonLandModel::moninobuk(double hu,   double ht,   double hq,     double displa, double z0m,
                                double z0h,  double z0q,  double obu,    double um,     double &ustar,
                                double &temp1,double &temp2,double &temp12m,double &temp22m,double &f10m,
                                double &fm,   double &fh,   double &fq)
{

    //// ======================================================================
    //// Original author : Yongjiu Dai, September 15, 1999
    ////
    //// calculation of friction velocity, relation for potential temperatur
    //// and humidity profiles of surface boundary layer.
    //// the scheme is based on the work of Zeng et al. (1998):
    //// Intercomparison of bulk aerodynamic algorithms for the computation
    //// of sea surface fluxes using TOGA CORE and TAO data. J. Climate, Vol. 11: 2628-2644
    //// ======================================================================
    //
    //
    //// ---------------------- dummy argument --------------------------------
    //
    //  double , INTENT(in) :: hu       // observational height of wind [m]
    //  double , INTENT(in) :: ht       // observational height of temperature [m]
    //  double , INTENT(in) :: hq       // observational height of humidity [m]
    //  double , INTENT(in) :: displa   // displacement height [m]
    //  double , INTENT(in) :: z0m      // roughness length, momentum [m]
    //  double , INTENT(in) :: z0h      // roughness length, sensible heat [m]
    //  double , INTENT(in) :: z0q      // roughness length, latent heat [m]
    //  double , INTENT(in) :: obu      // monin-obukhov length (m)
    //  double , INTENT(in) :: um       // wind speed including the stablity effect [m/s]
    //
    //  double , INTENT(out) :: ustar   // friction velocity [m/s]
    //  double , INTENT(out) :: temp1   // relation for potential temperature profile
    //  double , INTENT(out) :: temp2   // relation for specific humidity profile
    //  double , INTENT(out) :: temp12m // relation for temperature at 2m
    //  double , INTENT(out) :: temp22m // relation for specific humidity at 2m
    //  double , INTENT(out) :: f10m    // integral of profile function for momentum at 10m
    //  double , INTENT(out) :: fm      // integral of profile function for momentum
    //  double , INTENT(out) :: fh      // integral of profile function for heat
    //  double , INTENT(out) :: fq      // integral of profile function for moisture

    //------------------------ local variables ------------------------------

    double  zldis;  // reference height "minus" zero displacement heght [m]
    //double  psi;    // stability function for unstable case
    double  zetam;  // transition point of flux-gradient relation (wind profile)
    double  zetat;  // transition point of flux-gradient relation (temp. profile)
    double  zeta;   // dimensionless height used in Monin-Obukhov theory

    //-----------------------------------------------------------------------
    // adjustment factors for unstable (moz < 0) or stable (moz > 0) conditions.

    // wind profile
    zldis=hu-displa;
    zeta=zldis/obu;
    zetam=1.574;
    if (zeta < -zetam)          // zeta < -1
    {
        fm    =   std::log(-zetam*obu/z0m) - psi(1,-zetam) + psi(1,z0m/obu)+ 1.14*(pow((-zeta),0.333)-pow((zetam),0.333));
        ustar = vonkar*um/fm;
    }
    else
    {
        if (zeta < 0.0)         // -1 <= zeta < 0
        {
            fm    = std::log(zldis/z0m) - psi(1,zeta) + psi(1,z0m/obu);
            ustar = vonkar*um/fm;
        }
        else if (zeta <= 1.0)       //  0 <= ztea <= 1
        {
            fm    = std::log(zldis/z0m) + 5.*zeta - 5.*z0m/obu ;
            ustar = vonkar*um/fm;
        }
        else                            //  1 < zeta, phi=5+zeta
        {
            fm    = std::log(obu/z0m) + 5.0 - 5.0*z0m/obu + (5.0*std::log(zeta)+zeta-1.0);
            ustar = vonkar*um/fm;
        }
    }

    // for 10 meter wind-velocity
    zldis=10.+z0m;
    zeta=zldis/obu;
    zetam=1.574;
    if (zeta < -zetam)          // zeta < -1
    {
        f10m  =  std::log(-zetam*obu/z0m) - psi(1,-zetam)
                 + psi(1,z0m/obu) + 1.14*(pow((-zeta),0.333)-pow((zetam),0.333));
    }
    else
    {
        if (zeta < 0.0)         // -1 <= zeta < 0
            f10m  = std::log(zldis/z0m) - psi(1,zeta) + psi(1,z0m/obu);
        else if (zeta <= 1.0)       //  0 <= ztea <= 1
            f10m  = std::log(zldis/z0m) + 5.*zeta - 5.*z0m/obu;
        else                           //  1 < zeta, phi=5+zeta
            f10m  = std::log(obu/z0m) + 5. - 5.*z0m/obu + (5.*std::log(zeta)+zeta-1.);
    }

    // temperature profile
    zldis=ht-displa;
    zeta=zldis/obu;
    zetat=0.465;
    if (zeta < -zetat)          // zeta < -1
    {
        fh    =   std::log(-zetat*obu/z0h)-psi(2,-zetat)
                  + psi(2,z0h/obu) + 0.8*(pow(zetat,(-0.333))-pow((-zeta),(-0.333)));
        temp1 = vonkar/fh ;
    }
    else
    {
        if (zeta < 0.0)         // -1 <= zeta < 0
        {
            fh    = std::log(zldis/z0h) - psi(2,zeta) + psi(2,z0h/obu) ;
            temp1 = vonkar/fh ;
        }
        else if (zeta <= 1.0)       //  0 <= ztea <= 1
        {
            fh    = std::log(zldis/z0h) + 5.*zeta - 5.*z0h/obu ;
            temp1 = vonkar/fh;
        }
        else 				//  1 < zeta, phi=5+zeta
        {
            fh    = std::log(obu/z0h) + 5. - 5.*z0h/obu + (5.*std::log(zeta)+zeta-1.) ;
            temp1 = vonkar/fh;
        }
    }
    // for 2 meter screen temperature
    zldis=2.+z0h ; // ht-displa
    zeta=zldis/obu;
    zetat=0.465;
    if (zeta < -zetat)         // zeta < -1
        temp12m =   vonkar/( std::log(-zetat*obu/z0h)-psi(2,-zetat)
                             + psi(2,z0h/obu) + 0.8*(pow(zetat,(-0.333))-pow((-zeta),(-0.333))));
    else
    {
        if (zeta < 0.0)         // -1 <= zeta < 0
            temp12m = vonkar/( std::log(zldis/z0h) - psi(2,zeta) + psi(2,z0h/obu) );
        else if (zeta <= 1.0)      //  0 <= ztea <= 1
            temp12m = vonkar/( std::log(zldis/z0h) + 5.*zeta - 5.*z0h/obu );
        else                           //  1 < zeta, phi=5+zeta
            temp12m = vonkar/( std::log(obu/z0h) + 5. - 5.*z0h/obu + (5.*std::log(zeta)+zeta-1.) );

    }

    // humidity profile
    zldis=hq-displa;
    zeta=zldis/obu;
    zetat=0.465;
    if (zeta < -zetat)          // zeta < -1
    {
        fq    =   std::log(-zetat*obu/z0q) - psi(2,-zetat)
                  + psi(2,z0q/obu) + 0.8*(pow(zetat,(-0.333))-pow((-zeta),(-0.333))) ;
        temp2 = vonkar/fq;
    }
    else if (zeta < 0.0)        // -1 <= zeta < 0
    {
        fq    = std::log(zldis/z0q) - psi(2,zeta) + psi(2,z0q/obu);
        temp2 = vonkar/fq;
    }
    else if (zeta <= 1.0)        //  0 <= ztea <= 1
    {
        fq    = std::log(zldis/z0q) + 5.*zeta - 5.*z0q/obu;
        temp2 = vonkar/fq;
    }
    else                            //  1 < zeta, phi=5+zeta
    {
        fq    = std::log(obu/z0q) + 5. - 5.0*z0q/obu + (5.*std::log(zeta)+zeta-1.) ;
        temp2 = vonkar/fq;
    }

    // for 2 meter screen humidity
    zldis=2.+z0h;
    zeta=zldis/obu;
    zetat=0.465;
    if (zeta < -zetat)          // zeta < -1
        temp22m =   vonkar/(std::log(-zetat*obu/z0q)-psi(2,-zetat)
                            + psi(2,z0q/obu) + 0.8*(pow(zetat,(-0.333))-pow((-zeta),(-0.333))));
    else
    {
        if (zeta < 0.0)          // -1 <= zeta < 0
            temp22m = vonkar/(std::log(zldis/z0q)-psi(2,zeta)+psi(2,z0q/obu));
        else if (zeta <= 1.)       //  0 <= zeta <= 1
            temp22m = vonkar/(std::log(zldis/z0q)+5.*zeta-5.*z0q/obu);
        else                            // 1 < zeta, phi=5+zeta
            temp22m = vonkar/(std::log(obu/z0q)+5.0-5.0*z0q/obu+(5.0*std::log(zeta)+zeta-1.0));
    }

}
void CommonLandModel::moninobukini(double ur, double th,  double thm,  double thv,double dth,
                                   double dqh,double dthv,double zldis,double z0m,double &um,double &obu)
{

// ======================================================================
// Original author : Yongjiu Dai, September 15, 1999
//
// initialzation of Monin-Obukhov length,
// the scheme is based on the work of Zeng et al. (1998):
// Intercomparison of bulk aerodynamic algorithms for the computation
// of sea surface fluxes using TOGA CORE and TAO data. J. Climate, Vol. 11: 2628-2644
// ======================================================================


// Dummy argument
    //double , INTENT(in) :: ur    // wind speed at reference height [m/s]
    //double , INTENT(in) :: thm   // intermediate variable (tm+0.0098*ht)
    //double , INTENT(in) :: th    // potential temperature [kelvin]
    //double , INTENT(in) :: thv   // virtual potential temperature (kelvin)
    //double , INTENT(in) :: dth   // diff of virtual temp. between ref. height and surface
    //double , INTENT(in) :: dthv  // diff of vir. poten. temp. between ref. height and surface
    //double , INTENT(in) :: dqh   // diff of humidity between ref. height and surface
    //double , INTENT(in) :: zldis // reference height "minus" zero displacement heght [m]
    //double , INTENT(in) :: z0m   // roughness length, momentum [m]

    //double , INTENT(out) :: um   // wind speed including the stablity effect [m/s]
    //double , INTENT(out) :: obu  // monin-obukhov length (m)

// Local
    double  wc ;    // convective velocity [m/s]
    double  rib ;   // bulk Richardson number
    double  zeta;   // dimensionless height used in Monin-Obukhov theory
    double  ustar;  // friction velocity [m/s]

//-----------------------------------------------------------------------
// Initial values of u* and convective velocity

    ustar=0.06;
    wc=0.5;
    if (dthv >= 0.0)
        um=max(ur,0.1);
    else
        um=std::sqrt(ur*ur+wc*wc);

    rib=grav*zldis*dthv/(thv*um*um);

    if (rib >= 0)     // neutral or stable
    {
        zeta = rib*std::log(zldis/z0m)/(1.0-5.0*min(rib,0.19));
        zeta = min(2.0,max(zeta,1.e-6));
    }
    else                   // unstable
    {
        zeta = rib*std::log(zldis/z0m);
        zeta = max(-100.,min(zeta,-1.e-6));
    }
    obu=zldis/zeta;

}
void CommonLandModel::netsolar (int    itypwat  ,   double sigf,      double albg[2][2],  double albv[2][2],double alb[2][2],
                                double ssun[2][2],  double ssha[2][2],double sols,        double soll,      double solsd,
                                double solld,                                                                //in
                                double &parsun,    double &parsha,      double &sabvsun,   double &sabvsha,      //out
                                double &sabg,        double &sabvg )

{
//=======================================================================
// Net solar absorbed by surface
// Original author : Yongjiu Dai, 09/15/1999; 09/11/2001
//=======================================================================
    /*
    // Dummy argument
      integer, INTENT(in) :: itypwat  // land water type (99-sea)
      real(r8), dimension(1:2,1:2), INTENT(in) :: &
            albg,   &// albedo, ground [-]
            albv,   &// albedo, vegetation [-]
            alb,    &// averaged albedo [-]
            ssun,   &// sunlit canopy absorption for solar radiation
            ssha     // shaded canopy absorption for solar radiation

      real(r8), INTENT(in) :: &
            sigf,   &// fraction of veg cover, excluding snow-buried veg [-]
            sols,   &// atm vis direct beam solar rad onto srf [W/m2]
            soll,   &// atm nir direct beam solar rad onto srf [W/m2]
            solsd,  &// atm vis diffuse solar rad onto srf [W/m2]
            solld    // atm nir diffuse solar rad onto srf [W/m2]

      real(r8), INTENT(out) :: &
            parsun, &// PAR absorbed by sunlit vegetation [W/m2]
            parsha, &// PAR absorbed by shaded vegetation [W/m2]
            sabvsun,&// solar absorbed by sunlit vegetation [W/m2]
            sabvsha,&// solar absorbed by shaded vegetation [W/m2]
            sabg,   &// solar absorbed by ground  [W/m2]
            sabvg    // solar absorbed by ground + vegetation [W/m2]

    //=======================================================================
    */
    sabvsun = 0.0;
    sabvsha = 0.0;
    parsun = 0.0;
    parsha = 0.0;

    sabg = 0.0;
    sabvg = 0.0;

    if ((sols+soll+solsd+solld)>0.0)
    {
        if (itypwat<4)       //non lake and ocean
        {
            // Radiative fluxes onto surface
            parsun  = ssun[0][0]*sols + ssun[0][1]*solsd;
            parsha  = ssha[0][0]*sols + ssha[0][1]*solsd;
            sabvsun = sols*ssun[0][0] + soll*ssun[1][0]
                      + solsd*ssun[0][1] + solld*ssun[1][1];
            sabvsha = sols*ssha[0][0] + soll*ssha[1][0]
                      + solsd*ssha[0][1] + solld*ssha[1][1];
            sabvsun = sigf*sabvsun;
            sabvsha = sigf*sabvsha;
            sabvg = sols *(1.-alb[0][0]) + soll *(1.-alb[1][0])
                    + solsd*(1.-alb[0][1]) + solld*(1.-alb[1][1]);
            sabg = sabvg - sabvsun - sabvsha;
        }
        else                     //lake or ocean
        {
            sabvg = sols *(1.-alb[0][0]) + soll *(1.-alb[1][0])
                    + solsd*(1.-alb[0][1]) + solld*(1.-alb[1][1]);
            sabg = sabvg;
        }
    }
}
void CommonLandModel::newsnow (int    itypwat,                int   maxsnl,  double dtime,  double tm,double tg,
                               double pg,                     double tcrit,                                                   //in
                               double zi[MAXSNL+MAXSOILL+1],  double z[MAXSNL+MAXSOILL],     double dz[MAXSNL+MAXSOILL],
                               double tss[MAXSNL+MAXSOILL],   double wliq[MAXSNL+MAXSOILL],  double wice[MAXSNL+MAXSOILL],
                               double fiold[MAXSNL+MAXSOILL], int    &snl,                     double &sag,                     //inout
                               double &scv,			         double &snowdp,
                               double &pg_rain,               double &pg_snow)                                                                 //out
{
    // ------------------------ Dummy Argument ------------------------------

    //integer, INTENT(in) :: maxsnl  // maximum number of snow layers
    //integer, INTENT(in) :: itypwat // land water type (0=soil, 1=urban and built-up,
    //                               // 2=wetland, 3=land ice, 4=deep lake, 5=shallow lake)
    //double, INTENT(in) :: dtime  // model time step [second]
    //double, INTENT(in) :: tm     // temperature at agcm reference height [kelvin]
    //double, INTENT(in) :: tg     // ground surface temperature [k]
    //double, INTENT(in) :: pg     // water onto ground including canopy runoff [kg/(m2 s)]
    //double, INTENT(in) :: tcrit  // critical temp. to determine rain or snow

    //double, INTENT(inout) ::    zi(maxsnl:0)   // interface level below a "z" level (m)
    //double, INTENT(inout) ::     z(maxsnl+1:0) // layer depth (m)
    //double, INTENT(inout) ::    dz(maxsnl+1:0) // layer thickness (m)
    //double, INTENT(inout) ::   tss(maxsnl+1:0) // soil + snow layer temperature [K]
    //double, INTENT(inout) ::  wliq(maxsnl+1:0) // liquid water (kg/m2)
    //double, INTENT(inout) ::  wice(maxsnl+1:0) // ice lens (kg/m2)
    //double, INTENT(inout) :: fiold(maxsnl+1:0) // fraction of ice relative to the total water
    // integer, INTENT(inout) :: snl               // number of snow layers
    //double, INTENT(inout) :: sag               // non dimensional snow age [-]
    //double, INTENT(inout) :: scv               // snow mass (kg/m2)
    //double, INTENT(inout) :: snowdp            // snow depth (m)

    //double, INTENT(out) :: pg_rain  // liquid water onto ground [kg/(m2 s)]
    //double, INTENT(out) :: pg_snow  // ice onto ground [kg/(m2 s)]

    // ----------------------- Local  Variables -----------------------------

    double bifall;     // bulk density of newly fallen dry snow [kg/m3]
    double flfall;     // fraction of liquid water within falling precip.
    double dz_snowf ;  // layer thickness rate change due to precipitation [m/s]
    int  newnode;     // signification when new snow node is set, (1=yes, 0=non)

    //-----------------------------------------------------------------------
    // the upper limit of air temperature is set for snowfall, this cut-off
    // was selected based on Fig. 1, Plate 3-1, of Snow Hydrology (1956).
    // the percentage of liquid water by mass, which is arbitrarily set to
    // vary linearly with air temp, from 0% at 273.16 to 40% max at 275.16.

    double scv_old = scv;
    newnode = 0;

    //cout<<snowdp*1000<<endl;

    if (tm>(tfrz+tcrit))
    {
        flfall = 1.0;        // fraction of liquid water within falling precip.
        pg_snow = 0.0;       // ice onto ground (mm/s)
        pg_rain = pg;       // liquid water onto ground (mm/s)
        dz_snowf = 0.0;      // rate of snowfall, snow depth/s (m/s)
    }
    else
    {
        if (tm<=tfrz)
            flfall=0.0;
        else
        {
            if (tm<=(tfrz+2.0))
                flfall=-54.632+0.2*tm;
            else
                flfall=0.4;
        }


        // use Alta relationship, Anderson(1976); LaChapelle(1961),
        // U.S.Department of Agriculture Forest Service, Project F,
        // Progress Rep. 1, Alta Avalanche Study Center:Snow Layer Densification.

        if (tm>(tfrz+2.0))
            bifall=169.0;
        else
        {
            if (tm>(tfrz-15.0))
                bifall=50.+1.7*pow((tm-tfrz+15.0),1.5);
            else
                bifall=50.0;
        }

        pg_snow = pg*(1.0-flfall);


        pg_rain = pg*flfall;

        dz_snowf = pg_snow/bifall;
        snowdp = snowdp + dz_snowf*dtime;
        scv = scv + pg_snow*dtime;      // snow water equivalent (mm)

        //cout<<snowdp*1000<<"  "<<dz_snowf*dtime*1000<<endl;


        if ((itypwat==2)&&(tg>tfrz)) // snowfall on warmer wetland
        {
            scv=0.0;
            snowdp=0.0;
            sag=0.0;
        }
    }


    // when the snow accumulation exceeds 10 mm, initialize a snow layer
    if ((snl==maxsnl-1)&&(pg_snow>0.0)&&(snowdp>=0.01))
    {
        snl = snl -1;
        newnode = 1;
        dz[maxsnl-1]  = snowdp;             // meter
        z[maxsnl-1]   = -0.5*dz[maxsnl-1];
        zi[maxsnl-1]= -dz[maxsnl-1];
        // currently, the water temperature for the precipitation is simply set
        // as the surface air temperature
        sag = 0.0;                    // snow age
        tss[maxsnl-1]  = min(tfrz, tm);     // K
        wice[maxsnl-1] = scv;               // kg/m2
        wliq[maxsnl-1] = 0.0;                // kg/m2
        fiold[maxsnl-1] = 1.0;
        //**      write(6,*) 'snow layer is built'
    }

    // the change of ice partial density of surface node due to precipitation
    // only ice part of snowfall is added here, the liquid part will be added latter
    if ((snl<maxsnl-1)&&(newnode==0))
    {
        wice[snl+1] = wice[snl+1] +dtime*pg_snow;
        dz[snl+1] = dz[snl+1] +dz_snowf*dtime;
        z[snl+1] = zi[snl+2] - 0.5*dz[snl+1];
        zi[snl+1] = zi[snl+2] - dz[snl+1];
    }
}
double CommonLandModel::orb_coszen(double calday,double lon,double lat)
{

//-------------------------------------------------------------------------------
// FUNCTION to return the cosine of the solar zenith angle. Assumes 365.0 days/year.
// Compute earth/orbit parameters using formula suggested by
// Duane Thresher. Use formulas from Berger, Andre 1978: Long-Term Variations of Daily
// Insolation and Quaternary Climatic Changes. J. of the Atmo. Sci. 35:2362-2367.
//
// Original version:  Erik Kluzek, Oct/1997, Brian Kauffman, Jan/98
// CCSM2.0 standard
// yongjiu dai (07/23/2002)
//-------------------------------------------------------------------------------

//use precision
//implicit none

//real(r8), intent(in) :: calday        //Julian cal day (1.xx to 365.xx)
//real(r8), intent(in) :: lat           //Centered latitude (radians)
//real(r8), intent(in) :: lon           //Centered longitude (radians)
//real(r8) :: orb_coszen

// --- Local variables ---
    double   declin;                       //Solar declination (radians)
    double   eccf ;                        //Earth-sun distance factor (ie. (1/r)**2)
    double   lambm ;                       //Lambda m, mean long of perihelion (rad)
    double   lmm  ;                        //Intermediate argument involving lambm
    double   lamb ;                        //Lambda, the earths long of perihelion
    double   invrho ;                      //Inverse normalized sun/earth distance
    double   sinl ;                        //Sine of lmm
    double   dayspy=365.0;                 //days per year
    double   ve=80.5;                      //Calday of vernal equinox assumes Jan 1 = calday 1
    double   eccen = 1.672393084E-2;       //Eccentricity
    double   obliqr= 0.409214646;          //Earths obliquity in radians
    double   lambm0= -3.2625366E-2;        //Mean long of perihelion at the vernal equinox (radians)
    double   mvelpp= 4.92251015;           //moving vernal equinox longitude of
    //perihelion plus pi (radians)
    double   orb_coszen;
//-------------------------------------------------------------------------------

    lambm = lambm0 + (calday - ve)*2.*PI/dayspy;
    lmm = lambm  - mvelpp;

    sinl = std::sin(lmm);
    lamb = lambm + eccen*(2.*sinl + eccen*(1.25*std::sin(2.*lmm)
                                           + eccen*((13.0/12.0)*std::sin(3.*lmm) - 0.25*sinl)));
    invrho = (1. + eccen*std::cos(lamb - mvelpp)) / (1. - eccen*eccen);

    declin = std::asin(std::sin(obliqr)*std::sin(lamb));
    eccf = invrho*invrho;

    orb_coszen = std::sin(lat)*std::sin(declin)
                 - std::cos(lat)*std::cos(declin)*std::cos(calday*2.0*PI+lon);

    return orb_coszen;
}
double  CommonLandModel::psi(int k,double zeta)
{
    //=======================================================================
    // stability function for rib < 0

    // use precision
    // implicit none

    //long    k;
    double  psi1;   // stability function for unstable case
    //double  zeta;  // dimensionless height used in Monin-Obukhov theory
    double  chik;  //

    chik = pow((1.-16.*zeta),0.25);
    if (k == 1)
        psi1 = 2.0*std::log((1.+chik)*0.5)+std::log((1.+chik*chik)*0.5)-2.0*std::atan(chik)+2.0*std::atan(1.0);
    else
        psi1 = 2.*std::log((1.0+chik*chik)*0.5);
    return psi1;
}
void CommonLandModel::qsadv(double T,double p,double &es,double &esdT,double &qs,double &qsdT)
{

//=======================================================================
// Original author : Yongjiu Dai, September 15, 1999
//
//      Description: computes saturation mixing ratio and change in saturation
//                   mixing ratio with respect to temperature
//
//        Reference: polynomial approximations from:
//                   Piotr J. Flatau,et al,1992: polynomial fits to saturation
//                   vapor pressure. Journal of Applied meteorology,31,1507-1513.
//
//-----------------------------------------------------------------------
//  use precision
// implicit none

// dummy arguments
//  double, INTENT(in)  :: T        // temperature (K)
//  double, INTENT(in)  :: p        // surface atmospheric pressure (pa)
//
//  double, INTENT(out) :: es       // vapor pressure (pa)
//  double, INTENT(out) :: esdT     // d(es)/d(T)
//  double, INTENT(out) :: qs       // humidity (kg/kg)
//  double, INTENT(out) :: qsdT     // d(qs)/d(T)
//

//-----------------------------------------------------------------------
// local
    double td,vp,vp1,vp2;
    double a0,a1,a2,a3,a4,a5,a6,a7,a8;
    double b0,b1,b2,b3,b4,b5,b6,b7,b8;

    double c0,c1,c2,c3,c4,c5,c6,c7,c8;
    double d0,d1,d2,d3,d4,d5,d6,d7,d8;
    double T_limit;

// for water vapor (temperature range 0C-100C)
    a0 = 6.11213476    ;
    a1 =  0.444007856   ;
    a2 = 0.143064234e-01;
    a3 = 0.264461437e-03 ;
    a4 =  0.305903558e-05 ;
    a5 = 0.196237241e-07 ;
    a6 = 0.892344772e-10 ;
    a7 = -0.373208410e-12 ;
    a8 = 0.209339997e-15;

// for derivative:water vapor
    b0 = 0.444017302   ;
    b1 =  0.286064092e-01 ;
    b2 =  0.794683137e-03 ;
    b3 =  0.121211669e-04 ;
    b4 =  0.103354611e-06 ;
    b5 =  0.404125005e-09 ;
    b6 = -0.788037859e-12 ;
    b7 = -0.114596802e-13 ;
    b8 =  0.381294516e-16 ;

// for ice (temperature range -75C-0C)
    c0 = 6.11123516      ;
    c1 = 0.503109514     ;
    c2 = 0.188369801e-01 ;
    c3 = 0.420547422e-03 ;
    c4 = 0.614396778e-05 ;
    c5 = 0.602780717e-07 ;
    c6 = 0.387940929e-09 ;
    c7 = 0.149436277e-11 ;
    c8 = 0.262655803e-14 ;

//  for derivative:ice
    d0 = 0.503277922     ;
    d1 = 0.377289173e-01 ;
    d2 = 0.126801703e-02 ;
    d3 = 0.249468427e-04 ;
    d4 = 0.313703411e-06 ;
    d5 = 0.257180651e-08 ;
    d6 = 0.133268878e-10 ;
    d7 = 0.394116744e-13 ;
    d8 = 0.498070196e-16 ;

//=======================================================================
//根据fortran版本的更新(2008)
    T_limit     = T - 273.15; //273.16
    if (T_limit>100.0)
        T_limit = 100.0;
    if (T_limit<-75.0)
        T_limit = -75.0;
    td = T_limit;

    if (td >= 0.0)
    {
        es   = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 + td*(a5 + td*(a6 + td*(a7 + td*a8)))))));
        esdT = b0 + td*(b1 + td*(b2 + td*(b3 + td*(b4 + td*(b5 + td*(b6 + td*(b7 + td*b8)))))));
    }
    else
    {
        es   = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 + td*(c5 + td*(c6 + td*(c7 + td*c8)))))));
        esdT = d0 + td*(d1 + td*(d2 + td*(d3 + td*(d4 + td*(d5 + td*(d6 + td*(d7 + td*d8)))))));
    }
    es    = es    * 100.0;            // pa
    esdT  = esdT  * 100.0;            // pa/K

    vp    = 1.0   / (p - 0.378*es);
    vp1   = 0.622 * vp;
    vp2   = vp1   * vp;

    qs    = es    * vp1;             // kg/kg
    qsdT  = esdT  * vp2 * p;         // 1 / K

}

// ****************************************************************************
// 计算可见和近红波段的直接和散射辐射
// according to SiB2 model code
// ****************************************************************************
// Author: Chunlin Huang
// Date    August 29, 2006


// input:
// calday        //Julian cal day (1.xx to 365.xx)
// lat           //Centered latitude (radians)
// lon           //Centered longitude (radians)
//  swdown          downward short wave radiation

//output:
//double   sols     ; // atm vis direct beam solar rad onto srf [W/m2]
//double   soll     ; // atm nir direct beam solar rad onto srf [W/m2]
//double   solsd    ; // atm vis diffuse solar rad onto srf [W/m2]
//double   solld    ; // atm nir diffuse solar rad onto srf [W/m2]

void CommonLandModel::radiation(int jday,int msec,double lon,double lat, double swdown, double &sols, double &soll, double &solsd, double &solld)
{
    double sunang, cloud;

    double 	calday = (double)jday + (double)(msec/86400.0);

    sunang = orb_coszen(calday,lon,lat);

    cloud = (1160.*sunang - swdown) / (963. * sunang);
    cloud = max(cloud,0.);
    cloud = min(cloud,1.);
    cloud = max(0.58,cloud);

    double difrat = 0.0604 / ( sunang-0.0223 ) + 0.0683;
    if ( difrat < 0. ) difrat = 0.0;
    if ( difrat > 1. ) difrat = 1.0;

    difrat = difrat + ( 1. - difrat ) * cloud;
    double vnrat = ( 580. - cloud*464. ) / ( ( 580. - cloud*499. )
                   + ( 580. - cloud*464. ) );

    sols  = (1.-difrat)*vnrat*swdown;
    soll  = (1.-difrat)*(1.-vnrat)*swdown;
    solsd = difrat*vnrat*swdown;
    solld = difrat*(1.-vnrat)*swdown;

    //cout<<swdown<<"  "<<sols<<"  "<<soll<<"  "<<solsd<<"  "<<solld<<endl;


}
void CommonLandModel::seafluxes (double   oro,    double   hu,    double   ht,    double   hq,    double   us,
                                 double   vs,     double   tm,    double   qm,    double   rhoair,double   psrf,
                                 double   tssea,  double   &taux,  double  & tauy,  double  & fsena, double   &fevpa,
                                 double   &fseng,  double   &fevpg, double   &tref,  double   &qref,  double   &z0ma,
                                 double   &zol,    double   &rib,   double  &ustar,  double   &qstar, double   &tstar,
                                 double   &u10m,   double   &v10m,  double   &f10m,  double   &fm,    double   &fh,
                                 double   &fq,     double   &cgrndl,double   &cgrnds)
{

//=======================================================================
// this is the main subroutine to execute the calculation of thermal processes
// and surface fluxes
//
// Original author : Yongjiu Dai, 09/15/1999; 08/30/2002
//=======================================================================

    //use precision
    //use phycon_module, only : cpair,rgas,vonkar,grav
    //implicit none

//----------------------- Dummy argument --------------------------------

    //real(r8), INTENT(in) :: &
    //      oro,      &// ocean(0)/seaice(2)/ flag

    //      // atmospherical variables and agcm reference height
    //      hu,       &// agcm reference height of wind [m]
    //      ht,       &// agcm reference height of temperature [m]
    //      hq,       &// agcm reference height of humidity [m]
    //      us,       &// wind component in eastward direction [m/s]
    //      vs,       &// wind component in northward direction [m/s]
    //      tm,       &// temperature at agcm reference height [kelvin]
    //      qm,       &// specific humidity at agcm reference height [kg/kg]
    //      rhoair,   &// density air [kg/m3]
    //      psrf,     &// atmosphere pressure at the surface [pa] [not used]

    //      tssea      // sea surface temperature [K]

    //real(r8), INTENT(out) :: &
    //      taux,     &// wind stress: E-W [kg/m/s**2]
    //      tauy,     &// wind stress: N-S [kg/m/s**2]
    //      fsena,    &// sensible heat from agcm reference height to atmosphere [W/m2]
    //      fevpa,    &// evaporation from agcm reference height to atmosphere [mm/s]
    //      fseng,    &// sensible heat flux from ground [W/m2]
    //      fevpg,    &// evaporation heat flux from ground [mm/s]

    //      tref,     &// 2 m height air temperature [kelvin]
    //      qref,     &// 2 m height air humidity
    //      z0ma,     &// effective roughness [m]
    //      zol,      &// dimensionless height (z/L) used in Monin-Obukhov theory
    //      rib,      &// bulk Richardson number in surface layer
    //      ustar,    &// friction velocity [m/s]
    //      tstar,    &// temperature scaling parameter
    //      qstar,    &// moisture scaling parameter
    //      u10m,     &// 10m u-velocity
    //      v10m,     &// 10m v-velocity
    //      f10m,     &// integral of profile function for momentum at 10m
    //      fm,       &// integral of profile function for momentum
    //      fh,       &// integral of profile function for heat
    //      fq,       &// integral of profile function for moisture
    //      cgrndl,   &// deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
    //      cgrnds     // deriv of soil latent heat flux wrt soil temp [w/m**2/k]

//------------------------ LOCAL VARIABLES ------------------------------
    int i;
    int niters, // maximum number of iterations for surface temperature
        iter,      // iteration index
        nmozsgn;     // number of times moz changes sign

    double
    beta,      // coefficient of conective velocity [-]
    displax,   // zero-displacement height [m]
    dth,       // diff of virtual temp. between ref. height and surface
    dqh,       // diff of humidity between ref. height and surface
    dthv,      // diff of vir. poten. temp. between ref. height and surface
    eg,        // water vapor pressure at temperature T [Pa]
    degdT,     // d(eg)/dT
    obu,       // monin-obukhov length (m)
    obuold,    // monin-obukhov length from previous iteration
    qsatg,     // ground saturated specific humidity [kg/kg]
    qsatgdT,   // d(qsatg)/dT
    ram,       // aerodynamical resistance [s/m]
    rah,       // thermal resistance [s/m]
    raw,       // moisture resistance [s/m]
    raih,      // temporary variable [kg/m2/s]
    raiw,      // temporary variable [kg/m2/s]
    temp1,     // relation for potential temperature profile
    temp2,     // relation for specific humidity profile
    temp12m,   // relation for temperature at 2m
    temp22m,   // relation for specific humidity at 2m
    thm,       // intermediate variable (tm+0.0098*ht)
    th,        // potential temperature (kelvin)
    thv,       // virtual potential temperature (kelvin)
    thvstar,   // virtual potential temperature scaling parameter
    um,        // wind speed including the stablity effect [m/s]
    ur,        // wind speed at reference height [m/s]
    visa,      // kinematic viscosity of dry air [m2/s]
    wc,        // convective velocity [m/s]
    wc2,       // wc**2
    xt,        //
    xq,        //
    zii,       // convective boundary height [m]
    zldis,     // reference height "minus" zero displacement heght [m]
    z0mg,      // roughness length over ground, momentum [m]
    z0hg,      // roughness length over ground, sensible heat [m]
    z0qg;        // roughness length over ground, latent heat [m]

    const double zsice = 0.04;  // sea ice aerodynamic roughness length [m]

//-----------------------------------------------------------------------
// potential temperatur at the reference height
    beta = 1.0;      // -  (in computing W_*)
    zii = 1000.0 ;   // m  (pbl height)

//-----------------------------------------------------------------------
//     Compute sensible and latent fluxes and their derivatives with respect
//     to ground temperature using ground temperatures from previous time step.
//-----------------------------------------------------------------------
// Initialization variables
    nmozsgn = 0;
    obuold = 0.0;

    qsadv(tssea,psrf,eg,degdT,qsatg,qsatgdT);

// potential temperatur at the reference height
    thm = tm + 0.0098*ht ;             // intermediate variable equivalent to
    // tm*(pgcm/psrf)**(rgas/cpair)
    th = tm*pow((100000./psrf),(rgas/cpair)); // potential T
    thv = th*(1.+0.61*qm);             // virtual potential T
    ur = max(0.1,std::sqrt(us*us+vs*vs));   // limit set to 0.1

    dth   = thm-tssea;
    dqh   = qm-qsatg;
    dthv  = dth*(1.+0.61*qm)+0.61*th*dqh;
    zldis = hu-0.;

    if ((int)oro==0)         // ocean
        // Kinematic viscosity of dry air (m2/s)- Andreas (1989) CRREL Rep. 89-11
    {
        visa=1.326e-5*(1.+6.542e-3*tm + 8.301e-6*tm*tm - 4.84e-9*tm*tm*tm);

        // loop to obtain initial and good ustar and zo
        ustar=0.06;
        wc=0.5;
        if (dthv>=0.0)
            um=max(ur,0.1);
        else
            um=std::sqrt(ur*ur+wc*wc);

        for (i=1; i<=5; i++)
        {
            z0mg=0.013*ustar*ustar/grav+0.11*visa/ustar;
            ustar=vonkar*um/std::log(zldis/z0mg);
        }
    }

    else if ((int)oro==2)   // sea ice
    {
        z0mg = zsice;
        z0qg = z0mg;
        z0hg = z0mg;
    }

    moninobukini(ur,th,thm,thv,dth,dqh,dthv,zldis,z0mg,um,obu);

// Evaluated stability-dependent variables using moz from prior iteration
    niters=10;
    displax = 0.0;

    //----------------------------------------------------------------
    // do iter = 1, niters         // begin stability iteration
    //----------------------------------------------------------------

    for (iter=1; iter<=niters; iter++)
    {
        if ((int)oro==0)  // ocean
        {
            z0mg=0.013*ustar*ustar/grav + 0.11*visa/ustar;
            xq=2.67*pow((ustar*z0mg/visa),0.25) - 2.57;
            xt= xq;
            z0qg=z0mg/std::exp(xq);
            z0hg=z0mg/std::exp(xt);
        }
        moninobuk(hu,ht,hq,displax,z0mg,z0hg,z0qg,obu,um,
                  ustar,temp1,temp2,temp12m,temp22m,f10m,fm,fh,fq);

        tstar = temp1*dth;
        qstar = temp2*dqh;

        thvstar=tstar+0.61*th*qstar;
        zol=zldis*vonkar*grav*thvstar/(ustar*ustar*thv);
        if (zol >= 0.0)      // stable
            zol = min(2.,max(zol,1.e-6));
        else                     // unstable
            zol = max(-100.,min(zol,-1.e-6));

        obu = zldis/zol;

        if (zol >= 0.)
            um = max(ur,0.1);
        else
        {
            wc = pow((-grav*ustar*thvstar*zii/thv),(1./3.));
            wc2 = beta*beta*(wc*wc);
            um = std::sqrt(ur*ur+wc2);
        }
        if ((obuold*obu)< 0.0)
            nmozsgn = nmozsgn+1;
        if (nmozsgn >= 4)
            break;

        obuold = obu;
    }

    //----------------------------------------------------------------
    //enddo                       // end stability iteration
    //----------------------------------------------------------------

// Get derivative of fluxes with repect to ground temperature
    ram    = 1./(ustar*ustar/um);
    rah    = 1./(temp1*ustar) ;
    raw    = 1./(temp2*ustar) ;

    raih   = rhoair*cpair/rah;
    raiw   = rhoair/raw;
    cgrnds = raih;
    cgrndl = raiw*qsatgdT;

    rib = min(5.,zol*ustar*ustar/(vonkar*temp1*um*um));

// surface fluxes of momentum, sensible and latent
// using ground temperatures from previous time step
    taux   = -rhoair*us/ram ;
    tauy   = -rhoair*vs/ram;

    fseng  = -raih*dth;
    fevpg  = -raiw*dqh ;
    fsena  = fseng;
    fevpa  = fevpg;

// 2 m height air temperature
    tref   = thm + temp1*dth * (1./temp12m - 1./temp1);
    qref   = qm + temp2*dqh * (1./temp22m - 1./temp2);
    z0ma   = z0mg;

// 10 m wind
    u10m = us/max(0.1,ur) * ustar/vonkar * f10m;
    v10m = vs/max(0.1,ur) * ustar/vonkar * f10m;
}
void CommonLandModel::snowage(double dtime,double tg,double scv,double scvold,double &sag )
{

    //=======================================================================
    // Original version: Robert Dickinson
    // Update snow cover and snow age, based on BATS code
    //=======================================================================

    // use precision
    // use phycon_module, only : tfrz
    // implicit none

    //-------------------------- Dummy Argument -----------------------------

    // real(r8), INTENT(in) :: dtime  // time step [second]
    // real(r8), INTENT(in) :: tg     // temperature of soil at surface [K]
    // real(r8), INTENT(in) :: scv    // snow cover, water equivalent [mm]
    // real(r8), INTENT(in) :: scvold // snow cover for previous time step [mm]
    // real(r8), INTENT(inout) :: sag // non dimensional snow age [-]

    //-------------------------- Local variables ----------------------------

    double  age1;   // snow aging factor due to crystal growth [-]
    double  age2;   // snow aging factor due to surface growth [-]
    double  age3;   // snow aging factor due to accum of other particles [-]
    double  arg;    // temporary variable used in snow age calculation [-]
    double  arg2;   // temporary variable used in snow age calculation [-]
    double  dela;   // temporary variable used in snow age calculation [-]
    double  dels ;  // temporary variable used in snow age calculation [-]
    double  sge;    // temporary variable used in snow age calculation [-]

    //-----------------------------------------------------------------------
    if (scv <= 0.0)
        sag = 0.0;
    //
    // Over antarctica
    //
    else if (scv > 800.0)
        sag = 0.0;
    //
    // Away from antarctica
    //
    else
    {
        age3  = 0.3;
        arg   = 5.e3*(1./tfrz-1./tg);
        arg2  = min(0.,10.*arg);
        age2  = std::exp(arg2);
        age1  = std::exp(arg);
        dela  = 1.e-6*dtime*(age1+age2+age3);
        dels  = 0.1*max(0.0,scv-scvold);
        sge   = (sag+dela)*(1.0-dels);
        sag   = max(0.0,sge);
    }

}
void CommonLandModel::snowcompaction(int    lb,                   double dtime,               int    imelt[MAXSNL+MAXSOILL], double fiold[MAXSNL+MAXSOILL], double tss[MAXSNL+MAXSOILL],
                                     double wliq[MAXSNL+MAXSOILL], double wice[MAXSNL+MAXSOILL],double dz[MAXSNL+MAXSOILL] )

{
    //=======================================================================
    // Original author: Yongjiu Dai, September 15, 1999
    //
    // three of metamorphisms of changing snow characteristics are implemented,
    // i.e., destructive, overburden, and melt. The treatments of the former two
    // are from SNTHERM.89 and SNTHERM.99 (1991, 1999). The contribution due to
    // melt metamorphism is simply taken as a ratio of snow ice fraction after
    // the melting versus before the melting.
    //
    //=======================================================================

    // use precision
    // use phycon_module, only : denice, denh2o, tfrz
    // implicit none
    //
    //-------------------------- Dummy argument -----------------------------
    //
    //integer, INTENT(in) :: &
    //      lb             // lower bound of array

    //integer, INTENT(in) :: &
    //      imelt(lb : 0)  // signifies if node in melting (imelt = 1)

    //real(r8), INTENT(in) :: &
    //      dtime,        &// time step [second]
    //      fiold(lb : 0),&// fraction of ice relative to
    //                     // the total water content at the previous time step
    //      tss (lb : 0), &// nodal temperature [K]
    //      wice(lb : 0), &// ice lens [kg/m2]
    //      wliq(lb : 0)   // liquid water [kg/m2]

    //real(r8), INTENT(inout) :: &
    //      dz  (lb : 0)   // layer thickness [m]
    //
    //----------------------- local variables ------------------------------
    int i;            // Numeber of doing loop

    double
//		c1,            // = 2.777e-7 [m2/(kg s)]
    c2,            // = 21e-3 [m3/kg]
    c3,            // = 2.777e-6 [1/s]
    c4,            // = 0.04 [1/K]
    c5,            // = 2.0
    c6,            // = 5.15e-7.
    c7,            // = 4.
    dm,            // Upper Limit on Destructive Metamorphism Compaction [kg/m3]
    eta0;            // The Viscosity Coefficient Eta0 [kg-s/m2]

    double
    burden,        // pressure of overlying snow [kg/m2]
    ddz1,          // Rate of settling of snowpack due to destructive metamorphism.
    ddz2,          // Rate of compaction of snowpack due to overburden.
    ddz3,          // Rate of compaction of snowpack due to melt [1/s]
    dexpf,         // expf=std::exp(-c4*(273.15-tss)).
    fi,            // Fraction of ice relative to the total water content
    // at the current time step
    td,            // tss - tfrz [K]
    pdzdtc,        // Nodal rate of change in fractional-thickness
    // due to compaction [fraction/s]
    vol_id,          // vol_id (1 - vol_ice - vol_liq)
    wx,            // water mass (ice+liquid) [kg/m2]
    bi;              // partitial density of ice [kg/m3]

    c2 =  23.e-3;
    c3 =  2.777e-6;
    c4 =  0.04;
    c5 =  2.0;
    c6 =  5.15e-7;
    c7 =  4.0;
    dm =  100.0;
    eta0 = 9.0e+5;

    //=======================================================================

    burden = 0.0;

    for (i=lb; i<MAXSNL; i++)
    {
        wx = wice[i] + wliq[i];
        vol_id = 1.0- (wice[i]/denice + wliq[i]/denh2o)/dz[i];
        // cout<<i<<"  "<<wice[i]<<" "<<wliq[i]<<" "<<dz[i]<<endl;
        //cout<<i<<" "<<wx<<" "<<vol_id<<endl;

        // Disallow compaction for water saturated node and lower ice lens node.
        if ((vol_id <= 0.001)||(wice[i] <= 0.1))
        {
            burden = burden+wx;
            continue;
        }

        bi = wice[i] / dz[i];
        fi = wice[i] / wx;
        td = tfrz-tss[i];

        dexpf = std::exp(-c4*td);

        // Settling as a result of destructive metamorphism
        ddz1 = -c3*dexpf;

        if (bi > dm)
            ddz1 = ddz1*std::exp(-46.0e-3*(bi-dm));

        // Liquid water term
        if (wliq[i] > (0.01*dz[i]))
            ddz1=ddz1*c5;

        ddz2 = -burden*std::exp(-0.08*td - c2*bi)/eta0;

        // Compaction occurring during melt
        if (imelt[i] == 1)
            ddz3 = - 1./dtime * max(0.,(fiold[i] - fi)/fiold[i]);
        else
            ddz3 = 0.0;
//		cout<<i<<" "<<ddz1<<"  "<<ddz2<<"  "<<ddz3<<endl;

        // Time rate of fractional change in dz (units of s-1)
        pdzdtc = ddz1+ddz2+ddz3;
        //cout<<i<<"  "<<pdzdtc<<endl;

        // The change in dz due to compaction
        dz[i] = dz[i]*(1.+pdzdtc*dtime);
        //cout<<dz[i]<<endl;

        // Pressure of overlying snow
        burden = burden+wx;
    }
}
void CommonLandModel::snowfraction (double fveg,double z0m,double snowdp,double &wt,double &sigf,double &fsno)
{
//=======================================================================
// Original author : Yongjiu Dai, September 15, 1999
// Provide snow cover fraction
//=======================================================================

//  use precision
//  implicit none

// dummy arguments
//  real(r8), INTENT(in) :: snowdp // snow depth [m]
//  real(r8), INTENT(in) :: z0m    // aerodynamic roughness length [m]
//  real(r8), INTENT(in) :: fveg   // fractional vegetation cover [-]

//  real(r8), INTENT(out) :: wt    // fraction of vegetation covered with snow [-]
//  real(r8), INTENT(out) :: sigf  // fraction of veg cover, excluding snow-covered veg [-]
//  real(r8), INTENT(out) :: fsno  // fraction of soil covered by snow [-]

//-----------------------------------------------------------------------
    if (fveg > 0.001)
    {

// Fraction of vegetation buried (covered) by snow
        wt = 0.1*snowdp/z0m;
        wt = wt/(1.+wt);

// Fraction of vegetation cover free of snow
        sigf = (1.-wt)*fveg;

// Fraction of soil covered by snow
        fsno = snowdp/(0.1+snowdp);
    }
    else
    {
        wt = 0.0;
        sigf = 0.0;
        fsno = snowdp/(0.1+snowdp);
    }
    if (sigf < 0.001)
        sigf = 0.0;
    if (sigf > 0.999)
        sigf = 1.0;

}
void CommonLandModel::snowlayerscombine(int lb,                      int     &snl,                  double z[MAXSNL+MAXSOILL],  double  dz[MAXSNL+MAXSOILL],double zi[MAXSNL+MAXSOILL+1],
                                        double wliq[MAXSNL+MAXSOILL],double  wice[MAXSNL+MAXSOILL],double tss[MAXSNL+MAXSOILL],double &scv,
                                        double &snowdp )
{

    //=======================================================================
    // Original author : Yongjiu Dai, September 15, 1999
    //
    // checks for elements which are below prescribed minimum for thickness or mass.
    // If snow element thickness or mass is less than a prescribed minimum,
    // it is combined with neighboring element to be best combine with,
    // and executes the combination of mass and energy in clm_combo.f90
    //
    //=======================================================================

    //use precision
    //implicit none

    //-------------------------- Dummy argument -----------------------------
    //integer, INTENT(in) :: lb               // lower bound of array

    // numbering from 1 (bottom) mss (surface)
    //real(r8), INTENT(inout) :: wice(lb:1)   // ice lens [kg/m2]
    //real(r8), INTENT(inout) :: wliq(lb:1)   // liquid water {kg/m2]
    //real(r8), INTENT(inout) :: tss (lb:1)   // nodel temperature [K]
    //real(r8), INTENT(inout) :: dz  (lb:1)   // layer thickness [m]
    //real(r8), INTENT(inout) :: z   (lb:1)   // node depth [m]
    //real(r8), INTENT(inout) :: zi  (lb-1:1) // depth of layer interface [m]
    //real(r8), INTENT(inout) :: snowdp       // snow depth [m]
    //real(r8), INTENT(inout) :: scv          // snow mass - water equivalent [kg/m2]
    //integer, INTENT(inout) :: snl           // Number of snow

    //----------------------- Local variables ------------------------------
//	double  drr;           // thickness of the combined [m]
    double  dzmin[5];      // minimum of snow layer 1 (top) to msn0 (bottom)
    double  zwice;         // total ice mass in snow
    double  zwliq;         // total liquid water in snow

    int  i       ;       // number of do looping
    int  j       ;       // node index
    int  k       ;       // number of do looping
    int  l       ;       // node index
    int  msn_old ;       // number of snow layer 1 (top) to msn0 (bottom)
    int  mssi   ;        // node index
    int  neibor  ;       // adjacent node selected for combination

    dzmin[0] = 0.010;
    dzmin[1] = 0.015;
    dzmin[2] = 0.025;
    dzmin[3] = 0.055;
    dzmin[4] = 0.115;

    //-----------------------------------------------------------------------
    // check the mass of ice lens of snow, when the total less than a small value,
    // combine it with the underlying neighbor
//	int kk;
//	int tt;
//	int ii;
//	int jj;
//	int ll;

    msn_old = snl;
    for (j = msn_old+1; j<MAXSNL; j++)
    {
        if (wice[j] <= 0.1)
        {
            wliq[j+1] = wliq[j+1] + wliq[j];
            wice[j+1] = wice[j+1] + wice[j];
            //cout<<tt+1<<"  "<<wliq[tt+1]<<" "<<wice[tt+1]<<endl;

            // shift all elements above this down one.
            if ((j>snl+1)&&( snl<MAXSNL-2))
            {
                for (i=j; i>=snl+2; i--)
                {
                    tss[i]  = tss[i-1];
                    wliq[i] = wliq[i-1];
                    wice[i] = wice[i-1];
                    dz[i]   = dz[i-1];
                }
            }

            snl = snl + 1;
            //          write(6,*) 'one snow layer is gone'
        }
        //cout<<MAXSNL+j-1<<" "<<tss[MAXSNL+j-1]<<" "<<wliq[MAXSNL+j-1]
        //                <<" "<<wice[MAXSNL+j-1]<<" "<<dz[MAXSNL+j-1]<<endl;
    }

    if (snl==MAXSNL-1)
    {
        scv = 0.;
        snowdp = 0.;
        return;
    }
    else
    {
        scv = 0.;
        snowdp = 0.;
        zwice = 0.;
        zwliq = 0.;
        for (j = snl+1; j<MAXSNL; j++)
        {
            scv = scv + wice[j] + wliq[j];
            snowdp = snowdp + dz[j];
            zwice = zwice + wice[j];
            zwliq = zwliq + wliq[j];
        }
    }
    //-----------------------------------------------------------------------
    // check the snow depth

    if (snowdp < 0.01)     //all snow gone
    {
        //cout<<"snow depth is less than 0.01"<<endl;
        //getchar();
        snl = MAXSNL-1;
        scv = zwice;
        if (scv <= 0.)
            snowdp = 0.;

        // the liquid water assumed ponding on soil surface
        wliq[MAXSNL] = wliq[MAXSNL] + zwliq;
        //       write(6,'(17h all snow is gone)')
        return;
    }
    else                        // snow layers combined
    {
        // two or more layers
        if (snl<MAXSNL-2)
        {
            msn_old = snl;
            mssi = 0;
            for (i = msn_old+1; i<MAXSNL; i++)
            {
                // If top node is removed, combine with bottom neighbor
                //cout<<i<<" "<<dz[MAXSNL+i-1]<<"  "<<dzmin[mssi-1]<<endl;
                if (dz[i] < dzmin[mssi])
                {
                    if (i == snl+1)
                        neibor = i + 1;
                    // If the bottom neighbor is not snow, combine with the top neighbor
                    else if (i == MAXSNL-1)
                        neibor = i - 1;
                    // If none of the above special cases apply, combine with the thinnest neighbor
                    else
                    {
                        neibor = i + 1;
                        if ((dz[i-1]+dz[i]) < (dz[i]+dz[i+1]))
                            neibor = i-1;
                    }


                    // Node l and j are combined and stored as node j.

                    if (neibor> i)
                    {
                        j = neibor;
                        l = i;
                    }
                    else
                    {
                        j = i;
                        l = neibor;
                    }

                    //cout<<i<<" "<<neibor<<endl;
                    //cout<<j<<" "<<dz[MAXSNL+j-1]<<" "<<wliq[MAXSNL+j-1]<<" "<<wice[MAXSNL+j-1]<<" "<<tss[MAXSNL+j-1]<<endl;
                    //cout<<l<<" "<<dz[MAXSNL+l-1]<<" "<<wliq[MAXSNL+l-1]<<" "<<wice[MAXSNL+l-1]<<" "<<tss[MAXSNL+l-1]<<endl;

                    combo ( dz[j], wliq[j], wice[j], tss[j],
                            dz[l], wliq[l], wice[l], tss[l] );

                    //cout<<j<<" "<<dz[MAXSNL+j-1]<<" "<<wliq[MAXSNL+j-1]<<" "<<wice[MAXSNL+j-1]<<" "<<tss[MAXSNL+j-1]<<endl;
                    //cout<<l<<" "<<dz[MAXSNL+l-1]<<" "<<wliq[MAXSNL+l-1]<<" "<<wice[MAXSNL+l-1]<<" "<<tss[MAXSNL+l-1]<<endl;

                    // Now shift all elements above this down one.

                    if ((j-1)> (snl+1))
                    {
                        for (k=j-1; k>=snl+2; k--)
                        {

                            tss[k]  = tss[k-1];
                            wice[k] = wice[k-1];
                            wliq[k] = wliq[k-1];
                            dz[k]   = dz[k-1];
                        }
                    }
                    snl = snl + 1;

                    //    write(6,'(7h Nodes ,i4,4h and,i4,14h combined into,i4)') l,j,j

                    if (snl >=MAXSNL-2)
                        break;
                    // The layer thickness great than the prescibed minimum value
                }
                else
                    mssi = mssi + 1 ;

            }
        }
    }

    // Reset the node depth and the depth of layer interface

    for (k = MAXSNL-1; k>=snl+1; k--)
    {
        z[k] = zi[k+1] - 0.5*dz[k];
        zi[k] = zi[k+1] - dz[k];
    }
}
void CommonLandModel::snowlayersdivide (int lb,                      int    &snl,                  double z[MAXSNL+MAXSOILL],  double dz[MAXSNL+MAXSOILL],double zi[MAXSNL+MAXSOILL+1],
                                        double wliq[MAXSNL+MAXSOILL],double wice[MAXSNL+MAXSOILL],double tss[MAXSNL+MAXSOILL])
{
    //=======================================================================
    // Original author : Yongjiu Dai, September 15, 1999
    //
    // subdivides snow layer when its thickness exceed the prescribed maximum
    //=======================================================================

    //use precision
    //implicit none

    //-------------------------- Dummy argument -----------------------------

    // integer, INTENT(in) :: lb              // lower bound of array
    // integer, INTENT(inout) :: snl          // Number of snow
    //real(r8), INTENT(inout) :: wice(lb:0)   // ice lens [kg/m2]
    //real(r8), INTENT(inout) :: wliq(lb:0)   // liquid water [kg/m2]
    //real(r8), INTENT(inout) :: tss (lb:0)   // Nodel temperature [K]
    //real(r8), INTENT(inout) :: dz  (lb:0)   // Layer thickness [m]
    //real(r8), INTENT(inout) :: z   (lb:0)   // Node depth [m]
    //real(r8), INTENT(inout) :: zi  (lb-1:0) // Depth of layer interface [m]

    //----------------------- Local variables ------------------------------

    // numbering from 1 (surface) msno (bottom)
    double  drr;      // thickness of the combined [m]
    double  dzsno[MAXSNL]; // Snow layer thickness [m]
    double  swice[MAXSNL]; // Partial volume of ice [m3/m3]
    double  swliq[MAXSNL]; // Partial volume of liquid water [m3/m3]
    double  tsno[MAXSNL];  // Nodel temperature [K]

    int k;            // number of do looping
    int msno;         // number of snow layer 1 (top) to msno (bottom)

    double  zwice,zwliq,propor;

    //-----------------------------------------------------------------------
    msno = MAXSNL-snl-1;
    for (k=0; k<msno; k++)
    {
        dzsno[k] = dz  [k+snl+1];
        swice[k] = wice[k+snl+1];
        swliq[k] = wliq[k+snl+1];
        tsno[k]  = tss [k+snl+1];
        //cout<<dzsno[k]<<" "<<swice[k]<<" "<<swliq[k]<<" "<<tsno[k]<<endl;
    }

    if (msno == 1)
    {
        if (dzsno[0]> 0.03)
        {
            msno = 2;
            // Specified a new snow layer
            dzsno[0] = dzsno[0]/2.0;
            swice[0] = swice[0]/2.0;
            swliq[0] = swliq[0]/2.0;
            dzsno[1] = dzsno[0];
            swice[1] = swice[0];
            swliq[1] = swliq[0];
            tsno[1]  = tsno[0];
            // cout<<dzsno[1]<<" "<<swice[1]<<" "<<swliq[1]<<" "<<tsno[1]<<endl;
            //        write(6,*)'Subdivided Top Node into two layer (1/2)'
        }
    }

    if (msno > 1)
    {
        if (dzsno[0]> 0.02)
        {
            drr = dzsno[0] - 0.02;
            propor = drr/dzsno[0];
            zwice = propor*swice[0];
            zwliq = propor*swliq[0];
            propor = 0.02/dzsno[0];
            swice[0] = propor*swice[0];
            swliq[0] = propor*swliq[0];
            dzsno[0] = 0.02;
            combo(dzsno[1],swliq[1],swice[1],tsno[1],
                  drr,     zwliq,   zwice,   tsno[0]);
            //cout<<dzsno[1]<<"  "<<swliq[1]<<"  "<<swice[1]<<"  "<<tsno[1]<<endl;
            //cout<<drr<<"  "<<zwliq<<"  "<<zwice<<"  "<<tsno[0]<<endl;

            //        write(6,*) 'Subdivided Top Node &
            //                    20 mm combined into underlying neighbor'

            if ((msno <= 2)&&(dzsno[1] > 0.07))
            {
                // subdivided a new layer
                msno = 3;
                dzsno[1] = dzsno[1]/2.0;
                swice[1] = swice[1]/2.0;
                swliq[1] = swliq[1]/2.0;

                dzsno[2] = dzsno[1];
                swice[2] = swice[1];
                swliq[2] = swliq[1];
                tsno[2]  = tsno[1];
            }
        }
    }
    if (msno > 2)
    {
        if (dzsno[1]> 0.05)
        {
            drr = dzsno[1] - 0.05;
            propor = drr/dzsno[1];
            zwice = propor*swice[1];
            zwliq = propor*swliq[1];

            propor = 0.05/dzsno[1];
            swice[1] = propor*swice[1];
            swliq[1] = propor*swliq[1];
            dzsno[1] = 0.05;

            combo(dzsno[2],swliq[2],swice[2],tsno[2],
                  drr,     zwliq,   zwice,   tsno[1]);
            //cout<<dzsno[2]<<"  "<<swliq[2]<<"  "<<swice[2]<<"  "<<tsno[2]<<endl;
            //cout<<drr<<"  "<<zwliq<<"  "<<zwice<<"  "<<tsno[1]<<endl;
            //        write(6,*)'Subdivided 50 mm from the subsface layer &
            //                   &and combined into underlying neighbor'

            if ((msno <= 3)&&(dzsno[2]>0.18))
            {
                // subdivided a new layer
                msno =  4;
                dzsno[2] = dzsno[2]/2.0;
                swice[2] = swice[2]/2.0;
                swliq[2] = swliq[2]/2.0;

                dzsno[3] = dzsno[2];
                swice[3] = swice[2];
                swliq[3] = swliq[2];
                tsno[3]  = tsno[2];
            }
        }
    }

    if (msno > 3)
    {
        if (dzsno[2]> 0.11)
        {
            drr = dzsno[2] - 0.11;
            propor = drr/dzsno[2];
            zwice = propor*swice[2];
            zwliq = propor*swliq[2];

            propor = 0.11/dzsno[2];
            swice[2] = propor*swice[2];
            swliq[2] = propor*swliq[2];
            dzsno[2] = 0.11;

            combo(dzsno[3],swliq[3],swice[3],tsno[3],
                  drr,     zwliq,   zwice,   tsno[2]);

            //        write(6,*)'Subdivided 110 mm from the third Node &
            //                   &and combined into underlying neighbor'

            if ((msno <= 4)&&(dzsno[3]> 0.41))
            {
                // subdivided a new layer
                msno = 5;
                dzsno[3] = dzsno[3]/2.0;
                swice[3] = swice[3]/2.0;
                swliq[3] = swliq[3]/2.0;

                dzsno[4] = dzsno[3];
                swice[4] = swice[3];
                swliq[4] = swliq[3];
                tsno[4]  = tsno[3];
            }
        }
    }

    if (msno > 4)
    {
        if (dzsno[3] > 0.23)
        {
            drr = dzsno[3] - 0.23;
            propor = drr/dzsno[3];
            zwice = propor*swice[3];
            zwliq = propor*swliq[3];

            propor = 0.23/dzsno[3];
            swice[3] = propor*swice[3];
            swliq[3] = propor*swliq[3];
            dzsno[3] = 0.23;

            combo(dzsno[4],swliq[4],swice[4],tsno[4],
                  drr,     zwliq,   zwice,   tsno[3]);

            //cout<<dzsno[4]<<"  "<<swliq[4]<<"  "<<swice[4]<<"  "<<tsno[4]<<endl;
            //cout<<drr<<"  "<<zwliq<<"  "<<zwice<<"  "<<tsno[3]<<endl;

            //        write(6,*)'Subdivided 230 mm from the fourth Node &
            //                   'and combined into underlying neighbor'
        }
    }

    snl  = MAXSNL-msno-1;

    for (k=snl+1; k<MAXSNL; k++)
    {
        dz[k]   = dzsno[k-snl-1];
        wice[k] = swice[k-snl-1];
        wliq[k] = swliq[k-snl-1];
        tss[k]  = tsno [k-snl-1];
        //cout<<k<<"  "<<dz[k+MAXSNL-1]<<" "<<wice[k+MAXSNL-1]<<" "<<wliq[k+MAXSNL-1]<<" "<<tss[k+MAXSNL-1]<<endl;
    }

    for (k=MAXSNL-1; k>=snl+1; k--)
    {
        z[k]    = zi[k+1] - 0.5*dz[k];
        zi[k] = zi[k+1] - dz[k];
        //cout<<k<<" "<<z[k+MAXSNL-1]<<"  "<<zi[k+MAXSNL-1]<<endl;
    }
}

void CommonLandModel::snowwater (int lb, double dtime, double ssi, double wimp,double pg_rain,
                                 double qseva, double qsdew, double qsubl, double qfros, double  dz[MAXSNL+MAXSOILL],
                                 double wice[MAXSNL+MAXSOILL],  double  wliq[MAXSNL+MAXSOILL], double &qout_snowb)
{

    //-----------------------------------------------------------------------
    // Original author : Yongjiu Dai, September 15, 1999
    //
    // Water flow wihtin snow is computed by an explicit and non-physical based scheme,
    // which permits a part of liquid water over the holding capacity (a tentative value
    // is used, i.e., equal to 0.033*porosity) to percolate into the underlying layer,
    // except the case of that the porosity of one of the two neighboring layers is
    // less than 0.05, the zero flow is assumed. The water flow out of the bottom
    // snow pack will participate as the input of the soil water and runoff.
    //-----------------------------------------------------------------------

    //use precision
    //use phycon_module, only : denice, denh2o  // physical constant
    //implicit none

    //----------------------- dummy argument --------------------------------
    //integer, INTENT(in) :: &
    //      lb          // lower bound of array

    //real(r8), INTENT(in) :: &
    //      dtime,    // time step (s)
    //      ssi,      // irreducible water saturation of snow
    //      wimp,     // water impremeable if porosity less than wimp
    //      dz(lb:0), // layer thickness (m)

    //      pg_rain,  // rainfall after removal of interception (mm h2o/s)
    //      qseva,    // ground surface evaporation rate (mm h2o/s)
    //      qsdew,    // ground surface dew formation (mm h2o /s) [+]
    //      qsubl,    // sublimation rate from snow pack (mm h2o /s) [+]
    //      qfros       // surface dew added to snow pack (mm h2o /s) [+]

    //real(r8), INTENT(inout) :: &
    //      wice(lb:0),&// ice lens (kg/m2)
    //      wliq(lb:0)  // liquid water (kg/m2)
    //
    //real(r8), INTENT(out) :: &
    //      qout_snowb  // rate of water out of snow bottom (mm/s)

    //----------------------- local variables --------------------------------
    //int j ;        // k do loop/array indices

    double
    qin,       // water flow into the elmement (mm/s)
    qout,      // water flow out of the elmement (mm/s)
//		zwice,     // the sum of ice mass of snow cover (kg/m2)
    wgdif,     // ice mass after minus sublimation
    vol_liq[MAXSNL],// partitial volume of liquid water in layer
    vol_ice[MAXSNL],// partitial volume of ice lens in layer
    eff_porosity[MAXSNL]; // effective porosity = porosity - vol_ice

    //=======================================================================
    // renew the mass of ice lens (wice) and liquid (wliq) in the surface snow layer,
    // resulted by sublimation (frost) / evaporation (condense)

    wgdif = wice[lb] + (qfros - qsubl)*dtime;
    wice[lb] = wgdif;
    if (wgdif < 0.)
    {
        wice[lb] = 0.0;
        wliq[lb] = wliq[lb] + wgdif;
    }
    wliq[lb] = wliq[lb] + (pg_rain + qsdew - qseva)*dtime;
    wliq[lb] = max(0., wliq[lb]);

    // Porosity and partitial volume
    for (int j = lb; j<MAXSNL; j++)
    {
        vol_ice[j] = min(1., wice[j]/(dz[j]*denice));
        eff_porosity[j] = 1. - vol_ice[j];
        vol_liq[j] = min(eff_porosity[j], wliq[j]/(dz[j]*denh2o));
    }

    // Capillary forces within snow are usually two or more orders of magnitude
    // less than those of gravity. Only gravity term are considered.
    // the genernal expression for water flow is "K * ss**3", however,
    // no effective paramterization for "K". Thus, a very simple consideration
    // (not physical based) is introduced:
    // when the liquid water of layer exceeds the layer's holding
    // capacity, the excess meltwater adds to the underlying neighbor layer.

    qin = 0.;
    for (int j = lb; j<MAXSNL; j++)
    {
        wliq[j] = wliq[j] + qin;
        if (j<=MAXSNL-2)
        {
            // no runoff over snow surface, just ponding on surface
            if ((eff_porosity[j]<wimp)||(eff_porosity[j+1]<wimp))
                qout = 0.0;
            else
            {
                qout = max(0.0,(vol_liq[j]-ssi*eff_porosity[j])*dz[j]);
                qout = min(qout,(1.-vol_ice[j+1]-vol_liq[j+1])*dz[j+1]);
                qout = max(0.0, qout);
            }
        }
        else
            qout = max(0.,(vol_liq[j]-ssi*eff_porosity[j])*dz[j]);

        qout    = qout*1000.0;
        wliq[j] = wliq[j] - qout;
        qin     = qout;

    }

    qout_snowb = qout/dtime;
}
void CommonLandModel::SOCEAN(bool    dosst,     double dtime,    double  &oro,   double  hu,    double  ht,
                             double  hq,       double     us,    double  vs,     double  tm,    double  qm,
                             double  rhoair,   double  psrf,     double  sabg,   double frl,   double  &tssea,
                             double  tssub[4],  double &scv,     double & taux,   double & tauy,  double & fsena,
                             double & fevpa,    double & lfevpa,   double & fseng,  double & fevpg, double & tref,
                             double & qref,     double & z0ma,     double & zol,    double & rib,   double & ustar,
                             double & qstar,    double & tstar,    double & u10m,   double & v10m,  double & f10m,
                             double & fm,       double & fh,       double & fq,     double & emis,  double & olrg)
{
//-----------------------------------------------------------------------
//           Simple Ocean Model
// 1. calculate sea surface fluxes, based on CLM
// 2. calculate sea surface albedos and seaice/snow temperatures
//                as in NCAR CCM3.6.16
// Original authors : yongjiu dai and xin-zhong liang (08/30/2001)
//-----------------------------------------------------------------------

    //use precision
    //use phycon_module, only : tfrz, hvap, hsub, stefnc
    //implicit none

//------------------------------Arguments--------------------------------



    //logical,  INTENT(IN) :: dosst   // true to update sst/ice/snow before calculation
    //real(r8), INTENT(in) :: dtime   // time-step (s)
    //real(r8), INTENT(in) :: hu      // agcm reference height of wind [m]
    //real(r8), INTENT(in) :: ht      // agcm reference height of temperature [m]
    //real(r8), INTENT(in) :: hq      // agcm reference height of humidity [m]
    //real(r8), INTENT(in) :: us      // wind component in eastward direction [m/s]
    //real(r8), INTENT(in) :: vs      // wind component in northward direction [m/s]
    //real(r8), INTENT(in) :: tm      // temperature at agcm reference height [kelvin]
    //real(r8), INTENT(in) :: qm      // specific humidity at agcm reference height [kg/kg]
    //real(r8), INTENT(in) :: rhoair  // density air [kg/m3]
    //real(r8), INTENT(in) :: psrf    // atmosphere pressure at the surface [pa] [not used]
    //real(r8), INTENT(in) :: sabg    // surface solar absorbed flux [W/m2]
    //real(r8), INTENT(in) :: frl     // downward longwave radiation [W/m2]

    //real(r8), INTENT(inout) :: oro  // ocean(0)/seaice(2)/ flag
    //real(r8), INTENT(inout) :: scv  // snow water equivalent depth (mm)
    //real(r8), INTENT(inout) :: tssub(plsice) // surface/sub-surface temperatures [K]
    //real(r8), INTENT(out) :: tssea  // sea surface temperature [K]

    //real(r8), INTENT(out) :: taux   // wind stress: E-W [kg/m/s**2]
    //real(r8), INTENT(out) :: tauy   // wind stress: N-S [kg/m/s**2]
    //real(r8), INTENT(out) :: fsena  // sensible heat from reference height to atmosphere [W/m2]
    //real(r8), INTENT(out) :: fevpa  // evaporation from refence height to atmosphere [mm/s]
    //real(r8), INTENT(out) :: lfevpa // laten heat from reference height to atmosphere [W/m2]
    //real(r8), INTENT(out) :: fseng  // sensible heat flux from ground [W/m2]
    //real(r8), INTENT(out) :: fevpg  // evaporation heat flux from ground [mm/s]

    //real(r8), INTENT(out) :: tref   // 2 m height air temperature [kelvin]
    //real(r8), INTENT(out) :: qref   // 2 m height air humidity
    //real(r8), INTENT(out) :: z0ma   // effective roughness [m]
    //real(r8), INTENT(out) :: zol    // dimensionless height (z/L) used in Monin-Obukhov theory
    //real(r8), INTENT(out) :: rib    // bulk Richardson number in surface layer
    //real(r8), INTENT(out) :: ustar  // friction velocity [m/s]
    //real(r8), INTENT(out) :: tstar  // temperature scaling parameter
    //real(r8), INTENT(out) :: qstar  // moisture scaling parameter
    //real(r8), INTENT(out) :: u10m   // 10m u-velocity
    //real(r8), INTENT(out) :: v10m   // 10m v-velocity
    //real(r8), INTENT(out) :: f10m   // integral of profile function for momentum at 10m
    //real(r8), INTENT(out) :: fm     // integral of profile function for momentum
    //real(r8), INTENT(out) :: fh     // integral of profile function for heat
    //real(r8), INTENT(out) :: fq     // integral of profile function for moisture
    //real(r8), INTENT(out) :: emis   // averaged bulk surface emissivity
    //real(r8), INTENT(out) :: olrg   // longwave up flux at surface [W/m2]

//-----------------------------------------------------------------------
    const int psrfty = 7;  // Number of surface types
    const int plsice = 4;  // number of seaice levels

    int    isrfty  ;   // surface type index (1-7)
    double  cgrndl ;   // deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
    double  cgrnds ;   // deriv of soil latent heat flux wrt soil temp [w/m**2/k]
    double  dshf   ;   // Ts partial derivative for sensible heat flux
    double  dlhf   ;   // Ts partial derivative for latent heat flux
    double  fnt    ;   // net surface flux for input conditions [W/m2]
    double  dfntdt ;   // net surface flux ts partial derivative [W/m2]
    double  tsbsf[plsice]  ;   // Non-adjusted srfc/sub-srfc temperatures
    double  snowh  ;   // snow depth (liquid water equivalent) [m]
    double  sicthk ;   // sea-ice thickness [m]

    double emisi  = 1.0   ;   // (0.97) surface emissivity for ice or snow [-]
    double emisw  = 1.0   ;   // (0.97) surface emissivity for water [-]
    double tsice  = 271.36;   // freezing point of sea ice [K]
    double thsice = 2.0   ;   // initial thickness of sea ice [m]
    double snsice = 0.005 ;   // initial snow water equivalent over sea ice [m]

    int j;

//-----------------------------------------------------------------------

    snowh = scv/1000.0;

    if (dosst)
    {
// update sea temperatures and sea ice distribution
// as well as snow cover over sea ice
        if (((int)oro==2)&&(tssea>tsice))
        {
            oro = 0.0;         // old sea ice melt out
            snowh = 0.0;
            scv = 0.0;
            sicthk = 0.0;
            for (j=1; j<plsice; j++)
                tssub[j-1] = tssea;
        }
        else if (((int)oro==0)&&(tssea<=tsice))
        {
            oro = 2.0;         // new sea ice formed
            snowh = snsice;
            scv = snowh*1000.0;
            sicthk = thsice;
            for (j=1; j<=plsice; j++)
                tssub[j] = tssea;
        }
    }

    tssea = tssub[0];

// compute surface fluxes, derviatives, and exchange coefficiants
    seafluxes (oro   ,hu,    ht,      hq,      us,
               vs,    tm,    qm,      rhoair,  psrf,
               tssea, taux,  tauy,    fsena,   fevpa,
               fseng, fevpg, tref,    qref,    z0ma,
               zol,   rib,   ustar,   qstar ,  tstar,
               u10m,  v10m,  f10m,    fm,      fh,
               fq,    cgrndl,cgrnds);

    if ((int)oro==0)            // ocean
    {
        lfevpa = fevpa*hvap;
        olrg = stefnc*emisw*pow(tssea,4 )+ (1.-emisw)*frl;
        emis = emisw;
    }
    else if ((int)oro==2)      // sea ice
    {
        lfevpa = fevpa*hsub;

        // net surface flux and derivate at current surface temperature
        dshf = cgrnds;
        dlhf = hsub*cgrndl;
        olrg = stefnc*emisi*pow(tssea,4) + (1.-emisi)*frl;
        fnt = sabg + frl - olrg - fsena - lfevpa;
        dfntdt = -(dshf + dlhf) - stefnc*emisi*4.*pow(tssea,3);

        // initialize surface/subsurface temperatures for srftsb

        for (j=1; j<=plsice; j++)
            tsbsf[j-1] = tssub[j-1];

// set sea ice surface type
        isrfty = 2;

        // diffusion calculation for temperature
        srftsb(isrfty,dtime,fnt,dfntdt,snowh,tsbsf);

        for (j=1; j<=plsice; j++)
        {
            tsbsf[j-1] = min(tsbsf[j-1],tfrz);
            tssub[j] = tsbsf[j-1];
        }
        tssea = tssub[0];

        olrg = stefnc*emisi*pow(tssea,4) + (1.-emisi)*frl;
        emis = emisi;
    }

}
void CommonLandModel::soilwater(int    nl_soil,              double dtime,               double wimp,             double smpmin,            double porsl[MAXSOILL],
                                double phi0[MAXSOILL],       double bsw[MAXSOILL],       double hksati[MAXSOILL], double z[MAXSNL+MAXSOILL],double dz[MAXSNL+MAXSOILL],
                                double zi[MAXSNL+MAXSOILL+1],double tss[MAXSNL+MAXSOILL],double vol_liq[MAXSOILL],double vol_ice[MAXSOILL], double eff_porosity[MAXSOILL],
                                double qinfl,                double etr,                 double rootr[MAXSOILL],  double dwat[MAXSOILL],    double hk[MAXSOILL],
                                double dhkdw[MAXSOILL] )
{

//=======================================================================
// Original author : Yongjiu Dai, September 15, 1999
//
// Soil moisture is predicted from a 10-layer model (as with soil temperatures),
// in which the vertical soil moisture transport is governed by infiltration,
// runoff, gradient diffusion, gravity, and root extraction by canopy transpiration;
// the net water applied to the surface layer is by snowmelt, precipitation,
// the throughfall of canopy dew, and minus surface runoff and evaporation.
//
// The vertical water flow in an unsaturated porous media is described by Darcy's law,
// and the hydraulic conductivity and the soil negative potential vary with
// soil water content and soil texture based on the work of Clapp and Hornberger (1978)
// and Cosby et al. (1984). The equation is integrated over the layer thickness,
// in which the time rate of change in water mass must equal the net flow across
// the bounding interface, plus the rate of internal source or sink.
// The terms of water flow across the layer interfaces are linearly expanded
// by using first-order Taylor expansion. The equations result in a tridiagonal
// system equation.
//
//             Note: length units here are all millimeter
//                   (in temperature subroutine uses same  soil layer
//                   structure required but lengths are m)
//
//                   Richards equation:
//
//                   d wat     d     d wat d psi
//                   ----- =  -- [ k(----- ----- - 1) ] + S
//                     dt     dz       dz  d wat
//
//                   where: wat = volume of water per volume of soil (mm**3/mm**3)
//                          psi = soil matrix potential (mm)
//                          dt  = time step (s)
//                          z   = depth (mm)
//                          dz  = thickness (mm)
//                          qin = inflow at top (mm h2o /s)
//                          qout= outflow at bottom (mm h2o /s)
//                          s   = source/sink flux (mm h2o /s)
//                          k   = hydraulic conductivity (mm h2o /s)
//
//                                 d qin                 d qin
//           qin[n+1] = qin[n] +  -------- d wat(j-1) + --------- d wat(j)
//                                 d wat(j-1)            d wat(j)
//                          ==================|=================
//                                            < qin
//
//                          d wat(j)/dt * dz = qin[n+1] - qout[n+1] + S(j)
//
//                                            > qout
//                          ==================|=================
//                                    d qout               d qout
//             qout[n+1] = qout[n] + --------- d wat(j) + --------- d wat(j+1)
//                                    d wat(j)             d wat(j+1)
//
//
//                   Solution: linearize k and psi about d wat and use tridiagonal
//                             system of equations to solve for d wat,
//                             where for layer j
//
//                   r_j = a_j [d wat_j-1] + b_j [d wat_j] + c_j [d wat_j+1]
//
//=======================================================================

// use precision
// use phycon_module, only : tfrz
// implicit none
//
//------------------------- dummy argument ------------------------------
// integer, INTENT(in) :: nl_soil  // soil layers

// real(r8), INTENT(in) :: &
//       dtime,            &// time step [s]
//       wimp,             &//water impremeable if porosity less than wimp
//       smpmin,           &//restriction for min of soil poten. (mm)

//       phi0(1:nl_soil),  &// saturated soil suction [mm]
//       bsw(1:nl_soil),   &// Clapp-Hornberger "B"
//       porsl(1:nl_soil), &// saturated volumetric soil water content(porosity)
//       hksati(1:nl_soil),&// hydraulic conductivity at saturation [mm h2o/s]
//       z(1:nl_soil),     &// layer depth [mm]
//       dz(1:nl_soil),    &// layer thickness [mm]
//       zi(0:nl_soil),    &// interface level below a "z" level [mm]
//
//       qinfl,            &// net water input from top [mm h2o/s]
//       etr,              &// actual transpiration [mm h2o/s]
//       rootr(1:nl_soil), &// root resistance of a layer, all layers add to 1.0
//eff_porosity(1:nl_soil), &// effective porosity
//     vol_liq(1:nl_soil), &// soil water per unit volume [mm/mm]
//     vol_ice(1:nl_soil), &// soil water per unit volume [mm/mm]
//     tss(1:nl_soil)       // soil temperature

// real(r8), INTENT(out) :: &
//       dwat (1:nl_soil), &// change of soil water [m3/m3]
//       hk   (1:nl_soil), &// hydraulic conductivity [mm h2o/s]
//       dhkdw(1:nl_soil)   // d(hk)/d(vol_liq)

//-------------------------- local variables -----------------------------
    double
    amx[MAXSOILL],    // "a" left off diagonal of tridiagonal matrix
        bmx[MAXSOILL],    // "b" diagonal column for tridiagonal matrix
        cmx[MAXSOILL],    // "c" right off diagonal tridiagonal matrix
        den,               // used in calculating qin, qout
        dqidw0,            // d(qin)/d(vol_liq(i-1))
        dqidw1,            // d(qin)/d(vol_liq(i))
        dqodw1,            // d(qout)/d(vol_liq(i))
        dqodw2,            // d(qout)/d(vol_liq(i+1))
        dsmpdw[MAXSOILL], // d(smp)/d(vol_liq)
        num,               // used in calculating qin, qout
        qin,               // flux of water into soil layer [mm h2o/s]
        qout,              // flux of water out of soil layer [mm h2o/s]
        rmx[MAXSOILL],    // "r" forcing term of tridiagonal matrix
        s_node,            // soil wetness
        s1,                // "s" at interface of layer
        s2,                // k*s**(2b+2)
        smp[MAXSOILL];      // soil matrix potential [mm]

    int j;                // do loop indices

//=======================================================================
// Set zero to hydraulic conductivity if effective porosity 5% in any of
// two neighbour layers or liquid content (theta) less than 0.001
    for (j=0; j<nl_soil; j++)
    {
        //cout<<j<<"  "<<eff_porosity[j]<<" "<<vol_liq[j]<<" "<<porsl[j]<<" "<<bsw[j]<<endl;
        if ((eff_porosity[j] < wimp)||(eff_porosity[min(nl_soil-1,j+1)]< wimp)||(vol_liq[j] <= 1.e-3))
        {
            hk[j] = 0.0;
            dhkdw[j] = 0.0;
        }
        else
        {
            s1 = 0.5*(vol_liq[j]/porsl[j]+vol_liq[min(nl_soil-1,j+1)]/porsl[min(nl_soil-1,j+1)]);
            s2 = hksati[j]*pow(s1,(2.*bsw[j]+2.));
            hk[j] = s1*s2  ;
            dhkdw[j] = (2.*bsw[j]+3.)*s2*0.5/porsl[j];
            if (j == nl_soil-1)
                dhkdw[j] = dhkdw[j] * 2.0;
        }
        //cout<<dhkdw[j]<<","<<hk[j]<<endl;
    }


// Evaluate hydraulic conductivity, soil matrix potential,
//     d(smp)/d(vol_liq), and d(hk)/d(vol_liq).
//
    for (j=0; j<nl_soil; j++)
    {
        if (tss[MAXSNL+j]>tfrz)
        {
            if (porsl[j]<1.e-6)    // bed rock
            {
                s_node = 0.001;
                smp[j] = -phi0[j];
                dsmpdw[j] = 0.0;
            }
            else
            {
                s_node = min(1.,vol_liq[j]/porsl[j]);
                //cout<<s_node<<endl;
                s_node = max(s_node, 0.001);
                smp[j] = -phi0[j]*pow(s_node,(-bsw[j]));
                //cout<<smp[j]<<" "<<-phi0[j]<<" "<<bsw[j]<<endl;
                smp[j] = max(smpmin, smp[j]);        // Limit soil suction
                dsmpdw[j] = -bsw[j]*smp[j]/(s_node*porsl[j]);
                //cout<<smp[j]<<" "<<dsmpdw[j]<<endl;
            }
        }
        else
// When ice is present, the matric potential is only related to temperature
// by (Fuchs et al., 1978: Soil Sci. Soc. Amer. J. 42(3):379-385)
// Unit 1 Joule = 1 (kg m2/s2), J/kg /(m/s2) ==> m ==> 1e3 mm
        {
            smp[j] = 1.e3 * 0.3336e6/9.80616*(tss[MAXSNL+j]-tfrz)/tss[MAXSNL+j];
            smp[j] = max(smpmin, smp[j]);        // Limit soil suction
            dsmpdw[j] = 0.;
        }

        // cout<<j+1<<"  "<<smp[j]<<"  "<<dsmpdw[j]<<endl;
    }
//-----------------------------------------------------------------------
// Set up r, a, b, and c vectors for tridiagonal solution
// Node j=1

    j      = 0;
    qin    = qinfl;
    den    = (z[MAXSNL+j+1]-z[MAXSNL+j]);
    num    = (smp[j+1]-smp[j]) - den;
    qout   = -hk[j]*num/den;
    dqodw1 = -(-hk[j]*dsmpdw[j]   + num*dhkdw[j])/den;
    dqodw2 = -( hk[j]*dsmpdw[j+1] + num*dhkdw[j])/den;

    rmx[j] =  qin - qout - etr*rootr[j];
    amx[j] =  0.;
    bmx[j] =  dz[MAXSNL+j]/dtime + dqodw1;
    cmx[j] =  dqodw2;
//-----------------------------------------------------------------------
// Nodes j=2 to j=nl_soil-1

    for (j=1; j<nl_soil-1; j++)
    {
        den    = (z[MAXSNL+j] - z[MAXSNL+j-1]);
        num    = (smp[j]-smp[j-1]) - den;
        qin    = -hk[j-1]*num/den;
        dqidw0 = -(-hk[j-1]*dsmpdw[j-1] + num*dhkdw[j-1])/den;
        dqidw1 = -( hk[j-1]*dsmpdw[j]   + num*dhkdw[j-1])/den;

        den    = (z[MAXSNL+j+1]-z[MAXSNL+j]);
        num    = (smp[j+1]-smp[j]) - den;
        qout   = -hk[j]*num/den;
        dqodw1 = -(-hk[j]*dsmpdw[j]   + num*dhkdw[j])/den;
        dqodw2 = -( hk[j]*dsmpdw[j+1] + num*dhkdw[j])/den;
        rmx[j] =  qin - qout - etr*rootr[j];
        amx[j] = -dqidw0;
        bmx[j] =  dz[MAXSNL+j]/dtime - dqidw1 + dqodw1;
        cmx[j] =  dqodw2;
    }

//-----------------------------------------------------------------------
// node j=nl_soil
    j      = nl_soil-1;
    den    = (z[MAXSNL+j] - z[MAXSNL+j-1]);
    num    = (smp[j]-smp[j-1]) - den;
    qin    = -hk[j-1]*num/den;
    dqidw0 = -(-hk[j-1]*dsmpdw[j-1] + num*dhkdw[j-1])/den;
    dqidw1 = -( hk[j-1]*dsmpdw[j]   + num*dhkdw[j-1])/den;

    qout   =  hk[j];
    dqodw1 =  dhkdw[j];

    rmx[j] =  qin - qout - etr*rootr[j];
    amx[j] = -dqidw0;
    bmx[j] =  dz[MAXSNL+j]/dtime - dqidw1 + dqodw1;
    cmx[j] =  0.;


//-----------------------------------------------------------------------
// Solve for dwat
    tridia (nl_soil ,amx ,bmx ,cmx ,rmx ,dwat);
    //for(int i=0;i<nl_soil;i++)
    //  cout<<dwat[i]<<endl;
}
void CommonLandModel::sortin( double eyy[6], double pco2y[6], double range,double gammas,int ic )
{
    if ( ic < 3 )
    {
        pco2y[0] = gammas + 0.5*range;
        pco2y[1] = gammas + range*( 0.5 - 0.3*sign(1.0,eyy[0]) );
        pco2y[2] = pco2y[0] - (pco2y[0]-pco2y[1])/(eyy[0]-eyy[1]+1.e-10)*eyy[0];

        double pmin = min( pco2y[0], pco2y[1] );
        double emin = min(   eyy[0],   eyy[1] );
        if ( emin > 0.0 && pco2y[2] > pmin ) pco2y[2] = gammas;
    }

    else
    {
        int n = ic - 1;
        int i=0;
        for (int j=1; j<=n; j++)
        {
            double a = eyy[j];
            double b = pco2y[j];
            for (i = j-1; i>=0; i--)
            {
                if (eyy[i] <= a ) break;//goto L200;
                eyy[i+1] = eyy[i];
                pco2y[i+1] = pco2y[i];
                //if (i==0) i=-1;
            }
            //i = -1;
            //L200:
            eyy[i+1] = a;
            pco2y[i+1] = b;
        }

        double pco2b = 0.0;
        int is    = 0;
        for (int ix = 0; ix<=n; ix++)
        {
            if ( eyy[ix] < 0.0 ) pco2b = pco2y[ix];
            if ( eyy[ix] < 0.0 ) is = ix;
        }

        int i1 = is-1;
        i1 = max(0, i1);
        i1 = min(n-2, i1);
        int i2 = i1 + 1;
        int i3 = i1 + 2;
        int isp   = is + 1;
        isp = min( isp, n );
        is = isp - 1;

        double pco2yl = pco2y[is] - (pco2y[is]-pco2y[isp])/(eyy[is]-eyy[isp])*eyy[is];

        // method using a quadratic fit
        double ac1 = eyy[i1]*eyy[i1] - eyy[i2]*eyy[i2];
        double ac2 = eyy[i2]*eyy[i2] - eyy[i3]*eyy[i3];
        double bc1 = eyy[i1] - eyy[i2];
        double bc2 = eyy[i2] - eyy[i3];
        double cc1 = pco2y[i1] - pco2y[i2];
        double cc2 = pco2y[i2] - pco2y[i3];
        double bterm = (cc1*ac2-cc2*ac1)/(bc1*ac2-ac1*bc2);
        double aterm = (cc1-bc1*bterm)/ac1;
        double cterm = pco2y[i2] - aterm*eyy[i2]*eyy[i2] - bterm*eyy[i2];
        double pco2yq = cterm;
        pco2yq = max( pco2yq, pco2b );
        pco2y[ic] = ( pco2yl+pco2yq)/2.0;
    }

    pco2y[ic] = max ( pco2y[ic], 0.01 );
    //return;
}
void CommonLandModel::srftsb(int isrfty,double dtime,double fnt,double dfntdt,double snowh,double tsbsf[4])
{

//-----------------------------------------------------------------------
// Compute surface and subsurface temperatures over sea-ice surfaces.
//
// Sea ice temperatures are specified in 'plsice' layers of fixed
// thickness and thermal properties.  The forecast temperatures are
// determined from a backward/implicit diffusion calculation using
// linearized sensible/latent heat fluxes. The bottom ocean temperature
// is fixed at -2C, allowing heat flux exchange with underlying ocean.
//
// Sub-surface layers are indexed 1 at the surface, increasing downwards
// to plsice.  Layers have mid-points and interfaces between layers.
//
// Temperatures are defined at mid-points, while fluxes between layers
// and the top/bottom media are defined at layer interfaces.
//
//-----------------------------------------------------------------------

    // use precision
    // use phycon_module, only: tkice, tkair
    // implicit none

//------------------------------Arguments--------------------------------

    //integer, parameter :: psrfty = 7  // Number of surface types
    //integer, parameter :: plsice = 4  // number of seaice levels

    //integer, INTENT(in) :: isrfty  // surface type index (1 - 7)
    //real(r8), INTENT(in) :: dtime  // time step (s)
    //real(r8), INTENT(in) :: fnt    // top surface/atmosphere net energy flux
    //real(r8), INTENT(in) :: dfntdt // ts partial derivative of net sfc flux
    //real(r8), INTENT(in) :: snowh  // snow depth (liquid water equivalent) [m]

    //real(r8), INTENT(inout) :: tsbsf[plsice] // surface/sub-surface tmps

//---------------------------Local variables-----------------------------

    int j, jndx;           // sub-surface layer index
    const int psrfty = 7;  // Number of surface types
    const int plsice = 4;  // number of seaice levels


    double  cmass[plsice];  // specific heat of soil (J/kg/K)
    double  rho[plsice];    // mass densty of sub-sfc mat (kg/m3)
    double  tk[plsice];     // thermal conductivity (watts/m/K)
    double  diag[plsice];   // diagonal matrix elements
    double  htsrc[plsice];  // external heat source (W/m3)
    double  rhs[plsice];    // rhs of tri-diagonal matrix equation
    double  sbdiag[plsice]; // sub-diagonal matrix elements
    double  spdiag[plsice]; // super-diagonal matrix elements
    double  tin[plsice];    // initial sub-surface temperatures
    double  z[plsice+1];    // interface geometrical depth (m)
//   double  ws[plsice];     // working storage for mtdlss

    double  cmty  ;  // layer mass heat capacity
    double  fbt   ;  // ocean heat flux into sea-ice
    double  rhty  ;  // layer mass density
    double  thck  ;  // layer thickness
    double  tkty  ;  // layer thermal conductivity
    double  cmsnow;  // Snow mass heat capacity
    double  crt   ;  // cmass*rho*rdtime
    double  delz  ;  // layer thickness
    double  delzmn;  // thick from mid-point to lyr above mid-point
    double  delzpl;  // thick from mid-point to lyr below mid-point
    double  fmns  ;  // 1/(delz*delzmn)
    double  fpls  ;  // 1/(delz*delzpl)
    double  msnow ;  // mass path of snow
    double  mlice ;  // mass path of ice
    double  rdtime;  // inverse model time step
    double  rhsnow;  // snow mass density
    double  rztop ;  // 1/ztop
    double  tkbot ;  // bottom layer top interf thermal conduct
    double  tkmns ;  // layer bottom interface thermal conduct
    double  tkpls ;  // layer top interface thermal conductivity
    double  tksnow;  // snow thermal conducitivity
    double  tktop ;  // top layer bottom interface thermal conduct
    double  tmp   ;  // crt - dfntdt(i)*rztop
    double  zbot  ;  // bottom layer thickness
    double  zm    ;  // present layer mid-point depth
    double  zmmn  ;  // layer above mid-point depth
    double  zmpl  ;  // layer below mid-point depth
    double  zsnow ;  // snow geometric depth
    double  ztop  ;  // top layer thickness
    bool    scvr  ;  // true if surface snow covered

//--------------------------Data Statements------------------------------
// specified (and invariant) thermal properties for surface types


    const double   cmair  = 1.00e+3; // mass specific heat of air [J/kg/K]
    const double   cmice  = 2.07e+3; // mass specific heat of ice [J/kg/K]
    const double   frcair = 0.90;   // fraction of air assumed in mix of ice
    const double   rhair  = 1.25;   // mass density of surface air [kg/m3]
    const double   rhice  = 9.20e+2; // mass density of ice [kg/m3]
    const double   snwedp = 10.0;   // snow:water equivalent depth factor [-]

    //mass specific heat (J/kg/K)
    const double cmtype[psrfty][plsice]= {{4.20e+3,4.20e+3,4.20e+3,4.20e+3},
        {2.07e+3,2.07e+3,2.07e+3,2.07e+3},
        {2.07e+3,2.07e+3,2.07e+3,2.07e+3},
        {1.04e+3,1.04e+3,1.04e+3,1.04e+3},
        {7.20e+2,7.20e+2,7.20e+2,7.20e+2},
        {5.60e+2,5.60e+2,5.60e+2,5.60e+2},
        {4.16e+2,4.16e+2,4.16e+2,4.16e+2}
    };
    // mass density (kg/m3)
    const double rhtype[psrfty][plsice]= {{1.00e+3,1.00e+3,1.00e+3,1.00e+3},
        {9.20e+2,9.20e+2,9.20e+2,9.20e+2},
        {9.20e+2,9.20e+2,9.20e+2,9.20e+2},
        {2.50e+3,2.50e+3,2.50e+3,2.50e+3},
        {2.50e+3,2.50e+3,2.50e+3,2.50e+3},
        {2.50e+3,2.50e+3,2.50e+3,2.50e+3},
        {2.50e+3,2.50e+3,2.50e+3,2.50e+3}
    };

//layer thicknesses (m)
    const double thckly[psrfty][plsice]= {{2.0 , 5.0,   10.0,  33.0},
        {0.5,  0.5,   0.5,   0.5},
        {0.25, 0.5,   0.5,   8.5},
        {0.5,  0.366, 1.369, 6.99},
        {0.09, 0.39,  1.459, 7.45},
        {0.08, 0.435, 1.628, 8.31},
        {0.12, 0.492, 1.841, 9.4 }
    };


    const double tktype[psrfty][plsice]= {{15.0,  15.0,  15.0,  15.0},
        {2.2,   2.2,   2.2,   2.2},
        {2.2,   2.2,   2.2,   2.2},
        {1.408, 1.408 ,1.408, 1.408},
        {1.104, 1.104, 1.104, 1.104},
        {1.071, 1.071, 1.071, 1.071},
        {1.019, 1.019, 1.019, 1.019}
    };


    rdtime = 1./dtime;

// calculate snow properties
    cmsnow = (1.-frcair)*cmice + frcair*cmair;
    rhsnow = (1.-frcair)*rhice + frcair*rhair;
    tksnow = (1.-frcair)*tkice + frcair*tkair;

// no external heat source
    for (j=0; j<plsice; j++)
        htsrc[j] = 0.0;

// define logical for snow covered surfaces:
    scvr = (snowh>0.0);

// define thermal properities for each sub/surface layer, starting
// with the top layer
    jndx    = isrfty;
    thck = thckly[jndx-1][0];
    cmty = cmtype[jndx-1][0];
    rhty = rhtype[jndx-1][0];
    tkty = tktype[jndx-1][0];

// initialize fields for no snow cover
    z[0]     = 0.0;
    z[1]     = thck;
    cmass[0] = cmty;
    rho[0]   = rhty;
    tk[0]    = tkty;

// modify layer 1 fields for snow cover if present
// snow equivlnt depth times snow liquid water depth gives the physical
// depth of snow for thermal conduction computation; snow is mixed
// uniformly by mass with the top surface layer
    if (scvr)
    {
        zsnow    = snowh*snwedp;
        msnow    = rhsnow*zsnow;
        mlice    = rhty*thck;
        rho[0]   = (msnow*rhsnow + mlice*rhty)/(msnow+mlice);
        cmass[0] = (msnow*cmsnow + mlice*cmty)/(msnow+mlice);
        tk[0]    = (msnow*tksnow + mlice*tkty)/(msnow+mlice);
        z[1]     = (msnow+mlice) / rho[0];
    }

// set surface thermal properties for the lower sub/surface layers:
    for (j=2; j<=plsice; j++)
    {
        jndx       = isrfty;
        thck       = thckly[jndx-1][j-1];
        cmass[j-1] = cmtype[jndx-1][j-1];
        rho[j-1]   = rhtype[jndx-1][j-1];
        tk[j-1]    = tktype[jndx-1][j-1];
        z[j]     = z[j-1] + thck;
    }

// define set of linear equations for temperature
    for (j=0; j<plsice; j++)
        tin[j] = tsbsf[j];


// if sea ice, compute heat flux from underlying ocean, assumed to be at
// the temperature of -2C
    fbt = 0.0;
    if (isrfty==2)
    {
        zbot = 0.5*(z[plsice] - z[plsice-1]);
        fbt = -tk[plsice-1]*(271.16 - tin[plsice-1])/zbot;
    }
// set up linear equations
    sbdiag[0]      = 0.0;
    spdiag[plsice-1] = 0.0;

// single layer
    if (plsice==1)
    {
        rztop = 1./(z[1] - z[0]);
        crt = (cmass[0]*rho[0]*rdtime);
        diag[0] = crt - dfntdt*rztop;
        rhs[0] = diag[0]*tin[0] + fnt*rztop - fbt*rztop + htsrc[0];
    }

// more than one layer: top layer first
    else if (plsice>1)
    {
        crt       = cmass[0]*rho[0]*rdtime;
        ztop      = z[1] - z[0];
        rztop     = 1./ztop;
        tktop     = 0.5*(tk[0] + tk[1]);
        zmpl      = 0.5*(z[2] + z[1]);
        zm        = 0.5*(z[1] + z[0]);
        delzpl    = zmpl - zm;
        fpls      = 1./(ztop*delzpl);
        tmp       = crt - dfntdt*rztop;
        diag[0]   = tmp + tktop*fpls;
        spdiag[0] = -tktop*fpls;
        rhs[0]    = tmp*tin[0] + fnt*rztop + htsrc[0];

// intermediate layers
        for (j=2; j<=plsice-1; j++)
        {
            crt       = cmass[j-1]*rho[j-1]*rdtime;
            delz      = z[j] - z[j-1];
            zmpl      = 0.5*(z[j+1] + z[j]);
            zm        = 0.5*(z[j]   + z[j-1]);
            zmmn      = 0.5*(z[j-1] + z[j-2]);
            delzpl    = zmpl - zm;
            delzmn    = zm - zmmn;
            fpls      = 1./(delz*delzpl);
            fmns      = 1./(delz*delzmn);
            tkpls     = 0.5*(tk[j]+tk[j-1]);
            tkmns     = 0.5*(tk[j-1]+tk[j-2]);
            sbdiag[j-1] = -tkmns*fmns;
            diag[j-1]   = crt + (tkpls*fpls + tkmns*fmns);
            spdiag[j-1] = -tkpls*fpls;
            rhs[j-1]    = crt*tin[j-1] + htsrc[j-1];
        }

// bottom layer
        crt       = cmass[plsice-1]*rho[plsice-1]*rdtime;
        zbot      = z[plsice] - z[plsice-1];
        zm        = 0.5*(z[plsice]   + z[plsice-1]);
        zmmn      = 0.5*(z[plsice-1] + z[plsice-2]);
        delzmn    = zm - zmmn;
        tkbot     = 0.5*(tk[plsice-2] + tk[plsice-1]);
        fmns      = 1./(zbot*delzmn);
        sbdiag[plsice-1] = -tkbot*fmns;
        diag[plsice-1] = crt + (tkbot*fmns);
        rhs[plsice-1] = crt*tin[plsice-1] - fbt/zbot + htsrc[plsice-1];
    }
    //cout<<sbdiag[0]<<" "<<sbdiag[1]<<" "<<sbdiag[2]<<" "<<sbdiag[3]<<endl;
    //cout<<diag[0]<<" "<<diag[1]<<" "<<diag[2]<<" "<<diag[3]<<endl;
    //cout<<spdiag[0]<<" "<<spdiag[1]<<" "<<spdiag[2]<<" "<<spdiag[3]<<endl;
    //cout<<rhs[0]<<" "<<rhs[1]<<" "<<rhs[2]<<" "<<rhs[3]<<endl;

    if (plsice==1)
        tsbsf[0] = rhs[0]/diag[0];
    else
        tridia(plsice,sbdiag,diag,spdiag,rhs,tsbsf);
}
void CommonLandModel::stomata(double   vmax25, double    effcon, double    slti, double    hlti, double    shti,
                              double   hhti,   double    trda,   double    trdm, double    trop, double    gradm,
                              double   binter, double    tm,     double    psrf, double    po2m, double    pco2m,
                              double    pco2a, double    ea,     double    ei,   double    tlef, double    par,
                              double    rb,    double    ra,     double    rstfac,double    cint[3],double &   assim,
                              double &   respc, double &   rst)
{

//=======================================================================
//
//     calculation of canopy photosynthetic rate using the integrated
//     model relating assimilation and stomatal conductance.
//
//     Original author: Yongjiu Dai, 08/11/2001
//
//     units are converted from mks to biological units in this routine.
//
//                          units
//                         -------
//
//      pco2m, pco2a, pco2i, po2m                : pascals
//      co2a, co2s, co2i, h2oa, h2os, h2oa       : mol mol-1
//      vmax25, respcp, assim, gs, gb, ga        : mol m-2 s-1
//      effcon                                   : mol co2 mol quanta-1
//      1/rb, 1/ra, 1/rst                        : m s-1
//
//                       conversions
//                      -------------
//
//      1 mol h2o           = 0.018 kg
//      1 mol co2           = 0.044 kg
//      h2o (mol mol-1)     = ea / psrf ( pa pa-1 )
//      h2o (mol mol-1)     = q*mm/(q*mm + 1)
//      gs  (co2)           = gs (h2o) * 1./1.6
//      gs  (mol m-2 s-1 )  = gs (m s-1) * 44.6*tf/t*p/po
//      par (mol m-2 s-1 )  = par(w m-2) * 4.6*1.e-6
//      mm  (molair/molh2o) = 1.611
//
//----------------------------------------------------------------------

// use precision
// implicit none

//real(r8),intent(in) :: &
//     effcon,       &// quantum efficiency of RuBP regeneration (mol CO2 / mol quanta)
//     vmax25,       &// maximum carboxylation rate at 25 C at canopy top
//                    // the range : 30.e-6 <-> 100.e-6 (mol co2 m-2 s-1)

//     slti,         &// slope of low temperature inhibition function      (0.2)
//     hlti,         &// 1/2 point of low temperature inhibition function  (288.16)
//     shti,         &// slope of high temperature inhibition function     (0.3)
//     hhti,         &// 1/2 point of high temperature inhibition function (313.16)
//     trda,         &// temperature coefficient in gs-a model             (1.3)
//     trdm,         &// temperature coefficient in gs-a model             (328.16)
//     trop,         &// temperature coefficient in gs-a model             (298.16)
//     gradm,        &// conductance-photosynthesis slope parameter
//     binter         // conductance-photosynthesis intercept

//real(r8),intent(in) :: &
//     tm,           &// atmospheric air temperature (K)
//     psrf,         &// surface atmospheric pressure (pa)
//     po2m,         &// O2 concentration in atmos. (20900 pa)
//     pco2m,        &// CO2 concentration in atmos. (35 pa)
//     pco2a,        &// CO2 concentration in canopy air space (pa)
//     ea,           &// canopy air space vapor pressure (pa)
//     ei,           &// saturation h2o vapor pressure in leaf stomata (pa)
//     tlef,         &// leaf temperature (K)
//     par,          &// photosynthetic active radiation (W m-2)

//     rb,           &// boundary resistance from canopy to cas (s m-1)
//     ra,           &// aerodynamic resistance from cas to refence height (s m-1)
//     rstfac         // canopy resistance stress factors to soil moisture

//real(r8),intent(in), dimension(3) :: &
//     cint           // scaling up from leaf to canopy

//real(r8),intent(out) :: &// ATTENTION : all for canopy not leaf
//     assim,        &// canopy assimilation rate (mol m-2 s-1)
//     respc,        &// canopy respiration (mol m-2 s-1)
//     rst            // canopy stomatal resistance (s m-1)

//-------------------- local --------------------------------------------

    double
    c3,           // c3 vegetation : 1; 0 for c4
    c4,           // c4 vegetation : 1; 0 for c3
    qt,           // (tleaf - 298.16) / 10
    gammas,       // co2 compensation point (pa)
    kc,           // Michaelis-Menten constant for co2
    ko,           // Michaelis-Menten constant for o2
    rrkk,         // kc (1+o2/ko)

    templ,        // intermediate value
    temph,        // intermediate value

    vm,           // maximum catalytic activity of Rubison (mol co2 m-2 s-1)
    jmax25,       // potential rate of whole-chain electron transport at 25 C
    jmax,         // potential rate of whole-chain electron transport (mol electron m-2 s-1)
    epar,         // electron transport rate (mol electron m-2 s-1)
    respcp,       // respiration fraction of vmax (mol co2 m-2 s-1)
    bintc,        // residual stomatal conductance for co2 (mol co2 m-2 s-1)

    rgas,         // universal gas contant (8.314 J mol-1 K-1)
    tprcor,       // coefficient for unit transfer
    gbh2o,        // one side leaf boundary layer conductance (mol m-2 s-1)
    gah2o,        // aerodynamic conductance between cas and reference height (mol m-2 s-1)
    gsh2o,        // canopy conductance (mol m-2 s-1)

    atheta,       // wc, we coupling parameter
    btheta,       // wc & we, ws coupling parameter
    omss,         // intermediate calcuation for oms
    omc,          // rubisco limited assimilation (omega-c: mol m-2 s-1)
    ome,          // light limited assimilation (omega-e: mol m-2 s-1)
    oms,          // sink limited assimilation (omega-s: mol m-2 s-1)
    omp,          // intermediate calcuation for omc, ome

    co2m,         // co2 concentration in atmos (mol mol-1)
    co2a,         // co2 concentration at cas (mol mol-1)
    co2s,         // co2 concentration at canopy surface (mol mol-1)
    co2st,        // co2 concentration at canopy surface (mol mol-1)
//      co2i,         // internal co2 concentration (mol mol-1)
    pco2in,       // internal co2 concentration at the new iteration (pa)
    pco2i,        // internal co2 concentration (pa)
    es,           // canopy surface h2o vapor pressure (pa)

    sqrtin,       // intermediate calculation for quadratic
    assmt,        // net assimilation with a positive limitation (mol co2 m-2 s-1)
    assimn,       // net assimilation (mol co2 m-2 s-1)
    hcdma,        // a-1
    aquad,        // a: ax^2 + bx + c = 0
    bquad,        // b: ax^2 + bx + c = 0
    cquad ;      // c: ax^2 + bx + c = 0

    double
    eyy[6],      // differnce of pco2i at two iteration step
        pco2y[6],    //
        range;        //

    int ic,ic1;

//=======================================================================

    c3 = 0.0;
    if ( effcon>0.07 ) c3 = 1.0;
    c4 = 1. - c3 ;

//-----------------------------------------------------------------------
// dependence on leaf temperature
//     gammas - CO2 compensation point in the absence of day respiration
//     ko     - Michaelis-Menton constant for carboxylation by Rubisco
//     kc     - Michaelis-Menton constant for oxygenation by Rubisco
//-----------------------------------------------------------------------

    qt = 0.1*( tlef - trop );

    kc = 30.0 * pow(2.1,qt);
    ko = 30000.  * pow(1.2,qt) ;
    gammas = 0.5 * po2m / (2600. *pow(0.57,qt)) * c3;        // = 0. for c4 plant ???

    rrkk = kc * ( 1. + po2m/ko ) * c3;

//----------------------------------------------------------------------
// maximun capacity
// vm     - maximum catalytic activity of Rubisco in the presence of
//          saturating level of RuP2 and CO2 (mol m-2s-1)
// jmax   - potential rate of whole-chain electron transport (mol m-2s-1)
// epar   - electron transport rate for a given absorbed photon radiation
// respc  - dark resipration (mol m-2s-1)
// omss   - capacity of the leaf to export or utilize the products of photosynthesis.
// binter - coefficient from observation, 0.01 for c3 plant, 0.04 for c4 plant
//-----------------------------------------------------------------------

    vm = vmax25 * pow(2.1,qt);        // (mol m-2 s-1)
    templ = 1. + std::exp(slti*(hlti-tlef));
    temph = 1. + std::exp(shti*(tlef-hhti));
    vm = vm / temph * rstfac * c3 + vm / (templ*temph) * rstfac * c4;
    vm = vm * cint[0];

    rgas = 8.314;                 // universal gas constant (J mol-1 K-1)
//---> jmax25 = 2.39 * vmax25 - 14.2e-6        // (mol m-2 s-1)
    jmax25 = 2.1 * vmax25;        // (mol m-2 s-1)
    jmax = jmax25 * std::exp( 37.e3 * (tlef - trop) / (rgas*trop*tlef) ) *
           ( 1. + std::exp( (710.*trop-220.e3)/(rgas*trop) ) ) /
           ( 1. + std::exp( (710.*tlef-220.e3)/(rgas*tlef) ) ) ;
    // 37000  (J mol-1)
    // 220000 (J mol-1)
    // 710    (J K-1)

    jmax = jmax * rstfac;
    jmax = jmax * cint[1];

    epar = min(4.6e-6 * par * effcon, 0.25*jmax);

    respcp = 0.015 * c3 + 0.025 * c4;
    respc = respcp * vmax25 * pow(2.0,qt) / ( 1. + std::exp( trda*(tlef-trdm )) ) * rstfac;
//     respc = 0.75e-6 * std::exp(er(?)*(1./71.02 - 1./(tlef-46.02)) )
//     respc = 0.7e-6 * 2.0**qt / ( 1. + std::exp( trda*(tlef-trdm )) ) * rstfac
    respc = respc * cint[0];
    omss = ( vmax25/2. ) * pow(1.8,qt) / templ * rstfac * c3
           + ( vmax25/5. ) * pow(1.8,qt) * rstfac * c4 ;
    omss = omss * cint[0];

    bintc = binter * max( 0.1, rstfac );
    bintc = bintc * cint[2];

//-----------------------------------------------------------------------
    tprcor = 44.6*273.16*psrf/1.013e+5;

// one side leaf boundary layer conductance for water vapor [=1/(2*rb)]
// ATTENTION: rb in CLM is for one side leaf, but for SiB2 rb for
// 2-side leaf, so the gbh2o shold be " 0.5/rb * tprcor/tlef "
//     gbh2o  = 0.5/rb * tprcor/tlef                    // mol m-2 s-1
    gbh2o  = 1./rb * tprcor/tlef  ;                  // mol m-2 s-1
// rb is for single leaf, but here the flux is for canopy, thus
//修复数组越界问题(wlx)
    gbh2o  = gbh2o * cint[2];

//  aerodynamic condutance between canopy and reference height atmosphere
    gah2o  = 1.0/ra * tprcor/tm;                     // mol m-2 s-1

//-----------------------------------------------------------------------
//     first guess is midway between compensation point and maximum
//     assimilation rate. // pay attention on this iteration

    co2m = pco2m/psrf ;                              // mol mol-1
    co2a = pco2a/psrf;

    range = pco2m * ( 1. - 1.6/gradm ) - gammas;

    for (ic=1; ic<=6; ic++)
    {
        pco2y[ic-1] = 0.0;
        eyy[ic-1] = 0.0;
    }


    for (ic=1; ic<=6; ic++)
    {
        sortin(eyy, pco2y, range, gammas, ic-1);      // legacy from SiB2
        pco2i =  pco2y[ic-1];

//-----------------------------------------------------------------------
//                      NET ASSIMILATION
//     the leaf assimilation (or gross photosynthesis) rate is described
//     as the minimum of three limiting rates:
//     omc: the efficiency of the photosynthetic enzyme system (Rubisco-limited);
//     ome: the amount of PAR captured by leaf chlorophyll;
//     oms: the capacity of the leaf to export or utilize the products of photosynthesis.
//     to aviod the abrupt transitions, two quadratic equations are used:
//             atheta*omp^2 - omp*(omc+ome) + omc*ome = 0
//         btheta*assim^2 - assim*(omp+oms) + omp*oms = 0
//-----------------------------------------------------------------------

        atheta = 0.877;
        btheta = 0.95;

        omc = vm   * ( pco2i-gammas ) / ( pco2i + rrkk ) * c3 + vm * c4 ;
        ome = epar * ( pco2i-gammas ) / ( pco2i+2.*gammas ) * c3 + epar * c4 ;
        oms   = omss * c3 + omss*pco2i * c4;

        sqrtin= max( 0., ((ome+omc)*(ome+omc) - 4.*atheta*ome*omc ) ) ;
        omp   = ( ( ome+omc ) - std::sqrt( sqrtin ) ) / ( 2.*atheta );
        sqrtin= max( 0., ( (omp+oms)*(omp+oms) - 4.*btheta*omp*oms ) );
        assim = ( ( oms+omp ) - std::sqrt( sqrtin ) ) / ( 2.*btheta );

        assimn= ( assim - respc) ;                       // mol m-2 s-1
//-----------------------------------------------------------------------
//                      STOMATAL CONDUCTANCE
//
//  (1)   pathway for co2 flux
//                                                  co2m
//                                                   o
//                                                   |
//                                                   |
//                                                   <  |
//                                        1.37/gsh2o >  |  Ac-Rd-Rsoil
//                                                   <  v
//                                                   |
//                                     <--- Ac-Rd    |
//     o------/\/\/\/\/\------o------/\/\/\/\/\------o
//    co2i     1.6/gsh2o     co2s    1.37/gbh2o     co2a
//                                                   | ^
//                                                   | | Rsoil
//                                                   | |
//
//  (2)   pathway for water vapor flux
//
//                                                  em
//                                                   o
//                                                   |
//                                                   |
//                                                   <  ^
//                                           1/gsh2o >  | Ea
//                                                   <  |
//                                                   |
//                                     ---> Ec       //
//     o------/\/\/\/\/\------o------/\/\/\/\/\------o
//     ei       1/gsh2o      es       1/gbh2o       ea
//                                                   | ^
//                                                   | | Eg
//                                                   | |
//
//  (3)   the relationship between net assimilation and tomatal conductance :
//        gsh2o = m * An * [es/ei] / [pco2s/p] + b
//        es = [gsh2o *ei + gbh2o * ea] / [gsh2o + gbh2o]
//        ===>
//        a*gsh2o^2 + b*gsh2o + c = 0
//
//-----------------------------------------------------------------------

//-->  co2a = co2m - 1.37/max(0.446,gah2o) * (assimn - 0.)     // mol mol-1
        co2s = co2a - 1.37*assimn/gbh2o ;                 // mol mol-1

        co2st = min( co2s, co2a );
        co2st = max( co2st,1.e-5 );

        assmt = max( 1.e-12, assimn );
        hcdma = ei*co2st / ( gradm*assmt );

        aquad = hcdma ;
        bquad = gbh2o*hcdma - ei - bintc*hcdma;
        cquad = -gbh2o*( ea + hcdma*bintc ) ;

        sqrtin= max( 0., ( bquad*bquad - 4.*aquad*cquad ) );
        gsh2o = ( -bquad + std::sqrt( sqrtin ) ) / (2.*aquad);
        es  = ( gsh2o-bintc ) * hcdma;                   // pa
        es  = min( es, ei );
        es  = max( es, 1.e-2);
        gsh2o = es/hcdma + bintc;                        // mol m-2 s-1
        pco2in = ( co2s - 1.6 * assimn / gsh2o )*psrf ;  // pa
        eyy[ic-1] = pco2i - pco2in ;                       // pa
//-----------------------------------------------------------------------

        if ( abs(eyy[ic-1])< 0.1)
            break;
    }


// convert gsh2o (mol m-2 s-1) to resistance rst ( s m-1)
    rst   = min( 1.e6, 1./(gsh2o*tlef/tprcor) ) ;    // s m-1
}
void CommonLandModel::subsurfacerunoff (int    nl_soil,               double dtime,       double pondmx,         double dzmm[MAXSNL+MAXSOILL],double wliq[MAXSNL+MAXSOILL],
                                        double eff_porosity[MAXSOILL],double hk[MAXSOILL],double dhkdw[MAXSOILL],double dwat[MAXSOILL],double &rsubst)

{
    //=======================================================================
    // Original author : Yongjiu Dai, 07/29/2002
    //=======================================================================

    //use precision
    //implicit none

    //-----------------------Arguments---------------------------------------

    // integer, INTENT(in) :: nl_soil   // number of soil layers
    // real(r8), INTENT(in) :: &
    //        dtime,           &// time step (s)
    //        pondmx,          &// ponding depth (mm)
    //        dzmm(1:nl_soil), &// layer thickness (mm)
    //eff_porosity(1:nl_soil), &// effective porosity = porosity - vol_ice
    //          hk(1:nl_soil), &// hydraulic conductivity (mm h2o/s)
    //       dhkdw(1:nl_soil), &// d(hk)/d(vol_liq)
    //        dwat(1:nl_soil)   // change in soil water

    // real(r8), INTENT(inout) :: wliq(1:nl_soil)  // liquid water (kg/m2)
    // real(r8), INTENT(out) :: rsubst             // subsurface runoff (mm h2o/s)

    //-----------------------Local Variables---------------------------------

    double  xs;           // excess soil water above saturation
    int     i;            // loop counter

    //-----------------------End Variable List-------------------------------

    rsubst = 0.0;
    xs = 0.0;

    // determine water in excess of saturation and
    // add excess water back to underlying soil layers
    // any left over put into runoff

    for (i=nl_soil-1; i>=1; i--)
    {
        wliq[MAXSNL+i] = wliq[MAXSNL+i] + xs;
        //cout<<MAXSNL+i<<"  "<<wliq[MAXSNL+i]<<endl;
        xs = max(0., wliq[MAXSNL+i]-eff_porosity[i]*dzmm[MAXSNL+i]);    // [mm]
        //cout<<xs<<"  "<<eff_porosity[i]<<"  "<<dzmm[i]<<"  "<<wliq[MAXSNL+i]<<"  "<<(wliq[MAXSNL+i]-eff_porosity[i]*dzmm[i])<<endl;
        wliq[MAXSNL+i] = min(eff_porosity[i]*dzmm[MAXSNL+i], wliq[MAXSNL+i]);
        //cout<<MAXSNL+i<<"  "<<xs<<"  "<<eff_porosity[i]<<"  "<<dzmm[i]<<"  "<<wliq[MAXSNL+i]<<endl;
    }
    wliq[MAXSNL] = wliq[MAXSNL] + xs;
    xs = max(0., wliq[MAXSNL]-(pondmx+eff_porosity[0]*dzmm[MAXSNL]));
    wliq[MAXSNL] = min(pondmx+eff_porosity[0]*dzmm[MAXSNL], wliq[MAXSNL]);

    // sub-surface runoff and drainage
    //cout<<xs<<" "<<hk[nl_soil-1]<<"  "<<dhkdw[nl_soil-1]<<" "<<dwat[nl_soil-1]<<endl;
    rsubst = rsubst + xs/dtime
             + hk[nl_soil-1] + dhkdw[nl_soil-1]*dwat[nl_soil-1]; // [mm/s]

}
void CommonLandModel::surfacerunoff(int    nl_soil,          double wtfact,           double wimp,                  double bsw[MAXSOILL],      double porsl[MAXSOILL],
                                    double phi0[MAXSOILL],   double hksati[MAXSOILL], double z[MAXSNL+MAXSOILL],    double dz[MAXSNL+MAXSOILL],double zi[MAXSNL+MAXSOILL+1],
                                    double vol_liq[MAXSOILL],double vol_ice[MAXSOILL],double eff_porosity[MAXSOILL],double gwat,               double &rsur)

{
    //=======================================================================
    // the original code was provide by Robert E. Dickinson based on following clues:
    // a water table level determination level added including highland and
    // lowland levels and fractional area of wetland (water table above the surface.
    // Runoff is parametrized from the lowlands in terms of precip incident on
    // wet areas and a base flow, where these are estimated using ideas from TOPMODEL.
    //
    // Author : Yongjiu Dai, 07/29/2002
    //=======================================================================

    //use precision
    //implicit none

    //-----------------------Arguments---------------------------------------

    // integer, INTENT(in) :: nl_soil   // number of soil layers
    // real(r8), INTENT(in) :: &
    //       wtfact,        &// fraction of model area with high water table
    //       wimp,          &// water impremeable if porosity less than wimp
    //      bsw(1:nl_soil), &// Clapp-Hornberger "B"
    //    porsl(1:nl_soil), &// saturated volumetric soil water content(porosity)
    //     phi0(1:nl_soil), &// saturated soil suction (mm)
    //   hksati(1:nl_soil), &// hydraulic conductivity at saturation (mm h2o/s)
    //        z(1:nl_soil), &// layer depth (m)
    //       dz(1:nl_soil), &// layer thickness (m)
    //       zi(0:nl_soil), &// interface level below a "z" level (m)
    //  vol_liq(1:nl_soil), &// partial volume of liquid water in layer [-]
    //  vol_ice(1:nl_soil), &// partial volume of ice lens in layer [-]
    //eff_porosity(1:nl_soil), &// effective porosity = porosity - vol_ice
    //       gwat            // net water input from top

    // real(r8), INTENT(out) :: rsur    // surface runoff (mm h2o/s)

    //-----------------------Local Variables---------------------------------

    double fcov;         // fractional area with water table at surface
//	double fz;           // coefficient for water table depth
    double zwt ;         // water table depth
    double deficit1;     // water deficit from 10-layer soil moisture
    double deficit2;     // water deficit from 60-layer
    double dzfine;       // layer thickness for fine mesh soil layers
    double zfine[100];    // depth for fine mesh soil layers
    double s1,su,v;      // local variables to calculate qinmax
    double qinmax;       // maximum infiltration capability
    double watsatm;     // averaged soil porosity (-)
    double sucsatm;      // averaged minimum soil suction (mm)
    double bswm;         // averaged Clapp and Hornberger "b"

    int i  ;           // loop counter

    //-----------------------End Variable List-------------------------------

    // average soil properties

    watsatm = 0.0;
    sucsatm = 0.0;
    bswm    = 0.0;

    for (i=0; i<nl_soil; i++)
    {
        watsatm = watsatm + porsl[i];
        sucsatm = sucsatm + phi0[i];
        bswm    = bswm    + bsw[i];
    }

    watsatm = watsatm / nl_soil;
    sucsatm = sucsatm / nl_soil;
    bswm    = bswm    / nl_soil;

    // Determine water table
    deficit1 = 0.0;
    for (i=0; i<nl_soil; i++)
    {
        deficit1 = deficit1+ (porsl[i]-(vol_ice[i]+vol_liq[i]))*dz[MAXSNL+i];
    }
//correction 04/06/2007
    dzfine = 2.0*zi[MAXSNL+nl_soil]/100.0;
    for (i =0; i<100; i++)
        zfine[i] = (i+1)*dzfine;

    deficit2 = 0.0;
    zwt = 2.0*zi[MAXSNL+nl_soil] - 0.001;
    for (i=0; i<100; i++)
    {
        //a = (zwt-zfine[i])*1000.0/sucsatm;
        //a = max(1.e-3, a);
        //belta = pow(a,(-1./bswm));
        //deficit2 = deficit2 + belta*watsatm*dzfine;
        deficit2=deficit2+watsatm*(1.0-pow(1.0+(zwt-zfine[i])*1000.0/sucsatm,-1.0/bswm))*dzfine;
        if (fabs(deficit2-deficit1)<=0.1)
        {
            zwt = zfine[i];
            break;
        }
    }
    // Saturation fraction
    fcov = wtfact;

    if (eff_porosity[0]>=wimp)
        fcov = wtfact*min(1.,std::exp(-zwt));
    // Maximum infiltration capacity
    s1 = max(0.01,vol_liq[0])/max(wimp, eff_porosity[0]);
    su = max(0.,(s1-fcov)/(1.-fcov));
    v = -bsw[0]*phi0[0]/(0.5*dz[MAXSNL]*1000.);
    qinmax = (1.+v*(su-1.))*hksati[0];
    if (eff_porosity[0]<wimp)
        qinmax = 0.;
    // Surface runoff
    rsur = fcov*gwat + (1.-fcov)*max(0.,gwat-qinmax);
}
void CommonLandModel::THERMAL(int   itypwat ,          int lb      ,                    int   nl_soil,    		   double dtime  ,          	 double trsmx0 ,
                              double zlnd  ,           double zsno    ,			     double csoilc ,			   double dewmx  ,  			 double capr   ,
                              double cnfac   ,		   double csol[MAXSOILL]   ,         double porsl[MAXSOILL] ,      double phi0[MAXSOILL],		 double bsw[MAXSOILL] ,
                              double dkmg[MAXSOILL]   ,double dkdry[MAXSOILL]  ,         double dksatu[MAXSOILL],	   double lai    ,  			 double sai    ,
                              double z0m     ,         double displa  ,                 double sqrtdi ,      		   double rootfr[MAXSOILL] ,	 double effcon ,
                              double vmax25  ,         double slti    ,                 double hlti   ,  			   double shti   ,  			 double hhti   ,
                              double trda    ,         double trdm    ,                 double trop   ,      		   double gradm  ,		    	 double binter ,
                              double extkn   ,         double hu      ,                 double ht     ,   		   double hq     ,  			 double us     ,
                              double vs      ,         double tm      ,                 double qm     ,  			   double rhoair ,		     double psrf   ,
                              double pco2m   ,         double po2m    ,                 double coszen ,  			   double parsun ,   			 double parsha ,
                              double sabvsun ,         double sabvsha ,                 double sabg   ,  			   double frl    ,			     double extkb  ,
                              double extkd   ,         double thermk  ,                 double fsno   ,      		   double sigf   ,			     double dz[15] ,
                              double z[MAXSNL+MAXSOILL],double zi[MAXSNL+MAXSOILL+1] ,    double &tlsun  ,          	   double &tlsha  ,			     double tss[MAXSNL+MAXSOILL],
                              double wice[MAXSNL+MAXSOILL],double wliq[MAXSNL+MAXSOILL],  double &ldew   ,              double &scv    ,		    	 double &snowdp ,
                              int  imelt[MAXSNL+MAXSOILL], double &taux    ,              double &tauy   ,              double &fsena  ,			     double &fevpa  ,
                              double &lfevpa  ,         double &fsenl   ,                 double &fevpl  ,              double &etr    ,  			 double &fseng  ,
                              double &fevpg   ,         double &olrg    ,                 double &fgrnd  ,              double rootr[MAXSOILL] ,      double &qseva  ,
                              double &qsdew   ,         double &qsubl   ,                 double &qfros  ,              double &sm     ,			     double &tref   ,
                              double &qref    ,         double &trad    ,                 double &rst    ,              double &assim  ,		    	 double &respc  ,
                              double &errore  ,         double &emis    ,                 double &z0ma   ,              double &zol    ,			     double &rib    ,
                              double &ustar   ,         double &qstar   ,                 double &tstar  ,              double &u10m   ,			     double &v10m   ,
                              double &f10m    ,         double &fm      ,                 double &fh     ,              double &fq)
{
    ////=======================================================================
    //// this is the main subroutine to execute the calculation
    //// of thermal processes and surface fluxes
    ////
    //// Original author : Yongjiu Dai, 09/15/1999; 08/30/2002
    ////
    //// FLOW DIAGRAM FOR THERMAL.F90
    ////
    //// THERMAL ===> qsadv
    ////              groundfluxes
    ////              eroot                      |dewfraction
    ////              leaftemone |               |qsadv
    ////              leaftemtwo |  ---------->  |moninobukini
    ////                                         //moninobuk
    ////                                         |stomata
    ////
    ////              groundTem     ---------->   meltf
    ////
    ////=======================================================================


    int i,j;

    double
    cgrnd,        // deriv. of soil energy flux wrt to soil temp [w/m2/k]
    cgrndl,       // deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
    cgrnds,       // deriv of soil latent heat flux wrt soil temp [w/m**2/k]
    degdT,        // d(eg)/dT
    dqgdT,        // d(qg)/dT
    dlrad,        // downward intwave radiation blow the canopy [W/m2]
    eg,           // water vapor pressure at temperature T [pa]
    egsmax,       // max. evaporation which soil can provide at one time step
    egidif,       // the excess of evaporation over "egsmax"
    emg,          // ground emissivity (0.97 for snow,
    // glaciers and water surface; 0.96 for soil and wetland)
    //errore,       // energy balnce error [w/m2]
    etrc,         // maximum possible transpiration rate [mm/s]
    fac,          // soil wetness of surface layer
    fact[MAXSNL+MAXSOILL], // used in computing tridiagonal matrix
    fsun,         // fraction of sunlit canopy
    hr,           // relative humidity
    htvp,         // latent heat of vapor of water (or sublimation) [j/kg]
    olru,         // olrg excluding dwonwelling reflection [W/m2]
    olrb,         // olrg assuming blackbody emission [W/m2]
    psit,         // negative potential of soil
    par,          // PAR absorbed by canopy [W/m2]
    qg,           // ground specific humidity [kg/kg]
    qsatg,        // saturated humidity [kg/kg]
    qsatgdT,      // d(qsatg)/dT
    qred,         // soil surface relative humidity
    rstfac,       // factor of soil water stress
    sabv,         // solar absorbed by canopy [W/m2]
    thm,          // intermediate variable (tm+0.0098*ht)
    th,           // potential temperature (kelvin)
    thv,          // virtual potential temperature (kelvin)
    tl,           // leaf temperature
    tg,           // ground surface temperature [K]
    tssbef[MAXSNL+MAXSOILL], // soil/snow temperature before update
    tinc,         // temperature difference of two time step
    ur,           // wind speed at reference height [m/s]
    ulrad,        // upward intwave radiation above the canopy [W/m2]
    wice0[MAXSNL+MAXSOILL],// ice mass from previous time-step
    wliq0[MAXSNL+MAXSOILL],// liquid mass from previous time-step
    wx,           // patitial volume of ice and water of surface layer
    xmf;            // total latent heat of phase change of ground water

    double z0ma_g,zol_g,rib_g,ustar_g,qstar_g,tstar_g;
    double f10m_g, fm_g,fh_g,fq_g,temp1,temp2,temp12m,temp22m,um,obu;

    //double *fact = new double[maxsnl+nl_soil];
    //double *tssbef = new double[maxsnl+nl_soil];
    //double *wice0 = new double[maxsnl+nl_soil];
    //double *wliq0 = new double[maxsnl+nl_soil];

    //=======================================================================
    // [1] Initial set and propositional variables
    //=======================================================================
    // fluxes
    taux   = 0.0 ;
    tauy   = 0.0 ;
    fsena  = 0.0 ;
    fevpa  = 0.0 ;
    lfevpa = 0.0 ;
    fsenl  = 0.0 ;
    fevpl  = 0.0 ;
    etr    = 0.0 ;
    fseng  = 0.0 ;
    fevpg  = 0.0 ;
    dlrad  = 0.0 ;
    ulrad  = 0.0 ;
    cgrnds = 0.0 ;
    cgrndl = 0.0 ;
    cgrnd  = 0.0 ;
    tref   = 0.0 ;
    qref   = 0.0 ;
    rst    = 2.0e+4;
    assim  = 0.0 ;
    respc  = 0.0 ;

    emis   = 0.0 ;
    z0ma   = 0.0 ;
    zol    = 0.0 ;
    rib    = 0.0 ;
    ustar  = 0.0 ;
    qstar  = 0.0 ;
    tstar  = 0.0 ;   //rootr  = 0.0 ;

    // temperature and water mass from previous time step
    tg = tss[lb];


    for (i=lb; i<MAXSNL+nl_soil; i++)
    {
        tssbef[i] = tss[i];
        wice0[i] = wice[i];
        wliq0[i] = wliq[i];
    }

    // emissivity
    emg = 0.96;
    if ((scv>0.0)||(itypwat==3))
        emg = 0.97;

    // latent heat, assumed that the sublimation occured only as wliq=0
    htvp = hvap;
    if ((wliq[lb]<=0.0)&&(wice[lb]>0.0))
        htvp = hsub;

    // potential temperatur at the reference height
    thm = tm + 0.0098*ht;              // intermediate variable equivalent to
    // tm*(pgcm/psrf)**(rgas/cpair)
    th = tm*pow((100000./psrf),(rgas/cpair)); // potential T
    thv = th*(1.+0.61*qm);             // virtual potential T
    ur = max(0.1,std::sqrt(us*us+vs*vs));   // limit set to 0.1
    //=======================================================================
    // [2] specific humidity and its derivative at ground surface
    //=======================================================================

    qred = 1.0;
    qsadv(tg,psrf,eg,degdT,qsatg,qsatgdT);
    //cout<<"tg:  "<<tg<<endl;

    if (itypwat<=1)           // soil ground
    {
        wx   = (wliq[MAXSNL]/denh2o + wice[MAXSNL]/denice)/dz[MAXSNL];

        if (porsl[0]<1.e-6)    // bed rock
            fac  = 0.001;
        else
        {
            fac  = min(1.,wx/porsl[0]);
            fac  = max( fac, 0.001 );
        }
        psit = -phi0[0] * pow(fac,(- bsw[0]) );   // psit = max(smpmin, psit)
        hr   = std::exp(psit/roverg/tg);
        qred = (1.-fsno)*hr + fsno;
    }
    qg = qred*qsatg;
    dqgdT = qred*qsatgdT;

    if ((qsatg > qm)&&(qm > (qred*qsatg)))
    {
        qg = qm;
        dqgdT = 0.0;
    }

    //=======================================================================
    // [3] Compute sensible and latent fluxes and their derivatives with respect
    //     to ground temperature using ground temperatures from previous time step.
    //=======================================================================
    if (sigf <= 0.999)
    {
        groundfluxes(zlnd,   zsno,   hu,    ht,     hq,
                     us,     vs,     tm,    qm,     rhoair,
                     psrf,   ur,     thm,   th,     thv,
                     tg,     qg,     dqgdT, htvp,   fsno,
                     sigf,   cgrnd,  cgrndl,cgrnds, taux,
                     tauy,   fsena,  fevpa, fseng,  fevpg,
                     tref,   qref,	  z0ma_g,zol_g,  rib_g,
                     ustar_g,qstar_g,tstar_g,f10m_g,fm_g,
                     fh_g,   fq_g);
    }
    //=======================================================================
    // [4] Canopy temperature, fluxes from the canopy
    //=======================================================================
    par = parsun + parsha;
    sabv = sabvsun + sabvsha;

    if (sigf >= 0.001)

        // soil water strees factor on stomatal resistance
    {
        eroot (lb,nl_soil,trsmx0,porsl,bsw,phi0,rootfr,
               dz,tss,wliq,rootr,etrc,rstfac);
        // fraction of sunlit and shaded leaves of canopy
        fsun = ( 1. - std::exp(-extkb*lai) ) / max( extkb*lai, 1.e-6 );
        if ((coszen<=0.0)||(sabv<1))
            fsun = 0.0;
        if (fsun<=0.1)
        {
            fsun = 0.0;
            tl = tlsha;
            leaftemone (dtime ,csoilc ,dewmx  ,htvp   ,lai    ,
                        sai     ,displa   ,sqrtdi ,z0m    ,effcon ,vmax25 ,
                        slti    ,hlti     ,shti   ,hhti   ,trda   ,trdm   ,
                        trop    ,gradm    ,binter ,extkn  ,extkb  ,extkd  ,
                        hu      ,ht       ,hq     ,us     ,vs     ,thm    ,
                        th      ,thv      ,qm     ,psrf   ,rhoair ,par    ,
                        sabv    ,frl      ,thermk ,rstfac ,po2m   ,pco2m  ,
                        sigf    ,etrc     ,tg     ,qg     ,dqgdT  ,emg    ,
                        tl      ,ldew     ,taux   ,tauy   ,fseng  ,fevpg  ,
                        cgrnd   ,cgrndl   ,cgrnds ,tref   ,qref   ,rst    ,
                        assim   ,respc    ,fsenl  ,fevpl  ,etr    ,dlrad  ,
                        ulrad   ,z0ma     ,zol    ,rib    ,ustar  ,qstar  ,
                        tstar   ,f10m     ,fm     ,fh     ,fq);
            tlsun = tl;
            tlsha = tl;
        }

        else
        {
            leaftemtwo (dtime ,csoilc  ,dewmx  ,htvp   ,lai    ,
                        sai     ,displa   ,sqrtdi  ,z0m    ,effcon ,vmax25 ,
                        slti    ,hlti     ,shti    ,hhti   ,trda   ,trdm   ,
                        trop    ,gradm    ,binter  ,extkn  ,extkb  ,extkd  ,
                        hu      ,ht       ,hq      ,us     ,vs     ,thm    ,
                        th      ,thv      ,qm      ,psrf   ,rhoair ,parsun ,
                        parsha  ,sabvsun  ,sabvsha ,frl    ,fsun   ,thermk ,
                        rstfac  ,po2m     ,pco2m   ,sigf   ,etrc   ,tg     ,
                        qg      ,dqgdT    ,emg     ,tlsun  ,tlsha  ,ldew   ,
                        taux    ,tauy     ,fseng   ,fevpg  ,cgrnd  ,cgrndl ,
                        cgrnds  ,tref     ,qref    ,rst    ,assim  ,respc  ,
                        fsenl   ,fevpl    ,etr     ,dlrad  ,ulrad  ,z0ma   ,
                        zol     ,rib      ,ustar   ,qstar  ,tstar  ,f10m   ,
                        fm      ,fh       ,fq);
        }
    }



    // equate canopy temperature to air over bareland.
    // required as sigf=0 carried over to next time step
    if (sigf < 0.001)
    {
        tlsun = tm;
        tlsha = tm;
        ldew = 0.0;
    }
    //=======================================================================
    // [5] Gound temperature
    //=======================================================================
    groundtem (itypwat,lb,nl_soil,dtime,capr,
               cnfac,csol,porsl,dkmg,dkdry,
               dksatu,sigf,dz,z,zi,
               tss,wice,wliq,scv,snowdp,
               frl,dlrad,sabg,fseng,fevpg,
               cgrnd,htvp,emg,imelt,
               sm,xmf,fact);
    //=======================================================================
    // [6] Correct fluxes to present soil temperature
    //=======================================================================

    tg = tss[lb];
    tinc = tss[lb] - tssbef[lb];
    fseng = fseng + tinc*cgrnds ;
    fevpg = fevpg + tinc*cgrndl;

    // calculation of evaporative potential; flux in kg m-2 s-1.
    // egidif holds the excess energy if all water is evaporated
    // during the timestep.  this energy is later added to the sensible heat flux.

    egsmax = (wice[lb]+wliq[lb]) / dtime;

    egidif = max( 0.0, fevpg - egsmax );
    fevpg = min ( fevpg, egsmax );
    fseng = fseng + htvp*egidif;

    // total fluxes to atmosphere
    fsena = fsenl + fseng;
    fevpa = fevpl + fevpg;
    lfevpa= hvap*fevpl + htvp*fevpg;   // w/m2 (accouting for sublimation)

    qseva = 0.0;
    qsubl = 0.0;
    qfros = 0.0;
    qsdew = 0.0;

    if (fevpg >= 0.0)
        // not allow for sublimation in melting (melting ==> evap. ==> sublimation)
    {
        qseva = min(wliq[lb]/dtime, fevpg);
        qsubl = fevpg - qseva;
    }
    else if (tg < tfrz)
        qfros = fabs(fevpg);
    else
        qsdew = fabs(fevpg);




    // ground heat flux
    fgrnd =  sabg + dlrad + (1.-sigf)*emg*frl
             - emg*stefnc*pow(tssbef[lb],3)*(tssbef[lb] + 4.*tinc)
             - (fseng+fevpg*htvp);

    // outgoing int-wave radiation from canopy + ground
    olrg =   ulrad + (1.-sigf)*(1.-emg)*frl
             + (1.-sigf)*emg*stefnc * pow(tssbef[lb],4)
             + 4.*emg*stefnc*pow(tssbef[lb],3)*tinc;

    // for conservation we put the increase of ground intwave to outgoing

    // averaged bulk surface emissivity
    olrb = stefnc*pow(tssbef[lb],3)*((1.-sigf)*tssbef[lb] + 4.*tinc);
    olru = ulrad + emg*olrb;
    olrb = ulrad + olrb;
    emis = olru / olrb;

    // radiative temperature
    trad = pow((olrg/stefnc),0.25);

    // additonal variables required by WRF and RSM model
    if (sigf < 0.001)
    {
        ustar = ustar_g;
        tstar = tstar_g;
        qstar = qstar_g;
        rib   = rib_g;
        zol   = zol_g;
        z0ma  = z0ma_g;
        f10m  = f10m_g;
        fm    = fm_g;
        fh    = fh_g;
        fq    = fq_g;
    }

    else if (sigf <= 0.7)
    {
        z0ma  = sigf*z0ma  + (1.-sigf)*z0ma_g;
        // assumed um ~= ur here
        um = ur;
        ustar =   std::sqrt(max(1.e-6,std::sqrt(taux*taux+tauy*tauy))/rhoair);
        tstar = - fsena/(cpair*ustar*rhoair);
        qstar = - fevpa/(ustar*rhoair);
        zol = (hu-displa)*vonkar*grav*(tstar+0.61*th*qstar)/(ustar*ustar*thv);
        if (zol >= 0.0)  //stable
            zol = min(2.,max(zol,1.e-6));
        else                  //unstable
            zol = max(-100.,min(zol,-1.e-6));

        obu = (hu-displa)/zol;
        moninobuk(hu,ht,hq,displa,z0ma,
                  z0ma,z0ma,obu,um,ustar,
                  temp1,temp2,temp12m,temp22m,f10m,
                  fm,fh,fq);
        //rib = min(5.,zol*pow(vonkar,3)*ustar*ustar/(temp1*um*um));

        // correct 27, December, 2006
        rib = min(5.,zol*ustar*ustar/(vonkar*temp1*um*um));
    }
    else
    {
        ustar = ustar;
        tstar = tstar;
        qstar = qstar;
        rib   = rib;
        zol   = zol;
        z0ma  = z0ma;
        f10m  = f10m;
        fm    = fm;
        fh    = fh;
        fq    = fq;
    }
    u10m  = us/ur * ustar/vonkar * f10m;
    v10m  = vs/ur * ustar/vonkar * f10m;

    //=======================================================================
    // [7] energy balance error
    //=======================================================================
    errore = sabv + sabg + frl - olrg - fsena - lfevpa - xmf;
    for (j=lb; j<MAXSNL+nl_soil; j++)
        errore -=  (tss[j]-tssbef[j])/fact[j];
    // if(abs(errore)>.2)then
    // write(6,*) 'THERMAL.F90 : energy  balance violation'
    // write(6,100) errore,sabv,sabg,frl,olrg,fsenl,fseng,hvap*fevpl,htvp*fevpg,xmf
    // endif
    //delete [] fact;
    //delete [] tssbef;
    //delete [] wice0 ;
    //delete [] wliq0 ;

}
void CommonLandModel::ticktime (double dtime, int idate[3])
{

//=======================================================================
// Original author: Yongjiu Dai, September 15, 1999
//
//      Notes: Greenwich time
//
//=======================================================================

//  implicit none

//  real, INTENT(in) :: dtime          // model time step (second)
//  integer, INTENT(inout) :: idate(3) // year and julian day, sec, respectively
    int  maxday ;                    // days of one year (365 or 366)

//-----------------------------------------------------------------------

    idate[2] = idate[2] + int(dtime);
    if (idate[2]>=86400)
    {
        idate[2] = idate[2] - 86400;
        if ((idate[0]%4==0)&&((idate[0]%100)!=0)||(idate[0]%400==0))
            maxday = 366;
        else
            maxday = 365;
        idate[1] = idate[1] + 1;
        if (idate[1]>maxday)
        {
            idate[0] = idate[0] + 1;
            idate[1] = 1;
        }
    }
}
void CommonLandModel::tridia (int n, double a[], double b[], double c[], double r[],double u[])
{
    //use precision
    //implicit none
    //integer, intent(in) :: n        !length of diagonal element vector
    //real(r8), intent(in) :: a(1:n)  !subdiagonal elements
    //real(r8), intent(in) :: b(1:n)  !diagonal elements
    //real(r8), intent(in) :: c(1:n)  !superdiagonal elements
    //real(r8), intent(in) :: r(1:n)  !right hand side
    //real(r8), intent(out) :: u(1:n) !solution vector


    double  bet;
    double *gam = new double [n];
    int j;
// -----------------------------------------------------------------


    bet = b[0];
    u[0] = r[0] / bet;
    for (j=2; j<=n; j++)
    {
        gam[j-1] = c[j-2]/bet;
        bet = b[j-1] -a[j-1]*gam[j-1];
        u[j-1] = (r[j-1]-a[j-1]*u[j-2])/bet;
    }
    for (j=n-1; j>=1; j--)
    {
        u[j-1] = u[j-1] - gam[j]*u[j];
    }

    delete []gam;
}
void CommonLandModel::twostream (double   chil,   double ref[2][2],  double tran[2][2], double green,       double lai,
                                 double   sai,    double coszen,     double albg[2][2], double albv[2][2],  double tranc[2][2],
                                 double   &thermk, double &extkb,      double &extkd,      double ssun[2][2],  double ssha[2][2] )
{
    //-----------------------------------------------------------------------
    //
    //     calculation of canopy albedos via two stream approximation (direct
    //     and diffuse ) and partition of incident solar
    //
    // Original author: Yongjiu Dai, June 11, 2001
    //
    //-----------------------------------------------------------------------

    //use precision
    //implicit none

    // parameters
    //  real(r8), intent(in) :: &
    //          // static parameters associated with vegetation type
    //            chil,         // leaf angle distribution factor
    //            ref(2,2),     // leaf reflectance (iw=iband, il=life and dead)
    //            tran(2,2),    // leaf transmittance (iw=iband, il=life and dead)
    //
    //          // time-space varying vegetation parameters
    //            green,        // green leaf fraction
    //            lai,          // leaf area index of exposed canopy (snow-free)
    //            sai             // stem area index
    //
    //// environmental variables
    //  real(r8), intent(in) :: &
    //            coszen,       // consine of solar zenith angle
    //            albg(2,2)       // albedos of ground
    //
    //// output
    //  real(r8), intent(out) :: &
    //            albv(2,2),    // albedo, vegetation [-]
    //            tranc(2,2),   // canopy transmittances for solar radiation
    //            thermk,       // canopy gap fraction for tir radiation
    //            extkb,        // (k, g(mu)/mu) direct solar extinction coefficient
    //            extkd,        // diffuse and scattered diffuse PAR extinction coefficient
    //            ssun(2,2),    // sunlit canopy absorption for solar radiation
    //            ssha(2,2)       // shaded canopy absorption for solar radiation,
    //                            // normalized by the incident flux

    //-------------------------- local -----------------------------------
    double
    phi1,          // (phi-1)
    phi2,          // (phi-2)
    scat,          // (omega)
    proj,          // (g(mu))
    zmu,           // (int(mu/g(mu))
    zmu2,          // (zmu * zmu)
    as,            // (a-s(mu))
    upscat,        // (omega-beta)
    betao,         // (beta-0)
    psi,           // (h)

    be,            // (b)
    ce,            // (c)
    de,            // (d)
    fe,            // (f)

    power1,        // (h*lai)
    power2,        // (k*lai)
    power3,        //

    sigma,         //
    s1,            //
    s2,            //
    p1,            //
    p2,            //
    p3,            //
    p4,            //
    f1,            //
    f2,            //
    h1,            //
    h4,            //
    m1,            //
    m2,            //
    m3,            //
    n1,            //
    n2,            //
    n3,            //

    hh1,           // (h1/sigma)
    hh2,           // (h2)
    hh3,           // (h3)
    hh4,           // (h4/sigma)
    hh5,           // (h5)
    hh6,           // (h6)
    hh7,           // (h7)
    hh8,           // (h8)
    hh9,           // (h9)
    hh10,          // (h10)

    eup[2][2],      // (integral of i_up*std::exp(-kx) )
    edown[2][2];       // (integral of i_down*std::exp(-kx) )

    int iw;                //

    //-----------------------------------------------------------------------
    // projected area of phytoelements in direction of mu and
    // average inverse diffuse optical depth per unit leaf area

    phi1 = 0.5 - 0.633 * chil - 0.33 * chil * chil;
    phi2 = 0.877 * ( 1. - 2. * phi1 );

    proj = phi1 + phi2 * coszen;
    extkb = (phi1 + phi2 * coszen) / coszen;

    extkd = 0.719;

    if ((fabs(phi1)>1.e-6)&&(fabs(phi2)>1.e-6))
        zmu = 1. / phi2 * ( 1.0 - phi1 / phi2 * std::log ( ( phi1 + phi2 ) / phi1 ) );
    else if ((fabs(phi1)<=1.e-6))
        zmu = 1./0.877;
    else if ((fabs(phi2)<=1.e-6))
        zmu = 1./(2.*phi1);

    zmu2 = zmu * zmu;

    power3 = (lai+sai) / zmu;
    power3 = min( 50., power3 );
    power3 = max( 1.e-5, power3 );
    thermk = std::exp(-power3);

    //do iw = 1, 2    // WAVE_BAND_LOOP
    for (iw=0; iw<2; iw++)
    {

        //-----------------------------------------------------------------------
        //     calculate average scattering coefficient, leaf projection and
        //     other coefficients for two-stream model.
        //-----------------------------------------------------------------------

        scat = green * (tran[iw][0] + ref[iw][0] ) +
               ( 1. - green ) * ( tran[iw][1] + ref[iw][1] );

        as = scat / 2. * proj / ( proj + coszen * phi2 ) ;
        as = as * ( 1. - coszen * phi1 / ( proj + coszen * phi2 ) *
                    std::log ( ( proj + coszen * phi2 + coszen * phi1 ) / ( coszen * phi1 ) ) );

        upscat = green * tran[iw][0] + ( 1.0 - green ) * tran[iw][1]  ;
        upscat = 0.5 * ( scat + ( scat - 2. * upscat ) * pow((( 1. - chil ) / 2.0), 2));
        betao = ( 1. + zmu * extkb ) / ( scat * zmu * extkb ) * as ;

        //-----------------------------------------------------------------------
        //     intermediate variables identified in appendix of SE-85.
        //-----------------------------------------------------------------------

        be = 1. - scat + upscat ;
        ce = upscat ;
        de = scat * zmu * extkb * betao;
        fe = scat * zmu * extkb * ( 1. - betao ) ;

        psi = std::sqrt(be*be - ce*ce)/zmu;
        power1 = min( psi*lai, 50. ) ;
        power2 = min( extkb*lai, 50. ) ;
        s1 = std::exp( - power1 );
        s2 = std::exp( - power2 )  ;

        //-----------------------------------------------------------------------
        //     calculation of direct albedos and canopy transmittances.
        //     albv [iw][0]     ( i-up )
        //     tranc(iw,irad)  ( i-down )
        //-----------------------------------------------------------------------

        p1 = be + zmu * psi;
        p2 = be - zmu * psi;
        p3 = be + zmu * extkb;
        p4 = be - zmu * extkb;

        f1 = 1. - albg[iw][1]*p1/ce;
        f2 = 1. - albg[iw][1]*p2/ce;

        h1 = - ( de * p4 + ce * fe );
        h4 = - ( fe * p3 + ce * de );

        sigma = ( zmu * extkb ) *( zmu * extkb ) + ( ce*ce - be*be )  ;

        if (fabs(sigma)>1.e-10)    //<======
        {
            hh1 = h1 / sigma;
            hh4 = h4 / sigma;

            m1 = f1 * s1;
            m2 = f2 / s1;
            m3 = ( albg[iw][0] - ( hh1 - albg[iw][1] * hh4 ) ) * s2;

            n1 = p1 / ce;
            n2 = p2 / ce;
            n3 = - hh4;

            hh2 = (m3*n2 - m2*n3) / (m1*n2 - m2*n1);
            hh3 = (m3*n1 - m1*n3) / (m2*n1 - m1*n2);

            hh5 = hh2 * p1 / ce;
            hh6 = hh3 * p2 / ce;

            albv[iw][0] = hh1 + hh2 + hh3;
            tranc[iw][0] = hh4 * s2 + hh5 * s1 + hh6 / s1  ;

            eup[iw][0] = hh1 * (1. - s2*s2) / (2.*extkb)
                         + hh2 * (1. - s1*s2) / (extkb + psi)
                         + hh3 * (1. - s2/s1) / (extkb - psi);

            edown[iw][0] = hh4 * (1. - s2*s2) / (2.*extkb)
                           + hh5 * (1. - s1*s2) / (extkb + psi)
                           + hh6 * (1. - s2/s1) / (extkb - psi);
        }

        else                               //<======
        {
            m1 = f1 * s1;
            m2 = f2 / s1;
            m3 = h1 / zmu2 * ( lai + 1. / (2.*extkb) ) * s2
                 + albg[iw][1] / ce * ( - h1 / (2.*extkb) / zmu2 *
                                        ( p3*lai + p4 / (2.*extkb) ) - de ) * s2
                 + albg[iw][0] * s2;

            n1 = p1 / ce;
            n2 = p2 / ce;
            n3 = 1./ce * ( h1*p4 / (4.*extkb*extkb) / zmu2 + de);

            hh2 = (m3*n2 - m2*n3) / (m1*n2 - m2*n1);
            hh3 = (m3*n1 - m1*n3) / (m2*n1 - m1*n2);

            hh5 = hh2 * p1 / ce;
            hh6 = hh3 * p2 / ce;

            albv[iw][0] =  - h1 / (2.*extkb*zmu2) + hh2 + hh3;
            tranc[iw][0] = 1./ce * ( -h1 / (2.*extkb*zmu2) *
                                     ( p3*lai + p4 / (2.*extkb) ) - de ) * s2
                           + hh5 * s1 + hh6 / s1;

            eup[iw][0] = (hh2 - h1/(2.*extkb*zmu2)) * (1. - s2*s2) / (2.*extkb)
                         + hh3 * (lai - 0.)
                         + h1/(2.*extkb*zmu2) * ( lai*s2*s2 - (1. - s2*s2)/(2.*extkb) );

            edown[iw][0] = (hh5 - (h1*p4/(4.*extkb*extkb*zmu) + de)/ce) *
                           (1. - s2*s2) / (2.*extkb)
                           + hh6 * (lai - 0.)
                           + h1*p3/(ce*4.*extkb*extkb*zmu2) *
                           ( lai*s2*s2 - (1. - s2*s2)/(2.*extkb) );

        }                              //<======

        ssun[iw][0] = (1.-scat) * ( 1.-s2 + 1. / zmu * (eup[iw][0] + edown[iw][0]) ) ;
        ssha[iw][0] = scat * (1.-s2)
                      + ( albg[iw][1]*tranc[iw][0] + albg[iw][0]*s2 - tranc[iw][0] ) - albv[iw][0]
                      - ( 1. - scat ) / zmu * ( eup[iw][0] + edown[iw][0] ) ;

        //-----------------------------------------------------------------------
        //     calculation of diffuse albedos and canopy transmittances
        //     albv [iw][1] ( i-up )
        //     tranc[iw][1] ( i-down )
        //-----------------------------------------------------------------------

        m1 = f1 * s1;
        m2 = f2 / s1;
        m3 = 0.;

        n1 = p1 / ce;
        n2 = p2 / ce;
        n3 = 1.;

        hh7 = -m2 / (m1*n2 - m2*n1);
        hh8 = -m1 / (m2*n1 - m1*n2);

        hh9 = hh7 * p1 / ce;
        hh10 = hh8 * p2 / ce;

        albv[iw][1] =  hh7 + hh8 ;
        tranc[iw][1] = hh9 * s1 + hh10 / s1;

        if (fabs(sigma)>1.e-10)
        {
            eup[iw][1]   = hh7 * (1. - s1*s2) / (extkb + psi)
                           + hh8 * (1. - s2/s1) / (extkb - psi);
            edown[iw][1] = hh9 * (1. - s1*s2) / (extkb + psi)
                           + hh10 * (1. - s2/s1) / (extkb - psi);
        }
        else
        {
            eup[iw][1]   = hh7 * (1. - s1*s2) / ( extkb + psi) + hh8 * (lai - 0.);
            edown[iw][1] = hh9 * (1. - s1*s2) / ( extkb + psi) + hh10 * (lai - 0.);
        }
        ssun[iw][1] = (1.-scat) / zmu * (eup[iw][1] + edown[iw][1]);
        ssha[iw][1] = tranc[iw][1] * ( albg[iw][1] -1. ) - ( albv[iw][1] - 1. )
                      - ( 1. - scat ) / zmu * ( eup[iw][1] + edown[iw][1] ) ;

        //enddo           // WAVE_BAND_LOOP
    }
}
void CommonLandModel::vec2xy(int  lat_points,          int    lon_points,           int  numpatch,                   int ixy_patch[MAXPATCH],
                             int  jxy_patch[MAXPATCH], double wtxy_patch[MAXPATCH], int itypwat[MAXPATCH],           int  nfcon,
                             int  nforc,               int    nfldv,                double  fcon[MAXPATCH][NFCON],   double forcxy[MAXLONPOINT][MAXLATPOINT][NFORC],
                             double fldv[MAXPATCH][NFLDV],                          double fldxy_r[MAXLONPOINT][MAXLATPOINT][NFLDV])
{

    // ----------------------------------------------------------------------
    // perfrom grid-average from subgrid 1d vector
    // subgrid to grid average mapping: average a subgrid input vector [fldv]
    // of length to a 2-d [lon_points] x [lat_points] output array [fldxy_r]
    //
    // Created by Yongjiu Dai
    //--------------//-------------------------------------------------------
    // 01: taux     // wind stress: E-W [kg/m/s2]
    // 02: tauy     // wind stress: N-S [kg/m/s2]
    // 03: fsena    // sensible heat from canopy height to atmosphere [W/m2]
    // 04: lfevpa   // latent heat flux from canopy height to atmosphere [W/m2]
    // 05: fevpa    // evapotranspiration from canopy to atmosphere [mm/s]
    // 06: fsenl    // sensible heat from leaves [W/m2]
    // 07: fevpl    // evaporation+transpiration from leaves [mm/s]
    // 08: etr      // transpiration rate [mm/s]
    // 09: fseng    // sensible heat flux from ground [W/m2]
    // 10: fevpg    // evaporation heat flux from ground [mm/s]
    // 11: fgrnd    // ground heat flux [W/m2]
    // 12: sabvsun  // solar absorbed by sunlit canopy [W/m2]
    // 13: sabvsha  // solar absorbed by shaded [W/m2]
    // 14: sabg     // solar absorbed by ground  [W/m2]
    // 15: olrg     // outgoing long-wave radiation from ground+canopy [W/m2]
    // 16: rnet     // net radiation [W/m2]
    // 17: xerr     // the error of water banace [mm/s]
    // 18: zerr     // the error of energy balance [W/m2]
    // 19: rsur     // surface runoff [mm/s]
    // 20: rnof     // total runoff [mm/s]
    // 21: assim    // canopy assimilation rate [mol m-2 s-1]
    // 22: respc    // respiration (plant+soil) [mol m-2 s-1]
    //--------------//-------------------------------------------------------
    // 23:32: tss   // soil temperature [K]
    // 33:42: wliq  // liquid water in soil layers [kg/m2]
    // 43:52: wice  // ice lens in soil layers [kg/m2]
    //--------------//-------------------------------------------------------
    // 53: tg       // ground surface temperature [K]
    // 54: tlsun    // sunlit leaf temperature [K]
    // 55: tlsha    // shaded leaf temperature [K]
    // 56: ldew     // depth of water on foliage [mm]
    // 57: scv      // snow cover, water equivalent [mm]
    // 58: snowdp   // snow depth [meter]
    // 59: fsno     // fraction of snow cover on ground
    // 60: sigf     // fraction of veg cover, excluding snow-covered veg [-]
    // 61: green    // leaf greenness
    // 62: lai      // leaf area index
    // 63: sai      // stem area index
    // 64: alb(1,1) // averaged albedo [visible, direct]
    // 65: alb(1,2) // averaged albedo [visible, diffuse]
    // 66: alb(2,1) // averaged albedo [near-infrared, direct]
    // 67: alb(2,2) // averaged albedo [near-infrared,diffuse]
    // 68: emis     // averaged bulk surface emissivity
    // 69: z0ma     // effective roughness [m]
    //--------------//-------------------------------------------------------
    // 70: trad     // radiative temperature of surface [K]
    // 71: ustar    // u* in similarity theory [m/s]
    // 72: tstar    // t* in similarity theory [kg/kg]
    // 73: qstar    // q* in similarity theory [kg/kg]
    // 74: zol      // dimensionless height (z/L) used in Monin-Obukhov theory
    // 75: rib      // bulk Richardson number in surface layer
    // 76: fm       // integral of profile function for momentum
    // 77: fh       // integral of profile function for heat
    // 78: fq       // integral of profile function for moisture
    //--------------//-------------------------------------------------------
    // 79: tref     // 2 m height air temperature [kelvin]
    // 80: qref     // 2 m height air specific humidity [kg/kg]
    // 81: u10m     // 10m u-velocity [m/s]
    // 82: v10m     // 10m v-velocity [m/s]
    // 83: f10m     // integral of profile function for momentum at 10m [-]
    //--------------//-------------------------------------------------------
    // 84: us       // wind in eastward direction [m/s]
    // 85: vs       // wind in northward direction [m/s]
    // 86: tm       // temperature at reference height [kelvin]
    // 87: qm       // specific humidity at reference height [kg/kg]
    // 88: prc      // convective precipitation [mm/s]
    // 89: prl      // large scale precipitation [mm/s]
    // 90: pbot     // atmospheric pressure at the surface [pa]
    // 91: frl      // atmospheric infrared (longwave) radiation [W/m2]
    // 92: solar    // downward solar radiation at surface [W/m2]
    //--------------//-------------------------------------------------------

    // use precision
    // use phycon_module, only: vonkar, stefnc, cpair, rgas, grav
    // implicit none

    // arguments:
    //integer, intent(in) :: lon_points            // number of longitude points on model grid
    //integer, intent(in) :: lat_points            // number of latitude points on model grid

    //integer, intent(in) :: numpatch              // total number of patches of grids
    //integer, intent(in) :: ixy_patch(numpatch)   // patch longitude index
    //integer, intent(in) :: jxy_patch(numpatch)   // patch latitude index
    //integer, intent(in) :: itypwat(numpatch)     // land water type
    //real(r8), intent(in) :: wtxy_patch(numpatch) // patch weight

    //integer, INTENT(in) :: nfcon                 // number of time constant variables
    //integer, INTENT(in) :: nforc                 // number of forcing variables
    //integer, intent(in) :: nfldv                 // number of output variables
    //real(r8), INTENT(in) :: fcon(numpatch,nfcon) // time constant variables
    //real(r8), intent(in) :: fldv(numpatch,nfldv) // output variables
    //real(r8), intent(in) :: forcxy(lon_points,lat_points,nforc) // xy gridded forcing
    //real(r8), intent(out) :: fldxy_r(lon_points,lat_points,nfldv) // xy gridded output

    // local variables
    int  i,j,k,L,tt;                            // indices
    double  a;                                  //
    double  sumwt[MAXLONPOINT][MAXLATPOINT][NFLDV]; // sum of wt
    double  maxwt[MAXLONPOINT][MAXLATPOINT];       // maximum of patches
    int     maxfp[MAXLONPOINT][MAXLATPOINT];        // index of maximum patches
    double  rhoair,thm,th,thv,ur,displa,zldis,hu,ht,hq;
    double  z0m,z0h,z0q,us,vs,tm,qm,pbot,psrf;
    double  obu,temp1,temp2,temp12m,temp22m;
    double  um,thvstar,beta,zii,wc,wc2;
    int ivt;

    //----------------------------------------------------------------------



    for (i = 0; i<NFLDV; i++)
    {
        for (j = 0; j<lat_points; j++)
            for (k=0; k<lon_points; j++)
            {
                fldxy_r[k][j][i] = 0.0;
                sumwt[k][j][i]   = 0.0;
                maxwt[k][j]      = 0.0;
            }
    }

    for (k=0; k<numpatch; k++)
    {
        i = ixy_patch[k];
        j = jxy_patch[k];
        maxfp[i][j] = k ;
    }

    // Get the number of patch with the largest fraction of grid within 1d clm array
    for (k=0; k<numpatch; k++)
    {
        i = ixy_patch[k];
        j = jxy_patch[k];
        a = maxwt[i][j];
        maxwt[i][j] = max(a,wtxy_patch[k]);
        if ((maxwt[i][j]-a)>0.0)
            maxfp[i][j] = k ;
    }

    // Mapping the 1d [numpatch] to grid [lon_points x lat_points]
    for (L = 1; L<=69; L++)
    {

        for (k=0; k<numpatch; k++)
        {
            i = ixy_patch[k];
            j = jxy_patch[k];

            if (L<=22)            // fluxes
                //[1-22] Grid averages by area-weight over grid patches
            {
                fldxy_r[i][j][L-1] = fldxy_r[i][j][L-1] + wtxy_patch[k]*fldv[k][L-1];
                sumwt[i][j][L-1] = sumwt[i][j][L-1] + wtxy_patch[k];
            }
            else if (L<=52)      // soil temperature and water
                //[23-32 33-42 43-52] Area-weight over grid patches but excluding lake and ocean patches
            {
                if (itypwat[k]<=3) // lake and ocean excluded
                {
                    fldxy_r[i][j][L-1] = fldxy_r[i][j][L-1] + wtxy_patch[k]*fldv[k][L-1];
                    sumwt[i][j][L-1] = sumwt[i][j][L-1] + wtxy_patch[k];
                }
                else                        // clm state variables
                    //[53-69] Grid averages by area-weight over grid patches
                {
                    fldxy_r[i][j][L-1] = fldxy_r[i][j][L-1] + wtxy_patch[k]*fldv[k][L-1];
                    sumwt[i][j][L-1] = sumwt[i][j][L-1] + wtxy_patch[k];
                }
            }
        }
    }

    for (L = 1; L<=69; L++)
        for (j=0; j<lat_points; j++)
            for (i=1; i<lon_points; i++)
            {
                //              if(sumwt[i][j][L-1].gt.1.0 .or. sumwt[i][j][L-1].lt.0.0)then
                //                 write(6,*) 'summation of fraction patches = ', sumwt[i][j][L-1], i,j,j
                //                 call abort
                //              endif
                if (sumwt[i][j][L-1]>0.0)
                    fldxy_r[i][j][L-1] = fldxy_r[i][j][L-1]/sumwt[i][j][L-1];
                else
                    fldxy_r[i][j][L-1] = -9999.0;
            }


    //[70-78] Retrieve through averaged fluxes
    //     do L = 70, 78
    //        do k = 1, numpatch
    //           i = ixy_patch[k]
    //           j = jxy_patch[k]
    //           if(k.eq.maxfp(i,j))then  // take values as that at the largest patches
    //                 fldxy_r[i][j][L-1] = fldv[k][L-1]
    //                 sumwt[i][j][L-1] = 1.0
    //           endif
    //        enddo
    //     enddo

    for (L = 70; L<=NFLDV; L++)
    {
        for (k=1; k<numpatch; k++)
        {
            i = ixy_patch[k];
            j = jxy_patch[k];
            sumwt[i][j][L-1] = sumwt[i][j][L-1] + wtxy_patch[k];
        }
    }

    for (j=0; j<lat_points; j++)
        for (i=1; i<lon_points; i++)
        {
            if (sumwt[i][j][69]>0.0)        //For land only defined
            {
                z0m = fldxy_r[i][j][68];
                z0h = fldxy_r[i][j][69];
                z0q = fldxy_r[i][j][69];
                displa = 2./3.*z0m/0.07;
                hu = max(forcxy[i][j][15],5.+displa);
                ht = max(forcxy[i][j][16],5.+displa);
                hq = max(forcxy[i][j][17],5.+displa);
                zldis = hu-displa;
                us = forcxy[i][j][2];
                vs = forcxy[i][j][3];
                tm = forcxy[i][j][4];
                qm = forcxy[i][j][5];
                pbot = forcxy[i][j][8];
                psrf = forcxy[i][j][9];

                rhoair = (pbot-0.378*qm*pbot/(0.622+0.378*qm))/(rgas*tm);

                fldxy_r[i][j][69] = pow((fldxy_r[i][j][14]/stefnc),0.25);
                fldxy_r[i][j][70] = std::sqrt(max(1.e-6,std::sqrt(fldxy_r[i][j][0]*fldxy_r[i][j][0]+fldxy_r[i][j][1]*fldxy_r[i][j][1]))/rhoair);
                fldxy_r[i][j][71] = -fldxy_r[i][j][2]/(rhoair*fldxy_r[i][j][70])/cpair;
                fldxy_r[i][j][72] = -fldxy_r[i][j][4]/(rhoair*fldxy_r[i][j][70]);

                thm = tm + 0.0098*ht;
                th = tm*pow((100000./psrf),(rgas/cpair));
                thv = th*(1.+0.61*qm);

                fldxy_r[i][j][73] =   zldis*vonkar*grav
                                      * (fldxy_r[i][j][71]+0.61*th*fldxy_r[i][j][72])
                                      / (fldxy_r[i][j][70]*fldxy_r[i][j][70]*thv);

                if (fldxy_r[i][j][73]>=0.0)  //stable
                    fldxy_r[i][j][73] = min(2.,max(fldxy_r[i][j][73],1.e-6));
                else                              //unstable
                    fldxy_r[i][j][73] = max(-100.,min(fldxy_r[i][j][73],-1.e-6));

                beta = 1.0;
                zii = 1000.0;
                thvstar=fldxy_r[i][j][71]+0.61*th*fldxy_r[i][j][71];
                ur = std::sqrt(us*us+vs*vs);
                if (fldxy_r[i][j][73]>=0.0)
                    um = max(ur,0.1);
                else
                {
                    wc = pow((-grav*fldxy_r[i][j][70]*thvstar*zii/thv),(1.0/3.0));
                    wc2 = beta*beta*(wc*wc);
                    um = max(0.1,std::sqrt(ur*ur+wc2));
                }
                obu = zldis/fldxy_r[i][j][73];

                moninobuk(hu,ht,hq,displa,z0m,z0h,z0q,
                          obu,um,fldxy_r[i][j][70],temp1,temp2,temp12m,temp22m,
                          fldxy_r[i][j][82],fldxy_r[i][j][75],fldxy_r[i][j][76],fldxy_r[i][j][77]);

                fldxy_r[i][j][74] = fldxy_r[i][j][73]*pow(vonkar,3)*pow(fldxy_r[i][j][73],2)/(temp1*um*um);
                fldxy_r[i][j][74] = min(5.,fldxy_r[i][j][74]);
            }

            else
            {
                for (tt= 70; tt<=78; tt++)
                    fldxy_r[i][j][tt-1] = -9999.0;
            }
        }


    //[79-83] Grass land only (for matching the routine meteorological obs.)
    for (L=79; L<=83; L++)
    {
        for (k=0; k<numpatch; k++)
        {
            i = ixy_patch[k];
            j = jxy_patch[k];
            ivt = (int)(fcon[k][4]);
            if ((ivt==7)||(ivt==0))
                fldxy_r[i][j][L-1] = fldv[k][L-1];
        }
    }

    //[84-92] Meteorological forcing
    for (j=0; j<lat_points; j++)
        for (i=0; i<lon_points; i++)
        {
            if (sumwt[i][j][83]>=0.0)     //For land only defined
            {
                fldxy_r[i][j][83] = forcxy[i][j][2];
                fldxy_r[i][j][84] = forcxy[i][j][3];
                fldxy_r[i][j][85] = forcxy[i][j][4];
                fldxy_r[i][j][86] = forcxy[i][j][5];
                fldxy_r[i][j][87] = forcxy[i][j][6];
                fldxy_r[i][j][88] = forcxy[i][j][7];
                fldxy_r[i][j][89] = forcxy[i][j][8];
                fldxy_r[i][j][90] = forcxy[i][j][14];
                fldxy_r[i][j][90] = forcxy[i][j][10]+forcxy[i][j][11]+forcxy[i][j][12]+forcxy[i][j][13];
            }
            else
            {
                for (tt= 84; tt<=92; tt++)
                    fldxy_r[i][j][tt-1] = -9999.0;
            }
        }


}
void CommonLandModel::WATER (int    itypwat ,           int    lb,                    int    nl_soil ,             double dtime  ,               double z[MAXSNL+MAXSOILL],
                             double dz[MAXSNL+MAXSOILL],double zi[MAXSNL+MAXSOILL+1], double bsw[MAXSOILL] ,       double porsl[MAXSOILL]  ,     double phi0[MAXSOILL] ,
                             double hksati[MAXSOILL] ,  double rootr[MAXSOILL]   ,    double tss[MAXSNL+MAXSOILL] ,double wliq[MAXSNL+MAXSOILL] ,double  wice[MAXSNL+MAXSOILL] ,
                             double pg_rain,            double sm      ,              double etr    ,              double qseva  ,               double  qsdew   ,
                             double qsubl  ,            double qfros   ,              double &rsur   ,              double &rnof   ,               double   wtfact  ,
                             double pondmx ,            double ssi     ,              double wimp   ,              double smpmin  )
{
    //=======================================================================
    // this is the main subroutine to execute the calculation of
    // hydrological processes
    //
    // Original author : Yongjiu Dai, 09/15/1999; 08/30/2002
    //
    // FLOW DIAGRAM FOR WATER.F90
    //
    // WATER ===> snowwater
    //            surfacerunoff
    //            soilwater
    //            subsurfacerunoff
    //
    //=======================================================================

    // use precision
    // use phycon_module, only : denice, denh2o, tfrz
    // implicit none

    //-----------------------Argument---------- ------------------------------
    //integer, INTENT(in) :: &
    //      lb               , // lower bound of array
    //      nl_soil          , // upper bound of array
    //      itypwat             // land water type (0=soil, 1=urban or built-up, 2=wetland,
    //                          // 3=land ice, 4=deep lake, 5=shallow lake)
    //real(r8), INTENT(in) :: &
    //      dtime            , // time step (s)
    //      wtfact           , // fraction of model area with high water table
    //      pondmx           , // ponding depth (mm)
    //      ssi              , // irreducible water saturation of snow
    //      wimp             , // water impremeable if porosity less than wimp
    //      smpmin           , // restriction for min of soil poten. (mm)
    //      z (lb:nl_soil)   , // layer depth (m)
    //      dz(lb:nl_soil)   , // layer thickness (m)
    //      zi(lb-1:nl_soil) , // interface level below a "z" level (m)

    //      bsw(1:nl_soil)   , // Clapp-Hornberger "B"
    //      porsl(1:nl_soil) , // saturated volumetric soil water content(porosity)
    //      phi0(1:nl_soil)  , // saturated soil suction (mm)
    //      hksati(1:nl_soil), // hydraulic conductivity at saturation (mm h2o/s)
    //      rootr(1:nl_soil) , // root resistance of a layer, all layers add to 1.0

    //      tss(lb:nl_soil)  , // soil/snow skin temperature (K)
    //      pg_rain          , // rainfall after removal of interception (mm h2o/s)
    //      sm               , // snow melt (mm h2o/s)
    //      etr              , // actual transpiration (mm h2o/s)
    //      qseva            , //ground surface evaporation rate (mm h2o/s)
    //      qsdew            , //ground surface dew formation (mm h2o /s) [+]
    //      qsubl            , //sublimation rate from snow pack (mm h2o /s) [+]
    //      qfros               //surface dew added to snow pack (mm h2o /s) [+]

    //real(r8), INTENT(inout) :: &
    //      wice(lb:nl_soil) , // ice lens (kg/m2)
    //      wliq(lb:nl_soil)    // liquid water (kg/m2)

    //real(r8), INTENT(out) :: &
    //      rsur             , // surface runoff (mm h2o/s)
    //      rnof                // total runoff (mm h2o/s)
    //
    //-----------------------Local Variables------------------------------
    //
    int i;                 // loop counter

    double
    eff_porosity[MAXSOILL], // effective porosity = porosity - vol_ice
                 hk[MAXSOILL]     , // hydraulic conductivity (mm h2o/s)
                 dhkdw[MAXSOILL]  , // d(hk)/d(vol_liq)
                 dwat[MAXSOILL]   , // change in soil water
                 gwat              , // net water input from top
                 qinfl             , // infiltration rate (mm h2o/s)
                 rsubst            , // subsurface runoff (mm h2o/s)
                 vol_liq[MAXSOILL], // partitial volume of liquid water in layer
                 vol_ice[MAXSOILL], // partitial volume of ice lens in layer
                 zmm [MAXSNL+MAXSOILL]   , // layer depth (mm)
                 dzmm[MAXSNL+MAXSOILL]   , // layer thickness (mm)
                 zimm[MAXSNL+MAXSOILL+1];      // interface level below a "z" level (mm)

    //=======================================================================
    // [1] update the liquid water within snow layer and the water onto soil
    //=======================================================================

    //for(int i=0;i<MAXSNL+MAXSOILL;i++)
    //	cout<<wliq[i]<<endl;
    //getchar();
    gwat=0;
    if (lb>=MAXSNL)
        gwat = pg_rain + sm - qseva;
    else
        snowwater(lb,dtime,ssi,wimp,pg_rain,
                  qseva,qsdew,qsubl,qfros,dz,
                  wice,wliq,gwat);
    //=======================================================================
    // [2] surface runoff and infiltration
    //=======================================================================
    if (itypwat<=1)  // soil ground only
    {
        // porosity of soil, partitial volume of ice and liquid
        for (i = 0; i<nl_soil; i++)
        {
            vol_ice[i] = min(porsl[i], wice[MAXSNL+i]/(dz[MAXSNL+i]*denice));
            eff_porosity[i] = porsl[i]-vol_ice[i];
            vol_liq[i] = min(eff_porosity[i], wliq[MAXSNL+i]/(dz[MAXSNL+i]*denh2o));
            //cout<<vol_ice[i]<<"  "<<eff_porosity[i]<<"  "<<vol_liq[i]<<endl;
        }

        // surface runoff including water table
        if (gwat > 0.0)
        {
            surfacerunoff (nl_soil,wtfact,wimp, bsw,porsl,
                           phi0,   hksati,z,            dz, zi,
                           vol_liq,vol_ice,eff_porosity,gwat,rsur);
            //cout<<rsur<<endl;
        }
        else
            rsur = 0.0;

        // infiltration into surface soil layer
        qinfl = gwat - rsur ;
        //=======================================================================
        // [3] determine the change of soil water
        //=======================================================================

        // convert length units from m to mm
        for (int k=0; k<MAXSNL+MAXSOILL; k++)
        {
            zmm[k] = z[k]*1000.0;
            dzmm[k] = dz[k]*1000.0;
            zimm[k] = zi[k]*1000.0;
        }
        zimm[MAXSNL+MAXSOILL] = zi[MAXSNL+MAXSOILL]*1000.0;
        //for(int i=0;i<MAXSNL+MAXSOILL;i++)
        //cout<<wliq[i]<<endl;
        soilwater (nl_soil,dtime,wimp,smpmin,porsl,
                   phi0,bsw,hksati,zmm,dzmm,
                   zimm,tss,vol_liq,vol_ice,eff_porosity,
                   qinfl,etr,rootr,dwat,hk,dhkdw);

        // update the mass of liquid water
        for (i=0; i<nl_soil; i++)
        {
            wliq[MAXSNL+i] = max(0.0,wliq[MAXSNL+i]+dwat[i]*dzmm[MAXSNL+i]);
            //to avoid the little value of wliq(wlx,20090917)
            if (wliq[MAXSNL+i]<1e-7) wliq[MAXSNL+i]=0;
            //wliq[MAXSNL+i] = min(wliq[MAXSNL+i],dz[MAXSNL+i]*eff_porosity[i]*denh2o);
            //cout<<MAXSNL+i<<"  "<<dwat[i]<<"  "<<wliq[MAXSNL+i]<<"  "<<hk[i]<<"  "<<dhkdw[i]<<endl;
        }

        //=======================================================================
        // [4] subsurface runoff and the corrections
        //=======================================================================

        subsurfacerunoff (nl_soil,dtime,pondmx,dzmm,wliq,
                          eff_porosity,hk,dhkdw,dwat,rsubst);

        // total runoff
        //cout<<rsubst<<"  "<<rsur<<endl;
        rnof = rsubst + rsur ;

        // renew the ice and liquid mass due to condensation
        if (lb >= MAXSNL)
        {
            wliq[MAXSNL] = wliq[MAXSNL] + qsdew*dtime;
            wice[MAXSNL] = wice[MAXSNL] + (qfros-qsubl)*dtime;
        }


        /*       for (i=0;i<nl_soil;i++)
               {
                   wliq[MAXSNL+i] = min(wliq[MAXSNL+i],dz[MAXSNL+i]*eff_porosity[i]*denh2o);
                   wliq[MAXSNL+i] = max(wliq[MAXSNL+i],0.0);

                   ice_poros      = eff_porosity[i] - wliq[MAXSNL+i]/(dz[MAXSNL+i]*denh2o);
                   wice[MAXSNL+i] = min(wice[MAXSNL+i],dz[MAXSNL+i]*ice_poros*denice);
                   wice[MAXSNL+i] = max(wice[MAXSNL+i],0.0);
                   //getchar();
               }
        */

    }
    //=======================================================================
    // [6] arbitrarily hydrological processes in wetland and glacier
    //=======================================================================

    else
    {
        if (itypwat==2)       // WETLAND
        {
            rsur=0.0;
            qinfl=gwat;
            rsubst=0.0;
            rnof=0.0;
            for (i=0; i<nl_soil; i++)
            {
                if (tss[MAXSNL+i]>=tfrz)
                {
                    wice[MAXSNL+i] = 0.0;
                    wliq[MAXSNL+i] = 1000.*dz[MAXSNL+i];
                }
                else
                {
                    wice[MAXSNL+i] = 1000.*dz[MAXSNL+i];
                    wliq[MAXSNL+i] = 0.0;
                }
            }
        }
        if (itypwat==3)       // LAND ICE
        {
            rsur=gwat;
            qinfl=0.;
            rsubst=0.;
            rnof=rsur;
            for (int k=0; k<nl_soil; k++)
            {
                wice[MAXSNL+k] = 1000.0*dz[MAXSNL+k];
                wliq[MAXSNL+k] = 0.0;
            }
        }
    }

    //-----------------------------------------------------------------------

}
}
