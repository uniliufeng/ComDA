#include "cldas.h"
#include "exception.h"
#include "random.hpp"
#include "radiancetransfermodel.h"
#include "qhmodel.h"
#include "enkf.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>

using namespace ldas;

Cldas_Observation::Cldas_Observation(const int r,const int c):m_row(r),m_col(c)
{
    //******************************************************
    //                 SSMI
    //******************************************************
    ssmi_data_19v = new short * [m_row];
    ssmi_data_19h = new short * [m_row];
    ssmi_data_37v = new short * [m_row];
    ssmi_data_37h = new short * [m_row];
    for (int i=0; i<m_row; i++)
    {
        ssmi_data_19v[i] = new short [m_col];
        ssmi_data_19h[i] = new short [m_col];
        ssmi_data_37v[i] = new short [m_col];
        ssmi_data_37h[i] = new short [m_col];
    }
    //******************************************************
    //                 AMSR
    //******************************************************
    amsr_data_6v   = new short * [m_row];
    amsr_data_6h   = new short * [m_row];
    amsr_data_10v  = new short * [m_row];
    amsr_data_10h  = new short * [m_row];
    amsr_data_18v  = new short * [m_row];
    amsr_data_18h  = new short * [m_row];
    amsr_data_37v  = new short * [m_row];
    amsr_data_37h  = new short * [m_row];
    for (int i=0; i<m_row; i++)
    {
        amsr_data_6v[i]   = new short [m_col];
        amsr_data_6h[i]   = new short [m_col];
        amsr_data_10v[i]  = new short [m_col];
        amsr_data_10h[i]  = new short [m_col];
        amsr_data_18v[i]  = new short [m_col];
        amsr_data_18h[i]  = new short [m_col];
        amsr_data_37v[i]  = new short [m_col];
        amsr_data_37h[i]  = new short [m_col];
    }
}
Cldas_Observation::~Cldas_Observation()
{
    for (int i=0; i<m_row; i++)
    {
        delete [] ssmi_data_19v[i];
        delete [] ssmi_data_19h[i];
        delete [] ssmi_data_37v[i];
        delete [] ssmi_data_37h[i];
    }
    delete [] ssmi_data_19v;
    delete [] ssmi_data_19h;
    delete [] ssmi_data_37v;
    delete [] ssmi_data_37h;

    for (int i=0; i<m_row; i++)
    {
        delete [] amsr_data_6v[i]  ;
        delete [] amsr_data_6h[i]  ;
        delete [] amsr_data_10v[i] ;
        delete [] amsr_data_10h[i] ;
        delete [] amsr_data_18v[i] ;
        delete [] amsr_data_18h[i] ;
        delete [] amsr_data_37v[i] ;
        delete [] amsr_data_37h[i] ;
    }
    delete [] amsr_data_6v   ;
    delete [] amsr_data_6h   ;
    delete [] amsr_data_10v  ;
    delete [] amsr_data_10h  ;
    delete [] amsr_data_18v  ;
    delete [] amsr_data_18h  ;
    delete [] amsr_data_37v  ;
    delete [] amsr_data_37h  ;
}
void Cldas_Observation::open(const std::string& basename,const DATA_TYPE dt)
{
    std::string fn;
    std::ifstream obs_file;
    switch(dt)
    {
    case AMSRE:
        fn=basename+".6V";
        obs_file.open(fn.c_str(),std::ios_base::binary);
        if(!obs_file)
            throw Exception("Cldas_Observation","open","can't open passive microwave data file:"+fn);
        for(int i=0; i<m_row; i++)
        {
            obs_file.read((char *) amsr_data_6v[i],sizeof(short int)*m_col);
            if (is_bigendian()) for(int j=0; j<m_col; j++)SWAP_16(amsr_data_6v[i][j]);
        }
        obs_file.close();
        fn=basename+".6H";
        obs_file.open(fn.c_str(),std::ios_base::binary);
        if(!obs_file)
            throw Exception("Cldas_Observation","open","can't open passive microwave data file:"+fn);
        for(int i=0; i<m_row; i++)
        {
            obs_file.read((char *) amsr_data_6h[i],sizeof(short int)*m_col);
            if (is_bigendian()) for(int j=0; j<m_col; j++)SWAP_16(amsr_data_6h[i][j]);
        }
        obs_file.close();
        fn=basename+".10V";
        obs_file.open(fn.c_str(),std::ios_base::binary);
        if(!obs_file)
            throw Exception("Cldas_Observation","open","can't open passive microwave data file:"+fn);
        for(int i=0; i<m_row; i++)
        {
            obs_file.read((char *) amsr_data_10v[i],sizeof(short int)*m_col);
            if (is_bigendian()) for(int j=0; j<m_col; j++)SWAP_16(amsr_data_10v[i][j]);
        }
        obs_file.close();
        fn=basename+".10H";
        obs_file.open(fn.c_str(),std::ios_base::binary);
        if(!obs_file)
            throw Exception("Cldas_Observation","open","can't open passive microwave data file:"+fn);
        for(int i=0; i<m_row; i++)
        {
            obs_file.read((char *) amsr_data_10h[i],sizeof(short int)*m_col);
            if (is_bigendian()) for(int j=0; j<m_col; j++)SWAP_16(amsr_data_10h[i][j]);
        }
        obs_file.close();
        fn=basename+".18V";
        obs_file.open(fn.c_str(),std::ios_base::binary);
        if(!obs_file)
            throw Exception("Cldas_Observation","open","can't open passive microwave data file:"+fn);
        for(int i=0; i<m_row; i++)
        {
            obs_file.read((char *) amsr_data_18v[i],sizeof(short int)*m_col);
            if (is_bigendian()) for(int j=0; j<m_col; j++)SWAP_16(amsr_data_18v[i][j]);
        }
        obs_file.close();
        fn=basename+".18H";
        obs_file.open(fn.c_str(),std::ios_base::binary);
        if(!obs_file)
            throw Exception("Cldas_Observation","open","can't open passive microwave data file:"+fn);
        for(int i=0; i<m_row; i++)
        {
            obs_file.read((char *) amsr_data_18h[i],sizeof(short int)*m_col);
            if (is_bigendian()) for(int j=0; j<m_col; j++)SWAP_16(amsr_data_18h[i][j]);
        }
        obs_file.close();
        fn=basename+".37V";
        obs_file.open(fn.c_str(),std::ios_base::binary);
        if(!obs_file)
            throw Exception("Cldas_Observation","open","can't open passive microwave data file:"+fn);
        for(int i=0; i<m_row; i++)
        {
            obs_file.read((char *) amsr_data_37v[i],sizeof(short int)*m_col);
            if (is_bigendian()) for(int j=0; j<m_col; j++)SWAP_16(amsr_data_37v[i][j]);
        }
        obs_file.close();
        fn=basename+".37H";
        obs_file.open(fn.c_str(),std::ios_base::binary);
        if(!obs_file)
            throw Exception("Cldas_Observation","open","can't open passive microwave data file:"+fn);
        for(int i=0; i<m_row; i++)
        {
            obs_file.read((char *) amsr_data_37h[i],sizeof(short int)*m_col);
            if (is_bigendian()) for(int j=0; j<m_col; j++)SWAP_16(amsr_data_37h[i][j]);
        }
        obs_file.close();
        break;
    case SSMI:
    	fn=basename+".37V";
        obs_file.open(fn.c_str(),std::ios_base::binary);
        if(!obs_file)
            throw Exception("Cldas_Observation","open","can't open passive microwave data file:"+fn);
        for(int i=0; i<m_row; i++)
        {
            obs_file.read((char *) ssmi_data_37v[i],sizeof(short int)*m_col);
            if (is_bigendian()) for(int j=0; j<m_col; j++)SWAP_16(ssmi_data_37v[i][j]);
        }
        obs_file.close();
    	fn=basename+".37H";
        obs_file.open(fn.c_str(),std::ios_base::binary);
        if(!obs_file)
            throw Exception("Cldas_Observation","open","can't open passive microwave data file:"+fn);
        for(int i=0; i<m_row; i++)
        {
            obs_file.read((char *) ssmi_data_37h[i],sizeof(short int)*m_col);
            if (is_bigendian()) for(int j=0; j<m_col; j++)SWAP_16(ssmi_data_37h[i][j]);
        }
        obs_file.close();
        fn=basename+".19V";
        obs_file.open(fn.c_str(),std::ios_base::binary);
        if(!obs_file)
            throw Exception("Cldas_Observation","open","can't open passive microwave data file:"+fn);
        for(int i=0; i<m_row; i++)
        {
            obs_file.read((char *) ssmi_data_19v[i],sizeof(short int)*m_col);
            if (is_bigendian()) for(int j=0; j<m_col; j++)SWAP_16(ssmi_data_19v[i][j]);
        }
        obs_file.close();
    	fn=basename+".19H";
        obs_file.open(fn.c_str(),std::ios_base::binary);
        if(!obs_file)
            throw Exception("Cldas_Observation","open","can't open passive microwave data file:"+fn);
        for(int i=0; i<m_row; i++)
        {
            obs_file.read((char *) ssmi_data_19h[i],sizeof(short int)*m_col);
            if (is_bigendian()) for(int j=0; j<m_col; j++)SWAP_16(ssmi_data_19h[i][j]);
        }
        obs_file.close();

        break;
    default:
        break;
    }
}
Cldas::Cldas()
{
    //ctor
}

Cldas::~Cldas()
{
    //dtor
    if (m_config.row_num>0)
		delete m_observe;
	if (m_ensemble>0)
	{
		for(int i=0;i<m_ensemble;i++)
		{
			for(int j=0;j<m_config.total_patches;j++)
			{
				delete []fvar_en[i][j];
				delete []flux_en[i][j];
			}
			delete []fvar_en[i];
			delete []flux_en[i];
		}
		delete []fvar_en;
		delete []flux_en;
	}
	if (m_config.total_patches>0)
	{
		for(int i=0;i<m_config.total_patches;i++)
		{
			delete []fvar_mean[i];
			delete []flux_mean[i];
			delete []snl_en[i];
		}
		delete []fvar_mean;
		delete []flux_mean;
		delete []snl_en;
	}
	if (m_config.row_num>0)
	{
		for(int i=0;i<m_config.row_num;i++)
		{
			delete []sand_grid[i];
			delete []clay_grid[i];
		}
		delete []sand_grid;
		delete []clay_grid;
	}
}

Cldas::Cldas(const Cldas& other)
{
    //copy ctor
}

Cldas& Cldas::operator=(const Cldas& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}

void Cldas::config(const char* fn)
{
	CommonLandModel::config(fn);
	m_observe=new Cldas_Observation(m_config.row_num,m_config.col_num);
	m_ensemble=m_config.ensemble_size;
	fvar_mean=new double*[m_config.total_patches];
	flux_mean=new double*[m_config.total_patches];
	snl_en=new int*[m_config.total_patches];
	for(int i=0;i<m_config.total_patches;i++)
	{
		fvar_mean[i]=new double[NFVAR];
		snl_en[i]=new int[m_ensemble];
		flux_mean[i]=new double[OUTFLUX];
	}

	fvar_en=new double**[m_ensemble];
	flux_en=new double**[m_ensemble];

	sand_grid=new float*[m_config.row_num];
	clay_grid=new float*[m_config.row_num];
	for(int i=0;i<m_config.row_num;i++)
	{
		sand_grid[i]=new float[m_config.col_num];
		clay_grid[i]=new float[m_config.col_num];
	}
	std::ifstream sandfile(m_config.sand_path.c_str());
	for(int i=0;i<m_config.row_num;i++)
	{
		sandfile.read((char*)sand_grid[i],sizeof(float)*m_config.col_num);
		if (is_bigendian()) for(int j=0; j<m_config.col_num; j++) SWAP_FLOAT(sand_grid[i][j]);
	}
	sandfile.close();
	std::ifstream clayfile(m_config.clay_path.c_str());
	for(int i=0;i<m_config.row_num;i++)
	{
		clayfile.read((char*)clay_grid[i],sizeof(float)*m_config.col_num);
		if (is_bigendian()) for(int j=0; j<m_config.col_num; j++) SWAP_FLOAT(clay_grid[i][j]);
	}
	clayfile.close();

	for(int i=0;i<m_ensemble;i++)
	{
		fvar_en[i]=new double*[m_config.total_patches];
		flux_en[i]=new double*[m_config.total_patches];
		for(int j=0;j<m_config.total_patches;j++)
		{
			fvar_en[i][j]=new double[NFVAR];
			flux_en[i][j]=new double[OUTFLUX];
            for(int k=0; k<NFVAR; k++)
            {
                fvar_en[i][j][k]=parameter.fvar[j][k];
            }
		}
	}
}
void Cldas::run()
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
        sst.clear();
        sst.str("");
        sst<<m_config.forcing_path<<ym<<"/"<<ydh<<".forc";
        forcing.open(sst.str().c_str());
        if (month!=parameter.time.getMonth())
        {
            sst.clear();
            sst.str("");
            sst<<m_config.vegetation_path<<"LAI_"<<ym<<".900s";
            parameter.open(sst.str().c_str(),MODIS_LAI,m_config.row_num*m_config.col_num);
            month=parameter.time.getMonth();
        }
        // Calendar for NEXT time step
        parameter.time.addSeconds(m_config.time_step);

        //********************************************************
        // step 2 if obseravtion exist, read observation data
        //********************************************************

        //********************************************************
        // step 2.1  set observation file name
        //********************************************************
        sst.clear();
        sst.str("");
        sst<<m_config.ssmi_path<<"ssmi"<<ydh<<".19V";
        ifstream ssmi_in(sst.str().c_str());
        bool obs_ssmi_flag=ssmi_in;
        ssmi_in.close();

        sst.clear();
        sst.str("");
        sst<<m_config.amsr_path<<"AMSR"<<ydh<<".6V";
        ifstream amsr_in(sst.str().c_str());
        bool obs_amsr_flag=amsr_in;
        amsr_in.close();

        // no passive microwave observation
        obs_flag = false;
        // have ssmi no amsr
        if (obs_ssmi_flag && !obs_amsr_flag)
        {
            obs_flag = true;
            sensor = 1;
        }
        else if (!obs_ssmi_flag && obs_amsr_flag)
            // have amsr no ssmi
        {
            obs_flag = true;
            sensor = 2;
        }
        else if (obs_ssmi_flag && obs_amsr_flag)
            // have ssmi and amsr
        {
            obs_flag = true;
            sensor = 3;
        }

        //********************************************************
        // step 2.2 open and read current observation data
        //********************************************************
        if (obs_flag)
        {
            std::cout<<"Observe available in current time:"<<t_time<<std::endl;
            /*
            if (sensor==1)
                std::cout<<"SSMI data available in current time:  "<<date_int<<std::endl;
            if (sensor==2)
                std::cout<<"AMSR data available in current time:  "<<date_int<<std::endl;
            if (sensor==3)
                std::cout<<"SSMI and AMSR data available in current time:  "<<date_int<<std::endl;
            */
            //read ssmi data
            if ((sensor==1)||(sensor==3))
            {
            	sst.clear();
            	sst.str("");
            	sst<<m_config.ssmi_path<<"ssmi"<<ydh;
            	m_observe->open(sst.str(),SSMI);
            }
            //read amsr data
            if ((sensor==2)||(sensor==3))
            {
            	sst.clear();
            	sst.str("");
            	sst<<m_config.amsr_path<<"AMSR"<<ydh;
            	m_observe->open(sst.str(),AMSRE);
            }
        }

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
                    fvar_out[ipatch][i] = fvar_mean[ipatch][30+i];
                // output water content
                for (int i=0; i<15; i++)
                    fvar_out[ipatch][i+15] = fvar_mean[ipatch][45+i];
                // output ice content
                for (int i=0; i<15; i++)
                    fvar_out[ipatch][i+30] = fvar_mean[ipatch][60+i];

                // ground surface temperature
                fvar_out[ipatch][45] = fvar_mean[ipatch][75];
                // snow cover
                fvar_out[ipatch][46] = fvar_mean[ipatch][80];
                // snow depth
                fvar_out[ipatch][47] = fvar_mean[ipatch][81];
            }
            sst.str("");
            sst<<m_config.output_path<<ym<<"/"<<ydh<<".grid";
            output(sst.str().c_str(),VARNUM,fvar_out);
            // output flux
            sst.str("");
            sst<<m_config.output_path<<ym<<"/"<<ydh<<".flux";
            output(sst.str().c_str(),FLUXNUM,flux_mean);
        }

        //*************************************************************************
        // step 5    write restart file
        //*************************************************************************
        if ((istep%m_config.restart_span==0)||(istep==m_config.total_steps))
        {
            std::stringstream ss;
            ss<<m_config.restart_path<<ydh<<".patch";
            for(int i=0;i<m_config.total_patches;i++)
			{
				for(int j=0;j<NFVAR;j++)
				{
					parameter.fvar[i][j]=fvar_mean[i][j];
				}
			}
            parameter.save(ss.str().c_str());
        }
        //*************************************************************************
        cout<<"caculation is finished in current time"<<endl;
        //*************************************************************************
    }
}

void Cldas::run(const int ipatch)
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
    for(int ien=0;ien<m_ensemble;ien++)
		fvar_en[ien][ipatch][86]=parameter.lai[location];
    //********************************************************************
    // step 3.3  calculate forecast variable of patch
    //********************************************************************
    double air_tmp=forc_data[4];
    double rain=forc_data[6];
    Random random(std::time(NULL));
    for(int i=0;i<OUTFLUX;i++)
		flux_mean[ipatch][i]=0;
	for(int i=0;i<NFVAR;i++)
		fvar_mean[ipatch][i]=0;
    for(int ien=0;ien<m_ensemble;ien++)
	{
		forc_data[4]=air_tmp+random.normt_rnd(0,1,-1,1)*3;
		forc_data[6]=rain*(1 + 0.10*random.normt_rnd(0,1,-1,1));
        CommonLandModel::run(constant.ftune,constant.fcon[ipatch],forc_data,fvar_en[ien][ipatch],flux_en[ien][ipatch]);
        snl_en[ipatch][ien]=snowlayer();
        for(int i=0; i<OUTFLUX; i++)
            flux_mean[ipatch][i]+=flux_en[ien][ipatch][i];
		for(int i=0;i<NFVAR;i++)
			fvar_mean[ipatch][i]+=fvar_en[ien][ipatch][i];
    }
    for(int i=0;i<OUTFLUX;i++)
		flux_mean[ipatch][i]/=m_ensemble;
	for(int i=0;i<NFVAR;i++)
		fvar_mean[ipatch][i]/=m_ensemble;

	//进行观测、同化
	bool ass_flag=(obs_flag && !(parameter.lai[location]>1.5));
	int row=grid.row(ipatch);
	int column=grid.column(ipatch);
	int sensor_type=sensor;
	if (ass_flag)
	{
		switch(sensor_type)
		{
		case 1:
			ass_flag=(m_observe->ssmi_data_19h[row][column]>0);
			break;
		case 2:
			ass_flag=(m_observe->amsr_data_6h[row][column]>0);
			break;
		case 3:
			ass_flag=(m_observe->ssmi_data_19h[row][column]>0 || m_observe->amsr_data_6h[row][column]>0);
			if (m_observe->amsr_data_6h[row][column]<=0 && m_observe->ssmi_data_19h[row][column]>0)
				sensor_type=1;
			if (m_observe->amsr_data_6h[row][column]>0 && m_observe->ssmi_data_19h[row][column]<=0)
				sensor_type=2;
			break;
		default:
			break;
		}
	}
	LAND_TYPE landtype=LAND_UNUSED;
	if (ass_flag)
	{
		if (constant.fcon[ipatch][2]==4 || constant.fcon[ipatch][2]==5)
			landtype=LAND_WATER;
		else if (constant.fcon[ipatch][2]==0)
		{
			if (constant.fcon[ipatch][3]<10 || (constant.fcon[ipatch][3]>15 && constant.fcon[ipatch][2]!=21))
			{
				landtype=LAND_SOIL;
				//frozen soil and snow judgement
				bool is_frosen_soil=true;
				for(int i=0;i<m_ensemble;i++)
				{
					if (fvar_en[i][ipatch][35]>273.16)
					{
						is_frosen_soil=false;
						break;
					}
				}
				if (is_frosen_soil)
				{
					landtype=LAND_FROSEN_SOIL;
					bool is_snow=true;
					for(int i=0;i<m_ensemble;i++)
					{
						if(snl_en[ipatch][i]<=0)
						{
							is_snow=false;
							break;
						}
					}
					if (is_snow)
						landtype=LAND_SNOW;
				}
			}
		}
	}
    Matrix<double> obs_total,obs,HXEn,XfEn,R,Q,XaEn;
    if (ass_flag)
    {
        switch(sensor_type)
        {
        case 1:
            obs_total.setDim(4,1);
            obs_total.set(1,1,0.1*m_observe->ssmi_data_19h[row][column]);
            obs_total.set(2,1,0.1*m_observe->ssmi_data_19v[row][column]);
            obs_total.set(3,1,0.1*m_observe->ssmi_data_37h[row][column]);
            obs_total.set(4,1,0.1*m_observe->ssmi_data_37v[row][column]);
            obs.setDim(2,1);
            obs.set(1,1,obs_total.get(1,1));
            obs.set(2,1,obs_total.get(2,1));
            HXEn.setDim(2,m_ensemble);
            R.setDim(2,2);
            R.genZeros();
            R.set(1,1,R_ssmi_19v);
            R.set(2,2,R_ssmi_19h);
            break;
        case 2:
            obs_total.setDim(8,1);
            obs_total.set(1,1,0.1*m_observe->amsr_data_6h[row][column]);
            obs_total.set(2,1,0.1*m_observe->amsr_data_6v[row][column]);
            obs_total.set(3,1,0.1*m_observe->amsr_data_10h[row][column]);
            obs_total.set(4,1,0.1*m_observe->amsr_data_10v[row][column]);
            obs_total.set(5,1,0.1*m_observe->amsr_data_18h[row][column]);
            obs_total.set(6,1,0.1*m_observe->amsr_data_18v[row][column]);
            obs_total.set(7,1,0.1*m_observe->amsr_data_37h[row][column]);
            obs_total.set(8,1,0.1*m_observe->amsr_data_37v[row][column]);
            obs.setDim(4,1);
            obs.set(1,1,obs_total.get(1,1));
            obs.set(2,1,obs_total.get(2,1));
            obs.set(3,1,obs_total.get(3,1));
            obs.set(4,1,obs_total.get(4,1));
            HXEn.setDim(4,m_ensemble);
            R.setDim(4,4);
            R.genZeros();
            R.set(1,1,R_amsr_6v);
            R.set(2,2,R_amsr_6h);
            R.set(3,3,R_amsr_10v);
            R.set(4,4,R_amsr_10h);
            break;
        case 3:
            obs_total.setDim(12,1);

            obs_total.set(1,1,0.1*m_observe->ssmi_data_19h[row][column]);
            obs_total.set(2,1,0.1*m_observe->ssmi_data_19v[row][column]);
            obs_total.set(3,1,0.1*m_observe->ssmi_data_37h[row][column]);
            obs_total.set(4,1,0.1*m_observe->ssmi_data_37v[row][column]);

            obs_total.set(5, 1,0.1*m_observe->amsr_data_6h[row][column]);
            obs_total.set(6, 1,0.1*m_observe->amsr_data_6v[row][column]);
            obs_total.set(7, 1,0.1*m_observe->amsr_data_10h[row][column]);
            obs_total.set(8, 1,0.1*m_observe->amsr_data_10v[row][column]);
            obs_total.set(9, 1,0.1*m_observe->amsr_data_18h[row][column]);
            obs_total.set(10,1,0.1*m_observe->amsr_data_18v[row][column]);
            obs_total.set(11,1,0.1*m_observe->amsr_data_37h[row][column]);
            obs_total.set(12,1,0.1*m_observe->amsr_data_37v[row][column]);

            obs.setDim(6,1);
            obs.set(1,1,obs_total.get(1,1));
            obs.set(2,1,obs_total.get(2,1));
            obs.set(3,1,obs_total.get(5,1));
            obs.set(4,1,obs_total.get(6,1));
            obs.set(5,1,obs_total.get(7,1));
            obs.set(6,1,obs_total.get(8,1));

            HXEn.setDim(6,m_ensemble);
            R.setDim(6,6);
            R.genZeros();
            R.set(1,1,R_ssmi_19v);
            R.set(2,2,R_ssmi_19h);
            R.set(3,3,R_amsr_6v);
            R.set(4,4,R_amsr_6h);
            R.set(5,5,R_amsr_10v);
            R.set(6,6,R_amsr_10h);
            break;
        default:
            break;
        }
        //************************************************************************
        // assimilation passive brightness temperature for different land type
        //***********************************************************************
        Matrix<double> SoilThick(10,1);
        Matrix<double> T2En(1,m_ensemble);
        Matrix<double> icEn(1,m_ensemble);

        double porosity,sand,clay,thick_layer1;
		porosity = constant.fcon[ipatch][15];
        sand     = sand_grid[row][column]*0.01;
        clay     = clay_grid[row][column]*0.01;
        EnKF enkf;
        QHModel qh;
        SoilParameter sp;
        switch(landtype)
        {
		case LAND_SOIL:
            // soil thick
            for (int i=0; i<10; i++)
                SoilThick.set(i+1,1,fvar_en[0][ipatch][20+i]);

            thick_layer1       = SoilThick.get(1,1);

            sp.sand=sand;
            sp.clay=clay;
            sp.porosity=porosity;

            XfEn.setDim(6,m_ensemble);
            Q.setDim(6,6);
            Q.set(1,1,0.008);
            Q.set(2,2,0.005);
            Q.set(3,3,0.003);
            Q.set(4,4,0.10);
            Q.set(5,5,0.08);
            Q.set(6,6,0.05);
            for (int ien=0; ien<m_ensemble; ien++)
            {
                XfEn.set(1,ien+1,fvar_en[ien][ipatch][35]);
                XfEn.set(2,ien+1,fvar_en[ien][ipatch][36]);
                XfEn.set(3,ien+1,fvar_en[ien][ipatch][37]);
                XfEn.set(4,ien+1,fvar_en[ien][ipatch][50]);
                XfEn.set(5,ien+1,fvar_en[ien][ipatch][51]);
                XfEn.set(6,ien+1,fvar_en[ien][ipatch][52]);
                sp.Me=fvar_en[ien][ipatch][50]/(denh2o*thick_layer1);
                if (sp.Me<0)
					sp.Me=0.1;
				else if (sp.Me>sp.porosity)
					sp.Me=sp.porosity;
                sp.Te=fvar_en[ien][ipatch][35];
                switch(sensor_type)
				{
				case 1:
                    sp.frequency=ssmi_freq_19;
                    sp.theta=ssmi_theta;
                    qh.soilparameter(sp);
                    qh.Q(0.3);
                    qh.h(0.5);
                    qh.run();
                    HXEn.set(1,ien+1,qh.tb().v);
                    HXEn.set(2,ien+1,qh.tb().h);
                    break;
                case 2:
                    sp.frequency=amsr_freq_6;
                    sp.theta=amsr_theta;
                    qh.soilparameter(sp);
                    qh.Q(0.3);
                    qh.h(0.5);
                    qh.run();
                    HXEn.set(1,ien+1,qh.tb().v);
                    HXEn.set(2,ien+1,qh.tb().h);
                    sp.frequency=amsr_freq_10;
                    qh.soilparameter(sp);
                    qh.Q(0.3);
                    qh.h(0.5);
                    qh.run();
                    HXEn.set(3,ien+1,qh.tb().v);
                    HXEn.set(4,ien+1,qh.tb().h);
                    break;
				case 3:
                    sp.frequency=ssmi_freq_19;
                    sp.theta=ssmi_theta;
                    qh.soilparameter(sp);
                    qh.Q(0.3);
                    qh.h(0.5);
                    qh.run();
                    HXEn.set(1,ien+1,qh.tb().v);
                    HXEn.set(2,ien+1,qh.tb().h);

                    sp.frequency=amsr_freq_6;
                    sp.theta=amsr_theta;
                    qh.soilparameter(sp);
                    qh.Q(0.3);
                    qh.h(0.5);
                    qh.run();
                    HXEn.set(3,ien+1,qh.tb().v);
                    HXEn.set(4,ien+1,qh.tb().h);
                    sp.frequency=amsr_freq_10;
                    qh.soilparameter(sp);
                    qh.Q(0.3);
                    qh.h(0.5);
                    qh.run();
                    HXEn.set(5,ien+1,qh.tb().v);
                    HXEn.set(6,ien+1,qh.tb().h);
                    break;
				default:
					break;
				}
            }


            enkf.update(XfEn,HXEn,obs,R,Q);
            XaEn=enkf.GetXaEn();
            // check soil water range
            for(int i=1; i<=3; i++)
            {
                for(int j=1; j<=m_ensemble; j++)
                {
                    double sm = XaEn.get(3+i,j)/(SoilThick.get(i,1)*denh2o);
                    if(sm<0.0)
                        XaEn.set(i+3,j,0);
                    else if(sm>porosity)
                        XaEn.set(i+3,j,porosity*denh2o*SoilThick.get(i,1));
                }
            }
            for(int i=0;i<m_ensemble;i++)
			{
				fvar_en[i][ipatch][35]=XaEn.get(1,i+1);
				fvar_en[i][ipatch][36]=XaEn.get(2,i+1);
				fvar_en[i][ipatch][37]=XaEn.get(3,i+1);
				fvar_en[i][ipatch][50]=XaEn.get(4,i+1);
				fvar_en[i][ipatch][51]=XaEn.get(5,i+1);
				fvar_en[i][ipatch][52]=XaEn.get(6,i+1);
			}
			fvar_mean[ipatch][35]=enkf.GetXa().get(1,1);
			fvar_mean[ipatch][36]=enkf.GetXa().get(2,1);
			fvar_mean[ipatch][37]=enkf.GetXa().get(3,1);
			fvar_mean[ipatch][50]=enkf.GetXa().get(4,1);
			fvar_mean[ipatch][51]=enkf.GetXa().get(5,1);
			fvar_mean[ipatch][52]=enkf.GetXa().get(6,1);
            break;
        case LAND_FROSEN_SOIL:
        	//todo
            break;
		case LAND_SNOW:
			break;
		case LAND_WATER:
			break;
		default:
			break;
        }
    }
}
