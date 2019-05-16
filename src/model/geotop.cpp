#include "geotop.h"
using namespace ldas;
using namespace ldas::geotop;
Geotop::Geotop()
{
    //ctor
    adt=new GeotopParameter;
}

Geotop::Geotop(int Dt, double JD0, int year0, double TH,double ST, double Dt_output_discharge,
               double Dt_output_pixel, double Dt_output_basin, short nsky, double channel_thres,
               short format_out, short point_sim,short recover)
{
    this->Dt=Dt;
    this->JD0;
    this->year0;
    this->TH=TH;
    this->ST=ST;
    this->Dt_output_discharge=Dt_output_discharge;
    this->Dt_output_pixel=Dt_output_pixel;
    this->Dt_output_basin=Dt_output_basin;
    this->nsky=nsky;
    this->channel_thres=channel_thres;
    this->format_out=format_out;
    this->point_sim=point_sim;
    this->recover=recover;
}

void Geotop::parameters(int Dt, double JD0, int year0, double TH,double ST, double Dt_output_discharge,
                        double Dt_output_pixel, double Dt_output_basin, short nsky, double channel_thres,
                        short format_out, short point_sim,short recover)
{
    this->Dt=Dt;
    this->JD0;
    this->year0;
    this->TH=TH;
    this->ST=ST;
    this->Dt_output_discharge=Dt_output_discharge;
    this->Dt_output_pixel=Dt_output_pixel;
    this->Dt_output_basin=Dt_output_basin;
    this->nsky=nsky;
    this->channel_thres=channel_thres;
    this->format_out=format_out;
    this->point_sim=point_sim;
    this->recover=recover;
}

Geotop::~Geotop()
{
    //dtor
    dealloc_all(adt->T, adt->S, adt->L, adt->W, adt->C, adt->P, adt->E, adt->N, adt->G, adt->M);
    delete adt;
}

Geotop::Geotop(const Geotop& other)
{
    //copy ctor
}

void Geotop::config(const char* fn)
{
    UV=(T_INIT *)malloc(sizeof(T_INIT));
    if(!UV) t_error("UV was not allocated");

    adt->open(fn);
}

void Geotop::init()
{

    double Dt_output;
//---------------------------------------------------base parameters
    //printf("ENTER THE INTEGRATION INTERVAL [s]: %f\n",V->co[1]);
    adt->P->Dt=Dt;	//THE INTEGRATION INTERVAL [s]
    //printf("ENTER THE Decimal julian day at the beginning of simulation (0.0 - 365.99): %f\n",V->co[2]);
    adt->P->JD0=JD0;
    //printf("ENTER THE YEAR at the beginning of simulation (0.0 - 365.99): %f\n",V->co[3]);
    adt->P->year0=year0;
    //printf("ENTER THE NUMBER OF DAYS OF SIMULATION: %f\n",V->co[4]);
    adt->I->TH=TH;
    adt->I->TH*=24; //TH in hour
    //printf("ENTER THE Standard time to which all the output data are referred (difference respect UMT, in hour): %f\n",V->co[5]);
    adt->P->ST=ST;

    //printf("ENTER THE NUMBER OF Dt AFTER WHICH THE DISCHARGE IS PRINTED: %f\n",V->co[6]);
    Dt_output=Dt_output_discharge;
    if(Dt_output>0 && Dt_output*3600<adt->P->Dt) Dt_output=adt->P->Dt/3600.0;
    adt->I->n_discharge=(long)(Dt_output*3600/(long)adt->P->Dt);
    adt->I->i_discharge=0;/*counter for the output of a pixel*/
    if(adt->I->n_discharge>0)
    {
        adt->P->state_discharge=1;
    }
    else
    {
        adt->P->state_discharge=0;
    }

    //printf("ENTER THE NUMBER OF Dt AFTER WHICH THE OUTPUT FOR A SPECIFIED PIXEL ARE PRINTED: %f\n",V->co[7]);
    Dt_output=Dt_output_pixel;		//(double)V->co[7];
    if(Dt_output>0 && Dt_output*3600<adt->P->Dt) Dt_output=adt->P->Dt/3600.0;
    adt->I->n_pixel=(long)(Dt_output*3600/(long)adt->P->Dt);
    adt->I->i_pixel=0;/*counter for the output of a pixel*/
    if(adt->I->n_pixel>0)
    {
        adt->P->state_pixel=1;
    }
    else
    {
        adt->P->state_pixel=0;
    }

    //printf("ENTER THE NUMBER OF Dt AFTER WHICH THE OUTPUT FOR THE BASIN ARE PRINTED: %f\n",V->co[8]);
    Dt_output=Dt_output_basin;		//(double)V->co[8];
    if(Dt_output>0 && Dt_output*3600<adt->P->Dt) Dt_output=adt->P->Dt/3600.0;
    adt->I->n_basin=(long)(Dt_output*3600/(long)adt->P->Dt);
    adt->I->i_basin=0;/*counter for the output of a pixel*/
    if(adt->I->n_basin>0)
    {
        adt->P->state_basin=1;
    }
    else
    {
        adt->P->state_basin=0;
    }

    //printf("Multiplying factor decreasing the dem resolution for the calculation of the sky view factor: %f\n",V->co[9]);
    adt->P->nsky=nsky;
    if(adt->P->nsky<1) t_error("Multiplying factor for sky view factor calculation must be greater than or equal to 1: %f");
    //printf("Thershold for the definition of channel network (in pixel number draining): %f\n",V->co[10]);
    adt->P->channel_thres=channel_thres;
    //printf("OUTPUT MAPS in fluidturtle format (=1), GRASS ASCII (=2), ESRI ASCII (=3): %f\n",V->co[11]);
    adt->P->format_out=format_out;
    //printf("=1 for one point simulation, =0 for distributed simulation. The parameter file is different in these cases: %f\n",V->co[12]);
    adt->P->point_sim=point_sim;
    //printf("=1 if you want to recover a simulation, 0 otherwise: %f\n",V->co[13]);
    adt->P->recover=recover;
//----------------------------------------------------------------
}

Geotop& Geotop::operator=(const Geotop& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}

void Geotop::step()
{
    updates_times(adt->I, adt->P);
    //驱动数据处理
    meteo_distr(adt->M, adt->E, adt->W, adt->T, adt->N, adt->I->time, adt->P);
    if(adt->P->en_balance==1) energy_balance(adt->I, adt->P, adt->L, adt->T, adt->S, adt->M, adt->W, adt->E, adt->N, adt->G);
    if(adt->P->wat_balance==1) water_balance(adt->alldata());
}

void Geotop::run()
{
    do
    {
        step();
        output();
        adt->I->time+=adt->P->Dt;//Increase TIME
    }
    while(adt->I->time<(adt->I->TH*3600.0)); //end of time-cycle
}

void Geotop::output()
{
    adt->save();
}

GeotopParameter::GeotopParameter()
{
    T=new TOPO;
    S=new SOIL;
    P=new PAR;
    C=new CHANNEL;
    E=new ENERGY;
    N=new SNOW;
    G=new GLACIER;
    M=new METEO;
    I=new TIMES;
    W=new WATER;
    L=new LAND;
    adt=new ALLDATA;
}
GeotopParameter::~GeotopParameter()
{
    /* prevent second free
    delete T;
    delete S;
    delete P;
    delete C;
    delete E;
    delete N;
    delete G;
    delete M;
    delete W;
    delete L;
    */
    delete I;
    delete adt;
}
void GeotopParameter::open(const char* fn)
{
    char* path[2];
    path[0]="not_used";
    strcpy(path[1],fn);
    //path[1]=fn;
    get_all_input(2, path, T, S, L, M, W, C, P, E, N, G, I);
}
void GeotopParameter::save()
{
    write_output(I, W, C, P,T,L, S, E, N, G, M);
}
ALLDATA* GeotopParameter::alldata()
{
    adt->T=T;
    adt->S=S;
    adt->P=P;
    adt->C=C;
    adt->E=E;
    adt->N=N;
    adt->G=G;
    adt->M=M;
    adt->I=I;
    adt->W=W;
    adt->L=L;
    return adt;
}
GeotopParameter::GeotopParameter(const GeotopParameter& other)
{
    //copy ctor
}

GeotopParameter& GeotopParameter::operator=(const GeotopParameter& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator

    E->Rn_mean=new_doublematrix(Nr,Nc);

    if(rhs.P->distr_stat==1)
    {
        E->Rn_min=new_doublematrix(Nr,Nc);
        E->Rn_max=new_doublematrix(Nr,Nc);
        E->LW_max=new_doublematrix(Nr,Nc);
        E->LW_min=new_doublematrix(Nr,Nc);
        E->SW_max=new_doublematrix(Nr,Nc);
    }

    //E->LWin_mean=new_doublematrix(Nr,Nc);
    //E->LW_mean=new_doublematrix(Nr,Nc);
    //E->SW_mean=new_doublematrix(Nr,Nc);
    //if(par->distr_stat==1)egy->SW_max=new_doublematrix(Nr,Nc);

    E->ET_mean=new_doublematrix(Nr,Nc);
    if(rhs.P->distr_stat==1)E->ET_max=new_doublematrix(Nr,Nc);
    if(rhs.P->distr_stat==1)E->ET_min=new_doublematrix(Nr,Nc);

    if(rhs.P->output_H>0)E->H_mean=new_doublematrix(Nr,Nc);
    if(rhs.P->distr_stat==1 && rhs.P->output_H>0)E->H_max=new_doublematrix(Nr,Nc);
    if(rhs.P->distr_stat==1 && rhs.P->output_H>0)E->H_min=new_doublematrix(Nr,Nc);

    if (rhs.P->output_G>0)E->SEB_mean=new_doublematrix(Nr,Nc);
    if(rhs.P->distr_stat==1 && rhs.P->output_G>0)E->G_max=new_doublematrix(Nr,Nc);
    if(rhs.P->distr_stat==1 && rhs.P->output_G>0)E->G_min=new_doublematrix(Nr,Nc);
    if(rhs.P->output_G>0)E->G_snowsoil=new_doublematrix(Nr,Nc);

    if(rhs.P->output_Ts>0)E->Ts_mean=new_doublematrix(Nr,Nc);
    if(rhs.P->output_Ts>0&&rhs.P->distr_stat==1)E->Ts_max=new_doublematrix(Nr,Nc);
    if(rhs.P->output_Ts>0&&rhs.P->distr_stat==1)E->Ts_min=new_doublematrix(Nr,Nc);

    if(rhs.P->output_Rswdown>0)E->Rswdown_mean=new_doublematrix(Nr,Nc);
    if(rhs.P->output_Rswdown>0&&rhs.P->distr_stat==1)E->Rswdown_max=new_doublematrix(Nr,Nc);
    if(rhs.P->output_Rswdown>0)E->Rswbeam=new_doublematrix(Nr,Nc);


    E->Rn_mean=new_doublematrix(Nr,Nc);
    //E->Rn_max=new_doublematrix(Nr,Nc);
    //E->Rn_min=new_doublematrix(Nr,Nc);
    E->LW_in=new_doublematrix(Nr,Nc);
    E->LW=new_doublematrix(Nr,Nc);
    //E->LW_max=new_doublematrix(Nr,Nc);
    //E->LW_min=new_doublematrix(Nr,Nc);
    E->SW=new_doublematrix(Nr,Nc);

    //E->ET_mean=new_doublematrix(Nr,Nc);
    //E->ET_max=new_doublematrix(Nr,Nc);
    //E->ET_min=new_doublematrix(Nr,Nc);
    //E->H_mean=new_doublematrix(Nr,Nc);
    //E->H_max=new_doublematrix(Nr,Nc);
    //E->H_min=new_doublematrix(Nr,Nc);
    //E->SEB_mean=new_doublematrix(Nr,Nc);
    //E->G_max=new_doublematrix(Nr,Nc);
    //E->G_min=new_doublematrix(Nr,Nc);
    //E->G_snowsoil=new_doublematrix(Nr,Nc);
    //E->Ts_mean=new_doublematrix(Nr,Nc);
    //E->Ts_max=new_doublematrix(Nr,Nc);
    //E->Ts_min=new_doublematrix(Nr,Nc);
    //E->Rswdown_mean=new_doublematrix(Nr,Nc);
    //E->Rswdown_max=new_doublematrix(Nr,Nc);
    if(rhs.P->output_meteo>0 )E->Ta_mean=new_doublematrix(Nr,Nc);
    if(rhs.P->output_meteo>0&&rhs.P->distr_stat==1)E->Ta_max=new_doublematrix(Nr,Nc);
    if(rhs.P->output_meteo>0&&rhs.P->distr_stat==1)E->Ta_min=new_doublematrix(Nr,Nc);
    //E->Rswbeam=new_doublematrix(Nr,Nc);

    if(rhs.P->JD_plots->co[1]>=0)E->Hgplot=new_doublematrix(Nr,Nc);
    if(rhs.P->JD_plots->co[1]>=0)E->LEgplot=new_doublematrix(Nr,Nc);
    if(rhs.P->JD_plots->co[1]>=0)E->Hvplot=new_doublematrix(Nr,Nc);
    if(rhs.P->JD_plots->co[1]>=0)E->LEvplot=new_doublematrix(Nr,Nc);

    if(rhs.P->JD_plots->co[1]>=0)E->SWinplot=new_doublematrix(Nr,Nc);
    if(rhs.P->JD_plots->co[1]>=0)E->SWgplot=new_doublematrix(Nr,Nc);
    if(rhs.P->JD_plots->co[1]>=0)E->SWvplot=new_doublematrix(Nr,Nc);

    if(rhs.P->JD_plots->co[1]>=0)E->LWinplot=new_doublematrix(Nr,Nc);
    if(rhs.P->JD_plots->co[1]>=0)E->LWgplot=new_doublematrix(Nr,Nc);
    if(rhs.P->JD_plots->co[1]>=0)E->LWvplot=new_doublematrix(Nr,Nc);

    if(rhs.P->JD_plots->co[1]>=0)E->Tgplot=new_doublematrix(Nr,Nc);
    if(rhs.P->JD_plots->co[1]>=0)E->Tvplot=new_doublematrix(Nr,Nc);
    if(rhs.P->JD_plots->co[1]>=0)E->Tsplot=new_doublematrix(Nr,Nc);

    //E->SWin=new_doublematrix(Nr,Nc);
    //E->LWin=new_doublematrix(Nr,Nc);

    int r,c;
    for(r=1; r<=Nr; r++)
    {
        for(c=1; c<=Nc; c++)
        {
            E->Rn_mean->co[r][c]=rhs.E->Rn_mean->co[r][c];

            if(rhs.P->distr_stat==1 )
            {
                E->Rn_max->co[r][c]=rhs.E->Rn_max->co[r][c];
                E->Rn_min->co[r][c]=rhs.E->Rn_min->co[r][c];
                E->LW_max->co[r][c]=rhs.E->LW_max->co[r][c];
                E->LW_min->co[r][c]=rhs.E->LW_min->co[r][c];
                E->SW_max->co[r][c]=rhs.E->SW_max->co[r][c];
            }
            E->LW_in->co[r][c]=rhs.E->LW_in->co[r][c];
            E->LW->co[r][c]=rhs.E->LW->co[r][c];
            E->SW->co[r][c]=rhs.E->SW->co[r][c];

            if(rhs.P->distr_stat==1)E->ET_mean->co[r][c]=rhs.E->ET_mean->co[r][c];
            if(rhs.P->distr_stat==1)E->ET_max->co[r][c]=rhs.E->ET_max->co[r][c];
            if(rhs.P->distr_stat==1)E->ET_min->co[r][c]=rhs.E->ET_min->co[r][c];
            if(rhs.P->output_H>0)E->H_mean->co[r][c]=rhs.E->H_mean->co[r][c];
            if(rhs.P->distr_stat==1 && rhs.P->output_H>0)E->H_max->co[r][c]=rhs.E->H_max->co[r][c];
            if(rhs.P->distr_stat==1 && rhs.P->output_H>0)E->H_min->co[r][c]=rhs.E->H_min->co[r][c];
            if(rhs.P->output_G>0)E->SEB_mean->co[r][c]=rhs.E->SEB_mean->co[r][c];
            if(rhs.P->distr_stat==1 && rhs.P->output_G>0)E->G_max->co[r][c]=rhs.E->G_max->co[r][c];
            if(rhs.P->distr_stat==1 && rhs.P->output_G>0)E->G_min->co[r][c]=rhs.E->G_min->co[r][c];
            if(rhs.P->output_G>0)E->G_snowsoil->co[r][c]=rhs.E->G_snowsoil->co[r][c];
            if(rhs.P->output_Ts>0)E->Ts_mean->co[r][c]=rhs.E->Ts_mean->co[r][c];
            if(rhs.P->output_Ts>0&&rhs.P->distr_stat==1)E->Ts_max->co[r][c]=rhs.E->Ts_max->co[r][c];
            if(rhs.P->output_Ts>0&&rhs.P->distr_stat==1)E->Ts_min->co[r][c]=rhs.E->Ts_min->co[r][c];
            if(rhs.P->output_Rswdown>0)E->Rswdown_mean->co[r][c]=rhs.E->Rswdown_mean->co[r][c];
            if(rhs.P->output_Rswdown>0&&rhs.P->distr_stat==1)E->Rswdown_max->co[r][c]=rhs.E->Rswdown_max->co[r][c];
            if(rhs.P->output_meteo>0)E->Ta_mean->co[r][c]=rhs.E->Ta_mean->co[r][c];
            if(rhs.P->output_meteo>0&&rhs.P->distr_stat==1)E->Ta_max->co[r][c]=rhs.E->Ta_max->co[r][c];
            if(rhs.P->output_meteo>0&&rhs.P->distr_stat==1)E->Ta_min->co[r][c]=rhs.E->Ta_min->co[r][c];
            if(rhs.P->output_Rswdown>0)E->Rswbeam->co[r][c]=rhs.E->Rswbeam->co[r][c];

            if(rhs.P->JD_plots->co[1]>=0)E->Hgplot->co[r][c]=rhs.E->Hgplot->co[r][c];
            if(rhs.P->JD_plots->co[1]>=0)E->LEgplot->co[r][c]=rhs.E->LEgplot->co[r][c];
            if(rhs.P->JD_plots->co[1]>=0)E->Hvplot->co[r][c]=rhs.E->Hvplot->co[r][c];
            if(rhs.P->JD_plots->co[1]>=0)E->LEvplot->co[r][c]=rhs.E->LEvplot->co[r][c];

            if(rhs.P->JD_plots->co[1]>=0)E->SWinplot->co[r][c]=rhs.E->SWinplot->co[r][c];
            if(rhs.P->JD_plots->co[1]>=0)E->SWgplot->co[r][c]=rhs.E->SWgplot->co[r][c];
            if(rhs.P->JD_plots->co[1]>=0)E->SWvplot->co[r][c]=rhs.E->SWvplot->co[r][c];

            if(rhs.P->JD_plots->co[1]>=0)E->LWinplot->co[r][c]=rhs.E->LWinplot->co[r][c];
            if(rhs.P->JD_plots->co[1]>=0)E->LWgplot->co[r][c]=rhs.E->LWgplot->co[r][c];
            if(rhs.P->JD_plots->co[1]>=0)E->LWvplot->co[r][c]=rhs.E->LWvplot->co[r][c];

            if(rhs.P->JD_plots->co[1]>=0)E->Tgplot->co[r][c]=rhs.E->Tgplot->co[r][c];
            if(rhs.P->JD_plots->co[1]>=0)E->Tvplot->co[r][c]=rhs.E->Tvplot->co[r][c];
            if(rhs.P->JD_plots->co[1]>=0)E->Tsplot->co[r][c]=rhs.E->Tsplot->co[r][c];

            //E->SWin->co[r][c]=rhs.E->SWin->co[r][c];
            //E->LWin->co[r][c]=rhs.E->LWin->co[r][c];

        }
    }

    E->nDt_shadow=new_longmatrix(Nr,Nc);
    E->nDt_sun=new_longmatrix(Nr,Nc);
    initialize_longmatrix(E->nDt_shadow,0);
    initialize_longmatrix(E->nDt_sun,0);


    E->hsun=rhs.E->hsun;
    E->dsun=rhs.E->dsun;
    E->Dlay= new_doublevector( Nl + rhs.P->snowlayer_max + rhs.P->glaclayer_max );
    E->wliq= new_doublevector( Nl + rhs.P->snowlayer_max + rhs.P->glaclayer_max );
    E->wice= new_doublevector( Nl + rhs.P->snowlayer_max + rhs.P->glaclayer_max );
    E->Temp= new_doublevector( Nl + rhs.P->snowlayer_max + rhs.P->glaclayer_max );
    E->deltaw= new_doublevector( Nl + rhs.P->snowlayer_max + rhs.P->glaclayer_max );
    E->SWlayer= new_doublevector( rhs.P->snowlayer_max+1 );

    E->soil_transp_layer= new_doublevector(rhs.L->root_fraction->nch);
    initialize_doublevector(E->soil_transp_layer, 0.);

    E->dFenergy= new_doublevector( Nl + rhs.P->snowlayer_max + rhs.P->glaclayer_max );
    E->Kth0=new_doublevector( Nl + rhs.P->snowlayer_max + rhs.P->glaclayer_max );
    E->Kth1=new_doublevector( Nl + rhs.P->snowlayer_max + rhs.P->glaclayer_max );
    E->Fenergy=new_doublevector( Nl + rhs.P->snowlayer_max + rhs.P->glaclayer_max );
    E->Newton_dir=new_doublevector( Nl + rhs.P->snowlayer_max + rhs.P->glaclayer_max );
    E->T0=new_doublevector( Nl + rhs.P->snowlayer_max + rhs.P->glaclayer_max );
    E->T1=new_doublevector( Nl + rhs.P->snowlayer_max + rhs.P->glaclayer_max );
    E->Tstar=new_doublevector(Nl);
    E->THETA=new_doublevector(Nl);

    int l = 0;
    double z=0;
    do
    {
        l++;
        z += rhs.S->pa->co[1][jdz][l];
    }
    while(l<Nl && z < z_evap);

    E->soil_evap_layer_bare = new_doublevector(l);
    E->soil_evap_layer_veg = new_doublevector(l);
    initialize_doublevector(E->soil_evap_layer_bare, 0.);
    initialize_doublevector(E->soil_evap_layer_bare, 0.);







////////////////////////////////////////

    S->type=new_longmatrix(rhs.S->type->nrh,rhs.S->type->nch);
    for(int r=1; r<=rhs.S->type->nrh; r++)
    {
        for(int c=1; c<=rhs.S->type->nch; c++)
        {
            S->type->co[r][c]=rhs.S->type->co[r][c];
        }
    }
    //S->type=copylong_doublematrix(rhs.S->type);
    S->pa=new_doubletensor(rhs.S->pa->ndh, rhs.S->pa->nrh, rhs.S->pa->nch);/** contains the soil parameters from the __soil.txt, pa[SoilType][Property][Layer] */
    initialize_doubletensor(S->pa, 0.0);
    //DOUBLETENSOR *P;/**Psi,P[layer][Row][Cols]*/
    S->P=new_doubletensor0(rhs.S->P->ndh, rhs.S->P->nrh, rhs.S->P->nch);
    initialize_doubletensor(S->P,0.0);
    S->Ptot=new_doubletensor(rhs.S->Ptot->ndh, rhs.S->Ptot->nrh, rhs.S->Ptot->nch);
    initialize_doubletensor(S->Ptot,0.0);
    S->T=new_doubletensor(rhs.S->T->ndh, rhs.S->T->nrh, rhs.S->T->nch);
    initialize_doubletensor(S->T,0.0);	/**Temperature*/
    S->thice=new_doubletensor(rhs.S->thice->ndh, rhs.S->thice->nrh, rhs.S->thice->nch);
    initialize_doubletensor(S->thice,0.0);/**Volume metric ice fraction*/
    S->th=new_doubletensor(rhs.S->th->ndh, rhs.S->th->nrh, rhs.S->th->nch);
    initialize_doubletensor(S->th,0.0);	/**Volume metric water fraction*/
    S->Jinf=NULL;
    S->Tv=new_doublematrix(rhs.S->Tv->nrh,rhs.S->Tv->nch);
    initialize_doublematrix(S->Tv,0.0);/**vegetation temperature*/
    S->bc=NULL;
    S->ET=new_doubletensor(rhs.S->ET->ndh, rhs.S->ET->nrh, rhs.S->ET->nch);
    initialize_doubletensor(S->ET,0.0);
    if(rhs.P->state_pixel==1)
    {
        S->T_av=new_doublematrix(Nl, rhs.P->chkpt->nrh);
        S->th_av=new_doublematrix(Nl, rhs.P->chkpt->nrh);
        S->thice_av=new_doublematrix(Nl, rhs.P->chkpt->nrh);
    }
    for(int i=1; i<=Nl; i++)
    {
        for(int j=1; j<=rhs.P->chkpt->nrh; j++)
        {
            S->T_av->co[i][j]=rhs.S->T_av->co[i][j];
            S->th_av->co[i][j]=rhs.S->th_av->co[i][j];
            S->thice_av->co[i][j]=rhs.S->thice_av->co[i][j];
        }
    }
    for(int i=rhs.S->pa->ndl; i<=rhs.S->pa->ndh; i++)
    {
        for(int j=rhs.S->pa->nrl; j<=rhs.S->pa->nrh; j++)
        {
            for(int k=rhs.S->pa->ncl; k<=rhs.S->pa->nch; k++)
            {
                S->pa->co[i][j][k] = rhs.S->pa->co[i][j][k];
            }
        }

    }

    for(int i=rhs.S->P->ndl; i<=rhs.S->P->ndh; i++)
    {
        for(int j=rhs.S->P->nrl; j<=rhs.S->P->nrh; j++)
        {
            for(int k=rhs.S->P->ncl; k<=rhs.S->P->nch; k++)
            {
                S->P->co[i][j][k] = rhs.S->P->co[i][j][k];
            }
        }

    }

    for(int i=rhs.S->Ptot->ndl; i<=rhs.S->Ptot->ndh; i++)
    {
        for(int j=rhs.S->Ptot->nrl; j<=rhs.S->Ptot->nrh; j++)
        {
            for(int k=rhs.S->Ptot->ncl; k<=rhs.S->Ptot->nch; k++)
            {
                S->Ptot->co[i][j][k] = rhs.S->Ptot->co[i][j][k];

            }
        }

    }

    for(int i=rhs.S->T->ndl; i<=rhs.S->T->ndh; i++)
    {
        for(int j=rhs.S->T->nrl; j<=rhs.S->T->nrh; j++)
        {
            for(int k=rhs.S->T->ncl; k<=rhs.S->T->nch; k++)
            {

                S->T->co[i][j][k] = rhs.S->T->co[i][j][k];

            }
        }

    }

    for(int i=rhs.S->thice->ndl; i<=rhs.S->thice->ndh; i++)
    {
        for(int j=rhs.S->thice->nrl; j<=rhs.S->thice->nrh; j++)
        {
            for(int k=rhs.S->thice->ncl; k<=rhs.S->thice->nch; k++)
            {

                S->thice->co[i][j][k] = rhs.S->thice->co[i][j][k];

            }
        }

    }

    for(int i=rhs.S->th->ndl; i<=rhs.S->th->ndh; i++)
    {
        for(int j=rhs.S->th->nrl; j<=rhs.S->th->nrh; j++)
        {
            for(int k=rhs.S->th->ncl; k<=rhs.S->th->nch; k++)
            {

                S->th->co[i][j][k] = rhs.S->th->co[i][j][k];

            }
        }

    }

    for(int i=rhs.S->ET->ndl; i<=rhs.S->ET->ndh; i++)
    {
        for(int j=rhs.S->ET->nrl; j<=rhs.S->ET->nrh; j++)
        {
            for(int k=rhs.S->ET->ncl; k<=rhs.S->ET->nch; k++)
            {

                S->ET->co[i][j][k] = rhs.S->ET->co[i][j][k];
            }
        }

    }


////////////////////////////////////////
    //std::cout<<rhs.T->Z0->nrh<<"  "<<rhs.T->Z0->nch<<std::endl;
    T->Z0=new_doublematrix(rhs.T->Z0->nrh,rhs.T->Z0->nch);         /*elevetions of each pixel (DEM)*/
    for(int i=1; i<=rhs.T->Z0->nrh; i++)
    {
        for(int j=1; j<=rhs.T->Z0->nch; j++)
        {
            //std::cout<<rhs.T->Z0->co[i][j]<<std::endl;
            T->Z0->co[i][j]=rhs.T->Z0->co[i][j];
        }
    }


    if(rhs.T->Z1!=0)
    {
        T->Z1=new_doublematrix(rhs.T->Z1->nrh,rhs.T->Z1->nch);
        for(int i=1; i<=rhs.T->Z1->nrh; i++)
        {
            for(int j=1; j<=rhs.T->Z1->nch; j++)
            {
                T->Z1->co[i][j]=rhs.T->Z1->co[i][j];
            }
        }
    }

    T->Z0dp=new_doublematrix(rhs.T->Z0dp->nrh,rhs.T->Z0dp->nch);		  /*DEM depitted*/
    for(int i=1; i<=rhs.T->Z0dp->nrh; i++)
    {
        for(int j=1; j<=rhs.T->Z0dp->nch; j++)
        {
            T->Z0dp->co[i][j]=rhs.T->Z0dp->co[i][j];
        }
    }

    if(rhs.T->Z0ext!=0)
    {
        T->Z0ext=new_doublematrix(rhs.T->Z0ext->nrh,rhs.T->Z0ext->nch);     /*DEM extended (to avoid curvature problems)*/
        for(int i=1; i<=rhs.T->Z0ext->nrh; i++)
        {
            for(int j=1; j<=rhs.T->Z0ext->nch; j++)
            {
                T->Z0ext->co[i][j]=rhs.T->Z0ext->co[i][j];
            }
        }
    }


    T->Z=new_doubletensor(rhs.T->Z->ndh,rhs.T->Z->nrh,rhs.T->Z->nch);
    for(int i=1; i<=rhs.T->Z->ndh; i++)
    {
        for(int j=1; j<=rhs.T->Z->nrh; j++)
        {
            for(int k=1; k<=rhs.T->Z->nch; k++)
            {
                T->Z->co[i][j][k]=rhs.T->Z->co[i][j][k];
            }
        }

    }


    T->sky=new_doublematrix(rhs.T->sky->nrh,rhs.T->sky->nch);       /*view factor (of the sky) for each pixel*/
    for(int i=1; i<=rhs.T->sky->nrh; i++)
    {
        for(int j=1; j<=rhs.T->sky->nch; j++)
        {
            T->sky->co[i][j]=rhs.T->sky->co[i][j];
        }
    }

    T->pixel_type=new_shortmatrix(rhs.T->pixel_type->nrh,rhs.T->pixel_type->nch);
    for(int i=1; i<=rhs.T->pixel_type->nrh; i++)
    {
        for(int j=1; j<=rhs.T->pixel_type->nch; j++)
        {
            T->pixel_type->co[i][j]=rhs.T->pixel_type->co[i][j];
        }
    }


    T->DD=new_shortmatrix(rhs.T->DD->nrh,rhs.T->DD->nch);
    for(int i=1; i<=rhs.T->DD->nrh; i++)
    {
        for(int j=1; j<=rhs.T->DD->nch; j++)
        {
            T->DD->co[i][j]=rhs.T->DD->co[i][j];
        }
    }
    /*Drainage Directions for each pixel; ex matr_ev->slopes*/

    if(rhs.T->DDup!=0)
    {
        T->DDup=new_longmatrix(rhs.T->DDup->nrh,rhs.T->DDup->nch);
        for(int i=1; i<=rhs.T->DDup->nrh; i++)
        {
            for(int j=1; j<=rhs.T->DDup->nch; j++)
            {
                T->DDup->co[i][j]=rhs.T->DDup->co[i][j];
            }
        }
    }

    if(rhs.T->DDdown!=0)
    {
        T->DDdown=new_longvector(rhs.T->DDdown->nh);
        for(int i=1; i<=rhs.T->DDdown->nh; i++)
            T->DDdown->co[i]=rhs.T->DDdown->co[i];
    }

    if(rhs.T->i_DD!=0)
    {
        T->i_DD=new_doublematrix(rhs.T->i_DD->nrh,rhs.T->i_DD->nch);       /*view factor (of the sky) for each pixel*/
        for(int i=1; i<=rhs.T->i_DD->nrh; i++)
        {
            for(int j=1; j<=rhs.T->i_DD->nch; j++)
            {
                T->i_DD->co[i][j]=rhs.T->i_DD->co[i][j];
            }
        }
    }

    if(rhs.T->top_index!=0)
    {
        T->top_index=new_doublematrix(rhs.T->top_index->nrh,rhs.T->top_index->nch);       /*view factor (of the sky) for each pixel*/
        for(int i=1; i<=rhs.T->top_index->nrh; i++)
        {
            for(int j=1; j<=rhs.T->top_index->nch; j++)
            {
                T->top_index->co[i][j]=rhs.T->top_index->co[i][j];
            }
        }
    }


    if(rhs.T->area!=0)
    {
        T->area=new_doublematrix(rhs.T->area->nrh,rhs.T->area->nch);       /*view factor (of the sky) for each pixel*/
        for(int i=1; i<=rhs.T->area->nrh; i++)
        {
            for(int j=1; j<=rhs.T->area->nch; j++)
            {
                T->area->co[i][j]=rhs.T->area->co[i][j];
            }
        }
    }

    ///DOUBLEMATRIX *aspect;     /*aspect; ex: matr_ev->azimuth*/
    T->aspect=new_doublematrix(rhs.T->aspect->nrh,rhs.T->aspect->nch);       /*view factor (of the sky) for each pixel*/
    for(int i=1; i<=rhs.T->aspect->nrh; i++)
    {
        for(int j=1; j<=rhs.T->aspect->nch; j++)
        {
            T->aspect->co[i][j]=rhs.T->aspect->co[i][j];
        }
    }


///    DOUBLEMATRIX *slopes;     /*slope of the pixels; ex: matr_ev->slopes*/
    T->slopes=new_doublematrix(rhs.T->slopes->nrh,rhs.T->slopes->nch);       /*view factor (of the sky) for each pixel*/
    for(int i=1; i<=rhs.T->slopes->nrh; i++)
    {
        for(int j=1; j<=rhs.T->slopes->nch; j++)
        {
            T->slopes->co[i][j]=rhs.T->slopes->co[i][j];
        }
    }


    T->Zm=new_doublematrix(rhs.T->Zm->nrh,rhs.T->Zm->nch);       /*view factor (of the sky) for each pixel*/
    for(int i=1; i<=rhs.T->Zm->nrh; i++)
    {
        for(int j=1; j<=rhs.T->Zm->nch; j++)
        {
            T->Zm->co[i][j]=rhs.T->Zm->co[i][j];
        }
    }


    T->curv_m=new_doublematrix(rhs.T->curv_m->nrh,rhs.T->curv_m->nch);       /*view factor (of the sky) for each pixel*/
    for(int i=1; i<=rhs.T->curv_m->nrh; i++)
    {
        for(int j=1; j<=rhs.T->curv_m->nch; j++)
        {
            T->curv_m->co[i][j]=rhs.T->curv_m->co[i][j];
        }
    }



    T->slope_m=new_doublematrix(rhs.T->slope_m->nrh,rhs.T->slope_m->nch);
    for(int i=1; i<=rhs.T->slope_m->nrh; i++)
    {
        for(int j=1; j<=rhs.T->slope_m->nch; j++)
        {
            T->slope_m->co[i][j]=rhs.T->slope_m->co[i][j];
        }
    }



    T->slopeaz_m=new_doublematrix(rhs.T->slopeaz_m->nrh,rhs.T->slopeaz_m->nch);
    for(int i=1; i<=rhs.T->slopeaz_m->nrh; i++)
    {
        for(int j=1; j<=rhs.T->slopeaz_m->nch; j++)
        {
            T->slopeaz_m->co[i][j]=rhs.T->slopeaz_m->co[i][j];
        }
    }



    int i;
    if(rhs.T->horizon_height!=0)
    {
        T->horizon_height=(double ****)malloc(T->Z0->nrh*sizeof(double***));
        for(int r=1; r<=T->Z0->nrh; r++)
        {
            T->horizon_height[r-1]=(double ***)malloc(T->Z0->nch*sizeof(double**));
            for(int c=1; c<=rhs.T->Z0->nch; c++)
            {

                if(rhs.T->horizon_height[r-1][c-1]!=0)
                {
                    T->horizon_height[r-1][c-1]=alloc2(4,2);
                    for(int j=1; j<=4; j++)
                    {
                        T->horizon_height[r-1][c-1][j-1][0]=rhs.T->horizon_height[r-1][c-1][j-1][0];
                        T->horizon_height[r-1][c-1][j-1][1]=rhs.T->horizon_height[r-1][c-1][j-1][1];
                    }

                }
            }
        }
    }


    T->i_cont=(long ***)malloc((Nl+1)*sizeof(long**));
    for(int l=0; l<=Nl; l++)
    {
        T->i_cont[l]=(long **)malloc((Nr+1)*sizeof(long*));
        for(int r=1; r<=Nr; r++)
        {
            T->i_cont[l][r]=(long *)malloc((Nc+1)*sizeof(long));
        }
    }

    long cont=0;
    //std::cout<<rhs.T->lrc_cont<<std::endl;
    if(rhs.T->lrc_cont!=0)
    {
        //T->lrc_cont=new_longmatrix((Nl+1)*rhs.P->total_pixel + rhs.C->r->nh , 3);

        T->lrc_cont=new_longmatrix(rhs.T->lrc_cont->nrh, 3);
        //std::cout<<rhs.T->lrc_cont->nrh<<std::endl;
        initialize_longmatrix(T->lrc_cont, 0);
        std::cout<<rhs.T->lrc_cont->co[1][1]<<std::endl;
        for(int i=1; i<=rhs.T->lrc_cont->nrh; i++)
        {
            T->lrc_cont->co[i][1]=rhs.T->lrc_cont->co[i][1];
            T->lrc_cont->co[i][2]=rhs.T->lrc_cont->co[i][2];
            T->lrc_cont->co[i][3]=rhs.T->lrc_cont->co[i][3];
        }
    }

    for(int r=1; r<=Nr; r++)
    {
        for(int c=1; c<=Nc; c++)
        {
            for(l=0; l<=Nl; l++)
            {
                T->i_cont[l][r][c]=rhs.T->i_cont[l][r][c];
            }
        }
    }




    T->j_cont=(long **)malloc((Nr+1)*sizeof(long*));
    for(int r=1; r<=Nr; r++)
    {
        T->j_cont[r]=(long *)malloc((Nc+1)*sizeof(long));
    }

    //T->rc_cont=new_longmatrix(rhs.P->total_pixel,2);
    if(rhs.T->rc_cont!=0)
    {
        T->rc_cont=new_longmatrix(rhs.T->rc_cont->nrh,2);
        initialize_longmatrix(T->rc_cont, 0);
        for(int i=1; i<=rhs.T->rc_cont->nrh; i++)
        {
            T->rc_cont->co[i][1]=rhs.T->rc_cont->co[i][1];
            T->rc_cont->co[i][2]=rhs.T->rc_cont->co[i][1];
        }
    }

    for(int r=1; r<=Nr; r++)
    {
        for(int c=1; c<=Nc; c++)
        {
            T->j_cont[r][c]=rhs.T->j_cont[r][c];
        }
    }


    ///DOUBLETENSOR *Z;


    T->Lp=new_longvector(rhs.T->Lp->nh);
    for(int i=1; i<=rhs.T->Lp->nh; i++)
        T->Lp->co[i]=rhs.T->Lp->co[i];

    T->Li=new_longvector(rhs.T->Li->nh);
    for(int i=1; i<=rhs.T->Li->nh; i++)
        T->Li->co[i]=rhs.T->Li->co[i];


    if(rhs.T->Up != 0)
    {
        T->Up=new_longvector(rhs.T->Up->nh);
        for(int i=1; i<=rhs.T->Up->nh; i++)
            T->Up->co[i]=rhs.T->Up->co[i];
    }

    if(rhs.T->Ui!=0)
    {
        T->Ui=new_longvector(rhs.T->Ui->nh);
        for(int i=1; i<=rhs.T->Ui->nh; i++)
            T->Ui->co[i]=rhs.T->Ui->co[i];
    }

////////////////////////////////////////Land

    L->LC=new_doublematrix(rhs.L->LC->nrh,rhs.L->LC->nch);
    for(int i=1; i<=rhs.L->LC->nrh; i++)
    {
        for(int j=1; j<=rhs.L->LC->nch; j++)
        {
            L->LC->co[i][j]=rhs.L->LC->co[i][j];
        }
    }


    //std::cout<<rhs.L->albedo<<std::endl;
    if(rhs.L->albedo!=0)
    {
        L->albedo=new_doublematrix(rhs.L->albedo->nrh,rhs.L->albedo->nch);
        for(int i=1; i<=rhs.L->albedo->nrh; i++)
        {
            for(int j=1; j<=rhs.L->albedo->nch; j++)
            {
                L->albedo->co[i][j]=rhs.L->albedo->co[i][j];
            }
        }
    }


    L->shadow=new_shortmatrix(rhs.L->shadow->nrh,rhs.L->shadow->nch);
    for(int i=1; i<=rhs.L->shadow->nrh; i++)
    {
        for(int j=1; j<=rhs.L->shadow->nch; j++)
        {
            L->shadow->co[i][j]=rhs.L->shadow->co[i][j];
        }
    }


    if(rhs.L->cont!=0)
    {
        L->cont=new_longmatrix(rhs.L->cont->nrh,rhs.L->cont->nch);
        for(int i=1; i<=rhs.L->cont->nrh; i++)
        {
            for(int j=1; j<=rhs.L->cont->nch; j++)
            {
                L->cont->co[i][j]=rhs.L->cont->co[i][j];
            }
        }
    }

    if(rhs.L->ty!=0)
    {
        L->ty=new_doublematrix(rhs.L->ty->nrh,rhs.L->ty->nch);
        for(int i=1; i<=rhs.L->ty->nrh; i++)
        {
            for(int j=1; j<=rhs.L->ty->nch; j++)
            {
                L->ty->co[i][j]=rhs.L->ty->co[i][j];
            }
        }
    }


    if(rhs.L->vegparp!=0)
    {
        L->vegparp=new_longmatrix(rhs.L->vegparp->nrh,rhs.L->vegparp->nch);
        for(int i=1; i<=rhs.L->vegparp->nrh; i++)
        {
            for(int j=1; j<=rhs.L->vegparp->nch; j++)
            {
                L->vegparp->co[i][j]=rhs.L->vegparp->co[i][j];
            }
        }
    }

    //change
    ///double ***vegpars;
    ///double **vegparv;
    L->vegpars=rhs.L->vegpars;


    L->vegparv=(double **)malloc((rhs.P->n_landuses+1)*sizeof(double*));
    int offset=2;
    for(int i=1; i<=rhs.P->n_landuses; i++)
    {
        //vegparp[vegetationtype][vegetationvariable] represents the column of the veg file, so vegpars[vegetationtype][time][vegparp]
        L->vegparv[i]=alloc1(jdvegprop+offset);
        for(int c=1; c<=jdvegprop; c++)
        {
            L->vegparv[i][L->vegparp->co[i][c]]=NoV;
        }

    }



    L->vegpar=new_doublevector(rhs.L->vegpar->nh);
    for(int i=1; i<=rhs.L->vegpar->nh; i++)
        L->vegpar->co[i]=rhs.L->vegpar->co[i];


    if(rhs.L->root_fraction!=0)
    {
        L->root_fraction=new_doublematrix(rhs.L->root_fraction->nrh,rhs.L->root_fraction->nch);
        for(int i=1; i<=rhs.L->root_fraction->nrh; i++)
        {
            for(int j=1; j<=rhs.L->root_fraction->nch; j++)
            {
                L->root_fraction->co[i][j]=rhs.L->root_fraction->co[i][j];
            }
        }
    }




////////////////////////////////////////Channel

    C->r=new_longvector(rhs.C->r->nh);
    for(int i=1; i<=rhs.C->r->nh; i++)
        C->r->co[i]=rhs.C->r->co[i];

    C->c=new_longvector(rhs.C->c->nh);
    for(int i=1; i<=rhs.C->c->nh; i++)
        C->c->co[i]=rhs.C->c->co[i];

    C->ch=new_longmatrix(rhs.C->ch->nrh,rhs.C->ch->nch);
    for(int i=1; i<=rhs.C->ch->nrh; i++)
    {
        for(int j=1; j<=rhs.C->ch->nch; j++)
        {
            C->ch->co[i][j]=rhs.C->ch->co[i][j];
        }
    }


    C->Qsup=new_doublevector(rhs.C->Qsup->nh);
    for(int i=1; i<=rhs.C->Qsup->nh; i++)
        C->Qsup->co[i]=rhs.C->Qsup->co[i];


    C->Qsub=new_doublevector(rhs.C->Qsub->nh);
    for(int i=1; i<=rhs.C->Qsub->nh; i++)
        C->Qsub->co[i]=rhs.C->Qsub->co[i];



    C->h_sup=new_doublevector(rhs.C->h_sup->nh);
    for(int i=1; i<=rhs.C->h_sup->nh; i++)
        C->h_sup->co[i]=rhs.C->h_sup->co[i];


    C->dh_sup=new_doublevector(rhs.C->dh_sup->nh);
    for(int i=1; i<=rhs.C->dh_sup->nh; i++)
        C->dh_sup->co[i]=rhs.C->dh_sup->co[i];

    C->hsupav=new_doublematrix(rhs.C->hsupav->nrh,rhs.C->hsupav->nch);
    for(int i=1; i<=rhs.C->hsupav->nrh; i++)
    {
        for(int j=1; j<=rhs.C->hsupav->nch; j++)
        {
            C->hsupav->co[i][j]=rhs.C->hsupav->co[i][j];
        }
    }

    C->length=new_doublevector(rhs.C->length->nh);
    for(int i=1; i<=rhs.C->length->nh; i++)
        C->length->co[i]=rhs.C->length->co[i];

    C->Q_out=rhs.C->Q_out;

////////////////////////////////////////
    /// DOUBLEMATRIX *weights_Kriging; /*dimension=npixel*nstations Kriging weights to model the rainfall distribution*/
    W->weights_Kriging=new_doublematrix(rhs.W->weights_Kriging->nrh,rhs.W->weights_Kriging->nch);
    for(int i=1; i<=rhs.W->weights_Kriging->nrh; i++)
    {
        for(int j=1; j<=rhs.W->weights_Kriging->nch; j++)
        {
            W->weights_Kriging->co[i][j]=rhs.W->weights_Kriging->co[i][j];
        }
    }


    ///DOUBLEMATRIX *PrecTot;    /*total(snow+rain) precipitation in mm (in a Dt)*/
    W->PrecTot=new_doublematrix(rhs.W->PrecTot->nrh,rhs.W->PrecTot->nch);
    for(int i=1; i<=rhs.W->PrecTot->nrh; i++)
    {
        for(int j=1; j<=rhs.W->PrecTot->nch; j++)
        {
            W->PrecTot->co[i][j]=rhs.W->PrecTot->co[i][j];
        }
    }

    ///DOUBLEMATRIX *Pnet;
    W->Pnet=new_doublematrix(rhs.W->Pnet->nrh,rhs.W->Pnet->nch);
    for(int i=1; i<=rhs.W->Pnet->nrh; i++)
    {
        for(int j=1; j<=rhs.W->Pnet->nch; j++)
        {
            W->Pnet->co[i][j]=rhs.W->Pnet->co[i][j];
        }
    }


    ///DOUBLEMATRIX *wcan_rain;       /*liquid precipitation intercepted by vegetation in mm*/
    W->wcan_rain=new_doublematrix(rhs.W->wcan_rain->nrh,rhs.W->wcan_rain->nch);
    for(int i=1; i<=rhs.W->wcan_rain->nrh; i++)
    {
        for(int j=1; j<=rhs.W->wcan_rain->nch; j++)
        {
            W->wcan_rain->co[i][j]=rhs.W->wcan_rain->co[i][j];
        }
    }

    /// DOUBLEMATRIX *wcan_snow;       /*snow intercepted by vegetation in mm*/
    W->wcan_snow=new_doublematrix(rhs.W->wcan_snow->nrh,rhs.W->wcan_snow->nch);
    for(int i=1; i<=rhs.W->wcan_snow->nrh; i++)
    {
        for(int j=1; j<=rhs.W->wcan_snow->nch; j++)
        {
            W->wcan_snow->co[i][j]=rhs.W->wcan_snow->co[i][j];
        }
    }


    //DOUBLEMATRIX *PrTOT_mean;  /*Total precipitation [mm](on nDt_output_basin Dt time intervals) total precipitation (rain+snow) precipitation*/
    if(rhs.W->PrTOT_mean!=0)
    {
        W->PrTOT_mean=new_doublematrix(rhs.W->PrTOT_mean->nrh,rhs.W->PrTOT_mean->nch);
        for(int i=1; i<=rhs.W->PrTOT_mean->nrh; i++)
        {
            for(int j=1; j<=rhs.W->PrTOT_mean->nch; j++)
            {
                W->PrTOT_mean->co[i][j]=rhs.W->PrTOT_mean->co[i][j];
            }
        }
    }

    ///DOUBLEMATRIX *PrSNW_mean; /*output of total precipitation and interception*/
    if(rhs.W->PrSNW_mean!=0)
    {
        W->PrSNW_mean=new_doublematrix(rhs.W->PrSNW_mean->nrh,rhs.W->PrSNW_mean->nch);
        for(int i=1; i<=rhs.W->PrSNW_mean->nrh; i++)
        {
            for(int j=1; j<=rhs.W->PrSNW_mean->nch; j++)
            {
                W->PrSNW_mean->co[i][j]=rhs.W->PrSNW_mean->co[i][j];
            }
        }
    }


    ///DOUBLEMATRIX *hsupav;
    W->hsupav=new_doublematrix(rhs.W->hsupav->nrh,rhs.W->hsupav->nch);
    for(int i=1; i<=rhs.W->hsupav->nrh; i++)
    {
        for(int j=1; j<=rhs.W->hsupav->nch; j++)
        {
            W->hsupav->co[i][j]=rhs.W->hsupav->co[i][j];
        }
    }

    ///DOUBLEMATRIX *dh_sup;
    W->dh_sup=new_doublematrix(rhs.W->dh_sup->nrh,rhs.W->dh_sup->nch);
    for(int i=1; i<=rhs.W->dh_sup->nrh; i++)
    {
        for(int j=1; j<=rhs.W->dh_sup->nch; j++)
        {
            W->dh_sup->co[i][j]=rhs.W->dh_sup->co[i][j];
        }
    }

    ///DOUBLEMATRIX *error;
    if(rhs.W->error!=0)
    {
        W->error=new_doublematrix(rhs.W->error->nrh,rhs.W->error->nch);
        for(int i=1; i<=rhs.W->error->nrh; i++)
        {
            for(int j=1; j<=rhs.W->error->nch; j++)
            {
                W->error->co[i][j]=rhs.W->error->co[i][j];
            }
        }
    }

    //UMFPACK_REAL_TRIPLET *Jtriplet;
    //UMFPACK_REAL_MATRIX *Jmatrix;

    ///DOUBLEVECTOR *Lx;
    W->Lx=new_doublevector(rhs.W->Lx->nh);
    for(int i=1; i<=rhs.W->Lx->nh; i++)
        W->Lx->co[i]=rhs.W->Lx->co[i];

    ///DOUBLEVECTOR *Ux;
    if(rhs.W->Ux!=0)
    {
        W->Ux=new_doublevector(rhs.W->Ux->nh);
        for(int i=1; i<=rhs.W->Ux->nh; i++)
            W->Ux->co[i]=rhs.W->Ux->co[i];
    }

    ///DOUBLEVECTOR *H0;/*potential head for each active node in the basin at the beginning of iteration, vec length is Nl*Nr*Nc */
    W->H0=new_doublevector(rhs.W->H0->nh);
    for(int i=1; i<=rhs.W->H0->nh; i++)
        W->H0->co[i]=rhs.W->H0->co[i];

    ///DOUBLEVECTOR *H1;/*potential head for each active node in the basin, vec length is Nl*Nr*Nc */
    W->H1=new_doublevector(rhs.W->H1->nh);
    for(int i=1; i<=rhs.W->H1->nh; i++)
        W->H1->co[i]=rhs.W->H1->co[i];

    ///DOUBLEVECTOR *dH;
    W->dH=new_doublevector(rhs.W->dH->nh);
    for(int i=1; i<=rhs.W->dH->nh; i++)
        W->dH->co[i]=rhs.W->dH->co[i];


    ///DOUBLEVECTOR *B;
    W->B=new_doublevector(rhs.W->B->nh);
    for(int i=1; i<=rhs.W->B->nh; i++)
        W->B->co[i]=rhs.W->B->co[i];


    ///DOUBLEVECTOR *f;
    W->f=new_doublevector(rhs.W->f->nh);
    for(int i=1; i<=rhs.W->f->nh; i++)
        W->f->co[i]=rhs.W->f->co[i];

    //DOUBLEVECTOR *df;
    W->df=new_doublevector(rhs.W->df->nh);
    for(int i=1; i<=rhs.W->df->nh; i++)
        W->df->co[i]=rhs.W->df->co[i];


////////////////////////////////////////Times
    I->iter=rhs.I->iter;    /*current iteration of a time interval for egy-mass balance*/
    I->n_iter=rhs.I->n_iter;  /*n_iter=number of iterations to do for each time-step*/
    I->TH=rhs.I->TH;     /*TH=The numbers of days of simulation*/
    I->i_pixel=rhs.I->i_pixel;  /*counter for the output of a pixel*/
    I->n_pixel=rhs.I->n_pixel;  /*nDt_output_pixel=number of Dt after which the output of a pixel are printed*/
    I->i_basin=rhs.I->i_basin;  /*counter for the output of a pixel*/
    I->n_basin=rhs.I->n_basin;  /*THE NUMBER OF Dt AFTER WHICH THE OUTPUT FOR THE BASIN ARE PRINTED*/
    I->i_discharge=rhs.I->i_discharge;/*counter for the output of a pixel*/
    I->n_discharge=rhs.I->n_discharge;/*THE NUMBER OF Dt AFTER WHICH THE DISCHARGE IS PRINTED*/
    I->i_plot=rhs.I->i_plot;
    I->n_plot=rhs.I->n_plot;
    I->nt_plot=rhs.I->nt_plot;
    I->d_plot=rhs.I->d_plot;
    I->JD=rhs.I->JD;      /*day=current Julian day during the simulation*/
    I->day=rhs.I->day;       /*current day of the month during the simulation*/
    I->month=rhs.I->month;       /*current month of the year during the simulation*/
    I->year=rhs.I->year;     /*current year*/
    I->hour=rhs.I->hour;       /*current hour of day during the simulation*/
    I->min=rhs.I->min;       /*current minute of hour during the simulation*/
    I->time=rhs.I->time;    /*time=current time from the begin of simulation [s]*/
    I->egy=rhs.I->egy;  /*the time of egy subroutine [s]*/
    I->vert_wb=rhs.I->vert_wb; /*the time of water-vertical-distribution subroutine [s]*/
    I->horiz_wb=rhs.I->horiz_wb;/*the time of water-horizontal-distribution subroutine [s]*/
    I->writeout=rhs.I->writeout;/*the time of write-output subroutine [s]*/

////////////////////////////////////////Par
    P->Dt=rhs.P->Dt;      /*Dt=the integration time interval [s]*/
    P->JD0=rhs.P->JD0;		/*Decimal julian day at the beginning of simulation (0.0-365.99)*/
    P->year0=rhs.P->year0;		/*Year at the beginning of simulation*/
    P->ST=rhs.P->ST;		/*Standard time to which all the output data are referred (difference respect UMT, in hour)*/
    P->print=rhs.P->print;         /*1 IF YOU WANT TO PRINT MATRICES WITH INTERMEDIATE RESULTS, 0 OTHERWISE*/
    P->monin_obukhov=rhs.P->monin_obukhov;  /**#2      monin_obuhkov    1 stability and instability considered, 2 stability not considered, 3 instability not considered, 4 always neutrality,
			 5 Oerlemans&Grisogono (2002) on glacier*/

    P->gamma_m=rhs.P->gamma_m;   /**5 Exponent of the law of uniform motion on the surface*/
    P->T_rain=rhs.P->T_rain;    /*TEMPERATURE ABOVE WICH ALL PRECIPITAION IS RAIN [C]*/
    P->T_snow=rhs.P->T_snow;    /*TEMPERATURE BELOW WICH ALL PRECIPITAION IS SNOW [C]*/
    P->aep=rhs.P->aep;       /*ALBEDO EXTINCTION PARAMETER [m] if snow depth < aep, albedo is interpolated between soil and snow*/
    P->avo=rhs.P->avo;       /*NEW SNOW VISIBLE BAND REFLECTANCE*/
    P->airo=rhs.P->airo;      /*NEW NEAR INFRARED BAND REFLECTANCE*/
    P->Sr=rhs.P->Sr;		  /*WATER FRACTION RETAINED BY CAPILLARY FORCES IN SNOW IRREDUCIBLE WATER SATURATION [-] - from Colbeck (0.02 - 0.07)*/
    P->rho_ice=rhs.P->rho_ice;     /*Ice density [kg/mc]*/
    P->total_pixel=rhs.P->total_pixel;    /*The number of the valid pixel of the whole basin*/
    P->total_channel=rhs.P->total_channel;	/**The total number of channel-pixel in the basin */
    P->snowlayer_max=rhs.P->snowlayer_max; /*MAXIMUM NUMBER OF SNOW LAYERS*/
    P->snowlayer_inf=rhs.P->snowlayer_inf; /*layer of illimitate thickness (beginning from below)*/


    P->Dmin=new_doublevector(rhs.P->Dmin->nh);
    for(int i=1; i<=rhs.P->Dmin->nh; i++)
        P->Dmin->co[i]=rhs.P->Dmin->co[i];

    P->Dmax=new_doublevector(rhs.P->Dmax->nh);
    for(int i=1; i<=rhs.P->Dmax->nh; i++)
        P->Dmax->co[i]=rhs.P->Dmax->co[i];


    P->Sr_glac=rhs.P->Sr_glac;
    P->glaclayer_max=rhs.P->glaclayer_max;



    if(rhs.P->Dmin_glac!=0)
    {
        P->Dmin_glac=new_doublevector(rhs.P->Dmin_glac->nh);
        for(int i=1; i<=rhs.P->Dmin_glac->nh; i++)
            P->Dmin_glac->co[i]=rhs.P->Dmin_glac->co[i];
    }

    if(rhs.P->Dmax_glac!=0)
    {
        P->Dmax_glac=new_doublevector(rhs.P->Dmax_glac->nh);
        for(int i=1; i<=rhs.P->Dmax_glac->nh; i++)
            P->Dmax_glac->co[i]=rhs.P->Dmax_glac->co[i];
    }

    P->state_turb=rhs.P->state_turb;

    /** 9 block - BASE AND ADVANCED PARAMETERS
    #1      state_lwrad	 Which formula for incoming longwave radiation:
    	1 (Brutsaert, 1975), 2 (Satterlund, 1979), 3 (Idso, 1981), 4(Idso+Hodges),
    	5 (Koenig-Langlo & Augstein, 1994), 6 (Andreas & Ackley, 1982), 7 (Konzelmann, 1994),
    	8 (Prata, 1996), 9 (Dilley 1998)*/
    P->state_lwrad=rhs.P->state_lwrad;

    P->imp=rhs.P->imp;/**1 Impedence factor for (partially) frozen soil  1.Reduction factor of the hydraulic conductivity in partially frozen soil (K=K_no_ice*10^(impedence*Q), where Q is the ice ratio)
	suggested values: from 0 to 7, best value 4*/
    P->f_bound_Richards=rhs.P->f_bound_Richards;

    P->epsilon_snow=rhs.P->epsilon_snow;/**SNOW LONGWAVE EMISSIVITY [-]*/

    P->output_Txy=rhs.P->output_Txy;
    P->output_TETAxy=rhs.P->output_TETAxy;
    P->output_TETAICExy=rhs.P->output_TETAICExy;
    P->output_PSIxy=rhs.P->output_PSIxy;
    P->output_snow=rhs.P->output_snow;
    P->output_glac=rhs.P->output_glac;
    P->output_h_sup=rhs.P->output_h_sup;
    P->output_Rn=rhs.P->output_Rn;
    P->output_G=rhs.P->output_G;
    P->output_H=rhs.P->output_H;
    P->output_ET=rhs.P->output_ET;
    P->output_Ts=rhs.P->output_Ts;
    P->output_P=rhs.P->output_P;
    P->output_Wr=rhs.P->output_Wr;
    P->output_balancesn=rhs.P->output_balancesn;
    P->output_balancegl=rhs.P->output_balancegl;
    P->output_Rswdown=rhs.P->output_Rswdown;
    P->output_meteo=rhs.P->output_meteo;




    /**
     * chkpt COORDINATES of the points for which the simulations is printed
     * rc Rows and Cols of the points for which the simulations is printed
     */
    //DOUBLEMATRIX *chkpt;
    P->chkpt=new_doublematrix(rhs.P->chkpt->nrh,rhs.P->chkpt->nch);
    for(int i=1; i<=rhs.P->chkpt->nrh; i++)
    {
        for(int j=1; j<=rhs.P->chkpt->nch; j++)
        {
            P->chkpt->co[i][j]=rhs.P->chkpt->co[i][j];
        }
    }

    //LONGMATRIX *rc;
    P->rc=new_longmatrix(rhs.P->rc->nrh,rhs.P->rc->nch);
    for(int i=1; i<=rhs.P->rc->nrh; i++)
    {
        for(int j=1; j<=rhs.P->rc->nch; j++)
        {
            P->rc->co[i][j]=rhs.P->rc->co[i][j];
        }
    }


    /**
     * state_px_coord	 1 IF ALL COORDINATES ARE IN FORMAT (EAST,NORTH), 0 IF IN FORMAT ROW AND COLUMS (r,c)
     */
    P->state_px_coord=rhs.P->state_px_coord;

    P->integr_scale_rain=rhs.P->integr_scale_rain;
    P->variance_rain=rhs.P->variance_rain;

    P->recover=rhs.P->recover;		/*1 if you want to recover a simulation, 0 otherwise*/

    P->Vmin=rhs.P->Vmin;/**MINIMUM WIND VELOCITY (m/s)*/

    P->snowcorrfact=rhs.P->snowcorrfact;/*factor multiplying the snow precipitation*/
    P->raincorrfact=rhs.P->raincorrfact;/*factor multiplying the rain precipitation*/

    P->RHmin=rhs.P->RHmin;/**MINIMUM RELATIVE HUMIDITY VALUE (%)*/

    P->format_out=rhs.P->format_out;	/*OUTPUT MAPS in fluidturtle format (=1), GRASS ASCII (=2), ESRI ASCII (=3)*/
    P->nsky=rhs.P->nsky;			/*Multiplying factor decreasing the dem resolution for the calculation of the sky view factor*/
    P->channel_thres=rhs.P->channel_thres;/*Value of the threshold for definition of pixel channel [1.0-1000.0] (used if the network map is not provided)*/

    /**
     * saving_points saving points expressed as numbers of day after the beginning
     */
    ///DOUBLEVECTOR *saving_points;
    P->saving_points=new_doublevector(rhs.P->saving_points->nh);
    for(int i=1; i<=rhs.P->saving_points->nh; i++)
        P->saving_points->co[i]=rhs.P->saving_points->co[i];

    //double Vis; //visibility in km (>5 km)
    //double Lozone; //thickness of the stratospheric ozone layer (in cm normal conditions)

    P->point_sim=rhs.P->point_sim;/* =0 distributed simulation, =1 point simulation (the parameter files are different in the two cases) */

    P->snow_maxpor=rhs.P->snow_maxpor;/*MAXIMUM SNOW POROSITY ALLOWED*/
    P->snow_density_cutoff=rhs.P->snow_density_cutoff;
    P->drysnowdef_rate=rhs.P->drysnowdef_rate;/*SNOW COMPACTION (% per hour) DUE TO DESTRUCTIVE METAMORPHISM for SNOW DENSITY<snow_density_cutoff and DRY SNOW */
    P->wetsnowdef_rate=rhs.P->wetsnowdef_rate;/*ENHANCEMENT FACTOR IN PRESENCE OF WET SNOW*/
    P->snow_viscosity=rhs.P->snow_viscosity; /*SNOW VISCOSITY COEFFICIENT (kg s m^-2) at T=0 C and snow density=0*/

    P->latitude=rhs.P->latitude;/**latitude in radiance*/
    P->longitude=rhs.P->longitude;

    P->z0_snow=rhs.P->z0_snow;/*roughness length over snow (mm)*/
    P->n_landuses=rhs.P->n_landuses;/**The Maximum number of land use types*/

    /**JD_plots It is possible to select some days for which energy balance and
     * meteo data are plotted with a very short time step.
     */
    //DOUBLEVECTOR *JD_plots;
    P->JD_plots=new_doublevector(rhs.P->JD_plots->nh);
    for(int i=1; i<=rhs.P->JD_plots->nh; i++)
        P->JD_plots->co[i]=rhs.P->JD_plots->co[i];

    P->micromet=rhs.P->micromet;/**Use Micromet (=1), otherwise (=0)*/
    P->blowing_snow=rhs.P->blowing_snow;/**PBSM	 Use PBSM (snow horizontal transport) (=1), otherwise (=0)*/

    ///LONGVECTOR *r_points;
    ///LONGVECTOR *c_points;
    if(rhs.P->r_points!=0)
    {
        P->r_points=new_longvector(rhs.P->r_points->nh);
        for(int i=1; i<=rhs.P->r_points->nh; i++)
            P->r_points->co[i]=rhs.P->r_points->co[i];
    }

    if(rhs.P->c_points!=0)
    {
        P->c_points=new_longvector(rhs.P->c_points->nh);
        for(int i=1; i<=rhs.P->c_points->nh; i++)
            P->c_points->co[i]=rhs.P->c_points->co[i];
    }

    P->psimin=rhs.P->psimin;
    P->stmin=rhs.P->stmin;

    /**wat_balance	water balance 	 1 calculate, 0 do not calculate
     * en_balance	energy balance   1 calculate, 0 do not calculate
     */
    P->wat_balance=rhs.P->wat_balance;
    P->en_balance=rhs.P->en_balance;

    P->distr_stat=rhs.P->distr_stat;

    P->Dpsi=rhs.P->Dpsi;
    P->dtmin=rhs.P->dtmin;

    P->PsiInf=rhs.P->PsiInf;
    P->TolPsiInf=rhs.P->TolPsiInf;

    P->q1=rhs.P->q1;
    P->q2=rhs.P->q2;


    P->transect=NULL;
    P->vtrans=NULL;


    if(rhs.P->cont_trans!=0)
    {
        P->cont_trans=new_longvector(rhs.P->cont_trans->nh);
        for(int i=1; i<=rhs.P->cont_trans->nh; i++)
            P->cont_trans->co[i]=rhs.P->cont_trans->co[i];
    }


    if(rhs.P->ibeg!=0)
    {
        P->ibeg=new_longvector(rhs.P->ibeg->nh);
        for(int i=1; i<=rhs.P->ibeg->nh; i++)
            P->ibeg->co[i]=rhs.P->ibeg->co[i];
    }

    P->nLC=rhs.P->nLC;

    P->snow_fetch=rhs.P->snow_fetch;/*MINIMUM HORIZONTAL FETCH THAT ALLOWS FULL SNOW TRANSPORT*/

    P->ifill=rhs.P->ifill;
    P->iobsint=rhs.P->iobsint;
    P->dn=rhs.P->dn;
    P->curve_len_scale=rhs.P->curve_len_scale;/**__parameter.txt block 3 16*/
    P->slopewt=rhs.P->slopewt;
    P->curvewt=rhs.P->curvewt;
    P->topoflag=rhs.P->topoflag;/**if=1 top->Zm->co[r][c] = top->Z0->co[r][c] + DEPTH(r, c, snow->lnum, snow->Dzl);*/

    P->LRflag=rhs.P->LRflag;/**LRflag=1,has a lapse rates file;LRflag=0,without a lapse rates file*/


    P->vegflag=new_shortvector(rhs.P->vegflag->nh);
    for(int i=1; i<=rhs.P->vegflag->nh; i++)
        P->vegflag->co[i]=rhs.P->vegflag->co[i];



    P->harm_or_arit_mean=rhs.P->harm_or_arit_mean;
    P->MaxiterTol=rhs.P->MaxiterTol;/**Max iterations for the integration of Richards' equation*/
    P->MaxErrWb=rhs.P->MaxErrWb;/**MaximalErrorRichards*/
    P->TolVWb=rhs.P->TolVWb;/**2  Absolute Tolerance for the integration of Richards' equation*/

    P->incr_V_factor=rhs.P->incr_V_factor;/*Factor multiplying the averaged wind speed in order to run the blowing snow model*/

    P->alpha_snow=rhs.P->alpha_snow;
    P->tol_energy=rhs.P->tol_energy;
    P->maxiter_energy=rhs.P->maxiter_energy;
    P->maxiter_canopy=rhs.P->maxiter_canopy;
    P->maxiter_Ts=rhs.P->maxiter_Ts;
    P->maxiter_Loc=rhs.P->maxiter_Loc;
    P->maxiter_Businger=rhs.P->maxiter_Businger;
    P->stabcorr_incanopy=rhs.P->stabcorr_incanopy;

    P->state_pixel=rhs.P->state_pixel;/* if n_pixel>0, state_pixel=1,else state_pixel=0 */
    P->state_basin=rhs.P->state_basin;/* if n_basin>0, state_basin=1,else state_basin=0 */
    P->state_discharge=rhs.P->state_discharge;/* if n_discharge>0, state_discharge=1,else state_discharge=0 */

    P->DtminWb=rhs.P->DtminWb;
    P->nredDtWb=rhs.P->nredDtWb;

    P->TolCG=rhs.P->TolCG;/**Initial forcing term of Newton method */
    P->MaxiterCorr=rhs.P->MaxiterCorr;
    P->UpdateK=rhs.P->UpdateK;

    P->channel_network=rhs.P->channel_network;
    P->bedrock=rhs.P->bedrock;/** indicate whether the bedrock file is input,Yes=1,No=0 */

    P->thres_hsup=rhs.P->thres_hsup;/**Threshold on h_sup [mm] above which C_m is independent from h_sup*/
    P->thres_hchannel=rhs.P->thres_hchannel;

    //double cm_hsup_0;
    //double max_Courant_sup_sub;

    P->Kch_b=rhs.P->Kch_b;/**hydraulic conductivity of the sediments divided by sediment thickness*/
    P->w_dx=rhs.P->w_dx;/**channel width / pixel width*/

    P->RelTolVWb=rhs.P->RelTolVWb;/**RelativeErrorRichards*/

    P->snow_smin=rhs.P->snow_smin;/*MINIMUM SLOPE [degree] TO ADJUST PRECIPITATION REDUCTION*/
    P->snow_smax=rhs.P->snow_smax;/*MAXIMUM SLOPE [degree] TO ADJUST PRECIPITATION REDUCTION*/
    P->snow_curv=rhs.P->snow_curv;/*SHAPE PARAMETER FOR PRECIPITATION REDUCTION,correction of snow precipitation in case of steep slopes (contribution by Stephan Gruber)*/

    P->Zboundary=rhs.P->Zboundary;/**0 temperature amplitude depth in the soil [mm]*/
    P->Tboundary=rhs.P->Tboundary;/**temperature at the 0 temperature amplitude depth [C]*/

    P->Ks_channel=rhs.P->Ks_channel;/**Cm coefficient for the channel flow (the same gamma for surface flow is used)*/
    P->depr_channel=rhs.P->depr_channel;/**Depression of the channel bed with respect to the neighbouring slopes*/
////////////////////////////////////////SNOW
    ///SHORTMATRIX  *type;/**/
    N->type=new_shortmatrix(rhs.N->type->nrh,rhs.N->type->nch);
    for(int i=1; i<=rhs.N->type->nrh; i++)
    {
        for(int j=1; j<=rhs.N->type->nch; j++)
        {
            N->type->co[i][j]=rhs.N->type->co[i][j];
        }
    }
    ///LONGMATRIX	 *lnum;/*snow layer number*/
    N->lnum=new_longmatrix(rhs.N->lnum->nrh,rhs.N->lnum->nch);
    for(int i=1; i<=rhs.N->lnum->nrh; i++)
    {
        for(int j=1; j<=rhs.N->lnum->nch; j++)
        {
            N->lnum->co[i][j]=rhs.N->lnum->co[i][j];
        }
    }

    ///DOUBLETENSOR *Dzl;/**snow thickness for each layer*/
    ///DOUBLETENSOR *w_liq;/**liquid water content for each snow layer*/
    ///DOUBLETENSOR *w_ice;/**ice content for each snow layer*/
    ///DOUBLETENSOR *T;/*Temperature for each snow layer*/

    N->Dzl=new_doubletensor(rhs.N->Dzl->ndh, rhs.N->Dzl->nrh, rhs.N->Dzl->nch);/** contains the soil parameters from the __soil.txt, pa[SoilType][Property][Layer] */
    initialize_doubletensor(N->Dzl, 0.0);

    N->w_liq=new_doubletensor(rhs.N->w_liq->ndh, rhs.N->w_liq->nrh, rhs.N->w_liq->nch);/** contains the soil parameters from the __soil.txt, pa[SoilType][Property][Layer] */
    initialize_doubletensor(N->w_liq, 0.0);

    N->w_ice=new_doubletensor(rhs.N->w_ice->ndh, rhs.N->w_ice->nrh, rhs.N->w_ice->nch);/** contains the soil parameters from the __soil.txt, pa[SoilType][Property][Layer] */
    initialize_doubletensor(N->w_ice, 0.0);

    N->T=new_doubletensor(rhs.N->T->ndh, rhs.N->T->nrh, rhs.N->T->nch);/** contains the soil parameters from the __soil.txt, pa[SoilType][Property][Layer] */
    initialize_doubletensor(N->T, 0.0);

    for(int i=1; i<=rhs.N->Dzl->ndh; i++)
    {
        for(int j=1; j<=rhs.N->Dzl->nrh; j++)
        {
            for(int k=1; k<=rhs.N->Dzl->nch; k++)
            {
                N->Dzl->co[i][j][k]=rhs.N->Dzl->co[i][j][k];
            }
        }

    }
    for(int i=1; i<=rhs.N->w_liq->ndh; i++)
    {
        for(int j=1; j<=rhs.N->w_liq->nrh; j++)
        {
            for(int k=1; k<=rhs.N->w_liq->nch; k++)
            {
                N->w_liq->co[i][j][k]=rhs.N->w_liq->co[i][j][k];

            }
        }

    }
    for(int i=1; i<=rhs.N->w_ice->ndh; i++)
    {
        for(int j=1; j<=rhs.N->w_ice->nrh; j++)
        {
            for(int k=1; k<=rhs.N->w_ice->nch; k++)
            {

                N->w_ice->co[i][j][k]=rhs.N->w_ice->co[i][j][k];

            }
        }

    }
    for(int i=1; i<=rhs.N->T->ndh; i++)
    {
        for(int j=1; j<=rhs.N->T->nrh; j++)
        {
            for(int k=1; k<=rhs.N->T->nch; k++)
            {
                N->T->co[i][j][k]=rhs.N->T->co[i][j][k];
            }
        }

    }



    //DOUBLEMATRIX *nondimens_age;
    if(rhs.N->nondimens_age!=0)
    {
        N->nondimens_age=new_doublematrix(rhs.N->nondimens_age->nrh,rhs.N->nondimens_age->nch);
        for(int i=1; i<=rhs.N->nondimens_age->nrh; i++)
        {
            for(int j=1; j<=rhs.N->nondimens_age->nch; j++)
            {
                N->nondimens_age->co[i][j]=rhs.N->nondimens_age->co[i][j];
            }
        }
    }

    //DOUBLEMATRIX *dimens_age;
    N->dimens_age=new_doublematrix(rhs.N->dimens_age->nrh,rhs.N->dimens_age->nch);
    for(int i=1; i<=rhs.N->dimens_age->nrh; i++)
    {
        for(int j=1; j<=rhs.N->dimens_age->nch; j++)
        {
            N->dimens_age->co[i][j]=rhs.N->dimens_age->co[i][j];
        }
    }

    //DOUBLEMATRIX *max;
    N->max=new_doublematrix(rhs.N->max->nrh,rhs.N->max->nch);
    for(int i=1; i<=rhs.N->max->nrh; i++)
    {
        for(int j=1; j<=rhs.N->max->nch; j++)
        {
            N->max->co[i][j]=rhs.N->max->co[i][j];
        }
    }

    //DOUBLEMATRIX *average;
    N->average=new_doublematrix(rhs.N->average->nrh,rhs.N->average->nch);
    for(int i=1; i<=rhs.N->average->nrh; i++)
    {
        for(int j=1; j<=rhs.N->average->nch; j++)
        {
            N->average->co[i][j]=rhs.N->average->co[i][j];
        }
    }

    //DOUBLEMATRIX *MELTED;
    if(rhs.N->MELTED!=0)
    {
        N->MELTED=new_doublematrix(rhs.N->MELTED->nrh,rhs.N->MELTED->nch);
        for(int i=1; i<=rhs.N->MELTED->nrh; i++)
        {
            for(int j=1; j<=rhs.N->MELTED->nch; j++)
            {
                N->MELTED->co[i][j]=rhs.N->MELTED->co[i][j];
            }
        }
    }


    if(rhs.N->SUBL!=0)
    {
        N->SUBL=new_doublematrix(rhs.N->SUBL->nrh,rhs.N->SUBL->nch);
        for(int i=1; i<=rhs.N->SUBL->nrh; i++)
        {
            for(int j=1; j<=rhs.N->SUBL->nch; j++)
            {
                N->SUBL->co[i][j]=rhs.N->SUBL->co[i][j];
            }
        }
    }

    //DOUBLEMATRIX *t_snow;
    if(rhs.N->t_snow!=0)
    {
        N->t_snow=new_doublematrix(rhs.N->t_snow->nrh,rhs.N->t_snow->nch);
        for(int i=1; i<=rhs.N->t_snow->nrh; i++)
        {
            for(int j=1; j<=rhs.N->t_snow->nch; j++)
            {
                N->t_snow->co[i][j]=rhs.N->t_snow->co[i][j];
            }
        }
    }

    //DOUBLEMATRIX *DDF;
    if(rhs.N->DDF!=0)
    {
        N->DDF=new_doublematrix(rhs.N->DDF->nrh,rhs.N->DDF->nch);
        for(int i=1; i<=rhs.N->DDF->nrh; i++)
        {
            for(int j=1; j<=rhs.N->DDF->nch; j++)
            {
                N->DDF->co[i][j]=rhs.N->DDF->co[i][j];
            }
        }
    }

    //DOUBLEMATRIX *DDF1;
    if(rhs.N->DDF1!=0)
    {
        N->DDF1=new_doublematrix(rhs.N->DDF1->nrh,rhs.N->DDF1->nch);
        for(int i=1; i<=rhs.N->DDF1->nrh; i++)
        {
            for(int j=1; j<=rhs.N->DDF1->nch; j++)
            {
                N->DDF1->co[i][j]=rhs.N->DDF1->co[i][j];
            }
        }
    }

    //DOUBLEMATRIX *DDFvar;
    if(rhs.N->DDFvar!=0)
    {
        N->DDFvar=new_doublematrix(rhs.N->DDFvar->nrh,rhs.N->DDFvar->nch);
        for(int i=1; i<=rhs.N->DDFvar->nrh; i++)
        {
            for(int j=1; j<=rhs.N->DDFvar->nch; j++)
            {
                N->DDFvar->co[i][j]=rhs.N->DDFvar->co[i][j];
            }
        }
    }

    //LONGMATRIX   *DDFcont;
    if(rhs.N->DDFcont!=0)
    {
        N->DDFcont=new_longmatrix(rhs.N->DDFcont->nrh,rhs.N->DDFcont->nch);
        for(int i=1; i<=rhs.N->DDFcont->nrh; i++)
        {
            for(int j=1; j<=rhs.N->DDFcont->nch; j++)
            {
                N->DDFcont->co[i][j]=rhs.N->DDFcont->co[i][j];
            }
        }
    }

    //DOUBLEMATRIX *DDFTmin;
    if(rhs.N->DDFTmin!=0)
    {
        N->DDFTmin=new_doublematrix(rhs.N->DDFTmin->nrh,rhs.N->DDFTmin->nch);
        for(int i=1; i<=rhs.N->DDFTmin->nrh; i++)
        {
            for(int j=1; j<=rhs.N->DDFTmin->nch; j++)
            {
                N->DDFTmin->co[i][j]=rhs.N->DDFTmin->co[i][j];
            }
        }
    }

    //DOUBLEMATRIX *DDFmeltTL0;
    if(rhs.N->DDFmeltTL0!=0)
    {
        N->DDFmeltTL0=new_doublematrix(rhs.N->DDFmeltTL0->nrh,rhs.N->DDFmeltTL0->nch);
        for(int i=1; i<=rhs.N->DDFmeltTL0->nrh; i++)
        {
            for(int j=1; j<=rhs.N->DDFmeltTL0->nch; j++)
            {
                N->DDFmeltTL0->co[i][j]=rhs.N->DDFmeltTL0->co[i][j];
            }
        }
    }

    //DOUBLEMATRIX *DDFmelt;
    if(rhs.N->DDFmelt!=0)
    {
        N->DDFmelt=new_doublematrix(rhs.N->DDFmelt->nrh,rhs.N->DDFmelt->nch);
        for(int i=1; i<=rhs.N->DDFmelt->nrh; i++)
        {
            for(int j=1; j<=rhs.N->DDFmelt->nch; j++)
            {
                N->DDFmelt->co[i][j]=rhs.N->DDFmelt->co[i][j];
            }
        }
    }


    //DOUBLEMATRIX *rho_newsnow;
    if(rhs.N->rho_newsnow!=0)
    {
        N->rho_newsnow=new_doublematrix(rhs.N->rho_newsnow->nrh,rhs.N->rho_newsnow->nch);
        for(int i=1; i<=rhs.N->rho_newsnow->nrh; i++)
        {
            for(int j=1; j<=rhs.N->rho_newsnow->nch; j++)
            {
                N->rho_newsnow->co[i][j]=rhs.N->rho_newsnow->co[i][j];
            }
        }
    }

    //DOUBLEMATRIX *Qsub;
    if(rhs.N->Qsub!=0)
    {
        N->Qsub=new_doublematrix(rhs.N->Qsub->nrh,rhs.N->Qsub->nch);
        for(int i=1; i<=rhs.N->Qsub->nrh; i++)
        {
            for(int j=1; j<=rhs.N->Qsub->nch; j++)
            {
                N->Qsub->co[i][j]=rhs.N->Qsub->co[i][j];
            }
        }
    }

    ///DOUBLEMATRIX *Wtrans;
    if(rhs.N->Wtrans!=0)
    {
        N->Wtrans=new_doublematrix(rhs.N->Wtrans->nrh,rhs.N->Wtrans->nch);
        for(int i=1; i<=rhs.N->Wtrans->nrh; i++)
        {
            for(int j=1; j<=rhs.N->Wtrans->nch; j++)
            {
                N->Wtrans->co[i][j]=rhs.N->Wtrans->co[i][j];
            }
        }
    }

    if(rhs.N->Qtrans!=0)
    {
        N->Qtrans=new_doublematrix(rhs.N->Qtrans->nrh,rhs.N->Qtrans->nch);
        for(int i=1; i<=rhs.N->Qtrans->nrh; i++)
        {
            for(int j=1; j<=rhs.N->Qtrans->nch; j++)
            {
                N->Qtrans->co[i][j]=rhs.N->Qtrans->co[i][j];
            }
        }
    }

    ///DOUBLEMATRIX *Qtrans_x;
    if(rhs.N->Qtrans_x!=0)
    {
        N->Qtrans_x=new_doublematrix(rhs.N->Qtrans_x->nrh,rhs.N->Qtrans_x->nch);
        for(int i=1; i<=rhs.N->Qtrans_x->nrh; i++)
        {
            for(int j=1; j<=rhs.N->Qtrans_x->nch; j++)
            {
                N->Qtrans_x->co[i][j]=rhs.N->Qtrans_x->co[i][j];
            }
        }
    }

    ///DOUBLEMATRIX *Qtrans_y;
    if(rhs.N->Qtrans_y!=0)
    {
        N->Qtrans_y=new_doublematrix(rhs.N->Qtrans_y->nrh,rhs.N->Qtrans_y->nch);
        for(int i=1; i<=rhs.N->Qtrans_y->nrh; i++)
        {
            for(int j=1; j<=rhs.N->Qtrans_y->nch; j++)
            {
                N->Qtrans_y->co[i][j]=rhs.N->Qtrans_y->co[i][j];
            }
        }
    }

    ///DOUBLEMATRIX *Wtot;
    if(rhs.N->Wtot!=0)
    {
        N->Wtot=new_doublematrix(rhs.N->Wtot->nrh,rhs.N->Wtot->nch);
        for(int i=1; i<=rhs.N->Wtot->nrh; i++)
        {
            for(int j=1; j<=rhs.N->Wtot->nch; j++)
            {
                N->Wtot->co[i][j]=rhs.N->Wtot->co[i][j];
            }
        }
    }

    ///DOUBLEMATRIX *Wsubl_cum;
    if(rhs.N->Wsubl_cum!=0)
    {
        N->Wsubl_cum=new_doublematrix(rhs.N->Wsubl_cum->nrh,rhs.N->Wsubl_cum->nch);
        for(int i=1; i<=rhs.N->Wsubl_cum->nrh; i++)
        {
            for(int j=1; j<=rhs.N->Wsubl_cum->nch; j++)
            {
                N->Wsubl_cum->co[i][j]=rhs.N->Wsubl_cum->co[i][j];
            }
        }
    }

    ///DOUBLEMATRIX *Wsusp_cum;
    if(rhs.N->Wsusp_cum!=0)
    {
        N->Wsusp_cum=new_doublematrix(rhs.N->Wsusp_cum->nrh,rhs.N->Wsusp_cum->nch);
        for(int i=1; i<=rhs.N->Wsusp_cum->nrh; i++)
        {
            for(int j=1; j<=rhs.N->Wsusp_cum->nch; j++)
            {
                N->Wsusp_cum->co[i][j]=rhs.N->Wsusp_cum->co[i][j];
            }
        }
    }

    ///DOUBLEMATRIX *Wtrans_cum;
    if(rhs.N->Wtrans_cum!=0)
    {
        N->Wtrans_cum=new_doublematrix(rhs.N->Wtrans_cum->nrh,rhs.N->Wtrans_cum->nch);
        for(int i=1; i<=rhs.N->Wtrans_cum->nrh; i++)
        {
            for(int j=1; j<=rhs.N->Wtrans_cum->nch; j++)
            {
                N->Wtrans_cum->co[i][j]=rhs.N->Wtrans_cum->co[i][j];
            }
        }
    }

    if(rhs.N->Wsubgrid_cum!=0)
    {
        N->Wsubgrid_cum=new_doublematrix(rhs.N->Wsubgrid_cum->nrh,rhs.N->Wsubgrid_cum->nch);
        for(int i=1; i<=rhs.N->Wsubgrid_cum->nrh; i++)
        {
            for(int j=1; j<=rhs.N->Wsubgrid_cum->nch; j++)
            {
                N->Wsubgrid_cum->co[i][j]=rhs.N->Wsubgrid_cum->co[i][j];
            }
        }
    }

    ///DOUBLEMATRIX *out_bs;
    N->out_bs=new_doublematrix(rhs.N->out_bs->nrh,rhs.N->out_bs->nch);
    for(int i=1; i<=rhs.N->out_bs->nrh; i++)
    {
        for(int j=1; j<=rhs.N->out_bs->nch; j++)
        {
            N->out_bs->co[i][j]=rhs.N->out_bs->co[i][j];
        }
    }

    ///DOUBLEMATRIX *ListonSWE;
    if(rhs.N->ListonSWE!=0)
    {
        N->ListonSWE=new_doublematrix(rhs.N->ListonSWE->nrh,rhs.N->ListonSWE->nch);
        for(int i=1; i<=rhs.N->ListonSWE->nrh; i++)
        {
            for(int j=1; j<=rhs.N->ListonSWE->nch; j++)
            {
                N->ListonSWE->co[i][j]=rhs.N->ListonSWE->co[i][j];
            }
        }
    }

    ///DOUBLEMATRIX *softSWE;
    if(rhs.N->softSWE!=0)
    {
        N->softSWE=new_doublematrix(rhs.N->softSWE->nrh,rhs.N->softSWE->nch);
        for(int i=1; i<=rhs.N->softSWE->nrh; i++)
        {
            for(int j=1; j<=rhs.N->softSWE->nch; j++)
            {
                N->softSWE->co[i][j]=rhs.N->softSWE->co[i][j];
            }
        }
    }

    ///DOUBLEMATRIX *softSWE1;
    if(rhs.N->softSWE1!=0)
    {
        N->softSWE1=new_doublematrix(rhs.N->softSWE1->nrh,rhs.N->softSWE1->nch);
        for(int i=1; i<=rhs.N->softSWE1->nrh; i++)
        {
            for(int j=1; j<=rhs.N->softSWE1->nch; j++)
            {
                N->softSWE1->co[i][j]=rhs.N->softSWE1->co[i][j];
            }
        }
    }

    ///DOUBLEMATRIX *Dplot;
    if(rhs.N->Dplot!=0)
    {
        N->Dplot=new_doublematrix(rhs.N->Dplot->nrh,rhs.N->Dplot->nch);
        for(int i=1; i<=rhs.N->Dplot->nrh; i++)
        {
            for(int j=1; j<=rhs.N->Dplot->nch; j++)
            {
                N->Dplot->co[i][j]=rhs.N->Dplot->co[i][j];
            }
        }
    }


    ///DOUBLEVECTOR *CR1;
    ///DOUBLEVECTOR *CR2;
    ///DOUBLEVECTOR *CR3;
    ///LONGVECTOR *change_dir_wind;

    N->CR1=new_doublevector(rhs.N->CR1->nh);
    for(int i=1; i<=rhs.N->CR1->nh; i++)
        N->CR1->co[i]=rhs.N->CR1->co[i];

    N->CR2=new_doublevector(rhs.N->CR2->nh);
    for(int i=1; i<=rhs.N->CR2->nh; i++)
        N->CR2->co[i]=rhs.N->CR2->co[i];

    N->CR3=new_doublevector(rhs.N->CR3->nh);
    for(int i=1; i<=rhs.N->CR3->nh; i++)
        N->CR3->co[i]=rhs.N->CR3->co[i];

    if(rhs.N->change_dir_wind!=0)
    {
        N->change_dir_wind=new_longvector(rhs.N->change_dir_wind->nh);
        for(int i=1; i<=rhs.N->change_dir_wind->nh; i++)
            N->change_dir_wind->co[i]=rhs.N->change_dir_wind->co[i];
    }

    ///DOUBLEMATRIX *Psnow;/* The total snow fall to the ground,Dripp + through */
    N->Psnow=new_doublematrix(rhs.N->Psnow->nrh,rhs.N->Psnow->nch);
    for(int i=1; i<=rhs.N->Psnow->nrh; i++)
    {
        for(int j=1; j<=rhs.N->Psnow->nch; j++)
        {
            N->Psnow->co[i][j]=rhs.N->Psnow->co[i][j];
        }
    }

////////////////////////////////////////Glacier


    if(rhs.G->lnum!=0)
    {
        G->lnum=new_longmatrix(rhs.G->lnum->nrh,rhs.G->lnum->nch);
        for(int i=1; i<=rhs.G->lnum->nrh; i++)
        {
            for(int j=0; j<=rhs.G->lnum->nch; j++)
            {
                G->lnum->co[i][j]=rhs.G->lnum->co[i][j];
            }
        }
    }


    if(rhs.G->Dzl!=0)
    {
        G->Dzl=new_doubletensor(rhs.G->Dzl->ndh, rhs.G->Dzl->nrh, rhs.G->Dzl->nch);/** contains the soil parameters from the __soil.txt, pa[SoilType][Property][Layer] */
        initialize_doubletensor(G->Dzl, 0.0);

        for(int i=1; i<=rhs.G->Dzl->ndh; i++)
        {
            for(int j=1; j<=rhs.G->Dzl->nrh; j++)
            {
                for(int k=1; k<=rhs.G->Dzl->nch; k++)
                {
                    G->Dzl->co[i][j][k]=rhs.G->Dzl->co[i][j][k];
                }
            }
        }
    }

    if(rhs.G->w_liq!=0)
    {
        G->w_liq=new_doubletensor(rhs.G->w_liq->ndh, rhs.G->w_liq->nrh, rhs.G->w_liq->nch);
        initialize_doubletensor(G->w_liq, 0.0);

        for(int i=1; i<=rhs.G->w_liq->ndh; i++)
        {
            for(int j=1; j<=rhs.G->w_liq->nrh; j++)
            {
                for(int k=1; k<=rhs.G->w_liq->nch; k++)
                {
                    G->w_liq->co[i][j][k]=rhs.G->w_liq->co[i][j][k];
                }
            }
        }
    }

    if(rhs.G->w_ice!=0)
    {
        G->w_ice=new_doubletensor(rhs.G->w_ice->ndh, rhs.G->w_ice->nrh, rhs.G->w_ice->nch);
        initialize_doubletensor(G->w_ice, 0.0);

        for(int i=1; i<=rhs.G->w_ice->ndh; i++)
        {
            for(int j=1; j<=rhs.G->w_ice->nrh; j++)
            {
                for(int k=1; k<=rhs.G->w_ice->nch; k++)
                {
                    G->w_ice->co[i][j][k]=rhs.G->w_ice->co[i][j][k];
                }
            }
        }
    }

    if(rhs.G->T!=0)
    {
        G->T=new_doubletensor(rhs.G->T->ndh, rhs.G->T->nrh, rhs.G->T->nch);
        initialize_doubletensor(G->T, 0.0);

        for(int i=1; i<=rhs.G->T->ndh; i++)
        {
            for(int j=1; j<=rhs.G->T->nrh; j++)
            {
                for(int k=1; k<=rhs.G->T->nch; k++)
                {
                    G->T->co[i][j][k]=rhs.G->T->co[i][j][k];
                }
            }
        }
    }


    ///DOUBLEMATRIX *MELTED;
    if(rhs.G->MELTED!=0)
    {
        G->MELTED=new_doublematrix(rhs.G->MELTED->nrh,rhs.G->MELTED->nch);
        for(int i=1; i<=rhs.G->MELTED->nrh; i++)
        {
            for(int j=1; j<=rhs.G->MELTED->nch; j++)
            {
                G->MELTED->co[i][j]=rhs.G->MELTED->co[i][j];
            }
        }
    }

    ///DOUBLEMATRIX *SUBL;
    if(rhs.G->SUBL!=0)
    {
        G->SUBL=new_doublematrix(rhs.G->SUBL->nrh,rhs.G->SUBL->nch);
        for(int i=1; i<=rhs.G->SUBL->nrh; i++)
        {
            for(int j=1; j<=rhs.G->SUBL->nch; j++)
            {
                G->SUBL->co[i][j]=rhs.G->SUBL->co[i][j];
            }
        }
    }
////////////////////////////////////////METEO

    M->st=rhs.M->st;

    M->LRp=new_longvector(rhs.M->LRp->nh);
    for(int i=1; i<=rhs.M->LRp->nh; i++)
        M->LRp->co[i]=rhs.M->LRp->co[i];


    M->LRs=rhs.M->LRs;
    M->LRv=alloc1(5);
    for(int i=1; i<=3; i++)
    {
        M->LRv[M->LRp->co[i]]=NoV;
    }


    int nmetcols=0;
    for (int i=0; i<dim1l(rhs.M->column[0]); i++)
    {
        if (rhs.M->column[0][i]!=-1) nmetcols+=1;
    }

    M->column=alloc_long2(rhs.M->st->Z->nh);
    M->data=(double ***)malloc(rhs.M->st->Z->nh*sizeof(double**));
    M->horizon=(double ***)malloc(rhs.M->st->Z->nh*sizeof(double**));
    M->var=alloc2(rhs.M->st->Z->nh,nmet);

    int iii;
    for(int i=1; i<=rhs.M->st->Z->nh; i++)
    {

        M->column[i-1]=alloc_long1(nmet);
        for (int j=0; j<=nmet; j++)
        {
            M->column[i-1][j]=rhs.M->column[i-1][j];
        }

        iii=0;
        do
        {
            if(iii==0)
            {
                M->data[i-1]=(double **)malloc(sizeof(double*));
            }
            else
            {
                M->data[i-1]=(double **)realloc(M->data[i-1],(iii+1)*sizeof(double*));
            }
            M->data[i-1][iii]=(double *)malloc((nmetcols+1)*sizeof(double));

            for(int k=0; k<=nmetcols; k++)
            {
                M->data[i-1][iii][k]=rhs.M->data[i-1][iii][k];
            }
            iii++;
        }
        while(rhs.M->data[i-1][iii][0]!=end_vector);


        M->data[i-1]=(double **)realloc(M->data[i-1],(iii+1)*sizeof(double*));
        M->data[i-1][iii]=(double *)malloc(sizeof(double));
        M->data[i-1][iii][0]=end_vector;

        M->horizon[i-1]=alloc2(4,2);
        for(int j=1; j<=4; j++)
        {
            M->horizon[i-1][j-1][0]=rhs.M->horizon[i-1][j-1][0];
            M->horizon[i-1][j-1][1]=rhs.M->horizon[i-1][j-1][0];
        }
    }









    ///DOUBLEMATRIX *Tgrid;/**模型运行到某一时刻的驱动数据的空间分布，可有micromet插值得到，也可直接输入模型*/
    M->Tgrid=new_doublematrix(rhs.M->Tgrid->nrh,rhs.M->Tgrid->nch);
    for(int i=1; i<=rhs.M->Tgrid->nrh; i++)
    {
        for(int j=1; j<=rhs.M->Tgrid->nch; j++)
        {
            M->Tgrid->co[i][j]=rhs.M->Tgrid->co[i][j];
        }
    }


    M->Pgrid=new_doublematrix(rhs.M->Pgrid->nrh,rhs.M->Pgrid->nch);
    for(int i=1; i<=rhs.M->Pgrid->nrh; i++)
    {
        for(int j=1; j<=rhs.M->Pgrid->nch; j++)
        {
            M->Pgrid->co[i][j]=rhs.M->Pgrid->co[i][j];
        }
    }

    //DOUBLEMATRIX *Vgrid;
    M->Vgrid=new_doublematrix(rhs.M->Vgrid->nrh,rhs.M->Vgrid->nch);
    for(int i=1; i<=rhs.M->Vgrid->nrh; i++)
    {
        for(int j=1; j<=rhs.M->Vgrid->nch; j++)
        {
            M->Vgrid->co[i][j]=rhs.M->Vgrid->co[i][j];
        }
    }

    //DOUBLEMATRIX *Vdir;
    M->Vdir=new_doublematrix(rhs.M->Vdir->nrh,rhs.M->Vdir->nch);
    for(int i=1; i<=rhs.M->Vdir->nrh; i++)
    {
        for(int j=1; j<=rhs.M->Vdir->nch; j++)
        {
            M->Vdir->co[i][j]=rhs.M->Vdir->co[i][j];
        }
    }

    ///DOUBLEMATRIX *RHgrid;
    M->RHgrid=new_doublematrix(rhs.M->RHgrid->nrh,rhs.M->RHgrid->nch);
    for(int i=1; i<=rhs.M->RHgrid->nrh; i++)
    {
        for(int j=1; j<=rhs.M->RHgrid->nch; j++)
        {
            M->RHgrid->co[i][j]=rhs.M->RHgrid->co[i][j];
        }
    }

    ///DOUBLEMATRIX *Vspdmean;
    if(rhs.M->Vspdmean!=0)
    {
        M->Vspdmean=new_doublematrix(rhs.M->Vspdmean->nrh,rhs.M->Vspdmean->nch);
        for(int i=1; i<=rhs.M->Vspdmean->nrh; i++)
        {
            for(int j=1; j<=rhs.M->Vspdmean->nch; j++)
            {
                M->Vspdmean->co[i][j]=rhs.M->Vspdmean->co[i][j];
            }
        }
    }

    ///DOUBLEMATRIX *Vdirmean;
    if(rhs.M->Vdirmean!=0)
    {
        M->Vdirmean=new_doublematrix(rhs.M->Vdirmean->nrh,rhs.M->Vdirmean->nch);
        for(int i=1; i<=rhs.M->Vdirmean->nrh; i++)
        {
            for(int j=1; j<=rhs.M->Vdirmean->nch; j++)
            {
                M->Vdirmean->co[i][j]=rhs.M->Vdirmean->co[i][j];
            }
        }
    }

    ///DOUBLEMATRIX *RHmean;
    if(rhs.M->RHmean!=0)
    {
        M->RHmean=new_doublematrix(rhs.M->RHmean->nrh,rhs.M->RHmean->nch);
        for(int i=1; i<=rhs.M->RHmean->nrh; i++)
        {
            for(int j=1; j<=rhs.M->RHmean->nch; j++)
            {
                M->RHmean->co[i][j]=rhs.M->RHmean->co[i][j];
            }
        }
    }

    ///DOUBLEMATRIX *Taplot;
    if(rhs.M->Taplot!=0)
    {
        M->Taplot=new_doublematrix(rhs.M->Taplot->nrh,rhs.M->Taplot->nch);
        for(int i=1; i<=rhs.M->Taplot->nrh; i++)
        {
            for(int j=1; j<=rhs.M->Taplot->nch; j++)
            {
                M->Taplot->co[i][j]=rhs.M->Taplot->co[i][j];
            }
        }
    }

    ///DOUBLEMATRIX *Vspdplot;
    if(rhs.M->Vspdplot!=0)
    {
        M->Vspdplot=new_doublematrix(rhs.M->Vspdplot->nrh,rhs.M->Vspdplot->nch);
        for(int i=1; i<=rhs.M->Vspdplot->nrh; i++)
        {
            for(int j=1; j<=rhs.M->Vspdplot->nch; j++)
            {
                M->Vspdplot->co[i][j]=rhs.M->Vspdplot->co[i][j];
            }
        }
    }

    ///DOUBLEMATRIX *Vdirplot;
    if(rhs.M->Vdirplot!=0)
    {
        M->Vdirplot=new_doublematrix(rhs.M->Vdirplot->nrh,rhs.M->Vdirplot->nch);
        for(int i=1; i<=rhs.M->Vdirplot->nrh; i++)
        {
            for(int j=1; j<=rhs.M->Vdirplot->nch; j++)
            {
                M->Vdirplot->co[i][j]=rhs.M->Vdirplot->co[i][j];
            }
        }
    }

    ///DOUBLEMATRIX *RHplot;
    if(rhs.M->RHplot!=0)
    {
        M->RHplot=new_doublematrix(rhs.M->RHplot->nrh,rhs.M->RHplot->nch);
        for(int i=1; i<=rhs.M->RHplot->nrh; i++)
        {
            for(int j=1; j<=rhs.M->RHplot->nch; j++)
            {
                M->RHplot->co[i][j]=rhs.M->RHplot->co[i][j];
            }
        }
    }

    M->V=rhs.M->V;
    M->RH=rhs.M->RH;
    ///DOUBLEMATRIX *Tday;
    if(rhs.M->Tday!=0)
    {
        M->Tday=new_doublematrix(rhs.M->Tday->nrh,rhs.M->Tday->nch);
        for(int i=1; i<=rhs.M->Tday->nrh; i++)
        {
            for(int j=1; j<=rhs.M->Tday->nch; j++)
            {
                M->Tday->co[i][j]=rhs.M->Tday->co[i][j];
            }
        }
    }

    ///DOUBLEMATRIX *Tvar;
    if(rhs.M->Tvar!=0)
    {
        M->Tvar=new_doublematrix(rhs.M->Tvar->nrh,rhs.M->Tvar->nch);
        for(int i=1; i<=rhs.M->Tvar->nrh; i++)
        {
            for(int j=1; j<=rhs.M->Tvar->nch; j++)
            {
                M->Tvar->co[i][j]=rhs.M->Tvar->co[i][j];
            }
        }
    }

    M->nstsrad=rhs.M->nstsrad;/**nstsrad找到的具有段波辐射的第一个气象占的序号 1~气象站数目*/
    M->nstlrad=rhs.M->nstlrad;/**nstsrad找到的具有长波辐射的第一个气象占的序号 1~气象站数目*/
    M->nstcloud=rhs.M->nstcloud;/**nstsrad找到的具有云数据的第一个气象占的序号 1~气象站数目*/


    std::cout<<"Copy Constructor2"<<std::endl;
    return *this;
}

