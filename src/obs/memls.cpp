#include "memls.h"
using namespace ldas;
using namespace std;
MEMLS::MEMLS()
{
    //ctor
}

MEMLS::~MEMLS()
{
    //dtor
}

MEMLS::MEMLS(const MEMLS& other)
{
    //copy ctor
    m_snow=other.parameter();
}

MEMLS::MEMLS(const SnowParameter& sp):m_snow(sp)
{
    num.setDim(1,m_snow.col);
    //Ti=m_snow.temperature;
    //Wi=m_snow.wetness;
    //roi=m_snow.density;
    //di=m_snow.depth;
    //pci=m_snow.pci;
    for(int i=1; i<=m_snow.col; i++)
    {
        num.set(1,i,m_snow.col);
    }

    //pci.set(1,1,0.15/1000);
}
MEMLS& MEMLS::operator=(const MEMLS& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    m_snow=rhs.parameter();
    return *this;
}
int MEMLS::epsice ()
{
    try
    {
        Matrix <double> pp = ( 300.0 / m_snow.temperature )-1.0;
        const double B1 = 0.0207;
        const double b = 335.25;
        const double B2 = 1.16e-11;
        Matrix <double> db = exp(-10.02+0.0364 * (m_snow.temperature-273.0));
        Matrix <double> beta1 = (B1* exp(b/m_snow.temperature));
        Matrix <double> beta2 = dotP(m_snow.temperature,(exp(b/m_snow.temperature)-1.0)^2);
        Matrix <double>	beta = dotD(beta1,beta2) + B2 * (pow(m_snow.frequency,2.0)) + db;
        Matrix <double> alpha =dotP((0.00504 + 0.0062*pp),(exp(-22.1 * pp)));
        eice = (alpha/m_snow.frequency) + (beta*m_snow.frequency);
    }
    catch(ldas::Exception e)
    {
        e.what ();
        return 1;
    }
    return 0;
}

void MEMLS::epsr()
{
    epsi=m_snow.density;
    epsi.genZeros();
    Matrix <double> vfi = m_snow.density/0.917;
    double ehb = 0.99913;
    double esb = 1.4759;
    int i;
    for(i=1; i<=(int)m_snow.density.column(); i++)
    {
        if(m_snow.density.get(1,i)<=0.4)
            epsi.set(1,i,(1.0 + 1.5995 * m_snow.density.get(1,i) + 1.861 * pow(m_snow.density.get(1,i),3)));
        else //if ( epsi.get(1,i)>0.4)
            epsi.set(1,i,pow(((1.0 - vfi.get(1,i)) * ehb + vfi.get(1,i) * esb),3));
    }
}

void MEMLS::fresnelrc()
{
    int N=epsi.column();
    Matrix <double> FHFV_temp(1,N);
    FHFV_temp.genZeros();
    FH = FHFV_temp;
    FV = FHFV_temp;
    Matrix <double> epsr(1,N+1);
    for(int i=1; i<=N; i++)
    {
        epsr.set(1,i,epsi.get(1,i));
    }
    epsr.set(1,N+1,1.0);
    double epsn;
    double tein;
    double sinq;
    double qeps;
    double wurz;
    double wsub;
    double nd;
    for(int i=1; i<=N; i++)
    {
        epsn = epsr.get(1,i)/epsr.get(1,i+1);
        tein = tei.get(1,i+1);
        sinq = pow(std::sin(tein),2);
        qeps = sinq/epsn;
        wurz = std::sqrt(1.0-qeps);
        wsub = epsn-sinq;
        nd = std::sqrt(epsn);
        FH.set(1,i,((nd*wurz-std::cos(tein))/(nd*wurz+std::cos(tein))));
        FV.set(1,i,((wurz-nd*std::cos(tein))/(wurz+nd*std::cos(tein))));
    }
}

void MEMLS::ro2epsd ()
{
    epsice();
    epsr();
    Matrix <double> f = m_snow.density/0.917;
    const double ei = 3.185;
    long N = m_snow.density.column();
    Matrix <double> A(1,N);
    A += 0.3;
    int i;
    for(i=1; i<=N; i++)
    {
        if (f.get(1,i)<0.55)
            A.set(1,i,(0.476 - 0.64 * f.get(1,i)));
        if (f.get(1,i)<= 0.333)
            A.set(1,i,(0.1 + 0.5 * f.get(1,i)));
    }
    Matrix <double> epsp = epsi;
    Matrix <double> A3 = 1.0 - 2.0 * A;
    Matrix <double> ea = dotP(epsp,(1.0 - A)) + A;
    Matrix <double> ea3 = dotP(epsp,(-A3 + 1.0)) +A3;
    Matrix <double> K1 = (dotD(ea,ea + A * (ei-1.0)))^2;
    Matrix <double> K3 = (dotD(ea3,ea3 + A3 * (ei-1.0)))^2;
    Matrix <double> Ksq = (2.0 * K1 + K3) / 3.0;
    epsii = dotP(sqrt(epsi),dotP(eice,dotP(Ksq,f)));
}

void MEMLS::mixmod ()
{
    const double Aa = 0.005;
    const double Ab = 0.4975;
    const double Ac = 0.4975;
    const double euw = 4.9;
    const double esw = 88.045;
    const double frw = 0.11109;

    Matrix <double> esa = dotD((esw - epsi),3.0*(1.0+Aa*(esw/epsi-1.0)));
    Matrix <double> esb = dotD((esw - epsi),3.0*(1.0+Ab*(esw/epsi-1.0)));
    Matrix <double> esc = dotD((esw - epsi),3.0*(1.0+Ac*(esw/epsi-1.0)));
    Matrix <double> eua = dotD((euw - epsi),3.0*(1.0+Aa*(euw/epsi-1.0)));
    Matrix <double> eub = dotD((euw - epsi),3.0*(1.0+Ab*(euw/epsi-1.0)));
    Matrix <double> euc = dotD((euw - epsi),3.0*(1.0+Ac*(euw/epsi-1.0)));

    Matrix <double> fa = 1.0 + Aa * (esw - euw) / (epsi + Aa * (euw - epsi));
    Matrix <double> fb = 1.0 + Ab * (esw - euw) / (epsi + Ab * (euw - epsi));
    Matrix <double> fc = 1.0 + Ac * (esw - euw) / (epsi + Ac * (euw - epsi));

    Matrix <double> eea = esa -eua;
    Matrix <double> eeb = esb -eub;
    Matrix <double> eec = esc -euc;

    Matrix <double> fwa = frw / fa;
    Matrix <double> fwb = frw / fb;
    Matrix <double> fwc = frw / fc;

    Matrix <double> depsia = eua + dotD(eea,(1.0 + ((fwa * m_snow.frequency)^2)));
    Matrix <double> depsib = eub + dotD(eeb,1.0 + ((fwb * m_snow.frequency)^2));
    Matrix <double> depsic = euc + dotD(eec,1.0 + ((fwc * m_snow.frequency)^2));
    Matrix <double> depsi = dotP(m_snow.wetness,depsia + depsib + depsic);

    Matrix <double> depsiia = dotP((fwa * m_snow.frequency),dotD(eea,1.0 + ((fwa * m_snow.frequency)^2)));
    Matrix <double> depsiib = dotP((fwb * m_snow.frequency),dotD(eeb,1.0 + ((fwb * m_snow.frequency)^2)));
    Matrix <double> depsiic = dotP((fwc * m_snow.frequency),dotD(eec,1.0 + ((fwc * m_snow.frequency)^2)));
    Matrix <double> depsii = dotP(m_snow.wetness,depsiia + depsiib + depsiic);

    epsi = epsi + depsi;
    epsii = epsii + depsii;
}

void MEMLS::abscoeff ()
{
    const double c = 2.99793;
    ComplexMatrix <double> epsi_ii(1,m_snow.col,epsi,epsii);
    gai = (4*PI_NUMBER*10*m_snow.frequency)*((sqrt(epsi_ii)).imag())/c;
    //注意涉及到了虚数运算
}

void MEMLS::pfadi ()
{
    long N = m_snow.depth.column();
    dei = m_snow.depth;
    //由于tei为N+1列，所以不能使用dei=di.DotD(Cos(tei));
    for(int i=1; i<=N; i++)
    {
        dei.set(1,i, m_snow.depth.get(1,i)/std::cos(tei.get(1,i)));
    }
}
void MEMLS::fresnelc ()
{
    long N = epsi.column();

    siv = epsi;
    siv.genZeros();
    sih = epsi;
    sih.genZeros();
    Matrix <double> epsi_temp = tei;
    for(int i=1; i<=N; i++)
    {
        epsi_temp.set(1,i,epsi.get(1,i));
    }
    epsi_temp.set(1,N+1,1.0);

    for(int n=1; n<=N; n++)
    {
        double epso = epsi_temp.get(1,n+1);
        double epsu = epsi_temp.get(1,n);
        double tein = tei.get(1,n+1);
        sih.set(1,n,pow(((std::sqrt(epso)*std::cos(tein) - std::sqrt(epsu - epso * pow(std::sin(tein),2)))/(std::sqrt(epso)*std::cos(tein) + std::sqrt(epsu - epso * pow(std::sin(tein),2)))),2));
        siv.set(1,n,pow(((epsu*std::cos(tein) - std::sqrt(epso)*std::sqrt(epsu - epso*pow(std::sin(tein),2)))/(epsu*std::cos(tein) + std::sqrt(epso)*std::sqrt(epsu - epso*pow(std::sin(tein),2)))),2));
    }
}

void MEMLS::slred ()
{
    const double cc = 0.299793;
    const double FIC = 4 * PI_NUMBER * m_snow.frequency / cc;
    double fc = 4.712;
    const double repo = 0;
    int N = m_snow.density.column();
    double theta = tei.get(1,N+1);
    Matrix <double> ns = sqrt(epsi);
    Matrix <double> tei_temp(1,N);
    int m;
    for(m=1; m<=N; m++)
        tei_temp.set(1,m,tei.get(1,m));
    Matrix <double> fi = FIC*(dotP(m_snow.depth,dotP(ns,cos(tei_temp))));

    //if (repo==1)
    //	cout<<"repo==1"<<num<<di<<roi<<Ti<<pci<<ns<<gai<<epsi<<epsii<<fi<<tei*180/pi<<endl;

    Matrix <int> A(1,N);
    A.genZeros();
    int M_ColNum=0;
    for (m=1; m<=(int)fi.column(); m++)
    {
        if(fi.get(1,m)<fc)
        {
            M_ColNum++;
            A.set(1,m,1);
        }
    }
    A.set(1,1,1);
    if (M_ColNum > 0)
    {
        //cout<<M_ColNum<<" coherent layers of "<<N <<" detected: "<<freq<<" GHz"<<endl;

        int pl = 0;
        int sc = 0;
        int scmax = 0;
        int ml = 0;
        int mlo = 0;

        for(m=2; m<=N; m++)
        {
            if (A.get(1,m)==1 && pl==1)
            {
                if(ml==0)
                {
                    sc = sc +1;
                    ml = 1;
                    A.set(1,m-1,sc);
                }
                A.set(1,m,sc);
                scmax = sc;
            }
            else
            {
                if(pl==1)
                    pl = 0;
                else
                {
                    if(A.get(1,m) == 1)
                    {
                        pl = 1;
                        ml = 0;
                    }
                }
            }
        }
        //  combine succeeding coherent layers by weighting with the phase
        if (scmax>0)
        {
            for(m=2; m<=scmax; m++)
            {
                int B_ColNum = 0;
                int tal = 0;
                double fitot = 0;
                int j;
                for(j=1; j<N; j++)
                {
                    if(A.get(1,j) == m)
                        fitot += fi.get(1,j);
                    tal = j;
                    B_ColNum++;
                }
                Matrix <double> fitv(1,B_ColNum);
                int fitv_ = 1;
                for(j=1; j<N; j++)
                {
                    if(A.get(1,j) == m)
                    {
                        fitv.set(1,fitv_,fi.get(1,j)/fitot);
                        fitv_ ++;
                    }
                }
                double di_tal=0;
                double dei_tal=0;
                double roi_tal=0;
                double Ti_tal=0;
                double Wi_tal=0;
                double pci_tal=0;
                double ns_tal=0;
                double gai_tal=0;
                double epsi_tal=0;
                double epsii_tal=0;
                double tei_tal=0;
                for(j=1; j<N; j++)
                {
                    if(A.get(1,j) == m)
                    {
                        di_tal+=m_snow.depth.get(1,j);
                        dei_tal += dei.get(1,j);
                        roi_tal += (m_snow.density.get(1,j)*fitv.get(1,j));
                        Ti_tal += (m_snow.temperature.get(1,j)*fitv.get(1,j));
                        Wi_tal += (m_snow.wetness.get(1,j)*fitv.get(1,j));
                        pci_tal += (m_snow.pci.get(1,j)*fitv.get(1,j));
                        ns_tal += (ns.get(1,j)*fitv.get(1,j));
                        gai_tal += (gai.get(1,j)*fitv.get(1,j));
                        epsi_tal += (epsi.get(1,j)*fitv.get(1,j));
                        epsii_tal += (epsii.get(1,j)*fitv.get(1,j));
                        tei_tal += (tei.get(1,j)*fitv.get(1,j));
                    }
                }
                m_snow.depth.set(1,tal,di_tal);
                dei.set(1,tal,dei_tal);
                m_snow.density.set(1,tal,roi_tal);
                m_snow.temperature.set(1,tal,Ti_tal);
                m_snow.wetness.set(1,tal,Wi_tal);
                m_snow.pci.set(1,tal,pci_tal);
                ns.set(1,tal,ns_tal);
                gai.set(1,tal,gai_tal);
                epsi.set(1,tal,epsi_tal);
                epsii.set(1,tal,epsii_tal);
                tei.set(1,tal,tei_tal);
                A.set(1,tal,1);
            }
        }

        int A_i=0;
        int j;
        for(j=1; j<=(int)A.column(); j++)
        {
            if(A.get(1,j)<2)
                A_i++;
        }
        if(A_i>0)
        {
            int k=1;
            for(j=1; j<=(int)num.column(); j++)
            {
                if(A.get(1,j)<2)
                {
                    num.set(1,k,num.get(1,j));
                    m_snow.density.set(1,k,m_snow.density.get(1,j));
                    tei.set(1,k,tei.get(1,j));
                    m_snow.depth.set(1,k,m_snow.depth.get(1,j));
                    m_snow.temperature.set(1,k,m_snow.temperature.get(1,j));
                    m_snow.pci.set(1,k,m_snow.pci.get(1,j));
                    m_snow.wetness.set(1,k,m_snow.wetness.get(1,j));
                    gai.set(1,k,gai.get(1,j));
                    ns.set(1,k,ns.get(1,j));
                    epsi.set(1,k,epsi.get(1,j));
                    epsii.set(1,k,epsii.get(1,j));
                    fi.set(1,k,fi.get(1,j));
                    dei.set(1,k,dei.get(1,j));
                }
                k++;
            }
        }

        N = m_snow.density.column();
        A = num;
        A.genZeros();
        for(j=1; j<=N; j++)
        {
            if(fi.get(1,j)<fc)
                A.set(1,j,1);
        }
        A.set(1,1,0);
        sih = m_snow.density;
        siv = m_snow.density;
        Matrix <double> X = m_snow.density;
        sih.genZeros();
        siv.genZeros();
        X.genZeros();

        //	  calculate interface reflection coefficients
        fresnelrc();

        for(j=1; j<=N; j++)
        {
            if(A.get(1,j)==0)
            {
                sih.set(1,j,pow(FH.get(1,j),2));
                siv.set(1,j,pow(FV.get(1,j),2));
            }
        }
        ////for layers of type 1 shi-1 = ...

        for(j=1; j<=N; j++)
        {
            if(A.get(1,j)==1)
            {
                X.set(1,j,2*FH.get(1,j)*FH.get(1,j-1)*std::cos(fi.get(1,j)));
                sih.set(1,j-1,(pow(FH.get(1,j),2)+pow(FH.get(1,j-1),2)+X.get(1,j))
                        /(1+pow(FH.get(1,j-1),2)*pow(FH.get(1,j-1),2)+X.get(1,j)));
                X.set(1,j,2*FV.get(1,j)*FV.get(1,j-1)*std::cos(fi.get(1,j)));
                siv.set(1,j-1,(pow(FV.get(1,j),2)+pow(FV.get(1,j-1),2)+X.get(1,j))
                        /(1+pow(FV.get(1,j-1),2)*pow(FV.get(1,j-1),2)+X.get(1,j)));
            }
        }
        N = m_snow.depth.column();

        //	  remove layers of type 1
        int col_new = 0;
        for(j=1; j<=N; j++)
        {
            if(A.get(1,j)==0)
            {
                col_new++;
            }
        }
        Matrix <double> ns_new(1,col_new);
        Matrix <int> num_new(1,col_new);
        Matrix <double> roi_new(1,col_new);
        Matrix <double> tei_new(1,col_new);
        Matrix <double> di_new(1,col_new);
        Matrix <double> Ti_new(1,col_new);
        Matrix <double> pci_new(1,col_new);
        Matrix <double> Wi_new(1,col_new);
        Matrix <double> gai_new(1,col_new);
        Matrix <double> FH_new(1,col_new);
        Matrix <double> FV_new(1,col_new);
        Matrix <double> sih_new(1,col_new);
        Matrix <double> siv_new(1,col_new);
        Matrix <double> fi_new(1,col_new);
        Matrix <double> epsi_new(1,col_new);
        Matrix <double> epsii_new(1,col_new);
        Matrix <double> dei_new(1,col_new);
        col_new=1;
        for(j=1; j<=N; j++)
        {
            if(A.get(1,j)==0)
            {
                ns_new.set(1,col_new,ns.get(1,j));
                num_new.set(1,col_new,num.get(1,j));
                roi_new.set(1,col_new,m_snow.density.get(1,j));
                tei_new.set(1,col_new,tei.get(1,j));
                di_new.set(1,col_new,m_snow.depth.get(1,j));
                Ti_new.set(1,col_new,m_snow.temperature.get(1,j));
                pci_new.set(1,col_new,m_snow.pci.get(1,j));
                Wi_new.set(1,col_new,m_snow.wetness.get(1,j));
                gai_new.set(1,col_new,gai.get(1,j));
                FH_new.set(1,col_new,FH.get(1,j));
                FV_new.set(1,col_new,FV.get(1,j));
                sih_new.set(1,col_new,sih.get(1,j));
                siv_new.set(1,col_new,siv.get(1,j));
                fi_new.set(1,col_new,fi.get(1,j));
                epsi_new.set(1,col_new,epsi.get(1,j));
                epsii_new.set(1,col_new,epsii.get(1,j));
                dei_new.set(1,col_new,dei.get(1,j));
                col_new++;
            }
        }
        N = di_new.column();
        Matrix <double> tei_temp(1,N+1);
        for(j=1; j<=N; j++)
        {
            tei_temp.set(1,j,tei_new.get(1,j));
        }
        tei_temp.set(1,N+1,theta);
        tei = tei_temp;
        ns=ns_new;
        num=num_new;
        m_snow.density=roi_new;
        m_snow.depth=di_new;
        m_snow.temperature=Ti_new;
        m_snow.pci=pci_new;
        m_snow.wetness=Wi_new;
        gai=gai_new;
        FH=FH_new;
        FV=FV_new;
        sih=sih_new;
        siv=siv_new;
        fi=fi_new;
        epsi=epsi_new;
        epsii=epsii_new;
        dei=dei_new;
    }
    //else
    //cout<<"no coherent layers detected: "<<freq<<"GHz"<<endl;
}

void MEMLS::sccoeff (
    Matrix <double>& gbih,
    Matrix <double>& gbiv,
    Matrix <double>& ga2i)
{
    //specular component of scattering coefficient
    //usually 0 can be important in new snow!
    double c = 2.99;
    double roair = 0.001293;
    double roice = 0.917;
    double dgb0h = 0.0;
    double dgb0v = 0.0;

    // aus der Theorie scattering coefficient
    double k = m_snow.frequency*(2*PI_NUMBER/0.299793);
    double eice_sccoeff = 3.18;
    Matrix <double> vfi = m_snow.density/roice;

    //choose the scattering algorithm that should be used
    const double sccho=11;
    double wahl = sccho;

    ro2epsd();
    mixmod();
    //6-flux scattering coefficient

    if (wahl == 1)
    {
        gs6 = dotD((130 * pow(m_snow.frequency/50,2.7))*(m_snow.pci^3),((m_snow.density^1.3)+ 0.001));
    }
    //	fit vom 26.8.97 auf alle Daten v-pol, > 11 GHz
    if (wahl ==2)
    {
        gs6 = dotP(0.0704 * (m_snow.pci^2.32)*pow(m_snow.frequency,1.68),(m_snow.density^(-0.63)));
    }
    //for spheres
    Matrix <double> epseff = (2.0-eice+3.0*vfi*(eice_sccoeff-1)+ sqrt((1.0+3.0*vfi*((eice_sccoeff-1))^2)+8*eice_sccoeff))/4.0;
    Matrix <double> sphe = (3/32)*pow(k,4)*dotP(dotP(dotP(((0.001*m_snow.pci)^3),vfi),(1.0-vfi)),(dotD((2.0*epseff+1.0)*(eice_sccoeff-1),(2.0*epseff+eice_sccoeff))^2));

    if (wahl==4)
    {
        gs6 = sphe;
    }
    //for shells
    epseff = 1.0+dotD((vfi*(eice_sccoeff-1)*(2+1/eice_sccoeff)),(3.0-vfi*(1-1/eice_sccoeff)));
    Matrix <double> shel = (2/3 + 1/(3*pow(eice_sccoeff,2)))*dotP(dotP(((0.001*m_snow.pci)^3)*pow(k,2),vfi),dotD((1.0-vfi)*pow((eice_sccoeff-1),2),16.0*epseff));
    if (wahl == 5)
    {
        gs6 = shel;
    }

    //as linearcombination
    if (wahl == 6)
    {
        double a = 0.1664;
        double b = 0.2545;
        gs6 = a *sphe+b*shel;
    }

    //fit vom 26.9.97
    if (wahl == 7)
    {
        gs6 = 73.21*pow((m_snow.frequency/50),2.68)*dotP((m_snow.pci^3),m_snow.density^(-1));
    }
    if (wahl == 8)
    {
        gs6 = 136 * (pow((m_snow.frequency/50),2.5))* dotD((m_snow.pci^2.85),(m_snow.density + 0.001));
    }

    //fit vom 4.11.97 (without density)
    if (wahl == 9)
    {
        gs6 = 564.0 * (m_snow.pci^3.0) * pow((m_snow.frequency/50),2.5);
    }

    //fit vom 4.11.97 (without density, uses corr. length from exp. fit!)
    if (wahl == 10)
    {
        gs6 = (3.16 * m_snow.pci + 295.0 * (m_snow.pci^2.5))* pow((m_snow.frequency/50),2.5);
    }

    //fit vom 4.11.97 (with density, uses corr. length from exp. fit!)
    if (wahl == 11)
    {
        gs6 = ((9.20 * m_snow.pci - 1.23 * m_snow.density + 0.54)^2.5) * pow((m_snow.frequency/50),2.5);
        //cout<<"gs6:  "<<gs6<<endl;
    }
    Matrix <double> omega = sqrt(dotD((epsi-1.0),epsi));

    //Born Approximation
    if (wahl == 12)
    {
        //Sorry!
        //Born approximation is not provided at present.
        //I will find time to do that.
        //By CHE Tao
    }
    else
    {
        gb6 = 0.5*dotP(gs6,1.0-omega);
        gc6 = 0.25*dotP(gs6,omega);
    }

    int N = gs6.column();
    Matrix <double> gtr(1,N);
    Matrix <double> ga2(1,N);
    gtr.genZeros();
    ga2.genZeros();
    //gs6.genZeros();
    // -> 2 Flux
    gtr = dotD((4.0*gc6),(gai+2.0*gc6));
    ga2i = dotP(gai,(1.0+gtr));
    gbih = (gb6 + dgb0h) + dotP(gtr,gc6);
    gbiv = (gb6 + dgb0v) + dotP(gtr,gc6);

}

void MEMLS::pfadc (const double& teta)
{
    //*****	here is some problem aobut the Matrix
    const unsigned long N = m_snow.depth.column();
    Matrix <double> ns = sqrt(epsi);
    Matrix <double> costetasn = sqrt(1.0-(((std::sin(teta))/ns)^2));
    Matrix <double> cosc = sqrt(1.0-((1.0/ns)^2));
    Matrix <double> costetasc = 0.5 * (1.0 + cosc);
    dei = dotD(m_snow.depth,costetasn);

    int tauscat_N = epsi.column()+1;
    Matrix <double> tauscat(1,tauscat_N);
    tauscat.genZeros();
    tscat = epsi;
    tscat.genZeros();
    Matrix <double> costeta = epsi;
    costeta.genZeros();

    for(int m=N; m>=1; m--)
    {
        tauscat.set(1,m,(tauscat.get(1,m+1)+dei.get(1,m)*gs6.get(1,m)/2.0));
        tscat.set(1,m,std::exp(-1*tauscat.get(1,m)));
        costeta.set(1,m,(tscat.get(1,m)*costetasn.get(1,m)+(1-tscat.get(1,m))*costetasc.get(1,m)));
    }
    costeta.acos();
    tei = costeta;
}

void MEMLS::polmix ()
{
    Matrix <double> tscat_temp(1,tscat.column()+1);
    for(int i=1; i<=(int)tscat.column(); i++)
    {
        tscat_temp.set(1,i,tscat.get(1,i));
    }
    tscat_temp.set(1,tscat_temp.column(),1.0);
    tscat = tscat_temp;
    Matrix <double> smean = 0.5 * (sih + siv);
    Matrix <double> deltas = 0.5 * dotP(tscat,(sih - siv));
    sih = smean + deltas;
    siv = smean - deltas;
}

void MEMLS::rt (                Matrix <double>& gbi                )
{
    Matrix <double> gamma = sqrt(dotP(ga2i,(ga2i + 2.0*gbi)));
    Matrix <double> t0i = exp(-1.0*dotP(gamma,dei));
    Matrix <double> r0i = t0i;
    r0i.genZeros();

    for(int i=1; i<=(int)gbi.column(); i++)
    {
        if(gbi.get(1,i)>0.00001)
        {
            r0i.set(1,i,gbi.get(1,i)/(ga2i.get(1,i)+gbi.get(1,i)+gamma.get(1,i)));
        }
    }
    Matrix <double> r02 = r0i^2;
    Matrix <double> t02 = t0i^2;
    ri = dotD(dotP(r0i,(1.0 - t02)),(1.0 - dotP(t02,r02)));
    ti = dotD(dotP(t0i,(1.0 - r02)),(1.0 - dotP(t02,r02)));
}

void MEMLS::layer (Matrix <double>& si)
{
    int N = ri.column();
    int i,j;
    Matrix <double> ei = 1.0 - ri - ti;

    if (N<1)
    {
        cout<<"ERROR: NO scattering layer"<<endl;
        return;
    }
    if (N==1)
    {
        Matrix <double> D_temp(1,1);
        double k1 = (ri.get(1,1)*(1.0-si.get(1,1))*m_snow.temperature_ground+ti.get(1,1)*(1.0-si.get(1,2))*m_snow.temperature_sky+ei.get(1,1)*m_snow.temperature.get(1,1))/(1.0-ri.get(1,1)*si.get(1,1));
        double k2 = ti.get(1,1)*si.get(1,1)/(1.0-ri.get(1,1)*si.get(1,1));
        double k3 = 1.0-ri.get(1,1)*si.get(1,2)-ti.get(1,1)*si.get(1,1)*k2;
        D_temp.set(1,1,(ti.get(1,1)*si.get(1,1)*k1+ti.get(1,1)*(1-si.get(1,1))*m_snow.temperature_ground+ri.get(1,1)*(1-si.get(1,2))*m_snow.temperature_sky+ei.get(1,1)*m_snow.temperature.get(1,1))/k3);
        D=D_temp;
    }
    else
    {

        Matrix <double> M1(N,N);
        M1.genZeros();
        Matrix <double> H(N-1,N-1);
        H.genZeros();
        for(i=1; i<=N; i++)
        {
            M1.set(i,i,ri.get(1,i)*si.get(1,i));
        }
        for(i=1; i<=N-1; i++)
        {
            H.set(i,i,ti.get(1,i)*(1.0-si.get(1,i+1)));
        }
        for(i=1; i<=N-1; i++)
        {
            for(j=2; j<=N; j++)
            {
                M1.set(i,j,M1.get(i,j)+H.get(i,j-1));
            }
        }

        Matrix <double> I(N,N);
        for(i=1; i<=N; i++)
        {
            I.set(i,i,1);
        }

        Matrix <double> M2(N,N);
        M2.genZeros();
        H.genZeros();
        for(i=1; i<=N; i++)
        {
            M2.set(i,i,ti.get(1,i)*si.get(1,i+1));
        }
        for(i=1; i<=N-1; i++)
        {
            H.set(i,i,ri.get(1,i+1)*(1.0-si.get(1,i+1)));
        }
        for(i=2; i<=N; i++)
        {
            for(j=1; j<=N-1; j++)
            {
                M2.set(i,j,M2.get(i,j)+H.get(i-1 ,j));
            }
        }

        Matrix <double> M3(N,N);
        M3.genZeros();
        H.genZeros();
        for(i=1; i<=N; i++)
        {
            M3.set(i,i,ti.get(1,i)*si.get(1,i));
        }
        for(i=1; i<=N-1; i++)
        {
            H.set(i,i,ri.get(1,i)*(1.0-si.get(1,i+1)));
        }
        for(i=1; i<=N-1; i++)
        {
            for(j=2; j<=N; j++)
            {
                M3.set(i,j,M3.get(i,j)+H.get(i,j-1));
            }
        }

        Matrix <double> M4(N,N)			  ;
        M4.genZeros();
        H.genZeros();
        for(i=1; i<=N; i++)
        {
            M4.set(i,i,ri.get(1,i)*si.get(1,i+1));
        }
        for(i=1; i<=N-1; i++)
        {
            H.set(i,i,ti.get(1,i+1)*(1.0-si.get(1,i+1)));
        }
        for(i=2; i<=N; i++)
        {
            for(j=1; j<=N-1; j++)
            {
                M4.set(i,j,M4.get(i,j)+H.get(i-1,j));
            }
        }

        Matrix <double> E;
        E = dotP(ei,m_snow.temperature);
        E.set(1,1,E.get(1,1)+ri.get(1,1)*(1.0-si.get(1,1))*m_snow.temperature_ground);
        E.set(1,N,E.get(1,N)+ti.get(1,N)*(1.0-si.get(1,N+1))*m_snow.temperature_sky);

        Matrix <double> F;
        F = dotP(ei,m_snow.temperature);
        F.set(1,1,F.get(1,1)+ti.get(1,1)*(1.0-si.get(1,1))*m_snow.temperature_ground);
        F.set(1,N,F.get(1,N)+ri.get(1,N)*(1.0-si.get(1,N+1))*m_snow.temperature_sky);

        Matrix <double> M5;
        M5 = M3 * ((I - M1).inv() * M2) + M4;

        Matrix <double> D_tempM;
        try
        {
            D = ((I - M5).inv()) * (((M3 * ((I - M1).inv())) * E.tran()) + F.tran());

        }
        catch(ldas::Exception e)
        {
            e.what();
            return;
        }
    }//else end,ok
}

void MEMLS::run()
{
    //int col = 7;
    //double freq=10.7;
    //double teta=10.0*pi/180;
    //double Tsky=0;
    //double Tgnd=273.0;
    double teta = m_snow.theta*PI_NUMBER/180.0;
    ro2epsd ();
    mixmod ();
    abscoeff ();
    Matrix <double> ns = sqrt(epsi);
    //初始化tei
    Matrix <double> tei_temp(1,m_snow.col+1);
    int i;
    for(i=1; i<=m_snow.col; i++)
    {
        tei_temp.set(1,i,std::asin(std::sin(teta)/ns.get(1,i)));
    }
    tei_temp.set(1,m_snow.col+1,teta);
    tei = tei_temp;
    //初始化tei结束

    pfadi ();
    fresnelc ();
    slred ();
    sccoeff (gbih,gbiv,ga2i);
    pfadc (teta);

    int N = sih.column();
    Matrix <double> sih_temp(1,N+1);
    Matrix <double> siv_temp(1,N+1);

    const double s0h=0;
    const double s0v=0;
    sih_temp.set(1,1,s0h);
    siv_temp.set(1,1,s0v);
    for(i=1; i<=N; i++)
    {
        sih_temp.set(1,i+1,sih.get(1,i));
        siv_temp.set(1,i+1,siv.get(1,i));
    }

    sih = sih_temp;
    siv = siv_temp;
    polmix ();
    rt (gbih);
    layer(sih);
    N=m_snow.density.column();

    //double TBh,TBv;
    m_tb.h=(1.0-sih.get(1,N+1))*D.get(N,1) + sih.get(1,N+1)*m_snow.temperature_sky;
    rt(gbiv);
    layer(siv);
    m_tb.v=(1.0-siv.get(1,N+1))*D.get(N,1) + siv.get(1,N+1)*m_snow.temperature_sky;
    //cout<<"absorption coefficient: "<<gai<<endl;
    //cout<<"scattering coefficient in h-p: "<<gbih<<endl;
    //cout<<"scattering coefficient in h-v: "<<gbiv<<endl;
    //cout<<"TBh = "<<TBh<<endl;
    //cout<<"TBv = "<<TBv<<endl;
}
