#include "landgrid.h"
using namespace ldas;
using namespace std;
LandPoint& LandPoint::operator=(const LandPoint& other)
{
    if (this == &other) return *this; // handle self assignment
    //assignment operator
    lattitude=other.lattitude;
    longitude=other.longitude;
    return *this;
}
bool LandPoint::operator==(const LandPoint& other)
{
    return ((lattitude==other.lattitude) && (longitude==other.longitude));
}
LandPoint::LandPoint(const LandPoint& other)
{
    lattitude=other.lattitude;
    longitude=other.longitude;
}
LandGrid::LandGrid():m_use_mask(false),total_patches(0)
{
}

LandGrid::~LandGrid()
{
    if (m_use_mask)
    {
        delete []m_mask;
        delete []m_patch;
    }
}

LandGrid::LandGrid(const LandPoint& p,const unsigned int row,const unsigned int col,const double lat_span,const double lon_span):total_patches(0)
{
    corner_point=p;
    row_num=row;
    col_num=col;
    lattitude_span=lat_span;
    longitude_span=lon_span;
}

LandGrid::LandGrid(const LandGrid& other)
{
    m_use_mask=other.use_mask();
    longitude_span=other.lon_span();
    lattitude_span=other.lat_span();
    row_num=other.rows();
    col_num=other.cols();
    corner_point=other.corner();
    total_patches=other.patches();
    if (m_use_mask)
    {
        //copy the m_mask array
        m_mask=new bool[row_num*col_num];
        int k=0;
        for(int i=0; i<row_num*col_num; i++)
        {
            m_mask[i]=other.mask(i);
            if (m_mask[i])
            {
                m_patch[k]=i;
                ++i;
            }
        }
        total_patches=k;
    }
}

LandGrid& LandGrid::operator=(const LandGrid& other)
{
    if (this == &other) return *this; // handle self assignment
    //assignment operator
    m_use_mask=other.use_mask();
    longitude_span=other.lon_span();
    lattitude_span=other.lat_span();
    row_num=other.rows();
    col_num=other.cols();
    corner_point=other.corner();
    total_patches=other.patches();
    if (m_use_mask)
    {
        //copy the m_mask array
        m_mask=new bool[row_num*col_num];
        int k=0;
        for(int i=0; i<row_num*col_num; i++)
        {
            m_mask[i]=other.mask(i);
            if (m_mask[i])
            {
                m_patch[k]=i;
                ++k;
            }
        }
        total_patches=k;
    }
    return *this;
}
void LandGrid::mask(const char* fn)
{
    std::ifstream in;
    in.open(fn,std::ios_base::binary);
    if (!in)
    {
        string error_msg="could not open file:";
        error_msg+=fn;
        throw Exception("LandGrid","mask(fn)",error_msg.c_str());
    }
    char value;
    m_mask=new bool[row_num*col_num];
    m_patch=new int[row_num*col_num];
    m_use_mask=true;
    int k=0;
    for(int i=0; i<row_num; i++)
        for(int j=0; j<col_num; j++)
        {
            in.read((char *) &value,sizeof(value));
            if (value==48)//ascii, 0
                m_mask[i*col_num+j]=false;
            else
            {
                m_mask[i*col_num+j]=true;
                m_patch[k]=i*col_num+j;
                ++k;
            }
        }
    total_patches=k;
    in.close();
    cout<<"mask patches:"<<total_patches<<endl;
}
int LandGrid::location(const int index) const
{
    if (!m_use_mask) return index;
    return m_patch[index];
}
LandPoint LandGrid::index(const unsigned int i) const
{
    LandPoint p;
    p.lattitude=corner_point.lattitude-lattitude_span*int(i/col_num);
    p.longitude=corner_point.longitude+longitude_span*(i-int(i/col_num)*col_num);
    if (p.longitude>180)
        p.longitude-=360;
    if (p.longitude>180 || p.longitude<-180 || p.lattitude>90 || p.lattitude <-90)
        throw Exception("LandGrid","index()","Invalid position!");
    return p;
}
LandPoint LandGrid::index(const unsigned int row,const unsigned int col) const
{
    return index(row*col_num+col);
}
int LandGrid::row(const int index) const
{
	return int(location(index)/col_num);
}
int LandGrid::column(const int index) const
{
	return location(index)%col_num;
}
