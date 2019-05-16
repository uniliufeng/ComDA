#ifndef _LDAS_LANDGRID_H
#define _LDAS_LANDGRID_H

#include "exception.h"
#include "constant.h"
#include <iostream>
#include <fstream>
namespace ldas
{
class LandPoint
{
public:
    LandPoint():lattitude(0),longitude(0)
    {
    }
    LandPoint(const double lat,const double lon):lattitude(lat),longitude(lon)
    {
    }
    LandPoint(const LandPoint& other);
    LandPoint& operator=(const LandPoint& other);
    bool operator==(const LandPoint& other);
    double lon_radiance() const
    {
        return longitude*PI_NUMBER/360;
    }
    double lat_radiance() const
    {
        return lattitude*PI_NUMBER/360;
    }
    double lattitude,longitude;
};
class LandGrid
{
public:
    /** Default constructor */
    LandGrid();
    /** Default destructor */
    virtual ~LandGrid();
    /** Copy constructor
     *  \param other Object to copy from
     */
    LandGrid(const LandGrid& other);
    LandGrid(const LandPoint& p,const unsigned int row,const unsigned int col,const double lat_span,const double lon_span);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    LandGrid& operator=(const LandGrid& other);
    LandPoint index(const unsigned int row,const unsigned int col) const;
    LandPoint index(const unsigned int i) const;
    ///return patch index's landpoint
    LandPoint operator[](const int i) const
    {
        return index(location(i));
    }
    virtual void mask(const char* fn);
    bool mask(const unsigned int i) const
    {
        if (!m_use_mask)
            throw Exception("LandGrid","mask(int)","mask undefined!");
        return m_mask[i];
    }
    int patches() const
    {
        return total_patches;
    }
    ///row and col are starting from 0
    bool mask(const unsigned int row,const unsigned int col) const
    {
        if (!m_use_mask)
            throw Exception("LandGrid","mask(int)","mask undefined!");
        return m_mask[row*col_num+col];
    }
    ///return the value in the grid with a z-shape.
    int location(const int index) const;
    int rows() const
    {
        return row_num;
    }
    void rows(const unsigned int row)
    {
        row_num=row;
    }
    int cols() const
    {
        return col_num;
    }
    void cols(const unsigned int col)
    {
        col_num=col;
    }
    LandPoint corner() const
    {
        return corner_point;
    }
    bool use_mask() const
    {
        return m_use_mask;
    }
    double lon_span() const
    {
        return longitude_span;
    }
    void lon_span(const double span)
    {
        longitude_span=span;
    }
    double lat_span() const
    {
        return lattitude_span;
    }
    void lat_span(const double span)
    {
        lattitude_span=span;
    }
    void corner(const LandPoint& p)
    {
        corner_point=p;
    }
    int row(const int index)const;
    int column(const int index)const;
protected:
private:
    bool* m_mask;
    int* m_patch;
    double longitude_span,lattitude_span;
    ///left,top corner point
    LandPoint corner_point;
    //double min_longitude,max_longitude,min_lattitude,max_lattitude;
    int total_patches;
    int row_num,col_num;
    bool m_use_mask;
};
}
#endif // _LDAS_LANDGRID_H
