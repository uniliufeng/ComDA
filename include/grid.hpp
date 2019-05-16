// grid.hpp
// <Li Xin> Feb 19, 2001
// <Li Xin> rewrite in Aug 31, 2004

#ifndef GRID_H
#define GRID_H


#include <cmath>
#include <iostream>
//#include <stdlib.h>
#include "space3d.hpp"

#ifndef PI
#define PI ( 4.0*atan(1.0) )
#endif


// enum eProjectionType {Gauss, UTM, Trans, Albert, Lambert, LatLon};
// enum eRefrenceType {Rec, LatLon, Image};
// enum eUnit {m, km, dd, dms, rad, pixel};



namespace ldas
{
  //template <class TYPE> class Grid;
  // template<> class Grid<unsigned char>;
  // template<> class Grid<short int>;
  // template<> class Grid<int>;
  // template<> class Grid<float>;
  // template<> class Grid<double>;


  // class definition
  //template <class TYPE>
  class Grid
  {
  protected:
    double**	data;				// data matrix
    double	noValueData;		// novalue data	
    int		rows;				// rows
    int		cols;				// columns
    double	resolution;			// resolution
    double	llX;				// lowleft X
    double	llY;				// lowleft Y
    double	urX;				// upright X
    double	urY;				// upright Y
		
    /* eDataType DataType;
       unsigned int BandCount;
       eStorageType StorageType;
       eProjectionType ProjectionType;
       eRefrenceType RefrenceType;
       eUnit RefrenceUnit;
       eUnit ValueUnit;
       long MinX, MaxX, MinY, MaxY; */
		
    // file types
    enum fileType
      {
        generalBinary			= 1,
        generalASCII			= 2,
        IDRISI					= 10,
        ArcinfoGrid				= 20,
      };

    // data types
    enum dataType
      {
        BYTE			= 1,
        SHORT			= 2,
        INT				= 3,
        FLOAT			= 4,
        DOUBLE			= 5,
      };
	
  public:
    // constructor
    Grid();
    Grid(const int m_rows, const int m_cols);
    Grid(const int m_rows, const int m_cols, const double m_resolution);
    Grid(const int m_rows, const int m_cols, const double m_resolution, 
         const double m_llX, const double m_llY);
    Grid(const int m_rows, const int m_cols, const double m_llX, const double m_llY,
         const double m_urX, const double m_urY);
    Grid(const Grid& g);
	
    // destructor
    ~Grid();
		
    // grid map attribute operators
    void SetDimension(const int m_rows, const int m_cols, const double m_resolution=0, 
                      const double m_llX=0, const double m_llY=0);
    void SetResolution(double m_resolution);
    void SetNoValueData(double m_noValueData);
    int GetRows() const;
    int GetCols() const;
    double GetResolution() const;
    double GetllX() const {return llX;}
    double GetllY() const {return llY;}

    // grid cell access
    double GetAt(int row, int col) const;
    void SetAt(int row, int col, double newElement);
		
    // get grid cell coordinates
    // point getPixelPos(int row, int col) const;
    // point getXY(int row, int col) const;
    // point getLatlong(int row, int col) const;
		
		
    // logical operators for grid cell
    bool onEdge(int row, int col); // if a grid cell is on the edges
    bool inGrid(int row, int col); // if a grid cell is within the grid boundaries

    // if the attributes of two grid maps are same 
    bool ifOverlay(const Grid& g1, const Grid& g2); 

    // assignment operators
    Grid operator = (const Grid& g);
    // Grid operator + (const Grid& g);
    // Grid operator - (const Grid& g);
    Grid operator += (const Grid& g);
    Grid operator -= (const Grid& g);
    Grid operator *= (const Grid& g);
    Grid operator /= (const Grid& g);
    Grid operator += (const double scalar);
    Grid operator -= (const double scalar);
    Grid operator *= (const double scalar);
    Grid operator /= (const double scalar);

    // arithematic operators
    friend Grid operator + (const Grid& g1, const Grid& g2);
    friend Grid operator + (const Grid& g1, const double scalar);
    friend Grid operator + (const double scalar, const Grid& g1);
    friend Grid operator - (const Grid& g1, const Grid& g2);
    friend Grid operator - (const Grid& g1, const double scalar);
    friend Grid operator - (const double scalar, const Grid& g1);
    // friend Grid operator * (const Grid& g1, const Grid& g2);
    friend Grid operator * (const Grid& g1, const double scalar);
    friend Grid operator * (const double scalar, const Grid& g1);
    friend Grid operator / (const Grid& g1, const Grid& g2);
    friend Grid operator / (const Grid& g1, const double scalar);
    friend Grid operator / (const double scalar, const Grid& g1);
		
    // logical operators
    friend Grid operator && (const Grid& g1, const Grid& g2);
    friend Grid operator || (const Grid& g1, const Grid& g2);

    // relational operators
    friend bool operator == (const Grid& g1, const Grid& g2);
    friend bool operator != (const Grid& g1, const Grid& g2);
    friend bool operator >  (const Grid& g1, const Grid& g2);
    friend bool operator <  (const Grid& g1, const Grid& g2);
    friend bool operator >= (const Grid& g1, const Grid& g2);
    friend bool operator <= (const Grid& g1, const Grid& g2);

    // mathematic operators
    // friend Grid abs  (const Grid& g);
    // friend Grid cos  (const Grid& g);
    friend Grid cosh (const Grid& g);
    friend Grid exp  (const Grid& g);
    friend Grid log  (const Grid& g);
    friend Grid log10 (const Grid& g);
    // friend Grid pow  (const Grid& g1, const TYPE scalar);
    // friend Grid sin  (const Grid& g);
    friend Grid sinh (const Grid& g);
    // friend Grid sqrt (const Grid& g);
    // friend Grid tan  (const Grid& g);
    friend Grid tanh (const Grid& g);

    // topographic analysis
    friend Grid slope(const Grid& g);
    friend Grid aspect(const Grid& g);
    friend Grid shade(const Grid& g, 
                      double solarElevation, double solarAzimuth, int traceDepth=6);
    friend Grid Viso(const Grid& g, int traceDepth=6);
    friend Grid sumShapeFactor(const Grid& g, int traceDepth=6);
				
    // import and export
    bool readGrid (int fileType, char* fileName, Grid& m_grid);
    bool writeGrid(int fileType,char* fileName, Grid m_grid);
		
		
    // file IO
    friend Grid readBinary(char* fileName, int dataType, const int m_rows, const int m_cols,
                           const double m_resolution, const double m_llX=0.0, const double m_llY=0.0, 
                           const double m_noValueData=(double)-9999);
		
    bool writeBinary(char* fileName);
		
    bool readASCII (char* fileName);
    bool writeASCII(char* fileName);
		
    bool readIdrisi (char* fileName);
    bool writeIdrisi(char* fileName);
		
    bool readArcinfoGrid (char* fileName);
    bool writeArcinfoGrid (char* fileName);
		
    bool readERDAS (char* fileName);
    bool writeERDAS(char* fileName);
		
    bool readENVI (char* fileName, Grid& m_grid);
		
    // default data type is float
    friend bool writeENVI(char* fileName, const Grid& g, int dataType=4);

    bool readHDFEOS (char* fileName);
    bool writeHDFEOS(char* fileName);

    double shapeFactor(int k, int l, double r, double m_z1, double m_z2, 
                       double m_s1, double m_s2, double m_a1, double m_a2);

		
  };


  //template <class TYPE>
  Grid::Grid()
  {
    this->rows			= 0;
    this->cols			= 0;
    this->resolution	= 0;
    this->llX			= 0.0;
    this->llY			= 0.0;
    this->urX			= (double)(cols-1);
    this->urY			= (double)(rows-1);
    this->noValueData	= (double)(-9999);
  }

  //template <class TYPE>
  Grid::Grid(const int m_rows, const int m_cols)
  {
    this->rows			= m_rows;
    this->cols			= m_cols;
    this->resolution	= 1.0;
    this->llX			= 0.0;
    this->llY			= 0.0;
    this->urX			= (double)(cols-1);
    this->urY			= (double)(rows-1);
    this->noValueData	= (double)(-9999);
		
    data = new double* [rows];
    for (int i=0; i<rows; i++)
      {
        data[i] = new double [cols];
      }
  }

  //template <class TYPE>
  Grid::Grid(const int m_rows, const int m_cols, const double m_resolution)
  {
    this->rows			= m_rows;
    this->cols			= m_cols;
    this->resolution	= m_resolution;
    this->llX			= 0.0;
    this->llY			= 0.0;
    this->urX			= (double)(cols-1) * resolution;
    this->urY			= (double)(rows-1) * resolution;
    this->noValueData	= (double)(-9999);
		
    data = new double* [rows];
    for (int i=0; i<rows; i++)
      {
        data[i] = new double [cols];
      }
  }
	
  //template <class TYPE>
  Grid::Grid(const int m_rows, const int m_cols, const double m_resolution,
             const double m_llX, const double m_llY)
  {
    this->rows			= m_rows;
    this->cols			= m_cols;
    this->resolution	= m_resolution;
    this->llX			= m_llX;
    this->llY			= m_llY;
    this->urX			= llX + (double)(cols-1) * resolution;
    this->urY			= llY + (double)(rows-1) * resolution;
    this->noValueData	= (double)(-9999);
		
    data = new double* [rows];
    for (int i=0; i<rows; i++)
      {
        data[i] = new double [cols];
      }
  }

  //template <class TYPE>
  Grid::Grid(const int m_rows, const int m_cols, const double m_llX, const double m_llY,
             const double m_urX, const double m_urY)
  {
    this->rows			= m_rows;
    this->cols			= m_cols;
    this->llX			= m_llX;
    this->llY			= m_llY;
    this->urX			= m_urX;
    this->urY			= m_urY;
    // deal with the reception caused by (urX-llx)!=(urY-llY)
    this->resolution	= (urX - llX) / cols;
    this->noValueData	= (double)(-9999);
		
    data = new double* [rows];
    for (int i=0; i<rows; i++)
      {
        data[i] = new double [cols];
      }
  }

  //template <class TYPE>
  Grid::Grid(const Grid& g)
  {
    this->rows			= g.rows;
    this->cols			= g.cols;
    this->resolution	= g.resolution;
    this->llX			= g.llX;
    this->llY			= g.llY;
    this->urX			= g.urX;
    this->urY			= g.urY;
    this->noValueData	= g.noValueData;
		
    data=new double*[rows];
    for(int i=0; i<rows; i++)
      {
        data[i] = new double[cols];
      }

    for(int row=0; row<rows; row++)
      {
        for(int col=0; col<cols; col++)
          {
            data[row][col]=g.GetAt(row, col);
          }
      }
  }

  //template <class TYPE>
  Grid::~Grid()
  {
    for (int i=0; i < rows; i++)
      delete [] data[i];
    delete [] data;
  }

  //template <class TYPE>
  void Grid::SetDimension(int m_rows, int m_cols, const double m_resolution, 
                          const double m_llX, const double m_llY)
  {
    this->rows			= m_rows;
    this->cols			= m_cols;
    this->resolution	= m_resolution;
    this->llX			= m_llX;
    this->llY			= m_llY;
    this->urX			= llX + (double)(cols-1) * resolution;
    this->urY			= llY + (double)(rows-1) * resolution;
    this->noValueData	= (double)(-9999);
		
    data = new double* [rows];
    for (int i=0; i < rows; i++)
      {
        data[i] = new double [cols];
      }
  }

  //template <class TYPE>
  void Grid::SetResolution(double m_resolution)
  {
    this->resolution = m_resolution;
  }

  //template <class TYPE>
  void Grid::SetNoValueData(double m_noValueData)
  {
    this->noValueData = m_noValueData;
  }

  //template <class double>
  int Grid::GetRows() const
  {
    return rows;
  }

  //template <class TYPE>
  int Grid::GetCols() const
  {
    return cols;
  }

  //template <class TYPE>
  double Grid::GetResolution() const
  {
    return resolution;
  }
	
  //template <class TYPE>
  double Grid::GetAt(int row, int col) const
  {
    return data[row][col];
  }
	
  //template <class TYPE>
  void Grid::SetAt(int row, int col, double newElement)
  {
    data[row][col] = newElement;
  }
	
  //template <class TYPE>
  bool Grid::onEdge(int row, int col)
  {
    if ( row==0 || row==rows-1 || col==0 || col==cols-1 )
      return true;
    else
      return false;
  }

  //template <class TYPE>
  bool Grid::inGrid(int row, int col)
  {
    if ( row<0 || row>=rows || col<0 || col>=cols )
      return false;
    else
      return true;
  }

  //template <class TYPE>
  bool Grid::ifOverlay(const Grid& g1, const Grid& g2)
  {
    if ( g1.rows==g2.rows && g1.cols==g2.cols && g1.resolution==g2.resolution 
         && g1.llX==g2.llX && g1.llY==g2.llY )
      return true;
    else
      return false;
  }

  //template <class TYPE>
  Grid Grid::operator = (const Grid& g)
  {
    rows		= g.rows;
    cols		= g.cols;
    llX			= g.llX;
    llY			= g.llY;
    urX			= g.urX;
    urY			= g.urY;
    resolution	= g.resolution;
    noValueData	= g.noValueData;
		
    data = new double* [rows];
    for (int i=0; i<rows; i++)
      {
        data[i] = new double [cols];
      }
		
    for (int i=0; i<rows; i++)
      {
        for (int j=0; j<cols; j++)
          {
            SetAt( i, j, (g.GetAt(i,j)) );
          }
      }
    return (*this);
  }


  //template <class TYPE>
  Grid operator + (const Grid& g1, const Grid& g2)
  {
    Grid gTmp(g1.rows, g1.cols, g1.resolution, g1.llX, g1.llY);

    for (int i=0; i<gTmp.rows; i++)
      {
        for (int j=0; j<gTmp.cols; j++)
          {
            gTmp.SetAt(i, j, g1.GetAt(i, j) + g2.GetAt(i, j) );
          }
      }

    return gTmp;
  }

	
  //template <class TYPE>
  Grid operator - (const Grid& g1, const Grid& g2)
  {
    Grid gTmp(g1.rows, g1.cols, g1.resolution, g1.llX, g1.llY);

    for (int i=0; i<gTmp.rows; i++)
      {
        for (int j=0; j<gTmp.cols; j++)
          {
            gTmp.SetAt(i, j, g1.GetAt(i, j) - g2.GetAt(i, j) );
          }
      }

    return gTmp;
  }


  //template <class TYPE>
  Grid operator * (const double scalar, const Grid& g1)
  {
    Grid gTmp(g1.rows, g1.cols, g1.resolution, g1.llX, g1.llY);

    for (int i=0; i<gTmp.rows; i++)
      {
        for (int j=0; j<gTmp.cols; j++)
          {
            gTmp.SetAt(i, j, scalar * g1.GetAt(i, j) );
          }
      }

    return gTmp;
  }


  //template <class TYPE>
  Grid operator / (const Grid& g1, const Grid& g2)
  {
    Grid gTmp(g1.rows, g1.cols, g1.resolution, g1.llX, g1.llY);

    for (int i=0; i<gTmp.rows; i++)
      {
        for (int j=0; j<gTmp.cols; j++)
          {
            gTmp.SetAt(i, j, g1.GetAt(i, j) / g2.GetAt(i, j) );
          }
      }

    return gTmp;
  }

	
  //template <class TYPE>
  Grid operator / (const Grid& g1, const double scalar)
  {
    Grid gTmp(g1.rows, g1.cols, g1.resolution, g1.llX, g1.llY);

    for (int i=0; i<gTmp.rows; i++)
      {
        for (int j=0; j<gTmp.cols; j++)
          {
            gTmp.SetAt(i, j, g1.GetAt(i, j) / scalar );
          }
      }

    return gTmp;
  }



	
	
  //template <class TYPE>
  Grid slope(const Grid& g)
  {
    Grid slopeGrid(g.rows, g.cols, g.resolution, g.llX, g.llY);
    slopeGrid.SetNoValueData(-1.0);
			
    double slope;
	
    for (int i=0; i<slopeGrid.rows; i++)
      {
        for (int j=0; j<slopeGrid.cols; j++)
          {
            if ( slopeGrid.onEdge(i, j) )
              {
                slope = (double)slopeGrid.noValueData;
                slopeGrid.SetAt(i, j, (double)slope);
              }
            else
              {
                // 以下算法基于“ERDAS FIELD GUIDE”322页
                double dx1 = (double)( g.GetAt(i-1,j+1) - g.GetAt(i-1,j-1) );
                double dx2 = (double)( g.GetAt(i,j+1)   - g.GetAt(i,j-1) );
                double dx3 = (double)( g.GetAt(i+1,j+1) - g.GetAt(i+1,j-1) );
                double dy1 = (double)( g.GetAt(i-1,j-1) - g.GetAt(i+1,j-1) );
                double dy2 = (double)( g.GetAt(i-1,j)	- g.GetAt(i+1,j) );
                double dy3 = (double)( g.GetAt(i-1,j+1) - g.GetAt(i+1,j+1) );
					
                double dx = (dx1+dx2+dx3) / (3*g.resolution);
                double dy = (dy1+dy2+dy3) / (3*g.resolution);
					
                double s  = sqrt( dx*dx+dy*dy ) / 2.0;
                slope = atan(s) * 180.0 / PI;
                slopeGrid.SetAt(i, j, (double)slope);
              }
          }
      }
		
    return slopeGrid;
  }

	
  // North = 0.0, clockwise
  //template<class double>
  Grid aspect(const Grid& g)
  {
    Grid aspectGrid(g.rows, g.cols, g.resolution, g.llX, g.llY);
    double aspect;

    aspectGrid.SetNoValueData(-1.0);

    for (int i=0; i<aspectGrid.rows; i++)
      {
        for (int j=0; j<aspectGrid.cols; j++)
          {
            if ( aspectGrid.onEdge(i, j) )
              {
                aspect = (double)aspectGrid.noValueData;
                aspectGrid.SetAt(i, j, (double)aspect);
              }
            else
              {
                // 以下算法基于“ERDAS FIELD GUIDE”322页
                double dx1 = (double)( g.GetAt(i-1,j+1) - g.GetAt(i-1,j-1) );
                double dx2 = (double)( g.GetAt(i,j+1)   - g.GetAt(i,j-1) );
                double dx3 = (double)( g.GetAt(i+1,j+1) - g.GetAt(i+1,j-1) );
                double dy1 = (double)( g.GetAt(i-1,j-1) - g.GetAt(i+1,j-1) );
                double dy2 = (double)( g.GetAt(i-1,j)	- g.GetAt(i+1,j) );
                double dy3 = (double)( g.GetAt(i-1,j+1) - g.GetAt(i+1,j+1) );

                double dx  = dx1+dx2+dx3;
                double dy  = dy1+dy2+dy3;
		  
                if (dy!=0.0)
                  {
                    aspect = atan(dx/dy);
                    if (dx>0.0 && dy>0.0) 
                      aspect += PI;
                    else if (dx>0.0 && dy<0.0)
                      aspect += 2*PI;
                    else if (dx<0.0 && dy>0.0)
                      aspect += PI;
                    else if (dx<0.0 && dy<0.0)
                      aspect += 0.0;
                    else if (dx==0.0 && dy>0.0)
                      aspect += PI;
                    else // (dx==0.0 && dy<0.0)
                      aspect += 0.0;
                  }
                else
                  {
                    if (dx>0.0)
                      aspect = 1.5*PI;
                    else if (dx<0.0)
                      aspect = 0.5*PI;
                    else
                      // when terrain is flat, set the slope as facing south
                      aspect = PI;
                  }
					
                if ( aspect>=0.0 )
                  aspect = aspect * 180 / PI;
					
                aspectGrid.SetAt(i, j, (double)aspect);
              }
          }
      }
		
    return aspectGrid;
  }

	
  template<class TYPE>
  Grid shade(const Grid& g, 
             double solarElevation, double solarAzimuth, int traceDepth)
  {
    // solar azimuth's definition: north=0, clockwise
    solarElevation	*= PI/180.0;
    solarAzimuth	*= PI/180.0;
		
    Grid shadeGrid(g.rows, g.cols, g.resolution, g.llX, g.llY);
    shadeGrid.SetNoValueData(-1.0);

    double shade;
    double p1, p2, p3, p, dh;
    int nextRow, nextCol;

    p1 = fabs( cos(solarAzimuth) );
    p2 = fabs( sin(solarAzimuth) );
    p3 = (p1>p2) ? p1 : p2;

    for (int i=0; i<shadeGrid.rows; i++)
      {
        for (int j=0; j<shadeGrid.cols; j++)
          {
            shade = 1.0;
            for (int l=1; l<=traceDepth; l++)
              {
                p = double(l)/p3; 
                if ( p*cos(solarAzimuth)>=0.0 )
                  nextRow = i - int(p*cos(solarAzimuth)+0.5);
                else
                  nextRow = i - int(p*cos(solarAzimuth)-0.5);
                if (p*sin(solarAzimuth)>=0.0)
                  nextCol	= j + int(p*sin(solarAzimuth)+0.5);
                else
                  nextCol	= j + int(p*sin(solarAzimuth)-0.5);
			  
                if ( shadeGrid.inGrid(nextRow, nextCol) )
                  {
                    dh = (double)(g.GetAt(nextRow, nextCol) - g.GetAt(i, j));
                    if ( atan(dh/(shadeGrid.resolution*p)) > solarElevation )
                      {
                        shade = 0.0;
                        break;
                      }
                  }
              }
            shadeGrid.SetAt(i, j, (double)shade);
          }
      }

    return shadeGrid;
  }


  template<class TYPE>
  Grid Viso(const Grid& g, int traceDepth)
  {
    Grid VisoGrid(g.rows, g.cols, g.resolution, g.llX, g.llY);
    VisoGrid.SetNoValueData(-1.0);
		
    double k;			// percentage of unobstructed area at a given direction
    double sinh, sinhh;
    double viso;
    double solarAzimuth;

    double p1, p2, p3, p, dh;
    int nextRow, nextCol;

    for (int i=0; i<VisoGrid.rows; i++)
      {
        for (int j=0; j<VisoGrid.cols; j++)
          {
            viso = 0;
            for (int direction=0; direction<8; direction++)
              {
                solarAzimuth = direction * PI / 4.0;
                p1 = fabs( cos(solarAzimuth) );
                p2 = fabs( sin(solarAzimuth) );
                p3 = (p1>p2) ? p1 : p2;
				
                sinh = 0.0;
				
                for (int l=1; l<=traceDepth; l++)
                  {
                    p = double(l)/p3;
                    if ( p*cos(solarAzimuth)>=0.0 )
                      nextRow = i - int(p*cos(solarAzimuth)+0.5);
                    else
                      nextRow = i - int(p*cos(solarAzimuth)-0.5);
                    if (p*sin(solarAzimuth)>=0.0)
                      nextCol	= j + int(p*sin(solarAzimuth)+0.5);
                    else
                      nextCol	= j + int(p*sin(solarAzimuth)-0.5);
												
                    if ( VisoGrid.inGrid(nextRow, nextCol) )
						{
							dh = (double)(g.GetAt(nextRow, nextCol) - g.GetAt(i, j));
							sinhh = dh / sqrt(p*p*VisoGrid.resolution*VisoGrid.resolution + dh*dh);
						}
						else
							sinhh = 0.0;
						
						if (sinhh > sinh) 
							sinh = sinhh;
						
						k = 1 - sinh;
					}
					viso += k;
				}
				viso /= 8.0;
				VisoGrid.SetAt(i, j, (double)viso);
			}
		}

		return VisoGrid;
	}


	// calculate the sum value of shape factor of a grid
  //template<class TYPE>
	Grid sumShapeFactor(const Grid& g, int traceDepth)
	{
		Grid sfGrid(g.rows, g.cols, g.resolution, g.llX, g.llY);
		sfGrid.SetNoValueData(-1.0);

		Grid slopeGrid = slope(g);
		Grid aspectGrid = aspect(g);

		double sf, sfSum;

		for(int i=0; i<g.rows; i++)
		{
			for(int j=0; j<g.cols; j++)
			{
				sfSum  = 0.0;
				for (int k=-traceDepth; k<=traceDepth; k++)
				{
					for (int l=-traceDepth; l<=traceDepth; l++)
					{
						int nextRow	= i+k;
						int nextCol	= j+l;
						if ( sfGrid.inGrid(nextRow, nextCol) )
						{
							sf = sfGrid.shapeFactor ( k, l, g.resolution, 
								g.GetAt(i, j), g.GetAt(nextRow, nextCol), 
								slopeGrid.GetAt(i, j)*PI/180,
								slopeGrid.GetAt(nextRow, nextCol)*PI/180, 
								aspectGrid.GetAt(i, j)*PI/180, 
								aspectGrid.GetAt(nextRow, nextCol)*PI/180 );
							if (sf>0.0)
								sfSum += sf;
						}
					}
				}
				sfGrid.SetAt(i, j, (double)sfSum);
			}
		}
		
		return sfGrid;
	}


  //template<class TYPE>
	double Grid::shapeFactor(int k, int l, double r,
		double m_z1, double m_z2, 
		double m_s1, double m_s2, double m_a1, double m_a2)
	{
		double sf;		// shape factor
		
		if ( k==0 && l==0)
			return 0;
		
		else
		{
			// vector between v1(row, col) and v2(nextRow, nextCol), v2-v1
			Vector3D v ( l*r, -k*r, m_z2-m_z1 );
			v	= v.Normalize();
			
			// normal vector 1
			double s1 = m_s1;
			double a1;
			if (m_a1<=PI/2)
				a1 = PI/2 - m_a1;
			else
				a1 = 2*PI + PI/2 - m_a1;
			Vector3D n1( sin(s1)*cos(a1), sin(s1)*sin(a1), cos(s1) );
			
			// normal vector 2
			double s2 = m_s2;
			double a2;
			if (m_a2<=PI/2)
				a2 = PI/2 - m_a2;
			else
				a2 = 2*PI + PI/2 - m_a2;
			Vector3D n2( sin(s2)*cos(a2), sin(s2)*sin(a2), cos(s2) );
			
			double cosf1	= Dot(v, n1);
                        Vector3D v_m(-v);
			double cosf2	= Dot(v_m, n2);
			
			double dis		= fabs((double)k)*fabs((double)k) 
                          + fabs((double)l)*fabs((double)l)
                          + (m_z2 - m_z1) * (m_z2 - m_z1) / (r * r);
			
			if (cosf1 >= 0.0 && cosf2 >=0.0)
			{
				sf = cosf1 * cosf2 / ( PI*dis );
			}
			else 
				sf = 0.0;
			
			return sf;
		}
	}

	
	template<class TYPE>
	Grid readBinary(char* fileName, int dataType, 
		const int m_rows, const int m_cols,
		const double m_resolution, const double m_llX, const double m_llY, 
		const double m_noValueData)
	{
		Grid gTmp;
		gTmp.SetDimension(m_rows, m_cols, m_resolution, m_llX, m_llY);
		gTmp.SetNoValueData(m_noValueData);

		unsigned char dataByte;
		short int     dataShort;
		int			  dataInt;
		float		  dataFloat;
		double		  dataDouble;

		FILE *binaryFile;
		if ( (binaryFile=fopen(fileName, "rb"))==NULL )
		{
			printf("can not open the general binary file\n");
			exit(1);
		}

		for (int i=0; i<gTmp.rows; i++)
		{
			for (int j=0; j<gTmp.cols; j++)
			{
				if (dataType==Grid::BYTE)
				{
					fread (&dataByte, sizeof(unsigned char), 1, binaryFile);
					gTmp.SetAt(i, j, (double)dataByte);
				}
				else if (dataType==Grid::SHORT)
				{
					fread (&dataShort, sizeof(short int), 1, binaryFile);
					gTmp.SetAt(i, j, (double)dataShort);
				}
				else if (dataType==Grid::INT)
				{
					fread (&dataInt, sizeof(int), 1, binaryFile);
					gTmp.SetAt(i, j, (double)dataInt);
				}
				else if (dataType==Grid::FLOAT)
				{
					fread (&dataFloat, sizeof(float), 1, binaryFile);
					gTmp.SetAt(i, j, (double)dataFloat);
				}
				else if (dataType==Grid::DOUBLE)
				{
					fread (&dataDouble, sizeof(double), 1, binaryFile);
					gTmp.SetAt(i, j, (double)dataDouble);
				}
				else
					break;
			}
		}

		fclose(binaryFile);
	
		return (gTmp);
	}


	template<class TYPE>
	bool writeENVI(char* fileName, const Grid& g, int dataType)
	{
		FILE *imageFile;
		
		unsigned char dataByte;
		short int     dataShort;
		int			  dataInt;
		float		  dataFloat;
		double		  dataDouble;
	
		if ( (imageFile=fopen(fileName, "wb"))==NULL )
		{
			printf("can not open the general image file for writing\n");
			exit(1);
		}

		for (int i=0; i<g.rows; i++)
		{
			for (int j=0; j<g.cols; j++)
			{
				if (dataType==Grid::BYTE)
				{
					dataByte = (unsigned char)g.GetAt(i, j);
					fwrite (&dataByte, sizeof(unsigned char), 1, imageFile);
				}
				else if (dataType==Grid::SHORT)
				{
					dataShort = (short int)g.GetAt(i, j);
					fwrite (&dataShort, sizeof(short int), 1, imageFile);
				}
				else if (dataType==Grid::INT)
				{
					dataInt = (int)g.GetAt(i, j);
					fwrite (&dataInt, sizeof(int), 1, imageFile);
				}
				else if (dataType==Grid::FLOAT)
				{
					dataFloat = (float)g.GetAt(i, j);
					fwrite (&dataFloat, sizeof(float), 1, imageFile);
				}
				else if (dataType==Grid::DOUBLE)
				{
					dataDouble = (double)g.GetAt(i, j);
					fwrite (&dataDouble, sizeof(double), 1, imageFile);
				}
				else
					break;
			}
		}

		fclose(imageFile);
	
		return (true);
	}






				

	/*
	template<class TYPE>
	bool Grid::readENVI(char* fileName, Grid& m_grid)
	{
		TRY
		{
			CFile FileRead(FileName,CFile::modeRead);
			Dint unsigned charsRead = FileRead.ReadHuge(pData,Size*sizeof(TYPE));
			FileRead.Close();
		}
		CATCH( CFileException, e)
		{
			#ifdef _DEBUG
			afxDump << "File could not be opened " << e->m_cause << "\n";
			#endif
		}
		END_CATCH
	}
	
	template<class TYPE>
	void Grid::WriteENVI(char* FileName)
	{
	
		TRY
		{
			CFile FileWrite(FileName,CFile::modeCreate | CFile::modeWrite);
			FileWrite.WriteHuge(pData,Size*sizeof(TYPE));
			FileWrite.Close();
		}
		CATCH( CFileException, e)
		{
			#ifdef _DEBUG
			afxDump << "File could not be opened " << e->m_cause << "\n";
		#endif
		}
		END_CATCH
	}

	*/	
	

}

#endif
