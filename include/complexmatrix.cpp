namespace ldas
{
template <class T> ComplexMatrix<T>::ComplexMatrix(const T* RealBuf, const T* ImagBuf,const unsigned long row, const unsigned long col)
{
    m_row=row;
    m_col=col;
    m_Tolerance=(T)(1.0e-8);
    m_RealBuf= new T[m_row*m_col];
    m_ImagBuf= new T[m_row*m_col];
    unsigned long i,j;

    for(i=1; i<=m_row; i++)
        for(j=1; j<=m_col; j++)
        {
            set(i,j,RealBuf[(i-1)*m_col+j-1],ImagBuf[(i-1)*m_col+j-1]);
        }
}
template <class T> ComplexMatrix<T>::ComplexMatrix()
{
    m_row=0;
    m_col=0;
    m_Tolerance=(T)(1.0e-8);
    m_RealBuf= NULL;
    m_ImagBuf= NULL;
}
template <class T> ComplexMatrix<T>::ComplexMatrix(const unsigned long row, const unsigned long col)
{
    m_row=row;
    m_col=col;
    m_Tolerance=(T)(1.0e-8);
    m_RealBuf= new T[m_row*m_col];
    m_ImagBuf= new T[m_row*m_col];
    unsigned long k;
    for(k=0; k<m_row*m_col; k++)
    {
        m_RealBuf[k]=(T)0.0;
        m_ImagBuf[k]=(T)0.0;
    }
}
template <class T> ComplexMatrix<T>::ComplexMatrix(const ComplexMatrix<T>& cmatrix)
{
    m_row=cmatrix.row();
    m_col=cmatrix.column();
    m_Tolerance=(T)(1.0e-8);
    m_RealBuf= new T[m_row*m_col];
    m_ImagBuf= new T[m_row*m_col];
    unsigned long row, col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            this->set(row,col,cmatrix(row,col));
        }
}
template <class T> ComplexMatrix<T>::ComplexMatrix(const char *FileName)
{
    ifstream in;
    in.open(FileName);
    if(!in)//不能打开文件
        throw Exception("ComplexMatrix","ComplexMatrix(const char *FileName)","Can not open data file");
    unsigned long row,col;
    in>>this->m_row;
    in>>this->m_col;
    this->m_RealBuf=new T[m_row*m_col];
    this->m_ImagBuf=new T[m_row*m_col];
    this->m_Tolerance =(T)1.0e-8;
    complex<T> tmp;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            in>>tmp;
            set(row,col,tmp);
        }
    in.close();
}
template <class T> ComplexMatrix<T>::ComplexMatrix(const unsigned long row,const unsigned long col,const Matrix<T> real,const Matrix<T> imag)
{
    if((row==imag.row())&&(row==real.row())&&(col==imag.column())&&(col==real.column()))
    {
        m_row=row;
        m_col=col;
        m_Tolerance=(T)(1.0e-8);
        m_RealBuf= new T[m_row*m_col];
        m_ImagBuf= new T[m_row*m_col];
        unsigned long i,j;
        for(i=1; i<=m_row; i++)
            for(j=1; j<=m_col; j++)
            {
                this->set(i,j,real.get(i,j),imag.get(i,j));
            }
    }
    else
    {
        throw Exception("ComplexMatrix","ComplexMatrix(unsigned long row,unsigned long col,Matrix<T> real,Matrix<T> imag)","Dimension mismatch");
    }
}
template <class T> void ComplexMatrix<T>::setDim(const unsigned long row,const unsigned long col)
{
    m_row=row;
    m_col=col;
    m_Tolerance=(T)(1.0e-8);
    m_RealBuf= new T[m_row*m_col];
    m_ImagBuf= new T[m_row*m_col];
    unsigned long k;
    for(k=0; k<m_row*m_col; k++)
    {
        m_RealBuf[k]=(T)0.0;
        m_ImagBuf[k]=(T)0.0;
    }
}

template <class T> void ComplexMatrix<T>::save(const char *FileName)
{
    ofstream out;
    out.open(FileName);
    if(!out)
        throw Exception("ComplexMatrix","save","Can not create data file");
    unsigned long row,col;
    out<<this->m_row<<"	"<<this->m_col<<endl;
    for(row=1; row<=m_row; row++)
    {
        for(col=1; col<=m_col; col++)
        {
            out.setf(ios_base::showpoint);
            out<<this->get(row,col)<<" ";
        }
        out<<endl;
    }
    out.close();
}
template <class T> ComplexMatrix<T>::~ComplexMatrix()
{
    delete[] m_RealBuf;
    delete[] m_ImagBuf;
}
template <class T> unsigned long ComplexMatrix<T>::row() const
{
    return m_row;
}
template <class T> unsigned long ComplexMatrix<T>::column() const
{
    return m_col;
}
template <class T> ComplexMatrix<T>& ComplexMatrix<T>::neg()
{
    unsigned long row, col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            this->set(row,col,(-this->get(row,col))) ;
        }
    return *this;
}
template <class T> T ComplexMatrix<T>::ValueAbs(T value) const
{
    if(value>=0)
        return value;
    else
        return -value;
}
template <class T> T ComplexMatrix<T>::MaxAbs(unsigned long& row, unsigned long& col, unsigned long k)
{
    T a=(T)0.0;
    T br,bi,b;
    unsigned long i, j;
    for(i=1; i<=m_row; i++)
        for(j=1; j<=m_col; j++)
        {
            br=this->get(i,j).real();
            bi=this->get(i,j).imag();
            b=br*br+bi*bi;
            if(a<b)
            {
                a=b;
                row=i;
                col=j;
            }
        }
    return a;
}
template <class T> T ComplexMatrix<T>::Adjust(T a)
{
    T b=floor(a*10000.0+0.5)/10000.0;
    if(ValueAbs(b-a)<this->m_Tolerance)//如果舍入误差很小,小于容许误差，则抛弃
        a=b;
    return a;
}
template <class T> void ComplexMatrix<T>::DoAdjust()
{
    unsigned long k;
    for(k=0; k<=(m_row*m_col); k++)
    {
        this->m_ImagBuf[k]=this->Adjust(this->m_ImagBuf[k]);
        this->m_RealBuf[k]=this->Adjust(this->m_RealBuf[k]);
    }
}
template <class T> ComplexMatrix<T> ComplexMatrix<T>::operator -() const
{
    ComplexMatrix<T> temp(*this);
    return temp.neg();
}
template <class T> ComplexMatrix<T>& ComplexMatrix<T>::tran()
{
    unsigned long row, col;
    ComplexMatrix<T> cm(*this);
    row=this->m_row;
    this->m_row =this->m_col;
    this->m_col =row;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            this->set(row,col,cm.get(col,row));
        }
    return *this;
}
template <class T> ComplexMatrix<T> ComplexMatrix<T>::tran() const
{
    ComplexMatrix<T> temp(*this);
    return temp.tran();
}
template <class T> ComplexMatrix<T> ComplexMatrix<T>::conj_tran() const
{
    ComplexMatrix<T> temp(*this);
    temp.conj();
    return temp.tran();
}
template <class T> void ComplexMatrix<T>::set(const unsigned long row,const unsigned long col, const T real, const T image)
{
    this->m_RealBuf[(row-1)*m_col+col-1]=real;
    this->m_ImagBuf[(row-1)*m_col+col-1]=image;
}
template <class T> void ComplexMatrix<T>::set(const unsigned long row, const unsigned long col,const complex<T>& c)
{
    this->m_RealBuf[(row-1)*m_col+col-1]=c.real();
    this->m_ImagBuf[(row-1)*m_col+col-1]=c.imag();
}
template <class T> complex<T> ComplexMatrix<T>::get(const unsigned long row, const unsigned long col) const
{
    if(row>m_row||col>m_col)//如果下标超过矩阵的规模，抛出异常
        throw Exception("ComplexMatrix","get","row or col exceed thd size of matrix");
    T real,image;
    real=this->m_RealBuf[(row-1)*m_col+col-1];
    image=this->m_ImagBuf[(row-1)*m_col+col-1];
    return complex<T>(real,image);
}
template <class T> ComplexMatrix<T>& ComplexMatrix<T>::conj()
{
    unsigned long k;
    for(k=0; k<(m_row*m_col); k++)
        m_ImagBuf[k]*=(T)(-1);
    return *this;
}
template <class T> Matrix<T> ComplexMatrix<T>::real() const
{
    Matrix<T> m(this->row(),this->column());
    unsigned long row, col;
    for(row=1; row<=this->m_row; row++)
        for(col=1; col<=this->m_col; col++)
        {
            m.set(row,col,this->get(row,col).real());
        }
    return m;
}
template <class T> Matrix<T> ComplexMatrix<T>::imag() const
{
    Matrix<T> m(this->row(),this->column());
    unsigned long row, col;
    for(row=1; row<=this->m_row; row++)
        for(col=1; col<=this->m_col; col++)
        {
            m.set(row,col,this->get(row,col).imag());
        }
    return m;
}
template <class T> Matrix<T> ComplexMatrix<T>::arg() const
{
    Matrix<T> m(this->row(),this->column());
    unsigned long row, col;
    for(row=1; row<=this->m_row; row++)
        for(col=1; col<=this->m_col; col++)
        {
            m.set(row,col,arg(this->get(row,col)));
        }
    return m;
}
template <class T> complex<T> ComplexMatrix<T>::operator()(const unsigned long row, const unsigned long col) const//根据下标，返回一个复数
{
    T real, image;
    real=this->m_RealBuf[(row-1)*m_col+col-1];
    image=this->m_ImagBuf[(row-1)*m_col+col-1];
    return complex<T>(real,image);//返回一个complex<T>类型，应用real和image初始化的值
}
template <class T> ComplexMatrix<T>& ComplexMatrix<T>::operator =(const ComplexMatrix& rcmatrix)
{
    m_col=rcmatrix.column();
    m_row=rcmatrix.row();
    delete[] m_ImagBuf;
    delete[] m_RealBuf;
    m_ImagBuf=new T[m_row*m_col];
    m_RealBuf=new T[m_row*m_col];
    unsigned long row, col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            this->set(row,col,rcmatrix(row,col));
        }
    return *this;
}
template <class T> ostream& operator<<(ostream& os,const ComplexMatrix<T>& cmatrix)
{
    unsigned long row, col;
    for(row=1; row<=cmatrix.m_row; row++)
    {
        for(col=1; col<=cmatrix.m_col; col++)
        {
            os.setf(ios::right);//设置输出为左对齐
            os<<"(";
            os.width(13);
            os<<cmatrix.m_RealBuf[(row-1)*cmatrix.m_col+col-1];//(row,col).real();
            if(cmatrix.m_ImagBuf[(row-1)*cmatrix.m_col+col-1]>=0)
            {
                os<<" + ";
                os.setf(ios::left);
                os.width(13);
                os<<cmatrix.ValueAbs(cmatrix.m_ImagBuf[(row-1)*cmatrix.m_col+col-1])<<"i"<<")  ";
            }
            else
            {
                os<<" - ";
                os.setf(ios::left);
                os.width(13);
                os<<cmatrix.ValueAbs(cmatrix.m_ImagBuf[(row-1)*cmatrix.m_col+col-1])<<"i"<<")  ";
            }
        }
        os<<endl;
    }
    return os;
}
template <class T> ComplexMatrix<T>& ComplexMatrix<T>::operator +=(const ComplexMatrix<T>& rcmatrix)
{
    //实部与实部相加，虚部与虚部相加
    unsigned long row, col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            this->set(row,col, (rcmatrix(row,col)+this->get(row,col)));
        }
    return *this;
}
template <class T> ComplexMatrix<T> ComplexMatrix<T>::operator +(const ComplexMatrix<T>& rcmatrix) const
{
    if(this->m_row!=rcmatrix.row()||this->m_col!=rcmatrix.column())
        throw Exception("ComplexMatrix","operator +","the two matrice are not in same size, can not perform addition operation.");
    ComplexMatrix<T> temp(*this);
    return temp+=rcmatrix;
}
template <class T> ComplexMatrix<T>& ComplexMatrix<T>::operator -=(const ComplexMatrix<T>& rcmatrix)
{
    if(this->m_row!=rcmatrix.row||this->m_col!=rcmatrix.column())
        throw Exception("ComplexMatrix","operator -=","the two matrice are not in same size, can not perform minus operation.");
    return ((*this)+=(-rcmatrix));
}
template <class T> ComplexMatrix<T> ComplexMatrix<T>::operator -(const ComplexMatrix<T>& rcmatrix) const
{
    if(this->m_row!=rcmatrix.row||this->m_col!=rcmatrix.column())
        throw Exception("ComplexMatrix","operator -","the two matrice are not in same size, can not perform minus operation.");
    ComplexMatrix<T> temp(*this);
    return (temp-=rcmatrix);
}
template <class T> ComplexMatrix<T>& ComplexMatrix<T>::operator *=(const T k)
{
    unsigned long row;
    unsigned long col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            this->set(row,col,(k*this->get(row,col)));
        }
    return *this;
}
template <class T> ComplexMatrix<T> ComplexMatrix<T>::operator *(const T k) const
{
    ComplexMatrix<T> temp(*this);
    return temp*=k;
}
template <class T> ComplexMatrix<T>& ComplexMatrix<T>::operator /=(const T k)
{
    if(this->ValueAbs(k)<=this->m_Tolerance)
        throw Exception("ComplexMatrix","operator /=","divided by zero");
    (*this)*=(1/k);
    return *this;
}
template <class T> ComplexMatrix<T> ComplexMatrix<T>::operator /(const T k) const
{
    if(this->ValueAbs(k)<=this->m_Tolerance)
        throw Exception("ComplexMatrix","operator /","divided by zero");
    ComplexMatrix<T> temp(*this);
    return temp/=k;
}
template <class T> ComplexMatrix<T>& ComplexMatrix<T>::operator *=(const ComplexMatrix<T>& rcmatrix)
{
    if(this->m_col!=rcmatrix.row())//如果不满足矩阵相乘地条件
        throw Exception("ComplexMatrix","operator *=","The sizes of two matrices are mismatch");
    unsigned long row, col,k;
    ComplexMatrix<T> temp(*this);
    delete[] m_RealBuf;
    delete[] m_ImagBuf;
    m_RealBuf=new T[m_row*rcmatrix.column()];
    m_ImagBuf= new T[m_row*rcmatrix.column()];
    this->m_col=rcmatrix.column();
    for(row=1; row<=m_row; row++)
        for(col=1; col<=rcmatrix.column(); col++)
        {
            complex<T> c(0.0,0.0);
            for(k=1; k<=rcmatrix.row(); k++)
                c+=temp(row,k)*rcmatrix(k,col);
            this->set(row,col,c);
        }
    this->DoAdjust();
    return *this;
}
template <class T> ComplexMatrix<T> ComplexMatrix<T>::operator *(const ComplexMatrix<T>& rcmatrix) const
{
    if(this->m_col!=rcmatrix.row())//如果不满足矩阵相乘地条件
        throw Exception("ComplexMatrix","operator *","The sizes of two matrices are mismatch");
    ComplexMatrix<T> temp(*this);
    return (temp*=rcmatrix);
}
template <class T> ComplexMatrix<T>& ComplexMatrix<T>::operator ^=(const T x)
{
    unsigned long row, col;
    T value;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            value = this->get(row,col);
            if(value<0.0)
                this->set(row,col,0.0);
            else
                this->set(row,col,pow(value,x));
        }
    return *this;
}
template <class T> ComplexMatrix<T> ComplexMatrix<T>::operator ^(const T x) const
{
    ComplexMatrix<T> temp(*this);
    return temp^=x;
}
//复数矩阵的数学运算
template <class T> ComplexMatrix<T>& ComplexMatrix<T>::sin()
{
    unsigned long row,col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            this->set(row,col,sin(this->get(row,col)));
        }
    return *this;
}
template <class T> ComplexMatrix<T>& ComplexMatrix<T>::cos()
{
    unsigned long row,col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            this->set(row,col,cos(this->get(row,col)));
        }
    return *this;
}
template <class T> ComplexMatrix<T>& ComplexMatrix<T>::tan()
{
    unsigned long row,col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            this->set(row,col,tan(this->get(row,col)));
        }
    return *this;
}
template <class T> ComplexMatrix<T>& ComplexMatrix<T>::exp()
{
    unsigned long row,col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            this->set(row,col,exp(this->get(row,col)));
        }
    return *this;
}
template <class T> ComplexMatrix<T>& ComplexMatrix<T>::log()
{
    unsigned long row,col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            this->set(row,col,log(this->get(row,col)));
        }
    return *this;
}
template <class T> ComplexMatrix<T>& ComplexMatrix<T>::log10()
{
    unsigned long row,col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            this->set(row,col,log10(this->get(row,col)));
        }
    return *this;
}
template <class T> ComplexMatrix<T>& ComplexMatrix<T>::sinh()
{
    unsigned long row,col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            this->set(row,col,sinh(this->get(row,col)));
        }
    return *this;
}
template <class T> ComplexMatrix<T>& ComplexMatrix<T>::cosh()
{
    unsigned long row,col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            this->set(row,col,cosh(this->get(row,col)));
        }
    return *this;
}
template <class T> ComplexMatrix<T>& ComplexMatrix<T>::tanh()
{
    unsigned long row,col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            this->set(row,col,tanh(this->get(row,col)));
        }
    return *this;
}
template <class T> ComplexMatrix<T>& ComplexMatrix<T>::sqrt()
{
    unsigned long row,col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            this->set(row,col,std::sqrt(this->get(row,col)));
        }
    return *this;
}
template <class T> Matrix<T> ComplexMatrix<T>::abs() const
{
    Matrix<T> m(this->m_row,this->m_col);
    unsigned long row,col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            m.set(row,col,abs(this->get(row,col)));
        }
    return m;
}
template <class T> Matrix<T> ComplexMatrix<T>::norm() const
{
    Matrix<T> m(this->m_row,this->m_col);
    unsigned long row,col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            m.set(row,col,norm(this->get(row,col)));
        }
    return m;
}
//数学函数定义完毕
template <class T> void ComplexMatrix<T>::swapRows(const unsigned long row1,const unsigned long row2)
{
    if(row1>m_col||row1==0||row2>m_col||row2==0)//异常处理
        throw Exception("ComplexMatrix","swapCols","One or tow subscripts are our of range");
    unsigned long col;
    complex<T> temp(0,0);
    for(col=1; col<=m_col; col++)
    {
        temp=this->get(row1,col);
        this->set(row1,col,this->get(row1,col));
        this->set(row2,col,temp);
    }
}
template <class T> void ComplexMatrix<T>::swapCols(const unsigned long col1,const unsigned long col2)
{
    if(col1>m_col||col1==0||col2>m_col||col2==0)//异常处理
        throw Exception("ComplexMatrix","swapCols","one or tow subscripts are our of range");
    unsigned long row;
    complex<T> temp(0,0);
    for(row=1; row<=m_row; row++)
    {
        temp=this->get(row,col1);
        this->set(row,col1,this->get(row,col2));
        this->set(row,col2,temp);
    }
}
template <class T>ComplexMatrix<T>& ComplexMatrix<T>::inv()
{
    if(m_row!=m_col)//如果不是方阵，抛出异常
        throw Exception("ComplexMatrix","inv","the matrix is not square, can not inverse.");
    Matrix<T> m(2*m_row,2*m_col);
    unsigned long row,col;
    for(row=1; row<=m_row; row++) //左上角
        for(col=1; col<=m_col; col++)
            m.set(row,col,this->get(row,col).real());
    for(row=m_row+1; row<=2*m_row; row++) //左下角
        for(col=1; col<=m_col; col++)
            m.set(row,col,-this->get(row-m_row,col).imag());
    for(row=1; row<=m_row; row++) //右上角
        for(col=m_col+1; col<=2*m_col; col++)
            m.set(row,col,this->get(row,col-m_col).imag());
    for(row=m_row+1; row<=2*m_row; row++)
        for(col=m_col+1; col<=2*m_col; col++)
            m.set(row,col,this->get(row-m_row,col-m_col).real());
    m.inv();
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
            this->set(row,col,m(row,col),m(row,col+m_col));
    return *this;
}
template <class T> ComplexMatrix<T> ComplexMatrix<T>::operator ~() const
{
    if(m_row!=m_col)//如果不是方阵，抛出异常
        throw Exception("ComplexMatrix","operator ~","the matrix is not square, can not inverse.");
    ComplexMatrix<T> temp(m_row,m_col);
    Matrix<T> m(2*m_row,2*m_col);
    unsigned long row,col;
    for(row=1; row<=m_row; row++) //左上角
        for(col=1; col<=m_col; col++)
            m.set(row,col,this->get(row,col).real());
    for(row=m_row+1; row<=2*m_row; row++) //左下角
        for(col=1; col<=m_col; col++)
            m.set(row,col,-this->get(row-m_row,col).imag());
    for(row=1; row<=m_row; row++) //右上角
        for(col=m_col+1; col<=2*m_col; col++)
            m.set(row,col,this->get(row,col-m_col).imag());
    for(row=m_row+1; row<=2*m_row; row++)
        for(col=m_col+1; col<=2*m_col; col++)
            m.set(row,col,this->get(row-m_row,col-m_col).real());
    m.inv();
    //cout<<m<<endl;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
            temp.set(row,col,m(row,col),m(row,col+m_col));
    return temp;
}
template <class T> void ComplexMatrix<T>::genI()
{
    unsigned long row,col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            if(row==col)
                this->set(row,col,(T)1.0,(T)0.0);
            else
                this->set(row,col,(T)0.0,(T)0.0);
        }
}
template <class T> bool ComplexMatrix<T>::is_square() const
{
    return(m_row==m_col);
}
template <class T> bool ComplexMatrix<T>::is_symmetry() const
{
    unsigned long row,col;
    bool result=true;
    for(row=1; row<=m_row; row++)
        for(col=row+1; col<=m_col; col++)
        {
            if(ValueAbs(get(row,col)-get(col,row))>m_Tolerance)
                return false;
        }
}
/////function template for ComplexMatrix
template <class T> Matrix<T> abs(ComplexMatrix<T>& cm)//返回矩阵每个元素的模（modulus），The modulus of a complex number a + bi is sqrt(a^2 + b^2)
{
    return cm.abs();
}
template <class T> Matrix<T> norm(ComplexMatrix<T>& cm)//返回矩阵每一个元素的范数（norm）,The norm of a complex number a + bi is a^2 + b^2
{
    return cm.norm();
}
template <class T> ComplexMatrix<T> sin(ComplexMatrix<T>& cm)//对所有的元素取正弦，返回新的矩阵
{
    ComplexMatrix<T> temp(cm);
    return temp.sin();
}
template <class T> ComplexMatrix<T> cos(ComplexMatrix<T>& cm)//对所有元素取余弦，返回新的矩阵
{
    ComplexMatrix<T> temp(cm);
    return temp.cos();
}
template <class T> ComplexMatrix<T> tan(ComplexMatrix<T>& cm)//对所有元素取正切，返回新的矩阵
{
    ComplexMatrix<T> temp(cm);
    return temp.tan();
}
template <class T> ComplexMatrix<T> sinh(ComplexMatrix<T>& cm)//对所有元素取双曲正弦，返回新的矩阵
{
    ComplexMatrix<T> temp(cm);
    return temp.sinh();
}
template <class T> ComplexMatrix<T> cosh(ComplexMatrix<T>& cm)//对所有元素取双曲余弦，返回新的矩阵
{
    ComplexMatrix<T> temp(cm);
    return temp.cosh();
}
template <class T> ComplexMatrix<T> tanh(ComplexMatrix<T>& cm)//对所有元素取双曲正切，返回新的矩阵
{
    ComplexMatrix<T> temp(cm);
    return temp.tanh();
}
template <class T> ComplexMatrix<T> log(ComplexMatrix<T>& cm)//对所有元素取自然对数，返回新的矩阵
{
    ComplexMatrix<T> temp(cm);
    return temp.log();
}
template <class T> ComplexMatrix<T> log10(ComplexMatrix<T>& cm)//对所有元素取10为底的对数，返回新的矩阵
{
    ComplexMatrix<T> temp(cm);
    return temp.log10();
}
template <class T> ComplexMatrix<T> exp(ComplexMatrix<T>& cm)//对所有矩阵进行以e为底的指数运算，返回新的矩阵
{
    ComplexMatrix<T> temp(cm);
    return temp.exp();
}
template <class T> ComplexMatrix<T> sqrt(ComplexMatrix<T>& cm)//对矩阵所有的元素取二次方根，返回新的矩阵
{
    ComplexMatrix<T> temp(cm);
    return temp.sqrt();
}
template <class T> ComplexMatrix<T>& ComplexMatrix<T>::dotP(const ComplexMatrix<T>& multiplicator)
{
    if((this->m_col-multiplicator.column())*(this->m_row-multiplicator.row())!=0)
        throw Exception("ComplexMatrix","dotP","the dimensions of two matrices are not equal");
    unsigned long row, col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            this->set(row,col,this->get(row,col)*multiplicator.get(row,col));
        }
    return *this;
}
template <class T> ComplexMatrix<T>& ComplexMatrix<T>::dotD(const ComplexMatrix<T>& divisor)
{
    if((this->m_col-divisor.column())*(this->m_row-divisor.row())!=0)
        throw Exception("Matrix","dotD","the dimensions of two matrices are not equal");
    unsigned long row, col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            this->set(row,col,this->get(row,col)/divisor.get(row,col));
        }
    return *this;
}
template <class T> ComplexMatrix<T> dotP(const ComplexMatrix<T>& multiplicand, const ComplexMatrix<T> multiplicator)
{
    ComplexMatrix<T> temp(multiplicand);
    return temp.dotP(multiplicator);
}
template <class T> ComplexMatrix<T> dotD(const ComplexMatrix<T>& dividend,const ComplexMatrix<T>& divisor)
{
    ComplexMatrix<T> temp(dividend);
    return temp.dotD(divisor);
}
}
