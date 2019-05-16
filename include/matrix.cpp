namespace ldas
{
template <class T> Matrix<T>::Matrix()
{
    m_col =0;
    m_row =0;
    m_Tolerance =(T)(1.0e-8);
    m_buf=NULL;
}
template <class T> Matrix<T>::Matrix(const unsigned long row, const unsigned long col)
{
    m_row=row;
    m_col=col;
    m_Tolerance=(T)(1.0e-8);
    unsigned long k;
    m_buf=new T[m_row*m_col];
    for(k=0; k<m_row*m_col; k++)
    {
        m_buf[k]=(T)0.0;
    }
}

template <class T> Matrix<T>::Matrix(const Matrix<T>& matrix)//拷贝构造函数
{
    unsigned long i,j;
    m_row=matrix.row();
    m_col=matrix.column();
    m_Tolerance=(T)(1.0e-8);
    m_buf=new T[m_row*m_col];
    for(i=1; i<=m_row; i++) //将矩阵matrix的元素值复制到当前矩阵内
        for(j=1; j<=m_col; j++)
            this->set(i,j,matrix.get(i,j));
}

template <class T> Matrix<T>::Matrix(T* const buf,const unsigned long row,const unsigned long col)
{
    unsigned long k;
    m_row=row;
    m_col=col;
    m_Tolerance=(T)(1.0e-8);
    m_buf=new T[row*col];
    for(k=0; k<m_row*m_col; k++)
        m_buf[k]=buf[k];
}
template <class T> Matrix<T>::Matrix(T** const buf,const unsigned long row,const unsigned long col)
{
    unsigned long i,j;
    m_row=row;
    m_col=col;
    m_Tolerance=(T)(1.0e-8);
    m_buf=new T[m_row*m_col];
    for(i=0; i<m_row; i++)
        for(j=0; j<m_col; j++)
            m_buf[i*m_col+j]=buf[i][j];
}
template <class T> Matrix<T>::Matrix(const char* FileName)
{
    ifstream in;
    in.open(FileName) ;
    if(!in)//如果文件不存在或者不能打开，抛出异常
        throw Exception("Matrix","Matrix(const char* FileName)","Can not open data file");
    in>>this->m_row;//读入行数
    in>>this->m_col;//读入列数
    m_Tolerance=(T)(1.0e-8);
    m_buf=new T[m_row*m_col];
    T temp;
    unsigned long row,col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            in>>temp;
            this->set(row,col,temp);
        }
    in.close();
}

template <class T> void Matrix<T>::setDim(const unsigned long rowNum,const unsigned long colNum)
{
    this->m_col=colNum;
    this->m_row=rowNum;
    this->m_buf =new T[rowNum*colNum];
    unsigned long k;
    for(k=0; k<m_row*m_col; k++)
        m_buf[k]=(T)0.0;
}
template <class T> Matrix<T>& Matrix<T>::exp()
{
    unsigned long row,col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
            this->set(row,col,std::exp(this->get(row,col)));
    return (*this);
}
template <class T> Matrix<T>& Matrix<T>::sin()//对所有的元素取正弦，修改当前矩阵，返回当前矩阵
{
    unsigned long row,col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
            this->set(row,col,sin(this->get(row,col)));
    return (*this);
}
template <class T> Matrix<T>& Matrix<T>::cos()//对所有元素取余弦，修改当前矩阵，返回当前矩阵
{
    unsigned long row,col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
            this->set(row,col,std::cos(this->get(row,col)));
    return (*this);
}
template <class T> Matrix<T>& Matrix<T>::tan()//对所有元素取正切，修改当前矩阵，返回当前矩阵
{
    unsigned long row,col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
            this->set(row,col,tan(this->get(row,col)));
    return (*this);
}
template <class T> Matrix<T>& Matrix<T>::sinh()//对所有元素取双曲正弦，修改当前矩阵，返回当前矩阵
{
    unsigned long row,col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
            this->set(row,col,sinh(this->get(row,col)));
    return (*this);
}
template <class T> Matrix<T>& Matrix<T>::cosh()//对所有元素取双曲余弦，修改当前矩阵，返回当前矩阵
{
    unsigned long row,col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
            this->set(row,col,cosh(this->get(row,col)));
    return (*this);
}
template <class T> Matrix<T>& Matrix<T>::tanh()//对所有元素取双曲正切，修改当前矩阵，返回当前矩阵
{
    unsigned long row,col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
            this->set(row,col,sinh(this->get(row,col)));
    return (*this);
}
template <class T> Matrix<T>& Matrix<T>::log()//对所有元素取自然对数，修改当前矩阵，返回当前矩阵
{
    unsigned long row,col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            if(this->get(row,col)<=0)
                throw Exception("Matrix","log","one or more element of matrix is less than or equal to zero");
        }
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
            this->set(row,col,log(this->get(row,col)));
    return (*this);
}
template <class T> Matrix<T>& Matrix<T>::log10()//对所有元素取10为底的对数，修改当前矩阵，返回当前矩阵
{
    unsigned long row,col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            if(this->get(row,col)<=0)
                throw Exception("Matrix","log10","one or more element of matrix is less than or equal to zero");
        }
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
            this->set(row,col,log10(this->get(row,col)));
    return (*this);
}
template <class T> Matrix<T>& Matrix<T>::sqrt()
{
    unsigned long row,col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            if(this->get(row,col)<0)
                throw Exception("Matrix","sqrt","one or more element of matrix is less than zero");
        }
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
            this->set(row,col,std::sqrt(this->get(row,col)));
    return (*this);
}
template <class T> Matrix<T>& Matrix<T>::asin()
{
    unsigned long row,col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            this->set(row,col,(T)asin((double)this->get(row,col)));
        }
    return *this;
}
template <class T> Matrix<T>& Matrix<T>::acos()
{
    unsigned long row,col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            this->set(row,col,(T)std::acos((double)this->get(row,col)));
        }
    return *this;
}
template <class T> Matrix<T>& Matrix<T>::atan()
{
    unsigned long row,col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            this->set(row,col,(T)atan((double)this->get(row,col)));
        }
    return *this;
}
template <class T> void Matrix<T>::save(const char* FileName) const
{
    unsigned long row, col;
    ofstream out;
    out.open(FileName);
    if(!out)
        throw Exception("Matrix","save","Can not create data file");
    out<<m_row<<" "<<m_col<<endl;
    for(row=1; row<=m_row; row++)
    {
        for(col=1; col<=m_col; col++)
        {
            //out.setf(ios_base::showpoint);
            out.width(15);
            out<<get(row,col)<<" ";
        }
        out<<endl;
    }
    out.close();
}
template <class T> Matrix<T>::~Matrix()
{
    delete[] this->m_buf;//释放内存空间
}
template <class T> T Matrix<T>::Adjust(T a) const//对小数点后五位四舍五入
{
    T b=floor(a*10000.0+0.5)/10000.0;
    if(ValueAbs(b-a)<m_Tolerance)//如果舍入误差很小,小于容许误差，则抛弃
        a=b;
    return a;
}
template <class T> void Matrix<T>::DoAdjust()
{
    unsigned long row,col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            set(row,col,Adjust((*this)(row,col)));
        }
}
template <class T> unsigned long Matrix<T>::row() const//返回矩阵的行数
{
    return m_row;
}
template <class T> unsigned long Matrix<T>::column() const//返回矩阵的列数
{
    return m_col;
}
template <class T> bool Matrix<T>::is_equal(const Matrix<T>& matrix) const//判断两个矩阵是否相等
{
    // is this matrix equal to another matrix
    if( (m_col==matrix.column()) && (m_row==matrix.row()) )//如果两个矩阵的行数列数相等,比较各个元素
    {
        long i,j;
        for(i=1; i<=m_row; i++)
            for(j=1; j<=m_col; j++)
            {
                if (abs(this->get(i,j)-matrix.get(i,j)) > this->m_Tolerance)
                    return false;
            }
        return true;
    }
    else
        return false;
}
template <class T> bool Matrix<T>::is_square() const//判断当前矩阵是否为方阵
{
    //is a square matrix?
    if(m_row==m_col)
        return true;
    else
        return false;
}

template <class T> bool Matrix<T>::is_symmetry() const//判断当前矩阵是否为对称矩阵, 只对方阵
{
    if(this->m_col!=this->m_row)//如果当前矩阵不是方阵，则抛出异常，只有方阵才满足这项计算的要求
    {
        throw Exception("Matrix","is_symmetry","The matrix is not square.");
        return false;
    }

    unsigned long row, col;
    for(row=1; row<=this->m_row; row++ )
    {
        for(col=(row+1); col<=this->m_col; col++)
        {
            if(ValueAbs(this->get(row,col)-this->get(col,row))>this->m_Tolerance)
                return false;
        }
    }
    return true;
}


template <class T> bool Matrix<T>::is_I() const
{
    if(!this->is_square())//如果当前矩阵不是方阵，抛出异常
        throw Exception("Matrix","is_I","The matrix is not square.");
    else
    {
        unsigned long row,col;
        for(row=1; row<=this->m_row; row++)
            for(col=1; col<=this->m_col; col++)
            {
                if(row==col)
                {
                    T temp, unit;
                    unit=1;
                    temp=*this(row,col)-unit;
                    if(temp<0)// less than zero
                    {
                        if(temp<(-this->m_Tolerance))
                            return false;
                    }
                    else // more than zero
                    {
                        if(temp>this->m_Tolerance)
                            return false;
                    }
                }
                else
                {
                    if(*this(row,col)>0)//more than zero
                    {
                        if(*this(row,col)>this->m_Tolerance)
                            return false;
                    }
                    else //less than zero
                    {
                        if(*this(row,col)<(-this->m_Tolerance))
                            return false;
                    }
                }
            }
        return true;
    }

}

template <class T>
Matrix<T> Matrix<T>::extractRowV(const unsigned long row) const
{
    Matrix<T> matrix(1, this->m_col);	// row vector

    if(row>this->m_row)	// exception manipulation
        throw Exception("Matrix","extractRowV","Subscripts are over matrix size");
    else
    {
        for(int i=1; i<=matrix.column(); i++)
        {
            matrix.set(1, i, this->get(row, i));
        }
        return matrix;
    }
}

template <class T>
Matrix<T> Matrix<T>::extractColV(const unsigned long col) const
{
    Matrix<T> matrix(this->m_row, 1);	// column vector

    if(col>this->m_col)	// exception manipulation
        throw Exception("Matrix","extractColV","Subscripts are over matrix size");
    else
    {
        for(int i=1; i<=(int)matrix.row(); i++)
        {
            matrix.set(i, 1, this->get(i, col));
        }
        return matrix;
    }
}
//将两个矩阵连接，矩阵行数不变，列数为（COL1+COL2)
template <class T>
Matrix<T> Matrix<T>::appendCol(Matrix<T>& m1,Matrix<T>& m2) const
{
    int col1,col2,row1,row2;
    col1 = m1.column();
    col2 = m2.column();
    row1 = m1.row();
    row2 = m2.row();
    Matrix<T> matrix;
    matrix.SetDim(row1,col1+col2);
    if(row1!=row2)
        throw Exception("Matrix","appendCol","matrix row number are not equal.");
    else
    {
        for(int i=1; i<=row1; i++)
        {
            for(int j=1; j<=col1; j++)
                matrix.set(i,j,m1.get(i,j));
            for(int k=1; k<=col2; k++)
                matrix.set(i,col1+k,m2.get(i,k));
        }
        *this=matrix;
        return *this;
    }
}
//将两个矩阵连接，矩阵列数不变，行数为（row1+row2)
template <class T>
Matrix<T> Matrix<T>::appendRow(Matrix<T>& m1,Matrix<T>& m2) const
{
    int col1,col2,row1,row2;
    col1 = m1.column();
    col2 = m2.column();
    row1 = m1.row();
    row2 = m2.row();
    Matrix<T> matrix;
    matrix.SetDim(row1+row2,col1);
    if(col1!=col2)
        throw Exception("Matrix","appendRow","matrix col number are not equal.");
    else
    {
        for(int i=1; i<=col1; i++)
        {
            for(int j=1; j<=row1; j++)
                matrix.set(j,i,m1.get(j,i));
            for(int k=1; k<=row2; k++)
                matrix.set(row1+k,i,m2.get(k,i));
        }
        *this=matrix;
        return *this;
    }
}
//将三个矩阵连接，矩阵行数不变，列数为（COL1+COL2＋COL3)
template <class T>
Matrix<T> Matrix<T>::appendCol(Matrix<T>& m1,Matrix<T>& m2, Matrix<T> & m3) const
{
    Matrix<T> matrix;
    matrix = appendCol(m1,m2);
    matrix = appendCol(matrix,m3);
    *this=matrix;
    return *this;
}
//将三个矩阵连接，矩阵列数不变，行数为（row1+row2＋row3)
template <class T>
Matrix<T> Matrix<T>::appendRow(Matrix<T>& m1,Matrix<T>& m2, Matrix<T> & m3) const
{
    Matrix<T> matrix;
    matrix = appendRow(m1,m2);
    matrix = appendRow(matrix,m3);
    *this=matrix;
    return *this;
}
//复制矩阵
template <class T>
Matrix<T> Matrix<T>::repmat(Matrix<T> &m1, int Row, int Col) const
{

    Matrix<T> matrix(m1);
    for(int i=1; i<Col; i++)
        matrix = appendCol(matrix,m1);

    Matrix<T> tmp(matrix);
    for(int j=1; j<Row; j++)
        matrix = appendRow(matrix,tmp);

    *this=matrix;
    return *this;
}

// 计算矩阵列的和
template <class T>
Matrix<T> Matrix<T>::sumCol(const Matrix<T> &m1) const
{
    int col,row;
    row = m1.row();
    col = m1.column();
    Matrix<T> sum(row,1);
    sum.genZeros();
    for(int i=1; i<=col; i++)
    {
        sum = sum + m1.extractColV(i);
    }
    *this = sum;
    return *this;
}

// 计算矩阵行的和
template <class T>
Matrix<T> Matrix<T>::sumRow(Matrix<T> &m1) const
{
    int col,row;
    row = m1.row();
    col = m1.column();
    Matrix<T> sum(1,col);
    sum.genZeros();
    for(int i=1; i<=row; i++)
    {
        sum = sum + m1.extractRowV(i);
    }
    *this = sum;
    return *this;
}

// 计算矩阵列平均值
template <class T>
Matrix<T> Matrix<T>::meanCol(const Matrix<T> &m1) const
{

    Matrix<T> mean;
    mean.sumCol(m1);
    mean = mean/m1.column();
    *this = mean;
    return *this;
}
// 计算带权重的矩阵列平均值
template <class T>
Matrix<T> Matrix<T>::meanCol(Matrix<T> &m1, Matrix<T> &weight) const
{
    int col,row;
    row = m1.row();
    col = m1.column();
    Matrix<T> mean(row,1);
    mean.genZeros();
    for(int i=1; i<=col; i++)
    {
        mean = mean + m1.extractColV(i)*weight.get(1,i);
    }
    *this = mean;
    return *this;
}

// 计算矩阵行平均值
template <class T>
Matrix<T> Matrix<T>::meanRow(Matrix<T> &m1) const
{
    Matrix<T> mean;
    mean.sumRow(m1);
    mean = mean/m1.row();
    *this = mean;
    return *this;
}

//计算两个矩阵的协方差矩阵，认为矩阵的每一个列向量为一个随机变量集，返回新的矩阵。
template <class T>
Matrix<T> Matrix<T>::cov2(Matrix<T>& m1,Matrix<T>& m2) const
{
    Matrix<T> mean1,mean2, x1,x2,tmp1,tmp2;

    int col1,col2,row1,row2;
    col1 = m1.column();
    col2 = m2.column();
    row1 = m1.row();
    row2 = m2.row();

    mean1 = meanCol(m1);
    mean2 = meanCol(m2);

    Matrix<double> covariance(row1,row2);
    covariance.genZeros();

    for(int i=1; i<=col1; i++)
    {
        tmp1 = m1.extractColV(i)-mean1;
        tmp2 = m2.extractColV(i)-mean2;
        covariance = covariance +tmp1*tmp2.t();
    }
    covariance = covariance/col1;
    *this = covariance;
    return *this;
}
//计算带权重的两个矩阵的协方差矩阵，认为矩阵的每一个列向量为一个随机变量集，返回新的矩阵。
template <class T>
Matrix<T> Matrix<T>::cov2(Matrix<T>& m1,Matrix<T>& m2,Matrix<T> & weight) const
{
    Matrix<T> mean1,mean2, x1,x2,tmp1,tmp2;

    int col1,col2,row1,row2;
    col1 = m1.column();
    col2 = m2.column();
    row1 = m1.row();
    row2 = m2.row();

    mean1 = meanCol(m1);
    mean2 = meanCol(m2);
    Matrix<double> covariance(row1,row2);
    covariance.genZeros();
    for(int i=1; i<=col1; i++)
    {
        tmp1 = m1.extractColV(i)-mean1;
        tmp2 = m2.extractColV(i)-mean2;
        covariance = covariance +weight.get(1,i)*tmp1*tmp2.t();
    }
    *this = covariance;
    return *this;
}
//计算两个矩阵的方差矩阵，认为矩阵的每一个列向量为一个随机变量集，返回新的矩阵。
template <class T>
Matrix<T> Matrix<T>::var2(Matrix<T>& m1,Matrix<T>& m2) const
{
    Matrix<T> tmp;
    double value;

    int col1,col2,row1,row2;
    col1 = m1.column();
    col2 = m2.column();
    row1 = m1.row();
    row2 = m2.row();



    Matrix<double> variance(row1,1);
    variance.genZeros();

    if(ValueAbs((row1-row2)*(col1-col2))>this->m_Tolerance)	// exception manipulation
        throw Exception("Matrix","var2","Subscripts are over matrix size");
    for(int i=1; i<=col1; i++)
    {
        tmp = m1.extractColV(i)-m2.extractColV(i);
        for(int j=1; j<=row1; j++)
        {
            value = tmp.get(j,1);
            tmp.set(j,1,value*value);
        }
        variance = variance + tmp;
    }
    variance = variance/col1;
    *this = variance;
    return *this;
}


template <class T>
void Matrix<T>::setRowV(const unsigned long row, const Matrix<T>& rowElement)
{
    if(row>m_row)	// exception manipulation
        throw Exception("Matrix","setRowV","Subscripts are over matrix size");
    else
    {
        unsigned long i;
        for(i=1; i<=m_col; i++)
        {
            set(row, i, rowElement.get(1, i));
        }
    }
}

template <class T>
void Matrix<T>::setColV(const unsigned long col, const Matrix<T>& colElement)
{
    if(col>m_col)	// exception manipulation
        throw Exception("Matrix","setColV","Subscripts are over matrix size");
    else
    {
        unsigned long i;
        for(i=1; i<=m_row; i++)
        {
            set(i, col, colElement.get(i, 1));
        }
    }
}

template <class T> T Matrix<T>::get(const unsigned long row, const unsigned long col) const
{
    //这里要注意异常处理，可能存在要求返回的矩阵元素的下标超过矩阵的规模
    if(row>this->m_row || col>this->m_col || row<1 || col <1)//如果下标超出矩阵的范围，抛出异常
        throw Exception("Matrix","get","Subscripts are over matrix size");
    else
    {
        return m_buf[(row-1)*m_col+col-1];//从缓冲区读出数值 (row-1)*m_col+col-1
    }
}

template <class T> void Matrix<T>::set(unsigned long row,unsigned long col,const T value)
{
    if(row>this->m_row|| col > this->m_col || row<1 || col <1)//如果下标超过矩阵的规模
        throw Exception("Matrix","set","Subscripts are over matrix size");
    else
    {
        this->m_buf[(row-1)*m_col+col-1]=value;
    }
}
/*
template <class T> void Matrix<T>::tran()
{
    Matrix<T> matrix(this->m_col,this->m_row);// 一个和当前矩阵规模相同的矩阵
    unsigned long row, col;//矩阵的行数和列数
    for(row=1; row<=matrix.row(); row++)
        for(col=1; col<=matrix.column(); col++)
        {
            matrix.set(row,col,get(col,row));
        }
    *this=matrix;
    return *this;
}*/

template <class T>
Matrix<T> Matrix<T>::tran() const
{
    Matrix<T> matrix(this->m_col,this->m_row);// 一个和当前矩阵规模相同的矩阵
    unsigned long row, col;//矩阵的行数和列数
    for(row=1; row<=matrix.row(); row++)
        for(col=1; col<=matrix.column(); col++)
        {
            matrix.set(row,col,get(col,row));
        }
    return matrix;
}

template <class T> Matrix<T>& Matrix<T>::neg()
{
    unsigned long row,col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            this->set(row,col,(-(*this)(row,col)));
        }
    return *this;
}
template <class T> void Matrix<T>::swapCols(const unsigned long col1,const unsigned long col2)
{
    if(col1>this->m_col||col2>this->m_col)//如果列数超过了矩阵的规模，抛出异常。
        throw Exception("Matrix","swapCols","At lease one column number is over the matrix size");
    else
    {
        T temp;
        unsigned long row;
        for(row=1; row<=this->m_row; row++)
        {
            temp=get(row,col1);
            set(row,col1,get(row,col2));
            set(row,col2,temp);
        }
    }
}
template <class T> void Matrix<T>::swapRows(const unsigned long row1,const unsigned long row2)
{
    if(row1>m_row||row2>m_col)//如果行数超过了矩阵的规模，抛出异常
        throw Exception("Matrix","swapRows","At least one row number is over the matrix size");
    else
    {
        unsigned long col;
        T temp;
        for(col=1; col<=this->m_col; col++)
        {
            temp=get(row1,col);
            set(row1,col,get(row2,col));
            set(row2,col,temp);
        }
    }
}
template <class T> T Matrix<T>::operator()(const unsigned long row, const unsigned long col) const
{
    if(row>this->m_row|| col>this->m_col)//如果下标超过矩阵规模
        throw Exception("Matrix","Operator()","Subscripts are over the matrix size");
    else
    {
        return this->m_buf[(row-1)*m_col+col-1];
    }
}

template <class T> Matrix<T>& Matrix<T>::operator =(const Matrix<T>& rmatrix)
{
    m_col=rmatrix.column();
    m_row=rmatrix.row();
    unsigned long row, col;
    delete[] m_buf;
    m_buf = new T[m_col*m_row];
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            set(row,col,(rmatrix.get(row,col)));
        }
    return(*this);
}



template <class T>
Matrix<T>& Matrix<T>::operator =(const T& scalar)
{
    unsigned long row, col;
    for(row=1; row<=this->m_row; row++)
        for(col=1; col<=this->m_col; col++)
        {
            set(row,col,scalar);
        }
    return(*this);
}

template <class T> T Matrix<T>::ValueAbs(const T value) const//计算绝对值，私有函数
{
    if(value>=0)
        return value;
    else
        return -value;
}

template <class T> void Matrix<T>::genI()//产生单位矩阵
{
    unsigned long row, col;
    if(!this->is_square())//如果不是方阵，不能生成单位矩阵，抛出异常
        throw Exception("Matrix","genI","The current matrix is not square, can't generate identity matrix.");
    else
    {
        for(row=1; row<=this->m_row; row++)
            for(col=1; col<=this->m_col; col++)
            {
                if(row==col)
                    this->set(row,col,1.0);
                else
                    this->set(row,col,0.0);
            }
    }
}
template <class T> void Matrix<T>::genOnes()//产生所有元素都为1的矩阵
{
    unsigned long row, col;
    for(row=1; row<=this->m_row; row++ )
        for(col=1; col<=this->m_col; col++ )
        {
            this->set(row,col,1.0);
        }
}
template <class T> void Matrix<T>::genZeros()//产生所有元素都为零的矩阵
{
    unsigned long row, col;
    for(row=1; row<= this->m_row; row++ )
        for(col=1; col<= this->m_col; col++)
        {
            set(row,col,0.0);
        }
}
template <class T> T Matrix<T>::RowMaxAbs(unsigned long& max_row,unsigned long cur_col, unsigned long start_row)//用在计算逆矩阵的函数中
{
    //私有函数
    unsigned long i;
    T temp;
    temp=get(start_row,cur_col);//  *this(start_row,cur_col);//初始化变量
    max_row=start_row;
    for(i=start_row; i<=m_row; i++)
    {
        if(ValueAbs(temp)<=ValueAbs(get(i,cur_col)))
        {
            temp=get(i,cur_col);
            max_row=i;
        }
    }
    return temp;
}
template <class T> T Matrix<T>::MaxAbs(unsigned long &row,unsigned long &col,unsigned long k) const//计算的第k行和第k列之后的主元的位置和值，私有函数
{
    //私有函数
    T temp=0.0;
    unsigned long i,j;
    for(i=k; i<=m_row; i++)
        for(j=k; j<=m_col; j++)
        {
            if(temp<ValueAbs(get(i,j)))//如果有比零大的元素，临时变量就等于这个元素
            {
                row=i;
                col=j;
                temp=ValueAbs(get(i,j));
            }
        }
    return temp;
}
template <class T> T Matrix<T>::DetAndRank(unsigned long &rank) const
{
    if(this->m_col!=this->m_row||this->m_col==0)//如果这个矩阵不是方阵或者矩阵为空，抛出异常
        throw Exception("Matrix","DetAndRank","The matrix is not square or its dimension is zero.");
    else
    {
        unsigned long i,j,k,is,js;//其中is和js用来存放选中主元的位置
        T f,detv,q,d;
        Matrix<T> temp(*this);
        rank=0;//矩阵的秩初值设置为零
        f=1.0;
        detv=1.0;
        for(k=1; k<m_row; k++)
        {
            q=MaxAbs(is,js,k);//q是k行k列之后的主元的数值，is和js是选中主元的位置。
            if(q<=m_Tolerance)//如果整个矩阵最大主元是零，行列式的值
                return 0.0;
            rank++;//
            if(is!=k)//如果主元不在当前行，交换行
            {
                f=-f;
                temp.swapRows(k,is);
            }
            if(js!=k)//如果主元不在当前列，交换列
            {
                f=-f;
                temp.swapCols(k,js);
            }
            q=temp.get(k,k);
            detv*=q;
            for(i=k+1; i<=m_row; i++)
            {
                d=temp.get(i,k)/q;
                for(j=k+1; j<=m_row; j++)
                {
                    temp.set(i,j,(temp.get(i,j)-d*temp.get(k,j)));
                }
            }
        }//k循环结束了，呵呵
        q=temp.get(m_row,m_col);
        if(ValueAbs(q)>m_Tolerance)
        {
            rank++;
            return Adjust(f*detv*q);
        }
        else
        {
            return 0.0;
        }
    }
}
template <class T> T Matrix<T>::det() const//计算矩阵的行列式的值
{
    T temp;
    unsigned long rank;
    temp=DetAndRank(rank);
    return Adjust(temp);
}
template <class T> unsigned long Matrix<T>::rank() const//计算矩阵的秩
{
    unsigned long rank;
    this->DetAndRank(rank);
    return rank;
}

template <class T>
void Matrix<T>::lu(Matrix<T>& L,Matrix<T>& U)
{
    //注意，输入参数L,U的大小应该和当前矩阵相同
    //unsigned long k,i,j,jj;
    T det;
    /*T temp;*/
    Matrix<T> matrix(*this);
    if(!this->is_square())//如果当前矩阵不是方阵，抛出异常，程序结束
        throw Exception("Matrix","lu","the matrix is not square.");
    det=matrix.det();//计算当前矩阵的行列式的数值，如果行列式的数值为零，那么退出程序
    if(ValueAbs(det)<this->m_Tolerance)//如果行列式的值等于零，抛出异常
        throw Exception("Matrix","lu","the matrix is singular.");
    L=*this;
    U=*this;
    L.genZeros();
    U.genZeros();
    unsigned long i,r,k;
    T temp;
    for(i=1; i<=this->m_row; i++)
    {
        U.set(1,i,this->get(1,i));
        //L.set(i,1,this->get(i,1)/U.get(1,1));
    }
    L.set(1,1,1);
    for(i=2; i<=this->m_row; i++)
    {
        //U.set(1,i,this->get(1,i));
        L.set(i,1,this->get(i,1)/U.get(1,1));
        L.set(i,i,1);
    }
    for(r=2; r<=this->m_col; r++)
    {
        for(i=r; i<=this->m_col; i++)
        {
            temp=0;
            for(k=1; k<=r-1; k++)
                temp+=L(r,k)*U(k,i);
            U.set(r,i,(this->get(r,i)-temp));
        }
        for(i=r+1; i<=m_col; i++)
        {
            if(r!=m_row)
            {
                temp=0;
                for(k=1; k<=(r-1); k++)
                    temp+=L(i,k)*U(k,r);
                //cout<<U(r,r)<<endl;
                L.set(i,r,(this->get(i,r)-temp)/U(r,r));
            }
        }
    }
}

template <class T>
Matrix<T> Matrix<T>::chold() const
{
    if ( !this->is_symmetry() )
        throw Exception("Matrix","Chold","the matrix is not symmetric matrix");

    Matrix<T> matrix(*this);
    int n = this->row();
    int i, j;
    T s;

    for (i=1; i<=n; i++)
    {
        for (j=i; j<=n; j++)
        {
            s = this->get(i,j);
            for (int m=i-1; m>=1; m--)
                s -= matrix.get(i,m)*matrix.get(j,m);
            if (i==j)
            {
                if (s<m_Tolerance)
                    //throw Exception("Matrix","Chold","not positive definite");
                    matrix.set(i,i,m_Tolerance);

                else
                    matrix.set(i, i, std::sqrt(s));
            }
            else
                matrix.set(j, i, s/matrix.get(i,i));
        }
    }

    for (j=2; j<=n; j++)
        for (i=1; i<j; i++)
            matrix.set(i, j, T(0));

    return matrix;
}

template <class T> Matrix<T>& Matrix<T>::inv()//应用高斯-约当列主元消去法，计算当前矩阵的逆矩阵，
{

    if(!this->is_square())//如果当前矩阵不是方阵，则抛出异常，不能进行求逆运算
        throw Exception("Matrix","Inv","The current matrix is not square, can't calculate inverse.");
    if(ValueAbs(det())<=this->m_Tolerance)//如果矩阵是奇异矩阵，抛出异常，不能进行求逆运算
        throw Exception("Matrix","Inv","the current matrix is singular, can't calculate inverse");
    else
    {
        unsigned long i,j,k,dim;
        T det;//记录矩阵的行列式的值
        dim=m_row+1;
        unsigned long* const Ip=new unsigned long[dim];//用来存放主元位置
        //Ip=new[dim];//初始化，可以作为数组使用

        Matrix<T> matrix(*this);//临时矩阵
        Matrix<T> m(*this);
        //matrix=*this;//将当前矩阵设置到临时矩阵中
        //m=*this;
        m.genZeros();
        det=1.0;//(1)
        for(k=1; k<=m_col; k++)
        {
            //(2)按列选主元
            T c0,h;
            unsigned long max_row;
            c0=matrix.RowMaxAbs(max_row,k,k);
            //	cout<<c0<<" "<<max_row<<endl;
            Ip[k]=max_row;//记录
            if(c0<m_Tolerance)//(3)奇异矩阵，抛出异常
                Exception("Matrix","Inv","the current matrix is singular, can't calculate inverse");
            if(max_row!=k)//(4)如果最大行不等于K,交换行，否则，下一步
            {
                matrix.swapRows(max_row,k);
                det=-det;
            }
            det=det*c0;//(5)
            //(6)
            //matrix. (k,k)=1/c0;
            matrix.set(k,k,1/c0);
            h=matrix(k,k);
            for(i=1; i<=m_row; i++)
            {
                if(i!=k)
                {
                    m.set(i,k,(-matrix(i,k)*h));
                    matrix.set(i,k,m(i,k));
                }
            }
            //end(6)
            //(7)消元计算
            for(i=1; i<=m_row; i++)
                for(j=1; j<=m_col; j++)
                {
                    if((i!=k)&&(j!=k))//不是主元
                    {
                        matrix.set(i,j,(matrix(i,j)+m(i,k)*matrix(k,j)));
                    }
                }
            //end(7)
            //(8)计算主行
            for(j=1; j<=m_col; j++)
            {
                if(j!=k)
                {
                    matrix.set(k,j,(matrix(k,j)*h));
                }
            }
            //(8)end
        }//k循环结束
        unsigned long t;
        for(k=(m_col-1); k>=1; k--)
        {
            t=Ip[k];
            if(k!=t)
            {
                matrix.swapCols(t,k);
            }
        }
        matrix.DoAdjust();
        (*this)=matrix;
        return *this;
    }
}//inv()结束
template <class T> bool Matrix<T>::is_singular() const
{
    T det=this->det() ;
    if(this->ValueAbs(det)<=this->m_Tolerance)
        return true;
    else
        return false;
}
template <class T> Matrix<T> Matrix<T>::cov() const
{
    Matrix<T> cov(m_col,m_col);//矩阵自动初始化为零矩阵
    T* mean= new T[m_col];
    unsigned long row, col,k;
    for(col=1; col<=m_col; col++) //计算每一个随即变量的平均值（期望）
    {
        mean[col-1]=(T)0.0;
        for(row=1; row<=m_row; row++)
            mean[col-1]+=this->get(row,col);
        mean[col-1]/=m_row;
    }
    T temp;
    for(row=1; row<=m_col; row++)
        for(col=row; col<=m_col; col++)
        {
            temp=0;
            for(k=1; k<=m_row; k++)
            {
                temp+=(this->get(k,row)-mean[row-1])*(this->get(k,col)-mean[col-1]);
            }
            cout<<endl;
            temp/=(m_row-1);
            cov.set(row,col,temp);
        }
    for(row=1; row<=cov.row(); row++)
        for(col=1; col<row; col++)
        {
            cov.set(row,col,cov.get(col,row));
        }
    delete[] mean;
    return cov;
}
template <class T> Matrix<T> Matrix<T>::operator ~() const//对当前矩阵进行逆运算，返回一个新矩阵
{
    Matrix<T> temp(*this);
    return (temp.inv());
}
template <class T> Matrix<T>& Matrix<T>::operator ~() //对当前矩阵进行逆运算，返回一个新矩阵
{
    return (this->inv());
}
template <class T> bool Matrix<T>::operator==(const Matrix<T>& rm) const
{
    return is_equal(rm);
}
template <class T> bool Matrix<T>::operator!=(const Matrix<T>& rm) const
{
    return !is_equal(rm);
}
template <class T> Matrix<T>& Matrix<T>::operator -=(const T a)
{
    return (*this+=(-a));
}
template <class T> Matrix<T>& Matrix<T>::operator +=(const T a)//当前矩阵的所有元素都加上常数a
{
    unsigned long row,col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            this->set(row,col,((*this)(row,col)+a));
        }
    return *this;
}
template <class T> Matrix<T> Matrix<T>::operator +(const T a) const//当前矩阵的所有元素都加上常数a，返回一个新的矩阵
{
    Matrix<T> matrix(this->m_row,this->m_col);
    matrix=*this;
    matrix+=a;
    return matrix;
}
template <class T> Matrix<T>& Matrix<T>::operator +=(const Matrix<T>& matrix)//当前矩阵加上矩阵matrix,返回当前矩阵
{
    unsigned long row,col;
    if((m_row!=matrix.row())||(m_col!=matrix.column()))//如果两个矩阵的规模不同，不能进行相加运算，抛出异常
        throw Exception("Matrix","Operator +=","the two matrice are not in same size, can not perform addtion operation.");
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            this->set(row,col,(matrix.get(row,col)+(*this).get(row,col)));
        }
    return *this;
}
template <class T> Matrix<T> Matrix<T>::operator +(const Matrix<T>& matrix) const//当前矩阵加上matrix,返回一个新的矩阵
{
    if(m_row!=matrix.row()||m_col!=matrix.column())
        throw Exception("Matrix","Operator +","the two matrice are not in same size, can not perform addtion operation.");
    Matrix<T> temp(this->m_row,this->m_col);
    temp=*this;
    temp+=matrix;
    return temp;
}
template <class T> Matrix<T> Matrix<T>::operator -(void) const//当前矩阵的所有的元素都取复制负值，返回当前矩阵
{
    Matrix<T> matrix(*this);
    matrix.neg();
    return matrix;
}
template <class T> Matrix<T>& Matrix<T>::operator -=(const Matrix<T>& matrix)//当前矩阵减去矩阵matrix，返回当前矩阵
{
    unsigned long row,col;
    if((m_row!=matrix.row())||(m_col!=matrix.column()))//如果两个矩阵的规模不同，不能进行相加运算，抛出异常
        throw Exception("Matrix","Operator +=","the two matrice are not in same size, can not perform addtion operation.");
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            this->set(row,col,((*this).get(row,col)-matrix.get(row,col)));
        }
    return *this;
}
template <class T> Matrix<T> Matrix<T>::operator -(const T a) const//当前矩阵的所有元素都减去常数a，返回一个新的矩阵
{
    Matrix<T> matrix(*this);
    return(matrix+(-a));
}
template <class T> Matrix<T> Matrix<T>::operator -(const Matrix<T>& matrix) const//当前矩阵减去矩阵matriz,返还一个新的矩阵
{
    Matrix<T> temp(*this);
    if(m_row!=matrix.row()||m_col!=matrix.column())//如果两个矩阵的大小不一样，抛出异常，函数结束
        throw Exception("Matrix","Operator-","the two matices are not in same size.");
    temp-=matrix;
    temp.DoAdjust();
    //cout<<"here"<<endl;
    return temp;
}
template <class T> Matrix<T>& Matrix<T>::operator *=(const T a)//矩阵的每一个元素乘以常数a，返回当前矩阵
{
    unsigned long row, col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            this->set(row,col,((*this)(row,col)*a));
        }
    DoAdjust();//误差调整
    return *this;
}
template <class T> Matrix<T> Matrix<T>::operator *(const T a) const//矩阵数乘常数a，返回一个新的矩阵
{
    Matrix<T> matrix(*this);
    matrix*=a;
    return matrix;
}
template <class T> Matrix<T>& Matrix<T>::operator *=(const Matrix<T>& matrix)//当前矩阵乘以矩阵matrix,返回当前矩阵
{
    if(this->m_col!=matrix.row())//如果矩阵的规模不匹配，抛出异常，程序结束
        throw Exception("Matrix","operator *=","the size of two matrices is mismatch.");
    unsigned k,j,i;
    T t;
    Matrix<T> temp(*this);//临时矩阵，当前矩阵的数值复制到这个矩阵之中
    m_col=matrix.column();//当前矩阵的列数修改成乘数矩阵的列数
    delete[] this->m_buf;
    this->m_buf=new T[this->m_row*this->m_col];
    this->genZeros();//当前矩阵的所有元素的值设置为零
    //n=m_row;
    for(k=1; k<=m_row; k++)
        for(j=1; j<=m_col; j++)
        {
            t=0.0;
            for(i=1; i<=matrix.row(); i++)
                t+=temp.get(k,i)*matrix.get(i,j);
            this->set(k,j,t);
        }
    this->DoAdjust();
    return *this;
}
template <class T> Matrix<T> Matrix<T>::operator *(const Matrix<T>& matrix) const//当前矩阵乘以矩阵matrix,返回一个新的矩阵
{
    if(this->m_col!=matrix.row())//如果矩阵的规模不匹配，抛出异常，程序结束
        throw Exception("Matrix","operator *","the size of two matrices is mismatch.");
    Matrix<T> temp(*this);
    temp*=matrix;
    return temp;
}
template <class T> Matrix<T>& Matrix<T>::operator /=(const T a)
{
    if(ValueAbs(a)<=this->m_Tolerance)//如果a为零，抛出异常，终止程序
        throw Exception("Matrix","operator /=","Divided by zero");
    T c=1/a;
    (*this)*=c;
    return *this;
}
template <class T> Matrix<T> Matrix<T>::operator /(const T a) const
{
    if(ValueAbs(a)<=this->m_Tolerance)//如果a为零，抛出异常，终止程序
        throw Exception("Matrix","operator /","Divided by zero");
    Matrix<T> matrix(*this);
    matrix/=a;
    return matrix;
}
template <class T> Matrix<T>& Matrix<T>::operator ^=(const double x)
{
    unsigned long row, col;
    T temp1,temp2;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            temp1=(*this)(row,col);
            if(temp1>0.0)
                temp2=(T)pow((double)temp1,x);
            else
                temp2 = 0;
            this->set(row,col,temp2);
        }
    return *this;
}
template <class T> Matrix<T> Matrix<T>::operator ^(const double x) const
{
    Matrix<T> matrix(*this);
    matrix^=x;
    return matrix;
}
template <class T> T Matrix<T>::trace() const
{
    if(!this->is_square())//如果矩阵不是方阵，抛出异常，结束
        throw Exception("Matrix","Trace","the matrix is not square, can not calculate trace");
    T temp;
    unsigned long row;
    temp=(T)0.0;
    for(row=1; row<=m_row; row++)
        temp+=this->get(row,row);
    return temp;
}
template <class T> ostream& operator <<(ostream& os,const Matrix<T>& matrix)
{
    unsigned long row, col;
    for(row=1; row<=matrix.m_row; row++)
    {
        for(col=1; col<=matrix.m_col; col++)
        {
            //os.setf(ios::left);
            os.width(15);
            os<<matrix.m_buf[(row-1)*matrix.m_col+col-1]<<"  ";
        }
        os<<endl;
    }
    return os;
}

//*************************************************
//            function template
//*************************************************
//连接两个矩阵，行数不变，列数相加
template <class T> Matrix<T> appendCol(Matrix<T> &m1, Matrix<T> & m2)
{
    Matrix<T> temp;
    return temp.appendCol(m1,m2);
}
//连接两个矩阵，列数不变，行数相加
template <class T> Matrix<T> appendRow(Matrix<T> &m1, Matrix<T> & m2)
{
    Matrix<T> temp;
    return temp.appendRow(m1,m2);
}
//连接三个矩阵，行数不变，列数相加
template <class T> Matrix<T> appendCol(Matrix<T> &m1, Matrix<T> & m2, Matrix<T> & m3)
{
    Matrix<T> temp;
    return temp.appendCol(m1,m2,m3);
}
//连接三个矩阵，列数不变，行数相加
template <class T> Matrix<T> appendRow(Matrix<T> &m1, Matrix<T> & m2, Matrix<T> & m3)
{
    Matrix<T> temp;
    return temp.appendRow(m1,m2,m3);
}
//复制矩阵
template <class T>	Matrix<T> repmat(Matrix<T> &m1, int Row, int Col)
{
    Matrix<T> temp;
    return temp.repmat(m1,Row,Col);
}
// 矩阵列向量之和
template <class T> Matrix<T> sumCol(Matrix<T> & m1)
{
    Matrix<T> temp;
    return temp.sumCol(m1);
}

// 矩阵行向量之和
template <class T> Matrix<T> sumRow(Matrix<T> & m1)
{
    Matrix<T> temp;
    return temp.sumRow(m1);
}
// 矩阵列向量平均值
template <class T> Matrix<T> meanCol(const Matrix<T> & m1)
{
    Matrix<T> temp;
    return temp.meanCol(m1);
}
// 矩阵带权重列向量平均值
template <class T> Matrix<T> meanCol(Matrix<T> & m1, Matrix<T> &weight)
{
    Matrix<T> temp;
    return temp.meanCol(m1,weight);
}
// 矩阵行向量平均值
template <class T> Matrix<T> meanRow(Matrix<T> & m1)
{
    Matrix<T> temp;
    return temp.meanRow(m1);
}
//计算带权重的两个矩阵的协方差矩阵，认为矩阵的每一个列向量为一个随机变量集
template <class T>
Matrix<T> cov2(Matrix<T>& m1,Matrix<T>& m2)
{
    Matrix<T> temp;
    return temp.cov2(m1,m2);
}

//计算两个矩阵的协方差矩阵，认为矩阵的每一个列向量为一个随机变量集
template <class T>
Matrix<T> cov2(Matrix<T>& m1,Matrix<T>& m2,Matrix<T> & weight)
{
    Matrix<T> temp;
    return temp.cov2(m1,m2,weight);
}
//计算两个矩阵方差
template <class T>
Matrix<T> var2(Matrix<T>& m1,Matrix<T>& m2)
{
    Matrix<T> temp;
    return temp.var2(m1,m2);
}

template <class T> Matrix<T> sin(Matrix<T>& m)//对所有的元素取正弦，返回新的矩阵
{
    Matrix<T> temp(m);
    return temp.sin();
}
template <class T> Matrix<T> cos(Matrix<T>& m)//对所有元素取余弦，返回新的矩阵
{
    Matrix<T> temp(m);
    return m.cos();
}
template <class T> Matrix<T> tan(Matrix<T>& m)//对所有元素取正切，返回新的矩阵
{
    Matrix<T> temp(m);
    return temp.tan();
}
template <class T> Matrix<T> sinh(Matrix<T>& m)//对所有元素取双曲正弦，返回新的矩阵
{
    Matrix<T> temp(m);
    return temp.sinh();
}
template <class T> Matrix<T> cosh(Matrix<T>& m)//对所有元素取双曲余弦，返回新的矩阵
{
    Matrix<T> temp(m);
    return temp.cosh();
}
template <class T> Matrix<T> tanh(Matrix<T>& m)//对所有元素取双曲正切，返回新的矩阵
{
    Matrix<T> temp(m);
    return temp.tanh();
}
template <class T> Matrix<T> log(Matrix<T>& m)//对所有元素取自然对数，返回新的矩阵
{
    Matrix<T> temp(m);
    return temp.log();
}
template <class T> Matrix<T> log10(const Matrix<T>& m)//对所有元素取10为底的对数，返回新的矩阵
{
    Matrix<T> temp(m);
    return temp.log10();
}
template <class T> Matrix<T> exp(const Matrix<T>& m)//对所有矩阵进行以e为底的指数运算，返回新的矩阵
{
    Matrix<T> temp(m);
    return temp.exp();
}
template <class T> Matrix<T> sqrt(const Matrix<T>& m)//对矩阵所有的元素取二次方根，返回新的矩阵
{
    Matrix<T> temp(m);
    return temp.sqrt();
}
template <class T> Matrix<T>& Matrix<T>::dotP(const Matrix<T>& multiplicator)
{
    if((this->m_col-multiplicator.column())*(this->m_row-multiplicator.row())!=0)
        throw Exception("Matrix","dotP","the dimensions of two matrices are not equal");
    unsigned long row, col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            this->set(row,col,this->get(row,col)*multiplicator.get(row,col));
        }
    return *this;
}
template <class T> Matrix<T>& Matrix<T>::dotD(const Matrix<T>& divisor)
{
    if((this->m_col-divisor.column())*(this->m_row-divisor.row())!=0)
        throw Exception("Matrix","dotD","the dimensions of two matrices are not equal");
    unsigned long row, col;
    for(row=1; row<=m_row; row++)
        for(col=1; col<=m_col; col++)
        {
            if(this->ValueAbs(divisor.get(row,col))<=this->m_Tolerance)
                throw Exception("Matrix","dotD","Divided by zero");
            else
                this->set(row,col,this->get(row,col)/divisor.get(row,col));
        }
    return *this;
}
template <class T> Matrix<T> dotP(const Matrix<T>& multiplicand, const Matrix<T>& multiplicator)
{
    Matrix<T> temp(multiplicand);
    return temp.dotP(multiplicator);
}
template <class T> Matrix<T> dotD(const Matrix<T>& dividend,const Matrix<T>& divisor)
{
    Matrix<T> temp(dividend);
    return temp.dotD(divisor);
}
}
