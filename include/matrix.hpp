#ifndef __LDAS_MATRIX_HPP
#define __LDAS_MATRIX_HPP

#include <iostream>
#include <fstream>
#include <cmath>
#include "exception.h"

namespace ldas
{
using namespace std;

template <class T> class Matrix;
///模板类友元的预声明，约束模板友元函数
template <class T> std::ostream& operator <<(std::ostream& os,const Matrix<T>& matrix);
///the add,minus,time,divid of a number and a matrix
template <class T> Matrix<T> operator + (const T a, const Matrix<T>& m);
template <class T> Matrix<T> operator - (const T a, const Matrix<T>& m);
template <class T> Matrix<T> operator * (const T a, const Matrix<T>& m);
template <class T> Matrix<T> operator / (const T a, const Matrix<T>& m);
///mathematic operations of real matrix
///对所有的元素取正弦，返回新的矩阵
template <class T> Matrix<T> sin(Matrix<T>& m);
///对所有元素取余弦，返回新的矩阵
template <class T> Matrix<T> cos(Matrix<T>& m);
///对所有元素取正切，返回新的矩阵
template <class T> Matrix<T> tan(Matrix<T>& m);
///对所有元素取双曲正弦，返回新的矩阵
template <class T> Matrix<T> sinh(Matrix<T>& m);
///对所有元素取双曲余弦，返回新的矩阵
template <class T> Matrix<T> cosh(Matrix<T>& m);
///对所有元素取双曲正切，返回新的矩阵
template <class T> Matrix<T> tanh(Matrix<T>& m);
///对所有元素取自然对数，返回新的矩阵
template <class T> Matrix<T> log(Matrix<T>& m);
///对所有元素取10为底的对数，返回新的矩阵
template <class T> Matrix<T> log10(Matrix<T>& m);
///对所有矩阵进行以e为底的指数运算，返回新的矩阵
template <class T> Matrix<T> exp(const Matrix<T>& m);
///对矩阵所有的元素取二次方根，返回新的矩阵
template <class T> Matrix<T> sqrt(const Matrix<T>& m);
///计算矩阵的反正弦函数，返回一个新的矩阵
template <class T> Matrix<T> asin(Matrix<T>& m);
///计算矩阵的反余弦函数，返回一个新的矩阵
template <class T> Matrix<T> acos(Matrix<T>& m);
///计算矩阵的反正切函数，返回一个新的矩阵
template <class T> Matrix<T> atan(Matrix<T>& m);
template <class T> Matrix<T> meanCol(const Matrix<T> & m1);
///dot product and dot divide
template <class T> Matrix<T> dotP(const Matrix<T>& multiplicand, const Matrix<T>& multiplicator);
template <class T> Matrix<T> dotD(const Matrix<T>& dividend, const Matrix<T>& divisor);

///矩阵模板类
template<class T>
class Matrix
{
public:
    /** \brief 默认构造函数
    *
    * 构造一个空的矩阵，注意在使用之前一定要通过SetDim函数或者赋值操作，设定矩阵的维数
    */
    Matrix();
    /** \brief 构造函数
    *
    * 使用行数和列数构造一个空的矩阵，默认值为0
    */
    Matrix(const unsigned long row,const unsigned long col);
    ///拷贝构造函数，生成的矩阵和matrix完全一致
    Matrix(const Matrix<T>& matrix);
    /** \brief 一维数组构造函数
    *
    * 用一维数组初始化矩阵
    * \param buf 一维数组
    * \param row 行数
    * \param col 列数
    */
    Matrix(T* const buf,const unsigned long row,const unsigned long col);
    /// 用一个二维数组初始化
    Matrix(T** const buf,const unsigned long row,const unsigned long col);
    /** \brief 直接以文件形式打开矩阵
    *
    * 通过存放在文件中的数据构造矩阵，数据文件有特定的格式
    * 文件最开始部分表示unsigned long大小的行、列数
    * 之后是具体的数据
    * 通用性不强，不同硬件平台有可能字节大小不同，也存在CPU的BIG/LITTLE endian问题
    */
    Matrix(const char* FileName);
    ///析构函数，清空内存
    ~Matrix();
    ///获得矩阵的行数
    unsigned long row() const;
    ///获得矩阵的列数
    unsigned long column() const;
    ///判断矩阵是否为对称矩阵(symmetry)
    bool is_symmetry() const;
    ///判断矩阵是否为方阵
    bool is_square() const;
    ///判断两个矩阵是否相等
    bool is_equal(const Matrix<T>& matrix) const;
    ///判断矩阵是否是单位矩阵
    bool is_I() const;
    ///反断矩阵是否为奇异矩阵（退化矩阵）
    bool is_singular() const;

    ///提取某一列向量
    Matrix<T> extractColV(const unsigned long col) const;
    ///提取某一行向量
    Matrix<T> extractRowV(const unsigned long row) const;

    ///设置某一列
    void setColV(const unsigned long col, const Matrix<T>& colElement);
    ///设置某一行
    void setRowV(const unsigned long row, const Matrix<T>& rowElement);
    /** \brief 提取某一元素的值
    *
    * \param row 行数 row>=1
    * \param col 列数 col>=1
    */
    T get(const unsigned long row, const unsigned long col) const;
    /** \brief 设置某一元素的值
    *
    * \param row 行数 row>=1
    * \param col 列数 col>=1
    * \param value 赋值
    */
    void set(const unsigned long row, const unsigned long col,const T value);
    ///交换某两行
    void swapRows(const unsigned long row1, const unsigned long row2);
    ///交换某两列
    void swapCols(const unsigned long col1, const unsigned long col2);
    ///将当前矩阵设置成为单位矩阵
    void genI();
    ///将当前矩阵的所有元素都设置成为1
    void genOnes();
    ///将当前矩阵的所有元素的值都是之成为0
    void genZeros();

    /** \brief 求逆矩阵
    *
    * 计算矩阵的逆矩阵,修改当前矩阵，返回当前矩阵
    */
    Matrix<T>& inv();
    /** \brief 求转置矩阵
    *
    * 计算矩阵的转置矩阵，不产生新的矩阵，返回自身
    void tran();*/

    /** \brief 求转置矩阵
    *
    * 计算矩阵的转置矩阵，产生新的矩阵
    */
    Matrix<T> tran() const;
    /** \brief 求负矩阵
    *
    * 计算当前矩阵的负矩阵， 不产生新的矩阵，返回当前矩阵
    */
    Matrix<T>& neg();
    /// 计算矩阵的行列式的值
    T det() const;
    /// 计算矩阵的秩
    unsigned long rank() const;
    ///如果矩阵是方阵，计算矩阵的迹
    T trace() const;
    /** \brief 计算当前矩阵的协方差矩阵
    *
    * 认为当前矩阵的每一个列向量为一个随机变量集，返回新的矩阵
    */
    Matrix<T> cov() const;

    ///返回矩阵元素的值
    T operator ()(const unsigned long row, const unsigned long col) const;
    ///判断两个矩阵是否相等
    bool operator==(const Matrix<T>& r)const;
    ///判断两个矩阵是否不等
    bool operator!=(const Matrix<T>& r)const;
    ///计算逆矩阵，不产生新的矩阵，修改当前矩阵
    Matrix<T>& operator~();
    ///计算逆矩阵，返回新矩阵
    Matrix<T> operator~() const;
    ///从另一个矩阵中赋值
    Matrix<T>& operator=(const Matrix<T>& rmatrix);
    ///给矩阵赋常数
    Matrix<T>& operator=(const T& scalar);
    ///矩阵的每个元素都减去一个常数a
    Matrix<T> operator -(const T a) const;
    ///两个矩阵相减，返回一个新的矩阵
    Matrix<T> operator -(const Matrix<T>& matrix) const;
    ///当前矩阵的每一个元素都减去常数a，返回当前矩阵
    Matrix<T>& operator -=(const T a);
    ///当前矩阵减去矩阵matrix,返回当前矩阵
    Matrix<T>& operator -=(const Matrix<T>& matrix);
    ///计算当前矩阵的负矩阵，产生新的矩阵
    Matrix<T> operator -(void) const;
    ///矩阵的每一个元素都加上常数a，返回一个新的矩阵
    Matrix<T> operator +(const T a) const;
    ///两个矩阵相加，返回一个新的矩阵
    Matrix<T> operator +(const Matrix<T>& matrix) const;
    ///矩阵的每一个元素都加上常数a，返回当前矩阵
    Matrix<T>& operator +=(const T a);
    ///当前矩阵加上矩阵matrix，返回当前矩阵
    Matrix<T>& operator +=(const Matrix<T>& matrix);
    ///当前矩阵乘以常数a，返回一个新的矩阵
    Matrix<T> operator *(const T a) const;
    ///当前矩阵乘以矩阵matrix，返回一个新的矩阵
    Matrix<T> operator *(const Matrix<T>& matrix) const;
    ///当前矩阵乘以常数a,返回当前矩阵
    Matrix<T>& operator *=(const T a);
    ///当前矩阵乘以矩阵matrix,返回当前矩阵
    Matrix<T>& operator *=(const Matrix<T>& matrix);
    ///当前矩阵的所有元素除以常数a,返回一个新的矩阵
    Matrix<T> operator /(const T a) const;
    ///当前矩阵的所有元素都除以常数a,返回当前矩阵
    Matrix<T>& operator /=(const T a);
    ///当前矩阵的每一个元素都取y^x,返回当前矩阵
    Matrix<T>& operator ^=(const double x);
    ///当前矩阵的每一个元素都取y^x，返回新的矩阵
    Matrix<T> operator ^(const double x) const;

    ///输出操作符重载
    friend ostream& operator << <>(ostream& os,const Matrix<T>& matrix);
    ///加法重载
    friend Matrix<T> operator+ (const T a, const Matrix<T>& m)
    {
        Matrix<T> temp(m);
        return (temp+=a);
    }
    ///减法重载
    friend Matrix<T> operator- (const T a, const Matrix<T>& m)
    {
        Matrix<T> temp(m);
        unsigned long row,col;
        for(row=1; row<=m.row(); row++)
            for(col=1; col<=m.column(); col++)
                temp.set(row,col,(a-m.get(row,col)));
        return temp;
    }
    ///乘法重载
    friend Matrix<T> operator* (const T a, const Matrix<T>& m)
    {
        Matrix<T> temp(m);
        return (temp*=a);
    }
    ///除法重载
    friend Matrix<T> operator/ (const T a, const Matrix<T>& m)
    {
        unsigned long row, col;
        for(row=1; row<=m.row(); row++)
            for(col=1; col<=m.column(); col++)
            {
                if(fabs(m.get(row,col))<=m.m_Tolerance)
                    throw Exception("Matrix Template operator","/","divide by zero");
            }
        Matrix<T> temp(m);
        for(row=1; row<=m.row(); row++)
            for(col=1; col<=m.column(); col++)
            {
                temp.set(row,col,a/temp.get(row,col));
            }
        return temp;
    }

    /** \brief 直接以文件形式保存矩阵
    *
    * 数据文件有特定的格式，可以直接用构造函数读取
    * 文件最开始部分表示unsigned long大小的行、列数
    * 之后是具体的数据
    * 通用性不强，不同硬件平台有可能字节大小不同，也存在CPU的BIG/LITTLE endian问题
    */
    void save(const char* FileName) const;
    ///设定维数
    void setDim(const unsigned long rowNum, const unsigned long colNum);

    /// Matrix decomposition, 对矩阵进行三角分解，使用Doolittle算法
    void lu(Matrix<T>& L, Matrix<T>& U);
    /// Cholesky decomposition
    Matrix<T> chold() const;
    ///将两个矩阵连接，矩阵行数不变，列数为（COL1+COL2)
    Matrix<T> appendCol(Matrix<T>& m1,Matrix<T>& m2) const;
    ///将两个矩阵连接，矩阵列数不变，行数为（row1+row2)
    Matrix<T> appendRow(Matrix<T>& m1,Matrix<T>& m2) const;
    ///将三个矩阵连接，矩阵行数不变，列数为（COL1+COL2＋COL3)
    Matrix<T> appendCol(Matrix<T>& m1,Matrix<T>& m2,Matrix<T> & m3) const;
    ///将三个矩阵连接，矩阵列数不变，行数为（row1+row2＋row3)
    Matrix<T> appendRow(Matrix<T>& m1,Matrix<T>& m2,Matrix<T> & m3) const;
    ///计算两个矩阵的协方差矩阵，认为矩阵的每一个列向量为一个随机变量集，返回新的矩阵。
    Matrix<T> cov2(Matrix<T>& m1,Matrix<T>& m2) const;
    ///计算带权重的两个矩阵的协方差矩阵，认为矩阵的每一个列向量为一个随机变量集，返回新的矩阵。
    Matrix<T> cov2(Matrix<T>& m1,Matrix<T>& m2, Matrix<T> & weight) const;
    ///计算带权重的两个矩阵的方差矩阵，认为矩阵的每一个列向量为一个随机变量集，返回新的矩阵。
    Matrix<T> var2(Matrix<T>& m1,Matrix<T>& m2) const;

    ///复制矩阵
    Matrix<T> repmat(Matrix<T> &m1, int Row, int Col) const;
    ///计算矩阵列的和
    Matrix<T> sumCol(const Matrix<T> &m1) const;
    ///计算矩阵行的和
    Matrix<T> sumRow(Matrix<T> &m1) const;
    ///计算矩阵列平均值
    Matrix<T> meanCol(const Matrix<T> &m1) const;
    ///计算带权重的矩阵列向量平均值
    Matrix<T> meanCol(Matrix<T> &m1, Matrix<T> &weight) const;
    ///计算矩阵行平均值
    Matrix<T> meanRow(Matrix<T> &m1) const;

    ///对所有的元素取正弦，修改当前矩阵，返回当前矩阵
    Matrix<T>& sin();
    ///对所有元素取余弦，修改当前矩阵，返回当前矩阵
    Matrix<T>& cos();
    ///对所有元素取正切，修改当前矩阵，返回当前矩阵
    Matrix<T>& tan();
    ///对所有元素取双曲正弦，修改当前矩阵，返回当前矩阵
    Matrix<T>& sinh();
    ///对所有元素取双曲余弦，修改当前矩阵，返回当前矩阵
    Matrix<T>& cosh();
    ///对所有元素取双曲正切，修改当前矩阵，返回当前矩阵
    Matrix<T>& tanh();
    ///对所有元素取自然对数，修改当前矩阵，返回当前矩阵，如果某个矩阵元素小于等于零，抛出异常
    Matrix<T>& log();
    ///对所有元素取10为底的对数，修改当前矩阵，返回当前矩阵，如果某个矩阵元素小于等于零，抛出异常
    Matrix<T>& log10();
    ///对所有矩阵进行以e为底的指数运算，修改当前矩阵，返回当前矩阵
    Matrix<T>& exp();
    ///对所有元素取平方根，修改当前矩阵，返回当前矩阵，如果某个矩阵元素小于零，抛出异常
    Matrix<T>& sqrt();
    ///对当前矩阵所有元素取反正弦，修改当前矩阵，返回当前矩阵
    Matrix<T>& asin();
    ///对当前矩阵所有元素取反余弦，修改当前矩阵，返回当前矩阵
    Matrix<T>& acos();
    ///对当前矩阵所有元素取反正切，修改当前矩阵，返回当前矩阵
    Matrix<T>& atan();
    /** \brief matrix dot procuct operation
    *
    * 两个矩阵点乘，要求两个矩阵行数和列数相同，如果维数不同，抛出异常
    */
    Matrix<T>& dotP(const Matrix<T>& multiplicator);
    /** matrix dot divide operation
    *
    * 两个矩阵点除，要求两个矩阵行数和列数相同，如果维数不同，抛出异常，如果被除数存在零元素，抛出异常
    */
    Matrix<T>& dotD(const Matrix<T>& divisor);//
private:
    T* m_buf;//存放数据的缓冲区
    T m_Tolerance;//数值比较时的容忍误差
    unsigned long m_row;//矩阵的行数
    unsigned long m_col;//矩阵的列数
    T RowMaxAbs(unsigned long& max_row, unsigned long cur_col,unsigned long start_row);//计算第k列的主元的位置和值
    T ValueAbs(const T value)const;//计算绝对值
    T MaxAbs(unsigned long& row, unsigned long& col, unsigned long k) const;//计算第k行和第k列之后的主元的位置和值
    T DetAndRank(unsigned long &rank) const;//计算行列式的值和秩
    T Adjust(T a) const;//对数字小数点后五位进行四舍五入,取近似
    void DoAdjust();//对矩阵元素进行拉近
};
}//end namespace

//#ifdef __inline__
#include "matrix.cpp"
//#endif

#endif
