#ifndef __LDAS_COMPLEXMATRIX_HPP
#define __LDAS_COMPLEXMATRIX_HPP

#include <iostream>
#include <fstream>
#include <cmath>
#undef max
#undef min
#include <complex>
#include "exception.h"

namespace ldas
{
using namespace std;

template <class T> class Matrix;
template <class T> class ComplexMatrix;

template <class T> ostream& operator <<(ostream& os,const ComplexMatrix<T>& cmatrix);
template <class T> Matrix<T> abs(ComplexMatrix<T>& cm);//返回矩阵每个元素的模（modulus），The modulus of a complex number a + bi is sqrt(a^2 + b^2)
template <class T> Matrix<T> norm(ComplexMatrix<T>& cm);//返回矩阵每一个元素的范数（norm）,The norm of a complex number a + bi is a^2 + b^2
template <class T> ComplexMatrix<T> sin(ComplexMatrix<T>& cm);//对所有的元素取正弦，返回新的矩阵
template <class T> ComplexMatrix<T> cos(ComplexMatrix<T>& cm);//对所有元素取余弦，返回新的矩阵
template <class T> ComplexMatrix<T> tan(ComplexMatrix<T>& cm);//对所有元素取正切，返回新的矩阵
template <class T> ComplexMatrix<T> sinh(ComplexMatrix<T>& cm);//对所有元素取双曲正弦，返回新的矩阵
template <class T> ComplexMatrix<T> cosh(ComplexMatrix<T>& cm);//对所有元素取双曲余弦，返回新的矩阵
template <class T> ComplexMatrix<T> tanh(ComplexMatrix<T>& cm);//对所有元素取双曲正切，返回新的矩阵
template <class T> ComplexMatrix<T> log(ComplexMatrix<T>& cm);//对所有元素取自然对数，返回新的矩阵
template <class T> ComplexMatrix<T> log10(ComplexMatrix<T>& cm);//对所有元素取10为底的对数，返回新的矩阵
template <class T> ComplexMatrix<T> exp(const ComplexMatrix<T>& cm);//对所有矩阵进行以e为底的指数运算，返回新的矩阵
template <class T> ComplexMatrix<T> sqrt(ComplexMatrix<T>& cm);//对矩阵所有的元素取二次方根，返回新的矩阵
///dot product and dot divide
template <class T> ComplexMatrix<T> dotP(const ComplexMatrix<T>& multiplicand, const ComplexMatrix<T>& multiplicator);
template <class T> ComplexMatrix<T> dotD(const ComplexMatrix<T>& dividend, const ComplexMatrix<T>& divisor);

///复数矩阵模板类
template <class T>
class ComplexMatrix
{
public:
    /** \brief 默认构造函数
    *
    * 构造一个空的复数矩阵，注意在使用之前一定要通过SetDim函数或者赋值操作，设定矩阵的维数
    */
    ComplexMatrix();
    /** \brief 构造函数
    *
    * 使用行数和列数构造一个空的矩阵，实部和虚部都为零的复数矩阵
    */
    ComplexMatrix(const unsigned long row, const unsigned long col);
    ///拷贝构造函数
    ComplexMatrix(const ComplexMatrix& cmatrix);
    /** \brief 直接以文件形式打开矩阵
    *
    * 通过存放在文件中的数据构造矩阵，数据文件有特定的格式
    * 文件最开始部分表示unsigned long大小的行、列数
    * 之后是具体的数据
    * 通用性不强，不同硬件平台有可能字节大小不同，也存在CPU的BIG/LITTLE endian问题
    */
    ComplexMatrix(const char *FileName);
    ///通过两个一维数组和复数矩阵的行数和列数，构造一个复数矩阵
    ComplexMatrix(const T* RealBuf,const T* ImagBuf, const unsigned long row, const unsigned long col);
    ///通过实部矩阵和虚部矩阵构造复数矩阵
    ComplexMatrix(unsigned long row, unsigned long col, Matrix<T> real,Matrix<T> imag);
    ///析构函数，释放内存
    ~ComplexMatrix();

    ///返回矩阵的行数
    unsigned long row() const;
    ///返回矩阵的列数
    unsigned long column() const;

    ///根据下标给矩阵的元素赋值
    void set(const unsigned long row, const unsigned long col, const T real, const T image);
    ///根据下标，给矩阵的元素赋值
    void set(const unsigned long row, const unsigned long col, const complex<T>& c);
    ///根据下标返回矩阵元素的值，返回类型：复数
    complex<T> get(const unsigned long row, const unsigned long col) const;

    ///当前矩阵的实部和虚部都乘以-1，修改当前矩阵,返回当前矩阵
    ComplexMatrix<T>& neg();
    ///对当前矩阵进行转置操作，修改当前矩阵,返回当前矩阵
    ComplexMatrix<T>& tran();
    ///对当前矩阵进行转置，产生一个新的矩阵，不改变原来矩阵，返回新的矩阵
    ComplexMatrix<T> tran() const;
    ///对当前矩阵取共轭操作，返回当前矩阵
    ComplexMatrix<T>& conj();
    ///对当前矩阵共轭，再取转置，返回新的矩阵
    ComplexMatrix<T> conj_tran() const;

    ///返回当前矩阵的实部
    Matrix<T> real() const;
    ///返回当前矩阵的虚部
    Matrix<T> imag() const;
    ///返回当前矩阵所有元素的幅角（argument）,For a complex number a + bi, the argument is equal to arctan(b/a).
    Matrix<T> arg() const;

    /** \brief 返回矩阵每个元素的模（modulus）
    *
    * The modulus of a complex number a + bi is sqrt(a^2 + b^2)
    */
    Matrix<T> abs() const;
    /** \brief 返回矩阵每一个元素的范数（norm）
    *
    * The norm of a complex number a + bi is a^2 + b^2
    */
    Matrix<T> norm() const;
    ///对所有的元素取正弦，修改当前矩阵，返回当前矩阵
    ComplexMatrix<T>& sin();
    ///对所有元素取余弦，修改当前矩阵，返回当前矩阵
    ComplexMatrix<T>& cos();
    ///对所有元素取正切，修改当前矩阵，返回当前矩阵
    ComplexMatrix<T>& tan();
    ///对所有元素取双曲正弦，修改当前矩阵，返回当前矩阵
    ComplexMatrix<T>& sinh();
    ///对所有元素取双曲余弦，修改当前矩阵，返回当前矩阵
    ComplexMatrix<T>& cosh();
    ///对所有元素取双曲正切，修改当前矩阵，返回当前矩阵
    ComplexMatrix<T>& tanh();
    ///对所有元素取自然对数，修改当前矩阵，返回当前矩阵
    ComplexMatrix<T>& log();
    ///对所有元素取10为底的对数，修改当前矩阵，返回当前矩阵
    ComplexMatrix<T>& log10();
    ///对所有矩阵进行以e为底的指数运算，修改当前矩阵，返回当前矩阵
    ComplexMatrix<T>& exp();
    ///对矩阵所有的元素取二次方根，修改当前矩阵，返回当前矩阵
    ComplexMatrix<T>& sqrt();

    ///根据下标，返回一个复数
    complex<T> operator ()(const unsigned long row, const unsigned long col) const;
    ///当前矩阵等于right,返回当前矩阵
    ComplexMatrix<T>& operator =(const ComplexMatrix<T>& rcmatrix);
    ///当前矩阵加上right,返回当前矩阵
    ComplexMatrix<T>& operator +=(const ComplexMatrix<T>& rcmatrix);
    ///当前矩阵加上right,返回一个新的矩阵
    ComplexMatrix<T> operator +(const ComplexMatrix<T>& rcmatrix) const;
    ///矩阵求负，产生新的矩阵，返回新的矩阵
    ComplexMatrix<T> operator -() const;
    ///当前矩阵减去rcmatrix,修改当前矩阵，返回当前矩阵
    ComplexMatrix<T>& operator -=(const ComplexMatrix<T>& rcmatrix);
    ///当前矩阵减去rcmatrix,产生新的矩阵，返回新的矩阵
    ComplexMatrix<T> operator -(const ComplexMatrix<T>& rcmatrix) const;
    ///矩阵数乘k,修改当前矩阵，返回当前矩阵
    ComplexMatrix<T>& operator *=(const T k);
    ///矩阵数乘k，返回新的矩阵
    ComplexMatrix<T> operator *(const T k) const;
    ///矩阵的每一个元素都除以k,修改当前矩阵，返回当前矩阵
    ComplexMatrix<T>& operator /=(const T k);
    ///矩阵的每一个元素都除以k,返回一个新的矩阵
    ComplexMatrix<T> operator /(const T k) const;
    ///当前矩阵乘以矩阵rcmatrix,修改当前矩阵，返回当前矩阵
    ComplexMatrix<T>& operator *=(const ComplexMatrix<T>& rcmatrix);
    ///当前矩阵乘以矩阵rcmatrix，产生新的矩阵，返回新的矩阵
    ComplexMatrix<T> operator *(const ComplexMatrix<T>& rcmatrix) const;
    ///矩阵的每一个元素都进行参数为y^x的幂运算，修改当前矩阵，返回当前矩阵
    ComplexMatrix<T>& operator ^=(const T x);
    ///矩阵的每一个元素都进行按照y^x格式进行幂运算，返回一个新的矩阵
    ComplexMatrix<T> operator ^(const T x) const;
    ///输出符重载
    friend ostream& operator << <>(ostream& os,const ComplexMatrix<T>& cmatrix);
    ///矩阵求逆，返回新的矩阵
    ComplexMatrix<T> operator ~() const;

    ///当前矩阵交换两行
    void swapRows(const unsigned long row1, const unsigned long row2);
    ///当前矩阵交换两列
    void swapCols(const unsigned long col1, const unsigned long col2);
    ///矩阵求逆，修改当前矩阵，返回当前矩阵
    ComplexMatrix<T>& inv();
    ///将当前矩阵设置成为但为矩阵，主元实部为1，虚部为零，其它实部虚部都为零
    void genI();

    ///当前矩阵是否为方阵
    bool is_square() const;
    ///当前矩阵是否为对称矩阵
    bool is_symmetry() const;
    /** \brief 直接以文件形式保存矩阵
    *
    * 数据文件有特定的格式，可以直接用构造函数读取
    * 文件最开始部分表示unsigned long大小的行、列数
    * 之后是具体的数据
    * 通用性不强，不同硬件平台有可能字节大小不同，也存在CPU的BIG/LITTLE endian问题
    */
    void save(const char *FileName);
    ///设定矩阵的行列数
    void setDim(const unsigned long rowNum, const unsigned long colNum);
    ///the dot product of complex matrix
    ComplexMatrix<T>& dotP(const ComplexMatrix<T>& multiplicator);
    ///dot divide operation of complex matrix
    ComplexMatrix<T>& dotD(const ComplexMatrix<T>& divisor);
private:
    unsigned long m_row;//存放矩阵的行数
    unsigned long m_col;//存放矩阵的列数
    T m_Tolerance;//存放矩阵的容许误差，小于等于这个数认为为零
    T* m_RealBuf;//矩阵实部的内存缓冲区
    T* m_ImagBuf;//矩阵虚部的内存缓冲区
    T ValueAbs(T value) const;//计算元素的绝对值
    T MaxAbs(unsigned long& row, unsigned long& rol, unsigned long k);//计算第k行和第k列之后的主元和位置
    T Adjust(T a);//
    void DoAdjust();//应用拉近技术，调正数值
};

}//end namespace

#include "complexmatrix.cpp"
#endif
