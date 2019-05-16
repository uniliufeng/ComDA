#include <cpptest.h>
#include "matrix.hpp"
#include "complexmatrix.hpp"

class MatrixTest : public Test::Suite
{
public:
    MatrixTest( )
    {
        TEST_ADD (MatrixTest::MatrixTest1);
        TEST_ADD (MatrixTest::MatrixTest2);
        TEST_ADD (MatrixTest::MatrixTest3);
    }
private:
    void MatrixTest1()
    {
        //default constuction test
        ldas::Matrix<float> m;
        //setDim test
        m.setDim(2,2);
        //value copy test
        m=2;
        TEST_ASSERT(m.get(1,1)==2);
    }
    void MatrixTest2()
    {
        float** f=new float*[5];
        for(int i=0; i<5; i++)
        {
            f[i]=new float[7];
            for(int j=0; j<7; j++)
            {
                f[i][j]=i+j;
            }
        }
        //二维数组初始化测试
        ldas::Matrix<float> m(f,5,7);
        TEST_ASSERT(m(5,7)==10);
        for(int i=0; i<5; i++) delete []f[i];
        delete []f;
        // test open/save to file function
        m.save("test.dat");
        ldas::Matrix<float> m1("test.dat");
        //operator test
        TEST_ASSERT(m==m1);
        m1.set(1,1,m1.get(1,1)+10);
        TEST_ASSERT(m!=m1);
    }
    void MatrixTest3()
    {
        ldas::ComplexMatrix<float> cm;
        std::cout<<cm;
        cm.setDim(3,5);
        cm.set(1,1,2,3);
        //test file open/save function
        cm.save("test.dat");
        ldas::ComplexMatrix<float> cm1("test.dat");
        TEST_ASSERT(cm1.get(1,1)==cm.get(1,1));
    }
};

int main ( )
{
    MatrixTest tests;
    Test::TextOutput output(Test::TextOutput::Verbose);
    return tests.run(output);
}
