#ifndef __LDAS_EXCEPTION_H
#define __LDAS_EXCEPTION_H
#include <iostream>

namespace ldas
{
///异常类
class Exception
{
public:
    /** \brief 用类名、方法名和出错信息描述异常
    *
    * \param class_name 类名称
    * \param method_name 方法
    * \param description 描述信息
    */
    Exception(const char* class_name,const char* method_name,const char* description);
    Exception(const std::string class_name,const std::string method_name,const std::string description);
    ~Exception();
    ///输出出错信息
    void what();
private:
    const char* m_info;//存放异常信息
    const char* m_className;//发生异常的类的名字
    const char* m_methodName;//发生异常的方法名字
};
}
#endif // __LDAS_EXCEPTION_H
