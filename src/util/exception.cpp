#include "exception.h"

namespace ldas
{
Exception::Exception(const char* class_name,const char* method_name,const char* description):m_info(description),m_className(class_name),m_methodName(method_name)
{
}

Exception::Exception(const std::string class_name,const std::string method_name,const std::string description):m_info(description.c_str()),m_className(class_name.c_str()),m_methodName(method_name.c_str())
{
}

Exception::~Exception()
{
}
void Exception::what()
{
    std::cout<<"Exception information:"<<std::endl;
    std::cout<<"class name: "<<m_className<<std::endl;
    std::cout<<"method name: "<<m_methodName<<std::endl;
    std::cout<<m_info<<std::endl;
}
}
