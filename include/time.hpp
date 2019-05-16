#ifndef __LDAS_TIME_HPP
#define __LDAS_TIME_HPP

#include <ctime>
#include <iostream>

namespace ldas
{
class TimeSpan;
/** \brief 时间类，用于处理各种时间，默认采用UTC时间 */
class Time
{
private:
    int m_year;
    int m_month;
    int m_day;
    int m_hour;
    int m_minute;
    int m_second;
public:
    Time();
    ///拷贝自另一个时间类
    Time(const Time&);

    /** \brief 通过指定具体的日期构造时间
    *
    * \param year 年份
    * \param month 月份, 1..12
    * \param day 月份中的天, 1..31
    * \param hour 小时, 0..23
    * \param minute 分钟, 0..59
    * \param second 秒钟, 0..59
    * \param zone 整型时区，若需要包行分钟的时区，可使用字符串构建函数
    */
    Time(int year, int month, int day, int hour=0, int minute=0, int second=0,int zone=0);
    /** \brief 根据一个time_t类型的数据，构造一个时间，可以基于当地时间或UTC时间
    *
    * \param time time_t参数
    * \param greenwich 是否采用UTC时间
    */
    Time(time_t time,bool greenwich=true);
    /** \brief 通过字符串构造时间
    * 默认采用ISO 8601标准
    *
    * 如：2010
    * 2010-10
    * 2010-10-13
    * 2010-10-13T12:10
    */
    Time(const std::string& timestr);
    virtual ~Time() {}

    /** \brief Create time with a year, day in the year and seconds in the day.
    *
    * \param year
    * \param days the julian day in the year
    * \param seconds seconds in the day
    */
    void julian(const int year, const int days,const int seconds=0, const int zone=0);

    /** \brief 对当前时间进行小时的加减操作
    *
    * \param hours 小时
    */
    void addHours(const long hours);

    /** \brief 以天为单位增加/减少时间
    *
    * \param days 小时
    */
    void addDays(const long days);
    /** \brief 以分钟为单位增加/减少时间
    *
    * \param minutes 分钟数
    */
    void addMinutes(const long minutes);
    /** \brief 以分钟为单位增加/减少时间
    *
    * \param seconds 秒数
    */
    void addSeconds(const long seconds);

    /** \brief 返回当前时间的秒数值，取值范围在0～59
    */
    int getSecond() const
    {
        return m_second;
    }
    /** \brief return the seconds in a day */
    int getSeconds() const;

    /** \brief 返回当前时间的分钟值，取值范围在0～59
    */
    int getMinute() const
    {
        return m_minute;
    }

    /** \brief 返回当前时间的小时值，取值范围在0～23
    */
    int getHour() const
    {
        return m_hour;
    }

    /** \brief 返回当前时间在月份中的天数值，取值范围在1～31
    */
    int getDay() const
    {
        return m_day;
    }

    /** \brief 返回当前时间的月份，取值范围在1～12
    */
    int getMonth() const
    {
        return m_month;
    }

    /** \brief 返回当前时间的年
    */
    int getYear() const
    {
        return m_year;
    }
    /** \brief return the julian day: the days in the year, 1..366
    */
    int getDays() const;

    /** \brief 是否闰年
    *
    * \return 返回当前年是否为闰年
    */
    bool is_leap_year() const;
    /** \brief 是否闰年
    *
    * \param year 年份
    * \return 返回指定的年份是否为闰年
    */
    bool is_leap_year(const int year) const;

    /** 判断是否合法，不对年份进行判断 */
    bool is_valid() const;

    ///赋值操作
    Time& operator =(const Time&);
    ///当前时间加上一个TimeSpan，返回一个新的Time
    Time operator +(const TimeSpan& timeSpan) const;
    ///当前时间减去一个TimeSpan，返回一个新的Time
    Time operator -(const TimeSpan& timeSpan) const;
    ///前面的一个Time减去后面的一个Time，返回一个TimeSpan
    TimeSpan operator -(const Time& time) const;
    ///当前Time加上一个TimeSpan
    Time& operator +=(const TimeSpan&);
    /** \brief 当前Time减去一个TimeSpan
    *
    * \param span TimeSpan
    * \return 返回Time
    */
    Time& operator -=(const TimeSpan& span);
    ///判断两个时间是否相等，相等返回true
    bool operator ==(const Time& time) const;
    ///判断两个时间是否不等，不等返回true
    bool operator !=(const Time& time) const;
    ///判断前一个时间是否大于等于后一个时间，是返回true
    bool operator >=(const Time& time) const;
    ///判断前一个时间是否大于后一个时间，是返回true
    bool operator >(const Time& time) const;
    ///判断前一个时间是否小于等于后一个时间，是返回true
    bool operator <=(const Time& time) const;
    ///判断前一个时间是否小于后一个时间，是返回true
    bool operator <(const Time& time) const;
    ///输出操作
    friend std::ostream& operator<<(std::ostream& os, const Time& time);
};

///时间段，代表了一个时间跨度
class TimeSpan
{
private:
    long m_day;
    int m_hour;
    int m_minute;
    int m_second;
public:
    ///构造一个时间间隔为零秒的TimeSpan对象
    TimeSpan();
    ///用一个已有的TimeSpan对象构造一个新的TimeSpan对象
    TimeSpan(const TimeSpan& timeSpan);
    /** \brief 用time_t类型，即秒数初始化，构造一个TimeSpan对象
    * 只取年内的天数和时间值
    *
    */
    TimeSpan(time_t time,bool greenwich=true);

    /** \brief 通过指定具体的日数、小时数、分钟数、秒数构造时间段
    *
    * \param day 天数
    * \param hour 小时
    * \param minute 分钟
    * \param second 秒钟
    */
    TimeSpan(const long day, const int hour=0, const int minute=0, const int second=0);
    ///返回天数
    long getDay() const;
    ///返回当前天中的小时数，范围-23～23
    int getHour() const;
    ///返回当前小时中的分数，范围为-59～59
    int getMinute() const;
    ///返回当前分钟的秒数，范围为-59~59
    int getSecond() const;

    ///返回TimeSpan中完整的小时数目
    long getHours() const;
    ///返回TimeSpan中完整分钟的数目
    long getMinutes() const;
    ///返回TimeSpan中完整秒的数目
    long getSeconds() const;

    ///赋值操作
    TimeSpan& operator =(const TimeSpan& timeSpanSrc);
    ///两个TimeSpan相加
    TimeSpan operator +(const TimeSpan& timeSpan) const;
    ///两个TimeSpan相减
    TimeSpan operator -(const TimeSpan& timeSpan) const;
    ///前面的一个TimeSpan加上后面的一个TimeSpan
    TimeSpan& operator +=(const TimeSpan& timeSpan);
    ///前面的一个TimeSpan减去后面的一个TimeSpan
    TimeSpan& operator -=(const TimeSpan& timeSpan);
    ///判断两个TimeSpan是否相等
    bool operator ==(const TimeSpan& timeSpan) const;
    ///判断两个TimeSpan是否不等
    bool operator !=(const TimeSpan& timeSpan) const;
    ///判断前面一个TimeSpan是否大于等于后面的TimeSpan
    bool operator >=(const TimeSpan& timeSpan) const;
    ///判断前面一个TimeSpan是否小于等于后面的TimeSpan
    bool operator <=(const TimeSpan& timeSpan) const;
    ///判断前面的一个TimeSpan是否大于后面的一个TimeSpan
    bool operator >(const TimeSpan& timeSpan) const;
    ///判断前面的一个TimeSpan是否小于后面的一个TimeSpan
    bool operator <(const TimeSpan& timeSpan) const;
    ///输出操作
    friend std::ostream& operator<<(std::ostream& os,TimeSpan& timeSpan);
};

}

#endif

