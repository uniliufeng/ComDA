#include "time.hpp"
#include <sstream>

namespace ldas
{
Time::Time():m_year(0),m_month(1),m_day(1),m_hour(0),m_minute(0),m_second(0)
{
}
Time::Time(int year, int month, int day, int hour, int min, int sec,int zone)
    :m_year(year),m_month(month),m_day(day),m_hour(hour),m_minute(min),m_second(sec)
{
    addHours(-zone);
}
void Time::julian(const int year, const int days,const int seconds, const int zone)
{
    m_year=year;
    m_month=1;
    m_day=0;
    m_hour=0;
    m_minute=0;
    m_second=0;
    addDays(days);
    addSeconds(seconds);
    addHours(-zone);
}
Time::Time(time_t time,bool greenwich)
{
    tm* ptm;
    if (greenwich)
        ptm=gmtime(&time);
    else
        ptm=localtime(&time);
    m_second=ptm->tm_sec;
    m_minute=ptm->tm_min;
    m_hour=ptm->tm_hour;
    m_day=ptm->tm_mday;
    m_month=ptm->tm_mon;
    m_year=ptm->tm_year;
}

Time::Time(const Time& src)
{
    m_year   = src.getYear();
    m_month  = src.getMonth();
    m_day    = src.getDay();
    m_hour = src.getHour();
    m_minute = src.getMinute();
    m_second = src.getSecond();
}

Time::Time(const std::string& src)
{
    int year,month,day,hour,minute,second,zoneh,zonem;
    year=hour=minute=second=zoneh=zonem=0;
    month=day=1;
    //ISO 8601
    size_t found,prior;
    std::string str;
    std::istringstream ist;
    found=src.find_first_of("-");
    prior=0;
    if (found==std::string::npos)
    {
        ist.str(src);
        ist>>year;
    }
    else
    {
        str=src.substr(prior,found);
        ist.str(str);
        ist>>year;
        prior=found;
        found=src.find_first_of("-",found+1);
        if (found==std::string::npos)
        {
            str=src.substr(prior+1,found);
            ist.clear();
            ist.str(str);
            ist>>month;
            std::cout<<month<<std::endl;
        }
        else
        {
            str=src.substr(prior+1,found-1);
            ist.clear();
            ist.str(str);
            ist>>month;
            prior=found;
            found=src.find_first_of(" T",found+1);
            if (found==std::string::npos)
            {
                str=src.substr(prior+1,found);
                ist.clear();
                ist.str(str);
                ist>>day;
            }
            else
            {
                str=src.substr(prior+1,found-1);
                ist.clear();
                ist.str(str);
                ist>>day;
                prior=found;
                found=src.find_first_of(":",found+1);
                if (found==std::string::npos)
                {
                    str=src.substr(prior+1,found);
                    ist.clear();
                    ist.str(str);
                    ist>>hour;
                }
                else
                {
                    str=src.substr(prior+1,found-1);
                    ist.clear();
                    ist.str(str);
                    ist>>hour;
                    prior=found;
                    found=src.find_first_of(":",found+1);
                    //deal with Complete date plus hours and minutes: YYYY-MM-DDThh:mmTZD (eg 1997-07-16T19:20+01:00)
                    if (src.find_first_of("+-",prior+1)<found)
                    {
                        str=src.substr(src.find_first_of("+-",prior+1));
                        ist.clear();
                        ist.str(str);
                        ist>>zoneh;
                        ist>>zonem;
                    }
                    if (found==std::string::npos)
                    {
                        str=src.substr(prior+1,found);
                        ist.clear();
                        ist.str(str);
                        ist>>minute;
                    }
                    else
                    {
                        str=src.substr(prior+1,found-1);
                        ist.clear();
                        ist.str(str);
                        ist>>minute;
                        prior=found;
                        found=src.find_first_of("+-Z",found+1);
                        if (found==std::string::npos)
                        {
                            str=src.substr(prior+1,found);
                            ist.clear();
                            ist.str(str);
                            ist>>second;
                        }
                        else if (src[found]=='Z')
                        {
                            str=src.substr(prior+1,found-1);
                            ist.clear();
                            ist.str(str);
                            ist>>second;
                        }
                        else
                        {
                            str=src.substr(prior+1,found-1);
                            ist.clear();
                            ist.str(str);
                            ist>>second;
                            prior=found;
                            found=src.find_first_of(":",found+1);
                            if (found==std::string::npos)
                            {
                                str=src.substr(prior,found);
                                ist.clear();
                                ist.str(str);
                                ist>>zoneh;
                            }
                            else
                            {
                                str=src.substr(prior,found-1);
                                ist.clear();
                                ist.str(str);
                                ist>>zoneh;
                                str=src.substr(found);
                                ist.clear();
                                ist.str(str);
                                ist>>zonem;
                            }
                        }
                    }
                }
            }
        }
    }

    m_year=year;
    m_month=month;
    m_day=day;
    m_hour=hour;
    m_minute=minute;
    m_second=second;
    addMinutes(-zonem);
    addHours(-zoneh);
}
int Time::getSeconds() const
{
    return m_second+m_minute*60+m_hour*3600;
}
int Time::getDays() const
{
    if (m_month==1) return m_day;
    int month_day[12]= {31, 28, 31, 30, 31, 30,
                        31, 31, 30, 31, 30, 31
                       };
    if (is_leap_year())
        month_day[1] = 29;
    int day=m_day;
    for(int i=0; i<m_month-1; i++)
    {
        day+=month_day[i];
    }
    return day;
}
bool Time::is_leap_year() const
{
    if(((m_year%4)==0 && (m_year%100)!=0)||
            (m_year%400==0))
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool Time::is_leap_year(const int year) const
{
    if(((year%4)==0 && (year%100)!=0)||
            (year%400==0))
    {
        return true;
    }
    else
    {
        return false;
    }
}
void Time::addSeconds(const long seconds)
{
    if (seconds==0) return;
    int second=(seconds+m_second)%60;
    int minute=(seconds+m_second)/60;
    if (second<0)
    {
        --minute;
        second+=60;
    }
    m_second=second;
    addMinutes(minute);
}
void Time::addMinutes(const long minutes)
{
    if (minutes==0) return;
    int minute=(minutes+m_minute)%60;
    int hour=(minutes+m_minute)/60;
    if (minute<0)
    {
        --hour;
        minute+=60;
    }
    m_minute=minute;
    addHours(hour);
}
void Time::addHours(const long hours)
{
    if (hours==0) return;
    //处理小时数
    int hour=(hours+m_hour)%24;
    int day=(hours+m_hour)/24;
    if (hour<0)
    {
        --day;
        hour+=24;
    }
    m_hour=hour;
    addDays(day);
}
void Time::addDays(const long days)
{
    if (days==0) return;
    int month_day[12]= {31, 28, 31, 30, 31, 30,
                        31, 31, 30, 31, 30, 31
                       };
    //从年开始计算的天数
    int month_days[12];
    if(is_leap_year())
    {
        month_day[1] = 29;
    }
    month_days[0]=31;
    for(int i=0; i<11; i++)
    {
        month_days[i+1]=month_days[i]+month_day[i+1];
    }
    //是否跳出当前年范围
    int ydays;
    if (m_month==1)
        ydays=m_day+days;
    else
        ydays=month_days[m_month-2]+m_day+days;
    //是否大于四年？
    int f_years=ydays/(365*3+366);
    ydays=ydays%(365*3+366);
    m_year+=f_years*4;
    //保证在负4年内
    if (ydays<0)
    {
        m_year-=4;
        ydays+=365*3+366;
    }
    int year_days=365;
    if (is_leap_year()) year_days++;
    if (ydays>year_days)
    {
    	for(int i=1;i<4;i++)
		{
			if(is_leap_year(m_year+i))
				year_days+=366;
			else
				year_days+=365;
			if (ydays<=year_days)
			{
				m_year+=i;
				ydays-=year_days-365;
				if (is_leap_year()) ydays++;
				break;
			}
		}
    }
	if(is_leap_year())
        month_day[1] = 29;
    else
		month_day[1] = 28;
    for(int i=0; i<11; i++)
    {
        month_days[i+1]=month_days[i]+month_day[i+1];
    }
    if (ydays>0)
    {
        if (ydays<=31)
        {
            m_month=1;
            m_day=ydays;
        }
        else
        {
            for(int i=0; i<11; i++)
            {
                if (ydays>month_days[i] && ydays<=month_days[i+1])
                {
                    m_month=i+2;
                    m_day=ydays-month_days[i];
                }
            }
        }
    }
    else //ydays=0
    {
        --m_year;
        m_day=31;
        m_month=12;
    }
}

bool Time::is_valid() const
{
    if (m_month<1 || m_month>12)
        return false;
    else if (m_day<1 || m_day>31)
        return false;
    else if (m_hour<0 || m_hour>23)
        return false;
    else if (m_minute<0 || m_minute>59)
        return false;
    else if (m_second<0 || m_second>59)
        return false;
    else
    {
        int month_day[12]= {31, 28, 31, 30, 31, 30,
                            31, 31, 30, 31, 30, 31
                           };
        if (is_leap_year())
            month_day[1]=29;
        if (m_day>month_day[m_month-1])
            return false;
        else
            return true;
    }
}

std::ostream& operator<<(std::ostream& os, const Time& time)
{
    os<<time.getYear()<<"-";
    os.width(2);
    os.fill('0');
    os<<time.getMonth()<<"-";
    os.width(2);
    os.fill('0');
    os<<time.getDay()<<" ";
    os.width(2);
    os.fill('0');
    os<<time.getHour()<<":";
    os.width(2);
    os.fill('0');
    os<<time.getMinute()<<":";
    os.width(2);
    os.fill('0');
    os<<time.getSecond();
    return os;
}

Time& Time::operator =(const Time& src)
{
    m_year   = src.getYear();
    m_month  = src.getMonth();
    m_day    = src.getDay();
    m_hour = src.getHour();
    m_minute = src.getMinute();
    m_second = src.getSecond();
    return *this;
}

Time& Time::operator +=(const TimeSpan& span)
{
    addSeconds(span.getSeconds());
    return *this;
}

Time& Time::operator -=(const TimeSpan& span)
{
    addSeconds(-span.getSeconds());
    return *this;
}

Time Time::operator -(const TimeSpan& timeSpan) const
{
    Time temp(*this);
    return temp-=timeSpan;
}

Time Time::operator +(const TimeSpan& timeSpan) const
{
    Time temp(*this);
    return temp+=timeSpan;
}

TimeSpan Time::operator -(const Time& time) const
{
    int year=m_year-time.getYear();
    int day=year*365+year/4;
    int min_year=std::min(m_year,time.getYear());
    if (year!=0 && (year%4>0))
    {
        //年份中存在闰年，需要进行特别处理
        for(int i=0; i<(year%4); i++)
        {
            if (is_leap_year(min_year+i))
            {
                day++;
                break;
            }
        }
    }
    day+=m_day-time.getDay();
    //处理月份中的天数
    int month_day[12]= {31, 28, 31, 30, 31, 30,
                        31, 31, 30, 31, 30, 31
                       };
    if (is_leap_year())
        month_day[1] = 29;
    if (m_month>1)
        for(int i=1; i<=m_month-1; i++)
            day+=month_day[i-1];
    if (!is_leap_year(time.getYear()))
        month_day[1]=28;
    if (time.getMonth()>1)
        for(int i=1; i<time.getMonth()-1; i++)
            day-=month_day[i-1];

    int second=m_second-time.getSecond();
    int minute=m_minute;
    if (second<0)
    {
        minute--;
        second+=60;
    }
    minute-=time.getMinute();
    int hour=m_hour;
    if (minute<0)
    {
        hour--;
        minute+=60;
    }
    hour-=time.getHour();
    if (hour<0)
    {
        day--;
        hour+=24;
    }
    return TimeSpan(day,hour,minute,second);
}

bool Time::operator ==(const Time& time) const
{
    if (m_year==time.getYear() && m_month==time.getMonth() && m_day==time.getDay() && m_hour==time.getHour() && m_minute==time.getMinute() && m_second==time.getSecond())
        return true;
    return false;
}

bool Time::operator !=(const Time& time) const
{
    if (m_year==time.getYear() && m_month==time.getMonth() && m_day==time.getDay() && m_hour==time.getHour() && m_minute==time.getMinute() && m_second==time.getSecond())
        return false;
    return true;
}

bool Time::operator >(const Time& time) const
{
    if (m_year>time.getYear())
        return true;
    else if (m_year<time.getYear())
        return false;
    else if (m_month>time.getMonth())
        return true;
    else if (m_month<time.getMonth())
        return false;
    else if (m_day>time.getDay())
        return true;
    else if (m_day<time.getDay())
        return false;
    else if (m_hour>time.getHour())
        return true;
    else if (m_hour<time.getHour())
        return false;
    else if (m_minute>time.getMinute())
        return true;
    else if (m_minute<time.getMinute())
        return false;
    else if (m_second>time.getSecond())
        return true;
    else
        return false;
}

bool Time::operator >=(const Time& time) const
{
    if (m_year>time.getYear())
        return true;
    else if (m_year<time.getYear())
        return false;
    else if (m_month>time.getMonth())
        return true;
    else if (m_month<time.getMonth())
        return false;
    else if (m_day>time.getDay())
        return true;
    else if (m_day<time.getDay())
        return false;
    else if (m_hour>time.getHour())
        return true;
    else if (m_hour<time.getHour())
        return false;
    else if (m_minute>time.getMinute())
        return true;
    else if (m_minute<time.getMinute())
        return false;
    else if (m_second>time.getSecond())
        return true;
    else if (m_second<time.getSecond())
        return false;
    else
        return true;
}

bool Time::operator <(const Time& time) const
{
    return !(*this>=time);
}

bool Time::operator <=(const Time& time) const
{
    return !(*this>time);
}

// implementation of class TimeSpan
TimeSpan::TimeSpan():m_day(0),m_hour(0),m_minute(0),m_second(0)
{
}

TimeSpan::TimeSpan(const TimeSpan& span)
{
    m_day=span.getDay();
    m_hour=span.getHour();
    m_minute=span.getMinute();
    m_second=span.getSecond();
}

TimeSpan::TimeSpan(time_t time,bool greenwich)
{
    tm* ptm;
    if (greenwich)
        ptm=gmtime(&time);
    else
        ptm=localtime(&time);
    m_second=ptm->tm_sec;
    m_minute=ptm->tm_min;
    m_hour=ptm->tm_hour;
    m_day=ptm->tm_yday;
}

TimeSpan::TimeSpan(const long day,const int hour,const int minute,const int second):m_day(day),m_hour(hour),m_minute(minute),m_second(second)
{
}

long TimeSpan::getDay() const
{
    return m_day;
}
int TimeSpan::getHour() const
{
    return m_hour;
}
int TimeSpan::getMinute() const
{
    return m_minute;
}
int TimeSpan::getSecond() const
{
    return m_second;
}


long TimeSpan::getHours() const
{
    return (m_day*24+m_hour);
}

long TimeSpan::getMinutes() const
{
    return (m_day*24+m_hour)*60+m_minute;
}

long TimeSpan::getSeconds() const
{
    return ((m_day*24+m_hour)*60+m_minute)*60+m_second;
}

TimeSpan& TimeSpan::operator =(const TimeSpan& src)
{
    m_day=src.getDay();
    m_minute=src.getMinute();
    m_second=src.getSecond();
    m_hour=src.getHour();
    return *this;
}

TimeSpan& TimeSpan::operator +=(const TimeSpan &src)
{
    m_day+=src.getDay();
    m_minute+=src.getMinute();
    m_second+=src.getSecond();
    m_hour+=src.getHour();
    m_minute+=(m_second/60);
    m_second=(m_second%60);
    m_hour+=(m_minute/60);
    m_minute=(m_minute%60);
    m_day+=(m_hour/24);
    m_hour=(m_hour%24);
    return *this;
}

TimeSpan TimeSpan::operator +(const TimeSpan& timeSpan) const
{
    TimeSpan temp(*this);
    return temp+=timeSpan;
}

TimeSpan& TimeSpan::operator -=(const TimeSpan &src)
{
    m_day-=src.getDay();
    m_minute-=src.getMinute();
    m_second-=src.getSecond();
    m_hour-=src.getHour();
    m_minute+=(m_second/60);
    m_second=(m_second%60);
    m_hour+=(m_minute/60);
    m_minute=(m_minute%60);
    m_day+=(m_hour/24);
    m_hour=(m_hour%24);
    return *this;
}

TimeSpan TimeSpan::operator -(const TimeSpan& timeSpan) const
{
    TimeSpan temp(*this);
    return temp-=timeSpan;
}

bool TimeSpan::operator ==(const TimeSpan& src) const
{
    return (getSeconds()==src.getSeconds());
}

bool TimeSpan::operator !=(const TimeSpan& src) const
{
    return (getSeconds()!=src.getSeconds());
}

bool TimeSpan::operator <(const TimeSpan& src) const
{
    return (getSeconds()<src.getSeconds());
}

bool TimeSpan::operator <=(const TimeSpan& src) const
{
    return (getSeconds()<=src.getSeconds());
}

bool TimeSpan::operator >(const TimeSpan& src) const
{
    return (getSeconds()>src.getSeconds());
}

bool TimeSpan::operator >=(const TimeSpan& src) const
{
    return (getSeconds()>=src.getSeconds());
}

std::ostream& operator<<(std::ostream& os, TimeSpan& src)
{
    os<<src.getDay()<<" "<<src.getHour()<<":"<<src.getMinute()<<":"<<src.getSecond();
    return os;
}

}




