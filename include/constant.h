#ifndef __LDAS_CONSTANT_H
#define __LDAS_CONSTANT_H

#ifndef M_E
#define M_E        2.71828182845904523536028747135      /* e */
#endif

#ifndef PI
#define PI (std::atan(1.0)*4)
#endif

//#ifndef LIGHT_SPEED
//#define LIGHT_SPEED	2.9979e8
//#endif

namespace ldas
{
const double EULER_NUMBER=2.71828182845904523536028747135;      /* e */
const double PI_NUMBER=3.1415926535897932384626433832795;
const double LIGHT_SPEED=2.9979e8;
const int MAXBUFFER=255;
const int NULL_VALUE=-9999;

const int __unused_int_one = 1;
#define is_bigendian() ( (*(char*)&__unused_int_one) == 0 )
#define SWAP_FLOAT(b)   { unsigned &a = *(unsigned *)(&b); a = (a>>24)|((a>>8)&0xFF00)|((a<<8)&0xFF0000)|(a<<24); }
#define SWAP_16(a) a = ((unsigned short)(a)>>8)|((unsigned short)(a)<<8);
#define SWAP_32(a) swap_32((unsigned long&)a);
inline void swap_32(unsigned long &a)
{
    a = ((unsigned)(a)>>24)|(((unsigned)(a)>>8)&0xFF00)|(((unsigned)(a)<<8)&0xFF0000)|((unsigned)(a)<<24);
}
enum DATA_TYPE {MODIS_LAI,AMSRE,SSMI};

/** 1: Fresnel reflection coeff. approximated by R(incident_angle)
 2: Fresnel reflection coeff. approximated by R(specular_angle)
 3: Fresnel reflection coeff. approximated by R(Transition)                                   c
 4: calculate 1,2,3 simultaneously
*/
enum FRESNEL_TYPE {AIEM_INCIDENT_ANGLE=1,AIEM_SPECULAR_ANGLE,AIEM_TRANSITION,AIEM_SIMULTANEOUS};
/** select what type of surface correlation
 1  Gaussian correlation function
 2  exponential correlation function
 3  transformed exponential correlation
 4  power-law spectrum correlation function
 5  x-power correlated surface
 6  x-exponential  correlated surface
 7  exponential-like  correlated surface
 */
enum CORRELATION_TYPE {AIEM_GAUSSIAN=1,AIEM_EXPONENT,AIEM_TRANS_EXPONENT,AIEM_POWER_LAW,AIEM_XPOWER,AIEM_XEXPONENT,AIEM_EXPNENT_LIKE};
/** AIEM MODE*/
enum AIEM_MODE {AIEM_ACTIVE=0,AIEM_PASSIVE};
}
#endif // __LDAS_CONSTANT_H
