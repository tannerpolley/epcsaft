#pragma once
// Verification-only binary128 Base; include before cppad.hpp. No Engine ABI change.
#include <quadmath.h>
#include <array>
#include <cstdint>
#include <cstring>
#include <cppad/base_require.hpp>
namespace CppAD {
// CppAD promises eight-byte allocation alignment. Store the binary128 object
// representation in eight-byte-aligned words and copy to an aligned arithmetic
// temporary; never place a raw sixteen-byte-aligned scalar in its arena.
class Quad {
    std::array<std::uint64_t,2> words{};
public:
    Quad(__float128 v=0){std::memcpy(words.data(),&v,sizeof v);}
    __float128 raw() const {__float128 v;std::memcpy(&v,words.data(),sizeof v);return v;}
    explicit operator double() const {return double(raw());}
    explicit operator long double() const {return (long double)raw();}
    explicit operator int() const {return int(raw());}
    Quad operator-() const {return -raw();}
    Quad operator+() const {return *this;}
#define EPCSAFT_QUAD_BINARY(op) friend Quad operator op(Quad a,Quad b){return a.raw() op b.raw();}
    EPCSAFT_QUAD_BINARY(+) EPCSAFT_QUAD_BINARY(-) EPCSAFT_QUAD_BINARY(*) EPCSAFT_QUAD_BINARY(/)
#undef EPCSAFT_QUAD_BINARY
#define EPCSAFT_QUAD_COMPARE(op) friend bool operator op(Quad a,Quad b){return a.raw() op b.raw();}
    EPCSAFT_QUAD_COMPARE(==) EPCSAFT_QUAD_COMPARE(!=) EPCSAFT_QUAD_COMPARE(<) EPCSAFT_QUAD_COMPARE(<=) EPCSAFT_QUAD_COMPARE(>) EPCSAFT_QUAD_COMPARE(>=)
#undef EPCSAFT_QUAD_COMPARE
    Quad& operator+=(Quad b){return *this=*this+b;} Quad& operator-=(Quad b){return *this=*this-b;}
    Quad& operator*=(Quad b){return *this=*this*b;} Quad& operator/=(Quad b){return *this=*this/b;}
};
inline std::ostream& operator<<(std::ostream& os,const Quad& v){char b[128];quadmath_snprintf(b,sizeof b,"%.36Qg",v.raw());return os<<b;}
inline Quad CondExpOp(CompareOp c,const Quad& l,const Quad& r,const Quad& y,const Quad& n){return CondExpTemplate(c,l,r,y,n);}
CPPAD_COND_EXP_REL(Quad)
CPPAD_AZMUL(Quad)
inline bool EqualOpSeq(const Quad& x,const Quad& y){return x==y;}
inline bool IdenticalCon(const Quad&){return true;}
inline bool IdenticalZero(const Quad& x){return x==0;}
inline bool IdenticalOne(const Quad& x){return x==1;}
inline bool IdenticalEqualCon(const Quad& x,const Quad& y){return x==y;}
inline int Integer(const Quad& x){return int(x);}
inline bool GreaterThanZero(const Quad& x){return x>0;}
inline bool GreaterThanOrZero(const Quad& x){return x>=0;}
inline bool LessThanZero(const Quad& x){return x<0;}
inline bool LessThanOrZero(const Quad& x){return x<=0;}
inline bool abs_geq(const Quad& x,const Quad& y){return fabsq(x.raw())>=fabsq(y.raw());}
inline Quad abs(const Quad& x){return fabsq(x.raw());}
inline Quad sign(const Quad& x){return x>0?1:x<0?-1:0;}
#define EPCSAFT_QUAD_MATH(name) inline Quad name(const Quad& x){return name##q(x.raw());}
EPCSAFT_QUAD_MATH(acos) EPCSAFT_QUAD_MATH(asin) EPCSAFT_QUAD_MATH(atan) EPCSAFT_QUAD_MATH(cos) EPCSAFT_QUAD_MATH(cosh)
EPCSAFT_QUAD_MATH(exp) EPCSAFT_QUAD_MATH(fabs) EPCSAFT_QUAD_MATH(log) EPCSAFT_QUAD_MATH(log10) EPCSAFT_QUAD_MATH(sin)
EPCSAFT_QUAD_MATH(sinh) EPCSAFT_QUAD_MATH(sqrt) EPCSAFT_QUAD_MATH(tan) EPCSAFT_QUAD_MATH(tanh) EPCSAFT_QUAD_MATH(asinh)
EPCSAFT_QUAD_MATH(acosh) EPCSAFT_QUAD_MATH(atanh) EPCSAFT_QUAD_MATH(erf) EPCSAFT_QUAD_MATH(erfc) EPCSAFT_QUAD_MATH(expm1) EPCSAFT_QUAD_MATH(log1p)
#undef EPCSAFT_QUAD_MATH
inline Quad pow(const Quad& x,const Quad& y){return powq(x.raw(),y.raw());}
template<> class numeric_limits<Quad>{public:
    static Quad min(){return ldexpq(1,-16382);} static Quad max(){return ldexpq(2-ldexpq(1,-112),16383);}
    static Quad epsilon(){return ldexpq(1,-112);} static Quad quiet_NaN(){return nanq("");} static Quad infinity(){return __float128(__builtin_inf());}
    static const int digits10=33,max_digits10=36;
};
template<> struct to_string_struct<Quad>{std::string operator()(const Quad& v){char b[128];quadmath_snprintf(b,sizeof b,"%.36Qg",v.raw());return b;}};
}
#include <cppad/cppad.hpp>
