#include "fun_inicjaliz.h"

scalar funInitT1( double Psi, const dimensionedScalar& epsH, double par, const double x,const double y, const double z )
{
    scalar ret = 0.;
    ret = 0.5*(scalar(1.) + Foam::tanh( Psi/(2.*scalar(par)*epsH.value()) ));
    return ret;
}

scalar funInitT2( double Psi, const dimensionedScalar& epsH, double par, const double x,const double y, const double z )
{
    scalar ret = 0.;
    ret = Foam::exp(-400.*(x*x + y*y));
    return ret;
}

scalar funInitT3( double Psi, const dimensionedScalar& epsH, double par, const double x,const double y, const double z )
{
    scalar ret = 0.;
    scalar r = Foam::sqrt(x*x + (y-0.25)*(y-0.25));
    return Foam::exp(-100.* r*r );
}

/*
 * alpha(psi) circular interface
 */
scalar funInitT4( double Psi, const dimensionedScalar& epsH, double par, const double x,const double y, const double z )
{
    const static scalar R = 0.15;
    scalar eps = epsH.value();
    scalar ret = 0.;
    scalar r = Foam::sqrt(x*x + (y-0.25)*(y-0.25));
    return 1-0.5*(1.0 + Foam::tanh((r-R)/2/eps));
}

scalar funInitPsi1( double x, double y, double z )
{
    scalar ret = 0.;
    ret = x;
    return ret;
}

scalar funInitPsi2( double x, double y, double z )
{
    scalar ret = 0.;
    ret = Foam::sqrt(x*x + y*y);
    return ret;
}

vector funInitV1( double x, double y, double z )
{
    vector v(y*10,- x*10 , 0);
    return v;
}

vector funInitV2( double x, double y, double z )
{
    using Foam::constant::mathematical::pi;
    using Foam::sin;
    using Foam::cos;
    using Foam::pow;
    const static double v0 = 2; // powinno byc 1
    const static double L = 1;  // dl domeny
    const static double k = pi/L;

    x += 0.5;
    y += 0.5;

    return vector(-v0*sin(k*x)*sin(k*x)*sin(2*k*y), v0*sin(k*y)*sin(k*y)*sin(2*k*x), 0);
}
// odwrucony wektor predkosci
vector funInitV3( double x, double y, double z )
{
    using Foam::constant::mathematical::pi;
    using Foam::sin;
    using Foam::cos;
    using Foam::pow;
    const static double v0 = 2;
    const static double L = 1; // dl domeny
    const static double k = pi/L;

    x += 0.5;
    y += 0.5;

    return vector(v0*sin(k*x)*sin(k*x)*sin(2*k*y), -v0*sin(k*y)*sin(k*y)*sin(2*k*x), 0);
}
