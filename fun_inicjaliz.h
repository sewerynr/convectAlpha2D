#ifndef FUN_INICJALIZ_H
#define FUN_INICJALIZ_H

#include "fvCFD.H"
#include "fvOptions.H"

scalar funInitT1( double Psi, const dimensionedScalar& epsH, double par, const double x,const double y, const double z );

scalar funInitT2( double Psi, const dimensionedScalar& epsH, double par, const double x,const double y, const double z );

scalar funInitT3( double Psi, const dimensionedScalar& epsH, double par, const double x,const double y, const double z );\

scalar funInitT4( double Psi, const dimensionedScalar& epsH, double par, const double x,const double y, const double z );

scalar funInitPsi1( double x, double y, double z );

scalar funInitPsi1_06( double x, double y, double z );

scalar funInitPsi2( double x, double y, double z );

vector funInitV1( double x, double y, double z );

vector funInitV2( double x, double y, double z );

vector funInitV3( double x, double y, double z );

#endif // FUN_INICJALIZ_H
