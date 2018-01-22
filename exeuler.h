#ifndef EXEULER_H
#define EXEULER_H

#include <stdlib.h>
#include "funkcje.h"

#include "fvCFD.H"
#include "fvOptions.H"

void exEuler(const volScalarField& C, Time& runTime, const fvMesh& mesh, dimensionedScalar dtau,const volScalarField& Psi, volScalarField& PsiZero,
           volScalarField& T, volScalarField& Told, const dimensionedScalar& epsH, const double& eps, const bool& limitFieldT, const int& ilekrcz, double gamma, const bool& mapFunLog);

#endif // EXEULER_H
