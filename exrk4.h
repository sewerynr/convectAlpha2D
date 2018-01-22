#ifndef EXRK4_H
#define EXRK4_H

#include <stdlib.h>
#include "funkcje.h"

#include "fvCFD.H"
#include "fvOptions.H"

void exRK3(const volScalarField& C, Time& runTime, const fvMesh& mesh, dimensionedScalar dtau,const volScalarField& Psi, volScalarField& PsiZero,
           volScalarField& T, volScalarField& Told, const dimensionedScalar& epsH, const double& eps, const bool& limitFieldT, const int& ilekrcz, double gamma, const bool& mapFunLog, int ilePkt);

#endif // EXRK4_H
