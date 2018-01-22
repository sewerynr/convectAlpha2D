#ifndef FUNKCJE_H
#define FUNKCJE_H

#include "fvCFD.H"
#include "fvOptions.H"

void InitPsiXYZ(volScalarField& Psi, const fvMesh& mesh, scalar (*funIntP)(const double x,const double y, const double z) ) ;

void InitT(volScalarField& T, const volScalarField& Psi, const fvMesh& mesh, const dimensionedScalar epsH,
           scalar (*funInitT)(double psi, const dimensionedScalar& epsH, double par, const double x,const double y, const double z), double par
           );

void limitT(volScalarField& T);

void limitT1(const fvMesh& mesh, volScalarField& T);

surfaceScalarField createPhiFieldEx(const Time& runTime, const fvMesh& mesh, const volScalarField& T, const volScalarField& Psi, const volScalarField& PsiZero );

surfaceScalarField createPhiFieldExC(const volScalarField& C, const Time& runTime, const fvMesh& mesh, const volScalarField& T, const volScalarField& Psi, const volScalarField& PsiZero );

surfaceScalarField createPhiFieldImp(const volScalarField& C, const Time& runTime, const fvMesh& mesh, const volScalarField& T, const volScalarField& Psi, const volScalarField& PsiZero );

volScalarField createKField(string nazwa,const Time& runTime, const fvMesh& mesh);

void fileOpener(FILE** f, string outpath);     // wskaźnik wskazuje plik ale jest jak zmienna jezeli chce zmienic wartosc na ktory wskazuje to musze pobrac adres do wskaznika **

void PsiZero1(const volScalarField& T, volScalarField& PsiZero, const double& eps, const dimensionedScalar& epsH, double gamma);

void PsiZero2(const volScalarField& T, volScalarField& PsiZero, const double& eps, const dimensionedScalar& epsH, double gamma);

void InitUXYZ( volVectorField& U, const fvMesh& mesh, vector (*funInitV) ( const double x, const double y, const double z) );

#endif // FUNKCJE_H
