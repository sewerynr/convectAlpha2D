
#include <stdlib.h>
#include "funkcje.h"

void InitPsiXYZ(volScalarField& Psi, const fvMesh& mesh, scalar (*funIntP)(const double x,const double y, const double z) )
{
    forAll(mesh.cellCentres(),cellI)
    {
        Psi[cellI] = funIntP( mesh.cellCentres()[cellI].x(),mesh.cellCentres()[cellI].y(),mesh.cellCentres()[cellI].z() );
    }

    forAll(mesh.boundary(), patchi)  // mesh.boundary() daje liste adresow do war. brzeg.
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        fvPatchScalarField& psiPatch = Psi.boundaryFieldRef()[patchi];

        forAll(patch, faceId) // petla po centrach objetosci danego patcha
        {
            psiPatch[faceId] = funIntP( patch.Cf()[faceId].x(),  patch.Cf()[faceId].y(), patch.Cf()[faceId].z());
        }
    }
}

void InitT(volScalarField& T, const volScalarField& Psi, const fvMesh& mesh, const dimensionedScalar epsH,
           scalar (*funInitT)(double psi, const dimensionedScalar& epsH, double par,  double x, double y,  double z) , double par
           )
{
    forAll(mesh.cellCentres(),cellI)
    {
        double p = Psi[cellI];
        T[cellI] = funInitT(p, epsH, par, mesh.cellCentres()[cellI].x(),mesh.cellCentres()[cellI].y(),mesh.cellCentres()[cellI].z() );
    }

    forAll(mesh.boundary(), patchi)  // mesh.boundary() daje liste adresow do war. brzeg.
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        fvPatchScalarField& TPatch = T.boundaryFieldRef()[patchi];
        const fvPatchScalarField& psiPatch = Psi.boundaryField()[patchi];

        forAll(patch, faceId) // petla po centrach objetosci danego patcha
        {
            double pBon = psiPatch[faceId];
            TPatch[faceId] = funInitT( pBon, epsH, par, mesh.cellCentres()[faceId].x(),mesh.cellCentres()[faceId].y(),mesh.cellCentres()[faceId].z() );
        }
    }
}

void InitUXYZ( volVectorField& U, const fvMesh& mesh, vector (*funInitV) (const double x, const double y, const double z) )
{
    forAll (mesh.cellCentres(), I)
    {
        U[I] = funInitV( mesh.cellCentres()[I].x(), mesh.cellCentres()[I].y(), mesh.cellCentres()[I].z());
    }
    forAll (mesh.boundary(), patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];
        fvPatchVectorField& UPach = U.boundaryFieldRef()[patchI];

        forAll(patch, faceId)
        {
            UPach[faceId] = funInitV( mesh.cellCentres()[faceId].x(), mesh.cellCentres()[faceId].y(), mesh.cellCentres()[faceId].z());
        }
    }
    U.correctBoundaryConditions();
}

void limitT1(const fvMesh& mesh, volScalarField& T)
{
    forAll(mesh.cellCentres(),cellI)
    {
        if (T[cellI] > 1.)
        {
            T[cellI] = 1.;
        }
        else if (T[cellI] < 0.)
        {
            T[cellI] = 0.;
        }
    }

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        fvPatchScalarField& TPatch = T.boundaryFieldRef()[patchi];

        forAll(patch, faceId)
        {
            if (TPatch[faceId] > 1.)
            {
                TPatch[faceId] = 1.;
            }
            else if (TPatch[faceId] < 0.)
            {
                TPatch[faceId] = 0.;
            }
        }
    }
}

void limitT(volScalarField& T)
{
    T = Foam::min(Foam::max(T, scalar(0) ), scalar(1) );
}

surfaceScalarField createPhiFieldEx(const Time& runTime, const fvMesh& mesh, const volScalarField& T, const volScalarField& Psi, const volScalarField& PsiZero )
{
    dimensionedScalar SMALL_NUMBER("small", dimless, SMALL);
    surfaceScalarField phiR  // cv surface!
        (
          IOobject
          (
              "phiR",
              runTime.timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE
          ),
          linearInterpolate( T*(scalar(1.) - T)*(mag(fvc::grad(PsiZero))- scalar(1.) ) /( mag(fvc::grad(PsiZero)) + SMALL_NUMBER ) ) *  fvc::snGrad(PsiZero) * mesh.magSf()
        );
    return phiR;
}

surfaceScalarField createPhiFieldExC(const volScalarField& C, const Time& runTime, const fvMesh& mesh, const volScalarField& T, const volScalarField& Psi, const volScalarField& PsiZero )
{
    dimensionedScalar SMALL_NUMBER("small", dimless, SMALL);
    surfaceScalarField phiR  // cv surface!
        (
          IOobject
          (
              "phiR",
              runTime.timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE
          ), linearInterpolate( C*T*( scalar(1.) - T )*( mag(fvc::grad(PsiZero))- scalar(1.) ) * (fvc::grad(PsiZero) /( mag(fvc::grad(PsiZero)) + SMALL_NUMBER ))  ) & mesh.Sf()
        );
    return phiR;
}

surfaceScalarField createPhiFieldImp(const volScalarField& C, const Time& runTime, const fvMesh& mesh, const volScalarField& T, const volScalarField& Psi, const volScalarField& PsiZero )
{
    dimensionedScalar SMALL_NUMBER("small", dimless, SMALL);
    surfaceScalarField phiR  // cv surface!
        (
          IOobject
          (
              "phiR",
              runTime.timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE
          ),
          linearInterpolate( C*( scalar(1.) - T )*( mag(fvc::grad(PsiZero))- scalar(1.) ) * (fvc::grad(PsiZero) /( mag(fvc::grad(PsiZero)) + SMALL_NUMBER ))  )  & mesh.Sf()
//          linearInterpolate( (scalar(1.)-T)*(mag(fvc::grad(PsiZero))-scalar(1.)) /( mag(fvc::grad(Psi)) + SMALL_NUMBER) ) *  fvc::snGrad(Psi) * mesh.magSf()
        );
    return phiR;
}

volScalarField createKField(string nazwa, const Time& runTime, const fvMesh& mesh)
{
    volScalarField k
    (
        IOobject
        (
            nazwa,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(nazwa, dimless, scalar(0))
    );
    return k;
}


void fileOpener(FILE** f, string outpath)
{
    *f = fopen(outpath.c_str(), "r");
    if(NULL == f)
       {
         fprintf(stderr, "!!! fopen() _nr_1 failed.\n");
       }
}

void PsiZero1(const volScalarField& T, volScalarField& PsiZero, const double& eps, const dimensionedScalar& epsH, double gamma)
{
    PsiZero = epsH*Foam::log((T+eps)/(1-T+eps));
}

void PsiZero2(const volScalarField& T, volScalarField& PsiZero, const double& eps, const dimensionedScalar& epsH, double gamma)
{
    dimensionedScalar mn( "mn", dimLength , 1. ); // ?????????????????
    PsiZero = mn * Foam::pow((T + eps),gamma)/( Foam::pow((T + eps), gamma) + Foam::pow((1 - T + eps), gamma) );
}
