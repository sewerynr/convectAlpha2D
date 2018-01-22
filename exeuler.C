#include "exeuler.h"


void exEuler(const volScalarField& C, Time& runTime, const fvMesh& mesh, dimensionedScalar dtau,const volScalarField& Psi, volScalarField& PsiZero,
           volScalarField& T, volScalarField& Told, const dimensionedScalar& epsH, const double& eps, const bool& limitFieldT, const int& ilekrcz, double gamma, const bool& mapFunLog)
{
    PsiZero = epsH*Foam::log((T+eps)/(1-T+eps));
    surfaceScalarField phiR = createPhiFieldExC(C, runTime, mesh, T, Psi, PsiZero);

    Info<< "Explicit Euler !!! "  << endl;
    Info<< "Tsize =  "<< T.size()  << endl;
    scalar one(1);
    dimensionedScalar SMALL_NUMBER("small", dimless, SMALL);

    Told == T;     //        Told.deepCopy(T);   // to samo co: Told == T

    for(unsigned int i = 0; i <= ilekrcz; ++i)
    {
        T == Told + fvc::div(phiR)*dtau;
        scalar norm1 = 0.0;
        norm1 = Foam::sum(Foam::mag(T-Told)).value();
        norm1 = norm1 / T.size();
        Told == T;

        if (limitFieldT)
            limitT(T);

        PsiZero = epsH*Foam::log((T+eps)/(1-T+eps)) ;
        phiR = linearInterpolate( C*T*( scalar(1.) - T )*( mag(fvc::grad(PsiZero))- scalar(1.) ) * (fvc::grad(Psi) /( mag(fvc::grad(Psi)) + SMALL_NUMBER ))  )  & mesh.Sf();
        Info<< "Norma 1 = " << norm1 << endl;
        runTime.write(); // wszystko z AUTO_WRITE
        ++runTime;
    }

    Info<< "Norma 1 = " << phiR << endl;
}
