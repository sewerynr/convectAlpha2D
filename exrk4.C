#include "exrk4.h"
#include <iostream>
#include <fstream>
#include "basicexceptions.h"

dimensionedScalar SMALL_NUMBER("small", dimless, SMALL);

void updateGradPsi(const fvMesh& mesh, const volScalarField & Psi, volScalarField & gradPsi)
{
    gradPsi == mag(fvc::grad(Psi));
    forAll(gradPsi, cellId)
    {
        if(gradPsi[cellId]< SMALL_NUMBER.value())
        {
            gradPsi[cellId] = 1.;
        }
    }
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        fvPatchScalarField& gradPsiPatch = gradPsi.boundaryFieldRef()[patchi];
        forAll(patch, faceId) // petla po centrach objetosci danego patcha
        {
            if(gradPsiPatch[faceId]< SMALL_NUMBER.value())
            {
                gradPsiPatch[faceId] = 1.;
            }
        }
    }
}

void LimitGradPsi(const fvMesh& mesh, const volScalarField & Psi, volScalarField & gradPsi, double dx)
{
    forAll( mesh.cellCentres(), cellId)
    {
        if( (Psi[cellId] > 14*dx) || (Psi[cellId] < -14*dx) )
        {
            gradPsi[cellId] = 1.;
        }
    }

    forAll(mesh.boundary(), patchi)  // mesh.boundary() daje liste adresow do war. brzeg.
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        fvPatchScalarField& gradPsiPatch = gradPsi.boundaryFieldRef()[patchi];
        const fvPatchScalarField& PsiPatch = Psi.boundaryField()[patchi];

        forAll(patch, faceId) // petla po centrach objetosci danego patcha
        {
            if ( (PsiPatch[faceId]  > 14*dx) || (PsiPatch[faceId] < -14*dx) )
            {
                gradPsiPatch[faceId] = 1.;
//                Info << PsiPatch[faceId] << endl;
            }
        }
    }
}

// przekazuje PsiZero (rzeczywiste psi) jako Psi wiec wszystok licze z aktualnego psi nie z psi (rzeczywistego PsiZero)
void exRK3(const volScalarField& C, Time& runTime, const fvMesh& mesh, dimensionedScalar dtau,const volScalarField& Psi, volScalarField& PsiZero,
           volScalarField& T, volScalarField& Told, const dimensionedScalar& epsH, const double& eps, const bool& limitFieldT, const int& ilekrcz, double gamma, const bool& mapFunLog, int ilePkt)
    {
    std::fstream f;
    string nazwa, sciezka, metcalk;
    metcalk = "exRK3";
    std::ostringstream strs;
    strs << eps;
    std::string wielkosceps = strs.str();

    std::ostringstream strs2;
    strs2 << ilePkt;
    std::string ilePktStr = strs2.str();

    sciezka = "/home/sr/OpenFOAM/sr-4.1/run/convectAlphaSTest2D/";
    try
    {
    f.open(sciezka + "norma_" + metcalk + "_eps_" + wielkosceps + "_" + ilePktStr + ".vtk" , std::ios::out);
    if( f.good())
        Info << "plik otwarty" << endl;
    else
        throw FileException();
    }
    catch(BasicException& ex)
    { Info << ex.wyjatek << endl; }

    Info<< "Explicit RK3 !!! "  << endl;
    scalar one(1);

//    if (limitFieldT)
//        limitT(T);

    if (mapFunLog)
        PsiZero1(T, PsiZero, eps, epsH, gamma);
    else
        PsiZero2(T, PsiZero, eps, epsH, gamma);

    surfaceScalarField phiR = createPhiFieldExC(C, runTime, mesh, T, Psi, PsiZero);

    volScalarField gradPsi
    (
        IOobject
        (
            "gradPsi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mag(fvc::grad(PsiZero))
    );

    volScalarField gradTanalit
    (
        IOobject
        (
            "gradTanalit",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mag( ( one - Foam::tanh(Psi/(2*epsH))*Foam::tanh(Psi/(2*epsH)) )*(1/(4*epsH)) )
    );

    volScalarField magGradT
    (
        IOobject
        (
            "magGradT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mag(fvc::grad(T))
    );

    volScalarField GradT
    (
        IOobject
        (
            "GradT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        T*(1-T)*mag(fvc::grad(PsiZero))/epsH
    );

    volScalarField LaplaceT
    (
        IOobject
        (
            "LaplaceT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        T*(1-T)*(fvc::laplacian(PsiZero))/epsH
    );

    volScalarField LaplaceTStraignt
    (
        IOobject
        (
            "LaplaceTStraignt",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::laplacian(T)
    );

    volScalarField TAnalit
    (
        IOobject
        (
            "TAnalit",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        T
    );
    volScalarField TSt
    (
        IOobject
        (
            "TSt",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        T
    );

    volScalarField LaplaceTW1D
        (
            IOobject
            (
                "LaplaceTW1D",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            T*(1.-T)*( fvc::laplacian(PsiZero) + (1./epsH)*( fvc::grad(PsiZero)&fvc::grad(PsiZero) )*(1.-2.*T) )/epsH
         );

    volScalarField LaplaceAnalit
        (
            IOobject
            (
                "LaplaceAnalit",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
                - ( Foam::tanh(Psi/(2.0*epsH)) - Foam::tanh(Psi/(2.0*epsH))*Foam::tanh(Psi/(2.0*epsH))*Foam::tanh(Psi/(2.0*epsH)) )/(4*epsH*epsH)
         );

    surfaceScalarField surfgradT
    (
        IOobject
        (
            "surfgradT",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
       linearInterpolate( T*(1-T)*(fvc::grad(PsiZero))/epsH ) & mesh.Sf()
    );

    volScalarField GradTW
    (
        IOobject
        (
            "GradTW",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        T*(1.-T)*mag(fvc::grad(PsiZero))/epsH
    );
    TAnalit = 0.5*(1.+Foam::tanh(Psi/(2.*epsH)));



    Told == T;

    updateGradPsi(mesh, PsiZero, gradPsi);
    LimitGradPsi(mesh, PsiZero, gradPsi, 1./ilePkt);

    //updateGradPsi(mesh, PsiZero, gradPsi);
    volScalarField k1 = createKField("k1", runTime, mesh);

    k1 =  T + fvc::div(phiR)*dtau;

    limitT(k1);
    PsiZero1(k1, PsiZero, eps, epsH, gamma);
    gradPsi = mag(fvc::grad(PsiZero));

    updateGradPsi(mesh, PsiZero, gradPsi);
    LimitGradPsi(mesh, PsiZero, gradPsi, 1./ilePkt);

    phiR = linearInterpolate( C*k1*( scalar(1.) - k1 )*( gradPsi- scalar(1.) ) * (fvc::grad(PsiZero) /( gradPsi )) ) & mesh.Sf();

    volScalarField k2  = createKField("k2", runTime, mesh);

    k2 = 3./4* T + 1./4* k1 + 1./4* fvc::div(phiR)*dtau;

    limitT(k2);
    PsiZero1(k2, PsiZero, eps, epsH, gamma);
    gradPsi = mag(fvc::grad(PsiZero));

    updateGradPsi(mesh, PsiZero, gradPsi);
    LimitGradPsi(mesh, PsiZero, gradPsi, 1./ilePkt);

    phiR = linearInterpolate( C*k2*( scalar(1.) - k2 )*( gradPsi- scalar(1.) ) * (fvc::grad(PsiZero) /( gradPsi )) ) & mesh.Sf();
    T =  1./3* T + 2./3* k2 + 2./3*fvc::div(phiR)*dtau;

    limitT(T);
    PsiZero1(T, PsiZero, eps, epsH, gamma);
    gradPsi = mag(fvc::grad(PsiZero));
    updateGradPsi(mesh, PsiZero, gradPsi);
    LimitGradPsi(mesh, PsiZero, gradPsi, 1./ilePkt);

//    if (limitFieldT)
//        limitT(T);

    double norm1c = Foam::sum(Foam::mag(T-Told)).value() / T.size();
    Info << "Norma 1 = " << norm1c << endl;

//    if (mapFunLog)
//        PsiZero1(T, PsiZero, eps, epsH, gamma);
//    else
//        PsiZero2(T, PsiZero, eps, epsH, gamma);

//    updateGradPsi(mesh, PsiZero, gradPsi);
//    LimitGradPsi(mesh, PsiZero, gradPsi, 1./ilePkt);

    for(int i = 1; i <= ilekrcz; ++i)
    {
        Info<<"reinitialization iter "<< i << endl;

        Told == T;

        updateGradPsi(mesh, PsiZero, gradPsi);
        LimitGradPsi(mesh, PsiZero, gradPsi, 1./ilePkt);

        //updateGradPsi(mesh, PsiZero, gradPsi);
        volScalarField k1 = createKField("k1", runTime, mesh);

        k1 =  T + fvc::div(phiR)*dtau;

        limitT(k1);
        PsiZero1(k1, PsiZero, eps, epsH, gamma);
        gradPsi = mag(fvc::grad(PsiZero));

        updateGradPsi(mesh, PsiZero, gradPsi);
        LimitGradPsi(mesh, PsiZero, gradPsi, 1./ilePkt);

        phiR = linearInterpolate( C*k1*( scalar(1.) - k1 )*( gradPsi- scalar(1.) ) * (fvc::grad(PsiZero) /( gradPsi )) ) & mesh.Sf();

        volScalarField k2  = createKField("k2", runTime, mesh);

        k2 = 3./4* T + 1./4* k1 + 1./4* fvc::div(phiR)*dtau;

        limitT(k2);
        PsiZero1(k2, PsiZero, eps, epsH, gamma);
        gradPsi = mag(fvc::grad(PsiZero));

        updateGradPsi(mesh, PsiZero, gradPsi);
        LimitGradPsi(mesh, PsiZero, gradPsi, 1./ilePkt);

        phiR = linearInterpolate( C*k2*( scalar(1.) - k2 )*( gradPsi- scalar(1.) ) * (fvc::grad(PsiZero) /( gradPsi )) ) & mesh.Sf();
        T =  1./3* T + 2./3* k2 + 2./3*fvc::div(phiR)*dtau;
        limitT(T);
        PsiZero1(T, PsiZero, eps, epsH, gamma);
        gradPsi = mag(fvc::grad(PsiZero));

        updateGradPsi(mesh, PsiZero, gradPsi);
        LimitGradPsi(mesh, PsiZero, gradPsi, 1./ilePkt);

    //    if (limitFieldT)
    //        limitT(T);

        double norm1c = Foam::sum(Foam::mag(T-Told)).value() / T.size();
        Info << "Norma 1 = " << norm1c << endl;

//        if (mapFunLog)
//            PsiZero1(T, PsiZero, eps, epsH, gamma);
//        else
//            PsiZero2(T, PsiZero, eps, epsH, gamma);

//        updateGradPsi(mesh, PsiZero, gradPsi);
//        LimitGradPsi(mesh, PsiZero, gradPsi, 1./ilePkt);
        surfgradT = linearInterpolate( T*(1-T)*(fvc::grad(PsiZero))/epsH ) & mesh.Sf();
        GradT = T*(1-T)*mag(fvc::grad(PsiZero))/epsH;

        LaplaceT = fvc::div(surfgradT);
        LaplaceTStraignt = fvc::laplacian(T);

        norm1c = Foam::sum(Foam::mag(T-Told)).value() / T.size();
        Info << "Norma 1 = " << norm1c << endl;

        magGradT = mag(fvc::grad(T));
        GradTW = T*(1.-T)*gradPsi/epsH;
        LaplaceTStraignt = fvc::laplacian(T);
        LaplaceTW1D = T*(1.-T)*( fvc::laplacian(PsiZero) + (1./epsH)*( fvc::grad(PsiZero)&fvc::grad(PsiZero) )*(1.-2.*T) )/epsH;
        gradTanalit =  mag( ( one - Foam::tanh(Psi/(2*epsH))*Foam::tanh(Psi/(2*epsH)) )*(1/(4*epsH)) );

        double norm1 = Foam::sum(Foam::mag(TAnalit-T)).value() / T.size();
        double norm2 = Foam::pow(Foam::sum(Foam::pow(TAnalit-T, 2)).value(), 0.5 ) / T.size();
        double norm3 = Foam::max(Foam::mag(TAnalit-T)).value();
//        Info << "Norma 1 analit = " << norm1 << endl;
        double norm1gr = Foam::sum(Foam::mag(gradTanalit-GradTW)).value() / T.size();
        double norm2gr = Foam::pow(Foam::sum(Foam::pow((gradTanalit-GradTW), 2)).value(), 0.5 ) / T.size();
        double norm3gr = Foam::max(Foam::mag(gradTanalit-GradTW)).value();

        double norm1lap = Foam::sum(Foam::mag(LaplaceAnalit-LaplaceTW1D)).value() / T.size();
        double norm2lap = Foam::pow(Foam::sum(Foam::pow((LaplaceAnalit-LaplaceTW1D), 2)).value(), 0.5 ) / T.size();
        double norm3lap = Foam::max(Foam::mag(LaplaceAnalit-LaplaceTW1D)).value();

        f << i << " " << norm1 << " " << norm2 << " " << norm3 << " " << norm1gr << " " << norm2gr << " " << norm3gr << " " << norm1lap <<  " " << norm2lap << " " << norm3lap << " " << norm1c <<std::endl;

        runTime.write();
        ++runTime;
    }

    f.close();

    try
    {
    f.open(sciezka + "gradAlpha_" + metcalk + "_eps_" + wielkosceps + "_" + ilePktStr + ".vtk" , std::ios::out);
    if( f.good())
        Info << "plik otwarty" << endl;
    else
        throw FileException();
    }
    catch(BasicException& ex)
    { Info << ex.wyjatek << endl; }

    double dx = 1./ilePkt;
    double norm1_grad_pkt = 0.;
    double norm1_grad_pkt_1TW = 0.;
    forAll(mesh.cellCentres(), cellI )
    {
        if ( (mesh.cellCentres()[cellI].y() > 0) && (mesh.cellCentres()[cellI].y() < dx) )
        {
            norm1_grad_pkt_1TW = Foam::mag( gradTanalit[cellI] - GradTW[cellI] ) / ( Foam::mag( gradTanalit[cellI] ) + eps );
            norm1_grad_pkt = Foam::mag( gradTanalit[cellI] - magGradT[cellI] ) / ( Foam::mag( gradTanalit[cellI] ) + eps );
            f << mesh.cellCentres()[cellI].x() << " "<< T[cellI] << " " << TAnalit[cellI] << " "  << gradTanalit[cellI] << " " << magGradT[cellI] << " " <<GradTW[cellI] << " " << LaplaceAnalit[cellI] << " " << LaplaceTW1D[cellI] << " " << LaplaceTStraignt[cellI] << " " << TSt[cellI] << " " << norm1_grad_pkt_1TW << " " << norm1_grad_pkt << std::endl;
        }
    }
    f.close();

   std::ofstream newFile(sciezka + "AnalizaZbierz.vtk", std::ios_base::app);

       if(newFile.is_open())
       {
           double norm1 = Foam::sum(Foam::mag(TAnalit-T)).value() / T.size();
           double norm2 = Foam::pow(Foam::sum(Foam::pow(TAnalit-T, 2)).value(), 0.5 ) / T.size();
           double norm3 = Foam::max(Foam::mag(TAnalit-T)).value();

           newFile << ilePktStr << " " << norm1 << " " << norm1 << " " << norm2 << " " << norm2 << " " << norm3 << " " << norm3 << std::endl;
       }
       else
       {
           //You're in trouble now Mr!
       }
       newFile.close();

}
