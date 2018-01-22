/*---------------------------------------------------------------------------*\
| OpenFOAM: The Open Source CFD Toolbox
\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include <fstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "funkcje.h"
#include "exeuler.h"
#include "exrk4.h"
#include "implicit.h"
#include "fun_inicjaliz.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#include "setRootCase.H"   // sprawdza katalog kejsa
#include "createTime.H"
#include "createMesh.H"
    simpleControl simple(mesh);
#include "createFields.H"

//    InitUXYZ(U, mesh, funInitV2); //2D
//    InitUXYZ(U, mesh, funInitV1); //1D

#include "createPhi.H"

#include "createFvOptions.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CourantNo.H"
#include "ustawienia.h"

    dimensionedScalar epsH( "epsH", dimLength, Foam::pow(Foam::average(mesh.cellVolumes()), one/3.) /2. );
    epsH = epsH*Foam::sqrt(k);
//    dimensionedScalar dtau( "dTau", explicitSolver ? dimLength : dimLength/dimTime , epsH.value() / pardatu );
    dimensionedScalar dtau( "dTau", explicitSolver ? dimTime : dimless, epsH.value() / pardatu );

//    InitPsiXYZ(Psi, mesh, funInitPsi2); //2D
    InitPsiXYZ(Psi, mesh, funInitPsi1); //1D
    //InitT(T, Psi, mesh, epsH, funInitT3, par);      //bez inicjalizacji brzegu, Gaussian kropla
//    InitT(T, Psi, mesh, epsH, funInitT4, par);      // inijalizacja alpha(psi) 2D
    InitT(T, Psi, mesh, epsH, funInitT1, par);      //1D
    //    InitTB(T, Psi, mesh, epsH, funInitT3, par);
    runTime.write();
    Info << "dx =  " << Foam::pow(Foam::average(mesh.cellVolumes()), one/3.) << endl;
    Info << "epsH =  " << epsH.value() << endl;
    Info << "dtau =  " << dtau.value() << endl;
    Info << "ilosc Pkt =  " << ilePkt << endl;

    bool flipped = false;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        Info<< "while simple loop = " << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            Info<< "while simple loop ort cor = " << nl << endl;

            fvScalarMatrix TEqn
            (
                fvm::ddt(T)
                + fvm::div(phi, T)
            );

            TEqn.relax();
            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(T);
        }

        if(explicitSolver)
        {
            if(Euler)
            {
                exEuler(C, runTime, mesh, dtau, PsiZero, PsiZero, T, Told, epsH, eps, limitFieldT, ilKrCz, gamma, mapFunLog);
            }
            else
            {
                exRK3(C, runTime, mesh, dtau, Psi, PsiZero, T, Told, epsH, eps, limitFieldT, ilKrCz, gamma, mapFunLog, ilePkt);
            }
        }
        else
        {
            implicit(C, runTime, mesh, dtau, PsiZero, PsiZero, T, Told, epsH, eps, limitFieldT, ilKrCz, gamma, mapFunLog);
        }

        ++runTime;
        runTime.write();

//        zawroc przeplyw !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//        if (runTime.value() > runTime.endTime().value()/2 && !flipped)
//        {
//            InitUXYZ(U, mesh, funInitV3);
//            phi = fvc::flux(U);
//            flipped != flipped;
//        }

    }

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Only reinitialization !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//            if(explicitSolver)
//            {
//                if(Euler)
//                {
//                    exEuler(C, runTime, mesh, dtau, PsiZero, PsiZero, T, Told, epsH, eps, limitFieldT, ilKrCz, gamma, mapFunLog);
//                }
//                else
//                {
//                    exRK3(C, runTime, mesh, dtau, Psi, PsiZero, T, Told, epsH, eps, limitFieldT, ilKrCz, gamma, mapFunLog, ilePkt);
//                }
//            }
//            else
//            {
//                implicit(C, runTime, mesh, dtau, PsiZero, PsiZero, T, Told, epsH, eps, limitFieldT, ilKrCz, gamma, mapFunLog);
//            }

    Info<< "End\n" << endl;
    return 0;
}
