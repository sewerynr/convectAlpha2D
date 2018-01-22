#ifndef USTAWIENIA_H
#define USTAWIENIA_H


scalar one(1);  // stwórz obiekt klasy skalar o nazwie one zainicjalizowany 1
dimensionedScalar SMALL_NUMBER("small", dimless, SMALL);   //  aby nie dzielić przez zero zdefiniowana wartość "small" z przechowywana informacja o wymiarze i warosci SMALL (zdefiniowana w OF)

bool explicitSolver = runTime.controlDict().lookupOrDefault("explicit", true);
bool mapFunLog = runTime.controlDict().lookupOrDefault("mapFunLog", true);
bool limitFieldT = runTime.controlDict().lookupOrDefault("limitFieldT", false);
bool Euler = runTime.controlDict().lookupOrDefault("euler", true);

double gamma = runTime.controlDict().lookupOrDefault("gamma", 1e-6);
double eps = runTime.controlDict().lookupOrDefault("eps", 5e-16);
double k = runTime.controlDict().lookupOrDefault("k", 1.);
double par = runTime.controlDict().lookupOrDefault("par", 2.);

double pardatu = runTime.controlDict().lookupOrDefault("pardatu", 4.);

int ilKrCz = runTime.controlDict().lookupOrDefault("ilKrCz", 10);
int ilePkt = runTime.controlDict().lookupOrDefault("ilePkt", 128);
//double ParSiatki = runTime.controlDict().lookupOrDefault("ParSiatki", 1.);

#endif // USTAWIENIA_H
