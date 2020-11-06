#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "param.h"
#include "uEstructuras.h"
#include "uPropFisicas.h"
#include "uPropiedades.h"
#include "uCorrelaciones.h"

#define SIGMA 5.67e-8 // cosntante de Stefan Boltzaaman, WM^-2K^-4


int main()
{
  // tabla-aire-e.c.c 
  // i7/TANQUEs-4/TANQUE-4-REPO  5 nov 2020
  // genera una tabla de propiedades del aire
  // usa las librerias  modificadas y actualizadas
  //"uPropiedades.h"
  //"uPropFisicas.h"
  // + usa las estructura: struc propiedades pro
  // se compila con libSE
FILE * log;
log = fopen("tabla-aire-e.c.log", "w");

fprintf(log,"===========tabla-aire-e.c.c ====================\n");
fprintf(log,"LUS CARDON 2020 \n");
fprintf(log," Prueba de estructura propiedades.\n");

fprintf(log,"Programa tabla-aire-e.c.c \n");
fprintf(log,"En i7/TANQUEs-4/TANQUE-4-REPO. Compilado con LIBSE. \n");
fprintf(log,"Archivo: tabla-aire-e.c.log\n");

float T,TT;
float rho,mu,nu,k,Cp,gban,alpha,Pr;

struct propiedades aire;

// for en K para comparar
fprintf(stderr,"T\trho \tmu\tk\tCp\trho/nu\tgban\talpha \n");
fprintf(stderr,"[K]\t[kg/m^3]\t[10^-6Ns/m^2]\t[10^-3W/mK]\t[J/kgK]\t[10^3s/m^2]\t[1061/m^3K]\t[10-6m^2/s] \n");
for (TT=200 ; TT<=  400 ; TT=TT+10.)
{
  // Las funciones reuieren temperaturas en grados
  
T=TT-273.15;
propiedades_aire(T, &aire );

rho=aire.rho; 
mu=aire.mu; 
nu=aire.nu; 
k=aire.k;  
Cp=aire.Cp; 
gban=aire.gban; 
alpha=aire.alpha; 
Pr=aire.Pr;
fprintf(stderr,"%3.0f\t%1.4f\t%2.2f\t%2.2f\t%4.1f\t%3.1f\t%3.1f\t%2.2f \n",TT, rho, mu, k*1.E3, Cp, 1.*1.E-3/nu ,gban*1.E-6, alpha*1.E6);  
}

return 0;
}
