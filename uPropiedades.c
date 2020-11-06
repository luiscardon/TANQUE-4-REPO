#include "uEstructuras.h"
#include "uPropFisicas.h"
#include <stdio.h>
#include <math.h>

//typedef struct propiedades pro;

// uPropiedades.c
// antes propiedades.c
// Incorporada en TANQUES-4 14 de abril de 2016

// cambiar a uMateriales
// en la struct materiales (propiedades) incluir el nombre

// las correlaciones ponerlas en otrolado

// 5/nov/2020 incorporo propiedadesmde aire hecha
// originalmente en /u0-11correlaciones.c 

// El calculo de las propiedades de las substancias en funcin de la 
// temperatura se hace en uPropiedadesFisicas

float propiedades_aire(float T, pro *aire )
{
  // T, C 
  // las funciones convierten a Kelvin
  aire->rho=rho_aire(T);
  aire->Cp=Cp_aire(T);
  aire->nu=nu_aire(T);
  aire->mu=mu_aire(T);
  aire->alpha=alpha_aire(T);
  aire->k=k_aire(T);
  aire->Pr=Pr_aire(T);
  aire->gban=gban_aire(T);
  
}
float propiedades_agua(float T, pro * agua )
{
  // hacer 
  agua->rho=rho_agua(T);
  agua->Cp=Cp_agua(T);
  agua->nu=nu_agua(T);
  agua->mu=mu_agua(T);
  agua->alpha=alpha_agua(T);
  agua->k=k_agua(T); // a 40C Bejan
  agua->Pr=Pr_agua(T);
}
float print_propiedades(pro material )
{
  // imprime la propiedades termodinamicas de un material
  
fprintf(stderr,"#Propiedades de material:\n\n");
fprintf(stderr,"#Densidad              %f\n",material.rho);  
fprintf(stderr,"#Calor especifico      %f\n",material.Cp);  
fprintf(stderr,"#Viscosidad cinematica %f\n",material.nu);  
fprintf(stderr,"#Viscosidad  dinámica  %f\n",material.alpha);  
fprintf(stderr,"#Difusividad térmica   %f\n",material.k);  
fprintf(stderr,"#Numero de Prandlt     %f\n",material.Pr);  
fprintf(stderr,"#------------------------\n");  

}
/* Pasado a uPropFisicas
 * 
float Cp_agua(float T)
{ 
// calor específico del agua, J/kgK
// Qiu, pp. 21
return 4209.1 - 132.8e-2*T + 143.2e-4*T+T;
}
float k_agua(float T)
{
// conductividad térmica del agua, W/mK
return 0.5762 + 10.5e-4*T;
}
float Pr_agua(float T)
{
  // numero de Prandlt
  return 39.5345*pow(T,-0.144)-18.8396;
}
float beta_agua(float T)
{
  // coeficiente de expancion termica, K-1
  
  return (0.8*pow(T,0.5348) -1.9114)*1.e-4;
}
float nu_agua(float T)
{
  return 1.477e-6 * exp(-1.747e-2*T);
}
float rho_agua(float T)
{
  return 1000.;
}
*/
void print_propiedades_agua(float T)
{
  printf("Cp    %f    \n", Cp_agua(T) );
  printf("Pr    %f    \n", Pr_agua(T) );
  printf("k     %f    \n", k_agua(T) );
  printf("beta  %f    \n", beta_agua(T) );
  printf("nu    %f    \n", nu_agua(T) );
  printf("rho   %f    \n", rho_agua(T) );
  printf("===============================\n ");  
}
/*
  pasado a uPropFisicas 5 Nov 2020
 
//  Propiedades del Aire
// ==================================================
double Prandlt_aire(double T)
{
 // De las formulas de Mc Quillan
float alpha,nu;
T= T+273.15;
nu=nu_aire(T);
alpha=alpha_aire(T);
 return nu/alpha;      // 

}
double gban_aire(double T)
{
 // Ne la nota  Mc Quillan, Nseg/m^2
 // De las formulas de Mc Quillan
 // g beta/ alpha nu : m^-3 K^-1

float aux;
T= T+273.15;
 aux=6.8568E-3- 1.5079E-4*T + 1.5715E-6*T*T;
 aux=aux*aux;
 return 1.E6/aux;       
}

double mu_aire(double T)
{
 // Sutherland (Reid 1966). Nota de Mc Quillan, Nseg/m^2
 // Viscosidad dinamica
 T= T+273.15;
  return (1.4592*pow(T,1.5)/(109.10+T))*1.E-6;    

}
double Cp_aire(double T)
{
 // Sutherland (Reid 1966). Nota de Mc Quillan, Nseg/m^2
 // Calor específico, j/kg K
 T= T+273.15;
  return 1030.5 - 0.19975*T + 3.9743E-4*T*T;    

}
double nu_aire(double T)
{
 // Sutherland (Reid 1966). Nota de Mc Quillan, m^2/s
// Viscosidad cinematica, m^2s^-1
double nu;
T= T+273.15;
nu=2.4090E8/pow(T,1.5) + 2.6737E10/pow(T,2.5);
return  1./nu;    
}

double alpha_aire(double T)
{
 // Difusividad, m^2/s
 // Sutherland (Reid 1966). Nota de Mc Quillan
T= T+273.15;
 
return (-4.3274 + 4.119E-2* T + 1.5556E-4* T*T)* 1.E-6; 

}

double rho_aire(double T)
{
 // Densidad, Gas ideal + termino de coreciom, kg/m^3
 // Hacebfalta otra para corregirpro presion
T= T+273.15;
 
return 351.99/T + 344.84/T/T;       // Kg/m^3
  
}
double k_aire(double T)
{
T= T+273.15;
// Sutherland (Reid 1966). Nota de Mc Quillan
// Conductividad termica, W/M K

return 2.3340E-3*pow(T,1.5)/(164.54+T);
}

double Ra_aire(double T, double DT, double H)
{
  // Numero de raileigh
  // Tm, temperatura en grados Celcius
  // DT, diferencia de temperatura, K
  // H, altura, m
  
  float gban;
  T= T+273.15;
  gban=gban_aire(T);
  return gban* pow(H,3.)*DT;
}
*/
void print_propiedades_aire(float T)
{
  printf("Cp   %f   \n",Cp_aire(T) );
  printf("Pr   %f   \n",Pr_aire(T) );
  printf("k    %f   \n",k_aire(T) );
  printf("beta %f   \n",beta_aire(T) );
  printf("nu   %f   \n",nu_aire(T) );
  printf("mu   %f   \n",mu_aire(T) );
  printf("rho  %f   \n",rho_agua(T) );
 // printf("rho %f  \n",alpha_agua(T) );
  printf("gnan %f   \n",gban_aire(T) );
  
  printf("===============================\n ");  
}
// Otros materiales
// =====================================================
void p_Hormigon(pro material )
{

  material.rho=1200.; // kg/m3
  material.Cp=837.;   // J/kg C
  material.k=1.4;     // W/mC

  // no significativos
  material.nu=1.;
  material.mu=1.;
  material.alpha=1.;
  material.Pr=1.;
}

void p_BaldosaCeramica(pro material )
{
  material.rho=2000.; // kg/m3
  material.Cp=800.;   // J/kg C
  material.k=1.;     // W/mC

  // no significativos
  material.nu=1.;
  material.mu=1.;
  material.alpha=1.;
  material.Pr=1.;

}


void p_BaldosaPorcelana(pro material )
{
  material.rho=2300.; // kg/m3
  material.Cp=840.;   // J/kg C
  material.k=1.3;     // W/mC

  // no significativos
  material.nu=1.;
  material.mu=1.;
  material.alpha=1.;
  material.Pr=1.;

}



void p_BaldosaGres(pro material )
{
  material.rho=2500.; // kg/m3
  material.Cp=1000.;   // J/kg C
  material.k=2.3;     // W/mC

  // no significativos
  material.nu=1.;
  material.mu=1.;
  material.alpha=1.;
  material.Pr=1.;

}



void p_MorteroCalCemento(pro material )
{
  // Mortero de caly cemento
  material.rho=1900.; // kg/m3
  material.Cp=1000.;   // J/kg C
  material.k=0.7;     // W/mC

  // no significativos
  material.nu=1.;
  material.mu=1.;
  material.alpha=1.;
  material.Pr=1.;

}



void p_MorteroRevoque(pro material )
{
  // Mortero de revoque
  material.rho=1800.; // kg/m3
  material.Cp=0.;   // J/kg C
  material.k=1.163;     // W/mC

  // no significativos
  material.nu=1.;
  material.mu=1.;
  material.alpha=1.;
  material.Pr=1.;

}
void p_PoliestirenoExp(pro material )
{
  // Poliestireno expandido
  // a densidad 1000
  material.rho=16.; // kg/m3
  material.Cp=837.;   // J/kg C
  material.k=0.034;     // W/mC

  // no significativos
  material.nu=1.;
  material.mu=1.;
  material.alpha=1.;
  material.Pr=1.;

}
/*
void p_Hormigon(pro material )
void p_BaldosaCeramica(pro material )
void p_BaldosaPorcelana(pro material )
void p_BaldosaGres(pro material )
void p_MorteroCalCemento(pro material )
void p_MorteroRevoque(pro material )
void p_PoliestirenoExp(pro material )
*/
void p_suelo(pro * material )
{
  // propiedades del suelo
  // para en el paper de zhang20
  
  material->rho =2050.; // kg/m3
  material->Cp  =1840.;   // J/kg C
  material->k   = 0.52;     // W/mC
  material->alpha=0.0004963;

  // no significativos
  material->nu=1.;
  material->mu=1.;
  material->Pr=1.;
}
