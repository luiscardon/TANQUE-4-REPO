// uPropFisicas
#include "param.h"
#include "uPropFisicas.h"


double Prandlt_aire(double T)
{
  // hacer 
  return 0.707;      // 1

}
double mu_aire(double T)
{
  double aux;
// Mc Quillan, Culham and Yovanovich, 1984
// Ecuacion de Shutherland (Redid 1966)
T=T+273; 

aux= 1.4592*pow(T,3./2.)/(109.10 + T);
return aux*1.e-6;

// return 208e-7;     // Nseg/m^2

}
double nu_aire(double T)
{
 
  return rhoSnu_aire(T);     // m^2/s

}

double alpha_aire(double T)
{
// Mc Quillan, Culham and Yovanovich, 1984
// una atmosfera

return (-4.3274+4.1190e-2*T + 1.5556e-4*T*T)*1.e-6;     // m²/s

}

double rho_aire(double T)
{
// Mc Quillan, Culham and Yovanovich, 1984
// Marquardt(1963) cuadrados minimos no lineales

// T : temperatura en Celcius
// Rango 200-400K
// desvación 0.15% de datos tabulados
// una atmosfera

T=T+273; 
return   351.99/T + 344.84/(T*T); //kg/m^3
//  return 1.2;       // Kg/m^3

  
}
double k_aire(double T)
{
// Mc Quillan, Culham and Yovanovich, 1984
// Ecuacion de Shutherland (Redid 1966)
// una atmosfera
T=T+273; 
return 2.334e-3*pow(T,3./2.)/(164.54 +T); // W/mK  
//  return 0.0293;  //HACER
}
double cp_aire(double T)
{
// Mc Quillan, Culham and Yovanovich, 1984
// una atmosfera

T=T+273; 
return 1030.5 -0.19975*T + 3.9734e-4 *T*T; // J/kgK  

  
}


double Ra_aire(double T, double DT, double H)
{
  // Tm, temperatura en grados Celcius
  float g,beta,alpha, nu,Pr;
  float gbetaan;
  g=9.81;
  nu=nu_aire(T);           // HACER
  nu=1.96e-5;  // CAMBIAR
  alpha=alpha_aire(T);     // HACER
  beta=1./(T +273);
  Pr=Prandlt_aire(T);      // HACER
  Pr=0.71;                 //  CAMBIAR
  gbetaan=9.81*Pr*beta/(nu*nu);
  return gbetaan* pow(H,3.)*DT;
}
double rhoSnu_aire(double T)
{
// Mc Quillan, Culham and Yovanovich, 1984
// rho/mu
// una atmosfera
T=T+273; 

return  2.409e8/pow(T,3./2.) + 2.6737e10/pow(T,5./2.); // s/m^2
}
double gbSan_aire(double T)
{
double aux;  
// Mc Quillan, Culham and Yovanovich, 1984
// g beta /alpha nu
// una atmosfera

T=T+273; 
aux= (6.8568e-3 - 1.5079e-4*T + 1,5715e-6*T*T);
return  1.e6/(aux*aux); // m^-3K-1
}
