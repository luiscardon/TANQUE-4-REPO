#include <math.h>
// Propiedades del agua en funcion de la temperatura
double Cp_agua(double T)
{ 
// calor específico del agua, J/kgK
// Qiu, pp. 21
return 4209.1 - 132.8e-2*T + 143.2e-4*T+T;
}
double k_agua(double T)
{
// conductividad térmica del agua, W/mK
return 0.5762 + 10.5e-4*T;
}
double Pr_agua(double T)
{
  // numero de Prandlt
  return 39.5345*pow(T,-0.144)-18.8396;
}
/*double Prandlt_agua(double T)
{
  // numero de Prandlt
  return 39.5345*pow(T,-0.144)-18.8396;
}*/
double beta_agua(double T)
{
  // coeficiente de expancion termica, K-1
  
  return (0.8*pow(T,0.5348) -1.9114)*1.e-4;
}
double nu_agua(double T)
{
  return 1.477e-6 * exp(-1.747e-2*T);
}
double mu_agua(double T)
{ 
  return 1.;
}
double alpha_agua(double T)
{ 
  return 1.;
}


double rho_agua(double T)
{
  return 1000.;
}

// propiedades del aire en funciond de la temperatura
// ===============================================
/*double Prandlt_aire(double T)
{
  // hacer 
  return 0.707;      // 1

}
*/
double mu_aire(double T)
{
  double aux;
// Mc Quillan, Culham and Yovanovich, 1984
// Ecuacion de Shutherland (Redid 1966)
 // Viscosidad dinamica
T=T+273; 

aux= 1.4592*pow(T,3./2.)/(109.10 + T);

}



double alpha_aire(double T)
{
 // Difusividad termica, m^2/s
 // Sutherland (Reid 1966). Nota de Mc Quillan
 // 2020
T= T+273.15;
 
return (-4.3274 + 4.119E-2* T + 1.5556E-4* T*T)* 1.E-6; 

}


double rho_aire(double T)
{
 // Densidad, Gas ideal + termino de coreciom, kg/m^3
 // Hace falta otra para corregir por presion
 // Mc Quillan, Culham and Yovanovich, 1984
 // Marquardt(1963) cuadrados minimos no lineales

// T : temperatura en Celcius
// Rango 200-400K
// desvación 0.15% de datos tabulados
// una atmosfera
// 5 nov 2020
T= T+273.15; 
return 351.99/T + 344.84/T/T;       // Kg/m^3  
}

double nu_aire(double T)
{
 // Sutherland (Reid 1966). Nota de Mc Quillan, m^2/s
// Viscosidad cinematica, m^2s^-1
 // 2020
double nu;
T= T+273.15;
nu=2.4090E8/pow(T,1.5) + 2.6737E10/pow(T,2.5);
return  1./nu;    
}

double k_aire(double T)
{
// Mc Quillan, Culham and Yovanovich, 1984
// Ecuacion de Shutherland (Redid 1966)
// una atmosfera
// Conductividad termica, W/M K
// 5 nov 2020
T=T+273.15; 
return 2.334E-3*pow(T,1.5)/(164.54 +T); // W/mK  
//  return 0.0293;  //HACER
}


double Cp_aire(double T)
{
// Mc Quillan, Culham and Yovanovich, 1984
// Sutherland (Reid 1966). Nota de Mc Quillan, Nseg/m^2
 // Calor específico, j/kg K
 // una atmosfera
// 5 nov 2020

T=T+273; 
return 1030.5 -0.19975*T + 3.9734e-4 *T*T; // J/kgK  
 
}



double Pr_aire(double T)
{
 // De las formulas de Mc Quillan
 // 5 nov 2020
 double alpha,nu;
T= T+273.15;
nu=nu_aire(T);
alpha=alpha_aire(T);
 return nu/alpha;      // 
}

//  Propiedades del Aire
// ==================================================

double gban_aire(double T)
{
 // Ne la nota  Mc Quillan, Nseg/m^2
 // De las formulas de Mc Quillan
 // g beta/ alpha nu : m^-3 K^-1
 // 2020

double aux;
T= T+273.15;
 aux=6.8568E-3- 1.5079E-4*T + 1.5715E-6*T*T;
 aux=aux*aux;
 return 1.E6/aux;       
}
double beta_aire(double T)
{
T= T+273.15;
return 1./T;
  
}



double Ra_aire(double T, double DT, double H)
{
  // Numero de raileigh
  // Tm, temperatura en grados Celcius
  // DT, diferencia de temperatura, K
  // H, altura, m
 // 2020
  
  double gban;
  T= T+273.15;
  gban=gban_aire(T);
  return gban* pow(H,3.)*DT;
}

/*
double Ra_aire(double T, double DT, double H)
{
  // Tm, temperatura en grados Celcius
 // 2020
  double g,beta,alpha, nu,Pr;
  double gbetaan;
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
*/
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

double k_pvc(double T)
{
  // pvc Policloruro de vinilo
  
  return 0.25;
  
}
