// Propiedfades opticas de cubiertas
// Luis Cardon 24 de setiembre de 2013

#include "param.h"
#include "uEstructuras.h"
#include "uEstimacionRad.h"
#include "uColectores.h"
#include "uCubiertas.h"
#include "uSSS.h"
#include "uImprime.h"


// Duffie y Beckman 5.1
// Relaciones de Fresnell

float rfrac_n(float theta1, float theta2)
{
  // Luis Cardon - 24 de setiembre de 2013
  // Relación de Fresnell
  // Refracción de la luz no polarizada
  // Componente perpendicular de la radiacion polarizada
  // No funcionan  con incidencia normal, theta 1=0

  // Duffie y Beckman pp. 204, ec. 5.1.1
  
  // theta 1 : angulo de incidencia, grados
  // theta 2 : angulo de refracción, grados

  float thetamas,thetamenos,a,b;

 thetamenos=theta2-theta1;
 thetamas=theta2+theta1;

 thetamas=gradosarad(thetamas); 
 thetamenos=gradosarad(thetamenos);

 a=sin(thetamenos);
 b=sin(thetamas);

 return a*a/(b*b);
}

float rfrac_p(float theta1, float theta2)
{
  // Luis Cardon - 24 de setiembre de 2013
  // Relación de Fresnell
  // Refracción de la luz no polarizada
  // Componente paralela  de la radiacion polarizada
  // No funcionan  con incidencia normal, theta1=0

  // Duffie y Beckman pp. 204, ec. 5.1.2
  
  // theta 1 : angulo de incidencia, grados
  // theta 2 : angulo de refracción, grados

  float thetamas,thetamenos,a,b;
  
 thetamenos=theta2-theta1;
 thetamas=theta2+theta1;

 thetamas=gradosarad(thetamas); 
 thetamenos=gradosarad(thetamenos);

 a=tan(thetamenos);
 b=tan(thetamas);

 return a*a/(b*b);
}
float rfrac(float rn , float rp)
{
  // Luis Cardon - 24 de setiembre de 2013
  // Relación de Fresnell
  // Refracción de la luz no polarizada
  // Fracción reflejada de la  la radiacion polarizada
  // r= I_reflejada/I_incidente
  
  // Duffie y Beckman pp. 204, ec. 5.1.3
 return (rn+rp)/2.; 
}


// --------------------------------------------------------
// Sistemas de cubiertas iguales
// con absorcion
//

float tau_a2(float KL1, float KL2,float theta2)
{
  
  // absortancia de dos cubiertas 
  // mismo n, distinto KL

  // Ley de Bouguer
  // Duffie y Beckman, pp 208,  5.2.2
    
  // theta2 angulo de refracción, grados
  // L, espesor de la cubierta, m
  // K coeficiente de extincion, m^-1
  
  theta2=gradosarad(theta2);
// printf(" %f \n",K*L/cos(theta2));
  return exp(-(KL1+KL2)/cos(theta2));
}

// --------------------------------------------------------
// Sistemas de cubiertas iguales
// con absorcion
//

float tau_n(float rn)
{
  // Luis Cardon - 24 de setiembre de 2013
  // Transmitancia a la componente normal de la luz no polarizada
  // de una cubierta
  // Use la misma expresión para la componente paralela
  
  // Duffie y Beckman pp. 206, ec. 5.1.7
  return (1.-rn)/(1.+rn);

}

float tau_rfrac(float rfracn, float rfracp)
{
  // Luis Cardon - 24 de setiembre de 2013
  // Transmitancia  de la luz no polarizada
  // de una cubierta
  
  // Duffie y Beckman pp. 206, ec. 5.1.8
  return 0.5*(rfracn+rfracp);
}
// N cubiertas iguales

float tau_n_N(int N, float rp)
{
  // Luis Cardon - 24 de setiembre de 2013
  // Transmitancia a la componente normal de la luz no polarizada
  // de N  cubiertas iguales
  // Use la misma ecuacion para la componente paralela
  // Duffie y Beckman pp. 206, ec. 5.1.9
  
  return (1.-rp)/(1.+(2.*N-1)*rp);

}

float tau_r_N(float rpN, float rnN)
{
  // Luis Cardon - 24 de setiembre de 2013
  // Transmitancia a la componente normal de la luz no polarizada
  // de una cubierta
  // Duffie y Beckman pp. 206, ec. 5.1.9
  
  return (rpN+rnN)/2.;

}

// 
float tau_r(float taup, float taun)
{
    // Luis Cardon - 24 de setiembre de 2013
    // Transmitancia a la componente normal de la luz no polarizada
  // de una cubierta
  // Duffie y Beckman pp. 206
  return (taup+taun)/2.;

}
// --------------------------------------------------
// Cubieras con absorción
//

float tau_a(float L, float theta2, float K)
{
  
  // absortancia de un medio parcialmente transparente
  // Ley de Bouguer
  // Duffie y Beckman, pp 208,  5.2.2
    
  // theta2 angulo de refracción, grados
  // L, espesor de la cubierta, m
  // K coeficiente de extincion, m^-1
  
  theta2=gradosarad(theta2);
// printf(" %f \n",K*L/cos(theta2));
  return exp(-K*L/cos(theta2));
}

// --------------------------------------------------------
// Sistemas de cubiertas iguales
// con absorcion
//

float tau_a_n(float taua, float rn)
{
  // Luis Cardon - 24 de setiembre de 2013

  // Transmitancia  a la radiación normal de la luz no polarizada
  // de una cubierta simple con reflexion y absorcion
  
  // rn se clacula con las funciones rfrac()
  // Use la misma ecuación para la componente paralela
  // Duffie y Beckman pp. 208, ec. 5.3.1

  
  
  return taua* ( (1.-rn)/(1.+rn) )* (1.-rn*rn)/(1.-(rn*taua)*(rn*taua) );
}  

float rho_a_n(float taua, float rn, float taun)
{
    // Luis Cardon - 24 de setiembre de 2013

  // Reflectancia   a la radiación normal de la luz no polarizada
  // de una cubiert simple con reflexion y absorcion
  // Use la misma ecuación para la componente paralela
  
  // Duffie y Beckman pp. 208
  
  //  taua, absortancia de un medio parcialmente transparente 
  //  taun, transmitancia  a la radiación normal de la luz no
  //        polarizada de una cubierta simple con reflexion y
  //        absorcion, calcular con tau_a_n
  //  rn, refleccion de la luz no polarizada
  return rn * (1.+ taua*taun);
  
}  
float alpha_a_n(float taua, float rn)
{
    // Luis Cardon - 24 de setiembre de 2013

  // Absortancia  a la radiación normal de la luz no polarizada
  // de un sistema de cubiertas
  // Use la misma ecuación para la componente paralela
  
  //  taua, absortancia de un medio parcialmente transparente 
  //  rn, refleccion de la luz no polarizada
  
  return (1.-taua)* (1.- rn)/(1.-rn*taua);
  
}
float tau_a_r(float taun, float taup)
{
    // Luis Cardon - 24 de setiembre de 2013
  
    // transmitancia a la luz no polarizada
    // promedio de las componentes normal y paralela
  
    // Duffie y Beckman pp. 206, ec. 51.8

  return (taun+taup)*.5;
  
}
float rho_a_r(float rhon, float rhop)
{
    // Luis Cardon - 24 de setiembre de 2013
  
    // reflectancia  a la luz no polarizada
    // promedio de las componentes normal y paralela
  
    // Duffie y Beckman pp. 206, ec. ???
  
  return (rhon+rhop)*.5;
  
}
float alpha_a_r(float alphan, float alphap)
{
    // Luis Cardon - 24 de setiembre de 2013
  
    // transmitancia a la luz no polarizada
    // promedio de las componentes normal y paralela
  
    // Duffie y Beckman pp. 206, ec. ???  
  return (alphan+alphap)*.5;
  
}

// Sistemas de cubiertas con absorcion

float tau_s(float taun , float taup)
{
      // Luis Cardon - 24 de setiembre de 2013

  // Transmitancia  de un sistema de cubiertas

  return taun+taup;
  
}
float tau_simplif(float taua, float  taur)
{
   // Luis Cardon - 24 de setiembre de 2013

  // Reflectancia   a la radiación de la luz no polarizada
  // de una cubierta simple con reflexion y absorcion
  // ecuación simplificada  
  // Duffie y Beckman pp. 208, ec. 5.3.4
  
  //  taua, absortancia de un medio parcialmente transparente 
  //  taur, reflectancia a la luz inicailment no polarizada, ec. 5.1.8
  //        cubierta simple
  
  return taua*taur;
}
float alpha_simplif(float taua)
{
    // Luis Cardon - 24 de setiembre de 2013

  // Absortancia   a la radiación  de la luz no polarizada
  // de una cubierta simple con reflexion y absorcion
  // ecuación simplificada
  
  // Duffie y Beckman pp. 208, ec. 5.3.5
  
  //  taua, absortancia de un medio parcialmente transparente 
  
  return 1.-taua;
}
float rho_simplif(float taua,float tau)
{
    // Luis Cardon - 24 de setiembre de 2013

  // Reflectancia  a la radiación  de la luz no polarizada
  // de una cubiertA simple con reflexion y absorcion
  
  // Duffie y Beckman pp. 208, ec. 5.3.5
  
  //  taua, absortancia de un medio parcialmente transparente 
  //  tau, transmitancai de una cubierta simple,
  //       calculese con tau_simplif  
  return taua-tau;
}
// ------------------------------------------------------------
// Cubiertas dobles
//

// Cubiertas multiples 
float theta_g_eff(float beta)
{
  // Luis Cardon - 24 de setiembre de 2013
  // Angulo de incidencia efectivo para el albedo
  // DB pp. 215
  // Brandemuhel y Beckman
  return 90 - 0.5788*beta + 0.002693*beta*beta;
}
float theta_d_eff(float beta)
{
  // Luis Cardon - 24 de setiembre de 2013
  // Angulo de incidencia efectivo para radiacion difusa
  // DB pp. 215
  // Brandemuhel y Beckman
  
  return 59.7 - 0.1388*beta + 0.001497*beta*beta;
}
float ptaualpha(float tau, float alpha, float rhod)
{
  // Luis Cardon - 24 de setiembre de 2013
  // Producto transmitancia absortancia
  // DB pp. 216, ec. 5.5.1

  // alpha, absorción de la cubierta
  // tau, transmitancia de la cubierta
  // rhod, reflectancia del sistema de cubierta para la
  // la radiacion difusa incidente desde abajo
  // levemente diferente de la refelctancia difusa para la radiacion incidente
  
  return tau*alpha/(1.-(1.-alpha)*rhod);
}
float ptaualpha_s(float tau, float alpha)
{
  // Luis Cardon - 24 de setiembre de 2013
  // Producto transmitancia absortancia simplificado 
  // DB pp. 216, ec. 5.5.2 y 

  // alpha, absorción de la cubierta
  // tau, transmitancia de la cubierta
  // rhod, reflectancia del sistema de cubierta para la
  // la radiacion difusa incidente desde abajo
  // levemente diferente de la refelctancia difusa para la radiacion incidente
  
  return 1.01* tau*alpha;
}

float ta(int N, float  KL, float theta1,  float alpha, float n2)
{
  // Luis Cardon - 24 de setiembre de 2013
  // producto transmitancia absortancia
  // simplificado
  // Número de cubiertas

  float n1,stheta2,theta2;
  float taua, taun,taup,taur,tau;  
  float rn,rp,pta;
  n1=1.;

// N, numero de cubiertas
// alpha, absortancia de la placa
// KL, producto coeficiente de extincion
// n2, indice de rafracción de la placa
// n1, indice de refracción del aire
// stheta2, seno angulo de refracción 
// theta2, angulo de refracción, grados  
// rn  fraccion reflejada normal
// rp  fraccion reflejada paralela
// tau_a  transmitancia de la placa, ley de Bouguer
// tau_n  transmitancia de la cubierta, componente normal, ec. 5.1.7, 5.1.9  
// tau_p  transmitancia de la cubierta, componente normal
// tau_r  transmitancia promedio para
//        la radiacion inicialmente no polarizada, ec. 5.18
// tau transmitancia de una cubierta simple, ec.5.3.4
// pta  producto transmitancia absortancia simplificado,

// Calculo del angulo de refracción
stheta2=n1*sin(gradosarad(theta1))/n2;
theta2=asin(stheta2);
theta2=radagrados(theta2);

// fraccion reflejada, componente normal y paralela
rn=rfrac_n(theta1,theta2);
rp=rfrac_p(theta1,theta2);

// fracción transmitida luego del la absorcion, Ley de Bouguer
taua=tau_a(1.,theta2,KL);

// transmitancia 
taun= tau_n_N(N,rn);
taup= tau_n_N(N,rp);
taur=promedio(taun,taup); // 2020: hace el promedio
                       
// Transmitancia de una cubierta simple taur*taua
tau=tau_simplif(taua,taur);

// Producto transmitancia absortancia, ec 5.3.1
pta=ptaualpha_s(tau,alpha);

return pta;
}

FILE * ta_log(FILE *log, int N, float  KL, float theta1,  float alpha, float n2, float * pta)
{
  // Luis Cardon - 25 de setiembre de 2020
  // producto transmitancia absortancia
  // simplificado
  // Número de cubiertas

  float n1,stheta2,theta2;
  float taua, taun,taup,taur,tau;  
  float rn,rp;
  n1=1.;

// N, numero de cubiertas
// alpha, absortancia de la placa
// KL, producto coeficiente de extincion
// n2, indice de rafracción de la placa
// n1, indice de refracción del aire
// stheta2, seno angulo de refracción 
// theta2, angulo de refracción, grados  
// rn  fraccion reflejada normal
// rp  fraccion reflejada paralela
// tau_a  transmitancia de la placa, ley de Bouguer
// tau_n  transmitancia de la cubierta, componente normal, ec. 5.1.7, 5.1.9  
// tau_p  transmitancia de la cubierta, componente normal
// tau_r  transmitancia promedio para
//        la radiacion inicialmente no polarizada, ec. 5.18
// tau transmitancia de una cubierta simple, ec.5.3.4
// pta  producto transmitancia absortancia simplificado,

// Calculo del angulo de refracción
stheta2=n1*sin(gradosarad(theta1))/n2;
theta2=asin(stheta2);
theta2=radagrados(theta2);
// fraccion reflejada, componente normal y paralela
rn=rfrac_n(theta1,theta2);
rp=rfrac_p(theta1,theta2);

// fracción transmitida luego del la absorcion, Ley de Bouguer
taua=tau_a(1.,theta2,KL);

// transmitancia 
taun= tau_n_N(N,rn);
taup= tau_n_N(N,rp);
taur=promedio(taun,taup); // 2020: hace el promedio
                       
// Transmitancia de una cubierta simple taur*taua
tau=tau_simplif(taua,taur);

// Producto transmitancia absortancia, ec 5.3.1
*pta=ptaualpha_s(tau,alpha);
// ================================
fprintf(log,"ta_log:  theta2  :%f\n",theta2);
fprintf(log,"ta_log:  tau_a  :%f\n",taua);

fprintf(log,"ta_log:  rp  :%f\n",rp);
fprintf(log,"ta_log:  rn  :%f\n",rn);

fprintf(log,"ta_log:  taup  :%f\n",taup);
fprintf(log,"ta_log:  taun  :%f\n",taun);
fprintf(log,"ta_log:  1/2(taup+taun)  :%f\n",taur);
fprintf(log,"ta_log:  (tau alpha)_s  :%f\n",*pta);

return log;

}



float raz_alphan(float theta)
{
 // Razon de la absortancia a la absortancia normal
  // Ev 4.11.1 DV pp. 198
  // Beckman 1977
float theta2,theta3,theta4,theta5,theta6,theta7;
  theta2=theta*theta;
  theta3=theta2*theta;
  theta4=theta3*theta;
  theta5=theta4*theta;
  theta6=theta5*theta;
  theta7=theta6*theta;

  return 1. -1.5879e-3*theta + 2.7314e-4*theta2-2.3026e-5*theta3+
  9.0244e-7*theta4 -1.8e-8*theta5 + 1.7734e-10*theta6 -6.9937e-13 *theta7; 
}
float tau_theta(int N, float  KL, float theta1,  float alpha, float n2)
{
  
  float n1,stheta2,theta2;
  float taua, taun,taup,taur,tau;  
  float rn,rp,pta;
  n1=1.;
  
stheta2=n1*sin(gradosarad(theta1))/n2;
theta2=asin(stheta2);
theta2=radagrados(theta2);
rn=rfrac_n(theta1,theta2);
rp=rfrac_p(theta1,theta2);
taua=tau_a(1.,theta2,KL);
taun= tau_n_N(N,rn);
taup= tau_n_N(N,rp);
taur=tau_r(taun,taup);
tau=tau_simplif(taua,taur);
return tau;
}
float taos(int N,float KL, float theta,float beta,float alphann, float n2,
	    float * ptab,float *ptad,float * ptag)
{
  // Calcula el producto transmitancia absortancia
  // Determinacion individual de la dependencia 
  // angular de tao y de alpha
  
  // Energia solar absorbida por el colector
  // Ejemplo 5.9.1 DB pp. 222
  // Procedimiento 1
  
// Verificado
// ptab requiere entrada como dobles
/*
KL=0.037;
theta=17.;
beta=60.;
alphann=0.93;
*/

float alpha, theta1ed,  alphad;
// alpha, absortancia de la placa
// alphann, absortancia de la placa a incidencia normal
// theta, angulo de incidencia, grados
// ptab, producto transmitancia absortancia directa
// ptad, producto transmitancia absortancia difusa
// ptag, producto transmitancia absortancia albedo
// theta1ed, auxiliar
// N, numero de cubiertas
// KL, KL coeficiente de extincion
// n2, indice de refraccion del material de cubierta

// llama a: 
// raz_alphan(theta): razon alpha/alphan 
// ta(N, KL, theta,  alpha, n2); tao vs theta para cubiertas, uCubiertas.c
// theta_d_eff(beta), angulo de incidencia efectivo, uCubiertas.c
// theta_g_eff(beta), angulo de incidencia efectivo, uCubiertas.c

// calculo del producto transmitancia absortancia paralela
// la radiacion directa

alpha=alphann*raz_alphan(theta);
*ptab =ta(N, KL, theta,  alpha, n2); // como evitar n2 si tengo KL?

// calculo del producto transmitancia absortancia paralela
// la radiacion difusa

theta1ed=theta_d_eff(beta);
alphad=alphann*raz_alphan(theta1ed);
*ptad =ta(N, KL, theta1ed,  alphad,  n2); 

// calculo del producto transmitancia absortancia paralela
// la radiacion reflectada por el suelo, albedo

theta1ed=theta_g_eff(beta);
alphad=alphann*raz_alphan(theta1ed);
*ptag =ta(N, KL, theta1ed,  alphad, n2); 

return 0.;
}


// ========================================================

void  TyR_nn(int N, float tao, float rho, float * Tn, float * Rn)
{
int i;
float T1,R1;
// transmitancia de N cubiertas con absorcion
// formula recursiva
// usa las formulas de recurrencia DB ec. 5.37 y 5.38
//

// Luis Cardon - 30 de marzo de 2015

// datos
//     tao  : transmitancia del panel simple
//     rho  : reflectancia del panel simple
//     N    : numero de cubiertas

// Devuelve
//     Tn   : transmitancia de la cubierta
//     Rn   : reflectancai de la cubierta
// Agregar
//     An   : abortancia de la cubierta


T1=tao;
R1=rho;

*Tn=T1;
*Rn=R1;
if(N==1)
{  *Rn =  R1 + R1*T1*T1/(1.- *Rn*R1); 
   *Tn = T1* *Tn/(1.-R1* *Rn);
}else
{for(i=2;i<=N;i++)
 { *Rn=  R1 + R1*T1*T1/(1.- *Rn*R1); 
   *Tn = T1* *Tn/(1.-R1* *Rn);
 }  }
   return;
}

//===================================================
float Tmn(float Tm,float Tn, float Rm, float Rn)
{
// transmitancia de una cubierta doble de propiedades distintas
// con absorcion
// DB. pp 211 Ec. 5.3.7, 5.3.8 
// Luis Cardon - 30 de marzo de 2015

// datos
//     Tm   : transmitancia de la cubierta m
//     Rm   : reflectancia de la cubierta m
//     Tn   : transmitancia de la cubierta n
//     Rn   : reflectancia de la cubierta n

// Devuelve
//     : transmitancia de la cubierta

  return Tm*Tn/(1-Rm*Rn);
}
//==================================================
float Rmn( float Tm,float Tn, float Rm, float Rn)
{
// Reflectancia de una cubierta doble de propiedades distintas
// con absorcion
// DB. pp 211 Ec. 5.3.7, 5.3.8 
// Luis Cardon - 30 de marzo de 2015

// datos
//     T    : transmitancia de la cubierta doble. 
//     T    : debe calcularse antes que R

//     Tm   : transmitancia de la cubierta m
//     Rm   : reflectancia de la cubierta m
//     Tn   : transmitancia de la cubierta n
//     Rn   : reflectancia de la cubierta n

// Devuelve
//     : reflectancia  de la cubierta
return Rm + Rn*Tm*Tm/(1-Rm*Rn);
}
//=====================================================
float TyR_nn_2(int N, float n,float K, float L, float theta1, float *T, float *R)
{
// transmitancia y reflectancia de cubiertas absorbentes multiples
// de paneles iguales
// usa las formulas de recurrencia DB ec. 5.37 y 5.38
// Verificado con el problema 5.3.2 DB. pp. 210
// tau_DB=0.69 (formulas aproximadas) tau=0.67239

// Luis Cardon - 30 de marzo de 2015

// datos:
//       N numero de cubiertas
//       n indice  refraccion
//       K coeficientes de extincion
//       L  espesor de la cubierta
// devuelve:
//       Tn transmitancia
//       Rn reflectancia

float tau, rho;
float  theta1rad, sintheta1, theta2,sintheta2, theta2rad;
float rp,rn,taup,taun,taur,taua;
float rhon,rhop;
float Tn,Rn,Tp,Rp;
//Tn2=*Tn;
//Rn2=*Rn;
// Ver pp 210
// Cubierta simple, 60 grados,  ej. 5.3.1

// Angulo de refraccion
theta1rad= gradosarad(theta1);
sintheta2= sin(theta1rad)/n;
//printf("sintheta2 %f \n",sintheta2);
theta2rad=asin(sintheta2);
theta2= radagrados(theta2rad);
//printf("theta2 %f \n",theta2);
//
taua=tau_a(L,theta2,K);
// Componentes de polarizacion de la reflectancias
// para el panel simple
rp=rfrac_p(theta1,theta2);
rn=rfrac_n(theta1,theta2);
taup=tau_a_n(taua, rp);
taun=tau_a_n(taua, rn);
//tau=.5*(taup+taun);

rhon=rho_a_n(taua,rn,taun); 
rhop=rho_a_n(taua,rp,taup); 
//rho=0.5*(rhon+rhop);

TyR_nn(N,taup,rhop,&Tp,&Rp);
TyR_nn(N,taun,rhon,&Tn,&Rn);

//printf("%f  %f\n",Tn2, Rn2 );
*T=0.5*(Tn+Tp);
*R=0.5*(Rn+Rp);
}
//===================================================

float TyR_CC(int N, float nn,float Kn, float Ln, 
int M, float nm,float Km, float Lm,float theta1, float *T, float *R)
{
// transmitancia y reflectancia de cubiertas 
// absorbentes multiples de paneles diferentes
// usa las formulas de recurrencia DB ec. 5.37 y 5.38
// Verificado con el problema 5.3.2 DB. pp. 210
// tau_DB=0.69 (formulas aproximadas) tau=0.67239

// Luis Cardon - 30 de marzo de 2015

// datos:
//       N numero de cubiertas
//       n indice  refraccion
//       K coeficientes de extincion
//       L  espesor de la cubierta
// devuelve:
//       Tn transmitancia
//       Rn reflectancia

//       igual para M

//float taun, rhon, taum,rhom;

//fprintf(*log,"Datos\n");
//fprintf(*log,"N  :%d\n",N);
//fprintf(*log,"n  :%d\n",n);
//fprintf(*log,"K  :%f\n",K);
//fprintf(*log,"L  :%F\n",L);
//fprintf(*log,"-----------\n");
//fprintf(*log,"M  :%d\n",N);
//fprintf(*log,"n  :%d\n",n);
//fprintf(*log,"K  :%f\n",K);
//fprintf(*log,"L  :%F\n",L);
//fprintf(*log,"-----------\n");
//fprintf(*log,"Theta1  :%f\n",theta1);



float  theta1rad, sintheta1, theta2n, theta2, sintheta2, theta2rad;
//fprintf(*log,"Cubierta de N paneles:\n",);

float rp,rn,taup,taun,taur,taua;
float rhon,rhop;
float TNn,TNp, RNn,RNp ;
float TMn,TMp, RMn,RMp ;

// Calculo de la cubierta de M paneles iguales

// Angulo de refraccion
theta1rad= gradosarad(theta1);
sintheta2= sin(theta1rad)/nn;
//printf("sintheta2 %f \n",sintheta2);
theta2rad=asin(sintheta2);
theta2= radagrados(theta2rad);
//fprintf(*log,"Theta2  :%f\n",theta2);

//printf("theta2 %f \n",theta2);
//
taua=tau_a(Ln,theta2,Kn);
//fprintf(*log,"tau_a  :%f\n",taua);

// Componentes de polarizacion de la reflectancias
// para el panel simple
rp=rfrac_p(theta1,theta2);
rn=rfrac_n(theta1,theta2);
//fprintf(*log,"rp  :%f\n",rp);
//fprintf(*log,"rn  :%f\n",rn);

taup=tau_a_n(taua, rp);
taun=tau_a_n(taua, rn);
//fprintf(*log,"taup  :%f\n",taup);
//fprintf(*log,"taun  :%f\n",taun);

//tau=.5*(taup+taun);

rhon=rho_a_n(taua,rn,taun); 
rhop=rho_a_n(taua,rp,taup); 
//rho=0.5*(rhon+rhop);
//fprintf(*log,"rhop  :%f\n",rhop);
//fprintf(*log,"rhon  :%f\n",rhon);

TyR_nn(N,taup,rhop,&TNp,&RNp);
TyR_nn(N,taun,rhon,&TNn,&RNn);
printf("cubierta N\n");
printf("%f  %f\n",TNp, RNp );
printf("%f  %f\n",TNn, RNn );

//fprintf(*log,"TNp  RNp  :%f  %f\n",TNp,RNp);
//fprintf(*log,"TNn  RNn  :%f  %f\n",TNn,RNn);

// Calculo de la cubierta de M paneles iguales
//fprintf(*log,"Cubierta de N paneles:\n",);

// Angulo de refraccion
theta1rad= gradosarad(theta1);
sintheta2= sin(theta1rad)/nm;
//printf("sintheta2 %f \n",sintheta2);
theta2rad=asin(sintheta2);
theta2= radagrados(theta2rad);
//fprintf(*log,"Theta2  :%f\n",theta2);

//printf("theta2 %f \n",theta2);
//
taua=tau_a(Lm,theta2,Km);
//fprintf(*log,"tau_a  :%f\n",taua);
// Componentes de polarizacion de la reflectancias
// para el panel simple
rp=rfrac_p(theta1,theta2);
rn=rfrac_n(theta1,theta2);
//fprintf(*log,"rp  :%f\n",rp);
//fprintf(*log,"rn  :%f\n",rn);
taup=tau_a_n(taua, rp);
taun=tau_a_n(taua, rn);
//tau=.5*(taup+taun);
//fprintf(*log,"taup  :%f\n",taup);
//fprintf(*log,"taun  :%f\n",taun);

rhon=rho_a_n(taua,rn,taun); 
rhop=rho_a_n(taua,rp,taup); 
//rho=0.5*(rhon+rhop);
//fprintf(*log,"rhop  :%f\n",rhop);
//fprintf(*log,"rhon  :%f\n",rhon);

TyR_nn(M,taup,rhop,&TMp,&RMp);
TyR_nn(M,taun,rhon,&TMn,&RMn);
printf("cubierta M\n");

printf("%f  %f\n",TMp, RMp );
printf("%f  %f\n",TMn, RMn );
//fprintf(*log,"TNp  RNp  :%f  %f\n",TNp,RNp);
//fprintf(*log,"TNn  RNn  :%f  %f\n",TNn,RNn);


*T=0.5*( Tmn(TMp,TNp,RMp,RNp) + Tmn(TMn,TNn,RMn,RNn));
*R=0.5*( Rmn(TMp,TNp,RMp,RNp) + Rmn(TMn,TNn,RMn,RNn));

//printf(*log,"T   N %f  %f\n",*T, *R );

return;  
}
FILE * TyR_CC_log(FILE *log, int N, float nn,float Kn, float Ln, 
int M, float nm,float Km, float Lm,float theta1, float *T, float *R)
{
// Transmitancia y reflectancia de cubiertas absorbentes multiples
// de paneles diferentes
// usa las formulas de recurrencia DB ec. 5.37 y 5.38
// Verificado con el problema 5.3.2 DB. pp. 210
// tau_DB=0.69 (formulas aproximadas) tau=0.67239

// Luis Cardon - 30 de marzo de 2015

// datos:
//       N numero de cubiertas
//       n indice  refraccion
//       K coeficientes de extincion
//       L  espesor de la cubierta
// devuelve:
//       Tn transmitancia
//       Rn reflectancia

//       igual para M

//float taun, rhon, taum,rhom;
fprintf(log,"TyR_CC_log:  Calculo de la transmitancia T y Reflectancia de la cubierta:\n");
fprintf(log,"TyR_CC_log:  Datos\n");
fprintf(log,"TyR_CC_log:  N   %d\n",N);
fprintf(log,"TyR_CC_log:  n  :%f\n",nn);
fprintf(log,"TyR_CC_log:  K  :%f\n",Kn);
fprintf(log,"TyR_CC_log:  L  :%f\n",Ln);
fprintf(log,"TyR_CC_log:  KL  :%f\n",Ln*Kn);

fprintf(log,"TyR_CC_log:  -----------\n");
fprintf(log,"TyR_CC_log:  M  :%d\n",N);
fprintf(log,"TyR_CC_log:  n  :%f\n",nm);
fprintf(log,"TyR_CC_log:  K  :%f\n",Km);
fprintf(log,"TyR_CC_log:  L  :%f\n",Lm);
fprintf(log,"TyR_CC_log:  KL  :%f\n",Lm*Km);
fprintf(log,"TyR_CC_log:  -----------------------------\n");
fprintf(log,"TyR_CC_log:  theta1  :%f\n",theta1);



float  theta1rad, sintheta1, theta2n, theta2, sintheta2, theta2rad;
fprintf(log,"TyR_CC_log:  -----------------------------\n");
fprintf(log,"TyR_CC_log:  Cubierta de N=%d paneles:\n",N);

float rp,rn,taup,taun,taur,taua;
float rhon,rhop;
float TNn,TNp, RNn,RNp ;
float TMn,TMp, RMn,RMp ;

// Calculo de la cubierta de M paneles iguales

// Angulo de refraccion
theta1rad= gradosarad(theta1);
sintheta2= sin(theta1rad)/nn;
//printf("sintheta2 %f \n",sintheta2);
theta2rad=asin(sintheta2);
theta2= radagrados(theta2rad);
fprintf(log,"TyR_CC_log:  theta2  :%f\n",theta2);

//printf("theta2 %f \n",theta2);
//
taua=tau_a(Ln,theta2,Kn);
fprintf(log,"TyR_CC_log:  tau_a  :%f\n",taua);

// Componentes de polarizacion de la reflectancias
// para el panel simple
rp=rfrac_p(theta1,theta2);
rn=rfrac_n(theta1,theta2);
fprintf(log,"TyR_CC_log:  rp  :%f\n",rp);
fprintf(log,"TyR_CC_log:  rn  :%f\n",rn);

taup=tau_a_n(taua, rp);
taun=tau_a_n(taua, rn);
fprintf(log,"TyR_CC_log:  taup  :%f\n",taup);
fprintf(log,"TyR_CC_log:  taun  :%f\n",taun);

//tau=.5*(taup+taun);

rhon=rho_a_n(taua,rn,taun); 
rhop=rho_a_n(taua,rp,taup); 
//rho=0.5*(rhon+rhop);
fprintf(log,"TyR_CC_log:  rhop  :%f\n",rhop);
fprintf(log,"TyR_CC_log:  rhon  :%f\n",rhon);



TyR_nn(N,taup,rhop,&TNp,&RNp);
TyR_nn(N,taun,rhon,&TNn,&RNn);
printf("TyR_CC_log:  cubierta N\n");
printf("TyR_CC_log:  %f  %f\n",TNp, RNp );
printf("TyR_CC_log:  %f  %f\n",TNn, RNn );

fprintf(log,"TyR_CC_log:  TNp  RNp  :%f  %f\n",TNp,RNp);
fprintf(log,"TyR_CC_log:  TNn  RNn  :%f  %f\n",TNn,RNn);

fprintf(log,"TyR_CC_log:  TN  :%f  \n",.5*(TNp+TNn));
fprintf(log,"TyR_CC_log:  RN  :%f  \n",.5*(RNn+RNp));


// Calculo de la cubierta de M paneles iguales
fprintf(log,"TyR_CC_log:  -----------------------------\n");
fprintf(log,"TyR_CC_log:  Cubierta de M=%d paneles:\n",M);

// Angulo de refraccion
theta1rad= gradosarad(theta1);
sintheta2= sin(theta1rad)/nm;
//printf("sintheta2 %f \n",sintheta2);
theta2rad=asin(sintheta2);
theta2= radagrados(theta2rad);
fprintf(log,"TyR_CC_log:  theta2  :%f\n",theta2);

//printf("theta2 %f \n",theta2);
//
taua=tau_a(Lm,theta2,Km);
fprintf(log,"TyR_CC_log:  tau_a  :%f\n",taua);
// Componentes de polarizacion de la reflectancias
// para el panel simple
rp=rfrac_p(theta1,theta2);
rn=rfrac_n(theta1,theta2);
fprintf(log,"TyR_CC_log:  rp  :%f\n",rp);
fprintf(log,"TyR_CC_log:  rn  :%f\n",rn);
taup=tau_a_n(taua, rp);
taun=tau_a_n(taua, rn);
//tau=.5*(taup+taun);
fprintf(log,"TyR_CC_log:  taup  :%f\n",taup);
fprintf(log,"TyR_CC_log:  taun  :%f\n",taun);

rhon=rho_a_n(taua,rn,taun); 
rhop=rho_a_n(taua,rp,taup); 
//rho=0.5*(rhon+rhop);
fprintf(log,"TyR_CC_log:  rhop  :%f\n",rhop);
fprintf(log,"TyR_CC_log:  rhon  :%f\n",rhon);

TyR_nn(M,taup,rhop,&TMp,&RMp);
TyR_nn(M,taun,rhon,&TMn,&RMn);
printf("TyR_CC_log:  cubierta M\n");

printf("TyR_CC_log:  %f  %f\n",TMp, RMp );
printf("TyR_CC_log:  %f  %f\n",TMn, RMn );
fprintf(log,"TyR_CC_log:  TNp  RNp  :%f  %f\n",TNp,RNp);
fprintf(log,"TyR_CC_log:  TNn  RNn  :%f  %f\n",TNn,RNn);

fprintf(log,"TyR_CC_log:  TN  :%f  \n",.5*(TNp+TNn));
fprintf(log,"TyR_CC_log:  RN  :%f  \n",.5*(RNn+RNp));

fprintf(log,"TyR_CC_log: ---------------------\n \n");

*T=0.5*(Tmn(TMp,TNp,RMp,RNp)+Tmn(TMn,TNn,RMn,RNn));
*R=0.5*(Rmn(TMp,TNp,RMp,RNp)+Rmn(TMn,TNn,RMn,RNn));

fprintf(log,"TyR_CC_log:  T   R     :%f  %f\n\n", *T, *R );
fprintf(log,"TyR_CC_log:  fin TyR_CC_log ---------------------\n \n");
return log;
  
}

// ========================================================
float tau_cs(int N, float theta1, float n )
{
// tau_cs  - Luis Cardon - 27 marzo 2015
// Calcula la transmitancia de N cubiertas iguales sin absorcion
// Comprobado con Ej. 5.1.2 DB3. pp. 2017
// theta1 60, n= 1.526, N=2, tau 0.758780 vs 0.76

//  Datos:
//       N    : numero de cubiertas
//       thata: angulo de incidencia en grados
//       n    : indice de refraccion del panel

float   theta1rad, sintheta1, theta2,sintheta2, theta2rad;
float   rp,rn,taupN,taunN,taur;

// Angulo de refraccion
theta1rad= gradosarad(theta1);
sintheta2= sin(theta1rad)/n;
//printf("sintheta2 %f \n",sintheta2);
theta2rad=asin(sintheta2);
theta2= radagrados(theta2rad);
//printf("theta2 %f \n",theta2);

// Componentes de polarizacion de la reflectancias
// para el panel simple
rp=rfrac_p(theta1,theta2);
rn=rfrac_n(theta1,theta2);
// Componentes de polarizacon para la
// transmitancias de cunierta de N paneles 
taupN=tau_n_N(N,rp);
taunN=tau_n_N(N,rn);
// Transmitancias de cunierta de N paneles 
taur=tau_r_N(taupN,taunN);
return taur;  
}

float theta_2(float n1, float n2, float theta1)
{
float  theta1rad, sintheta1,  sintheta2, theta2rad, theta2;

theta1rad= gradosarad(theta1);
sintheta2= n1*sin(theta1rad)/n2;
theta2rad=asin(sintheta2);
//printf("sintheta2 %f \n",sintheta2);
theta2= radagrados(theta2rad);
return theta2;

}
float alpha_alphan(float theta, float alphan)
{
  // Calcula la absortancia en funcion del angulo de incidencia
  // y de la absortancia a incidencia normal
  // Correlacion de Duffie y Beckman, pp. 198
  
  // datos
  // theta, angulo de incidencia, grados
  // alphan, absortancia a incidencia normal
  
float ratio, t,t2,t4;
t=theta;
t2=t*t;
t4=t2*t2;
ratio= 1.- 1.5879e-3 * t  +
           2.7314e-4 * t2 -
           2.3026e-5 * t2*t +
           9.0244e-7 * t4 -
           1.8000e-8 * t4*t +
           1.7734e-10*t4*t2 -
           6.9937e-13*t4*t2*t;

return ratio*alphan;
}
// 
float promedio(float a, float b)
{
    // Luis Cardon - 25 de setiembre de 2020
  return (a+b)/2.;

}
