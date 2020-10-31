// UTILITARIOS - uEstimacionRad.c

// Luis Cardon 20 de setiembre de 2013
// derivado de fpp.c




#include "param.h"
#include "uEstructuras.h"
#include "uEstimacionRad.h"
#include "uColectores.h"
#include "uCubiertas.h"
#include "uSSS.h"
#include "uImprime.h"
#include <math.h>

// ------------------------------------------------------------------
// Calculo de la irradiancia  instantanea, W/m^2
float Gon(int n)
{
// Duffie y Beckman, 1980, pp  22.
// Duffie y Beckman, 2006, pp  9, ec. 1.4.1a
// Retorna radiacion extraterrestre sobre plano normal a la radiacion
// Datos
// n: día del año
// Sale:
// Gon radiacion extraterrestre, W/m^2


// ang angulo en grados
// Gon radiacion extraterrestre, W/m^2

float g0,ang,angrad;
ang=360.* n / 365.;             // DyB angulo en grados
angrad= ang *  PI/ 180.;
return GSC * ( 1. + 0.033 * cos(angrad) );
}
// ------------------------------------------------------------------
float ooGon(sgeo geo)
{
// Duffie y Beckman, 1980, pp  22.
// Duffie y Beckman, 2006, pp  9, ec. 1.4.1a
// Retorna radiacion extraterrestre sobre plano normal a la radiacion

// Nececita: n día del año

// Datos:
// geo: estructura de datos geograficos

// ang angulo en grados
// Gon radiacion extraterrestre, W/m^2

float g0,ang,angrad;
int n;
n=geo.n;

angrad= gradosarad(360.* n / 365.);
return GSC * ( 1. + 0.033 * cos(angrad) );
}
// ------------------------------------------------------------------

float Go(int n, float thetaz)
{
 // Luis Cardon - 20 de setiembre de 2013
  
// Duffie y Beckman, 1980, pp  22.
// Duffie y Beckman, 2006, pp 37, Ec. 1.10.1
// Retorna radiacion extraterrestre sobre plano horizontal
// Entra n:  día del año
// Entra thetaz: angulo cenital
//
  thetaz=gradosarad(thetaz);
  return Gon(n)*cos(thetaz);
}

// ------------------------------------------------------------------
// Caculos de la irradiacion horaria
float I_0(int n, float phi, float delta, float omega1, float omega2)
{
// Retorna I0, radiacion extraterrestre sobre superficie horizontal
// entre   omega_1 y omega_2
// Ecuacion  1.10.4, DB pp. 40 
// Comprobado ok con ejemplo 1.10.2 pp. 40 DyB 
  
float sinphi, cosphi, sindelta, cosdelta, sinomega1, sinomega2 ;
float phirad, deltarad,  omegarad1,omegarad2;
float angulo;
float a, domega;

domega=omega2-omega1;
domega=gradosarad(domega);

a=12.*3600.*GSC/PI;

phirad=gradosarad(phi);
deltarad=gradosarad(delta);
omegarad1=gradosarad(omega1);
omegarad2=gradosarad(omega2);
angulo=gradosarad(360.*n/365.);

sinphi=sin(phirad);
cosphi=cos(phirad);
sindelta=asin(deltarad);
cosdelta=cos(deltarad);
sinomega1=sin(omegarad1);
sinomega2=sin(omegarad2);


return a*(1.+ 0.033 *cos(angulo))*
(cosphi*cosdelta *(sinomega2-sinomega1)+ domega*sinphi*sindelta);
}
// Caculos de la irradiacion horaria
float I_0omega12(sgeo geo, float omega1, float omega2)
{
// Luis Cardon - 26 de mayo 2014
// de I_0, modifico para usar estructuras.

// Retorna I0, radiacion extraterrestre sobre superficie horizontal
// entre   omega_1 y omega_2
// Ecuacion  1.10.4, DB pp. 40 
// Comprobado ok con ejemplo 1.10.2 pp. 40 DyB 
  
float sinphi, cosphi, sindelta, cosdelta, sinomega1, sinomega2 ;
float phirad, deltarad,  omegarad1,omegarad2;
float angulo;
float a, domega;
float phi,beta,delta,gamma,hr;
int n;
phi=geo.phi;           // latitud
n=geo.n;               // febrero 20
beta=geo.beta;         // inclinacion de la superficie
delta=geo.delta;       // 11.57, ok
gamma=geo.gamma;


domega=omega2-omega1;
domega=gradosarad(domega);

a=12.*3600.*GSC/PI;

phirad=gradosarad(phi);
deltarad=gradosarad(delta);
omegarad1=gradosarad(omega1);
omegarad2=gradosarad(omega2);
angulo=gradosarad(360.*n/365.);

sinphi=sin(phirad);
cosphi=cos(phirad);
sindelta=sin(deltarad);
cosdelta=cos(deltarad);
sinomega1=sin(omegarad1);
sinomega2=sin(omegarad2);


return a*(1.+ 0.033 *cos(angulo))*
(cosphi*cosdelta *(sinomega2-sinomega1)+ domega*sinphi*sindelta);
}

float ooI0(sgeo geo, sgeohr geohr)
{
// Luis Cardon - 26 de mayo 2014
// de I_0, modifico para usar estructuras.
// dato: hr, calcula entre hr y hr+1

// Retorna I0, radiacion extraterrestre sobre superficie horizontal
// para hr + 1/2 
// calcula   omega_1 y omega_2 correspondiente hr y hr+1
// Ecuacion  1.10.4, DB pp. 40 
// Comprobado ok con ejemplo 1.10.2 pp. 40 DyB 
// Comprobado de nuevo 3.79MJ vs. 3795434  
float sinphi, cosphi, sindelta, cosdelta, sinomega1, sinomega2 ;
float phirad, deltarad,  omegarad1,omegarad2;
float angulo;
float a, domega;
float phi,beta,delta,gamma,hr;
float omega1,omega2;
int n;
phi=geo.phi;           // latitud
n=geo.n;               // febrero 20
beta=geo.beta;         // inclinacion de la superficie
delta=geo.delta;       // 11.57, ok
gamma=geo.gamma;

hr=geohr.hr;            // hora
omega1=anguloh_h(hr);
omega2=anguloh_h(hr+1);

// hay que calular omega1 y omega2 correspondiente

domega=omega2-omega1; // dirferencia en grados
domega=gradosarad(domega);

a=12.*3600.*GSC/PI;

phirad=gradosarad(phi);
deltarad=gradosarad(delta);
omegarad1=gradosarad(omega1);
omegarad2=gradosarad(omega2);
angulo=gradosarad(360.*n/365.);

sinphi=sin(phirad);
cosphi=cos(phirad);
sindelta=sin(deltarad);
cosdelta=cos(deltarad);
sinomega1=sin(omegarad1);
sinomega2=sin(omegarad2);

//return a*(1.+ 0.033 *cos(angulo))*
//(cosphi*sindelta *(sinomega2-sinomega1)+ domega*sinphi*sindelta);
// en DB hay un error es con sin o cos pp.40
// con cos reproduzco los calculos
return a*(1.+ 0.033 *cos(angulo))*
(cosphi*cosdelta *(sinomega2-sinomega1)+ domega*sinphi*sindelta);
}

float I_total(int n, float phi, float delta, float I, 
	     float beta, float gamma, float rhog, float hr,
	    float * Ibb, float * Idd, float * Igg )
{
  // I_total  discriminado
  // I_T Radiacion total sobre superficie inclinada
  // Luis Cardon  23 de setiembre de 2013

  // archivo: usEstimacionRad.c
  
// Radiacion sobre superficie inclinada: cielo isotropico
// en el per铆odo de una hora, +- 30 minutos de hr
  
// Duffie y Beckman, 1980, pp  90.
// Modelo isotr贸pico difuso
// Liu y Jordan 1963
  
  // hr hora para la cual se har谩n los calculos, 0-24 
  // I_T Radiacion total sobre superficie inclinada, J/m^2
  // Ib irradiacion directa horaria sobre superficie horizontal, J/m^2
  // Id irradiacion directa horaria sobre superficie horizontal, J/m^2
  // I irradiaci贸n total  horaria sobre superficie horizontal, J/m^2
  // I0 irradiaci贸n extraterrrestre total  horaria sobre superficie horizontal, J/m^2
  // beta, inclinacion de la superficie
  // rhog, albedo
  // indice de claridad::: de donde sacamos este valor???
  // se necesita I (dato medido) e I0
  
  // Comprobado con el ejemplo 2.15.1, pp. 90 
  // DB 1.18 MJ/m^2 aqui 1.17 MJ/m^2
  // Use la correlacion de Orgill y Hollands
  // La de Erbs parecetener un problema

  
float Ib,Id,I0,IT,tt;
float fracd,kt,Rbb,Rba; 
float omega, omega1, omega2;//,gamma;
float cosa, cosacim,theta,thetaz;
//gamma=0.;  // sur

omega=anguloh_h(hr);
tt=tiempo_hr(omega);

omega1=anguloh_h(hr-0.5);
omega2=anguloh_h(hr+0.5);
//printf("omega1  %f ,omega2  %f  omega %f  \n",omega1,omega2, omega);

// I0 irradiaci贸n extrat. total horaria sobre sup. horizontal, J/m^2
I0=I_0(n, phi, delta, omega1, omega2);

// Indice de claridad horario, DB. 2.9.3 pp. 72
kt=I/I0;
//printf("I_Tdis; kt=I/I_0  %f I0 %f \n",kt, I0*1.e-6);
//printf("Orgill y Hollands I_d/I %f \n",frac_d_o(kt));
//printf("Erbs %f %f \n",kt, frac_d(kt));

// frac_d(kt) fraccion Id/I, correlacion de  Erbs
// DB. 2.10.1 pp. 76
Id=I*frac_d(kt);
Ib=I -Id;

//printf(" ITdis :Id  %f Ib %f \n", Id*1.e-6, Ib*1.e-6 ); //ok

cosacim=costhetaz(phi, delta, omega);
thetaz=acos(cosacim);
thetaz=radagrados(thetaz);

//printf("%f   thetaz %f \n",cosacim,  thetaz); // ok

cosa= costheta(phi,delta, omega, gamma, beta);
theta=acos(cosa);
theta=radagrados(theta);

//printf("cos theta %f  theta %f \n", cosa , theta);

// razon IT_*/I 2.14.4, pp. 87
Rbb=Rb(theta,thetaz);
Rba=Rb_a(phi, delta, beta, gamma, omega1, omega2);

beta=gradosarad(beta);

//printf("Rb %f   Rb=a/b %f\n",Rbb, Rba);

// Componentes directa, isotropica difusa y albedo
*Ibb=Ib*Rbb;
*Idd=Id*(1.+cos(beta))*0.5;
*Igg=I*rhog*(1.-cos(beta))*0.5;

// I_T Radiacion total sobre superficie inclinada, J/m^2
IT=*Ibb +*Idd +*Igg;

//fprintf(stderr,"%f %f %f %f %f %f %f %f %f \n",tt,I0*1.e-6, Ib*1.e-6, Id*1.e-6, 
//	*Ibb*1.e-6, *Idd*1.e-6,
//	*Igg*1.e-6, IT*1.e-6, Rbb); 

return *Ibb +*Idd +*Igg;
}

float I_horizontal(int n, float phi, float delta, float I, 
	      float gamma,  float hr, float albedo,
	    float * Ib, float * Id, float * Ig)
{
  // I_horizontal  discriminado en directa y difusa.
  // Luis Cardon  6 de setiembre de 2020

  // archivo: usEstimacionRad.c
  
// Radiacion sobre superficie horizontal: cielo isotropico
// en el periodo de una hora, +- 30 minutos de hr
  
// Duffie y Beckman, 1980, pp  90.
// Modelo isotropico difuso
// Liu y Jordan 1963
  
  // hr hora para la cual se haran los calculos, 0-24 
  // Ib irradiacion directa horaria sobre superficie horizontal, J/m^2
  // Id irradiacion directa horaria sobre superficie horizontal, J/m^2
  // I irradiacion total  horaria sobre superficie horizontal, J/m^2

  // I0 irradiacion extraterrrestre total  horaria sobre superficie horizontal, J/m^2
  
  // indice de claridad::: de donde sacamos este valor???
  // se necesita I (dato medido) e I0
  
  // Comprobado con el ejemplo 2.15.1, pp. 90 
  // DB 1.18 MJ/m^2 aqui 1.17 MJ/m^2
  // Use la correlacion de Orgill y Hollands
  // La de Erbs parece tener un problema

  
float I0,IT,tt;
float fracd,kt,Rbb,Rba; 
float omega, omega1, omega2;//,gamma;
float cosa, cosacim,theta,thetaz;
//gamma=0.;  // sur

omega=anguloh_h(hr);
tt=tiempo_hr(omega);

omega1=anguloh_h(hr-0.5);
omega2=anguloh_h(hr+0.5);
//printf("omega1  %f ,omega2  %f  omega %f  \n",omega1,omega2, omega);

// I0 irradiacion extrat. total horaria sobre sup. horizontal, J/m^2
I0=I_0(n, phi, delta, omega1, omega2);

// Indice de claridad horario, DB. 2.9.3 pp. 72
kt=I/I0;
//printf("I_Tdis; kt=I/I_0  %f I0 %f \n",kt, I0*1.e-6);
//printf("Orgill y Hollands I_d/I %f \n",frac_d_o(kt));
//printf("Erbs %f %f \n",kt, frac_d(kt));

// frac_d(kt) fraccion Id/I, correlacion de  Erbs
// DB. 2.10.1 pp. 76
*Id=I*frac_d(kt);
*Ib=I -*Id;
*Ig=I*albedo;
return *Id + *Ib + *Ig;
}


// ------------------------------------------------------------------
// Calculos de la irradiacion diaria
// ------------------------------------------------------------------

float H_0(int n, float phi, float delta, float omegas)
{
  // Luis Cardon - 19 setiembre 2013
  
  // radiacion solar extraterrestre diaria
  // sobre superficie horizontal
  // Duffie  y Beckman, pp. 37, ec. 1.10.3

// Entra n, dia del año
// Entra phi, latitud, gardos
// Entra delta, declinacion, grados
// Entra omegas, angulo de la salida del sol, grados
  
// (360.*n/365.) es un angulo en grados ??
// ver mismo problema en calculo de Gon  

// resultados ok contra ejemplo de DyB pp 40
// DB 33.8 calculado 33.26

// H0 radiacion solar extraterreetres sobre sup. horizontal, J/M^2

// hay una nota en mi libro REVISAR sin o cos ??
  
  float H0, angulo;
phi=gradosarad(phi);
delta=gradosarad(delta);
omegas=gradosarad(omegas);
angulo=gradosarad(360.*n/365.);

H0=(24.* 3600.*GSC/PI)* ( 1.+ 0.033* cos(angulo) ) *
( cos(phi)* cos(delta)* sin(omegas) + omegas*sin(phi)*sin(delta));
return H0;
}
// ===========================================================
// ------------------------------------------------------------------
float r_t(float omega, float omegas)
{
  // Luis Cardon - 19 setiembre 2013
  // razon de radiacion horaria total a diaria total
  // r_t Collares-Pereira and Rabl correlation \cite{collares79}
  // Duffie y Beckman ec. 2.13.2 a, pp 82
  // DB 0.076  -> 0.076 verificado
  
  float a, b, r;
  float omegas60;
  omegas60=omegas-60.;
  omegas60=gradosarad(omegas60);
  omega=gradosarad(omega);
  omegas=gradosarad(omegas);
  
a= 0.409  + 0.5016*sin(omegas60);
b= 0.6609 - 0.4767*sin(omegas60);

r  = PI*(a + b*cos(omega) )/24.;
r=r*(cos(omega) - cos(omegas))/(sin(omegas)-omegas*cos(omegas) ); 
return r;
}
// ------------------------------------------------------------------
float oor_t(sgeohr geohr)
{
  // Luis Cardon - 29 de mayo de 2014
  // modifico para que calcula para la hora señalada
  
  // Luis Cardon - 19 setiembre 2013
  // razon de radiacion horaria total a diaria total
  // r_t Collares-Pereira and Rabl correlation \cite{collares79}
  // Duffie y Beckman ec. 2.13.2 a, pp 82
  // DB 0.076  -> 0.076 verificado
  
  float a, b, r;
  float omega,omegas,omegas60;
  float hr;
  hr=geohr.hr;
  omega=geohr.omega;
  omegas=geohr.omegas;
  
  omegas60=omegas-60.;
  omegas60=gradosarad(omegas60);
  omega=gradosarad(omega);
  omegas=gradosarad(omegas);
  
a= 0.409  + 0.5016*sin(omegas60);
b= 0.6609 - 0.4767*sin(omegas60);

r  = PI*(a + b*cos(omega) )/24.;
r=r*(cos(omega) - cos(omegas))/(sin(omegas)-omegas*cos(omegas) ); 
return r;
}
// ------------------------------------------------------------------

float rtt(float omega, float omegas)
{
// omega y omegas en grados
// las funciones trigonometricas de C requieren radianes  
// Duffie y Beckman, 1980, pp  79.
// Collares-Pereira and Rabl (1979a)

// Retorna rt, razon radiacion horaria total /radiacion diaria total, I/H
// Entra omega, angulo horario en grados para punto medio de la hora de calculo
// Entra omega_s, angulo horario del amanecer
// verificada con el problema 2.13.1 pp 83 DB
  
// 19 setiembre 2013

  // signo de omega_s
float omegas60;
float coso, cosos, cosom,cosoms;
float a,b;
omegas60=omegas-60;

omega=gradosarad(omega);
omegas=gradosarad(omegas);
omegas60=gradosarad(omegas60);

coso=cos(omega);
cosos=cos(omegas);

a= 0.409  + 0.5016 *sin(omegas60);
b= 0.6609 - 0.4767 *sin(omegas60);

return (PI/24.)*(a + b *coso) *(coso-cosos)/( sin(omegas)- (omegas*cosos));
}
// ------------------------------------------------------------------
float tao_b(float thetaz, float A, int clima)
{
// Luis Cardon 20 de setiembre de 2013
  
// Duffie y Beckman, 1980, pp  62.
// Duffie y Beckman,     , pp  68
// Hottel 1976

// Retorna taob, transmitancia atmosferica a la radiacion directa 
// Entra thetaz, angulo cenital, grados
// Entra A, altitud,   kilmetros

// clima=0 tropical
// clima=1 verano de media latitud
// clima=2 verano subartico
// clima=3 invierno de media latitud

// verificado con el probelma 2.8.1 DB pp. 69
// solo con clima 1
  
float cost, a0,a1,k,a0a,a1a,ka;

float r0[]={0.95, 0.97, 0.99, 1.03};
float r1[]={0.98, 0.99, 0.99, 1.01};
float rk[]={1.02, 1.02, 1.01, 1.00};

thetaz=gradosarad(thetaz);
cost=cos(thetaz);
a0a=0.4237 -0.00821 *(6.-A)*(6.-A);
a1a=0.5055 + 0.00595* (6.5-A)*(6.5-A);
ka=0.2711 + 0.01858*(2.5-A)*(2.5-A);

a0=a0a*r0[clima];
a1=a1a*r1[clima];
k=ka*rk[clima];
return a0 + a1 * exp(-k/cost);
}
// ------------------------------------------------------------------
float ootao_b(sgeo geo,sgeohr geohr)
{
// Luis Cardon 20 de setiembre de 2013
  
// Duffie y Beckman, 1980, pp  62.
// Duffie y Beckman,     , pp  68
// Hottel 1976

// Retorna taob, transmitancia atmosferica a la radiacion directa 
// Entra thetaz, angulo cenital, grados
// Entra A, altitud,   kilmetros

// clima=0 tropical
// clima=1 verano de media latitud
// clima=2 verano subartico
// clima=3 invierno de media latitud

// verificado con el probelma 2.8.1 DB pp. 69
// solo con clima 1
float A, thetaz;
int clima;
float hr;
float cost, a0,a1,k,a0a,a1a,ka;

float r0[4]={0.95, 0.97, 0.99, 1.03};  //23km
//float r0[]={0.92, 0.96, 0.98, 1.04}; 5km
float r1[4]={0.98, 0.99, 0.99, 1.01};
float rk[4]={1.02, 1.02, 1.01, 1.00};

clima=geo.clima;
A=geo.altitud;  // km

cost=oocosthetaz(geo,geohr);

//fimprimegeo(geo);
//fimprimegeohr(geohr);
 
//fprintf(stderr,"cost=%f \n",cost);
a0a=0.4237 - 0.00/821 * (6.0-A)*(6.0-A);
a1a=0.5055 + 0.00595 * (6.5-A)*(6.5-A);
ka= 0.2711 + 0.01858 * (2.5-A)*(2.5-A);

a0=a0a*r0[clima];
a1=a1a*r1[clima];
k=ka*rk[clima];

//fprintf(stderr,"cost=%f  a0a=%f a0=%f  a1a=%f a1=%f  k=%f  -k/cots=%f \n",cost,a0a,a0,a1a,a1,k,-k/cost);
return a0 + a1 * exp(-k/cost);
}
// ------------------------------------------------------------------


float tao_d(float taob)
{
// Luis Cardon 20 de setiembre de 2013
  
// Duffie y Beckman, 1980, pp  
// Duffie y Beckman,     , pp  70
// Liu Jordan 1969
  
// tao_d=I_d/I_0  radiacion difusa /radiacion extraterrestre
//   directa sobre sup horizaontal  

// Verificado con ejemplo 2.8.2 DB pp 70  

// Retorna tao_d, transmitancia atmosferica a la radiacion directa 
  return 0.271 - 0.294*taob;
}
float Icb(float Ion, float taob, float thetaz)
{
// Luis Cardon 20 de setiembre de 2013

// Duffie y Beckman, 1980, pp  63.
thetaz=gradosarad(thetaz);  
return Ion* taob * cos(thetaz);
}
float I_T(int n, float phi, float delta, float I, float beta, float rhog, float hr)
{
  // Luis Cardon  23 de setiembre de 2013

// Radiacion sobre superficie inclinada: cielo isotropico
// en el período de una hora, +- 30 minutos de hr
  
// Duffie y Beckman, 1980, pp  90.
// Modelo isotropico difuso
// Liu y Jordan 1963
// hr hora para la cual se haran los calculos 
  // Ib irradiacion directa horaria sobre superficie horizontal, J/m^2
  // Id irradiacion directa horaria sobre superficie horizontal, J/m^2
  // I irradiacion total  horaria sobre superficie horizontal, J/m^2
  // I0 irradiacion extraterrrestre total  horaria sobre superficie horizontal, J/m^2
  // beta, inclinacion de la superficie
  // rhog, albedo

  // indice de claridad::: de donde sacamos este valor???

  // se necesita I (dato medido) e I0
  // Comprobado con el ejemplo 2.15.1 
  // DB 1.18 MJ/m^2 aqui 1.17 MJ/m^2
  // Use la correlacion de Orgill y Hollands
  // La de Erbs parecetener un problema

  // I irradiacion total  horaria sobre superficie horizontal, J/m^2
  // dato meteorologico del lugar
  
  float Ib,Id,I0;
float fracd,kt,Rbb,Rba; 
float omega, omega1, omega2,gamma;
omega=anguloh_h(hr);
omega1=anguloh_h(hr-0.5);
omega2=anguloh_h(hr+0.5);
//printf("omega1  %f ,omega2  %f  omega %f  \n",omega1,omega2, omega);
I0=I_0(n, phi, delta, omega1, omega2);
//printf("I0 %f \n",I0);
kt=I/I0;
//printf("kt i= %f %f \n",kt, I0*1.e-6);
//printf("Orgill y Hollands %f \n",frac_d_o(kt));
//printf("Erbs %f \n",frac_d(kt));
Id=I*frac_d(kt);
Ib=I*(1.-frac_d(kt));
//printf("Id Ib %f %f \n", Id*1.e-6, Ib*1.e-6 );
//printf("-------------------------------\n");
gamma=0;
// ---------------
float cosa, cosacim,theta,thetaz;
cosacim=costhetaz(phi, delta, omega);
thetaz=acos(cosacim);
thetaz=radagrados(thetaz);
//printf("%f   thetaz %f \n",cosacim,  thetaz); // ok

cosa= costheta(phi,delta, omega, gamma, beta);
theta=acos(cosa);
theta=radagrados(theta);
//printf("cos theta %f  theta %f \n", cosa , theta);




Rbb=Rb(theta,thetaz);
Rba=Rb_a(phi, delta, beta, gamma, omega1, omega2);
printf("Rb %f   Rb=a/b %f\n",Rbb, Rba);
beta=gradosarad(beta);
return Ib*Rbb+ Id*(1.+cos(beta))*0.5 + I*rhog*(1.-cos(beta))*0.5;
}
// // ======================================================
//float I_Tdis(int n, float phi, float delta, float I, 
// 	     float beta, float rhog, float hr,
// 	    float * Ibb, float * Idd, float * Igg )
// {
//   // I_T discriminado
//   // I_T Radiacion total sobre superficie inclinada
//   // Luis Cardon  23 de setiembre de 2013
// 
//   // archivo: usEstimacionRad.c
//   
// // Radiacion sobre superficie inclinada: cielo isotropico
// // en el período de una hora, +- 30 minutos de hr
//   
// // Duffie y Beckman, 1980, pp  90.
// // Modelo isotropico difuso
// // Liu y Jordan 1963
//   
//   // hr hora para la cual se haran los calculos, 0-24 
//   // I_T Radiacon total sobre superficie inclinada, J/m^2
//   // Ib irradiacion directa horaria sobre superficie horizontal, J/m^2
//   // Id irradiacion directa horaria sobre superficie horizontal, J/m^2
//   // I irradiacion total  horaria sobre superficie horizontal, J/m^2
//   // I0 irradiacion extraterrrestre total  horaria sobre superficie horizontal, J/m^2
//   // beta, inclinacion de la superficie
//   // rhog, albedo
//   // indice de claridad::: de donde sacamos este valor???
//   // se necesita I (dato medido) e I0
//   
//   // Comprobado con el ejemplo 2.15.1, pp. 90 
//   // DB 1.18 MJ/m^2 aqui 1.17 MJ/m^2
//   // Use la correlacion de Orgill y Hollands
//   // La de Erbs parecetener un problema
// 
//   
// float Ib,Id,I0,IT,tt;
// float fracd,kt,Rbb,Rba; 
// float omega, omega1, omega2,gamma;
// float cosa, cosacim,theta,thetaz;
// gamma=0.;  // sur
// 
// omega=anguloh_h(hr);
// tt=tiempo_hr(omega);
// 
// omega1=anguloh_h(hr-0.5);
// omega2=anguloh_h(hr+0.5);
// //printf("omega1  %f ,omega2  %f  omega %f  \n",omega1,omega2, omega);
// 
// // I0 irradiacion extrat. total horaria sobre sup. horizontal, J/m^2
// I0=I_0(n, phi, delta, omega1, omega2);
// 
// // Indice de claridad horario, DB. 2.9.3 pp. 72
// kt=I/I0;
// //printf("I_Tdis; kt=I/I_0  %f I0 %f \n",kt, I0*1.e-6);
// //printf("Orgill y Hollands I_d/I %f \n",frac_d_o(kt));
// //printf("Erbs %f %f \n",kt, frac_d(kt));
// 
// // frac_d(kt) fraccion Id/I, correlacion de  Erbs
// // DB. 2.10.1 pp. 76
// Id=I*frac_d(kt);
// Ib=I -Id;
// 
// //printf(" ITdis :Id  %f Ib %f \n", Id*1.e-6, Ib*1.e-6 ); //ok
// 
// cosacim=costhetaz(phi, delta, omega);
// thetaz=acos(cosacim);
// thetaz=radagrados(thetaz);
// 
// //printf("%f   thetaz %f \n",cosacim,  thetaz); // ok
// 
// cosa= costheta(phi,delta, omega, gamma, beta);
// theta=acos(cosa);
// theta=radagrados(theta);
// 
// //printf("cos theta %f  theta %f \n", cosa , theta);
// 
// // razon IT_*/I 2.14.4, pp. 87
// Rbb=Rb(theta,thetaz);
// Rba=Rb_a(phi, delta, beta, gamma, omega1, omega2);
// 
// beta=gradosarad(beta);
// 
// //printf("Rb %f   Rb=a/b %f\n",Rbb, Rba);
// 
// // Componentes directa, isotropica difusa y albedo
// *Ibb=Ib*Rbb;
// *Idd=Id*(1.+cos(beta))*0.5;
// *Igg=I*rhog*(1.-cos(beta))*0.5;
// 
// // I_T Radiacion total sobre superficie inclinada, J/m^2
// IT=*Ibb +*Idd +*Igg;
// 
// //fprintf(stderr,"%f %f %f %f %f %f %f %f %f \n",tt,I0*1.e-6, Ib*1.e-6, Id*1.e-6, 
// //	*Ibb*1.e-6, *Idd*1.e-6,
// //	*Igg*1.e-6, IT*1.e-6, Rbb); 
// 
// return *Ibb +*Idd +*Igg;
// }

// ------------------------------------------------------------------
float frac_d(float kt)
{
  // Luis Cardon  203 de setiembre de 2013
  // Duffie y Beckman pp 76
  // Correlacion de Erbs 
  // 0.73 Orgill y Hollands frente a 0.76 de Erbs
  // Fraccion de la radiacion difusa a la total horaria
  // sobre superficie horizontal
  // Verificado con ejemplo 2.15.1, Duffie y Beckman,  pp  90.
  float fra;
  float kt2,kt3,kt4;
  kt2=kt*kt;
  kt3=kt2*kt;
  kt4=kt3*kt;
  if (kt >= 0. && kt <= 0.22)
    fra= 1.0 -0.09*kt;
    else if( kt > 0.22 && kt <= 0.8)
      fra=0.9511 - 0.1604*kt + 4.388*kt2- 16.638*kt3 +12.336*kt4;
      else if(kt > 0.8)
	fra=0.165;
	else
	  fra=0;
//  printf("Erbs: %f  %f\n",kt ,fra);
  return fra;
  
}
float frac_d_o(float kt)
{
  // Luis Cardon  23 de setiembre de 2013
  // Duffie y Beckman pp 76
  // Correlacion de Orgill y Hollands
  
  // Fraccion de la radiacion difusa a la total horaria
  // sobre superficie horizontal
  // Verificado con ejemplo 2.15.1, Duffie y Beckman,  pp  90.
  
  float fra;
  if (kt >=0 && kt <= 0.35)
    fra= 1.0 -0.249*kt;
    else if( kt > 0.35 && kt <= 0.75)
      fra=1.557 - 1.84*kt;
      else if(kt > 0.75)
	fra=0.177;
	
 
  return fra;
  
}

float costheta(float phi, float delta, float omega, 
	       float gamma, float beta)
{
  // Luis Cardon  20 de setiembre de 2013
  
// Duffie y Beckman, 1980, pp  11.
// Duffie y Beckman, 1980, pp  14.

// angulos ingresados en grados
  
// Retorna cos(theta) theta angulo de incidencia de la radiacion directa, 
// Entra phi,   latitud, grados
// Entra delta, declinacion, grados
// Entra omega, angulo horario, grados
// Entra gamma, angulo acimutal de la superficie, grados
// Entra beta,  pendiente, grados
// Verificada
phi=gradosarad(phi);
delta=gradosarad(delta);
omega=gradosarad(omega);
gamma=gradosarad(gamma);
beta=gradosarad(beta);
return sin(delta) *sin(phi) *cos(beta) - 
sin(delta) *cos(phi) *sin(beta)*cos(gamma) +
cos(delta)*cos(phi)*cos(beta)*cos(omega)+
cos(delta)*sin(phi)*sin(beta)*cos(gamma)*cos(omega)+
cos(delta)*sin(beta)*sin(gamma)*sin(omega);
}
// ------------------------------------------------------------------

float costhetaz(float phi, float delta, float omega)
{
  // Luis Cardon - 20 setiembre de 2013
  
// Duffie y Beckman, 1980, pp  11.
// Duffie y Beckman, 1980, pp  15.
// Comprobado con el ejemplo de DyB, 1.6.2 pp. 16 
// Retorna cos(thetaz)
// thetaz:  angulo de incidencia de la radiacion directa,
// para superficie horizontal: angulo cenital   
// phi, latitud
// delta, declinacion
// omega, angulo horario
// beta=0


phi=gradosarad(phi);
delta=gradosarad(delta);
omega=gradosarad(omega);


return cos(delta) * cos(phi)*cos(omega) +sin(delta)*sin(phi);
}

float oocosthetaz(sgeo geo,sgeohr geohr)
{
  // Luis Cardon - 20 setiembre de 2013
  
// Duffie y Beckman, 1980, pp  11.
// Duffie y Beckman, 1980, pp  15.
// Comprobado con el ejemplo de DyB, 1.6.2 pp. 16 
// Retorna cos(thetaz)
// thetaz:  angulo de incidencia de la radiacion directa,
// para superficie horizontal: angulo cenital   
// phi, latitud
// delta, declinacion
// omega, angulo horario
// beta=0

// ojo: puede trabajar con geohr.hr o geohr.omega !!
// ahora con geohr.omega  
float phi,delta,omega;
phi=geo.phi;
delta=geo.delta;
//omega=anguloh_h(geohr.hr);
omega=geohr.omega;
phi=gradosarad(phi);
delta=gradosarad(delta);
omega=gradosarad(omega);


return cos(delta) * cos(phi)*cos(omega) + sin(delta)*sin(phi);
}

// ------------------------------------------------------------------
float gradosarad(float grados)
{ 
  // convierte  radiantes a grados
  double Pi=3.141592653589793238462643383279502884;
  return  grados * 3.1415/ 180.;
}
// ------------------------------------------------------------------
float radagrados(float rad)
{ // convierte grados a radiantes
  double Pi=3.141592653589793238462643383279502884;
  return  rad*180. / Pi;
}
// ------------------------------------------------------------------
float delta_n(int n)
{
// Luis Cradon -20 setiembre 2013
// Cooper, 1969
// Duffie y Beckman, pp. 14
// verificado datos DyB pp. 115, error 0.1  
  
// n, dia del año
// delta, declinacion,   grados
  
float angulo;
angulo=360.*(284.+n)/365.;
angulo=gradosarad(angulo);
return  23.45* sin(angulo);
}
// ------------------------------------------------------------------
float Rb(float theta, float thetaz)
{
  // Luis Cardon 20 setiembre de 2013
  
//   Retorna Rb, razon de la radiacion directa sobre superfice 
//   inclinada a horizontal

float costheta, costhetaz;

theta=gradosarad(theta);
thetaz=gradosarad(thetaz);
//printf("Rb: %f %f \n",cos(thetarad), cos(thetazrad) );
return cos(theta)/cos(thetaz);
}
float Rb_a(float phi, float delta, float beta,
	   float gamma, float omega1, float omega2)
{
  // Luis Cardon 20 setiembre de 2013
  
// Retorna Rb,av, razon radiacion directa sobre superfice 
//   inclinada / horizaontal
// Duffie y Beckman, pp 89
float a,b;  
phi=gradosarad(phi);
delta=gradosarad(delta);
omega1=gradosarad(omega1);
omega2=gradosarad(omega2);
beta=gradosarad(beta);
gamma=gradosarad(gamma);

a= ( sin(delta)* sin(phi)*sin(beta) -
sin(delta)*cos(phi)*sin(beta)*cos(gamma) )*(omega2-omega1);
a= a + (cos(delta)*cos(phi)*cos(beta) +
cos(delta)*sin(phi)*sin(beta)*cos(gamma) )*
(sin(omega2)-sin(omega1));
a= a -(cos(delta)*sin(beta)*sin(gamma))*(cos(omega2)-cos(omega1));

b= (cos(phi)*cos(delta))*(sin(omega2)-sin(omega1)) 
+ (sin(phi)*sin(delta))*(omega2-omega1);
//printf("a %f b %f\n",a,b);
return a/b;
}
float ooRbab(sgeo geo, sgeohr geohr)
{
  // Luis Cardon 20 setiembre de 2013
  
// Retorna Rb,av, razon radiacion directa sobre superfice 
//   inclinada / horizaontal
// Duffie y Beckman, pp 89
float a,b;  
float phi,delta,hr,omega1,omega2,beta,gamma;
phi=gradosarad(geo.phi);
delta=gradosarad(geo.delta);
hr=geohr.hr;

omega1=anguloh_h(hr-0.5);
omega2=anguloh_h(hr+0.5);

omega1=gradosarad(omega1);
omega2=gradosarad(omega2);
beta=gradosarad(geo.beta);
gamma=gradosarad(geo.gamma);

a= ( sin(delta)* sin(phi)*sin(beta) -
sin(delta)*cos(phi)*sin(beta)*cos(gamma) )*(omega2-omega1);
a= a + (cos(delta)*cos(phi)*cos(beta) +
cos(delta)*sin(phi)*sin(beta)*cos(gamma) )*
(sin(omega2)-sin(omega1));
a= a -(cos(delta)*sin(beta)*sin(gamma))*(cos(omega2)-cos(omega1));

b= (cos(phi)*cos(delta))*(sin(omega2)-sin(omega1)) 
+ (sin(phi)*sin(delta))*(omega2-omega1);
//printf("a %f b %f\n",a,b);
return a/b;
}


float mH_dmHratio(float omegas, float KT)
{
// returns Manthly averaged diffuse fraction 
// $\overline{H_d}/\overline{H}$ 
// according to Collares-Pereira and Rabl 
// Dufie beckman 
float ratio;
ratio= 0.775 + 0.00653 *(omegas -90.)  - (0.505 + 0.00455*(omegas -90.))*cos(115.*KT -103. );
return ratio;
}
// ------------------------------------------------------------------
float rd(float omegas, float omega)
{
// Returns the total dayly to  to diffuse hourly ratio $r_d$
// from Liu and Jordan correlation,  \cite{liu60}:
float r_d;
r_d  = PI/24. *( cos(omega) - cos(omegas)/sin(omegas) - (2.*PI* omegas/360.)*cos(omegas)); 
return r_d;
}
float oord(sgeohr geohr)
{
// Returns the total dayly to  to diffuse hourly ratio $r_d$
// from Liu and Jordan correlation,  \cite{liu60}:
float r_d;
  float hr,omega,omegas;
  hr=geohr.hr;
  omega=geohr.omega;
  omegas=geohr.omegas;

r_d  = PI/24. *( cos(omega) - cos(omegas)/sin(omegas) - (2.*PI* omegas/360.)*cos(omegas)); 
return r_d;
}

// ------------------------------------------------------------------
float omega_s(float phi, float delta)
{
  // Angulo de la salida del sol en grados
  // Dufie y Beckman, pp 17, ec. 1.6.10
  // phi, latitud, grados
  // delta, declinacion, grados  
  // ok, comprobado con problema DB pp 18.
  // Luis Cardon 19 setiembre 2013
  
  float omegas, cosomegas;
phi=gradosarad(phi);
delta=gradosarad(delta);
//cosomegas= -sin(phi)* sin(delta)/cos(phi)*cos(delta);
cosomegas= -tan(phi)*tan(delta);
omegas=acos(cosomegas);
omegas=radagrados(omegas);
return omegas;
}

float omega_min(float min)
{
  // devuelve el angulo horario 

  // minutos antes o despues de mediodia
  return 0.25*min;
  
}
float omega_seg(float seg)
{
  // seg: segundos antes o despues de mediodia
  // negativo antes de mediodia
  // devuelve el angulo horario 
  // segundos antes o despues de mediodia
  return (15./3600.)*seg;
  
}
float omega_hr(float hr)
{
  // Luis Cardon - 19 setiembre 2013
  // devuelve el angulo horario 
  
  // horas antes o despues de mediodia
  // hr antes de mediodia debe ser negativo
  return 15.*hr;
}
float anguloh_h(float hora)
{
  float h, anguloh;
  // Luis Cardon - 19 setiembre 2013
  // entra la hora: 0-24
  // devuelve el angulo horario en grados
  
  // h  numero de horas antes o despues de mediodia
  // 9:30 -> -40.499 frente a 37.5 de DB pp. 16 (( facil de comprobar)
  // ojo 9:30 es 9.5 ok
  
  h= - (12. - hora);  // negativo  antes de mediodia
  anguloh= omega_hr(h);
  return anguloh;
}
float anguloh_min(float min)
{
  float h, anguloh;
  // entra la hora: 0-24
  // devuelve el angulo horario
  
  // h  numero de horas antes o despues de mediodia
  
  h= - (12.*60 - min);  // negativo  antes de mediodia
  anguloh= omega_min(h);
  return anguloh;
}
float anguloh_seg(float seg)
{
  float h, anguloh;
  // entra la hora: 0-24
  // devuelve el angulo horario
  
  // h  numero de horas antes o despues de mediodia
  
  h= - (12.*60*60 - seg);  // negativo  antes de mediodia
  anguloh= omega_seg(h);
  return anguloh;
}
float tiempo_seg(float angulo)
{
  // Luis Cardon -19 setiembre 2013
  // convierte el angulo horario a hora 0-24*60*60 seg
  
  // angulo:  angulo horario, grados
  // seg:     segundos antes o despues de mediodia 

  float seg;

  seg=angulo*3600./15.;
  
  return  12.*60.*60. + seg ;
}
float tiempo_min(float angulo)
{ float min;
  // probar 
  // angulo:  angulo horario, grados
  // devuelve segundos desde las cero horas
  // convierte el angulo horario a hora 0-24*60*60
  
  // min:     minutos antes o despues de mediodia 
  min=angulo*60./15.;
  
  return  12.*60. + min ;
}
float tiempo_hr(float angulo)
{
  // Luis Cardon -19 setiembre 2013
  // convierte el angulo horario a hora 0-24 hr
  
  // angulo:  angulo horario, grados
  // hr:     horas antes o despues de mediodia 

  float hr;

  hr=angulo/15.;
  
  return  12. + hr ;
}
/*float costheta2(float phi, float delta, float omega, 
	       float gamma, float beta)
{
// Duffie y Beckman, 1980, pp  11.
// Duffie y Beckman, , pp  15.

// angulos ingresados en grados
  
// Retorna cos(theta) theta angulo de incidencia de la radiacion directa, 
// Entra phi,   latitud, grados
// Entra delta, declinacion, grados
// Entra omega, angulo horario, grados
// Entra gamma, angulo acimutal de la superficie, grados
// Entra beta,  pendiente, grados

float phimbeta;
phimbeta=phi-beta;
phimbeta=gradosarad(phimbeta);

phi=gradosarad(phi);
delta=gradosarad(delta);
omega=gradosarad(omega);
gamma=gradosarad(gamma);
beta=gradosarad(beta);

return cos(phimbeta)*cos(delta)*cos(omega) +
sin(phimbeta)*sin(delta);
}

 /*
float O_theta_hr(struct s_geo geo, struct s_geohr geohr)
{
// Mal hecha por ahora
float hr,phi,delta,gamma,beta;
float omega,costeta,theta;

phi=geo.phi;
delta=geo.delta;
beta=geo.beta;
gamma=geo.gamma;

hr=geohr.hr;

omega=anguloh_h(hr);
costeta= costheta(phi,delta, omega, gamma, beta);
theta=acos(costeta);
theta=radagrados(theta);
return theta;
}
*/
srad  oohottel(sgeo geo, sgeohr geohr)
 {  // Requiere revision mayor.
 
   // Luis Cardon 5 junio de 2014
   // Metodo de Hottel-Liu Jordan
   // Calcula la radiacion horaria a 
   // partir de la radiacion horaria extraterrestre
   // y la transmitancia atmosferica de dia claro
   // DB. pp 68.
   // I0 radiacion extraterrestre sobre superficie horizonra, horaria, O(1)MJ/m²
   
   float I0,rt,rd,rho;
   float Icb,Icd,Ic,Ig;
   srad I;
   rho=0.6;
   
I0=ooI0(geo, geohr);
rt=oor_t(geohr);
rd=oord(geohr);
// 
Icb=I0*rt;
Icd=I0*rd;
Ic=Icb+Icd;
Ig=Ic*rho;

I.I=Ic;
I.Ib=Icb;
I.Id=Icd;
I.Ig=Ig;

fprintf(stderr, "oo hootel:\t%f\t %5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f \n", 
	                     geohr.hr, I0*1.e-6, rt, I.Ib*1.e-6,
	                     rd, I.Id*1.e-6,
	                     I.Ig*1.e-6,
	                     I.I*1.e-6);  
return I;    
   
}
void ooprintrad(sgeohr geohr,srad rad)
{
fprintf(stderr, "%f %5.2f    %5.2f %5.2f  %5.2f \n",geohr.hr,rad.I*1.e-6,rad.Ib*1.e-6,rad.Id*1.e-6,rad.Ig*1.e-6);  
}
srad ooSzero()
{ 
  srad S;
S.I=0.;
S.Ib=0.;
S.Id=0.;
S.Ig=0.;

  return S;
}
//===============================================================
srad  oohottelDB(sgeo geo, sgeohr geohr)
 {
   // Luis Cardon 5 junio de 2014
   // Metodo de Hottel-Liu Jordan
   // Calcula la radiacion horaria sobre superficie horizontal a  
   // partir de la radiacion horaria extraterrestre
   // y la transmitancia atmosferica de dia claro
   // DB. pp 68.
   // I0 radiacion extraterrestre sobre superficie normal, horaria, O(1)MJ/m²
   
   float I0,I0n, Gon,rt,rd,rho;
   float Icb,Icd,Ic,Ig;
   float costetaz, omega,taob,taod;
   srad I;

   rho=geo.albedo;
      // angulo horario, grados

   costetaz=oocosthetaz(geo,geohr);
   Gon=ooGon(geo);
   I0n=Gon*3600.;
   I0=I0n*costetaz;
   
taob=ootao_b(geo,geohr);
taod=tao_d(taob);

// 
Icb=I0*taob;
Icd=I0*taod;
Ic=Icb+Icd;
Ig=Ic*geo.albedo;

I.I=Ic;
I.Ib=Icb;
I.Id=Icd;
I.Ig=Ig;
/*
fprintf(stderr, "ooo hootel:\t%f\t%f\t %5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f \n", 
	                     geohr.hr, costetaz, I0*1.e-6, taob, I.Ib*1.e-6,
	                     taod, I.Id*1.e-6,
	                     I.Ig*1.e-6,
	                     I.I*1.e-6);  
*/
return I;    
   
}

srad ooITiso(sgeo geo, sgeohr geohr, srad I)
{
  
  // IooITiso  radiacion sobre superficie inclinada
  // modelo isotropico de Liu y Jorada 1963
  // en el período de una hora, +- 30 minutos de hr
  // Duffie y Beckman, 1980, pp  90.
  
  // Luis Cardon  9  de junio  de 2014

  // archivo: usEstimacionRad.c
  
// datos:
//    I.I : radiacion global sobre superficie horizontal, J/m^2
//    no usa las demas componentes (difusa y albedo), que recalcula
//    geohr.hr hora para la cual se haran los calculos, 0-24 
//    geo.beta, inclinacion de la superficie
//    geo.albedo, albedo

// Calculos auxiliares
//  indice de claridad::: de donde sacamos este valor???
//  I0 irradiacion extraterrrestre total  horaria sobre superficie horizontal, J/m^2


//  Retorna  
//  IT.I Radiacon  sobre superficie inclinada, J/m^2
//  IT.Ib irradiacion directa horaria sobre superficie inclinada, J/m^2
//  IT.Id irradiacion difusa horaria sobre superficie inclinada, J/m^2
//  IT.Ig irradiacion reflejada horaria sobre superficie inclinada, J/m^2

  
  // Comprobado con el ejemplo 2.15.1, pp. 90 
  // DB 1.18 MJ/m^2 aqui 1.17 MJ/m^2
  // Use la correlacion de Orgill y Hollands
  // La de Erbs parece tener un problema

srad IT;  
float Ib,Id,I0,I0n,tt,hr,rhog,Gon;
float fracd,kt,Rbb,Rba; 
float omega, omega1, omega2,gamma;
float cosa, cosacim,theta,thetaz;
float costetaz,beta;
hr=geohr.hr;
rhog=geo.albedo;

omega=anguloh_h(hr);
tt=tiempo_hr(omega);

omega1=anguloh_h(hr-0.5);
omega2=anguloh_h(hr+0.5);

// Calculo de la radiacion extraterrestre
   costetaz=oocosthetaz(geo,geohr);
   Gon=ooGon(geo);
   I0n=Gon*3600.;
   I0=I0n*costetaz;

// Calculo del indice de claridad horario, DB. 2.9.3 pp. 72
kt=I.I/I0;
//printf("I_Tdis; kt=I/I_0  %f I0 %f \n",kt, I0*1.e-6);
//printf("Orgill y Hollands I_d/I %f \n",frac_d_o(kt));
//printf("Erbs %f %f \n",kt, frac_d(kt));

// frac_d(kt) fraccion Id/I, correlacion de  Erbs
// DB. 2.10.1 pp. 76

// Calculo de las componentes difusa y reflejada sobre
// superficie horizontal
Id=I.I*frac_d(kt);
Ib=I.I - Id;

// calculo de los factores de vista
// razon IT_*/I 2.14.4, pp. 87
Rba=ooRbab(geo,geohr);
beta=gradosarad(geo.beta);

//printf("Rb %f   Rb=a/b %f\n",Rbb, Rba);

// Calculo de las componentes directa, isotropica difusa y albedo
IT.Ib  = Ib * Rba;
IT.Id  = Id * (1.+ cos(beta))*0.5;
IT.Ig  = I.I  * rhog*(1.-cos(beta))*0.5;

// Calculo de la radiacion total sobre superficie inclinada, J/m^2
IT.I= IT.Ib  + IT.Id ;


return IT;
}

// ------------------------------------------------------------------

 int ordinal(int dd, int mm, int yy)
{
  // programmingforstarters.blogspot.com
int julian;

julian=dd;
mm=mm-1;
switch(mm)
{
case 11 : julian=julian+30;
case 10 : julian=julian+31;
case 9 : julian=julian+30;
case 8 : julian=julian+31;
case 7 : julian=julian+31;
case 6: julian=julian+30;
case 5: julian=julian+31;
case 4 : julian=julian+30;
case 3 : julian=julian+31;
case 2 : if(yy%4==0)julian=julian+29;
             else julian=julian+28;
case 1 : julian=julian+31;
}
return julian;

  
}

srad ooITkt(sgeo geo, sgeohr geohr,int mes)
{
   // Luis Cardon - 24 de abril 2017
   
  // IooITkt  radiacion sobre superficie inclinada
  // Se estima la radiacion horaria con
  // la razon difusa/directa  en funcion de kt
  // se usa la correlacion de Erbs.
  
  // la radiacion horaria se calcula  a 
  // partir de la razon diaria horaria, r_t
  // en el período de una hora, +- 30 minutos de hr
  // Duffie y Beckman, 1980, pp  90.
  
  // Luis Cardon  9  de junio  de 2014

  // archivo: usEstimacionRad.c
  
// datos:
//    I.I : radiacion global sobre superficie horizontal, J/m^2
//    no usa las demas componentes (difusa y albedo), que recalcula
//    geohr.hr hora para la cual se haran los calculos, 0-24 
//    geo.beta, inclinacion de la superficie
//    geo.albedo, albedo

// Calculos auxiliares
//  indice de claridad::: de donde sacamos este valor???
//  I0 irradiacion extraterrrestre total  horaria sobre superficie horizontal, J/m^2


//  Retorna  
//  IT.I Radiacon  sobre superficie inclinada, J/m^2
//  IT.Ib irradiacion directa horaria sobre superficie inclinada, J/m^2
//  IT.Id irradiacion difusa horaria sobre superficie inclinada, J/m^2
//  IT.Ig irradiacion reflejada horaria sobre superficie inclinada, J/m^2

  
  // Comprobado con el ejemplo 2.15.1, pp. 90 
  // DB 1.18 MJ/m^2 aqui 1.17 MJ/m^2
  // Use la correlacion de Orgill y Hollands
  // La de Erbs parece tener un problema

srad IT;  
float I,Ib,Id,I0,I0n,tt,hr,rhog,Gon;
float fracd,kt,Rbb,Rba; 
float omega, omega1, omega2,gamma;
float cosa, cosacim,theta,thetaz;
float costetaz,beta;
hr=geohr.hr;
rhog=geo.albedo;

omega=anguloh_h(hr);
tt=tiempo_hr(omega);

omega1=anguloh_h(hr-0.5);
omega2=anguloh_h(hr+0.5);

// Calculo de la radiacion extraterrestre
   costetaz=oocosthetaz(geo,geohr);
   Gon=ooGon(geo);
   I0n=Gon*3600.;
   I0=I0n*costetaz;

// Calculo del indice de claridad horario, DB. 2.9.3 pp. 72
float H;           // este es el dato meteorologico disponible
H=geo.Hmm[mes];          // MJ julio
I=H*oor_t(geohr);   // 

kt=I/I0;

// Calculo de las componentes difusa y reflejada sobre
// superficie horizontal
Id=I*frac_d(kt);
Ib=I - Id;

// calculo de los factores de vista
// razon IT_*/I 2.14.4, pp. 87
Rba=ooRbab(geo,geohr);
beta=gradosarad(geo.beta);

//printf("Rb %f   Rb=a/b %f\n",Rbb, Rba);

// Calculo de las componentes directa, isotropica difusa y albedo
IT.Ib  = Ib * Rba;
IT.Id  = Id * (1.+ cos(beta))*0.5;
IT.Ig  = I  * rhog*(1.-cos(beta))*0.5;

// Calculo de la radiacion total sobre superficie inclinada, J/m^2
IT.I= IT.Ib  + IT.Id ;


return IT;
}
srad ooIkt(sgeo geo, sgeohr geohr,int mes)
{
   // Luis Cardon - 24 de abril 2017

  // ooIkt  radiacion sobre superficie horizontal
  // recortadade 
   
  // ooITkt  radiacion sobre superficie inclinada
  // Se estima la radiacion horaria con
  // la razon difusa/directa  en funcion de kt
  // se usa la correlacion de Erbs.
  
  // la radiacion horaria se calcula  a 
  // partir de la razon diaria horaria, r_t
  // en el período de una hora, +- 30 minutos de hr
  // Duffie y Beckman, 1980, pp  90.
  
  // Luis Cardon  4  de setiembre   de 2020

  // archivo: usEstimacionRad.c
  
// datos:
//    I.I : radiacion global sobre superficie horizontal, J/m^2
//    no usa las demas componentes (difusa y albedo), que recalcula
//    geohr.hr hora para la cual se haran los calculos, 0-24 
//    geo.beta, inclinacion de la superficie
//    geo.albedo, albedo

// Calculos auxiliares
//  indice de claridad::: de donde sacamos este valor???
//  I0 irradiacion extraterrrestre total  horaria sobre superficie horizontal, J/m^2


//  Retorna  
//  IT.I Radiacon  sobre superficie inclinada, J/m^2
//  IT.Ib irradiacion directa horaria sobre superficie inclinada, J/m^2
//  IT.Id irradiacion difusa horaria sobre superficie inclinada, J/m^2
//  IT.Ig irradiacion reflejada horaria sobre superficie inclinada, J/m^2

  
  // Comprobado con el ejemplo 2.15.1, pp. 90 
  // DB 1.18 MJ/m^2 aqui 1.17 MJ/m^2
  // Use la correlacion de Orgill y Hollands
  // La de Erbs parece tener un problema

srad IT;  
float I,Ib,Id,I0,I0n,tt,hr,rhog,Gon;
float fracd,kt,Rbb,Rba; 
float omega, omega1, omega2,gamma;
float cosa, cosacim,theta,thetaz;
float costetaz,beta;
hr=geohr.hr;
rhog=geo.albedo;

omega=anguloh_h(hr);
tt=tiempo_hr(omega);

omega1=anguloh_h(hr-0.5);
omega2=anguloh_h(hr+0.5);

// Calculo de la radiacion extraterrestre
   costetaz=oocosthetaz(geo,geohr);
   Gon=ooGon(geo);
   I0n=Gon*3600.;
   I0=I0n*costetaz;

// Calculo del indice de claridad horario, DB. 2.9.3 pp. 72
float H;           // este es el dato meteorologico disponible
H=geo.Hmm[mes];          // MJ julio
I=H*oor_t(geohr);   // 

kt=I/I0;

// Calculo de las componentes difusa y reflejada sobre
// superficie horizontal
Id=I*frac_d(kt);
Ib=I - Id;

// calculo de los factores de vista
// razon IT_*/I 2.14.4, pp. 87

//Rba=ooRbab(geo,geohr); suprimido 
//beta=gradosarad(geo.beta);

//printf("Rb %f   Rb=a/b %f\n",Rbb, Rba);

// Calculo de las componentes directa, isotropica difusa y albedo
IT.Ib  = Ib; // * Rba;
IT.Id  = Id; // * (1.+ cos(beta))*0.5;
IT.Ig  = I  * rhog; //*(1.-cos(beta))*0.5;

// Calculo de la radiacion total sobre superficie inclinada, J/m^2
IT.I= IT.Ib  + IT.Id ;


return IT;
}
// ======================================================
float I_Tdis(int n, float phi, float delta, float I, 
	     float beta, float rhog, float hr,
	    float * Ibb, float * Idd, float * Igg )
{
  //Traida de /LINSE...TANQUES-4
  // I_T discriminado
  // I_T Radiacion total sobre superficie inclinada
  // Luis Cardon  23 de setiembre de 2013

  // archivo: usEstimacionRad.c
  
// Radiacion sobre superficie inclinada: cielo isotropico
// en el período de una hora, +- 30 minutos de hr
  
// Duffie y Beckman, 1980, pp  90.
// Modelo isotropico difuso
// Liu y Jordan 1963
  
  // hr hora para la cual se haran los calculos, 0-24 
  // I_T Radiacon total sobre superficie inclinada, J/m^2
  // Ib irradiacion directa horaria sobre superficie horizontal, J/m^2
  // Id irradiacion directa horaria sobre superficie horizontal, J/m^2
  // I irradiacion total  horaria sobre superficie horizontal, J/m^2
  // I0 irradiacion extraterrrestre total  horaria sobre superficie horizontal, J/m^2
  // beta, inclinacion de la superficie
  // rhog, albedo
  // indice de claridad::: de donde sacamos este valor???
  // se necesita I (dato medido) e I0
  
  // Comprobado con el ejemplo 2.15.1, pp. 90 
  // DB 1.18 MJ/m^2 aqui 1.17 MJ/m^2
  // Use la correlacion de Orgill y Hollands
  // La de Erbs parecetener un problema

  
float Ib,Id,I0,IT,tt;
float fracd,kt,Rbb,Rba; 
float omega, omega1, omega2,gamma;
float cosa, cosacim,theta,thetaz;
gamma=0.;  // sur

omega=anguloh_h(hr);
tt=tiempo_hr(omega);

omega1=anguloh_h(hr-0.5);
omega2=anguloh_h(hr+0.5);
//printf("omega1  %f ,omega2  %f  omega %f  \n",omega1,omega2, omega);

// I0 irradiacion extrat. total horaria sobre sup. horizontal, J/m^2
I0=I_0(n, phi, delta, omega1, omega2);

// Indice de claridad horario, DB. 2.9.3 pp. 72
kt=I/I0;
//printf("I_Tdis; kt=I/I_0  %f I0 %f \n",kt, I0*1.e-6);
//printf("Orgill y Hollands I_d/I %f \n",frac_d_o(kt));
//printf("Erbs %f %f \n",kt, frac_d(kt));

// frac_d(kt) fraccion Id/I, correlacion de  Erbs
// DB. 2.10.1 pp. 76
Id=I*frac_d(kt);
Ib=I -Id;

//printf(" ITdis :Id  %f Ib %f \n", Id*1.e-6, Ib*1.e-6 ); //ok

cosacim=costhetaz(phi, delta, omega);
thetaz=acos(cosacim);
thetaz=radagrados(thetaz);

//printf("%f   thetaz %f \n",cosacim,  thetaz); // ok

cosa= costheta(phi,delta, omega, gamma, beta);
theta=acos(cosa);
theta=radagrados(theta);

//printf("cos theta %f  theta %f \n", cosa , theta);

// razon IT_*/I 2.14.4, pp. 87
Rbb=Rb(theta,thetaz);
Rba=Rb_a(phi, delta, beta, gamma, omega1, omega2);

beta=gradosarad(beta);

//printf("Rb %f   Rb=a/b %f\n",Rbb, Rba);

// Componentes directa, isotropica difusa y albedo
*Ibb=Ib*Rbb;
*Idd=Id*(1.+cos(beta))*0.5;
*Igg=I*rhog*(1.-cos(beta))*0.5;

// I_T Radiacion total sobre superficie inclinada, J/m^2
IT=*Ibb +*Idd +*Igg;

//fprintf(stderr,"%f %f %f %f %f %f %f %f %f \n",tt,I0*1.e-6, Ib*1.e-6, Id*1.e-6, 
//	*Ibb*1.e-6, *Idd*1.e-6,
//	*Igg*1.e-6, IT*1.e-6, Rbb); 

return *Ibb +*Idd +*Igg;
}
