#include "param.h"
#include "uEstructuras.h"
#include "uEstimacionRad.h"
#include "uColectores.h"
#include "uCubiertas.h"
#include "uSSS.h"

float SSS(sgeo geo, sgeohr geohr, scubierta cubierta, int mes)
{
// Calculo de la radiacion absorbida
// Ejemplo 5.9.1, pp. 222, verificado en pruebaCubiertas.c 

// Datos
// Irradiacion de dia claro para la hora: I
// este dato no entro en la lista  26 de mayo 2014

int n;  
float delta,phi,gamma,beta;
float hr,omega,theta;

float costeta;   
float I,rhog, Ibbb,Iddd,Ittt ;
float alpha,ptabd,ptadd,ptagd,rhod,Ss;
float KL,n2;
int N;

// datos
// lugar y hora
phi=geo.phi;           // latitud
n=geo.n;               // febrero 20
beta=geo.beta;         // inclinacion de la superficie
delta=geo.delta;       // 11.57, ok
gamma=geo.gamma;

hr=geohr.hr;            // hora
omega=geohr.omega;
rhog=0.6;          // reflectancia del piso, dato tipico de la locacion

KL=cubierta.KL;           //0.037;          // indice KL de extincion
alpha=cubierta.alpha;     //0.93;        // absortancia
N=cubierta.N;             // numero de cubiertas
n2=cubierta.n2;           //1.526;          // indice de refraccion de mat. de la cubierta

// datos de la cubierta e inclinacion 
hr=geohr.hr;            // hora
theta=geohr.theta;


// I irradiacion total  horaria sobre superficie horizontal, J/m^2
//I=1.04*1.e6;       // dato meteorologico del lugar

float H;           // este es el dato meteorologico disponible
H=geo.Hmm[mes];          // MJ julio
I=H*oor_t(geohr);   // 

// Calculo de la radiacion sobre superficie inclinada
I_Tdis(n, phi, delta, I, beta, rhog, hr, &Ibbb, &Iddd, &Ittt);

// Calculo de los productos transmitancia absortancia, uCubiertas.c
taos( 1, KL,  theta, beta, alpha, n2, &ptabd, &ptadd, &ptagd);

// Calculo de la radiacion absorbida sobre superficie inclinada cubierta
Ss=SSdis(&Ibbb, &Iddd, &Ittt, ptabd  ,ptadd  ,ptagd);

// Fin  de Radiacion sobre superficie inclinada: S");

//printf("SSS: %f\n",Ss);
return Ss;  
}
// =================================================================
srad  ooSSS(sgeo geo, sgeohr geohr, scubierta cubierta, srad Ibeta)
{
// Calculo de la radiacion absorbida
// Ejemplo 5.9.1, pp. 222, verificado en pruebaCubiertas.c 

// Datos
// Irradiacion de dia claro para la hora: I


srad S;
int n;  
float delta,phi,gamma,beta;
float hr,omega,theta;

float costeta;   
//float Icb,Icd,Ic,taod;
float rhog ;
float alpha,ptabd,ptadd,ptagd,rhod;
float KL,n2;
int N;

// datos
// lugar y hora
phi=geo.phi;           // latitud
n=geo.n;               // febrero 20
beta=geo.beta;         // inclinacion de la superficie
delta=geo.delta;       // 11.57, ok
gamma=geo.gamma;

hr=geohr.hr;            // hora
omega=geohr.omega;
rhog=geo.albedo;          // reflectancia del piso, dato tipico de la locacion

KL=cubierta.KL;           //0.037;          // indice KL de extincion
alpha=cubierta.alpha;     //0.93;        // absortancia
N=cubierta.N;             // numero de cubiertas
n2=cubierta.n2;           //1.526;          // indice de refraccion de mat. de la cubierta


// datos de la cubierta e inclinacion 
hr=geohr.hr;            // hora
theta=geohr.theta;





// Calculo de los productos transmitancia absortancia, uCubiertas.c
taos( N, KL,  theta, beta, alpha, n2, &ptabd, &ptadd, &ptagd);

S.Ib = Ibeta.Ib*ptabd;
S.Id = Ibeta.Id*ptadd;
S.Ig = Ibeta.Ig*ptagd;
S.I  = S.Ib+S.Id;
  


return S;  
}
