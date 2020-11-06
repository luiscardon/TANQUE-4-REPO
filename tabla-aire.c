#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "param.h"
#include "uEstructuras.h"
#include "uPropiedades.h"
#include "uPropFisicas.h"
#include "uCorrelaciones.h"

#define SIGMA 5.67e-8 // cosntante de Stefan Boltzaaman, WM^-2K^-4


int main()
{
  // tabla-aire.c 
  // i7/TANQUEs-4/TANQUE-4-REPO  5 nov 2020
  // genera una tabla de propiedades del aire
  // usa las librerias  modificadas y actualizadas
  //"uPropiedades.h"
  //"uPropFisicas.h"
  // se compila con libSE
FILE * log;
log = fopen("tabla-aire.log", "w");

fprintf(log,"===========tabla-aire.c ====================\n");
fprintf(log,"LUS CARDON 2020 \n");
fprintf(log," Genera una tabla de propiedades del aire.\n");

fprintf(log,"Programa tabla-aire.c \n");
fprintf(log,"En i7/TANQUEs-4/TANQUE-4-REPO. Compilado con LIBSE. \n");
fprintf(log,"Archivo: tabla-aire.log\n");

float T,TT;
float rho,mu,nu,k,Cp,gban,alpha,Pr;
T=27;


rho=rho_aire(T); //Ok
mu=mu_aire(T); //ok
nu=nu_aire(T); // ok
k=k_aire(T);  // ok
Cp=Cp_aire(T); // ok
gban=gban_aire(T); //ok
alpha=alpha_aire(T); //ok
Pr=Pr_aire(T);

// for en K para comparar
fprintf(stderr,"T\trho \tmu\tk\tCp\trho/nu\tgban\talpha \n");
fprintf(stderr,"[K]\t[kg/m^3]\t[10^-6Ns/m^2]\t[10^-3W/mK]\t[J/kgK]\t[10^3s/m^2]\t[1061/m^3K]\t[10-6m^2/s] \n");
for (TT=200 ; TT<=  400 ; TT=TT+10.)
{
  // Las funciones reuieren temperaturas en grados
  
T=TT-273.15;
rho=rho_aire(T); //Ok
mu=mu_aire(T); //ok
nu=nu_aire(T); // ok
k=k_aire(T);  // ok
Cp=Cp_aire(T); // ok
gban=gban_aire(T); //ok
alpha=alpha_aire(T); //ok
fprintf(stderr,"%3.0f\t%1.4f\t%2.2f\t%2.2f\t%4.1f\t%3.1f\t%3.1f\t%2.2f \n",TT, rho, mu, k*1.E3, Cp, 1.*1.E-3/nu ,gban*1.E-6, alpha*1.E6);  
}

return 0;
}

// =============================================
double hradplacasp(double epsilon1,double epsilon2, double T1,double T2)
{ 
  // Coeficiente de transferencia radiativo
  // entre placas planas
  // Entra temperaturas en  Celcius
  float s;
  s=SIGMA;
  T1=T1+273;
  T2=T2+273;
//  fprintf(log,"hradplacasp: s=%g \n",s);
//  fprintf(log,"hradplacasp: %f \n",epsilon1);
//  fprintf(log,"hradplacasp: %f \n",epsilon2);

  return SIGMA*(T1 +T2)*(T1*T1+ T2*T2)/(1./epsilon1 + 1./epsilon2 -1);
}
double hradcielo(double epsilonc, double Ts, double Tc, double Ta)
{
  // Coeficiente de transferencia radiativo
  // entre placas la placa y el cielo 
  // Entra temperaturas en  Celcius
  // Ts : temperatura del cielo
  // Tc : temperatura de la placa
  // Ta : temperatura ambiente
  // epsilonc : emisividad de la cubierta
  
  Ts=Ts+273;
  Tc=Tc+273;
  Ta=Ta+273;
//  fprintf(log,"hradcielo: %g \n",SIGMA);
//  fprintf(log,"hradcielo: %f \n",epsilonc);
  
  return SIGMA*epsilonc*(Tc + Ts)*(Tc*Tc+ Ts*Ts)*(Tc -Ts)/(Tc - Ta);
}
double hviento(double V, double L)
{
  // Coeficiente de transferencia convectivo
  // debido al viente
  // L : longitud característica
  // V : velocidad del viento

  // MEJORAR
  
  return 5.;
}
/*
float h_entreplacas(int i, float beta, float T[],float ep[], float D)
{
  // Coeficiente de transferencia convectivo + radiativo 
  // entre placas de una caja inclinada
  // Entra temperaturas en  Celcius
  // 
  // beta inclinacion en grados
  
fprintf(log,"\n h_entreplacas: %d-%d \n ",i,i+1);  
int N,itera;
float  Pr, Ra, Nu;
float V,L;
float hc,hr,hw,hs,h,Ut,R;
float ka;
float Tp, Ta, dT, dTi, Tm, Ts;

dTi=T[i]-T[i+1];
Tm=(T[i+1]+T[i])/2.;
Ra=Ra_aire(Tm,dTi,D);
Pr=Pr_aire(Tm);
Nu=Nu_Hollands(Ra,Pr,beta);
ka=k_aire(Tm);
fprintf(log,"Tm: %f \n",Tm);
fprintf(log,"dTi: %f \n",dTi);
fprintf(log,"Ra: %g \n",Ra);
fprintf(log,"Pr: %f \n",Pr);
fprintf(log,"Nu: %f \n",Nu);
fprintf(log,"ka: %f \n",ka);

hc=Nu*ka/D;
hr=hradplacasp(ep[i],ep[i+1], T[i],T[i+1]); 
fprintf(log,"hc: %f \n",hc);
fprintf(log,"hr: %f \n",hr);

h= (hc+hr);
fprintf(log,"h: %f  R_i=%f \n",h, 1./h);
return h;  
}
*/
float h_entreplacas2_log(FILE * log,int i, float beta,float T[],float ep[], float qr[], float qc[],float q[],float D)
{
  // Coeficiente de transferencia convectivo + radiativo 
  // entre placas de una caja inclinada
  // Entra temperaturas en  Celcius
  // 
  
fprintf(log,"\n -------------------------\n h_entreplacas2_log %d-%d\n ",i,i+1);  
int N,itera;
float  Pr, Ra, Nu;
float V,L;
float hc,hr,hw,hs,h,Ut,R;
float ka;
float Tp, Ta, dT, dTi, Tm, Ts;

// Hacer una funcion
dTi=T[i]-T[i+1];
Tm=(T[i+1]+T[i])/2.;
Ra=Ra_aire(Tm,dTi,D);
Pr=Pr_aire(Tm);
Nu=Nu_Hollands(Ra,Pr,beta);
ka=k_aire(Tm);
fprintf(log,"Tm: %f \n",Tm);
fprintf(log,"dTi: %f \n",dTi);
fprintf(log,"Ra: %g \n",Ra);
fprintf(log,"Pr: %f \n",Pr);
fprintf(log,"Nu: %f \n",Nu);
fprintf(log,"ka: %f \n",ka);

hc=Nu*ka/D;
hr=hradplacasp(ep[i],ep[i+1], T[i],T[i+1]); 

fprintf(log,"hc: %f \n",hc);
fprintf(log,"hr: %f \n",hr);

h= (hc+hr);
fprintf(log,"h: %f \n",h);

qc[i]=hc*(T[i]-T[i+1]);
qr[i]=hr*(T[i]-T[i+1]);
q[i]=h*(T[i]-T[i+1]);
fprintf(log,"h: %f  R_i=%f \n",h, 1./h);
fprintf(log,"h entre placas  : %f \n",h);
fprintf(log,"h[%d-%d]   : %f    R[%d-%d]  %f\n",i,i+1,h, i,i+1, 1./h);
fprintf(log,"flujos:  %d qc=%f qr=%f q=%f \n",i, qc[i],qr[i],q[i]);
fprintf(log,"\n Fin h_entreplacas2_log: %d-%d \n ",i,i+1);  
fprintf(log,"---------------------------------\n");
return h;
}

float hentrep_log(FILE * log, int i, float beta, float T[],float ep[], float D)
{
  // Coeficiente de transferencia convectivo + radiativo 
  // entre placas de una caja inclinada
  // Entra temperaturas en  Celcius
  // 
  // beta inclinacion en grados
  
fprintf(log,"\n ----------------------------\n hentrep_log: %d-%d \n ",i,i+1);  
int N,itera;
float  Pr, Ra, Nu;
float V,L;
float hc,hr,hw,hs,h,Ut,R;
float ka;
float Tp, Ta, dT, dTi, Tm, Ts;

dTi=T[i]-T[i+1];
Tm=(T[i+1]+T[i])/2.;
Ra=Ra_aire(Tm,dTi,D);
Pr=Pr_aire(Tm);
Nu=Nu_Hollands(Ra,Pr,beta);
ka=k_aire(Tm);
fprintf(log,"Tm: %f \n",Tm);
fprintf(log,"dTi: %f \n",dTi);
fprintf(log,"Ra: %g \n",Ra);
fprintf(log,"Pr: %f \n",Pr);
fprintf(log,"Nu: %f \n",Nu);
fprintf(log,"ka: %f \n",ka);

hc=Nu*ka/D;
hr=hradplacasp(ep[i],ep[i+1], T[i],T[i+1]); 
fprintf(log,"hc: %f \n",hc);
fprintf(log,"hr: %f \n",hr);

h= (hc+hr);
fprintf(log,"h: %f  R_i=%f \n",h, 1./h);
fprintf(log,"h entre placas  : %f \n",h);
fprintf(log,"h[%d-%d]   : %f    R[%d-%d]  %f\n",i,i+1,h, i,i+1, 1./h);
fprintf(log,"\n Fin hentrep_log: %d-%d \n ",i,i+1);  
fprintf(log,"---------------------------------\n");
return h;  
}
  
