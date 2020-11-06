#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "param.h"
#include "uEstructuras.h"
#include "uPropiedades.h"
#include "uPropFisicas.h"
#include "uCorrelaciones.h"

#define SIGMA 5.67e-8 // cosntante de Stefan Boltzaaman, WM^-2K^-4


double hradplacasp(double epsilon1,double epsilon2, double T1,double T2);
double hradcielo(double epsilonc, double Ts, double Tc, double Ta);
double hviento(double V, double L);
//float  h_entreplacas(int i, float beta, float T[],float ep[],float D);
float  h_entreplacas2_log(FILE * log,int i, float beta, float T[],float ep[], float qr[],float qc[],float q[],float D);
float  hentrep_log(FILE * log, int i, float beta, float T[],float ep[], float D);

int main()
{
  // R-corto.c  
  // Calculo de coeficiente de transferencia de tubo-suelo
  // se compila con libSE
FILE * log;
log = fopen("R-corto.log", "w");

fprintf(log,"===========UR-corto.c ====================\n");
fprintf(log,"LUS CARDON 2020 \n");
fprintf(log," Calcular el oeficiente global de perdidas.\n");

fprintf(log,"Programa R-corto.c \n");
fprintf(log,"En UTILITARIOS/.../U0-11. Compilado con LIBSE. \n");
fprintf(log,"CALCULOS Y RESULTADOS\n");
fprintf(log,"Archivo: R-corto.log\n");

/*
Datos
N   numero de cubiertas
*/
float Pr, Re, Nu, D;
float Q,V,L, A;
float hc,R;
float ka;
float T,TT;
float nu,nui,mu,alpha,rho,Cp,gban,k;
// D    :  Diametro del tubo, m
// L    :  Longitud, m
// V    :  velocidad media m/s
// Ta   :  temperatura del aire, C

// calculos intermedios
// hc  : coeficiente convectivo, W/m^2 
// ka  : conductividad del aire, W/mC

// variables auxiliares
// dT   
D=0.1;
Pr=0.7;

nu=15.E-6;
T=27;
Q=243./3600;

A=M_PI*D*D/4.;
V=Q/A;

Re=V*D/nu;


rho=rho_aire(T); //Ok
mu=mu_aire(T); //ok
nu=nu_aire(T); // ok
k=k_aire(T);  // ok
Cp=Cp_aire(T); // ok
gban=gban_aire(T); //ok
alpha=alpha_aire(T); //ok
Pr=Pr_aire(T);

for (TT=200 ; TT<=  400 ; TT=TT+10.)
{
  T=TT-273.15;
rho=rho_aire(T); //Ok
mu=mu_aire(T); //ok
nu=nu_aire(T); // ok
k=k_aire(T);  // ok
Cp=Cp_aire(T); // ok
gban=gban_aire(T); //ok
alpha=alpha_aire(T); //ok
fprintf(stderr,"%3.0f %1.4f %2.2f %2.2f %4.1f %3.1f %3.1f %2.2f \n",TT, rho, mu*1.E6, k*1.E3, Cp, 1.*1.E-3/nu ,gban*1.E-6, alpha*1.E6);  
}
Nu=Nu_Gnielinski(Re,Pr);

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
  
