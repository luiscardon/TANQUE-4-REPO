// UTILITARIOS - uColectores.c

// Luis Cardon 20 de setiembre de 2013
// derivado de fpp.c
#include "param.h"
#include "uEstructuras.h"
#include "uEstimacionRad.h"
#include "uColectores.h"
#include "uCubiertas.h"
#include "uSSS.h"
#include "uImprime.h"
#include "uCorrelaciones.h"
#include "uPropFisicas.h"



double angle_factor_man(double a, double b, double c);
// ------------------------------------------------------------------
// ==================================================================
float fpp(float ccr)
{
// Retorna el factor de flujo del colector
// Enra ccr, razon de capacitancia del colector
// Duffie y Beckman, 1980, pp 224.
return ccr*(1.-exp(-1./ccr));
}
// ------------------------------------------------------------------

float Qu(float Ac,float Fr, float S, float UL, float Ti, float Ta)
{
// Retorna el calor util
// Entra Ac, Area del colector
// Entra Fr, factor de remocion
// Entra S,
// Entra UL,
// Entra Ti,
// Entra Ta,

return Ac*Fr *(S-UL*(Ti-Ta));
}
float ooQuW(struct s_colector_g co, float S,float Ti, float Ta)
{
// Retorna el calor util
// Entra Ac, Area del colector
// Entra Fr, factor de remocion
// Entra S,
// Entra UL,
// Entra Ti,
// Entra Ta,
return co.Ac*co.FR*(S - co.UL*(Ti-Ta));
}
float ooQuJ(struct s_colector_g co, float S, float Ti, float Ta)
{
// Retorna el calor util
// Entra Ac, Area del colector
// Entra Fr, factor de remocion
// Entra S,
// Entra UL,
// Entra Ti,
// Entra Ta,


return co.Ac*co.FR*(S - co.UL*3600*(Ti-Ta));
}

// ------------------------------------------------------------------
/*
float Fp(float UL, float W, float D,float Cb,float Di,float hfi)
{
// Duffie y Beckman, 1980, pp 216.

// Retorna el factor de eficiencia de la placa
// Entra m, sqrt(UL/(k*delta))
// Entra W, distancia entre tubos
// Entra D, diámetro de los tubos
float a;
a=m*(W-D)/2.;
// ------------------------------------------------------------------

return (1./UL)/(W( 1/(UL*(D+(W-D)*F)+1./Cb+1/(PI*Di+hfi) )); 
}
*/
float Fr(float Fp, float Fpp)
{
// Duffie y Beckman, 1980, pp 216.

// Retorna el factor de remocion
// Entra factor de eficiencia de aleta, Fp
// Entra factor de  flujo, Fpp

return Fp*Fpp;
}
// ------------------------------------------------------------------

float S(float Rb, float Ib, float Id, float Ig, float ptaoalfab  ,
	float ptaoalfad  ,float ptaoalfag, float rhod, float beta)
{
    float cb,A,B;

  // Luis cardon - 24 de setiembre de 2013

  // Duffie y Beckman, 1980, pp 187.
// Verificada

// Retorna radiación incidente
// Entra Rb
// Entra Ib
// Entra Id
// Entra Ig
// Entra ptaoalfab
// Entra ptaoalfad
// Entra ptaoalfag
// Entra rhod
// Entra beta, 

beta=gradosarad(beta);
cb=cos(beta);
A=(1.+cb)*0.5;
B=(1.-cb)*0.5;

return Ib*Rb*ptaoalfab+ Id*ptaoalfad *A +rhod*(Ib+Id)*ptaoalfad *B; 
}
// ------------------------------------------------------------------
float SS(float Ib, float Id, float Ig, float ptaoalfab  ,
	float ptaoalfad  ,float ptaoalfag)
{
  // Luis cardon - 23 de setiembre de 2013
  
  // Utiliza los datos discriminados de I_Tdis
  // I_Tdis da los valores sobre superficie inclinada
  // Duffie y Beckman, 1980, pp 187.

// Retorna radiación incidente
// Entra Rb
// Entra Ib
// Entra Id
// Entra Ig
// Entra ptaoalfab
// Entra ptaoalfad
// Entra ptaoalfag
// Entra rhod

return Ib*ptaoalfab+ Id*ptaoalfad  + Ig*ptaoalfad ; 
}
// ------------------------------------------------------------------
float SSdis(float * Ib, float * Id, float * Ig, float ptaoalfab  ,
	float ptaoalfad  ,float ptaoalfag)
{
  // Luis cardon - 23 de setiembre de 2013

  // Archivo: uColectores.c
  
  // Utiliza los datos discriminados de I_Tdis
  // I_Tdis ya tiene los valores corregidos para superficie inclinada
  // Duffie y Beckman, 1980, pp 187.

// Retorna:
//         radiación absorbida sobre superficie inclinada.
//         radiacion directa, difusa y albedo absorbido  

  // Entra Ib : Ib*Rb
// Entra Id : Id*(1+cos beta)/2
// Entra Ig : Ig*(1-cos beta)/2
// Entra ptaoalfab : producto transmitancia absortancia
// Entra ptaoalfad : producto transmitancia absortancia
// Entra ptaoalfag : producto transmitancia absortancia
// Entra rhod

*Ib=*Ib*ptaoalfab;
*Id=*Id*ptaoalfad;
*Ig=*Ig*ptaoalfad;

return *Ib*ptaoalfab+ *Id*ptaoalfad  + *Ig*ptaoalfad ; 
}

// ------------------------------------------------------------------

float taoalfa (float tao, float alfa, float rhod, float n)
{
  // Luis Cardon 23 de setiembre de 2013
// Duffie y Beckman, , pp  216.

  // Retorna el producto transitancia absortancia del colector
// Entra n, numero de cubiertas
// Entra tao
// Entra alfa
// Entra rhod

return tao*alfa/(1.-(1.-alfa)*rhod);

}
// ------------------------------------------------------------------

// ------------------------------------------------------------------
// Aletas------------------------------------------------------------------
// ------------------------------------------------------------------

float Faleta(float U_L,float k, float delta, float W, float D)
{// retorna la eficiencia de aleta F
// Duffie Beckman pp. 258
float m ,aux;
m= sqrt(U_L/(k*delta));
aux=m*(W-D)/2.;
return tanh(aux)/aux;
}
float Faleta_kdelta(float U_L,float kd ,float W, float D)
{// retorna la eficiencia de aleta F
// Duffie Beckman pp. 258
float m ,aux;
m= sqrt(U_L/(kd));
aux=m*(W-D)/2.;
return tanh(aux)/aux;
}

float Faleta_mm(float mm)
{// retorna la eficiencia de aleta F
 // Duffie Beckman pp. 258
return tanh(mm)/mm;
}
float mm_aleta(float U_L,float k, float delta,float W, float D)
{
// retorna m
return  sqrt(U_L/(k*delta))*(W-D)/2.;
}
// ------------------------------------------------------------------
// Parametros de perfomance del coelector----------------------------
// ------------------------------------------------------------------

float Fprima(float F, float U_L,float W, float D, float D_i, float C_b, float h_fi)
{
// Factor de eficiencia del colector
// Colector de placa plana estandar
// Dufie y Beckman,  pp. 259, eq. 6.518
float R_L,R1,R2,R3;
R_L=1./U_L;
R1=1./(U_L*  ( D+ (W-D)*F ) );
R2=1./C_b;
R3=1./(PI*D_i*h_fi);
return R_L/( W*(R1+R2+R3) ) ;
}

// ------------------------------------------------------------------

float Col_cap(float mpunto, float Cp_fluido, float colA, float U_L, float Fp)
{// Capacitancia de Colector
return mpunto*Cp_fluido/(colA*U_L*Fp);
}
// ------------------------------------------------------------------
float F_pp(float colCap)
{
// Factor de flujo del coelctro
return colCap*(1. -exp(- 1./colCap));
}
// ------------------------------------------------------------------
float F_R(float Fp,float Fpp)
{
// Factor de remocion del colectro
return Fpp*Fp;
}
// ------------------------------------------------------------------

float T_fo(float Ta, float Tfi, float S, float U_L, float colCap )
{ 
// Temperatura de salida del colector
return (Tfi -Ta - S/U_L)* exp(-1./colCap) + Ta + S/U_L;
}
float Q_u(float colA, float FR, float S, float U_L, float Ti, float Ta)
{
// Ecuacion del colector
return colA*FR*(S-U_L*(Ti-Ta));
}
float Fp_colector_tipo_a(float W, float D, float h, float C_b ,float U_L, float F)
{
// colector tipo a, Duffie y Beckman pp 278
float Fp, R1, R2, R3;

R1= W*U_L /(PI*D*h);
R2=  W*U_L /C_b;
R3= W/(D +(W-D)*F);
Fp=1. /(R1+R2+R3);
return Fp;
}
// ------------------------------------------------------------------
float Fp_colector_tipo_b(float W, float D, float h, float C_b ,float U_L, float F)
{
// colector tipo a, Duffie y Beckman pp 278
float Fp, R1, R2, R3;

R1= W*U_L /(PI*D*h);
R2=  W*U_L /C_b;
R3= W/((W-D)*F);
Fp=1. /(R1 + (1. / (D/W + 1. / (R2+R3) ) ) );
return Fp;
}
// ------------------------------------------------------------------

float Fp_colector_tipo_c(float W, float D, float h, float C_b ,float U_L, float F)
{
// colector tipo a, Duffie y Beckman pp 278
float Fp, R1, R2, R3;

R1= W*U_L /(PI*D*h);
//R2=  W*U_L /C_b;
R3= W/((W-D)*F);
Fp=1. /(R1 + R3);
return Fp;
}/*
float Fp_colector_tipo_b(float W, float D, float h, float C_b ,float U_L, float F)
{
// colector tipo a, Duffie y Beckman pp 278
float Fp, R1, R2, R3;

R1= W*U_L /(PI*D*h);
//R2=  W*U_L /C_b;
R3= W/((W-D)*F);
Fp=1. /(R1 + R3);
return Fp;
}
*/
// ------------------------------------------------------------------
// Factores de forma-------------------------------------------------
// ------------------------------------------------------------------
double angle_factor_man(double a, double b, double c)
{
// Factor angular entre una superficie rectangular y una persona
//
  double X,Y,F;
c=1.8*c;
X=a/c;
Y=b/c;
//printf("X %f Y 5F\n",X,Y);
F=1./(4.*PI) * ( atan(Y/sqrt(1.+X*X)) * X/sqrt(1.+X*X) +
                  atan(X/sqrt(1.+Y*Y)) * Y/sqrt(1.+Y*Y) );
return F;
}
// ------------------------------------------------------------------
// Coeficientesde transferencia cubiertas----------------------------
// ------------------------------------------------------------------

double  Ut_klein80(int N, double Tp, double Ta, double epsilonp,
		   double epsilong, double hw, double beta )
{
// Coeficeinte global de transferencai de calor por la cubiertas
// Formula de Klein 80, DBnpp 211
// para beta >70, usar 70
double C,f, e, Ut,Uta,Utb;
double aux1,aux2, aux3;

if (beta > 70)beta = 70;
Ta=Ta+273.15;
Tp=Tp+273.15;
C=520.*(1.-0.000051*beta*beta);
f=(1.+ 0.089*hw - 0.1166*hw*epsilonp)*(1. + 0.07866*N);
e=0.43*(1.- 100./Tp);
aux1=C/Tp;
aux2= pow((Tp-Ta)/(N+f),e);
aux3=1./hw;
Uta= 1./ (N/(aux1*aux2) + aux3);
Utb=5.67e-8* (Tp*Tp +Ta*Ta) * (Tp+Ta) /( (1./ (epsilonp+0.00591*N*hw)) + ( (2.*N +f -1. +0.133*epsilonp)/epsilong ) - N);
Ut=Uta+Utb;
//printf("%f %f %f %f %f %F \n", f,C,e, Uta,Utb,Ut );
return Ut;  
}
// ------------------------------------------------------------------
double  Ut_agarwal81(int N, double Tp, double Ta, double epsilonp, double epsilong, double hw, double beta )
{
// Coeficeinte global de transferencai de calor por la cubiertas
// Formula de Agarwal y Larson  81, 
double C,f, e, Ut,Uta,Utb;
double aux1,aux2, aux3;

Ta=Ta+273.15;
Tp=Tp+273.15;
C=250.*(1.-0.0044*(beta-90));
f=(1.- 0.04*hw + 0.005*hw*hw)*(1. + 0.091*N);
//e=0.43*(1.- 100./Tp);
aux1=C/Tp;
aux2= pow((Tp-Ta)/(N+f),0.33);
aux3=1./hw;
Uta= 1./ (N/(aux1*aux2) + aux3);
Utb=5.67e-8* (Tp*Tp +Ta*Ta) * (Tp+Ta) /( (1./ (epsilonp+0.05*N*(1.-epsilonp)) ) + ( (2.*N + f -1.)/epsilong ) - N);
Ut=Uta+Utb;
//printf("%f %f %f %f %f %F \n", f,C,e, Uta,Utb,Ut );
return Ut;  
}
float litrosampunto(float litrosmin)
{ // devuelve kg/seg
 return litrosmin/60.; 
}
/*
float Col_cap(float mpunto, float Cp_fluido, float colA, float U_L, float Fp)
{// Capacitancia de Colector
return mpunto*Cp_fluido/(colA*U_L*Fp);
}
*/
float maleta(float UL,float k,float delta)
{
 return sqrt(UL/(k*delta)); 
}

// colector.c
float Tcolector(float Tcfi, float Ta, float S, struct s_colector_g micol)
{
  // respuesta instantanea del colector
  // S potencia instantanea, W/m^2
  float Tcfo;
  float FR,Ac,mp,Cp,UL;
  
  FR=micol.FR;
  Ac=micol.Ac;
  mp=micol.mpunto;
  Cp=micol.Cp;
  UL=micol.UL;
  
  Tcfo=Tcfi + FR*Ac*(S-UL*(Tcfi-Ta))/(mp*Cp);
  
 return Tcfo; 
} 
// ============================================
// colector-aux.h

void imprime_colector(struct s_colector_g a)
{ float rc;
 rc=a.mpunto*a.Cp/(a.Ac*a.UL*a.Fp);
 fprintf(stdout,"# Colector tipo\n");
 fprintf(stdout,"# Area del colector                     %f\n",a.Ac); 
 fprintf(stdout,"# Coeficiente global de pérdidas        %f\n",a.UL); 
 fprintf(stdout,"# Flujo másico                          %f\n",a.mpunto); 
 fprintf(stdout,"# Calor específico del fluido           %f\n",a.Cp); 
 fprintf(stdout,"# Eficiencia de aleta colector, F       %4.2f\n",a.F); 
 fprintf(stdout,"# Factor de eficiencia del colector, F' %4.2f\n",a.Fp); 
 fprintf(stdout,"# Factor de remoción, FR                %4.2f\n",a.FR); 
 fprintf(stdout,"# Factor de flujo del colector,  F\"     %4.2f\n",a.Fpp); 
 fprintf(stdout,"# Razon de capacitancia, mpuntoCp/AUF'  %f\n",rc);
 fprintf(stdout,"# dT activacion, dTON                   %f\n",a.dTON);
 fprintf(stdout,"# dT apagado, dTOFF                     %f\n",a.dTOFF);
 fprintf(stdout,"# Producto (tao alpha)                  %f\n",a.pta);
 fprintf(stdout,"# =========================================\n");
 
}
FILE * imprime_colector_log(FILE * log, struct s_colector_g a)
{ float rc;
 rc=a.mpunto*a.Cp/(a.Ac*a.UL*a.Fp);
 fprintf(log,"# Colector tipo\n");
 fprintf(log,"# Area del colector                     %f\n",a.Ac); 
 fprintf(log,"# Coeficiente global de pérdidas        %f\n",a.UL); 
 fprintf(log,"# Flujo másico                          %f\n",a.mpunto); 
 fprintf(log,"# Calor específico del fluido           %f\n",a.Cp); 
 fprintf(log,"# Eficiencia de aleta colector, F       %4.2f\n",a.F); 
 fprintf(log,"# Factor de eficiencia del colector, F' %4.2f\n",a.Fp); 
 fprintf(log,"# Factor de remoción, FR                %4.2f\n",a.FR); 
 fprintf(log,"# Factor de flujo del colector,  F\"     %4.2f\n",a.Fpp); 
 fprintf(log,"# Razon de capacitancia, mpuntoCp/AUF'  %f\n",rc);
 fprintf(log,"# dT activacion, dTON                   %f\n",a.dTON);
 fprintf(log,"# dT apagado, dTOFF                     %f\n",a.dTOFF);
 fprintf(log,"# Producto (tao alpha)                  %f\n",a.pta);
 fprintf(log,"# =========================================\n");
 
 return log;
}
void imprime_colector2(struct s_colector_g a,struct s_placa b)
{ float rc;
 rc=a.mpunto*a.Cp/(a.Ac*a.UL*a.Fp);
 fprintf(stdout,"# Colector tipo\n");
 fprintf(stdout,"# Area del colector                     %f\n",a.Ac); 
 fprintf(stdout,"# Coeficiente global de pérdidas        %f\n",a.UL); 
 fprintf(stdout,"# Flujo másico                          %f\n",a.mpunto); 
 fprintf(stdout,"# Calor específico del fluido           %f\n",a.Cp); 
 fprintf(stdout,"# Factor de aleta, m                    %4.3f\n",b.m); 
 fprintf(stdout,"# Eficiencia de aleta colector, F       %4.3f\n",a.F); 
 fprintf(stdout,"# Factor de eficiencia del colector, F' %4.3f\n",a.Fp); 
 fprintf(stdout,"# Factor de remoción, FR                %4.3f\n",a.FR); 
 fprintf(stdout,"# Factor de flujo del colector,  F\"    %4.3f\n",a.Fpp); 
 fprintf(stdout,"# Razon de capacitancia, mpuntoCp/AUF'  %f\n",a.rcC);
 fprintf(stdout,"# dT activacion, dTON                   %f\n",a.dTON);
 fprintf(stdout,"# dT apagado, dTOFF                     %f\n",a.dTOFF);
 fprintf(stdout,"# Producto (tao alpha)                  %f\n",a.pta);

 fprintf(stdout,"# =========================================\n");
 
}  

