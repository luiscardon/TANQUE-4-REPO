#include "param.h"
#include "uEstructuras.h"
#include "uEstimacionRad.h"
#include "uColectores.h"
#include "uCubiertas.h"
#include "uSSS.h"
#include "uImprime.h"

//  parcial-ERII-opticas-2c-2020.c
//  en TANQUES-4

// p1-parcial2015.c // en TANQUES-3
// pru-transmitancia-aMN-log.c  28-abril  2015
// Transmitancia de N x M cubiertas con absorcion


main()
{
  

//  datos:
//        angulo de incidencia  theta1, grados
//        indice de refraccion  n 1.526
FILE * log, * logT;
log = fopen("salta-optica.log", "w");
logT = fopen("salta-opticaT.log", "w");

fprintf(log,"ENERGIAS RENOVABLES II - LUS CARDON 2020 \n");
fprintf(log,"         PRIMER PARCIAL 2020 \n");
fprintf(log," Calcular los parametros opicos de una cubierta doble.\n");

fprintf(log,"Programa parcial-ERII-opticas-2c-2020.c \n");
fprintf(log,"En TANQUES-4. Compilado con LIBSE. \n");
fprintf(log,"CALCULOS Y RESULTADOS\n");
fprintf(log,"Arvhivo: salta-optica.log\n");


int N,M,i;
float  T=0., R=0.;
float beta;
float theta1, theta_d;
float theta2;
float nm,nn, Ln,Lm,Kn,Km;
float taua, tauav;
float pta,ptad,ptag;
float rhod;
float alpha, alphan;
float KLm,KLn,KL;
float TT,TTg,TTd;
int ncub;
float thetaeffd,thetaeffg;
ncub=2;
N=1;

nm=1.526;
Lm=0.004;
Km=32.;
KLm=Lm*Km;

M=1;
nn=1.526;
Ln=0.003;
Kn=4.;
KLn=Ln*Kn;

beta=24.7;
theta1= 5.01;

alphan=0.9; // absortanciade laplaca a incidencia nomal.
alpha=alpha_alphan(theta1,  alphan);

theta_d=60;
fprintf(log,"Datos Cubierta \n");

fprintf(log,"beta         :%f\n",beta);
fprintf(log,"N cubiertas: :%d\n",ncub);


fprintf(log,"\nCubierta N \n");
fprintf(log,"n             :%f\n",nn);
fprintf(log,"L             :%f\n",Ln);
fprintf(log,"K             :%f\n",Kn);
fprintf(log,"KL            :%f\n",KLn);

fprintf(log,"\nCubierta M \n");
fprintf(log,"n             :%f\n",nm);
fprintf(log,"L             :%f\n",Lm);
fprintf(log,"K             :%f\n",Km);
fprintf(log,"KL            :%f\n",KLn);
KL=KLm+KLn;
fprintf(log,"\nsumKL         :%f\n",KL);
fprintf(log,"\ntheta1       :%f\n",theta1);
fprintf(log,"alpha_n       :%f\n",alphan);
fprintf(log,"alpha         :%f\n",alpha);

fprintf(log,"=================================\n");
fprintf(log,"=======Resultados================\n");
fprintf(log,"--------------------------------\n");
fprintf(log,"Calculo de la reflectividad rho_d\n");

theta_d=60.;  // angulo deincidencia efectivo para  transmitancia difusa
              // proveniente del laplaca absorbedora
logT=TyR_CC_log(logT, M, nm, Km, Lm, N , nn, Kn, Ln, theta_d, &T, &R);
log=TyR_CC_log(log, M, nm, Km, Lm, N , nn, Kn, Ln, theta_d, &T, &R);
fprintf(log,"--------------------------------\n");
theta2=theta_2(1., nm, theta_d);
tauav=tau_a((Lm+Ln), theta2,Km);
taua=tau_a2(KLm,KLn, theta2);
rhod=taua-T;
TTd=T;
fprintf(log,"theta_d        :%f\n",theta_d);
fprintf(log,"theta2         :%f\n",theta2);
fprintf(log,"T(theta_d)     :%f\n",T);
fprintf(log,"--------------------------------\n");
fprintf(log,"Ln + Lm        :%f\n",Ln+Lm);
fprintf(log,"K max          :%f\n",Km);
fprintf(log,"KL=(KLn+KLm)*  :%f\n",KL);
fprintf(log,"taua suma      :%f\n",tauav);
fprintf(log,"taua nm=nn     :%f\n",taua);
fprintf(log,"rhod           :%f\n",rhod);


fprintf(log,"Calculo de la reflectividad rho_d finalizado\n");
fprintf(log,"======================================================\n");
// ==========================================================

fprintf(log,"\n Calculo de Transmitacicia de la radiacion directa \n para el angulo de incidencia \n");

logT=TyR_CC_log(logT, M, nm, Km, Lm, N , nn, Kn, Ln, theta1, &T, &R);
log=TyR_CC_log(log, M, nm, Km, Lm, N , nn, Kn, Ln, theta1, &T, &R);

TT=T;
theta2=theta_2(1., nm, theta_d);
tauav=tau_a((Lm+Ln), theta2,Km);
taua=tau_a2(KLm,KLn, theta2);

// Producto transmitancia abortancia directa
pta=ptaualpha(T, alpha, rhod);


fprintf(log,"theta1        :%f\n",theta1);
fprintf(log,"theta2         :%f\n",theta2);
fprintf(log,"T(%f)     :%f\n",theta1, T);
fprintf(log,"--------------------------------\n");
fprintf(log,"Ln + Lm        :%f\n",Ln+Lm);
fprintf(log,"K max          :%f\n",Km);
fprintf(log,"KL=(KLn+KLm)*  :%f\n",KL);
fprintf(log,"taua suma      :%f\n",tauav);
fprintf(log,"taua nm=nn     :%f\n",taua);
fprintf(log,"rhod           :%f\n",rhod);
fprintf(log,"alphan         :%f\n",alphan);
fprintf(log,"alpha          :%f\n",alpha);
fprintf(log,"--------------------------------\n");
fprintf(log,"producto (tau alpha)_beam\n");
fprintf(log,"\n ----Beam - Directa ----------\n");
fprintf(log,"theta1         :%f\n",theta1);
fprintf(log,"Td(%f)     :%f\n",theta1, T);
fprintf(log,"(tau alpha)b   :%f\n",pta);

fprintf(log,"Calculo del producto (tau alpha)_b  finalizado\n");
fprintf(log,"======================================================\n");

fprintf(log,"Calculo del producto (tau alpha)_b  simplificado\n");
fprintf(log,"Dos cubiertas iguales\n");

log=ta_log(log, 2, KLm+KLn ,  theta1,  alpha, nn, &pta);
fprintf(log,"(tau alpha)b   :%f\n",pta);

fprintf(log,"Calculo del producto (tau alpha)_b  simplificado finalizado\n");
fprintf(log,"======================================================\n");



fprintf(log,"\n Calculo de Transmitacicia de la radiacion difuse \n para el angulo efectivo \n");


// Producto transmitancia abortancia difusa
thetaeffd=theta_d_eff(beta);
logT=TyR_CC_log(logT, M, nm, Km, Lm, N , nn, Kn, Ln, thetaeffd, &T, &R);
log=TyR_CC_log(log, M, nm, Km, Lm, N , nn, Kn, Ln, thetaeffd, &T, &R);
TTd=T;
//TTd=T;
ptad=ptaualpha(TTd, alpha, rhod);
fprintf(log,"thetaeffd        :%f\n",thetaeffd);
fprintf(log,"theta2         :%f\n",theta2);
fprintf(log,"T(%f)     :%f\n",thetaeffd,T);
fprintf(log,"--------------------------------\n");
fprintf(log,"Ln + Lm        :%f\n",Ln+Lm);
fprintf(log,"K max          :%f\n",Km);
fprintf(log,"KL=(KLn+KLm)*  :%f\n",KL);
fprintf(log,"taua suma      :%f\n",tauav);
fprintf(log,"taua nm=nn     :%f\n",taua);
fprintf(log,"rhod           :%f\n",rhod);
fprintf(log,"alphan         :%f\n",alphan);
fprintf(log,"alpha          :%f\n",alpha);
fprintf(log,"--------------------------------\n");
fprintf(log,"producto (tau alpha)_beam\n");
fprintf(log,"\n ----Beam - Directa ----------\n");
fprintf(log,"theta1         :%f\n",thetaeffd);
fprintf(log,"Td(%f)     :%f\n",thetaeffd, T);
fprintf(log,"(tau alpha)b   :%f\n",pta);

fprintf(log,"Calculo del producto (tau alpha)_D  finalizado\n");
fprintf(log,"======================================================\n");

// Producto transmitancia absortancia ground
thetaeffg=theta_g_eff(beta);
logT=TyR_CC_log(logT, M, nm, Km, Lm, N , nn, Kn, Ln, thetaeffg, &T, &R);
ptag=ptaualpha(TTg, alpha, rhod);
TTg=T;
fprintf(log,"thetaeffg        :%f\n",thetaeffg);
fprintf(log,"theta2         :%f\n",theta2);
fprintf(log,"T(%f)     :%f\n",thetaeffd,T);
fprintf(log,"--------------------------------\n");
fprintf(log,"Ln + Lm        :%f\n",Ln+Lm);
fprintf(log,"K max          :%f\n",Km);
fprintf(log,"KL=(KLn+KLm)*  :%f\n",KL);
fprintf(log,"taua suma      :%f\n",tauav);
fprintf(log,"taua nm=nn     :%f\n",taua);
fprintf(log,"rhod           :%f\n",rhod);
fprintf(log,"alphan         :%f\n",alphan);
fprintf(log,"alpha          :%f\n",alpha);
fprintf(log,"--------------------------------\n");
fprintf(log,"producto (tau alpha)_beam\n");
fprintf(log,"\n ----Beam - Directa ----------\n");
fprintf(log,"theta1         :%f\n",thetaeffg);
fprintf(log,"T(%f)     :%f\n",thetaeffg,T);
fprintf(log,"(tau alpha)b   :%f\n",pta);

fprintf(log,"Calculo del producto (tau alpha)_G  finalizado\n");
fprintf(log,"======================================================\n");


fprintf(log,"RESUMEN\n");

fprintf(log,"Reflectividad rho_d\n");
fprintf(log,"theta_d     :%f\n",theta_d);
fprintf(log,"theta2     :%f\n",theta2);
fprintf(log,"T(%f)          :%f\n",theta_d,TT);
fprintf(log,"--------------------------------\n");
fprintf(log,"Ln + Lm        :%f\n",Ln+Lm);
fprintf(log,"K max          :%f\n",Km);
fprintf(log,"KL=(KLn+KLm)   :%f\n",KL);
fprintf(log,"taua suma      :%f\n",tauav);
fprintf(log,"taua nm=nn     :%f\n",taua);
fprintf(log,"rhod           :%f\n",rhod);
fprintf(log,"--------------------------------\n");


fprintf(log,"alphan         :%f\n",alphan);
fprintf(log,"alpha          :%f\n\n",alpha);
fprintf(log,"Producto (tau alpha)i\n\n");
fprintf(log,"\n ----i : Directa\n");
fprintf(log,"theta1         :%f\n",theta1);
fprintf(log,"Td(theta1)     :%f\n",TTd);
fprintf(log,"(tau alpha)b   :%f\n",pta);
fprintf(log,"\n ----i : Difusa\n");

fprintf(log,"thetaeffd      :%f\n",thetaeffd);
fprintf(log,"Td(%f)  :%f\n",thetaeffd,TTd);
fprintf(log,"(tau alpha)d   :%f\n",ptad);
fprintf(log,"\n ----i : Ground\n");
fprintf(log,"thetaeffg      :%f\n",thetaeffg);
fprintf(log,"Tg(%f)  :%f\n",thetaeffg, TTg);
fprintf(log,"(tau alpha)g   :%f\n",ptag);

return 0;
}
