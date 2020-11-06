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
	// R-tubo.c 
	// i7/TANQUEs-4/TANQUE-4-REPO  5 nov 2020
	// Calcula el coeficiente de perdidas al suelo de unntubo
	// usa las librerias  modificadas y actualizadas
	//"uPropiedades.h"
	//"uPropFisicas.h"
	// + usa las estructura: struc propiedades pro
	// se compila con libSE
	FILE * log;
	log = fopen("R-tubo.log", "w");
	
	fprintf(log,"===========R-tubo.c ====================\n");
	fprintf(log,"LUS CARDON 2020 \n");
	fprintf(log," Prueba de funcion Gnielisnsky y estructuras .\n");
	
	fprintf(log,"Programa R-tubo.c \n");
	fprintf(log,"En i7/TANQUEs-4/TANQUE-4-REPO. Compilado con LIBSE. \n");
	fprintf(log,"Archivo: R-tubo.log\n");
	
	float T,TT;
	float rho,mu,nu,k,Cp,gban,alpha,Pr;
	float ks,kp;
	float Re, Nu, h;
	float D,Di,ri,re,e,Q,V,A,L;
	float R,Rc,Rp,Rs;
	typedef struct propiedades pro;
	pro aire, suelo;
	
	T=27.;
	Di=0.11;  //Tigre Ramat pluvial IRAM
	e=0.0022;
	ri=Di/2;
	re=ri+e;
	L=4;
	A=M_PI*Di*Di/4.;
	Q=243./3600;
	V=Q/A;
	
	
	propiedades_aire(T, &aire );
	rho=aire.rho; 
	mu=aire.mu; 
	nu=aire.nu; 
	k=aire.k;  
	Cp=aire.Cp; 
	gban=aire.gban; 
	alpha=aire.alpha; 
	Pr=aire.Pr;
	
	kp=k_pvc(T);
	
	printf("T     %f  C    \n",T);
	printf("T     %f  K    \n",T+273.15);
	printf("Di    %f  m    \n",Di);
	printf("Q     %f  m^3/s\n", Q);
	printf("V     %f  m/s  \n",V);
	printf("k     %g  m    \n",k);
	printf("Cp    %f  m    \n",Cp);
	
	printf("k pvc %f  m    \n",kp);
	
	p_suelo( &suelo );
	ks=k.suelo;
	printf("k suelo %f  m    \n",ks);
	
	Re=V*D/nu;
	Nu=Nu_Gnielinski(Re,Pr);
	h=Nu*k/Di;
	Rc=1./(M_PI*Di*h);
	
	printf("Re   %f       \n",Re);
	printf("Nu   %f       \n",Nu);
	printf("h    %f\n",h);
	printf("Rc   %f\n",Rc);
	
	Rp=flog( (ri+e)/ri)/(2.*M_PI*kp);
	printf("Rp   %f\n",Rp);
	
	Rs=flog( (2.*re)/(ri+e) )/(2*M_PI*ks);
	printf("Rs   %f\n",Rs);
	R=Rc +Rs + Rp;
	printf("R   %f\n",R);
	
	return 0;
}

