#include "param.h"
#include "uEstructuras.h"
#include "uEstimacionRad.h"
#include "uColectores.h"
#include "uCubiertas.h"
#include "uSSS.h"
#include "uImprime.h"

void imprimegeo(sgeo geo)
{
printf("-------------------------------\n");
printf("geo.n       %d\n",geo.n);
printf("geo.delta   %f\n",geo.delta);
printf("geo.phi     %f\n",geo.phi);
printf("geo.gamma   %f\n",geo.gamma);
printf("geo.beta    %f\n",geo.beta);
printf("geo.altitud %f\n",geo.altitud);
printf("geo.clima   %d\n",geo.clima);
printf("geo.albedo  %f\n",geo.albedo);

}
void fimprimegeo(sgeo geo)
{

fprintf(stderr,"-------------------------------\n");
fprintf(stderr,"Nombre:      geo.nombre  %s\n",geo.nombre);
fprintf(stderr,"Dia   :      geo.n       %d [.]\n",geo.n);
fprintf(stderr,"Declinacion: geo.delta   %f [grados]\n",geo.delta);
fprintf(stderr,"Latitud:     geo.phi     %f [grados]\n",geo.phi);
fprintf(stderr,"Acimut:      geo.gamma   %f [grados]\n",geo.gamma);
fprintf(stderr,"Pendiente:   geo.beta    %f [grados]\n",geo.beta);
fprintf(stderr,"Altitud:     geo.altitud %f [km]\n",geo.altitud);
fprintf(stderr,"Clima        geo.clima   %d [.]\n",geo.clima);
fprintf(stderr,"Albedo       geo.albedo  %f [.]\n",geo.albedo);
}

void fimprimegeomas(sgeo geo, sgeohr geohr)
{

fprintf(stderr,"-------------------------------\n");
fprintf(stderr,"Nombre:      geo.nombre  %s\n",geo.nombre);
fprintf(stderr,"Dia   :      geo.n       %d [.]\n",geo.n);
fprintf(stderr,"Declinacion: geo.delta   %f [grados]\n",geo.delta);
fprintf(stderr,"Latitud:     geo.phi     %f [grados]\n",geo.phi);
fprintf(stderr,"Acimut:      geo.gamma   %f [grados]\n",geo.gamma);
fprintf(stderr,"Pendiente:   geo.beta    %f [grados]\n",geo.beta);
fprintf(stderr,"Altitud:     geo.altitud %f [km]\n",geo.altitud);
fprintf(stderr,"Clima        geo.clima   %d [.]\n",geo.clima);
fprintf(stderr,"Albedo       geo.albedo  %f [.]\n",geo.albedo);
fprintf(stderr,"Amanecer                 %f [hr]\n",tiempo_hr(geohr.omegas)-12);
}
void himprimegeomas(FILE * file, sgeo geo, sgeohr geohr)
{

fprintf(file,"-------------------------------\n");
fprintf(file,"Nombre:      geo.nombre  %s\n",geo.nombre);
fprintf(file,"Dia   :      geo.n       %d [.]\n",geo.n);
fprintf(file,"Declinacion: geo.delta   %f [grados]\n",geo.delta);
fprintf(file,"Latitud:     geo.phi     %f [grados]\n",geo.phi);
fprintf(file,"Acimut:      geo.gamma   %f [grados]\n",geo.gamma);
fprintf(file,"Pendiente:   geo.beta    %f [grados]\n",geo.beta);
fprintf(file,"Altitud:     geo.altitud %f [km]\n",geo.altitud);
fprintf(file,"Clima        geo.clima   %d [.]\n",geo.clima);
fprintf(file,"Albedo       geo.albedo  %f [.]\n",geo.albedo);
fprintf(file,"Amanecer                 %f [hr]\n",tiempo_hr(geohr.omegas)-12);
}

void imprimegeomas(sgeo geo, sgeohr geohr)
{

printf("-------------------------------\n");
printf("Nombre:      geo.nombre  %s\n",geo.nombre);
printf("Dia   :      geo.n       %d [.]\n",geo.n);
printf("Declinacion: geo.delta   %f [grados]\n",geo.delta);
printf("Latitud:     geo.phi     %f [grados]\n",geo.phi);
printf("Acimut:      geo.gamma   %f [grados]\n",geo.gamma);
printf("Pendiente:   geo.beta    %f [grados]\n",geo.beta);
printf("Altitud:     geo.altitud %f [km]\n",geo.altitud);
printf("Clima        geo.clima   %d [.]\n",geo.clima);
printf("Albedo       geo.albedo  %f [.]\n",geo.albedo);
printf("Amanecer    geohr.omegas %f [hr]\n",tiempo_hr(geohr.omegas)-12);
}

void imprimegeohr(sgeohr geohr)
{
printf("-------------------------------\n");
printf("geohr.hr     %f\n",geohr.hr);
printf("geohr.omega  %f\n",geohr.omega);
printf("geohr.omegas  %f\n",geohr.omegas);
printf("geohr.theta  %f\n",geohr.theta);
printf("geohr.alphas %f\n",geohr.alphas);
printf("geohr.gammas %f\n",geohr.gammas);
}
void fimprimegeohr(sgeohr geohr)
{
fprintf(stderr,"-------------------------------\n");
fprintf(stderr,"geohr.hr     %f\n",geohr.hr);
fprintf(stderr,"geohr.omega  %f\n",geohr.omega);
fprintf(stderr,"geohr.omegas  %f\n",geohr.omegas);
fprintf(stderr,"geohr.theta  %f\n",geohr.theta);
fprintf(stderr,"geohr.alphas %f\n",geohr.alphas);
fprintf(stderr,"geohr.gammas %f\n",geohr.gammas);
}


void imprimecubierta(scubierta cubierta)
{
printf("-------------------------------\n");
 printf("cubierta.KL    %f\n",cubierta.KL);
 printf("cubierta.alpha %f\n",cubierta.alpha);
 printf("cubierta.N     %d\n",cubierta.N);
 printf("cubierta.n2    %f\n",cubierta.n2);
}
