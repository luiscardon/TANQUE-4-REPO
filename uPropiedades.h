#include "uEstructuras.h"


// uPropiedades.c
// antes propiedades.c
// Incorporada en TANQUES-4 14 de abril de 2016

float propiedades_aire(float T, pro aire );
float propiedades_agua(float T, pro agua );

void p_Hormigon(pro material );
void p_BaldosaCeramica(pro material );
void p_BaldosaPorcelana(pro material );
void p_BaldosaGres(pro material );
void p_MorteroCalCemento(pro material );
void p_MorteroRevoque(pro material );
void p_PoliestirenoExp(pro material );


float print_propiedades(pro material );

float Cp_agua(float T);
float k_agua(float T);
float Pr_agua(float T);
float beta_agua(float T);
float nu_agua(float T);
float rho_agua(float T);  // constante

void  print_propiedades_agua(float T); // print_propiedades ?

