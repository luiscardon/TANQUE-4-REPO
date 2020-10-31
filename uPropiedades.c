#include "uEstructuras.h"
#include "uPropiedades.h"
#include "uPropFisicas.h"

// uPropiedades.c
// antes propiedades.c
// Incorporada en TANQUES-4 14 de abril de 2016

// cambiar a uMateriales
// en la struct materiales (propiedades) incluir el nombre

// las correlaciones ponerlas en otrolado


float propiedades_aire(float T, pro aire )
{
  aire.rho=1.;
  aire.cp=1.;
  aire.nu=1.;
  aire.mu=1.;
  aire.alpha=1.;
  aire.k=1.;
  aire.Pr=0.7;3.
}
float propiedades_agua(float T, pro agua )
{
  // hacer 
  agua.rho=1000.;
  agua.cp=4190;
  agua.nu=1.;
  agua.mu=1.;
  agua.alpha=1.;
  agua.k=0.63; // a 40C Bejan
  agua.Pr=10;
}
float print_propiedades(pro material )
{
  // imprime la propiedades termodinamicas de un material
  
fprintf(stderr,"#Propiedades de material:\n\n");
fprintf(stderr,"#Densidad              %f\n",material.rho);  
fprintf(stderr,"#Calor especifico      %f\n",material.cp);  
fprintf(stderr,"#Viscosidad cinematica %f\n",material.nu);  
fprintf(stderr,"#Viscosidad  dinámica  %f\n",material.alpha);  
fprintf(stderr,"#Difusividad térmica   %f\n",material.k);  
fprintf(stderr,"#Numero de Prandlt     %f\n",material.Pr);  
fprintf(stderr,"#------------------------\n");  

}

float Cp_agua(float T)
{ 
// calor específico del agua, J/kgK
// Qiu, pp. 21
return 4209.1 - 132.8e-2*T + 143.2e-4*T+T;
}
float k_agua(float T)
{
// conductividad térmica del agua, W/mK
return 0.5762 + 10.5e-4*T;
}
float Pr_agua(float T)
{
  // numero de Prandlt
  return 39.5345*pow(T,-0.144)-18.8396;
}
float beta_agua(float T)
{
  // coeficiente de expancion termica, K-1
  
  return (0.8*pow(T,0.5348) -1.9114)*1.e-4;
}
float nu_agua(float T)
{
  return 1.477e-6 * exp(-1.747e-2*T);
}
float rho_agua(float T)
{
  return 1000.;
}

void print_propiedades_agua(float T)
{
  printf("Cp %f \n ",Cp_agua(T) );
  printf("Pr %f \n ",Pr_agua(T) );
  printf("k %f \n ",k_agua(T) );
  printf("beta %f \n ",beta_agua(T) );
  printf("nu %f \n ",nu_agua(T) );
  printf("rho %f \n ",rho_agua(T) );
  printf("===============================\n ");  
}

void p_Hormigon(pro material )
{

  mateial.rho=1200.; // kg/m3
  material.cp=837.;   // J/kg C
  material.k=1.4;     // W/mC

  // no significativos
  material.nu=1.;
  material.mu=1.;
  material.alpha=1.;
  material.Pr=1.;

  
}

void p_BaldosaCeramica(pro material )
{
  mateial.rho=2000.; // kg/m3
  material.cp=800.;   // J/kg C
  material.k=1.;     // W/mC

  // no significativos
  material.nu=1.;
  material.mu=1.;
  material.alpha=1.;
  material.Pr=1.;

}


void p_BaldosaPorcelana(pro material )
{
  mateial.rho=2300.; // kg/m3
  material.cp=840.;   // J/kg C
  material.k=1.3;     // W/mC

  // no significativos
  material.nu=1.;
  material.mu=1.;
  material.alpha=1.;
  material.Pr=1.;

}



void p_BaldosaGres(pro material )
{
  mateial.rho=2500.; // kg/m3
  material.cp=1000.;   // J/kg C
  material.k=2.3;     // W/mC

  // no significativos
  material.nu=1.;
  material.mu=1.;
  material.alpha=1.;
  material.Pr=1.;

}



void p_MorteroCalCemento(pro material )
{
  // Mortero de caly cemento
  mateial.rho=1900.; // kg/m3
  material.cp=1000.;   // J/kg C
  material.k=0.7;     // W/mC

  // no significativos
  material.nu=1.;
  material.mu=1.;
  material.alpha=1.;
  material.Pr=1.;

}
void p_Hormigon(pro material )
void p_BaldosaCeramica(pro material )
void p_BaldosaPorcelana(pro material )
void p_BaldosaGres(pro material )
void p_MorteroCalCemento(pro material )
void p_MorteroRevoque(pro material )
void p_PoliestirenoExp(pro material )



void p_MorteroRevoque(pro material )
{
  // Mortero de revoque
  mateial.rho=1800.; // kg/m3
  material.cp=0.;   // J/kg C
  material.k=1.163;     // W/mC

  // no significativos
  material.nu=1.;
  material.mu=1.;
  material.alpha=1.;
  material.Pr=1.;

}
void p_PoliestirenoExp(pro material )
{
  // Poliestireno expandido
  // a densidad 1000
  mateial.rho=16.; // kg/m3
  material.cp=837.;   // J/kg C
  material.k=0.034;     // W/mC

  // no significativos
  material.nu=1.;
  material.mu=1.;
  material.alpha=1.;
  material.Pr=1.;

}



