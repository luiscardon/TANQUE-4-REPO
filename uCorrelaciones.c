#include "uCorrelaciones.h"
#include <math.h>

// uCorrelaciones
// reeplaza a correlaciones-h.h y a correlaciones.h
// incorporada en TANQUES-4

// Libreria de correlaciones para la transferencia de calor
// por conveccion
// Luis Cardon 14 de abril 2014

// 5/nov/2020, Gnielinski

double Num_Colburn(double Re,double Pr)
{
  // Correlacion de Colburn
  // Incropera pp.445
  // Flujo turbulento completamente desarrollado
  // en un tubo suave.
  // Propiedades evaluadas a Tm
  
  
  return 0.023 * pow(Re,4./5.) *pow(Pr,1./3.); 
}

double Num_Zhukauskas(double Re,double Pr,double Prs,double indice)
{
  double n,m,C;
  // Tabla de C y m para cilindros
  double cilindros[5][2]={ {0.75,0.4}, // 1<Re<40
                        {0.51,05},  // 40<Re<1000
                        {0.26,0.6}, // 10^3<Re<2e5
                        {0.076,0.7}, // 2e5<Re<2e6
                        {0.076,0.7}}; // 2e5<Re<2e6

//double cilindos[][]={ 0.75,0.4,0.51,05,0.26,0.6, 0.076,0.7, 0.076,0.7};
if(Re<40)  
{C=cilindros[0][1];
m=cilindros[0][2];
}if(Re<1000)
{C=cilindros[1][1];
m=cilindros[1][2];
}else if(Re<2.e5)
{C=cilindros[2][1];
m=cilindros[2][2];
}else if(Re<2.e6)
{C=cilindros[3][1];
m=cilindros[3][2];
}else
{C=cilindros[4][1];
m=cilindros[4][2];
}
n=0.47;
if(Pr>10)n=0.36;

return C * pow(Re,m)* pow(Pr,n)*pow( Pr/Prs,1./4.);
}

double Num_Churchill_Bernstein(double Re,double Pr)
{
// Correlacion de Churchill_Bernstein
// Incropera pp.470

  double up,down,fac;
if (Re*Pr>0.2)
{  up= 0.62 * pow(Re,1./2.)*pow(Pr,1./3.);
  down= pow(1 + pow(0.4/Pr,2./3.) , 1./4.);
  fac= pow(1. + pow(Re/282000.,5./8.)  , 4/5);
  return 0.3 + up*fac/down;
}else
{  return 0;}

}

double Num_Sieder_Tate(double Gz, double T,double Ts)
{ // Kreith Bhom pp. 382
  // para liquidos 
  // Incropera 444,
  // no menciona que sea para liquidos
  double mu,mus;
  mu=mu_aire(T);
  mus=mu_aire(Ts);
  mu=1;
  mus=1;
  
 return 1.86 * pow(Gz,1./3.) * pow (mu/mus,0.14); 
}
double Num_Hausen(double Gz)
{
 return 3.66 +  0.0668*Gz/(1 + 0.04*pow(Gz,2./3.));   
}
double Num_Dittus_Boelter(double Re,double Pr)
{
return 0.023*pow(Re,4./5.)*pow(Pr,0.4);
}
double Num_Swearingen_MXEligot(double Re,double Pr)
{
  // Kreith Bhom pp. 382
  
}
double Nu_Hollands(double Ra, float Pr, float beta)
{
  // Correlacion de Hollandas, 1976
  // para cavidades rectangulares inclinadas
  // requiere libreria matematica
  // datos
  // Ra : numero de Rayleigh, 
  // Pr : numero de Prandlt
  // beta: inclinacion de la caja
  
  // Ra se calcula en base a la separacion
  // entre placas
  
  float a,b,c,racosbeta,obeta;
  obeta=1.8*beta;
  obeta=obeta*M_PI/180.;
  
  beta=beta*M_PI/180.;
  
  racosbeta=Ra*cos(beta);
  a= 1708.*pow(sin(obeta),1.6)/(racosbeta);
  b=1.- 1708./racosbeta; if(b<0)b=0.;
  c= pow(racosbeta/5830., 1./3.) -1; if(c <0) c=0.;  
  return 1. + 1.44*(1.- a)*b +c;
}
float Nu_cilindro_vert(float Pr, float Ra)
{
  // Calcula el número de Nusselt exterior para cilindro vertical
  // podria aproximar con una placa plana

  // HACER !!!!
  
  return 5.;
}
float Num_McAdams(float Ra)
{ 
  // Correlacion de McAdams para
  // superficie superior de placa caliente
  // superficie inferior de placa fria
  
  // Incroprera pp. 498
  
 //  Dato: Ra_L
 //  calcular L=As/P
 
  // Devuelve Nu_L

  
  
  float Nu;
 if(1.e4<Ra && Ra <1.e7 )
   { Nu=0.54*pow(Ra,0.25);}
 else if (1.e7<Ra && Ra <1.e11 )
   {Nu=0.15*pow(Ra,0.333);}
  
 return Nu;
}

float Nu_Gnielinski(float Re, float Pr)
{
  // 
  float f,Nu;
  
  f=pow(1.182*log(Re)-1.64, -2.);

  if (Re< 2300)
    Nu= 3.66;
  else // (Nu < 5.E6)
  { Nu = (f/8)*(Re-1000.)*Pr;
    Nu=Nu / (1.+ 12.7*sqrt(f/8)*(pow(Pr, 2./3)-1.) );
  }
return Nu;      
}
