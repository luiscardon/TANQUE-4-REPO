#include "param.h"
#include "uEstructuras.h"
#include "uCorrelaciones.h"
#include "uPropFisicas.h"
#include "uTanques.h"



// =========================================================
// Funciones relacionadas al modelo del tanque estratificado
// y funciones de impresion
// Modelo de tanque estratificado de Kleinbach

// =========================================================
float FC(int i, int N, float Tco, float Ts[])
{  
  // Funcion de control de colector
  // DB, pp384
  // N = nuumero de estratos
  // notacion de DB  notacion de C
  // i=1..N          i=0 N-1
  
  // Tco temperatura del fluido de ingreso proveniente del colector
  float fc;
  
  if (i==0 && Tco > Ts[i]) // si i==0 Ts[i]= Ts[0]= Ts0
     { //printf("i=0=%d   y   Tco=%f  >  Ts[i]=%f\n",i, Tc0,Ts[i]);
       fc=1;
       }
  else if (Ts[i-1] >= Tco && Tco > Ts[i])
     { //printf("i != 0 = %d  y Ts[i-1]=%f  >= Tc0=%f  && Tc0 > Ts[i]= %f ",i, Ts[i-1],  Tc0,  Ts[i]);
       fc= 1;
       }
  else if (i==0-1 || i==N)   // fuera de límites
       fc=0.;
  else
       fc=0;

       return fc;
}
float FL(int i, int N, float TLr, float Ts[])
{
  // Funcion de control de retorno  
  // DB, pp384
  // TLr temperatura del fluido de retorno de la carga

  // N = nuumero de estratos
  // notacion de DB  notacion de C
  // i=1..N          i=0 N-1

  float fl;
  //printf("funcion FL i = %d\n",i);
  if(i==N-1  && TLr <= Ts[N-1] )
  {// printf("i = N-1 = %d  y TL%f  >=  Ts[i]= %f \n",i, TL, Ts[i] );
    fl= 1;}
  else if (Ts[i] >= TLr && TLr > Ts[i+1] )
  { // printf("i=0=%d Ts[i-1]=%f >= TL=%f  && TL > Ts[i]=%f \n",i,Ts[i-1],TL,Ts[i]);
    fl=1;}
  else if (i==-1 || i==N)  // fuera de limites
  {fl=0.;}
  else
  {  //printf("9999 funcion FL i = %d\n",i);
     fl=0.;}

    return fl;
}

float mpi(int i, int N, float mpC ,float fC[], float mpL, float fL[])
{
// calculo del flujo neto entre los tanques
// mixed flow rate
// mpi flujo neto del nodo i del nodo i-1
// no considera el flujo directo proveniente de la carga

// DBv pp. 335
// notese que las sumatorias dan 1 o 0, ya que
// solo uno de sus elemntos puede ser 1

  // N = numero de estratos
  // notacion de DB  notacion de C
  // i=1..N          i=0 N-1

int j;
float sumfc ,sumfl,mp;
sumfc=0;
sumfl=0;

//printf("%d %d\n",i,N);
//if(i==0)// || i==N) 
//{  mp=0;
//  printf("i==0  i=%d %d => mp=0\n",i,N);
//}
if (i >=0 && i <= N)  
{
  //printf("i >0 && i < N-1 i=%d  %d\n",i,N-1);
//  printf("calculo de Fc y FL\n");

for(j=0;j<=i-1;j++)
{//printf("j=%d ,  fC= %f\n",j,fC[j]);
  sumfc=sumfc + fC[j];
}
for(j=i;j<=N;j++) // Por que no considera el i mismo? -1q
{ //printf("j=%d ,  fL= %f\n",j,fL[j]);
  sumfl=sumfl + fL[j];
}  
//printf("i >0 && i < N-1 %d %d\n",i,N-1);
mp=mpC*sumfc - mpL * sumfl;
//printf("sumfC=%f mpC=%f sumfL=%f mpl=%f mp=%f\n", sumfc, mpC,sumfl,mpL, mp);
}
else if(i==N)   
{ // no hay estrato N, no hay flujo desde el estrato N   
  mp=0.;
  //printf("i=%d  N=%d  , \n",i,N);
//  printf("no hay estrato N, no hay flujo desde el estrato N\n");
  }
  else
  { //   printf("else i= %d\n",i);
  mp=9999.;}
//printf("desde mpi i=%d  mp=%f\n",i,mp);
return mp; 
}

float derivada_i(int i, int N ,float Ta, float Tc0,float TL, float Ts[], 
		 float fc[], float fl[], float mp[], float mpC, float mpL, struct datos_tanque tanque)
{
// Ta,   temperatura ambiente exterior al tanque
// Tc0,  temperatura del fluido prveniente de colector
// TL,   temperatura de fluido de retorno ed la carga
// Ts[i] temperatura del nod i
// fc[i] funcion
// fl[i] función 
// mp[i] flujo masico neto entrante al nodo i

// mpC   flujo másico proveniente de colector
// mpL   flujo másico proveniente de la carga

float U,A,Cp,m;
float sum=0.;
float UACp;
U=tanque.U;
Cp=tanque.cp;
A=tanque.Ale;
m=tanque.m;
UACp=(U*A/Cp);
//printf("Desde drivada, UA/Cp %f kg K/ seg\n" ,UACp);
//printf("i %d Ta %f Tc0 %f mpC %f mpL %f fc[i] %f fl[i] %f Ts[i] %f\n", i, Ta,Tc0,mpC, mpL, fc[i],fl[i],Ts[i]);

  sum=(U*A/Cp)*(Ta-Ts[i]);
  sum=sum+fc[i]*mpC*(Tc0-Ts[i]);
  sum=sum+fl[i]*mpL*(TL-Ts[i]);
//printf("sum %f\n",sum);
//printf("\n");
  if(mp[i]>0.)
    sum=sum + mp[i]*(Ts[i-1]-Ts[i]);

  if (mp[i]<0.)
    sum=sum+ mp[i+1]*(Ts[i]-Ts[i+1]);
  
//  printf("sum %f m %f sum/m %f K/seg \n",sum, m, sum/m);

  return sum/m;
  
}
void inicia_tanque(int N ,float Ta, float Tc,float TL, float Ts[], 
		 float fc[], float fl[], float mp[], float mpC, float mpL, struct datos_tanque tanque)
{int i;
 for(i=0;i<N;i++)
 {fc[i]=FC(i,N,Tc, Ts);
  fl[i]=FL(i,N,TL, Ts);
 }
   //
    // ====================================
 // calculo de los flujos mezclados  
 //    printf("Calculo de flujo mezclado\n");
 for(i=0;i<N;i++)
 {mp[i]=mpi( i, N,mpC , fc, mpL, fl);
  //printf("mp[i] %d %f \n",i,mp[i]); 
 } 
    printf("inicia_tanque: Calculo de los factores de control\n");
    printf("i Tc     TL     Ts[i]   FC            FL      mp[i]\n" );
    for(i=0;i<N;i++)
    printf("%d %4.2f  %4.2f  %4.2f  %3.1f  %3.1f  %f\n",i, Tc, TL, Ts[i], fc[i], fl[i],mp[i] );
   printf("====================================\n");
 
// imprime_mpi(N,mp,fc,fl);
}
//
// Funciones de impresion relacionadas
//
void imprime_mpi(int N, float mp[],float fc[],float fl[])
{int i;
  for(i=0;i<N;i++)
 {  printf("imprime_mpi: %d %f  %f  %f\n",i,mp[i],fc[i],fl[i]); 
 } 
}

// ================================================
// Funciones relacionadas al tanque de acumulacion
// #include "tanques-definicion.h"
// ================================================
// ================================================
/*
void tanque_geom(struct datos_tanque_g a, struct datos_tanque* tanque )
{ 
  // define un tanque
  // toma los datos geometricos basicos (e,k,D,L,rho,cp)
  // calcula
  //           parametros geometricos derivados
  //           parametros de transferencia
  
 float Ate; // area tapa y fondo exterior
 float Ale; // area lateral exterior
 float A;   // area total
 float m;   // masa edl fluido
 float U;   // coeficiente global de perdidas
 float Nu;  // Número ed Nusselt exterior
 float h;   // coeficiente convectivo
 
 float e;   // espesor de la aislacion
 float k;   // conductividad de la aislacion
 float D;   // diametro del tanque
 float rho; // desnidad del fluiod de trabajo
 float cp;  // calor especifico del fluido de trabajo
 float V;   // volumen del tanque
 float L;   // altura del tanque
 float R;   // resistencia termica de la asilacion
 float Pr;  // numero de Prandlt del aire
 float Ra;  // numero de Rayleigh exterior 
 //
 float ka;  // conductividad del aire
 
 // datos geometricos basicos
  e=  a.e; 
  k=  a.k;
  D=  a.D;
  L=  a.L;
  rho=a.rho;  
  cp= a.cp;

  //  Datos derivados
  Ate=M_PI*(D/2.)*(D/2.);
  Ale=M_PI*D*L;
  A=Ate+Ale;
  V=Ate*L;
  m=V*rho/6.;
  Nu=Nu_cilindro_vert(Pr,Ra);
  ka=0.024;
  h=Nu*ka/L; 
  R= 1./(h*Ale) + log((D+e)/D)/(2.*M_PI*k*L);
  U=1./R;
  
  tanque->Ale=Ale;
  tanque->Ate=Ate;
  tanque->A=A;
  tanque->rho=rho;
  tanque->cp=cp;
  tanque->m=m;
  tanque->U=U;
  tanque->V=V;
  tanque->R=R;
  tanque->h=h;
  tanque->Nu=Nu;
} 
*/
void tanque_geom_N(int N,struct datos_tanque_g a, struct datos_tanque* tanque )
{
  // define un estrato
  // toma los datos geometricos basicos (e,k,D,L,rho,cp)
  // calcula
  //           parametros geometricos derivados
  //           parametros de transferencia

  // derivada de tanque_geom_

 // se paso la estructura datos_tanque * por referencia
 // sus valores se modifican, es mas rapido
  
 // N requerido para que funcione como estrato
 // N numero de estratos 
 
  // datos geometricos  basicos: de  datos_tanque_g
  
 float e;   // espesor de la aislacion
 float k;   // conductividad de la aislacion
 float D;   // Diametro del tanque
 float L;   // Altura del tanque
 float rho; // Densidad del fluido del trabajo
 float cp;  // Calor específico del fluido de trabajo
 
  // datos derivados
  
 float Ate; // area tapa y fondo exterior
 float Ale; // area lateral exterior de cada estrato 
 float A;   // area total
 float m;   // masa edl fluido
 float U;   // coeficiente global de perdidas
 float Nu;  // Número ed Nusselt
 float h;   // coeficiente convectivo

 //  
 float V;   // Volumen del tanque
 float R;   // Resistencia termica de la aislacion
 float Pr;  // Numero de Prandlt del fluido de trabajo
 float Ra;  // Numero de Rayleigh, para el calculo de h 
 //
 float ka;  // conductividad del aire
 
  e=  a.e; 
  k=  a.k;
  D=  a.D;
  L=  a.L;
  rho=a.rho;  
  cp= a.cp;

  Ate=M_PI*(D/2.)*(D/2.);
  Ale=M_PI*D*L/N;
  A=Ate+Ale;
  V=Ate*L/N;
  m=V*rho; 
  Nu=Nu_cilindro_vert(Pr,Ra); // funcion no implementada devuelve 5
  ka=0.024;
  h=Nu*ka/L;                   // calculo en funcion de la altura total del tanque Ojo K del aire!
  R= 1./(h) + (D/ka)*log((D+e)/D);
  U=1./R;
  
  tanque->Ale=Ale;
  tanque->Ate=Ate;
  tanque->Ate=A;
  tanque->rho=rho;
  tanque->cp=cp;
  tanque->m=m;
  tanque->U=U;
  tanque->V=V;
  tanque->R=R;
  tanque->h=h;
  tanque->Nu=Nu;

} 
void tanque_geom(struct datos_tanque_g a, struct datos_tanque* tanque )
{
  // define un tanque
  // toma los datos geometricos basicos (e,k,D,L,rho,cp)
  // calcula
  //           parametros geometricos derivados
  //           parametros de transferencia

  // derivada de tanque_geom_

 // se paso la estructura datos_tanque * por referencia
 // sus valores se modifican, es mas rapido
  
 // N requerido para que funcione como estrato
 // N numero de estratos 
 
  // datos geometricos  basicos: de  datos_tanque_g
  
 float e;   // espesor de la aislacion
 float k;   // conductividad de la aislacion
 float D;   // Diametro del tanque
 float L;   // Altura del tanque
 float rho; // Densidad del fluido del trabajo
 float cp;  // Calor específico del fluido de trabajo
 
  // datos derivados
  
 float Ate; // area tapa y fondo exterior
 float Ale; // area lateral exterior de cada estrato 
 float A;   // area total
 float m;   // masa edl fluido
 float U;   // coeficiente global de perdidas
 float Nu;  // Número ed Nusselt
 float h;   // coeficiente convectivo

 //  
 float V;   // Volumen del tanque
 float R;   // Resistencia termica de la aislacion
 float Pr;  // Numero de Prandlt del fluido de trabajo
 float Ra;  // Numero de Rayleigh, para el calculo de h 
 //
 float ka;  // conductividad del aire
 
  e=  a.e; 
  k=  a.k;
  D=  a.D;
  L=  a.L;
  rho=a.rho;  
  cp= a.cp;

  Ate=M_PI*(D/2.)*(D/2.);
  Ale=M_PI*D*L;
  A=Ate+Ale;
  V=Ate*L;
  m=V*rho; 
  Nu=Nu_cilindro_vert(Pr,Ra); // funcion no implementada devuelve 5
  ka=0.024;
  h=Nu*ka/L;                   // calculo en funcion de la altura total del tanque Ojo K del aire!
  R= 1./(h) + (D/ka)*log((D+e)/D);
  U=1./R;
  
  tanque->Ale=Ale;
  tanque->Ate=Ate;
  tanque->A=A;
  tanque->rho=rho;
  tanque->cp=cp;
  tanque->m=m;
  tanque->U=U;
  tanque->V=V;
  tanque->R=R;
  tanque->h=h;
  tanque->Nu=Nu;

} 

void print_datos_o_tanque(struct datos_tanque tanque)
{
 // imprime los datos operativos del tanque de acumulacion
  
fprintf(stderr,"#\n#Datos operativos del tanque\n\n");
fprintf(stderr,"#Area de tapa y fondo           %f\n",tanque.Ate); 
fprintf(stderr,"#Area lateral                   %f\n",tanque.Ale);
fprintf(stderr,"#Volumen del tanque             %f\n",tanque.V);
fprintf(stderr,"#Densidad del fluido            %f\n",tanque.rho);
fprintf(stderr,"#Calor esecifico                %f\n",tanque.cp);
fprintf(stderr,"#Masa total                     %f\n",tanque.m);
fprintf(stderr,"#Coeficiente global de perdidas %f\n",tanque.U);
fprintf(stderr,"#Resistencia de la aislacion    %f\n",tanque.R);
fprintf(stderr,"#Coeficiente convectivo         %f\n",tanque.h);
fprintf(stderr,"#Numero de Nusselt externo      %f\n",tanque.Nu);
fprintf(stderr,"#---------------------------------------------\n");  
}
void print_datos_g_tanque(struct datos_tanque_g g)
{
  // imprime los datos de entrada minimos del tanque de acumulacion 
  
fprintf(stderr,"#\n#Datos geometricos y propiedades del tanque:\n\n");
fprintf(stderr,"#Espesor de la aislacion                    %f\n",g.e);  
fprintf(stderr,"#Conductividad del material de la aislacion %f\n",g.k);  
fprintf(stderr,"#Diametro interior                          %f\n",g.D);  
fprintf(stderr,"#Altura                                     %f\n",g.L);  
fprintf(stderr,"#Densidad del fluido                        %f\n",g.rho);  
fprintf(stderr,"#Calor especifico del fluido                %f\n",g.cp);  
fprintf(stderr,"#---------------------------------------------\n");  
} 

 void imprime_tanque_op(struct datos_tanque_op tanque)
{
 // imprime los datos operativos del tanque de acumulacion
  
fprintf(stderr,"#\n#Datos operativos del tanque\n\n");
fprintf(stderr,"#Temperatura agua alimentacion caliente   %f\n",tanque.Tc); 
fprintf(stderr,"#Temperatura agua alimentación fria       %f\n",tanque.TL);
fprintf(stderr,"#Temperatura del aire                     %f\n",tanque.Ta);
fprintf(stderr,"#Flujo masico entrada caliente            %f\n",tanque.mpC);
fprintf(stderr,"#Velocidad del fluido                     %f\n",tanque.vel);
fprintf(stderr,"#Flujo masico entrada fria                %f\n",tanque.mpL);
fprintf(stderr,"#Numero de estratos del modelo            %d\n",tanque.N);
fprintf(stderr,"#---------------------------------------------\n");  
}

// ==========================================================
// Tanques auxiliares
// #include "tanques-auxiliares.h"

void copia(int N, float a[] ,float b[])
{int i;
  for(i=0;i<N;i++)
    b[i]=a[i];
}    
void asigna(int N, float a ,float b[])
{int i;
  for(i=0;i<N;i++)
    b[i]=a;
}
void asigna_e(int N, int n,  float a , float c,float b[])
{int i;
  for(i=0;i<n;i++)
    b[i]=a;
  for(i=n;i<N;i++)
    b[i]=c;
  
}



void imprime_g(int N, float t, float Ts[])
{
// imprime para gnuplot
  int i;
for(i=0;i<N;i++)
    fprintf(stderr,"%d %f \n",N-i-1, Ts[i]);
    fprintf(stderr,"\n");
    fprintf(stderr,"\n");
}
void imprime_gnp_c(int N, float t, float Ts[], int cada, int * contador)
{
//   imprime para gnuplot con contador segun instrucciones

if(*contador==cada)
{  imprime_g(N,t,Ts); 
*contador=0;}

return;
}
void imprime_gnp_color(int N, float t, float Ts[], int cada, int * contador, int * color)
{
//   imprime para gnuplot con contador segun instrucciones

if(*contador==cada)
{  imprime_color(N,t,Ts,*color); 
*contador=0;
  *color=*color+3
  ;
}
return;
}
void imprime_gnp_color_normalizado(int N, float L,float t, float Ts[], int cada, int * contador, int * color)
{
//   imprime para gnuplot con contador segun instrucciones
//   llama a imprime_color_normalizado(N,L,t,Ts,*color)

if(*contador==cada)
{  imprime_color_normalizado(N,L,t,Ts,*color); 
*contador=0;
  *color=*color+1
  ;
}
return;
}

void imprime_color(int N, float t, float Ts[], int color)
{
// imprime para gnuplot
  int i;
fprintf(stderr,"#  %f \n",t);
  for(i=0;i<N;i++)
    fprintf(stderr,"%d %f %d\n",(N-i-1), Ts[i],color);
    fprintf(stderr,"\n");
    fprintf(stderr,"\n");
}
void imprime_color_normalizado(int N,float L, float t, float Ts[], int color)
{
// imprime para gnuplot
// imprime el perfil de temperatura en el tanque en el instante t
int i;
fprintf(stderr,"# tiempo %f \n",t);
  for(i=0;i<N;i++)
    fprintf(stderr,"%f %f %d\n",Ts[i],(L/N)*(N-i-1), color);
    fprintf(stderr,"\n");
    fprintf(stderr,"\n");
}
void imprimedat(FILE * fp,int N,float L, float t, float Ts[], int color)
{
// imprime para gnuplot
// imprime el perfil de temperatura en el tanque en el instante t
int i;
fprintf(fp,"# tiempo %f \n",t);
  for(i=0;i<N;i++)
    fprintf(fp,"%f %f %d\n",Ts[i],(L/N)*(N-i-1), color);
    fprintf(fp,"\n");
    fprintf(fp,"\n");
}

float conversion_flujo(float litros )
{
 // convierte litros por minuto a m3 por segundo
 //1min=60seg
 //1m3=1000litros 
  return litros/(1000.*60.);
}

//=======================================================

void imprime_gnp_cn_fc(int N, float L,float t, float Ts[], int cada, int * contador, int * color)
{
//   imprime para gnuplot con contador segun instrucciones
//   llama a imprime_color_normalizado(N,L,t,Ts,*color)

if(*contador==cada)
{  imprime_color_normalizado_fc(N,L,t,Ts,*color); 
*contador=0;
  *color=*color+1
  ;
}
return;
}



void imprime_color_normalizado_fc(int N,float L, float t, float fc[], int color)
{
// imprime para gnuplot
// imprime el perfil de temperatura en el tanque en el instante t
int i;
fprintf(stderr,"# tiempo %f \n",t);
  for(i=0;i<N;i++)
    fprintf(stderr,"%d %f %d\n",i,  t*fc[i],  color);
    fprintf(stderr,"\n");
    fprintf(stderr,"\n");
}

void print_control(  struct s_controlador  g)
{
  // imprime los datos de entrada minimos del tanque de acumulacion 
  
fprintf(stderr,"#\n#Control termostatico del tanque:\n\n");
fprintf(stderr,"#Temperatura minima:                %f\n",g.Tmin);  
fprintf(stderr,"#Temperatura maxima:                %f\n",g.Tmax);  
fprintf(stderr,"#Potencia auxiliar:                 %f\n",g.q);  
return;
}

float stepTTanqueMezcladoW(struct datos_tanque tanque, 
      float Ts, float Qu, float L, float Ta, float dt)
{
  // incremento de  temperatura del tanque mezclado
  // en el paso de tiempo
  // Qu potencia util, W
  // L potencia de la carga, W
  // Ta temperatura ambiente
  // dt, paso de tiempo
  float U,A,m,Cp;
  U=tanque.U;   // ojo revisar si corresponde al area total
  A=tanque.A;
  m=tanque.m;
  Cp=4190.;
//  printf("============== %f %f %f %f %f\n", U,A,(Ts - Ta),m,Cp);
//  printf("============== %g\n",dt*( Qu  - L - U*A*(Ts - Ta) )/(m*Cp));
  return  dt*( Qu  - L - U*A*(Ts - Ta) )/(m*Cp);
  
}
float stepTTanqueMezcladoJ(struct datos_tanque tanque, 
      float Ts, float Qu, float L, float Ta)
{
  // paso de temperatura del tanque mezclado
  // 
  // Qu potencia util, J en una hora
  // L potencia de la carga, J en una hora
  // Ta temoeratura ambiente
  // dt, paso de tiempo
  float U,A,m,Cp;
  U=tanque.U;   // ojo revisar si corresponde al area total
  A=tanque.A;
  m=tanque.m;
  Cp=4190.;
  return  ( Qu  - L - 3600*U*A*(Ts - Ta) )/(m*Cp);
  
}
float cargaLJ(float mpuntoL, float Ts, float Tsum, float es)
{
  // mpuntoL : flujo de la carga
  // TL      : temperatura de la carga
  // Tsum    : temperatura de suministro
  // esquema : esquema horario de utilizacion, minutos

 return mpuntoL*4190.*(Ts-Tsum)*es;
}

// factor de intercambio
float FRp(struct s_colector_g micolector, struct datos_tanque_op tanqueo, struct s_intercambiador interc )
{
  // Luis Cardon - 28 abril 2017
 // Factor colector-intercambiador de calor
 // DB. pp 427
float mpuntoCpC,mpuntoCpS,mpuntoCpmin;
float Ac,FR,mpuntoC,UL,CpC,CpS,mpuntoS,epsilon;
float aux;

Ac = micolector.Ac;
FR = micolector.FR;
UL = micolector.UL;
mpuntoC = micolector.mpunto;
CpC= micolector.Cp;
mpuntoCpC=mpuntoC*CpC;

mpuntoS=tanqueo.mpC;
CpS=4190.;
mpuntoCpS=mpuntoS*CpS;
mpuntoCpmin=minimo(mpuntoCpC,mpuntoCpS);

epsilon=interc.epsilon;

aux=1.+ ((Ac*FR*UL)/(mpuntoCpC))*(mpuntoCpC/(epsilon*mpuntoCpmin) -1.);
aux=1./aux;
return FR*aux;

}
float minimo(float a,float b)
{
 if(a >=b)return a;
 
 return b;
}

float colectorSerie(int N, struct s_colector_g micolector)
{
  
  // DB pp. 435
  // Db recalcula FRtaoalpha y FR UL
  // no obstante, parece que solo hace falta multiplicar FR por
  // el factor de correccion.
  float Ac,FR,UL,Cp,mpunto;
  float K,f;

  Ac=micolector.Ac;
  FR=micolector.FR;
  UL=micolector.UL;
  
  mpunto=micolector.mpunto;
  Cp=micolector.Cp;

  K=Ac*FR*UL/(mpunto*Cp);
  f=(1.-pow((1.-K),N))/(N*K);
  
  return FR=FR*f;
  
}

float dToff(struct s_colector_g micolector)
{
  // calcula la diferencia de temperatura de apagado
  // DB. pp 432
    float Ac,FR,UL,Cp,dTON,mpunto;
  
  Ac=micolector.Ac;
  FR=micolector.FR;
  UL=micolector.UL;  
  mpunto=micolector.mpunto;
  Cp=micolector.Cp;
  dTON=micolector.dTON;
 
  return Ac*FR*UL*dTON/(mpunto*Cp);

  
}
