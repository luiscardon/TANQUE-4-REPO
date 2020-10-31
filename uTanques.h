//#include "estructuras.h"

// tanque-klein.h

float FC(int i, int N, float Tco, float Ts[]);
float FL(int i, int N, float TLr, float Ts[]);
float mpi(int i, int N, float mpC ,float fC[], float mpL, float fL[]);
float derivada_i(int i, int N ,float Ta, float Tc0,float TL, float Ts[], 
	 float fc[], float fl[], float mp[], float mpC, float mpL,
	 struct datos_tanque tanque);
void inicia_tanque(int N ,float Ta, float Tc,float TL, float Ts[], 
	 float fc[], float fl[], float mp[], float mpC, float mpL,
	 struct datos_tanque tanque);
void imprime_mpi(int N, float mp[],float fc[],float fl[]);
void imprime_gnp_color(int N, float t, float Ts[], int cada,
		       int * contador, int * color);
void imprime_gnp_color(int N, float t, float Ts[], int cada,
		       int * contador, int * color);
void imprime_color(int N, float t, float Ts[], int color);

// tanques definicion

void tanque_geom_N(int N,struct datos_tanque_g a, struct datos_tanque* tanque);
void tanque_geom(struct datos_tanque_g a, struct datos_tanque* tanque );
void print_datos_o_tanque(struct datos_tanque tanque);
void print_datos_g_tanque(struct datos_tanque_g g);

// tanques auxiliares

// tanques-auxiliares.h

void copia(int N, float a[] ,float b[]);
void asigna(int N, float a ,float b[]);
void asigna_e(int N, int n,  float a , float c,float b[]);
void imprime_g(int N, float t, float Ts[]);
void imprime_gnp_c(int N, float t, float Ts[], int cada, int * contador);
void imprime_color(int N, float t, float Ts[], int color);
float conversion_flujo(float litros );
void imprime_color_normalizado(int N, float L,float t, float Ts[], int color);
void imprime_gnp_color_normalizado(int N, float L,float t, float Ts[], int cada,
				   int * contador, int * color);
void imprime_color_normalizado_fc(int N,float L, float t, float fc[], int color);
void imprime_gnp_cn_fc(int N, float L,float t, float Ts[], int cada, int * contador, int * color);
void imprimedat(FILE * fp,int N,float L, float t, float Ts[], int color);
void print_control(  struct s_controlador  g);

// =========================================
float colectorSerie(int N, struct s_colector_g micolector);
float dToff(struct s_colector_g micolector);
float FRp(struct s_colector_g colector, 
	  struct datos_tanque_op tanqueo, struct s_intercambiador interc );
float minimo(float a,float b);
float stepTTanqueMezcladoW(struct datos_tanque tanque, 
      float Ts, float Qu, float L, float Ta, float dt);
float stepTTanqueMezcladoJ(struct datos_tanque tanque, 
      float Ts, float Qu, float L, float Ta);
float cargaLJ(float mpuntoL, float Ts, float Tsum, float es);
