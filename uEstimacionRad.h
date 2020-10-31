
// UTILITARIOS - uEstimacionRad.c

// Luis Cardon 20 de setiembre de 2013
// derivado de fpp.c

// por que uEstimacionRad necesita uEstimacionRad. ??
// aqui creo que esta el problema


//#include "param.h"
//#include "uEstructuras.h"
//#include "uEstimacionRad.h"

// ------------------------------------------------------------------
// Cálculo de la irradiancia  instantánea, W/m^2
float Gon(int n);
float Go(int n, float thetaz);
float ooGon(sgeo geo);

float I_0(int n, float phi, float delta, float omega1, float omega2);
float I_0omega12(sgeo geo, float omega1, float omega2);
float ooI0(sgeo geo, sgeohr geohr);

float H_0(int n, float phi, float delta, float omegas);
float r_t(float omega, float omegas);
float r_thr(sgeohr geohr);
float oor_t(sgeohr geohr);

float rtt(float omega, float omegas);
// Metodo de Hottel
float tao_b(float thetaz, float A, int clima);
float ootao_b(sgeo geo,sgeohr geohr);

float tao_d(float taob);
float Icb(float Ion, float taob, float thetaz);
float I_T(int n, float phi, float delta, float I, 
	  float beta, float rhog, float hr);
float I_Tdis(int n, float phi, float delta, float I, 
	     float beta, float rhog, float hr,
	     float * Ibb, float * Idd, float * Igg );
float I_total(int n, float phi, float delta, float I, 
	     float beta, float gamma, float rhog, float hr,
	    float * Ibb, float * Idd, float * Igg );
srad ooITiso(sgeo geo, sgeohr geohr, srad I);

float frac_d(float kt);
float frac_d_o(float kt);

float costheta(float phi, float delta, float omega, 
	       float gamma, float beta);
float costhetaz(float phi, float delta, float omega);

float oocosthetaz(sgeo geo,sgeohr geohr);

float gradosarad(float grados);
float radagrados(float rad);
float delta_n(int n);
float Rb(float theta, float thetaz);
float Rb_a(float phi, float delta, float beta,
	   float gamma, float omega1, float omega2);
float ooRbab(sgeo geo, sgeohr geohr);
float mH_dmHratio(float omegas, float KT);
float rd(float omegas, float omega);
float oord(sgeohr geohr);
float omega_s(float phi, float delta);
float omega_min(float min);
float omega_seg(float seg);
float omega_hr(float hr);
float anguloh_h(float hora);
float anguloh_min(float min);
float anguloh_seg(float seg);
float tiempo_seg(float angulo);
float tiempo_min(float angulo);
float tiempo_hr(float angulo);
float costheta2(float phi, float delta, float omega, 
	       float gamma, float beta);
//float O_theta_hr(struct s_geo geo, struct s_geohr geohr);
srad  oohottel(sgeo geo, sgeohr geohr);
srad  oohottelDB(sgeo geo, sgeohr geohr);
void ooprintrad(sgeohr geohr, srad rad);
srad ooSzero();
int ordinal(int dd, int mm, int yy);
srad ooITkt(sgeo geo, sgeohr geohr,int mes);
srad ooIkt(sgeo geo, sgeohr geohr,int mes);
float I_horizontal(int n, float phi, float delta, float I, 
	      float gamma,  float hr, float albedo,
	    float * Ib, float * Id, float * Ig);
