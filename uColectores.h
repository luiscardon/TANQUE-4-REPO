

double Ut_agarwal81(int N, double Tp, double Ta, double epsilonp, double epsilong, double hw, double beta );
double Ut_klein80(int N, double Tp, double Ta, double epsilonp, double epsilong, double hw,double beta );
float fpp(float ccr);
float Qu(float Ac,float Fr, float S, float UL, float Ti, float T_a);
float ooQuW(struct s_colector_g co, float S,float Ti, float Ta);
float ooQuJ(struct s_colector_g co, float S,float Ti, float Ta);

//float Fp(float UL, float W, float D,float Cb,float Di,float hfi); // ojo
float Fr(float Fp, float Fpp);
float S(float Rb, float Ib, float Id, float Ig,
        float ptaoalfab  ,float ptaoalfad  ,float ptaoalfag, float  rhod, float beta);
float SS(float Ib, float Id, float Ig, float ptaoalfab  ,
	float ptaoalfad  ,float ptaoalfag);
float SSdis(float * Ib, float * Id, float * Ig, float ptaoalfab  ,
	float ptaoalfad  ,float ptaoalfag);

// colectores
float Faleta(float U_L, float k, float delta, float W, float D);
float Faleta_mm(float m);
float Faleta_kdelta(float U_L,float kd, float W, float D);
float mm_aleta(float U_L,float k, float delta,float W, float D);
float Fprima(float F, float U_L,float W, float D, float D_i, float C_b, float h_fi);
//float Col_cap(float mpunto, float colC, float colA, float U_L, float Fp);
float T_fo(float Ta, float Tfi, float S, float U_L, float colCap );
float F_R(float Fp,float Fpp);
float Col_cap(float mpunto, float Cp_fluido, float colA, float U_L, float Fp);
float Q_u(float colA, float FR, float S, float U_L, float Ti, float Ta);
float F_pp(float colCap);
//colectores tipo
float Fp_colector_tipo_a(float W, float D, float h, float C_b ,float U_L, float F);
float Fp_colector_tipo_b(float W, float D, float h, float C_b ,float U_L, float F);
float Fp_colector_tipo_c(float W, float D, float h, float C_b ,float U_L, float F);
float litrosampunto(float litrosmin);
float maleta(float UL,float k,float delta);

// colectores.h
float Tcolector(float Tcfi, float Ta, float S, struct s_colector_g micol);
// colectores-aux.h
// colector-aux.h
void imprime_colector(struct s_colector_g a);
void imprime_colector2(struct s_colector_g a,struct s_placa b);
FILE * imprime_colector_log(FILE * log, struct s_colector_g a);

