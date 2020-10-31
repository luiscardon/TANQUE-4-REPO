float rfrac_n(float theta1, float theta2);
float rfrac_p(float theta1, float theta2);
float rfrac(float rn , float rp);
float tau_n(float rfracn);
float tau_rfrac(float rfracn, float rfracp);
float tau_n_N(int N, float rflecp);
float tau_r_N(float rpN, float rnN);
float tau_r(float rp, float rn);
float tau_a(float L, float theta2, float K);
float tau_a2(float KL1, float KL2,float theta2);

float tau_a_n(float taua, float rn);
float rho_a_n(float taua, float rn, float taun);
float alpha_a_n(float taua, float rp);
float tau_a_r(float taun, float taup);
float rho_a_r(float rhon, float rhop);
float alpha_a_r(float alphan, float alphap);
float tau_s(float taun , float taup);
float tau_simplif(float taua, float  taur);
float alpha_simplif(float taua);
float rho_simplif(float taua,float tau);
float theta_g_eff(float beta);
float theta_d_eff(float beta);
float ptaualpha(float tau, float alpha, float rhod);
float ptaualpha_s(float tau, float alpha);
float ta(int N, float  KL, float theta1,  float alpha, float n2);

FILE * ta_log(FILE *log, int N, float  KL, float theta1,  float alpha, float n2, float * pta);

float raz_alphan(float theta);
float tau_theta(int N, float  KL, float theta1,  float alpha, float n2);
float taos(int N, float KL, float theta,float beta,float alphann, float n2,
	    float * ptab,float *ptabd,float * ptag);

void  TyR_nn(int N, float tao, float rho, float * Tn, float * Rn);
float Tmn(float Tm,float Tn, float Rm, float Rn);
float Rmn( float Tm,float Tn, float Rm, float Rn); // ojo reuiere T
float TyR_nn_2(int N, float n,float K, float L, float theta1, float *Tn, float *Rn);
float TyR_CC(int N, float nn,float Kn, float Ln, 
int M, float nm,float Km, float Lm,float theta1, float *T, float *R);
FILE * TyR_CC_log(FILE *log, int N, float nn,float Kn, float Ln, 
int M, float nm,float Km, float Lm,float theta1, float *T, float *R);



float tau_cs(int N, float theta1, float n);
float theta_2(float n1, float n2, float theta1);

float alpha_alphan(float theta, float alphan);
float promedio(float a, float b);
