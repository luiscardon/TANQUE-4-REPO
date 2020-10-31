struct u_rad
{
  float I;    //  radiacion global
  float Ib;   //  radiacion directa
  float Id;   //  radiacion difusa
  float Ig;   //  albedo
};

struct u_geo 
{ // propiedades del lugar  y la superficie
  // dependen del dia y del lugar, independientes de la hora
  char * nombre;
  int   n;       // numero de dia, adimensional
  float delta;   // declinacion, grados
  float phi;     // latitud, grados  
  float longitud;// longitud, grados 
  float altitud; // altitud, metros  
  float gamma;   // angulo acimutal de la superficie, grados
  float beta;    // pendiente, grados
  int   clima;   // tipo de clima, Hottel, 76, DB 2006, pp 69 0 a 3
  float albedo;
  float Hmm[12]; // radiacion diaria media mensual, J
  float Tmm[12]; // temperatura diaria media mensual
  float Tmmax[12];
  float Tmmin[12];
};
struct u_geohr
{
  // relaciones geometricas entre un plano con cualquier orientacion
  // particular relativa a la tierra, para cualquier tiempo y a la radiacion
  // solar incidente.
  // Dependen de la hora
  // Definiciones: DB, pp.12

  float hr;      // hora, 0-24
  float omega;   // angulo horario, grados
  float omegas;  // angulo horario de la salida del sol, grados  
  float theta;   // angulo cenital, grados
  float alphas;  // angulo de altitud solar, grados, complementario de theta, grados
  float gammas;  // angulo acimutal solar, grados

  
};

struct u_cubierta
{
int N;             // numero de cubiertas
float KL;          // indice KL de extincion, =0.037;
float alpha;       // absortancia, =0.93;
float n2;          // indice de refraccion de mat. de la cubierta,=1.526 
};

typedef struct u_rad srad;
typedef struct u_geo sgeo;
typedef struct u_geohr sgeohr;
typedef struct u_cubierta scubierta;

// ======================================================
// tanques

struct datos_tanque_op
{
// datos minimos 
 float Ts;   // Temperatura inicial (tanque mezclado)
 float Tc;   // Temperatura del agua caliente de entrada 
 float TL;   // Temperatura del agua fria  de entrada 
 float Ta;   // Temperatura ambiente
 float mpC;  // flujo masico de entrada del agua caliente 
 float mpL;  // flujo masico de entrada del agua fria
 float vel;  // velocidad del fluido
 int   N;    // numero de estratos del tanque
 };
struct datos_tanque_g //basicos
{
// datos minimos 
 float e;   // espesor de la aislacion
 float k;   // conductividad del material de la aislacion
 float D;   // diametro interior 
 float L;   // altura
 float rho; // densidad del fluido
 float cp;  // calor esecifico de fluido
 float q;   // potencia auxiliar, W
 float recuperacion;
 };
 
 struct datos_tanque
 {
 // datos derivados
 float Ate; // area tapa y ondo exterior
 float Ale; // area lateral exterior
 float A;   // area total
 float V  ; // volumen del tanque
 float rho; // densidad del fluido
 float cp;  // calor especifico de fluido
 float m;   // masa del fluido
 float U;   // coeficiente global de perdidas
 float UA;   // coeficiente global de perdidas
 float R;   // resistencia termica de la aislacion
 float h;   // coeficiente convectivo externo
 float Nu;  // Numero de Nusselt externo
   };

struct s_controlador
{
  float Tmin; //
  float Tmax;
  float q;
};

struct propiedades
{
  float rho;
  float cp;
  float nu;
  float mu;
  float alpha;
  float k;
  float Pr;
};

typedef struct propiedades pro;

// typedef
// ===============================================
// Estructuras colector
// estructuras-col.h

struct s_placa
{  // placa
   float delta;     // espesor 
   float W;         // ancho aleta
   float epsilon;   // Emisividad 
   float k;         // Conductividad
   // tubos
   float D;  // Diametro externo 
   float Di; // Diametro interno 
   float Cb; // Conductancia de juntura
   float m;  // sqrt (UL/(k delta))
   float hi; // coeficiente de transferencia interno  
};

struct s_cubierta
{ // cubierta simple
  int N;       // numero
  float ep;    // emisividad
  float eg;    // emisividad 
  float e;     // espesor
  float delta; // separacion
  float beta;  // inclinacion
  float Ta;    // T aire
  float Tpm;   // T media placa
  float hw;    // H viento
};

struct s_aislacion
{
  float e;
  float k;
};

struct s_colector_g
{ // 
  // parametros globales del colector  
  float Ac;
  float mpunto; 
  float Cp;
  float UL;  
  float FR;   // Factor de remocion
  float F;    // Eficiencia de aleta
  float Fp;   // Factor de eficiencia
  float Fpp;  // Factor de flujo
  float rcC;  // Razon de capacitancia
  float dTON;
  float dTOFF;
  float pta;  // producto transmitancia absortancia promedio
  
  
};
struct s_colector_d
{ // parametros de datalle del colector
  float W;
  float D;
  float delta;    // espesor de la placa
  float k;        // conductividad de la placa
  float epsilon;  // emisividad de la placa
  float l;        // espaciamiento entre tubos
  float hi;       // coeficiente de transferencia interno  
};

struct s_intercambiador
{
  float epsilon;  // efectividad  de intercambiador
  
};
