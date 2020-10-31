
// uCorrelaciones
// reeplaza a correlaciones-h.h y a correlaciones.h
// incorporada en TANQUES-4

// Libreria de correlaciones para la transferencia de calor
// por conveccion
// Luis Cardon 14 de abril 2014

// Porque Num?
// Ojo no todas estan terminadas
double Num_Colburn(double Re,double Pr);
double Num_Zhukauskas(double Re,double Pr,double Prs,double indice);
double Num_Churchill_Bernstein(double Re,double Pr);
double Num_Dittus_Boelter(double Re,double Pr); // HACER
double Num_Sieder_Tate(double Gz, double Tm,double Ts); // OJO prop
double Num_Hausen(double Gz);

double Nu_Hollands(double Ra, float Pr, float beta);
float  Nu_cilindro_vert(float Pr, float Ra); //HACER
float  Num_McAdams(float Ra);

// propiedades deben canbiar de lugar

