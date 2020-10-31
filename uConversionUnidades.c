// funciones utiles de conversion de unidades
#include "param.h"


float litrosmin_metros3seg(float litrosmin)
{
 // datos:
 // litrosmin: flujo litros/minuto
 // devuelve:
 // metris3/seg: m^3/seg
 
 
 return litrosmin/(1000.*60.);
}
float metros3_litros(float metros3)
{
 return metros3 * 1000.; 
}