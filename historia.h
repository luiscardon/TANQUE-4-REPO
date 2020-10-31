  // 14 abril 2016
  // TANQUE-4 estrato-9
  // unifico todos los tanques en uTanques
  // unifico todas las estructuras en uEstructuras
  // Creo uPropFisicas
  // Creo uCorrelaciones
  // creo makefile para generar una libreria estatica
  
  // estrato-6.c
  // comienzo - 18 setiembre 2013
  // de estrato-5.c
  // de estrato3.c 
  // Luis Cardon mayo 2013?
  // ejecutable ./es, funciona?
  // revision 10 setiembre 2012
  // estrato-1.c  Funciona 11 setiembre 2013
  // estrato-2.c  Funciona 11 setiembre 2013
  // funcion inicia_tanque()
  // estrato-3.c  funciones a sub-estrato3.c y sub-estrato3.h
  
  // compilar: gcc -o es3 estrato3.c sub-estrato3.c -lm
  
  // datos para gnuplot desde la salida error: ./es3 2> dat
  // ok, 11 setiembre 2013
  // estrato-4.c agrego imprime_tanque_op, agrego imprime:gnp_c con contador
  // 16 setiembre 2013:  los flujos no parecen mover la termoclina,
  //                     solo la desciende la temperatura
  // agrego imprime_gnp_color y imprime_color
  // plot "dat" using 2:1:3 with lines palette
  // plot "dat" using 2:1:(column (-2)) with lines  lc variable 
  // load "itd.txt"
  //? por que los dibujos en color no salen en color en el despues de gimp?

  // estrato-5.c
  // de estrato-4.
  // 16 setiembre 2013:  TANQUES-1 
  // compilar con make, correr: ./est 2> dat 1> datos
  // 17 de setiembre, fisica correcta ... 18 parece que no  
  // pruebas:
  // resultados antisimetricos en los  casos de carga y descarga 
  // resultados graficos en /A-MP/ER/tanques.tex
  // 18 setiembre
  // mp[0] incorrecto, carga6.dat, dibujos est5-n6-carga.pdf ok ?
  // primera corrida: carga ok, descarga con problemas
  // segunda corrida, problema en el nodo superior...??? corregidos if i==0{} else if
  // ok
  // falta mejorar la salida para 0 < y < 1
  
  // 18 setiembre . TANQUES-2/estrato-6.c incorpora colector


  // problema tanque estratificado
  // Duffie y Beckman 1980 , pp 331  
  // Modelo de Kleinbach de tanque estratificado
  // verificar los retoques en el calculo de mp_i

