TARGET = rtubo

$(TARGET): R-tubo.o libSE.a
	gcc $^ -o $@ -lm

R-corto.o: R-tubo.c  
	gcc -c $< -o $@ -lm

libSE.a: lib1.o lib2.o lib3.o lib4.o lib5.o lib6.o lib7.o lib8.o lib9.o
	ar rcs $@ $^ 

lib1.o: uCorrelaciones.c uCorrelaciones.h
	gcc -c -o $@ $< -lm

lib2.o: uPropFisicas.c uPropFisicas.h
	gcc -c -o $@ $< -lm

lib3.o: uEstimacionRad.c uEstimacionRad.h param.h uEstructuras.h 
	gcc -c -o $@ $< -lm


lib4.o: uColectores.c uColectores.h param.h uEstructuras.h
	gcc -c -o $@ $< -lm

lib5.o: uCubiertas.c  uCubiertas.h param.h uEstructuras.h
	gcc -c -o $@ $< -lm


lib6.o: uImprime.c uImprime.h param.h uEstructuras.h
	gcc -c -o $@ $< -lm

lib7.o: uSSS.c uSSS.h param.h uEstructuras.h
	gcc -c -o $@ $< -lm

lib8.o: uTanques.c  uTanques.h param.h uEstructuras.h
	gcc -c -o $@ $< -lm

lib9.o: uPropiedades.c  uPropiedades.h param.h uEstructuras.h uPropFisicas.h
	gcc -c -o $@ $< -lm

clean:
	rm -f *.o *.a $(TARGET)