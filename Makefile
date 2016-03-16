
GSL_INC = -I/home/ngan/libraries/gsl-1.16/include
GSL_LIB = -L/home/ngan/libraries/gsl-1.16/lib

CC = gcc

OPTIONS = -Wall -O3

driver: driver.o scf_potential.o Phi_nl.o gegenbauer.o ind.o
	$(CC) $(OPTIONS) -o driver driver.o scf_potential.o Phi_nl.o gegenbauer.o ind.o -lgsl -lgslcblas $(GSL_LIB)

driver.o: driver.c scf.h 
	$(CC) $(OPTIONS) -c driver.c

ind.o: ind.c
	$(CC) $(OPTIONS) -c ind.c

scf_potential.o : scf_potential.c
	$(CC) $(OPTIONS) -c scf_potential.c $(GSL_INC)
    
Phi_nl.o: Phi_nl.c
	$(CC) $(OPTIONS) -c Phi_nl.c $(GSL_INC)
	
gegenbauer.o: gegenbauer.c
	$(CC) $(OPTIONS) -c gegenbauer.c
	
clean:
	rm -f driver scf scf_spherical *.o *~
	
	
