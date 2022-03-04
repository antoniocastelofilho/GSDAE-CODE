DESTDIR=
%OBJS1= gsdae.o eqcamp.o
%OBJS2= gsdae.o exexp.o
%OBJS3= gsdae.o exesf.o
OBJS4= gsdae.o exvdp.o
%OBJS5= gsdae.o exedo.o
BINS=exvdp
%BINS= eqcamp exexp exesf exvdp exedo
#CC= gcc
#CFLAGS= -Wall -O3
CC= cc
CFLAGS= -g
LIBS = -lm 
%eqcamp: ${OBJS1}
	${CC} ${CFLAGS} ${LDFLAGS} -o eqcamp ${OBJS1} ${LIBS}

%exexp: ${OBJS2}
	${CC} ${CFLAGS} ${LDFLAGS} -o exexp ${OBJS2} ${LIBS}

%exesf: ${OBJS3}
	${CC} ${CFLAGS} ${LDFLAGS} -o exesf ${OBJS3} ${LIBS}

exvdp: ${OBJS4}
	${CC} ${CFLAGS} ${LDFLAGS} -o exvdp ${OBJS4} ${LIBS}

%exedo: ${OBJS5}
	${CC} ${CFLAGS} ${LDFLAGS} -o exedo ${OBJS5} ${LIBS}
