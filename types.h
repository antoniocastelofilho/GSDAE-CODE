/* ******************************************************* */
/*                                                         */
/*                    CODIGO GSDAE                         */
/*                   ARQUIVO types.h                       */          
/*                                                         */
/* INDICADO PARA INTEGRAR EQUACOES ALGEBRICO-DIFERENCIAIS  */
/*  DE QUALQUER ORDEM E INDICE ZERO OU UM PODENDO CONTER   */
/*                   SINGULARIDADES                        */
/*                                                         */
/* ******************************************************* */
/*                                                         */
/* VERSAO 1 - implementada por Renato Cesar da Silva       */
/*            em 07/96                                     */
/*                                                         */
/* VERSAO 2 - implementada por Antonio Castelo Filho       */
/*            em 12/96                                     */
/*                                                         */
/* VERSAO 3 - atualizada   por Antonio Castelo Filho      */
/*            em 01/97                                    */
/*                                                        */
/* ******************************************************* */



/* ******************************************************* */
/* Verificando se TYPES foi definido                       */
/* ******************************************************* */
#ifndef TYPES


/* ******************************************************* */
/* Definindo os tipos de dados usados em GSDAE             */
/* ******************************************************* */

/* definindo um numero real como double  */
typedef double   real;    

/* definindo um vetor de numeros inteiros  */
typedef int     *vint;

/* definindo uma matriz de numeros inteiros */
typedef int    **mint;

/* definindo uma matriz tridimensional de numeros inteiros */
typedef int  ***mmint;   

/* definindo um vetor de numeros reais */
typedef real    *vreal;

/* definindo uma matriz de numeros reais */
typedef real   **mreal;

/* definindo uma matriz tridimensional de numeros reais  */
typedef real  ***mmreal;   


/* ***************************************************** */
/* definindo a estrutura parameter que armazenara todos  */
/* os dados utilizados em GSDAE                          */
/* ***************************************************** */

typedef struct parameter  parameter; 

struct parameter {
  /* variavel de comprimento de arco */
  real   s;
  /* variavel x */
  real   cx;
  /* variavel y */
  mreal  cy;
  /* derivada de x */
  real   pdcx; 
  /* derivada de y */
  mreal  pdcy;
  /* jacobina de F */
  vreal  DFx;
  mmreal DFy;
  /* dimensao do vetor y */
  int    n;     
  /* ordem da EAD */
  int    o;
  /* rotina que define a EAD */
  void   (*F)(int,int,real,mreal,vreal);
  /* rotina que define a jacobiana */
  void   (*DF)(int,int,real,mreal,vreal,mmreal);
  /* ordem do metodo */
  int    k;
  /* posto de DFy[o] */
  int    rank;
  /* tolerancias */
  real   atolx;
  mreal  atoly;
  real   rtolx;
  mreal  rtoly;
  vreal  ftol;
  /* tamanho do passo */
  real   h;
  real   h0;
  real   hmin;
  real   hmax;
  real   cdmax;
  /* variaveis para armazenar os valores anteriores */
  real   sold;
  real   hold;
  int    kold;
  real   cjold;
  /* coeficientes do metodo */
  int    ns;
  real   cj;
  real   factor;
  real   alfas;
  real   alfa0;
  real   ck;
  real   S;
  real   SP;
  vreal  psi;
  vreal  alfa;
  vreal  beta;
  vreal  gama;
  vreal  sigma;
  vreal  phix;
  mmreal phiy;
  /* controladores */
  int    ifase;
  real   dir;
  /* variaveis para armazenar a funcao avaliada */
  vreal  u;
  vreal  delta;
  vreal  deltax;
  vreal  deltah;
  vreal  deltahx;
  /* matrizes para decomposicao QR */
  mreal  DH;
  mreal  Q;
  /* permutacoes */
  vint   p;
  vint   q;
  vint   paux;
  vint   qaux;
  /* contadores */
  int    nstep;       
  int    nNwf;       
  int    naF;        
  int    naDH;       
  int    ndQR;       
  int    nff;
  int    nstart;
  int    nDH;         
  int    aDH;
  int    cfalhas;
  int    nflhs;      
  int    irank;
  /* variaveis auxiliares */
  real   xend;
  real   x;
  mreal  y;
  real   pcx;
  mreal  pcy;
  real   cxx;
  mreal  cyx;
  real   pdcxold; 
  real   dx;
  mreal  dy;
  real   ccx;
  mreal  ccy;
  mreal  yx;
  real   wtx;
  mreal  wty;
  real   taux;
  mreal  tauy;
  real   Ex;
  mreal  Ey;
};


/* ******************************************************* */
/* Definindo as macro-funcoes utilizadas em GSDAE          */
/* ******************************************************* */

/* definindo a funcao maximo entre dois valores */
#define MAX2(x,y) ( ((x) > (y)) ? (x) : (y) ) 

/* definindo a funcao minimo entre dois valores */
#define MIN2(x,y) ( ((x) > (y)) ? (y) : (x) )

/* definindo a funcao maximo entre tres valores */
#define MAX3(x,y,z) ( MAX2(MAX2(x,y),z) ) 

/* definindo a funcao minimo entre tres valores */
#define MIN3(x,y,z) ( MIN2(MIN2(x,y),z) )

/* definindo a funcao sinal de um inteiro  */
#define SIGN(x) ( ((x) >= (0)) ? (1) : (-1) )

/* definindo a funcao sinal de um real  */
#define FSIGN(x) ( ((x) >= (0.0)) ? (1.0) : (-1.0) )

/* definindo a funcao quadrado */
#define SQR(x) ( (x) * (x) )

/* definindo TYPES */
#define TYPES

/* fim do if */
#endif
