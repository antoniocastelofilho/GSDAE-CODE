/* ****************************************************** */
/*                                                        */
/*                    CODIGO GSDAE                        */
/*                   ARQUIVO gsdae.h                      */          
/*                                                        */
/* INDICADO PARA INTEGRAR EQUACOES ALGEBRICO-DIFERENCIAIS */
/*  DE QUALQUER ORDEM E INDICE ZERO OU UM PODENDO CONTER  */
/*                   SINGULARIDADES                       */
/*                                                        */
/* ****************************************************** */
/*                                                        */
/* VERSAO 1 - implementada por Renato Cesar da Silva      */
/*            em 07/96                                    */
/*                                                        */
/* VERSAO 2 - implementada por Antonio Castelo Filho      */
/*            em 12/96                                    */
/*                                                        */
/* VERSAO 3 - atualizada   por Antonio Castelo Filho      */
/*            em 01/97                                    */
/*                                                        */
/* ****************************************************** */


/* ****************************************************** */
/*                  incluindo bibliotecas                 */
/* ****************************************************** */
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

/* ****************************************************** */
/*   incluindo os tipos de dados e macro-funcoes          */
/* ****************************************************** */
#include "types.h"


/* ****************************************************** */
/*   declarando a variavel global que armazena todos os   */
/*   dados para a execucao da rotina GSDAE e definindo-a  */
/*   como nao alocada                                     */
/* ****************************************************** */
parameter *par = NULL;


/* ****************************************************** */
/*   declarando todas as rotinas em gsdae.c               */
/* ****************************************************** */

int 
GSDAE (
int    n,           
int    o,           
real   h,           
real   hmin,        
real   hmax,        
real   cdmax,       
real  *s,           
real   send,        
real  *x,           
mreal  y,           
real   atolx,       
mreal  atoly,       
real   rtolx,       
mreal  rtoly,       
vreal  ftol,        
vint   infoinput,   
vint   infooutput   
);

int 
CSDAE (
int    n,           
int    o,           
real   h,          
real   hmin,        
real   hmax,        
real   cdmax,       
real  *s,           
real  *x,           
real   xend,        
mreal  y,           
real   atolx,       
mreal  atoly,       
real   rtolx,       
mreal  rtoly,       
vreal  ftol,        
vint   infoinput,   
vint   infooutput   
);

void 
STATISTICS (
real *s,
int  *npas,
int  *nreject,
int  *nsuc,
int  *nfunc,
int  *njac,
int  *nqr,
int  *nstart,
int  *nfnew
);

void 
STATUS (
int   status,        
char *mens
);

void 
weightvector (
int    n,
int    o,
int    r,
real   cx, 
mreal  cy,
vint   q,
real   atolx,
real   rtolx,
mreal  atoly,
mreal  rtoly,
real  *wtx,
mreal  wty
);

real 
weightnorm (
int    n, 
int    o,
int    r, 
real   cx, 
mreal  cy,
vint   q,
real   wtx, 
mreal  wty
);

void 
SETH (
int   n,          
int   o,
int   r,
real  h,
real  dpx,
mreal dpy,
real  x,        
mreal y,       
vint  p,
vint  q,
vreal delta, 
vreal deltah,
void  (*F)(int,int,real,mreal,vreal)
);

void 
SETDH (
int     n,          
int     o,
int     r,
real    cj,
real    h,
real    px,
mreal   py,
real    dpx,
mreal   dpy,
vint    p,
vint    q,
vreal   DFx,
mmreal  DFy,
void    (*DF)(int,int,real,mreal,vreal,mmreal),
mreal   DH 
);

void 
SETDHAPPROX (
int     n,
int     o,
int     r,
int     dim,
real    h,
real    cj,
real    uround,
real    x,
real    dx,
real    wtx,
mreal   y,
mreal   dy,
mreal   wty,
vint    p,
vint    q,
vreal   delta,
vreal   deltaaux,
vreal   deltah,
vreal   deltahaux,
mreal   DH,
void   (*F)(int,int,real,mreal,vreal) 
);

void 
predictor (
int     n, 
int     o, 
int     r, 
int     k, 
vreal   gama, 
vreal   phix,
mmreal  phiy,
real   *pcx,
mreal   pcy,
real   *pdcx,
mreal   pdcy,
vint    q,
vreal   u
);

void  
coefficient ( 
int     n,
int     o,
int     r,
int    *k,
real    h,
vreal   alfa,
vreal   beta,
vreal   gama,
vreal   sigma,
vreal   psi,
real   *alfas,
real   *cj,
real   *cjold,
real   *factor,
int    *aDH,
real   *ck,
vreal   phix,
mmreal  phiy,
int    *ns,
real   *hold,
int    *kold 
);

void 
update (
int    n,
int    o, 
int    kold,
int    kp1,
int    kp2,
real   Ex, 
mreal  Ey,  
vint   q,
vreal  phix,
mmreal phiy
);
 
int 
masterstep (
int     n,
int     o,
int     r,
real   *h,
real   *s,
int    *k,
int    *nNwf,
int    *aDH,
int     nDH,
real   *x,
mreal   y,
real   *cx,
mreal   cy,
real   *cxx,
mreal   cyx,
real   *pcx,
mreal   pcy,
real   *pdcx,
mreal   pdcy,
vreal   delta,
vint    p,
vint    q,
mreal   Q,
mreal   DH,
vreal   alfa,
vreal   beta,
vreal   gama,
vreal   sigma,
vreal   psi,
real   *alfas,
real   *cj,
real   *cjold,
real   *factor,
real   *ck,
vreal   DFx,
mmreal  DFy,
vreal   u,
real   *ccx,
mreal   ccy,
mreal   yx,
vreal   deltax,
real   *SP,
real   *S ,
vreal   phix,
mmreal  phiy,
real    hmin,
real    cdmax,
real    atolx,
real    rtolx,
mreal   atoly,
mreal   rtoly,
vreal   ftol,
int    *nff,
real   *wtx,
mreal   wty,
real   *Ex,
mreal   Ey,
int    *ifase,
int    *ns,
real   *hold,
int    *kold,
void   (*F)(int,int,real,mreal,vreal),
void   (*DF)(int,int,real,mreal,vreal,mmreal),
real   *sold,
real   *taux,
mreal   tauy,
int    *naF,
int    *naDH,
int    *ndQR,
vreal  deltah,
vreal  deltahx
);

int 
controlstep ( 
int     o,
int     n,
int     r,
int    *kold,
int    *k,
real    ck,
real   *hold,
real   *h,
real   *sold,
real   *s,
real   *cx,
mreal   cy,
real   *x,
mreal   y,
real   *cxx,
mreal   cyx,
vint    p,
vint    q,
vreal   beta,
vreal   sigma,
vreal   phix,
mmreal  phiy,
vreal   psi,
int    *nflhs,
int    *cfalhas,
real    atolx,
real    rtolx,
mreal   atoly,
mreal   rtoly,
real    wtx,
mreal   wty,
real    Ex,
mreal   Ey,
int    *ifase,
real    hmin,
real    hmax,
int     ns
);

void  
interpolator (
int     n,
int     o,
int     r,
real    hint,
real   *xint,
mreal   yint,
real   *dxint,
mreal   dyint,
vint    p,
vint    q,
int     kold,
vreal   phix,
mmreal  phiy,
vreal   psi
);

void  
interpolators (
int     n,
int     o,
int     r,
real   *send,
real   *x,
mreal   y,
real   *dx,
mreal   dy,
vint    p,
vint    q,
int     kold,
vreal   phix,
mmreal  phiy,
vreal   psi,
real    s
);

void  
interpolatorx (
int     n,
int     o,
int     r,
real    s,
real    x,
real    pdcx,
int     kold,
vreal   phix,
mmreal  phiy,
vreal   psi,
vint    p,
vint    q,
real   *send,
real   *xend,
mreal   yend,
real   *dxend,
mreal   dyend,
real    tol
);

void  
interpolatorsing (
int     n,
int     o,
int     r,
real    sold,
real    s,
real    xold,
real    x,
real    pdcxold,
real    pdcx,
vint    p,
vint    q,
int     kold,
vreal   phix,
mmreal  phiy,
vreal   psi,
real    tol,
real   *sint,
real   *xint,
mreal   yint,
real   *dxint,
mreal   dyint
);

void 
SETB (
int     n,
int     o,
int     r,
mreal   y,
vint    p,
vint    q,
vreal   DFx,
mmreal  DFy,
mreal   B
);

int 
settau (
int     n,
int    *o,
int    *r,
int     irank,
real    cx,
mreal   cy,
vreal   delta,
vreal   deltaaux,
vreal   DFx,
mmreal  DFy,
real    cdmax,
real    dir,
mreal   Q,
mreal   B,
vint    p,
vint    q,
vint    paux,
vint    qaux,
vreal   ftol,
real    atolx,
real    rtolx,
real   *taux,
mreal   tauy,
int     nDH,
int    *naF,
int    *naDH,
int    *nQR,
void    (*F)(int,int,real,mreal,vreal),
void    (*DF)(int,int,real,mreal,vreal,mmreal)
);

int 
rankneighbourhood (
int     n,
int     o,
int     r,
int     dir,
int     ord,
real    cx,
mreal   cy,
vreal   delta,
vreal   deltaaux,
vreal   DFx,
mmreal  DFy,
mreal   Q,
mreal   B,
vint    p,
vint    q,
int     nDH,
int    *naDH,
int    *nQR,
void   (*F)(int,int,real,mreal,vreal),
void   (*DF)(int,int,real,mreal,vreal,mmreal)
);

void 
DFAPPROX (
int     n,
int     o,
int     dir,
real    x,
mreal   y,
vreal   delta,
vreal   deltaaux,
vreal   DFx,
mmreal  DFy,
void    (*F)(int,int,real,mreal,vreal) 
);

void 
firststep (
int     n,
int     o,
int     r,
real    cx,
mreal   cy,
vreal   phix,
mmreal  phiy,
vint    p,
vint    q,
real    taux,
mreal   tauy,
real   *cj,
real   *cjold,
real   *factor,
real   *hold,
real    h,
int    *k,
int    *kold,
int    *ns,
vreal   psi,
int    *ifase,
real   *pdcx
);

int 
functionnorm (
int   n,
vreal x,
vreal tol
);

real 
POWER (
real x,
int  y
);

real 
ROOT ( 
real x,
int  y
);

void    
PIVOT (
int   n,
int   k,
mreal A,
mreal Q
); 

void 
GIVENS ( 
int   m,
int   n,
int   i,
int   k,
mreal A,
mreal Q,
real  s1,
real  s2
);

void 
QR (
int    m,
int    n,
mreal  A,
mreal  Q,
int    pivot,
real  *cond
);

void 
NEWTON (
int    n,
mreal  Q,
mreal  A,
vreal  u,
vreal  delta,
real  ac
);

real    
PIVOT2 (
int   n,
int   k,
mreal A,
vint  p,
vint  q
); 

void 
GIVENS2 ( 
int   n,
int   i,
int   k,
mreal A,
mreal Q,
vint  p,
vint  q,
real  s1,
real  s2
);

int 
QR2 (
int    n,
mreal  A,
mreal  Q,
vint   p,
vint   q
);

void 
SOLVESYSTEM (
int    n,
mreal  A,
mreal  Q,
vint   p,
vint   q,
vreal  x,
vreal  y
);

void 
ALLOCPAR (
int  n,
int  o,
mreal *y,
mreal *atoly,
mreal *rtoly,
vreal *ftol,
vint  *infoinput,
vint  *infooutput,
void (*F)(int,int,real,mreal,vreal),
void (*DF)(int,int,real,mreal,vreal,mmreal)
);

void 
FREEPAR ( 
int    n,
int    o,
mreal *y,
mreal *atoly,
mreal *rtoly,
vreal *ftol,
vint  *infoinput,
vint  *infooutput
);

vint   
ALLOCVINT (
int n
);

vint   
REALLOCVINT ( 
int  m,
int  n,
vint v
);

vint  
FREEVINT ( 
int  n,
vint v
);

mint   
ALLOCMINT ( 
int m,
int n
);

mint  
FREEMINT ( 
int  m,
int  n,
mint v
);

mmint  
ALLOCMMINT ( 
int m,
int n,
int k
);

mmint 
FREEMMINT (
int   m,
int   n,
int   k,
mmint v
);

vreal  
ALLOCVREAL (
int n
);

vreal   
REALLOCVREAL (
int   m,
int   n,
vreal v
);

vreal 
FREEVREAL (
int   n,
vreal v
);

mreal  
ALLOCMREAL (
int m,
int n
);

mreal 
FREEMREAL ( 
int   m,
int   n,
mreal v
);

mmreal  
ALLOCMMREAL ( 
int m,
int n,
int k
);

mmreal 
FREEMMREAL (
int    m,
int    n,
int    k,
mmreal v
);
