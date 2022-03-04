#include "gsdae.h"

void main ( void );
void FESF ( int, int, real, mreal, vreal );
void DFESF ( int, int, real, mreal, vreal, mmreal );

void main ( )
{
  int        n;
  int        q;
  int        r;
  real       h;
  real       hmin;
  real       hmax;
  real       cdmax;
  real       x;
  int        m;
  real       ATOLX;
  real       RTOLX;
  mreal      ATOLY;
  mreal      RTOLY;
  vreal      FTOL;
  mreal      y;
  vint       info;
  char       msg[250];
  real       s,sout,soutp,ds;

  int        erro,i,j,k,t;
  real       len;
  int        npas,nreject,nsuc,nfnew,nfunc,njac,nqr;

  n = 1;
  r = 1;
  q = 0;

  ALLOCPAR(n,q,r,&y,&ATOLY,&RTOLY,&FTOL,&info,FESF,DFESF);

  hmin   = 1.0e-16;
  hmax   = 0.0;
  cdmax  = 1.0e50;
  ATOLX = 1.0e-15;
  RTOLX = 1.0e-8;
  for (i = 1; i <= n; i++) {
    for (j = 0; j <= q; j++) {
      ATOLY[j][i] = 1.0e-15; 
      RTOLY[j][i] = 1.0e-8;
    }
    FTOL[i] = 1.0e-6;
  }
  info[1] = 0;
  info[2] = 1;
  h         = 1.0e-15;
  s         = 0.0;
  sout      = 6.5;
  x         = 0.0;
  y[0][1]   = 1.0;

  m = 10;

  ds = (sout-s)/(real) m;

  erro = 0;
  soutp  = s;

  printf("s = %14.16lf \n",s); 
  printf("x = %14.16lf \n",x); 
  for (i = 0; i <= q-1; i++)
    for (j = 1; j <= n; j++)
      printf("y[%d][%d] = %14.16lf \n",i,j,y[i][j]); 
  for (j = 1; j <= r; j++)
    printf("y[%d][%d] = %14.16lf \n",q,j,y[q][j]); 
  printf("\n\n"); 

  for (t = 1; (t <= m) && (erro >= 0); t++) {

    soutp += ds;

    do {

      erro = GSDAE (n,q,r,h,hmin,hmax,cdmax,&s,soutp,&x,y,
                    ATOLX,ATOLY,RTOLX,RTOLY,FTOL,info);

      printf("s = %14.16lf \n",s); 
      printf("x = %14.16lf \n",x); 
      for (i = 0; i <= q-1; i++)
        for (j = 1; j <= n; j++)
          printf("y[%d][%d] = %14.16lf \n",i,j,y[i][j]); 
      for (j = 1; j <= r; j++)
        printf("y[%d][%d] = %14.16lf \n",q,j,y[q][j]); 
      printf("\n\n"); 

    } while (erro == 1);

  }

  ERROR(erro,msg);
  STATISTIC(&len,&npas,&nreject,&nsuc,&nfunc,&njac,&nqr,&nfnew);

  printf("\n\nError menssage\n\n");
  printf("%s\n",msg);
  printf("\n\nStatistc\n\n");
  printf("Arc-lenght : %lf\n",len);
  printf("Number of Steps : %d\n",npas);
  printf("Number of Rejected Steps : %d\n",nreject);
  printf("Number of Sucess Steps : %d\n",nsuc);
  printf("Number of Newton Step Fail : %d\n",nfnew);
  printf("Number of Evaluation of the Function : %d\n",nfunc);
  printf("Number of Evaluation of the Jacobian : %d\n",njac);
  printf("Number of QR Decomposition : %d\n",nqr);

  FREEPAR(n,q,r,&y,&ATOLY,&RTOLY,&FTOL,&info);

  return;
}


void FESF( q, n, x, y, delta )
int    q;
int    n;
real   x;
mreal  y;
vreal delta ;
{
  delta[1] = x*x + y[0][1]*y[0][1] - 1.0;
}


void DFESF( q, n, x, y, DFx, DFy )
int    q;
int    n;
real   x;
mreal  y;
vreal  DFx;
mmreal DFy;
{
  DFx[1]       = 2.0*x;
  DFy[0][1][1] = 2.0*y[0][1];
}
