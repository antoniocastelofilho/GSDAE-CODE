/* ****************************************************** */
/*                                                        */
/*                    CODIGO GSDAE                        */
/*                   ARQUIVO gsdae.c                      */          
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
/* Incluindo o arquivo gsdae.h que contem                 */
/* os tipos de dados, macros e prototipo das              */
/* funcoes em gsdae.c                                     */
/* ****************************************************** */
#include "gsdae.h"



/* ****************************************************** */
/*                                                        */
/*               Rotinas GSDAE e CSDAE                    */
/*                                                        */
/*  Estas rotinas integram uma solucao geral de uma EAD   */
/*                                                        */
/*                F(x,y,y',..,y^(o)) = 0                  */
/*                                                        */
/*  onde y e suas derivadas sao representadas pelo vetor  */
/*  y[0],y[1],..,y[o], sendo "o" a ordem a EAD.           */
/*                                                        */
/*  A solucao geral representada por                      */
/*          c(s) = (x(s),y[0](s),..,y[o](s))              */
/*  e parametrizada por comprimento de arco, e pode       */
/*  conter singularidades transversais.                   */
/*                                                        */
/*  O sistema diferencial associado a esta EAD e dado     */
/*  por :                                                 */
/*                     |  F(c(s))         = 0             */
/*    H(c(s),c'(s)) = -|  w(c(s)) c'(s)   = 0             */
/*                     |  ||c'(s)||**2 -1 = 0             */
/*                                                        */
/*  O metodo numerico utilizado e o BDF com coeficientes  */
/*  principais fixos adaptados para o sistema diferencial */
/*  e para inicializar o metodo de Euler explicito para   */
/*  o problema equivalente :                              */
/*                                                        */
/*                  A(c(s))c'(s) = 0                      */
/*                                                        */
/*  onde                                                  */
/*             (   Fx   Fy[0] Fy[1] ... Fy[o-1] Fy[o] )   */
/*             (  y[o]    0    0    ...   -I      0   )   */
/*   A(c(s)) = ( y[o-1]   0    0    ...    0      0   )   */
/*             (   .      .    .     .     .      .   )   */
/*             (  y[1]   -I    0    ...    0      0   )   */
/*                                                        */
/*  As estrategias de controle de ordem e passo sao os    */
/*  as mesmas do codigo DASSL.                            */
/*                                                        */
/*  As rotinas GSDAE e CSDAE sao semelhantes quanto a     */
/*  entrada, saida de dados, metodos e estrategias. Elas  */
/*  diferem apenas quanto ao ponto final de integracao.   */
/*  A rotina GSDAE integra a solucao da condicao inicial  */
/*  c(s), ate um valor do parametro s definido por send   */
/*  ou uma singularidade ou erro.                         */
/*  A rotina CSDAE integra a solucao da condicao inicial  */
/*  c(s), ate um valor do parametro x definido por xend   */
/*  ou uma singularidade ou erro.                         */
/*                                                        */
/*                                                        */
/*  Descricao das rotinas :                               */
/*                                                        */
/*  Todos os parametros descritos a seguir sao validos    */
/*  para ambas as rotinas a menos que seja mensionado o   */
/*  contrario.                                            */
/*                                                        */
/*  Parametros das rotinas                                */
/*                                                        */
/*  n          : dimensao da EAD                          */
/*  o          : ordem da EAD                             */
/*  h          : passo de integracao                      */
/*  hmin       : passo minimo de integracao               */
/*  hmax       : passo maximo de integracao               */
/*  cdmax      : condicao maxima para a derivada          */
/*  s          : c(s) = (x(s),y(s))                       */
/*  x          : c(s) = (x(s),y(s))                       */
/*  y          : c(s) = (x(s),y(s))                       */
/*  send       : ponto final de integracao (rotina GSDAE) */
/*  xend       : ponto final de integracao (rotina CSDAE) */
/*  atolx      : erro absoluto para x                     */
/*  atoly      : erro absoluto para y                     */
/*  rtolx      : erro realtivo para x                     */
/*  rtoly      : erro realtivo para y                     */
/*  ftol       : erro absoluto para a funcao              */
/*  infoinput  : informacoes para entrada de dados        */
/*  infooutput : informacoes para saida de dados          */
/*                                                        */
/*                                                        */
/*  Parametro global                                      */
/*                                                        */
/*  par  : variavel que armazena todos os dados utiliza-  */
/*         dos nas rotinas                                */
/*                                                        */
/*         esta variavel e alocada na rotina ALLOCPAR e   */
/*         liberada na rotina FREEPAR                     */
/*                                                        */
/*                                                        */
/*  As rotinas GSDAE e CSDAE sao funces que retornam um   */
/*  valor inteiro. O valor retornado esta entre -16 e 4   */
/*  nesta versao. Este valor representa o estado da       */
/*  rotina no momento ou um erro detectado.               */
/*  O estado da rotina em condicoes de prosseguir sao     */
/*  representados por inteiros entre 0 e 5, e os erros    */
/*  detectados por inteiros negativos entre -16 e -1.     */
/*  Este valor representa:                                */
/*                                                        */
/*   0 : ponto regular o valor send (GSDAE) xend (CSDAE)  */
/*       foi atingido                                     */
/*                                                        */
/*   1 : singularidade transversal                        */
/*                                                        */
/*   2 : ponto regular e o posto de DFy[[o] diminuiu      */
/*                                                        */
/*   4 : ponto regular e a ordem da EAD diminuiu          */
/*                                                        */
/*   3 : singularidade transversal e o posto de DFy[o]    */
/*       diminuiu                                         */
/*                                                        */
/*   5 : singularidade transversal e a ordem da EAD       */
/*       diminuiu                                         */
/*                                                        */
/*  -1 : o parametro par nao foi alocado                  */
/*                                                        */
/*  -2 : erro na entrada de dados                         */
/*                                                        */
/*  -3 : ponto inicial inadequado - nao satisfaz a EAD    */
/*       com a tolerancia desejada                        */
/*                                                        */
/*  -4 : posto maior que o informado                      */
/*                                                        */
/*  -5 : EAD com ordem zero e posto de DFy[0] zero        */
/*                                                        */
/*  -6 : o posto de DFy[o] varia em uma vizinhaca do      */
/*       ponto                                            */
/*                                                        */
/*  -7 : a ordem diminuiu e o posto de DFy[o] varia em    */
/*       uma vizinhaca do ponto                           */
/*                                                        */
/*  -8 : singularidade nao tranversal                     */
/*                                                        */
/*  -9 : singularidade nao transversal com posto  de      */
/*       DFy[o] inferior                                  */
/*                                                        */
/* -10 : singularidade nao transversal com ordem inferior */
/*                                                        */
/* -11 : ponto inadequado - nao satisfaz a EAD com a      */
/*       tolerancia desejada                              */
/*                                                        */
/* -12 : h < hmin                                         */
/*                                                        */
/* -13 : condicao da jacobiana maior que a adimitida      */
/*                                                        */
/* -14 : numero de falhas na rotina masterstep maior que  */
/*       o admitida                                       */
/*                                                        */
/* -15 : houve falha na rotina e nehuma providencia foi   */
/*       tomada                                           */
/*                                                        */
/* -16 : ocorreu uma singularidade transversal na ultima  */
/*       chamada da rotina CSDAE e nehuma providencia foi */
/*       tomada                                           */
/*                                                        */
/*  Se a rotina GSDAE for chamada com o codigo de erro    */
/*  -15, o programa e abortado.                           */
/*  Se a rotina CSDAE for chamada com o codigo de erro    */
/*  -15 ou -16, o programa e abortado.                    */
/*                                                        */
/*  Parametros de entrada e saida                         */
/*                                                        */
/*  s          : variavel independente da parametrizacao  */
/*               da solucao geral c - valor real          */
/*                                                        */
/*  x          : variavel dependente x = x(s) que na EAD  */
/*               original e a variavel independente -     */
/*               valor real                               */
/*                                                        */
/*  y          : variavel dependente y = y(s) que na EAD  */
/*               original e a variavel dependente -       */
/*               matriz de numeros reais de dimensao o    */
/*               por n                                    */
/*                                                        */
/*               y[0] representa a variavel y, y[1]       */
/*               representa a primeira derivada de y,     */
/*               y[2] representa a segunda derivada de y, */
/*               y[o] representa a derivada de ordem o    */
/*               de y                                     */
/*                                                        */
/*                                                        */
/*  Parametros de entrada                                 */
/*                                                        */
/*  n          : dimensao da EAD - inteiro maior que zero */
/*                                                        */
/*  o          : ordem da EAD - inteiro maior ou igual a  */
/*               zero                                     */
/*               Se o = 0 a equacao e puramente algebrica */
/*                                                        */
/*  send       : valor final de integracao do parametro   */
/*               s - valor real (rotina GSDAE)            */
/*                                                        */
/*  xend       : valor final de integracao do parametro   */
/*               x - valor real (rotina CSDAE)            */
/*                                                        */
/*  h,hmin,hmax: passo atual, passo minimo e passo maximo */
/*               de integracao - valores reais maior ou   */
/*               igual a zero                             */
/*                                                        */
/*               se infoinput[1] = 0 (primeria chamada) e */
/*                                                        */
/*               se h = hmin = 0.0, sao definidos por     */
/*               default hmin = 1.0e-16, h = 1.0e-15      */
/*                                                        */
/*               se h = 0.0 e hmin != 0.0, h e definido   */
/*               por default h = 10.0*hmin                */
/*                                                        */
/*               se infoinput[1] = 1 (nao for a primeria  */
/*               chamada) e                               */
/*                                                        */
/*               se 0.0 < hmin < par->hmin, hmin e rede-  */
/*               finido pelo valor de entrada             */
/*                                                        */
/*               se 0.0 < h < par->h, h e redefinido pelo */
/*               valor de entrada como o minimo entre h   */
/*               par->h                                   */
/*                                                        */
/*  cdmax      : condicao maxima admitida para a matriz   */
/*               jacobiana - valor real maior ou igual a  */
/*               zero                                     */
/*                                                        */
/*               se cdmax = 0.0, cdmax e definido por     */
/*               default cdmax = 1.0e6                    */
/*                                                        */
/*  atolx,atoly: erro absoluto para c = (x,y) - atolx e   */
/*               um valor real maior ou igual a zero,     */
/*               atoly e uma matriz de dimensao o por n   */
/*               de valores reais maiores ou iguais a     */
/*               zero                                     */
/*                                                        */
/*  rtolx,rtoly: erro relativo para c = (x,y) - rtolx e   */
/*               um valor real maior ou igual a zero,     */
/*               rtoly e uma matriz de dimensao o por n   */
/*               de valores reais maiores ou iguais a     */
/*               zero                                     */
/*                                                        */
/*  ftol       : erro absoluto para as coordenadas da     */
/*               funcao que define a EAD - vetor de       */
/*               dimensao n de valores reais              */
/*                                                        */
/*               se ftol[1] = 0.0 o todo ponto e consi-   */
/*               derado satisfazendo a EAD no teste da    */
/*               rotina functionnorm                      */
/*                                                        */
/*               se infoinput[3] = 0 os valores de atolx, */
/*               atoly, rtolx, rtoly e ftol sao definidos */
/*               por default rtolx = rtoly[i][j] = 0.0,   */
/*               atolx = atoly[i][j] = 1.0e-15 e          */
/*               ftol[[i] = 1.0e-12                       */
/*                                                        */
/*               se infoinput[3] = 1 os valores de atolx, */
/*               rtolx, atoly, rtoly e ftol sao definidos */
/*               pelo usuario com escalares, isto e,      */
/*               atoly[i][j] = atolx, rtoly[i][j] = rtolx */
/*               e ftol[i] = ftol[1]                      */
/*                                                        */
/*               se infoinput[3] = 2 todos os valores de  */
/*               atolx, rtolx, atoly, rtoly e ftol sao    */
/*               definidos pelo usuario                   */
/*                                                        */
/*  infoinput  : informacoes sobre a entrada de dados na  */
/*               rotina GSDAE - vetor de dimensao 2n+10   */
/*               de valores inteiros                      */
/*                                                        */
/* infoinput[1]: infoinput[1] = 0 indica ao GSDAE que e a */
/*               primeira vez que a GSDAE e chamada       */
/*                                                        */
/*               infoinput[1] = 1 indica ao GSDAE que nao */
/*               e a primeira vez que a GSDAE e chamada   */
/*                                                        */
/*               quando infoinput[1] = 0, a rotina GSDAE  */
/*               redefine infoinput[1] = 1                */
/*                                                        */
/* infoinput[2]: infoinput[2] = 0 indica a rotina GSDAE   */
/*               ou CSDAE que e a matriz jacobiana deve   */
/*               ser aproximada                           */
/*                                                        */
/*               infoinput[2] = 1 indica a rotina GSDAE   */
/*               ou CSDAE que deve ser utilizada a matriz */
/*               jacobiana definida pelo usuario na ro-   */
/*               tina DF                                  */
/*                                                        */
/* infoinput[3]: infoinput[3] = 0 indica a rotina que as  */
/*               tolernacias atol, rtol e ftol devem ser  */
/*               definidas por default                    */
/*                                                        */
/*               infoinput[3] = 1 indica a rotina que as  */
/*               tolernacias atol, rtol e ftol foram de-  */
/*               finidas pelo usuario como escalares      */
/*                                                        */
/*               infoinput[3] = 2 indica a rotina que as  */
/*               tolernacias atol, rtol e ftol foram de-  */
/*               finidas pelo usuario como vetoriais      */
/*                                                        */
/* infoinput[4]: infoinput[4] = 0 indica a rotina que     */
/*               o posto de DFy[o] deve ser inicializado  */
/*               como sendo maximo, isto e, posto = n, e  */
/*               as permutacoes de coordendas da funcao e */
/*               das variaveis y[0],..,y[o] devem ser     */
/*               inicializadas como sendo a identidade    */
/*                                                        */
/*               infoinput[4] > 0 indica a rotina que     */
/*               o posto de DFy[o] e infoinput[4] e as    */
/*               permutacoes de coordenadas da funcao     */
/*               estao em infoinput[10+i] (i = 1..n) e as */
/*               permutacoes das variaveis y[0],..,y[o]   */
/*               estao em infoinput[10+n+i] (i = 1..n)    */
/*                                                        */
/* infoinput[i]: i = 11..10+n armazena as permutacoes de  */
/*               coordenadas da funcao que define a EAD   */
/*               quando infoinput[0] > 0                  */
/*                                                        */
/* infoinput[i]: i = 11+n..10+2n armazena as permutacoes  */
/*               das variaveis y[0],..,y[o]               */
/*               quando infoinput[0] > 0                  */
/*                                                        */
/* infoinput[i]: i = 0,5..10 nao sao utilizadas nesta     */
/*               versao                                   */
/*                                                        */
/*                                                        */
/*  Parametros de saida                                   */
/*                                                        */
/* infooutput   : informacoes sobre a saida de dados da   */
/*                rotina - vetor de dimensao 2n+10 de     */
/*                valores inteiros                        */
/*                                                        */
/* infooutput[1]: informa o estado da integracao na saida */
/*                da rotina - contem o ultimo valor       */
/*                retornado pela rotina                   */
/*                                                        */
/* infooutput[2]: informa o posto de DFy[o]               */
/*                                                        */
/* infooutput[3]: informa a ordem da equacao              */
/*                                                        */
/* infoinput[i]: i = 11..10+n informa  as permutacoes de  */
/*               coordenadas da funcao que define a EAD   */
/*               quando infoinput[0] > 0                  */
/*                                                        */
/* infoinput[i]: i = 11+n..10+2n informa  as permutacoes  */
/*               das variaveis y[0],..,y[o]               */
/*               quando infoinput[0] > 0                  */
/*                                                        */
/* infoinput[i]: i = 0,5..10 nao sao utilizadas nesta     */
/*               versao                                   */
/*                                                        */
/* ****************************************************** */

int 
GSDAE (
int    n,           /* dimensao da EAD                   */
int    o,           /* ordem da EAD                      */
real   h,           /* passo de integracao               */
real   hmin,        /* passo minimo de integracao        */
real   hmax,        /* passo maximo de integracao        */
real   cdmax,       /* condicao maxima para a derivada   */
real  *s,           /* c(s) = (x(s),y(s))                */
real   send,        /* ponto final de integracao         */
real  *x,           /* c(s) = (x(s),y(s))                */
mreal  y,           /* c(s) = (x(s),y(s))                */
real   atolx,       /* erro absoluto para x              */
mreal  atoly,       /* erro absoluto para y              */
real   rtolx,       /* erro realtivo para x              */
mreal  rtoly,       /* erro realtivo para y              */
vreal  ftol,        /* erro absoluto para a funcao       */
vint   infoinput,   /* informacoes para entrada de dados */
vint   infooutput   /* informacoes para saida de dados   */
)
{
  /* declarando variaveis locais estaticas para controle de erro */
  static int error, nerror;

  /* declarando variaveis locais de controle */
  int  converg, success, success0, success1, fail;

  /* declarando o contador local de passos */
  int  nstep;

  /* declarando controladores de lacos */
  int  i, j, k;


  /* verificando se foi alocado espaco para os dados */
  if (par == NULL) {

    /* dados nao alocadaos */
    return (-1);

  }

  /* verificando se e a primeira chamada do gsdae */
  if (infoinput[1] == 0) {

    /* e a primeira chamada do gsdae */
    /* checar e inicializar os dados */

    /* definindo como nao sendo a primeira chamada */
    infoinput[1] = 1;

    /* verificar a dimensao e a ordem da EAD */
    if ((o < 0) || (n <= 0)) {

      /* erro na entrada de dados */
      return (-2);

    } 

    /* inicializar a dimensao e a ordem */
    par->n     = n;
    par->o     = o;

    /* inicializar s e c(s) = (x(s),y(s)) */
    par->s     = *s;
    par->cx    = *x;
    for (j = 1; j <= par->n; j++) 
      for (i = 0; i <= par->o; i++) 
        par->cy[i][j] = y[i][j];

    /* verificar o tamanho do passo */
    if (h == 0.0) {
      
      /* definindo h */
      if (hmin == 0.0) {

	/* definindo hmin e h */
	par->hmin = 1.0e-16;
	par->h    = 1.0e-15;

      } else {

	/* definindo h */
	par->hmin = fabs(hmin);
	par->h    = 10.0*par->hmin;

      }

    } else if (fabs(h) > fabs(hmin)) {

      /* definindo h e hmin */
      par->h     = fabs(h);
      par->hmin  = fabs(hmin);

    } else {

      /* erro na entrada de dados */
      return (-2);

    }


    /* definindo hmax */
    par->hmax  = fabs(hmax);

    /* verificando hmax */
    if ((par->hmax > 0.0) && (par->hmax < par->h)) {

      /* erro na entrada de dados */
      return (-2);

    }

    /* definindo a condicao maxima */
    par->cdmax = fabs(cdmax);

    if (par->cdmax == 0.0) {

      /* definindo a condicao maxima com 1.0e6 */
      par->cdmax = 1.0e6;

    } else if (par->cdmax < 1.0e2) {

      /* erro na entrada de dados */
      return (-2);

    }

    /* definindo o tipo de avaliacao da jacobiana */
    if (infoinput[2] == 0) {

      /* a jacobina sera aproximada */
      par->nDH = 0;

    } else if (infoinput[2] == 1) {

      /* a jacobiana e definida na rotina DF */
      par->nDH = 1;

    } else {

      /* erro na entrada de dados */
      return (-2);

    }

    /* definindo as tolerancias */
    if (infoinput[3] == 0) {

      /* definindo os erros absolutos com 1.0e-15 e */
      /* os erros relativos com 0.0 e               */
      /* o erro na avaliacao da funcao com 1.0e-12  */

      /* definindo o erro absoluto e relativo para x  */
      par->atolx = 1.0e-15;
      par->rtolx = 0.0;
      for (j = 1; j <= par->n; j++) {
        /* definindo o erro absoluto para a funcao  */
        par->ftol[j] = 1.0e-12;
        /* definindo o erro absoluto e relativo para y  */
        for (i = 0; i <= par->o; i++) {
          par->atoly[i][j] = 1.0e-15;
          par->rtoly[i][j] = 0.0;
        }
      } 

    } else if (infoinput[3] == 1) {

      /* os erros sao definidos sendo todos iguais */
      /* segundo atolx, rtolx e ftol[1]            */
      /* caso ftol[1] = 0.0 nao e checada a funcao */
      par->atolx   = fabs(atolx);
      par->rtolx   = fabs(rtolx);
      par->ftol[1] = fabs(ftol[1]);

      /* verificando o erro relativo e absoluto */
      if (par->atolx > par->rtolx) {

        /* erro na entrada de dados */
        return (-2);

      }

      /* definindo as demais tolerancias */
      for (j = 1; j <= par->n; j++) {
        /* definindo o erro absoluto para a funcao  */
        par->ftol[j] = par->ftol[1];
        /* definindo o erro absoluto e relativo para y  */
        for (i = 0; i <= par->o; i++) {
          par->atoly[i][j] = par->atolx;
          par->rtoly[i][j] = par->rtolx;
        }
      } 

    } else if (infoinput[3] == 2) {

      /* os erros sao definidos pelo usuario podendo      */
      /* serem diferentes para cada variavel e coordenada */
      /* da funcao                                        */
      /* caso ftol[i] = 0.0 nao e checada esta coordenada */
      /* da funcao                                        */

      /* definindo o erro absoluto e relativo para x  */
      par->atolx   = fabs(atolx);
      par->rtolx   = fabs(rtolx);

      /* verificando o erro relativo e absoluto */
      if (par->atolx > par->rtolx) {

        /* erro na entrada de dados */
        return (-2);

      }

      for (j = 1; j <= par->n; j++) {

        /* definindo o erro absoluto para a funcao  */
        par->ftol[j] = fabs(ftol[j]);

        /* definindo o erro absoluto e relativo para y  */
        for (i = 0; i <= par->o; i++) {

          par->atoly[i][j] = fabs(atoly[i][j]);
          par->rtoly[i][j] = fabs(rtoly[i][j]);

          /* verificando o erro relativo e absoluto */
          if (par->atoly[i][j] > par->rtoly[i][j]) {

            /* erro na entrada de dados */
            return (-2);

          }

        }
      } 

    } else {

      /* erro na entrada de dados */
      return (-2);

    }

    /* definindo que as tolerancias nao devem ser mudadas */
    infoinput[3] = 0;

    /* verificando se o posto foi definido */
    if (infoinput[4] == 0) {

      /* definindo o posto sendo maximo e as */
      /* permutacoes sendo a identidade      */

      /* definindo a opcao informada */
      par->irank = 0;

      /* definindo o posto */
      par->rank  = par->n;

      /* definindo as permutacoes */
      for (i = 1; i <= n; i++) {
        par->p[i] = i;
        par->q[i] = i;
      }

      /* definindo fail no calculo do vetor tangente */
      fail = 0;

    } else if ((infoinput[4] > 0) && (infoinput[4] <= n)) {

      /* o posto e as permutacoes foram definidos */
      /* pelo usuario                             */

      /* definindo a opcao informada */
      par->irank = 1;

      /* definindo o posto */
      par->rank  = infoinput[4];

      /* definindo as permutacoes */
      for (i = 1; i <= n; i++) {
        par->p[i] = infoinput[10+i];
        par->q[i] = infoinput[10+n+i];
      }

      /* verificando a consistencia dos dados  */
      /* da permutacao de coordendas da funcao */
      for (i = 1; i <= n; i++) {

        /* definindo que ha falha */
        fail = 1;

        /* verificando se i pertence a p */
        for (j = 1; (j <= n) && fail; j++)
          if (par->p[j] == i) 
            fail = 0;

        /* verificando se houve falha */
        if (fail) {

          /* erro na entrada de dados */
          return (-2);

        }

      }

      /* verificando a consistencia dos dados     */
      /* da permutacao de coordendas de variaveis */
      for (i = 1; i <= n; i++) {

        /* definindo que ha falha */
        fail = 1;

        /* verificando se i pertence a q */
        for (j = 1; (j <= n) && fail; j++)
          if (par->q[j] == i) 
            fail = 0;

        /* verificando se houve falha */
        if (fail) {

          /* erro na entrada de dados */
          return (-2);

        }

      }

      /* definindo fail no calculo do vetor tangente */
      fail = 1;

    } else {

      /* erro na entrada de dados */
      return (-2);

    }

    /* inicializando os contadores de passos,     */
    /* passos falhos, falhas no metodo de Newton, */
    /* numero de avaliacao da funcao, numero de   */
    /* avaliacao da jacobiana, numero de          */
    /* decomposicoes QR                           */
    par->nstep  = 0;
    par->nff    = 0;
    par->nflhs  = 0;
    par->nNwf   = 0;
    par->naF    = 0;
    par->naDH   = 0;
    par->ndQR   = 0;
    par->nstart = 0;

    /* definindo a direcao de integracao de acordo */
    /* com o ponto final                           */
    if ((*s) <= send) {

      par->dir  = 1.0; 

    } else {

      par->dir  = -1.0;
      par->h   *= -1.0;

    }

    /* armazenando o passo inicial */
    par->h0    = par->h;

    /* definindo o vetor tangente ao ponto inicial e     */
    /* caso haja queda de posto definindo o novo posto e */
    /* permutacoes                                       */
    success0 = settau(par->n,&(par->o),&(par->rank),par->irank,
                      par->cx,par->cy,par->u,par->deltax,par->DFx,par->DFy,
                      par->cdmax,par->dir,par->Q,par->DH,par->p,par->q,
                      par->paux,par->qaux,par->ftol,par->atolx,par->rtolx,
                      &(par->x),par->y,par->nDH,&(par->naF),
                      &(par->naDH),&(par->ndQR),par->F,par->DF);
    (par->nstart)++;

    /* escrevendo no vetor de saida de comunicacao */
    /* com o usuario                               */

    /* escrevendo o tipo de ponto inicial ou o */
    /* erro do dado inicial                    */
    infooutput[1] = success0;

    /* escrevendo o posto e a ordem */
    infooutput[2] = par->rank;
    infooutput[3] = par->o;

    /* escrevendo as permutacoes */
    for (i = 1; i <= n; i++) {
      infooutput[10+i]   = par->p[i];
      infooutput[10+n+i] = par->q[i];
    }

    /* redefinindo a opcao para futuras quedas de posto */
    par->irank = 0;

    /* inicializando os dados para o laco principal      */
    /* os dados sao: coeficientes e diferencas divididas */
    /* para os polinomios preditor e corretor            */
    firststep(par->n,par->o,par->rank,par->cx,par->cy,par->phix,par->phiy,
              par->p,par->q,par->x,par->y,&(par->cj),&(par->cjold),
              &(par->factor),&(par->hold),par->h,&(par->k),&(par->kold),
              &(par->ns),par->psi,&(par->ifase),&(par->dx)); 

    /* definindo o contador de erros para o laco principal */
    nerror = 0;

    /* abortando no caso de erro na definicao do vetor tangente */
    if (success0 != 0) {

      /* retornado o erro no ponto inicial */
      return (success0);

    }

  } else if (infoinput[1] == 1) {

    /* nao e a primeira chamada do gsdae */

    /* verificando se houve erro anterior */
    if (infooutput[1] == -15) {

      /* houve erro nas duas ultimas chamadas do GSDAE e */
      /* nenhuma providencia foi tomada                  */
      /* abortando a execucao                            */
      printf("The two last step was terminated with a negative value and\n");
      printf("no appropriate action was taken");
      exit (-15);

    } else if (infooutput[1] < 0) {

      /* houve erro na ultima chamada do GSDAE e */
      /* nenhuma providencia foi tomada          */
      infooutput[1] = -15;
      return (-15);

    }

    /* verificando se nao houve mudanca no espacamento minimo*/
    if ((fabs(hmin) > 0.0) && (fabs(hmin) != par->hmin)) {

      par->hmin = fabs(hmin);

    }

    /* verificando se nao houve mudanca no espacamento maximo*/
    if ((fabs(hmax) > 0.0) && (fabs(hmax) != par->hmax)) {

      par->hmax = fabs(hmax);

    }

    /* verificando se nao houve mudanca na condicao maxima */
    if (fabs(cdmax) != par->cdmax) {

      par->cdmax = fabs(cdmax);

    }

    if (infoinput[3] == 2) {

      /* os erros sao redefinidos pelo usuario podendo    */
      /* serem diferentes para cada variavel e coordenada */
      /* da funcao                                        */
      /* caso ftol[i] = 0.0 nao e checada esta coordenada */
      /* da funcao                                        */

      /* definindo o erro absoluto e relativo para x  */
      par->atolx   = fabs(atolx);
      par->rtolx   = fabs(rtolx);

      /* verificando o erro relativo e absoluto */
      if (par->atolx > par->rtolx) {

        /* erro na entrada de dados */
        return (-2);

      }

      for (j = 1; j <= par->n; j++) {

        /* definindo o erro absoluto para a funcao  */
        par->ftol[j] = fabs(ftol[j]);

        /* definindo o erro absoluto e relativo para y  */
        for (i = 0; i <= par->o; i++) {

          par->atoly[i][j] = fabs(atoly[i][j]);
          par->rtoly[i][j] = fabs(rtoly[i][j]);

          /* verificando o erro relativo e absoluto */
          if (par->atoly[i][j] > par->rtoly[i][j]) {

            /* erro na entrada de dados */
            return (-2);

          }

        }
      } 

    }

  } else {

    /* erro na entrada de dados */
    return (-2);

  }

  /* definindo a nao convergencia */
  converg = 0;

  /* inicializando o contador de passos */
  nstep   = 0;

  /* definindo que nao ha erro */
  error   = 0;

  /* verificando a convergencia:                      */
  /* se o ponto final nao e posterior ao ponto atual  */
  /* basta interpolar a solucao                       */
  if ((par->s-send)*par->dir >= 0.0) {

    /* definindo a convergencia */
    converg = 1;

  }   

  /* enquanto nao houver convergencia ou ocorrer algum erro */
  while (!converg) {  

    /* passo principal na integracao da EAD */
    /* este passo define o proximo ponto    */
    success = masterstep(par->n,par->o,par->rank,&(par->h),&(par->s),
              &(par->k),&(par->nNwf),&(par->aDH),par->nDH,
              &(par->x),par->y,&(par->cx),par->cy,
              &(par->cxx), par->cyx,&(par->pcx),par->pcy,&(par->pdcx),
              par->pdcy,par->delta,par->p,par->q,par->Q,par->DH,
	      par->alfa,par->beta,par->gama,par->sigma,par->psi,
	      &(par->alfas),&(par->cj),&(par->cjold),&(par->factor),
              &(par->ck),par->DFx,par->DFy,par->u,&(par->ccx),par->ccy,
	      par->yx,par->deltax,&(par->SP),&(par->S),
	      par->phix,par->phiy,par->hmin,par->cdmax,
              par->atolx, par->rtolx, par->atoly,par->rtoly,par->ftol,
              &(par->nff),&(par->wtx),par->wty,&(par->Ex),par->Ey, 
	      &(par->ifase),&(par->ns),&(par->hold),&(par->kold),
              par->F,par->DF,&(par->sold),&(par->taux),par->tauy,
              &(par->naF),&(par->naDH),&(par->ndQR), 
              par->deltah,par->deltahx);

    /* verificar se houve sucesso no passo masterstep */
    if (success == 1) {

      /* houve sucesso no passo masterstep */

      /* passo que controla a ordem e o tamanho do passo */
      /* de integracao                                   */
      success1 = controlstep(par->o,par->n,par->rank,
                            &(par->kold),&(par->k),par->ck,&(par->hold),
                            &(par->h),&(par->sold),&(par->s),&(par->cx),
                            par->cy,&(par->x),par->y,&(par->cxx),par->cyx,
			    par->p,par->q,par->beta,par->sigma,par->phix,
			    par->phiy,par->psi,&(par->nflhs),&(par->cfalhas),
			    par->atolx,par->rtolx,par->atoly,par->rtoly,
			    par->wtx,par->wty,par->Ex,par->Ey,&(par->ifase),
                            par->hmin,par->hmax,par->ns);
               
      /* verificar se o passo foi aceito na rotina controlstep */
      if (success1 < 0) {  

        /* falha na determinacao da amplitude do passo */
        success = -1;

      } else if (success1 == 1) {  

        /* definindo a nao ocorrencia de erro nas rotinas masterstep */
        /* e controlstep                                             */
        error = 0;

      } else {

        /* o passo nao foi aceito                              */
        /* repetir o passo com a nova ordem e tamanho de passo */
        success = 0;

      }

    }

    /* verificando se houve sucesso nas rotinas masterstep e controlstep */
    if (success == 1) {

      /* houve sucesso nas rotinas masterstep e controlstep */

      /* verificando a ordem do metodo no passo dado */
      /* se a ordem e 1 a derivada fica inalterada   */
      if (par->kold > 1) {

        /* interpole a solucao para obter a derivada o ponto */
        /* anterior                                          */
        interpolator(par->n,par->o,par->rank,par->sold-par->s,
                     &(par->cxx),par->cyx,&(par->dx),par->dy,
                     par->p,par->q,par->kold,par->phix,par->phiy,
                     par->psi);

      }
  
      /* verificando a ocorrencia de uma singularidade transversal */
      if ((par->pdcx == 0.0) || (par->dx*par->pdcx < 0.0)) {
  
        /* ocorreu uma singularidade transversal entre o */
        /* ultimo passo                                  */
  
	/* verificando se o ponto anterior tambem era uma */
	/* singularidade transversal                      */
	if (par->dx != 0.0) {

          /* calcule uma aproximacao para a singularidade  */
          /* tansversal                                    */
          interpolatorsing(par->n,par->o,par->rank,par->sold,par->s,
                           par->cxx,par->cx,par->dx,par->pdcx,
                           par->p,par->q,par->kold,
                           par->phix,par->phiy,par->psi,par->rtolx,
                           s,x,y,&(par->dx),par->dy); 


          /* definindo a derivada no ponto */
          par->dx = par->pdcx;

          /* incrementando o contador local de passos */
          nstep++;

          /* atualizando o contador de passos global */
          par->nstep += nstep;
  
          /* retornar a ocorrencia de uma singularidade trnasversal */
          return (1);
        }

      }

      /* definindo a derivada no ponto */
      par->dx = par->pdcx;

    } 

    if (success < 0) {

      /* houve falha na rotina masterstep ou controlstep */
      /* reinicializa o metodo com h inicial             */

      /* restaura phi */
      for (k = par->ns+1; k <= par->k+1; k++) {  
        par->phix[k] /= par->beta[k];
        for (i = 0; i <= o; i++)
          for (j = 1; j <= n; j++)
            par->phiy[k][i][j] /= par->beta[k];
      }

      /* restaura psi */  
      for (i = 2; i <= par->k+1; i++)  
        par->psi[i-1] = par->psi[i]-par->h; 

      /* retornando ao ponto anterior */
      par->cx = par->cxx;
      for (i = 0; i <= o; i++)   
        for (j = 1; j <= n; j++) 
          par->cy[i][j] = par->cyx[i][j];
   
      /* verificar se e a primeira falha */
      if (!error) {

        /* definindo o passo como sendo o passo inicial */
        par->h = par->h0;

        /* recalculando o posto, as permutacoes e a */
	/* tangente no ponto                        */
        success0 = settau(par->n,&(par->o),&(par->rank),par->irank,
                          par->cx,par->cy,par->u,par->deltax,
                          par->DFx,par->DFy,par->cdmax,par->dir,par->Q,
                          par->DH,par->p,par->q,par->paux,par->qaux,
                          par->ftol,par->atolx,par->rtolx,
                          &(par->x),par->y,par->nDH,&(par->naF),
                          &(par->naDH),&(par->ndQR),par->F,par->DF);
        (par->nstart)++;

        /* escrevendo os novos dados */
        infooutput[1] = success0;
        infooutput[2] = par->rank;
        infooutput[3] = par->o;
        for (i = 1; i <= n; i++) {
          infooutput[10+i]   = par->p[i];
          infooutput[10+n+i] = par->q[i];
        }


      }

      /* se e a segunda ocorrencia consecutiva de erro */
      /* no passo masterstep, ou se ocorreram mais que */
      /* quatro falhas, ou se nao houve sucesso na     */
      /* definicao da tangente, o programa e abortado  */
      if (error || (nerror > 4) || (success0 < 0)) {

        /* copiando o ultimo ponto */
        *x         = par->cx;
        for (i = 0; i <= o-1; i++)
          for (j = 1; j <= n; j++)
            y[i][j] = par->cy[i][j]; 
        for (j = 1; j <= par->rank; j++)
          y[o][j] = par->cy[o][j]; 
        if (par->o > 0)
          for (i = par->rank+1; i <= par->n; i++)
            y[o][i] = par->dy[o-1][i]/par->dx; 

        /* incrementando o contador local de passos */
        nstep++;

        /* atualizando o contador de passos global */
        par->nstep += nstep;

        /* retornar o erro ocorrido */
        if (success0 < 0) {

          if (success0 == -3)
	    success0 = -11;

          infooutput[1] = success0;

          return (success0);

        } else {

          infooutput[1] = success;

          return (success);

        }

      }

      /* inicializar os dados para um recomeco */ 
      firststep(par->n,par->o,par->rank,par->cx,par->cy,par->phix,par->phiy,
              par->p,par->q,par->x,par->y,&(par->cj),&(par->cjold),
              &(par->factor),&(par->hold),par->h,&(par->k),&(par->kold),
              &(par->ns),par->psi,&(par->ifase),&(par->dx)); 

      /* definindo a ocorrencia de erro */
      error  = 1;

      /* incrementando o contador de erros */
      nerror++;

      /* informando ao usuario a mudanca de situacao */ 
      if (success0 > 0) {

        /* incrementando o contador local de passos */
        nstep++;

        /* atualizando o contador de passos global */
        par->nstep += nstep;
  
        /* retornar a ocorrencia */
        return (success0);

      }

    }
             
     /* verificando a convergencia */
    if ((par->s-send)*par->dir >= 0.0) {

      /* definindo a convergencia */
      converg = 1;

    }   
  
    /* incrementando o contador local de passos */
    nstep++;
       
  } 
  /* fim while */

  /* atualizando o contador de passos global */
  par->nstep += nstep;

  /* interpolando a solucao em send */
  *s = send;
  interpolators(par->n,par->o,par->rank,s,x,y,&(par->dx),par->dy,
		par->p,par->q,par->kold,par->phix,par->phiy,
		par->psi,par->s); 

  /* definindo a derivada no ponto */
  par->dx = par->pdcx;

  /* retornando sucesso no passo de integracao (s = send) */
  return (0);                   
}  
/* fim GSDAE */




int 
CSDAE (
int    n,           /* dimensao da EAD                   */
int    o,           /* ordem da EAD                      */
real   h,           /* passo de integracao               */
real   hmin,        /* passo minimo de integracao        */
real   hmax,        /* passo maximo de integracao        */
real   cdmax,       /* condicao maxima para a derivada   */
real  *s,           /* c(s) = (x(s),y(s))                */
real  *x,           /* c(s) = (x(s),y(s))                */
real   xend,        /* ponto final de integracao         */
mreal  y,           /* c(s) = (x(s),y(s))                */
real   atolx,       /* erro absoluto para x              */
mreal  atoly,       /* erro absoluto para y              */
real   rtolx,       /* erro realtivo para x              */
mreal  rtoly,       /* erro realtivo para y              */
vreal  ftol,        /* erro absoluto para a funcao       */
vint   infoinput,   /* informacoes para entrada de dados */
vint   infooutput   /* informacoes para saida de dados   */
)
{
  /* declarando variaveis locais estaticas para controle de erro */
  static int error, nerror;

  /* declarando variaveis locais de controle */
  int  converg, success, success0, success1, fail;

  /* declarando o contador local de passos */
  int  nstep;

  /* declarando controladores de lacos */
  int  i, j, k;


  /* verificando se foi alocado espaco para os dados */
  if (par == NULL) {

    /* dados nao alocadaos */
    return (-1);

  }

  /* verificando se e a primeira chamada do gsdae */
  if (infoinput[1] == 0) {

    /* e a primeira chamada do gsdae */
    /* checar e inicializar os dados */

    /* definindo como nao sendo a primeira chamada */
    infoinput[1] = 1;

    /* verificar a dimensao e a ordem da EAD */
    if ((o < 0) || (n <= 0)) {

      /* erro na entrada de dados */
      return (-2);

    } 

    /* inicializar a dimensao e a ordem */
    par->n     = n;
    par->o     = o;

    /* inicializar s e c(s) = (x(s),y(s)) */
    par->s     = *s;
    par->cx    = *x;
    for (j = 1; j <= par->n; j++) 
      for (i = 0; i <= par->o; i++) 
        par->cy[i][j] = y[i][j];

    /* verificar o tamanho do passo */
    if (h == 0.0) {
      
      /* definindo h */
      if (hmin == 0.0) {

	/* definindo hmin e h */
	par->hmin = 1.0e-16;
	par->h    = 1.0e-15;

      } else {

	/* definindo h */
	par->hmin = fabs(hmin);
	par->h    = 10.0*par->hmin;

      }

    } else if (fabs(h) > fabs(hmin)) {

      /* definindo h e hmin */
      par->h     = fabs(h);
      par->hmin  = fabs(hmin);

    } else {

      /* erro na entrada de dados */
      return (-2);

    }


    /* definindo hmax */
    par->hmax  = fabs(hmax);

    /* verificando hmax */
    if ((par->hmax > 0.0) && (par->hmax < par->h)) {

      /* erro na entrada de dados */
      return (-2);

    }

    /* definindo a condicao maxima */
    par->cdmax = fabs(cdmax);

    if (par->cdmax == 0.0) {

      /* definindo a condicao maxima com 1.0e6 */
      par->cdmax = 1.0e6;

    } else if (par->cdmax < 1.0e2) {

      /* erro na entrada de dados */
      return (-2);

    }

    /* definindo o tipo de avaliacao da jacobiana */
    if (infoinput[2] == 0) {

      /* a jacobina sera aproximada */
      par->nDH = 0;

    } else if (infoinput[2] == 1) {

      /* a jacobiana e definida na rotina DF */
      par->nDH = 1;

    } else {

      /* erro na entrada de dados */
      return (-2);

    }

    /* definindo as tolerancias */
    if (infoinput[3] == 0) {

      /* definindo os erros absolutos com 1.0e-15 e */
      /* os erros relativos com 0.0 e               */
      /* o erro na avaliacao da funcao com 1.0e-12  */

      /* definindo o erro absoluto e relativo para x  */
      par->atolx = 1.0e-15;
      par->rtolx = 0.0;
      for (j = 1; j <= par->n; j++) {
        /* definindo o erro absoluto para a funcao  */
        par->ftol[j] = 1.0e-12;
        /* definindo o erro absoluto e relativo para y  */
        for (i = 0; i <= par->o; i++) {
          par->atoly[i][j] = 1.0e-15;
          par->rtoly[i][j] = 0.0;
        }
      } 

    } else if (infoinput[3] == 1) {

      /* os erros sao definidos sendo todos iguais */
      /* segundo atolx, rtolx e ftol[1]            */
      /* caso ftol[1] = 0.0 nao e checada a funcao */
      par->atolx   = fabs(atolx);
      par->rtolx   = fabs(rtolx);
      par->ftol[1] = fabs(ftol[1]);

      /* verificando o erro relativo e absoluto */
      if (par->atolx > par->rtolx) {

        /* erro na entrada de dados */
        return (-2);

      }

      /* definindo as demais tolerancias */
      for (j = 1; j <= par->n; j++) {
        /* definindo o erro absoluto para a funcao  */
        par->ftol[j] = par->ftol[1];
        /* definindo o erro absoluto e relativo para y  */
        for (i = 0; i <= par->o; i++) {
          par->atoly[i][j] = par->atolx;
          par->rtoly[i][j] = par->rtolx;
        }
      } 

    } else if (infoinput[3] == 2) {

      /* os erros sao definidos pelo usuario podendo      */
      /* serem diferentes para cada variavel e coordenada */
      /* da funcao                                        */
      /* caso ftol[i] = 0.0 nao e checada esta coordenada */
      /* da funcao                                        */

      /* definindo o erro absoluto e relativo para x  */
      par->atolx   = fabs(atolx);
      par->rtolx   = fabs(rtolx);

      /* verificando o erro relativo e absoluto */
      if (par->atolx > par->rtolx) {

        /* erro na entrada de dados */
        return (-2);

      }

      for (j = 1; j <= par->n; j++) {

        /* definindo o erro absoluto para a funcao  */
        par->ftol[j] = fabs(ftol[j]);

        /* definindo o erro absoluto e relativo para y  */
        for (i = 0; i <= par->o; i++) {

          par->atoly[i][j] = fabs(atoly[i][j]);
          par->rtoly[i][j] = fabs(rtoly[i][j]);

          /* verificando o erro relativo e absoluto */
          if (par->atoly[i][j] > par->rtoly[i][j]) {

            /* erro na entrada de dados */
            return (-2);

          }

        }
      } 

    } else {

      /* erro na entrada de dados */
      return (-2);

    }

    /* definindo que as tolerancias nao devem ser mudadas */
    infoinput[3] = 0;

    /* verificando se o posto foi definido */
    if (infoinput[4] == 0) {

      /* definindo o posto sendo maximo e as */
      /* permutacoes sendo a identidade      */

      /* definindo a opcao informada */
      par->irank = 0;

      /* definindo o posto */
      par->rank  = par->n;

      /* definindo as permutacoes */
      for (i = 1; i <= n; i++) {
        par->p[i] = i;
        par->q[i] = i;
      }

      /* definindo fail no calculo do vetor tangente */
      fail = 0;

    } else if ((infoinput[4] > 0) && (infoinput[4] <= n)) {

      /* o posto e as permutacoes foram definidos */
      /* pelo usuario                             */

      /* definindo a opcao informada */
      par->irank = 1;

      /* definindo o posto */
      par->rank  = infoinput[4];

      /* definindo as permutacoes */
      for (i = 1; i <= n; i++) {
        par->p[i] = infoinput[10+i];
        par->q[i] = infoinput[10+n+i];
      }

      /* verificando a consistencia dos dados  */
      /* da permutacao de coordendas da funcao */
      for (i = 1; i <= n; i++) {

        /* definindo que ha falha */
        fail = 1;

        /* verificando se i pertence a p */
        for (j = 1; (j <= n) && fail; j++)
          if (par->p[j] == i) 
            fail = 0;

        /* verificando se houve falha */
        if (fail) {

          /* erro na entrada de dados */
          return (-2);

        }

      }

      /* verificando a consistencia dos dados     */
      /* da permutacao de coordendas de variaveis */
      for (i = 1; i <= n; i++) {

        /* definindo que ha falha */
        fail = 1;

        /* verificando se i pertence a q */
        for (j = 1; (j <= n) && fail; j++)
          if (par->q[j] == i) 
            fail = 0;

        /* verificando se houve falha */
        if (fail) {

          /* erro na entrada de dados */
          return (-2);

        }

      }

    } else {

      /* erro na entrada de dados */
      return (-2);

    }

    /* inicializando os contadores de passos,     */
    /* passos falhos, falhas no metodo de Newton, */
    /* numero de avaliacao da funcao, numero de   */
    /* avaliacao da jacobiana, numero de          */
    /* decomposicoes QR                           */
    par->nstep  = 0;
    par->nff    = 0;
    par->nflhs  = 0;
    par->nNwf   = 0;
    par->naF    = 0;
    par->naDH   = 0;
    par->ndQR   = 0;
    par->nstart = 0;

    /* definindo o vetor tangente ao ponto inicial e     */
    /* caso haja queda de posto definindo o novo posto e */
    /* permutacoes                                       */
    success0 = settau(par->n,&(par->o),&(par->rank),par->irank,
                      par->cx,par->cy,par->u,par->deltax,par->DFx,par->DFy,
                      par->cdmax,par->dir,par->Q,par->DH,par->p,par->q,
                      par->paux,par->qaux,par->ftol,par->atolx,par->rtolx,
                      &(par->x),par->y,par->nDH,&(par->naF),
                      &(par->naDH),&(par->ndQR),par->F,par->DF);
    (par->nstart)++;

    /* escrevendo no vetor de saida de comunicacao */
    /* com o usuario                               */

    /* redefinindo a opcao para futuras quedas de posto */
    par->irank = 0;

    /* escrevendo o tipo de ponto inicial ou o */
    /* erro do dado inicial                    */
    infooutput[1] = success0;

    /* escrevendo o posto e a ordem */
    infooutput[2] = par->rank;
    infooutput[3] = par->o;

    /* escrevendo as permutacoes */
    for (i = 1; i <= n; i++) {
      infooutput[10+i]   = par->p[i];
      infooutput[10+n+i] = par->q[i];
    }

    /* definindo a direcao de integracao de acordo */
    /* com o ponto final                           */
    if ((xend - (*x))*par->x >= 0.0) {

      par->dir  = 1.0; 

    } else {

      par->dir  = -1.0;
      par->h   *= -1.0;

    }

    /* armazenando o passo inicial */
    par->h0    = par->h;

    /* inicializando os dados para o laco principal      */
    /* os dados sao: coeficientes e diferencas divididas */
    /* para os polinomios preditor e corretor            */
    firststep(par->n,par->o,par->rank,par->cx,par->cy,par->phix,par->phiy,
              par->p,par->q,par->x,par->y,&(par->cj),&(par->cjold),
              &(par->factor),&(par->hold),par->h,&(par->k),&(par->kold),
              &(par->ns),par->psi,&(par->ifase),&(par->dx)); 

    /* definindo o contador de erros para o laco principal */
    nerror = 0;

    /* abortando no caso de erro na definicao do vetor tangente */
    if (success0 != 0) {

      /* retornado o erro no ponto inicial */
      return (success0);

    }

  } else if (infoinput[1] == 1) {

    /* nao e a primeira chamada do gsdae */

    /* verificando se houve erro anterior */
    if (infooutput[1] == -15) {

      /* houve erro nas duas ultimas chamadas do GSDAE e */
      /* nenhuma providencia foi tomada                  */
      /* abortando a execucao                            */
      printf("The two last step was terminated with a negative value and\n");
      printf("no appropriate action was taken\n");
      exit (-15);

    } else if (infooutput[1] == -16) {

      /* ocorreu uma singularidade transversal na penultima */
      /* chamada do CSDAE e nenhuma providencia foi tomada  */
      /* abortando a execucao                               */
      printf("The last step was terminated with a transversal \n");
      printf("singularity and no appropriate action was taken\n");
      exit (-16);

    } else if (infooutput[1] < 0) {

      /* houve erro na ultima chamada do GSDAE e */
      /* nenhuma providencia foi tomada          */
      infooutput[1] = -15;
      return (-15);

    } else if ((infooutput[1] == 1) ||
               (infooutput[1] == 3) ||
               (infooutput[1] == 5)) {

      /* ocorreu uma singularidade transversal na ultima   */
      /* chamada do CSDAE e nenhuma providencia foi tomada */
      return (-16);

    }

    /* verificando se nao houve mudanca no espacamento minimo*/
    if ((fabs(hmin) > 0.0) && (fabs(hmin) != par->hmin)) {

      par->hmin = fabs(hmin);

    }

    /* verificando se nao houve mudanca no espacamento maximo*/
    if ((fabs(hmax) > 0.0) && (fabs(hmax) != par->hmax)) {

      par->hmax = fabs(hmax);

    }

    /* verificando se nao houve mudanca na condicao maxima */
    if (fabs(cdmax) != par->cdmax) {

      par->cdmax = fabs(cdmax);

    }

    if (infoinput[3] == 2) {

      /* os erros sao redefinidos pelo usuario podendo    */
      /* serem diferentes para cada variavel e coordenada */
      /* da funcao                                        */
      /* caso ftol[i] = 0.0 nao e checada esta coordenada */
      /* da funcao                                        */

      /* definindo o erro absoluto e relativo para x  */
      par->atolx   = fabs(atolx);
      par->rtolx   = fabs(rtolx);

      /* verificando o erro relativo e absoluto */
      if (par->atolx > par->rtolx) {

        /* erro na entrada de dados */
        return (-2);

      }

      for (j = 1; j <= par->n; j++) {

        /* definindo o erro absoluto para a funcao  */
        par->ftol[j] = fabs(ftol[j]);

        /* definindo o erro absoluto e relativo para y  */
        for (i = 0; i <= par->o; i++) {

          par->atoly[i][j] = fabs(atoly[i][j]);
          par->rtoly[i][j] = fabs(rtoly[i][j]);

          /* verificando o erro relativo e absoluto */
          if (par->atoly[i][j] > par->rtoly[i][j]) {

            /* erro na entrada de dados */
            return (-2);

          }

        }
      } 

    }

  } else {

    /* erro na entrada de dados */
    return (-2);

  }

  /* definindo a nao convergencia */
  converg = 0;

  /* inicializando o contador de passos */
  nstep   = 0;

  /* definindo que nao ha erro */
  error   = 0;

  /* verificando a convergencia:                      */
  /* se o ponto final nao e posterior ao ponto atual  */
  /* basta interpolar a solucao                       */
  if ((xend-par->cx)*par->dir >= 0.0) {

    /* definindo a convergencia */
    converg = 1;

  }   

  /* enquanto nao houver convergencia ou ocorrer algum erro */
  while (!converg) {  

    /* passo principal na integracao da EAD */
    /* este passo define o proximo ponto    */
    success = masterstep(par->n,par->o,par->rank,&(par->h),&(par->s),
              &(par->k),&(par->nNwf),&(par->aDH),par->nDH,
              &(par->x),par->y,&(par->cx),par->cy,
              &(par->cxx), par->cyx,&(par->pcx),par->pcy,&(par->pdcx),
              par->pdcy,par->delta,par->p,par->q,par->Q,par->DH,
	      par->alfa,par->beta,par->gama,par->sigma,par->psi,
	      &(par->alfas),&(par->cj),&(par->cjold),&(par->factor),
              &(par->ck),par->DFx,par->DFy,par->u,&(par->ccx),par->ccy,
	      par->yx,par->deltax,&(par->SP),&(par->S),
	      par->phix,par->phiy,par->hmin,par->cdmax,
              par->atolx, par->rtolx, par->atoly,par->rtoly,par->ftol,
              &(par->nff),&(par->wtx),par->wty,&(par->Ex),par->Ey, 
	      &(par->ifase),&(par->ns),&(par->hold),&(par->kold),
              par->F,par->DF,&(par->sold),&(par->taux),par->tauy,
              &(par->naF),&(par->naDH),&(par->ndQR), 
              par->deltah,par->deltahx);

    /* verificar se houve sucesso no passo masterstep */
    if (success == 1) {

      /* houve sucesso no passo masterstep */

      /* passo que controla a ordem e o tamanho do passo */
      /* de integracao                                   */
      success1 = controlstep(par->o,par->n,par->rank,
                            &(par->kold),&(par->k),par->ck,&(par->hold),
                            &(par->h),&(par->sold),&(par->s),&(par->cx),
                            par->cy,&(par->x),par->y,&(par->cxx),par->cyx,
			    par->p,par->q,par->beta,par->sigma,par->phix,
			    par->phiy,par->psi,&(par->nflhs),&(par->cfalhas),
			    par->atolx,par->rtolx,par->atoly,par->rtoly,
			    par->wtx,par->wty,par->Ex,par->Ey,&(par->ifase),
                            par->hmin,par->hmax,par->ns);
               
      /* verificar se o passo foi aceito na rotina controlstep */
      if (success1 < 0) {  

        /* falha na determinacao da amplitude do passo */
        success = -1;

      } else if (success1 == 1) {  

        /* definindo a nao ocorrencia de erro nas rotinas masterstep */
        /* e controlstep                                             */
        error = 0;

      } else {

        /* o passo nao foi aceito                              */
        /* repetir o passo com a nova ordem e tamanho de passo */
        success = 0;

      }

    }

    /* verificando se houve sucesso nas rotinas masterstep e controlstep */
    if (success == 1) {

      /* houve sucesso nas rotinas masterstep e controlstep */

      /* verificando a ordem do metodo no passo dado */
      /* se a ordem e 1 a derivada fica inalterada   */
      if (par->kold > 1) {

        /* interpole a solucao para obter a derivada o ponto */
        /* anterior                                          */
        interpolator(par->n,par->o,par->rank,par->sold-par->s,
                     &(par->cxx),par->cyx,&(par->dx),par->dy,
                     par->p,par->q,par->kold,par->phix,par->phiy,
                     par->psi);

      }
  
      /* verificando a ocorrencia de uma singularidade transversal */
      if ((par->pdcx == 0.0) || (par->dx*par->pdcx < 0.0)) {
  
        /* ocorreu uma singularidade transversal entre o */
        /* ultimo passo                                  */
  
	/* verificando se o ponto anterior tambem era uma */
	/* singularidade transversal                      */
	if (par->dx != 0.0) {

          /* calcule uma aproximacao para a singularidade  */
          /* tansversal                                    */
          interpolatorsing(par->n,par->o,par->rank,par->sold,par->s,
                           par->cxx,par->cx,par->dx,par->pdcx,
                           par->p,par->q,par->kold,
                           par->phix,par->phiy,par->psi,par->rtolx,
                           s,x,y,&(par->dx),par->dy); 


          /* definindo a derivada no ponto */
          par->dx = par->pdcx;

          /* incrementando o contador local de passos */
          nstep++;

          /* atualizando o contador de passos global */
          par->nstep += nstep;
  
          /* retornar a ocorrencia de uma singularidade trnasversal */
          return (1);
        }

      }

      /* definindo a derivada no ponto */
      par->dx = par->pdcx;

    } 

    if (success < 0) {

      /* houve falha na rotina masterstep ou controlstep */
      /* reinicializa o metodo com h inicial             */

      /* restaura phi */
      for (k = par->ns+1; k <= par->k+1; k++) {  
        par->phix[k] /= par->beta[k];
        for (i = 0; i <= o; i++)
          for (j = 1; j <= n; j++)
            par->phiy[k][i][j] /= par->beta[k];
      }

      /* restaura psi */  
      for (i = 2; i <= par->k+1; i++)  
        par->psi[i-1] = par->psi[i]-par->h; 

      /* retornando ao ponto anterior */
      par->cx = par->cxx;
      for (i = 0; i <= o; i++)   
        for (j = 1; j <= n; j++) 
          par->cy[i][j] = par->cyx[i][j];

      /* verificar se e a primeira falha */
      if (!error) {

        /* definindo o passo como sendo o passo inicial */
        par->h = par->h0;

        /* recalculando o posto, as permutacoes e a */
	/* tangente no ponto                        */
        success0 = settau(par->n,&(par->o),&(par->rank),par->irank,
                          par->cx,par->cy,par->u,par->deltax,
                          par->DFx,par->DFy,par->cdmax,par->dir,par->Q,
                          par->DH,par->p,par->q,par->paux,par->qaux,
                          par->ftol,par->atolx,par->rtolx,
                          &(par->x),par->y,par->nDH,&(par->naF),
                          &(par->naDH),&(par->ndQR),par->F,par->DF);
        (par->nstart)++;

        /* escrevendo os novos dados */
        infooutput[1] = success0;
        infooutput[2] = par->rank;
        infooutput[3] = par->o;
        for (i = 1; i <= n; i++) {
          infooutput[10+i]   = par->p[i];
          infooutput[10+n+i] = par->q[i];
        }


      }

      /* se e a segunda ocorrencia consecutiva de erro */
      /* no passo masterstep, ou se ocorreram mais que */
      /* quatro falhas, ou se nao houve sucesso na     */
      /* definicao da tangente, o programa e abortado  */
      if (error || (nerror > 4) || (success0 < 0)) {

        /* copiando o ultimo ponto */
        *x         = par->cx;
        for (i = 0; i <= o-1; i++)
          for (j = 1; j <= n; j++)
            y[i][j] = par->cy[i][j]; 
        for (j = 1; j <= par->rank; j++)
          y[o][j] = par->cy[o][j]; 
        if (par->o > 0)
          for (i = par->rank+1; i <= par->n; i++)
            y[o][i] = par->dy[o-1][i]/par->dx; 

        /* incrementando o contador local de passos */
        nstep++;

        /* atualizando o contador de passos global */
        par->nstep += nstep;

        /* retornar o erro ocorrido */
        if (success0 < 0) {

          if (success0 == -3)
	    success0 = -11;

          infooutput[1] = success0;

          return (success0);

        } else {

          infooutput[1] = success;

          return (success);

        }

      }

      /* inicializar os dados para um recomeco */ 
      firststep(par->n,par->o,par->rank,par->cx,par->cy,par->phix,par->phiy,
                par->p,par->q,par->x,par->y,&(par->cj),&(par->cjold),
                &(par->factor),&(par->hold),par->h,&(par->k),&(par->kold),
                &(par->ns),par->psi,&(par->ifase),&(par->dx)); 

      /* definindo a ocorrencia de erro */
      error  = 1;

      /* incrementando o contador de erros */
      nerror++;

      /* informando ao usuario a mudanca de situacao */ 
      if (success0 > 0) {

        /* incrementando o contador local de passos */
        nstep++;

        /* atualizando o contador de passos global */
        par->nstep += nstep;
  
        /* retornar a ocorrencia */
        return (success0);

      }

    }
             
     /* verificando a convergencia */
    if ((xend-par->cx)*par->dir >= 0.0) {

      /* definindo a convergencia */
      converg = 1;

    }   
  
    /* incrementando o contador local de passos */
    nstep ++;
       
  } 
  /* fim while */

  /* atualizando o contador de passos global */
  par->nstep += nstep;

  /* interpolando a solucao em xend */
  *x = xend;
  interpolatorx(par->n,par->o,par->rank,par->s,par->cx,par->pdcx,
                par->kold,par->phix,par->phiy,par->psi,
                par->p,par->q,s,x,y,&(par->dx),par->dy,par->atolx); 

  /* definindo a derivada no ponto */
  par->dx = par->pdcx;

  /* retornando sucesso no passo de integracao (s = send) */
  return (0);                   
}  
/* fim CSDAE */



/*******************************************************/
/* rotina que retorna contadores para o usuario        */
/*******************************************************/

void 
STATISTICS (
real *s,
int  *nstep,
int  *nreject,
int  *nsuc,
int  *nfunc,
int  *njac,
int  *nqr,
int  *nstart,
int  *nfnew
)
{
  *s       = par->s;
  *nstep   = par->nstep;
  *nreject = par->nflhs;
  *nsuc    = par->nstep-par->nflhs;
  *nfunc   = par->naF;
  *njac    = par->naDH;
  *nqr     = par->ndQR;
  *nstart  = par->nstart;
  *nfnew   = par->nNwf;

  return;
}
                                         
  




/**********************************************************/
/* Esta rotina retorna uma menssagem de erro de acordo    */
/* com o valor da variavel status que e um parametro de   */
/* entrada. O significado da menssagem de acordo com      */
/* status e :                                             */
/*   0 : ponto regular (o valor send foi atingido)        */
/*   1 : singularidade transversal                        */
/*   2 : ponto regular e o posto de DFy[[o] diminuiu      */
/*   3 : singularidade transversal e o posto de DFy[o]    */
/*       diminuiu                                         */
/*   4 : ponto regular e a ordem da EAD diminuiu          */
/*   5 : singularidade transversal e a ordem da EAD       */
/*       diminuiu                                         */
/*  -1 : o parametro par nao foi alocado                  */
/*  -2 : erro na entrada de dados                         */
/*  -3 : ponto inicial inadequado - nao satisfaz a EAD    */
/*       com a tolerancia desejada                        */
/*  -4 : posto maior que o informado                      */
/*  -5 : EAD com ordem zero e posto de DFy[0] zero        */
/*  -6 : o posto de DFy[o] varia em uma vizinhaca do      */
/*       ponto                                            */
/*  -7 : a ordem diminuiu e o posto de DFy[o] varia em    */
/*       uma vizinhaca do ponto                           */
/*  -8 : singularidade nao tranversal                     */
/*  -9 : singularidade nao transversal com posto  de      */
/*       DFy[o] inferior                                  */
/* -10 : singularidade nao transversal com ordem inferior */
/* -11 : ponto inadequado - nao satisfaz a EAD com a      */
/*       tolerancia desejada                              */
/* -12 : h < hmin                                         */
/* -13 : condicao da jacobiana maior que a adimitida      */
/* -14 : as tolerancias exigidas nao foram atingidas      */
/* -15 : houve falha na rotina GSDAE e nehuma providencia */
/*       foi tomada                                       */
/* -16 : ocorreu uma singularidade na ultima chamada da   */
/*       rotina CSDAE e nehuma providencia foi tomada     */
/**********************************************************/

void 
STATUS (
int   status,        
char *mens
)
{
  if (status == -1) {

    /* Erro na entrada de dados */
    sprintf(mens,"Incompleat execution : parameters not allocated");

  } else if (status == -2) {

    /* Erro na entrada de dados */
    sprintf(mens,"Incompleat execution : error in data input");

  } else if (status == -3) {

    /* Tolerancia para F nao atingida */
    sprintf(mens,"Incompleat execution : initial point do not satisfies the DAE");

  } else if (status == -4) {

    /* posto maior que o informado */
    sprintf(mens,"Incompleat execution : rank(DFy[0]) is greater than the actual");
    sprintf(mens,"Incompleat execution : order 0 and rank(DFy[0]) < n-1");

  } else if (status == -5) {

    /* ordem e posto zero */
    sprintf(mens,"Incompleat execution : order 0 and rank(DFy[0]) = 0");

  } else if (status == -6) {

    /* o posto de DFy[o] varia na vizinhanca do ponto */
    sprintf(mens,"Incompleat execution : rank(DFy[o]) is not constant in neighbourhood of the point");

  } else if (status == -7) {

    /* ordem diminuiu e o posto de DFy[o] varia na vizinhanca do ponto */
    sprintf(mens,"Incompleat execution : order diminish and rank(DFy[o]) is not constant in neighbourhood of the point");

  } else if (status == -8) {

    /* singularidade nao transversal */
    sprintf(mens,"Incompleat execution : non transversal singularity");

  } else if (status == -9) {

    /* singularidade nao transversal e posto de DFy[o] diminuiu */
    sprintf(mens,"Incompleat execution : non transversal singularity and rank(DFy[o]) diminish");

  } else if (status == -10) {

    /* singularidade nao transversal e a ordem diminuiu */
    sprintf(mens,"Incompleat execution : non transversal singularity and the order diminish");

  } else if (status == -11) {

    /* Tolerancia para F nao atingida */
    sprintf(mens,"Incompleat execution : point do not satisfies the DAE");

  } else if (status == -12) {

    /* h < hmin  */
    sprintf(mens,"Incompleat execution : h < hmin");

  } else if (status == -13) {

    /* Condicao maior que a admitida */
    sprintf(mens,"Incompleat execution : big condition matrix number");

  } else if (status == -14) {

    /* Erro na rotina masterstep */
    sprintf(mens,"Incompleat execution : the corrector could not converge");

  } else if (status == -15) {

    /* houve falha na rotina GSDAE e nehuma providencia foi tomada  */
    sprintf(mens,"Incompleat execution : the last step was terminated with a negative value and no appropriate action was taken");

  } else if (status == -16) {

    /* ocorreu uma singularidade na ultima chamada da rotina CSDAE */
    /* e nehuma providencia foi tomada                             */
    sprintf(mens,"Incompleat execution : the last step was terminated with a transversal singularity and no appropriate action was taken");

  } else if (status == 0) {

    /* Execucao completa - valor final foi atingido  */
    sprintf(mens,"Compleat execution : no error");

  } else if (status == 1) {

    /* singularidade transversal */
    sprintf(mens,"Transversal singularity");

  } else if (status == 2) {

    /* ponto regular e o posto de DFy[o] diminuiu */
    sprintf(mens,"Regular point and rank(DFy[o]) diminish");

  } else if (status == 3) {

    /* singularidade transversal e o posto de DFy[o] diminuiu */
    sprintf(mens,"Transversal singularity and rank(DFy[o]) diminish");

  } else if (status == 4) {

    /* ponto regular e a ordem da EAD diminuiu */
    sprintf(mens,"Regular point and order < o");

  } else if (status == 5) {

    /* singularidade transversal e a ordem da EAD diminuiu */
    sprintf(mens,"Transversal singularity and order < o");

  } else {

    /* valor invalido para status */
    sprintf(mens,"Invalid value");

  }

  return;
}


/* ***************************************************************** */
/* calcula o vetor peso para o ponto c = (cx, cy)                    */
/* para o uso na rotina weightnorm                                   */
/* ***************************************************************** */

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
)
{  
  int i,j;

  *wtx = rtolx*fabs(cx)+atolx;

  for (i = 0; i <= o-1; i++)   
    for (j = 1; j <= n; j++)   
       wty[i][q[j]] = rtoly[i][q[j]]*fabs(cy[i][q[j]])+atoly[i][q[j]];

  for (j = 1; j <= r; j++)   
    wty[o][q[j]] = rtoly[o][q[j]]*fabs(cy[o][q[j]])+atoly[o][q[j]];

  return;
} 
/* fim weightvector */


/* *************************************** */
/* calcula a norma peso de c = (cx,cy) com */
/* relacao ao vetor peso wtc = (wtx,wty)   */         
/* *************************************** */

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
)
{
  int  i,j;
  int  neq;
  real vmax,aux;
  real norm;

  neq  = o*n+r+1;
  norm = 0.0;  

  vmax = fabs(cx/wtx); 
  for (i = 0; i <= o-1; i++) 
    for (j = 1; j <= n; j++)    
      if ((aux = fabs(cy[i][q[j]]/wty[i][q[j]])) >= vmax) 
        vmax = aux;
  for (j = 1; j <= r; j++)    
    if ((aux = fabs(cy[o][q[j]]/wty[o][q[j]])) >= vmax) 
      vmax = aux;

  if (vmax == 0.0) {

    return (0.0);

  } else {

    norm  = (cx/wtx)/vmax;  
    norm *= norm; 
    for (i = 0; i <= o-1; i++)
      for (j = 1; j <= n; j++) {
        aux    = (cy[i][q[j]]/wty[i][q[j]])/vmax;
        norm  += aux*aux;
      }
    for (j = 1; j <= r; j++) {
      aux    = (cy[o][q[j]]/wty[o][q[j]])/vmax;
      norm  += aux*aux;
    }
    norm = vmax*sqrt(norm/neq);

  }

  return (norm);   
} 
/* fim weightnorm */  



/* *********************************************************** */
/*  Rotina para a construcao de H(c) = (F(c),w(c)c',(c',c')-1) */
/* *********************************************************** */

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
)
{ 
  int  i,j;
  real norm, aux;
      
  /* deltah[i] = F(c(s))[i] (i = 1..n) */

  F(o,n,x,y,delta);  
  for (i = 1; i <= n; i++)
    deltah[i] = delta[p[i]];

  /* deltah[i] = w(c(s))c'(s)[i]  (i = n+i..o*n+r) */

  /* deltah[n+i] = y[o][i]*x'-y'[o-1][i] (i = 1..r) */
  norm = dpx*dpx;
  if (o > 0) {
    for (i = 1; i <= r; i++) {  
      aux         = dpy[o-1][q[i]];
      norm       += aux*aux;
      deltah[n+i] = h*(y[o][q[i]]*dpx-aux);
    }
    for (i = r+1; i <= n; i++)  {
      aux   = dpy[o-1][q[i]]; 
      norm += aux*aux;
    }
  }   

  /* deltah[(o-j)n+r+i] = y[j][i]*x'-y'[j-1][i] (j =o-1..1, i = 1..n) */
  for (j = o-1; j >= 1; j--) 
    for (i = 1; i <= n; i++) {
      aux                  = dpy[j-1][q[i]];
      norm                += aux*aux;
      deltah[(o-j)*n+r+i]  = h*(y[j][q[i]]*dpx-aux);
     } 
  for (i = 1; i <= r; i++)
    norm += dpy[o][q[i]]*dpy[o][q[i]];

  /* deltah[on+r+1] = (c'(s),c'(s)) - 1  */
  deltah[o*n+r+1] = h*(norm-1.0);
  
  return;
} 
/* fim SETH */



/* ********************************************************* */
/* Rotina para a construcao de DH(c0) onde c0 e o ponto      */
/* predito                                                   */
/* ********************************************************* */

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
)
{
  int  i,j,k,l;
  real cjaux;
  
  cjaux = cj*h;
      
  /* avaliacao de DF */
  DF(o,n,px,py,DFx,DFy); 

  /*  construcao de DH[i][j] (i = 1..n, j = 1..(o+1)n)  */
  for (i = 1; i <= r; i++)
    for (j = 1; j <= r; j++)
      DH[i][j] = DFy[o][q[i]][p[j]];
  for (i = r+1; i <= n; i++)
    for (j = 1; j <= r; j++)
      DH[i][j] = 0.0;
  for (k = o-1; k >= 0; k--)
    for (i = 1; i <= n; i++)
      for (j = 1; j <= n; j++)
        DH[i][(o-k-1)*n+r+j] = DFy[k][q[j]][p[i]];

  /*  DH[i][on+r+1] = DFx[i] (i = 1..n) */
  for (i = 1; i <= n; i++)
    DH[i][o*n+r+1] = DFx[p[i]];

  /*  DH[n+i][on+r+1] = cj y0[o][i] (i = 1..r) */
  for (i = 1; i <= r; i++)
    DH[n+i][o*n+r+1] = cjaux*py[o][q[i]];
  
  /*  DH[(o-k)n+r+i][on+r+1] = cj y0[k][i] (i = 1..r, k = o-1..1) */
  for (k = o-1; k >= 1; k--)
    for (i = 1; i <= n; i++)
      DH[(o-k)*n+r+i][o*n+r+1] = cjaux*py[k][q[i]];

  /*  construcao de DH[i][j] (i = n+r..on, j = r..(o+1)n) */
  for (l = o-1; l >= 1; l --)    
    for (k = o-1; k >= 0; k--)   
      for (i = 1; i <= n; i++)
        for (j = 1; j <= n; j++) 
          if (l == k) {
            if (i == j) 
              DH[(o-l)*n+r+i][(o-1-k)*n+r+j] = h*dpx;
            else  
              DH[(o-l)*n+r+i][(o-1-k)*n+r+j] = 0.0;
          } else if (l == k+1) {
            if (i == j) 
              DH[(o-l)*n+r+i][(o-1-k)*n+r+j] = -cjaux;
            else  
              DH[(o-l)*n+r+i][(o-1-k)*n+r+j] = 0.0;
          } else {
            DH[(o-l)*n+r+i][(o-1-k)*n+r+j] = 0.0;
          }
             
  /*  construcao de DH[i][j] (i = n+1..n+r, j = 1..o*n+r) */
  for (i = 1; i <= r; i++) { 
    for (j = 1; j <= r; j++)  
      if (i == j) {
        DH[n+i][j]   = h*dpx;
        DH[n+i][r+j] = -cjaux;
      } else {  
        DH[n+i][j]   = 0.0;
        DH[n+i][r+j] = 0.0;
      }
    for (j = n+r+1; j <= o*n+r; j++)  
      DH[i][j] = 0.0; 
  }    
    
  /*  construcao de DH[i][j] (i = n+1+r..(o-1)n, j = 1..r) */
  for (i = n+r+1; i <= o*n+r; i++)
    for (j = 1; j <= r; j ++)
      DH[i][j] = 0.0;

  /*  construcao de DH[(o+1)*n+1][i] (i = 1..on) */
  for (i = 1; i <= r; i++)
    DH[o*n+r+1][i] = 2.0*cjaux*dpy[o][i]; 

  for (k = o-1; k >= 0; k--)
    for (i = 1; i <= n; i++)
      DH[o*n+r+1][(o-k-1)*n+r+i] = 2.0*cjaux*dpy[k][q[i]]; 
    
  /*  construcao de DH[(o+1)*n+1][(o+1)*n+1] */
  DH[o*n+r+1][o*n+r+1] = 2.0*cjaux*dpx;

/*
  printf("\n\n\n\n");
  printf("cjaux = %lf\n",cjaux);
  printf("px = %lf, dpx = %lf\n",px,dpx);
  for (k = 0; k < o; k++)
    for (i = 1; i <= n; i++) 
      printf("py[%d][%d] = %lf, dpy[%d][%d] = %lf\n",k,i,py[k][i],k,i,dpy[k][i]);
  for (i = 1; i <= r; i++) 
    printf("py[%d][%d] = %lf, dpy[%d][%d] = %lf\n",o,i,py[k][i],o,i,dpy[k][i]);
  printf("\n\n");

  for (i = 1; i <= o*n+r+1; i++) {
    for (j = 1; j <= o*n+r+1; j++) {
      printf("%lf ",DH[i][j]);
    }
    printf("\n");
  }
  exit(0);
*/
      

  return;
} 
/* fim SETDH */



/* ********************************************************* */
/* Rotina para a construcao de uma aproximacao para DH(c0)   */
/* onde c0 e o ponto predito                                 */
/* ********************************************************* */

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
)
{
  int  i, j, k;  /* variaveis auxiliares   */
  real del;      /* armazenar o incremento */
  real save;     /* armazena o ponto       */
  real dsave;    /* armazena a derivada    */

  /* calculando a menor constante */ 
  uround   = sqrt(uround);

  /* calculando a derivada aproximada com relacao a x */

  /* calculando o incremento */
  del   = uround*MAX3(fabs(h*dx),fabs(x),fabs(wtx));
  del  *= FSIGN(h*dx);
  del   = (x+del)-x;

  /* salvando o ponto e a derivada */
  save  = x;
  dsave = dx;

  /* incrementando o ponto e a derivada */
  x    += del;
  dx   += del*cj;
  del   = 1.0/del;

  /* avaliando H */
  SETH(n,o,r,h,dx,dy,x,y,p,q,deltaaux,deltahaux,F);

  /* calculando a derivada aproximada */
  for (i = 1; i <= dim; i++)
    DH[i][dim] = (deltahaux[i]-deltah[i])*del;

  /* restaurando o ponto e a derivada */
  x     = save;
  dx    = dsave;

  /* calculando a derivada aproximada com relacao a y */
  for (k = o; k >= 0; k--)
    for (i = 1; i <= n; i++) {

      /* calculando o incremento */
      del       = uround*MAX3(fabs(h*dy[k][i]),fabs(y[k][i]),fabs(wty[k][i]));
      del      *= FSIGN(h*dy[k][i]);
      del       = (y[k][i]+del)-y[k][i];

      /* salvando o ponto e a derivada */
      save      = y[k][i];
      dsave     = dy[k][i];

      /* incrementando o ponto e a derivada */
      y[k][i]  += del;
      dy[k][i] += del*cj;
      del       = 1.0/del;

      /* avaliando H */
      SETH(n,o,r,h,dx,dy,x,y,p,q,deltaaux,deltahaux,F);

      /* calculando a derivada aproximada */
      for (j = 1; j <= dim; j++)
        DH[j][(o-k)*n+i] = (deltahaux[j]-deltah[j])*del;

      /* restaurando o ponto e a derivada */
      y[k][i]  = save;
      dy[k][i] = dsave;

    }

  return;
} 
/* fim SETDHAPPROX */



/* ******************************************************** */
/* rotina para realizar a predicao do ponto pc = (pcx,pcy), */
/* e de pdc = (pdcx,pdcy)  atraves do polinomio preditor    */    
/* ******************************************************** */

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
)
{   
  int i,j,l;

  /* calculo de pcx e de pdcx */
  *pcx = *pdcx = 0.0;
  for (i = 1; i <= k+1; i++) {
    (*pcx)  += phix[i];
    (*pdcx) += gama[i]*phix[i];
  }

  /* calculo de pcy e de pdcy */
  for (i = 0; i < o; i++) {   
    for (j = 1; j <= n; j++) {  
      pcy[i][q[j]] = pdcy[i][q[j]] = 0.0;
      for (l = 1; l <= k+1; l++) {  
         pcy[i][q[j]]  += phiy[l][i][q[j]];
         pdcy[i][q[j]] += gama[l]*phiy[l][i][q[j]];
      }
    }
  }  

  for (j = 1; j <= r; j++) {  
    pcy[o][q[j]] = pdcy[o][q[j]] = 0.0;
    for (l = 1; l <= k+1; l++) {  
       pcy[o][q[j]]  += phiy[l][o][q[j]];
       pdcy[o][q[j]] += gama[l]*phiy[l][o][q[j]];
    }
  }  

  /* calculo do vetor u utilizado no proc. NEWTON p/ correcao do ponto */
  for (j = 1; j <= r; j++) 
    u[j] = pcy[i][q[j]];
  for (i = o-1; i >= 0; i--) 
    for (j = 1; j <= n; j++) 
      u[(o-i)*n+j] = pcy[i][q[j]];
  u[o*n+r+1] = *pcx;

  return;
}  
/* fim predictor */
   


/* *********************************************************/
/* Esta rotina calcula os coeficientes para o calculo dos  */
/* polinomios preditor e corretor utilizando diferencas    */
/* divididas modificadas. Estes coeficientes serao         */
/* utilizados nas rotinas predictor, masterstep, update e  */
/* controlstep, e sao inicializados na rotina firststep.   */
/* *********************************************************/

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
)
{  
  int  i,j,l ;      /* controladores de loop     */
  int  kp1,kp2,km1; /*  valores de k+1, k+2, k-1 */ 
  real hx1,hx2;     /* variaveis auxiliares */
  real cjlast ;     /* var. p/ guardar o ultimo valor de cj  */ 
  real lambda;      /* raio de convergencia p/ a matriz jacobiana */  
  real alfa0;    
 
  kp1  = (*k)+1;
  kp2  = (*k)+2;
  km1  = (*k)-1;
  
  if ((h != (*hold)) || ((*k) != (*kold))) *ns =0;
  *ns = MIN2((*ns)+1,(*kold)+2);

  if (kp1 >= (*ns)) {

    /* inicializacao dos coef. alfa, beta, gama, sigma */
    alfa[1]  = 1.0;
    beta[1]  = 1.0;
    gama[1]  = 0.0;
    sigma[1] = 1.0;
    hx1      = h;   

    /*  UP DATE de psi, beta, alfa, sigma e gama */
    for (i = 2; i <= kp1; i++) { 
      hx2      = psi[i-1];   
      psi[i-1] = hx1;
      beta[i]  = beta[i-1]*psi[i-1]/hx2;
      hx1      = hx2+h;
      alfa[i]  = h/hx1;
      sigma[i] = (i-1)*sigma[i-1]*alfa[i];
      gama[i]  = gama[i-1]+alfa[i-1]/h;
    } 
  
    /* calculo de psi[k+1] */ 
    psi[kp1] = hx1;    
 
  } 
  /* fim if */

  /* calculo de alfas e alfa0 */
  (*alfas) = (alfa0) = 0.0;
  for (i = 1; i < kp1; i++) { 
    *alfas -= 1.0/(real)i ; 
    alfa0  -= alfa[i];
  }  

  /* calculo dos coeficientes principais */

  /* cj eh inicializado qdo k = 1  */
  cjlast = *cj;     
  *cj    = -(*alfas)/h; 

  /* calculo dos coeficientes p/ o erro com passo variavel */ 

  /* valor utilizado no teste de aceitacao do passo */
  *ck = alfa[kp1]+(*alfas)-alfa0;
  *ck = MAX2(fabs(*ck),alfa[kp1]); 

  /* teste p/ verificar se nova jacobiana eh necessaria  */

  /* calulo de lambda que estima o raio de convergencia */
  /* do metodo de Newton modificado                     */
  hx1    = 0.75/1.25;
  hx2    = 1.0/hx1;
  lambda = (*cj)/(*cjold);

  /* teste para o raio de convergencia */
  if ((lambda < hx1) || (lambda > hx2)) { 
    /* nova jacobiana eh necessaria */ 
    *aDH = 1;        
  } else {
    /* nova jacobiana nao eh necessaria */ 
    *aDH = 0;
  }
  
  /* verificando de e necessario reinicializar factor */
  if (*cj != cjlast) {
    /* inicializando factor */
    *factor = 100.0;
  }
 
  /* troca de phi por phi*,  phi = (phix,phiy) */
  if (kp1 >= (*ns)+1) {  
    /* troca de phix por phix estrela */
    for (i = (*ns)+1; i <= kp1; i++)
      phix[i] *= beta[i]; 
   
    /* troca de phiy por phiy estrela */
    for (l = (*ns)+1; l <= kp1; l++)  
      for (i = 0; i <= o; i++)
        for (j = 1; j <= n; j++) 
          phiy[l][i][j] *= beta[l];
  } 

  return;
}  
/* fim coefficient */



/* ********************************** */
/* UPDATE de phi, phi = (phix, phiy ) */
/* ********************************** */

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
)
{   
  int  i,j,l ;    /* controladores de loop  */
  
  /* calculando as novas diferencas divididas */

  /* se a ordem e menor que 5 */
  if (kold < 5) { 
   
    /* calculando as diferencas divididas de ordem k+2 */
    phix[kp2] = Ex;
    for (i = 0; i <= o; i ++)
      for (j = 1; j <= n; j ++)  
        phiy[kp2][i][q[j]] = Ey[i][q[j]];

  } 

  /* calculando as diferencas divididas de ordem k+1 */
  phix[kp1] += Ex ;
  for (i = 0; i <= o; i ++)
    for (j = 1; j <= n; j ++) 
      phiy[kp1][i][q[j]] += Ey[i][q[j]];
   
  /* calculando as diferencas divididas de ordem 1..k */
  for (l = 2; l <= kp1; l++)   
    phix[kp1-l+1] += phix[kp1-l+2];
  for (l = 2; l <= kp1; l++)   
    for (i = 0; i <= o; i++)
      for (j = 1; j <= n ; j++)
        phiy[kp1-l+1][i][q[j]] += phiy[kp1-l+2][i][q[j]];

  return;
} 
/* fim update */

 

/* ************************************************************* */
/* Esta rotina (rotina principal do GSDAE) e responsavel pela    */
/* integracao da EAD. Ela retorna o proximo ponto caso a         */
/* integracao seja bem sucedida ou a causa do erro no passo de   */
/* integracao.                                                   */
/* A integracao e feita utilizando o metodo BDF com coeficientes */
/* principais fixos. No calculo do ponto e utilizado um metodo   */
/* de Newton modificado.                                         */
/* As estrategias utilizadas nesta rotina sao semelhantes as da  */
/* rotina DASTP do codigo DASSL, so que adaptadas para este tipo */
/* de equacao.                                                   */
/* masterstep e uma funcao que retorna um valor inteiro que      */
/* significa:                                                    */
/*   1 : a rotina obteve sucesso no passo de integracao          */
/* -12 : se h < hmin                                             */
/* -13 : se o numero de condicao foi excedido                    */
/* -14 : se foram realizadas 20 tentativas sem exito ou devido   */
/*       as tolerancias ou ao numero de condicao da jacobiana    */
/* ************************************************************* */
 
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
)
{   
  real d;       /* armazena o valor de || cm+1 - cm ||                  */  
  real d0;      /* armazena o valor de || c1 - c0 ||, onde c0 e o valor */
                /* predito e c1 a primeira correcao obtida pelo metodo  */
	        /* de Newton modificado                                 */ 
  int  nint ;   /* numero de iteracoes do corretor2 caso haja neces-  */
                /*-sidade de reducao da amplitude do passo.             */
                /* observe que o maximo de iteracoes nesse caso eh 20   */
  int ncor;     /* numero de iteracoes para a resolucao do sistema de */
                /* equacoes atraves do metodo de NEWTON */
                /* o numero maximo de iteracoes eh 4.    */
  int  i,j,l;
  real tolerancia;
  real ro;
  real pnrm;      /* norma do ponto predito */
  real cond;
  real ac;      /* fator de aceleracao p/ o met de Newton */
  int  dim;

  /* definindo a tolerancia, a dimensao e o contador de correcoes */
  tolerancia = 2.3e-16;
  dim        = o*n+r+1;
  nint       = 0; 
   
  /* atualizando o vetor peso */
  weightvector(n,o,r,*cx,cy,q,atolx,rtolx,atoly,rtoly,wtx,wty);
  
  /* calculo dos coeficientes para o polinomio preditor e corretor */
  coefficient(n,o,r,k,*h,alfa,beta,gama,sigma,psi,alfas, 
              cj,cjold,factor,aDH,ck,phix,phiy,ns,hold,kold);

  /* calculo de pc = (pcx,pcy) e de pdc = (pdcx,pdcy) atraves do  */
  /* polinomio preditor                                           */
  predictor(n,o,r,*k,gama,phix,phiy,pcx,pcy,pdcx,pdcy,q,u);

  /* armazena os valores anteriores de c = (cx, cy) -- c(n)  */
  *cxx = *cx;
  *cx  = *pcx; 
  for (i = 0; i <= o; i++) 
    for (j = 1; j <= n; j++) {
      cyx[i][j] = cy[i][j];
      cy[i][j]  = pcy[i][j]; 
    }

  /* avalia a funcao no ponto predito */ 
  SETH(n,o,r,*h,*pdcx,pdcy,*pcx,pcy,p,q,deltax,deltahx,F); 
  (*naF) ++; 

  /* calcula a norma do ponto predito */
  pnrm = weightnorm(n,o,r,*pcx,pcy,q,*wtx,wty);
  
  /* copiando o valor de F no ponto predito */
  for (i = 1; i <= dim; i++) 
    deltah[i] = deltahx[i];

  /* avaliando a jacobiana se necessario */
  ncor = 0;
  if ((*aDH == 1) || ((*ifase == 0) && (*k == 1))) {  

    if (nDH == 0) {
      SETDHAPPROX(n,o,r,dim,*h,*cj,tolerancia,*pcx,*pdcx,*wtx,
                  pcy,pdcy,wty,p,q,delta,deltax,deltah,deltahx,DH,F); 
    } else {
      SETDH(n,o,r,*cj,*h,*pcx,pcy,*pdcx,pdcy,p,q,DFx,DFy,DF,DH); 
    }
    (*naDH) ++;
    *cjold  = *cj;
    *factor = 100.0;

    /* decomposicao QR de DH */
    QR(dim,dim,DH,Q,1,&cond);
    (*ndQR) ++;

    /* teste da condicao da jacobiana */
    if (cond > cdmax) { 
      /* ponto predito nao aceito */
      ncor = 4;
    }
       
    *aDH = 0;
  } 
  
  /* inicializando o vetor erro */
  *Ex = 0.0;
  for (i = 0; i <= o; i++)
    for (j = 1; j <= n; j ++)
      Ey[i][j] = 0.0;
                 
  /* laco para o calculo do ponto corrigido */
  do { 

    /* laco para a correcao pelo metodo de Newton modificado */ 
    /* numero maximo de iteracoes = 4                        */
    while (ncor < 4) {   
           
      ncor ++;
            
      /* fator de aceleracao do metodo de Newton modificado */
      ac = 2.0/(1.0+(*cj)/(*cjold));
  
      /* calculo de QRu = deltah, QR = DH, DH <- R */
      NEWTON(dim,Q,DH,u,deltah,ac) ; 

      /* calculo de v = cn(i+1) - cn(i) */ 
      for (i = o-1; i >= 0; i--)
        for (j = 1; j <= n; j++)  
          y[i][q[j]] = u[(o-i-1)*n+r+j];
      for (j = 1; j <= r; j++)  
        y[o][q[j]] = u[j];
      *x = u[o*n+r+1];

      /* atualizacao de c = (cx, cy0, cy1, ..., cyo), */
      /* sua derivada e o vetor de erro               */               
      *cx   -= (*x);
      *pdcx -= (*cj)*(*x);
      *Ex   -= (*x); 
      for (i = 0; i < o; i++) 
        for (j = 1; j <= n; j++) {
          cy[i][q[j]]   -= y[i][q[j]];
          pdcy[i][q[j]] -= (*cj)*y[i][q[j]];
          Ey[i][q[j]]   -= y[i][q[j]];
        }
      for (j = 1; j <= r; j++) {
        cy[o][q[j]]   -= y[o][q[j]];
        pdcy[o][q[j]] -= (*cj)*y[o][q[j]];
        Ey[o][q[j]]   -= y[o][q[j]];
      }

      /* calculo da norma peso no ponto corrigido */
      d = weightnorm(n,o,r,*x,y,q,*wtx,wty);
        
      /* verificando se o ponto corrigido deve ser aceito */ 
      if (d <= (tolerancia*100.0*pnrm)) { 

        /* se ftol[1] != 0 avalie a funcao */
        if (ftol[1] != 0.0) {

          /* Avaliando a funcao F no ponto aceito */
          SETH(n,o,r,*h,*pdcx,pdcy,*cx,cy,p,q,delta,deltah,F);
          (*naF) ++;

          /* calcula a norma peso de F */
          if (functionnorm(n,delta,ftol)) {
            /* ponto corrigido aceito */
            return (1);
          }

        } else {
          /* ponto corrigido aceito */
          return (1);
        }

      } else if (ncor == 1) {  

        d0 = d;

        if ((*factor)*d0 <= 0.33) {     
  
          /* se ftol[1] != 0 avalie a funcao */
          if (ftol[1] != 0.0) {
  
            /* Avaliando a funcao F no ponto aceito */
            SETH(n,o,r,*h,*pdcx,pdcy,*cx,cy,p,q,delta,deltah,F);
            (*naF) ++;
  
            /* calcula a norma peso de F */
            if (functionnorm(n,delta,ftol)) {
              /* ponto corrigido aceito */
              return (1);
            }
  
          } else {
            /* ponto corrigido aceito */
            return (1);
          }

        } else { 

          /* ponto corrigido nao aceito */
          SETH(n,o,r,*h,*pdcx,pdcy,*cx,cy,p,q,delta,deltah,F);
          (*naF) ++;

        }

      } else {    

        /* estima o raio de convergencia para iteracao de NEWTON */
        ro      = d/d0 ;    
        ro      = ROOT(ro,ncor-1); 
        *factor = ro/(1.0-ro);  

        if (ro > 0.9) {  

          /* ponto corrigido nao aceito                       */
          /* necessario reduzir o passo ou avaliar a derivada */
          /* no ponto predito                                 */
          ncor = 4;    

        } else if ((*factor)*d <= 0.33) {  

          /* se ftol[1] != 0 avalie a funcao */
          if (ftol[1] != 0.0) {
  
            /* Avaliando a funcao F no ponto aceito */
            SETH(n,o,r,*h,*pdcx,pdcy,*cx,cy,p,q,delta,deltah,F);
            (*naF) ++;
  
            /* calcula a norma peso de F */
            if (functionnorm(n,delta,ftol)) {
              /* ponto corrigido aceito */
              return (1);
            }

          } else {
            /* ponto corrigido aceito */
            return (1);
          }

        } else {

          /* ponto corrigido nao aceito */
          SETH (n,o,r,*h,*pdcx,pdcy,*cx,cy,p,q,delta,deltah,F);
          (*naF) ++;  

        }
  
      } 

    } 
    /* fim while */

    /* verificando o motivo da nao convergencia      */
    if (*aDH == 1) {    

      /* A matriz jacobiana nao foi avaliada no ponto predito    */
      /* Avalia a jacobiana e tenta-se novamente                 */

      /* calcula a norma do ponto predito */
      pnrm = weightnorm(n,o,r,*pcx,pcy,q,*wtx,wty);

      /* copiando o valor de F no ponto predito */
      for (i = 1; i <= dim; i++) 
        deltah[i] = deltahx[i];
                   
      /* calcula a jacobina no ponto predito */
      if (nDH == 0) {
        SETDHAPPROX(n,o,r,dim,*h,*cj,tolerancia,*pcx,*pdcx,*wtx,
                    pcy,pdcy,wty,p,q,delta,deltax,deltah,deltahx,DH,F);
      } else {
        SETDH(n,o,r,*cj,*h,*pcx,pcy,*pdcx,pdcy,p,q,DFx,DFy,DF,DH); 
      }
      (*naDH) ++; 
      *cjold  = *cj;
      *factor = 100.0;

      ncor = 0;
      /* decomposicao QR de DH */
      QR(dim,dim,DH,Q,1,&cond);
      (*ndQR) ++;
      if (cond > cdmax) { 
        /* matriz mal condicionada */
        ncor = 4;
        (*nff)++;
      } 
      *aDH = 0;   

      /* inicializacao do vetor erro */
      *cjold = *cj ;
      *Ex    = 0.0;
      for (i = 0; i <= o; i++)
        for (j = 1; j <= n; j ++)
           Ey[i][j] = 0.0;
                
    } else  {  

      /* O metodo de Newton nao convergiu. O ponto anterior e  */
      /* restaurado e uma nova tentativa com a amplitude do    */
      /* reduzida por um fator de um quarto e feita            */

      /* incrementa o contador de falhas */
      *nNwf += 1;     
      *ifase = 1;

      /* restaura psi e phi */
      for (l = (*ns)+1; l <= (*k)+1; l++) {  
        phix[l] = phix[l]/beta[l];
        for (i = 0; i <= o; i++)
          for (j = 1; j <= n; j++)
            phiy[l][i][q[j]] = phiy[l][i][q[j]]/beta[l];
      }
      for (i = 1; i <= *k; i++)  
	psi[i] = psi[i+1]-(*h);

      /* reduz o passo */
      (*h) *= 0.25; 

      if (fabs(*h) < hmin) { 
        /* passo muito pequeno */
        *cx = *cxx;
        for (i = 0; i <= o; i++)
          for (j = 1; j <= n; j++)
            cy[i][q[j]] = cyx[i][q[j]];

        return (-12); 
      }
         
      /* calculo dos coeficientes para o nova tentativa */
      coefficient(n,o,r,k,*h,alfa,beta,gama,sigma,psi,alfas, 
                  cj,cjold,factor,aDH,ck,phix,phiy,ns,hold,kold);

      /* calculo do ponto predito */
      predictor(n,o,r,*k,gama,phix,phiy,pcx,pcy,pdcx,pdcy,q,u);

      /* avalia a funcao no ponto predito  */
      SETH(n,o,r,*h,*pdcx,pdcy,*pcx,pcy,p,q,deltax,deltahx,F);          
      (*naF) ++; 

      /* calcula a norma do ponto predito */
      pnrm = weightnorm(n,o,r,*pcx,pcy,q,*wtx,wty);

      /* copiando o valor de F no ponto predito */
      for (i = 1; i <= dim; i++) 
        deltah[i] = deltahx[i];

      /* avaliacao da jacobiana DH */
      if (nDH == 0) {
        SETDHAPPROX(n,o,r,dim,*h,*cj,tolerancia,*pcx,*pdcx,*wtx,
                    pcy,pdcy,wty,p,q,delta,deltax,deltah,deltahx,DH,F);
      } else {
        SETDH(n,o,r,*cj,*h,*pcx,pcy,*pdcx,pdcy,p,q,DFx,DFy,DF,DH);
      }
      (*naDH)++;
      *aDH    = 0;
      *cjold  = *cj;
      *factor = 100.0;

      ncor = 0;
      /* decomposicao QR de DH */
      QR(dim,dim,DH,Q,1,&cond);
      (*ndQR) ++;
      if (cond > cdmax) {
        /* matriz mal condicionada */
        ncor = 4;
        (*nff)++;
      }

      /* inicializacao do erro e primeiro ponto */
      *cx = *pcx; 
      *Ex = 0.0;
      for (i = 0; i<= o; i++) {
        for (j = 1; j <= n; j++) {
          cy[i][q[j]] = pcy[i][q[j]]; 
          Ey[i][q[j]] = 0.0;
        }
      } 
 
    } 

    /* incrementando o contador de tentativas de correcao */
    nint++;  
   
  } while (nint <= 20); 
  /* fim do-while */ 

  /* falha na correcao - copiando o ponto anterior */
  *cx = *cxx;
  for (i = 0; i <= o; i++)
    for (j = 1; j <= n; j++)
      cy[i][q[j]] = cyx[i][q[j]];

  /* nao convergiu : excedeu numero maximo de iteracoes */
  if (cond > cdmax)
    return (-13);
  else  
    return (-14);

} 
/* fim masterstep */



/*************************************************************** */
/* Rotina para a aceitacao ou rejeicao do passo e criterio para  */
/* aumentar ou diminuir a ordem do metodo e a amplitude do passo */
/*************************************************************** */

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
)
{  
  int  i,j,l;                  
  int  knew,kdiff,kp1,kp2,km1;
  real terk,terkm1,terkm2,terkp1;
  real erk,erkm1,erkm2,erkp1;
  real enorm,err,est;
  real factor;
 
  /* armazena a ordem anterior */
  knew  = *k;
  kp1   = (*k)+1;
  kp2   = (*k)+2;
  km1   = (*k)-1;
  kdiff = (*k)-(*kold);
    
  /* construcao dos termos para o teste de estabilidade */

  /* calculo do termo de ordem k */
  enorm = weightnorm(n,o,r,Ex,Ey,q,wtx,wty);
  erk   = enorm*sigma[kp1]; 
  terk  = erk*kp1; 
  est   = erk;
 
  if (*k > 1) {  

    /* calculo do termo de ordem k-1 */
    *x = phix[kp1] + Ex;
    for (i = 0; i <= o; i++)
      for (j = 1; j <= n; j++)
        y[i][j] = phiy[kp1][i][j] + Ey[i][j] ;
    erkm1   = weightnorm(n,o,r,*x,y,q,wtx,wty);
    erkm1  *= sigma[*k];
    terkm1  = erkm1*(*k);

    if (*k > 2) { 

      /* calculo do termo de ordem k-2 */
      *x += phix[*k];
      for (i = 0; i <= o; i++)
        for (j = 1; j <= n; j++)
          y[i][j] += phiy[*k][i][j];
      erkm2   = weightnorm(n,o,r,*x,y,q,wtx,wty);  
      erkm2  *= sigma[km1];
      terkm2  = erkm2*km1;

      if (MAX2(terkm1,terkm2) <= terk) {

        /* calculo de est e knew */
        knew = (*k)-1;
        est  = erkm1;

      }

    } else if (terkm1 <= 0.5*terk) { 

      /* calculo de est e knew */
      knew = (*k)-1;
      est  = erkm1;

    }

  }

  err = ck*enorm;  

  /* teste para a aceitacao do passo */
  if (err > 1.0) {  

    /* passo nao foi aceito */
    *ifase = 1;
    (*cfalhas)++;
    (*nflhs)++;

    /* restaura phi */
    for (l = ns+1; l <= kp1; l++) {  
      phix[l] /= beta[l];
      for (i = 0; i <= o; i++)
        for (j = 1; j <= n; j++)
          phiy[l][i][j] /= beta[l];
    }

    /* restaura psi */  
    for (i = 2; i <= kp1; i++)  
      psi[i-1] = psi[i]-(*h) ; 
       
    /* definindo a nova ordem do metodo */
    *k = knew;

    if (*cfalhas == 1) {

      /* E a primeira falha na convergencia  */
      /* calculo do fator de escala de h     */
      factor  = 1.0/(2.0*est+0.0001);
      factor  = 0.9*ROOT(factor,(*k)+1);
      factor  = MAX2(0.25,MIN2(0.9,factor));
      (*h)   *= factor;

    } else if (*cfalhas == 2) {   

      /* E a segunda falha consecutiva na convergencia   */
      /* o passo eh reduzido por um fator de um quarto */
      (*h) *= 0.25;   
               
    } else { 

      /* Houve mais de duas falhas consecutivas */
      /* a ordem eh reduzida a unidade */
      /* a amplitude do passo eh reduzida por um fator de um quarto */ 
      (*k)  = 1;      
      (*h) *= 0.25;  

    }

    /* verificando a amplitude do passo */ 
    if (fabs(*h) < hmin) 
      return (-9);  

    return (0);

  } else  { 

    /* o passo foi aceito */
    *cfalhas = 0;    
    *kold    = *k;
    *hold    = *h;
    *sold    = *s;
    *s      += *h;

    /* teste para o fim da fase inicial */
    if ((knew == km1) || (*k == 5)) 
      *ifase = 1;

    if (*ifase == 0) {

      /* fase inicial */
      /* aumente a ordem e a amplitude do passo */
      (*k)++;
      (*h) *= 2.0;

    } else {

      if (knew == km1) { 

        /* diminua a ordem */
        est  = erkm1;
        (*k)--;

      } else if ((*k < 5) && (kp1 < ns) && (kdiff != 1))  {

        /* calculo do termo de ordem k+1 */
        *x = Ex - phix[kp2];
        for (i = 0; i <= o; i++)
          for (j = 1; j <= n; j++)
            y[i][j] = Ey[i][j] - phiy[kp2][i][j];
        terkp1 = weightnorm(n,o,r,*x,y,q,wtx,wty); 
        erkp1  = terkp1/kp2;
 
        /* teste para aumentar ou diminuir a ordem */
        if (*k == 1) {  
 
          if (terkp1 < 0.5*terk) {  
 
            /* aumente a ordem */
            est  = erkp1;
            (*k)++ ;
 
          }
 
        } else if (terkm1 < MIN2(terk,terkp1)) {
  
           /* diminua a ordem */
           est  = erkm1;
           (*k)--;     
  
        } else if (terkp1 < terk) { 
 
          /* aumente a ordem */
          est  = erkp1;          
          (*k)++;     
  
        } 

      }

      /* calculo do fator de escala de h     */
      factor = 1.0/(2.0*est+0.0001);
      factor = ROOT(factor,(*k)+1);
 
      if (factor >= 2.0) { 

        /* dobra a amplitude do passo     */ 
        (*h) *= 2.0 ;
        if ((hmax != 0.0) && (fabs(*h) > hmax)) 
          *h = SIGN(*h)*hmax;

      } else if (factor <= 1.0) { 

        /* o passo e reduzido pelo fator factor */
        factor  = MAX2(0.5,MIN2(0.9,factor));
        (*h)   *= factor;                 

      }

    }

    /* realizando o update das diferencas divididas */
    update(n,o,*kold,kp1,kp2,Ex,Ey,q,phix,phiy);

    return (1);  
  }

} 
/* fim controlstep */        



/***************************************************/
/* Avaliacao do polinomio interpolador nos pontos  */
/* c(s(n)),c(s(n-1)),..,c(s(n-kold)) no ponto      */
/* s(n)-h.                                          */
/***************************************************/

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
)
{  
  int  i,j,l;
  real c,d,gamma;

  /* inicializando com o ponto c(s(n)) e a derivada nula */
  *xint  = phix[1];
  *dxint = 0.0;
  for (i = 0; i <= o-1; i++)   
    for (j = 1; j <= n; j++) {
      yint[i][j]  = phiy[1][i][j];
      dyint[i][j] = 0.0;
    }
  for (j = 1; j <= r; j++) { 
    yint[o][j]  = phiy[1][o][j];
    dyint[o][j] = 0.0;
  }
 
  /* caluculo utilizando as diferencas divididas modificadas */
  c     = 1.0;
  d     = 0.0;
  gamma = hint/psi[1];
  for (l = 2; l <= kold+1; l++) { 
    d       = d*gamma+c/psi[l-1];
    c      *= gamma;
    gamma   = (hint+psi[l-1])/psi[l];
    *xint  += c*phix[l];
    *dxint += d*phix[l];
    for (i = 0; i <= o-1; i++)  
      for (j = 1; j <= n; j++) {   
        yint[i][j]  += c*phiy[l][i][j];
        dyint[i][j] += c*phiy[l][i][j];
      }
    for (j = 1; j <= r; j++) {  
      yint[o][j]  += c*phiy[l][o][j];
      dyint[o][j] += c*phiy[l][o][j];
    }
  }

  /* aproximacao da variavel desconhecida */
  if (o > 0) 
    for (i = r+1; i <= n; i++)
      yint[o][i] = dyint[o-1][i]/(*dxint);

  return;
} 
/* fim interpolator */



/***********************************************/
/* Avaliando o polinomio interpolador em send  */
/***********************************************/

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
)
{  
  real h;

  h = (*send)-s;

  /* avaliando o polinomio interpolador em send */
  interpolator(n,o,r,h,x,y,dx,dy,p,q,kold,phix,phiy,psi);

  return;
} 
/* fim interpolators */
  


/*************************************************/
/* Avaliando o polinomio interpolador em xend.   */
/* O calculo de xend e feito utilizando o metodo */
/* de Newton com tolerancia atolx                */
/*************************************************/

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
)
{  
  int  cont;
  real hend, xout;

  /* inicializando o xend */
  cont   = 0;
  xout   = *xend;
  *send  = s;
  *xend  = x;
  *dxend = pdcx;

  /* encontra send */
  do { 

    /* aproximacao de send */
    *send -= ((*xend)-xout)/(*dxend);
    hend   = (*send)-s;
    /* calculo do polinomio interpolador em send */
    interpolator(n,o,r,hend,xend,yend,dxend,dyend,p,q,kold,phix,phiy,psi);
    cont ++;

  } while ((fabs(xout-(*xend)) > tol*fabs(xout)) && 
           (cont < 500) && 
           (fabs(hend) > tol*fabs(s))); 

  /* verificaco da convergencia do metodo */
  if (cont == 500) {
    /*
    printf("\n\nInterpoladorx : nao convergiu\n\n");
    */
  }
  
  return; 
} 
/* fim interpolatorx */
  


/*******************************************************/
/* Avaliando o polinomio interpolador na singularidade */
/* transversal. O calculo de s (c(s) : ponto singular) */
/* e feito utilizando o metodo Regula-Falsi a partir   */
/* de x' nos dois ultimos pontos.                      */
/*******************************************************/

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
)
{  
  int  cont;
  real hint,dx;
  real sa,sb,xa,xb,pa,pb;

  /* xa : x avaliado em sa  */
  /* xb : x avaliado em sb  */
  /* pa : x' avaliado em sa */
  /* pb : x' avaliado em sb */
  sa     = sold;
  sb     = s;
  xa     = xold;
  xb     = x;
  pa     = pdcxold;
  pb     = pdcx;

  /* xint : x avaliado em sint  */
  /* dx   : x' avaliado em sint */
  *sint  = sb;
  *xint  = xb;
  dx     = pb;
  cont   = 0;

  if (kold > 1) {

    while ((cont < 500) && 
           (fabs(sb-sa) > tol*fabs(sb)) && 
           (fabs(dx) > tol)) { 

      /* calculo de sint pelo metodo Regula-Falsi */
      *sint = sa - (sb-sa)*pa/(pb-pa);
      hint  = (*sint)-s;
  
      /* avaliando o polinomio interpolador em sint */
      interpolator(n,o,r,hint,xint,yint,&dx,dyint,p,q,kold,phix,phiy,psi);

      /* verificando os estremos do intervalo */
      if (dx*pb < 0.0) {
        sa = *sint;
        xa = *xint;
        pa = dx;
      } else {
        sb = *sint;
        xb = *xint;
        pb = dx;
      }

      cont++;
    }
    /* fim do-while */

  } else {

  }

  /* verificacao da convergencia do metodo */
  if (cont == 500) {
    /*
    printf("\n\nInterpolatorsing : nao convergiu em 500 passos\n\n");
    */
  }

  return; 
} 
/* fim interpolatorsing */



/**************************************************************/
/* Calculo de B = (DFx+DFy[0]y[1]+...+DFy[o-1]y[o]  DFy[o])^t */
/**************************************************************/

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
)
{
  int  i,j,k;  /* variaveis auxiliares */

  /* B[1] = DFx + DFy[0] y[1] + .. + DFy[o-1] y[o] */
  for (i = 1; i <= n; i++) {
    B[1][i] = DFx[p[i]];
    for (k = 0; k <= o-2; k++) 
      for (j = 1; j <= n; j++) 
        B[1][i] += DFy[k][q[j]][p[i]]*y[k+1][q[j]];
    if (o > 0)  
      for (j = 1; j <= r; j++) 
        B[1][i] += DFy[o-1][q[j]][p[i]]*y[o][q[j]];
  }

  /* B[i+1][j] = DFy[o][i][j] (i = 1..n, j = 1..r) */
  for (i = 1; i <= r; i++) 
    for (j = 1; j <= r; j++)
      B[j+1][i] = DFy[o][q[j]][p[i]];
  for (i = r+1; i <= n; i++)
    for (j = r+1; j <= n; j++)
      B[j+1][i] = 0.0;

  /* B[i+1][j] = DFy[o-1][i][j] (i = 1..n, j = r+1..n) */
  if (o > 0)  
    for (i = 1; i <= n; i++) 
      for (j = r+1; j <= n; j++)
        B[j+1][i] = DFy[o-1][q[j]][p[i]];

  return;
} 
/* fim SETB */



/* ************************************************************* */
/* Esta rotina e responsavel pelo calculo da tangente no ponto   */
/* para inicializacao dos parametros para o primero passo do     */
/* metodo.                                                       */
/* Ela verifica se a condicao inicial satisfaz a EAD, calcula o  */
/* posto de DFy[o] e caso este seja menor que n calcula as       */
/* permutacoes de variaveis e linhas da funcao que define a EAD. */
/* Ela retorna um valor inteiro para controle do passo de        */
/* integracao e do usuario. Este inteiro significa :             */
/*   0 : ponto regular                                           */
/*   1 : singularidade tansversal                                */
/*   2 : diminuiu o posto de DFy[[o]                             */
/*   3 : singularidade transversal com posto de DFy[o] inferior  */
/*   4 : diminuiu a ordem da EAD                                 */
/*   5 : singularidade transversal com ordem inferior            */
/*  -3 : ponto inadequado - nao satisfaz a EAD com a tolerancia  */
/*       desejada                                                */
/*  -4 : posto maior que o informado                             */
/*  -5 : ordem zero e posto zero                                 */
/*  -6 : o posto de DFy[o] varia em uma vizinhaca do ponto       */
/*  -7 : ordem menor e o posto de DFy[o] varia em uma vizinhaca  */
/*       do ponto                                                */
/*  -8 : singularidade transversal                               */
/*  -9 : singularidade transversal com posto de DFy[o] inferior  */
/*  -10: singularidade transversal com ordem inferior            */
/* ************************************************************* */

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
)
{  
  real  cond;  /* condicao da matriz B              */
  real  norm;  /* norma de tau                      */
  int   raux1; /* posto de DFy[o]                   */
  int   raux2; /* posto de DFy[o] em uma vizinhanca */
  int   oaux;  /* ordem da EAD                      */
  int   i,j;   /* variaveis auxiliares              */

  /* tau = (taux,tauy) define a tangente no ponto inicial */
  /* tau pertence ao nucleo de  B^t = (DF(c),w(c)c')^t    */
  /* Q : decomposicao QR de B^t                           */


  /* avalida a funcao que define a EAD */
  F(*o,n,cx,cy,delta);
  (*naF)++;

  /* se foi definida tolerancias para a avaliacao da EAD */
  if (ftol[1] != 0.0) {

    /* verifica se o ponto satisfaz a EAD */
    if (!functionnorm(n,delta,ftol)) {

      /* ponto nao satisfaz a EAD */
      return (-3);

    } 

  } else {

    /* nao foi definida tolerancias para a avaliacao da EAD */
    /* definindo tolerancias iguais a rtolx                 */
    for (i = 1; i <= n; i++) 
      if (fabs(delta[i]) > rtolx) {

        /* ponto nao satisfaz a EAD */
        return (-3);

      }
  }

  /* verificando se o posto, a ordem e as permutacoes foram */
  /* predefinidos pelo usuario                              */
  if (!irank) {

    /* o posto nao foi definido pelo usuario               */
    /* definindo o posto, a ordem e as permutacoes para    */
    /* a inicializacao                                     */

    /* avaliando a jacobiana */
    if (nDH == 1) {

      /* jacobiana exata */
      DF(*o,n,cx,cy,DFx,DFy);

    } else {
    
      /* jacobiana aproximada */
      DFAPPROX(n,*o,dir,cx,cy,delta,deltaaux,DFx,DFy,F); 

    }
    (*naDH)++;
 
    /* verificando o posto de DFy[o] */
    /* B = DFy[o]                    */
    for (i = 1; i <= n; i++)
      for (j = 1; j <= n; j++)
        B[i][j] = DFy[*o][j][i];

    /* decomposicao QR de DFy[o]                            */
    /* r define o posto de DFy[o]                           */
    /* p a permutacao de linhas e q a permutacao de colunas */
    raux1 = QR2(n,B,Q,p,q);
    (nQR)++;

    /* copiando a ordem */
    oaux = *o;

    /* verificando se o posto de DFy[o] e diferennte do informado */
    if (raux1 > *r) {
  
      /* posto de DFy[o] maior que o informado */
      /* retornando erro na entrada de dados   */
      return (-4);

    } else if (raux1 < *r) {
    
      /* posto de DFy[o] menor que o informado */

      /* verifincando se o posto e constante em uma vizinhanca */
      raux2 =  rankneighbourhood(n,*o,raux1,dir,*o,cx,cy,delta,deltaaux, 
                                 DFx,DFy,Q,B,p,q,nDH,naDH,nQR,F,DF);

      /* reavaliando a jacobiana */
      if (nDH == 1) {

        /* jacobiana exata */
        DF(*o,n,cx,cy,DFx,DFy);

      } else {
    
        /* jacobiana aproximada */
        DFAPPROX(n,*o,dir,cx,cy,delta,deltaaux,DFx,DFy,F); 

      }
      (*naDH)++;

      /* verificando se o posto se manteve na vizinhanca */
      if (raux1 == raux2) {
  
        /* o posto se manteve na vizinhanca e e menor que o mensionado */

        /* verificando se a ordem se manteve */
        if (raux1 == 0) {

          /* a ordem diminuiu */

          /* verificando se a equacao e algebrica */
          if (*o == 0) {

            /* posto e ordem zero                  */
            /* retornando erro na entrada de dados */
            return (-5);

          }

          /* copiando a nova ordem */
          oaux = (*o)-1;

          /* definindo o posto sendo igual a dimensao */
          *r = n;

          /* verificando o posto de DFy[ord] */
          /* B = DFy[ord]                    */
          for (i = 1; i <= n; i++)
            for (j = 1; j <= n; j++)
              B[i][j] = DFy[oaux][j][i];

          /* decomposicao QR de DFy[oaux] */
          raux1 = QR2(n,B,Q,p,q);
          (nQR)++;

          /* verificando o posto de DFy[oaux] */
          if (raux1 < n-1) {

            /* verifincando se o posto e constante em uma vizinhanca */
            raux2 =  rankneighbourhood(n,*o,raux1,dir,*o,cx,cy,delta,
				       deltaaux,DFx,DFy,Q,B,p,q,nDH,
				       naDH,nQR,F,DF);

            /* reavaliando a jacobiana */
            if (nDH == 1) {

              /* jacobiana exata */
              DF(*o,n,cx,cy,DFx,DFy);

            } else {
    
              /* jacobiana aproximada */
              DFAPPROX(n,*o,dir,cx,cy,delta,deltaaux,DFx,DFy,F); 

            }
            (*naDH)++;

            /* verificando se o posto se manteve na vizinhanca */
            if ((raux1 != raux2) && (*r != raux2)) {

              /* a ordem diminuiu e o posto de DFy[o] varia na vizinhanca */
              /* singularidade nao transversal                            */
              return (-7);

            }

          }

          /* copiando o posto */
          raux1 = n;

        } 


      } else {

        /* possibilidade se singularidade transversal */

        if (raux2 != *r) {

	   /* posto de DFy[o] varia na vizinhanca */
           /* singularidade nao transversal       */
	   return (-6);

        }

        /* redefinindo o posto */ 
        raux1 = *r;

      }

    }

  } else {

    /* o posto a ordem e as permutacoes foram predefinidos */
    /* pelo usuario                                        */

    /* avaliando a jacobiana */
    if (nDH == 1) {

      /* jacobiana exata */
      DF(*o,n,cx,cy,DFx,DFy);

    } else {
    
      /* jacobiana aproximada */
      DFAPPROX(n,*o,dir,cx,cy,delta,deltaaux,DFx,DFy,F); 

    }
    (*naDH)++;

    /* definindo a ordem e o posto sendo iguais aos informados */
    oaux  = *o;
    raux1 = *r;

  }

  /* calcula a matriz que define a tangente no ponto inicial */
  SETB(n,oaux,raux1,cy,p,q,DFx,DFy,B);

  /* decomposicao QR de B */
  QR(n+1,n,B,Q,0,&cond);
  (nQR)++;

  /* testando se o posto de B e maximo */ 
  if (cond > cdmax) {

    /* o posto de B nao e maximo     */
    /* singularidade nao transversal */

    if (oaux != *o) {

      /* a ordem nao se manteve */

      /* copiando a ordem e o posto */
      *o = oaux;
      *r = raux1;

      return (-10);

    }

    if (raux1 != *r) {

      /* o posto nao se manteve */

      /* copiando a ordem e o posto */
      *o = oaux;
      *r = raux1;

      return (-9);

    } 

    return (-8);

  }

  /* o posto de B e maximo    */
  /* calculando o nucleo de B */

  /* calcula taux e tauy[o][i] (i = 1..r) */
  *taux = Q[n+1][1];
  norm = (*taux)*(*taux);
  for (i = 1; i <= raux1; i++) { 
    tauy[oaux][q[i]]  = Q[n+1][i+1];
    norm           += tauy[oaux][q[i]]*tauy[oaux][q[i]];
  }

  /* tauy[o-1][i] (i = r+1..n), tauy[o-1][i] (i = 1..r) */
  if (oaux > 0) {
    for (i = raux1+1; i <= n; i++) { 
      tauy[oaux-1][q[i]]  = Q[n+1][i+1];
      norm             += tauy[oaux-1][q[i]]*tauy[oaux-1][q[i]];
    }
    for (i = 1; i <= raux1; i++) {
      tauy[oaux-1][q[i]]  = cy[oaux][q[i]]*(*taux);
      norm             += tauy[oaux-1][q[i]]*tauy[oaux-1][q[i]];
    }
  }

  /* calcula tauy[j] j = 0,..,o-2 */
  for (j = 0; j <= oaux-2; j++)
    for (i = 1; i <= n; i++) {
      tauy[j][q[i]]  = (*taux)*cy[j+1][q[i]];
      norm          += tauy[j][q[i]]*tauy[j][q[i]];
    }
  
  /* norma de tau */
  norm = sqrt(norm);

  /* calculo de tau normalizado */ 
  *taux /= norm;
  for (j = 0; j <= oaux; j++)
    for (i = 1; i <= n; i++) 
      tauy[j][i] /= norm;


  /* retornando sucesso */
  if (fabs(*taux) < atolx) {

    /* singularidade transversal */

    if (oaux != *o) {

      /* a ordem nao se manteve */

      /* copiando a ordem e o posto */
      *o = oaux;
      *r = raux1;

      return (5);

    }

    if (raux1 != *r) {

      /* o posto nao se manteve */

      /* copiando a ordem e o posto */
      *o = oaux;
      *r = raux1;

      return (3);

    } 

    return (1);

  }

  /* ponto regular */

  if (oaux != *o) {

    /* a ordem nao se manteve */

    /* copiando a ordem e o posto */
    *o = oaux;
    *r = raux1;

    return (4);

  }

  if (raux1 != *r) {

    /* o posto nao se manteve */

    /* copiando a ordem e o posto */
    *o = oaux;
    *r = raux1;

    return (2);

  } 

  return (0); 
} 
/* fim settau */



/****************************************************************/
/* rotina que calucula o posto de DFy[ord] em uma vizinhanca de */
/* um ponto                                                     */
/****************************************************************/

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
)
{
  int  i,j,k,l; /* variaveis auxiliares */
  int  raux;    /* posto de DFy[ord]    */
  real save;    /* variavel para armazenar o valor original */ 

  /* salvando e perturbando a variavel x */
  save  = cx;
  cx   += 1.0e-6;

  /* avaliando a jacobiana */
  if (nDH == 1) {

    /* jacobiana exata */
    par->DF(o,n,cx,cy,DFx,DFy);

  } else {
    
    /* jacobiana aproximada */
    DFAPPROX(n,o,dir,cx,cy,delta,deltaaux,DFx,DFy,F); 

  }
  (*naDH)++;
 
  /* verificando o posto de DFy[ord] */
  /* B = DFy[ord]                    */
  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++)
      B[i][j] = DFy[ord][j][i];

  /* decomposicao QR de DFy[ord] */
  raux = QR2(n,B,Q,p,q);
  (nQR)++;

  /* voltando a variavel x ao valor original */
  cx = save;

  /* verificando se houve mudanca no posto de DFy[ord] */
  if (r != raux) {
   
    /* retornando o novo posto */
    return (raux);

  }

  /* selecionando a variavel y */
  for (k = 0; k <= o; k++) 
    for (l = 1; l <= n; l++) {

      /* salvando e perturbando a variavel y[k][l] */
      save      = cy[k][l];
      cy[k][l] += 1.0e-6;

      /* avaliando a jacobiana */
      if (nDH == 1) {

        /* jacobiana exata */
        par->DF(o,n,cx,cy,DFx,DFy);

      } else {
    
        /* jacobiana aproximada */
        DFAPPROX(n,o,dir,cx,cy,delta,deltaaux,DFx,DFy,F); 

      }
      (*naDH)++;
 
      /* verificando o posto de DFy[ord] */
      /* B = DFy[ord]                    */
      for (i = 1; i <= n; i++)
        for (j = 1; j <= n; j++)
          B[i][j] = DFy[ord][j][i];

      /* decomposicao QR de DFy[ord] */
      raux = QR2(n,B,Q,p,q);
      (nQR)++;
    
      /* voltando a variavel y[k][l] ao valor original */
      cy[k][l] = save;

      /* verificando se houve mudanca no posto de DFy[ord] */
      if (r != raux) {
   
        /* retornando o novo posto */
        return (raux);

      }

    }


  /* salvando e perturbando a variavel x */
  save  = cx;
  cx   -= 1.0e-6;

  /* avaliando a jacobiana */
  if (nDH == 1) {

    /* jacobiana exata */
    par->DF(o,n,cx,cy,DFx,DFy);

  } else {
    
    /* jacobiana aproximada */
    DFAPPROX(n,o,dir,cx,cy,delta,deltaaux,DFx,DFy,F); 

  }
  (*naDH)++;
 
  /* verificando o posto de DFy[ord] */
  /* B = DFy[ord]                    */
  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++)
      B[i][j] = DFy[ord][j][i];

  /* decomposicao QR de DFy[ord] */
  raux = QR2(n,B,Q,p,q);
  (nQR)++;

  /* voltando a variavel x ao valor original */
  cx = save;

  /* verificando se houve mudanca no posto de DFy[ord] */
  if (r != raux) {
   
    /* retornando o novo posto */
    return (raux);

  }

  /* selecionando a variavel y */
  for (k = 0; k <= o; k++) 
    for (l = 1; l <= n; l++) {

      /* salvando e perturbando a variavel y[k][l] */
      save      = cy[k][l];
      cy[k][l] -= 1.0e-6;

      /* avaliando a jacobiana */
      if (nDH == 1) {

        /* jacobiana exata */
        par->DF(o,n,cx,cy,DFx,DFy);

      } else {
    
        /* jacobiana aproximada */
        DFAPPROX(n,o,dir,cx,cy,delta,deltaaux,DFx,DFy,F); 

      }
      (*naDH)++;
 
      /* verificando o posto de DFy[ord] */
      /* B = DFy[ord]                    */
      for (i = 1; i <= n; i++)
        for (j = 1; j <= n; j++)
          B[i][j] = DFy[ord][j][i];

      /* decomposicao QR de DFy[ord] */
      raux = QR2(n,B,Q,p,q);
      (nQR)++;
    
      /* voltando a variavel y[k][l] ao valor original */
      cy[k][l] = save;

      /* verificando se houve mudanca no posto de DFy[ord] */
      if (r != raux) {
   
        /* retornando o novo posto */
        return (raux);

      }

    }


  /* retornando o posto inalterado */
  return (raux);
}



/***********************************************************/
/* rotina que calucula uma aproximacao para DF em um ponto */
/***********************************************************/

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
)
{
  int  i, j, k;   /* variaveis auxiliares        */
  real del;       /* incremento                  */
  real save;      /* salva o valor anterior      */
  real uround;    /* menor constante considerada */

  
  /* calculando a menor constante */
  uround = dir*sqrt(1.0e-15);

  /* calculo da derivada aproximada com relacao a x */

  /* calculando o incremento */
  del   = uround*fabs(x);
  del   = (x+del)-x;

  /* salvando o valor anterior */
  save  = x;

  /* incrementando a variavel */
  x    += del;
  del   = 1.0/del;

  /* avaliando a funcao */
  F(o,n,x,y,deltaaux);

  /* calculando uma aproximacao para a variavel */
  for (i = 1; i <= n; i++)
    DFx[i] = (deltaaux[i]-delta[i])*del;

  /* retornando ao valor anterior */
  x     = save;

  /* calculo da derivada aproximada com relacao a x */
  for (k = o; k >= 0; k--)
    for (i = 1; i <= n; i++) {

      /* calculando o incremento */
      del       = uround*fabs(y[k][i]);
      del       = (y[k][i]+del)-y[k][i];

      /* salvando o valor anterior */
      save      = y[k][i];

      /* incrementando a variavel */
      y[k][i]  += del;
      del       = 1.0/del;

      /* avaliando a funcao */
      F(o,n,x,y,deltaaux);

      /* calculando uma aproximacao para a variavel */
      for (j = 1; j <= n; j++)
        DFy[k][i][j] = (deltaaux[j]-delta[j])*del;

      /* retornando ao valor anterior */
      y[k][i]  = save;

    }

  return;
} 
/* fim DFAPPROX */



/* ************************************************************* */
/* Esta rotina inicializa todos os dados para a aplicacao do     */
/* metodo BDF.                                                   */
/* ************************************************************* */

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
)
{  
  int   i,j; /* variaveis auxiliares                */

  /* inicializando phi */
  phix[1] = cx;
  phix[2] = h*taux;
  for (i = 0; i <= o; i++) {
    for (j = 1; j <= n; j++) { 
      phiy[1][i][q[j]] = cy[i][q[j]] ;
      phiy[2][i][q[j]] = h*tauy[i][q[j]];
    }
  }  
        
  /* inicializando psi */
  psi[1] = h;

  /* inicializacao para realizar o primeiro passo */
  *cj     = 1.0/h ;
  *cjold  = *cj;
  *hold   = 0.0;
  *kold   = 0;
  *k      = 1;
  *ns     = 0; 
  *ifase  = 0;
  *pdcx   = taux;
  *factor = 100.0;

  return; 
} 
/* fim firststep */



int 
functionnorm (
int   n,
vreal x,
vreal tol
)
{
  int i;

  for (i = 1; i <= n; i++)
    if (fabs(x[i]) > tol[i]) return (0);

  return (1);
}



/****************************************************/
/* rotina que retorna x^y                           */
/****************************************************/

real 
POWER (
real x,
int  y
)
{
  real z;

  if (x == 0.0) {
    if (y == 0) 
      return (1.0); 
    else 
      return (0.0);
  } else {
    if ((abs(y) % 2) == 0) 
      z = exp(abs(y)*log(fabs(x)));
    else {
      if (x > 0.0) 
        z = exp(abs(y)*log(x)); 
      else 
        z = -exp(abs(y)*log(-x));
    }
    if (y > 0) return (z); else return (1.0/z);
  }
}



/****************************************************/
/* rotina que retorna x^(-y)                        */
/****************************************************/

real 
ROOT ( 
real x,
int  y
)
{
  real z;

  if ((abs(y) % 2) == 0) 
    z = exp(log(fabs(x))/abs(y));
  else {
    if (x > 0.0) 
      z = exp(log(x)/abs(y)); 
    else 
      z = -exp(log(-x)/abs(y));
  }
  if (y > 0) 
    return (z); 
  else 
    return (1.0/z);
}


void    
PIVOT (
int   n,
int   k,
mreal A,
mreal Q
) 
{
  int   i,imax;
  real  max;
  vreal aux;

  max  = fabs(A[k][k]);
  imax = k;

  for (i = k+1; i <= n; i++) 
    if (fabs(A[i][k]) > max) {
      imax = i;
      max  = fabs(A[i][k]);
    }

  if (imax != k) {
    aux     = A[k];
    A[k]    = A[imax];
    A[imax] = aux;
    aux     = Q[k];
    Q[k]    = Q[imax];
    Q[imax] = aux;
  }

  return;
}



/****************************************************************************/
/* Esta rotina aplica a matriz de rotacao de Givens nas matrizes A e Q para */
/* a decomposicao A = QR.                                                   */
/****************************************************************************/

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
)
{
  int   j;   /* variavel auxiliar para controle de lacos  */
  real  s,t; /* variaveis auxiliares para armazenar dados */

  if (fabs(s2)+fabs(s1) > 0.0) {
    if (fabs(s2) >= fabs(s1))
      s  = sqrt(1.0+(s1/s2)*(s1/s2))*fabs(s2);
    else
      s  = sqrt(1.0+(s2/s1)*(s2/s1))*fabs(s1);
    s1 = s1/s; s2 = s2/s; 
    /* aplicando a matriz de rotacao em A */
    for (j = 1; j <= n; j++) {
      s       =  s1*A[i][j]+s2*A[k][j];
      t       = -s2*A[i][j]+s1*A[k][j];
      A[i][j] =  s; 
      A[k][j] =  t;
    }
    /* aplicando a matriz de rotacao em Q */
    for (j = 1; j <= m; j++) {
      s       =  s1*Q[i][j]+s2*Q[k][j];
      t       = -s2*Q[i][j]+s1*Q[k][j];
      Q[i][j] =  s;
      Q[k][j] =  t;
    }
  }

  return;
}



/**************************************************************************/
/* Esta rotina retorna a decomposicao QR de uma matriz A (n+1)x(n) usando */
/* o metodo de Givens. A matriz R volta em A e em Q volta a matriz Q^t da */
/* decomposicao A = QR. Esta rotina tambem calcula uma estimativa para a  */
/* condicao da matriz A.                                                  */
/**************************************************************************/

void 
QR (
int    m,
int    n,
mreal  A,
mreal  Q,
int    pivot,
real  *cond
)
{
  int   i,j; /* variaveis auxiliares para controle de lacos       */

  /* inicializando Q com a identidade */
  for (i = 1; i <= m; i++) {
    for (j = 1; j <= m; j++) Q[i][j] = 0.0;
    Q[i][i] = 1.0;
  } 

  /* decomposicao QR por Givens */
  for (i = 1; i <= m-1; i++) {
    /* escolha do pivo */
    if (pivot)
      PIVOT(m,i,A,Q);
    for (j = i+1; j <= m; j++) {
      /* escolha da matriz de rotacao */
      GIVENS(m,n,i,j,A,Q,A[i][i],A[j][i]);
    }
  }

  /* calculo de uma estimativa para a condicao da matriz A */
  if (n == 1) {
   *cond = fabs(1.0/A[1][1]) ;
  } else {
    for (i = 2, (*cond) = 0.0; i <= n; i++)
      for (j = 1; j <= i-1; j++) *cond = MAX2(*cond,fabs(A[j][i]/A[i][i]));
  }

  return;
}



/**********************************************************************/
/* Esta rotina calcula um passo do metodo de Newton modificado para   */
/* sistemas indeterminados do tipo F(x) = 0, que e dado por           */
/*                  x_{k+1} = x_k - A^+ F(x_k)                        */
/* onde A^+ e a inversa de Moore-Penrouse de A = DF(x_0).             */
/*                                                                    */
/* Tendo a decomposicao A^t = QR, entao :                             */
/*                      A^+ = Q^t / R^(-t) \                          */
/*                                \  0^t   /                          */
/**********************************************************************/

void 
NEWTON (
int    n,
mreal  Q,
mreal  A,
vreal  u,
vreal  delta,
real   ac
)
{
  int   i,j; /* variaveis auxiliares para controle de lacos         */
  real  s;   /* variavel auxiliar para calculo da norma de A^+ F(u) */

  /* calculo de u = Q^t y */
  for (i = 1; i <= n; i++) {
    for (j = 1, s = 0.0; j <= n; j++) s += Q[i][j]*delta[j];
    u[i] = ac*s;
  }
  
  /* calculo de A b = y, onde y <- b */
  for (i = n; i >= 1; i--) {
    for (j = i+1; j <= n; j++) u[i] -= A[i][j]*u[j];
    u[i] /= A[i][i];
  } 

  return;
}


real    
PIVOT2 (
int   n,
int   k,
mreal A,
vint  p,
vint  q
) 
{
  int  i,j,imax,jmax;
  real max;

  max  = fabs(A[p[k]][q[k]]);
  imax = jmax = k;

  for (i = k; i <= n; i++) 
    for (j = k; j <= n; j++) 
      if (fabs(A[p[i]][q[j]]) > max) {
        imax = i;
        jmax = j;
        max  = fabs(A[p[i]][q[j]]);
      }

  if (imax != k) {
    i       = p[k];
    p[k]    = p[imax];
    p[imax] = i;
  }

  if (jmax != k) {
    j       = q[k];
    q[k]    = q[jmax];
    q[jmax] = j;
  }

  return (max);
}


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
)
{
  int   j;   /* variavel auxiliar para controle de lacos  */
  real  s,t; /* variaveis auxiliares para armazenar dados */

  if (fabs(s2)+fabs(s1) > 0.0) {
    if (fabs(s2) >= fabs(s1))
      s  = sqrt(1.0+(s1/s2)*(s1/s2))*fabs(s2);
    else
      s  = sqrt(1.0+(s2/s1)*(s2/s1))*fabs(s1);
    s1 = s1/s; s2 = s2/s; 
    /* aplicando a matriz de rotacao em A e Q */
    for (j = 1; j <= n; j++) {
      s             =  s1*A[p[i]][q[j]]+s2*A[p[k]][q[j]];
      t             = -s2*A[p[i]][q[j]]+s1*A[p[k]][q[j]];
      A[p[i]][q[j]] =  s; 
      A[p[k]][q[j]] =  t;
      s             =  s1*Q[p[i]][q[j]]+s2*Q[p[k]][q[j]];
      t             = -s2*Q[p[i]][q[j]]+s1*Q[p[k]][q[j]];
      Q[p[i]][q[j]] =  s;
      Q[p[k]][q[j]] =  t;
    }
  }

  return;
}


int 
QR2 (
int    n,
mreal  A,
mreal  Q,
vint   p,
vint   q
)
{
  int   i,k; 

  for (i = 1; i <= n; i++) {
    for (k = 1; k <= n; k++) Q[i][k] = 0.0;
    p[i]    = i;
    q[i]    = i;
    Q[i][i] = 1.0;
  } 

  for (i = 1; i < n; i++) {
    if (PIVOT2(n,i,A,p,q) < 1.0e-15) 
      return (i-1);
    for (k = i+1; k <= n; k++) 
      GIVENS2(n,i,k,A,Q,p,q,A[p[i]][q[i]],A[p[k]][q[i]]);
  }

  if (fabs(A[p[n]][q[n]]) < 1.0e-15) 
    return (n-1);

  return (n);
}


void 
SOLVESYSTEM (
int    n,
mreal  A,
mreal  Q,
vint   p,
vint   q,
vreal  x,
vreal  y
)
{
  int   i,j; 
  real  s;   

  for (j = 1; j <= n; j++) x[j] = y[j];

  for (i = 1; i <= n; i++) 
    for (j = 1, y[p[i]] = 0.0; j <= n; j++) 
      y[p[i]] += Q[p[i]][q[j]]*x[q[j]];
  
  for (i = n; i >= 1; i--) {
    for (j = i+1, s = 0.0; j <= n; j++) s += A[p[i]][q[j]]*x[q[j]];
    x[q[i]] = (y[p[i]]-s)/A[p[i]][q[i]];
  } 
  
  return;
}


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
)
{
  *y     = (mreal) ALLOCMREAL(o,n);
  if (*y == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  *atoly = (mreal) ALLOCMREAL(o,n);
  if (*atoly == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  *rtoly = (mreal) ALLOCMREAL(o,n);
  if (*rtoly == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  *ftol = (vreal) ALLOCVREAL(n);
  if (*ftol == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  *infoinput  = (vint)  ALLOCVINT(2*n+10);
  if (*infoinput == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  *infooutput  = (vint)  ALLOCVINT(2*n+10);
  if (*infooutput == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  par = (parameter *) malloc(sizeof(parameter));
  if (par == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  /* inicializando a ordem e a dimensao com 0 */
  par->n = 0;
  par->o = 0;

  /* definindo as funcoes */
  par->F  = F;
  par->DF = DF;

  /* aloca vetores inteiros */
  par->p = (vint) ALLOCVINT(n);
  if ( par->p == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  par->q = (vint) ALLOCVINT(n);
  if ( par->q == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  par->paux = (vint) ALLOCVINT(n);
  if ( par->paux == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  par->qaux = (vint) ALLOCVINT(n);
  if ( par->qaux == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }


  /* aloca vetores reais */
  par->DFx = (vreal) ALLOCVREAL(n);
  if ( par->DFx == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  par->u  = (vreal) ALLOCVREAL((o+1)*n+1);
  if ( par->u == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  par->delta = (vreal) ALLOCVREAL(n+1);
  if ( par->delta == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  par->deltax = (vreal) ALLOCVREAL(n+1);
  if ( par->deltax == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }
  
  par->deltah = (vreal) ALLOCVREAL((o+1)*n+1);
  if ( par->delta == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }
  
  par->deltahx = (vreal) ALLOCVREAL((o+1)*n+1);
  if ( par->delta == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  par->psi = (vreal) ALLOCVREAL(8);
  if ( par->psi == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  par->alfa = (vreal) ALLOCVREAL(8);
  if ( par->alfa == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  par->beta = (vreal) ALLOCVREAL(8);
  if ( par->beta == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  par->gama = (vreal) ALLOCVREAL(8);
  if ( par->gama == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  par->sigma = (vreal) ALLOCVREAL(8);
  if ( par->sigma == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  par->phix = (vreal) ALLOCVREAL(8);
  if ( par->phix == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }
 
  par->ftol = (vreal) ALLOCVREAL(n);
  if ( par->ftol == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }
 
  /* aloca matrizes bidimensionais */
 
  par->y = (mreal) ALLOCMREAL(o,n);
  if ( par->y == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  par->tauy = (mreal) ALLOCMREAL(o,n);
  if ( par->tauy == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  par->DH = (mreal) ALLOCMREAL((o+1)*n+1,(o+1)*n+1);
  if ( par->DH == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  par->cy = (mreal) ALLOCMREAL(o,n);
  if ( par->cy == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  par->cyx = (mreal) ALLOCMREAL(o,n);
  if ( par->cyx == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  par->pcy = (mreal) ALLOCMREAL(o,n);
  if ( par->pcy == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  par->pdcy = (mreal) ALLOCMREAL(o,n);
  if ( par->pdcy == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  par->dy = (mreal) ALLOCMREAL(o,n);
  if (par->dy == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  par->ccy = (mreal) ALLOCMREAL(o,n);
  if ( par->ccy == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  par->yx = (mreal) ALLOCMREAL(o,n);
  if ( par->yx == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  par->Q = (mreal) ALLOCMREAL((o+1)*n+1,(o+1)*n+1);
  if ( par->Q == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  par->atoly = (mreal) ALLOCMREAL(o,n);
  if ( par->atoly == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  par->rtoly = (mreal) ALLOCMREAL(o,n);
  if ( par->rtoly == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  par->wty = (mreal) ALLOCMREAL(o,n);
  if ( par->wty == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  par->Ey = (mreal) ALLOCMREAL(o,n);
  if ( par->Ey == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }
   
  /* aloca matrizes tridimensionais */
   
  par->DFy  = (mmreal) ALLOCMMREAL(o,n,n);
  if ( par->DFy == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  par->phiy = (mmreal) ALLOCMMREAL(8,o,n);
  if ( par->phiy == NULL) {
    printf("ALLOCPAR : nao alocado\n");
    exit(1);
    return;
  }

  return;
}




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
)
{
  extern parameter *par;

  *y           = (mreal) FREEMREAL(o,n,*y);
  *atoly       = (mreal) FREEMREAL(o,n,*atoly);
  *rtoly       = (mreal) FREEMREAL(o,n,*rtoly);
  *ftol        = (vreal) FREEVREAL(n,*ftol);
  *infoinput   = (vint)  FREEVINT(2*n+10,*infoinput);
  *infooutput  = (vint)  FREEVINT(2*n+10,*infooutput);

  /* desaloca vetores inteiros */
  par->p       = (vint)  FREEVINT(n,par->p);
  par->q       = (vint)  FREEVINT(n,par->q);
  par->paux    = (vint)  FREEVINT(n,par->paux);
  par->qaux    = (vint)  FREEVINT(n,par->qaux);

  /* desaloca vetores reais */
  par->DFx     = (vreal) FREEVREAL(n,par->DFx);
  par->u       = (vreal) FREEVREAL((o+1)*n+1,par->u);
  par->delta   = (vreal) FREEVREAL(n+1,par->delta);
  par->deltax  = (vreal) FREEVREAL(n+1,par->deltax);
  par->deltah  = (vreal) FREEVREAL((o+1)*n+1,par->delta);
  par->deltahx = (vreal) FREEVREAL((o+1)*n+1,par->deltax);
  par->psi     = (vreal) FREEVREAL(8,par->psi);
  par->alfa    = (vreal) FREEVREAL(8,par->alfa);
  par->beta    = (vreal) FREEVREAL(8,par->beta);
  par->gama    = (vreal) FREEVREAL(8,par->gama);
  par->sigma   = (vreal) FREEVREAL(8,par->sigma);
  par->phix    = (vreal) FREEVREAL(8,par->phix);
  par->ftol    = (vreal) FREEVREAL(n,par->ftol);
   
  /* desaloca matrizes bidimensionais */
  par->y       = (mreal) FREEMREAL(o,n,par->y);
  par->tauy    = (mreal) FREEMREAL(o,n,par->tauy); 
  par->DH      = (mreal) FREEMREAL((o+1)*n+1,(o+1)*n+1,par->DH);
  par->cy      = (mreal) FREEMREAL(o,n,par->cy);
  par->cyx     = (mreal) FREEMREAL(o,n,par->cyx);
  par->pcy     = (mreal) FREEMREAL(o,n,par->pcy);
  par->pdcy    = (mreal) FREEMREAL(o,n,par->pdcy);
  par->dy      = (mreal) FREEMREAL(o,n,par->dy);
  par->ccy     = (mreal) FREEMREAL(o-1,n,par->ccy);
  par->yx      = (mreal) FREEMREAL(o,n,par->yx);
  par->Q       = (mreal) FREEMREAL((o+1)*n+1,(o+1)*n+1,par->Q);
  par->atoly   = (mreal) FREEMREAL(o,n,par->atoly);
  par->rtoly   = (mreal) FREEMREAL(o,n,par->rtoly);
  par->wty     = (mreal) FREEMREAL(o,n,par->wty);
  par->Ey      = (mreal) FREEMREAL(o,n,par->Ey);
  
  /* desaloca matrizes tridimensionais */
  par->DFy     = (mmreal) FREEMMREAL(o,n,n,par->DFy);
  par->phiy    = (mmreal) FREEMMREAL(8,o,n,par->phiy);
  
  free(par);

  return;
}


/************************************************************/
/* Funcao para alocacao de vetor de inteiros                */
/************************************************************/
/* Aloca um vetor de inteiros de dimensao n, inserindo      */
/* na primeira posicao (v[0]) a dimensao n.                 */
/************************************************************/
/* PARAMETROS DE ENTRADA: a dimensao do vetor a ser alocado */
/* VALOR RETORNADO: ponteiro para o vetor                   */
/* FUNCOES ATIVADAS: malloc()                               */
/************************************************************/

vint   
ALLOCVINT (
int n
)
{
  vint v;  /* ponteiro para o vetor */
  int  i;  /* variavel auxiliar     */

  /* alocando o vetor */
  v = (vint) malloc((n+1)*sizeof(int));
  if (v == NULL) {
    printf("VECTORS : nao alocado\n");
    return (NULL);
  }
  for (i = 0; i <= n; i++) v[i] = 0;

  /* retornando o ponteiro */
  return (v);
}

/**********************************************************/
/* Funcao para redimensionar um conjunto de inteiros      */
/**********************************************************/
/* Redimensiona um conjunto de inteiros com m elementos   */
/* para conter n elementos, ou seja, realoca um vetor de  */
/* inteiros de dimensao m para um vetor de dimensao n     */
/* para dimensao n.                                       */
/**********************************************************/
/* PARAMETROS DE ENTRADA: ponteiro para o conjunto a ser  */
/*                        redimensionado e nova dimensao  */
/* VALOR RETORNADO: ponteiro para o novo conjunto         */ 
/* FUNCOES ATIVADAS: malloc()                             */
/**********************************************************/

vint   
REALLOCVINT ( 
int  m,
int  n,
vint v
)
{
  vint w;  /* ponteiro para o novo vetor                 */
  int  i;  /* variaveis auxiliares                       */

  /* alocando o novo vetor                                */
  w = (vint) malloc((n+1)*sizeof(int));
  if (w == NULL) {
    printf("VECTORS : nao alocado\n");
    return (NULL);
  }

  /* copiando o vetor no novo vetor */
  if (m < n) {
    for (i = 1; i <= m; i++)   w[i] = v[i];
    for (i = m+1; i <= n; i++) w[i] = 0;
  } else {
    for (i = 1; i <= n; i++) w[i] = v[i];
  }

  /* liberando o vetor antigo*/
  free(v);

  /* retornando o ponteiro para o novo vetor */
  return (w);
}


/*********************************************************/
/* Funcao para liberacao de um vetor de inteiros         */
/*********************************************************/
/* PARAMETROS DE ENTRADA: ponteiro para vetor e dimensao */
/* VALOR RETORNADO: ponteiro para o vetor                */
/* FUNCOES ATIVADAS: free()                             */
/*********************************************************/

vint  
FREEVINT ( 
int  n,
vint v
)
{
  /* liberando o vetor */
  free(v);

  /* retornando o ponteiro */
  return (NULL);
}


/**************************************************/
/* Funcao para alocacao de uma matriz de inteiros */
/**************************************************/
/* PARAMETROS DE ENTRADA: dimensoes da matriz     */
/* VALOR RETORNADO: ponteiro para matriz          */
/* FUNCOES ATIVADAS: malloc()                     */
/**************************************************/

mint   
ALLOCMINT ( 
int m,
int n
)
{
  mint v;    /* ponteiro para o vetor */
  int  i,j;  /* variaveis auxiliares  */

  /* alocando as linhas da matriz */
  v = (mint) malloc((m+1)*sizeof(vint));
  if (v == NULL) {
    printf("VECTORS : nao alocado\n");
    return (NULL);
  }

  /* alocando as colunas da matriz */
  for (i = 0; i <= m; i++) {
    v[i] = (vint) malloc((n+1)*sizeof(int));
    if (v[i] == NULL) {
      printf("VECTORS : nao alocado\n");
      return (NULL);
    }
    for (j = 0; j <= n; j++) v[i][j] = 0;
  }

  /* retornando o ponteiro */
  return (v);
}


/***********************************************************/
/* Funcao para liberar uma matriz de inteiros              */
/***********************************************************/
/* PARAMETROS DE ENTRADA: ponteiro para matriz e dimensoes */
/* VALOR RETORNADO: ponteiro para matriz                   */
/* FUNCOES ATIVADAS: free()                               */
/***********************************************************/

mint  
FREEMINT ( 
int  m,
int  n,
mint v
)
{
  int      i;  /* variavel auxiliar */

  /* liberando as linhas da matriz */
  for (i = 0; i <= m ; i++) free(v[i]);

  /* liberando a matriz */
  free(v);

  /* retornando o ponteiro */
  return (NULL);
}


/***********************************************************/
/* Esta rotina aloca uma matriz tridimensional de inteiros */
/***********************************************************/

mmint  
ALLOCMMINT ( 
int m,
int n,
int k
)
{
  mmint v;      /* ponteiro para a matriz */
  int   i,j,l;  /* variaveis auxiliares   */

  /* alocando a primeira dimensao da matriz */
  v = (mmint) malloc((m+1)*sizeof(mint));
  if (v == NULL) {
    printf("VECTORS : nao alocado\n");
    return (NULL);
  }

  /* alocando a segunda dimensao da matriz */
  for (i = 0; i <= m; i++) {
    v[i] = (mint) malloc((n+1)*sizeof(vint));
    if (v[i] == NULL) {
      printf("VECTORS : nao alocado\n");
      return (NULL);
    }
  }

  /* alocando a terceira dimensao da matriz */
  for (i = 0; i <= m; i++) {
    for (j = 0; j <= n; j++) {
      v[i][j] = (vint) malloc((k+1)*sizeof(int));
      if (v[i][j] == NULL) {
        printf("VECTORS : nao alocado\n");
        return (NULL);
      }
      for (l = 0; l <= k; l++) v[i][j][l] = 0;
    }
  }

  /* retornando o ponteiro */
  return (v);
}


/************************************************************/
/* Esta rotina libera uma matriz tridimensional de inteiros */
/************************************************************/

mmint 
FREEMMINT (
int   m,
int   n,
int   k,
mmint v
)
{
  int      i,j;  /* variaveis auxiliares */

  /* liberando a terceira dimensao da matriz */
  for (i = 0; i <= m; i++) 
    for (j = 0; j <= n; j++) free(v[i][j]);

  /* liberando a segunda dimensao da matriz */
  for (i = 0; i <= m; i++) free(v[i]);

  /* liberando a primeira dimensao da matriz */
  free(v);

  /* retornando o ponteiro */
  return (NULL);
}



/********************************************/
/* Funcao para alocar um vetor de reais     */
/********************************************/
/* PARAMETROS DE ENTRADA: dimensao do vetor */
/* VALOR RETORNADO: ponteiro para o vetor   */
/* FUNCOES ATIVADAS: malloc()               */
/********************************************/

vreal  
ALLOCVREAL (
int n
)
{
  vreal v;  /* ponteiro para o vetor */
  int  i;   /* variavel auxiliar     */

  /* alocando o vetor */
  v = (vreal) malloc((n+1)*sizeof(real));
  if (v == NULL) {
    printf("VECTORS : nao alocado\n");
    return (NULL);
  }
  for (i = 0; i <= n; i++) v[i] = 0.0;

  /* retornando o ponteiro */
  return (v);
}


/**********************************************************/
/* Funcao para redimensionar um conjunto de reais         */
/**********************************************************/
/* Redimensiona um conjunto de reais com m elementos      */
/* para conter n elementos, ou seja, realoca um vetor de  */
/* inteiros de dimensao m para um vetor de dimensao n     */
/* para dimensao n.                                       */
/**********************************************************/
/* PARAMETROS DE ENTRADA: ponteiro para o conjunto a ser  */
/*                         redimensionado e nova dimensao */
/* VALOR RETORNADO: ponteiro para o novo conjunto         */ 
/* FUNCOES ATIVADAS: malloc()                             */
/**********************************************************/

vreal   
REALLOCVREAL (
int   m,
int   n,
vreal v
)
{
  vreal w;  /* ponteiro para o novo vetor                 */
  int   i;  /* variaveis auxiliares                       */

  /* alocando o novo vetor                                */
  w = (vreal) malloc((n+1)*sizeof(real));
  if (w == NULL) {
    printf("VECTORS : nao alocado\n");
    return (NULL);
  }

  /* copiando o vetor no novo vetor */
  if (m < n) {
    for (i = 1; i <= m; i++)   w[i] = v[i];
    for (i = m+1; i <= n; i++) w[i] = 0.0;
  } else {
    for (i = 1; i <= n; i++) w[i] = v[i];
  }

  /* liberando o vetor antigo*/
  free(v);

  /* retornando o ponteiro para o novo vetor */
  return (w);
}


/*********************************************************/
/* Funcao para liberar um vetor de reais                 */
/*********************************************************/
/* PARAMETROS DE ENTRADA: ponteiro para vetor e dimensao */
/* VALOR RETORNADO: ponteiro para o vetor                */
/* FUNCOES ATIVADAS: free()                             */
/*********************************************************/

vreal 
FREEVREAL (
int   n,
vreal v
)
{
  /* liberando o vetor */
  free(v);

  /* retornando o ponteiro */
  return (NULL);
}


/**********************************************/
/* Funcao para alocar uma matriz de reais     */
/**********************************************/
/* PARAMETROS DE ENTRADA: dimensoes da matriz */
/* VALOR RETORNADO: ponteiro para a matriz    */
/* FUNCOES ATIVADAS: free()                  */
/**********************************************/

mreal  
ALLOCMREAL (
int m,
int n
)
{
  mreal v;   /* ponteiro para a matriz */
  int   i,j; /* variaveis auxiliares   */

  /* alocando as linhas da matriz */
  v = (mreal) malloc((m+1)*sizeof(vreal));
  if (v == NULL) {
    printf("VECTORS : nao alocado\n");
    return (NULL);
  }

  /* alocando as colunas da matriz */
  for (i = 0; i <= m; i++) {
    v[i] = (vreal) malloc((n+1)*sizeof(real));
    if (v[i] == NULL) {
      printf("VECTORS : nao alocado\n");
      return (NULL);
    }
    for (j = 0; j <= n; j++) v[i][j] = 0.0;
  }

  /* retornando o ponteiro */
  return (v);
}


/***********************************************************/
/* Funcao para liberar uma matriz de reais                 */
/***********************************************************/
/* PARAMETROS DE ENTRADA: ponteiro para matriz e dimensoes */
/* VALOR RETORNADO: ponteiro para matriz                   */
/* FUNCOES ATIVADAS: free()                               */
/***********************************************************/

mreal 
FREEMREAL ( 
int   m,
int   n,
mreal v
)
{
  int      i;  /* variavel auxiliar */

  /* liberando as linhas da matriz */
  for (i = 0; i <= m; i++) free(v[i]);

  /* libreando a matriz */
  free(v);

  /* retornando o ponteiro */
  return (NULL);
}


/********************************************************/
/* Esta rotina aloca uma matriz tridimensional de reais */
/********************************************************/

mmreal  
ALLOCMMREAL ( 
int m,
int n,
int k
)
{
  mmreal v;      /* ponteiro para a matriz */
  int    i,j,l;  /* variaveis auxiliares   */

  /* alocando a primeira dimensao da matriz */
  v = (mmreal) malloc((m+1)*sizeof(mreal));
  if (v == NULL) {
    printf("VECTORS : nao alocado\n");
    return (NULL);
  }

  /* alocando a segunda dimensao da matriz */
  for (i = 0; i <= m; i++) {
    v[i] = (mreal) malloc((n+1)*sizeof(vreal));
    if (v[i] == NULL) {
      printf("VECTORS : nao alocado\n");
      return (NULL);
    }
  }

  /* alocando a terceira dimensao da matriz */
  for (i = 0; i <= m; i++) {
    for (j = 0; j <= n; j++) {
      v[i][j] = (vreal) malloc((k+1)*sizeof(real));
      if (v[i][j] == NULL) {
        printf("VECTORS : nao alocado\n");
        return (NULL);
      }
      for (l = 0; l <= k; l++) v[i][j][l] = 0.0;
    }
  }

  /* retornando o ponteiro */
  return (v);
}


/*********************************************************/
/* Esta rotina libera uma matriz tridimensional de reais */
/*********************************************************/

mmreal 
FREEMMREAL (
int    m,
int    n,
int    k,
mmreal v
)
{
  int      i,j;  /* variaveis auxiliares */

  /* liberando a terceira dimensao da matriz */
  for (i = 0; i <= m; i++) 
    for (j = 0; j <= n; j++) free(v[i][j]);

  /* liberando a segunda dimensao da matriz */
  for (i = 0; i <= m; i++) free(v[i]);

  /* liberando a primeira dimensao da matriz */
  free(v);

  /* retornando o ponteiro */
  return (NULL);
}
