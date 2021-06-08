#include <stdio.h>
#include <math.h>
#include <string.h> 
#include <stdarg.h>
#include <stdlib.h>

#define LX 512 
#define LY 128  
#define F 9
#define A 32

#define M_PI 3.14159265358979323846
#define cs_sq (1. / 3.)
#define wC (4. / 9.)
#define wS (1. / 9.)
#define wL (1. / 36.)

enum LetVel {_C, _N, _S, _W, _E, _NE, _NW, _SW, _SE};
double cx[F] = {0, 0, 0, -1 ,1, 1, -1, -1, 1};
double cy[F] = {0, 1, -1, 0, 0, 1,  1, -1, -1};
double w[F] = {wC, wS, wS, wS, wS, wL, wL, wL, wL};

double cells[LX * LY * F];
double tmp_cells[LX * LY * F];

#define AR(x, y, f) ((f) + F * (x) + F * LX * (y)) 
#define ARM(x, y) ((x) + LX * (y))

double Phi[LX * LY];
double Rho[LX * LY];
double Ux[LX * LY];
double Uy[LX * LY];
int map[LX*LY];

double Re, Ma;
double Pe, phiH, phiL;

double tau, alphaLB;
double uLB, vLB;
double dx, dt;
double rhoLB_IN, uxLB_IN, uyLB_IN;
double phi0, u0, v0;

double fEq(int i, double rho, double ux, double uy) {
  double cu = (cx[i]*ux + cy[i]*uy) / cs_sq; 
    return rho*w[i]*(1. + cu + 0.5*cu*cu - 0.5*(ux*ux+uy*uy)/cs_sq);;
}

void collideAndStream(double *c, double *tc)
{
  int x, y, i, xm, ym, xp, yp;

  double rho, ux, uy, fstar[F];

  for (y = 0; y < LY; y++)
    for (x = 0; x < LX; x++)
    {
      // 0. BC

      /*if ((y == LY - 1) && (x != 0) && (x != LX - 1))
        for (i = 0; i < F; i++)
          c[AR(x, y, i)] = fEq(i, rhoLB_IN, uxLB_IN, uyLB_IN);*/
      if ((x == 0) && !map[ARM(x, y)]) {
        for(i = 0; i < F; i++) {
          c[AR(x, y, i)] = fEq(i, rhoLB_IN, uxLB_IN, uyLB_IN);
        }
      }
      if((x == LX-1) && !map[ARM(x, y)]) {
        for(i = 0; i < F; i++) {
          c[AR(x, y, i)] = c[AR(x - 1, y, i)];
        }
      }
      // 1. makro
      rho = 0;
      ux = 0;
      uy = 0;
      for (i = 0; i < F; i++)
      {
        rho += c[AR(x, y, i)];
        ux += c[AR(x, y, i)] * cx[i];
        uy += c[AR(x, y, i)] * cy[i];
      }
      ux /= rho;
      uy /= rho;
      // save macro to array
      Rho[ARM(x, y)] = rho;
      Ux[ARM(x, y)] = ux;
      Uy[ARM(x, y)] = uy;

      // 2. kolizja
      if (!map[ARM(x, y)])
      {
        double fi_E = c[AR(x, y, _E)];
        double fi_W = c[AR(x, y, _W)];
        double fi_S = c[AR(x, y, _S)];
        double fi_N = c[AR(x, y, _N)];
        double fi_SE = c[AR(x, y, _SE)];
        double fi_SW = c[AR(x, y, _SW)];
        double fi_NE = c[AR(x, y, _NE)];
        double fi_NW = c[AR(x, y, _NW)];

        double feq_E = fEq(_E, rho, ux, uy);
        double feq_W = fEq(_W, rho, ux, uy);
        double feq_S = fEq(_S, rho, ux, uy);
        double feq_N = fEq(_N, rho, ux, uy);
        double feq_SE = fEq(_SE, rho, ux, uy);
        double feq_SW = fEq(_SW, rho, ux, uy);
        double feq_NE = fEq(_NE, rho, ux, uy);
        double feq_NW = fEq(_NW, rho, ux, uy);

        double C = .18; // Smagorinsky constant
        double Q = (fi_E - feq_E + fi_W - feq_W + fi_NE - feq_NE + fi_SE - feq_SE + fi_NW - feq_NW + fi_SW - feq_SW) *
                       (fi_E - feq_E + fi_W - feq_W + fi_NE - feq_NE + fi_SE - feq_SE + fi_NW - feq_NW + fi_SW - feq_SW) +
                   (fi_S - feq_S + fi_N - feq_N + fi_NE - feq_NE + fi_SE - feq_SE + fi_NW - feq_NW + fi_SW - feq_SW) *
                       (fi_S - feq_S + fi_N - feq_N + fi_NE - feq_NE + fi_SE - feq_SE + fi_NW - feq_NW + fi_SW - feq_SW) +
                   2.f * (fi_NE - feq_NE - (fi_SE - feq_SE + fi_NW - feq_NW) + fi_SW - feq_SW) *
                       (fi_NE - feq_NE - (fi_SE - feq_SE + fi_NW - feq_NW) + fi_SW - feq_SW);

        double tau_t = .5 * (sqrt(tau * tau + C * C * sqrt(648. * Q) / rho) + tau);

        for (i = 0; i < F; i++)
          fstar[i] = (1. - 1. / tau_t) * c[AR(x, y, i)] + 1. / tau_t * fEq(i, rho, ux, uy);
      }
      else
      {
        // Bounce-back BC ~ wall
        fstar[_E] = c[AR(x, y, _W)];
        fstar[_W] = c[AR(x, y, _E)];
        fstar[_S] = c[AR(x, y, _N)];
        fstar[_N] = c[AR(x, y, _S)];
        fstar[_NE] = c[AR(x, y, _SW)];
        fstar[_NW] = c[AR(x, y, _SE)];
        fstar[_SE] = c[AR(x, y, _NW)];
        fstar[_SW] = c[AR(x, y, _NE)];
      }
      // 3. & 4. periodic BC & Streaming
      // BC
      xp = (x == LX - 1) ? 0 : x + 1;
      yp = (y == LY - 1) ? 0 : y + 1;
      xm = (x == 0) ? LX - 1 : x - 1;
      ym = (y == 0) ? LY - 1 : y - 1;

      // Streaming
      tc[AR(x, y, _C)] = fstar[_C];
      tc[AR(x, yp, _N)] = fstar[_N];
      tc[AR(x, ym, _S)] = fstar[_S];
      tc[AR(xm, y, _W)] = fstar[_W];
      tc[AR(xp, y, _E)] = fstar[_E];
      tc[AR(xp, yp, _NE)] = fstar[_NE];
      tc[AR(xp, ym, _SE)] = fstar[_SE];
      tc[AR(xm, yp, _NW)] = fstar[_NW];
      tc[AR(xm, ym, _SW)] = fstar[_SW];
    }
}
double PHI(double x){
   return 0.5*(1. + sin(2*M_PI*x));
}

void setInitialConditions(double *c, double rhoi, double uxi, double uyi){
  int x,y,i;
  for (y = 0; y < LY; y++){
    for (x = 0; x < LX ; x++){ 

      for (i = 0; i < F; i++){
        c[AR(x,y,i)] = fEq(i, rhoi, uxi, uyi);
        Rho[ARM(x,y)] = rhoi;
        Ux[ARM(x,y)] = uxi;
        Uy[ARM(x,y)] = uyi;
      }

      map[ARM(x,y)] = 0;
      /*if ( ( x==0)  || (x == LX-1) || (y == 0)){
        map[ARM(x,y)] = 1;
        Rho[ARM(x,y)] = 0.;
        Ux[ARM(x,y)] = 0.;
        Uy[ARM(x,y)] = 0.;
      }*/
      if ((y == 0) || (y == LY-1))
      {
        map[ARM(x, y)] = 1;
        Rho[ARM(x, y)] = 0.;
        Ux[ARM(x, y)] = 0.;
        Uy[ARM(x, y)] = 0.;
      }

      if((x >= .5*(LX*.5-A)) && (x < .5*(LX*.5+A)) && (y > .5*(LY-A)) && (y < .5*(LY+A))) {
        map[ARM(x, y)] = 1;
        Rho[ARM(x, y)] = 0.;
        Ux[ARM(x, y)] = 0.;
        Uy[ARM(x, y)] = 0.;
      }
    }
  }
}

void dumpStateVTK(char *fn, int iter){
 FILE *fp;
char fname[50];
 int x,y,dx,i;
 
    sprintf(fname,"%s%05d.vtk",fn,iter);
    fp = fopen(fname, "w");
  
 fprintf(fp,"# vtk DataFile Version 2.0\n");
 fprintf(fp,"2D-ADE data file \n");
 fprintf(fp,"ASCII\n");
 fprintf(fp,"DATASET RECTILINEAR_GRID\n");
 fprintf(fp,"DIMENSIONS %d %d %d\n", LX, LY, 1);
 fprintf(fp,"X_COORDINATES %d int\n", LX);
 for (x = 0; x < LX; x++) 
   fprintf(fp, "%d ", x);
 fprintf(fp,"\n");
 fprintf(fp,"Y_COORDINATES %d int\n", LY);
 for (x = 0; x < LY; x++)
   fprintf(fp, "%d ", x);
 fprintf(fp,"\n");
 fprintf(fp,"Z_COORDINATES 1 int\n");
 fprintf(fp, "0\n");
 fprintf(fp,"POINT_DATA %d \n", LX*LY);
 fprintf(fp,"SCALARS density double 1\n");
 fprintf(fp,"LOOKUP_TABLE default\n"); 
 for (y = 0; y < LY; y++)
   for (x = 0; x < LX; x++) {
   fprintf(fp, "%e \n", 1. -Rho[ARM(x,y)]); 
  }
  fprintf(fp,"VECTORS velocity double \n");
  for (y = 0; y < LY; y++) 
   for (x = 0; x < LX; x++) {

   fprintf(fp, "%e %e %e\n", Ux[ARM(x,y)], Uy[ARM(x,y)], 0.);
   }


    fclose(fp);
}

int main(int argc, char **argv){
  int iter = 0;
  int ITERMAX = 1e6;
  double rhoiLB, uxiLB, uyiLB, nuLB; 
  rhoiLB = 1.0;
  uxiLB = 0.0;
  uyiLB = 0.0;

  Re = 1e5;
  Ma = 0.1;

  rhoLB_IN = 1.;
  uxLB_IN = Ma*sqrt(cs_sq);
  uyLB_IN = 0.0; 
  
  nuLB = 0*uxLB_IN*LY/Re;
 
  tau = nuLB/cs_sq + .5;
  dx = 1. / (double)LX;

  setInitialConditions(cells, rhoiLB, uxiLB, uyiLB);
  dumpStateVTK("state",iter);

  int i, j;
  double f1 = 0.0, f2 = 0.0;

  do{
    if (iter % 2)
      collideAndStream(tmp_cells, cells);
    else
      collideAndStream(cells, tmp_cells);
    iter++;
    if (iter % 1000 == 0) dumpStateVTK("state",iter);
/*
    for (i = 0; i < LY; i++){
      for(j = 0; j < LY; j++) {
          f1 += Ux[ARM(i,j)]*Ux[ARM(i,j)] + Uy[ARM(i,j)]*Uy[ARM(i,j)];
      }
    }
    if (abs(f1-f2) < 1e-19 && iter >= 5000) break;
    f2 = f1;
    f1 = 0.0;*/
  }while (iter < ITERMAX);

  return 0;
}