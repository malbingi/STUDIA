#include <stdio.h>
#include <math.h>
#include <string.h> 
#include <stdarg.h>
#include <stdlib.h>

#define LX 128 
#define LY 128  
#define F 9

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

void collideAndStream(double *c, double *tc) {
  int i,x,y,xp,xm,yp,ym;
  double rho,ux,uy, fstar[F];
  //#pragma omp parallel for
  for (y = 0; y < LY; y++)
    for (x = 0; x < LX; x++){
      rho = 0; ux = 0; uy = 0;
      if (map[ARM(x,y)] != 1) {
        if (y==LY-1 ){
          for (i = 0; i < F; i++)
            c[AR(x,y,i)] = fEq(i, rhoLB_IN, uxLB_IN, uyLB_IN);
        }
        for (i = 0; i < F; i++){
            rho += c[AR(x,y,i)];
            ux += c[AR(x,y,i)]*cx[i];
            uy += c[AR(x,y,i)]*cy[i];
        }
        ux /= rho;
        uy /= rho;
        Rho[ARM(x,y)] = rho;
        Ux[ARM(x,y)] = ux;
        Uy[ARM(x,y)] = uy;
        for (i = 0; i < F ; i++){
          fstar[i]= (1. - 1. / tau) * c[AR(x, y, i)] + 1. / tau * fEq(i, rho, ux, uy);  
        }
      }
      else {
        fstar[_E] = c[AR(x,y, _W)];
        fstar[_W] = c[AR(x,y, _E)];
        fstar[_N] = c[AR(x,y, _S)];
        fstar[_S] = c[AR(x,y, _N)];
        fstar[_NE] = c[AR(x,y, _SW)];
        fstar[_NW] = c[AR(x,y, _SE)];
        fstar[_SW] = c[AR(x,y, _NE)];
        fstar[_SE] = c[AR(x,y, _NW)];
      }
      // dlaczego przeciwny?????
      xp = (x == LX - 1) ? 0 : x + 1; 
      yp = (y == LY - 1) ? 0 : y + 1;
      xm = (x == 0) ? LX - 1 : x - 1; 
      ym = (y == 0) ? LY - 1 : y - 1; 
      tc[AR(x, y, _C)] = fstar[_C];
      tc[AR(x,yp, _N)] = fstar[_N];
      tc[AR(x,ym, _S)] = fstar[_S];
      tc[AR(xm,y, _W)] = fstar[_W]; 
      tc[AR(xp,y, _E)] = fstar[_E];
      tc[AR(xp,yp, _NE)] = fstar[_NE];
      tc[AR(xp,ym, _SE)] = fstar[_SE];
      tc[AR(xm,yp, _NW)] = fstar[_NW];  
      tc[AR(xm,ym, _SW)] = fstar[_SW];
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
      if ( ( x==0)  || (x == LX-1) || (y == 0) ){
        map[ARM(x,y)] = 1;
        //Rho[ARM(x,y)] = 0;
        //Ux[ARM(x,y)] = 0;
        //Uy[ARM(x,y)] = 0;
      }
      map[ARM(x,y)] = 0;
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

  Re = 1000.0;
  Ma = 0.1;

  rhoLB_IN = 1.;
  uxLB_IN = Ma*sqrt(cs_sq);
  uyLB_IN = 0.0; 
  
  nuLB = uxLB_IN*LY/Re;
 
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
    if (iter % 10000 == 0) dumpStateVTK("state",iter);
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