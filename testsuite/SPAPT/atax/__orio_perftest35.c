/*
ACOPY_x:False
ACOPY_y:True
RT_K:32
T1_I:32
T1_J:512
T1_K:64
U_K:3
U_J:1
U_I:22
U1_I:14
T2_K:2048
T2_J:2048
T2_I:2048
VEC2:True
VEC1:False
RT_I:1
SCR:False
RT_J:1
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <time.h>

#define N 10000
#include "decl.h"
double orio_t_start, orio_t_end, orio_t = (double)LONG_MAX;

#include "init.c"

#include "validation.c"



extern double getClock(); 

int main(int argc, char *argv[]) {
  
#ifdef MAIN_DECLARATIONS
  MAIN_DECLARATIONS()
#endif  
  init_input_vars();


  int orio_i;

  /*
   Coordinate: [2, 6, 3, 6, 6, 6, 0, 1, 13, 21, 0, 2, 0, 0, 2, 0, 0, 1] 
  */
  
  
  for (orio_i=0; orio_i<ORIO_REPS; orio_i++) {
    orio_t_start = getClock();
    
    

#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

int i,j, k;
int it, jt, kt;
int ii, jj, kk;
int iii, jjj, kkk;

double* tmp=(double*) malloc(nx*sizeof(double));

/*@ begin Loop(

  transform Composite(
    unrolljam = (['i'],[U1_I]),
    vector = (VEC1, ['ivdep','vector always'])
  )
  for (i= 0; i<=ny-1; i++)
    y[i] = 0.0;

  transform Composite(
    tile = [('i',T1_I,'ii'),('j',T1_J,'jj'),('k',T1_K,'kk'),
            (('ii','i'),T2_I,'iii'),(('jj','j'),T2_J,'jjj'),(('kk','k'),T2_K,'kkk')],
    arrcopy = [(ACOPY_y, 'y[k]', [(T1_K if T1_K>1 else T2_K)],'_copy'),
               (ACOPY_x, 'x[j]', [(T1_J if T1_J>1 else T2_J)],'_copy')],
    unrolljam = (['k','j','i'],[U_K,U_J,U_I]),
    scalarreplace = (SCR, 'double'),
    regtile = (['i','j','k'],[RT_I,RT_J,RT_K]),
    vector = (VEC2, ['ivdep','vector always'])
  )
  for (i = 0; i<=nx-1; i++) {
    tmp[i] = 0;
    for (j = 0; j<=ny-1; j++)
      tmp[i] = tmp[i] + A[i*ny+j]*x[j];
    for (k = 0; k<=ny-1; k++)
      y[k] = y[k] + A[i*ny+k]*tmp[i];
  }
) @*/
{
  int i;
  for (i=0; i<=ny-14; i=i+14) {
    y[i]=0.0;
    y[(i+1)]=0.0;
    y[(i+2)]=0.0;
    y[(i+3)]=0.0;
    y[(i+4)]=0.0;
    y[(i+5)]=0.0;
    y[(i+6)]=0.0;
    y[(i+7)]=0.0;
    y[(i+8)]=0.0;
    y[(i+9)]=0.0;
    y[(i+10)]=0.0;
    y[(i+11)]=0.0;
    y[(i+12)]=0.0;
    y[(i+13)]=0.0;
  }
  for (i=ny-((ny-(0))%14); i<=ny-1; i=i+1) 
    y[i]=0.0;
}
{
  double x_copy;
  double y_copy[64];
  for (iii=0; iii<=nx-1; iii=iii+2048) 
    for (ii=iii; ii<=min(nx-1,iii+2016); ii=ii+32) {
      int i;
      for (i=ii; i<=min(nx-1,ii+31)-21; i=i+22) {
        tmp[i]=0;
        tmp[(i+1)]=0;
        tmp[(i+2)]=0;
        tmp[(i+3)]=0;
        tmp[(i+4)]=0;
        tmp[(i+5)]=0;
        tmp[(i+6)]=0;
        tmp[(i+7)]=0;
        tmp[(i+8)]=0;
        tmp[(i+9)]=0;
        tmp[(i+10)]=0;
        tmp[(i+11)]=0;
        tmp[(i+12)]=0;
        tmp[(i+13)]=0;
        tmp[(i+14)]=0;
        tmp[(i+15)]=0;
        tmp[(i+16)]=0;
        tmp[(i+17)]=0;
        tmp[(i+18)]=0;
        tmp[(i+19)]=0;
        tmp[(i+20)]=0;
        tmp[(i+21)]=0;
        for (jjj=0; jjj<=ny-1; jjj=jjj+2048) {
          for (jj=jjj; jj<=min(ny-1,jjj+1536); jj=jj+512) {
            register int cbv_1;
            cbv_1=min(ny-1,jj+511);
#pragma ivdep
#pragma vector always
            for (j=jj; j<=cbv_1; j=j+1) {
              tmp[i]=tmp[i]+A[i*ny+j]*x[j];
              tmp[(i+1)]=tmp[(i+1)]+A[(i+1)*ny+j]*x[j];
              tmp[(i+2)]=tmp[(i+2)]+A[(i+2)*ny+j]*x[j];
              tmp[(i+3)]=tmp[(i+3)]+A[(i+3)*ny+j]*x[j];
              tmp[(i+4)]=tmp[(i+4)]+A[(i+4)*ny+j]*x[j];
              tmp[(i+5)]=tmp[(i+5)]+A[(i+5)*ny+j]*x[j];
              tmp[(i+6)]=tmp[(i+6)]+A[(i+6)*ny+j]*x[j];
              tmp[(i+7)]=tmp[(i+7)]+A[(i+7)*ny+j]*x[j];
              tmp[(i+8)]=tmp[(i+8)]+A[(i+8)*ny+j]*x[j];
              tmp[(i+9)]=tmp[(i+9)]+A[(i+9)*ny+j]*x[j];
              tmp[(i+10)]=tmp[(i+10)]+A[(i+10)*ny+j]*x[j];
              tmp[(i+11)]=tmp[(i+11)]+A[(i+11)*ny+j]*x[j];
              tmp[(i+12)]=tmp[(i+12)]+A[(i+12)*ny+j]*x[j];
              tmp[(i+13)]=tmp[(i+13)]+A[(i+13)*ny+j]*x[j];
              tmp[(i+14)]=tmp[(i+14)]+A[(i+14)*ny+j]*x[j];
              tmp[(i+15)]=tmp[(i+15)]+A[(i+15)*ny+j]*x[j];
              tmp[(i+16)]=tmp[(i+16)]+A[(i+16)*ny+j]*x[j];
              tmp[(i+17)]=tmp[(i+17)]+A[(i+17)*ny+j]*x[j];
              tmp[(i+18)]=tmp[(i+18)]+A[(i+18)*ny+j]*x[j];
              tmp[(i+19)]=tmp[(i+19)]+A[(i+19)*ny+j]*x[j];
              tmp[(i+20)]=tmp[(i+20)]+A[(i+20)*ny+j]*x[j];
              tmp[(i+21)]=tmp[(i+21)]+A[(i+21)*ny+j]*x[j];
            }
          }
        }
        for (kkk=0; kkk<=ny-1; kkk=kkk+2048) {
          for (kk=kkk; kk<=min(ny-1,kkk+1984); kk=kk+64) {
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y_copy[(k-kk)]=y[k];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y_copy[(k-kk)]=y[k];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y_copy[(k-kk)]=y[k];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y_copy[(k-kk)]=y[k];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y_copy[(k-kk)]=y[k];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y_copy[(k-kk)]=y[k];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y_copy[(k-kk)]=y[k];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y_copy[(k-kk)]=y[k];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y_copy[(k-kk)]=y[k];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y_copy[(k-kk)]=y[k];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y_copy[(k-kk)]=y[k];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y_copy[(k-kk)]=y[k];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y_copy[(k-kk)]=y[k];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y_copy[(k-kk)]=y[k];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y_copy[(k-kk)]=y[k];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y_copy[(k-kk)]=y[k];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y_copy[(k-kk)]=y[k];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y_copy[(k-kk)]=y[k];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y_copy[(k-kk)]=y[k];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y_copy[(k-kk)]=y[k];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y_copy[(k-kk)]=y[k];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y_copy[(k-kk)]=y[k];
            register int cbv_2;
            cbv_2=min(ny-1,kk+63)-31;
#pragma ivdep
#pragma vector always
            for (kt=kk; kt<=cbv_2; kt=kt+32) {
              y_copy[(kt-kk)]=y_copy[(kt-kk)]+A[i*ny+kt]*tmp[i];
              y_copy[(kt-kk+1)]=y_copy[(kt-kk+1)]+A[i*ny+kt+1]*tmp[i];
              y_copy[(kt-kk+2)]=y_copy[(kt-kk+2)]+A[i*ny+kt+2]*tmp[i];
              y_copy[(kt-kk+3)]=y_copy[(kt-kk+3)]+A[i*ny+kt+3]*tmp[i];
              y_copy[(kt-kk+4)]=y_copy[(kt-kk+4)]+A[i*ny+kt+4]*tmp[i];
              y_copy[(kt-kk+5)]=y_copy[(kt-kk+5)]+A[i*ny+kt+5]*tmp[i];
              y_copy[(kt-kk+6)]=y_copy[(kt-kk+6)]+A[i*ny+kt+6]*tmp[i];
              y_copy[(kt-kk+7)]=y_copy[(kt-kk+7)]+A[i*ny+kt+7]*tmp[i];
              y_copy[(kt-kk+8)]=y_copy[(kt-kk+8)]+A[i*ny+kt+8]*tmp[i];
              y_copy[(kt-kk+9)]=y_copy[(kt-kk+9)]+A[i*ny+kt+9]*tmp[i];
              y_copy[(kt-kk+10)]=y_copy[(kt-kk+10)]+A[i*ny+kt+10]*tmp[i];
              y_copy[(kt-kk+11)]=y_copy[(kt-kk+11)]+A[i*ny+kt+11]*tmp[i];
              y_copy[(kt-kk+12)]=y_copy[(kt-kk+12)]+A[i*ny+kt+12]*tmp[i];
              y_copy[(kt-kk+13)]=y_copy[(kt-kk+13)]+A[i*ny+kt+13]*tmp[i];
              y_copy[(kt-kk+14)]=y_copy[(kt-kk+14)]+A[i*ny+kt+14]*tmp[i];
              y_copy[(kt-kk+15)]=y_copy[(kt-kk+15)]+A[i*ny+kt+15]*tmp[i];
              y_copy[(kt-kk+16)]=y_copy[(kt-kk+16)]+A[i*ny+kt+16]*tmp[i];
              y_copy[(kt-kk+17)]=y_copy[(kt-kk+17)]+A[i*ny+kt+17]*tmp[i];
              y_copy[(kt-kk+18)]=y_copy[(kt-kk+18)]+A[i*ny+kt+18]*tmp[i];
              y_copy[(kt-kk+19)]=y_copy[(kt-kk+19)]+A[i*ny+kt+19]*tmp[i];
              y_copy[(kt-kk+20)]=y_copy[(kt-kk+20)]+A[i*ny+kt+20]*tmp[i];
              y_copy[(kt-kk+21)]=y_copy[(kt-kk+21)]+A[i*ny+kt+21]*tmp[i];
              y_copy[(kt-kk+22)]=y_copy[(kt-kk+22)]+A[i*ny+kt+22]*tmp[i];
              y_copy[(kt-kk+23)]=y_copy[(kt-kk+23)]+A[i*ny+kt+23]*tmp[i];
              y_copy[(kt-kk+24)]=y_copy[(kt-kk+24)]+A[i*ny+kt+24]*tmp[i];
              y_copy[(kt-kk+25)]=y_copy[(kt-kk+25)]+A[i*ny+kt+25]*tmp[i];
              y_copy[(kt-kk+26)]=y_copy[(kt-kk+26)]+A[i*ny+kt+26]*tmp[i];
              y_copy[(kt-kk+27)]=y_copy[(kt-kk+27)]+A[i*ny+kt+27]*tmp[i];
              y_copy[(kt-kk+28)]=y_copy[(kt-kk+28)]+A[i*ny+kt+28]*tmp[i];
              y_copy[(kt-kk+29)]=y_copy[(kt-kk+29)]+A[i*ny+kt+29]*tmp[i];
              y_copy[(kt-kk+30)]=y_copy[(kt-kk+30)]+A[i*ny+kt+30]*tmp[i];
              y_copy[(kt-kk+31)]=y_copy[(kt-kk+31)]+A[i*ny+kt+31]*tmp[i];
              y_copy[(kt-kk)]=y_copy[(kt-kk)]+A[(i+1)*ny+kt]*tmp[(i+1)];
              y_copy[(kt-kk+1)]=y_copy[(kt-kk+1)]+A[(i+1)*ny+kt+1]*tmp[(i+1)];
              y_copy[(kt-kk+2)]=y_copy[(kt-kk+2)]+A[(i+1)*ny+kt+2]*tmp[(i+1)];
              y_copy[(kt-kk+3)]=y_copy[(kt-kk+3)]+A[(i+1)*ny+kt+3]*tmp[(i+1)];
              y_copy[(kt-kk+4)]=y_copy[(kt-kk+4)]+A[(i+1)*ny+kt+4]*tmp[(i+1)];
              y_copy[(kt-kk+5)]=y_copy[(kt-kk+5)]+A[(i+1)*ny+kt+5]*tmp[(i+1)];
              y_copy[(kt-kk+6)]=y_copy[(kt-kk+6)]+A[(i+1)*ny+kt+6]*tmp[(i+1)];
              y_copy[(kt-kk+7)]=y_copy[(kt-kk+7)]+A[(i+1)*ny+kt+7]*tmp[(i+1)];
              y_copy[(kt-kk+8)]=y_copy[(kt-kk+8)]+A[(i+1)*ny+kt+8]*tmp[(i+1)];
              y_copy[(kt-kk+9)]=y_copy[(kt-kk+9)]+A[(i+1)*ny+kt+9]*tmp[(i+1)];
              y_copy[(kt-kk+10)]=y_copy[(kt-kk+10)]+A[(i+1)*ny+kt+10]*tmp[(i+1)];
              y_copy[(kt-kk+11)]=y_copy[(kt-kk+11)]+A[(i+1)*ny+kt+11]*tmp[(i+1)];
              y_copy[(kt-kk+12)]=y_copy[(kt-kk+12)]+A[(i+1)*ny+kt+12]*tmp[(i+1)];
              y_copy[(kt-kk+13)]=y_copy[(kt-kk+13)]+A[(i+1)*ny+kt+13]*tmp[(i+1)];
              y_copy[(kt-kk+14)]=y_copy[(kt-kk+14)]+A[(i+1)*ny+kt+14]*tmp[(i+1)];
              y_copy[(kt-kk+15)]=y_copy[(kt-kk+15)]+A[(i+1)*ny+kt+15]*tmp[(i+1)];
              y_copy[(kt-kk+16)]=y_copy[(kt-kk+16)]+A[(i+1)*ny+kt+16]*tmp[(i+1)];
              y_copy[(kt-kk+17)]=y_copy[(kt-kk+17)]+A[(i+1)*ny+kt+17]*tmp[(i+1)];
              y_copy[(kt-kk+18)]=y_copy[(kt-kk+18)]+A[(i+1)*ny+kt+18]*tmp[(i+1)];
              y_copy[(kt-kk+19)]=y_copy[(kt-kk+19)]+A[(i+1)*ny+kt+19]*tmp[(i+1)];
              y_copy[(kt-kk+20)]=y_copy[(kt-kk+20)]+A[(i+1)*ny+kt+20]*tmp[(i+1)];
              y_copy[(kt-kk+21)]=y_copy[(kt-kk+21)]+A[(i+1)*ny+kt+21]*tmp[(i+1)];
              y_copy[(kt-kk+22)]=y_copy[(kt-kk+22)]+A[(i+1)*ny+kt+22]*tmp[(i+1)];
              y_copy[(kt-kk+23)]=y_copy[(kt-kk+23)]+A[(i+1)*ny+kt+23]*tmp[(i+1)];
              y_copy[(kt-kk+24)]=y_copy[(kt-kk+24)]+A[(i+1)*ny+kt+24]*tmp[(i+1)];
              y_copy[(kt-kk+25)]=y_copy[(kt-kk+25)]+A[(i+1)*ny+kt+25]*tmp[(i+1)];
              y_copy[(kt-kk+26)]=y_copy[(kt-kk+26)]+A[(i+1)*ny+kt+26]*tmp[(i+1)];
              y_copy[(kt-kk+27)]=y_copy[(kt-kk+27)]+A[(i+1)*ny+kt+27]*tmp[(i+1)];
              y_copy[(kt-kk+28)]=y_copy[(kt-kk+28)]+A[(i+1)*ny+kt+28]*tmp[(i+1)];
              y_copy[(kt-kk+29)]=y_copy[(kt-kk+29)]+A[(i+1)*ny+kt+29]*tmp[(i+1)];
              y_copy[(kt-kk+30)]=y_copy[(kt-kk+30)]+A[(i+1)*ny+kt+30]*tmp[(i+1)];
              y_copy[(kt-kk+31)]=y_copy[(kt-kk+31)]+A[(i+1)*ny+kt+31]*tmp[(i+1)];
              y_copy[(kt-kk)]=y_copy[(kt-kk)]+A[(i+2)*ny+kt]*tmp[(i+2)];
              y_copy[(kt-kk+1)]=y_copy[(kt-kk+1)]+A[(i+2)*ny+kt+1]*tmp[(i+2)];
              y_copy[(kt-kk+2)]=y_copy[(kt-kk+2)]+A[(i+2)*ny+kt+2]*tmp[(i+2)];
              y_copy[(kt-kk+3)]=y_copy[(kt-kk+3)]+A[(i+2)*ny+kt+3]*tmp[(i+2)];
              y_copy[(kt-kk+4)]=y_copy[(kt-kk+4)]+A[(i+2)*ny+kt+4]*tmp[(i+2)];
              y_copy[(kt-kk+5)]=y_copy[(kt-kk+5)]+A[(i+2)*ny+kt+5]*tmp[(i+2)];
              y_copy[(kt-kk+6)]=y_copy[(kt-kk+6)]+A[(i+2)*ny+kt+6]*tmp[(i+2)];
              y_copy[(kt-kk+7)]=y_copy[(kt-kk+7)]+A[(i+2)*ny+kt+7]*tmp[(i+2)];
              y_copy[(kt-kk+8)]=y_copy[(kt-kk+8)]+A[(i+2)*ny+kt+8]*tmp[(i+2)];
              y_copy[(kt-kk+9)]=y_copy[(kt-kk+9)]+A[(i+2)*ny+kt+9]*tmp[(i+2)];
              y_copy[(kt-kk+10)]=y_copy[(kt-kk+10)]+A[(i+2)*ny+kt+10]*tmp[(i+2)];
              y_copy[(kt-kk+11)]=y_copy[(kt-kk+11)]+A[(i+2)*ny+kt+11]*tmp[(i+2)];
              y_copy[(kt-kk+12)]=y_copy[(kt-kk+12)]+A[(i+2)*ny+kt+12]*tmp[(i+2)];
              y_copy[(kt-kk+13)]=y_copy[(kt-kk+13)]+A[(i+2)*ny+kt+13]*tmp[(i+2)];
              y_copy[(kt-kk+14)]=y_copy[(kt-kk+14)]+A[(i+2)*ny+kt+14]*tmp[(i+2)];
              y_copy[(kt-kk+15)]=y_copy[(kt-kk+15)]+A[(i+2)*ny+kt+15]*tmp[(i+2)];
              y_copy[(kt-kk+16)]=y_copy[(kt-kk+16)]+A[(i+2)*ny+kt+16]*tmp[(i+2)];
              y_copy[(kt-kk+17)]=y_copy[(kt-kk+17)]+A[(i+2)*ny+kt+17]*tmp[(i+2)];
              y_copy[(kt-kk+18)]=y_copy[(kt-kk+18)]+A[(i+2)*ny+kt+18]*tmp[(i+2)];
              y_copy[(kt-kk+19)]=y_copy[(kt-kk+19)]+A[(i+2)*ny+kt+19]*tmp[(i+2)];
              y_copy[(kt-kk+20)]=y_copy[(kt-kk+20)]+A[(i+2)*ny+kt+20]*tmp[(i+2)];
              y_copy[(kt-kk+21)]=y_copy[(kt-kk+21)]+A[(i+2)*ny+kt+21]*tmp[(i+2)];
              y_copy[(kt-kk+22)]=y_copy[(kt-kk+22)]+A[(i+2)*ny+kt+22]*tmp[(i+2)];
              y_copy[(kt-kk+23)]=y_copy[(kt-kk+23)]+A[(i+2)*ny+kt+23]*tmp[(i+2)];
              y_copy[(kt-kk+24)]=y_copy[(kt-kk+24)]+A[(i+2)*ny+kt+24]*tmp[(i+2)];
              y_copy[(kt-kk+25)]=y_copy[(kt-kk+25)]+A[(i+2)*ny+kt+25]*tmp[(i+2)];
              y_copy[(kt-kk+26)]=y_copy[(kt-kk+26)]+A[(i+2)*ny+kt+26]*tmp[(i+2)];
              y_copy[(kt-kk+27)]=y_copy[(kt-kk+27)]+A[(i+2)*ny+kt+27]*tmp[(i+2)];
              y_copy[(kt-kk+28)]=y_copy[(kt-kk+28)]+A[(i+2)*ny+kt+28]*tmp[(i+2)];
              y_copy[(kt-kk+29)]=y_copy[(kt-kk+29)]+A[(i+2)*ny+kt+29]*tmp[(i+2)];
              y_copy[(kt-kk+30)]=y_copy[(kt-kk+30)]+A[(i+2)*ny+kt+30]*tmp[(i+2)];
              y_copy[(kt-kk+31)]=y_copy[(kt-kk+31)]+A[(i+2)*ny+kt+31]*tmp[(i+2)];
              y_copy[(kt-kk)]=y_copy[(kt-kk)]+A[(i+3)*ny+kt]*tmp[(i+3)];
              y_copy[(kt-kk+1)]=y_copy[(kt-kk+1)]+A[(i+3)*ny+kt+1]*tmp[(i+3)];
              y_copy[(kt-kk+2)]=y_copy[(kt-kk+2)]+A[(i+3)*ny+kt+2]*tmp[(i+3)];
              y_copy[(kt-kk+3)]=y_copy[(kt-kk+3)]+A[(i+3)*ny+kt+3]*tmp[(i+3)];
              y_copy[(kt-kk+4)]=y_copy[(kt-kk+4)]+A[(i+3)*ny+kt+4]*tmp[(i+3)];
              y_copy[(kt-kk+5)]=y_copy[(kt-kk+5)]+A[(i+3)*ny+kt+5]*tmp[(i+3)];
              y_copy[(kt-kk+6)]=y_copy[(kt-kk+6)]+A[(i+3)*ny+kt+6]*tmp[(i+3)];
              y_copy[(kt-kk+7)]=y_copy[(kt-kk+7)]+A[(i+3)*ny+kt+7]*tmp[(i+3)];
              y_copy[(kt-kk+8)]=y_copy[(kt-kk+8)]+A[(i+3)*ny+kt+8]*tmp[(i+3)];
              y_copy[(kt-kk+9)]=y_copy[(kt-kk+9)]+A[(i+3)*ny+kt+9]*tmp[(i+3)];
              y_copy[(kt-kk+10)]=y_copy[(kt-kk+10)]+A[(i+3)*ny+kt+10]*tmp[(i+3)];
              y_copy[(kt-kk+11)]=y_copy[(kt-kk+11)]+A[(i+3)*ny+kt+11]*tmp[(i+3)];
              y_copy[(kt-kk+12)]=y_copy[(kt-kk+12)]+A[(i+3)*ny+kt+12]*tmp[(i+3)];
              y_copy[(kt-kk+13)]=y_copy[(kt-kk+13)]+A[(i+3)*ny+kt+13]*tmp[(i+3)];
              y_copy[(kt-kk+14)]=y_copy[(kt-kk+14)]+A[(i+3)*ny+kt+14]*tmp[(i+3)];
              y_copy[(kt-kk+15)]=y_copy[(kt-kk+15)]+A[(i+3)*ny+kt+15]*tmp[(i+3)];
              y_copy[(kt-kk+16)]=y_copy[(kt-kk+16)]+A[(i+3)*ny+kt+16]*tmp[(i+3)];
              y_copy[(kt-kk+17)]=y_copy[(kt-kk+17)]+A[(i+3)*ny+kt+17]*tmp[(i+3)];
              y_copy[(kt-kk+18)]=y_copy[(kt-kk+18)]+A[(i+3)*ny+kt+18]*tmp[(i+3)];
              y_copy[(kt-kk+19)]=y_copy[(kt-kk+19)]+A[(i+3)*ny+kt+19]*tmp[(i+3)];
              y_copy[(kt-kk+20)]=y_copy[(kt-kk+20)]+A[(i+3)*ny+kt+20]*tmp[(i+3)];
              y_copy[(kt-kk+21)]=y_copy[(kt-kk+21)]+A[(i+3)*ny+kt+21]*tmp[(i+3)];
              y_copy[(kt-kk+22)]=y_copy[(kt-kk+22)]+A[(i+3)*ny+kt+22]*tmp[(i+3)];
              y_copy[(kt-kk+23)]=y_copy[(kt-kk+23)]+A[(i+3)*ny+kt+23]*tmp[(i+3)];
              y_copy[(kt-kk+24)]=y_copy[(kt-kk+24)]+A[(i+3)*ny+kt+24]*tmp[(i+3)];
              y_copy[(kt-kk+25)]=y_copy[(kt-kk+25)]+A[(i+3)*ny+kt+25]*tmp[(i+3)];
              y_copy[(kt-kk+26)]=y_copy[(kt-kk+26)]+A[(i+3)*ny+kt+26]*tmp[(i+3)];
              y_copy[(kt-kk+27)]=y_copy[(kt-kk+27)]+A[(i+3)*ny+kt+27]*tmp[(i+3)];
              y_copy[(kt-kk+28)]=y_copy[(kt-kk+28)]+A[(i+3)*ny+kt+28]*tmp[(i+3)];
              y_copy[(kt-kk+29)]=y_copy[(kt-kk+29)]+A[(i+3)*ny+kt+29]*tmp[(i+3)];
              y_copy[(kt-kk+30)]=y_copy[(kt-kk+30)]+A[(i+3)*ny+kt+30]*tmp[(i+3)];
              y_copy[(kt-kk+31)]=y_copy[(kt-kk+31)]+A[(i+3)*ny+kt+31]*tmp[(i+3)];
              y_copy[(kt-kk)]=y_copy[(kt-kk)]+A[(i+4)*ny+kt]*tmp[(i+4)];
              y_copy[(kt-kk+1)]=y_copy[(kt-kk+1)]+A[(i+4)*ny+kt+1]*tmp[(i+4)];
              y_copy[(kt-kk+2)]=y_copy[(kt-kk+2)]+A[(i+4)*ny+kt+2]*tmp[(i+4)];
              y_copy[(kt-kk+3)]=y_copy[(kt-kk+3)]+A[(i+4)*ny+kt+3]*tmp[(i+4)];
              y_copy[(kt-kk+4)]=y_copy[(kt-kk+4)]+A[(i+4)*ny+kt+4]*tmp[(i+4)];
              y_copy[(kt-kk+5)]=y_copy[(kt-kk+5)]+A[(i+4)*ny+kt+5]*tmp[(i+4)];
              y_copy[(kt-kk+6)]=y_copy[(kt-kk+6)]+A[(i+4)*ny+kt+6]*tmp[(i+4)];
              y_copy[(kt-kk+7)]=y_copy[(kt-kk+7)]+A[(i+4)*ny+kt+7]*tmp[(i+4)];
              y_copy[(kt-kk+8)]=y_copy[(kt-kk+8)]+A[(i+4)*ny+kt+8]*tmp[(i+4)];
              y_copy[(kt-kk+9)]=y_copy[(kt-kk+9)]+A[(i+4)*ny+kt+9]*tmp[(i+4)];
              y_copy[(kt-kk+10)]=y_copy[(kt-kk+10)]+A[(i+4)*ny+kt+10]*tmp[(i+4)];
              y_copy[(kt-kk+11)]=y_copy[(kt-kk+11)]+A[(i+4)*ny+kt+11]*tmp[(i+4)];
              y_copy[(kt-kk+12)]=y_copy[(kt-kk+12)]+A[(i+4)*ny+kt+12]*tmp[(i+4)];
              y_copy[(kt-kk+13)]=y_copy[(kt-kk+13)]+A[(i+4)*ny+kt+13]*tmp[(i+4)];
              y_copy[(kt-kk+14)]=y_copy[(kt-kk+14)]+A[(i+4)*ny+kt+14]*tmp[(i+4)];
              y_copy[(kt-kk+15)]=y_copy[(kt-kk+15)]+A[(i+4)*ny+kt+15]*tmp[(i+4)];
              y_copy[(kt-kk+16)]=y_copy[(kt-kk+16)]+A[(i+4)*ny+kt+16]*tmp[(i+4)];
              y_copy[(kt-kk+17)]=y_copy[(kt-kk+17)]+A[(i+4)*ny+kt+17]*tmp[(i+4)];
              y_copy[(kt-kk+18)]=y_copy[(kt-kk+18)]+A[(i+4)*ny+kt+18]*tmp[(i+4)];
              y_copy[(kt-kk+19)]=y_copy[(kt-kk+19)]+A[(i+4)*ny+kt+19]*tmp[(i+4)];
              y_copy[(kt-kk+20)]=y_copy[(kt-kk+20)]+A[(i+4)*ny+kt+20]*tmp[(i+4)];
              y_copy[(kt-kk+21)]=y_copy[(kt-kk+21)]+A[(i+4)*ny+kt+21]*tmp[(i+4)];
              y_copy[(kt-kk+22)]=y_copy[(kt-kk+22)]+A[(i+4)*ny+kt+22]*tmp[(i+4)];
              y_copy[(kt-kk+23)]=y_copy[(kt-kk+23)]+A[(i+4)*ny+kt+23]*tmp[(i+4)];
              y_copy[(kt-kk+24)]=y_copy[(kt-kk+24)]+A[(i+4)*ny+kt+24]*tmp[(i+4)];
              y_copy[(kt-kk+25)]=y_copy[(kt-kk+25)]+A[(i+4)*ny+kt+25]*tmp[(i+4)];
              y_copy[(kt-kk+26)]=y_copy[(kt-kk+26)]+A[(i+4)*ny+kt+26]*tmp[(i+4)];
              y_copy[(kt-kk+27)]=y_copy[(kt-kk+27)]+A[(i+4)*ny+kt+27]*tmp[(i+4)];
              y_copy[(kt-kk+28)]=y_copy[(kt-kk+28)]+A[(i+4)*ny+kt+28]*tmp[(i+4)];
              y_copy[(kt-kk+29)]=y_copy[(kt-kk+29)]+A[(i+4)*ny+kt+29]*tmp[(i+4)];
              y_copy[(kt-kk+30)]=y_copy[(kt-kk+30)]+A[(i+4)*ny+kt+30]*tmp[(i+4)];
              y_copy[(kt-kk+31)]=y_copy[(kt-kk+31)]+A[(i+4)*ny+kt+31]*tmp[(i+4)];
              y_copy[(kt-kk)]=y_copy[(kt-kk)]+A[(i+5)*ny+kt]*tmp[(i+5)];
              y_copy[(kt-kk+1)]=y_copy[(kt-kk+1)]+A[(i+5)*ny+kt+1]*tmp[(i+5)];
              y_copy[(kt-kk+2)]=y_copy[(kt-kk+2)]+A[(i+5)*ny+kt+2]*tmp[(i+5)];
              y_copy[(kt-kk+3)]=y_copy[(kt-kk+3)]+A[(i+5)*ny+kt+3]*tmp[(i+5)];
              y_copy[(kt-kk+4)]=y_copy[(kt-kk+4)]+A[(i+5)*ny+kt+4]*tmp[(i+5)];
              y_copy[(kt-kk+5)]=y_copy[(kt-kk+5)]+A[(i+5)*ny+kt+5]*tmp[(i+5)];
              y_copy[(kt-kk+6)]=y_copy[(kt-kk+6)]+A[(i+5)*ny+kt+6]*tmp[(i+5)];
              y_copy[(kt-kk+7)]=y_copy[(kt-kk+7)]+A[(i+5)*ny+kt+7]*tmp[(i+5)];
              y_copy[(kt-kk+8)]=y_copy[(kt-kk+8)]+A[(i+5)*ny+kt+8]*tmp[(i+5)];
              y_copy[(kt-kk+9)]=y_copy[(kt-kk+9)]+A[(i+5)*ny+kt+9]*tmp[(i+5)];
              y_copy[(kt-kk+10)]=y_copy[(kt-kk+10)]+A[(i+5)*ny+kt+10]*tmp[(i+5)];
              y_copy[(kt-kk+11)]=y_copy[(kt-kk+11)]+A[(i+5)*ny+kt+11]*tmp[(i+5)];
              y_copy[(kt-kk+12)]=y_copy[(kt-kk+12)]+A[(i+5)*ny+kt+12]*tmp[(i+5)];
              y_copy[(kt-kk+13)]=y_copy[(kt-kk+13)]+A[(i+5)*ny+kt+13]*tmp[(i+5)];
              y_copy[(kt-kk+14)]=y_copy[(kt-kk+14)]+A[(i+5)*ny+kt+14]*tmp[(i+5)];
              y_copy[(kt-kk+15)]=y_copy[(kt-kk+15)]+A[(i+5)*ny+kt+15]*tmp[(i+5)];
              y_copy[(kt-kk+16)]=y_copy[(kt-kk+16)]+A[(i+5)*ny+kt+16]*tmp[(i+5)];
              y_copy[(kt-kk+17)]=y_copy[(kt-kk+17)]+A[(i+5)*ny+kt+17]*tmp[(i+5)];
              y_copy[(kt-kk+18)]=y_copy[(kt-kk+18)]+A[(i+5)*ny+kt+18]*tmp[(i+5)];
              y_copy[(kt-kk+19)]=y_copy[(kt-kk+19)]+A[(i+5)*ny+kt+19]*tmp[(i+5)];
              y_copy[(kt-kk+20)]=y_copy[(kt-kk+20)]+A[(i+5)*ny+kt+20]*tmp[(i+5)];
              y_copy[(kt-kk+21)]=y_copy[(kt-kk+21)]+A[(i+5)*ny+kt+21]*tmp[(i+5)];
              y_copy[(kt-kk+22)]=y_copy[(kt-kk+22)]+A[(i+5)*ny+kt+22]*tmp[(i+5)];
              y_copy[(kt-kk+23)]=y_copy[(kt-kk+23)]+A[(i+5)*ny+kt+23]*tmp[(i+5)];
              y_copy[(kt-kk+24)]=y_copy[(kt-kk+24)]+A[(i+5)*ny+kt+24]*tmp[(i+5)];
              y_copy[(kt-kk+25)]=y_copy[(kt-kk+25)]+A[(i+5)*ny+kt+25]*tmp[(i+5)];
              y_copy[(kt-kk+26)]=y_copy[(kt-kk+26)]+A[(i+5)*ny+kt+26]*tmp[(i+5)];
              y_copy[(kt-kk+27)]=y_copy[(kt-kk+27)]+A[(i+5)*ny+kt+27]*tmp[(i+5)];
              y_copy[(kt-kk+28)]=y_copy[(kt-kk+28)]+A[(i+5)*ny+kt+28]*tmp[(i+5)];
              y_copy[(kt-kk+29)]=y_copy[(kt-kk+29)]+A[(i+5)*ny+kt+29]*tmp[(i+5)];
              y_copy[(kt-kk+30)]=y_copy[(kt-kk+30)]+A[(i+5)*ny+kt+30]*tmp[(i+5)];
              y_copy[(kt-kk+31)]=y_copy[(kt-kk+31)]+A[(i+5)*ny+kt+31]*tmp[(i+5)];
              y_copy[(kt-kk)]=y_copy[(kt-kk)]+A[(i+6)*ny+kt]*tmp[(i+6)];
              y_copy[(kt-kk+1)]=y_copy[(kt-kk+1)]+A[(i+6)*ny+kt+1]*tmp[(i+6)];
              y_copy[(kt-kk+2)]=y_copy[(kt-kk+2)]+A[(i+6)*ny+kt+2]*tmp[(i+6)];
              y_copy[(kt-kk+3)]=y_copy[(kt-kk+3)]+A[(i+6)*ny+kt+3]*tmp[(i+6)];
              y_copy[(kt-kk+4)]=y_copy[(kt-kk+4)]+A[(i+6)*ny+kt+4]*tmp[(i+6)];
              y_copy[(kt-kk+5)]=y_copy[(kt-kk+5)]+A[(i+6)*ny+kt+5]*tmp[(i+6)];
              y_copy[(kt-kk+6)]=y_copy[(kt-kk+6)]+A[(i+6)*ny+kt+6]*tmp[(i+6)];
              y_copy[(kt-kk+7)]=y_copy[(kt-kk+7)]+A[(i+6)*ny+kt+7]*tmp[(i+6)];
              y_copy[(kt-kk+8)]=y_copy[(kt-kk+8)]+A[(i+6)*ny+kt+8]*tmp[(i+6)];
              y_copy[(kt-kk+9)]=y_copy[(kt-kk+9)]+A[(i+6)*ny+kt+9]*tmp[(i+6)];
              y_copy[(kt-kk+10)]=y_copy[(kt-kk+10)]+A[(i+6)*ny+kt+10]*tmp[(i+6)];
              y_copy[(kt-kk+11)]=y_copy[(kt-kk+11)]+A[(i+6)*ny+kt+11]*tmp[(i+6)];
              y_copy[(kt-kk+12)]=y_copy[(kt-kk+12)]+A[(i+6)*ny+kt+12]*tmp[(i+6)];
              y_copy[(kt-kk+13)]=y_copy[(kt-kk+13)]+A[(i+6)*ny+kt+13]*tmp[(i+6)];
              y_copy[(kt-kk+14)]=y_copy[(kt-kk+14)]+A[(i+6)*ny+kt+14]*tmp[(i+6)];
              y_copy[(kt-kk+15)]=y_copy[(kt-kk+15)]+A[(i+6)*ny+kt+15]*tmp[(i+6)];
              y_copy[(kt-kk+16)]=y_copy[(kt-kk+16)]+A[(i+6)*ny+kt+16]*tmp[(i+6)];
              y_copy[(kt-kk+17)]=y_copy[(kt-kk+17)]+A[(i+6)*ny+kt+17]*tmp[(i+6)];
              y_copy[(kt-kk+18)]=y_copy[(kt-kk+18)]+A[(i+6)*ny+kt+18]*tmp[(i+6)];
              y_copy[(kt-kk+19)]=y_copy[(kt-kk+19)]+A[(i+6)*ny+kt+19]*tmp[(i+6)];
              y_copy[(kt-kk+20)]=y_copy[(kt-kk+20)]+A[(i+6)*ny+kt+20]*tmp[(i+6)];
              y_copy[(kt-kk+21)]=y_copy[(kt-kk+21)]+A[(i+6)*ny+kt+21]*tmp[(i+6)];
              y_copy[(kt-kk+22)]=y_copy[(kt-kk+22)]+A[(i+6)*ny+kt+22]*tmp[(i+6)];
              y_copy[(kt-kk+23)]=y_copy[(kt-kk+23)]+A[(i+6)*ny+kt+23]*tmp[(i+6)];
              y_copy[(kt-kk+24)]=y_copy[(kt-kk+24)]+A[(i+6)*ny+kt+24]*tmp[(i+6)];
              y_copy[(kt-kk+25)]=y_copy[(kt-kk+25)]+A[(i+6)*ny+kt+25]*tmp[(i+6)];
              y_copy[(kt-kk+26)]=y_copy[(kt-kk+26)]+A[(i+6)*ny+kt+26]*tmp[(i+6)];
              y_copy[(kt-kk+27)]=y_copy[(kt-kk+27)]+A[(i+6)*ny+kt+27]*tmp[(i+6)];
              y_copy[(kt-kk+28)]=y_copy[(kt-kk+28)]+A[(i+6)*ny+kt+28]*tmp[(i+6)];
              y_copy[(kt-kk+29)]=y_copy[(kt-kk+29)]+A[(i+6)*ny+kt+29]*tmp[(i+6)];
              y_copy[(kt-kk+30)]=y_copy[(kt-kk+30)]+A[(i+6)*ny+kt+30]*tmp[(i+6)];
              y_copy[(kt-kk+31)]=y_copy[(kt-kk+31)]+A[(i+6)*ny+kt+31]*tmp[(i+6)];
              y_copy[(kt-kk)]=y_copy[(kt-kk)]+A[(i+7)*ny+kt]*tmp[(i+7)];
              y_copy[(kt-kk+1)]=y_copy[(kt-kk+1)]+A[(i+7)*ny+kt+1]*tmp[(i+7)];
              y_copy[(kt-kk+2)]=y_copy[(kt-kk+2)]+A[(i+7)*ny+kt+2]*tmp[(i+7)];
              y_copy[(kt-kk+3)]=y_copy[(kt-kk+3)]+A[(i+7)*ny+kt+3]*tmp[(i+7)];
              y_copy[(kt-kk+4)]=y_copy[(kt-kk+4)]+A[(i+7)*ny+kt+4]*tmp[(i+7)];
              y_copy[(kt-kk+5)]=y_copy[(kt-kk+5)]+A[(i+7)*ny+kt+5]*tmp[(i+7)];
              y_copy[(kt-kk+6)]=y_copy[(kt-kk+6)]+A[(i+7)*ny+kt+6]*tmp[(i+7)];
              y_copy[(kt-kk+7)]=y_copy[(kt-kk+7)]+A[(i+7)*ny+kt+7]*tmp[(i+7)];
              y_copy[(kt-kk+8)]=y_copy[(kt-kk+8)]+A[(i+7)*ny+kt+8]*tmp[(i+7)];
              y_copy[(kt-kk+9)]=y_copy[(kt-kk+9)]+A[(i+7)*ny+kt+9]*tmp[(i+7)];
              y_copy[(kt-kk+10)]=y_copy[(kt-kk+10)]+A[(i+7)*ny+kt+10]*tmp[(i+7)];
              y_copy[(kt-kk+11)]=y_copy[(kt-kk+11)]+A[(i+7)*ny+kt+11]*tmp[(i+7)];
              y_copy[(kt-kk+12)]=y_copy[(kt-kk+12)]+A[(i+7)*ny+kt+12]*tmp[(i+7)];
              y_copy[(kt-kk+13)]=y_copy[(kt-kk+13)]+A[(i+7)*ny+kt+13]*tmp[(i+7)];
              y_copy[(kt-kk+14)]=y_copy[(kt-kk+14)]+A[(i+7)*ny+kt+14]*tmp[(i+7)];
              y_copy[(kt-kk+15)]=y_copy[(kt-kk+15)]+A[(i+7)*ny+kt+15]*tmp[(i+7)];
              y_copy[(kt-kk+16)]=y_copy[(kt-kk+16)]+A[(i+7)*ny+kt+16]*tmp[(i+7)];
              y_copy[(kt-kk+17)]=y_copy[(kt-kk+17)]+A[(i+7)*ny+kt+17]*tmp[(i+7)];
              y_copy[(kt-kk+18)]=y_copy[(kt-kk+18)]+A[(i+7)*ny+kt+18]*tmp[(i+7)];
              y_copy[(kt-kk+19)]=y_copy[(kt-kk+19)]+A[(i+7)*ny+kt+19]*tmp[(i+7)];
              y_copy[(kt-kk+20)]=y_copy[(kt-kk+20)]+A[(i+7)*ny+kt+20]*tmp[(i+7)];
              y_copy[(kt-kk+21)]=y_copy[(kt-kk+21)]+A[(i+7)*ny+kt+21]*tmp[(i+7)];
              y_copy[(kt-kk+22)]=y_copy[(kt-kk+22)]+A[(i+7)*ny+kt+22]*tmp[(i+7)];
              y_copy[(kt-kk+23)]=y_copy[(kt-kk+23)]+A[(i+7)*ny+kt+23]*tmp[(i+7)];
              y_copy[(kt-kk+24)]=y_copy[(kt-kk+24)]+A[(i+7)*ny+kt+24]*tmp[(i+7)];
              y_copy[(kt-kk+25)]=y_copy[(kt-kk+25)]+A[(i+7)*ny+kt+25]*tmp[(i+7)];
              y_copy[(kt-kk+26)]=y_copy[(kt-kk+26)]+A[(i+7)*ny+kt+26]*tmp[(i+7)];
              y_copy[(kt-kk+27)]=y_copy[(kt-kk+27)]+A[(i+7)*ny+kt+27]*tmp[(i+7)];
              y_copy[(kt-kk+28)]=y_copy[(kt-kk+28)]+A[(i+7)*ny+kt+28]*tmp[(i+7)];
              y_copy[(kt-kk+29)]=y_copy[(kt-kk+29)]+A[(i+7)*ny+kt+29]*tmp[(i+7)];
              y_copy[(kt-kk+30)]=y_copy[(kt-kk+30)]+A[(i+7)*ny+kt+30]*tmp[(i+7)];
              y_copy[(kt-kk+31)]=y_copy[(kt-kk+31)]+A[(i+7)*ny+kt+31]*tmp[(i+7)];
              y_copy[(kt-kk)]=y_copy[(kt-kk)]+A[(i+8)*ny+kt]*tmp[(i+8)];
              y_copy[(kt-kk+1)]=y_copy[(kt-kk+1)]+A[(i+8)*ny+kt+1]*tmp[(i+8)];
              y_copy[(kt-kk+2)]=y_copy[(kt-kk+2)]+A[(i+8)*ny+kt+2]*tmp[(i+8)];
              y_copy[(kt-kk+3)]=y_copy[(kt-kk+3)]+A[(i+8)*ny+kt+3]*tmp[(i+8)];
              y_copy[(kt-kk+4)]=y_copy[(kt-kk+4)]+A[(i+8)*ny+kt+4]*tmp[(i+8)];
              y_copy[(kt-kk+5)]=y_copy[(kt-kk+5)]+A[(i+8)*ny+kt+5]*tmp[(i+8)];
              y_copy[(kt-kk+6)]=y_copy[(kt-kk+6)]+A[(i+8)*ny+kt+6]*tmp[(i+8)];
              y_copy[(kt-kk+7)]=y_copy[(kt-kk+7)]+A[(i+8)*ny+kt+7]*tmp[(i+8)];
              y_copy[(kt-kk+8)]=y_copy[(kt-kk+8)]+A[(i+8)*ny+kt+8]*tmp[(i+8)];
              y_copy[(kt-kk+9)]=y_copy[(kt-kk+9)]+A[(i+8)*ny+kt+9]*tmp[(i+8)];
              y_copy[(kt-kk+10)]=y_copy[(kt-kk+10)]+A[(i+8)*ny+kt+10]*tmp[(i+8)];
              y_copy[(kt-kk+11)]=y_copy[(kt-kk+11)]+A[(i+8)*ny+kt+11]*tmp[(i+8)];
              y_copy[(kt-kk+12)]=y_copy[(kt-kk+12)]+A[(i+8)*ny+kt+12]*tmp[(i+8)];
              y_copy[(kt-kk+13)]=y_copy[(kt-kk+13)]+A[(i+8)*ny+kt+13]*tmp[(i+8)];
              y_copy[(kt-kk+14)]=y_copy[(kt-kk+14)]+A[(i+8)*ny+kt+14]*tmp[(i+8)];
              y_copy[(kt-kk+15)]=y_copy[(kt-kk+15)]+A[(i+8)*ny+kt+15]*tmp[(i+8)];
              y_copy[(kt-kk+16)]=y_copy[(kt-kk+16)]+A[(i+8)*ny+kt+16]*tmp[(i+8)];
              y_copy[(kt-kk+17)]=y_copy[(kt-kk+17)]+A[(i+8)*ny+kt+17]*tmp[(i+8)];
              y_copy[(kt-kk+18)]=y_copy[(kt-kk+18)]+A[(i+8)*ny+kt+18]*tmp[(i+8)];
              y_copy[(kt-kk+19)]=y_copy[(kt-kk+19)]+A[(i+8)*ny+kt+19]*tmp[(i+8)];
              y_copy[(kt-kk+20)]=y_copy[(kt-kk+20)]+A[(i+8)*ny+kt+20]*tmp[(i+8)];
              y_copy[(kt-kk+21)]=y_copy[(kt-kk+21)]+A[(i+8)*ny+kt+21]*tmp[(i+8)];
              y_copy[(kt-kk+22)]=y_copy[(kt-kk+22)]+A[(i+8)*ny+kt+22]*tmp[(i+8)];
              y_copy[(kt-kk+23)]=y_copy[(kt-kk+23)]+A[(i+8)*ny+kt+23]*tmp[(i+8)];
              y_copy[(kt-kk+24)]=y_copy[(kt-kk+24)]+A[(i+8)*ny+kt+24]*tmp[(i+8)];
              y_copy[(kt-kk+25)]=y_copy[(kt-kk+25)]+A[(i+8)*ny+kt+25]*tmp[(i+8)];
              y_copy[(kt-kk+26)]=y_copy[(kt-kk+26)]+A[(i+8)*ny+kt+26]*tmp[(i+8)];
              y_copy[(kt-kk+27)]=y_copy[(kt-kk+27)]+A[(i+8)*ny+kt+27]*tmp[(i+8)];
              y_copy[(kt-kk+28)]=y_copy[(kt-kk+28)]+A[(i+8)*ny+kt+28]*tmp[(i+8)];
              y_copy[(kt-kk+29)]=y_copy[(kt-kk+29)]+A[(i+8)*ny+kt+29]*tmp[(i+8)];
              y_copy[(kt-kk+30)]=y_copy[(kt-kk+30)]+A[(i+8)*ny+kt+30]*tmp[(i+8)];
              y_copy[(kt-kk+31)]=y_copy[(kt-kk+31)]+A[(i+8)*ny+kt+31]*tmp[(i+8)];
              y_copy[(kt-kk)]=y_copy[(kt-kk)]+A[(i+9)*ny+kt]*tmp[(i+9)];
              y_copy[(kt-kk+1)]=y_copy[(kt-kk+1)]+A[(i+9)*ny+kt+1]*tmp[(i+9)];
              y_copy[(kt-kk+2)]=y_copy[(kt-kk+2)]+A[(i+9)*ny+kt+2]*tmp[(i+9)];
              y_copy[(kt-kk+3)]=y_copy[(kt-kk+3)]+A[(i+9)*ny+kt+3]*tmp[(i+9)];
              y_copy[(kt-kk+4)]=y_copy[(kt-kk+4)]+A[(i+9)*ny+kt+4]*tmp[(i+9)];
              y_copy[(kt-kk+5)]=y_copy[(kt-kk+5)]+A[(i+9)*ny+kt+5]*tmp[(i+9)];
              y_copy[(kt-kk+6)]=y_copy[(kt-kk+6)]+A[(i+9)*ny+kt+6]*tmp[(i+9)];
              y_copy[(kt-kk+7)]=y_copy[(kt-kk+7)]+A[(i+9)*ny+kt+7]*tmp[(i+9)];
              y_copy[(kt-kk+8)]=y_copy[(kt-kk+8)]+A[(i+9)*ny+kt+8]*tmp[(i+9)];
              y_copy[(kt-kk+9)]=y_copy[(kt-kk+9)]+A[(i+9)*ny+kt+9]*tmp[(i+9)];
              y_copy[(kt-kk+10)]=y_copy[(kt-kk+10)]+A[(i+9)*ny+kt+10]*tmp[(i+9)];
              y_copy[(kt-kk+11)]=y_copy[(kt-kk+11)]+A[(i+9)*ny+kt+11]*tmp[(i+9)];
              y_copy[(kt-kk+12)]=y_copy[(kt-kk+12)]+A[(i+9)*ny+kt+12]*tmp[(i+9)];
              y_copy[(kt-kk+13)]=y_copy[(kt-kk+13)]+A[(i+9)*ny+kt+13]*tmp[(i+9)];
              y_copy[(kt-kk+14)]=y_copy[(kt-kk+14)]+A[(i+9)*ny+kt+14]*tmp[(i+9)];
              y_copy[(kt-kk+15)]=y_copy[(kt-kk+15)]+A[(i+9)*ny+kt+15]*tmp[(i+9)];
              y_copy[(kt-kk+16)]=y_copy[(kt-kk+16)]+A[(i+9)*ny+kt+16]*tmp[(i+9)];
              y_copy[(kt-kk+17)]=y_copy[(kt-kk+17)]+A[(i+9)*ny+kt+17]*tmp[(i+9)];
              y_copy[(kt-kk+18)]=y_copy[(kt-kk+18)]+A[(i+9)*ny+kt+18]*tmp[(i+9)];
              y_copy[(kt-kk+19)]=y_copy[(kt-kk+19)]+A[(i+9)*ny+kt+19]*tmp[(i+9)];
              y_copy[(kt-kk+20)]=y_copy[(kt-kk+20)]+A[(i+9)*ny+kt+20]*tmp[(i+9)];
              y_copy[(kt-kk+21)]=y_copy[(kt-kk+21)]+A[(i+9)*ny+kt+21]*tmp[(i+9)];
              y_copy[(kt-kk+22)]=y_copy[(kt-kk+22)]+A[(i+9)*ny+kt+22]*tmp[(i+9)];
              y_copy[(kt-kk+23)]=y_copy[(kt-kk+23)]+A[(i+9)*ny+kt+23]*tmp[(i+9)];
              y_copy[(kt-kk+24)]=y_copy[(kt-kk+24)]+A[(i+9)*ny+kt+24]*tmp[(i+9)];
              y_copy[(kt-kk+25)]=y_copy[(kt-kk+25)]+A[(i+9)*ny+kt+25]*tmp[(i+9)];
              y_copy[(kt-kk+26)]=y_copy[(kt-kk+26)]+A[(i+9)*ny+kt+26]*tmp[(i+9)];
              y_copy[(kt-kk+27)]=y_copy[(kt-kk+27)]+A[(i+9)*ny+kt+27]*tmp[(i+9)];
              y_copy[(kt-kk+28)]=y_copy[(kt-kk+28)]+A[(i+9)*ny+kt+28]*tmp[(i+9)];
              y_copy[(kt-kk+29)]=y_copy[(kt-kk+29)]+A[(i+9)*ny+kt+29]*tmp[(i+9)];
              y_copy[(kt-kk+30)]=y_copy[(kt-kk+30)]+A[(i+9)*ny+kt+30]*tmp[(i+9)];
              y_copy[(kt-kk+31)]=y_copy[(kt-kk+31)]+A[(i+9)*ny+kt+31]*tmp[(i+9)];
              y_copy[(kt-kk)]=y_copy[(kt-kk)]+A[(i+10)*ny+kt]*tmp[(i+10)];
              y_copy[(kt-kk+1)]=y_copy[(kt-kk+1)]+A[(i+10)*ny+kt+1]*tmp[(i+10)];
              y_copy[(kt-kk+2)]=y_copy[(kt-kk+2)]+A[(i+10)*ny+kt+2]*tmp[(i+10)];
              y_copy[(kt-kk+3)]=y_copy[(kt-kk+3)]+A[(i+10)*ny+kt+3]*tmp[(i+10)];
              y_copy[(kt-kk+4)]=y_copy[(kt-kk+4)]+A[(i+10)*ny+kt+4]*tmp[(i+10)];
              y_copy[(kt-kk+5)]=y_copy[(kt-kk+5)]+A[(i+10)*ny+kt+5]*tmp[(i+10)];
              y_copy[(kt-kk+6)]=y_copy[(kt-kk+6)]+A[(i+10)*ny+kt+6]*tmp[(i+10)];
              y_copy[(kt-kk+7)]=y_copy[(kt-kk+7)]+A[(i+10)*ny+kt+7]*tmp[(i+10)];
              y_copy[(kt-kk+8)]=y_copy[(kt-kk+8)]+A[(i+10)*ny+kt+8]*tmp[(i+10)];
              y_copy[(kt-kk+9)]=y_copy[(kt-kk+9)]+A[(i+10)*ny+kt+9]*tmp[(i+10)];
              y_copy[(kt-kk+10)]=y_copy[(kt-kk+10)]+A[(i+10)*ny+kt+10]*tmp[(i+10)];
              y_copy[(kt-kk+11)]=y_copy[(kt-kk+11)]+A[(i+10)*ny+kt+11]*tmp[(i+10)];
              y_copy[(kt-kk+12)]=y_copy[(kt-kk+12)]+A[(i+10)*ny+kt+12]*tmp[(i+10)];
              y_copy[(kt-kk+13)]=y_copy[(kt-kk+13)]+A[(i+10)*ny+kt+13]*tmp[(i+10)];
              y_copy[(kt-kk+14)]=y_copy[(kt-kk+14)]+A[(i+10)*ny+kt+14]*tmp[(i+10)];
              y_copy[(kt-kk+15)]=y_copy[(kt-kk+15)]+A[(i+10)*ny+kt+15]*tmp[(i+10)];
              y_copy[(kt-kk+16)]=y_copy[(kt-kk+16)]+A[(i+10)*ny+kt+16]*tmp[(i+10)];
              y_copy[(kt-kk+17)]=y_copy[(kt-kk+17)]+A[(i+10)*ny+kt+17]*tmp[(i+10)];
              y_copy[(kt-kk+18)]=y_copy[(kt-kk+18)]+A[(i+10)*ny+kt+18]*tmp[(i+10)];
              y_copy[(kt-kk+19)]=y_copy[(kt-kk+19)]+A[(i+10)*ny+kt+19]*tmp[(i+10)];
              y_copy[(kt-kk+20)]=y_copy[(kt-kk+20)]+A[(i+10)*ny+kt+20]*tmp[(i+10)];
              y_copy[(kt-kk+21)]=y_copy[(kt-kk+21)]+A[(i+10)*ny+kt+21]*tmp[(i+10)];
              y_copy[(kt-kk+22)]=y_copy[(kt-kk+22)]+A[(i+10)*ny+kt+22]*tmp[(i+10)];
              y_copy[(kt-kk+23)]=y_copy[(kt-kk+23)]+A[(i+10)*ny+kt+23]*tmp[(i+10)];
              y_copy[(kt-kk+24)]=y_copy[(kt-kk+24)]+A[(i+10)*ny+kt+24]*tmp[(i+10)];
              y_copy[(kt-kk+25)]=y_copy[(kt-kk+25)]+A[(i+10)*ny+kt+25]*tmp[(i+10)];
              y_copy[(kt-kk+26)]=y_copy[(kt-kk+26)]+A[(i+10)*ny+kt+26]*tmp[(i+10)];
              y_copy[(kt-kk+27)]=y_copy[(kt-kk+27)]+A[(i+10)*ny+kt+27]*tmp[(i+10)];
              y_copy[(kt-kk+28)]=y_copy[(kt-kk+28)]+A[(i+10)*ny+kt+28]*tmp[(i+10)];
              y_copy[(kt-kk+29)]=y_copy[(kt-kk+29)]+A[(i+10)*ny+kt+29]*tmp[(i+10)];
              y_copy[(kt-kk+30)]=y_copy[(kt-kk+30)]+A[(i+10)*ny+kt+30]*tmp[(i+10)];
              y_copy[(kt-kk+31)]=y_copy[(kt-kk+31)]+A[(i+10)*ny+kt+31]*tmp[(i+10)];
              y_copy[(kt-kk)]=y_copy[(kt-kk)]+A[(i+11)*ny+kt]*tmp[(i+11)];
              y_copy[(kt-kk+1)]=y_copy[(kt-kk+1)]+A[(i+11)*ny+kt+1]*tmp[(i+11)];
              y_copy[(kt-kk+2)]=y_copy[(kt-kk+2)]+A[(i+11)*ny+kt+2]*tmp[(i+11)];
              y_copy[(kt-kk+3)]=y_copy[(kt-kk+3)]+A[(i+11)*ny+kt+3]*tmp[(i+11)];
              y_copy[(kt-kk+4)]=y_copy[(kt-kk+4)]+A[(i+11)*ny+kt+4]*tmp[(i+11)];
              y_copy[(kt-kk+5)]=y_copy[(kt-kk+5)]+A[(i+11)*ny+kt+5]*tmp[(i+11)];
              y_copy[(kt-kk+6)]=y_copy[(kt-kk+6)]+A[(i+11)*ny+kt+6]*tmp[(i+11)];
              y_copy[(kt-kk+7)]=y_copy[(kt-kk+7)]+A[(i+11)*ny+kt+7]*tmp[(i+11)];
              y_copy[(kt-kk+8)]=y_copy[(kt-kk+8)]+A[(i+11)*ny+kt+8]*tmp[(i+11)];
              y_copy[(kt-kk+9)]=y_copy[(kt-kk+9)]+A[(i+11)*ny+kt+9]*tmp[(i+11)];
              y_copy[(kt-kk+10)]=y_copy[(kt-kk+10)]+A[(i+11)*ny+kt+10]*tmp[(i+11)];
              y_copy[(kt-kk+11)]=y_copy[(kt-kk+11)]+A[(i+11)*ny+kt+11]*tmp[(i+11)];
              y_copy[(kt-kk+12)]=y_copy[(kt-kk+12)]+A[(i+11)*ny+kt+12]*tmp[(i+11)];
              y_copy[(kt-kk+13)]=y_copy[(kt-kk+13)]+A[(i+11)*ny+kt+13]*tmp[(i+11)];
              y_copy[(kt-kk+14)]=y_copy[(kt-kk+14)]+A[(i+11)*ny+kt+14]*tmp[(i+11)];
              y_copy[(kt-kk+15)]=y_copy[(kt-kk+15)]+A[(i+11)*ny+kt+15]*tmp[(i+11)];
              y_copy[(kt-kk+16)]=y_copy[(kt-kk+16)]+A[(i+11)*ny+kt+16]*tmp[(i+11)];
              y_copy[(kt-kk+17)]=y_copy[(kt-kk+17)]+A[(i+11)*ny+kt+17]*tmp[(i+11)];
              y_copy[(kt-kk+18)]=y_copy[(kt-kk+18)]+A[(i+11)*ny+kt+18]*tmp[(i+11)];
              y_copy[(kt-kk+19)]=y_copy[(kt-kk+19)]+A[(i+11)*ny+kt+19]*tmp[(i+11)];
              y_copy[(kt-kk+20)]=y_copy[(kt-kk+20)]+A[(i+11)*ny+kt+20]*tmp[(i+11)];
              y_copy[(kt-kk+21)]=y_copy[(kt-kk+21)]+A[(i+11)*ny+kt+21]*tmp[(i+11)];
              y_copy[(kt-kk+22)]=y_copy[(kt-kk+22)]+A[(i+11)*ny+kt+22]*tmp[(i+11)];
              y_copy[(kt-kk+23)]=y_copy[(kt-kk+23)]+A[(i+11)*ny+kt+23]*tmp[(i+11)];
              y_copy[(kt-kk+24)]=y_copy[(kt-kk+24)]+A[(i+11)*ny+kt+24]*tmp[(i+11)];
              y_copy[(kt-kk+25)]=y_copy[(kt-kk+25)]+A[(i+11)*ny+kt+25]*tmp[(i+11)];
              y_copy[(kt-kk+26)]=y_copy[(kt-kk+26)]+A[(i+11)*ny+kt+26]*tmp[(i+11)];
              y_copy[(kt-kk+27)]=y_copy[(kt-kk+27)]+A[(i+11)*ny+kt+27]*tmp[(i+11)];
              y_copy[(kt-kk+28)]=y_copy[(kt-kk+28)]+A[(i+11)*ny+kt+28]*tmp[(i+11)];
              y_copy[(kt-kk+29)]=y_copy[(kt-kk+29)]+A[(i+11)*ny+kt+29]*tmp[(i+11)];
              y_copy[(kt-kk+30)]=y_copy[(kt-kk+30)]+A[(i+11)*ny+kt+30]*tmp[(i+11)];
              y_copy[(kt-kk+31)]=y_copy[(kt-kk+31)]+A[(i+11)*ny+kt+31]*tmp[(i+11)];
              y_copy[(kt-kk)]=y_copy[(kt-kk)]+A[(i+12)*ny+kt]*tmp[(i+12)];
              y_copy[(kt-kk+1)]=y_copy[(kt-kk+1)]+A[(i+12)*ny+kt+1]*tmp[(i+12)];
              y_copy[(kt-kk+2)]=y_copy[(kt-kk+2)]+A[(i+12)*ny+kt+2]*tmp[(i+12)];
              y_copy[(kt-kk+3)]=y_copy[(kt-kk+3)]+A[(i+12)*ny+kt+3]*tmp[(i+12)];
              y_copy[(kt-kk+4)]=y_copy[(kt-kk+4)]+A[(i+12)*ny+kt+4]*tmp[(i+12)];
              y_copy[(kt-kk+5)]=y_copy[(kt-kk+5)]+A[(i+12)*ny+kt+5]*tmp[(i+12)];
              y_copy[(kt-kk+6)]=y_copy[(kt-kk+6)]+A[(i+12)*ny+kt+6]*tmp[(i+12)];
              y_copy[(kt-kk+7)]=y_copy[(kt-kk+7)]+A[(i+12)*ny+kt+7]*tmp[(i+12)];
              y_copy[(kt-kk+8)]=y_copy[(kt-kk+8)]+A[(i+12)*ny+kt+8]*tmp[(i+12)];
              y_copy[(kt-kk+9)]=y_copy[(kt-kk+9)]+A[(i+12)*ny+kt+9]*tmp[(i+12)];
              y_copy[(kt-kk+10)]=y_copy[(kt-kk+10)]+A[(i+12)*ny+kt+10]*tmp[(i+12)];
              y_copy[(kt-kk+11)]=y_copy[(kt-kk+11)]+A[(i+12)*ny+kt+11]*tmp[(i+12)];
              y_copy[(kt-kk+12)]=y_copy[(kt-kk+12)]+A[(i+12)*ny+kt+12]*tmp[(i+12)];
              y_copy[(kt-kk+13)]=y_copy[(kt-kk+13)]+A[(i+12)*ny+kt+13]*tmp[(i+12)];
              y_copy[(kt-kk+14)]=y_copy[(kt-kk+14)]+A[(i+12)*ny+kt+14]*tmp[(i+12)];
              y_copy[(kt-kk+15)]=y_copy[(kt-kk+15)]+A[(i+12)*ny+kt+15]*tmp[(i+12)];
              y_copy[(kt-kk+16)]=y_copy[(kt-kk+16)]+A[(i+12)*ny+kt+16]*tmp[(i+12)];
              y_copy[(kt-kk+17)]=y_copy[(kt-kk+17)]+A[(i+12)*ny+kt+17]*tmp[(i+12)];
              y_copy[(kt-kk+18)]=y_copy[(kt-kk+18)]+A[(i+12)*ny+kt+18]*tmp[(i+12)];
              y_copy[(kt-kk+19)]=y_copy[(kt-kk+19)]+A[(i+12)*ny+kt+19]*tmp[(i+12)];
              y_copy[(kt-kk+20)]=y_copy[(kt-kk+20)]+A[(i+12)*ny+kt+20]*tmp[(i+12)];
              y_copy[(kt-kk+21)]=y_copy[(kt-kk+21)]+A[(i+12)*ny+kt+21]*tmp[(i+12)];
              y_copy[(kt-kk+22)]=y_copy[(kt-kk+22)]+A[(i+12)*ny+kt+22]*tmp[(i+12)];
              y_copy[(kt-kk+23)]=y_copy[(kt-kk+23)]+A[(i+12)*ny+kt+23]*tmp[(i+12)];
              y_copy[(kt-kk+24)]=y_copy[(kt-kk+24)]+A[(i+12)*ny+kt+24]*tmp[(i+12)];
              y_copy[(kt-kk+25)]=y_copy[(kt-kk+25)]+A[(i+12)*ny+kt+25]*tmp[(i+12)];
              y_copy[(kt-kk+26)]=y_copy[(kt-kk+26)]+A[(i+12)*ny+kt+26]*tmp[(i+12)];
              y_copy[(kt-kk+27)]=y_copy[(kt-kk+27)]+A[(i+12)*ny+kt+27]*tmp[(i+12)];
              y_copy[(kt-kk+28)]=y_copy[(kt-kk+28)]+A[(i+12)*ny+kt+28]*tmp[(i+12)];
              y_copy[(kt-kk+29)]=y_copy[(kt-kk+29)]+A[(i+12)*ny+kt+29]*tmp[(i+12)];
              y_copy[(kt-kk+30)]=y_copy[(kt-kk+30)]+A[(i+12)*ny+kt+30]*tmp[(i+12)];
              y_copy[(kt-kk+31)]=y_copy[(kt-kk+31)]+A[(i+12)*ny+kt+31]*tmp[(i+12)];
              y_copy[(kt-kk)]=y_copy[(kt-kk)]+A[(i+13)*ny+kt]*tmp[(i+13)];
              y_copy[(kt-kk+1)]=y_copy[(kt-kk+1)]+A[(i+13)*ny+kt+1]*tmp[(i+13)];
              y_copy[(kt-kk+2)]=y_copy[(kt-kk+2)]+A[(i+13)*ny+kt+2]*tmp[(i+13)];
              y_copy[(kt-kk+3)]=y_copy[(kt-kk+3)]+A[(i+13)*ny+kt+3]*tmp[(i+13)];
              y_copy[(kt-kk+4)]=y_copy[(kt-kk+4)]+A[(i+13)*ny+kt+4]*tmp[(i+13)];
              y_copy[(kt-kk+5)]=y_copy[(kt-kk+5)]+A[(i+13)*ny+kt+5]*tmp[(i+13)];
              y_copy[(kt-kk+6)]=y_copy[(kt-kk+6)]+A[(i+13)*ny+kt+6]*tmp[(i+13)];
              y_copy[(kt-kk+7)]=y_copy[(kt-kk+7)]+A[(i+13)*ny+kt+7]*tmp[(i+13)];
              y_copy[(kt-kk+8)]=y_copy[(kt-kk+8)]+A[(i+13)*ny+kt+8]*tmp[(i+13)];
              y_copy[(kt-kk+9)]=y_copy[(kt-kk+9)]+A[(i+13)*ny+kt+9]*tmp[(i+13)];
              y_copy[(kt-kk+10)]=y_copy[(kt-kk+10)]+A[(i+13)*ny+kt+10]*tmp[(i+13)];
              y_copy[(kt-kk+11)]=y_copy[(kt-kk+11)]+A[(i+13)*ny+kt+11]*tmp[(i+13)];
              y_copy[(kt-kk+12)]=y_copy[(kt-kk+12)]+A[(i+13)*ny+kt+12]*tmp[(i+13)];
              y_copy[(kt-kk+13)]=y_copy[(kt-kk+13)]+A[(i+13)*ny+kt+13]*tmp[(i+13)];
              y_copy[(kt-kk+14)]=y_copy[(kt-kk+14)]+A[(i+13)*ny+kt+14]*tmp[(i+13)];
              y_copy[(kt-kk+15)]=y_copy[(kt-kk+15)]+A[(i+13)*ny+kt+15]*tmp[(i+13)];
              y_copy[(kt-kk+16)]=y_copy[(kt-kk+16)]+A[(i+13)*ny+kt+16]*tmp[(i+13)];
              y_copy[(kt-kk+17)]=y_copy[(kt-kk+17)]+A[(i+13)*ny+kt+17]*tmp[(i+13)];
              y_copy[(kt-kk+18)]=y_copy[(kt-kk+18)]+A[(i+13)*ny+kt+18]*tmp[(i+13)];
              y_copy[(kt-kk+19)]=y_copy[(kt-kk+19)]+A[(i+13)*ny+kt+19]*tmp[(i+13)];
              y_copy[(kt-kk+20)]=y_copy[(kt-kk+20)]+A[(i+13)*ny+kt+20]*tmp[(i+13)];
              y_copy[(kt-kk+21)]=y_copy[(kt-kk+21)]+A[(i+13)*ny+kt+21]*tmp[(i+13)];
              y_copy[(kt-kk+22)]=y_copy[(kt-kk+22)]+A[(i+13)*ny+kt+22]*tmp[(i+13)];
              y_copy[(kt-kk+23)]=y_copy[(kt-kk+23)]+A[(i+13)*ny+kt+23]*tmp[(i+13)];
              y_copy[(kt-kk+24)]=y_copy[(kt-kk+24)]+A[(i+13)*ny+kt+24]*tmp[(i+13)];
              y_copy[(kt-kk+25)]=y_copy[(kt-kk+25)]+A[(i+13)*ny+kt+25]*tmp[(i+13)];
              y_copy[(kt-kk+26)]=y_copy[(kt-kk+26)]+A[(i+13)*ny+kt+26]*tmp[(i+13)];
              y_copy[(kt-kk+27)]=y_copy[(kt-kk+27)]+A[(i+13)*ny+kt+27]*tmp[(i+13)];
              y_copy[(kt-kk+28)]=y_copy[(kt-kk+28)]+A[(i+13)*ny+kt+28]*tmp[(i+13)];
              y_copy[(kt-kk+29)]=y_copy[(kt-kk+29)]+A[(i+13)*ny+kt+29]*tmp[(i+13)];
              y_copy[(kt-kk+30)]=y_copy[(kt-kk+30)]+A[(i+13)*ny+kt+30]*tmp[(i+13)];
              y_copy[(kt-kk+31)]=y_copy[(kt-kk+31)]+A[(i+13)*ny+kt+31]*tmp[(i+13)];
              y_copy[(kt-kk)]=y_copy[(kt-kk)]+A[(i+14)*ny+kt]*tmp[(i+14)];
              y_copy[(kt-kk+1)]=y_copy[(kt-kk+1)]+A[(i+14)*ny+kt+1]*tmp[(i+14)];
              y_copy[(kt-kk+2)]=y_copy[(kt-kk+2)]+A[(i+14)*ny+kt+2]*tmp[(i+14)];
              y_copy[(kt-kk+3)]=y_copy[(kt-kk+3)]+A[(i+14)*ny+kt+3]*tmp[(i+14)];
              y_copy[(kt-kk+4)]=y_copy[(kt-kk+4)]+A[(i+14)*ny+kt+4]*tmp[(i+14)];
              y_copy[(kt-kk+5)]=y_copy[(kt-kk+5)]+A[(i+14)*ny+kt+5]*tmp[(i+14)];
              y_copy[(kt-kk+6)]=y_copy[(kt-kk+6)]+A[(i+14)*ny+kt+6]*tmp[(i+14)];
              y_copy[(kt-kk+7)]=y_copy[(kt-kk+7)]+A[(i+14)*ny+kt+7]*tmp[(i+14)];
              y_copy[(kt-kk+8)]=y_copy[(kt-kk+8)]+A[(i+14)*ny+kt+8]*tmp[(i+14)];
              y_copy[(kt-kk+9)]=y_copy[(kt-kk+9)]+A[(i+14)*ny+kt+9]*tmp[(i+14)];
              y_copy[(kt-kk+10)]=y_copy[(kt-kk+10)]+A[(i+14)*ny+kt+10]*tmp[(i+14)];
              y_copy[(kt-kk+11)]=y_copy[(kt-kk+11)]+A[(i+14)*ny+kt+11]*tmp[(i+14)];
              y_copy[(kt-kk+12)]=y_copy[(kt-kk+12)]+A[(i+14)*ny+kt+12]*tmp[(i+14)];
              y_copy[(kt-kk+13)]=y_copy[(kt-kk+13)]+A[(i+14)*ny+kt+13]*tmp[(i+14)];
              y_copy[(kt-kk+14)]=y_copy[(kt-kk+14)]+A[(i+14)*ny+kt+14]*tmp[(i+14)];
              y_copy[(kt-kk+15)]=y_copy[(kt-kk+15)]+A[(i+14)*ny+kt+15]*tmp[(i+14)];
              y_copy[(kt-kk+16)]=y_copy[(kt-kk+16)]+A[(i+14)*ny+kt+16]*tmp[(i+14)];
              y_copy[(kt-kk+17)]=y_copy[(kt-kk+17)]+A[(i+14)*ny+kt+17]*tmp[(i+14)];
              y_copy[(kt-kk+18)]=y_copy[(kt-kk+18)]+A[(i+14)*ny+kt+18]*tmp[(i+14)];
              y_copy[(kt-kk+19)]=y_copy[(kt-kk+19)]+A[(i+14)*ny+kt+19]*tmp[(i+14)];
              y_copy[(kt-kk+20)]=y_copy[(kt-kk+20)]+A[(i+14)*ny+kt+20]*tmp[(i+14)];
              y_copy[(kt-kk+21)]=y_copy[(kt-kk+21)]+A[(i+14)*ny+kt+21]*tmp[(i+14)];
              y_copy[(kt-kk+22)]=y_copy[(kt-kk+22)]+A[(i+14)*ny+kt+22]*tmp[(i+14)];
              y_copy[(kt-kk+23)]=y_copy[(kt-kk+23)]+A[(i+14)*ny+kt+23]*tmp[(i+14)];
              y_copy[(kt-kk+24)]=y_copy[(kt-kk+24)]+A[(i+14)*ny+kt+24]*tmp[(i+14)];
              y_copy[(kt-kk+25)]=y_copy[(kt-kk+25)]+A[(i+14)*ny+kt+25]*tmp[(i+14)];
              y_copy[(kt-kk+26)]=y_copy[(kt-kk+26)]+A[(i+14)*ny+kt+26]*tmp[(i+14)];
              y_copy[(kt-kk+27)]=y_copy[(kt-kk+27)]+A[(i+14)*ny+kt+27]*tmp[(i+14)];
              y_copy[(kt-kk+28)]=y_copy[(kt-kk+28)]+A[(i+14)*ny+kt+28]*tmp[(i+14)];
              y_copy[(kt-kk+29)]=y_copy[(kt-kk+29)]+A[(i+14)*ny+kt+29]*tmp[(i+14)];
              y_copy[(kt-kk+30)]=y_copy[(kt-kk+30)]+A[(i+14)*ny+kt+30]*tmp[(i+14)];
              y_copy[(kt-kk+31)]=y_copy[(kt-kk+31)]+A[(i+14)*ny+kt+31]*tmp[(i+14)];
              y_copy[(kt-kk)]=y_copy[(kt-kk)]+A[(i+15)*ny+kt]*tmp[(i+15)];
              y_copy[(kt-kk+1)]=y_copy[(kt-kk+1)]+A[(i+15)*ny+kt+1]*tmp[(i+15)];
              y_copy[(kt-kk+2)]=y_copy[(kt-kk+2)]+A[(i+15)*ny+kt+2]*tmp[(i+15)];
              y_copy[(kt-kk+3)]=y_copy[(kt-kk+3)]+A[(i+15)*ny+kt+3]*tmp[(i+15)];
              y_copy[(kt-kk+4)]=y_copy[(kt-kk+4)]+A[(i+15)*ny+kt+4]*tmp[(i+15)];
              y_copy[(kt-kk+5)]=y_copy[(kt-kk+5)]+A[(i+15)*ny+kt+5]*tmp[(i+15)];
              y_copy[(kt-kk+6)]=y_copy[(kt-kk+6)]+A[(i+15)*ny+kt+6]*tmp[(i+15)];
              y_copy[(kt-kk+7)]=y_copy[(kt-kk+7)]+A[(i+15)*ny+kt+7]*tmp[(i+15)];
              y_copy[(kt-kk+8)]=y_copy[(kt-kk+8)]+A[(i+15)*ny+kt+8]*tmp[(i+15)];
              y_copy[(kt-kk+9)]=y_copy[(kt-kk+9)]+A[(i+15)*ny+kt+9]*tmp[(i+15)];
              y_copy[(kt-kk+10)]=y_copy[(kt-kk+10)]+A[(i+15)*ny+kt+10]*tmp[(i+15)];
              y_copy[(kt-kk+11)]=y_copy[(kt-kk+11)]+A[(i+15)*ny+kt+11]*tmp[(i+15)];
              y_copy[(kt-kk+12)]=y_copy[(kt-kk+12)]+A[(i+15)*ny+kt+12]*tmp[(i+15)];
              y_copy[(kt-kk+13)]=y_copy[(kt-kk+13)]+A[(i+15)*ny+kt+13]*tmp[(i+15)];
              y_copy[(kt-kk+14)]=y_copy[(kt-kk+14)]+A[(i+15)*ny+kt+14]*tmp[(i+15)];
              y_copy[(kt-kk+15)]=y_copy[(kt-kk+15)]+A[(i+15)*ny+kt+15]*tmp[(i+15)];
              y_copy[(kt-kk+16)]=y_copy[(kt-kk+16)]+A[(i+15)*ny+kt+16]*tmp[(i+15)];
              y_copy[(kt-kk+17)]=y_copy[(kt-kk+17)]+A[(i+15)*ny+kt+17]*tmp[(i+15)];
              y_copy[(kt-kk+18)]=y_copy[(kt-kk+18)]+A[(i+15)*ny+kt+18]*tmp[(i+15)];
              y_copy[(kt-kk+19)]=y_copy[(kt-kk+19)]+A[(i+15)*ny+kt+19]*tmp[(i+15)];
              y_copy[(kt-kk+20)]=y_copy[(kt-kk+20)]+A[(i+15)*ny+kt+20]*tmp[(i+15)];
              y_copy[(kt-kk+21)]=y_copy[(kt-kk+21)]+A[(i+15)*ny+kt+21]*tmp[(i+15)];
              y_copy[(kt-kk+22)]=y_copy[(kt-kk+22)]+A[(i+15)*ny+kt+22]*tmp[(i+15)];
              y_copy[(kt-kk+23)]=y_copy[(kt-kk+23)]+A[(i+15)*ny+kt+23]*tmp[(i+15)];
              y_copy[(kt-kk+24)]=y_copy[(kt-kk+24)]+A[(i+15)*ny+kt+24]*tmp[(i+15)];
              y_copy[(kt-kk+25)]=y_copy[(kt-kk+25)]+A[(i+15)*ny+kt+25]*tmp[(i+15)];
              y_copy[(kt-kk+26)]=y_copy[(kt-kk+26)]+A[(i+15)*ny+kt+26]*tmp[(i+15)];
              y_copy[(kt-kk+27)]=y_copy[(kt-kk+27)]+A[(i+15)*ny+kt+27]*tmp[(i+15)];
              y_copy[(kt-kk+28)]=y_copy[(kt-kk+28)]+A[(i+15)*ny+kt+28]*tmp[(i+15)];
              y_copy[(kt-kk+29)]=y_copy[(kt-kk+29)]+A[(i+15)*ny+kt+29]*tmp[(i+15)];
              y_copy[(kt-kk+30)]=y_copy[(kt-kk+30)]+A[(i+15)*ny+kt+30]*tmp[(i+15)];
              y_copy[(kt-kk+31)]=y_copy[(kt-kk+31)]+A[(i+15)*ny+kt+31]*tmp[(i+15)];
              y_copy[(kt-kk)]=y_copy[(kt-kk)]+A[(i+16)*ny+kt]*tmp[(i+16)];
              y_copy[(kt-kk+1)]=y_copy[(kt-kk+1)]+A[(i+16)*ny+kt+1]*tmp[(i+16)];
              y_copy[(kt-kk+2)]=y_copy[(kt-kk+2)]+A[(i+16)*ny+kt+2]*tmp[(i+16)];
              y_copy[(kt-kk+3)]=y_copy[(kt-kk+3)]+A[(i+16)*ny+kt+3]*tmp[(i+16)];
              y_copy[(kt-kk+4)]=y_copy[(kt-kk+4)]+A[(i+16)*ny+kt+4]*tmp[(i+16)];
              y_copy[(kt-kk+5)]=y_copy[(kt-kk+5)]+A[(i+16)*ny+kt+5]*tmp[(i+16)];
              y_copy[(kt-kk+6)]=y_copy[(kt-kk+6)]+A[(i+16)*ny+kt+6]*tmp[(i+16)];
              y_copy[(kt-kk+7)]=y_copy[(kt-kk+7)]+A[(i+16)*ny+kt+7]*tmp[(i+16)];
              y_copy[(kt-kk+8)]=y_copy[(kt-kk+8)]+A[(i+16)*ny+kt+8]*tmp[(i+16)];
              y_copy[(kt-kk+9)]=y_copy[(kt-kk+9)]+A[(i+16)*ny+kt+9]*tmp[(i+16)];
              y_copy[(kt-kk+10)]=y_copy[(kt-kk+10)]+A[(i+16)*ny+kt+10]*tmp[(i+16)];
              y_copy[(kt-kk+11)]=y_copy[(kt-kk+11)]+A[(i+16)*ny+kt+11]*tmp[(i+16)];
              y_copy[(kt-kk+12)]=y_copy[(kt-kk+12)]+A[(i+16)*ny+kt+12]*tmp[(i+16)];
              y_copy[(kt-kk+13)]=y_copy[(kt-kk+13)]+A[(i+16)*ny+kt+13]*tmp[(i+16)];
              y_copy[(kt-kk+14)]=y_copy[(kt-kk+14)]+A[(i+16)*ny+kt+14]*tmp[(i+16)];
              y_copy[(kt-kk+15)]=y_copy[(kt-kk+15)]+A[(i+16)*ny+kt+15]*tmp[(i+16)];
              y_copy[(kt-kk+16)]=y_copy[(kt-kk+16)]+A[(i+16)*ny+kt+16]*tmp[(i+16)];
              y_copy[(kt-kk+17)]=y_copy[(kt-kk+17)]+A[(i+16)*ny+kt+17]*tmp[(i+16)];
              y_copy[(kt-kk+18)]=y_copy[(kt-kk+18)]+A[(i+16)*ny+kt+18]*tmp[(i+16)];
              y_copy[(kt-kk+19)]=y_copy[(kt-kk+19)]+A[(i+16)*ny+kt+19]*tmp[(i+16)];
              y_copy[(kt-kk+20)]=y_copy[(kt-kk+20)]+A[(i+16)*ny+kt+20]*tmp[(i+16)];
              y_copy[(kt-kk+21)]=y_copy[(kt-kk+21)]+A[(i+16)*ny+kt+21]*tmp[(i+16)];
              y_copy[(kt-kk+22)]=y_copy[(kt-kk+22)]+A[(i+16)*ny+kt+22]*tmp[(i+16)];
              y_copy[(kt-kk+23)]=y_copy[(kt-kk+23)]+A[(i+16)*ny+kt+23]*tmp[(i+16)];
              y_copy[(kt-kk+24)]=y_copy[(kt-kk+24)]+A[(i+16)*ny+kt+24]*tmp[(i+16)];
              y_copy[(kt-kk+25)]=y_copy[(kt-kk+25)]+A[(i+16)*ny+kt+25]*tmp[(i+16)];
              y_copy[(kt-kk+26)]=y_copy[(kt-kk+26)]+A[(i+16)*ny+kt+26]*tmp[(i+16)];
              y_copy[(kt-kk+27)]=y_copy[(kt-kk+27)]+A[(i+16)*ny+kt+27]*tmp[(i+16)];
              y_copy[(kt-kk+28)]=y_copy[(kt-kk+28)]+A[(i+16)*ny+kt+28]*tmp[(i+16)];
              y_copy[(kt-kk+29)]=y_copy[(kt-kk+29)]+A[(i+16)*ny+kt+29]*tmp[(i+16)];
              y_copy[(kt-kk+30)]=y_copy[(kt-kk+30)]+A[(i+16)*ny+kt+30]*tmp[(i+16)];
              y_copy[(kt-kk+31)]=y_copy[(kt-kk+31)]+A[(i+16)*ny+kt+31]*tmp[(i+16)];
              y_copy[(kt-kk)]=y_copy[(kt-kk)]+A[(i+17)*ny+kt]*tmp[(i+17)];
              y_copy[(kt-kk+1)]=y_copy[(kt-kk+1)]+A[(i+17)*ny+kt+1]*tmp[(i+17)];
              y_copy[(kt-kk+2)]=y_copy[(kt-kk+2)]+A[(i+17)*ny+kt+2]*tmp[(i+17)];
              y_copy[(kt-kk+3)]=y_copy[(kt-kk+3)]+A[(i+17)*ny+kt+3]*tmp[(i+17)];
              y_copy[(kt-kk+4)]=y_copy[(kt-kk+4)]+A[(i+17)*ny+kt+4]*tmp[(i+17)];
              y_copy[(kt-kk+5)]=y_copy[(kt-kk+5)]+A[(i+17)*ny+kt+5]*tmp[(i+17)];
              y_copy[(kt-kk+6)]=y_copy[(kt-kk+6)]+A[(i+17)*ny+kt+6]*tmp[(i+17)];
              y_copy[(kt-kk+7)]=y_copy[(kt-kk+7)]+A[(i+17)*ny+kt+7]*tmp[(i+17)];
              y_copy[(kt-kk+8)]=y_copy[(kt-kk+8)]+A[(i+17)*ny+kt+8]*tmp[(i+17)];
              y_copy[(kt-kk+9)]=y_copy[(kt-kk+9)]+A[(i+17)*ny+kt+9]*tmp[(i+17)];
              y_copy[(kt-kk+10)]=y_copy[(kt-kk+10)]+A[(i+17)*ny+kt+10]*tmp[(i+17)];
              y_copy[(kt-kk+11)]=y_copy[(kt-kk+11)]+A[(i+17)*ny+kt+11]*tmp[(i+17)];
              y_copy[(kt-kk+12)]=y_copy[(kt-kk+12)]+A[(i+17)*ny+kt+12]*tmp[(i+17)];
              y_copy[(kt-kk+13)]=y_copy[(kt-kk+13)]+A[(i+17)*ny+kt+13]*tmp[(i+17)];
              y_copy[(kt-kk+14)]=y_copy[(kt-kk+14)]+A[(i+17)*ny+kt+14]*tmp[(i+17)];
              y_copy[(kt-kk+15)]=y_copy[(kt-kk+15)]+A[(i+17)*ny+kt+15]*tmp[(i+17)];
              y_copy[(kt-kk+16)]=y_copy[(kt-kk+16)]+A[(i+17)*ny+kt+16]*tmp[(i+17)];
              y_copy[(kt-kk+17)]=y_copy[(kt-kk+17)]+A[(i+17)*ny+kt+17]*tmp[(i+17)];
              y_copy[(kt-kk+18)]=y_copy[(kt-kk+18)]+A[(i+17)*ny+kt+18]*tmp[(i+17)];
              y_copy[(kt-kk+19)]=y_copy[(kt-kk+19)]+A[(i+17)*ny+kt+19]*tmp[(i+17)];
              y_copy[(kt-kk+20)]=y_copy[(kt-kk+20)]+A[(i+17)*ny+kt+20]*tmp[(i+17)];
              y_copy[(kt-kk+21)]=y_copy[(kt-kk+21)]+A[(i+17)*ny+kt+21]*tmp[(i+17)];
              y_copy[(kt-kk+22)]=y_copy[(kt-kk+22)]+A[(i+17)*ny+kt+22]*tmp[(i+17)];
              y_copy[(kt-kk+23)]=y_copy[(kt-kk+23)]+A[(i+17)*ny+kt+23]*tmp[(i+17)];
              y_copy[(kt-kk+24)]=y_copy[(kt-kk+24)]+A[(i+17)*ny+kt+24]*tmp[(i+17)];
              y_copy[(kt-kk+25)]=y_copy[(kt-kk+25)]+A[(i+17)*ny+kt+25]*tmp[(i+17)];
              y_copy[(kt-kk+26)]=y_copy[(kt-kk+26)]+A[(i+17)*ny+kt+26]*tmp[(i+17)];
              y_copy[(kt-kk+27)]=y_copy[(kt-kk+27)]+A[(i+17)*ny+kt+27]*tmp[(i+17)];
              y_copy[(kt-kk+28)]=y_copy[(kt-kk+28)]+A[(i+17)*ny+kt+28]*tmp[(i+17)];
              y_copy[(kt-kk+29)]=y_copy[(kt-kk+29)]+A[(i+17)*ny+kt+29]*tmp[(i+17)];
              y_copy[(kt-kk+30)]=y_copy[(kt-kk+30)]+A[(i+17)*ny+kt+30]*tmp[(i+17)];
              y_copy[(kt-kk+31)]=y_copy[(kt-kk+31)]+A[(i+17)*ny+kt+31]*tmp[(i+17)];
              y_copy[(kt-kk)]=y_copy[(kt-kk)]+A[(i+18)*ny+kt]*tmp[(i+18)];
              y_copy[(kt-kk+1)]=y_copy[(kt-kk+1)]+A[(i+18)*ny+kt+1]*tmp[(i+18)];
              y_copy[(kt-kk+2)]=y_copy[(kt-kk+2)]+A[(i+18)*ny+kt+2]*tmp[(i+18)];
              y_copy[(kt-kk+3)]=y_copy[(kt-kk+3)]+A[(i+18)*ny+kt+3]*tmp[(i+18)];
              y_copy[(kt-kk+4)]=y_copy[(kt-kk+4)]+A[(i+18)*ny+kt+4]*tmp[(i+18)];
              y_copy[(kt-kk+5)]=y_copy[(kt-kk+5)]+A[(i+18)*ny+kt+5]*tmp[(i+18)];
              y_copy[(kt-kk+6)]=y_copy[(kt-kk+6)]+A[(i+18)*ny+kt+6]*tmp[(i+18)];
              y_copy[(kt-kk+7)]=y_copy[(kt-kk+7)]+A[(i+18)*ny+kt+7]*tmp[(i+18)];
              y_copy[(kt-kk+8)]=y_copy[(kt-kk+8)]+A[(i+18)*ny+kt+8]*tmp[(i+18)];
              y_copy[(kt-kk+9)]=y_copy[(kt-kk+9)]+A[(i+18)*ny+kt+9]*tmp[(i+18)];
              y_copy[(kt-kk+10)]=y_copy[(kt-kk+10)]+A[(i+18)*ny+kt+10]*tmp[(i+18)];
              y_copy[(kt-kk+11)]=y_copy[(kt-kk+11)]+A[(i+18)*ny+kt+11]*tmp[(i+18)];
              y_copy[(kt-kk+12)]=y_copy[(kt-kk+12)]+A[(i+18)*ny+kt+12]*tmp[(i+18)];
              y_copy[(kt-kk+13)]=y_copy[(kt-kk+13)]+A[(i+18)*ny+kt+13]*tmp[(i+18)];
              y_copy[(kt-kk+14)]=y_copy[(kt-kk+14)]+A[(i+18)*ny+kt+14]*tmp[(i+18)];
              y_copy[(kt-kk+15)]=y_copy[(kt-kk+15)]+A[(i+18)*ny+kt+15]*tmp[(i+18)];
              y_copy[(kt-kk+16)]=y_copy[(kt-kk+16)]+A[(i+18)*ny+kt+16]*tmp[(i+18)];
              y_copy[(kt-kk+17)]=y_copy[(kt-kk+17)]+A[(i+18)*ny+kt+17]*tmp[(i+18)];
              y_copy[(kt-kk+18)]=y_copy[(kt-kk+18)]+A[(i+18)*ny+kt+18]*tmp[(i+18)];
              y_copy[(kt-kk+19)]=y_copy[(kt-kk+19)]+A[(i+18)*ny+kt+19]*tmp[(i+18)];
              y_copy[(kt-kk+20)]=y_copy[(kt-kk+20)]+A[(i+18)*ny+kt+20]*tmp[(i+18)];
              y_copy[(kt-kk+21)]=y_copy[(kt-kk+21)]+A[(i+18)*ny+kt+21]*tmp[(i+18)];
              y_copy[(kt-kk+22)]=y_copy[(kt-kk+22)]+A[(i+18)*ny+kt+22]*tmp[(i+18)];
              y_copy[(kt-kk+23)]=y_copy[(kt-kk+23)]+A[(i+18)*ny+kt+23]*tmp[(i+18)];
              y_copy[(kt-kk+24)]=y_copy[(kt-kk+24)]+A[(i+18)*ny+kt+24]*tmp[(i+18)];
              y_copy[(kt-kk+25)]=y_copy[(kt-kk+25)]+A[(i+18)*ny+kt+25]*tmp[(i+18)];
              y_copy[(kt-kk+26)]=y_copy[(kt-kk+26)]+A[(i+18)*ny+kt+26]*tmp[(i+18)];
              y_copy[(kt-kk+27)]=y_copy[(kt-kk+27)]+A[(i+18)*ny+kt+27]*tmp[(i+18)];
              y_copy[(kt-kk+28)]=y_copy[(kt-kk+28)]+A[(i+18)*ny+kt+28]*tmp[(i+18)];
              y_copy[(kt-kk+29)]=y_copy[(kt-kk+29)]+A[(i+18)*ny+kt+29]*tmp[(i+18)];
              y_copy[(kt-kk+30)]=y_copy[(kt-kk+30)]+A[(i+18)*ny+kt+30]*tmp[(i+18)];
              y_copy[(kt-kk+31)]=y_copy[(kt-kk+31)]+A[(i+18)*ny+kt+31]*tmp[(i+18)];
              y_copy[(kt-kk)]=y_copy[(kt-kk)]+A[(i+19)*ny+kt]*tmp[(i+19)];
              y_copy[(kt-kk+1)]=y_copy[(kt-kk+1)]+A[(i+19)*ny+kt+1]*tmp[(i+19)];
              y_copy[(kt-kk+2)]=y_copy[(kt-kk+2)]+A[(i+19)*ny+kt+2]*tmp[(i+19)];
              y_copy[(kt-kk+3)]=y_copy[(kt-kk+3)]+A[(i+19)*ny+kt+3]*tmp[(i+19)];
              y_copy[(kt-kk+4)]=y_copy[(kt-kk+4)]+A[(i+19)*ny+kt+4]*tmp[(i+19)];
              y_copy[(kt-kk+5)]=y_copy[(kt-kk+5)]+A[(i+19)*ny+kt+5]*tmp[(i+19)];
              y_copy[(kt-kk+6)]=y_copy[(kt-kk+6)]+A[(i+19)*ny+kt+6]*tmp[(i+19)];
              y_copy[(kt-kk+7)]=y_copy[(kt-kk+7)]+A[(i+19)*ny+kt+7]*tmp[(i+19)];
              y_copy[(kt-kk+8)]=y_copy[(kt-kk+8)]+A[(i+19)*ny+kt+8]*tmp[(i+19)];
              y_copy[(kt-kk+9)]=y_copy[(kt-kk+9)]+A[(i+19)*ny+kt+9]*tmp[(i+19)];
              y_copy[(kt-kk+10)]=y_copy[(kt-kk+10)]+A[(i+19)*ny+kt+10]*tmp[(i+19)];
              y_copy[(kt-kk+11)]=y_copy[(kt-kk+11)]+A[(i+19)*ny+kt+11]*tmp[(i+19)];
              y_copy[(kt-kk+12)]=y_copy[(kt-kk+12)]+A[(i+19)*ny+kt+12]*tmp[(i+19)];
              y_copy[(kt-kk+13)]=y_copy[(kt-kk+13)]+A[(i+19)*ny+kt+13]*tmp[(i+19)];
              y_copy[(kt-kk+14)]=y_copy[(kt-kk+14)]+A[(i+19)*ny+kt+14]*tmp[(i+19)];
              y_copy[(kt-kk+15)]=y_copy[(kt-kk+15)]+A[(i+19)*ny+kt+15]*tmp[(i+19)];
              y_copy[(kt-kk+16)]=y_copy[(kt-kk+16)]+A[(i+19)*ny+kt+16]*tmp[(i+19)];
              y_copy[(kt-kk+17)]=y_copy[(kt-kk+17)]+A[(i+19)*ny+kt+17]*tmp[(i+19)];
              y_copy[(kt-kk+18)]=y_copy[(kt-kk+18)]+A[(i+19)*ny+kt+18]*tmp[(i+19)];
              y_copy[(kt-kk+19)]=y_copy[(kt-kk+19)]+A[(i+19)*ny+kt+19]*tmp[(i+19)];
              y_copy[(kt-kk+20)]=y_copy[(kt-kk+20)]+A[(i+19)*ny+kt+20]*tmp[(i+19)];
              y_copy[(kt-kk+21)]=y_copy[(kt-kk+21)]+A[(i+19)*ny+kt+21]*tmp[(i+19)];
              y_copy[(kt-kk+22)]=y_copy[(kt-kk+22)]+A[(i+19)*ny+kt+22]*tmp[(i+19)];
              y_copy[(kt-kk+23)]=y_copy[(kt-kk+23)]+A[(i+19)*ny+kt+23]*tmp[(i+19)];
              y_copy[(kt-kk+24)]=y_copy[(kt-kk+24)]+A[(i+19)*ny+kt+24]*tmp[(i+19)];
              y_copy[(kt-kk+25)]=y_copy[(kt-kk+25)]+A[(i+19)*ny+kt+25]*tmp[(i+19)];
              y_copy[(kt-kk+26)]=y_copy[(kt-kk+26)]+A[(i+19)*ny+kt+26]*tmp[(i+19)];
              y_copy[(kt-kk+27)]=y_copy[(kt-kk+27)]+A[(i+19)*ny+kt+27]*tmp[(i+19)];
              y_copy[(kt-kk+28)]=y_copy[(kt-kk+28)]+A[(i+19)*ny+kt+28]*tmp[(i+19)];
              y_copy[(kt-kk+29)]=y_copy[(kt-kk+29)]+A[(i+19)*ny+kt+29]*tmp[(i+19)];
              y_copy[(kt-kk+30)]=y_copy[(kt-kk+30)]+A[(i+19)*ny+kt+30]*tmp[(i+19)];
              y_copy[(kt-kk+31)]=y_copy[(kt-kk+31)]+A[(i+19)*ny+kt+31]*tmp[(i+19)];
              y_copy[(kt-kk)]=y_copy[(kt-kk)]+A[(i+20)*ny+kt]*tmp[(i+20)];
              y_copy[(kt-kk+1)]=y_copy[(kt-kk+1)]+A[(i+20)*ny+kt+1]*tmp[(i+20)];
              y_copy[(kt-kk+2)]=y_copy[(kt-kk+2)]+A[(i+20)*ny+kt+2]*tmp[(i+20)];
              y_copy[(kt-kk+3)]=y_copy[(kt-kk+3)]+A[(i+20)*ny+kt+3]*tmp[(i+20)];
              y_copy[(kt-kk+4)]=y_copy[(kt-kk+4)]+A[(i+20)*ny+kt+4]*tmp[(i+20)];
              y_copy[(kt-kk+5)]=y_copy[(kt-kk+5)]+A[(i+20)*ny+kt+5]*tmp[(i+20)];
              y_copy[(kt-kk+6)]=y_copy[(kt-kk+6)]+A[(i+20)*ny+kt+6]*tmp[(i+20)];
              y_copy[(kt-kk+7)]=y_copy[(kt-kk+7)]+A[(i+20)*ny+kt+7]*tmp[(i+20)];
              y_copy[(kt-kk+8)]=y_copy[(kt-kk+8)]+A[(i+20)*ny+kt+8]*tmp[(i+20)];
              y_copy[(kt-kk+9)]=y_copy[(kt-kk+9)]+A[(i+20)*ny+kt+9]*tmp[(i+20)];
              y_copy[(kt-kk+10)]=y_copy[(kt-kk+10)]+A[(i+20)*ny+kt+10]*tmp[(i+20)];
              y_copy[(kt-kk+11)]=y_copy[(kt-kk+11)]+A[(i+20)*ny+kt+11]*tmp[(i+20)];
              y_copy[(kt-kk+12)]=y_copy[(kt-kk+12)]+A[(i+20)*ny+kt+12]*tmp[(i+20)];
              y_copy[(kt-kk+13)]=y_copy[(kt-kk+13)]+A[(i+20)*ny+kt+13]*tmp[(i+20)];
              y_copy[(kt-kk+14)]=y_copy[(kt-kk+14)]+A[(i+20)*ny+kt+14]*tmp[(i+20)];
              y_copy[(kt-kk+15)]=y_copy[(kt-kk+15)]+A[(i+20)*ny+kt+15]*tmp[(i+20)];
              y_copy[(kt-kk+16)]=y_copy[(kt-kk+16)]+A[(i+20)*ny+kt+16]*tmp[(i+20)];
              y_copy[(kt-kk+17)]=y_copy[(kt-kk+17)]+A[(i+20)*ny+kt+17]*tmp[(i+20)];
              y_copy[(kt-kk+18)]=y_copy[(kt-kk+18)]+A[(i+20)*ny+kt+18]*tmp[(i+20)];
              y_copy[(kt-kk+19)]=y_copy[(kt-kk+19)]+A[(i+20)*ny+kt+19]*tmp[(i+20)];
              y_copy[(kt-kk+20)]=y_copy[(kt-kk+20)]+A[(i+20)*ny+kt+20]*tmp[(i+20)];
              y_copy[(kt-kk+21)]=y_copy[(kt-kk+21)]+A[(i+20)*ny+kt+21]*tmp[(i+20)];
              y_copy[(kt-kk+22)]=y_copy[(kt-kk+22)]+A[(i+20)*ny+kt+22]*tmp[(i+20)];
              y_copy[(kt-kk+23)]=y_copy[(kt-kk+23)]+A[(i+20)*ny+kt+23]*tmp[(i+20)];
              y_copy[(kt-kk+24)]=y_copy[(kt-kk+24)]+A[(i+20)*ny+kt+24]*tmp[(i+20)];
              y_copy[(kt-kk+25)]=y_copy[(kt-kk+25)]+A[(i+20)*ny+kt+25]*tmp[(i+20)];
              y_copy[(kt-kk+26)]=y_copy[(kt-kk+26)]+A[(i+20)*ny+kt+26]*tmp[(i+20)];
              y_copy[(kt-kk+27)]=y_copy[(kt-kk+27)]+A[(i+20)*ny+kt+27]*tmp[(i+20)];
              y_copy[(kt-kk+28)]=y_copy[(kt-kk+28)]+A[(i+20)*ny+kt+28]*tmp[(i+20)];
              y_copy[(kt-kk+29)]=y_copy[(kt-kk+29)]+A[(i+20)*ny+kt+29]*tmp[(i+20)];
              y_copy[(kt-kk+30)]=y_copy[(kt-kk+30)]+A[(i+20)*ny+kt+30]*tmp[(i+20)];
              y_copy[(kt-kk+31)]=y_copy[(kt-kk+31)]+A[(i+20)*ny+kt+31]*tmp[(i+20)];
              y_copy[(kt-kk)]=y_copy[(kt-kk)]+A[(i+21)*ny+kt]*tmp[(i+21)];
              y_copy[(kt-kk+1)]=y_copy[(kt-kk+1)]+A[(i+21)*ny+kt+1]*tmp[(i+21)];
              y_copy[(kt-kk+2)]=y_copy[(kt-kk+2)]+A[(i+21)*ny+kt+2]*tmp[(i+21)];
              y_copy[(kt-kk+3)]=y_copy[(kt-kk+3)]+A[(i+21)*ny+kt+3]*tmp[(i+21)];
              y_copy[(kt-kk+4)]=y_copy[(kt-kk+4)]+A[(i+21)*ny+kt+4]*tmp[(i+21)];
              y_copy[(kt-kk+5)]=y_copy[(kt-kk+5)]+A[(i+21)*ny+kt+5]*tmp[(i+21)];
              y_copy[(kt-kk+6)]=y_copy[(kt-kk+6)]+A[(i+21)*ny+kt+6]*tmp[(i+21)];
              y_copy[(kt-kk+7)]=y_copy[(kt-kk+7)]+A[(i+21)*ny+kt+7]*tmp[(i+21)];
              y_copy[(kt-kk+8)]=y_copy[(kt-kk+8)]+A[(i+21)*ny+kt+8]*tmp[(i+21)];
              y_copy[(kt-kk+9)]=y_copy[(kt-kk+9)]+A[(i+21)*ny+kt+9]*tmp[(i+21)];
              y_copy[(kt-kk+10)]=y_copy[(kt-kk+10)]+A[(i+21)*ny+kt+10]*tmp[(i+21)];
              y_copy[(kt-kk+11)]=y_copy[(kt-kk+11)]+A[(i+21)*ny+kt+11]*tmp[(i+21)];
              y_copy[(kt-kk+12)]=y_copy[(kt-kk+12)]+A[(i+21)*ny+kt+12]*tmp[(i+21)];
              y_copy[(kt-kk+13)]=y_copy[(kt-kk+13)]+A[(i+21)*ny+kt+13]*tmp[(i+21)];
              y_copy[(kt-kk+14)]=y_copy[(kt-kk+14)]+A[(i+21)*ny+kt+14]*tmp[(i+21)];
              y_copy[(kt-kk+15)]=y_copy[(kt-kk+15)]+A[(i+21)*ny+kt+15]*tmp[(i+21)];
              y_copy[(kt-kk+16)]=y_copy[(kt-kk+16)]+A[(i+21)*ny+kt+16]*tmp[(i+21)];
              y_copy[(kt-kk+17)]=y_copy[(kt-kk+17)]+A[(i+21)*ny+kt+17]*tmp[(i+21)];
              y_copy[(kt-kk+18)]=y_copy[(kt-kk+18)]+A[(i+21)*ny+kt+18]*tmp[(i+21)];
              y_copy[(kt-kk+19)]=y_copy[(kt-kk+19)]+A[(i+21)*ny+kt+19]*tmp[(i+21)];
              y_copy[(kt-kk+20)]=y_copy[(kt-kk+20)]+A[(i+21)*ny+kt+20]*tmp[(i+21)];
              y_copy[(kt-kk+21)]=y_copy[(kt-kk+21)]+A[(i+21)*ny+kt+21]*tmp[(i+21)];
              y_copy[(kt-kk+22)]=y_copy[(kt-kk+22)]+A[(i+21)*ny+kt+22]*tmp[(i+21)];
              y_copy[(kt-kk+23)]=y_copy[(kt-kk+23)]+A[(i+21)*ny+kt+23]*tmp[(i+21)];
              y_copy[(kt-kk+24)]=y_copy[(kt-kk+24)]+A[(i+21)*ny+kt+24]*tmp[(i+21)];
              y_copy[(kt-kk+25)]=y_copy[(kt-kk+25)]+A[(i+21)*ny+kt+25]*tmp[(i+21)];
              y_copy[(kt-kk+26)]=y_copy[(kt-kk+26)]+A[(i+21)*ny+kt+26]*tmp[(i+21)];
              y_copy[(kt-kk+27)]=y_copy[(kt-kk+27)]+A[(i+21)*ny+kt+27]*tmp[(i+21)];
              y_copy[(kt-kk+28)]=y_copy[(kt-kk+28)]+A[(i+21)*ny+kt+28]*tmp[(i+21)];
              y_copy[(kt-kk+29)]=y_copy[(kt-kk+29)]+A[(i+21)*ny+kt+29]*tmp[(i+21)];
              y_copy[(kt-kk+30)]=y_copy[(kt-kk+30)]+A[(i+21)*ny+kt+30]*tmp[(i+21)];
              y_copy[(kt-kk+31)]=y_copy[(kt-kk+31)]+A[(i+21)*ny+kt+31]*tmp[(i+21)];
            }
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y[k]=y_copy[(k-kk)];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y[k]=y_copy[(k-kk)];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y[k]=y_copy[(k-kk)];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y[k]=y_copy[(k-kk)];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y[k]=y_copy[(k-kk)];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y[k]=y_copy[(k-kk)];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y[k]=y_copy[(k-kk)];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y[k]=y_copy[(k-kk)];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y[k]=y_copy[(k-kk)];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y[k]=y_copy[(k-kk)];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y[k]=y_copy[(k-kk)];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y[k]=y_copy[(k-kk)];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y[k]=y_copy[(k-kk)];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y[k]=y_copy[(k-kk)];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y[k]=y_copy[(k-kk)];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y[k]=y_copy[(k-kk)];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y[k]=y_copy[(k-kk)];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y[k]=y_copy[(k-kk)];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y[k]=y_copy[(k-kk)];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y[k]=y_copy[(k-kk)];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y[k]=y_copy[(k-kk)];
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y[k]=y_copy[(k-kk)];
          }
        }
      }
      for (i=(min(nx-1,ii+31))-(((min(nx-1,ii+31))-(ii)+1)%22)+1; i<=min(nx-1,ii+31); i=i+1) {
        tmp[i]=0;
        for (jjj=0; jjj<=ny-1; jjj=jjj+2048) 
          for (jj=jjj; jj<=min(ny-1,jjj+1536); jj=jj+512) {
            register int cbv_3;
            cbv_3=min(ny-1,jj+511);
#pragma ivdep
#pragma vector always
            for (j=jj; j<=cbv_3; j=j+1) 
              tmp[i]=tmp[i]+A[i*ny+j]*x[j];
          }
        for (kkk=0; kkk<=ny-1; kkk=kkk+2048) 
          for (kk=kkk; kk<=min(ny-1,kkk+1984); kk=kk+64) {
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y_copy[(k-kk)]=y[k];
            register int cbv_4;
            cbv_4=min(ny-1,kk+63)-31;
#pragma ivdep
#pragma vector always
            for (kt=kk; kt<=cbv_4; kt=kt+32) {
              y_copy[(kt-kk)]=y_copy[(kt-kk)]+A[i*ny+kt]*tmp[i];
              y_copy[((kt+1)-kk)]=y_copy[((kt+1)-kk)]+A[i*ny+(kt+1)]*tmp[i];
              y_copy[((kt+2)-kk)]=y_copy[((kt+2)-kk)]+A[i*ny+(kt+2)]*tmp[i];
              y_copy[((kt+3)-kk)]=y_copy[((kt+3)-kk)]+A[i*ny+(kt+3)]*tmp[i];
              y_copy[((kt+4)-kk)]=y_copy[((kt+4)-kk)]+A[i*ny+(kt+4)]*tmp[i];
              y_copy[((kt+5)-kk)]=y_copy[((kt+5)-kk)]+A[i*ny+(kt+5)]*tmp[i];
              y_copy[((kt+6)-kk)]=y_copy[((kt+6)-kk)]+A[i*ny+(kt+6)]*tmp[i];
              y_copy[((kt+7)-kk)]=y_copy[((kt+7)-kk)]+A[i*ny+(kt+7)]*tmp[i];
              y_copy[((kt+8)-kk)]=y_copy[((kt+8)-kk)]+A[i*ny+(kt+8)]*tmp[i];
              y_copy[((kt+9)-kk)]=y_copy[((kt+9)-kk)]+A[i*ny+(kt+9)]*tmp[i];
              y_copy[((kt+10)-kk)]=y_copy[((kt+10)-kk)]+A[i*ny+(kt+10)]*tmp[i];
              y_copy[((kt+11)-kk)]=y_copy[((kt+11)-kk)]+A[i*ny+(kt+11)]*tmp[i];
              y_copy[((kt+12)-kk)]=y_copy[((kt+12)-kk)]+A[i*ny+(kt+12)]*tmp[i];
              y_copy[((kt+13)-kk)]=y_copy[((kt+13)-kk)]+A[i*ny+(kt+13)]*tmp[i];
              y_copy[((kt+14)-kk)]=y_copy[((kt+14)-kk)]+A[i*ny+(kt+14)]*tmp[i];
              y_copy[((kt+15)-kk)]=y_copy[((kt+15)-kk)]+A[i*ny+(kt+15)]*tmp[i];
              y_copy[((kt+16)-kk)]=y_copy[((kt+16)-kk)]+A[i*ny+(kt+16)]*tmp[i];
              y_copy[((kt+17)-kk)]=y_copy[((kt+17)-kk)]+A[i*ny+(kt+17)]*tmp[i];
              y_copy[((kt+18)-kk)]=y_copy[((kt+18)-kk)]+A[i*ny+(kt+18)]*tmp[i];
              y_copy[((kt+19)-kk)]=y_copy[((kt+19)-kk)]+A[i*ny+(kt+19)]*tmp[i];
              y_copy[((kt+20)-kk)]=y_copy[((kt+20)-kk)]+A[i*ny+(kt+20)]*tmp[i];
              y_copy[((kt+21)-kk)]=y_copy[((kt+21)-kk)]+A[i*ny+(kt+21)]*tmp[i];
              y_copy[((kt+22)-kk)]=y_copy[((kt+22)-kk)]+A[i*ny+(kt+22)]*tmp[i];
              y_copy[((kt+23)-kk)]=y_copy[((kt+23)-kk)]+A[i*ny+(kt+23)]*tmp[i];
              y_copy[((kt+24)-kk)]=y_copy[((kt+24)-kk)]+A[i*ny+(kt+24)]*tmp[i];
              y_copy[((kt+25)-kk)]=y_copy[((kt+25)-kk)]+A[i*ny+(kt+25)]*tmp[i];
              y_copy[((kt+26)-kk)]=y_copy[((kt+26)-kk)]+A[i*ny+(kt+26)]*tmp[i];
              y_copy[((kt+27)-kk)]=y_copy[((kt+27)-kk)]+A[i*ny+(kt+27)]*tmp[i];
              y_copy[((kt+28)-kk)]=y_copy[((kt+28)-kk)]+A[i*ny+(kt+28)]*tmp[i];
              y_copy[((kt+29)-kk)]=y_copy[((kt+29)-kk)]+A[i*ny+(kt+29)]*tmp[i];
              y_copy[((kt+30)-kk)]=y_copy[((kt+30)-kk)]+A[i*ny+(kt+30)]*tmp[i];
              y_copy[((kt+31)-kk)]=y_copy[((kt+31)-kk)]+A[i*ny+(kt+31)]*tmp[i];
            }
            for (k=kk; k<=min(ny-1,kk+63); k=k+1) 
              y[k]=y_copy[(k-kk)];
          }
      }
    }
}
/*@ end @*/



    
    orio_t_end = getClock();
    orio_t = orio_t_end - orio_t_start;
    printf("{'[2, 6, 3, 6, 6, 6, 0, 1, 13, 21, 0, 2, 0, 0, 2, 0, 0, 1]' : %g}\n", orio_t);
    
    if (orio_i==0) {
      
      if (!isValid()) {
         fprintf(stdout,"validation function isValid returned an error code.\n");
         return 1;
      }
    }
  }
  
  
  
  return 0;
}
