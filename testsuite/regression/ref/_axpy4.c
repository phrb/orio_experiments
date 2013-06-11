
void axpy4(int n, double *y, double a1, double *x1, double a2, double *x2, double a3, double *x3,
           double a4, double *x4) {

    register int i;

/*@ begin Align (x1[],x2[],x3[],x4[],y[]) @*/
#pragma disjoint (*x1,*x2,*x3,*x4,*y) 
if ((((int)(x1)|(int)(x2)|(int)(x3)|(int)(x4)|(int)(y)) & 0xF) == 0) {
  __alignx(16,x1); 
  __alignx(16,x2); 
  __alignx(16,x3); 
  __alignx(16,x4); 
  __alignx(16,y); 

  /*@ begin Loop (
    transform Unroll(ufactor=20, parallelize=True)
      for (i=0; i<=n-1; i++)
        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
        ) @*/
  {
    int i;
  #pragma omp parallel for private(i)
    for (i=0; i<=n-20; i=i+20) {
      y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
      y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
      y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
      y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
      y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
      y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
      y[(i+6)]=y[(i+6)]+a1*x1[(i+6)]+a2*x2[(i+6)]+a3*x3[(i+6)]+a4*x4[(i+6)];
      y[(i+7)]=y[(i+7)]+a1*x1[(i+7)]+a2*x2[(i+7)]+a3*x3[(i+7)]+a4*x4[(i+7)];
      y[(i+8)]=y[(i+8)]+a1*x1[(i+8)]+a2*x2[(i+8)]+a3*x3[(i+8)]+a4*x4[(i+8)];
      y[(i+9)]=y[(i+9)]+a1*x1[(i+9)]+a2*x2[(i+9)]+a3*x3[(i+9)]+a4*x4[(i+9)];
      y[(i+10)]=y[(i+10)]+a1*x1[(i+10)]+a2*x2[(i+10)]+a3*x3[(i+10)]+a4*x4[(i+10)];
      y[(i+11)]=y[(i+11)]+a1*x1[(i+11)]+a2*x2[(i+11)]+a3*x3[(i+11)]+a4*x4[(i+11)];
      y[(i+12)]=y[(i+12)]+a1*x1[(i+12)]+a2*x2[(i+12)]+a3*x3[(i+12)]+a4*x4[(i+12)];
      y[(i+13)]=y[(i+13)]+a1*x1[(i+13)]+a2*x2[(i+13)]+a3*x3[(i+13)]+a4*x4[(i+13)];
      y[(i+14)]=y[(i+14)]+a1*x1[(i+14)]+a2*x2[(i+14)]+a3*x3[(i+14)]+a4*x4[(i+14)];
      y[(i+15)]=y[(i+15)]+a1*x1[(i+15)]+a2*x2[(i+15)]+a3*x3[(i+15)]+a4*x4[(i+15)];
      y[(i+16)]=y[(i+16)]+a1*x1[(i+16)]+a2*x2[(i+16)]+a3*x3[(i+16)]+a4*x4[(i+16)];
      y[(i+17)]=y[(i+17)]+a1*x1[(i+17)]+a2*x2[(i+17)]+a3*x3[(i+17)]+a4*x4[(i+17)];
      y[(i+18)]=y[(i+18)]+a1*x1[(i+18)]+a2*x2[(i+18)]+a3*x3[(i+18)]+a4*x4[(i+18)];
      y[(i+19)]=y[(i+19)]+a1*x1[(i+19)]+a2*x2[(i+19)]+a3*x3[(i+19)]+a4*x4[(i+19)];
    }
    for (i=n-((n-(0))%20); i<=n-1; i=i+1) 
      y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
  }
  /*@ end @*/
  
} else {

  /*@ begin Loop (
    transform Unroll(ufactor=20, parallelize=True)
      for (i=0; i<=n-1; i++)
        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
        ) @*/
  {
    int i;
  #pragma omp parallel for private(i)
    for (i=0; i<=n-20; i=i+20) {
      y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
      y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
      y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
      y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
      y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
      y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
      y[(i+6)]=y[(i+6)]+a1*x1[(i+6)]+a2*x2[(i+6)]+a3*x3[(i+6)]+a4*x4[(i+6)];
      y[(i+7)]=y[(i+7)]+a1*x1[(i+7)]+a2*x2[(i+7)]+a3*x3[(i+7)]+a4*x4[(i+7)];
      y[(i+8)]=y[(i+8)]+a1*x1[(i+8)]+a2*x2[(i+8)]+a3*x3[(i+8)]+a4*x4[(i+8)];
      y[(i+9)]=y[(i+9)]+a1*x1[(i+9)]+a2*x2[(i+9)]+a3*x3[(i+9)]+a4*x4[(i+9)];
      y[(i+10)]=y[(i+10)]+a1*x1[(i+10)]+a2*x2[(i+10)]+a3*x3[(i+10)]+a4*x4[(i+10)];
      y[(i+11)]=y[(i+11)]+a1*x1[(i+11)]+a2*x2[(i+11)]+a3*x3[(i+11)]+a4*x4[(i+11)];
      y[(i+12)]=y[(i+12)]+a1*x1[(i+12)]+a2*x2[(i+12)]+a3*x3[(i+12)]+a4*x4[(i+12)];
      y[(i+13)]=y[(i+13)]+a1*x1[(i+13)]+a2*x2[(i+13)]+a3*x3[(i+13)]+a4*x4[(i+13)];
      y[(i+14)]=y[(i+14)]+a1*x1[(i+14)]+a2*x2[(i+14)]+a3*x3[(i+14)]+a4*x4[(i+14)];
      y[(i+15)]=y[(i+15)]+a1*x1[(i+15)]+a2*x2[(i+15)]+a3*x3[(i+15)]+a4*x4[(i+15)];
      y[(i+16)]=y[(i+16)]+a1*x1[(i+16)]+a2*x2[(i+16)]+a3*x3[(i+16)]+a4*x4[(i+16)];
      y[(i+17)]=y[(i+17)]+a1*x1[(i+17)]+a2*x2[(i+17)]+a3*x3[(i+17)]+a4*x4[(i+17)];
      y[(i+18)]=y[(i+18)]+a1*x1[(i+18)]+a2*x2[(i+18)]+a3*x3[(i+18)]+a4*x4[(i+18)];
      y[(i+19)]=y[(i+19)]+a1*x1[(i+19)]+a2*x2[(i+19)]+a3*x3[(i+19)]+a4*x4[(i+19)];
    }
    for (i=n-((n-(0))%20); i<=n-1; i=i+1) 
      y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
  }
  /*@ end @*/
  
} 
/*@ end @*/

}

