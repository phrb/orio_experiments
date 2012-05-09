void VecMDot(int n, int nv, double *x, double *y, double *r) {

  register int i,j;
  /*@ begin PerfTuning (
        def performance_params {
          param TC[] = range(16,33,16);
          param UIF[] = range(1,3);
        }
        def build {
          arg build_command = 'nvcc -arch=sm_20';
        }
        def input_params {
          param N[] = [10];
          param NV[] = [16];
        }
        def input_vars {
          decl dynamic double x[N] = random;
          decl dynamic double y[NV*N] = random;
          decl dynamic double r[NV] = 0;
        }
        def performance_counter {
          arg method = 'basic timer';
          arg repetitions = 1;
        }
  ) @*/

  int n=N;
  int nv=NV;

  /*@ begin Loop(transform CUDA(threadCount=TC, unrollInner=UIF)

  for (j=0; j<=nv-1; j++)
    for (i=0; i<=n-1; i++)
      r[j]+=x[i]*y[j*n+i];

  ) @*/

  for (j=0; j<=nv-1; j++)
    for (i=0; i<=n-1; i++)
      r[j]+=x[i]*y[j*n+i];

  /*@ end @*/
  /*@ end @*/
}