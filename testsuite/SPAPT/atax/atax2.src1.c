/*@ begin PerfTuning (
  def build
  {
  arg build_command = 'timeout 20m gcc -O2 -fopenmp -DDYNAMIC';
  arg libs = '-lm';
  }

  def performance_counter
  {
  arg repetitions = 30;
  }

  def search
  {
    arg algorithm = 'DLMT';
    arg total_runs = 75;
    arg dlmt_federov_sampling = 2;
    arg dlmt_extra_experiments = 10;
    arg dlmt_steps = 2;
    # arg dlmt_quadratic = '["T1_I", "T1_J", "T1_K", "T2_I", "T2_J", "T2_K", "U1_I", "U_I", "U_J", "U_K", "RT_I", "RT_J", "RT_K"]';
    # arg dlmt_cubic = '["T1_I", "T1_J", "T1_K", "T2_I", "T2_J", "T2_K", "U1_I", "U_I", "U_J", "U_K"]';
    # arg dlmt_linear = '["T1_I", "T1_J", "T1_K", "ACOPY_x", "ACOPY_y", "SCR", "VEC1", "VEC2"]';
    arg dlmt_linear = '["T1_I", "T1_J", "T1_K", "T2_I", "T2_J", "T2_K", "ACOPY_x", "ACOPY_y", "U1_I", "U_I", "U_J", "U_K", "RT_I", "RT_J", "RT_K", "SCR", "VEC1", "VEC2"]';
    arg dlmt_inverse = '["T1_I", "T1_J", "T1_K", "T2_I", "T2_J", "T2_K", "U1_I", "U_I", "U_J", "U_K", "RT_I", "RT_J", "RT_K"]';
    # arg dlmt_quadratic = '["RT_I", "RT_J", "RT_K", "T1_I", "T1_J", "T1_K", "T2_I", "T2_J", "T2_K"]';
    # arg dlmt_cubic = '["T1_I", "T1_J", "T1_K", "U1_I", "U_I", "U_J", "U_K"]';
    # arg dlmt_inverse = '["U1_I", "U_I", "U_J", "U_K"]';
    # arg dlmt_interactions = '["T1_I:T1_J", "T1_J:T1_K", "T1_I:T1_K", "T2_I:T2_J", "T2_J:T2_K", "T2_I:T2_K"]';
  }

  def performance_params
  {
    # Cache tiling
    param T1_I[] = [1,16,32,64,128,256,512];
    param T1_J[] = [1,16,32,64,128,256,512];
    param T1_K[] = [1,16,32,64,128,256,512];
    param T2_I[] = [1,64,128,256,512,1024,2048];
    param T2_J[] = [1,64,128,256,512,1024,2048];
    param T2_K[] = [1,64,128,256,512,1024,2048];

    # Array copy
    param ACOPY_x[] = [False,True];
    param ACOPY_y[] = [False,True];

    # Unroll-jam
    param U1_I[] = range(1,31);
    param U_I[]  = range(1,31);
    param U_J[]  = range(1,31);
    param U_K[]  = range(1,31);

    # Register tiling
    param RT_I[] = [1,8,32];
    param RT_J[] = [1,8,32];
    param RT_K[] = [1,8,32];

    # Scalar replacement
    param SCR[]  = [False,True];

    # Vectorization
    param VEC1[] = [False,True];
    param VEC2[] = [False,True];

    # Parallelization
    # param OMP[] = [False,True];
    # openmp = (OMP, 'omp parallel for private(iii,jjj,kkk,ii,jj,kk,i,j,k,y_copy,x_copy)')

    # Constraints
    constraint tileI = ((T2_I == 1) or (T2_I % T1_I == 0));
    constraint tileJ = ((T2_J == 1) or (T2_J % T1_J == 0));
    constraint tileK = ((T2_K == 1) or (T2_K % T1_K == 0));

    constraint reg_capacity = (RT_I*RT_J + RT_I*RT_K + RT_J*RT_K <= 150);
    constraint unroll_limit = ((U_I == 1) or (U_J == 1) or (U_K == 1));

  }

  def input_params
  {
  param N[] = [10000];
  }

  def input_vars
  {
  arg decl_file = 'decl.h';
  arg init_file = 'init.c';
  }

  def validation
  {
  arg validation_file = 'validation.c';
  }
) @*/

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
/*@ end @*/

/*@ end @*/
