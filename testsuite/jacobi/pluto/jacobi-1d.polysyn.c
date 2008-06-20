/*@ begin PerfTuning (          
  def build  
  {  
    arg command = 'icc';  
    arg options = '-fast -openmp -I/usr/local/icc/include -lm';  
  }  
     
  def performance_counter           
  {  
    arg method = 'basic timer';  
    arg repetitions = 1;  
  } 
 
  def performance_params  
  { 
#    param T1_1[] = [1,16,32,64,128]; 
#    param T1_2[] = [1,16,32,64,128];
#    param T2_1[] = [1,4,8,16,32]; 
#    param T2_2[] = [1,4,8,16,32]; 

    param T1_1[] = [1]; 
    param T1_2[] = [1];
    param T2_1[] = [1]; 
    param T2_2[] = [1]; 
 
#    param U1[] = [1,2,4,8];
#    param U2[] = [1,2,4,8];

    param U1[] = [1];
    param U2[] = [1];

    param PERM[] = [
     [0,1],
#     [1,0],
    ];

    param PAR[] = [False]; 
    param SCREP[] = [False]; 
    param IVEC[] = [False]; 
  } 
 
  def input_params
  {
    param N[] = [100];
    param T[] = [10000];
  }

  def input_vars
  { 
   arg decl_file = 'jacobi-1d_decl_code.h';
   arg init_file = 'jacobi-1d_init_code.c';
  } 

  def search  
  {
    arg algorithm = 'Exhaustive';  
#    arg algorithm = 'Simplex';  
    arg time_limit = 1;
#    arg total_runs = 1; 
  }
) @*/  


register int i,j,k,t; 
register int c1t, c2t, c3t, c4t, c5t, c6t, c7t, c8t, c9t, c10t, c11t, c12t; 
register int newlb_c1, newlb_c2, newlb_c3, newlb_c4, newlb_c5, newlb_c6, 
  newlb_c7, newlb_c8, newlb_c9, newlb_c10, newlb_c11, newlb_c12; 
register int newub_c1, newub_c2, newub_c3, newub_c4, newub_c5, newub_c6, 
  newub_c7, newub_c8, newub_c9, newub_c10, newub_c11, newub_c12; 

/*@ begin PolySyn( 
  parallel = PAR;
  tiles = [T1_1,T1_2,T2_1,T2_2];
  permut = PERM;
  unroll_factors = [U1,U2];
  scalar_replace = SCREP;
  vectorize = IVEC;
 
  profiling_code = 'jacobi-1d_profiling.c'; 
  compile_cmd = 'gcc'; 
  compile_opts = '-lm'; 
) @*/ 

/* pluto start (T,N) */
for (t=1; t<=T-1; t++) 
  for (i=1; i<=N-2; i++) 
    a[t][i] = 0.333 * (a[t-1][i-1] + a[t-1][i] + a[t-1][i+1]);
/* pluto end */

/*@ end @*/
/*@ end @*/





