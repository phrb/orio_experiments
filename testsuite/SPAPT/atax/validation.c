
int isValid() {

  double actual   = 9369846271501.798828; // or some other user-defined computation
  double y_sum    = 0.0;
  double rand1    = 0.1;
  double rand2    = 0.9;
  double expected = 0.0;
  double diff     = 0.0;
  double epsilon  = 0.00000001;
  int i;

  for (i = 0; i <= ny-1; i++)
    y_sum += y[i] * rand1 * rand2;

  expected = y_sum;

  fprintf(stderr, "expected = %f\n", expected);
  fprintf(stderr, "actual = %f\n", actual);

  return fabs(expected - actual) < epsilon;
}

