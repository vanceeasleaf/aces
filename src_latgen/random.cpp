// Park/Miller RNG

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "random.h"

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836

#define IA1 1366
#define IC1 150889
#define IM1 714025
#define IA2 8121
#define IC2 28411
#define IM2 134456
#define IA3 7141
#define IC3 54773
#define IM3 259200

/* ---------------------------------------------------------------------- */

RanPark::RanPark(int seed_init)
{
  if (seed_init <= 0){printf("\nInvalid seed for Park random # generator\n"); exit(1);}
  seed = seed_init;
  save = 0;
}

/* ----------------------------------------------------------------------
   uniform RN 
------------------------------------------------------------------------- */

double RanPark::uniform()
{
  int k = seed/IQ;
  seed = IA*(seed-k*IQ) - IR*k;
  if (seed < 0) seed += IM;
  double ans = AM*seed;
  return ans;
}

/* ----------------------------------------------------------------------
   gaussian RN 
------------------------------------------------------------------------- */

double RanPark::gaussian()
{
  double first,v1,v2,rsq,fac;

  if (!save) {
    int again = 1;
    while (again) {
      v1 = 2.0*uniform()-1.0;
      v2 = 2.0*uniform()-1.0;
      rsq = v1*v1 + v2*v2;
      if (rsq < 1.0 && rsq != 0.0) again = 0;
    }
    fac = sqrt(-2.0*log(rsq)/rsq);
    second = v1*fac;
    first = v2*fac;
    save = 1;
  } else {
    first = second;
    save = 0;
  }
  return first;
}

/* ---------------------------------------------------------------------- */

void RanPark::reset(int seed_init)
{
  if (seed_init <= 0){printf("\nInvalid seed for Park random # generator\n"); exit(2);}
  seed = seed_init;
  save = 0;
}

/* ---------------------------------------------------------------------- */

int RanPark::state()
{
  return seed;
}
