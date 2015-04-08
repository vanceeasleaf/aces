#ifndef A3B_H
#define A3B_H

#include "lattice.h"

class A3B : public lattice {
public:
  
  A3B();
  ~A3B();

private:
  void A3B_A15();
  void A3B_D019();
  void A3B_D022();
  void A3B_D03();
  void A3B_D09();
  void A3B_L12();
  void A3B_L60();

  double ca;
  int ip1, ip2;
};
#endif
