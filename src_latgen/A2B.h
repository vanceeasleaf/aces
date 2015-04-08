#ifndef A2B_H
#define A2B_H

#include "lattice.h"

class A2B : public lattice {
public:
  
  A2B();
  ~A2B();

private:
  void A2B_C1();
  void A2B_C15();
  void A2B_C32();

  double ca;
  int ip1, ip2;
};
#endif
