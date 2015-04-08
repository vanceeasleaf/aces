#ifndef HCP_H
#define HCP_H

#include "lattice.h"

class HCP : public lattice {
public:
  
  HCP();
  ~HCP();

private:
  void HCP001();
  void HCP100();
  void HCP110();
  void HCPm10();
  void Graphene();
  void Graphite();
  double ca;
};
#endif
