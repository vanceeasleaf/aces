#ifndef FCC_H
#define FCC_H

#include "lattice.h"

class FCC : public lattice {
public:
  
  FCC();
  ~FCC();

private:
  void FCC001();
  void FCC110();
  void FCC111();
  void Primitive();
  void DiamondPrim();
  void Diamond001();
  void Diamond110();
  void Diamond111();
};
#endif
