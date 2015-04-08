#ifndef BCC_H
#define BCC_H

#include "lattice.h"

class BCC : public lattice {
public:
  
  BCC();
  ~BCC();

private:
  void BCC001();
  void BCC110();
  void BCC111();
  void Primitive();
};
#endif
