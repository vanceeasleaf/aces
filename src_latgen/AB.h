#ifndef AB_H
#define AB_H

#include "lattice.h"

class AB : public lattice {
public:
  
  AB();
  ~AB();

private:
  void AB_B1();
  void AB_B2();
  void AB_B3();
  void AB_B4();
  void AB_L10();
  void AB_NiAs();
  void AB_Perov();

  double ba, ca;
};
#endif
