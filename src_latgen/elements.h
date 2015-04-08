#ifndef ELEMENTS_H
#define ELEMENTS_H

#include <string>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

class ChemElements {
public:
  int Name2Num(const char *);
  void Num2Name(const int, char *);
  double Name2Mass(const char *);
  double Num2Mass(const int);

private:
  static const int NumMax;
  static const double weight[];
  static const char   symbol[][3];
};

#endif
