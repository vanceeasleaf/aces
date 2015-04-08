#ifndef RANPARK_H
#define RANPARK_H

class RanPark{
 public:
  RanPark(int);
  double uniform();
  double gaussian();
  void reset(int);
  int state();

 private:
  int seed,save;
  double second;
};

#endif
