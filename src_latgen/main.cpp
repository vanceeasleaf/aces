#include "driver.h"

int main (int narg, char **arg)
{
  Driver *driver = new Driver;

  driver->modify();
  driver->write();

  delete driver;

return 0;
}
