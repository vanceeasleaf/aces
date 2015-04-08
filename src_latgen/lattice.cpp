#include "lattice.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include <list>
#include <map>

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

/* ----------------------------------------------------------------------
   constructor does nothing
------------------------------------------------------------------------- */
lattice::lattice()
{
  initialized = 0;
  perp_x = perp_y = perp_z = 0;
  nlayer = nucell = ntype = 0;

  name = NULL;
  attyp = NULL;
  atpos = NULL;

  layer = NULL;
  numlayer = NULL;
  h        = NULL;

  memory = new Memory;
}

/* ----------------------------------------------------------------------
   deconstructor destroy any allocated memory
------------------------------------------------------------------------- */
lattice::~lattice()
{
  if (atpos) memory->destroy(atpos);
  if (attyp) memory->destroy(attyp);
  if (name ) memory->destroy(name);
  if (layer) memory->destroy(layer);
  if (numlayer) memory->destroy(numlayer);
  if (h       ) memory->destroy(h);
  
  delete memory;
}

/* ----------------------------------------------------------------------
   Display lattice basic info
------------------------------------------------------------------------- */
void lattice::display()
{
  if (!initialized) return;
  setup();

  printf("\n"); for (int i=0; i<28; i++) printf("=");
  printf(" Lattice Info "); for (int i=0; i<28; i++) printf("="); printf("\n");
  printf("Lattice name......................: %s\n", name);
  printf("Number of atom per unit cell......: %d\n", nucell);
  printf("Number of atom types per unit cell: %d\n", ntype);
  printf("Number of layers in each unit cell: %d\n", nlayer);
  printf("Number of atoms  in each layer    :");
  for (int i=0; i<nlayer; i++) printf(" %d", numlayer[i]); printf("\n");
  printf("Expected height above each layer  :");
  for (int i=0; i<nlayer; i++) printf(" %g", h[i]); printf("\n");
  for (int i=0; i<14; i++) printf("-----");
  printf("\nLattice vectors:\n");
  for (int i=0; i<3; i++){
    for(int j=0; j<3; j++) printf("%lf ", latvec[i][j]*alat);
    printf("\n");
  }
  for (int i=0; i<14; i++) printf("-----");
  printf("\nBasis (Fractional coordinate & type):\n");
  for (int i=0; i<nucell; i++){
    printf("%lf %lf %lf %d\n", atpos[i][0], atpos[i][1], atpos[i][2], attyp[i]);
  }
  for (int i=0; i<14; i++) printf("====="); printf("\n\n");
}

/*------------------------------------------------------------------------------
 * Method to re-orient the lattice, hope to follow the rules of LAMMPS
 *----------------------------------------------------------------------------*/
void lattice::OrientLattice()
{
  double LV[3][3], box[3], oldy[3];
  box[0] = box[1] = box[2] = 0.;
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
      LV[i][j] = latvec[i][j];
      box[i] += latvec[i][j]*latvec[i][j];
      latvec[i][j] = 0.;
    }
    box[i] = sqrt(box[i]);
  }
  // new A
  latvec[0][0] = box[0];
  // new B
  for (int i=0; i<3; i++) latvec[1][0] += LV[0][i]*LV[1][i]/box[0];
  latvec[1][1] = sqrt(box[1]*box[1]-latvec[1][0]*latvec[1][0]);
  // new C
  for (int i=0; i<3; i++) latvec[2][0] += LV[0][i]*LV[2][i]/box[0];

  for (int i=0; i<3; i++) oldy[i] = LV[1][i] - latvec[1][0]*LV[0][i]/box[0];
  double yl=0.;
  for (int i=0; i<3; i++) yl += oldy[i]*oldy[i];
  latvec[2][1] = (oldy[0]*LV[2][0] + oldy[1]*LV[2][1] + oldy[2]*LV[2][2])/sqrt(yl);
  latvec[2][2] = sqrt(box[2]*box[2]-latvec[2][0]*latvec[2][0]-latvec[2][1]*latvec[2][1]);

return;
}

/*------------------------------------------------------------------------------
 * Method to count # of words in a string, without destroying the string
 *----------------------------------------------------------------------------*/
int lattice::count_words(const char *line)
{
  int n = strlen(line) + 1;
  char *copy = (char *) memory->smalloc(n*sizeof(char),"copy");
  strcpy(copy,line);

  char *ptr;
  if (ptr = strchr(copy,'#')) *ptr = '\0';

  if (strtok(copy," \t\n\r\f") == NULL) {
    memory->sfree(copy);
    return 0;
  }
  n = 1;
  while (strtok(NULL," \t\n\r\f")) n++;

  memory->sfree(copy);
  return n;
}

/*------------------------------------------------------------------------------
 * Method to setup the dependent parameters
 *----------------------------------------------------------------------------*/
void lattice::setup()
{
  if (initialized == 0) return;
  // get the layer ID
  std::list<double> zlist;
  std::list<double>::iterator it;
  std::map<double,int> zmap;

  zlist.clear();
  for (int i=0; i<nucell; i++) zlist.push_back(atpos[i][2]);
  zlist.sort(); zlist.unique();

  layer = memory->create(layer, nucell, "lattice:setup:layers");
  int il = 0;
  for (it = zlist.begin(); it != zlist.end(); it++) zmap[*it] = il++;
  for (int i=0; i<nucell; i++) layer[i] = zmap[atpos[i][2]];
  nlayer = (int) zmap.size();
  zlist.clear(); zmap.clear();

  // get the height above each layer
  if (h) memory->destroy(h);
  if (numlayer) memory->destroy(numlayer);
  h = memory->create(h, nlayer, "lattice->setup:h");
  numlayer = memory->create(numlayer, nlayer, "lattice->setup:numlayer");

  for (int i=0; i<nlayer; i++) numlayer[i] = 0;
  for (int i=0; i<nucell; i++) numlayer[layer[i]]++;

  int i0 = 0, i1 =0;
  double pos0, pos1;
  for (int i=0; i<nucell; i++) if (layer[i] == 0) {i0 = i; break;}
  pos0 = 0.;
  for (int i=0; i<3; i++) pos0 += atpos[i0][i] * latvec[i][2];

  for (int il=1; il<nlayer; il++){
    for (int i=0; i<nucell; i++) if (layer[i] == il) {i1 = i; break;}
    pos1 = 0.;
    for (int i=0; i<3; i++) pos1 += atpos[i1][i] * latvec[i][2];
    h[il-1] = (pos1 - pos0)*alat;
    pos0 = pos1;
  }

  pos1 = 0.;
  for (int i=0; i<3; i++) pos1 += (atpos[i0][i]+double(i/2)) * latvec[i][2];
  h[nlayer-1] = (pos1-pos0)*alat;

  // get the length of each vector
  lx = sqrt(DotProd(latvec[0], latvec[0]))*alat;
  ly = sqrt(DotProd(latvec[1], latvec[1]))*alat;
  lz = sqrt(DotProd(latvec[2], latvec[2]))*alat;

  // get the norm vector of each pair
  double Nxy[3], Nxz[3], Nyz[3];
  Cross(latvec[0], latvec[1], Nxy);
  Cross(latvec[0], latvec[2], Nxz);
  Cross(latvec[1], latvec[2], Nyz);
  double Sxy = sqrt(DotProd(&Nxy[0],&Nxy[0]));
  double Sxz = sqrt(DotProd(&Nxz[0],&Nxz[0]));
  double Syz = sqrt(DotProd(&Nyz[0],&Nyz[0]));
  double vola = DotProd(&Nxy[0], latvec[2])*alat;
  hx = vola/Syz;
  hy = vola/Sxz;
  hz = vola/Sxy;

  // check if one basic vector is perpendicular to the remaing two
  double AdotB = DotProd(latvec[0], latvec[1]);
  double AdotC = DotProd(latvec[0], latvec[2]);
  double BdotC = DotProd(latvec[1], latvec[2]);

  if ( (AdotB+AdotC) < 1.e-8 ) perp_x = 1;
  if ( (AdotB+BdotC) < 1.e-8 ) perp_y = 1;
  if ( (AdotC+BdotC) < 1.e-8 ) perp_z = 1;

return;
}

/*------------------------------------------------------------------------------
 * Method to rotate the lattice. angles are in unit of 2*pi
 *----------------------------------------------------------------------------*/
void lattice::RotateLattice(double *angles, double rotated[][3])
{
  const double tPI = 8.*atan(1.);
  double cosa = cos(angles[0]*tPI);
  double sina = sin(angles[0]*tPI);
  double cosb = cos(angles[1]*tPI);
  double sinb = sin(angles[1]*tPI);
  double cosg = cos(angles[2]*tPI);
  double sing = sin(angles[2]*tPI);

  double rot[3][3];
  rot[0][0] = cosb*cosg;
  rot[0][1] = cosb*sing;
  rot[0][2] = -sinb;
  rot[1][0] = sina * sinb * cosg - cosa * sing;
  rot[1][1] = sina * sinb * sing + cosa * cosg;
  rot[1][2] = sina * cosb;
  rot[2][0] = cosa * sinb * cosg + sina * sing;
  rot[2][1] = cosa * sinb * sing - sina * cosg;
  rot[2][2] = cosa * cosb;

  for (int i=0; i<3; i++)
  for (int j=0; j<3; j++) rotated[i][j] = 0.;

  for (int i=0; i<3; i++)
  for (int j=0; j<3; j++)
  for (int k=0; k<3; k++) rotated[i][j] += rot[i][k] * latvec[k][j];

return;
}

/*------------------------------------------------------------------------------
 * Method to rotate the lattice. angles are in unit of 2*pi
 *----------------------------------------------------------------------------*/
void lattice::RotateLattice(double *angles)
{
  const double tPI = 8.*atan(1.);
  double cosa = cos(angles[0]*tPI);
  double sina = sin(angles[0]*tPI);
  double cosb = cos(angles[1]*tPI);
  double sinb = sin(angles[1]*tPI);
  double cosg = cos(angles[2]*tPI);
  double sing = sin(angles[2]*tPI);

  double rot[3][3];
  rot[0][0] = cosb*cosg;
  rot[0][1] = cosb*sing;
  rot[0][2] = -sinb;
  rot[1][0] = sina * sinb * cosg - cosa * sing;
  rot[1][1] = sina * sinb * sing + cosa * cosg;
  rot[1][2] = sina * cosb;
  rot[2][0] = cosa * sinb * cosg + sina * sing;
  rot[2][1] = cosa * sinb * sing - sina * cosg;
  rot[2][2] = cosa * cosb;

  double rotated[3][3];
  for (int i=0; i<3; i++)
  for (int j=0; j<3; j++) rotated[i][j] = 0.;

  for (int i=0; i<3; i++)
  for (int j=0; j<3; j++)
  for (int k=0; k<3; k++) rotated[i][j] += rot[i][k] * latvec[k][j];

  for (int i=0; i<3; i++)
  for (int j=0; j<3; j++) latvec[i][j] = rotated[i][j];

return;
}

/*------------------------------------------------------------------------------
 * Private method to get the cross product of two vectors; C[3] = A[3] X B[3]
 *----------------------------------------------------------------------------*/
void lattice::Cross(double *A, double *B, double *C)
{
   C[0] = A[1] * B[2] - A[2] * B[1];
   C[1] = A[2] * B[0] - A[0] * B[2];
   C[2] = A[0] * B[1] - A[1] * B[0];

return;
}

/*------------------------------------------------------------------------------
 * Private method to get the dot product of two vectors of dim 3
 *----------------------------------------------------------------------------*/
double lattice::DotProd(double *A, double *B)
{
return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

/*----------------------------------------------------------------------------*/
