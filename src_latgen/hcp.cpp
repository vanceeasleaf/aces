#include "hcp.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

#define MAXLINE 256

using namespace std;

/*-----------------------------------------------------------------------
 * To select the orientation of the lattice
 *----------------------------------------------------------------------- */
HCP::HCP() : lattice()
{
  char str[MAXLINE];
  alat = 1.; ca = sqrt(8./3.);
  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Please input the lattice constant of the HCP lattice [1]:");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) alat = atof(strtok(str, " \t\n\r\f"));
  printf("Please input the value of c/a ratio or c (negative) [1.633]: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) ca = atof(strtok(str, " \t\n\r\f"));
  if (ca < 0.) ca = -ca/alat;
  printf("The lattice constants of your HCP: a = %g, c/a = %g.\n", alat, ca);

  int orient = 1;
  printf("Please select the orientation of the HCP lattice:\n");
  printf("   1. (001);         4. (-110); \n");
  printf("   2. (100);         5. Graphene;\n");
  printf("   3. (110);         6. Graphite;\n");
  printf("Your choice [%d]: ", orient);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) orient = atoi(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d", orient);
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  
  // initialize according to orientation
  initialized = 0;
  switch (orient){
  case 1:
    HCP001();
    break;
  case 2:
    HCP100();
    break;
  case 3:
    HCP110();
    break;
  case 4:
    HCPm10();
    break;
  case 5:
    Graphene();
    break;
  case 6:
    Graphite();
    break;
  default:
    break;
  }

}

/*-----------------------------------------------------------------------
 * Deconstructor does nothing
 *----------------------------------------------------------------------- */
HCP::~HCP()
{

}

/*-----------------------------------------------------------------------
 * Initialize for (001) orientation
 *----------------------------------------------------------------------- */
void HCP::HCP001()
{
  char str[MAXLINE];
  int surftype = 5;
  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Please select the type of HCP(001) surface:\n");
  printf("   1. U along x,  60 degree;\n");
  printf("   2. V along y,  60 degree;\n");
  printf("   3. U along x, 120 degree;\n");
  printf("   4. V along y, 120 degree;\n");
  printf("   5. Rectangle, long along x;\n");
  printf("   6. Rectangle, long along y;\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = atoi(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d", surftype);
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }
  name = memory->create(name,9,"HCP:name");
  strcpy(name, "HCP(001)");

  // initialize according to surface type
  switch (surftype){
  case 1:
    nucell = 2;
    ntype  = 1;
    
    latvec[0][0] = 1.;
    latvec[1][0] = 0.5;
    latvec[1][1] = sqrt(0.75);
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "HCP001_atpos");
    attyp = memory->create(attyp, nucell, "HCP:attyp");
    
    for (int i=0; i<nucell; i++) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 1./3.;
    atpos[1][1] = 1./3.;
    atpos[1][2] = 0.5;

    initialized = 1;
    break;
  case 2:
    nucell = 2;
    ntype  = 1;
    
    latvec[0][0] = sqrt(0.75);
    latvec[0][1] = 0.5;
    latvec[1][1] = 1.;
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "HCP001_atpos");
    attyp = memory->create(attyp, nucell, "HCP:attyp");
    
    for (int i=0; i<nucell; i++) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 1./3.;
    atpos[1][1] = 1./3.;
    atpos[1][2] = 0.5;

    initialized = 1;
    break;
  case 3:
    nucell = 2;
    ntype  = 1;
    
    latvec[0][0] =  1.;
    latvec[1][0] = -0.5;
    latvec[1][1] =  sqrt(0.75);
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "HCP001_atpos");
    attyp = memory->create(attyp, nucell, "HCP:attyp");
    
    for (int i=0; i<nucell; i++) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 2./3.;
    atpos[1][1] = 1./3.;
    atpos[1][2] = 0.5;

    initialized = 1;
    break;
  case 4:
    nucell = 2;
    ntype  = 1;
    
    latvec[0][0] = sqrt(0.75);
    latvec[0][1] = -0.5;
    latvec[1][1] = 1.;
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "HCP001_atpos");
    attyp = memory->create(attyp, nucell, "HCP:attyp");
    
    for (int i=0; i<nucell; i++) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 2./3.;
    atpos[1][1] = 1./3.;
    atpos[1][2] = 0.5;

    initialized = 1;
    break;
  case 5:
    nucell = 4;
    ntype  = 1;

    latvec[0][0] = sqrt(3.0);
    latvec[1][1] = 1.;
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "HCP001_atpos");
    attyp = memory->create(attyp, nucell, "HCP:attyp");

    for (int i=0; i<nucell; i++) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;

    atpos[2][0] = 1./6.;
    atpos[2][1] = 0.5;
    atpos[2][2] = 0.5;

    atpos[3][0] = 2./3.;
    atpos[3][1] = 0.0;
    atpos[3][2] = 0.5;

    initialized = 1;
    break;
  case 6:
    nucell = 4;
    ntype  = 1;

    latvec[0][0] = 1.;
    latvec[1][1] = sqrt(3.0);
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "HCP001_atpos");
    attyp = memory->create(attyp, nucell, "HCP:attyp");

    for (int i=0; i<nucell; i++) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;

    atpos[2][0] = 0.5;
    atpos[2][1] = 1./6.;
    atpos[2][2] = 0.5;

    atpos[3][0] = 0.0;
    atpos[3][1] = 2./3.;
    atpos[3][2] = 0.5;

    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/*-----------------------------------------------------------------------
 * Initialize for (110) orientation
 *----------------------------------------------------------------------- */
void HCP::HCP100()
{
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }
  name = memory->create(name,9,"HCP:name");
  strcpy(name, "HCP(100)");

  // initialize
  nucell = 4;
  ntype  = 1;
    
  latvec[0][0] = 1.;
  latvec[1][1] = ca;
  latvec[2][2] = sqrt(3.);

  atpos = memory->create(atpos, nucell, 3, "HCP100_atpos");
  attyp = memory->create(attyp, nucell, "HCP:attyp");
    
  for (int i=0; i<nucell; i++) attyp[i] = 1;
  atpos[0][0] = 0.50;
  atpos[0][1] = 0.75;
  atpos[0][2] = 0.00;

  atpos[1][0] = 0.50;
  atpos[1][1] = 0.25;
  atpos[1][2] = 1./3.;

  atpos[2][0] = 0.00;
  atpos[2][1] = 0.75;
  atpos[2][2] = 0.50;

  atpos[3][0] = 0.00;
  atpos[3][1] = 0.25;
  atpos[3][2] = 5./6.;

  initialized = 1;

return;
}

/*-----------------------------------------------------------------------
 * Initialize for (110) orientation
 *----------------------------------------------------------------------- */
void HCP::HCP110()
{
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }
  name = memory->create(name,9,"HCP:name");
  strcpy(name, "HCP(110)");

  // initialize
  nucell = 4;
  ntype  = 1;
    
  latvec[0][0] = ca;
  latvec[1][1] = sqrt(3.);
  latvec[2][2] = 1.;

  atpos = memory->create(atpos, nucell, 3, "HCP110_atpos");
  attyp = memory->create(attyp, nucell, "HCP:attyp");
    
  for (int i=0; i<nucell; i++) attyp[i] = 1;
  atpos[0][0] = 0.25;
  atpos[0][1] = 1./3.;
  atpos[0][2] = 0.00;

  atpos[1][0] = 0.75;
  atpos[1][1] = 2./3.;
  atpos[1][2] = 0.;

  atpos[2][0] = 0.75;
  atpos[2][1] = 1./6.;
  atpos[2][2] = 0.50;

  atpos[3][0] = 0.25;
  atpos[3][1] = 5./6.;
  atpos[3][2] = 0.50;

  initialized = 1;

return;
}

/*-----------------------------------------------------------------------
 * Initialize for (-110) orientation
 *----------------------------------------------------------------------- */
void HCP::HCPm10()
{
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }
  name = memory->create(name,10,"HCP:name");
  strcpy(name, "HCP(-110)");

  // initialize
  nucell = 4;
  ntype  = 1;
    
  latvec[0][0] = 1.;
  latvec[1][1] = ca;
  latvec[2][2] = sqrt(3.);

  atpos = memory->create(atpos, nucell, 3, "HCPm10_atpos");
  attyp = memory->create(attyp, nucell, "HCP:attyp");
    
  for (int i=0; i<nucell; i++) attyp[i] = 1;
  atpos[0][0] = 0.50;
  atpos[0][1] = 0.25;
  atpos[0][2] = 0.00;

  atpos[1][0] = 0.00;
  atpos[1][1] = 0.75;
  atpos[1][2] = 1./6.;

  atpos[2][0] = 0.00;
  atpos[2][1] = 0.25;
  atpos[2][2] = 0.50;

  atpos[3][0] = 0.50;
  atpos[3][1] = 0.75;
  atpos[3][2] = 2./3.;

  initialized = 1;

return;
}

/*-----------------------------------------------------------------------
 * Graphene with different orientations
 *----------------------------------------------------------------------- */
void HCP::Graphene()
{
  char str[MAXLINE];
  int surftype = 5;
  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("   1. Primitive, U along x,  60 degree;\n");
  printf("   2. Primitive, V along y,  60 degree;\n");
  printf("   3. Primitive, U along x, 120 degree;\n");
  printf("   4. Primitive, V along y, 120 degree;\n");
  printf("   5. Rectangle, long along x;\n");
  printf("   6. Rectangle, long along y;\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = atoi(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d", surftype);
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }
  name = memory->create(name,9,"HCP:name");
  strcpy(name, "Graphene");

  // initialize according to surface type
  switch (surftype){
  case 1:
    nucell = 2;
    ntype  = 1;
    
    latvec[0][0] = 1.;
    latvec[1][0] = 0.5;
    latvec[1][1] = sqrt(0.75);
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "Graphene_atpos");
    attyp = memory->create(attyp, nucell, "HCP:attyp");
    
    for (int i=0; i<nucell; i++) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 1./3.;
    atpos[1][1] = 1./3.;
    atpos[1][2] = 0.;

    initialized = 1;
    break;
  case 2:
    nucell = 2;
    ntype  = 1;
    
    latvec[0][0] = sqrt(0.75);
    latvec[0][1] = 0.5;
    latvec[1][1] = 1.;
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "Graphene_atpos");
    attyp = memory->create(attyp, nucell, "HCP:attyp");
    
    for (int i=0; i<nucell; i++) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 1./3.;
    atpos[1][1] = 1./3.;
    atpos[1][2] = 0.;

    initialized = 1;
    break;
  case 3:
    nucell = 2;
    ntype  = 1;
    
    latvec[0][0] =  1.;
    latvec[1][0] = -0.5;
    latvec[1][1] =  sqrt(0.75);
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "Graphene_atpos");
    attyp = memory->create(attyp, nucell, "HCP:attyp");
    
    for (int i=0; i<nucell; i++) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 1./3.;
    atpos[1][1] = 2./3.;
    atpos[1][2] = 0.;

    initialized = 1;
    break;
  case 4:
    nucell = 2;
    ntype  = 1;
    
    latvec[0][0] = sqrt(0.75);
    latvec[0][1] = -0.5;
    latvec[1][1] = 1.;
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "Graphene_atpos");
    attyp = memory->create(attyp, nucell, "HCP:attyp");
    
    for (int i=0; i<nucell; i++) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 1./3.;
    atpos[1][1] = 2./3.;
    atpos[1][2] = 0.;

    initialized = 1;
    break;
  case 5:
    nucell = 4;
    ntype  = 1;

    latvec[0][0] = sqrt(3.0);
    latvec[1][1] = 1.;
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "Graphene_atpos");
    attyp = memory->create(attyp, nucell, "HCP:attyp");

    for (int i=0; i<nucell; i++) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 1./3.;
    atpos[1][1] = 0.;
    atpos[1][2] = 0.;

    atpos[2][0] = 0.5;
    atpos[2][1] = 0.5;
    atpos[2][2] = 0.;

    atpos[3][0] = 5./6.;
    atpos[3][1] = 0.5;
    atpos[3][2] = 0.;

    initialized = 1;
    break;
  case 6:
    nucell = 4;
    ntype  = 1;

    latvec[0][0] = 1.;
    latvec[1][1] = sqrt(3.0);
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "Graphene_atpos");
    attyp = memory->create(attyp, nucell, "HCP:attyp");

    for (int i=0; i<nucell; i++) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 0.;
    atpos[1][1] = 1./3.;
    atpos[1][2] = 0.;

    atpos[2][0] = 0.5;
    atpos[2][1] = 0.5;
    atpos[2][2] = 0.;

    atpos[3][0] = 0.5;
    atpos[3][1] = 5./6.;
    atpos[3][2] = 0.;

    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/*-----------------------------------------------------------------------
 * Graphite with different orientations
 *----------------------------------------------------------------------- */
void HCP::Graphite()
{
  char str[MAXLINE];
  int surftype = 5;
  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("   1. Primitive, U along x,  60 degree;\n");
  printf("   2. Primitive, V along y,  60 degree;\n");
  printf("   3. Primitive, U along x, 120 degree;\n");
  printf("   4. Primitive, V along y, 120 degree;\n");
  printf("   5. Rectangle, long along x;\n");
  printf("   6. Rectangle, long along y;\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = atoi(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d", surftype);
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }
  name = memory->create(name,9,"HCP:name");
  strcpy(name, "Graphite");

  // initialize according to surface type
  switch (surftype){
  case 1:
    nucell = 4;
    ntype  = 1;
    
    latvec[0][0] = 1.;
    latvec[1][0] = 0.5;
    latvec[1][1] = sqrt(0.75);
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "Graphite_atpos");
    attyp = memory->create(attyp, nucell, "HCP:attyp");
    
    for (int i=0; i<nucell; i++) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 1./3.;
    atpos[1][1] = 1./3.;
    atpos[1][2] = 0.;

    atpos[2][0] = 0.;
    atpos[2][1] = 0.;
    atpos[2][2] = 0.5;

    atpos[3][0] = 2./3.;
    atpos[3][1] = 2./3.;
    atpos[3][2] = 0.5;

    initialized = 1;
    break;
  case 2:
    nucell = 4;
    ntype  = 1;
    
    latvec[0][0] = sqrt(0.75);
    latvec[0][1] = 0.5;
    latvec[1][1] = 1.;
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "Graphite_atpos");
    attyp = memory->create(attyp, nucell, "HCP:attyp");
    
    for (int i=0; i<nucell; i++) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 1./3.;
    atpos[1][1] = 1./3.;
    atpos[1][2] = 0.;

    atpos[2][0] = 0.;
    atpos[2][1] = 0.;
    atpos[2][2] = 0.5;

    atpos[3][0] = 2./3.;
    atpos[3][1] = 2./3.;
    atpos[3][2] = 0.5;

    initialized = 1;
    break;
  case 3:
    nucell = 4;
    ntype  = 1;
    
    latvec[0][0] =  1.;
    latvec[1][0] = -0.5;
    latvec[1][1] =  sqrt(0.75);
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "Graphite_atpos");
    attyp = memory->create(attyp, nucell, "HCP:attyp");
    
    for (int i=0; i<nucell; i++) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 1./3.;
    atpos[1][1] = 2./3.;
    atpos[1][2] = 0.;

    atpos[2][0] = 0.;
    atpos[2][1] = 0.;
    atpos[2][2] = 0.5;

    atpos[3][0] = 2./3.;
    atpos[3][1] = 1./3.;
    atpos[3][2] = 0.5;

    initialized = 1;
    break;
  case 4:
    nucell = 4;
    ntype  = 1;
    
    latvec[0][0] = sqrt(0.75);
    latvec[0][1] = -0.5;
    latvec[1][1] = 1.;
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "Graphite_atpos");
    attyp = memory->create(attyp, nucell, "HCP:attyp");
    
    for (int i=0; i<nucell; i++) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 1./3.;
    atpos[1][1] = 2./3.;
    atpos[1][2] = 0.;

    atpos[2][0] = 0.;
    atpos[2][1] = 0.;
    atpos[2][2] = 0.5;

    atpos[3][0] = 2./3.;
    atpos[3][1] = 1./3.;
    atpos[3][2] = 0.5;

    initialized = 1;
    break;
  case 5:
    nucell = 8;
    ntype  = 1;

    latvec[0][0] = sqrt(3.0);
    latvec[1][1] = 1.;
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "Graphite_atpos");
    attyp = memory->create(attyp, nucell, "HCP:attyp");

    for (int i=0; i<nucell; i++) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 1./3.;
    atpos[1][1] = 0.;
    atpos[1][2] = 0.;

    atpos[2][0] = 0.5;
    atpos[2][1] = 0.5;
    atpos[2][2] = 0.;

    atpos[3][0] = 5./6.;
    atpos[3][1] = 0.5;
    atpos[3][2] = 0.;

    atpos[4][0] = 0.;
    atpos[4][1] = 0.;
    atpos[4][2] = 0.5;

    atpos[5][0] = 2./3.;
    atpos[5][1] = 0.;
    atpos[5][2] = 0.5;

    atpos[6][0] = 0.5;
    atpos[6][1] = 0.5;
    atpos[6][2] = 0.5;

    atpos[7][0] = 1./6.;
    atpos[7][1] = 0.5;
    atpos[7][2] = 0.5;

    initialized = 1;
    break;
  case 6:
    nucell = 8;
    ntype  = 1;

    latvec[0][0] = 1.;
    latvec[1][1] = sqrt(3.0);
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "Graphite_atpos");
    attyp = memory->create(attyp, nucell, "HCP:attyp");

    for (int i=0; i<nucell; i++) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 0.;
    atpos[1][1] = 1./3.;
    atpos[1][2] = 0.;

    atpos[2][0] = 0.5;
    atpos[2][1] = 0.5;
    atpos[2][2] = 0.;

    atpos[3][0] = 0.5;
    atpos[3][1] = 5./6.;
    atpos[3][2] = 0.;

    atpos[4][0] = 0.;
    atpos[4][1] = 0.;
    atpos[4][2] = 0.5;

    atpos[5][0] = 0.;
    atpos[5][1] = 2./3.;
    atpos[5][2] = 0.5;

    atpos[6][0] = 0.5;
    atpos[6][1] = 0.5;
    atpos[6][2] = 0.5;

    atpos[7][0] = 0.5;
    atpos[7][1] = 1./6.;
    atpos[7][2] = 0.5;

    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/*-----------------------------------------------------------------------*/
