#include "fcc.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

#define MAXLINE 256

using namespace std;

/* ----------------------------------------------------------------------
   To select the orientation of the lattice
------------------------------------------------------------------------- */
FCC::FCC() : lattice()
{
  char str[MAXLINE];
  alat = 1.;
  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Please input the lattice constant of the FCC/Diamond lattice [1.]: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) alat = atof(strtok(str, " \t\n\r\f"));

  int orient = 3;
  printf("Please selection the orientation of the FCC lattice:\n");
  printf("   1. (001);                5. Diamond primitive;\n");
  printf("   2. (110);                6. Diamond (001);\n");
  printf("   3. (111);                7. Diamond (110);\n");
  printf("   4. Primitive cell;       8. Diamond (111);\n");
  for (int i=0; i<14; i++) printf("-----");
  printf("\nYour choice [%d]: ", orient);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) orient = atoi(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d", orient);
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  
  // initialize according to orientation
  initialized = 0;
  switch (orient){
  case 1:
    FCC001();
    break;
  case 2:
    FCC110();
    break;
  case 3:
    FCC111();
    break;
  case 4:
    Primitive();
    break;
  case 5:
    DiamondPrim();
    break;
  case 6:
    Diamond001();
    break;
  case 7:
    Diamond110();
    break;
  case 8:
    Diamond111();
    break;
  default:
    break;
  }

}

/* ----------------------------------------------------------------------
   Deconstructor does nothing
------------------------------------------------------------------------- */
FCC::~FCC()
{

}

/* ----------------------------------------------------------------------
   Initialize for (001) orientation
------------------------------------------------------------------------- */
void FCC::FCC001()
{
  char str[MAXLINE];
  int surftype = 2;
  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Please selection the type of FCC(001) surface:\n");
  printf("   1. primitive, horizental orientation;\n");
  printf("   2. conventional orientation;\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = atoi(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d", surftype);
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }
  name = memory->create(name,9,"FCC:name");
  strcpy(name, "FCC(001)");

  // initialize according to surface type
  switch (surftype){
  case 1:
    nucell = 2;
    ntype  = 1;
    
    latvec[0][0] = 1./sqrt(2.);
    latvec[1][1] = 1./sqrt(2.);
    latvec[2][2] = 1.;

    atpos = memory->create(atpos, nucell, 3, "FCC001_atpos");
    attyp = memory->create(attyp, nucell, "FCC:attyp");
  
    
    for (int i=0; i<nucell; i++) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.5;

    initialized = 1;
    break;
  case 2:
    nucell = 4;
    ntype  = 1;
    
    latvec[0][0] = 1.;
    latvec[1][1] = 1.;
    latvec[2][2] = 1.;

    atpos = memory->create(atpos, nucell, 3, "FCC001_atpos");
    attyp = memory->create(attyp, nucell, "FCC:attyp");
    
    for (int i=0; i<nucell; i++) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;

    atpos[2][0] = 0.;
    atpos[2][1] = 0.5;
    atpos[2][2] = 0.5;

    atpos[3][0] = 0.5;
    atpos[3][1] = 0.;
    atpos[3][2] = 0.5;

    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/* ----------------------------------------------------------------------
   Initialize for (110) orientation
------------------------------------------------------------------------- */
void FCC::FCC110()
{
  char str[MAXLINE];
  int surftype = 1;
  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Please selection the type of FCC(110) surface:\n");
  printf("   1. orthogonal, long side along x\n");
  printf("   2. orthogonal, long side along y\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = atoi(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d", surftype);
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }
  name = memory->create(name,9,"FCC:name");
  strcpy(name, "FCC(110)");

  // initialize according to surface type
  switch (surftype){
  case 2:
    nucell = 2;
    ntype  = 1;
    
    latvec[0][0] = 0.5*sqrt(2.);
    latvec[1][1] = 1.;
    latvec[2][2] = 0.5*sqrt(2.);

    atpos = memory->create(atpos, nucell, 3, "FCC110_atpos");
    attyp = memory->create(attyp, nucell, "FCC:attyp");
    
    for (int i=0; i<nucell; i++) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.5;

    initialized = 1;
    break;
  case 1:
    nucell = 2;
    ntype  = 1;
    
    latvec[0][0] = 1.;
    latvec[1][1] = sqrt(2.)*0.5;
    latvec[2][2] = sqrt(2.)*0.5;

    atpos = memory->create(atpos, nucell, 3, "FCC110_atpos");
    attyp = memory->create(attyp, nucell, "FCC:attyp");
    
    for (int i=0; i<nucell; i++) attyp[i] = 1;

    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.5;

    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/* ----------------------------------------------------------------------
   Initialize for (111) orientation
------------------------------------------------------------------------- */
void FCC::FCC111()
{
  char str[MAXLINE];
  int surftype = 5;
  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Please selection the type of FCC(111) surface:\n");
  printf("   1. hexgonal U along x, 60 deg;\n");
  printf("   2. hexgonal V along y, 60 deg;\n");
  printf("   3. hexgonal U along x, 120 deg;\n");
  printf("   4. hexgonal V along y, 120 deg;\n");
  printf("   5. orthogonal long side along x;\n");
  printf("   6. orthogonal long side along y;\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = atoi(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d", surftype);
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }
  name = memory->create(name,9,"FCC:name");
  strcpy(name, "FCC(111)");

  // initialize according to surface type
  switch (surftype){
  case 1:
    nucell = 3;
    ntype  = 1;
    
    latvec[0][0] = sqrt(0.5);
    latvec[1][0] = sqrt(0.125);
    latvec[1][1] = sqrt(0.375);
    latvec[2][2] = sqrt(3.);

    atpos = memory->create(atpos, nucell, 3, "FCC111_atpos");
    attyp = memory->create(attyp, nucell, "FCC:attyp");
    
    for (int i=0; i<nucell; i++) attyp[i] = 1;

    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 1./3.;
    atpos[1][1] = 1./3.;
    atpos[1][2] = 1./3.;

    atpos[2][0] = 2./3.;
    atpos[2][1] = 2./3.;
    atpos[2][2] = 2./3.;

    initialized = 1;
    break;
  case 2:
    nucell = 3;
    ntype  = 1;

    latvec[0][0] = sqrt(0.375);
    latvec[0][1] = sqrt(0.125);
    latvec[1][1] = sqrt(0.5);
    latvec[2][2] = sqrt(3.);

    atpos = memory->create(atpos, nucell, 3, "FCC111_atpos");
    attyp = memory->create(attyp, nucell, "FCC:attyp");

    for (int i=0; i<nucell; i++) attyp[i] = 1;

    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 1./3.;
    atpos[1][1] = 1./3.;
    atpos[1][2] = 1./3.;

    atpos[2][0] = 2./3.;
    atpos[2][1] = 2./3.;
    atpos[2][2] = 2./3.;

    initialized = 1;
    break;
  case 3:
    nucell = 3;
    ntype  = 1;

    latvec[0][0] =  sqrt(0.5);
    latvec[1][0] = -sqrt(0.125);
    latvec[1][1] =  sqrt(0.375);
    latvec[2][2] =  sqrt(3.);

    atpos = memory->create(atpos, nucell, 3, "FCC111_atpos");
    attyp = memory->create(attyp, nucell, "FCC:attyp");

    for (int i=0; i<nucell; i++) attyp[i] = 1;

    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 2./3.;
    atpos[1][1] = 1./3.;
    atpos[1][2] = 1./3.;

    atpos[2][0] = 1./3.;
    atpos[2][1] = 2./3.;
    atpos[2][2] = 2./3.;

    initialized = 1;
    break;
  case 4:
    nucell = 3;
    ntype  = 1;

    latvec[0][0] =  sqrt(0.375);
    latvec[0][1] = -sqrt(0.125);
    latvec[1][1] =  sqrt(0.5);
    latvec[2][2] =  sqrt(3.);

    atpos = memory->create(atpos, nucell, 3, "FCC111_atpos");
    attyp = memory->create(attyp, nucell, "FCC:attyp");

    for (int i=0; i<nucell; i++) attyp[i] = 1;

    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 1./3.;
    atpos[1][1] = 2./3.;
    atpos[1][2] = 1./3.;

    atpos[2][0] = 2./3.;
    atpos[2][1] = 1./3.;
    atpos[2][2] = 2./3.;

    initialized = 1;
    break;
  case 5:
    nucell = 6;
    ntype  = 1;

    latvec[0][0] =  sqrt(1.5);
    latvec[1][1] =  sqrt(0.5);
    latvec[2][2] =  sqrt(3.);

    atpos = memory->create(atpos, nucell, 3, "FCC111_atpos");
    attyp = memory->create(attyp, nucell, "FCC:attyp");

    for (int i=0; i<nucell; i++) attyp[i] = 1;

    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;

    atpos[2][0] = 1./6.;
    atpos[2][1] = 0.5;
    atpos[2][2] = 1./3.;

    atpos[3][0] = 2./3.;
    atpos[3][1] = 1.;
    atpos[3][2] = 1./3.;

    atpos[4][0] = 1./3.;
    atpos[4][1] = 0.;
    atpos[4][2] = 2./3.;

    atpos[5][0] = 5./6.;
    atpos[5][1] = 0.5;
    atpos[5][2] = 2./3.;

    initialized = 1;
    break;
  case 6:
    nucell = 6;
    ntype  = 1;

    latvec[0][0] = sqrt(0.5);
    latvec[1][1] = sqrt(1.5);
    latvec[2][2] = sqrt(3.);

    atpos = memory->create(atpos, nucell, 3, "FCC111_atpos");
    attyp = memory->create(attyp, nucell, "FCC:attyp");

    for (int i=0; i<nucell; i++) attyp[i] = 1;

    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;

    atpos[2][0] = 0.5;
    atpos[2][1] = 1./6.;
    atpos[2][2] = 1./3.;

    atpos[3][0] = 1.;
    atpos[3][1] = 2./3.;
    atpos[3][2] = 1./3.;

    atpos[4][0] = 0.;
    atpos[4][1] = 1./3.;
    atpos[4][2] = 2./3.;

    atpos[5][0] = 0.5;
    atpos[5][1] = 5./6.;
    atpos[5][2] = 2./3.;

    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/* ----------------------------------------------------------------------
   Initialize for (001) orientation
------------------------------------------------------------------------- */
void FCC::Primitive()
{
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }
  name = memory->create(name,9,"FCC:name");
  strcpy(name, "FCC-prim");

  nucell = 1;
  ntype  = 1;
  
  latvec[0][1] =  0.5;
  latvec[0][2] =  0.5;
  latvec[1][0] =  0.5;
  latvec[1][2] =  0.5;
  latvec[2][0] =  0.5;
  latvec[2][1] =  0.5;

  atpos = memory->create(atpos,nucell,3,"Primitive:atpos");
  attyp = memory->create(attyp,nucell,"Primitive:attyp");
  
  for (int i=0; i<nucell; i++) attyp[i] = 1;
  atpos[0][0] = 0.;
  atpos[0][1] = 0.;
  atpos[0][2] = 0.;

  initialized = 1;
return;
}

/* ----------------------------------------------------------------------
   Initialize Diamond primitive cell
------------------------------------------------------------------------- */
void FCC::DiamondPrim()
{
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }
  name = memory->create(name,9,"FCC:name");
  strcpy(name, "Dia-prim");

  nucell = 2;
  ntype  = 1;
  
  latvec[0][1] =  0.5;
  latvec[0][2] =  0.5;
  latvec[1][0] =  0.5;
  latvec[1][2] =  0.5;
  latvec[2][0] =  0.5;
  latvec[2][1] =  0.5;

  atpos = memory->create(atpos, nucell, 3, "Primitive:atpos");
  attyp = memory->create(attyp, nucell, "Primitive:attyp");
  
  for (int i=0; i<nucell; i++) attyp[i] = 1;
  atpos[0][0] = 0.;
  atpos[0][1] = 0.;
  atpos[0][2] = 0.;

  atpos[1][0] = 0.25;
  atpos[1][1] = 0.25;
  atpos[1][2] = 0.25;

  initialized = 1;
return;
}

/* ----------------------------------------------------------------------
   Initialize Diamond (001)
------------------------------------------------------------------------- */
void FCC::Diamond001()
{
  char str[MAXLINE];
  int surftype = 1;
  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Please selection the type of Diamond (001):\n");
  printf("   1. Conventional cubic;\n");
  printf("   2. Orthogonal, x//[110], y//[-110];\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = atoi(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d", surftype);
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");

  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }
  name = memory->create(name,13,"Dia001:name");
  strcpy(name, "Diamond(001)");

  switch (surftype){
  case 1:
    nucell = 8;
    ntype  = 1;
    
    latvec[0][0] = 1.;
    latvec[1][1] = 1.;
    latvec[2][2] = 1.;
  
    atpos = memory->create(atpos,nucell,3,"Dia001:atpos");
    attyp = memory->create(attyp,nucell,"Dia001:atpos");
    
    for (int i=0; i<nucell; i++) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
  
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
  
    atpos[2][0] = 0.25;
    atpos[2][1] = 0.25;
    atpos[2][2] = 0.25;
  
    atpos[3][0] = 0.75;
    atpos[3][1] = 0.75;
    atpos[3][2] = 0.25;
  
    atpos[4][0] = 0.5;
    atpos[4][1] = 0.0;
    atpos[4][2] = 0.5;
  
    atpos[5][0] = 0.0;
    atpos[5][1] = 0.5;
    atpos[5][2] = 0.5;
  
    atpos[6][0] = 0.75;
    atpos[6][1] = 0.25;
    atpos[6][2] = 0.75;
  
    atpos[7][0] = 0.25;
    atpos[7][1] = 0.75;
    atpos[7][2] = 0.75;
  
    initialized = 1;
    break;
  case 2:
    nucell = 4;
    ntype  = 1;
  
    latvec[0][0] = latvec[1][1] = sqrt(0.5);
    latvec[2][2] = 1.;
    
    atpos = memory->create(atpos,nucell,3,"atpos");
    attyp = memory->create(attyp,nucell,"attyp");
    
    attyp[0] =  1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] =  1;
    atpos[1][0] = 0.50;
    atpos[1][1] = 0.50;
    atpos[1][2] = 0.50;
    
    attyp[2] =  1;
    atpos[2][0] = 0.50;
    atpos[2][1] = 0.00;
    atpos[2][2] = 0.75;
    
    attyp[3] =  1;
    atpos[3][0] = 0.00;
    atpos[3][1] = 0.50;
    atpos[3][2] = 0.25;
    
    initialized = 1;
    break;
  default:
    break;
  }

return;
}

/* ----------------------------------------------------------------------
   Initialize Diamond (110)
------------------------------------------------------------------------- */
void FCC::Diamond110()
{
  char str[MAXLINE];
  int surftype = 1;
  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Please selection the type of Diamond (110):\n");
  printf("   1. x//[001], y//[110];\n");
  printf("   2. x//[110], y//[001];\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = atoi(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d", surftype);
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");

  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }
  name = memory->create(name,13,"Dia:name");
  strcpy(name, "Diamond(110)");

  switch (surftype){
  case 1:
    nucell = 4;
    ntype  = 1;
    
    latvec[0][0] = 1.;
    latvec[1][1] = latvec[2][2] = sqrt(0.5);
    
    atpos = memory->create(atpos,nucell,3,"atpos");
    attyp = memory->create(attyp,nucell,"attyp");
    
    attyp[0] =  1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] =  1;
    atpos[1][0] = 0.75;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
    
    attyp[2] =  1;
    atpos[2][0] = 0.5;
    atpos[2][1] = 0.5;
    atpos[2][2] = 0.5;
    
    attyp[3] =  1;
    atpos[3][0] = 0.25;
    atpos[3][1] = 0.;
    atpos[3][2] = 0.5;
    
    initialized = 1;
    break;

  case 2:
    nucell = 4;
    ntype  = 1;
    
    latvec[1][1] = 1.;
    latvec[0][0] = latvec[2][2] = sqrt(0.5);
    
    atpos = memory->create(atpos,nucell,3,"atpos");
    attyp = memory->create(attyp,nucell,"attyp");
    
    attyp[0] =  1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] =  1;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.75;
    atpos[1][2] = 0.;
    
    attyp[2] =  1;
    atpos[2][0] = 0.5;
    atpos[2][1] = 0.5;
    atpos[2][2] = 0.5;
    
    attyp[3] =  1;
    atpos[3][0] = 0.;
    atpos[3][1] = 0.25;
    atpos[3][2] = 0.5;
    
    initialized = 1;
    break;
  default:
    break;
  }

return;
}

/* ----------------------------------------------------------------------
 * Initialize Diamond (111)
 *----------------------------------------------------------------------- */
void FCC::Diamond111()
{
  char str[MAXLINE];
  int surftype = 1;
  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Please selection the type of Diamond (111):\n");
  printf("   1. x//[110], gamma = 120 deg;\n");
  printf("   2. y//[110], gamma = 120 deg;\n");
  printf("   3. x//[110], gamma = 60 deg;\n");
  printf("   4. y//[110], gamma = 60 deg;\n");
  printf("   5. Orthogonal, long along x;\n");
  printf("   6. Orthogonal, long along y;\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = atoi(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d", surftype);
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");

  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }
  name = memory->create(name,13,"Dia:name");
  strcpy(name, "Diamond(111)");

  // some constants
  const double one12 = 1./12., fiv12 = 5./12., sev12 = 7./12., ele12 = 11./12.;
  const double one6 = 1./6., fiv6 = 5./6.;
  const double one3 = 1./3., two3 = 2./3.;

  switch (surftype){
  case 1:
    nucell =  6;
    ntype  =  1;
    
    latvec[0][0] = sqrt(0.5);
    
    latvec[1][0] = -sqrt(0.125);
    latvec[1][1] =  sqrt(0.375);
    latvec[2][2] =  sqrt(3.);
    
    atpos = memory->create(atpos,nucell,3,"atpos");
    attyp = memory->create(attyp,nucell,"attyp");
    
    attyp[0] =  1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] =  1;
    atpos[1][0] = two3;
    atpos[1][1] = one3;
    atpos[1][2] = two3;
    
    attyp[2] =  1;
    atpos[2][0] = one3;
    atpos[2][1] = two3;
    atpos[2][2] = one3;
    
    attyp[3] =  1;
    atpos[3][0] = 0.;
    atpos[3][1] = 0.;
    atpos[3][2] = 0.25;
    
    attyp[4] =  1;
    atpos[4][0] = two3;
    atpos[4][1] = one3;
    atpos[4][2] = ele12;
    
    attyp[5] =  1;
    atpos[5][0] = one3;
    atpos[5][1] = two3;
    atpos[5][2] = sev12;
    
    initialized = 1;
    break;

  case 2:
    nucell =  6;
    ntype  =  1;
    
    latvec[0][0] =  sqrt(0.375);
    latvec[0][1] = -sqrt(0.125);
    latvec[1][1] =  sqrt(0.5);
    latvec[2][2] =  sqrt(3.);
    
    atpos = memory->create(atpos,nucell,3,"atpos");
    attyp = memory->create(attyp,nucell,"attyp");
    
    attyp[0] =  1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] =  1;
    atpos[1][0] = two3;
    atpos[1][1] = one3;
    atpos[1][2] = two3;
    
    attyp[2] =  1;
    atpos[2][0] = one3;
    atpos[2][1] = two3;
    atpos[2][2] = one3;
    
    attyp[3] =  1;
    atpos[3][0] = 0.;
    atpos[3][1] = 0.;
    atpos[3][2] = 0.25;
    
    attyp[4] =  1;
    atpos[4][0] = two3;
    atpos[4][1] = one3;
    atpos[4][2] = ele12;
    
    attyp[5] =  1;
    atpos[5][0] = one3;
    atpos[5][1] = two3;
    atpos[5][2] = sev12;
    
    initialized = 1;
    break;
  case 3:
    nucell =    6;
    ntype  =  1;
    
    latvec[0][0] = sqrt(0.5);
    latvec[1][0] = sqrt(0.125);
    latvec[1][1] = sqrt(0.375);
    latvec[2][2] = sqrt(3.);
    
    atpos = memory->create(atpos,nucell,3,"atpos");
    attyp = memory->create(attyp,nucell,"attyp");
    
    attyp[0] =  1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] =  1;
    atpos[1][0] = two3;
    atpos[1][1] = two3;
    atpos[1][2] = one3;
    
    attyp[2] =  1;
    atpos[2][0] = one3;
    atpos[2][1] = one3;
    atpos[2][2] = two3;
    
    attyp[3] =  1;
    atpos[3][0] = 0.;
    atpos[3][1] = 0.;
    atpos[3][2] = 0.25;
    
    attyp[4] =  1;
    atpos[4][0] = two3;
    atpos[4][1] = two3;
    atpos[4][2] = sev12;
    
    attyp[5] =  1;
    atpos[5][0] = one3;
    atpos[5][1] = one3;
    atpos[5][2] = ele12;
    
    initialized = 1;
    break;

  case 4:
    nucell =    6;
    ntype  =  1;
    
    latvec[0][0] = sqrt(0.375);
    latvec[0][1] = sqrt(0.125);
    latvec[1][1] = sqrt(0.5);
    latvec[2][2] = sqrt(3.);
    
    atpos = memory->create(atpos,nucell,3,"atpos");
    attyp = memory->create(attyp,nucell,"attyp");
    
    attyp[0] =  1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] =  1;
    atpos[1][0] = two3;
    atpos[1][1] = two3;
    atpos[1][2] = one3;
    
    attyp[2] =  1;
    atpos[2][0] = one3;
    atpos[2][1] = one3;
    atpos[2][2] = two3;
    
    attyp[3] =  1;
    atpos[3][0] = 0.;
    atpos[3][1] = 0.;
    atpos[3][2] = 0.25;
    
    attyp[4] =  1;
    atpos[4][0] = two3;
    atpos[4][1] = two3;
    atpos[4][2] = sev12;
    
    attyp[5] =  1;
    atpos[5][0] = one3;
    atpos[5][1] = one3;
    atpos[5][2] = ele12;
    
    initialized = 1;
    break;

  case 6:
    nucell = 12;
    ntype  = 1;
    
    latvec[0][0] = sqrt(0.5);
    latvec[1][1] = sqrt(1.5);
    latvec[2][2] = sqrt(3.);
    
    atpos = memory->create(atpos,nucell,3,"atpos");
    attyp = memory->create(attyp,nucell,"attyp");
    
    attyp[ 0] =  1;
    atpos[ 0][0] = 0.;
    atpos[ 0][1] = 0.;
    atpos[ 0][2] = 0.;
    
    attyp[ 1] =  1;
    atpos[ 1][0] = 0.;
    atpos[ 1][1] = two3;
    atpos[ 1][2] = one3;
    
    attyp[ 2] =  1;
    atpos[ 2][0] = 0.5;
    atpos[ 2][1] = 0.5;
    atpos[ 2][2] = 0.;
    
    attyp[ 3] =  1;
    atpos[ 3][0] = 0.;
    atpos[ 3][1] = one3;
    atpos[ 3][2] = two3;
    
    attyp[ 4] =  1;
    atpos[ 4][0] = 0.5;
    atpos[ 4][1] = one6;
    atpos[ 4][2] = one3;
    
    attyp[ 5] =  1;
    atpos[ 5][0] = 0.5;
    atpos[ 5][1] = fiv6;
    atpos[ 5][2] = two3;
    
    attyp[ 6] =  1;
    atpos[ 6][0] = 0.;
    atpos[ 6][1] = 0.;
    atpos[ 6][2] = 0.25;
    
    attyp[ 7] =  1;
    atpos[ 7][0] = 0.;
    atpos[ 7][1] = two3;
    atpos[ 7][2] = sev12;
    
    attyp[ 8] =  1;
    atpos[ 8][0] = 0.;
    atpos[ 8][1] = one3;
    atpos[ 8][2] = ele12;
    
    attyp[ 9] =  1;
    atpos[ 9][0] = 0.5;
    atpos[ 9][1] = 0.5;
    atpos[ 9][2] = 0.25;
    
    attyp[10] =  1;
    atpos[10][0] = 0.5;
    atpos[10][1] = one6;
    atpos[10][2] = sev12;
    
    attyp[11] =  1;
    atpos[11][0] = 0.5;
    atpos[11][1] = fiv6;
    atpos[11][2] = ele12;
    
    initialized = 1;
    break;

  case 5:
    nucell = 12;
    ntype  = 1;
    
    latvec[0][0] = sqrt(1.5);
    latvec[1][1] = sqrt(0.5);
    latvec[2][2] = sqrt(3.);
    
    atpos = memory->create(atpos,nucell,3,"atpos");
    attyp = memory->create(attyp,nucell,"attyp");
    
    attyp[ 0] =  1;
    atpos[ 0][0] = 0.;
    atpos[ 0][1] = 0.;
    atpos[ 0][2] = 0.;
    
    attyp[ 1] =  1;
    atpos[ 1][1] = 0.;
    atpos[ 1][0] = two3;
    atpos[ 1][2] = one3;
    
    attyp[ 2] =  1;
    atpos[ 2][1] = 0.5;
    atpos[ 2][0] = 0.5;
    atpos[ 2][2] = 0.;
    
    attyp[ 3] =  1;
    atpos[ 3][1] = 0.;
    atpos[ 3][0] = one3;
    atpos[ 3][2] = two3;
    
    attyp[ 4] =  1;
    atpos[ 4][1] = 0.5;
    atpos[ 4][0] = one6;
    atpos[ 4][2] = one3;
    
    attyp[ 5] =  1;
    atpos[ 5][1] = 0.5;
    atpos[ 5][0] = fiv6;
    atpos[ 5][2] = two3;
    
    attyp[ 6] =  1;
    atpos[ 6][1] = 0.;
    atpos[ 6][0] = 0.;
    atpos[ 6][2] = 0.25;
    
    attyp[ 7] =  1;
    atpos[ 7][1] = 0.;
    atpos[ 7][0] = two3;
    atpos[ 7][2] = sev12;
    
    attyp[ 8] =  1;
    atpos[ 8][1] = 0.;
    atpos[ 8][0] = one3;
    atpos[ 8][2] = ele12;
    
    attyp[ 9] =  1;
    atpos[ 9][1] = 0.5;
    atpos[ 9][0] = 0.5;
    atpos[ 9][2] = 0.25;
    
    attyp[10] =  1;
    atpos[10][1] = 0.5;
    atpos[10][0] = one6;
    atpos[10][2] = sev12;
    
    attyp[11] =  1;
    atpos[11][1] = 0.5;
    atpos[11][0] = fiv6;
    atpos[11][2] = ele12;
    
    initialized = 1;
    break;
  default:
    break;
  }

return;
}
/* ---------------------------------------------------------------------- */
