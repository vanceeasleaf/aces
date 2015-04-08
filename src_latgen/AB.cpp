#include "AB.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

#define MAXLINE 256

using namespace std;

/* ----------------------------------------------------------------------
   To select the orientation of the lattice
------------------------------------------------------------------------- */
AB::AB() : lattice()
{
  char str[MAXLINE];
  alat = 1.; ca = 1.;
  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  int lattype = 1;
  printf("Please select the type of your lattice:\n");
  printf("   1. B1 (NaCl);        4. L10 (CuAu);\n");
  printf("   2. B2 (CsCl);        5. B81 (a-NiAs);\n");
  printf("   3. B3 (Zincblende);  6. B4 (Wurtzite);\n");
  for (int i=0; i<14; i++) printf("-----"); printf("\n");
  printf("   7. Perovskite;\n");
  for (int i=0; i<14; i++) printf("-----"); printf("\n");
  printf("Your choice [%d]: ", lattype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) lattype = atoi(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d\n", lattype);

  printf("Please input the lattice constant of your AB crystal [%g]: ", alat);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) alat = atof(strtok(str, " \t\n\r\f"));
  if (lattype == 4 || lattype == 5 || lattype == 6){
    if (lattype == 6) ca = sqrt(8./3.);
    printf("Please input the c/a or c (negative) of your crystal [%g]: ", ca);
    if (count_words(fgets(str,MAXLINE,stdin)) > 0) ca = atof(strtok(str, " \t\n\r\f"));
    if (ca < 0.) ca = -ca/alat;
  }
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  
  // initialize according to orientation
  initialized = 0;
  switch (lattype){
  case 1:
    AB_B1();
    break;
  case 2:
    AB_B2();
    break;
  case 3:
    AB_B3();
    break;
  case 4:
    AB_L10();
    break;
  case 5:
    AB_NiAs();
    break;
  case 6:
    AB_B4();
    break;
  case 7:
    AB_Perov();
    break;
  default:
    break;
  }

}

/* ----------------------------------------------------------------------
   Deconstructor does nothing
------------------------------------------------------------------------- */
AB::~AB()
{

}

/* ----------------------------------------------------------------------
   Initialize for B1 lattice
------------------------------------------------------------------------- */
void AB::AB_B1()
{
  char str[MAXLINE];
  int surftype = 1;
  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Please selection the type of AB-B1 (NaCl) surface:\n");
  printf("   1. (001), small;\n");
  printf("   2. (001), conventional;\n");
  printf("   3. (110), long along y;\n");
  printf("   4. (111), long along x, orthogonal;\n");
  printf("   5. primitive cell;\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = atoi(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d", surftype);
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }

  // initialize according to surface type
  switch (surftype){
  case 1:
    name = memory->create(name,11,"AB:name");
    strcpy(name, "AB-B1(001)");

    ntype  = 2;
    nucell = 4;
    
    latvec[0][0] = 1./sqrt(2.);
    latvec[1][1] = 1./sqrt(2.);
    latvec[2][2] = 1.;

    atpos = memory->create(atpos, nucell, 3, "AB:atpos");
    attyp = memory->create(attyp, nucell, "AB:attyp");
    
    attyp[0] = 2;
    atpos[0][0] = 0.5;
    atpos[0][1] = 0.5;
    atpos[0][2] = 0.;
    
    attyp[1] = 1;
    atpos[1][0] = 0.;
    atpos[1][1] = 0.;
    atpos[1][2] = 0.;
    
    attyp[2] = 2;
    atpos[2][0] = 0.;
    atpos[2][1] = 0.;
    atpos[2][2] = 0.5;
    
    attyp[3] = 1;
    atpos[3][0] = 0.5;
    atpos[3][1] = 0.5;
    atpos[3][2] = 0.5;

    initialized = 1;
    break;
  case 2:
    name = memory->create(name,11,"AB:name");
    strcpy(name, "AB-B1(001)");

    ntype  = 2;
    nucell = 8;
    
    latvec[0][0] = 1.;
    latvec[1][1] = 1.;
    latvec[2][2] = 1.;

    atpos = memory->create(atpos, nucell, 3, "AB:atpos");
    attyp = memory->create(attyp, nucell, "AB:attyp");

    attyp[0] = 2;
    atpos[0][0] = 0.5;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = 2;
    atpos[1][0] = 0.;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
    
    attyp[2] = 1;
    atpos[2][0] = 0.;
    atpos[2][1] = 0.;
    atpos[2][2] = 0.;
    
    attyp[3] = 1;
    atpos[3][0] = 0.5;
    atpos[3][1] = 0.5;
    atpos[3][2] = 0.;
    
    attyp[4] = 2;
    atpos[4][0] = 0.5;
    atpos[4][1] = 0.5;
    atpos[4][2] = 0.5;
    
    attyp[5] = 2;
    atpos[5][0] = 0.;
    atpos[5][1] = 0.;
    atpos[5][2] = 0.5;
    
    attyp[6] = 1;
    atpos[6][0] = 0.;
    atpos[6][1] = 0.5;
    atpos[6][2] = 0.5;
    
    attyp[7] = 1;
    atpos[7][0] = 0.5;
    atpos[7][1] = 0.;
    atpos[7][2] = 0.5;
    
    initialized = 1;
    break;
  case 3:
    name = memory->create(name,11,"AB:name");
    strcpy(name, "AB-B1(110)");

    ntype  = 2;
    nucell = 4;
    
    latvec[0][0] = 1.;
    latvec[1][1] = 1./sqrt(2.);
    latvec[2][2] = 1./sqrt(2.);

    atpos = memory->create(atpos, nucell, 3, "AB:atpos");
    attyp = memory->create(attyp, nucell, "AB:attyp");
    
    attyp[0] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = 2;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.;
    atpos[1][2] = 0.;
    
    attyp[2] = 2;
    atpos[2][0] = 0;
    atpos[2][1] = 0.5;
    atpos[2][2] = 0.5;
    
    attyp[3] = 1;
    atpos[3][0] = 0.5;
    atpos[3][1] = 0.5;
    atpos[3][2] = 0.5;

    initialized = 1;
    break;
  case 4:
    name = memory->create(name,11,"AB:name");
    strcpy(name, "AB-B1(111)");

    ntype  = 2;
    nucell = 12;
    
    latvec[0][0] = sqrt(6.)/2.;
    latvec[1][1] = 1./sqrt(2.);
    latvec[2][2] = sqrt(3.);

    atpos = memory->create(atpos, nucell, 3, "AB:atpos");
    attyp = memory->create(attyp, nucell, "AB:attyp");
    
    attyp[0] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = 1;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
    
    attyp[2] = 2;
    atpos[2][0] = 1./3.;
    atpos[2][1] = 0.;
    atpos[2][2] = 1./6.;
    
    attyp[3] = 2;
    atpos[3][0] = 5./6.;
    atpos[3][1] = 0.5;
    atpos[3][2] = 1./6.;
    
    attyp[4] = 1;
    atpos[4][0] = 2./3.;
    atpos[4][1] = 0.;
    atpos[4][2] = 1./3.;
    
    attyp[5] = 1;
    atpos[5][0] = 1./6.;
    atpos[5][1] = 0.5;
    atpos[5][2] = 1./3.;
    
    attyp[6] = 2;
    atpos[6][0] = 0.;
    atpos[6][1] = 0.;
    atpos[6][2] = 0.5;
    
    attyp[7] = 2;
    atpos[7][0] = 0.5;
    atpos[7][1] = 0.5;
    atpos[7][2] = 0.5;
    
    attyp[8] = 1;
    atpos[8][0] = 1./3.;
    atpos[8][1] = 0.;
    atpos[8][2] = 2./3.;
    
    attyp[9] = 1;
    atpos[9][0] = 5./6.;
    atpos[9][1] = 0.5;
    atpos[9][2] = 2./3.;
    
    attyp[10] = 2;
    atpos[10][0] = 2./3.;
    atpos[10][1] = 0.;
    atpos[10][2] = 5./6.;
    
    attyp[11] = 2;
    atpos[11][0] = 1./6.;
    atpos[11][1] = 0.5;
    atpos[11][2] = 5./6.;

    initialized = 1;
    break;
  case 5:
    name = memory->create(name,16,"AB:name");
    strcpy(name, "AB-B1-primitive");

    ntype  = 2;
    nucell = 2;
    
    latvec[0][1] = 0.5;
    latvec[0][2] = 0.5;
    latvec[1][0] = 0.5;
    latvec[1][2] = 0.5;
    latvec[2][0] = 0.5;
    latvec[2][1] = 0.5;

    atpos = memory->create(atpos, nucell, 3, "AB:atpos");
    attyp = memory->create(attyp, nucell, "AB:attyp");
    
    attyp[0] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = 2;
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
   Initialize for B2
------------------------------------------------------------------------- */
void AB::AB_B2()
{
  char str[MAXLINE];
  int surftype = 1;
  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Please selection the type of AB-B2 (CsCl) surface:\n");
  printf("   1. (100);\n");
  printf("   2. (110);\n");
  printf("   3. (111);\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = atoi(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d", surftype);
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }

  // initialize according to surface type
  switch (surftype){
  case 1:
    name = memory->create(name,11,"AB:name");
    strcpy(name, "AB-B2(001)");

    ntype  = 2;
    nucell = 2;
    
    latvec[0][0] = 1.;
    latvec[1][1] = 1.;
    latvec[2][2] = 1.;

    atpos = memory->create(atpos, nucell, 3, "AB:atpos");
    attyp = memory->create(attyp, nucell, "AB:attyp");
    
    attyp[0]    = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    attyp[1]    = 2;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.5;

    initialized = 1;
    break;
  case 2:
    name = memory->create(name,11,"AB:name");
    strcpy(name, "AB-B2(110)");

    ntype  = 2;
    nucell = 4;
    
    latvec[0][0] = 1.;
    latvec[1][1] = sqrt(2.);
    latvec[2][2] = sqrt(2.);

    atpos = memory->create(atpos, nucell, 3, "AB:atpos");
    attyp = memory->create(attyp, nucell, "AB:attyp");
    
    attyp[0] = 2;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = 1;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
    
    attyp[2] = 2;
    atpos[2][0] = 0.;
    atpos[2][1] = 0.5;
    atpos[2][2] = 0.5;
    
    attyp[3] = 1;
    atpos[3][0] = 0.5;
    atpos[3][1] = 0.;
    atpos[3][2] = 0.5;

    initialized = 1;
    break;
  case 3:
    name = memory->create(name,11,"AB:name");
    strcpy(name, "AB-B2(111)");

    ntype  = 2;
    nucell = 12;
    
    latvec[0][0] = sqrt(2.);
    latvec[1][1] = sqrt(6.);
    latvec[2][2] = sqrt(3.);

    atpos = memory->create(atpos, nucell, 3, "AB:atpos");
    attyp = memory->create(attyp, nucell, "AB:attyp");
    
    attyp[0] = 2;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = 2;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
    
    attyp[2] = 1;
    atpos[2][0] = 0.;
    atpos[2][1] = 2./3.;
    atpos[2][2] = 1./6.;
    
    attyp[3] = 1;
    atpos[3][0] = 0.5;
    atpos[3][1] = 1./6.;
    atpos[3][2] = 1./6.;
    
    attyp[4] = 2;
    atpos[4][0] = 0.;
    atpos[4][1] = 1./3.;
    atpos[4][2] = 1./3.;
    
    attyp[5] = 2;
    atpos[5][0] = 0.5;
    atpos[5][1] = 5./6.;
    atpos[5][2] = 1./3.;
    
    attyp[6] = 1;
    atpos[6][0] = 0.;
    atpos[6][1] = 0.;
    atpos[6][2] = 0.5;
    
    attyp[7] = 1;
    atpos[7][0] = 0.5;
    atpos[7][1] = 0.5;
    atpos[7][2] = 0.5;
    
    attyp[8] = 2;
    atpos[8][0] = 0.;
    atpos[8][1] = 2./3.;
    atpos[8][2] = 2./3.;
    
    attyp[9] = 2;
    atpos[9][0] = 0.5;
    atpos[9][1] = 1./6.;
    atpos[9][2] = 2./3.;
    
    attyp[10] = 1;
    atpos[10][0] = 0.;
    atpos[10][1] = 1./3.;
    atpos[10][2] = 5./6.;
    
    attyp[11] = 1;
    atpos[11][0] = 0.5;
    atpos[11][1] = 5./6.;
    atpos[11][2] = 5./6.;

    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/* ----------------------------------------------------------------------
   Initialize for B3 lattice
------------------------------------------------------------------------- */
void AB::AB_B3()
{
  char str[MAXLINE];
  int surftype = 1;
  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Please selection the type of AB-B3 surface:\n");
  printf("   1. (100), small;\n");
  printf("   2. (100), conventional;\n");
  printf("   3. (110), long along y, orthogonal;\n");
  printf("   4. (111), long along y, orthogonal;\n");
  printf("   5. primitive;\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = atoi(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d", surftype);
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }

  // initialize according to surface type
  switch (surftype){
  case 1:
    name = memory->create(name,11,"AB:name");
    strcpy(name, "AB-B3(100)");

    ntype  = 2;
    nucell = 4;
    
    latvec[0][0] = 1./sqrt(2.);
    latvec[1][1] = 1./sqrt(2.);
    latvec[2][2] = 1.;

    atpos = memory->create(atpos, nucell, 3, "AB:atpos");
    attyp = memory->create(attyp, nucell, "AB:attyp");
    
    attyp[0] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = 2;
    atpos[1][0] = 0.;
    atpos[1][1] = 0.50;
    atpos[1][2] = 0.25;
    
    attyp[2] = 1;
    atpos[2][0] = 0.50;
    atpos[2][1] = 0.50;
    atpos[2][2] = 0.50;
    
    attyp[3] = 2;
    atpos[3][0] = 0.50;
    atpos[3][1] = 0.00;
    atpos[3][2] = 0.75;

    initialized = 1;
    break;
  case 2:
    name = memory->create(name,11,"AB:name");
    strcpy(name, "AB-B3(100)");

    ntype  = 2;
    nucell = 8;
    
    latvec[0][0] = 1.;
    latvec[1][1] = 1.;
    latvec[2][2] = 1.;

    atpos = memory->create(atpos, nucell, 3, "AB:atpos");
    attyp = memory->create(attyp, nucell, "AB:attyp");
    
    attyp[0] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = 1;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
    
    attyp[2] = 2;
    atpos[2][0] = 0.25;
    atpos[2][1] = 0.25;
    atpos[2][2] = 0.25;
    
    attyp[3] = 2;
    atpos[3][0] = 0.75;
    atpos[3][1] = 0.75;
    atpos[3][2] = 0.25;
    
    attyp[4] = 1;
    atpos[4][0] = 0.;
    atpos[4][1] = 0.5;
    atpos[4][2] = 0.5;
    
    attyp[5] = 1;
    atpos[5][0] = 0.5;
    atpos[5][1] = 0.;
    atpos[5][2] = 0.5;
    
    attyp[6] = 2;
    atpos[6][0] = 0.75;
    atpos[6][1] = 0.25;
    atpos[6][2] = 0.75;
    
    attyp[7] = 2;
    atpos[7][0] = 0.25;
    atpos[7][1] = 0.75;
    atpos[7][2] = 0.75;

    initialized = 1;
    break;
  case 3:
    name = memory->create(name,11,"AB:name");
    strcpy(name, "AB-B3(110)");

    ntype  = 2;
    nucell = 4;
    
    latvec[0][0] = 1.;
    latvec[1][1] = 1./sqrt(2.);
    latvec[2][2] = 1./sqrt(2.);

    atpos = memory->create(atpos, nucell, 3, "AB:atpos");
    attyp = memory->create(attyp, nucell, "AB:attyp");
    
    attyp[0] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    attyp[1] = 2;
    atpos[1][0] = 0.75;
    atpos[1][1] = 0.50;
    atpos[1][2] = 0.;
    
    attyp[2] = 2;
    atpos[2][0] = 0.25;
    atpos[2][1] = 0.;
    atpos[2][2] = 0.50;
    
    attyp[3] = 1;
    atpos[3][0] = 0.5;
    atpos[3][1] = 0.5;
    atpos[3][2] = 0.5;
    
    initialized = 1;
    break;
  case 4:
    name = memory->create(name,11,"AB:name");
    strcpy(name, "AB-B3(111)");

    ntype  = 2;
    nucell = 12;
    
    latvec[0][0] = sqrt(2.)*0.5;
    latvec[1][1] = sqrt(6.)*0.5;
    latvec[2][2] = sqrt(3.);

    atpos = memory->create(atpos, nucell, 3, "AB:atpos");
    attyp = memory->create(attyp, nucell, "AB:attyp");
    
    attyp[0] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = 1;
    atpos[1][0] = 0.50;
    atpos[1][1] = 0.50;
    atpos[1][2] = 0.;
    
    attyp[2] = 2;
    atpos[2][0] = 0.;
    atpos[2][1] = 0.;
    atpos[2][2] = 0.25;
    
    attyp[3] = 2;
    atpos[3][0] = 0.5;
    atpos[3][1] = 0.5;
    atpos[3][2] = 0.25;
    
    attyp[4] = 1;
    atpos[4][0] = 0.;
    atpos[4][1] = 2./3.;
    atpos[4][2] = 1./3.;
    
    attyp[5] = 1;
    atpos[5][0] = 0.5;
    atpos[5][1] = 1./6.;
    atpos[5][2] = 1./3.;
    
    attyp[6] = 2;
    atpos[6][0] = 0.;
    atpos[6][1] = 2./3.;
    atpos[6][2] = 7./12.;
    
    attyp[7] = 2;
    atpos[7][0] = 0.5;
    atpos[7][1] = 1./6.;
    atpos[7][2] = 7./12.;
    
    attyp[8] = 1;
    atpos[8][0] = 0.;
    atpos[8][1] = 1./3.;
    atpos[8][2] = 2./3.;
    
    attyp[9] = 1;
    atpos[9][0] = 0.5;
    atpos[9][1] = 5./6.;
    atpos[9][2] = 2./3.;
    
    attyp[10] = 2;
    atpos[10][0] = 0.5;
    atpos[10][1] = 5./6.;
    atpos[10][2] = 11./12.;
    
    attyp[11] = 2;
    atpos[11][0] = 0.;
    atpos[11][1] = 1./3.;
    atpos[11][2] = 11./12.;

    initialized = 1;
    break;
  case 5:
    name = memory->create(name,16,"AB:name");
    strcpy(name, "AB-B3-primitive");

    ntype  = 2;
    nucell = 2;
    
    latvec[0][1] = 0.5;
    latvec[0][2] = 0.5;
    latvec[1][0] = 0.5;
    latvec[1][2] = 0.5;
    latvec[2][0] = 0.5;
    latvec[2][1] = 0.5;

    atpos = memory->create(atpos, nucell, 3, "AB:atpos");
    attyp = memory->create(attyp, nucell, "AB:attyp");
    
    attyp[0] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    attyp[1] = 2;
    atpos[1][0] = 0.25;
    atpos[1][1] = 0.25;
    atpos[1][2] = 0.25;
    
    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/* ----------------------------------------------------------------------
   Initialize for L10 lattice
------------------------------------------------------------------------- */
void AB::AB_L10()
{
  char str[MAXLINE];
  int surftype = 1;
  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Please selection the type of AB-L10 surface:\n");
  printf("   1. (001);\n");
  printf("   2. (100);\n");
  printf("   3. (110), long along y, orthogonal;\n");
  printf("   4. primitive cell;\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = atoi(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d", surftype);
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }

  // initialize according to surface type
  switch (surftype){
  case 1:
    name = memory->create(name,12,"AB:name");
    strcpy(name, "AB-L10(001)");

    ntype  = 2;
    nucell = 4;
    
    latvec[0][0] = 1.;
    latvec[1][1] = 1.;
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "AB:atpos");
    attyp = memory->create(attyp, nucell, "AB:attyp");

    attyp[0] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = 1;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
    
    attyp[2] = 2;
    atpos[2][0] = 0.5;
    atpos[2][1] = 0.;
    atpos[2][2] = 0.5;
    
    attyp[3] = 2;
    atpos[3][0] = 0.;
    atpos[3][1] = 0.5;
    atpos[3][2] = 0.5;

    initialized = 1;
    break;
  case 2:
    name = memory->create(name,14,"AB:name");
    strcpy(name, "AB-L10(100)");

    ntype  = 2;
    nucell = 4;
    
    latvec[0][0] = 1.;
    latvec[1][1] = ca;
    latvec[2][2] = 1.;

    atpos = memory->create(atpos, nucell, 3, "AB:atpos");
    attyp = memory->create(attyp, nucell, "AB:attyp");
    
    attyp[0] = 2;
    atpos[0][0] = 0.5;
    atpos[0][1] = 0.5;
    atpos[0][2] = 0.;
    
    attyp[1] = 1;
    atpos[1][0] = 0.;
    atpos[1][1] = 0.;
    atpos[1][2] = 0.;
    
    attyp[2] = 2;
    atpos[2][0] = 0.;
    atpos[2][1] = 0.5;
    atpos[2][2] = 0.5;
    
    attyp[3] = 1;
    atpos[3][0] = 0.5;
    atpos[3][1] = 0.;
    atpos[3][2] = 0.5;

    initialized = 1;
    break;
  case 3:
    name = memory->create(name,12,"AB:name");
    strcpy(name, "AB-L10(110)");

    ntype  = 2;
    nucell = 4;
    
    latvec[0][0] = ca;
    latvec[1][1] = sqrt(2.);
    latvec[2][2] = sqrt(2.)/2.;

    atpos = memory->create(atpos, nucell, 3, "AB:atpos");
    attyp = memory->create(attyp, nucell, "AB:attyp");
    
    attyp[0] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = 1;
    atpos[1][0] = 0.;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
    
    attyp[2] = 2;
    atpos[2][0] = 0.5;
    atpos[2][1] = 0.75;
    atpos[2][2] = 0.5;
    
    attyp[3] = 2;
    atpos[3][0] = 0.5;
    atpos[3][1] = 0.25;
    atpos[3][2] = 0.5;

    initialized = 1;
    break;
  case 4:
    name = memory->create(name,17,"AB:name");
    strcpy(name, "AB-L10-primitive");

    ntype  = 2;
    nucell = 2;
    
    latvec[0][0] = 1./sqrt(2.);
    latvec[1][1] = 1./sqrt(2.);
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "AB:atpos");
    attyp = memory->create(attyp, nucell, "AB:attyp");
    
    attyp[0] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = 2;
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
   Initialize for a-NiAs lattice
------------------------------------------------------------------------- */
void AB::AB_NiAs()
{
  char str[MAXLINE];
  int surftype = 1;
  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Please selection the type of AB-a-NiAs surface:\n");
  printf("   1. (001), conventional;\n");
  printf("   2. (001), orthogonal;\n");
  printf("   3. (100);\n");
  printf("   4. (110), orthogonal, long along y;\n");
  printf("   5. (1-10), orthogonal, long along y;\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = atoi(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d", surftype);
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }

  // initialize according to surface type
  switch (surftype){
  case 1:
    name = memory->create(name,13,"AB:name");
    strcpy(name, "AB-NiAs(001)");

    ntype  = 2;
    nucell = 4;
    
    latvec[0][0] = 1.;
    latvec[1][0] = -0.5;
    latvec[1][1] = sqrt(3.)/2.;
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "AB:atpos");
    attyp = memory->create(attyp, nucell, "AB:attyp");

    attyp[0] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = 2;
    atpos[1][0] = 1./3.;
    atpos[1][1] = 2./3.;
    atpos[1][2] = 0.25;
    
    attyp[2] = 1;
    atpos[2][0] = 0.;
    atpos[2][1] = 0.;
    atpos[2][2] = 0.5;
    
    attyp[3] = 2;
    atpos[3][0] = 2./3.;
    atpos[3][1] = 1./3.;
    atpos[3][2] = 0.75;

    initialized = 1;
    break;
  case 2:
    name = memory->create(name,13,"AB:name");
    strcpy(name, "AB-NiAs(001)");

    ntype  = 2;
    nucell = 8;
    
    latvec[0][0] = 1.;
    latvec[1][1] = sqrt(3.);
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "AB:atpos");
    attyp = memory->create(attyp, nucell, "AB:attyp");

    attyp[0] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = 1;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
    
    attyp[2] = 2;
    atpos[2][0] = 0.;
    atpos[2][1] = 2./3.;
    atpos[2][2] = 0.25;
    
    attyp[3] = 2;
    atpos[3][0] = 0.5;
    atpos[3][1] = 1./6.;
    atpos[3][2] = 0.25;
    
    attyp[4] = 1;
    atpos[4][0] = 0.;
    atpos[4][1] = 0.;
    atpos[4][2] = 0.5;
    
    attyp[5] = 1;
    atpos[5][0] = 0.5;
    atpos[5][1] = 0.5;
    atpos[5][2] = 0.5;
    
    attyp[6] = 2;
    atpos[6][0] = 0.;
    atpos[6][1] = 1./3.;
    atpos[6][2] = 0.75;
    
    attyp[7] = 2;
    atpos[7][0] = 0.5;
    atpos[7][1] = 5./6.;
    atpos[7][2] = 0.75;
    
    initialized = 1;
    break;
  case 3:
    name = memory->create(name,13,"AB:name");
    strcpy(name, "AB-NiAs(100)");

    ntype  = 2;
    nucell = 8;
    
    latvec[0][0] = 1.;
    latvec[1][1] = ca;
    latvec[2][2] = sqrt(3.);

    atpos = memory->create(atpos, nucell, 3, "AB:atpos");
    attyp = memory->create(attyp, nucell, "AB:attyp");

    attyp[0] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = 1;
    atpos[1][0] = 0.;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
    
    attyp[2] = 2;
    atpos[2][0] = 0.5;
    atpos[2][1] = 0.25;
    atpos[2][2] = 1./6.;
    
    attyp[3] = 2;
    atpos[3][0] = 0.;
    atpos[3][1] = 0.75;
    atpos[3][2] = 1./3.;
    
    attyp[4] = 1;
    atpos[4][0] = 0.5;
    atpos[4][1] = 0.;
    atpos[4][2] = 0.5;
    
    attyp[5] = 1;
    atpos[5][0] = 0.5;
    atpos[5][1] = 0.5;
    atpos[5][2] = 0.5;
    
    attyp[6] = 2;
    atpos[6][0] = 0.;
    atpos[6][1] = 0.25;
    atpos[6][2] = 2./3.;
    
    attyp[7] = 2;
    atpos[7][0] = 0.5;
    atpos[7][1] = 0.75;
    atpos[7][2] = 5./6.;

    initialized = 1;
    break;
  case 4:
    name = memory->create(name,13,"AB:name");
    strcpy(name, "AB-NiAs(110)");

    ntype  = 2;
    nucell = 8;
    
    latvec[0][0] = ca;
    latvec[1][1] = sqrt(3.);
    latvec[2][2] = 1.;

    atpos = memory->create(atpos, nucell, 3, "AB:atpos");
    attyp = memory->create(attyp, nucell, "AB:attyp");

    attyp[0] = 2;
    atpos[0][0] = 0.25;
    atpos[0][1] = 1./3.;
    atpos[0][2] = 0.;
    
    attyp[1] = 2;
    atpos[1][0] = 0.75;
    atpos[1][1] = 2./3.;
    atpos[1][2] = 0.;
    
    attyp[2] = 1;
    atpos[2][0] = 0.;
    atpos[2][1] = 0.;
    atpos[2][2] = 0.;
    
    attyp[3] = 1;
    atpos[3][0] = 0.5;
    atpos[3][1] = 0.;
    atpos[3][2] = 0.;
    
    attyp[4] = 2;
    atpos[4][0] = 0.25;
    atpos[4][1] = 5./6.;
    atpos[4][2] = 0.5;
    
    attyp[5] = 2;
    atpos[5][0] = 0.75;
    atpos[5][1] = 1./6.;
    atpos[5][2] = 0.5;
    
    attyp[6] = 1;
    atpos[6][0] = 0.;
    atpos[6][1] = 0.5;
    atpos[6][2] = 0.5;
    
    attyp[7] = 1;
    atpos[7][0] = 0.5;
    atpos[7][1] = 0.5;
    atpos[7][2] = 0.5;

    initialized = 1;
    break;
  case 5:
    name = memory->create(name,14,"AB:name");
    strcpy(name, "AB-NiAs(1-10)");

    ntype  = 2;
    nucell = 8;
    
    latvec[0][0] = 1.;
    latvec[1][1] = ca;
    latvec[2][2] = sqrt(3.);

    atpos = memory->create(atpos, nucell, 3, "AB:atpos");
    attyp = memory->create(attyp, nucell, "AB:attyp");

    attyp[0] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = 1;
    atpos[1][0] = 0.;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
    
    attyp[2] = 2;
    atpos[2][0] = 0.5;
    atpos[2][1] = 0.75;
    atpos[2][2] = 1./6.;
    
    attyp[3] = 2;
    atpos[3][0] = 0.;
    atpos[3][1] = 0.25;
    atpos[3][2] = 1./3.;
    
    attyp[4] = 1;
    atpos[4][0] = 0.5;
    atpos[4][1] = 0.;
    atpos[4][2] = 0.5;
    
    attyp[5] = 1;
    atpos[5][0] = 0.5;
    atpos[5][1] = 0.5;
    atpos[5][2] = 0.5;
    
    attyp[6] = 2;
    atpos[6][0] = 0.;
    atpos[6][1] = 0.75;
    atpos[6][2] = 2./3.;
    
    attyp[7] = 2;
    atpos[7][0] = 0.5;
    atpos[7][1] = 0.25;
    atpos[7][2] = 5./6.;

    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/* ----------------------------------------------------------------------
   Initialize for B4 (Wurtzite) lattice
------------------------------------------------------------------------- */
void AB::AB_B4()
{
  char str[MAXLINE];
  int surftype = 1;
  double u = 0.375;
  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Please indicate the relative position of the two lattices along c [%g]: ", u);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) u = atof(strtok(str, " \t\n\r\f"));
  printf("\nPlease selection the type of AB-B4 surface:\n");
  printf("   1. (001), gamma =  60 degree;\n");
  printf("   2. (001), gamma = 120 degree;\n");
  printf("   3. (001), orthogonal, long x;\n");
  printf("   4. (001), orthogonal, long y;\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = atoi(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d", surftype);
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }

  // initialize according to surface type
  switch (surftype){
  case 1:
    name = memory->create(name,11,"AB:name");
    strcpy(name, "AB-B4(001)");

    ntype  = 2;
    nucell = 4;
    
    latvec[0][0] = 1.;
    latvec[1][0] = 0.5;
    latvec[1][1] = sqrt(0.75);
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "AB:atpos");
    attyp = memory->create(attyp, nucell, "AB:attyp");
    
    attyp[0]    = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    attyp[1]    = 1;
    atpos[1][0] = 1./3.;
    atpos[1][1] = 1./3.;
    atpos[1][2] = 0.5;

    attyp[2]    = 2;
    atpos[2][0] = 0.;
    atpos[2][1] = 0.;
    atpos[2][2] = u;

    attyp[3]    = 2;
    atpos[3][0] = 1./3.;
    atpos[3][1] = 1./3.;
    atpos[3][2] = 0.5+u;

    initialized = 1;
    break;
  case 2:
    name = memory->create(name,11,"AB:name");
    strcpy(name, "AB-B4(001)");

    ntype  = 2;
    nucell = 4;
    
    latvec[0][0] = 1.;
    latvec[1][0] =-0.5;
    latvec[1][1] = sqrt(0.75);
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "AB:atpos");
    attyp = memory->create(attyp, nucell, "AB:attyp");
    
    attyp[0]    = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    attyp[1]    = 1;
    atpos[1][0] = 2./3.;
    atpos[1][1] = 1./3.;
    atpos[1][2] = 0.5;

    attyp[2]    = 2;
    atpos[2][0] = 0.;
    atpos[2][1] = 0.;
    atpos[2][2] = u;

    attyp[3]    = 2;
    atpos[3][0] = 2./3.;
    atpos[3][1] = 1./3.;
    atpos[3][2] = 0.5+u;

    initialized = 1;
    break;
  case 3:
    name = memory->create(name,11,"AB:name");
    strcpy(name, "AB-B4(001)");

    nucell = 8;
    ntype  = 2;

    latvec[0][0] = sqrt(3.0);
    latvec[1][1] = 1.;
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "AB:atpos");
    attyp = memory->create(attyp, nucell, "AB:attyp");

    attyp[0]    = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    attyp[1]    = 1;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;

    attyp[2]    = 1;
    atpos[2][0] = 1./6.;
    atpos[2][1] = 0.5;
    atpos[2][2] = 0.5;

    attyp[3]    = 1;
    atpos[3][0] = 2./3.;
    atpos[3][1] = 0.0;
    atpos[3][2] = 0.5;

    attyp[4]    = 2;
    atpos[4][0] = 0.;
    atpos[4][1] = 0.;
    atpos[4][2] = 0.+u;

    attyp[5]    = 2;
    atpos[5][0] = 0.5;
    atpos[5][1] = 0.5;
    atpos[5][2] = 0.+u;

    attyp[6]    = 2;
    atpos[6][0] = 1./6.;
    atpos[6][1] = 0.5;
    atpos[6][2] = 0.5+u;

    attyp[7]    = 2;
    atpos[7][0] = 2./3.;
    atpos[7][1] = 0.0;
    atpos[7][2] = 0.5+u;

    initialized = 1;
    break;
  case 4:
    name = memory->create(name,11,"AB:name");
    strcpy(name, "AB-B4(001)");
    nucell = 8;
    ntype  = 2;

    latvec[0][0] = 1.;
    latvec[1][1] = sqrt(3.0);
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "AB:atpos");
    attyp = memory->create(attyp, nucell, "AB:attyp");

    attyp[0]    = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    attyp[1]    = 1;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;

    attyp[2]    = 1;
    atpos[2][0] = 0.5;
    atpos[2][1] = 1./6.;
    atpos[2][2] = 0.5;

    attyp[3]    = 1;
    atpos[3][0] = 0.0;
    atpos[3][1] = 2./3.;
    atpos[3][2] = 0.5;

    attyp[4]    = 2;
    atpos[4][0] = 0.;
    atpos[4][1] = 0.;
    atpos[4][2] = u;

    attyp[5]    = 2;
    atpos[5][0] = 0.5;
    atpos[5][1] = 0.5;
    atpos[5][2] = u;

    attyp[6]    = 2;
    atpos[6][0] = 0.5;
    atpos[6][1] = 1./6.;
    atpos[6][2] = 0.5+u;

    attyp[7]    = 2;
    atpos[7][0] = 0.0;
    atpos[7][1] = 2./3.;
    atpos[7][2] = 0.5+u;

    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/* -------------------------------------------------------------------- */
