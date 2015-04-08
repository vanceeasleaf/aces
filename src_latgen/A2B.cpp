#include "A2B.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

#define MAXLINE 256

using namespace std;

/* ----------------------------------------------------------------------
   To select the orientation of the lattice
------------------------------------------------------------------------- */
A2B::A2B() : lattice()
{
  char str[MAXLINE];
  alat = 1.; ca = 1.;
  int ctype = 1;
  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Please select the composition of your lattice:\n");
  printf("   1. A2B;\n");
  printf("   2. AB2;\n");
  printf("Your choice [%d]: ", ctype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) ctype = atoi(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d\n", ctype);
  if (ctype == 1){ip1 = 1; ip2=2;}
  else {ip1=2; ip2=1;}

  int lattype = 1;
  printf("Please select the type of your lattice:\n");
  printf("   1. C1 (Fluorite);\n");
  printf("   2. C15 (Cu2Mg);\n");
  printf("   3. C32 (AlB2);\n");
  printf("Your choice [%d]: ", lattype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) lattype = atoi(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d\n", lattype);

  printf("Please input the lattice constant of the A3B crystal [%g]: ", alat);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) alat = atof(strtok(str, " \t\n\r\f"));
  if (lattype == 3){
    printf("Please input the c/a or c (negative) of your crystal [%g]: ", ca);
    if (count_words(fgets(str,MAXLINE,stdin)) > 0) ca = atof(strtok(str, " \t\n\r\f"));
    if (ca < 0.) ca = -ca/alat;
  }
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  
  // initialize according to orientation
  initialized = 0;
  switch (lattype){
  case 1:
    A2B_C1();
    break;
  case 2:
    A2B_C15();
    break;
  case 3:
    A2B_C32();
    break;
  default:
    break;
  }

}

/* ----------------------------------------------------------------------
   Deconstructor does nothing
------------------------------------------------------------------------- */
A2B::~A2B()
{

}

/* ----------------------------------------------------------------------
   Initialize for C1 (Fluorite) lattice
------------------------------------------------------------------------- */
void A2B::A2B_C1()
{
  char str[MAXLINE];
  int surftype = 1;
  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Please selection the type of A2B-C1 surface:\n");
  printf("   1. (001), small;\n");
  printf("   2. (001), conventional;\n");
  printf("   3. (110), long along y;\n");
  printf("   4. (111), long along x, orthogonal;\n");
  printf("   5. primitive cell;\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = atoi(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d\n", surftype);
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }

  // initialize according to surface type
  switch (surftype){
  case 1:
    name = memory->create(name,12,"A2B:name");
    strcpy(name, "A2B-C1(001)");

    ntype  = 2;
    nucell = 6;
    
    latvec[0][0] = 1./sqrt(2.);
    latvec[1][1] = 1./sqrt(2.);
    latvec[2][2] = 1.;

    atpos = memory->create(atpos, nucell, 3, "A2B:atpos");
    attyp = memory->create(attyp, nucell, "A2B:attyp");

    attyp[0] = ip1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = ip2;
    atpos[1][0] = 0.;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.25;
    
    attyp[2] = ip2;
    atpos[2][0] = 0.5;
    atpos[2][1] = 0.;
    atpos[2][2] = 0.25;
    
    attyp[3] = ip1;
    atpos[3][0] = 0.5;
    atpos[3][1] = 0.5;
    atpos[3][2] = 0.5;
    
    attyp[4] = ip2;
    atpos[4][0] = 0.5;
    atpos[4][1] = 0.;
    atpos[4][2] = 0.75;
    
    attyp[5] = ip2;
    atpos[5][0] = 0.;
    atpos[5][1] = 0.5;
    atpos[5][2] = 0.75;

    initialized = 1;
    break;
  case 2:
    name = memory->create(name,12,"A2B:name");
    strcpy(name, "A2B-C1(001)");

    ntype  = 2;
    nucell = 12;
    
    latvec[0][0] = 1.;
    latvec[1][1] = 1.;
    latvec[2][2] = 1.;

    atpos = memory->create(atpos, nucell, 3, "A2B:atpos");
    attyp = memory->create(attyp, nucell, "A2B:attyp");

    attyp[0] = ip1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = ip1;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
    
    attyp[2] = ip2;
    atpos[2][0] = 0.25;
    atpos[2][1] = 0.25;
    atpos[2][2] = 0.25;
    
    attyp[3] = ip2;
    atpos[3][0] = 0.75;
    atpos[3][1] = 0.75;
    atpos[3][2] = 0.25;
    
    attyp[4] = ip2;
    atpos[4][0] = 0.75;
    atpos[4][1] = 0.25;
    atpos[4][2] = 0.25;
    
    attyp[5] = ip2;
    atpos[5][0] = 0.25;
    atpos[5][1] = 0.75;
    atpos[5][2] = 0.25;
    
    attyp[6] = ip1;
    atpos[6][0] = 0.;
    atpos[6][1] = 0.5;
    atpos[6][2] = 0.5;
    
    attyp[7] = ip1;
    atpos[7][0] = 0.5;
    atpos[7][1] = 0.;
    atpos[7][2] = 0.5;
    
    attyp[8] = ip2;
    atpos[8][0] = 0.75;
    atpos[8][1] = 0.25;
    atpos[8][2] = 0.75;
    
    attyp[9] = ip2;
    atpos[9][0] = 0.25;
    atpos[9][1] = 0.75;
    atpos[9][2] = 0.75;
    
    attyp[10] = ip2;
    atpos[10][0] = 0.25;
    atpos[10][1] = 0.25;
    atpos[10][2] = 0.75;
    
    attyp[11] = ip2;
    atpos[11][0] = 0.75;
    atpos[11][1] = 0.75;
    atpos[11][2] = 0.75;

    initialized = 1;
    break;
  case 3:
    name = memory->create(name,12,"A2B:name");
    strcpy(name, "A2B-C1(110)");

    ntype  = 2;
    nucell = 6;
    
    latvec[0][0] = 1.;
    latvec[1][1] = 1./sqrt(2.);
    latvec[2][2] = 1./sqrt(2.);

    atpos = memory->create(atpos, nucell, 3, "A2B:atpos");
    attyp = memory->create(attyp, nucell, "A2B:attyp");

    attyp[0] = ip2;
    atpos[0][0] = 0.75;
    atpos[0][1] = 0.5;
    atpos[0][2] = 0.;
    
    attyp[1] = ip2;
    atpos[1][0] = 0.25;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
    
    attyp[2] = ip1;
    atpos[2][0] = 0.;
    atpos[2][1] = 0.;
    atpos[2][2] = 0.;
    
    attyp[3] = ip2;
    atpos[3][0] = 0.25;
    atpos[3][1] = 0.;
    atpos[3][2] = 0.5;
    
    attyp[4] = ip2;
    atpos[4][0] = 0.75;
    atpos[4][1] = 0.;
    atpos[4][2] = 0.5;
    
    attyp[5] = ip1;
    atpos[5][0] = 0.5;
    atpos[5][1] = 0.5;
    atpos[5][2] = 0.5;

    initialized = 1;
    break;
  case 4:
    name = memory->create(name,12,"A2B:name");
    strcpy(name, "A2B-C1(111)");

    ntype  = 2;
    nucell = 18;
    
    latvec[0][0] = sqrt(6.)/2.;
    latvec[1][1] = sqrt(2.)/2.;
    latvec[2][2] = sqrt(3.);

    atpos = memory->create(atpos, nucell, 3, "A2B:atpos");
    attyp = memory->create(attyp, nucell, "A2B:attyp");

    attyp[0] = ip1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = ip1;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
    
    attyp[2] = ip2;
    atpos[2][0] = 5./6.;
    atpos[2][1] = 0.5;
    atpos[2][2] = 1./12.;
    
    attyp[3] = ip2;
    atpos[3][0] = 1./3.;
    atpos[3][1] = 0.;
    atpos[3][2] = 1./12.;
    
    attyp[4] = ip2;
    atpos[4][0] = 0.;
    atpos[4][1] = 0.;
    atpos[4][2] = 0.25;
    
    attyp[5] = ip2;
    atpos[5][0] = 0.5;
    atpos[5][1] = 0.5;
    atpos[5][2] = 0.25;
    
    attyp[6] = ip1;
    atpos[6][0] = 5./6.;
    atpos[6][1] = 0.5;
    atpos[6][2] = 1./3.;
    
    attyp[7] = ip1;
    atpos[7][0] = 1./3.;
    atpos[7][1] = 0.;
    atpos[7][2] = 1./3.;
    
    attyp[8] = ip2;
    atpos[8][0] = 1./6.;
    atpos[8][1] = 0.5;
    atpos[8][2] = 5./12.;
    
    attyp[9] = ip2;
    atpos[9][0] = 2./3.;
    atpos[9][1] = 0.;
    atpos[9][2] = 5./12.;
    
    attyp[10] = ip2;
    atpos[10][0] = 1./3.;
    atpos[10][1] = 0.;
    atpos[10][2] = 7./12.;
    
    attyp[11] = ip2;
    atpos[11][0] = 5./6.;
    atpos[11][1] = 0.5;
    atpos[11][2] = 7./12.;
    
    attyp[12] = ip1;
    atpos[12][0] = 2./3.;
    atpos[12][1] = 0.;
    atpos[12][2] = 2./3.;
    
    attyp[13] = ip1;
    atpos[13][0] = 1./6.;
    atpos[13][1] = 0.5;
    atpos[13][2] = 2./3.;
    
    attyp[14] = ip2;
    atpos[14][0] = 0.5;
    atpos[14][1] = 0.5;
    atpos[14][2] = 0.75;
    
    attyp[15] = ip2;
    atpos[15][0] = 0.;
    atpos[15][1] = 0.;
    atpos[15][2] = 0.75;
    
    attyp[16] = ip2;
    atpos[16][0] = 2./3.;
    atpos[16][1] = 0.;
    atpos[16][2] = 11./12.;
    
    attyp[17] = ip2;
    atpos[17][0] = 1./6.;
    atpos[17][1] = 0.5;
    atpos[17][2] = 11./12.;

    initialized = 1;
    break;
  case 5:
    name = memory->create(name,17,"A2B:name");
    strcpy(name, "A2B-C1-primitive");

    ntype  = 2;
    nucell = 3;
    
    latvec[0][1] = 0.5;
    latvec[0][2] = 0.5;
    latvec[1][0] = 0.5;
    latvec[1][2] = 0.5;
    latvec[2][0] = 0.5;
    latvec[2][1] = 0.5;

    atpos = memory->create(atpos, nucell, 3, "A2B:atpos");
    attyp = memory->create(attyp, nucell, "A2B:attyp");

    attyp[0] = ip1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = ip2;
    atpos[1][0] = 0.25;
    atpos[1][1] = 0.25;
    atpos[1][2] = 0.25;
    
    attyp[2] = ip2;
    atpos[2][0] = 0.75;
    atpos[2][1] = 0.75;
    atpos[2][2] = 0.75;
    
    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/* ----------------------------------------------------------------------
   Initialize for C15 (Cu2Mg) lattice
------------------------------------------------------------------------- */
void A2B::A2B_C15()
{
  char str[MAXLINE];
  int surftype = 1;
  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Please selection the type of A2B-C15 surface:\n");
  printf("   1. (001), small;\n");
  printf("   2. (001), conventional;\n");
  printf("   3. (110), long along y;\n");
  printf("   4. (111), long along x, orthogonal;\n");
  printf("   5. primitive cell;\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = atoi(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d\n", surftype);
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }

  // initialize according to surface type
  switch (surftype){
  case 1:
    name = memory->create(name,13,"A2B:name");
    strcpy(name, "A2B-C15(001)");

    ntype  = 2;
    nucell = 12;
    
    latvec[0][0] = 1./sqrt(2.);
    latvec[1][1] = 1./sqrt(2.);
    latvec[2][2] = 1.;

    atpos = memory->create(atpos, nucell, 3, "A2B:atpos");
    attyp = memory->create(attyp, nucell, "A2B:attyp");

    attyp[0] = ip1;
    atpos[0][0] = 0.5;
    atpos[0][1] = 0.5;
    atpos[0][2] = 0.;
    
    attyp[1] = ip2;
    atpos[1][0] = 0.;
    atpos[1][1] = 0.25;
    atpos[1][2] = 0.125;
    
    attyp[2] = ip2;
    atpos[2][0] = 0.;
    atpos[2][1] = 0.75;
    atpos[2][2] = 0.125;
    
    attyp[3] = ip1;
    atpos[3][0] = 0.5;
    atpos[3][1] = 0.;
    atpos[3][2] = 0.25;
    
    attyp[4] = ip2;
    atpos[4][0] = 0.25;
    atpos[4][1] = 0.5;
    atpos[4][2] = 0.375;
    
    attyp[5] = ip2;
    atpos[5][0] = 0.75;
    atpos[5][1] = 0.5;
    atpos[5][2] = 0.375;
    
    attyp[6] = ip1;
    atpos[6][0] = 0.;
    atpos[6][1] = 0.;
    atpos[6][2] = 0.5;
    
    attyp[7] = ip2;
    atpos[7][0] = 0.5;
    atpos[7][1] = 0.75;
    atpos[7][2] = 0.625;
    
    attyp[8] = ip2;
    atpos[8][0] = 0.5;
    atpos[8][1] = 0.25;
    atpos[8][2] = 0.625;
    
    attyp[9] = ip1;
    atpos[9][0] = 0.;
    atpos[9][1] = 0.5;
    atpos[9][2] = 0.75;
    
    attyp[10] = ip2;
    atpos[10][0] = 0.75;
    atpos[10][1] = 0.;
    atpos[10][2] = 0.875;
    
    attyp[11] = ip2;
    atpos[11][0] = 0.25;
    atpos[11][1] = 0.;
    atpos[11][2] = 0.875;

    initialized = 1;
    break;
  case 2:
    name = memory->create(name,13,"A2B:name");
    strcpy(name, "A2B-C15(001)");

    ntype  = 2;
    nucell = 24;
    
    latvec[0][0] = 1.;
    latvec[1][1] = 1.;
    latvec[2][2] = 1.;

    atpos = memory->create(atpos, nucell, 3, "A2B:atpos");
    attyp = memory->create(attyp, nucell, "A2B:attyp");

    attyp[0] = ip1;
    atpos[0][0] = 0.5;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = ip1;
    atpos[1][0] = 0.;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
    
    attyp[2] = ip2;
    atpos[2][0] = 0.125;
    atpos[2][1] = 0.125;
    atpos[2][2] = 0.125;
    
    attyp[3] = ip2;
    atpos[3][0] = 0.375;
    atpos[3][1] = 0.375;
    atpos[3][2] = 0.125;
    
    attyp[4] = ip2;
    atpos[4][0] = 0.625;
    atpos[4][1] = 0.625;
    atpos[4][2] = 0.125;
    
    attyp[5] = ip2;
    atpos[5][0] = 0.875;
    atpos[5][1] = 0.875;
    atpos[5][2] = 0.125;
    
    attyp[6] = ip1;
    atpos[6][0] = 0.25;
    atpos[6][1] = 0.75;
    atpos[6][2] = 0.25;
    
    attyp[7] = ip1;
    atpos[7][0] = 0.75;
    atpos[7][1] = 0.25;
    atpos[7][2] = 0.25;
    
    attyp[8] = ip2;
    atpos[8][0] = 0.875;
    atpos[8][1] = 0.625;
    atpos[8][2] = 0.375;
    
    attyp[9] = ip2;
    atpos[9][0] = 0.625;
    atpos[9][1] = 0.875;
    atpos[9][2] = 0.375;
    
    attyp[10] = ip2;
    atpos[10][0] = 0.125;
    atpos[10][1] = 0.375;
    atpos[10][2] = 0.375;
    
    attyp[11] = ip2;
    atpos[11][0] = 0.375;
    atpos[11][1] = 0.125;
    atpos[11][2] = 0.375;
    
    attyp[12] = ip1;
    atpos[12][0] = 0.5;
    atpos[12][1] = 0.5;
    atpos[12][2] = 0.5;
    
    attyp[13] = ip1;
    atpos[13][0] = 0.;
    atpos[13][1] = 0.;
    atpos[13][2] = 0.5;
    
    attyp[14] = ip2;
    atpos[14][0] = 0.125;
    atpos[14][1] = 0.625;
    atpos[14][2] = 0.625;
    
    attyp[15] = ip2;
    atpos[15][0] = 0.625;
    atpos[15][1] = 0.125;
    atpos[15][2] = 0.625;
    
    attyp[16] = ip2;
    atpos[16][0] = 0.875;
    atpos[16][1] = 0.375;
    atpos[16][2] = 0.625;
    
    attyp[17] = ip2;
    atpos[17][0] = 0.375;
    atpos[17][1] = 0.875;
    atpos[17][2] = 0.625;
    
    attyp[18] = ip1;
    atpos[18][0] = 0.75;
    atpos[18][1] = 0.75;
    atpos[18][2] = 0.75;
    
    attyp[19] = ip1;
    atpos[19][0] = 0.25;
    atpos[19][1] = 0.25;
    atpos[19][2] = 0.75;
    
    attyp[20] = ip2;
    atpos[20][0] = 0.625;
    atpos[20][1] = 0.375;
    atpos[20][2] = 0.875;
    
    attyp[21] = ip2;
    atpos[21][0] = 0.125;
    atpos[21][1] = 0.875;
    atpos[21][2] = 0.875;
    
    attyp[22] = ip2;
    atpos[22][0] = 0.875;
    atpos[22][1] = 0.125;
    atpos[22][2] = 0.875;
    
    attyp[23] = ip2;
    atpos[23][0] = 0.375;
    atpos[23][1] = 0.625;
    atpos[23][2] = 0.875;

    initialized = 1;
    break;
  case 3:
    name = memory->create(name,13,"A2B:name");
    strcpy(name, "A2B-C15(110)");

    ntype  = 2;
    nucell = 12;
    
    latvec[0][0] = 1.;
    latvec[1][1] = 1./sqrt(2.);
    latvec[2][2] = 1./sqrt(2.);

    atpos = memory->create(atpos, nucell, 3, "A2B:atpos");
    attyp = memory->create(attyp, nucell, "A2B:attyp");

    attyp[0] = ip2;
    atpos[0][0] = 0.875;
    atpos[0][1] = 0.75;
    atpos[0][2] = 0.;
    
    attyp[1] = ip2;
    atpos[1][0] = 0.875;
    atpos[1][1] = 0.25;
    atpos[1][2] = 0.;
    
    attyp[2] = ip1;
    atpos[2][0] = 0.5;
    atpos[2][1] = 0.;
    atpos[2][2] = 0.;
    
    attyp[3] = ip1;
    atpos[3][0] = 0.25;
    atpos[3][1] = 0.5;
    atpos[3][2] = 0.;
    
    attyp[4] = ip2;
    atpos[4][0] = 0.125;
    atpos[4][1] = 0.;
    atpos[4][2] = 0.25;
    
    attyp[5] = ip2;
    atpos[5][0] = 0.625;
    atpos[5][1] = 0.5;
    atpos[5][2] = 0.25;
    
    attyp[6] = ip2;
    atpos[6][0] = 0.375;
    atpos[6][1] = 0.25;
    atpos[6][2] = 0.5;
    
    attyp[7] = ip2;
    atpos[7][0] = 0.375;
    atpos[7][1] = 0.75;
    atpos[7][2] = 0.5;
    
    attyp[8] = ip1;
    atpos[8][0] = 0.;
    atpos[8][1] = 0.5;
    atpos[8][2] = 0.5;
    
    attyp[9] = ip1;
    atpos[9][0] = 0.75;
    atpos[9][1] = 0.;
    atpos[9][2] = 0.5;
    
    attyp[10] = ip2;
    atpos[10][0] = 0.625;
    atpos[10][1] = 0.5;
    atpos[10][2] = 0.75;
    
    attyp[11] = ip2;
    atpos[11][0] = 0.125;
    atpos[11][1] = 0.;
    atpos[11][2] = 0.75;

    initialized = 1;
    break;
  case 4:
    name = memory->create(name,13,"A2B:name");
    strcpy(name, "A2B-C15(111)");

    ntype  = 2;
    nucell = 36;
    
    latvec[0][0] = sqrt(2.)/2.;
    latvec[1][1] = sqrt(6.)/2.;
    latvec[2][2] = sqrt(3.);

    atpos = memory->create(atpos, nucell, 3, "A2B:atpos");
    attyp = memory->create(attyp, nucell, "A2B:attyp");

    attyp[0] = ip2;
    atpos[0][0] = 0.5;
    atpos[0][1] = 2./3.;
    atpos[0][2] = 0.;
    
    attyp[1] = ip2;
    atpos[1][0] = 0.25;
    atpos[1][1] = 11./12.;
    atpos[1][2] = 0.;
    
    attyp[2] = ip2;
    atpos[2][0] = 0.;
    atpos[2][1] = 1./6.;
    atpos[2][2] = 0.;
    
    attyp[3] = ip2;
    atpos[3][0] = 0.75;
    atpos[3][1] = 5./12.;
    atpos[3][2] = 0.;
    
    attyp[4] = ip2;
    atpos[4][0] = 0.25;
    atpos[4][1] = 5./12.;
    atpos[4][2] = 0.;
    
    attyp[5] = ip2;
    atpos[5][0] = 0.75;
    atpos[5][1] = 11./12.;
    atpos[5][2] = 0.;
    
    attyp[6] = ip1;
    atpos[6][0] = 0.5;
    atpos[6][1] = 1./6.;
    atpos[6][2] = 0.125;
    
    attyp[7] = ip1;
    atpos[7][0] = 0.;
    atpos[7][1] = 2./3.;
    atpos[7][2] = 0.125;
    
    attyp[8] = ip2;
    atpos[8][0] = 0.;
    atpos[8][1] = 0.;
    atpos[8][2] = 1./6.;
    
    attyp[9] = ip2;
    atpos[9][0] = 0.5;
    atpos[9][1] = 0.5;
    atpos[9][2] = 1./6.;
    
    attyp[10] = ip1;
    atpos[10][0] = 0.5;
    atpos[10][1] = 5./6.;
    atpos[10][2] = 5./24.;
    
    attyp[11] = ip1;
    atpos[11][0] = 0.;
    atpos[11][1] = 1./3.;
    atpos[11][2] = 5./24.;
    
    attyp[12] = ip2;
    atpos[12][0] = 0.75;
    atpos[12][1] = 5./6.;
    atpos[12][2] = 1./3.;
    
    attyp[13] = ip2;
    atpos[13][0] = 0.5;
    atpos[13][1] = 1./3.;
    atpos[13][2] = 1./3.;
    
    attyp[14] = ip2;
    atpos[14][0] = 0.75;
    atpos[14][1] = 1./12.;
    atpos[14][2] = 1./3.;
    
    attyp[15] = ip2;
    atpos[15][0] = 0.25;
    atpos[15][1] = 5./6.;
    atpos[15][2] = 1./3.;
    
    attyp[16] = ip2;
    atpos[16][0] = 0.;
    atpos[16][1] = 5./6.;
    atpos[16][2] = 1./3.;
    
    attyp[17] = ip2;
    atpos[17][0] = 0.25;
    atpos[17][1] = 1./12.;
    atpos[17][2] = 1./3.;
    
    attyp[18] = ip1;
    atpos[18][0] = 0.5;
    atpos[18][1] = 5./6.;
    atpos[18][2] = 11./24.;
    
    attyp[19] = ip1;
    atpos[19][0] = 0.;
    atpos[19][1] = 1./3.;
    atpos[19][2] = 11./24.;
    
    attyp[20] = ip2;
    atpos[20][0] = 0.;
    atpos[20][1] = 2./3.;
    atpos[20][2] = 0.5;
    
    attyp[21] = ip2;
    atpos[21][0] = 0.5;
    atpos[21][1] = 1./6.;
    atpos[21][2] = 0.5;
    
    attyp[22] = ip1;
    atpos[22][0] = 0.;
    atpos[22][1] = 0.;
    atpos[22][2] = 13./24.;
    
    attyp[23] = ip1;
    atpos[23][0] = 0.5;
    atpos[23][1] = 0.5;
    atpos[23][2] = 13./24.;
    
    attyp[24] = ip2;
    atpos[24][0] = 0.5;
    atpos[24][1] = 0.;
    atpos[24][2] = 2./3.;
    
    attyp[25] = ip2;
    atpos[25][0] = 0.;
    atpos[25][1] = 0.5;
    atpos[25][2] = 2./3.;
    
    attyp[26] = ip2;
    atpos[26][0] = 0.25;
    atpos[26][1] = 0.25;
    atpos[26][2] = 2./3.;
    
    attyp[27] = ip2;
    atpos[27][0] = 0.75;
    atpos[27][1] = 0.75;
    atpos[27][2] = 2./3.;
    
    attyp[28] = ip2;
    atpos[28][0] = 0.25;
    atpos[28][1] = 0.75;
    atpos[28][2] = 2./3.;
    
    attyp[29] = ip2;
    atpos[29][0] = 0.75;
    atpos[29][1] = 0.25;
    atpos[29][2] = 2./3.;
    
    attyp[30] = ip1;
    atpos[30][0] = 0.5;
    atpos[30][1] = 0.5;
    atpos[30][2] = 19./24.;
    
    attyp[31] = ip1;
    atpos[31][0] = 0.;
    atpos[31][1] = 0.;
    atpos[31][2] = 19./24.;
    
    attyp[32] = ip2;
    atpos[32][0] = 0.;
    atpos[32][1] = 1./3.;
    atpos[32][2] = 5./6.;
    
    attyp[33] = ip2;
    atpos[33][0] = 0.5;
    atpos[33][1] = 5./6.;
    atpos[33][2] = 5./6.;
    
    attyp[34] = ip1;
    atpos[34][0] = 0.;
    atpos[34][1] = 2./3.;
    atpos[34][2] = 0.875;
    
    attyp[35] = ip1;
    atpos[35][0] = 0.5;
    atpos[35][1] = 1./6.;
    atpos[35][2] = 0.875;

    initialized = 1;
    break;
  case 5:
    name = memory->create(name,18,"A2B:name");
    strcpy(name, "A2B-C15-primitive");

    ntype  = 2;
    nucell = 6;
    
    latvec[0][1] = 0.5;
    latvec[0][2] = 0.5;
    latvec[1][0] = 0.5;
    latvec[1][2] = 0.5;
    latvec[2][0] = 0.5;
    latvec[2][1] = 0.5;

    atpos = memory->create(atpos, nucell, 3, "A2B:atpos");
    attyp = memory->create(attyp, nucell, "A2B:attyp");

    attyp[0] = ip1;
    atpos[0][0] = 0.5;
    atpos[0][1] = 0.5;
    atpos[0][2] = 0.5;

    attyp[1] = ip1;
    atpos[1][0] = 0.;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.5;

    attyp[2] = ip1;
    atpos[2][0] = 0.5;
    atpos[2][1] = 0.;
    atpos[2][2] = 0.5;

    attyp[3] = ip1;
    atpos[3][0] = 0.5;
    atpos[3][1] = 0.5;
    atpos[3][2] = 0.;

    attyp[4] = ip2;
    atpos[4][0] = 0.125;
    atpos[4][1] = 0.125;
    atpos[4][2] = 0.125;

    attyp[5] = ip2;
    atpos[5][0] = 0.875;
    atpos[5][1] = 0.875;
    atpos[5][2] = 0.875;

    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/* ----------------------------------------------------------------------
   Initialize for C32 (AlB2) lattice
------------------------------------------------------------------------- */
void A2B::A2B_C32()
{
  char str[MAXLINE];
  int surftype = 1;
  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Please selection the type of A2B-C32 surface:\n");
  printf("   1. (001), conventional;\n");
  printf("   2. (001), orthogonal;\n");
  printf("   3. (100);\n");
  printf("   4. (110), long along y;\n");
  printf("   5. (1-10), long along y;\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = atoi(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d\n", surftype);
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }

  // initialize according to surface type
  switch (surftype){
  case 1:
    name = memory->create(name,13,"A2B:name");
    strcpy(name, "A2B-C32(001)");

    ntype  = 2;
    nucell = 3;
    
    latvec[0][0] = 1.;
    latvec[1][0] = -0.5;
    latvec[1][1] = sqrt(3.)/2.;
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "A2B:atpos");
    attyp = memory->create(attyp, nucell, "A2B:attyp");

    attyp[0] = ip2;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = ip1;
    atpos[1][0] = 1./3.;
    atpos[1][1] = 2./3.;
    atpos[1][2] = 0.5;
    
    attyp[2] = ip1;
    atpos[2][0] = 2./3.;
    atpos[2][1] = 1./3.;
    atpos[2][2] = 0.5;

    initialized = 1;
    break;
  case 2:
    name = memory->create(name,13,"A2B:name");
    strcpy(name, "A2B-C32(001)");

    ntype  = 2;
    nucell = 6;
    
    latvec[0][0] = 1.;
    latvec[1][1] = sqrt(3.);
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "A2B:atpos");
    attyp = memory->create(attyp, nucell, "A2B:attyp");

    attyp[0] = ip2;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = ip2;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
    
    attyp[2] = ip1;
    atpos[2][0] = 0.;
    atpos[2][1] = 2./3.;
    atpos[2][2] = 0.5;
    
    attyp[3] = ip1;
    atpos[3][0] = 0.5;
    atpos[3][1] = 1./6.;
    atpos[3][2] = 0.5;
    
    attyp[4] = ip1;
    atpos[4][0] = 0.;
    atpos[4][1] = 1./3.;
    atpos[4][2] = 0.5;
    
    attyp[5] = ip1;
    atpos[5][0] = 0.5;
    atpos[5][1] = 5./6.;
    atpos[5][2] = 0.5;

    initialized = 1;
    break;
  case 3:
    name = memory->create(name,13,"A2B:name");
    strcpy(name, "A2B-C32(100)");

    ntype  = 2;
    nucell = 6;
    
    latvec[0][0] = 1.;
    latvec[1][1] = ca;
    latvec[2][2] = sqrt(3.);

    atpos = memory->create(atpos, nucell, 3, "A2B:atpos");
    attyp = memory->create(attyp, nucell, "A2B:attyp");

    attyp[0] = ip2;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = ip1;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 1./6.;
    
    attyp[2] = ip1;
    atpos[2][0] = 0.;
    atpos[2][1] = 0.5;
    atpos[2][2] = 1./3.;
    
    attyp[3] = ip2;
    atpos[3][0] = 0.5;
    atpos[3][1] = 0.;
    atpos[3][2] = 0.5;
    
    attyp[4] = ip1;
    atpos[4][0] = 0.;
    atpos[4][1] = 0.5;
    atpos[4][2] = 2./3.;
    
    attyp[5] = ip1;
    atpos[5][0] = 0.5;
    atpos[5][1] = 0.5;
    atpos[5][2] = 5./6.;

    initialized = 1;
    break;
  case 4:
    name = memory->create(name,13,"A2B:name");
    strcpy(name, "A2B-C32(110)");

    ntype  = 2;
    nucell = 6;
    
    latvec[0][0] = ca;
    latvec[1][1] = sqrt(3.);
    latvec[2][2] = 1.;

    atpos = memory->create(atpos, nucell, 3, "A2B:atpos");
    attyp = memory->create(attyp, nucell, "A2B:attyp");

    attyp[0] = ip2;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = ip1;
    atpos[1][0] = 0.5;
    atpos[1][1] = 1./3.;
    atpos[1][2] = 0.;
    
    attyp[2] = ip1;
    atpos[2][0] = 0.5;
    atpos[2][1] = 2./3.;
    atpos[2][2] = 0.;
    
    attyp[3] = ip2;
    atpos[3][0] = 0.;
    atpos[3][1] = 0.5;
    atpos[3][2] = 0.5;
    
    attyp[4] = ip1;
    atpos[4][0] = 0.5;
    atpos[4][1] = 5./6.;
    atpos[4][2] = 0.5;
    
    attyp[5] = ip1;
    atpos[5][0] = 0.5;
    atpos[5][1] = 1./6.;
    atpos[5][2] = 0.5;

    initialized = 1;
    break;
  case 5:
    name = memory->create(name,14,"A2B:name");
    strcpy(name, "A2B-C32(1-10)");

    ntype  = 2;
    nucell = 6;
    
    latvec[0][0] = 1.;
    latvec[1][1] = ca;
    latvec[2][2] = sqrt(3.);

    atpos = memory->create(atpos, nucell, 3, "A2B_C32_(1-10)_atpos");
    attyp = memory->create(attyp, nucell, "A2B:attyp");

    attyp[0] = ip2;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = ip1;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 1./6.;
    
    attyp[2] = ip1;
    atpos[2][0] = 0.;
    atpos[2][1] = 0.5;
    atpos[2][2] = 1./3.;
    
    attyp[3] = ip2;
    atpos[3][0] = 0.5;
    atpos[3][1] = 0.;
    atpos[3][2] = 0.5;
    
    attyp[4] = ip1;
    atpos[4][0] = 0.;
    atpos[4][1] = 0.5;
    atpos[4][2] = 2./3.;
    
    attyp[5] = ip1;
    atpos[5][0] = 0.5;
    atpos[5][1] = 0.5;
    atpos[5][2] = 5./6.;

    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/* ------------------------------------------------------------------- */
