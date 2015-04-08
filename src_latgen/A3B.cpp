#include "A3B.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

#define MAXLINE 256

using namespace std;

/* ----------------------------------------------------------------------
   To select the orientation of the lattice
------------------------------------------------------------------------- */
A3B::A3B() : lattice()
{
  char str[MAXLINE];
  alat = 1.; ca = 1.;
  int ctype = 1;
  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Please select the composition of your lattice:\n");
  printf("   1. A3B;\n");
  printf("   2. AB3;\n");
  printf("Your choice [%d]: ", ctype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) ctype = atoi(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d\n", ctype);
  if (ctype == 1){ip1 = 1; ip2=2;}
  else {ip1=2; ip2=1;}

  int lattype = 1;
  printf("Please select the type of your lattice:\n");
  printf("   1. A15;         5. D09;\n");
  printf("   2. D019;        6. L12;\n");
  printf("   3. D022;        7. L60;\n");
  printf("   4. D03;\n");
  printf("Your choice [%d]: ", lattype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) lattype = atoi(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d\n", lattype);

  printf("Please input the lattice constant of the A3B crystal [%g]: ", alat);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) alat = atof(strtok(str, " \t\n\r\f"));
  if (lattype == 2 || lattype == 3 || lattype == 7){
    printf("Please input the c/a or c (negative) of your crystal [%g]: ", ca);
    if (count_words(fgets(str,MAXLINE,stdin)) > 0) ca = atof(strtok(str, " \t\n\r\f"));
    if (ca < 0.) ca = -ca/alat;
  }
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  
  // initialize according to orientation
  initialized = 0;
  switch (lattype){
  case 1:
    A3B_A15();
    break;
  case 2:
    A3B_D019();
    break;
  case 3:
    A3B_D022();
    break;
  case 4:
    A3B_D03();
    break;
  case 5:
    A3B_D09();
    break;
  case 6:
    A3B_L12();
    break;
  case 7:
    A3B_L60();
    break;
  default:
    break;
  }

}

/* ----------------------------------------------------------------------
   Deconstructor does nothing
------------------------------------------------------------------------- */
A3B::~A3B()
{

}

/* ----------------------------------------------------------------------
   Initialize for A15 lattice
------------------------------------------------------------------------- */
void A3B::A3B_A15()
{
  char str[MAXLINE];
  int surftype = 1;
  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Please selection the type of A3B-A15 surface:\n");
  printf("   1. (001) conventional orientation;\n");
  printf("   2. (110), long along y;\n");
  printf("   3. (111), long along x, orthogonal;\n");
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
    name = memory->create(name,13,"A3B:name");
    strcpy(name, "A3B-A15(001)");

    ntype  = 2;
    nucell = 8;
    
    latvec[0][0] = 1.;
    latvec[1][1] = 1.;
    latvec[2][2] = 1.;

    atpos = memory->create(atpos, nucell, 3, "A3B:atpos");
    attyp = memory->create(attyp, nucell, "A3B:attyp");
    
    attyp[0]    = ip2;
    atpos[0][0] = 0.00;
    atpos[0][1] = 0.00;
    atpos[0][2] = 0.00;

    attyp[1]    = ip1;
    atpos[1][0] = 0.50;
    atpos[1][1] = 0.25;
    atpos[1][2] = 0.00;

    attyp[2]    = ip1;
    atpos[2][0] = 0.50;
    atpos[2][1] = 0.75;
    atpos[2][2] = 0.00;

    attyp[3]    = ip1;
    atpos[3][0] = 0.00;
    atpos[3][1] = 0.50;
    atpos[3][2] = 0.25;

    attyp[4]    = ip2;
    atpos[4][0] = 0.50;
    atpos[4][1] = 0.50;
    atpos[4][2] = 0.50;

    attyp[5]    = ip1;
    atpos[5][0] = 0.25;
    atpos[5][1] = 0.00;
    atpos[5][2] = 0.50;

    attyp[6]    = ip1;
    atpos[6][0] = 0.75;
    atpos[6][1] = 0.00;
    atpos[6][2] = 0.50;

    attyp[7]    = ip1;
    atpos[7][0] = 0.00;
    atpos[7][1] = 0.50;
    atpos[7][2] = 0.75;

    initialized = 1;
    break;
  case 2:
    name = memory->create(name,13,"A3B:name");
    strcpy(name, "A3B-A15(110)");

    ntype  = 2;
    nucell = 16;
    
    latvec[0][0] = 1.;
    latvec[1][1] = sqrt(2.);
    latvec[2][2] = sqrt(2.);

    atpos = memory->create(atpos, nucell, 3, "A3B:atpos");
    attyp = memory->create(attyp, nucell, "A3B:attyp");
    
    attyp[0]    = ip2;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    attyp[1]    = ip2;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
    
    attyp[2]    = ip1;
    atpos[2][0] = 0.5;
    atpos[2][1] = 0.125;
    atpos[2][2] = 0.125;
    
    attyp[3]    = ip1;
    atpos[3][0] = 0.;
    atpos[3][1] = 0.375;
    atpos[3][2] = 0.125;
    
    attyp[4]    = ip1;
    atpos[4][0] = 0.75;
    atpos[4][1] = 0.75;
    atpos[4][2] = 0.25;
    
    attyp[5]    = ip1;
    atpos[5][0] = 0.25;
    atpos[5][1] = 0.75;
    atpos[5][2] = 0.25;
    
    attyp[6]    = ip1;
    atpos[6][0] = 0.5;
    atpos[6][1] = 0.375;
    atpos[6][2] = 0.375;
    
    attyp[7]    = ip1;
    atpos[7][0] = 0.;
    atpos[7][1] = 0.125;
    atpos[7][2] = 0.375;
    
    attyp[8]    = ip2;
    atpos[8][0] = 0.;
    atpos[8][1] = 0.5;
    atpos[8][2] = 0.5;
    
    attyp[9]    = ip2;
    atpos[9][0] = 0.5;
    atpos[9][1] = 0.;
    atpos[9][2] = 0.5;

    attyp[10]    = ip1;
    atpos[10][0] = 0.5;
    atpos[10][1] = 0.625;
    atpos[10][2] = 0.625;
    
    attyp[11]    = ip1;
    atpos[11][0] = 0.;
    atpos[11][1] = 0.875;
    atpos[11][2] = 0.625;
    
    attyp[12]    = ip1;
    atpos[12][0] = 0.25;
    atpos[12][1] = 0.25;
    atpos[12][2] = 0.75;
    
    attyp[13]    = ip1;
    atpos[13][0] = 0.75;
    atpos[13][1] = 0.25;
    atpos[13][2] = 0.75;
    
    attyp[14]    = ip1;
    atpos[14][0] = 0.5;
    atpos[14][1] = 0.875;
    atpos[14][2] = 0.875;
    
    attyp[15]    = ip1;
    atpos[15][0] = 0.;
    atpos[15][1] = 0.625;
    atpos[15][2] = 0.875;

    initialized = 1;
    break;
  case 3:
    name = memory->create(name,13,"A3B:name");
    strcpy(name, "A3B-A15(111)");

    ntype  = 2;
    nucell = 24;
    
    latvec[0][0] = sqrt(6.);
    latvec[1][1] = sqrt(2.);
    latvec[2][2] = sqrt(3.)/2.;

    atpos = memory->create(atpos, nucell, 3, "A3B:atpos");
    attyp = memory->create(attyp, nucell, "A3B:attyp");
    
    attyp[0]    = ip2;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1]    = ip2;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;

    attyp[2]    = ip1;
    atpos[2][0] = 1./24.;
    atpos[2][1] = 0.375;
    atpos[2][2] = 1./6.;
    
    attyp[3]    = ip1;
    atpos[3][0] = 13./24.;
    atpos[3][1] = 0.875;
    atpos[3][2] = 1./6.;
    
    attyp[4]    = ip1;
    atpos[4][0] = 7./24.;
    atpos[4][1] = 0.375;
    atpos[4][2] = 1./6.;
    
    attyp[5]    = ip1;
    atpos[5][0] = 19./24.;
    atpos[5][1] = 0.875;
    atpos[5][2] = 1./6.;
    
    attyp[6]    = ip1;
    atpos[6][0] = 2./3.;
    atpos[6][1] = 0.25;
    atpos[6][2] = 1./6.;
    
    attyp[7]    = ip1;
    atpos[7][0] = 1./6.;
    atpos[7][1] = 0.75;
    atpos[7][2] = 1./6.;
    
    attyp[8]    = ip2;
    atpos[8][0] = 1./3.;
    atpos[8][1] = 0.;
    atpos[8][2] = 1./3.;
    
    attyp[9]    = ip2;
    atpos[9][0] = 5./6.;
    atpos[9][1] = 0.5;
    atpos[9][2] = 1./3.;
    
    attyp[10]    = ip1;
    atpos[10][0] = 0.5;
    atpos[10][1] = 0.25;
    atpos[10][2] = 0.5;
    
    attyp[11]    = ip1;
    atpos[11][0] = 0.;
    atpos[11][1] = 0.75;
    atpos[11][2] = 0.5;
    
    attyp[12]    = ip1;
    atpos[12][0] = 0.125;
    atpos[12][1] = 0.125;
    atpos[12][2] = 0.5;
    
    attyp[13]    = ip1;
    atpos[13][0] = 0.625;
    atpos[13][1] = 0.625;
    atpos[13][2] = 0.5;
    
    attyp[14]    = ip1;
    atpos[14][0] = 0.875;
    atpos[14][1] = 0.125;
    atpos[14][2] = 0.5;
    
    attyp[15]    = ip1;
    atpos[15][0] = 0.375;
    atpos[15][1] = 0.625;
    atpos[15][2] = 0.5;
    
    attyp[16]    = ip2;
    atpos[16][0] = 2./3.;
    atpos[16][1] = 0.;
    atpos[16][2] = 2./3.;
    
    attyp[17]    = ip2;
    atpos[17][0] = 1./6.;
    atpos[17][1] = 0.5;
    atpos[17][2] = 2./3.;

    attyp[18]    = ip1;
    atpos[18][0] = 5./6.;
    atpos[18][1] = 0.75;
    atpos[18][2] = 5./6.;
    
    attyp[19]    = ip1;
    atpos[19][0] = 17./24.;
    atpos[19][1] = 0.375;
    atpos[19][2] = 5./6.;
    
    attyp[20]    = ip1;
    atpos[20][0] = 5./24.;
    atpos[20][1] = 0.875;
    atpos[20][2] = 5./6.;
    
    attyp[21]    = ip1;
    atpos[21][0] = 23./24.;
    atpos[21][1] = 0.375;
    atpos[21][2] = 5./6.;
    
    attyp[22]    = ip1;
    atpos[22][0] = 11./24.;
    atpos[22][1] = 0.875;
    atpos[22][2] = 5./6.;
    
    attyp[23]    = ip1;
    atpos[23][0] = 1./3.;
    atpos[23][1] = 0.25;
    atpos[23][2] = 5./6.;
    
    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/* ----------------------------------------------------------------------
   Initialize for D019
------------------------------------------------------------------------- */
void A3B::A3B_D019()
{
  char str[MAXLINE];
  int surftype = 1;
  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Please selection the type of A3B-D019 surface:\n");
  printf("   1. (001), conventional;\n");
  printf("   2. (001), orthogonal, long along y;\n");
  printf("   3. (100), orthogonal;\n");
  printf("   4. (1-10), conventional;\n");
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
    name = memory->create(name,14,"A3B:name");
    strcpy(name, "A3B-D019(001)");

    ntype  = 2;
    nucell = 8;
    
    latvec[0][0] = 1.;
    latvec[1][0] = -0.5;
    latvec[1][1] = sqrt(3.)*0.5;
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "A3B:atpos");
    attyp = memory->create(attyp, nucell, "A3B:attyp");
    
    attyp[0]    = ip2;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    attyp[1]    = ip1;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.;
    atpos[1][2] = 0.;

    attyp[2]    = ip1;
    atpos[2][0] = 0.5;
    atpos[2][1] = 0.5;
    atpos[2][2] = 0.;

    attyp[3]    = ip1;
    atpos[3][0] = 0.;
    atpos[3][1] = 0.5;
    atpos[3][2] = 0.;

    attyp[4]    = ip2;
    atpos[4][0] = 1./3.;
    atpos[4][1] = 2./3.;
    atpos[4][2] = 0.5;

    attyp[5]    = ip1;
    atpos[5][0] = 1./3.;
    atpos[5][1] = 1./6.;
    atpos[5][2] = 0.5;

    attyp[6]    = ip1;
    atpos[6][0] = 5./6.;
    atpos[6][1] = 1./6.;
    atpos[6][2] = 0.5;

    attyp[7]    = ip1;
    atpos[7][0] = 5./6.;
    atpos[7][1] = 2./3.;
    atpos[7][2] = 0.5;

    initialized = 1;
    break;
  case 2:
    name = memory->create(name,14,"A3B:name");
    strcpy(name, "A3B-D019(001)");

    ntype  = 2;
    nucell = 16;
    
    latvec[0][0] = 1.;
    latvec[1][1] = sqrt(3.);
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "A3B:atpos");
    attyp = memory->create(attyp, nucell, "A3B:attyp");
    
    attyp[0]    = ip2;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.5;
    atpos[0][2] = 0.;
    
    attyp[1]    = ip2;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.;
    atpos[1][2] = 0.;
    
    attyp[2]    = ip1;
    atpos[2][0] = 0.25;
    atpos[2][1] = 0.25;
    atpos[2][2] = 0.;
    
    attyp[3]    = ip1;
    atpos[3][0] = 0.75;
    atpos[3][1] = 0.75;
    atpos[3][2] = 0.;
    
    attyp[4]    = ip1;
    atpos[4][0] = 0.25;
    atpos[4][1] = 0.75;
    atpos[4][2] = 0.;
    
    attyp[5]    = ip1;
    atpos[5][0] = 0.75;
    atpos[5][1] = 0.25;
    atpos[5][2] = 0.;
    
    attyp[6]    = ip1;
    atpos[6][0] = 0.;
    atpos[6][1] = 0.;
    atpos[6][2] = 0.;
    
    attyp[7]    = ip1;
    atpos[7][0] = 0.5;
    atpos[7][1] = 0.5;
    atpos[7][2] = 0.;
    
    attyp[8]    = ip2;
    atpos[8][0] = 0.;
    atpos[8][1] = 1./6.;
    atpos[8][2] = 0.5;
    
    attyp[9]    = ip2;
    atpos[9][0] = 0.5;
    atpos[9][1] = 2./3.;
    atpos[9][2] = 0.5;

    attyp[10]    = ip1;
    atpos[10][0] = 0.75;
    atpos[10][1] = 11./12.;
    atpos[10][2] = 0.5;
    
    attyp[11]    = ip1;
    atpos[11][0] = 0.;
    atpos[11][1] = 2./3.;
    atpos[11][2] = 0.5;
    
    attyp[12]    = ip1;
    atpos[12][0] = 0.5;
    atpos[12][1] = 1./6.;
    atpos[12][2] = 0.5;
    
    attyp[13]    = ip1;
    atpos[13][0] = 0.25;
    atpos[13][1] = 11./12.;
    atpos[13][2] = 0.5;
    
    attyp[14]    = ip1;
    atpos[14][0] = 0.75;
    atpos[14][1] = 5./12.;
    atpos[14][2] = 0.5;
    
    attyp[15]    = ip1;
    atpos[15][0] = 0.25;
    atpos[15][1] = 5./12.;
    atpos[15][2] = 0.5;

    initialized = 1;
    break;
  case 3:
    name = memory->create(name,14,"A3B:name");
    strcpy(name, "A3B-D019(100)");

    ntype  = 2;
    nucell = 16;
    
    latvec[0][0] = 1.;
    latvec[1][1] = ca;
    latvec[2][2] = sqrt(3.);

    atpos = memory->create(atpos, nucell, 3, "A3B:atpos");
    attyp = memory->create(attyp, nucell, "A3B:attyp");
    
    attyp[0]    = ip1;
    atpos[0][0] = 0.75;
    atpos[0][1] = 0.5;
    atpos[0][2] = 0.;
    
    attyp[1]    = ip1;
    atpos[1][0] = 0.25;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
    
    attyp[2]    = ip1;
    atpos[2][0] = 0.75;
    atpos[2][1] = 0.;
    atpos[2][2] = 1./6.;
    
    attyp[3]    = ip1;
    atpos[3][0] = 0.25;
    atpos[3][1] = 0.;
    atpos[3][2] = 1./6.;

    attyp[4]    = ip2;
    atpos[4][0] = 0.5;
    atpos[4][1] = 0.5;
    atpos[4][2] = 0.25;
    
    attyp[5]    = ip1;
    atpos[5][0] = 0.;
    atpos[5][1] = 0.5;
    atpos[5][2] = 0.25;
    
    attyp[6]    = ip2;
    atpos[6][0] = 0.;
    atpos[6][1] = 0.;
    atpos[6][2] = 5./12.;
    
    attyp[7]    = ip1;
    atpos[7][0] = 0.5;
    atpos[7][1] = 0.;
    atpos[7][2] = 5./12.;
    
    attyp[8]    = ip1;
    atpos[8][0] = 0.25;
    atpos[8][1] = 0.5;
    atpos[8][2] = 0.5;
    
    attyp[9]    = ip1;
    atpos[9][0] = 0.75;
    atpos[9][1] = 0.5;
    atpos[9][2] = 0.5;
    
    attyp[10]    = ip1;
    atpos[10][0] = 0.25;
    atpos[10][1] = 0.;
    atpos[10][2] = 2./3.;
    
    attyp[11]    = ip1;
    atpos[11][0] = 0.75;
    atpos[11][1] = 0.;
    atpos[11][2] = 2./3.;
    
    attyp[12]    = ip2;
    atpos[12][0] = 0.;
    atpos[12][1] = 0.5;
    atpos[12][2] = 0.75;
    
    attyp[13]    = ip1;
    atpos[13][0] = 0.5;
    atpos[13][1] = 0.5;
    atpos[13][2] = 0.75;
    
    attyp[14]    = ip2;
    atpos[14][0] = 0.5;
    atpos[14][1] = 0.;
    atpos[14][2] = 11./12.;

    attyp[15]    = ip1;
    atpos[15][0] = 0.;
    atpos[15][1] = 0.;
    atpos[15][2] = 11./12.;
    
    initialized = 1;
    break;
  case 4:
    name = memory->create(name,14,"A3B:name");
    strcpy(name, "A3B-D019(110)");

    ntype  = 2;
    nucell = 16;
    
    latvec[0][0] = 1.;
    latvec[1][1] = ca;
    latvec[2][2] = sqrt(3.);

    atpos = memory->create(atpos, nucell, 3, "A3B:atpos");
    attyp = memory->create(attyp, nucell, "A3B:attyp");
    
    attyp[0] = ip1;
    atpos[0][0] = 0.75;
    atpos[0][1] = 0.25;
    atpos[0][2] = 0.;
    
    attyp[1] = ip1;
    atpos[1][0] = 0.25;
    atpos[1][1] = 0.25;
    atpos[1][2] = 0.;
    
    attyp[2] = ip1;
    atpos[2][0] = 0.75;
    atpos[2][1] = 0.75;
    atpos[2][2] = 1./6.;
    
    attyp[3] = ip1;
    atpos[3][0] = 0.25;
    atpos[3][1] = 0.75;
    atpos[3][2] = 1./6.;
    
    attyp[4] = ip2;
    atpos[4][0] = 0.5;
    atpos[4][1] = 0.25;
    atpos[4][2] = 0.25;
    
    attyp[5] = ip1;
    atpos[5][0] = 0.;
    atpos[5][1] = 0.25;
    atpos[5][2] = 0.25;
    
    attyp[6] = ip1;
    atpos[6][0] = 0.5;
    atpos[6][1] = 0.75;
    atpos[6][2] = 5./12.;
    
    attyp[7] = ip2;
    atpos[7][0] = 0.;
    atpos[7][1] = 0.75;
    atpos[7][2] = 5./12.;
    
    attyp[8] = ip1;
    atpos[8][0] = 0.75;
    atpos[8][1] = 0.25;
    atpos[8][2] = 0.5;
    
    attyp[9] = ip1;
    atpos[9][0] = 0.25;
    atpos[9][1] = 0.25;
    atpos[9][2] = 0.5;
    
    attyp[10] = ip1;
    atpos[10][0] = 0.25;
    atpos[10][1] = 0.75;
    atpos[10][2] = 2./3.;
    
    attyp[11] = ip1;
    atpos[11][0] = 0.75;
    atpos[11][1] = 0.75;
    atpos[11][2] = 2./3.;
    
    attyp[12] = ip2;
    atpos[12][0] = 0.;
    atpos[12][1] = 0.25;
    atpos[12][2] = 0.75;
    
    attyp[13] = ip1;
    atpos[13][0] = 0.5;
    atpos[13][1] = 0.25;
    atpos[13][2] = 0.75;
    
    attyp[14] = ip1;
    atpos[14][0] = 0;
    atpos[14][1] = 0.75;
    atpos[14][2] = 11./12.;
    
    attyp[15] = ip2;
    atpos[15][0] = 0.5;
    atpos[15][1] = 0.75;
    atpos[15][2] = 11./12.;

    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/* ----------------------------------------------------------------------
   Initialize for D022 lattice
------------------------------------------------------------------------- */
void A3B::A3B_D022()
{
  char str[MAXLINE];
  int surftype = 1;
  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Please selection the type of A3B-D022 surface:\n");
  printf("   1. (001);\n");
  printf("   2. (100);\n");
  printf("   3. (110), long along y, orthogonal;\n");
  printf("   4. primitive cell;\n");
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
    name = memory->create(name,14,"A3B:name");
    strcpy(name, "A3B-D022(001)");

    ntype  = 2;
    nucell = 8;
    
    latvec[0][0] = 1.;
    latvec[1][1] = 1.;
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "A3B:atpos");
    attyp = memory->create(attyp, nucell, "A3B:attyp");
    
    attyp[0]    = ip2;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1]    = ip1;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
    
    attyp[2]    = ip1;
    atpos[2][0] = 0.;
    atpos[2][1] = 0.5;
    atpos[2][2] = 0.25;
    
    attyp[3]    = ip1;
    atpos[3][0] = 0.5;
    atpos[3][1] = 0.;
    atpos[3][2] = 0.25;
    
    attyp[4]    = ip2;
    atpos[4][0] = 0.5;
    atpos[4][1] = 0.5;
    atpos[4][2] = 0.5;

    attyp[5]    = ip1;
    atpos[5][0] = 0.;
    atpos[5][1] = 0.;
    atpos[5][2] = 0.5;
    
    attyp[6]    = ip1;
    atpos[6][0] = 0.;
    atpos[6][1] = 0.5;
    atpos[6][2] = 0.75;
    
    attyp[7]    = ip1;
    atpos[7][0] = 0.5;
    atpos[7][1] = 0.;
    atpos[7][2] = 0.75;
    
    initialized = 1;
    break;
  case 2:
    name = memory->create(name,14,"A3B:name");
    strcpy(name, "A3B-D022(100)");

    ntype  = 2;
    nucell = 8;
    
    latvec[0][0] = 1.;
    latvec[1][1] = ca;
    latvec[2][2] = 1.;

    atpos = memory->create(atpos, nucell, 3, "A3B:atpos");
    attyp = memory->create(attyp, nucell, "A3B:attyp");
    
    attyp[0]    = ip2;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1]    = ip1;
    atpos[1][0] = 0.;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
    
    attyp[2]    = ip1;
    atpos[2][0] = 0.5;
    atpos[2][1] = 0.25;
    atpos[2][2] = 0.;
    
    attyp[3]    = ip1;
    atpos[3][0] = 0.5;
    atpos[3][1] = 0.75;
    atpos[3][2] = 0.;
    
    attyp[4]    = ip2;
    atpos[4][0] = 0.5;
    atpos[4][1] = 0.5;
    atpos[4][2] = 0.5;

    attyp[5]    = ip1;
    atpos[5][0] = 0.5;
    atpos[5][1] = 0.;
    atpos[5][2] = 0.5;
    
    attyp[6]    = ip1;
    atpos[6][0] = 0.;
    atpos[6][1] = 0.25;
    atpos[6][2] = 0.5;
    
    attyp[7]    = ip1;
    atpos[7][0] = 0.;
    atpos[7][1] = 0.75;
    atpos[7][2] = 0.5;

    initialized = 1;
    break;
  case 3:
    name = memory->create(name,14,"A3B:name");
    strcpy(name, "A3B-D022(110)");

    ntype  = 2;
    nucell = 16;
    
    latvec[0][0] = sqrt(2.);
    latvec[1][1] = ca;
    latvec[2][2] = sqrt(2.);

    atpos = memory->create(atpos, nucell, 3, "A3B:atpos");
    attyp = memory->create(attyp, nucell, "A3B:attyp");
    
    attyp[0]    = ip2;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1]    = ip2;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
    
    attyp[2]    = ip1;
    atpos[2][0] = 0.;
    atpos[2][1] = 0.5;
    atpos[2][2] = 0.;
    
    attyp[3]    = ip1;
    atpos[3][0] = 0.5;
    atpos[3][1] = 0.;
    atpos[3][2] = 0.;
    
    attyp[4]    = ip1;
    atpos[4][0] = 0.25;
    atpos[4][1] = 0.25;
    atpos[4][2] = 0.25;
    
    attyp[5]    = ip1;
    atpos[5][0] = 0.75;
    atpos[5][1] = 0.75;
    atpos[5][2] = 0.25;
    
    attyp[6]    = ip1;
    atpos[6][0] = 0.25;
    atpos[6][1] = 0.75;
    atpos[6][2] = 0.25;
    
    attyp[7]    = ip1;
    atpos[7][0] = 0.75;
    atpos[7][1] = 0.25;
    atpos[7][2] = 0.25;
    
    attyp[8]    = ip2;
    atpos[8][0] = 0.5;
    atpos[8][1] = 0.;
    atpos[8][2] = 0.5;
    
    attyp[9]    = ip2;
    atpos[9][0] = 0.;
    atpos[9][1] = 0.5;
    atpos[9][2] = 0.5;

    attyp[10]    = ip1;
    atpos[10][0] = 0.5;
    atpos[10][1] = 0.5;
    atpos[10][2] = 0.5;
    
    attyp[11]    = ip1;
    atpos[11][0] = 0.;
    atpos[11][1] = 0.;
    atpos[11][2] = 0.5;
    
    attyp[12]    = ip1;
    atpos[12][0] = 0.25;
    atpos[12][1] = 0.25;
    atpos[12][2] = 0.75;
    
    attyp[13]    = ip1;
    atpos[13][0] = 0.75;
    atpos[13][1] = 0.25;
    atpos[13][2] = 0.75;
    
    attyp[14]    = ip1;
    atpos[14][0] = 0.25;
    atpos[14][1] = 0.75;
    atpos[14][2] = 0.75;
    
    attyp[15]    = ip1;
    atpos[15][0] = 0.75;
    atpos[15][1] = 0.75;
    atpos[15][2] = 0.75;

    initialized = 1;
    break;
  case 4:
    name = memory->create(name,19,"A3B:name");
    strcpy(name, "A3B-D022-primitive");

    ntype  = 2;
    nucell = 4;
    
    latvec[0][0] = -0.5;
    latvec[0][1] =  0.5;
    latvec[0][2] =  0.5;
    latvec[1][0] =  0.5;
    latvec[1][1] = -0.5;
    latvec[1][2] =  0.5;
    latvec[2][0] =  0.5;
    latvec[2][1] =  0.5;
    latvec[2][2] = -0.5;

    atpos = memory->create(atpos, nucell, 3, "A3B:atpos");
    attyp = memory->create(attyp, nucell, "A3B:attyp");
    
    attyp[0]    = ip2;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1]    = ip1;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
    
    attyp[2]    = ip1;
    atpos[2][0] = 0.75;
    atpos[2][1] = 0.25;
    atpos[2][2] = 0.5;
    
    attyp[3]    = ip1;
    atpos[3][0] = 0.25;
    atpos[3][1] = 0.75;
    atpos[3][2] = 0.5;
    
    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/* ----------------------------------------------------------------------
   Initialize for D03 lattice
------------------------------------------------------------------------- */
void A3B::A3B_D03()
{
  char str[MAXLINE];
  int surftype = 1;
  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Please selection the type of A3B-D03 surface:\n");
  printf("   1. (001), small;\n");
  printf("   2. (001), conventional;\n");
  printf("   3. (110), long along y, orthogonal;\n");
  printf("   4. (111), long along y, orthogonal;\n");
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
    name = memory->create(name,13,"A3B:name");
    strcpy(name, "A3B-D03(001)");

    ntype  = 2;
    nucell = 8;
    
    latvec[0][0] = 1./sqrt(2.);
    latvec[1][1] = 1./sqrt(2.);
    latvec[2][2] = 1.;

    atpos = memory->create(atpos, nucell, 3, "A3B:atpos");
    attyp = memory->create(attyp, nucell, "A3B:attyp");

    attyp[0] = ip1;
    atpos[0][0] = 0.5;
    atpos[0][1] = 0.5;
    atpos[0][2] = 0.;
    
    attyp[1] = ip2;
    atpos[1][0] = 0.;
    atpos[1][1] = 0.;
    atpos[1][2] = 0.;
    
    attyp[2] = ip1;
    atpos[2][0] = 0.;
    atpos[2][1] = 0.5;
    atpos[2][2] = 0.25;
    
    attyp[3] = ip1;
    atpos[3][0] = 0.5;
    atpos[3][1] = 0.;
    atpos[3][2] = 0.25;
    
    attyp[4] = ip1;
    atpos[4][0] = 0.;
    atpos[4][1] = 0.;
    atpos[4][2] = 0.5;
    
    attyp[5] = ip2;
    atpos[5][0] = 0.5;
    atpos[5][1] = 0.5;
    atpos[5][2] = 0.5;
    
    attyp[6] = ip1;
    atpos[6][0] = 0.5;
    atpos[6][1] = 0.;
    atpos[6][2] = 0.75;
    
    attyp[7] = ip1;
    atpos[7][0] = 0.;
    atpos[7][1] = 0.5;
    atpos[7][2] = 0.75;  

    initialized = 1;
    break;
  case 2:
    name = memory->create(name,14,"A3B:name");
    strcpy(name, "A3B-D03(100)");

    ntype  = 2;
    nucell = 16;
    
    latvec[0][0] = 1.;
    latvec[1][1] = 1.;
    latvec[2][2] = 1.;

    atpos = memory->create(atpos, nucell, 3, "A3B:atpos");
    attyp = memory->create(attyp, nucell, "A3B:attyp");
    
    attyp[0] = ip1;
    atpos[0][0] = 0.5;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = ip1;
    atpos[1][0] = 0.;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
    
    attyp[2] = ip2;
    atpos[2][0] = 0.;
    atpos[2][1] = 0.;
    atpos[2][2] = 0.;
    
    attyp[3] = ip2;
    atpos[3][0] = 0.5;
    atpos[3][1] = 0.5;
    atpos[3][2] = 0.;
    
    attyp[4] = ip1;
    atpos[4][0] = 0.75;
    atpos[4][1] = 0.25;
    atpos[4][2] = 0.25;
    
    attyp[5] = ip1;
    atpos[5][0] = 0.25;
    atpos[5][1] = 0.75;
    atpos[5][2] = 0.25;
    
    attyp[6] = ip1;
    atpos[6][0] = 0.25;
    atpos[6][1] = 0.25;
    atpos[6][2] = 0.25;
    
    attyp[7] = ip1;
    atpos[7][0] = 0.75;
    atpos[7][1] = 0.75;
    atpos[7][2] = 0.25;
    
    attyp[8] = ip1;
    atpos[8][0] = 0.5;
    atpos[8][1] = 0.5;
    atpos[8][2] = 0.5;
    
    attyp[9] = ip1;
    atpos[9][0] = 0.;
    atpos[9][1] = 0.;
    atpos[9][2] = 0.5;
    
    attyp[10] = ip2;
    atpos[10][0] = 0.;
    atpos[10][1] = 0.5;
    atpos[10][2] = 0.5;
    
    attyp[11] = ip2;
    atpos[11][0] = 0.5;
    atpos[11][1] = 0.;
    atpos[11][2] = 0.5;
    
    attyp[12] = ip1;
    atpos[12][0] = 0.75;
    atpos[12][1] = 0.75;
    atpos[12][2] = 0.75;
    
    attyp[13] = ip1;
    atpos[13][0] = 0.75;
    atpos[13][1] = 0.25;
    atpos[13][2] = 0.75;
    
    attyp[14] = ip1;
    atpos[14][0] = 0.25;
    atpos[14][1] = 0.75;
    atpos[14][2] = 0.75;
    
    attyp[15] = ip1;
    atpos[15][0] = 0.25;
    atpos[15][1] = 0.25;
    atpos[15][2] = 0.75;

    initialized = 1;
    break;
  case 3:
    name = memory->create(name,13,"A3B:name");
    strcpy(name, "A3B-D03(110)");

    ntype  = 2;
    nucell = 8;
    
    latvec[0][0] = 1.;
    latvec[1][1] = sqrt(2.);
    latvec[2][2] = sqrt(2.);

    atpos = memory->create(atpos, nucell, 3, "A3B:atpos");
    attyp = memory->create(attyp, nucell, "A3B:attyp");
    
    attyp[0] = ip1;
    atpos[0][0] = 0.5;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = ip1;
    atpos[1][0] = 0.75;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
    
    attyp[2] = ip1;
    atpos[2][0] = 0.25;
    atpos[2][1] = 0.5;
    atpos[2][2] = 0.;
    
    attyp[3] = ip2;
    atpos[3][0] = 0.;
    atpos[3][1] = 0.;
    atpos[3][2] = 0.;
    
    attyp[4] = ip1;
    atpos[4][0] = 0.;
    atpos[4][1] = 0.5;
    atpos[4][2] = 0.5;
    
    attyp[5] = ip1;
    atpos[5][0] = 0.25;
    atpos[5][1] = 0.;
    atpos[5][2] = 0.5;
    
    attyp[6] = ip1;
    atpos[6][0] = 0.75;
    atpos[6][1] = 0.;
    atpos[6][2] = 0.5;
    
    attyp[7] = ip2;
    atpos[7][0] = 0.5;
    atpos[7][1] = 0.5;
    atpos[7][2] = 0.5;

    initialized = 1;
    break;
  case 4:
    name = memory->create(name,13,"A3B:name");
    strcpy(name, "A3B-D03(111)");

    ntype  = 2;
    nucell = 24;
    
    latvec[0][0] = 1./sqrt(2.);
    latvec[1][1] = sqrt(6.)/2.;
    latvec[2][2] = sqrt(3.);

    atpos = memory->create(atpos, nucell, 3, "A3B:atpos");
    attyp = memory->create(attyp, nucell, "A3B:attyp");
    
    attyp[0] = ip2;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = ip2;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
    
    attyp[2] = ip1;
    atpos[2][0] = 0.5;
    atpos[2][1] = 1./6.;
    atpos[2][2] = 1./12.;
    
    attyp[3] = ip1;
    atpos[3][0] = 0.;
    atpos[3][1] = 2./3.;
    atpos[3][2] = 1./12.;
    
    attyp[4] = ip1;
    atpos[4][0] = 0.5;
    atpos[4][1] = 5./6.;
    atpos[4][2] = 1./6.;
    
    attyp[5] = ip1;
    atpos[5][0] = 0.;
    atpos[5][1] = 1./3.;
    atpos[5][2] = 1./6.;
    
    attyp[6] = ip1;
    atpos[6][0] = 0.5;
    atpos[6][1] = 0.5;
    atpos[6][2] = 0.25;
    
    attyp[7] = ip1;
    atpos[7][0] = 0.;
    atpos[7][1] = 0.;
    atpos[7][2] = 0.25;
    
    attyp[8] = ip2;
    atpos[8][0] = 0.;
    atpos[8][1] = 2./3.;
    atpos[8][2] = 1./3.;
    
    attyp[9] = ip2;
    atpos[9][0] = 0.5;
    atpos[9][1] = 1./6.;
    atpos[9][2] = 1./3.;
    
    attyp[10] = ip1;
    atpos[10][0] = 0.5;
    atpos[10][1] = 5./6.;
    atpos[10][2] = 5./12.;
    
    attyp[11] = ip1;
    atpos[11][0] = 0.;
    atpos[11][1] = 1./3.;
    atpos[11][2] = 5./12.;
    
    attyp[12] = ip1;
    atpos[12][0] = 0.5;
    atpos[12][1] = 0.5;
    atpos[12][2] = 0.5;
    
    attyp[13] = ip1;
    atpos[13][0] = 0.;
    atpos[13][1] = 0.;
    atpos[13][2] = 0.5;
    
    attyp[14] = ip1;
    atpos[14][0] = 0.;
    atpos[14][1] = 2./3.;
    atpos[14][2] = 7./12.;
    
    attyp[15] = ip1;
    atpos[15][0] = 0.5;
    atpos[15][1] = 1./6.;
    atpos[15][2] = 7./12.;
    
    attyp[16] = ip2;
    atpos[16][0] = 0.;
    atpos[16][1] = 1./3.;
    atpos[16][2] = 2./3.;
    
    attyp[17] = ip2;
    atpos[17][0] = 0.5;
    atpos[17][1] = 5./6.;
    atpos[17][2] = 2./3.;
    
    attyp[18] = ip1;
    atpos[18][0] = 0.;
    atpos[18][1] = 0.;
    atpos[18][2] = 0.75;
    
    attyp[19] = ip1;
    atpos[19][0] = 0.5;
    atpos[19][1] = 0.5;
    atpos[19][2] = 0.75;
    
    attyp[20] = ip1;
    atpos[20][0] = 0.;
    atpos[20][1] = 2./3.;
    atpos[20][2] = 5./6.;
    
    attyp[21] = ip1;
    atpos[21][0] = 0.5;
    atpos[21][1] = 1./6.;
    atpos[21][2] = 5./6.;
    
    attyp[22] = ip1;
    atpos[22][0] = 0.5;
    atpos[22][1] = 5./6.;
    atpos[22][2] = 11./12.;
    
    attyp[23] = ip1;
    atpos[23][0] = 0.;
    atpos[23][1] = 1./3.;
    atpos[23][2] = 11./12.;

    initialized = 1;
    break;
  case 5:
    name = memory->create(name,18,"A3B:name");
    strcpy(name, "A3B-D03-primitive");

    ntype  = 2;
    nucell = 4;
    
    latvec[0][1] = 0.5;
    latvec[0][2] = 0.5;
    latvec[1][0] = 0.5;
    latvec[1][2] = 0.5;
    latvec[2][0] = 0.5;
    latvec[2][1] = 0.5;

    atpos = memory->create(atpos, nucell, 3, "A3B:atpos");
    attyp = memory->create(attyp, nucell, "A3B:attyp");
    
    attyp[0] = ip2;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = ip1;
    atpos[1][0] = 0.25;
    atpos[1][1] = 0.25;
    atpos[1][2] = 0.25;
    
    attyp[2] = ip1;
    atpos[2][0] = 0.5;
    atpos[2][1] = 0.5;
    atpos[2][2] = 0.5;
    
    attyp[3] = ip1;
    atpos[3][0] = 0.75;
    atpos[3][1] = 0.75;
    atpos[3][2] = 0.75;
    
    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/* ----------------------------------------------------------------------
   Initialize for D09 lattice
------------------------------------------------------------------------- */
void A3B::A3B_D09()
{
  char str[MAXLINE];
  int surftype = 1;
  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Please selection the type of A3B-D09 surface:\n");
  printf("   1. (001), conventional;\n");
  printf("   2. (110), orthogonal, long along y;\n");
  printf("   3. (111), orthogonal, long along y;\n");
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
    name = memory->create(name,13,"A3B:name");
    strcpy(name, "A3B-D09(001)");

    ntype  = 2;
    nucell = 4;
    
    latvec[0][0] = 1.;
    latvec[1][1] = 1.;
    latvec[2][2] = 1.;

    atpos = memory->create(atpos, nucell, 3, "A3B:atpos");
    attyp = memory->create(attyp, nucell, "A3B:attyp");

    attyp[0] = ip1;
    atpos[0][0] = 0.5;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = ip1;
    atpos[1][0] = 0.;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
    
    attyp[2] = ip2;
    atpos[2][0] = 0.;
    atpos[2][1] = 0.;
    atpos[2][2] = 0.;
    
    attyp[3] = ip1;
    atpos[3][0] = 0.;
    atpos[3][1] = 0.;
    atpos[3][2] = 0.5;

    initialized = 1;
    break;
  case 2:
    name = memory->create(name,13,"A3B:name");
    strcpy(name, "A3B-D09(110)");

    ntype  = 2;
    nucell = 8;
    
    latvec[0][0] = 1.;
    latvec[1][1] = sqrt(2.);
    latvec[2][2] = sqrt(2.);

    atpos = memory->create(atpos, nucell, 3, "A3B:atpos");
    attyp = memory->create(attyp, nucell, "A3B:attyp");
    
    attyp[0] = ip1;
    atpos[0][0] = 0.5;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = ip2;
    atpos[1][0] = 0.;
    atpos[1][1] = 0.;
    atpos[1][2] = 0.;
    
    attyp[2] = ip1;
    atpos[2][0] = 0.;
    atpos[2][1] = 0.25;
    atpos[2][2] = 0.25;
    
    attyp[3] = ip1;
    atpos[3][0] = 0.;
    atpos[3][1] = 0.75;
    atpos[3][2] = 0.25;
    
    attyp[4] = ip1;
    atpos[4][0] = 0.5;
    atpos[4][1] = 0.50;
    atpos[4][2] = 0.50;
    
    attyp[5] = ip2;
    atpos[5][0] = 0.;
    atpos[5][1] = 0.50;
    atpos[5][2] = 0.50;
    
    attyp[6] = ip1;
    atpos[6][0] = 0.;
    atpos[6][1] = 0.75;
    atpos[6][2] = 0.75;
    
    attyp[7] = ip1;
    atpos[7][0] = 0.;
    atpos[7][1] = 0.25;
    atpos[7][2] = 0.75;

    initialized = 1;
    break;
  case 3:
    name = memory->create(name,13,"A3B:name");
    strcpy(name, "A3B-D09(111)");

    ntype  = 2;
    nucell = 24;
    
    latvec[0][0] = sqrt(2.);
    latvec[1][1] = sqrt(6.);
    latvec[2][2] = sqrt(3.);

    atpos = memory->create(atpos, nucell, 3, "A3B:atpos");
    attyp = memory->create(attyp, nucell, "A3B:attyp");

    attyp[0] = ip2;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = ip2;
    atpos[1][0] = 0.50;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
    
    attyp[2] = ip1;
    atpos[2][0] = 0.25;
    atpos[2][1] = 5./12.;
    atpos[2][2] = 1./6.;
    
    attyp[3] = ip1;
    atpos[3][0] = 0.25;
    atpos[3][1] = 11./12.;
    atpos[3][2] = 1./6.;
    
    attyp[4] = ip1;
    atpos[4][0] = 0.75;
    atpos[4][1] = 11./12.;
    atpos[4][2] = 1./6.;
    
    attyp[5] = ip1;
    atpos[5][0] = 0.75;
    atpos[5][1] = 5./12.;
    atpos[5][2] = 1./6.;
    
    attyp[6] = ip1;
    atpos[6][0] = 0.;
    atpos[6][1] = 1./6.;
    atpos[6][2] = 1./6.;
    
    attyp[7] = ip1;
    atpos[7][0] = 0.50;
    atpos[7][1] = 2./3.;
    atpos[7][2] = 1./6.;
    
    attyp[8] = ip2;
    atpos[8][0] = 0.;
    atpos[8][1] = 1./3.;
    atpos[8][2] = 1./3.;
    
    attyp[9] = ip2;
    atpos[9][0] = 0.50;
    atpos[9][1] = 5./6.;
    atpos[9][2] = 1./3.;
    
    attyp[10] = ip1;
    atpos[10][0] = 0.50;
    atpos[10][1] = 0.;
    atpos[10][2] = 0.5;
    
    attyp[11] = ip1;
    atpos[11][0] = 0.25;
    atpos[11][1] = 0.75;
    atpos[11][2] = 0.5;
    
    attyp[12] = ip1;
    atpos[12][0] = 0.75;
    atpos[12][1] = 0.25;
    atpos[12][2] = 0.5;
    
    attyp[13] = ip1;
    atpos[13][0] = 0.25;
    atpos[13][1] = 0.25;
    atpos[13][2] = 0.5;
    
    attyp[14] = ip1;
    atpos[14][0] = 0.75;
    atpos[14][1] = 0.75;
    atpos[14][2] = 0.5;
    
    attyp[15] = ip1;
    atpos[15][0] = 0.;
    atpos[15][1] = 0.5;
    atpos[15][2] = 0.5;
    
    attyp[16] = ip2;
    atpos[16][0] = 0.;
    atpos[16][1] = 2./3.;
    atpos[16][2] = 2./3.;
    
    attyp[17] = ip2;
    atpos[17][0] = 0.50;
    atpos[17][1] = 1./6.;
    atpos[17][2] = 2./3.;
    
    attyp[18] = ip1;
    atpos[18][0] = 0.;
    atpos[18][1] = 5./6.;
    atpos[18][2] = 5./6.;
    
    attyp[19] = ip1;
    atpos[19][0] = 0.50;
    atpos[19][1] = 1./3.;
    atpos[19][2] = 5./6.;
    
    attyp[20] = ip1;
    atpos[20][0] = 0.25;
    atpos[20][1] = 1./12.;
    atpos[20][2] = 5./6.;
    
    attyp[21] = ip1;
    atpos[21][0] = 0.75;
    atpos[21][1] = 7./12.;
    atpos[21][2] = 5./6.;
    
    attyp[22] = ip1;
    atpos[22][0] = 0.25;
    atpos[22][1] = 7./12.;
    atpos[22][2] = 5./6.;
    
    attyp[23] = ip1;
    atpos[23][0] = 0.75;
    atpos[23][1] = 1./12.;
    atpos[23][2] = 5./6.;

    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/* ----------------------------------------------------------------------
   Initialize for L12 lattice
------------------------------------------------------------------------- */
void A3B::A3B_L12()
{
  char str[MAXLINE];
  int surftype = 1;
  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Please selection the type of A3B-L12 surface:\n");
  printf("   1. (001), conventional;\n");
  printf("   2. (110), orthogonal, long along y;\n");
  printf("   3. (111), orthogonal, long along y;\n");
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
    name = memory->create(name,13,"A3B:name");
    strcpy(name, "A3B-L12(001)");

    ntype  = 2;
    nucell = 4;
    
    latvec[0][0] = 1.;
    latvec[1][1] = 1.;
    latvec[2][2] = 1.;

    atpos = memory->create(atpos, nucell, 3, "A3B:atpos");
    attyp = memory->create(attyp, nucell, "A3B:attyp");

    attyp[0] = ip2;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] = ip1;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
    
    attyp[2] = ip1;
    atpos[2][0] = 0.;
    atpos[2][1] = 0.5;
    atpos[2][2] = 0.5;
    
    attyp[3] = ip1;
    atpos[3][0] = 0.5;
    atpos[3][1] = 0.;
    atpos[3][2] = 0.5;

    initialized = 1;
    break;
  case 2:
    name = memory->create(name,13,"A3B:name");
    strcpy(name, "A3B-L12(110)");

    ntype  = 2;
    nucell = 8;
    
    latvec[0][0] = 1.;
    latvec[1][1] = sqrt(2.);
    latvec[2][2] = sqrt(2.);

    atpos = memory->create(atpos, nucell, 3, "A3B:atpos");
    attyp = memory->create(attyp, nucell, "A3B:attyp");
    
    attyp[0] = ip1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.50;
    atpos[0][2] = 0.;
    
    attyp[1] = ip2;
    atpos[1][0] = 0.;
    atpos[1][1] = 0.;
    atpos[1][2] = 0.;
    
    attyp[2] = ip1;
    atpos[2][0] = 0.5;
    atpos[2][1] = 0.25;
    atpos[2][2] = 0.25;
    
    attyp[3] = ip1;
    atpos[3][0] = 0.5;
    atpos[3][1] = 0.75;
    atpos[3][2] = 0.25;
    
    attyp[4] = ip1;
    atpos[4][0] = 0.;
    atpos[4][1] = 0.;
    atpos[4][2] = 0.50;
    
    attyp[5] = ip2;
    atpos[5][0] = 0.;
    atpos[5][1] = 0.50;
    atpos[5][2] = 0.50;
    
    attyp[6] = ip1;
    atpos[6][0] = 0.5;
    atpos[6][1] = 0.75;
    atpos[6][2] = 0.75;
    
    attyp[7] = ip1;
    atpos[7][0] = 0.5;
    atpos[7][1] = 0.25;
    atpos[7][2] = 0.75;

    initialized = 1;
    break;
  case 3:
    name = memory->create(name,13,"A3B:name");
    strcpy(name, "A3B-L12(111)");

    ntype  = 2;
    nucell = 24;
    
    latvec[0][0] = sqrt(2.);
    latvec[1][1] = sqrt(6.);
    latvec[2][2] = sqrt(3.);

    atpos = memory->create(atpos, nucell, 3, "A3B:atpos");
    attyp = memory->create(attyp, nucell, "A3B:attyp");

    attyp[0] = ip1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.5;
    atpos[0][2] = 0.;
    
    attyp[1] = ip1;
    atpos[1][0] = 0.25;
    atpos[1][1] = 0.75;
    atpos[1][2] = 0.;
    
    attyp[2] = ip1;
    atpos[2][0] = 0.75;
    atpos[2][1] = 0.25;
    atpos[2][2] = 0.;
    
    attyp[3] = ip1;
    atpos[3][0] = 0.50;
    atpos[3][1] = 0.;
    atpos[3][2] = 0.;
    
    attyp[4] = ip1;
    atpos[4][0] = 0.25;
    atpos[4][1] = 0.25;
    atpos[4][2] = 0.;
    
    attyp[5] = ip1;
    atpos[5][0] = 0.75;
    atpos[5][1] = 0.75;
    atpos[5][2] = 0.;
    
    attyp[6] = ip2;
    atpos[6][0] = 0.;
    atpos[6][1] = 0.;
    atpos[6][2] = 0.;
    
    attyp[7] = ip2;
    atpos[7][0] = 0.50;
    atpos[7][1] = 0.5;
    atpos[7][2] = 0.;
    
    attyp[8] = ip1;
    atpos[8][0] = 0.75;
    atpos[8][1] = 1./12.;
    atpos[8][2] = 1./3.;
    
    attyp[9] = ip1;
    atpos[9][0] = 0.25;
    atpos[9][1] = 1./12.;
    atpos[9][2] = 1./3.;
    
    attyp[10] = ip1;
    atpos[10][0] = 0.75;
    atpos[10][1] = 7./12.;
    atpos[10][2] = 1./3.;
    
    attyp[11] = ip1;
    atpos[11][0] = 0.;
    atpos[11][1] = 5./6.;
    atpos[11][2] = 1./3.;
    
    attyp[12] = ip1;
    atpos[12][0] = 0.50;
    atpos[12][1] = 1./3.;
    atpos[12][2] = 1./3.;
    
    attyp[13] = ip1;
    atpos[13][0] = 0.25;
    atpos[13][1] = 7./12.;
    atpos[13][2] = 1./3.;
    
    attyp[14] = ip2;
    atpos[14][0] = 0.;
    atpos[14][1] = 1./3.;
    atpos[14][2] = 1./3.;
    
    attyp[15] = ip2;
    atpos[15][0] = 0.50;
    atpos[15][1] = 5./6.;
    atpos[15][2] = 1./3.;
    
    attyp[16] = ip1;
    atpos[16][0] = 0.25;
    atpos[16][1] = 11./12.;
    atpos[16][2] = 2./3.;
    
    attyp[17] = ip1;
    atpos[17][0] = 0.75;
    atpos[17][1] = 5./12.;
    atpos[17][2] = 2./3.;
    
    attyp[18] = ip1;
    atpos[18][0] = 0.25;
    atpos[18][1] = 5./12.;
    atpos[18][2] = 2./3.;
    
    attyp[19] = ip1;
    atpos[19][0] = 0.75;
    atpos[19][1] = 11./12.;
    atpos[19][2] = 2./3.;
    
    attyp[20] = ip1;
    atpos[20][0] = 0.;
    atpos[20][1] = 1./6.;
    atpos[20][2] = 2./3.;
    
    attyp[21] = ip1;
    atpos[21][0] = 0.50;
    atpos[21][1] = 2./3.;
    atpos[21][2] = 2./3.;
    
    attyp[22] = ip2;
    atpos[22][0] = 0.;
    atpos[22][1] = 2./3.;
    atpos[22][2] = 2./3.;
    
    attyp[23] = ip2;
    atpos[23][0] = 0.50;
    atpos[23][1] = 1./6.;
    atpos[23][2] = 2./3.;

    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/* ----------------------------------------------------------------------
   Initialize for L60 lattice
------------------------------------------------------------------------- */
void A3B::A3B_L60()
{
  char str[MAXLINE];
  int surftype = 1;
  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Please selection the type of A3B-L60 surface:\n");
  printf("   1. (001), conventional;\n");
  printf("   2. (100), conventional;\n");
  printf("   3. (110), orthogonal, long along y;\n");
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
    name = memory->create(name,13,"A3B:name");
    strcpy(name, "A3B-L60(001)");

    ntype  = 2;
    nucell = 4;
    
    latvec[0][0] = 1.;
    latvec[1][1] = 1.;
    latvec[2][2] = ca;

    atpos = memory->create(atpos, nucell, 3, "A3B:atpos");
    attyp = memory->create(attyp, nucell, "A3B:attyp");

    attyp[0] = ip1;
    atpos[0][0] = 0.5;
    atpos[0][1] = 0.5;
    atpos[0][2] = 0.;
    
    attyp[1] = ip2;
    atpos[1][0] = 0.;
    atpos[1][1] = 0.;
    atpos[1][2] = 0.;
    
    attyp[2] = ip1;
    atpos[2][0] = 0.5;
    atpos[2][1] = 0.;
    atpos[2][2] = 0.5;
    
    attyp[3] = ip1;
    atpos[3][0] = 0.;
    atpos[3][1] = 0.5;
    atpos[3][2] = 0.5;

    initialized = 1;
    break;
  case 2:
    name = memory->create(name,13,"A3B:name");
    strcpy(name, "A3B-L60(100)");

    ntype  = 2;
    nucell = 4;
    
    latvec[0][0] = 1.;
    latvec[1][1] = ca;
    latvec[2][2] = 1.;

    atpos = memory->create(atpos, nucell, 3, "A3B:atpos");
    attyp = memory->create(attyp, nucell, "A3B:attyp");
    
    attyp[0] = ip1;
    atpos[0][0] = 0.5;
    atpos[0][1] = 0.5;
    atpos[0][2] = 0.;
    
    attyp[1] = ip2;
    atpos[1][0] = 0.;
    atpos[1][1] = 0.;
    atpos[1][2] = 0.;
    
    attyp[2] = ip1;
    atpos[2][0] = 0.5;
    atpos[2][1] = 0.;
    atpos[2][2] = 0.5;
    
    attyp[3] = ip1;
    atpos[3][0] = 0.;
    atpos[3][1] = 0.5;
    atpos[3][2] = 0.5;

    initialized = 1;
    break;
  case 3:
    name = memory->create(name,13,"A3B:name");
    strcpy(name, "A3B-L60(110)");

    ntype  = 2;
    nucell = 8;
    
    latvec[0][0] = ca;
    latvec[1][1] = sqrt(2.);
    latvec[2][2] = sqrt(2.);

    atpos = memory->create(atpos, nucell, 3, "A3B:atpos");
    attyp = memory->create(attyp, nucell, "A3B:attyp");

    attyp[0] = ip1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.50;
    atpos[0][2] = 0.;
    
    attyp[1] = ip2;
    atpos[1][0] = 0.;
    atpos[1][1] = 0.;
    atpos[1][2] = 0.;

    attyp[2] = ip1;
    atpos[2][0] = 0.5;
    atpos[2][1] = 0.25;
    atpos[2][2] = 0.25;
    
    attyp[3] = ip1;
    atpos[3][0] = 0.5;
    atpos[3][1] = 0.75;
    atpos[3][2] = 0.25;
    
    attyp[4] = ip1;
    atpos[4][0] = 0.;
    atpos[4][1] = 0.;
    atpos[4][2] = 0.50;
    
    attyp[5] = ip2;
    atpos[5][0] = 0.;
    atpos[5][1] = 0.50;
    atpos[5][2] = 0.50;
    
    attyp[6] = ip1;
    atpos[6][0] = 0.5;
    atpos[6][1] = 0.75;
    atpos[6][2] = 0.75;
    
    attyp[7] = ip1;
    atpos[7][0] = 0.5;
    atpos[7][1] = 0.25;
    atpos[7][2] = 0.75;

    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/* -------------------------------------------------------------------- */
