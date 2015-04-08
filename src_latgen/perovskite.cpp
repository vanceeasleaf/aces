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
void AB::AB_Perov()
{
  char str[MAXLINE];
  int orient = 1;
  printf("Please selection the orientation of the Perovskite lattice:\n");
  printf("   1. Conventional (001);        2. (001), x//[110], y//[-110];\n");
  printf("   3. (110), x//[110];           4. (110), x//[001];\n");
  printf("   5. (111), gamma = 120deg;     6. (111), orthogonal, x//[100];\n");
  printf("   7. (111), orth., x//[110];\n");
  for (int i=0; i<14; i++) printf("-----");
  printf("\nYour choice [%d]: ", orient);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) orient = atoi(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d", orient);
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }

  // some constants
  const double one12 = 1./12., fiv12 = 5./12., sev12 = 7./12., ele12 = 11./12.;
  const double one6 = 1./6., fiv6 = 5./6.;
  const double one3 = 1./3., two3 = 2./3.;

  // initialize according to orientation
  initialized = 0;
  switch (orient){
  case 1:
    name = memory->create(name,11,"Perov:name");
    strcpy(name, "Perov(001)");

    nucell = 5;
    ntype  = 3;
    
    latvec[0][0] = latvec[1][1] = latvec[2][2] = 1.;
    
    atpos = memory->create(atpos,nucell,3,"atpos");
    attyp = memory->create(attyp,nucell,"attyp");
    
    attyp[0] =  3;
    atpos[0][0] = 0.5;
    atpos[0][1] = 0.5;
    atpos[0][2] = 0.5;
    
    attyp[1] =  1;
    atpos[1][0] = 0.0;
    atpos[1][1] = 0.0;
    atpos[1][2] = 0.0;
    
    attyp[2] =  2;
    atpos[2][0] = 0.0;
    atpos[2][1] = 0.5;
    atpos[2][2] = 0.5;
    
    attyp[3] =  2;
    atpos[3][0] = 0.5;
    atpos[3][1] = 0.0;
    atpos[3][2] = 0.5;
    
    attyp[4] =  2;
    atpos[4][0] = 0.5;
    atpos[4][1] = 0.5;
    atpos[4][2] = 0.0;
    
    initialized = 1;
    break;

  case 2:
    name = memory->create(name,11,"Perov:name");
    strcpy(name, "Perov(001)");

    nucell = 10;
    ntype  = 3;
    
    latvec[0][0] = latvec[1][1] = sqrt(2.);
    latvec[2][2] = 1.;
    
    atpos = memory->create(atpos,nucell,3,"atpos");
    attyp = memory->create(attyp,nucell,"attyp");
    
    attyp[0] =  3;
    atpos[0][0] = 0.00;
    atpos[0][1] = 0.50;
    atpos[0][2] = 0.50;
    
    attyp[1] =  3;
    atpos[1][0] = 0.50;
    atpos[1][1] = 0.00;
    atpos[1][2] = 0.50;
    
    attyp[2] =  1;
    atpos[2][0] = 0.00;
    atpos[2][1] = 0.00;
    atpos[2][2] = 0.00;
    
    attyp[3] =  1;
    atpos[3][0] = 0.50;
    atpos[3][1] = 0.50;
    atpos[3][2] = 0.00;
    
    attyp[4] =  2;
    atpos[4][0] = 0.25;
    atpos[4][1] = 0.25;
    atpos[4][2] = 0.50;
    
    attyp[5] =  2;
    atpos[5][0] = 0.75;
    atpos[5][1] = 0.75;
    atpos[5][2] = 0.50;
    
    attyp[6] =  2;
    atpos[6][0] = 0.25;
    atpos[6][1] = 0.75;
    atpos[6][2] = 0.50;
    
    attyp[7] =  2;
    atpos[7][0] = 0.75;
    atpos[7][1] = 0.25;
    atpos[7][2] = 0.50;
    
    attyp[8] =  2;
    atpos[8][0] = 0.00;
    atpos[8][1] = 0.50;
    atpos[8][2] = 0.00;
    
    attyp[9] =  2;
    atpos[9][0] = 0.50;
    atpos[9][1] = 0.00;
    atpos[9][2] = 0.00;
    
    initialized = 1;
    break;

  case 3:
    name = memory->create(name,11,"Perov:name");
    strcpy(name, "Perov(110)");

    nucell = 10;
    ntype  = 3;
    
    latvec[1][1] = 1.;
    latvec[0][0] = latvec[2][2] = sqrt(2.);
    
    atpos = memory->create(atpos,nucell,3,"atpos");
    attyp = memory->create(attyp,nucell,"attyp");
    
    attyp[0] =  3;
    atpos[0][0] = 0.50;
    atpos[0][1] = 0.50;
    atpos[0][2] = 0.00;
    
    attyp[1] =  3;
    atpos[1][0] = 0.00;
    atpos[1][1] = 0.50;
    atpos[1][2] = 0.50;
    
    attyp[2] =  1;
    atpos[2][0] = 0.00;
    atpos[2][1] = 0.00;
    atpos[2][2] = 0.00;
    
    attyp[3] =  1;
    atpos[3][0] = 0.50;
    atpos[3][1] = 0.00;
    atpos[3][2] = 0.50;
    
    attyp[4] =  2;
    atpos[4][0] = 0.75;
    atpos[4][1] = 0.50;
    atpos[4][2] = 0.25;
    
    attyp[5] =  2;
    atpos[5][0] = 0.25;
    atpos[5][1] = 0.50;
    atpos[5][2] = 0.75;
    
    attyp[6] =  2;
    atpos[6][0] = 0.25;
    atpos[6][1] = 0.50;
    atpos[6][2] = 0.25;
    
    attyp[7] =  2;
    atpos[7][0] = 0.75;
    atpos[7][1] = 0.50;
    atpos[7][2] = 0.75;
    
    attyp[8] =  2;
    atpos[8][0] = 0.50;
    atpos[8][1] = 0.00;
    atpos[8][2] = 0.00;
    
    attyp[9] =  2;
    atpos[9][0] = 0.00;
    atpos[9][1] = 0.00;
    atpos[9][2] = 0.50;
    
    initialized = 1;
    break;

  case 4:
    name = memory->create(name,11,"Perov:name");
    strcpy(name, "Perov(110)");

    nucell = 10;
    ntype  = 3;
    
    latvec[0][0] = 1.;
    latvec[1][1] = latvec[2][2] = sqrt(2.);
    
    atpos = memory->create(atpos,nucell,3,"atpos");
    attyp = memory->create(attyp,nucell,"attyp");
    
    attyp[0] =  3;
    atpos[0][0] = 0.50;
    atpos[0][1] = 0.50;
    atpos[0][2] = 0.00;
    
    attyp[1] =  3;
    atpos[1][0] = 0.50;
    atpos[1][1] = 0.00;
    atpos[1][2] = 0.50;
    
    attyp[2] =  1;
    atpos[2][0] = 0.00;
    atpos[2][1] = 0.00;
    atpos[2][2] = 0.00;
    
    attyp[3] =  1;
    atpos[3][0] = 0.00;
    atpos[3][1] = 0.50;
    atpos[3][2] = 0.50;
    
    attyp[4] =  2;
    atpos[4][0] = 0.50;
    atpos[4][1] = 0.75;
    atpos[4][2] = 0.25;
    
    attyp[5] =  2;
    atpos[5][0] = 0.50;
    atpos[5][1] = 0.25;
    atpos[5][2] = 0.75;
    
    attyp[6] =  2;
    atpos[6][0] = 0.50;
    atpos[6][1] = 0.25;
    atpos[6][2] = 0.25;
    
    attyp[7] =  2;
    atpos[7][0] = 0.50;
    atpos[7][1] = 0.75;
    atpos[7][2] = 0.75;
    
    attyp[8] =  2;
    atpos[8][0] = 0.00;
    atpos[8][1] = 0.50;
    atpos[8][2] = 0.00;
    
    attyp[9] =  2;
    atpos[9][0] = 0.00;
    atpos[9][1] = 0.00;
    atpos[9][2] = 0.50;
    
    initialized = 1;
    break;

  case 5:
    name = memory->create(name,11,"Perov:name");
    strcpy(name, "Perov(111)");
    nucell =   15;
    ntype  =  3;
    
    latvec[0][0] = 1.;
    latvec[1][0] = -0.5;
    latvec[1][1] = sqrt(3.)*0.5;
    latvec[2][2] = sqrt(6.)*0.5;
    
    atpos = memory->create(atpos,nucell,3,"atpos");
    attyp = memory->create(attyp,nucell,"attyp");
    
    attyp[ 0] =  3;
    atpos[ 0][0] = one3;
    atpos[ 0][1] = two3;
    atpos[ 0][2] = one6;
    
    attyp[ 1] =  3;
    atpos[ 1][0] = 0.;
    atpos[ 1][1] = 0.;
    atpos[ 1][2] = 0.5;
    
    attyp[ 2] =  3;
    atpos[ 2][0] = two3;
    atpos[ 2][1] = one3;
    atpos[ 2][2] = fiv6;
    
    attyp[ 3] =  1;
    atpos[ 3][0] = 0.;
    atpos[ 3][1] = 0.;
    atpos[ 3][2] = 0.;
    
    attyp[ 4] =  1;
    atpos[ 4][0] = two3;
    atpos[ 4][1] = one3;
    atpos[ 4][2] = one3;
    
    attyp[ 5] =  1;
    atpos[ 5][0] = one3;
    atpos[ 5][1] = two3;
    atpos[ 5][2] = two3;
    
    attyp[ 6] =  2;
    atpos[ 6][0] = 0.;
    atpos[ 6][1] = 0.5;
    atpos[ 6][2] = 0.;
    
    attyp[ 7] =  2;
    atpos[ 7][0] = two3;
    atpos[ 7][1] = fiv6;
    atpos[ 7][2] = one3;
    
    attyp[ 8] =  2;
    atpos[ 8][0] = one3;
    atpos[ 8][1] = one6;
    atpos[ 8][2] = two3;
    
    attyp[ 9] =  2;
    atpos[ 9][0] = 0.5;
    atpos[ 9][1] = 0.5;
    atpos[ 9][2] = 0.;
    
    attyp[10] =  2;
    atpos[10][0] = one6;
    atpos[10][1] = fiv6;
    atpos[10][2] = one3;
    
    attyp[11] =  2;
    atpos[11][0] = fiv6;
    atpos[11][1] = one6;
    atpos[11][2] = two3;
    
    attyp[12] =  2;
    atpos[12][0] = 0.5;
    atpos[12][1] = 0.;
    atpos[12][2] = 0.;
    
    attyp[13] =  2;
    atpos[13][0] = one6;
    atpos[13][1] = one3;
    atpos[13][2] = one3;
    
    attyp[14] =  2;
    atpos[14][0] = fiv6;
    atpos[14][1] = two3;
    atpos[14][2] = two3;
    
    initialized = 1;
    break;

  case 6:
    name = memory->create(name,11,"Perov:name");
    strcpy(name, "Perov(111)");

    nucell = 30;
    ntype  = 3;
    
    latvec[0][0] = 1.;
    latvec[1][1] = sqrt(3.);
    latvec[2][2] = sqrt(6.)*0.5;
    
    atpos = memory->create(atpos,nucell,3,"atpos");
    attyp = memory->create(attyp,nucell,"attyp");
    
    attyp[ 0] =  3;
    atpos[ 0][0] = 0.0;
    atpos[ 0][1] = two3;
    atpos[ 0][2] = one6;
    
    attyp[ 1] =  3;
    atpos[ 1][0] = 0.0;
    atpos[ 1][1] = 0.0;
    atpos[ 1][2] = 0.5;
    
    attyp[ 2] =  3;
    atpos[ 2][0] = 0.5;
    atpos[ 2][1] = one6;
    atpos[ 2][2] = one6;
    
    attyp[ 3] =  3;
    atpos[ 3][0] = 0.5;
    atpos[ 3][1] = 0.5;
    atpos[ 3][2] = 0.5;
    
    attyp[ 4] =  3;
    atpos[ 4][0] = 0.0;
    atpos[ 4][1] = one3;
    atpos[ 4][2] = fiv6;
    
    attyp[ 5] =  3;
    atpos[ 5][0] = 0.5;
    atpos[ 5][1] = fiv6;
    atpos[ 5][2] = fiv6;
    
    attyp[ 6] =  1;
    atpos[ 6][0] = 0.0;
    atpos[ 6][1] = 0.0;
    atpos[ 6][2] = 0.0;
    
    attyp[ 7] =  1;
    atpos[ 7][0] = 0.5;
    atpos[ 7][1] = 0.5;
    atpos[ 7][2] = 0.0;
    
    attyp[ 8] =  1;
    atpos[ 8][0] = 0.0;
    atpos[ 8][1] = one3;
    atpos[ 8][2] = one3;
    
    attyp[ 9] =  1;
    atpos[ 9][0] = 0.5;
    atpos[ 9][1] = fiv6;
    atpos[ 9][2] = one3;
    
    attyp[10] =  1;
    atpos[10][0] = 0.0;
    atpos[10][1] = two3;
    atpos[10][2] = two3;
    
    attyp[11] =  1;
    atpos[11][0] = 0.5;
    atpos[11][1] = one6;
    atpos[11][2] = two3;
    
    attyp[12] =  2;
    atpos[12][0] = 0.25;
    atpos[12][1] = 0.25;
    atpos[12][2] = 0.0;
    
    attyp[13] =  2;
    atpos[13][0] = 0.75;
    atpos[13][1] = 0.75;
    atpos[13][2] = 0.0;
    
    attyp[14] =  2;
    atpos[14][0] = 0.25;
    atpos[14][1] = sev12;
    atpos[14][2] = one3;
    
    attyp[15] =  2;
    atpos[15][0] = 0.25;
    atpos[15][1] = ele12;
    atpos[15][2] = two3;
    
    attyp[16] =  2;
    atpos[16][0] = 0.75;
    atpos[16][1] = one12;
    atpos[16][2] = one3;
    
    attyp[17] =  2;
    atpos[17][0] = 0.75;
    atpos[17][1] = fiv12;
    atpos[17][2] = two3;
    
    attyp[18] =  2;
    atpos[18][0] = 0.0;
    atpos[18][1] = 0.5;
    atpos[18][2] = 0.0;
    
    attyp[19] =  2;
    atpos[19][0] = 0.5;
    atpos[19][1] = 0.0;
    atpos[19][2] = 0.0;
    
    attyp[20] =  2;
    atpos[20][0] = 0.0;
    atpos[20][1] = fiv6;
    atpos[20][2] = one3;
    
    attyp[21] =  2;
    atpos[21][0] = 0.5;
    atpos[21][1] = one3;
    atpos[21][2] = one3;
    
    attyp[22] =  2;
    atpos[22][0] = 0.0;
    atpos[22][1] = one6;
    atpos[22][2] = two3;
    
    attyp[23] =  2;
    atpos[23][0] = 0.5;
    atpos[23][1] = two3;
    atpos[23][2] = two3;
    
    attyp[24] =  2;
    atpos[24][0] = 0.25;
    atpos[24][1] = 0.75;
    atpos[24][2] = 0.0;
    
    attyp[25] =  2;
    atpos[25][0] = 0.75;
    atpos[25][1] = 0.25;
    atpos[25][2] = 0.0;
    
    attyp[26] =  2;
    atpos[26][0] = 0.25;
    atpos[26][1] = one12;
    atpos[26][2] = one3;
    
    attyp[27] =  2;
    atpos[27][0] = 0.75;
    atpos[27][1] = sev12;
    atpos[27][2] = one3;
    
    attyp[28] =  2;
    atpos[28][0] = 0.25;
    atpos[28][1] = fiv12;
    atpos[28][2] = two3;
    
    attyp[29] =  2;
    atpos[29][0] = 0.75;
    atpos[29][1] = ele12;
    atpos[29][2] = two3;
    
    initialized = 1;
    break;

  case 7:
    name = memory->create(name,11,"Perov:name");
    strcpy(name, "Perov(111)");

    nucell = 30;
    ntype  = 3;
    
    latvec[1][1] = 1.;
    latvec[0][0] = sqrt(3.);
    latvec[2][2] = sqrt(6.)*0.5;
    
    atpos = memory->create(atpos,nucell,3,"atpos");
    attyp = memory->create(attyp,nucell,"attyp");
    
    attyp[ 0] =  3;
    atpos[ 0][1] = 0.0;
    atpos[ 0][0] = two3;
    atpos[ 0][2] = one6;
    
    attyp[ 1] =  3;
    atpos[ 1][1] = 0.0;
    atpos[ 1][0] = 0.0;
    atpos[ 1][2] = 0.5;
    
    attyp[ 2] =  3;
    atpos[ 2][1] = 0.5;
    atpos[ 2][0] = one6;
    atpos[ 2][2] = one6;
    
    attyp[ 3] =  3;
    atpos[ 3][1] = 0.5;
    atpos[ 3][0] = 0.5;
    atpos[ 3][2] = 0.5;
    
    attyp[ 4] =  3;
    atpos[ 4][1] = 0.0;
    atpos[ 4][0] = one3;
    atpos[ 4][2] = fiv6;
    
    attyp[ 5] =  3;
    atpos[ 5][1] = 0.5;
    atpos[ 5][0] = fiv6;
    atpos[ 5][2] = fiv6;
    
    attyp[ 6] =  1;
    atpos[ 6][1] = 0.0;
    atpos[ 6][0] = 0.0;
    atpos[ 6][2] = 0.0;
    
    attyp[ 7] =  1;
    atpos[ 7][1] = 0.5;
    atpos[ 7][0] = 0.5;
    atpos[ 7][2] = 0.0;
    
    attyp[ 8] =  1;
    atpos[ 8][1] = 0.0;
    atpos[ 8][0] = one3;
    atpos[ 8][2] = one3;
    
    attyp[ 9] =  1;
    atpos[ 9][1] = 0.5;
    atpos[ 9][0] = fiv6;
    atpos[ 9][2] = one3;
    
    attyp[10] =  1;
    atpos[10][1] = 0.0;
    atpos[10][0] = two3;
    atpos[10][2] = two3;
    
    attyp[11] =  1;
    atpos[11][1] = 0.5;
    atpos[11][0] = one6;
    atpos[11][2] = two3;
    
    attyp[12] =  2;
    atpos[12][1] = 0.25;
    atpos[12][0] = 0.25;
    atpos[12][2] = 0.0;
    
    attyp[13] =  2;
    atpos[13][1] = 0.75;
    atpos[13][0] = 0.75;
    atpos[13][2] = 0.0;
    
    attyp[14] =  2;
    atpos[14][1] = 0.25;
    atpos[14][0] = sev12;
    atpos[14][2] = one3;
    
    attyp[15] =  2;
    atpos[15][1] = 0.25;
    atpos[15][0] = ele12;
    atpos[15][2] = two3;
    
    attyp[16] =  2;
    atpos[16][1] = 0.75;
    atpos[16][0] = one12;
    atpos[16][2] = one3;
    
    attyp[17] =  2;
    atpos[17][1] = 0.75;
    atpos[17][0] = fiv12;
    atpos[17][2] = two3;
    
    attyp[18] =  2;
    atpos[18][1] = 0.0;
    atpos[18][0] = 0.5;
    atpos[18][2] = 0.0;
    
    attyp[19] =  2;
    atpos[19][1] = 0.5;
    atpos[19][0] = 0.0;
    atpos[19][2] = 0.0;
    
    attyp[20] =  2;
    atpos[20][1] = 0.0;
    atpos[20][0] = fiv6;
    atpos[20][2] = one3;
    
    attyp[21] =  2;
    atpos[21][1] = 0.5;
    atpos[21][0] = one3;
    atpos[21][2] = one3;
    
    attyp[22] =  2;
    atpos[22][1] = 0.0;
    atpos[22][0] = one6;
    atpos[22][2] = two3;
    
    attyp[23] =  2;
    atpos[23][1] = 0.5;
    atpos[23][0] = two3;
    atpos[23][2] = two3;
    
    attyp[24] =  2;
    atpos[24][1] = 0.25;
    atpos[24][0] = 0.75;
    atpos[24][2] = 0.0;
    
    attyp[25] =  2;
    atpos[25][1] = 0.75;
    atpos[25][0] = 0.25;
    atpos[25][2] = 0.0;
    
    attyp[26] =  2;
    atpos[26][1] = 0.25;
    atpos[26][0] = one12;
    atpos[26][2] = one3;
    
    attyp[27] =  2;
    atpos[27][1] = 0.75;
    atpos[27][0] = sev12;
    atpos[27][2] = one3;
    
    attyp[28] =  2;
    atpos[28][1] = 0.25;
    atpos[28][0] = fiv12;
    atpos[28][2] = two3;
    
    attyp[29] =  2;
    atpos[29][1] = 0.75;
    atpos[29][0] = ele12;
    atpos[29][2] = two3;
    
    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/* --------------------------------------------------------------------*/
