#include "user.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

#define MAXLINE 256

using namespace std;

/* ----------------------------------------------------------------------
   To select the orientation of the lattice
------------------------------------------------------------------------- */
USER::USER() : lattice()
{
  initialized = 0;
  name = memory->create(name,10,"user:name");
  strcpy(name, "USER_LATT");
  char str[MAXLINE];
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Would you like to read the lattice info from a file?(y/n)[n]: ");
  fgets(str,MAXLINE,stdin);
  char *flag = strtok(str," \t\n\r\f");
  if (flag != NULL && (strcmp(flag,"y")==0 || strcmp(flag,"Y")==0)){
    FILE *fp;
    do printf("\nPlease input the name of the file that contains the lattice info: ");
    while (count_words(fgets(str,MAXLINE,stdin)) < 1);
    char *fname =  strtok(str," \t\n\r\f");
    fp = fopen(fname,"r");
    if (fp == NULL){
      printf("Error: File %s not found!\n", fname);
      return;
    }
    // if read from file, the file would be similar to vasp POSCAR (direct)
    //  alat
    //  xx xy xz
    //  yx yy yz
    //  zx zy zz
    //  ntype1 ntype2 ntype3 ...
    //  sx1 sy1 sz1
    //  ...
    fgets(str,MAXLINE,fp); if (feof(fp)){fclose(fp); return;}
    alat = atof(strtok(str, " \t\n\r\f"));
    for (int i=0; i<3; i++){
      fgets(str,MAXLINE,fp); if (feof(fp)){fclose(fp); return;}
      latvec[i][0] = atof(strtok(str,  " \t\n\r\f"));
      latvec[i][1] = atof(strtok(NULL, " \t\n\r\f"));
      latvec[i][2] = atof(strtok(NULL, " \t\n\r\f"));
    }
    fgets(str,MAXLINE,fp); if (feof(fp)){fclose(fp); return;}
    ntype = count_words(str);

    int *ntm = new int[ntype];
    char *ptr = strtok(str," \t\n\r\f");
    nucell = ntm[0] = atoi(ptr);

    for (int i=1; i<ntype; i++){
      ptr = strtok(NULL," \t\n\r\f");
      ntm[i] = atoi(ptr);
      nucell += ntm[i];
    }

    atpos = memory->create(atpos, nucell, 3, "USER_atpos");
    attyp = memory->create(attyp, nucell, "USER:attyp");

    int iatom =0;
    for (int ip=0; ip<ntype; ip++){
      for (int i=0; i<ntm[ip]; i++){
        fgets(str,MAXLINE,fp); if (feof(fp)){fclose(fp); return;}
        atpos[iatom][0] = atof(strtok(str,  " \t\n\r\f"));
        atpos[iatom][1] = atof(strtok(NULL, " \t\n\r\f"));
        atpos[iatom][2] = atof(strtok(NULL, " \t\n\r\f"));
        attyp[iatom++] = ip+1;
      }
    }
    fclose(fp);
    delete []ntm;
  } else {
    // ask for lattice constant
    alat = 1.;
    printf("\nPlease input the lattice constant of the USER lattice [1.0]: ");
    if (count_words(fgets(str,MAXLINE,stdin))>0) alat = atof(strtok(str, " \t\n\r\f"));

    // ask for lattice vectors
    for (int i=0; i<3; i++){
      while (1){
        printf("Please input the lattice vector A%d: ", i+1);
        if (count_words(fgets(str,MAXLINE,stdin))<3) continue;
        latvec[i][0] = atof(strtok(str,  " \t\n\r\f"));
        latvec[i][1] = atof(strtok(NULL, " \t\n\r\f"));
        latvec[i][2] = atof(strtok(NULL, " \t\n\r\f"));

        break;
      }
    }
    // ask for # of atoms and # of atom types
    nucell = ntype = 1;
    printf("Please input the number of atoms per unit cell [1]: ");
    if (count_words(fgets(str,MAXLINE,stdin)) > 0) nucell = atoi(strtok(str, " \t\n\r\f"));
    if (nucell < 1) return;
  
    if (nucell != 1){
      printf("Please input the number of atom types in cell [1]: ");
      if (count_words(fgets(str,MAXLINE,stdin))>0) ntype = atoi(strtok(str, " \t\n\r\f"));
      if (ntype < 1) return;
      if (ntype > nucell) ntype = nucell;
    }
      
    atpos = memory->create(atpos, nucell, 3, "USER_atpos");
    attyp = memory->create(attyp, nucell, "USER:attyp");
    // ask for atom coordinates and types
    for (int i=0; i<nucell; i++){
      do printf("Please input [type xs ys zs] for atom %d: ", i+1);
      while (count_words(fgets(str,MAXLINE,stdin)) < 4);
      attyp[i]    = atoi(strtok(str,  " \t\n\r\f"));
      atpos[i][0] = atof(strtok(NULL, " \t\n\r\f"));
      atpos[i][1] = atof(strtok(NULL, " \t\n\r\f"));
      atpos[i][2] = atof(strtok(NULL, " \t\n\r\f"));
    }
  }
  for (int i=0; i<14; i++) printf("====="); printf("\n");
  initialized = 1;

return;
}

/* ----------------------------------------------------------------------
   Deconstructor does nothing
------------------------------------------------------------------------- */
USER::~USER()
{

}

/* ------------------------------------------------------------------- */
