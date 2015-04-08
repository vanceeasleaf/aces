#include "driver.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <time.h>
#include <math.h>
#include "version.h"

#define ZERO   1.e-8
#define MAXLINE 512
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
/* ----------------------------------------------------------------------
   constructor to initialize
------------------------------------------------------------------------- */
Driver::Driver()
{
  nx = ny = nz = natom = ntype = 0;

  alat = 0.;
  latt = NULL;
  atpos = NULL;
  attyp = NULL;
  random = NULL;
  element = NULL;
  typeID = numtype = NULL;
  xmap = ymap = zmap = umap = NULL;

  memory = new Memory();

  for (int i=0; i<3; i++)
  for (int j=0; j<3; j++) latvec[i][j] = 0.;

  ShowVersion();
  MainMenu();

return;
}

int Driver::ShowMenu(const int flag)
{
  int ltype = 1;
  char str[MAXLINE];

  // print out the menu
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  if (flag>0) printf("Please select the lattice type for lattice: %c\n", flag+'A'-1);
  else printf("Please select the lattice type of your system:\n");
  printf(" 1. FCC/Diamond;               |  4. A3B;\n");
  printf(" 2. BCC;                       |  5. A2B;\n");
  printf(" 3. HCP/Graphene;              |  6. AB & ABXn;\n");
  for (int i=0; i<31; i++) printf("-"); printf("+"); for (int i=0; i<38; i++) printf("-"); printf("\n");
  if (flag != 0){
    printf(" 7. User defined;              |  0. Exit.\n");
  } else {
    printf(" 7. User defined;              |  8. Multi-layer.\n");
    for (int i=0; i<14; i++) printf("-----"); printf("\n");
#ifdef Poly
    printf(" 9. Polycrystal;               | ");
#endif
    printf(" 0. Exit.\n");
  }
  for (int i=0; i<14; i++) printf("-----");
  printf("\nYour choice [1]: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) ltype = atoi(strtok(str," \t\n\r\f"));
  printf("You selected: %d\n", ltype);
  for (int i=0; i<14; i++) printf("====="); printf("\n");

  switch (ltype){
  case 1: latt = new FCC(); break;
  case 2: latt = new BCC(); break;
  case 3: latt = new HCP(); break;
  case 4: latt = new A3B(); break;
  case 5: latt = new A2B(); break;
  case 6: latt = new AB(); break;
  case 7: latt = new USER(); break;
  case 8: if (flag == 0) FormLayers(); break;
#ifdef Poly
  case 9: if (flag == 0) PolyCrystal(); break;
#endif
  default: exit(1);}

  int rflag = 8 - ltype;
  if (flag == 0 && rflag > 0){
    // re-orient the lattice
    printf("Would you like to re-orient the unit cell to comply with LAMMPS? (y/n)[n]: ");
    if (count_words(fgets(str,MAXLINE,stdin)) > 0){
      char *ptr = strtok(str," \t\n\r\f");
      if (strcmp(ptr,"y")==0 || strcmp(ptr,"Y")==0) latt->OrientLattice();
    }
  }

return rflag;
}

void Driver::MainMenu()
{
  if ( ShowMenu(0) > 0 ){
    if (latt->initialized == 0) exit(2);

    name = memory->create(name, strlen(latt->name)+1, "driver->MainMenu:name");
    strcpy(name, latt->name);
    alat = latt->alat;
  
    latt->display();

    generate();
  }
return;
}

/* ----------------------------------------------------------------------
   deconstructor to free memory
------------------------------------------------------------------------- */
Driver::~Driver()
{
  memory->destroy(atpos);
  memory->destroy(attyp);
  memory->destroy(xmap);
  memory->destroy(ymap);
  memory->destroy(zmap);
  memory->destroy(umap);
  memory->destroy(name);

  type2num.clear();

  if (latt != NULL) delete latt;
  if (element != NULL) delete element;

  memory->destroy(typeID);
  memory->destroy(numtype);
  if (random != NULL) delete random;

  delete memory;
}

/* ----------------------------------------------------------------------
   method to generate atomic configurations
------------------------------------------------------------------------- */
void Driver::generate()
{
  char str[MAXLINE];
  int leading_dir = 1;
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Please input the extensions in x, y, and z directions: ");
  if (count_words(fgets(str,MAXLINE,stdin)) < 3) exit(2);
  nx = atoi(strtok(str,  " \t\n\r\f"));
  ny = atoi(strtok(NULL, " \t\n\r\f"));
  nz = atoi(strtok(NULL, " \t\n\r\f"));
  nucell = latt->nucell;
  natom = nx*ny*nz*nucell;
  if (natom < 1) exit(3);

  printf("Your system would be %d x %d x %d with %d atoms.\n",nx,ny,nz,natom);
  printf("Please indicate which direction should goes fast(1:x; other: z)[1]: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) leading_dir = atoi(strtok(str, " \t\n\r\f"));
  for (int i=0; i<14; i++) printf("====="); printf("\n");

  atpos = memory->create(atpos, natom, 3, "driver->generate:atpos");
  attyp = memory->create(attyp, natom, "driver->generate:attyp");
  xmap = memory->create(xmap, natom, "driver->generate:xmap");
  ymap = memory->create(ymap, natom, "driver->generate:ymap");
  zmap = memory->create(zmap, natom, "driver->generate:zmap");
  umap = memory->create(umap, natom, "driver->generate:umap");

  int iatom = 0;
  if ( leading_dir == 1){
    for (int k=0; k<nz; k++){
      for (int j=0; j<ny; j++){
        for (int i=0; i<nx; i++){
          for (int u=0; u<latt->nucell; u++){
            xmap[iatom] = i;
            ymap[iatom] = j;
            zmap[iatom] = k;
            umap[iatom] = u;
            attyp[iatom] = latt->attyp[u];
            atpos[iatom][0] = latt->atpos[u][0] + double(i);
            atpos[iatom][1] = latt->atpos[u][1] + double(j);
            atpos[iatom][2] = latt->atpos[u][2] + double(k);
            iatom++;
          }
        }
      }
    }
  }else{
    for (int i=0; i<nx; i++){
      for (int j=0; j<ny; j++){
        for (int k=0; k<nz; k++){
          for (int u=0; u<latt->nucell; u++){
            xmap[iatom] = i;
            ymap[iatom] = j;
            zmap[iatom] = k;
            umap[iatom] = u;
            attyp[iatom] = latt->attyp[u];
            atpos[iatom][0] = latt->atpos[u][0] + double(i);
            atpos[iatom][1] = latt->atpos[u][1] + double(j);
            atpos[iatom][2] = latt->atpos[u][2] + double(k);
            iatom++;
          }
        }
      }
    }
  }
  // convert fractional coordinate to cartesian
  double tmp[3];
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = latt->latvec[i][j]*latt->alat;
  }
  for (int i=0; i<natom; i++){
    for (int idim=0; idim<3; idim++) tmp[idim] = atpos[i][idim];
    atpos[i][0] = tmp[0]*latvec[0][0] + tmp[1]*latvec[1][0] + tmp[2]*latvec[2][0];
    atpos[i][1] = tmp[0]*latvec[0][1] + tmp[1]*latvec[1][1] + tmp[2]*latvec[2][1];
    atpos[i][2] = tmp[0]*latvec[0][2] + tmp[1]*latvec[1][2] + tmp[2]*latvec[2][2];
  }
  for (int i=0; i<3; i++){
    latvec[0][i] *= double(nx);
    latvec[1][i] *= double(ny);
    latvec[2][i] *= double(nz);
  }
  // find the total # of types and # of atoms for each type
  typescan();

return;
}

/* ----------------------------------------------------------------------
   method to find the total # of types and # of atoms for each type
------------------------------------------------------------------------- */
void Driver::typescan()
{
  // allocate memory
  int typmax = 10;
  if (typeID != NULL) memory->destroy(typeID);
  if (numtype!= NULL) memory->destroy(numtype);
  typeID  = memory->create(typeID,  typmax, "driver->typescan:typeID");
  numtype = memory->create(numtype, typmax, "driver->typescan:numtype");
  for (int i=0; i<typmax; i++) numtype[i] = 0;

  ntype = 0;
  // now to identify the total number of types
  for (int i=0; i<natom; i++){
    int id = lookup(attyp[i]);
    if (id < 0){
      if (ntype == typmax){
        typmax += 5;
        typeID  = memory->grow(typeID, typmax, "driver->typescan:typeID");
        numtype = memory->grow(numtype,typmax, "driver->typescan:numtype");
      }
      typeID[ntype] = attyp[i];
      id            = ntype++;
    }
    numtype[id]++;
  }
return;
}

/* ----------------------------------------------------------------------
   method to reset the atomic type ID
------------------------------------------------------------------------- */
void Driver::ResetTypeID()
{
  char str[MAXLINE];
  printf("\n"); for (int i=0; i<14; i++) printf("=====");
  printf("\nThere are %d atomic types in system, and their IDs are:\nIndex : ", ntype);
  for (int i=0; i<ntype; i++) printf("%4d", i+1); printf("\nTypeID: ");
  for (int i=0; i<ntype; i++) printf("%4d", typeID[i]); printf("\nNatTyp: ");
  for (int i=0; i<ntype; i++) printf("%4d", numtype[i]);
  printf("\nPlease input the new type IDs in sequence, enter to skip: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0){
    int newID[ntype], num=0;
    char *ptr;
    ptr = strtok(str, " \t\n\r\f");
    while (ptr != NULL){
      newID[num] = atoi(ptr);
      if (++num >= ntype) break;
      ptr = strtok(NULL, " \t\n\r\f");
    }
    if (num == ntype){
      printf("\nThe new type IDs are:\nIndex : ");
      for (int i=0; i<ntype; i++) printf("%4d", i+1); printf("\nTypeID: ");
      for (int i=0; i<ntype; i++) printf("%4d", newID[i]);
      printf("\nIs this what you want? (y/n)[y]: ");
      int nw = count_words(fgets(str,MAXLINE,stdin));
      char *flag = strtok(str, " \t\n\r\f");
      if ( nw < 1 || ((strcmp(flag,"n") != 0) && (strcmp(flag,"N") != 0))){
        for (int i=0; i<natom; i++){
          int id = lookup(attyp[i]);
          attyp[i] = newID[id];
        }

        typescan();
      }
    }
  }
  for (int i=0; i<14; i++) printf("====="); printf("\n");

return;
}

/* ----------------------------------------------------------------------
 * method to map atomic type to real element
 * ---------------------------------------------------------------------- */
void Driver::MapElement()
{
  char str[MAXLINE];
  printf("\n"); for (int i=0; i<14; i++) printf("=====");
  printf("\nThere are %d atomic types in system, and their IDs are:\nIndex : ", ntype);
  for (int i=0; i<ntype; i++) printf("%4d", i+1); printf("\nTypeID: ");
  for (int i=0; i<ntype; i++) printf("%4d", typeID[i]); printf("\nNatTyp: ");
  for (int i=0; i<ntype; i++) printf("%4d", numtype[i]);
  printf("\nPlease input the element symbols in sequence, enter to skip: ");
  if (count_words(fgets(str,MAXLINE,stdin)) >= ntype){

    if (element == NULL) element = new ChemElements();
    type2num.clear();

    char *ptr = strtok(str, " \t\n\r\f");
    for (int i=0; i<ntype; i++){
      int ip = typeID[i];
      type2num[ip] = element->Name2Num(ptr);

      ptr = strtok(NULL, " \t\n\r\f");
    }
    printf("Element:");
    for (int i=0; i<ntype; i++){
      char ename[3];
      int ip = typeID[i]; element->Num2Name(type2num[ip],ename);
      printf("%4s", ename);
    } printf("\n");
  }
  for (int i=0; i<14; i++) printf("====="); printf("\n");

return;
}

/* ----------------------------------------------------------------------
   method to find the ID of the current atomic type
------------------------------------------------------------------------- */
int Driver::lookup(int ip)
{
  for (int i=0; i<ntype; i++){if (ip == typeID[i]) return i;}

  return -1;
}

/* ----------------------------------------------------------------------
   method to write out atomic configuraton and mapping info
------------------------------------------------------------------------- */
void Driver::write()
{
  if (natom < 1) return;
  FILE *fp;
  char str[MAXLINE], *posfile, *mapfile, *lmpfile;
  int flag_lmp_data = 1;
  if (latvec[0][1]*latvec[0][1]+latvec[0][2]*latvec[0][2]+latvec[1][2]*latvec[1][2] > 1.e-6) flag_lmp_data = 0;

  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Please input the filename of the output xyz file [atomcfg.xyz]: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0){
    int n = strlen(str) + 1;
    posfile = new char[n];
    strcpy(posfile, strtok(str," \t\n\r\f"));
  } else {
    posfile = new char[12];
    strcpy(posfile, "atomcfg.xyz");
  }
  if (flag_lmp_data){
    printf("Please input the filename of the lammps atomic file [data.pos]: ");
    if (count_words(fgets(str,MAXLINE,stdin)) > 0){
      int n = strlen(str) + 1;
      lmpfile = new char[n];
      strcpy(lmpfile, strtok(str," \t\n\r\f"));
    } else {
      lmpfile = new char[9];
      strcpy(lmpfile, "data.pos");
    }
  }

  if (xmap){
    printf("Please input the filename of the output map file [map.in]: ");
    if (count_words(fgets(str,MAXLINE,stdin)) > 0){
      int n = strlen(str) + 1;
      mapfile = new char[n];
      strcpy(mapfile, strtok(str," \t\n\r\f"));
    } else {
      mapfile = new char[7];
      strcpy(mapfile, "map.in");
    } 
  }

  printf("\nThe atomic configuration will be written to: %s", posfile);
  if (flag_lmp_data) printf(" and %s\n", lmpfile); else printf("\n");
  if (xmap) printf("The FFT map information  will be written to file: %s\n", mapfile);
  if (natom < 3){
    printf("\nThe basis vectors of your system is:\n");
    for (int i=0; i<3; i++) printf("  %16.16e %16.16e %16.16e\n", latvec[i][0], latvec[i][1], latvec[i][2]);
  }
  for (int i=0; i<14; i++) printf("====="); printf("\n");

  // write the xyz position file 
  fp = fopen(posfile, "w");
  fprintf(fp, "%d\n", natom);
  fprintf(fp, "%s cell with dimension %dx%dx%d and a= %g\n", name, nx, ny, nz, alat);
  int nr = 3;
  if (natom < nr) nr = natom;
  if (type2num.size() == ntype){
    char ename[3];
    for (int i=0; i<nr; i++){
      int ip = attyp[i]; element->Num2Name(type2num[ip], ename);
      fprintf(fp,"%2s %16.16e %16.16e %16.16e crystal_vector %d %16.16e %16.16e %16.16e\n", ename, atpos[i][0],
      atpos[i][1], atpos[i][2], i+1, latvec[i][0], latvec[i][1], latvec[i][2]);
    }
    for (int i=nr; i<natom; i++){
      int ip = attyp[i]; element->Num2Name(type2num[ip], ename);
      fprintf(fp,"%2s %16.16e %16.16e %16.16e\n", ename, atpos[i][0], atpos[i][1], atpos[i][2]);
    }
  } else {
    for (int i=0; i<nr; i++){
      fprintf(fp,"%d %16.16e %16.16e %16.16e crystal_vector %d %16.16e %16.16e %16.16e\n", attyp[i], atpos[i][0],
      atpos[i][1], atpos[i][2], i+1, latvec[i][0], latvec[i][1], latvec[i][2]);
    }
    for (int i=nr; i<natom; i++) fprintf(fp,"%d %16.16e %16.16e %16.16e\n", attyp[i], atpos[i][0], atpos[i][1], atpos[i][2]);
  }
  fclose(fp);
  delete []posfile;

  // write the lammps atomic style file
  if (flag_lmp_data){
    fp = fopen(lmpfile,"w");
    fprintf(fp, "# %s cell with dimension %d x %d x %d and a = %g\n", name, nx, ny, nz, alat);
    fprintf(fp, "%10d  atoms\n", natom);
    fprintf(fp, "%10d  atom types\n\n", ntype);
    fprintf(fp, " 0. %20.14f  xlo xhi\n", latvec[0][0]);
    fprintf(fp, " 0. %20.14f  ylo yhi\n", latvec[1][1]);
    fprintf(fp, " 0. %20.14f  zlo zhi\n", latvec[2][2]);
    if ( latvec[1][0]*latvec[1][0] + latvec[2][0]*latvec[2][0] + latvec[2][1]*latvec[2][1] > 1.e-8 )
      fprintf(fp, "%20.14f %20.14f %20.14f xy xz yz\n", latvec[1][0], latvec[2][0], latvec[2][1]);

    // write atomic mass info (g/mol) if element mapping is done
    if (type2num.size() == ntype){
      fprintf(fp, "\nMasses\n\n");

      for (std::map<int,int>::iterator it = type2num.begin(); it != type2num.end(); it++){
        int ip = (*it).first; int num = (*it).second;
        fprintf(fp,"%d %g\n", ip, element->Num2Mass(num));
      }
    }

    fprintf(fp, "\nAtoms\n\n");
  
    for (int i=0; i<natom; i++) fprintf(fp,"%d %d %20.14f %20.14f %20.14f\n", i+1, attyp[i], atpos[i][0], atpos[i][1], atpos[i][2]);
    fclose(fp);
    delete []lmpfile;
  }

  // write the map file, useful to fix_phonon only.
  if (xmap){
     fp = fopen(mapfile, "w");
     fprintf(fp,"%d %d %d %d\n", nx, ny, nz, nucell);
     fprintf(fp,"Map file for %dx%dx%d %s cell.\n",nx,ny,nz, name);
     for (int i=0; i<natom; i++)
       fprintf(fp,"%d %d %d %d %d\n", xmap[i], ymap[i], zmap[i], umap[i], i+1);
     fclose(fp);

     delete []mapfile;
  }

return;
}

/* ----------------------------------------------------------------------
   method to modify the resultant model
------------------------------------------------------------------------- */
void Driver::modify()
{
  if (natom < 1) return;
  char str[MAXLINE];
  int ncycle = 1;
  while (ncycle){
    int job=0;
    // to display the menu for modification
    printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
    printf("Please select the modification you want to do:\n");
    printf("  1. Create substitutional solid solution;\n");
    printf("  2. Reset atomic types;\n");
    printf("  3. Map atomic type to elements;\n");

    if (ncycle == 1) printf("  0. Nothing.\n");
    else printf("  0. Done.\n");
    printf("Your choice [0]: ");

    if (count_words(fgets(str,MAXLINE,stdin)) >0) job = atoi(strtok(str, " \t\n\r\f"));

    printf("You selected: %d", job);
    printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
   
    switch (job){ 
    case 1: solidsol();    break;
    case 2: ResetTypeID(); break;
    case 3: MapElement();  break;
    default: return;
    }

    ncycle++;
  }
return;
}

/* ----------------------------------------------------------------------
   private method to create substitutional solid solution
------------------------------------------------------------------------- */
void Driver::solidsol()
{
  char str[MAXLINE];
  int job = 0;
  printf("\n"); for (int i=0; i<14; i++) printf("====="); printf("\n");
  printf("Please select the region to create the solid solution:\n");
  printf("  1. Limit solid solution inside/outside a block region;\n");
  printf("  2. Limit solid solution inside/outside a spherical region;\n");
  printf("  3. Union of 1 & 2;\n");
  printf("  4. Intersection of 1 & 2;\n");
  printf("  5. All atoms in the box;\n");
  printf("  0. Return;\nYour choice [%d]: ", job);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) job = atoi(strtok(str," \t\n\r\f"));
  printf("You selected: %d\n", job);
  if (job < 1 || job > 5){
    for (int i=0; i<14; i++) printf("====="); printf("\n");
    return;
  }

  double block[6], cpos[3], radius;
  int flag, logand = 1, nsel, atsel[natom];
  for (int id = 0; id<natom; id++) atsel[id] = 1;
  if (job==1 || job==3 || job==4) flag |= 1;
  if (job==2 || job==3 || job==4) flag |= 2;
  if (job==3) logand = 0;

  if (flag & 1){ // block region needed
    printf("\nThe basis vectors that define the current simulation box is:\n");
    for (int idim=0; idim<3; idim++) printf("  A%d: %lg %lg %lg\n", idim, latvec[idim][0], latvec[idim][1], latvec[idim][2]); 
    printf("Please input the bounds of the block region, in the format of `xlo xhi ylo yhi zlo zhi pos`,\n");
    printf("where `pos` is either `i` or `o`, indicating inside or outside the block. If no limit in certain\n");
    printf("direction, use `NULL`. For non-orthogonal box, this might not work well.\nPlease input them now: ");
    if (count_words(fgets(str,MAXLINE,stdin)) >= 7){
      char *ptr = strtok(str," \t\n\r\f");
      for (int i=0; i<6; i++){
        if (strcmp(ptr, "NULL") == 0) block[i] = double(i%2)*latvec[i/2][i/2];
        else block[i] = atof(ptr);

        ptr = strtok(NULL, " \n\t\r\f");
      }
      int inside = 1;
      if (strcmp(ptr,"o")==0 || strcmp(ptr,"O")==0){
        inside = 0;
        printf("\nAtoms outside block region bounded by [%g %g], [%g %g], [%g %g]\n",
        block[0], block[1], block[2], block[3], block[4], block[5]);
      } else {
        printf("\nAtoms inside block region bounded by [%g %g], [%g %g], [%g %g]\n",
        block[0], block[1], block[2], block[3], block[4], block[5]);
      }
      printf("will be used to create the substitutional solid solution.\n");

      for (int i=0; i<natom; i++){
        if (inside){
          if (atpos[i][0] - block[0] < ZERO || atpos[i][0] - block[1] >-ZERO||
              atpos[i][1] - block[2] < ZERO || atpos[i][1] - block[3] >-ZERO||
              atpos[i][2] - block[4] < ZERO || atpos[i][2] - block[5] >-ZERO) atsel[i] = 0;
        } else {
          if (atpos[i][0] - block[0] >-ZERO && atpos[i][0] - block[1] < ZERO &&
              atpos[i][1] - block[2] >-ZERO && atpos[i][1] - block[3] < ZERO &&
              atpos[i][2] - block[4] >-ZERO && atpos[i][2] - block[5] < ZERO ) atsel[i] = 0;
        }
      }
    } else {
      printf("\nNo input read, operation terminated!\n");
      for (int i=0; i<14; i++) printf("====="); printf("\n");
    }
  }

  nsel = 0;
  for (int i=0; i<natom; i++) nsel += atsel[i];

  if (flag & 2){ // spherical region needed
    printf("\nThe basis vectors that define the current simulation box is:\n");
    for (int idim=0; idim<3; idim++) printf("A%d: %lg %lg %lg\n", idim, latvec[idim][0], latvec[idim][1], latvec[idim][2]); 
    printf("There are %d atoms in current selection. Please input the necessary parameters that define the\n", nsel);
    printf("spherical region in the format of `x y z r pos`, where `pos` is either `i` or `o`, indicating\n");
    printf("inside or outside the block. Please input them now: ");
    if (count_words(fgets(str,MAXLINE,stdin)) >= 5){
      char *ptr = strtok(str," \t\n\r\f");
      for (int i=0; i<3; i++){
        cpos[i] = atof(ptr);
        ptr = strtok(NULL, " \n\t\r\f");
      }
      radius = atof(ptr);

      int inside = 1;
      ptr = strtok(NULL, " \n\t\r\f");
      if (strcmp(ptr,"o")==0 || strcmp(ptr,"O")==0){
        inside = 0;
        printf("\nAtoms outside sphere centered at [%g %g %g] with radius of %g\n", cpos[0], cpos[1], cpos[2], radius);
      } else {
        printf("\nAtoms inside sphere centered at [%g %g %g] with radius of %g\n", cpos[0], cpos[1], cpos[2], radius);
      }
      printf(" will be used to create the substitutional solid solution.\n");
      radius *= radius;

      for (int i=0; i<natom; i++){
        double r2 = 0.;
        for (int idim=0; idim<3; idim++){
          double dx = atpos[i][idim]-cpos[idim];
          r2 += dx*dx;
        }
        if (logand){
          if (inside==1 && r2>radius) atsel[i] = 0;
          if (inside==0 && r2<radius) atsel[i] = 0;
        } else {
          if (inside==1 && r2<radius) atsel[i] = 1;
          if (inside==0 && r2>radius) atsel[i] = 1;
        }
      }
    }
  }
  nsel = 0;
  int nsel_type[ntype];
  for (int i=0; i<ntype; i++) nsel_type[i] = 0;

  for (int i=0; i<natom; i++){
    nsel += atsel[i];
    int ip = lookup(attyp[i]);
    nsel_type[ip] += atsel[i];
  }

  printf("\nTotal number of atoms in the system: %d\n", natom);
  printf("Total number of atoms in selection : %d\n", nsel);
  printf("Total number of atomimc types      : %d\n", ntype);
  printf("Atomic type number for each type   :");
  for (int i=0; i<ntype; i++) printf(" %d", typeID[i]);
  printf("\nNumber of atoms for each  type     :");
  for (int i=0; i<ntype; i++) printf(" %d", numtype[i]); printf("\n");
  printf("\nNumber of atoms in selection for each type:");
  for (int i=0; i<ntype; i++) printf(" %d", nsel_type[i]); printf("\n");

  int ipsrc, idsrc=-1, numsub;
  printf("\nPlease input the atomic type to be substituted: ");
  while (count_words(fgets(str,MAXLINE,stdin)) < 1);
  ipsrc = atoi(strtok(str, " \t\n\r\f"));
  idsrc = lookup(ipsrc);
  if (idsrc < 0){
    printf("\nInput atomic type not found, operation terminated!\n");
    for (int i=0; i<14; i++) printf("====="); printf("\n");
    return;
  }

  printf("Total # of atoms with type %d in selection is %d.\n", ipsrc, nsel_type[idsrc]);
  if (nsel_type[idsrc] < 1){
    printf("Not enough atoms to create substitutional solid solution.\n");
    for (int i=0; i<14; i++) printf("====="); printf("\n");
    return;
  }

  double frac = -1.;
  printf("Please input the fraction or total # of atoms to be replaced: ");
  while (count_words(fgets(str,MAXLINE,stdin)) < 1);
  frac = atof(strtok(str, " \t\n\r\f"));
  if (frac < 0. || int(frac) > nsel_type[idsrc]){
    printf("Not enough atoms to create substitutional solid solution.\n");
    for (int i=0; i<14; i++) printf("====="); printf("\n");
    return;
  }

  if (frac <=1.) numsub = MIN(int(frac*double(nsel_type[idsrc])+0.5), nsel_type[idsrc]);
  else numsub = int(frac);
  printf("There will be %d atoms with type %d to be replaced.\n", numsub, ipsrc);
  //if (numsub < 1) return;

  int ipdes, iddes = 1;
  do {
    printf("Please input the atomic type to be assigned: ");
    while (count_words(fgets(str,MAXLINE,stdin)) < 1);
    ipdes = atoi(strtok(str, " \t\n\r\f"));
    iddes = lookup(ipdes);
    if (iddes >= 0){
      iddes = -1;
      printf("***Note: assigned type already exist, continue? (y/n)[y]: ");
      if (count_words(fgets(str,MAXLINE,stdin)) > 0){
        char *ptr = strtok(str, " \t\n\r\f");
        if (strcmp(ptr,"n")==0 || strcmp(ptr,"N")==0) iddes = 1;
      }
    }
  } while (iddes > 0);

  // create random number generator if not created; current time as seed;
  if (random == NULL){
    time_t ctime;
    time(&ctime);
    random = new RanPark((int) ctime);
  }

  // now to generate the solid solution
  int isub =0;
  while (isub < numsub){
    int i = int(random->uniform()*double(natom));
    if (atsel[i] == 0) continue;

    if (attyp[i] == ipsrc){
      attyp[i] = ipdes;
      isub++;
    }
  }

  // reset type info
  typescan();

  // atomic type info after replacement
  int newsel_type[ntype];
  for (int i=0; i<ntype; i++) newsel_type[i] = 0;

  for (int i=0; i<natom; i++){
    nsel += atsel[i];
    int ip = lookup(attyp[i]);
    newsel_type[ip] += atsel[i];
  }
  printf("\nSystem info after creation of solid solution:\n");
  printf("  Total number of atoms in the system: %d\n", natom);
  printf("  Total number of atoms in selection : %d\n", nsel);
  printf("  Total number of atomimc types      : %d\n", ntype);
  printf("  Atomic type number for each type   :");
  for (int i=0; i<ntype; i++) printf(" %d", typeID[i]);
  printf("\n  Number of atoms for each  type     :");
  for (int i=0; i<ntype; i++) printf(" %d", numtype[i]); printf("\n");
  printf("\n  Number of atoms in selection for each type:");
  for (int i=0; i<ntype; i++) printf(" %d", newsel_type[i]); printf("\n");

  for (int i=0; i<14; i++) printf("====="); printf("\n");
return;
}

/*------------------------------------------------------------------------------
 * Method to form multilayers
 *----------------------------------------------------------------------------*/
void Driver::FormLayers()
{
  int nlat = 1, idum;
  int *mynx, *myny;
  char str[MAXLINE];

  printf("\n\n>>>>>>======  To form multilayers with multiple lattices  ======<<<<<<\n");
  printf("NOTE: The 3rd axis of these lattices must be perpendicular to the other 2!\n");
  printf("\nPlease input the number of lattices in your multi-layer system: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) nlat = atoi(strtok(str," \t\n\r\f"));
  if (nlat < 1) return;

  lattice *latts[nlat];
  for (int i=0; i<nlat; i++){
    if (! ShowMenu(i+1) ){nlat = i; break;}

    latts[i] = latt;
    latt = NULL;
  }

  for (int i=0; i<nlat; i++){
    printf("\n>>>>>>   Lattice info for lattice: %c - %s    <<<<<", 'A'+i, latts[i]->name);
    if (latts[i]->perp_z == 0)
      printf("\nWARNING: A3 is not perpendicular to A1 and A2, this lattice cannot be used to form layers!\n");

    latts[i]->display();
  }

  idum = 0;
  printf("\nYou have defined %d lattices: ", nlat);
  for (int i=0; i<nlat; i++) printf(" %c = %s;", 'A'+i, latts[i]->name); printf("\n");
  
  mynx = memory->create(mynx, nlat, "mynx");
  myny = memory->create(myny, nlat, "myny");
  for (int ilat=0; ilat<nlat; ilat++){
    printf("Please input the lateral extensions (nx & ny) for lattice %c: ", 'A'+ilat);
    while (1){
      if ( count_words(fgets(str,MAXLINE,stdin)) == 2 ){
        mynx[ilat] = atoi(strtok(str, " \t\n\r\f"));
        myny[ilat] = atoi(strtok(NULL," \t\n\r\f"));
       if (mynx[ilat] > 0 && myny[ilat] > 0) break;
      }
    }
  }

  nx = ny = nz = 0;
  for (int j=0; j<3; j++){
    latvec[0][j] = latvec[1][j] = 0.;
  }

  printf("\nThe surface vectors for each lattice will be:\n");
  for (int i=0; i<nlat; i++){
    printf("  %c: [", i+'A');
    for (int j=0; j<3; j++){
      double xi = mynx[i]*latts[i]->latvec[0][j]*latts[i]->alat;
      printf("%g ", xi);
      latvec[0][j] += xi;
    }
    printf("] [");
    for (int j=0; j<3; j++){
      double yi = myny[i]*latts[i]->latvec[1][j]*latts[i]->alat;
      printf("%g ", yi);
      latvec[1][j] += yi;
    }
    printf("]\n");
  }
  for (int j=0; j<3; j++){
    latvec[0][j] /= double(nlat);
    latvec[1][j] /= double(nlat);
  }
  printf("Please input your desired surface vectors [%g %g %g, %g %g %g]: ",
    latvec[0][0], latvec[0][1], latvec[0][2], latvec[1][0], latvec[1][1], latvec[1][2]);
  if ( count_words(fgets(str,MAXLINE,stdin)) == 6 ){
    char *ptr = strtok(str," \n\r\t\f");
    for (int i=0; i<2; i++)
    for (int j=0; j<3; j++){ latvec[i][j] = atof(ptr); ptr = strtok(NULL, " \n\r\t\f");}
  }

  double lx, ly, lx0=0., ly0=0.;
  for (int j=0; j<3; j++){
    lx0 += latvec[0][j]*latvec[0][j];
    ly0 += latvec[1][j]*latvec[1][j];
  }
  lx0 = sqrt(lx0); ly0 = sqrt(ly0);

  for (int i=0; i<nlat; i++){
    lx = ly = 0.;
    for (int j=0; j<3; j++){
      double xi = mynx[i]*latts[i]->latvec[0][j]*latts[i]->alat;
      double yi = myny[i]*latts[i]->latvec[1][j]*latts[i]->alat;
      lx += xi*xi;
      ly += yi*yi;
    }
    lx = sqrt(lx); ly = sqrt(ly);
    printf("Lateral misfit for lattice %c is: [%lg %lg]\n", 'A'+i, (lx-lx0)/lx0, (ly-ly0)/ly0);
  }

  double H = 0.;
  int zprev[nlat], ntprev[nlat], zflag = 0, iatom = 0;
  for (int i=0; i<nlat; i++) zprev[i] = 0;

  ntprev[0] = 0;
  for (int i=1; i<nlat; i++) ntprev[i] = ntprev[i-1] + latts[i-1]->ntype;

  char realized[MAXLINE]; strcpy(realized, "");

  // prepare for map info if only one lattice is selected, suggesting layered structure
  int fmap = 0;
  if (nlat == 1) fmap = 1;

  int istr = 0; // start layer ID
  int first = 1, flag_no_interlayer = 0;
  double Hlast = 0., Hfirst = 0., Hextra;
  double shift[2];
  shift[0] = shift[1] = 0.; // define shift of a lattice, in unit of basis vectors

  printf("\nPlease input the layer sequences, for example, if you have two lattices: A and B,\n");
  printf("and you want to have 4 layers A, 5 layers B and then 3 layers A, input A4 B5 A3.\n");
  printf("If extra distance between different lattices is needed, just insert a number between\n");
  printf("them, for example: A4 0.5 B5 0.4 A4 B5...; multiple numbers will add multiple distances.\n");
  printf("If you want to form the 2nd A layers from its first layer in the unit cell, use\n");
  printf("lower case 'a' instead of 'A'. Now, input your sequences: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) {
    char *ptr; if (ptr = strchr(str,'#')) *ptr = '\0';
    ptr = strtok(str," \n\r\t\f");
    while (ptr != NULL){
      zflag = 0; Hextra = 0.;
      int ilat = ptr[0] - 'A';
      if (ilat > nlat){ ilat = ptr[0] - 'a'; zflag = 1;}
      if (ilat >=0 && ilat < nlat){
        latt = latts[ilat];

        strcat(realized," ");strcat(realized, ptr);
        if (zflag) zprev[ilat] = 0;
        zprev[ilat] += istr;

        ptr[0] = ' ';
        int nl_new  = atoi(ptr);
        int ntm_new = 0;
        for (int i=0; i<nl_new; i++) ntm_new += latt->numlayer[(i+zprev[ilat])%latt->nlayer];

        double Hbelow = latt->h[(latt->nlayer-1+zprev[ilat])%latt->nlayer];
        if (first){
          first = 0;
          Hlast = 0.;
          Hfirst = Hbelow;
          H = -Hfirst;
        }
        if (nl_new > 0){
          if (flag_no_interlayer) flag_no_interlayer = 0;
          else H += MAX(Hlast,Hbelow);
        }

        ntm_new *= (mynx[ilat]*myny[ilat]);
        natom += ntm_new;
        atpos = memory->grow(atpos, natom, 3, "atpos");
        attyp = memory->grow(attyp, natom, "attyp");
        if (fmap){
          xmap = memory->grow(xmap, natom, "xmap");
          ymap = memory->grow(ymap, natom, "ymap");
          zmap = memory->grow(zmap, natom, "zmap");
          umap = memory->grow(umap, natom, "umap");

          nx = mynx[ilat];
          ny = myny[ilat];
          nucell = 0;
        }

        for (int k=0; k<nl_new; k++){
          int il = (k+zprev[ilat])%latt->nlayer;
          for (int ia=0; ia<latt->nucell; ia++){
            if ( latt->layer[ia] == il ){
              for (int i=0; i<mynx[ilat]; i++)
              for (int j=0; j<myny[ilat]; j++){
                atpos[iatom][0] = (double(i)+latt->atpos[ia][0]+shift[0])/double(mynx[ilat]);
                atpos[iatom][1] = (double(j)+latt->atpos[ia][1]+shift[1])/double(myny[ilat]);
                atpos[iatom][2] = H;

                if (fmap){
                  xmap[iatom] = i;
                  ymap[iatom] = j;
                  zmap[iatom] = 0;
                  umap[iatom] = nucell;
                }

                attyp[iatom++]  = latt->attyp[ia] + ntprev[ilat];

              }
              nucell++;
            }
          }
          nz++;

          Hlast = latt->h[il];
          H += Hlast;
        }
        if (nl_new > 0) H -= Hlast;
        zprev[ilat] += nl_new%latt->nlayer;
        shift[0] = shift[1] = 0.;
        istr = 0;

      } else if (strcmp(ptr, "-s") == 0){ // to define the start layer of each lattice, must be defined before this lattice
        ptr = strtok(NULL, " \n\r\t\f");
        if (ptr){
          istr = atoi(ptr);
          strcat(realized," -s ");strcat(realized, ptr);
        }

      } else if (strcmp(ptr, "-S") == 0){ // to define the shift in xy direction of the lattice, in unit of unit cell basis vectors, must be defined before the lattice
        char *s0 = strtok(NULL, " \n\r\t\f");
        char *s1 = strtok(NULL, " \n\r\t\f");
        if (s0 && s1){
          shift[0] = atof(s0);
          shift[1] = atof(s1);

          strcat(realized," -S ");
          strcat(realized, s0);
          strcat(realized, s1);
        }

      } else if (strcmp(ptr, "-z") == 0){
        flag_no_interlayer = 1;
        strcat(realized," ");strcat(realized, ptr);

      } else {
        Hextra = atof(ptr);
        strcat(realized," ");strcat(realized, ptr);
        H += Hextra;
      }

      ptr = strtok(NULL, " \n\r\t\f");
    }
  }

  printf("\nThe layer sequences realized is: %s\n", realized);
  printf("In total, %d layers and %d atoms are created.\n", nz, iatom);
  if (fmap) nz = 1;

  if (flag_no_interlayer == 0) H += MAX(Hlast,Hfirst);
  latt = NULL; alat = 1.;
  latvec[2][0] = latvec[2][1] = 0.; latvec[2][2] = H;

  strcpy(str,"Multilayer: ");
  for (int i=0; i<nlat; i++){
    char info[MAXLINE];
    sprintf(info, "%dx%d-%s ", mynx[i], myny[i], latts[i]->name);
    strcat(str,info);
  }
  name = memory->create(name, strlen(str)+1, "driver->FormLayers"); strcpy(name, str);

  double tmp[2];
  for (int i=0; i<natom; i++){
    for (int idim=0; idim<2; idim++) tmp[idim] = atpos[i][idim];
    atpos[i][0] = tmp[0]*latvec[0][0] + tmp[1]*latvec[1][0];
    atpos[i][1] = tmp[0]*latvec[0][1] + tmp[1]*latvec[1][1];
  }

  for (int i=0; i<nlat; i++) delete latts[i];

  memory->destroy(mynx);
  memory->destroy(myny);
  // find the total # of types and # of atoms for each type
  typescan();

return;
}

/*------------------------------------------------------------------------------
 * Method to count # of words in a string, without destroying the string
 *----------------------------------------------------------------------------*/
int Driver::count_words(const char *line)
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
 * Private method to show the version info
 *----------------------------------------------------------------------------*/
void Driver::ShowVersion()
{
  printf("\nLatGen  version 1.%d, compiled on %s.\n", VERSION, __DATE__);
}
/*----------------------------------------------------------------------------*/
