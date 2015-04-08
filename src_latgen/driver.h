#ifndef DRIVER_H
#define DRIVER_H

#include "lists.h"
#include "random.h"
#include "memory.h"
#include "elements.h"
#include <map>

class Driver{
public:
  Driver();
  ~Driver();
  
  void generate();  // method to generate the atomic configuration
  void modify();    // to modify the resultant model

  void write();     // method to write the atomic configuration (xyz) and mapping info

  void FormLayers();

#ifdef Poly
  void PolyCrystal();
#endif

private:
  int ShowMenu(int);
  void MainMenu();
  void ShowVersion();

  lattice *latt;

  char *name;
  double alat;
  int nx, ny, nz, nucell;         // size in three dimension and # of atoms per unit cell
  int natom, ntype;               // otal number of atoms and atom types
  int *attyp, *numtype, *typeID;  // array to store the atomic types for all
  int *xmap, *ymap, *zmap, *umap; // arrays to store the mapping info
  double **atpos, latvec[3][3];   // arrays to store the atomic positions and lattice info
  void typescan();                // to scan the total number of atomic types in system
  int lookup(int);                // to find the ID of an atomic type

  RanPark *random;                // class object to create random numbers
  Memory *memory;

  // private modification methods
  void solidsol(void);            // method to create subsutitutional solid solution
  void ResetTypeID(void);         // method to reset the atomic type ID
  void MapElement(void);          // method to map atomic type to real elements

  ChemElements * element;
  std::map<int,int> type2num;

  int count_words(const char *);
};

#endif
