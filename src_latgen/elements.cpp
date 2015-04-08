#include "elements.h"
#include "string.h"

int ChemElements::Name2Num(const char *ename)
{
  int AtNum;
  for (AtNum = 1; AtNum <= NumMax; AtNum++) if (strcmp(ename, symbol[AtNum]) == 0) break;
  if (AtNum > NumMax) AtNum = 0;

return AtNum;
}

void ChemElements::Num2Name(const int AtNum, char * ename)
{
  int num = AtNum;
  if (num < 0||num > NumMax) num = 0;
  strcpy(ename, symbol[num]);
}

double ChemElements::Name2Mass(const char * ename)
{
  int AtNum = Name2Num(ename);

return weight[AtNum];
}

double ChemElements::Num2Mass(const int AtNum)
{
  int num = AtNum;
  if (num < 0||num > NumMax) num = 0;
  return weight[num];
}

/* -----------------------------------------------------------------------------
 * Atomic weight data taken from:
 * Pure Appl. Chem., Vol. 83, No. 2, pp. 359–396, 2011.
 * Atomic weights of the elements 2009 (IUPAC Technical Report)
 * ---------------------------------------------------------------------------*/
const int ChemElements::NumMax = 112;
const double ChemElements::weight[] = {0.,                                  // Null
    1.00797,    4.0026,     6.939,      9.012182,    10.811,                // H - B
   12.01115,   14.0067,    15.9994,    18.9984032,   20.17976,   22.989769, // C - Na
   24.30506,   26.9815386, 28.086,     30.973762,    32.064,                // Mg - S
   35.453,     39.948,     39.0983,    40.078,       44.955912,             // Cl - Sc
   47.867,     50.9415,    51.9961,    54.938045,    55.845,                // Ti - Fe
   58.933195,  58.6934,    63.546,     65.38,        69.723,                // Co - Ga
   72.63,      74.9216,    78.96,      79.904,       83.798,                // Ge - Kr
   85.4678,    87.62,      88.90585,   91.224,       92.90638,              // Rb - Nb
   95.96,      98.9062,   101.07,     102.9055,     106.42,                 // Mo - Pd
  107.8682,   112.411,    114.818,    118.71,       121.76,                 // Ag - Sb
  127.6,      126.90447,  131.293,    132.9054519,  137.327,                // Te - Ba
  138.90547,  140.116,    140.90765,  144.242,      147.,                   // La - Pm
  150.36,     151.964,    157.25,     158.92535,    162.5,                  // Sm - Dy
  164.93032,  167.259,    168.93421,  173.054,      174.9668,               // Ho - Lu
  178.49,     180.94788,  183.84,     186.207,      190.23,                 // Hf - Os
  192.217,    195.084,    196.966569, 200.59,       204.383,                // Ir - Tl
  207.2,      208.9804,   209.,       210.,         222.,                   // Pb - Rn
  223.,       226.025,    227.028,    232.03806,    231.03588,              // Fr - Pa
  238.02891,  237.048,    244.,       243.,         247.,                   // U - Cm
  247.,       251.,       252.,       257.,         258.,                   // Bk - Md
  259.,       260.,       261.11,     262.11,       263.12,                 // No - Sg
  262.12,     265.,       266.,       269.,         272.,      285.         // Bh - Cn
};

const char ChemElements::symbol[][3] = { "X",
    "H",  "He", "Li", "Be", "B",
    "C",  "N",  "O",  "F",  "Ne", "Na",
    "Mg", "Al", "Si", "P",  "S",
    "Cl", "Ar", "K",  "Ca", "Sc",
    "Ti", "V",  "Cr", "Mn", "Fe",
    "Co", "Ni", "Cu", "Zn", "Ga",
    "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y",  "Zr", "Nb",
    "Mo", "Tc", "Ru", "Rh", "Pd",
    "Ag", "Cd", "In", "Sn", "Sb",
    "Te", "I",  "Xe", "Cs", "Ba",
    "La", "Ce", "Pr", "Nd", "Pm",
    "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W",  "Re", "Os",
    "Ir", "Pt", "Au", "Hg", "Tl",
    "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa",
    "U",  "Np", "Pu", "Am", "Cm",
    "Bk", "Cf", "Es", "Fm", "Md",
    "No", "Lr", "Rf", "Db", "Sg",
    "Bh", "Hs", "Mt", "Ds", "Rg", "Cn"
};
