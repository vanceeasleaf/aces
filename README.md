[![Build Status](https://travis-ci.org/vanceeasleaf/aces.svg?branch=master)](https://travis-ci.org/vanceeasleaf/aces)

Automatical Computational Experiment System(ACES)
========

@Authour Yang Zhou 
@Fudan University Computational Condensed Matter Group(CCMG)

Docs at [http://vanceeasleaf.github.io/aces/](http://vanceeasleaf.github.io/aces/)

## Introduction

A set of thermal conductivity calculation algorithm including NEMD,Muller-Plathe,Green-Kubo,NEGF,NMA,BTE are realized.

Besides,these calculators are controlled by a universal system which helps you to generate input files and structure files for lammps,VASP,shengbte,almode etc. and use pbs to manage your tasks.Then automatical analysis including data formating and graph generation are carried out.

ACES is a wrapper for many computational codes for thermal conductivity because people want to use all the tools with similar workflows.

A set of thermal conductivity calculation algorithm including NEMD, Muller-Plathe, Green-Kubo, NEGF, NMA, BTE are realized.  Besides, these calculators are controlled by a universal system which helps you to generate input files and structure files for lammps, VASP, ShengBTE, alamode etc. and use pbs to manage your tasks. Then automatically analysis including data formatting and graph generation are carried out.
Meanwhile it’s very easy for you to extend ACES if you have some new codes to integrate. For example, you can write a wrapper for Quantum Expresso easily to use it as an engine to structure optimization, band structure calculation or phonon calculation to have the same work flow with VASP.

