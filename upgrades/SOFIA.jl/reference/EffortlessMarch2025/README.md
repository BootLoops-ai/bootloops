# Effortless

Effortless is a Mathematica package for efficiently constructing odd (algebraic) letters for Feynman integrals with multiple square roots as leading singularities. 
The code takes as inputs the list of rational letters and the list of possible square roots, and it outputs a list of independent algebraic letters. 

The algorithm is described in Section 3.4 of the following publication
- A. Matijašić, *Singularity structure of Feynman integrals with applications to six-particle scattering processes*, [PhD thesis, Munich U., 2024](https://edoc.ub.uni-muenchen.de/34154/)

## Installation and usage

Several functions are implemented using finite field methods implemented in FiniteFlow to enhance efficiency. Functions that use FiniteFlow are characterized by the suffix “FF” at the end of their names. To use those functions, one needs to first install FiniteFlow and all of its dependencies, as outlined in the [FiniteFlow documentation](https://github.com/peraro/finiteflow). 

For all functions using FiniteFlow, there exists a corresponding Mathematica-only implementation. Consequently, users can utilize the package without installing any additional dependencies by importing `Effortless.m` in a Mathematica notebook. 

The usage of the package is demonstrated in [Examples/Examples.nb](Examples/Examples.nb). 


