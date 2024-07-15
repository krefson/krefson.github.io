---
author:
- Keith Refson, STFC Rutherford Appleton Laboratory
bibliography:
- dfpt.bib
- dft.bib
title: Phonon and related calculations using CASTEP
---

Copyright © Science and Technology Facilities Council, 2009

### Introduction

This guide introduces the concepts, keywords and techniques needed to
set up and run CASTEP calculations of lattice dynamics and vibrational
spectroscopy. The material covers the various lattice-dynamical and
related methods implemented in CASTEP, how to set up a calculation, and
presents simple examples of the major types of calculation. It is
assumed that the reader is familiar with the general CASTEP input and
output files and keywords to the level of, for example the tutorials on
geometry optimisation which may be found at
[www.castep.org](http://www.castep.org). It describes how to analyse the
results and generate graphical output using the CASTEP tools, but does
not cover the modelling of IR, Raman or inelastic neutron or X-ray
spectra, which are large subjects and beyond the scope of this document.

It does not describe the compilation and installation of CASTEP and its
tools, nor does it describe the operational details of invoking and
running CASTEP. Computational clusters, HPC computing environments and
batch systems vary to a considerable degree, and the reader is referred
to their cluster or computer centre documentation. It does discuss
general aspects of managing and running large phonon calculations
including restarts and parallelism in
section <a href="#sec:large-calculations" data-reference-type="ref"
data-reference="sec:large-calculations">6</a>.

There are many useful textbooks on the theory of vibrations in crystals,
or *lattice dynamics* as the subject is usually known. A good beginner’s
guide is “Introduction to Lattice Dynamics” by Martin Dove (Dove 1993).
More advanced texts are available by Srivistava (Srivastava 1990),
Maradudin (Maradudin and Horton 1980) and many others. The following
section presents only a brief summary to introduce the notation.

#### Theory of Lattice Dynamics

Consider a crystal with a unit cell containing $N$ atoms, labelled
$\kappa$, and $a$ labels the primitive cells in the lattice. The crystal
is initially in mechanical equilibrium, with Cartesian co-ordinates
$R_{{\kappa,\alpha}}$, ($\alpha=1..3$ denotes the Cartesian $x$,$y$ or
$z$ direction).
${{\mathbf{u}}_{{\kappa,\alpha},a}}= x_{{\kappa,\alpha,a}} - R_{{\kappa,\alpha,a}}$
denotes the displacement of an atom from its equilibrium position.
Harmonic lattice dynamics is based on a Taylor expansion of total energy
about structural equilibrium co-ordinates.

$$\label{eq:taylor}
     E = E_0 + \sum_{{\kappa,\alpha,a}} \frac{\partial E}{\partial{{\mathbf{u}}_{{\kappa,\alpha},a}}}.{{\mathbf{u}}_{{\kappa,\alpha},a}}+ {\frac{1}{2}}
    \sum_{{\kappa,\alpha,a},{\kappa^\prime,\alpha^\prime},a^\prime}  {{\mathbf{u}}_{{\kappa,\alpha},a}}.{\Phi^{\kappa,\kappa^\prime}_{\alpha,\alpha^\prime}(a,a^\prime)}. {{\mathbf{u}}_{{\kappa^\prime,\alpha^\prime},a^\prime}}+ ...$$

where ${{\mathbf{u}}_{{\kappa,\alpha},a}}$ is the vector of atomic
displacements from equilibrium and
${\Phi^{\kappa,\kappa^\prime}_{\alpha,\alpha^\prime}(a)}$ is the matrix
of *force constants*

$$\label{eq:fcmat}
 {\Phi^{\kappa,\kappa^\prime}_{\alpha,\alpha^\prime}(a)}\equiv {\Phi^{\kappa,\kappa^\prime}_{\alpha,\alpha^\prime}(a,0)}= \frac{\partial^2 E}{\partial {{\mathbf{u}}_{{\kappa,\alpha},a}}\partial {{\mathbf{u}}_{{\kappa^\prime,\alpha^\prime},0}}} \;.$$

At equilibrium the forces
$F_{{\kappa,\alpha}} = -\frac{\partial E}{\partial{{\mathbf{u}}_{{\kappa,\alpha}}}}$
are all zero so the first-order term vanishes. In the *Harmonic
Approximation* the 3$^{\text{rd}}$ and higher order terms are neglected
[^1]. Assume *Born von Karman* periodic boundary conditions and
substitute a plane-wave guess for the solution

$${{\mathbf{u}}_{{\kappa,\alpha}}}= {\mathbf{\varepsilon}_{m{\kappa,\alpha}{\mathbf{q}}}}\exp( i {\mathbf{q}}. {\mathbf{R}}_{{\kappa,\alpha}} - \omega_{m} t)$$

with a phonon wavevector ${\mathbf{q}}$ and a *polarization vector*
${\mathbf{\varepsilon}_{m{\kappa,\alpha}{\mathbf{q}}}}$. This yields an
eigenvalue equation

$$\label{eq:dmat-eigen}
{D^{\kappa,\kappa^\prime}_{\alpha,\alpha^\prime}(\mathbf{q})}{\mathbf{\varepsilon}_{m{\kappa,\alpha}{\mathbf{q}}}}= \omega_{m,{\mathbf{q}}}^2 {\mathbf{\varepsilon}_{m{\kappa,\alpha}{\mathbf{q}}}}\; .$$

The **dynamical matrix** is defined as

$$\label{eq:dmat}
{D^{\kappa,\kappa^\prime}_{\alpha,\alpha^\prime}(\mathbf{q})}= \frac{1}{\sqrt{M_\kappa M_{\kappa^\prime}}} {C^{\kappa,\kappa^\prime}_{\alpha,\alpha^\prime}(\mathbf{q})}=
\frac{1}{\sqrt{M_\kappa M_{\kappa^\prime}}} \sum_{a} {\Phi^{\kappa,\kappa^\prime}_{\alpha,\alpha^\prime}(a)}
    e^{-i \mathbf{q}.{\mathbf{r}}_a}$$

where the index $a$ runs over all lattice images of the primitive cell.
That is, the dynamical matrix is the mass-reduced Fourier transform of
the force constant matrix.

The eigenvalue equation
<a href="#eq:dmat-eigen" data-reference-type="ref"
data-reference="eq:dmat-eigen">[eq:dmat-eigen]</a> can be solved by
standard numerical methods.
${D^{\kappa,\kappa^\prime}_{\alpha,\alpha^\prime}(\mathbf{q})}$ is
Hermitian by construction. The vibrational frequencies at each mode are
obtained as the square roots of the eigenvalues, and the eigenvectors
give the pattern of atomic displacements belonging to each mode.

The central question of *ab-initio* lattice dynamics is therefore how to
determine the force constants
${\Phi^{\kappa,\kappa^\prime}_{\alpha,\alpha^\prime}}$ which are the
*second derivatives* of the total energy with respect to two atomic
displacements. The *first* derivatives are the forces, and may be
determined straightforwardly using the Hellman–Feynmann Theorem as the
derivative of the quantum mechanical energy with respect to an atomic
displacement, $\lambda$

$$\begin{aligned}
 E &=  {\left < \psi \right |}{\hat{H}}{\left | \psi \right >}\qquad \text{with}\qquad {\hat{H}}= \nabla^2 + V_{\text{SCF}}\\
F  &=  - \frac{dE}{d\lambda} \\
   &=  - {\left < {\frac{d \psi}{d \lambda}}\right |}{\hat{H}}{\left | \psi \right >}-
    {\left < \psi \right |}{\hat{H}}{\left | {\frac{d \psi}{d \lambda}}\right >}- {\left < \psi \right |}{\frac{d V}{d \lambda}}{\left | \psi \right >}\label{eq:hellman} \; .
\end{aligned}$$

If ${\left < \psi \right |}$ represents the ground state of ${\hat{H}}$
then the first two terms vanish because
${\left < \psi \right |}{\hat{H}}{\left | {\frac{d \psi}{d \lambda}}\right >}= \epsilon_n {\left < \psi \right |}{\left | {\frac{d \psi}{d \lambda}}\right >}= 0$.
Differentiating again yields

$$\frac{d^2 E}{d \lambda^2} = {\left < {\frac{d \psi}{d \lambda}}\right |}{\frac{d V}{d \lambda}}{\left | \psi \right >}+
    {\left < \psi \right |}{\frac{d V}{d \lambda}}{\left | {\frac{d \psi}{d \lambda}}\right >}- {\left < \psi \right |}{\frac{d^2 V}{d \lambda^2}}{\left | \psi \right >}\;.$$

Unlike equation <a href="#eq:hellman" data-reference-type="ref"
data-reference="eq:hellman">[eq:hellman]</a> the terms involving the
derivatives of the wavefunctions do not vanish. This means that it is
necessary to somehow compute the *electronic response* of the system to
the displacement of an atom to perform *ab-initio* lattice dynamics
calculations. This may be accomplished either by a finite-displacement
method, *i.e.* performing two calculations which differ by a small
displacement of an atomic co-ordinate and evaluating a numerical
derivative, or by using perturbation theory to evaluate the response
wavefunction ${\frac{d \psi}{d \lambda}}$. (In Kohn-Sham DFT, as
implemented in CASTEP we actually evaluate the first-order response
orbitals.)

A good general reference on *ab-initio* lattice dynamics methods is the
review paper by Baroni *et al.*(Baroni et al. 2001).

#### Prerequisites

##### Geometry Optimisation

A successful phonon calculation almost always requires a preceding
geometry optimisation (except for small, high symmetry system where all
atoms lie on crystallographic high-symmetry positions and the forces are
zero by symmetry). It is *not* necessary to perform a variable-cell
optimisation - the lattice dynamics is well defined at any stress or
pressure, and phonons in high-pressure or strained systems are
frequently of scientific interest. The two most convenient ways of
achieving this are to

- Set up the phonon run as a continuation of the geometry optimisation
  by setting the parameters keyword `continuation : <geom-seed.check>`.

- Add the keyword `write_cell_structure : TRUE` to the geometry
  optimization run and modify the resulting \<seedname\>-out.cell to use
  as the input for a new run.

The importance of a high-quality structure optimisation can not be
overemphasized - the energy expansion in
equation <a href="#eq:taylor" data-reference-type="ref"
data-reference="eq:taylor">[eq:taylor]</a> makes the explicit assumption
that the system is in mechanical equilibrium and that all atomic forces
are zero.

This is sufficiently important that before performing a phonon
calculation CASTEP will compute the residual forces to determine if the
geometry is converged. If any component of the force exceeds
geom_force_tol it will print an error message and abort the run. Should
a run fail with this message, it may be because the geometry
optimisation run did not in fact succeed, or because some parameter
governing the convergence (*e.g.* the cutoff energy) differs in the
phonon run compared to the geometry optimisation. In that case the
correct procedure would be to re-optimise the geometry using the same
parameters as needed for the phonon run. Alternatively if the geometry
error and size of the force residual are tolerable, then the value of
geom_force_tol may be increased in the .param file of the phonon
calculation which will allow the run to proceed.

How accurate a geometry optimization is needed? Accumulated practical
experience suggests that substantially tighter tolerances are needed to
generate reasonable quality phonons than are needed for structural or
energetics calculations. For many crystalline systems a geometry force
convergence tolerance set using parameter geom_force_tol of 0.01 eV/Å is
typically needed. For “soft” materials containing weak bonds such as
molecular crystals or in the presence of hydrogen bonds, an even smaller
value is frequently necessary. Only careful convergence testing of the
geometry and resulting frequencies can determine the value to use. To
achieve a high level of force convergence, it is obviously essential
that the forces be evaluated to at least the same precision. This will
in turn govern the choice of electronic k-point sampling, and probably
require a smaller than default SCF convergence tolerance,
elec_energy_tol. See
section <a href="#sec:convergence" data-reference-type="ref"
data-reference="sec:convergence">2.4</a> for further discussion of
geometry optimisation and convergence for phonon calculations.

If a lattice dynamics calculation is performed at the configuration
which minimises the energy the force constant matrix
${\Phi^{\kappa,\kappa^\prime}_{\alpha,\alpha^\prime}}$ is positive
definite, and all of its eigenvalues are positive. Consequently the
vibrational frequencies which are the square roots of the eigenvalues
are real numbers. If on the other hand the system is not at a minimum
energy equilibrium configuration
${\Phi^{\kappa,\kappa^\prime}_{\alpha,\alpha^\prime}}$ is not
necessarily positive definite, the eigenvalues may be negative. In that
case the frequencies are imaginary and do not correspond to a physical
vibrational mode. For convenience CASTEP prints such values using a a
minus sign “$-$”. These are sometimes loosely referred to as “negative”
frequencies, but bear in mind that this convention really indicates
negative *eigenvalues* and *imaginary* frequencies. [^2]

##### Use of Symmetry

CASTEP lattice dynamics uses symmetry to reduce the number of
perturbation calculations to the irreducible set, and applies the
symmetry transformations to generate the full dynamical matrix. Combined
with the usual reduction in k-points in the IBZ, savings in CPU time by
use of symmetry can easily be a factor of 10 or more. It is therefore
important to maximise the symmetry available to the calculation, and to
either specify the symmetry operations in the .cell file using the
`%block SYMMETRY_OPS` keyword or generate them using keyword
symmetry_generate.

CASTEP generates the dynamical matrix using perturbations along
Cartesian X, Y and Z directions. To maximise the use of symmetry and
minimise the number of perturbations required, the symmetry axes of the
crystal or molecule should be aligned along Cartesian X, Y, and Z. If
the cell axes are specified using `%block LATTICE_ABC`, CASTEP attempts
to optimally orient the cell in many cases, but if using
`%block LATTICE_CART` it is the responsibility of the user. This is one
of the few cases where the absolute orientation makes any difference in
a CASTEP calculation. An optimal orientation will use the fewest
perturbations and lowest CPU time and will also exhibit best convergence
with respect to, for example electronic k-point set.

There is a further consideration related to symmetry; structure, cell
and atomic co-ordinates must be specified to a reasonably high accuracy
in the input files. CASTEP uses stochastic methods to analyse the effect
of symmetry operations on the dynamical matrix, and this analysis can
fail to detect symmetries or otherwise misbehave unless symmetry
operations, lattice vectors and atomic co-ordinates are consistent to a
reasonable degree of precision. It is recommended that symmetry-related
lattice vectors and internal co-ordinates are consistent to at least
10$^{-6}$Å  which can only be achieved if the values in the .cell file
are expressed to this number of decimal places. This is particularly
significant in the case of trigonal or hexagonal systems specified with
`%block LATTICE_CART` where equality of the $a$ and $b$ cell vector
lengths is only as precise as the number of significant figures used to
represent the components.

Two features of CASTEP’s .cell file input may be helpful. First is the
ability to input cell vectors and atomic positions using (a limited set
of) mathematical syntax. See
figure <a href="#fig:example-gamma" data-reference-type="ref"
data-reference="fig:example-gamma">2.1</a> for an example. Second, if
the cell file keyword snap_to_symmetry is present, CASTEP will adjust
co-ordinates and vectors to satisfy symmetry to high precision.

Alternatively there is a utility program in the academic distribution
named symmetry_snap which implements the same functionality. This is
invoked simply as

> symmetry_snap *seedname*

which reads *seedname*.cell and outputs the symmetrised version to
*seedname*-symm.cell

### Running phonon calculations

CASTEP contains implementations of several methods for computing the
electronic response orbitals and the dynamical matrix (Refson, Tulip,
and Clark 2006). There are implementation restrictions and the choice of
the most suitable one depends on the type of calculation, the
Hamiltonian, and the property calculatsions required. See
table <a href="#tbl:captable" data-reference-type="ref"
data-reference="tbl:captable">2.1</a> for details. For straightforward
semi-local DFT calculations (LDA, PBE etc.) the *density-functional
perturbation theory method* (DFPT) method is preferred
(section <a href="#sec:dfpt-gamma" data-reference-type="ref"
data-reference="sec:dfpt-gamma">2.1</a>), as this is not only the most
efficient, but also allows the calculation of infra red and Raman
intensities for modelling of spectra. DFPT is not yet implemented for
ultrasoft pseudopotentials, or for Hubbard U or some
dispersion-corrected DFT methods (as of release 24.1), and in these
cases the *finite displacement* method may be used
(section <a href="#sec:fd" data-reference-type="ref"
data-reference="sec:fd">2.3</a>). If a density of states or finely
sampled dispersion curve along high symmetry directions is needed, then
either the DFPT method with *Fourier interpolation*
(section <a href="#sec:ddos" data-reference-type="ref"
data-reference="sec:ddos">2.2</a>) is the most suitable or the
*finite-displacement supercell*
(section <a href="#sec:supercell" data-reference-type="ref"
data-reference="sec:supercell">2.3.4</a>; if ultrasoft pseudopotentials
are required). A summary of the recommended approach for various
property calculations is given in
table <a href="#tbl:method-selection" data-reference-type="ref"
data-reference="tbl:method-selection">2.2</a>.

<div id="tbl:captable">

|                         | DFPT (Phonon) | DFPT (E-field) | FD (Phonon) |
|:------------------------|:-------------:|:--------------:|:-----------:|
| USP                     |       ✘       |       ✘        |      ✓      |
| NCP                     |       ✓       |       ✓        |      ✓      |
| LDA, GGA                |       ✓       |       ✓        |      ✓      |
| MGGA                    |       ✘       |       ✘        |      ✓      |
| DFT+U                   |       ✘       |       ✘        |      ✓      |
| NCM/SOC                 |       ✘       |       ✘        |      ✓      |
| PBE0, Hybrid XC         |       ✘       |       ✘        |      ✓      |
| DFT+D: TS, D2           |       ✓       |       ✓        |      ✓      |
| DFT+D: D3,D4, MBD\*,XDM |       ✘       |       ✓        |      ✓      |

CASTEP phonons implementation matrix. In general DFPT is avalable for
semilocal DFPT, but not other Hamiltonians. The combination of DFPT and
USPs is not implemented

</div>

<div class="threeparttable">

<div id="tbl:method-selection">

| Target Property                                 | Preferred method                                       |
|:------------------------------------------------|:-------------------------------------------------------|
| IR spectrum                                     | DFPT at q=0 with NCP potentials                        |
|                                                 | FD at q=0 with NCP potentials and DFPT e-field         |
|                                                 | FD at q=0 with USP potentials and Berry Phase $Z^{*}$  |
| Raman spectrum                                  | DFPT at q=0 with NCP potentials (2n+1 theorem)         |
| Born Effective Charges ($Z^{*}$)                | DFPT E-field using NCP potentials                      |
|                                                 | FD at q=0 with USP potentials and Berry Phase $Z^{*}$  |
| Dielectric Permittivity ($\epsilon_\infty$)     | DFPT E-field using NCP potentials                      |
| Nonlinear optical susceptibility ($\chi^{(2)}$) | DFPT E-field (2n+1 theorem) using NCP potentials       |
| Phonon dispersion or DOS                        | DFPT plus interpolation with NCP potentials            |
|                                                 | FD plus interpolation using USP or NCP potentials      |
|                                                 | FD with supercell using USP or NCP potentials          |
| Vibrational Thermodynamics                      | same as DOS                                            |

Available and recommended methods for different property calculations.

</div>

</div>

<div class="tablenotes">

a Will be released in version 25.1. See section
<a href="#sec:berry-fd" data-reference-type="ref"
data-reference="sec:berry-fd">2.3.2.1</a>.

</div>

#### A DFPT phonon calculation at the $\Gamma$-point

Many of the principles of setting up and running phonon calculations can
be illustrated in the simplest case - computing phonon frequencies at
the $q=(0,0,0)$, often referred to as the $\Gamma$ point. This is a very
common calculation, as it forms the basis of modelling of infra-red or
Raman spectra and group theoretical analysis and assignment of the
modes.

##### Input files

The setup for a CASTEP phonon calculation requires a few additional
keywords in the *seedname*.cell file. Like any type of calculation, the
unit cell must be specified using either of `%block LATTICE_ABC` or
`%block LATTICE_CART`, and the atomic co-ordinates using
`%block POSITIONS_FRAC` or `%block POSITIONS_ABS`. Figure
<a href="#fig:example-gamma" data-reference-type="ref"
data-reference="fig:example-gamma">2.1</a> shows a complete input file
for a calculation on boron nitride in the hexagonal Wurtzite structure.
The additional keywords. phonon_kpoint_list is used to specify that a
single phonon wavevector (0,0,0) is to be computed[^3]. More wavevectors
could be specified using additional lines in this block.[^4]

<figure id="fig:example-gamma">
<pre><code>%block LATTICE_ABC
   2.536000   2.536000     4.199000  ! a, b, c
  90.000000  90.000000   120.000000  ! alpha, beta, gamma
%endblock LATTICE_ABC

%block POSITIONS_FRAC
  B   2/3   1/3    0.000000  ! Atom co-ordinates in crystallographic
  N   2/3   1/3    0.374536  ! fractional system.
  B   1/3   2/3    0.500000
  N   1/3   2/3    0.874536
%endblock POSITIONS_FRAC

%block SPECIES_POT
B   B_00.recpot              ! File names of pseudopotential to use for B
N   N_00.recpot
%endblock SPECIES_POT

kpoints_mp_spacing 0.07      ! Generate MP grid of electronic k-points.
symmetry_generate            ! Find symmetry operations of crystal structure

%block PHONON_KPOINT_LIST
    0.0   0.0   0.0    1.0   ! Wavevector of phonon(s) to compute ( qx qy qz, weight)
%endblock PHONON_KPOINT_LIST
</code></pre>
<pre><code>task             : PHONON    ! Choose a phonon calculation
xc_functional    : LDA
opt_strategy     : SPEED     ! Optimise for speed over memory saving
cut_off_energy   : 700.0 eV  ! Plane-wave cutoff for this system
fix_occupancy    : TRUE      ! Mandatory for a DFPT phonon calculation
elec_method      : DM        ! Use fast DM solver even for insulating system
phonon_sum_rule  : TRUE      ! Enforce acoustic sum rule on calculated D</code></pre>
<figcaption>Annotated example cell and parameter files for BN in the
Wurtzite structure. Note that <span>phonon_method</span> defaults to
<span>DFPT</span>.</figcaption>
</figure>

<figure>
<pre><code> ==============================================================================
 +                           Vibrational Frequencies                          +
 +                           -----------------------                          +
 +                                                                            +
 + Performing frequency calculation at  1 wavevector  (q-pt )                 +
 + ========================================================================== +
 +                                                                            +
 + -------------------------------------------------------------------------- +
 +  q-pt=    1 (  0.000000  0.000000  0.000000)     1.0000000000              +
 + -------------------------------------------------------------------------- +
 +  Acoustic sum rule correction &lt;   8.094974 cm-1 applied                    +
 +     N      Frequency irrep.    ir intensity active            raman active + 
 +                (cm-1)         ((D/A)**2/amu)                               + 
 +                                                                            +
 +     1      -0.049170   a          0.0000000  N                       N     + 
 +     2      -0.034768   b          0.0000000  N                       N     + 
 +     3      -0.034768   b          0.0000000  N                       N     + 
 +     4     475.116083   c          0.0000000  N                       Y     + 
 +     5     475.116083   c          0.0000000  N                       Y     + 
 +     6     952.000075   c          0.0000000  N                       Y     + 
 +     7     952.000075   c          0.0000000  N                       Y     + 
 +     8     963.032787   d          0.0000000  N                       N     + 
 +     9    1016.312039   a         28.5895083  Y                       Y     + 
 +    10    1051.124801   b         25.8079204  Y                       Y     + 
 +    11    1051.124801   b         25.8079204  Y                       Y     + 
 +    12    1123.472590   d          0.0000000  N                       N     + 
 + .......................................................................... +
 +        Character table from group theory analysis of eigenvectors          +
 +                           Point Group =  25, C6v                           +
 + .......................................................................... +
 +  Rep  Mul |    E   6   3   2 m_v m_v                                       +
 +           | ------------------------                                       +
 + a A1    2 |    1   1   1   1   1   1                                       +
 + b E1    2 |    2   1  -1  -2   0   0                                       +
 + c E2    2 |    2  -1  -1   2   0   0                                       +
 + d B1    2 |    1  -1   1  -1  -1   1                                       +</code></pre>
<figcaption>Output in the .castep file generated by the preceding input
files. The columns show the mode frequencies, a letter labelling of the
irreducible representation of the mode, the infrared absorptivity and
Raman activity, and whether or not the mode is Raman and/or IR active.
Not all of these are present in every calculation depending on the use
of symmetry and input options.</figcaption>
</figure>

<figure id="fig:example-gamma-out">
<pre><code> + -------------------------------------------------------------------------- +
 +  q-pt=    1 (  0.000000  0.000000  0.000000)     0.0000000000              +
 +  q-&gt;0 along (  0.000000  0.000000  1.000000)                               +
 + -------------------------------------------------------------------------- +
 +  Acoustic sum rule correction &lt;   8.094974 cm-1 applied                    +
 +     N      Frequency irrep.    ir intensity active            raman active + 
 +                (cm-1)         ((D/A)**2/amu)                               + 
 +                                                                            +
 +     1      -0.037851   a          0.0000000  N                       N     + 
 +     2      -0.034768   b          0.0000000  N                       N     + 
 +     3      -0.034768   b          0.0000000  N                       N     + 
 +     4     475.116083   c          0.0000000  N                       Y     + 
 +     5     475.116083   c          0.0000000  N                       Y     + 
 +     6     952.000075   c          0.0000000  N                       Y     + 
 +     7     952.000075   c          0.0000000  N                       Y     + 
 +     8     963.032787   d          0.0000000  N                       N     + 
 +     9    1051.124801   b         25.8079204  Y                       Y     + 
 +    10    1051.124801   b         25.8079204  Y                       Y     + 
 +    11    1123.472590   d          0.0000000  N                       N     + 
 +    12    1262.933291   a         28.5895083  Y                       Y     + 
 + .......................................................................... +
 +        Character table from group theory analysis of eigenvectors          +
 +                           Point Group =  25, C6v                           +
 + (Due to LO/TO splitting this character table may not contain some symmetry +
 +  operations of the full crystallographic point group.  Additional          +
 +  representations may be also be present corresponding to split LO modes.   +
 +  A conventional analysis can be generated by specifying an additional null +
 +  (all zero) field direction or one along any unique crystallographic axis  +
 +  in %BLOCK PHONON_GAMMA_DIRECTIONS in &lt;seedname&gt;.cell.)                    +
 + .......................................................................... +
 +  Rep  Mul |    E   6   3   2 m_v m_v                                       +
 +           | ------------------------                                       +
 + a A1    2 |    1   1   1   1   1   1                                       +
 + b E1    2 |    2   1  -1  -2   0   0                                       +
 + c E2    2 |    2  -1  -1   2   0   0                                       +
 + d B1    2 |    1  -1   1  -1  -1   1                                       +
 ==============================================================================</code></pre>
<figcaption>Run output (continued): If the LO-TO-splitting calculation
is active the <span class="math inline">\(\Gamma\)</span> point
frequency table is repeated, excluding and including the non-analytic
contribution which generates the LO/TO splitting in one or more
directions (See Gonze and Lee <span class="citation"
data-cites="GonzeL97">(Gonze and Lee 1997)</span>). The group theory
analysis usually only corresponds to the full crystallographic point
group in the case without symmetry breaking by LO/TO
splitting.</figcaption>
</figure>

##### Reading the output

Figure <a href="#fig:example-gamma-out" data-reference-type="ref"
data-reference="fig:example-gamma-out">2.2</a> shows part of the
phonon-relevant output extracted from the .castep file obtained by
running the input files of the previous section. There are several
blocks of output, one per direction chosen for
${\mathbf{q}}\rightarrow 0$ in the LO-TO splitting terms. Note that
CASTEP has added a calculation without LO-TO splitting (the first
block), even though this was not explicitly requested. Within each block
frequencies are listed one per line. Also on the line are (a) a label
indicating the irreducible representation of the mode from a group
theory analysis, (b) the computed absorptivity of the mode in a powder
(or otherwise orientationally averaged) infrared experiment, (c) whether
the mode is IR active, (d) the Raman activity (if computed) and (e)
whether the mode is Raman active.

As with other calculations the amount of information written to the
.castep file is controlled by the value of the parameter iprint. The
levels of output are

1.  Minimal output as in
    figure <a href="#fig:example-gamma-out" data-reference-type="ref"
    data-reference="fig:example-gamma-out">2.2</a>. No progress info.

2.  As for `iprint : 0` but including a reassuring progress report of
    q-points and perturbations.

3.  More detailed output including details of k-points and symmetry for
    each perturbation, cycle-by cycle DFPT minimiser report, printing of
    dynamical matrices (and force constant matrices).

In addition to the user-readable output in the .castep file[^5], every
phonon calculation generates an additional output file with the suffix
.phonon which is intended for postprocessing analysis by other programs.
This includes not only the frequencies but also the eigenvectors
${\mathbf{\varepsilon}_{m{\kappa,\alpha}{\mathbf{q}}}}$ resulting from
diagonalising the dynamical matrix. These eigenvectors are orthonormal
by construction, and the relationship between the eigenvectors and
atomic displacements is given by

$$u_{\kappa,\alpha,a} = \frac{1}{\sqrt{M_{\kappa}}} Q_{m} {\mathbf{\varepsilon}_{m{\kappa,\alpha}{\mathbf{q}}}}\exp(i ({\mathbf{q}}\cdot {\mathbf{R}}_{\kappa,\alpha}(a) - \omega_m t))$$

where $Q_M$ is the amplitude of mode $m$ and the other notation is as
set out in the introduction.

##### Raman activity calculations

In addition to the infrared absorptivity computed as a default part of a
$\Gamma$-point phonon calculation, CASTEP is also capable of computing
the Raman activity tensors in the case of non-resonant scattering for an
insulating system[^6]. This calculation uses an extension of DFPT known
as the “2n+1 theorem” (Baroni et al. 2001; Miwa 2011) which requires
substantially more computational effort than a bare phonon calculation
and is therefore not enabled by default. To compute Raman activity set

> `calculate_raman : TRUE`

in the .param file, which calculates and prints the Raman atomic
susceptibility tensors and the mode susceptibility tensors to the
.castep file. These may be extracted for use when modelling the Raman
scattering in a polarised single-crystal geometry.[^7]

As with the Born effective charges and dynamical matrices, there is a
sum rule on the values of the atomic polar tensors, which may not be
well satisfied because of numerical approximations. A correction will be
applied if parameter

> `raman_sum_rule : TRUE`

is set.

A 3d-average Raman activity (see Porezag and Pedersen (Porezag and
Pederson 1996)) is computed from the susceptibility tensors and printed
as an additional column in the usual frequency output given in
figure <a href="#fig:example-gamma-out" data-reference-type="ref"
data-reference="fig:example-gamma-out">2.2</a>. This may be used for
simple modelling of a powder spectrum, for example using the dos.pl tool
(see section <a href="#sec:dos-pl" data-reference-type="ref"
data-reference="sec:dos-pl">5.1</a>).

#### Dispersion and density of states

Phonons of a nonzero wavevector play an important role in the
thermophysical properties of crystalline solids and the physics of many
solid-state phase transitions. Proving the mechanical stability of a
crystal structure by testing for real frequencies requires a vibrational
calculation over the full Brillouin Zone. And dispersion curves and
densities of states are frequently required for comparison with
inelastic neutron and X-ray scattering experiments.

One of the benefits of density functional perturbation theory is that
CASTEP can calculate vibrational modes at ${\mathbf{q}}\ne \mathbf{0}$
as easily as at ${\mathbf{q}}= \mathbf{0}$ (but at the cost of an
increase in CPU time due to the decreased symmetry). It is possible to
perform this calculation simply by providing a list of q-points in the
.cell file ( using one of the blocks `%block PHONON_KPOINT_LIST`, or
`%block PHONON_KPOINT_PATH` or keyword phonon_kpoint_mp_grid). However
to generate a reasonable quality dispersion plot or density of states
will usually require hundreds or thousands of q-points. The time
required by such a calculation would be several hundred times that for a
single-point energy and therefore infeasibly large.

Fortunately there is a way to achieve the same result at a far smaller
computational cost. This method exploits the fact that the interatomic
interactions in a solid have a finite range and decay rapidly to zero.
Specifically, the elements of the force constant matrix in
Eq. <a href="#eq:fcmat" data-reference-type="ref"
data-reference="eq:fcmat">[eq:fcmat]</a> decrease as $1 / r^{5}$ with
interatomic distance. Consequently the dynamical matrix defined by
Eq. <a href="#eq:dmat" data-reference-type="ref"
data-reference="eq:dmat">[eq:dmat]</a> and its eigenvalues
$\omega^{2}({\mathbf{q}})$ are smoothly varying with wavevector
${\mathbf{q}}$. Fourier interpolation is used to generate dynamical
matrices on an arbitrarily fine grid or linear path in reciprocal space
from a set of DFPT calculations on a much coarser grid. For a full
description of the method see references (Baroni et al. 2001; Giannozzi
et al. 1991; Gonze and Lee 1997).

A complication arises in the case of polar solids where the
dipole-dipole interaction generated upon displacing an atom leads to a
longer ranged force constant matrix which decays only as $1 / r^{3}$.
CASTEP models this term analytically using Born effective charges and
dielectric permittivity calculated using an electric field response DFPT
calculation (see
section <a href="#sec:born-charges" data-reference-type="ref"
data-reference="sec:born-charges">[sec:born-charges]</a> and Ref. (Gonze
and Lee 1997)). It is therefore able to perform the Fourier
interpolation only for the part of the force constant matrix which
varies as $1 / r^{5}$, and does not require a finer grid than in the
case of non-polar solids.

##### Setting up an interpolation calculation

In the .cell file choose the q-points of the coarse grid of points at
which to perform DFPT calculations. This may be specified as a
`%block PHONON_KPOINT_LIST` containing the reduced set of points in the
irreducible Brillouin Zone, but is is almost always more straightforward
to use the alternative keywords

> phonon_kpoint_mp_grid : *p q r*  
> phonon_kpoint_mp_offset : INCLUDE_GAMMA

to specify the grid parameters and offset[^8]. The grid parameters *p*,
*q* and *r* should normally be chosen to give a roughly uniform sampling
of reciprocal space taking length of the reciprocal lattice vectors into
account, and should be compatible with the symmetry of the simulation
cell. Alternatively a grid may be specified using the minimum spacing,
for example

> phonon_kpoint_mp_grid_spacing 0.1 1/ang

Normally the grid should be chosen to contain the $\Gamma$ point, which
usually gives better convergence properties of the interpolation (*i.e.*
convergence at smaller *p,q,r*) than otherwise. (This is the opposite of
the convergence behaviour of electronic k-point sampling - see
Ref. (Probert and Payne 2003)). The special keyword value
`phonon_kpoint_mp_offset : INCLUDE_GAMMA` avoids the need to work out
the offset explicitly.

The choice of grid parameters *p*, *q* and *r* will be influenced by the
nature of the system under study. Ionically bonded systems tend to have
fairly short-ranged force-constant matrices and need relatively coarse
grids for convergence. For example, sodium chloride in the rocksalt
structure with a primitive lattice parameter 3.75Å a
$4 \times 4 \times 4$ grid is reasonably close to convergence. This
corresponds to a truncation of the force constant matrix at a distance
greater then $7.5$Å. On the other hand covalent systems tend to have
fairly long-ranged force-constant matrices and need finer grids for
convergence. Silicon is a good example of this. The primitive lattice
parameter is similar to sodium chloride at $3.81$Å, but the dispersion
curve is not fully converged until an $8 \times 8 \times 8$ grid is
used. A more detailed examination of this point may be found in
reference (Ackland, Warren, and Clark 1997).

Additional keywords in the .param file control the interpolation. Most
importantly

> `phonon_fine_method : INTERPOLATE`

instructs CASTEP to perform the Fourier interpolation step following the
usual DFPT calculation.

The target set of wavevectors for the interpolation is set up using
additional keywords or blocks in the *.cell* file. Either

> phonon_fine_kpoint_mp_grid *p q r*  
> phonon_fine_kpoint_mp_offset *o$_1$ o$_2$ o$_3$*

or

> phonon_fine_kpoint_mp_spacing 0.03 1/ang .

will perform interpolation onto a regular (possibly offset) grid. In
fact only the points in the irreducible wedge of the Brillouin Zone are
included, and a suitable weight is computed so that the weighted average
is identical to a uniform sampling of the BZ. This is the usual method
for computing a phonon density of states.

If a set of dispersion curves along high symmetry directions is
required, an empty cell keyword block

> `%block PHONON_FINE_KPOINT_PATH`  
> `%endblock PHONON_FINE_KPOINT_PATH`

will internally generate a list of ${\mathbf{q}}$-points sampling a
default path through the Brillouin-Zone according to the symmetry of the
calculation. Alternatively a list of ${\mathbf{q}}$-points explicitly
specifying the vertices of the path may be input as

>     %block PHONON_FINE_KPOINT_PATH
>     0.0 0.0 0.0
>     0.5 0.5 0.0
>     0.5 0.5 0.5
>     break
>     0.0 0.0 0.0
>     0.5 0.5 0.5
>     ...
>     %endblock PHONON_FINE_KPOINT_PATH

which traverses the directions between the vertices specified, except in
the presence of the break keyword where the path jumps without including
any intermediate points. The fineness of sampling along the path is set
by an additional keyword

> `phonon_fine_kpoint_path_spacing : 0.03 1/ang` .

As a final alternative, a simple list of points can be directly input to
model any sampling you choose using

>     %block PHONON_FINE_KPOINT_LIST
>     0.0 0.0 0.0
>     0.1 0.2 0.3
>     ...
>     %endblock PHONON_FINE_KPOINT_LIST

This is the approach used for files generated by Accelrys Materials
Studio.

The same keywords are used in the closely related method of a finite
difference calculation; see the example .cell file of
figure <a href="#fig:al-sc-cell" data-reference-type="ref"
data-reference="fig:al-sc-cell">2.4</a>.

##### Continuation

CASTEP stores the results of all phonon calculations - dynamical
matrices and force constant matrices - in the binary *seedname*.check
file written at the end of a successful run. This can be used to change
some values of parameters relating to interpolation, and to change or
indeed replace the entire fine phonon k-point set *without* any need to
repeat the expensive DFPT (or supercell) part of the calculation. In
fact a single calculation of the force constant matrix is sufficient to
compute a DOS at a variety of sampling densities plus arbitrarily smooth
dispersion curves.

Setting up a continuation calculation is simple. Just add the
continuation keyword in the usual way to a renamed copy of the param
file

> continuation *orig-seedname*.check

and make any changes fine k-point sampling parameters in a copy of the
*seedname*.cell file. It is recommended that you work with renamed
copies to avoid overwriting the original and valuable *seedname*.check
file. Running CASTEP on the new set of input files is exactly the same
as running a new calculation, except that the result will be generated
much more quickly.

When setting up a continuation run, take care not to change the standard
phonon k-point set, the electronic k-point set or any electronic
structure parameters such as elec_energy_tol. If a mismatch is detected
with the values stored in the checkpoint file, CASTEP will discard the
saved dynamical matrix data and restart the full calculation from the
beginning.

Continuation calculations are also used in conjunction with
*checkpointing* where a partially complete calculation is written to a
file, also in .check format and which may be used to restart an
interrupted calculation. See
section <a href="#sec:checkpointing" data-reference-type="ref"
data-reference="sec:checkpointing">6.2</a> for a description.

##### Control of interpolation scheme

The final representation of the Force constant matrix derived from the
dynamical matrices is actually a periodic representation and is
equivalent to the ${\mathbf{q}}=0$ dynamical matrix of a (fictitious)
$p \times q \times r$ supercell. (See
section <a href="#sec:supercell" data-reference-type="ref"
data-reference="sec:supercell">2.3.4</a> for more explanation.) CASTEP
must determine a mapping between elements of the periodic dynamical
matrix and aperiodic force constant matrix using a minimum-image
convention for ionic site pairs and impose a cutoff scheme in real space
to exclude (supercell-) periodic images. In fact CASTEP implements two
distinct schemes.

The cumulant scheme  (Parlinski, Li, and Kawazoe 1997) includes all
image force constants with a suitable weighting to avoid multiple
counting of images. This is achieved by including image force constants
in any direction if they lie within half the distance to the nearest
periodic repeat of the fictitious supercell lattice in that direction.
If an atom-atom pair vector lies *exactly* half way to the supercell
repeat so that image force constants occur, for example, at both
$\mathbf{L}$ and $-\mathbf{L}$, all images at the same distance are
included with a suitable weighting factor to preserve the symmetry of
the cumulant force constant matrix. (See Refs. (Parlinski, Li, and
Kawazoe 1997) and (Ye et al. 2004) for a more detailed explanation).
This scheme is selected by specifying the param file keyword

> phonon_fine_cutoff_method CUMULANT

and is in fact the default method in CASTEP.

CASTEP also implements a simple spherical cutoff, controlled by the
parameter $R_c$ and specified by parameter

> phonon_force_constant_cutoff 10.0 ang

The value $R_c$ should be chosen to satisfy
$2 R_c < min(p L_1,q L_2,r L_3) L$ where $L_n$ is the cell edge of the
simulation cell and $p,q,r$ are the (coarse) grid of phonon wavevectors.
It is usually easiest to specify a value of zero, in which case CASTEP
chooses the largest allowable value automatically. This scheme is most
suitable for bulk materials of cubic symmetry. The spherical scheme is
chosen using keyword

> phonon_fine_cutoff_method SPHERICAL .

Within the default method a smaller cutoff volume can be decreased by a
radius scaling factor, *e.g.*

> `phonon_force_constant_cut_scale : 1.0` .

This may be useful for testing the effect of long-ranged contributions
to the IFC matrix. However any departure of this parameter’s value from
1 does not preserve the superior convergence properties of the cumulant
scheme and will in general require a larger supercell than the exact
cumulant method.

##### Acoustic Sum Rule correction

The vibrational Hamiltonian is invariant to a uniform translation of the
system in space. This symmetry is the origin of the well-known result
that any crystal has three acoustic vibrational modes at
${\mathbf{q}}=0$ with a frequency of zero. Any book on lattice dynamics
will discuss the so-called *acoustic sum rule* (or ASR) which has
mathematical expressions for the force constant and $\Gamma$-point
dynamical matrices

$$\begin{aligned}
 \sum_{\kappa,a} {\Phi^{\kappa,\kappa^\prime}_{\alpha,\alpha^\prime}(a)}&= 0\\
 \sum_{\kappa} {D^{\kappa,\kappa^\prime}_{\alpha,\alpha^\prime}(\mathbf{q}=0)}&= 0
 \label{eq:asr}
\end{aligned}$$

In plane-wave calculations the translational invariance is broken as
atoms translate with respect to the fixed FFT grid, so the sum rule is
never exactly satisfied. Consequently it is sometimes observed even in
an otherwise apparently very well converged calculation that the three
acoustic modes at $\Gamma$ depart significantly from zero frequency.
Depending on the XC functional[^9] these frequencies may reach or exceed
50 cm$^{\text{-1}}$. One solution is simply to increase the FFT grid
parameters using parameter fine_grid_scale to increase the density grid
and not the wavefunction grid. However this can be prohibitively costly
in computer time and memory.

Provided that the amplitude of the symmetry violation is not too large,
it is possible to apply a transformation to the computed dynamical or
force constant matrix so that it exactly satisfies the ASR. CASTEP
implements a scheme which projects out the acoustic mode eigenvectors
and adjusts their frequency to zero, while having minimal impact on the
optic mode frequencies. This scheme is controlled by parameters (in the
.param file)

> `phonon_sum_rule : TRUE / FALSE`  
> `phonon_sum_rule_method : REAL / RECIP / REAL-RECIP / NONE`

The first of these simply activates or deactivates the correction. The
second chooses which of the variants of the ASR in
equation <a href="#eq:asr" data-reference-type="ref"
data-reference="eq:asr">[eq:asr]</a> to enforce, the real-space force
constant matrix, the reciprocal-space dynamical matrix, or both.
(`phonon_sum_rule_method : NONE` is a synonym for
`phonon_sum_rule : FALSE`). The real-space method is only applicable to
interpolation or supercell/finite displacement calculations, but the
reciprocal-space method can be used for any type of phonon calculation.

Both variants change the acoustic mode frequencies away from but near
${\mathbf{q}}=0$, the realspace method implicitly, and the reciprocal
space method explicitly, by determining the correction to
${D^{\kappa,\kappa^\prime}_{\alpha,\alpha^\prime}(\mathbf{q}=0)}$ and
subtracting it from
${D^{\kappa,\kappa^\prime}_{\alpha,\alpha^\prime}(\mathbf{q})}$ as
suggested by Gonze (Gonze and Lee 1997). This usually results in
acoustic mode behaviour which is indistinguishable from a very well
converged calculation using a very fine FFT grid[^10].

In addition to the sum rule on frequencies there is another for the Born
effective charges (see
section <a href="#sec:born-charges" data-reference-type="ref"
data-reference="sec:born-charges">[sec:born-charges]</a>)

$$\sum_{\kappa} {Z^{*}_{\kappa,\beta,\alpha} }= 0$$

which is activated by parameter

> `born_charge_sum_rule : TRUE`

The default behaviour of CASTEP is that neither sum rule is enforced. If
this was not requested in the original run, then it may be added in
post-processing fashion in a continuation run, as per
section <a href="#sec:continuation" data-reference-type="ref"
data-reference="sec:continuation">2.2.2</a>. Only the raw dynamical
and/or force constant matrices are stored in the checkpoint file
*without* the effect of ASR enforcement, which is only applied at the
printout stage. Therefore the effect can be turned off, or altered by a
post-processing calculation as well as turned on.

##### Density of States

The definition of the phonon density of states requires an integration
of the calculated frequencies $\omega_{i{\mathbf{q}}}$ over the
Brillouin-Zone

$$\label{eq:dos}
g(\omega) = \int_{\text{BZ}}d^3{\mathbf{q}}\sum_i \delta(\omega -
\omega_{i{\mathbf{q}}})  .$$

A simple computational approximation would be to compute
$\omega_{i{\mathbf{q}}}$ on a regular, discrete grid of ${\mathbf{q}}$,
replace the integral with a sum over the grid and convolute the result
with some peak shape function $h(\omega)$ such as a Gaussian

$$\label{eq:broadendos}
g(\omega) \approx \int d\omega^{\prime} h(\omega - \omega^{\prime}) \sum_{{\mathbf{q}}} \sum_i \delta(\omega^{\prime} -
\omega_{i{\mathbf{q}}})  .$$

In the CASTEP toolset, this is implemented by performing a Fourier
interpolation calculation onto a suitably fine ${\mathbf{q}}$-point set,
specifying either phonon_fine_kpoint_mp_grid or
phonon_fine_kpoint_mp_spacing. The resulting .phonon file may then be
analysed using the dos.pl tool (see
section <a href="#sec:dos-pl" data-reference-type="ref"
data-reference="sec:dos-pl">5.1</a>) which implements
equation <a href="#eq:broadendos" data-reference-type="ref"
data-reference="eq:broadendos">[eq:broadendos]</a>.

However the broadening smooths out any sharp features of the DOS, and
fails to reproduce the smoothly-curved or sharply peaked segments
typical of theoretical densities-of-states. This is discussed in more
detail in references (Yates et al. 2007) and  (Pickard 1997). A more
faithful rendering may be obtained using the so-called *adaptive
broadening* approach of Yates et. al (Yates et al. 2007) which uses the
*gradients* of the phonon branch dispersion to narrow or widen the
broadening for flat or steep branches respectively. It is activated by
the parameters keyword

> `phonon_calculate_dos : TRUE`

which computes an adaptively-broadened DOS during the Fourier
interpolation stage of the calculation and writes an output file named
\<seedname\>.phonon_dos containing the tabulated DOS, plus the
per-atomic-species resolved DOS in additional columns. The dos.pl script
is able to read and present .phonon_dos files just as with .phonon.

The range and resolution may be modified using additional parameters,
for example

> `phonon_dos_spacing : 0.01 THz`  
> `phonon_dos_limit : 10.0 THz`

Example output for the phonons of diamond is shown in
figure <a href="#fig:c2-pdos" data-reference-type="ref"
data-reference="fig:c2-pdos">2.3</a>.

<figure id="fig:c2-pdos">
<div class="center">
<img src="diamond-DOS.svg" style="width:60.0%" />
</div>
<figcaption>Phonon DOS for Diamond, with an <span
class="math inline">\(18^3\)</span> fine grid comparing Gaussian
broadened (red) vs adaptively broadened (black) methods. Artefactual
wiggles are clearly visible in the Gaussian-broadened case, and the
height of the sharp peaks is not well reproduced.</figcaption>
</figure>

#### Finite Displacement

In addition to the DFPT method of computing force constants, CASTEP
implements schemes based on numerical differentiation of forces when
atoms are displaced by a small amount from their equilibrium positions.
This method is useful for cases where DFPT has not been implemented,
which as of CASTEP release 24.1 includes ultrasoft pseudopotentials,
Hubbard U and exact exchange and hybrid functionals, and some of the
newer classes of dispersion correction.

There are three variants of this scheme.

##### Primitive Cell Finite Displacement

The basic finite displacement method is selected by setting parameter

> `phonon_method : FINITEDISPLACEMENT`

In contrast to DFPT such displacements are necessarily periodic with the
simulation cell, and therefore only ${\mathbf{q}}=0$ phonons are
commensurate with this condition. As in the case of DFPT lattice
dynamics the phonon wavevectors are specified by
`%block PHONON_KPOINT_LIST`, `%block PHONON_KPOINT_PATH` or
phonon_kpoint_mp_grid in the .cell file but only the $\Gamma$ point,
$(0,0,0)$ is meaningful. CASTEP will print a warning in the output file
and ignore any non-zero value in the list.

CASTEP proceeds by shifting each atom by a small amount, then performing
a SCF calculation to evaluate the forces on the perturbed configuration.
Both positive and negative displacements are performed in each direction
so that the corresponding force constants can be evaluated using the
accurate “central difference” method of numerical differentiation.

$$\label{eq:finite-diff}
 \frac{d^{2} E_{0}}{d u_{\kappa,\alpha}
        d u_{\kappa^{\prime},\alpha^{\prime}}} = \frac{d F_{\kappa,\alpha}}{d u} \approx \frac{F^{+}_{\kappa,\alpha} - F^{-}_{\kappa,\alpha}}{2 u} \;.$$

Equation <a href="#eq:finite-diff" data-reference-type="ref"
data-reference="eq:finite-diff">[eq:finite-diff]</a> demonstrates that a
single pair of displaced calculation yields an entire row of the
dynamical matrix. As with DFPT calculations, only the minimal set of
perturbations is performed and the space-group symmetry is used to build
the complete dynamical matrix.

The SCF calculations on the perturbed configuration are efficient,
typically taking only only a few cycles in CASTEP 5.0 or later. This
efficiency is achieved by first making a good guess for the electron
density of the perturbed system based on the ground state of the
unperturbed system, and applying a displacement of an atomic-like
density of the pseudo-atom in question. Second, the SCF is started using
the Kohn-Sham orbitals of the unperturbed state as the initial guess for
the perturbed configuration. To exploit this efficiency it is essential
to use the density-mixing (Davidson) minimiser, selected by
`elec_method : DM` in the .param file. (As the all-bands method has no
means of initialising the density).

The default displacement used is $0.01$ bohr This can be changed if
necessary by setting a parameter,

> `phonon_finite_disp : 0.02 ang`

in the .param file. Except for the differences discussed above, input
and output formats are the same as for DFPT calculations.

##### Born Charges, Permittivity and LO/TO splitting with FD

FD phonon calculations are useful where atomic-displacement response
DFPT is not implemented, notably for ultrasoft pseudopotentials (USPs)
and post-DFT exchange and correlation including LDA+U, SOC, hybrid
functionals (see table <a href="#tbl:captable" data-reference-type="ref"
data-reference="tbl:captable">2.1</a>). However properties including IR
spectral intensity, LO/TO splitting and the effectiveness of Fourier
interpolation of dynamical matrices depend on the Born effective charges
and dielectric permittivity (see
section <a href="#sec:efield" data-reference-type="ref"
data-reference="sec:efield">3</a>), whose calculation by electric field
response DFPT is not implemented.

In some cases it may be convenient to read in externally computed or
approximate values of Born charges and permittivity to be used to
calculate IR spectra, LO/TO splitting etc. For example, a DFT+U,
meta-GGA, or hybrid calculation might make use of values calculated
using PBE-DFPT and NCPs as a fair approximation, or even computed using
a different DFT code.

This may be accomplished by specifying the name of an external file
containing the values (in this case BORN.DAT) in the .param file

>     %block DEVEL_CODE
>     PHONON:READ_EXTERNAL_BORN=BORN.DAT:ENDPHONON
>     %endblock DEVEL_CODE

The file should be formatted as for the Phonopy code
(<https://phonopy.github.io/phonopy/input-files.html#born-optional>). A
file of this format may be written at the end of a DFT, NCP E-field
response calculation by

>     %block DEVEL_CODE
>     PHONON:WRITE_EXTERNAL_BORN=BORN.DAT:ENDPHONON
>     %endblock DEVEL_CODE

###### Born charges from Berry Phase Polarization

As of CASTEP release 25.1, it will be possible to calculate Born
effective charges (but not dielectric permittivity) using a
finite-difference numerical differentiation of the Berry-phase
polarization. This is performed as an adjunct to the FD phonon
calculation, and will be automatically selected in the case of an FD
phonon calculation with USPs.

##### Finite Displacement using non-diagonal supercells

As with the DFPT method, calculation of phonon dispersion and DOS using
the finite displacement method is achieved using Fourier interpolation
of dynamical matrixes. First, an approximation of the full
force-constant-matrix
${\Phi^{\kappa,\kappa^\prime}_{\alpha,\alpha^\prime}}$ is calculated by
explicitly computing a set of
${D^{\kappa,\kappa^\prime}_{\alpha,\alpha^\prime}(\mathbf{q})}$ with
${\mathbf{q}}$ sampling a coarse grid of ${\mathbf{q}}$ and using the
inverse Fourier relation to
equation <a href="#eq:dmat" data-reference-type="ref"
data-reference="eq:dmat">[eq:dmat]</a>. Then
equation <a href="#eq:dmat" data-reference-type="ref"
data-reference="eq:dmat">[eq:dmat]</a> is used to generate the dynamical
matrices at all of the (fine) phonon wavevectors required.

The key step of computing every
${D^{\kappa,\kappa^\prime}_{\alpha,\alpha^\prime}(\mathbf{q})}$ on the
“coarse” regular grid is achieved by constructing a list of supercells
commensurate with each of the coarse ${\mathbf{q}}$-points in turn (so
that ${\mathbf{q}}\cdot T_{\text{SC}} = 2 n \pi$ where $T_{\text{SC}}$
is a lattice vector of that supercell). This ${\mathbf{q}}$ maps to a
$\Gamma$ point calculation on the corresponding supercell despite being
incommensurate with the primitive cell, and
${D^{\kappa,\kappa^\prime}_{\alpha,\alpha^\prime}(\mathbf{q})}$ for the
primitive cell at ${\mathbf{q}}$ may be extracted using
equation <a href="#eq:dmat" data-reference-type="ref"
data-reference="eq:dmat">[eq:dmat]</a>. It was shown by Lloyd-Williams
and Monserrat (Lloyd-Williams and Monserrat 2015) that by using
so-called non-diagonal supercells such a set need contain only the least
common multiple of the coarse grid subdivisions ($p, q, r$).
Consequently the number of distinct supercell calculations requires is
much smaller than the number of ${\mathbf{q}}$-vectors of the grid, and
most importantly increases with the linear dimension of the grid instead
of the product.

A non-diagonal supercell calculation is selected by setting the
parameters

> `phonon_method : FINITEDISPLACEMENT`  
> `phonon_fine_method : INTERPOLATE`,

the cell keywords

> `phonon_kpoint_mp_grid : p q r`

and one of the specifications of phonon_fine_kpoints exactly as in the
case of a DFPT calculation. When the phonon calculation begins, CASTEP
generates a list of supercells using the algorithm of Lloyd-Williams and
Monserrat (Lloyd-Williams and Monserrat 2015). For each of these, CASTEP
generates the corresponding supercell and performs a finite-displacement
phonon calculation. The dynamical matrices
${D^{\kappa,\kappa^\prime}_{\alpha,\alpha^\prime}(\mathbf{q})}$ are then
extracted for each ${\mathbf{q}}$ on the coarse KPOINT grid, and the
calculation proceeds exactly as if the DFPT method has been used.

The choice of coarse ${\mathbf{q}}$-point sampling set by
phonon_kpoint_mp_grid requires one additional consideration compared to
the DFPT case to give best computational performance. The
supercell-generation algorithm of Lloyd Williams and Monserrat
guarantees that the largest supercell in the generated list contains
$N_\text{SC} = \text{LCM}  (p, q, r)$ replications of the primitive
cell. For example, an $4 \times 4 \times 4$ sampling of a cubic cell
gives $N_\text{SC} = 4$. However a $2 \times 3 \times 4$ sampling of an
orthorhombic primitive cell yields $N_\text{SC} = 12$, the LCM of 4 and
3. In that case it might be more efficient to use a $2 \times 4
\times 4$ sampling which only requires $N_\text{SC} = 8$ for the largest
of the generated supercells.

As well as the ${\mathbf{q}}$-point sampling, FD interpolation
calculations require one further parameter choice, namely to select the
*electronic* k-point sampling for the supercell calculations. This may
most conveniently be set using the cell keyword

> `supercell_kpoint_mp_spacing : 0.1 1/ang`.

An explicit setting of grid sampling divisions with
supercell_kpoint_mp_grid would be inappropriate as the correct choice
will be different for each supercell, and unknown until run-time. By
default CASTEP chooses a default spacing which is consistent with the
primitive cell spacing, in an attempt to achieve a consistent level of
convergence across supercells.

Apart from the above two considerations, the setting up and execution of
the FD interpolation follows exactly the same lines as with DFPT. In
particular sum-rule, DOS and thermodynamics calculations and
interpolation scheme tweaks apply identically.

##### Legacy Finite Displacement/Supercell

<figure id="fig:al-sc-cell">
<pre><code>%block LATTICE_CART
          0     2.02475      2.02475
    2.02475           0      2.02475
    2.02475     2.02475            0
%endblock LATTICE_CART

%block POSITIONS_ABS
Al          0       0       0
%endblock POSITIONS_ABS

%block SPECIES_POT
Al Al_00.usp
%endblock SPECIES_POT

kpoint_mp_grid           14 14 14
supercell_kpoint_mp_grid  2  2  2

symmetry_generate

%block PHONON_FINE_KPOINT_PATH
0.0 0.0 0.0
0.5 0.5 0.0
1.0 1.0 1.0
0.5 0.5 0.5
0.5 0.5 0.0
0.5 0.75 0.25
0.5 0.5 0.5
%endblock PHONON_FINE_KPOINT_PATH

%block PHONON_SUPERCELL_MATRIX
-3  3  3
 3 -3  3
 3  3 -3
%endblock PHONON_SUPERCELL_MATRIX</code></pre>
<figcaption>Example cell file for aluminium supercell phonon
calculation. This calculation computes a set of dispersion curves along
high-symmetry directions. Note that this calculation is not fully
converged with supercell size - there is a noticeable change in
frequency on some of the branches on increasing the supercell matrix
entries from 3 to 4.</figcaption>
</figure>

<figure id="fig:-sc-paramal">
<pre><code>task                            : PHONON
phonon_fine_method              : SUPERCELL
phonon_calc_lo_to_splitting     : FALSE
phonon_sum_rule                 : TRUE
calculate_born_charges          : FALSE
phonon_force_constant_ellipsoid : 1.0

cut_off_energy    : 150 eV
metals_method     : DM
smearing_width    : 0.04 eV
nextra_bands      : 4
spin_polarized    : FALSE

opt_strategy      : SPEED
num_dump_cycles   : 0
xc_functional     : LDA</code></pre>
<figcaption>Example param file for aluminium supercell phonon
calculation. It is not strictly necessary to turn off the LO/TO
splitting calculation - CASTEP will warn that this is not possible and
turn it off anyway. The explicit request for a non-spin polarized
calculation is necessary for fcc Al, because CASTEP chooses
spin-polarized by default due to the odd number of electrons.
</figcaption>
</figure>

The limitation of the primitive-cell finite displacement approach to
${\mathbf{q}}=0$ may also be overcome by combining the method with the
use of a supercell. This method, sometimes known as the “direct
method” (Kunc and Martin 1982, sec. 19.2) relies on the short-ranged
decay of the force constant matrix with interatomic distance and makes
the assumption that force constants for separations larger than some
value, $R_c$ are negligibly small and can be treated as zero. A
supercell can be constructed to contain an imaginary sphere of radius
$R_c$ beyond which force constants may be neglected. Then the *dynamical
matrix* of a supercell satisfying $L > 2 R_{c}$ is identical to the
*force constant matrix*, *i.e.*
${C^{\kappa,\kappa^\prime}_{\alpha,\alpha^\prime}(\mathbf{q}=0)}= {\Phi^{\kappa,\kappa^\prime}_{\alpha,\alpha^\prime}}$.
Therefore complete knowledge of the force constant matrix to a
reasonable approximation may be derived from a single calculation of the
${\mathbf{q}}=0$ dynamical matrix of a supercell containing several
primitive cells. From this set of force constants the dynamical matrices
at any phonon wavevector may be computed using
equation <a href="#eq:dmat" data-reference-type="ref"
data-reference="eq:dmat">[eq:dmat]</a> in exactly the same way as used
in an interpolation calculation.

**Note:** This method is mostly superceded by the non-diagonal
supercell/finite displacement approach of
section <a href="#sec:ndsc" data-reference-type="ref"
data-reference="sec:ndsc">2.3.3</a>. That method gives an equivalent
calculation of the force constant matrix, usually with a smaller
computational cost, because of its superior scaling with cell size.

In a CASTEP calculation the supercell must be chosen and explicitly
specified in the input files. It is defined by a matrix, which
multiplies the ordinary simulation cell vectors and is specified as a
$3\times3$ matrix in the .cell file of the form

>     %block PHONON_SUPERCELL_MATRIX
>      4  0  0
>      0  4  0
>      0  0  4
>     %endblock PHONON_SUPERCELL_MATRIX

This typical example using a diagonal matrix creates a supercell
expanded along each lattice vector **a**, **b** and **c** by factor of
4.

The supercell to be used in a phonon calculation needs to be chosen with
care.

- It must be large enough to contain a sphere of radius $R_c$, the
  typical range of a force constant matrix. In a simple metal a value of
  $R_c$ as small as 6Å may be satisfactory, but more complex and
  structured systems will need a larger value[^11]. In a covalent
  semiconductor the required $R_c$ may be 10Å or larger.

- However a supercell can quickly grow to generate a very large
  calculation indeed. For example a $4 \times 4 \times 4$ supercell of
  even a 2 atom primitive cell contains 128 atoms. Consequently it is
  rarely feasible to use this method on a uniprocessor or desktop
  computer, and a substantial parallel computer is usually required.

- The *shape* of the supercell should usually be as near cubic as
  possible, irrespective of the shape of the primitive cell, to optimise
  the supercell size/$R_c$ ratio. This can be achieved using an
  off-diagonal supercell matrix, as in the example input files of
  figure <a href="#fig:al-sc-cell" data-reference-type="ref"
  data-reference="fig:al-sc-cell">2.4</a> where a cubic supercell
  containing 108 atoms ($3 \times 4 = 12$ unit cells) is generated from
  the rhombohedral primitive cell of aluminium containing a single atom.

- If the system is highly anisotropic, perhaps a slab model of a surface
  then a uniform supercell is clearly not appropriate and a suitable
  supercell must be designed from a consideration of the exact nature of
  the model.

The electronic Brillouin-Zone integrals for the supercell calculation
use the special k-points method, which are specified in the .cell file
using an separate but analogous set of keywords to those pertaining to
the primitive cell sampling. Specifically,

- `%block SUPERCELL_KPOINT_LIST` allows an exact specification of
  k-points and weights

- supercell_kpoint_mp_grid_spacing d chooses a Monkhorst-Pack grid with
  the specified spacing.

- supercell_kpoint_mp_grid p q r (and optionally
  supercell_kpoint_mp_offset ) allow an exact specification of a
  Monkhorst-Pack grid.

- Finally if no supercell kpoint keywords or blocks are given, a grid is
  chosen to generate a similar sampling density to the primitive cell
  calculation.

Once the force constant matrix has been calculated using the supercell,
the remainder of the lattice dynamics proceeds exactly as in the case of
a Fourier interpolation calculation. The keywords controlling the
interpolation scheme and cutoff radius, the fine phonon kpoint set and
the acoustic sum rule enforcement work in exactly the same way. See
section <a href="#sec:interpolation-setup" data-reference-type="ref"
data-reference="sec:interpolation-setup">2.2.1</a> for details.

#### Convergence

Phonon and dielectric response calculations give rise to a number of
issues with convergence, in addition to those encountered in
ground-state calculations, and a systematic and step-by-step approach is
necessary to achieve well-converged results. CASTEP assists the user by
setting default values for many convergence parameters based on the
particular task which incorporate the experience of many calculations.
Therefore it is a good practice not to specify convergence parameters
explicitly unless specific convergence test results are known -
otherwise a well-chosen default could be overridden with an untested
value. (The example input of
figure <a href="#fig:example-gamma" data-reference-type="ref"
data-reference="fig:example-gamma">2.1</a> contains no phonon-specific
convergence parameters.)

It was emphasised previously
(section <a href="#sec:geometry" data-reference-type="ref"
data-reference="sec:geometry">1.2.1</a>) that a well-converged geometry
optimisation is a prerequisite for a phonon calculation. In turn this
mandates a reasonably high level of convergence of plane-wave cutoff,
SCF convergence and electronic k-point sampling. It is typical to run
phonon calculations at a `basis_precision : FINE` level of plane-wave
cutoff, (The header of each .usp or .recpot pseudopotential file
contains a translation into eV units.) DFPT electric field response
calculations can require finer electronic k-point sampling than suffices
for a ground state calculation, so the effect of the kpoint_mp_grid on
the dielectric permittivity should be tested.

##### Convergence of DFPT solver

The second stage of a phonon or E-field calculation is the variational
DFPT solver, and there are a number of associated parameters to control
the convergence. Parameters

> phonon_energy_tol 1.0e-5 eV/ang\*\*2

and

> efield_energy_tol 1.0e-5 ang\*\*3

govern the exit criterion for the DFPT self-consistency loop. The above
default values are usually sufficient for frequencies converged to
$< 1 \text{cm}^{-1}$ and permittivities to two or three decimal places
and rarely need to be changed.

However the converged results depend not only on the DFPT SCF
calculation, but also, and rather strongly on the degree of convergence
of the *ground state* electron density and wavefunctions. The
ground-state wavefunctions enter the DFPT equations both directly and as
a consequence of the orthogonality condition between ground-state and
first-order response orbitals (Refson, Tulip, and Clark 2006; Gonze and
Lee 1997). Specifically, the error in the second-order energy or force
constants is variational, and therefore depends *quadratically* on the
error in the DFPT first-order response orbitals. However it is
*non-variational* and *linear* with respect to the error in the
ground-state orbitals. The practical consequence is that there is an
error in DFPT results which varies as the *square root* of the
ground-state convergence parameter elec_energy_tol which therefore must
be very small for good accuracy. A heuristic rule of thumb is that
elec_energy_tol $\approx$ phonon_energy_tol$^2$ is necessary to converge
the second order energy to the value of phonon_energy_tol. The default
setting of elec_energy_tol is based on this heuristic and should only be
changed after careful testing of the effect on frequencies or dielectric
properties.

##### Convergence of finite displacement forces

The criteria for choosing the “ground-state” convergence parameters for
a finite-displacement phonon calculation (in either primitive or
supercell mode) are naturally almost identical to those governing DFPT
calculations. In this case there is no additional self-consistent
electronic calculation; instead there is a sequence of additional
ground-state calculations at displaced geometries.

Finite-displacement calculations require very well converged *forces* to
be computed in the ground state SCF calculation. This is because the
numerical evaluation of the second derivative depends on differences
between the (small) forces at perturbed configurations
(equation <a href="#eq:finite-diff" data-reference-type="ref"
data-reference="eq:finite-diff">[eq:finite-diff]</a>). Clearly the
numerical derivative has lower relative precision than the argument (the
forces), which must therefore be evaluated to a rather high precision.
Unlike the ground-state *energy*, which is variational with respect to
the orbitals, the forces are not and the error in the forces is linear
in the error in the orbitals. The main parameter governing their
accuracy is again

> elec_energy_tol

which is again set by default from the value of phonon_energy_tol$^2$.

It is normally sufficient to use the default value of elec_energy_tol
chosen as the default in a phonon calculation. However there is also a
way of *directly* achieving the requested force tolerance, by setting
parameter

> elec_force_tol 1e-3 eV/ang .

This is not used by default as it adds some overhead to the SCF
calculation resulting in longer run times. However it does guarantee
that the forces really are converged to the required accuracy, unlike
setting elec_energy_tol.

##### Convergence of Interpolation and Supercell

Fourier interpolation and supercell calculations add yet another
criterion which must be satisfied to achieve well-converged results -
the range of the force-constant matrix in real space.
Sections <a href="#sec:interpolation-setup" data-reference-type="ref"
data-reference="sec:interpolation-setup">2.2.1</a>
and <a href="#sec:supercell" data-reference-type="ref"
data-reference="sec:supercell">2.3.4</a> discuss some of the criteria
applicable to these types of calculation.

Unfortunately convergence testing of the size of the supercell for
FD/supercell calculations can be prohibitively expensive as the the
volume and number of atoms increase as the cube of the linear dimension
under test. Consequently the CPU time will increase with the sixth power
or higher! There is no easy solution to this problem and the reader is
cautioned not to necessarily take on trust that FD/supercell phonon
calculations published in the literature are fully converged!

However some progress may be made by performing a calculation using the
largest feasible supercell. If

> `phonon_fine_cutoff_method : SPHERICAL`

is selected, it is then possible to vary the cutoff radius using
parameter

> phonon_force_constant_CUTOFF

to lower values than the maximum of half the largest box size. These
tests may be performed using the “continuation” method
(section <a href="#sec:continuation" data-reference-type="ref"
data-reference="sec:continuation">2.2.2</a>) or the “phonons” utility
(section <a href="#sec:phonons-tool" data-reference-type="ref"
data-reference="sec:phonons-tool">5.4</a>) without repeating the
expensive supercell calculation.

Another convergence test is automatically performed if the spherical
cutoff method is selected. The frequencies resulting from the
interpolation with the spherical cutoff are compared to those from an
“exact” interpolation at wavevectors commensurate with the supercell.
The results are written to the .castep file. This provides a
quantitative measure of the interpolation error.

A similar cutoff scaling test might be possible for the
`phonon_fine_cutoff_method : DEFAULT` case, although it is considerably
less useful. Parameter

> phonon_force_constant_cut_scale 0.9

will scale the range of the cutoff by the factor specified. However the
default value of 1.0 is “special” in including constants separated by
exactly half the supercell Wigner Seitz cell. Consequently convergence
is not smooth in this parameter and while 0.95 might be underconverged,
1.0 could be very well converged. (This is also the reason for its
superior convergence performance over the spherical cutoff scheme.)

### Dielectric properties

<span id="sec:born-charges" label="sec:born-charges"></span>

Density functional perturbation theory is not limited to atomic
displacement perturbations, and may also be used to calculate other
response properties with respect to an electric field perturbation of an
insulating system[^12]. These are the dielectric permittivity

$$\epsilon_{\alpha,\beta}^{\infty} = \delta_{\alpha,\beta} - \frac{4 \pi}{\Omega}\frac{\partial^2 E}{\partial \varepsilon_{\alpha} \partial \varepsilon_{\beta}} \;,$$

the “molecular” polarizability

$$\alpha_{\alpha,\beta}^{\infty} =  - \frac{1}{\Omega}\frac{\partial^2 E}{\partial \varepsilon_{\alpha} \partial \varepsilon_{\beta}} \;,$$

and the Born effective charges

$$Z^{*}_{\alpha,\beta} = \Omega \frac{\partial P_{\beta}}{\partial u_{\kappa,\alpha}} = \frac{\partial F_{\kappa,\alpha}}{\partial  \varepsilon_{\beta}} = 
\frac{\partial^2 E}{\partial u_{\kappa,\alpha}\partial  \varepsilon_{\beta}}$$

which play a strong role in lattice dynamics of crystals, and in
particular governs the frequency-dependent dielectric response in the
infra-red region (Gonze and Lee 1997).

An electric field response calculation is selected using the task
keyword in the .param file.

> `task : EFIELD`

which computes the optical frequency dielectric permittivity tensor and
the low frequency (ionic lattice) response to a time-varying field in
the regime of the phonon modes and the Born charges. Alternatively

> `task : PHONON+EFIELD`

performs both dielectric and phonon tasks.

The low and near-infrared frequency contributions to the permittivity
and polarizability and the Born effective charges are printed to the
.castep output file. An extract for the same Wurtzite BN calculation as
earlier is shown in
figure <a href="#fig:efield-out" data-reference-type="ref"
data-reference="fig:efield-out">3.1</a>. An additional output file
*seedname*.efield is also written and contains the frequency-dependent
permittivity over the entire range, with a spacing determined by
parameter efield_freq_spacing, and a Lorentzian broadening governed by a
fixed $Q$, efield_oscillator_q.

The low-frequency contribution of the phonons to the dielectric
polarizability and permittivity is only well-defined when all mode
frequencies are real and positive. In the presence of any imaginary mode
or one of zero frequency the low-frequency dielectric and polarizability
tensors are not calculated and are reported as “infinity”. By default
the three lowest frequency modes are assumed to be acoustic modes and
not included in the calculation. To support molecule-in-supercell
calculations, the parameter efield_ignore_molec_modes may be set to
molecule, which excludes the 6 lowest frequency modes from the
dielectric calculation. Allowed values are CRYSTAL, MOLECULE and
LINEAR_MOLECULE which ignore 3, 6 and 5 modes respectively.

In fact CASTEP usually performs an electric field response calculation
even for a `task : PHONON` calculation because the permittivity tensor
and Born charges are required to calculate the LO/TO splitting terms.
Conversely a pure `task : EFIELD` calculation also performs a
$\Gamma$-point phonon calculation which is needed to compute the ionic
contribution to the permittivity and polarizability. Consequently the
only real difference between any of the tasks PHONON, EFIELD and
PHONON+EFIELD lies in what is printed to the output file as the same
computations are performed in each case. Only if one or other of those
properties is specifically turned off with one of the parameters

> `phonon_calc_lo_to_splitting : FALSE`  
> `efield_calc_ion_permittivity : FALSE`

is a “pure” phonon or E-field calculation ever performed.

It is well documented that the LDA tends to overestimate dielectric
permittivities - by over 10% in the case of Si or Ge (Levine and Allan
1989). It is possible to include an ad-hoc correction term to model the
missing self-energy, by applying the so-called “scissors operator”,
which consists of a rigid upshift of all conduction band states. This
was incorporated into DFPT electric field response calculations by Gonze
and Lee (Gonze and Lee 1997). The parameters keyword

> excited_state_scissors 1.0 eV

is used to model the effect of a 1 eV (in this example) upshift of
conduction band states and will include the effects on dielectric
permittivity and Born charges (but not phonons). The value to use must
be determined from high-level calculations or empirically.

<figure id="fig:efield-out">
<pre><code> ===============================================================================
        Optical Permittivity (f-&gt;infinity)             DC Permittivity (f=0)
        ----------------------------------             ---------------------
         4.50788     0.00000     0.00000         6.64363     0.00000     0.00000
         0.00000     4.50788     0.00000         0.00000     6.64363     0.00000
         0.00000     0.00000     4.63846         0.00000     0.00000     7.15660
 ===============================================================================</code></pre>
<pre><code> ===============================================================================
                                    Polarisabilities (A**3)
                Optical (f-&gt;infinity)                       Static  (f=0)
                ---------------------                       -------------
         6.52844     0.00000     0.00000        10.50324     0.00000     0.00000
         0.00000     6.52844     0.00000         0.00000    10.50324     0.00000
         0.00000     0.00000     6.77146         0.00000     0.00000    11.45793
 ===============================================================================
 ===================================================
                   Born Effective Charges
                   ----------------------
   B       1         1.84636     0.00000     0.00000
                     0.00000     1.84636     0.00000
                     0.00000     0.00000     1.94523
   B       2         1.84636     0.00000     0.00000
                     0.00000     1.84636     0.00000
                     0.00000     0.00000     1.94523
   N       1        -1.85371     0.00000     0.00000
                     0.00000    -1.85371     0.00000
                     0.00000     0.00000    -1.94009
   N       2        -1.85371     0.00000     0.00000
                     0.00000    -1.85371     0.00000
                     0.00000     0.00000    -1.94009
 ===================================================</code></pre>
<figcaption>Extract from the <span>.castep</span> output file generated
from the hexagonal BN run of figure <a href="#fig:example-gamma"
data-reference-type="ref" data-reference="fig:example-gamma">2.1</a>,
with <code>task : EFIELD</code>. The Born effective charges are laid out
with the columns representing the X,Y,Z electric field directions and
the rows the X,Y,Z displacement directions.</figcaption>
</figure>

#### Non-linear optical susceptibility

In addition to the linear response properties calculated with task
EFIELD or phonon+efield the non-linear dielectric susceptibility may be
computed if the parameter

> `efield_calculate_nonlinear : CHI2`

is set. The calculation uses the “2n+1 theorem” extension of
DFPT (Baroni et al. 2001; Miwa 2011) to compute the static response,
$\chi^{(2)}$. The results are reported in the .castep file in reduced
tensor form as the “d-matrix” (Boyd 2003). See
figure <a href="#fig:nlo" data-reference-type="ref"
data-reference="fig:nlo">3.2</a>.

This is not activated by default, because it is substantially more
expensive than a baseline E-field linear response calculation. (It
requires the calculation of three sets of response functions with a full
k-point set in P1 symmetry, even for crystals with higher, even cubic
symmetry.)

<figure id="fig:nlo">
<pre><code> ===========================================================================
      Nonlinear Optical Susceptibility (pm/V)
      ---------------------------------------
         1.13621    -1.13625     0.00001     0.00002    -7.73417     0.00002
         0.00002    -0.00003    -0.00008    -7.73421     0.00002    -1.13625
        -7.73417    -7.73421   -30.12197    -0.00008     0.00001     0.00002
 ===========================================================================</code></pre>
<p>. <span id="fig:nlo" label="fig:nlo"></span></p>
<figcaption>Nonlinear optical susceptibility <span
class="math inline">\(\chi^{(2)}\)</span> expressed as a d-matrix for
LiNbO3</figcaption>
</figure>

### Thermodynamic properties

One of the important motivations for lattice dynamical calculations of
crystalline solids is that the harmonic approximation gives access to
thermodynamic properties including the zero-point energy and the free
energy as a function of temperature. CASTEP lattice dynamics
calculations can be followed by a thermodynamics calculation to
calculate the zero-point energy and temperature dependent free energy,
entropy, and specific heat

$$\begin{aligned}
 F &= E_{\text{elec}} + k_{B} T \sum_{{\mathbf{q}},i} ln \left [ 2 \sinh \left (\frac{\hbar \omega_{{\mathbf{q}},i}}{ k_{B} T} \right ) \right ],\\
 S &= \frac{1}{2T} \sum_{{\mathbf{q}},i} \hbar \omega_{{\mathbf{q}},i}  \coth \left (\frac{\hbar \omega_{{\mathbf{q}},i}}{ k_{B} T} \right )
       - k_{B} \sum_{{\mathbf{q}},i} ln \left [ 2 \sinh \left (\frac{\hbar \omega_{{\mathbf{q}},i}}{ k_{B} T} \right ) \right ],\\
 C_V &=  \sum_{{\mathbf{q}},i} k_{B} \left ( \frac{\hbar \omega_{{\mathbf{q}},i}}{ k_{B} T} \right )^2
         \exp{\left (\frac{\hbar \omega_{{\mathbf{q}},i}}{ k_{B} T}\right
         )}  /%
         \left [\exp{\left (\frac{\hbar \omega_{{\mathbf{q}},i}}{ k_{B}
         T}\right )} -1 \right ]^2 .
\end{aligned}$$

The thermodynamics calculation follows a previous phonon calculation. It
is selected by setting the parameter

> `task : THERMODYNAMICS`

This may be used in a continuation run from a previous phonon
calculation, where the value of the continuation parameter is the name
of the previous .check file. Alternatively it may be configured as a new
run from scratch by setting the remainder of the parameters exactly as
if this were a phonon task, and the phonon calculation will be performed
first. Only the phonons defined on the “fine” set of phonon kpoints will
be used to compute the free energy as it is normally expected that a
thermodynamics calculation will follow an interpolation or supercell
calculation. However provided that the standard and fine sets of phonon
k-points are identical, it may also be used following a standard phonon
calculation.

This task will compute and print the free energy, entropy and specific
heat plus the vibrational atomic displacement parameters (“ADP”s) in the
range of temperatures specified by the parameters thermo_t_start and
thermo_t_stop. The number of temperatures is set by one or other of the
parameters thermo_t_spacing or thermo_t_npoints. All temperatures are
absolute and the default unit is K. The results of the calculation are
written to the .castep file.

The harmonic approximation free energy is only defined if all
frequencies are greater than or equal to zero. Any zero or imaginary
frequencies are automatically omitted from the calculation and a warning
message is printed. It is the responsibility of the end user to check
that the computed free energy is not rendered meaningless by the
presence of an imaginary mode.

### Plotting and analysis tools

Additional tools and programs are available to help automate
post-processing analysis and plotting of the results of phonon
calculations.

Two programs dispersion.pl and dos.pl analyse the .castep or .phonon
files, read the phonon information and can generate plots of phonon
dispersion curves across the Brillouin Zone, and of phonon densities of
states respectively [^13]. These programs are written in the PERL
language which is almost universally available on modern operating
systems. Generation of the plots is handled by either of the
[xmgrace](http://plasma-gate.weizmann.ac.il/Grace/) or
[gnuplot](http://gnuplot.info/) graphics programs (the perl programs
generate an xmgrace script and invoke it automatically). [^14]

#### dos.pl

The dos.pl program can read any of a .castep, a .phonon, or a
.phonon_dos file and generate a phonon density of states plot. Arguments
are given in “unix style” following a minus sign (dash). The command

> dos.pl -xg -w 10 \<seedname\>.phonon

will generate a phonon density of states using a Gaussian broadening
with a FWHM of 10 cm$^{-1}$, and invoke xmgrace to generate the plot. If
a .phonon file is specified the DOS is constructed as a weighted average
over all q-points present in the file. Alternatively an
adaptively-weighted dos computed by CASTEP is read from a .phonon_dos
file if given as command-line argument.

The options -np, -ps, -eps change the default output and write an
xmgrace file, a PostScript or Encapsulated PostScript file respectively.
These files are written to standard output so use a shell redirect ($>$)
to name your plot file. The option -ir option weights the computed DOS
by the computed infra-red intensity in the .castep or .phonon file. This
simple algorithm is not a fully realistic model of an infrared powder
spectrum, and should be regarded as a simple approximation only. [^15]
Likewise the -raman option extracts Raman intensities and computes the
weighted spectrum. The computed spectrum does include the
frequency-dependent and Stokes thermal factors, and is therefore a more
realistic model spectrum than in the infrared case. The -lorentz option
switches from Gaussian to a Lorentzian broadening and the -temp T option
sets the temperature for the Stokes thermal population term in Kelvin.

Figure <a href="#fig:bn-dos-ir" data-reference-type="ref"
data-reference="fig:bn-dos-ir">5.1</a> shows an example derived from the
output produced by using dos.pl on the example run of
figure <a href="#fig:example-gamma-out" data-reference-type="ref"
data-reference="fig:example-gamma-out">2.2</a>.

<figure id="fig:bn-dos-ir">
<div class="center">
<img src="BN_irdos.svg" style="width:60.0%" />
</div>
<figcaption>Example output from dos.pl based on the run of figure <a
href="#fig:example-gamma-out" data-reference-type="ref"
data-reference="fig:example-gamma-out">2.2</a>. Infrared spectrum and
DOS curves based on just the TO modes or TO plus LO have been combined
into into one plot with a slightly shifted baseline, scaled, and legends
added.</figcaption>
</figure>

#### dispersion.pl

The dispersion.pl program can read either a .castep or a .phonon file
and generate a dispersion curve plot using xmgrace. Unlike the behaviour
of dos.pl there is an important difference between the behaviour when
reading these two different output files concerning the detection of
branch crossings. A high-quality dispersion plot requires that phonon
branches are drawn as continuous lines even when two branches cross in
between the computed wavevectors. Dispersion.pl contains an algorithm
based on matching of eigenvectors at adjacent points to determine branch
connectivity. Only the .phonon file contains eigenvector information, so
only in this case is crossing detection enabled. The algorithm can be
time consuming and take several minutes to complete in large cases, so
patience is sometimes required. The -nj option (“no-join”) disables the
crossing detection even when the input file is a .phonon one.

The options -np, -ps, -eps behave exactly as for dos.pl. One useful
output option is -symmetry \<symm\> which attempts to label the high
symmetry points using the conventional Brillouin zone notation of
Bradley and Cracknell. The symmetry keywords cubic, fcc, bcc,
tetragonal, tetragonal, tetragonal-I, orthorhombic, hexagonal, trigonal,
trigonal-h (and minor variants) are understood.

Figure <a href="#fig:rbbr-dispersion" data-reference-type="ref"
data-reference="fig:rbbr-dispersion">5.2</a> demonstrates the effect of
the flags and the branch joining algorithm. The plots were produced from
a Fourier interpolation calculation of fcc RbBr using the commands

> dispersion.pl -xg -symmetry fcc RbBr.phonon

and

> dispersion.pl -xg -symmetry fcc -nj RbBr.phonon

respectively.

<figure id="fig:rbbr-dispersion">
<table>
<tbody>
<tr class="odd">
<td style="text-align: center;"><img src="RbBr-join.svg"
alt="image" /></td>
<td style="text-align: center;"><img src="RbBr-nj.svg"
alt="image" /></td>
</tr>
</tbody>
</table>
<figcaption>Phonon dispersion curve plots of RbBr generated using the
<span>dispersion.pl</span> script and xmgrace. The Brillouin zone
labelling is generated using the <span>-symmetry fcc</span> option. The
left-hand plot was generated using the default branch crossing detection
algorithm, which was disabled using the <span>-nj</span> option for the
right-hand plot. The algorithm has discriminated between modes which do
cross and the four genuine avoided crossings in the left-hand
plot</figcaption>
</figure>

#### mode_follow

mode_follow is one of the Fortran tools suite in the CASTEP source,
which is compiled using the command make tools. As the name implies its
function is to generate new “frozen phonon” configurations based on
perturbation by a mode generated from a previous phonon calculation,
which it outputs by writing one or more new .cell files. In fact it has
two modes of operation.

1.  To generate a sequence of .cell files perturbed by the a frozen
    phonon at a range of amplitudes, which may be used to explore the
    energy profile along the mode it is invoked as

    > mode follow -mode *mode_num* -namp *num_amplitudes* -amp
    > *max_amplitude* -qpoint *qx qy qz*
    > $<$*seedname*$>$\|$<$*seedname*$>$.phonon

    which reads the .phonon and corresponding .cell files and generates
    a sequence of $N+1$ files *seedname*-i.cell containing structures
    perturbed by the selected eigenvector scaled by a non-dimensional
    amplitude factor $f=A i/N, i=0 .. N$. The arguments are

    *mode_num*  
    is the integer index number selecting which mode to use (default 1)

    *qx qy qz*  
    is the q-point to extract from the .phonon file (default (0,0,0))

    *max_amplitude*  
    is $A$ the non-dimensional scale of the eigenvectors used to create
    the (largest) displacement

    *num_amplitudes*  
    is $N$, one less than the number of configurations to generate
    (default 2).

    *seedname*  
    is the seed name of a previous, successful phonon run.

    The file *seedname*.phonon must exist and be readable. If
    mode_follow is invoked without the .phonon extension, it will also
    attempt to read *seedname*.cell file if it exists and will copy most
    other cell parameters and settings to its output .cell files.

    As an alternative to specifying the arguments on the command line
    mode_follow will also attempt to read them from a file named
    *seedname*.mode-param if it exists. This file should contain a
    Fortran namelist named freeze, whose entries are the identical to
    the command-line argument names. For example

    >     &freeze
    >     mode=4
    >     num_amplitudes=5,
    >     /

    To produce a frozen phonon configuration for a nonzero
    ${\mathbf{q}}$-vector it is also necessary to generate a supercell
    which is commensurate with a frozen phonon at wavevector
    ${\mathbf{q}}$. This supercell may be specified in the .cell file
    using the usual phonon_supercell_matrix block, or by the entry
    SUPERCELL in namelist freeze in the .mode-param file. (There is no
    corresponding command line argument). For example a .mode-param file
    requesting a zone-boundary phonon might contain

    >     &freeze
    >     mode=4
    >     num_amplitudes=5,
    >     qpoint=0,0.5,0
    >     supercell=1,0,0,0,2,0,0,0,1
    >     /

2.  Mode_follow’s second mode of operation is to generate output files
    with the structure perturbed by a frozen phonon at the same
    amplitude but a progressive sequence of phases, which could be used
    for an animation of the mode. In that case the num_amplitudes
    argument should be omitted, and the alternative nframes argument
    given instead (either on the command line or in the .mode-param
    file). This will generate a sequence of $N$ frames with phases
    separated by $2 \pi/N$. Otherwise the arguments and behaviour are
    identical.

    One of the scripts in the cteprouts package, *e.g.* cell2xtl,
    cell2pdb, cell2xyz may be then be used to convert the .cell files to
    a form suitable for visualisation.

#### phonons

The program phonons, one of the CASTEP tools suite is a general purpose
phonon post-processing tool. It can read all of the dynamical matrix or
force constant matrix data from a .check or .castep_bin file generated
in any phonon calculation and re-generate the final phonon output with
changes to one or several “finalisation” options, without needing to
repeat the expensive “electronic” DFPT or supercell parts of a lattice
dynamics calculation. For example, an acoustic sum-rule correction may
be applied to a calculation where this was not chosen initially.

phonons is invoked exactly as is CASTEP and reads a .cell and a .param
file exactly as CASTEP does. These may be identical to or near copies of
the original calculation, but the .param file must contain the
continuation keyword which must point to the .check or .castep_bin file
which contains the dynamical matrix data. If the run is successful it
will write a log file with the extension .phonon_out and a new .phonon
file. In this respect a run of phonons on a continuation deck is very
similar to re-running castep on the same deck. However it will not
attempt to perform any “electronic” calculation and will ignore any
attempt to try. For example, if the parameter elec_energy_tol was
changed castep would discard the saved dynamical matrix data and restart
from the beginning. phonons will ignore this and process the saved data
as a continuation run.

This post-processing can be used for a number of tasks, including

- Turning on or off or modifying an acoustic sum-rule correction by
  changing the parameters keywords phonon_sum_rule or
  phonon_sum_rule_method.

- turning on or off the inclusion of LO/TO splitting terms by changing
  parameters keyword phonon_calc_lo_to_splitting or changing the set of
  directions for ${\mathbf{q}}\rightarrow 0$ by adding or changing cell
  `%block PHONON_GAMMA_DIRECTIONS`.

- adding or omitting the low-frequency ionic term in an E-field
  calculation by changing parameter efield_calc_ion_permittivity. Note
  that an attempt to turn this on will only succeed if there result of a
  $\Gamma$-point phonon calculation is already stored in the .check
  file. phonons will not attempt the electronic calculation necessary to
  compute this if it is not.

- changing the set of fine phonon k-points used as the target of a
  supercell or interpolation calculation by adding or changing
  `%block PHONON_FINE_KPOINT_LIST`, `%block PHONON_FINE_KPOINT_PATH`,
  phonon_fine_kpoint_mp_grid or similar in the .cell file. This permits
  the calculation of *both* a set of dispersion curves *and* a DOS from
  the same electronic run (DFPT/Interpolation or supercell).

- taking the result of a phonon calculation on a Monkhorst Pack grid of
  standard (not fine) phonon kpoints and performing Fourier
  interpolation as a post-processing step.

- switching interpolation methods between spherical and anisotropic
  schemes.

All of the above could also be performed using castep rather than
phonons provided care is taken not to change any parameters which
control the properties of the “electronic” part of the calculation.
However phonons can also perform some additional processing which castep
can not, most notably isotopically substituted lattice dynamics
calculations.

### Running Large Calculations

Phonon calculations even on small crystalline systems typically require
many times the CPU resources of a ground state calculation. DFPT
calculations of phonon dispersion compute dynamical matrices at a number
of phonon wavevectors, each of which contains calculations of several
perturbations. Each perturbation will typically require a large k-point
set due to symmetry breaking by the perturbation. If the supercell
method is used, converged calculations require a system of a typical
size of a few hundred atoms, and many perturbations, although the
k-point set used is smaller. Consequently, calculations on systems of
scientific interest frequently require departmental, university or
national-level supercomputing facilities, usually parallel cluster class
machines.

Much of the advice for effective use of cluster or supercomputer class
resources is the same as for ground-state or other types of CASTEP
calculations, but there are a few special considerations for phonon
calculations, set out below. Among the particularly relevant general
items are the choice of memory/speed tradeoff; usually the best approach
is to select the highest speed option `opt_strategy_bias : 3`[^16] which
retains wavefunction coefficients in memory rather than paging to disk.
For large calculations it is very frequently the case that that the
memory requirement (in particular of the wavefunctions) is the most
important consideration in choosing a parallel distribution. If a run
fails due to exceeding the available memory per node, the processor
count requested should be increased to distribute the wavefunction
arrays across a larger set of processors, reducing the memory/processor
requirement.

If increasing the degree of parallel distribution is not possible,
opt_strategy_bias can be reduced to 0[^17], which will page
wavefunctions to disk. In that case it is vital to ensure that the
temporary scratch files are written to high-speed disk (either local or
a high-performance filesystem). This is usually controlled by setting
the environment variable CASTEP_PAGE_TMPDIR to point to a directory on
an appropriate filesystem[^18].

#### Parallel execution

CASTEP implements a parallel strategy based on a hierarchical
distribution of wavefunction data by k-points, plane-waves, bands and
OpenMP across processors[^19]. In a phonon calculation this is used to
speed up the execution within each perturbation and q-point which are
still executed serially in sequence. Normally an efficient distribution
is chosen automatically providing that the data_distribution parameter
is not changed from the default value MIXED.

K-points, plane-wave, band and task farm parallelism are all implemented
using the [Message-Passing
Interface](https://hpc-tutorials.llnl.gov/mpi) (MPI) system, and CASTEP
must usually be launched by starting the executable using the mpiexec or
similar commands, viz

> mpiexec -n 1024 castep.mpi \<seedname\>

which will start 1024 MPI processes and distribute the calculation
across them.

By contrast [OpenMP](https://www.openmp.org/) parallelism is is
activated by setting the environment variable
`CASTEP_NUM_THREADS : <n>`, which may be done before the mpiexec command
to activate hybrid MPI+OpenMP mode.

##### k-point and g-vector parallelism

To best exploit the k-point component of the parallel distribution, the
total number of the processors requested should be a multiple of the
number of electronic k-points or have a large common divisor. The
parallel distribution is printed to the .castep file, where in this
example four k-points are used:

>     Calculation parallelised over   32 nodes.
>     K-points are distributed over    4 groups, each containing    8 nodes.

For non-phonon CASTEP calculations it is sufficient to choose a
processor count which is a multiple of $N_{\text{k}}$, and that the
degree of plane-wave-parallelism is not so large that efficiency is
lost. However the choice of processor number in a phonon calculation is
severely complicated by the fact that the number of electronic k-points
in the irreducible Brillouin Zone changes during the run as the
perturbations and phonon q-points have different symmetries. It is not
convenient to compute the number for any perturbation individually
without a detailed analysis, and some compromise between all of the
perturbations should be chosen. To assist in this choice a utility
program, phonon_kpoints is provided. This reads the configuration of the
proposed calculation from the .cell file and is simply invoked by

> phonon_kpoints *seedname*

It then determines and prints the k-point counts, and provides a “figure
of merit” for a range of possible processor counts. On most parallel
architectures the efficiency of the plane-wave parallelism becomes
unacceptable if there are fewer than around 200 plane-waves per node. It
is usually possible to choose a processor count which allows a highly
parallel run while keeping the number of plane-waves per node
considerably higher than this.

##### band parallelism

From CASTEP version 24.1 band parallelism is implemented for FD
calculations, but not yet DFPT. This may be enabled using a “devel-code”
string in the .param file[^20]

>     %block DEVEL_CODE
>     PARALLEL:NBANDS=8:ENDPARALLEL
>     %endblock DEVEL_CODE

will attempt to set up 8-way band parallelism in addition to k-point and
g-vector.

##### perturbation parallelism

From CASTEP release 24.1, a further level of parallelism called
task-farming may be used to distribute perturbations across processors.
This is set up and run in the same manner as for PIMD or NEB
calculations by setting parameters file keyword `num_farms : <n>`. The
value of $n$ should not be too large, as the performance will be limited
by load-balance issues, and in any case never greater than the number of
perturbations. This will produce $n$ output files named

> \<seedname\>\_farm00\<n\>.castep

instead of the usual, single .castep file. Of these,
\<seedname\>\_farm001.castep will contain the calculated final
frequencies, Born effective charges etc. but a single instance of the
.phonon, .efield, .phonon_dos and .check files are written as with other
parallel schemes.

##### hybrid OpenMP parallelism

In addition to the above-described parallel distribution strategies -
all based on MPI parallelism, CASTEP also offers a degree of OpenMP
parallelism[^21]. This can speed up some operations such as matrix
diagonalisations and is activated by setting the environment variable
`CASTEP_NUM_THREADS : <n>` where $n$ in the range 2-16 is most effective
(only a modest speed-up is accessible using OpenMP). However this can
give a useful gain in addition to MPI parallelism especially in the case
where compute nodes must be underpopulated because of memory
requirements.

#### Checkpointing and Restarting

Even with a parallel computer, it is frequently the case that a
calculation can not be completed in a single run. Many machines have a
maximum time limit on a batch queue which may be too short. On desktop
machines, run time may be limited by reliability and uptime limitations.
CASTEP is capable of periodically writing “checkpoint” files containing
a complete record of the state of the calculation and of restarting and
completing a calculation from such a checkpoint file. In particular
dynamical matrices from complete q-points, and partial dynamical
matrices from each perturbation are saved and can be used in a restart
calculation. To enable the writing of periodic checkpoint files, set the
parameter

> backup_interval 3600

which will write a checkpoint file named *seedname*.check every hour
(the time is specified in seconds) or on completion of the next
perturbation thereafter. To restart a calculation, set the parameter

> `continuation : default`

in the .param file before resubmitting the job. This will attempt to
read *seedname*.check and restart the calculation from there.
Alternatively the full filename of a checkpoint file may be given as
argument to the continuation keyword to read an explicitly named file.

At the end of the calculation a checkpoint file *seedname*.check is
always written. As with the intermediate checkpoint files this contains
a (now complete) record of the dynamical matrices or force constant
matrix resulting from phonon calculation. This may be analysed in a
post-processing phase using the phonons utility.

### Advanced Topics

#### Molecular vibrational calculations

CASTEP is primarily a solid state code with periodic boundary
conditions, and is not necessarily the first choice for performing
vibrational spectroscopy calculations on molecules. Nevertheless it can
sometimes be convenient to use it in this mode (for example, if it is
desired to compare a molecular crystal with an isolated molecule at
exactly the same level of theory).

With certain limitations this can be done in CASTEP using the “molecule
in large box” method. The idea is to place the molecule in a vacuum by
constructing a large supercell. If the supercell is sufficiently large
that the periodic images of the molecule do not interact then the
frequencies and eigenvectors in the molecular limit can be recovered.
Because the zero charge density in a large volume must still be
represented using a plane-wave basis sat, such calculations can become
expensive, but it is nevertheless essential to perform careful cell size
convergence tests.

A molecule in free space possesses rotational as well as translational
symmetry, which results in three (or two for a linear molecule) free
librational modes of zero frequency. A periodic CASTEP calculation does
not possess this symmetry and the librational mode frequencies will have
nonzero values, either real and positive or imaginary. If a molecular
electric field calculation is desired the parameter

> efield_ignore_molec_modes

may be specified with the keyword value of MOLECULE or LINEAR_MOLECULE.
This omits the 6 (or 5) lowest modes from the infra-red polarizability
calculation, assuming them to be the librations.

As with crystalline systems it is desirable to maximise the use of
symmetry by optimally orienting the molecule. In this case an additional
criterion applies - the simulation cell vectors should be organised for
maximal compatibility with the molecular symmetry. Any incompatible
symmetry operations will not be found if using the SYMMETRY_GENERATE
keyword. For example a tetrahedral molecule is best modelled in a cubic
cell with the 3-fold axes oriented along the cube diagonals. By contrast
a hexagonal molecule such as benzene is better modelled using a
hexagonal simulation cell with the molecule commensurately oriented.

A limitation of this method is that CASTEP implements only the
crystallographic point group symmetry operations, and does not include
all of the molecular point groups. It is not always possible to take
full advantage of molecular symmetry, and 5-fold rotation axes and
icosahedral groups can not be represented. In such cases some
degeneracies will be lifted and only approximately satisfied, and the
group theoretical analysis of the eigenvectors will not be correct for
the molecular point group.

Finally, it is sometimes desirable to compute the vibrational spectrum
of an ion rather than a neutral molecule. This adds a number of
complications of different degrees of severity. First is the well-known
result that semi-local DFT (LDA and GGA) severely underbinds anions, and
the description of the molecular states in some systems may be
substantially incorrect. Second, the introduction of a non-zero
molecular charge in a periodic cell necessitates the addition of a
compensating charge to avoid a divergent Coulomb energy. CASTEP
implicitly adds a uniform charge distribution which integrates over the
cell to the negative of the sum of ionic and electronic charge in the
system. This model while effective, gives rise to additional
electrostatic terms which decrease only as $1/L^{3}$. Consequently it is
impossible to completely converge such a calculation with respect to
cell size. It is beyond the scope of this document to describe the
techniques which may be used to recover the “infinite” volume limit. The
reader is referred to ref (Leslie and Gillan 1985) and the extensive
literature which cites this paper for further reading. An example of a
successful approach involving an extrapolation of frequencies of a
molecular anion to the infinite cell limit is given in Parker *et
al.* (Parker et al. 2006).

#### Isotopic substitution calculations

A powerful technique in experimental spectroscopy is *isotopic
substitution* where the sample is modified by the substitution of a
common isotope of some element by a heavier or lighter one. This can
provide valuable site-specific spectroscopic information, particularly
if the substitution can be performed in a site-specific manner. Force
constants and dynamical matrices do not depend on nuclear mass, and are
therefore one set of computations will suffice irrespective of isotopic
substitution. The phonons tool
(section <a href="#sec:phonons-tool" data-reference-type="ref"
data-reference="sec:phonons-tool">5.4</a>) allows recomputation of the
frequencies for different isotopic substitutions without recomputing
dynamical matrices.

The simplest case is complete substitution. The new .cell file should
contain a block specifying the new mass, *e.g.*

>      %block SPECIES_MASS
>        H     2.014101778
>      %endblock SPECIES_MASS

Running the phonon code will compute the vibrational frequencies and
eigenvectors for a completely deuterated sample.

Slightly more complicated is the case where only one or a few of a
number of sites is to be substituted. In that case the above
prescription is inadequate as it would change all of the sites
containing the species in question. It is therefore necessary to use
CASTEP’s sub-species labelling capabilities. For example the original
cell containing a methane molecule might contain

    %block SPECIES_POT
    C C_00.recpot
    H H_04.recpot
    %endblock SPECIES_POT

    %block POSITIONS_ABS
    C    0.0         0.0          0.0
    H    0.0         0.0          1.09
    H    0.9316592   0.0         -0.36333333
    H   -0.4658296   0.8068405   -0.36333333
    H   -0.4658296  -0.8068405   -0.36333333
    %endblock POSITIONS_ABS

The continuation cell for input to phonons should contain

    %block SPECIES_POT
    C C_00.recpot
    H H_04.recpot
    H:D H_04.recpot
    %endblock SPECIES_POT

    %block POSITIONS_ABS
    C    0.0         0.0          0.0
    H:D  0.0         0.0          1.09
    H    0.9316592   0.0         -0.36333333
    H   -0.4658296   0.8068405   -0.36333333
    H   -0.4658296  -0.8068405   -0.36333333
    %endblock POSITIONS_ABS

    %block SPECIES_MASS
       H:D     2.014101778
       H       1.007825032 
       C       12.0107 
    %endblock SPECIES_MASS

The new atom label may contain any alphanumeric extension following the
colon up to a maximum of 8 characters.

This defines a new system with one hydrogen atom of a methane molecule
replaced by deuterium. Notice that this has broken the initial symmetry
of the cell. The phonons program generates a new system of reduced
symmetry with three distinct atom types, and copies the dynamical matrix
data from the original methane system from the .check file over to the
corresponding atoms in the new system. It then diagonalises the
dynamical matrix, applies whatever post-processing is specified and
writes the frequencies and eigenvectors as usual.

This straightforward approach fails to take into account isotopic shifts
in bond lengths and geometry, and is therefore an approximate one. Since
isotopic shifts in bond lengths depend in turn on the vibrational
frequency, within the quasiharmionic approximation a self-consistent
approach under which bond lengths and frequencies are simultaneously
adjusted to self-consistency would be required. This would also require
new DFPT electronic calculations at each stage and is beyond the scope
of this utility.

#### Constrained lattice dynamics

It is sometimes desirable to compute the vibrational frequencies of a
restricted region of an *ab initio* model. Consider the case of a
molecule adsorbed on a surface. If only the frequencies of the molecular
vibrations are needed it would be desirable to compute a subset of the
modes were this possible. While it is not practical (as it would require
prior knowledge of the eigenvectors of the full calculation), an
alternative approach is to apply constraints to certain atoms.

CASTEP implements a technique known as *constrained lattice dynamics*,
also known as the *partial hessian* method. Nominated atoms are assumed
to be “frozen”, and the corresponding entries of the dynamical matrix
are set to zero. The model is effectively one whereby the atoms in the
region of the system deemed “irrelevant” are assigned a mass of
infinity. It is not necessary to perform any computations for
perturbation of these atoms, and a considerable saving of computational
effort may be achieved.

Frozen atoms are specified a CASTEP .cell file using the same syntax as
applies to geometry operations. The block

    %block IONIC_CONSTRAINTS
    1 Si  1   1   1   1
    2 Si  2   1   1   1
    3 C   4   1   1   1
    %endblock IONIC_CONSTRAINTS 

constrains silicon atoms numbered 1-2 and C atom number 4 to be fixed.
These do not move during geometry optimisation, and their perturbations
are not considered during a lattice dynamics calculation. In fact the
constrained lattice dynamics method does not make full use of the
generality of CASTEP’s linear constraints block (see for example the
tutorial on MD at <http://www.castep.org>) but only identifies atoms
which are fully constrained not to move. As in the example above, there
should be a single line for each atom which creates a uniquely numbered
constraint. This should contain a “1” in all of the x, y, z positions.

Except in rather specialised geometries the presence of fixed atom
constraints is incompatible with most symmetry operations, and therefore
symmetry should usually be turned off during a constrained lattice
dynamics calculation. There is also an incompatibility with the acoustic
sum rule (section <a href="#sec:asr" data-reference-type="ref"
data-reference="sec:asr">2.2.4</a>) as constraining the atoms breaks the
translational invariance of the Hamiltonian. Acoustic sum rule
correction is therefore disabled automatically if constraints are
present.

### Keyword Reference

The cell and param file keywords used specifically for phonon and
related calculations are listed here alphabetically with brief
descriptions.

#### CELL file keywords

PHONON_FINE_KPOINT_LIST  
*List of phonon wavevectors on the fine grid* (Block)  
Phonon frequencies are calculated on a coarse set of wavevectors using
DFPT and interpolated onto this finer list of points.

PHONON_FINE_KPOINT_MP_GRID  
*Fine MP-grid of phonon wavevectors* (Integer Vector)  
Phonon frequencies are calculated on a coarse set of wavevectors using
DFPT and interpolated onto this finer grid of wavevectors.

PHONON_FINE_KPOINT_MP_OFFSET  
*Origin offset of the fine phonon MP grid* (Real Vector)  
The offset of the fine MP grid at which the phonons calculated using
DFPT are interpolated.

PHONON_FINE_KPOINT_MP_SPACING  
*The spacing of points on the fine MP set for phonons* (Physical)  
This specifies the minimum spacing between points on a Monkhorst-Pack
grid that phonons will be interpolated onto from the coarser phonon
grid.

PHONON_FINE_KPOINT_PATH  
*Path of phonon wavevectors on a fine scale* (Block)  
Phonon frequencies are calculated on a coarse set of wavevectors using
DFPT and interpolated onto this finer path.

PHONON_FINE_KPOINT_PATH_SPACING  
*The fine spacing of points on a path at which phonons are calculated*
(Physical)  
The spacing of k-points along a path (specified by
PHONON_FINE_KPOINT_PATH) at which phonons will be interpolated from a
coarser grid

PHONON_GAMMA_DIRECTIONS  
*Phonon gamma-point LO/TO splitting* (Block)  
This is a list of directions along which ${\mathbf{q}}\rightarrow 0$
will be calculated for the non-analytic LO/TO term in a phonon
calculation at ${\mathbf{q}}=0$. Fractional coordinates must be used.
Default value: The k-point before gamma in the k-point list, or the one
after or (0.1, 0, 0)

PHONON_KPOINTS_LIST  
*Alias for* PHONON_KPOINT_LIST (Block)  
Default value: Determined from PHONON_KPOINT_LIST

PHONON_KPOINTS_PATH  
*Alias for* PHONON_KPOINT_PATH (Block)  
Default value: Determined from PHONON_KPOINT_PATH

PHONON_KPOINTS_PATH_SPACING  
*Alias for* PHONON_KPOINT_PATH_SPACING (Physical)  
Default value: Determined from PHONON_KPOINT_PATH_SPACING

PHONON_KPOINT_LIST  
*List of phonon k-points* (Block)  
A list of discrete k-points at which phonon frequencies and eigenvectors
will be calculated. Default value: SCF k-points are used if an
alternative PHONON_KPOINT\_\* specifier is not given.

PHONON_KPOINT_MP_GRID  
*Phonon wavevector Monkhorst-Pack grid* (Integer Vector)  
The phonon wavevectors defined by a Monkhorst-Pack grid. Symmetry (if
specified) will be used to generate the wavevector list and weights.

PHONON_KPOINT_MP_OFFSET  
*Phonon wavevector Monkhorst Pack grid offset* (Real Vector)  
The offset of the origin of the Monkhorst-Pack set for phonons in
fractional coordinates, or the keyword INCLUDE_GAMMA. Default value:
INCLUDE_GAMMA

PHONON_KPOINT_MP_SPACING  
*Phonon wavevector Monkhorst-Pack grid density* (Physical)  
The density of wavevectors on a a Monkhorst-Pack grid for phonon
calculations. Units of inverse length should be specified. Default
value: 0.1 A$^{-1}$.

PHONON_KPOINT_PATH  
*Phonon dispersion k-point path* (Block)  
The path continuous through the BZ on which phonon dispersion is
calculated. This is specified in fractional coordinates. Default value:
None.

PHONON_KPOINT_PATH_SPACING  
*Phonon dispersion path spacing* (Physical)  
The maximum spacing between kpoints along the path specified by
PHONON_KPOINT_PATH. Units of inverse length must be specified. Default
value: 0.1 A$^{-1}$.

PHONON_SUPERCELL_MATRIX  
*Supercell matrix for finite difference phonon calculations* (Block)  
The supercelling matrix for force constant matrix calculations. The
supercell matrix is specified by a 3x3 integer matrix which gives the
supercell used in finite-difference phonon calculations.

SUPERCELL_KPOINTS_LIST  
*SCF k-points for FD phonon supercell* (Block)  
A list of k-points in the Brillouin zone (with associated weights) used
for BZ integration during a supercell FD phonon calculation. The k-point
weights must sum to 1. Default value: Generated from
SUPERCELL_KPOINTS_MP_SPACING and the crystal symmetry.

SUPERCELL_KPOINTS_MP_GRID  
*SCF Monkhorst-Pack grid for FD phonon supercell calculation* (Integer
Vector)  
The k-points defined by a Monkhorst-Pack grid when doing a finite
displacement phonon calculation. Symmetry (if specified) will be used to
generate the k-point list and weights. Default value: Generated from
SUPERCELL_KPOINTS_MP_SPACING.

SUPERCELL_KPOINTS_MP_OFFSET  
*SCF Monkhorst Pack grid offset for a FD phonon supercell calculation*
(Real Vector)  
The offset of the origin of the Monkhorst-Pack set in fractional
coordinates when performing a finite displacement phonon calculation.
Default value: (0, 0, 0).

SUPERCELL_KPOINTS_MP_SPACING  
*SCF Monkhorst-Pack grid density for FD phonon supercell calculation*
(Physical)  
The k-point density of a Monkhorst-Pack grid for a supercell FD phonon
calculation. Units of inverse length should be specified. Default value:
0.1 A$^{-1}$.

#### PARAM file keywords

NUM_BACKUP_ITER  
*md/geom iterations between backups* (Integer)  
Specifies the number of iterations between backups of all data for
restarts, for a geometry optimization or molecular dynamics run. Allowed
values: (any integer) $>$ 0 Default value : 5

BACKUP_INTERVAL  
*seconds between backups* (Integer)  
Specifies the interval, in seconds, between backups of all data for
restarts, for a geometry optimization/molecular dynamics/phonon run - if
less than or equal to zero then no timed backups. Allowed values: (any)
Default value : 0

BORN_CHARGE_SUM_RULE  
*enforce Born charge sum rule* (Logical)  
Selects whether to explicitly correct the Born effective charge tensor
to enforce the sum rule that effective charges sum to zero.

CALCULATE_BORN_CHARGES  
*calculate Born effective charges* (Logical)  
Selects whether to compute Born effective charge tensors as part of a
phonon or E-field linear-response calculation. Allowed values: TRUE or
FALSE Default value : TRUE

CALCULATE_RAMAN  
*calculate Raman intensities* (Logical)  
Selects whether to compute Raman intensities as part of a phonon or
E-field linear-response calculation. Allowed values: TRUE or FALSE
Default value : FALSE

EFIELD_CALC_ION_PERMITTIVITY  
*calculate zero-frequency permittivity* (Logical)  
Specifies whether or not to compute the zero-frequency dielectric
permittivity based on the ionic response to electric fields. This
requires a gamma-point phonon calculation in addition to the EFIELD
linear response one. Allowed values: TRUE or FALSE Default value : TRUE

EFIELD_CALCULATE_NONLINEAR  
*calculate non-linear optical susceptibility* (String)  
Select which non-linear optical property to calculate during TASK=EFIELD
calculation. Allowed values: NONE, CHI2 Default value : NONE

EFIELD_CONVERGENCE_WIN  
*convergence tolerance window in EFIELD* (Integer)  
The LR convergence criteria must be met for EFIELD_CONVERGENCE_WIN
iterations before acceptance. Allowed values: (any integer) $\ge$ 2
Default value : 2

EFIELD_DFPT_METHOD  
*E-field DFPT solver method* (String)  
Selects the solver for E-field density functional perturbation theory..
Allowed values: ALLBANDS(=VARIATIONAL) or DM(=GREEN) to select Gonze
variational or Baroni Green function with DM solver. Default value :
ALLBANDS

EFIELD_ENERGY_TOL  
*E(2) convergence tolerance in EFIELD* (Physical)  
Tolerance for accepting convergence of the field constants during PHONON
calculation. The difference between max and min E(2) values over
EFIELD_CONVERGENCE_WIN iterations must be less than this. NB This is an
INTENSIVE parameter and has units of volume. Allowed values: (any) $>$
0.0 Default value : $10^{-5}$ A$^3$

EFIELD_FREQ_SPACING  
*Spacing of frequencies in permittivity calculation* (Physical)  
Spacing of frequencies in calculation of frequency-dependent
permittivity. Allowed values: (any) $>$ 0.0 Default value : 1.0
cm$^{-1}$

EFIELD_IGNORE_MOLEC_MODES  
*Ignore lowest modes in permittivity calculation* (String)  
Ignore the lowest lying (3,5,6) modes when computing the ionic
contribution to the permittivity and polarizability. Allowed values:
CRYSTAL, MOLECULE, LINEAR_MOLECULE Default value : CRYSTAL

EFIELD_MAX_CG_STEPS  
*max. number of cg steps in EFIELD* (Integer)  
The maximum number of conjugate gradient steps in EFIELD calculation
before performing a SD reset. Allowed values: (any integer) $\ge$ 0
Default value : 0

EFIELD_MAX_CYCLES  
*maximum cycles in EFIELD* (Integer)  
The maximum number of SCF cycles in EFIELD calculation regardless of
convergence. Allowed values: (any integer) $>$ 0 Default value : 50

EFIELD_OSCILLATOR_Q  
*Q-factor for line-shape broadening* (Real)  
Oscillator Q-factor for line-shape broadening in calculation of
frequency- dependent permittivity. Allowed values: (any) $>$ 0.0 Default
value : 50.0

EFIELD_UNIT  
*unit of electric field in output* (String)  
Controls the units used for all electric field in output - many
different units are supported. Default value : eV/A/E

ELEC_METHOD  
*treatment of metals or finite temperature insulator* (String)  
The treatment of metals or finite temperature insulator to be used. An
alias for METALS_METHOD. Allowed values: NONE (=ALLBANDS), DM, EDFT
Default value : DM

EXCITED_STATE_SCISSORS  
*“scissors” operator band-gap correction* (Physical)  
Effectively adds an offset to conduction-band eigenvalues as empirical
correction for LDA/GGA underestimation of band-gaps. Allowed values:
(any) Default value : 0.0

FIX_OCCUPANCY  
*treat system as an insulator* (Logical)  
Determines if the system is treated as an insulator or a metal. Allowed
values: TRUE or FALSE Default value : FALSE

GEOM_FORCE_TOL  
*geometry optimization force convergence tolerance* (Physical)  
Tolerance for accepting convergence of the maximum \|ionic force\|
during geometry optimization. Allowed values: (any) $>$ 0.0 Default
value : 0.05 eV/A

PHONON_CALCULATE_DOS  
*density of states calculation* (Logical)  
Determines whether or not the phonon density of states will be
calculated. Allowed values: TRUE or FALSE Default value : FALSE

PHONON_DOS_SPACING  
*density of states calculation* (Physical)  
The resolution at which a phonon density-of-states will be calculated.
Allowed values: (any) $> 0.0$ Default value : 10.0 cm$^{-1}$.

PHONON_DOS_LIMIT  
*density of states calculation* (Logical)  
The largest phonon to be included in a phonon density-of-states
calculation. Allowed values: (any) $>$ PHONON_DOS_SPACING Default value
: 5000.0 cm$^{-1}$.

PHONON_CALC_LO_TO_SPLITTING  
*gamma-point phonon LO/TO correction* (Logical)  
Selects whether to compute non-analytic contribution to dynamical matrix
from long-ranged electric field effects responsible for LO/TO splitting.
This requires calculation of the dielectric permittivity by E-field
linear-response and the Born effective charges. Allowed values: TRUE or
FALSE Default value : TRUE

PHONON_CONVERGENCE_WIN  
*convergence tolerance window in LR* (Integer)  
The LR convergence criteria must be met for PHONON_CONVERGENCE_WIN
iterations before acceptance. Allowed values: (any integer) $\ge$ 2
Default value : 2

PHONON_ENERGY_TOL  
*E(2) convergence tolerance in LR* (Physical)  
Tolerance for accepting convergence of the force constants during PHONON
calculation. The difference between max and min E(2) values over
PHONON_CONVERGENCE_WIN iterations must be less than this. Allowed
values: (any) $>$ 0.0 Default value : $10^{-5}$ eV/A$^2$

PHONON_FINE_METHOD  
*fine phonon calculation method* (String)  
Selects which calculation method to use for phonon calculation on a fine
grid. Allowed values: NONE, SUPERCELL, INTERPOLATE Default value :
SUPERCELL if TASK=THERMODYNAMICS else NONE

PHONON_FINITE_DISP  
*finite displacement amplitude* (Physical)  
The amplitude of the ionic perturbation to be used in a finite
displacement phonon calculation. Allowed values: (any) $>$ 0.0 Default
value : 0.01 $a_0$

PHONON_FORCE_CONSTANT_CUTOFF  
*Cutoff for force constant matrix* (Physical)  
The cutoff for the force constant matrix in a phonon calculation on a
fine grid with supercell method. Allowed values: (any) $\ge$ 0.0 Default
value : 0.0

PHONON_FINE_CUTOFF_METHOD  
*Selects which method to use to extract non-periodic force constant
matrix from periodic supercell.* (String)  
With the CUMULANT method, all contributions from the periodic supercell
are summed with a suitable weighting factor to avoid double counting of
image contributions.  
The SPHERICAL method, uses a minimum image convention with a spherical
cutoff given by PHONON_FORCE_CONSTANT_CUTOFF.  
Allowed values: CUMULANT and SPHERICAL. Default value : CUMULANT.

PHONON_FORCE_CONSTANT_CUT_SCALE  
*Scaling factor for aspherical force constant matrix cutoff* (Real)  
The range of force constant terms included is up to s times halfway to
the Wigner Seitz cell boundary. This parameter supplies the value of s.
Allowed values: 0.0 $\le$ (any) $\ge$ 1.0 Default value : 0.0

PHONON_FORCE_CONSTANT_ELLIPSOID  
*Ellipsoid size for force constant matrix* (Real)  
Alias for PHONON_FORCE_CONSTANT_CUT_SCALE (deprecated).

PHONON_MAX_CG_STEPS  
*max. number of cg steps in LR* (Integer)  
The maximum number of conjugate gradient steps in PHONON calculation
before performing a SD reset. Allowed values: (any integer) $\ge$ 0
Default value : 0

PHONON_MAX_CYCLES  
*maximum cycles in LR* (Integer)  
The maximum number of SCF cycles in PHONON calculation regardless of
convergence. Allowed values: (any integer) $\ge$ 0 or if TASK=PHONON etc
$\ge$ PHONON_CONVERGENCE_WIN Default value : 50

PHONON_METHOD  
*phonon calculation method* (String)  
Selects which calculation method to use for phonons. Allowed values:
DFPT, LINEARRESPONSE, FINITEDISPLACEMENT Default value : set by
PHONON_FINE_METHOD

PHONON_DFPT_METHOD  
*phonon DFPT solver method* (String)  
Selects the solver for phonon density functional perturbation theory..
Allowed values: ALLBANDS(=VARIATIONAL) or DM(=GREEN) to select Gonze
variational or Baroni Green function with DM solver. Default value : DM
if `FIX_OCCUPANCY : FALSE`, otherwise ALLBANDS

PHONON_PRECONDITIONER  
*scheme to use in LR* (String)  
The preconditioning scheme used by the CG minimiser in LR. Allowed
values: RTPA, TPA, PS, NONE Default value : TPA

PHONON_SUM_RULE  
*enforce phonon sum rule* (Logical)  
Selects whether to explicitly correct the dynamical matrix to enforce
the acoustic q=0 phonon sum rule, i.e. that 3 modes have zero frequency
at q=0. Allowed values: TRUE or FALSE Default value : FALSE

PHONON_SUM_RULE_METHOD  
*select method to enforce phonon sum rule* (String)  
Selects a method to use when enforcing acoustic phonon sum rule. Allowed
values: NONE : No sum-rule correction will be applied. RECIPROCAL :
Correct dynamical matrix D(q) at each q using D(q=0). REALSPACE :
Correct the real-space force constant matrix C(R). REAL-RECIP : Correct
C(R) in realspace followed by D(q) in reciprocal space. MOLECULAR :
Correct D(0) using rotational as well as translational sum-rule.
MOLECULAR-1D : Correct D(0) for a linear molecule using rotational as
well as translational sum-rule. Default value : RECIPROCAL

PHONON_USE_KPOINT_SYMMETRY  
*reduced or full kpoint set in LR* (Logical)  
Selects which k-point set to use For each phonon q-vector in LR: T =\>
use the irreducible k-point set of the (reduced) symmetry, F =\> use the
complete fully symmetric k-point set. Allowed values: TRUE or FALSE
Default value : TRUE

PHONON_WRITE_FORCE_CONSTANTS  
*Write out real-space force constant matrix* (Logical)  
Selects whether to write out the real-space force constant matrix from a
phonon supercell or interpolation calculation (to the
\<seedname\>.castep file) for the case of PHONON_FINE_METHOD /= NONE.
Allowed values: TRUE or FALSE Default value : FALSE

PHONON_WRITE_DYNAMICAL  
*Write out reciprocal space dynamical matrix* (Logical)  
Selects whether to write out the reciprocal space dynamical matrices
from a phonon calculation (to the \<seedname\>.castep file). /= NONE.
Allowed values: TRUE or FALSE Default value : FALSE

THERMO_T_NPOINTS  
*Number of points in temperature interval* (Integer)  
The number of points in the temperature interval for the thermodynamics
calculation. Allowed values: (any integer) $\ge$ 1 Default value : 2 if
THERMO_T_STOP /= THERMO_T_START or 1 otherwise

THERMO_T_SPACING  
*Temperature spacing* (Physical)  
The spacing between temperature values for the thermodynamics
calculation. Allowed values: (any) $>$ (-epsilon) Default value :
THERMO_T_STOP - THERMO_T_START

THERMO_T_START  
*Starting temperature* (Physical)  
The desired starting temperature for the thermodynamics calculation.
Allowed values: (any) $>$ 0.0 Default value : 298 K

THERMO_T_STOP  
*Final temperature* (Physical)  
The desired final temperature for the thermodynamics calculation.
Allowed values: (any) $\ge$ THERMO_T_START Default value : 298

### Bibliography

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-AcklandWC97" class="csl-entry">

Ackland, G J, M C Warren, and S J Clark. 1997. “Practical Methods in Ab
Initio Lattice Dynamics.” *Journal of Physics: Condensed Matter* 9 (37):
7861–72. <http://stacks.iop.org/0953-8984/9/7861>.

</div>

<div id="ref-BalanSMC01" class="csl-entry">

Balan, E., A. M. Saitta, F. Mauri, and G. Calas. 2001. “First-Principles
Modeling of the Infra-Red Spectrum of Kaolinite.” *American
Mineralogist* 86: 1321–30.

</div>

<div id="ref-BaroniDDG01" class="csl-entry">

Baroni, S., S. De Gironcoli, A. Dal Corso, and P. Giannozzi. 2001.
“Phonons and Related Crystal Properties from Density- Functional
Perturbation Theory.” *Rev. Mod. Phys.* 73: 515–62.

</div>

<div id="ref-Boyd03" class="csl-entry">

Boyd, Robert W. 2003. *Nonlinear Optics*. San Diego, Calif.: Academic
Press. <http://www.myilibrary.com?id=105034>.

</div>

<div id="ref-ChengRC20" class="csl-entry">

Cheng, Y. Q., and A. J. Ramirez-Cuesta. 2020. “Calculation of the
Thermal Neutron Scattering Cross-Section of Solids Using OCLIMAX.” *J.
Chem. Theory Comput.* 16 (8): 5212–17.
<https://doi.org/10.1021/acs.jctc.0c00569>.

</div>

<div id="ref-Dove93" class="csl-entry">

Dove, M. T. 1993. *Introduction to Lattice Dynamics*. Cambridge
University Press.

</div>

<div id="ref-FairJVJLRP22" class="csl-entry">

Fair, R., A. Jackson, D. Voneshen, D. Jochym, D. Le, K. Refson, and T.
Perring. 2022. “Euphonic: Inelastic Neutron Scattering Simulations from
Force Constants and Visualization Tools for Phonon Properties.” *J Appl
Cryst* 55 (6): 1689–1703. <https://doi.org/10.1107/S1600576722009256>.

</div>

<div id="ref-GianozziGPB91" class="csl-entry">

Giannozzi, Paolo, Stefano de Gironcoli, Pasquale Pavone, and Stefano
Baroni. 1991. “Ab Initio Calculation of Phonon Dispersions in
Semiconductors.” *Phys. Rev. B* 43 (9): 7231–42.
<https://doi.org/10.1103/PhysRevB.43.7231>.

</div>

<div id="ref-GonzeL97" class="csl-entry">

Gonze, X., and C. Lee. 1997. “Dynamical Matrices, Born Effective
Charges, Dielectric Permittivity Tensors, and Interatomic Force
Constants from\> Density-Functional Perturbation Theory.” *Phys. Rev. B*
55: 10355–68.

</div>

<div id="ref-KendrickB16" class="csl-entry">

Kendrick, John, and Andrew D. Burnett. 2016. “PDielec: The Calculation
of Infrared and Terahertz Absorption for Powdered Crystals.” *Journal of
Computational Chemistry* 37 (16): 1491–1504.
<https://doi.org/10.1002/jcc.24344>.

</div>

<div id="ref-KuncM82" class="csl-entry">

Kunc, K., and Richard M. Martin. 1982. “Ab Initio Force Constants of
GaAs: A New Approach to Calculation of Phonons and Dielectric
Properties.” *Phys. Rev. Lett.* 48 (6): 406–9.
<https://doi.org/10.1103/PhysRevLett.48.406>.

</div>

<div id="ref-leslie:85" class="csl-entry">

Leslie, M., and M. J. Gillan. 1985. “The Energy and Elastic Dipole
Tensor of Defects in Ionic Crystals Calculated by the Supercell Method.”
*Journal of Physics: Condensed Matter* 18: 973–82.

</div>

<div id="ref-LevineA89" class="csl-entry">

Levine, Zachary H., and Douglas C. Allan. 1989. “Linear Optical Response
in Silicon and Germanium Including Self-Energy Effects.” *Phys. Rev.
Lett.* 63 (16): 1719–22. <https://doi.org/10.1103/PhysRevLett.63.1719>.

</div>

<div id="ref-LWM15" class="csl-entry">

Lloyd-Williams, Jonathan H., and Bartomeu Monserrat. 2015. “Lattice
Dynamics and Electron-Phonon Coupling Calculations Using Nondiagonal
Supercells.” *Physical Review B* 92 (18): 184301.
<https://doi.org/10.1103/PhysRevB.92.184301>.

</div>

<div id="ref-Maradudin80" class="csl-entry">

Maradudin, A. A., and G. K. Horton, eds. 1980. *Dynamical Properties of
Solids*. Vol. 1–7. Elsevier.

</div>

<div id="ref-Miwa2011" class="csl-entry">

Miwa, Kazutoshi. 2011. “Prediction of Raman Spectra with Ultrasoft
Pseudopotentials.” *Physical Review B* 84 (9): 094304.
<https://doi.org/10.1103/PhysRevB.84.094304>.

</div>

<div id="ref-ParkerRWBHY06" class="csl-entry">

Parker, S. F., K. Refson, K. P. J. Williams, D. A. Braden, B. S. Hudson,
and K. Yvon. 2006. “Spectroscopic and Ab Initio Characterization of the
ReH9(2-) Ion.” *Inorg. Chem.* 45 (26): 10951–57.

</div>

<div id="ref-ParlinskiLK97" class="csl-entry">

Parlinski, K., Z. Q. Li, and Y. Kawazoe. 1997. “First-Principles
Determination of the Soft Mode in Cubic $ZrO2$.” *Phys. Rev. Lett.* 78
(21): 4063–66. <https://doi.org/10.1103/PhysRevLett.78.4063>.

</div>

<div id="ref-Pickard97" class="csl-entry">

Pickard, Chris J. 1997. “Ab Initio Electron Energy Loss Spectroscopy.”
PhD Thesis, University of Cambridge.

</div>

<div id="ref-PorezagP:96" class="csl-entry">

Porezag, Dirk, and Mark R. Pederson. 1996. “Infrared Intensities and
Raman-Scattering Activities Within Density-Functional Theory.” *Phys.
Rev. B* 54 (11): 7830–36. <https://doi.org/10.1103/PhysRevB.54.7830>.

</div>

<div id="ref-ProbertP03" class="csl-entry">

Probert, M. I. J., and M. C. Payne. 2003. “Improving the Convergence of
Defect Calculations in Supercells: An Ab Initio Study of the Neutral
Silicon Vacancy.” *Phys. Rev. B* 67 (7): 075204.
<https://doi.org/10.1103/PhysRevB.67.075204>.

</div>

<div id="ref-RamirezCuesta04" class="csl-entry">

Ramirez-Cuesta, A. J. 2004. “<span class="nocase">aCLIMAX</span> 4.0.1,
the New Version of the Software for Analyzing and Interpreting INS
Spectra.” *Computer Physics Communications* 157 (3): 226–38.
[https://doi.org/DOI:
10.1016/S0010-4655(03)00520-4](https://doi.org/DOI: 10.1016/S0010-4655(03)00520-4).

</div>

<div id="ref-RefsonTC06" class="csl-entry">

Refson, K, P. R. Tulip, and S. J. Clark. 2006. “Variational
Density-Functional Perturbation Theory for Dielectrics and Lattice
Dynamics.” *Physical Review B* 73: 155114.

</div>

<div id="ref-Srivastava90" class="csl-entry">

Srivastava, G. P. 1990. *The Physics of Phonons*. Institute of Physics
Publishing.

</div>

<div id="ref-YatesWVS07" class="csl-entry">

Yates, Jonathan R, Xinjie Wang, David Vanderbilt, and Ivo Souza. 2007.
“Spectral and Fermi Surface Properties from Wannier Interpolation.”
*Physical Review B* 75 (19): 195121.
<https://doi.org/10.1103/PhysRevB.75.195121>.

</div>

<div id="ref-YeLWH04" class="csl-entry">

Ye, Lin-Hui, Bang-Gui Liu, Ding-Sheng Wang, and Rushan Han. 2004. “Ab
Initio Phonon Dispersions of Single-Wall Carbon Nanotubes.” *Phys. Rev.
B* 69 (23): 235409. <https://doi.org/10.1103/PhysRevB.69.235409>.

</div>

</div>

[^1]: Thus the harmonic approximation is noticeably in error for very
    asymmetric potentials, which for example can lead to disagreement
    with experimental frequencies on the order of 5% for the librational
    modes of small molecular crystals, or OH bonds. It also neglects
    phonon-phonon interactions and therefore does not contain any line
    broadening.

[^2]: There are several occasions when a consideration of imaginary
    frequencies may be useful. A lattice dynamics calculation can be
    used to test whether a purported equilibrium crystal structure
    really is at equilibrium. The presence of imaginary frequencies
    anywhere in the Brillouin Zone would indicate that the structure is
    unstable to a distortion along the corresponding eigenvector.
    Another occasion would be to investigate a transition state. By
    definition a transition state geometry is at a stationary saddle
    point of the energy hypersurface and must therefore exhibit
    precisely one imaginary frequency. A final example is the case of a
    system exhibiting a “soft-mode” phase transition. In such systems it
    can be useful to explore how the frequency approaches zero as some
    external variable such as pressure or lattice constant is varied, in
    order to precisely locate the phase transition.

[^3]: In general phonon wavevectors ${\mathbf{q}}$ are specified using
    one of the variants of the phonon_kpoint\_\* cell keywords for
    lists, grids or paths. The example given uses a list of length 1 to
    specify a single point.

[^4]: This example will use the DFPT method which is the default if
    phonon_method is not present.

[^5]: Further optional output writing to .castep of the dynamical
    matrixes and force constant matrices is enabled at any iprint level
    by the parameters phonon_write_dynamical and
    phonon_write_force_constants respectively.

[^6]: The Placzek theory applies to the case of insulators only, and no
    general formulation to compute the activity for Raman scattering of
    conducting systems is available

[^7]: If the parameter `raman_method : FINITEDISPLACEMENT` is set the
    calculation uses an older, finite-difference implementation (see
    section <a href="#sec:fd" data-reference-type="ref"
    data-reference="sec:fd">2.3</a>) to compute the mode displacement
    derivatives of the polarizability tensor
    (section <a href="#sec:efield" data-reference-type="ref"
    data-reference="sec:efield">3</a>) using the approach of Porezag and
    Pedersen (Porezag and Pederson 1996). The first stage is to perform
    a full phonon calculation at ${\mathbf{q}}=0$, to determine the mode
    eigenvectors and identify the Raman-active modes. Then CASTEP loops
    over the active modes only computing the Raman tensor, activity and
    depolarization for each. This method is deprecated as obsolete and
    computationally expensive.

[^8]: The offset may be specified as three fractions of the grid
    spacing, so $1/2p \; 1/2q \; 1/2r$ shifts the grid exactly half a
    subdivision in each of the three directions. As of release 24.1 the
    INCLUDE_GAMMA keyword will automate this calculation and in fact is
    the default so that omitting phonon_kpoint_mp_offset entirely will
    do the right thing.

[^9]: The effect is usually small for LDA, intermediate for PBE and
    largest for PW91, as the gradient terms are most sensitive to the
    ASR violation.

[^10]: There is one known case where the reciprocal-space ASR correction
    fails to give the correct behaviour. Two-dimensional planar or
    layered systems with weak inter-layer bonding such as graphene or
    graphite exhibit a quadratic behaviour of an acoustic mode near
    $\Gamma$. The reciprocal space scheme fails to reproduce the
    correct, quadratic behaviour and linearises the dispersion. The
    realspace method generates correct behaviour in these cases.

[^11]: A non-negligible force constant at longer range may also occur in
    a metallic system in the presence of a Kohn anomaly.

[^12]: Because the response of a metallic system to an applied field is
    the generation of a current flow, the dielectric permittivity
    diverges and the Born charges become zero. CASTEP will check that
    the system has a band gap before proceeding with an E-field response
    calculation and abort if there is none.

[^13]: dispersion.pl and dos.pl are also able perform very similar
    analysis and plotting of *electronic* eigenvalues from .castep or
    .bands files and generate band structure and electronic DOS plots.

[^14]: On linux systems xmgrace can usually be found as a contributed
    package and installed using the system package manager. Xmgrace is
    also available for Microsoft Windows systems as part of the “cygwin”
    suite of programs (<http://www.cygwin.com>) along with shells, the
    PERL interpreter and an X-windows server (Xmgrace is an X-windows
    program and requires a running X server to display).

[^15]: A more sophisticated model of infra-red spectra was introduced by
    Balan and Mauri ((Balan et al. 2001; Kendrick and Burnett 2016) and
    subsequent works). Modelling of inelastic neutron spectra is
    discussed by Fair (Fair et al. 2022), Ramirez-Cuesta (Ramirez-Cuesta
    2004; Cheng and Ramirez-Cuesta 2020)

[^16]: equivalent to `opt_strategy : SPEED`

[^17]: equivalent to `opt_strategy : MEMORY`

[^18]: If CASTEP_PAGE_TMPDIR is not set, the code falls back to the
    value of CASTEP_TMPDIR.

[^19]: k-points, plane-waves, bands correspond to indices of the
    wavefunction data, so parallelism in any of these will distribute
    the wavefunctions across MPI processes and reduce the memory/process
    requirement.

[^20]: This activation mechanism will be replaced by a “first-class”
    parameter to set parallelism in future releases of CASTEP. N.B.
    DEVEL_CODE is the only %block allowed in the .param file.

[^21]: See [eCSE
    report](http://www.archer.ac.uk/community/eCSE/eCSE01-017/eCSE01-017_Final_Report_technical.pdf)
    for a description of the implementation and benchmark results
