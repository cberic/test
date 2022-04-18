#### Example input for XP-PCM calculation on a single atom or ion

#### "xenon" for CsCl compression chamber; "krypton" for NaCl compression chamber
solvent = "cyclohexane"

#### User defined dielectric permittivity of the solvent or the crystal
dielectric = 2.0 #7.2 #1.0025

#### Choose "Bondi" or "Rahm" vdW radii for atoms, or "ionic" for ions
radiustype = "ionic"
scalingfactors = [1.0, 0.95, 0.925]
#scalingfactors = [1.0, 0.95, 0.925, 0.9, 0.85, 0.8, 0.75, 0.725, 0.7]

#### Additional parameters for radiustype = "ionic"
#lattice = "NaCl"
#𝛼ᵣ = 1.747565  # Madelung constant, https://en.wikipedia.org/wiki/Madelung_constant
#𝑟₀ = 2.820 * 1.88973  # https://en.wikipedia.org/wiki/Born–Mayer_equation
#lattice = "noLattice" #"CsCl"
#𝛼ᵣ = 1.762675  # Madelung constant, https://en.wikipedia.org/wiki/Madelung_constant
#𝑟₀ = 3.571 * 1.88973 # closest distance between two ions of opposite charge in bohr
model = "pointcharges"  # evoking compensating, spherically-distributed point charges
                        # do not write "charge" in the keyword variable below
#model = "Born"  # Born model for computing the lattice Madelung energy
impose_equilibrium = true  # imposing p=0 at the 1st scaling factor

#### The mean area in Å² of the tesserae by which the surface of the cavity
#### is partitioned. Suggested value = 0.075.
tesserae = 0.075

#### Empirical Pauli repulsion parameter; recommended values: 3, 6 or 9
#### Larger 𝜂 leads to higher calculated pressures
𝜂 = 9

#### Gaussian 09/16 parameters
#### If "gen" is used, provide the custom basis set in a separate file and pass
#### the file name to the variable gen_filename
nproc = 4
mem = "4gb"
keywords = "pbe1pbe gen scf=(Damp,fermi,conver=6) int=finegrid"
charge = -1
multiplicity = 1
gen_filename = "svp"

#### Keep the coordinates (in Å) within the triple """ block.
#### Do not include comments or other text in the """ block.
#### Blank lines and spaces are ok.
geometries = """

Cl  0.0  0.0  0.0

"""
