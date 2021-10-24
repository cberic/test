#### Example input for XP-PCM calculation on a single atom or ion

#### Only "cyclohexane", "benzene", and "argon" are currently implemented
solvent = "cyclohexane"

#### User defined dielectric permittivity of the solvent
#### Comment it out to use the default value
dielectric = 1.0025

#### Choose "Bondi" or "Rahm" vdW radii for atoms, or "Rahm_ionic" for ions
#### with the suggested scaling factors for the corresponding radius type
#radiustype = "bondi"
#scalingfactors = [1.2, 1.15, 1.1, 1.05, 1.0, 0.975, 0.95]
#radiustype = "rahm"
#scalingfactors = [1.3, 1.25, 1.2, 1.15, 1.1, 1.05, 1.0, 0.975, 0.95, 0.925, 0.90]
radiustype = "rahm_ionic"
scalingfactors = [1.3, 1.25, 1.2, 1.15, 1.1, 1.05, 1.0, 0.95]

#### Additional parameters for radiustype = "rahm_ionic"
#### NaCl, simple cubic lattice
#𝛼ᵣ = 1.747565  # Madelung constant, https://en.wikipedia.org/wiki/Madelung_constant
#𝑟₀ = 2.820 * 1.88973  # https://en.wikipedia.org/wiki/Born–Mayer_equation
#### CsCl, body-centered cubic lattice
𝛼ᵣ = 1.762675  # Madelung constant, https://en.wikipedia.org/wiki/Madelung_constant
𝑟₀ = 3.571 * 1.88973 # closest distance between two ions of opposite charge in bohr


#### The mean area in Å² of the tesserae by which the surface of the cavity
#### is partitioned. Suggested value = 0.075.
tesserae = 0.075

#### Empirical Pauli repulsion parameter; recommended values: 3, 6 or 9
#### Larger 𝜂 leads to higher calculated pressures
𝜂 = 2

#### Gaussian 09/16 parameters
#### If "gen" is used, provide the custom basis set in a separate file named "gen"
nproc = 4
mem = "4gb"
keywords = "pbe1pbe 6-31+g* scf=(Damp,fermi,conver=6) int=finegrid"
#keywords = "pbe1pbe/6-311++(2d,2p)  scf=(Damp,fermi,conver=6) int=finegrid"    # Gaussian keywords; add more if needed
charge = -1
multiplicity = 1

#### Keep the coordinates (in Å) within the triple """ block.
#### Do not include comments or other text in the """ block.
#### Blank lines and spaces are ok.
geometries = """

Cl  0.0  0.0  0.0

"""
