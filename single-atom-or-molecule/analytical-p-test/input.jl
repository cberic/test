#### Example input for XP-PCM calculation on a single atom or ion

#### Only "cyclohexane", "benzene", and "argon" are currently implemented
solvent = "cyclohexane"

#### User defined dielectric permittivity of the solvent or the crystal
#### Comment it out to use the default value
dielectric = 1.64^2 #square of the refractive index

#### Choose "Bondi" or "Rahm" vdW radii for atoms, or "Rahm_ionic" for ions
#### with the suggested scaling factors for the corresponding radius type
#radiustype = "bondi"
#scalingfactors = [1.2, 1.15, 1.1, 1.05, 1.0, 0.975, 0.95]
#radiustype = "rahm"
#scalingfactors = [1.3, 1.25, 1.2, 1.15, 1.1, 1.05, 1.0, 0.975, 0.95, 0.925, 0.90]
radiustype = "rahm_ionic"
scalingfactors = [1.0, 0.975, 0.95, 0.925, 0.9, 0.8, 0.75, 0.725]

#### Additional parameters for radiustype = "rahm_ionic"
#lattice = "NaCl"
#𝛼ᵣ = 1.747565  # Madelung constant, https://en.wikipedia.org/wiki/Madelung_constant
#𝑟₀ = 2.820 * 1.88973  # https://en.wikipedia.org/wiki/Born–Mayer_equation
#lattice = "CsCl"
#𝛼ᵣ = 1.762675  # Madelung constant, https://en.wikipedia.org/wiki/Madelung_constant
#𝑟₀ = 3.571 * 1.88973 # closest distance between two ions of opposite charge in bohr
lattice = "noLattice"
𝛼ᵣ = 0.0
𝑟₀ = 3.571 * 1.88973 # closest distance between two ions of opposite charge in bohr

#### The mean area in Å² of the tesserae by which the surface of the cavity
#### is partitioned. Suggested value = 0.075.
tesserae = 0.075

#### Empirical Pauli repulsion parameter; recommended values: 3, 6 or 9
#### Larger 𝜂 leads to higher calculated pressures
𝜂 = 5

#### Gaussian 09/16 parameters
#### If "gen" is used, provide the custom basis set in a separate file named "gen"
nproc = 4
mem = "4gb"
keywords = "pbe1pbe gen scf=(Damp,fermi,conver=6) int=finegrid"
#keywords = "pbe1pbe/6-311++(2d,2p)  scf=(Damp,fermi,conver=6) int=finegrid"    # Gaussian keywords; add more if needed
charge = -1
multiplicity = 1

#### Keep the coordinates (in Å) within the triple """ block.
#### Do not include comments or other text in the """ block.
#### Blank lines and spaces are ok.
geometries = """

Cl  0.0  0.0  0.0

"""
