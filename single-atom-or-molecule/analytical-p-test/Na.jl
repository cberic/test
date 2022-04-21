#### Example input for XP-PCM calculation on a single atom or ion

#### "xenon" for CsCl compression chamber; "krypton" for NaCl compression chamber
solvent = "xenon"

#### User defined dielectric permittivity of the solvent or the crystal
#dielectric = 7.2 #1.0025

#### Choose "Bondi" or "Rahm" vdW radii for atoms, or "Rahm_ionic" for ions
radiustype = "rahm_ionic"
scalingfactors = [1.0, 0.95, 0.925, 0.9, 0.85, 0.8, 0.75, 0.725, 0.7]

#### Additional parameters for radiustype = "rahm_ionic"
#lattice = "NaCl"
#𝛼ᵣ = 1.747565  # Madelung constant, https://en.wikipedia.org/wiki/Madelung_constant
#𝑟₀ = 2.820 * 1.88973  # https://en.wikipedia.org/wiki/Born–Mayer_equation
lattice = "noLattice" #"CsCl"
𝛼ᵣ = 1.762675  # Madelung constant, https://en.wikipedia.org/wiki/Madelung_constant
𝑟₀ = 3.571 * 1.88973 # closest distance between two ions of opposite charge in bohr

#### The mean area in Å² of the tesserae by which the surface of the cavity
#### is partitioned. Suggested value = 0.075.
tesserae = 0.075

#### Empirical Pauli repulsion parameter; recommended values: 3, 6 or 9
#### Larger 𝜂 leads to higher calculated pressures
𝜂 = 5

#### Gaussian 09/16 parameters
#### If "gen" is used, provide the custom basis set in a separate file and pass
#### the file name to the variable gen_filename
nproc = 4
mem = "4gb"
keywords = "pbe1pbe gen scf=(Damp,fermi,conver=6) int=finegrid"
charge = 1
multiplicity = 1
gen_filename = "ano-rcc-Na"

#### Keep the coordinates (in Å) within the triple """ block.
#### Do not include comments or other text in the """ block.
#### Blank lines and spaces are ok.
geometries = """

Na  0.0  0.0  0.0

"""
