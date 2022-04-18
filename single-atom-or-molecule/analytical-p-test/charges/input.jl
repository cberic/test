# Example input for XP-PCM calculat, n on a single atom .
# Both numerical and analytical pressures will be computed

# Only "cyclohexane", "benzene", and "argon" are currently implemented
solvent = "argon"

# Change the dielectric permittivity of the solvent to be close to 1 for
# calculations on charged systems. Leave the below commented to use the default value.
dielectric = 40
# Choose either the Bondi or Rahm vdW radii and the corresponding set of scaling factors #or ionic
#radiustype = "bondi"
#scalingfactors = [1.2, 1.15, 1.1, 1.05, 1.0, 0.975, 0.95]
radiustype = "ionic"
scalingfactors = [1.00, 0.975, 0.95, 0.925, 0.90, 0.875, 0.85, 0.825, 0.80, 0.775, 0.750, 0.725, 0.700, 0.675, 0.650]

# The mean area in Ang^2 of the tesserae by which the surfaces of the cavity
# is partitioned. Suggested value = 0.075.
tesserae = 0.075

# Empirical Pauli repulsion parameter; recommended values: 3, 6 or 9 #remember that these values are recommended for atoms, for spec. cations you might need higher values.
# Larger 𝜂 leads to higher calculated pressures
𝜂 = 3

# Gaussian 09/16 parameters
# If "gen" is used, provide the custom basis set in a separate file named "gen"
nproc = 4
mem = "4gb"
keywords = "pbe1pbe scf=(Damp,fermi,conver=6) int=finegrid gen sp Charge NoSymm"
#keywords = "pbe1pbe/aug-cc-pvtz"    # Gaussian keywords; add more if needed
charge = -1
multiplicity = 1

# Keep the coordinates (in Angstrom) within the triple """ block.
# Do not include comments or other text in the """ block.
# Blank lines and spaces are ok.
geometries = """

Cl 0.0 0.0 0.0

"""
