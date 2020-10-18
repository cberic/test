# Example input for XP-PCM calculation on a single atom or molecule

# Only "cyclohexane", "benzene", and "argon" are currently implemented
solvent = "cyclohexane"

# The default dielectric permittivity of the solvent may be set close to 1
# for the calculations on charged systems. Comment it out to use the default value.
#dielectric = 1.1

# Choose either the Bondi or Rahm vdW radii and the corresponding set of 
# scaling factors
#radiustype = "bondi"
#scalingfactors = [1.2, 1.15, 1.1, 1.05, 1.0, 0.975, 0.95]
radiustype = "rahm"
scalingfactors = [1.3, 1.25, 1.2, 1.15, 1.1, 1.05, 1.0, 0.975, 0.95]

# The mean area in Ang^2 of the tesserae by which the surfaces of the cavity
# is partitioned. Suggested value = 0.075.
tesserae = 0.075

# Empirical Pauli repulsion parameter; recommended values: 3, 6 or 9
# Larger 𝜂 leads to higher calculated pressures
𝜂 = 6

# Gaussian 09/16 parameters
nproc = 1
mem = "1gb"
keywords = "b3lyp 6-31g* int=finegrid"    # Gaussian keywords; add more if needed
charge = 0
multiplicity = 1

# Keep the coordinates (in Angstrom) within the triple """ block.
# Do not include comments or other text in the """ block.
# Blank lines and spaces are ok.
geometries = """

He 0.0 0.0 0.0

"""
