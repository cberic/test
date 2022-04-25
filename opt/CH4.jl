# example input for XP-PCM calculations

# path to the directory of modified Gaussian links
exedir = "/scratch/bochen/Gaussian16/xppcm-links/exe-dir"

solvent = "cyclohexane"      # "cyclohexane", "benzene", or "argon"
cavity = "vdw"               # "vdw" or "ses" (default) cavity

tesserae = 0.075             # the mean area in Å² of the tesserae by
                             # which the surfaces of the cavity
                             # is partitioned. default value = 0.075

# The default dielectric permittivity of the solvent may be set close to 1
# for calculations on charged systems.
#dielectric = 1.0025
#radiustype = "rahm" # The default is "bondi"

# 𝜂
# Empirical Pauli repulsion parameter; recommended values: 3, 6 or 9
# Larger 𝜂 leads to higher calculated pressures
𝜂 = 3

# scalingfactors
# Scaling factors of the vdW atomic radii for constructing the cavity.
# Include the values inside the () and separate them by ,
#scalingfactors = (1.2, 1.15, 1.1, 1.05, 1.0, 0.975, 0.95)
scalingfactors = (1.2, 1.15)

# Gaussian 09/16 parameters
nproc = 2     # change to total cpus if ismultithreading = false
mem = "2gb"   # memory per core; change to total memory if ismultithreading = false
keywords = "pbepbe/cc-pvdz"    # Gaussian keywords; add more if needed
charge = 0
multiplicity = 1

# Keep the coordinates within the triple """ block.
# Separate each geometry by one or more blank lines.
# Do not include comments or other text in the """ block.
# Leading or trailing spaces on each line are ok.
# Atoms may be specified by element symbols or atomic numbers.
geometries = """
 C       0.000000      0.000000      0.000000
 H       0.000000      0.000000      1.105570
 H       1.042341      0.000000     -0.368523
 H      -0.521171     -0.902694     -0.368523
 H      -0.521171      0.902694     -0.368523

"""
