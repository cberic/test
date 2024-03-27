# example input for XP-PCM optimizations

# path to the directory of modified Gaussian links
exedir = "/scratch/bochen/Gaussian16/xppcm-links/exe-dir"

solvent = "cyclohexane"      # "cyclohexane", "benzene", or "argon"
cavity = "custom"            # "vdw" or "ses" (default) cavity, or "custom"

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
nproc = 4     # change to total cpus if ismultithreading = false
mem = "4gb"   # memory per core; change to total memory if ismultithreading = false
keywords = "pbepbe/cc-pvdz"    # Gaussian keywords; add more if needed
charge = 0
multiplicity = 1

# Geometry
# Can be provided in Cartesian or Z-matrix format. 
# When using Z-matrix, `atomlist` needs to be provided.

#= cartesian = """
          6                   0.000000    0.000000    0.000000
          1                   0.000000    0.000000    1.105570
          1                   1.042341    0.000000   -0.368523
          1                  -0.521171   -0.902694   -0.368523
          1                  -0.521171    0.902694   -0.368523
""" =#

atomlist = "C H H  H H  "
zmatrix = """
C
H  1 r1
H  1 r1  2 a
H  1 r1  2 a  3  d1
H  1 r1  2 a  3  d2
Variables:
r1=1.1055698
Constant:
a=109.471220
d1=120.0
d2=-120.0
"""

# Custom cavity
# When `cavity = custom`, list the sphere locations (on which atoms) and radii.
spherespec = [2  2.04 
             4  2.06
             3  2.08]
