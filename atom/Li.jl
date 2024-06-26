# Example input for XP-PCM calculations on a single atom

solvent = "cyclohexane"      # "cyclohexane", "benzene", or "argon"

# The default dielectric permittivity of the solvent may be set close to 1
# for calculations on charged systems.
#dielectric = 1.0025

radiustype = "bondi"        # "bondi" (default) or "rahm"

# 𝜂
# Empirical Pauli repulsion parameter; recommended values: 3, 6 or 9
# Larger 𝜂 leads to higher calculated pressures
𝜂 = 6

# scalingfactors
# Scaling factors of the vdW atomic radii for constructing the cavity.
# Include the values inside () or [] and separate them by comma.
# Some examples below
#scalingfactors = (1.2,)  # Don't forget the comma
#scalingfactors = 1.2:-0.1:0.9  # from 1.2 to 0.9 with stepsize -0.1
scalingfactors = [1.2, 1.15, 1.1, 1.05, 1.0, 0.975, 0.95]

# The mean area in Å² of the tesserae by which the surface of the cavity
# is partitioned. Suggested value = 0.075.
tesserae = 0.075

# Gaussian 09/16 parameters
nproc = 4     # change to total cpus if ismultithreading = false
mem = "4gb"   # memory per core; change to total memory if ismultithreading = false
keywords = "pbepbe def2svp"    # Gaussian keywords; add more if needed
charge = 0
multiplicity = 2

# Geometry
# Keep the coordinates within the triple """ block.
# Do not include comments or other text in the """ block.
# Leading or trailing spaces on each line are ok.
# Atoms may be specified by element symbols or atomic numbers.
geometry = """
Li    0.0    0.0    0.0

"""
