# example input for XP-PCM calculations
# If not sure about the options below, only modify the geometry and Gaussian 
# realted keywords for your calculation.

# multithreading
# "true" means each Ger job is assigned only 1 cpu core and multiple jobs will 
# be running at the same time. "false" means all cpu cores work on one job at 
# a time. Vc and Gcav jobs are set to be always multithreaded. It is advised to
# use multithreading if the total number of structures are much larger than the 
# total number of cpu cores assigned to the job.
ismultithreading = false

# restart
# Run a new calculation or restart from previously interupted Ger jobs.
# The program will read all existing .log files and figure out which jobs
# are not finished, and then restart the unfinished jobs.
restart = false

# molecularity
# chose "uni" or "bi"-molecular reaction. This option only affect how the 
# relative energies are computed. "uni", the energies are computed relative
# to the first structure. "bi", the energies are computed relative to the sum
# of the energies of the first two structures.
molecularity = "uni"

solvent = "cyclohexane"      # "cyclohexane", "benzene", or "argon"
cavity = "vdw"               # "vdw" or "ses" cavity
sphere = "hard"              # fixed cavity ("hard") or varied cavity ("soft") in Gcav jobs
tesserae = 0.075             # the mean area in Å² of the tesserae by
                             # which the surfaces of the cavity
                             # is partitioned.

# The default dielectric permittivity of the solvent may be set close to 1
# for calculations on charged systems.
#dielectric = 1.0025

radiustype = "bondi"        # "bondi" (default) or "rahm"

# 𝜂
# Empirical Pauli repulsion parameter; recommended values: 3, 6 or 9
# Larger 𝜂 leads to higher calculated pressures
𝜂 = 3

# scalingfactors
# Scaling factors of the vdW atomic radii for constructing the cavity.
# Include the values inside () or [] and separate them by ,
scalingfactors = (1.2, 1.15, 1.1, 1.05, 1.0, 0.975, 0.95)

# Gaussian 09/16 parameters
nproc = 4     # change to total cpus if ismultithreading = false
mem = "2gb"   # memory per core; change to total memory if ismultithreading = false
keywords = "pbepbe def2svp"    # Gaussian keywords; add more if needed
charge = 0
multiplicity = 1

# Keep the coordinates within the triple """ block.
# Separate each geometry by one or more blank lines.
# Do not include comments or other text in the """ block.
# Leading or trailing spaces on each line are ok.
# Atoms may be specified by element symbols or atomic numbers.
geometries = """
H    -0.37 0.0 0.0
H    0.37 0.0 0.0

H    -0.365 0.0 0.0
H    0.365 0.0 0.0

H    -0.36 0.0 0.0
H    0.36 0.0 0.0

"""
