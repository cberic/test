# example input for XP-PCM calculations

# multithreading
# Turn multithreading "on" (highly recommended) or "off" for Ger jobs. 
# On means each Ger job is assigned only 1 cpu core and multiple jobs will 
# be running at the same time. Off means all cpu cores work on one job at 
# a time. Tests show that multithreading is much faster. Vc and Gcav jobs 
# are set to be always multithreaded.
const multithreading = "on"

# restart
# "yes" or "no". Notice the lower case and the double quotation.
# No means to start a new job from begining.
# Yes means to restart from previously interupted Ger jobs. The program will 
# read all existing Ger.log files and figure out which jobs are finished and 
# which are not, and then restart the Ger jobs for the unfinished jobs.
const restart = "no"

const molecularity = "bi"         # "uni" or "bi"-molecular reaction
const solvent = "cyclohexane"      # "cyclohexane", "benzene", or "argon"
const cavity = "ses"               # "vdw" or "ses" cavity for Gcav jobs
const sphere = "hard"              # soft sphere model not yet implemented
const tesserae = 0.075             # the mean area in Ang^2 of the tesserae by 
                                   # which the surfaces of the cavity
                                   # is partitioned. default value = 0.075 

# The default dielectric permittivity of the solvent may be set close to 1
# for calculations on charged systems.
const dielectric = 1.0025

# 𝜂
# Empirical Pauli repulsion parameter; recommended values: 3, 6 or 9
# Larger 𝜂 leads to higher calculated pressures
const 𝜂 = 3

# scalingfactors
# Scaling factors of the vdW atomic radii for constructing the cavity.
# Include the values inside the [] and separate them by ,
const scalingfactors = [1.2, 1.15, 1.1, 1.05, 1.0, 0.975, 0.95]

# Gaussian 09/16 parameters
const nproc = 1     # change to total cpus if multithreading is off
const mem = "1gb"   # memory per core; change to total memory if multithreading is off
const keywords = "b3lyp 6-31g* int=finegrid"    # Gaussian keywords; add more if needed
const charge = 0
const multiplicity = 1

# Keep the coordinates within the triple """ block.
# Separate each geometry by one or more blank lines.
# Do not include comments or other text in the """ block.
# Leading or trailing spaces on each line are ok.
# Atoms may be specified by element symbols or atomic numbers.
const geometries = """
H -0.37 0.0 0.0
H 0.37 0.0 0.0

H -0.37 0.0 0.0
H 0.37 0.0 0.0

H 1.0 0.0 0.0
H 1.74 0.0 0.0
H -1.0 0.0 0.0
H -1.74 0.0 0.0

H 0.5 0.0 0.0
H 1.24 0.0 0.0
H -0.5 0.0 0.0
H -1.24 0.0 0.0

"""
