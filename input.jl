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
const cavity = "vdw"               # "vdw" or "ses" cavity for Gcav jobs
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

6   1.095015   -0.347434    0.000511
6   0.717910    0.959938   -0.000290
6  -0.718022    0.959859    0.000057
6  -1.094973   -0.347557    0.000276
8   0.000067   -1.160967   -0.000523
1   2.050270   -0.849914    0.000841
1   1.373596    1.819548   -0.000523
1  -1.373805    1.819396    0.000095
1  -2.050175   -0.850135    0.000447

6   1.095015   -0.347434    0.000511
6   0.717910    0.959938   -0.000290
6  -0.718022    0.959859    0.000057
6  -1.094973   -0.347557    0.000276
8   0.000067   -1.160967   -0.000523
1   2.050270   -0.849914    0.000841
1   1.373596    1.819548   -0.000523
1  -1.373805    1.819396    0.000095
1  -2.050175   -0.850135    0.000447

6 -1.323136 -1.239946 -1.584710
6 0.011915 -1.416941 -1.777599
6 0.559132 -0.108870 -1.998123
6 -0.488675 0.755656 -1.921291
8 -1.648396 0.082060 -1.671714
1 -2.145120 -1.911022 -1.387243
1 0.545943 -2.356795 -1.761906
1 1.592761 0.147738 -2.183047
1 -0.581482 1.826433 -2.019278
6 -0.211607 1.136702 1.621912
6 -0.471642 -0.141065 2.010073
6 0.802935 -0.793282 2.109373
6 1.732374 0.141615 1.773151
8 1.130687 1.328956 1.474417
1 -0.832729 1.994795 1.416287
1 -1.445402 -0.570402 2.199408
1 0.994425 -1.819478 2.390777
1 2.809366 0.132984 1.702070



"""
