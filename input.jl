# example input for XP-PCM calculations

# multithreading
# Turn multithreading "on" (highly recommended) or "off" for Ger jobs. 
# On means each Ger job is assigned only 1 cpu core and multiple jobs will 
# be running at the same time. Off means all cpu cores work on one job at 
# a time. Tests show that multithreading is much faster. Vc and Gcav jobs 
# are set to be always multithreaded.
ismultithreading = false

# restart
# Run a new calculation or restart from previously interupted Ger jobs.
# The program will read all existing Ger.log files and figure out which jobs
# are not finished, and then restart the Ger jobs for the unfinished jobs.
restart = true

molecularity = "uni"          # "uni" or "bi"-molecular reaction
solvent = "cyclohexane"      # "cyclohexane", "benzene", or "argon"
cavity = "vdw"               # "vdw" or "ses" (default) cavity
sphere = "hard"              # fixed cavity ("hard", which is the default)
                             # or varied cavity ("soft") in Gcav job
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
scalingfactors = (1.2, 1.15, 1.1, 1.05, 1.0, 0.975, 0.95)

# Gaussian 09/16 parameters
nproc = 4     # change to total cpus if ismultithreading = false
mem = "4gb"   # memory per core; change to total memory if ismultithreading = false
keywords = "uwb97xd def2svp guess=mix"    # Gaussian keywords; add more if needed
charge = 0
multiplicity = 1

# Keep the coordinates within the triple """ block.
# Separate each geometry by one or more blank lines.
# Do not include comments or other text in the """ block.
# Leading or trailing spaces on each line are ok.
# Atoms may be specified by element symbols or atomic numbers.
geometries = """
6                     1.1       0.71888   0. 
6                     1.1      -0.71888   0. 
6                    -1.1       0.71888   0. 
6                    -1.1      -0.71888   0. 
1                     1.30832   1.26099   0.92589 
1                     1.30832   1.26099  -0.92589 
1                     1.30832  -1.26099   0.92589 
1                     1.30832  -1.26099  -0.92589 
1                    -1.30832   1.26099   0.92589 
1                    -1.30832   1.26099  -0.92589 
1                    -1.30832  -1.26099   0.92589 
1                    -1.30832  -1.26099  -0.92589 

"""
