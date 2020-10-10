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

# 𝜂
# Empirical Pauli repulsion parameter; recommended values: 3, 6 or 9
# Larger 𝜂 leads to higher calculated pressures

const 𝜂 = 3

# scalingfactors
# Scaling factors of the vdW atomic radii for constructing the cavity.
# Include the values inside the [] and separate them by ,

const scalingfactors = [1.2, 1.15, 1.1, 1.05, 1.0, 0.975, 0.95]

# fitting
# Program for nonlinear regression fitting for the Murnaghan equation of state.
# Choose from "julia" (recommended), "python", or "mathematica".
# For "julia", install the LsqFit package. For "python", install the lmfit package.
# For "python" or "mathematica", the script will simply call "python" or "math",
# so make sure the call refers to the correct python or mathematica build.

# "Pkg.add("LsqFit")"
const fitting = "julia"

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
const geometries = """

 6 0.799019 0.201672 1.365900
 6 0.799567 -0.089594 -1.377464
 1 1.093133 0.353236 2.405045
 1 1.094169 -0.159706 -2.425121
 6 -1.389967 1.163858 -0.790502
 6 -0.798958 -0.089745 -1.377790
 1 -1.780654 1.935601 -1.443272
 1 -1.093113 -0.159931 -2.425569
 6 1.390104 1.164105 -0.789905
 1 1.780937 1.935925 -1.442497
 6 1.390177 1.303816 0.528256
 1 1.781057 2.195358 1.004390
 6 1.457286 -1.274900 -0.639945
 1 0.998787 -2.217936 -0.934871

 1 2.489541 -1.324618 -0.989990
 6 1.455581 -1.112858 0.894310
 1 2.487001 -1.090702 1.249635
 1 0.994209 -1.972452 1.378623
 6 -1.390625 1.303587 0.527657
 6 -0.799643 0.201553 1.365570
 1 -1.781862 2.195068 1.003610
 1 -1.094214 0.353089 2.404589
 6 -1.455800 -1.113084 0.893720
 1 -2.487389 -1.091060 1.248563
 1 -0.994534 -1.972598 1.378273
 6 -1.456765 -1.275167 -0.640532
 1 -2.488848 -1.325118 -0.991048
 1 -0.997928 -2.218113 -0.935219

 6 0.799239 0.200450 1.366113
 6 0.799579 -0.091184 -1.377450
 1 1.093427 0.352831 2.405110
 1 1.094367 -0.160924 -2.425091
 6 -1.384683 1.164852 -0.790981
 6 -0.798967 -0.091334 -1.377779
 1 -1.767626 1.940189 -1.444028
 1 -1.093306 -0.161147 -2.425541
 6 1.384823 1.165098 -0.790383
 1 1.767913 1.940509 -1.443255
 6 1.385956 1.304346 0.527804
 1 1.769931 2.199204 1.003294
 6 1.461595 -1.273803 -0.639631
 1 1.008024 -2.218920 -0.936296
 1 2.494640 -1.317658 -0.988241
 6 1.456782 -1.113509 0.894724
 1 2.487368 -1.092358 1.252540
 1 0.993785 -1.973416 1.376941
 6 -1.386406 1.304119 0.527204
 6 -0.799865 0.200331 1.365781
 1 -1.770739 2.198916 1.002516
 1 -1.094512 0.352685 2.404652
 6 -1.457003 -1.113734 0.894130
 1 -2.487759 -1.092716 1.251465
 1 -0.994111 -1.973562 1.376588
 6 -1.461074 -1.274070 -0.640223
 1 -2.493950 -1.318156 -0.989306
 1 -1.007163 -2.219098 -0.936650

 6 0.799561 0.199254 1.366271
 6 0.799538 -0.092771 -1.377379
 1 1.093694 0.352453 2.405152
 1 1.094532 -0.162140 -2.424996
 6 -1.379200 1.165932 -0.791509
 6 -0.798925 -0.092922 -1.377708
 1 -1.754286 1.944827 -1.444906
 1 -1.093469 -0.162363 -2.425447
 6 1.379341 1.166177 -0.790912
 1 1.754575 1.945144 -1.444137
 6 1.382102 1.304765 0.527316
 1 1.759561 2.202741 1.002167
 6 1.465966 -1.272630 -0.639269
 1 1.017397 -2.219734 -0.937730
 1 2.499732 -1.310547 -0.986450
 6 1.457996 -1.114226 0.895171
 1 2.487713 -1.094249 1.255574
 1 0.993215 -1.974410 1.375184
 6 -1.382551 1.304538 0.526716
 6 -0.800187 0.199134 1.365938
 1 -1.760369 2.202456 1.001392
 1 -1.094781 0.352306 2.404693
 6 -1.458218 -1.114451 0.894576
 1 -2.488105 -1.094606 1.254497
 1 -0.993540 -1.974556 1.374831
 6 -1.465446 -1.272897 -0.639864
 1 -2.499045 -1.311043 -0.987518
 1 -1.016541 -2.219916 -0.938089



"""
