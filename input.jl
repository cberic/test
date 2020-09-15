# example input for XP-PCM calculations

const molecularity = "uni"         # uni (default) or bi -molecular reaction
const solvent = "cyclohexane"      # cyclohexane (default), benzene, or argon

const cavity = "vdw"               # vdw (default) or ses; for cavitation 
                                   # energy calculation only

const tesserae = 0.075             # the mean area in Ang^2 of the tesserae by 
                                   # which the surfaces of the cavity     		                
                                   # is partitioned. default value = 0.075 

const 𝜂 = 3                        # empirical Pauli repulsion parameter
                                   # 3 (default), 6 or 9

const scalingfactors = [1.2, 1.15, 1.1, 1.05, 1.0, 0.975, 0.95]   # (default)

const fitting = "python"            # Murnaghan equation of state fitting method
                                   # julia (default), python, or mathematica
# Gaussian 16 parameters
const nproc = 8
const mem = "1gb"
const method = "wb97xd"
const basis = "def2tzvp"
const charge = 0
const multiplicity = 1

# Keep the coordinates within the triple """ block.
# Separate each geometry by one or more blank lines.
# Do not include comments or other text in the """ block.
# Leading or trainling spaces are ok.
const geometries = """

    C    1.10712900   -1.48465400    0.00000000  
    C   -0.00001100   -0.72877000    0.00000000  
    C    0.00000000    0.72875800    0.00000000  
    C   -1.10712300    1.48466300    0.00000000  
    H    2.09998900   -1.03986800    0.00000000  
    H    1.05995700   -2.56940600    0.00000000  
    H   -0.97811200   -1.21104800    0.00000000  
    H    0.97811000    1.21101800    0.00000000  
    H   -2.09999700    1.03991100    0.00000000  
    H   -1.05991600    2.56941400    0.00000000  
     

C    1.10712900   -1.48465400    0.00000000
C   -0.00001100   -0.72877000    0.00000000
C    0.00000000    0.72875800    0.00000000
C   -1.10712300    1.48466300    0.00000000
H    2.09998900   -1.03986800    0.00000000
H    1.05995700   -2.56940600    0.00000000
H   -0.97811200   -1.21104800    0.00000000
H    0.97811000    1.21101800    0.00000000
H   -2.09999700    1.03991100    0.00000000
H   -1.05991600    2.56941400    2.00000000


     

  C    1.10712900   -1.48465400    0.00000000
C   -0.00001100   -0.72877000    0.00000000
C    0.00000000    0.72875800    0.00000000
C   -1.10712300    1.48466300    0.00000000
   H    2.09998900   -1.03986800    0.00000000
H    1.05995700   -2.56940600    0.00000000
H   -0.97811200   -1.21104800    0.00000000
H    0.97811000    1.21101800    0.00000000
H   -2.09999700    1.03991100    0.00000000
    H   -1.05991600    2.56941400    3.00000000

  

  
"""
