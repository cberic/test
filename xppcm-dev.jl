include("input.jl")
include("solvent.jl")
include("atomicradii.jl")
include("geometries.jl")
include("io.jl")
include("rungaussian.jl")
include("arithmetic.jl")

#------------------------------------------------------------------------------
# Step 1: cavity volume 𝑉𝑐(𝑓) Gaussian jobs and solvent property calculations
#------------------------------------------------------------------------------

writegjf("Vc")          # write .gjf files for cavity volume "Vc" calculation 

rungaussian("Vc")       # run Gaussian jobs

get𝑉𝑐()                 # extract cavity volume 𝑉𝑐 from Gaussian output

calculate𝑠()            # calculate linear scaling 𝑠 from 𝑉𝑐 date

average𝑠()              # calculated the average of 𝑠 over all structures

calculate𝜀()            # calculate dielectric permitivity 𝜀

calculate𝑍()            # calculate Pauli repulsion barrier 𝑍

calculate𝑉ₘ()           # calculate the molar volume of solvent 𝑉ₘ

#------------------------------------------------------------------------------
# Step 2: electronic structure Gaussian jobs and pressure calculations 
#------------------------------------------------------------------------------

writegjf("Ger")          # write .gjf files for cavity volume "Ger" calculation 

rungaussian("Ger")       # run Gaussian jobs

get𝐺𝑒𝑟()                 # extract 𝐺𝑒𝑟 from Gaussian output

calculate𝑝()             # calculate pressure 𝑝

#------------------------------------------------------------------------------
# Step 3: cavitation energy Gaussian jobs
#------------------------------------------------------------------------------

writegjf("Gcav")         # write .gjf files for cavitation energy "Gcav" calculation 

rungaussian("Gcav")      # run Gaussian jobs

get𝐸𝑐𝑎𝑣()                 # extract 𝐺𝑐𝑎𝑣 from Gaussian output

calculate𝐺𝑐𝑎𝑣()           # calculate cavitation energy 𝐺𝑐𝑎𝑣

calculate𝐺𝑡𝑜𝑡()            # calculate total energy 𝐺𝑡𝑜𝑡


writeproperties()       # write properties.dat file
