#------------------------------------------------------------------------------
# Usage:
#       julia xppcm.jl 
#------------------------------------------------------------------------------

# read input parameters
include("input.jl")
include("atomicradii.jl")



# function to write Gaussian input files
function output(geometry)
    # Split the long string geometry into short stings, 
    # and stored in an array called structure
    structure = split(geometry, r"\n\n+(?!$)")
    for i in 1:length(structure)
        content = """
	    %chk=file-$i.chk
            %nproc=$nproc
            %mem=$mem
            #p $method $basis
             
            title

            $charge $multiplicity

            $(structure[i])
         
            """
        open("file-$i".gjf, "w") do f  
            write(f, content)
        end
    end
end

# Physical parameters of the solvent
function solventparameters(solvent)
    if solvent == "cyclohexane"
        𝜀[1] = 2.0165;        # dielectric constant
        𝜌[1] = 0.7781;        # density (gr/cm^3)
        𝑀 = 84.1595;          # molar mass (gr/mol)
        𝑛ₑ = 36;             # number of valence electrons 
        𝑟ₛₒₗ = 2.815;         # molecular radius (Ang.)
    elseif solvent == "benzene"
        𝜀[1] = 2.2706;        # dielectric constant
        𝜌[1] = 0.8756;        # density (gr/cm^3)
        𝑀 = 78.1118;          # molar mass (gr/mol)
        𝑛ₑ = 30;             # number of valence electrons 
        𝑟ₛₒₗ = 2.63;         # molecular radius (Ang.)
    elseif solvent == "argon"
        𝜀[1] = 1.43;        # dielectric constant
        𝜌[1] = 1.3954;        # density (gr/cm^3)
        𝑀 = 39.948;          # molar mass (gr/mol)
        𝑛ₑ = 8;             # number of valence electrons 
        #r_solv = 1.875;         # molecular radius (Ang.)
        𝑟ₛₒₗ = 1.705;         # molecular radius (Ang.)
    else 
        println("Solvent not implemented. Try cyclohexane, benzene, or argon")
    end
end

#------------------------------------------------------------------------------
# Step 1: Computinthe cavity (SES) volume V_{c}(f) as a function of 
#         the cavity scaling factor f: V_{c}(f0)...V_c(fn) 
#-----------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
# Step 2: Computing the permittivity \eps(f) and density RhoS(f) of the 
#          external medium (see Eq.s (9) and (10) of SI in Ref. (2))
#------------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Step 3: Computing the electronic energy G_er(f) as a function of the 
#         cavity scaling factor f. 
#-----------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Step 4: Computing the pressure p=-dG_er/dVc
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Step 5: Computing the mean value of the pressure 
#------------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Step 6 and 7: Cavitation, G_cav (p), and total G_tot(p) energies.
#-----------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Step 8: Computing the reaction energy profile, DG_tot(ns;p), as a function
#         of the pressure p:
#         DG_tot(ns;p)= G_(tot)(ns;p)-[G_tot(ref;p)+G_tot(ref;p)]
#-----------------------------------------------------------------------------


function tidyinput(nproc,mem,method,basis,charge,multiplicity;
    molecularity = "uni",
    solvent = "cyclohexane",
    cavity = "vdw",
    tesserae = 0.075,
    𝜂 = 3,
    scaling = [1.2, 1.15, 1.1, 1.05, 1.0, 0.975, 0.95],
    fitting = "julia")
Dict(
"nproc" => string(nproc),
"mem" => mem,
"method" => method,
"basis" => basis,
"charge" => string(charge),
"multiplicity" => string(multiplicity),
"molecularity" => molecularity,
"solvent" => solvent,
"cavity" => cavity,
"tesserae" => string(tesserae),
"𝜂" => 𝜂,
"scaling" => scaling,
"fitting" => fitting,
)
end