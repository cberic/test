using Statistics
using Printf
using DelimitedFiles
using LsqFit

include("input.jl")


#------------------------------------------------------------------------------
# solvent.jl
#------------------------------------------------------------------------------
# returns a 5-element tuple of (dielectric constant 𝜀, valence electron density 𝜌,
# molar mass 𝑀, number of valence electrons, molecular radius (Ang.) )
function solventparameters(s = solvent)
    if s == "cyclohexane"
        return (2.0165, 0.7781, 84.1595, 36, 2.815)
    elseif s == "benzene"
        return (2.2706, 0.8756, 78.1118, 30, 2.63)
    elseif s == "argon"
        return (1.43, 1.3954, 39.948, 8, 1.705)
    else 
        println("Solvent not implemented. Try cyclohexane, benzene, or argon")
    end
end

#------------------------------------------------------------------------------
# atomicradii.jl
#------------------------------------------------------------------------------
function atomicradii()
    return Dict(
    "H" => 1.20,    "1" => 1.20,
    "B" => 1.92,    "5" => 1.92,  
    "C" => 1.70,    "6" => 1.70,
    "N" => 1.55,    "7" => 1.55,
    "O" => 1.52,    "8" => 1.52,
    "F" => 1.47,    "9" => 1.47)
end

#------------------------------------------------------------------------------
# geometries.jl
#------------------------------------------------------------------------------
function tidygeometries(geom = geometries)
    # remove leading and trailing spaces/blank lines 
    # and add a space to the beginning
    # "*" is string concatenation operator
    a = " " * strip(geom)

    # change block separators to "NEXT"
    # r"" represents a regular expression
    b = replace(a, r"\s*\n+\s*\n+\s*" => "NEXT")
    
    # remove leading and trailing spaces of each coordinate line 
    # and add a space to the begininig of coordinate line
    c = replace(b, r"\s*\n\s*" => "\n ")
    
    # change the block separators back to "\n\n" and add a space
    return replace(c, "NEXT" => "\n\n ")
end


# 1D array of geometry blocks
function geometryblocks(geom = geometries)
    return split(tidygeometries(geom), "\n\n", keepempty=false)
end


function numberofstructures(geom = geometries)
    return length(geometryblocks(geom))
end


# 1D array whose elements are 1D arrays of lines in each geometry block
function structurelines(geom = geometries)
    blocks = geometryblocks(geom)
    nos = numberofstructures(geom)
    # lines is an array of length nos of 1D arrays
    # lines[i] is a 1D array of length numberofatoms
    lines = Array{Array{String}}(undef, nos)
    for i in 1:nos
        lines[i] = split(blocks[i], "\n", keepempty=false)
    end
    return lines
end


function numberofatoms(geom = geometries)
    return length(structurelines(geom)[1])
end


# nos * noa 2D array of atom lists
function atomlist(geom = geometries)
    nos = numberofstructures(geom)
    lines = structurelines(geom)
    noa = numberofatoms(geom)
    atoms = Array{String}(undef, nos,noa)
    for i in 1:nos
        for j in 1:noa
            # get the 1st field of each coordinate line
            atoms[i,j] = split(lines[i][j], keepempty=false, limit=2)[1]
        end
    end
    return atoms
end


# nos * noa 2D array of xyz coordinates
function coordinatelines(geom = geometries)
    nos = numberofstructures(geom)
    lines = structurelines(geom)
    noa = numberofatoms(geom)
    coordlines = Array{String}(undef, nos,noa)
    for i in 1:nos
        for j in 1:noa
            # get the 2nd-last fields of each coordinate line & add a leading space
            coordlines[i,j] = " " * split(lines[i][j], keepempty=false, limit=2)[2]
        end
    end
    return coordlines
end


# nos * noa * 3 3D array of each coordinates; not correct
#function atomcoordinates(geom = geometries)
    # split into a long 1D array of length nos * nos * 4
#    array = split(geom)
#    nos = numberofstructures(geom)
#    noa = numberofatoms(geom)
#    coordinates = Array{String}(undef, nos,noa,3)
#    for i in 1:nos
#        for j in 1:noa
#            for k in 1:3
#                coordinates[i,j,k] = array[(i-1)*40 + (j-1)*4 + k+1]
#            end
#        end
#    end
#    return coordinates
#end


#------------------------------------------------------------------------------
# io.jl
#------------------------------------------------------------------------------
function writegjf(jobtype)
    if jobtype == "Vc"        # Initial cavitation volume jobs
        gjfvc()
    elseif jobtype == "Ger"   # Electronic energy SCRF jobs
        gjfger()
    elseif jobtype == "Gcav"  # Cavitation energy jobs
        gjfgcav()
    end
end


# write gjf files for Vc jobs
function gjfvc(geom = geometries, 𝑓 = scalingfactors)
    a = length(𝑓)
    geomblocks = geometryblocks(geom)
    nos = numberofstructures(geom)
    noa = numberofatoms(geom)
    atoms = atomlist(geom)
    coordlines = coordinatelines(geom)
    sp = solventparameters()
    𝑟ₐ = atomicradii()
    
    run(`mkdir -p tmp`)
    
    Threads.@threads for i in 1:nos  # use multithreading
        for j in 1:a
            # writing mode for the first and appending mode for other 𝑓
            open("tmp/structure-$i-Vc.gjf", "$(j == 1 ? "w" : "a")") do file
                write(file, """
                    %kjob l301
                    %nproc=1
                    %mem=1000mb
                    #p $keywords
                    # scrf=(iefpcm,solvent=$solvent,read) nosym guess=only pop=none
                    
                    title
                    
                    $charge $multiplicity
                    $(geomblocks[i])
                    
                    qrep pcmdoc geomview nodis nocav g03defaults tsare=$tesserae
                    nsfe=$noa
                    nvesolv=$(sp[4]) solvmw=$(sp[3]) rsolv=$(sp[5])
                    eps=$(sp[1]) rhos=$(sp[2])
                    
                    """)
                
                for k in 1:noa
                    write(file, " $(coordlines[i,k])    $(𝑟ₐ[atoms[i,k]])    $(𝑓[j])\n")
                end
                
                write(file, "\n")
                
                if j != a  # do not write --link1-- for the last scaling factor
                    write(file, "--link1--\n")
                end
            end
        end
    end
end


# write gjf files for Ger calculations
function gjfger(geom = geometries, 𝑓 = scalingfactors)
    a = length(𝑓)
    geomblocks = geometryblocks(geom)
    nos = numberofstructures(geom)
    noa = numberofatoms(geom)
    atoms = atomlist(geom)
    coordlines = coordinatelines(geom)
    sp = solventparameters()
    𝑟ₐ = atomicradii()
    𝜀 = calculate𝜀()    # data of 𝜀 and 𝑍 needed for the gjf files
    𝑍 = calculate𝑍()
    
    Threads.@threads for i in 1:nos  # use multithreading
        for j in 1:a
            # writing mode for the first and appending mode for other 𝑓
            open("tmp/structure-$i-Ger.gjf", "$(j == 1 ? "w" : "a")") do file
                write(file, """
                    $(j == 1 ? "" : "%kjob l502\n")%chk=structure-$i-Ger.chk
                    %nproc=$nproc
                    %mem=$mem
                    #p $keywords $(j == 1 ? "" : "guess=read")
                    # scrf=(iefpcm,solvent=$solvent,read) nosym 6d 10f
                    
                    title
                    
                    $charge $multiplicity
                    $(geomblocks[i])
                    
                    qrep pcmdoc geomview nodis nocav g03defaults tsare=$tesserae
                    nsfe=$noa
                    nvesolv=$(sp[4]) solvmw=$(sp[3]) rsolv=$(sp[5])
                    eps=$(𝜀[j]) rhos=$(𝑍[j])
                    
                    """)
                
                for k in 1:noa
                    write(file, " $(coordlines[i,k])    $(𝑟ₐ[atoms[i,k]])    $(𝑓[j])\n")
                end
                
                write(file, "\n")
                # do not write "--link1--" for the last scaling factor
                if j != a
                    write(file, "--link1--\n")
                end
            end
        end
    end
end


# write gjf files for Gcav calculations
function gjfgcav(cav = cavity, geom = geometries, 𝑓 = scalingfactors)
    a = length(𝑓)
    geomblocks = geometryblocks(geom)
    nos = numberofstructures(geom)
    noa = numberofatoms(geom)
    atoms = atomlist(geom)
    sp = solventparameters()
    𝑟ₐ = atomicradii()
    𝑉ₘ = calculate𝑉ₘ()    # molar volume 𝑉ₘ of the solvent
    
    Threads.@threads for i in 1:nos  # use multithreading
        for j in 1:a
            # writing mode for the first and appending mode for other 𝑓
            open("tmp/structure-$i-Gcav.gjf", "$(j == 1 ? "w" : "a")") do file
                # also testing for vdw cavity to add the noaddsph keyword
                write(file, """
                    %kjob l301
                    %nproc=1
                    %mem=1000mb
                    #p $keywords
                    # scrf=(iefpcm,solvent=$solvent,read) nosym guess=only pop=none
                    
                    title
                    
                    $charge $multiplicity
                    $(geomblocks[i])
                    
                    norep pcmdoc geomview nodis cav g03defaults tsare=$tesserae
                    nsfe=$noa $(cav == "vdw" ? "noaddsph" : "") 
                    Vmol=$(𝑉ₘ[j]) rsolv=$(sp[5])
                    
                    """)
                
                for k in 1:noa
                    write(file, " $k    $(𝑟ₐ[atoms[i,k]] * 𝑓[1])    1.0\n")
                end
                
                write(file, "\n")
                # do not write "--link1--" after the last structure
                if j != a
                    write(file, "--link1--\n")
                end
            end
        end
    end
end


# extract volume 𝑉𝑐 data from gaussian output files
function get𝑉𝑐(geom = geometries, 𝑓 = scalingfactors)
    nos = numberofstructures(geom)
    a = length(𝑓)
    𝑉𝑐 = Array{Float64}(undef, nos,a)    # nos * a 2D array
    Threads.@threads for i in 1:nos
        j = 1    # j ranges from 1:length(𝑓)
        open("tmp/structure-$i-Vc.log") do file
            for line in eachline(file)
                if occursin("GePol: Cavity volume", line)
                    𝑉𝑐[i,j] = parse(Float64, split(line)[5])
                    j += 1    # j ranges from 1:length(𝑓)
                end 
            end
        end
    end
    return 𝑉𝑐
end


# extract electronic energy 𝐺𝑒𝑟 data from gaussian output files
function get𝐺𝑒𝑟(geom = geometries, 𝑓 = scalingfactors)
    nos = numberofstructures(geom)
    a = length(𝑓)
    𝐺𝑒𝑟 = Array{Float64}(undef, nos,a)    # 2D array with dimensions nos * a
    Threads.@threads for i in 1:nos
        j = 1    # j ranges from 1:length(𝑓)
        open("tmp/structure-$i-Ger.log") do file
            for line in eachline(file)
                if occursin("SCF Done", line)
                    𝐺𝑒𝑟[i,j] = parse(Float64, split(line)[5])
                    j += 1    # j ranges from 1:length(𝑓)
                end
            end
        end
    end
    return 𝐺𝑒𝑟
end


# extract cavity volume 𝑉𝑐𝑎𝑣 & non-electrostatic cavitation energy 𝐸𝑐𝑎𝑣
function get𝑉𝑐𝑎𝑣𝐸𝑐𝑎𝑣(geom = geometries, 𝑓 = scalingfactors)
    nos = numberofstructures(geom)
    a = length(𝑓)
    𝑉𝑐𝑎𝑣 = Array{Float64}(undef, nos,a)    # nos * a 2D array
    𝐸𝑐𝑎𝑣 = Array{Float64}(undef, nos,a)    # nos * a 2D array
    Threads.@threads for i in 1:nos
        j = 1    # j ranges from 1:length(𝑓)
        open("tmp/structure-$i-Gcav.log") do file
            for line in eachline(file)
                if occursin("GePol: Cavity volume", line)
                    𝑉𝑐𝑎𝑣[i,j] = parse(Float64, split(line)[5])
                elseif occursin("PCM non-electrostatic energy", line)
                    𝐸𝑐𝑎𝑣[i,j] = parse(Float64, split(line)[5])
                    j += 1    # j ranges from 1:length(𝑓)
                end
            end
        end
    end
    return 𝑉𝑐𝑎𝑣, 𝐸𝑐𝑎𝑣
end


# use the Printf package to write the properties.dat file
function writeproperties(𝑉𝑐 = 𝑉𝑐, 𝐺𝑒𝑟 = 𝐺𝑒𝑟, 𝐸𝑐𝑎𝑣 = 𝐸𝑐𝑎𝑣, geom = geometries, 𝑓 = scalingfactors)
    nos = numberofstructures(geom)
    a = length(𝑓)
    𝑠 = calculate𝑠()
    𝑠̄ = average𝑠()
    𝜀 = calculate𝜀()
    𝑍 = calculate𝑍()
    𝑉ₘ = calculate𝑉ₘ()
    𝑝 = calculate𝑝()
    𝑝̄ = average𝑝()
    𝐺𝑐𝑎𝑣 = calculate𝐺𝑐𝑎𝑣()
    𝐺𝑡𝑜𝑡 = calculate𝐺𝑡𝑜𝑡()
    Δ𝐺𝑡𝑜𝑡 = calculateΔ𝐺𝑡𝑜𝑡()
    open("properties.dat", "w") do file
        for i in 1:nos
            write(file, "structure $i\n")
            write(file, "#    𝑓       𝑉𝑐(𝑓)      𝑠(𝑓)         𝑠̄(𝑓,𝑛ₛ)      𝜀(𝑠̄)        𝑍(𝑠̄)        𝑉ₘ(𝑠̄)      𝐺𝑒𝑟(𝑓)          𝑝(𝑓)     𝑝̄(𝑠̄)      𝐸𝑐𝑎𝑣(𝑓)      𝐺𝑐𝑎𝑣(𝑓)      𝐺𝑡𝑜𝑡(𝑓)      Δ𝐺𝑡𝑜𝑡\n")
            for j in 1:a
                @printf(file, "%d    %.2f    %7.3f    %.6f    %.6f    %.6f    %.6f    %7.3f    %.8f    %.3f    %.3f    %.8f    %.8f    %.8f    %.2f\n", 
                               j,    𝑓[j],   𝑉𝑐[i,j], 𝑠[i,j],  𝑠̄[j],   𝜀[j],   𝑍[j],  𝑉ₘ[j], 𝐺𝑒𝑟[i,j], 𝑝[i,j], 𝑝̄[j], 𝐸𝑐𝑎𝑣[i,j], 𝐺𝑐𝑎𝑣[i,j], 𝐺𝑡𝑜𝑡[i,j], Δ𝐺𝑡𝑜𝑡[i,j])
            end
            write(file, "\n")
        end
    end
end


#------------------------------------------------------------------------------
# rungaussian.jl
#------------------------------------------------------------------------------
# check if g16 or g09 is installed and loaded
function gaussianversion()
    script = raw"""
        #!/usr/bin/bash
        which g16 > /dev/null 2>&1
        a=$?
        which g09 > /dev/null 2>&1
        b=$?
        if [[ $a == 0 ]]; then
            echo g16
        elif [[ $b == 0 ]]; then
            echo g09
        else
            echo "missing"
        fi"""
    open("checkgaussian.sh", "w") do file
        write(file, script)
    end
    version = strip(read(`bash checkgaussian.sh`, String))
    rm("checkgaussian.sh")
    if version === Missing
        println("g16 and g09 not found")
    else
        return version
    end    
end


function rungaussian(jobtype, geom = geometries, multi = multithreading, gau = gaussian)
    nos = numberofstructures(geom)
    cd("tmp")
    if jobtype == "Vc" || jobtype == "Gcav"
        Threads.@threads for i in 1:nos
            run(`$gau structure-$i-$jobtype.gjf`)
        end
    elseif jobtype == "Ger" && multi == "on"
        Threads.@threads for i in 1:nos
            run(`$gau structure-$i-$jobtype.gjf`)
        end
    elseif jobtype == "Ger" && multi == "off"
        for i in 1:nos
            run(`$gau structure-$i-$jobtype.gjf`)
        end
    end
    cd("..")
end


#------------------------------------------------------------------------------
# algebra.jl
#------------------------------------------------------------------------------
# linear scaling factor 𝑠, as the cubic root of the volume scaling
function calculate𝑠(𝑉𝑐 = 𝑉𝑐)  # 𝑉𝑐 is a global variable
    #𝑉𝑐 = get𝑉𝑐()
    # @. is a macro for vectoried operation/broadcasting
    return @. cbrt(𝑉𝑐/𝑉𝑐[:,1])    # nos * a 2D array
end


# average of 𝑠 over all structures at the same scalingfactor 𝑓
function average𝑠()
    𝑠 = calculate𝑠()
    # the mean function from the Staticstics package
    return mean(𝑠, dims=1)    # 1D array of length a
end


# dielectric permitivity 𝜀 = 1 + (𝜀₀-1)/𝑠̄³
function calculate𝜀()
    𝜀₀ = solventparameters()[1]
    𝑠̄ = average𝑠()
    return @. 1 + (𝜀₀ - 1) / 𝑠̄^3    # 1D array of length a
end


# Pauli repulson barrier 𝑍 = 𝑍₀/𝑠̄⁽³⁺𝜂⁾, where 𝑍₀ = 𝜌₀
function calculate𝑍(𝜂=𝜂)
    𝜌₀ = solventparameters()[2]
    𝑠̄ = average𝑠()
    return @. 𝜌₀ / 𝑠̄^(3+𝜂)    # 1D array of length a
end


# molar volume of solvent 𝑉ₘ = (𝑀/𝜌₀) * 𝑠̄³
function calculate𝑉ₘ()
    (𝜌₀, 𝑀) = solventparameters()[2:3]
    𝑠̄ = average𝑠()
    return @. (𝑀/𝜌₀) * 𝑠̄^3    # 1D array of length a
end


# Murnaghan equation of state fitting for pressure 𝑝 calculation
# three methods - julia (default), python, or mathematica
function murnaghan_eos(ft = fitting)
    if ft == "julia"
        return juliafitting()
    elseif ft =="python"
        return pythonfitting()
    elseif ft == "mathematica"
        return mathematicafitting()
    end
end


# using LsqFit
function juliafitting(𝑉𝑐 = 𝑉𝑐, 𝐺𝑒𝑟 = 𝐺𝑒𝑟, geom = geometries)
    nos = numberofstructures(geom)
    abc_parameters = Array{Float64}(undef, nos,3)
    Threads.@threads for i in 1:nos
    # python: y = (a/b)*(1/x)**b+(a-c)*x; y is Ger, x is Vc
    # mathematica: a*x ((1/b)*(t[[1, 1]]/x)^(b + 1) + 1) - c*x
    # a=p[1], b=p[2], c=p[3], x is Vc
        @. model(x, p) = (p[1]/p[2])*x^(-p[2]) + (p[1]-p[3])*x
        xdata = 𝑉𝑐[i,:] ./ 𝑉𝑐[i,1]
        ydata = 𝐺𝑒𝑟[i,:] .- 𝐺𝑒𝑟[i,1]
        p0 = [0.0, 5.0, 0.0]
        fit = curve_fit(model, xdata, ydata, p0)
        abc_parameters[i,1] = fit.param[1]/𝑉𝑐[i,1]
        abc_parameters[i,2] = fit.param[2]
        abc_parameters[i,3] = fit.param[3]/𝑉𝑐[i,1]
    end
    return abc_parameters    # nos * 3 2D array
end


# Using the lmfit package (to install: `pip install lmfit`)
# Make sure to use the correct python verion and
# the path to python inside the function is correct.
function pythonfitting(𝑉𝑐 = 𝑉𝑐, 𝐺𝑒𝑟 = 𝐺𝑒𝑟, geom = geometries)
    nos = numberofstructures(geom)
    #𝑉𝑐 = get𝑉𝑐()
    #𝐺𝑒𝑟 = get𝐺𝑒𝑟()
    abc_parameters = Array{Float64}(undef, nos,3)
    Threads.@threads for i in 1:nos
        script = """#!/usr/bin/env python
            #<examples/doc_model1.py>
            from numpy import sqrt, pi, exp, linspace, loadtxt, savetxt
            from lmfit import Model
            
            ### import matplotlib.pyplot as plt
            
            #data = loadtxt('Vc-Ger.dat')
            #x = data[:, 0]/data[0,0]
            #y = data[:, 1]-data[0,1]
            x = $(𝑉𝑐[i,:] ./ 𝑉𝑐[i,1])
            y = $(𝐺𝑒𝑟[i,:] .- 𝐺𝑒𝑟[i,1])
            
            def murnaghan(x, a, b, c):
                "1-d gaussian: murnaghan(x, a,b, c)"
                return (a/b)*(1/x)**b+(a-c)*x
            
            gmod = Model(murnaghan)
            result = gmod.fit(y, x=x, a=0, b=5, c=0)
            
            print(result.fit_report())
            
            #<end examples/doc_model1.py>"""
        
        # write fitting results to structure-$i-fitting.out files
        results = read(pipeline(`echo $script`, `/scratch/bochen/Python-2.7.18/bin/python`), String)
        #results = read(pipeline(`echo $script`, `python3`), String)
        open("tmp/structure-$i-fitting.out", "w") do file
            write(file, results)
        end
        
        # use an awk script to extract the a, b, c parameters
        awkscript = raw"/a:/ {a=$2}; /b:/ {b=$2}; /c:/ {c=$2}; END {print a,b,c}"
        abc = read(`awk $awkscript "tmp/structure-$i-fitting.out"`, String)
        abc_parameters[i,1] = parse(Float64, split(abc)[1])/𝑉𝑐[i,1]
        abc_parameters[i,2] = parse(Float64, split(abc)[2])
        abc_parameters[i,3] = parse(Float64, split(abc)[3])/𝑉𝑐[i,1]
    end
    return abc_parameters   # nos * 3 2D array
end


#= Not working; to be fully implemented.
function mathematicafitting(𝑉𝑐 = 𝑉𝑐, 𝐺𝑒𝑟 = 𝐺𝑒𝑟, geom = geometries)
    nos = numberofstructures(geom)
    #𝑉𝑐 = get𝑉𝑐()
    #𝐺𝑒𝑟 = get𝐺𝑒𝑟()
    abc_parameters = Array{Float64}(undef, nos,3)
    Threads.@threads for i in 1:nos
        # write data to Vc-Ger.dat file which is then read by Mathematica
        # the DelimitedFiles package is used
        open("tmp/Vc-Ger.dat", "w") do file
            writedlm(file, [𝑉𝑐[i,:] 𝐺𝑒𝑟[i,:]])
        end
        script = """###!/usr/local/bin/WolframScript -script

                    (* generate high-precision samples of a mixed distribution *)

                    t = Import["tmp/Vc-Ger.dat", "Table"]

                    Print[murnaghan-eos]

                    << NonlinearRegression`

                    Print[NonlinearRegress[t, 
                    t[[1, 2]] + a*x ((1/b)*(t[[1, 1]]/x)^(b + 1) + 1) - c*x, {a, b, 
                    c}, x, RegressionReport -> {BestFitParameters, ParameterCITable, 
                    ANOVATable, BestFit, PredictedResponse}]]"""
        results = read(`math -script $script`, String)
        #results = read(pipeline(`echo $script`, `math`), String)
        open("tmp/structure-$i-fitting.out", "w") do file
            write(file, results)
        end
        # extraction of abc parameters to be implemented
    end
    rm("tmp/Vc-Ger.dat")
    return abc_parameters    # 2D array of dimension nos * 3
end
=#


function calculate𝑝(𝑉𝑐 = 𝑉𝑐)
    #𝑉𝑐 = get𝑉𝑐()    # nos * a 2D array
    abc = murnaghan_eos()    # nos * 3 2D array
    𝑎 = abc[:,1]   # 1D array of length nos
    𝑏 = abc[:,2]
    𝑐 = abc[:,3]
    # $p[$n_fact]=($a * (($volume[1] / $volume[$n_fact]) ** ($b + 1) -1 ) + $c)*$Eh_o_Ang3_to_GPa;
    #𝑝[i,j] = (𝑎[i] * ((𝑉𝑐[i,1] / 𝑉𝑐[i,j]) ^ (𝑏[i] + 1) - 1) + 𝑐[i]) * 4359.74417
    #𝑝[i,j] = (𝑎[i] * (1 - (𝑉𝑐[i,1] / 𝑉𝑐[i,j]) ^ 𝑏[i] ) + 𝑐[i]) * 4359.74417
    
    # nos * a 2D array; 1 hartree/Å³ = 4359.74417 GPa
    return @. (𝑎 * ( (𝑉𝑐[:,1]/𝑉𝑐)^(𝑏+1) - 1 ) + 𝑐) * 4359.74417 
    #return @. (𝑎 * (1 - (𝑉𝑐[:,1]/𝑉𝑐)^𝑏) + 𝑐) * 4359.74417
end


# average of 𝑝 over all structures at the same scalingfactor 𝑓
function average𝑝()
    𝑝 = calculate𝑝()    # nos * a 2D array
    return mean(𝑝, dims=1)   # 1 * a 2D array
end


function calculate𝐺𝑐𝑎𝑣(𝐸𝑐𝑎𝑣 = 𝐸𝑐𝑎𝑣, 𝑉𝑐𝑎𝑣 = 𝑉𝑐𝑎𝑣)
    𝑝̄ = average𝑝()      # 1 * a 2D array
    #𝐸𝑐𝑎𝑣 and 𝑉𝑐𝑎𝑣 are nos * a 2D arrays; 1 GPa*Å³ = 2.293712569e-4 Hartree
    return @. 𝐸𝑐𝑎𝑣 + 𝑝̄ * 𝑉𝑐𝑎𝑣 * 2.293712569e-4    # nos * a 2D array
end


function calculate𝐺𝑡𝑜𝑡(𝐺𝑒𝑟 = 𝐺𝑒𝑟)
    return 𝐺𝑒𝑟 .+ calculate𝐺𝑐𝑎𝑣()    # nos * a 2D array
end


function calculateΔ𝐺𝑡𝑜𝑡()
    𝐺𝑡𝑜𝑡 = calculate𝐺𝑡𝑜𝑡()
    Δ𝐺𝑡𝑜𝑡 = Array{Float64}(undef, size(𝐺𝑡𝑜𝑡))  # nos * a 2D array
    for i in 1:length(𝐺𝑡𝑜𝑡[1,:])
        @. Δ𝐺𝑡𝑜𝑡[:,i] = (𝐺𝑡𝑜𝑡[:,i] - 𝐺𝑡𝑜𝑡[1,i]) * 627.509  # 1 hartree = 627.509 kcal/mol
    end
    return Δ𝐺𝑡𝑜𝑡
end

#= need to figure out how to calculate barrier; use the TS structure or the maximum?
function calculateΔ𝑉activation()
    𝑝̄ = average𝑝()    # 1 * a 2D array
    Δ𝐺𝑡𝑜𝑡 = calculateΔ𝐺𝑡𝑜𝑡()   # nos * a 2D array
    Δ𝐺𝑡𝑜𝑡activation = Δ𝐺𝑡𝑜𝑡[50,:]    # 1D array of a
    slope, intercept = [ones(length(𝑝̄)) vec(𝑝̄)] \ Δ𝐺𝑡𝑜𝑡activation
    return slope * 4.184    # 1 kcal mol⁻¹ / GPa = 4.184 cm³/mol; 4.184 * 10^3 / 10^9 * 10^6
end
=#

# restart from previously interupted Ger jobs
# The idea is to first generate a set of input filenames, and remove those finished
# from the set. Then restart Ger jobs for the rest of filenames in the set.
function restartger(geom=geometries, 𝑓 = scalingfactors, multi = multithreading)
    nos = numberofstructures(geom)
    a = length(𝑓)
    all = [1:nos;]    # Int64 array containing all job numbers
    
    cd("tmp")
    
    # awk script to return the finished job numbers
    script = raw"for file in *Ger.log; do i=${file#*structure-}; i=${i%-Ger.log}; grep 'SCF Done' $file | wc -l; echo $i; done | paste - - | awk '/" * "$a" * raw"/ {print $2}'"
    open("restart.sh", "w") do file
        write(file, "$script")
    end
    
    string = read(`bash restart.sh`, String)
    rm("restart.sh")
    
    finished = parse.(Int64, split(string))    # the finished job numbers
    unfinished = setdiff(all, finished)  # remove the finished job numbers from all
    
    if multi == "on"
        Threads.@threads for i in unfinished
            run(`g16 structure-$i-Ger.gjf`)
        end
    elseif multi == "off"
        for i in unfinished
            run(`g16 structure-$i-Ger.gjf`)
        end
    end
    cd("..")
end

#------------------------------------------------------------------------------
# main procedure
#------------------------------------------------------------------------------
#function main(restart = "no")
const gaussian = gaussianversion()
    if restart == "no"
        writegjf("Vc")
        rungaussian("Vc")
        const 𝑉𝑐 = get𝑉𝑐()
        writegjf("Ger")          # write .gjf files for cavity volume "Ger" calculation 
        rungaussian("Ger")
    elseif restart == "yes"
        const 𝑉𝑐 = get𝑉𝑐()
        restartger()
    else
        println("restart only accepts \"yes\" or \"no\"")
    end

    const 𝐺𝑒𝑟 = get𝐺𝑒𝑟()

    writegjf("Gcav")         # write .gjf files for cavitation energy "Gcav" calculation 

    rungaussian("Gcav")      # run Gaussian jobs

    const (𝑉𝑐𝑎𝑣,𝐸𝑐𝑎𝑣) = get𝑉𝑐𝑎𝑣𝐸𝑐𝑎𝑣()

    writeproperties()        # write properties.dat file
#end
#------------------------------------------------------------------------------
# Step 1: cavity volume 𝑉𝑐(𝑓) Gaussian jobs and solvent property calculations
#------------------------------------------------------------------------------

#writegjf("Vc")          # write .gjf files for cavity volume "Vc" calculation 

#rungaussian("Vc")       # run Gaussian jobs

#const 𝑉𝑐 = get𝑉𝑐()       # extract cavity volume 𝑉𝑐 from Gaussian output

#calculate𝑠()            # calculate linear scaling 𝑠 from 𝑉𝑐 date

#average𝑠()              # calculated the average of 𝑠 over all structures

#calculate𝜀()            # calculate dielectric permitivity 𝜀

#calculate𝑍()            # calculate Pauli repulsion barrier 𝑍

#calculate𝑉ₘ()           # calculate the molar volume of solvent 𝑉ₘ

#------------------------------------------------------------------------------
# Step 2: electronic structure Gaussian jobs and pressure calculations 
#------------------------------------------------------------------------------

#    writegjf("Ger")      # write .gjf files for cavity volume "Ger" calculation

#    rungaussian("Ger")   # run Gaussian jobs

#    const 𝐺𝑒𝑟 = get𝐺𝑒𝑟()      # extract 𝐺𝑒𝑟 from Gaussian output

#calculate𝑝()             # calculate pressure 𝑝

#------------------------------------------------------------------------------
# Step 3: cavitation energy Gaussian jobs
#------------------------------------------------------------------------------

#writegjf("Gcav")         # write .gjf files for cavitation energy "Gcav" calculation 

#rungaussian("Gcav")      # run Gaussian jobs

#const 𝑉𝑐𝑎𝑣 = get𝑉𝑐𝑎𝑣𝐸𝑐𝑎𝑣()[1]
#const 𝐸𝑐𝑎𝑣 = get𝑉𝑐𝑎𝑣𝐸𝑐𝑎𝑣()[2]    # extract 𝐺𝑐𝑎𝑣 from Gaussian output

#calculate𝐺𝑐𝑎𝑣()           # calculate cavitation energy 𝐺𝑐𝑎𝑣

#calculate𝐺𝑡𝑜𝑡()            # calculate total energy 𝐺𝑡𝑜𝑡

#------------------------------------------------------------------------------
# print results
#------------------------------------------------------------------------------

#writeproperties()        # write properties.dat file
