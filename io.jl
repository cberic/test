function writegjf(jobtype)
    if jobtype == "Vc"
        gjfvc()
    elseif jobtype == "Ger"
        gjfger()
    elseif jobtype == "Gcav"
        gjfgcav()
    end
end


function gjfvc(geom = geometries, 𝑓 = scalingfactors)
    a = length(𝑓)
    st = structures(geom)
    nos = numberofstructures(geom)
    noa = numberofatoms(geom)
    atoms = atomlist(geom)
    lines = coordinatelines(geom)
    sp = solventparameters()
    r = atomicradii()
    
    run(`mkdir -p tmp`)  # make a tmp dir in the working folder

    # job killed after L301; 
    Threads.@threads for i in 1:nos  # loop through structures; use nproc=1 and multithreading
        for j in 1:a   # loop through scaling factors
            # writing mode for the first scalingfactor and appending mode for the rest
            open("tmp/structure-$i-Vc.gjf", "$(j == 1 ? "w" : "a")") do file 
                println(file, "%kjob l301")
                println(file, "%chk=structure-$i-Vc-$(𝑓[j]).chk")
                println(file, "%nproc=1")
                println(file, "%mem=200mb")
                println(file, "#p $method $basis")
                println(file, "# scrf=(iefpcm,solvent=$solvent,read) nosym guess=only pop=none")
                println(file)
                println(file, "title")
                println(file)
                println(file, "$charge $multiplicity")
                println(file, "$(st[i])")
                println(file)
                println(file, "qrep pcmdoc geomview nodis nocav g03defaults tsare=$tesserae")
                println(file, "nsfe=$noa")
                println(file, "nvesolv=$(sp[4]) solvmw=$(sp[3]) rsolv=$(sp[5])")
                println(file, "eps=$(sp[1]) rhos=$(sp[2])")         
                println(file)
                
                for k in 1:noa
                    println(file, " $(lines[i][k])  $(r[atoms[k]])  $(𝑓[j])")
                end
                println(file)
                if j != a
                    println(file, "--link1--")
                end
            end
        end
    end
end


function gjfger(geom = geometries, 𝑓 = scalingfactors)
    a = length(𝑓)
    st = structures(geom)
    nos = numberofstructures(geom)
    noa = numberofatoms(geom)
    atoms = atomlist(geom)
    lines = coordinatelines(geom)
    sp = solventparameters()
    r = atomicradii()
    𝜀 = calculate𝜀()    # needs data of 𝜀 and 𝜌 at different scalingfactors
    𝑍 = calculate𝑍()
    # job killed after L502
    Threads.@threads for i in 1:nos  # loop through structures
        for j in 1:a   # loop through scaling factors
            # writing mode for the first scalingfactor and appending mode for the rest
            open("tmp/structure-$i-Ger.gjf", "$(j == 1 ? "w" : "a")") do file
                write(file, """
                %kjob l502
                %chk=structure-$i-Ger-$(𝑓[j]).chk
                %nproc=$nproc
                %mem=$mem
                #p $method $basis
                # scrf=(iefpcm,solvent=$solvent,read) nosym 6d 10f

                title

                $charge $multiplicity
                $(st[i])

                qrep pcmdoc geomview nodis nocav g03defaults tsare=$tesserae
                nsfe=$noa
                nvesolv=$(sp[4]) solvmw=$(sp[3]) rsolv=$(sp[5])
                eps=$(𝜀[j]) rhos=$(𝑍[j])                

                """)

                for k in 1:noa
                    write(file, " $(lines[i][k])  $(r[atoms[k]])  $(𝑓[j])\n")
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


function gjfgcav(cav = cavity, geom = geometries, 𝑓 = scalingfactors)
    a = length(𝑓)
    st = structures(geom)
    nos = numberofstructures(geom)
    noa = numberofatoms(geom)
    atoms = atomlist(geom)
    #lines = coordinatelines(geom)
    sp = solventparameters()
    r = atomicradii()
    𝑉ₘ = calculate𝑉ₘ()    # needs data of 𝜀 and 𝜌 at different scalingfactors

    # job killed after L301
    Threads.@threads for i in 1:nos  # loop through structures
        for j in 1:a   # loop through scaling factors
            # writing mode for the first scalingfactor and appending mode for the rest
            open("tmp/structure-$i-Gcav.gjf", "$(j == 1 ? "w" : "a")") do file
                write(file, """
                %kjob l301
                %chk=structure-$i-Gcav-$(𝑓[j]).chk
                %nproc=$nproc
                %mem=$mem
                #p $method $basis
                # scrf=(iefpcm,solvent=$solvent,read) nosym guess=only pop=none

                title

                $charge $multiplicity
                $(st[i])

                norep pcmdoc geomview nodis cav g03defaults $(cav == "vdw" ? "noaddsph" : "" ) tsare=$tesserae
                nsfe=$noa
                Vmol=$(𝑉ₘ[j]) rsolv=$(sp[5])

                """)

                for k in 1:noa
                    write(file, " $k  $(r[atoms[k]] * 𝑓[1])  1.0\n")
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


# extract volume 𝑉𝑐 from gaussian output files 
# 𝑉𝑐 as a function of scaling factor 𝑓
function get𝑉𝑐(geom = geometries, 𝑓 = scalingfactors)
    nos = numberofstructures(geom)
    a = length(𝑓)
    𝑉𝑐 = Array{Float64}(undef, nos,a)    # 2D array of nos and a
    for i in 1:nos
        j = 1    # j indexes the length(𝑓)
        open("tmp/structure-$i-Vc.log") do file
            for line in eachline(file)
                if occursin("GePol: Cavity volume", line)
                    𝑉𝑐[i,j] = parse(Float64, split(line)[5])
                    j += 1    # j indexes the length(𝑓)
                end 
            end
        end
    end
    return 𝑉𝑐
end


# extract electronic energy 𝐺𝑒𝑟 from gaussian output files
function get𝐺𝑒𝑟(geom = geometries, 𝑓 = scalingfactors)
    nos = numberofstructures(geom)
    a = length(𝑓)
    𝐺𝑒𝑟 = Array{Float64}(undef, nos,a)    # 2D array of nos and a
    for i in 1:nos
        j = 1    # j should index the length(𝑓)
        open("tmp/structure-$i-Ger.log") do file
            for line in eachline(file)
                if occursin("SCF Done", line)
                    𝐺𝑒𝑟[i,j] = parse(Float64, split(line)[5])
                    j += 1    # j should index the length(𝑓)
                end
            end
        end
    end
    return 𝐺𝑒𝑟
end


# extract non-electrostatic cavitation energy 𝐸𝑐𝑎𝑣 from gaussian output files
function get𝐸𝑐𝑎𝑣(geom = geometries, 𝑓 = scalingfactors)
    nos = numberofstructures(geom)
    a = length(𝑓)
    𝐸𝑐𝑎𝑣 = Array{Float64}(undef, nos,a)    # 2D array of nos and a
    for i in 1:nos
        j = 1    # j should index the length(𝑓)
        open("tmp/structure-$i-Gcav.log") do file
            for line in eachline(file)
                if occursin("PCM non-electrostatic energy", line)
                    𝐸𝑐𝑎𝑣[i,j] = parse(Float64, split(line)[5])
                    j += 1    # j should index the length(𝑓)
                end
            end
        end
    end
    return 𝐸𝑐𝑎𝑣
end


# write the properties.dat file; the Printf package is used
using Printf
function writeproperties(geom = geometries, 𝑓 = scalingfactors)
    nos = numberofstructures(geom)
    a = length(𝑓)
    𝑉𝑐 = get𝑉𝑐()
    𝑠 = calculate𝑠()
    𝑠̄ = average𝑠()
    𝜀 = calculate𝜀()
    𝑍 = calculate𝑍()
    𝑉ₘ = calculate𝑉ₘ()
    𝐺𝑒𝑟 = get𝐺𝑒𝑟()
    𝑝 = calculate𝑝()
    𝑝̄ = average𝑝()
    𝐸𝑐𝑎𝑣 = get𝐸𝑐𝑎𝑣()
    𝐺𝑐𝑎𝑣 = calculate𝐺𝑐𝑎𝑣()
    𝐺𝑡𝑜𝑡 = calculate𝐺𝑡𝑜𝑡()
    open("properties.dat", "w") do file
        for i in 1:nos
            write(file, "structure $i\n")
            write(file, "#    𝑓       𝑉𝑐(𝑓)      𝑠(𝑓)         𝑠̄(𝑓,𝑛ₛ)      𝜀(𝑠̄)        𝑍(𝑠̄)        𝑉ₘ(𝑠̄)      𝐺𝑒𝑟(𝑓)            𝑝(𝑓)      𝑝̄(𝑠̄)      𝐸𝑐𝑎𝑣(𝑓)        𝐺𝑐𝑎𝑣(𝑓)        𝐺𝑡𝑜𝑡(𝑓)\n")
            for j in 1:a
                @printf(file, "%d    %.2f    %7.3f    %.6f    %.6f    %.6f    %.6f    %7.3f    %.8f    %7.3f    %6.3f    %.8f      %.8f      %.8f\n", 
                               j,    𝑓[j],   𝑉𝑐[i,j], 𝑠[i,j],  𝑠̄[j],   𝜀[j],   𝑍[j],  𝑉ₘ[j], 𝐺𝑒𝑟[i,j], 𝑝[i,j], 𝑝̄[j], 𝐸𝑐𝑎𝑣[i,j], 𝐺𝑐𝑎𝑣[i,j], 𝐺𝑡𝑜𝑡[i,j])
            end
            write(file, "\n")
        end
    end
end