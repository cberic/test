function tidygeometries(geom = geometries)
    # remove trailing and tailing spaces/blank lines 
    # and add a space to the begininig
    a = " " * strip(geom)

    # change block separators to "NEXT"
    b = replace(a, r"\s*\n+\s*\n+\s*" => "NEXT")
    
    # remove trailing and tailing spaces within each coordinate line 
    # and add a space to the begininig of coordinate line
    c = replace(b, r"\s*\n\s*" => "\n ")
    
    # change block separators back to "\n\n" and add a space
    return replace(c, "NEXT" => "\n\n ")
end

#= the function below generates an array of blocks
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
=#
function structures(geom = geometries)
    return split(tidygeometries(geom), "\n\n", keepempty=false)
end

function numberofstructures(geom = geometries)
    return length(structures(geom))
end

#= the function below generates an array of arrays of lines
 C    1.10712900   -1.48465400    0.00000000
=#
function structurelines(geom = geometries)
    blocks = structures(geom)
    nos = numberofstructures(geom)
    lines = Array{Array{String}}(undef, nos)
    for i in 1:nos
        lines[i] = split(blocks[i], "\n", keepempty=false)
    end
    return lines
end

function numberofatoms(geom = geometries)
    return length(structurelines(geom)[1])
end

#= the function below generates an array of atoms (no coordinates)
 C
=#
function atomlist(geom = geometries)
    lines1 = structurelines(geom)[1]
    noa = numberofatoms(geom)
    atoms = Array{String}(undef, noa)
    for i in 1:noa
        atoms[i] = split(lines1[i], keepempty=false, limit=2)[1]
    end
    return atoms
end

#= the function below generates an array of arrays of coordinate lines (no atom label)
 1.10712900   -1.48465400    0.00000000
=#
function coordinatelines(geom = geometries)
    nos = numberofstructures(geom)
    lines = structurelines(geom)
    noa = numberofatoms(geom)
    for i in 1:nos
        for j in 1:noa
            lines[i][j] = " " * split(lines[i][j], keepempty=false, limit=2)[2]
        end
    end
    return lines
end

function coordinatelines2(geom = geometries)
    nos = numberofstructures(geom)
    lines = structurelines(geom)
    noa = numberofatoms(geom)
    array = Array{String}(undef, nos,noa)
    for i in 1:nos
        for j in 1:noa
            array[i,j] = " " * split(lines[i][j], keepempty=false, limit=2)[2]
        end
    end
    return array
end

# slower function with much more allocations
#function atomlist2(geometries)
#    noa = numberofatoms(geometries)
#    atoms = Array{String}(undef, noa) 
#    for i in 1:noa
#        atoms[i] = split(geometries)[4*i-3]
#    end
#    return atoms
#end


function atomcoordinates(geom = geometries)
    array = split(geom)

    nos = numberofstructures(geom)
    noa = numberofatoms(geom)

    # 3D array storing each coordinate
    coordinates = Array{String}(undef, nos,noa,3) 
    for i in 1:nos
        for j in 1:noa
            for k in 1:3
                coordinates[i,j,k] = array[(i-1)*40 + (j-1)*4 + k+1]
            end
        end
    end
    return coordinates
    #return coordinates[m,n,:]
    #return string(coordinates[m,n,1], " ", coordinates[m,n,2], " ", coordinates[m,n,3])
end