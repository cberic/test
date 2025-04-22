using GLMakie,LinearAlgebra, DelaunayTriangulation, PlyIO

# Define van der Waals radii for common elements
VDW_RADII = Dict(
    "H" => 1.2, "C" => 1.7, "N" => 1.55, "O" => 1.52, "S" => 1.8
)

# Step 1: Load molecular structure (atoms and coordinates)
struct Atom
    element::String
    position::Vector{Float64}
    radius::Float64
end

function load_molecule(atoms::Vector{Tuple{String,Vector{Float64}}})
    return [Atom(a[1], a[2], get(VDW_RADII, a[1], 1.5)) for a in atoms]
end

# Step 2: Generate covering points using Fibonacci lattice
function fibonacci_sphere(atom::Atom, num_points::Int)
    phi = (1 + sqrt(5)) / 2  # Golden ratio
    points = Vector{Vector{Float64}}()
    for i in 0:num_points-1
        theta = 2π * (i / phi)
        z = 1 - (2 * i / (num_points - 1))
        r = sqrt(1 - z^2)
        x = r * cos(theta)
        y = r * sin(theta)
        push!(points, atom.position .+ atom.radius .* [x, y, z])
    end
    return points
end

# Step 3: Identify exposed surface points
function filter_exposed_points(atoms::Vector{Atom}, points::Vector{Vector{Float64}})
    exposed = Vector{Vector{Float64}}()
    for p in points
        is_exposed = true
        for atom in atoms
            if norm(p - atom.position) < atom.radius - 0.01  # Small tolerance
                is_exposed = false
                break
            end
        end
        if is_exposed
            push!(exposed, p)
        end
    end
    return exposed
end

# Step 4: Perform Delaunay triangulation on exposed points
function triangulate_surface(points::Vector{Vector{Float64}})
    points_matrix = hcat(points...)'  # Convert to Nx3 matrix
    triangulation = triangulate(points_matrix)
    return triangulation
end

# Step 5: Output the surface as a PLY file
function save_to_ply(filename::String, points::Vector{Vector{Float64}}, triangles)
    # Convert points to the correct dictionary format
    vertices = [Dict("x" => p[1], "y" => p[2], "z" => p[3]) for p in points]
    
    # Convert triangles to a list of dictionaries
    faces = [Dict("vertex_indices" => collect(tri .- 1)) for tri in triangles]  # Convert to 0-based indexing

    # Create PLY structure
    ArrayProperty()
    ply = PlyIO.Ply(ASCII=true, elements=[
        PlyIO.PlyElement("vertex", vertices),
        PlyIO.PlyElement("face", faces)
    ])
    #rand() randn()
    # Save to file
    PlyIO.save_ply(ply, filename)
end

# Function to generate sphere mesh
function sphere(center, radius; res=30)
    θ = range(0, 2π, length=res)
    ϕ = range(0, π, length=res)
    
    x = [radius * sin(ϕ[j]) * cos(θ[i]) + center[1] for i in 1:res, j in 1:res]
    y = [radius * sin(ϕ[j]) * sin(θ[i]) + center[2] for i in 1:res, j in 1:res]
    z = [radius * cos(ϕ[j]) + center[3] for i in 1:res, j in 1:res]
    
    return x, y, z
end

# Function to visualize vdW spheres and points
function visualize_vdw(atoms, all_points)
    fig = Figure(resolution=(800, 600))
    ax = Axis3(fig[1, 1], title="vdW Spheres & Points", aspect=:data)

    # Plot vdW spheres
    for atom in atoms
        x, y, z = sphere(atom.position, atom.radius)
        surface!(ax, x, y, z, color=:lightblue, transparency=true)
    end

    # Plot all generated points
    xs = [p[1] for p in all_points]
    ys = [p[2] for p in all_points]
    zs = [p[3] for p in all_points]
    scatter!(ax, xs, ys, zs, markersize=3, color=:red)

    display(fig)
end

# Example usage
atoms = [
    ("C", [0.0, 0.0, 0.0]),
    ("O", [1.2, 0.0, 0.0]),
    ("H", [-0.6, 1.0, 0.0]),
    ("H", [-0.6, -1.0, 0.0])
]
molecule = load_molecule(atoms)

all_points = vcat([fibonacci_sphere(a, 1000) for a in molecule]...)
exposed_points = filter_exposed_points(molecule, all_points)

triangles = triangulate_surface(exposed_points)
#save_to_ply("vanderwaals_surface.ply", exposed_points, triangles)

println("Molecular vdW surface saved as 'vanderwaals_surface.ply'.")

# Visualize
visualize_vdw(molecule, exposed_points)