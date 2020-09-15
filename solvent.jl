# returns a 5-element tuple of solvent properties: 
# (dielectric constant 𝜀, density 𝜌, molar mass 𝑀, number of valence electrons, molecular radius (Ang.) )
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