using LsqFit

# nonlinear least squares fitting for the Murnaghan equation of state
function nonlinearlsq(𝑉𝑐 = 𝑉𝑐, 𝐺𝑒𝑟 = 𝐺𝑒𝑟, geom = geometries)
    nos = numberofstructures(geom)
    array = Array{Float64}(undef, nos,3)
    Threads.@threads for i in 1:nos
    # python: y = (a/b)*(1/x)**b+(a-c)*x; y is Ger, x is Vc
    # mathematica: a*x ((1/b)*(t[[1, 1]]/x)^(b + 1) + 1) - c*x
    # a=p[1], b=p[2], c=p[3], x is Vc
        @. model(x, p) = (p[1]/p[2])*x^(-p[2]) + (p[1]-p[3])*x
        xdata = 𝑉𝑐[i,:] ./ 𝑉𝑐[i,1]
        ydata = 𝐺𝑒𝑟[i,:] .- 𝐺𝑒𝑟[i,1]
        p0 = [0.0, 5.0, 0.0]
        fit = curve_fit(model, xdata, ydata, p0)
        array[i,1] = fit.param[1]/𝑉𝑐[i,1]
        array[i,2] = fit.param[2]
        array[i,3] = fit.param[3]/𝑉𝑐[i,1]
    end
    return array
end

# linear least squares fitting for Δ𝐺𝑡𝑜𝑡-𝑝̄ curve to obtain volume of activation
function linearlsq(xdata = 𝑝̄, ydata = Δ𝐺𝑡𝑜𝑡‡)
    # python: y = (a/b)*(1/x)**b+(a-c)*x; y is Ger, x is Vc
    # mathematica: a*x ((1/b)*(t[[1, 1]]/x)^(b + 1) + 1) - c*x
    # a=p[1], b=p[2], c=p[3], x is Vc
    @. model(x, p) = p[1] * x + p[2]
    p0 = [-1.0, 5.0]
    fit = curve_fit(model, xdata, ydata, p0)
end