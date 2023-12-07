
# define necessary functions
function integral2(t, x)
    I = zeros(length(x))
    for i in 2:length(x)
        I[i] = I[i-1] + (x[i-1] + x[i])/2*(t[i]-t[i-1])
    end
    ddxdt2 = (I[3]-I[2])/(t[3]-t[2])
    I[1] = I[2] - ddxdt2*(t[2]-t[1])
    return I
end
    
function integral2(t, x, bounds)
    _, it1 = findmin(abs.(t .- bounds[1]))
    _, it2 = findmin(abs.(t .- bounds[2]))
    I = integral2(t, x)
    return I[it2] - I[it1]
end
    
#using Roots
#using SpecialFunctions
    
    # system parameters
#dPdz = -2.0
#rho = 10.0 #Moved to main.jl
mu = rho*nu
r1 = R1
r2 = R2

    u0 = (r)->0 #(r-r1)*(r2-r)
    
    # steady state solution
    c1 = dPdz*(r1^2-r2^2)/(4*mu*log(r2/r1))
    c2 = -dPdz*r1^2/4/mu - c1*log(r1)
    Vz_ss = (r)->dPdz*r^2/4/mu + c1*log(r) + c2
    
    # transient solution
    rs = range(r1, stop = r2, length = integralresolution)
    J0 = besselj0
    Y0 = bessely0
    char_eq2 = (x)->J0(x*r1)*Y0(x*r2)-J0(x*r2)*Y0(x*r1)
    lambda = Roots.find_zeros(char_eq2, [0, 100])
    R = (n)->(r)->J0(lambda[n]*r)/J0(lambda[n]*r2) - Y0(lambda[n]*r)/Y0(lambda[n]*r2)
    T = (n)->(t)->exp(-mu/rho*lambda[n]^2*t)
    norm = (n)->integral2(rs, rs.*R(n).(rs).^2, [r1,r2])
    c = (n)->integral2(rs, (u0.(rs)-Vz_ss.(rs)).*R(n).(rs).*rs, [r1,r2])/norm(n)
    #Added precompute loop
    n = 7 #num eigs
    cs = []
    Rs = []
    Ts = []
    progressbar("Precomputing Vz",0)
    for i in 1:n
        push!(cs,integral2(rs, (u0.(rs)-Vz_ss.(rs)).*R(i).(rs).*rs, [r1,r2])/integral2(rs, rs.*R(i).(rs).^2, [r1,r2]))
        push!(Rs,(r)->J0(lambda[i]*r)/J0(lambda[i]*r2) - Y0(lambda[i]*r)/Y0(lambda[i]*r2))
        push!(Ts,(t)->exp(-mu/rho*lambda[i]^2*t))
        progressbar("Precomputing Vz",i/n)
    end
    Vz_t = (n)->(r,t)->begin
    val = 0.0
    for i in 1:n
        Rn = Rs[i]
        Tn = Ts[i]
        val += cs[i]*Rn(r)*Tn(t)#c(i)*R(i)(r)*T(i)(t)
    end
    return val
    end
    
    # complete solution
    Vz2 = (r,t)-> Vz_ss(r) + Vz_t(n)(r,t)