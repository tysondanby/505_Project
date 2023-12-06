using Plots,SpecialFunctions, Roots, ForwardDiff,Printf
#Problem Parameters--------------------------------------------------------------------------------------------------------------------------------
R1 = 1
R2 = 2
nu = .1
function omega(t)
    return 0.1*sin(2*pi*t)#.1*t
end
tspan = (0,3)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Visualization Parameters----------------------------------------------------------------------------------------------------------------------------------------------------------
nframes = 90
number_eigenvalues = 30
integralresolution = 10000
Rresolution = 250
FPS = 30
Vθ_3Ddfilename = "3Danimation.gif"
Vθ_2Dfilename = "2Danimation.gif"
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

function plotanulus(rs,vs,time,vθlims)

    #PARAMS
    θresolution = 30
    coldv = vθlims[1]
    hotv = vθlims[2]
    #
    xs = []
    ys = []
    fs=[]
    lims = (-1.5*rs[end],1.5*rs[end])
    zrange = (vθlims[1]*2,vθlims[2]*2)
    
    p=plot(xs,ys,fs,xlims = lims,ylims = lims,zlims = zrange,linecolor = RGBA(1,0,0,1),legend = false,camera = (15,50),plot_title="t=$(@sprintf("%.2f", time))"*"s",plot_titlelocation = :left)

    #Plot rings
    θs = collect(0:2*pi/θresolution:2*pi)
    for i = 1:1:length(rs)
        
        colorR = (vs[i]-coldv)/(hotv-coldv)
        colorB = 1-colorR
        newxs = @.sin(θs)*rs[i]
        newys = @.cos(θs)*rs[i]
        newfs = @.cos(θs)/cos(θs) *vs[i]
        plot!(newxs,newys,newfs,linecolor = RGBA(colorR,0,colorB,1))
    end

    #plot spokes
    colorS = 15/length(θs)
    for i = 1:1:length(θs)
        newxs = @.sin(θs[i])*rs
        newys = @.cos(θs[i])*rs
        plot!(newxs,newys,vs,linecolor = RGBA(colorS,colorS,colorS,1))
    end
    return p
end

function integrate(f,lowb,upb,steps)#Midpoint method
    I = 0
    for i =1:2:(2*steps-1)
        x=lowb + i*(upb-lowb)/(steps*2)
        I = I + f(x)*(upb-lowb)/steps
    end
    return I
end

function eigFun(eig,r)
    return besselj1(eig*r)/besselj1(eig*R2) - bessely1(eig*r)/bessely1(eig*R2)#besselj1(eig*r)*bessely1(eig*R2)-besselj1(eig*R2)*bessely1(eig*r)#
end

function char_eq(mu)
    return besselj1(mu*R1)*bessely1(mu*R2)-besselj1(mu*R2)*bessely1(mu*R1)
end

function findeigs(nEigs)
    #PARAMS
    step = 0.01 #rootfinding step
    #
    eigs = []
    lo = 0
    hi = step
    while length(eigs) < nEigs
        newzeros=find_zeros(char_eq, lo,hi)
        lo = lo + step
        hi = hi + step
        eigs = [eigs;newzeros]
    end
    #println(eigs)
    yns = []
    ynprimeR2s = zeros(length(eigs))
    sqrnorms = zeros(length(eigs))
    for i = 1:1:length(eigs)
        newyn = function yn(x)
                    return eigFun(eigs[i],x)
                end
        push!(yns,newyn)
        ynprimeR2s[i] = ForwardDiff.derivative(newyn, R2)
        function sqrnormint(x)
            return newyn(x)^2*x
        end
        sqrnorms[i] = integrate(sqrnormint,R1,R2,integralresolution)
    end

    return eigs,yns,ynprimeR2s,sqrnorms
end

eigs,yns,ynprimeR2s,sqrnorms = findeigs(number_eigenvalues)



function Vth(r,t)

    vth = 0
    for i = 1:1:length(eigs)
        yn = yns[i]
        function vthbarint(tau)
            return -1*omega(t-tau)*nu*R2^2*ynprimeR2s[i]*exp(-nu*eigs[i]^2*tau)
        end
        vthbar = integrate(vthbarint,0,t,integralresolution)
        #function test(x)
        #    return x*yn(x)*f(x)
        #end
        #testint = integrate(test,R1,R2,100*i)
        #vth=vth + (testint/sqrnorm)*yn(r) #TEST
        vth=vth + (vthbar/sqrnorms[i])*yn(r)
    end
    return vth
end

function Vth(r::Vector,t)

    vth = zeros(length(r))
    for i = 1:1:length(eigs)
        yn = yns[i]
        function vthbarint(tau)
            return -1*omega(t-tau)*nu*R2^2*ynprimeR2s[i]*exp(-nu*eigs[i]^2*tau)
        end
        vthbar = integrate(vthbarint,0,t,integralresolution)

        for j = 1:1:length(r)
            vth[j]=vth[j] + (vthbar/sqrnorms[i])*yn(r[j])
        end
    end
    return vth
end

function simulate()
    rs = collect(R1:(R2-R1)/(Rresolution-1):R2)
    ts = collect(tspan[1]:(tspan[2]-tspan[1])/(nframes-1):tspan[2])
    vθs = []
    minvθ = 0
    maxvθ = 0
    for frame = 1:1:nframes
        vs = zeros(length(rs))
        vs = Vth(rs,ts[frame])
        if maximum(vs) > maxvθ
            maxvθ = maximum(vs)
        end
        if minimum(vs) < minvθ
            minvθ = minimum(vs)
        end
        push!(vθs,vs)
    end
    return rs,ts,vθs,minvθ,maxvθ
end

rs,ts,vθs,minvθ,maxvθ = simulate()
#animate
vθanim3d = @animate for frame = 1:1:nframes
    plotanulus(rs,vθs[frame],ts[frame],(minvθ,maxvθ))
end
vθanim2d = @animate for frame = 1:1:nframes
    plot(rs,vθs[frame],xlims = (0,R2),ylims=(minvθ,maxvθ),legend = false,plot_title="t=$(@sprintf("%.2f", ts[frame]))"*"s",plot_titlelocation = :left)
end


gif(vθanim3d, Vθ_3Ddfilename, fps=FPS)
gif(vθanim2d, Vθ_2Ddfilename, fps=FPS)

