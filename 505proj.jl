using Plots,SpecialFunctions, Roots, ForwardDiff,Printf
#Problem Parameters
R1 = 1
R2 = 2
nu = .1
tspan = (0,3)
nframes = 90
number_eigenvalues = 30
integralresolution = 10000
Rresolution = 250
FPS = 30
threedfilename = "3Danimation.gif"
twodfilename = "2Danimation.gif"
zmax = 0.6
maxv = 0.2 #Max Vth value expected (for colors)
zmin = -0.21

function omega(t)
    return 0.1*sin(2*pi*t)#.1*t
end

#Test Function
function f(x)
    if x > 1.5
        return 1
    else
        return 0
    end
end

function plotanulus(rs,vs,time)
    #timetrunc = round(time*100)/100
    #PARAMS
    resolution = 30
    coldv = zmin
    hotv = maxv
    #
    ths = collect(0:2*pi/resolution:2*pi)
    xs = []
    ys = []
    fs=[]
    lims = (-1.5*rs[end],1.5*rs[end])
    top = zmax#maximum(vs) + 1*(maximum(vs)-minimum(vs))
    bot = zmin#minimum(vs) - 1*(maximum(vs)-minimum(vs))
    zrange = (bot,top)
    
    p=plot(xs,ys,fs,xlims = lims,ylims = lims,zlims = zrange,linecolor = RGBA(1,0,0,1),legend = false,camera = (15,50),plot_title="t=$(@sprintf("%.2f", time))"*"s",plot_titlelocation = :left)

    #Plot rings
    for i = 1:1:length(rs)
        value = vs[i]
        if (vs[i] >= 0) && (vs[i] <= hotv)
        elseif vs[i] < 0
            value = 0
        else
            value = hotv
            #println(vs[i])
        end
        colorR = (value-coldv)/(hotv-coldv)
        colorB = 1-colorR
        newxs = @.sin(ths)*rs[i]
        newys = @.cos(ths)*rs[i]
        newfs = @.cos(ths)/cos(ths) *vs[i]
        plot!(newxs,newys,newfs,linecolor = RGBA(colorR,0,colorB,1))
    end

    #plot spokes
    colorS = 15/length(ths)
    for i = 1:1:length(ths)
        newxs = @.sin(ths[i])*rs
        newys = @.cos(ths[i])*rs
        plot!(newxs,newys,vs,linecolor = RGBA(colorS,colorS,colorS,1))
    end
    #Add time:
    #annotate!(0,0,0,"t=$time"*"s")
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

rs = collect(R1:(R2-R1)/(Rresolution-1):R2)
ts = collect(tspan[1]:(tspan[2]-tspan[1])/(nframes-1):tspan[2])
vss = []

anim3d = @animate for frame = 1:1:nframes
    local vs = zeros(length(rs))
    #for i = 1:1:length(vs)
    #    vs[i] = Vth(rs[i],ts[frame])
    #end
    vs = Vth(rs,ts[frame])
    push!(vss,vs)
    plotanulus(rs,vs,ts[frame])
end
anim2d = @animate for frame = 1:1:nframes
    plot(rs,vss[frame],xlims = (0,R2),ylims=(zmin,zmax),legend = false,plot_title="t=$(@sprintf("%.2f", ts[frame]))"*"s",plot_titlelocation = :left)
end


gif(anim3d, threedfilename, fps=FPS)
gif(anim2d, twodfilename, fps=FPS)

