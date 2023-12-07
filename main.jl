using Plots,SpecialFunctions, Roots, ForwardDiff,Printf
#Problem Parameters--------------------------------------------------------------------------------------------------------------------------------
R1 = 1
R2 = 2
nu = .1
rho = 1.0
dPdz = -2.0
function omega(t)
    return 0.2*sin(2*pi*t)#.1*t
end
tspan = (0,3)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Visualization Parameters----------------------------------------------------------------------------------------------------------------------------------------------------------
nframes = 30
number_eigenvalues = 30
integralresolution = 10000
Rresolution = 200
FPS = 10
Vθ_3Dfilename = "Vθ_3Danimation.gif"
Vθ_2Dfilename = "Vθ_2Danimation.gif"
Vz_3Dfilename = "Vz_3Danimation.gif"
Vz_2Dfilename = "Vz_2Danimation.gif"
V_3Dfilename = "V_3Danimation.gif"
nprofiles = 8

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



function progressbar(s,x)
    n = 35
    percentsimulated = round(x*1000)/10
    bar = ""
    for i = 1:1:n
        if x >= (i/n)
            bar = bar*"█"
        else
            bar = bar*"░"
        end
    end
    println("\n\n\n\n\n\n\n\n\n\n\n\n"*s*": "*bar*" $percentsimulated"*"%"*"\n\n")
end

include("Vz_stuff.jl")

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

function plotvprofile(rs,vz,vθ,time,zrange,maxv)

    #PARAMS
    θresolution = 30
    coldv = 0
    hotv = maxv
    az = 15 *pi/180
    #
    xs = []
    ys = []
    fs=[]
    lims = (-1.1*rs[end],1.1*rs[end])
    p=plot(xs,ys,fs,xlims = lims,ylims = lims,zlims = (2*zrange[1],2*zrange[2]),linecolor = RGBA(1,0,0,1),legend = false,camera = (az*180/pi,40),plot_title="t=$(@sprintf("%.2f", time))"*"s",plot_titlelocation = :left)

    #Plot Inner Ring
    θs = collect(0:2*pi/θresolution:2*pi)
    inxs = @.sin(θs)*R1
    inys = @.cos(θs)*R1
    newfs = zeros(length(θs))
    plot!(inxs,inys,newfs,linecolor = RGBA(0,0,0,1))
    #Plot Outer Ring
    outxs = @.sin(θs)*R2
    outys = @.cos(θs)*R2
    plot!(outxs,outys,newfs,linecolor = RGBA(0,0,0,1))
    #Plot spokes
    #plot!([R1,R2],[0,0],[0,0],linecolor = RGBA(.5,.5,.5,1))
    #plot!([-R1,-R2],[0,0],[0,0],linecolor = RGBA(.5,.5,.5,1))
    #plot!([0,0],[R1,R2],[0,0],linecolor = RGBA(.5,.5,.5,1))
    #plot!([0,0],[-R1,-R2],[0,0],linecolor = RGBA(.5,.5,.5,1))
    #Plot Velocity Vectors on Spokes, start with back.
    for i = 1:1:length(rs)
        mag = sqrt(vθ[end-i+1]^2+vz[end-i+1]^2)
        colorR = (mag-coldv)/(hotv-coldv)
        colorB = 1-colorR
        for th = (2*pi)/nprofiles:(2*pi)/nprofiles:(2*pi)
            if (th <= pi+az) && (th >= az)
                xpt = [rs[end-i+1]*cos(th),rs[end-i+1]*cos(th)-vθ[end-i+1]*sin(th)]
                ypt = [rs[end-i+1]*sin(th),rs[end-i+1]*sin(th)+vθ[end-i+1]*cos(th)]
                plot!(xpt,ypt,[0,vz[end-i+1]],linecolor = RGBA(colorR,0,colorB,1))
            end
        end
    end
    #Plot Velocity Vectors on Spokes, end with front.
    for i = 1:1:length(rs)
        mag = sqrt(vθ[i]^2+vz[i]^2)
        colorR = (mag-coldv)/(hotv-coldv)
        colorB = 1-colorR
        for th = (2*pi)/nprofiles:(2*pi)/nprofiles:(2*pi)
            if (th > (pi+az)) || (th < az)
                xpt = [rs[i]*cos(th),rs[i]*cos(th)-vθ[i]*sin(th)]
                ypt = [rs[i]*sin(th),rs[i]*sin(th)+vθ[i]*cos(th)]
                plot!(xpt,ypt,[0,vz[i]],linecolor = RGBA(colorR,0,colorB,1))
            end
        end
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
    progressbar("Finding Eigenvalues",0)
    while length(eigs) < nEigs
        newzeros=find_zeros(char_eq, lo,hi)
        lo = lo + step
        hi = hi + step
        eigs = [eigs;newzeros]
        if newzeros != []
            progressbar("Finding Eigenvalues",length(eigs)/nEigs)
        end
    end
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
        progressbar("Precomputing Vθ",i/nEigs)
    end

    return eigs,yns,ynprimeR2s,sqrnorms
end

eigs,yns,ynprimeR2s,sqrnorms = findeigs(number_eigenvalues)

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

function Vz(r::Vector,t)#TODO: add this function from andrew's code
    return @. Vz2(r,t)
end

function simulate()
    rs = collect(R1:(R2-R1)/(Rresolution-1):R2)
    ts = collect(tspan[1]:(tspan[2]-tspan[1])/(nframes-1):tspan[2])
    vθs = []
    vzs = []
    minvθ = 0
    maxvθ = 0
    minvz = 0
    maxvz = 0
    maxv = 0
    for frame = 1:1:nframes
        #Actual Calculations
        vθ = zeros(length(rs))#TODO: These two lines probably arent needed.
        vz = zeros(length(rs))
        vθ = Vth(rs,ts[frame])
        vz = Vz(rs,ts[frame])
        push!(vθs,vθ)
        push!(vzs,vz)
        #Progress Bar
        progressbar("Simulating",frame/nframes)
        
        #Track Extrema
        if maximum(vθ) > maxvθ
            maxvθ = maximum(vθ)
        end
        if minimum(vθ) < minvθ
            minvθ = minimum(vθ)
        end
        if maximum(vz) > maxvz
            maxvz = maximum(vz)
        end
        if minimum(vz) < minvz
            minvz = minimum(vz)
        end
        if maximum(@. (vθ^2 + vz^2)) > maxv^2
            maxv = maximum(@. sqrt(vθ^2 + vz^2))
        end
    end
    return rs,ts,vθs,minvθ,maxvθ,vzs,minvz,maxvz,maxv
end

rs,ts,vθs,minvθ,maxvθ,vzs,minvz,maxvz,maxv = simulate()
#animate
vθanim3d = @animate for frame = 1:1:nframes
    plotanulus(rs,vθs[frame],ts[frame],(minvθ,maxvθ))
    progressbar("Animating",.2*frame/nframes)
end
vθanim2d = @animate for frame = 1:1:nframes
    plot(rs,vθs[frame],xlims = (0,R2),ylims=(minvθ,maxvθ),legend = false,plot_title="t=$(@sprintf("%.2f", ts[frame]))"*"s",plot_titlelocation = :left)
    progressbar("Animating",.2*frame/nframes+.2)
end
vzanim3d = @animate for frame = 1:1:nframes
    plotanulus(rs,vzs[frame],ts[frame],(minvz,maxvz))
    progressbar("Animating",.2*frame/nframes+.4)
end
vzanim2d = @animate for frame = 1:1:nframes
    plot(rs,vzs[frame],xlims = (0,R2),ylims=(minvz,maxvz),legend = false,plot_title="t=$(@sprintf("%.2f", ts[frame]))"*"s",plot_titlelocation = :left)
    progressbar("Animating",.2*frame/nframes+.6)
end
vanim3d = @animate for frame = 1:1:nframes
    plotvprofile(rs,vzs[frame],vθs[frame],ts[frame],(minvz,maxvz),maxv)
    progressbar("Animating",.2*frame/nframes+.8)
end

dir = "Animations/"
gif(vθanim3d, dir*Vθ_3Dfilename, fps=FPS)
gif(vθanim2d, dir*Vθ_2Dfilename, fps=FPS)
gif(vzanim3d, dir*Vz_3Dfilename, fps=FPS)
gif(vzanim2d, dir*Vz_2Dfilename, fps=FPS)
gif(vanim3d, dir*V_3Dfilename, fps=FPS)