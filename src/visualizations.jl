function store_path(path, ts, us)
    if typeof(path.t[1]) == Float64
        open("files/paths.txt", "a") do io
            DF.writedlm(io, path.t')
            DF.writedlm(io, path)
        end 
        open("files/controls.txt", "a") do io
            DF.writedlm(io, ts')
            DF.writedlm(io, us)
        end 
    else
        tPath = zeros(length(path.t))
        Vinfs = zeros(length(path.t))
        gammas = zeros(length(path.t))
        thetaDots = zeros(length(path.t))
        thetas = zeros(length(path.t))
        posxs = zeros(length(path.t))
        posys = zeros(length(path.t))
        tsU = zeros(length(ts))
        us1 = zeros(length(ts))
        us2 = zeros(length(ts))
        for i in 1:length(path.t)
            tPath[i] = path.t[i].value
            Vinfs[i] = path.u[i][1].value
            gammas[i] = path.u[i][2].value
            thetaDots[i] = path.u[i][3].value
            thetas[i] = path.u[i][4].value
            posxs[i] = path.u[i][5].value
            posys[i] = path.u[i][6].value
        end
        for i in 1:length(ts)
            tsU[i] = ts[i].value
            us1[i] = us[1][i].value
            us2[i] = us[2][i].value
        end
        open("files/paths.txt", "a") do io
            DF.writedlm(io, tPath')
            DF.writedlm(io, Vinfs')
            DF.writedlm(io, gammas')
            DF.writedlm(io, thetaDots')
            DF.writedlm(io, thetas')
            DF.writedlm(io, posxs')
            DF.writedlm(io, posys')
        end 
        open("files/controls.txt", "a") do io
            DF.writedlm(io, tsU')
            DF.writedlm(io, us1')
            DF.writedlm(io, us2')
        end 
    end
end 

function clear_paths()
    open("files/paths.txt", "w") do io
        DF.writedlm(io, "")
    end
    open("files/controls.txt", "w") do io
        DF.writedlm(io, "")
    end
end

function visualize_paths(skip, delay)
    A = readdlm("files/paths.txt")  
    B = readdlm("files/controls.txt")  
    i = 1
    j = 1
    while i < length(A[:,1])
        tPath = truncate(A[i,:])
        Vinfs = truncate(A[i+1,:])
        gammas = truncate(A[i+2,:])
        thetaDots = truncate(A[i+3,:])
        thetas = truncate(A[i+4,:])
        aoas = thetas .- gammas
        posxs = truncate(A[i+5,:])
        posys = truncate(A[i+6,:])
        tU = truncate(B[j,:])
        u1s = truncate(B[j+1,:])
        u2s = truncate(B[j+2,:])
    
        p1 = plot(tPath, Vinfs, xlabel = "Time (s)", ylabel = "State", label = "Velocity", legend = :topleft)
        plot!(tPath, gammas, label = "Flightpath Angle")
        plot!(tPath, thetas, label = "Pitch Angle")
        plot!(tPath, aoas, label = "Angle of Attack")
        plot!(tPath, posys, label = "Y Position")

        # plot!(t, posx, label = "X Position")
        # p2 = plot(t, u_top, xlabel = "Time (s)", ylabel = "Thrust (N)", label = "Top", legend = :topleft)

        p2 = plot(tU, u1s, xlabel = "Time (s)", ylabel = "Input", label = "u1", legend = :topleft)
        plot!(tU, u2s, label = "u2")
        p = plot(p1, p2, layout = 2)
        display(p)

        sleep(delay)

        i += 7*skip
        j += 3*skip
    end
end

function truncate(vector)
    i = 1
    while (i <= length(vector))&&(typeof(vector[i]) == Float64)
        i += 1
    end
    return vector[1:i-1]
end







function trajectory_tracker(model, initial, delay)
    A = readdlm("files/paths.txt")  
    for i in 1:length(A[:,1])
        designVariables = A[i,:]
        us = [designVariables[1:Int((end-1)/2)], designVariables[Int((end-1)/2)+1:end-1]]
        ts = range(0, stop = designVariables[end], length = length(us[1]))
        tSpan = [0, ts[end]]
        path = simulate(initial, ts, us, model, tSpan)
        plot_simulation(path, ts, us)
        sleep(delay)
    end
end

function trim_tracker(model, tSpan, delay)
    A = readdlm("files/iterations.txt")  
    for i in 1:length(A[:,1])
        designVariables = A[i,:]
        u = designVariables[1:2]
        x = designVariables[3:end] 
        ts = range(0, stop = tSpan[2], length = 100)
        us = [u[1]*ones(length(ts)), u[2]*ones(length(ts))]
        path = simulate(x, ts, us, model, tSpan)
        plot_simulation(path, ts, us)
        sleep(delay)
    end
end