
function MM1c(λ, μ, T, c, seed)
    #This function simulates a capacitated M/M/1 queue 
    #Used for Figures 1, 2, and 3
    #λ is the arrival rate
    #μ is the service rate
    #T is the total duration of the simulation
    #c is the capacity for the number in system

    t = 0   #current time in the simulation
    currentnoinsystem = 0   #current number in the system
    NoinSystem = 0   #initialise the number in system series
    EventTimes = 0.0   #initialise the event time series
    Throughput = 0   #initialise the throughput series
    Arrivals = Vector{Float64}()   #initialise empty vector to contain arrival times
    WaitingTimes = 0.0   #initialise waiting times series

    ArrivalEvent = Exponential(1/λ)   
    EitherEvent = Exponential(1/(λ+μ))
    DepartureEvent = Exponential(1/μ)

    Random.seed!(seed)
    while t < T

        if currentnoinsystem == 0    #only an arrival event can happen
            T1 = rand(ArrivalEvent)   #T1 will be the next event time (inter event time)
            currentnoinsystem += 1
            Throughput = vcat(Throughput, Throughput[end])
            Arrivals = vcat(Arrivals, t+T1)
            WaitingTimes = vcat(WaitingTimes, WaitingTimes[end])
        elseif currentnoinsystem == c    #only a departure event can happen
            T1 = rand(DepartureEvent)
            currentnoinsystem -= 1
            Throughput = vcat(Throughput, Throughput[end]+1)
            WaitingTimes = vcat(WaitingTimes, t+T1 - Arrivals[Int(Throughput[end])])
        else   #either an arrival or a departure can happen next
            T1 = rand(EitherEvent)
            p = rand()   #is the event an arrival or a departure?
            if p < (λ/(λ+μ))
                currentnoinsystem += 1
                Throughput = vcat(Throughput, Throughput[end])
                Arrivals = vcat(Arrivals, t+T1)
                WaitingTimes = vcat(WaitingTimes, WaitingTimes[end])
            else
                currentnoinsystem -= 1
                Throughput = vcat(Throughput, Throughput[end]+1)
                WaitingTimes = vcat(WaitingTimes, t+T1 - Arrivals[Int(Throughput[end])])
            end

        end
        
        t += T1   #progress the simulation clock
        EventTimes = vcat(EventTimes, t)
        NoinSystem = vcat(NoinSystem, currentnoinsystem)
    end

    EventTimes[end] = T
    NoinSystem[end] = NoinSystem[end-1]
    Throughput[end] = Throughput[end-1]
    WaitingTimes[end] = WaitingTimes[end-1]
    series = hcat(EventTimes, NoinSystem, Throughput, WaitingTimes)

    return series
end

function Figure1()
    series = MM1c(2,2,20,5,3)  
    noinsystem_plot = plot(series[:,1], series[:,2], linetype=:steppost, lw=2, color=:plum4, xlabel=L"t", ylabel=L"\textrm{number \ in \ system}", legend=:none, tickfontsize=5, xguidefontsize=8, yguidefontsize=8, dpi=600, margin=3mm)
    throughput_plot = plot(series[:,1], series[:,3], linetype=:steppost, lw=2, color=:plum4, xlabel=L"t", ylabel=L"\textrm{throughput}", legend=:none, tickfontsize=5, xguidefontsize=8, yguidefontsize=8, dpi=600, margin=3mm)
    waittime_plot = plot(series[3:end,1], series[3:end,4], linetype=:steppost, lw=2, color=:plum4, xlabel=L"t", ylabel=L"\textrm{latest \ sojourn \ time}", legend=:none, tickfontsize=5, xguidefontsize=8, yguidefontsize=8, dpi=600)
    final_plot = plot(noinsystem_plot, throughput_plot, waittime_plot, layout=(1,3), size=(800,200), margin=3mm)
    return final_plot
end

function distEuclid_shapelet_series(shapelet, series)
    #This function calculates the distance between a shapelet and a series using Euclidean distance in the time series setting. No z-normalisation
    #Used for Figure 2
    ℓ = length(shapelet)
    smallest_dist = sum((shapelet .- series[1:ℓ]).^2)
    position = 1
    for i in 2:(length(series)-ℓ+1)
        test_dist = 0
        j = 1
        while test_dist < smallest_dist && j <= ℓ
            test_dist += (shapelet[j] - series[i+j-1])^2
            j += 1
        end
        if test_dist < smallest_dist
            smallest_dist = test_dist
            position = i
        end
    end
    return smallest_dist
end

function Figure2()
    shapelet = MM1c(1,1,3,5,1)
    twoseries = [MM1c(1,1,20,5,17), MM1c(1,1,20,5,16)]
    m = twoseries[1][end,1]
    ℓ = shapelet[end,1]

    distance_samplings = Vector{Array{Float64}}(undef, 0)
    sampling_range = 0.1:0.01:0.5
    for series in twoseries
        pos = 0
        distance_freq = zeros(length(sampling_range))
        for freq in sampling_range
            seriesRegTimes = 0:freq:m
            seriesRegValues = zeros(length(seriesRegTimes))
            for i in eachindex(seriesRegTimes)
                place = argmax(j -> series[:,1][j], findall(x -> x<=seriesRegTimes[i], series[:,1]))
                seriesRegValues[i] = series[place,2]
            end
            seriesReg = hcat(seriesRegTimes, seriesRegValues)
            shapeletRegTimes = 0:freq:ℓ
            shapeletRegValues = zeros(length(shapeletRegTimes))
            for i in eachindex(shapeletRegTimes)
                place = argmax(j -> shapelet[:,1][j], findall(x -> x<=shapeletRegTimes[i], shapelet[:,1]))
                shapeletRegValues[i] = shapelet[place,2]
            end
            shapeletReg = hcat(shapeletRegTimes, shapeletRegValues)

            pos += 1
            distance_freq[pos] = distEuclid_shapelet_series(shapeletReg[:,2], seriesReg[:,2])
        end
        distances_series = hcat(sampling_range, distance_freq)
        push!(distance_samplings, distances_series)
    end

    shapelet[:,1] .+= 21 #to plot shapelet on the right hand side after the series

    series_plot = plot([twoseries[i][:,1] for i in 1:2], [twoseries[i][:,2] for i in 1:2], xlim=(0,25), ylim=(-0.2,6), xticks=0:5:20, linetype=:steppost, linewidth=2, legend=:topright, label=[L"y_1" L"y_2"], color=[:plum4 :sienna1], xlabel=L"t", ylabel=L"y(t)", tickfontsize=6, xguidefontsize=10, yguidefontsize=10, legendfontsize=8, margins=3mm, right_margin=10mm, size=(800,200), dpi=600)
    series_plot = plot!(shapelet[:,1], shapelet[:,2], linetype=:steppost, linewidth=2, color=:mediumseagreen, label=L"s")
    series_plot = plot!(twinx(),shapelet[:,1], shapelet[:,2], xlim=(0,25), ylim=(-0.2,6), ylabel=L"s(t)", tickfontsize=6, xticks=:none, label=:none, linetype=:steppost, linewidth=1, color=:mediumseagreen)

    sampfreq_plot = scatter([distance_samplings[i][:,1] for i in 1:2], [distance_samplings[i][:,2] for i in 1:2], ms=3, markerstrokewidth=0.5, thickness_scaling=1, color=[:plum4 :sienna1], legend=:topright, label=[L"Z_1" L"Z_2"], xlabel=L"\mathrm{sampling \ frequency}", ylabel=L"\mathrm{dist}(S,Z)", tickfontsize=6, xguidefontsize=10, yguidefontsize=10, legendfontsize=8, margins=3mm, size=(800,200), dpi=600)
    sampfreq_plot = plot!([distance_samplings[i][:,1] for i in 1:2], [distance_samplings[i][:,2] for i in 1:2], ms=3, color=[:plum4 :sienna1], label=:none)

    final_plot = plot(series_plot, sampfreq_plot, layout=(2,1), size=(800,450), margins=3mm)
    return final_plot
end

function Figure3()
    series = MM1c(1,1,9,5,1)
    shapelet = MM1c(1,1,2,5,1)
    shapelet[:,1] .+= 2
    series[:,1] .-= 0.2
    series[1,1] = 0.0
    series[end,1] = 9.0

    Times_combined = unique(sort(vcat(series[:,1], shapelet[:,1])))
    series_combinedvalues = zeros(length(Times_combined))
    shapelet_combinedvalues = zeros(length(Times_combined))
    shapeletindices = Vector{Int64}()
    for i in eachindex(series_combinedvalues)
        place1 = argmax(j -> series[:,1][j], findall(x -> x<=Times_combined[i], series[:,1]))
        series_combinedvalues[i] = series[place1,2]
        shapelet_combinedvalues[i] = series[place1,2]
        if Times_combined[i] >= shapelet[1,1] && Times_combined[i] < shapelet[end,1]
            place2 = argmax(j -> shapelet[:,1][j], findall(x -> x<=Times_combined[i], shapelet[:,1]))
            shapelet_combinedvalues[i] = shapelet[place2,2]
            shapeletindices = vcat(shapeletindices, i)
        end
    end
    series_combined = hcat(Times_combined, series_combinedvalues)
    shapelet_combined = hcat(Times_combined, shapelet_combinedvalues)
    only_shapelet_combined = shapelet_combined[shapeletindices,:]
    only_shapelet_combined = [only_shapelet_combined; shapelet[end,1:2]']

    series_shapelet_plot = plot(Times_combined, series_combined[:,2], fillalpha=0.2, linetype=:steppost, color=:plum4, linewidth=2, xlabel=L"t", ylabel=L"y(t)", label=:none, xguidefontsize=10, tickfontsize=6, yguidefontsize=10, legendfontsize=8, margins=2mm, ylim=(-0.2,6), size=(800,200), dpi=600)
    series_shapelet_plot = plot!(Times_combined, series_combined[:,2], linetype=:steppost, color=:plum4, linewidth=2, label=L"y")
    series_shapelet_plot = plot!(only_shapelet_combined[:,1], only_shapelet_combined[:,2], linetype=:steppost, color=:mediumseagreen, label=L"s", linewidth=2)
    series_shapelet_plot = plot!(Times_combined, series_combined[:,2], fillrange=shapelet_combined[:,2], fillcolor=:mediumseagreen, fillalpha=0.2,  linetype=:steppost, color=:mediumseagreen, alpha=0.2, linewidth=0.1, label=L"\Vert s - y([2,4]) \Vert_1")

    seriesRegTimes = 0:0.2:9
    seriesRegValues = zeros(length(seriesRegTimes))
    shapeletindices = Vector{Int64}()
    for i in eachindex(seriesRegTimes)
        place = argmax(j -> series[:,1][j], findall(x -> x<=seriesRegTimes[i], series[:,1]))
        seriesRegValues[i] = series[place,2]
        if seriesRegTimes[i] >= shapelet[1,1] && seriesRegTimes[i] <=shapelet[end,1]
            shapeletindices = vcat(shapeletindices, i)
        end
    end

    shapeletRegTimes = shapelet[1,1]:0.2:shapelet[end,1]
    shapeletRegValues = zeros(length(shapeletRegTimes))
    for i in eachindex(shapeletRegTimes)
        place = argmax(j -> shapelet[:,1][j], findall(x -> x<=shapeletRegTimes[i], shapelet[:,1]))
        shapeletRegValues[i] = shapelet[place,2]
    end

    regsampling_plot = scatter(seriesRegTimes, seriesRegValues, markerstrokewidth=0.5, ms=3, color=:plum4, xlabel=L"i", ylabel=L"z_i", label=L"Z", xticks=(0:2:8, ["1","11","21","31","41"]), xguidefontsize=10, tickfontsize=6, yguidefontsize=10, legendfontsize=8, margins=2mm, ylim=(-0.2,6), size=(800,200), dpi=600)
    regsampling_plot = scatter!(shapeletRegTimes, shapeletRegValues, color=:mediumseagreen, markerstrokewidth=0.5, ms=3, label=L"S")
    regsampling_plot = plot!([[shapeletRegTimes[i], shapeletRegTimes[i]] for i in eachindex(shapeletRegTimes)], [[shapeletRegValues[i], seriesRegValues[shapeletindices][i]] for i in eachindex(shapeletRegTimes)], lw=2, color=:mediumseagreen, alpha=0.5, label=:none, linestyle=:dash)
    regsampling_plot = plot!([shapeletRegTimes[1], shapeletRegTimes[1]], [shapeletRegValues[1], shapeletRegValues[1]], lw=2, color=:mediumseagreen, alpha=0.5, label=L"d(S,z_{11:21})", linestyle=:dash)

    final_plot = plot(regsampling_plot, series_shapelet_plot, layout=(2,1), link=:x, size=(800,450))
    return final_plot
end

function Figure4()
    LRWcontrolled_s₁ = load_object("res/Breakdowns/LRWcontrolled_s1.jld2")
    MRWcontrolled_s₂ = load_object("res/Breakdowns/MRWcontrolled_s2.jld2")

    traindists_s₁_LRW = load_object("res/Breakdowns/controlled_traindists_s1_LRW.jld2") 
    traindists_s₁_MRW = load_object("res/Breakdowns/controlled_traindists_s1_MRW.jld2") 
    traindists_s₂_LRW = load_object("res/Breakdowns/controlled_traindists_s2_LRW.jld2") 
    traindists_s₂_MRW = load_object("res/Breakdowns/controlled_traindists_s2_MRW.jld2")
    testdists_s₁_LRW = load_object("res/Breakdowns/controlled_testdists_s1_LRW.jld2")
    testdists_s₁_MRW = load_object("res/Breakdowns/controlled_testdists_s1_MRW.jld2")
    testdists_s₂_LRW = load_object("res/Breakdowns/controlled_testdists_s2_LRW.jld2")
    testdists_s₂_MRW = load_object("res/Breakdowns/controlled_testdists_s2_MRW.jld2")

    WFshapelets = plot(LRWcontrolled_s₁[1][:,1], LRWcontrolled_s₁[1][:,2] .- minimum(LRWcontrolled_s₁[1][:,2]), linetype=:steppost, lw=2, color=:sienna1, xlabel=L"t", ylabel=L"s(t) - \min_u s(u)", label=L"s_1", tickfontsize=6, xguidefontsize=9, yguidefontsize=9, legendfontsize=9, legend=:topleft, dpi=600)
    WFshapelets = plot!(MRWcontrolled_s₂[1][:,1], MRWcontrolled_s₂[1][:,2] .- minimum(MRWcontrolled_s₂[1][:,2]), linetype=:steppost, lw=2, color=:plum4, label=L"s_2")
    WFscatter_train = scatter(traindists_s₁_LRW, traindists_s₂_LRW, color=:sienna1, alpha=0.5, ms=2, markerstrokewidth=0.1, label="LRW", xlabel=L"\tilde{\mathrm{dist}}(s_1, y)", ylabel=L"\tilde{\mathrm{dist}}(s_2, y)", xlim=(-4,130),ylim=(-4,130), yticks=0:40:120, xticks=0:40:120, legend=:topright, tickfontsize=6, xguidefontsize=9, yguidefontsize=9, legendfontsize=7, titlefontsize=9, dpi=600) #xlim=(0,15000), ylim=(0,8200), xticks=(0:5000:15000, 0:7), yticks=(0:1000:6000, 0:6)
    WFscatter_train = annotate!(65, 135, text(L"\textbf{training}", :black, 9))
    WFscatter_train = scatter!(traindists_s₁_MRW, traindists_s₂_MRW, color=:plum4, alpha=0.5, ms=2, markerstrokewidth=0.1, label="MRW")
    WFscatter_test = scatter(testdists_s₁_LRW, testdists_s₂_LRW, color=:sienna1, alpha=0.5, ms=2, markerstrokewidth=0.1, label="LRW", xlabel=L"\tilde{\mathrm{dist}}(s_1, y)", ylabel=L"\tilde{\mathrm{dist}}(s_2, y)", xlim=(-4,130),ylim=(-4,130), yticks=0:40:120, xticks=0:40:120, legend=:topright, tickfontsize=6, xguidefontsize=9, yguidefontsize=9, legendfontsize=7, titlefontsize=9, dpi=600)
    WFscatter_test = annotate!(65, 135, text(L"\textbf{testing}", :black, 9))
    WFscatter_test = scatter!(testdists_s₁_MRW, testdists_s₂_MRW, color=:plum4, alpha=0.5, ms=2, markerstrokewidth=0.1, label="MRW")
    l = @layout [a{0.4w} b{0.3w} c{0.3w}]
    final_plot = plot(WFshapelets, WFscatter_train, WFscatter_test, layout=l, size=(800,200), margin=3mm)
    return final_plot
end

function Figure5()
    LRWrandom_s₁ = load_object("res/Breakdowns/LRWrandom_s1.jld2")
    MRWrandom_s₂ = load_object("res/Breakdowns/MRWrandom_s2.jld2")

    r_traindists_s₁_LRW = load_object("res/Breakdowns/random_traindists_s1_LRW.jld2") 
    r_traindists_s₁_MRW = load_object("res/Breakdowns/random_traindists_s1_MRW.jld2") 
    r_traindists_s₂_LRW = load_object("res/Breakdowns/random_traindists_s2_LRW.jld2") 
    r_traindists_s₂_MRW = load_object("res/Breakdowns/random_traindists_s2_MRW.jld2")
    r_testdists_s₁_LRW = load_object("res/Breakdowns/random_testdists_s1_LRW.jld2")
    r_testdists_s₁_MRW = load_object("res/Breakdowns/random_testdists_s1_MRW.jld2")
    r_testdists_s₂_LRW = load_object("res/Breakdowns/random_testdists_s2_LRW.jld2")
    r_testdists_s₂_MRW = load_object("res/Breakdowns/random_testdists_s2_MRW.jld2")

    WFshapelets = plot(LRWrandom_s₁[1][:,1], LRWrandom_s₁[1][:,2] .- minimum(LRWrandom_s₁[1][:,2]), linetype=:steppost, lw=2, color=:sienna1, xlabel=L"t", ylabel=L"s(t) - \min_u s(u)", label=L"s_1", xticks=0:10:80, tickfontsize=6, xguidefontsize=9, yguidefontsize=9, legendfontsize=9, legend=:topleft, dpi=600)
    WFshapelets = plot!(MRWrandom_s₂[1][:,1], MRWrandom_s₂[1][:,2] .- minimum(MRWrandom_s₂[1][:,2]), linetype=:steppost, lw=2, color=:plum4, label=L"s_2")
    WFscatter_train = scatter(r_traindists_s₁_LRW, r_traindists_s₂_LRW, color=:sienna1, alpha=0.5, ms=2, markerstrokewidth=0.1, label="LRW", xlabel=L"\tilde{\mathrm{dist}}(s_1, y)", ylabel=L"\tilde{\mathrm{dist}}(s_2, y)", xlim=(-8,400), ylim=(-20, 840), legend=:topleft, tickfontsize=6, xguidefontsize=9, yguidefontsize=9, legendfontsize=7, titlefontsize=9, dpi=600)
    WFscatter_train = annotate!(200, 850, text(L"\textbf{training}", :black, 9))
    WFscatter_train = scatter!(r_traindists_s₁_MRW, r_traindists_s₂_MRW, color=:plum4, alpha=0.5, ms=2, markerstrokewidth=0.1, label="MRW")
    WFscatter_test =scatter(r_testdists_s₁_LRW, r_testdists_s₂_LRW, color=:sienna1, alpha=0.5, ms=2, markerstrokewidth=0.1, label="LRW", xlabel=L"\tilde{\mathrm{dist}}(s_1, y)", ylabel=L"\tilde{\mathrm{dist}}(s_2, y)", xlim=(-8,400), ylim=(-20, 840), legend=:topleft, tickfontsize=6, xguidefontsize=9, yguidefontsize=9, legendfontsize=7, titlefontsize=9, dpi=600)
    WFscatter_test =annotate!(200, 850, text(L"\textbf{testing}", :black, 9))
    WFscatter_test = scatter!(r_testdists_s₁_MRW, r_testdists_s₂_MRW, color=:plum4, alpha=0.5, ms=2, markerstrokewidth=0.1, label="MRW")
    l = @layout [a{0.4w} b{0.3w} c{0.3w}]
    final_plot = plot(WFshapelets, WFscatter_train, WFscatter_test, layout=l, size=(800,200), margin=3mm)
    return final_plot
end


function Figure6()
    t = 0.01:0.01:10
    arrivalrates = plot(t, (sin.(2*pi*t/10)./2).+0.5, color=:mediumseagreen, linewidth=2, xlabel=L"t", ylabel=L"\lambda(t)", label=L"\textrm{true}", yticks=0:0.5:1, xticks=0:5:10, dpi=600, size=(300,300),legend=:topright, tickfontsize=6, xguidefontsize=8, yguidefontsize=8, legendfontsize=8, margins=2mm)
    arrivalrates = plot!([0,5],[1,1], color=:plum4, linewidth=2,label="")
    arrivalrates = plot!([5,10],[0,0], color=:plum4, linewidth=2,label=L"\textrm{model \ 1}")
    arrivalrates = plot!([0,(10/(2*pi))*asin(0.5)],[0.5,0.5], color=:sienna1, linestyle=:dash, linewidth=2,label="")
    arrivalrates = plot!([(10/(2*pi))*asin(0.5), 5-(10/(2*pi))*asin(0.5)],[1,1], color=:sienna1, linestyle=:dash, linewidth=2,label="")
    arrivalrates = plot!([5-(10/(2*pi))*asin(0.5), 5+(10/(2*pi))*asin(0.5)],[0.5,0.5], color=:sienna1, linestyle=:dash, linewidth=2,label="")
    arrivalrates = plot!([5+(10/(2*pi))*asin(0.5), 10-(10/(2*pi))*asin(0.5)],[0,0], color=:sienna1, linestyle=:dash, linewidth=2,label="")
    arrivalrates = plot!([10-(10/(2*pi))*asin(0.5), 10],[0.5,0.5], color=:sienna1, linestyle=:dash, linewidth=2,label=L"\textrm{model \ 2}")

    GR.setarrowsize(0.5)
    rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
    Queuediagram = plot(xlim=(0,4), ylim=(0,14), grid=:false, xaxis=:false,yaxis=:false,xticks=:none,yticks=:none, size=(100,300), dpi=600)
    Queuediagram = plot!(rectangle(4,2,0,2), lw=0.5, opacity=0.7, color=:mediumseagreen, label=:none)
    Queuediagram = plot!(rectangle(4,2,0,6), lw=0.5, opacity=0.7, color=:mediumseagreen, label=:none)
    Queuediagram = plot!(rectangle(4,2,0,10), lw=0.5, opacity=0.7, color=:mediumseagreen, label=:none)
    Queuediagram = annotate!(2, 11, text(L"\textrm{Station} \ A", :black, 8))
    Queuediagram = annotate!(2, 7, text(L"\textrm{Station} \ B", :black, 8))
    Queuediagram = annotate!(2, 3, text(L"\textrm{Station} \ C", :black, 8))
    Queuediagram = plot!([2,2],[14,12.05],arrow=true,color=:black,linewidth=1,label="")
    Queuediagram = plot!([2,2],[10,8.05],arrow=true,color=:black,linewidth=1,label="")
    Queuediagram = plot!([2,2],[6,4.05],arrow=true,color=:black,linewidth=1,label="")
    Queuediagram = plot!([2,2],[2,0.05],arrow=true,color=:black,linewidth=1,label="")
    Queuediagram = annotate!(1.25, 13, text(L"\lambda(t)", :black, 8))
    Queuediagram = annotate!(2.75, 1, text(L"\textrm{exit}", :black, 8))

    l = @layout [a{0.2w} b{0.8w}]
    final_plot = plot(Queuediagram, arrivalrates, layout = l, size=(600,250), dpi=600)
    return final_plot
end


function Figure7()
    Model1_s₁ = load_object("res/Validation/Model1_s1.jld2")
    Model1_s₂ = load_object("res/Validation/Model1_s2.jld2")

    m1_traindists_s₁_true = load_object("res/Validation/m1_traindists_s1_true.jld2")
    m1_traindists_s₁_model = load_object("res/Validation/m1_traindists_s1_model.jld2")
    m1_traindists_s₂_true = load_object("res/Validation/m1_traindists_s2_true.jld2")
    m1_traindists_s₂_model = load_object("res/Validation/m1_traindists_s2_model.jld2")
    m1_testdists_s₁_true = load_object("res/Validation/m1_testdists_s1_true.jld2")
    m1_testdists_s₁_model = load_object("res/Validation/m1_testdists_s1_model.jld2")
    m1_testdists_s₂_true = load_object("res/Validation/m1_testdists_s2_true.jld2")
    m1_testdists_s₂_model = load_object("res/Validation/m1_testdists_s2_model.jld2")


    m1_shapelets = plot(Model1_s₁[1][:,1], Model1_s₁[1][:,2] .- minimum(Model1_s₁[1][:,2]), linetype=:steppost, linewidth=2, xlabel=L"t", ylabel=L"s(t) - \min_u s(u)", color=:mediumseagreen, label=L"s_1", xticks=0:2:10, xguidefontsize=8, yguidefontsize=8, legendfontsize=8, xtickfontsize=5, ytickfontsize=5, legend=:topleft, dpi=600)
    m1_shapelets = plot!(Model1_s₂[1][:,1], Model1_s₂[1][:,2] .- minimum(Model1_s₂[1][:,2]), linetype=:steppost, linewidth=2, color=:plum4, label=L"s_2")
    m1_scatter_train = scatter(m1_traindists_s₁_true, m1_traindists_s₂_true, color=:mediumseagreen, alpha=0.5, ms=2, markerstrokewidth=0.1, label=L"\textrm{true}", xlim=(-0.1,5.8), ylim=(-0.35,10), yticks=0:2:10, xlabel=L"\tilde{\mathrm{dist}}(s_1, y)", ylabel=L"\tilde{\mathrm{dist}}(s_2, y)", legend=:bottomleft, xguidefontsize=8, yguidefontsize=8, legendfontsize=4, xtickfontsize=5, ytickfontsize=5, titlefontsize=8, dpi=600)
    m1_scatter_train = annotate!(2.85, 10.5, text(L"\textbf{training}", :black, 8))
    m1_scatter_train = scatter!(m1_traindists_s₁_model, m1_traindists_s₂_model, color=:plum4, alpha=0.5, ms=2, markerstrokewidth=0.1, label=L"\textrm{model \ 1}")
    m1_scatter_test = scatter(m1_testdists_s₁_true, m1_testdists_s₂_true, color=:mediumseagreen, alpha=0.5, ms=2, markerstrokewidth=0.1, label=L"\textrm{true}", xlim=(0.75,7.8), ylim=(1.5,10), xlabel=L"\tilde{\mathrm{dist}}(s_1, y)", ylabel=L"\tilde{\mathrm{dist}}(s_2, y)", legend=:bottomright, xguidefontsize=8, yguidefontsize=8, legendfontsize=4, xtickfontsize=5, ytickfontsize=5, titlefontsize=8, dpi=600)
    m1_scatter_test = annotate!(4.275, 10.5, text(L"\textbf{testing}", :black, 8))
    m1_scatter_test = scatter!(m1_testdists_s₁_model, m1_testdists_s₂_model, color=:plum4, alpha=0.5, ms=2, markerstrokewidth=0.1, label=L"\textrm{model \ 2}")

    true_dists = Matrix{Float64}(hcat(m1_testdists_s₁_true, m1_testdists_s₂_true)')
    trueGaussian = fit(MvNormal, true_dists)
    model1_dists = Matrix{Float64}(hcat(m1_testdists_s₁_model, m1_testdists_s₂_model)')
    model1Gaussian = fit(MvNormal, model1_dists)
    X = range(1.5, 5.5, length=100)
    Y = range(3, 9, length=100)
    ftrue1(x,y) = pdf(trueGaussian, [x,y])
    fmodel1(x,y) = pdf(model1Gaussian, [x,y])
    contour(X, Y, ftrue1, linewidth=2, levels=10, color=cgrad([:white, :mediumseagreen], [0.3, 0.7]), cbar=false, xlabel=L"\tilde{\mathrm{dist}}(s_1, y)", ylabel=L"\tilde{\mathrm{dist}}(s_2, y)", yticks=4:2:8, xguidefontsize=8, yguidefontsize=8, legendfontsize=6, xtickfontsize=5, ytickfontsize=5, dpi=600)
    contour!(X, Y, fmodel1, linewidth=2, levels=10, color=cgrad([:white, :plum4], [0.3, 0.7]))
    Model1contours = plot!([3,3.01],[6.95,6.96], linewidth=0.5, color=:mediumseagreen, label=L"\textrm{true}", legendfontsize=4)
    Model1contours = plot!([3.95,3.95],[5.75,5.76], linewidth=0.5, color=:plum4, label=L"\textrm{model \ 1}", legendfontsize=4)

    final_plot = plot(m1_shapelets, m1_scatter_train, m1_scatter_test, Model1contours, layout=(2,2), margins=1mm, size=(600,300), dpi=600)
    return final_plot
end

function Figure8()
    Model2_s₁ = load_object("res/Validation/Model2_s1.jld2")
    Model2_s₂ = load_object("res/Validation/Model2_s2.jld2")

    m2_traindists_s₁_true = load_object("res/Validation/m2_traindists_s1_true.jld2")
    m2_traindists_s₁_model = load_object("res/Validation/m2_traindists_s1_model.jld2")
    m2_traindists_s₂_true = load_object("res/Validation/m2_traindists_s2_true.jld2")
    m2_traindists_s₂_model = load_object("res/Validation/m2_traindists_s2_model.jld2")
    m2_testdists_s₁_true = load_object("res/Validation/m2_testdists_s1_true.jld2")
    m2_testdists_s₁_model = load_object("res/Validation/m2_testdists_s1_model.jld2")
    m2_testdists_s₂_true = load_object("res/Validation/m2_testdists_s2_true.jld2")
    m2_testdists_s₂_model = load_object("res/Validation/m2_testdists_s2_model.jld2")

    m2_shapelets = plot(Model2_s₁[1][:,1], Model2_s₁[1][:,2] .- minimum(Model2_s₁[1][:,2]), linewidth=2, linetype=:steppost, xlabel=L"t", ylabel=L"s(t) - \min_u s(u)", color=:mediumseagreen, label=L"s_1", xticks=0:2:12, xlim=(0,12), xguidefontsize=8, yguidefontsize=8, legendfontsize=8, xtickfontsize=5, ytickfontsize=5, legend=:topleft, dpi=600)
    m2_shapelets = plot!(Model2_s₂[1][:,1], Model2_s₂[1][:,2] .- minimum(Model2_s₂[1][:,2]), linewidth=2, linetype=:steppost, color=:sienna1, label=L"s_2")
    m2_scatter_train = scatter(m2_traindists_s₁_true, m2_traindists_s₂_true, color=:mediumseagreen, alpha=0.5, ms=2, markerstrokewidth=0.1, label=L"\textrm{true}", xlim=(-0.2,12.5), ylim=(-0.3,9), xlabel=L"\tilde{\mathrm{dist}}(s_1, y)", ylabel=L"\tilde{\mathrm{dist}}(s_2, y)", xguidefontsize=8, yguidefontsize=8, legendfontsize=4, xtickfontsize=5, ytickfontsize=5, titlefontsize=8, legend=:bottomleft, dpi=600)
    m2_scatter_train = annotate!(6.15, 9.5, text(L"\textbf{training}", :black, 8))
    m2_scatter_train = scatter!(m2_traindists_s₁_model, m2_traindists_s₂_model, color=:sienna1, alpha=0.5, ms=2, markerstrokewidth=0.1, label=L"\textrm{model \ 2}")
    m2_scatter_test = scatter(m2_testdists_s₁_true, m2_testdists_s₂_true, color=:mediumseagreen, alpha=0.5, ms=2, markerstrokewidth=0.1, label=L"\textrm{true}", xlim=(2.6,15.1), ylim=(1.2,10), xticks=4:2:14, xlabel=L"\tilde{\mathrm{dist}}(s_1, y)", ylabel=L"\tilde{\mathrm{dist}}(s_2, y)", xguidefontsize=8, yguidefontsize=8, legendfontsize=4, xtickfontsize=5, ytickfontsize=5, titlefontsize=8, legend=:bottomright, dpi=600)
    m2_scatter_test = annotate!(8.85, 10.5, text(L"\textbf{testing}", :black, 8))
    m2_scatter_test = scatter!(m2_testdists_s₁_model, m2_testdists_s₂_model, color=:sienna1, alpha=0.5, ms=2, markerstrokewidth=0.1, label=L"\textrm{model \ 2}")

    true_dists = Matrix{Float64}(hcat(m2_testdists_s₁_true, m2_testdists_s₂_true)')
    trueGaussian = fit(MvNormal, true_dists)
    model2_dists = Matrix{Float64}(hcat(m2_testdists_s₁_model, m2_testdists_s₂_model)')
    model2Gaussian = fit(MvNormal, model2_dists)
    X = range(4.5, 10.5, length=100)
    Y = range(2.5, 8, length=100)
    ftrue2(x,y) = pdf(trueGaussian, [x,y])
    fmodel2(x,y) = pdf(model2Gaussian, [x,y])
    contour(X, Y, ftrue2, linewidth=2, levels=10, color=cgrad([:white, :mediumseagreen], [0.3, 0.7]), cbar=false, xlabel=L"\tilde{\mathrm{dist}}(s_1, y)", ylabel=L"\tilde{\mathrm{dist}}(s_2, y)", xticks=6:2:10, yticks=4:2:8, xguidefontsize=8, yguidefontsize=8, xtickfontsize=5, ytickfontsize=5, legendfontsize=6, dpi=600)
    contour!(X, Y, fmodel2, linewidth=2, levels=10, color=cgrad([:white, :sienna1], [0.3, 0.7]))
    Model2contours = plot!([7,7],[5.2,5.2], linewidth=0.5, color=:mediumseagreen, label=L"\textrm{true}", legendfontsize=4)
    Model2contours = plot!([7,7],[5.55,5.55], linewidth=0.5,color=:sienna1, label=L"\textrm{model \ 2}", legendfontsize=4)

    final_plot = plot(m2_shapelets, m2_scatter_train, m2_scatter_test, Model2contours, layout=(2,2), margins=1mm, size=(600,300), dpi=600)
    return final_plot
end

function Figure9()
    IndepA_s₁ = load_object("res/Multivariate/IndepA_s1.jld2")
    IndepA_s₂ = load_object("res/Multivariate/IndepA_s2.jld2")
    IndepB_s₁ = load_object("res/Multivariate/IndepB_s1.jld2")
    IndepB_s₂ = load_object("res/Multivariate/IndepB_s2.jld2")
    IndepC_s₁ = load_object("res/Multivariate/IndepC_s1.jld2")
    IndepC_s₂ = load_object("res/Multivariate/IndepC_s2.jld2")

    testdistsA_s₁_system1 = load_object("res/Multivariate/testdistsA_s1_system1.jld2")
    testdistsA_s₁_system2 = load_object("res/Multivariate/testdistsA_s1_system2.jld2")
    testdistsA_s₂_system1 = load_object("res/Multivariate/testdistsA_s2_system1.jld2")
    testdistsA_s₂_system2 = load_object("res/Multivariate/testdistsA_s2_system2.jld2")
    testdistsB_s₁_system1 = load_object("res/Multivariate/testdistsB_s1_system1.jld2")
    testdistsB_s₁_system2 = load_object("res/Multivariate/testdistsB_s1_system2.jld2")
    testdistsB_s₂_system1 = load_object("res/Multivariate/testdistsB_s2_system1.jld2")
    testdistsB_s₂_system2 = load_object("res/Multivariate/testdistsB_s2_system2.jld2")
    testdistsC_s₁_system1 = load_object("res/Multivariate/testdistsC_s1_system1.jld2")
    testdistsC_s₁_system2 = load_object("res/Multivariate/testdistsC_s1_system2.jld2")
    testdistsC_s₂_system1 = load_object("res/Multivariate/testdistsC_s2_system1.jld2")
    testdistsC_s₂_system2 = load_object("res/Multivariate/testdistsC_s2_system2.jld2")


    IndepA_shapelets = plot(IndepA_s₁[1][:,1], IndepA_s₁[1][:,2].-minimum(IndepA_s₁[1][:,2]), linetype=:steppost, lw=2, color=:sienna1, xlabel=L"t", ylabel=L"s(t) - \min_u s(u)", label=L"s_{A1}", ylim=(-0.3,6.3), xlim=(-0.3, 12), yticks=0:2:6, xticks=0:2:12, tickfontsize=5, xguidefontsize=8, yguidefontsize=8, legendfontsize=6, legend=:topright, size=(600,500), dpi=600)
    IndepA_shapelets = plot!(IndepA_s₂[1][:,1], IndepA_s₂[1][:,2].-minimum(IndepA_s₂[1][:,2]), linetype=:steppost, lw=2, color=:plum4, label=L"s_{A2}")
    IndepA_testscatter = scatter(testdistsA_s₁_system1, testdistsA_s₂_system1, color=:sienna1, alpha=0.5, ms=2, markerstrokewidth=0.1, label=L"\textrm{system \ 1}", xlabel=L"\tilde{\mathrm{dist}}(s_{A1}, y_A)", ylabel=L"\tilde{\mathrm{dist}}(s_{A2}, y_A)", legend=:topleft, tickfontsize=5, xguidefontsize=8, yguidefontsize=8, legendfontsize=4, titlefontsize=8, size=(600,500), dpi=600)
    IndepA_testscatter = scatter!(testdistsA_s₁_system2, testdistsA_s₂_system2, color=:plum4, alpha=0.5, ms=2, markerstrokewidth=0.1, label=L"\textrm{system \ 2}")

    IndepB_shapelets = plot(IndepB_s₁[1][:,1], IndepB_s₁[1][:,2].-minimum(IndepB_s₁[1][:,2]), linetype=:steppost, lw=2, color=:sienna1, xlabel=L"t", ylabel=L"s(t) - \min_u s(u)", label=L"s_{B1}", ylim=(-0.3,6.3), xlim=(-0.3, 12), yticks=0:2:6, xticks=0:2:12, tickfontsize=5, xguidefontsize=8, yguidefontsize=8, legendfontsize=6, legend=:topright, size=(600,500), dpi=600)
    IndepB_shapelets = plot!(IndepB_s₂[1][:,1], IndepB_s₂[1][:,2].-minimum(IndepB_s₂[1][:,2]), linetype=:steppost, lw=2, color=:plum4, label=L"s_{B2}")
    IndepB_testscatter = scatter(testdistsB_s₁_system1, testdistsB_s₂_system1, color=:sienna1, alpha=0.5, ms=2, markerstrokewidth=0.1, label=L"\textrm{system \ 1}", xlabel=L"\tilde{\mathrm{dist}}(s_{B1}, y_B)", ylabel=L"\tilde{\mathrm{dist}}(s_{B2}, y_B)", legend=:topright, tickfontsize=5, xguidefontsize=8, yguidefontsize=8, legendfontsize=4, titlefontsize=8, size=(600,500), dpi=600)
    IndepB_testscatter = scatter!(testdistsB_s₁_system2, testdistsB_s₂_system2, color=:plum4, alpha=0.5, ms=2, markerstrokewidth=0.1, label=L"\textrm{system \ 2}")

    IndepC_shapelets = plot(IndepC_s₁[1][:,1], IndepC_s₁[1][:,2].-minimum(IndepC_s₁[1][:,2]), linetype=:steppost, lw=2, color=:sienna1, xlabel=L"t", ylabel=L"s(t) - \min_u s(u)", label=L"s_{C1}", ylim=(-0.3,6.3), xlim=(-0.3, 12), yticks=0:2:6, xticks=0:2:12, tickfontsize=5, xguidefontsize=8, yguidefontsize=8, legendfontsize=6, legend=:topleft, size=(600,500), dpi=600)
    IndepC_shapelets = plot!(IndepC_s₂[1][:,1], IndepC_s₂[1][:,2].-minimum(IndepC_s₂[1][:,2]), linetype=:steppost, lw=2, color=:plum4, label=L"s_{C2}")
    IndepC_testscatter = scatter(testdistsC_s₁_system1, testdistsC_s₂_system1, color=:sienna1, alpha=0.5, ms=2, markerstrokewidth=0.1, label=L"\textrm{system \ 1}", xlabel=L"\tilde{\mathrm{dist}}(s_{C1}, y_C)", ylabel=L"\tilde{\mathrm{dist}}(s_{C2}, y_C)", legend=:topright, tickfontsize=5, xguidefontsize=8, yguidefontsize=8, legendfontsize=4, titlefontsize=8, size=(600,500), dpi=600)
    IndepC_testscatter = scatter!(testdistsC_s₁_system2, testdistsC_s₂_system2, color=:plum4, alpha=0.5, ms=2, markerstrokewidth=0.1, label=L"\textrm{system \ 2}")

    final_plot = plot(IndepA_shapelets, IndepB_shapelets, IndepC_shapelets, IndepA_testscatter, IndepB_testscatter, IndepC_testscatter, layout=(2,3), size=(600,300), margin=2mm)
    return final_plot
end

function Figure10()
    Multi_s₁ = load_object("res/Multivariate/Multi_s1.jld2")
    Multi_s₂ = load_object("res/Multivariate/Multi_s2.jld2")

    testdistsmulti_s₁_system1 = load_object("res/Multivariate/testdistsmulti_s1_system1.jld2")
    testdistsmulti_s₁_system2 = load_object("res/Multivariate/testdistsmulti_s1_system2.jld2")
    testdistsmulti_s₂_system1 = load_object("res/Multivariate/testdistsmulti_s2_system1.jld2")
    testdistsmulti_s₂_system2 = load_object("res/Multivariate/testdistsmulti_s2_system2.jld2")

    S1 = plot(Multi_s₁[1][1][:,1], Multi_s₁[1][1][:,2].-minimum(Multi_s₁[1][1][:,2]), linetype=:steppost, lw=1.5, xlabel=L"t", ylim=(-0.3,6.2), xlim=(-0.3, 12), xticks=0:2:12, yticks=0:2:6, ylabel=L"s(t) - \min_u s(u)", color=:sienna1, label=L"s_A", xguidefontsize=8, yguidefontsize=8, tickfontsize=6, legendfontsize=6, legend=:right, dpi=600)
    S1 = plot!(Multi_s₁[1][2][:,1], Multi_s₁[1][2][:,2].-minimum(Multi_s₁[1][2][:,2]), linetype=:steppost, linestyle=:dash, lw=1.5, color=:sienna1, label=L"s_B")
    S1 = plot!(Multi_s₁[1][3][:,1], Multi_s₁[1][3][:,2].-minimum(Multi_s₁[1][3][:,2]), linetype=:steppost, linestyle=:dot, lw=2, color=:sienna1, label=L"s_C")
    S1 = annotate!(6.0, 6.5, text(L"\mathbf{\textit{s}}_1", :black, 12))
    S2 = plot(Multi_s₂[1][1][:,1], Multi_s₂[1][1][:,2].-minimum(Multi_s₂[1][1][:,2]), linetype=:steppost, lw=1.5, xlabel=L"t", ylim=(-0.3,6.2), xlim=(-0.3, 12), xticks=0:2:12, yticks=0:2:6, ylabel=L"s(t) - \min_u s(u)", color=:plum4, label=L"s_A", xguidefontsize=8, yguidefontsize=8, tickfontsize=6, legendfontsize=6, legend=:topright, dpi=600)
    S2 = plot!(Multi_s₂[1][2][:,1], Multi_s₂[1][2][:,2].-minimum(Multi_s₂[1][2][:,2]), linetype=:steppost, linestyle=:dash, lw=1.5, color=:plum4, label=L"s_B")
    S2 = plot!(Multi_s₂[1][3][:,1], Multi_s₂[1][3][:,2].-minimum(Multi_s₂[1][3][:,2]), linetype=:steppost, linestyle=:dot, lw=2, color=:plum4, label=L"s_C")
    S2 = annotate!(6.0, 6.5, text(L"\mathbf{\textit{s}}_2", :black, 12))
    scatter_test = scatter(testdistsmulti_s₁_system1, testdistsmulti_s₂_system1, color=:sienna1, alpha=0.5, ms=2, markerstrokewidth=0.1, label=L"\textrm{system \ 1}", xlabel=L"\tilde{\mathrm{dist}}(\mathbf{\textit{s}}_1,\mathbf{\textit{y}})", ylabel=L"\tilde{\mathrm{dist}}(\mathbf{\textit{s}}_2,\mathbf{\textit{y}})", legend=:topright, xguidefontsize=8, yguidefontsize=8, legendfontsize=5, tickfontsize=6, titlefontsize=8, dpi=600)
    scatter_test = scatter!(testdistsmulti_s₁_system2, testdistsmulti_s₂_system2, color=:plum4, alpha=0.5, ms=2, markerstrokewidth=0.1, label=L"\textrm{system \ 2}")
    
    final_plot = plot(S1, S2, scatter_test, layout=(1,3), margin=3mm, size=(800,200), dpi=600)
    return final_plot
end


function Figure11()
    series = [0.0 1.0; 4.0 2.0; 4.8 1.5; 8.0 2.0; 11.2 2.0]
    shapelet = [2.5 3.0; 6.0 3.0]
    fullshapelet = [0.75 2.0; 2.5 3.0; 6.0 2.0; 7.2 2.0]

    Times_combined = unique(sort(vcat(series[:,1], shapelet[:,1])))
    series_combinedvalues = zeros(length(Times_combined))
    shapelet_combinedvalues = zeros(length(Times_combined))
    shapeletindices = Vector{Int64}()
    for i in eachindex(series_combinedvalues)
        place1 = argmax(j -> series[:,1][j], findall(x -> x<=Times_combined[i], series[:,1]))
        series_combinedvalues[i] = series[place1,2]
        shapelet_combinedvalues[i] = series[place1,2]
        if Times_combined[i] >= shapelet[1,1] && Times_combined[i] < shapelet[end,1]
            place2 = argmax(j -> shapelet[:,1][j], findall(x -> x<=Times_combined[i], shapelet[:,1]))
            shapelet_combinedvalues[i] = shapelet[place2,2]
            shapeletindices = vcat(shapeletindices, i)
        end
    end
    series_combined = hcat(Times_combined, series_combinedvalues)
    shapelet_combined = hcat(Times_combined, shapelet_combinedvalues)
    only_shapelet_combined = shapelet_combined[shapeletindices,:]
    only_shapelet_combined = [only_shapelet_combined; shapelet[end,:]']

    proofdiagram = plot(Times_combined, series_combined[:,2], linetype=:steppost, color=:plum4, alpha=0.5, linewidth=2, label=:none, ylim=(0,5), xticks=([0.75, 1.5, 2.5, 4.0, 4.8, 6.0], [L" ",L" ",L"x+u_j",L"t_1",L"t_2",L"x+u_{j+1}"]), xguidefontsize=12, xtickfontsize=12, yguidefontsize=8, legendfontsize=9, tickfontsize=7, margins=3mm, size=(700,300), dpi=600)
    proofdiagram = plot!(Times_combined, series_combined[:,2], linetype=:steppost, color=:plum4, linewidth=2, label=L"y")
    proofdiagram = plot!(fullshapelet[:,1], fullshapelet[:,2], linetype=:steppost, color=:mediumseagreen, linestyle=:dash, label=L"s", linewidth=3)
    proofdiagram = plot!(only_shapelet_combined[:,1], only_shapelet_combined[:,2], linetype=:steppost, color=:mediumseagreen, label=L"s(u_j)", linewidth=3)
    proofdiagram = plot!(Times_combined, series_combined[:,2], fillrange=shapelet_combined[:,2], fillcolor=:mediumseagreen, fillalpha=0.2,  linetype=:steppost, color=:mediumseagreen, alpha=0.2, linewidth=0.1, label=L"\Vert s(u_j) - y([x+u_j, x+u_{j+1}]) \Vert_1")
    proofdiagram = [plot!([series[i,1],series[i,1]],[0,5],color=:plum4,linestyle=:dot, linewidth=2, label=:none) for i in 2:3]
    proofdiagram = [plot!([shapelet[i,1],shapelet[i,1]],[0,5],color=:mediumseagreen,linestyle=:dot, linewidth=2, label=:none) for i in 1:size(shapelet,1)]
    proofdiagram = plot!([0.75,1.5],[0.15,0.15],arrow=true,color=:black,linewidth=1,label="")
    proofdiagram = plot!([2.5,3.25],[2.8,2.8],arrow=true,color=:black,linewidth=1,label="")
    proofdiagram = plot!([6.0,6.75],[2.8,2.8],arrow=true,color=:black,linewidth=1,label="")
    proofdiagram = annotate!(1.125, 0.35, text(L"\delta", :black, 12))
    proofdiagram = annotate!(2.875, 2.6, text(L"\delta", :black, 12))
    proofdiagram = annotate!(6.375, 2.6, text(L"\delta", :black, 12))
    proofdiagram = annotate!(10.0, -0.35, text(L"\textrm{time}", :black, 10))
    proofdiagram = annotate!(0.75, -0.26, text(L"x", :black, 12))
    proofdiagram = annotate!(1.5, -0.21, text(L"x^\prime", :black, 12))
    return proofdiagram
end

function Figure12()
    series = [0.0 3.0; 4.0 1.0; 11.0 1.0]
    shapelet = [1.0 5.0; 3.0 2.0; 6.0 4.0; 7.5 4.0]

    Times_combined = unique(sort(vcat(series[:,1], shapelet[:,1])))
    series_combinedvalues = zeros(length(Times_combined))
    shapelet_combinedvalues = zeros(length(Times_combined))
    shapeletindices = Vector{Int64}()
    for i in eachindex(series_combinedvalues)
        place1 = argmax(j -> series[:,1][j], findall(x -> x<=Times_combined[i], series[:,1]))
        series_combinedvalues[i] = series[place1,2]
        shapelet_combinedvalues[i] = series[place1,2]
        if Times_combined[i] >= shapelet[1,1] && Times_combined[i] < shapelet[end,1]
            place2 = argmax(j -> shapelet[:,1][j], findall(x -> x<=Times_combined[i], shapelet[:,1]))
            shapelet_combinedvalues[i] = shapelet[place2,2]
            shapeletindices = vcat(shapeletindices, i)
        end
    end
    series_combined = hcat(Times_combined, series_combinedvalues)
    shapelet_combined = hcat(Times_combined, shapelet_combinedvalues)
    only_shapelet_combined = shapelet_combined[shapeletindices,:]
    only_shapelet_combined = [only_shapelet_combined; shapelet[end,:]']

    proofdiagram = plot(series[:,1], series[:,2], linetype=:steppost, color=:plum4, linewidth=2, ylim=(0,6), label=L"y", xticks=([1.0, 3.0, 4.0, 6.0, 7.5], [L"w^t_3",L"w^t_5",L"w^t_4",L"w^t_1",L"w^t_2"]), xguidefontsize=12, xtickfontsize=12, yguidefontsize=8, legendfontsize=9, tickfontsize=7, margins=3mm, size=(700,300), dpi=600)
    proofdiagram = plot!(shapelet[:,1], shapelet[:,2], linetype=:steppost, color=:mediumseagreen, linestyle=:dash, label=L"s", linewidth=3)
    proofdiagram = plot!(Times_combined, series_combined[:,2], fillrange=shapelet_combined[:,2], fillcolor=:mediumseagreen, fillalpha=0.2,  linetype=:steppost, color=:mediumseagreen, alpha=0.2, linewidth=0.1, label=L"\Vert s - y([t, t + \ell]) \Vert_1")
    proofdiagram = [plot!([series[i,1],series[i,1]],[0,6],color=:plum4,linestyle=:dot, linewidth=2, label=:none) for i in 2:2]
    proofdiagram = [plot!([shapelet[i,1],shapelet[i,1]],[0,6],color=:mediumseagreen,linestyle=:dot, linewidth=2, label=:none) for i in 1:size(shapelet,1)]
    proofdiagram = plot!([1.0,2.9],[0.15,0.15],arrow=true,color=:black,linewidth=1,label="")
    proofdiagram = plot!([3.0,3.9],[0.15,0.15],arrow=true,color=:black,linewidth=1,label="")
    proofdiagram = plot!([4.0,5.9],[0.15,0.15],arrow=true,color=:black,linewidth=1,label="")
    proofdiagram = plot!([6.0,7.4],[0.15,0.15],arrow=true,color=:black,linewidth=1,label="")
    proofdiagram = annotate!(2.0, 0.45, text(L"\lambda_3", :black, 12))
    proofdiagram = annotate!(3.5, 0.45, text(L"\lambda_5", :black, 12))
    proofdiagram = annotate!(5.0, 0.45, text(L"\lambda_4", :black, 12))
    proofdiagram = annotate!(6.75, 0.45, text(L"\lambda_1", :black, 12))
    proofdiagram = annotate!(8.2, 0.45, text(L"\lambda_2 = 0", :black, 12))
    proofdiagram = annotate!(1.0, 6.2, text(L"t", :black, 12))
    proofdiagram = annotate!(7.5, 6.2, text(L"t + \ell", :black, 12))
    proofdiagram = annotate!(10.0, -0.4, text(L"\textrm{time}", :black, 10))
    return proofdiagram
end