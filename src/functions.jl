function splitseriesWF_weeks(data, T, warmup)   
    #split data into trajectories of length T, discarding a warmup period
    #data is a vector, T and warmup are integer/floats
    collected_series = Vector{Matrix{Float64}}()
    data = vcat([0 0], data)   #start from 0
    n = floor(data[end,1]/T)
    for i in 1:n
        start = T*(i-1) + warmup
        finish = i*T
        series_start = argmin(data[:,1] .<= start) - 1
        series_end = argmin(data[:,1] .<= finish) - 1

        series = data[series_start:series_end,:]
        series[:,1] .-= start
        series[1,1] = 0.0
        series[:,2] .-= series[1,2]
        if series[end,1] != T
            series = vcat(series, [T-warmup series[end,2]])
        end
        push!(collected_series, series)
    end
    return collected_series
end

function splitseries(traj, T)  #validation
    collected_series = Vector{Matrix{Float64}}()
    traj = vcat([0 0], traj)  
    n = floor(traj[end,1]/T)
    for i in 2:2:n   #discarding warm up periods
        start = T*(i-1) - 20
        finish = i*T
        series_start = argmin(traj[:,1] .<= start) - 1
        series_end = argmin(traj[:,1] .<= finish) - 1
        series = traj[series_start:series_end,:]
        series[:,1] .-= start
        series[1,1] = 0.0
        if series[end,1] != T+20
            series = vcat(series, [T+20 series[end,2]])
        end
        push!(collected_series, series)
    end
    return collected_series
end

function splitmultiseriesV2(trajs, T)  #multi
    collected_series = Vector{Vector{Matrix{Float64}}}()
    d = length(trajs)
    n = floor(trajs[1][end,1]/(T+20))
    for i in 1:d
        multiseries = Vector{Matrix{Float64}}()
        traj = vcat([0 0], trajs[i])   #start from 0
        for i in 1:n
            start = (T+20)*(i-1) + 20
            finish = start + T
            series_start = argmin(traj[:,1] .<= start) - 1
            series_end = argmin(traj[:,1] .<= finish) - 1
            series = traj[series_start:series_end,:]
            series[:,1] .-= start
            series[1,1] = 0.0
            if series[end,1] != T
                series = vcat(series, [T series[end,2]])
            end
            push!(multiseries, series)
        end
        push!(collected_series, multiseries)
    end
    return collected_series
end

function update_distance(series1, series2, currentmin; LocInv::LocInvMethod=Yes())
    #updates the minimum distance as a shapelet moves along a series.
    #series 1 is the shapelet, series 2 is the segment of the series being tested. 
    #currentmin holds the current minimum. The distance calculation is abandoned as soon as this is exceeded
    Times_combined = unique(sort(vcat(series1[:,1], series2[:,1])))
    length_Times_combined = length(Times_combined)
    series1_combinedvalues = Vector{Float64}(undef, length_Times_combined)
    series2_combinedvalues = Vector{Float64}(undef, length_Times_combined)
    j = 1
    for i in 1:(size(series1, 1)-1)
        while Times_combined[j] < series1[i+1, 1]
            series1_combinedvalues[j] = series1[i, 2]
            j += 1
        end
    end
    series1_combinedvalues[end] = series1[end,2]
    j = 1
    for i in 1:(size(series2, 1)-1)
        while Times_combined[j] < series2[i+1, 1]
            series2_combinedvalues[j] = series2[i, 2]
            j += 1
        end
    end
    series2_combinedvalues[end] = series2[end,2]

    if LocInv isa Yes
        differences = (series1_combinedvalues .- series2_combinedvalues)[1:(end-1)]
        ordering = sortperm(differences)
        widths = [Times_combined[i+1] - Times_combined[i] for i in 1:(length(Times_combined)-1)]
        #ordered_widths = widths[ordering]
        median_position = argmin(cumsum(widths[ordering]) .< series1[end,1]/2)
        c = differences[ordering][median_position]
        differences = abs.((series2_combinedvalues.+c) .- series1_combinedvalues)[1:(end-1)]
        distance = 0.0
        i = 1
        while distance < currentmin && i < length_Times_combined
            distance += widths[i] * differences[i]
            i += 1
        end
        return minimum([distance, currentmin])
    end
    differences = abs.(series1_combinedvalues - series2_combinedvalues)
    distance = 0.0
    i = 1
    while distance < currentmin && i <= length(differences)
        distance += (Times_combined[i+1] - Times_combined[i]) * differences[i]
        i += 1
    end
    return minimum([distance, currentmin])
end

function dist_shapelet_series(shapelet, series; LocInv::LocInvMethod=Yes())
    #calculate the distance between a shapelet and a series
    l = shapelet[end,1]  
    m = series[end,1]  
    snap_points = Array{Float64}(undef, 0)
    for i in shapelet[:,1]
        for j in series[:,1]
            if i <= j && l-i+j <= m 
                snap_points = vcat(snap_points, j-i)
            end
        end
    end
    snap_points = unique(sort(snap_points))
    no_points = length(snap_points)
    distance = 1000000.0
    j = 1
    for i in 1:(size(series, 1) - 1)
        while j <= no_points && snap_points[j] <= series[i+1,1]
            series_end = argmin(series[i:end, 1] .< snap_points[j]+l) + i - 2
            series_segment = series[i:series_end,:]
            series_segment[:,1] .-= snap_points[j]
            series_segment[1,1] = 0.0
            if series_segment[end,1] != l
                series_segment = vcat(series_segment, [l series_segment[end,2]])
            end
            distance = update_distance(shapelet, series_segment, distance, LocInv=LocInv)
            j += 1
        end
    end
    return round(distance, digits=8)
end

function OptimalSplit(class1dists, class2dists)
    #Given vectors of distances of a shapelet to class 1 series and class 2 series, return the associated information gain with the optimal splitting point and the difference between the mean distance to class 1 series and the mean distance to class 2 series
    combined_unique = unique(sort([class1dists; class2dists]))
    combined = [class1dists; class2dists]
    initialproportion = length(class1dists)/length(combined)
    initialentropy = entropy([initialproportion, 1 - initialproportion])
    best_gain = 0
    best_split = 0
    for i in 1:(length(combined_unique)-1)
        split = combined_unique[i]
        class1propD₁ = sum(class1dists .<= split) / sum(combined .<= split)
        class1propD₂ = sum(class1dists .> split) / sum(combined .> split)
        gain = initialentropy - (((sum(combined .<= split) / length(combined))*entropy([class1propD₁, 1 - class1propD₁])) + ((sum(combined .> split) / length(combined))*entropy([class1propD₂, 1 - class1propD₂])))
        if gain > best_gain
            best_gain = gain
            best_split = split
        end
    end
    mean_dist = abs(mean(class1dists) - mean(class2dists))
    return best_gain, best_split, mean_dist
end

function OptimalGain(class1dists, class2dists)  
    #Given vectors of distances of a shapelet to class 1 series and class 2 series, return the optimal information gain
    #This is the same as OptimalSplit() but only returning the information gain
    combined_unique = unique(sort([class1dists; class2dists]))
    combined = [class1dists; class2dists]
    initialproportion = length(class1dists)/length(combined)
    initialentropy = entropy([initialproportion, 1 - initialproportion])
    best_gain = 0
    for i in 1:(length(combined_unique)-1)
        split = combined_unique[i]
        class1propD₁ = sum(class1dists .<= split) / sum(combined .<= split)
        class1propD₂ = sum(class1dists .> split) / sum(combined .> split)
        gain = initialentropy - (((sum(combined .<= split) / length(combined))*entropy([class1propD₁, 1 - class1propD₁])) + ((sum(combined .> split) / length(combined))*entropy([class1propD₂, 1 - class1propD₂])))
        if gain > best_gain
            best_gain = gain
        end
    end
    return best_gain
end

function TestShapelet(shapelet, class1series, class2series, current_best_gain; LocInv::LocInvMethod=Yes())
    #test if a shapelet can be better than the current best, by iteratively calculating distances to series and updating the best possible information gain. As soon as we have calculated enough distances to know that we can't improve on the information gain from the current best shapelet, we abandon this test shapelet
    best_possible = 100000.0
    
    n1 = length(class1series)
    n2 = length(class2series)
    no_series = n1 + n2

    class1dists = dist_shapelet_series(shapelet, class1series[1], LocInv=LocInv)
    class2dists = dist_shapelet_series(shapelet, class2series[1], LocInv=LocInv)

    class1_counter = 1   #counts how many class 1 dists have been calculated
    class2_counter = 1   #counts how many class 2 dists have been calculated

    while best_possible >= current_best_gain && (class1_counter + class2_counter) < no_series
        if class1_counter < n1
            new_class1dist = dist_shapelet_series(shapelet, class1series[class1_counter+1], LocInv=LocInv)
            class1dists = vcat(class1dists, new_class1dist)
            class1_counter += 1
        end
        if class2_counter < n2
            new_class2dist = dist_shapelet_series(shapelet, class2series[class2_counter+1], LocInv=LocInv)
            class2dists = vcat(class2dists, new_class2dist)
            class2_counter += 1
        end

        best_possible = bestPossible(class1dists, class2dists, n1, n2)
        
    end
    return best_possible, class1dists, class2dists
    #If all distances were calculated, then best_possible is the actual infortmation gain
end

function bestPossible(class1dists, class2dists, n1, n2)
    #Given the already-claculated distances to class 1 and class 2 series, and the total number of class 1 and class 2 series, calculate the best possible information gain, considering the not-yet-calculated distances
    alldists = vcat(class1dists, class2dists)
    mindist = minimum(alldists)
    maxdist = maximum(alldists)
    
    class1dists_try1 = vcat(class1dists, repeat([mindist], n1-length(class1dists)))
    class2dists_try1 = vcat(class2dists, repeat([maxdist], n2-length(class2dists)))
    bestPossible_try1 = OptimalGain(class1dists_try1, class2dists_try1)

    class1dists_try2 = vcat(class1dists, repeat([maxdist], n1-length(class1dists)))
    class2dists_try2 = vcat(class2dists, repeat([mindist], n2-length(class2dists)))
    bestPossible_try2 = OptimalGain(class1dists_try2, class2dists_try2)

    return maximum([bestPossible_try1, bestPossible_try2])
end

function FindShapelet(class1series, class2series, l, lrange, τ; shapeletrange=(0,class1series[1][end,1]-l), LocInv::LocInvMethod=Yes(), Search::SearchOver=All())  #current_gain=0.0 
    #l is the target length of shapelets to look for.
    #lrange is the full set of lengths to consider
    #τ is the time steps for shapelet start times
    #shapeletrange is the range of starting positions of shapelets. Default is over the whole trajectory length (0,m-l)
    #class1series and class2series are vectors of vectors
    if Search isa All
        searchseries = vcat(class1series, class2series)
    elseif Search isa Class1
        searchseries = class1series
    elseif Search isa Class2
        searchseries = class2series
    end
    #m = searchseries[1][end,1]  #find trajectory lengths (assume all the same)
    best_gain = 0.0 
    shapelet_meandist = 0.0  
    shapelet_threshold = 0.0
    best_shapelet = vcat(searchseries[1][1:(argmin(searchseries[1][:,1] .<= l)-1),:], [l searchseries[1][argmin(searchseries[1][:,1] .< l)-1, 2]])
    perfect_gain = OptimalGain(zeros(length(class1series)), ones(length(class2series)))  #if we find a shapelet with perfect information gain, we will stop the search
    @showprogress "Shapelets Search: " for i in 1:length(searchseries)
        for j in shapeletrange[1]:τ:shapeletrange[2]
            shapelet_start = argmin(searchseries[i][:,1] .<= j) - 1
            shapelet_end = argmin(searchseries[i][:,1] .<= j + l) - 1
            shapelet_end = ifelse(shapelet_end == 0, length(searchseries[i][:,1]), shapelet_end) ###
            shapelet = searchseries[i][shapelet_start:shapelet_end,:]

            shapelet[:,1] .-= j
            shapelet[1,1] = 0.0
            shapelet = vcat(shapelet, [l shapelet[end,2]])

            test = TestShapelet(shapelet, class1series, class2series, best_gain; LocInv=LocInv)
            test_gain = test[1]
            if Search isa Class1
                if test_gain > best_gain && mean(test[2]) < mean(test[3])
                    println("new shapelet found via I")
                    localsearch = LocalShapeletSearch(class1series, class2series, lrange, τ, LocInv=LocInv, series=searchseries[i], position=j, current_gain=best_gain)
                    best_shapelet = localsearch[1]
                    best_gain = localsearch[2]
                    shapelet_threshold = localsearch[3]
                    shapelet_meandist = localsearch[4]

                elseif test_gain == best_gain
                    localsearch = LocalShapeletSearch(class1series, class2series, lrange, τ, LocInv=LocInv, series=searchseries[i], position=j, current_gain=best_gain)
                    if localsearch[2] > best_gain
                        best_shapelet = localsearch[1]
                        best_gain = localsearch[2]
                        shapelet_threshold = localsearch[3]
                        shapelet_meandist = localsearch[4]
                        println("new shapelet found via I")
                    elseif localsearch[4] > shapelet_meandist
                        best_shapelet = localsearch[1]
                        best_gain = localsearch[2]
                        shapelet_threshold = localsearch[3]
                        shapelet_meandist = localsearch[4]
                        println("new shapelet found via mean dist")
                    end
                end
            elseif Search isa Class2
                if test_gain > best_gain && mean(test[3]) < mean(test[2])
                    println("new shapelet found via I")
                    localsearch = LocalShapeletSearch(class1series, class2series, lrange, τ, LocInv=LocInv, series=searchseries[i], position=j, current_gain=best_gain)
                    best_shapelet = localsearch[1]
                    best_gain = localsearch[2]
                    shapelet_threshold = localsearch[3]
                    shapelet_meandist = localsearch[4]

                elseif test_gain == best_gain
                    localsearch = LocalShapeletSearch(class1series, class2series, lrange, τ, LocInv=LocInv, series=searchseries[i], position=j, current_gain=best_gain)
                    if localsearch[2] > best_gain
                        best_shapelet = localsearch[1]
                        best_gain = localsearch[2]
                        shapelet_threshold = localsearch[3]
                        shapelet_meandist = localsearch[4]
                        println("new shapelet found via I")
                    elseif localsearch[4] > shapelet_meandist
                        best_shapelet = localsearch[1]
                        best_gain = localsearch[2]
                        shapelet_threshold = localsearch[3]
                        shapelet_meandist = localsearch[4]
                        println("new shapelet found via mean dist")
                    end
                end
            else
                if test_gain > best_gain
                    best_gain = test_gain
                    println("new shapelet found via I")
                    localsearch = LocalShapeletSearch(class1series, class2series, lrange, τ, LocInv=LocInv, series=searchseries[i], position=j, current_gain=best_gain)
                    best_shapelet = localsearch[1]
                    best_gain = localsearch[2]
                    shapelet_threshold = localsearch[3]
                    shapelet_meandist = localsearch[4]
    
                elseif test_gain == best_gain
                    localsearch = LocalShapeletSearch(class1series, class2series, lrange, τ, LocInv=LocInv, series=searchseries[i], position=j, current_gain=best_gain)
                    if localsearch[2] > best_gain
                        best_shapelet = localsearch[1]
                        best_gain = localsearch[2]
                        shapelet_threshold = localsearch[3]
                        shapelet_meandist = localsearch[4]
                        println("new shapelet found via I")
                    elseif localsearch[4] > shapelet_meandist
                        best_shapelet = localsearch[1]
                        best_gain = localsearch[2]
                        shapelet_threshold = localsearch[3]
                        shapelet_meandist = localsearch[4]
                        println("new shapelet found via mean dist")
                    end
                end
            end
            if best_gain == perfect_gain
                return best_shapelet, best_gain, shapelet_threshold, shapelet_meandist
            end
        end
    end
    return best_shapelet, best_gain, shapelet_threshold, shapelet_meandist#, shapelet_series, shapelet_positions
end

function LocalShapeletSearch(class1series, class2series, lrange, τ; LocInv::LocInvMethod=Yes(), series, position, current_gain)
    m = series[end,1]  #find time series length (assume all the same)
    best_gain = current_gain 
    shapelet_meandist = 0.0  
    shapelet_threshold = 0.0
    best_shapelet = series
    perfect_gain = OptimalGain(zeros(length(class1series)), ones(length(class2series)))
    @showprogress "Local Search: " for l in lrange
        for j in maximum([0,position-abs(l - mean(lrange))]):τ:minimum([position+abs(l - mean(lrange)),m-l]) 
            shapelet_start = argmin(series[:,1] .<= j) - 1
            shapelet_end = argmin(series[:,1] .< j + l) - 1
            shapelet = series[shapelet_start:shapelet_end,:]
            shapelet[:,1] .-= j
            shapelet[1,1] = 0.0
            shapelet = vcat(shapelet, [l shapelet[end,2]])

            test = TestShapelet(shapelet, class1series, class2series, best_gain; LocInv=LocInv)
            test_gain = test[1]

            if test_gain > best_gain
                best_shapelet = shapelet
                class1dists = test[2]
                class2dists = test[3]
                test_split = OptimalSplit(class1dists, class2dists)

                shapelet_threshold = test_split[2]
                shapelet_meandist = test_split[3]
                best_gain = test_gain
                println("new shapelet found via I")
            elseif test_gain == best_gain
                class1dists = test[2]
                class2dists = test[3]
                test_split = OptimalSplit(class1dists, class2dists)
                if test_split[3] > shapelet_meandist   
                    best_shapelet = shapelet
                    shapelet_threshold = test_split[2]
                    shapelet_meandist = test_split[3]
                    println("new shapelet found via mean dist")
                end
            end
            if best_gain == perfect_gain
                return best_shapelet, best_gain, shapelet_threshold, shapelet_meandist
            end
        end
    end
    return best_shapelet, best_gain, shapelet_threshold, shapelet_meandist
end



#Functions for multivariate shapelets
function multi_update_distance(series1, series2, currentmin; LocInv::LocInvMethod=Yes())
    distance = 0.0
    for d in 1:length(series1) 
        Times_combined = unique(sort(vcat(series1[d][:,1], series2[d][:,1])))
        length_Times_combined = length(Times_combined)
        series1_combinedvalues = zeros(length_Times_combined)
        series2_combinedvalues = zeros(length_Times_combined)
        j = 1
        for i in 1:(size(series1[d], 1)-1)
            while Times_combined[j] < series1[d][i+1, 1]
                series1_combinedvalues[j] = series1[d][i, 2]
                j += 1
            end
        end
        series1_combinedvalues[end] = series1[d][end,2]
        j = 1
        for i in 1:(size(series2[d], 1)-1)
            while Times_combined[j] < series2[d][i+1, 1]
                series2_combinedvalues[j] = series2[d][i, 2]
                j += 1
            end
        end
        series2_combinedvalues[end] = series2[d][end,2]
        if LocInv isa Yes
            differences = (series1_combinedvalues .- series2_combinedvalues)[1:(end-1)]
            ordering = sortperm(differences)
            widths = [Times_combined[i+1] - Times_combined[i] for i in 1:(length(Times_combined)-1)]
            median_position = argmin(cumsum(widths[ordering]) .< series1[d][end,1]/2)
            c = differences[ordering][median_position]
            differences = abs.((series2_combinedvalues.+c) .- series1_combinedvalues)[1:(end-1)]
            
            i = 1
            while distance < currentmin && i < length_Times_combined
                distance += widths[i] * differences[i]
                i += 1
            end
            
        end
        
    end
    return minimum([distance, currentmin])
end

function multidep_dist_shapelet_series(multishapelet, series; LocInv::LocInvMethod=Yes())
    #calculate the distance between a multivariate shapelet and a multivariate series
    l = multishapelet[1][end,1]   
    m = series[1][end,1]    
    distances = Vector{Float64}()
    snap_points = Array{Float64}(undef, 0)
    mindist_position = 0.0
    for d in 1:length(multishapelet)
        for i in multishapelet[d][:,1]
            for j in series[d][:,1]
                if i <= j && l-i+j <= m 
                    snap_points = vcat(snap_points, j-i)
                end
            end
        end
    end
    snap_points = unique(sort(snap_points))
    no_points = length(snap_points)
    distance = 1000000.0

    for j in snap_points
        multiseries_segment = Vector{Matrix{Float64}}()
        for d in 1:length(series)
            series_start = argmin(series[d][:, 1] .<= j) -1
            series_end = argmin(series[d][:, 1] .< j+l) - 1
            series_segment = series[d][series_start:series_end,:]
            series_segment[:,1] .-= j
            series_segment[1,1] = 0.0
            series_segment = vcat(series_segment, [l series_segment[end,2]])
            push!(multiseries_segment, series_segment)
        end
        testdistance = multi_update_distance(multishapelet, multiseries_segment, distance, LocInv=LocInv)
        if testdistance < distance
            mindist_position = j
        end
        distance = testdistance
    end
    return round(distance, digits=8), mindist_position
end

function MultiTestShapelet(multishapelet, class1series, class2series, current_best_gain; LocInv::LocInvMethod=Yes())
    #multivariate version of TestShapelet()
    best_possible = 100000.0
    
    n1 = length(class1series[1])
    n2 = length(class2series[1])
    no_series = n1 + n2

    class1dists = multidep_dist_shapelet_series(multishapelet, [class1series[d][1] for d in 1:length(class1series)], LocInv=LocInv)[1]
    class2dists = multidep_dist_shapelet_series(multishapelet, [class2series[d][1] for d in 1:length(class2series)], LocInv=LocInv)[1]

    class1_counter = 1   #counts how many class 1 dists have been calculated
    class2_counter = 1   #counts how many class 2 dists have been calculated

    while best_possible >= current_best_gain && (class1_counter + class2_counter) < no_series
        if class1_counter < n1
            new_class1dist = multidep_dist_shapelet_series(multishapelet, [class1series[d][class1_counter+1] for d in 1:length(class1series)], LocInv=LocInv)[1]
            class1dists = vcat(class1dists, new_class1dist)
            class1_counter += 1
        end
        if class2_counter < n2
            new_class2dist = multidep_dist_shapelet_series(multishapelet, [class2series[d][class2_counter+1] for d in 1:length(class2series)], LocInv=LocInv)[1]
            class2dists = vcat(class2dists, new_class2dist)
            class2_counter += 1
        end

        best_possible = bestPossible(class1dists, class2dists, n1, n2)
        
    end
    return best_possible, class1dists, class2dists
    #If all distances were calculated, then best_possible is the actual infortmation gain
end

function MultiDepShapelet(class1series, class2series, l, lrange, τ; LocInv::LocInvMethod=Yes(), Search::SearchOver=All())  #current_gain=0.0 
    #l is the length of shapelets to look for. lrange is the plus/minus around this to consider
    #class1series and class2series are length d vectors of length n vectors of series
    if Search isa All
        searchseries = vcat(class1series[d], class2series[d])  #needs fixing. 
    elseif Search isa Class1
        searchseries = class1series
    elseif Search isa Class2
        searchseries = class2series
    end
    n = length(searchseries[1])
    m = searchseries[1][1][end,1]  #find time series length (assume all the same)
    best_gain = 0.0  #changed from = 0.0, and added current_gain as an argument. #to store the information gains of the best shapelets
    shapelet_meandist = 0.0  
    shapelet_threshold = 0.0
    best_shapelet = searchseries[1][1]
    perfect_gain = OptimalGain(zeros(n), ones(n))
    #shapelets = Vector{Matrix{Float64}}()   #to store the best shapelets
    @showprogress "Shapelets Search: " for i in 1:n
        #series_consider = length(findall(x -> x<=m-l, searchseries[i][:,1]))
        for j in 0:τ:(m-l)
            multishapelet = Vector{Matrix{Float64}}()
            for d in 1:length(searchseries)
                shapelet_start = argmin(searchseries[d][i][:, 1] .<= j) -1
                shapelet_end = argmin(searchseries[d][i][:, 1] .< j+l) -1
                shapelet = searchseries[d][i][shapelet_start:shapelet_end,:]
                shapelet[1,1] = j
                shapelet[:,1] .-= j
                shapelet = vcat(shapelet, [l shapelet[end,2]])
                push!(multishapelet, shapelet)
            end
            test = MultiTestShapelet(multishapelet, class1series, class2series, best_gain; LocInv=LocInv)
            test_gain = test[1]

            if test_gain > best_gain
                best_gain = test_gain
                println("new shapelet found via I")
                localsearch = MultiLocalShapeletSearch(class1series, class2series, lrange, LocInv=LocInv, series=[searchseries[d][i] for d in 1:length(searchseries)], position=j, current_gain=best_gain)
                best_shapelet = localsearch[1]
                best_gain = localsearch[2]
                shapelet_threshold = localsearch[3]
                shapelet_meandist = localsearch[4]

            elseif test_gain == best_gain
                localsearch = MultiLocalShapeletSearch(class1series, class2series, lrange, LocInv=LocInv, series=[searchseries[d][i] for d in 1:length(searchseries)], position=j, current_gain=best_gain)
                if localsearch[2] > best_gain
                    best_shapelet = localsearch[1]
                    best_gain = localsearch[2]
                    shapelet_threshold = localsearch[3]
                    shapelet_meandist = localsearch[4]
                    println("new shapelet found via I")
                elseif localsearch[4] > shapelet_meandist
                    best_shapelet = localsearch[1]
                    best_gain = localsearch[2]
                    shapelet_threshold = localsearch[3]
                    shapelet_meandist = localsearch[4]
                    #println("new shapelet found via mean dist")
                end

            end
            if best_gain == perfect_gain
                return best_shapelet, best_gain, shapelet_threshold, shapelet_meandist
            end
        end
    end
    return best_shapelet, best_gain, shapelet_threshold, shapelet_meandist#, shapelet_series, shapelet_positions
end

function MultiLocalShapeletSearch(class1series, class2series, lrange; LocInv::LocInvMethod=Yes(), series, position, current_gain)
    m = series[1][end,1]  #find time series length (assume all the same)
    best_gain = current_gain   #to store the information gains of the best shapelets
    shapelet_meandist = 0.0  
    shapelet_threshold = 0.0
    best_shapelet = series
    perfect_gain = OptimalGain(zeros(length(class1series[1])), ones(length(class2series[1])))
    for l in lrange
        for j in maximum([0,position-1]):0.2:minimum([m-l,position+1]) 
            multishapelet = Vector{Matrix{Float64}}()
            for d in 1:length(series)
                shapelet_start = argmin(series[d][:, 1] .<= j) -1
                shapelet_end = argmin(series[d][:, 1] .< j+l) -1
                shapelet = series[d][shapelet_start:shapelet_end,:]
                shapelet[1,1] = j
                shapelet[:,1] .-= j
                shapelet = vcat(shapelet, [l shapelet[end,2]])
                push!(multishapelet, shapelet)
            end
            test = MultiTestShapelet(multishapelet, class1series, class2series, best_gain; LocInv=LocInv)
            test_gain = test[1]

            if test_gain > best_gain
                best_shapelet = multishapelet
                class1dists = test[2]
                class2dists = test[3]
                test_split = OptimalSplit(class1dists, class2dists)

                shapelet_threshold = test_split[2]
                shapelet_meandist = test_split[3]
                best_gain = test_gain
                println("new shapelet found via I")
            elseif test_gain == best_gain
                class1dists = test[2]
                class2dists = test[3]
                test_split = OptimalSplit(class1dists, class2dists)
                if test_split[3] > shapelet_meandist    #test_meandist > shapelet_meandist
                    best_shapelet = multishapelet
                    shapelet_threshold = test_split[2]
                    shapelet_meandist = test_split[3]
                end
            end
            if best_gain == perfect_gain
                return best_shapelet, best_gain, shapelet_threshold, shapelet_meandist
            end
        end
    end
    return best_shapelet, best_gain, shapelet_threshold, shapelet_meandist
end