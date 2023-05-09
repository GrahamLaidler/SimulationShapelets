

####functions for distance calculations and optimal information gain####

function update_distance(series1, series2, currentmin; LocInv::LocInvMethod=Yes())
    #updates the minimum distance as a shapelet moves along a series.
    #series 1 is the shapelet, series 2 is the segment of the series being tested. These are matrices - column 1 contains the change times and column 2 the values
    #currentmin holds the current minimum. The distance calculation is abandoned as soon as this is exceeded, and we stay with currentmin
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
        widths = [Times_combined[i+1] - Times_combined[i] for i in 1:(length_Times_combined-1)]
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
    ℓ = shapelet[end,1]  
    m = series[end,1]  
    V = Array{Float64}(undef, 0) #the set containing shapelet start times for which a change in the shapelet will coincide with a change in the series. \mathcal{V} in the theorems 
    for i in shapelet[:,1]
        for j in series[:,1]
            if i <= j && ℓ-i+j <= m 
                V = vcat(V, j-i)
            end
        end
    end
    V = unique(sort(V))
    no_points = length(V)
    distance = 1000000.0  #this value must exceed the minimum distance
    j = 1
    for i in 1:(size(series, 1) - 1)
        while j <= no_points && V[j] <= series[i+1,1]
            series_end = argmin(series[i:end, 1] .< V[j]+ℓ) + i - 2
            series_segment = series[i:series_end,:]
            series_segment[:,1] .-= V[j]
            series_segment[1,1] = 0.0
            if series_segment[end,1] != ℓ
                series_segment = vcat(series_segment, [ℓ series_segment[end,2]])
            end
            distance = update_distance(shapelet, series_segment, distance, LocInv=LocInv)
            j += 1
        end
    end
    return round(distance, digits=8)
end

function OptimalSplit(class1dists, class2dists)
    #Given vectors of distances of a shapelet to class 1 series and class 2 series, return the information gain associated with the optimal split point, the optimal split point, and the difference between the mean distance to class 1 series and the mean distance to class 2 series
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

####functions used for admissible entropy pruning####

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
    #test if a shapelet can be better than the current best, by iteratively calculating distances to series and updating the best possible information gain.
    #As soon as we have calculated enough distances to know that we can't improve on the information gain from the current best shapelet, we can abandon this test shapelet
    #shapelet is a matrix - column 1 contains the change times and column 2 the values
    #class1series and class2series are vectors of matrices (of the same form as shapelet)
    #current_best_gain is the optimal information gain from the best shapelet found so far
    best_possible = 100000.0
    
    n₁ = length(class1series)  
    n₂ = length(class2series) 
    n = n₁ + n₂

    class1dists = dist_shapelet_series(shapelet, class1series[1], LocInv=LocInv)  #calculate the distance to the first class 1 series
    class2dists = dist_shapelet_series(shapelet, class2series[1], LocInv=LocInv)  #calculate the distance to the first class 2 series

    class1_counter = 1   #counts how many class 1 dists have been calculated
    class2_counter = 1   #counts how many class 2 dists have been calculated

    while best_possible >= current_best_gain && (class1_counter + class2_counter) < n
        if class1_counter < n₁
            new_class1dist = dist_shapelet_series(shapelet, class1series[class1_counter+1], LocInv=LocInv)
            class1dists = vcat(class1dists, new_class1dist)
            class1_counter += 1
        end
        if class2_counter < n₂
            new_class2dist = dist_shapelet_series(shapelet, class2series[class2_counter+1], LocInv=LocInv)
            class2dists = vcat(class2dists, new_class2dist)
            class2_counter += 1
        end

        best_possible = bestPossible(class1dists, class2dists, n₁, n₂)
        
    end
    return best_possible, class1dists, class2dists
    #If all distances were calculated, then best_possible is the actual infortmation gain
end

function bestPossible(class1dists, class2dists, n₁, n₂)
    #Given the already-claculated distances to class 1 and class 2 series, and the total number of class 1 and class 2 series, calculate the best possible information gain, considering the possible arrangements of the not-yet-calculated distances
    alldists = vcat(class1dists, class2dists)
    mindist = minimum(alldists)
    maxdist = maximum(alldists)
    
    class1dists_option1 = vcat(class1dists, repeat([mindist], n₁-length(class1dists)))
    class2dists_option1 = vcat(class2dists, repeat([maxdist], n₂-length(class2dists)))
    bestPossible_option1 = OptimalGain(class1dists_option1, class2dists_option1)

    class1dists_option2 = vcat(class1dists, repeat([maxdist], n₁-length(class1dists)))
    class2dists_option2 = vcat(class2dists, repeat([mindist], n₂-length(class2dists)))
    bestPossible_option2 = OptimalGain(class1dists_option2, class2dists_option2)

    return maximum([bestPossible_option1, bestPossible_option2])
end

####Shapelet finding functions####

function FindShapelet(class1series, class2series, ℓ, lrange, τ; shapeletrange=(0,class1series[1][end,1]-ℓ), LocInv::LocInvMethod=Yes(), Search::SearchOver=All()) 
    #ℓ is the target length of shapelets to look for.
    #lrange is the full set of lengths to consider (\mathcal{I} in paper)
    #τ is the time steps for shapelet extractions
    #shapeletrange is the range of starting positions of shapelets. Default is over the whole trajectory length (0,m-ℓ)
    #class1series and class2series are vectors of matrices - for each matrix, column 1 contains the change times and column 2 the values
    #LocInv = Yes() uses location invariant distance function (default), otherwise use LocInv = No()
    #if Search = Class1() or Class2(), we find the optimal characteristic shapelet of class 1 or 2, respectively. If Search = All(), we find the overall optimal shapelet from either class.
    if Search isa All
        searchseries = vcat(class1series, class2series)
    elseif Search isa Class1
        searchseries = class1series
    elseif Search isa Class2
        searchseries = class2series
    end
    best_gain = 0.0 
    shapelet_meandist = 0.0  
    shapelet_threshold = 0.0
    best_shapelet = vcat(searchseries[1][1:(argmin(searchseries[1][:,1] .<= ℓ)-1),:], [ℓ searchseries[1][argmin(searchseries[1][:,1] .< ℓ)-1, 2]])
    perfect_gain = OptimalGain(zeros(length(class1series)), ones(length(class2series)))  #if we find a shapelet which gives perfect discrimination, we will stop the search
    @showprogress "Shapelets Search: " for i in eachindex(searchseries)
        for j in shapeletrange[1]:τ:shapeletrange[2]
            ####extract candidate shapelet####
            shapelet_start = argmin(searchseries[i][:,1] .<= j) - 1
            shapelet_end = argmin(searchseries[i][:,1] .<= j + ℓ) - 1
            shapelet_end = ifelse(shapelet_end == 0, length(searchseries[i][:,1]), shapelet_end) ###
            shapelet = searchseries[i][shapelet_start:shapelet_end,:]
            shapelet[:,1] .-= j
            shapelet[1,1] = 0.0
            shapelet = vcat(shapelet, [ℓ shapelet[end,2]])
            ####test candidate shapelet with admissible entropy pruning####
            test = TestShapelet(shapelet, class1series, class2series, best_gain; LocInv=LocInv)
            test_gain = test[1]
            if test_gain >= best_gain  #if the candidate is at least as good as the current best shapelet, move to a local shapelet search
                if Search isa Class1 || Search isa Class2  #if looking specifically for characteristic shapelets of class 1 or class 2, check the subset proportions to validate this
                    Optsplit = OptimalSplit(test[2], test[3])
                    γ_star = Optsplit[2]
                    subset_split = vcat(test[2], test[3]) .<= γ_star
                    class1_split = test[2] .<= γ_star
                    class2_split = test[3] .<= γ_star
                    prop1near = sum(class1_split)/sum(subset_split)
                    prop1far = (length(class1_split) - sum(class1_split))/(length(subset_split) - sum(subset_split))
                    prop2near = sum(class2_split)/sum(subset_split)
                    prop2far = (length(class2_split) - sum(class2_split))/(length(subset_split) - sum(subset_split))
                    if (Search isa Class1 && prop1near > prop1far) || (Search isa Class2 && prop2near > prop2far)
                        localsearch = LocalShapeletSearch(class1series, class2series, lrange, τ, LocInv=LocInv, series=searchseries[i], position=j, current_gain=best_gain)
                        best_shapelet = localsearch[1]
                        best_gain = localsearch[2]
                        shapelet_threshold = localsearch[3]
                        shapelet_meandist = localsearch[4]
                    end
                else  #if Search isa All
                    localsearch = LocalShapeletSearch(class1series, class2series, lrange, τ, LocInv=LocInv, series=searchseries[i], position=j, current_gain=best_gain)
                    best_shapelet = localsearch[1]
                    best_gain = localsearch[2]
                    shapelet_threshold = localsearch[3]
                    shapelet_meandist = localsearch[4]
                end
            end

            if best_gain == perfect_gain
                return best_shapelet, best_gain, shapelet_threshold, shapelet_meandist
            end
        end
    end
    return best_shapelet, best_gain, shapelet_threshold, shapelet_meandist
end

function LocalShapeletSearch(class1series, class2series, lrange, τ; LocInv::LocInvMethod=Yes(), series, position, current_gain)
    #For use within FindShapelet()
    m = series[end,1]  #find time series length (assume all the same)
    best_gain = current_gain 
    shapelet_meandist = 0.0  
    shapelet_threshold = 0.0
    best_shapelet = series
    perfect_gain = OptimalGain(zeros(length(class1series)), ones(length(class2series)))
    @showprogress "Local Search: " for ℓ in lrange
        for j in maximum([0,position-abs(ℓ - mean(lrange))]):τ:minimum([position+abs(ℓ - mean(lrange)),m-ℓ]) 
            ####extract candidate shapelet####
            shapelet_start = argmin(series[:,1] .<= j) - 1
            shapelet_end = argmin(series[:,1] .< j + ℓ) - 1
            shapelet = series[shapelet_start:shapelet_end,:]
            shapelet[:,1] .-= j
            shapelet[1,1] = 0.0
            shapelet = vcat(shapelet, [ℓ shapelet[end,2]])
            ####test candidate shapelet with admissible entropy pruning####
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
                if test_split[3] > shapelet_meandist   #break ties by maximum class separation  ..test_meandist > shapelet_meandist
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



####Equivalent functions for multivariate shapelets####

function multi_update_distance(series1, series2, currentmin; LocInv::LocInvMethod=Yes())
    #updates the minimum multivariate distance as a shapelet moves along a series.
    #series 1 is the shapelet, series 2 is the segment of the series being tested. These are vectors of matrices. The first vector element is the first dimension etc. For each matrix, column 1 contains the change times and column 2 the values. 
    #currentmin holds the current minimum. The distance calculation is abandoned as soon as this is exceeded, and we stay with currentmin
    distance = 0.0
    for d in eachindex(series1) 
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
        else
            differences = abs.(series1_combinedvalues - series2_combinedvalues)
            i = 1
            while distance < currentmin && i <= length(differences)
                distance += (Times_combined[i+1] - Times_combined[i]) * differences[i]
                i += 1
            end
        end
    end
    return minimum([distance, currentmin])
end

function multi_dist_shapelet_series(multishapelet, series; LocInv::LocInvMethod=Yes())
    #calculate the distance between a multivariate shapelet and a multivariate series
    ℓ = multishapelet[1][end,1]   
    m = series[1][end,1]    
    V = Array{Float64}(undef, 0)  #the set containing shapelet start times for which a change in the shapelet will coincide with a change in the series. \mathcal{V} in the theorems 
    for d in eachindex(multishapelet)
        for i in multishapelet[d][:,1]
            for j in series[d][:,1]
                if i <= j && ℓ-i+j <= m 
                    V = vcat(V, j-i)
                end
            end
        end
    end
    V = unique(sort(V))
    distance = 1000000.0

    for j in V
        multiseries_segment = Vector{Matrix{Float64}}()
        for d in eachindex(series)
            series_start = argmin(series[d][:, 1] .<= j) -1
            series_end = argmin(series[d][:, 1] .< j+ℓ) - 1
            series_segment = series[d][series_start:series_end,:]
            series_segment[:,1] .-= j
            series_segment[1,1] = 0.0
            series_segment = vcat(series_segment, [ℓ series_segment[end,2]])
            push!(multiseries_segment, series_segment)
        end
        testdistance = multi_update_distance(multishapelet, multiseries_segment, distance, LocInv=LocInv)
        distance = testdistance
    end
    return round(distance, digits=8) 
end

function MultiTestShapelet(multishapelet, class1series, class2series, current_best_gain; LocInv::LocInvMethod=Yes())
    #multivariate version of TestShapelet()
    best_possible = 100000.0
    
    n₁ = length(class1series[1])
    n₂ = length(class2series[1])
    n = n₁ + n₂

    class1dists = multi_dist_shapelet_series(multishapelet, [class1series[d][1] for d in eachindex(class1series)], LocInv=LocInv)#[1]
    class2dists = multi_dist_shapelet_series(multishapelet, [class2series[d][1] for d in eachindex(class2series)], LocInv=LocInv)#[1]

    class1_counter = 1   #counts how many class 1 dists have been calculated
    class2_counter = 1   #counts how many class 2 dists have been calculated

    while best_possible >= current_best_gain && (class1_counter + class2_counter) < n
        if class1_counter < n₁
            new_class1dist = multi_dist_shapelet_series(multishapelet, [class1series[d][class1_counter+1] for d in eachindex(class1series)], LocInv=LocInv)#[1]
            class1dists = vcat(class1dists, new_class1dist)
            class1_counter += 1
        end
        if class2_counter < n₂
            new_class2dist = multi_dist_shapelet_series(multishapelet, [class2series[d][class2_counter+1] for d in eachindex(class2series)], LocInv=LocInv)#[1]
            class2dists = vcat(class2dists, new_class2dist)
            class2_counter += 1
        end

        best_possible = bestPossible(class1dists, class2dists, n₁, n₂)
        
    end
    return best_possible, class1dists, class2dists
    #If all distances were calculated, then best_possible is the actual infortmation gain
end

function MultiFindShapelet(class1series, class2series, ℓ, lrange, τ; LocInv::LocInvMethod=Yes(), Search::SearchOver=All())  #current_gain=0.0 
    #ℓ is the target length of shapelets to look for.
    #lrange is the full set of lengths to consider (\mathcal{I} in paper)
    #class1series and class2series are length d (dimension) vectors of length n (number of series) vectors of matrices - for each matrix, column 1 contains the change times and column 2 the values
    #τ is the time steps for shapelet extractions
    #LocInv = Yes() uses location invariant distance function (default), otherwise use LocInv = No()
    #if Search = Class1() or Class2(), we find the optimal characteristic shapelet of class 1 or 2, respectively. If Search = All(), we find the overall optimal shapelet from either class.
    
    if Search isa All
        searchseries = vcat(class1series[d], class2series[d])  #needs fixing. 
    elseif Search isa Class1
        searchseries = class1series
    elseif Search isa Class2
        searchseries = class2series
    end
    n = length(searchseries[1])
    m = searchseries[1][1][end,1]  #find time series length (assume all the same)
    best_gain = 0.0 
    shapelet_meandist = 0.0  
    shapelet_threshold = 0.0
    best_shapelet = searchseries[1][1]
    perfect_gain = OptimalGain(zeros(n), ones(n))
    @showprogress "Shapelets Search: " for i in 1:n
        for j in 0:τ:(m-ℓ)
            multishapelet = Vector{Matrix{Float64}}()
            for d in eachindex(searchseries)
                shapelet_start = argmin(searchseries[d][i][:, 1] .<= j) -1
                shapelet_end = argmin(searchseries[d][i][:, 1] .< j+ℓ) -1
                shapelet = searchseries[d][i][shapelet_start:shapelet_end,:]
                shapelet[1,1] = j
                shapelet[:,1] .-= j
                shapelet = vcat(shapelet, [ℓ shapelet[end,2]])
                push!(multishapelet, shapelet)
            end
            test = MultiTestShapelet(multishapelet, class1series, class2series, best_gain; LocInv=LocInv)
            test_gain = test[1]

            if test_gain >= best_gain  #if the candidate is at least as good as the current best shapelet, move to a local shapelet search
                if Search isa Class1 || Search isa Class2  #if looking specifically for characteristic shapelets of class 1 or class 2, check the subset proportions to validate this
                    Optsplit = OptimalSplit(test[2], test[3])
                    γ_star = Optsplit[2]
                    subset_split = vcat(test[2], test[3]) .<= γ_star
                    class1_split = test[2] .<= γ_star
                    class2_split = test[3] .<= γ_star
                    prop1near = sum(class1_split)/sum(subset_split)
                    prop1far = (length(class1_split) - sum(class1_split))/(length(subset_split) - sum(subset_split))
                    prop2near = sum(class2_split)/sum(subset_split)
                    prop2far = (length(class2_split) - sum(class2_split))/(length(subset_split) - sum(subset_split))
                    if (Search isa Class1 && prop1near > prop1far) || (Search isa Class2 && prop2near > prop2far)
                        localsearch = MultiLocalShapeletSearch(class1series, class2series, lrange, τ, LocInv=LocInv, series=[searchseries[d][i] for d in eachindex(searchseries)], position=j, current_gain=best_gain)
                        best_shapelet = localsearch[1]
                        best_gain = localsearch[2]
                        shapelet_threshold = localsearch[3]
                        shapelet_meandist = localsearch[4]
                    end
                else  #if Search isa All
                    localsearch = MultiLocalShapeletSearch(class1series, class2series, lrange, τ, LocInv=LocInv, series=[searchseries[d][i] for d in eachindex(searchseries)], position=j, current_gain=best_gain)
                    best_shapelet = localsearch[1]
                    best_gain = localsearch[2]
                    shapelet_threshold = localsearch[3]
                    shapelet_meandist = localsearch[4]
                end
            end

            if best_gain == perfect_gain
                return best_shapelet, best_gain, shapelet_threshold, shapelet_meandist
            end
        end
    end
    return best_shapelet, best_gain, shapelet_threshold, shapelet_meandist#, shapelet_series, shapelet_positions
end

function MultiLocalShapeletSearch(class1series, class2series, lrange, τ; LocInv::LocInvMethod=Yes(), series, position, current_gain)
    #For use within MultiFindShapelet()
    m = series[1][end,1]  #find time series length (assume all the same)
    best_gain = current_gain  
    shapelet_meandist = 0.0  
    shapelet_threshold = 0.0
    best_shapelet = series
    perfect_gain = OptimalGain(zeros(length(class1series[1])), ones(length(class2series[1])))
    for ℓ in lrange
        for j in maximum([0,position-1]):(τ/50):minimum([m-ℓ,position+1]) 
            ####extract candidate shapelet####
            multishapelet = Vector{Matrix{Float64}}()
            for d in eachindex(series)
                shapelet_start = argmin(series[d][:, 1] .<= j) -1
                shapelet_end = argmin(series[d][:, 1] .< j+ℓ) -1
                shapelet = series[d][shapelet_start:shapelet_end,:]
                shapelet[1,1] = j
                shapelet[:,1] .-= j
                shapelet = vcat(shapelet, [ℓ shapelet[end,2]])
                push!(multishapelet, shapelet)
            end
            ####test candidate shapelet with admissible entropy pruning####
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
                if test_split[3] > shapelet_meandist    ##break ties by maximum class separation.. test_meandist > shapelet_meandist
                    best_shapelet = multishapelet
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