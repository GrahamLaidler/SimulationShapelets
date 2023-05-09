include("src/functionswrapper.jl")

#load data
system1series = load_object("dat/Multivariate/system1series.jld2");
system2series = load_object("dat/Multivariate/system2series.jld2");

#create training and testing sets
system1series_train = [system1series[i][1:100] for i in eachindex(system1series)];
system2series_train = [system2series[i][1:100] for i in eachindex(system2series)];
system1series_test = [system1series[i][101:1100] for i in eachindex(system1series)];
system2series_test = [system2series[i][101:1100] for i in eachindex(system2series)];


####Find optimal univariate shapelets in each dimension####
IndepA_s₁ = FindShapelet(system1series_train[1], system2series_train[1], 11, 10:12, 2, LocInv=Yes(), Search=Class1())
IndepA_s₂ = FindShapelet(system1series_train[1], system2series_train[1], 11, 10:12, 2, LocInv=Yes(), Search=Class2())

IndepB_s₁ = FindShapelet(system1series_train[2], system2series_train[2], 11, 10:12, 2, LocInv=Yes(), Search=Class1())
IndepB_s₂ = FindShapelet(system1series_train[2], system2series_train[2], 11, 10:12, 2, LocInv=Yes(), Search=Class2())

IndepC_s₁ = FindShapelet(system1series_train[3], system2series_train[3], 11, 10:12, 2, LocInv=Yes(), Search=Class1())
IndepC_s₂ = FindShapelet(system1series_train[3], system2series_train[3], 11, 10:12, 2, LocInv=Yes(), Search=Class2())

save_object("res/Multivariate/IndepA_s1.jld2", IndepA_s₁)
save_object("res/Multivariate/IndepA_s2.jld2", IndepA_s₂)
save_object("res/Multivariate/IndepB_s1.jld2", IndepB_s₁)
save_object("res/Multivariate/IndepB_s2.jld2", IndepB_s₂)
save_object("res/Multivariate/IndepC_s1.jld2", IndepC_s₁)
save_object("res/Multivariate/IndepC_s2.jld2", IndepC_s₂)

#Calculate distances from the optimal shapelets to the training and testing trajectories for scatter plots
testdistsA_s₁_system1 = [dist_shapelet_series(IndepA_s₁[1], system1series_test[1][i], LocInv=Yes()) for i in eachindex(system1series_test[1])];
testdistsA_s₁_system2 = [dist_shapelet_series(IndepA_s₁[1], system2series_test[1][i], LocInv=Yes()) for i in eachindex(system2series_test[1])];
testdistsA_s₂_system1 = [dist_shapelet_series(IndepA_s₂[1], system1series_test[1][i], LocInv=Yes()) for i in eachindex(system1series_test[1])];
testdistsA_s₂_system2 = [dist_shapelet_series(IndepA_s₂[1], system2series_test[1][i], LocInv=Yes()) for i in eachindex(system2series_test[1])];
testdistsB_s₁_system1 = [dist_shapelet_series(IndepB_s₁[1], system1series_test[2][i], LocInv=Yes()) for i in eachindex(system1series_test[2])];
testdistsB_s₁_system2 = [dist_shapelet_series(IndepB_s₁[1], system2series_test[2][i], LocInv=Yes()) for i in eachindex(system2series_test[2])];
testdistsB_s₂_system1 = [dist_shapelet_series(IndepB_s₂[1], system1series_test[2][i], LocInv=Yes()) for i in eachindex(system1series_test[2])];
testdistsB_s₂_system2 = [dist_shapelet_series(IndepB_s₂[1], system2series_test[2][i], LocInv=Yes()) for i in eachindex(system2series_test[2])];
testdistsC_s₁_system1 = [dist_shapelet_series(IndepC_s₁[1], system1series_test[3][i], LocInv=Yes()) for i in eachindex(system1series_test[3])];
testdistsC_s₁_system2 = [dist_shapelet_series(IndepC_s₁[1], system2series_test[3][i], LocInv=Yes()) for i in eachindex(system2series_test[3])];
testdistsC_s₂_system1 = [dist_shapelet_series(IndepC_s₂[1], system1series_test[3][i], LocInv=Yes()) for i in eachindex(system1series_test[3])];
testdistsC_s₂_system2 = [dist_shapelet_series(IndepC_s₂[1], system2series_test[3][i], LocInv=Yes()) for i in eachindex(system2series_test[3])];
save_object("res/Multivariate/testdistsA_s1_system1.jld2", testdistsA_s₁_system1)
save_object("res/Multivariate/testdistsA_s1_system2.jld2", testdistsA_s₁_system2)
save_object("res/Multivariate/testdistsA_s2_system1.jld2", testdistsA_s₂_system1)
save_object("res/Multivariate/testdistsA_s2_system2.jld2", testdistsA_s₂_system2)
save_object("res/Multivariate/testdistsB_s1_system1.jld2", testdistsB_s₁_system1)
save_object("res/Multivariate/testdistsB_s1_system2.jld2", testdistsB_s₁_system2)
save_object("res/Multivariate/testdistsB_s2_system1.jld2", testdistsB_s₂_system1)
save_object("res/Multivariate/testdistsB_s2_system2.jld2", testdistsB_s₂_system2)
save_object("res/Multivariate/testdistsC_s1_system1.jld2", testdistsC_s₁_system1)
save_object("res/Multivariate/testdistsC_s1_system2.jld2", testdistsC_s₁_system2)
save_object("res/Multivariate/testdistsC_s2_system1.jld2", testdistsC_s₂_system1)
save_object("res/Multivariate/testdistsC_s2_system2.jld2", testdistsC_s₂_system2)


####Find optimal multivariate shapelets####
Multi_s₁ = MultiFindShapelet(system1series_train, system2series_train, 11, 10:12, 10, LocInv=Yes(), Search=Class1())
Multi_s₂ = MultiFindShapelet(system1series_train, system2series_train, 11, 10:12, 10, LocInv=Yes(), Search=Class2())
save_object("res/Multivariate/Multi_s1.jld2", Multi_s₁)
save_object("res/Multivariate/Multi_s2.jld2", Multi_s₂)

#Calculate distances from the optimal shapelets to the training and testing trajectories for scatter plots
testdistsmulti_s₁_system1 = [multi_dist_shapelet_series(Multi_s₁[1], [system1series_test[d][i] for d in eachindex(system1series_test)], LocInv=Yes()) for i in eachindex(system1series_test[1])]
testdistsmulti_s₁_system2 = [multi_dist_shapelet_series(Multi_s₁[1], [system2series_test[d][i] for d in eachindex(system2series_test)], LocInv=Yes()) for i in eachindex(system2series_test[1])]
testdistsmulti_s₂_system1 = [multi_dist_shapelet_series(Multi_s₂[1], [system1series_test[d][i] for d in eachindex(system1series_test)], LocInv=Yes()) for i in eachindex(system1series_test[1])]
testdistsmulti_s₂_system2 = [multi_dist_shapelet_series(Multi_s₂[1], [system2series_test[d][i] for d in eachindex(system2series_test)], LocInv=Yes()) for i in eachindex(system2series_test[1])]
save_object("res/Multivariate/testdistsmulti_s1_system1.jld2", testdistsmulti_s₁_system1)
save_object("res/Multivariate/testdistsmulti_s1_system2.jld2", testdistsmulti_s₁_system2)
save_object("res/Multivariate/testdistsmulti_s2_system1.jld2", testdistsmulti_s₂_system1)
save_object("res/Multivariate/testdistsmulti_s2_system2.jld2", testdistsmulti_s₂_system2)