include("src/functionswrapper.jl")

#########Experiment with controlled breakdowns#############
seriesLRWcontrolled = load_object("dat/Breakdowns/LRW_controlled.jld2");
seriesMRWcontrolled = load_object("dat/Breakdowns/MRW_controlled.jld2");

#claculate average throughput and confidence intervals for LRW and MRW systems
LRWmeanTP = mean([seriesLRWcontrolled[i][end,2] for i in 1:length(seriesLRWcontrolled)])
MRWmeanTP = mean([seriesMRWcontrolled[i][end,2] for i in 1:length(seriesMRWcontrolled)])

LRWsd = std([seriesLRWcontrolled[i][end,2] for i in 1:length(seriesLRWcontrolled)]);
MRWsd = std([seriesMRWcontrolled[i][end,2] for i in 1:length(seriesMRWcontrolled)]);
LRWmeanTP95CI = (LRWmeanTP-(1.96*LRWsd/sqrt(length(seriesLRWcontrolled))), LRWmeanTP+(1.96*LRWsd/sqrt(length(seriesLRWcontrolled))))
MRWmeanTP95CI = (MRWmeanTP-(1.96*MRWsd/sqrt(length(seriesMRWcontrolled))), MRWmeanTP+(1.96*MRWsd/sqrt(length(seriesMRWcontrolled))))

#convert time units from minutes to hours
for i in 1:length(seriesLRWcontrolled)
    seriesLRWcontrolled[i][:,1] ./= 60
    seriesMRWcontrolled[i][:,1] ./= 60
end

#create training and testing sets
seriesLRW_train = seriesLRWcontrolled[1:100];
seriesMRW_train = seriesMRWcontrolled[1:100];
seriesLRW_test = seriesLRWcontrolled[101:end];
seriesMRW_test = seriesMRWcontrolled[101:end];

#Find optimal shapelets
LRWcontrolled_s₁ = FindShapelet(seriesLRW_train, seriesMRW_train, 45, 40:5:50, 2, shapeletrange = (26,34), Center=Opt(), Search=Class1())
MRWcontrolled_s₂ = FindShapelet(seriesLRW_train, seriesMRW_train, 45, 40:5:50, 2, shapeletrange = (24,32), Center=Opt(), Search=Class2())
save_object("res/Breakdowns/LRWcontrolled_s1.jld2", LRWcontrolled_s₁)
save_object("res/Breakdowns/MRWcontrolled_s2.jld2", MRWcontrolled_s₂)

traindists_s₁_LRW = [dist_shapelet_series(LRWcontrolled_s₁[1], seriesLRW_train[i], Center=Opt()) for i in 1:length(seriesLRW_train)];
traindists_s₁_MRW = [dist_shapelet_series(LRWcontrolled_s₁[1], seriesMRW_train[i], Center=Opt()) for i in 1:length(seriesMRW_train)];
traindists_s₂_LRW = [dist_shapelet_series(MRWcontrolled_s₂[1], seriesLRW_train[i], Center=Opt()) for i in 1:length(seriesLRW_train)];
traindists_s₂_MRW = [dist_shapelet_series(MRWcontrolled_s₂[1], seriesMRW_train[i], Center=Opt()) for i in 1:length(seriesMRW_train)];
testdists_s₁_LRW = [dist_shapelet_series(LRWcontrolled_s₁[1], seriesLRW_test[i], Center=Opt()) for i in 1:length(seriesLRW_test)];
testdists_s₁_MRW = [dist_shapelet_series(LRWcontrolled_s₁[1], seriesMRW_test[i], Center=Opt()) for i in 1:length(seriesMRW_test)];
testdists_s₂_LRW = [dist_shapelet_series(MRWcontrolled_s₂[1], seriesLRW_test[i], Center=Opt()) for i in 1:length(seriesLRW_test)];
testdists_s₂_MRW = [dist_shapelet_series(MRWcontrolled_s₂[1], seriesMRW_test[i], Center=Opt()) for i in 1:length(seriesMRW_test)];
save_object("res/Breakdowns/controlled_traindists_s1_LRW.jld2", traindists_s₁_LRW)
save_object("res/Breakdowns/controlled_traindists_s1_MRW.jld2", traindists_s₁_MRW)
save_object("res/Breakdowns/controlled_traindists_s2_LRW.jld2", traindists_s₂_LRW)
save_object("res/Breakdowns/controlled_traindists_s2_MRW.jld2", traindists_s₂_MRW)
save_object("res/Breakdowns/controlled_testdists_s1_LRW.jld2", testdists_s₁_LRW)
save_object("res/Breakdowns/controlled_testdists_s1_MRW.jld2", testdists_s₁_MRW)
save_object("res/Breakdowns/controlled_testdists_s2_LRW.jld2", testdists_s₂_LRW)
save_object("res/Breakdowns/controlled_testdists_s2_MRW.jld2", testdists_s₂_MRW)



##########Experiment with random breakdowns##########
seriesLRWrandom = load_object("dat/Breakdowns/LRW_random.jld2");
seriesMRWrandom = load_object("dat/Breakdowns/MRW_random.jld2");


#convert time units from minutes to hours
for i in 1:length(seriesLRWrandom)
    seriesLRWrandom[i][:,1] ./= 60
    seriesMRWrandom[i][:,1] ./= 60
end

#create training and testing sets
seriesLRWrandom_train = seriesLRWrandom[1:100];
seriesMRWrandom_train = seriesMRWrandom[1:100];
seriesLRWrandom_test = seriesLRWrandom[101:end];
seriesMRWrandom_test = seriesMRWrandom[101:end];

#Find optimal shapelets
LRWrandom_s₁ = FindShapelet(seriesLRWrandom_train, seriesMRWrandom_train, 80, 70:5:90, 10, Center=Opt(), Search=Class1())
MRWrandom_s₂ = FindShapelet(seriesLRWrandom_train, seriesMRWrandom_train, 80, 70:5:90, 10, Center=Opt(), Search=Class2())
save_object("res/Breakdowns/LRWrandom_s1.jld2", LRWrandom_s₁)
save_object("res/Breakdowns/MRWrandom_s2.jld2", MRWrandom_s₂)


r_traindists_s₁_LRW = [dist_shapelet_series(LRWrandom_s₁[1], seriesLRWrandom_train[i], Center=Opt()) for i in 1:length(seriesLRWrandom_train)];
r_traindists_s₁_MRW = [dist_shapelet_series(LRWrandom_s₁[1], seriesMRWrandom_train[i], Center=Opt()) for i in 1:length(seriesMRWrandom_train)];
r_traindists_s₂_LRW = [dist_shapelet_series(MRWrandom_s₂[1], seriesLRWrandom_train[i], Center=Opt()) for i in 1:length(seriesLRWrandom_train)];
r_traindists_s₂_MRW = [dist_shapelet_series(MRWrandom_s₂[1], seriesMRWrandom_train[i], Center=Opt()) for i in 1:length(seriesMRWrandom_train)];
r_testdists_s₁_LRW = [dist_shapelet_series(LRWrandom_s₁[1], seriesLRWrandom_test[i], Center=Opt()) for i in 1:length(seriesLRWrandom_test)];
r_testdists_s₁_MRW = [dist_shapelet_series(LRWrandom_s₁[1], seriesMRWrandom_test[i], Center=Opt()) for i in 1:length(seriesMRWrandom_test)];
r_testdists_s₂_LRW = [dist_shapelet_series(MRWrandom_s₂[1], seriesLRWrandom_test[i], Center=Opt()) for i in 1:length(seriesLRWrandom_test)];
r_testdists_s₂_MRW = [dist_shapelet_series(MRWrandom_s₂[1], seriesMRWrandom_test[i], Center=Opt()) for i in 1:length(seriesMRWrandom_test)];
save_object("res/Breakdowns/random_traindists_s1_LRW.jld2", r_traindists_s₁_LRW)
save_object("res/Breakdowns/random_traindists_s1_MRW.jld2", r_traindists_s₁_MRW)
save_object("res/Breakdowns/random_traindists_s2_LRW.jld2", r_traindists_s₂_LRW)
save_object("res/Breakdowns/random_traindists_s2_MRW.jld2", r_traindists_s₂_MRW)
save_object("res/Breakdowns/random_testdists_s1_LRW.jld2", r_testdists_s₁_LRW)
save_object("res/Breakdowns/random_testdists_s1_MRW.jld2", r_testdists_s₁_MRW)
save_object("res/Breakdowns/random_testdists_s2_LRW.jld2", r_testdists_s₂_LRW)
save_object("res/Breakdowns/random_testdists_s2_MRW.jld2", r_testdists_s₂_MRW)