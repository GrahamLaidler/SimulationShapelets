include("src/functionswrapper.jl")


trueseries = load_object("dat/Validation/trueseries.jld2");
model1series = load_object("dat/Validation/model1series.jld2");
model2series = load_object("dat/Validation/model2series.jld2");




#create training and testing sets
trueseries_train = trueseries[1:200];
trueseries_test = trueseries[201:5200];
model1series_train = model1series[1:200];
model1series_test = model1series[201:5200];
model2series_train = model2series[1:200];
model2series_test = model2series[201:5200];

Model1_s₁ = FindShapelet(trueseries_train, model1series_train, 10, 8:12, 2, Center=Opt(), Search=Class1())
Model1_s₂ = FindShapelet(trueseries_train, model1series_train, 10, 8:12, 2, Center=Opt(), Search=Class2())
Model2_s₁ = FindShapelet(trueseries_train, model2series_train, 10, 8:12, 2, Center=Opt(), Search=Class1())
Model2_s₂ = FindShapelet(trueseries_train, model2series_train, 10, 8:12, 2, Center=Opt(), Search=Class2())
save_object("res/Validation/Model1_s1.jld2", Model1_s₁)
save_object("res/Validation/Model1_s2.jld2", Model1_s₂)
save_object("res/Validation/Model2_s1.jld2", Model2_s₁)
save_object("res/Validation/Model2_s2.jld2", Model2_s₂)


m1_traindists_s₁_true = [dist_shapelet_series(Model1_s₁[1], trueseries_train[i], Center=Opt()) for i in 1:length(trueseries_train)];
m1_traindists_s₁_model = [dist_shapelet_series(Model1_s₁[1], model1series_train[i], Center=Opt()) for i in 1:length(model1series_train)];
m1_traindists_s₂_true = [dist_shapelet_series(Model1_s₂[1], trueseries_train[i], Center=Opt()) for i in 1:length(trueseries_train)];
m1_traindists_s₂_model = [dist_shapelet_series(Model1_s₂[1], model1series_train[i], Center=Opt()) for i in 1:length(model1series_train)];
m1_testdists_s₁_true = [dist_shapelet_series(Model1_s₁[1], trueseries_test[i], Center=Opt()) for i in 1:length(trueseries_test)];
m1_testdists_s₁_model = [dist_shapelet_series(Model1_s₁[1], model1series_test[i], Center=Opt()) for i in 1:length(model1series_test)];
m1_testdists_s₂_true = [dist_shapelet_series(Model1_s₂[1], trueseries_test[i], Center=Opt()) for i in 1:length(trueseries_test)];
m1_testdists_s₂_model = [dist_shapelet_series(Model1_s₂[1], model1series_test[i], Center=Opt()) for i in 1:length(model1series_test)];
save_object("res/Validation/m1_traindists_s1_true.jld2", m1_traindists_s₁_true)
save_object("res/Validation/m1_traindists_s1_model.jld2", m1_traindists_s₁_model)
save_object("res/Validation/m1_traindists_s2_true.jld2", m1_traindists_s₂_true)
save_object("res/Validation/m1_traindists_s2_model.jld2", m1_traindists_s₂_model)
save_object("res/Validation/m1_testdists_s1_true.jld2", m1_testdists_s₁_true)
save_object("res/Validation/m1_testdists_s1_model.jld2", m1_testdists_s₁_model)
save_object("res/Validation/m1_testdists_s2_true.jld2", m1_testdists_s₂_true)
save_object("res/Validation/m1_testdists_s2_model.jld2", m1_testdists_s₂_model)

m2_traindists_s₁_true = [dist_shapelet_series(Model2_s₁[1], trueseries_train[i], Center=Opt()) for i in 1:length(trueseries_train)];
m2_traindists_s₁_model = [dist_shapelet_series(Model2_s₁[1], model2series_train[i], Center=Opt()) for i in 1:length(model2series_train)];
m2_traindists_s₂_true = [dist_shapelet_series(Model2_s₂[1], trueseries_train[i], Center=Opt()) for i in 1:length(trueseries_train)];
m2_traindists_s₂_model = [dist_shapelet_series(Model2_s₂[1], model2series_train[i], Center=Opt()) for i in 1:length(model2series_train)];
m2_testdists_s₁_true = [dist_shapelet_series(Model2_s₁[1], trueseries_test[i], Center=Opt()) for i in 1:length(trueseries_test)];
m2_testdists_s₁_model = [dist_shapelet_series(Model2_s₁[1], model2series_test[i], Center=Opt()) for i in 1:length(model2series_test)];
m2_testdists_s₂_true = [dist_shapelet_series(Model2_s₂[1], trueseries_test[i], Center=Opt()) for i in 1:length(trueseries_test)];
m2_testdists_s₂_model = [dist_shapelet_series(Model2_s₂[1], model2series_test[i], Center=Opt()) for i in 1:length(model2series_test)];
save_object("res/Validation/m2_traindists_s1_true.jld2", m2_traindists_s₁_true)
save_object("res/Validation/m2_traindists_s1_model.jld2", m2_traindists_s₁_model)
save_object("res/Validation/m2_traindists_s2_true.jld2", m2_traindists_s₂_true)
save_object("res/Validation/m2_traindists_s2_model.jld2", m2_traindists_s₂_model)
save_object("res/Validation/m2_testdists_s1_true.jld2", m2_testdists_s₁_true)
save_object("res/Validation/m2_testdists_s1_model.jld2", m2_testdists_s₁_model)
save_object("res/Validation/m2_testdists_s2_true.jld2", m2_testdists_s₂_true)
save_object("res/Validation/m2_testdists_s2_model.jld2", m2_testdists_s₂_model)