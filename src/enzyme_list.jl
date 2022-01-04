using DataFrames, CSV

dir = @__DIR__

home_dir = joinpath(split(dir, '/')[1:end-1]...)
enzyme_list = CSV.read("/$home_dir/data/re_enzyme_list.csv", DataFrame)
enzyme_list[!,:enzyme] = convert.(String, enzyme_list[!,:enzyme] )
