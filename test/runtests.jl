using Test
using DiffusionGarnet
using GeoParams
using DelimitedFiles
using JLD2
using HDF5
using LinearAlgebra: norm

function runtests()
    files = readdir(@__DIR__)
    test_files = filter(startswith("test_"), files)

    for f in test_files
        if !isdir(f)
            include(f)
        end
    end
    return
end

runtests()