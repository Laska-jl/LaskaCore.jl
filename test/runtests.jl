using Test
using LaskaCore

include("Aqua.jl")

@testset "findmax" begin
    alleqspv = RelativeSpikeVector([ones(100) for _ in 1:10], 30)
    alleqv = [ones(100) for _ in 1:10]

    @inferred LaskaCore.minval(alleqspv)
    @test LaskaCore.minval(alleqspv) == 1.0
    @test LaskaCore.minval(alleqspv) isa eltype(alleqspv)

    @inferred LaskaCore.maxval(alleqspv)
    @test LaskaCore.maxval(alleqspv) == 1.0
    @test LaskaCore.maxval(alleqspv) isa eltype(alleqspv)


    @inferred LaskaCore.minval(alleqv)
    @test LaskaCore.minval(alleqv) == 1.0
    @test LaskaCore.minval(alleqv) isa eltype(eltype(alleqspv))

    @inferred LaskaCore.maxval(alleqv)
    @test LaskaCore.maxval(alleqv) == 1.0
    @test LaskaCore.maxval(alleqv) isa eltype(eltype(alleqspv))

end

@testset "isiviolations" begin
    testvec1 = [0, 2, 4, 6, 7, 8, 10, 11]
    @test LaskaCore.isiviolations(testvec1, 1) == [5, 6, 8]
    @test LaskaCore.isiviolations(testvec1, 2) == 2:length(testvec1) |> collect
    @test LaskaCore.isiviolations(testvec1, 0.5) == Int64[]

    testvecspike = SpikeVector(testvec1, 1)
    @test LaskaCore.isiviolations(testvecspike, 1) == [5, 6, 8]
    @test LaskaCore.isiviolations(testvecspike, 2) == 2:length(testvec1) |> collect
    @test LaskaCore.isiviolations(testvecspike, 0.5) == Int64[]

    testvecrel = RelativeSpikeVector([deepcopy(testvec1) for _ in 1:10], 1)
    @test all([LaskaCore.isiviolations(testvecrel, 1)[i] == [5, 6, 8] for i in 1:10])
    @test all([LaskaCore.isiviolations(testvecrel, 2)[i] == 2:length(testvec1) |> collect for i in 1:10])
    @test all([LaskaCore.isiviolations(testvecrel, 0.5)[i] == Int64[] for i in 1:10])
end
