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

