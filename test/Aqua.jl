using Aqua

@testset "aqua deps compat" begin
    Aqua.test_deps_compat(LaskaCore)
end


@testset "aqua unbound_args" begin
    Aqua.test_unbound_args(LaskaCore)
end

@testset "aqua undefined exports" begin
    Aqua.test_undefined_exports(LaskaCore)
end


# Gives false positives from dependencies
#
# @testset "aqua test ambiguities" begin
#     Aqua.test_ambiguities([MyModule, Core, Base])
# end

@testset "aqua piracy" begin
    Aqua.test_piracies(LaskaCore)
end

@testset "aqua project extras" begin
    Aqua.test_project_extras(LaskaCore)
end

@testset "aqua state deps" begin
    Aqua.test_stale_deps(LaskaCore)
end
