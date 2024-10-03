using TestItemRunner

@run_package_tests verbose = true

@testitem "Simple Test" begin
    
    @test 1 + 1 == 2

end