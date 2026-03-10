using StaticArrays

struct TestStruct{N}
    mat::SMatrix{N,N,Float64}
end

t = TestStruct{2}(@SMatrix [1.0 2.0; 3.0 4.0])
println("Success: ", t)

struct TestStructWithL{N}
    # This is expected to fail if N is a TypeVar
    mat::SMatrix{N,N,Float64,N * N}
end
