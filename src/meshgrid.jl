# copied and pasted directly from
# https://github.com/ChrisRackauckas/VectorizedRoutines.jl/blob/master/src/matlab.jl
#
# using VectorizedRoutines.jl directly causes Pkg versioning issues with SpecialFunctions.jl
# Not sure why, so I just put the meshgrid code directly in here.

"""
meshgrid(vx)
Computes an (x,y)-grid from the vectors (vx,vx).
For more information, see the MATLAB documentation.

Copied and pasted directly from [VectorizedRoutines.jl](https://github.com/ChrisRackauckas/VectorizedRoutines.jl/blob/master/src/matlab.jl).
Using VectorizedRoutines.jl directly causes Pkg versioning issues with SpecialFunctions.jl
"""
meshgrid(v::AbstractVector) = meshgrid(v, v)

"""
meshgrid(vx,vy)
Computes an (x,y)-grid from the vectors (vx,vy).
For more information, see the MATLAB documentation.

Copied and pasted directly from [VectorizedRoutines.jl](https://github.com/ChrisRackauckas/VectorizedRoutines.jl/blob/master/src/matlab.jl).
Using VectorizedRoutines.jl directly causes Pkg versioning issues with SpecialFunctions.jl
"""
function meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T}) where {T}
    m, n = length(vy), length(vx)
    return repeat(reshape(vx, 1, n), m, 1), repeat(reshape(vy, m, 1), 1, n)
end

"""
meshgrid(vx,vy,vz)
Computes an (x,y,z)-grid from the vectors (vx,vy,vz).
For more information, see the MATLAB documentation.

Copied and pasted directly from [VectorizedRoutines.jl](https://github.com/ChrisRackauckas/VectorizedRoutines.jl/blob/master/src/matlab.jl).
Using VectorizedRoutines.jl directly causes Pkg versioning issues with SpecialFunctions.jl
"""
function meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T},
                  vz::AbstractVector{T}) where {T}
    m, n, o = length(vy), length(vx), length(vz)
    vx = reshape(vx, 1, n, 1)
    vy = reshape(vy, m, 1, 1)
    vz = reshape(vz, 1, 1, o)
    om = ones(Int, m)
    on = ones(Int, n)
    oo = ones(Int, o)
    (vx[om, :, oo], vy[:, on, oo], vz[om, on, :])
end
