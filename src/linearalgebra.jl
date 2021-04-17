
# Multiplication

function mul!(C::StridedVecOrMat{<:AbstractQuantity{T}}, 
              A::StridedMatrix{<:AbstractQuantity{T}}, 
              B::StridedVecOrMat{<:AbstractQuantity{T}},
              alpha::Number, beta::Number) where {T<:Base.HWNumber}
    _mul!(C, A, B, alpha, beta)
end

function mul!(C::StridedVecOrMat{<:AbstractQuantity{T}}, 
              A::AdjOrTransAbsMat{<:AbstractQuantity{T}, <:StridedMatrix}, 
              B::StridedVecOrMat{<:AbstractQuantity{T}},
              alpha::Number, beta::Number) where {T<:Base.HWNumber}
    _mul!(C, A, B, alpha, beta)
end

function _mul!(C, A, B, alpha, beta)
    if unit(beta) != NoUnits
        throw(DimensionError("beta", 1.0))
    elseif unit(eltype(C)) != unit(eltype(A)) * unit(eltype(B)) * unit(alpha)
        throw(DimensionError("A * B .* α", "C"))
    end
    C0 = ustrip(C)
    A0 = ustrip(A)
    B0 = ustrip(B)
    mul!(C0, A0, B0)
    _linearalgebra_count()
    return C
end

function dot(A::StridedArray{<:AbstractQuantity{T}}, 
             B::StridedArray{<:AbstractQuantity{T}}) where {T<:Base.HWNumber}
    A0 = ustrip(A)
    B0 = ustrip(B)
    C0 = dot(A0, B0)
    _linearalgebra_count()
    C = C0 * unit(eltype(A)) * unit(eltype(B))
    return C
end

# Division

function (\)(A::StridedMatrix{<:AbstractQuantity{T}}, 
             B::StridedVecOrMat{<:AbstractQuantity{T}}) where {T<:Base.HWNumber}
    A0 = ustrip(A)
    B0 = ustrip(B)
    C0 = A0 \ B0
    _linearalgebra_count()
    u = unit(eltype(B)) / unit(eltype(A))
    Tu = typeof(one(eltype(C0)) * u)
    return reinterpret(Tu, C0)
end

function (/)(A::StridedVecOrMat{<:AbstractQuantity{T}}, 
             B::StridedVecOrMat{<:AbstractQuantity{T}}) where {T<:Base.HWNumber}
    A0 = ustrip(A)
    B0 = ustrip(B)
    C0 = A0 / B0
    _linearalgebra_count()
    u = unit(eltype(A)) / unit(eltype(B))
    Tu = typeof(one(eltype(C0)) * u)
    return reinterpret(Tu, C0)
end

function inv(A::StridedMatrix{<:AbstractQuantity{T}}) where {T<:Base.HWNumber}
    C0 = inv(ustrip(A))
    _linearalgebra_count()
    u = inv(unit(eltype(A)))
    Tu = typeof(one(eltype(C0)) * u)
    return reinterpret(Tu, C0)
end

function pinv(A::StridedMatrix{<:AbstractQuantity{T}}; kw...) where {T<:Base.HWNumber}
    C0 = pinv(ustrip(A); kw...)
    _linearalgebra_count()
    u = inv(unit(eltype(A)))
    Tu = typeof(one(eltype(C0)) * u)
    return reinterpret(Tu, C0)
end

# This function is re-defined during testing, to check we hit the fast path:
_linearalgebra_count() = nothing

