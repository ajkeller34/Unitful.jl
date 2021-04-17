"""
    promote_unit(::Units, ::Units...)
Given `Units` objects as arguments, this function returns a `Units` object appropriate
for the result of promoting quantities which have these units. This function is kind
of like `promote_rule`, except that it doesn't take types. It also does not return a tuple,
but rather just a [`Unitful.Units`](@ref) object (or it throws an error).

Although we had used `promote_rule` for `Units` objects in prior versions of Unitful,
this was always kind of a hack; it doesn't make sense to promote units directly for
a variety of reasons.
"""
function promote_unit end

# Generic methods
@inline promote_unit(x) = _promote_unit(x)
@inline _promote_unit(x::Units) = x

@inline promote_unit(x,y) = _promote_unit(x,y)

promote_unit(x::Units, y::Units, z::Units, t::Units...) =
    promote_unit(_promote_unit(x,y), z, t...)

@inline _promote_unit(x::T, y::T) where {T <: FreeUnits} = T()
# Use configurable fall-back mechanism for FreeUnits
@inline _promote_unit(x::FreeUnits{N1,D}, y::FreeUnits{N2,D}) where {N1,N2,D} =
    upreferred(dimension(x))

@inline _promote_unit(x::ContextUnits{N,D,P,A}, y::ContextUnits{N,D,P,A}) where {N,D,P,A} =
    ContextUnits{N,D,P,A}() # ambiguity reasons
# same units, but promotion context disagrees
@inline _promote_unit(x::ContextUnits{N,D,P1,A}, y::ContextUnits{N,D,P2,A}) where {N,D,P1,P2,A} =
    ContextUnits{N,D,promote_unit(P1(), P2()),A}()
# different units, but promotion context agrees
@inline _promote_unit(x::ContextUnits{N1,D,P}, y::ContextUnits{N2,D,P}) where {N1,N2,D,P} =
    ContextUnits(P(), P())
# different units, promotion context disagrees, fall back to FreeUnits
@inline _promote_unit(x::ContextUnits{N1,D,P1}, y::ContextUnits{N2,D,P2}) where {N1,N2,D,P1,P2} =
    promote_unit(FreeUnits(x), FreeUnits(y))

# ContextUnits beat FreeUnits
@inline _promote_unit(x::ContextUnits{N,D,P,A}, y::FreeUnits{N,D,A}) where {N,D,P,A} = x
@inline _promote_unit(x::ContextUnits{N,D,P,A1}, y::FreeUnits{N,D,A2}) where {N,D,P,A1,A2} =
    ContextUnits(P(), P())
@inline _promote_unit(x::ContextUnits{N1,D,P}, y::FreeUnits{N2,D}) where {N1,N2,D,P} =
    ContextUnits(P(), P())
@inline _promote_unit(x::FreeUnits, y::ContextUnits) = promote_unit(y,x)

# FixedUnits beat everything
@inline _promote_unit(x::T, y::T) where {T <: FixedUnits} = T()
@inline _promote_unit(x::FixedUnits{M,D}, y::Units{N,D}) where {M,N,D} = x
@inline _promote_unit(x::Units, y::FixedUnits) = promote_unit(y,x)

# Different units but same dimension are not fungible for FixedUnits
@inline _promote_unit(x::FixedUnits{M,D}, y::FixedUnits{N,D}) where {M,N,D} =
    error("automatic conversion prohibited.")

# If we didn't handle it above, the dimensions mismatched.
@inline _promote_unit(x::Units, y::Units) = throw(DimensionError(x,y))

####
#    Base.promote_rule

# quantity, quantity (different dims)
Base.promote_rule(::Type{Quantity{S1,D1,U1}},
        ::Type{Quantity{S2,D2,U2}}) where {S1,D1,U1,S2,D2,U2} =
    Quantity{promote_type(S1,S2)}

# quantity, quantity (same dims, different units)
function Base.promote_rule(::Type{Quantity{S1,D,U1}},
        ::Type{Quantity{S2,D,U2}}) where {S1,S2,D,U1,U2}

    p = promote_unit(U1(), U2())
    c1 = convfact(p, U1())
    c1′ = affinetranslation(U1())
    c2 = convfact(p, U2())
    c2′ = affinetranslation(U2())
    numtype = promote_type(S1, S2,
        promote_type(typeof(c1), typeof(c2), typeof(c1′), typeof(c2′)))
    if !isunitless(p)
        if U1 <: ContextUnits && U2 <: ContextUnits
            up1 = upreferred(U1())
            if up1 === upreferred(U2())
                return Quantity{numtype,D,typeof(ContextUnits(p,up1))}
            else
                return Quantity{numtype,D,typeof(p)}
            end
        elseif U1 <: ContextUnits || U2 <: ContextUnits
            return Quantity{numtype,D,typeof(ContextUnits(p,p))}
        else
            return Quantity{numtype,D,typeof(p)}
        end
    else
        return numtype
    end
end

# number, quantity
function Base.promote_rule(::Type{Quantity{S,D,U}}, ::Type{T}) where {S,T <: Number,D,U}
    if D == NoDims
        promote_type(S,T,typeof(convfact(NoUnits,U())))
    else
        Quantity{promote_type(S,T)}
    end
end

Base.promote_rule(::Type{S}, ::Type{T}) where {S<:AbstractIrrational,S2,T<:Quantity{S2}} =
    promote_type(promote_type(S, real(S2)), T)

Base.promote_rule(::Type{Quantity{S}}, ::Type{T}) where {S,T <: Number} =
    Quantity{promote_type(S,T)}

# With only one of these, you can get a segmentation fault because you # fall back to the
# number, quantity promote_rule above and there is an infinite recursion.
Base.promote_rule(::Type{Quantity{T}}, ::Type{Quantity{S,D,U}}) where {T,S,D,U} =
    Quantity{promote_type(T,S)}

Base.promote_rule(::Type{Quantity{S,D,U}}, ::Type{Quantity{T}}) where {T,S,D,U} =
    Quantity{promote_type(T,S)}

Base.promote_typejoin(::Type{T}, ::Type{Quantity{T,D,U}}) where {T,D,U} = Quantity{T}
