# This is a generated function to avoid determining the dimensions of a given
# set of units each time a new quantity is made.
@generated function _Quantity(x::Number, y::Units)
    u = y()
    du = dimension(u)
    dx = dimension(x)
    d = du*dx
    :(Quantity{typeof(x), $d, typeof($u)}(x))
end

"""
    Quantity(x::Number, y::Units)

Create a `Quantity` with numerical value `x` and units `y`.

# Example

```jldoctest
julia> Quantity(5, u"m")
5 m
```
"""
Quantity(x::Number, y::Units) = _Quantity(x, y)
Quantity(x::Number, y::Units{()}) = x

*(x::Number, y::Units, z::Units...) = Quantity(x,*(y,z...))
*(x::Units, y::Number) = *(y,x)

*(x::AbstractQuantity, y::Units, z::Units...) = Quantity(x.val, *(unit(x),y,z...))
*(x::AbstractQuantity, y::AbstractQuantity) = Quantity(x.val*y.val, unit(x)*unit(y))

# Next two lines resolves some method ambiguity:
*(x::Bool, y::T) where {T <: AbstractQuantity} =
    ifelse(x, y, ifelse(signbit(y), -zero(y), zero(y)))
*(x::AbstractQuantity, y::Bool) = Quantity(x.val*y, unit(x))

*(y::Number, x::AbstractQuantity) = *(x,y)
function *(x::AbstractQuantity, y::Number)
    x isa AffineQuantity &&
        throw(AffineError("an invalid operation was attempted with affine quantities: $x*$y"))
    return Quantity(x.val*y, unit(x))
end

*(A::Units, B::AbstractArray) = broadcast(*, A, B)
*(A::AbstractArray, B::Units) = broadcast(*, A, B)

/(A::AbstractArray, B::Units) = broadcast(/, A, B)

# Division (units)
/(x::AbstractQuantity, y::Units) = Quantity(x.val, unit(x) / y)
/(x::Units, y::AbstractQuantity) = Quantity(1/y.val, x / unit(y))
/(x::Number, y::Units) = Quantity(x,inv(y))
/(x::Units, y::Number) = (1/y) * x

//(x::AbstractQuantity, y::Units) = Quantity(x.val, unit(x) / y)
//(x::Units, y::AbstractQuantity) = Quantity(1//y.val, x / unit(y))
//(x::Number, y::Units) = Rational(x)/y
//(x::Units, y::Number) = (1//y) * x

/(x::AbstractQuantity, y::AbstractQuantity) = Quantity(/(x.val, y.val), unit(x) / unit(y))
/(x::AbstractQuantity, y::Number) = Quantity(/(x.val, y), unit(x) / unit(y))
/(x::Number, y::AbstractQuantity) = Quantity(/(x, y.val), unit(x) / unit(y))
//(x::AbstractQuantity, y::AbstractQuantity) = Quantity(//(x.val, y.val), unit(x) / unit(y))
//(x::AbstractQuantity, y::Number) = Quantity(//(x.val, y), unit(x) // unit(y))
//(x::Number, y::AbstractQuantity) = Quantity(//(x, y.val), unit(x) / unit(y))

# ambiguity resolution
//(x::AbstractQuantity, y::Complex) = Quantity(//(x.val, y), unit(x))

for f in (:fld, :cld)
    @eval begin
        function ($f)(x::AbstractQuantity, y::AbstractQuantity)
            z = uconvert(unit(y), x)        # TODO: use promote?
            ($f)(z.val,y.val)
        end

        ($f)(x::Number, y::AbstractQuantity) = Quantity(($f)(x, ustrip(y)), unit(x) / unit(y))

        ($f)(x::AbstractQuantity, y::Number) = Quantity(($f)(ustrip(x), y), unit(x))
    end
end

function div(x::AbstractQuantity, y::AbstractQuantity, r...)
    z = uconvert(unit(y), x)        # TODO: use promote?
    div(z.val,y.val, r...)
end

function div(x::Number, y::AbstractQuantity, r...)
    Quantity(div(x, ustrip(y), r...), unit(x) / unit(y))
end

function div(x::AbstractQuantity, y::Number, r...)
    Quantity(div(ustrip(x), y, r...), unit(x))
end

for f in (:mod, :rem)
    @eval function ($f)(x::AbstractQuantity, y::AbstractQuantity)
        z = uconvert(unit(y), x)        # TODO: use promote?
        Quantity(($f)(z.val,y.val), unit(y))
    end
end

_affineerror(f, args...) =
    throw(AffineError("an invalid operation was attempted with affine quantities: $f($(join(args, ", ")))"))

for f in (:div, :rem, :divrem)
    for r = (RoundNearest, RoundNearestTiesAway, RoundNearestTiesUp,
             RoundToZero, RoundUp, RoundDown)
        @eval begin
            $f(x::AffineQuantity, y::AffineQuantity, ::typeof($r)) = _affineerror($f, x, y, $r)
            $f(x::AffineQuantity, y::AbstractQuantity, ::typeof($r)) = _affineerror($f, x, y, $r)
            $f(x::AbstractQuantity, y::AffineQuantity, ::typeof($r)) = _affineerror($f, x, y, $r)
        end
    end
end
for f = (:div, :cld, :fld, :rem, :mod)
    @eval begin
        $f(x::AffineQuantity, y::AffineQuantity) = _affineerror($f, x, y)
        $f(x::AffineQuantity, y::AbstractQuantity) = _affineerror($f, x, y)
        $f(x::AbstractQuantity, y::AffineQuantity) = _affineerror($f, x, y)
    end
end

Base.mod2pi(x::DimensionlessQuantity) = mod2pi(uconvert(NoUnits, x))
Base.mod2pi(x::AbstractQuantity{S, NoDims, <:Units{(Unitful.Unit{:Degree, NoDims}(0, 1//1),),
    NoDims}}) where S = mod(x, 360°)

# Addition / subtraction
for op in [:+, :-]
    @eval ($op)(x::AbstractQuantity{S,D,U}, y::AbstractQuantity{T,D,U}) where {S,T,D,U} =
        Quantity(($op)(x.val, y.val), U())

    @eval function ($op)(x::AbstractQuantity{S,D,SU}, y::AbstractQuantity{T,D,TU}) where {S,T,D,SU,TU}
        ($op)(promote(x,y)...)
    end

    @eval ($op)(x::AbstractQuantity, y::AbstractQuantity) = throw(DimensionError(x,y))
    @eval ($op)(x::AbstractQuantity) = Quantity(($op)(x.val), unit(x))
end

function +(x::AffineQuantity{S,D}, y::AbstractQuantity{T,D}) where {S,T,D}
    pu = promote_unit(unit(x), unit(y))     # units for the final result.

    # Get x on an absolute scale. FreeUnits in the line below prevents
    # promote(x′, y) from yielding affine quantities. If x had `ContextUnits` and
    # the promotion units were affine units, x′+y would error without this.
    x′ = Quantity(x.val - affinetranslation(unit(x)), FreeUnits(absoluteunit(x)))

    # Likewise if y were not affine but y had ContextUnits and the promotion units were
    # affine, x′+y could also fail.
    y′ = Quantity(y.val, FreeUnits(unit(y)))

    return uconvert(pu, x′+y′)  # we get back the promotion context in the end
end
+(x::AbstractQuantity, y::AffineQuantity) = +(y,x)

# Disallow addition of affine quantities
+(x::AffineQuantity, y::AffineQuantity) = throw(AffineError(
   "an invalid operation was attempted with affine quantities: $x + $y"))

# Specialize subtraction of affine quantities
-(x::AffineQuantity, y::AffineQuantity) = -(promote(x,y)...)
function -(x::T, y::T) where T <: AffineQuantity
    return Quantity(x.val - y.val, absoluteunit(unit(x)))
end

# Disallow subtracting an affine quantity from a quantity
-(x::AbstractQuantity, y::AffineQuantity) =
    throw(AffineError("an invalid operation was attempted with affine quantities: $x - $y"))

# Needed until LU factorization is made to work with unitful numbers
function inv(x::StridedMatrix{T}) where {T <: AbstractQuantity}
    m = inv(ustrip(x))
    iq = eltype(m)
    reinterpret(Quantity{iq, inv(dimension(T)), typeof(inv(unit(T)))}, m)
end

for x in (:istriu, :istril)
    @eval ($x)(A::AbstractMatrix{T}) where {T <: AbstractQuantity} = ($x)(ustrip(A))
end

# Other mathematical functions

# `fma` and `muladd`
# The idea here is that if the numeric backing types are not the same, they
# will be promoted to be the same by the generic `fma(::Number, ::Number, ::Number)`
# method. We then catch the possible results and handle the units logic with one
# performant method.

for (_x,_y) in [(:fma, :_fma), (:muladd, :_muladd)]
    # Catch some signatures pre-promotion
    @eval @inline ($_x)(x::Number, y::AbstractQuantity, z::AbstractQuantity) = ($_y)(x,y,z)
    @eval @inline ($_x)(x::AbstractQuantity, y::Number, z::AbstractQuantity) = ($_y)(x,y,z)
    @eval @inline ($_x)(x::AbstractQuantity, y::AbstractQuantity, z::AbstractQuantity) = ($_y)(x,y,z)

    # Post-promotion
    @eval @inline ($_x)(x::AbstractQuantity{A}, y::AbstractQuantity{B}, z::AbstractQuantity{C}) where {
        A <: Number, B <: Number, C <: Number} = ($_y)(x,y,z)

    # It seems like most of this is optimized out by the compiler, including the
    # apparent runtime check of dimensions, which does not appear in @code_llvm.
    @eval @inline function ($_y)(x,y,z)
        dimension(x) * dimension(y) != dimension(z) && throw(DimensionError(x*y,z))
        uI = unit(x)*unit(y)
        uF = promote_unit(uI, unit(z))
        c = ($_x)(ustrip(x), ustrip(y), ustrip(uconvert(uI, z)))
        uconvert(uF, Quantity(c, uI))
    end
end

sqrt(x::AbstractQuantity) = Quantity(sqrt(x.val), sqrt(unit(x)))
cbrt(x::AbstractQuantity) = Quantity(cbrt(x.val), cbrt(unit(x)))

for _y in (:sin, :cos, :tan, :asin, :acos, :atan, :sinh, :cosh, :tanh, :asinh, :acosh, :atanh,
           :sinpi, :cospi, :sinc, :cosc, :cis)
    @eval ($_y)(x::DimensionlessQuantity) = ($_y)(uconvert(NoUnits, x))
end

atan(y::AbstractQuantity{T1,D,U1}, x::AbstractQuantity{T2,D,U2}) where {T1,T2,D,U1,U2} =
    atan(promote(y,x)...)
atan(y::AbstractQuantity{T,D,U}, x::AbstractQuantity{T,D,U}) where {T,D,U} = atan(y.val,x.val)
atan(y::AbstractQuantity, x::AbstractQuantity) = throw(DimensionError(x,y))

for (f, F) in [(:min, :<), (:max, :>)]
    @eval @generated function ($f)(x::AbstractQuantity, y::AbstractQuantity)    #TODO
        xdim = x.parameters[2]
        ydim = y.parameters[2]
        if xdim != ydim
            return :(throw(DimensionError(x,y)))
        end

        isa(x.parameters[3](), FixedUnits) &&
            isa(y.parameters[3](), FixedUnits) &&
            x.parameters[3] !== y.parameters[3] &&
            error("automatic conversion prohibited.")

        xunits = x.parameters[3].parameters[1]
        yunits = y.parameters[3].parameters[1]

        factx = mapreduce((x,y)->broadcast(*,x,y), xunits) do x
            vcat(basefactor(x)...)
        end
        facty = mapreduce((x,y)->broadcast(*,x,y), yunits) do x
            vcat(basefactor(x)...)
        end

        tensx = mapreduce(tensfactor, +, xunits)
        tensy = mapreduce(tensfactor, +, yunits)

        convx = *(factx..., (10.0)^tensx)
        convy = *(facty..., (10.0)^tensy)

        :($($F)(x.val*$convx, y.val*$convy) ? x : y)
    end
end

abs(x::AbstractQuantity) = Quantity(abs(x.val), unit(x))
abs2(x::AbstractQuantity) = Quantity(abs2(x.val), unit(x)*unit(x))
angle(x::AbstractQuantity{<:Complex}) = angle(x.val)

copysign(x::AbstractQuantity, y::Number) = Quantity(copysign(x.val,y/unit(y)), unit(x))
flipsign(x::AbstractQuantity, y::Number) = Quantity(flipsign(x.val,y/unit(y)), unit(x))

for (i,j) in zip((:<, :isless), (:_lt, :_isless))
    @eval ($i)(x::AbstractQuantity, y::AbstractQuantity) = ($j)(x,y)
    @eval ($i)(x::AbstractQuantity, y::Number) = ($i)(promote(x,y)...)
    @eval ($i)(x::Number, y::AbstractQuantity) = ($i)(promote(x,y)...)

    # promotion might not yield Quantity types
    @eval @inline ($j)(x::AbstractQuantity{T1}, y::AbstractQuantity{T2}) where {T1,T2} = ($i)(promote(x,y)...)
    # If it does yield Quantity types, we'll get back here,
    # since at least the numeric part can be promoted.
    @eval @inline ($j)(x::AbstractQuantity{T,D,U}, y::AbstractQuantity{T,D,U}) where {T,D,U} = ($i)(x.val,y.val)
    @eval @inline ($j)(x::AbstractQuantity{T,D,U1}, y::AbstractQuantity{T,D,U2}) where {T,D,U1,U2} = ($i)(promote(x,y)...)
    @eval @inline ($j)(x::AbstractQuantity{T,D1,U1}, y::AbstractQuantity{T,D2,U2}) where {T,D1,D2,U1,U2} = throw(DimensionError(x,y))
end

Base.rtoldefault(::Type{<:AbstractQuantity{T,D,U}}) where {T,D,U} = Base.rtoldefault(T)
isapprox(x::AbstractQuantity{T,D,U}, y::AbstractQuantity{T,D,U}; atol=zero(Quantity{real(T),D,U}), kwargs...) where {T,D,U} =
    isapprox(x.val, y.val; atol=uconvert(unit(y), atol).val, kwargs...)
function isapprox(x::AbstractQuantity, y::AbstractQuantity; kwargs...)
    dimension(x) != dimension(y) && return false
    return isapprox(promote(x,y)...; kwargs...)
end
isapprox(x::AbstractQuantity, y::Number; kwargs...) = isapprox(promote(x,y)...; kwargs...)
isapprox(x::Number, y::AbstractQuantity; kwargs...) = isapprox(y, x; kwargs...)

function isapprox(x::AbstractArray{<:AbstractQuantity{T1,D,U1}},
        y::AbstractArray{<:AbstractQuantity{T2,D,U2}}; rtol::Real=Base.rtoldefault(T1,T2,0),
        atol=zero(Quantity{T1,D,U1}), norm::Function=norm) where {T1,D,U1,T2,U2}

    d = norm(x - y)
    if isfinite(d)
        return d <= atol + rtol*max(norm(x), norm(y))
    else
        # Fall back to a component-wise approximate comparison
        return all(ab -> isapprox(ab[1], ab[2]; rtol=rtol, atol=atol), zip(x, y))
    end
end
isapprox(x::AbstractArray{S}, y::AbstractArray{T};
    kwargs...) where {S <: AbstractQuantity,T <: AbstractQuantity} = false
function isapprox(x::AbstractArray{S}, y::AbstractArray{N};
    kwargs...) where {S <: AbstractQuantity,N <: Number}
    if dimension(N) == dimension(S)
        isapprox(map(x->uconvert(NoUnits,x),x),y; kwargs...)
    else
        false
    end
end
isapprox(y::AbstractArray{N}, x::AbstractArray{S};
    kwargs...) where {S <: AbstractQuantity,N <: Number} = isapprox(x,y; kwargs...)

for cmp in [:(==), :isequal]
    @eval $cmp(x::AbstractQuantity{S,D,U}, y::AbstractQuantity{T,D,U}) where {S,T,D,U} = $cmp(x.val, y.val)
    @eval function $cmp(x::AbstractQuantity, y::AbstractQuantity)
        dimension(x) != dimension(y) && return false
        $cmp(promote(x,y)...)
    end

    @eval function $cmp(x::AbstractQuantity, y::Number)
        $cmp(promote(x,y)...)
    end
    @eval $cmp(x::Number, y::AbstractQuantity) = $cmp(y,x)
end
<=(x::AbstractQuantity, y::AbstractQuantity) = <(x,y) || x==y

_dimerr(f) = error("$f can only be well-defined for dimensionless ",
        "numbers. For dimensionful numbers, different input units yield physically ",
        "different results.")
isinteger(x::AbstractQuantity) = _dimerr(isinteger)
isinteger(x::DimensionlessQuantity) = isinteger(uconvert(NoUnits, x))

_rounderr() = error("specify the type of the quantity to convert to ",
    "when rounding quantities. Example: round(typeof(1u\"m\"), 137u\"cm\").")

# convenience methods
round(u::Units, q::AbstractQuantity, r::RoundingMode=RoundNearest; kwargs...) =
    Quantity(round(ustrip(u, q), r; kwargs...), u)
round(::Type{T}, u::Units, q::AbstractQuantity, r::RoundingMode=RoundNearest;
        kwargs...) where {T<:Number} =
    round(Quantity{T, dimension(u), typeof(u)}, q, r; kwargs...)

# workhorse methods
round(x::AbstractQuantity, r::RoundingMode=RoundNearest; kwargs...) =
    _rounderr()
round(x::DimensionlessQuantity; kwargs...) = round(uconvert(NoUnits, x); kwargs...)
round(x::DimensionlessQuantity, r::RoundingMode; kwargs...) =
    round(uconvert(NoUnits, x), r; kwargs...)
round(::Type{T}, x::AbstractQuantity, r=RoundingMode=RoundNearest;
    kwargs...) where {T<:Number} = _dimerr(:round)
round(::Type{T}, x::DimensionlessQuantity, r::RoundingMode=RoundNearest;
    kwargs...) where {T<:Number} = round(T, uconvert(NoUnits, x), r; kwargs...)
function round(::Type{T}, x::AbstractQuantity;
        kwargs...) where {S, T <: Quantity{S}}
    u = unit(T)
    unitless = ustrip(u, x)
    return Quantity{S, dimension(T), typeof(u)}(round(unitless; kwargs...))
end
function round(::Type{T}, x::AbstractQuantity, r::RoundingMode;
        kwargs...) where {S, T <: Quantity{S}}
    u = unit(T)
    unitless = ustrip(u, x)
    return Quantity{S, dimension(T), typeof(u)}(round(unitless, r; kwargs...))
end

# that should actually be fixed in Base ↓
for (f,r) = ((:trunc, :RoundToZero), (:floor, :RoundDown), (:ceil, :RoundUp))
    @eval $f(x::AbstractQuantity; kwargs...) = round(x, $r; kwargs...)
    @eval $f(::Type{T}, x::AbstractQuantity; kwargs...) where {T<:Number} =
        round(T, x, $r; kwargs...)
    @eval $f(u::Units, x::AbstractQuantity; kwargs...) = round(u, x, $r; kwargs...)
end

zero(x::AbstractQuantity) = Quantity(zero(x.val), unit(x))
zero(x::AffineQuantity) = Quantity(zero(x.val), absoluteunit(x))
zero(x::Type{<:AbstractQuantity{T}}) where {T} = throw(ArgumentError("zero($x) not defined."))
zero(x::Type{<:AbstractQuantity{T,D}}) where {T,D} = zero(T) * upreferred(D)
zero(x::Type{<:AbstractQuantity{T,D,U}}) where {T,D,U<:ScalarUnits} = zero(T)*U()
zero(x::Type{<:AbstractQuantity{T,D,U}}) where {T,D,U<:AffineUnits} = zero(T)*absoluteunit(U())

one(x::AbstractQuantity) = one(x.val)
one(x::AffineQuantity) =
    throw(AffineError("no multiplicative identity for affine quantity $x."))
oneunit(x::AffineQuantity) = Quantity(one(x.val), absoluteunit(x))
oneunit(x::Type{<:AbstractQuantity{T,D,U}}) where {T,D,U<:AffineUnits} = Quantity(one(T), absoluteunit(U()))
get_T(::Type{<:AbstractQuantity{T}}) where T = T
get_T(::Type{<:AbstractQuantity{T,D}}) where {T,D} = T
get_T(::Type{<:AbstractQuantity{T,D,U}}) where {T,D,U} = T
one(x::Type{<:AbstractQuantity}) = one(get_T(x))
one(x::Type{<:AffineQuantity}) =
    throw(AffineError("no multiplicative identity for affine quantity type $x."))

isreal(x::AbstractQuantity) = isreal(x.val)
isfinite(x::AbstractQuantity) = isfinite(x.val)
isinf(x::AbstractQuantity) = isinf(x.val)
isnan(x::AbstractQuantity) = isnan(x.val)
@static if VERSION ≥ v"1.7.0-DEV.119"
    isunordered(x::AbstractQuantity) = isunordered(x.val)
end

eps(x::T) where {T<:AbstractQuantity} = T(eps(x.val))
eps(x::Type{T}) where {T<:AbstractQuantity} = eps(Unitful.numtype(T))

unsigned(x::AbstractQuantity) = Quantity(unsigned(x.val), unit(x))

for f in (:exp, :exp10, :exp2, :expm1, :log, :log10, :log1p, :log2)
    @eval ($f)(x::DimensionlessQuantity) = ($f)(uconvert(NoUnits, x))
end

real(x::AbstractQuantity) = Quantity(real(x.val), unit(x))
imag(x::AbstractQuantity) = Quantity(imag(x.val), unit(x))
conj(x::AbstractQuantity) = Quantity(conj(x.val), unit(x))

@inline norm(x::AbstractQuantity, p::Real=2) =
    p == 0 ? (x==zero(x) ? typeof(abs(x))(0) : typeof(abs(x))(1)) : abs(x)

"""
    sign(x::AbstractQuantity)
Returns the sign of `x`.
"""
sign(x::AbstractQuantity) = sign(x.val)

"""
    signbit(x::AbstractQuantity)
Returns the sign bit of the underlying numeric value of `x`.
"""
signbit(x::AbstractQuantity) = signbit(x.val)

prevfloat(x::AbstractQuantity{T}) where {T <: AbstractFloat} = Quantity(prevfloat(x.val), unit(x))
nextfloat(x::AbstractQuantity{T}) where {T <: AbstractFloat} = Quantity(nextfloat(x.val), unit(x))

function frexp(x::AbstractQuantity{T}) where {T <: AbstractFloat}
    a,b = frexp(x.val)
    a*unit(x), b
end

for f in (:float, :BigFloat, :Float64, :Float32, :Float16)
    @eval begin
    """
        $($f)(x::AbstractQuantity)
    Convert the numeric backing type of `x` to a floating-point representation.
    Returns a `Quantity` with the same units.
    """
    (Base.$f)(x::AbstractQuantity) = Quantity($f(x.val), unit(x))
    end
end
   
"""
    Integer(x::AbstractQuantity)
Convert the numeric backing type of `x` to an integer representation.
Returns a `Quantity` with the same units.
"""
Integer(x::AbstractQuantity) = Quantity(Integer(x.val), unit(x))

"""
    Rational(x::AbstractQuantity)
Convert the numeric backing type of `x` to a rational number representation.
Returns a `Quantity` with the same units.
"""
Rational(x::AbstractQuantity) = Quantity(Rational(x.val), unit(x))

typemin(::Type{<:AbstractQuantity{T,D,U}}) where {T,D,U} = typemin(T)*U()
typemin(x::AbstractQuantity{T}) where {T} = typemin(T)*unit(x)

typemax(::Type{<:AbstractQuantity{T,D,U}}) where {T,D,U} = typemax(T)*U()
typemax(x::AbstractQuantity{T}) where {T} = typemax(T)*unit(x)

Base.literal_pow(::typeof(^), x::AbstractQuantity, ::Val{v}) where {v} =
    Quantity(Base.literal_pow(^, x.val, Val(v)),
             Base.literal_pow(^, unit(x), Val(v)))

# All of these are needed for ambiguity resolution
^(x::AbstractQuantity, y::Integer) = Quantity((x.val)^y, unit(x)^y)
^(x::AbstractQuantity, y::Rational) = Quantity((x.val)^y, unit(x)^y)
^(x::AbstractQuantity, y::Real) = Quantity((x.val)^y, unit(x)^y)

Base.rand(r::Random.AbstractRNG, ::Random.SamplerType{<:AbstractQuantity{T,D,U}}) where {T,D,U} =
    rand(r, T) * U()
Base.ones(Q::Type{<:AbstractQuantity}, dims::NTuple{N,Integer}) where {N} =
    fill!(Array{Q,N}(undef, map(Base.to_dim, dims)), oneunit(Q))
Base.ones(Q::Type{<:AbstractQuantity}, dims::Tuple{}) = fill!(Array{Q}(undef), oneunit(Q))
Base.ones(a::AbstractArray, Q::Type{<:AbstractQuantity}) = fill!(similar(a,Q), oneunit(Q))
