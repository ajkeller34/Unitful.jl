using Unitful
using Test, LinearAlgebra, Random, ConstructionBase
import Unitful: DimensionError, AffineError
import Unitful: LogScaled, LogInfo, Level, Gain, MixedUnits, Decibel
import Unitful: FreeUnits, ContextUnits, FixedUnits, AffineUnits, AffineQuantity

import Unitful:
    nm, μm, mm, cm, m, km, inch, ft, mi,
    ac,
    mg, g, kg,
    Ra, °F, °C, K,
    rad, °,
    ms, s, minute, hr, d, yr, Hz,
    J, A, N, mol, V,
    mW, W,
    dB, dB_rp, dB_p, dBm, dBV, dBSPL, Decibel,
    Np, Np_rp, Np_p, Neper

import Unitful: 𝐋, 𝐓, 𝐍, 𝚯

import Unitful:
    Length, Area, Volume,
    Luminosity,
    Time, Frequency,
    Mass,
    Current,
    Temperature, AbsoluteScaleTemperature, RelativeScaleTemperature,
    Action,
    Power,
    MassFlow,
    MolarFlow,
    VolumeFlow

import Unitful: LengthUnits, AreaUnits, MassUnits, TemperatureUnits

using Dates:
    Dates,
    Nanosecond, Microsecond, Millisecond, Second, Minute, Hour, Day, Week,
    Month, Year,
    CompoundPeriod

const colon = Base.:(:)

@testset "Construction" begin
    @test isa(NoUnits, FreeUnits)
    @test typeof(𝐋) === Unitful.Dimensions{(Unitful.Dimension{:Length}(1),)}
    @test 𝐋*𝐋 === 𝐋^2
    @test typeof(1.0m) === Unitful.Quantity{Float64, 𝐋,
        Unitful.FreeUnits{(Unitful.Unit{:Meter, 𝐋}(0,1),), 𝐋, nothing}}
    @test typeof(1m^2) === Unitful.Quantity{Int, 𝐋^2,
            Unitful.FreeUnits{(Unitful.Unit{:Meter, 𝐋}(0,2),), 𝐋^2, nothing}}
    @test typeof(1ac) === Unitful.Quantity{Int, 𝐋^2,
            Unitful.FreeUnits{(Unitful.Unit{:Acre, 𝐋^2}(0,1),), 𝐋^2, nothing}}
    @test typeof(ContextUnits(m,μm)) ===
        ContextUnits{(Unitful.Unit{:Meter, 𝐋}(0,1),), 𝐋, typeof(μm), nothing}
    @test typeof(1.0*ContextUnits(m,μm)) === Unitful.Quantity{Float64, 𝐋,
        ContextUnits{(Unitful.Unit{:Meter, 𝐋}(0,1),), 𝐋, typeof(μm), nothing}}
    @test typeof(1.0*FixedUnits(m)) === Unitful.Quantity{Float64, 𝐋,
        FixedUnits{(Unitful.Unit{:Meter, 𝐋}(0,1),), 𝐋, nothing}}
    @test 3mm != 3*(m*m)                        # mm not interpreted as m*m
    @test (3+4im)*V === V*(3+4im) === (3V+4V*im)  # Complex quantity construction
    @test !isreal(Base.complex(3.0/m, 4.0/m))
    @test !isreal(Base.complex((3.0+4.0im)/m))
    @test Base.reim(Base.complex((3.0+4.0im)/m)) == (3.0/m, 4.0/m)
    @test Base.complex(1m, 1.5m) == Base.complex(1.0m, 1.5m)
    @test Base.widen(Base.complex(Float32(3.0)/m)) == Base.complex(Float64(3.0)/m)
    @test Base.complex(1.0/m) == (1.0/m + (0.0/m)*im)
    @test Base.complex(1.0/m + (1.5/m)*im) == (1.0/m + (1.5/m)*im)
    @test 3*NoUnits === 3
    @test 3*(FreeUnits(m)/FreeUnits(m)) === 3
    @test 3*(ContextUnits(m)/ContextUnits(m)) === 3
    @test 3*(FixedUnits(m)/FixedUnits(m)) === 3
    @test ContextUnits(mm) === ContextUnits(mm,m)
    @test Quantity(3, NoUnits) === 3
    @test FreeUnits(ContextUnits(mm,m)) === FreeUnits(mm)
    @test FreeUnits(FixedUnits(mm)) === FreeUnits(mm)
    @test isa(FreeUnits(m), FreeUnits)
    @test isa(ContextUnits(m), ContextUnits)
    @test isa(ContextUnits(m,mm), ContextUnits)
    @test isa(FixedUnits(m), FixedUnits)
    @test ContextUnits(m, FixedUnits(mm)) === ContextUnits(m, mm)
    @test ContextUnits(m, ContextUnits(mm, mm)) === ContextUnits(m, mm)
    @test_throws DimensionError ContextUnits(m,kg)
    @test ConstructionBase.constructorof(typeof(1.0m))(2) === 2m
end

@testset "LinearAlgebra functions" begin
    CNT = Ref(0)
    Unitful._linearalgebra_count() = (CNT[] += 1; nothing)
    @testset "> Matrix multiplication: *" begin
        M = rand(3,3) .* u"m"
        M_ = view(M,:,1:3)
        v = rand(3) .* u"V"
        v_ = view(v, 1:3)

        CNT[] = 0

        @test unit(first(M * M)) == u"m*m"
        @test M * M == M_ * M == M * M_ == M_ * M_

        @test unit(first(M * v)) == u"m*V"
        @test M * v == M_ * v == M * v_ == M_ * v_

        VERSION >= v"1.3" && @test CNT[] == 10

        @test unit(first(v' * M)) == u"m*V"
        @test v' * M == v_' * M == v_' * M == v_' * M_

        VERSION >= v"1.3" && @test CNT[] == 15

        @test unit(v' * v) == u"V*V"
        @test v' * v == v_' * v == v_' * v == v_' * v_

        VERSION >= v"1.3" && @test CNT[] == 20

        # Mixed with & without units
        N = rand(3,3)
        w = rand(3)

        CNT[] = 0
        
        @test unit(first(M * N)) == u"m"
        @test unit(first(N * M)) == u"m"

        @test unit(first(M * w)) == u"m"
        @test unit(first(N * v)) == u"V"

        @show CNT[] # not specialised yet

    end
    @testset "> Matrix multiplication: mul!" begin
        A = rand(3,3) .* u"m"
        B = rand(3,3) .* u"m"
        C = fill(zero(eltype(A*B)), 3, 3)
        CNT[] = 0

        mul!(C, A, B)
        if VERSION >= v"1.3" # the 5-arm mul! exists
            mul!(C, A, B, true, true)
            mul!(C, A, B, 3, 7) # not specialised yet

            @show CNT[]
        end
    end
end

@testset "Types" begin
    @test Base.complex(Quantity{Float64,NoDims,NoUnits}) ==
        Quantity{Complex{Float64},NoDims,NoUnits}
end

@testset "Conversion" begin
    @testset "> Unitless ↔ unitful conversion" begin
        @test_throws DimensionError convert(typeof(3m), 1)
        @test_throws DimensionError convert(Quantity{Float64, typeof(𝐋)}, 1)
        @test_throws DimensionError convert(Float64, 3m)
        @test @inferred(3m/unit(3m)) === 3
        @test @inferred(3.0g/unit(3.0g)) === 3.0
        # 1-arg ustrip
        @test @inferred(ustrip(3*FreeUnits(m))) === 3
        @test @inferred(ustrip(3*ContextUnits(m,μm))) === 3
        @test @inferred(ustrip(3*FixedUnits(m))) === 3
        @test @inferred(ustrip(3)) === 3
        @test @inferred(ustrip(3.0m)) === 3.0
        # ustrip with type and unit arguments
        @test @inferred(ustrip(m, 3.0m)) === 3.0
        @test @inferred(ustrip(m, 2mm)) === 1//500
        @test @inferred(ustrip(mm, 3.0m)) === 3000.0
        @test @inferred(ustrip(NoUnits, 3.0m/1.0m)) === 3.0
        @test @inferred(ustrip(NoUnits, 3.0m/1.0cm)) === 300.0
        @test @inferred(ustrip(cm, missing)) === missing
        @test @inferred(ustrip(NoUnits, missing)) === missing
        @test_throws DimensionError ustrip(NoUnits, 3.0m/1.0s)
        @test @inferred(ustrip(Float64, m, 2mm)) === 0.002
        @test @inferred(ustrip(Int, mm, 2.0m)) === 2000
        @test @inferred(ustrip(Float32, NoUnits, 5.0u"m"/2.0u"m")) === Float32(2.5)
        @test @inferred(ustrip(Int, NoUnits, 3.0u"m"/1.0u"cm")) === 300
        # convert
        @test convert(typeof(1mm/m), 3) == 3000mm/m
        @test convert(typeof(1mm/m), 3*NoUnits) == 3000mm/m
        @test convert(typeof(1*ContextUnits(mm/m, NoUnits)), 3) == 3000mm/m
        @test convert(typeof(1*FixedUnits(mm/m)), 3) == 3000*FixedUnits(mm/m)
        @test convert(Int, 1*FreeUnits(m/mm)) === 1000
        @test convert(Int, 1*FixedUnits(m/mm)) === 1000
        @test convert(Int, 1*ContextUnits(m/mm, NoUnits)) === 1000

        # w/ units distinct from w/o units
        @test 1m != 1
        @test 1 != 1m
        @test (3V+4V*im) != (3+4im)

        # Issue 26
        @unit altL "altL" altLiter 1000*cm^3 true
        @test convert(Float64, 1altL/cm^3) === 1000.0
        # Issue 327
        @test uconvert(u"√cm", 1u"√m") == 10u"√cm"
    end
    @testset "> Unitful ↔ unitful conversion" begin
        @testset ">> Numeric conversion" begin
            @test @inferred(float(3m)) === 3.0m
            @test @inferred(Float32(3m)) === 3.0f0m
            @test @inferred(Integer(3.0A)) === 3A
            @test Rational(3.0m) === (Int64(3)//1)*m
            @test typeof(convert(typeof(0.0°), 90°)) == typeof(0.0°)
            @test (3.0+4.0im)*V == (3+4im)*V
            @test im*V == Complex(0,1)*V
        end
        @testset ">> Intra-unit conversion" begin
            # an essentially no-op uconvert should not disturb numeric type
            @test @inferred(uconvert(g,1g)) === 1g
            @test @inferred(uconvert(m,0x01*m)) === 0x01*m
            @test @inferred(convert(Quantity{Float64, 𝐋}, 1m)) === 1.0m
            @test 1kg === 1kg
            @test typeof(1m)(1m) === 1m

            @test 1*FreeUnits(kg) == 1*ContextUnits(kg,g)
            @test 1*ContextUnits(kg,g) == 1*FreeUnits(kg)

            @test 1*FreeUnits(kg) == 1*FixedUnits(kg)
            @test 1*FixedUnits(kg) == 1*FreeUnits(kg)

            @test 1*ContextUnits(kg,g) == 1*FixedUnits(kg)
            @test 1*FixedUnits(kg) == 1*ContextUnits(kg,g)

            # No auto conversion when working with FixedUnits exclusively.
            @test_throws ErrorException 1*FixedUnits(kg) == 1000*FixedUnits(g)
        end
        @testset ">> Inter-unit conversion" begin
            @test 1g == 0.001kg                   # Issue 56
            @test 0.001kg == 1g                   # Issue 56
            @test 1kg == 1000g
            @test !(1kg === 1000g)
            @test 1inch == (254//100)*cm
            @test 1ft == 12inch
            @test 1/mi == 1//(5280ft)
            @test 1minute == 60s
            @test 1hr == 60minute
            @test 1d == 24hr
            @test 1yr == 365.25d
            @test 1J == 1kg*m^2/s^2
            @test typeof(1cm)(1m) === 100cm
            @test (3V+4V*im) != (3m+4m*im)
            @test_throws DimensionError uconvert(m, 1kg)
            @test_throws DimensionError uconvert(m, 1*ContextUnits(kg,g))
            @test_throws DimensionError uconvert(ContextUnits(m,mm), 1kg)
            @test_throws DimensionError uconvert(m, 1*FixedUnits(kg))
            @test uconvert(g, 1*FixedUnits(kg)) == 1000g         # manual conversion okay
            @test (1kg, 2g, 3mg, missing) .|> g === (1000g, 2g, (3//1000)g, missing)
            # Issue 79:
            @test isapprox(upreferred(Unitful.ɛ0), 8.85e-12u"F/m", atol=0.01e-12u"F/m")
            # Issue 261:
            @test 1u"rps" == 360°/s
            @test 1u"rps" == 2π/s
            @test 1u"rpm" == 360°/minute
            @test 1u"rpm" == 2π/minute
        end
    end
end

include("dates.jl")

@testset "Temperature and affine quantities" begin
    @testset "Affine transforms and quantities" begin
        @test 1°C isa RelativeScaleTemperature
        @test !isa(1°C, AbsoluteScaleTemperature)
        @test 1K isa AbsoluteScaleTemperature
        @test !isa(1K, RelativeScaleTemperature)

        @test_throws AffineError °C*°C
        @test_throws AffineError °C*K
        @test_throws AffineError (0°C)*(0°C)
        @test_throws AffineError (1°C)/(1°C)
        @test_throws AffineError °C^2
        let x = 2
            @test_throws AffineError °C^x
        end
        @test_throws AffineError inv(°C)
        @test_throws AffineError inv(0°C)
        @test_throws AffineError sqrt(°C)
        @test_throws AffineError sqrt(0°C)
        @test_throws AffineError cbrt(°C)
        @test_throws AffineError cbrt(0°C)
        @test_throws AffineError 32°F + 1°F
        @test_throws AffineError (32°F) * 2
        @test_throws AffineError 2 * (32°F)
        @test_throws AffineError (32°F) / 2
        @test_throws AffineError 2 / (32°F)

        for f = (:div, :rem, :divrem)
            @eval for r = (RoundNearest, RoundNearestTiesAway, RoundNearestTiesUp,
                           RoundToZero, RoundUp, RoundDown)
                @test_throws AffineError $f(32°F, 2°F, r)
                @test_throws AffineError $f(32°F, 2K, r)
                @test_throws AffineError $f(32K, 2°F, r)
            end
        end
        for f = (:div, :cld, :fld, :rem, :mod, :divrem, :fldmod)
            @eval begin
                @test_throws AffineError $f(32°F, 2°F)
                @test_throws AffineError $f(32°F, 2K)
                @test_throws AffineError $f(32K, 2°F)
            end
        end

        @test zero(100°C) === 0K
        @test zero(typeof(100°C)) === 0K
        @test oneunit(100°C) === 1K
        @test oneunit(typeof(100°C)) === 1K
        @test_throws AffineError one(100°C)
        @test_throws AffineError one(typeof(100°C))

        @test 0°C isa AffineQuantity{T, 𝚯} where T    # is "relative temperature"
        @test 0°C isa Temperature                             # dimensional correctness
        @test °C isa AffineUnits{N, 𝚯} where N
        @test °C isa TemperatureUnits

        @test @inferred(uconvert(°F, 0°C))  === (32//1)°F   # Some known conversions...
        @test @inferred(uconvert(°C, 32°F)) === (0//1)°C    #  ⋮
        @test @inferred(uconvert(°C, 212°F)) === (100//1)°C #  ⋮
        @test @inferred(uconvert(°C, 0x01*°C)) === 0x01*°C  # Preserve numeric type

        # The next test is a little funky but checks the `affineunit` functionality
        @test @inferred(uconvert(°F,
            0*Unitful.affineunit(27315K//100 + 5K//9))) === (33//1)°F
    end
    @testset "Temperature differences" begin
        @test @inferred(uconvert(Ra, 0K)) === 0Ra//1
        @test @inferred(uconvert(K, 1Ra)) === 5K//9
        @test @inferred(uconvert(μm/(m*Ra), 9μm/(m*K))) === 5μm/(m*Ra)//1

        @test @inferred(uconvert(FreeUnits(Ra), 4.2K)) ≈ 7.56Ra
        @test @inferred(unit(uconvert(FreeUnits(Ra), 4.2K))) === FreeUnits(Ra)
        @test @inferred(uconvert(FreeUnits(Ra), 4.2*ContextUnits(K))) ≈ 7.56Ra
        @test @inferred(unit(uconvert(FreeUnits(Ra), 4.2*ContextUnits(K)))) === FreeUnits(Ra)
        @test @inferred(unit(uconvert(ContextUnits(Ra), 4.2K))) === ContextUnits(Ra)

        let cc = ContextUnits(°C, °C), kc = ContextUnits(K, °C), rac = ContextUnits(Ra, °C)
            @test 100°C + 1K === (7483//20)K
            @test 100cc + 1K === (101//1)cc
            @test 100cc + 1K == (101//1)°C
            @test 1K + 100cc === (101//1)cc
            @test 1K + 100cc == (101//1)°C
            @test 100°C + 1Ra === (67267//180)K
            @test 100°C - 212°F === (0//1)K
            @test 100°C - 211°F === (5//9)K
            @test 100°C - 1°C === 99K
            @test 100°C - 32°F === (100//1)K
            @test 10cc + 2.0K/hr * 60minute + 3.0K/hr * 60minute === 15.0cc
            @test 10cc + 5kc === (15//1)cc
            @test 10°C + 5kc === (15//1)cc
            @test 10°C + (9//5)rac === (11//1)cc
        end
    end
    @testset "Promotion" begin
        @test_throws ErrorException Unitful.preferunits(°C)
        @test @inferred(eltype([1°C, 1K])) <: Quantity{Rational{Int}, 𝚯, typeof(K)}
        @test @inferred(eltype([1.0°C, 1K])) <: Quantity{Float64, 𝚯, typeof(K)}
        @test @inferred(eltype([1°C, 1°F])) <: Quantity{Rational{Int}, 𝚯, typeof(K)}
        @test @inferred(eltype([1.0°C, 1°F])) <: Quantity{Float64, 𝚯, typeof(K)}

        # context units should be identifiable as affine
        @test ContextUnits(°C, °F) isa AffineUnits

        let fc = ContextUnits(°F, °C), cc = ContextUnits(°C, °C)
            @test @inferred(promote(1fc, 1cc)) === ((-155//9)cc, (1//1)cc)
            @test @inferred(eltype([1cc, 1°C])) <: Quantity{Rational{Int}, 𝚯, typeof(cc)}
        end
    end
end

# preferred units work on AbstractQuantity
struct QQQ <: Unitful.AbstractQuantity{Float64,𝐋,typeof(cm)}
    val::Float64
end
Unitful.uconvert(U::Unitful.Units, q::QQQ) = uconvert(U, Quantity(q.val, cm))

@testset "Promotion" begin
    @testset "> Unit preferences" begin
        # Only because we favor SI, we have the following:
        @test @inferred(upreferred(N)) === kg*m/s^2
        @test @inferred(upreferred(dimension(N))) === kg*m/s^2
        @test @inferred(upreferred(g)) === kg
        @test @inferred(upreferred(FreeUnits(g))) === FreeUnits(kg)

        # Test special units behaviors
        @test @inferred(upreferred(ContextUnits(g,mg))) === ContextUnits(mg,mg)
        @test @inferred(upreferred(FixedUnits(kg))) === FixedUnits(kg)
        @test @inferred(upreferred(upreferred(1.0ContextUnits(kg,g)))) ===
            1000.0ContextUnits(g,g)
        @test @inferred(upreferred(unit(1g |> ContextUnits(g,mg)))) === ContextUnits(mg,mg)
        @test @inferred(upreferred(1g |> ContextUnits(g,mg))) == 1000mg

        @test @inferred(upreferred(1N)) === 1*kg*m/s^2
        @test ismissing(upreferred(missing))

        # preferred units work on AbstractQuantity
        @test @inferred(upreferred(QQQ(10))) == 0.1m
    end
    @testset "> promote_unit" begin
        @test Unitful.promote_unit(FreeUnits(m)) === FreeUnits(m)
        @test Unitful.promote_unit(ContextUnits(m,mm)) === ContextUnits(m,mm)
        @test Unitful.promote_unit(FixedUnits(kg)) === FixedUnits(kg)
        @test Unitful.promote_unit(ContextUnits(m,mm), ContextUnits(km,mm)) ===
            ContextUnits(mm,mm)
        @test Unitful.promote_unit(FreeUnits(m), ContextUnits(mm,km)) ===
            ContextUnits(km,km)
        @test Unitful.promote_unit(FixedUnits(kg), ContextUnits(g,g)) === FixedUnits(kg)
        @test Unitful.promote_unit(ContextUnits(g,g), FixedUnits(kg)) === FixedUnits(kg)
        @test Unitful.promote_unit(FixedUnits(kg), FreeUnits(g)) === FixedUnits(kg)
        @test Unitful.promote_unit(FreeUnits(g), FixedUnits(kg)) === FixedUnits(kg)
        @test_throws DimensionError Unitful.promote_unit(m,kg)

        # FixedUnits throw a promotion error
        @test_throws ErrorException Unitful.promote_unit(FixedUnits(m), FixedUnits(mm))

        # Only because we favor SI, we have the following:
        @test Unitful.promote_unit(m,km) === m
        @test Unitful.promote_unit(m,km,cm) === m
        @test Unitful.promote_unit(ContextUnits(m,mm), ContextUnits(km,cm)) ===
            FreeUnits(m)
    end
    @testset "> Simple promotion" begin
        # promotion should do nothing to units alone
        # promote throws an error if no types are be changed
        @test_throws ErrorException promote(m, km)
        @test_throws ErrorException promote(ContextUnits(m, km), ContextUnits(mm, km))
        @test_throws ErrorException promote(FixedUnits(m), FixedUnits(km))

        # promote the numeric type
        @test @inferred(promote(1.0m, 1m)) === (1.0m, 1.0m)
        @test @inferred(promote(1m, 1.0m)) === (1.0m, 1.0m)
        @test @inferred(promote(1.0g, 1kg)) === (0.001kg, 1.0kg)
        @test @inferred(promote(1g, 1.0kg)) === (0.001kg, 1.0kg)
        @test @inferred(promote(1.0m, 1kg)) === (1.0m, 1.0kg)
        @test @inferred(promote(1kg, 1.0m)) === (1.0kg, 1.0m)
        @test_broken @inferred(promote(1.0m, 1)) === (1.0m, 1.0)         # issue 52
        @test @inferred(promote(π, 180°)) === (float(π), float(π))       # issue 168
        @test @inferred(promote(180°, π)) === (float(π), float(π))       # issue 168

        # prefer no units for dimensionless numbers
        @test @inferred(promote(1.0mm/m, 1.0km/m)) === (0.001,1000.0)
        @test @inferred(promote(1.0cm/m, 1.0mm/m, 1.0km/m)) === (0.01,0.001,1000.0)
        @test @inferred(promote(1.0rad,1.0°)) === (1.0,π/180.0)

        # Quantities with promotion context
        # Context overrides free units
        nm2μm = ContextUnits(nm,μm)
        μm2μm = ContextUnits(μm,μm)
        μm2mm = ContextUnits(μm,mm)
        @test @inferred(promote(1.0nm2μm, 1.0m)) === (0.001μm2μm, 1e6μm2μm)
        @test @inferred(promote(1.0m, 1.0μm2μm)) === (1e6μm2μm, 1.0μm2μm)
        @test ===(upreferred.(unit.(promote(1.0nm2μm, 2nm2μm)))[1], ContextUnits(μm,μm))
        @test ===(upreferred.(unit.(promote(1.0nm2μm, 2nm2μm)))[2], ContextUnits(μm,μm))

        # Context agreement
        @test @inferred(promote(1.0nm2μm, 1.0μm2μm)) ===
            (0.001μm2μm, 1.0μm2μm)
        @test @inferred(promote(1μm2μm, 1.0nm2μm, 1.0m)) ===
            (1.0μm2μm, 0.001μm2μm, 1e6μm2μm)
        @test @inferred(promote(1μm2μm, 1.0nm2μm, 1.0s)) ===
            (1.0μm2μm, 1.0nm2μm, 1.0s)
        # Context disagreement: fall back to free units
        @test @inferred(promote(1.0nm2μm, 1.0μm2mm)) === (1e-9m, 1e-6m)
    end
    @testset "> Promotion during array creation" begin
        @test typeof([1.0m,1.0m]) == Array{typeof(1.0m),1}
        @test typeof([1.0m,1m]) == Array{typeof(1.0m),1}
        @test typeof([1.0m,1cm]) == Array{typeof(1.0m),1}
        @test typeof([1kg,1g]) == Array{typeof(1kg//1),1}
        @test typeof([1.0m,1]) == Array{Quantity{Float64},1}
        @test typeof([1.0m,1kg]) == Array{Quantity{Float64},1}
        @test typeof([1.0m/s 1; 1 0]) == Array{Quantity{Float64},2}
    end
    @testset "> Issue 52" begin
        x,y = 10m, 1
        px,py = promote(x,y)

        # promoting the second time should not change the types
        @test_throws ErrorException promote(px, py)
    end
    @testset "> Some internal behaviors" begin
        # quantities
        @test Unitful.numtype(Quantity{Float64}) <: Float64
        @test Unitful.numtype(Quantity{Float64, 𝐋}) <: Float64
        @test Unitful.numtype(typeof(1.0kg)) <: Float64
        @test Unitful.numtype(1.0kg) <: Float64
    end
end

@testset "Unit string parsing" begin
    @test uparse("m") == m
    @test uparse("m,s") == (m,s)
    @test uparse("1.0") == 1.0
    @test uparse("m/s") == m/s
    @test uparse("N*m") == N*m
    @test uparse("1.0m/s") == 1.0m/s
    @test uparse("m^-1") == m^-1
    @test uparse("dB/Hz") == dB/Hz
    @test uparse("3.0dB/Hz") == 3.0dB/Hz

    # Invalid unit strings
    @test_throws Meta.ParseError uparse("N m")
    @test_throws ArgumentError uparse("abs(2)")
    @test_throws ArgumentError uparse("(1,2)")
    @test_throws ArgumentError uparse("begin end")

    # test ustrcheck_bool
    @test_throws ArgumentError uparse("basefactor") # non-Unit symbols
    # ustrcheck_bool(::Quantity)
    @test uparse("h") == Unitful.h
    @test uparse("π") == π              # issue 112
end

@testset "Unit and dimensional analysis" begin
    @test @inferred(unit(1m^2)) === m^2
    @test @inferred(unit(typeof(1m^2))) === m^2
    @test @inferred(unit(Float64)) === NoUnits
    @test @inferred(unit(Union{typeof(1m^2),Missing})) === m^2
    @test @inferred(unit(Union{Float64,Missing})) === NoUnits
    @test @inferred(unit(missing)) === missing
    @test @inferred(unit(Missing)) === missing
    @test @inferred(dimension(1m^2)) === 𝐋^2
    @test @inferred(dimension(1*ContextUnits(m,km)^2)) === 𝐋^2
    @test @inferred(dimension(typeof(1m^2))) === 𝐋^2
    @test @inferred(dimension(Float64)) === NoDims
    @test @inferred(dimension(m^2)) === 𝐋^2
    @test @inferred(dimension(1m/s)) === 𝐋/𝐓
    @test @inferred(dimension(m/s)) === 𝐋/𝐓
    @test @inferred(dimension(1u"mol")) === 𝐍
    @test @inferred(dimension(μm/m)) === NoDims
    @test @inferred(dimension(missing)) === missing
    @test @inferred(dimension(Missing)) === missing
    @test dimension.([1u"m", 1u"s"]) == [𝐋, 𝐓]
    @test dimension.([u"m", u"s"]) == [𝐋, 𝐓]
    @test (𝐋/𝐓)^2 === 𝐋^2 / 𝐓^2
    @test isa(m, LengthUnits)
    @test isa(ContextUnits(m,km), LengthUnits)
    @test isa(FixedUnits(m), LengthUnits)
    @test !isa(m, AreaUnits)
    @test !isa(ContextUnits(m,km), AreaUnits)
    @test !isa(FixedUnits(m), AreaUnits)
    @test !isa(m, MassUnits)
    @test isa(m^2, AreaUnits)
    @test !isa(m^2, LengthUnits)
    @test isa(1m, Length)
    @test isa(1*ContextUnits(m,km), Length)
    @test isa(1*FixedUnits(m), Length)
    @test !isa(1m, LengthUnits)
    @test !isa(1m, Area)
    @test !isa(1m, Luminosity)
    @test isa(1ft, Length)
    @test isa(1m^2, Area)
    @test !isa(1m^2, Length)
    @test isa(1inch^3, Volume)
    @test isa(1/s, Frequency)
    @test isa(1kg, Mass)
    @test isa(1s, Time)
    @test isa(1A, Current)
    @test isa(1K, Temperature)
    @test isa(1u"cd", Luminosity)
    @test isa(2π*rad*1.0m, Length)
    @test isa(u"h", Action)
    @test isa(3u"dBm", Power)
    @test isa(3u"dBm*Hz*s", Power)
    @test isa(1kg/s, MassFlow)
    @test isa(1mol/s, MolarFlow)
    @test isa(1m^3/s, VolumeFlow)
end

@testset "Mathematics" begin
    @testset "> Comparisons" begin
        # make sure we are just picking one of the arguments, without surprising conversions
        # happening to the units...
        @test min(1FreeUnits(hr), 1FreeUnits(s)) === 1FreeUnits(s)
        @test min(1FreeUnits(hr), 1ContextUnits(s,ms)) === 1ContextUnits(s,ms)
        @test min(1ContextUnits(s,ms), 1FreeUnits(hr)) === 1ContextUnits(s,ms)
        @test min(1ContextUnits(s,minute), 1ContextUnits(hr,s)) ===
            1ContextUnits(s,minute)
        @test max(1ContextUnits(ft, mi), 1ContextUnits(m,mm)) ===
            1ContextUnits(m,mm)
        @test min(1FreeUnits(hr), 1FixedUnits(s)) === 1FixedUnits(s)
        @test min(1FixedUnits(s), 1FreeUnits(hr)) === 1FixedUnits(s)
        @test min(1FixedUnits(s), 1ContextUnits(ms,s)) === 1ContextUnits(ms,s)
        @test min(1ContextUnits(ms,s), 1FixedUnits(s)) === 1ContextUnits(ms,s)

        # automatic conversion prohibited
        @test_throws ErrorException min(1FixedUnits(s), 1FixedUnits(hr))
        @test min(1FixedUnits(s), 1FixedUnits(s)) === 1FixedUnits(s)

        # now move on, presuming that's working well.
        @test max(1ft, 1m) == 1m
        @test max(10J, 1kg*m^2/s^2) === 10J
        @test max(1J//10, 1kg*m^2/s^2) === 1kg*m^2/s^2
        @test @inferred(0m < 1.0m)
        @test @inferred(2.0m < 3.0m)
        @test @inferred(2.0m .< 3.0m)
        @test !(@inferred 3.0m .< 3.0m)
        @test @inferred(2.0m <= 3.0m)
        @test @inferred(3.0m <= 3.0m)
        @test @inferred(2.0m .<= 3.0m)
        @test @inferred(3.0m .<= 3.0m)
        @test @inferred(1μm/m < 1)
        @test @inferred(1 > 1μm/m)
        @test @inferred(1μm/m < 1mm/m)
        @test @inferred(1mm/m > 1μm/m)
        @test_throws DimensionError 1m < 1kg
        @test_throws DimensionError 1m < 1
        @test_throws DimensionError 1 < 1m
        @test_throws DimensionError 1mm/m < 1m
        @test Base.rtoldefault(typeof(1.0u"m")) === Base.rtoldefault(typeof(1.0))
        @test Base.rtoldefault(typeof(1u"m")) === Base.rtoldefault(Int)
    end
    @testset "> Addition and subtraction" begin
        @test @inferred(+(1A)) == 1A                    # Unary addition
        @test @inferred(3m + 3m) == 6m                  # Binary addition
        @test @inferred(-(1kg)) == (-1)*kg              # Unary subtraction
        @test @inferred(3m - 2m) == 1m                  # Binary subtraction
        @test @inferred(zero(1m)) === 0m                # Additive identity
        @test @inferred(zero(typeof(1m))) === 0m
        @test @inferred(zero(typeof(1.0m))) === 0.0m
        @test_throws ArgumentError zero(Quantity{Int})
        @test zero(Quantity{Int, 𝐋}) == 0m
        @test zero(Quantity{Int, 𝐋}) isa Quantity{Int}
        @test @inferred(π/2*u"rad" + 90u"°") ≈ π        # Dimless quantities
        @test @inferred(π/2*u"rad" - 90u"°") ≈ 0        # Dimless quantities
        @test_throws DimensionError 1+1m                # Dim mismatched
        @test_throws DimensionError 1-1m
    end
    @testset "> Multiplication" begin
        @test @inferred(FreeUnits(g)*FreeUnits(kg)) === FreeUnits(g*kg)
        @test @inferred(ContextUnits(m,mm)*ContextUnits(kg,g)) ===
            ContextUnits(m*kg,mm*g)
        @test @inferred(ContextUnits(m,mm)*FreeUnits(g)) ===
            ContextUnits(m*g,mm*kg)
        @test @inferred(FreeUnits(g)*ContextUnits(m,mm)) ===
            ContextUnits(m*g,mm*kg)
        @test @inferred(FixedUnits(g)*FixedUnits(kg)) === FixedUnits(g*kg)
        @test @inferred(FixedUnits(g)*FreeUnits(kg)) === FixedUnits(g*kg)
        @test @inferred(FixedUnits(g)*ContextUnits(kg,kg)) === FixedUnits(g*kg)
        @test @inferred(*(1s)) === 1s                     # Unary multiplication
        @test @inferred(3m * 2cm) === 3cm * 2m            # Binary multiplication
        @test @inferred((3m)*m) === 3*(m*m)               # Associative multiplication
        @test @inferred(true*1kg) === 1kg                 # Boolean multiplication (T)
        @test @inferred(false*1kg) === 0kg                # Boolean multiplication (F)
        @test typeof(one(eltype([1.0s, 1kg]))) <: Float64 # issue 159, multiplicative identity
    end
    @testset "> Division" begin
        @test 360° / 2 === 180.0°            # Issue 110
        @test 360° // 2 === 180°//1
        @test 2m // 5s == (2//5)*(m/s)       # Units propagate through rationals
        @test (2//3)*m // 5 == (2//15)*m     # Quantity // Real
        @test 5.0m // s === 5.0m/s           # Quantity // Unit. Just pass units through
        @test s//(5m) === (1//5)*s/m         # Unit // Quantity. Will fail if denom is float
        @test (m//2) === 1//2 * m            # Unit // Real
        @test (2//m) === (2//1) / m          # Real // Unit
        @test (m//s) === m/s                 # Unit // Unit
        @test m / missing === missing        # Unit / missing
        @test missing / m === missing        # Missing / Unit (// is not defined for Missing)
        @test @inferred(div(10m, -3cm)) === -333
        @test @inferred(div(10m, 3)) === 3m
        @test @inferred(div(10, 3m)) === 3/m
        @test @inferred(fld(10m, -3cm)) === -334
        @test @inferred(fld(10m, 3)) === 3m
        @test @inferred(fld(10, 3m)) === 3/m
        @test @inferred(cld(10m, 3)) === 4m
        @test @inferred(cld(10, 3m)) === 4/m
        @test rem(10m, -3cm) == 1.0cm
        @test mod(10m, -3cm) == -2.0cm
        @test mod(1hr+3minute+5s, 24s) == 17s
        @test mod2pi(360°) === 0°           # 2pi is 360°
        @test mod2pi(0.5pi*u"m/dm") ≈ pi    # just testing the dimensionless fallback
        @test @inferred(inv(s)) === s^-1
        @test inv(ContextUnits(m,km)) === ContextUnits(m^-1,km^-1)
        @test inv(FixedUnits(m)) === FixedUnits(m^-1)
    end
    @testset "> Exponentiation" begin
        @test @inferred(m^3/m) === m^2
        @test @inferred(𝐋^3/𝐋) === 𝐋^2
        @test @inferred(sqrt(4m^2)) === 2.0m
        @test sqrt(4m^(2//3)) === 2.0m^(1//3)
        @test @inferred(sqrt(𝐋^2)) === 𝐋
        @test @inferred(sqrt(m^2)) === m
        @test @inferred(cbrt(8m^3)) === 2.0m
        @test cbrt(8m) === 2.0m^(1//3)
        @test @inferred(cbrt(𝐋^3)) === 𝐋
        @test @inferred(cbrt(m^3)) === m
        @test (2m)^3 === 8*m^3
        @test (8m)^(1//3) === 2.0*m^(1//3)
        @test @inferred(cis(90°)) ≈ im

        # Test inferrability of integer literal powers
        _pow_m3(x) = x^-3
        _pow_0(x) = x^0
        _pow_3(x) = x^3
        _pow_2_3(x) = x^(2//3)

        @test_throws ErrorException @inferred(_pow_2_3(m))
        @test_throws ErrorException @inferred(_pow_2_3(𝐋))
        @test_throws ErrorException @inferred(_pow_2_3(1.0m))

        @test @inferred(_pow_m3(m)) == m^-3
        @test @inferred(_pow_0(m)) == NoUnits
        @test @inferred(_pow_3(m)) == m^3

        @test @inferred(_pow_m3(𝐋)) == 𝐋^-3
        @test @inferred(_pow_0(𝐋)) == NoDims
        @test @inferred(_pow_3(𝐋)) == 𝐋^3

        @test @inferred(_pow_m3(1.0m)) == 1.0m^-3
        @test @inferred(_pow_0(1.0m)) == 1.0
        @test @inferred(_pow_3(1.0m)) == 1.0m^3
    end
    @testset "> Trigonometry" begin
        @test @inferred(sin(0.0rad)) == 0.0
        @test @inferred(cos(π*rad)) == -1
        @test @inferred(tan(π*rad/4)) ≈ 1
        @test @inferred(csc(π*rad/2)) == 1
        @test @inferred(sec(0.0*rad)) == 1
        @test @inferred(cot(π*rad/4)) ≈ 1
        @test @inferred(sin(90°)) == 1
        @test @inferred(cos(0.0°)) == 1
        @test @inferred(tan(45°)) == 1
        @test @inferred(csc(90°)) == 1
        @test @inferred(sec(0°)) == 1
        @test @inferred(cot(45°)) == 1
        @test @inferred(asin(1m/1000mm)) == 90°
        @test @inferred(acos(-1mm/1000μm)) == π*rad
        @test @inferred(atan(2000sqrt(3)ms/2.0s)) == 60°
        @test @inferred(acsc(5.0Hz*0.2s)) == π/2
        @test @inferred(asec(1m/1nm)) ≈ π/2
        @test @inferred(acot(2sqrt(3)s/2000ms)) ≈ 30°

        @test @inferred(sinh(0.0rad)) == 0.0
        @test @inferred(sinh(1J/N/m) + cosh(1rad)) ≈ MathConstants.e
        @test @inferred(tanh(1m/1µm)) == 1
        @test @inferred(csch(0.0°)) == Inf
        @test @inferred(sech(0K/Ra)) == 1
        @test @inferred(coth(1e3m*mm^-1)) == 1
        @test @inferred(asinh(0.0mg/kg)) == 0
        @test @inferred(acosh(1mm/1000μm)) == 0
        @test @inferred(atanh(0W*s/J)) == 0
        @test @inferred(acsch(hr/yr * 0)) == Inf
        @test @inferred(asech(1.0m/1000.0mm)) == 0
        @test @inferred(acoth(1km/1000m)) == Inf

        @test @inferred(sinpi(rad/2)) == 1
        @test @inferred(cospi(1rad)) == -1
        @test @inferred(sinc(1rad)) === 0
        @test @inferred(cosc(1ft/3inch)) === 0.25

        @test @inferred(atan(m*sqrt(3),1m)) ≈ 60°
        @test @inferred(atan(m*sqrt(3),1.0m)) ≈ 60°
        @test @inferred(atan(m*sqrt(3),1000mm)) ≈ 60°
        @test @inferred(atan(m*sqrt(3),1e+3mm)) ≈ 60°
        @test_throws DimensionError atan(m*sqrt(3),1e+3s)
        @test @inferred(angle((3im)*V)) ≈ 90°
    end
    @testset "> Exponentials and logarithms" begin
        for f in (exp, exp10, exp2, expm1, log, log10, log1p, log2)
            @test f(1.0 * u"m/dm") ≈ f(10)
        end
    end
    @testset "> Matrix inversion" begin
        @test inv([1 1; -1 1]u"nm") ≈ [0.5 -0.5; 0.5 0.5]u"nm^-1"
    end
    @testset "> Is functions" begin
        @test_throws ErrorException isinteger(1.0m)
        @test isinteger(1.4m/mm)
        @test !isinteger(1.4mm/m)
        @test isfinite(1.0m)
        @test !isfinite(Inf*m)
        @test isnan(NaN*m)
        @test !isnan(1.0m)
        @static if VERSION ≥ v"1.7.0-DEV.119"
            @test isunordered(NaN*m)
            @test !isunordered(Inf*m)
            @test !isunordered(1.0*m)
        end
    end
    @testset "> Floating point tests" begin
        @test isapprox(1.0u"m",(1.0+eps(1.0))u"m")
        @test isapprox(1.0u"μm/m",1e-6)
        @test !isapprox(1.0u"μm/m",1e-7)
        @test !isapprox(1.0u"m",5)
        @test frexp(1.5m) == (0.75m, 1.0)
        @test unit(nextfloat(0.0m)) == m
        @test unit(prevfloat(0.0m)) == m

        # NaN behavior
        @test NaN*m != NaN*m
        @test isequal(NaN*m, NaN*m)

        @test isapprox(1.0u"m", 1.1u"m"; atol=0.2u"m")
        @test !isapprox(1.0u"m", 1.1u"m"; atol=0.05u"m")
        @test isapprox(1.0u"m", 1.1u"m"; atol=200u"mm")
        @test !isapprox(1.0u"m", 1.1u"m"; atol=50u"mm")
        @test isapprox(1.0u"m", 1.1u"m"; rtol=0.2)
        @test !isapprox(1.0u"m", 1.1u"m"; rtol=0.05)

        # Test eps
        @test eps(1.0u"s") == eps(1.0)u"s"
        @test eps(typeof(1.0u"s")) == eps(Float64)

        # Test promotion behavior
        @test !isapprox(1.0u"m", 1.0u"s")
        @test isapprox(1.0u"m", 1000.0u"mm")
        @test_throws ErrorException isapprox(1.0*FixedUnits(u"m"), 1000.0*FixedUnits(u"mm"))
    end
end

@testset "Fast mathematics" begin
    @testset "> fma and muladd" begin
        m2cm = ContextUnits(m,cm)
        m2mm = ContextUnits(m,mm)
        mm2cm = ContextUnits(mm,cm)
        cm2cm = ContextUnits(cm,cm)
        fm, fmm = FixedUnits(m), FixedUnits(mm)
        @test @inferred(fma(2.0, 3.0m, 1.0m)) === 7.0m               # llvm good
        @test @inferred(fma(2.0, 3.0m2cm, 1.0mm2cm)) === 600.1cm2cm
        @test @inferred(fma(2.0, 3.0fm, 1.0fm)) === 7.0fm
        @test @inferred(fma(2.0, 3.0m, 35mm)) === 6.035m             # llvm good
        @test @inferred(fma(2.0, 3.0m2cm, 35mm2cm)) === 603.5cm2cm
        @test_throws ErrorException fma(2.0, 3.0fm, 35fmm)  #automatic conversion prohibited
        @test @inferred(fma(2.0m, 3.0, 35mm)) === 6.035m             # llvm good
        @test @inferred(fma(2.0m2cm, 3.0, 35mm2cm)) === 603.5cm2cm
        @test @inferred(fma(2.0m, 1.0/m, 3.0)) === 5.0               # llvm good
        @test @inferred(fma(2.0m2cm, 1.0/m2mm, 3.0)) === 5.0
        @test @inferred(fma(2.0cm, 1.0/s, 3.0mm/s)) === .023m/s      # llvm good
        @test @inferred(fma(2.0cm2cm, 1.0/s, 3.0mm/s)) === 2.3cm2cm/s
        @test @inferred(fma(2m, 1/s, 3m/s)) === 5m/s                 # llvm good
        @test @inferred(fma(2, 1.0μm/m, 1)) === 1.000002             # llvm good
        @test @inferred(fma(1.0mm/m, 1.0mm/m, 1.0mm/m)) === 0.001001 # llvm good
        @test @inferred(fma(1.0mm/m, 1.0, 1.0)) ≈ 1.001              # llvm good
        @test @inferred(fma(1.0, 1.0μm/m, 1.0μm/m)) === 2.0μm/m      # llvm good
        @test @inferred(fma(2, 1.0, 1μm/m)) === 2.000001             # llvm BAD
        @test @inferred(fma(2, 1μm/m, 1mm/m)) === 501//500000    # llvm BAD
        @test @inferred(muladd(2.0, 3.0m, 1.0m)) === 7.0m
        @test @inferred(muladd(2.0, 3.0m, 35mm)) === 6.035m
        @test @inferred(muladd(2.0m, 3.0, 35mm)) === 6.035m
        @test @inferred(muladd(2.0m, 1.0/m, 3.0)) === 5.0
        @test @inferred(muladd(2.0cm, 1.0/s, 3.0mm/s)) === .023m/s
        @test @inferred(muladd(2m, 1/s, 3m/s)) === 5m/s
        @test @inferred(muladd(2, 1.0μm/m, 1)) === 1.000002
        @test @inferred(muladd(1.0mm/m, 1.0mm/m, 1.0mm/m)) === 0.001001
        @test @inferred(muladd(1.0mm/m, 1.0, 1.0)) ≈ 1.001
        @test @inferred(muladd(1.0, 1.0μm/m, 1.0μm/m)) === 2.0μm/m
        @test @inferred(muladd(2, 1.0, 1μm/m)) === 2.000001
        @test @inferred(muladd(2, 1μm/m, 1mm/m)) === 501//500000
        @test_throws DimensionError fma(2m, 1/m, 1m)
        @test_throws DimensionError fma(2, 1m, 1V)
        @test muladd(1s, 1.0mol/s, 2.0mol) === 3.0mol               # issue 138
    end
    @testset "> @fastmath" begin
        one32 = one(Float32)*m
        eps32 = eps(Float32)*m
        eps32_2 = eps32/2

        # Note: Cannot use local functions since these are not yet optimized
        fm_ieee_32(x) = x + eps32_2 + eps32_2
        fm_fast_32(x) = @fastmath x + eps32_2 + eps32_2
        @test fm_ieee_32(one32) == one32
        @test (fm_fast_32(one32) == one32 ||
            fm_fast_32(one32) == one32 + eps32 > one32)

        one64 = one(Float64)*m
        eps64 = eps(Float64)*m
        eps64_2 = eps64/2

        # Note: Cannot use local functions since these are not yet optimized
        fm_ieee_64(x) = x + eps64_2 + eps64_2
        fm_fast_64(x) = @fastmath x + eps64_2 + eps64_2
        @test fm_ieee_64(one64) == one64
        @test (fm_fast_64(one64) == one64 ||
            fm_fast_64(one64) == one64 + eps64 > one64)

        # check updating operators
        fm_ieee_64_upd(x) = (r=x; r+=eps64_2; r+=eps64_2)
        fm_fast_64_upd(x) = @fastmath (r=x; r+=eps64_2; r+=eps64_2)
        @test fm_ieee_64_upd(one64) == one64
        @test (fm_fast_64_upd(one64) == one64 ||
            fm_fast_64_upd(one64) == one64 + eps64 > one64)

        for T in (Float32, Float64, BigFloat)
            _zero = convert(T, 0)*m
            _one = convert(T, 1)*m + eps(T)*m
            _two = convert(T, 2)*m + 1m//10
            _three = convert(T, 3)*m + 1m//100

            @test isapprox((@fastmath +_two), +_two)
            @test isapprox((@fastmath -_two), -_two)
            @test isapprox((@fastmath _zero+_one+_two), _zero+_one+_two)
            @test isapprox((@fastmath _zero-_one-_two), _zero-_one-_two)
            @test isapprox((@fastmath _one*_two*_three), _one*_two*_three)
            @test isapprox((@fastmath _one/_two/_three), _one/_two/_three)
            @test isapprox((@fastmath rem(_two, _three)), rem(_two, _three))
            @test isapprox((@fastmath mod(_two, _three)), mod(_two, _three))
            @test (@fastmath cmp(_two, _two)) == cmp(_two, _two)
            @test (@fastmath cmp(_two, _three)) == cmp(_two, _three)
            @test (@fastmath cmp(_three, _two)) == cmp(_three, _two)
            @test (@fastmath _one/_zero) == convert(T, Inf)
            @test (@fastmath -_one/_zero) == -convert(T, Inf)
            @test isnan(@fastmath _zero/_zero) # must not throw

            for x in (_zero, _two, convert(T, Inf)*m, convert(T, NaN)*m)
                @test (@fastmath isfinite(x))
                @test !(@fastmath isinf(x))
                @test !(@fastmath isnan(x))
                @test !(@fastmath issubnormal(x))
            end
        end

        for T in (Complex{Float32}, Complex{Float64}, Complex{BigFloat})
            _zero = convert(T, 0)*m
            _one = convert(T, 1)*m + im*eps(real(convert(T,1)))*m
            _two = convert(T, 2)*m + im*m//10
            _three = convert(T, 3)*m + im*m//100

            @test isapprox((@fastmath +_two), +_two)
            @test isapprox((@fastmath -_two), -_two)
            @test isapprox((@fastmath _zero+_one+_two), _zero+_one+_two)
            @test isapprox((@fastmath _zero-_one-_two), _zero-_one-_two)
            @test isapprox((@fastmath _one*_two*_three), _one*_two*_three)
            @test isapprox((@fastmath _one/_two/_three), _one/_two/_three)
            @test (@fastmath _three == _two) == (_three == _two)
            @test (@fastmath _three != _two) == (_three != _two)
            @test isnan(@fastmath _one/_zero)  # must not throw
            @test isnan(@fastmath -_one/_zero) # must not throw
            @test isnan(@fastmath _zero/_zero) # must not throw

            for x in (_zero, _two, convert(T, Inf)*m, convert(T, NaN)*m)
                @test (@fastmath isfinite(x))
                @test !(@fastmath isinf(x))
                @test !(@fastmath isnan(x))
                @test !(@fastmath issubnormal(x))
            end
        end


        # real arithmetic
        for T in (Float32, Float64, BigFloat)
            half = 1m/convert(T, 2)
            third = 1m/convert(T, 3)

            for f in (:+, :-, :abs, :abs2, :conj, :inv, :sign, :sqrt)
                @test isapprox((@eval @fastmath $f($half)), (@eval $f($half)))
                @test isapprox((@eval @fastmath $f($third)), (@eval $f($third)))
            end
            for f in (:+, :-, :*, :/, :%, :(==), :!=, :<, :<=, :>, :>=,
                      :atan, :hypot, :max, :min)
                @test isapprox((@eval @fastmath $f($half, $third)),
                               (@eval $f($half, $third)))
                @test isapprox((@eval @fastmath $f($third, $half)),
                               (@eval $f($third, $half)))
            end
            for f in (:minmax,)
                @test isapprox((@eval @fastmath $f($half, $third))[1],
                               (@eval $f($half, $third))[1])
                @test isapprox((@eval @fastmath $f($half, $third))[2],
                               (@eval $f($half, $third))[2])
                @test isapprox((@eval @fastmath $f($third, $half))[1],
                               (@eval $f($third, $half))[1])
                @test isapprox((@eval @fastmath $f($third, $half))[2],
                               (@eval $f($third, $half))[2])
            end

            half = 1°/convert(T, 2)
            third = 1°/convert(T, 3)
            for f in (:cos, :sin, :tan)
                @test isapprox((@eval @fastmath $f($half)), (@eval $f($half)))
                @test isapprox((@eval @fastmath $f($third)), (@eval $f($third)))
            end
        end

        # complex arithmetic
        for T in (Complex{Float32}, Complex{Float64}, Complex{BigFloat})
            half = (1+1im)V/T(2)
            third = (1-1im)V/T(3)

            # some of these functions promote their result to double
            # precision, but we want to check equality at precision T
            rtol = Base.rtoldefault(real(T))

            for f in (:+, :-, :abs, :abs2, :conj, :inv, :sign, :sqrt)
                @test isapprox((@eval @fastmath $f($half)), (@eval $f($half)), rtol=rtol)
                @test isapprox((@eval @fastmath $f($third)), (@eval $f($third)), rtol=rtol)
            end
            for f in (:+, :-, :*, :/, :(==), :!=)
                @test isapprox((@eval @fastmath $f($half, $third)),
                               (@eval $f($half, $third)), rtol=rtol)
                @test isapprox((@eval @fastmath $f($third, $half)),
                               (@eval $f($third, $half)), rtol=rtol)
            end

            _d = 90°/T(2)
            @test isapprox((@fastmath cis(_d)), cis(_d))
        end

        # mixed real/complex arithmetic
        for T in (Float32, Float64, BigFloat)
            CT = Complex{T}
            half = 1V/T(2)
            third = 1V/T(3)
            chalf = (1+1im)V/CT(2)
            cthird = (1-1im)V/CT(3)

            for f in (:+, :-, :*, :/, :(==), :!=)
                @test isapprox((@eval @fastmath $f($chalf, $third)),
                               (@eval $f($chalf, $third)))
                @test isapprox((@eval @fastmath $f($half, $cthird)),
                               (@eval $f($half, $cthird)))
                @test isapprox((@eval @fastmath $f($cthird, $half)),
                               (@eval $f($cthird, $half)))
                @test isapprox((@eval @fastmath $f($third, $chalf)),
                               (@eval $f($third, $chalf)))
            end

            @test isapprox((@fastmath third^3), third^3)
            @test isapprox((@fastmath chalf/third), chalf/third)
            @test isapprox((@fastmath chalf^3), chalf^3)
        end
    end
end

@testset "Rounding" begin
    @test_throws ErrorException floor(3.7m)
    @test_throws ErrorException ceil(3.7m)
    @test_throws ErrorException trunc(3.7m)
    @test_throws ErrorException round(3.7m)
    @test_throws ErrorException floor(Integer, 3.7m)
    @test_throws ErrorException ceil(Integer, 3.7m)
    @test_throws ErrorException trunc(Integer, 3.7m)
    @test_throws ErrorException round(Integer, 3.7m)
    @test_throws ErrorException floor(Int, 3.7m)
    @test_throws ErrorException ceil(Int, 3.7m)
    @test_throws ErrorException trunc(Int, 3.7m)
    @test_throws ErrorException round(Int, 3.7m)
    @test floor(1.0314m/mm) === 1031.0
    @test floor(1.0314m/mm; digits=1) === 1031.4
    @test ceil(1.0314m/mm) === 1032.0
    @test ceil(1.0314m/mm; digits=1) === 1031.4
    @test trunc(-1.0314m/mm) === -1031.0
    @test trunc(-1.0314m/mm; digits=1) === -1031.4
    @test round(1.0314m/mm) === 1031.0
    @test round(1.0314m/mm; digits=1) === 1031.4
    @test floor(Integer, 1.0314m/mm) === Integer(1031.0)
    @test ceil(Integer, 1.0314m/mm) === Integer(1032.0)
    @test trunc(Integer, -1.0314m/mm) === Integer(-1031.0)
    @test round(Integer, 1.0314m/mm) === Integer(1031.0)
    @test floor(Int16, 1.0314m/mm) === Int16(1031.0)
    @test ceil(Int16, 1.0314m/mm) === Int16(1032.0)
    @test trunc(Int16, -1.0314m/mm) === Int16(-1031.0)
    @test round(Int16, 1.0314m/mm) === Int16(1031.0)
    @test floor(typeof(1mm), 1.0314m) === 1031mm
    @test floor(typeof(1.0mm), 1.0314m) === 1031.0mm
    @test floor(typeof(1.0mm), 1.0314m; digits=1) === 1031.4mm
    @test ceil(typeof(1mm), 1.0314m) === 1032mm
    @test ceil(typeof(1.0mm), 1.0314m) === 1032.0mm
    @test ceil(typeof(1.0mm), 1.0314m; digits=1) === 1031.4mm
    @test trunc(typeof(1mm), -1.0314m) === -1031mm
    @test trunc(typeof(1.0mm), -1.0314m) === -1031.0mm
    @test trunc(typeof(1.0mm), -1.0314m; digits=1) === -1031.4mm
    @test round(typeof(1mm), 1.0314m) === 1031mm
    @test round(typeof(1.0mm), 1.0314m) === 1031.0mm
    @test round(typeof(1.0mm), 1.0314m; digits=1) === 1031.4mm
    @test round(u"inch", 1.0314m) === 41.0u"inch"
    @test round(Int, u"inch", 1.0314m) === 41u"inch"
    @test round(typeof(1m), 137cm) === 1m
    @test round(137cm/m) === 1//1
    @test round(u"m", -125u"cm", sigdigits=2) === -1.2u"m"
    @test round(u"m", (125//1)u"cm", sigdigits=2) === 1.2u"m"
    @test round(u"m", -125u"cm", RoundNearestTiesUp, sigdigits=2) === -1.2u"m"
    @test round(u"m", (125//1)u"cm", RoundNearestTiesUp, sigdigits=2) === 1.3u"m"
    @test floor(u"m", -125u"cm", sigdigits=2) === -1.3u"m"
    @test floor(u"m", (125//1)u"cm", sigdigits=2) === 1.2u"m"
    @test ceil(u"m", -125u"cm", sigdigits=2) === -1.2u"m"
    @test ceil(u"m", (125//1)u"cm", sigdigits=2) === 1.3u"m"
    @test trunc(u"m", -125u"cm", sigdigits=2) === -1.2u"m"
    @test trunc(u"m", (125//1)u"cm", sigdigits=2) === 1.2u"m"
end

@testset "Sgn, abs, &c." begin
    @test @inferred(abs(3V+4V*im)) == 5V
    @test @inferred(norm(3V+4V*im)) == 5V
    @test @inferred(abs2(3V+4V*im)) == 25V^2
    @test @inferred(abs(-3m)) == 3m
    @test @inferred(abs2(-3m)) == 9m^2
    @test @inferred(sign(-3.3m)) == -1.0
    @test @inferred(signbit(0.0m)) == false
    @test @inferred(signbit(-0.0m)) == true
    @test @inferred(copysign(3.0m, -4.0s)) == -3.0m
    @test @inferred(copysign(3.0m, 4)) == 3.0m
    @test @inferred(flipsign(3.0m, -4)) == -3.0m
    @test @inferred(flipsign(-3.0m, -4)) == 3.0m
    @test @inferred(real(3m)) == 3.0m
    @test @inferred(real((3+4im)V)) == 3V
    @test @inferred(imag(3m)) == 0m
    @test @inferred(imag((3+4im)V)) == 4V
    @test @inferred(conj(3m)) == 3m
    @test @inferred(conj((3+4im)V)) == (3-4im)V
    @test @inferred(typemin(1.0m)) == -Inf*m
    @test @inferred(typemax(typeof(1.0m))) == Inf*m
    @test @inferred(typemin(0x01*m)) == 0x00*m
    @test @inferred(typemax(typeof(0x01*m))) == 0xff*m
    @test @inferred(rand(typeof(1u"m"))) isa typeof(1u"m")
    @test @inferred(rand(MersenneTwister(0), typeof(1u"m"))) isa typeof(1u"m")
end

@testset "Collections" begin
    @testset "> Ranges" begin
        @testset ">> Some of test/ranges.jl, with units" begin
            @test @inferred(size(10m:1m:0m)) == (0,)
            @test length(1m:.2m:2m) == 6
            @test length(1.0m:.2m:2.0m) == 6
            @test length(2m:-.2m:1m) == 6
            @test length(2.0m:-.2m:1.0m) == 6
            @test @inferred(length(2m:.2m:1m)) == 0
            @test length(2.0m:.2m:1.0m) == 0

            @test length(1m:2m:0m) == 0
            L32 = range(Int32(1)*m, stop=Int32(4)*m, length=4)
            L64 = range(Int64(1)*m, stop=Int64(4)*m, length=4)
            @test L32[1] == 1m && L64[1] == 1m
            @test L32[2] == 2m && L64[2] == 2m
            @test L32[3] == 3m && L64[3] == 3m
            @test L32[4] == 4m && L64[4] == 4m

            r = 5m:-1m:1m
            @test @inferred(r[1])==5m
            @test r[2]==4m
            @test r[3]==3m
            @test r[4]==2m
            @test r[5]==1m

            @test length(.1m:.1m:.3m) == 3
            # @test length(1.1m:1.1m:3.3m) == 3
            @test @inferred(length(1.1m:1.3m:3m)) == 2
            @test length(1m:1m:1.8m) == 1

            @test (1m:2m:13m)[2:6] == 3m:2m:11m
            @test typeof((1m:2m:13m)[2:6]) == typeof(3m:2m:11m)
            @test (1m:2m:13m)[2:3:7] == 3m:6m:13m
            @test typeof((1m:2m:13m)[2:3:7]) == typeof(3m:6m:13m)
        end
        @testset ">> StepRange" begin
            r = @inferred(colon(1m, 1m, 5m))
            @test isa(r, StepRange{typeof(1m)})
            @test @inferred(length(r)) === 5
            @test @inferred(step(r)) === 1m
            @test @inferred(first(range(1mm, step=2m, length=4))) === 1mm
            @test @inferred(step(range(1mm, step=2m, length=4))) === 2000mm
            @test @inferred(last(range(1mm, step=2m, length=4))) === 6001mm
            @test @inferred(first(range(1m, step=2mm, length=4))) === (1//1)m
            @test @inferred(step(range(1m, step=2mm, length=4))) === (1//500)m
            @test @inferred(last(range(1m, step=2mm, length=4))) === (503//500)m
            @test_throws DimensionError(1m, 2V) range(1m, step=2V, length=5)
            @test_throws ArgumentError 1m:0m:5m
        end
        @testset ">> StepRangeLen" begin
            @test isa(@inferred(colon(1.0m, 1m, 5m)), StepRangeLen{typeof(1.0m)})
            @test @inferred(length(1.0m:1m:5m)) === 5
            @test @inferred(step(1.0m:1m:5m)) === 1.0m
            @test @inferred(length(0:10°:360°)) == 37 # issue 111
            @test @inferred(length(0.0:10°:2pi)) == 37 # issue 111 fallout
            @test @inferred(last(0°:0.1:360°)) === 6.2 # issue 111 fallout
            @test @inferred(first(range(1mm, step=0.1mm, length=50))) === 1.0mm # issue 111
            @test @inferred(step(range(1mm, step=0.1mm, length=50))) === 0.1mm # issue 111
            @test @inferred(last(range(0, step=10°, length=37))) == 2pi
            @test @inferred(last(range(0°, step=2pi/36, length=37))) == 2pi
            @test step(range(1.0m, step=1m, length=5)) === 1.0m
            @test @inferred(first(range(1.0mm, step=2.0m, length=4))) === 1.0mm
            @test @inferred(step(range(1.0mm, step=2.0m, length=4))) === 2000.0mm
            @test @inferred(last(range(1.0mm, step=2.0m, length=4))) === 6001.0mm
            @test @inferred(first(range(1.0m, step=2.0mm, length=4))) === 1.0m
            @test @inferred(step(range(1.0m, step=2.0mm, length=4))) === 0.002m
            @test @inferred(last(range(1.0m, step=2.0mm, length=4))) === 1.006m
            @test_throws DimensionError range(1.0m, step=1.0V, length=5)
            @test_throws ArgumentError 1.0m:0.0m:5.0m
            @test (-2.0Hz:1.0Hz:2.0Hz)/1.0Hz == -2.0:1.0:2.0  # issue 160
            @test (range(0, stop=2, length=5) * u"°")[2:end] ==
                range(0.5, stop=2, length=4) * u"°"  # issue 241
        end
        @testset ">> LinSpace" begin
            # Not using Compat.range for these because kw args don't infer in julia 0.6.2
            @test isa(@inferred(range(1.0m, stop=3.0m, length=5)),
                StepRangeLen{typeof(1.0m), Base.TwicePrecision{typeof(1.0m)}})
            @test isa(@inferred(range(1.0m, stop=10m, length=5)),
                StepRangeLen{typeof(1.0m), Base.TwicePrecision{typeof(1.0m)}})
            @test isa(@inferred(range(1m, stop=10.0m, length=5)),
                StepRangeLen{typeof(1.0m), Base.TwicePrecision{typeof(1.0m)}})
            @test isa(@inferred(range(1m, stop=10m, length=5)),
                StepRangeLen{typeof(1.0m), Base.TwicePrecision{typeof(1.0m)}})
            @test_throws Unitful.DimensionError range(1m, stop=10, length=5)
            @test_throws Unitful.DimensionError range(1, stop=10m, length=5)
            r = range(1m, stop=3m, length=3)
            @test r[1:2:end] == range(1m, stop=3m, length=2)
        end
        @testset ">> Range → Array" begin
            @test isa(collect(1m:1m:5m), Array{typeof(1m),1})
            @test isa(collect(1m:2m:10m), Array{typeof(1m),1})
            @test isa(collect(1.0m:2m:10m), Array{typeof(1.0m),1})
            @test isa(collect(range(1.0m, stop=10.0m, length=5)),
                Array{typeof(1.0m),1})
        end
        @testset ">> unit multiplication" begin
            @test @inferred((1:5)*mm) === 1mm:1mm:5mm
            @test @inferred(mm*(1:5)) === 1mm:1mm:5mm
            @test @inferred((1:2:5)*mm) === 1mm:2mm:5mm
            @test @inferred((1.0:2.0:5.01)*mm) === 1.0mm:2.0mm:5.0mm
            r = @inferred(range(0.1, step=0.1, length=3) * 1.0s)
            @test r[3] === 0.3s
            @test *(1:5, mm, s^-1) === 1mm*s^-1:1mm*s^-1:5mm*s^-1
            @test *(1:5, mm, s^-1, mol^-1) === 1mm*s^-1*mol^-1:1mm*s^-1*mol^-1:5mm*s^-1*mol^-1
        end
    end
    @testset "> Arrays" begin
        @testset ">> Array multiplication" begin
            # Quantity, quantity
            @test @inferred([1m, 2m]' * [3m, 4m])    == 11m^2
            @test @inferred([1m, 2m]' * [3/m, 4/m])  == 11
            @test typeof([1m, 2m]' * [3/m, 4/m])     == Int
            @test typeof([1m, 2V]' * [3/m, 4/V])     == Int
            @test @inferred([1V,2V]*[0.1/m, 0.4/m]') == [0.1V/m 0.4V/m; 0.2V/m 0.8V/m]

            # Probably broken as soon as we stopped using custom promote_op methods
            @test_broken @inferred([1m, 2V]' * [3/m, 4/V])  == [11]
            @test_broken @inferred([1m, 2V] * [3/m, 4/V]') ==
                [3 4u"m*V^-1"; 6u"V*m^-1" 8]

            # Quantity, number or vice versa
            @test @inferred([1 2] * [3m,4m])         == [11m]
            @test typeof([1 2] * [3m,4m])            == Array{typeof(1u"m"),1}
            @test @inferred([3m 4m] * [1,2])         == [11m]
            @test typeof([3m 4m] * [1,2])            == Array{typeof(1u"m"),1}

            @test @inferred([1,2] * [3m,4m]')    == [3m 4m; 6m 8m]
            @test typeof([1,2] * [3m,4m]')       == Array{typeof(1u"m"),2}
            @test @inferred([3m,4m] * [1,2]')    == [3m 6m; 4m 8m]
            @test typeof([3m,4m] * [1,2]')       == Array{typeof(1u"m"),2}

            # re-allow vector*(1-row matrix), PR 20423
            @test @inferred([1,2] * [3m 4m])     == [3m 4m; 6m 8m]
            @test typeof([1,2] * [3m 4m])        == Array{typeof(1u"m"),2}
            @test @inferred([3m,4m] * [1 2])     == [3m 6m; 4m 8m]
            @test typeof([3m,4m] * [1 2])        == Array{typeof(1u"m"),2}
        end
        @testset ">> Element-wise multiplication" begin
            @test @inferred([1m, 2m, 3m] * 5)            == [5m, 10m, 15m]
            @test typeof([1m, 2m, 3m] * 5)               == Array{typeof(1u"m"),1}
            @test @inferred([1m, 2m, 3m] .* 5m)          == [5m^2, 10m^2, 15m^2]
            @test typeof([1m, 2m, 3m] * 5m)              == Array{typeof(1u"m^2"),1}
            @test @inferred(5m .* [1m, 2m, 3m])          == [5m^2, 10m^2, 15m^2]
            @test typeof(5m .* [1m, 2m, 3m])             == Array{typeof(1u"m^2"),1}
            @test @inferred(Matrix{Float64}(I, 2, 2)*V)  == [1.0V 0.0V; 0.0V 1.0V]
            @test @inferred(V*Matrix{Float64}(I, 2, 2))  == [1.0V 0.0V; 0.0V 1.0V]
            @test @inferred(Matrix{Float64}(I, 2, 2).*V) == [1.0V 0.0V; 0.0V 1.0V]
            @test @inferred(V.*Matrix{Float64}(I, 2, 2)) == [1.0V 0.0V; 0.0V 1.0V]
            @test @inferred([1V 2V; 0V 3V].*2)           == [2V 4V; 0V 6V]
            @test @inferred([1V, 2V] .* [true, false])   == [1V, 0V]
            @test @inferred([1.0m, 2.0m] ./ 3)           == [1m/3, 2m/3]
            @test @inferred([1V, 2.0V] ./ [3m, 4m])      == [1V/(3m), 0.5V/m]

            @test @inferred([1, 2]kg)                  == [1, 2] * kg
            @test @inferred([1, 2]kg .* [2, 3]kg^-1)   == [2, 6]
            @test @inferred([1, 2]/kg)                 == [1/kg, 2/kg]
        end
        @testset ">> Array addition" begin
            @test @inferred([1m, 2m] + [3m, 4m])     == [4m, 6m]
            @test @inferred([1m, 2m] + [1m, 1cm])    == [2m, 201m//100]
            @test @inferred([1m] + [1cm])            == [(101//100)*m]

            # issue 127
            b = [0.0, 0.0m]
            @test b + b == b
            @test b .+ b == b
            @test eltype(b+b) === Quantity{Float64}

            # Dimensionless quantities
            @test @inferred([1mm/m] + [1.0cm/m])     == [0.011]
            @test typeof([1mm/m] + [1.0cm/m])        == Array{Float64,1}
            @test @inferred([1mm/m] + [1cm/m])       == [11//1000]
            @test typeof([1mm/m] + [1cm/m])          == Array{Rational{Int},1}
            @test @inferred([1mm/m] + [2])           == [2001//1000]
            @test typeof([1mm/m] + [2])              == Array{Rational{Int},1}
            @test_throws DimensionError [1m] + [2V]
            @test_throws DimensionError [1] + [1m]
        end
        @testset ">> Element-wise addition" begin
            @test @inferred(5m .+ [1m, 2m, 3m])      == [6m, 7m, 8m]
        end
        @testset ">> Element-wise comparison" begin
            @test @inferred([0.0m, 2.0m] .< [3.0m, 2.0μm]) == BitArray([true,false])
            @test @inferred([0.0m, 2.0m] .> [3.0m, 2.0μm]) == BitArray([false,true])
            @test @inferred([0.0m, 0.0μm] .<= [0.0mm, 0.0mm]) == BitArray([true, true])
            @test @inferred([0.0m, 0.0μm] .>= [0.0mm, 0.0mm]) == BitArray([true, true])
            @test @inferred([0.0m, 0.0μm] .== [0.0mm, 0.0mm]) == BitArray([true, true])

            # Want to make sure we play nicely with StaticArrays
            for j in (<, <=, >, >=, ==)
                @test @inferred(Base.promote_op(j, typeof(1.0m), typeof(1.0μm))) == Bool
            end
        end
        @testset ">> isapprox on arrays" begin
            @test !isapprox([1.0m], [1.0V])
            @test isapprox([1.0μm/m], [1e-6])
            @test isapprox([1cm, 200cm], [0.01m, 2.0m])
            @test !isapprox([1.0], [1.0m])
            @test !isapprox([1.0m], [1.0])
        end
        @testset ">> Unit stripping" begin
            @test @inferred(ustrip([1u"m", 2u"m"])) == [1,2]
            @test_deprecated ustrip([1,2])
            @test ustrip.([1,2]) == [1,2]
            @test typeof(ustrip([1u"m", 2u"m"])) <: Base.ReinterpretArray{Int,1}

            # With target type
            @test @inferred(ustrip(u"m", [1, 2]u"m")) == [1,2]
            @test @inferred(ustrip(u"km", [1, 2]u"m")) == [1//1000, 2//1000]
            @test typeof(ustrip(u"m", [1, 2]u"m")) <: Base.ReinterpretArray{Int,1}
            @test typeof(ustrip(u"m/ms", [1, 2]*(u"km/s"))) <: Base.ReinterpretArray{Int,1}

            # Structured matrices
            @test typeof(ustrip(adjoint([1,2]u"m"))) <: Adjoint{Int}
            @test typeof(ustrip(transpose([1 2; 3 4]u"m"))) <: Transpose{Int}
            @test typeof(ustrip(Diagonal([1,2]u"m"))) <: Diagonal{Int}
            @test typeof(ustrip(Bidiagonal([1,2,3]u"m", [1,2]u"m", :U))) <:
                Bidiagonal{Int}
            @test typeof(ustrip(Tridiagonal([1,2]u"m", [3,4,5]u"m", [6,7]u"m"))) <:
                Tridiagonal{Int}
            @test typeof(ustrip(SymTridiagonal([1,2,3]u"m", [4,5]u"m"))) <:
                SymTridiagonal{Int}

            @test typeof(ustrip(u"m", adjoint([1,2]u"m"))) <: Adjoint{Int}
            @test typeof(ustrip(u"m", Diagonal([1,2]u"m"))) <: Diagonal{Int}
        end
        @testset ">> Linear algebra" begin
            @test istril([1 1; 0 1]u"m") == false
            @test istriu([1 1; 0 1]u"m") == true
        end

        @testset ">> Array initialization" begin
            Q = typeof(1u"m")
            @test @inferred(zeros(Q, 2)) == [0, 0]u"m"
            @test @inferred(zeros(Q, (2,))) == [0, 0]u"m"
            @test @inferred(zeros(Q)[]) == 0u"m"
            @test @inferred(fill!(similar([1.0, 2.0, 3.0], Q), zero(Q))) == [0, 0, 0]u"m"
            @test @inferred(ones(Q, 2)) == [1, 1]u"m"
            @test @inferred(ones(Q, (2,))) == [1, 1]u"m"
            @test @inferred(ones(Q)[]) == 1u"m"
            @test @inferred(fill!(similar([1.0, 2.0, 3.0], Q), oneunit(Q))) == [1, 1, 1]u"m"
            @test size(rand(Q, 2)) == (2,)
            @test size(rand(Q, 2, 3)) == (2,3)
            @test eltype(@inferred(rand(Q, 2))) == Q
            @test_throws ArgumentError zero([1u"m", 1u"s"])
            @test zero(Quantity{Int,𝐋}[1u"m", 1u"mm"]) == [0, 0]u"m"
        end
    end
end

@testset "Display" begin
    withenv("UNITFUL_FANCY_EXPONENTS" => false) do
        @static if VERSION ≥ v"1.6.0-DEV.770"
            @test string(typeof(1.0m/s)) ==
                "Quantity{Float64, 𝐋 𝐓^-1, FreeUnits{(m, s^-1), 𝐋 𝐓^-1, nothing}}"
            @test string(typeof(m/s)) ==
                "FreeUnits{(m, s^-1), 𝐋 𝐓^-1, nothing}"
        else
            @test string(typeof(1.0m/s)) ==
                "Quantity{Float64,𝐋 𝐓^-1,FreeUnits{(m, s^-1),𝐋 𝐓^-1,nothing}}"
            @test string(typeof(m/s)) ==
                "FreeUnits{(m, s^-1),𝐋 𝐓^-1,nothing}"
        end
        @test string(dimension(1u"m/s")) == "𝐋 𝐓^-1"
        @test string(NoDims) == "NoDims"
    end
end

struct Foo <: Number end
Base.show(io::IO, x::Foo) = print(io, "1")
Base.show(io::IO, ::MIME"text/plain", ::Foo) = print(io, "42.0")

@testset "Show quantities" begin
    withenv("UNITFUL_FANCY_EXPONENTS" => false) do
        @test repr(1.0 * u"m * s * kg^-1") == "1.0 m s kg^-1"
        @test repr("text/plain", 1.0 * u"m * s * kg^-1") == "1.0 m s kg^-1"
        @test repr(Foo() * u"m * s * kg^-1") == "1 m s kg^-1"
        @test repr("text/plain", Foo() * u"m * s * kg^-1") == "42.0 m s kg^-1"

        # Complex quantities
        @test repr((1+2im) * u"m/s") == "(1 + 2im) m s^-1"
        @test repr("text/plain", (1+2im) * u"m/s") == "(1 + 2im) m s^-1"

        # Angular degree printing #253
        @test sprint(show, 1.0°)       == "1.0°"
        @test repr("text/plain", 1.0°) == "1.0°"

        # Concise printing of ranges
        @test repr((1:10)*u"kg/m^3") == "(1:10) kg m^-3"
        @test repr((1.0:0.1:10.0)*u"kg/m^3") == "(1.0:0.1:10.0) kg m^-3"
        @test repr((1:10)*°) == "(1:10)°"
    end
    withenv("UNITFUL_FANCY_EXPONENTS" => true) do
        @test repr(1.0 * u"m * s * kg^(-1//2)") == "1.0 m s kg⁻¹ᐟ²"
    end
    withenv("UNITFUL_FANCY_EXPONENTS" => nothing) do
        @test repr(1.0 * u"m * s * kg^(-1//2)") ==
            (Sys.isapple() ? "1.0 m s kg⁻¹ᐟ²" : "1.0 m s kg^-1/2")
    end
end

@testset "DimensionError message" begin
    function errorstr(e)
        b = IOBuffer()
        Base.showerror(b,e)
        String(take!(b))
    end
    @test errorstr(DimensionError(1u"m",2)) ==
        "DimensionError: 1 m and 2 are not dimensionally compatible."
    @test errorstr(DimensionError(1u"m",NoDims)) ==
        "DimensionError: 1 m and NoDims are not dimensionally compatible."
    @test errorstr(DimensionError(u"m",2)) ==
        "DimensionError: m and 2 are not dimensionally compatible."
end

@testset "Logarithmic quantities" begin
    @testset "> Explicit construction" begin
        @testset ">> Level" begin
            # Outer constructor
            @test Level{Decibel,1}(2) isa Level{Decibel,1,Int}
            @test_throws DimensionError Level{Decibel,1}(2V)
            @test_throws DimensionError Level{Decibel,1V}(2)
            @test_throws InexactError Level{Decibel,1}(1+1im)
            @test_throws DomainError Level{Decibel,1+0im}(2)
            @test_throws DomainError Level{Decibel,(1+0im)V}(2V)
            @test_throws DomainError Level{Decibel,(1+1im)V}(2V)

            # Inner constructor
            @test Level{Decibel,1,Int}(2) === Level{Decibel,1}(2)
            @test Level{Decibel,1,Int}(2) === Level{Decibel,1,Int}(2.0)
            @test Level{Decibel,1,Int}(2) === Level{Decibel,1,Int}(2+0im)
            @test_throws DimensionError Level{Decibel,1,typeof(2V)}(2V)
            @test_throws DimensionError Level{Decibel,1V,Int}(2)
            @test_throws TypeError Level{Decibel,1,Complex{Int}}(1+1im)
            @test_throws TypeError Level{Decibel,1V,typeof((1+1im)V)}((1+1im)V)
            @test_throws DomainError Level{Decibel,1+0im,Int}(2)
            @test_throws DomainError Level{Decibel,(1+0im)V,typeof(2V)}(2V)
        end

        @testset ">> Gain" begin
            @test Gain{Decibel}(1) isa Gain{Decibel,:?,Int}
            @test Gain{Decibel,:rp}(1) isa Gain{Decibel,:rp,Int}
            @test_throws MethodError Gain{Decibel}(1V)
            @test_throws MethodError Gain{Decibel,:?}(1V)
            @test_throws TypeError Gain{Decibel,:?,typeof(1V)}(1V)
        end
    end

    @testset "> Implicit construction" begin
        @testset ">> Level" begin
            @test 20*dBm == (@dB 100mW/mW) == (@dB 100mW/1mW) == dB(100mW,mW) == dB(100mW,1mW)
            @test 20*dBV == (@dB 10V/V) == (@dB 10V/1V) == dB(10V,V) == dB(10V,1V)
            @test_throws ArgumentError @dB 10V/V true
        end

        @testset ">> Gain" begin
            @test_throws LoadError @eval @dB 10
            @test 20*dB === dB*20
        end

        @testset ">> MixedUnits" begin
            @test dBm === MixedUnits{Level{Decibel, 1mW}}()
            @test dBm/Hz === MixedUnits{Level{Decibel, 1mW}}(Hz^-1)
        end
    end

    @testset "> Unit and dimensional analysis" begin
        @testset ">> Level" begin
            @test dimension(1dBm) === dimension(1mW)
            @test dimension(typeof(1dBm)) === dimension(1mW)
            @test dimension(1dBV) === dimension(1V)
            @test dimension(typeof(1dBV)) === dimension(1V)
            @test dimension(1dB) === NoDims
            @test dimension(typeof(1dB)) === NoDims
            @test dimension(@dB 3V/2.14V) === dimension(1V)
            @test dimension(typeof(@dB 3V/2.14V)) === dimension(1V)
            @test logunit(1dBm) === dBm
            @test logunit(typeof(1dBm)) === dBm
        end

        @testset ">> Gain" begin
            @test logunit(3dB) === dB
            @test logunit(3dB_rp) === dB_rp
            @test logunit(3dB_p) === dB_p
            @test logunit(typeof(3dB)) === dB
            @test logunit(typeof(3dB_rp)) === dB_rp
            @test logunit(typeof(3dB_p)) === dB_p
        end

        @testset ">> Quantity{<:Level}" begin
            @test dimension(1dBm/Hz) === dimension(1mW/Hz)
            @test dimension(typeof(1dBm/Hz)) === dimension(1mW/Hz)
            @test dimension(1dB/Hz) === dimension(Hz^-1)
            @test dimension(typeof(1dB/Hz)) === dimension(Hz^-1)
            @test dimension((@dB 3V/2.14V)/Hz) === dimension(1V/Hz)
            @test dimension(typeof((@dB 3V/2.14V)/Hz)) === dimension(1V/Hz)
        end
    end

    @testset "> Conversion" begin
        @test float(3dB) == 3.0dB
        @test float(@dB 3V/1V) === @dB 3.0V/1V

        @test uconvert(V, (@dB 3V/2.14V)) === 3V
        @test uconvert(V, (@dB 3V/1V)) === 3V
        @test uconvert(mW/Hz, 0dBm/Hz) == 1mW/Hz
        @test uconvert(mW/Hz, (@dB 1mW/mW)/Hz) === 1mW/Hz
        @test uconvert(dB, 1Np) ≈ 8.685889638065037dB
        @test uconvert(dB, 10dB*m/mm) == 10000dB

        @test convert(typeof(1.0dB), 1Np) ≈ 8.685889638065037dB
        @test convert(typeof(1.0dBm), 1W) == 30.0dBm
        @test_throws DimensionError convert(typeof(1.0dBm), 1V)
        @test convert(typeof(3dB), 3dB) === 3dB
        @test convert(typeof(3.0dB), 3dB) === 3.0dB
        @test convert(Float64, 0u"dBFS") === 1.0

        @test_throws ErrorException convert(Float64, u"10dB")
        @test convert(Float64, u"10dB_p") === 10.0
        @test convert(Float64, u"20dB_rp") === 10.0

        @test isapprox(uconvertrp(NoUnits, 6.02dB), 2.0, atol=0.001)
        @test uconvertrp(NoUnits, 1Np) ≈ MathConstants.e
        @test uconvertrp(Np, MathConstants.e) == 1Np
        @test uconvertrp(NoUnits, 1) == 1
        @test uconvertrp(NoUnits, 20dB) == 10
        @test uconvertrp(dB, 10) == 20dB
        @test isapprox(uconvertp(NoUnits, 3.01dB), 2.0, atol=0.001)
        @test uconvertp(NoUnits, 1Np) == (MathConstants.e)^2
        @test uconvertp(Np, (MathConstants.e)^2) == 1Np
        @test uconvertp(NoUnits, 1) == 1
        @test uconvertp(NoUnits, 20dB) == 100
        @test uconvertp(dB, 100) == 20dB

        @test linear(@dB(1mW/mW)/Hz) === 1mW/Hz
        @test linear(@dB(1.4V/2.8V)/s) === 1.4V/s

        @test convert(Gain{Decibel}, 0dB_rp) === 0dB_rp
        @test convert(Gain{Neper}, 10dB) ≈ 1.151292546Np
        @test convert(Gain{Decibel,:?}, 0dB_rp) === 0dB
        @test convert(Gain{Neper,:?}, 10dB) ≈ 1.151292546Np
        @test convert(Gain{Decibel,:?,Float32}, 0dB_rp) === 0.0f0*dB
        @test convert(Gain{Neper,:rp,Float32}, 10dB) === 1.1512926f0*Np_rp
    end

    @testset "> Equality" begin
        @testset ">> Level" begin
            @test big(3.0)dBm == big(3.0)dBm
            @test isequal(big(3.0)dBm, big(3.0)dBm)
            @test_broken hash(big(3.0)dBm) == hash(big(3.0)dBm)

            @test @dB(3.0V/2.0V) == @dB(3V/V)
            @test isequal(@dB(3.0V/2.0V), @dB(3V/V))
            @test_broken hash(@dB(3.0V/2.0V)) == hash(@dB(3V/V))
        end

        @testset ">> Gain" begin
            @test 3dB == (3//1)dB
            @test isequal(3dB, (3//1)dB)
            @test_broken hash(3dB) == hash((3//1)dB)

            @test big(3)dB == big(3)dB
            @test isequal(big(3)dB, big(3)dB)
            @test_broken hash(big(3)dB) == hash(big(3)dB)

            @test 0.0dB == -0.0dB
            @test !isequal(0.0dB, -0.0dB)
            @test hash(0.0dB) != hash(-0.0dB)
        end

        @test !(20dBm == 20dB)
        @test !(20dB == 20dBm)
        @test !(20dBm == 20dBV)
        @test !(20dBV == 20dBm)
    end

    @testset "> Addition and subtraction" begin
        @testset ">> Level" begin
            @test isapprox(10dBm + 10dBm, 13dBm; atol=0.02dBm)
            @test !isapprox(10dBm + 10dBm, 13dBm; atol=0.00001dBm)
            @test isapprox(13dBm, 20mW; atol = 0.1mW)
            @test @dB(10mW/mW) + 1mW === 11mW
            @test 1mW + @dB(10mW/mW) === 11mW
            @test @dB(10mW/mW) + @dB(90mW/mW) === @dB(100mW/mW)
            @test (@dB 10mW/3mW) + (@dB 11mW/2mW) === 21mW
            @test (@dB 10mW/3mW) + 2mW === 12mW
            @test (@dB 10mW/3mW) + 1W === 101u"kg*m^2/s^3"//100
            @test 20dB + 20dB == 40dB
            @test 20dB + 20.2dB == 40.2dB
            @test 1Np + 1.5Np == 2.5Np
            @test_throws DimensionError (1dBm + 1dBV)
            @test_throws DimensionError (1dBm + 1V)
        end

        @testset ">> Gain" begin
            for op in (:+, :*)
                @test @eval ($op)(20dB, 10dB)      === 30dB
                @test @eval ($op)(20dB_rp, 10dB)   === 30dB_rp
                @test @eval ($op)(20dB, 10dB_rp)   === 30dB_rp
                @test @eval ($op)(20dB_p, 10dB)    === 30dB_p
                @test @eval ($op)(20dB, 10dB_p)    === 30dB_p
                @test @eval ($op)(20dB_rp, 10dB_p) === 30dB
                @test @eval ($op)(20dB_p, 10dB_rp) === 30dB
                @test_throws ErrorException @eval ($op)(1dB, 1Np) # no promotion
                @test_throws ErrorException @eval ($op)(1dB_rp, 1Np)
            end
            for op in (:-, :/)
                @test @eval ($op)(20dB, 10dB)      === 10dB
                @test @eval ($op)(20dB_rp, 10dB)   === 10dB_rp
                @test @eval ($op)(20dB, 10dB_rp)   === 10dB_rp
                @test @eval ($op)(20dB_p, 10dB)    === 10dB_p
                @test @eval ($op)(20dB, 10dB_p)    === 10dB_p
                @test @eval ($op)(20dB_rp, 10dB_p) === 10dB
                @test @eval ($op)(20dB_p, 10dB_rp) === 10dB
                @test_throws ErrorException @eval ($op)(1dB, 1Np) # no promotion
                @test_throws ErrorException @eval ($op)(1dB_rp, 1Np)
            end
        end

        @testset ">> Level, meet Gain" begin
            for op in (:+, :*)
                @test @eval ($op)(10dBm, 30dB)    == 40dBm
                @test @eval ($op)(30dB, 10dBm)    == 40dBm
                @test @eval ($op)(10dBm, 30dB_rp) == 40dBm
                @test @eval ($op)(30dB_rp, 10dBm) == 40dBm
                @test @eval ($op)(10dBm, 30dB_p)  == 40dBm
                @test @eval ($op)(30dB_p, 10dBm)  == 40dBm
                @test @eval ($op)(0Np, 3dBm)      == 3dBm
            end
            for op in (:-, :/)
                @test @eval ($op)(10dBm, 30dB)    == -20dBm
                @test @eval ($op)(10dBm, 30dB_rp) == -20dBm
                @test @eval ($op)(10dBm, 30dB_p) == -20dBm
                @test @eval isapprox(($op)(10dBm, 1Np), 1.314dBm; atol=0.001dBm)
                @test @eval isapprox(($op)(10dBm, 1Np_rp), 1.314dBm; atol=0.001dBm)
                @test @eval isapprox(($op)(10dBm, 1Np_p), 1.314dBm; atol=0.001dBm)

                # cannot subtract Levels from Gains
                @test_throws ArgumentError @eval ($op)(10dB, 30dBm)
                @test_throws ArgumentError @eval ($op)(10dB_rp, 30dBm)
                @test_throws ArgumentError @eval ($op)(10dB_p, 30dBm)
                @test_throws ArgumentError @eval ($op)(1Np, 10dBm)
                @test_throws ArgumentError @eval ($op)(1Np_rp, 10dBm)
                @test_throws ArgumentError @eval ($op)(1Np_p, 10dBm)
            end
        end
    end

    @testset "> Multiplication and division" begin
        @testset ">> Level" begin
            @test (0dBm) * 10 == (10dBm)
            @test @dB(10V/V)*10 == 100V
            @test @dB(10V/V)/20 == 0.5V
            @test 10*@dB(10V/V) == 100V
            @test 10/@dB(10V/V) == 1V^-1
            @test (0dBm) * (1W) == 1*mW*W
            @test (1W) * (0dBm) == 1*W*mW
            @test 100*((0dBm)/s) == (20dBm)/s
            @test isapprox((3.01dBm)*(3.01dBm), 4mW^2, atol=0.01mW^2)
            @test typeof((1dBm * big"2").val.val) == BigFloat
            @test 10dBm/10Hz == 1mW/Hz
            @test 10Hz/10dBm == 1Hz/mW
            @test true*3dBm == 3dBm
            @test false*3dBm == -Inf*dBm
            @test 3dBm*true == 3dBm
            @test 3dBm*false == -Inf*dBm
            @test (0dBV)*(1Np) ≈ 8.685889638dBV
            @test dBm/5 ≈ 0.2dBm
            @test linear((@dB 3W/W)//3) === 1W//1
            @test (@dB 3W/W)//(@dB 3W/W) === 1//1
        end

        @testset ">> Gain" begin
            @test 3dB * 2 === 6dB
            @test 2 * 3dB === 6dB
            @test 3dB_rp * 2 === 6dB_rp
            @test 2 * 3dB_rp === 6dB_rp
            @test 3dB_p * 2 === 6dB_p
            @test 2 * 3dB_p === 6dB_p
            @test 20dB * 1mW == 100mW
            @test 1mW * 20dB == 100mW
            @test 3dB * 2.1 ≈ 6.3dB
            @test 3dB * false == 0*dB
            @test false * 3dB == 0*dB
            @test 1V * 20dB == 10V
            @test 20dB * 1V == 10V
            @test 10J * 10dB == 100J
            @test 10W/m^3 * 10dB == 100W/m^3
        end

        @testset ">> MixedUnits" begin
            @test_throws ArgumentError dB*dB
            @test_throws ArgumentError dB/Np
            @test dB/dB === NoUnits
            @test (dB*m)/(dB*s) === m/s
            @test m*dB === dB*m
            @test [1,2,3]u"dB" == u"dB"*[1,2,3]
        end
    end

    @testset "> Comparisons" begin
        @test 3dB < 5dB
        @test 3dBm < 5dBm
        @test_throws MethodError 3dB < 5dBm
    end

    @testset "> zero, one" begin
        @test zero(3dB) === 0dB
        @test zero(3dB_rp) === 0dB_rp
        @test zero(typeof(3dB)) === 0dB
        @test one(3dB) === 0dB
        @test one(3dB_rp) === 0dB_rp
        @test one(typeof(3dB)) === 0dB

        @test zero(3dBm) === (-Inf)*dBm
        @test zero(typeof(3dBm)) === (-Inf)*dBm
        @test one(3dBm) === 1.0
        @test one(typeof(3dBm)) === 1.0
        @test one(@dB 3mW/1mW) === 1
    end

    @testset "> Unit stripping" begin
        @test ustrip(500.0Np) === 500.0
        @test ustrip(20dB/Hz) === 20
        @test ustrip(20dB) === 20
        @test ustrip(20dB_rp) === 20
        @test ustrip(13dBm) ≈ 13
        @test ustrip(missing) === missing
    end

    @testset "> Display" begin
        withenv("UNITFUL_FANCY_EXPONENTS" => false) do
            @test repr(3u"dB/Hz") == "[3 dB] Hz^-1"
            @test repr("text/plain", 3u"dB/Hz") == "[3 dB] Hz^-1"
        end
        @test Unitful.abbr(3u"dBm") == "dBm"
        @test Unitful.abbr(@dB 3V/1.241V) == "dB (1.241 V)"
        @test string(360°) == "360°"
    end

    @testset "> Thanks for signing up for Log Facts!" begin
        @test_throws ErrorException 20dB == 100
        @test 20dBm ≈ 100mW
        @test 20dBV ≈ 10V
        @test 40dBV ≈ 100V

        # Maximum sound pressure level is a full swing of atmospheric pressure
        @test isapprox(uconvert(dBSPL, 1u"atm"), 194dBSPL, atol=0.1dBSPL)
    end
end

@testset "Output ordered by unit exponent" begin
    ordered = Unitful.sortexp(typeof(u"J*mol^-1*K^-1").parameters[1])
    @test typeof(ordered[1]) <: Unitful.Unit{:Joule,<:Any}
    @test typeof(ordered[2]) <: Unitful.Unit{:Kelvin,<:Any}
    @test typeof(ordered[3]) <: Unitful.Unit{:Mole,<:Any}

    ordered = Unitful.sortexp(typeof(u"mol*J^-1*K^-1").parameters[1])
    @test typeof(ordered[1]) <: Unitful.Unit{:Mole,<:Any}
    @test typeof(ordered[2]) <: Unitful.Unit{:Joule,<:Any}
    @test typeof(ordered[3]) <: Unitful.Unit{:Kelvin,<:Any}
end

# Test that the @u_str macro will find units in other modules.
module ShadowUnits
    using Unitful
    @unit m "m" MyMeter 1u"m" false
    Unitful.register(ShadowUnits)
end

@test (@test_logs (:warn, r"found in multiple") eval(:(typeof(u"m")))) ==
    Unitful.FreeUnits{(Unitful.Unit{:MyMeter, 𝐋}(0, 1//1),), 𝐋, nothing}

# Test that the @u_str macro will not find units in modules which are
# not loaded before the u_str invocation.
module FooUnits
    using Unitful
    @unit foo "foo" MyFoo 1u"m" false
    Unitful.register(FooUnits)
end
# NB: The following is LoadError in 1.0 but an ErrorException in 0.7.
@test_throws Exception eval(:(module ShouldUseFooUnits
                                  using Unitful
                                  foo() = 1u"foo"
                              end))
# Test that u_str works when FooUnits is correctly loaded.
module DoesUseFooUnits
    using Unitful, ..FooUnits
    foo() = 1u"foo"
end
@test DoesUseFooUnits.foo() === 1u"foo"

# Tests for unit extension modules in unit parsing
@test_throws ArgumentError uparse("foo", unit_context=Unitful)
@test uparse("foo", unit_context=FooUnits) === u"foo"
@test uparse("foo", unit_context=[Unitful, FooUnits]) === u"foo"
@test uparse("foo", unit_context=[FooUnits, Unitful]) === u"foo"

# Test for #272
module OnlyUstrImported
    import Unitful: @u_str
    u = u"m"
end
@test OnlyUstrImported.u === m

# Test to make sure user macros are working properly
module TUM
    using Unitful
    using Test

    @dimension f "f" FakeDim12345
    @derived_dimension FakeDim212345 f^2
    @refunit fu "fu" FakeUnit12345 f false
    @unit fu2 "fu2" FakeUnit212345 1fu false
end

@testset "User macros" begin
    @test typeof(TUM.f) == Unitful.Dimensions{(Unitful.Dimension{:FakeDim12345}(1//1),)}
    @test 1(TUM.fu) == 1(TUM.fu2)
    @test isa(1(TUM.fu), TUM.FakeDim12345)
    @test isa(TUM.fu, TUM.FakeDim12345Units)
    @test isa(1(TUM.fu)^2, TUM.FakeDim212345)
    @test isa(TUM.fu^2, TUM.FakeDim212345Units)
end


struct Num <: Real
   x::Float64
end
Base.:+(a::Num, b::Num) = Num(a.x + b.x)
Base.:-(a::Num, b::Num) = Num(a.x - b.x)
Base.:*(a::Num, b::Num) = Num(a.x * b.x)
Base.promote_rule(::Type{Num}, ::Type{<:Real}) = Num
Base.ArithmeticStyle(::Type{Num}) = Base.ArithmeticRounds()
Base.OrderStyle(::Type{Num}) = Base.Unordered()

@testset "Custom types" begin
    # Test that @generated functions work with Quantities + custom types (#231)
    @test uconvert(u"°C", Num(373.15)u"K") == Num(100)u"°C"
end

@testset "Traits" begin
    @testset "> ArithmeticStyle" begin
        @test Base.ArithmeticStyle(1m) === Base.ArithmeticWraps()
        @test Base.ArithmeticStyle(1.0m) === Base.ArithmeticRounds()
        @test Base.ArithmeticStyle((1//1)m) === Base.ArithmeticUnknown()
        @test Base.ArithmeticStyle(Num(1)m) === Base.ArithmeticRounds()
    end

    @testset "> OrderStyle" begin
        @test Base.OrderStyle(1m) === Base.Ordered()
        @test Base.OrderStyle((1+1im)m) === Base.Unordered()
        @test Base.OrderStyle(Num(1)m) === Base.Unordered()
    end
end

# Test precompiled Unitful extension modules
load_path = mktempdir()
load_cache_path = mktempdir()
try
    write(joinpath(load_path, "ExampleExtension.jl"),
          """
          module ExampleExtension
          using Unitful

          @unit year "year" JulianYear 365u"d" true

          function __init__()
              Unitful.register(ExampleExtension)
          end
          end
          """)
    pushfirst!(LOAD_PATH, load_path)
    pushfirst!(DEPOT_PATH, load_cache_path)
    @eval using ExampleExtension
    # Delay u"year" expansion until test time
    @eval @test uconvert(u"d", 1u"year") == 365u"d"
finally
    rm(load_path, recursive=true)
    rm(load_cache_path, recursive=true)
end
