```@meta
DocTestSetup = quote
    using Unitful
end
```
# Logarithmic scales

!!! note
    Logarithmic scales are new to Unitful and should be considered experimental.

Unitful provides a way to use logarithmically-scaled quantities as of v0.4.0. Some
compromises have been made in striving for logarithmic quantities to be both usable and
consistent. In the following discussion, for pedagogical purposes, we will assume prior
familiarity with the definitions of `dB` and `dBm`.

## Constructing logarithmic quantities

Left- or right-multiplying a pure number by a logarithmic "unit", whether dimensionful or
dimensionless, is short-hand for constructing a logarithmic quantity.

```jldoctest
julia> 3u"dB"
3 dB

julia> 3u"dBm"
3.0 dBm

julia> u"dB"*3 === 3u"dB"
true
```

Currently implemented are `dB`, `B`, `dBm`, `dBV`, `dBu`, `dBμV`, `dBSPL`, `dBFS`, `cNp`,
`Np`.

One can also construct logarithmic quantities using the `@dB`, `@B`, `@cNp`, `@Np` macros to
use an arbitrary reference level:

```jldoctest
julia> using Unitful: mW, V

julia> @dB 10mW/mW
10.0 dBm

julia> @dB 10V/V
20.0 dBV

julia> @dB 3V/4V
-2.498774732165999 dB (4 V)

julia> @Np ℯ*V/V    # ℯ = 2.71828...
1.0 Np (1 V)
```

These macros are exported by default since empirically macros are defined less often than
variables and generic functions. When using the macros, the levels are constructed at parse
time. The scales themselves are callable as functions if you need to construct a level that
way (they are not exported):

```jldoctest
julia> using Unitful: dB, mW, V

julia> dB(10mW,mW)
10.0 dBm
```

In calculating the logarithms, the log function appropriate to the scale in question is used
(`log10` for decibels, `log` for nepers).

There is an important difference in these two approaches to constructing logarithmic
quantities. When we construct `0dBm`, the power in `mW` is calculated and stored,
resulting in a lossy floating-point conversion. This can be avoided by constructing
`0 dBm` as `@dB 1mW/mW`.

It is important to keep in mind that the reference level is just used to calculate the
logarithms, and nothing more. When there is ambiguity about what to do, we fall back
to the underlying linear quantities, paying no mind to the reference levels:

```jldoctest
julia> using Unitful: mW

julia> (@dB 10mW/1mW) + (@dB 10mW/2mW)
20 mW
```

With a few exceptions, dimensionful logarithmic units, such as `dBm` behave just like the underlying
linear unit for purposes of arithmetic (i.e. arithmetic operations commute with `linear`). However,
note that they will still be displayed logarithmically. In contrast, arithmetic on dimensionless
logarithmic units (i.e. gains/attenuations) such as `dB` behaves logarithmically. This will be explored in more detail below.

The `@dB` and `@Np` macros will fail if either a dimensionless number or a ratio of
dimensionless numbers is used. This is because the ratio could be of power quantities or of
root-power quantities, leading to ambiguities. After all, usually it is the ratio that is
dimensionless, not the numerator and denominator that make up the ratio. In some cases
it may nonetheless be convenient to have a dimensionless reference level. By providing an
extra `Bool` argument to these macros, you can explicitly choose whether the resulting ratio
should be considered a "root-power" or "power" ratio. You can only do this for dimensionless
numbers:

```jldoctest
julia> @dB 10/1 true   # is a root-power (amplitude) ratio
20.0 dBFS

julia> @dB 10/1 false  # is not a root-power ratio; is a power ratio
10.0 dB (power ratio with reference 1)
```

Note that `dBFS` is defined to represent amplitudes relative to 1 in `dB`, hence the
custom display logic.

Also, you can of course use functions instead of macros:

```jldoctest
julia> using Unitful: dB, mW

julia> dB(10,1,true)
20.0 dBFS

julia> dB(10mW,mW,true)
ERROR: ArgumentError: when passing a final Bool argument, this can only be used with dimensionless numbers.
[...]
```

### Logarithmic quantities with no reference level specified

Logarithmic quantities with no reference level specified typically represent some amount of
gain or attenuation, i.e. a ratio which is dimensionless. These can be constructed as,
for example, `10*dB`, which displays similarly (`10 dB`). The type of this kind of
logarithmic quantity is:

```@docs
    Unitful.Gain
```

One might expect that any gain / attenuation factor should be convertible to a scalar,
that is, to `x == y/z` if you had `10*log10(x)` dB. However, it turns out that in dB, a ratio
of powers is defined as `10*log10(y/z)`, but a ratio of voltages or other root-power
quantities is defined as `20*log10(y/z)`. Clearly, converting back from decibels to a scalar
is ambiguous, and so we have not implemented automatic promotion to avoid incorrect
results. You can use [`Unitful.uconvertp`](@ref) to interpret a `Gain` as a ratio of power
quantities (hence the `p` in `uconvertp`), or [`Unitful.uconvertrp`](@ref) to interpret as
a ratio of root-power (field) quantities.

### "Dimensionful" logarithmic quantities?

In this package, quantities with units like `dBm` are considered to have the dimension of
power, even though the expression `P(dBm) = 10*log10(P/1mW)` is dimensionless and formed
from a dimensionless ratio. Practically speaking, these kinds of logarithmic quantities are
fungible whenever they share the same dimensions, so it is more convenient to adopt this
convention (people refer to `dBm/Hz` as a power spectral density, etc.) Presumably, one
would like to have `10dBm isa Unitful.Power` for dispatch too. Therefore, in the following
discussion, we will shamelessly (okay, with some shame) speak of dimensionful logarithmic
quantities, or `Level`s for short:

```@docs
    Unitful.Level
```

Actually, the defining characteristic of a `Level` is that it has a reference level,
which may or may not be dimensionful. It usually is, but is not in the case of e.g. `dBFS`.

Finally, for completeness we note that both `Level` and `Gain` are subtypes of `LogScaled`:

```@docs
    Unitful.LogScaled
```

## Addition and Multiplication rules

For dimensionaless logarithmic units, addition behaves as one might expect:

```jldoctest
julia> 10u"dB" + 10u"dB"
20 dB
```

I.e. the gains add. However, as hinted at above, dimensionful logarithmic units,
behave as their corresponding linear unit:

```
julia> 10u"dBm" + 10u"dBm"
13.010299956639813 dBm

julia> linear(10u"dBm") + linear(10u"dBm")
20.0 mW

julia> uconvert(u"dBm", ans)
13.010299956639813 dBm
```

Note that this may seem strange from an arithmetic perspective, as written, but
the arithmetic is entirely consistent. It can be helpful to think of the arithmetic
as being performed on the linear units, with the logarithmic units simply being a
display hint (although the quantity being stored is indeed the displayed logarithmic
value).

Multiplication by a scalar is consistent with the addition rules above:

```jldoctest
julia> 3u"dB" * 2
6 dB

julia> 3u"dB" + 3u"dB"
6 dB

julia> 2 * 0u"dB"
0 dB

julia> 0u"dBm" * 2
3.010299956639812 dBm

julia> 0u"dBm" + 0u"dBm"
3.010299956639812 dBm
```

Logarithmic quantities can only be multiplied by scalar, linear units, or quantities,
but not logarithmic "units" or quantities.  When a logarithmic quantity is multiplied by a
linear quantity, the logarithmic quantity is linearized and multiplication proceeds as
usual:

```jldoctest
julia> (0u"dBm") * (1u"W")
1.0 mW W
```

The previous example returns a floating point value because in constructing the level
`0 dBm`, the power in `mW` is calculated and stored, entailing a floating point
conversion. This can be avoided by constructing `0 dBm` as `@dB 1mW/mW`:

```jldoctest
julia> (@dB 1u"mW"/u"mW") * (1u"W")
1 mW W
```

We refer to a quantity with both logarithmic "units" and linear units as a mixed quantity.
For mixed quantities, the numeric value associates with the logarithmic unit, and the
quantity is displayed in a way that makes this explicit:

```jldoctest
julia> (0u"dBm")/u"Hz"
[0.0 dBm] Hz^-1

julia> (0u"dB")/u"Hz"
[0 dB] Hz^-1

julia> 0u"dB/Hz"
[0 dB] Hz^-1
```

Since dimensionful logarithmic units still behave as their corresponding linear unit,
working with dimensionful units is entirely consistnet.

The behavior of addition and multiplication are summarized in the following table, with entries marked by
† indicate prohibited operations. This table is populated automatically whenever the docs
are built.

```@eval
using Latexify, Unitful
head = ["100", "20dB", "1Np", "10.0dBm", "10.0dBV", "1mW"]
side = ["+"; "**" .* head .* "**"]
quantities = uparse.(head)
tab = fill("", length(head), length(head))
for col = eachindex(head), row = 1:col
    try
        tab[row, col] = string(quantities[row] + quantities[col])
    catch
        tab[row, col] = "†"
    end
end
mdtable(tab, latex=false, head=head, side=side)
```

```@eval
using Latexify, Unitful
head = ["10", "Hz^-1", "dB", "dBm", "1/Hz", "1mW", "3dB", "3dBm"]
side = ["*"; "**" .* head .* "**"]
quantities = uparse.(head)
tab = fill("", length(head), length(head))
for col = eachindex(head), row = 1:col
    try
        tab[row, col] = string(quantities[row] * quantities[col])
    catch
        if quantities[row] === u"1/Hz" && quantities[col] === u"3dB"
            tab[row, col] = "† ‡"
        else
            tab[row, col] = "†"
        end
    end
end
mdtable(tab, latex=false, head=head, side=side)
```

‡: `1/Hz * 3dB` could be allowed, technically, but we throw an error if it's unclear whether
a quantity is a root-power or power quantity:

```jldoctest
julia> u"1/Hz" * u"3dB"
ERROR: undefined behavior. Please file an issue with the code needed to reproduce.
```

On the other hand, if it can be determined that a power quantity or root-power quantity
is being multiplied by a gain, then the gain is interpreted as a power ratio or root-power
ratio, respectively:

```jldoctest
julia> 1u"mW" * 20u"dB"
100.0 mW

julia> 1u"V" * 20u"dB"
10.0 V
```

## Mixed Arithmetic

One final question to answer is how arithmetic behaves when it involves both dimensionless
and dimensionful logarithmic units. The answer here is that in mixed arithmetic, both
dimensionless and dimensionful units are treated logarithmically. This is done for
convenience and can break commutativity and associativity, so should be probably avoided in generic
code.

```
julia> 10u"dBm" + 20u"dB"
30.0 dBm

julia> (10u"dBm" + 10u"dBm") + 20u"dB"
33.01029995663981 dBm

julia> 10u"dBm" + (10u"dBm" + 20u"dB")
30.043213737826427 dBm

julia> 10u"dBm" * 20u"dB"
ERROR: ArgumentError: Multiplying a level by a Gain is disallowed. Use addition, or `linear` depending on context.

julia> linear(10u"dBm") * 20u"dB"
1000.0 mW

julia> linear(10u"dBm") + 20u"dB"
ERROR: ArgumentError: Adding a gain to a linear quantity is disallowed. Use multiplication or convert to `Level` first
```

## Conversion

As alluded to earlier, conversions can be tricky because so-called logarithmic units are not
units in the conventional sense.

You may use [`linear`](@ref) to convert to a linear scale when you have a `Level` or
`Quantity{<:Level}` type. There is a fallback for `Number`, which just returns the number.

```jldoctest
julia> linear(@dB 10u"mW"/u"mW")
10 mW

julia> linear(20u"dBm/Hz")
100.0 mW Hz^-1

julia> linear(30u"W")
30 W

julia> linear(12)
12
```

Linearizing a `Quantity{<:Gain}` or a `Gain` to a real number is ambiguous, because the real
number may represent a ratio of powers or a ratio of root-power (field) quantities. We
implement [`Unitful.uconvertp`](@ref) and [`Unitful.uconvertrp`](@ref) which may be
thought of as disambiguated `uconvert` functions. There is a one argument version that
assumes you are converting to a unitless number. These functions can take either a `Gain`
or a `Real` so that they may be used somewhat generically.

```jldoctest
julia> uconvertrp(NoUnits, 20u"dB")
10.0

julia> uconvertp(NoUnits, 20u"dB")  
100.0

julia> uconvertp(u"dB", 100)
20.0 dB

julia> uconvertp(u"Np", ℯ^2)
1.0 Np

julia> uconvertrp(u"Np", ℯ)
1//1 Np
```

## Notation

This package displays logarithmic quantities using shorthand like `dBm` where available.
This should probably not be done in polite company. To quote "Guide for the Use of the
International System of Units (SI)," NIST Special Pub. 811 (2008):

> The rules of Ref. [5: IEC 60027-3] preclude, for example, the use of the symbol dBm to
> indicate a reference level of power of 1 mW. This restriction is based on the rule of Sec.
> 7.4, which does not permit attachments to unit symbols.

The authorities say the reference level should always specified. In practice, this hasn't
stopped the use of `dBm` and the like on commercially available test equipment. Dealing with
these units is unavoidable in practice. When no shorthand exists, we follow NIST's advice in
displaying logarithmic quantities:

> When such data are presented in a table or in a figure, the following condensed notation
> may be used instead: -0.58 Np (1 μV/m); 25 dB (20 μPa).

## Custom logarithmic scales

```@docs
    Unitful.@logscale
```

## API

```@docs
    Unitful.linear
    Unitful.reflevel
    Unitful.uconvertp
    Unitful.uconvertrp
```
