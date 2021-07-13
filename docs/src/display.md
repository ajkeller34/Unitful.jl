# How units are displayed

By default, exponents on units or dimensions are indicated using Unicode superscripts on
macOS and without superscripts on other operating systems. You can set the environment
variable `UNITFUL_FANCY_EXPONENTS` to either `true` or `false` to force using or not using
the exponents. You can also set the `:fancy_exponent` IO context property to either `true`
or `false` to force using or not using the exponents.

## Specifying precision

If you need to set precision of displayed value of your unitful quantity, you can use the
`round` function as follows:

```julia
julia> d = 1.7438748921932u"mm"
1.7438748921932 mm

julia> round(typeof(d), d; sigdigits=3)
1.74 mm
```


```@docs
Unitful.BracketStyle
Unitful.abbr
Unitful.prefix
Unitful.show(::IO, ::Quantity)
Unitful.show(::IO, ::Unitful.Unitlike)
Unitful.showrep
Unitful.showval
Unitful.superscript
```
