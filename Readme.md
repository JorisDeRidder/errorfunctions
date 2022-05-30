
This crate allows to compute:

1. The error function `erf(z)` for complex and real arguments `z`:

$$ {\rm erf}(z) = \frac{2}{\sqrt{\pi}} \int_0^z e^{-t^2} dt $$

2. The complementary error function `erfc(z)`for complex and real arguments `z`:

$$ {\rm erfc}(z) = 1 - {\rm erf}(z) $$

3. The imaginary error function `erfi(z)` for complex and real arguments `z`:

$$ {\rm erfi}(z) = -i\  {\rm erf}(i z) = \frac{2}{\sqrt{\pi}} \int_0^z e^{t^2} dt $$

4. Dawson's function `dawson(z)` for complex and real arguments `z`:

$$ {\rm dawson}(z) = \frac{\sqrt{\pi}}{2} \ e^{-z^2} \ {\rm erfi}(z) = e^{-z^2} \int_0^z e^{t^2} dt $$

5. The Faddeeva function `w(z)` for complex and real arguments `z`:

$$ {\rm w}(z) = e^{-z^2}\ {\rm erfc}(-i z) = e^{-z^2} \ \left(1 + \frac{2i}{\sqrt{\pi}} \int_0^z e^{t^2} dt \right) $$

6. The scaled complementary error function `erfcx(z)` for complex and real arguments `z`:

$$ {\rm erfcx}(z) = e^{z^2} \ {\rm erfc}(z) = {\rm w}(i z) $$ 

7. The imaginary part of the Faddeeva function `w_im(x)` for real arguments `x`:

$$ {\rm w}\\_{\rm im}(x) = {\rm Im}({\rm w}(x)) = e^{-x^2} {\rm erfi}(x) $$

The implementation of this crate is a port of Steven G. Johnson's 
[Faddeeva](http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package) C/C++ library in Rust. 
The functions are computed in an efficient way up to machine precision for `Complex<f64>` or `f64` arguments. 
The functions handle NaN's and (positive and negative) infinities as arguments correctly.



## Examples

Computing the error functions for a complex argument can be done as in the following example:

```
use num::complex::Complex;
use errorfunctions::ComplexErrorFunctions;

fn main() {

let z = Complex::<f64>::new(1.21, -0.93);

println!("z = {}", z);
println!("erf(z)    = {}",    z.erf());
println!("erfc(z)   = {}",   z.erfc());
println!("erfcx(z)  = {}",  z.erfcx());
println!("erfi(z)   = {}",   z.erfi());
println!("w(z)      = {}",      z.w());
println!("dawson(z) = {}", z.dawson());
}
```

Computing the error functions for a real argument can be done as in the following example:

```
use errorfunctions::RealErrorFunctions;

fn main() {

let x: f64 = 0.934;

println!("x = {}", x);
println!("erf(x)    = {}",    x.erf());
println!("erfc(x)   = {}",   x.erfc()); 
println!("erfcx(x)  = {}",  x.erfcx()); 
println!("erfi(x)   = {}",   x.erfi());
println!("Im(w(x))  = {}",   x.w_im()); 
println!("dawson(x) = {}", x.dawson());
}
```

If, for some reason, you don't need machine precision, you can specify the desired 
[relative error](https://en.wikipedia.org/wiki/Approximation_error) as follows:

```
use num::complex::Complex;
use errorfunctions::*;

fn main() {

let z = Complex::<f64>::new(1.21, -0.93);
let relerror = 1.0e-3;

println!("z = {}", z);
println!("erf(z)    = {}",    erf_with_relerror(z, relerror));
println!("erfc(z)   = {}",   erfc_with_relerror(z, relerror));
println!("erfcx(z)  = {}",  erfcx_with_relerror(z, relerror));
println!("erfi(z)   = {}",   erfi_with_relerror(z, relerror));
println!("w(z)      = {}",      w_with_relerror(z, relerror));
println!("dawson(z) = {}", dawson_with_relerror(z, relerror));
}
```


## Tests

The extenstive set of unit tests in the [Faddeeva](http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package) 
code was also ported to Rust and is included in this crate.


## Credits

Since this is a close to literal translation in Rust of Steven G. Johnson's code, credit should go to him.





