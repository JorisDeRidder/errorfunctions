
use crate::auxilliary::*;


pub trait RealErrorFunctions {
    fn erfcx(self) -> Self;
    fn erf(self) -> Self;
    fn w_im(self) -> Self;
    fn erfi(self) -> Self;
    fn erfc(self) -> Self;
    fn dawson(self) -> Self;
}


impl RealErrorFunctions for f64 {

    // Compute erfcx(x) = exp(x^2) erfc(x) function, for real x
    // This function combines a few different ideas.
    // First, for x > 50, it uses a continued-fraction expansion (same as for the Faddeeva function,
    // but with algebraic simplifications for z=i*x).
    // Second, for 0 <= x <= 50, it uses Chebyshev polynomial approximations, but with two twists:
    //   a) It maps x to y = 4 / (4+x) in [0,1].  This simple transformation, inspired by a similar
    //      transformation in the octave-forge/specfun erfcx by Soren Hauberg, results in much faster
    //      Chebyshev convergence than other simple transformations I have examined.
    //   b) Instead of using a single Chebyshev polynomial for the entire [0,1] y interval,
    //      we break the interval up into 100 equal subintervals, with a switch/lookup table,
    //      and use much lower degree Chebyshev polynomials in each subinterval. This greatly
    //      improves performance in my tests.
    //
    // For x < 0, we use the relationship erfcx(-x) = 2 exp(x^2) - erfc(x), with the usual checks
    // for overflow etcetera.
    //
    // Performance-wise, it seems to be substantially faster than either the SLATEC DERFC function
    // (or an erfcx function derived therefrom) or Cody's CALERF function (from netlib.org/specfun),
    // while retaining near machine precision in accuracy.
    //
    fn erfcx(self) -> Self {
        if self >= 0.0 {
            if self > 50.0 {
                // continued-fraction expansion is faster
                const ISPI: f64 = 0.56418958354775628694807945156;                // 1 / sqrt(pi)
                if self > 5.0e7 {
                    // 1-term expansion, important to avoid overflow
                    return ISPI / self;
                } else {
                    //  5-term expansion (rely on compiler for CSE), simplified from:
                    //  ispi / (x+0.5/(x+1/(x+1.5/(x+2/x))))  */
                    let selfsqr = self*self;
                    return ISPI * (selfsqr * (selfsqr + 4.5) + 2.0) / (self * (selfsqr * (selfsqr + 5.0) + 3.75));
                }
            } else {
                return erfcx_y100(400.0 / (4.0 + self));
            }
        } else {
            if self < -26.7 {
                return 2.0 * Self::MAX;       // Mimicking HUGE_VAL in C++: a value larger than what can be represented with Self.
            } else {
                if self < -6.1 {
                    return 2.0 * (self * self).exp();
                } else {
                    return 2.0 * (self * self).exp() - erfcx_y100(400.0 / (4.0 - self));
                }
            }
        }
    }




    // The error functions
    fn erf(self) -> Self {
        let mx2 = -self * self;
        if mx2 < -750.0 {
            if self >= 0.0 {
                1.0
            } else {
                -1.0
            }
        } else {
            if self >= 0.0 {
                if self < 0.08 {
                    // Use Taylor series for small |x|, to avoid cancellation inaccuracy
                    //   erf(x) = 2/sqrt(pi) * x * (1 - x^2/3 + x^4/10 - x^6/42 + x^8/216 + ...)
                    self * (1.1283791670955125739 + mx2 * (0.37612638903183752464 + mx2 * (0.11283791670955125739 + mx2 * (0.026866170645131251760 + mx2 * 0.0052239776254421878422))))
                } else {
                    1.0 - mx2.exp() * self.erfcx()
                }
            } else {
                if self > -0.08 {
                    // Use the same Taylor series as above.
                    self * (1.1283791670955125739 + mx2 * (0.37612638903183752464 + mx2 * (0.11283791670955125739 + mx2 * (0.026866170645131251760 + mx2 * 0.0052239776254421878422))))
                } else {
                    mx2.exp() * (-self).erfcx() - 1.0
                }
            }
        }
    }





    // Compute a scaled Dawson integral
    //        w_im(x) = 2*Dawson(x)/sqrt(pi)
    // equivalent to the imaginary part w(x) for real x.
    // Uses methods similar to the erfcx calculation above: continued fractions for large |x|,
    // a lookup table of Chebyshev polynomials for smaller |x|, and finally a Taylor expansion for |x|<0.01.
    //
    fn w_im(self) -> Self {
        if self >= 0.0 {
            if self > 45.0 {
                // continued-fraction eselfpansion is faster
                const ISPI: f64 = 0.56418958354775628694807945156;             // 1 / sqrt(pi)
                if self > 5e7 {
                    // 1-term expansion, important to avoid overflow
                    return ISPI / self;
                } else {
                    // 5-term expansion (rely on compiler for CSE), simplified from:
                    //       ispi / (x-0.5/(x-1/(x-1.5/(x-2/x))))
                    let selfsqr = self * self;
                    return ISPI * (selfsqr * (selfsqr - 4.5) + 2.0)
                        / (self * (selfsqr * (selfsqr - 5.0) + 3.75));
                }
            } else {
                return w_im_y100(100.0 / (1.0 + self), self);
            }
        } else {
            // = -w_im(-x)
            if self < -45.0 {
                // continued-fraction expansion is faster
                const ISPI: f64 = 0.56418958354775628694807945156;             // 1 / sqrt(pi)
                if self < -5e7 {
                    // 1-term expansion, important to avoid overflow
                    return ISPI / self;
                } else {
                    // 5-term expansion (rely on compiler for CSE), simplified from:
                    //       ispi / (self-0.5/(self-1/(self-1.5/(self-2/self))))
                    let selfsqr = self * self;
                    return ISPI * (selfsqr * (selfsqr - 4.5) + 2.0)
                        / (self * (selfsqr * (selfsqr - 5.0) + 3.75));
                }
            } else {
                return -w_im_y100(100.0 / (1.0 - self), -self);
            }
        }
    }





    // The imaginary error function erfi(x) = -i erf(ix)
    
    fn erfi(self) -> Self {
        if self * self > 720.0 {
            if self > 0.0 {
                return f64::INFINITY;
            } else {
                return f64::NEG_INFINITY;
            }
        } else {
            return (self * self).exp() * self.w_im();
        }
    }



    // The complementary error function erfc(x) = 1 - erf(x)
    
    fn erfc(self) -> Self {
        if self * self > 750.0 {                                                // underflow
            if self >= 0.0 {
                return 0.0;
            } else {
                return 2.0;
            }
        } else {
            if self >= 0.0 {
                return (-self * self).exp() * self.erfcx();
            } else {
                return 2.0 - (-self * self).exp() * (-self).erfcx();
            }
        }
    }



    // Dawson's function Dawson(x) = sqrt(pi)/2  *  exp(-x^2) * erfi(x)
    
    fn dawson(self) -> Self {
        const SPI2: f64 = 0.8862269254527580136490837416705725913990;           // sqrt(pi)/2
        return SPI2 * self.w_im();
    }

} // end implementation of trait Faddeeva for type f64







#[cfg(test)]
mod tests {

    use super::*;
    use crate::auxilliary::relerr;
    use crate::complexerrorfunctions::*;
    use num::complex::Complex;
    type Complex64 = Complex<f64>;

    #[test]
    fn test_erf_real() {
        let tolerance = 1.0e-13;
        for i in 0..10000 {
            let x = 10.0_f64.powf(-300.0 + i as f64 * 600.0 / (10000.0 - 1.0));
            let computed_erf = x.erf();
            let expected_erf = Complex64::new(x, x*1.0e-20).erf().re;
            // println!( "{}: x={} --> erf(x) computed: {}  vs expected:  {}", i, x, computed_erf, expected_erf);
            let relative_error = relerr(computed_erf, expected_erf);
            assert!(relative_error < tolerance);

            let computed_erf = (-x).erf();
            let expected_erf = erf_with_relerror(Complex64::new(-x, x*1.0e-20), 0.0).re;
            // println!( "{}: x={} --> erf(x) computed: {}  vs expected:  {}", i, -x, computed_erf, expected_erf);
            let relative_error = relerr(computed_erf, expected_erf);
            assert!(relative_error < tolerance);
        }

        let x = f64::INFINITY;
        let computed_erf = x.erf();
        let expected_erf = erf_with_relerror(Complex64::new(x, 0.0), 0.0).re;
        let relative_error = relerr(computed_erf, expected_erf);
        assert!(relative_error < tolerance);

        let x = f64::NEG_INFINITY;
        let computed_erf = x.erf();
        let expected_erf = erf_with_relerror(Complex64::new(x, 0.0), 0.0).re;
        let relative_error = relerr(computed_erf, expected_erf);
        assert!(relative_error < tolerance);

        let x = f64::NAN;
        let computed_erf = x.erf();
        let expected_erf = erf_with_relerror(Complex64::new(x, 0.0), 0.0).re;
        let relative_error = relerr(computed_erf, expected_erf);
        assert!(relative_error < tolerance);
    }




    #[test]
    fn test_erfi_real() {
        let tolerance = 1.0e-13;
        for i in 0..10000 {
            let x = 10.0_f64.powf(-300.0 + i as f64 * 600.0 / (10000.0 - 1.0));
            let computed_erfi = x.erfi();
            let expected_erfi = erfi_with_relerror(Complex64::new(x, 0.0), 0.0).re;
            let relative_error = relerr(computed_erfi, expected_erfi);
            assert!(relative_error < tolerance);

            let computed_erfi = (-x).erfi();
            let expected_erfi = erfi_with_relerror(Complex64::new(-x, 0.0), 0.0).re;
            // println!( "{}: x={} --> erfi(x) computed: {}  vs expected:  {}", i, -x, computed_erfi, expected_erfi);
            let relative_error = relerr(computed_erfi, expected_erfi);
            assert!(relative_error < tolerance);
        }

        let x = f64::INFINITY;
        let computed_erfi = x.erfi();
        let expected_erfi = erfi_with_relerror(Complex64::new(x, 0.0), 0.0).re;
        let relative_error = relerr(computed_erfi, expected_erfi);
        assert!(relative_error < tolerance);

        let x = f64::NEG_INFINITY;
        let computed_erfi = x.erfi();
        let expected_erfi = erfi_with_relerror(Complex64::new(x, 0.0), 0.0).re;
        let relative_error = relerr(computed_erfi, expected_erfi);
        assert!(relative_error < tolerance);

        let x = f64::NAN;
        let computed_erfi = x.erfi();
        let expected_erfi = erfi_with_relerror(Complex64::new(x, 0.0), 0.0).re;
        let relative_error = relerr(computed_erfi, expected_erfi);
        assert!(relative_error < tolerance);
    }





    #[test]
    fn test_erfcx_real() {
        let tolerance = 1.0e-13;
        for i in 0..10000 {
            let x = 10.0_f64.powf(-300.0 + i as f64 * 600.0 / (10000.0 - 1.0));
            let computed_erfcx = x.erfcx();
            let expected_erfcx = erfcx_with_relerror(Complex64::new(x, 0.0), 0.0).re;
            let relative_error = relerr(computed_erfcx, expected_erfcx);
            if relative_error >= tolerance {
                println!( "erfcx_real: {}: x={} --> erfcx(x) computed: {}  vs expected:  {}", i, x, computed_erfcx, expected_erfcx);
            }
            assert!(relative_error < tolerance);

            let computed_erfcx = (-x).erfcx();
            let expected_erfcx = erfcx_with_relerror(Complex64::new(-x, 0.0), 0.0).re;
            let relative_error = relerr(computed_erfcx, expected_erfcx);
            if relative_error >= tolerance {
                println!( "erfcx_real: {}: x={} --> erfcx(x) computed: {}  vs expected:  {}", i, -x, computed_erfcx, expected_erfcx);
            }
            assert!(relative_error < tolerance);
        }

        let x = f64::INFINITY;
        let computed_erfcx = x.erfcx();
        let expected_erfcx = erfcx_with_relerror(Complex64::new(x, 0.0), 0.0).re;
        println!( "x={} --> erfcx(x) computed: {}  vs expected:  {}", x, computed_erfcx, expected_erfcx);
        let relative_error = relerr(computed_erfcx, expected_erfcx);
        assert!(relative_error < tolerance);

        let x = f64::NEG_INFINITY;
        let computed_erfcx = x.erfcx();
        let expected_erfcx = erfcx_with_relerror(Complex64::new(x, 0.0), 0.0).re;
        println!( "x={} --> erfcx(x) computed: {}  vs expected:  {}", x, computed_erfcx, expected_erfcx);
        let relative_error = relerr(computed_erfcx, expected_erfcx);
        assert!(relative_error < tolerance);

        let x = f64::NAN;
        let computed_erfcx = x.erfcx();
        let expected_erfcx = erfcx_with_relerror(Complex64::new(x, 0.0), 0.0).re;
        println!( "x={} --> erfcx(x) computed: {}  vs expected:  {}", x, computed_erfcx, expected_erfcx);
        let relative_error = relerr(computed_erfcx, expected_erfcx);
        assert!(relative_error < tolerance);
    }




    #[test]
    fn test_erfc_real() {
        let tolerance = 1.0e-13;
        for i in 0..10000 {
            let x = 10.0_f64.powf(-300.0 + i as f64 * 600.0 / (10000.0 - 1.0));
            let computed_erfc = x.erfc();
            let expected_erfc = erfc_with_relerror(Complex64::new(x, x*1.0e-20), 0.0).re;
            // println!( "{}: x={} --> erfc(x) computed: {}  vs expected:  {}", i, x, computed_erfc, expected_erfc);
            let relative_error = relerr(computed_erfc, expected_erfc);
            assert!(relative_error < tolerance);

            let computed_erfc = (-x).erfc();
            let expected_erfc = erfc_with_relerror(Complex64::new(-x, x*1.0e-20), 0.0).re;
            // println!( "{}: x={} --> erfc(x) computed: {}  vs expected:  {}", i, -x, computed_erfc, expected_erfc);
            let relative_error = relerr(computed_erfc, expected_erfc);
            assert!(relative_error < tolerance);
        }

        let x = f64::INFINITY;
        let computed_erfc = x.erfc();
        let expected_erfc = erfc_with_relerror(Complex64::new(x, 0.0), 0.0).re;
        let relative_error = relerr(computed_erfc, expected_erfc);
        assert!(relative_error < tolerance);

        let x = f64::NEG_INFINITY;
        let computed_erfc = x.erfc();
        let expected_erfc = erfc_with_relerror(Complex64::new(x, 0.0), 0.0).re;
        let relative_error = relerr(computed_erfc, expected_erfc);
        assert!(relative_error < tolerance);

        let x = f64::NAN;
        let computed_erfc = x.erfc();
        let expected_erfc = erfc_with_relerror(Complex64::new(x, 0.0), 0.0).re;
        let relative_error = relerr(computed_erfc, expected_erfc);
        assert!(relative_error < tolerance);
    }





    #[test]
    fn test_dawson_real() {
        let tolerance = 1.0e-13;
        for i in 0..10000 {
            let x = 10.0_f64.powf(-300.0 + i as f64 * 600.0 / (10000.0 - 1.0));
            let computed_dawson = x.dawson();
            let expected_dawson = dawson_with_relerror(Complex64::new(x, x*1.0e-20), 0.0).re;
            let relative_error = relerr(computed_dawson, expected_dawson);
            if relative_error >= tolerance {
                println!("i: {},  x={},  computed y: {},  expected y: {}", i, x, computed_dawson, expected_dawson);
            }
            assert!(relative_error < tolerance);

            let computed_dawson = (-x).dawson();
            let expected_dawson = dawson_with_relerror(Complex64::new(-x, x*1.0e-20), 0.0).re;
            let relative_error = relerr(computed_dawson, expected_dawson);
            if relative_error >= tolerance {
                println!("i: {},  x={},  computed y: {},  expected y: {}", i, x, computed_dawson, expected_dawson);
            }
            assert!(relative_error < tolerance);
        }

        let x = f64::INFINITY;
        let computed_dawson = x.dawson();
        let expected_dawson = dawson_with_relerror(Complex64::new(x, 0.0), 0.0).re;
        let relative_error = relerr(computed_dawson, expected_dawson);
        if relative_error >= tolerance {
            println!("x={},  computed y: {},  expected y: {}", x, computed_dawson, expected_dawson);
        }
        assert!(relative_error < tolerance);

        let x = f64::NEG_INFINITY;
        let computed_dawson = x.dawson();
        let expected_dawson = dawson_with_relerror(Complex64::new(x, 0.0), 0.0).re;
        let relative_error = relerr(computed_dawson, expected_dawson);
        if relative_error >= tolerance {
            println!("x={},  computed y: {},  expected y: {}", x, computed_dawson, expected_dawson);
        }
        assert!(relative_error < tolerance);

        let x = f64::NAN;
        let computed_dawson = x.dawson();
        let expected_dawson = dawson_with_relerror(Complex64::new(x, 0.0), 0.0).re;
        let relative_error = relerr(computed_dawson, expected_dawson);
        if relative_error >= tolerance {
            println!("x={},  computed y: {},  expected y: {}", x, computed_dawson, expected_dawson);
        }
        assert!(relative_error < tolerance);
    }

} 

