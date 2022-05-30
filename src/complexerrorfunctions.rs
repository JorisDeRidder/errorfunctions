
use num::complex::Complex;

use std::f64::consts::PI;
use crate::realerrorfunctions::RealErrorFunctions;
use crate::auxilliary::*;

type Complex64 = Complex<f64>;


// The Faddeeva function

pub fn w_with_relerror(z: Complex64, mut relerr: f64) -> Complex64 {
    if z.re == 0.0 {
        return Complex64::new((z.im).erfcx(), z.re); // give correct sign of 0 in w.im
    } else if z.im == 0.0 {
        return Complex64::new((-z.re * z.re).exp(), (z.re).w_im());
    }

    let a: f64;
    let a2: f64;
    let c: f64;
    if relerr <= f64::EPSILON {
        relerr = f64::EPSILON;
        a = 0.518321480430085929872; // pi / sqrt(-log(eps*0.5))
        c = 0.329973702884629072537; // (2/pi) * a;
        a2 = 0.268657157075235951582; // a^2
    } else {
        if relerr > 0.1 {
            relerr = 0.1; // not sensible to compute < 1 digit
        }
        a = PI / (-(relerr * 0.5).log2()).sqrt();
        c = (2.0 / PI) * a;
        a2 = a * a;
    }

    let x: f64 = z.re.abs();
    let y: f64 = z.im;
    let ya: f64 = y.abs();

    let ret: Complex64; // return value

    let mut sum1 = 0.0;
    let mut sum2 = 0.0;
    let mut sum3 = 0.0;
    let mut sum4 = 0.0;
    let mut sum5 = 0.0;

    // Use a continued fraction for large |z|, as it is faster.
    // As pointed out by M. Zaghloul, the continued fraction seems to give a large relative error in
    //    Re w(z) for |x| ~ 6 and small |y|, so use algorithm 816 in this region.

    if ya > 7.0 || (x > 6.0 && (ya > 0.1 || (x > 8.0 && ya > 1e-10) || x > 28.0)) {
        // Poppe & Wijers suggest using a number of terms
        //      nu = 3 + 1442 / (26*rho + 77)
        // where rho = sqrt((x/x0)^2 + (y/y0)^2) where x0=6.3, y0=4.4.
        // They only use this expansion for rho >= 1, but rho a little less than 1 seems okay too.
        // Instead, I did my own fit to a slightly different function that avoids the hypotenuse calculation, using NLopt to minimize
        // the sum of the squares of the errors in nu with the constraint that the estimated nu be >= minimum nu to attain machine precision.
        // I also separate the regions where nu == 2 and nu == 1.

        const ISPI: f64 = 0.56418958354775628694807945156; // 1 / sqrt(pi)
        let xs: f64 = if y < 0.0 { -z.re } else { z.re }; // Compute for -z if y < 0
        if x + ya > 4000.0 {
            // nu <= 2
            if x + ya > 1.0e7 {
                // nu == 1, w(z) = i/sqrt(pi) / z
                // scale to avoid overflow
                if x > ya {
                    let yax = ya / xs;
                    let denom = ISPI / (xs + yax * ya);
                    ret = Complex64::new(denom * yax, denom);
                } else {
                    if ya.is_infinite() && ya.is_sign_positive() {
                        if x.is_nan() || y < 0.0 {
                            return Complex64::new(f64::NAN, f64::NAN);
                        } else {
                            return Complex64::new(0.0, 0.0);
                        }
                    } else {
                        let xya = xs / ya;
                        let denom = ISPI / (xya * xs + ya);
                        ret = Complex64::new(denom, denom * xya);
                    }
                }
            } else {
                // nu == 2, w(z) = i/sqrt(pi) * z / (z*z - 0.5)
                let dr = xs * xs - ya * ya - 0.5;
                let di = 2.0 * xs * ya;
                let denom = ISPI / (dr * dr + di * di);
                ret = Complex64::new(denom * (xs * di - ya * dr), denom * (xs * dr + ya * di));
            }
        } else {
            // compute nu(z) estimate and do general continUed fraction
            const C0: f64 = 3.9;
            const C1: f64 = 11.398;
            const C2: f64 = 0.08254;
            const C3: f64 = 0.1421;
            const C4: f64 = 0.2023; // fit
            let nu = (C0 + C1 / (C2 * x + C3 * ya + C4)).floor();
            let mut wr = xs;
            let mut wi = ya;
            let mut nu = 0.5 * (nu - 1.0);
            while nu > 0.4 {
                // w <- z - nu/w:
                let denom = nu / (wr * wr + wi * wi);
                wr = xs - wr * denom;
                wi = ya + wi * denom;
                nu -= 0.5;
            }
            {
                // w(z) = i/sqrt(pi) / w:
                let denom = ISPI / (wr * wr + wi * wi);
                ret = Complex64::new(denom * wi, denom * wr);
            }
        }

        if y < 0.0 {
            // use w(z) = 2.0*exp(-z*z) - w(-z), but be careful of overflow in exp(-z*z) = exp(-(xs*xs-ya*ya) -2*i*xs*ya)
            return 2.0 * Complex64::new((ya - xs) * (xs + ya), 2.0 * xs * y).cexp() - ret;
        } else {
            return ret;
        }
    } else {
        // Note: The test that seems to be suggested in the paper is x < sqrt(-log(DBL_MIN)), about 26.6, since otherwise exp(-x^2)
        // underflows to zero and sum1,sum2,sum4 are zero.  However, long before this occurs, the sum1,sum2,sum4 contributions are
        // negligible in double precision; I find that this happens for x > about 6, for all y.  On the other hand, I find that the case
        // where we compute all of the sums is faster (at least with the precomputed expa2n2 table) until about x=10.  Furthermore, if we
        // try to compute all of the sums for x > 20, I find that we sometimes run into numerical problems because underflow/overflow
        // problems start to appear in the various coefficients of the sums, below.  Therefore, we use x < 10 here.

        if x < 10.0 {
            let mut prod2ax = 1.0;
            let mut prodm2ax = 1.0;
            let expx2;

            if y.is_nan() {
                return Complex64::new(y, y);
            }

            // Somewhat ugly copy-and-paste duplication here, but I see significant speedups from using the special-case code
            // with the precomputed exponential, and the x < 5e-4 special case is needed for accuracy.

            if relerr == f64::EPSILON {
                // Use precomputed exp(-a2*(n*n)) table
                if x < 5e-4 {
                    // compute sum4 and sum5 together as sum5-sum4
                    let x2 = x * x;
                    expx2 = 1.0 - x2 * (1.0 - 0.5 * x2); // exp(-x*x) via Taylor

                    // compute exp(2*a*x) and exp(-2*a*x) via Taylor, to double precision
                    let ax2 = 1.036642960860171859744 * x; // 2*a*x
                    let exp2ax =
                        1.0 + ax2 * (1.0 + ax2 * (0.5 + 0.166666666666666666667 * ax2));
                    let expm2ax =
                        1.0 - ax2 * (1.0 - ax2 * (0.5 - 0.166666666666666666667 * ax2));
                    let mut n: i32 = 1;
                    loop {
                        let coef = EXPA2N2[n as usize - 1] * expx2
                            / (a2 * (f64::from(n) * f64::from(n)) + y * y);
                        prod2ax *= exp2ax;
                        prodm2ax *= expm2ax;
                        sum1 += coef;
                        sum2 += coef * prodm2ax;
                        sum3 += coef * prod2ax;

                        // really = sum5 - sum4
                        sum5 += coef
                            * (2.0 * a)
                            * f64::from(n)
                            * sinh_taylor((2.0 * a) * f64::from(n) * x);

                        // test convergence via sum3
                        if coef * prod2ax < relerr * sum3 {
                            break;
                        }
                        n += 1;
                    }
                } else {
                    // x > 5e-4, compute sum4 and sum5 separately
                    expx2 = (-x * x).exp();
                    let exp2ax = ((2.0 * a) * x).exp();
                    let expm2ax = 1.0 / exp2ax;
                    let mut n: i32 = 1;
                    loop {
                        let coef = EXPA2N2[n as usize - 1] * expx2
                            / (a2 * (f64::from(n) * f64::from(n)) + y * y);
                        prod2ax *= exp2ax;
                        prodm2ax *= expm2ax;
                        sum1 += coef;
                        sum2 += coef * prodm2ax;
                        sum4 += (coef * prodm2ax) * (a * n as f64);
                        sum3 += coef * prod2ax;
                        sum5 += (coef * prod2ax) * (a * n as f64);
                        // test convergence via sum5, since this sum has the slowest decay
                        if coef * prod2ax * a * f64::from(n) < relerr * sum5 {
                            break;
                        }
                        n += 1;
                    }
                }
            } else {
                // relerr != f64::EPSILON, compute exp(-a2*(n*n)) on the fly
                let exp2ax = (2.0 * a * x).exp();
                let expm2ax = 1.0 / exp2ax;
                if x < 5e-4 {
                    // compute sum4 and sum5 together as sum5-sum4
                    let x2 = x * x;
                    expx2 = 1.0 - x2 * (1.0 - 0.5 * x2); // exp(-x*x) via Taylor
                    let mut n: i32 = 1;
                    loop {
                        let coef = (-a2 * f64::from(n) * f64::from(n)).exp() * expx2
                            / (a2 * (f64::from(n) * f64::from(n)) + y * y);
                        prod2ax *= exp2ax;
                        prodm2ax *= expm2ax;
                        sum1 += coef;
                        sum2 += coef * prodm2ax;
                        sum3 += coef * prod2ax;

                        // Really = sum5 - sum4
                        sum5 += coef
                            * (2.0 * a)
                            * f64::from(n)
                            * sinh_taylor((2.0 * a) * f64::from(n) * x);

                        // Test convergence via sum3
                        if coef * prod2ax < relerr * sum3 {
                            break;
                        }
                        n += 1;
                    }
                } else {
                    // x > 5e-4, compute sum4 and sum5 separately
                    expx2 = (-x * x).exp();
                    let mut n: i32 = 1;
                    loop {
                        let coef = (-a2 * (f64::from(n) * f64::from(n))).exp() * expx2
                            / (a2 * (f64::from(n) * f64::from(n)) + y * y);
                        prod2ax *= exp2ax;
                        prodm2ax *= expm2ax;
                        sum1 += coef;
                        sum2 += coef * prodm2ax;
                        sum4 += (coef * prodm2ax) * (a * n as f64);
                        sum3 += coef * prod2ax;
                        sum5 += (coef * prod2ax) * (a * n as f64);
                        // Test convergence via sum5, since this sum has the slowest decay
                        if (coef * prod2ax) * (a * n as f64) < relerr * sum5 {
                            break;
                        }
                        n += 1;
                    }
                }
            }

            let expx2erfcxy: f64 = if y > -6.0 {
                expx2 * y.erfcx()
            } else {
                2.0 * (y * y - x * x).exp()
            }; // avoid spurious overflow for large negative y
            if y > 5.0 {
                // imaginary terms cancel
                let sinxy = (x * y).sin();
                ret = Complex64::new(
                    (expx2erfcxy - c * y * sum1) * (2.0 * x * y).cos()
                        + (c * x * expx2) * sinxy * sinc(x * y, sinxy),
                    0.0,
                );
            } else {
                let xs = z.re;
                let sinxy = (xs * y).sin();
                let sin2xy = (2.0 * xs * y).sin();
                let cos2xy = (2.0 * xs * y).cos();
                let coef1 = expx2erfcxy - c * y * sum1;
                let coef2 = c * xs * expx2;
                ret = Complex64::new(
                    coef1 * cos2xy + coef2 * sinxy * sinc(xs * y, sinxy),
                    coef2 * sinc(2.0 * xs * y, sin2xy) - coef1 * sin2xy,
                );
            }
        } else {
            // x >= 10: only sum3 & sum5 contribute (see above note)
            if x.is_nan() {
                return Complex64::new(x, x);
            }
            if y.is_nan() {
                return Complex64::new(y, y);
            }

            ret = Complex64::new((-x * x).exp(), 0.0); // |y| < 1e-10, so we only need exp(-x*x) term

            // Round instead of ceil as in original paper; note that x/a > 1 here
            let n0 = (x / a + 0.5).floor(); // sum in both directions, starting at n0
            let dx = a * n0 - x;
            sum3 = (-dx * dx).exp() / (a2 * (n0 * n0) + y * y);
            sum5 = a * n0 * sum3;
            let exp1 = (4.0 * a * dx).exp();
            let mut exp1dn = 1.0;
            let mut dn: i32 = 1;
            while dn < n0 as i32 {
                // loop over n0-dn and n0+dn terms
                let np = n0 + dn as f64;
                let nm = n0 - dn as f64;
                let mut tp = (-sqr(a * f64::from(dn) + dx)).exp();
                exp1dn *= exp1;
                let mut tm = tp * exp1dn; // trick to get tm from tp
                tp /= a2 * (np * np) + y * y;
                tm /= a2 * (nm * nm) + y * y;
                sum3 += tp + tm;
                sum5 += a * (np * tp + nm * tm);
                if a * (np * tp + nm * tm) < relerr * sum5 {
                    return ret + Complex64::new( (0.5 * c) * y * (sum2 + sum3), (0.5 * c) * (sum5 - sum4).copysign(z.re) );
                }
                dn += 1;
            }

            loop {
                // loop over n0+dn terms only (since n0-dn <= 0)
                let np = n0 + dn as f64;
                dn += 1;
                let tp = (-sqr(a * f64::from(dn) + dx)).exp() / (a2 * (np * np) + y * y);
                sum3 += tp;
                sum5 += a * np * tp;
                if a * np * tp < relerr * sum5 {
                    return ret + Complex64::new( (0.5 * c) * y * (sum2 + sum3), (0.5 * c) * (sum5 - sum4).copysign(z.re) );
                }
            }
        }
    }

    return ret + Complex64::new( (0.5 * c) * y * (sum2 + sum3), (0.5 * c) * (sum5 - sum4).copysign(z.re) );
}








// The complex scaled complementary error function

pub fn erfcx_with_relerror(z: Complex64, relerr: f64) -> Complex64 {
    w_with_relerror(Complex64::new(-z.im, z.re), relerr)
}







// The complex error function

pub fn erf_with_relerror(z: Complex64, relerr: f64) -> Complex64 {
    let x = z.re;
    let y = z.im;

    // If the imaginary part is 0, make sure the sign of 0 is preserved
    if y == 0.0 {
        return Complex64::new(x.erf(), y);
    }

    // If the real part is 0, make sure the sign of 0 is preserved.
    // Handle the case of y -> Inf limit manually, since exp(y^2) -> Inf
    // but Im[w(y)] -> 0, so IEEE will give us a NaN when it should be Inf.
    if x == 0.0 {
        let imag = if y * y > 720.0 {
            if y > 0.0 {
                f64::INFINITY
            } else {
                f64::NEG_INFINITY
            }
        } else {
            (y * y).exp() * y.w_im()
        };

        return Complex64::new(x, imag);
    }

    let m_re_z2: f64 = (y - x) * (x + y);                    // Re(-z^2), being careful of overflow
    let m_im_z2: f64 = -2.0 * x * y;                         // Im(-z^2)
    if m_re_z2 < -750.0 {
        // underflow
        if x >= 0.0 {
            return Complex64::new(1.0, 0.0);
        } else {
            return Complex64::new(-1.0, 0.0);
        }
    }

    // Handle positive and negative x via different formulas, using the mirror symmetries of w,
    // to avoid overflow/underflow problems from multiplying exponentially large and small quantities.
    if x >= 0.0 {
        if x < 8.0e-2 {
            if y.abs() < 1.0e-2 {
                // Use Taylor series for small |z|, to avoid cancellation inaccuracy
                //   erf(z) = 2/sqrt(pi) * z * (1 - z^2/3 + z^4/10 - z^6/42 + z^8/216 + ...)

                let mz2 = Complex64::new(m_re_z2, m_im_z2); // -z^2
                return z * (1.1283791670955125739 + mz2 * (0.37612638903183752464 + mz2 * (0.11283791670955125739 + mz2 * (0.026866170645131251760 + mz2 * 0.0052239776254421878422))));
            } else if m_im_z2.abs() < 5.0e-3 && x < 5.0e-3 {
                // For small |x| and small |xy|, use Taylor series to avoid cancellation inaccuracy:
                //     erf(x+iy) = erf(iy)
                //        + 2*exp(y^2)/sqrt(pi) *
                //          [ x * (1 - x^2 * (1+2y^2)/3 + x^4 * (3+12y^2+4y^4)/30 + ...
                //            - i * x^2 * y * (1 - x^2 * (3+2y^2)/6 + ...) ]
                //   where:
                //      erf(iy) = exp(y^2) * Im[w(y)]

                let x2: f64 = x * x;
                let y2: f64 = y * y;
                let expy2: f64 = y2.exp();
                return Complex64::new(
                    expy2 * x * (1.1283791670955125739
                            - x2 * (0.37612638903183752464 + 0.75225277806367504925 * y2)
                            + x2 * x2 * (0.11283791670955125739 + y2 * (0.45135166683820502956 + 0.15045055561273500986 * y2))),
                    expy2 * (y.w_im() - x2 * y * (1.1283791670955125739 - x2 * (0.56418958354775628695 + 0.37612638903183752464 * y2)))
                );
            }
        }

        // Don't use the complex exp function, since that will produce spurious NaN values when multiplying w in an overflow situation.
        return 1.0 - m_re_z2.exp() * Complex64::new(m_im_z2.cos(), m_im_z2.sin()) * w_with_relerror(Complex64::new(-y, x), relerr);
    } else {
        // x < 0
        if x > -8.0e-2 {
            // duplicate from above to avoid abs(x) call
            if y.abs() < 1.0e-2 {
                // Use Taylor series for small |z|, to avoid cancellation inaccuracy
                //   erf(z) = 2/sqrt(pi) * z * (1 - z^2/3 + z^4/10 - z^6/42 + z^8/216 + ...)

                let mz2: Complex64 = Complex64::new(m_re_z2, m_im_z2); // -z^2
                return z * (1.1283791670955125739 + mz2 * (0.37612638903183752464 + mz2 * (0.11283791670955125739 + mz2 * (0.026866170645131251760 + mz2 * 0.0052239776254421878422))));
            } else if m_im_z2.abs() < 5.0e-3 && x > -5.0e-3 {
                // For small |x| and small |xy|, use Taylor series to avoid cancellation inaccuracy:
                //     erf(x+iy) = erf(iy)
                //        + 2*exp(y^2)/sqrt(pi) *
                //          [ x * (1 - x^2 * (1+2y^2)/3 + x^4 * (3+12y^2+4y^4)/30 + ...
                //            - i * x^2 * y * (1 - x^2 * (3+2y^2)/6 + ...) ]
                //   where:
                //      erf(iy) = exp(y^2) * Im[w(y)]

                let x2: f64 = x * x;
                let y2: f64 = y * y;
                let expy2: f64 = y2.exp();
                return Complex64::new(
                    expy2 * x * (1.1283791670955125739
                            - x2 * (0.37612638903183752464 + 0.75225277806367504925 * y2)
                            + x2 * x2 * (0.11283791670955125739 + y2 * (0.45135166683820502956 + 0.15045055561273500986 * y2))),
                    expy2 * (y.w_im() - x2 * y * (1.1283791670955125739 - x2 * (0.56418958354775628695 + 0.37612638903183752464 * y2)))
                );
            }
        } else if x.is_nan() {
            if y == 0.0 {
                return Complex64::new(f64::NAN, 0.0);
            } else {
                return Complex64::new(f64::NAN, f64::NAN);
            }
        }

        // Don't use the complex exp function, since that will produce spurious NaN values when multiplying w in an overflow situation.
        // Although in principle complex multiplication should be associative, in numerical applications this is not necessarily the case.
        // In the following I therefore deliberately write a * (b * c) rather than a * b * c. If not, one of the test cases will break
        // and return (NaN, -Inf) rather than (-Inf, -Inf).

        m_re_z2.exp() * (Complex64::new(m_im_z2.cos(), m_im_z2.sin()) * w_with_relerror(Complex64::new(y, -x), relerr)) - 1.0
    }
}






// The imaginary error function

pub fn erfi_with_relerror(z: Complex64, relerr: f64) -> Complex64 {
    let e: Complex64 = erf_with_relerror(Complex64::new(-z.im, z.re), relerr);
    Complex64::new(e.im, -e.re)
}





// The complex complementary error function

pub fn erfc_with_relerror(z: Complex64, relerr: f64) -> Complex64 {
    let x: f64 = z.re;
    let y: f64 = z.im;

    if x == 0.0 {
        // Handle y -> Inf limit manualyy, since exp(y^2) -> Inf but Im[w(y)] -> 0, so
        // IEEE will give us a NaN when it should be Inf

        if y * y > 720.0 {
            if y > 0.0 {
                return Complex64::new(1.0, f64::NEG_INFINITY);
            } else {
                return Complex64::new(1.0, f64::INFINITY);
            }
        } else {
            return Complex64::new(1.0, -(y * y).exp() * y.w_im());
        }
    }

    if y == 0.0 {
        if x * x > 750.0 {
            // underflow
            if x >= 0.0 {
                return Complex64::new(0.0, -y); // Preserve sign of 0
            } else {
                return Complex64::new(2.0, -y); // Preserve sign of 0
            }
        } else {
            if x >= 0.0 {
                return Complex64::new((-x*x).exp() * x.erfcx(), -y);
            } else {
                return Complex64::new(2.0 - (-x*x).exp() * (-x).erfcx(), -y);
            }
        }
    }

    let m_re_z2: f64 = (y - x) * (x + y);      // Re(-z^2), being careful of overflow
    let m_im_z2: f64 = -2.0 * x * y;           // Im(-z^2)
    if m_re_z2 < -750.0 {
        // underflow
        if x >= 0.0 {
            return Complex64::new(0.0, 0.0);
        } else {
            return Complex64::new(2.0, 0.0);
        }
    }

    if x >= 0.0 {
        return Complex64::new(m_re_z2, m_im_z2).cexp() * w_with_relerror(Complex64::new(-y, x), relerr);
    } else {
        return 2.0 - Complex64::new(m_re_z2, m_im_z2).cexp() * w_with_relerror(Complex64::new(y, -x), relerr);
    }
}





// Dawson's function 

pub fn dawson_with_relerror(z: Complex64, relerr: f64) -> Complex64 {
    const SPI2: f64 = 0.8862269254527580136490837416705725913990; // sqrt(pi)/2
    let x = z.re;
    let y = z.im;

    // Handle axes separately for speed & proper handling of x or y = Inf or NaN

    if y == 0.0 {
        return Complex64::new(SPI2 * x.w_im(), -y); // preserve sign of 0
    }

    if x == 0.0 {
        let y2 = y * y;
        if y2 < 2.5e-5 {
            // Taylor expansion
            return Complex64::new( x, y * (1. + y2 * (0.6666666666666666666666666666666666666667 + y2 * 0.26666666666666666666666666666666666667)));
        } else {
            // Preserve sign of 0
            if y >= 0.0 {
                return Complex64::new(x, SPI2 * (y2.exp() - y.erfcx()));
            } else {
                return Complex64::new(x, SPI2 * ((-y).erfcx() - y2.exp()));
            }
        }
    }

    let m_re_z2 = (y - x) * (x + y);             // Re(-z^2), being careful of overflow
    let m_im_z2 = -2.0 * x * y;                  // Im(-z^2)
    let mz2 = Complex64::new(m_re_z2, m_im_z2);  // -z^2

    // Handle positive and negative x via different formulas, using the mirror symmetries of w, to avoid overflow/underflow
    // problems from multiplying exponentially large and small quantities.

    if y >= 0.0 {
        if y < 5e-3 {
            if x.abs() < 5e-3 {
                // Use Taylor series for small |z|, to avoid cancellation inaccuracy: dawson(z) = z - 2/3 z^3 + 4/15 z^5 + ...
                return z * (1. + mz2 * (0.6666666666666666666666666666666666666667 + mz2 * 0.2666666666666666666666666666666666666667));
            } else if m_im_z2.abs() < 5.0e-3 {
                // (**) For small |y| and small |xy|, use Taylor series to avoid cancellation inaccuracy:
                //     dawson(x + iy) = D + y^2 (D + x - 2Dx^2) + y^4 (D/2 + 5x/6 - 2Dx^2 - x^3/3 + 2Dx^4/3)
                //                        + iy [ (1-2Dx) + 2/3 y^2 (1 - 3Dx - x^2 + 2Dx^3)
                //                        + y^4/15 (4 - 15Dx - 9x^2 + 20Dx^3 + 2x^4 - 4Dx^5) ] + ...
                // where D = dawson(x)
                //
                // However, for large |x|, 2Dx -> 1 which gives cancellation problems in this series (many of the leading terms cancel).
                // So, for large |x|, we need to substitute a continued-fraction expansion for D:
                //     dawson(x) = 0.5 / (x-0.5/(x-1/(x-1.5/(x-2/(x-2.5/(x...))))))
                //
                // The 6 terms shown here seems to be the minimum needed to be accurate as soon as the simpler Taylor expansion above starts
                // breaking down.  Using this 6-term expansion, factoring out the denominator, and simplifying with Maple, we obtain:
                //     Re dawson(x + iy) * (-15 + 90x^2 - 60x^4 + 8x^6) / x = 33 - 28x^2 + 4x^4 + y^2 (18 - 4x^2) + 4 y^4
                //     Im dawson(x + iy) * (-15 + 90x^2 - 60x^4 + 8x^6) / y = -15 + 24x^2 - 4x^4 + 2/3 y^2 (6x^2 - 15) - 4 y^4
                //
                // Finally, for |x| > 5e7, we can use a simpler 1-term continued-fraction expansion for the real part,
                // and a 2-term expansion for the imaginary part. This avoids overflow problems for huge |x|.  This yields:
                //     Re dawson(x + iy) = [1 + y^2 (1 + y^2/2 - (xy)^2/3)] / (2x)
                //     Im dawson(x + iy) = y [ -1 - 2/3 y^2 + y^4/15 (2x^2 - 4) ] / (2x^2 - 1)

                let x2 = x * x;
                if x2 > 1600.0 {
                    // |x| > 40
                    let y2 = y * y;
                    if x2 > 25.0e14 {
                        // |x| > 5e7
                        let xy2 = (x * y) * (x * y);
                        return Complex64::new(
                            (0.5 + y2 * (0.5 + 0.25 * y2 - 0.16666666666666666667 * xy2)) / x,
                            y * (-1.0
                                + y2 * (-0.66666666666666666667
                                    + 0.13333333333333333333 * xy2
                                    - 0.26666666666666666667 * y2))
                                / (2.0 * x2 - 1.0),
                        );
                    }
                    return (1. / (-15. + x2 * (90. + x2 * (-60. + 8. * x2))))
                        * Complex64::new(
                            x * (33. + x2 * (-28. + 4. * x2) + y2 * (18. - 4. * x2 + 4. * y2)),
                            y * (-15. + x2 * (24. - 4. * x2) + y2 * (4. * x2 - 10. - 4. * y2)),
                        );
                } else {
                    let d = SPI2 * x.w_im();
                    let y2 = y * y;
                    return Complex64::new(
                        d + y2 * (d + x - 2. * d * x2)
                            + y2 * y2
                                * (d * (0.5 - x2 * (2. - 0.66666666666666666667 * x2))
                                    + x * (0.83333333333333333333
                                        - 0.33333333333333333333 * x2)),
                        y * (1. - 2. * d * x
                            + y2 * 0.66666666666666666667 * (1. - x2 - d * x * (3. - 2. * x2))
                            + y2 * y2
                                * (0.26666666666666666667
                                    - x2 * (0.6 - 0.13333333333333333333 * x2)
                                    - d * x
                                        * (1.
                                            - x2 * (1.3333333333333333333
                                                - 0.26666666666666666667 * x2)))),
                    );
                }
            }
        }
        let res: Complex64 = mz2.cexp() - w_with_relerror(z, relerr);
        return SPI2 * Complex64::new(-res.im, res.re);
    } else {
        // y < 0
        if y > -5.0e-3 {
            // duplicate from above to avoid abs(x) call
            if x.abs() < 5.0e-3 {
                // Use Taylor series for small |z|, to avoid cancellation inaccuracy: dawson(z) = z - 2/3 z^3 + 4/15 z^5 + ...
                return z * (1. + mz2 * (0.6666666666666666666666666666666666666667 + mz2 * 0.2666666666666666666666666666666666666667));
            } else if m_im_z2.abs() < 5.0e-3 {
                // For small |y| and small |xy|, use Taylor series to avoid cancellation inaccuracy.
                // See explanation above at (**).
                let x2 = x * x;
                if x2 > 1600.0 {
                    // |x| > 40
                    let y2 = y * y;
                    if x2 > 25.0e14 {
                        // |x| > 5e7
                        let xy2 = (x * y) * (x * y);
                        return Complex64::new(
                            (0.5 + y2 * (0.5 + 0.25 * y2 - 0.16666666666666666667 * xy2)) / x,
                            y * (-1.0
                                + y2 * (-0.66666666666666666667
                                    + 0.13333333333333333333 * xy2
                                    - 0.26666666666666666667 * y2))
                                / (2.0 * x2 - 1.0),
                        );
                    }
                    return (1. / (-15. + x2 * (90. + x2 * (-60. + 8. * x2))))
                        * Complex64::new(
                            x * (33. + x2 * (-28. + 4. * x2) + y2 * (18. - 4. * x2 + 4. * y2)),
                            y * (-15. + x2 * (24. - 4. * x2) + y2 * (4. * x2 - 10. - 4. * y2))
                        );
                } else {
                    let d = SPI2 * x.w_im();
                    let y2 = y * y;
                    return Complex64::new(
                        d + y2 * (d + x - 2. * d * x2)
                            + y2 * y2
                                * (d * (0.5 - x2 * (2. - 0.66666666666666666667 * x2))
                                    + x * (0.83333333333333333333
                                        - 0.33333333333333333333 * x2)),
                        y * (1. - 2. * d * x
                            + y2 * 0.66666666666666666667 * (1. - x2 - d * x * (3. - 2. * x2))
                            + y2 * y2
                                * (0.26666666666666666667
                                    - x2 * (0.6 - 0.13333333333333333333 * x2)
                                    - d * x
                                        * (1.
                                            - x2 * (1.3333333333333333333
                                                - 0.26666666666666666667 * x2)))),
                    );
                }
            }
        } else if y.is_nan() {
            if x == 0.0 {
                return Complex64::new(0.0, f64::NAN);
            } else {
                return Complex64::new(f64::NAN, f64::NAN);
            }
        }
        let res: Complex64 = w_with_relerror(-z, relerr) - mz2.cexp();
        return SPI2 * Complex64::new(-res.im, res.re);
    }
}





pub trait ComplexErrorFunctions {
    fn erfcx(self) -> Self;
    fn erf(self) -> Self;
    fn w(self) -> Self;
    fn erfi(self) -> Self;      
    fn erfc(self) -> Self;
    fn dawson(self) -> Self;
}



impl ComplexErrorFunctions for Complex64 {

    fn erfcx(self) -> Self {
        erfcx_with_relerror(self, 0.0)
    }

    fn erf(self) -> Self {
        erf_with_relerror(self, 0.0)
    }

    fn w(self) -> Self {
        w_with_relerror(self, 0.0)
    }

    fn erfi(self) -> Self {
        erfi_with_relerror(self, 0.0)
    }

    fn erfc(self) -> Self {
        erfc_with_relerror(self, 0.0)
    }

    fn dawson(self) -> Self {
        dawson_with_relerror(self, 0.0)
    }
}









#[cfg(test)]
mod tests {

    use super::*;
    use crate::auxilliary::relerr;

    #[test]
    fn test_w_complex() {
        // Function arguments z.

        let z: [Complex64; 57] = [
            Complex64::new(624.2, -0.26123),
            Complex64::new(-0.4, 3.),
            Complex64::new(0.6, 2.),
            Complex64::new(-1., 1.),
            Complex64::new(-1., -9.),
            Complex64::new(-1., 9.),
            Complex64::new(-0.0000000234545, 1.1234),
            Complex64::new(-3., 5.1),
            Complex64::new(-53.0, 30.1),
            Complex64::new(0.0, 0.12345),
            Complex64::new(11.0, 1.0),
            Complex64::new(-22.0, -2.0),
            Complex64::new(9.0, -28.0),
            Complex64::new(21.0, -33.0),
            Complex64::new(1e5, 1e5),
            Complex64::new(1e14, 1e14),
            Complex64::new(-3001.0, -1000.0),
            Complex64::new(1e160, -1e159),
            Complex64::new(-6.01, 0.01),
            Complex64::new(-0.7, -0.7),
            Complex64::new(2.611780000000000e+01, 4.540909610972489e+03),
            Complex64::new(0.8e7, 0.3e7),
            Complex64::new(-20.0, -19.8081),
            Complex64::new(1e-16, -1.1e-16),
            Complex64::new(2.3e-8, 1.3e-8),
            Complex64::new(6.3, -1e-13),
            Complex64::new(6.3, 1e-20),
            Complex64::new(1e-20, 6.3),
            Complex64::new(1e-20, 16.3),
            Complex64::new(9.0, 1e-300),
            Complex64::new(6.01, 0.11),
            Complex64::new(8.01, 1.01e-10),
            Complex64::new(28.01, 1e-300),
            Complex64::new(10.01, 1e-200),
            Complex64::new(10.01, -1e-200),
            Complex64::new(10.01, 0.99e-10),
            Complex64::new(10.01, -0.99e-10),
            Complex64::new(1e-20, 7.01),
            Complex64::new(-1.0, 7.01),
            Complex64::new(5.99, 7.01),
            Complex64::new(1.0, 0.0),
            Complex64::new(55.0, 0.0),
            Complex64::new(-0.1, 0.0),
            Complex64::new(1e-20, 0.0),
            Complex64::new(0.0, 5e-14),
            Complex64::new(0.0, 51.0),
            Complex64::new(f64::INFINITY, 0.0),
            Complex64::new(f64::NEG_INFINITY, 0.0),
            Complex64::new(0.0, f64::INFINITY),
            Complex64::new(0.0, f64::NEG_INFINITY),
            Complex64::new(f64::INFINITY, f64::INFINITY),
            Complex64::new(f64::INFINITY, f64::NEG_INFINITY),
            Complex64::new(f64::NAN, f64::NAN),
            Complex64::new(f64::NAN, 0.0),
            Complex64::new(0.0, f64::NAN),
            Complex64::new(f64::NAN, f64::INFINITY),
            Complex64::new(f64::INFINITY, f64::NAN),
        ];

        // w(z) computed with WolframAlpha. WolframAlpha is problematic for some of the above inputs, so I had to use the
        // continued-fraction expansion in WolframAlpha in some cases, or switch to Maple

        let expected_w: [Complex64; 57] = [
            Complex64::new( -3.78270245518980507452677445620103199303131110e-7, 0.000903861276433172057331093754199933411710053155),
            Complex64::new( 0.1764906227004816847297495349730234591778719532788, -0.02146550539468457616788719893991501311573031095617),
            Complex64::new( 0.2410250715772692146133539023007113781272362309451, 0.06087579663428089745895459735240964093522265589350),
            Complex64::new( 0.30474420525691259245713884106959496013413834051768, -0.20821893820283162728743734725471561394145872072738),
            Complex64::new( 7.317131068972378096865595229600561710140617977e34, 8.321873499714402777186848353320412813066170427e34),
            Complex64::new( 0.0615698507236323685519612934241429530190806818395, -0.00676005783716575013073036218018565206070072304635),
            Complex64::new( 0.3960793007699874918961319170187598400134746631, -5.593152259116644920546186222529802777409274656e-9),
            Complex64::new( 0.08217199226739447943295069917990417630675021771804, -0.04701291087643609891018366143118110965272615832184),
            Complex64::new( 0.00457246000350281640952328010227885008541748668738, -0.00804900791411691821818731763401840373998654987934),
            Complex64::new(0.8746342859608052666092782112565360755791467973338452, 0.),
            Complex64::new( 0.00468190164965444174367477874864366058339647648741, 0.0510735563901306197993676329845149741675029197050),
            Complex64::new( -0.0023193175200187620902125853834909543869428763219, -0.025460054739731556004902057663500272721780776336),
            Complex64::new( 9.11463368405637174660562096516414499772662584e304, 3.97101807145263333769664875189354358563218932e305),
            Complex64::new( -4.4927207857715598976165541011143706155432296e281, -2.8019591213423077494444700357168707775769028e281),
            Complex64::new( 2.820947917809305132678577516325951485807107151e-6, 2.820947917668257736791638444590253942253354058e-6),
            Complex64::new( 2.82094791773878143474039725787438662716372268e-15, 2.82094791773878143474039725773333923127678361e-15),
            Complex64::new( -0.0000563851289696244350147899376081488003110150498, -0.000169211755126812174631861529808288295454992688),
            Complex64::new( -5.586035480670854326218608431294778077663867e-162, 5.586035480670854326218608431294778077663867e-161),
            Complex64::new( 0.00016318325137140451888255634399123461580248456, -0.095232456573009287370728788146686162555021209999),
            Complex64::new( 0.69504753678406939989115375989939096800793577783885, -1.8916411171103639136680830887017670616339912024317),
            Complex64::new( 0.0001242418269653279656612334210746733213167234822, 7.145975826320186888508563111992099992116786763e-7),
            Complex64::new( 2.318587329648353318615800865959225429377529825e-8, 6.182899545728857485721417893323317843200933380e-8),
            Complex64::new( -0.0133426877243506022053521927604277115767311800303, -0.0148087097143220769493341484176979826888871576145),
            Complex64::new( 1.00000000000000012412170838050638522857747934, 1.12837916709551279389615890312156495593616433e-16),
            Complex64::new( 0.9999999853310704677583504063775310832036830015, 2.595272024519678881897196435157270184030360773e-8),
            Complex64::new( -1.4731421795638279504242963027196663601154624e-15, 0.090727659684127365236479098488823462473074709),
            Complex64::new( 5.79246077884410284575834156425396800754409308e-18, 0.0907276596841273652364790985059772809093822374),
            Complex64::new( 0.0884658993528521953466533278764830881245144368, 1.37088352495749125283269718778582613192166760e-22),
            Complex64::new( 0.0345480845419190424370085249304184266813447878, 2.11161102895179044968099038990446187626075258e-23),
            Complex64::new( 6.63967719958073440070225527042829242391918213e-36, 0.0630820900592582863713653132559743161572639353),
            Complex64::new( 0.00179435233208702644891092397579091030658500743634, 0.0951983814805270647939647438459699953990788064762),
            Complex64::new( 9.09760377102097999924241322094863528771095448e-13, 0.0709979210725138550986782242355007611074966717),
            Complex64::new( 7.2049510279742166460047102593255688682910274423e-304, 0.0201552956479526953866611812593266285000876784321),
            Complex64::new( 3.04543604652250734193622967873276113872279682e-44, 0.0566481651760675042930042117726713294607499165),
            Complex64::new( 3.04543604652250734193622967873276113872279682e-44, 0.0566481651760675042930042117726713294607499165),
            Complex64::new( 0.5659928732065273429286988428080855057102069081e-12, 0.056648165176067504292998527162143030538756683302),
            Complex64::new( -0.56599287320652734292869884280802459698927645e-12, 0.0566481651760675042929985271621430305387566833029),
            Complex64::new( 0.0796884251721652215687859778119964009569455462, 1.11474461817561675017794941973556302717225126e-22),
            Complex64::new( 0.07817195821247357458545539935996687005781943386550, -0.01093913670103576690766705513142246633056714279654),
            Complex64::new( 0.04670032980990449912809326141164730850466208439937, 0.03944038961933534137558064191650437353429669886545),
            Complex64::new( 0.36787944117144232159552377016146086744581113103176, 0.60715770584139372911503823580074492116122092866515),
            Complex64::new(0.0, 0.010259688805536830986089913987516716056946786526145), 
            Complex64::new( 0.99004983374916805357390597718003655777207908125383, -0.11208866436449538036721343053869621153527769495574),
            Complex64::new( 0.99999999999999999999999999999999999999990000, 1.12837916709551257389615890312154517168802603e-20),
            Complex64::new(0.999999999999943581041645226871305192054749891144158, 0.0),
            Complex64::new(0.0110604154853277201542582159216317923453996211744250, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(f64::INFINITY, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(f64::NAN, f64::NAN),
            Complex64::new(f64::NAN, f64::NAN),
            Complex64::new(f64::NAN, f64::NAN),
            Complex64::new(f64::NAN, 0.0),
            Complex64::new(f64::NAN, f64::NAN),
            Complex64::new(f64::NAN, f64::NAN),
        ];

        let tolerance = 1.0e-13;
        for n in 0..z.len() {
            let computed_w = z[n].w();
            // println!( "{}: z={:?} --> w(z) computed: {:?}  vs expected:  {:?}", n, z[n], computed_w, expected_w[n]);
            let relative_error_re = relerr(computed_w.re, expected_w[n].re);
            let relative_error_im = relerr(computed_w.im, expected_w[n].im);
            assert!(relative_error_re < tolerance);
            assert!(relative_error_im < tolerance);
        }
    }



    #[test]
    fn test_erf_complex() {

        let z: [Complex64; 41] = [
            Complex64::new(1.0, 2.0),
            Complex64::new(-1.0, 2.0),
            Complex64::new(1.0, -2.0),
            Complex64::new(-1.0, -2.0),
            Complex64::new(9.0, -28.0),
            Complex64::new(21.0, -33.0),
            Complex64::new(1.0e3, 1.0e3),
            Complex64::new(-3001.0, -1000.0),
            Complex64::new(1.0e160, -1.0e159),
            Complex64::new(5.1e-3, 1.0e-8),
            Complex64::new(-4.9e-3, 4.95e-3),
            Complex64::new(4.9e-3, 0.5),
            Complex64::new(4.9e-4, -0.5e1),
            Complex64::new(-4.9e-5, -0.5e2),
            Complex64::new(5.1e-3, 0.5),
            Complex64::new(5.1e-4, -0.5e1),
            Complex64::new(-5.1e-5, -0.5e2),
            Complex64::new(1e-6, 2e-6),
            Complex64::new(0.0, 2e-6),
            Complex64::new(0.0, 2.0),
            Complex64::new(0.0, 20.0),
            Complex64::new(0.0, 200.0),
            Complex64::new(f64::INFINITY, 0.0),
            Complex64::new(f64::NEG_INFINITY, 0.0),
            Complex64::new(0.0, f64::INFINITY),
            Complex64::new(0.0, f64::NEG_INFINITY),
            Complex64::new(f64::INFINITY, f64::INFINITY),
            Complex64::new(f64::INFINITY, f64::NEG_INFINITY),
            Complex64::new(f64::NAN, f64::NAN),
            Complex64::new(f64::NAN, 0.0),
            Complex64::new(0.0, f64::NAN),
            Complex64::new(f64::NAN, f64::INFINITY),
            Complex64::new(f64::INFINITY, f64::NAN),
            Complex64::new(1e-3, f64::NAN),
            Complex64::new(7e-2, 7e-2),
            Complex64::new(7e-2, -7e-4),
            Complex64::new(-9e-2, 7e-4),
            Complex64::new(-9e-2, 9e-2),
            Complex64::new(-7e-4, 9e-2),
            Complex64::new(7e-2, 0.9e-2),
            Complex64::new(7e-2, 1.1e-2)
        ];
        
        let expected_erf: [Complex64; 41] = [
            Complex64::new(-0.5366435657785650339917955593141927494421, -5.049143703447034669543036958614140565553),
            Complex64::new(0.5366435657785650339917955593141927494421, -5.049143703447034669543036958614140565553),
            Complex64::new(-0.5366435657785650339917955593141927494421, 5.049143703447034669543036958614140565553),
            Complex64::new(0.5366435657785650339917955593141927494421, 5.049143703447034669543036958614140565553),
            Complex64::new(0.3359473673830576996788000505817956637777e304, -0.1999896139679880888755589794455069208455e304),
            Complex64::new(0.3584459971462946066523939204836760283645e278, 0.3818954885257184373734213077678011282505e280),
            Complex64::new(0.9996020422657148639102150147542224526887, 0.00002801044116908227889681753993542916894856),
            Complex64::new(-1.0, 0.0),
            Complex64::new(1.0, 0.0),
            Complex64::new(0.005754683859034800134412990541076554934877, 0.1128349818335058741511924929801267822634e-7),
            Complex64::new(-0.005529149142341821193633460286828381876955, 0.005585388387864706679609092447916333443570),
            Complex64::new(0.007099365669981359632319829148438283865814, 0.6149347012854211635026981277569074001219),
            Complex64::new(0.3981176338702323417718189922039863062440e8, -0.8298176341665249121085423917575122140650e10),
            Complex64::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
            Complex64::new(0.007389128308257135427153919483147229573895, 0.6149332524601658796226417164791221815139),
            Complex64::new(0.4143671923267934479245651547534414976991e8, -0.8298168216818314211557046346850921446950e10),
            Complex64::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
            Complex64::new(0.1128379167099649964175513742247082845155e-5, 0.2256758334191777400570377193451519478895e-5),
            Complex64::new(0.0, 0.2256758334194034158904576117253481476197e-5),
            Complex64::new(0.0, 18.56480241457555259870429191324101719886),
            Complex64::new(0.0, 0.1474797539628786202447733153131835124599e173),
            Complex64::new(0.0, f64::INFINITY),
            Complex64::new(1.0,0.0),
            Complex64::new(-1.0,0.0),
            Complex64::new(0.0,f64::INFINITY),
            Complex64::new(0.0,f64::NEG_INFINITY),
            Complex64::new(f64::NAN,f64::NAN),
            Complex64::new(f64::NAN,f64::NAN),
            Complex64::new(f64::NAN,f64::NAN),
            Complex64::new(f64::NAN,0.0),
            Complex64::new(0.0,f64::NAN),
            Complex64::new(f64::NAN,f64::NAN),
            Complex64::new(f64::NAN,f64::NAN),
            Complex64::new(f64::NAN,f64::NAN),
            Complex64::new(0.07924380404615782687930591956705225541145, 0.07872776218046681145537914954027729115247),
            Complex64::new(0.07885775828512276968931773651224684454495, -0.0007860046704118224342390725280161272277506),
            Complex64::new(-0.1012806432747198859687963080684978759881, 0.0007834934747022035607566216654982820299469),
            Complex64::new(-0.1020998418798097910247132140051062512527, 0.1010030778892310851309082083238896270340),
            Complex64::new(-0.0007962891763147907785684591823889484764272, 0.1018289385936278171741809237435404896152),
            Complex64::new(0.07886408666470478681566329888615410479530, 0.01010604288780868961492224347707949372245),
            Complex64::new(0.07886723099940260286824654364807981336591, 0.01235199327873258197931147306290916629654)
        ];

        let tolerance = 1.0e-13;
        for n in 0..z.len() {
            let computed_erf = z[n].erf();
            let relative_error_re = relerr(computed_erf.re, expected_erf[n].re);
            let relative_error_im = relerr(computed_erf.im, expected_erf[n].im);
            assert!(relative_error_re < tolerance);
            assert!(relative_error_im < tolerance);
        }
    }



    #[test]
    fn test_erfi_complex() {
        let tolerance = 1.0e-13;
        let z = Complex64::new(1.234, 0.5678);
        let computed_erfi = erfi_with_relerror(z, 0.0);
        let expected_erfi = Complex64::new(1.081032284405373149432716643834106923212, 1.926775520840916645838949402886591180834);
        let relative_error_re = relerr(computed_erfi.re, expected_erfi.re);
        let relative_error_im = relerr(computed_erfi.im, expected_erfi.im);
        assert!(relative_error_re < tolerance);
        assert!(relative_error_im < tolerance);
    }



    #[test]
    fn test_erfcx_complex() {
        let tolerance = 1.0e-13;
        let z = Complex64::new(1.234, 0.5678);
        let computed_erfcx = erfcx_with_relerror(z, 0.0);
        let expected_erfcx = Complex64::new(0.3382187479799972294747793561190487832579, -0.1116077470811648467464927471872945833154);
        let relative_error_re = relerr(computed_erfcx.re, expected_erfcx.re);
        let relative_error_im = relerr(computed_erfcx.im, expected_erfcx.im);
        assert!(relative_error_re < tolerance);
        assert!(relative_error_im < tolerance);
    }




    #[test]
    fn test_erfc_complex() {

        let z: [Complex64; 30] = [
            Complex64::new(1.0, 2.0),
            Complex64::new(-1.0, 2.0),
            Complex64::new(1.0, -2.0),
            Complex64::new(-1.0, -2.0),
            Complex64::new(9.0, -28.0),
            Complex64::new(21.0, -33.0),
            Complex64::new(1.0e3, 1.0e3),
            Complex64::new(-3001.0, -1000.0),
            Complex64::new(1.0e160, -1.0e159),
            Complex64::new(5.1e-3, 1.0e-8),
            Complex64::new(0.0, 2.0e-6),
            Complex64::new(0.0, 2.0),
            Complex64::new(0.0, 20.0),
            Complex64::new(0.0, 200.0),
            Complex64::new(2.0e-6, 0.0),
            Complex64::new(2.0, 0.0),
            Complex64::new(20.0, 0.0),
            Complex64::new(200.0, 0.0),
            Complex64::new(f64::INFINITY, 0.0),
            Complex64::new(f64::NEG_INFINITY, 0.0),
            Complex64::new(0.0, f64::INFINITY),
            Complex64::new(0.0,f64::NEG_INFINITY),
            Complex64::new(f64::INFINITY, f64::INFINITY),
            Complex64::new(f64::INFINITY, f64::NEG_INFINITY),
            Complex64::new(f64::NAN, f64::NAN),
            Complex64::new(f64::NAN, 0.0),
            Complex64::new(0.0, f64::NAN),
            Complex64::new(f64::NAN, f64::INFINITY),
            Complex64::new(f64::INFINITY, f64::NAN),
            Complex64::new(88.0, 0.0)
        ];
        
        let expected_erfc: [Complex64; 30] = [
            Complex64::new(1.536643565778565033991795559314192749442, 5.049143703447034669543036958614140565553),
            Complex64::new(0.4633564342214349660082044406858072505579, 5.049143703447034669543036958614140565553),
            Complex64::new(1.536643565778565033991795559314192749442, -5.049143703447034669543036958614140565553),
            Complex64::new(0.4633564342214349660082044406858072505579, -5.049143703447034669543036958614140565553),
            Complex64::new(-0.3359473673830576996788000505817956637777e304, 0.1999896139679880888755589794455069208455e304),
            Complex64::new(-0.3584459971462946066523939204836760283645e278, -0.3818954885257184373734213077678011282505e280),
            Complex64::new(0.0003979577342851360897849852457775473112748, -0.00002801044116908227889681753993542916894856),
            Complex64::new(2.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.9942453161409651998655870094589234450651, -0.1128349818335058741511924929801267822634e-7),
            Complex64::new(1.0, -0.2256758334194034158904576117253481476197e-5),
            Complex64::new(1.0, -18.56480241457555259870429191324101719886),
            Complex64::new(1.0, -0.1474797539628786202447733153131835124599e173),
            Complex64::new(1.0, f64::NEG_INFINITY),
            Complex64::new(0.9999977432416658119838633199332831406314, 0.0),
            Complex64::new(0.004677734981047265837930743632747071389108, 0.0),
            Complex64::new(0.5395865611607900928934999167905345604088e-175, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(2.0, 0.0),
            Complex64::new(1.0, f64::NEG_INFINITY),
            Complex64::new(1.0, f64::INFINITY),
            Complex64::new(f64::NAN, f64::NAN),
            Complex64::new(f64::NAN, f64::NAN),
            Complex64::new(f64::NAN, f64::NAN),
            Complex64::new(f64::NAN, 0.0),
            Complex64::new(1.0, f64::NAN),
            Complex64::new(f64::NAN, f64::NAN),
            Complex64::new(f64::NAN, f64::NAN),
            Complex64::new(0.0, 0.0)
        ];

        let tolerance = 1.0e-13;
        for n in 0..z.len() {
            let computed_erfc = erfc_with_relerror(z[n], 0.);
            let relative_error_re = relerr(computed_erfc.re, expected_erfc[n].re);
            let relative_error_im = relerr(computed_erfc.im, expected_erfc[n].im);
            assert!(relative_error_re < tolerance);
            assert!(relative_error_im < tolerance);
        }
    }




    #[test]
    fn test_dawson_complex() {
        let z: [Complex64; 48] = [
            Complex64::new(2.0, 1.0),
            Complex64::new(-2.0, 1.0),
            Complex64::new(2.0, -1.0),
            Complex64::new(-2.0, -1.0),
            Complex64::new(-28.0, 9.0),
            Complex64::new(33.0, -21.0),
            Complex64::new(1.0e3, 1.0e3),
            Complex64::new(-1000.0, -3001.0),
            Complex64::new(1.0e-8, 5.1e-3),
            Complex64::new(4.95e-3, -4.9e-3),
            Complex64::new(5.1e-3, 5.1e-3),
            Complex64::new(0.5, 4.9e-3),
            Complex64::new(-0.5e1, 4.9e-4),
            Complex64::new(-0.5e2, -4.9e-5),
            Complex64::new(0.5e3, 4.9e-6),
            Complex64::new(0.5, 5.1e-3),
            Complex64::new(-0.5e1, 5.1e-4),
            Complex64::new(-0.5e2, -5.1e-5),
            Complex64::new(1e-6, 2e-6),
            Complex64::new(2e-6, 0.0),
            Complex64::new(2.0, 0.0),
            Complex64::new(20.0, 0.0),
            Complex64::new(200.0, 0.0),
            Complex64::new(0.0, 4.9e-3),
            Complex64::new(0.0, -5.1e-3),
            Complex64::new(0.0, 2e-6),
            Complex64::new(0.0, -2.0),
            Complex64::new(0.0, 20.0),
            Complex64::new(0.0, -200.0),
            Complex64::new(f64::INFINITY, 0.0),
            Complex64::new(f64::NEG_INFINITY, 0.0),
            Complex64::new(0.0, f64::INFINITY),
            Complex64::new(0.0, f64::NEG_INFINITY),
            Complex64::new(f64::INFINITY, f64::INFINITY),
            Complex64::new(f64::INFINITY,f64::NEG_INFINITY),
            Complex64::new(f64::NAN,f64::NAN),
            Complex64::new(f64::NAN, 0.0),
            Complex64::new(0.0, f64::NAN),
            Complex64::new(f64::NAN, f64::INFINITY),
            Complex64::new(f64::INFINITY, f64::NAN),
            Complex64::new(39.0, 6.4e-5),
            Complex64::new(41.0, 6.09e-5),
            Complex64::new(4.9e7, 5.0e-11),
            Complex64::new(5.1e7, 4.8e-11),
            Complex64::new(1.0e9, 2.4e-12),
            Complex64::new(1.0e11, 2.4e-14),
            Complex64::new(1.0e13, 2.4e-16),
            Complex64::new(1.0e300, 2.4e-303)
        ];

        let expected_dawson: [Complex64; 48] = [
            Complex64::new(0.1635394094345355614904345232875688576839, -0.1531245755371229803585918112683241066853),
            Complex64::new(-0.1635394094345355614904345232875688576839, -0.1531245755371229803585918112683241066853),
            Complex64::new(0.1635394094345355614904345232875688576839, 0.1531245755371229803585918112683241066853),
            Complex64::new(-0.1635394094345355614904345232875688576839, 0.1531245755371229803585918112683241066853),
            Complex64::new(-0.01619082256681596362895875232699626384420, -0.005210224203359059109181555401330902819419),
            Complex64::new(0.01078377080978103125464543240346760257008, 0.006866888783433775382193630944275682670599),
            Complex64::new(-0.5808616819196736225612296471081337245459, 0.6688593905505562263387760667171706325749),
            Complex64::new(f64::INFINITY, f64::NEG_INFINITY),
            Complex64::new(0.1000052020902036118082966385855563526705e-7, 0.005100088434920073153418834680320146441685),
            Complex64::new(0.004950156837581592745389973960217444687524, -0.004899838305155226382584756154100963570500),
            Complex64::new(0.005100176864319675957314822982399286703798, 0.005099823128319785355949825238269336481254),
            Complex64::new(0.4244534840871830045021143490355372016428, 0.002820278933186814021399602648373095266538),
            Complex64::new(-0.1021340733271046543881236523269967674156, -0.00001045696456072005761498961861088944159916),
            Complex64::new(-0.01000200120119206748855061636187197886859, 0.9805885888237419500266621041508714123763e-8),
            Complex64::new(0.001000002000012000023960527532953151819595, -0.9800058800588007290937355024646722133204e-11),
            Complex64::new(0.4244549085628511778373438768121222815752, 0.002935393851311701428647152230552122898291),
            Complex64::new(-0.1021340732357117208743299813648493928105, -0.00001088377943049851799938998805451564893540),
            Complex64::new(-0.01000200120119126652710792390331206563616, 0.1020612612857282306892368985525393707486e-7),
            Complex64::new(0.1000000000007333333333344266666666664457e-5, 0.2000000000001333333333323199999999978819e-5),
            Complex64::new(0.1999999999994666666666675199999999990248e-5, 0.0),
            Complex64::new(0.3013403889237919660346644392864226952119, 0.0),
            Complex64::new(0.02503136792640367194699495234782353186858, 0.0),
            Complex64::new(0.002500031251171948248596912483183760683918, 0.0),
            Complex64::new(0.0, 0.004900078433419939164774792850907128053308),
            Complex64::new(0.0, -0.005100088434920074173454208832365950009419),
            Complex64::new(0.0, 0.2000000000005333333333341866666666676419e-5),
            Complex64::new(0.0, -48.16001211429122974789822893525016528191),
            Complex64::new(0.0, 0.4627407029504443513654142715903005954668e174),
            Complex64::new(0.0, f64::NEG_INFINITY),
            Complex64::new(0.0, 0.0),
            Complex64::new(-0.0, 0.0),
            Complex64::new(0.0, f64::INFINITY),
            Complex64::new(0.0, f64::NEG_INFINITY),
            Complex64::new(f64::NAN, f64::NAN),
            Complex64::new(f64::NAN, f64::NAN),
            Complex64::new(f64::NAN, f64::NAN),
            Complex64::new(f64::NAN, 0.0),
            Complex64::new(0.0, f64::NAN),
            Complex64::new(f64::NAN, f64::NAN),
            Complex64::new(f64::NAN, f64::NAN),
            Complex64::new(0.01282473148489433743567240624939698290584, -0.2105957276516618621447832572909153498104e-7),
            Complex64::new(0.01219875253423634378984109995893708152885, -0.1813040560401824664088425926165834355953e-7),
            Complex64::new(0.1020408163265306334945473399689037886997e-7, -0.1041232819658476285651490827866174985330e-25),
            Complex64::new(0.9803921568627452865036825956835185367356e-8, -0.9227220299884665067601095648451913375754e-26),
            Complex64::new(0.5000000000000000002500000000000000003750e-9, -0.1200000000000000001800000188712838420241e-29),
            Complex64::new(5.00000000000000000000025000000000000000000003e-12, -1.20000000000000000000018000000000000000000004e-36),
            Complex64::new(5.00000000000000000000000002500000000000000000e-14, -1.20000000000000000000000001800000000000000000e-42),
            Complex64::new(5.0e-301, 0.0) 
        ];

        let tolerance = 1.0e-13;
        for n in 0..z.len() {
            let computed_dawson = dawson_with_relerror(z[n], 0.);
            let relative_error_re = relerr(computed_dawson.re, expected_dawson[n].re);
            let relative_error_im = relerr(computed_dawson.im, expected_dawson[n].im);
            assert!(relative_error_re < tolerance);
            assert!(relative_error_im < tolerance);
        }
    }

} 

