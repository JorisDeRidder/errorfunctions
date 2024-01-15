// Run with: cargo run --release --example example

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
