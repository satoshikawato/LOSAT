use super::{exp_impl, ln_impl, log10_impl, pow_impl, Emulated, Fma, Fused};

pub fn fma_emulated(a: f64, b: f64, c: f64) -> f64 {
    Emulated::fma(a, b, c)
}
pub fn fma_emulated_full(a: f64, b: f64, c: f64) -> f64 {
    Emulated::full(a, b, c)
}
pub fn exp_emulated(x: f64) -> f64 {
    exp_impl::<Emulated>(x)
}
pub fn ln_emulated(x: f64) -> f64 {
    ln_impl::<Emulated>(x)
}
pub fn log10_emulated(x: f64) -> f64 {
    log10_impl::<Emulated>(x)
}
pub fn pow_emulated(x: f64, y: f64) -> f64 {
    pow_impl::<Emulated>(x, y)
}
pub fn exp_fused(x: f64) -> f64 {
    exp_impl::<Fused>(x)
}
pub fn ln_fused(x: f64) -> f64 {
    ln_impl::<Fused>(x)
}
pub fn log10_fused(x: f64) -> f64 {
    log10_impl::<Fused>(x)
}
pub fn pow_fused(x: f64, y: f64) -> f64 {
    pow_impl::<Fused>(x, y)
}
