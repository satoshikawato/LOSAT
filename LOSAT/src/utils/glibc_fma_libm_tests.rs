use super::variants as v;
use super::{exp, ln, log10, pow};

/// (function, x, y, result) as bits, function 0 `exp`, 1 `ln`, 2 `log10`, 3 `pow`: arguments
/// where glibc 2.35 and 2.39 on x86-64 with FMA give `result` and glibc with its FMA variants
/// off (the same algorithm as musl's, so WASI's) gives other bits.
#[rustfmt::skip]
const PINNED: [(u8, u64, u64, u64); 80] = [
    (0, 0xbfb7346e152d2b50, 0x0000000000000000, 0x3fed3a1b212e1efd), // exp[-1,1]
    (0, 0x3fe95f2fc2688c1e, 0x0000000000000000, 0x4001ad85e637996d), // exp[-1,1]
    (0, 0xbfe9b99b5117812a, 0x0000000000000000, 0x3fdca51949f00589), // exp[-1,1]
    (0, 0xc008608684fe4de0, 0x0000000000000000, 0x3fa8514291f81ddb), // exp[-60,0]
    (0, 0xc03ea6da0f9bb927, 0x0000000000000000, 0x3d2b73c49f291119), // exp[-60,0]
    (0, 0xc04a250d7859b4ac, 0x0000000000000000, 0x3b379ff8df28f1c1), // exp[-60,0]
    (0, 0xc04b7466a52670f0, 0x0000000000000000, 0x3afb8598c6df27a7), // exp[-745,710]
    (0, 0xc04e5001b3d27720, 0x0000000000000000, 0x3a773532ba6abab1), // exp[-745,710]
    (0, 0x408292cbb12e1cd6, 0x0000000000000000, 0x758615d8d0597efd), // exp[-745,710]
    (0, 0xc086253befc3fa92, 0x0000000000000000, 0x000c5d08786e7dd3), // exp[-745,-708]
    (0, 0xc086247de02b3c29, 0x0000000000000000, 0x000d90d1b2f6a67f), // exp[-745,-708]
    (0, 0xc08630542322e779, 0x0000000000000000, 0x000316d8552aba69), // exp[-745,-708]
    (0, 0xc005af8431352386, 0x0000000000000000, 0x3fb1057fba5e6c13), // exp[bits]
    (0, 0xc04ec2ddd59eb197, 0x0000000000000000, 0x3a62ebe01ba5b68b), // exp[bits]
    (0, 0x404456c5103f2292, 0x0000000000000000, 0x4399bcbf6fa78141), // exp[bits]
    (1, 0x3ff173bfb5d375f0, 0x0000000000000000, 0x3fb63d5c42031633), // ln[0.5,2]
    (1, 0x3ff2203196002bec, 0x0000000000000000, 0x3fbfef4d89eeeb90), // ln[0.5,2]
    (1, 0x3ff50394db9b6778, 0x0000000000000000, 0x3fd17245f941fc78), // ln[0.5,2]
    (1, 0x3feee9f1015f2ff7, 0x0000000000000000, 0xbfa1ae31e42e66fb), // ln[near1]
    (1, 0x3ff069c9c9d1c7db, 0x0000000000000000, 0x3f9a1c7e3882fb17), // ln[near1]
    (1, 0x3ff09fa57208a91e, 0x0000000000000000, 0x3fa393a2cd7cabbf), // ln[near1]
    (1, 0x413fc35fc9b9d5f0, 0x0000000000000000, 0x402d18ea4fd58f10), // ln[1,1e7]
    (1, 0x40f9633d3eacb7e9, 0x0000000000000000, 0x40271aa38fef55bb), // ln[1,1e7]
    (1, 0x40f206ccd46ef0cf, 0x0000000000000000, 0x40266b52564646c8), // ln[1,1e7]
    (1, 0x729b79abeb172e68, 0x0000000000000000, 0x40818feb5823828d), // ln[posbits]
    (1, 0x3d61dd35777df42b, 0x0000000000000000, 0xc03c4f0f7e9132e0), // ln[posbits]
    (1, 0x3bb38ab5519c1979, 0x0000000000000000, 0xc047778e9d465d6d), // ln[posbits]
    (1, 0x000976060b0abca6, 0x0000000000000000, 0xc086275fe65325c5), // ln[subnormal]
    (1, 0x3fece4e6c6ce0e67, 0x0000000000000000, 0xbfba2305570b4413), // ln[bits]
    (1, 0x3b89ce2fc16490f6, 0x0000000000000000, 0xc0485e235c8054cc), // ln[bits]
    (1, 0x54fc61b6b5fcf89e, 0x0000000000000000, 0x406d2f0f392a0aa9), // ln[bits]
    (2, 0x401e4458d91a5290, 0x0000000000000000, 0x3fec2005febfa766), // log10[posbits]
    (2, 0x3ffa3fe2cec11f1b, 0x0000000000000000, 0x3fcb852f1f8b31de), // log10[posbits]
    (2, 0x38083b5e0cf0cd5c, 0x0000000000000000, 0xc04306782824ebdc), // log10[posbits]
    (2, 0x3fe7684d532711f3, 0x0000000000000000, 0xbfc161c57c4cf001), // log10[1e-300,1e3]
    (2, 0x3ff97f340663b4e1, 0x0000000000000000, 0x3fc9e72b000ac001), // log10[1e-300,1e3]
    (2, 0x3fe6835ae0cdc81e, 0x0000000000000000, 0xbfc38c26b0345000), // log10[1e-300,1e3]
    (2, 0x3fb894eb91f41dc5, 0x0000000000000000, 0xbff04833fd669cf2), // log10[bits]
    (2, 0x44a334e05ee3f129, 0x0000000000000000, 0x4036a8159013bc5f), // log10[bits]
    (2, 0x40727763e100434c, 0x0000000000000000, 0x4003c39645736ba2), // log10[bits]
    (3, 0x40135012c00000dc, 0x3ff8bce027903aa4, 0x4026d0e1007ec7a7), // pow[0,10]x[-5,5]
    (3, 0x401e613a32570e30, 0x4005d8a80d052564, 0x406fba46c387441f), // pow[0,10]x[-5,5]
    (3, 0x4023da2f44fa623a, 0x400c131266523fbc, 0x40a897b395667ccf), // pow[0,10]x[-5,5]
    (3, 0x400ea5eb0693542b, 0xbfa1b636d3a50080, 0x3fee8c14724b6659), // pow[0,10]x[-5,5]
    (3, 0x3ff2bd6cd5900238, 0x3fe0f1193b7fb1d6, 0x3ff16587e3f39be7), // pow[0,4]x[0,1]
    (3, 0x400d2b0e306f1031, 0x3fe4e31006bda2ac, 0x40029cbdec3e1a9d), // pow[0,4]x[0,1]
    (3, 0x3ffc9ab45a50f7fa, 0x3fc3c0db53373ff8, 0x3ff180334679ac79), // pow[0,4]x[0,1]
    (3, 0x4000372b47465044, 0x3fae92b556245c80, 0x3ff0b080a4bcc283), // pow[0,4]x[0,1]
    (3, 0x4024000000000000, 0x405f800000000000, 0x5a17a2ecc414a03f), // pow10^int
    (3, 0x4024000000000000, 0x405f800000000000, 0x5a17a2ecc414a03f), // pow10^int
    (3, 0x4024000000000000, 0x405f800000000000, 0x5a17a2ecc414a03f), // pow10^int
    (3, 0x4024000000000000, 0x405f800000000000, 0x5a17a2ecc414a03f), // pow10^int
    (3, 0x3feffda9549186e8, 0xc0f15ccfe8c6db48, 0x41c39b28f0000f03), // pow[near1]x[big]
    (3, 0x3feffb0c1bcfdce7, 0x4117e74593989f50, 0x2a93c6d03463cbcf), // pow[near1]x[big]
    (3, 0x3fefff84d6b72c19, 0x40e2d94b6f0aed90, 0x3fba86534d13fc21), // pow[near1]x[big]
    (3, 0x3ff00242ac054991, 0x412adc2e0866e670, 0x6bb7909b8b62dc17), // pow[near1]x[big]
    (3, 0xbffd5aafcb95bcf2, 0xc058c00000000000, 0xba840da9e7914d69), // pow[neg]x[int]
    (3, 0xbff1a261f339bae2, 0x4024000000000000, 0x400528713c1e5757), // pow[neg]x[int]
    (3, 0xbff60ac87808c33a, 0xc046800000000000, 0xbea265a0852bf0db), // pow[neg]x[int]
    (3, 0xbffac9c1474afc0a, 0x4008000000000000, 0xc012c5cdb12e4433), // pow[neg]x[int]
    (3, 0x300a8d6a7a6a0dd1, 0xbfdc074dc018c9e0, 0x46e47adabdcedd37), // pow[posbits]x[-2,2]
    (3, 0x51ff461ae6bf320e, 0x3fe089adcdf84630, 0x4944406bad06ff6d), // pow[posbits]x[-2,2]
    (3, 0x54b56930780fc9c6, 0xbff18c4d7623a760, 0x29255c51b09dad81), // pow[posbits]x[-2,2]
    (3, 0x58814ada472ea9ea, 0xbfe20860cdac6170, 0x321632ed94d182f7), // pow[posbits]x[-2,2]
    (3, 0x0002a8b242943bb4, 0x3fd3b506dc3c7710, 0x2c36b26193f158e9), // pow[subn]x[-1,1]
    (3, 0x00080d8a74ecc080, 0xbf9688e2935d3240, 0x4156d2fdd889844b), // pow[subn]x[-1,1]
    (3, 0x000ec724df16f69c, 0xbfd2decfcea072a0, 0x52c4acde4ad4d013), // pow[subn]x[-1,1]
    (3, 0x000ff1f02d7edff4, 0xbfe9921dbff5bee2, 0x72f97678a909b4f9), // pow[subn]x[-1,1]
    (3, 0x46698bef4eb624f5, 0x3f9cf6d49d358f15, 0x401e89750613d523), // pow[bits]x[bits]
    (3, 0x24df72154c46d6ec, 0x3ff136ff4a783591, 0x22d11b37cb48b749), // pow[bits]x[bits]
    (3, 0x2801c12f0e55700f, 0x3f34b1cbaa38d9c1, 0x3fed6d8242cc45a3), // pow[bits]x[bits]
    (3, 0x1821408f368c1705, 0xbf4a699f2074e54e, 0x3ff6d687cd32677d), // pow[bits]x[bits]
    (3, 0x3ff9f95f9d1b139c, 0x40921fcc90e0ec26, 0x729c035093405175), // pow[ovf]
    (3, 0x3ffd969a00f66d59, 0x4090d874d00a9f83, 0x7bb2de524f617177), // pow[ovf]
    (3, 0x3ffbc0a346964152, 0x4091bbfc0a276b1c, 0x784bf8f8a0a60da1), // pow[ovf]
    (3, 0x3ffb5c4386785a90, 0x4093fc0d93a96f72, 0x7dcf890f34c08f57), // pow[ovf]
    (3, 0x3fda8ee00ea9a966, 0x407b161b1bebc578, 0x1d90d6437bd90ae9), // pow[unf]
    (3, 0x3fd87ccca77b5bc0, 0x407ea46956b476f2, 0x1576152909f85fb6), // pow[unf]
    (3, 0x3fdeab6e8bc8116e, 0x407d49fc164ceaee, 0x20d96c56eab1fe33), // pow[unf]
    (3, 0x3fe96ebcd599258e, 0x407b1665889bbd2c, 0x36f4bab57ab5f09f), // pow[unf]
];

fn call(f: u8, x: f64, y: f64) -> [f64; 3] {
    match f {
        0 => [exp(x), v::exp_emulated(x), v::exp_fused(x)],
        1 => [ln(x), v::ln_emulated(x), v::ln_fused(x)],
        2 => [log10(x), v::log10_emulated(x), v::log10_fused(x)],
        _ => [pow(x, y), v::pow_emulated(x, y), v::pow_fused(x, y)],
    }
}

fn bits(x: f64) -> u64 {
    if x.is_nan() {
        0x7ff8_0000_0000_0000
    } else {
        x.to_bits()
    }
}

#[test]
fn every_target_writes_glibcs_fma_bits_where_other_libms_differ() {
    for (f, x, y, want) in PINNED {
        let (x, y) = (f64::from_bits(x), f64::from_bits(y));
        for (i, got) in call(f, x, y).into_iter().enumerate() {
            assert_eq!(got.to_bits(), want, "function {f} x={x:e} y={y:e} way {i}");
        }
    }
}

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    fn unif(&mut self, a: f64, b: f64) -> f64 {
        a + (b - a) * ((self.next() >> 11) as f64 / (1u64 << 53) as f64)
    }
}

/// Arguments of every range the functions treat apart: (function, x, y).
/// 2^-4 to 2^-53: the distance from 1 of an argument of `log`'s interval near 1.
fn near_one(r: &mut Rng) -> f64 {
    2f64.powi(-4 - (r.next() % 50) as i32)
}

fn arguments(n: usize) -> Vec<(u8, f64, f64)> {
    let mut r = Rng(12345);
    let mut out = Vec::new();
    for _ in 0..n {
        let any = f64::from_bits(r.next());
        out.extend([
            (0, r.unif(-1.0, 1.0), 0.0),
            (0, r.unif(-60.0, 0.0), 0.0),
            (0, r.unif(-745.2, 709.8), 0.0),
            (0, r.unif(-745.2, -708.0), 0.0),
            (0, any, 0.0),
            (0, r.unif(-1.0, 1.0) * 1e-17, 0.0),
            (1, r.unif(0.5, 2.0), 0.0),
            (1, r.unif(0.937, 1.0647), 0.0),
            (1, 1.0 + r.unif(-1.0, 1.0) * near_one(&mut r), 0.0),
            (1, r.unif(1.0, 1e7), 0.0),
            (1, f64::from_bits(r.next() >> 1), 0.0),
            (1, f64::from_bits(r.next() >> 12), 0.0),
            (1, any, 0.0),
            (2, f64::from_bits(r.next() >> 1), 0.0),
            (2, any, 0.0),
            (3, r.unif(0.0, 10.0), r.unif(-5.0, 5.0)),
            (3, 10.0, ((r.next() % 700) as i32 - 400) as f64),
            (3, r.unif(0.999, 1.001), r.unif(-1e6, 1e6)),
            (3, -r.unif(0.0, 3.0), ((r.next() % 200) as i64 - 100) as f64),
            (3, f64::from_bits(r.next() >> 1), r.unif(-2.0, 2.0)),
            (3, f64::from_bits(r.next() >> 12), r.unif(-1.0, 1.0)),
            (3, any, f64::from_bits(r.next())),
            (3, r.unif(1.5, 3.0), r.unif(900.0, 1900.0)),
            (3, r.unif(0.1, 0.9), r.unif(300.0, 800.0)),
        ]);
    }
    let special = [
        0.0,
        -0.0,
        1.0,
        -1.0,
        2.0,
        -2.0,
        0.5,
        10.0,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NAN,
        709.782712893384,
        709.7827128933841,
        -745.1332191019411,
        -745.1332191019412,
        -708.3964185322641,
        1e-300,
        5e-324,
        f64::MAX,
        f64::MIN_POSITIVE,
        3.0,
        -3.0,
        0.1,
        1e300,
        1024.0,
        -1075.0,
        1e-20,
        f64::from_bits(0x3be0_0000_0000_0000),
        f64::from_bits(0x43e0_0000_0000_0000),
    ];
    for x in special {
        out.extend([(0, x, 0.0), (1, x, 0.0), (2, x, 0.0)]);
        for y in special {
            out.push((3, x, y));
        }
    }
    out
}

#[test]
fn the_emulated_fma_way_gives_the_fused_way_bits() {
    for (f, x, y) in arguments(20_000) {
        let [public, emulated, fused] = call(f, x, y);
        assert_eq!(bits(emulated), bits(fused), "function {f} x={x:e} y={y:e}");
        assert_eq!(bits(public), bits(fused), "function {f} x={x:e} y={y:e}");
    }
}

#[test]
fn the_emulated_fma_rounds_once() {
    let mut r = Rng(99);
    for i in 0..1_000_000 {
        let ea = (r.next() % 120) as i32 - 60;
        let eb = (r.next() % 120) as i32 - 60;
        let sign = |r: &mut Rng| if r.next() & 1 == 0 { 1.0 } else { -1.0 };
        let a = r.unif(1.0, 2.0) * 2f64.powi(ea) * sign(&mut r);
        let b = r.unif(1.0, 2.0) * 2f64.powi(eb) * sign(&mut r);
        let ph = a * b;
        let c = match i % 5 {
            0 => -ph,
            1 => -ph + r.unif(-1.0, 1.0) * ph.abs() * 2f64.powi(-((r.next() % 60) as i32)),
            2 => r.unif(-4.0, 4.0) * 2f64.powi(ea + eb + (r.next() % 140) as i32 - 70),
            3 => {
                // ph + c on a midpoint of the sum's precision
                let s = ph * r.unif(1.0, 3.0);
                let u = f64::from_bits(s.to_bits() + 1) - s;
                s - ph + u * 0.5
            }
            _ => a.mul_add(b, -ph) * r.unif(-3.0, 3.0) - ph * 0.5,
        };
        let want = bits(a.mul_add(b, c));
        assert_eq!(
            bits(v::fma_emulated(a, b, c)),
            want,
            "a={a:e} b={b:e} c={c:e}"
        );
        assert_eq!(
            bits(v::fma_emulated_full(a, b, c)),
            want,
            "a={a:e} b={b:e} c={c:e}"
        );
    }
}

#[cfg(all(target_os = "linux", target_env = "gnu", target_arch = "x86_64"))]
#[test]
#[expect(
    clippy::disallowed_methods,
    reason = "compares the port with glibc itself"
)]
fn glibc_with_fma_gives_the_same_bits() {
    if !std::arch::is_x86_feature_detected!("fma") {
        return; // glibc then takes its variants without FMA
    }
    for (f, x, y) in arguments(20_000) {
        let glibc = match f {
            0 => x.exp(),
            1 => x.ln(),
            2 => x.log10(),
            _ => x.powf(y),
        };
        assert_eq!(
            bits(call(f, x, y)[0]),
            bits(glibc),
            "function {f} x={x:e} y={y:e}"
        );
    }
}
