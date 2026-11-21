use super::arch::*;

struct Rng(u64);
impl Rng {
    fn bytes(&mut self) -> [u8; 16] {
        let mut b = [0u8; 16];
        for x in b.iter_mut() {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            // mostly small values and values near the i8 / i16 limits
            *x = match self.0 >> 61 {
                0 => 0,
                1 => 0xFF,
                2 => 0x7F,
                3 => 0x80,
                _ => (self.0 >> 8) as u8,
            };
        }
        b
    }
}

unsafe fn v(b: [u8; 16]) -> V {
    load(b.as_ptr())
}
unsafe fn b(x: V) -> [u8; 16] {
    let mut o = [0u8; 16];
    store(o.as_mut_ptr(), x);
    o
}
fn w(b: [u8; 16]) -> [i16; 8] {
    core::array::from_fn(|i| i16::from_le_bytes([b[2 * i], b[2 * i + 1]]))
}
fn from_w(x: [i16; 8]) -> [u8; 16] {
    let mut o = [0u8; 16];
    for (i, l) in x.iter().enumerate() {
        o[2 * i..2 * i + 2].copy_from_slice(&l.to_le_bytes());
    }
    o
}
fn mask(t: bool) -> i16 {
    if t {
        -1
    } else {
        0
    }
}

#[test]
fn neon_arch_layer_matches_the_scalar_definitions() {
    let mut rng = Rng(0x2545_F491_4F6C_DD1D);
    let kcross: [u8; 16] = [
        128, 128, 128, 128, 128, 128, 128, 128, 6, 7, 6, 7, 6, 7, 6, 7,
    ];
    let kcross8: [u8; 16] = [
        128, 128, 128, 128, 128, 128, 128, 128, 7, 7, 7, 7, 7, 7, 7, 7,
    ];
    for _ in 0..20_000 {
        let (x, y, z) = (rng.bytes(), rng.bytes(), rng.bytes());
        let (xw, yw) = (w(x), w(y));
        // SAFETY: NEON is part of every aarch64 build (`arch` layer comment); the block passes
        // pointers only through `v` and `b`, which load and store one `[u8; 16]`.
        unsafe {
            let (vx, vy, vz) = (v(x), v(y), v(z));
            let m = x
                .iter()
                .enumerate()
                .map(|(i, &c)| u32::from(c >> 7) << i)
                .sum::<u32>();
            assert_eq!(mask8(vx), m, "mask8 {x:?}");
            let p: [u8; 16] = core::array::from_fn(|i| *x.get(z[i] as usize).unwrap_or(&0));
            assert_eq!(b(pick(vx, vz)), p, "pick");
            let t: Vec<u8> = x.iter().chain(y.iter()).copied().collect();
            let l: [u8; 16] = core::array::from_fn(|i| *t.get(z[i] as usize).unwrap_or(&0));
            assert_eq!(b(lut32(vx, vy, vz)), l, "lut32");
            for (n, got) in [(1, bsl1(vx)), (2, bsl2(vx)), (4, bsl4(vx)), (8, bsl8(vx))] {
                let e: [u8; 16] = core::array::from_fn(|i| if i >= n { x[i - n] } else { 0 });
                assert_eq!(b(got), e, "bsl{n}");
            }
            for (k, got) in [(8, shlq8(vx)), (16, shlq16(vx)), (32, shlq32(vx))] {
                let lo = u64::from_le_bytes(x[..8].try_into().unwrap()) << k;
                let hi = u64::from_le_bytes(x[8..].try_into().unwrap()) << k;
                let mut e = [0u8; 16];
                e[..8].copy_from_slice(&lo.to_le_bytes());
                e[8..].copy_from_slice(&hi.to_le_bytes());
                assert_eq!(b(got), e, "shlq{k}");
            }
            let e: [i16; 8] = core::array::from_fn(|i| xw[i].saturating_add(yw[i]));
            assert_eq!(b(adds16(vx, vy)), from_w(e), "adds16");
            let e: [i16; 8] = core::array::from_fn(|i| xw[i].wrapping_add(yw[i]));
            assert_eq!(b(add16(vx, vy)), from_w(e), "add16");
            let e: [i16; 8] =
                core::array::from_fn(|i| (xw[i] as u16).saturating_sub(yw[i] as u16) as i16);
            assert_eq!(b(subsu16(vx, vy)), from_w(e), "subsu16");
            let e: [i16; 8] = core::array::from_fn(|i| xw[i].max(yw[i]));
            assert_eq!(b(max16(vx, vy)), from_w(e), "max16");
            let e: [i16; 8] = core::array::from_fn(|i| mask(xw[i] > yw[i]));
            assert_eq!(b(gt16(vx, vy)), from_w(e), "gt16");
            let e: [u8; 16] =
                core::array::from_fn(|i| (x[i] as i8).saturating_add(y[i] as i8) as u8);
            assert_eq!(b(adds8(vx, vy)), e, "adds8");
            let e: [u8; 16] = core::array::from_fn(|i| x[i].wrapping_add(y[i]));
            assert_eq!(b(add8(vx, vy)), e, "add8");
            let e: [u8; 16] = core::array::from_fn(|i| x[i].saturating_sub(y[i]));
            assert_eq!(b(subsu8(vx, vy)), e, "subsu8");
            let e: [u8; 16] = core::array::from_fn(|i| (x[i] as i8).max(y[i] as i8) as u8);
            assert_eq!(b(max8(vx, vy)), e, "max8");
            let e: [u8; 16] =
                core::array::from_fn(|i| if (x[i] as i8) > (y[i] as i8) { 0xFF } else { 0 });
            assert_eq!(b(gt8(vx, vy)), e, "gt8");
            let e: [u8; 16] = core::array::from_fn(|i| x[i] & y[i]);
            assert_eq!(b(and(vx, vy)), e, "and");
            let e: [u8; 16] = core::array::from_fn(|i| x[i] | y[i]);
            assert_eq!(b(or(vx, vy)), e, "or");
            let e: [u8; 16] = core::array::from_fn(|i| !x[i] & y[i]);
            assert_eq!(b(andnot(vx, vy)), e, "andnot");
            let e: [u8; 16] =
                core::array::from_fn(|i| if i % 2 == 0 { x[i / 2] } else { y[i / 2] });
            assert_eq!(b(unpacklo8(vx, vy)), e, "unpacklo8");
            let e: [u8; 16] = core::array::from_fn(|i| {
                if i % 2 == 0 {
                    x[8 + i / 2]
                } else {
                    y[8 + i / 2]
                }
            });
            assert_eq!(b(unpackhi8(vx, vy)), e, "unpackhi8");
            let e: [u8; 16] = core::array::from_fn(|i| {
                let l = if i < 8 { xw[i] } else { yw[i - 8] };
                l.clamp(0, 255) as u8
            });
            assert_eq!(b(packus16(vx, vy)), e, "packus16");
            assert_eq!(b(bcast_last16(vx, vz)), from_w([xw[7]; 8]), "bcast_last16");
            assert_eq!(b(bcast_last8(vx, vz)), [x[15]; 16], "bcast_last8");
            let e: [i16; 8] = core::array::from_fn(|i| if i >= 4 { xw[3] } else { 0 });
            assert_eq!(b(cross16(vx, v(kcross))), from_w(e), "cross16");
            let e: [u8; 16] = core::array::from_fn(|i| if i >= 8 { x[7] } else { 0 });
            assert_eq!(b(cross8(vx, v(kcross8))), e, "cross8");
            assert_eq!(lane0(vx), xw[0] as i32, "lane0");
            assert_eq!(lane7(vx), xw[7] as i32, "lane7");
            assert_eq!(lane0_u8(vx), x[0] as i32, "lane0_u8");
            assert_eq!(lane15_u8(vx), x[15] as i32, "lane15_u8");
            assert_eq!(b(splat16(xw[2])), from_w([xw[2]; 8]), "splat16");
            assert_eq!(b(splat8(x[3] as i8)), [x[3]; 16], "splat8");
            assert_eq!(b(zero()), [0u8; 16], "zero");
        }
    }
}
