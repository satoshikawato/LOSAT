use super::ext::*;
use super::{load, store, V};

unsafe fn v(b: [u8; 16]) -> V {
    load(b.as_ptr())
}
unsafe fn b(x: V) -> [u8; 16] {
    let mut o = [0u8; 16];
    store(o.as_mut_ptr(), x);
    o
}
fn lanes<const N: usize, const W: usize>(b: [u8; 16]) -> [[u8; W]; N] {
    core::array::from_fn(|i| core::array::from_fn(|j| b[W * i + j]))
}

#[test]
fn neon_ext_layer_matches_the_scalar_definitions() {
    let mut s = 0x9E37_79B9_7F4A_7C15u64;
    let mut bytes = move || -> [u8; 16] {
        core::array::from_fn(|_| {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            match s >> 61 {
                0 => 0,
                1 => 0xFF,
                2 => 0x7F,
                3 => 0x80,
                _ => (s >> 8) as u8,
            }
        })
    };
    for _ in 0..20_000 {
        let (x, y) = (bytes(), bytes());
        let xw = lanes::<8, 2>(x).map(u16::from_le_bytes);
        let yw = lanes::<8, 2>(y).map(u16::from_le_bytes);
        let xd = lanes::<4, 4>(x).map(i32::from_le_bytes);
        let yd = lanes::<4, 4>(y).map(i32::from_le_bytes);
        let w16 =
            |l: [u16; 8]| -> [u8; 16] { core::array::from_fn(|i| l[i / 2].to_le_bytes()[i % 2]) };
        let w32 =
            |l: [i32; 4]| -> [u8; 16] { core::array::from_fn(|i| l[i / 4].to_le_bytes()[i % 4]) };
        // SAFETY: NEON is part of every aarch64 build (`arch` layer comment in `xdrop_simd.rs`)
        // and the `ext` functions take registers only; `v` and `b` touch one `[u8; 16]` each.
        unsafe {
            let (vx, vy) = (v(x), v(y));
            assert_eq!(
                b(max_u16(vx, vy)),
                w16(core::array::from_fn(|i| xw[i].max(yw[i]))),
                "max_u16"
            );
            assert_eq!(
                b(adds_u16(vx, vy)),
                w16(core::array::from_fn(|i| xw[i].saturating_add(yw[i]))),
                "adds_u16"
            );
            let e: [u8; 16] = core::array::from_fn(|i| x[i].max(y[i]));
            assert_eq!(b(max_u8(vx, vy)), e, "max_u8");
            let e: [u8; 16] = core::array::from_fn(|i| x[i].saturating_add(y[i]));
            assert_eq!(b(adds_u8(vx, vy)), e, "adds_u8");
            assert_eq!(any(vx), x.iter().any(|&c| c != 0), "any");
            let e: [u8; 16] = core::array::from_fn(|i| x[15 - i]);
            assert_eq!(b(rev8(vx)), e, "rev8");
            assert_eq!(b(splat32(xd[1])), w32([xd[1]; 4]), "splat32");
            let e: [i32; 4] = core::array::from_fn(|i| if xd[i] > yd[i] { -1 } else { 0 });
            assert_eq!(b(gt32(vx, vy)), w32(e), "gt32");
            let e: [u16; 8] = core::array::from_fn(|i| {
                let l = if i < 4 { xd[i] } else { yd[i - 4] };
                l.clamp(i16::MIN as i32, i16::MAX as i32) as i16 as u16
            });
            assert_eq!(b(narrow32(vx, vy)), w16(e), "narrow32");
            let e: [u8; 16] = core::array::from_fn(|i| {
                let l = if i < 8 { xw[i] } else { yw[i - 8] } as i16;
                l.clamp(i8::MIN as i16, i8::MAX as i16) as i8 as u8
            });
            assert_eq!(b(narrow16(vx, vy)), e, "narrow16");
            let e: [u8; 16] =
                core::array::from_fn(|i| if i < 8 { x[2 * i] } else { x[2 * (i - 8) + 1] });
            assert_eq!(b(split16(vx)), e, "split16");
            let e: [u8; 16] = core::array::from_fn(|i| if i < 8 { x[i] } else { y[i - 8] });
            assert_eq!(b(lo64(vx, vy)), e, "lo64");
            let e: [u8; 16] = core::array::from_fn(|i| if i < 8 { x[8 + i] } else { y[i] });
            assert_eq!(b(hi64(vx, vy)), e, "hi64");
            for (n, got) in [(1, up1b(vx)), (2, up2b(vx)), (4, up4b(vx)), (8, up8b(vx))] {
                let e: [u8; 16] = core::array::from_fn(|i| if i >= n { x[i - n] } else { 0 });
                assert_eq!(b(got), e, "up{n}b");
            }
        }
    }
}
