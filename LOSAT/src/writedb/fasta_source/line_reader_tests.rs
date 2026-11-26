use super::*;

/// `find_eol` is the first `\n` or `\r` (the byte-by-byte scan it replaced) for every byte at every
/// place of buffers across the blocks of the search, and on random buffers up to a whole
/// `BUFFER_SIZE` with line ends far apart, close together and absent.
#[test]
fn find_eol_is_the_first_line_end() {
    let scan = |bytes: &[u8]| {
        bytes
            .iter()
            .position(|&byte| byte == b'\n' || byte == b'\r')
    };
    for byte in 0..=u8::MAX {
        for len in [1, 2, 31, 32, 33, 61, 63, 64, 65, 127, 128, 129, 200, 300] {
            for at in 0..len {
                let mut bytes = vec![b'A'; len];
                bytes[at] = byte;
                assert_eq!(find_eol(&bytes), scan(&bytes), "{byte:#04x} {len} {at}");
            }
        }
    }
    let mut state: u64 = 0x2545_f491_4f6c_dd1d;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    for round in 0..200 {
        let len = match round % 4 {
            0 => (next() % 100) as usize,
            1 => (next() % 1000) as usize,
            _ => (next() % BUFFER_SIZE as u64) as usize,
        };
        let sparse = [8u64, 80, 4000, u64::MAX][round % 4];
        let bytes: Vec<u8> = (0..len)
            .map(|_| match next() % sparse {
                0 => [b'\n', b'\r'][(next() % 2) as usize],
                _ => next() as u8 | 0x20,
            })
            .collect();
        for from in [0, len / 3, len / 2] {
            assert_eq!(
                find_eol(&bytes[from..]),
                scan(&bytes[from..]),
                "{round} {from}"
            );
        }
    }
}
