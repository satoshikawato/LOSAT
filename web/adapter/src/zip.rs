//! A ZIP archive written as it goes (PKWARE APPNOTE 6.3.10): the outputs of outfmt 13 and 14, a
//! master file and one side file per report, as one file on stream 4 (docs/web/abi_v2.md §5). The
//! entries are stored (method 0, no compression), each with its CRC-32 and sizes in its local
//! header, and the same modification time (1980-01-01 00:00, the first time that MS-DOS dates
//! hold), so that the archive of a run is the same bytes every time. ZIP64 records come only where
//! the classic fields cannot hold a value: 65,535 entries or more, or an offset or the central
//! directory's size of 4 GiB - 1 or more.

use std::io::{self, Write};

const LOCAL_HEADER: u32 = 0x0403_4b50;
const CENTRAL_HEADER: u32 = 0x0201_4b50;
const END_OF_CENTRAL_DIRECTORY: u32 = 0x0605_4b50;
const ZIP64_END_OF_CENTRAL_DIRECTORY: u32 = 0x0606_4b50;
const ZIP64_END_LOCATOR: u32 = 0x0706_4b50;
/// The version needed to extract a stored entry (2.0), and one with ZIP64 fields (4.5).
const VERSION: u16 = 20;
const VERSION_ZIP64: u16 = 45;
/// 1980-01-01 in the MS-DOS date format (year - 1980, month, day), at 00:00:00.
const DOS_DATE: u16 = (1 << 5) | 1;
const DOS_TIME: u16 = 0;

/// The CRC-32 of ZIP (ISO-HDLC: the reflected polynomial 0xEDB88320).
fn crc32(bytes: &[u8]) -> u32 {
    const TABLE: [u32; 256] = {
        let mut table = [0u32; 256];
        let mut n = 0;
        while n < 256 {
            let mut c = n as u32;
            let mut k = 0;
            while k < 8 {
                c = if c & 1 != 0 {
                    0xEDB8_8320 ^ (c >> 1)
                } else {
                    c >> 1
                };
                k += 1;
            }
            table[n] = c;
            n += 1;
        }
        table
    };
    let mut crc = !0u32;
    for &byte in bytes {
        crc = TABLE[((crc ^ u32::from(byte)) & 0xFF) as usize] ^ (crc >> 8);
    }
    !crc
}

/// An entry written: what its central directory header repeats.
struct Entry {
    name: Vec<u8>,
    crc: u32,
    size: u32,
    offset: u64,
}

/// A ZIP archive that writes each entry to `out` as it is added and the central directory at
/// `finish`.
pub struct ZipStream<W: Write> {
    out: W,
    written: u64,
    entries: Vec<Entry>,
}

impl<W: Write> ZipStream<W> {
    pub fn new(out: W) -> Self {
        Self {
            out,
            written: 0,
            entries: Vec::new(),
        }
    }

    fn put(&mut self, bytes: &[u8]) -> io::Result<()> {
        self.out.write_all(bytes)?;
        self.written += bytes.len() as u64;
        Ok(())
    }

    /// Adds the file `name` with the contents `bytes` (stored). A file of 4 GiB - 1 bytes or more
    /// (a size that needs ZIP64's) is refused: a browser's memory cannot hold one.
    pub fn add(&mut self, name: &[u8], bytes: &[u8]) -> io::Result<()> {
        let size = u32::try_from(bytes.len())
            .ok()
            .filter(|&size| size != u32::MAX)
            .ok_or_else(|| io::Error::other("a file of 4 GiB or more in a ZIP archive"))?;
        let name_len = u16::try_from(name.len())
            .map_err(|_| io::Error::other("a file name of 64 KiB or more in a ZIP archive"))?;
        let crc = crc32(bytes);
        let offset = self.written;
        let mut header = Vec::with_capacity(30 + name.len());
        header.extend_from_slice(&LOCAL_HEADER.to_le_bytes());
        // The version that its central header gives (ZIP64 for an offset that needs it).
        let version = if offset >= u64::from(u32::MAX) {
            VERSION_ZIP64
        } else {
            VERSION
        };
        header.extend_from_slice(&version.to_le_bytes());
        header.extend_from_slice(&0u16.to_le_bytes()); // flags
        header.extend_from_slice(&0u16.to_le_bytes()); // stored
        header.extend_from_slice(&DOS_TIME.to_le_bytes());
        header.extend_from_slice(&DOS_DATE.to_le_bytes());
        header.extend_from_slice(&crc.to_le_bytes());
        header.extend_from_slice(&size.to_le_bytes()); // compressed
        header.extend_from_slice(&size.to_le_bytes()); // uncompressed
        header.extend_from_slice(&name_len.to_le_bytes());
        header.extend_from_slice(&0u16.to_le_bytes()); // extra field
        header.extend_from_slice(name);
        self.put(&header)?;
        self.put(bytes)?;
        self.entries.push(Entry {
            name: name.to_vec(),
            crc,
            size,
            offset,
        });
        Ok(())
    }

    /// Writes the central directory and the end records, and gives back the writer.
    pub fn finish(mut self) -> io::Result<W> {
        let start = self.written;
        let entries = std::mem::take(&mut self.entries);
        for entry in &entries {
            let zip64 = entry.offset >= u64::from(u32::MAX);
            let mut header = Vec::with_capacity(46 + entry.name.len() + 12);
            header.extend_from_slice(&CENTRAL_HEADER.to_le_bytes());
            let version = if zip64 { VERSION_ZIP64 } else { VERSION };
            header.extend_from_slice(&version.to_le_bytes()); // made by (MS-DOS)
            header.extend_from_slice(&version.to_le_bytes()); // needed
            header.extend_from_slice(&0u16.to_le_bytes()); // flags
            header.extend_from_slice(&0u16.to_le_bytes()); // stored
            header.extend_from_slice(&DOS_TIME.to_le_bytes());
            header.extend_from_slice(&DOS_DATE.to_le_bytes());
            header.extend_from_slice(&entry.crc.to_le_bytes());
            header.extend_from_slice(&entry.size.to_le_bytes());
            header.extend_from_slice(&entry.size.to_le_bytes());
            header.extend_from_slice(&(entry.name.len() as u16).to_le_bytes());
            header.extend_from_slice(&(if zip64 { 12u16 } else { 0 }).to_le_bytes());
            header.extend_from_slice(&0u16.to_le_bytes()); // comment
            header.extend_from_slice(&0u16.to_le_bytes()); // disk
            header.extend_from_slice(&0u16.to_le_bytes()); // internal attributes
            header.extend_from_slice(&0u32.to_le_bytes()); // external attributes
            let offset = if zip64 { u32::MAX } else { entry.offset as u32 };
            header.extend_from_slice(&offset.to_le_bytes());
            header.extend_from_slice(&entry.name);
            if zip64 {
                // The ZIP64 extended information: the local header's offset only.
                header.extend_from_slice(&1u16.to_le_bytes());
                header.extend_from_slice(&8u16.to_le_bytes());
                header.extend_from_slice(&entry.offset.to_le_bytes());
            }
            self.put(&header)?;
        }
        let size = self.written - start;
        let count = entries.len() as u64;
        let zip64 = count >= u64::from(u16::MAX)
            || start >= u64::from(u32::MAX)
            || size >= u64::from(u32::MAX);
        if zip64 {
            let at = self.written;
            let mut record = Vec::with_capacity(56 + 20);
            record.extend_from_slice(&ZIP64_END_OF_CENTRAL_DIRECTORY.to_le_bytes());
            record.extend_from_slice(&44u64.to_le_bytes()); // the size of the rest
            record.extend_from_slice(&VERSION_ZIP64.to_le_bytes());
            record.extend_from_slice(&VERSION_ZIP64.to_le_bytes());
            record.extend_from_slice(&0u32.to_le_bytes()); // this disk
            record.extend_from_slice(&0u32.to_le_bytes()); // the central directory's disk
            record.extend_from_slice(&count.to_le_bytes());
            record.extend_from_slice(&count.to_le_bytes());
            record.extend_from_slice(&size.to_le_bytes());
            record.extend_from_slice(&start.to_le_bytes());
            record.extend_from_slice(&ZIP64_END_LOCATOR.to_le_bytes());
            record.extend_from_slice(&0u32.to_le_bytes()); // the disk of the record
            record.extend_from_slice(&at.to_le_bytes());
            record.extend_from_slice(&1u32.to_le_bytes()); // disks
            self.put(&record)?;
        }
        let small = |value: u64, max: u64| if zip64 { max } else { value };
        let mut end = Vec::with_capacity(22);
        end.extend_from_slice(&END_OF_CENTRAL_DIRECTORY.to_le_bytes());
        end.extend_from_slice(&0u16.to_le_bytes()); // this disk
        end.extend_from_slice(&0u16.to_le_bytes()); // the central directory's disk
        let entries16 = small(count, u64::from(u16::MAX)) as u16;
        end.extend_from_slice(&entries16.to_le_bytes());
        end.extend_from_slice(&entries16.to_le_bytes());
        end.extend_from_slice(&(small(size, u64::from(u32::MAX)) as u32).to_le_bytes());
        end.extend_from_slice(&(small(start, u64::from(u32::MAX)) as u32).to_le_bytes());
        end.extend_from_slice(&0u16.to_le_bytes()); // comment
        self.put(&end)?;
        Ok(self.out)
    }

    /// The writer, without the central directory (a run that failed).
    pub fn into_inner(self) -> W {
        self.out
    }
}

/// The entries of a stored ZIP archive as (name, bytes), read through its central directory, with
/// every CRC-32 and size checked (the adapter's tests).
#[cfg(test)]
pub fn read(archive: &[u8]) -> Vec<(Vec<u8>, Vec<u8>)> {
    let u16_at = |at: usize| u16::from_le_bytes(archive[at..at + 2].try_into().unwrap());
    let u32_at = |at: usize| u32::from_le_bytes(archive[at..at + 4].try_into().unwrap());
    let u64_at = |at: usize| u64::from_le_bytes(archive[at..at + 8].try_into().unwrap());
    let end = archive.len() - 22;
    assert_eq!(u32_at(end), END_OF_CENTRAL_DIRECTORY);
    let (mut count, mut start) = (u64::from(u16_at(end + 10)), u64::from(u32_at(end + 16)));
    if count == u64::from(u16::MAX) {
        let locator = end - 20;
        assert_eq!(u32_at(locator), ZIP64_END_LOCATOR);
        let record = u64_at(locator + 8) as usize;
        assert_eq!(u32_at(record), ZIP64_END_OF_CENTRAL_DIRECTORY);
        (count, start) = (u64_at(record + 32), u64_at(record + 48));
    }
    let mut at = start as usize;
    let mut files = Vec::new();
    for _ in 0..count {
        assert_eq!(u32_at(at), CENTRAL_HEADER);
        let (crc, size) = (u32_at(at + 16), u32_at(at + 24) as usize);
        let (name_len, extra_len) = (u16_at(at + 28) as usize, u16_at(at + 30) as usize);
        let mut offset = u64::from(u32_at(at + 42));
        let name = archive[at + 46..at + 46 + name_len].to_vec();
        if offset == u64::from(u32::MAX) {
            offset = u64_at(at + 46 + name_len + 4);
        }
        let local = offset as usize;
        assert_eq!(u32_at(local), LOCAL_HEADER);
        assert_eq!(u32_at(local + 14), crc);
        assert_eq!(u32_at(local + 22) as usize, size);
        let data = local + 30 + u16_at(local + 26) as usize + u16_at(local + 28) as usize;
        let bytes = archive[data..data + size].to_vec();
        assert_eq!(crc32(&bytes), crc);
        files.push((name, bytes));
        at += 46 + name_len + extra_len;
    }
    files
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crc32_is_zips() {
        // The check value of CRC-32/ISO-HDLC, and zlib.crc32(b"") and (b"a").
        assert_eq!(crc32(b"123456789"), 0xCBF4_3926);
        assert_eq!(crc32(b""), 0);
        assert_eq!(crc32(b"a"), 0xE8B7_BE43);
    }

    // An archive of the same files is the same bytes; its entries read back with their names and
    // contents; an empty archive is the end record alone.
    #[test]
    fn entries_read_back_and_the_bytes_are_fixed() {
        let make = || {
            let mut zip = ZipStream::new(Vec::new());
            zip.add(b"r_1.json", b"{\n}\n").unwrap();
            zip.add(b"r_2.json", b"").unwrap();
            zip.add(b"r.json", b"{\"File\": \"r_1.json\"}").unwrap();
            zip.finish().unwrap()
        };
        let archive = make();
        assert_eq!(archive, make());
        assert_eq!(
            read(&archive),
            [
                (b"r_1.json".to_vec(), b"{\n}\n".to_vec()),
                (b"r_2.json".to_vec(), Vec::new()),
                (b"r.json".to_vec(), b"{\"File\": \"r_1.json\"}".to_vec()),
            ]
        );
        assert_eq!(&archive[..4], b"PK\x03\x04");
        let empty = ZipStream::new(Vec::new()).finish().unwrap();
        assert_eq!(empty.len(), 22);
        assert!(read(&empty).is_empty());
    }

    // An entry whose local header starts at 4 GiB - 1 or later has its offset in a ZIP64 extra field
    // of its central header, and the central directory there needs the ZIP64 end records (an
    // archive that starts 16 bytes before that point; only its own bytes are kept).
    #[test]
    fn offsets_of_4_gib_take_zip64_fields() {
        let base = u64::from(u32::MAX) - 16;
        let mut zip = ZipStream {
            out: Vec::new(),
            written: base,
            entries: Vec::new(),
        };
        zip.add(b"a", b"").unwrap(); // at base: 31 bytes
        zip.add(b"b", b"xyz").unwrap(); // at base + 31
        let bytes = zip.finish().unwrap();
        let u16_at = |at: usize| u16::from_le_bytes(bytes[at..at + 2].try_into().unwrap());
        let u32_at = |at: usize| u32::from_le_bytes(bytes[at..at + 4].try_into().unwrap());
        let u64_at = |at: usize| u64::from_le_bytes(bytes[at..at + 8].try_into().unwrap());
        // The local headers: version 2.0, then 4.5 for the one past the point.
        assert_eq!((u16_at(4), u16_at(31 + 4)), (VERSION, VERSION_ZIP64));
        // The central headers: the first with its offset, the second with 0xFFFFFFFF and the
        // ZIP64 extra field (id 1, 8 bytes, the offset).
        let first = 31 + 34;
        assert_eq!(u32_at(first), CENTRAL_HEADER);
        assert_eq!((u16_at(first + 30), u32_at(first + 42)), (0, (base) as u32));
        let second = first + 46 + 1;
        assert_eq!(u32_at(second), CENTRAL_HEADER);
        assert_eq!(
            (u16_at(second + 6), u16_at(second + 30)),
            (VERSION_ZIP64, 12)
        );
        assert_eq!(u32_at(second + 42), u32::MAX);
        let extra = second + 46 + 1;
        assert_eq!((u16_at(extra), u16_at(extra + 2)), (1, 8));
        assert_eq!(u64_at(extra + 4), base + 31);
        // The ZIP64 end record and its locator, then the end record with 0xFFFFFFFF offsets.
        let record = extra + 12;
        assert_eq!(u32_at(record), ZIP64_END_OF_CENTRAL_DIRECTORY);
        assert_eq!((u64_at(record + 24), u64_at(record + 48)), (2, base + 65));
        let locator = record + 56;
        assert_eq!(u32_at(locator), ZIP64_END_LOCATOR);
        assert_eq!(u64_at(locator + 8), base + record as u64);
        let end = locator + 20;
        assert_eq!(u32_at(end), END_OF_CENTRAL_DIRECTORY);
        assert_eq!((u16_at(end + 10), u32_at(end + 16)), (u16::MAX, u32::MAX));
        assert_eq!(bytes.len(), end + 22);
    }

    // 65,535 entries or more need the ZIP64 end records: the classic count field holds 65,534 at
    // most here (65,535 marks the ZIP64 count).
    #[test]
    fn many_entries_take_the_zip64_end_records() {
        for count in [65_534u32, 65_535, 70_000] {
            let mut zip = ZipStream::new(Vec::new());
            for n in 0..count {
                zip.add(format!("{n}").as_bytes(), &n.to_le_bytes())
                    .unwrap();
            }
            let archive = zip.finish().unwrap();
            let files = read(&archive);
            assert_eq!(files.len(), count as usize);
            assert_eq!(files[count as usize - 1].1, (count - 1).to_le_bytes());
            let zip64 = archive
                .windows(4)
                .rev()
                .take(100)
                .any(|window| window == ZIP64_END_LOCATOR.to_le_bytes());
            assert_eq!(zip64, count >= 65_535, "{count}");
        }
    }
}
