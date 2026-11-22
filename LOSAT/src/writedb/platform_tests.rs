use super::*;

// What NCBI's source picks for each platform (writedb_lmdb.hpp:51-57, mdb.c:3979-4003,
// writedb_files.cpp:346-360), whatever platform the test runs on.
#[test]
fn ncbi_s_platform_values_by_operating_system() {
    let windows = ncbi_platform(NcbiOs::Windows, PageSize::Host, 4, u64::MAX);
    assert_eq!(
        (
            windows.default_lmdb_map_size,
            windows.default_taxid_map_size,
            windows.file_length_is_map_size,
            windows.sizeof_long
        ),
        (500_000, 500_000, true, 4)
    );
    let unix = ncbi_platform(NcbiOs::Unix, PageSize::Fixed(16384), 8, u64::MAX);
    assert_eq!(
        (
            unix.default_lmdb_map_size,
            unix.default_taxid_map_size,
            unix.file_length_is_map_size,
            unix.sizeof_long,
            unix.page_size
        ),
        (700_000_000, 300_000_000, false, 8, PageSize::Fixed(16384))
    );
    assert_eq!(
        Platform::LINUX_X86_64,
        ncbi_platform(NcbiOs::Unix, PageSize::Fixed(4096), 8, u64::MAX)
    );
    assert_eq!(
        Platform::LINUX_X86_64.lmdb().unwrap(),
        LmdbPlatform::LINUX_X86_64
    );
}

// The row `Platform::host` takes on the platform the test runs on.
#[test]
fn the_host_is_ncbi_s_platform_of_this_target() {
    let host = Platform::host();
    #[cfg(windows)]
    assert_eq!(
        host,
        ncbi_platform(NcbiOs::Windows, PageSize::Host, 4, usize::MAX as u64)
    );
    #[cfg(all(unix, target_pointer_width = "64"))]
    assert_eq!(
        host,
        ncbi_platform(NcbiOs::Unix, PageSize::Host, 8, u64::MAX)
    );
    #[cfg(target_family = "wasm")]
    assert_eq!(host, Platform::LINUX_X86_64);
    let page = host.lmdb().unwrap().os_page_size;
    assert!(page.is_power_of_two() && page >= 4096, "{page}");
    // Linux x86-64 and Windows x64 have pages of 4 KiB; macOS on Apple silicon 16 KiB.
    #[cfg(any(
        all(target_os = "linux", target_arch = "x86_64"),
        all(windows, target_arch = "x86_64")
    ))]
    assert_eq!(page, 4096);
    #[cfg(all(target_os = "macos", target_arch = "aarch64"))]
    assert_eq!(page, 16384);
}

// The kernel's alignment is the `AT_PAGESZ` that `sysconf` returns (read here by rustix
// from /proc, which the test machines have).
#[cfg(any(target_os = "linux", target_os = "android"))]
#[test]
fn the_kernel_s_page_size_is_at_pagesz() {
    assert_eq!(kernel_page_size().unwrap(), rustix::param::page_size());
    assert_eq!(host_page_size().unwrap(), rustix::param::page_size());
}
