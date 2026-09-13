use std::{env, path::PathBuf, process::Command};

// NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:145-188
// TBlastThreads the_threads(GetNumberOfThreads());
// (*thread)->Run(); (*thread)->Join(&result);
// WASI implementation boundary: Rust 1.92's cdylib link omits startup objects.
// Its bundled crt1-reactor.o exports _initialize, which calls __wasi_init_tp
// and __wasm_call_ctors. Linking it as the entry also prevents LLD from wrapping
// each API/thread export with process constructors/destructors. This is Rust's
// runtime startup, never an NCBI runtime or build dependency.
fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rustc-check-cfg=cfg(losat_wasi_threads)");
    let target = env::var("TARGET").expect("Cargo TARGET");
    if target == "wasm32-wasip1-threads" {
        println!("cargo:rustc-cfg=losat_wasi_threads");
    }
    if target != "wasm32-wasip1" && target != "wasm32-wasip1-threads" {
        return;
    }
    let output = Command::new(env::var_os("RUSTC").expect("Cargo RUSTC"))
        .args(["--print", "sysroot"])
        .output()
        .expect("query Rust sysroot for the WASI reactor CRT");
    assert!(output.status.success(), "rustc --print sysroot failed");
    let sysroot = String::from_utf8(output.stdout).expect("UTF-8 Rust sysroot");
    let crt = PathBuf::from(sysroot.trim())
        .join("lib/rustlib")
        .join(target)
        .join("lib/self-contained/crt1-reactor.o");
    assert!(crt.is_file(), "WASI reactor CRT missing: {}", crt.display());
    println!("cargo:rustc-cdylib-link-arg={}", crt.display());
    println!("cargo:rustc-cdylib-link-arg=--entry=_initialize");
}
