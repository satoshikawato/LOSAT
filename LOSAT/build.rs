use std::{env, path::PathBuf, process::Command};

// The argument tables of `ncbi_app::args`, from the embedded `-xmlhelp` XML (line RG DECISIONS
// K4i-1).
#[path = "build/ncbi_args.rs"]
mod ncbi_args;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:145-188
// WASI implementation boundary: Rust 1.92's cdylib link omits startup objects.
// Its bundled crt1-reactor.o exports _initialize, which calls __wasi_init_tp
// and __wasm_call_ctors. Linking it as the entry also prevents LLD from wrapping
// each API/thread export with process constructors/destructors. This is Rust's
// runtime startup, never an NCBI runtime or build dependency.
fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=build/ncbi_args.rs");
    let manifest_dir =
        PathBuf::from(env::var_os("CARGO_MANIFEST_DIR").expect("Cargo CARGO_MANIFEST_DIR"));
    let out_dir = PathBuf::from(env::var_os("OUT_DIR").expect("Cargo OUT_DIR"));
    ncbi_args::write_tables(&manifest_dir, &out_dir);
    println!("cargo:rustc-check-cfg=cfg(losat_wasi_threads)");
    let target = env::var("TARGET").expect("Cargo TARGET");
    if target == "wasm32-wasip1-threads" {
        println!("cargo:rustc-cfg=losat_wasi_threads");
    }
    // LOSAT Web's ABI v1 (`web_api`, `web_v1_argv`, and what only they call): every wasm32 build
    // but one with the feature `no-web-abi-v1` (Cargo.toml), which the WASI command line selects so
    // that it links neither v1's exports nor clap (line RG item K4i-3d). One cfg for every such item.
    println!("cargo:rustc-check-cfg=cfg(losat_web_abi_v1)");
    if env::var("CARGO_CFG_TARGET_ARCH").as_deref() == Ok("wasm32")
        && env::var_os("CARGO_FEATURE_NO_WEB_ABI_V1").is_none()
    {
        println!("cargo:rustc-cfg=losat_web_abi_v1");
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
