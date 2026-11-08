//! The command `losat`. Its body is the library's `main_entry` (`src/cli/main_entry.rs`), which
//! cites NCBI's application entries (`blastn_app.cpp`'s `NcbiSys_main` -> `AppMain`).

fn main() -> anyhow::Result<()> {
    LOSAT::main_entry()
}
