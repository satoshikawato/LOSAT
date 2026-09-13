//! Search-scoped execution; NCBI owns search semantics, Rayon only schedules.
use anyhow::{ensure, Context, Result};

// NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:145-188
// TBlastThreads the_threads(GetNumberOfThreads());
// (*thread)->Run(); (*thread)->Join(&result);
// This borrow cannot outlive the one search that owns the pool.
pub struct SearchPool<'a> {
    requested: &'a usize,
    #[cfg(feature = "parallel")]
    pool: Option<&'a rayon::ThreadPool>,
}

impl SearchPool<'_> {
    // NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:147
    // TBlastThreads the_threads(GetNumberOfThreads());
    pub fn threads(&self) -> usize {
        *self.requested
    }

    // NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:177-188
    // (*thread)->Run(); (*thread)->Join(&result);
    pub fn enabled(&self) -> bool {
        #[cfg(feature = "parallel")]
        return self.pool.is_some();
        #[cfg(not(feature = "parallel"))]
        false
    }

    // NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:177-188
    // (*thread)->Run(); (*thread)->Join(&result);
    // Nested DP/linking/redo work uses this same pool, never the global pool.
    pub fn install<F, R>(&self, work: F) -> R
    where
        F: FnOnce() -> R + Send,
        R: Send,
    {
        #[cfg(feature = "parallel")]
        if let Some(pool) = self.pool {
            return pool.install(work);
        }
        work()
    }
}

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3152-3187,3205-3236
// arg_desc.SetConstraint(kArgNumThreads, new CArgAllowValuesGreaterThanOrEqual(1));
// NCBI warns when restricting threads. LOSAT's explicit-N contract instead
// rejects unsupported requests, without silently reducing the requested pool.
pub fn validate_threads(requested: usize) -> Result<()> {
    ensure!(requested > 0, "num_threads must be greater than zero");
    if let Some(raw) = std::env::var_os("LOSAT_WASI_THREAD_CAP") {
        let cap = raw
            .to_str()
            .context("LOSAT_WASI_THREAD_CAP must be a positive integer")?
            .parse::<usize>()
            .context("LOSAT_WASI_THREAD_CAP must be a positive integer")?;
        ensure!(cap > 0, "LOSAT_WASI_THREAD_CAP must be a positive integer");
        ensure!(
            requested <= cap,
            "requested {requested} threads exceeds LOSAT_WASI_THREAD_CAP={cap}"
        );
    }
    if requested == 1 {
        return Ok(());
    }
    ensure!(
        cfg!(all(
            feature = "parallel",
            any(
                not(target_arch = "wasm32"),
                all(losat_wasi_threads, feature = "wasm-threads")
            )
        )),
        "unsupported num_threads={requested}: this build does not support parallel search"
    );
    #[cfg(feature = "parallel")]
    ensure!(
        requested <= rayon::max_num_threads(),
        "requested {requested} threads exceeds Rayon maximum {}",
        rayon::max_num_threads()
    );
    Ok(())
}

// NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:145-188
// TBlastThreads the_threads(GetNumberOfThreads());
// (*thread)->Run(); (*thread)->Join(&result);
// Rayon-core 1.13.0 src/lib.rs:322-343 uses spawn_scoped and std::thread::scope;
// the pool drops before scope joins all workers, including partial spawn errors.
// Do not use use_current_thread: its caller registry survives pool destruction.
pub fn with_search_pool<F, R>(requested: usize, program: &str, work: F) -> Result<R>
where
    F: FnOnce(&SearchPool<'_>) -> Result<R> + Send,
    R: Send,
{
    validate_threads(requested)?;
    #[cfg(feature = "parallel")]
    if requested > 1 {
        return rayon::ThreadPoolBuilder::new()
            .num_threads(requested)
            .build_scoped(
                |thread| thread.run(),
                |pool| {
                    ensure!(
                        pool.current_num_threads() == requested,
                        "requested {requested} threads but constructed pool has {}",
                        pool.current_num_threads()
                    );
                    report_pool(program, requested, pool.current_num_threads());
                    let search = SearchPool {
                        requested: &requested,
                        pool: Some(pool),
                    };
                    search.install(|| work(&search))
                },
            )
            .with_context(|| format!("failed to build {program} pool with {requested} threads"))?;
    }
    report_pool(program, requested, 0);
    work(&SearchPool {
        requested: &requested,
        #[cfg(feature = "parallel")]
        pool: None,
    })
}

// NCBI reference: c++/src/algo/blast/core/blast_engine.c:1659-1668
// /* Use a local diagnostics structure ... shared between multiple threads */
pub fn diagnostics_enabled() -> bool {
    std::env::var("LOSAT_WASI_THREADS_DEBUG")
        .is_ok_and(|v| v == "1" || v.eq_ignore_ascii_case("true"))
}

// NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:147,177-188
// TBlastThreads the_threads(GetNumberOfThreads()); (*thread)->Join(&result);
fn report_pool(program: &str, requested: usize, pool: usize) {
    if diagnostics_enabled() {
        eprintln!("[losat-thread-pool] program={program} requested_threads={requested} pool_threads={pool} caller_participates=false measured_activity=null");
    }
}

// NCBI reference: c++/src/algo/blast/core/blast_engine.c:1659-1668
// /* Use a local diagnostics structure ... shared between multiple threads */
// Record stage selection separately from the already verified pool size.
pub fn report_stage(program: &str, stage: &str, work_items: usize, parallel: bool) {
    if diagnostics_enabled() {
        eprintln!("[losat-thread-stage] program={program} stage={stage} work_items={work_items} parallel_selected={parallel} measured_activity=null");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3152-3187
    // arg_desc.SetConstraint(kArgNumThreads, new CArgAllowValuesGreaterThanOrEqual(1));
    #[test]
    fn invalid_request_never_enters_search() {
        assert!(with_search_pool(0, "test", |_| -> Result<()> {
            panic!("invalid request entered search")
        })
        .is_err());
    }

    // NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:145-188
    // TBlastThreads the_threads(GetNumberOfThreads()); (*thread)->Run(); (*thread)->Join(&result);
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    #[test]
    fn repeated_pools_keep_nested_work_off_caller_and_preserve_errors() {
        let caller = std::thread::current().id();
        for n in [1, 2, 4, 2, 1] {
            let error = with_search_pool(n, "test", |search| -> Result<()> {
                assert_eq!(search.threads(), n);
                assert_eq!(std::thread::current().id() == caller, n == 1);
                if let Some(pool) = search.pool {
                    let workers = pool.broadcast(|_| std::thread::current().id());
                    assert_eq!(workers.len(), n);
                    assert!(workers.iter().all(|id| *id != caller));
                }
                anyhow::bail!("original search error")
            })
            .unwrap_err();
            assert_eq!(format!("{error:#}"), "original search error");
        }
    }
}
