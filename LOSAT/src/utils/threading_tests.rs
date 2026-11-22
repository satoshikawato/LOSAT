use super::*;

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3152-3187
#[test]
fn invalid_request_never_enters_search() {
    assert!(with_search_pool(0, "test", |_| -> Result<()> {
        panic!("invalid request entered search")
    })
    .is_err());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:145-188
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
#[test]
fn repeated_pools_include_caller_and_preserve_errors() {
    let caller = std::thread::current().id();
    for n in [1, 2, 4, 2, 1] {
        let error = with_search_pool(n, "test", |search| -> Result<()> {
            assert_eq!(search.threads(), n);
            assert_eq!(std::thread::current().id(), caller);
            if let Some(pool) = search.pool {
                let workers = pool.broadcast(|_| std::thread::current().id());
                assert_eq!(workers.len(), n);
                assert_eq!(workers.iter().filter(|id| **id == caller).count(), 1);
            }
            anyhow::bail!("original search error")
        })
        .unwrap_err();
        assert_eq!(format!("{error:#}"), "original search error");
        assert_eq!(rayon::current_thread_index(), None);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:172-180
#[test]
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
fn caller_registration_is_cleared_after_panic() {
    let failed = std::panic::catch_unwind(|| {
        let _ = with_search_pool(2, "test", |_| -> Result<()> { panic!("search panic") });
    });
    assert!(failed.is_err());
    assert_eq!(rayon::current_thread_index(), None);
    with_search_pool(4, "test", |_| Ok(())).unwrap();
    assert_eq!(rayon::current_thread_index(), None);
}

// The searches of a run share its pool: the same worker threads, no pool of their own; a
// search of another size cannot start inside it, and after the run a search builds its pool.
#[test]
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
fn searches_of_a_run_share_its_pool() {
    let workers = |search: &SearchPool<'_>| {
        let mut ids = search
            .pool
            .expect("a pool")
            .broadcast(|_| format!("{:?}", std::thread::current().id()));
        ids.sort();
        ids
    };
    let (first, second) = with_run_pool(3, "test", || {
        let first = with_search_pool(3, "test", |search| Ok(workers(search))).unwrap();
        let second = with_search_pool(3, "test", |search| Ok(workers(search))).unwrap();
        assert!(with_search_pool(2, "test", |_| Ok(())).is_err());
        assert!(!with_search_pool(1, "test", |search| Ok(search.enabled())).unwrap());
        (first, second)
    });
    assert_eq!(first.len(), 3);
    assert_eq!(first, second);
    // A search after the run builds its pool: new threads besides the caller.
    let caller = format!("{:?}", std::thread::current().id());
    let after = with_search_pool(3, "test", |search| Ok(workers(search))).unwrap();
    assert!(after
        .iter()
        .filter(|id| **id != caller)
        .all(|id| !first.contains(id)));
    // A run of one thread has no pool, and its searches none either.
    assert!(
        !with_run_pool(1, "test", || with_search_pool(1, "test", |search| Ok(
            search.enabled()
        )))
        .unwrap()
    );
    // A run whose pool cannot be built runs without it: its search fails as on its own.
    let too_many = rayon::max_num_threads() + 1;
    let mut ran = false;
    let error = with_run_pool(too_many, "test", || {
        ran = true;
        with_search_pool(too_many, "test", |_| Ok(()))
    })
    .unwrap_err();
    assert!(ran);
    assert_eq!(
        format!("{error:#}"),
        format!(
            "{:#}",
            with_search_pool(too_many, "test", |_| Ok(())).unwrap_err()
        )
    );
}
