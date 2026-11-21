use super::*;
use std::sync::Mutex;

#[derive(Default)]
struct Kept(Mutex<Vec<(Vec<u8>, Vec<u8>)>>);

impl FileReceiver for Kept {
    fn receive(&self, name: &[u8], bytes: &[u8]) -> io::Result<()> {
        self.0.lock().unwrap().push((name.to_vec(), bytes.to_vec()));
        Ok(())
    }
}

/// The names of `-logfile` of a receiver: NCBI's special names, a file in memory handed over
/// under its name, `cerr` for a name with a directory, which a run that searches refuses
/// (`log_refusal`; line WQA's D14, 2026-11-21).
#[test]
fn a_receiver_keeps_the_log_in_memory_and_refuses_directories() {
    let kept = Arc::new(Kept::default());
    let files = RunFiles::Receiver(kept.clone());
    let target = |name: Option<&str>| files.log_target(name.map(OsStr::new));
    let refusal = |name: Option<&str>| files.log_refusal(name.map(OsStr::new));
    assert!(matches!(target(None), DiagTarget::Stderr));
    assert!(matches!(target(Some("-")), DiagTarget::Stderr));
    assert!(matches!(target(Some("")), DiagTarget::Discard));
    assert!(matches!(target(Some("/dev/null")), DiagTarget::Discard));
    assert!(matches!(target(Some(".")), DiagTarget::Stderr));
    assert!(matches!(target(Some("..")), DiagTarget::Stderr));
    for name in [
        None,
        Some("-"),
        Some(""),
        Some("/dev/null"),
        Some("."),
        Some(".."),
        Some("run.log"),
    ] {
        assert_eq!(refusal(name), None, "{name:?}");
    }
    for name in ["d/run.log", "./run.log", "/tmp/run.log"] {
        assert!(matches!(target(Some(name)), DiagTarget::Stderr));
        assert_eq!(
            refusal(Some(name)).as_deref(),
            Some(
                format!(
                    "the -logfile name {name} has a directory, which LOSAT Web does not have: give a file name without one"
                )
                .as_str()
            )
        );
        // The command line writes the file (or to `cerr`) and refuses nothing.
        assert_eq!(RunFiles::Disk.log_refusal(Some(OsStr::new(name))), None);
    }
    let log = target(Some("run.log"));
    log.write_text(b"posted\n", &mut io::sink());
    files.deliver_log(Some(OsStr::new("run.log")), log).unwrap();
    // A target without a file hands nothing over.
    files
        .deliver_log(Some(OsStr::new("-")), DiagTarget::Stderr)
        .unwrap();
    assert_eq!(
        *kept.0.lock().unwrap(),
        vec![(b"run.log".to_vec(), b"posted\n".to_vec())]
    );
}
