# NCBI reference: c++/src/app/blast/blastn_app.cpp:59-67
# m_StopWatch.Start(); m_StopWatch.Elapsed();
# Measurement-only matched Rust dispatch clocks; no inner-loop instrumentation.
from pathlib import Path
import argparse, shutil, subprocess

p = argparse.ArgumentParser()
p.add_argument('versions', nargs='+')
a = p.parse_args()
e = Path(__file__).resolve().parent
for version in a.versions:
    name = version + 'body'
    dest = e / 'work' / (name + '-source')
    shutil.copytree(e / 'work' / (version + '-source'), dest)
    source = dest / 'LOSAT/src/main.rs'
    raw = source.read_bytes()
    needle = b'    match cli.command {'
    assert raw.count(needle) == 1
    raw = raw.replace(needle, b'    // NCBI blastn_app.cpp:59-67: m_StopWatch.Start(); m_StopWatch.Elapsed();\n    let body_start = std::time::Instant::now();\n' + needle)
    needle = b'    Ok(())'
    assert raw.count(needle) == 1
    raw = raw.replace(needle, b'    eprintln!("[BODY_SCOPE_SECONDS] {:.9}", body_start.elapsed().as_secs_f64());\n' + needle)
    source.write_bytes(raw)
    subprocess.run(['python3', str(e / 'build.py'), name, '--kinds', 'native', 'threaded'], check=True)
