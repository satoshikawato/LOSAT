# NCBI reference: c++/src/app/blast/blastn_app.cpp:59-67
# m_StopWatch.Start(); m_StopWatch.Elapsed();
# Diagnostic-only coarse scope. No per-HSP/cell counters are enabled.
from pathlib import Path
import argparse,shutil,subprocess
p=argparse.ArgumentParser();p.add_argument('version');a=p.parse_args();b=Path(__file__).resolve().parent;name=a.version+'body';dest=b/(name+'-source');shutil.copytree(b/(a.version+'-source'),dest,dirs_exist_ok=True)
p=dest/'LOSAT/src/main.rs';s=p.read_text();needle='    match cli.command {';assert s.count(needle)==1;s=s.replace(needle,'''    // NCBI reference: c++/src/app/blast/blastn_app.cpp:59-67
    // m_StopWatch.Start(); m_StopWatch.Elapsed();
    // Diagnostic scope: engine call including configuration, input, search and
    // output; CLI parse, host startup/compile and process teardown are excluded.
    let body_start = std::time::Instant::now();
'''+needle);s=s.replace('    Ok(())','    eprintln!("[BODY_SCOPE_SECONDS] {:.9}", body_start.elapsed().as_secs_f64());\n    Ok(())');p.write_text(s)
subprocess.run(['python3',str(b/'build.py'),name,'--kinds','native','threaded'],check=True)
