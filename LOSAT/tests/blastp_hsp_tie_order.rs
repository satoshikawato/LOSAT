//! Order of two HSPs of one query-subject pair whose scores tie after the composition-based
//! redo (BLASTP, `-comp_based_stats 2`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1574-1577
//! ```c++
//!         // Sort HSPs with e-values as first priority and scores as
//!         // tie-breakers, since that is the order we want to see them in
//!         // in Seq-aligns.
//!         Blast_HSPListSortByEvalue(hsp_list);
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1392-1395
//! ```c
//!     const double epsilon = 1.0e-180;
//!     if (evalue1 < epsilon && evalue2 < epsilon) {
//!         return 0;
//!     }
//! ```
//! The redo sorts the HSPs by their scaled scores (48512 and 48503 here) and rounds them to the
//! reported score (1516 for both) only afterwards. Both E-values are below 1.0e-180, so the
//! Seq-align conversion ties them and orders them by ScoreCompareHSPs: the smaller subject
//! offset first. The pair is protein 6_2347 against 6_1211 of the OrthoFinder-style benchmark
//! of 2026-10-13, one of its two such pairs; shorter pieces of the two proteins change the
//! composition-adjusted scores and no longer tie. The expected rows are NCBI BLAST+ 2.17.0's
//! (`blastp -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore
//! score" -evalue 0.001 -query q.fa -subject s.fa`).

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

const OUTFMT: &str = "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore score";

const NCBI_ROWS: &str = concat!(
    "6_2347\t6_1211\t60.773\t543\t1\t514\t1\t534\t0.0\t588\t1516\n",
    "6_2347\t6_1211\t57.741\t549\t604\t1120\t700\t1233\t0.0\t588\t1516\n",
);

const QUERY: &str = concat!(
    ">6_2347\n",
    "MANDSPAKSLVDIDLASLRDPAGIFELVEVVGNGTYGQVYKVGSFFFYDEEEEIKLEINM\n",
    "LKKYSHHRNIATYYGAFIKKSPPGHDDQLWLVMEFCGAGSITDLVKNTKGNSLKEDWIAY\n",
    "ISREILRGLAHLHAHHVIHRDIKGQNVLLTENAEVKLVDFGVSAQLDRTVGRRNTFIGTP\n",
    "YWMAPEVIACDENPDATYDYRSDLWSCGITAIEMAEGAPPLCDMHPMRALFLIPRNPPPR\n",
    "LKSKKWSKKFFSFIEGCLVKNYTQRPPTEQLLKHPFIRDQPNERQVRIQLKDHIDRTKKK\n",
    "RGEKDETEYEYSGSEEEEEESSEQEGEPSSIVNVPGESTLRRDFIRLQQENKERSEALRR\n",
    "QQLLQEQQLREQEEYKRQLLAERQKRIEQQKEQRRRLEEQQRREREMRRQQEREQRRREQ\n",
    "EEKRRMEEIERRRKEEEERRRAEEEKRRADREQEYIRRQLEEEQRHLEILQQQLLHEQAM\n",
    "LLEYKWRELEEQRQAQKLQRQLQLEQAYLLSLQQTQNQDNSKTAQQQSTKQPLSIEPDRE\n",
    "NPGKTPNLDVPVRTTSRSPVLSRRDSPHHHGNASHRDAVNRNTVSSAEPRLLWERVEKLV\n",
    "PRTTGNSGSGGSSSSSSSSNSSSQPSSHCGSGERFRARCESTRPFSKSSLKTHTRTGELR\n",
    "AVDDVRPPNKVTEYSSSSEDSDTTDEDDDEEVDQEAGEESTSGPEDSKAVSSKLSNGETE\n",
    "SVKTMVVHEGESDAGSTPCKDSTLIVRQSQSSNNMQKHKSSSSFTPFIDPRLLQVSPSTG\n",
    "SALNNVGYSNDARLAEALRADPSRKGSVVNVNPVNTRPQSDTPEIRKYKKRFNSEILCAA\n",
    "LWGVNLLVGTESGLMLLDRSGQGKVYPLINRRRFQQMDVLEGLNVLVTISGKKNKLRVYY\n",
    "LSWLRNKILHNDPEVEKKQGWTTVGDLEGCVHYKVVKYERIKFLVLALKNSVEVYAWAPK\n",
    "PYHKFMAFKSFGDLVHKPLLVDLTVEEGQRLKVIYGSYSGFHAVDVDSGAVYDIYLPTHI\n",
    "QTSIQSHAIIILPNTDGIELLVCYEDEGVYVNTYGRITKDVVLQWGEMPTSVAYIRSNQI\n",
    "MGWGEKAIEIRSVETGHLDGVFMHKRAQRLKFLCERNDKVFFASVRSGGSSQVYFMTLGR\n",
    "SNLLSW\n",
);

const SUBJECT: &str = concat!(
    ">6_1211\n",
    "MSRENTSRSLDSIDLAALRDPAGIFELVEVVGNGTYGQVYKGRHVKTGQLAAIKVMEVTE\n",
    "EEEEEIKLEINMLKSYSHHRNIATYYGAFVKKSPAGQDHQLWLVMEYCGAGSVTDLVKKT\n",
    "KGNCFKEDWIAYICREVLRGLSHLHSHHVIHRDIKGQNVLLTENAEVKLVDFGVSAQLDK\n",
    "TIGRRNTFIGTPYWMAPEVIACDENPEATYDYRSDLWSLGITALEMAEGAPPLCDMHPMR\n",
    "ALFLIPRNPPPKLKSKKWSKRFLSFVDSCLVKNHLHRPATDALLRHAFIRDLPNERHVRI\n",
    "MLKDHLDRTRKKRDKEGLEYEYSGSEDDDDEVTEEEGEPSSIVNVPGEWTLRREFLRLQT\n",
    "DERSKEGGHRGGQQRQQQHRQQLAEQERYKQQFLAERQRRIRQQQEQRQKLEDQQQRLQL\n",
    "TDGAFQWAELQEVFVGEEPDLHNSRIPLDEKNRRCDRRQDQTGRQRSADGSVRPLDSTAD\n",
    "QNAIAPNQLLQLKTQHFQPQAKLLPPPQHQQPSAGRRSHRTLPKDQTQLLSLQHDHRHRD\n",
    "KERQKPREPYMAAVHKLTANSGNKAAAAAAAAAAAPSPSRPANRQHTAMELKESSVASIL\n",
    "LAANLPSQVHSRLGSGSSSSPCTPAMQRAHFNQNANLRSSRDSPHSSSMSDMVISPLSPF\n",
    "APLSPTDNVFWDPAEEPPKVPERTTSKFYCKELVMAHKRTASSGSPVSAADKSRGSASQG\n",
    "GPPLSSNHQGVLKLESPLPQHRAQRKAENAASRLAGFKPVLGRLNKAAEYSSSSDDVSSD\n",
    "DEERIRSGLSNGKGKYFRGISDGIDLQLHRNLNQVSRLGSDDTYLLPDLLQKTSSSSPTR\n",
    "DALSQQRTLDFSYPMSPTAPGSPIFQYLHDVSSPSGKSTSSSSSFLSFIDPRLIPTSSPP\n",
    "HSPTALKSEPIKRENHRKGSVVNVNPTNIRPQSDTPEIRKYKKKFNAEVLCAGLWGVNLL\n",
    "VGTENGLWLLDRSGQGKVYSLISRRRFQQMDVLEGLNVLITISGKKNKLRLYYLSWLRNK\n",
    "ILHNDPEVEKRQGWTSVGDLEGCVHYKVVRYEKIKFLVVALKNAVEVYAWAPKPYHKFMA\n",
    "FKSFRSLPQKPLFVDLTVEEGQRLKVIYGSSSGFHAIDVDSGTPYDLYLPTHIQGLIHPH\n",
    "AIIILPNKAGTEVLVCYEDEGVYIDTYGRITKDTVLQWGEMPASVAYLQSNQVMGWGEKA\n",
    "IELRSANTGNLEGVFMHKKAQKLKFLCERNDKVGSHSKLEPRPPFDSKTVIVSDISTGNP\n",
    "RCADQTKSKVVCSKT\n",
);

/// The query and subject files, written once.
fn inputs() -> &'static (PathBuf, PathBuf) {
    static FILES: OnceLock<(PathBuf, PathBuf)> = OnceLock::new();
    FILES.get_or_init(|| {
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("blastp_hsp_tie_order");
        std::fs::create_dir_all(&dir).expect("work directory");
        let files = (dir.join("query.faa"), dir.join("subject.faa"));
        std::fs::write(&files.0, QUERY).expect("query file");
        std::fs::write(&files.1, SUBJECT).expect("subject file");
        files
    })
}

/// One BLASTP run of the inputs with no `LOSAT_*` variables except `LOSAT_BLASTX_WORKER_LOG`:
/// a build with the worker probe (`--all-features`) writes to that file as soon as a pool
/// starts, so the caller's value is kept, or a file of the tests is named.
fn blastp(threads: usize) -> String {
    const WORKER_LOG: &str = "LOSAT_BLASTX_WORKER_LOG";
    let (query, subject) = inputs();
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
    for (name, _) in std::env::vars_os() {
        let name_text = name.to_string_lossy();
        if name_text.starts_with("LOSAT_") && name_text != WORKER_LOG {
            command.env_remove(name);
        }
    }
    if std::env::var_os(WORKER_LOG).is_none() {
        command.env(WORKER_LOG, query.with_file_name("worker.log"));
    }
    let output = command
        .arg("blastp")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subject)
        .args(["-outfmt", OUTFMT, "-evalue", "0.001"])
        .args(["-num_threads", &threads.to_string()])
        .output()
        .expect("run LOSAT");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).expect("UTF-8 report")
}

#[test]
fn hsps_with_rounded_score_ties_below_the_evalue_epsilon_are_in_ncbi_order() {
    for threads in [1, 4] {
        assert_eq!(blastp(threads), NCBI_ROWS, "{threads} threads");
    }
}
