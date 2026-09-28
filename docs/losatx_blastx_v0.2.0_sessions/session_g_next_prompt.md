# LOSATX（BLASTX）v0.2.0 — Session G 総合認証・性能・独立監査

これはソフトウェアエンジニアリングです。配列は計算プログラムの固定テスト入力として扱い、生物学的解釈を行わない。

受入済み Session F とその後の実装コミットを起点に、Session G のローカル実装・検証・正式性能測定・認証資料と release handoff の準備を完遂してください。G はまだ開始されていない。タグ、release/package 公開、deployment、新たな commit/push は今回許可しない。未実施の必須 gate を合格扱いせず、必要な作業を継続する。製品・governance 承認や実行環境が本当に不足する境界では、独立して進められる作業を完了してから具体的な残件を示す。

## 固定作業場所と承継 identity

- worktree: `/tmp/losatx-blastx-v020`
- branch: `feature/losatx-blastx-v0.2.0`
- base: `2976bd5f427cc5448315b6605ab85ee0787545a5`
- F 受入時 HEAD: `eb7cb676f4b401169928a7e85dde42d9739738f9`。その後の実装コミットが現在の起点なので、古い HEAD へ戻さない。
- コミット対応表: `docs/losatx_blastx_session_f_checkpoint.json`
- F evidence: `docs/evidence/losatx_stage_f_parallel_wasi_web/run-20260927T001743Z/`
- frozen candidate SHA-256: `2757f2fb21dc7902d28e3ee9d1415a8d313b223ed1bb15fa23d4cb3e191a9a3f`
- frozen contract SHA-256: `be2d3c85800bdb077d3d46ebf03cdf9a5ab18ab762777b790aff68d397477187`
- frozen strict validator SHA-256: `f1f0934212d060af2b78d3791251b967a0bbab4abd3db7197b3e98f0ed7b050f`
- frozen MANIFEST SHA-256: `d2a836a0ff17318ea70873a5a1b09fc74b94e7638dc62004faf165795bf9ede5`
- frozen source archive SHA-256: `fdfeed2fbaf1333f9d58602e632414eb19ab1a4f5a19d916047359a05b6e936b`

作業ツリー・index・remote の実状態を読み取りで確認し、現在の commit SHA と source closure を固定する。別 worktree 作成、branch 切替、reset/checkout/stash、既存変更の巻戻し、過去の成果物の削除はしない。

## 開始前に読むものと gate

1. `AGENTS.md`、`verify-ncbi-parity-and-speed` と必須参照、総合計画 `docs/losatx_blastx_v0.2.0_plan.md` 第 7〜14 節、`docs/losatx_blastx_v0.2.0_sessions/session_g_certification.md`。
2. checkpoint JSON/説明、F の `SESSION_G_HANDOFF.md`、`GATE.md`、`REPRODUCE.md`、`completion_state.json`、`coverage.json`、`selection_rules_final.json`、`source_call_state.md`、`first_differences.md`、`independent_audit.md` と A〜E の必要な authority/coverage/raw 証拠。
3. F の strict acceptance `commands/strict_acceptance_v1` と最終 integrity `manifest_receipts/final` の実 receipt、stdout/stderr/status、auditor binding、22 refusal 記録を検証する。artifact integrity は実装受入を決定しない。
4. コミット後は凍結 F validator の旧 HEAD/index 条件が現在と一致しない。旧 validator/contract/candidate を書き換えず、拒否を PASS と読み替えず、Git を巻き戻さない。新コミットの全 485 対象 Git blob と working file が checkpoint/F closure に一致し、F raw 証拠・archive・authority・監査・refusal が有効であることを検証する、現在の HEAD/index に別途結び付けた G 承継 gate を作る。F の凍結 manifest validator はそのまま実行して整合性を確認する。
5. 開始時 inventory/index/closure/archive、保全対象と編集 owner を固定する。G の証拠・target/output は新規 namespace に保存し、F evidence を上書きしない。G 完了状態や handoff を gate 前に生成しない。

## 維持する authority と実装契約

NCBI C/C++ の唯一の authority は commit `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`。official BLASTX は `/home/kawato/micromamba/bin/blastx`、2.17.0+、SHA-256 `5c29dd472f7b18db37b688077ea4a48eea8f46ccfbeab64e11ca5efddc142c05` と A の runtime-library registry。未知 fingerprint を拒否する。

NCBI は source inspection と比較 oracle に限る。runtime/build/FFI/subprocess/fallback 依存にしない。Rust の変更は NCBI snippet・file・line を直前コメントに付け、candidate、HSP、cutoff、collector/heap、Kappa/link/filter、統計、丸めと report の timing/order/context を保つ。

BLASTX query code は `1,2,3,4,5,6,9,10,11,12,13,14,15,16,21,22,23,24,25,26,27,28,29,30,31,33` の 26 ID。33 を含め、32 を拒否する。protein subject に db_gencode はなく、5,000,000 残基超は現在明示的未対応。TBLASTN/TBLASTX の例外を BLASTX に拡張しない。

同じ Rust BLASTX engine と serial input-OID reduction を使う。plain `wasm32-wasip1` は serial、実 threading は `wasm32-wasip1-threads` + `wasm-threads`。Web は既存 ABI と alloc/store/release/result/error 所有権を維持し、保存した原 FASTA bytes を再解析する。generic FASTA handle に molecule-kind tag はない。CLI default outfmt 0 と Web default 6 を区別し、比較は同じ resolved options と固定 FASTA/label identity で行う。新 SDK、UI、gbdraw 改修は範囲外。

## Session G 必須作業

1. X01〜X18 → NCBI owner → Rust production owner → fixture/分岐到達 → raw evidence の双方向 ledger を検査し、M1〜M11 と既知の相互作用・境界を計画どおり認証する。M1 は 26 codes × formats 0/6/7 × native threads 1/4 = 156 code-format-thread 行、検索 option としては 52 code-thread 組である。再実行や format/入口を独立検索として重複計上しない。全 codon unit matrix は検索 matrix の代替にしない。新しい実行場所で必要な再現を行い、既存証拠は code/input/environment/acceptance が不変と証明できる範囲だけ再利用する。
2. Native threads 1/2/4/8、serial/threaded command WASI、serial/threaded reactor の memory/handles と実 Chromium を同じ resolved options で確認する。実 worker spawn/compute/join、pool resize/failure/回復、stale result/error clearing、反復メモリと worker 回収を証明する。ビルド成功や Node だけを browser/runtime 認証に使わない。
3. 同一候補 source identity で build/test/fmt/clippy、適用可能な artifact smoke、BLASTN/megablast、BLASTP、TBLASTX、TBLASTN の必要回帰を実施する。PR6 Gate A = 凍結 PR5 LOSAT bytes、Gate B = 六 official 検索の登録済み platform fingerprints を分離したまま維持する。
4. 配布 Linux x86_64、Windows x86_64、macOS arm64/x86_64 および計画で要求された追加対象の source/toolchain/flags/binary/oracle identity と native gate を確認する。未実施 platform を認証済みにしない。BLASTX の未説明 platform raw variance は hard-fail。新 bounded authority の受入は特性化・review・明示的製品/governance 承認が必要であり、既存 PD を流用しない。
5. パリティ確立後、第 11 節の B1〜B5 固定 workload（小型、R06 fragment batch、R01/R02/R05 全長 viral、R03/R04 code4、R07 bacterial）を正式測定する。各 case/mode で untimed warmup 1 回 + retained timed 3 回ちょうど。全 raw samples、中央値、三標本 min–max、wall/取得可能な CPU time・peak RSS、output hash、threads、input/source/binary/toolchain/environment identity を保存する。最速値選択や遅い標本の除外をしない。追加測定は明示要求または三標本が実際に判定不能の場合のみ、理由を記録する。
6. NCBI timing は事前 `makeblastdb -dbtype prot` の `-db` 検索で測る。DB 作成 command/version/input hash/経過時間は検索時間と別保存。BLASTX parity/hit-distribution oracle は同じ query_gencode の NCBI local `-subject`。timing DB output と分布 oracle を別名・別 metadata にし、local-subject を threaded NCBI timing と呼ばない。LOSAT local-subject と NCBI DB の違いを説明し、同条件 speedup と短絡しない。
7. 最初の正しい serial BLASTX を baseline とし、別 program を BLASTX baseline にしない。cold command と warm reactor の計測境界（load/compile/worker startup を含む範囲）を固定し、並列化で遅くなる結果も保持する。diagnostic worker-probe の時間を production 性能に使わない。共有部変更の性能退行を適切な同条件三標本で評価し、明確な悪化は原因解消または具体的な製品判断まで未解決とする。長時間 benchmark のプロセス poll は原則 10 分間隔、ユーザーへの説明は別途継続する。
8. README、CHANGELOG、scope/help、manifest、使用例、release handoff を最終候補と実証済み能力に揃える。現行文書に残る Session E の serial-only 説明を修正する。巨大な再生成可能出力を無条件に Git に入れず、fixture/provenance/hash/再現手順と認証記録を恒久的に保存し、外部一時 path だけで release 証拠を完結させない。

raw compare に sort、E-value tolerance、header/footer/空白除去、独自 canonicalization を導入しない。未知 fingerprint、missing、skip、重複、既知差分を PASS にしない。差分は最初の byte/field と最初の NCBI/Rust stage を特定し、同じ owner の全 offender を修正して必要 gate を再実行する。

## 受入 gate と独立監査

G の実装認証 gate は summary を信頼せず、実行 argv/cwd/env、input bytes、binary/target/features、stdout/stderr/exit、ABI/host buffers、worker と failure/recovery/memory 記録、source/archive/index、coverage、authority、保全 baseline、正式 benchmark の warmup/全三標本/集計/output hashes を raw から検証する。artifact integrity gate は独立させる。代表的な欠落・改変・未知 identity・stale audit を実際に拒否する negative mutation を記録し、G contract と validator 自身も固定する。

release-facing parity/performance 受入前に、本物の `ncbi_parity_auditor` custom agent へ独立 read-only 監査を依頼する。この監査の delegation は許可する。監査を自作の PASS で代替しない。candidate/contract/validator/raw/benchmark を同一 binding にし、指摘を解決する。変更後は失効した証拠・監査を使わず、必要な gate と監査を更新する。

X01〜X18 と第 14 節をすべて満たした場合だけ Session G / v0.2.0 gate PASS とする。有限 fixture を全入力の証明と呼ばない。未実施 platform/必須測定/監査が残る場合は正確に未完了を記録する。

## 最終報告と handoff

- 現在の commit、最終 source/candidate/contract/validator/manifest と各 production artifact SHA。
- 全 gate、M1〜M11、platform/target/入口の状態と正確な coverage、再利用証拠の根拠。
- B1〜B5 各 mode の全三標本、中央値・min–max、計測境界と出力一致。
- 初回差分・NCBI owner・修正・例外範囲・未対応・残件、実独立監査と refusal 証拠。
- REPRODUCE、認証記録、release handoff、英語 proposed commit title と短い summary。

最終状態と handoff は actual gate/監査 PASS の後に作り、最終 manifest を再検証する。今回のプロンプトは commit/push/tag/publication を許可しない。
