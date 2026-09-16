# 採用判断と最終結果

**H1aのserial Module再利用を採用。X1は速度基準を満たさず見送り。X2・M1・H1bは条件不成立で見送り。**

Node 26.8.2 / V8 14.6.202.34-node.28、i9-14900HX。主条件は通常のNode設定。
Wasm / Native ≤1.20は到達目標であり、採用の必須条件ではない。現状からの改善と非回帰で採否を判断する。 有効queryの指定fixtureによる限定評価であり、全入力・全platformの互換認証ではない。

## 最終結果

| 対象・区間 | baseline → 最終版 | 短縮 | 減った費用 | Wasm / Native | 範囲・未達事項 |
|---|---:|---:|---|---:|---|
| H1a・serial準備 | 0.007541 → 0.004700 s | 0.002841 s (37.68%) | 検証済みbytesへの二度目のcompile要求 | 対象外 | Node serial command・valid-nohit |
| H1a・起動込み | 0.051205 → 0.050399 s | 0.000806 s (1.57%) | 同上 | 未測定 | 全体5%短縮は主張しない |
| MjeNMV / MelaMJNV・threaded n1本体 | 22.406 s（baseline参照） | 0 s | 本体変更なし | 1.567 | 通常Node、目標≤1.20未達 |
| MjeNMV / MelaMJNV・threaded n8本体 | 35.077 s（baseline参照） | 0 s | 本体変更なし | 4.203 | 通常Node、目標≤1.20未達 |
| AP027280 self・threaded n1本体 | 37.147 s（baseline参照） | 0 s | 本体変更なし | 1.370 | 通常Node、目標≤1.20未達 |
| AP027280 self・threaded n8本体 | 61.850 s（baseline参照） | 0 s | 本体変更なし | 3.412 | 通常Node、目標≤1.20未達 |

H1aは検索kernelを変更しない。threaded本体の表は凍結baselineの参考値で、採用した本体改善量は0としている。H1適用後のthreaded起動込み時間を再計測した値ではない。
検査APIは元と同じJSONを返し、serial実行だけが検証済みModuleを再利用する。
新しいinstance/memoryの所有権、reactor初期化、終了コード、Worker失敗伝播は維持した。
ブラウザ配布artifactは更新しておらず、ブラウザ・Wasmtimeの速度は未測定。

## H1a：全条件と非回帰基準

各版warmup 1回、AB/BA/ABの3組。主区間は測定前に宣言したserial準備区間。
時間は3回の中央値、RSS基準は保守的に3回の最大peak RSSでも確認した。
時間悪化はmax(5%, 50 ms)、RSS増加はmax(10%, 16 MiB)以内。

| Fixture | 準備 baseline → H1 (ms) | 起動込み baseline → H1 (ms) | 時間差 H1−baseline (ms) | 最大RSS baseline → H1 (MiB) | 判定 |
|---|---:|---:|---:|---:|---|
| valid-nohit | 7.541 → 4.700 | 51.205 → 50.399 | -0.806 | 77.12 → 72.09 | PASS |
| single32 | 10.389 → 4.508 | 100.469 → 94.661 | -5.808 | 97.67 → 96.11 | PASS |
| word7 | 10.898 → 4.408 | 110.680 → 125.764 | +15.084 | 101.32 → 96.44 | PASS |
| EDL933.Sakai.losatn.megablast | 10.460 → 4.364 | 1300.815 → 1301.094 | +0.279 | 329.44 → 316.75 | PASS |
| MelaMJNV.PemoMJNVA.losatn.blastn | 10.621 → 4.511 | 1489.423 → 1480.943 | -8.480 | 119.70 → 115.34 | PASS |

BLASTNのMelaMJNV.PemoMJNVAという表示名はLC738874/LC738870を指す。全FASTA・argvはmanifestに記録した。 word7の15.084 ms悪化を除外していない。3組では変動の原因を断定しない。
全40回（warmup込み）でfresh oracleのraw bytesと終了状態が一致した。
二つのcompile API要求があったことから、機械語生成も二重だったとは主張しない。

## X1：独立実装を評価し、見送り

| 主入力 | n | baseline本体 (s) | X1本体 (s) | 悪化 | 判定 |
|---|---:|---:|---:|---:|---|
| MjeNMV / MelaMJNV | 1 | 22.406 | 22.688 | 1.26% | 見送り |
| MjeNMV / MelaMJNV | 8 | 35.077 | 37.823 | 7.83% | 見送り |
| AP027280 self | 1 | 37.147 | 42.024 | 13.13% | 見送り |
| AP027280 self | 8 | 61.850 | 65.333 | 5.63% | 見送り |

全32回のraw bytesは一致したが、全4条件で変更前より遅くなった。1.20倍未達を理由に棄却したわけではない。
統合版へRust変更を取り込んでいない。棄却後の全runtime統合計測も行っていない。
X1のNative速度・反復速度・ブラウザ速度は採用基準として未測定であり、合格扱いしない。

## Native absolute times and ordinary threaded Wasm

The same frozen Rust inputs and matched body interval are used. Native has one warmup plus three timed runs; Wasm uses the baseline side of the three alternating pairs. Native and Wasm runs are sequential, not simultaneous. Values below are medians in seconds.

| Input | n | Native body / process | Wasm body / process | Wasm / Native body |
|---|---:|---:|---:|---:|
| MjeNMV.MelaMJNV.tlosatx | 1 | 14.295 / 14.299 | 22.406 / 22.500 | 1.567 |
| MjeNMV.MelaMJNV.tlosatx | 8 | 8.346 / 8.352 | 35.077 / 35.177 | 4.203 |
| AP027280.AP027280.tlosatx | 1 | 27.108 / 27.112 | 37.147 / 37.237 | 1.370 |
| AP027280.AP027280.tlosatx | 8 | 18.127 / 18.131 | 61.850 / 61.953 | 3.412 |

Compile/module reuse and same-instance timings, setup cost and memory evidence are reported separately in [REUSE.md](REUSE.md). All24 per-session/case guards pass, and the linear-memory sequences are identical.

## 原因と残る候補

[DIAGNOSIS.md](DIAGNOSIS.md)に仕事量、NCBI owner、生成コード、条件付き候補を記録した。
n1/n2/n4/n8でlarge-gap訪問数は一致し、両linking段階の総数は約40.48億／99.08億だった。
通常設定ではn8が遅いが、両版に同じTurboFan固定flagを与えた補助計測では両入力で逆転が消えた。
これはコンパイル段階を主要因として強く支持する。個々のgroupと実行tierの厳密な対応、
全待ち時間の因果分解までは確定していない。

- X2：過大capacityだけでは物理メモリ・初期化時間の主要費用を証明できず見送り。
- M1：非圧縮一致走査で16塩基以上一致する呼出しは0.070–0.134%。SIMD導入の根拠不足で見送り。
- H1b：raw検査約4 ms、共有メモリguard約47 ms、guard後compile約2.5 ms。guardを維持し、大きな新規parserは追加しない。
- C1局所スケジューリング案：同一module反復では毎回直列prefixが増えるため、別の条件付き案として記録。未実装・未採用。

## 検証範囲

- B0：fresh oracleに対し15 fixtures × Native n1/n8・serial n1・threaded n1/n8の75条件がraw一致。
- X1 focused：35条件raw一致。長いAP027131/AP027133 code4はfresh DB oracleと5条件一致。
- X1状態：176合成cases、Native 19,017行、同点補強20,159行、Wasm 19,017行が各target内の旧新で一致。
- X1 trace：sum拒否経路を含む132行がNative/Wasmそれぞれ一致。段階出力はn1/n2/n4/n8の36 file pairsが一致。
- H1 host：19 tests PASS。正しい／不正なABI、別instanceのmemory独立性、終了0/23/70、Worker trap、実shared-memory拡張を含む。
- H1最終比較：15 fixturesのserial raw出力・終了状態が一致。四処理、短い閾値、word7、single32、有効no-hitを含む。
- H1実artifact：command/reactor × serial/threadedの検査CLI JSONが旧新でbyte一致。serial reactorの成功・拒否・回復も一致。
- H1反復：AB/BAの独立2 sessions、各case warmup 1＋測定3。compiled-module command/reactorとsame-instance reactorの192 invocationsがraw一致。

反復はこの有限windowに限る。既存reactor memory問題や無制限反復の安定性は解決・認証していない。
統計的に無効なTBLASTX queryのexit-status差も未解決。
PR5 Gate A／公式platform Gate Bの認証は本ラウンドでは実施していない。

## Rust・最終差分・失敗記録

- `V0-cargo-test`: exit 101。全出力は実験archiveの`V0-cargo-test.log`。
- `V0-cargo-clippy`: exit 0。全出力は実験archiveの`V0-cargo-clippy.log`。
- `V0-cargo-fmt`: exit 0。全出力は実験archiveの`V0-cargo-fmt.log`。
- リポジトリ指定の標準 `cargo test`（独立target）：625 passed / 0 failed / 3 ignored。`V0-cargo-test-standard.log`を参照。

最初のreleaseテストは隔離snapshotのfixture不足で停止した。元のfixtureを追加した後、
既存のrlib/cdylib出力名衝突とpanic=abort / unwind不一致でreleaseテストのコンパイルが失敗した。
release profileは変更せず、指定の標準cargo testを独立targetで実行した。

Rustの最終gateは変更していないbaselineの凍結入力に対して実行した。
ホスト変更はH1-sourceのhashと作業ツリーのhashを照合して引き渡す。
production Rust、Nodeホスト、検証、文書、実験artifactを分けてレビューした。

最初のNode host試験はsandboxのspawnSync EPERMで停止した。同じ失敗を変更前で確認し、
同じ候補テストを制限外で実行して19件PASS。失敗ログを保存し、成功sampleへ混ぜていない。
最初のX1準備はassertionで停止し、未使用build名X1は評価対象外。正式な実験名はX1v1。
最初のtrace座標は対象がなく0行だった。実在する座標211,234,217,240に直して比較した。

全コマンド・raw出力・warmup・測定sample・source/build/artifact hashは
[README.md](README.md)とbundle manifestから辿れる。過去の採用判断とexpected出力は変更していない。
