# Session S10 — W2：データ層

## INSTRUCTION PROMPT

LOSAT Web の段階 W2 を実行する。先に [セッション README](README.md) の共通規則を読み、それに従う。完了条件の正本は、総合計画書 §7 の S10 の行である。設計は計画 §5.4 と §5.6 にある。

目的：入力ファイルと実行結果を、Engine worker から独立して安全に保持する（設計書 §4、§6、§7）。

1. Data worker を作り、`DataGateway`（`src/ports/data.ts`）をその中で実装する。S01 の `MemoryDataGateway` は、同じ契約の Memory 実装として Data worker の中へ移す。
2. OPFS の実装を作る。同期アクセスハンドルは Data worker の中だけで使う。配置は `tmp/<session-token>/runs/<run-token>/…` とし、ファイル名に入力名を入れない。Memory と OPFS の両方に、同じ契約試験（書き込み、確定、破棄、読み出し、大きなブロック）をかける。OPFS が使えないときは Memory にする。経路はブラウザ名ではなく、機能を試して選ぶ。
3. 入力：File は参照だけを持ち、`File.slice` で読む。ABI v2 の `scan_*` でレコード表を作り、レコードごとの SHA-256 と DatasetRevision を作る。`register` の時点で、エンジンが解析した ID・長さとレコード表を照合し、食い違えば止める。
4. 実行の結果は、確定するまで「仮」とし、取消・失敗なら破棄する。Engine worker から Data worker へ結果のチャンクを直接送る経路（worker どうしの MessagePort）を作る。
5. 作業ごとに `navigator.locks` のロックを持ち続け、起動時には誰もロックを持っていない `tmp/*` だけを消す。Web Locks が使えない環境では消さない。
6. 容量不足（`QuotaExceededError`）が起きたら、その実行を理由付きで失敗にし、過去の結果を残す。使用量を、操作を妨げない表示で示す。
7. 試験：タブを強制終了した後の回収、2 つのタブで互いのデータを消さないこと、停止中のタブのデータを消さないこと、容量不足を注入した場合、を Playwright で確かめる。

完了条件は計画 §7 の S10 の行による。記録は `docs/evidence/losat_web_w2/README.md`。

## 終了・引き継ぎ

README の規則 8 に従う。次は [S11 — 領域の指定](session_s11_e2d_query_subject_loc.md)。
