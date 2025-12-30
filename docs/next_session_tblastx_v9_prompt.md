# 次セッション開始時のプロンプト

```
@docs/next_session_tblastx_v9.md に従って、TBLASTXのSEGフィルター実装を進めましょう。

NCBI BLASTのコードベースはこちら:
C:\Users\kawato\Documents\GitHub\ncbi-blast

現状の問題:
- LOSATはNZ_CP006932 self-comparisonで最大113,828 aaまで拡張（NCBIは2,260 aa）
- センチネル実装は完了したが、runaway extension問題は解決していない
- NCBI BLASTのSEGフィルター（-segオプション）を実装する必要がある

目標:
1. NCBI BLASTのSEG実装を調査
2. LOSATにSEGフィルターを実装
3. NZ_CP006932 self-comparisonで最大アラインメント長が2,500 aa以下になることを確認

がんばろう！
```


