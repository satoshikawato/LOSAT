# NCBI BLAST TBLASTX実装の詳細

このドキュメントは、NCBI BLASTのTBLASTX実装について、ソースコードから確認した内容を網羅的にまとめたものです。

## 目次

1. [Two-hit要件の実装](#two-hit要件の実装)
2. [cutoff_scoreとgap_trigger](#cutoff_scoreとgap_trigger)
3. [Extensionの実装](#extensionの実装)
4. [E-value計算とEffective Search Space](#e-value計算とeffective-search-space)
5. [Diagonal Suppression](#diagonal-suppression)
6. [Seed Score Threshold](#seed-score-threshold)
7. [One-hit vs Two-hit Extension](#one-hit-vs-two-hit-extension)
8. [重要な発見と疑問点](#重要な発見と疑問点)

---

## Two-hit要件の実装

### 基本動作

NCBI BLASTでは、`ewp->diag_table->multiple_hits`が`TRUE`の場合、`s_BlastAaWordFinder_TwoHit`が使用されます。

**参照**: `ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:214-232`

```c
if (ewp->diag_table->multiple_hits) {
    status = s_BlastAaWordFinder_TwoHit(subject, query, ...);
} else {
    status = s_BlastAaWordFinder_OneHit(subject, query, ...);
}
```

### Two-hit要件のチェック

**参照**: `ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:532-573`

```c
/* find the distance to the last hit on this diagonal */
last_hit = diag_array[diag_coord].last_hit - diag_offset;
diff = subject_offset - last_hit;

if (diff >= window) {
    /* We are beyond the window for this diagonal; start a new hit */
    diag_array[diag_coord].last_hit = subject_offset + diag_offset;
    continue;  // Extensionは実行されない
}

/* If the difference is less than the wordsize (i.e. last
   hit and this hit overlap), give up */
if (diff < wordsize) {
    continue;  // Extensionは実行されない
}

/* Check if the last hit hits current query. Because last_hit
   is never reset, it may contain a hit to the previous
   concatenated query */
if (query_offset - diff <
    query_info->contexts[curr_context].query_offset) {
    /* there was no last hit for this diagnol; start a new hit */
    diag_array[diag_coord].last_hit = subject_offset + diag_offset;
    continue;  // Extensionは実行されない
}

// Two-hit要件を満たした場合のみ、extensionが実行される
score = s_BlastAaExtendTwoHit(...);
```

### 重要な点

1. **Two-hit要件を満たさない場合**: `continue`でスキップされ、extensionは**実行されない**
2. **Two-hit要件を満たす条件**:
   - `diff < window` (window内に前のhitがある)
   - `diff >= wordsize` (重複していない)
   - 前のhitが現在のqueryに属している

### 疑問点

**現在のLOSATの実装では、Two-hit要件を満たさない場合でもextensionを実行していますが、これはNCBI BLASTの実装と一致していません。**

NCBI BLASTのコードを見る限り、Two-hit要件を満たさない場合は`continue`でスキップされ、extensionは実行されません。

---

## cutoff_scoreとgap_trigger

### cutoff_scoreの計算

**参照**: `ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:348-374`

```c
if (!gapped_calculation || sbp->matrix_only_scoring) {
    double cutoff_e = s_GetCutoffEvalue(program_number);
    Int4 query_length = query_info->contexts[context].query_length;
    
    kbp = kbp_array[context];
    BLAST_Cutoffs(&new_cutoff, &cutoff_e, kbp, 
                  MIN((Uint8)subj_length, 
                      (Uint8)query_length)*((Uint8)subj_length), 
                  TRUE, gap_decay_rate);

    /* Perform this check for compatibility with the old code */
    if (program_number != eBlastTypeBlastn)  
        new_cutoff = MIN(new_cutoff, gap_trigger);
} else {
    new_cutoff = gap_trigger;
}
new_cutoff *= (Int4)sbp->scale_factor;
new_cutoff = MIN(new_cutoff, 
                 hit_params->cutoffs[context].cutoff_score_max);
curr_cutoffs->cutoff_score = new_cutoff;
```

### CUTOFF_E_TBLASTX

**参照**: `ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:134-156`

```c
s_GetCutoffEvalue(EBlastProgramType program)
{
   switch(program) {
   case eBlastTypeTblastx:
      return CUTOFF_E_TBLASTX;  // 1e-300
   ...
   }
}
```

**参照**: `ncbi-blast/c++/include/algo/blast/core/blast_parameters.h:80`
```c
#define CUTOFF_E_TBLASTX 1e-300
```

### gap_triggerの計算

**参照**: `ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:340-346`

```c
Int4 gap_trigger = INT4_MAX;
if (sbp->kbp_std) {
    kbp = sbp->kbp_std[context];
    if (s_BlastKarlinBlkIsValid(kbp)) {
        gap_trigger = (Int4)((kOptions->gap_trigger * NCBIMATH_LN2 + 
                                 kbp->logK) / kbp->Lambda);
    }
}
```

**参照**: `ncbi-blast/c++/include/algo/blast/core/blast_options.h:137`
```c
#define BLAST_GAP_TRIGGER_PROT 22.0
```

### cutoff_scoreの適用

**参照**: `ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:587-591`

```c
/* if the hsp meets the score threshold, report it */
if (score >= cutoffs->cutoff_score)
    BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
                     query_offset, subject_offset, hsp_len,
                     score);
```

### 重要な点

1. **cutoff_scoreの計算**: `BLAST_Cutoffs`関数で、`CUTOFF_E_TBLASTX = 1e-300`から計算
2. **search space**: `MIN(subj_length, query_length) * subj_length`を使用
3. **gap_triggerの制限**: `program_number != eBlastTypeBlastn`の場合、`new_cutoff = MIN(new_cutoff, gap_trigger)`
4. **適用タイミング**: Extension実行後、E-valueチェック前

### 重要な理解

1. **cutoff_scoreのsearch space**: `MIN(subj_length, query_length) * subj_length`を使用
   - **目的**: Ungapped extensionの事前フィルタリング用
   - **E-value計算との違い**: E-value計算では`effective_db_length * (query_length - length_adjustment)`を使用
   - **結論**: これは意図的な違い。cutoff_scoreは事前フィルタリング用、E-valueは最終的な統計評価用

2. **gap_triggerの制限**: `gap_trigger`は`BLAST_GAP_TRIGGER_PROT = 22.0`（bit score）から計算され、Karlin-Altschulパラメータを用いて最終的な値に変換される
   - **現在のLOSAT**: `gap_trigger=46`が`cutoff_score`を制限
   - **適切性**: 使用されるスコア行列やパラメータに依存する

---

## Extensionの実装

### Two-hit Extension

**参照**: `ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:1088-1158`

```c
static Int4
s_BlastAaExtendTwoHit(Int4 ** matrix,
                      const BLAST_SequenceBlk * subject,
                      const BLAST_SequenceBlk * query,
                      Int4 s_left_off, Int4 s_right_off, 
                      Int4 q_right_off, Int4 dropoff, 
                      Int4 * hsp_q, Int4 * hsp_s, Int4 * hsp_len,
                      Boolean use_pssm, Int4 word_size, 
                      Boolean * right_extend, Int4 * s_last_off)
{
    // 1. Word内で最良のスコア位置を見つける (lines 1108-1119)
    for (i = 0; i < word_size; i++) {
        score += matrix[q[q_right_off + i]][s[s_right_off + i]];
        if (score > left_score) {
            left_score = score;
            right_d = i + 1;
        }
    }
    q_right_off += right_d;
    s_right_off += right_d;

    // 2. Left extensionを実行 (lines 1127-1135)
    left_score = s_BlastAaExtendLeft(matrix, subject, query,
                                     s_right_off - 1, q_right_off - 1,
                                     dropoff, &left_d, 0);

    // 3. Left extensionがfirst hitに到達した場合のみright extensionを実行 (line 1138)
    if (left_d >= (s_right_off - s_left_off)) {
        *right_extend = TRUE;
        right_score = s_BlastAaExtendRight(matrix, subject, query,
                                           s_right_off, q_right_off,
                                           dropoff, &right_d, left_score,
                                           s_last_off);
    }

    // 4. 常に結果を返す (line 1157)
    return MAX(left_score, right_score);
}
```

### 重要な点

1. **Left extensionは常に実行される**: Two-hit要件を満たした場合、left extensionは必ず実行される
2. **Right extensionの条件**: Left extensionがfirst hitに到達した場合のみ実行される
3. **戻り値**: `MAX(left_score, right_score)` - right extensionが実行されない場合、`left_score`が返される

### 重要な理解

**Two-hit要件を満たさない場合の処理**

NCBI BLASTのコードでは、Two-hit要件を満たさない場合は`continue`でスキップされ、`s_BlastAaExtendTwoHit`は呼び出されません。つまり、extensionは実行されません。

**Two-hit要件を満たした場合の動作**

Two-hit要件を満たした場合、`s_BlastAaExtendTwoHit`が呼び出されます。この関数内では：

1. **Left extensionは常に実行される**: Two-hit要件を満たした場合でも、left extensionがfirst hitに到達しない場合があります
2. **Right extensionの条件**: Left extensionがfirst hitに到達した場合のみ実行される
3. **戻り値**: `MAX(left_score, right_score)` - right extensionが実行されない場合、`left_score`が返される
4. **HSP保存条件**: スコアが`cutoff_score`以上であればHSPとして保存される

つまり、「Two-hitを満たさない場合でも、left extensionの結果がcutoff_scoreを超えればHSPとして保存される」という情報は、`s_BlastAaExtendTwoHit`関数内での動作を説明している可能性があります。しかし、NCBI BLASTのコードを見る限り、Two-hit要件を満たさない場合はextension自体が実行されません。

---

## E-value計算とEffective Search Space

### Effective Search Spaceの計算

**参照**: `ncbi-blast/c++/src/algo/blast/core/blast_setup.c:730-848`

```c
if (Blast_SubjectIsTranslated(program_number))
    db_length = db_length/3;  // TBLASTXの場合、amino acid lengthに変換

for (index = query_info->first_context;
     index <= query_info->last_context;
     index++) {
    Int4 query_length = query_info->contexts[index].query_length;
    
    // Length adjustmentの計算
    BLAST_ComputeLengthAdjustment(kbp->K, kbp->logK,
                                  alpha/kbp->Lambda, beta,
                                  query_length, db_length,
                                  db_num_seqs, &length_adjustment);

    // Effective search spaceの計算
    Int8 effective_db_length = db_length - ((Int8)db_num_seqs * length_adjustment);
    if (effective_db_length <= 0)
        effective_db_length = 1;

    effective_search_space = effective_db_length *
                            (query_length - length_adjustment);
    
    query_info->contexts[index].eff_searchsp = effective_search_space;
    query_info->contexts[index].length_adjustment = length_adjustment;
}
```

### 重要な点

1. **各contextごとに計算**: 各query frameとsubject frameの組み合わせごとに異なるeffective search spaceを計算
2. **Length adjustment**: Karlin-Altschul統計パラメータから計算
3. **TBLASTXの場合**: `db_length = db_length/3`でamino acid lengthに変換

### 疑問点

**cutoff_scoreの計算で使用されるsearch spaceと、E-value計算で使用されるsearch spaceが異なる**

- **cutoff_score**: `MIN(subj_length, query_length) * subj_length`
- **E-value**: `effective_db_length * (query_length - length_adjustment)`

この不一致は意図的か、それとも実装の違いか？

---

## Diagonal Suppression

### diag_maskの計算

**参照**: `ncbi-blast/c++/src/algo/blast/core/blast_extend.c:42-67`

```c
static BLAST_DiagTable*
s_BlastDiagTableNew (Int4 qlen, Boolean multiple_hits, Int4 window_size)
{
    Int4 diag_array_length = 1;
    /* What power of 2 is just longer than the query? */
    while (diag_array_length < (qlen+window_size))
    {
        diag_array_length = diag_array_length << 1;
    }
    diag_table->diag_array_length = diag_array_length;
    diag_table->diag_mask = diag_array_length-1;
    ...
}
```

### diag_coordの計算

**参照**: `ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:351`

```c
diag_coord = (query_offset - subject_offset) & diag_mask;
```

### last_hitの更新

**参照**: `ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:593-606`

```c
/* If an extension to the right happened, reset the last hit
   so that future hits to this diagonal must start over. */
if (right_extend) {
    diag_array[diag_coord].flag = 1;
    diag_array[diag_coord].last_hit =
        s_last_off - (wordsize - 1) + diag_offset;
}
/* Otherwise, make the present hit into the previous hit for
   this diagonal */
else {
    diag_array[diag_coord].last_hit =
        subject_offset + diag_offset;
}
```

### 重要な点

1. **diag_mask**: `diag_array_length - 1`（2のべき乗 - 1）
2. **diag_coord**: `(query_offset - subject_offset) & diag_mask`
3. **last_hitの更新**: Extension実行後、スコアチェック後に更新

---

## Seed Score Threshold

**参照**: `ncbi-blast/c++/include/algo/blast/core/blast_options.h:115`

```c
#define BLAST_WORD_THRESHOLD_TBLASTX 13
```

**参照**: `ncbi-blast/c++/src/algo/blast/core/blast_options.c:1204`

TBLASTXの場合、base threshold (11)に2を加えて13を使用。

---

## One-hit vs Two-hit Extension

### 選択条件

**参照**: `ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:214-232`

```c
if (ewp->diag_table->multiple_hits) {
    status = s_BlastAaWordFinder_TwoHit(...);
} else {
    status = s_BlastAaWordFinder_OneHit(...);
}
```

### multiple_hitsの設定

TBLASTXの場合、`window_size > 0`であれば`multiple_hits = TRUE`となり、Two-hit要件が使用されます。

### 重要な点

1. **TBLASTXは常にTwo-hit要件を使用**: `window_size = 40`のため、`multiple_hits = TRUE`
2. **One-hit extensionは使用されない**: TBLASTXでは`s_BlastAaWordFinder_OneHit`は呼ばれない

---

## 重要な発見と疑問点

### 確認済みの実装

1. ✅ **Two-hit要件**: `diff < window && diff >= wordsize`でチェック
2. ✅ **cutoff_score**: `CUTOFF_E_TBLASTX = 1e-300`から計算、`gap_trigger`で制限
3. ✅ **Extension**: Left extensionは常に実行、right extensionは条件付き
4. ✅ **Effective search space**: 各contextごとに計算
5. ✅ **Diagonal suppression**: `diag_coord = (query_offset - subject_offset) & diag_mask`
6. ✅ **Seed score threshold**: `BLAST_WORD_THRESHOLD_TBLASTX = 13`

### 確認済みの実装詳細

1. ✅ **Two-hit要件を満たさない場合のextension**: 
   - **NCBI BLASTのコード**: `continue`でスキップ（extension実行されない）
   - **確認**: `aa_ungapped.c:538-543`, `aa_ungapped.c:549-551`, `aa_ungapped.c:566-573`
   - **結論**: Two-hit要件を満たさない場合は、`s_BlastAaExtendTwoHit`は呼び出されない

2. ✅ **cutoff_scoreのsearch space**:
   - **cutoff_score計算**: `MIN(subj_length, query_length) * subj_length` (`blast_parameters.c:360-363`)
   - **E-value計算**: `effective_db_length * (query_length - length_adjustment)`
   - **結論**: 意図的な違い。cutoff_scoreはungapped extensionの事前フィルタリング用、E-valueは最終的な統計評価用

3. ✅ **gap_triggerの計算**:
   - **計算方法**: `BLAST_GAP_TRIGGER_PROT = 22.0`（bit score）から、Karlin-Altschulパラメータを用いて計算
   - **制限**: `program_number != eBlastTypeBlastn`の場合、`new_cutoff = MIN(new_cutoff, gap_trigger)`
   - **現在のLOSAT**: `gap_trigger=46`が`cutoff_score`を制限
   - **適切性**: 使用されるスコア行列やパラメータに依存

4. ❓ **低アイデンティティヒットの検出**:
   - **LOSAT**: 21-28%の比率（NCBI BLAST比）
   - **NCBI BLAST**: 42-63%の比率
   - **問題**: なぜこれほど差があるのか？要調査

### 次の調査項目

1. ✅ **Two-hit要件を満たさない場合の処理**: 確認済み - NCBI BLASTではextensionは実行されない
2. ✅ **cutoff_scoreの計算方法**: 確認済み - `BLAST_Cutoffs`関数で計算、search spaceは`MIN(subj_length, query_length) * subj_length`
3. ✅ **gap_triggerの適用方法**: 確認済み - `program_number != eBlastTypeBlastn`の場合に制限
4. ❓ **低アイデンティティヒットの検出**: 要調査 - なぜLOSATは21-28%の比率しか検出できないのか？
5. ❓ **gap_trigger=46の適切性**: 要調査 - 91.8%のヒットが`cutoff_score`でフィルタリングされている原因
6. ❓ **cutoff_scoreでフィルタリングされているヒットのスコア分布**: 要分析 - 低アイデンティティヒットが失われている原因を特定

---

## 参照ファイル

### NCBI BLASTソースコード

- `ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c`: Two-hit要件とextensionの実装
- `ncbi-blast/c++/src/algo/blast/core/blast_parameters.c`: cutoff_scoreとgap_triggerの計算
- `ncbi-blast/c++/src/algo/blast/core/blast_setup.c`: Effective search spaceの計算
- `ncbi-blast/c++/src/algo/blast/core/blast_extend.c`: Diagonal tableの初期化
- `ncbi-blast/c++/include/algo/blast/core/blast_options.h`: 定数定義
- `ncbi-blast/c++/include/algo/blast/core/blast_parameters.h`: 定数定義

### LOSAT実装

- `LOSAT/src/algorithm/tblastx/utils.rs`: メイン実行ロジック
- `LOSAT/src/algorithm/tblastx/extension.rs`: Extensionロジック
- `LOSAT/src/algorithm/tblastx/constants.rs`: 定数定義

---

## 更新履歴

- 2024-12-28: 初版作成。NCBI BLASTの実装を網羅的にまとめ、疑問点を明確化。

