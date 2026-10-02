Wall seconds / CPU seconds (user+sys) per tblastx label, from frozen-regressions/runs.json

| label | green 35201795582 | red 36745543144 | red 36882087513 | wall ratio red1/green | red2/green |
|---|---|---|---|---|---|
| native/tblastx/d01_nz_self_code4 | 1497s / cpu 1055 | 2327s / cpu 1503 | 1898s / cpu 1347 | 1.55 | 1.27 |
| threaded/tblastx/d01_nz_self_code4 | 850s / cpu 1286 | 1313s / cpu 1809 | 1125s / cpu 1672 | 1.55 | 1.32 |
| native/tblastx/d02_ap027132_nz_code4 | 1737s / cpu 1412 | 2389s / cpu 2017 | 2356s / cpu 1797 | 1.38 | 1.36 |
| threaded/tblastx/d02_ap027132_nz_code4 | 884s / cpu 1678 | 1262s / cpu 2376 | 1149s / cpu 2173 | 1.43 | 1.30 |
| native/tblastx/d03_ap027078_ap027131_code4 | 736s / cpu 580 | 1013s / cpu 898 | 835s / cpu 768 | 1.38 | 1.13 |
| threaded/tblastx/d03_ap027078_ap027131_code4 | 321s / cpu 676 | 513s / cpu 1020 | 506s / cpu 898 | 1.60 | 1.57 |
| native/tblastx/d04_ap027131_ap027133_code4 | 174s / cpu 161 | 232s / cpu 208 | 279s / cpu 207 | 1.33 | 1.60 |
| threaded/tblastx/d04_ap027131_ap027133_code4 | 104s / cpu 217 | 128s / cpu 284 | 148s / cpu 279 | 1.23 | 1.42 |
| native/tblastx/d05_ap027133_ap027132_code4 | 869s / cpu 773 | 1309s / cpu 1186 | 1140s / cpu 1045 | 1.51 | 1.31 |
| threaded/tblastx/d05_ap027133_ap027132_code4 | 534s / cpu 1076 | 792s / cpu 1603 | 692s / cpu 1400 | 1.48 | 1.30 |
| native/tblastx/d06_ap027131_ap027133_db4 | 36s / cpu 36 | 51s / cpu 45 | 48s / cpu 48 | 1.43 | 1.34 |
| threaded/tblastx/d06_ap027131_ap027133_db4 | 37s / cpu 56 | 45s / cpu 71 | 45s / cpu 71 | 1.21 | 1.23 |
| native/tblastx/p01_ap027280_self | 48s / cpu 43 | 63s / cpu 56 | 59s / cpu 52 | 1.31 | 1.23 |
| threaded/tblastx/p01_ap027280_self | 39s / cpu 61 | 47s / cpu 76 | 54s / cpu 86 | 1.19 | 1.38 |
| native/tblastx/p02_mje_mela | 21s / cpu 21 | 31s / cpu 29 | 29s / cpu 27 | 1.46 | 1.36 |
| threaded/tblastx/p02_mje_mela | 22s / cpu 35 | 28s / cpu 44 | 31s / cpu 47 | 1.27 | 1.39 |
| native/tblastx/p03_mela_pemojnva | 3s / cpu 3 | 6s / cpu 6 | 6s / cpu 6 | 1.86 | 1.81 |
| threaded/tblastx/p03_mela_pemojnva | 6s / cpu 9 | 7s / cpu 11 | 8s / cpu 12 | 1.24 | 1.34 |
| native/tblastx/p04_pemojnva_pesemjnv | 23s / cpu 22 | 27s / cpu 26 | 33s / cpu 31 | 1.18 | 1.43 |
| threaded/tblastx/p04_pemojnva_pesemjnv | 31s / cpu 47 | 37s / cpu 63 | 41s / cpu 67 | 1.18 | 1.32 |
| native/tblastx/p05_pesemjnv_pemojnva | 61s / cpu 59 | 69s / cpu 63 | 78s / cpu 72 | 1.12 | 1.27 |
| threaded/tblastx/p05_pesemjnv_pemojnva | 68s / cpu 113 | 99s / cpu 157 | 110s / cpu 175 | 1.46 | 1.61 |
| native/tblastx/p06_pemojnva_lvmjnv | 452s / cpu 299 | 438s / cpu 299 | 567s / cpu 373 | 0.97 | 1.25 |
| threaded/tblastx/p06_pemojnva_lvmjnv | 515s / cpu 1179 | 596s / cpu 1329 | 690s / cpu 1560 | 1.16 | 1.34 |
| native/tblastx/p07_lvmjnv_trcumjnv | 6s / cpu 3 | 5s / cpu 3 | 6s / cpu 4 | 0.95 | 1.12 |
| threaded/tblastx/p07_lvmjnv_trcumjnv | 7s / cpu 6 | 8s / cpu 6 | 10s / cpu 7 | 1.07 | 1.40 |
| native/tblastx/p08_trcumjnv_mellatmjnv | 22s / cpu 14 | 27s / cpu 18 | 28s / cpu 18 | 1.22 | 1.27 |
| threaded/tblastx/p08_trcumjnv_mellatmjnv | 24s / cpu 28 | 26s / cpu 34 | 31s / cpu 38 | 1.09 | 1.27 |
| native/tblastx/p09_mellatmjnv_meenmjnv | 96s / cpu 62 | 121s / cpu 84 | 126s / cpu 84 | 1.26 | 1.31 |
| threaded/tblastx/p09_mellatmjnv_meenmjnv | 84s / cpu 106 | 108s / cpu 129 | 116s / cpu 138 | 1.28 | 1.38 |
| native/tblastx/p10_meenmjnv_mejomjnv | 157s / cpu 107 | 238s / cpu 145 | 229s / cpu 143 | 1.52 | 1.46 |
| threaded/tblastx/p10_meenmjnv_mejomjnv | 134s / cpu 242 | 156s / cpu 246 | 177s / cpu 308 | 1.16 | 1.32 |
| native/tblastx/p11_avclpv_psclpv | 1909s / cpu 1248 | 1413s / cpu 997 | 2060s / cpu 1472 | 0.74 | 1.08 |
| threaded/tblastx/p11_avclpv_psclpv | 3358s / cpu 6617 | 3600s **TIMEOUT** | 3600s **TIMEOUT** | 1.07 | 1.07 |
| native/tblastx/p12_lc738874_lc738875_default | 6s / cpu 6 | 7s / cpu 7 | 12s / cpu 8 | 1.21 | 1.92 |
| threaded/tblastx/p12_lc738874_lc738875_default | 11s / cpu 17 | 14s / cpu 17 | 17s / cpu 20 | 1.36 | 1.63 |
| native/tblastx/p13_mela_mje_reverse | 18s / cpu 18 | 35s / cpu 23 | 32s / cpu 23 | 1.98 | 1.79 |
| threaded/tblastx/p13_mela_mje_reverse | 18s / cpu 27 | 31s / cpu 34 | 34s / cpu 35 | 1.72 | 1.86 |
| native/tblastx/p14_ap027131_ap027133_query4 | 38s / cpu 38 | 73s / cpu 49 | 74s / cpu 50 | 1.93 | 1.97 |
| threaded/tblastx/p14_ap027131_ap027133_query4 | 37s / cpu 57 | 68s / cpu 73 | 52s / cpu 75 | 1.83 | 1.40 |

threaded/native wall ratio per case (green / red1 / red2):
  d01_nz_self_code4: 0.57 / 0.56 / 0.59
  d02_ap027132_nz_code4: 0.51 / 0.53 / 0.49
  d03_ap027078_ap027131_code4: 0.44 / 0.51 / 0.61
  d04_ap027131_ap027133_code4: 0.60 / 0.55 / 0.53
  d05_ap027133_ap027132_code4: 0.61 / 0.60 / 0.61
  d06_ap027131_ap027133_db4: 1.02 / 0.87 / 0.93
  p01_ap027280_self: 0.82 / 0.74 / 0.92
  p02_mje_mela: 1.04 / 0.90 / 1.06
  p03_mela_pemojnva: 1.82 / 1.21 / 1.35
  p04_pemojnva_pesemjnv: 1.36 / 1.36 / 1.25
  p05_pesemjnv_pemojnva: 1.12 / 1.45 / 1.41
  p06_pemojnva_lvmjnv: 1.14 / 1.36 / 1.22
  p07_lvmjnv_trcumjnv: 1.31 / 1.47 / 1.63
  p08_trcumjnv_mellatmjnv: 1.08 / 0.97 / 1.09
  p09_mellatmjnv_meenmjnv: 0.88 / 0.89 / 0.92
  p10_meenmjnv_mejomjnv: 0.85 / 0.65 / 0.77
  p11_avclpv_psclpv: 1.76 / 2.55 / 1.75
  p12_lc738874_lc738875_default: 1.72 / 1.93 / 1.46
  p13_mela_mje_reverse: 1.03 / 0.90 / 1.07
  p14_ap027131_ap027133_query4: 0.99 / 0.94 / 0.71
