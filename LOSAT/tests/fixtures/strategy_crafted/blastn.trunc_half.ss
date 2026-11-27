Blast4-request ::= {
  body queue-search {
    program "blastn",
    service "megablast",
    queries bioseq-set {
      seq-set {
        seq {
          id {
            local str "Query_1"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">qn1 query nucleotide one"
                }
              }
            },
            title "qn1 query nucleotide one"
          },
          inst {
            repr raw,
            mol na,
            length 480,
            seq-data ncbi2na '049131660FBCB5C3BA3A0475B9C0240FF96E7DB02E3A1856
00C6A3A953A7B0E0A9EA2073E7A5660CE27A02D7CE0B38BFD758CD457E43F3A3F9E8D5C288E905
528EBA0868C6F933A480EADEE8E4C2CD50BB3BBB0EB88F448D0DDD693E803C14E19C1FE201F08D
F51C4EBAB4EA2049378A90410CB56CA6A831'H
          }
        },
        seq {
          id {
            local str "Query_