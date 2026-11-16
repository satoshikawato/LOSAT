Blast4-request ::= {
  body queue-search {
    program "blastp",
    service "plain",
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
                  data str ">q1 one"
                }
              }
            },
            title "q1 one"
          },
          inst {
            repr raw,
            mol aa,
            length 120,
            seq-data ncbieaa "MWGIKSVTGYTMTSRLDLQLTSCYQKMVIECMLHDNEMITKMANTLVP
PGECNVEDFDEMACQKNHHPRYSNWVRLSDASPPNEIAVRMGSYSCDYLRKQTQKSYSPAPMMGEGMYLYGN"
          }
        }
      }
    },
    subject database "xpost/db/pR3/pR3",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, 1 }
      },
      {
        name "MaskAtHash",
        value boolean FALSE
      },
      {
        name "SegFiltering",
        value boolean FALSE
      },
      {
        name "WordThreshold",
        value integer 11
      },
      {
        name "WindowSize",
        value integer 40
      },
      {
        name "HitlistSize",
        value integer 500
      },
      {
        name "UngappedMode",
        value boolean FALSE
      },
      {
        name "CompositionBasedStats",
        value integer 2
      },
      {
        name "SmithWatermanMode",
        value boolean FALSE
      }
    }
  }
}
Blast4-request ::= {
  body queue-search {
    program "blastp",
    service "plain",
    queries bioseq-set {
      seq-set {
        seq {
          id {
            local str "Query_2"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">q2 a\303\251b z"
                }
              }
            },
            title "q2 a##b z"
          },
          inst {
            repr raw,
            mol aa,
            length 120,
            seq-data ncbieaa "MWGIKSVTGYTMTSRLDLQLTSCYQKMVIECMLHDNEMITKMANTLVP
PGECNVEDFDEMACQKNHHPRYSNWVRLSDASPPNEIAVRMGSYSCDYLRKQTQKSYSPAPMMGEGMYLYGN"
          }
        }
      }
    },
    subject database "xpost/db/pR3/pR3",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, 1 }
      },
      {
        name "MaskAtHash",
        value boolean FALSE
      },
      {
        name "SegFiltering",
        value boolean FALSE
      },
      {
        name "WordThreshold",
        value integer 11
      },
      {
        name "WindowSize",
        value integer 40
      },
      {
        name "HitlistSize",
        value integer 500
      },
      {
        name "UngappedMode",
        value boolean FALSE
      },
      {
        name "CompositionBasedStats",
        value integer 2
      },
      {
        name "SmithWatermanMode",
        value boolean FALSE
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 11449
      }
    }
  }
}
