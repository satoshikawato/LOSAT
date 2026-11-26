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
            title "q1"
          },
          inst {
            repr raw,
            mol na,
            length 300,
            seq-data ncbi2na '2D372C8246E6A8BF35F71B036DE048ABC5BF017A0D5200F026A81F5BFFF8208F77A57A37FB62C1B8C4BA0A6D734F715906E2F92AD0246227E12208DB5DD3A4DE29D5FE7B0A3BA2EDD8B694'H
          }
        }
      }
    },
    subject database "tests/db_search/arch_snucl",
    algorithm-options {
      {
        name "WindowSize",
        value integer -40
      },
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, 1 }
      },
      {
        name "UngappedMode",
        value boolean FALSE
      },
      {
        name "HitlistSize",
        value integer 500
      }
    }
  }
}
