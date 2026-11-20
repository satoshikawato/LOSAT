Blast4-request ::= {
  body queue-search {
    program "blastn",
    service "megablast",
    queries seq-loc-list {
      int {
        from 0,
        to 479,
        strand both,
        id genbank {
          accession "AB000001",
          version 1
        }
      },
      int {
        from 0,
        to 89,
        strand both,
        id genbank {
          accession "AB000002",
          version 1
        }
      }
    },
    subject database "tests/db_search/arch_snucl",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, 1 }
      },
      {
        name "MaskAtHash",
        value boolean TRUE
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
