Blast4-request ::= {
  body queue-search {
    program "blastp",
    service "plain",
    queries seq-loc-list {
      int {
        from 0,
        to 139,
        strand unknown,
        id swissprot {
          name "TST1_HUMAN",
          accession "P99901",
          release "reviewed",
          version 1
        }
      },
      int {
        from 0,
        to 44,
        strand unknown,
        id swissprot {
          name "TST2_HUMAN",
          accession "P99902",
          release "reviewed",
          version 1
        }
      }
    },
    subject database "tests/db_search/arch_sprot",
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
