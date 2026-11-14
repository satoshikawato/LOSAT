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
                  data str ">q1"
                }
              }
            },
            title "q1"
          },
          inst {
            repr raw,
            mol na,
            length 900,
            seq-data ncbi2na 'D7E64DED8F80C92585ADCEF1D3FEE6718A6A5D6FC5236B91
DC5AF6044C3319361EBA005D4621E221E5387CD47F6230DC859ED5DC9527379A57520342B4F010
98D2A84A378A9F41DCAAC624FE2E3D91A6D4B6FF8177D1DE0168E6974E0E493F291BEA5C534444
4444444444444444444444444444440444444444444444D6262CAAEC7FDD55DB2C9708BF056465
3A7E250D9E229E28FBF3D29120A2480A9E533FAA536D109DC166DE9B18E1C19DDA7F795A093716
D64BDBB116291B7AC237E816AF68520D3CBF6BCCAC889BCF25065FEDC39D52028B223FA4496F53
5239D335F681'H
          }
        }
      }
    },
    subject database "tests/db_search/dbmask/zh_k5",
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
    },
    program-options {
      {
        name "DbFilteringAlgorithmKey",
        value string "mydust"
      },
      {
        name "SubjectMaskingType",
        value integer 2
      }
    }
  }
}
