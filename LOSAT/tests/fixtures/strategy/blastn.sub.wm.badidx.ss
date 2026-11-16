Blast4-request ::= {
  body queue-search {
    program "blastn",
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
                  data str ">wm_q1 query from wm_mq00"
                }
              }
            },
            title "wm_q1 query from wm_mq00"
          },
          inst {
            repr raw,
            mol na,
            length 560,
            seq-data ncbi2na '83DF8F4E991EE6E973DBB4B9C4FDD399012A6C99416E98F9
D909B96FE13113E0218D2A9ABAFD738004DE7B63808C18D21E509DA6678F8E7A6C6F5F423E0A9B
C408CE271B820CD1694F365239FBE50CF8206E98DE9999F4D2B6C7813F7FBD4D868421F945BAC9
123CFFB8395B0E98CFA66F73D9603FB69796B2B82CB80BA20234DD8EF96F16251E903A63A5E0'H
          }
        }
      }
    },
    subject sequences {
      {
        id {
          local str "Subject_1"
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">wm_s1 subject from wm_mq00"
              }
            }
          },
          title "wm_s1 subject from wm_mq00"
        },
        inst {
          repr raw,
          mol na,
          length 600,
          seq-data ncbi2na '09B96FE13113E0218D2A9ABAFD738004DE7B63808C18D21E50
9DA6678F8E7A6C6F5F423E0A9BC408CE271B820CD1694F365239FBE50CF8206E98DE9999F4D2B6
C7813F7FBD4D868421F945BAC9123CFFB8395B0E98CFA66F73D9603FB69796B2B82CB80BA20234
DD8EF96F16251E903A63A5E06EF8F01B82DB6BFDA9431AD8801688034B164CCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCED'H
        }
      }
    },
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
        name "WindowMaskerDatabase",
        value string "tests/fasta/strategy/wm_badidx.obinary"
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
