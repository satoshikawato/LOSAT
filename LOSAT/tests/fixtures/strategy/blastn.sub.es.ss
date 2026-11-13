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
                  data str ">e empty"
                }
              }
            },
            title "e empty"
          },
          inst {
            repr virtual,
            mol na,
            length 0
          }
        },
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
                  data str ">s two"
                }
              }
            },
            title "s two"
          },
          inst {
            repr raw,
            mol na,
            length 2,
            seq-data ncbi2na '10'H
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
                data str ">sn1 subject nucleotide alpha"
              }
            }
          },
          title "sn1 subject nucleotide alpha"
        },
        inst {
          repr raw,
          mol na,
          length 488,
          seq-data ncbi2na 'EFE7BED928FD49E4F6C028CEE8E811D6E72090FFE5B9F6C0B8
E86158031A8EA54E9E4382A7A881CF9E95983309E80B5F382FE2FF5D633515F90FCE4FE7970A33
A4054A3ED821A31BE4CAD243AB7BA293093742ACEEEC1BAE20D12343775A4F9017053867079880
FC237D4413AEAD3A88124DE2A41BAE1A7F82C9'H
        }
      },
      {
        id {
          local str "Subject_2"
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">sn2 subject nucleotide beta reverse"
              }
            }
          },
          title "sn2 subject nucleotide beta reverse"
        },
        inst {
          repr raw,
          mol na,
          length 477,
          seq-data ncbi2na 'EFADCF303F2746F3BC7C95D279EF754E1514ECBA28DF0BF744
BC9B4E9C3FD4F94888B0FEE0D117B111F11FA8C7CE4D4685BFDE146CF61F550D7A8B94DC7CA8D8
9025103E426E8CDA380150070131A1FD49D23C9A9A90CB9D4E57D3C494FA54C7B3FEA36D4D1F47
095A403F9FF2468BBD4D8C6FCD071831C0'H
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
