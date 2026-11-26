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
            local str "Query_2"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">qn2 unrelated nucleotide (no hits expected)"
                }
              }
            },
            title "qn2 unrelated nucleotide (no hits expected)"
          },
          inst {
            repr raw,
            mol na,
            length 90,
            seq-data ncbi2na 'D39CDCE9103A13CFFE82B6F2ED62ADD32352FC197E8F80'H
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
    },
    program-options {
      {
        name "LCaseMask",
        value query-mask {
          locations {
            packed-int {
              {
                from 0,
                to 19,
                id local str "Query_1"
              },
              {
                from 0,
                to 19,
                id local str "Query_1"
              }
            }
          },
          frame plus1
        }
      }
    }
  }
}
