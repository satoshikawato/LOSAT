Blast4-request ::= {
  body queue-search {
    program "tblastn",
    service "rpsblast",
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
                  data str ">stopq"
                }
              }
            },
            title "stopq"
          },
          inst {
            repr raw,
            mol na,
            length 600,
            seq-data ncbi2na 'E32E38CB2C30CB8CB2C30E38C30CB2CB0CB8E32CB2CB0CB0
CB0C30C30CB0E38E32C32C32E38C32C38CB8E32C32E38C32E30C32C38C38E30C32C30C30E38CB8
E30CB2C32E38C38E30E30C30CB0C30E38C32CB0E32C32C38CB8C38E30C38CB8C30C38CB2CB0C30
C38E30CB8E38C30C30CB2C30CB8E32E30CB8CB0E30CB0E30CB0CB0C32CB2C32E38CB0E30C30C32
E30CB2E38E32E38CB0'H
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
                  data str ">YP_009220476.1 bt strand=+ flank=0"
                }
              }
            },
            title "YP_009220476.1 bt strand=+ flank=0"
          },
          inst {
            repr raw,
            mol na,
            length 637,
            seq-data ncbi2na '6E83933BAAEFC6794FF26ABC17C31DBEF3C932FB870CF81F
04001C8400230A39F356E8B488C3417C53C3AEC998A814FF07FB03A86F10118071A80BF1BA7AD2
A84261262F2E08002A8CF37E17675DAE3AA1445741B810A5837020FA0F31958B330C1102663346
BA1DAFD9ED97A7C0E616B620FF79430BFBFEC7F09412EF0008168D047FAFF68985FE38A3323847
0484737FE2CE2BF2D5A418FEF42BEB8821D9C0'H
          }
        },
        seq {
          id {
            local str "Query_3"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">YP_009220478.1 bt strand=- flank=0"
                }
              }
            },
            title "YP_009220478.1 bt strand=- flank=0"
          },
          inst {
            repr raw,
            mol na,
            length 905,
            seq-data ncbi2na 'C0A9D3F34BC9D7F0B7D2D327FCFF3DB0E5D8D460817CCD17
D9F6BF047F528F5328DCAB32D6067633FFC61C4A37A7A3B09638A4D80DC7D75521395B7701A8C0
26CD48D77DF5DB5D75377D54737003DB43B6F36D80DCEDCED34D36DF454271C8A87E9BD8B68AE8
AA8BE8B2BB29AA8BE8BE9BA9AE9AA82405775281AC52F5D5515BBC7D4BBC1E3AF9D81D821F7D75
DFE5B7DC69D7DF6D9858EA5CF75A7465A58A1E6575BB4D75FF7F77F35FFEBFFFBF82CB6EBFCF2B
7A13033865451C458609310328D0E304108134F3CC60FFFE4743CD994800C027FEC8A49CFD5651
19675BE1A7BDF4C0'H
          }
        },
        seq {
          id {
            local str "Query_4"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">YP_009220482.1 bt strand=+ flank=1"
                }
              }
            },
            title "YP_009220482.1 bt strand=+ flank=1"
          },
          inst {
            repr raw,
            mol na,
            length 624,
            seq-data ncbi2na '9B02790B3843BA0C036C30414068E07A5CEBACA2EF7D3DC9
4D1E6F24BA773135A4E5DC04220B32C4DD13C1892FC59C249D52D1B1774004AC27C28CE4DE111D
23F4E9F1BA1BD10DB10015F09154EE478609EBEA64FDE7C725D668A4C5860C2B1A8AA62F6A5E00
13D5C9A018D4B616D4233B4814CCD1133FD0820C686C8EA468390DECFA208168100F6CEAA1AC01
8BDA21730FED630B73373EDA6F9D4E'H
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, 1 }
      },
      {
        name "QueryGeneticCode",
        value integer 1
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
        name "UngappedMode",
        value boolean FALSE
      },
      {
        name "HitlistSize",
        value integer 500
      },
      {
        name "CompositionBasedStats",
        value integer 0
      },
      {
        name "SmithWatermanMode",
        value boolean FALSE
      }
    }
  }
}
