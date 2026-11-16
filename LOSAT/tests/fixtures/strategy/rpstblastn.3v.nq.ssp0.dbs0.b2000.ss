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
        }
      }
    },
    subject database "rps_mini3v",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, 1 }
      },
      {
        name "DbLength",
        value big-integer 0
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 0
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
Blast4-request ::= {
  body queue-search {
    program "tblastn",
    service "rpsblast",
    queries bioseq-set {
      seq-set {
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
        },
        seq {
          id {
            local str "Query_5"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">YP_009220483.1 bt strand=- flank=1"
                }
              }
            },
            title "YP_009220483.1 bt strand=- flank=1"
          },
          inst {
            repr raw,
            mol na,
            length 864,
            seq-data ncbi2na '941F74780F1A241A0C36D7610CF92E36B1E8AB1DE1AD5755
0CCFF9621EC136490FD7BFE6EC41F83B4D1AA0823CD7E717FE3267DCEF0C01E5C30FD44E3643B8
C008B709E9331649E80EEB50343E568792E78E42BDD0658DF3E8FA4CEF65350D0BF2E3AEDC9128
BBFFA62FC1E5DD5FABD7CA02C85A75F75B9EF2FF0CFFDDE39E3DD60B02CBE5C38C8351F632ED9F
F470763B667730BE90135FB425931AFE0FD4BF49A5AE7B667E4B673425B6824D0B362FA8A299A4
1A304E92E3B4E0F37C8B2D41D421B40B7EF5D17D7537EF2ED4EF627AAE3EFCE75682132F'H
          }
        },
        seq {
          id {
            local str "Query_6"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">YP_009220649.1 bt strand=+ flank=2"
                }
              }
            },
            title "YP_009220649.1 bt strand=+ flank=2"
          },
          inst {
            repr raw,
            mol na,
            length 712,
            seq-data ncbi2na '9384D15957DD6DB5115036D3B18732D0C8E69F7E7230C102
2B2F3F38C70C90BF0B989608376278DAAEFB7AEEDA306EB4648FDD585237FBAB1F33AC03307ADC
02EB4E06BDD633B4043E5C8413C4216D62B4D8B7B4F080109D21FC2380C4620B70E8511220B1CF
DD2B85001C24034BD3AE471D53C41AEA0E41D65C803FEF0A2811806B2DB51055D6B385ABF0CCFC
1DB45EFC7D8AD4C380E8E6F8FFCE9E7FFE09304BFBAAE23039833A17B11FD8F30A2E906567'H
          }
        }
      }
    },
    subject database "rps_mini3v",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, 1 }
      },
      {
        name "DbLength",
        value big-integer 560
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 0
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
      },
      {
        name "MatrixName",
        value string "BLOSUM62"
      },
      {
        name "GapOpeningCost",
        value integer 11
      },
      {
        name "GapExtensionCost",
        value integer 1
      }
    }
  }
}
Blast4-request ::= {
  body queue-search {
    program "tblastn",
    service "rpsblast",
    queries bioseq-set {
      seq-set {
        seq {
          id {
            local str "Query_7"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">YP_009220484.1 bt strand=- flank=2"
                }
              }
            },
            title "YP_009220484.1 bt strand=- flank=2"
          },
          inst {
            repr raw,
            mol na,
            length 718,
            seq-data ncbi2na '0FC14294C32224438E7AF37C35588E3119A2001F096609AA
74FF7F383F2A25F752BFDF7CF6DFFBF7EDAFAFFF4367EFBD0FEE0FEABED8FD8E99F7D7FCF5CDE4
243335F74CFF3816B61B1FFE3CB63E66EB9EF9871D8D8202F26D2586C877B367B7DF7D77D94240
7B71BEBB4877834C77DF5DF5DF7D358E0ACD92182FA6071F8E27238E3B23B6394DDFCFACB6F1F8
53BDF7D84DB7DF5D0FF3D22663F0FBE78C3862CFAE12C71D443999760924765568D4F728172360
'H
          }
        },
        seq {
          id {
            local str "Query_8"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">YP_009220489.1 bt strand=+ flank=0"
                }
              }
            },
            title "YP_009220489.1 bt strand=+ flank=0"
          },
          inst {
            repr raw,
            mol na,
            length 612,
            seq-data ncbi2na '0FC34E86FD74C40D863E3CD30CD8233A185F41A905C4285A
3A21E27ED7A23A6082828671B9BDDC9053A383E82C520815946E73179A7313F34A6489D65D709A
03537D2205B9971562F9C5E3916256B1A7E289EBE3B1782B937F7D173556CE18E8D36174DB44A4
8398A1871514533857D20054408250BF9873E441DB42D74B24E1A34E56E60FC3000DD00A802024
7287E969B043F4BCF0D4C011'H
          }
        },
        seq {
          id {
            local str "Query_9"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">YP_009220491.1 bt strand=- flank=0"
                }
              }
            },
            title "YP_009220491.1 bt strand=- flank=0"
          },
          inst {
            repr raw,
            mol na,
            length 627,
            seq-data ncbi2na '4B1221FFFAAE0533E301BB82FEEC60D5DF53B7DFE000302F
E0CFBF8F83A84E10AB75067E5C0B1DF973DFFD153B18B145743ADF11FFE331D76D3EC40E2D3D6D
15D3CCBE7B3081906511FB73EAFF7DAFDF536B6CBA8842F505D57EB539E26521763AB5AF2B7DAC
7CCAAF14C582C5F33158DB9962565F9BF5C7C44BF9EC292389BFD7B4600D8DE5F2EAD36EF6DD9D
F7DFFFF4B7E347A2407FE0183FA25330'H
          }
        },
        seq {
          id {
            local str "Query_10"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">YP_009220493.1 bt strand=+ flank=1"
                }
              }
            },
            title "YP_009220493.1 bt strand=+ flank=1"
          },
          inst {
            repr raw,
            mol na,
            length 884,
            seq-data ncbi2na '2CF85E8E8A88934B3781FAC6603DE888A6E002861645AC4A
04FB437FA1648E67F25E661847355B1D10500037DE99AECB6089CD3AF4F2280C48670A07E50000
0C4CFE83981A75660F7E8F72AA7B3858922AE376A56CCE9F50E8A47FA26833845E49DA3C47A82A
CC8D0FA70CCF801DF6200D78E16C838C3869FA0D53A1F1174E9FC5D6ED3385992FF3BA60EA81FD
3941FC5268BAE1B68FAABD4FC33E4DB32DEF054DD3A624EBEB78252B83D3DC1DDA639133710C74
CC8AF74009073B62E51857F50270A33E3A552208C8A1F473E13BFC5C8ABC505512E805F48E08E9
ED63'H
          }
        }
      }
    },
    subject database "rps_mini3v",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, 1 }
      },
      {
        name "DbLength",
        value big-integer 560
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 85680
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
      },
      {
        name "MatrixName",
        value string "BLOSUM62"
      },
      {
        name "GapOpeningCost",
        value integer 11
      },
      {
        name "GapExtensionCost",
        value integer 1
      }
    }
  }
}
Blast4-request ::= {
  body queue-search {
    program "tblastn",
    service "rpsblast",
    queries bioseq-set {
      seq-set {
        seq {
          id {
            local str "Query_11"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">YP_009220496.1 bt strand=- flank=1"
                }
              }
            },
            title "YP_009220496.1 bt strand=- flank=1"
          },
          inst {
            repr raw,
            mol na,
            length 845,
            seq-data ncbi2na '76D2C8154A9B830E39C5F7A74D2837E6840600B2B9C655FC
89FEBFDE86FB3ACC04FDF5FE2155FFFB74C363800C5319384FB3E78F661DA38D0BFDE513CF6B67
DD0DF6D302189DFE81693DD76C5A5DD7CDB94C776CD54E32B8A8F321F7DCFF54FCD93431A21EF4
F804C7CAAFBE65AB13FFB1D67FC0923D02529B37D8223802C516DD4B9E3FC580DB9303D58F36FC
BD7EDED3A760DF77B84F75F45BBFDFF346F6410C71B32ED74675D7DB553CDF7532CBF5D3E3FF5D
83FDD4C07B79828D96C5C6DCC756CE3AE1959F2636FD8DF3620C157D4CC5E600'H
          }
        },
        seq {
          id {
            local str "Query_12"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">YP_009220497.1 bt strand=+ flank=2"
                }
              }
            },
            title "YP_009220497.1 bt strand=+ flank=2"
          },
          inst {
            repr raw,
            mol na,
            length 909,
            seq-data ncbi2na '9A21AC8A0E9DD4200A334536B06D2E599E22274ECF75074C
E41346309E0C6F5E8D80EC8B6D098A000D29036D8C23C8D50105DE3374887B410ED1860E91E6FC
97FD82315F82A3AD0820CFA2305D437AE1F4F8E868001F061382D5C9C7AE637C24A48A87FBA68B
5304AC4021E9E0273CEAF9A4D359CC7B3DB8305F958D447E30E38E1025986C7304754CCE9EB1E4
03A63408F180F5A3D33240BA29E24C5B62A1B08E17D0CC0BA6A98473C105C30634D574B07A737A
624FD628648EE4A45FEFBE0FB043EFD32296138B2DECCE8600801E0989F45F81F1588C82A96390
39EDC623FB36FB6F40'H
          }
        },
        seq {
          id {
            local str "Query_13"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">YP_009220507.1 bt strand=- flank=2"
                }
              }
            },
            title "YP_009220507.1 bt strand=- flank=2"
          },
          inst {
            repr raw,
            mol na,
            length 706,
            seq-data ncbi2na '60609EA5BE3068FDFACF621CFEEBE471FBE55C456799EBEE
A8DAE7B265AE9249A79649777DBEF9E762CE37C919F3D375B6D72E08E51AF31E841FEFDCCF45CC
DFCD9BF74081F253F65926FB6E695192D3E5382D03E348EF3BEEBCFC3237166C653128A69DBE0E
F62C7E73F669F87BBA307FA35B4F4A0F49F7DDAEC0263B9579E12D534D4318392B5EA36BB76436
5307AEDA71EC039868C8953CD141558AC224962CC7D57BCFE9BB04ABC8C1526DCFDE93B0'H
          }
        }
      }
    },
    subject database "rps_mini3v",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, 1 }
      },
      {
        name "DbLength",
        value big-integer 560
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 99330
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
      },
      {
        name "MatrixName",
        value string "BLOSUM62"
      },
      {
        name "GapOpeningCost",
        value integer 11
      },
      {
        name "GapExtensionCost",
        value integer 1
      }
    }
  }
}
Blast4-request ::= {
  body queue-search {
    program "tblastn",
    service "rpsblast",
    queries bioseq-set {
      seq-set {
        seq {
          id {
            local str "Query_14"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">YP_009220508.1 bt strand=+ flank=0"
                }
              }
            },
            title "YP_009220508.1 bt strand=+ flank=0"
          },
          inst {
            repr raw,
            mol na,
            length 585,
            seq-data ncbi2na '96574BB363B74CF0491379661862490C77D737914512550E
3249266659A43E4D1B1379B596C1B5303D41DC6D76D643BF66F552A4C26C1E53965EE5267C72FE
160DC5A4628E2771D40EF8EE3D6BD5A550010B70B50000A3DC2600003401A016F42638A189625E
B0D861AB6A5A4B6003002C75EAF1CBF271841A48880C1844726556617F5C5045E08A804CE5C09C
150ACAB2ECC0'H
          }
        },
        seq {
          id {
            local str "Query_15"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">YP_009220514.1 bt strand=- flank=0"
                }
              }
            },
            title "YP_009220514.1 bt strand=- flank=0"
          },
          inst {
            repr raw,
            mol na,
            length 472,
            seq-data ncbi2na 'A02927578AFCF2F8C4EB0B179DFF84C99D3579CDF8B9FEFC
BB1FF073A23EA139F3955C73237CF26D0CF8FF3FF37FF5E6EF1EF9CA62B9EB2C7ADB6192C489E7
5EBBEF5D75EBD13CD77D06CF7DC95036777CDF75F07575512DD63A28F99C838DFAFBE3BE00C8AC
413CF8EFBF8DE54F9A920FF1EA4FDC96'H
          }
        },
        seq {
          id {
            local str "Query_16"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">YP_009220601.1 bt strand=+ flank=1"
                }
              }
            },
            title "YP_009220601.1 bt strand=+ flank=1"
          },
          inst {
            repr raw,
            mol na,
            length 1018,
            seq-data ncbi2na 'B2EAB3E03CE41D33CDD34DC18E9E96AD8BEB45A573BF96B0
3617178204FCD0420068A41011E5FFED011DDF8852B1A7BA46A4536C1939E882D320A39107BB7D
4D0A270C2EFF958ACB61BCC70F8E2A50F7F1617EE76D04A7E58C2AC26CCBCF9A6F23AC7DA34923
BD25535301FF16C450773ED03E9B403938FEC03211008A657F16D6FFAC38438C30ED37BBAAAE48
0ECE5967EDB87BCC000D0400800C2A82F2FB7E829A2367BA4224742601D71F1C3003D16F31A6A8
B33B9FD68DAD89447AC2F8F0E1CDD1C208C89D8634E4991F1FF4B9C36BA8ACE7508C0107DD863A
0D4ECB4E84ACF3CC9ACDB7DD1F958A6D846B22EBBD662E3829F1549485EBDC78EF9CB490'H
          }
        }
      }
    },
    subject database "rps_mini3v",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, 1 }
      },
      {
        name "DbLength",
        value big-integer 560
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 117970
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
      },
      {
        name "MatrixName",
        value string "BLOSUM62"
      },
      {
        name "GapOpeningCost",
        value integer 11
      },
      {
        name "GapExtensionCost",
        value integer 1
      }
    }
  }
}
