Blast4-request ::= {
  body queue-search {
    program "tblastx",
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
                  data str ">q_p_plus PemoMJNVA 10601-12200"
                }
              }
            },
            title "q_p_plus PemoMJNVA 10601-12200"
          },
          inst {
            repr raw,
            mol na,
            length 1600,
            seq-data ncbi2na '1EC0CCCCECDF4CF5AFB3A8ED3DB713CCFF3FF23ECC118E32
00A320C03A4008F77D47FCFD7D2080590098DAF7044FE5B14D3F82EC6224091D722B2823983678
4229A0E863E7ABA1AF1BD0588FFC70A54EDDD3047D0A209EE61EFF5BD0CF30112290FEA50A6194
2D80064C6BA2BF4D3688691A8BEB089220326ACE3EC16431EB20431733202262D5101655DB0CEB
C1B926077CF6407B01893F960B9E11E9C3989033DBC28B3818C1ADE186321FA0A0187EF4D4FF18
678007F7A48C40E11B98F355C5905C6B3EF7C1F306B78E983FE285D3A033F8C2A4F33E5849ACF7
B862FD30050813C98C0382EAD3D22D8262231868861263441D54D9DFED94E6411BFEC804220301
E49E16E961D6C09EFEE8D17020C41FDFEE38B265156FD00C342371338D8D134F283327598CD37D
D850C41A055C4ACE0F7DA2CF611CFFE88930C0BEF31041A8FCEF25B18D24E028308E6732D48254
DD659C13B138C060337DF17A680E03B7A5086911D63B42BBEA8200D8D230D0B4000EA2D94DA478
E2D9FA08BFA6092D0C4A80EEB82E94F285A7123EB3D962E832'H
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
                  data str ">q_p_minus PemoMJNVA 32001-33500"
                }
              }
            },
            title "q_p_minus PemoMJNVA 32001-33500"
          },
          inst {
            repr raw,
            mol na,
            length 1500,
            seq-data ncbi2na 'B06096DB60CDD59335BDD239B8CCF71A5F87F4374E2CF59C
FBFFDFFAF3A8E0C53EA28387244BDD0C8C3FDA2232FDF1E827B7855DF7AECE94A0C80F74D4D4FC
7B4048C4DA4C4F700BD0DBBC3B19F7FAE6AD036778E015B98618E0240DB13BD0030FF906466CC6
F8903FD43D4942EDD9360C45A692D4C2606D93FDB14AC93E1BD76642D03A97F2BA0131B80C4DDD
2A8186CDD08B0DFD8C3E6FC523F3A8C35533D9D35219C6818023A0CD6B7BDDBD37D52A5ED6BA67
8E18F758633D7C4ADFDDD72476B22D84032EDAD80BA6DB34D8DB932B3D3CF9EBAEBAEBAEB16516
9E8C5DDF83BEDC1252DA30DD8777B603F762DBDD4B68810FDDC09BF88D85DCEF7FDDE384F77063
FC5BF6F04864AFD2DA5B3DFB079B7E5B4DBC33CDF52C90422F69CE2AEC883F698DFBEEB802F33F
1CF258924CA6FFCC365489375E56509DD3B73159281B22D523DF3FDECD75DCA5369F68FF534C48
D7D27EFF8DB45527177779F5D4023236333DFC49F7D0296431CFF57A674F27583320DB4FF3D84E
'H
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
                  data str ">q_v_a LvMJNV 28501-29500"
                }
              }
            },
            title "q_v_a LvMJNV 28501-29500"
          },
          inst {
            repr raw,
            mol na,
            length 1000,
            seq-data ncbi2na '354308D17BFCCEC10AC22BD2EFB6FFD452EFCEFFFFCF2CFF
0C837F43F2CBB09C0BD7F8C3C3F42083E87FFFC31634CF1137CB7F3FF7FF3E490CED0BE09A4CD0
A13D34001E8021760825C4B8A40036319A87774D09F90DDE690A9B5844E0B14FEB0131309EBD62
0B860C4A3A67A9181BE26F7645B454A20822063D080A22285ED023D4B37B909BA365C60D88A2B7
85BA1280ACE4880CE21EBCE80C722193D9E14C3A174C831CC6B621F5C6075DC4011D881680C42C
6624DA6C3DCAD178067D321EA6E458860DB55DE627C50D881107FC4D09EED8082A81B6004E0C80
26922CEF5463E180EA55F8D4311024DA84D3A1A5C36341522215FC60C6B961'H
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
                  data str ">q_v_b LvMJNV 37801-38600"
                }
              }
            },
            title "q_v_b LvMJNV 37801-38600"
          },
          inst {
            repr raw,
            mol na,
            length 800,
            seq-data ncbi2na '0702CE1DB480302C0CE10812091B39B0433CF00D1034070C
09CC88B38088F3033240033040CF70380AAC30CF00013C0E33204073F4800408301FCF20010203
20300CCBFC333FF3FCDA40F47C4295731ECC51633403DFA50C50FCD04063F0C4FC3F63B67ADD00
DED2DDDFCC3734800412ED0FEACECEFC5234C8B238D404E3D00C4BE1ADF9030D9898C54B74EE76
DEFACC5CD8B0FE3A40CDCD0023E1C4429114F900343AF56BEE776C106B250BE714CC53A8F98C5C
CE8FCCCD8500EFF34CD10C885D9AF526F910D11E'H
          }
        }
      }
    },
    subject database "tests/db_search/xslice",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, 1 }
      },
      {
        name "LongestIntronLength",
        value integer 0
      },
      {
        name "MaskAtHash",
        value boolean FALSE
      },
      {
        name "SegFilteringWindow",
        value integer 12
      },
      {
        name "SegFilteringLocut",
        value real { 22, 10, -1 }
      },
      {
        name "SegFilteringHicut",
        value real { 25, 10, -1 }
      },
      {
        name "WordThreshold",
        value integer 13
      },
      {
        name "WindowSize",
        value integer 40
      },
      {
        name "QueryGeneticCode",
        value integer 1
      },
      {
        name "DbGeneticCode",
        value integer 1
      },
      {
        name "HitlistSize",
        value integer 500
      }
    }
  }
}
