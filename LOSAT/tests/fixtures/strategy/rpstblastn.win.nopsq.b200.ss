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
                  data str ">w0 rnd"
                }
              }
            },
            title "w0 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 181,
            seq-data ncbi2na 'A72ED1E644B013CD913FF06AE26A4F0734523B8E6BF5E54A
5049287EB78ADA01B57C8F36B440DC9AC74E8929E440'H
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
                  data str ">w1 rnd"
                }
              }
            },
            title "w1 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 184,
            seq-data ncbi2na 'F4B612A7967DF1FC28BA5D6CEBB963FAFFD588A648156582
DC1FBBC878F86130101DEE72263614FBAF986E7AB0D9'H
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -8 }
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
            local str "Query_3"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w2 rnd"
                }
              }
            },
            title "w2 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 212,
            seq-data ncbi2na 'BAAB1DA9AB209C9D85617B7FFA5EB90BB79BC4C953E15E95
18CD38FB0F2D229BB0AA004559C6F98BD4A8EE88A49407A36A2D43D7E5'H
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -8 }
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 21252
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
                  data str ">w3 rnd"
                }
              }
            },
            title "w3 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 211,
            seq-data ncbi2na '7D1D62F3D54647D4D7668D8CD7401E4C61C840AA873EA0A6
BE0317FC1DE906F94D2CBA2A0B64D460489702D5EC5B030627C830FBF4'H
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -8 }
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 25653
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
            local str "Query_5"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w4 rnd"
                }
              }
            },
            title "w4 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 200,
            seq-data ncbi2na 'F53978547DF5B9D06EEF00140FCF1872B7645F501F8FC5BC
A17432943C81B5E0EB786DE445C82F5B756AED95603DB27FEAE8'H
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -8 }
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 25653
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
            local str "Query_6"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w5 rnd"
                }
              }
            },
            title "w5 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 178,
            seq-data ncbi2na '6C9E58D6B07E1FB921CF1FF57B8F8421436239282FCFB4EF
6924A49206C3A1B39B33822A26B3A2302C3261DEF0'H
          }
        },
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
                  data str ">w6 rnd"
                }
              }
            },
            title "w6 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 144,
            seq-data ncbi2na '7112C2F4BAC61A596301FAE6416ACD39ED7D680B206A95E1
D7F74849BDA20923E629D217'H
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -8 }
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 24288
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
            local str "Query_8"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w7 rnd"
                }
              }
            },
            title "w7 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 94,
            seq-data ncbi2na '9387417437716906B945D882CD2AC92B484FD05BB842C520
'H
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
                  data str ">w8 rnd"
                }
              }
            },
            title "w8 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 154,
            seq-data ncbi2na 'ABC9E9B711403896922653F7E264606B59E9AEA4E6D8EE77
B1C314A814CAAD970E95C67F195D60'H
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -8 }
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 20746
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
            local str "Query_10"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w9 rnd"
                }
              }
            },
            title "w9 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 228,
            seq-data ncbi2na '47B0861581A0FB837B3FC98101786C09F11D7D71124B5FF7
5062B0E097AE0A2C9D19B2C7C1ABED5A2D44691B3D21373AF8049D03E305EAE7DE'H
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -8 }
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 14170
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
                  data str ">w10 rnd"
                }
              }
            },
            title "w10 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 94,
            seq-data ncbi2na 'E10317C98218420A094451A1A5CBD7B04064A39D6294B430
'H
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
                  data str ">w11 rnd"
                }
              }
            },
            title "w11 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 162,
            seq-data ncbi2na '7C2B2B3921B5F66E2E91091ACE46A7943CB8EE60A52A2EA8
B242BD696696F448AE162567C31AA929C0'H
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -8 }
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 28000
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
            local str "Query_13"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w12 rnd"
                }
              }
            },
            title "w12 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 217,
            seq-data ncbi2na '39682B8CCCA5061C7DB22B0E4955BED354392138557BED6D
7846182B8B5AFEDDA5292C45C5BA261F89A0C22DD931D327CF7E51671F8F40'H
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -8 }
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 14170
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
                  data str ">w13 long"
                }
              }
            },
            title "w13 long"
          },
          inst {
            repr raw,
            mol na,
            length 250,
            seq-data ncbi2na 'BDB8BBD703397FDAD447C7043B0DB8EC26B3C046ADD47D26
EBE1309DCB912D746B820F3EDACB36FA05AC40DFBC32B3A00DBFAB937EEA6A3CA8C230D9D0C5A0
'H
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -8 }
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 26659
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
            local str "Query_15"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w13n nrun"
                }
              }
            },
            title "w13n nrun"
          },
          inst {
            repr raw,
            mol na,
            length 210,
            seq-data ncbi4na 'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
FFFFFF'H
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -8 }
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 31500
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
            local str "Query_16"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w14 rnd"
                }
              }
            },
            title "w14 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 63,
            seq-data ncbi2na '6FCAC4E389C1E32E0996F0008DC9183C'H
          }
        },
        seq {
          id {
            local str "Query_17"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w15 rnd"
                }
              }
            },
            title "w15 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 67,
            seq-data ncbi2na '75C1FA0C07FF94045AF51280EA88B67C30'H
          }
        },
        seq {
          id {
            local str "Query_18"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w16 rnd"
                }
              }
            },
            title "w16 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 108,
            seq-data ncbi2na '05AC800982F7AC33BE2B09D9A68EFC6756D1F9BFC8245D2E
8569C5'H
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -8 }
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 31500
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
            local str "Query_19"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w17 rnd"
                }
              }
            },
            title "w17 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 111,
            seq-data ncbi2na 'DA4561EB420B628E2AB25E47C5811D6D6ADDEB9D1255A5F6
5A3DCFF4'H
          }
        },
        seq {
          id {
            local str "Query_20"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w18 rnd"
                }
              }
            },
            title "w18 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 153,
            seq-data ncbi2na '91A0CA506C7F153C963C15B9A23082758F9B95969AC5573B
C565D5F8840ADA326909DF66A62F80'H
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -8 }
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 11760
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
            local str "Query_21"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w19 rnd"
                }
              }
            },
            title "w19 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 203,
            seq-data ncbi2na 'C14A83A3FB058E62F755155848A77F299820D8E405048E71
7778654983E9DFF01157E588383FF7ACDC729AB3215A925431AC5C'H
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -8 }
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 13702
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
            local str "Query_22"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w20 rnd"
                }
              }
            },
            title "w20 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 192,
            seq-data ncbi2na 'C179FD5705EFB633DC9C2C6D3B32F5E82C28DD37B72FF58E
9D1B434C6A8A95F84BE0AFA82638623E3DA116F1CA17971D'H
          }
        },
        seq {
          id {
            local str "Query_23"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w21 rnd"
                }
              }
            },
            title "w21 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 199,
            seq-data ncbi2na '13B141DFE0120DF04CE60C429FCB02CCC7DF42A60BAFC4D5
9174EDF71A3297DDE80E6FC23FA27D495AF47155EA3C4532BF7C'H
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -8 }
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 24144
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
            local str "Query_24"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w22 rnd"
                }
              }
            },
            title "w22 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 75,
            seq-data ncbi2na '51C148504A8E0BF2CBCA695D7352A14C309B70'H
          }
        },
        seq {
          id {
            local str "Query_25"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w23 rnd"
                }
              }
            },
            title "w23 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 156,
            seq-data ncbi2na '911B57DF8A84F7B5D6AFDE81DD1067506D9CBB783B895FA7
0BA71A87E2753B91D77E50B05CD069'H
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -8 }
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 23276
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
            local str "Query_26"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w24 rnd"
                }
              }
            },
            title "w24 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 74,
            seq-data ncbi2na '40D33E86E6CA6B5DACDB9F2C1407F3CE5C6E10'H
          }
        },
        seq {
          id {
            local str "Query_27"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w25 rnd"
                }
              }
            },
            title "w25 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 129,
            seq-data ncbi2na '822E3602A34B62E11E7B0681B99828E9F819DC3377B549A3
108924D9263BD49E40'H
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -8 }
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 14000
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
            local str "Query_28"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w26 empty"
                }
              }
            },
            title "w26 empty"
          },
          inst {
            repr virtual,
            mol na,
            length 0
          }
        },
        seq {
          id {
            local str "Query_29"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w27 rnd"
                }
              }
            },
            title "w27 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 122,
            seq-data ncbi2na '65F9CFA38F1C2595BA54FA5FE0E28C840BEB180E548C2FCF
58FDA6A1779230'H
          }
        },
        seq {
          id {
            local str "Query_30"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w28 rnd"
                }
              }
            },
            title "w28 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 198,
            seq-data ncbi2na '5950E935670C3A34C0A2499F3169797631167C78DC799CA7
8E4F6243618AB5016A32E595466C656C4CC5B75BEA9792677BB0'H
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -8 }
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 13440
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
            local str "Query_31"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w29 rnd"
                }
              }
            },
            title "w29 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 164,
            seq-data ncbi2na 'FDEA392744F6A9D128E8C8E2DCE7978DA97FFB513EE93689
934B37E79F4C96D4A0AB26AB9FC02F9A3B'H
          }
        },
        seq {
          id {
            local str "Query_32"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w30 rnd"
                }
              }
            },
            title "w30 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 208,
            seq-data ncbi2na '833B05CA2C746D3B7BFFDB08ADEA668A85537DFA781D463A
FDDEE642116663ACCCB0FCC25F6736463BBDAF281DEDAB741CC522AB'H
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -8 }
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 0
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
            local str "Query_33"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w31 rnd"
                }
              }
            },
            title "w31 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 117,
            seq-data ncbi2na '3D18CAAB29FA491CE2E6CB5CA50EACA78C512AA4BD7B0911
64F7ECFF9300'H
          }
        },
        seq {
          id {
            local str "Query_34"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w32 rnd"
                }
              }
            },
            title "w32 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 65,
            seq-data ncbi2na 'BF190C964451709FA6AA4E133D8D354880'H
          }
        },
        seq {
          id {
            local str "Query_35"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w33 rnd"
                }
              }
            },
            title "w33 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 154,
            seq-data ncbi2na '90255CB45102946BD060CEAD861F04DF349035061A78F93C
18036E7A22B6FB1F3AC6DFF1825BB0'H
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -8 }
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 18833
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
            local str "Query_36"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w34 rnd"
                }
              }
            },
            title "w34 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 156,
            seq-data ncbi2na '273C082D4796DA1285D2A321F9C6DBA18E2090B76D680A1D
34B5D9E10D392B7DAC622E99918225'H
          }
        },
        seq {
          id {
            local str "Query_37"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w35 rnd"
                }
              }
            },
            title "w35 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 92,
            seq-data ncbi2na '1D944DDF67ECEC95BC8FECF52D11491BD3967F9C012194'H
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -8 }
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 14148
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
            local str "Query_38"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w36 rnd"
                }
              }
            },
            title "w36 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 127,
            seq-data ncbi2na '8211E5D04FC696FDE9CDABEBEFB324BFE85667A2EED89A53
B943B593D9A1CD98'H
          }
        },
        seq {
          id {
            local str "Query_39"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w37 rnd"
                }
              }
            },
            title "w37 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 222,
            seq-data ncbi2na '44D9E39B595EE0FEC81FEB9C9006E3A223A717899E746247
7160D529951EF78B3353801EA088227DDEF23922E219BD80F055D5CB7E03F090'H
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -8 }
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 18432
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
            local str "Query_40"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w38 rnd"
                }
              }
            },
            title "w38 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 235,
            seq-data ncbi2na '5628FCCD5218C27288FF1432F722F143EB5C2949226BAAAF
57F0291D8D74E5247D1DAAD0A3A49676144BDFB58A5AB41BC3157FFA928ED0CE01C290'H
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -8 }
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 13905
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
            local str "Query_41"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w39 rnd"
                }
              }
            },
            title "w39 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 147,
            seq-data ncbi2na '8F36C74E6EED7F8FC581CC5A15503A4225EF25BA6A1E3793
94A0519F06E6304CED76AFB2C0'H
          }
        },
        seq {
          id {
            local str "Query_42"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w40 rnd"
                }
              }
            },
            title "w40 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 235,
            seq-data ncbi2na '720B58A064A4AD71A7E3745C10F78B973AD3D5B98FBA1753
A7729E2EC66A37E4C09FDCDF6B4FFB8ED49E8E9F743EACB2D949FA20AB82119A924EE0'H
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -8 }
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 29000
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
            local str "Query_43"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w41 rnd"
                }
              }
            },
            title "w41 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 214,
            seq-data ncbi2na '13474A9474E3640328EB3C0509BC6A4512A5E95619C40516
FE3AFE98A5ECAD6A3538923139C5BB8135B75D6C7F01717BFADA324050D0'H
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -8 }
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 16896
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
            local str "Query_44"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w42 rnd"
                }
              }
            },
            title "w42 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 94,
            seq-data ncbi2na '2A9B0F7D79A0C485F45FA766394DC09A81CB56425B49BC60
'H
          }
        },
        seq {
          id {
            local str "Query_45"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w43 rnd"
                }
              }
            },
            title "w43 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 155,
            seq-data ncbi2na '9A33BD7157657907ED4294A40D09A28DFDDB213FABC7FBB7
22060862DCDEFE1042A904B3126828'H
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -8 }
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 26156
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
            local str "Query_46"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w44 rnd"
                }
              }
            },
            title "w44 rnd"
          },
          inst {
            repr raw,
            mol na,
            length 108,
            seq-data ncbi2na '8428387BB1AAE28734A69F7C035608DE68F554F908FA7366
6DFC00'H
          }
        },
        seq {
          id {
            local str "Query_47"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">w45 hit 0 95-133"
                }
              }
            },
            title "w45 hit 0 95-133"
          },
          inst {
            repr raw,
            mol na,
            length 127,
            seq-data ncbi2na '2CD8E0AD2F7F1A372DFAE0D9F24842A4268CB8D924723AC7
D614123BFDD708A4'H
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -8 }
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 14170
      }
    }
  }
}
