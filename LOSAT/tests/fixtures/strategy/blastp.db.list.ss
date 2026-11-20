Blast4-request ::= {
  body queue-search {
    program "blastp",
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
                  data str ">P69905.2 RecName: Full=Hemoglobin subunit alpha;
 AltName: Full=Alpha-globin; AltName: Full=Hemoglobin alpha chain; Contains:
 RecName: Full=Hemopressin [Homo sapiens] >P69906.2 RecName: Full=Hemoglobin
 subunit alpha; AltName: Full=Alpha-globin; AltName: Full=Hemoglobin alpha
 chain; Contains: RecName: Full=Hemopressin [Pan paniscus] >P69907.2 RecName:
 Full=Hemoglobin subunit alpha; AltName: Full=Alpha-globin; AltName:
 Full=Hemoglobin alpha chain; Contains: RecName: Full=Hemopressin [Pan
 troglodytes]"
                }
              }
            },
            title "P69905.2 RecName: Full=Hemoglobin subunit alpha; AltName:
 Full=Alpha-globin; AltName: Full=Hemoglobin alpha chain; Contains: RecName:
 Full=Hemopressin [Homo sapiens] >P69906.2 RecName: Full=Hemoglobin subunit
 alpha; AltName: Full=Alpha-globin; AltName: Full=Hemoglobin alpha chain;
 Contains: RecName: Full=Hemopressin [Pan paniscus] >P69907.2 RecName:
 Full=Hemoglobin subunit alpha; AltName: Full=Alpha-globin; AltName:
 Full=Hemoglobin alpha chain; Contains: RecName: Full=Hemopressin [Pan
 troglodytes]"
          },
          inst {
            repr raw,
            mol aa,
            length 142,
            seq-data ncbieaa "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFD
LSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASL
DKFLASVSTVLTSKYR"
          }
        }
      }
    },
    subject database "tests/db_search/taxfilter/tfp",
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
    },
    program-options {
      {
        name "TaxidList",
        value integer-list {
        }
      }
    }
  }
}
