qx1	sx1	89.744	39	4	0	4	120	40	156	1.25e-25	95.9
qx1	sx1	79.487	39	8	0	119	3	155	39	1.36e-18	72.5
qx1	sx1	58.333	36	15	0	120	13	156	49	1.25e-17	69.3
qx1	sx1	62.162	37	14	0	9	119	45	155	2.18e-12	51.9
qx1	sx1	60.000	35	14	0	17	121	53	157	1.07e-11	49.6
qx1	sx1	55.000	40	18	0	121	2	157	38	5.23e-11	47.3
qx1	sx1	25.000	16	12	0	66	113	150	103	4.1	11.1
qx1	sx2	88.000	25	3	0	28	102	190	116	7.18e-15	60.2
qx1	sx2	63.333	30	11	0	101	12	117	206	1.25e-09	42.8
qx1	sx2	48.148	27	14	0	32	112	186	106	4.83e-06	30.9
qx1	sx2	52.174	23	11	0	30	98	188	120	9.12e-06	29.9
qx1	sx2	48.000	25	13	0	105	31	113	187	9.12e-06	29.9
qx1	sx2	56.522	23	10	0	100	32	118	186	3.00e-04	24.9
qx2	sx2	100.000	36	0	0	120	13	175	68	1.64e-24	92.3
qx2	sx2	64.865	37	13	0	14	124	69	179	1.47e-11	49.2
qx2	sx1	85.000	20	3	0	120	61	79	138	3.03e-13	48.7
qx2	sx1	81.818	11	2	0	45	13	142	174	3.03e-13	21.2
qx2	sx1	70.588	17	5	0	74	124	125	75	9.12e-06	29.9
              label str "DefLine",
                  data str ">qx2 minus strand gene"
                }
              }
            },
            title "qx2 minus strand gene"
          },
          inst {
            repr raw,
            mol na,
            length 129,
            seq-data ncbi2na '6D8BF231F4D206E6674EFA639CECA3824338F54C5F4F1E53
EB5B560FF6AD058140'H
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
                  data str ">qx3 no hits"
                }
              }
            },
            title "qx3 no hits"
          },
          inst {
            repr raw,
            mol na,
            length 90,
            seq-data ncbi2na '9050B43EF2DFF11020771ACE787B7D4758BA590B1F7570'H
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
                data str ">sx1 plus strand gene"
              }
            }
          },
          title "sx1 plus strand gene"
        },
        inst {
          repr raw,
          mol na,
          length 218,
          seq-data ncbi2na 'DB2A39E2E5F72DA449BFC0A5BBD96D29B0AD3A385400D221A9
53A6BCE029EA00F3E7A5240CEDFA02DF0E0B38BFDB584D4D614720721440'H
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
                data str ">sx2 minus strand gene"
              }
            }
          },
          title "sx2 minus strand gene"
        },
        inst {
          repr raw,
          mol na,
          length 208,
          seq-data ncbi2na '7687CAC5E73F04AD4F8F13A3A3968F8074CFD3C21389E13BD9
CC8BA2CE012C1DD511FD30494FA16DDF37CAB6D4C17F688DCC4C4E'H
        }
      },
      {
        id {
          local str "Subject_3"
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">sx3 no hits"
              }
            }
          },
          title "sx3 no hits"
        },
        inst {
          repr raw,
          mol na,
          length 150,
          seq-data ncbi2na '1DFA657DAE6FE51B93E11C7026237415CE1BD917E1ABE83A41
A76D24B8EC9928D4C9E7097300'H
        }
      }
    },
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
