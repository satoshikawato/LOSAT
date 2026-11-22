Blast4-request ::= {
  body queue-search {
    program "blastp",
    service "plain",
    queries bioseq-set {
      seq-set {
        seq {
          id {
            local str "dup2"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">lcl|dup2 first"
                }
              }
            },
            title "first"
          },
          inst {
            repr raw,
            mol aa,
            length 150,
            seq-data ncbieaa "MASLGLLSTLRGILQGKRNADIENLISKLHYRVTTTTVFCFCLLVTTN
TLVGDPISCVYDSGHPKLNDKVVNTYCWIHTTFTLPNRPYRDMSAISSPGIGTHNLSDPVTYHAYYQWVPFMLFIQGL
LFYAPHYIWKTWEGGRMQSITADL"
          }
        }
      }
    },
    subject sequences {
      {
        id {
          ddbj {
            accession "BDT61857",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61857.1 innexin [Penaeus monodon majanivirus A]"
              }
            }
          },
          title "innexin [Penaeus monodon majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 378,
          seq-data ncbieaa "MASLGLLSTLRGILQGKRNADIENLISKLHYRVTTTTVFCFCLLVTTNTL
VGDPISCVYDSGHPKLNDKVVNTYCWIHTTFTLPNRPYRDMSAISSPGIGTHNLSDPVTYHAYYQWVPFMLFIQGLLF
YAPHYIWKTWEGGRMQSITADLSVPILNKEDQRKRVNQLADYLQNSLWNNNTYALKYLFCQCLSLVNCIVNMVLMNKF
LNGKFLTYGADVFSYSEAVDQESRTDPMIEVFPRVTKCSMSMFGSSGTLQKHDFLCVLAVNVINEKIYLAMWLWFILL
ISLSILVFFFHLAVFAVPVIRRLYLVYLLGSGDKARIVARTFNQMGDVFLLYLLKKNVECQSLNLLLTELCREPITGD
DESKKYQKYLRRESEV"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61858",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61858.1 wsv026-like protein [Penaeus monodon
 majanivirus A]"
              }
            }
          },
          title "wsv026-like protein [Penaeus monodon majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 1472,
          seq-data ncbieaa "MIEKDRIKWQKDSLPLLFPSEETAKAIGSNTFAVPSFEVYESKALLEVGR
CESLTRGGMDDCWVDGYVQPRFLLRPMSLITLQGEAVRLFSVQYYKHRGNLGQGDAKSKNAYGGGFIIGDGTGVGKSR
EIAGMIVTAILVEQYLYRKSESHKPPLVIWLTCSELLFANCKRAFAEVLTLANASKYSLRSMNDNGLTTIDLEGNDLF
IHFTTLKNFLADTNDTCDYPLPQPTVLFLTYNGLMANFEDLMEIFDKGIILPTAVFCDEFHKTKNISDKMKWVIQSRS
ERYDGDDSDHNSPSLFVAMRNTFCRNKRNKLQLTVADSVKLFVDHLKNTTFFVMSSATPFQNNQDLHMIDHIIRNIAP
RYHLLDQYNGTPTGMNSSEYSTLFLESIIKLLYNNGIYVSRTISMKGIRCAIVQKPISAANMYMINEISSYLAEMKCL
AKTAHSDVKVLGKKIDQIIKSKKWESPSALMSRLEEFGEAVNTGKCGEVAIRTGYRLVFAEWNSAWQYLYDGRRRKTT
YSNDDDDDDDDDDDDDIFSKRKKIKMDPTIQESLSSSSPPPPPPHQTEKGDTIDSDIISMSKVFSTRSCRLLSGQLKV
NLAACCVAINKNILLSLKNKSVLKYIETVRKQKEPQKMVITVDQTGDSYYSHITDELIKIPLLESGAKKKKKGKKKKN
AISVNRLSVFDKAPASHFIGNAYKTLLKIILLDTTYRITSLEKERVCFLLLPRLPPPAPLALLEENFLDSIEVGVGEG
RYIEVTKRKYSGRLMSEGTLEIIPINTSAKNTTRDFQLFKHSKSVDVAIVGRSGNTGFDFHDSKDNITRARRIHIIAD
LPHDAISFLQSIGRTHRNNQRSLPRYCLFLTNIPTEKRFVDTLESRVRDSKAGTYGDRYCQNTISLKKEGKGVDGKII
ASTTDGDSVFAGSGTDDNDDTAATAIAAVTNSVTSTITTTANTTSVIETKRPTNKDITTCTDSSHDDLSKQNKCIDKI
DFLENDFVMKVLGILLQLLAGEYNPLDIFLQIAKSGYRDGSLFIPLEGLGDSKIYLSIVEFAFECAAALLYESKIKCA
SQSFRISKVARATVSKIDRDKKRKTLDEILERLCYFYDSSIMDIGLRGNIVEEKREIKNLVSVAIQFWKTFSPLKMTK
FSHLLVAVKDIQRGNNISVNSISKSSSNSSGNNYAASRQRKMSIRDIFRRSAVLVFGYTVRDNIQLPPPTDVLTVSIH
DGTMKLIMNGVVSNIPMVAACNLLDALREKTPSLLLDLVNEITPYKKENSDDVESCLLTLSKLLVNKGLTYKQFRNYF
MTIPLSEQDMLQSLFSKVQEWLSTEERFSRMCMMRRSNCISAEYTNVIACPSDVSRLYARVDDNDVTNITCDSNFAID
CRLSLCNDNEGPARRSTLTNTTTTTTTTTTTTSTSTGVTTAARGDINIRLTPENTIFINNLFNTLLLTNARSDMVKFK
FWNSEIVGLPSVALVRYQ"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61861",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61861.1 wsv206-like protein [Penaeus monodon
 majanivirus A]"
              }
            }
          },
          title "wsv206-like protein [Penaeus monodon majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 234,
          seq-data ncbieaa "MTEYSSRWENRIFREDVKNYHIIEGDIFADRFRKKNICICQQCNCITLKP
HGLSADIVKHFGEYANSYGRRRAKPNTINTATLETRAEPGTVDFCEGTPCVANLFAQFFYGSSSNVQYLSRKYLDPHI
RKGVDNDTAEKRLMYFRSCLNNLMFDLRNKYTYIDTVVFPHRIGCGMAGGSWADYQNVLLQFADRLSEDNITVFEVRN
DAAIEATKRYHSRPVEGRWREREGITTD"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61864",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61864.1 hypothetical protein [Penaeus monodon
 majanivirus A]"
              }
            }
          },
          title "hypothetical protein [Penaeus monodon majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 154,
          seq-data ncbieaa "MPSLLAGLEDPSTPSEGTTVEQPLMSSGGAEGEEVERGEKMDPDNLSSPS
VENNEADQRSQPQPQSKLSPMPELKPLPRTRYQSPQISQCQSLQMTQTLSLKAIHSPTIRVFLKEDAQPIILQDSLPW
LDGRILEPEYDVFVRLTRRDNRPNPM"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61865",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61865.1 hypothetical protein [Penaeus monodon
 majanivirus A]"
              }
            }
          },
          title "hypothetical protein [Penaeus monodon majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 215,
          seq-data ncbieaa "MEEYSSESTESSMGSPVVTINRKRRRMKNPEKEYSVDKRKKCQMSKNDPR
RNEVGRSNSSDDTHSHNDDDDDDDDDDQGPYSRLSEDEATSPADVHHIDSDETDIEIYDDNNIENNNSNDYDDDNNDN
NNNNNNNDDGNIVHYRFRPFTAEAALAQQILSNTAGNICGILGKLEYHLSELCPKIRAITTYSNHFSKELEKLKDNME
KDLENSSST"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61866",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61866.1 wsv220-like protein [Penaeus monodon
 majanivirus A]"
              }
            }
          },
          title "wsv220-like protein [Penaeus monodon majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 731,
          seq-data ncbieaa "MSNKMTILYSELMSARENSIAAFEEAVKNISIYLLEEAERGSWGDDQNKL
KDLYDGKSEADGLEEDTENKNLGLYVPAGIDMRAWRAGDALGDYKNAYAARLIVNITFHHNKIAEILYTLIAELLLLL
EDNINDDGKTQLQEYGRLKPASVNETVKSLENVIREKNIEVDLKRFREIVLRLENDSRNFDRESRLSDWLLDNIQERY
PAGGGTTTTTTTSNNEYLCTIDDTTPLSTDTILSTLPSARRERPVRNIVGESSSAPPDRPWEDEREQTGYSIFSFRSV
WMSEYGDYPINLVNAIIERLLLRYVVPLRDVFHVCFHLKAHLTCARNVNATWYENATFAYGLPPVYSMRDTWLELENL
LNVYACVAKLFLNMYDLLSSSSHGFHQSDLTRTKEAYINTIELLECMPMYLFDSKWMDENFYSWPYTRRGVRQLPVRN
YLSENYLLRTVLVIPPNGIHPITKRKNNSGILMRLKVKAVEYITHLRTDMREIFDDGFVTIVKDNLQKTTNTTPRSDS
RNINGMNDFEARILCDSTFENNPYVVHPLYHLKIAHRPSNWSGKGLASIEEWDARYHLASESGTLSSSVTNKNNNNNN
NNNNNNNNNRSSWFWRSLPLFKPSLLKDSSSDSLPYKFWEAALRHIHSHQEFDRPMMNYSKNTEALEIMIKMVEVTEA
STEKKKKSDIEDILVQDKQLESYEKQEERRLRNVERSIFLTILCKNMLTTLMNVGGI"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61867",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61867.1 wsv216-like protein [Penaeus monodon
 majanivirus A]"
              }
            }
          },
          title "wsv216-like protein [Penaeus monodon majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 918,
          seq-data ncbieaa "MEGYEIFSFISVRLVVIVLLLLSLLLYPIIRMISDYKKQTDEIVSPKVDS
FVNTEISHIVKDIDDKEHANYIASKIRNRVYDCCNVYDNDDDDNERSDNVEDDYSSSKDTLLSIKQVRIRMSEKGENN
RIVSLSALYFLKASILRGRLPVSVRITFYRKDRPLCYMDTIIQSLGKPIVLKEISKCLIYGEIYRVRIESRRTVFFDD
FVLIDSPFPEMVPVLICQDAILLKCNYKTQTATPFKAVEYNSKLFSGVYLTDFPKLLPTAGVQTGRDNILEPIPGKTS
GHCLLYAKDPIRFDGSLSFFSDESDGTTNADLTLIKSSANPIEYYVLLKRDSPQPTRRTHANIKGWVINYRSSEGAEY
KQIFSLSASIEKYYPRSQIEIRPNEKLRDITPEDTVTVWNTDEKIGVDKPETGTATAAETITVNGSKLFSRNAKDKGI
SFPKDSSFFQRDPHIVLDRETGECVMVTLGFHEDVALNPRWASSPCHFTSRDKLLVINPQRKHHLSPFVTSVRIYRGI
PGGRSFTSQQGEALCSSVIVDRGDNLPSLHLPTISLPEERQRFCVHIQHKEPGTSRYSLIEHFPRSLPLETPLDVSSD
RIYRIYTRIPSSGIASPTTANQTYANVGVINETDRPYSDLANILDDDAEGDGDDRAKHRERRTCTIFAIGLVNTYSDG
TSIPDTIVKAIRSSLPEVVSITDINGTEVVRSKNIIVNNSSVLIILKNIQIDGSSQRWISSTPFGLGPCKTSGRCGIL
GNWVPLYTVDHVLPICKVFSYRYADDSTARGFTVYLSNTCANGRSVSESGISEIVALYLSPEIQGKRIVLNAVAINCN
ANGKDSNNSTDSKRTKIFSIPSKESIGGGAEGATSIIILDRVVYALMAVVITKDVERERRHFARNGLFTARSVTLSEE
FVRRHRNVMN"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61868",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61868.1 hypothetical protein [Penaeus monodon
 majanivirus A]"
              }
            }
          },
          title "hypothetical protein [Penaeus monodon majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 76,
          seq-data ncbieaa "MFKEKIRINKQSTVILFSFYSLKIHASTYLYKIILKETIIFLLTIIEDRG
VHIDPIILLPFLSSFLETSWVTEIMI"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61874",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61874.1 hypothetical protein [Penaeus monodon
 majanivirus A]"
              }
            }
          },
          title "hypothetical protein [Penaeus monodon majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 741,
          seq-data ncbieaa "MSDKQLLEVERVLLGASSTTDPQKGLNVITKGAIRQLLSDRERKLKEWSP
VPKPTPSTSSSTSIESAILSATKSNRADTNVADYKTCKDNSIEDNMKYIDDEFSDSDIDIMTYKVTKDMIMEEGKRRC
SYNNDDNNINNTTAATTSSNNDNDNIISDSNNFTGYQKRIDDEASGISNINITTYTGSLLKNLLETDKRYTIVSSRQT
LRDTPNDPHAQQQTPTGDQEQSVEQQIIRDKNTTAQQQQQIHIIQDTDQKPIRIYHLQKNSSAIALPLKPVTESSISS
SFSSAAATPPPPPPLPPPVPIETPAPALPLLLPPFPPPPITKPASAPVSAPALAPTSAPASALASVLLQPPPPPPPTP
PPPPSQPMNDADSDVNVVDLTRSTSPLPPPPPPPLAHPNNSTQASSNTPQIPEFQPAAGSLFAVPYTFLIHSNYPIPL
YLPYGSCYYHNSQYAPYNYSFRYLDNFYLPPPPPPPPPLPLPITPPFYSYRNPYSYPYFHSNHCHYCHRRTPMNIDRQ
QYDLNYYRDHHHHHNYHHHHHYHYNQEGKQEVEEGNEEIARERTEIPEEMEREEVNEQDEQQFQQQQYYHHHHHHQHQ
QQQQQQQQYRHRQYNNESEPRHGEYQYYGYNHYQQQRRSALPQRQRQRYFRNRSIPPYVNSYNYNYNNNDDNVGGWYI
DSNTSTMREKRSGSVFEHDINNMREKRIKHQNNNSNSNSNNNTNNSSSNNNNKNNNSDDNNNTANKK"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61876",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61876.1 baculoviral IAP repeat-containing
 protein [Penaeus monodon majanivirus A]"
              }
            }
          },
          title "baculoviral IAP repeat-containing protein [Penaeus monodon
 majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 181,
          seq-data ncbieaa "MQNIFICARLRAQATLSASAIYIERSLAMPVSHFLKTECRFIDTQYNSHI
RPTHLLLTNVYFYITVLIILGSKHFPNPIHQEYVMHQRRLETFRWWPKETAQSPHKLAEAGFVSSGLGDTVYCFHCGV
YLEQWEKNDMPWKEHAKWSPRCRYVKMVKGDKFIIRWGPLDIDSVDIYDTIER"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61877",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61877.1 hypothetical protein [Penaeus monodon
 majanivirus A]"
              }
            }
          },
          title "hypothetical protein [Penaeus monodon majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 162,
          seq-data ncbieaa "MLSKCNRSRNRKNNKIKKKRNVTIERMDLDNELRRWLLDGWKTGGPLNRF
RRKNTGSTNNEAPTTSQGKSKQNISSVNNTTNNTTNNSNNKKVNDFSDNSNSCSNISNSSYSSSSSSSSSSSSSSGSN
SGNVSSSSSSRRSGNNRTIVKTNDEKKIGLNNYY"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61886",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61886.1 hypothetical protein [Penaeus monodon
 majanivirus A]"
              }
            }
          },
          title "hypothetical protein [Penaeus monodon majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 233,
          seq-data ncbieaa "MVAPHNSGIATRSFVRISGSGPCSAAGTVYMDFSRVCGAAAIVAARHRPP
QAVMQAQHRRLSNSMDAAFLRKMAAIEDSFAESARVEIGGSEARRARREVSNLRAGARAGGGASARPAAGEGRGVPAR
KRTSALPARDRSPDGGQRKRLRSHPSEENGKAIAPPQAKAPKTIYIVVGMKGDDARLRTVAAGQPFLLSAESLRRLMD
QTNRLLQENQALSQRVRAYQLNESVLR"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61889",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61889.1 hypothetical protein [Penaeus monodon
 majanivirus A]"
              }
            }
          },
          title "hypothetical protein [Penaeus monodon majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 504,
          seq-data ncbieaa "PLASVHIEHTLLGFLFKMLARNIVLLLIACQWLTSHPHRHRVFGYPIEDT
SSPFILKSPWTIRFPSAEKLDFICGICSSSSNSSSSTVEYGLTVSTNVYMYQYPIRLQRSQNCKALSITFNPIAETGL
INVRLTEQNKNKDSNTTTSSSINSSSSSSSSSSSSSSRSRSTRSNIINNNIEVVEGADIPYSISSDEKDEKHKDVYIY
ALNFQNVDSRIRLLGSAADGFATVAEDYWSYGISGITCIFRSFPLSTDDDKSLPDEQTEITTNDNTDSVNLLQPLGTP
SSTLPITALAALGTFSSTSPPAPTIPVSLTAIELVPSPPVEPSVLSSSTVAPPNNEREKENTTNDYFYEYYIYDDDDN
ETRGRLSTVSPRTTTVAAENAADDDYGNVEDDYNDGEEDPENKDRDRPSASSATTTTTNTTTTTNTTTTTTTSKRTTT
TVLEGIGDNESGNRVDPVITETEKEYIVGDIGDDRVMLIMLILIVVMLVLELCIIATGLCILGY"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61891",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61891.1 hypothetical protein [Penaeus monodon
 majanivirus A]"
              }
            }
          },
          title "hypothetical protein [Penaeus monodon majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 82,
          seq-data ncbieaa "MGNLRIPVIRSTDDYDRDEYVDDGNNIVDDNNNNDVNTNDDDDDDDNNND
DDDDDDDDNDNNNNNSDDDDDDDKLDTIMQPK"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61892",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61892.1 hypothetical protein [Penaeus monodon
 majanivirus A]"
              }
            }
          },
          title "hypothetical protein [Penaeus monodon majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 137,
          seq-data ncbieaa "MAASPIMFIRHVLLIWMFTNSIPFLRQADGHRHTFIKLAPTTIHRLKIKT
ILDDNKVIAIWGFCDRPPVCEVCHDPPLCKHRNNAWTGSYRLDGYTQVPCDPSTQDPCNSSTRLHFHEGETLCGRSLG
RKARGAAAR"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61894",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61894.1 dUTPase [Penaeus monodon majanivirus A]"
              }
            }
          },
          title "dUTPase [Penaeus monodon majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 146,
          seq-data ncbieaa "MSELNNVLKFKKLTEHAFTPSKGSKLAAGFDLYSAYDLTIPANGKSLVKT
DIQIELPEGCYGRIAPRSGLSWKHHLDVGAGVIDRDYRGDIGVILFNHAKSDYKVKKGDRIAQLICEKIIYPDIEQVE
ELTKTERGEGGYGSTGKQ"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61897",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61897.1 wsv427-like protein [Penaeus monodon
 majanivirus A]"
              }
            }
          },
          title "wsv427-like protein [Penaeus monodon majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 639,
          seq-data ncbieaa "MFNNRIDTYTNRARDLLAEINNLEYVLMNNPSLLPSKEKLHTPKDEESQL
DETHESYDESVASSSSSSSSSIFYKPECDLSLANVPLSEMVNHLTPLRTTIGARVREYVDCYGYLKKGIKINDALKNI
WRDIIDFSFTASDNIVKETVTNRLVNVWKDTLKGVQDHQDVFSNITEKRRSSVLILLRVIERIIDANMQYSQHFVKDA
LLCRLIRGATYSLNVKNTGVLPTVMRYIPLCTETSLQPPPPPPPLSSSSSSSSSSSSASSSSLSTSTIPETDTTLEYR
IKGVYGYGEFHPRCEDDHKIALRDICSDYFLKLLFEFLLSPHIHSTWRETYLIVFISSALKLIRNHEPFILALGKKKY
QQEQQEQEEEQEQQQQNTLDEQQTEYNCKKYFDDNIDFTMRNPITAAGCHPFVFLPTATDFWSFCKNGQPIYNVRNVL
VNARPSSFPPICGGPAGISQGKQLAWERYCQRQGVCDHGCECIACNRRNFLVDRVRQRIGQISGEARVTDVTKWHQGT
ILMITFYAFTGLLAHRNYARDMDDKADFFKNLITAIEKSHYMRSKDKLITSGGGESDTSNILHTITFGVEEEFRQKLK
SPQKKPIFITASIDDDILEGDNISRVQLTALRNRRVSEFNCWF"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61900",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61900.1 wsv303-like protein [Penaeus monodon
 majanivirus A]"
              }
            }
          },
          title "wsv303-like protein [Penaeus monodon majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 851,
          seq-data ncbieaa "MNKRLEKIAHNTYTRDNLCIGISLDFGKNKNNVSFLPKNKKPRLIDKTIM
VKSLNETSQKQFTALRLREEFPPCSILIFSMGSVIIVGLRKLGLTPLSVLQAVSAIADTQEHPITISNVEIVNVVSTF
NLFRLNMQSLKEFFKENKIAYTHVPDSFPGLFFKVLVPIRHLKRGETIGGYYTRAAREREAGDADAVKKNFRGKTVLV
FQVGKCTILGACGGDDIHVIFKMLFGFFWYFIDDTVHLSAEELECGNRTYHIPPLSWYLMIDFFLRTVPGYQLKPNVA
LADLAVEDKPSVVNRRGASDLLDAFDRLSITKDRMLRRTPRENHGAWFINCLEVQPSTEGRRGRGGGGVSASAQKRQN
LFEKELGQLKTVQEKSTANQWGVTTQQCYRRLLRLHNLTREIGEETPTPIHAPVSHTARCIFEGFEGYLPYTEEELSH
PERIRIEMPPKDARVTKFHQRRPLFSKESKAREVRSNHRMKIRQLDQLEEIYDLVGKAANEGRCDEIAAFLGVQPYII
SGLLNRDISVPLQESSVFGRDYTHMTPAQIFFNKNLYKEYSNAEVANENRRGVLDFMINKDANTESIKQDHNNFNNIK
EEAEQEQENEEEMEVEIDIDDDGNDFKDDDIYDDTSNDDDDDDDDDDDDDINNHDGDDDDDDDDNRDNDTEIYLDIHG
MEKFLEEQVQHGMIFAPLSNYRTQHHQSLQTRLILEDLKRKTNMNMLQQQQFKQQQQQQQQQQKLNSTMLNVKSLYKK
KESFLHKECCECKYKNKYESNFNFNTVTIPYTNVNTTTSFTTTTTSTSTTSSTTTTTTTEKLKLKRSKCITNSDITEN
KRQKYTLPKWKRQEHTLPQKE"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61901",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61901.1 hypothetical protein [Penaeus monodon
 majanivirus A]"
              }
            }
          },
          title "hypothetical protein [Penaeus monodon majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 325,
          seq-data ncbieaa "MPDRYIIPQRRKRRGGGGRGERRRSTSTYKDEQIKNDYHISQGIIHHTYR
NDEINNIRNRDNWRLLKSQYNDNNNNNNHRHDHHHHQHPHHRHGWQQLIDDHRVYRNFRKCNKAKSWSNASDEHEGVV
EMISYLNSNNEHIQKQQQQRISLYDTVTVSTETTEKEICNSVATVLTLPMSSLSQPPLLPSLLPPPLLSSSFSSSSSF
SSSSPSYFSNIPEEGIVKVFSLSSRDNNLFATVARNVKTFSYHRDHILYLAVVAFIFLRESSLVLTPTDMKSMFNSPV
VESLLKKLIIILQYENKTVETCMIIYNVERVAYDVICDYYN"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61913",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61913.1 wsv136-like protein [Penaeus monodon
 majanivirus A]"
              }
            }
          },
          title "wsv136-like protein [Penaeus monodon majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 139,
          seq-data ncbieaa "MCRVKNEFSKRLENFLKTQETYKFFLLLFLITAVCLCASHFIYYEIVYNS
AHKPLEYLHQLSTSLPRVNPRKAAILIYGLQDDDHSSSSSLLLSSVDLLPPKMMDNLRKTLIDPRSFLMTSYGTHLGR
SDYVDAPVYMS"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61930",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61930.1 semaphorin-1A-like protein [Penaeus
 monodon majanivirus A]"
              }
            }
          },
          title "semaphorin-1A-like protein [Penaeus monodon majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 575,
          seq-data ncbieaa "MLPYGKEDVRLSCFKTCFLEYLITTLSIFKEAICLNNIKTVDSKSNSKSI
TDSDSDNQHFLKFKKINFSHSTSERRICLNMYIERDDQFCMEYVVYVLLAIDMFISMLDKTDDKHNNNNDNDDNDANR
DNYDGKKTRRSLIDISYINRDVLLAEDCKSCENSLYGFAAALLENRVKPFSFLFLSRNKKNFNDYYAGDPTSVFFTAD
EIKSISRRIKTSRKNERIKIILSKYRKKIKDDSTKTENGLVLRARATNILIRCDDDNDIFMDVCVHIVDKESYVWIKT
IEKQTKQLNKSLIDYLDALVKIDDSVCPDRLRVVNDVIAFKERMTDNTHKYFPTFKFSGSRKYAAISRAIETLANECE
PHASTPQYHHTTIKEHCSYVKDNYGSFEPFMNTILPLLINETITNIWFNRLFLLNTQREVATNNDIRSTGFLSFLKGV
EREADKTKVMDVLRKFSNVCVTLRRNRLPFPLREYDRPTKRQKNTAVADNRYNESQPLRDIQYIQTLSKGTVNVKSLS
SHRTPKISGDIVVGKYGADDIVNETDNNSNTSRRLANIDAGQFLQYFTPIAKGLFAQ"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61932",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61932.1 hypothetical protein [Penaeus monodon
 majanivirus A]"
              }
            }
          },
          title "hypothetical protein [Penaeus monodon majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 118,
          seq-data ncbieaa "MSSKACRSFFFFDVKISTIIPAQHSRSVQFLSALAFCGSTSFQQQNTISR
ANETFTRHLEHDGDFLDRLFFDIKISTIRPAQHSRSLSGSVQILLSHYHFFFCMNKTMFYFVSVYSIL"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61938",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61938.1 wsv308-like protein [Penaeus monodon
 majanivirus A]"
              }
            }
          },
          title "wsv308-like protein [Penaeus monodon majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 592,
          seq-data ncbieaa "MFPMKRKRIQGDYTRGQSRHTSTTVIPARDDYLKDSLLFETRLESFIKRM
IERQAALSIQNAMDISDIISPNVGVDKDVLKMQAVYSTCDPDTYIENLLDFLKGEFAGNLIENNIILRNMKERRLKSL
ETIRRCLLEDAMCGDFTYDIFRSCVTLSDKREGSINADNNIVSNNPYRYRYQHQHCTPVIDAHEAFIVTMDEQNKSNR
SIDPSHSGNDEYNKESVPLDSDDDDDDVVDDDDDDDNDDSVDYKNDDTDKLERESKGDTETDIINNVVRMLNAILYDE
ESTRGEIGKRNDPCYIYNHVKNIFKKSQDILQTGECLSCNHQPIEKFLFTLDHNSSTRSMEDLLAQLEGVINFTFPKG
ESQLSSIEEITVLRIFFKNILLLKTATLVNRYRLFLSFVVTYSISHNLNHNIRVYFQKNRAYNTHCMMDRTVNPNVDQ
RVVDKFLLSKLSSEESVFILSLLAHFLGPCLETFLCPANLYDAAINYIDSANSEWGKRYVDVSLSSKIKRILDQIATR
ETFLKERRKRKSINTVIAYMTENRDSDTLLNPLHCIISDIGNVVEIITDGFRDLFIRTTSALNIISTLISPRNK"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61939",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61939.1 wsv306-like protein [Penaeus monodon
 majanivirus A]"
              }
            }
          },
          title "wsv306-like protein [Penaeus monodon majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 426,
          seq-data ncbieaa "MHTAHIKPNCSAVDIFCTDDRDCHVQCGSSALSVSSHTIIPEKDDYSVVT
IYLCNKSTKRCEPNTVLKGFYDMAMHDLKGLLKRKSTLSRGLHSYLQMYDGKYFILARARFYESLSLHLHDNENLEMA
IRNIRSEHQRKRRTYLYEKIKYLEHEGAVTDKDEVKKMMDWLSTPALTSQNPIAENTTITMDVYKQQLFGHRMPHPPR
DTTVLTRAEEYDDRRGEYVQNVTAVLPRQHDRHRRGEHLLKPTTTTTTLERHEELSSNLQTKSPSRELQFHCNEEMGA
GIMMAAFRDGDGIILPVCICTHPLYLTGPTCKHRTYDSVIDYELWEKSGIPEFLTDPFNDYTKANDICKRTTLNTTAV
FDEREGVFRCQPLAAHIGQSLQFRGPYEPALILDRDYIDTFNANHQQGFKINVDYIDLLGKFWQ"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61945",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61945.1 hypothetical protein [Penaeus monodon
 majanivirus A]"
              }
            }
          },
          title "hypothetical protein [Penaeus monodon majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 251,
          seq-data ncbieaa "MSKSYSVIGIMNKLAKKENRDISSSFIYKNTRDEDVYIGKEKQDEDEKDD
DDTSSNVVVMQELKRVQTQRDKNDDMPSSDKVHEDDDDDDDHTPSRDNEHKDTRVVKVHDDERRECNEDDGRNENVNV
FDLQKATNTTNTNTITVTNIKIDDNDTTNLHDQDAFKIERRLYFINSLKRVLARRKNKPANIFSIHRKKEKDDVTYIS
AEDFAINNTKEKERTDCHKCLIHVFNRLVDCFELLIKLKTRSSVL"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61950",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61950.1 hypothetical protein [Penaeus monodon
 majanivirus A]"
              }
            }
          },
          title "hypothetical protein [Penaeus monodon majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 224,
          seq-data ncbieaa "MCNFINSIASTSIIENGKDFLRTDESNTVIIILSEDLIRNNSSASITEQQ
RLIEKTCVSTLLSTIVQKKYINKALCNLSFYIQNKMRTKREYCIVIKNPRHIEIYSYKSVMCSCDRDDEYGTDTEEEY
ICEWTTKRCINVPYIPPNYKKKIKPPYVFIFTGRHIPSILRNFSDDICAPLFVKIIYITRDDEFANVGLDDLLEHLPY
KNCTMKVLIQRLKSNNRS"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61954",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61954.1 hypothetical protein [Penaeus monodon
 majanivirus A]"
              }
            }
          },
          title "hypothetical protein [Penaeus monodon majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 284,
          seq-data ncbieaa "MNKQGEGISPSLIFKSNSGDIKRPIDLQYWLGIINARPSKIAKQSFPGDD
KYWQNEENLKNFLVNEKTVDRVYHFVRSEFDLNEFRLLLRTYYDRRPLYVEAVFFIHLMDGRLIHTLSNLMTLTFNVQ
HFICSVVGIYKDRCAILKFIKDDELDARATWNTALPLQHLCYRAAAKLWENGRTTLCKIPPPRIWETIKATVAIRKIE
EGYERRICMAEHLVYPVDRAAVRERGEHSTIAYSPREVYLNINYFREDDYDSEEENSDEYYNSLVKRRPPTNNTYRVL"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61955",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61955.1 hypothetical protein [Penaeus monodon
 majanivirus A]"
              }
            }
          },
          title "hypothetical protein [Penaeus monodon majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 266,
          seq-data ncbieaa "MNRAIDFKCYTNLLNGVDSSVAENVEKAIEWINTRLYRNGKVYIFPDEAK
DWKDEEVTKIAFDKDRVDRVYYFMGFMFDPTTFHIAFRTHYKSRLLYVEFVTLNLDCHDETMLILTYNAQRFVDTAVR
SLHKDRDTVLQHMLDDRVNVCRSWRNISPLKQLCSRTITDQKIDFLSENIPKCIKDTLIKDEEISKIESKRARLIEEV
YGNIRHCWLPGNYRVTHFAENHMWDDDNSDDEPGEYRYHEHDRRTIDTFKREILEILDII"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61956",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61956.1 hypothetical protein [Penaeus monodon
 majanivirus A]"
              }
            }
          },
          title "hypothetical protein [Penaeus monodon majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 353,
          seq-data ncbieaa "MDIDGGKGGGLFVWKPETRILDSDFELLTSGSPPKKICASRNDYSGIYNA
INMNPNKGENNRGLLSSVTPVASALYRRRTSSKAGSPHKLLKTPPKPLTKFILDGTTISDGDMSTEQNQNSPVGEGRP
GGTKRTSIRPDVTIRSILSAEPDRPIGSHRTSSQSECETIHPLGKINSPTDADDSERESKSNNSPTDSLCLSKEIEET
LTGIENRVEMDRGGRAEFEADRHVRSRRTFVPGPDQLNGLAELTFDDIYISEEDDDDDDDHETNHPLEMANSITRVET
DRGGRAEFEADRHVRSRRTFVPGPDQLNGLAELTFDDIYISEEDDDDDDHETNHPLEMANSITHGIVDG"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61959",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61959.1 baculoviral IAP repeat-containing
 protein [Penaeus monodon majanivirus A]"
              }
            }
          },
          title "baculoviral IAP repeat-containing protein [Penaeus monodon
 majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 261,
          seq-data ncbieaa "MNSYSVKCKSLRYESERLLTFINWPFGSLAKSLAWNGFYRKCPDDLVCYF
CGICINVEDYSSLIEKGSGYINMKNIHYRSCPFIYGNSSMNVTREQEILISEMMRRYCKRYAKRTRTLYKTPKKRFLS
FTNFPNRYCNVNRMAESGFYFTGELDAVRCFSCQIDIYGWQHNDDPISVHAQWSPDCEYIQAIVDDNFIDSVRDKRYD
TLIDPNVYQTLKDCGIPDYIIRHAHRKHIIQSGELFTDAVQCLSVIEKLYNFKDK"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61960",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61960.1 hypothetical protein [Penaeus monodon
 majanivirus A]"
              }
            }
          },
          title "hypothetical protein [Penaeus monodon majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 492,
          seq-data ncbieaa "MYDYIDYKLLFFGEDVIHRILQIRNQALFVINGIKKVQEHDNLLFSKEKI
YTVNNVMELASCTASDDDADADTNNPILMDIIGVLIYGILNEGWSESTMRCIENELVDKCLVHKYITEEEYRYCRKKY
NNSDNDKSKNDNIPYIIDFLWNMCYFNKLYRDRYLNKNMTLRLIEKSTSSMSRNALIKHDLSFFHLPTAAGIAIYEYL
THRDTNLRLVKEAAKTLVAEEKYLRYSRNIRNTIEDDTTNKKKNIRLAYVNKRSKTSCALLLCKSASKINYSKKDTKS
KLFNLQWLFNVTDIQNQKFKKCLGVIYDSIDVGDMLPNTSLSQKSKPLDWYYRFMNTGIVPFDQHNGCSEFRQYHREE
GNEEKEHIYHNIDEENTYVIMNINTNTIDRVTPTTLSNTDVNTVQIHNKSMKHYNTLPVDKKKNRTHRTETFTPTVPS
SSMMMSSSVSVYDIPKKATTSIYNCKSLMIQCHKCQRYNISPYVVKCNRCRC"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61961",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61961.1 hypothetical protein [Penaeus monodon
 majanivirus A]"
              }
            }
          },
          title "hypothetical protein [Penaeus monodon majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 207,
          seq-data ncbieaa "MREYGKIYYSSMVGLLNALTEAGYPISIHDLSDRDGIKSSYRDDGVIITH
ASEEALCLTLYCALGLLMEKGITVNLKRSNAAEYLLHFPTVKVEVRPNCETNTIVRGRIQRNDDRVAIGEEPVNTSDI
YNGGKAYVYQNKALYSHSSGHLKTVIIESKGRGHTANDTLSTVLYIDSYDRSRLPDIIEDIAKYPLKFFRAIILRGKLL"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61962",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61962.1 baculoviral IAP repeat-containing
 protein [Penaeus monodon majanivirus A]"
              }
            }
          },
          title "baculoviral IAP repeat-containing protein [Penaeus monodon
 majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 577,
          seq-data ncbieaa "MDISHPNSFLQQEEYRLTTFSEHETSPIDTTRVYRNMPFVSYDCLRYEEK
RLGTFIDWPHEWLKPSELAADGFYYLRKDDHCACIFCRGIVGAWEKTDTPRGEHERHFPHCPFIRGQPTGNVPIVQGN
ILARLPVLHASLDDCGSRRVEADACRFSSRHMAGSYPECRGSPDETGLYQHSGPKRSDFTTLDSRLKSYVEGWPILGL
HARDLAEAGFFYCGLSDHVRCFHCGKGLRNWISSDIPWKEHARWYPKCRFVLLTKGQDYINEVQQEHPPYTRTASTTI
YSGIDKAASGSSNSPTEVQSVSEQELDLLMSLDVVKCVLEMGYPSCIVRACLKDRVEQTGEPYFVIEPCLEDVIHRMS
VTENRVADNCSESTEREPSYGGTDISSSIFQQQAPIFPSTSETMGNRVASIQSLPSESTEVELLSQQQQVSILPSTPE
IIRQQASIFPTAPLNTTLSSSSSSSERENTETEANHILQIVEEIISQVSEEEKEGYEEKEGYEEKEEFKSRQSYEQLK
EELERMREDRTCKVCMDSEIAIVFLPCSHMVTCVSCAVSLRQCPICRRDIKHVCRAIVS"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61963",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61963.1 baculoviral IAP repeat-containing
 protein [Penaeus monodon majanivirus A]"
              }
            }
          },
          title "baculoviral IAP repeat-containing protein [Penaeus monodon
 majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 728,
          seq-data ncbieaa "MSDSGNMDNATDIDRDLTRMCEETHRANTYGDYDQYVHHPKIMAKAGFYE
LDKRLICFRCGLNIDIDDIGKTDDIVEVHQTKRPSCVFARNLPSLPRSRLVVASDLRYETKRLETFIDWPLGEVIEAR
DLAAAGFFYSRRNVECICFCCKVRISFKPCTRHGCCARHASTSQYITEKHLNSAYCDFAMEKVTGNVPIRINKIISKI
IGHQETDNEKNSCGNNEGRSRGASGEEDGGDSRESENDSSSRISAGNILSEMGIESYSGPEHENFVTPISRKYTYDDW
PLGDIMLPKKLIAAGFFYSGYGDHVNCFHCGVGIRNWEKGDDPWTQHAMWNPECHYVRTLMGQTFIDTAKTVKNVLFR
KPTRERTYSVTDEDWNAVMGLDITKDLLISLFPVAAVRDCVREQILVHKIPFITKHEGLYYVREKIAAYGNMDITLAK
KISRMSTEYIQNYHAVYSVMYFLMDLKRIFIMNQDVDRNVSIDLFKPDQSLFLFRQPSYNPAATPMYQYRLSFTIPHP
SHNSDFKFQQEYRIALEQLNRGLPVHVTNLYFRERDSISLFDIWNHLIRLSIIINDNNRSSGSSINSSNNNNNNNNNN
NNNNNSSSSSSSSSSSSNNNNNNNNNNNNNNNNNNDDDDDDDDDDDKNSYDVFSVIREHHQRQEDEERSSEIKEAPDS
VIGDRYLCKICLNSCVKTVTLPCRHMSCCGNCIVTQTKCPICRKKIKYTLEPYL"
        }
      },
      {
        id {
          ddbj {
            accession "BDT61964",
            version 1
          }
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">BDT61964.1 baculoviral IAP repeat-containing
 protein [Penaeus monodon majanivirus A]"
              }
            }
          },
          title "baculoviral IAP repeat-containing protein [Penaeus monodon
 majanivirus A]"
        },
        inst {
          repr raw,
          mol aa,
          length 719,
          seq-data ncbieaa "MSQNIPLHIDDQQYHVDFSGMTSEYNRRMTFGHESSSYVINPKILAKAGF
MRCGLSSKQYIVCVDCHLYIDIDTIMPWDDIVKIHKDKNPLCLLARNLPPLPRSKIFENRFDPATGLPADLRYEIKRL
ETFIDWPNKVVQAKDLAAAGFFYTRDTTCCICFSCHKYMNFSYDKPIVIRDMHLKKSPNCKFAQGKSVGNIPIRMSNI
LNSIVNQKRVAVTSATATATATATVDQEDLIKTRYEYAGPLHDYYKTTQSRLESFLNWPHRHRQKIEALVEAGFYYLG
ISDHVQCFHCGACLRNWEEDDDPWRVHAEWYPECYYVQLKKGKEFIDKVRLEVSRIPMWRHSSSTSQDTSEHDWNILL
GLDYTKRLILVQGFPVVAVRDALCEQILQTGIPFSREVDCTSAVRKKIESYEEKNNNKLLTAERIANYRTVYAVMRFI
QMVMGKSYNRSYPTGCIENSVYLYKQSAFNDSGDPNIPMYRHEISFNLPYPDHTYTGYHIESEYRMTLDRLNSYEITS
RISSILDQEQSASLTAIWGYLVRQTTTTKKTNTSASNVDRGDITNALDSVGEAAAAAAAAVDEEEELEVEEEEEELEV
EEEEEDIVVSELENQSAIGITIAPVPTMTPIPTVHTNGRDLIFHTRLYIISNANSRRVVPYQNFNSIIRGEPDIQHRQ
FNNRDEEAESKLSRSTGIEEGSMVNQTRTSIKTRGISVRSAWMPK"
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
    }
  }
}
