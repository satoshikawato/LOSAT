//! NCBI composition-adjustment helpers used by `blastp` redo post-processing.
//!
//! Reference: ncbi-blast/c++/include/algo/blast/composition_adjustment/composition_adjustment.h
//! Reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c
//! Reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/compo_mode_condition.c
//! Reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/unified_pvalues.c
//! Reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/matrix_frequency_data.c

use anyhow::{bail, Result};
use std::sync::OnceLock;

use crate::config::ScoringMatrix;
use crate::utils::matrix::{protein_score_direct, BLASTAA_SIZE};

use super::redo_alignment::{BlastCompoAdjustMode, EMatrixAdjustRule};

// NCBI reference: ncbi-blast/c++/include/algo/blast/composition_adjustment/composition_constants.h:39-44
// ```c
// #define COMPO_SCORE_MIN INT2_MIN
// #define COMPO_NUM_TRUE_AA 20
// #define COMPO_LARGEST_ALPHABET 28
// ```
const COMPO_SCORE_MIN: i32 = i16::MIN as i32;
pub(crate) const COMPO_NUM_TRUE_AA: usize = 20;
const COMPO_LARGEST_ALPHABET: usize = BLASTAA_SIZE;

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:77-80
// ```c
// static const double kFixedReBlosum62 = 0.44;
// static const double kMaximumXscore = -1.0;
// ```
#[allow(dead_code)]
const K_FIXED_RE_BLOSUM62: f64 = 0.44;
const K_MAXIMUM_X_SCORE: f64 = -1.0;

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:64-75
// ```c
// static const int kLambdaIterationLimit = 100;
// static const double kLambdaErrorTolerance = 0.0000001;
// static const double kLambdaFunctionTolerance = 0.00001;
// static const double kCompoAdjustErrTolerance = 0.00000001;
// static const int kCompoAdjustIterationLimit = 2000;
// ```
const K_LAMBDA_ITERATION_LIMIT: usize = 100;
const K_LAMBDA_ERROR_TOLERANCE: f64 = 1.0e-7;
const K_LAMBDA_FUNCTION_TOLERANCE: f64 = 1.0e-5;
pub(crate) const K_COMPO_ADJUST_ERR_TOLERANCE: f64 = 1.0e-8;
pub(crate) const K_COMPO_ADJUST_ITERATION_LIMIT: usize = 2000;

// NCBI reference: ncbi-blast/c++/include/algo/blast/composition_adjustment/composition_adjustment.h:42-45
// ```c
// enum { eGapChar = 0, eBchar = 2, eCchar = 3,  eDchar = 4,  eEchar = 5,
//        eIchar = 9, eLchar = 11,  eNchar = 13, eQchar = 15, eXchar = 21,
//        eZchar = 23,  eSelenocysteine = 24, eStopChar = 25,
//        eOchar = 26,  eJchar = 27};
// ```
const E_BCHAR: usize = 2;
const E_CCHAR: usize = 3;
const E_DCHAR: usize = 4;
const E_ECHAR: usize = 5;
const E_ICHAR: usize = 9;
const E_LCHAR: usize = 11;
const E_NCHAR: usize = 13;
const E_QCHAR: usize = 15;
const E_XCHAR: usize = 21;
const E_ZCHAR: usize = 23;
const E_SELENOCYSTEINE: usize = 24;
const E_STOP_CHAR: usize = 25;
const E_OCHAR: usize = 26;
const E_JCHAR: usize = 27;

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1097-1100
// ```c
// #define LambdaRatioLowerBound 0.5
// ```
const LAMBDA_RATIO_LOWER_BOUND: f64 = 0.5;

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/compo_mode_condition.c:36-42
// ```c
// #define QUERY_MATCH_DISTANCE_THRESHOLD 0.16
// #define LENGTH_RATIO_THRESHOLD 3.0
// #define ANGLE_DEGREE_THRESHOLD 70.0
// #define HIGH_PAIR_THRESHOLD 0.4
// #define LENGTH_LOWER_THRESHOLD 50
// ```
const QUERY_MATCH_DISTANCE_THRESHOLD: f64 = 0.16;
const LENGTH_RATIO_THRESHOLD: f64 = 3.0;
const ANGLE_DEGREE_THRESHOLD: f64 = 70.0;
const HIGH_PAIR_THRESHOLD: f64 = 0.4;
const LENGTH_LOWER_THRESHOLD: i32 = 50;

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/unified_pvalues.c:37-48
// ```c
// #define LENGTH_MIN             40
// #define LENGTH_MAX             4000
// #define LAMBDA_BIN_SIZE      0.001
// #define LAMBDA_BIN_TOTAL       565
// #define LOW_LAMBDA_BIN_CUT     35
// ```
const COMPO_MIN_LAMBDA: f64 = 0.034;
const LAMBDA_BIN_SIZE: f64 = 0.001;
const LAMBDA_BIN_TOTAL: usize = 565;
const LOW_LAMBDA_BIN_CUT: usize = 35;

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:44-54
// ```c
// static int trueCharPositions[COMPO_NUM_TRUE_AA] =
// {1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22};
// static int alphaConvert[COMPO_LARGEST_ALPHABET] =
//   {(-1), 0, (-1),  4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5, 1, 15,
//    16, 19, 17, (-1), 18, (-1), (-1), (-1), (-1), (-1)};
// ```
const TRUE_CHAR_POSITIONS: [usize; COMPO_NUM_TRUE_AA] = [
    1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22,
];
const ALPHA_CONVERT: [i32; COMPO_LARGEST_ALPHABET] = [
    -1, 0, -1, 4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5, 1, 15, 16, 19, 17, -1, 18, -1, -1, -1,
    -1, -1,
];

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/matrix_frequency_data.c:383-540
// ```c
// static const double BLOSUM62_JOINT_PROBS[COMPO_NUM_TRUE_AA][COMPO_NUM_TRUE_AA] = { ... };
// static const double BLOSUM62_bg[COMPO_NUM_TRUE_AA] = { ... };
// ```
const BLOSUM62_JOINT_PROBS: [[f64; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA] = [
    [
        2.1497573378347484e-02,
        2.3470224274721213e-03,
        1.9493235258876179e-03,
        2.1674844853066858e-03,
        1.5903351423026848e-03,
        1.9242657898716525e-03,
        2.9879059292799641e-03,
        5.8158526388051033e-03,
        1.1076584657559144e-03,
        3.1880644746334580e-03,
        4.4186245468471547e-03,
        3.3466571942021082e-03,
        1.3412107617355408e-03,
        1.6360627863999076e-03,
        2.1568959784943114e-03,
        6.2524987419815400e-03,
        3.7180506975672363e-03,
        4.0281679108936688e-04,
        1.2999956675626666e-03,
        5.0679056444508912e-03,
    ],
    [
        2.3470224274721213e-03,
        1.7757465118386322e-02,
        1.9786027128591904e-03,
        1.5865480081162602e-03,
        3.9365984789376245e-04,
        2.4858611089731411e-03,
        2.6933867548771758e-03,
        1.7221140903704937e-03,
        1.2407382229440791e-03,
        1.2435878276496955e-03,
        2.4193952633248727e-03,
        6.2339060289407083e-03,
        8.0309461712520876e-04,
        9.3181986323789834e-04,
        9.5783034332718700e-04,
        2.2660898636037261e-03,
        1.7802796534180537e-03,
        2.6571979312581875e-04,
        9.2634607111251918e-04,
        1.5810185245264004e-03,
    ],
    [
        1.9493235258876179e-03,
        1.9786027128591904e-03,
        1.4140291972553610e-02,
        3.7201973506001745e-03,
        4.3845466068066216e-04,
        1.5304436972610567e-03,
        2.2097156829738759e-03,
        2.8591871815612977e-03,
        1.4301072616183181e-03,
        9.9437221166923172e-04,
        1.3690958423974782e-03,
        2.4402105140841090e-03,
        5.2943633069226512e-04,
        7.5004227978192801e-04,
        8.6016459857770028e-04,
        3.1466019144814608e-03,
        2.2360795375444384e-03,
        1.6159545671597605e-04,
        7.0048422794024819e-04,
        1.2014015528772706e-03,
    ],
    [
        2.1674844853066858e-03,
        1.5865480081162602e-03,
        3.7201973506001745e-03,
        2.1274574617480089e-02,
        3.9909227141697264e-04,
        1.6481246723433428e-03,
        4.9158017471929655e-03,
        2.5221102126636373e-03,
        9.5384849402143984e-04,
        1.2347404942429857e-03,
        1.5202051791453383e-03,
        2.4453087721980561e-03,
        4.6429229320514104e-04,
        7.6023722413111566e-04,
        1.2373315413524663e-03,
        2.8035127901697272e-03,
        1.8961512776990257e-03,
        1.6218020183662784e-04,
        5.9842263937853702e-04,
        1.3158365660538270e-03,
    ],
    [
        1.5903351423026848e-03,
        3.9365984789376245e-04,
        4.3845466068066216e-04,
        3.9909227141697264e-04,
        1.1931352277704348e-02,
        3.0937204045913537e-04,
        3.8338775043186374e-04,
        7.6951976030099293e-04,
        2.2976387481074697e-04,
        1.0956590131781735e-03,
        1.5682982157153873e-03,
        5.0124929379033781e-04,
        3.7717165634097634e-04,
        5.1389991547056834e-04,
        3.6111795849154795e-04,
        1.0432626586831986e-03,
        9.3041313726939057e-04,
        1.4474923964368156e-04,
        3.4603772624580643e-04,
        1.3606607271146112e-03,
    ],
    [
        1.9242657898716525e-03,
        2.4858611089731411e-03,
        1.5304436972610567e-03,
        1.6481246723433428e-03,
        3.0937204045913537e-04,
        7.3292255467189687e-03,
        3.5385780499965817e-03,
        1.3683038039160171e-03,
        1.0489026828741754e-03,
        8.9102936026571569e-04,
        1.6174411456311808e-03,
        3.0968229715707327e-03,
        7.3993258722701268e-04,
        5.4255147972143906e-04,
        8.4668181752066874e-04,
        1.8931125300036275e-03,
        1.3796838284921874e-03,
        2.2737931366728891e-04,
        6.7584155312457842e-04,
        1.1660966117775285e-03,
    ],
    [
        2.9879059292799641e-03,
        2.6933867548771758e-03,
        2.2097156829738759e-03,
        4.9158017471929655e-03,
        3.8338775043186374e-04,
        3.5385780499965817e-03,
        1.6133927472163669e-02,
        1.9380952488713059e-03,
        1.3667885452189439e-03,
        1.2192061706431622e-03,
        2.0030316026648431e-03,
        4.1322603720305197e-03,
        6.7909745467514783e-04,
        8.5179405867513139e-04,
        1.4216207127018586e-03,
        2.9539180653600089e-03,
        2.0493063257644955e-03,
        2.6488552587183780e-04,
        8.7044186256788659e-04,
        1.6987763526262680e-03,
    ],
    [
        5.8158526388051033e-03,
        1.7221140903704937e-03,
        2.8591871815612977e-03,
        2.5221102126636373e-03,
        7.6951976030099293e-04,
        1.3683038039160171e-03,
        1.9380952488713059e-03,
        3.7804346453413303e-02,
        9.5813607255887238e-04,
        1.3849118546156933e-03,
        2.0864716056392773e-03,
        2.5392537741810947e-03,
        7.3281559749652399e-04,
        1.1976708695723554e-03,
        1.3641171883713547e-03,
        3.8342830901664762e-03,
        2.1858459940987062e-03,
        4.0740829083805248e-04,
        8.3467413018106177e-04,
        1.8218235950233687e-03,
    ],
    [
        1.1076584657559144e-03,
        1.2407382229440791e-03,
        1.4301072616183181e-03,
        9.5384849402143984e-04,
        2.2976387481074697e-04,
        1.0489026828741754e-03,
        1.3667885452189439e-03,
        9.5813607255887238e-04,
        9.2802502369336622e-03,
        5.8089627083019206e-04,
        9.8696608463236094e-04,
        1.1873625842258938e-03,
        3.8264639620910225e-04,
        8.1041076335565583e-04,
        4.7770135861914477e-04,
        1.1052034635193162e-03,
        7.4371746073077327e-04,
        1.5168037757411286e-04,
        1.5213771111755425e-03,
        6.4882907765797669e-04,
    ],
    [
        3.1880644746334580e-03,
        1.2435878276496955e-03,
        9.9437221166923172e-04,
        1.2347404942429857e-03,
        1.0956590131781735e-03,
        8.9102936026571569e-04,
        1.2192061706431622e-03,
        1.3849118546156933e-03,
        5.8089627083019206e-04,
        1.8441526588740136e-02,
        1.1382470627796603e-02,
        1.5655862274689192e-03,
        2.5081290988482057e-03,
        3.0458868657559346e-03,
        1.0068164685944146e-03,
        1.7225081689171561e-03,
        2.6953622613315018e-03,
        3.6183761166072852e-04,
        1.3821121844492116e-03,
        1.1972663837662637e-02,
    ],
    [
        4.4186245468471547e-03,
        2.4193952633248727e-03,
        1.3690958423974782e-03,
        1.5202051791453383e-03,
        1.5682982157153873e-03,
        1.6174411456311808e-03,
        2.0030316026648431e-03,
        2.0864716056392773e-03,
        9.8696608463236094e-04,
        1.1382470627796603e-02,
        3.7141460156350926e-02,
        2.4634345023228079e-03,
        4.9293545515183088e-03,
        5.4151301166464015e-03,
        1.4146090399381900e-03,
        2.4277107072013821e-03,
        3.3238031308707055e-03,
        7.3206640617832933e-04,
        2.2096734692836624e-03,
        9.4786263030457313e-03,
    ],
    [
        3.3466571942021082e-03,
        6.2339060289407083e-03,
        2.4402105140841090e-03,
        2.4453087721980561e-03,
        5.0124929379033781e-04,
        3.0968229715707327e-03,
        4.1322603720305197e-03,
        2.5392537741810947e-03,
        1.1873625842258938e-03,
        1.5655862274689192e-03,
        2.4634345023228079e-03,
        1.6113385590544604e-02,
        9.0876633395557617e-04,
        9.4875149773685364e-04,
        1.5773020912564391e-03,
        3.1016069999481111e-03,
        2.3467014804084987e-03,
        2.7198500003555514e-04,
        9.9908866586876396e-04,
        1.9360424083099779e-03,
    ],
    [
        1.3412107617355408e-03,
        8.0309461712520876e-04,
        5.2943633069226512e-04,
        4.6429229320514104e-04,
        3.7717165634097634e-04,
        7.3993258722701268e-04,
        6.7909745467514783e-04,
        7.3281559749652399e-04,
        3.8264639620910225e-04,
        2.5081290988482057e-03,
        4.9293545515183088e-03,
        9.0876633395557617e-04,
        4.0477309321969848e-03,
        1.1901770463553603e-03,
        4.0824445213456919e-04,
        8.5603787638552766e-04,
        1.0095451907679563e-03,
        1.9872537223131380e-04,
        5.7145288352831449e-04,
        2.3123361470140736e-03,
    ],
    [
        1.6360627863999076e-03,
        9.3181986323789834e-04,
        7.5004227978192801e-04,
        7.6023722413111566e-04,
        5.1389991547056834e-04,
        5.4255147972143906e-04,
        8.5179405867513139e-04,
        1.1976708695723554e-03,
        8.1041076335565583e-04,
        3.0458868657559346e-03,
        5.4151301166464015e-03,
        9.4875149773685364e-04,
        1.1901770463553603e-03,
        1.8277684015431908e-02,
        5.2528021756783813e-04,
        1.1939618185901600e-03,
        1.1624184369750680e-03,
        8.4917468952377874e-04,
        4.2392005745634370e-03,
        2.5763052227920180e-03,
    ],
    [
        2.1568959784943114e-03,
        9.5783034332718700e-04,
        8.6016459857770028e-04,
        1.2373315413524663e-03,
        3.6111795849154795e-04,
        8.4668181752066874e-04,
        1.4216207127018586e-03,
        1.3641171883713547e-03,
        4.7770135861914477e-04,
        1.0068164685944146e-03,
        1.4146090399381900e-03,
        1.5773020912564391e-03,
        4.0824445213456919e-04,
        5.2528021756783813e-04,
        1.9066033679132538e-02,
        1.6662567934883051e-03,
        1.3511005665728870e-03,
        1.4152209821874487e-04,
        4.5224391125285910e-04,
        1.2451325046931832e-03,
    ],
    [
        6.2524987419815400e-03,
        2.2660898636037261e-03,
        3.1466019144814608e-03,
        2.8035127901697272e-03,
        1.0432626586831986e-03,
        1.8931125300036275e-03,
        2.9539180653600089e-03,
        3.8342830901664762e-03,
        1.1052034635193162e-03,
        1.7225081689171561e-03,
        2.4277107072013821e-03,
        3.1016069999481111e-03,
        8.5603787638552766e-04,
        1.1939618185901600e-03,
        1.6662567934883051e-03,
        1.2585947097159817e-02,
        4.7004857686835334e-03,
        2.8731729176487776e-04,
        1.0299846310599138e-03,
        2.3587292053265561e-03,
    ],
    [
        3.7180506975672363e-03,
        1.7802796534180537e-03,
        2.2360795375444384e-03,
        1.8961512776990257e-03,
        9.3041313726939057e-04,
        1.3796838284921874e-03,
        2.0493063257644955e-03,
        2.1858459940987062e-03,
        7.4371746073077327e-04,
        2.6953622613315018e-03,
        3.3238031308707055e-03,
        2.3467014804084987e-03,
        1.0095451907679563e-03,
        1.1624184369750680e-03,
        1.3511005665728870e-03,
        4.7004857686835334e-03,
        1.2514818886617953e-02,
        2.8575770858467209e-04,
        9.4161039895612720e-04,
        3.6402328079338207e-03,
    ],
    [
        4.0281679108936688e-04,
        2.6571979312581875e-04,
        1.6159545671597605e-04,
        1.6218020183662784e-04,
        1.4474923964368156e-04,
        2.2737931366728891e-04,
        2.6488552587183780e-04,
        4.0740829083805248e-04,
        1.5168037757411286e-04,
        3.6183761166072852e-04,
        7.3206640617832933e-04,
        2.7198500003555514e-04,
        1.9872537223131380e-04,
        8.4917468952377874e-04,
        1.4152209821874487e-04,
        2.8731729176487776e-04,
        2.8575770858467209e-04,
        6.4699301717154852e-03,
        8.8744160259272527e-04,
        3.5578318710317554e-04,
    ],
    [
        1.2999956675626666e-03,
        9.2634607111251918e-04,
        7.0048422794024819e-04,
        5.9842263937853702e-04,
        3.4603772624580643e-04,
        6.7584155312457842e-04,
        8.7044186256788659e-04,
        8.3467413018106177e-04,
        1.5213771111755425e-03,
        1.3821121844492116e-03,
        2.2096734692836624e-03,
        9.9908866586876396e-04,
        5.7145288352831449e-04,
        4.2392005745634370e-03,
        4.5224391125285910e-04,
        1.0299846310599138e-03,
        9.4161039895612720e-04,
        8.8744160259272527e-04,
        1.0246100213822419e-02,
        1.5489827890922993e-03,
    ],
    [
        5.0679056444508912e-03,
        1.5810185245264004e-03,
        1.2014015528772706e-03,
        1.3158365660538270e-03,
        1.3606607271146112e-03,
        1.1660966117775285e-03,
        1.6987763526262680e-03,
        1.8218235950233687e-03,
        6.4882907765797669e-04,
        1.1972663837662637e-02,
        9.4786263030457313e-03,
        1.9360424083099779e-03,
        2.3123361470140736e-03,
        2.5763052227920180e-03,
        1.2451325046931832e-03,
        2.3587292053265561e-03,
        3.6402328079338207e-03,
        3.5578318710317554e-04,
        1.5489827890922993e-03,
        1.9631915140537640e-02,
    ],
];

const BLOSUM62_BG: [f64; COMPO_NUM_TRUE_AA] = [
    7.4216205067993410e-02,
    5.1614486141284638e-02,
    4.4645808512757915e-02,
    5.3626000838554413e-02,
    2.4687457167944848e-02,
    3.4259650591416023e-02,
    5.4311925684587502e-02,
    7.4146941452644999e-02,
    2.6212984805266227e-02,
    6.7917367618953756e-02,
    9.8907868497150955e-02,
    5.8155682303079680e-02,
    2.4990197579643110e-02,
    4.7418459742284751e-02,
    3.8538003320306206e-02,
    5.7229029476494421e-02,
    5.0891364550287033e-02,
    1.3029956129972148e-02,
    3.2281512313758580e-02,
    7.2919098205619245e-02,
];

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/matrix_freq_ratios.c:438-466
// ```c
// static const double BLOSUM62_FREQRATIOS[BLASTAA_SIZE][BLASTAA_SIZE] =
// {{0.00000000e+00, 0.00000000e+00, ...},
//  {0.00000000e+00, 3.90294070e+00, ...},
//  ...};
// ```
const BLOSUM62_START_FREQ_RATIOS: [[f64; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET] =
    include!("blosum62_start_freq_ratios.inc.rs");

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/unified_pvalues.c:52-205
// ```c
// static double P_lambda_table[LAMBDA_BIN_TOTAL] = { ... };
// ```
const P_LAMBDA_TABLE: [f64; LAMBDA_BIN_TOTAL] = [
    1.247750e-07,
    1.247750e-07,
    1.247750e-07,
    1.247750e-07,
    1.247750e-07,
    1.247750e-07,
    1.247750e-07,
    1.247750e-07,
    1.247750e-07,
    1.247750e-07,
    1.247750e-07,
    1.247750e-07,
    2.495500e-07,
    2.495500e-07,
    2.495500e-07,
    2.495500e-07,
    2.495500e-07,
    2.495500e-07,
    2.495500e-07,
    2.495500e-07,
    2.495500e-07,
    2.495500e-07,
    3.743250e-07,
    4.991000e-07,
    4.991000e-07,
    4.991000e-07,
    4.991000e-07,
    7.486490e-07,
    7.486490e-07,
    7.486490e-07,
    7.486490e-07,
    8.734240e-07,
    9.981990e-07,
    9.981990e-07,
    1.122974e-06,
    1.122974e-06,
    1.372523e-06,
    1.497298e-06,
    1.622073e-06,
    1.871622e-06,
    1.871622e-06,
    1.996397e-06,
    2.370721e-06,
    2.620270e-06,
    2.620270e-06,
    2.620270e-06,
    2.745045e-06,
    2.745045e-06,
    2.745045e-06,
    2.745045e-06,
    2.869820e-06,
    3.493695e-06,
    3.493695e-06,
    3.493695e-06,
    3.868019e-06,
    4.117568e-06,
    4.367117e-06,
    4.866216e-06,
    4.866216e-06,
    5.240540e-06,
    5.240540e-06,
    5.864415e-06,
    6.363514e-06,
    6.737838e-06,
    7.236937e-06,
    7.860812e-06,
    8.110361e-06,
    8.484685e-06,
    9.233334e-06,
    9.857209e-06,
    1.035631e-05,
    1.060586e-05,
    1.085541e-05,
    1.210316e-05,
    1.272703e-05,
    1.322613e-05,
    1.409955e-05,
    1.509775e-05,
    1.647027e-05,
    1.784279e-05,
    1.859144e-05,
    1.921532e-05,
    2.008874e-05,
    2.133649e-05,
    2.245946e-05,
    2.333288e-05,
    2.445585e-05,
    2.557882e-05,
    2.732567e-05,
    2.832387e-05,
    3.019549e-05,
    3.231666e-05,
    3.468738e-05,
    3.568558e-05,
    3.718288e-05,
    3.942882e-05,
    4.292251e-05,
    4.479413e-05,
    4.666575e-05,
    4.903647e-05,
    5.140719e-05,
    5.352836e-05,
    5.539998e-05,
    5.714683e-05,
    5.964232e-05,
    6.263691e-05,
    6.550673e-05,
    6.825177e-05,
    7.124636e-05,
    7.411618e-05,
    7.611258e-05,
    7.910717e-05,
    8.235131e-05,
    8.534590e-05,
    8.971302e-05,
    9.457923e-05,
    9.907112e-05,
    1.028144e-04,
    1.068072e-04,
    1.124220e-04,
    1.150423e-04,
    1.202828e-04,
    1.241509e-04,
    1.297657e-04,
    1.352558e-04,
    1.408707e-04,
    1.474838e-04,
    1.535977e-04,
    1.605851e-04,
    1.679469e-04,
    1.766811e-04,
    1.849162e-04,
    1.943991e-04,
    2.023847e-04,
    2.111190e-04,
    2.208514e-04,
    2.298352e-04,
    2.390685e-04,
    2.500487e-04,
    2.599059e-04,
    2.705118e-04,
    2.826149e-04,
    2.953419e-04,
    3.058230e-04,
    3.209207e-04,
    3.336477e-04,
    3.469986e-04,
    3.607238e-04,
    3.758215e-04,
    3.915431e-04,
    4.100098e-04,
    4.283517e-04,
    4.454458e-04,
    4.661584e-04,
    4.845003e-04,
    5.032165e-04,
    5.194372e-04,
    5.399003e-04,
    5.642314e-04,
    5.905589e-04,
    6.173855e-04,
    6.407184e-04,
    6.646751e-04,
    6.896300e-04,
    7.199503e-04,
    7.487733e-04,
    7.780954e-04,
    8.119093e-04,
    8.442260e-04,
    8.782895e-04,
    9.142246e-04,
    9.507836e-04,
    9.890894e-04,
    1.025898e-03,
    1.060710e-03,
    1.102010e-03,
    1.141065e-03,
    1.179870e-03,
    1.225413e-03,
    1.268086e-03,
    1.313254e-03,
    1.362665e-03,
    1.413199e-03,
    1.465479e-03,
    1.514516e-03,
    1.568668e-03,
    1.620325e-03,
    1.678595e-03,
    1.738861e-03,
    1.800375e-03,
    1.863511e-03,
    1.933760e-03,
    1.999641e-03,
    2.066271e-03,
    2.136269e-03,
    2.208639e-03,
    2.282256e-03,
    2.363235e-03,
    2.451451e-03,
    2.535050e-03,
    2.626385e-03,
    2.715973e-03,
    2.811676e-03,
    2.906879e-03,
    3.008695e-03,
    3.108390e-03,
    3.220687e-03,
    3.337351e-03,
    3.458507e-03,
    3.587275e-03,
    3.717165e-03,
    3.853419e-03,
    3.998407e-03,
    4.150882e-03,
    4.307349e-03,
    4.470305e-03,
    4.650105e-03,
    4.835520e-03,
    5.028047e-03,
    5.239291e-03,
    5.459394e-03,
    5.699835e-03,
    5.944019e-03,
    6.201679e-03,
    6.472939e-03,
    6.767532e-03,
    7.067490e-03,
    7.402760e-03,
    7.779704e-03,
    8.163261e-03,
    8.581007e-03,
    9.043672e-03,
    9.538778e-03,
    1.006046e-02,
    1.061571e-02,
    1.121351e-02,
    1.187918e-02,
    1.259427e-02,
    1.336213e-02,
    1.420611e-02,
    1.510486e-02,
    1.609220e-02,
    1.716439e-02,
    1.831930e-02,
    1.961309e-02,
    2.100770e-02,
    2.254767e-02,
    2.423026e-02,
    2.604985e-02,
    2.803414e-02,
    3.023304e-02,
    3.261836e-02,
    3.527531e-02,
    3.818206e-02,
    4.134298e-02,
    4.483305e-02,
    4.870456e-02,
    5.296474e-02,
    5.765964e-02,
    6.281384e-02,
    6.850944e-02,
    7.481781e-02,
    8.173895e-02,
    8.937667e-02,
    9.774769e-02,
    1.069522e-01,
    1.170782e-01,
    1.282548e-01,
    1.405027e-01,
    1.539673e-01,
    1.686191e-01,
    1.846574e-01,
    2.020771e-01,
    2.210459e-01,
    2.415374e-01,
    2.637095e-01,
    2.872893e-01,
    3.124426e-01,
    3.391945e-01,
    3.673232e-01,
    3.966229e-01,
    4.269345e-01,
    4.580085e-01,
    4.895437e-01,
    5.211873e-01,
    5.525534e-01,
    5.832051e-01,
    6.129409e-01,
    6.414552e-01,
    6.685366e-01,
    6.941538e-01,
    7.180993e-01,
    7.403272e-01,
    7.608830e-01,
    7.798100e-01,
    7.972520e-01,
    8.132446e-01,
    8.277589e-01,
    8.410031e-01,
    8.531792e-01,
    8.642644e-01,
    8.743145e-01,
    8.835475e-01,
    8.919562e-01,
    8.996102e-01,
    9.066574e-01,
    9.130520e-01,
    9.189736e-01,
    9.244128e-01,
    9.293928e-01,
    9.339722e-01,
    9.381913e-01,
    9.421451e-01,
    9.457657e-01,
    9.491230e-01,
    9.522394e-01,
    9.551589e-01,
    9.578458e-01,
    9.603525e-01,
    9.626876e-01,
    9.648771e-01,
    9.669361e-01,
    9.688645e-01,
    9.706878e-01,
    9.723798e-01,
    9.739466e-01,
    9.754310e-01,
    9.768135e-01,
    9.781311e-01,
    9.793539e-01,
    9.805091e-01,
    9.816204e-01,
    9.826387e-01,
    9.835868e-01,
    9.844858e-01,
    9.853414e-01,
    9.861397e-01,
    9.869005e-01,
    9.876234e-01,
    9.882862e-01,
    9.889344e-01,
    9.895357e-01,
    9.901063e-01,
    9.906435e-01,
    9.911373e-01,
    9.916119e-01,
    9.920576e-01,
    9.924838e-01,
    9.928730e-01,
    9.932470e-01,
    9.935931e-01,
    9.939352e-01,
    9.942430e-01,
    9.945419e-01,
    9.948185e-01,
    9.950818e-01,
    9.953442e-01,
    9.955881e-01,
    9.958246e-01,
    9.960354e-01,
    9.962523e-01,
    9.964441e-01,
    9.966252e-01,
    9.967933e-01,
    9.969595e-01,
    9.971159e-01,
    9.972685e-01,
    9.974102e-01,
    9.975482e-01,
    9.976761e-01,
    9.977899e-01,
    9.978979e-01,
    9.979996e-01,
    9.980933e-01,
    9.981938e-01,
    9.982869e-01,
    9.983711e-01,
    9.984526e-01,
    9.985349e-01,
    9.986089e-01,
    9.986793e-01,
    9.987437e-01,
    9.988027e-01,
    9.988587e-01,
    9.989131e-01,
    9.989656e-01,
    9.990185e-01,
    9.990667e-01,
    9.991120e-01,
    9.991549e-01,
    9.991993e-01,
    9.992385e-01,
    9.992718e-01,
    9.993088e-01,
    9.993393e-01,
    9.993723e-01,
    9.993982e-01,
    9.994278e-01,
    9.994549e-01,
    9.994825e-01,
    9.995078e-01,
    9.995336e-01,
    9.995536e-01,
    9.995732e-01,
    9.995946e-01,
    9.996155e-01,
    9.996358e-01,
    9.996535e-01,
    9.996697e-01,
    9.996842e-01,
    9.996983e-01,
    9.997115e-01,
    9.997263e-01,
    9.997404e-01,
    9.997520e-01,
    9.997612e-01,
    9.997721e-01,
    9.997834e-01,
    9.997917e-01,
    9.998010e-01,
    9.998097e-01,
    9.998197e-01,
    9.998272e-01,
    9.998347e-01,
    9.998427e-01,
    9.998504e-01,
    9.998564e-01,
    9.998629e-01,
    9.998688e-01,
    9.998745e-01,
    9.998795e-01,
    9.998834e-01,
    9.998882e-01,
    9.998937e-01,
    9.998978e-01,
    9.999021e-01,
    9.999067e-01,
    9.999094e-01,
    9.999134e-01,
    9.999174e-01,
    9.999209e-01,
    9.999240e-01,
    9.999275e-01,
    9.999312e-01,
    9.999336e-01,
    9.999374e-01,
    9.999390e-01,
    9.999415e-01,
    9.999435e-01,
    9.999452e-01,
    9.999480e-01,
    9.999506e-01,
    9.999517e-01,
    9.999541e-01,
    9.999555e-01,
    9.999570e-01,
    9.999581e-01,
    9.999595e-01,
    9.999608e-01,
    9.999623e-01,
    9.999636e-01,
    9.999655e-01,
    9.999665e-01,
    9.999680e-01,
    9.999691e-01,
    9.999702e-01,
    9.999712e-01,
    9.999718e-01,
    9.999733e-01,
    9.999739e-01,
    9.999751e-01,
    9.999761e-01,
    9.999771e-01,
    9.999782e-01,
    9.999792e-01,
    9.999801e-01,
    9.999807e-01,
    9.999813e-01,
    9.999816e-01,
    9.999826e-01,
    9.999831e-01,
    9.999834e-01,
    9.999836e-01,
    9.999838e-01,
    9.999841e-01,
    9.999846e-01,
    9.999854e-01,
    9.999856e-01,
    9.999859e-01,
    9.999870e-01,
    9.999874e-01,
    9.999882e-01,
    9.999885e-01,
    9.999890e-01,
    9.999892e-01,
    9.999892e-01,
    9.999893e-01,
    9.999895e-01,
    9.999900e-01,
    9.999904e-01,
    9.999907e-01,
    9.999909e-01,
    9.999914e-01,
    9.999917e-01,
    9.999923e-01,
    9.999923e-01,
    9.999925e-01,
    9.999928e-01,
    9.999932e-01,
    9.999934e-01,
    9.999937e-01,
    9.999939e-01,
    9.999942e-01,
    9.999948e-01,
    9.999950e-01,
    9.999950e-01,
    9.999952e-01,
    9.999954e-01,
    9.999955e-01,
    9.999957e-01,
    9.999957e-01,
    9.999959e-01,
    9.999959e-01,
    9.999960e-01,
    9.999960e-01,
    9.999962e-01,
    9.999963e-01,
    9.999965e-01,
    9.999965e-01,
    9.999967e-01,
    9.999967e-01,
    9.999968e-01,
    9.999973e-01,
    9.999973e-01,
    9.999974e-01,
    9.999975e-01,
    9.999977e-01,
    9.999977e-01,
    9.999977e-01,
    9.999977e-01,
    9.999979e-01,
    9.999979e-01,
    9.999979e-01,
    9.999982e-01,
    9.999983e-01,
    9.999984e-01,
    9.999985e-01,
    9.999985e-01,
    9.999985e-01,
    9.999985e-01,
    9.999985e-01,
    9.999985e-01,
    9.999985e-01,
    9.999985e-01,
    9.999987e-01,
    9.999987e-01,
    9.999988e-01,
    9.999988e-01,
    9.999988e-01,
    9.999988e-01,
    9.999989e-01,
    9.999989e-01,
    9.999989e-01,
    9.999989e-01,
    9.999992e-01,
];

// NCBI reference: ncbi-blast/c++/include/algo/blast/composition_adjustment/composition_adjustment.h:49-58
// ```c
// typedef struct Blast_AminoAcidComposition {
//     double prob[COMPO_LARGEST_ALPHABET];
//     int numTrueAminoAcids;
// } Blast_AminoAcidComposition;
// ```
#[derive(Debug, Clone, PartialEq)]
pub struct BlastAminoAcidComposition {
    pub prob: [f64; COMPO_LARGEST_ALPHABET],
    pub num_true_amino_acids: i32,
}

impl BlastAminoAcidComposition {
    pub fn empty() -> Self {
        Self {
            prob: [0.0; COMPO_LARGEST_ALPHABET],
            num_true_amino_acids: 0,
        }
    }
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/composition_adjustment/composition_adjustment.h:76-86
// ```c
// typedef struct Blast_MatrixInfo {
//     char * matrixName;
//     int    **startMatrix;
//     double **startFreqRatios;
//     int      rows;
//     int      cols;
//     int      positionBased;
//     double   ungappedLambda;
// } Blast_MatrixInfo;
// ```
#[derive(Debug, Clone, PartialEq)]
pub struct BlastMatrixInfo {
    pub matrix: ScoringMatrix,
    pub position_based: bool,
    pub cols: usize,
    pub ungapped_lambda: f64,
    pub start_matrix: [[i32; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET],
    pub start_freq_ratios: [[f64; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET],
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/composition_adjustment/composition_adjustment.h:109-117
// ```c
// typedef struct Blast_CompositionWorkspace {
//     double ** mat_b;
//     double ** mat_final;
//     double * first_standard_freq;
//     double * second_standard_freq;
// } Blast_CompositionWorkspace;
// ```
#[derive(Debug, Clone)]
pub(crate) struct BlastCompositionWorkspace {
    mat_b: [[f64; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA],
    mat_final: [[f64; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA],
    first_standard_freq: [f64; COMPO_NUM_TRUE_AA],
    second_standard_freq: [f64; COMPO_NUM_TRUE_AA],
}

impl BlastCompositionWorkspace {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1271-1332
    // ```c
    // Blast_CompositionWorkspace * Blast_CompositionWorkspaceNew(void)
    // int Blast_CompositionWorkspaceInit(Blast_CompositionWorkspace * NRrecord,
    //                                    const char *matrixName)
    // {
    //     if (0 == Blast_GetJointProbsForMatrix(NRrecord->mat_b,
    //                                           NRrecord->first_standard_freq,
    //                                           NRrecord->second_standard_freq,
    //                                           matrixName)) {
    //         return 0;
    //     }
    // }
    // ```
    pub(crate) fn new_blosum62() -> Self {
        Self {
            mat_b: BLOSUM62_JOINT_PROBS,
            mat_final: [[0.0; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA],
            first_standard_freq: BLOSUM62_BG,
            second_standard_freq: BLOSUM62_BG,
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1025-1151
// ```c
// s_ScalePSSM(...); s_ScaleSquareMatrix(...);
// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdjustedProteinMatrix {
    pub scores: [[i32; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET],
}

impl AdjustedProteinMatrix {
    #[inline]
    pub fn score(&self, query_residue: u8, subject_residue: u8) -> i32 {
        self.scores[query_residue as usize][subject_residue as usize]
    }

    #[inline]
    pub fn is_positive(&self, query_residue: u8, subject_residue: u8) -> bool {
        self.score(query_residue, subject_residue) > 0
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1414-1530
// ```c
// int Blast_AdjustScores(..., EMatrixAdjustRule *matrix_adjust_rule, ...,
//                        double *pvalueForThisPair, int compositionTestIndex,
//                        double *ratioToPassBack)
// ```
#[derive(Debug, Clone)]
pub struct BlastAdjustScoresResult {
    pub matrix_adjust_rule: EMatrixAdjustRule,
    pub adjusted_matrix: AdjustedProteinMatrix,
    pub pvalue_for_this_pair: Option<f64>,
    pub lambda_ratio: f64,
}

// NCBI c++/src/algo/blast/composition_adjustment/composition_adjustment.c:61,1204-1248:
// static const int kCompositionMargin = 20;
// left = start; for (i = left; i > 0; i--) if (subject_data[i-1] == eStopChar)
//     { if (i + kCompositionMargin < left) left = i + kCompositionMargin; break; }
// if (i == 0) left = 0;
// right = finish; for (i = right; i < length; i++)
//     if (subject_data[i] == eStopChar)
//     { if (i - kCompositionMargin > right) right = i - kCompositionMargin; break; }
// if (i == length) right = length;
#[allow(dead_code)] // TBLASTN calls this after its translated-range callback is ported.
pub(crate) fn blast_get_composition_range(
    subject_data: &[u8],
    start: usize,
    finish: usize,
) -> Result<(usize, usize)> {
    if start > finish || finish > subject_data.len() {
        anyhow::bail!("invalid NCBI translated composition interval");
    }
    let mut left = 0;
    for i in (1..=start).rev() {
        if subject_data[i - 1] == 25 {
            left = if i + 20 < start { i + 20 } else { start };
            break;
        }
    }
    let mut right = subject_data.len();
    for i in finish..subject_data.len() {
        if subject_data[i] == 25 {
            right = if i.saturating_sub(20) > finish {
                i - 20
            } else {
                finish
            };
            break;
        }
    }
    Ok((left, right))
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1170-1194
// ```c
// void Blast_ReadAaComposition(Blast_AminoAcidComposition * composition,
//                              int alphsize, const Uint1 * sequence, int length)
// {
//     ...
// }
// ```
pub fn read_aa_composition(sequence: &[u8]) -> BlastAminoAcidComposition {
    let mut composition = BlastAminoAcidComposition::empty();
    for &residue in sequence {
        let residue_index = residue as usize;
        if residue_index < COMPO_LARGEST_ALPHABET
            && (ALPHA_CONVERT[residue_index] >= 0 || residue == 24)
        {
            composition.prob[residue_index] += 1.0;
            composition.num_true_amino_acids += 1;
        }
    }

    if composition.prob[24] > 0.0 {
        composition.prob[3] += composition.prob[24];
        composition.prob[24] = 0.0;
    }

    if composition.num_true_amino_acids > 0 {
        let denom = composition.num_true_amino_acids as f64;
        for prob in &mut composition.prob {
            *prob /= denom;
        }
    }
    composition
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2216-2231
// ```c
// self->ungappedLambda = sbp->kbp_ideal->Lambda / scale_factor;
// status = s_GetStartFreqRatios(self->startFreqRatios, matrixName);
// Blast_Int4MatrixFromFreq(self->startMatrix, self->cols,
//                          self->startFreqRatios, self->ungappedLambda);
// ```
pub fn build_matrix_info(matrix: ScoringMatrix, ungapped_lambda: f64) -> Result<BlastMatrixInfo> {
    if matrix != ScoringMatrix::Blosum62 {
        bail!("composition adjustment matrix info is only ported for BLOSUM62");
    }

    let start_freq_ratios = BLOSUM62_START_FREQ_RATIOS;
    let mut start_matrix = [[0i32; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET];
    int4_matrix_from_freq(
        &mut start_matrix,
        COMPO_LARGEST_ALPHABET,
        &start_freq_ratios,
        ungapped_lambda,
    );

    Ok(BlastMatrixInfo {
        matrix,
        position_based: false,
        cols: COMPO_LARGEST_ALPHABET,
        ungapped_lambda,
        start_matrix,
        start_freq_ratios,
    })
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:31-48
// ```c
// double ** Nlm_DenseMatrixNew(int nrows, int ncols)
// ```
fn dense_matrix_new(nrows: usize, ncols: usize) -> Vec<Vec<f64>> {
    vec![vec![0.0; ncols]; nrows]
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:52-74
// ```c
// double ** Nlm_LtriangMatrixNew(int n)
// ```
fn ltriang_matrix_new(n: usize) -> Vec<Vec<f64>> {
    (0..n).map(|row| vec![0.0; row + 1]).collect()
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:120-143
// ```c
// void Nlm_FactorLtriangPosDef(double ** A, int n)
// ```
fn factor_ltriang_pos_def(a: &mut [Vec<f64>], n: usize) {
    for i in 0..n {
        for j in 0..i {
            let mut temp = a[i][j];
            for k in 0..j {
                temp -= a[i][k] * a[j][k];
            }
            a[i][j] = temp / a[j][j];
        }
        let mut temp = a[i][i];
        for k in 0..i {
            temp -= a[i][k] * a[i][k];
        }
        a[i][i] = temp.sqrt();
    }
}

// EXPERIMENT (LOSAT_X_COMPFAST / LOSAT_X_COMPSHADOW): 0 = original loops, 1 = reordered loops,
// 2 = reordered loops checked bit-for-bit against the original on every call.
// No NCBI counterpart: reads LOSAT_X_COMPFAST and LOSAT_X_COMPSHADOW once. It chooses, for each function below that has a reordered variant, between the original loops, the reordered loops, or both with a bitwise comparison; it does not change any value NCBI computes.
pub(crate) fn x_comp_mode() -> u8 {
    use std::sync::OnceLock;
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_COMPSHADOW").is_some() {
            2
        } else if std::env::var_os("LOSAT_X_COMPFAST").is_some() {
            1
        } else {
            0
        }
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:725-727
// ```c
//     n  = alphsize * alphsize;
//     mA = 2 * alphsize - 1;
//     m  = constrain_rel_entropy ? mA + 1 : mA;
// ```
// The stride is `m` for `alphsize == 20` and `constrain_rel_entropy` (the size of `W`).
/// Column stride of the column-major copy of the Newton matrix.
const X_WC_STRIDE: usize = 2 * COMPO_NUM_TRUE_AA;

// EXPERIMENT: `Nlm_FactorLtriangPosDef` with the loops interchanged.
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:120-143
// ```c
// for (i = 0;  i < n;  i++) {
//     for (j = 0;  j < i;  j++) {
//         temp = A[i][j];
//         for (k = 0;  k < j;  k++) {
//             temp -= A[i][k] * A[j][k];
//         }
//         A[i][j] = temp/A[j][j];
//     }
//     temp = A[i][i];
//     for (k = 0;  k < i;  k++) {
//         temp -= A[i][k] * A[i][k];
//     }
//     A[i][i] = sqrt(temp);
// }
// ```
// Every element (i, j) still receives exactly the same floating-point
// operations in the same order: the subtractions for k = 0, 1, ..., j-1 (each
// one `temp -= A[i][k] * A[j][k]` with the final values of column k), then the
// division by A[j][j] (or the square root on the diagonal). Only the order in
// which *different* elements are advanced changes: after column k is final it
// is subtracted from all later columns, whose elements are independent of each
// other, so the innermost loop has no loop-carried dependency.
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:143-156
// ```c
//     for (i = 0;  i < n;  i++) {
//         for (j = 0;  j < i;  j++) {
//             temp = A[i][j];
//             for (k = 0;  k < j;  k++) {
//                 temp -= A[i][k] * A[j][k];
//             }
//             A[i][j] = temp/A[j][j];
//         }
//         temp = A[i][i];
//         for (k = 0;  k < i;  k++) {
//             temp -= A[i][k] * A[i][k];
//         }
//         A[i][i] = sqrt(temp);
//     }
// ```
// (The line numbers in the older comment above are for a different revision; these are the pinned ones.)
// The loop nest is interchanged: column k is finished and then subtracted from every later
// column. Element (i, j) still receives its subtractions for k = 0..j-1 in order, then the division
// or square root.
fn x_factor_ltriang_pos_def_cols(a: &mut [Vec<f64>], n: usize, wc: &mut [f64]) {
    let s = X_WC_STRIDE;
    debug_assert!(n <= s && wc.len() >= s * s && s % 2 == 0);
    // Column-major copy. Slots above a column's diagonal and below row n are
    // scratch: they are written by the padded loops below but never read as
    // matrix elements.
    for k in 0..n {
        let col = &mut wc[k * s..(k + 1) * s];
        col[..k].fill(0.0);
        for (i, slot) in col[k..n].iter_mut().enumerate() {
            *slot = a[k + i][k];
        }
        col[n..].fill(0.0);
    }
    let end = (n + 1) & !1;
    for k in 0..n {
        let (done, rest) = wc.split_at_mut((k + 1) * s);
        let ck = &mut done[k * s..(k + 1) * s];
        let d = ck[k].sqrt();
        ck[k] = d;
        for value in &mut ck[k + 1..n] {
            *value /= d;
        }
        let ck = &*ck;
        // A multiplier that is +0.0 subtracts `y * 0.0` from its whole target
        // column.  For a finite `y` that product is +0.0, or -0.0 when `y` has
        // its sign bit set; `x - (+0.0)` is `x` for every `x` (NaN, the
        // infinities and both zeros included) and `x - (-0.0)` is `x` except
        // that it turns -0.0 into +0.0.  So such a column is left alone, apart
        // from that one case in the few rows whose `y` is negative.  The Newton
        // matrix `A D A^T` has these multipliers throughout its leading diagonal
        // block (optimize_target_freq.c:119-147 only ever adds to `W[j][j]`,
        // `W[i + alphsize - 1][j]` and the second diagonal); its negative `y`
        // are in the relative-entropy row.
        const NEG_ROWS: usize = 4;
        let mut neg_rows = [0usize; NEG_ROWS];
        let mut neg_count = 0usize;
        // Integer tests on the bit patterns (all exponent bits set = infinity
        // or NaN; the top bit = negative), without branches.
        const EXPONENT: u64 = 0x7FF0_0000_0000_0000;
        let mut not_finite = 0u32;
        let mut signs = 0u64;
        for &v in &ck[k + 1..n] {
            let bits = v.to_bits();
            not_finite += u32::from(bits & EXPONENT == EXPONENT);
            signs |= bits;
        }
        let finite = not_finite == 0;
        if finite && (signs >> 63) != 0 {
            for (offset, &v) in ck[k + 1..n].iter().enumerate() {
                if v.is_sign_negative() {
                    if neg_count < NEG_ROWS {
                        neg_rows[neg_count] = k + 1 + offset;
                    }
                    neg_count += 1;
                }
            }
        }
        let plain = finite && neg_count <= NEG_ROWS;
        // The effect of the skipped update on column `column` (see above).
        let settle = |rest: &mut [f64], column: usize| {
            for &row in &neg_rows[..neg_count] {
                if row >= column {
                    let cell = &mut rest[(column - k - 1) * s + row];
                    if cell.to_bits() == (-0.0f64).to_bits() {
                        *cell = 0.0;
                    }
                }
            }
        };
        // Trailing update, four columns at a time. Rows start at an even
        // index <= j so every column pair of rows is handled by one 2-lane
        // step; the extra rows are scratch slots of those columns.
        let mut j = k + 1;
        while j + 4 <= n {
            let (s0, s1, s2, s3) = (ck[j], ck[j + 1], ck[j + 2], ck[j + 3]);
            if plain && (s0.to_bits() | s1.to_bits() | s2.to_bits() | s3.to_bits()) == 0 {
                if neg_count != 0 {
                    for column in j..j + 4 {
                        settle(rest, column);
                    }
                }
                j += 4;
                continue;
            }
            let i0 = j & !1;
            let len = end - i0;
            let block = &mut rest[(j - k - 1) * s..(j - k + 3) * s];
            let (c0, block) = block.split_at_mut(s);
            let (c1, block) = block.split_at_mut(s);
            let (c2, c3) = block.split_at_mut(s);
            let y = &ck[i0..i0 + len];
            let c0 = &mut c0[i0..i0 + len];
            let c1 = &mut c1[i0..i0 + len];
            let c2 = &mut c2[i0..i0 + len];
            let c3 = &mut c3[i0..i0 + len];
            for i in 0..len {
                let yi = y[i];
                c0[i] -= yi * s0;
                c1[i] -= yi * s1;
                c2[i] -= yi * s2;
                c3[i] -= yi * s3;
            }
            j += 4;
        }
        while j < n {
            let sj = ck[j];
            if plain && sj.to_bits() == 0 {
                if neg_count != 0 {
                    settle(rest, j);
                }
                j += 1;
                continue;
            }
            let i0 = j & !1;
            let len = end - i0;
            let cj = &mut rest[(j - k - 1) * s + i0..(j - k - 1) * s + i0 + len];
            for (x, &y) in cj.iter_mut().zip(&ck[i0..i0 + len]) {
                *x -= y * sj;
            }
            j += 1;
        }
    }
    for (i, row) in a.iter_mut().enumerate().take(n) {
        for (k, value) in row.iter_mut().enumerate().take(i + 1) {
            *value = wc[k * s + i];
        }
    }
}

// EXPERIMENT: `Nlm_SolveLtriangPosDef`, forward substitution by columns.
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:147-171
// ```c
// for (i = 0;  i < n;  i++) {
//     temp = x[i];
//     for (j = 0;  j < i;  j++) {
//         temp -= L[i][j] * x[j];
//     }
//     x[i] = temp/L[i][i];
// }
// for (j = n - 1;  j >= 0;  j--) {
//     x[j] /= L[j][j];
//     for (i = 0;  i < j;  i++) {
//         x[i] -= L[j][i] * x[j];
//     }
// }
// ```
// x[i] still receives `-= L[i][j] * x[j]` for j = 0, ..., i-1 in that order
// and is then divided by L[i][i]; the second loop is unchanged.
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:171-184
// ```c
//     for (i = 0;  i < n;  i++) {
//         temp = x[i];
//         for (j = 0;  j < i;  j++) {
//             temp -= L[i][j] * x[j];
//         }
//         x[i] = temp/L[i][i];
//     }
//     /* Now x = z.  Back solve the system L\T y = z */
//     for (j = n - 1;  j >= 0;  j--) {
//         x[j] /= L[j][j];
//         for (i = 0;  i < j;  i++) {
//             x[i] -= L[j][i] * x[j];
//         }
//     }
// ```
// Forward substitution is done by columns; the back substitution is the reference loop.
fn x_solve_ltriang_pos_def_cols(x: &mut [f64], n: usize, l: &[Vec<f64>], wc: &[f64]) {
    let s = X_WC_STRIDE;
    for j in 0..n {
        let cj = &wc[j * s..j * s + n];
        let xj = x[j] / cj[j];
        x[j] = xj;
        for (xi, &lij) in x[j + 1..n].iter_mut().zip(&cj[j + 1..n]) {
            *xi -= lij * xj;
        }
    }
    for j in (0..n).rev() {
        let row = &l[j];
        x[j] /= row[j];
        let xj = x[j];
        for (xi, &lji) in x[..j].iter_mut().zip(&row[..j]) {
            *xi -= lji * xj;
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:147-171
// ```c
// void Nlm_SolveLtriangPosDef(double x[], int n, double ** L)
// ```
fn solve_ltriang_pos_def(x: &mut [f64], n: usize, l: &[Vec<f64>]) {
    for i in 0..n {
        let mut temp = x[i];
        for (j, &xj) in x.iter().enumerate().take(i) {
            temp -= l[i][j] * xj;
        }
        x[i] = temp / l[i][i];
    }
    for j in (0..n).rev() {
        x[j] /= l[j][j];
        for i in 0..j {
            x[i] -= l[j][i] * x[j];
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:175-193
// ```c
// double Nlm_EuclideanNorm(const double v[], int n)
// ```
fn euclidean_norm(v: &[f64]) -> f64 {
    let mut sum = 1.0;
    let mut scale = 0.0;
    for &value in v {
        if value != 0.0 {
            let abs_value = value.abs();
            if scale < abs_value {
                sum = 1.0 + sum * (scale / abs_value) * (scale / abs_value);
                scale = abs_value;
            } else {
                sum += (abs_value / scale) * (abs_value / scale);
            }
        }
    }
    scale * sum.sqrt()
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:197-205
// ```c
// void Nlm_AddVectors(double y[], int n, double alpha, const double x[])
// ```
fn add_vectors(y: &mut [f64], alpha: f64, x: &[f64]) {
    for (left, right) in y.iter_mut().zip(x.iter()) {
        *left += alpha * *right;
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:209-223
// ```c
// double Nlm_StepBound(const double x[], int n, const double step_x[], double max)
// ```
fn step_bound(x: &[f64], step_x: &[f64], max: f64) -> f64 {
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:224-238
    // ```c
    // double
    // Nlm_StepBound(const double x[], int n, const double step_x[], double max)
    // {
    //     int i;                 /* iteration index */
    //     double alpha = max;    /* current largest permitted step */
    //     for (i = 0; i < n; i++) {
    //         double alpha_i;    /* a step to the boundary for the current i */
    //         alpha_i = -x[i] / step_x[i];
    //         if (alpha_i >= 0 && alpha_i < alpha) {
    //             alpha = alpha_i;
    //         }
    //     }
    //     return alpha;
    // ```
    // Dispatch point of LOSAT_X_COMPFAST / LOSAT_X_COMPSHADOW: the reference path is
    // `step_bound_reference`, the port of this function.
    let mode = x_comp_mode();
    if mode != 0 && x.len() <= 400 && x.len() == step_x.len() {
        let fast = x_step_bound(x, step_x, max);
        if mode == 2 {
            let reference = step_bound_reference(x, step_x, max);
            assert!(
                fast.to_bits() == reference.to_bits(),
                "x_step_bound differs from the reference"
            );
        }
        return fast;
    }
    step_bound_reference(x, step_x, max)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:230-238
// ```c
//     for (i = 0; i < n; i++) {
//         double alpha_i;    /* a step to the boundary for the current i */
//         alpha_i = -x[i] / step_x[i];
//         if (alpha_i >= 0 && alpha_i < alpha) {
//             alpha = alpha_i;
//         }
//     }
//     return alpha;
// ```
// The quotients are computed first and then scanned in the original order with the same test.
// EXPERIMENT: `Nlm_StepBound`. The quotients do not depend on each other;
// compute them first, then scan them in the original order.
fn x_step_bound(x: &[f64], step_x: &[f64], max: f64) -> f64 {
    let mut ratio = [0.0f64; 400];
    for ((out, &x_i), &step_i) in ratio.iter_mut().zip(x).zip(step_x) {
        *out = -x_i / step_i;
    }
    let mut alpha = max;
    for &alpha_i in &ratio[..x.len()] {
        if alpha_i >= 0.0 && alpha_i < alpha {
            alpha = alpha_i;
        }
    }
    alpha
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:224-238
// ```c
// double
// Nlm_StepBound(const double x[], int n, const double step_x[], double max)
// {
//     int i;                 /* iteration index */
//     double alpha = max;    /* current largest permitted step */
//     for (i = 0; i < n; i++) {
//         double alpha_i;    /* a step to the boundary for the current i */
//         alpha_i = -x[i] / step_x[i];
//         if (alpha_i >= 0 && alpha_i < alpha) {
//             alpha = alpha_i;
//         }
//     }
//     return alpha;
// ```
// The original `step_bound`, renamed; it is the port of this function.
fn step_bound_reference(x: &[f64], step_x: &[f64], max: f64) -> f64 {
    let mut alpha = max;
    for (&x_i, &step_i) in x.iter().zip(step_x.iter()) {
        let alpha_i = -x_i / step_i;
        if alpha_i >= 0.0 && alpha_i < alpha {
            alpha = alpha_i;
        }
    }
    alpha
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:329-337
// ```c
// double Blast_MatrixEntropy(double ** matrix, int alphsize, const double row_prob[],
//                            const double col_prob[], double Lambda)
// ```
fn matrix_entropy(
    matrix: &[Vec<f64>],
    alphsize: usize,
    row_prob: &[f64; COMPO_NUM_TRUE_AA],
    col_prob: &[f64; COMPO_NUM_TRUE_AA],
    lambda: f64,
) -> f64 {
    let mut entropy = 0.0;
    for i in 0..alphsize {
        for j in 0..alphsize {
            let nat_score = lambda * matrix[i][j];
            entropy += nat_score.exp() * nat_score * row_prob[i] * col_prob[j];
        }
    }
    entropy
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:358-376
// ```c
// double Blast_TargetFreqEntropy(double ** target_freq)
// ```
fn target_freq_entropy(target_freq: &[[f64; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA]) -> f64 {
    let mut row_prob = [0.0; COMPO_NUM_TRUE_AA];
    let mut col_prob = [0.0; COMPO_NUM_TRUE_AA];
    for (i, row) in target_freq.iter().enumerate() {
        for (j, &freq) in row.iter().enumerate() {
            row_prob[i] += freq;
            col_prob[j] += freq;
        }
    }
    let mut entropy = 0.0;
    for (i, row) in target_freq.iter().enumerate() {
        for (j, &freq) in row.iter().enumerate() {
            entropy += freq * (freq / row_prob[i] / col_prob[j]).ln();
        }
    }
    entropy
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:385-397
// ```c
// void Blast_CalcFreqRatios(double ** ratios, int alphsize,
//                           double row_prob[], double col_prob[])
// ```
fn calc_freq_ratios_dense(
    ratios: &mut [Vec<f64>],
    alphsize: usize,
    row_prob: &[f64],
    col_prob: &[f64],
) {
    for i in 0..alphsize {
        if row_prob[i] > 0.0 {
            for j in 0..alphsize {
                if col_prob[j] > 0.0 {
                    ratios[i][j] /= row_prob[i] * col_prob[j];
                }
            }
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:401-415
// ```c
// void Blast_FreqRatioToScore(double ** matrix, int rows, int cols, double Lambda)
// ```
fn freq_ratio_to_score_dense(matrix: &mut [Vec<f64>], rows: usize, cols: usize, lambda: f64) {
    for row in matrix.iter_mut().take(rows) {
        for value in row.iter_mut().take(cols) {
            if *value == 0.0 {
                *value = COMPO_SCORE_MIN as f64;
            } else {
                *value = value.ln() / lambda;
            }
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:401-415
// ```c
// void Blast_FreqRatioToScore(double ** matrix, int rows, int cols, double Lambda)
// ```
fn freq_ratio_to_score_square(
    matrix: &mut [[f64; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET],
    rows: usize,
    cols: usize,
    lambda: f64,
) {
    for row in matrix.iter_mut().take(rows) {
        for value in row.iter_mut().take(cols) {
            if *value == 0.0 {
                *value = COMPO_SCORE_MIN as f64;
            } else {
                *value = value.ln() / lambda;
            }
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:117-240
// ```c
// void Blast_CalcLambdaFullPrecision(double * plambda, int *piterations,
//                                    double **score, int alphsize,
//                                    const double row_prob[],
//                                    const double col_prob[],
//                                    double lambda_tolerance,
//                                    double function_tolerance,
//                                    int max_iterations)
// ```
fn calc_lambda_full_precision_from_scores(
    score: &[Vec<f64>],
    alphsize: usize,
    row_prob: &[f64; COMPO_NUM_TRUE_AA],
    col_prob: &[f64; COMPO_NUM_TRUE_AA],
) -> (f64, usize) {
    let mut f = 4.0;
    let mut left = 0.0;
    let mut right = 1.0;
    let mut x = 0.367879441171f64;
    let mut is_newton = false;
    let mut max_score = COMPO_SCORE_MIN as f64;
    let mut avg_score = 0.0;

    for i in 0..alphsize {
        if row_prob[i] == 0.0 {
            continue;
        }
        for (j, &col_p) in col_prob.iter().enumerate().take(alphsize) {
            if col_p == 0.0 {
                continue;
            }
            if max_score < score[i][j] {
                max_score = score[i][j];
            }
            avg_score += row_prob[i] * col_p * score[i][j];
        }
    }
    if max_score <= 0.0 || avg_score >= 0.0 {
        return (-1.0, K_LAMBDA_ITERATION_LIMIT);
    }

    let mut k = 0usize;
    while k < K_LAMBDA_ITERATION_LIMIT {
        let lambda = -x.ln();
        let was_newton = is_newton;
        let x_pow_max_score = (-max_score * lambda).exp();
        f = -x_pow_max_score;
        let mut slope = max_score * f / x;
        let fold = f;

        for i in 0..alphsize {
            if row_prob[i] == 0.0 {
                continue;
            }
            for (j, &col_p) in col_prob.iter().enumerate().take(alphsize) {
                if col_p == 0.0 {
                    continue;
                }
                let ff = if max_score != score[i][j] {
                    let diff_score = max_score - score[i][j];
                    let ff = row_prob[i] * col_p * (-lambda * diff_score).exp();
                    slope += diff_score * ff / x;
                    ff
                } else {
                    row_prob[i] * col_p
                };
                f += ff;
            }
        }

        if f > 0.0 {
            left = x;
        } else if f < 0.0 {
            right = x;
        } else {
            break;
        }
        if right - left <= 2.0 * left * (1.0 - right) * K_LAMBDA_ERROR_TOLERANCE
            && (f / x_pow_max_score).abs() <= K_LAMBDA_FUNCTION_TOLERANCE
        {
            x = (left + right) / 2.0;
            break;
        }
        // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:290-320
        // ```c
        // if ((was_newton && fabs(f) > .9 * fabs(fold)) || slope >= 0) {
        //     x = (left + right)/2;
        // } else {
        //     ...
        //     if (y <= left || y >= right) {
        //         x = (left + right)/2;
        //     } else {
        //         is_newton = 1;
        //         x = y;
        //     }
        // }
        // ```
        if (was_newton && f.abs() > 0.9 * fold.abs()) || slope >= 0.0 {
            x = (left + right) / 2.0;
        } else {
            let p = -f / slope;
            let y = x + p;
            if y <= left || y >= right {
                x = (left + right) / 2.0;
            } else {
                is_newton = true;
                x = y;
                if p.abs() <= K_LAMBDA_ERROR_TOLERANCE * x * (1.0 - x)
                    && (f / x_pow_max_score).abs() <= K_LAMBDA_FUNCTION_TOLERANCE
                {
                    break;
                }
            }
        }
        k += 1;
    }

    if k < K_LAMBDA_ITERATION_LIMIT {
        (-x.ln(), k)
    } else {
        (-1.0, k)
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:136-147
// ```c
// double Blast_GetRelativeEntropy(const double A[], const double B[])
// {
//     ...
//     return sqrt(value);
// }
// ```
pub fn blast_get_relative_entropy(
    a: &[f64; COMPO_NUM_TRUE_AA],
    b: &[f64; COMPO_NUM_TRUE_AA],
) -> f64 {
    let mut value = 0.0;
    for i in 0..COMPO_NUM_TRUE_AA {
        let temp = (a[i] + b[i]) / 2.0;
        if temp > 0.0 {
            if a[i] > 0.0 {
                value += a[i] * (a[i] / temp).ln() / 2.0;
            }
            if b[i] > 0.0 {
                value += b[i] * (b[i] / temp).ln() / 2.0;
            }
        }
    }
    if value < 0.0 {
        0.0
    } else {
        value.sqrt()
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:89-107
// ```c
// void Blast_ApplyPseudocounts(double * probs20, int number_of_observations,
//                              const double * background_probs20, int pseudocounts)
// ```
fn apply_pseudocounts(
    probs20: &mut [f64; COMPO_NUM_TRUE_AA],
    number_of_observations: i32,
    background_probs20: &[f64; COMPO_NUM_TRUE_AA],
    pseudocounts: i32,
) {
    let mut sum = 0.0;
    for &value in probs20.iter() {
        sum += value;
    }
    if sum == 0.0 {
        sum = 1.0;
    }
    let weight = pseudocounts as f64 / (number_of_observations as f64 + pseudocounts as f64);
    for i in 0..COMPO_NUM_TRUE_AA {
        probs20[i] = (1.0 - weight) * probs20[i] / sum + weight * background_probs20[i];
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:428-441
// ```c
// static void s_GatherLetterProbs(double * outputLetterProbs,
//                                 const double * inputLetterProbs, int alphsize)
// ```
fn gather_letter_probs(
    output_letter_probs: &mut [f64; COMPO_NUM_TRUE_AA],
    input_letter_probs: &[f64; COMPO_LARGEST_ALPHABET],
) {
    for c in 0..COMPO_LARGEST_ALPHABET {
        if ALPHA_CONVERT[c] >= 0 {
            output_letter_probs[ALPHA_CONVERT[c] as usize] = input_letter_probs[c];
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:444-458
// ```c
// static void s_UnpackLetterProbs(double std_probs[], int alphsize, const double probs[])
// ```
fn unpack_letter_probs(
    std_probs: &mut [f64; COMPO_LARGEST_ALPHABET],
    probs: &[f64; COMPO_NUM_TRUE_AA],
) {
    for c in 0..COMPO_LARGEST_ALPHABET {
        if ALPHA_CONVERT[c] >= 0 {
            std_probs[c] = probs[ALPHA_CONVERT[c] as usize];
        } else {
            std_probs[c] = 0.0;
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:467-478
// ```c
// static void s_SetPairAmbigProbsToSum(double probs[], int alphsize)
// ```
fn set_pair_ambig_probs_to_sum(probs: &mut [f64; COMPO_LARGEST_ALPHABET]) {
    probs[E_BCHAR] = probs[E_DCHAR] + probs[E_NCHAR];
    probs[E_ZCHAR] = probs[E_ECHAR] + probs[E_QCHAR];
    probs[E_JCHAR] = probs[E_ICHAR] + probs[E_LCHAR];
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:485-536
// ```c
// int Blast_EntropyOldFreqNewContext(...)
// ```
fn entropy_old_freq_new_context(
    target_freq: &[[f64; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA],
    row_prob: &[f64; COMPO_NUM_TRUE_AA],
    col_prob: &[f64; COMPO_NUM_TRUE_AA],
) -> Result<Option<f64>> {
    let mut scores = dense_matrix_new(COMPO_NUM_TRUE_AA, COMPO_NUM_TRUE_AA);
    let mut old_col_prob = [0.0; COMPO_NUM_TRUE_AA];
    let mut old_row_prob = [0.0; COMPO_NUM_TRUE_AA];
    for i in 0..COMPO_NUM_TRUE_AA {
        for j in 0..COMPO_NUM_TRUE_AA {
            old_row_prob[i] += target_freq[i][j];
            old_col_prob[j] += target_freq[i][j];
            scores[i][j] = target_freq[i][j];
        }
    }
    calc_freq_ratios_dense(&mut scores, COMPO_NUM_TRUE_AA, &old_row_prob, &old_col_prob);
    freq_ratio_to_score_dense(&mut scores, COMPO_NUM_TRUE_AA, COMPO_NUM_TRUE_AA, 1.0);
    let (lambda, iter_count) =
        calc_lambda_full_precision_from_scores(&scores, COMPO_NUM_TRUE_AA, row_prob, col_prob);
    if iter_count < K_LAMBDA_ITERATION_LIMIT && lambda > 0.0 {
        Ok(Some(matrix_entropy(
            &scores,
            COMPO_NUM_TRUE_AA,
            row_prob,
            col_prob,
            lambda,
        )))
    } else {
        Ok(None)
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:538-604
// ```c
// void Blast_TrueAaToStdTargetFreqs(double ** StdFreq, int StdAlphsize, double ** freq)
// ```
fn true_aa_to_std_target_freqs(
    std_freq: &mut [[f64; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET],
    freq: &[[f64; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA],
) {
    let mut sum = 0.0;
    for row in freq {
        for &value in row {
            sum += value;
        }
    }
    for a in 0..COMPO_LARGEST_ALPHABET {
        if ALPHA_CONVERT[a] < 0 {
            for b in 0..COMPO_LARGEST_ALPHABET {
                std_freq[a][b] = 0.0;
            }
        } else {
            let small_a = ALPHA_CONVERT[a] as usize;
            for b in 0..COMPO_LARGEST_ALPHABET {
                if ALPHA_CONVERT[b] < 0 {
                    std_freq[a][b] = 0.0;
                } else {
                    let small_b = ALPHA_CONVERT[b] as usize;
                    std_freq[a][b] = freq[small_a][small_b] / sum;
                }
            }
            std_freq[a][E_BCHAR] = std_freq[a][E_DCHAR] + std_freq[a][E_NCHAR];
            std_freq[a][E_ZCHAR] = std_freq[a][E_ECHAR] + std_freq[a][E_QCHAR];
            std_freq[a][E_JCHAR] = std_freq[a][E_ICHAR] + std_freq[a][E_LCHAR];
        }
    }
    for b in 0..COMPO_LARGEST_ALPHABET {
        std_freq[E_BCHAR][b] = std_freq[E_DCHAR][b] + std_freq[E_NCHAR][b];
        std_freq[E_ZCHAR][b] = std_freq[E_ECHAR][b] + std_freq[E_QCHAR][b];
        std_freq[E_JCHAR][b] = std_freq[E_ICHAR][b] + std_freq[E_LCHAR][b];
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:607-639
// ```c
// static double s_CalcAvgScore(...)
// static double s_CalcXScore(...)
// ```
fn calc_avg_score(
    scores: &[[f64; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET],
    row: usize,
    col_probs: &[f64; COMPO_LARGEST_ALPHABET],
) -> f64 {
    let mut score = 0.0;
    for col in 0..COMPO_LARGEST_ALPHABET {
        if ALPHA_CONVERT[col] >= 0 {
            score += scores[row][col] * col_probs[col];
        }
    }
    score
}

fn calc_x_score(
    scores: &[[f64; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET],
    col: usize,
    row_probs: &[f64; COMPO_LARGEST_ALPHABET],
) -> f64 {
    let mut score = 0.0;
    for row in 0..COMPO_LARGEST_ALPHABET {
        if ALPHA_CONVERT[row] >= 0 {
            score += scores[row][col] * row_probs[row];
        }
    }
    score.min(K_MAXIMUM_X_SCORE)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:708-711
// ```c
// static long Nint(double x)
// {
//     x += (x >= 0. ? 0.5 : -0.5);
// ```
fn nint(x: f64) -> i32 {
    (x + if x >= 0.0 { 0.5 } else { -0.5 }) as i32
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:648-690
// ```c
// static void s_SetXUOScores(double ** M, int alphsize, ...)
// ```
fn set_xuo_scores(
    scores: &mut [[f64; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET],
    row_probs: &[f64; COMPO_LARGEST_ALPHABET],
    col_probs: &[f64; COMPO_LARGEST_ALPHABET],
) {
    let mut score_xx = 0.0;
    for i in 0..COMPO_LARGEST_ALPHABET {
        if ALPHA_CONVERT[i] >= 0 {
            let avg_ix = calc_avg_score(scores, i, col_probs);
            scores[i][E_XCHAR] = avg_ix.min(K_MAXIMUM_X_SCORE);
            score_xx += avg_ix * row_probs[i];
            scores[E_XCHAR][i] = calc_x_score(scores, i, row_probs);
        }
    }
    scores[E_XCHAR][E_XCHAR] = score_xx.min(K_MAXIMUM_X_SCORE);
    scores[E_BCHAR][E_XCHAR] = calc_avg_score(scores, E_BCHAR, col_probs).min(K_MAXIMUM_X_SCORE);
    scores[E_XCHAR][E_BCHAR] = calc_x_score(scores, E_BCHAR, row_probs);
    scores[E_ZCHAR][E_XCHAR] = calc_avg_score(scores, E_ZCHAR, col_probs).min(K_MAXIMUM_X_SCORE);
    scores[E_XCHAR][E_ZCHAR] = calc_x_score(scores, E_ZCHAR, row_probs);
    scores[E_JCHAR][E_XCHAR] = calc_avg_score(scores, E_JCHAR, col_probs).min(K_MAXIMUM_X_SCORE);
    scores[E_XCHAR][E_JCHAR] = calc_x_score(scores, E_JCHAR, row_probs);

    scores[E_SELENOCYSTEINE] = scores[E_CCHAR];
    for row in scores.iter_mut() {
        row[E_SELENOCYSTEINE] = row[E_CCHAR];
    }
    scores[E_OCHAR] = scores[E_XCHAR];
    for row in scores.iter_mut() {
        row[E_OCHAR] = row[E_XCHAR];
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:701-734
// ```c
// static long Nint(double x)
// static void s_RoundScoreMatrix(int **matrix, int rows, int cols, double **floatScoreMatrix)
// ```
fn round_score_matrix(
    scores: &[[f64; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET],
) -> [[i32; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET] {
    let mut matrix = [[0i32; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET];
    for row in 0..COMPO_LARGEST_ALPHABET {
        for col in 0..COMPO_LARGEST_ALPHABET {
            matrix[row][col] = if scores[row][col] < i32::MIN as f64 {
                i32::MIN
            } else {
                nint(scores[row][col])
            };
        }
    }
    matrix
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:945-958
// ```c
// void Blast_Int4MatrixFromFreq(int **matrix, int alphsize,
//                               double ** freq, double Lambda)
// {
//     ...
//     Blast_FreqRatioToScore(dMatrix, 1, alphsize, Lambda);
//     s_RoundScoreMatrix(&matrix[i], 1, alphsize, dMatrix);
// }
// ```
fn int4_matrix_from_freq(
    matrix: &mut [[i32; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET],
    alphsize: usize,
    freq: &[[f64; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET],
    lambda: f64,
) {
    let mut scores = *freq;
    freq_ratio_to_score_square(&mut scores, alphsize, alphsize, lambda);
    let rounded = round_score_matrix(&scores);
    for row in 0..alphsize {
        matrix[row][..alphsize].copy_from_slice(&rounded[row][..alphsize]);
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:760-793
// ```c
// static int s_ScoresStdAlphabet(int ** Matrix, int Alphsize, double ** target_freq,
//                                int ** StartMatrix, const double row_prob[],
//                                const double col_prob[], double Lambda)
// ```
fn scores_std_alphabet(
    target_freq: &[[f64; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA],
    matrix_info: &BlastMatrixInfo,
    row_prob: &[f64; COMPO_NUM_TRUE_AA],
    col_prob: &[f64; COMPO_NUM_TRUE_AA],
) -> AdjustedProteinMatrix {
    let scores = scores_std_alphabet_unrounded(target_freq, matrix_info, row_prob, col_prob);
    scores_std_alphabet_round(&scores, matrix_info)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:779-788
// ```c
//     s_UnpackLetterProbs(RowProb, Alphsize, row_prob);
//     s_SetPairAmbigProbsToSum(RowProb, Alphsize);
//     s_UnpackLetterProbs(ColProb, Alphsize, col_prob);
//     s_SetPairAmbigProbsToSum(ColProb, Alphsize);
//     Blast_TrueAaToStdTargetFreqs(Scores, Alphsize, target_freq);
//     Blast_CalcFreqRatios(Scores, Alphsize, RowProb, ColProb);
//     Blast_FreqRatioToScore(Scores, Alphsize, Alphsize, Lambda);
//     s_SetXUOScores(Scores, Alphsize, RowProb, ColProb);
// ```
// This function ports lines 779-788 of `s_ScoresStdAlphabet`; `scores_std_alphabet_round` ports lines 790-796.
/// The part of `s_ScoresStdAlphabet` before the scores are rounded.
fn scores_std_alphabet_unrounded(
    target_freq: &[[f64; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA],
    matrix_info: &BlastMatrixInfo,
    row_prob: &[f64; COMPO_NUM_TRUE_AA],
    col_prob: &[f64; COMPO_NUM_TRUE_AA],
) -> [[f64; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET] {
    let mut row_prob_std = [0.0; COMPO_LARGEST_ALPHABET];
    let mut col_prob_std = [0.0; COMPO_LARGEST_ALPHABET];
    unpack_letter_probs(&mut row_prob_std, row_prob);
    set_pair_ambig_probs_to_sum(&mut row_prob_std);
    unpack_letter_probs(&mut col_prob_std, col_prob);
    set_pair_ambig_probs_to_sum(&mut col_prob_std);

    let mut scores = [[0.0; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET];
    true_aa_to_std_target_freqs(&mut scores, target_freq);
    calc_freq_ratios(&mut scores, &row_prob_std, &col_prob_std);
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:403-413
    // ```c
    // Blast_FreqRatioToScore(double ** matrix, int rows, int cols, double Lambda)
    // {
    //     int i;
    //     for (i = 0;  i < rows;  i++) {
    //         int j;
    //         for (j = 0;  j < cols;  j++) {
    //             if (0.0 == matrix[i][j]) {
    //                 matrix[i][j] = COMPO_SCORE_MIN;
    //             } else {
    //                 matrix[i][j] = log(matrix[i][j])/Lambda;
    //             }
    // ```
    // `Blast_FreqRatioToScore`: every non-zero entry becomes `log(entry) / Lambda`. The logs of all non-zero
    // entries are computed first in one slice call, then divided by Lambda in the same entry order.
    // ln of every non-zero ratio (through x_logclone::log_slice: libm unless
    // LOSAT_X_LOGCLONE is set), then / Lambda, as Blast_FreqRatioToScore does.
    let mut ratios = [0.0f64; COMPO_LARGEST_ALPHABET * COMPO_LARGEST_ALPHABET];
    let mut logs = [0.0f64; COMPO_LARGEST_ALPHABET * COMPO_LARGEST_ALPHABET];
    let mut count = 0usize;
    for row in &scores {
        for &value in row.iter() {
            if value != 0.0 {
                ratios[count] = value;
                count += 1;
            }
        }
    }
    crate::utils::x_logclone::log_slice(&ratios[..count], &mut logs[..count]);
    let mut next = 0usize;
    for row in &mut scores {
        for value in row.iter_mut() {
            if *value == 0.0 {
                *value = COMPO_SCORE_MIN as f64;
            } else {
                *value = logs[next] / matrix_info.ungapped_lambda;
                next += 1;
            }
        }
    }
    set_xuo_scores(&mut scores, &row_prob_std, &col_prob_std);
    scores
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:790-796
// ```c
//     s_RoundScoreMatrix(Matrix, Alphsize, Alphsize, Scores);
//     Nlm_DenseMatrixFree(&Scores);
//     for (i = 0;  i < Alphsize;  i++) {
//         Matrix[i][eStopChar] = StartMatrix[i][eStopChar];
//         Matrix[eStopChar][i] = StartMatrix[eStopChar][i];
//     }
// ```
/// The rounding and the stop-codon scores of `s_ScoresStdAlphabet`.
fn scores_std_alphabet_round(
    scores: &[[f64; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET],
    matrix_info: &BlastMatrixInfo,
) -> AdjustedProteinMatrix {
    let mut rounded = round_score_matrix(scores);
    for index in 0..COMPO_LARGEST_ALPHABET {
        rounded[index][E_STOP_CHAR] = matrix_info.start_matrix[index][E_STOP_CHAR];
        rounded[E_STOP_CHAR][index] = matrix_info.start_matrix[E_STOP_CHAR][index];
    }
    AdjustedProteinMatrix { scores: rounded }
}

// No NCBI counterpart: reads LOSAT_X_ADJMEMO_STATS once; it does not change any value NCBI computes.
fn x_adjmemo_enabled() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_X_ADJMEMO_STATS").is_some())
}
// No NCBI counterpart: call counter and set of hashed inputs for the LOSAT_X_ADJMEMO_STATS report; it does not change any value NCBI computes.
static X_ADJMEMO: std::sync::Mutex<Option<(u64, std::collections::HashSet<u64>)>> =
    std::sync::Mutex::new(None);
// No NCBI counterpart: counts calls of `Blast_AdjustScores` that reach the matrix adjustment, and distinct inputs; it does not change any value NCBI computes.
fn x_adjmemo_record(
    query_composition: &BlastAminoAcidComposition,
    subject_composition: &BlastAminoAcidComposition,
    rule: EMatrixAdjustRule,
) {
    if !x_adjmemo_enabled() {
        return;
    }
    use std::hash::{Hash, Hasher};
    let mut h = std::collections::hash_map::DefaultHasher::new();
    for v in query_composition.prob.iter() {
        v.to_bits().hash(&mut h);
    }
    query_composition.num_true_amino_acids.hash(&mut h);
    for v in subject_composition.prob.iter() {
        v.to_bits().hash(&mut h);
    }
    subject_composition.num_true_amino_acids.hash(&mut h);
    (rule as u8).hash(&mut h);
    let key = h.finish();
    let mut guard = X_ADJMEMO.lock().expect("adjmemo lock");
    let entry = guard.get_or_insert_with(|| (0, std::collections::HashSet::new()));
    entry.0 += 1;
    entry.1.insert(key);
}
// No NCBI counterpart: counter of the shadow comparisons; it does not change any value NCBI computes.
/// Calls compared with the reference by LOSAT_X_NEWTONEXACTSHADOW.
pub(crate) static X_NEWTONEXACT_SHADOWED: std::sync::atomic::AtomicU64 =
    std::sync::atomic::AtomicU64::new(0);

// No NCBI counterpart: prints the shadow counter; it does not change any value NCBI computes.
/// Shadow-mode summary (printed by `main` at exit).
pub fn x_newtonexact_print_stats() {
    if super::x_newton_exact::mode() == 2 {
        eprintln!(
            "[X_SHADOW] NEWTONEXACT calls_compared={} (status and all 400 target frequencies bit-identical)",
            X_NEWTONEXACT_SHADOWED.load(std::sync::atomic::Ordering::Relaxed)
        );
    }
    // No NCBI counterpart: LOSAT_X_NEWTONLANES shadow and LOSAT_X_STATS counters; it does not change any
    // value NCBI computes.
    super::x_newton_lanes::print_stats();
}

// No NCBI counterpart: prints the LOSAT_X_ADJMEMO_STATS counters; it does not change any value NCBI computes.
pub fn x_adjmemo_print_stats() {
    if !x_adjmemo_enabled() {
        return;
    }
    let guard = X_ADJMEMO.lock().expect("adjmemo lock");
    if let Some((calls, set)) = guard.as_ref() {
        eprintln!(
            "[X_STATS] ADJMEMO matrix_adjust_calls={} distinct_inputs={}",
            calls,
            set.len()
        );
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:95-126
// ```c
// static void ScaledSymmetricProductA(double ** W, const double diagonal[], int alphsize)
// ```
fn scaled_symmetric_product_a(w: &mut [Vec<f64>], diagonal: &[f64], alphsize: usize) {
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:126-137
    // ```c
    //     for (i = 0;  i < alphsize;  i++) {
    //         for (j = 0;  j < alphsize;  j++) {
    //             double dd;     /* an individual diagonal element */
    //             dd = diagonal[i * alphsize + j];
    //             W[j][j] += dd;
    //             if (i > 0) {
    //                 W[i + alphsize - 1][j] += dd;
    //                 W[i + alphsize - 1][i + alphsize - 1] += dd;
    //             }
    //         }
    // ```
    // Dispatch point of LOSAT_X_COMPFAST / LOSAT_X_COMPSHADOW: the reference path is
    // `scaled_symmetric_product_a_reference`, the port of this function.
    match x_comp_mode() {
        1 if alphsize <= COMPO_NUM_TRUE_AA => {
            return x_scaled_symmetric_product_a(w, diagonal, alphsize)
        }
        2 if alphsize <= COMPO_NUM_TRUE_AA => {
            let mut check: Vec<Vec<f64>> = w.to_vec();
            x_scaled_symmetric_product_a(&mut check, diagonal, alphsize);
            scaled_symmetric_product_a_reference(w, diagonal, alphsize);
            let m = 2 * alphsize - 1;
            for i in 0..m {
                for k in 0..=i {
                    assert!(
                        check[i][k].to_bits() == w[i][k].to_bits(),
                        "x_scaled_symmetric_product_a differs at ({i},{k})"
                    );
                }
            }
            return;
        }
        _ => {}
    }
    scaled_symmetric_product_a_reference(w, diagonal, alphsize)
}

// EXPERIMENT: `ScaledSymmetricProductA` with register accumulators.
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:108-125
// ```c
// for (i = 0;  i < alphsize;  i++) {
//     for (j = 0;  j < alphsize;  j++) {
//         double dd;
//         dd = diagonal[i * alphsize + j];
//         W[j][j] += dd;
//         if (i > 0) {
//             W[i + alphsize - 1][j] += dd;
//             W[i + alphsize - 1][i + alphsize - 1] += dd;
//         }
//     }
// }
// ```
// W[j][j] still accumulates dd for i = 0, 1, ... in order starting from 0.0,
// W[i+alphsize-1][i+alphsize-1] accumulates for j = 0, 1, ... in order starting
// from 0.0, and every W[i+alphsize-1][j] is the single sum 0.0 + dd.
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:126-137
// ```c
//     for (i = 0;  i < alphsize;  i++) {
//         for (j = 0;  j < alphsize;  j++) {
//             double dd;     /* an individual diagonal element */
//             dd = diagonal[i * alphsize + j];
//             W[j][j] += dd;
//             if (i > 0) {
//                 W[i + alphsize - 1][j] += dd;
//                 W[i + alphsize - 1][i + alphsize - 1] += dd;
//             }
//         }
// ```
// (The line numbers in the older comment above are for a different revision; these are the pinned ones.)
fn x_scaled_symmetric_product_a(w: &mut [Vec<f64>], diagonal: &[f64], alphsize: usize) {
    let n = alphsize;
    let m = 2 * n - 1;
    for (row_idx, row) in w.iter_mut().enumerate().take(m) {
        for value in row.iter_mut().take(row_idx + 1) {
            *value = 0.0;
        }
    }
    let mut col_sum = [0.0f64; COMPO_NUM_TRUE_AA];
    for i in 0..n {
        let row = &diagonal[i * n..(i + 1) * n];
        for (acc, &dd) in col_sum[..n].iter_mut().zip(row) {
            *acc += dd;
        }
    }
    for j in 0..n {
        w[j][j] = col_sum[j];
    }
    for i in 1..n {
        let row = &diagonal[i * n..(i + 1) * n];
        let w_row = &mut w[i + n - 1];
        let mut acc = 0.0f64;
        for (slot, &dd) in w_row[..n].iter_mut().zip(row) {
            *slot = 0.0 + dd;
            acc += dd;
        }
        w_row[i + n - 1] = acc;
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:112-113,121-125
// ```c
// static void
// ScaledSymmetricProductA(double ** W, const double diagonal[], int alphsize)
// ...
//     for (rowW = 0;  rowW < m;  rowW++) {
//         for (colW = 0;  colW <= rowW;  colW++) {
//             W[rowW][colW] = 0.0;
//         }
//     }
// ```
// The original `scaled_symmetric_product_a`, renamed; it is the port of this function.
fn scaled_symmetric_product_a_reference(w: &mut [Vec<f64>], diagonal: &[f64], alphsize: usize) {
    let m = 2 * alphsize - 1;
    for (row_idx, row) in w.iter_mut().enumerate().take(m) {
        for value in row.iter_mut().take(row_idx + 1) {
            *value = 0.0;
        }
    }
    for i in 0..alphsize {
        for j in 0..alphsize {
            let dd = diagonal[i * alphsize + j];
            w[j][j] += dd;
            if i > 0 {
                w[i + alphsize - 1][j] += dd;
                w[i + alphsize - 1][i + alphsize - 1] += dd;
            }
        }
    }
}

// EXPERIMENT: `MultiplyByA` with the same additions per element in the same
// order, arranged so that independent elements advance together.
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:155-164
// ```c
// for (i = 0;  i < alphsize;  i++) {
//     for (j = 0;  j < alphsize;  j++) {
//         y[j] += alpha * x[i * alphsize + j];
//     }
// }
// for (i = 1;  i < alphsize;  i++) {
//     for (j = 0;  j < alphsize;  j++) {
//         y[i + alphsize - 1] += alpha * x[i * alphsize + j];
//     }
// }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:170-179
// ```c
//     for (i = 0;  i < alphsize;  i++) {
//         for (j = 0;  j < alphsize;  j++) {
//             y[j] += alpha * x[i * alphsize + j];
//         }
//     }
//     for (i = 1;  i < alphsize;  i++) {
//         for (j = 0;  j < alphsize;  j++) {
//             y[i + alphsize - 1] += alpha * x[i * alphsize + j];
//         }
//     }
// ```
// (The line numbers in the older comment above are for a different revision; these are the pinned ones.)
fn x_multiply_by_a(beta: f64, y: &mut [f64], alphsize: usize, alpha: f64, x: &[f64]) {
    if beta == 0.0 {
        for value in y.iter_mut() {
            *value = 0.0;
        }
    } else if beta != 1.0 {
        for value in y.iter_mut() {
            *value *= beta;
        }
    }
    let n = alphsize;
    {
        let head = &mut y[..n];
        for i in 0..n {
            let row = &x[i * n..(i + 1) * n];
            for (yj, &xv) in head.iter_mut().zip(row) {
                *yj += alpha * xv;
            }
        }
    }
    // y[i + n - 1] accumulates x[i*n + j] for j = 0..n in order; four rows at
    // a time, each with its own accumulator.
    let mut i = 1;
    while i + 4 <= n {
        let r0 = &x[i * n..(i + 1) * n];
        let r1 = &x[(i + 1) * n..(i + 2) * n];
        let r2 = &x[(i + 2) * n..(i + 3) * n];
        let r3 = &x[(i + 3) * n..(i + 4) * n];
        let (mut a0, mut a1, mut a2, mut a3) = (y[i + n - 1], y[i + n], y[i + n + 1], y[i + n + 2]);
        for j in 0..n {
            a0 += alpha * r0[j];
            a1 += alpha * r1[j];
            a2 += alpha * r2[j];
            a3 += alpha * r3[j];
        }
        y[i + n - 1] = a0;
        y[i + n] = a1;
        y[i + n + 1] = a2;
        y[i + n + 2] = a3;
        i += 4;
    }
    while i < n {
        let row = &x[i * n..(i + 1) * n];
        let mut acc = y[i + n - 1];
        for &xv in row {
            acc += alpha * xv;
        }
        y[i + n - 1] = acc;
        i += 1;
    }
}

// EXPERIMENT: `MultiplyByAtranspose`, same two additions per element.
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:195-204
// ```c
// for (i = 0;  i < alphsize;  i++) {
//     for (j = 0;  j < alphsize;  j++) {
//         k = i * alphsize + j;
//         y[k] += alpha * x[j];
//         if (i > 0) {
//             y[k] += alpha * x[i + alphsize - 1];
//         }
//     }
// }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:214-223
// ```c
//     for (i = 0;  i < alphsize;  i++) {
//         for (j = 0;  j < alphsize;  j++) {
//             k = i * alphsize + j;
//             y[k] += alpha * x[j];
//             if (i > 0) {
//                 y[k] += alpha * x[i + alphsize - 1];
//             }
//         }
//     }
// ```
// (The line numbers in the older comment above are for a different revision; these are the pinned ones.)
fn x_multiply_by_at(beta: f64, y: &mut [f64], alphsize: usize, alpha: f64, x: &[f64]) {
    if beta == 0.0 {
        for value in y.iter_mut() {
            *value = 0.0;
        }
    } else if beta != 1.0 {
        for value in y.iter_mut() {
            *value *= beta;
        }
    }
    let n = alphsize;
    let head = &x[..n];
    for i in 0..n {
        let row = &mut y[i * n..(i + 1) * n];
        if i > 0 {
            let xi = x[i + n - 1];
            for (yk, &xj) in row.iter_mut().zip(head) {
                let t = *yk + alpha * xj;
                *yk = t + alpha * xi;
            }
        } else {
            for (yk, &xj) in row.iter_mut().zip(head) {
                *yk += alpha * xj;
            }
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:138-165
// ```c
// static void MultiplyByA(double beta, double y[], int alphsize, double alpha, const double x[])
// ```
fn multiply_by_a(beta: f64, y: &mut [f64], alphsize: usize, alpha: f64, x: &[f64]) {
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:170-179
    // ```c
    //     for (i = 0;  i < alphsize;  i++) {
    //         for (j = 0;  j < alphsize;  j++) {
    //             y[j] += alpha * x[i * alphsize + j];
    //         }
    //     }
    //     for (i = 1;  i < alphsize;  i++) {
    //         for (j = 0;  j < alphsize;  j++) {
    //             y[i + alphsize - 1] += alpha * x[i * alphsize + j];
    //         }
    //     }
    // ```
    // Dispatch point of LOSAT_X_COMPFAST / LOSAT_X_COMPSHADOW: the reference path is
    // `multiply_by_a_reference`, the port of this function.
    match x_comp_mode() {
        1 => return x_multiply_by_a(beta, y, alphsize, alpha, x),
        2 => {
            let mut check = y.to_vec();
            x_multiply_by_a(beta, &mut check, alphsize, alpha, x);
            multiply_by_a_reference(beta, y, alphsize, alpha, x);
            assert!(
                check
                    .iter()
                    .zip(y.iter())
                    .all(|(a, b)| a.to_bits() == b.to_bits()),
                "x_multiply_by_a differs from the reference"
            );
            return;
        }
        _ => {}
    }
    multiply_by_a_reference(beta, y, alphsize, alpha, x)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:154-157,170-179
// ```c
// static void
// MultiplyByA(double beta, double y[], int alphsize,
//             double alpha, const double x[])
// {
// ...
//     for (i = 0;  i < alphsize;  i++) {
//         for (j = 0;  j < alphsize;  j++) {
//             y[j] += alpha * x[i * alphsize + j];
//         }
//     }
//     for (i = 1;  i < alphsize;  i++) {
//         for (j = 0;  j < alphsize;  j++) {
//             y[i + alphsize - 1] += alpha * x[i * alphsize + j];
//         }
//     }
// ```
// The original `multiply_by_a`, renamed; it is the port of this function.
fn multiply_by_a_reference(beta: f64, y: &mut [f64], alphsize: usize, alpha: f64, x: &[f64]) {
    if beta == 0.0 {
        for value in y.iter_mut() {
            *value = 0.0;
        }
    } else if beta != 1.0 {
        for value in y.iter_mut() {
            *value *= beta;
        }
    }
    for i in 0..alphsize {
        for j in 0..alphsize {
            y[j] += alpha * x[i * alphsize + j];
        }
    }
    for i in 1..alphsize {
        for j in 0..alphsize {
            y[i + alphsize - 1] += alpha * x[i * alphsize + j];
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:176-205
// ```c
// static void MultiplyByAtranspose(double beta, double y[], int alphsize,
//                                  double alpha, const double x[])
// ```
fn multiply_by_at(beta: f64, y: &mut [f64], alphsize: usize, alpha: f64, x: &[f64]) {
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:214-223
    // ```c
    //     for (i = 0;  i < alphsize;  i++) {
    //         for (j = 0;  j < alphsize;  j++) {
    //             k = i * alphsize + j;
    //             y[k] += alpha * x[j];
    //             if (i > 0) {
    //                 y[k] += alpha * x[i + alphsize - 1];
    //             }
    //         }
    //     }
    // ```
    // Dispatch point of LOSAT_X_COMPFAST / LOSAT_X_COMPSHADOW: the reference path is
    // `multiply_by_at_reference`, the port of this function.
    match x_comp_mode() {
        1 => return x_multiply_by_at(beta, y, alphsize, alpha, x),
        2 => {
            let mut check = y.to_vec();
            x_multiply_by_at(beta, &mut check, alphsize, alpha, x);
            multiply_by_at_reference(beta, y, alphsize, alpha, x);
            assert!(
                check
                    .iter()
                    .zip(y.iter())
                    .all(|(a, b)| a.to_bits() == b.to_bits()),
                "x_multiply_by_at differs from the reference"
            );
            return;
        }
        _ => {}
    }
    multiply_by_at_reference(beta, y, alphsize, alpha, x)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:195-197,214-223
// ```c
// static void
// MultiplyByAtranspose(double beta, double y[], int alphsize,
//                      double alpha, const double x[])
// ...
//     for (i = 0;  i < alphsize;  i++) {
//         for (j = 0;  j < alphsize;  j++) {
//             k = i * alphsize + j;
//             y[k] += alpha * x[j];
//             if (i > 0) {
//                 y[k] += alpha * x[i + alphsize - 1];
//             }
//         }
//     }
// ```
// The original `multiply_by_at`, renamed; it is the port of this function.
fn multiply_by_at_reference(beta: f64, y: &mut [f64], alphsize: usize, alpha: f64, x: &[f64]) {
    if beta == 0.0 {
        for value in y.iter_mut() {
            *value = 0.0;
        }
    } else if beta != 1.0 {
        for value in y.iter_mut() {
            *value *= beta;
        }
    }
    for i in 0..alphsize {
        for j in 0..alphsize {
            let k = i * alphsize + j;
            y[k] += alpha * x[j];
            if i > 0 {
                y[k] += alpha * x[i + alphsize - 1];
            }
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:216-231
// ```c
// static void ResidualsLinearConstraints(...)
// ```
fn residuals_linear_constraints(
    ra: &mut [f64],
    alphsize: usize,
    x: &[f64],
    row_sums: &[f64; COMPO_NUM_TRUE_AA],
    col_sums: &[f64; COMPO_NUM_TRUE_AA],
) {
    ra[..alphsize].copy_from_slice(&col_sums[..alphsize]);
    for i in 1..alphsize {
        ra[i + alphsize - 1] = row_sums[i];
    }
    multiply_by_a(1.0, ra, alphsize, -1.0, x);
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:241-262
// ```c
// static void DualResiduals(double resids_x[], int alphsize, double ** grads,
//                           const double z[], int constrain_rel_entropy)
// ```
fn dual_residuals(
    resids_x: &mut [f64],
    alphsize: usize,
    grads: &[Vec<f64>],
    z: &[f64],
    constrain_rel_entropy: bool,
) {
    if constrain_rel_entropy {
        let eta = z[2 * alphsize - 1];
        for i in 0..resids_x.len() {
            resids_x[i] = -grads[0][i] + eta * grads[1][i];
        }
    } else {
        for i in 0..resids_x.len() {
            resids_x[i] = -grads[0][i];
        }
    }
    multiply_by_at(1.0, resids_x, alphsize, 1.0, z);
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:282-323
// ```c
// static void CalculateResiduals(...)
// ```
fn calculate_residuals(
    resids_x: &mut [f64],
    alphsize: usize,
    resids_z: &mut [f64],
    values: &[f64; 2],
    grads: &[Vec<f64>],
    row_sums: &[f64; COMPO_NUM_TRUE_AA],
    col_sums: &[f64; COMPO_NUM_TRUE_AA],
    x: &[f64],
    z: &[f64],
    constrain_rel_entropy: bool,
    relative_entropy: f64,
) -> f64 {
    dual_residuals(resids_x, alphsize, grads, z, constrain_rel_entropy);
    let norm_resids_x = euclidean_norm(resids_x);
    residuals_linear_constraints(resids_z, alphsize, x, row_sums, col_sums);
    let norm_resids_z = if constrain_rel_entropy {
        resids_z[2 * alphsize - 1] = relative_entropy - values[1];
        euclidean_norm(&resids_z[..(2 * alphsize)])
    } else {
        euclidean_norm(&resids_z[..(2 * alphsize - 1)])
    };
    (norm_resids_x * norm_resids_x + norm_resids_z * norm_resids_z).sqrt()
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:336-344
// ```c
// typedef struct ReNewtonSystem {
//     ...
// } ReNewtonSystem;
// ```
struct ReNewtonSystem {
    alphsize: usize,
    constrain_rel_entropy: bool,
    // EXPERIMENT: column-major copy of the factor (see x_factor_ltriang_pos_def_cols).
    // No NCBI counterpart: a column-major copy of the lower triangle of `W` after the factorisation, and a flag that says whether it is current; it does not change any value NCBI computes.
    x_wc: Vec<f64>,
    x_wc_valid: bool,
    w: Vec<Vec<f64>>,
    dinv: Vec<f64>,
    grad_re: Vec<f64>,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:376-394
// ```c
// static ReNewtonSystem * ReNewtonSystemNew(int alphsize)
// ```
fn re_newton_system_new(alphsize: usize) -> ReNewtonSystem {
    ReNewtonSystem {
        alphsize,
        constrain_rel_entropy: true,
        // No NCBI counterpart: the copy is allocated only when LOSAT_X_COMPFAST or LOSAT_X_COMPSHADOW is set; it does not change any value NCBI computes.
        x_wc: if x_comp_mode() != 0 {
            vec![0.0; X_WC_STRIDE * X_WC_STRIDE]
        } else {
            Vec::new()
        },
        x_wc_valid: false,
        w: ltriang_matrix_new(2 * alphsize),
        dinv: vec![0.0; alphsize * alphsize],
        grad_re: vec![0.0; alphsize * alphsize],
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:408-470
// ```c
// static void FactorReNewtonSystem(...)
// ```
fn factor_re_newton_system(
    newton_system: &mut ReNewtonSystem,
    x: &[f64],
    z: &[f64],
    grads: &[Vec<f64>],
    constrain_rel_entropy: bool,
    workspace: &mut [f64],
) {
    let alphsize = newton_system.alphsize;
    let n = alphsize * alphsize;
    let m = if constrain_rel_entropy {
        2 * alphsize
    } else {
        2 * alphsize - 1
    };

    newton_system.constrain_rel_entropy = constrain_rel_entropy;
    if constrain_rel_entropy {
        let eta = z[m - 1];
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:507-510
        // ```c
        //         eta = z[m - 1];
        //         for (i = 0;  i < n;  i++) {
        //             Dinv[i] = x[i] / (1 - eta);
        //         }
        // ```
        // Dispatch point of LOSAT_X_COMPFAST: the reference path is the loop in the `else` branch. The
        // divisor `1 - eta` is the same value for every element, so computing it once gives the same quotients.
        if x_comp_mode() == 1 {
            // Same quotient per element; the divisor is the same value each time.
            let divisor = 1.0 - eta;
            for (out, &xi) in newton_system.dinv[..n].iter_mut().zip(&x[..n]) {
                *out = xi / divisor;
            }
        } else {
            for i in 0..n {
                newton_system.dinv[i] = x[i] / (1.0 - eta);
            }
        }
    } else {
        newton_system.dinv.copy_from_slice(x);
    }

    scaled_symmetric_product_a(&mut newton_system.w, &newton_system.dinv, alphsize);
    if constrain_rel_entropy {
        newton_system.grad_re.copy_from_slice(&grads[1]);
        newton_system.w[m - 1][m - 1] = 0.0;
        for i in 0..n {
            workspace[i] = newton_system.dinv[i] * newton_system.grad_re[i];
            newton_system.w[m - 1][m - 1] += newton_system.grad_re[i] * workspace[i];
        }
        multiply_by_a(
            0.0,
            &mut newton_system.w[m - 1][..m - 1],
            alphsize,
            1.0,
            workspace,
        );
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:534-534
    // ```c
    //     Nlm_FactorLtriangPosDef(W, m);
    // ```
    // Dispatch point of LOSAT_X_COMPFAST / LOSAT_X_COMPSHADOW: the reference path is the
    // `factor_ltriang_pos_def` call at the end of this function.
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:143-149
    // ```c
    //     for (i = 0;  i < n;  i++) {
    //         for (j = 0;  j < i;  j++) {
    //             temp = A[i][j];
    //             for (k = 0;  k < j;  k++) {
    //                 temp -= A[i][k] * A[j][k];
    //             }
    //             A[i][j] = temp/A[j][j];
    // ```
    newton_system.x_wc_valid = false;
    if m <= X_WC_STRIDE {
        match x_comp_mode() {
            1 => {
                x_factor_ltriang_pos_def_cols(&mut newton_system.w, m, &mut newton_system.x_wc);
                newton_system.x_wc_valid = true;
                return;
            }
            2 => {
                let mut check = newton_system.w.clone();
                x_factor_ltriang_pos_def_cols(&mut check, m, &mut newton_system.x_wc);
                factor_ltriang_pos_def(&mut newton_system.w, m);
                for i in 0..m {
                    for k in 0..=i {
                        assert!(
                            check[i][k].to_bits() == newton_system.w[i][k].to_bits(),
                            "x_factor_ltriang_pos_def_cols differs at ({i},{k})"
                        );
                    }
                }
                newton_system.x_wc_valid = true;
                return;
            }
            _ => {}
        }
    }
    factor_ltriang_pos_def(&mut newton_system.w, m);
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:480-524
// ```c
// static void SolveReNewtonSystem(double x[], double z[],
//                                 const ReNewtonSystem * newton_system, double workspace[])
// ```
fn solve_re_newton_system(
    x: &mut [f64],
    z: &mut [f64],
    newton_system: &ReNewtonSystem,
    workspace: &mut [f64],
) {
    let alphsize = newton_system.alphsize;
    let n = alphsize * alphsize;
    let ma = 2 * alphsize - 1;
    let m = if newton_system.constrain_rel_entropy {
        ma + 1
    } else {
        ma
    };

    for i in 0..n {
        workspace[i] = x[i] * newton_system.dinv[i];
    }
    multiply_by_a(1.0, z, alphsize, -1.0, workspace);
    if newton_system.constrain_rel_entropy {
        for i in 0..n {
            z[m - 1] -= newton_system.grad_re[i] * workspace[i];
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:587-587
    // ```c
    //     Nlm_SolveLtriangPosDef(z, m, W);
    // ```
    // Dispatch point of LOSAT_X_COMPFAST / LOSAT_X_COMPSHADOW: the reference path is the
    // `solve_ltriang_pos_def` call in the `else` branch.
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:171-177
    // ```c
    //     for (i = 0;  i < n;  i++) {
    //         temp = x[i];
    //         for (j = 0;  j < i;  j++) {
    //             temp -= L[i][j] * x[j];
    //         }
    //         x[i] = temp/L[i][i];
    //     }
    // ```
    if newton_system.x_wc_valid {
        if x_comp_mode() == 2 {
            let mut check = z[..m].to_vec();
            x_solve_ltriang_pos_def_cols(&mut check, m, &newton_system.w, &newton_system.x_wc);
            solve_ltriang_pos_def(z, m, &newton_system.w);
            assert!(
                check
                    .iter()
                    .zip(z.iter())
                    .all(|(a, b)| a.to_bits() == b.to_bits()),
                "x_solve_ltriang_pos_def_cols differs from the reference"
            );
        } else {
            x_solve_ltriang_pos_def_cols(z, m, &newton_system.w, &newton_system.x_wc);
        }
    } else {
        solve_ltriang_pos_def(z, m, &newton_system.w);
    }
    if newton_system.constrain_rel_entropy {
        for i in 0..n {
            x[i] += newton_system.grad_re[i] * z[m - 1];
        }
    }
    multiply_by_at(1.0, x, alphsize, 1.0, z);
    for i in 0..n {
        x[i] *= newton_system.dinv[i];
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:538-565
// ```c
// static void EvaluateReFunctions(...)
// ```
fn evaluate_re_functions(
    values: &mut [f64; 2],
    grads: &mut [Vec<f64>],
    alphsize: usize,
    x: &[f64],
    q: &[f64],
    scores: &[f64],
    constrain_rel_entropy: bool,
) {
    values[0] = 0.0;
    values[1] = 0.0;
    for k in 0..alphsize * alphsize {
        let mut temp = (x[k] / q[k]).ln();
        values[0] += x[k] * temp;
        grads[0][k] = temp + 1.0;
        if constrain_rel_entropy {
            temp += scores[k];
            values[1] += x[k] * temp;
            grads[1][k] = temp + 1.0;
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:578-598
// ```c
// static void ComputeScoresFromProbs(double scores[], int alphsize, ...)
// ```
fn compute_scores_from_probs(
    scores: &mut [f64],
    alphsize: usize,
    target_freqs: &[f64],
    row_freqs: &[f64; COMPO_NUM_TRUE_AA],
    col_freqs: &[f64; COMPO_NUM_TRUE_AA],
) {
    for i in 0..alphsize {
        for j in 0..alphsize {
            let k = i * alphsize + j;
            scores[k] = (target_freqs[k] / (row_freqs[i] * col_freqs[j])).ln();
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:602-713
// ```c
// int Blast_OptimizeTargetFrequencies(...)
// ```
fn optimize_target_frequencies(
    x: &mut [f64],
    alphsize: usize,
    q: &[f64],
    row_sums: &[f64; COMPO_NUM_TRUE_AA],
    col_sums: &[f64; COMPO_NUM_TRUE_AA],
    constrain_rel_entropy: bool,
    relative_entropy: f64,
) -> i32 {
    // EXPERIMENT (LOSAT_X_NEWTONEXACT): the same iteration, bit for bit, in
    // fixed-size storage (see x_newton_exact.rs).
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:686-687,693-696
    // ```c
    // int
    // Blast_OptimizeTargetFrequencies(double x[],
    // ...
    //                                 int constrain_rel_entropy,
    //                                 double relative_entropy,
    //                                 double tol,
    //                                 int maxits)
    // ```
    // Dispatch point of LOSAT_X_NEWTONEXACT / LOSAT_X_NEWTONEXACTSHADOW: the reference path is
    // `optimize_target_frequencies_reference`, the port of this function. The exact module handles only
    // `alphsize == 20` with `constrain_rel_entropy`.
    if constrain_rel_entropy && alphsize == COMPO_NUM_TRUE_AA {
        // EXPERIMENT (LOSAT_X_NEWTONLANES / LOSAT_X_NEWTONLANESSHADOW): a result computed ahead
        // (four problems per vector, x_newton_lanes.rs) for exactly this input.
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:686-687,693-696
        // ```c
        // int
        // Blast_OptimizeTargetFrequencies(double x[],
        // ...
        //                                 int constrain_rel_entropy,
        //                                 double relative_entropy,
        //                                 double tol,
        //                                 int maxits)
        // ```
        // No NCBI counterpart: the stored result was computed by the operations of this function on an
        // input with the same bits (q, row and column sums, relative entropy); the shadow mode runs the
        // port of this function as well and compares status and all 400 values by bits. It does not
        // change any value NCBI computes.
        let x_lanes_mode = super::x_newton_lanes::mode();
        if x_lanes_mode != 0 {
            if let Some(status) =
                super::x_newton_lanes::take(x, q, row_sums, col_sums, relative_entropy)
            {
                if x_lanes_mode == 2 {
                    let mut x_reference = vec![0.0f64; x.len()];
                    let status_reference = optimize_target_frequencies_reference(
                        &mut x_reference,
                        alphsize,
                        q,
                        row_sums,
                        col_sums,
                        constrain_rel_entropy,
                        relative_entropy,
                    );
                    assert_eq!(
                        status, status_reference,
                        "LOSAT_X_NEWTONLANESSHADOW: convergence status differs"
                    );
                    for (k, (a, b)) in x.iter().zip(x_reference.iter()).enumerate() {
                        assert!(
                            a.to_bits() == b.to_bits(),
                            "LOSAT_X_NEWTONLANESSHADOW: x[{k}] differs: {a:e} vs {b:e}"
                        );
                    }
                    super::x_newton_lanes::SHADOWED
                        .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                }
                return status;
            }
        }
        match super::x_newton_exact::mode() {
            1 => {
                // No NCBI counterpart: call counter and timer of the xstats report; it does not change any value NCBI computes.
                let x_t0 = crate::utils::xstats::now();
                crate::utils::xstats::add(&crate::utils::xstats::COMPO_ADJ_CALLS, 1);
                let status = super::x_newton_exact::optimize_target_frequencies_exact(
                    x,
                    q,
                    row_sums,
                    col_sums,
                    relative_entropy,
                );
                // No NCBI counterpart: timer of the xstats report; it does not change any value NCBI computes.
                crate::utils::xstats::add_ns(&crate::utils::xstats::COMPO_NS, x_t0);
                return status;
            }
            2 => {
                // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:686-687,693-696
                // ```c
                // int
                // Blast_OptimizeTargetFrequencies(double x[],
                // ...
                //                                 int constrain_rel_entropy,
                //                                 double relative_entropy,
                //                                 double tol,
                //                                 int maxits)
                // ```
                // Shadow mode runs the exact module and the reference port on the same input and compares status
                // and every target frequency by bits.
                X_NEWTONEXACT_SHADOWED.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                let mut x_fast = vec![0.0f64; x.len()];
                let status_fast = super::x_newton_exact::optimize_target_frequencies_exact(
                    &mut x_fast,
                    q,
                    row_sums,
                    col_sums,
                    relative_entropy,
                );
                let status = optimize_target_frequencies_reference(
                    x,
                    alphsize,
                    q,
                    row_sums,
                    col_sums,
                    constrain_rel_entropy,
                    relative_entropy,
                );
                assert_eq!(
                    status_fast, status,
                    "LOSAT_X_NEWTONEXACTSHADOW: convergence status differs"
                );
                for (k, (a, b)) in x_fast.iter().zip(x.iter()).enumerate() {
                    assert!(
                        a.to_bits() == b.to_bits(),
                        "LOSAT_X_NEWTONEXACTSHADOW: x[{k}] differs: {a:e} vs {b:e}"
                    );
                }
                return status;
            }
            _ => {}
        }
    }
    optimize_target_frequencies_reference(
        x,
        alphsize,
        q,
        row_sums,
        col_sums,
        constrain_rel_entropy,
        relative_entropy,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:686-687,693-696
// ```c
// int
// Blast_OptimizeTargetFrequencies(double x[],
// ...
//                                 int constrain_rel_entropy,
//                                 double relative_entropy,
//                                 double tol,
//                                 int maxits)
// ```
// The port below with `alphsize == 20` and `constrain_rel_entropy`, for the bitwise tests of
// `x_newton_lanes`.
#[cfg(test)]
pub(super) fn x_optimize_reference_for_test(
    x: &mut [f64],
    q: &[f64],
    row_sums: &[f64; COMPO_NUM_TRUE_AA],
    col_sums: &[f64; COMPO_NUM_TRUE_AA],
    relative_entropy: f64,
) -> i32 {
    optimize_target_frequencies_reference(
        x,
        COMPO_NUM_TRUE_AA,
        q,
        row_sums,
        col_sums,
        true,
        relative_entropy,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:686-687,693-697,745-748
// ```c
// int
// Blast_OptimizeTargetFrequencies(double x[],
// ...
//                                 int constrain_rel_entropy,
//                                 double relative_entropy,
//                                 double tol,
//                                 int maxits)
// {
// ...
//     ComputeScoresFromProbs(old_scores, alphsize, q, row_sums, col_sums);
//     /* Use q as the initial value for x */
//     memcpy(x, q, n * sizeof(double));
// ```
// The original `optimize_target_frequencies`, renamed; it is the port of this function.
fn optimize_target_frequencies_reference(
    x: &mut [f64],
    alphsize: usize,
    q: &[f64],
    row_sums: &[f64; COMPO_NUM_TRUE_AA],
    col_sums: &[f64; COMPO_NUM_TRUE_AA],
    constrain_rel_entropy: bool,
    relative_entropy: f64,
) -> i32 {
    let n = alphsize * alphsize;
    let ma = 2 * alphsize - 1;
    let m = if constrain_rel_entropy { ma + 1 } else { ma };
    // No NCBI counterpart: call counter and timer of the xstats report; it does not change any value NCBI computes.
    let x_t0 = crate::utils::xstats::now();
    crate::utils::xstats::add(&crate::utils::xstats::COMPO_ADJ_CALLS, 1);

    let mut newton_system = re_newton_system_new(alphsize);
    let mut resids_x = vec![0.0; n];
    let mut resids_z = vec![0.0; ma + 1];
    let mut z = vec![0.0; ma + 1];
    let mut old_scores = vec![0.0; n];
    let mut workspace = vec![0.0; n];
    let mut grads = dense_matrix_new(2, n);
    let mut values = [0.0; 2];

    compute_scores_from_probs(&mut old_scores, alphsize, q, row_sums, col_sums);
    x.copy_from_slice(q);
    let mut rnorm = 0.0;
    let mut its = 0usize;
    while its <= K_COMPO_ADJUST_ITERATION_LIMIT {
        evaluate_re_functions(
            &mut values,
            &mut grads,
            alphsize,
            x,
            q,
            &old_scores,
            constrain_rel_entropy,
        );
        rnorm = calculate_residuals(
            &mut resids_x,
            alphsize,
            &mut resids_z,
            &values,
            &grads,
            row_sums,
            col_sums,
            x,
            &z,
            constrain_rel_entropy,
            relative_entropy,
        );
        if !(rnorm > K_COMPO_ADJUST_ERR_TOLERANCE) {
            break;
        }
        its += 1;
        if its <= K_COMPO_ADJUST_ITERATION_LIMIT {
            // No NCBI counterpart: Newton iteration counter of the xstats report; it does not change any value NCBI computes.
            crate::utils::xstats::add(&crate::utils::xstats::NEWTON_ITERS, 1);
            factor_re_newton_system(
                &mut newton_system,
                x,
                &z,
                &grads,
                constrain_rel_entropy,
                &mut workspace,
            );
            solve_re_newton_system(
                &mut resids_x,
                &mut resids_z[..m],
                &newton_system,
                &mut workspace,
            );
            let mut alpha = step_bound(x, &resids_x, 1.0 / 0.95);
            alpha *= 0.95;
            add_vectors(x, alpha, &resids_x);
            add_vectors(&mut z[..m], alpha, &resids_z[..m]);
        }
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:786-792
    // ```c
    // converged = 0;
    // if (its <= maxits && rnorm <= tol) {
    //     if (!constrain_rel_entropy || z[m - 1] < 1) {
    //         converged = 1;
    //     }
    // }
    // ```
    // No NCBI counterpart: timer of the xstats report; it does not change any value NCBI computes.
    crate::utils::xstats::add_ns(&crate::utils::xstats::COMPO_NS, x_t0);
    let converged = its <= K_COMPO_ADJUST_ITERATION_LIMIT
        && rnorm <= K_COMPO_ADJUST_ERR_TOLERANCE
        && (!constrain_rel_entropy || z[m - 1] < 1.0);
    if converged {
        0
    } else {
        1
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1336-1408
// ```c
// int Blast_CompositionMatrixAdj(int ** matrix, int alphsize, ...)
// ```
fn composition_matrix_adjust(
    matrix_info: &BlastMatrixInfo,
    matrix_adjust_rule: EMatrixAdjustRule,
    query_length: i32,
    subject_length: i32,
    query_prob_std: &[f64; COMPO_LARGEST_ALPHABET],
    subject_prob_std: &[f64; COMPO_LARGEST_ALPHABET],
    pseudocounts: i32,
    specified_re: f64,
    workspace: &mut BlastCompositionWorkspace,
) -> Result<Option<AdjustedProteinMatrix>> {
    let mut row_probs = [0.0; COMPO_NUM_TRUE_AA];
    let mut col_probs = [0.0; COMPO_NUM_TRUE_AA];
    gather_letter_probs(&mut row_probs, query_prob_std);
    gather_letter_probs(&mut col_probs, subject_prob_std);

    let mut desired_re = 0.0;
    match matrix_adjust_rule {
        EMatrixAdjustRule::UnconstrainedRelEntropy => {}
        EMatrixAdjustRule::RelEntropyOldMatrixNewContext => {
            if let Some(value) =
                entropy_old_freq_new_context(&workspace.mat_b, &row_probs, &col_probs)?
            {
                desired_re = value;
            }
        }
        EMatrixAdjustRule::RelEntropyOldMatrixOldContext => {
            desired_re = target_freq_entropy(&workspace.mat_b);
        }
        EMatrixAdjustRule::UserSpecifiedRelEntropy => {
            desired_re = specified_re;
        }
        _ => bail!(
            "unsupported matrix-adjust rule {:?} for composition matrix adjustment",
            matrix_adjust_rule
        ),
    }

    apply_pseudocounts(
        &mut row_probs,
        query_length,
        &workspace.first_standard_freq,
        pseudocounts,
    );
    apply_pseudocounts(
        &mut col_probs,
        subject_length,
        &workspace.second_standard_freq,
        pseudocounts,
    );

    let mut x = vec![0.0; COMPO_NUM_TRUE_AA * COMPO_NUM_TRUE_AA];
    let mut q = vec![0.0; COMPO_NUM_TRUE_AA * COMPO_NUM_TRUE_AA];
    for i in 0..COMPO_NUM_TRUE_AA {
        for j in 0..COMPO_NUM_TRUE_AA {
            q[i * COMPO_NUM_TRUE_AA + j] = workspace.mat_b[i][j];
        }
    }
    let status = optimize_target_frequencies(
        &mut x,
        COMPO_NUM_TRUE_AA,
        &q,
        &row_probs,
        &col_probs,
        desired_re > 0.0,
        desired_re,
    );
    if status != 0 {
        return Ok(None);
    }
    for i in 0..COMPO_NUM_TRUE_AA {
        for j in 0..COMPO_NUM_TRUE_AA {
            workspace.mat_final[i][j] = x[i * COMPO_NUM_TRUE_AA + j];
        }
    }
    Ok(Some(scores_std_alphabet(
        &workspace.mat_final,
        matrix_info,
        &row_probs,
        &col_probs,
    )))
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/unified_pvalues.c:208-260
// ```c
// double Blast_CompositionPvalue(double lambda) { ... }
// double Blast_Overall_P_Value(double p_comp, double p_alignment) { ... }
// ```
pub fn blast_composition_pvalue(lambda: f64) -> f64 {
    let initial_u_score = (lambda - COMPO_MIN_LAMBDA) / LAMBDA_BIN_SIZE;
    if initial_u_score < LOW_LAMBDA_BIN_CUT as f64 {
        return P_LAMBDA_TABLE[LOW_LAMBDA_BIN_CUT];
    }
    if initial_u_score > (LAMBDA_BIN_TOTAL - 1) as f64 {
        return 1.0;
    }
    let bin = initial_u_score as usize;
    if bin == LAMBDA_BIN_TOTAL - 1 {
        return P_LAMBDA_TABLE[LAMBDA_BIN_TOTAL - 1];
    }
    let delta = initial_u_score - bin as f64;
    (1.0 - delta) * P_LAMBDA_TABLE[bin] + delta * P_LAMBDA_TABLE[bin + 1]
}

pub fn blast_overall_p_value(p_comp: f64, p_alignment: f64) -> f64 {
    let product = p_comp * p_alignment;
    if product > 0.0 {
        product * (1.0 - product.ln())
    } else {
        product
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1414-1530
// ```c
// int Blast_AdjustScores(...);
// ```
pub(crate) fn blast_adjust_scores(
    matrix_info: &BlastMatrixInfo,
    query_composition: &BlastAminoAcidComposition,
    query_length: i32,
    subject_composition: &BlastAminoAcidComposition,
    subject_length: i32,
    composition_adjust_mode: BlastCompoAdjustMode,
    composition_test_index: i32,
    workspace: &mut BlastCompositionWorkspace,
    calc_lambda: fn(&[f64], i32, i32, f64) -> Result<f64>,
) -> Result<Option<BlastAdjustScoresResult>> {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1433-1439
    // ```c
    // if (query_composition->numTrueAminoAcids == 0 ||
    //     subject_composition->numTrueAminoAcids == 0) {
    //     return 1;
    // }
    // ```
    if query_composition.num_true_amino_acids == 0 || subject_composition.num_true_amino_acids == 0
    {
        return Ok(None);
    }

    let mut permuted_query_probs = [0.0; COMPO_NUM_TRUE_AA];
    let mut permuted_match_probs = [0.0; COMPO_NUM_TRUE_AA];
    gather_letter_probs(&mut permuted_query_probs, &query_composition.prob);
    gather_letter_probs(&mut permuted_match_probs, &subject_composition.prob);

    let pvalue_for_this_pair = if composition_test_index > 0 {
        Some(compute_pair_composition_pvalue(
            &permuted_query_probs,
            &permuted_match_probs,
        )?)
    } else {
        None
    };

    let mut matrix_adjust_rule = if matrix_info.position_based
        || composition_adjust_mode == BlastCompoAdjustMode::CompositionBasedStats
    {
        EMatrixAdjustRule::CompoScaleOldMatrix
    } else {
        blast_choose_matrix_adjust_rule(
            query_length,
            subject_length,
            &permuted_query_probs,
            &permuted_match_probs,
            matrix_info.matrix,
            composition_adjust_mode,
        )?
    };
    if matrix_adjust_rule != EMatrixAdjustRule::CompoScaleOldMatrix {
        // EXPERIMENT (LOSAT_X_ADJMEMO_STATS): how often the same inputs recur
        // (a memo of the adjusted matrix would be exact: the function is pure).
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1501-1507
        // ```c
        //     if (eCompoScaleOldMatrix != *matrix_adjust_rule) {
        //         /* Try matrix optimization, if it fails to converge, we
        //            fall back to traditional scaling below */
        //         int status =
        //             Blast_CompositionMatrixAdj(matrix,
        //                                        alphsize,
        //                                        *matrix_adjust_rule,
        // ```
        // No NCBI counterpart: LOSAT_X_ADJMEMO_STATS counts how many times the call that follows (the port of
        // `Blast_CompositionMatrixAdj`) is made and with how many distinct inputs; it does not change any value NCBI computes.
        x_adjmemo_record(query_composition, subject_composition, matrix_adjust_rule);
        if let Some(adjusted_matrix) = composition_matrix_adjust(
            matrix_info,
            matrix_adjust_rule,
            query_composition.num_true_amino_acids,
            subject_composition.num_true_amino_acids,
            &query_composition.prob,
            &subject_composition.prob,
            20,
            K_FIXED_RE_BLOSUM62,
            workspace,
        )? {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1509-1517
            // ```c
            // *ratioToPassBack = 1.0;
            // if (status <= 0)
            //     return status;
            // ```
            return Ok(Some(BlastAdjustScoresResult {
                matrix_adjust_rule,
                adjusted_matrix,
                pvalue_for_this_pair,
                lambda_ratio: 1.0,
            }));
        }
        matrix_adjust_rule = EMatrixAdjustRule::CompoScaleOldMatrix;
    }

    let lambda_ratio = compute_lambda_ratio(
        matrix_info,
        &query_composition.prob,
        &subject_composition.prob,
        composition_test_index > 0,
        calc_lambda,
    );
    let adjusted_matrix = build_scaled_square_matrix(
        matrix_info,
        &query_composition.prob,
        &subject_composition.prob,
        lambda_ratio,
    );

    Ok(Some(BlastAdjustScoresResult {
        matrix_adjust_rule,
        adjusted_matrix,
        pvalue_for_this_pair,
        lambda_ratio,
    }))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1485-1505
// ```c
//     if (matrixInfo->positionBased ||
//         composition_adjust_mode == eCompositionBasedStats) {
//         /* Use old-style composition-based statistics unconditionally. */
//         *matrix_adjust_rule =  eCompoScaleOldMatrix;
//     } else {
//         /* else call Yi-Kuo's code to choose mode for matrix adjustment. */
//         *matrix_adjust_rule =
//             Blast_ChooseMatrixAdjustRule(queryLength, subjectLength,
//                                          permutedQueryProbs,
//                                          permutedMatchProbs,
//                                          matrixInfo->matrixName,
//                                          composition_adjust_mode);
//     }  /* end else call Yi-Kuo's code to choose mode for matrix adjustment. */
//
//     if (eCompoScaleOldMatrix != *matrix_adjust_rule) {
//         /* Try matrix optimization, if it fails to converge, we
//            fall back to traditional scaling below */
//         int status =
//             Blast_CompositionMatrixAdj(matrix,
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1354-1355,1385-1399
// ```c
//     s_GatherLetterProbs(row_probs, stdaa_row_probs, alphsize);
//     s_GatherLetterProbs(col_probs, stdaa_col_probs, alphsize);
// ...
//     Blast_ApplyPseudocounts(row_probs, length1,
//                             NRrecord->first_standard_freq, pseudocounts);
//     Blast_ApplyPseudocounts(col_probs, length2,
//                             NRrecord->second_standard_freq, pseudocounts);
//
//     status =
//         Blast_OptimizeTargetFrequencies(&NRrecord->mat_final[0][0],
//                                         COMPO_NUM_TRUE_AA,
//                                         &iteration_count,
//                                         &NRrecord->mat_b[0][0],
//                                         row_probs, col_probs,
//                                         (desired_re > 0.0),
//                                         desired_re,
// ```
// EXPERIMENT (LOSAT_X_NEWTONLANES): the prefix of `blast_adjust_scores` and `composition_matrix_adjust`
// up to the `optimize_target_frequencies` call, with the same helpers, arguments and order, returning
// that call's input. `None` where the call would not reach the relative-entropy-constrained Newton
// iteration (a composition without true amino acids, position-based or `CompositionBasedStats`,
// rule `CompoScaleOldMatrix` or an unsupported rule, `desired_re <= 0`, an error). The p-value of the
// pair (`compositionTestIndex`) does not enter the input and is not computed.
// No NCBI counterpart: the input is used only to compute a result ahead of the real call, which uses
// it only when its own input has the same bits; it does not change any value NCBI computes.
/// The Newton input that `blast_adjust_scores` with these arguments would pass to
/// `optimize_target_frequencies`, if it reaches the constrained iteration.
pub(crate) fn x_newton_input(
    matrix_info: &BlastMatrixInfo,
    query_composition: &BlastAminoAcidComposition,
    query_length: i32,
    subject_composition: &BlastAminoAcidComposition,
    subject_length: i32,
    composition_adjust_mode: BlastCompoAdjustMode,
    workspace: &BlastCompositionWorkspace,
) -> Option<super::x_newton_lanes::NewtonInput> {
    if query_composition.num_true_amino_acids == 0 || subject_composition.num_true_amino_acids == 0
    {
        return None;
    }
    let mut permuted_query_probs = [0.0; COMPO_NUM_TRUE_AA];
    let mut permuted_match_probs = [0.0; COMPO_NUM_TRUE_AA];
    gather_letter_probs(&mut permuted_query_probs, &query_composition.prob);
    gather_letter_probs(&mut permuted_match_probs, &subject_composition.prob);
    if matrix_info.position_based
        || composition_adjust_mode == BlastCompoAdjustMode::CompositionBasedStats
    {
        return None;
    }
    let matrix_adjust_rule = blast_choose_matrix_adjust_rule(
        query_length,
        subject_length,
        &permuted_query_probs,
        &permuted_match_probs,
        matrix_info.matrix,
        composition_adjust_mode,
    )
    .ok()?;
    if matrix_adjust_rule == EMatrixAdjustRule::CompoScaleOldMatrix {
        return None;
    }
    // composition_matrix_adjust(matrix_info, rule, query num_true_amino_acids, subject
    // num_true_amino_acids, query prob, subject prob, 20, K_FIXED_RE_BLOSUM62, workspace)
    let mut row_probs = [0.0; COMPO_NUM_TRUE_AA];
    let mut col_probs = [0.0; COMPO_NUM_TRUE_AA];
    gather_letter_probs(&mut row_probs, &query_composition.prob);
    gather_letter_probs(&mut col_probs, &subject_composition.prob);
    let mut desired_re = 0.0;
    match matrix_adjust_rule {
        EMatrixAdjustRule::UnconstrainedRelEntropy => {}
        EMatrixAdjustRule::RelEntropyOldMatrixNewContext => {
            if let Some(value) =
                entropy_old_freq_new_context(&workspace.mat_b, &row_probs, &col_probs).ok()?
            {
                desired_re = value;
            }
        }
        EMatrixAdjustRule::RelEntropyOldMatrixOldContext => {
            desired_re = target_freq_entropy(&workspace.mat_b);
        }
        EMatrixAdjustRule::UserSpecifiedRelEntropy => {
            desired_re = K_FIXED_RE_BLOSUM62;
        }
        _ => return None,
    }
    if !(desired_re > 0.0) {
        return None;
    }
    apply_pseudocounts(
        &mut row_probs,
        query_composition.num_true_amino_acids,
        &workspace.first_standard_freq,
        20,
    );
    apply_pseudocounts(
        &mut col_probs,
        subject_composition.num_true_amino_acids,
        &workspace.second_standard_freq,
        20,
    );
    let mut q = [0.0; COMPO_NUM_TRUE_AA * COMPO_NUM_TRUE_AA];
    for i in 0..COMPO_NUM_TRUE_AA {
        for j in 0..COMPO_NUM_TRUE_AA {
            q[i * COMPO_NUM_TRUE_AA + j] = workspace.mat_b[i][j];
        }
    }
    Some(super::x_newton_lanes::NewtonInput {
        q,
        row: row_probs,
        col: col_probs,
        relative_entropy: desired_re,
    })
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1444-1481
// ```c
// if (compositionTestIndex > 0) {
//     ...
//     Blast_CalcLambdaFullPrecision(&lambdaForPair, &iter_count, scores, ...);
//     if (iter_count >= kLambdaIterationLimit) {
//         lambdaForPair = COMPO_MIN_LAMBDA;
//     }
//     *pvalueForThisPair = Blast_CompositionPvalue(lambdaForPair);
// }
// ```
fn compute_pair_composition_pvalue(
    permuted_query_probs: &[f64; COMPO_NUM_TRUE_AA],
    permuted_match_probs: &[f64; COMPO_NUM_TRUE_AA],
) -> Result<f64> {
    let mut scores = dense_matrix_new(COMPO_NUM_TRUE_AA, COMPO_NUM_TRUE_AA);
    for (i, row) in scores.iter_mut().enumerate().take(COMPO_NUM_TRUE_AA) {
        for (j, value) in row.iter_mut().enumerate().take(COMPO_NUM_TRUE_AA) {
            *value = protein_score_direct(ScoringMatrix::Blosum62, i, j) as f64;
        }
    }
    let (mut lambda_for_pair, iter_count) = calc_lambda_full_precision_from_scores(
        &scores,
        COMPO_NUM_TRUE_AA,
        permuted_query_probs,
        permuted_match_probs,
    );
    if iter_count >= K_LAMBDA_ITERATION_LIMIT {
        lambda_for_pair = COMPO_MIN_LAMBDA;
    }
    Ok(blast_composition_pvalue(lambda_for_pair))
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1105-1155
// ```c
// int
// Blast_CompositionBasedStats(...,
//                             double (*calc_lambda)(double*,int,int,double),
//                             int pValueAdjustment)
// {
//     ...
//     correctUngappedLambda =
//         calc_lambda(scoreArray, obs_min, obs_max, ss->ungappedLambda);
//     *LambdaRatio = correctUngappedLambda / ss->ungappedLambda;
//     ...
// }
// ```
fn compute_lambda_ratio(
    matrix_info: &BlastMatrixInfo,
    query_prob: &[f64; COMPO_LARGEST_ALPHABET],
    subject_prob: &[f64; COMPO_LARGEST_ALPHABET],
    pvalue_adjustment: bool,
    calc_lambda: fn(&[f64], i32, i32, f64) -> Result<f64>,
) -> f64 {
    let (score_probs, obs_min, obs_max) =
        get_matrix_score_probs(&matrix_info.start_matrix, query_prob, subject_prob);
    let ratio = match calc_lambda(&score_probs, obs_min, obs_max, matrix_info.ungapped_lambda) {
        Ok(lambda) => lambda / matrix_info.ungapped_lambda,
        Err(_) => -1.0,
    };

    let mut lambda_ratio = ratio;
    if !pvalue_adjustment {
        lambda_ratio = lambda_ratio.min(1.0);
    }
    lambda_ratio.max(LAMBDA_RATIO_LOWER_BOUND)
}

fn build_scaled_square_matrix(
    matrix_info: &BlastMatrixInfo,
    row_prob: &[f64; COMPO_LARGEST_ALPHABET],
    col_prob: &[f64; COMPO_LARGEST_ALPHABET],
    lambda_ratio: f64,
) -> AdjustedProteinMatrix {
    let scaled_lambda = matrix_info.ungapped_lambda / lambda_ratio;
    let mut scores = matrix_info.start_freq_ratios;
    freq_ratio_to_score_square(
        &mut scores,
        COMPO_LARGEST_ALPHABET,
        COMPO_LARGEST_ALPHABET,
        scaled_lambda,
    );
    set_xuo_scores(&mut scores, row_prob, col_prob);
    let mut rounded = round_score_matrix(&scores);
    for index in 0..COMPO_LARGEST_ALPHABET {
        rounded[index][25] = matrix_info.start_matrix[index][25];
        rounded[25][index] = matrix_info.start_matrix[25][index];
    }
    AdjustedProteinMatrix { scores: rounded }
}

fn get_matrix_score_probs(
    matrix: &[[i32; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET],
    query_prob: &[f64; COMPO_LARGEST_ALPHABET],
    subject_prob: &[f64; COMPO_LARGEST_ALPHABET],
) -> (Vec<f64>, i32, i32) {
    let mut obs_min = 0;
    let mut obs_max = 0;
    for row in 0..COMPO_LARGEST_ALPHABET {
        for &col in &TRUE_CHAR_POSITIONS {
            let score = matrix[row][col];
            if score > obs_max {
                obs_max = score;
            }
            if score < obs_min && score > COMPO_SCORE_MIN {
                obs_min = score;
            }
        }
    }

    let mut score_prob = vec![0.0; (obs_max - obs_min + 1) as usize];
    for row in 0..COMPO_LARGEST_ALPHABET {
        for &col in &TRUE_CHAR_POSITIONS {
            let score = matrix[row][col];
            if score >= obs_min {
                score_prob[(score - obs_min) as usize] += query_prob[row] * subject_prob[col];
            }
        }
    }
    (score_prob, obs_min, obs_max)
}

fn calc_freq_ratios(
    ratios: &mut [[f64; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET],
    row_prob: &[f64; COMPO_LARGEST_ALPHABET],
    col_prob: &[f64; COMPO_LARGEST_ALPHABET],
) {
    for row in 0..COMPO_LARGEST_ALPHABET {
        for col in 0..COMPO_LARGEST_ALPHABET {
            let product = row_prob[row] * col_prob[col];
            ratios[row][col] = if product > 0.0 {
                ratios[row][col] / product
            } else {
                0.0
            };
        }
    }
}

fn high_pair_frequencies(letter_probs: &[f64; COMPO_NUM_TRUE_AA], length: i32) -> bool {
    if length <= LENGTH_LOWER_THRESHOLD {
        return false;
    }
    let mut max = 0.0;
    let mut second = 0.0;
    for &prob in letter_probs {
        if prob > second {
            second = prob;
            if prob > max {
                second = max;
                max = prob;
            }
        }
    }
    (max + second) > HIGH_PAIR_THRESHOLD
}

fn high_pair_either_seq(
    query_prob: &[f64; COMPO_NUM_TRUE_AA],
    query_length: i32,
    match_prob: &[f64; COMPO_NUM_TRUE_AA],
    match_length: i32,
) -> bool {
    high_pair_frequencies(query_prob, query_length)
        || high_pair_frequencies(match_prob, match_length)
}

fn blast_choose_matrix_adjust_rule(
    length1: i32,
    length2: i32,
    prob_array1: &[f64; COMPO_NUM_TRUE_AA],
    prob_array2: &[f64; COMPO_NUM_TRUE_AA],
    matrix: ScoringMatrix,
    composition_adjust_mode: BlastCompoAdjustMode,
) -> Result<EMatrixAdjustRule> {
    match composition_adjust_mode {
        BlastCompoAdjustMode::NoCompositionBasedStats => Ok(EMatrixAdjustRule::DontAdjustMatrix),
        BlastCompoAdjustMode::CompositionBasedStats => Ok(EMatrixAdjustRule::CompoScaleOldMatrix),
        BlastCompoAdjustMode::ForceFullMatrixAdjust => {
            Ok(EMatrixAdjustRule::UserSpecifiedRelEntropy)
        }
        BlastCompoAdjustMode::CompositionMatrixAdjust => {
            if matrix != ScoringMatrix::Blosum62 {
                bail!("matrix-adjust rule selection is only ported for BLOSUM62");
            }

            let mut corr_factor = 0.0;
            for i in 0..COMPO_NUM_TRUE_AA {
                corr_factor +=
                    (prob_array1[i] - BLOSUM62_BG[i]) * (prob_array2[i] - BLOSUM62_BG[i]);
            }
            let d_m_mat = blast_get_relative_entropy(prob_array2, &BLOSUM62_BG);
            let d_q_mat = blast_get_relative_entropy(prob_array1, &BLOSUM62_BG);
            let d_m_q = blast_get_relative_entropy(prob_array2, prob_array1);

            let angle = if d_m_mat > 0.0 && d_q_mat > 0.0 {
                let cos_angle = ((d_m_mat * d_m_mat) + (d_q_mat * d_q_mat) - (d_m_q * d_m_q))
                    / (2.0 * d_m_mat * d_q_mat);
                cos_angle.clamp(-1.0, 1.0).acos().to_degrees()
            } else {
                0.0
            };

            let len_large = (length1.max(length2)) as f64;
            let len_small = (length1.min(length2)) as f64;
            let _ = corr_factor;
            if high_pair_either_seq(prob_array1, length1, prob_array2, length2) {
                Ok(EMatrixAdjustRule::UserSpecifiedRelEntropy)
            } else if d_m_q > QUERY_MATCH_DISTANCE_THRESHOLD
                && len_small > 0.0
                && len_large / len_small > LENGTH_RATIO_THRESHOLD
                && angle > ANGLE_DEGREE_THRESHOLD
            {
                Ok(EMatrixAdjustRule::CompoScaleOldMatrix)
            } else {
                Ok(EMatrixAdjustRule::UserSpecifiedRelEntropy)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::blast_stat::composition::compute_lambda_from_score_probs;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:551-572
    // ```c
    // s_CalcLambda(double probs[], int min_score, int max_score, double lambda0)
    // {
    //     return Blast_KarlinLambdaNR(&freq, lambda0);
    // }
    // ```
    fn test_compute_lambda_from_score_probs(
        score_probs: &[f64],
        min_score: i32,
        max_score: i32,
        lambda0: f64,
    ) -> Result<f64> {
        compute_lambda_from_score_probs(score_probs, min_score, max_score, lambda0)
            .map_err(anyhow::Error::msg)
    }

    // NCBI composition_adjustment/composition_adjustment.c:1204-1248:
    // if (i + kCompositionMargin < left) left = i + kCompositionMargin;
    // if (i - kCompositionMargin > right) right = i - kCompositionMargin;
    #[test]
    fn translated_composition_range_honors_stop_margin_strict_boundaries() {
        let mut data = vec![1u8; 100];
        assert_eq!(
            blast_get_composition_range(&data, 50, 60).unwrap(),
            (0, 100)
        );
        data[20] = 25;
        data[80] = 25;
        assert_eq!(
            blast_get_composition_range(&data, 50, 60).unwrap(),
            (41, 60)
        );
        assert_eq!(
            blast_get_composition_range(&data, 50, 50).unwrap(),
            (41, 60)
        );
        data[29] = 25;
        assert_eq!(
            blast_get_composition_range(&data, 50, 50).unwrap(),
            (50, 60)
        );
    }

    // NCBI composition_adjustment/redo_alignment.c:1234-1254;
    // composition_adjustment/composition_adjustment.c:1414-1530:
    // Blast_AdjustScores(matrix, query_composition, query.length,
    //                    &subject_composition, subject.length,
    //                    scaledMatrixInfo, compo_adjust_mode,
    //                    RE_pseudocounts, NRrecord, &matrix_adjust_rule, ...);
    #[test]
    fn pinned_tblastn_mode2_adjustment_rules_match_ncbi_call_inputs() {
        for case in [
            "seg_hard_query_20260924_default",
            "multi_query_20260924_default",
        ] {
            let trace = std::fs::read_to_string(format!(
                "{}/../docs/evidence/tlosan_stage_d/kappa_composition_matrix_scores_20260924/{case}.tsv",
                env!("CARGO_MANIFEST_DIR")
            ))
            .unwrap();
            let rows: Vec<_> = trace.lines().collect();
            assert_eq!(rows.len() % 6, 0);
            // NCBI core/blast_kappa.c:2229-2236:
            // self->ungappedLambda = sbp->kbp_ideal->Lambda / scale_factor;
            // Blast_Int4MatrixFromFreq(self->startMatrix, self->cols,
            //                          self->startFreqRatios, self->ungappedLambda);
            let matrix_row: Vec<_> = rows[1].split('\t').collect();
            assert_eq!(matrix_row[0], "K_MATRIX");
            let matrix_lambda: f64 = matrix_row[5].parse().unwrap();
            let matrix_info = build_matrix_info(ScoringMatrix::Blosum62, matrix_lambda).unwrap();
            let mut workspace = BlastCompositionWorkspace::new_blosum62();
            for block in rows.chunks_exact(6) {
                let call: Vec<_> = block[0].split('\t').collect();
                let matrix: Vec<_> = block[1].split('\t').collect();
                let result: Vec<_> = block[4].split('\t').collect();
                assert_eq!(matrix[0], "K_MATRIX");
                assert_eq!(matrix[1], call[1]);
                assert_eq!(matrix[2], "28");
                assert_eq!(matrix[3], "28");
                assert_eq!(matrix[4], "0");
                assert_eq!(
                    matrix[5].parse::<f64>().unwrap().to_bits(),
                    matrix_lambda.to_bits()
                );
                assert_eq!(matrix_info.start_matrix[1][1], matrix[6].parse().unwrap());
                assert_eq!(matrix_info.start_matrix[1][16], matrix[7].parse().unwrap());
                assert_eq!(matrix_info.start_matrix[22][22], matrix[8].parse().unwrap());
                assert_eq!(call[0], "K_COMP_CALL");
                assert_eq!(call[4], "2");
                assert_eq!(call[5], "20");
                assert_eq!(result[0], "K_COMP_RESULT");
                assert_eq!(call[1], result[1]);
                assert_eq!(result[2], "0");
                let compositions: Vec<_> = block[2..4]
                    .iter()
                    .enumerate()
                    .map(|(idx, row)| {
                        let fields: Vec<_> = row.split('\t').collect();
                        assert_eq!(fields[0], "K_COMP");
                        assert_eq!(fields[1], call[1]);
                        assert_eq!(fields[2], if idx == 0 { "query" } else { "subject" });
                        let mut comp = BlastAminoAcidComposition::empty();
                        comp.num_true_amino_acids = fields[3].parse().unwrap();
                        assert_eq!(fields.len(), 32);
                        for (value, hex) in comp.prob.iter_mut().zip(&fields[4..]) {
                            *value = f64::from_bits(u64::from_str_radix(hex, 16).unwrap());
                        }
                        comp
                    })
                    .collect();
                let adjusted = blast_adjust_scores(
                    &matrix_info,
                    &compositions[0],
                    call[2].parse().unwrap(),
                    &compositions[1],
                    call[3].parse().unwrap(),
                    BlastCompoAdjustMode::CompositionMatrixAdjust,
                    call[6].parse().unwrap(),
                    &mut workspace,
                    test_compute_lambda_from_score_probs,
                )
                .unwrap()
                .unwrap();
                assert_eq!(
                    adjusted.matrix_adjust_rule as i32,
                    result[3].parse().unwrap(),
                    "{case}: call {}",
                    call[1]
                );
                assert_eq!(
                    adjusted.lambda_ratio.to_bits(),
                    result[4].parse::<f64>().unwrap().to_bits(),
                    "{case}: call {}",
                    call[1]
                );
                let adjusted_row: Vec<_> = block[5].split('\t').collect();
                assert_eq!(adjusted_row[0], "K_ADJUSTED");
                assert_eq!(adjusted_row[1], call[1]);
                assert_eq!(adjusted_row.len(), 2 + 28 * 28);
                for (index, (observed, expected)) in adjusted
                    .adjusted_matrix
                    .scores
                    .iter()
                    .flatten()
                    .zip(&adjusted_row[2..])
                    .enumerate()
                {
                    assert_eq!(
                        *observed,
                        expected.parse::<i32>().unwrap(),
                        "{case}: call {}, matrix index {index}",
                        call[1]
                    );
                }
            }
        }
    }

    #[test]
    fn test_read_aa_composition_ignores_x_and_merges_u_into_c() {
        let seq = [1u8, 3, 21, 24];
        let composition = read_aa_composition(&seq);
        assert_eq!(composition.num_true_amino_acids, 3);
        assert!((composition.prob[1] - (1.0 / 3.0)).abs() < 1.0e-12);
        assert!((composition.prob[3] - (2.0 / 3.0)).abs() < 1.0e-12);
        assert_eq!(composition.prob[21], 0.0);
        assert_eq!(composition.prob[24], 0.0);
    }

    #[test]
    fn test_build_matrix_info_blosum62_rebuilds_start_matrix_from_freq_ratios() {
        let info = build_matrix_info(ScoringMatrix::Blosum62, 0.3176).unwrap();
        assert!((info.start_freq_ratios[1][1] - 3.90294070).abs() < 1.0e-8);
        assert!((info.start_freq_ratios[22][22] - 9.83220341).abs() < 1.0e-8);
        assert_eq!(info.start_matrix[1][1], 4);
        assert_eq!(info.start_matrix[1][16], -2);
        assert!(info.start_freq_ratios[1][1] > 0.0);
    }

    #[test]
    fn test_blast_composition_pvalue_clamps_left_tail() {
        assert_eq!(
            blast_composition_pvalue(COMPO_MIN_LAMBDA),
            P_LAMBDA_TABLE[LOW_LAMBDA_BIN_CUT]
        );
    }

    #[test]
    fn test_blast_overall_p_value_matches_fisher_style_formula() {
        let combined = blast_overall_p_value(0.2, 0.3);
        let product = 0.06f64;
        assert!((combined - (product * (1.0 - product.ln()))).abs() < 1.0e-12);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:686-687,693-696
    // ```c
    // int
    // Blast_OptimizeTargetFrequencies(double x[],
    // ...
    //                                 int constrain_rel_entropy,
    //                                 double relative_entropy,
    //                                 double tol,
    //                                 int maxits)
    // ```
    // The ignored test times the port of this function against `optimize_target_frequencies_exact`
    // and checks the last result bit for bit.
    /// Timing of the Newton kernels on realistic inputs (run with --ignored --nocapture;
    /// LOSAT_X_COMPFAST selects the reference variant, LOSAT_X_NEWTONEXACT_LEVEL the exact one).
    #[test]
    #[ignore]
    fn x_bench_newton_kernels() {
        let mut state = 0x2545_F491_4F6C_DD1Du64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        let workspace = BlastCompositionWorkspace::new_blosum62();
        let mut q = vec![0.0f64; COMPO_NUM_TRUE_AA * COMPO_NUM_TRUE_AA];
        for i in 0..COMPO_NUM_TRUE_AA {
            for j in 0..COMPO_NUM_TRUE_AA {
                q[i * COMPO_NUM_TRUE_AA + j] = workspace.mat_b[i][j];
            }
        }
        let cases: Vec<([f64; COMPO_NUM_TRUE_AA], [f64; COMPO_NUM_TRUE_AA])> = (0..256)
            .map(|_| {
                let mut row = [0.0f64; COMPO_NUM_TRUE_AA];
                let mut col = [0.0f64; COMPO_NUM_TRUE_AA];
                let (mut sr, mut sc) = (0.0, 0.0);
                for k in 0..COMPO_NUM_TRUE_AA {
                    row[k] = BLOSUM62_BG[k] * (0.3 + 1.4 * next());
                    col[k] = BLOSUM62_BG[k] * (0.3 + 1.4 * next());
                    sr += row[k];
                    sc += col[k];
                }
                for k in 0..COMPO_NUM_TRUE_AA {
                    row[k] /= sr;
                    col[k] /= sc;
                }
                apply_pseudocounts(&mut row, 300, &BLOSUM62_BG, 20);
                apply_pseudocounts(&mut col, 250, &BLOSUM62_BG, 20);
                (row, col)
            })
            .collect();
        let mut x_ref = vec![0.0f64; q.len()];
        let mut x_new = vec![0.0f64; q.len()];
        let reps = 8;
        let t0 = std::time::Instant::now();
        let mut status_sum = 0;
        for _ in 0..reps {
            for (row, col) in &cases {
                status_sum += optimize_target_frequencies_reference(
                    &mut x_ref,
                    COMPO_NUM_TRUE_AA,
                    &q,
                    row,
                    col,
                    true,
                    K_FIXED_RE_BLOSUM62,
                );
            }
        }
        let t_ref = t0.elapsed().as_secs_f64() / (reps * cases.len()) as f64 * 1e6;
        let t0 = std::time::Instant::now();
        let mut status_sum2 = 0;
        for _ in 0..reps {
            for (row, col) in &cases {
                status_sum2 += super::super::x_newton_exact::optimize_target_frequencies_exact(
                    &mut x_new,
                    &q,
                    row,
                    col,
                    K_FIXED_RE_BLOSUM62,
                );
            }
        }
        let t_new = t0.elapsed().as_secs_f64() / (reps * cases.len()) as f64 * 1e6;
        // bitwise check on the last case
        for (a, b) in x_ref.iter().zip(x_new.iter()) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
        eprintln!(
            "[BENCH] reference {:.1} us/call (status sum {}), exact {:.1} us/call (status sum {})",
            t_ref, status_sum, t_new, status_sum2
        );
        super::super::x_newton_exact::prof_print();
    }

    #[test]
    fn test_blast_adjust_scores_returns_adjusted_matrix_for_blosum62() {
        let matrix_info = build_matrix_info(ScoringMatrix::Blosum62, 0.3176).unwrap();
        let query = read_aa_composition(&[1u8, 1, 1, 4, 4, 4, 17, 17]);
        let subject = read_aa_composition(&[1u8, 4, 17, 1, 4, 17, 1, 4]);
        let mut workspace = BlastCompositionWorkspace::new_blosum62();
        let adjusted = blast_adjust_scores(
            &matrix_info,
            &query,
            8,
            &subject,
            8,
            BlastCompoAdjustMode::CompositionMatrixAdjust,
            0,
            &mut workspace,
            test_compute_lambda_from_score_probs,
        )
        .unwrap()
        .expect("adjusted matrix result");
        assert_eq!(
            adjusted.matrix_adjust_rule,
            EMatrixAdjustRule::UserSpecifiedRelEntropy
        );
        assert_eq!(adjusted.lambda_ratio, 1.0);
        assert_ne!(
            adjusted.adjusted_matrix.scores[1][1],
            matrix_info.start_matrix[1][1]
        );
    }

    #[test]
    fn test_blast_adjust_scores_returns_none_for_zero_composition_ncbi_style() {
        let matrix_info = build_matrix_info(ScoringMatrix::Blosum62, 0.3176).unwrap();
        let zero = read_aa_composition(&[21u8, 21, 21]);
        let subject = read_aa_composition(&[1u8, 4, 17]);
        let mut workspace = BlastCompositionWorkspace::new_blosum62();
        let adjusted = blast_adjust_scores(
            &matrix_info,
            &zero,
            3,
            &subject,
            3,
            BlastCompoAdjustMode::CompositionMatrixAdjust,
            0,
            &mut workspace,
            test_compute_lambda_from_score_probs,
        )
        .unwrap();
        assert!(adjusted.is_none());
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/composition_adjustment.c:682-690
    // ```c
    // M[eBchar][eXchar] = s_CalcXScore(M[eBchar], alphsize, 1, col_probs);
    // M[eXchar][eBchar] = s_CalcXScore(&M[0][eBchar], alphsize, cols, row_probs);
    // M[eZchar][eXchar] = s_CalcXScore(M[eZchar], alphsize, 1, col_probs);
    // ```
    #[test]
    fn test_set_xuo_scores_uses_row_average_for_pairwise_ambiguity_to_x() {
        let mut scores = [[0.0; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET];
        let mut row_probs = [0.0; COMPO_LARGEST_ALPHABET];
        let mut col_probs = [0.0; COMPO_LARGEST_ALPHABET];

        scores[E_BCHAR][E_DCHAR] = -10.0;
        scores[E_BCHAR][E_NCHAR] = -30.0;
        scores[E_DCHAR][E_BCHAR] = -7.0;
        scores[E_NCHAR][E_BCHAR] = -11.0;

        col_probs[E_DCHAR] = 0.75;
        col_probs[E_NCHAR] = 0.25;
        row_probs[E_DCHAR] = 0.25;
        row_probs[E_NCHAR] = 0.75;

        set_xuo_scores(&mut scores, &row_probs, &col_probs);

        assert!((scores[E_BCHAR][E_XCHAR] - (-15.0)).abs() < 1.0e-12);
        assert!((scores[E_XCHAR][E_BCHAR] - (-10.0)).abs() < 1.0e-12);
    }
}
