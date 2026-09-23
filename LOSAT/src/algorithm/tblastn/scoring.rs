//! Pinned NCBI protein matrix gap pairs used for TBLASTN CLI validation.

// NCBI reference: c++/src/algo/blast/core/blast_stat.c:2949-2994,3578-3638
// ```c
// matrix_info = MatrixInfoNew("BLOSUM62", blosum62_values, blosum62_prefs, BLOSUM62_VALUES_MAX);
// if (strcasecmp(matrix_info->name, matrix_name) == 0)
// {
//     values = matrix_info->values;
// }
// if (BLAST_Nint(values[index][0]) == gap_open &&
//     BLAST_Nint(values[index][1]) == gap_extend)
// ```
pub(super) struct MatrixParams {
    pub preferred: (i32, i32),
    pub allowed: &'static [(i32, i32)],
}

// NCBI reference: c++/src/algo/blast/core/blast_stat.c:183-198
// ```c
// static array_of_8 blosum45_values[BLOSUM45_VALUES_MAX] = {
//     {(double) INT2_MAX, (double) INT2_MAX, (double) INT2_MAX, 0.2291, 0.0924, 0.2514, 0.9113, -5.7, 0.641318, 9.611060, 9.611060},
//     {13, 3, (double) INT2_MAX, 0.207, 0.049, 0.14, 1.5, -22, 0.671128, 35.855900, 35.963900},
// ```
const BLOSUM45: &[(i32, i32)] = &[
    (13, 3),
    (12, 3),
    (11, 3),
    (10, 3),
    (16, 2),
    (15, 2),
    (14, 2),
    (13, 2),
    (12, 2),
    (19, 1),
    (18, 1),
    (17, 1),
    (16, 1),
];

// NCBI reference: c++/src/algo/blast/core/blast_stat.c:219-236
// ```c
// static array_of_8 blosum50_values[BLOSUM50_VALUES_MAX] = {
//     {(double) INT2_MAX, (double) INT2_MAX, (double) INT2_MAX, 0.2318, 0.112, 0.3362, 0.6895, -4.0, 0.609639, 5.388310, 5.388310},
//     {13, 3, (double) INT2_MAX, 0.212, 0.063, 0.19, 1.1, -16, 0.639287, 18.113800, 18.202800},
// ```
const BLOSUM50: &[(i32, i32)] = &[
    (13, 3),
    (12, 3),
    (11, 3),
    (10, 3),
    (9, 3),
    (16, 2),
    (15, 2),
    (14, 2),
    (13, 2),
    (12, 2),
    (19, 1),
    (18, 1),
    (17, 1),
    (16, 1),
    (15, 1),
];

// NCBI reference: c++/src/algo/blast/core/blast_stat.c:258-271
// ```c
// static array_of_8 blosum62_values[BLOSUM62_VALUES_MAX] = {
//     {(double) INT2_MAX, (double) INT2_MAX, (double) INT2_MAX, 0.3176, 0.134, 0.4012, 0.7916, -3.2, 0.623757, 4.964660, 4.964660},
//     {11, 2, (double) INT2_MAX, 0.297, 0.082, 0.27, 1.1, -10, 0.641766, 12.673800, 12.757600},
// ```
const BLOSUM62: &[(i32, i32)] = &[
    (11, 2),
    (10, 2),
    (9, 2),
    (8, 2),
    (7, 2),
    (6, 2),
    (13, 1),
    (12, 1),
    (11, 1),
    (10, 1),
    (9, 1),
];

// NCBI reference: c++/src/algo/blast/core/blast_stat.c:290-301
// ```c
// static array_of_8 blosum80_values[BLOSUM80_VALUES_MAX] = {
//     {(double) INT2_MAX, (double) INT2_MAX, (double) INT2_MAX, 0.3430, 0.177, 0.6568, 0.5222, -1.6, 0.564057, 1.918130, 1.918130},
//     {25, 2, (double) INT2_MAX, 0.342, 0.17, 0.66, 0.52, -1.6, 0.563956, 1.731000, 1.731300},
// ```
const BLOSUM80: &[(i32, i32)] = &[
    (25, 2),
    (13, 2),
    (9, 2),
    (8, 2),
    (7, 2),
    (6, 2),
    (11, 1),
    (10, 1),
    (9, 1),
];

// NCBI reference: c++/src/algo/blast/core/blast_stat.c:317-326
// ```c
// static array_of_8 blosum90_values[BLOSUM90_VALUES_MAX] = {
//     {(double) INT2_MAX, (double) INT2_MAX, (double) INT2_MAX, 0.3346, 0.190, 0.7547, 0.4434, -1.4 , 0.544178, 1.377760, 1.377760},
//     {9, 2, (double) INT2_MAX, 0.310, 0.12, 0.46, 0.67, -6 , 0.570267, 4.232290, 4.334170},
// ```
const BLOSUM90: &[(i32, i32)] = &[(9, 2), (8, 2), (7, 2), (6, 2), (11, 1), (10, 1), (9, 1)];

// NCBI reference: c++/src/algo/blast/core/blast_stat.c:340-357
// ```c
// static array_of_8 pam250_values[PAM250_VALUES_MAX] = {
//     {(double) INT2_MAX, (double) INT2_MAX, (double) INT2_MAX, 0.2252, 0.0868, 0.2223, 0.98, -5.0, 0.660059, 11.754300, 11.754300},
//     {15, 3, (double) INT2_MAX, 0.205, 0.049, 0.13, 1.6, -23, 0.687656, 34.578400, 34.928000},
// ```
const PAM250: &[(i32, i32)] = &[
    (15, 3),
    (14, 3),
    (13, 3),
    (12, 3),
    (11, 3),
    (17, 2),
    (16, 2),
    (15, 2),
    (14, 2),
    (13, 2),
    (21, 1),
    (20, 1),
    (19, 1),
    (18, 1),
    (17, 1),
];

// NCBI reference: c++/src/algo/blast/core/blast_stat.c:379-391
// ```c
// static array_of_8 pam30_values[PAM30_VALUES_MAX] = {
//     {(double) INT2_MAX, (double) INT2_MAX, (double) INT2_MAX, 0.3400, 0.283, 1.754, 0.1938, -0.3, 0.436164, 0.161818, 0.161818},
//     {7, 2, (double) INT2_MAX, 0.305, 0.15, 0.87, 0.35, -3, 0.479087, 1.014010, 1.162730},
// ```
const PAM30: &[(i32, i32)] = &[
    (7, 2),
    (6, 2),
    (5, 2),
    (10, 1),
    (9, 1),
    (8, 1),
    (15, 3),
    (14, 2),
    (14, 1),
    (13, 3),
];

// NCBI reference: c++/src/algo/blast/core/blast_stat.c:409-419
// ```c
// static array_of_8 pam70_values[PAM70_VALUES_MAX] = {
//     {(double) INT2_MAX, (double) INT2_MAX, (double) INT2_MAX, 0.3345, 0.229, 1.029, 0.3250,   -0.7, 0.511296, 0.633439, 0.633439},
//     {8, 2, (double) INT2_MAX, 0.301, 0.12, 0.54, 0.56, -5, 0.549019, 2.881650, 3.025710},
// ```
const PAM70: &[(i32, i32)] = &[
    (8, 2),
    (7, 2),
    (6, 2),
    (11, 1),
    (10, 1),
    (9, 1),
    (11, 2),
    (12, 3),
];

// NCBI reference: c++/src/algo/blast/core/blast_stat.c:578-581
// ```c
// static array_of_8 prot_idenity_values[PROT_IDENTITY_VALUES_MAX] = {
//     {(double) INT2_MAX, (double) INT2_MAX, (double) INT2_MAX, 0.28768, 0.282, 1.69, 0.1703, -0.3, 0.43828, 0.16804, 0.16804},
//     {15, 2, (double) INT2_MAX, 0.2835, 0.255, 1.49, 0.19, -1, 0.44502, 0.24613, 0.22743}
// ```
const IDENTITY: &[(i32, i32)] = &[(15, 2)];

// NCBI reference: c++/src/algo/blast/core/blast_stat.c:2952-2994,3374-3399
// ```c
// matrix_info = MatrixInfoNew("BLOSUM62", blosum62_values, blosum62_prefs, BLOSUM62_VALUES_MAX);
// if(pref_flags[i]==BLAST_MATRIX_BEST) {
//     (*gap_existence) = gapOpen_arr[i];
//     (*gap_extension) = gapExtend_arr[i];
//     break;
// }
// ```
pub(super) fn matrix_params(matrix: &str) -> Option<MatrixParams> {
    let (preferred, allowed) = match matrix.to_ascii_uppercase().as_str() {
        "BLOSUM45" => ((14, 2), BLOSUM45),
        "BLOSUM50" => ((13, 2), BLOSUM50),
        "BLOSUM62" => ((11, 1), BLOSUM62),
        "BLOSUM80" => ((10, 1), BLOSUM80),
        "BLOSUM90" => ((10, 1), BLOSUM90),
        "PAM250" => ((15, 2), PAM250),
        "PAM30" => ((9, 1), PAM30),
        "PAM70" => ((10, 1), PAM70),
        "IDENTITY" => ((15, 2), IDENTITY),
        _ => return None,
    };
    Some(MatrixParams { preferred, allowed })
}
// NCBI reference: c++/src/algo/blast/core/blast_options.c:1174-1208
// ```c
// if(strcasecmp(matrixName, "BLOSUM45") == 0) *threshold = 14;
// else if(strcasecmp(matrixName, "PAM30") == 0) *threshold = 16;
// if (Blast_SubjectIsTranslated(program_number) == TRUE) *threshold += 2;
// ```
pub(super) fn suggested_threshold(matrix: &str) -> f64 {
    let raw = match matrix.to_ascii_uppercase().as_str() {
        "BLOSUM45" | "PAM70" => 14.0,
        "BLOSUM80" => 12.0,
        "PAM30" => 16.0,
        "IDENTITY" => 27.0,
        _ => 11.0,
    };
    raw + 2.0
}

// NCBI reference: c++/src/algo/blast/core/blast_options.c:1211-1236
// ```c
// if(strcasecmp(matrixName, "BLOSUM45") == 0) *window_size = 60;
// else if(strcasecmp(matrixName, "BLOSUM80") == 0) *window_size = 25;
// else if(strcasecmp(matrixName, "PAM30") == 0) *window_size = 15;
// else if(strcasecmp(matrixName, "PAM70") == 0) *window_size = 20;
// else *window_size = kB62_windowsize;
// ```
pub(super) fn suggested_window_size(matrix: &str) -> usize {
    match matrix.to_ascii_uppercase().as_str() {
        "BLOSUM45" => 60,
        "BLOSUM80" => 25,
        "PAM30" => 15,
        "PAM70" => 20,
        _ => 40,
    }
}
