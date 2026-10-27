/* Comparison-only oracle driver (never linked into LOSAT): NCBI 2.17.0 (598d8ae6) BLAST_SpougeStoE, blast_stat.c:5175-5232, copied verbatim below, with
   boost_erf.c compiled from the NCBI tree. Driver reads the same lines as losat_spouge. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "boost_erf.h"
typedef int Int4; typedef long long Int8;
#define MAX(a,b) ((a)>(b)?(a):(b))
#define ASSERT(x)
typedef struct { double Lambda, K, logK, H; } Blast_KarlinBlk;
typedef struct { double Lambda, C, G, a, Alpha, Sigma, a_un, Alpha_un, b, Beta, Tau; Int8 db_length; } Blast_GumbelBlk;
double
BLAST_SpougeStoE(Int4 y_,
                 Blast_KarlinBlk* kbp,
                 Blast_GumbelBlk* gbp,
                 Int4 m_, Int4 n_)
{
    // the score and lambda may have been rescaled.  We will compute the scaling factor
    // and use it to scale a, alpha, and Sigma similarly.
    double scale_factor = kbp->Lambda / gbp->Lambda;

    // the pair-wise e-value must be scaled back to db-wise e-value
    double db_scale_factor = (gbp->db_length) ?
            (double)gbp->db_length/(double)n_ : 1.0;

    double lambda_    = kbp->Lambda;
    double k_         = kbp->K;
    double ai_hat_    = gbp->a * scale_factor;
    double bi_hat_    = gbp->b;
    double alphai_hat_= gbp->Alpha * scale_factor;
    double betai_hat_ = gbp->Beta;
    double sigma_hat_ = gbp->Sigma * scale_factor;
    double tau_hat_   = gbp->Tau;

    /* here we consider symmetric matrix only */
    double aj_hat_    = ai_hat_;
    double bj_hat_    = bi_hat_;
    double alphaj_hat_= alphai_hat_;
    double betaj_hat_ = betai_hat_;

    /* this is 1/sqrt(2.0*PI) */
    static double const_val = 0.39894228040143267793994605993438;

    double m_li_y, vi_y, sqrt_vi_y, m_F, P_m_F;
    double n_lj_y, vj_y, sqrt_vj_y, n_F, P_n_F;
    double c_y, p1, p2, area;
    double e_value;

    m_li_y = m_ - (ai_hat_*y_ + bi_hat_);
    vi_y = MAX(2.0*alphai_hat_/lambda_, alphai_hat_*y_+betai_hat_);
    sqrt_vi_y = sqrt(vi_y);
    m_F = m_li_y/sqrt_vi_y;
    P_m_F = ErfC(-m_F / sqrt(2.0)) / 2.0;
    p1 = m_li_y * P_m_F + sqrt_vi_y * const_val * exp(-0.5*m_F*m_F);

    n_lj_y = n_ - (aj_hat_*y_ + bj_hat_);
    vj_y = MAX(2.0*alphaj_hat_/lambda_, alphaj_hat_*y_+betaj_hat_);
    sqrt_vj_y = sqrt(vj_y);
    n_F = n_lj_y/sqrt_vj_y;
    P_n_F = ErfC(-n_F / sqrt(2.0)) / 2.0;
    p2 = n_lj_y * P_n_F + sqrt_vj_y * const_val * exp(-0.5*n_F*n_F);

    c_y = MAX(2.0*sigma_hat_/lambda_, sigma_hat_*y_+tau_hat_);
    area = p1 * p2 + c_y * P_m_F * P_n_F;

    e_value = area * k_ * exp(-lambda_ * y_) * db_scale_factor;
    ASSERT(e_value >= 0.0);

    return e_value;
}

double
SpougeDebug(Int4 y_,
                 Blast_KarlinBlk* kbp,
                 Blast_GumbelBlk* gbp,
                 Int4 m_, Int4 n_)
{
    // the score and lambda may have been rescaled.  We will compute the scaling factor
    // and use it to scale a, alpha, and Sigma similarly.
    double scale_factor = kbp->Lambda / gbp->Lambda;

    // the pair-wise e-value must be scaled back to db-wise e-value
    double db_scale_factor = (gbp->db_length) ?
            (double)gbp->db_length/(double)n_ : 1.0;

    double lambda_    = kbp->Lambda;
    double k_         = kbp->K;
    double ai_hat_    = gbp->a * scale_factor;
    double bi_hat_    = gbp->b;
    double alphai_hat_= gbp->Alpha * scale_factor;
    double betai_hat_ = gbp->Beta;
    double sigma_hat_ = gbp->Sigma * scale_factor;
    double tau_hat_   = gbp->Tau;

    /* here we consider symmetric matrix only */
    double aj_hat_    = ai_hat_;
    double bj_hat_    = bi_hat_;
    double alphaj_hat_= alphai_hat_;
    double betaj_hat_ = betai_hat_;

    /* this is 1/sqrt(2.0*PI) */
    static double const_val = 0.39894228040143267793994605993438;

    double m_li_y, vi_y, sqrt_vi_y, m_F, P_m_F;
    double n_lj_y, vj_y, sqrt_vj_y, n_F, P_n_F;
    double c_y, p1, p2, area;
    double e_value;

    m_li_y = m_ - (ai_hat_*y_ + bi_hat_);
    vi_y = MAX(2.0*alphai_hat_/lambda_, alphai_hat_*y_+betai_hat_);
    sqrt_vi_y = sqrt(vi_y);
    m_F = m_li_y/sqrt_vi_y;
    P_m_F = ErfC(-m_F / sqrt(2.0)) / 2.0;
    p1 = m_li_y * P_m_F + sqrt_vi_y * const_val * exp(-0.5*m_F*m_F);

    n_lj_y = n_ - (aj_hat_*y_ + bj_hat_);
    vj_y = MAX(2.0*alphaj_hat_/lambda_, alphaj_hat_*y_+betaj_hat_);
    sqrt_vj_y = sqrt(vj_y);
    n_F = n_lj_y/sqrt_vj_y;
    P_n_F = ErfC(-n_F / sqrt(2.0)) / 2.0;
    p2 = n_lj_y * P_n_F + sqrt_vj_y * const_val * exp(-0.5*n_F*n_F);

    c_y = MAX(2.0*sigma_hat_/lambda_, sigma_hat_*y_+tau_hat_);
    area = p1 * p2 + c_y * P_m_F * P_n_F;

    e_value = area * k_ * exp(-lambda_ * y_) * db_scale_factor;
{ double v=m_li_y; uint64_t b; memcpy(&b,&v,8); printf("  m_li_y %016llx %.17e\n",(unsigned long long)b,v); }
{ double v=vi_y; uint64_t b; memcpy(&b,&v,8); printf("  vi_y %016llx %.17e\n",(unsigned long long)b,v); }
{ double v=sqrt_vi_y; uint64_t b; memcpy(&b,&v,8); printf("  sqrt_vi_y %016llx %.17e\n",(unsigned long long)b,v); }
{ double v=m_F; uint64_t b; memcpy(&b,&v,8); printf("  m_F %016llx %.17e\n",(unsigned long long)b,v); }
{ double v=P_m_F; uint64_t b; memcpy(&b,&v,8); printf("  P_m_F %016llx %.17e\n",(unsigned long long)b,v); }
{ double v=p1; uint64_t b; memcpy(&b,&v,8); printf("  p1 %016llx %.17e\n",(unsigned long long)b,v); }
{ double v=n_lj_y; uint64_t b; memcpy(&b,&v,8); printf("  n_lj_y %016llx %.17e\n",(unsigned long long)b,v); }
{ double v=vj_y; uint64_t b; memcpy(&b,&v,8); printf("  vj_y %016llx %.17e\n",(unsigned long long)b,v); }
{ double v=sqrt_vj_y; uint64_t b; memcpy(&b,&v,8); printf("  sqrt_vj_y %016llx %.17e\n",(unsigned long long)b,v); }
{ double v=n_F; uint64_t b; memcpy(&b,&v,8); printf("  n_F %016llx %.17e\n",(unsigned long long)b,v); }
{ double v=P_n_F; uint64_t b; memcpy(&b,&v,8); printf("  P_n_F %016llx %.17e\n",(unsigned long long)b,v); }
{ double v=p2; uint64_t b; memcpy(&b,&v,8); printf("  p2 %016llx %.17e\n",(unsigned long long)b,v); }
{ double v=c_y; uint64_t b; memcpy(&b,&v,8); printf("  c_y %016llx %.17e\n",(unsigned long long)b,v); }
{ double v=area; uint64_t b; memcpy(&b,&v,8); printf("  area %016llx %.17e\n",(unsigned long long)b,v); }
{ double v=e_value; uint64_t b; memcpy(&b,&v,8); printf("  e_value %016llx %.17e\n",(unsigned long long)b,v); }
{ double v=exp(-lambda_ * y_), w=exp(-0.5*m_F*m_F), x=exp(-0.5*n_F*n_F), za=-m_F / sqrt(2.0), zb=-n_F / sqrt(2.0); uint64_t b1,b2,b3,b4,b5; memcpy(&b1,&v,8); memcpy(&b2,&w,8); memcpy(&b3,&x,8); memcpy(&b4,&za,8); memcpy(&b5,&zb,8); printf("  exp_ly %016llx\n  exp_m %016llx\n  exp_n %016llx\n  erfc_arg_m %016llx\n  erfc_arg_n %016llx\n",(unsigned long long)b1,(unsigned long long)b2,(unsigned long long)b3,(unsigned long long)b4,(unsigned long long)b5); }
    ASSERT(e_value >= 0.0);

    return e_value;
}


static void gumbel(Blast_GumbelBlk* g, Int8 d) {
    g->Lambda = 0.267; g->C = 0.669720; g->G = 11 + 1;
    g->a = 1.9; g->Alpha = 42.6028; g->Sigma = 43.6362; g->a_un = 0.7916; g->Alpha_un = 4.964660;
    g->b = 2.0 * g->G * (g->a_un - g->a); g->Beta = 2.0 * g->G * (g->Alpha_un - g->Alpha);
    g->Tau = 2.0 * g->G * (g->Alpha_un - g->Sigma); g->db_length = d;
}
int main(void) {
    Blast_KarlinBlk k; Blast_GumbelBlk g; char line[256];
    k.Lambda = 0.267; k.K = 0.041; k.H = 0.14; k.Lambda /= 32.0; k.logK = log(k.K);
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] == 'S') {
            int y, m, n; long long d; double e; uint64_t b;
            sscanf(line + 2, "%d %d %d %lld", &y, &m, &n, &d);
            gumbel(&g, d);
            e = BLAST_SpougeStoE(y, &k, &g, m, n); memcpy(&b, &e, 8);
            if (line[1] == 'D') SpougeDebug(y, &k, &g, m, n);
            printf("S %d %d %d %lld %016llx %.17e\n", y, m, n, d, (unsigned long long)b, e);
        } else if (line[0] == 'E') {
            uint64_t zb, rb; double z, r; char hex[32];
            sscanf(line + 1, "%31s", hex); zb = strtoull(hex, 0, 16); memcpy(&z, &zb, 8);
            r = ErfC(z); memcpy(&rb, &r, 8);
            printf("E %s %016llx %.17e\n", hex, (unsigned long long)rb, r);
        }
    }
    return 0;
}
