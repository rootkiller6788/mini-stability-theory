/* fts_analysis.c ? extension module */
#include "fts_analysis.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/* fts_analysis_op0: absolute value */
int fts_analysis_op0(int n, const double *x, double *out) {
    if (!x || !out || n < 1 || n > 4096) return -1;
    for (int j = 0; j < n; j++) {
        double v = x[j];
        if (!isfinite(v)) { out[j] = 0.0; continue; }
        out[j] = fabs;
        if (!isfinite(out[j])) out[j] = 0.0;
    }
    return 0;
}

/* fts_analysis_batch0: batched absolute value */
int fts_analysis_batch0(int n, int m, const double *x, double *out) {
    if (!x || !out || n < 1 || m < 1 || n*m > 65536) return -1;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double v = x[i*n + j];
            if (!isfinite(v)) { out[i*n + j] = 0.0; continue; }
            out[i*n + j] = fabs;
            if (!isfinite(out[i*n + j])) out[i*n + j] = 0.0;
        }
    }
    return 0;
}

/* fts_analysis_op1: square */
int fts_analysis_op1(int n, const double *x, double *out) {
    if (!x || !out || n < 1 || n > 4096) return -1;
    for (int j = 0; j < n; j++) {
        double v = x[j];
        if (!isfinite(v)) { out[j] = 0.0; continue; }
        out[j] = (v)*(v);
        if (!isfinite(out[j])) out[j] = 0.0;
    }
    return 0;
}

/* fts_analysis_batch1: batched square */
int fts_analysis_batch1(int n, int m, const double *x, double *out) {
    if (!x || !out || n < 1 || m < 1 || n*m > 65536) return -1;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double v = x[i*n + j];
            if (!isfinite(v)) { out[i*n + j] = 0.0; continue; }
            out[i*n + j] = (v)*(v);
            if (!isfinite(out[i*n + j])) out[i*n + j] = 0.0;
        }
    }
    return 0;
}

/* fts_analysis_op2: cube */
int fts_analysis_op2(int n, const double *x, double *out) {
    if (!x || !out || n < 1 || n > 4096) return -1;
    for (int j = 0; j < n; j++) {
        double v = x[j];
        if (!isfinite(v)) { out[j] = 0.0; continue; }
        out[j] = (v)*(v)*(v);
        if (!isfinite(out[j])) out[j] = 0.0;
    }
    return 0;
}

/* fts_analysis_batch2: batched cube */
int fts_analysis_batch2(int n, int m, const double *x, double *out) {
    if (!x || !out || n < 1 || m < 1 || n*m > 65536) return -1;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double v = x[i*n + j];
            if (!isfinite(v)) { out[i*n + j] = 0.0; continue; }
            out[i*n + j] = (v)*(v)*(v);
            if (!isfinite(out[i*n + j])) out[i*n + j] = 0.0;
        }
    }
    return 0;
}

/* fts_analysis_op3: logistic sigmoid */
int fts_analysis_op3(int n, const double *x, double *out) {
    if (!x || !out || n < 1 || n > 4096) return -1;
    for (int j = 0; j < n; j++) {
        double v = x[j];
        if (!isfinite(v)) { out[j] = 0.0; continue; }
        out[j] = 1.0/(1.0+exp(-(v)));
        if (!isfinite(out[j])) out[j] = 0.0;
    }
    return 0;
}

/* fts_analysis_batch3: batched logistic sigmoid */
int fts_analysis_batch3(int n, int m, const double *x, double *out) {
    if (!x || !out || n < 1 || m < 1 || n*m > 65536) return -1;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double v = x[i*n + j];
            if (!isfinite(v)) { out[i*n + j] = 0.0; continue; }
            out[i*n + j] = 1.0/(1.0+exp(-(v)));
            if (!isfinite(out[i*n + j])) out[i*n + j] = 0.0;
        }
    }
    return 0;
}

/* fts_analysis_op4: ReLU activation */
int fts_analysis_op4(int n, const double *x, double *out) {
    if (!x || !out || n < 1 || n > 4096) return -1;
    for (int j = 0; j < n; j++) {
        double v = x[j];
        if (!isfinite(v)) { out[j] = 0.0; continue; }
        out[j] = (v)>0?(v):0;
        if (!isfinite(out[j])) out[j] = 0.0;
    }
    return 0;
}

/* fts_analysis_batch4: batched ReLU activation */
int fts_analysis_batch4(int n, int m, const double *x, double *out) {
    if (!x || !out || n < 1 || m < 1 || n*m > 65536) return -1;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double v = x[i*n + j];
            if (!isfinite(v)) { out[i*n + j] = 0.0; continue; }
            out[i*n + j] = (v)>0?(v):0;
            if (!isfinite(out[i*n + j])) out[i*n + j] = 0.0;
        }
    }
    return 0;
}

/* fts_analysis_op5: hyperbolic tangent */
int fts_analysis_op5(int n, const double *x, double *out) {
    if (!x || !out || n < 1 || n > 4096) return -1;
    for (int j = 0; j < n; j++) {
        double v = x[j];
        if (!isfinite(v)) { out[j] = 0.0; continue; }
        out[j] = tanh(v);
        if (!isfinite(out[j])) out[j] = 0.0;
    }
    return 0;
}

/* fts_analysis_batch5: batched hyperbolic tangent */
int fts_analysis_batch5(int n, int m, const double *x, double *out) {
    if (!x || !out || n < 1 || m < 1 || n*m > 65536) return -1;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double v = x[i*n + j];
            if (!isfinite(v)) { out[i*n + j] = 0.0; continue; }
            out[i*n + j] = tanh(v);
            if (!isfinite(out[i*n + j])) out[i*n + j] = 0.0;
        }
    }
    return 0;
}

/* fts_analysis_op6: softplus */
int fts_analysis_op6(int n, const double *x, double *out) {
    if (!x || !out || n < 1 || n > 4096) return -1;
    for (int j = 0; j < n; j++) {
        double v = x[j];
        if (!isfinite(v)) { out[j] = 0.0; continue; }
        out[j] = log1p(exp(v));
        if (!isfinite(out[j])) out[j] = 0.0;
    }
    return 0;
}

/* fts_analysis_batch6: batched softplus */
int fts_analysis_batch6(int n, int m, const double *x, double *out) {
    if (!x || !out || n < 1 || m < 1 || n*m > 65536) return -1;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double v = x[i*n + j];
            if (!isfinite(v)) { out[i*n + j] = 0.0; continue; }
            out[i*n + j] = log1p(exp(v));
            if (!isfinite(out[i*n + j])) out[i*n + j] = 0.0;
        }
    }
    return 0;
}

/* fts_analysis_op7: Gaussian RBF */
int fts_analysis_op7(int n, const double *x, double *out) {
    if (!x || !out || n < 1 || n > 4096) return -1;
    for (int j = 0; j < n; j++) {
        double v = x[j];
        if (!isfinite(v)) { out[j] = 0.0; continue; }
        out[j] = exp(-(v)*(v));
        if (!isfinite(out[j])) out[j] = 0.0;
    }
    return 0;
}

/* fts_analysis_batch7: batched Gaussian RBF */
int fts_analysis_batch7(int n, int m, const double *x, double *out) {
    if (!x || !out || n < 1 || m < 1 || n*m > 65536) return -1;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double v = x[i*n + j];
            if (!isfinite(v)) { out[i*n + j] = 0.0; continue; }
            out[i*n + j] = exp(-(v)*(v));
            if (!isfinite(out[i*n + j])) out[i*n + j] = 0.0;
        }
    }
    return 0;
}

/* fts_analysis_op8: swish activation */
int fts_analysis_op8(int n, const double *x, double *out) {
    if (!x || !out || n < 1 || n > 4096) return -1;
    for (int j = 0; j < n; j++) {
        double v = x[j];
        if (!isfinite(v)) { out[j] = 0.0; continue; }
        out[j] = (v)/(1.0+exp(-(v)));
        if (!isfinite(out[j])) out[j] = 0.0;
    }
    return 0;
}

/* fts_analysis_batch8: batched swish activation */
int fts_analysis_batch8(int n, int m, const double *x, double *out) {
    if (!x || !out || n < 1 || m < 1 || n*m > 65536) return -1;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double v = x[i*n + j];
            if (!isfinite(v)) { out[i*n + j] = 0.0; continue; }
            out[i*n + j] = (v)/(1.0+exp(-(v)));
            if (!isfinite(out[i*n + j])) out[i*n + j] = 0.0;
        }
    }
    return 0;
}

double fts_analysis_mean(int n, const double *x) {
    if (!x || n < 1) return 0.0;
    double sum=0; for(int j=0;j<n;j++) sum+=x[j]; return sum/(double)n;
}

double fts_analysis_variance(int n, const double *x) {
    if (!x || n < 1) return 0.0;
    double m=0,s=0; for(int j=0;j<n;j++) m+=x[j]; m/=(double)n; for(int j=0;j<n;j++){double d=x[j]-m; s+=d*d;} return s/(double)(n-1);
}

int fts_analysis_minmax(int n, const double *x, double *out) {
    if (!x || n < 1) return -1;
    double mn=x[0],mx=x[0]; for(int j=1;j<n;j++){if(x[j]<mn)mn=x[j];if(x[j]>mx)mx=x[j];} out[0]=mn; out[1]=mx; return 0;
}

double fts_analysis_sum(int n, const double *x) {
    if (!x || n < 1) return 0.0;
    double s=0; for(int j=0;j<n;j++) s+=x[j]; return s;
}

double fts_analysis_norm_l2(int n, const double *x) {
    if (!x || n < 1) return 0.0;
    double s=0; for(int j=0;j<n;j++) s+=x[j]*x[j]; return sqrt(s);
}

double fts_analysis_norm_l1(int n, const double *x) {
    if (!x || n < 1) return 0.0;
    double s=0; for(int j=0;j<n;j++) s+=fabs(x[j]); return s;
}

double fts_analysis_norm_linf(int n, const double *x) {
    if (!x || n < 1) return 0.0;
    double m=fabs(x[0]); for(int j=1;j<n;j++){if(fabs(x[j])>m)m=fabs(x[j]);} return m;
}
