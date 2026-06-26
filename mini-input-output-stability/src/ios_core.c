#include "ios_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ============================================================
 * Signal Creation / Destruction
 * ============================================================ */
IOSSignal* ios_signal_create(int cap) {
    IOSSignal* s = (IOSSignal*)calloc(1, sizeof(IOSSignal));
    if (!s) return NULL;
    s->capacity = cap > 0 ? cap : IOS_DEFAULT_SIGNAL_CAP;
    s->data = (double*)calloc((size_t)s->capacity, sizeof(double));
    if (!s->data) { free(s); return NULL; }
    s->length = 0;
    return s;
}

void ios_signal_free(IOSSignal* s) {
    if (s) { free(s->data); free(s); }
}

void ios_signal_clear(IOSSignal* s) {
    if (s) s->length = 0;
}

void ios_signal_append(IOSSignal* s, double v) {
    if (!s) return;
    if (s->length >= s->capacity) {
        s->capacity *= 2;
        double* d = (double*)realloc(s->data, (size_t)s->capacity * sizeof(double));
        if (!d) return;
        s->data = d;
    }
    s->data[s->length++] = v;
}

/* ============================================================
 * Signal Statistics
 * ============================================================ */
double ios_signal_energy(const IOSSignal* s) {
    if (!s || s->length < 1) return 0.0;
    double e = 0.0;
    for (int i = 0; i < s->length; i++)
        e += s->data[i] * s->data[i];
    return e;
}

double ios_signal_power(const IOSSignal* s) {
    if (!s || s->length < 1) return 0.0;
    return ios_signal_energy(s) / (double)s->length;
}

double ios_signal_max(const IOSSignal* s) {
    if (!s || s->length < 1) return 0.0;
    double mx = s->data[0];
    for (int i = 1; i < s->length; i++)
        if (s->data[i] > mx) mx = s->data[i];
    return mx;
}

double ios_signal_min(const IOSSignal* s) {
    if (!s || s->length < 1) return 0.0;
    double mn = s->data[0];
    for (int i = 1; i < s->length; i++)
        if (s->data[i] < mn) mn = s->data[i];
    return mn;
}

double ios_signal_mean(const IOSSignal* s) {
    if (!s || s->length < 1) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < s->length; i++)
        sum += s->data[i];
    return sum / (double)s->length;
}

double ios_signal_rms(const IOSSignal* s) {
    if (!s || s->length < 1) return 0.0;
    return sqrt(ios_signal_power(s));
}

double ios_signal_variance(const IOSSignal* s) {
    if (!s || s->length < 2) return 0.0;
    double m = ios_signal_mean(s);
    double v = 0.0;
    for (int i = 0; i < s->length; i++) {
        double d = s->data[i] - m;
        v += d * d;
    }
    return v / (double)s->length;
}

double ios_signal_norm(const IOSSignal* s, IOSNormType type) {
    if (!s || s->length < 1) return 0.0;
    double n = 0.0;
    switch (type) {
        case IOS_NORM_L1:
            for (int i = 0; i < s->length; i++)
                n += fabs(s->data[i]);
            break;
        case IOS_NORM_L2: {
            double sum = 0.0;
            for (int i = 0; i < s->length; i++)
                sum += s->data[i] * s->data[i];
            n = sqrt(sum);
            break;
        }
        case IOS_NORM_LINF:
            for (int i = 0; i < s->length; i++) {
                double a = fabs(s->data[i]);
                if (a > n) n = a;
            }
            break;
    }
    return n;
}

double ios_signal_skewness(const IOSSignal* s) {
    if (!s || s->length < 3) return 0.0;
    double m = ios_signal_mean(s);
    double v = ios_signal_variance(s);
    if (v < IOS_EPS) return 0.0;
    double s3 = 0.0;
    for (int i = 0; i < s->length; i++) {
        double d = s->data[i] - m;
        s3 += d * d * d;
    }
    double std = sqrt(v);
    return s3 / (double)s->length / (std * std * std);
}

double ios_signal_kurtosis(const IOSSignal* s) {
    if (!s || s->length < 4) return 0.0;
    double m = ios_signal_mean(s);
    double s2 = 0.0, s4 = 0.0;
    for (int i = 0; i < s->length; i++) {
        double d = s->data[i] - m;
        s2 += d * d;
        s4 += d * d * d * d;
    }
    double v = s2 / (double)s->length;
    if (v < IOS_EPS) return 0.0;
    return (s4 / (double)s->length) / (v * v) - 3.0;
}

double ios_signal_median(IOSSignal* s) {
    if (!s || s->length < 1) return 0.0;
    double* copy = (double*)malloc((size_t)s->length * sizeof(double));
    if (!copy) return 0.0;
    memcpy(copy, s->data, (size_t)s->length * sizeof(double));
    for (int i = 0; i < s->length - 1; i++)
        for (int j = i + 1; j < s->length; j++)
            if (copy[j] < copy[i]) {
                double t = copy[i]; copy[i] = copy[j]; copy[j] = t;
            }
    double med = (s->length % 2 == 1) ? copy[s->length / 2]
               : (copy[s->length / 2 - 1] + copy[s->length / 2]) / 2.0;
    free(copy);
    return med;
}

/* ============================================================
 * Signal Correlation & Convolution
 * ============================================================ */
double ios_signal_inner_product(const IOSSignal* a, const IOSSignal* b) {
    if (!a || !b) return 0.0;
    int n = a->length < b->length ? a->length : b->length;
    double ip = 0.0;
    for (int i = 0; i < n; i++)
        ip += a->data[i] * b->data[i];
    return ip;
}

double ios_signal_corr(const IOSSignal* a, const IOSSignal* b, int lag) {
    if (!a || !b || a->length < 2) return 0.0;
    double ma = ios_signal_mean(a);
    double mb = ios_signal_mean(b);
    double num = 0.0, da = 0.0, db = 0.0;
    int n = a->length < b->length ? a->length : b->length;
    for (int i = 0; i < n; i++) {
        double ai = a->data[i] - ma;
        double bi = b->data[i] - mb;
        da += ai * ai;
        db += bi * bi;
        if (i + lag >= 0 && i + lag < n)
            num += ai * (b->data[i + lag] - mb);
    }
    double denom = sqrt(da * db);
    return (denom > IOS_EPS) ? num / denom : 0.0;
}

double ios_signal_cross_correlation(const IOSSignal* a, const IOSSignal* b) {
    return ios_signal_corr(a, b, 0);
}

IOSSignal* ios_signal_convolve(const IOSSignal* a, const IOSSignal* b) {
    if (!a || !b || a->length < 1 || b->length < 1) return NULL;
    IOSSignal* c = ios_signal_create(a->length + b->length);
    if (!c) return NULL;
    for (int i = 0; i < a->length; i++) {
        for (int j = 0; j < b->length; j++) {
            int k = i + j;
            while (c->length <= k) ios_signal_append(c, 0.0);
            c->data[k] += a->data[i] * b->data[j];
        }
    }
    return c;
}

/* ============================================================
 * Signal Generation (modifying existing IOSSignal*)
 * ============================================================ */
void ios_signal_fill(IOSSignal* s, double v) {
    if (s)
        for (int i = 0; i < s->length; i++)
            s->data[i] = v;
}

void ios_signal_ramp(IOSSignal* s, double slope) {
    if (s)
        for (int i = 0; i < s->length; i++)
            s->data[i] = slope * (double)i;
}

void ios_signal_sine(IOSSignal* s, double amp, double freq, double dt) {
    if (s)
        for (int i = 0; i < s->length; i++)
            s->data[i] = amp * sin(2.0 * IOS_PI * freq * (double)i * dt);
}

void ios_signal_step(IOSSignal* s, double amp, int start) {
    if (!s) return;
    for (int i = 0; i < s->length; i++)
        s->data[i] = (i >= start) ? amp : 0.0;
}

void ios_signal_random(IOSSignal* s, double amp) {
    if (s)
        for (int i = 0; i < s->length; i++)
            s->data[i] = amp * ((double)rand() / (double)RAND_MAX - 0.5) * 2.0;
}

/* ============================================================
 * Signal Operations
 * ============================================================ */
void ios_signal_differentiate(const IOSSignal* s, IOSSignal* ds, double dt) {
    if (!s || !ds || s->length < 2 || dt <= IOS_EPS) return;
    ios_signal_clear(ds);
    for (int i = 1; i < s->length; i++)
        ios_signal_append(ds, (s->data[i] - s->data[i - 1]) / dt);
}

void ios_signal_integrate(const IOSSignal* s, IOSSignal* is, double dt) {
    if (!s || !is || dt <= IOS_EPS) return;
    ios_signal_clear(is);
    double sum = 0.0;
    for (int i = 0; i < s->length; i++) {
        sum += s->data[i] * dt;
        ios_signal_append(is, sum);
    }
}

void ios_signal_normalize(IOSSignal* s) {
    if (!s || s->length < 2) return;
    double mn = ios_signal_min(s);
    double mx = ios_signal_max(s);
    double range = mx - mn;
    if (range < IOS_EPS) return;
    for (int i = 0; i < s->length; i++)
        s->data[i] = (s->data[i] - mn) / range;
}

void ios_signal_standardize(IOSSignal* s) {
    if (!s || s->length < 2) return;
    double m = ios_signal_mean(s);
    double v = ios_signal_variance(s);
    double std = sqrt(v);
    if (std < IOS_EPS) return;
    for (int i = 0; i < s->length; i++)
        s->data[i] = (s->data[i] - m) / std;
}

void ios_signal_moving_average(const IOSSignal* s, IOSSignal* out, int window) {
    if (!s || !out || window < 1 || window > s->length) return;
    ios_signal_clear(out);
    int half = window / 2;
    for (int i = 0; i < s->length; i++) {
        double sum = 0.0;
        int count = 0;
        for (int j = i - half; j <= i + half; j++) {
            if (j >= 0 && j < s->length) { sum += s->data[j]; count++; }
        }
        ios_signal_append(out, count > 0 ? sum / (double)count : s->data[i]);
    }
}

void ios_signal_exponential_smooth(const IOSSignal* s, IOSSignal* out, double alpha) {
    if (!s || !out || s->length < 1) return;
    ios_signal_clear(out);
    double smoothed = s->data[0];
    ios_signal_append(out, smoothed);
    for (int i = 1; i < s->length; i++) {
        smoothed = alpha * s->data[i] + (1.0 - alpha) * smoothed;
        ios_signal_append(out, smoothed);
    }
}

int ios_signal_clip(IOSSignal* s, double lo, double hi) {
    if (!s) return 0;
    int c = 0;
    for (int i = 0; i < s->length; i++) {
        if (s->data[i] < lo) { s->data[i] = lo; c++; }
        if (s->data[i] > hi) { s->data[i] = hi; c++; }
    }
    return c;
}

void ios_signal_sub(IOSSignal* a, const IOSSignal* b) {
    if (!a || !b) return;
    int n = a->length < b->length ? a->length : b->length;
    for (int i = 0; i < n; i++)
        a->data[i] -= b->data[i];
}

void ios_signal_resample(const IOSSignal* s, IOSSignal* out, int factor) {
    if (!s || !out || factor < 2) return;
    ios_signal_clear(out);
    for (int i = 0; i < s->length; i += factor)
        ios_signal_append(out, s->data[i]);
}

void ios_signal_interpolate(const IOSSignal* s, IOSSignal* out, int factor) {
    if (!s || !out || factor < 2) return;
    ios_signal_clear(out);
    for (int i = 0; i < s->length - 1; i++)
        for (int j = 0; j < factor; j++) {
            double t = (double)j / (double)factor;
            ios_signal_append(out, s->data[i] + t * (s->data[i + 1] - s->data[i]));
        }
}

int ios_signal_binarize(const IOSSignal* s, double th, IOSSignal* out) {
    if (!s || !out) return -1;
    ios_signal_clear(out);
    for (int i = 0; i < s->length; i++)
        ios_signal_append(out, s->data[i] > th ? 1.0 : 0.0);
    return 0;
}

/* ============================================================
 * Signal Analysis
 * ============================================================ */
int ios_signal_find_peak(const IOSSignal* s, double* peak_val, int* peak_idx) {
    if (!s || s->length < 1) return -1;
    *peak_val = s->data[0];
    *peak_idx = 0;
    for (int i = 1; i < s->length; i++)
        if (s->data[i] > *peak_val) {
            *peak_val = s->data[i];
            *peak_idx = i;
        }
    return 0;
}

int ios_find_crossings(const IOSSignal* s, double thresh, int* idx, int max_n) {
    if (!s || !idx || s->length < 2 || max_n < 1) return 0;
    int count = 0;
    for (int i = 1; i < s->length && count < max_n; i++) {
        if ((s->data[i - 1] - thresh) * (s->data[i] - thresh) < 0)
            idx[count++] = i;
    }
    return count;
}

int ios_signal_count_threshold(const IOSSignal* s, double th, int direction) {
    if (!s) return 0;
    int c = 0;
    for (int i = 1; i < s->length; i++) {
        if (direction > 0 && s->data[i - 1] <= th && s->data[i] > th) c++;
        if (direction < 0 && s->data[i - 1] >= th && s->data[i] < th) c++;
    }
    return c;
}

int ios_signal_zero_crossings(const IOSSignal* s) {
    if (!s || s->length < 2) return 0;
    int c = 0;
    for (int i = 1; i < s->length; i++)
        if (s->data[i - 1] * s->data[i] < 0) c++;
    return c;
}

double ios_signal_rise_time(const IOSSignal* s, double dt) {
    if (!s || s->length < 2) return -1.0;
    double fv = s->data[s->length - 1];
    if (fabs(fv) < IOS_EPS) return -1.0;
    for (int i = 0; i < s->length; i++)
        if (s->data[i] >= 0.9 * fv)
            return (double)i * dt;
    return -1.0;
}

double ios_signal_settling_time(const IOSSignal* s, double dt, double tol) {
    if (!s || s->length < 2) return -1.0;
    double fv = s->data[s->length - 1];
    if (fabs(fv) < IOS_EPS) return 0.0;
    for (int i = s->length - 1; i >= 0; i--)
        if (fabs(s->data[i] / fv - 1.0) > tol)
            return (double)i * dt;
    return 0.0;
}

double ios_signal_overshoot_pct(const IOSSignal* s) {
    if (!s || s->length < 2) return 0.0;
    double fv = s->data[s->length - 1];
    if (fabs(fv) < IOS_EPS) return 0.0;
    double pk = s->data[0];
    for (int i = 1; i < s->length; i++)
        if (s->data[i] > pk) pk = s->data[i];
    return fabs(fv) > IOS_EPS ? (pk / fv - 1.0) * 100.0 : 0.0;
}

double ios_signal_duty_cycle(const IOSSignal* s, double th) {
    if (!s || s->length < 1) return 0.0;
    int c = 0;
    for (int i = 0; i < s->length; i++)
        if (s->data[i] > th) c++;
    return (double)c / (double)s->length;
}

double ios_signal_snr(const IOSSignal* signal, const IOSSignal* noise) {
    if (!signal || !noise) return 0.0;
    double Ps = ios_signal_power(signal);
    double Pn = ios_signal_power(noise);
    return (Pn > IOS_EPS) ? 10.0 * log10(Ps / Pn) : 100.0;
}

/* ============================================================
 * State-Space Operations
 * ============================================================ */
IOSStateSpace* ios_ss_create(int n, int m, int p) {
    if (n < 1 || m < 1 || p < 1) return NULL;
    IOSStateSpace* sys = (IOSStateSpace*)calloc(1, sizeof(IOSStateSpace));
    if (!sys) return NULL;
    sys->n = n; sys->m = m; sys->p = p;
    sys->A = (double*)calloc((size_t)(n * n), sizeof(double));
    sys->B = (double*)calloc((size_t)(n * m), sizeof(double));
    sys->C = (double*)calloc((size_t)(p * n), sizeof(double));
    sys->D = (double*)calloc((size_t)(p * m), sizeof(double));
    if (!sys->A || !sys->B || !sys->C || !sys->D) {
        free(sys->A); free(sys->B); free(sys->C); free(sys->D); free(sys);
        return NULL;
    }
    return sys;
}

void ios_ss_free(IOSStateSpace* sys) {
    if (sys) { free(sys->A); free(sys->B); free(sys->C); free(sys->D); free(sys); }
}

void ios_ss_simulate(const IOSStateSpace* sys, const IOSSignal* u, IOSSignal* y) {
    if (!sys || !u || !y || sys->n < 1) return;
    ios_signal_clear(y);
    double* x = (double*)calloc((size_t)sys->n, sizeof(double));
    if (!x) return;
    for (int t = 0; t < u->length; t++) {
        double* xdot = (double*)calloc((size_t)sys->n, sizeof(double));
        if (!xdot) break;
        for (int i = 0; i < sys->n; i++) {
            for (int j = 0; j < sys->n; j++)
                xdot[i] += sys->A[i * sys->n + j] * x[j];
            for (int j = 0; j < sys->m; j++)
                xdot[i] += sys->B[i * sys->m + j] * u->data[t];
        }
        for (int i = 0; i < sys->n; i++)
            x[i] += xdot[i];
        free(xdot);
        for (int i = 0; i < sys->p; i++) {
            double out = 0.0;
            for (int j = 0; j < sys->n; j++)
                out += sys->C[i * sys->n + j] * x[j];
            for (int j = 0; j < sys->m; j++)
                out += sys->D[i * sys->m + j] * u->data[t];
            ios_signal_append(y, out);
        }
    }
    free(x);
}

void ios_ss_print(const IOSStateSpace* sys) {
    if (!sys) return;
    printf("SS[n=%d,m=%d,p=%d]\n", sys->n, sys->m, sys->p);
}

void ios_ss_copy(const IOSStateSpace* src, IOSStateSpace* dst) {
    if (!src || !dst || src->n != dst->n) return;
    memcpy(dst->A, src->A, (size_t)(src->n * src->n) * sizeof(double));
    memcpy(dst->B, src->B, (size_t)(src->n * src->m) * sizeof(double));
    memcpy(dst->C, src->C, (size_t)(src->p * src->n) * sizeof(double));
    memcpy(dst->D, src->D, (size_t)(src->p * src->m) * sizeof(double));
}

double ios_ss_eigen_max_real(const IOSStateSpace* sys) {
    if (!sys || sys->n < 1) return 0.0;
    double mx = sys->A[0];
    for (int i = 1; i < sys->n; i++) {
        double d = sys->A[i * sys->n + i];
        if (d > mx) mx = d;
    }
    return mx;
}

int ios_ss_is_siso(const IOSStateSpace* sys) {
    return (sys && sys->m == 1 && sys->p == 1) ? 1 : 0;
}

int ios_ss_is_mimo(const IOSStateSpace* sys) {
    return (sys && (sys->m > 1 || sys->p > 1)) ? 1 : 0;
}

/* ============================================================
 * Transfer Function Operations
 * ============================================================ */
IOSTransferFcn* ios_tf_create(int num_n, int den_n) {
    if (num_n < 1 || den_n < 1) return NULL;
    IOSTransferFcn* tf = (IOSTransferFcn*)calloc(1, sizeof(IOSTransferFcn));
    if (!tf) return NULL;
    tf->n_num = num_n; tf->n_den = den_n;
    tf->num = (double*)calloc((size_t)num_n, sizeof(double));
    tf->den = (double*)calloc((size_t)den_n, sizeof(double));
    if (!tf->num || !tf->den) {
        free(tf->num); free(tf->den); free(tf); return NULL;
    }
    return tf;
}

void ios_tf_free(IOSTransferFcn* tf) {
    if (tf) { free(tf->num); free(tf->den); free(tf); }
}

double ios_tf_eval(const IOSTransferFcn* tf, double s) {
    if (!tf) return 0.0;
    double num_re = 0.0, num_im = 0.0;
    double den_re = 0.0, den_im = 0.0;
    for (int i = 0; i < tf->n_num; i++) {
        int pwr = tf->n_num - 1 - i;
        if (pwr % 2 == 0) num_re += tf->num[i] * pow(s, pwr);
        else num_im += tf->num[i] * pow(s, pwr);
    }
    for (int i = 0; i < tf->n_den; i++) {
        int pwr = tf->n_den - 1 - i;
        if (pwr % 2 == 0) den_re += tf->den[i] * pow(s, pwr);
        else den_im += tf->den[i] * pow(s, pwr);
    }
    double den_mag_sq = den_re * den_re + den_im * den_im;
    if (den_mag_sq < IOS_EPS) return IOS_INF;
    return (num_re * den_re + num_im * den_im) / den_mag_sq;
}

int ios_tf_order(const IOSTransferFcn* tf) {
    return tf ? tf->n_den - 1 : -1;
}

double ios_tf_dcgain(const IOSTransferFcn* tf) {
    if (!tf || tf->n_den < 1 || tf->n_num < 1) return 0.0;
    /* DC gain = G(0) = b_{m-1} / a_{n-1} */
    double num_last = tf->num[tf->n_num - 1];
    double den_last = tf->den[tf->n_den - 1];
    return fabs(den_last) > IOS_EPS ? num_last / den_last : IOS_INF;
}

bool ios_tf_is_proper(const IOSTransferFcn* tf) {
    return tf && tf->n_den >= tf->n_num;
}

bool ios_tf_is_stable(const IOSTransferFcn* tf) {
    if (!tf || tf->n_den < 1) return false;
    if (tf->n_den == 1) return fabs(tf->den[0]) > IOS_EPS;
    if (tf->n_den == 2) {
        double a = tf->den[0], b = tf->den[1];
        if (fabs(a) < IOS_EPS) return false;
        return b / a > 0.0;
    }
    if (tf->n_den == 3) {
        double a0 = tf->den[0], a1 = tf->den[1], a2 = tf->den[2];
        if (fabs(a0) < IOS_EPS) return false;
        return a1 > 0.0 && a2 > 0.0 && a1 * a2 > a0;
    }
    double sign = tf->den[0];
    if (fabs(sign) < IOS_EPS) return false;
    for (int i = 1; i < tf->n_den; i++)
        if (tf->den[i] * sign <= 0.0) return false;
    return true;
}

void ios_print_tf(const IOSTransferFcn* tf) {
    if (!tf) { printf("TF[null]\n"); return; }
    printf("TF[%d/%d]: num=[", tf->n_num, tf->n_den);
    for (int i = 0; i < tf->n_num; i++) printf("%.3f ", tf->num[i]);
    printf("] den=[");
    for (int i = 0; i < tf->n_den; i++) printf("%.3f ", tf->den[i]);
    printf("]\n");
}

/* ============================================================
 * Frequency Response
 * ============================================================ */
double ios_frequency_response(const IOSTransferFcn* tf, double w, double* mag, double* phase) {
    if (!tf || !mag || !phase) return 0.0;
    double Re = 0.0, Im = 0.0;
    for (int i = 0; i < tf->n_num; i++) {
        int pwr = tf->n_num - 1 - i;
        double coeff = tf->num[i];
        if (pwr % 4 == 0) Re += coeff * pow(w, pwr);
        else if (pwr % 4 == 1) Im += coeff * pow(w, pwr);
        else if (pwr % 4 == 2) Re -= coeff * pow(w, pwr);
        else Im -= coeff * pow(w, pwr);
    }
    double den_Re = 0.0, den_Im = 0.0;
    for (int i = 0; i < tf->n_den; i++) {
        int pwr = tf->n_den - 1 - i;
        double coeff = tf->den[i];
        if (pwr % 4 == 0) den_Re += coeff * pow(w, pwr);
        else if (pwr % 4 == 1) den_Im += coeff * pow(w, pwr);
        else if (pwr % 4 == 2) den_Re -= coeff * pow(w, pwr);
        else den_Im -= coeff * pow(w, pwr);
    }
    double dn = den_Re * den_Re + den_Im * den_Im;
    if (dn < IOS_EPS) { *mag = 0.0; *phase = 0.0; return 0.0; }
    double G_re = (Re * den_Re + Im * den_Im) / dn;
    double G_im = (Im * den_Re - Re * den_Im) / dn;
    *mag = sqrt(G_re * G_re + G_im * G_im);
    *phase = atan2(G_im, G_re);
    return G_re;
}

void ios_bode_plot(const IOSTransferFcn* tf, double* w, double* mag, double* phase, int n) {
    if (!tf || !w || !mag || !phase || n < 1) return;
    for (int i = 0; i < n; i++) {
        w[i] = 0.01 * pow(10.0, 4.0 * (double)i / (double)n);
        ios_frequency_response(tf, w[i], &mag[i], &phase[i]);
    }
}

double ios_bandwidth(const IOSTransferFcn* tf) {
    if (!tf) return 0.0;
    double dc = ios_tf_dcgain(tf);
    if (fabs(dc) < IOS_EPS) return 0.0;
    for (int k = 0; k < IOS_FREQ_POINTS; k++) {
        double w = 0.01 * pow(10.0, 4.0 * (double)k / (double)IOS_FREQ_POINTS);
        double mag = 0.0, phase = 0.0;
        ios_frequency_response(tf, w, &mag, &phase);
        if (mag < fabs(dc) * 0.7071) return w;
    }
    return 100.0;
}

/* ============================================================
 * System-Level Operations (IOSSystem wrapper)
 * ============================================================ */
IOSSystem* ios_system_create(IOSStateSpace* ss) {
    if (!ss) return NULL;
    IOSSystem* sys = (IOSSystem*)calloc(1, sizeof(IOSSystem));
    if (!sys) return NULL;
    sys->sys = ss;
    sys->L2_gain = IOS_INF;
    sys->L1_gain = IOS_INF;
    sys->Linf_gain = IOS_INF;
    sys->bibo = IOS_BIBO_CRITICAL;
    sys->analyzed = false;
    return sys;
}

void ios_system_free(IOSSystem* sys) {
    if (sys) free(sys);
}

void ios_system_analyze(IOSSystem* sys) {
    if (!sys || !sys->sys) return;
    sys->bibo = ios_check_bibo_ss(sys->sys);
    sys->L2_gain = ios_L2_gain_ss(sys->sys);
    sys->L1_gain = ios_L1_induced_norm(sys->sys);
    sys->Linf_gain = ios_Linf_induced_norm(sys->sys);
    sys->passivity = ios_check_passivity_ss(sys->sys);
    sys->analyzed = true;
}

double ios_system_gain_margin(const IOSSystem* sys) {
    if (!sys || !sys->sys) return IOS_INF;
    return ios_system_norm_ratio(sys->sys);
}

void ios_print_system(const IOSSystem* sys) {
    if (!sys) return;
    printf("IOSystem: BIBO=%s L2=%.4f L1=%.4f Linf=%.4f analyzed=%s\n",
           sys->bibo == IOS_BIBO_STABLE ? "STABLE" :
           (sys->bibo == IOS_BIBO_UNSTABLE ? "UNSTABLE" : "CRITICAL"),
           sys->L2_gain, sys->L1_gain, sys->Linf_gain,
           sys->analyzed ? "Yes" : "No");
}

/* ============================================================
 * Convenience Signal Constructors (return new IOSSignal*)
 * ============================================================ */
IOSSignal* ios_signal_impulse(int len, double amp, int idx) {
    IOSSignal* s = ios_signal_create(len);
    if (s) {
        s->length = len;
        if (idx >= 0 && idx < len) s->data[idx] = amp;
    }
    return s;
}

IOSSignal* ios_signal_step_signal(int len, double amp, int start) {
    IOSSignal* s = ios_signal_create(len);
    if (s) {
        s->length = len;
        for (int i = start; i < len; i++)
            s->data[i] = amp;
    }
    return s;
}

IOSSignal* ios_signal_sine_wave(int len, double amp, double freq, double dt) {
    IOSSignal* s = ios_signal_create(len);
    if (s) {
        s->length = len;
        for (int i = 0; i < len; i++)
            s->data[i] = amp * sin(2.0 * IOS_PI * freq * (double)i * dt);
    }
    return s;
}

IOSSignal* ios_signal_random_noise(int len, double amp) {
    IOSSignal* s = ios_signal_create(len);
    if (s) {
        s->length = len;
        for (int i = 0; i < len; i++)
            s->data[i] = amp * ((double)rand() / (double)RAND_MAX - 0.5) * 2.0;
    }
    return s;
}

/* ============================================================
 * Transfer function helpers
 * ============================================================ */
double ios_tf_phase_at(const IOSTransferFcn* tf, double w) {
    double mag = 0.0, phase = 0.0;
    ios_frequency_response(tf, w, &mag, &phase);
    return phase;
}

double ios_tf_magnitude_at(const IOSTransferFcn* tf, double w) {
    double mag = 0.0, phase = 0.0;
    ios_frequency_response(tf, w, &mag, &phase);
    return mag;
}

void ios_tf_scale(IOSTransferFcn* tf, double k) {
    if (tf)
        for (int i = 0; i < tf->n_num; i++)
            tf->num[i] *= k;
}

void ios_tf_add(const IOSTransferFcn* a, const IOSTransferFcn* b, IOSTransferFcn* result) {
    (void)a; (void)b; (void)result;
}

void ios_tf_print_full(const IOSTransferFcn* tf) {
    ios_print_tf(tf);
}

/* ============================================================
 * System report and utility
 * ============================================================ */
void ios_system_report(const IOSSystem* s) {
    ios_print_system(s);
}

double ios_system_norm_ratio(const IOSStateSpace* sys) {
    if (!sys) return 0.0;
    double L2 = ios_L2_gain_ss(sys);
    double L1 = ios_L1_induced_norm(sys);
    return (L1 > IOS_EPS) ? L2 / L1 : 0.0;
}

double ios_dc_gain_ss(const IOSStateSpace* sys) {
    if (!sys || sys->n < 1) return 0.0;
    double g = 0.0;
    for (int i = 0; i < sys->p; i++)
        for (int j = 0; j < sys->m; j++) {
            double sum = sys->D[i * sys->m + j];
            for (int k = 0; k < sys->n; k++)
                sum -= sys->C[i * sys->n + k] * sys->B[k * sys->m + j]
                       / fmax(sys->A[k * sys->n + k], IOS_EPS);
            g += fabs(sum);
        }
    return g;
}

void ios_signal_all_stats(const IOSSignal* s) {
    if (!s) return;
    printf("Signal[%d]: min=%.4f max=%.4f mean=%.4f rms=%.4f energy=%.4f var=%.4f\n",
           s->length, ios_signal_min(s), ios_signal_max(s),
           ios_signal_mean(s), ios_signal_rms(s),
           ios_signal_energy(s), ios_signal_variance(s));
}
