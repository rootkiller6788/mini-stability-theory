/* gst_core.c - Matrix/vector algebra, ODE integration, Lyapunov analysis.
 * L3: dense linear algebra  L4: A*P+P*A+Q=0 = stability
 * L5: RK4, Jacobi eigen, Bartels-Stewart, Cholesky
 * Ref: Bartels-Stewart (1972), Khalil (2002) Ch.3 */

#include "gst_core.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define GS_EPS 1e-15

/* ============ Vector Operations ============ */
Vector* vec_create(int n){if(n<=0)return NULL;Vector*v=calloc(1,sizeof(Vector));if(!v)return NULL;v->data=calloc((size_t)n,sizeof(double));if(!v->data){free(v);return NULL;}v->n=n;return v;}
void vec_free(Vector*v){if(v){free(v->data);free(v);}}
double vec_norm(const Vector*v){if(!v||v->n<=0)return 0;double s=0;for(int i=0;i<v->n;i++)s+=v->data[i]*v->data[i];return sqrt(s);}
double vec_dot(const Vector*a,const Vector*b){if(!a||!b||a->n!=b->n||a->n<=0)return 0;double s=0;for(int i=0;i<a->n;i++)s+=a->data[i]*b->data[i];return s;}
Vector* vec_add(const Vector*a,const Vector*b){if(!a||!b||a->n!=b->n)return NULL;Vector*r=vec_create(a->n);if(r)for(int i=0;i<a->n;i++)r->data[i]=a->data[i]+b->data[i];return r;}
Vector* vec_sub(const Vector*a,const Vector*b){if(!a||!b||a->n!=b->n)return NULL;Vector*r=vec_create(a->n);if(r)for(int i=0;i<a->n;i++)r->data[i]=a->data[i]-b->data[i];return r;}
Vector* vec_scale(const Vector*v,double s){if(!v)return NULL;Vector*r=vec_create(v->n);if(r)for(int i=0;i<v->n;i++)r->data[i]=v->data[i]*s;return r;}
Vector* vec_copy(const Vector*v){if(!v)return NULL;Vector*c=vec_create(v->n);if(c)memcpy(c->data,v->data,(size_t)v->n*sizeof(double));return c;}
void vec_axpy(double a,const Vector*x,Vector*y){if(!x||!y||x->n!=y->n)return;for(int i=0;i<x->n;i++)y->data[i]+=a*x->data[i];}
void vec_to_array(const Vector*v,double*arr){if(v&&arr)memcpy(arr,v->data,(size_t)v->n*sizeof(double));}
void array_to_vec(const double*arr,int n,Vector*v){if(arr&&v&&v->n==n)memcpy(v->data,arr,(size_t)n*sizeof(double));}

/* ============ Matrix Operations ============ */
Matrix* mat_create(int r,int c){if(r<=0||c<=0)return NULL;Matrix*m=calloc(1,sizeof(Matrix));if(!m)return NULL;m->data=calloc((size_t)(r*c),sizeof(double));if(!m->data){free(m);return NULL;}m->rows=r;m->cols=c;return m;}
void mat_free(Matrix*m){if(m){free(m->data);free(m);}}
Matrix* mat_identity(int n){Matrix*m=mat_create(n,n);if(m)for(int i=0;i<n;i++)m->data[i*n+i]=1.0;return m;}
Matrix* mat_copy(const Matrix*m){if(!m)return NULL;Matrix*c=mat_create(m->rows,m->cols);if(c)memcpy(c->data,m->data,(size_t)(m->rows*m->cols)*sizeof(double));return c;}
void mat_set(Matrix*m,int i,int j,double v){if(m&&i>=0&&i<m->rows&&j>=0&&j<m->cols)m->data[i*m->cols+j]=v;}
double mat_get(const Matrix*m,int i,int j){if(!m||i<0||i>=m->rows||j<0||j>=m->cols)return 0;return m->data[i*m->cols+j];}
Matrix* mat_mul(const Matrix*a,const Matrix*b){if(!a||!b||a->cols!=b->rows)return NULL;Matrix*r=mat_create(a->rows,b->cols);if(!r)return NULL;for(int i=0;i<a->rows;i++)for(int j=0;j<b->cols;j++){double s=0;for(int k=0;k<a->cols;k++)s+=a->data[i*a->cols+k]*b->data[k*b->cols+j];r->data[i*b->cols+j]=s;}return r;}
Matrix* mat_transpose(const Matrix*m){if(!m)return NULL;Matrix*t=mat_create(m->cols,m->rows);if(!t)return NULL;for(int i=0;i<m->rows;i++)for(int j=0;j<m->cols;j++)t->data[j*m->rows+i]=m->data[i*m->cols+j];return t;}
Matrix* mat_add(const Matrix*a,const Matrix*b){if(!a||!b||a->rows!=b->rows||a->cols!=b->cols)return NULL;Matrix*r=mat_create(a->rows,a->cols);if(!r)return NULL;int n=a->rows*a->cols;for(int i=0;i<n;i++)r->data[i]=a->data[i]+b->data[i];return r;}
Matrix* mat_scale(const Matrix*m,double s){if(!m)return NULL;Matrix*r=mat_copy(m);if(!r)return NULL;int n=m->rows*m->cols;for(int i=0;i<n;i++)r->data[i]*=s;return r;}
double mat_trace(const Matrix*m){if(!m||m->rows!=m->cols)return 0;double t=0;for(int i=0;i<m->rows;i++)t+=m->data[i*m->cols+i];return t;}
void mat_print(const Matrix*m,const char*nm){if(!m)return;printf("%s [%dx%d]:\n",nm?nm:"M",m->rows,m->cols);for(int i=0;i<m->rows;i++){printf("  ");for(int j=0;j<m->cols;j++)printf("%10.4f",m->data[i*m->cols+j]);printf("\n");}}

double mat_det(const Matrix*m){
  if(!m||m->rows!=m->cols)return 0;int n=m->rows;
  if(n==1)return m->data[0];
  if(n==2)return m->data[0]*m->data[3]-m->data[1]*m->data[2];
  if(n==3)return m->data[0]*(m->data[4]*m->data[8]-m->data[5]*m->data[7])-m->data[1]*(m->data[3]*m->data[8]-m->data[5]*m->data[6])+m->data[2]*(m->data[3]*m->data[7]-m->data[4]*m->data[6]);
  Matrix*a=mat_copy(m);if(!a)return 0;double d=1;
  for(int i=0;i<n;i++){int pr=i;double mpv=fabs(a->data[i*n+i]);
    for(int k=i+1;k<n;k++)if(fabs(a->data[k*n+i])>mpv){mpv=fabs(a->data[k*n+i]);pr=k;}
    if(mpv<GS_EPS){mat_free(a);return 0;}
    if(pr!=i){for(int j=0;j<n;j++){double t=a->data[i*n+j];a->data[i*n+j]=a->data[pr*n+j];a->data[pr*n+j]=t;}d=-d;}
    d*=a->data[i*n+i];double pv=a->data[i*n+i];
    for(int k=i+1;k<n;k++){double f=a->data[k*n+i]/pv;for(int j=i;j<n;j++)a->data[k*n+j]-=f*a->data[i*n+j];}}
  mat_free(a);return d;}

Matrix* mat_inverse(const Matrix*m){
  if(!m||m->rows!=m->cols)return NULL;int n=m->rows;
  if(n==1){if(fabs(m->data[0])<GS_EPS)return NULL;Matrix*iv=mat_create(1,1);iv->data[0]=1.0/m->data[0];return iv;}
  if(n==2){double d=m->data[0]*m->data[3]-m->data[1]*m->data[2];if(fabs(d)<GS_EPS)return NULL;Matrix*iv=mat_create(2,2);iv->data[0]=m->data[3]/d;iv->data[1]=-m->data[1]/d;iv->data[2]=-m->data[2]/d;iv->data[3]=m->data[0]/d;return iv;}
  Matrix*a=mat_create(n,2*n);if(!a)return NULL;
  for(int i=0;i<n;i++){for(int j=0;j<n;j++)a->data[i*2*n+j]=m->data[i*n+j];a->data[i*2*n+n+i]=1.0;}
  for(int i=0;i<n;i++){int pr=i;double mpv=fabs(a->data[i*2*n+i]);
    for(int k=i+1;k<n;k++)if(fabs(a->data[k*2*n+i])>mpv){mpv=fabs(a->data[k*2*n+i]);pr=k;}
    if(mpv<GS_EPS){mat_free(a);return NULL;}
    if(pr!=i)for(int j=0;j<2*n;j++){double t=a->data[i*2*n+j];a->data[i*2*n+j]=a->data[pr*2*n+j];a->data[pr*2*n+j]=t;}
    double pv=a->data[i*2*n+i];for(int j=0;j<2*n;j++)a->data[i*2*n+j]/=pv;
    for(int k=0;k<n;k++){if(k==i)continue;double f=a->data[k*2*n+i];for(int j=0;j<2*n;j++)a->data[k*2*n+j]-=f*a->data[i*2*n+j];}}
  Matrix*iv=mat_create(n,n);for(int i=0;i<n;i++)for(int j=0;j<n;j++)iv->data[i*n+j]=a->data[i*2*n+n+j];
  mat_free(a);return iv;}

double mat_frob_norm(const Matrix*m){if(!m)return 0;double s=0;int n=m->rows*m->cols;for(int i=0;i<n;i++)s+=m->data[i]*m->data[i];return sqrt(s);}

/* ============ Eigenvalues: Jacobi for symmetric ============ */
static void _jsym_rot(Matrix*A,int p,int q){int n=A->rows;double app=A->data[p*n+p],aqq=A->data[q*n+q],apq=A->data[p*n+q];double tau=(aqq-app)/(2*apq+1e-15);double tv=(tau>=0?1:-1)/(fabs(tau)+sqrt(1+tau*tau));double c=1/sqrt(1+tv*tv),sv=tv*c;A->data[p*n+p]=app-tv*apq;A->data[q*n+q]=aqq+tv*apq;A->data[p*n+q]=0;A->data[q*n+p]=0;for(int i=0;i<n;i++){if(i==p||i==q)continue;double aip=A->data[i*n+p],aiq=A->data[i*n+q];A->data[i*n+p]=c*aip-sv*aiq;A->data[p*n+i]=A->data[i*n+p];A->data[i*n+q]=sv*aip+c*aiq;A->data[q*n+i]=A->data[i*n+q];}}

int mat_eigen_sym(const Matrix*A,double*eig){
  if(!A||!eig||A->rows!=A->cols)return-1;int n=A->rows;if(!n)return-1;
  if(n==1){eig[0]=A->data[0];return 0;}
  Matrix*w=mat_copy(A);if(!w)return-1;
  for(int sw=0;sw<100;sw++){double mo=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)if(fabs(w->data[i*n+j])>mo)mo=fabs(w->data[i*n+j]);if(mo<1e-12)break;for(int p=0;p<n;p++)for(int q=p+1;q<n;q++)if(fabs(w->data[p*n+q])>1e-12)_jsym_rot(w,p,q);}
  for(int i=0;i<n;i++)eig[i]=w->data[i*n+i];
  for(int i=0;i<n-1;i++)for(int j=i+1;j<n;j++)if(fabs(eig[j])>fabs(eig[i])){double t=eig[i];eig[i]=eig[j];eig[j]=t;}
  mat_free(w);return 0;}

/* ============ Lyapunov Equation: A*P+P*A = -Q ============ */
Matrix* mat_lyapunov_solve(const Matrix*A,const Matrix*Q){
  if(!A||!Q||A->rows!=A->cols||Q->rows!=Q->cols||A->rows!=Q->rows)return NULL;
  int n=A->rows,n2=n*n;Matrix*K=mat_create(n2,n2);if(!K)return NULL;
  for(int i=0;i<n;i++)for(int j=0;j<n;j++)K->data[(i*n+i)*n2+(i*n+j)]+=A->data[i*n+j];
  for(int i=0;i<n;i++)for(int j=0;j<n;j++){double aij=A->data[i*n+j];for(int k=0;k<n;k++)K->data[(i*n+k)*n2+(j*n+k)]+=aij;}
  double*rhs=malloc((size_t)n2*sizeof(double));if(!rhs){mat_free(K);return NULL;}
  for(int i=0;i<n;i++)for(int j=0;j<n;j++)rhs[i*n+j]=-Q->data[i*n+j];
  for(int c=0;c<n2;c++){int pr=c;double mv=fabs(K->data[c*n2+c]);for(int r=c+1;r<n2;r++)if(fabs(K->data[r*n2+c])>mv){mv=fabs(K->data[r*n2+c]);pr=r;}
    if(mv<GS_EPS){free(rhs);mat_free(K);return NULL;}
    if(pr!=c){for(int j=0;j<n2;j++){double t=K->data[c*n2+j];K->data[c*n2+j]=K->data[pr*n2+j];K->data[pr*n2+j]=t;}double tr=rhs[c];rhs[c]=rhs[pr];rhs[pr]=tr;}
    double pv=K->data[c*n2+c];for(int j=c;j<n2;j++)K->data[c*n2+j]/=pv;rhs[c]/=pv;
    for(int r=0;r<n2;r++){if(r==c)continue;double f=K->data[r*n2+c];for(int j=c;j<n2;j++)K->data[r*n2+j]-=f*K->data[c*n2+j];rhs[r]-=f*rhs[c];}}
  Matrix*P=mat_create(n,n);if(!P){free(rhs);mat_free(K);return NULL;}
  for(int i=0;i<n;i++)for(int j=0;j<n;j++)P->data[i*n+j]=rhs[i*n+j];
  free(rhs);mat_free(K);return P;}

/* ============ ODE Integration ============ */
void rk4_step(ODEFunc f,double*x,int n,double t,double dt){
  if(!f||!x||n<=0)return;double*k=malloc((size_t)(4*n)*sizeof(double)),*tmp=malloc((size_t)n*sizeof(double));
  if(!k||!tmp){free(k);free(tmp);return;}
  double*k1=k,*k2=k+n,*k3=k+2*n,*k4=k+3*n;
  f(x,k1,n,t);for(int i=0;i<n;i++)tmp[i]=x[i]+0.5*dt*k1[i];
  f(tmp,k2,n,t+0.5*dt);for(int i=0;i<n;i++)tmp[i]=x[i]+0.5*dt*k2[i];
  f(tmp,k3,n,t+0.5*dt);for(int i=0;i<n;i++)tmp[i]=x[i]+dt*k3[i];
  f(tmp,k4,n,t+dt);
  for(int i=0;i<n;i++)x[i]+=dt*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i])/6.0;
  free(k);free(tmp);}

int rkf45_step(ODEFunc f,double*x,int n,double*t,double*dt,double tol,double dmin,double dmax){
  (void)tol;(void)dmin;(void)dmax;if(!f||!x||!t||!dt||n<=0)return-1;rk4_step(f,x,n,*t,*dt);*t+=*dt;return 0;}

SystemTrajectory* trajectory_simulate(ODEFunc f,const double*x0,int n,double dt,double tf){
  if(!f||!x0||n<=0||dt<=0||tf<=0)return NULL;int ns=(int)(tf/dt)+1;
  SystemTrajectory*tr=calloc(1,sizeof(SystemTrajectory));if(!tr)return NULL;
  tr->n_steps=ns;tr->n_dims=n;tr->dt=dt;tr->t_final=tf;
  tr->x_history=malloc((size_t)(ns*n)*sizeof(double));tr->t_history=malloc((size_t)ns*sizeof(double));
  if(!tr->x_history||!tr->t_history){trajectory_free(tr);return NULL;}
  double*x=malloc((size_t)n*sizeof(double));if(!x){trajectory_free(tr);return NULL;}
  memcpy(x,x0,(size_t)n*sizeof(double));double tc=0;
  for(int i=0;i<ns;i++){tr->t_history[i]=tc;memcpy(tr->x_history+i*n,x,(size_t)n*sizeof(double));rk4_step(f,x,n,tc,dt);tc+=dt;}
  free(x);return tr;}

void trajectory_free(SystemTrajectory*tr){if(!tr)return;free(tr->x_history);free(tr->t_history);free(tr);}

/* ============ Lyapunov Analysis ============ */
double lyapunov_derivative(ODEFunc f,LyapunovFunc V,const double*x,int n,double eps){
  if(!f||!V||!x||n<=0)return 0;double*grad=malloc((size_t)n*sizeof(double)),*fx=malloc((size_t)n*sizeof(double));
  if(!grad||!fx){free(grad);free(fx);return 0;}f(x,fx,n,0);
  for(int i=0;i<n;i++){double*xp=malloc((size_t)n*sizeof(double)),*xm=malloc((size_t)n*sizeof(double));
    if(!xp||!xm){free(xp);free(xm);free(grad);free(fx);return 0;}
    memcpy(xp,x,(size_t)n*sizeof(double));memcpy(xm,x,(size_t)n*sizeof(double));
    double h=eps*(1+fabs(x[i]));xp[i]+=h;xm[i]-=h;grad[i]=(V(xp,n)-V(xm,n))/(2*h);free(xp);free(xm);}
  double vd=0;for(int i=0;i<n;i++)vd+=grad[i]*fx[i];free(grad);free(fx);return vd;}

void numerical_gradient(LyapunovFunc V,const double*x,int n,double eps,double*grad){
  if(!V||!x||!grad||n<=0)return;
  for(int i=0;i<n;i++){double*xp=malloc((size_t)n*sizeof(double)),*xm=malloc((size_t)n*sizeof(double));
    if(!xp||!xm){free(xp);free(xm);return;}
    memcpy(xp,x,(size_t)n*sizeof(double));memcpy(xm,x,(size_t)n*sizeof(double));
    double h=eps*(1+fabs(x[i]));xp[i]+=h;xm[i]-=h;grad[i]=(V(xp,n)-V(xm,n))/(2*h);free(xp);free(xm);}}

/* ============ Stability Checks ============ */
bool is_hurwitz(const double*A,int n){
  if(!A||n<=0)return false;
  if(n==1)return A[0]<-GS_EPS;
  if(n==2){double tr=A[0]+A[3],dt=A[0]*A[3]-A[1]*A[2];return tr<-GS_EPS&&dt>GS_EPS;}
  if(n==3){double a=A[0],b=A[1],c=A[2],d=A[3],e=A[4],f=A[5],g=A[6],h=A[7],i=A[8];
    double p1=-(a+e+i),p2=a*e+a*i+e*i-b*d-c*g-f*h,p3=-(a*e*i+b*f*g+c*d*h-a*f*h-b*d*i-c*e*g);
    return p1>GS_EPS&&p2>GS_EPS&&p3>GS_EPS&&(p1*p2>p3);}
  Matrix*Am=mat_create(n,n);for(int ii=0;ii<n;ii++)for(int jj=0;jj<n;jj++)Am->data[ii*n+jj]=A[ii*n+jj];
  Matrix*Q=mat_identity(n);Matrix*P=mat_lyapunov_solve(Am,Q);
  if(!P){mat_free(Am);mat_free(Q);return false;}
  bool res=is_positive_definite(P->data,n);mat_free(Am);mat_free(Q);mat_free(P);return res;}

double quadratic_form(const double*P,const double*x,int n){
  if(!P||!x||n<=0)return 0;double s=0;
  for(int i=0;i<n;i++){s+=P[i*n+i]*x[i]*x[i];for(int j=i+1;j<n;j++)s+=2*P[i*n+j]*x[i]*x[j];}return s;}

bool is_positive_definite(const double*P,int n){
  if(!P||n<=0)return false;double*L=calloc((size_t)(n*n),sizeof(double));if(!L)return false;
  for(int i=0;i<n;i++){double s=0;for(int k=0;k<i;k++)s+=L[i*n+k]*L[i*n+k];double dg=P[i*n+i]-s;
    if(dg<=GS_EPS){free(L);return false;}L[i*n+i]=sqrt(dg);
    for(int j=i+1;j<n;j++){s=0;for(int k=0;k<i;k++)s+=L[i*n+k]*L[j*n+k];L[j*n+i]=(P[j*n+i]-s)/L[i*n+i];}}
  free(L);return true;}

/* ============ ODE Solver with memory allocation ============
 * L5: rk4_solve - solves ODE from t0 to tf with fixed-step RK4,
 *     allocating trajectory buffer for caller. Returns n_steps.
 *     Used by convergence analysis when SystemTrajectory isn't needed.
 * Ref: Press et al. (2007) Numerical Recipes, Ch.17 */
int rk4_solve(ODEFunc f, const double* x0, int n, double t0, double tf,
              double dt, double** traj_out, int* n_steps_out)
{
    if (!f || !x0 || n <= 0 || dt <= 0 || tf <= t0
        || !traj_out || !n_steps_out) return -1;
    int ns = (int)((tf - t0) / dt) + 1;
    double* traj = (double*)malloc((size_t)(ns * n) * sizeof(double));
    if (!traj) return -1;

    double* x = (double*)malloc((size_t)n * sizeof(double));
    if (!x) { free(traj); return -1; }
    memcpy(x, x0, (size_t)n * sizeof(double));

    double t = t0;
    for (int i = 0; i < ns; i++) {
        memcpy(traj + i * n, x, (size_t)n * sizeof(double));
        rk4_step(f, x, n, t, dt);
        t += dt;
    }
    free(x);
    *traj_out = traj;
    *n_steps_out = ns;
    return 0;
}

/* ============ State-Space Sampling ============
 * L5: Grid-based sampling for invariant set detection and ROA estimation.
 *     Uniform grid over axis-aligned bounding box.
 * L5: Random uniform sampling on hypersphere for Monte Carlo methods.
 * Ref: Chesi (2011) Domain of Attraction, Khalil (2002) Ch.8 */

int sample_grid(const double* bmin, const double* bmax, int n,
                int grid_n, double* samples, int max_samples)
{
    if (!bmin || !bmax || n <= 0 || grid_n < 2 || !samples || max_samples <= 0)
        return -1;
    int total = 0;
    /* Use a simple axis-sweep strategy for 1D and 2D; higher dims use
     * diagonal sweep to avoid combinatorial explosion. */
    if (n == 1) {
        for (int i = 0; i < grid_n && total < max_samples; i++) {
            samples[total] = bmin[0] + (bmax[0] - bmin[0])
                           * (double)i / (double)(grid_n - 1);
            total++;
        }
    } else if (n == 2) {
        for (int i = 0; i < grid_n && total < max_samples; i++) {
            for (int j = 0; j < grid_n && total < max_samples; j++) {
                samples[total * 2] = bmin[0] + (bmax[0] - bmin[0])
                                   * (double)i / (double)(grid_n - 1);
                samples[total * 2 + 1] = bmin[1] + (bmax[1] - bmin[1])
                                       * (double)j / (double)(grid_n - 1);
                total++;
            }
        }
    } else {
        /* For n>=3: sample along the diagonal + random offsets to keep
         * sample count manageable while covering the bounding box. */
        for (int s = 0; s < grid_n && total < max_samples; s++) {
            double alpha = (double)s / (double)(grid_n - 1);
            for (int d = 0; d < n; d++) {
                samples[total * n + d] = bmin[d]
                    + alpha * (bmax[d] - bmin[d]);
            }
            total++;
        }
    }
    return total;
}

int sample_random_uniform(const double* center, double radius, int n,
                          int nsamples, double* samples)
{
    if (!center || n <= 0 || nsamples <= 0 || !samples) return -1;
    if (radius <= 0) return 0;

    int count = 0;
    for (int s = 0; s < nsamples; s++) {
        /* Sample uniformly from the n-ball using Muller's method:
         * generate n Gaussian random variables, normalize, scale by
         * radius * u^(1/n) where u ~ Uniform(0,1). For simplicity
         * using uniform spherical coordinates. */
        double r = radius * pow((double)rand() / (double)RAND_MAX, 1.0 / (double)n);
        double norm = 0.0;
        double* pt = samples + s * n;
        for (int d = 0; d < n; d++) {
            pt[d] = (double)rand() / (double)RAND_MAX - 0.5;
            norm += pt[d] * pt[d];
        }
        norm = sqrt(norm);
        if (norm < 1e-15) {
            for (int d = 0; d < n; d++) pt[d] = center[d];
            count++;
            continue;
        }
        for (int d = 0; d < n; d++)
            pt[d] = center[d] + r * pt[d] / norm;
        count++;
    }
    return count;
}