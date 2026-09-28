/* Odd/even residual-vector bounds regression using the production public API.
 * Diagonal independent measurements have an analytic KF/Joseph solution.
 * Each legacy positional group has constant normalized residuals, so robust
 * reweighting is inactive. No phase/code ordering is assumed by these fixtures.
 */
#include <assert.h>
#include <math.h>
#include "rtklib.h"

static void diagonal_case(int m)
{
    double *x=mat(m,1),*P=zeros(m,m),*H=eye(m),*v=mat(m,1),*R=zeros(m,m);
    double expected;
    int i,j;
    /* Exactly representable normalized residuals avoid roundoff in the
     * existing group-variance calculation activating robust reweighting. */
    for(i=0;i<m;i++) {
        x[i]=2.0+i*0.125;
        P[i+i*m]=3.0;
        R[i+i*m]=1.0;
        v[i]=i%2?-0.25:0.5;
    }
    assert(filter_vbakf(x,P,H,v,R,m,m)==0);
    for(i=0;i<m;i++) {
        expected=2.0+i*0.125+0.75*v[i];
        assert(fabs(x[i]-expected)<1e-11);
        for(j=0;j<m;j++) {
            expected=i==j?0.75:0.0;
            assert(fabs(P[i+j*m]-expected)<1e-11);
        }
    }
    free(x);free(P);free(H);free(v);free(R);
    printf("PASS m=%d: finite state, analytic covariance, all residuals processed\n",m);
}

/* Capture nonuniform even-m results for comparison against the original code.
 * Includes a large residual to exercise the existing robust weight policy.
 */
static void even_fingerprint(int m)
{
    double x[3]={1.0,2.0,3.0},P[9]={4.0,0,0,0,5.0,0,0,0,6.0};
    double *H=zeros(3,m),*v=mat(m,1),*R=zeros(m,m);
    int i;
    for(i=0;i<m;i++) {
        H[i%3+i*3]=1.0;
        R[i+i*m]=3.0+(i%4);
        v[i]=(i%7)-3.0;
    }
    v[m-1]=50.0;
    assert(filter_vbakf(x,P,H,v,R,3,m)==0);
    printf("m=%d",m);
    for(i=0;i<3;i++) printf(" %.17g",x[i]);
    for(i=0;i<9;i++) printf(" %.17g",P[i]);
    putchar('\n');
    free(H);free(v);free(R);
}

int main(int argc,char **argv)
{
    int sizes[]={1,2,3,20,21,24,95,96};
    unsigned int i;
    if(argc==2&&!strcmp(argv[1],"--even-fingerprint")) {
        even_fingerprint(20);even_fingerprint(24);return 0;
    }
    for(i=0;i<sizeof(sizes)/sizeof(sizes[0]);i++) diagonal_case(sizes[i]);
    return 0;
}
