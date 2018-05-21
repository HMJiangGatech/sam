#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "R.h"
#include "math.h"
#include "R_ext/BLAS.h"
#include "R_ext/Lapack.h"

void get_ols_norm(double *sub_X, double *y, int*nn, int *pp, double *ols_norm)
{
    int i,k,n,p;
    n = *nn;
    p = *pp;
    double *ols;
    ols = (double *) malloc(p*sizeof(double));
    
    *ols_norm = 0;
    for(k=0;k<p;k++){
        ols[k] = 0;
        for(i=0;i<n;i++)
            ols[k] = ols[k] + sub_X[k*n+i]*y[i];
        *ols_norm = *ols_norm + pow(ols[k],2);
    }
    *ols_norm = pow(*ols_norm,0.5);
    free(ols);
}

void get_sol_ext(double *sub_w, double *llambda, double *sub_X, double *r, int *nn, int *pp, double *w_norm)
{
    int n,p;
    double lambda;
    n = *nn;
    p = *pp;
    lambda = *llambda;
    
    int i,k;
    double *ols_w;
    
    ols_w = (double *) malloc(p*sizeof(double));
    double ratio;
    
    *w_norm = 0;
    for(k=0;k<p;k++)
    {
        ols_w[k] = 0;
        for(i=0;i<n;i++)
            ols_w[k] = ols_w[k] + sub_X[k*n+i]*r[i];
        
        *w_norm = *w_norm + pow(ols_w[k],2);
    }
    
    *w_norm = pow(*w_norm,0.5);
    
    //for(k=0;k<p;k++)
    //{
    //    printf("%f\t",ols_w[k]);
    //}
    //printf("\n");
    
    //printf("%f\n",*w_norm);
    
    if(*w_norm<=lambda){
        for(k=0;k<p;k++)
            sub_w[k] = 0;
        *w_norm = 0;
    }
    else
    {
        ratio = (*w_norm - lambda)/(*w_norm);
        for(k=0;k<p;k++)
            sub_w[k] = ols_w[k]*ratio;
        
        *w_norm = *w_norm - lambda;
        
        for(k=0;k<p;k++)
            for(i=0;i<n;i++){
                r[i] = r[i] - sub_X[k*n+i]*sub_w[k];
            }
    }
    
    
    /*
    for(k=0;k<p;k++)
    {
        printf("%f\t",sub_w[k]);
    }
    printf("\n");
    
    for(i=0;i<n;i++)
    {
        printf("%f\t",r[i]);
    }
    printf("\n");
     */
    free(ols_w);
}

void get_sol_int(double *sub_w, double *llambda, double *sub_X, double *r, int *nn, int *pp, double *w_norm, double *gap_tmp)
{
    int n,p;
    double lambda;
    n = *nn;
    p = *pp;
    lambda = *llambda;
    
    int i,k;
    double *ols_w;
    double tmp;
    ols_w = (double *) malloc(p*sizeof(double));
    double ratio;
    
    //for(k=0;k<p;k++)
    //{
    //    printf("%f\t",sub_w[k]);
    //}
    //printf("\n");
    
    
    double *sub_w_prev;
    sub_w_prev = (double *) malloc(p*sizeof(double));
    for(k=0;k<p;k++)
        sub_w_prev[k] = sub_w[k];
    
    for(k=0;k<p;k++)
        for(i=0;i<n;i++){
            r[i] = r[i] + sub_X[k*n+i]*sub_w[k];
        }
    
    *w_norm = 0;
    for(k=0;k<p;k++)
    {
        ols_w[k] = 0;
        for(i=0;i<n;i++)
            ols_w[k] = ols_w[k] + sub_X[k*n+i]*r[i];
        
        *w_norm = *w_norm + pow(ols_w[k],2);
    }
    
    //for(k=0;k<p;k++)
    //{
    //    printf("%f\t",ols_w[k]);
    //}
    //printf("\n");
    
    *w_norm = pow(*w_norm,0.5);
    
    //printf("%f\n",*w_norm);
    
    if(*w_norm<lambda){
        for(k=0;k<p;k++)
            sub_w[k] = 0;
        *w_norm = 0;
    }
    else
    {
        ratio = (*w_norm - lambda)/(*w_norm);
        for(k=0;k<p;k++)
            sub_w[k] = ols_w[k]*ratio;
        
        *w_norm = *w_norm - lambda;
        
        for(k=0;k<p;k++)
            for(i=0;i<n;i++){
                r[i] = r[i] - sub_X[k*n+i]*sub_w[k];
            }
    }
    
    
    tmp = 0;
    for(k=0;k<p;k++)
        tmp = tmp + pow(sub_w[k] - sub_w_prev[k],2);
    *gap_tmp = *gap_tmp + pow(tmp,0.5);
    
    /*
    for(k=0;k<p;k++)
      {
    printf("%f\t",sub_w[k]);
      }
    printf("\n");
     
    for(i=0;i<n;i++)
    {
      printf("%f\t",r[i]);
      }
    printf("\n");
    */
    free(ols_w);
    free(sub_w_prev);
}

void get_sol_final(double *sub_w, double *sub_w_final, double *sub_V, int *pp)
{
    int k,s,p;
    p = *pp;
    
    for(k=0;k<p;k++)
    {
        for(s=0;s<p;s++)
            sub_w_final[k] = sub_w_final[k] + sub_V[k*p+s]*sub_w[s];
    }
}

void grplasso(double *y, double *X, double *lambda, int *nnlambda, int *nn, int *dd, int *pp, double *ww, int *mmax_ite, double *tthol, int *iinput, int *df, double *sse, double *func_norm)
{
    
    int counter,n,d,p,m,max_ite,nlambda;
    int ite_ext,ite_int;
    int i,j,k,s;
    int input;
    
    
    int gap_ext,change_ext,back;
    double ilambda,thol,gap_int,gap_tmp0,gap_tmp1;
    double lambda_max,ols_norm;
    double w_norm;
    double *aw_norm;
    double *ols;
    
    nlambda = *nnlambda;
    n = *nn;
    d = *dd;
    p = *pp;
    m = d*p;
    max_ite = *mmax_ite;
    thol = *tthol;
    input = *iinput;
    
    double *U0, *D0, *V;
    double *w;
    
    D0 = (double *) malloc(p*sizeof(double));
    ols = (double *) malloc(p*sizeof(double));
    V = (double *) malloc((d*p*p)*sizeof(double));
    w = (double *) malloc(m*sizeof(double));
    
    char JOBU,JOBVT;
    int LDA,LDU,LDVT,LWORK,INFO;
    int *IWORK;
    double *WORK;
    JOBU = 'O';
    JOBVT = 'S';
    LDA = n;
    LDU = n;
    LDVT = p;
    
    lambda_max = 0;
    
    for(j=0;j<d;j++){
        LWORK = -1;
        WORK = (double *) malloc(1*sizeof(double));
        dgesvd_(&JOBU, &JOBVT, &n, &p, X+j*p*n, &LDA, D0, U0, &LDU, V+j*p*p, &LDVT, WORK, &LWORK, &INFO);
        
        LWORK = *WORK;
        free(WORK);
        WORK = (double *) malloc(LWORK*sizeof(double));
        
        dgesvd_(&JOBU, &JOBVT, &n, &p, X+j*p*n, &LDA, D0, U0, &LDU, V+j*p*p, &LDVT, WORK, &LWORK, &INFO);
        
        for(k=0;k<p;k++)
            for(s=0;s<p;s++)
                V[j*p*p+k*p+s] = V[j*p*p+k*p+s]/D0[s];
        /*
        for(k=0;k<p;k++)
            printf("%f\t",D0[k]);
        printf("\n");
        printf("\n");
        
        for(i=0;i<n;i++){
            for(k=0;k<p;k++)
                printf("%f\t",U0[k*n+i]);
            printf("\n");
        }
        
        printf("\n");
        printf("\n");
        
        for(i=0;i<p;i++){
            for(k=0;k<p;k++)
                printf("%f\t",V0[i*p+k]);
            printf("\n");
        }
        
        printf("\n");
        printf("\n");
         */
        
        if(input==0){
            get_ols_norm(X+j*p*n, y, &n, &p, &ols_norm);
            //printf("%f\t",ols_norm);
            if(ols_norm>lambda_max)
                lambda_max = ols_norm;
        }
        //printf("\n");
            
    }
    free(D0);
    free(WORK);
    free(ols);
    
    int *a_idx,*i_idx;
    int a_size,i_size;
    a_idx = (int *) malloc(d*sizeof(int));
    i_idx = (int *) malloc(d*sizeof(int));
    a_size = 0;
    i_size = d;
    for(j=0;j<d;j++){
        a_idx[j] = -1;
    }
    
    for(j=0;j<d;j++){
        i_idx[j] = 1;
    }
    
    double *r;
    r = (double *) malloc(n*sizeof(double));
    for(i=0;i<n;i++)
        r[i] = y[i];

    for(counter=0;counter<nlambda;counter++)
    {
        if(input==0)
            lambda[counter] = lambda[counter]*lambda_max;
        
        ilambda = lambda[counter];
        
        gap_ext = 1;
        ite_ext = 0;
        while((gap_ext>0) && (ite_ext<max_ite)){
            change_ext = 0;
            for(j=0;j<d;j++){
            
                if(i_idx[j]==1){
                    get_sol_ext(w+j*p,&ilambda,X+j*p*n,r,&n,&p,&w_norm);
                    if(w_norm!=0){
                        i_idx[j] = -1;
                        a_idx[a_size] = j;
                        a_size++;
                        change_ext = 1;
                    }
                }
            }
            if(change_ext==0)
                gap_ext=0;
            else
                gap_ext=1;
        
            /*
        for(j=0;j<d;j++)
            printf("%d\t",a_idx[j]);
        printf("\n");
        
        for(j=0;j<d;j++)
            printf("%d\t",i_idx[j]);
        printf("\n");
        
        printf("%d\n",a_size);
         */
        
        
            gap_int = 1;
            aw_norm = (double *) malloc(a_size*sizeof(double));
            ite_int = 0;
            while((gap_int>thol) && (ite_int<max_ite)){
                gap_tmp0 = 0;
                gap_tmp1 = 0;
                for(j=0;j<a_size;j++){
                    get_sol_int(w+a_idx[j]*p,&ilambda,X+a_idx[j]*p*n,r,&n,&p,aw_norm+j,&gap_tmp0);
                    gap_tmp1 = gap_tmp1 + aw_norm[j];
                }
            
                //printf("%f,%f",gap_tmp0,gap_tmp1);
                //gap_int = 0;
            
                gap_int = gap_tmp0/(gap_tmp1+1e-10);
                ite_int++;
            }
        
            /*
            for(j=0;j<a_size;j++){
                printf("%d\n",a_idx[j]);
                printf("%f\n",aw_norm[j]);
            }
            printf("\n");
            printf("\n");
             */
            
            back = 0;
            for(j=0;j<a_size;j++)
                if(aw_norm[j]==0){
                    back++;
                    i_idx[a_idx[j]] = 1;
                }
                else
                    a_idx[j-back] = a_idx[j];
            for(j=a_size - back;j<a_size;j++)
                a_idx[j] = -1;
            a_size = a_size - back;
        
            /*
            for(j=0;j<d;j++)
                printf("%d\t",a_idx[j]);
            printf("\n");
        
            for(j=0;j<d;j++)
                printf("%d\t",i_idx[j]);
            printf("\n");
        
            printf("%d\n",a_size);
            */
        
            free(aw_norm);
            ite_ext++;
        }
    
        for(j=0;j<a_size;j++){
            get_sol_final(w+a_idx[j]*p,ww+counter*m+a_idx[j]*p,V+a_idx[j]*p*p,&p);
            for(k=0;k<p;k++)
                func_norm[counter*d+a_idx[j]] = func_norm[counter*d+a_idx[j]] + pow(w[a_idx[j]*p+k],2);
            func_norm[counter*d+a_idx[j]] = pow(func_norm[counter*d+a_idx[j]],2);
            
        }
        df[counter] = a_size;
    
        for(i=0;i<n;i++)
            sse[counter] = sse[counter] + pow(r[i],2);
    }
    
    free(i_idx);
    free(a_idx);
    free(r);
    free(V);
    free(w);
}
