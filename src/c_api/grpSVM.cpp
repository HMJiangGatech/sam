#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "R.h"
#include "math.h"
#include "R_ext/BLAS.h"
#include "R_ext/Lapack.h"

void get_residual(double *r,double *Z, double* x, int *xa_idx, int *nn, int *dd, int *pp, int *mm)
{
    int i,j,k,g_idx,b_idx;
    int n,d,p,m;
    n = *nn;
    d = *dd;
    p = *pp;
    m = *mm;
    
    for(i=0;i<n;i++)
        r[i] = 1;

    for(j=0;j<d;j++){
        g_idx = j*p;
        if(xa_idx[j]==1){
            for(k=0;k<p;k++){
                b_idx = g_idx + k;
                for(i=0;i<n;i++)
                    r[i] = r[i] - Z[b_idx*n+i]*x[b_idx];
            }
        }
    }
    
    for(i=0;i<n;i++)
        r[i] = r[i] - Z[m*n+i]*x[m];
}

void get_dual(double *u, double *r, int *ua_idx, double *mmu, int *nn)
{
    int i,n;
    double mu;
    mu = *mmu;
    n = *nn;
    for(i=0;i<n;i++){
        u[i] = r[i]/mu;
        if(u[i]>=1){
            u[i] = 1;
            ua_idx[i] = 1;
        }
        else if(u[i]<=0){
            u[i] = 0;
            ua_idx[i] = 0;
        }
        else ua_idx[i] = 1;
    }
}

void get_grad_SVM(double *g, double *Z, double *u, int *ua_idx, int *mm, int *nn)
{
    int i,j,b_idx;
    int m,n;
    
    m = *mm;
    n = *nn;
    
    for(j=0;j<(m+1);j++){
        b_idx = j*n;
        g[j] = 0;
        for(i=0;i<n;i++){
            if(ua_idx[i]==1)
                g[j] = g[j] - Z[b_idx+i]*u[i];
        }
    }
}

void get_base_SVM(double *H0, double *u, double *r, int *ua_idx, double *mmu, int *nn)
{
    int i,n;
    double mu;
    mu = *mmu;
    n = *nn;
    *H0 = 0;
    for(i=0;i<n;i++){
        if(ua_idx[i]==1)
            *H0 = *H0 + u[i]*r[i] - mu*pow(u[i],2)/2;
        }
}
void grp_sth_SVM(double *sub_x, int *pp, double *iilambda0, double *gnorm)
{
    int k,p;
    double w_ratio,ilambda0;
    p = *pp;
    ilambda0 = *iilambda0;
    
    *gnorm = 0;
    for(k=0;k<p;k++)
        *gnorm = *gnorm + pow(sub_x[k],2);
    *gnorm = sqrt(*gnorm);
    
    if(*gnorm <= ilambda0){
        *gnorm = 0;
        for(k=0;k<p;k++)
            sub_x[k] = 0;
    }
    else {
        w_ratio = (*gnorm - ilambda0)/(*gnorm);
        *gnorm = *gnorm - ilambda0;
        for(k=0;k<p;k++)
            sub_x[k] = sub_x[k]*w_ratio;
    }
}

//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::plugins(openmp)]

extern "C" void grpSVM(double *Z, double *lambda, int *nnlambda, double *LL0, int *nn, int *dd, int *pp, double *aa0, double *xx, double *mmu, int *mmax_ite, double *tthol, double *aalpha, double *df, double *func_norm)
{
    
    int nlambda,n,d,p,m,max_ite;
    double a0,mu,thol,ilambda,L0,alpha;
    double L,ilambda0,gnorm,reg_norm;
    double H0,Q,H,Hx0,Hx1;
    double tmp;
    nlambda = *nnlambda;
    a0 = *aa0;
    n = *nn;
    d = *dd;
    p = *pp;
    m = d*p;
    max_ite = *mmax_ite;
    mu = *mmu;
    thol = *tthol;
    L0 = *LL0;
    alpha = *aalpha;
    
    int ite_ext,ite_int;
    int counter,i,j,k,s,g_idx,b_idx,lambda_idx0,lambda_idx1;
    
    double gap,gap_x,gap_y,tmp_x,gap_xx,gap_yy;
    int iter,tracking;
    //int gap_ext,change_ext;
    //double gap_int;
    
    int *xa_idx,*ya_idx,*ua_idx;
    xa_idx = (int *) malloc(d*sizeof(int));
    ya_idx = (int *) malloc(d*sizeof(int));
    ua_idx = (int *) malloc(n*sizeof(int));
    
    for(j=0;j<d;j++){
        xa_idx[j] = 0;
    }
    
    double *r,*u,*g,*x0,*x1,*y1,*y2,t1,t2;
    r = (double *) malloc(n*sizeof(double));
    u = (double *) malloc(n*sizeof(double));
    g = (double *) malloc((m+1)*sizeof(double));
    x0 = (double *) malloc((m+1)*sizeof(double));
    x1 = (double *) malloc((m+1)*sizeof(double));
    y1 = (double *) malloc((m+1)*sizeof(double));
    y2 = (double *) malloc((m+1)*sizeof(double));
    
    
    for(j=0;j<m;j++){
        x0[j] = 0;
    }
    x0[m] = a0;
    
    
    for(counter=0;counter<nlambda;counter++){
        
        for(j=0;j<(m+1);j++){
            y1[j] = x0[j];
        }
        
        for(j=0;j<d;j++){
            ya_idx[j] = xa_idx[j];
        }
    
        ilambda = lambda[counter];

        
        //printf("%f\n",ilambda);
    
        t1 = 1;
    
        get_residual(r,Z,y1,ya_idx,&n,&d,&p,&m);
        
    
        //for(i=0;i<n;i++)
        //    printf("%f\n",r[i]);
        //printf("\n");
    
    
        get_dual(u,r,ua_idx,&mu,&n);
        
        //for(i=0;i<n;i++)
        //    printf("%f\n",u[i]);
        //printf("\n");
    
    
        get_grad_SVM(g,Z,u,ua_idx,&m,&n);
    
        //for(j=0;j<(m+1);j++)
        //    printf("%f\n",g[j]);
        //printf("\n");
    
    
        L = L0;
    
        get_base_SVM(&H0,u,r,ua_idx,&mu,&n);
    
        //printf("%f\n",H0);
    
        tracking = 1;
    
        while(tracking==1){
        
            ilambda0 = ilambda/L;
            for(j=0;j<(m+1);j++)
                x1[j] = y1[j] - g[j]/L;
        
            for(j=0;j<d;j++){
                grp_sth_SVM(x1+j*p,&p,&ilambda0,&gnorm);
                if(gnorm!=0)
                    xa_idx[j] = 1;
                else
                    xa_idx[j] = 0;
            }
        
            //for(j=0;j<(m+1);j++)
            //    printf("%f\n",x1[j]);
            //printf("\n");
        
            //for(j=0;j<d;j++)
            //    printf("%d,",xa_idx[j]);
            //printf("\n");
        
        
            Q = H0;
            for(j=0;j<(m+1);j++){
                Q = Q + g[j]*(x1[j]-y1[j]) + L*pow(x1[j]-y1[j],2)/2;
            }
        
            //printf("%f\n",Q);
        
            get_residual(r,Z,x1,xa_idx,&n,&d,&p,&m);
        
            //for(i=0;i<n;i++)
            //    printf("%f\n",r[i]);
            //printf("\n");
        
            get_dual(u,r,ua_idx,&mu,&n);
        
            H = 0;
            for(i=0;i<n;i++){
                if(ua_idx[i]==1)
                    H = H + u[i]*r[i] - mu*pow(u[i],2)/2;
            }
        
            //printf("%f\n",H);
        
            if(Q>=H)
                L = L*alpha;
            else{
                L = L/alpha;
                tracking = 0;
            }
        }
    
        //if(iter == 1){
    
        //    printf("%f\n",Q);
        //    printf("%f\n",H);
        //}
    
        ilambda0 = ilambda/L;
        for(j=0;j<(m+1);j++)
            x1[j] = y1[j] - g[j]/L;
    
        reg_norm = 0;
        for(j=0;j<d;j++){
            grp_sth_SVM(x1+j*p,&p,&ilambda0,&gnorm);
            reg_norm = reg_norm + gnorm;
            if(gnorm!=0)
                xa_idx[j] = 1;
            else
                xa_idx[j] = 0;
        }
    
        //for(j=0;j<(m+1);j++)
        //    printf("%f\n",x1[j]);
        //printf("\n");
    
    
        t2 = (1+sqrt(1+4*pow(t1,2)))/2;
        //printf("%f\n",t2);
        tmp = (t1-1)/t2;
        //printf("%f\n",tmp);
    
        for(j=0;j<d;j++){
            g_idx = j*p;
            for(k=0;k<p;k++)
                y2[g_idx+k] = x1[g_idx+k] + tmp*(x1[g_idx+k]-x0[g_idx+k]);
            if(y2[g_idx]!=0)
                ya_idx[j] = 1;
            else
                ya_idx[j] = 0;
        }
        y2[m] = x1[m] + tmp*(x1[m]-x0[m]);
    
        
        gap_x = 0;
        gap_y = 0;
        gap_xx = 0;
        gap_yy = 0;
        for(j=0;j<(m+1);j++){
            tmp_x = x1[j] - x0[j];
            y2[j] = x1[j] + tmp*tmp_x;
            gap_x = gap_x + pow(tmp_x,2);
            gap_y = gap_y + pow(y2[j] - y1[j],2);
            gap_xx = gap_xx + pow(x1[j],2);
            gap_yy = gap_yy + pow(y2[j],2);
            y1[j] = y2[j];
            x0[j] = x1[j];
        }
            gap_x = sqrt(gap_x/(gap_xx+1e-16));
            gap_y = sqrt(gap_y/(gap_yy+1e-16));
        
        if(gap_x>gap_y)
            gap = gap_x;
        else
            gap = gap_y;
     
        //printf("%f\n",gap);
        
        //for(j=0;j<(m+1);j++)
        //    printf("%f\n",y2[j]-x1[j]);
        //printf("\n");
    
    
        //for(j=0;j<(m+1);j++)
        //    printf("%f\n",y2[j]);
        //printf("\n");
    
    
        for(j=0;j<(m+1);j++)
            x0[j] = x1[j];
    
        get_residual(r,Z,x0,xa_idx,&n,&d,&p,&m);
    
        //for(i=0;i<n;i++)
        //    printf("%f\n",r[i]);
        //printf("\n");
    
        get_dual(u,r,ua_idx,&mu,&n);
    
        //for(i=0;i<n;i++)
        //    printf("%f\n",u[i]);
        //printf("\n");
    
    
        Hx0 = ilambda*reg_norm;
        for(i=0;i<n;i++){
            if(ua_idx[i]==1)
                Hx0 = Hx0 + u[i]*r[i] - mu*pow(u[i],2)/2;
        }
    
        //printf("%f\n",Hx0);
    
        for(j=0;j<(m+1);j++)
            y1[j] = y2[j];
    
        //for(j=0;j<(m+1);j++)
        //    printf("%f\n",y1[j]-x1[j]);
        //printf("\n");
    
        t1 = t2;
    
        iter = 0;
    
    
        while((gap>thol)&&(iter<=max_ite)){
        
            get_residual(r,Z,y1,ya_idx,&n,&d,&p,&m);
        
            //for(i=0;i<n;i++)
            //    printf("%f\n",r[i]);
            //printf("\n");
        
        
            get_dual(u,r,ua_idx,&mu,&n);
            
            //for(i=0;i<n;i++)
            //    printf("%f\n",u[i]);
            //printf("\n");
        
        
            get_grad_SVM(g,Z,u,ua_idx,&m,&n);
        
            get_base_SVM(&H0,u,r,ua_idx,&mu,&n);
        
            if(L<L0){
        
                tracking = 1;
        
                while(tracking==1){
            
                    ilambda0 = ilambda/L;
                    for(j=0;j<(m+1);j++)
                        x1[j] = y1[j] - g[j]/L;
            
                    reg_norm = 0;
                    for(j=0;j<d;j++){
                        grp_sth_SVM(x1+j*p,&p,&ilambda0,&gnorm);
                        reg_norm = reg_norm + gnorm;
                        if(gnorm!=0)
                            xa_idx[j] = 1;
                        else
                            xa_idx[j] = 0;
                    }
            
                    //for(j=0;j<(m+1);j++)
                    //    printf("%f\n",x1[j]);
                    //printf("\n");
            
                    //for(j=0;j<d;j++)
                    //    printf("%d,",xa_idx[j]);
                    //printf("\n");
            
            
                    Q = H0;
                    for(j=0;j<(m+1);j++){
                        Q = Q + g[j]*(x1[j]-y1[j]) + L*pow(x1[j]-y1[j],2)/2;
                    }
            
                    //printf("%f\n",Q);
            
                    get_residual(r,Z,x1,xa_idx,&n,&d,&p,&m);
                    
                    //for(i=0;i<n;i++)
                    //    printf("%f\n",r[i]);
                    //printf("\n");
                    
                    
                    get_dual(u,r,ua_idx,&mu,&n);
                                
                    H = 0;
                    for(i=0;i<n;i++){
                        if(ua_idx[i]==1)
                            H = H + u[i]*r[i] - mu*pow(u[i],2)/2;
                    }
            
                    if(Q>H)
                        tracking = 0;
                    else
                        L = L/alpha;
                }
            }
            else{
                ilambda0 = ilambda/L0;
                for(j=0;j<(m+1);j++)
                    x1[j] = y1[j] - g[j]/L0;
            
                reg_norm = 0;
                for(j=0;j<d;j++){
                    grp_sth_SVM(x1+j*p,&p,&ilambda0,&gnorm);
                    reg_norm = reg_norm + gnorm;
                    if(gnorm!=0)
                        xa_idx[j] = 1;
                    else
                        xa_idx[j] = 0;
                }
            
                get_residual(r,Z,x1,xa_idx,&n,&d,&p,&m);
                
                //for(i=0;i<n;i++)
                //    printf("%f\n",r[i]);
                //printf("\n");
                
                
                get_dual(u,r,ua_idx,&mu,&n);
                
                //for(i=0;i<n;i++)
                //    printf("%f\n",u[i]);
                //printf("\n");
                
                H = 0;
                for(i=0;i<n;i++){
                    if(ua_idx[i]==1)
                        H = H + u[i]*r[i] - mu*pow(u[i],2)/2;
                }
            }
        
            //printf("%f\n",Q);
            //printf("%f\n",H);
        
        
            Hx1 = H + ilambda*reg_norm;
        
            //printf("%f\n",Hx1);
        
            t2 = (1+sqrt(1+4*pow(t1,2)))/2;
            tmp = (t1-1)/t2;
        
        
            if(Hx1<Hx0){
                Hx0 = Hx1;
                gap_x = 0;
                gap_y = 0;
                gap_xx = 0;
                gap_yy = 0;
                for(j=0;j<(m+1);j++){
                    tmp_x = x1[j] - x0[j];
                    y2[j] = x1[j] + tmp*tmp_x;
                    gap_x = gap_x + pow(tmp_x,2);
                    gap_y = gap_y + pow(y2[j] - y1[j],2);
                    gap_xx = gap_xx + pow(x1[j],2);
                    gap_yy = gap_yy + pow(y2[j],2);
                    y1[j] = y2[j];
                    x0[j] = x1[j];
                }
                gap_x = sqrt(gap_x/(gap_xx+1e-16));
                gap_y = sqrt(gap_y/(gap_yy+1e-16));
            
            
                for(j=0;j<d;j++){
                    g_idx = j*p;
                    for(k=0;k<p;k++)
                        if(y1[g_idx]!=0)
                            ya_idx[j] = 1;
                        else
                            ya_idx[j] = 0;
                }
            
                if(gap_x<gap_y)
                    gap = gap_y;
                else
                    gap = gap_x;
        
                //printf("%f\n",gap_x);
                //printf("%f\n",gap_y);
            }
            else{
                gap = 0;
                gap_yy = 0;
                for(j=0;j<(m+1);j++){
                    gap = gap + pow(y1[j] - x0[j],2);
                    gap_yy = gap_yy + pow(x0[j],2);
                    y1[j] = x0[j];
                }
                for(j=0;j<d;j++){
                    g_idx = j*p;
                    for(k=0;k<p;k++)
                        if(y1[g_idx]!=0)
                            ya_idx[j] = 1;
                        else
                            ya_idx[j] = 0;
                }
                gap = sqrt(gap/gap_yy);
                //printf("%f\n",gap);
            }
        
            t1 = t2;
        
            //printf("%f\n",gap);
        
            //if(iter<10){
            //  printf("%f\n",Hx0);
            //}
        
       
            //for(j=0;j<(m+1);j++)
            //  printf("%f\n",x0[j]);
            //printf("\n");
        
            //for(j=0;j<(m+1);j++)
            //    printf("%f\n",y1[j]);
            //printf("\n");
        
            iter = iter + 1;
        }
    
        //printf("%d\n",iter);
    
    
        //df[counter] = 0;
        lambda_idx0 = counter*(m+1);
        lambda_idx1 = counter*d;
        df[counter] = 0;
        for(j=0;j<d;j++){
            func_norm[lambda_idx1+j] = 0;
            if(xa_idx[j]==1){
                df[counter]++;
                g_idx = j*p;
                for(k=0;k<p;k++){
                    xx[lambda_idx0 + g_idx+k] = x0[g_idx+k];
                    func_norm[lambda_idx1+j] = func_norm[lambda_idx1+j] + pow(x0[g_idx+k],2);
                }
                func_norm[lambda_idx1+j] = sqrt(func_norm[lambda_idx1+j]);
            }
        }
        xx[lambda_idx0+m] = x0[m];
    
    
        //printf("%d\n",iter);
    
        /*
        for(j=0;j<(m+1);j++)
            printf("%f\n",x0[j]);
        printf("\n");
    
        for(j=0;j<(m+1);j++)
            printf("%f\n",y1[j]);
        printf("\n");
    
        for(j=0;j<d;j++)
            printf("%d\n",xa_idx[j]);
    
        printf("\n");
        for(j=0;j<d;j++)
            printf("%d\n",ya_idx[j]);
    
        printf("\n");
    
        printf("%f\n",gap);
         */
    
    }
    
    free(xa_idx);
    free(ya_idx);
    free(ua_idx);
    free(r);
    free(u);
    free(g);
    free(x0);
    free(x1);
    free(y1);
    free(y2);
}
