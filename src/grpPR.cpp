#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "R.h"
#include "math.h"
#include "R_ext/BLAS.h"
#include "R_ext/Lapack.h"

void grp_sth_PR(double *sub_x, int *pp, double *iilambda0, double *gnorm)
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

void get_EC(double *ec, double *A, double *x, int *xa_idx, int *nn, int *dd, int* pp, int *mm)
{
    int i,j,k,b_idx,g_idx;
    int n,d,p,m;
    n = *nn;
    d = *dd;
    p = *pp;
    m = *mm;
    
    for(i=0;i<n;i++)
        ec[i] = x[m];
    
    for(j=0;j<d;j++)
        if(xa_idx[j]==1){
            g_idx = j*p;
            for(k=0;k<p;k++){
                b_idx = n*(g_idx+k);
                for(i=0;i<n;i++)
                    ec[i] = ec[i] + A[b_idx+i]*x[g_idx+k];
            }
        }
    for(i=0;i<n;i++){
        ec[i] = exp(ec[i]);
    }
}

void get_grad_PR(double *g, double *A, double *ec, double *z, int *nn, int *mm)
{
    int i,j,b_idx;
    int n,m;
    n = *nn;
    m = *mm;
    
    for(j=0;j<(m+1);j++){
        g[j] = -z[j];
        b_idx = n*j;
        for(i=0;i<n;i++)
            g[j] = g[j] + ec[i]*A[b_idx+i];
    }
}

void get_base_PR(double *H0, double *z, double *x, double *ec, int *xa_idx, int *nn, int *dd, int *pp, int *mm)
{
    int i,j,k,g_idx;
    int n,m,d,p;
    n = *nn;
    m = *mm;
    d = *dd;
    p = *pp;
    
    *H0 = -z[m]*x[m];
    for(j=0;j<d;j++){
        if(xa_idx[j]==1){
            g_idx = j*p;
            for(k=0;k<p;k++)
                *H0 = *H0 - z[g_idx+k]*x[g_idx+k];
        }
    }

    for(i=0;i<n;i++)
        *H0 = *H0 + ec[i];
}

void grpPR(double *A, double *lambda, int *nnlambda, double *LL0, int *nn, int *dd, int *pp, double *xx, double *aa0, int *mmax_ite, double *tthol, double *aalpha, double *z, int *df, double *func_norm)
{
    
    int n,d,p,m,max_ite;
    double thol,ilambda,L0,alpha,nlambda;
    double L,ilambda0,gnorm,reg_norm;
    double H0,Q,H,Hx0,Hx1;
    double tmp;
    n = *nn;
    d = *dd;
    p = *pp;
    m = d*p;
    max_ite = *mmax_ite;
    thol = *tthol;
    L0 = *LL0;
    alpha = *aalpha;
    nlambda = *nnlambda;
    
    int ite_ext,ite_int;
    int i,j,k,s,g_idx,b_idx,lambda_idx0,lambda_idx1;
    
    double gap,gap_x,gap_y,tmp_x,gap_xx,gap_yy;
    int counter,iter,tracking;
    //int gap_ext,change_ext;
    //double gap_int;
    
    int *xa_idx,*ya_idx;
    xa_idx = (int *) malloc(d*sizeof(int));
    ya_idx = (int *) malloc(d*sizeof(int));
    
    double *ec,*g,*x0,*x1,*y1,*y2,t1,t2;
    ec = (double *) malloc(n*sizeof(double));
    g = (double *) malloc((m+1)*sizeof(double));
    x0 = (double *) malloc((m+1)*sizeof(double));
    x1 = (double *) malloc((m+1)*sizeof(double));
    y1 = (double *) malloc((m+1)*sizeof(double));
    y2 = (double *) malloc((m+1)*sizeof(double));
    
    for(j=0;j<d;j++){
        xa_idx[j] = 0;
    }
    
    for(j=0;j<m;j++){
        x0[j] = 0;
    }
    x0[m] = *aa0;
    
    for(counter=0;counter<nlambda;counter++){
    
        ilambda = lambda[counter];
        
        L = L0;
        
        for(j=0;j<(m+1);j++)
            y1[j] = x0[j];
    
        for(j=0;j<d;j++){
            ya_idx[j] = xa_idx[j];
        }

        t1 = 1;
    
        get_EC(ec,A,y1,ya_idx,&n,&d,&p,&m);
    
        /*
        for(i=0;i<n;i++)
            printf("%f\n",e[i]);
        printf("\n");
     
        for(i=0;i<n;i++)
            printf("%f\n",ep[i]);
        printf("\n");
        */

    
        get_grad_PR(g,A,ec,z,&n,&m);
    
        /*
    for(j=0;j<(m+1);j++)
        printf("%f\n",g[j]);
    printf("\n");
    */
    
        get_base_PR(&H0, z, y1, ec, ya_idx, &n, &d, &p, &m);

        //printf("%f\n",H0);
        //printf("\n");
    
        tracking = 1;
    
        while(tracking==1){
        
            ilambda0 = ilambda/L;
            for(j=0;j<(m+1);j++)
                x1[j] = y1[j] - g[j]/L;
        
            for(j=0;j<d;j++){
                grp_sth_PR(x1+j*p,&p,&ilambda0,&gnorm);
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
        
            get_EC(ec,A,x1,xa_idx,&n,&d,&p,&m);
        
            get_base_PR(&H, z, x1, ec, xa_idx, &n, &d, &p, &m);

            //printf("%f\n",H);
        
            if(Q>=H)
                L = L*alpha;
            else{
                L = L/alpha;
                tracking = 0;
            }
        }
    
        //printf("%f\n",Q);
        //printf("%f\n",H);
    
        ilambda0 = ilambda/L;
        for(j=0;j<(m+1);j++)
            x1[j] = y1[j] - g[j]/L;
    
        reg_norm = 0;
        for(j=0;j<d;j++){
            grp_sth_PR(x1+j*p,&p,&ilambda0,&gnorm);
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
        gap_x = sqrt(gap_x/(gap_xx+1e-10));
        gap_y = sqrt(gap_y/(gap_yy+1e-10));
    
        if(gap_x>gap_y)
            gap = gap_x;
        else
            gap = gap_y;

        //printf("%f\n",gap);
        //printf("%f\n",L);
    
        //for(j=0;j<(m+1);j++)
        //    printf("%f\n",y2[j]-x1[j]);
        //printf("\n");
    
    
        //for(j=0;j<(m+1);j++)
        //    printf("%f\n",y2[j]);
        //printf("\n");
    
        for(j=0;j<(m+1);j++)
            x0[j] = x1[j];
    
        get_EC(ec,A,x1,xa_idx,&n,&d,&p,&m);
    
        get_base_PR(&Hx0, z, x1, ec, xa_idx, &n, &d, &p, &m);

        Hx0 = Hx0 + ilambda*reg_norm;
    
        //printf("%f\n",Hx0);
    
        for(j=0;j<(m+1);j++)
            y1[j] = y2[j];
    
        //for(j=0;j<(m+1);j++)
        //    printf("%f\n",y1[j]-x1[j]);
        //printf("\n");
    
        t1 = t2;
    
        iter = 0;
    
        while((gap>thol)&&(iter<=max_ite)){
        
            get_EC(ec,A,y1,ya_idx,&n,&d,&p,&m);
        
            /*
         for(i=0;i<n;i++)
         printf("%f\n",e[i]);
         printf("\n");
         
         for(i=0;i<n;i++)
         printf("%f\n",ep[i]);
         printf("\n");
         */
        
        
            get_grad_PR(g,A,ec,z,&n,&m);
        
            /*
         for(j=0;j<(m+1);j++)
         printf("%f\n",g[j]);
         printf("\n");
         */
        
            get_base_PR(&H0, z, y1, ec, ya_idx, &n, &d, &p, &m);
        
            tracking = 1;
            
            while(tracking==1){
                
                ilambda0 = ilambda/L;
                for(j=0;j<(m+1);j++)
                    x1[j] = y1[j] - g[j]/L;
            
                reg_norm = 0;
                for(j=0;j<d;j++){
                    grp_sth_PR(x1+j*p,&p,&ilambda0,&gnorm);
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
            
                get_EC(ec,A,x1,xa_idx,&n,&d,&p,&m);
            
                get_base_PR(&H, z, x1, ec, xa_idx, &n, &d, &p, &m);
                
                if(Q>H)
                    tracking = 0;
                else
                    L = L/alpha;
            }
        
            //printf("%f\n",Q);
            //printf("%f\n",H);
            //printf("%f\n",L);
        
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
                gap_x = sqrt(gap_x/(gap_xx+1e-10));
                gap_y = sqrt(gap_y/(gap_yy+1e-10));
            
            
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
        
        
            //gap = 0;
        
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
    free(ec);
    free(g);
    free(x1);
    free(y1);
    free(y2);
}
