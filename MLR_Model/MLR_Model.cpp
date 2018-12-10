// MLR_Model.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"

#include <stdlib.h>



double	 *	 Mlr(int n, double   x[], int xm, double y[], int ym);                     //多元线性回归函数声明

double   *   MatrixOpp(double   *A,int   m,int   n);   /*矩阵求逆*/ //m 维数，n 样本数
double   *   MatrixInver(double   *A,int   m,int   n);   /*矩阵转置*/ 
double   Surplus(double   A[],int   m,int   n);   /*求矩阵行列式*/ 
double	 *	 MatrixMul(int n, double   A[], int xm, double B[], int ym);//矩阵相乘


int _tmain(int argc, _TCHAR* argv[])
{


	int   i,j; 
	double	 *b;
	double   x[3][4]={6,2,3,4,2,6,7,8,9,3,4,9};
	double   y[1][4]={1,2,3,4};


	b = Mlr(4,(double   *)x,3,(double   *)y,1);  // b=x*y'/(x*x')

	printf("MLR回归系数为：\n");
	for(i=0;   i <3;   i++) 
		{
			for(j=0;   j <1;   j++) 
				printf("%lf  ",b[i*1+j]);
			printf("\n");
		};


	free(b);

	return 0;
}



double   *   MatrixInver(double   A[],int   m,int   n)   /*矩阵转置*/ 
{ 
          int   i,j; 
          double   *B=NULL; 
          B=(double   *)malloc(m*n*sizeof(double)); 
        
          for(i=0;i <n;i++) 
          for(j=0;j <m;j++) 
          B[i*m+j]=A[j*n+i];         
          return   B; 
} 

double   Surplus(double   A[],int   m,int   n)   /*求矩阵行列式*/ 
{         
          int   i,j,k,p,r; 
          double   X,temp=1,temp1=1,s=0,s1=0; 
        
          if(n==2) 
          { 
                  for(i=0;i <m;i++) 
                  for(j=0;j <n;j++) 
                          if((i+j)%2)   temp1*=A[i*n+j]; 
                          else   temp*=A[i*n+j]; 
                  X=temp-temp1; 
          } 
          else 
          { 
                  for(k=0;k <n;k++) 
                  { 
                          for(i=0,j=k;i <m,j <n;i++,j++) 
                          temp*=A[i*n+j]; 
                          if(m-i) 
                          { 
                                  for(p=m-i,r=m-1;p> 0;p--,r--) 
                                  temp*=A[r*n+p-1]; 
                          } 
                          s+=temp; 
                          temp=1; 
                  } 
                
                  for(k=n-1;k >= 0;k--) 
                  { 
                          for(i=0,j=k;i <m,j >=0;i++,j--) 
                          temp1*=A[i*n+j]; 
                          if(m-i) 
                          {for(p=m-1,r=i;r <m;p--,r++) 
                          temp1*=A[r*n+p];} 
                          s1+=temp1; 
                          temp1=1; 
                  }                 
                  X=s-s1; 
          } 
          return   X; 
} 

double   *   MatrixOpp(double   A[],int   m,int   n)   /*矩阵求逆*/ 
{ 
          int   i,j,x,y,k; 
          double   *SP=NULL,*AB=NULL,*B=NULL,X,*C=NULL; 
          SP=(double   *)malloc(m*n*sizeof(double)); 
          AB=(double   *)malloc(m*n*sizeof(double)); 
          B=(double   *)malloc(m*n*sizeof(double)); 
		  C=(double   *)malloc(m*m*sizeof(double)); 
          X=Surplus(A,m,n); 
          X=1/X; 
        
          for(i=0;i <m;i++) 
          for(j=0;j <n;j++) 
          { 
                  for(k=0;k <m*n;k++) 
                  B[k]=A[k]; 
                  { 
                          for(x=0;x <n;x++) 
                          B[i*n+x]=0; 
                          for(y=0;y <m;y++) 
                          B[m*y+j]=0; 
                          B[i*n+j]=1; 
                          SP[i*n+j]=Surplus(B,m,n); 
                          AB[i*n+j]=X*SP[i*n+j]; 
                  } 
          } 
          C=MatrixInver(AB,m,n); 
        free(B);
		free(SP);
		free(AB);
          return   C; 
} 

double	 *	 MatrixMul(int n, double   A[], int xm, double B[], int ym)  //矩阵相乘
{
	int	 i,j,k;
	double   sum;
	double	 * C = NULL;
	
	C=(double *)malloc(xm*ym*sizeof(double));

	for(i=0; i<xm; ++i)
		for(j=0; j<ym; ++j)
		{
			sum=0;
			for(k=0; k<n; ++k)
			{
				sum=sum+A[i*n+k]*B[j*n+k];
			}
			C[i*ym+j]=sum;
		};
	return C;
}

double   *   Mlr(int n, double x[], int xm,double y[], int ym )      // b=x*y'/(x*x')
{
	double  *xx_arr;
	double  *xx_inv; 
	double  *xy_arr;
	double	 * result;

	xy_arr = MatrixMul(n,(double   *)x,xm,(double   *)y,ym);//矩阵相乘
	xx_arr = MatrixMul(n,(double   *)x,xm,(double   *)x,xm);//

	xx_inv = MatrixOpp((double   *)xx_arr,xm,xm);     //求逆 

	result = MatrixMul(xm,(double   *)xy_arr,ym,(double   *)xx_inv,xm); 


	free(xy_arr);
	free(xx_arr);
	free(xx_inv);

	return result;
}
