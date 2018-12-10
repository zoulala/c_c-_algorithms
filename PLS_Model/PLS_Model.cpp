// PLS_Model.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"

#include <stdlib.h>
#include <math.h>


/*----------------- 函数声明部分 --------------------*/


																					/*m 维数，n 样本数*/
double   Surplus(double   A[],int   m,int   n);												//求矩阵行列式
int	 *	 EigSort(double	  A[],int	m);														//特征值排序,返回序号数组
double	 *	 MatrixZscore(double	*A,int m, int n);										//矩阵标准化,并返回每行：均值+标准差；
double	 *	 MatrixMul(int n, double   A[], int xm, double B[], int ym);					//矩阵相乘
double   *   MatrixInver(double   *A,int   m,int   n);										//矩阵转置
double   *   MatrixOpp(double   *A,int  m,int   n); 										//矩阵求逆
void	 MatrixEig(double	*A,double	*B,double	*C,int	m,double eps,int	nit);		//矩阵特征值&特征向量,B存储特征值，C存储特征向量，eps收敛精度，nit迭代次数
double	 *	 MatrixSub(double   A[], double B[], int m, int n);								//矩阵相减
int eejcb(double a[],int n,double v[],double eps,int jt);
double   *   Mlr(int n, double x[], int xm,double y[], int ym );							//最小二乘，b=x*y'/(x*x')


/*----------------- 主函数入口 ------------------*/
int _tmain(int argc, _TCHAR* argv[])
{
	//int kk,p;
	//double eps=0.0001;
	//double a[]={10,2,5,2,4,3,5,3,-1};
	//double vlue[3];
	//double vector[9];
	//double ve[9];
	//MatrixEig(a, vlue, vector, 3, 0.0001, 20);
	//kk= eejcb(a,3,ve,eps,20);
	//for(p=0; p<9; p++)
	//	printf("%lf ",ve[p]);

	int i,j;
	double X[]={1,2,3,9,0,7,5,2,3,2.2,1,5.3};
	double Y[]={3,2,1,3.1};
	int m=3;
	int n=4;
	int r_pls=2;    //选取主元个数
	double *B;	//存放主元回归系数
	double *T=NULL;
	T=(double*)malloc(m*n*sizeof(double));//存放得分矩阵

	double * X_std_mean;
	double * Y_std_mean;
	X_std_mean = MatrixZscore(X,m, n);
	Y_std_mean = MatrixZscore(Y,1, n);

	double * E0=NULL;
	double * F0=NULL;
	E0=(double*)malloc(m*n*sizeof(double));
	F0=(double*)malloc(1*n*sizeof(double));
	for(i=0; i<n*m; i++)
		E0[i] = X[i];
	for(i=0; i<n*1; i++)
		F0[i] = Y[i];


	
	double	w_all[50];		//存放w特征向量矩阵
	double	ap_all[50];	//得分回归系数


	for(i=0; i<m; i++)
	{
		double *temp1,*temp2,*M;
		temp1 = MatrixMul(n,E0,m,F0,1);
		temp2 = MatrixMul(1,temp1,m,F0,n);
		M	  = MatrixMul(n,temp2,m,E0,m);//M=E0'*F0*F0'*E0;


		double eigvlue[3],eigvector[9];
		double eps=0.00001;
		int nit=20;
		MatrixEig(M, eigvlue, eigvector, m, eps, nit);


		double w[10];
		double *t;
		for(j=0; j<m; j++)
		{
			w[j]=eigvector[j];//w为最大特征值对应的特征向量
		}
		for(j=0; j<m; j++)
			w_all[i*m+j]=w[j];

		double *E0_tran;
		E0_tran = MatrixInver(E0, m, n);
		t	    = MatrixMul(m,E0_tran,n,w,1);//计算主成分
		for(j=0; j<n; j++)
			T[i*n+j]=t[j];


		double *t2;
		double *ap;
		t2	= MatrixMul(n,t,1,t,1);
		ap  = MatrixMul(n,E0,m,t,1);
		for(j=0; j<m; j++)
		{
			ap[j] = ap[j]/t2[0];//计算特征矩阵
		}
		for(j=0; j<m; j++)
			ap_all[i*m+j]=ap[j];

		double *Et,*Etran,*E1;
		Et	  = MatrixMul(1,t,n,ap,m);
		Etran = MatrixInver(Et,n,m);
		E1    = MatrixSub(E0,Etran,m,n);//计算残差矩阵 

		for(j=0; j<m*n; j++)
		{
			E0[j]=E1[j];
		}
	
	}

	

	for(i=0; i<m; i++)
	{
		for(j=0; j<n; j++)
			printf("%lf  ",T[i*n+j]);
		printf("\n");
	}


	double *TT;

	B=Mlr(n, (double*)T, m,(double*)Y,  1 );      // b=x*y'/(x*x')

	printf("主元回归系数：\n");
	for(i=0;i<m;i++)
	{
		printf("%lf  ",B[i]);
	}


	//计算降维矩阵p[]
	double eye[50];
	double temp[50],temp3[50],temp4[50],temp5[50],temp6[10];//变量维数m不能超过7
	double p[50];
	double *temp7;
	for(i=0; i<m*m; i++)
	{
		temp[i] = 0.0;
	}
	for(i=0; i<m; i++)
	{
		temp[i*m+i] = 1.0;
	}
	for(i=0; i<m*m; i++)
		eye[i]=temp[i];

	for(i=0; i<m; i++)
	{
		p[i]=w_all[i];
	}
	for(i=1; i<m; i++)
	{
		for(j=0; j<i; j++)
		{
			int r,k,l;
			for(r=0; r<m; r++)
				for(k=0; k<m; k++)
					temp3[r*m+k]=w_all[j*m+r]*ap_all[j*m+k];
			for(r=0; r<m*m; r++)
				temp4[r]=eye[r]-temp3[r];

			for(r=0; r<m; r++)
			{
				for(k=0; k<m; k++)
				{
					double sum=0;
					for(l=0; l<m; l++)
						sum=sum+temp[r*m+l]*temp4[l*m+k];
					temp5[r*m+k]=sum;
				}
			}

			for(r=0; r<m*m; r++)
				temp[r]=temp5[r];
		}

	
		for(j=0; j<m; j++)
			temp6[j]=w_all[i*m+j];
		temp7=MatrixMul(m, temp, m, temp6, 1);
		for(j=0; j<m; j++)
			p[i*m+j]=temp7[j];
	}

	return 0;
}


/*----------------- 函数定义部分 -------------------*/
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

double	 *	 MatrixMul(int n, double   A[], int xm, double B[], int ym)  //矩阵相乘,返回xm*ym的矩阵
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

int	 *	 EigSort(double	A[], int  m)
{
	double temp;
	int i,j,B_temp;
	int *B = NULL;
	B=(int   *)malloc(m*sizeof(int)); 
	
	B[0]=0;
	for(i=1; i<m; ++i)
	{
		B[i]=B[i-1]+1;	//B中存储序号
	}

	for(i=0; i<m; ++i)
	{
		for(j=i+1; j<m; ++j)
		{
			if(A[i]<A[j])
			{
				temp = A[i];
				A[i] = A[j];
				A[j] = temp;

				B_temp=B[i];
				B[i] = B[j];
				B[j] = B_temp;
			}
		}
	}

	return B;
}

double	 *	 MatrixZscore(double	A[],int m, int n)						//将A矩阵标准化，并返回每行：均值+标准差；
{
	int i,j;
	double sum,sum2;
	double mean[10],var[10];
	double *B=NULL;
	B=(double*)malloc(2*m*sizeof(double));
	
	for(i=0; i<m; i++)
	{
		sum=0.0;sum2=0.0;
		for(j=0; j<n; j++)
		{
			sum=sum+A[i*n+j];
		}
		mean[i] = sum/n;
		for(j=0; j<n; j++)
		{
			sum2=sum2+(A[i*n+j]-mean[i])*(A[i*n+j]-mean[i]);
		}
		var[i]=sum2/(n-1);
		for(j=0; j<n; j++)
		{
			A[i*n+j] = (A[i*n+j]-mean[i])/sqrt(var[i]);
		}
	
		B[i*2]=mean[i];
		B[i*2+1]=sqrt(var[i]);

	}

	return B;
}

void	MatrixEig(double	*A,double	*B,double	*C,int	m,double eps,int	nit)		//矩阵特征值&特征向量,B存储特征值，C存储特征向量，eps收敛精度，nit迭代次数.(雅格比(Jacobi)方法求实对称矩阵的全部特征值及特征向量)
{
	/*--初始化特征向量为单位对角矩阵--*/
	int i,j;
	for(i=0; i<m*m; i++)
	{
		C[i] = 0.0;
	}
	for(i=0; i<m; i++)
	{
		C[i*m+i] = 1.0;
	}


	/*--迭代开始--*/
	int nCount = 0;
	while(1)
	{
		/*--找到实对称矩阵A的非对角线上元素最大值--*/
		double nMax = A[1];
		int nRow = 0;
		int nCol = 1;

		for ( i = 0; i < m; i ++)			//行
		{
			for ( j = 0; j < m; j ++)		//列
			{
				double d = fabs(A[i*m+j]); 

				if((i!=j) && (d> nMax)) 
				{ 
					nMax = d;   
					nRow = i;   
					nCol = j; 
				} 
			}
		}

		if(nMax < eps)     //精度符合要求 
			break;  

		if(nCount > nit)       //迭代次数超过限制
			break;

		nCount++;
		

		double dbApp = A[nRow*m+nRow];
		double dbApq = A[nRow*m+nCol];
		double dbAqq = A[nCol*m+nCol];

		//计算旋转角度
		double dbAngle = 0.5*atan2(-2*dbApq,dbAqq-dbApp);
		double dbSinTheta = sin(dbAngle);
		double dbCosTheta = cos(dbAngle);
		double dbSin2Theta = sin(2*dbAngle);
		double dbCos2Theta = cos(2*dbAngle);

		A[nRow*m+nRow] = dbApp*dbCosTheta*dbCosTheta + 
			dbAqq*dbSinTheta*dbSinTheta + 2*dbApq*dbCosTheta*dbSinTheta;
		A[nCol*m+nCol] = dbApp*dbSinTheta*dbSinTheta + 
			dbAqq*dbCosTheta*dbCosTheta - 2*dbApq*dbCosTheta*dbSinTheta;
		A[nRow*m+nCol] = 0.5*(dbAqq-dbApp)*dbSin2Theta + dbApq*dbCos2Theta;
		A[nCol*m+nRow] = A[nRow*m+nCol];

		for(int i = 0; i < m; i ++) 
		{ 
			if((i!=nCol) && (i!=nRow)) 
			{ 
				int u = i*m + nRow;	//p  
				int w = i*m + nCol;	//q
				nMax = A[u]; 
				A[u]= A[w]*dbSinTheta + nMax*dbCosTheta; 
				A[w]= A[w]*dbCosTheta - nMax*dbSinTheta; 
			} 
		} 

		for (int j = 0; j < m; j ++)
		{
			if((j!=nCol) && (j!=nRow)) 
			{ 
				int u = nRow*m + j;	//p
				int w = nCol*m + j;	//q
				nMax = A[u]; 
				A[u]= A[w]*dbSinTheta + nMax*dbCosTheta; 
				A[w]= A[w]*dbCosTheta - nMax*dbSinTheta; 
			} 
		}

		//计算特征向量
		for(int i = 0; i < m; i ++) 
		{ 
			int u = i*m + nRow;		//p   
			int w = i*m + nCol;		//q
			nMax = C[u]; 
			C[u] = C[w]*dbSinTheta + nMax*dbCosTheta; 
			C[w] = C[w]*dbCosTheta - nMax*dbSinTheta; 
		} 
	}


	//特征值,即A主对角线上的元素 
	for(i=0; i<m; i++)
	{
		B[i]=A[i*m+i];
	}

	/*---特征值&特征向量 按大小排序---*/
	int * num;
	num = EigSort(B,m);//大到小排序，num存序号

	double *C_Sort;
	C_Sort = (double*)malloc(m*m*sizeof(double));
	for(i=0; i<m; i++)
	{
		for(j=0; j<m; j++)
		{
			C_Sort[i*m+j]=C[j*m+num[i]];//特征向量按特征值大小顺序排序，且每行是一个特征向量；
		}
	}
	for(i=0; i<m*m; i++)
	{
		C[i] = C_Sort[i];
	}

	free(C_Sort);
	return;
}

double	 *	 MatrixSub(double   A[], double B[], int m, int n)								//矩阵相减
{
	int	 i;
	double	 * C = NULL;
	
	C=(double *)malloc(m*n*sizeof(double));

	for(i=0; i<m*n; ++i)
	{
		C[i]=A[i]-B[i];
	}
	return C;
}

//求特征值特征向量的雅格比(Jacobi)方法2.
int eejcb(double a[],int n,double v[],double eps,int jt) 
//求实对称矩阵的特征值及特征向量的雅格比法 
//利用雅格比(Jacobi)方法求实对称矩阵的全部特征值及特征向量 
//返回值小于0表示超过迭代jt次仍未达到精度要求 
//返回值大于0表示正常返回 
//a-长度为n*n的数组，存放实对称矩阵，返回时对角线存放n个特征值 
//n-矩阵的阶数 
//u-长度为n*n的数组，返回特征向量(按列存储) 
//eps-控制精度要求 
//jt-整型变量，控制最大迭代次数 
{ 
	int i,j,p,q,u,w,t,s,l; 
	double fm,cn,sn,omega,x,y,d; 
	l=1; 

	/*--初始化特征向量为单位对角矩阵--*/
	for (i=0; i<=n-1; i++) 
	{ 
		v[i*n+i]=1.0; 
		for (j=0; j<=n-1; j++) 
		{ 
			if (i!=j) 
			{ 
				v[i*n+j]=0.0; 
			} 
		} 
	} 

	/*--迭代开始--*/
	while (1==1) 
	{ 
		/*--找到实对称矩阵A的非对角线上元素最大值--*/
		fm=0.0; 
		for (i=0; i<=n-1; i++) 
		{ 
			for (j=0; j<=n-1; j++) 
			{ 
				d=fabs(a[i*n+j]); 
				if ((i!=j)&&(d>fm)) 
				{ 
					fm=d; 
					p=i; 
					q=j; 
				} 
			} 
		} 

		/*--迭代结束条件--*/
		if (fm<eps) 
		{ 
			return(1); 
		} 
		if (l>jt) 
		{ 
			return(-1); 
		} 

		/*--计算旋转角度--*/		
		l=l+1; 
		u=p*n+q; 
		w=p*n+p; 
		t=q*n+p; 
		s=q*n+q; 
		x=-a[u]; 
		y=(a[s]-a[w])/2.0; 
		omega=x/sqrt(x*x+y*y); 
		if (y<0.0) 
		{ 
			omega=-omega; 
		} 
		sn=1.0+sqrt(1.0-omega*omega); 
		sn=omega/sqrt(2.0*sn); 
		cn=sqrt(1.0-sn*sn); 
		fm=a[w]; 
		a[w]=fm*cn*cn+a[s]*sn*sn+a[u]*omega; 
		a[s]=fm*sn*sn+a[s]*cn*cn-a[u]*omega; 
		a[u]=0.0; 
		a[t]=0.0; 

		for (j=0; j<=n-1; j++) 
		{ 
			if ((j!=p)&&(j!=q)) 
			{ 
				u=p*n+j; 
				w=q*n+j; 
				fm=a[u]; 
				a[u]=fm*cn+a[w]*sn; 
				a[w]=-fm*sn+a[w]*cn; 
			} 
		} 

		for (i=0; i<=n-1; i++) 
		{ 
			if ((i!=p)&&(i!=q)) 
			{ 
				u=i*n+p; 
				w=i*n+q; 
				fm=a[u]; 
				a[u]=fm*cn+a[w]*sn; 
				a[w]=-fm*sn+a[w]*cn; 
			} 
		} 

		//计算特征向量
		for (i=0; i<=n-1; i++) 
		{ 
			u=i*n+p; 
			w=i*n+q; 
			fm=v[u]; 
			v[u]=fm*cn+v[w]*sn; 
			v[w]=-fm*sn+v[w]*cn; 
		} 

	} 
	return(1); 
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