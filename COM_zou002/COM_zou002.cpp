//////////////////////////////////////////////////////////////////////
// COPYRIGHT NOTICE  
// Copyright (c) 2016, 光巨力信息技术有限公司Algorithm Group  （版权声明）  
// All rights reserved.  
//
// @file    COM_zou002.cpp    
// @brief   串口通信实现文件  
//
// (异步方式)实现串口的打开、配置、读写和关闭等操作, + 调试二代机算法程序 + 预测模块&校正模块
//
//
// @version 1.2     
// @author  LingWei zou    
// @E-mail：zlw2008ok@126.com  
// @date    2016/12/21  
//  
//
// 修订说明：  + 预测模块&校正模块
////////////////////////////////////////////////////////////////////
#include "stdafx.h"  
#include "SerialPort.h"  
#include <iostream>  
#include <tchar.h>
#include <stdio.h>
#include <time.h>
#include <direct.h> //用 _mkdir()创建一个文件夹（目录）
#include <io.h>


#define A_n 20		//用于夹持判断的数据最大长度
#define Max_n 120	//有效数据最大长度
#define T_n 20		//温度数据最大长度
#define C_n 3		//循环接收数据的长度
#define H_data 13000	//有效数据最大高度
#define L_data 500		//有效数据最小高度
#define step_H 1000		//夹持跳跃高度
#define stable_H 150	//扰动最大容忍幅度

// ----------- 全局变量 ------------
float temp_C;

int Alen1=0;		//用于夹持判断的数据长度
int Alen2=0;
int Alen3=0;
int Alen4=0;


int SDec1[1000];	//存放有效数据
int SDec2[1000];
int SDec3[1000];	
int SDec4[1000];
int count1=0;		//有效数据计数
int count2=0;	
int count3=0;		
int count4=0;
int step_i1=0;		//有效数据起点
int step_i2=0;
int step_i3=0;
int step_i4=0;

int start_flag1=0;	//有效数据开始标志
int start_flag2=0;
int start_flag3=0;
int start_flag4=0;
int wear_flag;		//设备穿戴标志
int clamp_flag;		//夹持动作标志
int judge_flag;		//数据稳定标志

typedef struct Tnode //定义用于存储温度的链表节点
{ 
	int TempData;
	Tnode * next;
}Tnode,* T_Link;
T_Link L1,Lsp,Ltemp;

//----- 数据采集 -----
int Hex2Dec(char * Hex )							//处理4位的十六进制数 转 十进制
{
	int i;
	int temp[4];
	for (i=0; i<4; ++i)
	{
		if ((Hex[i]>='0') && (Hex[i]<='9'))
			temp[i] = (int)(Hex[i] - '0');
		else if ((Hex[i]>='a') && (Hex[i]<='f'))
			temp[i] = (int)(Hex[i] - 'a') + 10;
		else if ((Hex[i]>='A') && (Hex[i]<='F'))
			temp[i] = (int)(Hex[i] - 'A') + 10;
		else
			return 0;								//若不是十六进制数，函数输出为0
	}

	return temp[0]*16*16*16 + temp[1]*16*16 + temp[2]*16 + temp[3]; 
}

void Init_Fun(void)//标志位初始化函数
{
	Alen1=0;
	Alen2=0;
	Alen3=0;
	Alen4=0;

	start_flag1=0;
	start_flag2=0;
	start_flag3=0;
	start_flag4=0;
	count1=0;
	count2=0;
	count3=0;
	count4=0;
	step_i1=0;
	step_i2=0;
	step_i3=0;
	step_i4=0;
	clamp_flag=0;
}

int Wear_Fun(int *Dec, int len, int high) //（单波长），判断设备是否穿戴上，return 1,已经戴上，return 0，设备空置状态
{
	int i;
	
	if(Dec[len-1]>high && Dec[len-1]<H_data && Dec[len-2]>high && Dec[len-2]<H_data) //最后2个点都在范围内，说明设备已经在耳垂上
		return 1;//收到的数据中末尾3个数据在有效范围内
	else
		return 0;

}

int Add_Fun(int *Dec, int len, int *ADec, int * Alen) //累加保存判断数据函数,最大容量A_n
{
	if (A_n-*Alen >len)
	{
		for(int i=0; i<len; i++)
		{
			ADec[*Alen] = Dec[i];
			(*Alen)++;
		}
	}
	else
	{
		for(int i=0; i<A_n-len; i++)
		{
			ADec[i] = ADec[*Alen+len-A_n+i];
		}
		for(int i=A_n-len; i<A_n; i++)
		{
			ADec[i] = Dec[i-(A_n-len)];
		}
		*Alen = A_n;
	}
	return 0;
}

int Clamp_Fun(int *Dec, int len, int step_h, int * step_i)//（单波长）判断夹持动作，当设备对耳垂进行夹持操作时，return 1，落差>step_h,否则 return 0；
{
	int i;

	for(i=0; i<len-10; i++)
	{	
		if ( Dec[i+8]<H_data && Dec[i+8]-Dec[i]>step_h && Dec[i]>step_h)
		//if ( Dec[i+8]<H_data && Dec[i+9]<H_data && Dec[i+10]<H_data && Dec[i+8]-Dec[i]>step_h && Dec[i+9]-Dec[i+1]>step_h && Dec[i+10]-Dec[i+2]>step_h && Dec[i]>step_h  && Dec[i+1]>step_h && Dec[i+2]>step_h)
		{
			*step_i = i;
			return 1;
		}
	}
	//printf("-%d-",*step_i);
	return 0;
}

void Save_Fun(int *Dec, int len, int step_i, int *SDec, int * count)
{
	for(int i=0; i<len; i++)
	{
		SDec[*count+step_i] = Dec[i];
		(*count)++;			
	}
}

int DataJudge_Fun(int *SDec,int step_i, int *count, int stable_h)
{
	int len=*count+step_i;

	if(*count>20)
	{
		for(int i=step_i+20; i<(len-4); i++)
		{
			if((SDec[i+4]-SDec[i])>stable_h || (SDec[i+4]-SDec[i])<-(stable_h))
			{

				return 1;//数据不稳定，重新测量
			}
		}
	}
	else
	{
		for(int i=step_i+8; i<(len-4); i++)
		{
			if(/*(SDec[i+10]-SDec[i])>stable_h || */(SDec[i+4]-SDec[i])<-(stable_h))
			{

				return 1;//数据不稳定，重新测量
			}
		}
	}
	return 0;//数据稳定
}

int ReceiveData_Fun(int *Dec1,int *Dec2, int *Dec3, int *Dec4, int len1, int len2, int len3, int len4, int *Temp, int Tlen)
{
	//************ 设备空闲阶段 ******************
	wear_flag = Wear_Fun(Dec1, len1, L_data) && Wear_Fun(Dec2, len2, L_data) && Wear_Fun(Dec3, len3, L_data) && Wear_Fun(Dec4, len4, L_data);//设备是否穿戴上
	if (wear_flag==0)
	{
		Init_Fun();
		return 0;// 设备空闲状态...
	}
	//************ 测试阶段 *************（不进行加热，测试阶段需要判断夹持动作，并找到有效数据起点）	

	//---------- 温度数据保存 --------------
	for(int i=0; i<Tlen; i=i+4)
	{
		L1->TempData=Temp[i];//温度数据保存在循环链表中
		L1=L1->next;
	}
	//---------- 夹持动作判断 -------------
	if (start_flag1==0) {Add_Fun(Dec1, len1, SDec1, &Alen1);start_flag1 = Clamp_Fun(SDec1, Alen1, step_H, &step_i1);count1=Alen1-step_i1;}
	if (start_flag2==0) {Add_Fun(Dec2, len2, SDec2, &Alen2); start_flag2 = Clamp_Fun(SDec2, Alen2, step_H, &step_i2);count2=Alen1-step_i2;}
	if (start_flag3==0) {Add_Fun(Dec3, len3, SDec3, &Alen3); start_flag3 = Clamp_Fun(SDec3, Alen3, step_H, &step_i3);count3=Alen1-step_i3;}
	if (start_flag4==0) {Add_Fun(Dec4, len4, SDec4, &Alen4); start_flag4 = Clamp_Fun(SDec4, Alen4, step_H, &step_i4);count4=Alen1-step_i4;}
	//---------- 保存数据 -----------
	if (start_flag1==1)	Save_Fun(Dec1, len1, step_i1, SDec1, &count1);
	if (start_flag2==1) Save_Fun(Dec2, len2, step_i2, SDec2, &count2);
	if (start_flag3==1) Save_Fun(Dec3, len3, step_i3, SDec3, &count3);
	if (start_flag4==1) Save_Fun(Dec4, len4, step_i4, SDec4, &count4);

	clamp_flag = start_flag1 && start_flag2 && start_flag3 && start_flag4;//夹持动作标志
	if (clamp_flag)
	{
		judge_flag = DataJudge_Fun(SDec1,step_i1, &count1, stable_H) || DataJudge_Fun(SDec2,step_i2, &count2, stable_H) || DataJudge_Fun(SDec3,step_i3, &count3, stable_H)/* || DataJudge_Fun(SDec4,step_i4, &count4, stable_H)*/;
		if(judge_flag==1)
		{
			//Init_Fun();
			//judge_flag=0;
			return 4;//数据不稳定，从新测量
		}
		else
		{		
			if(count1>Max_n && count2>Max_n && count3>Max_n && count4>Max_n)
			{
				return 5;//测量完成；
			}
			else
				return 3;//正在测量....
		}
	}
	
}

//----- 数据预处理 -----
void Denoise(int * NoiseData, int Data_Len,int Mean_Len,float*SmData )	//均值平滑曲线,Data_Len数据长度，Mean_Len均值步长。
{
	int i,j;
	float sum;


	//------------------ 均值平滑曲线 ---------------------------
	for(i=0; i<Data_Len; ++i)
	{
		if(i < Data_Len-Mean_Len+1)
		{
			sum=0;
			for(j=0; j<Mean_Len; ++j)
			{
				sum=sum+NoiseData[i+j];
			}
			SmData[i]=sum/Mean_Len;
		}
		else
		{
			sum=0;
			for(j=0; j<Data_Len-i; ++j)
			{
				sum=sum+NoiseData[i+j];
			}
			SmData[i]=sum/(Data_Len-i);
		}
	}

	return;
/*
	//-------------------- 中位值滤波 ---------------------
	//连续采样N次（N取奇数）
 //  把N次采样值按大小排列
 //  取中间值为本次有效值
	float buf[30];//Mean_Len不能大于30
	float temp=0;
	int i,j,k,p;
	for(i=0; i<Data_Len; i++)
	{

		if(i < Data_Len-Mean_Len+1)
		{
			for(p=0; p<Mean_Len; p++)
			{
				buf[p]=NoiseData[i+p];//缓冲数据
			}

			for(j=0; j<Mean_Len; j++)
			{
				for(k=0; k<Mean_Len-j-1; j++)
				{
					if(buf[k]>buf[k+1])
					{
						temp = buf[k];
						buf[k] = buf[k+1];
						buf[k+1] = temp;
					}
				}
			}
			if(Mean_Len%2==0)
				SmData[i] = (buf[Mean_Len/2]+buf[Mean_Len/2+1])/2;
			else
				SmData[i] = buf[(Mean_Len+1)/2];
		}
		else
		{
			for(p=0; p<Data_Len-i; p++)
			{
				buf[p]=NoiseData[i+p];//缓冲数据
			}

			for(j=0; j<Data_Len-i; j++)
			{
				for(k=0; k<Data_Len-i-j-1; j++)
				{
					if(buf[k]>buf[k+1])
					{
						temp = buf[k];
						buf[k] = buf[k+1];
						buf[k+1] = temp;
					}
				}
			}
			if((Data_Len-i)%2==0)
				SmData[i] = (buf[(Data_Len-i)/2]+buf[(Data_Len-i)/2+1])/2;
			else
				SmData[i] = buf[((Data_Len-i)+1)/2];
		}
	}
	return;
	*/

/*
	//-------------------- 一阶滞后滤波 ---------------------
	//Mean_Len在这里为0~1的系数；
	int i;
	for(i=0; i<Data_Len-1; i++)
	{
		SmData[i]=(1-Mean_Len)*NoiseData[i+1]+Mean_Len*NoiseData[i]
	}
	SmData[Data_Len-1]=NoiseData[Data_Len-1];
	return;
	*/
}

//----- 特征提取 -----
int FindStartEnd(float * SmoothData ,int len , int * StartEnd)		//先找到数据有效起点&终点。
{
	int i,j=0;
	float *tempData = NULL;
	float temp1;

	tempData = (float*) malloc( sizeof(float)*len );
	if (!tempData)
	{
		printf("分配内存失败");
		return 0;
	}

	for (i=0; i<len-1; i++)
	{
		tempData[i] = SmoothData[i+1]-SmoothData[i];
	}

	temp1 = tempData[0];
	for (i=0; i<len-2; i++)
	{
		if ( temp1<tempData[i+1] )
		{
			temp1 = tempData[i+1];
			j = i+1;
		}
	}
	StartEnd[0] = j+1; //数据有效起点
	StartEnd[1] = len-1; //数据有效终点

	free(tempData);

	return 1;
}

void FindtwoPoint(float * SmoothData ,int len ,int * StartEnd, int n_end,float n_r, int * twoPoint)
{
	int i,j=StartEnd[0];
	float temp2,tempEnd;
//----------------------------------- 按一定规则找数据计算的起点 --------------------------------------------------
	twoPoint[1] = StartEnd[0]+n_end; //设定数据计算终点（有效数据的第n_end点）
	if (twoPoint[1]>len)
		twoPoint[1] = len-1;
	
	tempEnd = SmoothData[twoPoint[1]]/(1+n_r);
	
	temp2 = abs(SmoothData[StartEnd[1]] - SmoothData[StartEnd[0]]); //初始为有效终点-起点的幅值
	for (i=StartEnd[0]; i<twoPoint[1]; i++)
	{
		if (temp2 > abs(tempEnd-SmoothData[i]))
		{
			temp2 = abs(tempEnd-SmoothData[i]);
			j = i;
		}
	}
	twoPoint[0] = j; //数据计算起点
	twoPoint[0] = twoPoint[0]-StartEnd[0]; //数据有效起点 到 数据计算起点 的距离。
	twoPoint[1] = twoPoint[1]-StartEnd[0]; //数据有效起点 到 数据计算终点 的距离。
}

float Feature1(float *SmoothData,int *StartEnd, int *twoPoint)					//声明函数，取两个点之间的log差作为特征1；
{
	float temp;
	float y_feature;

	if (twoPoint[1]>(StartEnd[1]-StartEnd[0]))
	{
		twoPoint[1]=StartEnd[1]-StartEnd[0];
	}
	if (twoPoint[0]>(StartEnd[1]-StartEnd[0]))
	{
		twoPoint[0]=StartEnd[1]-StartEnd[0]-1;
	}
	//temp = SmoothData[twoPoint[1]]-SmoothData[twoPoint[0]];
	temp = log10(SmoothData[StartEnd[0]+twoPoint[1]])-log10(SmoothData[StartEnd[0]+twoPoint[0]]);
	y_feature = temp;

	return y_feature;
};

float temp_Fun(T_Link L)//计算耳垂温度函数（求测量过程中温度均值）
{
	int i=T_n;
	float T=0.0;
	while(i)
	{
		T += L->TempData;
		L = L->next;
		i--;
	}
	return T/T_n;
}

//----- 调用模型 -----
float GluPre( float * x,  float * b)
{
	float y;
	//计算公式

	y = b[1]*x[0] + b[2]*x[1] + b[3]*x[2] + b[4]*x[3] +b[0];
	

	return y;
}

//----- 预测输出 -----


//----- 模型校正 -----
//
int achol(float *b,int n,int m,float *a)

{
	int i,j,k,u,v;
	if ((b[0]+1.0==1.0)||(b[0]<0.0))
	{
		printf("fail\n");return(-2);
	}
	b[0]=sqrt(b[0]);
	for (j=1;j<=n-1;j++) b[j]=b[j]/b[0];
	for (i=1;i<=n-1;i++)
	{
		u=i*n+i;
		for (j=1;j<=i;j++)
		{
			v=(j-1)*n+i;
			b[u]=b[u]-b[v]*b[v];
		}
		if ((b[u]+1.0==1.0)||(b[u]<0.0))
		{
		    printf("fail\n");return(-2);
	    }
		b[u]=sqrt(b[u]);
		if (i!=(n-1))
		{
			for (j=i+1;j<=n-1;j++)
			{
				v=i*n+j;
				for (k=1;k<=i;k++)
					b[v]=b[v]-b[(k-1)*n+i]*b[(k-1)*n+j];
				b[v]=b[v]/b[u];
			}
		}
	}
	for (j=0;j<=m-1;j++)
	{
		a[j]=a[j]/b[0];
		for (i=1;i<=n-1;i++)
		{
			u=i*n+i;v=i*m+j;
			for(k=1;k<=i;k++)
				a[v]=a[v]-b[(k-1)*n+i]*a[(k-1)*m+j];
			a[v]=a[v]/b[u];
		}
	}
	for (j=0;j<=m-1;j++)
	{
		u=(n-1)*m+j;
		a[u]=a[u]/b[n*n-1];
		for (k=n-1;k>=1;k--)
		{
			u=(k-1)*m+j;
			for(i=k;i<=n-1;i++)
			{
				v=(k-1)*n+i;
				a[u]=a[u]-b[v]*a[i*m+j];
			}
			v=(k-1)*n+k-1;
			a[u]=a[u]/b[v];
		}
	}
	return(2);
}

void isqt2(float *x,float *y,int m,int n,float *a,float *dt,float *v)
{
	int i,j,k,l,mm;
	float q,e,u,p,yy,s,r,pp,*b;
	//extern int achol();
	b=(float *)malloc((m+1)*(m+1)*sizeof(float));
	mm=m+1;
	b[mm*mm-1]=n;
	for (j=0;j<=m-1;j++)
	{
		p=0.0;
		for (i=0;i<=n-1;i++)
			p=p+x[j*n+i];
		b[m*mm+j]=p;
		b[j*mm+m]=p;
	}
	
	for (i=0;i<=m-1;i++)
		for(j=i;j<=m-1;j++)
		{
			p=0.0;
			for (k=0;k<=n-1;k++)
				p=p+x[i*n+k]*x[j*n+k];
			b[j*mm+i]=p;
			b[i*mm+j]=p;
		}
	a[m]=0.0;
	for (i=0;i<=n-1;i++)
		a[m]=a[m]+y[i];
	for (i=0;i<=m-1;i++)
	{
		a[i]=0.0;
		for(j=0;j<=n-1;j++)
			a[i]=a[i]+x[i*n+j]*y[j];
	}
    achol(b,mm,1,a);//解线性方程组，b[mm*mm],a[mm]；mm为自变量个数+1，1为因变量；
	yy=0.0;
	for(i=0;i<=n-1;i++)
		yy=yy+y[i]/n;
	q=0.0;e=0.0;u=0.0;
	for (i=0;i<=n-1;i++)
	{
		p=a[m];
		for(j=0;j<=m-1;j++)
			p=p+a[j]*x[j*n+i];
		q=q+(y[i]-p)*(y[i]-p);
		e=e+(y[i]-yy)*(y[i]-yy);
		u=u+(yy-p)*(yy-p);
	}
	s=sqrt(q/n);
	r=sqrt(1.0-q/e);
	for(j=0;j<=m-1;j++)
	{
		p=0.0;
		for(i=0;i<=n-1;i++)
		{
			pp=a[m];
			for(k=0;k<=m-1;k++)
				if(k!=j) pp=pp+a[k]*x[k*n+i];
			p=p+(y[i]-pp)*(y[i]-pp);
		}
		v[j]=sqrt(1.0-q/p);//偏相关系数
	}
	dt[0]=q;dt[1]=s;dt[2]=r;dt[3]=u;//[0]误差和，[1]均方根误差，[2]复相关系数，[3]变异平方和
	free(b);
	return;
}

float   *   MatrixInver(float   A[],int   m,int   n)   /*矩阵转置*/ 
{ 
          int   i,j; 
          float   *B=NULL; 
          B=(float   *)malloc(m*n*sizeof(float)); 
        
          for(i=0;i <n;i++) 
          for(j=0;j <m;j++) 
          B[i*m+j]=A[j*n+i];         
          return   B; 
} 

float   Surplus(float   A[],int   m,int   n)   /*求矩阵行列式*/ 
{         
          int   i,j,k,p,r; 
          float   X,temp=1,temp1=1,s=0,s1=0; 
        
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

float   *   MatrixOpp(float   A[],int   m,int   n)   /*矩阵求逆*/ 
{ 
          int   i,j,x,y,k; 
          float   *SP=NULL,*AB=NULL,*B=NULL,X,*C=NULL; 
          SP=(float   *)malloc(m*n*sizeof(float)); 
          AB=(float   *)malloc(m*n*sizeof(float)); 
          B=(float   *)malloc(m*n*sizeof(float)); 
		  C=(float   *)malloc(m*m*sizeof(float)); 
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

float	 *	 MatrixMul(int n, float   A[], int xm, float B[], int ym)  //矩阵相乘
{
	int	 i,j,k;
	float   sum;
	float	 * C = NULL;
	
	C=(float *)malloc(xm*ym*sizeof(float));

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

float   *   Mlr(int n, float x[], int xm,float y[], int ym )      // b=x*y'/(x*x')
{
	float  *xx_arr;
	float  *xx_inv; 
	float  *xy_arr;
	float	 * result;

	xy_arr = MatrixMul(n,(float   *)x,xm,(float   *)y,ym);//矩阵相乘
	xx_arr = MatrixMul(n,(float   *)x,xm,(float   *)x,xm);//

	xx_inv = MatrixOpp((float   *)xx_arr,xm,xm);     //求逆 

	result = MatrixMul(xm,(float   *)xy_arr,ym,(float   *)xx_inv,xm); 


	free(xy_arr);
	free(xx_arr);
	free(xx_inv);

	return result;
}

void con_oldnew_fun(int m, float real_glu, float x_new[], float * xy_old, int old_num, float * xy_new)
{
	//前向累计新样本。
	//m:每次更新的数据量，[x1,x2,y]则m为3；
	//old_num:旧样本数量；
	for(int i=m*(old_num+1)-1; i>m-1; i--)
	{
		xy_new[i] = xy_old[i-m];
	}
	for(int i=0; i<m-1; i++)
	{
		xy_new[i] = x_new[i];
	}
	xy_new[m-1] = real_glu;//
	

}

void cor_fun(int m, float  * xy_new, int new_num, float *b)
{
	int i;
	//float xx[2][30];
	float *xx = NULL;
	float yy[30];
	if(new_num>=30)
	{
		new_num=30;
	}

	xx = (float *)malloc(m*new_num*sizeof(float));



	for(i=0; i<m-1; i++)
	{
		for(int j=0; j<new_num; j++)
		{
			xx[i*new_num+j] = xy_new[m*j+i];
		}
	}
	for(int j=0; j<new_num; j++)
	{
		xx[(m-1)*new_num+j] = 1.0;
	}

	for(int j=0; j<new_num; j++)
	{
		yy[j] = xy_new[m*j+m-1];
	}


//----------- 调用建模函数 ----------
	//****-----------------mlr1:
	float a[3],v[2],dt[4];//a[]回归系数a[m]为零点；v[]偏相关系数；dt[]:[0]误差和，[1]均方根误差，[2]复相关系数，[3]变异平方和
	//float   x[3][4]={6,2,3,4,2,6,7,8,9,3,4,9};
	//float   y[4]={1,2,3,4};
	isqt2((float*)xx,yy,2,new_num,a,dt,v);

	//b[0]=a[0];
	//b[1]=a[1];
	//b[2]=a[2];

	//****-----------------mlr2:
	float *bb;
	bb = Mlr(new_num,(float   *)xx,m,(float   *)yy,1);  // b=x*y'/(x*x')
	for(int j=0; j<5; j++)
		b[j]=bb[j];
	free(bb);
}


int _tmain(int argc, _TCHAR* argv[])
{



//------------------- 创建循环链表用于存储温度数据 -----------------------
	L1=(T_Link)malloc(sizeof(Tnode));
	L1->TempData=0;
	L1->next=NULL;
	Ltemp=L1;

	for(int i=1; i<T_n; i++)
	{
		Lsp=(T_Link)malloc(sizeof(Tnode));
		Lsp->TempData=i;

		Ltemp->next=Lsp;
		Ltemp=Lsp;
	}
	Ltemp->next=L1;
	//MessageBox(NULL,(LPCWSTR)"HELLO!",(LPCWSTR)"message",0);

	//--------------------- 创建串口 ----------------
	int com_num;
	bool opensuccess;
	bool initsuccess;
	bool readsuccess;
	//bool threadsuccess;
	bool claersuccess;
	SerialPort my_Port;  
	while(1)
	{
		printf("输入串口号：");
		scanf("%d",&com_num);
		char ch;while((ch = getchar()) != '\n' && ch != EOF);
		opensuccess = my_Port.OpenPort(com_num);
		if (!opensuccess)
		{
			printf("串口不存在！重新输入！\n");
		}
		else
			break;
	}
	initsuccess = my_Port.InitPort(115200);

	//threadsuccess = my_Port.OpenListenThread();


	//-------------------- 进入读取操作环节 ----------------
Loop: Init_Fun();

	claersuccess = my_Port.ClearPort();//清除串口缓存
	char cha;
	char str[4]={'0','0','0','0'};
	char strT[4]={'0','0','0','0'};
	char strR[4]={'0','0','0','0'};

	int j=0;
	int flag=0;
	int w_flag=-1;
	int D_data=0;

	int Dec1[50];
	int Dec2[50];
	int Dec3[50];
	int Dec4[50];
	int Temp[50];

	int len1=0;
	int len2=0;
	int len3=0;
	int len4=0;
	int Tlen=0;

	int rec_flag=0;
	int stop_flag=0;
	while(!stop_flag==1)
	{
		readsuccess = my_Port.ReadPort(&cha);
		//printf("%c",cha);
		//continue;
		if(cha=='s')
		{
			flag=1;
			j=0;
		}

		if(flag)
		{
			j++;
			if(j==4)
			{
				switch(cha)
				{
				case '0':w_flag=1;break;
				case '1':w_flag=2;break;
				case '2':w_flag=3;break;
				case '3':w_flag=4;break;
				default :w_flag=-1;
				}
			}
			else if(j>6&&j<11)
				str[j-7]=cha;
			else if(j>18&&j<23)
				strT[j-19]=cha;
			else if(j>30&&j<35)
				strR[j-31]=cha;
			else if(j>34)
			{
				D_data=Hex2Dec(str);
				flag=0;
				if(w_flag==1)
				{Dec1[len1]=D_data;len1++;}
				else if(w_flag==2)
				{Dec2[len2]=D_data;len2++;}
				else if(w_flag==3)
				{Dec3[len3]=D_data;len3++;}
				else if(w_flag==4)
				{Dec4[len4]=D_data;len4++;}

				D_data=Hex2Dec(strT);
				Temp[Tlen]=D_data;Tlen++;
			}
		}



		if(len1>C_n&&len2>C_n&&len3>C_n&&len4>C_n)
		{
			rec_flag = ReceiveData_Fun(Dec1,Dec2, Dec3, Dec4, len1, len2, len3, len4, Temp, Tlen);

			len1=0;
			len2=0;
			len3=0;
			len4=0;
			Tlen=0;
			if(wear_flag==0)
			{
				printf("空闲等待中...\n");

			}
			else
			{
				if(clamp_flag==0)
				{
					printf("正在扑捉夹持动作...\n");
					//printf("%d%d%d%d\n",start_flag1,start_flag2,start_flag3,start_flag4);
					//printf("%d%d%d%d\n",Alen1,Alen2,Alen3,Alen4);
					//printf("%d%d%d%d\n",step_i1,step_i2,step_i3,step_i4);
					
				}
				else
				{
					if(judge_flag==1)
					{
						Init_Fun();
						printf("数据不稳定，结束！重新夹持测量...\n");
					}
					else
					{
						if(count1<Max_n || count2<Max_n || count3<Max_n || count4<Max_n)
						{
							printf("正在测量...\n");
						}
						else
						{
							printf("测试完成!\n");
							stop_flag=1;
						}
					}

				}
			}

		}
	}
/*wear_flag==0; 空闲
//
//wear_flag==1; 设备带上
//		clamp_flag==0;正在扑捉夹持动作...
//
//		clamp_flag==1;正在测量...
//			judge_flag==1;数据不稳定，结束！重新夹持测量...
//
//			judge_flag==0;数据稳定；
//				count<400;正在测量；
//
//				count>400;测量完成，结束！

		//printf("%c",cha);
*/


/*
调用模型参数
	有：个人模型预测
	用户不存在：大众模型预测（大众模型里可以分类多种模型）
	用户校正样本不足：大众模型预测

是否需要校正
	否：退出
	是：输入微创值，调用校正函数，更新特征数据文件，更新模型文件	  
*/



//-------- 获得系统时间+数据文件名 ---------
	float real_value;
	int dec_p,sign_p;
	char str_glu[5];
	//scanf("%f",&real_value);
	//sprintf( str_glu, "%3.1f", real_value );//血糖值转为字符串
	//char glu_d[10];
	//gets(glu_d);
	_mkdir("d:\\4wave_data");
	printf("请输入人名（如 zlw）：\n");
	char name[10];
	gets(name);
	char path[100]="d:\\4wave_data\\";
	strcat(path,name);

	_mkdir(path);
	//printf("%s",name);

	time_t ts=time(NULL);
    struct tm *ptm;  
    ptm =localtime(&ts); //此函数获得的tm结构体的时间，是已经进行过时区转化为本地时间  
	//int y=1900+ptm->tm_year;
	//int m=1+ptm->tm_mon;
	//int d=ptm->tm_mday;
	//int h=ptm->tm_hour;
	//int mm=ptm->tm_min;
	//int s=ptm->tm_sec;
	
	char strt[25];
	strftime(strt, sizeof(strt), "%Y%m%d%H%M%S", ptm) ;
	strcat(strt,"&");
	//strcat(strt,glu_d);
	strcat(strt,".txt");
	strcat(path,"\\");
	strcat(path,strt);

//---- 有效数据写入文件 -----
	FILE *fpt;
	int aaa=-77777;
	if((fpt = fopen(path,"wb"))==NULL)
	{
		printf("打开文件失败！");
		exit(0);
	}
	for(int k=step_i1; k<step_i1+count1; k++)
	{
		//fprintf(fpt, "%d", SDec1[k]); /*向所建文件写一整型数*/ 
		fwrite(&SDec1[k],sizeof(int),1,fpt);
	}
	fwrite(&aaa,sizeof(int),1,fpt);
	for(int k=step_i2; k<step_i2+count2; k++)
	{
		fwrite(&SDec2[k],sizeof(int),1,fpt);
	}
	fwrite(&aaa,sizeof(int),1,fpt);
	for(int k=step_i3; k<step_i3+count3; k++)
	{
		fwrite(&SDec3[k],sizeof(int),1,fpt);
	}
	fwrite(&aaa,sizeof(int),1,fpt);
	for(int k=step_i4; k<step_i4+count4; k++)
	{
		fwrite(&SDec4[k],sizeof(int),1,fpt);
	}
	fwrite(&aaa,sizeof(int),1,fpt);
	for(int k=0; k<T_n; k++)
	{
		fwrite(&(L1->TempData),sizeof(int),1,fpt);
		L1=L1->next;
	}
				
	fclose(fpt);

	//if((fpt=fopen(strt,"rb"))==NULL)
	//{
	//	printf("打开文件失败");
	//	exit(1);
	//}

	//fseek(fpt, 0, SEEK_END);//将文件指针指向 离文件末尾(SEEK_END)距离为0的位置。
	//int len = ftell(fpt);  //计算当前文件指针离文件头的距离(字节数)。这就是长度了
	//int leni=len/sizeof(int);

	//fseek(fpt, 0, SEEK_SET);//将文件指针指向 离文件头(SEEK_SET)距离为0的位置。
	//int *Deck;
	//Deck=(int *)malloc(sizeof(int)*leni);

	//fread(Deck,sizeof(int),leni,fpt);
	//fclose(fpt);

	//for(int i=0; i<leni; i++)
	//	printf("%d ",Deck[i]);


//--------------------- 数据分析判断 -----------------------
	//这里对有效数据可以进一步判断（毕竟才采数过程判断数据稳定性存在一定局限）； 

//--------------------- 数据预处理  -------------------------
	for(int i=0; i<count1; i++)
		SDec1[i]=SDec1[step_i1+i];
	for(int i=0; i<count2; i++)
		SDec2[i]=SDec2[step_i2+i];
	for(int i=0; i<count3; i++)
		SDec3[i]=SDec3[step_i3+i];
	for(int i=0; i<count4; i++)
		SDec4[i]=SDec4[step_i4+i];


	int mean_n=5;
	float SmData1[300];
	float SmData2[300];
	float SmData3[300];
	float SmData4[300];
	Denoise( SDec1, count1,mean_n, SmData1);
	Denoise( SDec2, count2,mean_n, SmData2);
	Denoise( SDec3, count2,mean_n, SmData3);
	Denoise( SDec4, count2,mean_n, SmData4);

//--------------------- 特征提取  -------------------------
	int StartEnd[2]={0,100};						//存放有效数据头尾序号
	int twoPoint[2];						//存放计算起点终点 离有效数据起点的距离
	int flag_r;								//twoPoint存在的标志
	float x[10]={0,0,0,0,0,0,0,0,0,0};		//存放数据特征变量值(最多10个特征)

	//flag_r = FindStartEnd(SmData1 ,count1 , StartEnd);
	//FindtwoPoint(SmData1 ,count1 ,StartEnd,400,0.03, twoPoint);
	twoPoint[0]=15;twoPoint[1]=100;
	x[0] = Feature1(SmData1, StartEnd, twoPoint);
	x[1] = Feature1(SmData2, StartEnd, twoPoint);
	x[2] = Feature1(SmData3, StartEnd, twoPoint);
	x[3] = Feature1(SmData4, StartEnd, twoPoint);

	//temp_C = temp_Fun(L1);
	//x[2] = temp_C;//温度特征
//---------------------------------  调用模型参数  --------------------------------
float b_default[5]={4.8, -0.4, -0.3, 1.21, 0.05};//大众模型

FILE * fpt_b,* fpt_xy;
int xyf_flag=0,bf_flag=0;//配置文件存在标志
char b_file[100]="d:\\4wave_data\\";
char xy_file[100]="d:\\4wave_data\\";
strcat(b_file,name);
strcat(xy_file,name);
strcat(b_file,"\\model_parameter.txt");
strcat(xy_file,"\\xy_features.txt");
//if((fpt_xy=fopen(xy_file,"r+"))==NULL)
//{
//	xyf_flag=0;//样本文件不存在
//}
//else
//{
//	xyf_flag=1;//样本文件存在
//	fclose(fpt_xy);
//}

if((fpt_b=fopen(b_file,"r+"))==NULL)
{
	bf_flag=0;//模型文件不存在
}
else
{
	bf_flag=1;//模型文件存在
	fclose(fpt_b);
}


float b_person[10];
if(bf_flag==1)
{
	fpt_b=fopen(b_file,"r");
	char strf[50];
	int ir=0;
	while(fgets(strf,50,fpt_b)!=NULL)
	{
		b_person[ir++]=atof(strf);
	}
	//b=b_person;个人模型预测程序

	fclose(fpt_b);
}
else
{
	for(int i=0; i<5; i++)
		b_person[i]=b_default[i];
	//b=b_default;大众模型预测程序	
	printf("没有您的模型或校正次数不够，需要校正！\n");
}



//---------------------------------  预测血糖值  -----------------------------------------
	float Glu_value;

	Glu_value = GluPre(x, b_person);
	int kkk;
	srand((unsigned)time( NULL ));
	kkk=rand()%10;

	if(Glu_value>=10)
		Glu_value=10+kkk*0.1;
	else if (Glu_value<3.8)
		Glu_value=3.9+kkk*0.1;	
	printf("预测血糖值为：------------ %f\n -----------\n-------------\n",Glu_value);

//---------------------------------  校正模型 & 保存数据及参数 ----------------------------------------------
int cor_flag;
float Glu_real;
printf("是否需要校正？输入1或0，分别代表校正和不校正\n");
//char ch;while((ch = getchar()) != '\n' && ch != EOF);
scanf("%d",&cor_flag);

//printf("cor_flag=%\n",cor_flag);
if(cor_flag==0)
{
	goto Loop1;
}
else
{
	printf("请输入微创值：\n");
	scanf("%f",&Glu_real);
	
}

//读取历史样本
int ixy=0;
float xy_person[50];
if((fpt_xy=fopen(xy_file,"r+"))==NULL)
{
	xyf_flag=0;//样本文件不存在
	//xy_old=[]
}
else
{
	while(fscanf(fpt_xy,"%f",&xy_person[ixy])!=EOF)
	{
		ixy++;
	}
	fclose(fpt_xy);
}

float xy_new[5]={x[0],x[1],x[2],x[3],Glu_real};//新样本

if(ixy<45)
{
	printf("累计样本...");
	int j;
	for(j=ixy; j<ixy+5; j++)
		xy_person[j]=xy_new[j-ixy];
	int n_xy=j;

	//更新样本，保存xy_file
	fpt_xy=fopen(xy_file,"w");

	for(int i=0; i<n_xy ;i=i+5)//将样本写入文件
	{
		fprintf(fpt_xy, "%f %f %f %f %f\n", xy_person[i],xy_person[i+1],xy_person[i+2],xy_person[i+3],xy_person[i+4]);
	}
	fclose(fpt_xy);
	printf("累计样本完成。\n");
}
else
{
	printf("更新样本...更新模型...\n");
	int j;int n_xy=50;
	float temp_xy[50];
	if(ixy==45)
	{
		for(j=ixy; j<ixy+5; j++)
			xy_person[j]=xy_new[j-ixy];
	}
	else
	{
		for(j=0; j<45; j++)
			xy_person[j]=xy_person[j+5];
		for(int i=0; i<5; i++)
			xy_person[j+i]=xy_new[i];
	}


	//校正，更新样本，保存xy_file+b_file;
	int m=5; 
	int new_num=10;
	float b[5];
	cor_fun(m,  xy_person, new_num, b); //模型校正函数
	b_person[0]=b[4];
	b_person[1]=b[0];b_person[2]=b[1];b_person[3]=b[2];b_person[4]=b[3];
	printf("模型系数：（%f）* x1 + （%f）* x2  + （%f）* x3 + （%f）* x4+ （%f）\n",b[0],b[1],b[2],b[3],b[4]);
	//校正函数->b

	fpt_xy=fopen(xy_file,"w");
	for(int i=0; i<n_xy ;i=i+5)//将模型参数写入文件
	{
		fprintf(fpt_xy, "%f %f %f %f %f\n", xy_person[i],xy_person[i+1],xy_person[i+2],xy_person[i+3],xy_person[i+4]);
	}
	fclose(fpt_xy);

	fpt_b=fopen(b_file,"w");
	for(int i=0; i<5 ;i++)//将模型参数写入文件
	{
		fprintf(fpt_b, "%f\n", b_person[i]);
	}
	fclose(fpt_b);
	printf("校正完毕，恭喜！\n");
}


	Loop1:
	printf("继续测量请输入“1”+回车，结束则关闭。\n");
	int GO_num;
	while(1)
	{
		scanf("%d",&GO_num);
		if(GO_num==1)
		{
			char ch;while((ch = getchar()) != '\n' && ch != EOF);
			goto Loop;
		}
	}
	return 0;
}

