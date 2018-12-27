// filter01.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include <stdio.h>
#include <math.h>
#include <malloc.h>

#include "kfft.h"

int main()
{
		int i,j;
    double pr[64],pi[64],fr[64],fi[64];
    for (i=0; i<=63; i++)
      { pr[i]=exp(-0.1*(i+0.5)); pi[i]=0.0;}
    printf("\n");
    for (i=0; i<=15; i++)
      { for (j=0; j<=3; j++)
          printf("%e   ",pr[4*i+j]);
        printf("\n");
      }
    printf("\n");
    kfft(pr,pi,64,6,fr,fi,0,1);//���ٸ���Ҷ�任��ʱ�� ʵ��pr������pi;   Ƶ��fr+fi; ��MATLAB����һ��
    for (i=0; i<=15; i++)
      { for (j=0; j<=3; j++)
          printf("%e   ",fr[4*i+j]);
        printf("\n");
      }
    printf("\n");
    for (i=0; i<=15; i++)
      { for (j=0; j<=3; j++)
          printf("%e   ",fi[4*i+j]);
        printf("\n");
      }
    printf("\n");
    for (i=0; i<=15; i++)
      { for (j=0; j<=3; j++)
          printf("%e   ",pr[4*i+j]);
        printf("\n");
      }
    printf("\n");
    for (i=0; i<=15; i++)
      { for (j=0; j<=3; j++)
          printf("%e   ",pi[4*i+j]);
        printf("\n");
      }
    printf("\n");
    kfft(fr,fi,64,6,pr,pi,1,1);//����Ҷ���任��pr��fr������Ϊԭʱ���źţ�pi��fi����Ϊ0��
    for (i=0; i<=15; i++)
      { for (j=0; j<=3; j++)
          printf("%e   ",pr[4*i+j]);
        printf("\n");
      }
    printf("\n");



	return 0;
}

