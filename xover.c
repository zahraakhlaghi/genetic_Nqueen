#include <stdio.h>
#include <stdlib.h>  /* for calloc */
#include "type.h"

enum xoverType { order1, PMX };
int lchrom;


int rnd(int low, int high);
int flip(double);
int muteX(POPULATION *p, int pa, int pb);

void selectTwoRandomPoints(int *fromIndex,int *toindex)
{
    int fIndex1 = rnd(0,lchrom-1);
    int tIndex1 =  rnd(0,lchrom-1);
    while (abs(fIndex1 - tIndex1)<=1)
    {
        fIndex1 = rnd(0,lchrom-1);
        tIndex1 =  rnd(0,lchrom-1);
    }

    int fIndex = (tIndex1 < fIndex1)?tIndex1:fIndex1;
    int tIndex = (fIndex1 > tIndex1)?fIndex1:tIndex1;

    *fromIndex = fIndex;
    *toindex= tIndex;
}
///////////////////////////////////////////////////////////////////////////////////
void swap_mutation(POPULATION *pop,IPTR ind)
{
    if (flip(pop->pMut))
    {
        lchrom = pop->lchrom;
        int fromIndex = 0, toIndex = lchrom;
        selectTwoRandomPoints(&fromIndex, &toIndex);
        int tmp = (ind->chrom)[fromIndex];
        (ind->chrom)[fromIndex] = (ind->chrom)[toIndex];
        (ind->chrom)[toIndex] = tmp;
    }
}

void insert_mutation(POPULATION *pop,IPTR ind)
{
    lchrom=pop->lchrom;
    int fromIndex = 0, toIndex = lchrom;
    selectTwoRandomPoints(&fromIndex, &toIndex);
    int tmp = (ind->chrom)[toIndex];
    for(int i= toIndex -1;i> fromIndex;i--)
        (ind->chrom)[i+1]=(ind->chrom)[i];
    (ind->chrom)[fromIndex+1] = tmp;
}

void inversion_mutation(POPULATION *pop,IPTR ind)
{
    lchrom=pop->lchrom;
    int fromIndex = 0, toIndex = lchrom;
    selectTwoRandomPoints(&fromIndex, &toIndex);
    int i = fromIndex,j=toIndex;
    while(j-i > 0)
    {
        int tmp=(ind->chrom)[i];
        (ind->chrom)[i]=(ind->chrom)[j];
        (ind->chrom)[j]=tmp;
        j = j - 1;
        i = i + 1;
    }
}
//////////////////////////////////////////////////////////////////////////////////

void crossover(POPULATION *p, IPTR p1, IPTR p2, IPTR c1, IPTR c2)
{
/* p1,p2,c1,c2,m1,m2,mc1,mc2 */
  int *pi1,*pi2,*ci1,*ci2;
  int xp, i;

  pi1 = p1->chrom;
  pi2 = p2->chrom;
  ci1 = c1->chrom;
  ci2 = c2->chrom;

  if(flip(p->pCross)){

    xp = rnd(0, p->lchrom - 1);
    for(i = 0; i < xp; i++){
      ci1[i] = muteX(p, pi1[i],pi2[i]);
      ci2[i] = muteX(p, pi2[i],pi1[i]);
    }
    for(i = xp; i < p->lchrom; i++){
      ci1[i] = muteX(p, pi2[i],pi1[i]);
      ci2[i] = muteX(p, pi1[i],pi2[i]);
    }
  } else {
    for(i = 0; i < p->lchrom; i++){
      ci1[i] = muteX(p, pi1[i],pi2[i]);
      ci2[i] = muteX(p, pi2[i],pi1[i]);
    }
  }
}

void copyItems(IPTR fromInd , int fromIndex , int toIndex , IPTR toInd)
{
    for(int i = fromIndex ; i < toIndex;i++)
    {
        (toInd->chrom)[i] = (fromInd->chrom)[i];
    }

}

int isExist(int *chrom , int fromIndex , int toIndex,int key)
{
    for(int i = fromIndex;i<toIndex;i++)
        if(chrom[i]==key)
            return i;
    return (0==1);
}

void fillOtherChromItems(enum xoverType type, int fromIndex, int toIndex, IPTR c, IPTR p)
{
    int j=0,k=0;
    switch (type)
    {
        //پر کردن سایر کروموزوم ها برای روش order crossover
        case order1:
            j = toIndex, k = toIndex;
            while (j < (lchrom) && (j >= toIndex || j <= fromIndex))
            {
                if (!isExist(c->chrom, fromIndex + 1, toIndex, (p->chrom)[k]))
                {
                    (c->chrom)[j] = (p->chrom)[k];
                    j++;
                    k++;
                }
                else
                {
                    k++;
                }
                if (k == lchrom)
                    k = 0;
                if (j == lchrom)
                    j = 0;
            }

            break;

            //پر کردن سایر کروموزوم ها برای روش PMX crossover
        case PMX:
            j = fromIndex;
            for (int i = fromIndex ; i<toIndex;i++)
            {
                if (!isExist((c->chrom), fromIndex , toIndex, (p->chrom)[i]))
                {
                    int l = 0;
                    j = isExist((p->chrom), 0, lchrom, (c->chrom)[i]);
                    while (j>=fromIndex && j<toIndex && l<lchrom)
                    {

                        j = isExist((p->chrom), 0, lchrom, (c->chrom)[j]);
                        l++;

                    }
                    (c->chrom)[j] = (p->chrom)[i];
                }

            }
            for (int i = 0;i<lchrom;i++)
                if ((c->chrom)[i] == -1)
                    (c->chrom)[i] = (p->chrom)[i];
            break;
    }


}

void order1_crossover(POPULATION *pop,IPTR p1,IPTR p2 ,IPTR c1,IPTR c2)
{
    for (int i = 0; i < lchrom;i++)
    {
        (c1->chrom)[i] = -1;
        (c2->chrom)[i] = -1;
    }
    lchrom = pop->lchrom;
    // با بررسی شرط احتمال crossover
    if (flip(pop->pCross))
    {

        //دو نقطه در کروموزوم به صورت تصادفی انتخاب می شود
        int lPoint=0, rPoint= lchrom;
        selectTwoRandomPoints(&lPoint, &rPoint);

        //مقادیر بین دو نقطه در فرزندان کپی می شوند.
        copyItems(p1, lPoint + 1, rPoint, c1);
        copyItems(p2, lPoint + 1, rPoint, c2);

        //پر کردن سایر بخش های باقی مانده در کروموزوم
        fillOtherChromItems(order1,lPoint, rPoint, c1, p2);
        fillOtherChromItems(order1,lPoint, rPoint, c2, p1);
    }
    else
    {
        copyItems(p1, 0, lchrom, c1);
        copyItems(p2, 0, lchrom, c2);
    }

}

void PMX_crossover(POPULATION *pop,IPTR p1,IPTR p2 ,IPTR c1,IPTR c2)
{
    lchrom=pop->lchrom;

    //انتخاب دو نقطه تصادفی
    int lPoint = 0, rPoint = lchrom;
    selectTwoRandomPoints(&lPoint, &rPoint);
    for (int i = 0; i < lchrom;i++)
    {
        (c1->chrom)[i] = -1;
        (c2->chrom)[i] = -1;
    }
    //کپی مقادیر بین دو نقطه در فرزندان
    copyItems(p1,lPoint+1,rPoint,c1);
    copyItems(p2,lPoint+1,rPoint,c2);

    int c1_after[100];
    int c2_after[100];
    for (int i = 0;i < pop->lchrom;i++)
    {
        c1_after[i] = (c1->chrom)[i];
        c2_after[i] = (c2->chrom)[i];
    }

    //پر کردن کروموزوم ها
    fillOtherChromItems(PMX,lPoint+1,rPoint,c1, p2);
    fillOtherChromItems(PMX,lPoint+1,rPoint,c2, p1);


    int c1_com[100];
    int c2_com[100];
    for (int i = 0;i < pop->lchrom;i++)
    {
        c1_com[i] = (c1->chrom)[i];
        c2_com[i] = (c2->chrom)[i];
    }
}

int muteX(POPULATION *p, int pa, int pb)
{
  return (flip(p->pMut) ? 1 - pa  : pa);
}



