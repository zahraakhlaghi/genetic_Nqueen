#include <stdio.h>
#include <math.h>    /* for pow(x, y) */
#include <stdlib.h>  /* for calloc */
#include "type.h"

double decode(IPTR pj, int index, int size);

double binToDec(int *chrom, int l);

double eval(POPULATION *p, IPTR pj)
/* Called from gen.c and init.c */
{
    int n = p->lchrom;

    int maxFitness = (n * (n - 1)) / 2;

    double horizontal_collisions = 0;
    for (int i = 0; i < p->lchrom; i++) {
        for (int j = 0; j < p->lchrom; j++)
            if (pj->chrom[i] == pj->chrom[j])
                horizontal_collisions++;
        horizontal_collisions = horizontal_collisions - 1;
    }
    horizontal_collisions = horizontal_collisions / 2;


    int *left, *right;
    left = (int *) calloc(2 * n, sizeof(int));
    right = (int *) calloc(2 * n, sizeof(int));

    for (int i = 0; i < 2 * n; i++) {
        right[i] = 0;
        left[i] = 0;
    }
    for (int i = 0; i < n; i++) {
        left[i + pj->chrom[i] - 1] += 1;
        right[n - i + pj->chrom[i] - 2] += 1;
    }

    double diagonal_collisions = 0;

    for (int i = 0; i < (2 * n - 1); i++) {
        int counter = 0;
        if (left[i] > 1)
            counter += left[i] - 1;
        if (right[i] > 1)
            counter += right[i] - 1;
        diagonal_collisions += counter / (n - abs(i - n + 1));
    }

    return maxFitness - (horizontal_collisions + diagonal_collisions);

}
////////////////////////////////////////////////////////////////////////

double eval2(POPULATION *p, IPTR pj)
{
    int n = p->lchrom;

    int maxFitness = (n * (n - 1)) / 2;

    double horizontal_collisions = 0;
    for (int i = 0; i < p->lchrom; i++) {
        for (int j = 0; j < p->lchrom; j++)
            if (pj->chrom[i] == pj->chrom[j])
                horizontal_collisions++;
        horizontal_collisions = horizontal_collisions - 1;
    }
    horizontal_collisions = horizontal_collisions / 2;


    double diagonal_collisions = 0;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (abs(i - j) == abs(pj->chrom[i] - pj->chrom[j])){
                diagonal_collisions++;
            }
        }
        diagonal_collisions--;
    }

    diagonal_collisions/=2;

    return maxFitness - (horizontal_collisions + diagonal_collisions);


}
//////////////////////////////////////////////////////////////////////////

int one_queen_penalty(int index, IPTR p, int pSize) {

    int col = index;
    int row = p->chrom[index];
    int penalty = 0;
    for (int i = 0; i < pSize; i++) {
        if (i == col)
            continue;
        if (p->chrom[i] < row && p->chrom[i] + abs(col - i) == row)
            penalty++;
        else if (p->chrom[i] > row && p->chrom[i] - abs(col - i) == row)
            penalty++;
    }
    return penalty;
}

int configuration_penalty(int size, IPTR pj) {

    int sum = 0;
    for (int i = 0; i < size; i++) {
        sum += one_queen_penalty(i, pj, size);
    }
    return sum;
}

double eval3(POPULATION *p, IPTR pj) {

    int penalty = configuration_penalty(p->lchrom, pj);

    if (penalty > 0)
        return (double)(1 / (double)penalty);
    else
        return 2;

}

///////////////////////////////////////////////////////////////////////////
double decode(IPTR pj, int index, int size) {
    return ((double) binToDec(&(pj->chrom[0]), size));
}

double binToDec(int *chrom, int l) {
    int i;
    double prod;

    prod = 0.0;

    for (i = 0; i < l; i++) {
        prod += (chrom[i] == 0 ? 0.0 : pow((double) 2.0, (double) i));
    }
    return prod;
}

void decToBin(int ad, int *barray, int size) {
    int i, t;

    t = ad;
    for (i = 0; i < size; i++) {
        barray[i] = t % 2;
        t = t / 2;
    }
}
