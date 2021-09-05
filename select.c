
#include "type.h"

double f_random();

int rnd(int low, int high);

int roulette(IPTR pop, double sumFitness, int popsize) {

    /* select a single individual by roulette wheel selection */
    //یک عدد تصادفی تولید میکند سپس فیتنس ها را یکی یکی اضافه میکن تا عددتصافی تولید شده در بازه فیتنس ان قرار گیرد
    double rand, partsum;
    int i;
    partsum = 0.0;
    i = 0;
    rand = f_random() * sumFitness;
    i = -1;
    do {
        i++;
        partsum += pop[i].fitness;
    } while (partsum < rand && i < popsize - 1);

    return i;
}

int tournament_selection(IPTR p, int pSize) {

    int k1 = 0, k2 = 0;

    do {
        k1 = rnd(0, pSize);
        k2 = rnd(0, pSize);

    } while (k1 == k2);

    return p[k1].fitness > p[k2].fitness ? k1 : k2;
}

