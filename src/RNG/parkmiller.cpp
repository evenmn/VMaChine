#include "parkmiller.h"
#include "rng.h"
#include <cmath>

ParkMiller::ParkMiller()  :
    RandomNumberGenerator() {
}

long     ParkMiller::iy = 0;
long     ParkMiller::iv[NTAB];
long     ParkMiller::seed = -1;
void ParkMiller::setSeed(long seed) {
    ParkMiller::seed = seed;
}

double ParkMiller::nextGaussian(double mean, double standardDeviation) {
    double standardNormalRandomNumber = sqrt( -2.0*log(1.0 - nextDouble()) ) * cos( 6.283185307 * nextDouble() );
    return standardDeviation*standardNormalRandomNumber + mean;
}

int ParkMiller::nextInt(int upperLimit) {
    return int(std::floor(nextDouble() * upperLimit));
}

double ParkMiller::nextDouble()
{
    int             j;
    long            k;
    double          temp;
    if (ParkMiller::seed <= 0 || !iy) {
        if (-(ParkMiller::seed) < 1) ParkMiller::seed=1;
        else ParkMiller::seed = -(ParkMiller::seed);
        for(j = NTAB + 7; j >= 0; j--) {
            k     = (ParkMiller::seed)/IQ;
            ParkMiller::seed = IA*(ParkMiller::seed - k*IQ) - IR*k;
            if(ParkMiller::seed < 0) ParkMiller::seed += IM;
            if(j < NTAB) iv[j] = ParkMiller::seed;
        }
        iy = iv[0];
    }
    k     = (ParkMiller::seed)/IQ;
    ParkMiller::seed = IA*(ParkMiller::seed - k*IQ) - IR*k;
    if(ParkMiller::seed < 0) ParkMiller::seed += IM;
    j     = int(iy/NDIV);
    iy    = iv[j];
    iv[j] = ParkMiller::seed;
    if((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}

