#include "vector.cpp"
#include <random>


std::random_device rd;
std::default_random_engine eng(rd());
std::uniform_real_distribution<float> distr(0, 1);

double generate() {
    return distr(eng);
}

void boxMuller (double stdev , double &x , double &y ) {
    double r1 = distr(eng);
    double r2 = distr(eng) ;
    x = sqrt(-2 * log (r1)) *cos(2*M_PI*r2) *stdev ;
    y = sqrt(-2 * log (r1)) *sin(2*M_PI*r2) *stdev ;
}

Vector hemisphere(){
    double r1 = distr(eng);
    double r2 = distr(eng) ;
    double x = cos(2*M_PI*r1)*sqrt(1. - r2);
    double y = sin(2*M_PI*r1)*sqrt(1. - r2);
    double z = sqrt(r2);
    return Vector(x,y,z);
}

Vector random_cos(const Vector& N) {
    double bv = fabs(N[0]);
    int bi = 0;
    if (N[1] < bv) bv = fabs(N[1]), bi = 1;
    if (N[2] < bv) bv = fabs(N[2]), bi = 2;
    
    auto T1 = N;
    T1[bi] = 0;
    int a = bi == 0;
    int b = 1 + (bi != 2);

    T1[b] = -T1[b];
    T1 = normalization(T1);

    Vector T2 = cross(N, T1);

    Vector u = hemisphere();

    return u[0] * T1 + u[1] * T2 + u[2] * N;
}

