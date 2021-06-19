#include "vector_3d.cpp"
#include <random>
#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"


std::random_device rd;
std::default_random_engine eng(rd());
std::uniform_real_distribution<float> distr(0, 1);

Vector sample_sphere(){
    double r1 = distr(eng);
    double r2 = distr(eng) ;
    double x = cos(2*M_PI*r1)*sqrt(1. - r2);
    double y = sin(2*M_PI*r1)*sqrt(1. - r2);
    double z = sqrt(r2);
    if (distr(eng) < .5) z = -z;
    return Vector(x,y,z);
}

void Color_transfer(std::vector<Vector>& I, std::vector<Vector>& M){
    int n = I.size();
    int nbiter = 30;
    for(int iter = 0; iter < nbiter; iter++){
        Vector v = sample_sphere();
        std::vector<std::pair<double, int> > projI(n);
        std::vector<std::pair<double, int> > projM(n);
        for(int i = 0; i<n; i++ ){
            projI[i] = {dot(I[i],v),i};
            projM[i] = {dot(M[i],v),i};
        }
        std::sort(projI.begin(), projI.end());
        std::sort(projM.begin(), projM.end());

        for(int i = 0; i<n; i++ ){
            int j = projI[i].second;
            I[j] += (projM[i].first - projI[i].first)*v;
        }
    }
}

int main() {
    int x,y,n;
    unsigned char *im1 = stbi_load("gabriel.jpeg", &x, &y, &n, 3);
    unsigned char *im2 = stbi_load("axes.jpeg", &x, &y, &n, 3);

    std::vector<Vector> I(x*y), M(x*y);
    for (int i = 0; i < x*y; i++) {
        I[i] = Vector(im1[3*i+0],im1[3*i+1],im1[3*i+2]);
        M[i] = Vector(im2[3*i+0],im2[3*i+1],im2[3*i+2]);
    }
    Color_transfer(I, M);
    for (int i = 0; i < x*y; i++) {
        im1[3*i+0] = I[i][0];
        im1[3*i+1] = I[i][1];
        im1[3*i+2] = I[i][2];
    }
    stbi_write_jpg("out.jpg", x, y, n, im1, 80);
}