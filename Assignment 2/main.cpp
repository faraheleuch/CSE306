#include <iostream>
#include "objects.cpp"
#include "svg.cpp"
#include "vector.cpp"
#include <random>


std::random_device rd;
std::default_random_engine eng(rd());
std::uniform_real_distribution<float> distr(0, 1);

double generate() {
    return distr(eng);
}


int main(int argc, char **argv){

    // part 1
    /*

    const int n = 100;
    std::vector<Vector> points;
    for (int i=0; i<n; i++) {
        points.push_back(Vector(generate(), generate()));
    }
    // std::vector<Vector> points = {Vector(0.1, 0.9), Vector(0.5, 1), Vector(0.7, 0.8), Vector(0.6, 0.55), Vector(0.9, 0.4), Vector(0.8, 0.2), Vector(0.3, 0.2)};
    //const int n = points.size();

    std::vector<double> lambda(n); // tu peux mettre autre chose stv
    for(int i = 0; i<points.size();i++) {
        Vector C = Vector(0.5,0.5);
        lambda[i]= exp(-norm(points[i]-C)/0.02);
    }
    
    auto w = LBFGS(points, lambda);
    auto res = compute_vor(points, w);
    save_svg(res, "image.svg");
    */


    int n_water = 50;
    int n_air = 200;
    int n = n_water + n_air;
    std::vector<Vector> points;
    std::vector<double> mass;
    for (int i=0; i<n_water; i++) {
        points.push_back(Vector(generate(), generate()));
        mass.push_back(200);
    }
    for (int i=0; i<n_air; i++) {
        points.push_back(Vector(generate(), generate()));
        mass.push_back(100);
    }
    std::vector<Vector> velocity(n);
    std::vector<double> lambda(n);
    
    int frames = 30;
    for (int f=0; f<frames; f++) {
        std::cout << f << std::endl;
        auto [x,v] = Galouet_Merigot(points, velocity, mass);
        points = x;
        velocity = v;

        auto weights = LBFGS(points, lambda);
        auto vor = compute_vor(points, weights);
        save_svg_animated(vor, "water.svg", f, frames, n_water);
    }
   

}