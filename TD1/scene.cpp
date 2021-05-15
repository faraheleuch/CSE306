#pragma once
#include <string>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include "vector.cpp"
#include "sphere.cpp"
#include "objects.cpp"
#include "random.cpp"

const double pi = acos(-1);

class Light {
    public :
    double I;
    Vector S;
    Light() {}
    Light(double Intensity, Vector Position) {
        I = Intensity;
        S = Position;
    }
};



class Scene {
    public:
    std::vector<Geometry*> elements ;
    Light light;

    Scene() {}
    void add_object(Geometry* O){
        O->index = elements.size();
        this->elements.push_back(O);
    }
    
    void set_light(Light l) {
        light = l;
    }


    Intersection intersection(const Ray& r) const {
        double closest = inf;
        Intersection res;
        for (int i = 0; i < elements.size(); i++){
            Intersection inter = elements[i]->intersect(r);
            if (inter.is_intersection){
                // std::cout << inter.distance << ' ';
                if (inter.distance < closest && inter.distance > 1e-5){
                    closest = inter.distance;
                    res = inter;
                }
            }
        };
        return res;
    }

    Vector Lambertian(Intersection inter) const {
        auto I = light.I;
        auto P = inter.position;
        auto albedo = elements[inter.index]->albedo;
        auto N = inter.normal;
        auto d = norm(light.S - P);
        auto w_i = (light.S - P) / d;
        double visibility = 0;
        Ray r (P, w_i);
        auto x = this->intersection(r); 
        if (!x.is_intersection || x.distance > d) {
            visibility = 1;
        }
        // visibility = 1;
        return (albedo/pi) * I / (4*pi*pow(d,2)) * visibility * dot(N,w_i);
        
    }

    Vector getColor (const Ray& ray , int ray_depth) const {
        if (ray_depth < 0 ) return Vector ( 0. , 0. , 0. ) ;
        Intersection x = this->intersection(ray); 
        if (x.is_intersection ) { 
            //
            Geometry* g = elements[x.index];
            if (g->mirror ) {
                // std::cout << "mirror" << std::endl;
                Ray reflected_ray = Ray(x.position, ray.direction - x.normal* 2*dot(ray.direction,x.normal)) ;
                return getColor(reflected_ray , ray_depth - 1) ;
            }else if (g->transparent) {
                // std::cout << "transp" << std::endl;
                // refraction occurs
                double n1 = 1;
                double n2 = g->refract_index;
                Vector w_i = ray.direction;

                // are we going out of the sphere?
                bool invert = dot(ray.direction, x.normal) > 0;

                if (invert) {
                    std::swap(n1,n2);
                }
                if (x.index == 0 || x.index == 1) {
                    // std::cout << x.index << x.position << std::endl;
                    //std::cout << (ray.direction[2] > 0) << " ";
                }

                double k_0 = pow(n1-n2,2)/pow(n1+n2,2);
                double R = k_0 + (1-k_0) * pow(1. - fabs(dot(x.normal,w_i)),5);
                double T = 1 - R;

                double u = generate();
                //std::cout << n1 << ' ' << n2 << ' ' << R << std::endl;
                if(u<R) {
                    Ray reflected_ray = Ray(x.position, ray.direction - x.normal* 2*dot(ray.direction,x.normal)) ;
                    return getColor(reflected_ray , ray_depth - 1) ;
                }

                
                Vector w_tT = (w_i - x.normal * dot(w_i,x.normal)) * n1/n2;
                Vector w_tN =  -x.normal * sqrt(1-pow(n1/n2,2) * (1-pow(dot(ray.direction,x.normal),2)));
                if (invert) w_tN = -w_tN;
                Ray refracted_ray = Ray(x.position, w_tT+w_tN);
                auto col = getColor(refracted_ray, ray_depth - 1);
                // std::cout << col << std::endl;
                return col;
            } else {/// diffuse    
                //return Lambertian(x);
                Vector Lo = Lambertian(x);
                Ray randomRay(x.position, random_cos(x.normal));
                Lo += elements[x.index]->albedo * getColor(randomRay, ray_depth -1);
                return Lo;

            } 
        }
        //std::cout << "no" << std::endl;
        return Vector(0,0,0);
    }
};
