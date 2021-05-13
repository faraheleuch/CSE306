#include <string>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include "vector.cpp"
#include "objects.cpp"


class Sphere : public Geometry {
    public:
        Vector center;
        double radius;

        explicit Sphere(Vector C = Vector(0,0,0), double R = 0, Vector A = Vector(0,0,0)) {
            center = C;
            radius = R;
            albedo = A;
        }

        Intersection intersect(const Ray &r) const {
            double t = 0;
            bool is_inter = true;
            Vector u = r.direction;
            Vector O = r.origin;
            Vector C = this->center;
            double R = this->radius;
            Vector OC = C-O;
            auto delta = pow(dot(u, OC), 2) - (pow(norm(OC),2) - pow(R, 2));
            
            if (delta < 0){         // no solutions
                is_inter = false;
                return Intersection();
            } else {
                auto t2 = dot(u, OC) + sqrt(delta);
                if (t2 < 0) { 
                    is_inter = false; // no intersection
                    return Intersection();
                }
            
                auto t1 = dot(u, OC) - sqrt(delta);

                if (t1 > 1e-5){
                    t = t1;
                } else {
                    t = t2;
                }

                Vector inter_point = O + u*t; //intersection point
                Vector normal = (inter_point - C) / norm(inter_point - C); //unit normal N at P 
                return Intersection(is_inter, inter_point, normal, t, index);
            }
        }

};