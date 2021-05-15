#pragma once
#include <string>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include "vector.cpp"

const double inf = 1e29;

class Camera{
    public:
        Vector center;
        double fov; // field of view
        double f;   // distance screen-camera center
        double width;
        double height; 

        explicit Camera(Vector center, double alpha, int W, int H)
        {
            this->center = center;
            fov = alpha;
            width = W;
            height = H;
        }

        Vector pixel(int x, int y) {
            double f = this->width / (2 * tan(this->fov / 2));
            return  Vector(this->center[0] + x + 0.5 - this->width / 2,
                            this->center[1] + y + 0.5 - this->height / 2,
                            this->center[2] - f);
        }
        
};

class Ray{
    public:
        Vector origin;
        Vector direction;

        explicit Ray(Vector O, Vector d) {
            origin = O;
            direction = normalization(d);
        }
};

class Intersection {
    public:

        bool is_intersection = 0;
        Vector position;
        Vector normal;
        double distance;
        int index;
        
        Intersection(){}

        Intersection(bool b, Vector P, Vector N, double d, int i){
            is_intersection = b;
            position = P;
            normal = N;
            distance = d;
            index = i;        
        }
};

class Geometry{
    public:
        Vector albedo = Vector(1,1,1);
        int index;
        bool mirror = false;
        bool transparent = false;
        double refract_index;
        virtual Intersection intersect(const Ray &r) const = 0;
        void setTransparent(double ri) {
            refract_index = ri;
            transparent = true;
        }
        void setMirror() {
            mirror = true;
        }
};
