#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"
#include <omp.h>


#include <chrono>

#include "mesh_loader.cpp"
#include "scene.cpp"

//test g++ -O3 -std=c++11 main.cpp -o a.out 
//./a.out


// #include "simple_obj_file_reader.hpp"

using namespace std::chrono;


int main(int argc, char **argv){
    auto start = high_resolution_clock::now();

    // Sphere s_right_1 = Sphere(Vector(21, 0, 0), 9, Vector(1, 1, 1), "transparent",1.5,true);
    Geometry* s_left = new Sphere(Vector(-21,0,0), 10, Vector(0,1,1));
    s_left->setMirror();
    Geometry* s_middle = new Sphere(Vector(0, 0, 0), 10, Vector(1, 1, 1));
    s_middle->setTransparent(1.5);
    Geometry* s_right_1 = new Sphere(Vector(21, 0, 0), 10, Vector(1, 0, 1)); // , "transparent",1.5);
    Geometry* s_right_2 = new Sphere(Vector(21, 0, 0), 9.5, Vector(1, 0, 1)); // , "transparent",1.5);
    s_right_1->setTransparent(1.5);
    s_right_2->setTransparent(1./1.5);
    Geometry* s_green = new Sphere(Vector(0, 0, -1000), 940, Vector(0, 1, 0));
    Geometry* s_blue = new Sphere(Vector(0, -1000, 0), 990, Vector(0, 0, 1));
    Geometry* s_magenta = new Sphere(Vector(0, 0, 1000), 940, Vector(1,0,1));
    Geometry* s_red = new Sphere(Vector(0, 1000, 0), 940,Vector(1, 0, 0));
    Geometry* s_cyan = new Sphere(Vector(-1000, 0, 0), 940, Vector(0, 1, 1));
    Geometry* s_yellow = new Sphere(Vector(1000, 0, 0), 940, Vector(1, 1, 0));
    
    Geometry* cat = new TriangleMesh();
    ((TriangleMesh *)cat)->readOBJ("cat_model/cat.obj");
	((TriangleMesh *)cat)->transform(Vector(0,-10,0), .6);
    ((TriangleMesh *)cat)->init();


    const double pi = acos(-1);
    Vector center = Vector(0,0,55); 
    double alpha = 60. / 180 * pi;                  
    const int W = 520;
    const int H = 520;
    Camera cam = Camera(center,alpha,W,H);
    Light L = Light(pow(10.,5.),Vector(-10,20,40));
    int max_path_length = 5;
    int samples = 32;

    Scene scene;
    // scene.add_object(s_left);
    // scene.add_object(s_middle);
    // scene.add_object(s_right_1);
    // scene.add_object(s_right_2);
    scene.add_object(s_blue);
    scene.add_object(s_magenta);
    scene.add_object(s_green);
    scene.add_object(s_red);
    scene.add_object(s_cyan);
    scene.add_object(s_yellow);
    scene.add_object(cat);
    scene.set_light(L);

    unsigned char image[W * H * 3];

    #pragma omp parallel for schedule(static,1)
    for (int i = 0; i < H; i++){
        for (int j = 0; j < W; j++){
            bool muller = false;         // if using muller box for antialiasing

            // std::vector<double> index = {1};
            // color = scene.getColor(r, max_path_length);
            // std::cout << *index.size() << std::endl;

            Vector cumul = Vector(0,0,0);
            for (int u = 0; u < samples; u++){
                auto direction = cam.pixel(j, H - i - 1) - center;
                double x,y; boxMuller(0.01,x,y);
                direction[0] += x;
                direction[1] += y;
                Ray r = Ray(center, direction);

                Vector color = scene.getColor(r,max_path_length);
                // std::cout << std::endl;
                // color = scene.getColor(r,max_path_length);
                // print(color);
                // std::cout << color[0] << ", " << color[1] << ", " << color[2] << std::endl;
                cumul = cumul + color;
            }
            
            Vector color = cumul/samples;
            // print(color);

            double power = 1. / 2.2;
            image[(i * W + j) * 3 + 0] = std::min(255., std::max(0., pow(color[0], power) * 255));
            image[(i * W + j) * 3 + 1] = std::min(255., std::max(0., pow(color[1], power) * 255));
            image[(i * W + j) * 3 + 2] = std::min(255., std::max(0., pow(color[2], power) * 255));
        }
    }

    stbi_write_jpg("image.jpg", W, H, 3, image, 0);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    duration = duration / 1000000;
    std::cout << "Time taken: " << duration.count() << " seconds" << std::endl;
    return 0;
}