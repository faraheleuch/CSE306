#pragma once
#include "vector.cpp"
#include <vector>
#include <math.h>

#include <stdio.h>


double inf = 1e88;

class Polygon{
    public:
        std::vector<Vector> vertices;

        Polygon(){};
        Polygon(std::vector<Vector>& V){
        vertices = V;
    }

    Vector Centroid(){
        Vector res(0,0);
        const int N = vertices.size();
        for(int i = 0; i< N; i++){
            res += vertices[i];
        }
        return res / N;
    }
};

 void print(Polygon &a){
    std::cout << "[" ;
    for (auto bubbles : a.vertices) {
        std::cout << bubbles << ',';
    };
    std::cout  << "]," << std::endl;
}

Vector intersect(Vector A, Vector B, std::pair<Vector,Vector> line){
    // returns the point of intersection between the Edge [A,B] and line (u,v)
    Vector u = line.first;
    Vector v = line.second;
    Vector N = Vector(v[1]-u[1],u[0]-v[0]);     
    double t = dot(u-A,N)/dot(B-A,N);
    if (t < 0 || t > 1){    // no intersection
        return Vector(inf,inf);
    }
    return A + t*(B-A);
};

bool inside(Vector P, std::pair<Vector,Vector> clipEdge){
    Vector u = clipEdge.first;
    Vector v = clipEdge.second;
    Vector N = Vector(v[1] - u[1], u[0] - v[0]); //outwards normal to clipEdge
    if (dot(P-u,N) <= 0) return true ;
    return false;
}

Polygon clip_poly(Polygon subjectPolygon, Polygon clipPolygon){
    Polygon outPolygon;
    for(int i = 0; i < clipPolygon.vertices.size(); i++) {
        std::pair<Vector, Vector> clipEdge(clipPolygon.vertices[i], clipPolygon.vertices[(i != 0) ? i - 1 : clipPolygon.vertices.size() - 1]);
        outPolygon = Polygon();
        for (size_t i = 0; i < subjectPolygon.vertices.size(); i++){
            Vector curVertex = subjectPolygon.vertices[i];
            Vector prevVertex = subjectPolygon.vertices[(i>0)?(i-1):subjectPolygon.vertices.size()-1];
            Vector intersection = intersect(prevVertex,curVertex,clipEdge);
            if (inside(curVertex,clipEdge)){
                if (!inside(prevVertex,clipEdge)) outPolygon.vertices.push_back(intersection);
                outPolygon.vertices.push_back(curVertex);
            }
            else if (inside(prevVertex,clipEdge)) outPolygon.vertices.push_back(intersection);
        }
        subjectPolygon = outPolygon;
    }
    return outPolygon;
};

class Triangle {
    public:
        std::vector<Vector> vtx;
    explicit Triangle(std::vector<Vector> vertices = {}){
        vtx = vertices;
    }
};

bool in_circle(Vector Pi, Triangle ABC) {
    Vector u = ABC.vtx[1] - ABC.vtx[0];
    Vector v = ABC.vtx[2] - ABC.vtx[0];
    Vector M = (ABC.vtx[0] + ABC.vtx[1])/2;
    Vector N = (ABC.vtx[0] + ABC.vtx[2])/2;
    Vector up = Vector(u[1], -u[0]);
    Vector vp = Vector(v[1], -v[0]);

    Vector K = (dot(u,M) * vp - up * dot(v,N))/(u[0]*v[1]-u[1]*v[0]);
    double r = abs(sqrt((pow(K[0]-ABC.vtx[0][0],2))+(pow(K[1]-ABC.vtx[0][1],2))));

    if (abs(sqrt(pow(Pi[0]-K[0],2)+(pow(Pi[1]-K[1],2)))) <= r) return true ;
    return false;
}

std::vector<Polygon> compute_vor(std::vector<Vector>& pts, std::vector<double>& w) {
    if (w.size() == 0) {
        w = std::vector<double>(pts.size());
    }
    std::vector<Polygon> ricky;
    for (int i=0; i<pts.size(); i++) {
        std::vector<Vector> tmp = {Vector(1,0),Vector(1,1),Vector(0,1), Vector(0,0)};
        Polygon pol_i = Polygon(tmp);
        for (int j=0; j<pts.size(); j++) {
            if (j == i) continue;
            Vector Pi = pts[i];
            Vector Pj = pts[j];
            // compute mediatrice pi <> pj
            Vector o = (Pi+Pj)/2 + (w[i]-w[j])/(2*pow(norm(Pj-Pi),2)) * (Pj-Pi);
            Vector d = normalization(Pi-Pj); // points toward i
            Vector N = Vector(d[1], -d[0]);

            // create big square
            Vector A = o + N * 10;
            Vector D = o - N * 10;
            Vector B = A + d *20 ;
            Vector C = D + d *20;

            tmp = {A,D,C,B};
            Polygon square = Polygon(tmp);
            //print(pol_i);
            //print(square);

            pol_i = clip_poly(pol_i, square);
            //print(pol_i);
        }
        ricky.push_back(pol_i);        
    }
    return ricky;
    
}

double integrale(Vector P, Polygon& pol) {
    int N = pol.vertices.size();
    
    double res = 0 ;
    for(int k = 0; k<N; k++){
        Vector a = pol.vertices[k] - P;
        Vector b = pol.vertices[(k+1)%N] - P;
        res += cross(a,b) * dot(a,b) / 6;
    }
    /*
    std::vector<double> x(N), y(N);

    for(int i = 0; i<N; i++){
        x[i] = pol.vertices[i][0];
        y[i] = pol.vertices[i][1];
    }
    double res = 0 ;
    for(int k = 1; k<=N; k++){
        int kn = k % N;
        res += (x[k-1] * y[kn] - x[kn]*y[k-1])*(pow(x[k-1],2) + x[k-1]*x[kn] + pow(x[kn],2) + pow(y[k-1],2) + y[k-1]*y[kn] + pow(y[kn],2)
        - 4*(P[0]*(x[k-1]+ x[kn])+ P[1]*(y[k-1]+y[kn])) + 6*dot(P,P));
    }
    */

    return res;
}


double area(const Polygon& p) {
    double res = 0;
    Vector o = p.vertices[0];
    for (int i=1; i+1<p.vertices.size(); i++) {
        Vector a = p.vertices[i] - o;
        Vector b = p.vertices[i+1] - o;
        res += cross(a,b);
    }
    return res / 2;
}
/*

#include "lbfgs.h"
static lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
    )
{
    auto [pts, lambda] = *(std::pair<std::vector<Vector>, std::vector<double> >*) instance;
    std::vector<double> w(n);
    for (int i=0; i<n; i++) w[i] = x[i];
    auto pols = compute_vor(pts, w);

    lbfgsfloatval_t fx = 0.0;
    for (int i=0; i<n; i++) {
        double A = area(pols[i]);
        fx -= - integrale(pts[i], pols[i]) - A*w[i] + lambda[i]*w[i];
        g[i] = - (-A + lambda[i]);
    }
    return fx;
}
static int progress(
    void *instance,
    const lbfgsfloatval_t *x,
    const lbfgsfloatval_t *g,
    const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm,
    const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step,
    int n,
    int k,
    int ls
    )
{
    printf("Iteration %d:\n", k);
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");
    return 0;
}


std::vector<double> LBFGS(std::vector<Vector>& pts, std::vector<double>& lambda)  {
    int n = pts.size();

    int i, ret = 0;
    lbfgsfloatval_t fx;
    lbfgsfloatval_t *x = lbfgs_malloc(n); // weights
    lbfgs_parameter_t param;

    for (i = 0; i < n; i++) {
        x[i] = 0;
    }

    std::pair<std::vector<Vector>, std::vector<double> > args = {pts, lambda};
    lbfgs(n, x, &fx, evaluate, progress, &args, &param);

    std::vector<double> w (n);
    for (int i=0; i<n; i++) w[i] = x[i];
    return w;
}
*/

std::vector<double> LBFGS(std::vector<Vector>& pts, std::vector<double>& lambda)  {
    int n = pts.size();
    std::vector<double> w (n);

    double eps = .1;
    for(int it = 0; it<20; it++) {
        double loss = 0;
        auto pols = compute_vor(pts, w);
        double A = 0;
        for(int i = 0; i < n; i++){
            double a = area(pols[i]);
            A += a;
            double dGi = - a + lambda[i];
            w[i] += eps * dGi;
            loss += abs(dGi);
        }
        // std::cout << A << "loss:" << loss << std::endl;
    }
    // std::cout << "end" << std::endl;
    
    return w;
}

std::pair< std::vector<Vector>, std::vector<Vector> > Galouet_Merigot(
    std::vector<Vector>& X, std::vector<Vector>& v, std::vector<double>& m
){
    const int N = X.size();
    std::vector<double> Uniform(N,1./N);
    std::vector<double> Vw = LBFGS(X,Uniform);
    std::vector<Polygon> Pol = compute_vor(X, Vw) ;
    std::vector<Vector> vnext(N); 
    std::vector<Vector> Xnext(N); 
    double eps = 0.004;
    double dt = 0.002;
    for(int i = 0; i<N; i++) {
        Vector C;
        if (Pol[i].vertices.size()) C = Pol[i].Centroid() ;
        else C = Vector(
            std::min(1.,std::max(0.,X[i][0])),
            std::min(1.,std::max(0.,X[i][1]))
        );

        Vector Fspring = 1./pow(eps,2) * (C - X[i]);
        Vector F = Fspring + Vector(0, -9.81);
        vnext[i] = v[i] + dt*F/m[i];
        Xnext[i] = X[i] + dt*v[i];
        //std::cout << Xnext[i] << Pol[i].vertices.size();
    }
    //std::cout << std::endl << std::endl;
    return {Xnext, vnext};
}