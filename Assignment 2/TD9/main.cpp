#include <vector>
#include <map>
#include "mesh_loader.cpp"

std::vector<Vector> Tutte(const TriangleMesh& mesh, int nb_iter=5) {

    int n0 = mesh.vertices.size();
    std::vector<std::vector<int> > nbs(n0);
    std::vector<bool> on_boundary(n0);
    int n;
    std::vector<int> b;
    {
        /* Find the boundary {bi}
        assuming the triangles are oriented
        */
        std::map<std::pair<int, int>, int> cnt;
        for (TriangleIndices tri: mesh.indices) {
            int idx[] = {tri.vtxi, tri.vtxj, tri.vtxk};
            for (int i=0; i<3; i++) {
                int a = idx[i];
                int b = idx[(i+1) % 3];
                cnt[{a,b}]++;
                cnt[{b,a}]--;
                nbs[a].push_back(b);
            }
        }
        std::map<int, int> dM;
        for (auto [ab,k]: cnt) if (k==1) {
            auto [a,b] = ab;
            nbs[b].push_back(a);
            dM[a] = b;
        }

        n = dM.size();
        for (int p = dM.begin()->first, i=0; i < n; p = dM[p], i++) {
            on_boundary[p] = true;
            b.push_back(p);
        }
    }
    // std::cout << mesh.indices.size() << std::endl;
    // std::cout << b.size() << std::endl;

    b.push_back(b[0]); // bn = b0
    double s = 0;
    for(int i = 0; i<n;i++) {
        s += norm(mesh.vertices[b[i+1]]-mesh.vertices[b[i]]);
    }
    std::vector<std::vector<Vector> > v(nb_iter+1);
    v[0] = mesh.vertices;

    double cs = 0;
    for (int i=0; i<n; i++) {
        double th_i = 2*M_PI * cs / s;
        v[0][b[i]] = Vector(cos(th_i), sin(th_i), 0);
        cs += norm(mesh.vertices[b[i+1]]-mesh.vertices[b[i]]);
    }

    for (int n = 0; n < nb_iter; n++) {
        v[n+1].resize(n0);
        for (int i=0; i<n0; i++) {
            if (!on_boundary[i]) {
                int K = nbs[i].size();
                Vector pos;
                for (int j: nbs[i]) pos += v[n][j];
                v[n+1][i] = pos / K;
            } else {
                v[n+1][i] = v[n][i];
            }
        }
    }

    return v[nb_iter];
}

int main() {
    TriangleMesh mesh;
    // mesh.readOBJ("cat_model/cat.obj");
    mesh.readOBJ("cuuube.obj");

    auto tutte = Tutte(mesh, 10);

    for (Vector v: tutte) {
        printf("v %f %f %f\n", v[0], v[1], v[2]);
    }
    for (TriangleIndices tri: mesh.indices) {
        printf("f %d %d %d\n", tri.vtxi+1, tri.vtxj+1, tri.vtxk+1);
    }
}