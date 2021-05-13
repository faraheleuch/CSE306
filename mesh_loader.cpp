#include <string>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <list>

#include "vector.cpp"
#include "objects.cpp"

class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
	};
	int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
	int uvi, uvj, uvk;  // indices within the uv coordinates array
	int ni, nj, nk;  // indices within the normals array
	int group;       // face group
};

class TriangleMesh;

class BoundingBox {
	public:

	Vector Bmin, Bmax;

    BoundingBox() {}
	BoundingBox(int start_index, int end_index, TriangleMesh& T);
    bool intersect(const Ray& r, double& t) const;
};

class BVH {
	public:
	BoundingBox bbox;
	int starting_triangle;
	int ending_triangle;
	BVH *child_left = nullptr;
	BVH *child_right = nullptr;
	
	BVH() {}
	BVH(int starting_triangle, int ending_triangle, TriangleMesh& T);
	Intersection intersect(const Ray& r, const TriangleMesh& t);
};


class TriangleMesh : public Geometry {
public:
  ~TriangleMesh() {}
	TriangleMesh() {};
	
	void readOBJ(const char* obj) {

		char matfile[255];
		char grp[255];

		FILE* f;
		f = fopen(obj, "r");
		int curGroup = -1;
		while (!feof(f)) {
			char line[255];
			if (!fgets(line, 255, f)) break;

			std::string linetrim(line);
			linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
			strcpy(line, linetrim.c_str());

			if (line[0] == 'u' && line[1] == 's') {
				sscanf(line, "usemtl %[^\n]\n", grp);
				curGroup++;
			}

			if (line[0] == 'v' && line[1] == ' ') {
				Vector vec;

				Vector col;
				if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
					col[0] = std::min(1., std::max(0., col[0]));
					col[1] = std::min(1., std::max(0., col[1]));
					col[2] = std::min(1., std::max(0., col[2]));

					vertices.push_back(vec);
					vertexcolors.push_back(col);

				} else {
					sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
					vertices.push_back(vec);
				}
			}
			if (line[0] == 'v' && line[1] == 'n') {
				Vector vec;
				sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
				normals.push_back(vec);
			}
			if (line[0] == 'v' && line[1] == 't') {
				Vector vec;
				sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
				uvs.push_back(vec);
			}
			if (line[0] == 'f') {
				TriangleIndices t;
				int i0, i1, i2, i3;
				int j0, j1, j2, j3;
				int k0, k1, k2, k3;
				int nn;
				t.group = curGroup;

				char* consumedline = line + 1;
				int offset;

				nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
				if (nn == 9) {
					if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
					if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
					if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
					if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
					if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
					if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
					if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
					if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
					if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
					indices.push_back(t);
				} else {
					nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
					if (nn == 6) {
						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
						if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
						if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
						if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
						indices.push_back(t);
					} else {
						nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
						if (nn == 3) {
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							indices.push_back(t);
						} else {
							nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
							if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
							if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
							indices.push_back(t);
						}
					}
				}

				consumedline = consumedline + offset;

				while (true) {
					if (consumedline[0] == '\n') break;
					if (consumedline[0] == '\0') break;
					nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
					TriangleIndices t2;
					t2.group = curGroup;
					if (nn == 3) {
						if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
						if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
						if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
						if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
						if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
						if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
						if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
						if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
						if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
						indices.push_back(t2);
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						k2 = k3;
					} else {
						nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
							if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
							if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
							if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
							if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
							if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
							consumedline = consumedline + offset;
							i2 = i3;
							j2 = j3;
							indices.push_back(t2);
						} else {
							nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
							if (nn == 2) {
								if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
								if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
								if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
								if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
								if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
								if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;								
								consumedline = consumedline + offset;
								i2 = i3;
								k2 = k3;
								indices.push_back(t2);
							} else {
								nn = sscanf(consumedline, "%u%n", &i3, &offset);
								if (nn == 1) {
									if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
									if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
									if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
									consumedline = consumedline + offset;
									i2 = i3;
									indices.push_back(t2);
								} else {
									consumedline = consumedline + 1;
								}
							}
						}
					}
				}

			}

		}
		fclose(f);

	}

	Intersection intersect(const Ray& r) const {
		return root->intersect(r, *this);
	}

	void transform(Vector offset, double scale) {
		for (auto& v: vertices) {
			v = v * scale + offset;
		}
	}

	void init() {
		root = new BVH(0, indices.size(), *this);
	}

	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;

	BVH* root;
};

/*
Intersection TriangleMesh::intersect(const Ray& r) const {

	double closest = inf;
	Intersection res = Intersection();

	for(int i=0; i < indices.size(); i++) {
		Vector A = vertices[indices[i].vtxi];
		Vector B = vertices[indices[i].vtxj];
		Vector C = vertices[indices[i].vtxk];
		Vector O = r.origin;
		Vector u = r.direction;
		auto e1 = B - A;
		auto e2 = C - A;
		auto N = cross(e1,e2);

		double beta = dot(e2, cross((A-O),u)) / dot(u,N);
		double gamma = - dot(e1, cross((A-O),u)) / dot(u,N);
		double alpha = 1 - beta - gamma;
		double t = dot(A-O,N) / dot(u,N);

		if((0 <= alpha && alpha <= 1) && (0 <= beta && beta <= 1) && (0 <= gamma && gamma <= 1)) {
			Vector P = alpha*A + beta*B + gamma*C;
			Intersection inter(1, P, normalization(N), t, index);

			if (inter.distance < closest && inter.distance > 1e-5){
				closest = inter.distance;
				res = inter;
			}
		}
	}
	return res;
}
*/

BoundingBox::BoundingBox(int start_index, int end_index, TriangleMesh& T) {
	double x_min, y_min, z_min;
	double x_max, y_max, z_max;
	
	x_min = y_min = z_min =  1e69;
	x_max = y_max = z_max = -1e69;
	for(int i = start_index; i < end_index; i++) {
		auto tri = T.indices[i];
		for (int j: {tri.vtxi, tri.vtxj, tri.vtxk}) {
			Vector v = T.vertices[j];
			// add point to BoundingBox
			x_max = std::max(x_max, v[0]);
			y_max = std::max(y_max, v[1]);
			z_max = std::max(z_max, v[2]);
			x_min = std::min(x_min, v[0]);
			y_min = std::min(y_min, v[1]);
			z_min = std::min(z_min, v[2]);
		}
	}
	Bmin = Vector(x_min, y_min, z_min);
	Bmax = Vector(x_max, y_max, z_max);
}

bool BoundingBox::intersect(const Ray& r, double& t) const {
	double  x_min = Bmin[0],
			y_min = Bmin[1],
			z_min = Bmin[2];
	double  x_max = Bmax[0],
			y_max = Bmax[1],
			z_max = Bmax[2];

	double t_0x = (x_min - r.origin[0])/r.direction[0];
	double t_0y = (y_min - r.origin[1])/r.direction[1];
	double t_0z = (z_min - r.origin[2])/r.direction[2];

	double t_1x = (x_max - r.origin[0])/r.direction[0];
	double t_1y = (y_max - r.origin[1])/r.direction[1];
	double t_1z = (z_max - r.origin[2])/r.direction[2];

	if (t_0x > t_1x) std::swap(t_0x,t_1x);
	if (t_0y > t_1y) std::swap(t_0y,t_1y);
	if (t_0z > t_1z) std::swap(t_0z,t_1z); 

	if (std::min({t_1x,t_1y,t_1z}) > std::max({t_0x,t_0y,t_0z})) {
		t = std::max({t_0x,t_0y,t_0z});
		return true;
	}
	
	return false;	
}

BVH::BVH(int starting_triangle, int ending_triangle, TriangleMesh& T) {
	this->bbox = BoundingBox(starting_triangle , ending_triangle, T) ;
	this->starting_triangle = starting_triangle ;
	this->ending_triangle = ending_triangle ;
	Vector diag = this->bbox.Bmax - this->bbox.Bmin;
	Vector middle_diag = this->bbox.Bmin + diag * .5;
	int longuest_axis = (diag[1] > diag[0]);
	if (diag[2] > diag[longuest_axis]) longuest_axis = 2;

	int pivot_index = starting_triangle;
	for (int i=starting_triangle; i<ending_triangle; i++) {
		Vector A = T.vertices[T.indices[i].vtxi];
		Vector B = T.vertices[T.indices[i].vtxj];
		Vector C = T.vertices[T.indices[i].vtxk];

		Vector barycenter = (A+B+C)/3;

		if (barycenter[longuest_axis] < middle_diag[longuest_axis]) {
			std::swap(T.indices[i], T.indices[pivot_index]);
			pivot_index++;
		}
	}
	if( pivot_index<=starting_triangle || pivot_index>=ending_triangle-1 || ending_triangle-starting_triangle<5 ) return;

		this->child_left = new BVH(starting_triangle , pivot_index, T) ;
		this->child_right = new BVH(pivot_index , ending_triangle, T) ;
}

Intersection BVH::intersect(const Ray& r, const TriangleMesh& T) {
	double t;
	Intersection res;
	double best_inter_distance = std::numeric_limits<double>::max();

	if (! this->bbox.intersect(r,t)) return res ;
	
	std::list<BVH*> nodes_to_visit;
	nodes_to_visit.push_front(this);
	
	while (!nodes_to_visit.empty()) {
		BVH* curNode = nodes_to_visit.back();
		nodes_to_visit.pop_back();
		if(curNode -> child_left) {
			if(curNode->child_left->bbox.intersect(r,t)) {
				if(t < best_inter_distance) {
					nodes_to_visit.push_back(curNode->child_left);
				}
			}
			if(curNode->child_right->bbox.intersect(r,t)) {
				if(t < best_inter_distance) {
					nodes_to_visit.push_back(curNode->child_right);
				}
			}
		} else {
			for (int i=this->starting_triangle; i < this->ending_triangle; i++) {
				

				Vector A = T.vertices[T.indices[i].vtxi];
				Vector B = T.vertices[T.indices[i].vtxj];
				Vector C = T.vertices[T.indices[i].vtxk];
				Vector O = r.origin;
				Vector u = r.direction;
				auto e1 = B - A;
				auto e2 = C - A;
				auto N = cross(e1,e2);

				double beta = dot(e2, cross((A-O),u)) / dot(u,N);
				double gamma = - dot(e1, cross((A-O),u)) / dot(u,N);
				double alpha = 1 - beta - gamma;
				double t = dot(A-O,N) / dot(u,N);

				if((0 <= alpha && alpha <= 1) && (0 <= beta && beta <= 1) && (0 <= gamma && gamma <= 1)) {
					Vector P = alpha*A + beta*B + gamma*C;
					Intersection inter(1, P, normalization(N), t, T.index);

					if (inter.distance < best_inter_distance && inter.distance > 1e-5){
						best_inter_distance = inter.distance;
						res = inter;
					}
				}
			}
		}
	}
	
	return res;
}