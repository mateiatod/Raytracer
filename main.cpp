#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define M_P

#include <iostream>
#include <vector>
#include <string>
#include <list>
#include <stdio.h>
#include <algorithm>
#include <random>
#include <cmath>
#include <chrono>
#include <omp.h>


std::random_device rd;
thread_local std::mt19937 gen(rd());

class Vector {
public:
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        coords[0] = x;
        coords[1] = y;
        coords[2] = z;
    }

    Vector& operator+=(const Vector& b) {
    coords[0] += b[0];
    coords[1] += b[1];
    coords[2] += b[2];
    return * this;
    }

    const double& operator[](int i) const {
        return coords[i];
    }

    double& operator[](int i) {
        return coords[i];
    }

    double norm2() const {
        return coords[0] * coords[0] + coords[1] * coords[1] + coords[2] * coords[2];
    }

    double norm() const {
        return sqrt(norm2());
    }

    void normalize() {
        double n = norm();
        coords[0] /= n;
        coords[1] /= n;
        coords[2] /= n;
    }

    Vector operator-(const Vector &b) const {
        return Vector(coords[0] - b[0], coords[1] - b[1], coords[2] - b[2]);
    }

    Vector operator*(const double a) const {
        return Vector(coords[0] * a, coords[1] * a, coords[2] * a);
    }
    
    Vector operator*(const Vector &a) const{
        return Vector(coords[0] * a[0], coords[1] * a[1], coords[2] * a[2]);
    }

    Vector operator/(const double a) const {
        return Vector(coords[0] / a, coords[1] / a, coords[2] / a);
    }

    double dot(const Vector &b) const {
        return coords[0] * b[0] + coords[1] * b[1] + coords[2] * b[2];
    }

    Vector cross(const Vector &b) const {
        return Vector(coords[1] * b[2] - coords[2] * b[1], coords[2] * b[0] - coords[0] * b[2], coords[0] * b[1] - coords[1] * b[0]);
    }

private:
    double coords[3];
};

Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

class Ray {
public:
    Ray(const Vector &origin = Vector(), const Vector &direction = Vector()) {
    this->origin = origin;
    this->direction = direction;
    }

    Vector origin, direction;
};

struct Intersection {
    Vector P, N, albedo;
    double distance;
    bool mirror;

    Intersection(Vector P = Vector(0, 0, 0), Vector N = Vector(0, 0, 0), 
                 double distance = std::numeric_limits<double>::max(), Vector albedo = Vector(0, 0, 0)) {
        this->P = P;
        this->N = N;
        this->distance = distance;
        this->albedo = albedo;
        this->mirror = false; // default value
    }
};

class Object{
    public:
        Object(Vector albedo, bool mirror){
        this->albedo = albedo;
        this->mirror = mirror;
        }

        Vector albedo;
        virtual bool intersect(const Ray &r, Intersection &curr_inter) const = 0;
        bool mirror;
};

class Sphere: public Object {
public:

    Sphere(const Vector &C, const double &R, const Vector &albedo, bool mirror = false) : C(C), R(R), ::Object(albedo, mirror) {};

    bool intersect(const Ray &r, Intersection &curr_inter) const {
        Vector u = r.direction;
        Vector O = r.origin;
        double b = u.dot(O - C);
        double c = (O - C).norm2() - R * R;
        double delta = b * b - c;
        if (delta >= 0) {
            double t1 = (u.dot(C - O) - sqrt(delta));
            double t2 = (u.dot(C - O) + sqrt(delta));
            if (t2 < 0) {
                return false;
            }
            if (t1 > 0) {
                curr_inter.distance = t1;
            }
            else{
                curr_inter.distance = t2;
            }
            curr_inter.P = O + u * curr_inter.distance;
            curr_inter.N = curr_inter.P - C;
            curr_inter.N.normalize();
            curr_inter.albedo = albedo;
            curr_inter.mirror = mirror;
            return true;
        }
        return false;
    }

    Vector C;
    double R;
};

struct TriangleIndices {
public:
	TriangleIndices(int vtx0 = -1, int vtx1 = -1, int vtx2 = -1, int normal0 = -1, int normal1 = -1, 
            int normal2 = -1, int uv0 = -1, int uv1 = -1, int uv2 = -1, int group = -1){
        vtxindices[0] = vtx0;
        vtxindices[1] = vtx1;
        vtxindices[2] = vtx2;

        normalindices[0] = normal0;
        normalindices[1] = normal1;
        normalindices[2] = normal2;

        uvindices[0] = uv0;
        uvindices[1] = uv1;
        uvindices[2] = uv2;

        group = group;
    }
    int vtxindices[3];     // Refers to 3 indices in the vertices array of the Mesh class
    int normalindices[3];  // Refers to 3 indices in the normal array of the Mesh class
    int uvindices[3];

	int group;       // face group
};

class BoundingBox {
    public:

        BoundingBox() {
            Bmin = Vector(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
            Bmax = Vector(std::numeric_limits<double>::min(), std::numeric_limits<double>::min(), std::numeric_limits<double>::min());
        }

        bool intersect(const Ray& r) const {
            Vector O = r.origin;
            Vector u = r.direction;

            double xmin = (Bmin[0] - O[0]) / u[0];
            double xmax = (Bmax[0] - O[0]) / u[0];

            double ymin = (Bmin[1] - O[1]) / u[1];
            double ymax = (Bmax[1] - O[1]) / u[1];

            double zmin = (Bmin[2] - O[2]) / u[2];
            double zmax = (Bmax[2] - O[2]) / u[2];

            double tmin = std::max({std::min(xmin, xmax), std::min(ymin, ymax), std::min(zmin, zmax)});
            double tmax = std::min({std::max(xmin, xmax), std::max(ymin, ymax), std::max(zmin, zmax)});

            if (tmax > 0 && tmax > tmin) 
                return true;
            return false;
        }

        Vector compute_diag() const {
            return Bmax - Bmin;
        }

        void include_triangle(TriangleIndices& triangle, std::vector<Vector> vertices) {
            for (int i = 0; i < 3; ++i) {
                const Vector& point = vertices[triangle.vtxindices[i]];
                Bmin = Vector(std::min(Bmin[0], point[0]), std::min(Bmin[1], point[1]), std::min(Bmin[2], point[2]));
                Bmax = Vector(std::max(Bmax[0], point[0]), std::max(Bmax[1], point[1]), std::max(Bmax[2], point[2]));
    }
        }

        int get_longest(Vector diag) const {
        if (diag[0] >= diag[1] && diag[0] >= diag[2]) return 0;
        else if (diag[1] >= diag[0] && diag[1] >= diag[2]) return 1;
        else return 2;
    }

        Vector Bmin, Bmax;
};

class BVHNode {
public:
    BoundingBox bbox;
    BVHNode* right, *left;
    size_t starting_triangle, ending_triangle;

    BVHNode(){
    right = nullptr;
    left = nullptr;
    starting_triangle = -1;
    ending_triangle = -1;
    }


    ~BVHNode() {
        delete right;
        delete left;
    }
};

class TriangleMesh: public Object {
public:

    TriangleMesh(Vector albedo, bool mirror): Object(albedo, mirror){
    this->albedo = albedo;
    this->mirror = mirror;
    }


    ~TriangleMesh() {
        delete root_node;
    }

    bool intersect_triangle(const Ray& ray, size_t i, Intersection &intersection) const {
        const TriangleIndices &curr_triangle = mesh_list[i];
        const Vector &A = vertices[curr_triangle.vtxindices[0]];
        const Vector &B = vertices[curr_triangle.vtxindices[1]];
        const Vector &C = vertices[curr_triangle.vtxindices[2]];

        Vector AO = A - ray.origin;
        Vector e1 = B - A;
        Vector e2 = C - A;
        Vector u = ray.direction;
        Vector vec = AO.cross(u);
        Vector N = e1.cross(e2);
        double triangle_normal = u.dot(N);
        if (std::abs(triangle_normal) < 0.000001){
            return false;
        }

        double beta = e2.dot(vec)/triangle_normal;
        if (beta < 0 || beta > 1){
            return false;
        }
        double gamma = -e1.dot(vec)/triangle_normal; 
        if (gamma < 0 || gamma + beta> 1){
            return false;
        }

        double distance = AO.dot(N)/triangle_normal;

        if (distance > 0.000001 && distance < intersection.distance) {
            intersection.N = e1.cross(e2);
            intersection.P = ray.origin + ray.direction * distance;
            intersection.distance = distance;
            intersection.N.normalize();
            return true;
        }
        return false;
    }

    Vector compute_barycenter(TriangleIndices &triangle){
        Vector barycenter(0, 0, 0);
        for (int j = 0; j < 3; ++j) {
            barycenter += vertices[triangle.vtxindices[j]];
        }
        barycenter = barycenter / 3.0;
        return barycenter;
    }

    BVHNode* build_tree(size_t starting_triangle, size_t ending_triangle, int min_triangles) {
        BVHNode* curr_node = new BVHNode();
        size_t num_triangles = ending_triangle - starting_triangle;
        fit_bbox(curr_node->bbox, starting_triangle, ending_triangle);

        curr_node->starting_triangle = starting_triangle;
        curr_node->ending_triangle = ending_triangle;
        fit_bbox(curr_node->bbox,starting_triangle, ending_triangle); // BBox from starting_triangle included to ending_triangle excluded
        

        if (num_triangles <= min_triangles) {
            curr_node->starting_triangle = starting_triangle;
            curr_node->ending_triangle = ending_triangle;
            return curr_node;
        }
        Vector diag = curr_node->bbox.compute_diag();
        Vector middle_diag = curr_node->bbox.Bmin + diag * 0.5;
        int longest_axis = curr_node->bbox.get_longest(diag);
        int pivot_index = starting_triangle;
        for (int i = starting_triangle; i < ending_triangle; ++i) {
            Vector barycenter = compute_barycenter(mesh_list[i]);
            // the swap below guarantees triangles whose barycenter are smaller than middle_diag are before "pivot_index"
            if (barycenter[longest_axis] < middle_diag[longest_axis]) {
                std::swap(mesh_list[i], mesh_list[pivot_index]);
                pivot_index++;
            }
        }
        // Stopping criterion
        if (pivot_index <= starting_triangle || pivot_index >= ending_triangle || ending_triangle - starting_triangle < 5) {
            return curr_node;
        }
        curr_node->left = build_tree(starting_triangle, pivot_index, min_triangles);
        curr_node->right = build_tree(pivot_index, ending_triangle, min_triangles);
        return curr_node;
    }

    void intersect_BVH(const Ray &ray, BVHNode *curr_node, Intersection &intersection) const {
        if (!curr_node) {
            return;
        }
        if (!curr_node->bbox.intersect(ray)){
            return;
        }
        if (curr_node->left == nullptr && curr_node->right == nullptr) {
            for (size_t i = curr_node->starting_triangle; i < curr_node->ending_triangle; ++i) {
                Intersection curr_inter = intersection;
                if (intersect_triangle(ray, i, curr_inter) && curr_inter.distance < intersection.distance) {
                    intersection = curr_inter;
                }
            }
        } else {
            intersect_BVH(ray, curr_node->left, intersection);
            intersect_BVH(ray, curr_node->right, intersection);
        }
    }

    bool intersect(const Ray& ray, Intersection &intersection) const {
        intersect_BVH(ray, root_node, intersection);

        if (intersection.distance < std::numeric_limits<double>::max()) {
            intersection.albedo = albedo;
            intersection.mirror = mirror;
            return true;
        } else {
            return false;
        }
    }

    void fit_bbox(BoundingBox &bbox, size_t starting_triangle, size_t ending_triangle) {
        for (size_t i = starting_triangle; i < ending_triangle; ++i) {
            TriangleIndices &curr_triangle = mesh_list[i];
            bbox.include_triangle(curr_triangle, vertices);
        }
    }

    void place(const double scalar = 1, const Vector offset = Vector(0, 0, 0)){
        for (int i = 0; i< vertices.size(); ++i){
            vertices[i] = vertices[i]*scalar + offset;
        }
    }

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
                    if (i0 < 0) t.vtxindices[0] = vertices.size() + i0; else	t.vtxindices[0] = i0 - 1;
                    if (i1 < 0) t.vtxindices[1] = vertices.size() + i1; else	t.vtxindices[1] = i1 - 1;
                    if (i2 < 0) t.vtxindices[2] = vertices.size() + i2; else	t.vtxindices[2] = i2 - 1;
                    if (j0 < 0) t.uvindices[0] = uvs.size() + j0; else	t.uvindices[0] = j0 - 1;
                    if (j1 < 0) t.uvindices[1] = uvs.size() + j1; else	t.uvindices[1] = j1 - 1;
                    if (j2 < 0) t.uvindices[2] = uvs.size() + j2; else	t.uvindices[2] = j2 - 1;
                    if (k0 < 0) t.normalindices[0] = normals.size() + k0; else	t.normalindices[0] = k0 - 1;
                    if (k1 < 0) t.normalindices[1] = normals.size() + k1; else	t.normalindices[1] = k1 - 1;
                    if (k2 < 0) t.normalindices[2] = normals.size() + k2; else	t.normalindices[2] = k2 - 1;
                    mesh_list.push_back(t);
                } else {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6) {
                        if (i0 < 0) t.vtxindices[0] = vertices.size() + i0; else	t.vtxindices[0] = i0 - 1;
                        if (i1 < 0) t.vtxindices[1] = vertices.size() + i1; else	t.vtxindices[1] = i1 - 1;
                        if (i2 < 0) t.vtxindices[2] = vertices.size() + i2; else	t.vtxindices[2] = i2 - 1;
                        if (j0 < 0) t.uvindices[0] = uvs.size() + j0; else	t.uvindices[0] = j0 - 1;
                        if (j1 < 0) t.uvindices[1] = uvs.size() + j1; else	t.uvindices[1] = j1 - 1;
                        if (j2 < 0) t.uvindices[2] = uvs.size() + j2; else	t.uvindices[2] = j2 - 1;
                        mesh_list.push_back(t);
                    } else {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3) {
                            if (i0 < 0) t.vtxindices[0] = vertices.size() + i0; else	t.vtxindices[0] = i0 - 1;
                            if (i1 < 0) t.vtxindices[1] = vertices.size() + i1; else	t.vtxindices[1] = i1 - 1;
                            if (i2 < 0) t.vtxindices[2] = vertices.size() + i2; else	t.vtxindices[2] = i2 - 1;
                            mesh_list.push_back(t);
                        } else {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0) t.vtxindices[0] = vertices.size() + i0; else	t.vtxindices[0] = i0 - 1;
                            if (i1 < 0) t.vtxindices[1] = vertices.size() + i1; else	t.vtxindices[1] = i1 - 1;
                            if (i2 < 0) t.vtxindices[2] = vertices.size() + i2; else	t.vtxindices[2] = i2 - 1;
                            if (k0 < 0) t.normalindices[0] = normals.size() + k0; else	t.normalindices[0] = k0 - 1;
                            if (k1 < 0) t.normalindices[1] = normals.size() + k1; else	t.normalindices[1] = k1 - 1;
                            if (k2 < 0) t.normalindices[2] = normals.size() + k2; else	t.normalindices[2] = k2 - 1;
                            mesh_list.push_back(t);
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
                        if (i0 < 0) t2.vtxindices[0] = vertices.size() + i0; else	t2.vtxindices[0] = i0 - 1;
                        if (i2 < 0) t2.vtxindices[1] = vertices.size() + i2; else	t2.vtxindices[1] = i2 - 1;
                        if (i3 < 0) t2.vtxindices[2] = vertices.size() + i3; else	t2.vtxindices[2] = i3 - 1;
                        if (j0 < 0) t2.uvindices[0] = uvs.size() + j0; else	t2.uvindices[0] = j0 - 1;
                        if (j2 < 0) t2.uvindices[1] = uvs.size() + j2; else	t2.uvindices[1] = j2 - 1;
                        if (j3 < 0) t2.uvindices[2] = uvs.size() + j3; else	t2.uvindices[2] = j3 - 1;
                        if (k0 < 0) t2.normalindices[0] = normals.size() + k0; else	t2.normalindices[0] = k0 - 1;
                        if (k2 < 0) t2.normalindices[1] = normals.size() + k2; else	t2.normalindices[1] = k2 - 1;
                        if (k3 < 0) t2.normalindices[2] = normals.size() + k3; else	t2.normalindices[2] = k3 - 1;
                        mesh_list.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    } else {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2) {
                            if (i0 < 0) t2.vtxindices[0] = vertices.size() + i0; else	t2.vtxindices[0] = i0 - 1;
                            if (i2 < 0) t2.vtxindices[1] = vertices.size() + i2; else	t2.vtxindices[1] = i2 - 1;
                            if (i3 < 0) t2.vtxindices[2] = vertices.size() + i3; else	t2.vtxindices[2] = i3 - 1;
                            if (j0 < 0) t2.uvindices[0] = uvs.size() + j0; else	t2.uvindices[0] = j0 - 1;
                            if (j2 < 0) t2.uvindices[1] = uvs.size() + j2; else	t2.uvindices[1] = j2 - 1;
                            if (j3 < 0) t2.uvindices[2] = uvs.size() + j3; else	t2.uvindices[2] = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            mesh_list.push_back(t2);
                        } else {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                            if (nn == 2) {
                                if (i0 < 0) t2.vtxindices[0] = vertices.size() + i0; else	t2.vtxindices[0] = i0 - 1;
                                if (i2 < 0) t2.vtxindices[1] = vertices.size() + i2; else	t2.vtxindices[1] = i2 - 1;
                                if (i3 < 0) t2.vtxindices[2] = vertices.size() + i3; else	t2.vtxindices[2] = i3 - 1;
                                if (k0 < 0) t2.normalindices[0] = normals.size() + k0; else	t2.normalindices[0] = k0 - 1;
                                if (k2 < 0) t2.normalindices[1] = normals.size() + k2; else	t2.normalindices[1] = k2 - 1;
                                if (k3 < 0) t2.normalindices[2] = normals.size() + k3; else	t2.normalindices[2] = k3 - 1;								
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                mesh_list.push_back(t2);
                            } else {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1) {
                                    if (i0 < 0) t2.vtxindices[0] = vertices.size() + i0; else	t2.vtxindices[0] = i0 - 1;
                                    if (i2 < 0) t2.vtxindices[1] = vertices.size() + i2; else	t2.vtxindices[1] = i2 - 1;
                                    if (i3 < 0) t2.vtxindices[2] = vertices.size() + i3; else	t2.vtxindices[2] = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    mesh_list.push_back(t2);
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
    BVHNode* root_node;
	std::vector<TriangleIndices> mesh_list;
	std::vector<Vector> vertices, normals, uvs, vertexcolors;
    Vector albedo;
    bool mirror;
    
};

const int W = 512;
const int H = 512;
const double alpha = 60. * M_PI / 180.;
const Vector camera_center(0, 0, 55);
const Vector source(-10, 20, 40);
const long long light_intensity = 3 * 1e8;
std::default_random_engine engine(10);
std::uniform_real_distribution<double> uniform(0.0, 1.0);

class Scene {
public:
    std::vector<const Object*> objects;
    void addGeometry(const Object *obj) {
        objects.push_back(obj);
    }

    bool intersect(const Ray &r, Intersection &best_inter) {
        bool intersection_found = false;
        double min_dist = std::numeric_limits<double>::max();
        for (auto &object : objects) {
            Intersection curr_inter;
            if (object->intersect(r, curr_inter) && curr_inter.distance < min_dist) {
                min_dist = curr_inter.distance;
                best_inter = curr_inter;
                intersection_found = true;
            }
        }
        return intersection_found;
    } //TODO

    int compute_visibility(const Vector &point, const Vector &source) {
        Vector u = source - point;
        double distance = u.norm();
        u.normalize();
        Intersection backwards_inter;
        Ray backwards_ray(point + u * 0.000001, u);
        bool intersected = intersect(backwards_ray, backwards_inter);

        if (!intersected) {
            return 0;
        } else if (backwards_inter.distance < distance) {
            return 0;
        }
        return 1; //TODO
    }

    Vector getColor(const Ray &ray, int n_bounces, int max_distance, const Vector source, const long long light_intensity) {
        if (n_bounces < 0)
            return Vector(0.0, 0.0, 0.0); // Terminates recursion at some point

        Intersection curr_inter;
        if (!intersect(ray, curr_inter)) {
            return Vector(0., 0., 0.);
        } 
        Vector N = curr_inter.N, P = curr_inter.P;
        if (curr_inter.mirror) {
            Vector u = ray.direction;
            Vector refl_direction = u - N * 2 * u.dot(N);
            Ray ray_refl(curr_inter.P + refl_direction * 0.0001, refl_direction);
            return getColor(ray_refl, n_bounces - 1, max_distance, source, light_intensity);
        } else {
            int visib = compute_visibility(curr_inter.P, source);
            Vector Lo = curr_inter.albedo * (light_intensity / (4.0 * M_PI * (source - curr_inter.P).norm2())) *
                        visib * std::max(0.0, curr_inter.N.dot(source - curr_inter.P)) / M_PI;

            double r1 = uniform(engine);
            double r2 = uniform(engine);
            double x = cos(2 * M_PI * r1) * sqrt(1 - r2);
            double y = sin(2 * M_PI * r1) * sqrt(1 - r2);
            double z = sqrt(r2);

            
            Vector T1;

            if (abs(N[0]) < abs(N[1]) && abs(N[0]) < abs(N[2])){
                T1 = Vector(0, -N[2], N[1]);
            }
            else{
                if(abs(N[1]) < abs(N[2])){
                    T1 = Vector(-N[2], 0, N[0]);
                }
                else{
                    T1 = Vector(-N[1], N[0], 0);
                }
            }
            T1.normalize();

            Vector random_vector =  T1 * x + N.cross(T1) * y + N * z;

            Lo += curr_inter.albedo * getColor(Ray(curr_inter.P + random_vector * 0.0001, random_vector), n_bounces - 1, max_distance, source, light_intensity);
            return Lo;
        }

    }

};


int main() {
    std::vector<unsigned char> image(W * H * 3, 0);

    Scene scene;
    scene.addGeometry(new Sphere(Vector(0, 16, 0), 5, Vector(0.1, 0.1, 0.1)));
    scene.addGeometry(new Sphere(Vector(7, 11, -4), 5, Vector(1., 1., 1.), true));
    scene.addGeometry(new Sphere(Vector(-7, 11, -4), 5, Vector(1., 1., 1.), true));
    scene.addGeometry(new Sphere(Vector(0, 0, -1000), 940, Vector(0, 1, 0)));
    scene.addGeometry(new Sphere(Vector(0, -1000, 0), 990, Vector(0, 0, 1)));
    scene.addGeometry(new Sphere(Vector(0, 1000, 0), 940, Vector(1, 0, 0)));
    scene.addGeometry(new Sphere(Vector(0, 0, 1000), 940, Vector(1, 0, 1)));
    scene.addGeometry(new Sphere(Vector(-1000, 0, 0), 940, Vector(1, 1, 0)));
    scene.addGeometry(new Sphere(Vector(1000, 0, 0), 940, Vector(0, 1, 1)));

    int max_distance = 5;
    int NB_PATHS = 64;

    TriangleMesh mesh(Vector(1., 1., 1.), false);
    mesh.readOBJ("cat.obj");
    mesh.place(0.6, Vector(0, -10, 0));
    mesh.root_node = mesh.build_tree(0, mesh.mesh_list.size(), 4);
    scene.addGeometry(&mesh);

    auto start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for schedule(dynamic, 1) 
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            Vector pixelColor(0, 0, 0);
            for (int k = 0; k < NB_PATHS; ++k) {
                double randomX = uniform(gen);
                double randomY = uniform(gen);
                // double randomX = 0 , randomY = 0;
                Vector dir = Vector(j - W / 2. + 0.5 + randomX * 0.5, -i + H / 2. + 0.5 + randomY * 0.5, -W/(2.*tan(alpha/2.)));
                dir.normalize();
                Ray r(camera_center, dir);
                pixelColor += scene.getColor(r, max_distance, max_distance, source, light_intensity);
                
            }

            Vector color = pixelColor / NB_PATHS;
            double gamma = 0;
            gamma = 2.2;
            image[(i * W + j) * 3 + 0] = std::min(255., std::pow(color[0], 1 / gamma));
            image[(i * W + j) * 3 + 1] = std::min(255., std::pow(color[1], 1 / gamma));
            image[(i * W + j) * 3 + 2] = std::min(255., std::pow(color[2], 1 / gamma));
        }
    }
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    double elapsedSeconds = elapsed.count();

    std::cout << "Elapsed time: " << elapsedSeconds << " seconds " << std::endl;

    stbi_write_png("image.png", W, H, 3, &image[0], 0);
    return 0;
}