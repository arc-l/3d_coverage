#ifndef SOLUTION_H
#define SOLUTION_H

#include "gurobi_c++.h"
#include "Vec.h"

#include<pcl/io/pcd_io.h>
#include<pcl/point_types.h>
#include<pcl/io/ply_io.h>
#include<pcl/io/vtk_lib_io.h>
#include<pcl/Vertices.h>

#include <list>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/AABB_triangle_primitive.h>

#include <random>
#include <iterator>
#include<vector>
#include<cmath>
#include<cstdlib>
#include<cstdio>
#include<set>
#include<ctime>
#include<cassert>
#include<algorithm>
#include <chrono> 
#include <iostream>   
#include <fstream>   
#include<functional>
#include<cstring>
#include<omp.h>

using namespace std;
using namespace std::chrono;

typedef CGAL::Simple_cartesian<double> K;
typedef K::Segment_3 Segment;
typedef K::FT FT;
typedef K::Ray_3 Ray;
typedef K::Line_3 Line;
typedef K::Point_3 Point;
typedef K::Triangle_3 Triangle;
typedef std::list<Triangle>::iterator Iterator;
typedef CGAL::AABB_triangle_primitive<K, Iterator> Primitive;
typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;

// const int R[5] = {205, 89, 53, 58, 116};
// const int G[5] = {109, 72, 107, 95, 61};
// const int B[5] = {1, 133, 59, 106, 100};

// orange(255, 140, 0)
// pink(255,105,180)
// green(0,128,0)
// purple(128,0,128)
// cyan(0,255,255)

// Orange, Green, Pink Purple, Cyan
const int R[8] = {205, 68,  238,  120, 122};
const int G[8] = {109, 164, 77,   80,  215};
const int B[8] = {1,   22,  192,  219, 236};

typedef pcl::PointXYZRGB PointXYZRGB;
// Default pcl version on Ubuntu is faulty now
// //****************************************
// // PCL v0.10.0 (default for ubuntu focal LTS) does not have proper constructor for pointxyzrgb
// //****************************************
// inline pcl::PointXYZRGB PointXYZRGB(float x, float y, float z, uint8_t r, uint8_t g, uint8_t b){
// 	pcl::PointXYZRGB ret;
// 	ret.x = x; ret.y = y; ret.z = z; ret.r = r; ret.g = g; ret.b = b;
// 	return ret;
// }

class Solution{
public:
	GRBEnv env = GRBEnv(true);
	pcl::PointCloud<pcl::PointNormal>::Ptr cloud;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr loc_cloud;
	pcl::PolygonMesh mesh;
	Tree tree;
	double height = .0;
	// vector<double> xval, yval; // used for default sensor locations to select from
	vector<vec3> possible_locs;
	int num_of_pts, num_of_locations;
    double maxx, minx, maxy, miny, maxz, minz;
	Solution(){
		// env = GRBEnv(true);
		env.set("LogFile", "mip1.log");
		env.start();
		env.set(GRB_IntParam_OutputFlag, false);
		cloud = pcl::PointCloud<pcl::PointNormal>::Ptr(new pcl::PointCloud<pcl::PointNormal>);
		loc_cloud = pcl::PointCloud<pcl::PointXYZRGB>::Ptr(new pcl::PointCloud<pcl::PointXYZRGB>);
		// load pcd file from 
		if(pcl::io::loadPCDFile<pcl::PointNormal> ("uniform_out.pcd", *cloud) == -1){
			cerr<<"error loading pcd files"<<endl; exit(1);
		}
		num_of_pts = cloud->size();
		// load ply file
		if(pcl::io::loadPolygonFilePLY("model.ply", mesh) == -1){
			cerr<<"error loading ply files"<<endl; exit(1);
		}
        maxx = std::numeric_limits<double>::min();
		minx = std::numeric_limits<double>::max();
		miny = minx, maxy = maxx, minz= minx, maxz = maxx;
		for(const auto& pt: *cloud){
			double x1 = pt.x, y1 = pt.y, z1 = pt.z;
			minz = min(minz, z1); maxz = max(z1, maxz);
			maxx = max(maxx, x1); minx = min(x1, minx);
			miny = min(miny, y1); maxy = max(y1, maxy);
		}
	}
	void sub_sample_coarse(int new_num_of_locations, int new_num_of_pts, vector<vector<float>> & inter);
	void comp_possible_sensor_locations(int num);
    void build_AABB_Tree();
	void sub_sample(int sam_num = 20000);
	void precompute_2_1();
	void compute_2(int num_of_sensors);
	void precompute_3_1();
	void compute_3(int num_of_sensors, double thresh, bool Coarse = false, int new_num_of_locations = 100, int new_num_of_pts = 1000);
	void precompute_4_1();
	void compute_4(int num_of_sensors, double ratio, bool Coarse = false, int new_num_of_locations = 100, int new_num_of_pts = 1000);
	void visualize(int problem_t = 2, double thresh = 0.0);
	void local_imprv_3();
	void local_imprv_4();
};

#endif