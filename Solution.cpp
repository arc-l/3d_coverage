#include "Solution.h"

void Solution::sub_sample_coarse(int new_num_of_locations, int new_num_of_pts, vector<vector<float>> & inter){
    // Using less possible locations to save time for ilp computation
    // though it is not very necessary
    // resample to new_num_of_points
    // and uses new_num_of_location locations
    if(new_num_of_locations > num_of_locations || new_num_of_pts > new_num_of_pts){
        cout << "Sub Sample failed" << endl;
        exit(1);
    }
    vector<double> v(num_of_pts), v2(num_of_locations);
    iota(begin(v), end(v), 0);
    iota(begin(v2), end(v2), 0);
    std::random_device rd;
    std::mt19937 g(rd());
    shuffle(begin(v), end(v), g);
    shuffle(begin(v2), end(v2), g);
    v.resize(new_num_of_pts);
    v2.resize(new_num_of_locations);
    sort(begin(v), end(v));
    sort(begin(v2), end(v2));
    num_of_pts = new_num_of_pts;
    num_of_locations = new_num_of_locations;
    for(int i=0;i<num_of_locations;i++){
        for(int j=0;j<num_of_pts;j++){
            inter[i][j] = inter[v2[i]][v[j]];
        }
        inter[i].resize(num_of_pts);
    }
    inter.resize(num_of_locations);
    possible_locs.resize(num_of_locations);
}

void Solution::comp_possible_sensor_locations(int num){
    //***********************************************************************
    // default version uses a grid, 
    // and assumes the ceiling is axis-aligned, 
    // and z = zmax - (zmax - zmin) / 10
    //***********************************************************************
    double div = 0.01;
    int x_seg, y_seg;
    double madiv = maxx - minx + maxy - miny, midiv = 0;
    while(madiv - midiv > 1e-6 && (madiv - midiv) / madiv > 1e-6){
        double div = (madiv + midiv) / 2;
        if(ceil((maxx - minx) / div)  * ceil((maxy - miny) / div) <= num) {
            midiv = div;
        }else madiv = div;
    }
    div = madiv;
    x_seg = (maxx - minx) / div + 1;
    y_seg = (maxy - miny) / div + 1;
    double z0 = maxz - (maxz - minz) / 10;
    vector<double> xval, yval;
    xval.resize(x_seg), yval.resize(y_seg);
    for(int i=0;i<x_seg;i++){
        xval[i] = minx + (maxx-minx) / (x_seg+2) * (i+1.5);
    }
    for(int i=0;i<y_seg;i++){
        yval[i] = miny + (maxy-miny) / (y_seg+2) * (i+1.5);
    }
    for(int i=0;i<x_seg;i++){
        for(int j=0;j<y_seg;j++){
            possible_locs.emplace_back(xval[i], yval[j], z0);
        }
    }
    num_of_locations = possible_locs.size();
}

void Solution::build_AABB_Tree(){
    pcl::PointCloud<pcl::PointXYZ> pc2;
    pcl::fromPCLPointCloud2(mesh.cloud, pc2);
    std::list<Triangle> triangles;
    // triangles.clear();
    for(auto &p : mesh.polygons){
        pcl::PointXYZ a = pc2.points[p.vertices[0]];
        pcl::PointXYZ b = pc2.points[p.vertices[1]];
        pcl::PointXYZ c = pc2.points[p.vertices[2]]; 
        triangles.emplace_back(Point(a.x, a.y, a.z), Point(b.x, b.y, b.z), Point(c.x, c.y, c.z));
    }
    tree.clear();
    tree.insert(begin(triangles), end(triangles));
}

void Solution::sub_sample(int sam_num){
    
    //***********************************************************************
    // take only those points visibility to at least one of the possible sensor locations
    // and sub_sample to sam_num of points
    //***********************************************************************

    if(cloud->size() < sam_num){
        cerr << "Sample point size larger than the original pointcloud" << endl;
        exit(1);
    }
    pcl::PointCloud<pcl::PointNormal>::Ptr newcloud_for_save(new pcl::PointCloud<pcl::PointNormal>());
    pcl::PointCloud<pcl::PointNormal>::Ptr newcloud_for_save_subsam(new pcl::PointCloud<pcl::PointNormal>());
    for(const auto &pt: *cloud){
        double x1 = pt.x, y1 = pt.y, z1 = pt.z;
        for(const vec3 &v: possible_locs){
            double x0 = v.x, y0 = v.y, z0 = v.z;
            Segment s(Point(x1 * .99999 + x0 * .00001, y1 * .99999 + y0 * .00001, z1 * .99999 + z0 * .00001), 
                Point(x0 * .99999 + x1 * .00001, y0 * .99999 + y1 * .00001, z0 * .99999 + z1 * .00001));
            if(s.is_degenerate()) continue;
            if(!tree.do_intersect(s)){
                auto v1 = s.to_vector();
                if(v1[0] * pt.normal_x + v1[1] * pt.normal_y + v1[2] * pt.normal_z > 0){
                    newcloud_for_save->push_back(pt);
                    break;
                }
            }
        }
    }
    vector<int> v(newcloud_for_save->size());
    iota(begin(v), end(v), 0);
    std::random_device rd;
    std::mt19937 g(rd());
    shuffle(begin(v), end(v), g);
    for(int i=0;i<sam_num;i++){
        newcloud_for_save_subsam->push_back((*newcloud_for_save)[v[i]]);
    }
    if(pcl::io::savePCDFile<pcl::PointNormal> ("uniform_out.pcd", *newcloud_for_save_subsam) == -1){
        cout<<"Error saving files"<<endl;
    }
    cloud = newcloud_for_save_subsam;
    num_of_pts = sam_num;
}

void Solution::visualize(int problem_t, double thresh){
    //***********************************
    // visualize coverd point clouds by saving it into a colored pointcloud
    // problem t: model number applied
    //***********************************
    if(pcl::io::loadPCDFile<pcl::PointXYZRGB>("locations_t" + to_string(problem_t) + ".pcd", *loc_cloud) == -1){
        cerr<<"Error loading location cloud"<<endl; exit(1);
    }
    sort(begin(*loc_cloud), end(*loc_cloud), [](auto i, auto j){return i.x + i.y < j.x + j.y;});
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr covered_col_cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
    for(const auto &pt: *cloud){
        double x1 = pt.x, y1 = pt.y, z1 = pt.z;
        double dis = std::numeric_limits<double>::max(), ma_quality = 0;
        int cnt = 0, id = -1;
        double quality = 0;
        for(const auto &pt2: *loc_cloud){
            double x0 = pt2.x, y0 = pt2.y, z0 = pt2.z;
            Segment s(Point(x1 * .99999 + x0 * .00001, y1 * .99999 + y0 * .00001, z1 * .99999 + z0 * .00001), 
                Point(x0 * .99999 + x1 * .00001, y0 * .99999 + y1 * .00001, z0 * .99999 + z1 * .00001));
            if(s.is_degenerate()) continue;
            if(!tree.do_intersect(s)){
                if(problem_t == 3){
                    auto v1 = s.to_vector();
                    double tmp = max(0.0, (v1[0] * pt.normal_x + v1[1] * pt.normal_y + v1[2] * pt.normal_z) / 
                        sqrt( pt.normal_x * pt.normal_x  + pt.normal_y * pt.normal_y + pt.normal_z * pt.normal_z) 
                        /  pow(v1.squared_length(), 1.5));
                    quality += tmp;
                    if(tmp > ma_quality) {
                        ma_quality = tmp; id = cnt;
                    }
                }else if (problem_t = 2){
                    if(s.squared_length() < dis){
                        dis = s.squared_length();
                        id = cnt;
                    }
                }else if(problem_t == 4){
                    auto v1 = s.to_vector();
                    double tmp = sqrt(v1.squared_length());
                    if(tmp < dis){
                        dis = tmp; id = cnt;
                    }
                }
            }
            cnt ++;
        }
        if(problem_t==2) {
            if(id != -1){
                covered_col_cloud->push_back(PointXYZRGB(x1, y1, z1, R[id], G[id], B[id]));
            }
        }else if(problem_t == 3){
            if(quality > thresh){
                covered_col_cloud->push_back(PointXYZRGB(x1, y1, z1, R[id], G[id], B[id]));
            }
        }else if(problem_t == 4){
            if(dis > thresh){
                covered_col_cloud->push_back(PointXYZRGB(x1, y1, z1, R[id], G[id], B[id]));
            }
        }
    }
    if(pcl::io::savePCDFile<pcl::PointXYZRGB>("covered_col_cloud_t" + to_string(problem_t) + ".pcd", *covered_col_cloud) == -1){
        cerr<<"Error saving covered cloud"<<endl; exit(1);
    }
}