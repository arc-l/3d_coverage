#include "Solution.h"

void Solution::precompute_2_1(){
    
    //***********************************
    // precompute for bare visiblity model
    // saving into precomputevis.bin
    //***********************************
    
    FILE *fp = fopen("precomputevis.bin", "wb");
    vector<vector<unsigned char>> to_write(num_of_locations);
    fwrite(&num_of_locations, sizeof(int), 1, fp);
    // Dont know why openmp is not working here
    // #pragma omp parellel for
    for(int i = 0; i < num_of_locations; i++){
        vector<bool> inter(num_of_pts, false);
        int cnt = 0;
        double x0 = possible_locs[i].x, y0 = possible_locs[i].y, z0 = possible_locs[i].z;
        for(const auto& pt: *cloud){
            double x1 = pt.x, y1 = pt.y, z1 = pt.z;
            Segment s(Point(x1 * .99999 + x0 * .00001, y1 * .99999 + y0 * .00001, z1 * .99999 + z0 * .00001), 
                Point(x0 * .99999 + x1 * .00001, y0 * .99999 + y1 * .00001, z0 * .99999 + z1 * .00001));
            if(s.is_degenerate()) continue;
            if(!tree.do_intersect(s)){
                inter[cnt] = true;
            }
            cnt++;
        }
        for(int k=0;k<num_of_pts;k+=8){
            unsigned char tmp = 0;
            for(int inc = 0; inc < 8 && k + inc < num_of_pts; inc++){
                tmp <<= 1;
                tmp |= inter[k + inc];
            }
            to_write[i].push_back(tmp);
        }
    }
    for(int i=0; i < num_of_locations; i++){
        fwrite(&possible_locs[i].x, sizeof(double), 1, fp);
        fwrite(&possible_locs[i].y, sizeof(double), 1, fp);
        fwrite(&possible_locs[i].z, sizeof(double), 1, fp);
        fwrite(&num_of_pts, sizeof(int), 1, fp);
        fwrite(to_write[i].data(), 1, to_write[i].size(), fp);
    }
    fclose(fp);
}

void Solution::compute_2(int num_of_sensors){
    //******************************************
    // Compute for bare visibility
    //******************************************
    
    possible_locs.clear();
    FILE* fp = fopen("precomputevis.bin", "rb");
    double x0, y0, z0;
    fread(&num_of_locations, 1, sizeof(int), fp);
    
    // Load precomputed results into ``inter''
    
    vector<vector<bool>> inter(num_of_sensors);
    for(int i=0;i<num_of_sensors;i++){
        fread(&x0, sizeof(double), 1, fp);
        fread(&y0, sizeof(double), 1, fp);
        fread(&z0, sizeof(double), 1, fp);
        fread(&num_of_pts, sizeof(int), 1, fp);
        possible_locs.emplace_back(x0, y0, z0);
        inter[i].resize(num_of_pts);
        for(int k=0;k<num_of_pts;k+=8){
            unsigned char c;
            fread(&c, 1, 1, fp);
            for(int inc = min(7, num_of_pts - k -1);~inc;inc--){
                inter[i][k + inc] = c & 1;
                c >>= 1;
            }
        }
    }
    fclose(fp);

    auto start_time = std::chrono::high_resolution_clock::now();
    // cout<<"Start initiating grb constraints and vars ..."<<endl;
    GRBModel model = GRBModel(env);
    GRBVar(*ys) = new GRBVar[num_of_pts];
    GRBVar(*zs) = new GRBVar[num_of_sensors];

    GRBLinExpr lsum = - num_of_sensors;
    for(int j=0;j<num_of_locations;j++){
        zs[j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "z_" + to_string(j));
        lsum += zs[j];
    }
    for(int i=0;i<num_of_pts;i++){
        ys[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y_" + to_string(i));
    }
    model.addConstr(lsum <= 0, "sensor_constr");
    GRBLinExpr sum = 0;
    for(int i=0;i<num_of_pts;i++){
        GRBLinExpr expr = ys[i];
        for(int j=0;j<num_of_locations;j++){
            if(inter[j][i]) expr -= zs[j];
        }
        model.addConstr(expr <= 0, "cover_" + to_string(i));
        sum += ys[i];
    }
    // cout<<"Finished setting up, start solving"<<endl;
    model.setObjective(sum, GRB_MAXIMIZE);
    model.optimize();
    auto stop_time = std::chrono::high_resolution_clock::now();
    cout << duration_cast<microseconds>(stop_time - start_time).count() << " microseconds"<<endl;
    cout<<"[";
    int cnt = 0;
    for(int i=0;i<num_of_locations;i++){
        if(zs[i].get(GRB_DoubleAttr_X) == 1.0){
            // cout<<"["<<possible_locs[i].x <<", "<<possible_locs[i].y<<", "<<possible_locs[i].z<<"], ";
            loc_cloud->push_back(PointXYZRGB(possible_locs[i].x, possible_locs[i].y, possible_locs[i].z, 180, 0, 0));
            cnt ++;
        }
    }
    cout<<"]"<<endl;
    int tot_p = 0;
    for(int i=0;i<num_of_pts;i++){
        if(ys[i].get(GRB_DoubleAttr_X) == 1.0){
            tot_p ++;
        }
    }
    cout<<"Total points "<<tot_p<<endl;
    delete[] ys, zs; 
    if(pcl::io::savePCDFile<pcl::PointXYZRGB>("locations_t1.pcd", *loc_cloud)==-1){
        cerr<<"error writing pcd file"<<endl; exit(1);
    }
}

