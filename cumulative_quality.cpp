#include "Solution.h"

void Solution::precompute_3_1(){
    //*********************************************
    // precompute for quality with respect to max(0, cos / r^2)
    //*********************************************
    FILE* fp = fopen("precomputevis_plus.bin", "wb");
    vector<vector<float>> inter(num_of_locations);
    fwrite(&num_of_locations, sizeof(int), 1, fp);
    #pragma omp parellel for
    for(int i=0;i<num_of_locations;i++){
        inter[i].resize(num_of_pts);
        int cnt = 0;
        float x0 = possible_locs[i].x, y0 = possible_locs[i].y, z0 = possible_locs[i].z;
        for(const auto& pt: *cloud){
            float x1 = pt.x, y1 = pt.y, z1 = pt.z;
            Segment s(Point(x1 * .99999 + x0 * .00001, y1 * .99999 + y0 * .00001, z1 * .99999 + z0 * .00001), 
                Point(x0 * .99999 + x1 * .00001, y0 * .99999 + y1 * .00001, z0 * .99999 + z1 * .00001));
            if(s.is_degenerate()) continue;
            if(tree.do_intersect(s)){
                inter[i][cnt] = 0.0;
            }else{
                auto v1 = s.to_vector();
                inter[i][cnt] = max(0.0, (v1[0] * pt.normal_x + v1[1] * pt.normal_y + v1[2] * pt.normal_z) / 
                        sqrt( pt.normal_x * pt.normal_x  + pt.normal_y * pt.normal_y + pt.normal_z * pt.normal_z) /  pow(v1.squared_length(), 1.5));
            }
            cnt++;
        }
    }
    for(int i=0;i<num_of_locations;i++){
        fwrite(&possible_locs[i].x, sizeof(possible_locs[i].x), 1, fp);
        fwrite(&possible_locs[i].y, sizeof(possible_locs[i].y), 1, fp);
        fwrite(&possible_locs[i].z, sizeof(possible_locs[i].z), 1, fp);
        fwrite(&num_of_pts, sizeof(int), 1, fp);
        fwrite(inter[i].data(), sizeof(float), num_of_pts, fp);
    }
    fclose(fp);	
}

void Solution::compute_3(int num_of_sensors, double thresh, bool Coarse, int new_num_of_locations, int new_num_of_pts){
    //****************************************
    // Compute for cumulative quality model
    //****************************************
    FILE* fp = fopen("precomputevis_plus.bin", "rb");
    fread(&num_of_locations, 1, sizeof(int), fp);
    vector<vector<float>> inter(num_of_locations);
    double x0, y0, z0;
    possible_locs.clear();
    for(int i=0;i<num_of_locations;i++){
        fread(&x0, sizeof(double), 1, fp);
        fread(&y0, sizeof(double), 1, fp);
        fread(&z0, sizeof(double), 1, fp);
        fread(&num_of_pts, sizeof(double), 1, fp);
        possible_locs.emplace_back(x0, y0, z0);
        inter[i].resize(num_of_pts);
        fread(inter[i].data(), sizeof(float), num_of_pts, fp);
    }
    fclose(fp);
    cout<<"Starting initiating grb constraints and vars ..."<<endl;
    if(Coarse){
        sub_sample_coarse(new_num_of_locations, new_num_of_pts, inter);
    }
    GRBModel model = GRBModel(env);
    // if(argc <= 3)model.set("TimeLimit", "600.0");
    // else model.set("TimeLimit", "60.0");
    // model.set("MIPGap", "0.1");
    GRBVar *ys = new GRBVar[num_of_pts];
    GRBVar *zs = new GRBVar[num_of_locations];
    auto start_time = std::chrono::high_resolution_clock::now();
    GRBLinExpr lsum = -num_of_sensors;
    for(int i=0;i<num_of_pts;i++){
        ys[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y_"+to_string(i));
    }
    for(int i=0;i<num_of_locations;i++){
        zs[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "z_"+to_string(i));
        lsum += zs[i];
    }
    model.addConstr(lsum <= 0, "sensor_constr");
    GRBLinExpr sum = 0;
    for(int i=0;i<num_of_pts;i++){
        GRBLinExpr expr = ys[i] * 10;
        for(int j=0;j<num_of_locations;j++){
            if(inter[j][i] / thresh > 0) expr -= zs[j] * inter[j][i] / thresh * 10;
        }
        model.addConstr(expr <= 0, "cover_" + to_string(i));
        sum += ys[i];
    }
    cout<<"Finished setting up, start solving"<<endl;
    model.setObjective(sum, GRB_MAXIMIZE);
    model.optimize();
    auto stop_time = std::chrono::high_resolution_clock::now();
    cout << duration_cast<microseconds>(stop_time - start_time).count() << " microseconds"<<endl;
    // cout<<"[";
    int cnt = 0;
    for(int i=0;i<num_of_locations;i++){
        if(zs[i].get(GRB_DoubleAttr_X) == 1.0){
            // cout<<"["<<xval[i] <<", "<<yval[j]<<", "<<z0<<"], ";
            loc_cloud->push_back(PointXYZRGB(possible_locs[i].x, possible_locs[i].y, 
                                                    possible_locs[i].z, 180, 0, 0));
            cnt ++;
        }
    }
    // cout<<"]"<<endl;
    int tot_p = 0;
    for(int i=0;i<num_of_pts;i++){
        if(ys[i].get(GRB_DoubleAttr_X) == 1.0){
            tot_p ++;
        }
    }
    // cout<<"Total points "<<tot_p<<endl;
    delete[] ys, zs;
    if(pcl::io::savePCDFile<pcl::PointXYZRGB>("locations_t2.pcd", *loc_cloud)==-1){
        cerr<<"error writing pcd file"<<endl; exit(1);
    }
}

void Solution::local_imprv_3(){
    if(pcl::io::loadPCDFile<pcl::PointXYZRGB> ("locations_t2.pcd", *loc_cloud) == -1){
        cerr<<"error loading location pcd file"<<endl; exit(1);
    }
    auto start_time = chrono::high_resolution_clock::now();
    int num_of_sensors = loc_cloud->size();
    vector<vec3> cent_loc;
    // vector<vec3> out_cloud;

    for(const auto& pt: *loc_cloud){
        cent_loc.emplace_back(pt.x, pt.y, pt.z);
    }
    double thresh;
    thresh = 1.0 / (maxz - minz)  / (maxz - minz) * sqrt(3) / 2 * 3 / 4;
    int round_num = 1;
    vector<vector<float>> contrib(num_of_sensors, vector<float>(num_of_pts));
    int pid = 0;
    for(const auto & pt: *cloud){
        double x1 = pt.x, y1 = pt.y, z1 = pt.z, dis = std::numeric_limits<double>::max();
        int cnt = 0;
        double sum = 0;
        for(const auto &pt2: cent_loc){
            double x0 = pt2.x, y0 = pt2.y, z0 = pt2.z;
            Segment s(Point(x1 * .99999 + x0 * .00001, y1 * .99999 + y0 * .00001, z1 * .99999 + z0 * .00001), 
                Point(x0 * .99999 + x1 * .00001, y0 * .99999 + y1 * .00001, z0 * .99999 + z1 * .00001));
            if(s.is_degenerate()) continue;
            if(!tree.do_intersect(s)){
                auto v1 = s.to_vector();
                // double tmp = int(max(0.0, (v1[0] * pt.normal_x + v1[1] * pt.normal_y + v1[2] * pt.normal_z) / 
                // 				sqrt( pt.normal_x * pt.normal_x  + pt.normal_y * pt.normal_y + pt.normal_z * pt.normal_z) /  pow(v1.squared_length(), 1.5)) / thresh * 10);
                double tmp = (max(0.0, (v1[0] * pt.normal_x + v1[1] * pt.normal_y + v1[2] * pt.normal_z) / 
                                sqrt( pt.normal_x * pt.normal_x  + pt.normal_y * pt.normal_y + pt.normal_z * pt.normal_z) /  pow(v1.squared_length(), 1.5)) / thresh * 10);

                if(tmp>1e-8){
                    contrib[cnt][pid] = tmp;
                    sum += tmp;
                }
            }
            cnt ++;
        }
        pid ++;
    }

    vector<float> cumu(num_of_pts, 0);
    vector<bool> is_covered(num_of_pts, false);
    int res = 0;
    for(int i=0; i<num_of_pts; i++){
        // cumu[i].second = i;
        for(int j=0; j<num_of_sensors;j++){
            cumu[i] += contrib[j][i];
        }
        is_covered[i] = (cumu[i] >= 10);
        if(is_covered[i]) res ++;
    }
    cout<<"Initially Covered "<< res<<endl;

    double initial =  (maxz - minz) / 20;
    double end_criteria = (maxz - minz) / 10000;
    function<bool(const int,  vec3 &)> func = [&](const int s_id, vec3 &loc){
        double step = initial;
        double x = loc.x, y = loc.y, z=loc.z, nx, ny, nz=loc.z;
        bool c = false;
        while(step > end_criteria){
            double dx[4] = {step, -step, 0, 0};
            double dy[4] = {0, 0, step, -step};
            // cout<<"step  = " <<step << "=======   quality = "<<r<<endl;
            while(1){
                bool t = false;
                auto pre = contrib[s_id];
                auto p_covered = is_covered;
                auto p_cumu = cumu;
                for(int i=0;i<4;i++){
                    double nx = x + dx[i];
                    double ny = y + dy[i];
                    // bool succ = true;
                    // double len = 0;
                    int gain = 0;
                    int p_id = 0;
                    for(auto &pt: *cloud){
                        double x0 = pt.x, y0 = pt.y, z0 = pt.z;
                        Segment s(Point(nx * .99999 + x0 * .00001, ny * .99999 + y0 * .00001, nz * .99999 + z0 * .00001), 
                            Point(x0 * .99999 + nx * .00001, y0 * .99999 + ny * .00001, z0 * .99999 + nz * .00001));
                        if(s.is_degenerate()) continue;
                        if(!tree.do_intersect(s)){
                            auto v1 = s.to_vector();
                            // contrib[s_id][p_id] = int(max(0.0, (v1[0] * pt.normal_x + v1[1] * pt.normal_y + v1[2] * pt.normal_z) / 
                            // 			sqrt( pt.normal_x * pt.normal_x  + pt.normal_y * pt.normal_y + pt.normal_z * pt.normal_z) /  pow(v1.squared_length(), 1.5)) / thresh * 10);
                            contrib[s_id][p_id] = (max(0.0, (v1[0] * pt.normal_x + v1[1] * pt.normal_y + v1[2] * pt.normal_z) / 
                                        sqrt( pt.normal_x * pt.normal_x  + pt.normal_y * pt.normal_y + pt.normal_z * pt.normal_z) /  pow(v1.squared_length(), 1.5)) / thresh * 10);
                        }else{
                            contrib[s_id][p_id] = 0;
                        }
                        cumu[p_id] += contrib[s_id][p_id] - pre[i];
                        is_covered[p_id] = (cumu[p_id] >= 10);
                        if(is_covered[p_id] ^ p_covered[p_id]) gain += (is_covered[p_id] ? 1: -1);
                        p_id ++;
                    }
                    if(gain > 0){
                        res += gain;
                        cout<<res<<endl;
                        x = nx; y= ny;
                        t = true; c = true;
                        break;
                    }else{
                        (contrib[s_id] = pre);
                        (is_covered = p_covered);
                        (cumu = p_cumu);
                    }
                    // if(succ){
                    // 	x = nx, y = ny, r = len;
                    // 	t = true;
                    // 	break;
                    // }
                }
                if(!t) break;
            }
            step /= 2;
        }
        loc = vec3(x, y, z);
        return c;
    };
    while(1){
        bool changed = false;
        for(int i=0;i<num_of_sensors;i++) changed |= func(i, cent_loc[i]);
        if(!changed) break;
    }
    auto stop_time = chrono::high_resolution_clock::now();
    cout<<duration_cast<microseconds>(stop_time - start_time).count() << " microseconds" <<endl;
    loc_cloud->clear();
    for(int i=0;i<num_of_sensors;i++){
        loc_cloud->push_back(PointXYZRGB(cent_loc[i].x, cent_loc[i].y, cent_loc[i].z, 180, 0, 0));
    }
    if(pcl::io::savePCDFile<pcl::PointXYZRGB>("locations_t2.pcd", *loc_cloud) == -1){
        cout<<"Error saving pcd file"<<endl;
        exit(1);
    }
}