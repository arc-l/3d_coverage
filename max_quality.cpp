#include "Solution.h"


void Solution::precompute_4_1(){
    //****************************************
    // precompute for max quality, distance and visiblity matter
    //****************************************
    FILE* fp = fopen("precomputevis_mq.bin", "wb");
    vector<vector<float>> inter(num_of_locations);
    fwrite(&num_of_locations, sizeof(int), 1, fp);
    #pragma omp parellel for
    for(int i=0;i<num_of_locations;i++){
        inter[i].resize(num_of_pts);
        int cnt = 0;
        double x0 = possible_locs[i].x, y0 = possible_locs[i].y, z0 = possible_locs[i].z;
        for(const auto& pt: *cloud){
            double x1 = pt.x, y1 = pt.y, z1 = pt.z;
            Segment s(Point(x1 * .99999 + x0 * .00001, y1 * .99999 + y0 * .00001, z1 * .99999 + z0 * .00001), 
                Point(x0 * .99999 + x1 * .00001, y0 * .99999 + y1 * .00001, z0 * .99999 + z1 * .00001));
            if(s.is_degenerate()) continue;
            if(tree.do_intersect(s)){
                inter[i][cnt] = 0.0;
            }else{
                auto v1 = s.to_vector();
                // inter[i][cnt] = max(0.0, (v1[0] * pt.normal_x + v1[1] * pt.normal_y + v1[2] * pt.normal_z) / 
                // 		sqrt( pt.normal_x * pt.normal_x  + pt.normal_y * pt.normal_y + pt.normal_z * pt.normal_z) /  pow(v1.squared_length(), 1.5));
                inter[i][cnt] = sqrt(v1.squared_length());
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

void Solution::compute_4(int num_of_sensors, double ratio, bool Coarse, int new_num_of_locations, int new_num_of_pts){
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
    auto start_time = std::chrono::high_resolution_clock::now();
    // cout<<"Starting initiating grb constraints and vars ..."<<endl;
    if(Coarse){
        sub_sample_coarse(new_num_of_locations, new_num_of_pts, inter);
    }
    double ma = maxz - minz + maxx - minx + maxy - miny;
    double mi = (maxz - minz)/2;
    
    
    while(ma - mi > mi * 0.01){
        GRBModel model = GRBModel(env);
        GRBVar(*ys) = new GRBVar[num_of_pts];
        GRBVar(*zs) = new GRBVar[num_of_locations];
        GRBLinExpr sum = int(ratio * num_of_pts);
        GRBLinExpr lsum = -num_of_sensors;
        for(int i=0;i<num_of_locations;i++){
            zs[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "z_"+to_string(i));
            lsum += zs[i];
        }
        for(int i=0;i<num_of_pts;i++){
            ys[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y_"+to_string(i));
            sum -= ys[i];
        }			
        model.addConstr(lsum <= 0, "sensor_constr");
        model.addConstr(sum <= 0, "sensor_constr");
        double mid = (ma + mi) / 2;
        for(int i=0;i<num_of_pts;i ++){
            GRBLinExpr expr = ys[i];
            for(int j=0;j<num_of_locations;j++){
                if(inter[j][i] > 1e-8 && inter[j][i] < mid)
                    expr -= zs[j];
            }
            model.addConstr(expr <= 0, "cover_" + to_string(i));
        }
        model.optimize();
        if(model.get(GRB_IntAttr_SolCount) > 0){
            loc_cloud -> clear();
            vector<int> locs;
            for(int i=0;i<num_of_locations;i++){
                if(zs[i].get(GRB_DoubleAttr_X) == 1.0){
                    loc_cloud->push_back(PointXYZRGB(possible_locs[i].x, possible_locs[i].y, possible_locs[i].z, 180, 0, 0));
                    locs.push_back(i);
                }
            }
            vector<double> distance;
            for(int i=0;i<num_of_pts;i++){
                float var = mid;
                for(int it: locs){
                    if(inter[it][i] > 1e-8 && inter[it][i] < mid) {
                        var = min(var, inter[it][i]);
                    }
                }
                if(var <= mid) distance.push_back(var);
            }
            nth_element(begin(distance), begin(distance) + int(ratio * num_of_pts) - 1, end(distance));
            ma = *(begin(distance) + int(ratio * num_of_pts) -1 );
        }else{
            mi = mid;
        }
        delete[] ys, zs;
    }
    auto stop_time = std::chrono::high_resolution_clock::now();
    cout << duration_cast<microseconds>(stop_time - start_time).count() << " microseconds"<<endl;
    cout << "R = " << ma<<endl;
    if(pcl::io::savePCDFile<pcl::PointXYZRGB>("locations_t1.5.pcd", *loc_cloud)==-1){
        cerr<<"error writing pcd file"<<endl; exit(1);
    }
}

void Solution::local_imprv_4(){
    if(pcl::io::loadPCDFile<pcl::PointXYZRGB> ("locations_t1.5.pcd", *loc_cloud) == -1){
        cerr<<"error loading location pcd file"<<endl; exit(1);
    }
    int num_of_sensors = loc_cloud->size();
    cout << "Input ratio: ";
    double rau;
    cin >> rau;
    auto start_time = chrono::high_resolution_clock::now();
    int tot_pts = rau * num_of_pts;
    vector<vector<vec3>> gp(num_of_sensors);
    vector<vec3> cent_loc;
    for(int i=0;i<num_of_sensors;i++){
        const auto & pt = (*loc_cloud)[i];
        cent_loc.emplace_back(pt.x, pt.y, pt.z);
    }
    vector<pair<double,pair<int,int>>> dist;
    vector<int> par(cloud->size(), -1);
    double radius = 0, initial, end_criteria;
    int round_num = 1;
    function<double(const vector<vec3> &, vec3 &)> func = [&](const vector<vec3> &pts, vec3 &loc){
        double step = initial;
        double x = loc.x, y = loc.y, z=loc.z, nx, ny, nz=loc.z;
        // compute distance
        double r = 0;
        for(auto &pt: pts) r = max(r, dis(pt, loc));
        while(step > end_criteria){
            double dx[4] = {step, -step, 0, 0};
            double dy[4] = {0, 0, step, -step};
            //cout<<"step  = " <<step << "=======   quality = "<<r<<endl;
            while(1){
                bool t = false;
                for(int i=0;i<4;i++){
                    double nx = x + dx[i];
                    double ny = y + dy[i];
                    bool succ = true;
                    double len = 0;
                    for(auto &it: pts){
                        double x0 = it.x, y0 = it.y, z0 = it.z;
                        Segment s(Point(nx * .99999 + x0 * .00001, ny * .99999 + y0 * .00001, nz * .99999 + z0 * .00001), 
                            Point(x0 * .99999 + nx * .00001, y0 * .99999 + ny * .00001, z0 * .99999 + nz * .00001));
                        if(s.is_degenerate()) continue;
                        if(tree.do_intersect(s)){
                            succ = false; break;
                        }
                        len = max(len, 1 / 0.99998 * sqrt(s.squared_length()));
                        if(len > r) {succ=false;break;}
                    }
                    if(succ){
                        x = nx, y = ny, r = len;
                        t = true;
                        break;
                    }
                }
                if(!t) break;
            }
            step /= 2;
        }
        loc = vec3(x, y, z);
        return r;
    };
    while(1){
        for(auto &v: gp) v.clear();
        dist.clear();
        int pid = 0;
        for(const auto & pt: *cloud){
            double x1 = pt.x, y1 = pt.y, z1 = pt.z, dis = std::numeric_limits<double>::max();
            int cnt = 0, id = -1;
            for(const auto &pt2: cent_loc){
                double x0 = pt2.x, y0 = pt2.y, z0 = pt2.z;
                Segment s(Point(x1 * .99999 + x0 * .00001, y1 * .99999 + y0 * .00001, z1 * .99999 + z0 * .00001), 
                    Point(x0 * .99999 + x1 * .00001, y0 * .99999 + y1 * .00001, z0 * .99999 + z1 * .00001));
                if(s.is_degenerate()) continue;
                if(!tree.do_intersect(s)){
                    auto v1 = s.to_vector();
                    double tmp = sqrt(s.squared_length());
                    if(dis > tmp) {
                        dis = tmp; id = cnt;
                    }
                }
                cnt ++;
            }
            if(id >=0 ){
                dist.emplace_back(dis, make_pair(id, pid));
                par[pid] = id;
            }
            pid ++;
        }
        nth_element(begin(dist), begin(dist) + tot_pts - 1, end(dist));
        radius = dist[tot_pts-1].first * 1 / 0.99998;
        // if(round_num++ == 1) 
        //     cout<<"Initial Radius = "<<radius<<endl;

        initial = radius / 20;
        end_criteria = radius / 10000;

        for(int i=0; i<tot_pts; i++){
            const auto& pt = (*cloud)[dist[i].second.second];
            gp[dist[i].second.first].emplace_back(pt.x, pt.y, pt.z);
        }
        
        double nr = 0;
        for(int i=0;i<num_of_sensors;i++) {nr = max(nr, func(gp[i], cent_loc[i]));}
        if( radius - nr < 1e-8) break;
        radius = nr;
    }
    auto stop_time = chrono::high_resolution_clock::now();
    cout<<"Result radius: "<<radius<<endl;
    cout << duration_cast<microseconds>(stop_time - start_time).count() << " microseconds"<<endl;
    loc_cloud->clear();
    for(int i=0;i<num_of_sensors;i++){
        loc_cloud->push_back(PointXYZRGB((float)cent_loc[i].x, (float)cent_loc[i].y, (float)cent_loc[i].z, (uint8_t)180, (uint8_t)0, (uint8_t)0));
    }
    if(pcl::io::savePCDFile<pcl::PointXYZRGB>("locations_t1.5.pcd", *loc_cloud) == -1){
        cerr<<"Error saving covered cloud"<<endl; exit(1);
    }
}