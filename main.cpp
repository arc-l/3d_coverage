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

#include"Vec.h"
#include"Solution.h"



using namespace std;
using namespace std::chrono; 


int main(int argc, char*argv[]){
	if(argc<3){
		printf("Format:\n\t"
			"2 0 preprocess piont close (only visible points from the ceiling will remain) \n\t"
			"2 1 precompute vis\n\t"
			"2 2 compute vis ilp model\n\t"
			"3 1 precompute cumulative model\n\t"
			"3 2 compute cumulative model\n\t"
			"3 2 1 compute cumulative model (coarse)\n\t"
			"4 1 precompute mq model\n\t"
			"4 2 compute mq model\n\t"
			"4 2 1 compute mq model (coarse)\n\t"
			"5 1 locally improve max quality model\n\t"
			"5 2 locally improve cumulative quality model\n"
			);
		return 0;
	}
	if(argv[1][0] == '2'){
		Solution S;
		if(argv[2][0] == '0'){
			S.build_AABB_Tree();
			S.sub_sample();
		}else if(argv[2][0] == '1'){
			S.build_AABB_Tree();
			S.comp_possible_sensor_locations(600);
			S.precompute_2_1();
		}else if(argv[2][0] == '2'){
			cout<<"Input num of sensors: ";
			int input_num;
			cin >> input_num;
			S.build_AABB_Tree();
			S.compute_2(input_num);
		}else{ // visualize
			S.build_AABB_Tree();
			S.visualize(2);
		}
		return 0;
	}
	else if(argv[1][0] == '3'){
		Solution S;
		if(argv[2][0] == '1'){
			S.build_AABB_Tree();
			S.comp_possible_sensor_locations(200);
			S.precompute_3_1();
		}else if(argv[2][0] == '2'){
			cout<<"Input num of sensors: ";
			int input_num;
			cin>>input_num;
			S.build_AABB_Tree();
			// set default threshold
			double threshold = 1 / pow(S.maxz - S.minz, 2) * 3. * sqrt(3) / 8.;
			if(argc > 3){
				S.compute_3(input_num, threshold = threshold, true, /*new_num_of_pts =*/ 1000, /*new_num_of_locaitons =*/ 100);
			}
			else S.compute_3(input_num, threshold = threshold);
		}else if(argv[2][0] == '3'){
			S.build_AABB_Tree();
			double threshold = 1 / pow(S.maxz - S.minz, 2) * 3. * sqrt(3) / 8.;
			S.visualize(3, threshold);
		}
	}
	else if(argv[1][0] == '4'){
		Solution S;
		if(argv[2][0] == '1'){
			S.build_AABB_Tree();
			S.comp_possible_sensor_locations(200);
			S.precompute_4_1();
		}else if(argv[2][0] == '2'){
			int input_num;
			double ratio;
			cout<<"Input num of sensors: ";
			cin>>input_num;
			cout<<"Input num of ratio: ";
			cin>>ratio;
			S.build_AABB_Tree();
			if(argc > 3){
				S.compute_4(input_num, ratio, true,  1000, 100);
			}
			else S.compute_4(input_num, ratio = ratio);
		}else if(argv[2][0] == '3'){
			double radius;
			cin >> radius;
			S.build_AABB_Tree();
			S.visualize(3, radius);
		}
	}else if(argv[1][0] == '5'){
		// Locally improve a solution
		if(argv[2][0] == '1'){
			// min radius model
			// load locations_t1.5.pcd
			Solution S;
			S.build_AABB_Tree();
			S.local_imprv_4();
		}
		else if(argv[2][0] == '2'){
			// Locally improve cumulative exposure model
			Solution S;
			S.build_AABB_Tree();
			S.local_imprv_3();	
		}
	}
	
	return 0;
}

