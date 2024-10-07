#include "clustering.h"


struct processedFps{
    Eigen::ArrayXXf data;
    Eigen::ArrayXf n_objects;
    std::list<std::list<int>> obj_indices;
};

struct processedFps pre_process(Eigen::ArrayXXf data){
    processedFps processed_data;
    processed_data.data = data;
    processed_data.n_objects = Eigen::ArrayXf::Zero(data.rows())+1;
    std::list<std::list<int>> obj_indices;
    for (int i = 0; i < data.rows(); i++){
        obj_indices.push_back({i});
    }
    return processed_data;
}

std::list<int> max_indices(struct processedFps fps, std::string n_ary){
    int n_objects = fps.data.rows();
    double max_sim = -3.08;
    int max1 = n_objects;
    int max2 = n_objects;
    int cols = fps.data.cols();
    #pragma omp parallel for
    for(int i = 0; i < n_objects; i++){
        for(int j = 0; j < n_objects; j++){
            if(i != j){
                Eigen::ArrayXf fp1 = fps.data.row(i);
                Eigen::ArrayXf fp2 = fps.data.col(j);
                int n = fps.n_objects[i] + fps.n_objects[j];
                double sim = gen_sim_dict(fp1+fp2, n, 1.0)[n_ary];
                if (sim > max_sim){
                    max_sim = sim;
                    max1 = i;
                    max2 = j;
                }
            }
        }
    }
    return {max1, max2};
}