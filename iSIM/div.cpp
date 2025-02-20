#include "div.h"

int get_new_index_n(Eigen::ArrayXXf total_data, Eigen::ArrayXf selected_condensed, int n, std::vector<int> select_from_n, std::string n_ary){
    int index = total_data.rows();
    double min_value = 3.08;
    for (int i : select_from_n){
        
        Eigen::ArrayXf fp = total_data.row(i);
        
        Eigen::ArrayXf col_sum = selected_condensed + fp;
        double sim = gen_sim_dict(col_sum, n+1, 1)[n_ary];
        if (sim < min_value){
            index = i;
            min_value = sim;
        }
    }
    return index;
}

int get_new_index_sqrt(Eigen::ArrayXXf total_data, Eigen::ArrayXf selected_condensed, int n, int k, std::vector<int> select_from_n, std::string n_ary){
    int index = total_data.rows();
    double min_value = 3.08;
    for (int i : select_from_n){
        Eigen::ArrayXf fp = total_data.row(i);
        Eigen::ArrayXf col_sum = selected_condensed + fp;
        double sim = gen_sim_dict(col_sum, n+1, k)[n_ary];
        if (sim < min_value){
            index = i;
            min_value = sim;
        }
    }
    return index;
}

int get_new_index_reverse(Eigen::ArrayXXf total_data, Eigen::ArrayXf selected_condensed, int n, std::vector<int> select_from_n, std::string n_ary){
    int index = total_data.rows();
    double min_value = 3.08;
    for (int i : select_from_n){
        Eigen::ArrayXf fp = total_data.row(i);
        Eigen::ArrayXf col_sum = selected_condensed - fp;
        double sim = gen_sim_dict(col_sum, n-1, 1)[n_ary];
        if (sim < min_value){
            index = i;
            min_value = sim;
        }
    }
    return index;
}

int get_new_indices_b_max(Eigen::ArrayXXf total_data, std::vector<int> selected_b, std::vector<int> select_from_b, std::string n_ary){
    double min_sim = 3.08;
    int idx = -1;
    for (int i : select_from_b){
        std::vector<double> comps;
        double max_sim = -3.08;
        for (int j : selected_b){
            Eigen::ArrayXXf fps = total_data({i, j}, Eigen::all);
            double sim = gen_sim_dict(fps, 1)[n_ary];
            comps.push_back(sim);
            if (sim>max_sim){
                max_sim = sim;
            }
        }
        if (max_sim < min_sim){
            min_sim = max_sim;
            idx = i;
        }
    }
    if (idx == -1){
        throw std::invalid_argument("No valid index found");
    }
    return idx;
}

std::vector<int> diversity(Eigen::ArrayXXf data, double percentage, std::string start, std::string n_ary, std::string method, int k){
    int n_total = data.rows();
    std::vector<int> select_from_n(n_total);
    std::iota (select_from_n.begin(), select_from_n.end(), 0);
    std::vector<int> selected_n;

    if (start=="medoid"){
        int seed = calculate_medoid(data, n_ary);
        selected_n.push_back(seed);
    }
    else if(start=="random"){
        std::random_device rand_dev;
        std::mt19937 generator(rand_dev());
        std::uniform_int_distribution<int> distr(0, n_total-1);
        int seed = distr(generator);
        selected_n.push_back(seed);
    }
    else if(start=="outlier"){
        int seed = calculate_outlier(data, n_ary);
        selected_n.push_back(seed);
    }
    else{
        throw std::invalid_argument("Select a correct starting point: medoid, random, outlier. ");
    }

    int n_max = n_total * percentage/100;
    select_from_n.erase(std::find(select_from_n.begin(), select_from_n.end(), selected_n[0]));

    Eigen::ArrayXf selected_condensed = data.row(selected_n[0]);
    for (int i=1; i<n_max; i++){
        int n = selected_n.size();
        int new_idx;
        if (method == "isim"){
            new_idx = get_new_index_n(data, selected_condensed, n, select_from_n, n_ary);
        }
        else if (method == "bmax"){
            new_idx = get_new_indices_b_max(data, selected_n, select_from_n, n_ary);
        }
        else if(method == "power"){
            new_idx = get_new_index_sqrt(data, selected_condensed, n, k, select_from_n, n_ary);
        }
        Eigen::ArrayXf new_fp = data.row(new_idx);
        selected_condensed = selected_condensed + new_fp;
        selected_n.push_back(new_idx);
        select_from_n.erase(std::find(select_from_n.begin(), select_from_n.end(), new_idx));
    }
    return selected_n;
}

std::vector<int> reverse_diversity(Eigen::ArrayXXf data, double percentage, std::string n_ary){
    int n_total = data.rows();
    std::vector<int> select_from_n(n_total);
    std::iota (select_from_n.begin(), select_from_n.end(), 0);
    int n_max = n_total * percentage/100;
    std::vector<int> deselected;
    Eigen::ArrayXf selected_condensed = data.colwise().sum();
    int n = n_total;
    for (int i=0; i< n_total-n_max; i++){
        int new_idx = get_new_index_reverse(data, selected_condensed, n, select_from_n, n_ary);
        Eigen::ArrayXf fp = data.row(new_idx);
        selected_condensed = selected_condensed - fp;
        deselected.push_back(new_idx);
        n -= 1;
        select_from_n.erase(std::find(select_from_n.begin(), select_from_n.end(), new_idx));
    }
    return select_from_n;
}