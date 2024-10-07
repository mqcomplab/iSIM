#include "isim.h"
#include "comp.h"


Eigen::ArrayXf calculate_a(const Eigen::ArrayXf col_sum){
    return col_sum*(col_sum-1.0)/2.0;
}

double pairwise_average(const Eigen::ArrayXXf fingerprints, const std::string n_ary){
    double avg = 0;
    int n_objects = fingerprints.rows(); // num rows
    #pragma omp parallel for reduction(+:avg)
    for (int i =0; i< n_objects - 1; i++){
        for (int j = i+1; j < n_objects; j++){
            avg += calculate_isim(fingerprints.row(i)+fingerprints.row(j), 2, n_ary);
        }
    }
    return avg/(n_objects*(n_objects-1.0)/2.0);
}

std::map<std::string, double> calculate_counters(const Eigen::ArrayXf col_sum, 
const int n_objects, const int k) { // data is column wise sum
    std::map<std::string, double> counters;
    Eigen::ArrayXf a_array = calculate_a(col_sum);
    Eigen::ArrayXf off_coincidence = n_objects - col_sum;
    Eigen::ArrayXf d_array = calculate_a(off_coincidence);
    Eigen::ArrayXf dis_array = off_coincidence*col_sum;
    double a = a_array.pow(1./k).sum();
    double d = d_array.pow(1./k).sum();
    double total_dis = dis_array.pow(1./k).sum();
    double total_sim = a + d;
    double p = total_sim + total_dis;

    counters["a"] = a;
    counters["d"] = d;
    counters["total_dis"] = total_dis;
    counters["total_sim"] = total_sim;
    counters["p"] = p;
    return counters;
}

 std::map<std::string, double> calculate_counters( const Eigen::ArrayXXf data,
 const int k){ // array of arrays sub array contains binary object
    int n_objects = data.rows(); // num rows
    Eigen::ArrayXf col_sum = data.colwise().sum();
    return calculate_counters(col_sum, n_objects, k);
}

double calculate_isim(const Eigen::ArrayXf col_sum, 
                     const int n_objects, std::string n_ary){
    if (n_ary == "RR"){
        double a = (col_sum*(col_sum-1.0)*0.5).sum();
        double p = n_objects*(n_objects-1.0)*col_sum.size()*0.5;
        return a/p;
    }
    else if(n_ary == "JT"){
        double a = (col_sum*(col_sum-1.0)/2.0).sum();
        Eigen::ArrayXf off_coincidence = n_objects - col_sum;
        double total_dis = (off_coincidence*col_sum).sum();
        return a/(a+total_dis);
    }
    else if(n_ary == "SM"){
        double a = (col_sum*(col_sum-1.0)/2.0).sum();
        Eigen::ArrayXf off_coincidence = n_objects - col_sum;
        double d = (off_coincidence*(off_coincidence-1.0)/2.0).sum();
        double p = n_objects*(n_objects-1.0)*col_sum.size()/2.0;
        return (a+d)/p;
    }
    else{
        throw std::invalid_argument("Invalid n_ary value.");
    }

}

double calculate_isim(const Eigen::ArrayXXf data, std::string n_ary){
    int n_objects = data.rows(); // num rows
    Eigen::ArrayXf col_sum = data.colwise().sum();

    return calculate_isim(col_sum, n_objects, n_ary);
}


Eigen::ArrayXf calculate_comp_sim(const  Eigen::ArrayXXf data, const std::string n_ary){
    int n_objects = data.rows() - 1; 
    Eigen::ArrayXf col_sum = data.colwise().sum();
    Eigen::ArrayXXf col_sum_m = col_sum.transpose().replicate(n_objects+1, 1);
    Eigen::ArrayXXf comp_matrix = col_sum_m - data;
    Eigen::ArrayXXf a = comp_matrix*(comp_matrix-1.0)/2.0;
    int m = data.cols();
    
    if (n_ary == "RR"){
        Eigen::ArrayXf comp_sim = a.rowwise().sum()/(m * n_objects * (n_objects-1.0)/2.0);
        return comp_sim;
    }
    else if (n_ary == "JT"){
        Eigen::ArrayXf comp_sims = a.rowwise().sum();
        comp_sims /= ((a+ comp_matrix*(n_objects - comp_matrix)).rowwise().sum());
        return comp_sims;
    }
    else if (n_ary == "SM"){
        Eigen::ArrayXf comp_sims = (a + (n_objects - comp_matrix) * (n_objects - comp_matrix - 1.0)/2.0).rowwise().sum();
        comp_sims /= (m * n_objects * (n_objects - 1.0)/2.0);
        return comp_sims;
    }
    else
        throw std::invalid_argument("Invalid n_ary value.");
}

int calculate_medoid(const Eigen::ArrayXXf data, const std::string n_ary){
    Eigen::ArrayXf comp_sims = calculate_comp_sim(data, n_ary);
    int pos;
    float min_sim = comp_sims.minCoeff(&pos);
    return pos;
}

int calculate_outlier(const Eigen::ArrayXXf data, const std::string n_ary){
    Eigen::ArrayXf comp_sims = calculate_comp_sim(data, n_ary);
    int pos;
    float max_sim = comp_sims.maxCoeff(&pos);
    return pos;
}


std::map<std::string, double> gen_sim_dict(const Eigen::ArrayXf col_sum,
const int n_objects, const int k){
    std::map<std::string, double> counters = calculate_counters(col_sum, n_objects, k);
    std::map<std::string, double> sim_dict;
    sim_dict["AC"] = asin(sqrt(counters["total_sim"]/counters["p"]))*2.0/M_PI;
    sim_dict["BUB"] = (sqrt(counters["a"] * counters["d"])) + counters["a"];
    sim_dict["BUB"] /= (sqrt(counters["a"] * counters["d"])) + counters["a"] + counters["total_dis"];
    sim_dict["Fai"] = (counters["a"] + 0.5 * counters["d"])/counters["p"];
    sim_dict["Gle"] = 2.0*counters["a"]/(2.0 * counters["a"] + counters["total_dis"]);
    sim_dict["Ja"] = 3.0 * counters["a"]/(3.0 * counters["a"] + counters["total_dis"]);
    sim_dict["JT"] = counters["a"]/(counters["a"] + counters["total_dis"]);
    sim_dict["RT"] = counters["total_sim"]/(counters["p"] + counters["total_dis"]);
    sim_dict["RR"] = counters["a"]/counters["p"];
    sim_dict["SM"] = counters["total_sim"]/counters["p"];
    sim_dict["SS1"] = counters["a"]/(counters["a"] + 2.0*counters["total_dis"]);
    sim_dict["SS2"] = 2.0*counters["total_sim"]/(counters["p"] + counters["total_sim"]);
    return sim_dict;
}

std::map<std::string, double> gen_sim_dict(const Eigen::ArrayXXf data, 
const int k){
    int n_objects = data.rows(); // num rows
    Eigen::ArrayXf col_sum = data.colwise().sum();
    return gen_sim_dict(col_sum, n_objects, k);
}
