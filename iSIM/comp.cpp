#include "isim.h"
#include "comp.h"


Eigen::ArrayXf calculate_a(const Eigen::ArrayXf col_sum){
    return col_sum*(col_sum-1.0)/2.0;
}

/**
 * @brief Calculates the pairwise average similarity.
 * 
 * @param fingerprints A 2D array of fingerprints (n_fps x n_features).
 * @param n_ary The type of similarity metric to use. RR, JT, or SM.
 * @return double The pairwise average similarity.
 */
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

/**
 * @brief Calculates various counters based on the column sum.
 * 
 * @param col_sum The column-wise sum of the data.
 * @param n_objects The number of objects (rows).
 * @param k The exponent used in the calculations.
 * @return std::map<std::string, double> A map containing the calculated counters.
 */
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

/**
 * @brief Calculates various counters based on the data.
 * 
 * @param data A 2D array of fingerprints (n_fps x n_features).
 * @param k The exponent used in the calculations.
 * @return std::map<std::string, double> A map containing the calculated counters.
 */
 std::map<std::string, double> calculate_counters( const Eigen::ArrayXXf data,
 const int k){ // array of arrays sub array contains binary object
    int n_objects = data.rows(); // num rows
    Eigen::ArrayXf col_sum = data.colwise().sum();
    return calculate_counters(col_sum, n_objects, k);
}

/**
 * @brief Calculates the instant similarity value based on the column sum.
 * 
 * @param col_sum The column-wise sum of the data.
 * @param n_objects The number of objects (rows).
 * @param n_ary The type of similarity metric to use. RR, JT, or SM.
 * @return double The calculated iSIM value.
 */
double calculate_isim(const Eigen::ArrayXf col_sum, 
                     const int n_objects, std::string n_ary){
    if (n_ary == "RR"){
        double a = (col_sum.cast<double>()*(col_sum.cast<double>()-1.0)*0.5).sum();
        double p = n_objects*(n_objects-1.0)*col_sum.size()*0.5;
        return a/p;
    }
    else if(n_ary == "JT"){
        double a = (col_sum.cast<double>()*(col_sum.cast<double>()-1.0)*0.5).sum();
        Eigen::ArrayXf off_coincidence = n_objects - col_sum;
        double total_dis = (off_coincidence.cast<double>()*col_sum.cast<double>()).sum();
        return a/(a+total_dis);
    }
    else if(n_ary == "SM"){
        double a = (col_sum.cast<double>()*(col_sum.cast<double>()-1.0)*0.5).sum();
        Eigen::ArrayXf off_coincidence = n_objects - col_sum;
        double d = (off_coincidence.cast<double>()*(off_coincidence.cast<double>()-1.0)/2.0).sum();
        double p = n_objects*(n_objects-1.0)*col_sum.size()/2.0;
        return (a+d)/p;
    }
    else{
        throw std::invalid_argument("Invalid n_ary value.");
    }

}

/**
 * @brief Calculates the instant similarity value based on the column sum.
 * 
 * @param data  A 2D Eigen::ArrayXXf of fingerprints (n_fps x n_features).
 * @param n_ary The type of similarity metric to use. RR, JT, or SM.
 * @return double The calculated iSIM value.
 */
double calculate_isim(const Eigen::ArrayXXf data, std::string n_ary){
    int n_objects = data.rows(); // num rows
    Eigen::ArrayXf col_sum = data.colwise().sum();

    return calculate_isim(col_sum, n_objects, n_ary);
}


/**
 * @brief Calculate the complementary similarity matrix based on the given data and method.
 * 
 * @param data A 2D Eigen::ArrayXXf of fingerprints (n_fps x n_features).
 * @param n_ary The type of similarity metric to use. RR, JT, or SM.
 * @return Eigen::ArrayXd The competition similarity vector calculated based on the specified method.
 * 
 * @throws std::invalid_argument If an invalid `n_ary` value is provided.
 */
Eigen::ArrayXd calculate_comp_sim(const  Eigen::ArrayXXf data, const std::string n_ary){
    int n_objects = data.rows();
    int remaining_objects = n_objects - 1;
    Eigen::ArrayXf col_sum = data.colwise().sum();
    Eigen::ArrayXd comp_sims(n_objects);
    #pragma omp parallel for
    for (int i = 0; i < n_objects; i++){
        Eigen::ArrayXf current_row = data.row(i); // note: this needs to be declared here so that col_sum_current is also a vector of lenght n_features
        Eigen::ArrayXf col_sum_current = col_sum - current_row;
        comp_sims(i) = calculate_isim(col_sum_current, remaining_objects, n_ary);
    }
    return comp_sims;
}

/**
 * @brief Calculate the medoid of the given data based on the specified similarity metric.
 * 
 * @param data A 2D Eigen::ArrayXXf of fingerprints (n_fps x n_features).
 * @param n_ary The type of similarity metric to use. RR, JT, or SM.
 * @return int The index of the medoid in the data.
 */
int calculate_medoid(const Eigen::ArrayXXf data, const std::string n_ary){
    Eigen::ArrayXd comp_sims = calculate_comp_sim(data, n_ary);
    int pos;
    comp_sims.minCoeff(&pos);
    return pos;
}

/**
 * @brief Calculate the outlier of the given data based on the specified similarity metric.
 * 
 * @param data A 2D Eigen::ArrayXXf of fingerprints (n_fps x n_features).
 * @param n_ary The type of similarity metric to use. RR, JT, or SM.
 * @return int The index of the outlier in the data.
 */
int calculate_outlier(const Eigen::ArrayXXf data, const std::string n_ary){
    Eigen::ArrayXd comp_sims = calculate_comp_sim(data, n_ary);
    int pos;
    comp_sims.maxCoeff(&pos);
    return pos;
}

/**
 * @brief Generates a dictionary of similarity metrics based on the column sum.
 * 
 * This function calculates various similarity metrics using the provided column sum,
 * number of objects, and exponent. The metrics are stored in a map with their respective
 * names as keys.
 * 
 * @param col_sum The column-wise sum of the data.
 * @param n_objects The number of objects (rows).
 * @param k The exponent used in the calculations.
 * @return std::map<std::string, double> A map containing the calculated similarity metrics.
 */
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

/**
 * @brief Generates a dictionary of similarity metrics based on the data.
 * 
 * @param data A 2D Eigen::ArrayXXf of fingerprints (n_fps x n_features).
 * @param k The exponent used in the calculations.
 * @return std::map<std::string, double> A map containing the calculated similarity metrics.
 */
std::map<std::string, double> gen_sim_dict(const Eigen::ArrayXXf data, 
const int k){
    int n_objects = data.rows(); // num rows
    Eigen::ArrayXf col_sum = data.colwise().sum();
    return gen_sim_dict(col_sum, n_objects, k);
}
