#include "real.h"

/**
 * @brief Calculate the pairwise similarity between two fingerprints using RR methods.
 * @param fp1 The first fingerprint as an Eigen::ArrayXf.
 * @param fp2 The second fingerprint as an Eigen::ArrayXf.
 * @return The pairwise similarity as a double.
 */
double pair_rr(const Eigen::ArrayXf fp1, const Eigen::ArrayXf fp2){
    return fp1.matrix().dot(fp2.matrix())/fp1.size();
}

/**
 * @brief Calculate the pairwise similarity between two fingerprints using JT methods.
 * @param fp1 The first fingerprint as an Eigen::ArrayXf.
 * @param fp2 The second fingerprint as an Eigen::ArrayXf.
 * @return The pairwise similarity as a double.
 */
double pair_jt(const Eigen::ArrayXf fp1, const Eigen::ArrayXf fp2){
    double fp1_fp2 = fp1.matrix().dot(fp2.matrix());
    double fp1_fp1 = fp1.matrix().dot(fp1.matrix());
    double fp2_fp2 = fp2.matrix().dot(fp2.matrix());
    return fp1_fp2/(fp1_fp1 + fp2_fp2 - fp1_fp2);
}

/**
 * @brief Calculate the pairwise similarity between two fingerprints using SM methods.
 * @param fp1 The first fingerprint as an Eigen::ArrayXf.
 * @param fp2 The second fingerprint as an Eigen::ArrayXf.
 * @return The pairwise similarity as a double.
 */
double pair_sm(const Eigen::ArrayXf fp1, const Eigen::ArrayXf fp2){
    double fp1_fp2 = fp1.matrix().dot(fp2.matrix());
    double fp1_minus = (1-fp1).matrix().dot((1-fp2).matrix());
    return (fp1_fp2+fp1_minus)/fp1.size();
}

/**
 * @brief Calculate the average pairwise similarity for all fingerprints in the dataset.
 * @param fingerprints A 2D Eigen Array containing all the fingerprints.
 * @param n_ary A string representing the n-ary type for similarity calculation. Options are "RR", "JT", or "SM".
 * @return The average pairwise similarity as a double.
 */
double pairwise_average_real(const Eigen::ArrayXXf fingerprints, std::string n_ary){
    double avg = 0;
    int n_objects = fingerprints.rows();
    #pragma omp parallel for reduction(+:avg)
    for (int i =0; i< n_objects - 1; i++){
        for (int j = i+1; j < n_objects; j++){
            if (n_ary=="RR"){
                avg += pair_rr(fingerprints.row(i), fingerprints.row(j));
            }
            else if (n_ary=="JT"){
                avg += pair_jt(fingerprints.row(i), fingerprints.row(j));
            }
            else if (n_ary=="SM"){
                avg += pair_sm(fingerprints.row(i), fingerprints.row(j));
            }
        }
    }

    return avg/(n_objects*(n_objects-1.0)/2.0);
}

/**
 * @brief Process the matrix data by summing each column.
 * @param data A 2D Eigen Array containing the matrix data.
 * @return An Eigen::ArrayXf containing the sum of each column.
 */
Eigen::ArrayXf process_matrix_data_sum(const Eigen::ArrayXXf data){
    return data.colwise().sum();
}

/**
 * @brief Process the matrix data by summing the square of each column.
 * @param data A 2D Eigen Array containing the matrix data.
 * @return An Eigen::ArrayXf containing the sum of the square of each column.
 */
Eigen::ArrayXf process_matrix_sq_data_sum(const Eigen::ArrayXXf data){
    return data.square().colwise().sum();
}

/**
 * @brief Process the matrix data to calculate the sum of the squared differences.
 * @param data A 2D Eigen Array containing the matrix data.
 * @return A double representing the sum of the squared differences.
 */
double process_matrix_ij(const Eigen::ArrayXXf data){
    Eigen::ArrayXf col_sum = data.colwise().sum();
    Eigen::ArrayXf sq_sum = data.square().colwise().sum();
    Eigen::ArrayXf ij_arr = 0.5*(col_sum.square() - sq_sum);
    return ij_arr.sum();
}

/**
 * @brief Calculate the iSIM value for a given set of fingerprints based on the specified n-ary method.
 * @param fingerprints A 2D Eigen Array containing the fingerprints.
 * @param n_ary A string representing the n-ary type for similarity calculation. Options are "RR", "JT", or "SM".
 * @return The iSIM value as a double.
 * @throw  std::invalid_argument if an invalid n_ary value is provided.
 */
double calculate_isim_real(const Eigen::ArrayXXf fingerprints, std::string n_ary){
    if (n_ary == "RR"){
        int n_objects = fingerprints.rows();
        double ij = process_matrix_ij(fingerprints);
        int m = fingerprints.cols();
        return 2.0*ij/(m*n_objects*(n_objects-1.0));
    }
    else if (n_ary == "JT"){
        int n_objects = fingerprints.rows();
        Eigen::ArrayXf sq_sum = process_matrix_sq_data_sum(fingerprints);
        double ij = process_matrix_ij(fingerprints);
        double inners = (n_objects-1)*sq_sum.sum();
        return ij/(inners-ij);
    }
    else if (n_ary == "SM"){
        int n_objects = fingerprints.rows();
        int m = fingerprints.cols();
        double ij = process_matrix_ij(fingerprints);
        Eigen::ArrayXXf flip_data = 1- fingerprints;
        double flip_ij = process_matrix_ij(flip_data);
        return 2*(ij + flip_ij)/(m*n_objects*(n_objects-1.0));
    }
    else{
        throw std::invalid_argument("Invalid n_ary value.");
    }

}

/**
 * @brief Calculate the complementary similarity for all fingerprints in the dataset based on the RR method.
 * @param fingerprints A 2D Eigen Array containing the fingerprints.
 * @return An Eigen::ArrayXf containing the complementary similarity for all fingerprints.
 */
Eigen::ArrayXf comp_sim_rr(const Eigen::ArrayXXf fingerprints){
    int n_objects = fingerprints.rows() - 1;
    Eigen::ArrayXf sum_sq = process_matrix_sq_data_sum(fingerprints);
    Eigen::ArrayXf col_sum = process_matrix_data_sum(fingerprints);
    Eigen::ArrayXf comp_sq_sum = (col_sum.transpose().replicate(n_objects+1, 1) - fingerprints).square().rowwise().sum();
    Eigen::ArrayXf comp_sum_sq = (sum_sq.transpose().replicate(n_objects+1, 1) - fingerprints.square()).rowwise().sum();
    Eigen::ArrayXf comp_sim = (comp_sq_sum - comp_sum_sq)/(fingerprints.cols() * n_objects * (n_objects-1.0));
    return comp_sim;
}

/**
 * @brief Calculate the complementary similarity for all fingerprints in the dataset based on the SM method.
 * @param fingerprints A 2D Eigen Array containing the fingerprints.
 * @return An Eigen::ArrayXf containing the complementary similarity for all fingerprints.
 */
Eigen::ArrayXf comp_sim_sm(const Eigen::ArrayXXf fingerprints){
    int n = fingerprints.rows()-1;
    Eigen::ArrayXXf had_matrix = fingerprints.square();
    Eigen::ArrayXXf flip_matrix = 1 - fingerprints;
    Eigen::ArrayXXf had_flip = flip_matrix.square();
    Eigen::ArrayXf col_sum = process_matrix_data_sum(fingerprints);
    Eigen::ArrayXf had_sum = process_matrix_data_sum(had_matrix);
    Eigen::ArrayXf flip_sum = process_matrix_data_sum(flip_matrix);
    Eigen::ArrayXf had_flip_sum = process_matrix_data_sum(had_flip);

    Eigen::ArrayXf comp_sq_sum = (col_sum.transpose().replicate(n+1, 1) - fingerprints).square().rowwise().sum();
    Eigen::ArrayXf comp_sum_sq = (had_sum.transpose().replicate(n+1, 1) - had_matrix).rowwise().sum();
    Eigen::ArrayXf comp_flip_sq_sum = (flip_sum.transpose().replicate(n+1, 1) - flip_matrix).square().rowwise().sum();
    Eigen::ArrayXf comp_flip_sum_sq = (had_flip_sum.transpose().replicate(n+1, 1) - had_flip).rowwise().sum();

    Eigen::ArrayXf comp_sim = (comp_sq_sum + comp_flip_sq_sum - comp_sum_sq - comp_flip_sum_sq)/(fingerprints.cols() * n * (n-1.0));
    return comp_sim;
}   

/**
 * @brief Calculate the complementary similarity for all fingerprints in the dataset based on the JT method.
 * @param fingerprints A 2D Eigen Array containing the fingerprints.
 * @return An Eigen::ArrayXf containing the complementary similarity for all fingerprints.
 */
Eigen::ArrayXf comp_sim_jt(const Eigen::ArrayXXf fingerprints){
    int n_objects = fingerprints.rows() - 1;
    Eigen::ArrayXf sum_sq = process_matrix_sq_data_sum(fingerprints);
    Eigen::ArrayXf col_sum = process_matrix_data_sum(fingerprints);
    Eigen::ArrayXf comp_sq_sum = (col_sum.transpose().replicate(n_objects+1, 1) - fingerprints).square().rowwise().sum();   
    Eigen::ArrayXf comp_sum_sq = (sum_sq.transpose().replicate(n_objects+1, 1) - fingerprints.square()).rowwise().sum();
    Eigen::ArrayXf comp_sim_denom = (n_objects-1)*comp_sum_sq - 0.5*(comp_sq_sum - comp_sum_sq);
    Eigen::ArrayXf comp_sim = 0.5*(comp_sq_sum - comp_sum_sq)/comp_sim_denom;
    return comp_sim;
}

/**
 * @brief Calculate the complementary similarity for all fingerprints in the dataset based on the specified n-ary method.
 * @param fingerprints A 2D Eigen Array containing the fingerprints.
 * @param n_ary A string representing the n-ary type for similarity calculation. Options are "RR", "JT", or "SM".
 * @return An Eigen::ArrayXf containing the complementary similarity for all fingerprints.
 * @throw std::invalid_argument if an invalid n_ary value is provided.
 * @todo Increase precision to double. 
 */
// TODO : Increase precision to double 
Eigen::ArrayXf calculate_comp_sim_real(const Eigen::ArrayXXf fingerprints, std::string n_ary){
    if (n_ary == "RR"){
        return comp_sim_rr(fingerprints);
    }
    else if (n_ary == "JT"){
        return comp_sim_jt(fingerprints);
    }
    else if (n_ary == "SM"){
        return comp_sim_sm(fingerprints);
    }
    else{
        throw std::invalid_argument("Invalid n_ary value.");
    }
}

/**
 * @brief Calculate the medoid fingerprint based on the complementary similarity.
 * @param data A 2D Eigen Array containing the fingerprints.
 * @param n_ary A string representing the n-ary type for similarity calculation. Options are "RR", "JT", or "SM".
 * @return The index of the medoid fingerprint as an integer.
 */
int calculate_medoid_real(const Eigen::ArrayXXf data, const std::string n_ary){
    Eigen::ArrayXf comp_sims = calculate_comp_sim_real(data, n_ary);
    int pos;
    comp_sims.minCoeff(&pos);
    return pos;
}

/**
 * @brief Calculate the outlier fingerprint based on the complementary similarity.
 * @param data A 2D Eigen Array containing the fingerprints.
 * @param n_ary A string representing the n-ary type for similarity calculation. Options are "RR", "JT", or "SM".
 * @return The index of the outlier fingerprint as an integer.
 */
int calculate_outlier_real(const Eigen::ArrayXXf data, const std::string n_ary){
    Eigen::ArrayXf comp_sims = calculate_comp_sim_real(data, n_ary);
    int pos;
    comp_sims.maxCoeff(&pos);
    return pos;
}