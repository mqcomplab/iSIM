// Function declarations
std::map<std::string, double> calculate_counters(const Eigen::ArrayXf col_sum, 
                                                const int n_objects, const int k);

std::map<std::string, double> calculate_counters( const Eigen::ArrayXXf data,
                                                 const int k);

double pairwise_average(const Eigen::ArrayXXf fingerprints, const std::string n_ary);

double calculate_isim(const Eigen::ArrayXf col_sum, 
                     const int n_objects, std::string n_ary);

Eigen::ArrayXf calculate_a(const Eigen::ArrayXf col_sum);

double calculate_isim(const Eigen::ArrayXXf data, std::string n_ary);

std::map<std::string, double> gen_sim_dict(const Eigen::ArrayXf col_sum,
const int n_objects, const int k);

std::map<std::string, double> gen_sim_dict(const Eigen::ArrayXXf data, 
const int k);

Eigen::ArrayXf calculate_comp_sim(const  Eigen::ArrayXXf data, const std::string n_ary);
int calculate_medoid(const Eigen::ArrayXXf data, const std::string n_ary);

int calculate_outlier(const Eigen::ArrayXXf data, const std::string n_ary);
