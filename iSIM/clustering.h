#include "isim.h"
#include "comp.h"

struct processedFps pre_process(Eigen::ArrayXXf data);
std::list<int> max_indices(struct processedFps fps, std::string n_ary);