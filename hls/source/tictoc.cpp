#include "tictoc.h"

tictoc::tictoc(int len){
	tocs.reserve(len);
}
void tictoc::reset(){
	tocs.clear();
}

void tictoc::exportTocs(std::ofstream &csvOut){
	int size = tocs.size();
	for (int i=0; i<(size-1); i++){
        csvOut << tocs[i] << ", ";
    }
    csvOut << tocs.back() << "\n";
}

void tictoc::rmWarmup(int n){
	tocs.erase(tocs.begin(), tocs.begin() + n);
}

double tictoc::mean(){
	return std::accumulate(tocs.begin(), tocs.end(), 0.0)/tocs.size();
}
double tictoc::stdev(){
	std::vector<double> diff(tocs.size());
	std::transform(tocs.begin(), tocs.end(), diff.begin(),
	               std::bind2nd(std::minus<double>(), mean()));
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);

	return  std::sqrt(sq_sum / (tocs.size()-1));
}

double tictoc::max(){
	return  *max_element(tocs.begin(), tocs.end());
}

double tictoc::min(){
	return  *min_element(tocs.begin(), tocs.end());
}
