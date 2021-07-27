#include <iostream>
#include <list>
#include <cstdlib>

// set up the overlap buffer, make the filter image. 
    int findGoodSize(std::list<int> factors, int best_so_far, int test_size, int min, int max)
    {
      if(factors.empty()) return best_so_far; 
      int m = factors.front();
      std::list<int> rfactors = factors;
      rfactors.pop_front();
      std::cerr << "fgs(..," << best_so_far << ", " << test_size << ", " << min << ", " << max << ")\n";
      for(int f = test_size; f < max; f *= m) {
	if((f >= min) && (f < best_so_far)) {
	  best_so_far = f;
	}
	best_so_far = findGoodSize(rfactors, best_so_far, f, min, max);
      }
      return best_so_far; 
    }

int main(int argc, char * argv[]) {
  int fs, vs, rm;

  fs = atoi(argv[1]);
  vs = atoi(argv[2]);
  rm = atoi(argv[3]);

  int target_min = vs + fs - 1;
  int target_max = target_min * 3 / 2; 
  
  std::cerr << "fs = " << fs << " vs = " << vs << " rm = " << rm << " " ;

  std::list<int> factors;
  factors.push_back(2); factors.push_back(3); factors.push_back(5); factors.push_back(7); 
  int good_size = findGoodSize(factors, (1 << 24), rm, target_min, target_max);
  
  std::cerr << " goodsize = " << good_size << "\n";
}
