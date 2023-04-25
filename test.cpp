#include <iostream> 
#include <map>
#include <vector>
#include <algorithm>

int main()

{

  typedef std::map< std::vector<int> , int > VecToInts; 

  std::vector<int> v1{10, 20, 30}; 
  std::vector<int> v2{ 5, 15, 25}; 

  VecToInts map; 
  map[v1] =1; 
  map[v2] =2; 

  int element_to_find = 10; 
  int cluster_ID = 0; 

  VecToInts::const_iterator iter = map.find(v1); 
  if( map.end() != iter) {
    std::cout << iter->second << std::endl;
    std::vector<int> myvector = iter->first; 

    if (std::find(myvector.begin(), myvector.end(), element_to_find) != myvector.end()){
      std::cout << "found: " <<  element_to_find <<std::endl;
    }
    else std::cout<< 0 << std::endl;
  }

  for (auto const& x : map)
    {
      std::vector<int> myvector = x.first; 
      if (std::find(myvector.begin(), myvector.end(), element_to_find) != myvector.end()){
	cluster_ID = x.second; 	
      }
    }
  std::cout << cluster_ID << std::endl;
  return 0; 
}
