#include "umap_functions.h"

/// @brief Given an unordered map that has RowVectors as values, return with as values the sum of the corresponding RowVector
/// @param umap 
/// @return 
unordered_map<string, double> umap_sum_each_key(unordered_map<string, RowVector> umap){
    unordered_map<string, double> sumap; 
    for (auto el:umap){
        sumap[el.first] = umap[el.first].Sum();
    }
    return sumap;
}

/// @brief Given an unordered map that has double as values, return the sum of the values
/// @param umap
/// @return 
double umap_sum_over_keys(unordered_map<string, double> umap){
    double sum=0;
    for (auto x:umap){
        sum += x.second;
    }
    return sum;
}

/// @brief Given an unordered map that has RowVectors as values, return the sum of the values in position p
/// @param umap 
/// @return 
double umap_sum_over_keys_1_lag(unordered_map<string, RowVector> umap, int p){
    double sum=0;
    for (auto x:umap){
        sum += x.second(p);
    }
    return sum;
}

