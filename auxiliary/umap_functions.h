// Functions to perform operations on unordered_map class

#ifndef UMAP_FUNCTIONS
#define UMAP_FUNCTIONS

#include <unordered_map>
#include <string>

#include "../newmat10/newmat.h"

using namespace std;

//Functions

unordered_map<string, double> umap_sum_each_key(unordered_map<string, RowVector> umap);
double umap_sum_over_keys(unordered_map<string, double> umap);
double umap_sum_over_keys_1_lag(unordered_map<string, RowVector> umap, int p);

#endif