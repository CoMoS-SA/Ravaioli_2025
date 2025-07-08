// Functions to read values from inputs and generate unordered_maps for households classes 

#ifndef READ_MH_INPUTS
#define READ_MH_INPUTS

#include <vector>
#include <unordered_map>
#include <string>

#include "../newmat10/newmat.h"
#include "../rapidjson/document.h"

using namespace std;

//Functions
RowVector read_input_vector(RowVector vector, const rapidjson::Value& input, int dim);
unordered_map<string, double> create_umap(vector<string> keys, RowVector values);
unordered_map<string, RowVector> create_umap_1_lag(vector<string> keys, RowVector values);
unordered_map<string, double> read_input_vector_into_umap(vector<string> keys, const rapidjson::Value& input);

#endif