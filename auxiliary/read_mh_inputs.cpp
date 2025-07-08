#include "read_mh_inputs.h"


/// @brief Read inside a RowVector the values from an input file
/// @param vector 
/// @param input 
/// @param dim 
/// @return 
RowVector read_input_vector(RowVector vector, const rapidjson::Value& input, int dim){
    vector.ReSize(dim);
    for(int i=1; i<=dim; i++){
        vector(i)=input[(i-1)].GetDouble();
    }
    return vector;
}

/// @brief Create an unordered_map from keys (vector of strings) and values (RowVector) 
/// @param keys
/// @param values 
/// @return 
unordered_map<string, double> create_umap(vector<string> keys, RowVector values){
    int dim = keys.size();
    unordered_map<string, double> umap;
    for (int i=1; i<=dim; i++) {
        umap[keys[i-1]] = values(i);
    }
    return umap;
}

/// @brief Create an unordered_map [keys:RowVector] from keys (vector of strings) and values (RowVector) including 1 lag 
/// @param keys
/// @param values 
/// @return 
unordered_map<string, RowVector> create_umap_1_lag(vector<string> keys, RowVector values){
    int dim = keys.size();
    unordered_map<string, RowVector> umap;
    for (int i=1; i<=dim; i++) {
        RowVector key_values(2);
        key_values = values(i);
        umap[keys[i-1]] = key_values;
    }
    return umap;
}

/// @brief Create an unordered_map [keys:RowVector] from keys (vector of strings) and values as input
/// @param keys 
/// @param input values read from input Document
/// @return 
unordered_map<string, double> read_input_vector_into_umap(vector<string> keys, const rapidjson::Value& input){
    int dim = keys.size();
    RowVector vector;
    vector=read_input_vector(vector, input, dim);
    unordered_map<string, double> umap;
    umap=create_umap(keys, vector);
    return umap;
}