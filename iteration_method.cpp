#include <iostream>
#include <vector>
#include <cmath>
#include "gauss.cpp"


void transforming(std::vector<std::vector<double>>&, std::vector<double>&);
double infinity_norm_of_vector(std::vector<double>);
double infinity_norm_of_matrix(std::vector<std::vector<double>>);
std::vector<double> simple_iteration_method(std::vector<std::vector<double>> , std::vector<double> , int);
std::vector<double> operator-(std::vector<double>, std::vector<double>);


void transforming(std::vector<std::vector<double>> &A_matrix, std::vector<double> &b_vector){
    int n = b_vector.size();
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            if(i == j)
                A_matrix[i][j] = 0;
            else
                A_matrix[i][j] /= A_matrix[i][i];
        }
        b_vector[i] /= A_matrix[i][i];
    }
}

double infinity_norm_of_vector(std::vector<double> vector) {
    double max_el = fabs(vector[0]);
    for(int i = 1; i < vector.size(); ++i){
       if (vector[i] > max_el){
           max_el = vector[i];
       }
    }
    return max_el;
}

double infinity_norm_of_matrix(std::vector<std::vector<double>> matrix){
    double max_sum = 0;
    for(int i = 0; i < matrix[0].size(); ++i){
        double temp_sum = 0;
        for(int j = 0; j < matrix[0].size(); ++j){
            temp_sum += fabs(matrix[i][j]);
        }
        if (temp_sum > max_sum)
            max_sum = temp_sum;
    }
    return max_sum;
}

std::vector<double> simple_iteration_method(std::vector<std::vector<double>> H_matrix, std::vector<double> g_vector, int iteration){
    int n = g_vector.size();
    std::vector<double> x_vector_current(n, 0);
    for(int k = 0; k < iteration; ++k){
        for(int i = 0; i < n; ++i) {
            double sum = 0;
            for (int j = 0; j < n; ++j) {
                sum += H_matrix[i][j] * x_vector_current[j];
            }
            x_vector_current[i] = sum + g_vector[i];
        }
    }
    return x_vector_current;
}

std::vector<double> operator-(std::vector<double> v1, std::vector<double> v2){
    int size1 = v1.size(), size2 = v2.size();
    for(int i = size1; i < size2; ++i){
        v1.push_back(0);
    }
    for(int i = size2; i < size1; ++i){
        v2.push_back(0);
    }
    for(int i = 0; i < v1.size(); ++i){
        v1[i] -= v2[i];
    }
    return v1;
}


int main() {
    int n;
    std::cout << "Enter number of rows and cols: " << std::endl;
    std::cin >> n;
    std::vector<std::vector<double>> arr(n, std::vector<double> (n));
    std::cout << "Enter A matrix values: " << std::endl;
    std::string string;
    std::getline(std::cin, string);     //first getline for reading /n symbol, which has been left by std::cout
    std::vector<double> values_temp;
    for(int i = 0; i < n; ++i){
        std::getline(std::cin, string);
        values_temp = split_string_into_array(string);
        for(int j = 0; j < n; ++j){
            arr[i][j] = values_temp[j];
        }
    }
    std::vector<double> b_vector;
    std::cout << "Enter b-vector:" << std::endl;
    getline(std::cin, string);
    b_vector = split_string_into_array(string);
    transforming(arr, b_vector);            //"A" matrix is "H" now, meanwhile "b_vector" is "g" vector now
    std::cout << "Infinity norm for H matrix is: " << infinity_norm_of_matrix(arr) << std::endl;
    std::vector <double> accurate_answer = Gauss(b_vector, arr, n, n);
    std::cout << "Accurate answer is: ";        //Gauss answer
    for(auto el : accurate_answer)
        std::cout << el << " ";
    std::cout << std::endl;
    std::vector<double> approximate_answer = simple_iteration_method(arr, b_vector, 7);
    std::cout << "A priori error estimate of ||x^(7) - x*||_{infinity} "
                 "= ||H||^k||x^(0)|| + ||H||^k/(1 - ||H||)*||g|| = " << pow(infinity_norm_of_matrix(arr), 7)/
                 (1 - infinity_norm_of_matrix(arr)) * infinity_norm_of_vector(b_vector) << std::endl;
    std::cout << "Approximate value of x^(7) is ";
    for(auto el : approximate_answer)
        std::cout << el << " ";
    std::cout << std::endl;
    std::cout << "Actual error is: " << infinity_norm_of_vector(approximate_answer - accurate_answer) << std::endl;
    std::cout << "A posterior error estimate of ||x^(7) - x*||_{infinity}: "
                << infinity_norm_of_matrix(arr)/(1 - infinity_norm_of_matrix(arr))
                *infinity_norm_of_vector(approximate_answer - simple_iteration_method(arr, b_vector, 6));
}