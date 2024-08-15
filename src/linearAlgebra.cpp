#include "../include/linearAlgebra.hpp"

//Display 1D vector
void displayVector(const vector<double> &vec)
{
    for (const auto &v : vec)
        cout << v << " ";
    cout << endl;
}

//Display 2D matrix
void displayMatrix(const vector<vector<double>> &mat)
{
    for (size_t i{0}; i < mat.size(); ++i)
    {
        for (size_t j{0}; j < mat[i].size(); ++j)
        {
            cout << mat[i][j] << " ";
        }
        cout << endl;
    }
}

//Generating 1D field with the uniform spacing
vector<double> linspace(double x0, double xn, size_t size)
{
    vector<double> space(size);
    double dx{0};
    
    dx = (xn - x0)/static_cast<double>(size-1);

    for (size_t i{0}; i < size; ++i)
    {
        space[i] = x0 + dx*static_cast<double>(i);
    }

    return space;
}

// Overloaded matrix multiplication functions 
vector<double> matMul(const vector<vector<double>>& mat, const vector<double> vec)
{
    int numRow = mat.size();
    int numCol = mat[0].size();
    int dimVec = vec.size();
    vector<double> temp(dimVec);

    if (numCol != dimVec) {
        cout << "Dimension of matrix and vector is not matched";
        return temp;
    }
    
    for (int i{0}; i < numRow; ++i) {
        for (int j{0}; j < numCol; ++j) {
            temp[i] = temp[i] + mat[i][j]*vec[j];
        }
    }
    
    return temp;
}

vector<double> slice(const vector<double> &vec, int m, int n) {
    int start = (m < 0) ? 0 : m;
    int end   = (n > int(vec.size())) ? vec.size() : n;

    if (start > end) {
        return {};
    }

    return vector<double>(vec.begin() + start, vec.begin() + end);

}

void replace(vector<double> &vec1, vector<double> &&vec2, size_t m)
/*vec1 : The oringinal vector, vec2 : The replacement vector, m : starting point, n : end point*/
{
    if (vec1.size() < m-1 + vec2.size()) {
        cout << "The target size is small for replacement" << endl;
    } else
    {
        copy(vec2.begin(), vec2.end(), vec1.begin() + m);
    }
}

void replace(vector<double> &vec1, vector<double> &vec2, size_t m)
/*vec1 : The oringinal vector, vec2 : The replacement vector, m : starting point, n : end point*/
{
    if (vec1.size() < m-1 + vec2.size()) {
        cout << "The target size is small for replacement" << endl;
    } else
    {
        copy(vec2.begin(), vec2.end(), vec1.begin() + m);
    }
}

double dotProd(const vector<double> &vec1, const vector<double> &vec2)
{
    size_t dim{vec1.size()};
    double temp{0};

    if (vec1.size() != vec2.size()) {
        cout << "Dimension of two vectors is not matched";
        return temp;
    }

    for (size_t i{0}; i < dim; ++i)
        temp = temp + vec1[i]*vec2[i];
    
    return temp;
}

double norm(const vector<double> &vec, string &&type)
{
    double temp{0};
    
    size_t dim{static_cast<size_t>(vec.size())};

    if (type == "norm2") {
        for (size_t i{0}; i < dim; ++i) {
            temp += std::pow((vec[i]), 2);
        }
        temp = std::pow(temp, 0.5);    
    }

    else if (type == "max") {
        temp = vec[0];
        for (size_t i{0}; i < dim; ++i) {
            if (temp <= vec[i]) {
                temp = vec[i];
            } else {
                continue;
            }
        }
    }

    return temp;

}

// Gauss-Jordan Elimination
void gaussJordan(vector<vector<double>>& mat) {
    int numRow = mat.size();
    int numCol = mat[0].size();

    for (int i = 0; i < numRow; ++i) {
        // Partial pivoting: 가장 큰 값을 가진 행을 현재 행으로 이동시킵니다.
        int maxRowIndex = i;
        double maxVal = mat[i][i];
        for (int k = i + 1; k < numCol; ++k) {
            if (abs(mat[k][i]) > abs(maxVal)) {
                maxVal = mat[k][i];
                maxRowIndex = k;
            }
        }
        swap(mat[i], mat[maxRowIndex]); // 행 교환

        // Normalize the row: 현재 행을 정규화합니다. 대각 요소를 1로 만듭니다.
        double pivot = mat[i][i];
        for (int j = i; j < numCol; ++j) {
            mat[i][j] /= pivot;
        }

        // Elimination step: 다른 행들을 수정하여 0으로 만듭니다.
        for (int k = 0; k < numRow; ++k) {
            if (k == i) continue;
            double factor = mat[k][i];
            for (int j = i; j < numCol; ++j) {
                mat[k][j] -= factor * mat[i][j];
            }
        }
    }
}

// Gauss-Jordan Elimination for Inverse
void gaussJordanInv(vector<vector<double>>& mat) {
    int numRow = mat.size();
    for (int i = 0; i < numRow; i++) {
        mat[i].resize(2 * numRow, 0);
        mat[i][numRow + i] = 1; // 단위 행렬을 추가
    }

    for (int i = 0; i < numRow; i++) {
        // Pivot 선택
        int maxRow = i;
        for (int k = i + 1; k < numRow; k++) {
            if (abs(mat[k][i]) > abs(mat[maxRow][i])) {
                maxRow = k;
            }
        }
        if (mat[maxRow][i] == 0.0) {
            cout << "The matrix is singular and cannot be inverted." << endl; // 특이행렬 판정
        }
        swap(mat[i], mat[maxRow]);

        double pivotValue = mat[i][i];
        for (int j = 0; j < 2 * numRow; j++) {
            mat[i][j] /= pivotValue;
        }

        for (int k = 0; k < numRow; k++) {
            if (k != i) {
                double factor = mat[k][i];
                for (int j = 0; j < 2 * numRow; j++) {
                    mat[k][j] -= factor * mat[i][j];
                }
            }
        }
    }

    for (int i = 0; i < numRow; i++) {
        for (int j = 0; j < numRow; j++) {
            mat[i][j] = mat[i][j + numRow];
        }
        mat[i].resize(numRow);
    }
}

vector<double> TDMA(const vector<double>& a, const vector<double>& b, const vector<double>& c, const vector<double>& d) {
    int n = d.size();
    vector<double> c_star(n, 0.0), d_star(n, 0.0), x(n, 0.0);

    //forward elimination
    c_star[0] = c[0] / b[0];
    d_star[0] = d[0] / b[0];
    for (int i = 1; i < n; i++) {
        double m = 1.0 / (b[i] - a[i] * c_star[i - 1]);
        c_star[i] = c[i] * m;
        d_star[i] = (d[i] - a[i] * d_star[i - 1]) * m;
    }

    //backward substitution
    x[n - 1] = d_star[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = d_star[i] - c_star[i] * x[i + 1];
    }

    return x;
}