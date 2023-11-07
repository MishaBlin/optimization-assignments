#include <bits/stdc++.h>

using namespace std;

class Matrix {
public:
    int r;
    int c;
    vector<vector<double>> matrix;

    Matrix(int r, int c) : r(r), c(c), matrix(r, vector<double>(c)) {}

    Matrix operator+(const Matrix &other) {
        if (r != other.r or c != other.c) {
            throw invalid_argument("Matrices with different dimensions");
        }
        Matrix res(r, c);
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                res.matrix[i][j] = this->matrix[i][j] + other.matrix[i][j];
            }
        }
        return res;
    }

    Matrix operator-(const Matrix &other) {
        if (r != other.r or c != other.c) {
            throw invalid_argument("Matrices with different dimensions");
        }
        Matrix res(r, c);
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                res.matrix[i][j] = this->matrix[i][j] - other.matrix[i][j];
            }
        }
        return res;
    }

    Matrix operator*(const Matrix &other) {
        if (c != other.r) {
            throw invalid_argument("Cannot perform multiplication for such dimensions");
        }
        Matrix res(r, other.c);
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < other.c; j++) {
                for (int k = 0; k < c; k++) {
                    res.matrix[i][j] += this->matrix[i][k] * other.matrix[k][j];
                }
            }
        }
        return res;
    }

    Matrix transpose() {
        Matrix res(c, r);
        for (int i = 0; i < c; i++) {
            for (int j = 0; j < r; j++) {
                res.matrix[i][j] = this->matrix[j][i];
            }
        }
        return res;
    }

    Matrix inverse() {
        if (r != c) {
            throw invalid_argument("Inverse matrix does not exist");
        }

        Matrix augmMatrix(r, 2 * c);
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                augmMatrix.matrix[i][j] = this->matrix[i][j];
            }
        }

        for (int j = 0; j < c; j++) {
            augmMatrix.matrix[j][r + j] = 1;
        }

        for (int i = 0; i < r; i++) {
            int mxElemRow = i;
            for (int j = i + 1; j < r; j++) {
                if (abs(augmMatrix.matrix[j][i]) > abs(augmMatrix.matrix[mxElemRow][i])) {
                    mxElemRow = j;
                }
            }
            if (mxElemRow != i) {
                swap(augmMatrix.matrix[mxElemRow], augmMatrix.matrix[i]);
            }
            if (!augmMatrix.matrix[i][i]) {
                cout << "The method is not applicable!" << endl;
                return {0, 0};
            }
            for (int j = 0; j < r; j++) {
                if (i == j) {
                    continue;
                }
                double temp = augmMatrix.matrix[j][i] / augmMatrix.matrix[i][i];
                for (int k = 0; k < 2 * c; k++) {
                    augmMatrix.matrix[j][k] -= augmMatrix.matrix[i][k] * temp;
                }
            }
        }

        Matrix res(r, c);
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                res.matrix[i][j] = augmMatrix.matrix[i][j + c] / augmMatrix.matrix[i][i];
            }
        }
        return res;
    }

    friend ostream &operator<<(ostream &, const Matrix &);

    friend istream &operator>>(istream &, Matrix &);
};

ostream &operator<<(ostream &out, const Matrix &matrix) {
    for (int i = 0; i < matrix.r; i++) {
        for (int j = 0; j < matrix.c; j++) {
            out << matrix.matrix[i][j] << ' ';
        }
        out << '\n';
    }
    return out;
}

istream &operator>>(istream &in, Matrix &matrix) {
    for (int i = 0; i < matrix.r; i++) {
        for (int j = 0; j < matrix.c; j++) {
            in >> matrix.matrix[i][j];
        }
    }
    return in;
}

bool
solution_existence_by_simplex(Matrix A, Matrix X0, vector<int> startNonBasic, vector<int> BasicVarIndices, Matrix b) {
    Matrix x(b.r, b.c);
    for (int i = 0; i < startNonBasic.size(); i++) {
        bool q = false;
        for (int j = 0; j < BasicVarIndices.size(); j++) {
            if (startNonBasic[i] == BasicVarIndices[j]) {
                x.matrix[i][0] = X0.matrix[j][0];
                q = true;
                break;
            }
        }
        if (!q) x.matrix[i][0] = 0;
    }

    for (int i = 0; i < A.r; i++) {
        double sum = 0;
        for (int j = 0; j < startNonBasic.size(); j++) {
            sum += x.matrix[j][0] * A.matrix[i][startNonBasic[j]];
        }
        if (round(sum) > b.matrix[i][0]) return false;
    }

    return true;
}

void simplex(bool problemKind, Matrix A, Matrix b, Matrix z) {

    // check if we have '>=' in our constraints

    bool flag = false;
    for (int i = 0; i < z.r; i++) {
        if (z.matrix[i][0] == 0) {
            for (int j = 0; j < A.r; j++) {
                if (A.matrix[j][i] < 0) flag = true;
            }
        }
    }
    if (flag) {
        cout << "The method is not applicable!\n";
        return;
    }

    for (int i = 0; i < b.r; i++) {
        if (b.matrix[i][0] < 0) {
            cout << "Negative number in the right-hand side of constraints!" << endl;
            exit(0);
        }
    }
    vector<int> BasicVarIndices;
    vector<int> nonBasicVarIndices;
    for (int i = 0; i < z.r; i++) {
        if (z.matrix[i][0])
            nonBasicVarIndices.push_back(i);
        else
            BasicVarIndices.push_back(i);
    }
    vector<int> startNonBasic = nonBasicVarIndices;
    while (true) {
        Matrix CB0(1, BasicVarIndices.size());
        for (int i = 0; i < BasicVarIndices.size(); i++) {
            CB0.matrix[0][i] = z.matrix[BasicVarIndices[i]][0];
        }
        Matrix B0(BasicVarIndices.size(), BasicVarIndices.size());
        for (int i = 0; i < BasicVarIndices.size(); i++) {
            for (int j = 0; j < BasicVarIndices.size(); j++) {
                B0.matrix[i][j] = A.matrix[i][BasicVarIndices[j]];
            }
        }
        Matrix B0_Inverse = B0.inverse();
        Matrix X0 = B0_Inverse * b;
        vector<double> z_c;
        for (int nonBasicVarIndice: nonBasicVarIndices) {
            Matrix col(BasicVarIndices.size(), 1);
            for (int j = 0; j < BasicVarIndices.size(); j++) {
                col.matrix[j][0] = A.matrix[j][nonBasicVarIndice];
            }
            Matrix x(1, 1);
            x = CB0 * B0_Inverse * col;
            x.matrix[0][0] -= z.matrix[nonBasicVarIndice][0];
            z_c.push_back(x.matrix[0][0]);
        }
        double mn = 1e9;
        int enteringVector = -1;
        for (int i = 0; i < z_c.size(); i++) {
            if (z_c[i] < mn) {
                enteringVector = nonBasicVarIndices[i];
                mn = z_c[i];
            }
        }

        if (mn >= 0) {

            bool method_has_solution = solution_existence_by_simplex(A, X0, startNonBasic, BasicVarIndices, b);
            if (!method_has_solution) {
                cout << "The problem does not have solution!\n";
                return;
            }

            for (int i: startNonBasic) {
                bool q = false;
                for (int j = 0; j < BasicVarIndices.size(); j++) {
                    if (i == BasicVarIndices[j]) {
                        cout << 'x' << i + 1 << ":"
                             << " " << X0.matrix[j][0] << endl;
                        q = true;
                        break;
                    }
                }
                if (!q)
                    cout << 'x' << i + 1 << ":"
                         << " " << 0 << endl;
            }
            Matrix Z = CB0 * X0;
            int coeff = (problemKind ? 1 : -1);
            cout << "z: " << coeff * Z.matrix[0][0] << endl;
            break;
        }
        Matrix EnteringVec(BasicVarIndices.size(), 1);
        for (int i = 0; i < BasicVarIndices.size(); i++) {
            EnteringVec.matrix[i][0] = A.matrix[i][enteringVector];
        }
        Matrix FC = B0_Inverse * EnteringVec;
        vector<double> x1;
        for (int i = 0; i < BasicVarIndices.size(); i++) {
            if (FC.matrix[i][0] != 0 && X0.matrix[i][0] / FC.matrix[i][0] >= 0) {
                x1.push_back(X0.matrix[i][0] / FC.matrix[i][0]);
            } else
                x1.push_back(1e9);
        }
        mn = 1e9;
        int leavingVector = -1;
        for (int i = 0; i < x1.size(); i++) {
            if (x1[i] < mn) {
                mn = x1[i];
                leavingVector = BasicVarIndices[i];
            }
        }
        for (int &BasicVarIndice: BasicVarIndices) {
            if (BasicVarIndice == leavingVector) {
                BasicVarIndice = enteringVector;
            }
        }
        for (int &nonBasicVarIndice: nonBasicVarIndices) {
            if (nonBasicVarIndice == enteringVector) {
                nonBasicVarIndice = leavingVector;
            }
        }
    }
}

bool solution_existence_by_interior_point(Matrix A, Matrix x, Matrix b) {
    Matrix A_multiplied_x = A * x;
    for (int i = 0; i < b.r; i++) {
        if (round(A_multiplied_x.matrix[i][0]) > b.matrix[i][0]) {
            return false;
        }
    }
    return true;
}

void interior_point(bool kind, int coefficientsNumber, Matrix A, Matrix z, double alpha, Matrix b, Matrix x) {
    Matrix D(coefficientsNumber, coefficientsNumber);

    while (true) {
        Matrix x_init = x;

        for (int i = 0; i < coefficientsNumber; i++) {
            D.matrix[i][i] = x_init.matrix[i][0];
        }

        Matrix A_hat = A * D;
        Matrix c_hat = D * z;

        Matrix T = A_hat.transpose();
        Matrix N = (A_hat * T).inverse();

        // check if matrix passed inverse function validation check (for identifying if method applicable or not)
        if (N.r == 0 && N.c == 0) return;

        Matrix I = Matrix(coefficientsNumber, coefficientsNumber);
        for (int i = 0; i < coefficientsNumber; i++) {
            I.matrix[i][i] = 1;
        }

        Matrix P = I - T * N * A_hat;

        Matrix c_p = P * c_hat;

        auto v = DBL_MIN;

        for (int i = 0; i < coefficientsNumber; i++) {
            double curr = c_p.matrix[i][0];

            if (curr < 0 && abs(curr) > v) {
                v = abs(curr);
            }
        }

        Matrix X_hat(coefficientsNumber, 1);

        for (int i = 0; i < coefficientsNumber; i++) {
            X_hat.matrix[i][0] = 1 + (alpha / v) * c_p.matrix[i][0];
        }

        x = D * X_hat;

        double sum = 0;
        for (int i = 0; i < coefficientsNumber; i++) {
            double diff = x.matrix[i][0] - x_init.matrix[i][0];
            sum += diff * diff;
        }

        if (sqrt(sum) < 0.00001) {
            break;
        }
    }

    // check if interior point has solution (look if our obtained x in feasible region)
    bool method_has_solution = solution_existence_by_interior_point(A, x, b);

    if (method_has_solution) {
        for (int i = 0; i < coefficientsNumber; i++) {
            if (z.matrix[i][0] != 0) {
                cout << "x" << i + 1 << ": " << x.matrix[i][0] << "\n";
            }
        }

        double optima = (z.transpose() * x).matrix[0][0];
        if (!kind) {
            optima *= -1;
        }

        cout << "z: " << optima;
    } else {
        cout << "The problem does not have solution\n";
    }
}


int main() {
// Input
    cout << "Enter 1 if this problem is a maximizing problem and 0 if it's a minimizing problem" << endl;
    bool problemKind;
    cin >> problemKind;
    cout << "Enter the number coefficients of objective function (including the slack variables)" << endl;
    int coefficientsNumber;
    cin >> coefficientsNumber;
    cout << "Enter the objective function coefficients" << endl;
    Matrix z(coefficientsNumber, 1);
    cin >> z;
    if (!problemKind) {
        for (int i = 0; i < z.r; i++) {
            z.matrix[i][0] *= -1;
        }
    }
    int numberOfBasicVar = 0;
    for (int i = 0; i < z.r; i++) {
        if (!z.matrix[i][0]) numberOfBasicVar++;
    }
    Matrix A(numberOfBasicVar, coefficientsNumber);
    cout << "The matrix of coefficients of constraint function" << endl;
    cin >> A;
    Matrix b(numberOfBasicVar, 1);
    cout << "The vector of right-hand side numbers" << endl;
    cin >> b;
    cout << "The approximation accuracy" << endl;
    int approxAccuracy;
    cin >> approxAccuracy;
    cout << fixed << setprecision(approxAccuracy);
    Matrix x(coefficientsNumber, 1);
    cout << "The initial trial solution" << endl;
    cin >> x;

// Interior point
    cout << "Interior point solution with alpha = 0.5:\n";
    interior_point(problemKind, coefficientsNumber, A, z, 0.5, b, x);

    cout << "\n\nInterior point solution with alpha = 0.9:\n";
    interior_point(problemKind, coefficientsNumber, A, z, 0.9, b, x);

// Simplex
    cout << "\n\nSolution with simplex method:\n";
    simplex(problemKind, A, b, z);
}
