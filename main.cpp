#include <bits/stdc++.h>

using namespace std;
class Matrix {
public:
    int r;
    int c;
    vector<vector<double>> matrix;
    Matrix(int r, int c) : r(r), c(c), matrix(r, vector<double>(c)) {}

    Matrix operator+(const Matrix& other) {
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

    Matrix operator-(const Matrix& other) {
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

    Matrix operator*(const Matrix& other) {
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
                exit(0);
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

    friend ostream& operator<<(ostream&, const Matrix&);
    friend istream& operator>>(istream&, Matrix&);
};

ostream& operator<<(ostream& out, const Matrix& matrix) {
    for (int i = 0; i < matrix.r; i++) {
        for (int j = 0; j < matrix.c; j++) {
            out << matrix.matrix[i][j] << ' ';
        }
        out << '\n';
    }
    return out;
}

istream& operator>>(istream& in, Matrix& matrix) {
    for (int i = 0; i < matrix.r; i++) {
        for (int j = 0; j < matrix.c; j++) {
            in >> matrix.matrix[i][j];
        }
    }
    return in;
}

int main() {
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
    for (int i = 0; i < b.r; i++) {
        if (b.matrix[i][0] < 0) {
            cout << "The method is not applicable!" << endl;
            cout << "Negative number in the right-hand side of constraints!" << endl;
            return 0;
        }
    }
    cout << fixed << setprecision(approxAccuracy);
    vector<int> BasicVarIndices;
    vector<int> nonBasicVarIndices;
    for (int i = 0; i < z.r; i++) {
        if (z.matrix[i][0])
            nonBasicVarIndices.push_back(i);
        else
            BasicVarIndices.push_back(i);
    }
    vector<int> startNonBasic = nonBasicVarIndices;
    while (1) {
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
        for (int i = 0; i < nonBasicVarIndices.size(); i++) {
            Matrix col(BasicVarIndices.size(), 1);
            for (int j = 0; j < BasicVarIndices.size(); j++) {
                col.matrix[j][0] = A.matrix[j][nonBasicVarIndices[i]];
            }
            Matrix x(1, 1);
            x = CB0 * B0_Inverse * col;
            x.matrix[0][0] -= z.matrix[nonBasicVarIndices[i]][0];
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
            for (int i = 0; i < startNonBasic.size(); i++) {
                bool q = 0;
                for (int j = 0; j < BasicVarIndices.size(); j++) {
                    if (startNonBasic[i] == BasicVarIndices[j]) {
                        cout << 'x' << startNonBasic[i] + 1 << ":"
                             << " " << X0.matrix[j][0] << endl;
                        q = 1;
                        break;
                    }
                }
                if (!q) cout << 'x' << startNonBasic[i] + 1 << ":"
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
        for (int i = 0; i < BasicVarIndices.size(); i++) {
            if (BasicVarIndices[i] == leavingVector) {
                BasicVarIndices[i] = enteringVector;
            }
        }
        for (int i = 0; i < nonBasicVarIndices.size(); i++) {
            if (nonBasicVarIndices[i] == enteringVector) {
                nonBasicVarIndices[i] = leavingVector;
            }
        }
    }
}
