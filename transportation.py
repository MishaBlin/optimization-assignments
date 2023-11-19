import numpy as np


def VogellsApproximation(matrix: list, supply: list, demand: list, n: int, m: int):
    cost_matrix = np.array(matrix)
    result = []

    penalty_row = np.array([0] * n)
    penalty_col = np.array([0] * m)

    while sum(supply) > 0 and sum(demand) > 0:
        for i in range(n):
            sorted_row = sorted(cost_matrix[i])
            penalty_row[i] = sorted_row[1] - sorted_row[0]

        for j in range(m):
            sorted_column = sorted(cost_matrix[:, j])
            penalty_col[j] = sorted_column[1] - sorted_column[0]

        max_cell_row = 0
        for i in range(n):
            if penalty_row[i] > penalty_row[max_cell_row]:
                max_cell_row = i
        max_cell_col = 0
        for j in range(m):
            if penalty_col[j] > penalty_col[max_cell_col]:
                max_cell_col = j
        min_i, min_j = max_cell_row, max_cell_col
        if penalty_row[max_cell_row] >= penalty_col[max_cell_col]:
            min_cost = float("inf")
            for i in range(m):
                if cost_matrix[max_cell_row][i] <= min_cost:
                    min_j = i
                    min_cost = cost_matrix[max_cell_row][i]
        else:
            min_cost = float("inf")
            col = cost_matrix[:, max_cell_col]
            for i in range(n):
                if col[i] <= min_cost:
                    min_i = i
                    min_cost = cost_matrix[i][max_cell_col]
        if demand[min_j] >= supply[min_i]:
            demand[min_j] -= supply[min_i]
            result.append(supply[min_i])
            supply[min_i] = 0
            cost_matrix[min_i] = [100000] * m
        else:
            supply[min_i] -= demand[min_j]
            result.append(demand[min_j])
            demand[min_j] = 0
            for i in range(n):
                cost_matrix[i][min_j] = 100000
    return result


def NorthWestCornerMethod(supply, matrix, demand):
    result = []
    n = len(supply)
    m = len(demand)
    i, j = 0, 0
    while i < n and j < m:
        if demand[j] >= supply[i]:
            result.append(supply[i])
            demand[j] -= supply[i]
            supply[i] = 0
            i += 1
            if demand[j] == 0:
                j += 1
        elif supply[i] > demand[j]:
            result.append(demand[j])
            supply[i] -= demand[j]
            demand[j] = 0
            j += 1
    return result


def RusselApproximationMethod(supply, matrix, demand):
    result = []
    n = len(supply)
    m = len(demand)
    while sum(supply) > 0:
        rows_max = [-1] * n
        columns_max = [-1] * m
        for i in range(n):
            rows_max[i] = max(matrix[i])
        for j in range(m):
            for i in range(n):
                columns_max[j] = max(columns_max[j], matrix[i][j])
        delta_i, delta_j = -1, -1
        delta = -1
        for i in range(n):
            for j in range(m):
                if (
                    columns_max[j] != -1
                    and rows_max[i] != -1
                    and matrix[i][j] != -1
                    and supply[i] != 0
                    and demand[j] != 0
                ):
                    delta_tmp = matrix[i][j] - columns_max[j] - rows_max[i]
                    if delta_tmp < delta or [delta_i, delta_j] == [-1, -1]:
                        delta = delta_tmp
                        delta_i = i
                        delta_j = j
        if demand[delta_j] >= supply[delta_i]:
            result.append(supply[delta_i])
            demand[delta_j] -= supply[delta_i]
            supply[delta_i] = 0
        elif supply[delta_i] > demand[delta_j]:
            result.append(demand[delta_j])
            supply[delta_i] -= demand[delta_j]
            demand[delta_j] = 0
        matrix[delta_i][delta_j] = -1
    return result


supply = []
matrix = []
demand = []
try:
    supply = [int(x) for x in input().split()]
    matrix = []
    for i in range(len(supply)):
        row = [int(x) for x in input().split()]
        matrix.append(row)
    demand = [int(x) for x in input().split()]
    if len(matrix) != len(supply):
        raise ValueError("method is not applicable!")
    if len(matrix[0]) != len(demand):
        raise ValueError("method is not applicable!")
except ValueError as e:
    print("Method is not applicable!")
    exit(0)

if sum(supply) != sum(demand):
    print("The problem is not balanced!")
    exit(0)

for i in range(len(supply)):
    print(" ".join(map(str, (matrix[i] + [supply[i]]))))
print(" ".join(map(str, demand)) + " -")

north_west_result = NorthWestCornerMethod(supply.copy(), matrix.copy(), demand.copy())
print(north_west_result)

russel_approximation_method = RusselApproximationMethod(
    supply.copy(), matrix.copy(), demand.copy()
)
print(russel_approximation_method)
