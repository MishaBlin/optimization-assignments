import numpy as np


def findDiff(matrix: list, n: int, m: int) -> tuple:
    row_diff = []
    col_diff = []
    for i in range(n):
        row = matrix[i][:]
        row.sort()
        if row[1] == 10**3 and row[0] != 10**3:
            row_diff.append(row[0])
            continue
        row_diff.append(row[1] - row[0])
    col = 0
    while col < m:
        column = []
        for i in range(n):
            column.append(matrix[i][col])
        column.sort()
        if column[1] == 10**3 and column[0] != 10**3:
            col_diff.append(column[0])
            col += 1
            continue
        col_diff.append(column[1] - column[0])
        col += 1
    return row_diff, col_diff


def findMaxDiff(row: list, col: list):
    max_row = max(row)
    max_col = max(col)
    max_row_index = row.index(max_row)
    max_col_index = col.index(max_col)
    return max_row, max_row_index, max_col, max_col_index


def VogelsApproximation(cost_matrix: list, supply: list, demand: list, n: int, m: int):
    result = []

    while sum(supply) > 0 and sum(demand) > 0:
        penalty_row, penalty_col = findDiff(cost_matrix, n, m)
        maxr, maxr_index, maxc, maxc_index = findMaxDiff(penalty_row, penalty_col)

        if maxr >= maxc:
            min_cell = min(cost_matrix[maxr_index])
            min_cell_i = maxr_index
            min_cell_j = cost_matrix[maxr_index].index(min_cell)
        else:
            col = [row[maxc_index] for row in cost_matrix]
            min_cell = min(col)
            min_cell_i = col.index(min_cell)
            min_cell_j = maxc_index

        if demand[min_cell_j] >= supply[min_cell_i]:
            result.append(supply[min_cell_i])
            demand[min_cell_j] -= supply[min_cell_i]
            supply[min_cell_i] = 0
            cost_matrix[min_cell_i] = [10**3 for i in range(m)]
        else:
            result.append(demand[min_cell_j])
            supply[min_cell_i] -= demand[min_cell_j]
            demand[min_cell_j] = 0
            for i in range(n):
                cost_matrix[i][min_cell_j] = 10**3
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


def RussellsApproximationMethod(supply, matrix, demand):
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


print("-" * 30)
print("Table:")
print("-" * 30)
for i in range(len(supply)):
    row = matrix[i] + [supply[i]]
    row_str = " ".join(f"{val:5}" for val in row)
    print(row_str)
demand_str = " ".join(f"{val:5}" for val in demand)
print(demand_str + f" {sum(supply):5}")
print("-" * 30)

north_west_result = NorthWestCornerMethod(supply.copy(), matrix.copy(), demand.copy())
print("North West Corner Method:" + str(north_west_result))

vogel_approximation_result = VogelsApproximation(
    matrix.copy(), supply.copy(), demand.copy(), len(supply), len(demand)
)
print("Vogel's Approximation Method:" + str(vogel_approximation_result))


russells_approximation_result = RussellsApproximationMethod(
    supply.copy(), matrix.copy(), demand.copy()
)
print("Russel's Approximation Method:" + str(russells_approximation_result))
