#include <iostream>
#include <vector>
#include <omp.h>

const int MATRIX_SIZE = 1000;
const int BLOCK_SIZE = 100;

// Функция для заполнения матрицы случайными числами
void fillMatrix(std::vector<std::vector<int>>& matrix) {
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            matrix[i][j] = rand() % 1001; // случайное число от 0 до 1000
        }
    }
}

// Функция для умножения двух матриц
void multiplyMatrixBlock(const std::vector<std::vector<int>>& matrixA,
    const std::vector<std::vector<int>>& matrixB,
    std::vector<std::vector<int>>& result,
    int rowStart, int rowEnd, int colStart, int colEnd) {
    for (int i = rowStart; i < rowEnd; ++i) {
        for (int j = colStart; j < colEnd; ++j) {
            result[i][j] = 0;
            for (int k = 0; k < MATRIX_SIZE; ++k) {
                result[i][j] += matrixA[i][k] * matrixB[k][j];
            }
        }
    }
}

int main() {
    // Инициализация матриц
    std::vector<std::vector<int>> matrixA(MATRIX_SIZE, std::vector<int>(MATRIX_SIZE));
    std::vector<std::vector<int>> matrixB(MATRIX_SIZE, std::vector<int>(MATRIX_SIZE));
    std::vector<std::vector<int>> result(MATRIX_SIZE, std::vector<int>(MATRIX_SIZE));

    fillMatrix(matrixA);
    fillMatrix(matrixB);

    double total_time = 0.0;
    const int num_runs = 10;

    for (int run = 0; run < num_runs; ++run) {
        double start_time = omp_get_wtime();

#pragma omp parallel for
        for (int i = 0; i < MATRIX_SIZE; i += BLOCK_SIZE) {
            for (int j = 0; j < MATRIX_SIZE; j += BLOCK_SIZE) {
                multiplyMatrixBlock(matrixA, matrixB, result, i, i + BLOCK_SIZE, j, j + BLOCK_SIZE);
            }
        }

        double end_time = omp_get_wtime();
        double run_time = end_time - start_time;
        total_time += run_time;
        std::cout << "Run " << run + 1 << " time: " << run_time << " seconds" << std::endl;
    }

    std::cout << "Average time over " << num_runs << " runs: " << total_time / num_runs << " seconds" << std::endl;

    return 0;
}

