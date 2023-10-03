#include <iostream>
#include <vector>
#include <omp.h>
#include <random>
#include <chrono>

// Функция для умножения матрицы на вектор
void matrixVectorMultiply(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vector, std::vector<double>& result) {
    int numRows = matrix.size();
    int numCols = matrix[0].size();

#pragma omp parallel for
    for (int i = 0; i < numRows; i++) {
        double sum = 0.0;
        for (int j = 0; j < numCols; j++) {
            sum += matrix[i][j] * vector[j];
        }
        result[i] = sum;
    }
}

int main() {
    // Задаем размеры матрицы и вектора
    int numRows = 1000;
    int numCols = 1000;
    int numThreads = 4; // Количество потоков

    // Инициализируем матрицу, вектор и результат
    std::vector<std::vector<double>> matrix(numRows, std::vector<double>(numCols));
    std::vector<double> vector(numCols);
    std::vector<double> result(numRows);

    // Заполняем матрицу и вектор случайными значениями в диапазоне от 0 до 1000
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distribution(0.0, 1000.0);

    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            matrix[i][j] = distribution(gen);
        }
    }

    for (int j = 0; j < numCols; j++) {
        vector[j] = distribution(gen);
    }

    // Устанавливаем количество потоков
    omp_set_num_threads(numThreads);

    // Замеряем время выполнения умножения 10 раз и усредняем результаты
    int numTrials = 10;
    double totalTime = 0.0;

    for (int trial = 0; trial < numTrials; trial++) {
        auto start = std::chrono::high_resolution_clock::now();
        matrixVectorMultiply(matrix, vector, result);
        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> duration = end - start;
        totalTime += duration.count();

        std::cout << "Time taken in trial " << trial + 1 << ": " << duration.count() << " seconds" << std::endl;
    }

    double averageTime = totalTime / numTrials;

    std::cout << "Average time taken over " << numTrials << " trials: " << averageTime << " seconds" << std::endl;

    return 0;
}
