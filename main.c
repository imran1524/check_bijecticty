#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int i = 256;
int j = 8;
int** temp;

// Function to generate a random binary value (0 or 1)
int generateBinaryValue() {
    return rand() % 2;
}
int**v_i_j;
int** allocate_2D_matrix(int rows, int cols);
int** generateSboxMatrix(int rows, int cols);
void free_2D_array(int** array_2D, int array_2D_length);
int* calculate_x_v(int* A, int n);
int calculate_F_D_j(int *x_k_j, int n, int cols);
int** generate_WHM(int order);
//Function definition of generating Walsh Hadamard Matrix of order k
int main() {
    int **ptr_sbox_matrix = allocate_2D_matrix(i, j);
    ptr_sbox_matrix = generateSboxMatrix(i, j);

    int offset = 0;
    //CALCULATING v_i_j
    v_i_j = allocate_2D_matrix(i, j);
    for (int rows = 0; rows < i; rows++) {
        offset = 1;
        int sum = 0;
        for (int cols = j - 1; cols >= 0; cols--) {
            sum = ptr_sbox_matrix[rows][cols] * 1 << offset - 1;
            v_i_j[rows][offset - 1] = sum;
            offset++;
        }
    }
    //PRINTING THE v_i_j(j)
    for (int cols = 1; cols < j; cols++) {
        for (int rows = 0; rows < i; rows++) {
            v_i_j[rows][cols] = v_i_j[rows][cols] + v_i_j[rows][cols - 1];
        }
    }

    int* temp = (int*) malloc(j * sizeof(int));
    int* x_v = (int*) malloc(i * sizeof(int));
    int* F_D_j = (int*) malloc(j * sizeof(int));

    for (int cols = 0; cols < j; cols++) {
        printf("j = %d\n", cols);
        for (int rows = 0; rows < i; rows++) {
            temp[rows] = v_i_j[rows][cols];
            printf("v_j[%d] = %d\n",rows, temp[rows]);
        }

        printf("2^%d-1-%d = %d\n", 8, cols, 1 << (8-1-cols));

        int sum = 0;
        //CALCULATING #x_v(j)
        x_v = calculate_x_v(temp, i);
        for (int index = 0; index < i; index++) {
            printf("x_v[%d] = %d\n", index, x_v[index]);
            sum = x_v[index] + sum;
        }
//            printf("%d\n", sum);
            printf("\n");

        //CALCULATING F(D(J))
        F_D_j[cols] = calculate_F_D_j(x_v, j, cols);
//        printf("F_D_j[%d] = %d\n", cols, F_D_j[cols]);
    }

    int** WH = allocate_2D_matrix(j,j);
    WH = generate_WHM(4);
    for(int row = 0; row < (1 << 3); row++){
        for(int col = 0; col < (1 << 3); col++){
//            printf("%d \n", WH[row][col]);
        }
    }

    free_2D_array(ptr_sbox_matrix, i);
    free_2D_array(v_i_j, i);
    free(temp);
    free(x_v);
    free_2D_array(WH, j);
    return 0;
} //END OF MAIN FUNCTION

//-----------------------------------------------------------------------------------------------------------
//DECLARATION OF FUNCTIONS
//FUNCTION DEFINITION OF 2D ALLOCATION
int** allocate_2D_matrix(int rows, int cols){
    int** array_2D = (int**) malloc(rows * sizeof(int));
    for(int index = 0; index < rows; index++){
        array_2D[index] = (int*) malloc(cols * sizeof(int));
    }
    return array_2D;
}

int** generateSboxMatrix(int row, int col){
    // Seed the random number generator with current time
//    srand(time(0));

    // Define the size of the matrix
    int rows = 256;
    int cols = 8;

    // Allocate memory for the matrix
    int** matrix = (int**)malloc(rows * sizeof(int*));
    for (int index = 0; index < rows; index++) {
        matrix[index] = (int*)malloc(cols * sizeof(int));
    }

    // Generate binary values for the matrix
    for (int row = 0; row < rows; row++) {
        for (int col = 0; col < cols; col++) {
            matrix[row][col] = generateBinaryValue();
        }
    }
    return matrix;
}
//CALCULATE WALSH HADAMARD MATRIX
//Function definition of generating Walsh Hadamard Matrix of order k

//FREE UP 2D-ARRAY
void free_2D_array(int** array_2D, int array_2D_length){
    size_t i;
    for(i = 0; i < array_2D_length; i++){
        if(array_2D != NULL){
            free(array_2D[i]);
        }
    }
    free(array_2D);
}

// Function to calculate the frequency of all elements in an array
int* calculate_x_v(int* A, int n){
    // create a count array of size `n` to store the count of all array elements
    int freq[n];
    int* x_v = (int*)malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) {
        freq[i] = 0;
    }

    // update frequency of each element
    for (int i = 0; i < n; i++) {
        freq[A[i]]++;
    }
    int counter = 0;
    // iterate through the array to print frequencies
    for (int i = 0; i < n; i++)
    {
        if (freq[i]) {
            //printf("%d: %d\n", i, freq[i]);
            x_v[i] = freq[i];
            counter++;
        }
    }
    return x_v;
}

int calculate_F_D_j(int *x_k_j, int n, int cols){
    int F_D_j = 0;
    int sum = 0;
    int temp;
//        printf("j = %d\n", cols );
        for(int k = 0; k < (2 << (cols + 1) - 1); k++){
//            printf("#x_k_j[%d]) = %d\n", k, x_k_j[k]);
//            printf("2^%d - 1 - %d = %d\n", n, cols, 1 << (n - 1 - cols));
            temp = abs(x_k_j[k]  - (1 << (n-1-cols)));
            sum = temp + sum;
            F_D_j = -sum;
//            printf("F_D_j[%d] = %d\n", k,temp);
//            printf("\n");
        }
//        printf("F((D(%d) = %d\n",cols,  F_D_j);
//        printf("\n");
    return F_D_j;
}

int** generate_WHM(int order){
    int n = 1 << order;

    int** hadamard = (int**) malloc(n * sizeof(int*));
    for (int row = 0; row < n; row++){
        hadamard[row] = (int*) malloc(n * sizeof(int));
    }

    hadamard[0][0] = 1;

    // Loop to copy elements to other quarters of the matrix
    for(int k = 1; k < n; k += k){

        for(int i = 0; i < k; i++){
            for(int j = 0; j < k; j++){
                hadamard[i+k][j] = hadamard[i][j];
                hadamard[i][j+k] = hadamard[i][j];
                hadamard[i+k][j+k] = -hadamard[i][j];
            }
        }
    }

    //Displaying final hadamard matrix
    for(int i = 0; i < n;i++){
        for(int j = 0; j < n; j++){
//            printf("%d ", hadamard[i][j]);
        }
//        printf("\n");
    }

    return hadamard;
}