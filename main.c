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
int **v_i_j;
int** allocate_2D_matrix(int rows, int cols);
int** generateSboxMatrix(int rows, int cols);
int** generate_WHM(int k);
void free_2D_array(int** array_2D, int array_2D_length);

int main() {

    int** ptr_sbox_matrix = allocate_2D_matrix(i, j);
    ptr_sbox_matrix = generateSboxMatrix(i, j);

    // Print the generated matrix
    printf("Generated Matrix:\n");
    for (int rows = 0; rows < i; rows++) {
        for (int cols = 0; cols < j; cols++) {
            printf("%d",  ptr_sbox_matrix[rows][cols]);
        }
        printf("\n");
    }

    //j = 0
    //v_0_0 = a_0_0 * 2^0
    //v_1_0 = a_1_0 * 2^0
    //...
    //v_254_0 = a_254_0 * 2_0
    //v_255_0 = a_255_0 * 2_0

    //j = 1
    //v_0_1 = (a_0_1 * 2^1) + (a_0_0 + 2^0)
    //v_1_1 = (a_1_1 * 2^1) + (a_1_0 + 2^0)
    //...
    //v_254_1 = (a_254_1 * 2^1) + (a_254_0 * 2^0)
    //v_255_1 = (a_255_1 * 2^1) + (a_255_0 * 2^0)

    //j = 2
    //v_0_1 = (a_0_2 * 2^2) + (a_0_1 * 2^1) + (a_0_0 + 2^0 )
    //v_1_1 = (a_1_2 * 2^2) + (a_1_1 * 2^1) + (a_1_0 + 2^0 )
    //...
    //v_254_1 = (a_254_1 * 2^2) + (a_254_0 * 2^1) + (a_254 * 2^0)
    //v_255_1 = (a_255_1 * 2^2) + (a_255_0 * 2^1) + (a_255 * 2^0)

    //...
    //j = 7
    //v_0_7 = (a_0_7 * 2^7) + (a_0_6 * 2^6) + (a_0_5 + 2^5) + (a_0_4 + 2^4) + (a_0_3 + 2^3) + (a_0_2 + 2^2) + (a_0_1 + 2^1) + (a_0_0 + 2^0)
    //v_1_7 = (a_1_7 * 2^7) + (a_1_6 * 2^6) + (a_1_5 + 2^5) + (a_1_4 + 2^4) + (a_1_3 + 2^3) + (a_1_2 + 2^2) + (a_1_1 + 2^1) + (a_1_0 + 2^0)
    //...
    //v_254_7 = (a_254_7 * 2^7) + (a_254_6 * 2^6) + (a_254_5 + 2^5) + (a_254_4 + 2^4) + (a_254_3 + 2^3) + (a_254_2 + 2^2) + (a_254_1 + 2^1) + (a_254_0 + 2^0)
    //v_255_7 = (a_255_7 * 2^7) + (a_255_6 * 2^6) + (a_255_5 + 2^5) + (a_255_4 + 2^4) + (a_255_3 + 2^3) + (a_255_2 + 2^2) + (a_255_1 + 2^1) + (a_255_0 + 2^0)

    // X_v_j
    // int* X_v = (int*) malloc(i * sizeof(int));
    int offset = 0;
    v_i_j = allocate_2D_matrix(i,j);
    int sum = 0;
    for(int rows = 0; rows < i; rows++){
        offset = 1;
        int sum = 0;
        for(int cols = j-1; cols >=0; cols--){
            sum = ptr_sbox_matrix[rows][cols] * 1 << offset-1;
            v_i_j[rows][offset-1] = sum;
            offset++;
        }
    }

    printf("cols = 0\n");
    for(int cols = 0; cols  < 1; cols++){
        for(int rows = 0; rows < 256; rows++){
            printf("v_i_j[%d][%d] = %d\n", cols, rows, v_i_j[rows][cols]);
        }
    }

    //PRINTING THE WEIGHTED VALUES OF THE BITS
    printf("v_i_j[%d][%d] = %d\n", 0, 0, v_i_j[0][0]);
    for(int cols = 1; cols < j; cols++){
        printf("cols = %d\n", cols);
        for(int rows = 0; rows < i; rows++){
            v_i_j[rows][cols] = v_i_j[rows][cols] + v_i_j[rows][cols-1];
            printf("v_i_j[%d][%d] = %d\n", rows, cols, v_i_j[rows][cols]);
        }
    }


//
//    for(int cols = 0; cols < j; cols++){
//        printf("j = %d\n", cols);
//        for(int rows = 0; rows < i; rows++){
//            printf("v_i_j[%d][%d] = %d\n", rows, cols, v_i_j[rows][cols]);
//        }
//    }

    //ACCUMULATIONS
    //CALCULATE BIJECTIVE FITNESS FUNCTION
    //Free S-box
    free_2D_array(ptr_sbox_matrix, i);

    return 0;
} //End of mail connection


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
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = generateBinaryValue();
        }
    }
    return matrix;

}
//CALCULATE WALSH HADAMARD MATRIX
//Function definition of generating Walsh Hadamard Matrix of order k
int** generate_WHM(int k) {
//    int **hadamard = memory_allocation_for_2D_array(GENE_SIZE, GENE_SIZE);
    int rows;
    int cols;
    int** hadamard = (int**) malloc(rows * sizeof(int*));
    for(int index = 0; index < rows; index++){
        hadamard[index] = (int*)malloc(cols * sizeof(int));
    }

    hadamard[0][0] = 1;
    for (int k = 1; k < 256; k += k) {
        // Loop to copy elements to
        // other quarters of the matrix
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < k; j++) {
                hadamard[i + k][j] = hadamard[i][j];
                hadamard[i][j + k] = hadamard[i][j];
                hadamard[i + k][j + k] = -hadamard[i][j];
            }
        }
    }

    // Displaying the final hadamard matrix
//    for (int i = 0; i < 256; i++) {
//        for (int j = 0; j < 256; j++) {
//            printf(" %d ",hadamard[i][j]);
//        }
//        printf("\n");
//    }

    return hadamard;
}

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

