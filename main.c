#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int i = 256;
int j = 8;
int** temp;

struct W_j_object{
    int* array;
    int array_size;
};

struct array_object{
    int* array;
    int array_size;
};

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

struct W_j_object W_j_1_plus;
struct W_j_object W_j_1_minus;
struct W_j_object W_j_2_plus;
struct W_j_object W_j_2_minus;
struct W_j_object W_j_3_plus;
struct W_j_object W_j_3_minus;

struct W_j_object calculate_W_1_plus(int* F_omaga, int size, int WH_max);
struct W_j_object calculate_W_1_minus(int* F_omega, int size, int WH_max);
struct W_j_object calculate_W_2_plus(int* F_omega, int size, int WH_max);
struct W_j_object calculate_W_2_minus(int* F_omega, int size, int WH_max);
struct W_j_object calculate_W_3_plus(int* F_omega, int size, int WH_max);
struct W_j_object calculate_W_3_minus(int* F_omega, int size, int WH_max);
//Function definition of generating Walsh Hadamard Matrix of order k
int main() {
    int** ptr_sbox_matrix = allocate_2D_matrix(i, j);
    ptr_sbox_matrix = generateSboxMatrix(i, j);
    int* f_j = (int*) malloc(i * sizeof(int));

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
    int* V = (int*) malloc(j + sizeof(int));

    //ASSIGNING f_j
    for(int cols = 0; cols < j; cols++){
        printf("j = %d\n", cols);
        for(int index = 0; index < i; index++){
            f_j[index] = ptr_sbox_matrix[index][7-cols];
//            printf("f_j[%d] = %d\n", index, f_j[index]);
        }
    }

    int** WH = allocate_2D_matrix(j,j);
    WH = generate_WHM(3);
    for(int row = 0; row < (1 << 3); row++){
        for(int col = 0; col < (1 << 3); col++){
            printf("%d ", WH[row][col]);
        }
        printf("\n");
    }

    //HADAMARD MATRIX®
    //    S_j_omega = matrix_multiplication(f, hadamard);
    int* F_omega = (int*)malloc(16 * sizeof(int));
//    F_omega = S_j_omega;

//    printf("Walsh-Hadamard transform for f Boolean function:\n");
    for(int omega = 0; omega < 8; omega++){
//        printf("%d ", S_j_omega[omega]);
    }

    for (int cols = 0; cols < j; cols++) {
        printf("j = %d\n", cols);


        for (int rows = 0; rows < i; rows++) {
            temp[rows] = v_i_j[rows][cols];
            //printf("v_j[%d] = %d\n", rows, temp[rows]);
        }

        //CALCULATING #x_v(j)
        int sum = 0;
        int v = 1 << cols+1;
        x_v = calculate_x_v(temp, i);
        for (int index = 0; index < v; index++) {
            //printf("#{x_v[%d]} = %d\n", index, x_v[index]);
            sum = x_v[index] + sum;
        }

        printf("2^(n - 1 - %d) = %d\n", cols, 1 << (8 - 1 - cols));
        int index_v = 0;
        for (int index = 0; index < v; index++) {
            printf("#{x_v[%d]} = %d\n", index, x_v[index]);
            if(x_v[index] > (1 << (8 - 1 - cols))){
                V[index_v] = index;
                //printf("V[%d] = %d\n", index_v, V[index_v]);
                index_v++;
            }
        }
        printf("\n");

        for(int index = 0; index < index_v; index++){
            printf("V[%d] = %d\n", index, V[index]);
        }

        //Calling W_1_plus
        printf("W_j_1_plus:\n");
        printf("W_j_1_plus.array_size = %d\n", W_j_1_plus.array_size);
        for(int index = 0; index < W_j_1_plus.array_size; index++){
            //printf("W_j_1_plus[%d] = %d\n", index, W_j_1_plus.array[index]);
        }
        printf("\n");
        printf("W_j_1_minus.array_size = %d\n", W_j_1_minus.array_size);
        printf("W_j_1_minus:\n");
        for(int index = 0; index < W_j_1_minus.array_size; index++){
            //printf("W_j_1_minus[%d] = %d\n", index, W_j_1_minus.array[index]);
        }

        //W_j_1 = calculate_union(W_j_1_plus.array, W_j_1_minus.array, W_j_1_plus.array_size, W_j_1_minus.array_size);
        //printf("\n");
        //printf("W_j_1.array_size = %d\n", W_j_1.array_size);
        //printf("W_j_1:\n");
        //for(int index = 0; index < W_j_1.array_size; index++){
        //  printf("W_j_1[%d] = %d\n", index, W_j_1.array[index]);
        //}
        printf("\n");
        //****************************************************************************************
        printf("W_j_1_plus:\n");
        printf("W_j_1_plus.array_size = %d\n", W_j_1_plus.array_size);
        for(int index = 0; index < W_j_1_plus.array_size; index++){
        //printf("W_j_1_plus[%d] = %d\n", index, W_j_1_plus.array[index]);
        }

        printf("\n");
        printf("W_j_2_plus:\n");
        printf("W_j_2_plus.array_size = %d\n", W_j_2_plus.array_size);
        for(int index = 0; index < W_j_2_plus.array_size; index++){
        //printf("W_j_2_plus[%d] = %d\n", index, W_j_2_plus.array[index]);
        }
        printf("\n");
        //W_j_1_2_plus
        //W_j_1_2_plus = calculate_union(W_j_1_plus.array, W_j_2_plus.array, W_j_1_plus.array_size, W_j_2_plus.array_size);
        //printf("W_j_1_2_plus:\n");
        //printf("W_j_1_2_plus.array_size = %d\n", W_j_1_2_plus.array_size);
        //for(int index = 0; index < W_j_1_2_plus.array_size; index++){
        //printf("W_j_1_2_plus[%d] = %d\n", index, W_j_1_2_plus.array[index]);
        //}
        //printf("\n");
        //*******************************************************************************
        //Calling W_j_1_minus
        printf("W_j_1_minus:\n");
        printf("W_j_1_minus.array_size = %d\n", W_j_1_minus.array_size);
        for(int index = 0; index < W_j_1_minus.array_size; index++){
            printf("W_j_1_minus[%d] = %d\n", index, W_j_1_minus.array[index]);
        }
        printf("\n");
        //Calling W_2_minus
        printf("W_j_2_minus:\n");
        printf("W_j_2_plus.array_size = %d\n", W_j_2_minus.array_size);
        for(int index = 0; index < W_j_2_minus.array_size; index++){
            printf("W_j2_minus[%d] = %d\n", index, W_j_2_minus.array[index]);
        }

        //W_j_1_2_minus
        //W_j_1_2_minus = calculate_union(W_j_1_minus.array, W_j_2_minus.array, W_j_1_minus.array_size, W_j_2_minus.array_size);
        //printf("\n");
        //printf("W_j_1_2_minus.array_size = %d\n", W_j_1_2_minus.array_size);
        //printf("W_j_1_2_minus\n");
        //for(int index = 0; index < W_j_1_2_minus.array_size; index++){
        //printf("W_j_1_2_minus[%d] = %d\n", index, W_j_1_2_minus.array[index]);
        //}
        // printf("\n");
        //**********************************************************************************
        printf("W_j_2_plus:\n");
        printf("W_j_2_plus.array_size = %d\n", W_j_2_plus.array_size);
        for(int index = 0; index < W_j_2_plus.array_size; index++){
            printf("W_j_2_plus[%d] = %d\n", index, W_j_2_plus.array[index]);
        }
        printf("\n");

        printf("W_j_3_plus:\n");
        printf("W_j_3_plus.array_size = %d\n", W_j_3_plus.array_size);
        for(int index = 0; index < W_j_3_plus.array_size; index++){
            printf("W_j_3_plus[%d] = %d\n", index, W_j_3_plus.array[index]);
        }

        //printf("\n");
        //printf("W_j_2_3_plus:\n");
        //W_j_2_3_plus = calculate_union(W_j_2_plus.array, W_j_3_plus.array, W_j_2_plus.array_size, W_j_3_plus.array_size);
        //printf("W_j_2_3_plus.array_size = %d\n", W_j_2_3_plus.array_size);
        //for(int index = 0; index < W_j_2_3_plus.array_size; index++){
        //printf("W_j_2_3_plus[%d] = %d\n", index, W_j_2_3_plus.array[index]);
        //}
        //printf("\n");
        //**************************************************************************************
        printf("W_j_2_minus:\n");
        printf("W_j_2_minus.array_size = %d\n", W_j_2_minus.array_size);
        for(int index = 0; index < W_j_2_minus.array_size; index++){
         printf("W_j_2_minus[%d] = %d\n", index, W_j_2_minus.array[index]);
        }
        printf("\n");
        printf("W_j_3_minus.array_size = %d\n", W_j_3_minus.array_size);
        printf("W_j_3_minus:\n");
        for(int index = 0; index < W_j_3_minus.array_size; index++){
            printf("W_j_3_minus[%d] = %d\n", index, W_j_3_minus.array[index]);
        }

        //W_j_2_3_minus = calculate_union(W_j_2_minus.array, W_j_3_minus.array,W_j_2_minus.array_size, W_j_3_minus.array_size);
        //printf("\n");
        //printf("W_j_2_3_minus.array_size = %d\n", W_j_2_3_minus.array_size);
        //printf("W_j_2_3_minus:\n");
        //for(int index = 0; index < W_j_2_3_minus.array_size; index++){
        //  printf("W_j_2_3_minus[%d] = %d\n", index, W_j_2_3_minus.array[index]);
        //}

        //CALCULATING F(D(J))
        F_D_j[cols] = calculate_F_D_j(x_v, j, cols);
        //  printf("F_D_j[%d] = %d\n", cols, F_D_j[cols]);
    }



    printf("\n");
//    int WH_max = calculate_WH_max(S_j_omega, 8);
    //CALCULATION OF W_j_plus, W_j_minus, W_1_plus, W_1_minus,  W_j_1_2_plus, W_j_1_2_minus, W_j_2
//    W_j_1_plus = calculate_W_1_plus(F_omega, 8, WH_max);
//    W_j_1_minus = calculate_W_1_minus(F_omega, 8, WH_max);
//    W_j_2_plus = calculate_W_2_plus(F_omega, 8, WH_max);
//    W_j_2_minus = calculate_W_2_minus(F_omega, 8, WH_max);
//    W_j_3_plus = calculate_W_3_plus(F_omega, 8, WH_max);
//    W_j_3_minus = calculate_W_3_minus(F_omega, 8, WH_max);


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
    size_t index;
    for(index = 0; index < array_2D_length; index++){
        if(array_2D != NULL){
            free(array_2D[index]);
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

struct W_j_object calculate_W_1_plus(int* F_omaga, int size, int WH_max){
    struct W_j_object W_1_plus;
    int index = 0;
    int counter = 0;
    W_1_plus.array = (int*)malloc(8 * sizeof(int));
    W_1_plus.array_size = 0;
//    printf("WH_max = %d\n", WH_max);
    for(int omega = 0; omega < size; omega++){
        if (F_omaga[omega] == WH_max) {
            W_1_plus.array[index++] = omega;
            counter++;
        }
    }
    W_1_plus.array_size = counter;
    return W_1_plus;
}

//w_1_2_plus
struct W_j_object calculate_W_1_minus(int* F_omaga, int size, int WH_max){
    int index = 0;
    int counter = 0;
    struct W_j_object W_1_minus;
    W_1_minus.array = (int*)malloc(8 * sizeof(int));
//    printf("-WH_max = %d\n", -WH_max);
    for(int omega = 0; omega < size; omega++){
        if (F_omaga[omega] == -WH_max) {
            W_1_minus.array[index++] = omega;
            counter++;
        }
    }
    W_1_minus.array_size = counter;
    return W_1_minus;
}

//w_2_plus
struct W_j_object calculate_W_2_plus(int* F_omega, int size, int WH_max){
    int index = 0;
    int counter = 0;
    struct W_j_object W_2_plus;
    W_2_plus.array = (int*) malloc(8 * sizeof(int));
//    printf("WH_max - 2 = %d\n", WH_max - 2);
    for(int omega = 0; omega < size; omega++){
        if(F_omega[omega] == WH_max - 2){
            W_2_plus.array[index++] = omega;
            counter++;
        }
    }
    W_2_plus.array_size = counter;
    return  W_2_plus;
}

//w_2_minus
struct W_j_object calculate_W_2_minus(int* F_omega, int size, int WH_max){
    int index = 0;
    int counter = 0;
    struct W_j_object W_2_minus;
//    printf("-WH_max + 2 = %d\n", -WH_max + 2);
    for(int omega = 0; omega < size; omega++){
        if(F_omega[omega] == -WH_max + 2){
            W_2_minus.array[index++] = omega;
            counter++;
        }
    }
    W_2_minus.array_size = counter;
    return W_2_minus;
}

//w_3_plus
struct W_j_object calculate_W_3_plus(int* F_omega, int size, int WH_max){
    int index = 0;
    int counter = 0;
    struct W_j_object W_3_plus;
//    printf("WH_max - 4 = %d\n", WH_max - 4);
    W_3_plus.array = (int*)malloc(8 * sizeof(int));
    for(int omega = 0; omega < size; omega++) {
        if (F_omega[omega] == WH_max - 4) {
            W_3_plus.array[index++] = omega;
            counter++;
        }
    }
    W_3_plus.array_size = counter;
    return W_3_plus;
}

//w_3_minus
struct W_j_object calculate_W_3_minus(int* F_omega, int size, int WH_max){
    int index = 0;
    int counter = 0;
    struct W_j_object W_3_minus;
//    printf("-WH_max + 4 = %d\n", -WH_max + 4);
    W_3_minus.array = (int*) malloc(8 * sizeof(int));
    for(int omega = 0; omega < size; omega++){
        if(F_omega[omega] == (-WH_max) + 4){
            W_3_minus.array[index++] = omega;
            counter++;
        }
    }
    W_3_minus.array_size = counter;
    return W_3_minus;
}