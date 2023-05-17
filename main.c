#include <stdio.h>
#include <stdlib.h>
#include <time.h>

struct W_j_object{
    int* array;
    int array_size;
};

struct array_object{
    int* array;
    int array_size;
};
int n = 8;
int i;
int j;
int** v_i_j;
int* v_j;
struct W_j_object W_j_1_plus;
struct W_j_object W_j_1_minus;
struct W_j_object W_j_2_plus;
struct W_j_object W_j_2_minus;
struct W_j_object W_j_3_plus;
struct W_j_object W_j_3_minus;

struct array_object W_j_1;
struct array_object W_j_1_2_plus;
struct array_object W_j_1_2_minus;
struct array_object W_j_2_3_plus;
struct array_object W_j_2_3_minus;

//FUNCTION SIGNATURES
int generateBinaryValue();
int** allocate_2D_matrix(int rows, int cols);
int** generateSboxMatrix(int rows, int cols);
int** D;
void free_2D_array(int** array_2D, int array_2D_length);
int* calculate_X_v(int* A, int n);
int calculate_F_D_j(int *x_k_j, int rows, int cols);
int** generate_WHM(int order);
int calculate_WH_max (int* walsh_spectrum, int size);
struct W_j_object calculate_W_1_plus(int* F_omaga, int size, int WH_max);
struct W_j_object calculate_W_1_minus(int* F_omega, int size, int WH_max);
struct W_j_object calculate_W_2_plus(int* F_omega, int size, int WH_max);
struct W_j_object calculate_W_2_minus(int* F_omega, int size, int WH_max);
struct W_j_object calculate_W_3_plus(int* F_omega, int size, int WH_max);
struct W_j_object calculate_W_3_minus(int* F_omega, int size, int WH_max);
int* matrix_multiplication(int *array, int **matrix);
struct array_object calculate_union(int* array_1, int* array_2, int array_size_1, int array_size_2);
int complement_f_j_i(int f_j_i);
void print_Sbox(int **S_box, int row_number, int col_number);
int** calculate_v_i_j(int** S_box, int row_number, int col_number);

int main() {
    i = 1 << n;
    j = n;
    int **ptr_sbox_matrix = allocate_2D_matrix(i, j);
    D = allocate_2D_matrix(i, j);
    ptr_sbox_matrix = generateSboxMatrix(i, j);
    int **f_j = allocate_2D_matrix(i, j);
    int *S_j_omega = (int *) malloc(i * sizeof(int));
    int offset = 0;
    v_i_j = allocate_2D_matrix(i, j);
    v_j = (int *) malloc(j * sizeof(int));
    int *X_v = (int *) malloc(i * sizeof(int));
    int *updated_X_v = (int *) malloc(i * sizeof(int));
    int *F_D_j = (int *) malloc(j * sizeof(int));
    int *f = (int *) malloc(i * sizeof(int));
    int *F_omega = (int *) malloc(16 * sizeof(int));
    int *V = (int *) malloc(j + sizeof(int));
    int WH_max;
    int *updated_D = (int *) malloc(j * sizeof(int));

    //GENERATING WALSH HADAMARD MATRIX OF ORDER j
    int **WH = allocate_2D_matrix(j, j);
    WH = generate_WHM(j);

    //PRINTING WALSH HADAMARD MATRIX
    for (int row = 0; row < (1 << j); row++) {
        for (int col = 0; col < (1 << j); col++) {
            //printf("%d ", WH[row][col]);
        }
        //printf("\n");
    }

    //ASSIGNING f_j where f_0 is MSB and f_(j-1) is MSB as represented by the array
    for (int col = 0; col < j; col++) {
        printf("j = %d\n", col);
        for (int row = 0; row < i; row++) {
            //CONSTRUCT D(j) FROM FROM THE RANDOM SBOX WHERE j = 0, 1, 2, ..., n-1
            D[row][col] = ptr_sbox_matrix[row][7 - col];
            f[row] = D[row][col];
            //printf("D[%d][%d] = %d\n", row, col, D[row][col]);

            //printf("f_j[%d][%d] = %d\n", row, col, f_j[row][col]);
            //printf("f[%d] = %d\n", row, f[row]);
        }

        //CALCULATION OF v_i_j USING FUNCTION
        v_i_j = calculate_v_i_j(D, i, col);
//        printf("j = %d\n", col);
        int sum = 0;
        for (int row = 0; row < i; row++) {
            v_j[row] = v_i_j[row][col];
            //printf("v_j[%d] = %d\n", row, v_j[row]);
        }

        //CALCULATION OF x_v and F(D(j)
        //CALCULATING OF #x_v(j)
        //X_v CALCULATES REPETITION OF THE VALUE v IN SUB-MATRIX D(j) FOR j = 0, 1, 2, ... , n-1
        X_v = calculate_X_v(v_j, i);
        int index_v = 0;
        int v = 1 << (col + 1) - 1;
        for (int index = 0; index <= v; index++) {
            //printf("#{X_v[%d]} = %d\n", index, X_v[index]);
            //printf("2^n-1-%d  = %d\n", col, 1 << (n - 1 - col));
            if (X_v[index] > (1 << (n - 1 - col))) {
                V[index_v++] = index;
                //printf("%d > %d\n", X_v[index], 1 << (n - 1 - col));
                //printf("v = %d\n", index);
            }
            //printf("\n");
            //sum = X_v[index] + sum;
        }
        //printf("index_v = %d\n", index_v);
        //CALCULATING BIJECTIVE FITNESS FUNCTION, F(D(j))
        F_D_j[col] = calculate_F_D_j(X_v, j, col);
        //printf("F_D_j[%d] = %d\n", col, F_D_j[col]);
        //printf("\n");

        v = (1 << (col + 1)) - 1;
        printf("v = %d\n", v);
        for (int index = 0; index <= v; index++) {
            //printf("#{X_v[%d]} = %d\n", index, X_v[index]);
            //printf("2^n-1-%d  = %d\n", col, 1 << (n - 1 - col));
            if (X_v[index] > (1 << (n - 1 - col))) {
                V[index_v++] = index;
                //printf("%d > %d\n", X_v[index], 1 << (n - 1 - col));
                //printf("v = %d\n", index);
            }
            //printf("\n");
            //sum = X_v[index] + sum;
        }

        for(int index = 0; index < index_v; index++) {
            printf("V[%d] = %d\n", index, V[index]);
        }


        //APPLYING WALSH HADAMARD MATRIX ON f_j TO CALCULATE WALSH HADAMARD SPECTRUM OF f_j
        S_j_omega = matrix_multiplication(f, WH);
        F_omega = S_j_omega;

        // printf("j = %d\n", col);
        //printf("Walsh-Hadamard transform for f Boolean function:\n");
        for (int omega = 0; omega < i; omega++) {
            //printf("%d ", S_j_omega[omega]);
            //printf("f[%d] = %d\n", omega,  f[omega]);
            //printf("S_j_omega[%d] = %d\n", omega, S_j_omega[omega]);
            //printf("f[%d] = %d\n", omega,  F_omega[omega]);
        }
        WH_max = calculate_WH_max(S_j_omega, i);
        //printf("WH_max = %d\n", WH_max);

        W_j_1_plus = calculate_W_1_plus(S_j_omega, i, WH_max);
        W_j_1_minus = calculate_W_1_minus(S_j_omega, i, WH_max);
        W_j_2_plus = calculate_W_2_plus(S_j_omega, i, WH_max);
        W_j_2_minus = calculate_W_2_minus(S_j_omega, i, WH_max);
        W_j_3_plus = calculate_W_3_plus(S_j_omega, i, WH_max);
        W_j_3_minus = calculate_W_3_minus(S_j_omega, i, WH_max);
        printf("\n");
    }

    //CALCULATION OF #X_v(j) AND BIJECTIVE FITNESS FUNCTION, F(D(j)
    for (int col = 0; col < j; col++) {
//        v_i_j = calculate_v_i_j(D, i, col);
//        //printf("j = %d\n", col);
//        int sum = 0;
//        for (int row = 0; row < i; row++) {
//            v_j[row] = v_i_j[row][col];
//            //printf("v_j[%d] = %d\n", row, v_j[row]);
//        }
#if 0
        //CALCULATING OF #x_v(j)
        //X_v CALCULATES REPETITION OF THE VALUE v IN SUB-MATRIX D(j) FOR j = 0, 1, 2, ... , n-1
        X_v = calculate_X_v(v_j, i);
        int index_v = 0;
        int v = 1 << (col + 1) - 1;
        for (int index = 0; index <= v; index++) {
            //printf("#{X_v[%d]} = %d\n", index, X_v[index]);
            //printf("2^n-1-%d  = %d\n", col, 1 << (n - 1 - col));
            if (X_v[index] > (1 << (n - 1 - col))) {
                V[index_v++] = index;
                //printf("%d > %d\n", X_v[index], 1 << (n - 1 - col));
                //printf("v = %d\n", index);
            }
            //printf("\n");
            //sum = X_v[index] + sum;
        }
        //printf("index_v = %d\n", index_v);
        //CALCULATING BIJECTIVE FITNESS FUNCTION, F(D(j))
        F_D_j[col] = calculate_F_D_j(X_v, j, col);
        //printf("F_D_j[%d] = %d\n", col, F_D_j[col]);
        //printf("\n");
#if 1
        //STEP 5: Set a=V(m),here V(m)is the m-th element in V int a;
        int a;
        // printf("index_v = %d\n", index_v);
        for(int m = 0; m < index_v; m++){
            a = V[m];
            printf("a = %d\n", a);
            for(col = 0; col < j; col++){
                int count_X_v_1 = 0;
                for(int row = 0; row < i; row++){
                    if(a == v_i_j[row][col]){
                        //printf("rows_with_1 = %d\n", row);
                       count_X_v_1++;
                        //printf("f_j[%d][%d] = %d\n", row, col, f_j[row][col]);
                        //printf("count_X_v_1 = %d\n", count_X_v_1);
                        if(count_X_v_1 > (1 << (n-1-col))){
                            D[row][col] = complement_f_j_i(f_j[row][col]);
                        }
                        //printf("D[%d][%d] = %d\n", row, col, D[row][col]);
                        //printf("\n");
                    }
                }
            }
            //printf("V[%d] = %d\n", m, V[m]);
        }
    }
//    printf("Printing Updated_X_v\n");
//    for(int col = 0; col < j; col++){
//        printf("j = %d\n", col);
//        for(int row = 0; row < i; row++){
//            updated_D[row] = D[row][col];
//        }
//        updated_X_v = calculate_X_v(updated_D, i);
//        for(int v = 0; v < (1 << col+1); v++){
//            printf("updated_X_v[%d] = %d\n", v, updated_X_v[v]);
//        }
//    }

        //######################
//PRINTING X_v(j) FOR THE UPDATED
    for (int rows = 0; rows < i; rows++) {
        offset = 1;
        int sum = 0;
        //COLUMN NUMBER STARTS FROM THE LSB: j = 0 IS LSB, j = n-1 IS THE MSB
        //ONLY CALCULATE THE SUM FOR EACH ROW WITH SET BITS
        for (int cols = 0; cols < n; cols++) {
            sum = D[rows][cols] * 1 << (offset - 1);
            v_i_j[rows][offset - 1] = sum;
            //printf("v_i_j[%d][%d] = %d\n", rows, offset-1, v_i_j[rows][offset-1]);
            offset++;
        }
    }
    for (int col = 1; col < j; col++) {
        for (int rows = 0; rows < i; rows++) {
            D[rows][col] = D[rows][col] + D[rows][col - 1];
            //printf("v_i_j[%d][%d] = %d\n", rows, col, v_i_j[rows][col]);
        }
    }

    //CALCULATION OF #X_v(j) AND BIJECTIVE FITNESS FUNCTION, F(D(j)
    for(int col = 0; col < j; col++) {
        //printf("j = %d\n", col);
        int sum = 0;
        for (int row = 0; row < i; row++) {
            v_j[row] = v_i_j[row][col];
            //printf("v_j[%d] = %d\n", row, v_j[row]);
        }

        //CALCULATING OF #x_v(j)
        //X_v CALCULATES REPETITION OF THE VALUE v IN SUB-MATRIX D(j) FOR j = 0, 1, 2, ... , n-1
        X_v = calculate_X_v(v_j, i);
        int index_v = 0;
        int v = 1 << (col + 1);
        for (int index = 0; index < v; index++) {
            //printf("#{X_v[%d]} = %d\n", index, X_v[index]);
            //printf("2^n-1-%d  = %d\n", col, 1 << (n - 1 - col));

            //printf("%d > %d\n", X_v[index], 1 << (n - 1 - col));
            //printf("v = %d\n", index);
        }

        //CALCULATING BIJECTIVE FITNESS FUNCTION, F(D(j))
//        F_D_j[col] = calculate_F_D_j(X_v, j, col);
        //printf("F_D_j[%d] = %d\n", col, F_D_j[col]);
        printf("\n");
    }
        //####################################################################
        for (int cols = 0; cols < j; cols++) {
            //printf("j = %d \n", cols);
            //printf("f_%d = \n", cols);

            for (int rows = 0; rows < i; rows++) {
                //printf("f_j_i[%d] = %d\n", rows, f[rows]);

                //TODO:
                //Check the condition

                //Add complement f_i_j
                int complement_f = complement_f_j_i(f[rows]);
                //        printf("complement_f_j_i[%d] = %d\n", rows, complement_f);
            }

            //    printf("\n");
            //    for(int omega = 0; omega < i; omega++){
            //        printf("%d ", F_omega[omega]);
            //    }

            //STEP 3: According to the Boolean function f_j in D(j),calculate the sets W_j_plus, W_j_minus, W_1_plus, W_1_minus, W_j_1_2_plus, W_j_1_2_minus, W_j_2
            //            printf("\n");
//            int WH_max = calculate_WH_max(S_j_omega, i);
//            printf("WH_max = %d\n", WH_max);

//            W_j_1_plus = calculate_W_1_plus(F_omega, i, WH_max);
//            W_j_1_minus = calculate_W_1_minus(F_omega, i, WH_max);
//            W_j_2_plus = calculate_W_2_plus(F_omega, i, WH_max);
//            W_j_2_minus = calculate_W_2_minus(F_omega, i, WH_max);
//            W_j_3_plus = calculate_W_3_plus(F_omega, i, WH_max);
//            W_j_3_minus = calculate_W_3_minus(F_omega, i, WH_max);

            //Calling W_1_plus
            //    printf("W_j_1_plus:\n");
            //    printf("W_j_1_plus.array_size = %d\n", W_j_1_plus.array_size);

            for (int index = 0; index < W_j_1_plus.array_size; index++) {
                //printf("W_j_1_plus[%d] = %d\n", index, W_j_1_plus.array[index]);
            }

            //printf("\n");
            //    printf("W_j_1_minus.array_size = %d\n", W_j_1_minus.array_size);
            //    printf("W_j_1_minus:\n");
            for (int index = 0; index < W_j_1_minus.array_size; index++) {
                //        printf("W_j_1_minus[%d] = %d\n", index, W_j_1_minus.array[index]);
            }

            //W_j_1 = calculate_union(W_j_1_plus.array, W_j_1_minus.array, W_j_1_plus.array_size, W_j_1_minus.array_size);
            //printf("\n");
            //    printf("W_j_1.array_size = %d\n", W_j_1.array_size);
            //    printf("W_j_1:\n");
            for (int index = 0; index < W_j_1.array_size; index++) {
                //        printf("W_j_1[%d] = %d\n", index, W_j_1.array[index]);
            }
            //    printf("**********************************************************************\n");
            //    printf("W_j_1_plus:\n");
            //    printf("W_j_1_plus.array_size = %d\n", W_j_1_plus.array_size);
            for (int index = 0; index < W_j_1_plus.array_size; index++) {
                //        printf("W_j_1_plus[%d] = %d\n", index, W_j_1_plus.array[index]);
            }

            //    printf("\n");
            //    printf("W_j_2_plus:\n");
            //    printf("W_j_2_plus.array_size = %d\n", W_j_2_plus.array_size);
            for (int index = 0; index < W_j_2_plus.array_size; index++) {
                //        printf("W_j_2_plus[%d] = %d\n", index, W_j_2_plus.array[index]);
            }

            //    printf("\n");
            //W_j_1_2_plus
            //    printf("W_j_1_2_plus:\n");
            //W_j_1_2_plus = calculate_union(W_j_1_plus.array, W_j_2_plus.array, W_j_1_plus.array_size,
            //W_j_2_plus.array_size);
            //    printf("W_j_1_2_plus.array_size = %d\n", W_j_1_2_plus.array_size);
            for (int index = 0; index < W_j_1_2_plus.array_size; index++) {
                //        printf("W_j_1_2_plus[%d] = %d\n", index, W_j_1_2_plus.array[index]);
            }
            //    printf("**********************************************************************\n");
            //*******************************************************************************
            //Calling W_j_1_minus
            //    printf("W_j_1_minus:\n");
            //    printf("W_j_1_minus.array_size = %d\n", W_j_1_minus.array_size);
            for (int index = 0; index < W_j_1_minus.array_size; index++) {
                //        printf("W_j_1_minus[%d] = %d\n", index, W_j_1_minus.array[index]);
            }

            //    printf("\n");
            //Calling W_2_minus
            //    printf("W_j_2_minus:\n");
            //    printf("W_j_2_minus.array_size = %d\n", W_j_2_minus.array_size);
            for (int index = 0; index < W_j_2_minus.array_size; index++) {
                //        printf("W_j_2_minus[%d] = %d\n", index, W_j_2_minus.array[index]);
            }
            //    printf("\n");

            //    printf("W_j_1_2_minus\n");
            //W_j_1_2_minus = calculate_union(W_j_1_minus.array, W_j_2_minus.array, W_j_1_minus.array_size,
            //W_j_2_minus.array_size);
            //    printf("W_j_1_2_minus.array_size = %d\n", W_j_1_2_minus.array_size);
            for (int index = 0; index < W_j_1_2_minus.array_size; index++) {
                //        printf("W_j_1_2_minus[%d] = %d\n", index, W_j_1_2_minus.array[index]);
            }
            //    printf("**********************************************************************\n");
            //**********************************************************************************
            //    printf("W_j_2_plus:\n");
            //    printf("W_j_2_plus.array_size = %d\n", W_j_2_plus.array_size);
            for (int index = 0; index < W_j_2_plus.array_size; index++) {
                //        printf("W_j_2_plus[%d] = %d\n", index, W_j_2_plus.array[index]);
            }
            //    printf("\n");
            //    printf("W_j_3_plus:\n");
            //    printf("W_j_3_plus.array_size = %d\n", W_j_3_plus.array_size);
            for (int index = 0; index < W_j_3_plus.array_size; index++) {
                //        printf("W_j_3_plus[%d] = %d\n", index, W_j_3_plus.array[index]);
            }

            //    printf("\n");
            //    printf("W_j_2_3_plus:\n");
            //W_j_2_3_plus = calculate_union(W_j_2_plus.array, W_j_3_plus.array, W_j_2_plus.array_size,
            //W_j_3_plus.array_size);
            //    printf("W_j_2_3_plus.array_size = %d\n", W_j_2_3_plus.array_size);
            for (int index = 0; index < W_j_2_3_plus.array_size; index++) {
                //        printf("W_j_2_3_plus[%d] = %d\n", index, W_j_2_3_plus.array[index]);
            }
            //    printf("**********************************************************************\n");
            //**************************************************************************************printf("W_j_2_minus:\n");
            //    printf("W_j_2_minus.array_size = %d\n", W_j_2_minus.array_size);
            for (int index = 0; index < W_j_2_minus.array_size; index++) {
                //        printf("W_j_2_minus[%d] = %d\n", index, W_j_2_minus.array[index]);
            }

            //    printf("\n");
            //    printf("W_j_3_minus.array_size = %d\n", W_j_3_minus.array_size);
            //    printf("W_j_3_minus:\n");
            for (int index = 0; index < W_j_3_minus.array_size; index++) {
                //        printf("W_j_3_minus[%d] = %d\n", index, W_j_3_minus.array[index]);
            }

            //W_j_2_3_minus = calculate_union(W_j_2_minus.array, W_j_3_minus.array, W_j_2_minus.array_size,
            //W_j_3_minus.array_size);
            //    printf("\n");
            //    printf("W_j_2_3_minus.array_size = %d\n", W_j_2_3_minus.array_size);
            //    printf("W_j_2_3_minus:\n");
            for (int index = 0; index < W_j_2_3_minus.array_size; index++) {
                //        printf("W_j_2_3_minus[%d] = %d\n", index, W_j_2_3_minus.array[index]);
            }
        }
#endif
        //printf("\n");

        //printf("Updated S-box\n");
        //print_Sbox(D, i, j);
#endif
    free_2D_array(v_i_j, i);
    free(v_j);
    free_2D_array(WH, i);
    v_j = NULL;

    free(S_j_omega);
    S_j_omega = NULL;

    free(f);
    f = NULL;

    free(F_D_j);
    F_D_j = NULL;

    free(V);
    V = NULL;

    free(W_j_1_plus.array);
    W_j_1_plus.array = NULL;

    free(W_j_1_minus.array);
    W_j_1_plus.array = NULL;

    free(W_j_2_plus.array);
    W_j_2_plus.array = NULL;

    free(W_j_2_minus.array);
    W_j_2_minus.array = NULL;

    free(W_j_3_plus.array);
    W_j_3_plus.array = NULL;

    free(X_v);
    X_v = NULL;

    free_2D_array(ptr_sbox_matrix, i);


    return 0;
} //END OF MAIN FUNCTION
}
//######################################################################################################################
//DECLARATION OF FUNCTIONS

// Function to generate a random binary value (0 or 1)
int generateBinaryValue() {
    return rand() % 2;
}

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
int* calculate_X_v(int* A, int row_number){
    // create a count array of size `n` to store the count of all array elements
    int freq[row_number];
    int* x_v = (int*)malloc(row_number * sizeof(int));
    for (int index = 0; index < row_number; index++) {
        freq[index] = 0;
    }

    // update frequency of each element
    for (int index = 0; index < row_number; index++) {
        freq[A[index]]++;
    }
    int counter = 0;
    // iterate through the array to print frequencies
    for (int index = 0; index < row_number; index++)
    {
        if (freq[index]) {
            //printf("%d: %d\n", i, freq[i]);
            x_v[index] = freq[index];
            counter++;
        }
    }

    return x_v;
}

int calculate_F_D_j(int *x_k_j, int rows, int cols){
    int F_D_j = 0;
    int sum = 0;
    int temp;
//        printf("j = %d\n", cols );
        for(int k = 0; k < (2 << (cols + 1) - 1); k++){
            //printf("#x_k_j[%d]) = %d\n", k, x_k_j[k]);
            //printf("2^%d - 1 - %d = %d\n", n, cols, 1 << (n - 1 - cols));
            temp = abs(x_k_j[k]  - (1 << (rows-1-cols)));
            sum = temp + sum;
            F_D_j = -sum;
            //printf("F_D_j[%d] = %d\n", k,temp);
            //printf("\n");
        }
//        printf("F((D(%d) = %d\n",cols,  F_D_j);
//        printf("\n");
    return F_D_j;
}

int** generate_WHM(int order){
    int n_WHM = 1 << order;
    int** hadamard = (int**) malloc(n_WHM * sizeof(int*));
    for (int row = 0; row < n_WHM; row++){
        hadamard[row] = (int*) malloc(n_WHM * sizeof(int));
    }
    hadamard[0][0] = 1;
    // Loop to copy elements to other quarters of the matrix
    for(int k = 1; k < n_WHM; k += k){
        for(int i = 0; i < k; i++){
            for(int j = 0; j < k; j++){
                hadamard[i+k][j] = hadamard[i][j];
                hadamard[i][j+k] = hadamard[i][j];
                hadamard[i+k][j+k] = -hadamard[i][j];
            }
        }
    }
    //Displaying final hadamard matrix
    for(int i = 0; i < n_WHM;i++){
        for(int j = 0; j < n_WHM; j++){
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
    W_1_plus.array = (int*)malloc(i * sizeof(int));
    W_1_plus.array_size = 0;
//    printf("WH_max = %d\n", WH_max);
    for(int omega = 0; omega < size; omega++){
        if (F_omaga[omega] == WH_max) {
            W_1_plus.array[index++] = omega;
            counter++;
        }
    }
    W_1_plus.array_size = counter;

    for(int omega = 0; omega < W_1_plus.array_size; omega++){
//        printf("W_1_plus[%d] = %d:\n", index, W_1_plus.array[omega]);
    }
    return W_1_plus;
}

//w_1_2_plus
struct W_j_object calculate_W_1_minus(int* F_omaga, int size, int WH_max){
    int index = 0;
    int counter = 0;
    struct W_j_object W_1_minus;
    W_1_minus.array = (int*)malloc(i * sizeof(int));
//    printf("-WH_max = %d\n", -WH_max);
    for(int omega = 0; omega < size; omega++){
        if (F_omaga[omega] == -WH_max) {
            W_1_minus.array[index++] = omega;
            counter++;
        }
    }
    W_1_minus.array_size = counter;

    for(int omega = 0; omega < W_1_minus.array_size; omega++){
//        printf("W_1_minus[%d] = %d:\n", omega, W_1_minus.array[omega]);
    }
    return W_1_minus;
}

//w_2_plus
struct W_j_object calculate_W_2_plus(int* F_omega, int size, int WH_max){
    int index = 0;
    int counter = 0;
    struct W_j_object W_2_plus;
    W_2_plus.array = (int*) malloc(i * sizeof(int));
//    printf("WH_max - 2 = %d\n", WH_max - 2);
    for(int omega = 0; omega < size; omega++){
        if(F_omega[omega] == WH_max - 2){
            W_2_plus.array[index++] = omega;
            counter++;
        }
    }
    W_2_plus.array_size = counter;

    for(int omega = 0; omega < W_2_plus.array_size; omega++){
//        printf("W_2_plus[%d] = %d:\n", omega, W_2_plus.array[omega]);
    }
    return  W_2_plus;
}

//w_2_minus
struct W_j_object calculate_W_2_minus(int* F_omega, int size, int WH_max){
    int index = 0;
    int counter = 0;
    struct W_j_object W_2_minus;
//    printf("-WH_max + 2 = %d\n", -WH_max + 2);
    W_2_minus.array = (int*)malloc(i * sizeof(int));
    for(int omega = 0; omega < size; omega++){
        if(F_omega[omega] == -WH_max + 2){
            W_2_minus.array[index++] = omega;
            counter++;
        }
    }
    W_2_minus.array_size = counter;
    for(int omega = 0; omega < W_2_minus.array_size; omega++){
//        printf("W_2_minus[%d] = %d:\n", omega, W_2_minus.array[omega]);
    }
    return W_2_minus;
}

//w_3_plus
struct W_j_object calculate_W_3_plus(int* F_omega, int size, int WH_max){
    int index = 0;
    int counter = 0;
    struct W_j_object W_3_plus;
//    printf("WH_max - 4 = %d\n", WH_max - 4);
    W_3_plus.array = (int*)malloc(i * sizeof(int));
    for(int omega = 0; omega < size; omega++) {
        if (F_omega[omega] == WH_max - 4) {
            W_3_plus.array[index++] = omega;
            counter++;
        }
    }
    W_3_plus.array_size = counter;
    for(int omega = 0; omega < W_3_plus.array_size; omega++){
//        printf("W_3_plus[%d] = %d\n", omega, W_3_plus.array[omega]);
    }
    return W_3_plus;
}

//w_3_minus
struct W_j_object calculate_W_3_minus(int* F_omega, int size, int WH_max){
    int index = 0;
    int counter = 0;
    struct W_j_object W_3_minus;
//    printf("-WH_max + 4 = %d\n", -WH_max + 4);
    W_3_minus.array = (int*) malloc(i * sizeof(int));
    for(int omega = 0; omega < size; omega++){
        if(F_omega[omega] == (-WH_max) + 4){
            W_3_minus.array[index++] = omega;
            counter++;
        }
    }
    W_3_minus.array_size = counter;
    for(int omega = 0; omega < W_3_minus.array_size; omega++){
//        printf("W_3_minus[%d] = %d\n", omega, W_3_minus.array[omega]);
    }
    return W_3_minus;
}

int* matrix_multiplication(int *array, int **matrix) {
    int *result = (int *) calloc(i, sizeof(int));
    for(int row = 0; row < i; row++) {
        int temp_sum = 0;
        int sum = 0;
        for(int col = 0; col < i; col++){
            temp_sum = array[col] * matrix[row][col];
            sum = sum + temp_sum;
        }
        result[row] = sum;
    }
    return result;
}

int calculate_WH_max (int* walsh_spectrum, int size){
    int WH_max = abs(walsh_spectrum[1]);
    for(int index = 1; index < size; index++){
        if(WH_max < abs(walsh_spectrum[index])){
            WH_max = abs(walsh_spectrum[index]);
        }
    }
    return WH_max;
}

struct array_object calculate_union(int* array_1, int* array_2, int array_size_1, int array_size_2){
    struct array_object union_array;
    union_array.array = (int*) malloc((array_size_1 + array_size_2) * sizeof(int));
    union_array.array[0] = array_1[0];
    union_array.array_size = 0;
    int index_i;
    int index_j;
    int index_k;
    // copy the element of set A in set C
    for(index_i = 0; index_i < array_size_1; index_i++){
        for(index_j = 0; index_j < index_k; index_j++){
            if(union_array.array[index_j] == array_1[index_i]){
                break;
            }
        }
        //if not repesated then store value in set c
        if(index_j == index_k){
            union_array.array[index_k] = array_1[index_i];
            index_k++;
//            printf("index_k = %d\n", index_k);
        }
    }
    // copy element of set B in set C
    for(index_i = 0; index_i < array_size_2; index_i++){
        for(index_j = 0; index_j < index_k; index_j++){
            if(union_array.array[index_j] == array_2[index_i]){
                break;
            }
        }
        // if element is not repeated then store in set C
        if(index_j == index_k){
            union_array.array[index_k] = array_2[index_i];
            index_k++;
//            printf("index_k = %d\n", index_k);
        }
    }
    union_array.array_size = index_k ;
    return union_array;
};

//Complement f_j_i
int complement_f_j_i(int f_j_i){
    return f_j_i == 0 ? 1 : 0;
}

void print_Sbox(int **S_box, int row_number, int col_number){
    for(int col = 0; col < col_number; col++){
        printf("j = %d\n", col);
        for(int row = 0; row < row_number; row++){
            printf("D[%d][%d] = %d\n", row, col, S_box[row][col]);
        }
    }
}

int** calculate_v_i_j(int** S_box, int row_number, int col_number) {
    //print_Sbox(D,i,j);
    //CALCULATING v_i_j: CALCULATE INTEGER VALUE OF f_j_i AT EACH ROW i FROM THE BINARY VALUES FROM COLUMN j
    //FOR EXAMPLE: i = 0, j = 7, 00000111 = 7
    int offset = 0;

    for (int row = 0; row < row_number; row++) {
        offset = 1;
        int sum = 0;
        //COLUMN NUMBER STARTS FROM THE LSB: j = 0 IS LSB, j = n-1 IS THE MSB
        //ONLY CALCULATE THE SUM FOR EACH ROW WITH SET BITS
        for (int bit_position = 0; bit_position < n; bit_position++) {
            sum = D[row][bit_position] * (1 << (offset - 1));
            v_i_j[row][offset - 1] = sum;
            //printf("v_i_j[%d][%d] = %d\n", rows, offset-1, v_i_j[rows][offset-1]);
            offset++;
        }
    }


    //CALCULATING INTEGER VALUE OF v_i_j FOR EACH ROW OF D(j) when j = 1, 2, 3, ... n-1
    //FOR EXAMPLE, j = 2: (v[0][1] * 2^1) + (v[0][0] * 2^0)
    //             j = 7: (v[0][7] * 2^7) +
    //                    (v[0][6] * 2^6) +
    //                    (v[0][5] * 2^5) +
    //                    (v[0][4] * 2^4) +
    //                    (v[0][3] * 2^3) +
    //                    (v[0][2] * 2^2) +
    //                    (v[0][1] * 2^1) +
    //                    (v[0][0] * 2^0) +

    for (int col = 1; col < j; col++) {
        for (int rows = 0; rows < i; rows++) {
            v_i_j[rows][col] = v_i_j[rows][col] + v_i_j[rows][col - 1];
            //printf("v_i_j[%d][%d] = %d\n", rows, col, v_i_j[rows][col]);
        }
    }

    //CALCULATION OF #X_v(j) AND BIJECTIVE FITNESS FUNCTION, F(D(j)
    for (int col = 0; col < j; col++) {
        //printf("j = %d\n", col);
        for (int row = 0; row < i; row++) {
            v_j[row] = v_i_j[row][col];
            //printf("v_j[%d] = %d\n", row, v_j[row]);
        }
    }

    return v_i_j;
}