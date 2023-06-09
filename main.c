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

struct D_object{
    int* array;
    int array_size;
    int** v_i_j;
    int v_i_j_size;
    int* X_v;
    int X_v_size;
    int* F_D_j;
    int F_D_j_size;
};
int n = 8;
int i;
int j;
int** v_i_j;
int* v_j;

struct W_j_object W_plus_j_1;
struct W_j_object W_minus_j_1;
struct W_j_object W_plus_j_2;
struct W_j_object W_minus_j_2;
struct W_j_object W_plus_j_3;
struct W_j_object W_minus_j_3;

struct array_object W_j_1;
struct array_object W_plus_j_1_2;
struct array_object W_minus_j_1_2;
struct array_object W_plus_j_2_3;
struct array_object W_minus_j_2_3;
struct array_object W;

//FUNCTION SIGNATURES
int generateBinaryValue();
int** allocate_2D_matrix(int rows, int cols);
int** generate_Sbox_Matrix(int row_number, int col_number);
void free_2D_array(int** array_2D, int array_2D_length);
int* calculate_X_v(int* A, int n);
int calculate_F_D_j(int *x_k_j, int rows, int cols);
int** generate_WHM(int order);
int calculate_WH_max (int* walsh_spectrum, int size);
int** ptr_sbox_matrix;
int** D;
struct W_j_object calculate_W_plus_1(int* F_omaga, int size, int WH_max);
struct W_j_object calculate_W_minus_1(int* F_omega, int size, int WH_max);
struct W_j_object calculate_W_plus_2(int* F_omega, int size, int WH_max);
struct W_j_object calculate_W_minus_2(int* F_omega, int size, int WH_max);
struct W_j_object calculate_W_plus_3(int* F_omega, int size, int WH_max);
struct W_j_object calculate_W_minus_3(int* F_omega, int size, int WH_max);
int* matrix_multiplication(int *array, int **matrix);
struct array_object calculate_union(int* array_1, int* array_2, int array_size_1, int array_size_2);
int complement_f_j_i(int f_j_i);
void print_Sbox(int **S_box, int row_number, int col_number);
int** calculate_v_i_j(int** S_box, int row_number, int col_number);
int calculate_L_omega (int i_value, int omega);
struct D_object calculate_V(int** D_input_for_V, int col_number);
int is_satisfy_theorem_3_without_omega(int* f, int row, int* W_plus_j_1_2_input, int* W_minus_j_1_2_input);
int is_satisfy_theorem_3_with_omega(int* f, int row, int* W_plus_j_1_2_input, int* W_minus_j_1_2_input, int *W_input);
int** make_D_bijective(int** D, int j_index);
//int** make_D_bijective(int** D, int col_number);
int main() {
    i = 1 << n;
    printf("");
    j = n;

    // GENERATION OF A RANDOM S_BOX
    ptr_sbox_matrix = allocate_2D_matrix(i, j);
    ptr_sbox_matrix = generate_Sbox_Matrix(i, j);
    D = allocate_2D_matrix(i, j);
    int** bijective_D = allocate_2D_matrix(i, j);
    for (int col = 0; col < 1; col++) {
        for (int row = 0; row < i; row++) {
            //CONSTRUCT D(j) FROM FROM THE RANDOM SBOX WHERE j = 0, 1, 2, ..., n-1
            D[row][col] = ptr_sbox_matrix[row][j - col - 1];
            //printf("D[%d][%d] = %d\n", row, col, D[row][col]);
        }
    }

    bijective_D = make_D_bijective(D, 0);

    for (int col = 0; col < 1; col++) {
        for (int row = 0; row < i; row++) {
            printf("bijective_D[%d][%d] = %d\n", row, col, bijective_D[row][col]);
        }
    }

    return 0;
} //END OF MAIN FUNCTION

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

int** generate_Sbox_Matrix(int row_number, int col_number){
    // Seed the random number generator with current time
//    srand(time(0));

    // Define the size of the matrix
    int rows = row_number;
    int cols = col_number;

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
        for(int index_i = 0; index_i < k; index_i++){
            for(int index_j = 0; index_j < k; index_j++){
                hadamard[index_i+k][index_j] = hadamard[index_i][index_j];
                hadamard[index_i][index_j+k] = hadamard[index_i][index_j];
                hadamard[index_i+k][index_j+k] = -hadamard[index_i][index_j];
            }
        }
    }
    //Displaying final hadamard matrix
    for(int index_i = 0; index_i < n_WHM;index_i++){
        for(int index_j = 0; index_j < n_WHM; index_j++){
//            printf("%d ", hadamard[i][j]);
        }
//        printf("\n");
    }

    return hadamard;
}

struct W_j_object calculate_W_plus_1(int* F_omega, int size, int WH_max){
    struct W_j_object W_plus_1;
    int index = 0;
    int counter = 0;
    W_plus_1.array = (int*)malloc(i * sizeof(int));
    W_plus_1.array_size = 0;
//    printf("WH_max = %d\n", WH_max);
    for(int omega = 0; omega < size; omega++){
        if (F_omega[omega] == WH_max) {
            W_plus_1.array[index++] = omega;
            counter++;
        }
    }
    W_plus_1.array_size = counter;
    //printf("W_plus_1:\n");
    for(int omega = 0; omega < W_plus_1.array_size; omega++){
        //printf("W_plus_1[%d] = %d:\n", index, W_plus_1.array[omega]);
    }
    return W_plus_1;
}

//w_1_2_plus
struct W_j_object calculate_W_minus_1(int* F_omaga, int size, int WH_max){
    int index = 0;
    int counter = 0;
    struct W_j_object W_minus_1;
    W_minus_1.array = (int*)malloc(i * sizeof(int));
//    printf("-WH_max = %d\n", -WH_max);
    for(int omega = 0; omega < size; omega++){
        if (F_omaga[omega] == -WH_max) {
            W_minus_1.array[index++] = omega;
            counter++;
        }
    }
    W_minus_1.array_size = counter;

    for(int omega = 0; omega < W_minus_1.array_size; omega++){
        //printf("W_minus_1[%d] = %d:\n", omega, W_minus_1.array[omega]);
    }
    return W_minus_1;
}

//w_2_plus
struct W_j_object calculate_W_plus_2(int* F_omega, int size, int WH_max){
    int index = 0;
    int counter = 0;
    struct W_j_object W_plus_2;
    W_plus_2.array = (int*) malloc(i * sizeof(int));
//    printf("WH_max - 2 = %d\n", WH_max - 2);
    for(int omega = 0; omega < size; omega++){
        if(F_omega[omega] == WH_max - 2){
            W_plus_2.array[index++] = omega;
            counter++;
        }
    }
    W_plus_2.array_size = counter;

    for(int omega = 0; omega < W_plus_2.array_size; omega++){
        //printf("W_plus_2[%d] = %d:\n", omega, W_plus_2.array[omega]);
    }
    return  W_plus_2;
}

//w_2_minus
struct W_j_object calculate_W_minus_2(int* F_omega, int size, int WH_max){
    int index = 0;
    int counter = 0;
    struct W_j_object W_minus_2;
//    printf("-WH_max + 2 = %d\n", -WH_max + 2);
    W_minus_2.array = (int*)malloc(i * sizeof(int));
    for(int omega = 0; omega < size; omega++){
        if(F_omega[omega] == -WH_max + 2){
            W_minus_2.array[index++] = omega;
            counter++;
        }
    }
    W_minus_2.array_size = counter;
    for(int omega = 0; omega < W_minus_2.array_size; omega++){
        //printf("W_minus_2[%d] = %d:\n", omega, W_minus_2.array[omega]);
    }
    return W_minus_2;
}

//w_3_plus
struct W_j_object calculate_W_plus_3(int* F_omega, int size, int WH_max){
    int index = 0;
    int counter = 0;
    struct W_j_object W_plus_3;
//    printf("WH_max - 4 = %d\n", WH_max - 4);
    W_plus_3.array = (int*)malloc(i * sizeof(int));
    for(int omega = 0; omega < size; omega++) {
        if (F_omega[omega] == WH_max - 4) {
            W_plus_3.array[index++] = omega;
            counter++;
        }
    }
    W_plus_3.array_size = counter;
    for(int omega = 0; omega < W_plus_3.array_size; omega++){
        //printf("W_plus_3[%d] = %d\n", omega, W_plus_3.array[omega]);
    }
    return W_plus_3;
}

//w_3_minus
struct W_j_object calculate_W_minus_3(int* F_omega, int size, int WH_max){
    int index = 0;
    int counter = 0;
    struct W_j_object W_minus_3;
//    printf("-WH_max + 4 = %d\n", -WH_max + 4);
    W_minus_3.array = (int*) malloc(i * sizeof(int));
    for(int omega = 0; omega < size; omega++){
        if(F_omega[omega] == (-WH_max) + 4){
            W_minus_3.array[index++] = omega;
            counter++;
        }
    }
    W_minus_3.array_size = counter;
    for(int omega = 0; omega < W_minus_3.array_size; omega++){
        //printf("W_minus_3[%d] = %d\n", omega, W_minus_3.array[omega]);
    }
    return W_minus_3;
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

struct array_object calculate_union(int* array_1, int* array_2, int array_size_1, int array_size_2) {
    struct array_object union_array;
    union_array.array = (int *) malloc((array_size_1 + array_size_2) * sizeof(int));
    union_array.array[0] = array_1[0];
    union_array.array_size = 0;
    int index_i;
    int index_j;
    int index_k;

    if ((array_size_1 == 0) && (array_size_2 == 0)) {
        union_array.array_size = 0;
        return union_array;
    } else {
        // copy the element of set A in set C
        for (index_i = 0; index_i < array_size_1; index_i++) {
            for (index_j = 0; index_j < index_k; index_j++) {
                if (union_array.array[index_j] == array_1[index_i]) {
                    break;
                }
            }
            //if not repeated then store value in set c
            if (index_j == index_k) {
                union_array.array[index_k] = array_1[index_i];
                index_k++;
//            printf("index_k = %d\n", index_k);
            }
        }
        // copy element of set B in set C
        for (index_i = 0; index_i < array_size_2; index_i++) {
            for (index_j = 0; index_j < index_k; index_j++) {
                if (union_array.array[index_j] == array_2[index_i]) {
                    break;
                }
            }
            // if element is not repeated then store in set C
            if (index_j == index_k) {
                union_array.array[index_k] = array_2[index_i];
                index_k++;
//            printf("index_k = %d\n", index_k);
            }
        }
        union_array.array_size = index_k;
        return union_array;
    };
}
//Complement f_j_i
int complement_f_j_i(int f_j_i){
    return f_j_i == 0 ? 1 : 0;
}

void print_Sbox(int **S_box, int row_number, int col_number){
    for(int col = 0; col < col_number; col++){
//        printf("j = %d\n", col);
        for(int row = 0; row < row_number; row++){
//            printf("D[%d][%d] = %d\n", row, col, S_box[row][col]);
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

int calculate_L_omega (int i_value, int omega){
    int L_omega = 0;
    int product = 0;
    int x_bit;
    int omega_bit;

    //x =     00000000
    //omega = 00000001
    //L_omega = 0.1 ^ 0.0 ^ 0.0 ^ 0^0 ^ 0.0 ^ 0.0 ^ 0.0 ^ 0.0 = 0

    //x =     00000001
    //omega = 00000001
    //L_omega = 1.1 ^ 0.0 ^ 0.0 ^ 0^0 ^ 0.0 ^ 0.0 ^ 0.0 ^ 0.0 = 1

    //x =     00000010
    //omega = 00000001
    //L_omega = 0.1 ^ 1.0 ^ 0.0 ^ 0^0 ^ 0.0 ^ 0.0 ^ 0.0 ^ 0.0 = 1

    for(int bit_position = 0;bit_position < 8; bit_position++){
        x_bit = (i_value >> bit_position) & 0x01;
        omega_bit = (omega >> bit_position) & 0x01;
        //printf("x_bit = %d\n", x_bit);
        //printf("omega_bit = %d\n", omega_bit);
       if(bit_position == 0){
           product = x_bit  * omega_bit ;
           //printf("%d * %d = %d\n", x_bit, omega_bit, product);

           L_omega = product;
       }else{
           //printf("%d ^ %d = %d\n", L_omega, product,L_omega ^ product);
           product = x_bit  * omega_bit ;
           //printf("%d * %d = %d\n", x_bit, omega_bit, product);
           L_omega = L_omega ^ product;
       }

        //printf("L_omega = %d\n", L_omega);
        //printf("\n");
    }

    return L_omega;
}



int is_satisfy_theorem_3_without_omega(int* f, int row, int* W_plus_j_1_2_input, int* W_minus_j_1_2_input){
    int is_condition_1_satisfied = 0;
    int is_condition_2_satisfied = 0;
    int is_satisfy;
    int L_omega;
    //TO SATISFY THEOREM 3:
    //ASSUMPTION F(D(j-1)) = 0 and F(D(j)) > 2^n-1-j
    //CONDITION # 1
    //TO SATISFY THE CONDITION WE NEED TO KNOW omega, i, j
    //j IS THE COLUMN WE GET FROM D(j)
    //i IS THE ROW
    //omega IS THE SET OF W^+_j,1,2

    for(int index = 0; index < W_plus_j_1_2.array_size; index++){
        for(int omega = 0; omega < i; omega++){
            if(omega == W_plus_j_1_2.array[index]){
                L_omega = calculate_L_omega(row, omega);
                if(f[row] == L_omega){
                    is_condition_1_satisfied = 1;
                    //printf("is_condition_1_satisfied = %d\n", is_condition_1_satisfied);
                };
            }
        }
    }
    //CONDITION # 2
    //TO SATISFY THE CONDITION WE NEED TO KNOW omega, i, j
    //j IS THE COLUMN WE GET FROM D(j)
    //i IS THE ROW
    //omega IS THE SET OF W^-_j,1,2
    for(int index = 0; index < W_minus_j_1_2.array_size; index++){
        //printf("W_minus_j_1_2.array[%d] = %d\n", index, W_minus_j_1_2.array[index]);
        for(int omega = 0; omega < i; omega++){
            if(omega == W_minus_j_1_2.array[index]){
//                printf("omega = %d\n", omega);
                L_omega = calculate_L_omega(row, omega);
//                printf("L_omega = %d\n", L_omega);
//                printf("f[%d] = %d\n", row, f[row]);
                if(f[row] != L_omega){
                    is_condition_2_satisfied = 1;
                    //printf("is_condition_1_satisfied = %d\n", is_condition_1_satisfied);
                };
            }
        }
    }

    if( is_condition_1_satisfied && is_condition_2_satisfied){
        is_satisfy = 1;
    }else {
        is_satisfy = 0;
    }

    return is_satisfy;
}

int is_satisfy_theorem_3_with_omega(int* f, int row, int* W_plus_j_1_2_input, int* W_minus_j_1_2_input, int *W_input){
    int is_condition_1_satisfied = 0;
    int is_condition_2_satisfied = 0;
    int is_satisfy;
    int L_omega;
    //TO SATISFY THEOREM 3:
    //ASSUMPTION F(D(j-1)) = 0 and F(D(j)) > 2^n-1-j
    //CONDITION # 1
    //TO SATISFY THE CONDITION WE NEED TO KNOW omega, i, j
    //j IS THE COLUMN WE GET FROM D(j)
    //i IS THE ROW
    //omega IS THE SET OF W^+_j,1,2

    for(int index = 0; index < W_plus_j_1_2.array_size; index++){
        for(int omega = 0; omega < i; omega++){
            if(omega == W_plus_j_1_2.array[index]){
                L_omega = calculate_L_omega(row, omega);
                if(f[row] == L_omega){
                    is_condition_1_satisfied = 1;
                    //printf("is_condition_1_satisfied = %d\n", is_condition_1_satisfied);
                };
            }
        }
    }
    //CONDITION # 2
    //TO SATISFY THE CONDITION WE NEED TO KNOW omega, i, j
    //j IS THE COLUMN WE GET FROM D(j)
    //i IS THE ROW
    //omega IS THE SET OF W^-_j,1,2
    for(int index = 0; index < W_minus_j_1_2.array_size; index++){
        //printf("W_minus_j_1_2.array[%d] = %d\n", index, W_minus_j_1_2.array[index]);
        for(int omega = 0; omega < i; omega++){
            if(omega == W_minus_j_1_2.array[index]){
//                printf("omega = %d\n", omega);
                L_omega = calculate_L_omega(row, omega);
//                printf("L_omega = %d\n", L_omega);
//                printf("f[%d] = %d\n", row, f[row]);
                if(f[row] != L_omega){
                    is_condition_2_satisfied = 1;
                    //printf("is_condition_1_satisfied = %d\n", is_condition_1_satisfied);
                };
            }
        }
    }
    if( is_condition_1_satisfied && is_condition_2_satisfied){
        is_satisfy = 1;
    }else {
        is_satisfy = 0;
    }
    return is_satisfy;
}

//CALCULATION OF V
struct D_object calculate_V(int** D_input_for_V, int col_number){
    struct D_object V;

    V.array = (int*) malloc(i * sizeof(int));
    V.array_size = 0;
    V.v_i_j = allocate_2D_matrix(i,j);
    V.v_i_j_size = i;

    V.X_v = (int*)malloc(i * sizeof(int));
    V.X_v_size = 0;

    V.F_D_j = (int*) malloc(i * sizeof(int));
    V.F_D_j_size = 0;

    //CALCULATION OF v_i_j USING FUNCTION
    int counter = 0;
    counter++;
//        printf("counter = %d\n", counter);
    V.v_i_j = calculate_v_i_j(D_input_for_V, i, col_number);
    // printf("j = %d\n", col);
    int sum = 0;
    for (int row = 0; row < i; row++) {
        v_j[row] =  V.v_i_j[row][col_number];
        //printf("v_j[%d] = %d\n", row, v_j[row]);
    }

    //CALCULATION OF x_v and F(D(j)
    //CALCULATING OF #x_v(j)
    //X_v CALCULATES REPETITION OF THE VALUE v IN SUB-MATRIX D(j) FOR j = 0, 1, 2, ... , n-1
    V.X_v = calculate_X_v(v_j, i);

    //CALCULATING BIJECTIVE FITNESS FUNCTION, F(D(j))
    V.F_D_j[col_number] = calculate_F_D_j(V.X_v, j, col_number);
//    printf("F_D_j[%d] = %d\n", col_number, V.F_D_j[col_number]);
    //printf("\n");

    //STEP # 1: FOR D(j) construct a set V = {v = v| v ∈ {0, 1, 2, 3, ..., 2^(j+1)-1} if #{X_v(j) > 2^n-1-j}
    int v = (1 << (col_number + 1)) - 1;
    //printf("v = %d\n", v);
    int index_V = 0;
    for (int index = 0; index <= v; index++) {
        //printf("#{X_v[%d]} = %d\n", index, X_v[index]);
        //printf("2^n-1-%d  = %d\n", col, 1 << (n - 1 - col));
        if (V.X_v[index] > (1 << (n - 1 - col_number))) {
            V.array[index_V++] = index;
            //printf("%d > %d\n", X_v[index], 1 << (n - 1 - col));
            //printf("v = %d\n", index);
            V.array_size = index_V;
        }
        //printf("\n");
        //sum = X_v[index] + sum;
    }
    return V;
}

//struct D_object calculate_V(int* X_v, int col_number){
//
//    struct D_object V;
//    //CALCULATION OF v_i_j USING FUNCTION
//    int counter = 0;
//    counter++;
//        printf("counter = %d\n", counter);
//    v_i_j = calculate_v_i_j(D, i, col_number);
//    // printf("j = %d\n", col);
//    int sum = 0;
//    for (int row = 0; row < i; row++) {
//        v_j[row] = v_i_j[row][col_number];
//        //printf("v_j[%d] = %d\n", row, v_j[row]);
//    }
//
//    //CALCULATION OF x_v and F(D(j)
//    //CALCULATING OF #x_v(j)
//    //X_v CALCULATES REPETITION OF THE VALUE v IN SUB-MATRIX D(j) FOR j = 0, 1, 2, ... , n-1
//    X_v = calculate_X_v(v_j, i);
//
//    //CALCULATING BIJECTIVE FITNESS FUNCTION, F(D(j))
//    F_D_j[col_number] = calculate_F_D_j(X_v, j, col_number);
//    //printf("F_D_j[%d] = %d\n", col, F_D_j[col]);
//    //printf("\n");
//
//
//    //CALCULATION OF V
//    struct array_object V;
//    V.array = (int*)malloc(i * sizeof(int));
//    V.array_size = 0;
//
//    //STEP # 1: FOR D(j) construct a set V = {v = v| v ∈ {0, 1, 2, 3, ..., 2^(j+1)-1} if #{X_v(j) > 2^n-1-j}
//    int v = (1 << (col_number + 1)) - 1;
//    //printf("v = %d\n", v);
//    int index_V = 0;
//    for (int index = 0; index <= v; index++) {
//        //printf("#{X_v[%d]} = %d\n", index, X_v[index]);
//        //printf("2^n-1-%d  = %d\n", col, 1 << (n - 1 - col));
//        if (X_v[index] > (1 << (n - 1 - col_number))) {
//            V.array[index_V++] = index;
//            //printf("%d > %d\n", X_v[index], 1 << (n - 1 - col));
//            //printf("v = %d\n", index);
//            V.array_size = index_V;
//        }
//        //printf("\n");
//        //sum = X_v[index] + sum;
//    }
//    return V;
//}


int** make_D_bijective(int** D_input, int j_index){
    int **f_j = allocate_2D_matrix(i, j);
    int *S_j_omega = (int *) malloc(i * sizeof(int));
    int offset = 0;
    v_i_j = allocate_2D_matrix(i, j);
    v_j = (int *) malloc(j * sizeof(int));
    int *X_v = (int *) malloc(i * sizeof(int));
    int *updated_X_v = (int *) malloc(i * sizeof(int));
    int *F_D_j = (int *) malloc(j * sizeof(int));
    int *f_j_i = (int *) malloc(i * sizeof(int));
    int *F_omega = (int *) malloc(16 * sizeof(int));
    int WH_max;
    int *updated_D = (int *) malloc(j * sizeof(int));
    int L_omega;
    struct D_object V;
    int is_bijective_flag;

    printf("Calling bijective function\n");

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
    printf("col_number_inside function = %d\n", j_index);
    //ASSIGNING f_j where f_0 is MSB and f_(j-1) is MSB as represented by the array
    //WE HAVE TO MAKE D BIJECTIVE FOR EACH col AND MOVE TO THE NEXT col.
    for (int col = 0; col <= j_index; col++) {
        printf("j = %d\n", col);
        for (int row = 0; row < i; row++) {
            f_j_i[row] = D_input[row][col];
//            printf("Input_D[%d][%d] = %d\n", row, col, D_input[row][col]);
//            printf("f_j[%d][%d] = %d\n", row, col, f_j[row][col]);
            //printf("f_j_i[%d] = %d\n", row, f_j_i[row]);
        }

        //CALCULATION OF v_i_j USING FUNCTION

        int counter = 0;
        counter++;

        //CALCULATION OF V
        V = calculate_V(D_input, j_index);


        //        printf("counter = %d\n", counter);
//        v_i_j = calculate_v_i_j(D, i, col);
        // printf("j = %d\n", col);
        int sum = 0;
        for (int row = 0; row < i; row++) {
            v_j[row] = V.v_i_j[row][col];
            //printf("v_j[%d] = %d\n", row, v_j[row]);
        }

        //CALCULATION OF x_v and F(D(j)
        //CALCULATING OF #x_v(j)
        //X_v CALCULATES REPETITION OF THE VALUE v IN SUB-MATRIX D(j) FOR j = 0, 1, 2, ... , n-1
        X_v = calculate_X_v(v_j, i);

        //CALCULATING BIJECTIVE FITNESS FUNCTION, F(D(j))
//        F_D_j[col] = calculate_F_D_j(X_v, j, col);
        printf("F_D_j[%d] = %d\n", col, V.F_D_j[col]);

        //printf("\n");
        //Calculate F_D(j-1) and F_F(j)

        //printf("index_V = %d\n");
        for (int index = 0; index < V.array_size; index++) {
            printf("V[%d] = %d\n", index, V.array[index]);
        }

        //STEP # 2: IF V IS ∅, OUTPUT D(j), OTHERWISE GO TO STEP # 3
        if (V.array_size == 0) {
            printf("D(%d) is bijective\n", col);
            for (int row = 0; row < i; row++) {
                //printf("D[%d][%d] = %d\n", row, col, D[row][col]);
            }
            continue; //IF D IS BIJECTIVE FOR j THEN TO GOT
        } else {
            printf("D(%d) is NOT bijective\n", col);
            printf("Going to Step 3\n", col);
            //STEP # 3. According to the Boolean function f_j in D(j), calculate the sets W^plus_j,1, W^minus_j,1,W^plus_j,2 AND W^minus_j,2.
            //APPLYING WALSH HADAMARD MATRIX ON f_j TO CALCULATE WALSH HADAMARD SPECTRUM OF f_j
            S_j_omega = matrix_multiplication(f_j_i, WH);

            // printf("j = %d\n", col);
            //printf("Walsh-Hadamard transform for f Boolean function:\n");
            //PRINTING WALSH-HADAMARD TRANSFORM
            for (int omega = 0; omega < i; omega++) {
                //printf("%d ", S_j_omega[omega]);
                //printf("f[%d] = %d\n", omega,  f[omega]);
                //printf("S_j_omega[%d] = %d\n", omega, S_j_omega[omega]);
                //printf("f[%d] = %d\n", omega,  F_omega[omega]);
            }

            //CALCULATING MAXIMUM VALUE FOR WALSH-HADAMARD TRANSFORM
            WH_max = calculate_WH_max(S_j_omega, i);
            printf("WH_max = %d\n", WH_max);

            //CALCULATION OF W_plus_j_1, W_minus_j_1, W_plus, W_plus_j_2, W_minus_j_2, W_plus_j_3 and W_minus_j_3
            W_plus_j_1 = calculate_W_plus_1(S_j_omega, i, WH_max);
            W_minus_j_1 = calculate_W_minus_1(S_j_omega, i, WH_max);
            W_plus_j_2 = calculate_W_plus_2(S_j_omega, i, WH_max);
            W_minus_j_2 = calculate_W_minus_2(S_j_omega, i, WH_max);
            W_plus_j_3 = calculate_W_plus_3(S_j_omega, i, WH_max);
            W_minus_j_3 = calculate_W_minus_3(S_j_omega, i, WH_max);

            if (W_plus_j_1.array_size == 0) {
//                printf("W_plus_j_1 is empty\n");
            } else {
//                printf("W_plus_j_1.array_size = %d\n", W_plus_j_1.array_size);
                for (int index = 0; index < W_plus_j_1.array_size; index++) {
//                    printf("W_plus_j_1[%d] = %d\n", index, W_plus_j_1.array[index]);
                }
            }

            if (W_minus_j_1.array_size == 0) {
//                printf("W_minus_j_1 is empty\n");
            } else {
//                printf("W_minus_j_1.array_size = %d\n", W_minus_j_1.array_size);
                for (int index = 0; index < W_minus_j_1.array_size; index++) {
//                    printf("W_minus_j_1[%d] = %d\n", index, W_minus_j_1.array[index]);
                }
            }

            if (W_plus_j_2.array_size == 0) {
//                printf("W_plus_j_2 is empty\n");
            } else {
//                printf("W_plus_j_2.array_size = %d\n", W_plus_j_2.array_size);
                for (int index = 0; index < W_plus_j_2.array_size; index++) {
//                    printf("W_j_2_plus[%d] = %d\n", index, W_plus_j_2.array[index]);
                }
            }

            if (W_minus_j_2.array_size == 0) {
//                printf("W_minus_j_2 is empty\n");
            } else {
//                printf("W_minus_j_2.array_size = %d\n", W_minus_j_2.array_size);
                for (int index = 0; index < W_minus_j_2.array_size; index++) {
//                    printf("W_j_2_minus[%d] = %d\n", index, W_minus_j_2.array[index]);
                }
            }

            if (W_plus_j_3.array_size == 0) {
//                printf("W_plus_j_3 is empty\n");
            } else {
//                printf("W_plus_j_3.array_size = %d\n", W_plus_j_3.array_size);
                for (int index = 0; index < W_plus_j_3.array_size; index++) {
//                    printf("W_plus_j_3_[%d] = %d\n", index, W_plus_j_3.array[index]);
                }
            }

            if (W_minus_j_3.array_size == 0) {
//                printf("W_minus_j_3 is empty\n");
            } else {
//                printf("W_minus_j_3.array_size = %d\n", W_minus_j_3.array_size);
                for (int index = 0; index < W_minus_j_3.array_size; index++) {
//                    printf("W_j_3_minus[%d] = %d\n", index, W_minus_j_3.array[index]);
                }
            }

            //CALCULATING W^+_j,1,2
            W_plus_j_1_2 = calculate_union(W_plus_j_1.array, W_plus_j_2.array, W_plus_j_1.array_size,
                                           W_plus_j_2.array_size);


            //CALCULATING W^+_j,1,2
            W_minus_j_1_2 = calculate_union(W_minus_j_1.array, W_minus_j_2.array, W_minus_j_1.array_size,
                                            W_minus_j_2.array_size);

            //printf("W_plus_j_1_2.array_size = %d\n", W_plus_j_1_2.array_size);
            if (W_plus_j_1_2.array_size == 0) {
//                printf("W_plus_j_1_2 is empty\n");
            } else {
                for (int W_plus_j_1_2_index = 0;
                     W_plus_j_1_2_index < W_plus_j_1_2.array_size; W_plus_j_1_2_index++) {
//                    printf("W_plus_j_1_2[%d] = %d\n", W_plus_j_1_2_index, W_plus_j_1_2.array[W_plus_j_1_2_index]);
                }
            }

            if (W_minus_j_1_2.array_size == 0) {
//                printf("W_minus_j_1_2 is empty\n");
            } else {
                for (int W_minus_j_1_2_index = 0;
                     W_minus_j_1_2_index < W_minus_j_1_2.array_size; W_minus_j_1_2_index++) {
//                    printf("W_minus_j_1_2[%d] = %d\n", W_minus_j_1_2_index, W_minus_j_1_2.array[W_minus_j_1_2_index]);
                }
            }

            W = calculate_union(W_plus_j_1.array, W_minus_j_1.array, W_plus_j_1.array_size, W_minus_j_1.array_size);

            printf("\n");
            int a;
            int is_theorem_3_satisfied_without_omega;
            int is_theorem_3_satisfied_with_omega;


            // printf("index_v = %d\n", index_v);
            //STEP 4: Set m = 1
            for (int m = 0; m < V.array_size; m++) {
                //STEP 5: Set a=V(m),here V(m)is the m-th element in V int a;
                printf("m = %d\n", m);
                a = V.array[m];
//                printf("a = %d\n", a);
                int count_X_v_1 = 0;
                //STEP 6: CHECK FOR v_i_j = a and v_i_j satisfies Theorem 3, complement the Boolean function value f_j_i, update and return to Step 1;
                // Otherwise, go to Step 7.
                int match_flag = 0;
                for (int row = 0; row < i; row++) {
//                    printf("v_j_i[%d][%d] = %d\n", row, col, V.v_i_j[row][col]);
                    if ((a == v_j[row]) && (match_flag == 0)) { //Checking if we can get v_i_j = a
                        match_flag = 1;
//                        printf("rows_with_1 = %d\n", row);
                        //CHECK IF v_i_j SATISFY THEOREM 3 WHERE i = row

                        is_theorem_3_satisfied_without_omega = is_satisfy_theorem_3_without_omega(f_j_i, row,
                                                                                                  W_plus_j_1_2.array,
                                                                                                  W_minus_j_1_2.array);

//                        printf("is_theorem_3_satisfied_without_omega = %d\n", is_theorem_3_satisfied_without_omega);
                        if (is_theorem_3_satisfied_without_omega == 1) {
//                            printf("Theorem 3 is satisfied without omega\n");
//                            printf("Complement f_i_j for row = %d\n", row);
//                            printf("Before complement:f[%d] = %d\n", row, f_j_i[row]);

                            //COMPLEMENT f_i_j
                            complement_f_j_i(f_j_i[row]);

                            //UPDATE D(j)
                            D[row][col] = f_j_i[row];
                            D = make_D_bijective(D,col);
//                            printf("After complement:f[%d] = %d\n", row, f[row]);
                            break;
                        }
                    } else if ((a == v_i_j[row][col]) && (match_flag == 1)) {
                        W_j_1 = calculate_union(W_plus_j_1.array, W_minus_j_1.array, W_plus_j_1.array_size,
                                                W_minus_j_1.array_size);

                        is_theorem_3_satisfied_with_omega = is_satisfy_theorem_3_with_omega(f_j_i, row, W_plus_j_1_2.array,
                                                                                           W_minus_j_1_2.array,
                                                                                           W.array);;
//                        printf("is_theorem_3_satisfied_wit_omega = %d\n", is_theorem_3_satisfied_wit_omega);
                        if (is_theorem_3_satisfied_with_omega == 1) {
//                            printf("Theorem 3 is satisfied with omega\n");
//                            printf("Complement f_i_j for row = %d\n", row);
//                            printf("Before complement:f[%d] = %d\n", row, f[row]);

                            //COMPLEMENT f_i_j
                            complement_f_j_i(f_j_i[row]);

                            //UPDATE D(j)
                            D[row][col] = f_j_i[row];
                            D = make_D_bijective(D, col);
//                            printf("After complement:f[%d] = %d\n", row, f[row]);
                        }
                        //printf("\n");
                        //    printf("W_j_1.array_size = %d\n", W_j_1.array_size);
                        //    printf("W_j_1:\n");


                        if (col > 0) {
//                            printf("j = %d\n", col);
//                            printf("F(D(%d)) = %d\n", col-1, F_D_j[col-1]);
//                            printf("F(D(%d)) = %d\n", col, F_D_j[col]);
                        }
                        count_X_v_1++;
                        //printf("f_j[%d][%d] = %d\n", row, col, f_j[row][col]);
                        //printf("count_X_v_1 = %d\n", count_X_v_1);
                        if (count_X_v_1 > (1 << (n - 1 - col))) {
                            D[row][col] = complement_f_j_i(f_j[row][col]);
                        }
                    }else{
//                        printf("Theorem 3 is NOT satisfied without omega and moving to Step 8\n");
                    }


//                        if (col == 0) {
//                            printf("j = %d\n", col);
//                            printf("F(D(%d)) = %d\n", col, F_D_j[col]);
//                        }
                    //CALCULATING W^+_j,1,2
                    // W^-_j,1,2, W^+_j,2,3, W^-_j,2,3
//                        //Calling W_1_plus
//                        //printf("\n");
//                        //    printf("W_j_1_minus.array_size = %d\n", W_j_1_minus.array_size);
//                        //    printf("W_j_1_minus:\n");
//                        for (int index = 0; index < W_minus_j_1.array_size; index++) {
//                            //        printf("W_j_1_minus[%d] = %d\n", index, W_j_1_minus.array[index]);
//                        }
//

                    //printf("D[%d][%d] = %d\n", row, col, D[row][col]);
                    //printf("\n");

                }//END OF FOR LOOP of row CHECKING WHICH ROW IS EQUAL TO
                if(m <= V.array_size){
                    printf("Breaking loop for m\n");
                    break;
                }
            }//END OF FOR LOOP OF m

            printf("After breaking m loop\n");
        }
//        printf("V[%d] = %d\n", 0, V.array[0]);
    }
    //printf("\n");
    //END OF FOR LOOP COL OF j
    //printf("Updated S-box\n");
    //print_Sbox(D, i, j);

    for(int row = 0; row < i; row++){
        if(V.array[0] == v_j[row]){
            f_j_i[row] = complement_f_j_i(f_j_i[row]);
            D[row][j_index] = f_j_i[row];
            printf("F_D_j = %d\n", V.F_D_j[0]);
            if(V.F_D_j[0] != 0){
                is_bijective_flag = 0;
                D = make_D_bijective(D, j_index);
            }else{
                is_bijective_flag = 1;
            }
        }
    }

    return D;

}
