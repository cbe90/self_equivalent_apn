/* Compile with: g++ -O2 main_det8_class1.c -o det8_class1
Authors: C. Beierle, G. Leander  -- Feb 2020

This program constructs and returns all 8-bit APN permutations F s.t. AF = FB with A, B corresponding to class 1, i.e.,
A = Comp(X^8+X^7+X^6+X^4+X^2+X+1)
B = Comp(X^8+X^5+X^4+X^3+1),
up to affine equivalence. It may return multiple representatives in one affine-equivalence class. 
*/

#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <stdio.h>
#include <time.h>
#include "commuting_matrices_A.h"
#include "commuting_matrices_B.h"

// dimension of the APN permutation to find
#define  N 8

time_t start;
time_t runtime;

// the matrix A
int A[(1<<N)] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100, 102, 104, 106, 108, 110, 112, 114, 116, 118, 120, 122, 124, 126, 128, 130, 132, 134, 136, 138, 140, 142, 144, 146, 148, 150, 152, 154, 156, 158, 160, 162, 164, 166, 168, 170, 172, 174, 176, 178, 180, 182, 184, 186, 188, 190, 192, 194, 196, 198, 200, 202, 204, 206, 208, 210, 212, 214, 216, 218, 220, 222, 224, 226, 228, 230, 232, 234, 236, 238, 240, 242, 244, 246, 248, 250, 252, 254, 215, 213, 211, 209, 223, 221, 219, 217, 199, 197, 195, 193, 207, 205, 203, 201, 247, 245, 243, 241, 255, 253, 251, 249, 231, 229, 227, 225, 239, 237, 235, 233, 151, 149, 147, 145, 159, 157, 155, 153, 135, 133, 131, 129, 143, 141, 139, 137, 183, 181, 179, 177, 191, 189, 187, 185, 167, 165, 163, 161, 175, 173, 171, 169, 87, 85, 83, 81, 95, 93, 91, 89, 71, 69, 67, 65, 79, 77, 75, 73, 119, 117, 115, 113, 127, 125, 123, 121, 103, 101, 99, 97, 111, 109, 107, 105, 23, 21, 19, 17, 31, 29, 27, 25, 7, 5, 3, 1, 15, 13, 11, 9, 55, 53, 51, 49, 63, 61, 59, 57, 39, 37, 35, 33, 47, 45, 43, 41};

// the matrix B
int B[(1<<N)] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100, 102, 104, 106, 108, 110, 112, 114, 116, 118, 120, 122, 124, 126, 128, 130, 132, 134, 136, 138, 140, 142, 144, 146, 148, 150, 152, 154, 156, 158, 160, 162, 164, 166, 168, 170, 172, 174, 176, 178, 180, 182, 184, 186, 188, 190, 192, 194, 196, 198, 200, 202, 204, 206, 208, 210, 212, 214, 216, 218, 220, 222, 224, 226, 228, 230, 232, 234, 236, 238, 240, 242, 244, 246, 248, 250, 252, 254, 57, 59, 61, 63, 49, 51, 53, 55, 41, 43, 45, 47, 33, 35, 37, 39, 25, 27, 29, 31, 17, 19, 21, 23, 9, 11, 13, 15, 1, 3, 5, 7, 121, 123, 125, 127, 113, 115, 117, 119, 105, 107, 109, 111, 97, 99, 101, 103, 89, 91, 93, 95, 81, 83, 85, 87, 73, 75, 77, 79, 65, 67, 69, 71, 185, 187, 189, 191, 177, 179, 181, 183, 169, 171, 173, 175, 161, 163, 165, 167, 153, 155, 157, 159, 145, 147, 149, 151, 137, 139, 141, 143, 129, 131, 133, 135, 249, 251, 253, 255, 241, 243, 245, 247, 233, 235, 237, 239, 225, 227, 229, 231, 217, 219, 221, 223, 209, 211, 213, 215, 201, 203, 205, 207, 193, 195, 197, 199};

// inverse of B
int B_inv[(1<<N)] = {0, 156, 1, 157, 2, 158, 3, 159, 4, 152, 5, 153, 6, 154, 7, 155, 8, 148, 9, 149, 10, 150, 11, 151, 12, 144, 13, 145, 14, 146, 15, 147, 16, 140, 17, 141, 18, 142, 19, 143, 20, 136, 21, 137, 22, 138, 23, 139, 24, 132, 25, 133, 26, 134, 27, 135, 28, 128, 29, 129, 30, 130, 31, 131, 32, 188, 33, 189, 34, 190, 35, 191, 36, 184, 37, 185, 38, 186, 39, 187, 40, 180, 41, 181, 42, 182, 43, 183, 44, 176, 45, 177, 46, 178, 47, 179, 48, 172, 49, 173, 50, 174, 51, 175, 52, 168, 53, 169, 54, 170, 55, 171, 56, 164, 57, 165, 58, 166, 59, 167, 60, 160, 61, 161, 62, 162, 63, 163, 64, 220, 65, 221, 66, 222, 67, 223, 68, 216, 69, 217, 70, 218, 71, 219, 72, 212, 73, 213, 74, 214, 75, 215, 76, 208, 77, 209, 78, 210, 79, 211, 80, 204, 81, 205, 82, 206, 83, 207, 84, 200, 85, 201, 86, 202, 87, 203, 88, 196, 89, 197, 90, 198, 91, 199, 92, 192, 93, 193, 94, 194, 95, 195, 96, 252, 97, 253, 98, 254, 99, 255, 100, 248, 101, 249, 102, 250, 103, 251, 104, 244, 105, 245, 106, 246, 107, 247, 108, 240, 109, 241, 110, 242, 111, 243, 112, 236, 113, 237, 114, 238, 115, 239, 116, 232, 117, 233, 118, 234, 119, 235, 120, 228, 121, 229, 122, 230, 123, 231, 124, 224, 125, 225, 126, 226, 127, 227};

// an ordering of the indices mod 2^n. It orders cycle-wise with respect to the cycles of B
int lex[(1<<N)] = {0,1, 2, 4, 8, 16, 32, 64, 128, 57, 114, 228, 241, 219, 143, 39, 78, 156,3, 6, 12, 24, 48, 96, 192, 185, 75, 150, 21, 42, 84, 168, 105, 210, 157,5, 10, 20, 40, 80, 160, 121, 242, 221, 131, 63, 126, 252, 193, 187, 79, 158,7, 14, 28, 56, 112, 224, 249, 203, 175, 103, 206, 165, 115, 230, 245, 211, 159,9, 18, 36, 72, 144, 25, 50, 100, 200, 169, 107, 214, 149, 19, 38, 76, 152,11, 22, 44, 88, 176, 89, 178, 93, 186, 77, 154, 13, 26, 52, 104, 208, 153,15, 30, 60, 120, 240, 217, 139, 47, 94, 188, 65, 130, 61, 122, 244, 209, 155,17, 34, 68, 136, 41, 82, 164, 113, 226, 253, 195, 191, 71, 142, 37, 74, 148,23, 46, 92, 184, 73, 146, 29, 58, 116, 232, 233, 235, 239, 231, 247, 215, 151,27, 54, 108, 216, 137, 43, 86, 172, 97, 194, 189, 67, 134, 53, 106, 212, 145,31, 62, 124, 248, 201, 171, 111, 222, 133, 51, 102, 204, 161, 123, 246, 213, 147,33, 66, 132, 49, 98, 196, 177, 91, 182, 85, 170, 109, 218, 141, 35, 70, 140,45, 90, 180, 81, 162, 125, 250, 205, 163, 127, 254, 197, 179, 95, 190, 69, 138,55, 110, 220, 129, 59, 118, 236, 225, 251, 207, 167, 119, 238, 229, 243, 223, 135,83, 166, 117, 234, 237, 227, 255, 199, 183, 87, 174, 101, 202, 173, 99, 198, 181};

int pos[15] = {1,3,5,7,9,11,15,17,23,27,31,33,45,55,83}; // one element of each cycle within the cycle structure of B (without 0). The next free position x will be selected from this list of positions (long cycles indexed first)

int repr[(1<<N)] = {0, 1, 1, 3, 1, 5, 3, 7, 1, 9, 5, 11, 3, 11, 7, 15, 1, 17, 9, 19, 5, 21, 11, 7, 3, 19, 11, 27, 7, 7, 15, 9, 1, 27, 17, 35, 9, 5, 19, 9, 5, 5, 21, 43, 11, 45, 7, 47, 3, 35, 19, 1, 11, 43, 27, 15, 7, 9, 7, 15, 15, 47, 9, 47, 1, 19, 27, 21, 17, 43, 35, 19, 9, 35, 5, 27, 19, 17, 9, 1, 5, 43, 5, 21, 21, 3, 43, 3, 11, 17, 45, 45, 7, 47, 47, 35, 3, 21, 35, 45, 19, 21, 1, 27, 11, 27, 43, 15, 27, 45, 15, 15, 7, 19, 9, 27, 7, 3, 15, 11, 15, 1, 47, 15, 9, 35, 47, 43, 1, 35, 19, 3, 27, 45, 21, 21, 17, 47, 43, 3, 35, 17, 19, 43, 9, 47, 35, 43, 5, 45, 27, 5, 19, 7, 17, 21, 9, 11, 1, 5, 5, 45, 43, 47, 5, 11, 21, 7, 21, 45, 3, 35, 43, 17, 3, 47, 11, 17, 17, 17, 45, 17, 45, 11, 7, 11, 47, 45, 47, 17, 35, 45, 3, 3, 21, 43, 35, 47, 45, 17, 19, 43, 21, 19, 1, 17, 27, 35, 11, 3, 27, 19, 43, 35, 15, 1, 27, 21, 45, 21, 15, 45, 15, 27, 7, 21, 19, 17, 9, 7, 27, 19, 7, 5, 3, 1, 15, 11, 11, 9, 15, 43, 1, 35, 47, 47, 15, 9, 9, 5, 35, 27, 47, 45, 43, 5}; // the cycle representatives of each element of A from 0-255.

// numbers of even Hamming weight up to 255 (needed for more efficient APN check)
int evens[127] = {3, 5, 6, 9, 10, 12, 15, 17, 18, 20, 23, 24, 27, 29, 30, 33, 34, 36, 39, 40, 43, 45, 46, 48, 51, 53, 54, 57, 58, 60, 63, 65, 66, 68, 71, 72, 75, 77, 78, 80, 83, 85, 86, 89, 90, 92, 95, 96, 99, 101, 102, 105, 106, 108, 111, 113, 114, 116, 119, 120, 123, 125, 126, 129, 130, 132, 135, 136, 139, 141, 142, 144, 147, 149, 150, 153, 154, 156, 159, 160, 163, 165, 166, 169, 170, 172, 175, 177, 178, 180, 183, 184, 187, 189, 190, 192, 195, 197, 198, 201, 202, 204, 207, 209, 210, 212, 215, 216, 219, 221, 222, 225, 226, 228, 231, 232, 235, 237, 238, 240, 243, 245, 246, 249, 250, 252, 255};

int solutions;
long long iterations;
int n_skipped;
int skipped[128];

int sbox[(1<<N)]; // it stores the (partial) S-box that we are going to construct (undefined values represented by -1)
int sbox_DDT[(1<<N)][(1<<N)]; // it stores the (partial) DDT of sbox

void printArray(int A[1<<N]) {
	for (int i=0;i<(1<<N);i++)
		printf("%2x ",i);
	printf("\n");
	for (int i=0;i<(1<<N);i++)
		if (A[i]!=-1) printf("%2x ",A[i]);
			else printf("-- ");
	printf("\n\n");
}

void printSkippedArray(int A[128]) {
	for (int i=0;i<32;i++)
		printf("%d ",A[i]);
	printf("\n\n");
}

void printArray2(int A[1<<N]) {
	printf("{0x%02x",A[0]);
	for (int i=1;i<(1<<N);i++)
		printf(",0x%02x",A[i]);
	printf("};\n");
}

void fprintArray2(FILE *fp, int A[1<<N]) {
	fprintf(fp, "{0x%02x",A[0]);
	for (int i=1;i<(1<<N);i++)
		fprintf(fp, ",0x%02x",A[i]);
	fprintf(fp, "};\n");
}

void printDDT(int A[1<<N][1<<N]) {
    for (int i=0; i<(1<<N); i++) {
        for (int j=0; j<(1<<N); j++) {
            printf("%d ",A[i][j]);
        }
        printf("\n");
    }
}

// checks if A does not contain -1
int isComplete(int A[1<<N]) {
    int val = A[(1<<N)-1];
    int i;
    if (val==-1) return 0;
    A[(1<<N)-1] = -1;
    for (i=0; A[i] != -1; i++);
    A[(1<<N)-1] = val;
	return i==((1<<N)-1);

}


// returns 1 if A is the smallest permutation up to conjugation with matrices in commuting_matrices. We use the lexicographic order of the lookup tables.
int is_smallest_in_class(int A[1<<N]){
    int E;
    for (int i=0;i<N_COMMUTING_MATRICES_B;i++) {
        for (int j=0; j<N_COMMUTING_MATRICES_A;j++) {
            for (int x=0;x<(1<<N);x++){
                if (A[lex[x]]==-1) goto NEXT_MAT; // We cannot say whether A is smaller than current matrix because of the undefined value
                if (A[commuting_matrices_B[i][lex[x]]] == -1) goto NEXT_MAT; // We cannot say whether A is smaller than current matrix because of the undefined value
                else {
                    E = commuting_matrices_A[j][A[commuting_matrices_B[i][lex[x]]]];
                    if (E < A[lex[x]]) return 0; // E is smaller, so A cannot be the smallest in class
                    if (E > A[lex[x]]) goto NEXT_MAT; // A is smaller, so consider next matrix
                } // consider next position
            }
            NEXT_MAT:;
        }
    }
    // A is the smallest of all, so return 1
    return 1;
}

// return 1 if y is not yet in the image of sbox. returns 0 if it is already taken (i.e., check whether y is contained in sbox or not). It only checks whether the cycle representative is already contained. The function does not work for y=0, since those are not included in the pos array
int isNotTaken(int y) {
	for (int i=0; i<15; i++) {
		if (repr[sbox[pos[i]]] == repr[y])
			return 0;
	}
	return 1;
}

// continue building DDT for new set point sbox[c]. Returns 0 if the S-box cannot be APN anymore (because one entry is set to 4), returns 1 otherwise
int addDDTInformation(int c) {
    int a;
    // for all input differences
	for (int i=0;i<127;i++) {
		a = evens[i];	// we only need to check those input differences corresponding to vectors with even Hamming Weight
                // if sbox[c + a] is set, increase DDT entry by 2
                if (sbox[c^a]!= -1) {
                    sbox_DDT[a][sbox[c]^sbox[c^a]]+=2;
                    // if set above 2, return that it cannot be APN anymore
                    if (sbox_DDT[a][sbox[c]^sbox[c^a]]>2) return 0;
                }

	}
	return 1;
}

// reducing DDT entries for new removed point sbox[c]. It reduces all entries that were increased during the APN check of addDDTInformation(c).
void removeDDTInformation(int c) {
	int a;
	for (int i=0;i<127;i++) {
		a = evens[i];
                if (sbox[c^a]!= -1) {
                    sbox_DDT[a][sbox[c]^sbox[c^a]]-=2;
                    if (sbox_DDT[a][sbox[c]^sbox[c^a]]==2) break; // important because addDDTInformation(c) did not increment any values above that point.
                }

	}
}

int maxDepth=0;

// recursive construction of the sbox.
void nextValue_with_filter(int depth, char* filename) {
	int xS,yS,j;
	iterations++;

	if (depth>maxDepth) maxDepth=depth; // to show maxDepth
	if ((time(NULL)-start)>=300) { // display the status every 300 seconds
		start=time(NULL);
		printf(" depth:%d\n",depth);
		printArray(sbox);
		printf("solutions so far:%d maxDepth:%d \n",solutions,maxDepth);
		printf("iterations:%lli  skipped: %d \nDepth skipped: \n",iterations,n_skipped);
		printSkippedArray(skipped);
		printf("running time: %li sec\n", time(NULL)-runtime);
		fflush(stdout);
	}

	//complete
	if (isComplete(sbox)) {
			solutions++;
			printf("found a new apn permutation: #%d\n",solutions);
			printArray2(sbox);
			fflush(stdout);
			FILE *fp = fopen(filename, "a");
			if (fp == NULL)
            {
                printf("Error opening file!\n");
		fflush(stdout);
            }
            else
            {
                fprintArray2(fp, sbox);
            }
            fclose(fp);
		return;
	}

	//not complete

	int x=pos[depth];	// get next free position in the look-up table
	for (int y=(1<<N)-1;y>=1;y--) {		// do not include y=0 since those values are already set in the beginning. The function isNotTaken(y) would not work for 0.
		if (isNotTaken(y)) { // we know the order of each element is 17
			xS = x;
			yS = y;
			for (int i=0; i<16; i++) { // we know the order of each element is 17
				// we set sbox[xS] = yS. Add this to the partial DDT information and check whether it can still be APN
				sbox[xS]=yS;
				if (!addDDTInformation(xS)) goto UNDO_WITH_FILTER; // undo if it can't be APN anymore
				xS=B[xS];
                    		yS=A[yS];
			}
			sbox[xS]=yS;
			if (!addDDTInformation(xS)) goto UNDO_WITH_FILTER; // undo if it can't be APN anymore
			
                	if (depth <= 2) { // if it is too deep, it is faster to omit the check for smallest representative.
                        	if (is_smallest_in_class(sbox)) nextValue_with_filter(depth+1, filename);
                        	else {n_skipped++; skipped[depth]++;}
                	}
                	else { // omit check for smallest representative
                    		nextValue_with_filter(depth+1, filename);
                	}
            
			//undo the changes (i.e., delete the set points of sbox from this step and undo the changes in the DDT)
			UNDO_WITH_FILTER:
			// we remove the points in reverse order to undo the DDT changes correctly
			do { 
				removeDDTInformation(xS);
				sbox[xS] = -1;
				xS = B_inv[xS];
			} while (xS != B_inv[x]);
		}
	}
}

void test_for_matrix(char* filename, int first, int second) {
	if (N!=8) {
		printf("invalid N\n");
		exit(-1);
	}
    int xS, yS;
    solutions=0;
    iterations=0;
    memset(sbox_DDT,0,(1<<N)*(1<<N)*sizeof(int));
    for (int i=0; i<128; i++) skipped[i]=0;
    FILE *fp = fopen(filename, "w");
    if (fp == NULL)
    {
            printf("Error opening file!\n");
    }
    fclose(fp);
    printf("Search S-boxes invariant for matrices: \n");
    printArray2(B);
    printArray2(A);
    printf("\n\n\n");
    fflush(stdout);
    for (int x=0;x<(1<<N);x++) {
        sbox[x]=-1;
    }

    // w.l.o.g., we can set sbox[0]=0
    sbox[0]=0;
    if (!addDDTInformation(0)) exit(0);

    sbox[1]=first;
    if (!addDDTInformation(1)) exit(0);
    xS=B[1];
    yS=A[first];
    // set the orbits accordingly
    for (int i=0;i<16;i++) {
                    sbox[xS]=yS;
                    if (!addDDTInformation(xS)) exit(0);
                    xS=B[xS];
                    yS=A[yS];
    }
    sbox[3]=second;
    if (!addDDTInformation(3)) exit(0);
    xS=B[3];
    yS=A[second];
    // set the orbits accordingly
    for (int i=0;i<16;i++) {
                    sbox[xS]=yS;
                    if (!addDDTInformation(xS)) exit(0);
                    xS=B[xS];
                    yS=A[yS];
    }
    if (!is_smallest_in_class(sbox)) {
	exit(0);
    } 

    // start the recursive search
    nextValue_with_filter(2, filename);
}

int main(int argc, char* argv[])
{
	if (argc != 4) {
        	printf("Usage: %s <outfile> <SBox[1]> <SBox[3]> ...\n", argv[0]);
        	return -1;
    	}
	srand (time(NULL));
	start = time(NULL);
	runtime = time(NULL);
	
	
	if ((atoi(argv[2]) == 0) or (atoi(argv[3]) == 0)) {
		printf("No permutation ..\n");
		return -1;
	}
	if (repr[atoi(argv[2])] == repr[atoi(argv[3])]) {
		printf("No permutation ..\n");
		return -1;
	}
	test_for_matrix(argv[1],atoi(argv[2]),atoi(argv[3]));
	printf("\nFinished. Total running time: %li sec\n", time(NULL)-runtime);
	exit(0);
}
