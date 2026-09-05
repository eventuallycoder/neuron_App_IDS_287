#include <stdio.h>
#include <stdlib.h>

int main() {
    printf("--- Simulating Scattered Memory (Fragmentation) ---\n");

    // 1. We malloc three separate integers.
    // malloc finds a free spot on the "Heap". These spots are rarely next to each other.
    // The computer puts "metadata" (bookkeeping info) between them.
    int *p1 = (int*)malloc(sizeof(int));
    int *p2 = (int*)malloc(sizeof(int));
    int *p3 = (int*)malloc(sizeof(int));

    *p1 = 10;
    *p2 = 20;
    *p3 = 30;

    // 2. Let's look at the addresses to see the GAPS.
    printf("Address of p1: %lld\n", (long long)p1);
    printf("Address of p2: %lld\n", (long long)p2);
    
    long long diff = (long long)p2 - (long long)p1;
    printf("Gap size     : %lld bytes\n", diff);
    printf("(We expected 4 bytes, but got %lld. The extra bytes are 'other memory things'.)\n\n", diff);

    // 3. THE FAIL: Trying to use pointer arithmetic to jump the gap.
    // If we pretend this is an array and just do p1 + 1...
    int *ghostPointer = p1 + 1; 

    printf("--- The Pointer Arithmetic Fail ---\n");
    printf("Target Value (p2) : %d\n", *p2);
    printf("Pointer Guess (p1+1): %d\n", *ghostPointer);

    if (ghostPointer == p2) {
        printf("RESULT: Success! (This almost never happens)\n");
    } else {
        printf("RESULT: Failed. p1+1 landed in the gap, not on p2.\n");
    }

    // Clean up
    free(p1);
    free(p2);
    free(p3);

    return 0;
}