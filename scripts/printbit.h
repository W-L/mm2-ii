#ifndef MM2II_PRINTBIT_H
#define MM2II_PRINTBIT_H

#include <stdio.h>
#include <inttypes.h>


void printbytes(const void *b, int size);

void printbits(const uint8_t byte);



void printbits(const uint8_t byte)
{
    char out[9];

    uint8_t mask = 1;
    uint8_t tmp = 0;

    int i = 0;
    for (; i < 8; ++i) {
        tmp = (byte >> i) & mask;
        out[7 - i] = tmp + '0';
    }

    out[8] = '\0';
    printf("%s", out);
}

void printbytes(const void *b, int size)
{
    uint8_t *bytes = (uint8_t *) b;
    uint8_t *byte = bytes;
    int i = 0;
    for (; i < size; ++i) {
        printf("%p:\t", byte);
        printbits(*byte);
        printf("\n");
        byte++;
    }
}


#endif //MM2II_PRINTBIT_H
