#include <string.h>
#include <stdio.h>

int main(void) {
    unsigned char buf1[33];
    unsigned char buf2[33];
    int ret;
    void *(*volatile const volatile_memset)(void *, int, size_t) = memset;

    memset(buf1, 5, 33);
    volatile_memset(buf2, 5, 33);
    ret = memcmp(buf1, buf2, 33);
    printf("%d\n", ret);

    return 0;
}
