#include <string.h>
#include <stdio.h>

int main(void) {
    unsigned char buf1[33];
    unsigned char buf2[33];
    int ret;

    memset(buf1, 5, 33);
    memset(buf2, 5, 33);
    ret = memcmp(buf1, buf2, 33);
    printf("%d\n", ret);

    return 0;
}
