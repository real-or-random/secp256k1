#include <string.h>

int main(void) {
    unsigned char buf1[33];
    unsigned char buf2[33];

    memset(buf1, 5, 33);
    memset(buf2, 5, 33);
    memcmp(buf1, buf2, 33);

    return 0;
}
