with (import <nixpkgs> {});
mkShell {
   buildInputs = [
       bash file pkgconfig autoconf automake libtool gmp valgrind clang gcc
   ];
}
