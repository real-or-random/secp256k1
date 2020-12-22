with (import <nixpkgs> {}).pkgsi686Linux;
mkShell {
   buildInputs = [
       bash file pkgconfig autoconf automake libtool gmp valgrind clang gcc
   ];
}
