FROM nixos/nix

COPY ci/shell.nix /tmp
COPY ci/shell-i686.nix /tmp

# Run dummy command "true" in the nix-shell just to get the packages prepared.
RUN nix-shell -I channel:nixos-unstable /tmp/shell.nix --command true
RUN nix-shell -I channel:nixos-unstable /tmp/shell-i686.nix --command true
