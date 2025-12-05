{
  description = "Dev shell for MolSim Praktikum";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs {
          inherit system;
        };

        llvm = pkgs.llvmPackages_20;
      in
      {
        devShells.default = pkgs.mkShell.override { stdenv = llvm.libcxxStdenv; } {
          packages = with pkgs; [
            (llvm.clang-tools.override { enableLibcxx = true; })
            llvm.libcxxClang
            gcc14
            cmake
            doxygen
            paraview
            vtk
            cmake-format
            llvm.openmp
            (texliveTeTeX.withPackages (ps: [ ps.newunicodechar ]))
            ghostscript_headless
            lcov
            hyperfine
          ];
          shellHook = ''
            export NIX_CFLAGS_COMPILE="$NIX_CFLAGS_COMPILE -B${llvm.libcxxClang.libcxx}/lib";
            export CLICOLOR=0;
          '';
        };
      }
    );
}
