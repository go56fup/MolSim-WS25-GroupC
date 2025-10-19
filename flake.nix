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
            llvm.libcxxClang
            llvm.llvm
            gcc15
            cmake
            doxygen
            paraview
            vtk
            cmake-format
            llvm.openmp
          ];
        };
      }
    );
}
