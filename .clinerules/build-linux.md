# Build Linux

To build opencopter on linux, first ensure that the `opencopter` conda environment is active. Then use the command
`dub build -c library-python312 -b debug-native-512 --compiler=ldc2 --force` to build the library.

