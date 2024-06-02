

# Write-Host "VCPKG_ROOT=$env:VCPKG_ROOT"
# Write-Host "CMAKE_ARGS=$env:CMAKE_ARGS"
# Write-Host "-DCMAKE_TOOLCHAIN_FILE=$env:VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake"

vcpkg install --triplet x64-windows gsl boost-headers
vcpkg integrate install

New-Item -Path "Env:\CMAKE_ARGS" -Value "-DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake"

# Write-Host "VCPKG_ROOT=$env:VCPKG_ROOT"
Write-Host "CMAKE_ARGS=$env:CMAKE_ARGS"
# Write-Host "-DCMAKE_TOOLCHAIN_FILE=$env:VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake"
