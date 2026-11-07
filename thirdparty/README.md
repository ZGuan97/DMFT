# Third-party software

## nlohmann/json

- Pinned version: 3.11.2 (MIT license, bundled as `LICENSE.MIT`).
- Upstream: https://github.com/nlohmann/json/tree/v3.11.2
- Headers copied unchanged from the local conda package
  `nlohmann_json-3.11.2-h313beb8_0` on 2026-11-07.
- Location: `nlohmann-json-3.11.2/include/nlohmann/`.
- Header-only, compatible with C++11. Add
  `thirdparty/nlohmann-json-3.11.2/include` to the include path; no library link
  or runtime conda environment is required.
- Used by the cRPA interface independently of ComCTQMC.

## ComCTQMC

- Distribution source: [doi:10.17632/nvcc595t2b.1](https://doi.org/10.17632/nvcc595t2b.1).
- Source: local `ComDMFT-v2.0/src/CTQMC/V201`.
- Bundled version: `V201`.
- Copied on: 2026-10-10.
- Upstream commit: unknown; the source directory has no Git metadata.

# Eigen

Vendored Eigen 3.4.0 headers (Eigen/ directory, unchanged) and upstream COPYING.*
files. Downloaded 2026-09-29 from the official tagged source archive:

https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz

SHA-256 of downloaded archive:
`8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72`

The pinned 3.4 series supports the C++11 numerical tests and the existing C++14
production build. Only headers and licenses are retained; upstream tests,
benchmarks, documentation and unsupported modules are not needed here.
Add eigen-3.4.0 to the include path. No separate Eigen library build is required.
The new lattice and density_matrix modules use Eigen's native CPU algorithms,
without EIGEN_USE_MKL_ALL. Other DMFT-dev modules still require MKL.
