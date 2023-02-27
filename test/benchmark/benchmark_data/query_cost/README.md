## Query Cost Benchmark

Source code and machine information used to determine the query costs in `include/chopper/layout/ibf_query_cost.hpp`.

Related PR: https://github.com/seqan/chopper/pull/138

We used the **mean** values.

The original benchmark is in
[Raptor](https://github.com/seqan/raptor/tree/fe59a68b7174cbdf34db31d72d5224c3f33ac14b/test/performance), which provides
the build infrastructure to repeat this test. An exact copy of Raptor's `bin_influence_benchmark.cpp` as it was used
for benchmarking can be found in this directory.

### Benchmark Call
```bash
bin_influence_benchmark --benchmark_repetitions=5 \
                        --benchmark_enable_random_interleaving=true \
                        --benchmark_out=results.json \
                        --benchmark_out_format=json \
                        --benchmark_min_time=0.01
```

### Operating System
```bash
$ uname -a # Nero
Linux nero 5.10.0-17-amd64 SMP Debian 5.10.136-1 (2022-08-13) x86_64 GNU/Linux
```

### Compiler version
```
$ g++-12.2 -v
Using built-in specs.
COLLECT_GCC=/gcc-12.2/bin/g++-12.2
COLLECT_LTO_WRAPPER=/gcc-12.2/libexec/gcc/x86_64-linux-gnu/12/lto-wrapper
Target: x86_64-linux-gnu
Configured with: /gcc/configure --build=x86_64-linux-gnu --disable-bootstrap --disable-libssp --disable-libstdcxx-pch --disable-libunwind-exceptions --disable-multilib --disable-werror --enable-cet=auto --enable-clocale=gnu --enable-__cxa_atexit --enable-default-pie --enable-default-ssp --enable-gnu-indirect-function --enable-gnu-unique-object --enable-install-libiberty --enable-languages=c,c++,lto --enable-linker-build-id --enable-lto --enable-objc-gc=auto --enable-plugin --enable-shared --enable-threads=posix --host=x86_64-linux-gnu --no-create --no-recursion --prefix=/gcc-12.2 --program-suffix=-12.2 --target=x86_64-linux-gnu --with-abi=m64 --with-arch-32=i686 --with-bugurl=http://gcc.gnu.org/bugzilla/ --with-gcc-major-version-only --with-linker-hash-style=gnu --without-cuda-driver --with-system-zlib --with-tune=generic
Thread model: posix
Supported LTO compression algorithms: zlib zstd
gcc version 12.2.0 (GCC)
```

### Compiler flags
Default Raptor flags.
```
CXX_FLAGS = -O3 -DNDEBUG -pedantic -Wall -Wextra -Werror -std=c++20 -fopenmp -fopenmp-simd -DSIMDE_ENABLE_OPENMP -Wno-psabi -march=native -flto=auto -s -pthread
```

### CPU flags
```
$ lscpu
Architecture:                    x86_64
CPU op-mode(s):                  32-bit, 64-bit
Byte Order:                      Little Endian
Address sizes:                   43 bits physical, 48 bits virtual
CPU(s):                          128
On-line CPU(s) list:             0-127
Thread(s) per core:              2
Core(s) per socket:              64
Socket(s):                       1
NUMA node(s):                    1
Vendor ID:                       AuthenticAMD
CPU family:                      23
Model:                           49
Model name:                      AMD EPYC 7702P 64-Core Processor
Stepping:                        0
Frequency boost:                 enabled
CPU MHz:                         1497.889
CPU max MHz:                     3353.5149
CPU min MHz:                     1500.0000
BogoMIPS:                        3999.87
Virtualization:                  AMD-V
L1d cache:                       2 MiB
L1i cache:                       2 MiB
L2 cache:                        32 MiB
L3 cache:                        256 MiB
NUMA node0 CPU(s):               0-127
Vulnerability Itlb multihit:     Not affected
Vulnerability L1tf:              Not affected
Vulnerability Mds:               Not affected
Vulnerability Meltdown:          Not affected
Vulnerability Mmio stale data:   Not affected
Vulnerability Retbleed:          Mitigation; untrained return thunk; SMT enabled with STIBP protection
Vulnerability Spec store bypass: Mitigation; Speculative Store Bypass disabled via prctl and seccomp
Vulnerability Spectre v1:        Mitigation; usercopy/swapgs barriers and __user pointer sanitization
Vulnerability Spectre v2:        Mitigation; Retpolines, IBPB conditional, STIBP always-on, RSB filling, PBRSB-eIBRS Not affected
Vulnerability Srbds:             Not affected
Vulnerability Tsx async abort:   Not affected
Flags:                           fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush mmx fxsr sse sse2 ht syscall nx mmxext fxsr_opt pdpe1gb rdtscp lm constant_tsc rep_good nopl nonstop_tsc cpuid extd_apicid aperfmperf pni pclmulqdq monitor ssse3 fma cx16 sse4_1 sse4_2 movbe popcnt aes xsave avx f16c rdrand lahf_lm cmp_legacy svm extapic cr8_legacy abm sse4a misalignsse 3dnowprefetch osvw ibs skinit wdt tce topoext perfctr_core perfctr_nb bpext perfctr_llc mwaitx cpb cat_l3 cdp_l3 hw_pstate sme ssbd mba sev ibrs ibpb stibp vmmcall sev_es fsgsbase bmi1 avx2 smep bmi2 cqm rdt_a rdseed adx smap clflushopt clwb sha_ni xsaveopt xsavec xgetbv1 xsaves cqm_llc cqm_occup_llc cqm_mbm_total cqm_mbm_local clzero irperf xsaveerptr rdpru wbnoinvd arat npt lbrv svm_lock nrip_save tsc_scale vmcb_clean flushbyasid decodeassists pausefilter pfthreshold avic v_vmsave_vmload vgif umip rdpid overflow_recov succor smca
```

### Memory
```
$ free -g
               total
Mem:             995
Swap:              7
```
