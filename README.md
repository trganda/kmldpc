# kmldpc

Author: marrow1203@gmail.com

This project was develop for lab works. Using a hacked kmeans algorithm for bind-detect and solve the ambiguity.

## Compile && Simulate

compile:

```shell
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make -j 6
```

run:

Before running this binary, you need to edit the config file (`config.toml`) under the `build/kmldpc`.

```shell
cd build/kmldpc
./kmldpc
```

## System Model

y = hx + w, where x is the modulated symbols, it's multiplied by h () and plus to additive gaussia noise.