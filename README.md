# Project 2 | CSE 5242 (AU2024) 

To compile the `db5242.c` file, use the following command:
```shell
gcc -O3 -mavx2 -o db5242 db5242.c
```

# (g) Profiling the code of binary search
## Compile

```shell
for N in 10 100 1000 10000 100000 1000000; do
  echo "Testing with N=$N"
  for R in 10 100 1000 10000; do
    echo "Testing with R=$R"
    ./db5242 $N 0 0 0 $R
  done
done
```

# (h) Profiling the code of band join
## Compile
```shell
```

