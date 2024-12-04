# Project 2 | CSE 5242 (AU2024) 

To compile the `db5242.c` file, use the following command:
```shell
gcc -O3 -mavx2 -o db5242 db5242.c

```
To proceed on using the script provided below, make sure to change the permission by using below.
```shell
chmod +x (script_name)
```


# (g) Profiling the code of binary search
## Compile

```shell
#!/bin/bash
./db5242 10 50 50 100000 1000000
./db5242 10 50 50 100000 500000
./db5242 10 50 50 100000 100000
./db5242 100 50 50 100000 500000
./db5242 100 50 50 100000 100000
./db5242 100 50 50 100000 50000
./db5242 1000 50 50 100000 100000
./db5242 1000 50 50 100000 50000
./db5242 1000 50 50 100000 10000
./db5242 10000 50 50 100000 50000
./db5242 10000 50 50 100000 10000
./db5242 100000 50 50 100000 10000
./db5242 100000 50 50 100000 5000
./db5242 100000 50 50 100000 1000
./db5242 1000000 50 50 100000 1000
./db5242 1000000 50 50 100000 500
./db5242 10000000 50 50 100000 50
./db5242 10000000 50 50 100000 10
```

# (h) Profiling the code of band join
## Compile
```shell
#!/bin/bash

for Z in 10 100 1000 10000 100000 1000000; do
        echo "Testing with Z=$Z"

        for i in {1..10}
        do
                ./db5242 8192 16384 7001 $Z
        done
done
```

