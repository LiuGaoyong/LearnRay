#!/bin/bash

which python
python ../launch.py                         \
    -n  4                                   \
    --exp-name  test_ray                    \
    --command "python highly_parallel.py"