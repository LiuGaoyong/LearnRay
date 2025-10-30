#!/bin/bash

for i in $(ls test-ray-*); do
    sbatch $i
done