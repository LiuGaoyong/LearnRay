import random
import time
from fractions import Fraction

import ray

# Let's start Ray
ray.init()
print("-" * 32)


@ray.remote
def pi4_sample(sample_count) -> Fraction:
    """pi4_sample runs sample_count experiments.

    returns the fraction of time it was inside the circle.
    """
    in_count = 0
    for _ in range(sample_count):
        x = random.random()
        y = random.random()
        if x * x + y * y <= 1:
            in_count += 1
    return Fraction(in_count, sample_count)


SAMPLE_COUNT = 1000 * 1000
start = time.time()
future = pi4_sample.remote(SAMPLE_COUNT)
pi4 = ray.get(future)
end = time.time()
dur = end - start
print("Done Once")
print(f"Running {SAMPLE_COUNT} tests took {dur} seconds")
print("-" * 32)


FULL_SAMPLE_COUNT = 100 * 1000 * 1000  # * 1000  # 100 billion samples!
BATCHES = int(FULL_SAMPLE_COUNT / SAMPLE_COUNT)
print(f"Doing {BATCHES} batches")
start = time.time()
results = []
for _ in range(BATCHES):
    results.append(pi4_sample.remote(SAMPLE_COUNT))
output = ray.get(results)
print(f"Running {FULL_SAMPLE_COUNT} tests took {time.time() - start} seconds")
print("-" * 32)
