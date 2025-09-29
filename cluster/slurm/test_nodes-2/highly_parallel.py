import random
import time
from fractions import Fraction

import ray
from loguru import logger


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


if __name__ == "__main__":
    import sys

    if len(sys.argv) <= 1:
        BATCHES = 100
    else:
        BATCHES = int(sys.argv[1])

    t0 = time.time()
    # ray.init()
    logger.info(f"Ray start by {time.time() - t0:.2}s.")
    logger.info("-" * 32)

    start = time.time()
    SAMPLE_COUNT = 1000 * 1000
    future = pi4_sample.remote(SAMPLE_COUNT)
    logger.info("Doing Once: ")
    pi4 = ray.get(future)
    end = time.time()
    dur = end - start
    logger.info(f"Running {SAMPLE_COUNT} tests took {dur:.2f} seconds")
    logger.info("-" * 32)

    logger.info(f"Doing {BATCHES} batches")
    start = time.time()
    results = []
    for _ in range(BATCHES):
        results.append(pi4_sample.remote(SAMPLE_COUNT))
    output = ray.get(results)
    FULL_SAMPLE_COUNT = BATCHES * SAMPLE_COUNT
    logger.info(
        f"Running {FULL_SAMPLE_COUNT} tests took "
        f"{time.time() - start:.2f} seconds."
    )
    logger.info("-" * 32)

    start = time.time()
    SAMPLE_COUNT = 1000 * 1000
    future = pi4_sample.remote(SAMPLE_COUNT)
    logger.info("Doing Once: ")
    pi4 = ray.get(future)
    end = time.time()
    dur = end - start
    logger.info(f"Running {SAMPLE_COUNT} tests took {dur:.2f} seconds")
    logger.info("-" * 32)
