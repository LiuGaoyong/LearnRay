import math
import random
import time

import ray
import ray.actor

ray.init()
t0 = time.time()


@ray.remote
class ProgressActor:
    def __init__(self, total_num_samples: int):
        self.total_num_samples = total_num_samples
        self.num_samples_completed_per_task = {}

    def report_progress(self, task_id: int, num_samples_completed: int) -> None:
        self.num_samples_completed_per_task[task_id] = num_samples_completed

    def get_progress(self) -> float:
        return (
            sum(self.num_samples_completed_per_task.values())
            / self.total_num_samples
        )


@ray.remote
def sampling_task(
    num_samples: int,
    task_id: int,
    progress_actor: ray.actor.ActorHandle,
) -> int:
    num_inside = 0
    for i in range(num_samples):
        x, y = random.uniform(-1, 1), random.uniform(-1, 1)
        if math.hypot(x, y) <= 1:
            num_inside += 1

        # Report progress every 1 million samples.
        if (i + 1) % 1_000_000 == 0:
            # This is async.
            progress_actor.report_progress.remote(task_id, i + 1)

    # Report the final progress.
    progress_actor.report_progress.remote(task_id, num_samples)
    return num_inside


# Change this to match your cluster scale.
NUM_SAMPLING_TASKS = 10
NUM_SAMPLES_PER_TASK = 10_000_000
TOTAL_NUM_SAMPLES = NUM_SAMPLING_TASKS * NUM_SAMPLES_PER_TASK

# Create the progress actor.
progress_actor = ProgressActor.remote(TOTAL_NUM_SAMPLES)
# Create and execute all sampling tasks in parallel.
results = [
    sampling_task.remote(NUM_SAMPLES_PER_TASK, i, progress_actor)  # type: ignore
    for i in range(NUM_SAMPLING_TASKS)
]

# Query progress periodically.
while True:
    progress = ray.get(progress_actor.get_progress.remote())  # type: ignore
    print(f"Progress: {int(progress * 100)}%")
    if progress == 1:
        break
    time.sleep(0.25)
# 使用 actor 报告运行状态，会额外占用一部分时间。
#   其总耗时具体依赖于调用的 actor 方法的耗时:
#       actor null          :  8.73s
#       actor sleep 0.0s    : 10.92s
#       actor sleep 0.3s    :  8.89s
#       actor sleep 0.5s    :  8.92s
#       actor sleep 1.0s    :  9.41s

# Get all the sampling tasks results.
total_num_inside = sum(ray.get(results))
pi = (total_num_inside * 4) / TOTAL_NUM_SAMPLES
print(f"Estimated value of π is: {pi}")
print(f"Finished by {time.time() - t0:.2f}s")
