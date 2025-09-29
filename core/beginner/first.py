import time
from typing import Any

import ray
import ray.actor

database = [
    "Learning",
    "Ray",
    "Flexible",
    "Distributed",
    "Python",
    "for",
    "Machine",
    "Learning",
]


def retrieve(item) -> tuple[Any, Any]:
    """Retrieve data from database."""
    time.sleep(item / 10.0)
    return item, database[item]


def print_runtime(input_data, start_time) -> None:
    """Print runtime and data."""
    print(f"Runtime: {time.time() - start_time:.2f} seconds, data:")
    print(*input_data, sep="\n")


start = time.time()
data = [retrieve(item) for item in range(8)]
print("# Serial Runner: ")
print_runtime(data, start)
print("-" * 32)
##################################################################
# The serial runner is generally not efficient.


# @ray.remote
# def retrieve_task(item) -> tuple[Any, Any]:
#     """Retrieve data from database by Ray Task."""
#     return retrieve(item)
#
# This is equivalent to the following code block.
retrieve_task = ray.remote(retrieve)

t0 = time.time()
ray.init()
print(f"Ray init time: {time.time() - t0:.2f} seconds")
start = time.time()
object_references = [retrieve_task.remote(item) for item in range(8)]
data = ray.get(object_references)
print("# Ray Runner: ")
print_runtime(data, start)
print("-" * 32)
###################################################################
# The Ray runner by task is much efficent.


db_objref = ray.put(database)


@ray.remote
def retrieve_task_use_db(item, db) -> tuple[Any, Any]:
    """Retrieve data from database by Ray Task & Object Store."""
    time.sleep(item / 10.0)
    return item, db[item]


start = time.time()
object_references = [retrieve_task.remote(item) for item in range(8)]
data = ray.get(object_references)
print("# Ray Runner & Object Store: ")
print_runtime(data, start)
print("-" * 32)
###################################################################
# use object store to store database. It can be more efficient.


@ray.remote
class DataTracker:  # noqa: D101
    def __init__(self) -> None:  # noqa: D107
        self._counts = 0

    def increment(self) -> None:  # noqa: D102
        self._counts += 1

    def counts(self) -> int:  # noqa: D102
        return self._counts


@ray.remote
def retrieve_tracker_task(item: int, tracker, db) -> tuple[Any, Any]:
    """Retrieve data from database by Ray Task & Object Store.

    This function will communite with tracker to count the number of tasks.
    """
    time.sleep(item / 10.0)
    tracker.increment.remote()
    return item, db[item]


tracker = DataTracker.remote()

object_references = [
    retrieve_tracker_task.remote(item, tracker, db_objref) for item in range(8)
]
data = ray.get(object_references)

print(data)
print(ray.get(tracker.counts.remote()))  # type: ignore
###############################################################


# ray.init()    : 初始化你的 Ray 集群。传入一个地址以连接到现有集群。
# @ray.remote   : 将函数转换为任务，将类转换为 actor。
# ray.put()     : 将值放入 Ray 的对象存储中。
# ray.get()     : 从对象存储中获取值。返回你存放的值或由任务或演员计算出的值。
# .remote()     : 在您的 Ray 集群上运行 actor 方法或任务，并用于实例化 actors。
# ray.wait()    : 返回两个对象引用列表，
#                    一个包含我们正在等待的已完成任务，
#                    另一个包含未完成的任务。
