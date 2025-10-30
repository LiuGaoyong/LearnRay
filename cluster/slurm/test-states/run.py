import time

import ray
from ray._private.state import available_resources_per_node

ray.init()
for _ in range(6):
    print("=" * 54)
    print(ray.cluster_resources())
    print(ray.available_resources())
    print(available_resources_per_node())
    print("=" * 54)
    time.sleep(30)
