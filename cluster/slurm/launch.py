#!usr/bin/env python3

import argparse
import subprocess
import sys
import time
from pathlib import Path

PORT = "$PORT"
JOB_NAME = "$JOB_NAME"
NUM_NODES = "$NUM_NODES"
NCPUS_PER_NODE = "$NCPUS_PER_NODE"
PARTITION_SUBMIT = "$PARTITION_SUBMIT"
COMMAND_PLACEHOLDER = "$COMMAND_PLACEHOLDER"
template_file = Path(__file__).parent.joinpath("./template.sh")


def parse_args() -> argparse.Namespace:  # noqa: D103
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--exp-name",
        type=str,
        required=True,
        help="The job name and path to logging file (exp_name.log).",
    )
    parser.add_argument(
        "--num-nodes",
        "-n",
        type=int,
        default=1,
        help="Number of nodes to use.",
    )

    parser.add_argument(
        "--cpu-per-node",
        "-cpu",
        type=int,
        default=0,
        help="Number of cpu for each nodes to use.",
    )

    parser.add_argument(
        "--partition",
        "-p",
        type=str,
        default="ihicnormal",
    )
    parser.add_argument(
        "--port",
        "-P",
        type=str,
        default="6379",
    )

    parser.add_argument(
        "--command",
        type=str,
        required=True,
        help="The command you wish to execute. For example: "
        " --command 'python test.py'. "
        "Note that the command must be a string.",
    )

    parser.add_argument(
        "--submit",
        action="store_true",
        help="Wether submit this jos.",
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    job_name = "{}_{}".format(
        args.exp_name,
        time.strftime("%m%d-%H%M", time.localtime()),
    )

    # ===== Modified the template script =====
    with open(template_file, "r") as f:
        text = f.read()
    text = text.replace(JOB_NAME, job_name)
    text = text.replace(NUM_NODES, str(args.num_nodes))
    text = text.replace(PARTITION_SUBMIT, args.partition)
    if int(args.cpu_per_node) <= 0:
        s = f"#SBATCH --cpus-per-task={NCPUS_PER_NODE}"
        text = text.replace(s, "")
    else:
        text = text.replace(NCPUS_PER_NODE, str(args.cpu_per_node))
    text = text.replace(COMMAND_PLACEHOLDER, str(args.command))
    text = text.replace(PORT, str(args.port))

    # ===== Save the script =====
    script_file = "{}.sh".format(job_name)
    with open(script_file, "w") as f:
        f.write(text)

    # ===== Submit the job =====
    if args.submit:
        print("Starting to submit job!")
        subprocess.Popen(["sbatch", script_file])
        print("Job submitted! ")
        sys.exit(0)
