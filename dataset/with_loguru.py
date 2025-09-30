import ray
from lark import logger
from loguru._logger import Core, Logger


@ray.remote
class MyLogger:
    def __init__(self, logfile: str | None = None) -> None:
        self.__logger = Logger(
            core=Core(),
            exception=None,
            depth=0,
            record=False,
            lazy=False,
            colors=False,
            raw=False,
            capture=True,
            patchers=[],
            extra={},
        )
        self.__logger.add("log.txt")

    def info(self, msg: str) -> None:
        self.__logger.info(msg)


if __name__ == "__main__":
    ds = ray.data.range(100)
    logger = MyLogger.remote()
    ds2 = ds.map(lambda x: {"out": ray.get(logger.info.remote(f"{x}"))})
    print(ds2.take_all())
