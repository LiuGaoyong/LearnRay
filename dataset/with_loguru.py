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
        if logfile is not None:
            self.__logger.add(logfile)

    def info(self, msg: str) -> None:
        self.__logger.info(msg)


if __name__ == "__main__":
    from pathlib import Path

    f = Path(__file__)
    name = f.name.split(".")[0]
    log = f.with_name(f"{name}_log.txt")
    logger = MyLogger.remote(log)

    ds = ray.data.range(100)
    ds2 = ds.map(lambda x: {"out": ray.get(logger.info.remote(f"{x}"))})  # type: ignore
    print(ds2.take_all())

    # 在Ray中使用loguru， 必须使用ray.remote把Logger类封装起来！！
