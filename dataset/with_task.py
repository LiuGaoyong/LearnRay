from math import nan

import ray
from ase.build import molecule
from ase.collections import g2
from xtb.interface import Calculator, Param


def spe(name: str) -> float:
    try:
        atoms = molecule(name)
        xtb = Calculator(
            param=Param.GFNFF,
            numbers=atoms.numbers,  # type: ignore
            positions=atoms.positions,  # type: ignore
        )
        xtb.set_verbosity(0)
        result = xtb.singlepoint()
        return result.get_energy()
    except Exception:
        return nan


if __name__ == "__main__":
    # print(g2.names, spe("CO"))

    @ray.remote
    def fun(name: str) -> float:
        return spe(name)

    ds = ray.data.from_items(g2.names).map(
        lambda row: {"e": fun.remote(row["item"])},
    )
    print(ds.take_batch())
