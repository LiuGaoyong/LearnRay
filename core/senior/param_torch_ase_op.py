from typing import Callable

import torch
from ase import Atoms
from ase.calculators.calculator import Calculator
from torch import Tensor

FUNC_TROCH_SPE_TYPE = Callable[[torch.Tensor], torch.Tensor]


def get_torch_ase_op(atoms: Atoms) -> FUNC_TROCH_SPE_TYPE:
    name = atoms.symbols.formula.format("hill")
    name = f"{atoms.calc.__class__.__name__}_{name}"
    name = f"get_torch_ase_op::{name}"
    if not isinstance(atoms.calc, Calculator):
        raise RuntimeError(
            f"atoms.calc must be a Calculator, but got {type(atoms.calc)}"
        )
    atoms.calc.reset()

    @torch.library.custom_op(name, mutates_args=(), device_types="cpu")
    def func(x: Tensor) -> Tensor:
        atoms.positions = x.numpy()
        assert isinstance(atoms.calc, Calculator)
        atoms.calc.results.clear()
        e = atoms.get_potential_energy()
        return torch.tensor(e)

    @func.register_fake
    def _(x: Tensor) -> Tensor:
        return torch.tensor(0.0)

    def backward(ctx, grad) -> torch.Tensor:
        forces = atoms.get_forces()
        return grad * torch.from_numpy(forces)

    torch.library.register_autograd(name, backward)
    return func


if __name__ == "__main__":
    from ase.build import molecule
    from ase.calculators.emt import EMT

    print(__name__)

    CO = molecule("CO")
    CO.calc = EMT()
    EMT.calculate
    print(CO.get_potential_energy())
    print(CO.get_forces())

    func = get_torch_ase_op(CO)

    x = torch.from_numpy(CO.positions)
    x.requires_grad_(True)
    print("x=", x)
    y = func(x)
    print("y=", y)
    y.backward()
    print("grad=", x.grad)
