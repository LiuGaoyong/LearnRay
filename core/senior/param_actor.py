import torch
from torch import nn
from torch.autograd import Function

from ase import Atoms
from ase.calculators.emt import EMT



class Torch(Function):



class TorchModuleEMT(nn.Module):
    def __init__(self, atoms: Atoms) -> None:
        self.__numbers = torch.from_numpy(atoms.numbers)
        self.__cell = torch.from_numpy(atoms.cell.array)
        self.__pbc = torch.from_numpy(atoms.pbc)
        self.__atoms = atoms

    def forward(self, positions: torch.Tensor) -> torch.Tensor: