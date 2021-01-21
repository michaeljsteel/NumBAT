
# types.py is a subroutine of NumBAT that contains various property datatypes. 

# Copyright (C) 2017 Bjorn Sturmberg, Kokou Dossou.

# NumBAT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

from enum import IntEnum, auto

class PointGroup(IntEnum):
  Unknown = 1
  C2V = auto()
  C3V = auto()
  C6V = auto()

class SymRep(IntEnum):
#E = 1
  A = 0
  B1 = auto()
  B2 = auto()
  B3 = auto()
  Unknown = auto()


class QAcMethod(IntEnum):
  NotSet=1
  Fixed=auto()
  Intrinsic=auto()
