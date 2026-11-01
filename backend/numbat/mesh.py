"""Mesh objects: ``from numbat.mesh import FemMesh``."""

from numbat_core.geometry.usermesh import UserGeometryBase
from numbat_core.mesh.fem import FemMesh
from numbat_core.mesh.mail import MailData

__all__ = ["FemMesh", "MailData", "UserGeometryBase"]
