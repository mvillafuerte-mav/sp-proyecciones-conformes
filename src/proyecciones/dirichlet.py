"""Proyección óptima de Chebyshev (Paso 4).

Pasos:
  1. Llevar el contorno de México al plano w.
  2. Resolver Δh = 0 en Ω con h = ln(R cos φ) + c en ∂Ω.
  3. Conjugada armónica h~.
  4. f' = exp(h + i h~) e integrar para obtener f.
"""


def resolver_dirichlet(w_frontera, valores_frontera):
    """Devuelve h (y de ser posible h~) como funciones evaluables en w."""
    # TODO Persona 2: elegir método (polinomio complejo, MFS o diferencias finitas)
    raise NotImplementedError


def proyectar_optima(lat, lon):
    """(lat, lon) -> (x, y) con la proyección óptima."""
    # TODO Persona 2
    raise NotImplementedError
