"""Factor de escala σ y distorsión δ(Ω) = max ln σ − min ln σ (Paso 3)."""


def sigma_mercator(lat):
    """σ = sec φ."""
    raise NotImplementedError


def sigma_lambert(lat, n, F):
    """σ = n F exp(−n q(φ)) / cos φ."""
    raise NotImplementedError


def distorsion(sigma):
    """δ a partir de los valores de σ evaluados en una malla dentro de Ω."""
    raise NotImplementedError


def error_max_escala(delta):
    """Error máximo de escala en %, con la escala centrada: 100 (e^{δ/2} − 1)."""
    raise NotImplementedError
