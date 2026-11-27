import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from proyecciones.conversion import convertir_csv, latlon_a_xy
from proyecciones.dirichlet import ProyeccionOptima, obtener_modelo


def test_csv_mayusculas_no_sobrescribe(tmp_path):
    entrada = tmp_path / 'clientes.CSV'
    texto = 'lat,lon\n19.4326,-99.1332\n'
    entrada.write_text(texto)
    subprocess.run([sys.executable, '-m', 'proyecciones.conversion', str(entrada)],
                   cwd=Path(__file__).resolve().parents[1] / 'src', check=True)
    assert entrada.read_text() == texto
    assert (tmp_path / 'clientes_xy.csv').exists()
    with pytest.raises(ValueError, match='diferente'):
        convertir_csv(entrada, entrada)


@pytest.mark.parametrize('nombre', ['mercator', 'lambert', 'recomendada', 'optima'])
@pytest.mark.parametrize('lat,lon', [(100, -100), (90, 0), (20, 181), (np.inf, 0)])
def test_coordenadas_invalidas(nombre, lat, lon):
    with pytest.raises(ValueError):
        latlon_a_xy(lat, lon, nombre)


def test_persistencia_dominio_y_conversion(tmp_path):
    original = obtener_modelo()
    original.guardar(tmp_path / 'modelo.npz')
    cargado = ProyeccionOptima.cargar(tmp_path / 'modelo.npz')
    cargado.guardar(tmp_path / 'segundo.npz')
    assert np.allclose(original.proyectar(20, -100), cargado.proyectar(20, -100))
    assert cargado.contiene(20, -100)[0]
    assert not cargado.contiene(32.3571137518, -117.124498078)[0]
    with pytest.warns(UserWarning, match='fuera'):
        latlon_a_xy(32.3571137518, -117.124498078, 'optima')


def test_modelo_compartido():
    assert obtener_modelo() is obtener_modelo('optima')


def test_interfaces_numericas_documentadas():
    p = obtener_modelo()
    a, b = 27034, 26828
    w = (p.w[a] + p.w[b]) / 2
    from proyecciones.dirichlet import latlon_de
    from proyecciones.datos import haversine
    direction = (p.w[b]-p.w[a])/abs(p.w[b]-p.w[a])
    lat, lon = latlon_de(w + np.array([-1, 1])*direction*1e-10)
    x, y = p.proyectar(lat, lon)
    reference = haversine(lat[0], lon[0], lat[1], lon[1])
    assert abs(np.hypot(x[1]-x[0], y[1]-y[0])-reference) < .01
