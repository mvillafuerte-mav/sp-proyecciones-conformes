# Proyecciones conformes óptimas para México

Situación problema de Fismat: wiki técnica sobre qué proyección conviene usar
para analizar datos geoespaciales de México (Mercator vs. Lambert de INEGI vs.
la proyección óptima de Chebyshev).

La consigna completa está en [`docs/consigna.md`](docs/consigna.md).

## Estructura

```
.
├── wiki/                 Documento final en LaTeX
│   ├── main.tex          Archivo principal (solo junta las secciones)
│   ├── secciones/        Una sección por archivo
│   ├── figuras/          Figuras finales que se usan en el documento
│   └── referencias.bib   Bibliografía (formato IEEE)
├── src/proyecciones/     Código en Python
│   ├── constantes.py     R, paralelos de INEGI, etc. (todos usamos este)
│   ├── clasicas.py       Mercator y Lambert
│   ├── distorsion.py     Factor de escala σ y distorsión δ
│   ├── dirichlet.py      Solver de Laplace + proyección óptima
│   ├── datos.py          Contorno de México y datos de INEGI
│   └── conversion.py     Herramienta para el cliente: (lat, lon) -> (x, y)
├── notebooks/            Pruebas y exploración (no va al documento)
├── data/raw/             Datos descargados (no se suben a git)
├── data/processed/       Datos limpios (no se suben a git)
├── resultados/           Figuras y tablas generadas por el código
├── tests/                Pruebas (comparación con pyproj, caso conocido)
└── docs/                 Consigna y notas del equipo
```

## Quién hace qué

| Persona | Rol | Secciones | Código |
|---|---|---|---|
| 1 | Teoría | Fundamentos, Pasos 1 y 2, teoría del Paso 4 | — |
| 2 | Solver numérico | Paso 4 numérico, método numérico | `dirichlet.py` |
| 3 | Datos y figuras | Paso 5, figuras, herramienta para el cliente | `datos.py`, `conversion.py` |
| 4 | Contexto, fuentes y edición | Resumen, intro, Paso 3, limitaciones, desarrollos, conclusiones, referencias | `distorsion.py`, `tests/` |

## Cómo trabajar

1. Clonar e instalar:
   ```bash
   git clone <url-del-repo>
   cd sp-proyecciones-conformes
   python -m venv .venv && source .venv/bin/activate
   pip install -r requirements.txt
   ```
2. Antes de empezar: `git pull`.
3. Cada quien edita **sus** archivos de `wiki/secciones/` para no pisarnos.
4. Commits chicos y con mensaje claro, ej. `paso 2: derivación de n`.
5. Todo número que aparezca en la wiki tiene que salir del código, no se copia a mano.

## Datos

Los datos no se suben al repo. Para descargarlos y procesarlos (tarda unos
minutos, hace una consulta a la API del INEGI por municipio):

```bash
cd src
python -m proyecciones.datos
```

Esto deja en `data/processed/`:
- `contorno_mexico.csv`: contorno de México continental (sin islas), a partir
  del Marco Geoestadístico del INEGI, simplificado a ~1 km.
- `cabeceras.csv`: cabecera de cada municipio con lat, lon, región
  (Norte/Centro/Sur por entidad) y si cae dentro del contorno continental.

Si sale un error de certificado SSL en Mac, correr una vez
`/Applications/Python\ 3.*/Install\ Certificates.command`.

## Resultados

Cada número de la wiki sale de alguno de estos scripts (desde `src/`):

```bash
python -m proyecciones.distorsion     # δ y errores de escala (Paso 3)
python -m proyecciones.verificacion   # comparación con PROJ y con GRS80
python -m proyecciones.dirichlet      # proyección óptima (Paso 4) y su validación
python -m proyecciones.efecto         # efecto en los datos: distancias, k-means, vecinos, áreas (Paso 5)
```

Las tablas quedan en `resultados/tablas/`.

## Compilar la wiki

Con Overleaf: subir la carpeta `wiki/` (o conectar el repo) y compilar `main.tex`.

En local:
```bash
cd wiki
latexmk -pdf main.tex
```

## Convenciones

- Tierra esférica, R = 6371.0088 km.
- Las funciones que reciben puntos (lat, lon) los toman en grados; los
  parámetros internos como φ1, φ2 van en radianes.
- Notación: φ latitud, λ longitud, q(φ) latitud isométrica, w = λ + i q,
  σ factor de escala, Ω región de México, δ(Ω) distorsión.
- Máximo 3 figuras, con pie y fuente.
- Referencias IEEE: si usas una fuente, agrégala a `referencias.bib` con DOI o enlace.
- Extensión: 2000–3000 palabras sin ecuaciones, código ni referencias.
