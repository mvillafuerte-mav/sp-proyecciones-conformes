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
│   ├── verificacion.py   Comparación con PROJ y con el elipsoide
│   ├── efecto.py         Efecto en distancias, k-means, vecinos y áreas
│   ├── figuras.py        Tabla comparativa y figuras
│   └── conversion.py     Herramienta para el cliente: (lat, lon) -> (x, y)
├── notebooks/            Pruebas y exploración (no va al documento)
├── data/raw/             Datos descargados (no se suben a git)
├── data/processed/       Datos limpios (no se suben a git)
├── resultados/           Tablas, figuras y modelos generados por el código
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
python -m proyecciones.figuras        # tabla comparativa y las 3 figuras de la wiki
```

Las tablas quedan en `resultados/tablas/`, las figuras en `wiki/figuras/` y la
tabla comparativa en `wiki/tablas/`.

## Herramienta para el cliente

Convierte (lat, lon) a (x, y) en metros. No necesita descargar datos: usa los
modelos ya ajustados de `resultados/modelos/`.

```bash
cd src
python -m proyecciones.conversion                               # ejemplo con 3 ciudades
python -m proyecciones.conversion clientes.csv clientes_xy.csv   # CSV con columnas lat, lon geográficas (EPSG:4326)
python -m proyecciones.conversion clientes.csv salida.csv --proyeccion optima_envolvente
```

Desde Python: `from proyecciones.conversion import latlon_a_xy`.

La entrada geográfica se interpreta en EPSG:4326. Se rechazan infinitos,
latitudes fuera de (-90, 90), longitudes fuera de [-180, 180] y la escritura
del CSV sobre su propio archivo de entrada. Las filas incompletas quedan vacías.
Los modelos NPZ incluyen el polígono de ajuste para avisar sobre extrapolación.
La cota teórica es de escala local sobre la esfera: las distancias entre pares
se verifican aparte. La evaluación de la óptima tiene saltos de hasta unos
6 mm en interfaces de celdas y no debe tratarse como exactamente conforme
a escala milimétrica.

En Windows, desde la raíz, ejecutar pruebas con
`.\.venv\Scripts\python.exe -m pytest -q`.
Para ejecutar módulos, entrar en `src` y usar
`..\.venv\Scripts\python.exe -m proyecciones.conversion`.
La revisión se realizó con Python 3.12.10; las versiones de las bibliotecas
quedaron registradas en `requirements-lock.txt`.

## Compilar la wiki

Con Overleaf: subir la carpeta `wiki/` (o conectar el repo) y compilar `main.tex`.

En local:
```bash
cd wiki
latexmk -xelatex main.tex
```

El documento usa Times New Roman y requiere XeLaTeX (por ejemplo,
`latexmk -xelatex main.tex`) o `tectonic main.tex`. La portada y el logo
proceden de la plantilla aportada, conservada en `wiki/recursos/`.
Los cuatro campos de integrantes se editan en `wiki/datos_portada.tex`.

## Convenciones

- Tierra esférica, R = 6371.0088 km.
- Las funciones que reciben puntos (lat, lon) los toman en grados; los
  parámetros internos como φ1, φ2 van en radianes.
- Notación: φ latitud, λ longitud, q(φ) latitud isométrica, w = λ + i q,
  σ factor de escala, Ω región de México, δ(Ω) distorsión.
- Máximo 3 figuras, con pie y fuente.
- Referencias IEEE: si usas una fuente, agrégala a `referencias.bib` con DOI o enlace.
- Extensión actualizada: hasta 4000 palabras sin ecuaciones, código ni referencias.

## Wiki revisada

El PDF compilado está en [wiki/main.pdf](wiki/main.pdf). La fuente principal
es [wiki/main.tex](wiki/main.tex) y los integrantes se completan en
[wiki/datos_portada.tex](wiki/datos_portada.tex). El límite actualizado es de
4000 palabras. Los teoremas auxiliares están en el anexo.

La Figura 1 combina la construcción estereográfica con las retículas. Para
regenerarla, instalar opcionalmente PyMuPDF (`python -m pip install pymupdf`)
y ejecutar `python scripts/generar_geometria.py` desde la raíz. El script
usa la retícula existente de `wiki/figuras/fig1_reticula.pdf`.
