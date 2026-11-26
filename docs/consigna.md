# Situación Problema

En esta unidad de formación abordarás situaciones problema que te permitirán desarrollar las subcompetencias asignadas a este materia y alcanzar los niveles de dominio esperados.

## Consultoría

Una consultoría de prestigio recibe a una empresa de analítica geoespacial. La empresa procesa millones de registros con coordenadas (latitud, longitud) en todo México: ubicaciones de comercios, rutas de reparto, sensores y clientes. Para hacer agrupamiento, búsqueda de vecinos cercanos, mapas de calor y optimización de rutas, necesita trabajar en coordenadas planas (x, y), donde las distancias euclidianas sean confiables.

El problema es que la Tierra no se puede aplanar sin deformarla. Hoy la empresa usa la proyección de Mercator "porque es la de los mapas web" y sus analistas sospechan que las distancias en el norte del país salen infladas respecto a las del sur. El INEGI usa oficialmente otra proyección, la Cónica Conforme de Lambert (EPSG:6361, paralelos estándar 17.5° N y 29.5° N). El cliente pregunta: ¿cuál conviene usar?, ¿cuánto error introduce cada una? y ¿existe una proyección óptima para México?

Una matemática sénior de la consultoría reconoce la estructura del problema. Una proyección conforme es una función analítica de una coordenada compleja sobre la esfera; su deformación local está dada por $|f'|$, y el logaritmo de $|f'|$ es una función armónica. Un teorema de Chebyshev (1856) dice que la mejor proyección conforme de una región es la que tiene escala constante sobre su frontera. Encontrarla se reduce a resolver un problema de Dirichlet para la ecuación de Laplace.

El equipo de consultores (ustedes) debe entregar al cliente una wiki técnica que explique el método, resuelva el problema y cuantifique el efecto de la proyección sobre sus análisis de datos. La wiki debe responder también las tres preguntas que la consultoría hace en sus entrevistas:

1. ¿Cómo ayuda el análisis complejo, y en particular el mapeo conforme, a resolver un problema que a primera vista es de geometría y de datos?
2. ¿Qué se gana respecto a elegir una proyección "de catálogo" o corregir los datos por fuerza bruta?
3. ¿Cuáles son los límites del método y qué se usa cuando deja de ser aplicable?

---

# Descripción general

En equipos de **3 integrantes**, elaborarán una **entrada de wiki** de aproximadamente 4 **páginas (1500–3000 palabras)** sobre una **aplicación real** en física, matemáticas, ciencia de datos o en la industria donde se utilicen de manera fundamental los métodos matemáticos estudiados durante el semestre.

La entrada debe ser clara, concisa y estar escrita para un público con formación matemática, pero no necesariamente especializado en el área específica de aplicación. Pueden elegir cualquier aplicación **real y contemporánea**, siempre que haya una base matemática clara relacionada con el contenido del curso.

# Objetivo y descripción general

El objetivo es que el equipo aprenda mapeos conformes resolviendo con sus propias manos un problema real, y que comunique la solución en una wiki clara, rigurosa y útil para el cliente.

El hilo conductor son tres ideas del curso. Primero, una función analítica con $f'(z) \neq 0$ preserva ángulos y escala igual en todas direcciones, con factor $|f'(z)|$. Segundo, si $f$ es analítica y no se anula su derivada, $\ln|f'(z)|$ es la parte real de la función analítica $\ln f'(z)$ y por lo tanto es armónica. Tercero, una función armónica queda determinada por sus valores en la frontera, y su conjugada armónica permite reconstruir $f$.

$$
\nabla^2 \ln|f'(w)| = 0, \qquad f'(w) = \exp\big(h(w) + i\,\tilde{h}(w)\big)
$$

Por eso el diseño del mejor mapa se reduce a resolver la ecuación de Laplace en la región que ocupa México, calcular la conjugada armónica e integrar.

El trabajo se realiza en equipos de hasta 3 integrantes. No hay reporte por separado: todo, desde la teoría hasta el cómputo y las recomendaciones al cliente, va en la wiki.

---

# Producto esperado: la wiki para el cliente

Una entrada de wiki (estilo Wikipedia) dirigida al equipo del cliente: científicos de datos y analistas con formación matemática, pero no especialistas en análisis complejo ni en cartografía. Debe contener las secciones siguientes.

## 1. Título y resumen

Título claro y específico, por ejemplo: "Proyecciones conformes óptimas para el análisis de datos geoespaciales en México". El resumen da la respuesta al cliente en un párrafo: qué proyección recomendar, cuánto error introduce cada opción y cuánto se gana con la proyección óptima.

## 2. Introducción y contexto del problema

- Por qué los datos geoespaciales necesitan coordenadas planas y qué algoritmos dependen de distancias euclidianas (agrupamiento, vecinos cercanos, interpolación espacial, rutas).
- Por qué no existe un mapa perfecto: la curvatura de la esfera impide preservar a la vez distancias, ángulos y áreas (Euler, Gauss).
- Qué significa que una proyección sea conforme y por qué es la elección natural cuando importan las formas y las distancias locales.
- Las variables del problema: latitud $\varphi$, longitud $\lambda$, radio terrestre $R$, factor de escala $\sigma$, la forma y extensión de la región.

## 3. Fundamentos teóricos

- Funciones analíticas, ecuaciones de Cauchy–Riemann y por qué un mapeo analítico con $f'(z) \neq 0$ preserva ángulos y escala por el factor $|f'(z)|$.
- Funciones armónicas, conjugada armónica y el problema de Dirichlet; principio del máximo.
- La esfera de Riemann y la proyección estereográfica como mapeo conforme.
- Coordenada compleja sobre la esfera: con la latitud isométrica $q(\varphi) = \ln\tan(\pi/4 + \varphi/2)$ y $w = \lambda + i\,q$, la métrica de la esfera es conformemente plana. Toda proyección conforme es entonces $z = f(w)$ con $f$ analítica, y su factor de escala es

$$
ds^2 = R^2\cos^2\varphi\,|dw|^2, \qquad \sigma(w) = \frac{|f'(w)|}{R\cos\varphi}
$$

- Se recomienda una figura generada por el equipo: la retícula de meridianos y paralelos de México bajo dos proyecciones, mostrando que los ángulos rectos se conservan.

## 4. Solución del problema del cliente

Derivación completa y explicada, paso a paso (se acepta el modelo esférico de la Tierra):

**Paso 1. Mercator.** Mostrar que $f(w) = R\,w$ define la proyección de Mercator, que es conforme y que su factor de escala es $\sigma = \sec\varphi$. Mostrar que Mercator es, salvo constantes, el logaritmo complejo de la proyección estereográfica.

**Paso 2. Cónica Conforme de Lambert.** Mostrar que $f(w) = R\,F\,e^{i n w}$ es conforme, calcular su factor de escala y obtener la constante del cono $n$ a partir de dos paralelos estándar $\varphi_1$ y $\varphi_2$ donde $\sigma$ toma el mismo valor:

$$
\sigma(\varphi) = \frac{n F e^{-n q(\varphi)}}{\cos\varphi}, \qquad n = \frac{\ln(\cos\varphi_1/\cos\varphi_2)}{q(\varphi_2) - q(\varphi_1)}
$$

**Paso 3. Medir la distorsión.** Definir la distorsión de una proyección sobre una región $\Omega$ como la oscilación de $\ln\sigma$. Calcularla para Mercator y para la Lambert de INEGI ($\varphi_1 = 17.5°$, $\varphi_2 = 29.5°$) sobre el territorio de México.

$$
\delta(\Omega) = \max_\Omega \ln\sigma - \min_\Omega \ln\sigma
$$

**Paso 4. La proyección óptima de Chebyshev.** Mostrar que $h = \ln|f'|$ es armónica y que pedir $\sigma$ constante en la frontera equivale a un problema de Dirichlet. Explicar con el principio del máximo por qué esa condición minimiza $\delta$.

$$
\nabla^2 h = 0 \text{ en } \Omega, \qquad h = \ln(R\cos\varphi) + c \text{ en } \partial\Omega
$$

Resolver numéricamente el problema en el plano $w$ con el contorno de México, calcular la conjugada armónica, reconstruir $f' = \exp(h + i\,\tilde{h})$ e integrar para obtener $f$. Comparar la distorsión $\delta$ de la proyección óptima con la de Mercator y la de Lambert.

**Paso 5. Efecto en los datos del cliente.** Con un conjunto de datos real de México (por ejemplo, las cabeceras municipales o las unidades económicas del DENUE de INEGI):

- comparar las distancias euclidianas en cada proyección contra la distancia geodésica (fórmula del haversine) y reportar la distribución del error relativo por región del país;
- repetir un agrupamiento k-means o una búsqueda de k vecinos cercanos en cada proyección y medir cuántas asignaciones cambian;
- explicar por qué las áreas se deforman por $\sigma^2$ y qué implica para mapas de densidad.

## 5. Implementación y verificación

- Describir el método numérico usado para el problema de Dirichlet (diferencias finitas, elemento finito o método de soluciones fundamentales) y verificarlo con un caso de solución conocida.
- Verificar que las proyecciones propias coinciden con una biblioteca estándar (por ejemplo, PROJ vía pyproj en Python, o Proj.jl en Julia) para Mercator y Lambert.
- Se recomienda Julia o Python. El código se entrega como apéndice o enlace a un repositorio.

## 6. Herramienta para el cliente

- Una tabla comparativa de las tres proyecciones: distorsión $\delta$, error máximo de escala en porcentaje y error típico en distancias del conjunto de datos.
- Una función o script que convierta (latitud, longitud) a (x, y) con la proyección recomendada.

## 7. Limitaciones y alcance de validez

- Efecto del modelo esférico frente al elipsoide de referencia que usa INEGI.
- La conformidad preserva ángulos pero no áreas: cuándo el cliente debería usar una proyección equivalente o calcular directamente sobre la esfera.
- El teorema de Riemann garantiza la existencia de mapeos conformes pero no da fórmulas; en dimensión 3 o mayor solo hay mapeos conformes de Möbius (teorema de Liouville).
- Costo práctico: una proyección "a la medida" no viene en el software estándar.

## 8. Desarrollos recientes y otras aplicaciones

Una sección breve que muestre al cliente que el método está vigente: parametrización conforme de superficies en gráficos por computadora y neuroimagen, mapeos conformes discretos, y al menos una aplicación en ingeniería o física (por ejemplo, la transformación de Joukowski en aerodinámica o Schwarz–Christoffel en electrostática).

## 9. Conclusiones y recomendaciones al cliente

Responder las tres preguntas de la consultoría para este caso y dar una recomendación concreta: qué proyección usar y en qué tipo de análisis el error deja de ser despreciable.

## 10. Referencias

- Mínimo 10 referencias académicas: libros, artículos revisados por pares, capítulos de libro o recursos de instituciones científicas. Al menos 3 deben ser artículos revisados por pares.
- Formato IEEE, con consistencia obligatoria y DOI o enlace verificable cuando exista.
- Puntos de partida: J. Milnor, "A problem in cartography", *Amer. Math. Monthly*, vol. 76, no. 10, pp. 1101–1112, 1969; M. Bermejo y J. Otero, "Minimum conformal mapping distortion according to Chebyshev's principle: a case study over peninsular Spain", *J. Geodesy*, vol. 79, pp. 124–134, 2005; J. P. Snyder, *Map Projections: A Working Manual*, USGS, 1987; los textos del curso (Boas, Marsden y Hoffman, Ahlfors).

---

# Requisitos de formato

| Elemento | Requisito |
|---|---|
| Extensión | 2000–4000 palabras, sin contar ecuaciones, código ni referencias |
| Formato | PDF, LaTeX o Markdown |
| Lenguaje | Técnico, claro y coherente, pensado para el cliente |
| Figuras | Máximo 3, con pie de figura y fuente: se recomiendan la retícula bajo dos proyecciones, el mapa de σ sobre México y una de errores en distancias |
| Código | Apéndice o enlace a repositorio |
| Equipos | Hasta 3 integrantes; misma calificación para todos; no se cambian una vez definidos |

Nota: la evaluación y retroalimentación de las evidencias se realiza en **Canvas**.
