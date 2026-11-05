# stokes-solar

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23981088.svg)](https://doi.org/10.5281/zenodo.23981088)

Código e figuras do artigo **"Espectropolarimetria solar: dos parâmetros de Stokes ao campo magnético fotosférico, com Python"** (em preparação para a Revista Brasileira de Ensino de Física).

## Conteúdo

| Pasta | O que tem |
|---|---|
| `stokespy/` | Módulo com as funções do artigo: Stokes e Mueller, polarímetro, Zeeman e Milne-Eddington, campo fraco |
| `notebooks/` | Um notebook por bloco do artigo, já executados |
| `figuras/` | Figuras em PNG (300 dpi) e PDF |
| `tests/` | Testes de sanidade (`python -m pytest tests`) |
| `tools/` | Scripts que geram os notebooks |

| Notebook | Seções do artigo | Figuras do artigo | Material suplementar | Colab |
|---|---|---|---|---|
| `01_stokes_mueller_polarimetro` | 2, 3 | 1 | Figuras S1, S2, S3 | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ChristofherdoVP/stokes-solar/blob/main/notebooks/01_stokes_mueller_polarimetro.ipynb) |
| `02_zeeman_milne_eddington` | 4 | 2, 3 | Figuras S4, S5 | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ChristofherdoVP/stokes-solar/blob/main/notebooks/02_zeeman_milne_eddington.ipynb) |
| `03_campo_fraco` | 5.1, 5.2 | 4 | Tabela S1 | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ChristofherdoVP/stokes-solar/blob/main/notebooks/03_campo_fraco.ipynb) |
| `04_mancha_sintetica` | 5.3 | 5 | | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ChristofherdoVP/stokes-solar/blob/main/notebooks/04_mancha_sintetica.ipynb) |

Os arquivos em `figuras/` seguem a numeração do artigo (`fig1` a `fig5`) e do material suplementar (`figS1` a `figS5`). A pasta `suplementar/` traz o material suplementar em Word e PDF.

## Convenções

- $\delta = \phi_y - \phi_x$; $V > 0$ para rotação anti-horária vista pelo observador (IAU).
- Comprimentos de onda em Å, campo em G, velocidades em km/s.
- Síntese de Milne-Eddington pela forma matricial $\mathbf{I}(0) = (S_0\,\mathbb{1} + \mu S_1\,\mathbf{K}^{-1})\,\mathbf{e}_0$.

## Uso

```bash
pip install -r requirements.txt
python -m pytest tests
jupyter lab notebooks/
```

## Como citar

Veja `CITATION.cff`. DOI do Zenodo: [10.5281/zenodo.23981088](https://doi.org/10.5281/zenodo.23981088).

Licença MIT.
