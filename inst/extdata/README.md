# Geographic reference for the metric MDS article

`europe-map.csv` is a cropped and simplified extract of Natural Earth's
`ne_50m_admin_0_countries` dataset, version 5.1.2. Natural Earth releases its
vector and raster map data into the public domain:

- <https://www.naturalearthdata.com/about/terms-of-use/>
- <https://github.com/nvkelso/natural-earth-vector/tree/v5.1.2>

The source GeoJSON has SHA-256 checksum
`3e458fc036ad0a66411f2c1e6cac49c5d7bfb81cb1123bc513b22511a2b7fdeb`.
It was cropped to longitude -12 through 31 and latitude 33 through 63, polygon
parts were extracted, and boundaries were simplified with a 0.025-degree
tolerance while preserving topology. The CSV retains polygon-ring groups plus
longitude and latitude; the metric MDS vignette performs its own projection.

`eurodist-geography.csv` uses the corresponding city coordinates from Natural
Earth's `ne_10m_populated_places_simple` dataset, version 5.1.2. Its source
GeoJSON has SHA-256 checksum
`fd3fa867a320cbd5c5b6bb5bc550afeec2939fb2cef688e508007282a55ac42f`.
Dataset spellings were matched to the 21 labels in `datasets::eurodist`.

Natural Earth does not include Hook of Holland in that city layer. Its
coordinate comes from GeoNames under the Creative Commons Attribution 4.0
license. The attribution and license apply to that row independently of the
package's BSD 2-clause license:

- <https://www.geonames.org/advanced-search.html?country=NL&q=&startRow=200>
- <https://www.geonames.org/>
- <https://creativecommons.org/licenses/by/4.0/>

The files were prepared on 2026-08-27. They are bundled so rendering the
article does not require network access or a spatial R package.
