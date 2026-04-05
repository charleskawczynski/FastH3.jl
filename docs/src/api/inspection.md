# Inspection

Functions for examining properties of H3 indexes.

Corresponds to the H3 C API [Inspection](https://h3geo.org/docs/api/inspection) section.

## Index Properties

```@docs
getResolution
getBaseCellNumber
isValidCell
isValidIndex
isPentagon
isResClassIII
maxFaceCount
getIcosahedronFaces
```

## String Conversion

```@docs
stringToH3
h3ToString
```

## Error Descriptions

```@docs
describeH3Error
```

## Types

```@docs
H3Index
H3_NULL
H3Error
LatLng
CellBoundary
CoordIJ
GeoLoop
GeoPolygon
ContainmentMode
```
