# Measurement

Functions for computing distances, areas, and edge lengths.

Corresponds to the H3 C API [Misc](https://h3geo.org/docs/api/misc) section.

## Great Circle Distance

```@docs
greatCircleDistanceRads
greatCircleDistanceKm
greatCircleDistanceM
```

## Unit Conversion

```@docs
degsToRads
radsToDegs
```

## Average Hexagon Area

```@docs
getHexagonAreaAvgKm2
getHexagonAreaAvgM2
```

## Average Edge Length

```@docs
getHexagonEdgeLengthAvgKm
getHexagonEdgeLengthAvgM
```

## Cell Count

```@docs
getNumCells
```

## Exact Cell Area

```@docs
FastH3.cellAreaRads2
FastH3.cellAreaKm2
FastH3.cellAreaM2
```

## Exact Edge Length

```@docs
FastH3.edgeLengthRads
FastH3.edgeLengthKm
FastH3.edgeLengthM
```
