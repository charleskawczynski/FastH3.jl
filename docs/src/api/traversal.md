# Traversal

Functions for finding cells near other cells and computing paths across the grid.

Corresponds to the H3 C API [Traversal](https://h3geo.org/docs/api/traversal) section.

## Grid Disk

```@docs
gridDisk
gridDiskDistances
gridDiskUnsafe
gridDiskDistancesUnsafe
gridRingUnsafe
maxGridDiskSize
```

## Grid Distance and Paths

```@docs
gridDistance
gridPathCells
gridPathCellsSize
```

## Local IJ Coordinates

```@docs
cellToLocalIj
localIjToCell
```

## Neighbor Detection

```@docs
areNeighborCells
```
