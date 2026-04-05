# Hierarchy

Functions for navigating the parent/child cell hierarchy and compact/uncompact operations.

Corresponds to the H3 C API [Hierarchy](https://h3geo.org/docs/api/hierarchy) section.

## Parent / Child Navigation

```@docs
cellToParent
cellToCenterChild
cellToChildren
cellToChildrenSize
cellToChildPos
childPosToCell
```

## Compact / Uncompact

```@docs
compactCells
uncompactCells
uncompactCellsSize
```

## Resolution 0 and Pentagons

```@docs
getRes0Cells
getPentagons
pentagonCount
```
