# Base cell lookup tables and access functions
# Translated from H3 library baseCells.c / baseCells.h

const INVALID_BASE_CELL = 127
const MAX_FACE_COORD = 2
const INVALID_ROTATIONS = -1

struct BaseCellData
    homeFijk::FaceIJK
    isPentagon::Int32
    cwOffsetPent::Tuple{Int32, Int32}
end

struct BaseCellRotation
    baseCell::Int32
    ccwRot60::Int32
end

# Neighboring base cell ID in each IJK direction.
# Indexed as baseCellNeighbors[baseCell + 1][dir + 1].
# 127 (INVALID_BASE_CELL) indicates no neighbor in that direction.
const baseCellNeighbors = NTuple{7, Int32}[
    (0, 1, 5, 2, 4, 3, 8),                          # base cell 0
    (1, 7, 6, 9, 0, 3, 2),                           # base cell 1
    (2, 6, 10, 11, 0, 1, 5),                         # base cell 2
    (3, 13, 1, 7, 4, 12, 0),                         # base cell 3
    (4, INVALID_BASE_CELL, 15, 8, 3, 0, 12),         # base cell 4 (pentagon)
    (5, 2, 18, 10, 8, 0, 16),                        # base cell 5
    (6, 14, 11, 17, 1, 9, 2),                        # base cell 6
    (7, 21, 9, 19, 3, 13, 1),                        # base cell 7
    (8, 5, 22, 16, 4, 0, 15),                        # base cell 8
    (9, 19, 14, 20, 1, 7, 6),                        # base cell 9
    (10, 11, 24, 23, 5, 2, 18),                      # base cell 10
    (11, 17, 23, 25, 2, 6, 10),                      # base cell 11
    (12, 28, 13, 26, 4, 15, 3),                      # base cell 12
    (13, 26, 21, 29, 3, 12, 7),                      # base cell 13
    (14, INVALID_BASE_CELL, 17, 27, 9, 20, 6),       # base cell 14 (pentagon)
    (15, 22, 28, 31, 4, 8, 12),                      # base cell 15
    (16, 18, 33, 30, 8, 5, 22),                      # base cell 16
    (17, 11, 14, 6, 35, 25, 27),                     # base cell 17
    (18, 24, 30, 32, 5, 10, 16),                     # base cell 18
    (19, 34, 20, 36, 7, 21, 9),                      # base cell 19
    (20, 14, 19, 9, 40, 27, 36),                     # base cell 20
    (21, 38, 19, 34, 13, 29, 7),                     # base cell 21
    (22, 16, 41, 33, 15, 8, 31),                     # base cell 22
    (23, 24, 11, 10, 39, 37, 25),                    # base cell 23
    (24, INVALID_BASE_CELL, 32, 37, 10, 23, 18),     # base cell 24 (pentagon)
    (25, 23, 17, 11, 45, 39, 35),                    # base cell 25
    (26, 42, 29, 43, 12, 28, 13),                    # base cell 26
    (27, 40, 35, 46, 14, 20, 17),                    # base cell 27
    (28, 31, 42, 44, 12, 15, 26),                    # base cell 28
    (29, 43, 38, 47, 13, 26, 21),                    # base cell 29
    (30, 32, 48, 50, 16, 18, 33),                    # base cell 30
    (31, 41, 44, 53, 15, 22, 28),                    # base cell 31
    (32, 30, 24, 18, 52, 50, 37),                    # base cell 32
    (33, 30, 49, 48, 22, 16, 41),                    # base cell 33
    (34, 19, 38, 21, 54, 36, 51),                    # base cell 34
    (35, 46, 45, 56, 17, 27, 25),                    # base cell 35
    (36, 20, 34, 19, 55, 40, 54),                    # base cell 36
    (37, 39, 52, 57, 24, 23, 32),                    # base cell 37
    (38, INVALID_BASE_CELL, 34, 51, 29, 47, 21),     # base cell 38 (pentagon)
    (39, 37, 25, 23, 59, 57, 45),                    # base cell 39
    (40, 27, 36, 20, 60, 46, 55),                    # base cell 40
    (41, 49, 53, 61, 22, 33, 31),                    # base cell 41
    (42, 58, 43, 62, 28, 44, 26),                    # base cell 42
    (43, 62, 47, 64, 26, 42, 29),                    # base cell 43
    (44, 53, 58, 65, 28, 31, 42),                    # base cell 44
    (45, 39, 35, 25, 63, 59, 56),                    # base cell 45
    (46, 60, 56, 68, 27, 40, 35),                    # base cell 46
    (47, 38, 43, 29, 69, 51, 64),                    # base cell 47
    (48, 49, 30, 33, 67, 66, 50),                    # base cell 48
    (49, INVALID_BASE_CELL, 61, 66, 33, 48, 41),     # base cell 49 (pentagon)
    (50, 48, 32, 30, 70, 67, 52),                    # base cell 50
    (51, 69, 54, 71, 38, 47, 34),                    # base cell 51
    (52, 57, 70, 74, 32, 37, 50),                    # base cell 52
    (53, 61, 65, 75, 31, 41, 44),                    # base cell 53
    (54, 71, 55, 73, 34, 51, 36),                    # base cell 54
    (55, 40, 54, 36, 72, 60, 73),                    # base cell 55
    (56, 68, 63, 77, 35, 46, 45),                    # base cell 56
    (57, 59, 74, 78, 37, 39, 52),                    # base cell 57
    (58, INVALID_BASE_CELL, 62, 76, 44, 65, 42),     # base cell 58 (pentagon)
    (59, 63, 78, 79, 39, 45, 57),                    # base cell 59
    (60, 72, 68, 80, 40, 55, 46),                    # base cell 60
    (61, 53, 49, 41, 81, 75, 66),                    # base cell 61
    (62, 43, 58, 42, 82, 64, 76),                    # base cell 62
    (63, INVALID_BASE_CELL, 56, 45, 79, 59, 77),     # base cell 63 (pentagon)
    (64, 47, 62, 43, 84, 69, 82),                    # base cell 64
    (65, 58, 53, 44, 86, 76, 75),                    # base cell 65
    (66, 67, 81, 85, 49, 48, 61),                    # base cell 66
    (67, 66, 50, 48, 87, 85, 70),                    # base cell 67
    (68, 56, 60, 46, 90, 77, 80),                    # base cell 68
    (69, 51, 64, 47, 89, 71, 84),                    # base cell 69
    (70, 67, 52, 50, 83, 87, 74),                    # base cell 70
    (71, 89, 73, 91, 51, 69, 54),                    # base cell 71
    (72, INVALID_BASE_CELL, 73, 55, 80, 60, 88),     # base cell 72 (pentagon)
    (73, 91, 72, 88, 54, 71, 55),                    # base cell 73
    (74, 78, 83, 92, 52, 57, 70),                    # base cell 74
    (75, 65, 61, 53, 94, 86, 81),                    # base cell 75
    (76, 86, 82, 96, 58, 65, 62),                    # base cell 76
    (77, 63, 68, 56, 93, 79, 90),                    # base cell 77
    (78, 74, 59, 57, 95, 92, 79),                    # base cell 78
    (79, 78, 63, 59, 93, 95, 77),                    # base cell 79
    (80, 68, 72, 60, 99, 90, 88),                    # base cell 80
    (81, 85, 94, 101, 61, 66, 75),                   # base cell 81
    (82, 96, 84, 98, 62, 76, 64),                    # base cell 82
    (83, INVALID_BASE_CELL, 74, 70, 100, 87, 92),    # base cell 83 (pentagon)
    (84, 69, 82, 64, 97, 89, 98),                    # base cell 84
    (85, 87, 101, 102, 66, 67, 81),                  # base cell 85
    (86, 76, 75, 65, 104, 96, 94),                   # base cell 86
    (87, 83, 102, 100, 67, 70, 85),                  # base cell 87
    (88, 72, 91, 73, 99, 80, 105),                   # base cell 88
    (89, 97, 91, 103, 69, 84, 71),                   # base cell 89
    (90, 77, 80, 68, 106, 93, 99),                   # base cell 90
    (91, 73, 89, 71, 105, 88, 103),                  # base cell 91
    (92, 83, 78, 74, 108, 100, 95),                  # base cell 92
    (93, 79, 90, 77, 109, 95, 106),                  # base cell 93
    (94, 86, 81, 75, 107, 104, 101),                 # base cell 94
    (95, 92, 79, 78, 109, 108, 93),                  # base cell 95
    (96, 104, 98, 110, 76, 86, 82),                  # base cell 96
    (97, INVALID_BASE_CELL, 98, 84, 103, 89, 111),   # base cell 97 (pentagon)
    (98, 110, 97, 111, 82, 96, 84),                  # base cell 98
    (99, 80, 105, 88, 106, 90, 113),                 # base cell 99
    (100, 102, 83, 87, 108, 114, 92),                # base cell 100
    (101, 102, 107, 112, 81, 85, 94),                # base cell 101
    (102, 101, 87, 85, 114, 112, 100),               # base cell 102
    (103, 91, 97, 89, 116, 105, 111),                # base cell 103
    (104, 107, 110, 115, 86, 94, 96),                # base cell 104
    (105, 88, 103, 91, 113, 99, 116),                # base cell 105
    (106, 93, 99, 90, 117, 109, 113),                # base cell 106
    (107, INVALID_BASE_CELL, 101, 94, 115, 104, 112),# base cell 107 (pentagon)
    (108, 100, 95, 92, 118, 114, 109),               # base cell 108
    (109, 108, 93, 95, 117, 118, 106),               # base cell 109
    (110, 98, 104, 96, 119, 111, 115),               # base cell 110
    (111, 97, 110, 98, 116, 103, 119),               # base cell 111
    (112, 107, 102, 101, 120, 115, 114),             # base cell 112
    (113, 99, 116, 105, 117, 106, 121),              # base cell 113
    (114, 112, 100, 102, 118, 120, 108),             # base cell 114
    (115, 110, 107, 104, 120, 119, 112),             # base cell 115
    (116, 103, 119, 111, 113, 105, 121),             # base cell 116
    (117, INVALID_BASE_CELL, 109, 118, 113, 121, 106),# base cell 117 (pentagon)
    (118, 120, 108, 114, 117, 121, 109),             # base cell 118
    (119, 111, 115, 110, 121, 116, 120),             # base cell 119
    (120, 115, 114, 112, 121, 119, 118),             # base cell 120
    (121, 116, 120, 119, 117, 113, 118),             # base cell 121
]

# Neighboring base cell rotations in each IJK direction.
# For each base cell, for each direction, the number of 60 degree CCW rotations
# to the coordinate system of the neighbor. -1 indicates no neighbor.
# Indexed as baseCellNeighbor60CCWRots[baseCell + 1][dir + 1].
const baseCellNeighbor60CCWRots = NTuple{7, Int32}[
    (0, 5, 0, 0, 1, 5, 1),    # base cell 0
    (0, 0, 1, 0, 1, 0, 1),    # base cell 1
    (0, 0, 0, 0, 0, 5, 0),    # base cell 2
    (0, 5, 0, 0, 2, 5, 1),    # base cell 3
    (0, -1, 1, 0, 3, 4, 2),   # base cell 4 (pentagon)
    (0, 0, 1, 0, 1, 0, 1),    # base cell 5
    (0, 0, 0, 3, 5, 5, 0),    # base cell 6
    (0, 0, 0, 0, 0, 5, 0),    # base cell 7
    (0, 5, 0, 0, 0, 5, 1),    # base cell 8
    (0, 0, 1, 3, 0, 0, 1),    # base cell 9
    (0, 0, 1, 3, 0, 0, 1),    # base cell 10
    (0, 3, 3, 3, 0, 0, 0),    # base cell 11
    (0, 5, 0, 0, 3, 5, 1),    # base cell 12
    (0, 0, 1, 0, 1, 0, 1),    # base cell 13
    (0, -1, 3, 0, 5, 2, 0),   # base cell 14 (pentagon)
    (0, 5, 0, 0, 4, 5, 1),    # base cell 15
    (0, 0, 0, 0, 0, 5, 0),    # base cell 16
    (0, 3, 3, 3, 3, 0, 3),    # base cell 17
    (0, 0, 0, 3, 5, 5, 0),    # base cell 18
    (0, 3, 3, 3, 0, 0, 0),    # base cell 19
    (0, 3, 3, 3, 0, 3, 0),    # base cell 20
    (0, 0, 0, 3, 5, 5, 0),    # base cell 21
    (0, 0, 1, 0, 1, 0, 1),    # base cell 22
    (0, 3, 3, 3, 0, 3, 0),    # base cell 23
    (0, -1, 3, 0, 5, 2, 0),   # base cell 24 (pentagon)
    (0, 0, 0, 3, 0, 0, 3),    # base cell 25
    (0, 0, 0, 0, 0, 5, 0),    # base cell 26
    (0, 3, 0, 0, 0, 3, 3),    # base cell 27
    (0, 0, 1, 0, 1, 0, 1),    # base cell 28
    (0, 0, 1, 3, 0, 0, 1),    # base cell 29
    (0, 3, 3, 3, 0, 0, 0),    # base cell 30
    (0, 0, 0, 0, 0, 5, 0),    # base cell 31
    (0, 3, 3, 3, 3, 0, 3),    # base cell 32
    (0, 0, 1, 3, 0, 0, 1),    # base cell 33
    (0, 3, 3, 3, 3, 0, 3),    # base cell 34
    (0, 0, 3, 0, 3, 0, 3),    # base cell 35
    (0, 0, 0, 3, 0, 0, 3),    # base cell 36
    (0, 3, 0, 0, 0, 3, 3),    # base cell 37
    (0, -1, 3, 0, 5, 2, 0),   # base cell 38 (pentagon)
    (0, 3, 0, 0, 3, 3, 0),    # base cell 39
    (0, 3, 0, 0, 3, 3, 0),    # base cell 40
    (0, 0, 0, 3, 5, 5, 0),    # base cell 41
    (0, 0, 0, 3, 5, 5, 0),    # base cell 42
    (0, 3, 3, 3, 0, 0, 0),    # base cell 43
    (0, 0, 1, 3, 0, 0, 1),    # base cell 44
    (0, 0, 3, 0, 0, 3, 3),    # base cell 45
    (0, 0, 0, 3, 0, 3, 0),    # base cell 46
    (0, 3, 3, 3, 0, 3, 0),    # base cell 47
    (0, 3, 3, 3, 0, 3, 0),    # base cell 48
    (0, -1, 3, 0, 5, 2, 0),   # base cell 49 (pentagon)
    (0, 0, 0, 3, 0, 0, 3),    # base cell 50
    (0, 3, 0, 0, 0, 3, 3),    # base cell 51
    (0, 0, 3, 0, 3, 0, 3),    # base cell 52
    (0, 3, 3, 3, 0, 0, 0),    # base cell 53
    (0, 0, 3, 0, 3, 0, 3),    # base cell 54
    (0, 0, 3, 0, 0, 3, 3),    # base cell 55
    (0, 3, 3, 3, 0, 0, 3),    # base cell 56
    (0, 0, 0, 3, 0, 3, 0),    # base cell 57
    (0, -1, 3, 0, 5, 2, 0),   # base cell 58 (pentagon)
    (0, 3, 3, 3, 3, 3, 0),    # base cell 59
    (0, 3, 3, 3, 3, 3, 0),    # base cell 60
    (0, 3, 3, 3, 3, 0, 3),    # base cell 61
    (0, 3, 3, 3, 3, 0, 3),    # base cell 62
    (0, -1, 3, 0, 5, 2, 0),   # base cell 63 (pentagon)
    (0, 0, 0, 3, 0, 0, 3),    # base cell 64
    (0, 3, 3, 3, 0, 3, 0),    # base cell 65
    (0, 3, 0, 0, 0, 3, 3),    # base cell 66
    (0, 3, 0, 0, 3, 3, 0),    # base cell 67
    (0, 3, 3, 3, 0, 0, 0),    # base cell 68
    (0, 3, 0, 0, 3, 3, 0),    # base cell 69
    (0, 0, 3, 0, 0, 3, 3),    # base cell 70
    (0, 0, 0, 3, 0, 3, 0),    # base cell 71
    (0, -1, 3, 0, 5, 2, 0),   # base cell 72 (pentagon)
    (0, 3, 3, 3, 0, 0, 3),    # base cell 73
    (0, 3, 3, 3, 0, 0, 3),    # base cell 74
    (0, 0, 0, 3, 0, 0, 3),    # base cell 75
    (0, 3, 0, 0, 0, 3, 3),    # base cell 76
    (0, 0, 0, 3, 0, 5, 0),    # base cell 77
    (0, 3, 3, 3, 0, 0, 0),    # base cell 78
    (0, 0, 1, 3, 1, 0, 1),    # base cell 79
    (0, 0, 1, 3, 1, 0, 1),    # base cell 80
    (0, 0, 3, 0, 3, 0, 3),    # base cell 81
    (0, 0, 3, 0, 3, 0, 3),    # base cell 82
    (0, -1, 3, 0, 5, 2, 0),   # base cell 83 (pentagon)
    (0, 0, 3, 0, 0, 3, 3),    # base cell 84
    (0, 0, 0, 3, 0, 3, 0),    # base cell 85
    (0, 3, 0, 0, 3, 3, 0),    # base cell 86
    (0, 3, 3, 3, 3, 3, 0),    # base cell 87
    (0, 0, 0, 3, 0, 5, 0),    # base cell 88
    (0, 3, 3, 3, 3, 3, 0),    # base cell 89
    (0, 0, 0, 0, 0, 0, 1),    # base cell 90
    (0, 3, 3, 3, 0, 0, 0),    # base cell 91
    (0, 0, 0, 3, 0, 5, 0),    # base cell 92
    (0, 5, 0, 0, 5, 5, 0),    # base cell 93
    (0, 0, 3, 0, 0, 3, 3),    # base cell 94
    (0, 0, 0, 0, 0, 0, 1),    # base cell 95
    (0, 0, 0, 3, 0, 3, 0),    # base cell 96
    (0, -1, 3, 0, 5, 2, 0),   # base cell 97 (pentagon)
    (0, 3, 3, 3, 0, 0, 3),    # base cell 98
    (0, 5, 0, 0, 5, 5, 0),    # base cell 99
    (0, 0, 1, 3, 1, 0, 1),    # base cell 100
    (0, 3, 3, 3, 0, 0, 3),    # base cell 101
    (0, 3, 3, 3, 0, 0, 0),    # base cell 102
    (0, 0, 1, 3, 1, 0, 1),    # base cell 103
    (0, 3, 3, 3, 3, 3, 0),    # base cell 104
    (0, 0, 0, 0, 0, 0, 1),    # base cell 105
    (0, 0, 1, 0, 3, 5, 1),    # base cell 106
    (0, -1, 3, 0, 5, 2, 0),   # base cell 107 (pentagon)
    (0, 5, 0, 0, 5, 5, 0),    # base cell 108
    (0, 0, 1, 0, 4, 5, 1),    # base cell 109
    (0, 3, 3, 3, 0, 0, 0),    # base cell 110
    (0, 0, 0, 3, 0, 5, 0),    # base cell 111
    (0, 0, 0, 3, 0, 5, 0),    # base cell 112
    (0, 0, 1, 0, 2, 5, 1),    # base cell 113
    (0, 0, 0, 0, 0, 0, 1),    # base cell 114
    (0, 0, 1, 3, 1, 0, 1),    # base cell 115
    (0, 5, 0, 0, 5, 5, 0),    # base cell 116
    (0, -1, 1, 0, 3, 4, 2),   # base cell 117 (pentagon)
    (0, 0, 1, 0, 0, 5, 1),    # base cell 118
    (0, 0, 0, 0, 0, 0, 1),    # base cell 119
    (0, 5, 0, 0, 5, 5, 0),    # base cell 120
    (0, 0, 1, 0, 1, 5, 1),    # base cell 121
]

# Resolution 0 base cell lookup table for each face.
# Indexed as faceIjkBaseCells[face+1, i+1, j+1, k+1].
# Each entry is a BaseCellRotation(baseCell, ccwRot60).
const faceIjkBaseCells = let
    R(b, r) = BaseCellRotation(Int32(b), Int32(r))
    a = Array{BaseCellRotation}(undef, 20, 3, 3, 3)

    function setface!(a, f, data)
        idx = 1
        for i in 1:3, j in 1:3, k in 1:3
            a[f, i, j, k] = data[idx]
            idx += 1
        end
    end

    # face 0
    setface!(a, 1, [
        R(16,0), R(18,0), R(24,0),  R(33,0), R(30,0), R(32,3),  R(49,1), R(48,3), R(50,3),
        R(8,0),  R(5,5),  R(10,5),  R(22,0), R(16,0), R(18,0),  R(41,1), R(33,0), R(30,0),
        R(4,0),  R(0,5),  R(2,5),   R(15,1), R(8,0),  R(5,5),   R(31,1), R(22,0), R(16,0),
    ])
    # face 1
    setface!(a, 2, [
        R(2,0),  R(6,0),  R(14,0),  R(10,0), R(11,0), R(17,3),  R(24,1), R(23,3), R(25,3),
        R(0,0),  R(1,5),  R(9,5),   R(5,0),  R(2,0),  R(6,0),   R(18,1), R(10,0), R(11,0),
        R(4,1),  R(3,5),  R(7,5),   R(8,1),  R(0,0),  R(1,5),   R(16,1), R(5,0),  R(2,0),
    ])
    # face 2
    setface!(a, 3, [
        R(7,0),  R(21,0), R(38,0),  R(9,0),  R(19,0), R(34,3),  R(14,1), R(20,3), R(36,3),
        R(3,0),  R(13,5), R(29,5),  R(1,0),  R(7,0),  R(21,0),  R(6,1),  R(9,0),  R(19,0),
        R(4,2),  R(12,5), R(26,5),  R(0,1),  R(3,0),  R(13,5),  R(2,1),  R(1,0),  R(7,0),
    ])
    # face 3
    setface!(a, 4, [
        R(26,0), R(42,0), R(58,0),  R(29,0), R(43,0), R(62,3),  R(38,1), R(47,3), R(64,3),
        R(12,0), R(28,5), R(44,5),  R(13,0), R(26,0), R(42,0),  R(21,1), R(29,0), R(43,0),
        R(4,3),  R(15,5), R(31,5),  R(3,1),  R(12,0), R(28,5),  R(7,1),  R(13,0), R(26,0),
    ])
    # face 4
    setface!(a, 5, [
        R(31,0), R(41,0), R(49,0),  R(44,0), R(53,0), R(61,3),  R(58,1), R(65,3), R(75,3),
        R(15,0), R(22,5), R(33,5),  R(28,0), R(31,0), R(41,0),  R(42,1), R(44,0), R(53,0),
        R(4,4),  R(8,5),  R(16,5),  R(12,1), R(15,0), R(22,5),  R(26,1), R(28,0), R(31,0),
    ])
    # face 5
    setface!(a, 6, [
        R(50,0), R(48,0), R(49,3),  R(32,0), R(30,3), R(33,3),  R(24,3), R(18,3), R(16,3),
        R(70,0), R(67,0), R(66,3),  R(52,3), R(50,0), R(48,0),  R(37,3), R(32,0), R(30,3),
        R(83,0), R(87,3), R(85,3),  R(74,3), R(70,0), R(67,0),  R(57,1), R(52,3), R(50,0),
    ])
    # face 6
    setface!(a, 7, [
        R(25,0), R(23,0), R(24,3),  R(17,0), R(11,3), R(10,3),  R(14,3), R(6,3),  R(2,3),
        R(45,0), R(39,0), R(37,3),  R(35,3), R(25,0), R(23,0),  R(27,3), R(17,0), R(11,3),
        R(63,0), R(59,3), R(57,3),  R(56,3), R(45,0), R(39,0),  R(46,3), R(35,3), R(25,0),
    ])
    # face 7
    setface!(a, 8, [
        R(36,0), R(20,0), R(14,3),  R(34,0), R(19,3), R(9,3),   R(38,3), R(21,3), R(7,3),
        R(55,0), R(40,0), R(27,3),  R(54,3), R(36,0), R(20,0),  R(51,3), R(34,0), R(19,3),
        R(72,0), R(60,3), R(46,3),  R(73,3), R(55,0), R(40,0),  R(71,3), R(54,3), R(36,0),
    ])
    # face 8
    setface!(a, 9, [
        R(64,0), R(47,0), R(38,3),  R(62,0), R(43,3), R(29,3),  R(58,3), R(42,3), R(26,3),
        R(84,0), R(69,0), R(51,3),  R(82,3), R(64,0), R(47,0),  R(76,3), R(62,0), R(43,3),
        R(97,0), R(89,3), R(71,3),  R(98,3), R(84,0), R(69,0),  R(96,3), R(82,3), R(64,0),
    ])
    # face 9
    setface!(a, 10, [
        R(75,0),  R(65,0),  R(58,3),  R(61,0),  R(53,3),  R(44,3),  R(49,3),  R(41,3),  R(31,3),
        R(94,0),  R(86,0),  R(76,3),  R(81,3),  R(75,0),  R(65,0),  R(66,3),  R(61,0),  R(53,3),
        R(107,0), R(104,3), R(96,3),  R(101,3), R(94,0),  R(86,0),  R(85,3),  R(81,3),  R(75,0),
    ])
    # face 10
    setface!(a, 11, [
        R(57,0),  R(59,0),  R(63,3),  R(74,0),  R(78,3),  R(79,3),  R(83,3),  R(92,3),  R(95,3),
        R(37,0),  R(39,3),  R(45,3),  R(52,0),  R(57,0),  R(59,0),  R(70,3),  R(74,0),  R(78,3),
        R(24,0),  R(23,3),  R(25,3),  R(32,3),  R(37,0),  R(39,3),  R(50,3),  R(52,0),  R(57,0),
    ])
    # face 11
    setface!(a, 12, [
        R(46,0), R(60,0), R(72,3),  R(56,0), R(68,3), R(80,3),  R(63,3), R(77,3), R(90,3),
        R(27,0), R(40,3), R(55,3),  R(35,0), R(46,0), R(60,0),  R(45,3), R(56,0), R(68,3),
        R(14,0), R(20,3), R(36,3),  R(17,3), R(27,0), R(40,3),  R(25,3), R(35,0), R(46,0),
    ])
    # face 12
    setface!(a, 13, [
        R(71,0),  R(89,0),  R(97,3),   R(73,0),  R(91,3),  R(103,3),  R(72,3),  R(88,3),  R(105,3),
        R(51,0),  R(69,3),  R(84,3),   R(54,0),  R(71,0),  R(89,0),   R(55,3),  R(73,0),  R(91,3),
        R(38,0),  R(47,3),  R(64,3),   R(34,3),  R(51,0),  R(69,3),   R(36,3),  R(54,0),  R(71,0),
    ])
    # face 13
    setface!(a, 14, [
        R(96,0),  R(104,0), R(107,3),  R(98,0),  R(110,3), R(115,3),  R(97,3),  R(111,3), R(119,3),
        R(76,0),  R(86,3),  R(94,3),   R(82,0),  R(96,0),  R(104,0),  R(84,3),  R(98,0),  R(110,3),
        R(58,0),  R(65,3),  R(75,3),   R(62,3),  R(76,0),  R(86,3),   R(64,3),  R(82,0),  R(96,0),
    ])
    # face 14
    setface!(a, 15, [
        R(85,0),  R(87,0),  R(83,3),   R(101,0), R(102,3), R(100,3),  R(107,3), R(112,3), R(114,3),
        R(66,0),  R(67,3),  R(70,3),   R(81,0),  R(85,0),  R(87,0),   R(94,3),  R(101,0), R(102,3),
        R(49,0),  R(48,3),  R(50,3),   R(61,3),  R(66,0),  R(67,3),   R(75,3),  R(81,0),  R(85,0),
    ])
    # face 15
    setface!(a, 16, [
        R(95,0),  R(92,0),  R(83,0),   R(79,0),  R(78,0),  R(74,3),   R(63,1),  R(59,3),  R(57,3),
        R(109,0), R(108,0), R(100,5),  R(93,1),  R(95,0),  R(92,0),   R(77,1),  R(79,0),  R(78,0),
        R(117,4), R(118,5), R(114,5),  R(106,1), R(109,0), R(108,0),  R(90,1),  R(93,1),  R(95,0),
    ])
    # face 16
    setface!(a, 17, [
        R(90,0),  R(77,0),  R(63,0),   R(80,0),  R(68,0),  R(56,3),   R(72,1),  R(60,3),  R(46,3),
        R(106,0), R(93,0),  R(79,5),   R(99,1),  R(90,0),  R(77,0),   R(88,1),  R(80,0),  R(68,0),
        R(117,3), R(109,5), R(95,5),   R(113,1), R(106,0), R(93,0),   R(105,1), R(99,1),  R(90,0),
    ])
    # face 17
    setface!(a, 18, [
        R(105,0), R(88,0),  R(72,0),   R(103,0), R(91,0),  R(73,3),   R(97,1),  R(89,3),  R(71,3),
        R(113,0), R(99,0),  R(80,5),   R(116,1), R(105,0), R(88,0),   R(111,1), R(103,0), R(91,0),
        R(117,2), R(106,5), R(90,5),   R(121,1), R(113,0), R(99,0),   R(119,1), R(116,1), R(105,0),
    ])
    # face 18
    setface!(a, 19, [
        R(119,0), R(111,0), R(97,0),   R(115,0), R(110,0), R(98,3),   R(107,1), R(104,3), R(96,3),
        R(121,0), R(116,0), R(103,5),  R(120,1), R(119,0), R(111,0),  R(112,1), R(115,0), R(110,0),
        R(117,1), R(113,5), R(105,5),  R(118,1), R(121,0), R(116,0),  R(114,1), R(120,1), R(119,0),
    ])
    # face 19
    setface!(a, 20, [
        R(114,0), R(112,0), R(107,0),  R(100,0), R(102,0), R(101,3),  R(83,1),  R(87,3),  R(85,3),
        R(118,0), R(120,0), R(115,5),  R(108,1), R(114,0), R(112,0),  R(92,1),  R(100,0), R(102,0),
        R(117,0), R(121,5), R(119,5),  R(109,1), R(118,0), R(120,0),  R(95,1),  R(108,1), R(114,0),
    ])

    a
end

# Resolution 0 base cell data table.
# Indexed as baseCellData[baseCell + 1].
const baseCellData = BaseCellData[
    BaseCellData(FaceIJK(1,  CoordIJK(1, 0, 0)), 0, (0, 0)),      # base cell 0
    BaseCellData(FaceIJK(2,  CoordIJK(1, 1, 0)), 0, (0, 0)),      # base cell 1
    BaseCellData(FaceIJK(1,  CoordIJK(0, 0, 0)), 0, (0, 0)),      # base cell 2
    BaseCellData(FaceIJK(2,  CoordIJK(1, 0, 0)), 0, (0, 0)),      # base cell 3
    BaseCellData(FaceIJK(0,  CoordIJK(2, 0, 0)), 1, (-1, -1)),    # base cell 4
    BaseCellData(FaceIJK(1,  CoordIJK(1, 1, 0)), 0, (0, 0)),      # base cell 5
    BaseCellData(FaceIJK(1,  CoordIJK(0, 0, 1)), 0, (0, 0)),      # base cell 6
    BaseCellData(FaceIJK(2,  CoordIJK(0, 0, 0)), 0, (0, 0)),      # base cell 7
    BaseCellData(FaceIJK(0,  CoordIJK(1, 0, 0)), 0, (0, 0)),      # base cell 8
    BaseCellData(FaceIJK(2,  CoordIJK(0, 1, 0)), 0, (0, 0)),      # base cell 9
    BaseCellData(FaceIJK(1,  CoordIJK(0, 1, 0)), 0, (0, 0)),      # base cell 10
    BaseCellData(FaceIJK(1,  CoordIJK(0, 1, 1)), 0, (0, 0)),      # base cell 11
    BaseCellData(FaceIJK(3,  CoordIJK(1, 0, 0)), 0, (0, 0)),      # base cell 12
    BaseCellData(FaceIJK(3,  CoordIJK(1, 1, 0)), 0, (0, 0)),      # base cell 13
    BaseCellData(FaceIJK(11, CoordIJK(2, 0, 0)), 1, (2, 6)),      # base cell 14
    BaseCellData(FaceIJK(4,  CoordIJK(1, 0, 0)), 0, (0, 0)),      # base cell 15
    BaseCellData(FaceIJK(0,  CoordIJK(0, 0, 0)), 0, (0, 0)),      # base cell 16
    BaseCellData(FaceIJK(6,  CoordIJK(0, 1, 0)), 0, (0, 0)),      # base cell 17
    BaseCellData(FaceIJK(0,  CoordIJK(0, 0, 1)), 0, (0, 0)),      # base cell 18
    BaseCellData(FaceIJK(2,  CoordIJK(0, 1, 1)), 0, (0, 0)),      # base cell 19
    BaseCellData(FaceIJK(7,  CoordIJK(0, 0, 1)), 0, (0, 0)),      # base cell 20
    BaseCellData(FaceIJK(2,  CoordIJK(0, 0, 1)), 0, (0, 0)),      # base cell 21
    BaseCellData(FaceIJK(0,  CoordIJK(1, 1, 0)), 0, (0, 0)),      # base cell 22
    BaseCellData(FaceIJK(6,  CoordIJK(0, 0, 1)), 0, (0, 0)),      # base cell 23
    BaseCellData(FaceIJK(10, CoordIJK(2, 0, 0)), 1, (1, 5)),      # base cell 24
    BaseCellData(FaceIJK(6,  CoordIJK(0, 0, 0)), 0, (0, 0)),      # base cell 25
    BaseCellData(FaceIJK(3,  CoordIJK(0, 0, 0)), 0, (0, 0)),      # base cell 26
    BaseCellData(FaceIJK(11, CoordIJK(1, 0, 0)), 0, (0, 0)),      # base cell 27
    BaseCellData(FaceIJK(4,  CoordIJK(1, 1, 0)), 0, (0, 0)),      # base cell 28
    BaseCellData(FaceIJK(3,  CoordIJK(0, 1, 0)), 0, (0, 0)),      # base cell 29
    BaseCellData(FaceIJK(0,  CoordIJK(0, 1, 1)), 0, (0, 0)),      # base cell 30
    BaseCellData(FaceIJK(4,  CoordIJK(0, 0, 0)), 0, (0, 0)),      # base cell 31
    BaseCellData(FaceIJK(5,  CoordIJK(0, 1, 0)), 0, (0, 0)),      # base cell 32
    BaseCellData(FaceIJK(0,  CoordIJK(0, 1, 0)), 0, (0, 0)),      # base cell 33
    BaseCellData(FaceIJK(7,  CoordIJK(0, 1, 0)), 0, (0, 0)),      # base cell 34
    BaseCellData(FaceIJK(11, CoordIJK(1, 1, 0)), 0, (0, 0)),      # base cell 35
    BaseCellData(FaceIJK(7,  CoordIJK(0, 0, 0)), 0, (0, 0)),      # base cell 36
    BaseCellData(FaceIJK(10, CoordIJK(1, 0, 0)), 0, (0, 0)),      # base cell 37
    BaseCellData(FaceIJK(12, CoordIJK(2, 0, 0)), 1, (3, 7)),      # base cell 38
    BaseCellData(FaceIJK(6,  CoordIJK(1, 0, 1)), 0, (0, 0)),      # base cell 39
    BaseCellData(FaceIJK(7,  CoordIJK(1, 0, 1)), 0, (0, 0)),      # base cell 40
    BaseCellData(FaceIJK(4,  CoordIJK(0, 0, 1)), 0, (0, 0)),      # base cell 41
    BaseCellData(FaceIJK(3,  CoordIJK(0, 0, 1)), 0, (0, 0)),      # base cell 42
    BaseCellData(FaceIJK(3,  CoordIJK(0, 1, 1)), 0, (0, 0)),      # base cell 43
    BaseCellData(FaceIJK(4,  CoordIJK(0, 1, 0)), 0, (0, 0)),      # base cell 44
    BaseCellData(FaceIJK(6,  CoordIJK(1, 0, 0)), 0, (0, 0)),      # base cell 45
    BaseCellData(FaceIJK(11, CoordIJK(0, 0, 0)), 0, (0, 0)),      # base cell 46
    BaseCellData(FaceIJK(8,  CoordIJK(0, 0, 1)), 0, (0, 0)),      # base cell 47
    BaseCellData(FaceIJK(5,  CoordIJK(0, 0, 1)), 0, (0, 0)),      # base cell 48
    BaseCellData(FaceIJK(14, CoordIJK(2, 0, 0)), 1, (0, 9)),      # base cell 49
    BaseCellData(FaceIJK(5,  CoordIJK(0, 0, 0)), 0, (0, 0)),      # base cell 50
    BaseCellData(FaceIJK(12, CoordIJK(1, 0, 0)), 0, (0, 0)),      # base cell 51
    BaseCellData(FaceIJK(10, CoordIJK(1, 1, 0)), 0, (0, 0)),      # base cell 52
    BaseCellData(FaceIJK(4,  CoordIJK(0, 1, 1)), 0, (0, 0)),      # base cell 53
    BaseCellData(FaceIJK(12, CoordIJK(1, 1, 0)), 0, (0, 0)),      # base cell 54
    BaseCellData(FaceIJK(7,  CoordIJK(1, 0, 0)), 0, (0, 0)),      # base cell 55
    BaseCellData(FaceIJK(11, CoordIJK(0, 1, 0)), 0, (0, 0)),      # base cell 56
    BaseCellData(FaceIJK(10, CoordIJK(0, 0, 0)), 0, (0, 0)),      # base cell 57
    BaseCellData(FaceIJK(13, CoordIJK(2, 0, 0)), 1, (4, 8)),      # base cell 58
    BaseCellData(FaceIJK(10, CoordIJK(0, 0, 1)), 0, (0, 0)),      # base cell 59
    BaseCellData(FaceIJK(11, CoordIJK(0, 0, 1)), 0, (0, 0)),      # base cell 60
    BaseCellData(FaceIJK(9,  CoordIJK(0, 1, 0)), 0, (0, 0)),      # base cell 61
    BaseCellData(FaceIJK(8,  CoordIJK(0, 1, 0)), 0, (0, 0)),      # base cell 62
    BaseCellData(FaceIJK(6,  CoordIJK(2, 0, 0)), 1, (11, 15)),    # base cell 63
    BaseCellData(FaceIJK(8,  CoordIJK(0, 0, 0)), 0, (0, 0)),      # base cell 64
    BaseCellData(FaceIJK(9,  CoordIJK(0, 0, 1)), 0, (0, 0)),      # base cell 65
    BaseCellData(FaceIJK(14, CoordIJK(1, 0, 0)), 0, (0, 0)),      # base cell 66
    BaseCellData(FaceIJK(5,  CoordIJK(1, 0, 1)), 0, (0, 0)),      # base cell 67
    BaseCellData(FaceIJK(16, CoordIJK(0, 1, 1)), 0, (0, 0)),      # base cell 68
    BaseCellData(FaceIJK(8,  CoordIJK(1, 0, 1)), 0, (0, 0)),      # base cell 69
    BaseCellData(FaceIJK(5,  CoordIJK(1, 0, 0)), 0, (0, 0)),      # base cell 70
    BaseCellData(FaceIJK(12, CoordIJK(0, 0, 0)), 0, (0, 0)),      # base cell 71
    BaseCellData(FaceIJK(7,  CoordIJK(2, 0, 0)), 1, (12, 16)),    # base cell 72
    BaseCellData(FaceIJK(12, CoordIJK(0, 1, 0)), 0, (0, 0)),      # base cell 73
    BaseCellData(FaceIJK(10, CoordIJK(0, 1, 0)), 0, (0, 0)),      # base cell 74
    BaseCellData(FaceIJK(9,  CoordIJK(0, 0, 0)), 0, (0, 0)),      # base cell 75
    BaseCellData(FaceIJK(13, CoordIJK(1, 0, 0)), 0, (0, 0)),      # base cell 76
    BaseCellData(FaceIJK(16, CoordIJK(0, 0, 1)), 0, (0, 0)),      # base cell 77
    BaseCellData(FaceIJK(15, CoordIJK(0, 1, 1)), 0, (0, 0)),      # base cell 78
    BaseCellData(FaceIJK(15, CoordIJK(0, 1, 0)), 0, (0, 0)),      # base cell 79
    BaseCellData(FaceIJK(16, CoordIJK(0, 1, 0)), 0, (0, 0)),      # base cell 80
    BaseCellData(FaceIJK(14, CoordIJK(1, 1, 0)), 0, (0, 0)),      # base cell 81
    BaseCellData(FaceIJK(13, CoordIJK(1, 1, 0)), 0, (0, 0)),      # base cell 82
    BaseCellData(FaceIJK(5,  CoordIJK(2, 0, 0)), 1, (10, 19)),    # base cell 83
    BaseCellData(FaceIJK(8,  CoordIJK(1, 0, 0)), 0, (0, 0)),      # base cell 84
    BaseCellData(FaceIJK(14, CoordIJK(0, 0, 0)), 0, (0, 0)),      # base cell 85
    BaseCellData(FaceIJK(9,  CoordIJK(1, 0, 1)), 0, (0, 0)),      # base cell 86
    BaseCellData(FaceIJK(14, CoordIJK(0, 0, 1)), 0, (0, 0)),      # base cell 87
    BaseCellData(FaceIJK(17, CoordIJK(0, 0, 1)), 0, (0, 0)),      # base cell 88
    BaseCellData(FaceIJK(12, CoordIJK(0, 0, 1)), 0, (0, 0)),      # base cell 89
    BaseCellData(FaceIJK(16, CoordIJK(0, 0, 0)), 0, (0, 0)),      # base cell 90
    BaseCellData(FaceIJK(17, CoordIJK(0, 1, 1)), 0, (0, 0)),      # base cell 91
    BaseCellData(FaceIJK(15, CoordIJK(0, 0, 1)), 0, (0, 0)),      # base cell 92
    BaseCellData(FaceIJK(16, CoordIJK(1, 0, 1)), 0, (0, 0)),      # base cell 93
    BaseCellData(FaceIJK(9,  CoordIJK(1, 0, 0)), 0, (0, 0)),      # base cell 94
    BaseCellData(FaceIJK(15, CoordIJK(0, 0, 0)), 0, (0, 0)),      # base cell 95
    BaseCellData(FaceIJK(13, CoordIJK(0, 0, 0)), 0, (0, 0)),      # base cell 96
    BaseCellData(FaceIJK(8,  CoordIJK(2, 0, 0)), 1, (13, 17)),    # base cell 97
    BaseCellData(FaceIJK(13, CoordIJK(0, 1, 0)), 0, (0, 0)),      # base cell 98
    BaseCellData(FaceIJK(17, CoordIJK(1, 0, 1)), 0, (0, 0)),      # base cell 99
    BaseCellData(FaceIJK(19, CoordIJK(0, 1, 0)), 0, (0, 0)),      # base cell 100
    BaseCellData(FaceIJK(14, CoordIJK(0, 1, 0)), 0, (0, 0)),      # base cell 101
    BaseCellData(FaceIJK(19, CoordIJK(0, 1, 1)), 0, (0, 0)),      # base cell 102
    BaseCellData(FaceIJK(17, CoordIJK(0, 1, 0)), 0, (0, 0)),      # base cell 103
    BaseCellData(FaceIJK(13, CoordIJK(0, 0, 1)), 0, (0, 0)),      # base cell 104
    BaseCellData(FaceIJK(17, CoordIJK(0, 0, 0)), 0, (0, 0)),      # base cell 105
    BaseCellData(FaceIJK(16, CoordIJK(1, 0, 0)), 0, (0, 0)),      # base cell 106
    BaseCellData(FaceIJK(9,  CoordIJK(2, 0, 0)), 1, (14, 18)),    # base cell 107
    BaseCellData(FaceIJK(15, CoordIJK(1, 0, 1)), 0, (0, 0)),      # base cell 108
    BaseCellData(FaceIJK(15, CoordIJK(1, 0, 0)), 0, (0, 0)),      # base cell 109
    BaseCellData(FaceIJK(18, CoordIJK(0, 1, 1)), 0, (0, 0)),      # base cell 110
    BaseCellData(FaceIJK(18, CoordIJK(0, 0, 1)), 0, (0, 0)),      # base cell 111
    BaseCellData(FaceIJK(19, CoordIJK(0, 0, 1)), 0, (0, 0)),      # base cell 112
    BaseCellData(FaceIJK(17, CoordIJK(1, 0, 0)), 0, (0, 0)),      # base cell 113
    BaseCellData(FaceIJK(19, CoordIJK(0, 0, 0)), 0, (0, 0)),      # base cell 114
    BaseCellData(FaceIJK(18, CoordIJK(0, 1, 0)), 0, (0, 0)),      # base cell 115
    BaseCellData(FaceIJK(18, CoordIJK(1, 0, 1)), 0, (0, 0)),      # base cell 116
    BaseCellData(FaceIJK(19, CoordIJK(2, 0, 0)), 1, (-1, -1)),    # base cell 117
    BaseCellData(FaceIJK(19, CoordIJK(1, 0, 0)), 0, (0, 0)),      # base cell 118
    BaseCellData(FaceIJK(18, CoordIJK(0, 0, 0)), 0, (0, 0)),      # base cell 119
    BaseCellData(FaceIJK(19, CoordIJK(1, 0, 1)), 0, (0, 0)),      # base cell 120
    BaseCellData(FaceIJK(18, CoordIJK(1, 0, 0)), 0, (0, 0)),      # base cell 121
]

# Fast pentagon lookup: true for base cells that are pentagons
const isBaseCellPentagonArr = let
    arr = falses(NUM_BASE_CELLS)
    for bc in (4, 14, 24, 38, 49, 58, 63, 72, 83, 97, 107, 117)
        arr[bc + 1] = true
    end
    arr
end

function _isBaseCellPentagon(baseCell::Int)::Bool
    if baseCell < 0 || baseCell >= NUM_BASE_CELLS
        return false
    end
    return isBaseCellPentagonArr[baseCell + 1]
end

function _isBaseCellPolarPentagon(baseCell::Int)::Bool
    return baseCell == 4 || baseCell == 117
end

function _faceIjkToBaseCell(h::FaceIJK)::Int
    return Int(faceIjkBaseCells[h.face + 1, h.coord.i + 1, h.coord.j + 1, h.coord.k + 1].baseCell)
end

function _faceIjkToBaseCellCCWrot60(h::FaceIJK)::Int
    return Int(faceIjkBaseCells[h.face + 1, h.coord.i + 1, h.coord.j + 1, h.coord.k + 1].ccwRot60)
end

function _baseCellToFaceIjk(baseCell::Int)::FaceIJK
    return baseCellData[baseCell + 1].homeFijk
end

function _baseCellToCCWrot60(baseCell::Int, face::Int)::Int
    if face < 0 || face > NUM_ICOSA_FACES
        return INVALID_ROTATIONS
    end
    for i in 1:3, j in 1:3, k in 1:3
        if faceIjkBaseCells[face + 1, i, j, k].baseCell == baseCell
            return Int(faceIjkBaseCells[face + 1, i, j, k].ccwRot60)
        end
    end
    return INVALID_ROTATIONS
end

function _baseCellIsCwOffset(baseCell::Int, testFace::Int)::Bool
    return baseCellData[baseCell + 1].cwOffsetPent[1] == testFace ||
           baseCellData[baseCell + 1].cwOffsetPent[2] == testFace
end

function _getBaseCellNeighbor(baseCell::Int, dir::Direction)::Int
    return Int(baseCellNeighbors[baseCell + 1][Int(dir) + 1])
end

function _getBaseCellNeighbor(baseCell::Int, dir::Int)::Int
    return Int(baseCellNeighbors[baseCell + 1][dir + 1])
end

function _getBaseCellDirection(originBaseCell::Int, neighboringBaseCell::Int)::Direction
    for dir_i in Int32(0):Int32(6)
        testBaseCell = _getBaseCellNeighbor(originBaseCell, Direction(dir_i))
        if testBaseCell == neighboringBaseCell
            return Direction(dir_i)
        end
    end
    return INVALID_DIGIT
end

function res0CellCount()::Int
    return NUM_BASE_CELLS
end
