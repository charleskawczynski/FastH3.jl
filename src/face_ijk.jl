# Face IJK coordinate system functions

const M_SQRT7 = 2.6457513110645905905016157536392604257102
const M_RSQRT7 = 0.37796447300922722721451653623418006081576

struct FaceIJK
    face::Int32
    coord::CoordIJK
end
FaceIJK() = FaceIJK(Int32(0), CoordIJK())
FaceIJK(face::Int, coord::CoordIJK) = FaceIJK(Int32(face), coord)

struct FaceOrientIJK
    face::Int32
    translate::CoordIJK
    ccwRot60::Int32
end

isResolutionClassIII(r::Int)::Bool = (r % 2) != 0
isResolutionClassIII(r::Int32)::Bool = (r % 2) != 0

@enum Overage::Int32 begin
    NO_OVERAGE = 0
    FACE_EDGE = 1
    NEW_FACE = 2
end

const INVALID_FACE = -1

# Face neighbor direction indices
const _IJ = 1
const _KI = 2
const _JK = 3

const faceCenterGeo = [
    LatLng(0.803582649718989942, 1.248397419617396099),
    LatLng(1.307747883455638156, 2.536945009877921159),
    LatLng(1.054751253523952054, -1.347517358900396623),
    LatLng(0.600191595538186799, -0.450603909469755746),
    LatLng(0.491715428198773866, 0.401988202911306943),
    LatLng(0.172745327415618701, 1.678146885280433686),
    LatLng(0.605929321571350690, 2.953923329812411617),
    LatLng(0.427370518328979641, -1.888876200336285401),
    LatLng(-0.079066118549212831, -0.733429513380867741),
    LatLng(-0.230961644455383637, 0.506495587332349035),
    LatLng(0.079066118549212831, 2.408163140208925497),
    LatLng(0.230961644455383637, -2.635097066257444203),
    LatLng(-0.172745327415618701, -1.463445768309359553),
    LatLng(-0.605929321571350690, -0.187669323777381622),
    LatLng(-0.427370518328979641, 1.252716453253507838),
    LatLng(-0.600191595538186799, 2.690988744120037492),
    LatLng(-0.491715428198773866, -2.739604450678486295),
    LatLng(-0.803582649718989942, -1.893195233972397139),
    LatLng(-1.307747883455638156, -0.604647643711872080),
    LatLng(-1.054751253523952054, 1.794075294689396615),
]

const faceCenterPoint = [
    Vec3d(0.2199307791404606, 0.6583691780274996, 0.7198475378926182),
    Vec3d(-0.2139234834501421, 0.1478171829550703, 0.9656017935214205),
    Vec3d(0.1092625278784797, -0.4811951572873210, 0.8697775121287253),
    Vec3d(0.7428567301586791, -0.3593941678278028, 0.5648005936517033),
    Vec3d(0.8112534709140969, 0.3448953237639384, 0.4721387736413930),
    Vec3d(-0.1055498149613921, 0.9794457296411413, 0.1718874610009365),
    Vec3d(-0.8075407579970092, 0.1533552485898818, 0.5695261994882688),
    Vec3d(-0.2846148069787907, -0.8644080972654206, 0.4144792552473539),
    Vec3d(0.7405621473854482, -0.6673299564565524, -0.0789837646326737),
    Vec3d(0.8512303986474293, 0.4722343788582681, -0.2289137388687808),
    Vec3d(-0.7405621473854481, 0.6673299564565524, 0.0789837646326737),
    Vec3d(-0.8512303986474292, -0.4722343788582682, 0.2289137388687808),
    Vec3d(0.1055498149613919, -0.9794457296411413, -0.1718874610009365),
    Vec3d(0.8075407579970092, -0.1533552485898819, -0.5695261994882688),
    Vec3d(0.2846148069787908, 0.8644080972654204, -0.4144792552473539),
    Vec3d(-0.7428567301586791, 0.3593941678278027, -0.5648005936517033),
    Vec3d(-0.8112534709140971, -0.3448953237639382, -0.4721387736413930),
    Vec3d(-0.2199307791404607, -0.6583691780274996, -0.7198475378926182),
    Vec3d(0.2139234834501420, -0.1478171829550704, -0.9656017935214205),
    Vec3d(-0.1092625278784796, 0.4811951572873210, -0.8697775121287253),
]

const faceAxesAzRadsCII = [
    (5.619958268523939882, 3.525563166130744542, 1.431168063737548730),
    (5.760339081714187279, 3.665943979320991689, 1.571548876927796127),
    (0.780213654393430055, 4.969003859179821079, 2.874608756786625655),
    (0.430469363979999913, 4.619259568766391033, 2.524864466373195467),
    (6.130269123335111400, 4.035874020941915804, 1.941478918548720291),
    (2.692877706530642877, 0.598482604137447119, 4.787272808923838195),
    (2.982963003477243874, 0.888567901084048369, 5.077358105870439581),
    (3.532912002790141181, 1.438516900396945656, 5.627307105183336758),
    (3.494305004259568154, 1.399909901866372864, 5.588700106652763840),
    (3.003214169499538391, 0.908819067106342928, 5.097609271892733906),
    (5.930472956509811562, 3.836077854116615875, 1.741682751723420374),
    (0.138378484090254847, 4.327168688876645809, 2.232773586483450311),
    (0.448714947059150361, 4.637505151845541521, 2.543110049452346120),
    (0.158629650112549365, 4.347419854898940135, 2.253024752505744869),
    (5.891865957979238535, 3.797470855586042958, 1.703075753192847583),
    (2.711123289609793325, 0.616728187216597771, 4.805518392002988683),
    (3.294508837434268316, 1.200113735041072948, 5.388903939827463911),
    (3.804819692245439833, 1.710424589852244509, 5.899214794638635174),
    (3.664438879055192436, 1.570043776661997111, 5.758833981448388027),
    (2.361378999196363184, 0.266983896803167583, 4.455774101589558636),
]

# Face neighbor table: faceNeighbors[face+1][dir] where dir: 1=central, 2=IJ, 3=KI, 4=JK
const faceNeighbors = [
    # face 0
    (FaceOrientIJK(0, CoordIJK(0, 0, 0), 0), FaceOrientIJK(4, CoordIJK(2, 0, 2), 1), FaceOrientIJK(1, CoordIJK(2, 2, 0), 5), FaceOrientIJK(5, CoordIJK(0, 2, 2), 3)),
    # face 1
    (FaceOrientIJK(1, CoordIJK(0, 0, 0), 0), FaceOrientIJK(0, CoordIJK(2, 0, 2), 1), FaceOrientIJK(2, CoordIJK(2, 2, 0), 5), FaceOrientIJK(6, CoordIJK(0, 2, 2), 3)),
    # face 2
    (FaceOrientIJK(2, CoordIJK(0, 0, 0), 0), FaceOrientIJK(1, CoordIJK(2, 0, 2), 1), FaceOrientIJK(3, CoordIJK(2, 2, 0), 5), FaceOrientIJK(7, CoordIJK(0, 2, 2), 3)),
    # face 3
    (FaceOrientIJK(3, CoordIJK(0, 0, 0), 0), FaceOrientIJK(2, CoordIJK(2, 0, 2), 1), FaceOrientIJK(4, CoordIJK(2, 2, 0), 5), FaceOrientIJK(8, CoordIJK(0, 2, 2), 3)),
    # face 4
    (FaceOrientIJK(4, CoordIJK(0, 0, 0), 0), FaceOrientIJK(3, CoordIJK(2, 0, 2), 1), FaceOrientIJK(0, CoordIJK(2, 2, 0), 5), FaceOrientIJK(9, CoordIJK(0, 2, 2), 3)),
    # face 5
    (FaceOrientIJK(5, CoordIJK(0, 0, 0), 0), FaceOrientIJK(10, CoordIJK(2, 2, 0), 3), FaceOrientIJK(14, CoordIJK(2, 0, 2), 3), FaceOrientIJK(0, CoordIJK(0, 2, 2), 3)),
    # face 6
    (FaceOrientIJK(6, CoordIJK(0, 0, 0), 0), FaceOrientIJK(11, CoordIJK(2, 2, 0), 3), FaceOrientIJK(10, CoordIJK(2, 0, 2), 3), FaceOrientIJK(1, CoordIJK(0, 2, 2), 3)),
    # face 7
    (FaceOrientIJK(7, CoordIJK(0, 0, 0), 0), FaceOrientIJK(12, CoordIJK(2, 2, 0), 3), FaceOrientIJK(11, CoordIJK(2, 0, 2), 3), FaceOrientIJK(2, CoordIJK(0, 2, 2), 3)),
    # face 8
    (FaceOrientIJK(8, CoordIJK(0, 0, 0), 0), FaceOrientIJK(13, CoordIJK(2, 2, 0), 3), FaceOrientIJK(12, CoordIJK(2, 0, 2), 3), FaceOrientIJK(3, CoordIJK(0, 2, 2), 3)),
    # face 9
    (FaceOrientIJK(9, CoordIJK(0, 0, 0), 0), FaceOrientIJK(14, CoordIJK(2, 2, 0), 3), FaceOrientIJK(13, CoordIJK(2, 0, 2), 3), FaceOrientIJK(4, CoordIJK(0, 2, 2), 3)),
    # face 10
    (FaceOrientIJK(10, CoordIJK(0, 0, 0), 0), FaceOrientIJK(5, CoordIJK(2, 2, 0), 3), FaceOrientIJK(6, CoordIJK(2, 0, 2), 3), FaceOrientIJK(15, CoordIJK(0, 2, 2), 3)),
    # face 11
    (FaceOrientIJK(11, CoordIJK(0, 0, 0), 0), FaceOrientIJK(6, CoordIJK(2, 2, 0), 3), FaceOrientIJK(7, CoordIJK(2, 0, 2), 3), FaceOrientIJK(16, CoordIJK(0, 2, 2), 3)),
    # face 12
    (FaceOrientIJK(12, CoordIJK(0, 0, 0), 0), FaceOrientIJK(7, CoordIJK(2, 2, 0), 3), FaceOrientIJK(8, CoordIJK(2, 0, 2), 3), FaceOrientIJK(17, CoordIJK(0, 2, 2), 3)),
    # face 13
    (FaceOrientIJK(13, CoordIJK(0, 0, 0), 0), FaceOrientIJK(8, CoordIJK(2, 2, 0), 3), FaceOrientIJK(9, CoordIJK(2, 0, 2), 3), FaceOrientIJK(18, CoordIJK(0, 2, 2), 3)),
    # face 14
    (FaceOrientIJK(14, CoordIJK(0, 0, 0), 0), FaceOrientIJK(9, CoordIJK(2, 2, 0), 3), FaceOrientIJK(5, CoordIJK(2, 0, 2), 3), FaceOrientIJK(19, CoordIJK(0, 2, 2), 3)),
    # face 15
    (FaceOrientIJK(15, CoordIJK(0, 0, 0), 0), FaceOrientIJK(16, CoordIJK(2, 0, 2), 1), FaceOrientIJK(19, CoordIJK(2, 2, 0), 5), FaceOrientIJK(10, CoordIJK(0, 2, 2), 3)),
    # face 16
    (FaceOrientIJK(16, CoordIJK(0, 0, 0), 0), FaceOrientIJK(17, CoordIJK(2, 0, 2), 1), FaceOrientIJK(15, CoordIJK(2, 2, 0), 5), FaceOrientIJK(11, CoordIJK(0, 2, 2), 3)),
    # face 17
    (FaceOrientIJK(17, CoordIJK(0, 0, 0), 0), FaceOrientIJK(18, CoordIJK(2, 0, 2), 1), FaceOrientIJK(16, CoordIJK(2, 2, 0), 5), FaceOrientIJK(12, CoordIJK(0, 2, 2), 3)),
    # face 18
    (FaceOrientIJK(18, CoordIJK(0, 0, 0), 0), FaceOrientIJK(19, CoordIJK(2, 0, 2), 1), FaceOrientIJK(17, CoordIJK(2, 2, 0), 5), FaceOrientIJK(13, CoordIJK(0, 2, 2), 3)),
    # face 19
    (FaceOrientIJK(19, CoordIJK(0, 0, 0), 0), FaceOrientIJK(15, CoordIJK(2, 0, 2), 1), FaceOrientIJK(18, CoordIJK(2, 2, 0), 5), FaceOrientIJK(14, CoordIJK(0, 2, 2), 3)),
]

# adjacentFaceDir[face+1][otherface+1]: direction from origin face to dest face, or -1 if not adjacent
const adjacentFaceDir = [
    [0, _KI, -1, -1, _IJ, _JK, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [_IJ, 0, _KI, -1, -1, -1, _JK, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [-1, _IJ, 0, _KI, -1, -1, -1, _JK, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [-1, -1, _IJ, 0, _KI, -1, -1, -1, _JK, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [_KI, -1, -1, _IJ, 0, -1, -1, -1, -1, _JK, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [_JK, -1, -1, -1, -1, 0, -1, -1, -1, -1, _IJ, -1, -1, -1, _KI, -1, -1, -1, -1, -1],
    [-1, _JK, -1, -1, -1, -1, 0, -1, -1, -1, _KI, _IJ, -1, -1, -1, -1, -1, -1, -1, -1],
    [-1, -1, _JK, -1, -1, -1, -1, 0, -1, -1, -1, _KI, _IJ, -1, -1, -1, -1, -1, -1, -1],
    [-1, -1, -1, _JK, -1, -1, -1, -1, 0, -1, -1, -1, _KI, _IJ, -1, -1, -1, -1, -1, -1],
    [-1, -1, -1, -1, _JK, -1, -1, -1, -1, 0, -1, -1, -1, _KI, _IJ, -1, -1, -1, -1, -1],
    [-1, -1, -1, -1, -1, _IJ, _KI, -1, -1, -1, 0, -1, -1, -1, -1, _JK, -1, -1, -1, -1],
    [-1, -1, -1, -1, -1, -1, _IJ, _KI, -1, -1, -1, 0, -1, -1, -1, -1, _JK, -1, -1, -1],
    [-1, -1, -1, -1, -1, -1, -1, _IJ, _KI, -1, -1, -1, 0, -1, -1, -1, -1, _JK, -1, -1],
    [-1, -1, -1, -1, -1, -1, -1, -1, _IJ, _KI, -1, -1, -1, 0, -1, -1, -1, -1, _JK, -1],
    [-1, -1, -1, -1, -1, _KI, -1, -1, -1, _IJ, -1, -1, -1, -1, 0, -1, -1, -1, -1, _JK],
    [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, _JK, -1, -1, -1, -1, 0, _IJ, -1, -1, _KI],
    [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, _JK, -1, -1, -1, _KI, 0, _IJ, -1, -1],
    [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, _JK, -1, -1, -1, _KI, 0, _IJ, -1],
    [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, _JK, -1, -1, -1, _KI, 0, _IJ],
    [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, _JK, _IJ, -1, -1, _KI, 0],
]

# Overage distance table (1-indexed: index = res+1)
const maxDimByCIIres = [
    2, -1, 14, -1, 98, -1, 686, -1, 4802, -1,
    33614, -1, 235298, -1, 1647086, -1, 11529602
]

# Unit scale distance table (1-indexed: index = res+1)
const unitScaleByCIIres = [
    1, -1, 7, -1, 49, -1, 343, -1, 2401, -1,
    16807, -1, 117649, -1, 823543, -1, 5764801
]

function _geoToClosestFace(g::LatLng)::Tuple{Int,Float64}
    v3d = _geoToVec3d(g)
    face = 0
    sqd = 5.0
    for f in 0:(NUM_ICOSA_FACES-1)
        sqdT = _pointSquareDist(faceCenterPoint[f+1], v3d)
        if sqdT < sqd
            face = f
            sqd = sqdT
        end
    end
    return (face, sqd)
end

function _geoToHex2d(g::LatLng, res::Int)::Tuple{Int,Vec2d}
    face, sqd = _geoToClosestFace(g)

    r = acos(clamp(1.0 - sqd * 0.5, -1.0, 1.0))

    if r < EPSILON
        return (face, Vec2d(0.0, 0.0))
    end

    theta = _posAngleRads(faceAxesAzRadsCII[face+1][1] -
                          _posAngleRads(_geoAzimuthRads(faceCenterGeo[face+1], g)))

    if isResolutionClassIII(res)
        theta = _posAngleRads(theta - M_AP7_ROT_RADS)
    end

    r = tan(r)
    r *= INV_RES0_U_GNOMONIC
    for _ in 1:res
        r *= M_SQRT7
    end

    return (face, Vec2d(r * cos(theta), r * sin(theta)))
end

function _geoToFaceIjk(g::LatLng, res::Int)::FaceIJK
    face, v = _geoToHex2d(g, res)
    coord = _hex2dToCoordIJK(v)
    return FaceIJK(face, coord)
end

function _hex2dToGeo(v::Vec2d, face::Int, res::Int, substrate::Int)::LatLng
    r = _v2dMag(v)

    if r < EPSILON
        return faceCenterGeo[face+1]
    end

    theta = atan(v.y, v.x)

    for _ in 1:res
        r *= M_RSQRT7
    end

    if substrate != 0
        r *= M_ONETHIRD
        if isResolutionClassIII(res)
            r *= M_RSQRT7
        end
    end

    r *= RES0_U_GNOMONIC
    r = atan(r)

    if substrate == 0 && isResolutionClassIII(res)
        theta = _posAngleRads(theta + M_AP7_ROT_RADS)
    end

    theta = _posAngleRads(faceAxesAzRadsCII[face+1][1] - theta)

    return _geoAzDistanceRads(faceCenterGeo[face+1], theta, r)
end

function _faceIjkToGeo(h::FaceIJK, res::Int)::LatLng
    v = _ijkToHex2d(h.coord)
    return _hex2dToGeo(v, Int(h.face), res, 0)
end

function _faceIjkPentToVerts(fijk::FaceIJK, res::Int)::Tuple{NTuple{NUM_PENT_VERTS,FaceIJK},Int}
    vertsCII = (
        CoordIJK(2, 1, 0), CoordIJK(1, 2, 0), CoordIJK(0, 2, 1),
        CoordIJK(0, 1, 2), CoordIJK(1, 0, 2),
    )
    vertsCIII = (
        CoordIJK(5, 4, 0), CoordIJK(1, 5, 0), CoordIJK(0, 5, 4),
        CoordIJK(0, 1, 5), CoordIJK(4, 0, 5),
    )

    verts = isResolutionClassIII(res) ? vertsCIII : vertsCII

    coord = _downAp3r(_downAp3(fijk.coord))

    if isResolutionClassIII(res)
        coord = _downAp7r(coord)
        res += 1
    end

    face = fijk.face
    fijkVerts = let coord = coord, verts = verts
        ntuple(Val(NUM_PENT_VERTS)) do v
            FaceIJK(face, _ijkNormalize(_ijkAdd(coord, verts[v])))
        end
    end
    return (fijkVerts, res)
end

function _faceIjkToVerts(fijk::FaceIJK, res::Int)::Tuple{NTuple{NUM_HEX_VERTS,FaceIJK},Int}
    vertsCII = (
        CoordIJK(2, 1, 0), CoordIJK(1, 2, 0), CoordIJK(0, 2, 1),
        CoordIJK(0, 1, 2), CoordIJK(1, 0, 2), CoordIJK(2, 0, 1),
    )
    vertsCIII = (
        CoordIJK(5, 4, 0), CoordIJK(1, 5, 0), CoordIJK(0, 5, 4),
        CoordIJK(0, 1, 5), CoordIJK(4, 0, 5), CoordIJK(5, 0, 1),
    )

    verts = isResolutionClassIII(res) ? vertsCIII : vertsCII

    coord = _downAp3r(_downAp3(fijk.coord))

    if isResolutionClassIII(res)
        coord = _downAp7r(coord)
        res += 1
    end

    face = fijk.face
    fijkVerts = let coord = coord, verts = verts
        ntuple(Val(NUM_HEX_VERTS)) do v
            FaceIJK(face, _ijkNormalize(_ijkAdd(coord, verts[v])))
        end
    end
    return (fijkVerts, res)
end

function _adjustOverageClassII(fijk::FaceIJK, res::Int, pentLeading4::Int, substrate::Int)::Tuple{Overage,FaceIJK}
    overage = NO_OVERAGE
    ijk = fijk.coord
    face = fijk.face

    maxDim = maxDimByCIIres[res+1]
    if substrate != 0
        maxDim *= 3
    end

    if substrate != 0 && ijk.i + ijk.j + ijk.k == maxDim
        overage = FACE_EDGE
    elseif ijk.i + ijk.j + ijk.k > maxDim
        overage = NEW_FACE

        if ijk.k > 0
            if ijk.j > 0  # jk quadrant
                fijkOrient = faceNeighbors[face+1][_JK+1]
            else  # ki quadrant
                fijkOrient = faceNeighbors[face+1][_KI+1]

                if pentLeading4 != 0
                    origin = CoordIJK(maxDim, 0, 0)
                    tmp = _ijkSub(ijk, origin)
                    tmp = _ijkRotate60cw(tmp)
                    ijk = _ijkAdd(tmp, origin)
                end
            end
        else  # ij quadrant
            fijkOrient = faceNeighbors[face+1][_IJ+1]
        end

        face = fijkOrient.face

        for _ in 1:fijkOrient.ccwRot60
            ijk = _ijkRotate60ccw(ijk)
        end

        transVec = _ijkScale(fijkOrient.translate, unitScaleByCIIres[res+1] * (substrate != 0 ? 3 : 1))
        ijk = _ijkNormalize(_ijkAdd(ijk, transVec))

        if substrate != 0 && ijk.i + ijk.j + ijk.k == maxDim
            overage = FACE_EDGE
        end
    end

    return (overage, FaceIJK(face, ijk))
end

function _adjustPentVertOverage(fijk::FaceIJK, res::Int)::Tuple{Overage,FaceIJK}
    overage = NO_OVERAGE
    while true
        overage, fijk = _adjustOverageClassII(fijk, res, 0, 1)
        overage != NEW_FACE && break
    end
    return (overage, fijk)
end

function _faceIjkPentToCellBoundary(h::FaceIJK, res::Int, start::Int, length::Int)::CellBoundary
    adjRes = res
    fijkVerts, adjRes = _faceIjkPentToVerts(h, adjRes)

    additionalIteration = length == NUM_PENT_VERTS ? 1 : 0

    numVerts = Int32(0)
    verts_arr = ntuple(_ -> LatLng(), MAX_CELL_BNDRY_VERTS)
    lastFijk = FaceIJK()

    for vert in start:(start+length-1+additionalIteration)
        v = mod(vert, NUM_PENT_VERTS) + 1  # 1-indexed

        _, fijk = _adjustPentVertOverage(fijkVerts[v], adjRes)

        if isResolutionClassIII(res) && vert > start
            tmpFijk = fijk
            orig2d0 = _ijkToHex2d(lastFijk.coord)

            currentToLastDir = adjacentFaceDir[tmpFijk.face+1][lastFijk.face+1]
            fijkOrient = faceNeighbors[tmpFijk.face+1][currentToLastDir+1]

            tmpCoord = tmpFijk.coord
            tmpFace = fijkOrient.face

            for _ in 1:fijkOrient.ccwRot60
                tmpCoord = _ijkRotate60ccw(tmpCoord)
            end

            transVec = _ijkScale(fijkOrient.translate, unitScaleByCIIres[adjRes+1] * 3)
            tmpCoord = _ijkNormalize(_ijkAdd(tmpCoord, transVec))

            orig2d1 = _ijkToHex2d(tmpCoord)

            maxDim_val = maxDimByCIIres[adjRes+1]
            v0 = Vec2d(3.0 * maxDim_val, 0.0)
            v1 = Vec2d(-1.5 * maxDim_val, 3.0 * M_SQRT3_2 * maxDim_val)
            v2 = Vec2d(-1.5 * maxDim_val, -3.0 * M_SQRT3_2 * maxDim_val)

            faceDir = adjacentFaceDir[tmpFace+1][fijk.face+1]
            if faceDir == _IJ
                edge0, edge1 = v0, v1
            elseif faceDir == _JK
                edge0, edge1 = v1, v2
            else
                edge0, edge1 = v2, v0
            end

            inter = _v2dIntersect(orig2d0, orig2d1, edge0, edge1)
            numVerts += Int32(1)
            verts_arr = Base.setindex(verts_arr, _hex2dToGeo(inter, Int(tmpFace), adjRes, 1), numVerts)
        end

        if vert < start + NUM_PENT_VERTS
            vec = _ijkToHex2d(fijk.coord)
            numVerts += Int32(1)
            verts_arr = Base.setindex(verts_arr, _hex2dToGeo(vec, Int(fijk.face), adjRes, 1), numVerts)
        end

        lastFijk = fijk
    end

    return CellBoundary(numVerts, verts_arr)
end

function _faceIjkToCellBoundary(h::FaceIJK, res::Int, start::Int, length::Int)::CellBoundary
    adjRes = res
    fijkVerts, adjRes = _faceIjkToVerts(h, adjRes)
    centerFace = h.face

    additionalIteration = length == NUM_HEX_VERTS ? 1 : 0

    numVerts = Int32(0)
    verts_arr = ntuple(_ -> LatLng(), MAX_CELL_BNDRY_VERTS)
    lastFace = -1
    lastOverage = NO_OVERAGE

    for vert in start:(start+length-1+additionalIteration)
        v = mod(vert, NUM_HEX_VERTS) + 1  # 1-indexed

        overage, fijk = _adjustOverageClassII(fijkVerts[v], adjRes, 0, 1)

        if isResolutionClassIII(res) && vert > start &&
           fijk.face != lastFace && lastOverage != FACE_EDGE
            lastV = mod(vert - 1, NUM_HEX_VERTS) + 1
            orig2d0 = _ijkToHex2d(fijkVerts[lastV].coord)
            orig2d1 = _ijkToHex2d(fijkVerts[v].coord)

            maxDim_val = maxDimByCIIres[adjRes+1]
            v0 = Vec2d(3.0 * maxDim_val, 0.0)
            v1 = Vec2d(-1.5 * maxDim_val, 3.0 * M_SQRT3_2 * maxDim_val)
            v2 = Vec2d(-1.5 * maxDim_val, -3.0 * M_SQRT3_2 * maxDim_val)

            face2 = (lastFace == centerFace) ? fijk.face : lastFace
            faceDir = adjacentFaceDir[centerFace+1][face2+1]
            if faceDir == _IJ
                edge0, edge1 = v0, v1
            elseif faceDir == _JK
                edge0, edge1 = v1, v2
            else
                edge0, edge1 = v2, v0
            end

            inter = _v2dIntersect(orig2d0, orig2d1, edge0, edge1)
            isIntersectionAtVertex = _v2dAlmostEquals(orig2d0, inter) ||
                                     _v2dAlmostEquals(orig2d1, inter)
            if !isIntersectionAtVertex
                numVerts += Int32(1)
                verts_arr = Base.setindex(verts_arr, _hex2dToGeo(inter, Int(centerFace), adjRes, 1), numVerts)
            end
        end

        if vert < start + NUM_HEX_VERTS
            vec = _ijkToHex2d(fijk.coord)
            numVerts += Int32(1)
            verts_arr = Base.setindex(verts_arr, _hex2dToGeo(vec, Int(fijk.face), adjRes, 1), numVerts)
        end

        lastFace = fijk.face
        lastOverage = overage
    end

    return CellBoundary(numVerts, verts_arr)
end
