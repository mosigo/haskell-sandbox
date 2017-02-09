module AminoAcid
        ( AminoAcid (..)
        , Radical (..)
        , Hydrated (..)
        , HydratedAminoAcid
        , AtomType (..)
        , atomCoordinates
        , connections
        ) where

import Linear.V3 (V3)

data AminoAcid a = AminoAcid { nitro       :: a
                             , carbonAlpha :: a
                             , carbon      :: a
                             , oxi2        :: a
                             , oxi         :: Maybe a
                             , radical     :: Radical a
                             }
  deriving Show

data Radical a = Alanine a
               | Glysine
               | Valine a a a
               | Isoleucine a a a a
               | Leucine a a a a
               | Methionine a a a a
               | Phenylalanine a a a a a a a
               | Tyrosine a a a a a a a a
               | Tryptophan a a a a a a a a a a
               | Serine a a
               | Threonine a a a
               | Asparagine a a a a
               | Glutamine a a a a a
               | Cysteine a a
               | Proline a a a
               | Arginine a a a a a a a
               | Histidine a a a a a a
               | Lysine a a a a a
               | AsparticAcid a a a a
               | GlutamicAcid a a a a a
  deriving Show

data Hydrated a = Hydrated a [a]
  deriving Show

data AtomType = N | CA | C | O | OXT
              | CB
              | CD | CD1 | CD2
              | CE | CE1 | CE2 | CE3
              | CG | CG1 | CG2
              | CZ | CZ2 | CZ3
              | CH | CH2
              | SD
              | SG
              | OD | OD1 | OD2
              | OH | OG | OG1
              | OE | OE1 | OE2
              | ND1 | ND2
              | NE | NE1 | NE2
              | NH1 | NH2
              | NZ
              | H | H1 | H2 | H3 | HXT
              | HA | HA2 | HA3
              | HB | HB1 | HB2 | HB3 | HD11 | HD12 | HD13
              | HG | HG1 | HG2 | HG3 | HG11 | HG12 | HG13 | HG21 | HG22 | HG23
              | HD | HD1 | HD2 | HD3 | HD21 | HD22 | HD23
              | HE | HE1 | HE2 | HE3 | HE21 | HE22
              | HZ | HZ1 | HZ2 | HZ3
              | HH | HH2 | HH11 | HH12 | HH21 | HH22
  deriving (Show, Read, Eq)

type HydratedAminoAcid = AminoAcid (Hydrated (V3 Float))

atomCoordinates :: HydratedAminoAcid -> AtomType -> V3 Float
atomCoordinates aminoAcid atomType =
  case atomType of
    N   -> nitroCoords
    CA  -> carbonAlphaCoords
    C   -> carbonCoords
    O   -> oxi2Coords
    OXT -> case oxi aminoAcid of
      Just (Hydrated oxiCoords _) -> oxiCoords
      Nothing -> error "No atom type OXT"
    HXT -> case oxi aminoAcid of
      Just (Hydrated _ oxiHydrogensCoords) -> head oxiHydrogensCoords
      Nothing -> error "No atom type HXT"
    HA  -> head carbonAlphaHydrogensCoords
    HA2 -> carbonAlphaHydrogensCoords !! 1
    HA3 -> carbonAlphaHydrogensCoords !! 2
    H   -> head nitroHydrogensCoords
    H1  -> head nitroHydrogensCoords
    H2  -> nitroHydrogensCoords !! 1
    H3  -> nitroHydrogensCoords !! 2
    _   -> radicalAtomCoordinates (radical aminoAcid) atomType
  where Hydrated nitroCoords       nitroHydrogensCoords       = nitro       aminoAcid
        Hydrated carbonAlphaCoords carbonAlphaHydrogensCoords = carbonAlpha aminoAcid
        Hydrated carbonCoords      _                          = carbon      aminoAcid
        Hydrated oxi2Coords        _                          = oxi2        aminoAcid

radicalAtomCoordinates :: Radical (Hydrated (V3 Float)) -> AtomType -> V3 Float
radicalAtomCoordinates (Alanine (Hydrated atomCoords hydrogensCoords)) atomType =
  case atomType of
    CB  -> atomCoords
    HB1 -> head hydrogensCoords
    HB2 -> hydrogensCoords !! 1
    HB3 -> hydrogensCoords !! 2
    _   -> error $ "No atom type " ++ show atomType ++ "in Alanine"
radicalAtomCoordinates Glysine atomType = error $ "No atom type " ++ show atomType ++ "in Glysine, because it hasn't radical"
radicalAtomCoordinates (Valine cb cg1 cg2) atomType =
  let (Hydrated cbCoords cbHydrogensCoords)   = cb
      (Hydrated cg1Coords cg1HydrogensCoords) = cg1
      (Hydrated cg2Coords cg2HydrogensCoords) = cg2
  in case atomType of
    CB   -> cbCoords
    CG1  -> cg1Coords
    CG2  -> cg2Coords
    HB   -> head cbHydrogensCoords
    HG11 -> head cg1HydrogensCoords
    HG12 -> cg1HydrogensCoords !! 1
    HG13 -> cg1HydrogensCoords !! 2
    HG21 -> head cg2HydrogensCoords
    HG22 -> cg2HydrogensCoords !! 1
    HG23 -> cg2HydrogensCoords !! 2
    _    -> error $ "No atom type " ++ show atomType ++ "in Valine"
radicalAtomCoordinates (Isoleucine cb cg1 cg2 cd1) atomType =
  let (Hydrated cbCoords cbHydrogensCoords)   = cb
      (Hydrated cg1Coords cg1HydrogensCoords) = cg1
      (Hydrated cg2Coords cg2HydrogensCoords) = cg2
      (Hydrated cd1Coords cd1HydrogensCoords) = cd1
  in case atomType of
    CB   -> cbCoords
    CG1  -> cg1Coords
    CG2  -> cg2Coords
    CD1  -> cd1Coords
    HB   -> head cbHydrogensCoords
    HG12 -> head cg1HydrogensCoords
    HG13 -> cg1HydrogensCoords !! 1
    HG21 -> head cg2HydrogensCoords
    HG22 -> cg2HydrogensCoords !! 1
    HG23 -> cg2HydrogensCoords !! 2
    HD11 -> head cd1HydrogensCoords
    HD12 -> cd1HydrogensCoords !! 1
    HD13 -> cd1HydrogensCoords !! 2
    _    -> error $ "No atom type " ++ show atomType ++ "in Isoleucine"
radicalAtomCoordinates (Leucine cb cg cd1 cd2) atomType =
  let (Hydrated cbCoords cbHydrogensCoords)   = cb
      (Hydrated cgCoords cgHydrogensCoords)   = cg
      (Hydrated cd1Coords cd1HydrogensCoords) = cd1
      (Hydrated cd2Coords cd2HydrogensCoords) = cd2
  in case atomType of
    CB   -> cbCoords
    CG   -> cgCoords
    CD1  -> cd1Coords
    CD2  -> cd2Coords
    HB2  -> head cbHydrogensCoords
    HB3  -> cbHydrogensCoords !! 1
    HG   -> head cgHydrogensCoords
    HD11 -> head cd1HydrogensCoords
    HD12 -> cd1HydrogensCoords !! 1
    HD13 -> cd1HydrogensCoords !! 2
    HD21 -> head cd2HydrogensCoords
    HD22 -> cd2HydrogensCoords !! 1
    HD23 -> cd2HydrogensCoords !! 2
    _    -> error $ "No atom type " ++ show atomType ++ "in Leucine"
radicalAtomCoordinates (Methionine cb cg sd ce) atomType =
  let (Hydrated cbCoords cbHydrogensCoords) = cb
      (Hydrated cgCoords cgHydrogensCoords) = cg
      (Hydrated sdCoords _                ) = sd
      (Hydrated ceCoords ceHydrogensCoords) = ce
  in case atomType of
    CB   -> cbCoords
    CG   -> cgCoords
    SD   -> sdCoords
    CE   -> ceCoords
    HE1  -> head ceHydrogensCoords
    HE2  -> ceHydrogensCoords !! 1
    HE3  -> ceHydrogensCoords !! 2
    HG2  -> head cgHydrogensCoords
    HG3  -> cgHydrogensCoords !! 1
    HB2  -> head cbHydrogensCoords
    HB3  -> cbHydrogensCoords !! 1
    _    -> error $ "No atom type " ++ show atomType ++ "in Methionine"
radicalAtomCoordinates (Phenylalanine cb cg cd1 ce1 cz ce2 cd2) atomType =
  let (Hydrated cbCoords cbHydrogensCoords) = cb
      (Hydrated cgCoords _                ) = cg
      (Hydrated cd1Coords cd1HydrogensCoords) = cd1
      (Hydrated ce1Coords ce1HydrogensCoords) = ce1
      (Hydrated czCoords czHydrogensCoords) = cz
      (Hydrated ce2Coords ce2HydrogensCoords) = ce2
      (Hydrated cd2Coords cd2HydrogensCoords) = cd2
  in case atomType of
    CB   -> cbCoords
    CG   -> cgCoords
    CD1  -> cd1Coords
    CE1  -> ce1Coords
    CZ   -> czCoords
    CE2  -> ce2Coords
    CD2  -> cd2Coords
    HB2  -> head cbHydrogensCoords
    HB3  -> cbHydrogensCoords !! 1
    HD1  -> head cd1HydrogensCoords
    HE1  -> head ce1HydrogensCoords
    HZ   -> head czHydrogensCoords
    HE2  -> head ce2HydrogensCoords
    HD2  -> head cd2HydrogensCoords
    _    -> error $ "No atom type " ++ show atomType ++ "in Phenylalanine"
radicalAtomCoordinates (Tyrosine cb cg cd1 ce1 cz ce2 cd2 oh) atomType =
  let (Hydrated cbCoords cbHydrogensCoords) = cb
      (Hydrated cgCoords _                ) = cg
      (Hydrated cd1Coords cd1HydrogensCoords) = cd1
      (Hydrated ce1Coords ce1HydrogensCoords) = ce1
      (Hydrated czCoords czHydrogensCoords) = cz
      (Hydrated ce2Coords ce2HydrogensCoords) = ce2
      (Hydrated cd2Coords cd2HydrogensCoords) = cd2
      (Hydrated ohCoords ohHydrogensCoords) = oh
  in case atomType of
    CB   -> cbCoords
    CG   -> cgCoords
    CD1  -> cd1Coords
    CE1  -> ce1Coords
    CZ   -> czCoords
    CE2  -> ce2Coords
    CD2  -> cd2Coords
    OH   -> ohCoords
    HB2  -> head cbHydrogensCoords
    HB3  -> cbHydrogensCoords !! 1
    HD1  -> head cd1HydrogensCoords
    HE1  -> head ce1HydrogensCoords
    HZ   -> head czHydrogensCoords
    HE2  -> head ce2HydrogensCoords
    HD2  -> head cd2HydrogensCoords
    HH   -> head ohHydrogensCoords
    _    -> error $ "No atom type " ++ show atomType ++ "in Tyrosine"
radicalAtomCoordinates (Tryptophan cb cg cd1 ne1 ce2 cd2 ce3 cz3 ch2 cz2) atomType =
  let (Hydrated cbCoords cbHydrogensCoords) = cb
      (Hydrated cgCoords _                ) = cg
      (Hydrated cd1Coords cd1HydrogensCoords) = cd1
      (Hydrated ne1Coords ne1HydrogensCoords) = ne1
      (Hydrated ce2Coords _                 ) = ce2
      (Hydrated cd2Coords _                 ) = cd2
      (Hydrated ce3Coords ce3HydrogensCoords) = ce3
      (Hydrated cz3Coords cz3HydrogensCoords) = cz3
      (Hydrated ch2Coords ch2HydrogensCoords) = ch2
      (Hydrated cz2Coords cz2HydrogensCoords) = cz2
  in case atomType of
    CB   -> cbCoords
    CG   -> cgCoords
    CD1  -> cd1Coords
    NE1  -> ne1Coords
    CE2  -> ce2Coords
    CD2  -> cd2Coords
    CE3  -> ce3Coords
    CZ3  -> cz3Coords
    CH2  -> ch2Coords
    CZ2  -> cz2Coords
    HB2  -> head cbHydrogensCoords
    HB3  -> cbHydrogensCoords !! 1
    HD1  -> head cd1HydrogensCoords
    HE1  -> head ne1HydrogensCoords
    HZ2  -> head cz2HydrogensCoords
    HE3  -> head ce3HydrogensCoords
    HZ3  -> head cz3HydrogensCoords
    HH2  -> head ch2HydrogensCoords
    _    -> error $ "No atom type " ++ show atomType ++ "in Tryptophan"
radicalAtomCoordinates (Serine cb og) atomType =
  let (Hydrated cbCoords cbHydrogensCoords) = cb
      (Hydrated ogCoords ogHydrogensCoords) = og
  in case atomType of
    CB  -> cbCoords
    HB2 -> head cbHydrogensCoords
    HB3 -> cbHydrogensCoords !! 1
    OG  -> ogCoords
    HG  -> head ogHydrogensCoords
    _   -> error $ "No atom type " ++ show atomType ++ "in Serine"
radicalAtomCoordinates (Threonine cb og1 cg2) atomType =
  let (Hydrated cbCoords cbHydrogensCoords)   = cb
      (Hydrated og1Coords og1HydrogensCoords) = og1
      (Hydrated cg2Coords cg2HydrogensCoords) = cg2
  in case atomType of
    CB   -> cbCoords
    OG1  -> og1Coords
    CG2  -> cg2Coords
    HB   -> head cbHydrogensCoords
    HG1  -> head og1HydrogensCoords
    HG21 -> head cg2HydrogensCoords
    HG22 -> cg2HydrogensCoords !! 1
    HG23 -> cg2HydrogensCoords !! 2
    _    -> error $ "No atom type " ++ show atomType ++ "in Threonine"
radicalAtomCoordinates (Asparagine cb cg od1 nd2) atomType =
  let (Hydrated cbCoords cbHydrogensCoords)   = cb
      (Hydrated cgCoords _                ) = cg
      (Hydrated od1Coords _               ) = od1
      (Hydrated nd2Coords nd2HydrogensCoords) = nd2
  in case atomType of
    CB   -> cbCoords
    CG   -> cgCoords
    OD1  -> od1Coords
    ND2  -> nd2Coords
    HB2  -> head cbHydrogensCoords
    HB3  -> cbHydrogensCoords !! 1
    HD21 -> head nd2HydrogensCoords
    HD22 -> nd2HydrogensCoords !! 1
    _    -> error $ "No atom type " ++ show atomType ++ "in Asparagine"
radicalAtomCoordinates (Glutamine cb cg cd oe1 ne2) atomType =
  let (Hydrated cbCoords cbHydrogensCoords) = cb
      (Hydrated cgCoords cgHydrogensCoords) = cg
      (Hydrated cdCoords _                ) = cd
      (Hydrated oe1Coords _               ) = oe1
      (Hydrated ne2Coords ne2HydrogensCoords) = ne2
  in case atomType of
    CB   -> cbCoords
    CG   -> cgCoords
    CD   -> cdCoords
    OE1  -> oe1Coords
    NE2  -> ne2Coords
    HB2  -> head cbHydrogensCoords
    HB3  -> cbHydrogensCoords !! 1
    HG2  -> head cgHydrogensCoords
    HG3  -> cgHydrogensCoords !! 1
    HE21 -> head ne2HydrogensCoords
    HE22 -> ne2HydrogensCoords !! 1
    _    -> error $ "No atom type " ++ show atomType ++ "in Glutamine"
radicalAtomCoordinates (Cysteine cb sg) atomType =
  let (Hydrated cbCoords cbHydrogensCoords) = cb
      (Hydrated sgCoords sgHydrogensCoords) = sg
  in case atomType of
    CB  -> cbCoords
    HB2 -> head cbHydrogensCoords
    HB3 -> cbHydrogensCoords !! 1
    SG  -> sgCoords
    HG  -> head sgHydrogensCoords
    _   -> error $ "No atom type " ++ show atomType ++ "in Cysteine"
radicalAtomCoordinates (Proline cb cg cd) atomType =
  let (Hydrated cbCoords cbHydrogensCoords) = cb
      (Hydrated cgCoords cgHydrogensCoords) = cg
      (Hydrated cdCoords cdHydrogensCoords) = cd
  in case atomType of
    CB   -> cbCoords
    CG   -> cgCoords
    CD   -> cdCoords
    HB2  -> head cbHydrogensCoords
    HB3  -> cbHydrogensCoords !! 1
    HG2  -> head cgHydrogensCoords
    HG3  -> cgHydrogensCoords !! 1
    HD2  -> head cdHydrogensCoords
    HD3  -> cdHydrogensCoords !! 1
    _    -> error $ "No atom type " ++ show atomType ++ "in Proline"
radicalAtomCoordinates (Arginine cb cg cd ne cz nh1 nh2) atomType =
  let (Hydrated cbCoords cbHydrogensCoords) = cb
      (Hydrated cgCoords cgHydrogensCoords) = cg
      (Hydrated cdCoords cdHydrogensCoords) = cd
      (Hydrated neCoords neHydrogensCoords) = ne
      (Hydrated czCoords _                ) = cz
      (Hydrated nh1Coords nh1HydrogensCoords) = nh1
      (Hydrated nh2Coords nh2HydrogensCoords) = nh2
  in case atomType of
    CB   -> cbCoords
    CG   -> cgCoords
    CD   -> cdCoords
    NE   -> neCoords
    CZ   -> czCoords
    NH1  -> nh1Coords
    NH2  -> nh2Coords
    HB2  -> head cbHydrogensCoords
    HB3  -> cbHydrogensCoords !! 1
    HG2 -> head cgHydrogensCoords
    HG3 -> cgHydrogensCoords !! 1
    HD2 -> head cdHydrogensCoords
    HD3 -> cdHydrogensCoords !! 1
    HE  -> head neHydrogensCoords
    HH11 -> head nh1HydrogensCoords
    HH12 -> nh1HydrogensCoords !! 1
    HH21 -> head nh2HydrogensCoords
    HH22 -> nh2HydrogensCoords !! 1
    _    -> error $ "No atom type " ++ show atomType ++ "in Arginine"
radicalAtomCoordinates (Histidine cb cg nd1 ce1 ne2 cd2) atomType =
  let (Hydrated cbCoords cbHydrogensCoords) = cb
      (Hydrated cgCoords _                ) = cg
      (Hydrated nd1Coords nd1HydrogensCoords) = nd1
      (Hydrated ce1Coords ce1HydrogensCoords) = ce1
      (Hydrated ne2Coords ne2HydrogensCoords) = ne2
      (Hydrated cd2Coords cd2HydrogensCoords) = cd2
  in case atomType of
    CB   -> cbCoords
    CG   -> cgCoords
    ND1  -> nd1Coords
    CE1  -> ce1Coords
    NE2  -> ne2Coords
    CD2  -> cd2Coords
    HB2  -> head cbHydrogensCoords
    HB3  -> cbHydrogensCoords !! 1
    HD1  -> head nd1HydrogensCoords
    HE1  -> head ce1HydrogensCoords
    HE2  -> head ne2HydrogensCoords
    HD2  -> head cd2HydrogensCoords
    _    -> error $ "No atom type " ++ show atomType ++ "in Histidine"
radicalAtomCoordinates (Lysine cb cg cd ce nz) atomType =
  let (Hydrated cbCoords cbHydrogensCoords) = cb
      (Hydrated cgCoords cgHydrogensCoords) = cg
      (Hydrated cdCoords cdHydrogensCoords) = cd
      (Hydrated ceCoords ceHydrogensCoords) = ce
      (Hydrated nzCoords nzHydrogensCoords) = nz
  in case atomType of
    CB   -> cbCoords
    CG   -> cgCoords
    CD   -> cdCoords
    CE   -> ceCoords
    NZ   -> nzCoords
    HB2  -> head cbHydrogensCoords
    HB3  -> cbHydrogensCoords !! 1
    HG2  -> head cgHydrogensCoords
    HG3  -> cgHydrogensCoords !! 1
    HD2  -> head cdHydrogensCoords
    HD3  -> cdHydrogensCoords !! 1
    HE2  -> head ceHydrogensCoords
    HE3  -> ceHydrogensCoords !! 1
    HZ1  -> head nzHydrogensCoords
    HZ2  -> nzHydrogensCoords !! 1
    HZ3  -> nzHydrogensCoords !! 2
    _    -> error $ "No atom type " ++ show atomType ++ "in Lysine"
radicalAtomCoordinates (AsparticAcid cb cg od1 od2) atomType =
  let (Hydrated cbCoords cbHydrogensCoords) = cb
      (Hydrated cgCoords _                ) = cg
      (Hydrated od1Coords _               ) = od1
      (Hydrated od2Coords od2HydrogensCoords) = od2
  in case atomType of
    CB   -> cbCoords
    CG   -> cgCoords
    OD1  -> od1Coords
    OD2  -> od2Coords
    HB2  -> head cbHydrogensCoords
    HB3  -> cbHydrogensCoords !! 1
    HD2 -> head od2HydrogensCoords
    _    -> error $ "No atom type " ++ show atomType ++ "in AsparticAcid"
radicalAtomCoordinates (GlutamicAcid cb cg cd oe1 oe2) atomType =
  let (Hydrated cbCoords cbHydrogensCoords) = cb
      (Hydrated cgCoords cgHydrogensCoords) = cg
      (Hydrated cdCoords _                ) = cd
      (Hydrated oe1Coords _               ) = oe1
      (Hydrated oe2Coords oe2HydrogensCoords) = oe2
  in case atomType of
    CB   -> cbCoords
    CG   -> cgCoords
    CD   -> cdCoords
    OE1  -> oe1Coords
    OE2  -> oe2Coords
    HB2  -> head cbHydrogensCoords
    HB3  -> cbHydrogensCoords !! 1
    HG2  -> head cgHydrogensCoords
    HG3  -> cgHydrogensCoords !! 1
    HE2 -> head oe2HydrogensCoords
    _    -> error $ "No atom type " ++ show atomType ++ "in GlutamicAcid"

connections :: HydratedAminoAcid -> [(AtomType, AtomType)]
connections aminoacid = [(N, CA), (CA, C), (C, O)] ++ caConnections aminoacid ++ radicalConnections (radical aminoacid)

caConnections :: HydratedAminoAcid -> [(AtomType, AtomType)]
caConnections (AminoAcid _ _ _ _ _ Glysine) = [(CA, HA2), (CA, HA3)]
caConnections _                             = [(CA, HA)]

radicalConnections :: Radical (Hydrated (V3 Float)) -> [(AtomType, AtomType)]
radicalConnections Alanine {} = [(CA, CB), (CB, HB1), (CB, HB2), (CB, HB3)]
radicalConnections Glysine {} = []
radicalConnections Valine {} = [(CA, CB), (CB, HB), (CB, CG1), (CB, CG2), (CG1, HG11), (CG1, HG12), (CG1, HG13), (CG2, HG21), (CG2, HG22), (CG2, HG23)]
radicalConnections Isoleucine {} = [(CA, CB), (CB, HB), (CB, CG1), (CB, CG2), (CG1, CD1), (CG1, HG12), (CG1, HG13), (CG2, HG21), (CG2, HG22), (CG2, HG23), (CD1, HD11), (CD1, HD12), (CD1, HD12)]
radicalConnections Leucine {} = [(CA, CB), (CB, CG), (CG, CD1), (CG, CD2), (CB, HB2), (CB, HB3), (CG, HG), (CD1, HD11), (CD1, HD12), (CD1, HD13), (CD2, HD21), (CD2, HD22), (CD2, HD23)]
radicalConnections Methionine {} = [(CA, CB), (CB, CG), (CG, SD), (SD, CE), (CB, HB2), (CB, HB3), (CG, HG2), (CG, HG3), (CE, HE1), (CE, HE2), (CE, HE3)]
radicalConnections Phenylalanine {} = [(CA, CB), (CB, CG), (CG, CD1), (CD1, CE1), (CE1, CZ), (CZ, CE2), (CE2, CD2), (CD2, CG), (CB, HB2), (CD, HB3), (CD1, HD1), (CE1, HE1), (CZ, HZ), (CE2, HE2), (CD2, HD2)]
radicalConnections Tyrosine {} = [(CA, CB), (CB, CG), (CG, CD1), (CD1, CE1), (CE1, CZ), (CZ, CE2), (CE2, CD2), (CD2, CG), (CB, HB2), (CD, HB3), (CD1, HD1), (CE1, HE1), (CZ, OH), (CE2, HE2), (CD2, HD2), (OH, HH)]
radicalConnections Tryptophan {} = [(CA, CB), (CB, CG), (CG, CD1), (CD1, NE1), (NE1, CE2), (CE2, CD2), (CD2, CG), (CD2, CE3), (CE3, CZ3), (CZ3, CH2), (CH2, CZ2), (CZ2, CE2), (CB, HB2), (CD, HB3), (CD1, HD1), (NE1, HE1), (CZ2, HZ2), (CH2, HH2), (CZ3, HZ3), (CE3, HE3)]
radicalConnections Serine {} = [(CA, CB), (CB, OG), (CB, HB2), (CB, HB3), (OG, HG)]
radicalConnections Threonine {} = [(CA, CB), (CB, HB), (CB, OG1), (CB, CG2), (OG1, HG1), (CG2, HG21), (CG2, HG22), (CG2, HG23)]
radicalConnections Asparagine {} = [(CA, CB), (CB, HB2), (CB, HB3), (CB, CG), (CG, OD1), (CG, ND2), (ND2, HD21), (ND2, HD22)]
radicalConnections Glutamine {} = [(CA, CB), (CB, HB2), (CB, HB3), (CB, CG), (CG, HG2), (CG, HG3), (CG, CD), (CD, OE1), (CD, NE2), (NE2, HE21), (NE2, HE22)]
radicalConnections Cysteine {} = [(CA, CB), (CB, SG), (CB, HB2), (CB, HB3), (SG, HG)]
radicalConnections Proline {} = [(CA, CB), (CB, CG), (CG, CD), (CD, N), (CB, HB2), (CB, HB3), (CG, HG2), (CG, HG3), (CD, HD2), (CD, HD3)]
radicalConnections Arginine {} = [(CA, CB), (CB, CG), (CG, CD), (CD, NE), (NE, CZ), (CZ, NH2), (CB, HB2), (CB, HB3), (CG, HG2), (CG, HG3), (CD, HD2), (CD, HD3), (NE, HE), (NH1, HH11), (NH1, HH12), (NH2, HH21), (NH2, HH22)]
radicalConnections Histidine {} = [(CA, CB), (CB, CG), (CG, ND1), (ND1, CE1), (CE1, NE2), (NE2, CD2), (CD2, CG), (CB, HB2), (CB, HB3), (ND1, HD1), (CE1, HE1), (NE2, HE2), (CD2, HD2)]
radicalConnections Lysine {} = [(CA, CB), (CB, CG), (CG, CD), (CD, CE), (CE, NZ), (CB, HB2), (CB, HB3), (CG, HG2), (CG, HG3), (CD, HD2), (CD, HD3), (CE, HE2), (CE, HE3), (NZ, HZ1), (NZ, HZ2), (NZ, HZ3)]
radicalConnections AsparticAcid {} = [(CA, CB), (CB, HB2), (CB, HB3), (CB, CG), (CG, OD1), (CG, OD2), (OD2, HD2)]
radicalConnections GlutamicAcid {} = [(CA, CB), (CB, HB2), (CB, HB3), (CB, CG), (CG, HG2), (CG, HG3), (CG, CD), (CD, OE1), (CD, OE2), (OE2, HE2)]
