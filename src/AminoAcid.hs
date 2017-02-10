module AminoAcid
        ( AminoAcid (..)
        , Radical (..)
        , Hydrated (..)
        , HydratedAminoAcid
        , AtomType (..)
        , atomCoordinates
        , atomHydrogens
        , connections
        , atoms
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

data Radical a = Alanine a                      -- CB
               | Glysine                        -- nothing
               | Valine a a a                   -- CB
               | Isoleucine a a a a             -- CB CG1 CG2 CD1
               | Leucine a a a a                -- CB CG CD1 CD2
               | Methionine a a a a             -- CB CG SD CE
               | Phenylalanine a a a a a a a    -- CB CG CD1 CD2 CE1 CE2 CZ
               | Tyrosine a a a a a a a a       -- CB CG CD1 CD2 CE1 CE2 CZ OH
               | Tryptophan a a a a a a a a a a -- CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2
               | Serine a a                     -- CB OG
               | Threonine a a a                -- CB OG1 CG2
               | Asparagine a a a a             -- CB CG OD1 ND2
               | Glutamine a a a a a            -- CB CG CD OE1 NE2
               | Cysteine a a                   -- CB SG
               | Proline a a a                  -- CB CG CD
               | Arginine a a a a a a a         -- CB CG CD NE CZ NH1 NH2
               | Histidine a a a a a a          -- CB CG ND1 CD2 CE1 NE2
               | Lysine a a a a a               -- CB CG CD CE NZ
               | AsparticAcid a a a a           -- CB CG OD1 OD2
               | GlutamicAcid a a a a a         -- CB CG CD OE1 OE2
  deriving Show

data Hydrated a = Hydrated a [a]
  deriving Show

data AtomType = N   | CA  | C   | O  | OXT
              | CB
              | CG  | CG1 | CG2
              | CD  | CD1 | CD2
              | CE  | CE1 | CE2 | CE3
              | CZ  | CZ2 | CZ3
              | CH2
              | SG
              | SD
              | OG  | OG1
              | OD1 | OD2
              | OE1 | OE2
              | OH
              | ND1 | ND2
              | NE  | NE1 | NE2
              | NZ
              | NH1 | NH2
  deriving (Show, Read, Eq)

type HydratedAminoAcid = AminoAcid (Hydrated (V3 Float))

atoms :: HydratedAminoAcid -> [AtomType]
atoms aminoacid = [N, CA, C, O] ++ case radical aminoacid of
  Alanine       {} -> [CB]
  Glysine       {} -> []
  Valine        {} -> [CB]
  Isoleucine    {} -> [CB, CG1, CG2, CD1]
  Leucine       {} -> [CB, CG, CD1, CD2]
  Methionine    {} -> [CB, CG, SD, CE]
  Phenylalanine {} -> [CB, CG, CD1, CD2, CE1, CE2, CZ]
  Tyrosine      {} -> [CB, CG, CD1, CD2, CE1, CE2, CZ, OH]
  Tryptophan    {} -> [CB, CG, CD1, CD2, NE1, CE2, CE3, CZ2, CZ3, CH2]
  Serine        {} -> [CB, OG]
  Threonine     {} -> [CB, OG1, CG2]
  Asparagine    {} -> [CB, CG, OD1, ND2]
  Glutamine     {} -> [CB, CG, CD, OE1, NE2]
  Cysteine      {} -> [CB, SG]
  Proline       {} -> [CB, CG, CD]
  Arginine      {} -> [CB, CG, CD, NE, CZ, NH1, NH2]
  Histidine     {} -> [CB, CG, ND1, CD2, CE1, NE2]
  Lysine        {} -> [CB, CG, CD, CE, NZ]
  AsparticAcid  {} -> [CB, CG, OD1, OD2]
  GlutamicAcid  {} -> [CB, CG, CD, OE1, OE2]

atomCoordinates :: HydratedAminoAcid -> AtomType -> V3 Float
atomCoordinates aminoAcid atomType =
  case atomType of
    N   -> nitroCoords
    CA  -> carbonAlphaCoords
    C   -> carbonCoords
    O   -> oxi2Coords
    OXT -> case oxi aminoAcid of
      Just (Hydrated oxiCoords _) -> oxiCoords
      Nothing -> error "No atom OXT"
    _   -> radicalAtomCoordinates (radical aminoAcid) atomType
  where Hydrated nitroCoords       _ = nitro       aminoAcid
        Hydrated carbonAlphaCoords _ = carbonAlpha aminoAcid
        Hydrated carbonCoords      _ = carbon      aminoAcid
        Hydrated oxi2Coords        _ = oxi2        aminoAcid

radicalAtomCoordinates :: Radical (Hydrated (V3 Float)) -> AtomType -> V3 Float
radicalAtomCoordinates (Alanine (Hydrated atomCoords _)) atomType =
  case atomType of
    CB  -> atomCoords
    _   -> error $ "No atom type " ++ show atomType ++ "in Alanine"
radicalAtomCoordinates Glysine atomType = error $ "No atom type " ++ show atomType ++ "in Glysine, because it hasn't radical"
radicalAtomCoordinates (Valine cb cg1 cg2) atomType =
  let (Hydrated cbCoords  _) = cb
      (Hydrated cg1Coords _) = cg1
      (Hydrated cg2Coords _) = cg2
  in case atomType of
    CB   -> cbCoords
    CG1  -> cg1Coords
    CG2  -> cg2Coords
    _    -> error $ "No atom type " ++ show atomType ++ "in Valine"
radicalAtomCoordinates (Isoleucine cb cg1 cg2 cd1) atomType =
  let (Hydrated cbCoords  _) = cb
      (Hydrated cg1Coords _) = cg1
      (Hydrated cg2Coords _) = cg2
      (Hydrated cd1Coords _) = cd1
  in case atomType of
    CB   -> cbCoords
    CG1  -> cg1Coords
    CG2  -> cg2Coords
    CD1  -> cd1Coords
    _    -> error $ "No atom type " ++ show atomType ++ "in Isoleucine"
radicalAtomCoordinates (Leucine cb cg cd1 cd2) atomType =
  let (Hydrated cbCoords  _) = cb
      (Hydrated cgCoords  _) = cg
      (Hydrated cd1Coords _) = cd1
      (Hydrated cd2Coords _) = cd2
  in case atomType of
    CB   -> cbCoords
    CG   -> cgCoords
    CD1  -> cd1Coords
    CD2  -> cd2Coords
    _    -> error $ "No atom type " ++ show atomType ++ "in Leucine"
radicalAtomCoordinates (Methionine cb cg sd ce) atomType =
  let (Hydrated cbCoords _) = cb
      (Hydrated cgCoords _) = cg
      (Hydrated sdCoords _) = sd
      (Hydrated ceCoords _) = ce
  in case atomType of
    CB   -> cbCoords
    CG   -> cgCoords
    SD   -> sdCoords
    CE   -> ceCoords
    _    -> error $ "No atom type " ++ show atomType ++ "in Methionine"
radicalAtomCoordinates (Phenylalanine cb cg cd1 cd2 ce1 ce2 cz) atomType =
  let (Hydrated cbCoords  _) = cb
      (Hydrated cgCoords  _) = cg
      (Hydrated cd1Coords _) = cd1
      (Hydrated cd2Coords _) = cd2
      (Hydrated ce1Coords _) = ce1
      (Hydrated ce2Coords _) = ce2
      (Hydrated czCoords  _) = cz
  in case atomType of
    CB   -> cbCoords
    CG   -> cgCoords
    CD1  -> cd1Coords
    CE1  -> ce1Coords
    CZ   -> czCoords
    CE2  -> ce2Coords
    CD2  -> cd2Coords
    _    -> error $ "No atom type " ++ show atomType ++ "in Phenylalanine"
radicalAtomCoordinates (Tyrosine cb cg cd1 cd2 ce1 ce2 cz oh) atomType =
  let (Hydrated cbCoords  _) = cb
      (Hydrated cgCoords  _) = cg
      (Hydrated cd1Coords _) = cd1
      (Hydrated cd2Coords _) = cd2
      (Hydrated ce1Coords _) = ce1
      (Hydrated ce2Coords _) = ce2
      (Hydrated czCoords  _) = cz
      (Hydrated ohCoords  _) = oh
  in case atomType of
    CB   -> cbCoords
    CG   -> cgCoords
    CD1  -> cd1Coords
    CE1  -> ce1Coords
    CZ   -> czCoords
    CE2  -> ce2Coords
    CD2  -> cd2Coords
    OH   -> ohCoords
    _    -> error $ "No atom type " ++ show atomType ++ "in Tyrosine"
radicalAtomCoordinates (Tryptophan cb cg cd1 cd2 ne1 ce2 ce3 cz2 cz3 ch2) atomType =
  let (Hydrated cbCoords  _) = cb
      (Hydrated cgCoords  _) = cg
      (Hydrated cd1Coords _) = cd1
      (Hydrated cd2Coords _) = cd2
      (Hydrated ne1Coords _) = ne1
      (Hydrated ce2Coords _) = ce2
      (Hydrated ce3Coords _) = ce3
      (Hydrated cz2Coords _) = cz2
      (Hydrated cz3Coords _) = cz3
      (Hydrated ch2Coords _) = ch2
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
    _    -> error $ "No atom type " ++ show atomType ++ "in Tryptophan"
radicalAtomCoordinates (Serine cb og) atomType =
  let (Hydrated cbCoords _) = cb
      (Hydrated ogCoords _) = og
  in case atomType of
    CB  -> cbCoords
    OG  -> ogCoords
    _   -> error $ "No atom type " ++ show atomType ++ "in Serine"
radicalAtomCoordinates (Threonine cb og1 cg2) atomType =
  let (Hydrated cbCoords  _) = cb
      (Hydrated og1Coords _) = og1
      (Hydrated cg2Coords _) = cg2
  in case atomType of
    CB   -> cbCoords
    OG1  -> og1Coords
    CG2  -> cg2Coords
    _    -> error $ "No atom type " ++ show atomType ++ "in Threonine"
radicalAtomCoordinates (Asparagine cb cg od1 nd2) atomType =
  let (Hydrated cbCoords  _) = cb
      (Hydrated cgCoords  _) = cg
      (Hydrated od1Coords _) = od1
      (Hydrated nd2Coords _) = nd2
  in case atomType of
    CB   -> cbCoords
    CG   -> cgCoords
    OD1  -> od1Coords
    ND2  -> nd2Coords
    _    -> error $ "No atom type " ++ show atomType ++ "in Asparagine"
radicalAtomCoordinates (Glutamine cb cg cd oe1 ne2) atomType =
  let (Hydrated cbCoords  _) = cb
      (Hydrated cgCoords  _) = cg
      (Hydrated cdCoords  _) = cd
      (Hydrated oe1Coords _) = oe1
      (Hydrated ne2Coords _) = ne2
  in case atomType of
    CB   -> cbCoords
    CG   -> cgCoords
    CD   -> cdCoords
    OE1  -> oe1Coords
    NE2  -> ne2Coords
    _    -> error $ "No atom type " ++ show atomType ++ "in Glutamine"
radicalAtomCoordinates (Cysteine cb sg) atomType =
  let (Hydrated cbCoords _) = cb
      (Hydrated sgCoords _) = sg
  in case atomType of
    CB  -> cbCoords
    SG  -> sgCoords
    _   -> error $ "No atom type " ++ show atomType ++ "in Cysteine"
radicalAtomCoordinates (Proline cb cg cd) atomType =
  let (Hydrated cbCoords _) = cb
      (Hydrated cgCoords _) = cg
      (Hydrated cdCoords _) = cd
  in case atomType of
    CB   -> cbCoords
    CG   -> cgCoords
    CD   -> cdCoords
    _    -> error $ "No atom type " ++ show atomType ++ "in Proline"
radicalAtomCoordinates (Arginine cb cg cd ne cz nh1 nh2) atomType =
  let (Hydrated cbCoords  _) = cb
      (Hydrated cgCoords  _) = cg
      (Hydrated cdCoords  _) = cd
      (Hydrated neCoords  _) = ne
      (Hydrated czCoords  _) = cz
      (Hydrated nh1Coords _) = nh1
      (Hydrated nh2Coords _) = nh2
  in case atomType of
    CB   -> cbCoords
    CG   -> cgCoords
    CD   -> cdCoords
    NE   -> neCoords
    CZ   -> czCoords
    NH1  -> nh1Coords
    NH2  -> nh2Coords
    _    -> error $ "No atom type " ++ show atomType ++ "in Arginine"
radicalAtomCoordinates (Histidine cb cg nd1 cd2 ce1 ne2) atomType =
  let (Hydrated cbCoords  _) = cb
      (Hydrated cgCoords  _) = cg
      (Hydrated nd1Coords _) = nd1
      (Hydrated cd2Coords _) = cd2
      (Hydrated ce1Coords _) = ce1
      (Hydrated ne2Coords _) = ne2
  in case atomType of
    CB   -> cbCoords
    CG   -> cgCoords
    ND1  -> nd1Coords
    CE1  -> ce1Coords
    NE2  -> ne2Coords
    CD2  -> cd2Coords
    _    -> error $ "No atom type " ++ show atomType ++ "in Histidine"
radicalAtomCoordinates (Lysine cb cg cd ce nz) atomType =
  let (Hydrated cbCoords _) = cb
      (Hydrated cgCoords _) = cg
      (Hydrated cdCoords _) = cd
      (Hydrated ceCoords _) = ce
      (Hydrated nzCoords _) = nz
  in case atomType of
    CB   -> cbCoords
    CG   -> cgCoords
    CD   -> cdCoords
    CE   -> ceCoords
    NZ   -> nzCoords
    _    -> error $ "No atom type " ++ show atomType ++ "in Lysine"
radicalAtomCoordinates (AsparticAcid cb cg od1 od2) atomType =
  let (Hydrated cbCoords  _) = cb
      (Hydrated cgCoords  _) = cg
      (Hydrated od1Coords _) = od1
      (Hydrated od2Coords _) = od2
  in case atomType of
    CB   -> cbCoords
    CG   -> cgCoords
    OD1  -> od1Coords
    OD2  -> od2Coords
    _    -> error $ "No atom type " ++ show atomType ++ "in AsparticAcid"
radicalAtomCoordinates (GlutamicAcid cb cg cd oe1 oe2) atomType =
  let (Hydrated cbCoords  _) = cb
      (Hydrated cgCoords  _) = cg
      (Hydrated cdCoords  _) = cd
      (Hydrated oe1Coords _) = oe1
      (Hydrated oe2Coords _) = oe2
  in case atomType of
    CB   -> cbCoords
    CG   -> cgCoords
    CD   -> cdCoords
    OE1  -> oe1Coords
    OE2  -> oe2Coords
    _    -> error $ "No atom type " ++ show atomType ++ "in GlutamicAcid"

atomHydrogens :: HydratedAminoAcid -> AtomType -> [V3 Float]
atomHydrogens aminoacid atomType =
  case atomType of
    N   -> nitroHydrogensCoords
    CA  -> carbonAlphaHydrogensCoords
    C   -> carbonHydrogensCoords
    O   -> oxi2HydrogensCoords
    OXT -> case oxi aminoacid of
      Just (Hydrated _ oxiHydrogensCoords) -> oxiHydrogensCoords
      Nothing -> error "No atom OXT"
    _   -> radicalAtomHydrogens (radical aminoacid) atomType
  where Hydrated _ nitroHydrogensCoords       = nitro       aminoacid
        Hydrated _ carbonAlphaHydrogensCoords = carbonAlpha aminoacid
        Hydrated _ carbonHydrogensCoords      = carbon      aminoacid
        Hydrated _ oxi2HydrogensCoords        = oxi2        aminoacid

radicalAtomHydrogens :: Radical (Hydrated (V3 Float)) -> AtomType -> [V3 Float]
radicalAtomHydrogens (Alanine (Hydrated _ hydrogens)) atomType =
  case atomType of
    CB  -> hydrogens
    _   -> error $ "No atom type " ++ show atomType ++ "in Alanine"
radicalAtomHydrogens Glysine atomType = error $ "No atom type " ++ show atomType ++ "in Glysine, because it hasn't radical"
radicalAtomHydrogens (Valine cb cg1 cg2) atomType =
  let (Hydrated _ cbHydrogens ) = cb
      (Hydrated _ cg1Hydrogens) = cg1
      (Hydrated _ cg2Hydrogens) = cg2
  in case atomType of
    CB   -> cbHydrogens
    CG1  -> cg1Hydrogens
    CG2  -> cg2Hydrogens
    _    -> error $ "No atom type " ++ show atomType ++ "in Valine"
radicalAtomHydrogens (Isoleucine cb cg1 cg2 cd1) atomType =
  let (Hydrated _ cbHydrogens ) = cb
      (Hydrated _ cg1Hydrogens) = cg1
      (Hydrated _ cg2Hydrogens) = cg2
      (Hydrated _ cd1Hydrogens) = cd1
  in case atomType of
    CB   -> cbHydrogens
    CG1  -> cg1Hydrogens
    CG2  -> cg2Hydrogens
    CD1  -> cd1Hydrogens
    _    -> error $ "No atom type " ++ show atomType ++ "in Isoleucine"
radicalAtomHydrogens (Leucine cb cg cd1 cd2) atomType =
  let (Hydrated _ cbHydrogens ) = cb
      (Hydrated _ cgHydrogens ) = cg
      (Hydrated _ cd1Hydrogens) = cd1
      (Hydrated _ cd2Hydrogens) = cd2
  in case atomType of
    CB   -> cbHydrogens
    CG   -> cgHydrogens
    CD1  -> cd1Hydrogens
    CD2  -> cd2Hydrogens
    _    -> error $ "No atom type " ++ show atomType ++ "in Leucine"
radicalAtomHydrogens (Methionine cb cg sd ce) atomType =
  let (Hydrated _ cbHydrogens) = cb
      (Hydrated _ cgHydrogens) = cg
      (Hydrated _ sdHydrogens) = sd
      (Hydrated _ ceHydrogens) = ce
  in case atomType of
    CB   -> cbHydrogens
    CG   -> cgHydrogens
    SD   -> sdHydrogens
    CE   -> ceHydrogens
    _    -> error $ "No atom type " ++ show atomType ++ "in Methionine"
radicalAtomHydrogens (Phenylalanine cb cg cd1 cd2 ce1 ce2 cz) atomType =
  let (Hydrated _ cbHydrogens ) = cb
      (Hydrated _ cgHydrogens ) = cg
      (Hydrated _ cd1Hydrogens) = cd1
      (Hydrated _ cd2Hydrogens) = cd2
      (Hydrated _ ce1Hydrogens) = ce1
      (Hydrated _ ce2Hydrogens) = ce2
      (Hydrated _ czHydrogens ) = cz
  in case atomType of
    CB   -> cbHydrogens
    CG   -> cgHydrogens
    CD1  -> cd1Hydrogens
    CE1  -> ce1Hydrogens
    CZ   -> czHydrogens
    CE2  -> ce2Hydrogens
    CD2  -> cd2Hydrogens
    _    -> error $ "No atom type " ++ show atomType ++ "in Phenylalanine"
radicalAtomHydrogens (Tyrosine cb cg cd1 cd2 ce1 ce2 cz oh) atomType =
  let (Hydrated _ cbHydrogens ) = cb
      (Hydrated _ cgHydrogens ) = cg
      (Hydrated _ cd1Hydrogens) = cd1
      (Hydrated _ cd2Hydrogens) = cd2
      (Hydrated _ ce1Hydrogens) = ce1
      (Hydrated _ ce2Hydrogens) = ce2
      (Hydrated _ czHydrogens ) = cz
      (Hydrated _ ohHydrogens ) = oh
  in case atomType of
    CB   -> cbHydrogens
    CG   -> cgHydrogens
    CD1  -> cd1Hydrogens
    CE1  -> ce1Hydrogens
    CZ   -> czHydrogens
    CE2  -> ce2Hydrogens
    CD2  -> cd2Hydrogens
    OH   -> ohHydrogens
    _    -> error $ "No atom type " ++ show atomType ++ "in Tyrosine"
radicalAtomHydrogens (Tryptophan cb cg cd1 cd2 ne1 ce2 ce3 cz2 cz3 ch2) atomType =
  let (Hydrated _ cbHydrogens ) = cb
      (Hydrated _ cgHydrogens ) = cg
      (Hydrated _ cd1Hydrogens) = cd1
      (Hydrated _ cd2Hydrogens) = cd2
      (Hydrated _ ne1Hydrogens) = ne1
      (Hydrated _ ce2Hydrogens) = ce2
      (Hydrated _ ce3Hydrogens) = ce3
      (Hydrated _ cz2Hydrogens) = cz2
      (Hydrated _ cz3Hydrogens) = cz3
      (Hydrated _ ch2Hydrogens) = ch2
  in case atomType of
    CB   -> cbHydrogens
    CG   -> cgHydrogens
    CD1  -> cd1Hydrogens
    NE1  -> ne1Hydrogens
    CE2  -> ce2Hydrogens
    CD2  -> cd2Hydrogens
    CE3  -> ce3Hydrogens
    CZ3  -> cz3Hydrogens
    CH2  -> ch2Hydrogens
    CZ2  -> cz2Hydrogens
    _    -> error $ "No atom type " ++ show atomType ++ "in Tryptophan"
radicalAtomHydrogens (Serine cb og) atomType =
  let (Hydrated _ cbHydrogens) = cb
      (Hydrated _ ogHydrogens) = og
  in case atomType of
    CB  -> cbHydrogens
    OG  -> ogHydrogens
    _   -> error $ "No atom type " ++ show atomType ++ "in Serine"
radicalAtomHydrogens (Threonine cb og1 cg2) atomType =
  let (Hydrated _ cbHydrogens ) = cb
      (Hydrated _ og1Hydrogens) = og1
      (Hydrated _ cg2Hydrogens) = cg2
  in case atomType of
    CB   -> cbHydrogens
    OG1  -> og1Hydrogens
    CG2  -> cg2Hydrogens
    _    -> error $ "No atom type " ++ show atomType ++ "in Threonine"
radicalAtomHydrogens (Asparagine cb cg od1 nd2) atomType =
  let (Hydrated _ cbHydrogens ) = cb
      (Hydrated _ cgHydrogens ) = cg
      (Hydrated _ od1Hydrogens) = od1
      (Hydrated _ nd2Hydrogens) = nd2
  in case atomType of
    CB   -> cbHydrogens
    CG   -> cgHydrogens
    OD1  -> od1Hydrogens
    ND2  -> nd2Hydrogens
    _    -> error $ "No atom type " ++ show atomType ++ "in Asparagine"
radicalAtomHydrogens (Glutamine cb cg cd oe1 ne2) atomType =
  let (Hydrated _ cbHydrogens ) = cb
      (Hydrated _ cgHydrogens ) = cg
      (Hydrated _ cdHydrogens ) = cd
      (Hydrated _ oe1Hydrogens) = oe1
      (Hydrated _ ne2Hydrogens) = ne2
  in case atomType of
    CB   -> cbHydrogens
    CG   -> cgHydrogens
    CD   -> cdHydrogens
    OE1  -> oe1Hydrogens
    NE2  -> ne2Hydrogens
    _    -> error $ "No atom type " ++ show atomType ++ "in Glutamine"
radicalAtomHydrogens (Cysteine cb sg) atomType =
  let (Hydrated _ cbHydrogens) = cb
      (Hydrated _ sgHydrogens) = sg
  in case atomType of
    CB  -> cbHydrogens
    SG  -> sgHydrogens
    _   -> error $ "No atom type " ++ show atomType ++ "in Cysteine"
radicalAtomHydrogens (Proline cb cg cd) atomType =
  let (Hydrated _ cbHydrogens) = cb
      (Hydrated _ cgHydrogens) = cg
      (Hydrated _ cdHydrogens) = cd
  in case atomType of
    CB   -> cbHydrogens
    CG   -> cgHydrogens
    CD   -> cdHydrogens
    _    -> error $ "No atom type " ++ show atomType ++ "in Proline"
radicalAtomHydrogens (Arginine cb cg cd ne cz nh1 nh2) atomType =
  let (Hydrated _ cbHydrogens ) = cb
      (Hydrated _ cgHydrogens ) = cg
      (Hydrated _ cdHydrogens ) = cd
      (Hydrated _ neHydrogens ) = ne
      (Hydrated _ czHydrogens ) = cz
      (Hydrated _ nh1Hydrogens) = nh1
      (Hydrated _ nh2Hydrogens) = nh2
  in case atomType of
    CB   -> cbHydrogens
    CG   -> cgHydrogens
    CD   -> cdHydrogens
    NE   -> neHydrogens
    CZ   -> czHydrogens
    NH1  -> nh1Hydrogens
    NH2  -> nh2Hydrogens
    _    -> error $ "No atom type " ++ show atomType ++ "in Arginine"
radicalAtomHydrogens (Histidine cb cg nd1 cd2 ce1 ne2) atomType =
  let (Hydrated _ cbHydrogens ) = cb
      (Hydrated _ cgHydrogens ) = cg
      (Hydrated _ nd1Hydrogens) = nd1
      (Hydrated _ cd2Hydrogens) = cd2
      (Hydrated _ ce1Hydrogens) = ce1
      (Hydrated _ ne2Hydrogens) = ne2
  in case atomType of
    CB   -> cbHydrogens
    CG   -> cgHydrogens
    ND1  -> nd1Hydrogens
    CE1  -> ce1Hydrogens
    NE2  -> ne2Hydrogens
    CD2  -> cd2Hydrogens
    _    -> error $ "No atom type " ++ show atomType ++ "in Histidine"
radicalAtomHydrogens (Lysine cb cg cd ce nz) atomType =
  let (Hydrated _ cbHydrogens) = cb
      (Hydrated _ cgHydrogens) = cg
      (Hydrated _ cdHydrogens) = cd
      (Hydrated _ ceHydrogens) = ce
      (Hydrated _ nzHydrogens) = nz
  in case atomType of
    CB   -> cbHydrogens
    CG   -> cgHydrogens
    CD   -> cdHydrogens
    CE   -> ceHydrogens
    NZ   -> nzHydrogens
    _    -> error $ "No atom type " ++ show atomType ++ "in Lysine"
radicalAtomHydrogens (AsparticAcid cb cg od1 od2) atomType =
  let (Hydrated _ cbHydrogens ) = cb
      (Hydrated _ cgHydrogens ) = cg
      (Hydrated _ od1Hydrogens) = od1
      (Hydrated _ od2Hydrogens) = od2
  in case atomType of
    CB   -> cbHydrogens
    CG   -> cgHydrogens
    OD1  -> od1Hydrogens
    OD2  -> od2Hydrogens
    _    -> error $ "No atom type " ++ show atomType ++ "in AsparticAcid"
radicalAtomHydrogens (GlutamicAcid cb cg cd oe1 oe2) atomType =
  let (Hydrated _ cbHydrogens ) = cb
      (Hydrated _ cgHydrogens ) = cg
      (Hydrated _ cdHydrogens ) = cd
      (Hydrated _ oe1Hydrogens) = oe1
      (Hydrated _ oe2Hydrogens) = oe2
  in case atomType of
    CB   -> cbHydrogens
    CG   -> cgHydrogens
    CD   -> cdHydrogens
    OE1  -> oe1Hydrogens
    OE2  -> oe2Hydrogens
    _    -> error $ "No atom type " ++ show atomType ++ "in GlutamicAcid"

connections :: HydratedAminoAcid -> [(AtomType, AtomType)]
connections aminoacid = [(N, CA), (CA, C), (C, O)] ++ caConnections aminoacid ++ radicalConnections (radical aminoacid)

caConnections :: HydratedAminoAcid -> [(AtomType, AtomType)]
caConnections aminoacid = case radical aminoacid of
  Glysine {} -> []
  _          -> [(CA, CB)]

radicalConnections :: Radical (Hydrated (V3 Float)) -> [(AtomType, AtomType)]
radicalConnections Alanine {} = []
radicalConnections Glysine {} = []
radicalConnections Valine {} = [(CB, CG1), (CB, CG2)]
radicalConnections Isoleucine {} = [(CB, CG1), (CB, CG2), (CG1, CD1)]
radicalConnections Leucine {} = [(CB, CG), (CG, CD1), (CG, CD2)]
radicalConnections Methionine {} = [(CB, CG), (CG, SD), (SD, CE)]
radicalConnections Phenylalanine {} = [(CB, CG), (CG, CD1), (CD1, CE1), (CE1, CZ), (CZ, CE2), (CE2, CD2), (CD2, CG)]
radicalConnections Tyrosine {} = [(CB, CG), (CG, CD1), (CD1, CE1), (CE1, CZ), (CZ, CE2), (CE2, CD2), (CD2, CG)]
radicalConnections Tryptophan {} = [(CB, CG), (CG, CD1), (CD1, NE1), (NE1, CE2), (CE2, CD2), (CD2, CG), (CD2, CE3), (CE3, CZ3), (CZ3, CH2), (CH2, CZ2), (CZ2, CE2)]
radicalConnections Serine {} = [(CB, OG)]
radicalConnections Threonine {} = [(CB, OG1), (CB, CG2)]
radicalConnections Asparagine {} = [(CB, CG), (CG, OD1), (CG, ND2)]
radicalConnections Glutamine {} = [(CB, CG), (CG, CD), (CD, OE1), (CD, NE2)]
radicalConnections Cysteine {} = [(CB, SG)]
radicalConnections Proline {} = [(CB, CG), (CG, CD), (CD, N)]
radicalConnections Arginine {} = [(CB, CG), (CG, CD), (CD, NE), (NE, CZ), (CZ, NH2)]
radicalConnections Histidine {} = [(CB, CG), (CG, ND1), (ND1, CE1), (CE1, NE2), (NE2, CD2), (CD2, CG)]
radicalConnections Lysine {} = [(CB, CG), (CG, CD), (CD, CE), (CE, NZ)]
radicalConnections AsparticAcid {} = [(CB, CG), (CG, OD1), (CG, OD2)]
radicalConnections GlutamicAcid {} = [(CB, CG), (CG, CD), (CD, OE1), (CD, OE2)]
