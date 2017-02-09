module AminoAcidPDB
        ( pdbAminoParser

        ) where

import Linear.V3 (V3 (..))
import AminoAcid (AminoAcid (..), Radical (..), Hydrated (..), HydratedAminoAcid, AtomType (..))
import Parser.PDB (PDBLine (..), pdbLinesParser)
import Text.Parsec.Text (Parser)
import Data.Text (unpack)
import Data.Maybe (fromMaybe, catMaybes)
import Debug.Trace (trace, traceShow)

pdbAminoParser :: Parser [HydratedAminoAcid]
pdbAminoParser = pdbLinesToAminos <$> pdbLinesParser

pdbLinesToAminos :: [PDBLine] -> [HydratedAminoAcid]
pdbLinesToAminos list = reverse $ pdbLinesToAmino <$> groupPDBLinesByAminos' list1
  where list1 = filter atomLineFilter list

atomLineFilter :: PDBLine -> Bool
atomLineFilter AtomLine {} = True
atomLineFilter _ = False

groupPDBLinesByAminos' :: [PDBLine] -> [[PDBLine]]
groupPDBLinesByAminos' = groupPDBLinesByAminos []

groupPDBLinesByAminos :: [[PDBLine]] -> [PDBLine] -> [[PDBLine]]
groupPDBLinesByAminos res []       = res
groupPDBLinesByAminos [] (x : xs)  = groupPDBLinesByAminos [[x]] xs
groupPDBLinesByAminos res@(lastList : resultTail) (x : xs)
  | atomResSeq x == atomResSeq (head lastList) = groupPDBLinesByAminos ((x : lastList) : resultTail) xs
  | otherwise                                  = groupPDBLinesByAminos ([x] : res) xs

pdbLineToAtom :: PDBLine -> (AtomType, V3 Float)
pdbLineToAtom line = (atomType, coordinates)
  where atomType = read $ unpack $ atomName line
        coordinates = V3 (atomX line) (atomY line) (atomZ line)

pdbLinesToAmino :: [PDBLine] -> HydratedAminoAcid
pdbLinesToAmino list =
  AminoAcid { nitro       = Hydrated (atomByName N) (findHydrogens [H, H1, H2, H3] convertedList)
            , carbonAlpha = Hydrated (atomByName CA) (findHydrogens [HA, HA2, HA3] convertedList) 
            , carbon      = Hydrated (atomByName C)  []
            , oxi2        = Hydrated (atomByName O)  []
            , oxi         = Nothing -- TODO: FIX ME
            , radical     = pdbLinesToRadical (unpack . atomResName $ head list) convertedList
            }
  where convertedList = pdbLineToAtom <$> list
        atomByName = findAtom convertedList

findHydrogens :: [AtomType] -> [(AtomType, V3 Float)] -> [V3 Float]
findHydrogens hydrogens atoms = catMaybes $ (`lookup` atoms) <$> hydrogens
-- findNitroHydrogens atoms = case lookup H atoms of
--   Just coords -> [coords]
--   Nothing     -> [findAtom atoms H1, findAtom atoms H2, findAtom atoms H3]

findAtom :: [(AtomType, V3 Float)] -> AtomType -> V3 Float
findAtom atoms atomType = fromMaybe
  (error $ "No atom " ++ show atomType ++ " in " ++ show (fst <$> atoms))
  (lookup atomType atoms)

pdbLinesToRadical :: String -> [(AtomType, V3 Float)] -> Radical (Hydrated (V3 Float))
pdbLinesToRadical "ALA" atoms = Alanine (Hydrated (findAtom atoms CB) [findAtom atoms HB1, findAtom atoms HB2, findAtom atoms HB3])
pdbLinesToRadical "GLY" _     = Glysine
pdbLinesToRadical "VAL" atoms = Valine
  (Hydrated (findAtom atoms CB) [findAtom atoms HB])
  (Hydrated (findAtom atoms CG1) [findAtom atoms HG11, findAtom atoms HG12, findAtom atoms HG13])
  (Hydrated (findAtom atoms CG2) [findAtom atoms HG21, findAtom atoms HG22, findAtom atoms HG23])
pdbLinesToRadical "ILE" atoms = Isoleucine
  (Hydrated (findAtom atoms CB) [findAtom atoms HB])
  (Hydrated (findAtom atoms CG1) [findAtom atoms HG12, findAtom atoms HG13])
  (Hydrated (findAtom atoms CG2) [findAtom atoms HG21, findAtom atoms HG22, findAtom atoms HG23])
  (Hydrated (findAtom atoms CD1) [findAtom atoms HD11, findAtom atoms HD12, findAtom atoms HD13])
pdbLinesToRadical "LEU" atoms = Leucine
  (Hydrated (findAtom atoms CB) [findAtom atoms HB2, findAtom atoms HB3])
  (Hydrated (findAtom atoms CG) [findAtom atoms HG])
  (Hydrated (findAtom atoms CD1) [findAtom atoms HD11, findAtom atoms HD12, findAtom atoms HD13])
  (Hydrated (findAtom atoms CD2) [findAtom atoms HD21, findAtom atoms HD22, findAtom atoms HD23])
pdbLinesToRadical "MET" atoms = Methionine
  (Hydrated (findAtom atoms CB) [findAtom atoms HB2, findAtom atoms HB3])
  (Hydrated (findAtom atoms CG) [findAtom atoms HG2, findAtom atoms HG3])
  (Hydrated (findAtom atoms SD) [])
  (Hydrated (findAtom atoms CE) [findAtom atoms HE1, findAtom atoms HE2, findAtom atoms HE3])
pdbLinesToRadical "PHE" atoms = Phenylalanine
  (Hydrated (findAtom atoms CB) [findAtom atoms HB2, findAtom atoms HB3])
  (Hydrated (findAtom atoms CG) [])
  (Hydrated (findAtom atoms CD1) [findAtom atoms HD1])
  (Hydrated (findAtom atoms CE1) [findAtom atoms HE1])
  (Hydrated (findAtom atoms CZ) [findAtom atoms HZ])
  (Hydrated (findAtom atoms CE2) [findAtom atoms HE2])
  (Hydrated (findAtom atoms CD2) [findAtom atoms HD2])
pdbLinesToRadical "TYR" atoms = Tyrosine
  (Hydrated (findAtom atoms CB) [findAtom atoms HB2, findAtom atoms HB3])
  (Hydrated (findAtom atoms CG) [])
  (Hydrated (findAtom atoms CD1) [findAtom atoms HD1])
  (Hydrated (findAtom atoms CE1) [findAtom atoms HE1])
  (Hydrated (findAtom atoms CZ) [findAtom atoms HZ])
  (Hydrated (findAtom atoms CE2) [findAtom atoms HE2])
  (Hydrated (findAtom atoms CD2) [findAtom atoms HD2])
  (Hydrated (findAtom atoms OH) [findAtom atoms HH])
pdbLinesToRadical "TRP" atoms = Tryptophan
  (Hydrated (findAtom atoms CB) [findAtom atoms HB2, findAtom atoms HB3])
  (Hydrated (findAtom atoms CG) [])
  (Hydrated (findAtom atoms CD1) [findAtom atoms HD1])
  (Hydrated (findAtom atoms NE1) [findAtom atoms HE1])
  (Hydrated (findAtom atoms CE2) [])
  (Hydrated (findAtom atoms CD2) [])
  (Hydrated (findAtom atoms CE3) [findAtom atoms HE3])
  (Hydrated (findAtom atoms CZ3) [findAtom atoms HZ3])
  (Hydrated (findAtom atoms CH2) [findAtom atoms HH2])
  (Hydrated (findAtom atoms CZ2) [findAtom atoms HZ2])
pdbLinesToRadical "SER" atoms = Serine
  (Hydrated (findAtom atoms CB) [findAtom atoms HB2, findAtom atoms HB3])
  (Hydrated (findAtom atoms OG) [findAtom atoms HG])
pdbLinesToRadical "THR" atoms = Threonine
  (Hydrated (findAtom atoms CB) [findAtom atoms HB])
  (Hydrated (findAtom atoms OG1) [findAtom atoms HG1])
  (Hydrated (findAtom atoms CG2) [findAtom atoms HG21, findAtom atoms HG22, findAtom atoms HG23])
pdbLinesToRadical "ASN" atoms = Asparagine
  (Hydrated (findAtom atoms CB) [findAtom atoms HB2, findAtom atoms HB3])
  (Hydrated (findAtom atoms CG) [])
  (Hydrated (findAtom atoms OD1) [])
  (Hydrated (findAtom atoms ND2) [findAtom atoms HD21, findAtom atoms HD22])
pdbLinesToRadical "GLN" atoms = Glutamine
  (Hydrated (findAtom atoms CB) [findAtom atoms HB2, findAtom atoms HB3])
  (Hydrated (findAtom atoms CG) [findAtom atoms HG2, findAtom atoms HG3])
  (Hydrated (findAtom atoms CD) [])
  (Hydrated (findAtom atoms OE1) [])
  (Hydrated (findAtom atoms NE2) [findAtom atoms HE21, findAtom atoms HE22])
pdbLinesToRadical "CYS" atoms = Cysteine
  (Hydrated (findAtom atoms CB) [findAtom atoms HB2, findAtom atoms HB3])
  (Hydrated (findAtom atoms SG) [findAtom atoms HG])
pdbLinesToRadical "PRO" atoms = Proline
  (Hydrated (findAtom atoms CB) [findAtom atoms HB2, findAtom atoms HB3])
  (Hydrated (findAtom atoms CG) [findAtom atoms HG2, findAtom atoms HG3])
  (Hydrated (findAtom atoms CD) [findAtom atoms HD2, findAtom atoms HD3])
pdbLinesToRadical "ARG" atoms = Arginine
  (Hydrated (findAtom atoms CB) [findAtom atoms HB2, findAtom atoms HB3])
  (Hydrated (findAtom atoms CG) [findAtom atoms HG2, findAtom atoms HG3])
  (Hydrated (findAtom atoms CD) [findAtom atoms HD2, findAtom atoms HD3])
  (Hydrated (findAtom atoms NE) [findAtom atoms HE])
  (Hydrated (findAtom atoms CZ) [])
  (Hydrated (findAtom atoms NH1) [findAtom atoms HH11, findAtom atoms HH12])
  (Hydrated (findAtom atoms NH2) [findAtom atoms HH21, findAtom atoms HH22])
pdbLinesToRadical "HIS" atoms = Histidine
  (Hydrated (findAtom atoms CB) [findAtom atoms HB2, findAtom atoms HB3])
  (Hydrated (findAtom atoms CG) [])
  (Hydrated (findAtom atoms ND1) [findAtom atoms HD1])
  (Hydrated (findAtom atoms CE1) [findAtom atoms HE1, findAtom atoms HE2])
  (Hydrated (findAtom atoms NE2) [findAtom atoms HE2])
  (Hydrated (findAtom atoms CD2) [findAtom atoms HD2])
pdbLinesToRadical "LYS" atoms = Lysine
  (Hydrated (findAtom atoms CB) [findAtom atoms HB2, findAtom atoms HB3])
  (Hydrated (findAtom atoms CG) [findAtom atoms HG2, findAtom atoms HG3])
  (Hydrated (findAtom atoms CD) [findAtom atoms HD2, findAtom atoms HD3])
  (Hydrated (findAtom atoms CE) [findAtom atoms HE2, findAtom atoms HE3])
  (Hydrated (findAtom atoms NZ) [findAtom atoms HZ1, findAtom atoms HZ2, findAtom atoms HZ3])
pdbLinesToRadical "ASP" atoms = AsparticAcid
  (Hydrated (findAtom atoms CB) [findAtom atoms HB2, findAtom atoms HB3])
  (Hydrated (findAtom atoms CG) [])
  (Hydrated (findAtom atoms OD1) [])
  (Hydrated (findAtom atoms OD2) [findAtom atoms HD2])
pdbLinesToRadical "GLU" atoms = GlutamicAcid
  (Hydrated (findAtom atoms CB) [findAtom atoms HB2, findAtom atoms HB3])
  (Hydrated (findAtom atoms CG) [findAtom atoms HG2, findAtom atoms HG3])
  (Hydrated (findAtom atoms CD) [])
  (Hydrated (findAtom atoms OE1) [])
  (Hydrated (findAtom atoms OE2) [findAtom atoms HE2])
pdbLinesToRadical amino _ = error $ "Unknown amino acid code => " ++ amino
