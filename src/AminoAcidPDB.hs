module AminoAcidPDB
        ( pdbAminoParser
        ) where

import           AminoAcid        (AminoAcid (..), AtomType (..), Hydrated (..),
                                   HydratedAminoAcid, Radical (..))
import           Data.Map         (Map, fromList, keys, lookup)
import           Data.Maybe       (catMaybes, fromMaybe)
import           Data.Text        (unpack)
import           Linear.V3        (V3 (..))
import           Parser.PDB       (PDBLine (..), pdbLinesParser)
import           Prelude          hiding (lookup)
import           Text.Parsec.Text (Parser)

pdbAminoParser :: Parser [HydratedAminoAcid]
pdbAminoParser = pdbLinesToAminos <$> pdbLinesParser

pdbLinesToAminos :: [PDBLine] -> [HydratedAminoAcid]
pdbLinesToAminos list = reverse $ pdbLinesToAmino <$> groupPDBLinesByAminos atomLines
  where atomLines = filter onlyAtoms list

pdbLinesToAmino :: [PDBLine] -> HydratedAminoAcid
pdbLinesToAmino list =
  AminoAcid { nitro       = hydratedFrom atomsMap N
            , carbonAlpha = hydratedFrom atomsMap CA
            , carbon      = hydratedFrom atomsMap C
            , oxi2        = hydratedFrom atomsMap O
            , oxi         = Hydrated <$> lookup (show OXT) atomsMap <*> pure (findHydrogensForAtom atomsMap OXT)
            , radical     = pdbLinesToRadical (unpack . atomResName $ head list) atomsMap
            }
  where atomsMap = fromList $ pdbLineToAtom <$> list

pdbLinesToRadical :: String -> Map String (V3 Float) -> Radical (Hydrated (V3 Float))
pdbLinesToRadical "ALA" atoms = hydratedFrom atoms <$> Alanine CB
pdbLinesToRadical "GLY" _     = Glysine
pdbLinesToRadical "VAL" atoms = hydratedFrom atoms <$> Valine CB CG1 CG2
pdbLinesToRadical "ILE" atoms = hydratedFrom atoms <$> Isoleucine CB CG1 CG2 CD1
pdbLinesToRadical "LEU" atoms = hydratedFrom atoms <$> Leucine CB CG CD1 CD2
pdbLinesToRadical "MET" atoms = hydratedFrom atoms <$> Methionine CB CG SD CE
pdbLinesToRadical "PHE" atoms = hydratedFrom atoms <$> Phenylalanine CB CG CD1 CD2 CE1 CE2 CZ
pdbLinesToRadical "TYR" atoms = hydratedFrom atoms <$> Tyrosine CB CG CD1 CD2 CE1 CE2 CZ OH
pdbLinesToRadical "TRP" atoms = hydratedFrom atoms <$> Tryptophan CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2
pdbLinesToRadical "SER" atoms = hydratedFrom atoms <$> Serine CB OG
pdbLinesToRadical "THR" atoms = hydratedFrom atoms <$> Threonine CB OG1 CG2
pdbLinesToRadical "ASN" atoms = hydratedFrom atoms <$> Asparagine CB CG OD1 ND2
pdbLinesToRadical "GLN" atoms = hydratedFrom atoms <$> Glutamine CB CG CD OE1 NE2
pdbLinesToRadical "CYS" atoms = hydratedFrom atoms <$> Cysteine CB SG
pdbLinesToRadical "PRO" atoms = hydratedFrom atoms <$> Proline CB CG CD
pdbLinesToRadical "ARG" atoms = hydratedFrom atoms <$> Arginine CB CG CD NE CZ NH1 NH2
pdbLinesToRadical "HIS" atoms = hydratedFrom atoms <$> Histidine CB CG ND1 CE1 NE2 CD2
pdbLinesToRadical "LYS" atoms = hydratedFrom atoms <$> Lysine CB CG CD CE NZ
pdbLinesToRadical "ASP" atoms = hydratedFrom atoms <$> AsparticAcid CB CG OD1 OD2
pdbLinesToRadical "GLU" atoms = hydratedFrom atoms <$> GlutamicAcid CB CG CD OE1 OE2
pdbLinesToRadical amino _ = error $ "Unknown amino acid code => " ++ amino

hydratedFrom ::  Map String (V3 Float) -> AtomType -> Hydrated (V3 Float)
hydratedFrom atoms C = Hydrated (findAtom atoms C) []
hydratedFrom atoms O = Hydrated (findAtom atoms O) []
hydratedFrom atoms atomType = Hydrated (findAtom atoms atomType) (findHydrogensForAtom atoms atomType)

findHydrogensForAtom :: Map String (V3 Float) -> AtomType -> [V3 Float]
findHydrogensForAtom atoms atomType = findHydrogens atoms hydrogens
  where (_:atomCode) = show atomType
        hydrogens    = (('H' : atomCode) ++) <$> ["", "1", "2", "3"]

findHydrogens :: Map String (V3 Float) -> [String] -> [V3 Float]
findHydrogens atoms hydrogens = catMaybes $ (`lookup` atoms) <$> hydrogens

findAtom :: Map String (V3 Float) -> AtomType -> V3 Float
findAtom atoms atomType = fromMaybe
  (error $ "No atom " ++ show atomType ++ " in " ++ show (keys atoms))
  (lookup (show atomType) atoms)

onlyAtoms :: PDBLine -> Bool
onlyAtoms AtomLine {} = True
onlyAtoms _ = False

groupPDBLinesByAminos :: [PDBLine] -> [[PDBLine]]
groupPDBLinesByAminos = groupPDBLinesByAminosHelper []

groupPDBLinesByAminosHelper :: [[PDBLine]] -> [PDBLine] -> [[PDBLine]]
groupPDBLinesByAminosHelper res []       = res
groupPDBLinesByAminosHelper [] (x : xs)  = groupPDBLinesByAminosHelper [[x]] xs
groupPDBLinesByAminosHelper res@(lastList : resultTail) (x : xs)
  | atomResSeq x == atomResSeq (head lastList) = groupPDBLinesByAminosHelper ((x : lastList) : resultTail) xs
  | otherwise                                  = groupPDBLinesByAminosHelper ([x] : res) xs

pdbLineToAtom :: PDBLine -> (String, V3 Float)
pdbLineToAtom line = (atomType, coordinates)
  where atomType    = unpack $ atomName line
        coordinates = V3 (atomX line) (atomY line) (atomZ line)
