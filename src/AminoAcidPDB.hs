module AminoAcidPDB
        ( pdbAminoParser

        ) where

import Linear.V3 (V3 (..))
import AminoAcid (AminoAcid (..), Radical (..), Hydrated (..), HydratedAminoAcid, AtomType (..))
import Parser.PDB (PDBLine (..), pdbLinesParser)
import Text.Parsec.Text (Parser)
import Data.Text (unpack)
import Data.Maybe (fromMaybe, catMaybes)
import Debug.Trace (trace)

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

pdbLineToAtom :: PDBLine -> (String, V3 Float)
pdbLineToAtom line = (atomType, coordinates)
  where atomType = unpack $ atomName line
        coordinates = V3 (atomX line) (atomY line) (atomZ line)

pdbLinesToAmino :: [PDBLine] -> HydratedAminoAcid
pdbLinesToAmino list =
  AminoAcid { nitro       = createHydrated convertedList N
            , carbonAlpha = createHydrated convertedList CA
            , carbon      = createHydrated convertedList C
            , oxi2        = createHydrated convertedList O
            , oxi         = Nothing -- TODO: FIX ME
            , radical     = pdbLinesToRadical (unpack . atomResName $ head list) convertedList
            }
  where convertedList = pdbLineToAtom <$> list

createHydrated ::  [(String, V3 Float)] -> AtomType -> Hydrated (V3 Float)
createHydrated atoms C = Hydrated (findAtom atoms C) []
createHydrated atoms O = Hydrated (findAtom atoms O) []
createHydrated atoms atomType = Hydrated (findAtom atoms atomType) (findHydrogensForAtom atoms atomType)

findHydrogensForAtom :: [(String, V3 Float)] -> AtomType -> [V3 Float]
findHydrogensForAtom atoms atomType = trace (show hydrogens) (findHydrogens atoms hydrogens)
  where (_:atomCode) = show atomType
        hydrogens    = (('H' : atomCode) ++) <$> ["", "1", "2", "3"]

findHydrogens :: [(String, V3 Float)] -> [String] -> [V3 Float]
findHydrogens atoms hydrogens = catMaybes $ (`lookup` atoms) <$> hydrogens

findAtom :: [(String, V3 Float)] -> AtomType -> V3 Float
findAtom atoms atomType = fromMaybe
  (error $ "No atom " ++ show atomType ++ " in " ++ show (fst <$> atoms))
  (lookup (show atomType) atoms)

pdbLinesToRadical :: String -> [(String, V3 Float)] -> Radical (Hydrated (V3 Float))
pdbLinesToRadical "ALA" atoms = createHydrated atoms <$> Alanine CB
pdbLinesToRadical "GLY" _     = Glysine
pdbLinesToRadical "VAL" atoms = createHydrated atoms <$> Valine CB CG1 CG2
pdbLinesToRadical "ILE" atoms = createHydrated atoms <$> Isoleucine CB CG1 CG2 CD1
pdbLinesToRadical "LEU" atoms = createHydrated atoms <$> Leucine CB CG CD1 CD2
pdbLinesToRadical "MET" atoms = createHydrated atoms <$> Methionine CB CG SD CE
pdbLinesToRadical "PHE" atoms = createHydrated atoms <$> Phenylalanine CB CG CD1 CD2 CE1 CE2 CZ
pdbLinesToRadical "TYR" atoms = createHydrated atoms <$> Tyrosine CB CG CD1 CD2 CE1 CE2 CZ OH
pdbLinesToRadical "TRP" atoms = createHydrated atoms <$> Tryptophan CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2
pdbLinesToRadical "SER" atoms = createHydrated atoms <$> Serine CB OG
pdbLinesToRadical "THR" atoms = createHydrated atoms <$> Threonine CB OG1 CG2
pdbLinesToRadical "ASN" atoms = createHydrated atoms <$> Asparagine CB CG OD1 ND2
pdbLinesToRadical "GLN" atoms = createHydrated atoms <$> Glutamine CB CG CD OE1 NE2
pdbLinesToRadical "CYS" atoms = createHydrated atoms <$> Cysteine CB SG
pdbLinesToRadical "PRO" atoms = createHydrated atoms <$> Proline CB CG CD
pdbLinesToRadical "ARG" atoms = createHydrated atoms <$> Arginine CB CG CD NE CZ NH1 NH2
pdbLinesToRadical "HIS" atoms = createHydrated atoms <$> Histidine CB CG ND1 CE1 NE2 CD2
pdbLinesToRadical "LYS" atoms = createHydrated atoms <$> Lysine CB CG CD CE NZ
pdbLinesToRadical "ASP" atoms = createHydrated atoms <$> AsparticAcid CB CG OD1 OD2
pdbLinesToRadical "GLU" atoms = createHydrated atoms <$> GlutamicAcid CB CG CD OE1 OE2
pdbLinesToRadical amino _ = error $ "Unknown amino acid code => " ++ amino
