module AminoAcidPDB
        ( pdbAminoParser

        ) where

import Linear.V3 (V3 (..))
import AminoAcid (AminoAcid (..), Radical (..), Hydrated (..), HydratedAminoAcid, AtomType (..))
import Parser.PDB (PDBLine (..), pdbLinesParser)
import Text.Parsec.Text (Parser)
import Data.Text (unpack)
import Data.Maybe (fromMaybe, catMaybes)

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
  AminoAcid { nitro       = createHydrated N convertedList
            , carbonAlpha = createHydrated CA convertedList
            , carbon      = createHydrated C convertedList
            , oxi2        = createHydrated O convertedList
            , oxi         = Nothing -- TODO: FIX ME
            , radical     = pdbLinesToRadical (unpack . atomResName $ head list) convertedList
            }
  where convertedList = pdbLineToAtom <$> list

createHydrated :: AtomType -> [(String, V3 Float)] -> Hydrated (V3 Float)
createHydrated C atoms = Hydrated (findAtom atoms C) []
createHydrated O atoms = Hydrated (findAtom atoms O) []
createHydrated atomType atoms = Hydrated (findAtom atoms atomType) (findHydrogensForAtom atoms atomType)

findHydrogensForAtom :: [(String, V3 Float)] -> AtomType -> [V3 Float]
findHydrogensForAtom atoms atomType = findHydrogens hydrogens atoms
  where (_:atomCode) = show atomType
        hydrogens    = (atomCode ++) <$> ["", "1", "2", "3"]

findHydrogens :: [String] -> [(String, V3 Float)] -> [V3 Float]
findHydrogens hydrogens atoms = catMaybes $ (`lookup` atoms) <$> hydrogens

findAtom :: [(String, V3 Float)] -> AtomType -> V3 Float
findAtom atoms atomType = fromMaybe
  (error $ "No atom " ++ show atomType ++ " in " ++ show (fst <$> atoms))
  (lookup (show atomType) atoms)

pdbLinesToRadical :: String -> [(String, V3 Float)] -> Radical (Hydrated (V3 Float))
pdbLinesToRadical "ALA" atoms = Alanine (createHydrated CB atoms)
pdbLinesToRadical "GLY" _     = Glysine
pdbLinesToRadical "VAL" atoms = Valine
  (createHydrated CB atoms)
  (createHydrated CG1 atoms)
  (createHydrated CG2 atoms)
pdbLinesToRadical "ILE" atoms = Isoleucine
  (createHydrated CB atoms)
  (createHydrated CG1 atoms)
  (createHydrated CG2 atoms)
  (createHydrated CD1 atoms)
pdbLinesToRadical "LEU" atoms = Leucine
  (createHydrated CB atoms)
  (createHydrated CG atoms)
  (createHydrated CD1 atoms)
  (createHydrated CD2 atoms)
pdbLinesToRadical "MET" atoms = Methionine
  (createHydrated CB atoms)
  (createHydrated CG atoms)
  (createHydrated SD atoms)
  (createHydrated CE atoms)
pdbLinesToRadical "PHE" atoms = Phenylalanine
  (createHydrated CB atoms)
  (createHydrated CG atoms)
  (createHydrated CD1 atoms)
  (createHydrated CD2 atoms)
  (createHydrated CE1 atoms)
  (createHydrated CE2 atoms)
  (createHydrated CZ atoms)
pdbLinesToRadical "TYR" atoms = Tyrosine
  (createHydrated CB atoms)
  (createHydrated CG atoms)
  (createHydrated CD1 atoms)
  (createHydrated CD2 atoms)
  (createHydrated CE1 atoms)
  (createHydrated CE2 atoms)
  (createHydrated CZ atoms)
  (createHydrated OH atoms)
pdbLinesToRadical "TRP" atoms = Tryptophan
  (createHydrated CB atoms)
  (createHydrated CG atoms)
  (createHydrated CD1 atoms)
  (createHydrated CD2 atoms)
  (createHydrated NE1 atoms)
  (createHydrated CE2 atoms)
  (createHydrated CE3 atoms)
  (createHydrated CZ2 atoms)
  (createHydrated CZ3 atoms)
  (createHydrated CH2 atoms)
pdbLinesToRadical "SER" atoms = Serine
  (createHydrated CB atoms)
  (createHydrated OG atoms)
pdbLinesToRadical "THR" atoms = Threonine
  (createHydrated CB atoms)
  (createHydrated OG1 atoms)
  (createHydrated CG2 atoms)
pdbLinesToRadical "ASN" atoms = Asparagine
  (createHydrated CB atoms)
  (createHydrated CG atoms)
  (createHydrated OD1 atoms)
  (createHydrated ND2 atoms)
pdbLinesToRadical "GLN" atoms = Glutamine
  (createHydrated CB atoms)
  (createHydrated CG atoms)
  (createHydrated CD atoms)
  (createHydrated OE1 atoms)
  (createHydrated NE2 atoms)
pdbLinesToRadical "CYS" atoms = Cysteine
  (createHydrated CB atoms)
  (createHydrated SG atoms)
pdbLinesToRadical "PRO" atoms = Proline
  (createHydrated CB atoms)
  (createHydrated CG atoms)
  (createHydrated CD atoms)
pdbLinesToRadical "ARG" atoms = Arginine
  (createHydrated CB atoms)
  (createHydrated CG atoms)
  (createHydrated CD atoms)
  (createHydrated NE atoms)
  (createHydrated CZ atoms)
  (createHydrated NH1 atoms)
  (createHydrated NH2 atoms)
pdbLinesToRadical "HIS" atoms = Histidine
  (createHydrated CB atoms)
  (createHydrated CG atoms)
  (createHydrated ND1 atoms)
  (createHydrated CE1 atoms)
  (createHydrated NE2 atoms)
  (createHydrated CD2 atoms)
pdbLinesToRadical "LYS" atoms = Lysine
  (createHydrated CB atoms)
  (createHydrated CG atoms)
  (createHydrated CD atoms)
  (createHydrated CE atoms)
  (createHydrated NZ atoms)
pdbLinesToRadical "ASP" atoms = AsparticAcid
  (createHydrated CB atoms)
  (createHydrated CG atoms)
  (createHydrated OD1 atoms)
  (createHydrated OD2 atoms)
pdbLinesToRadical "GLU" atoms = GlutamicAcid
  (createHydrated CB atoms)
  (createHydrated CG atoms)
  (createHydrated CD atoms)
  (createHydrated OE1 atoms)
  (createHydrated OE2 atoms)
pdbLinesToRadical amino _ = error $ "Unknown amino acid code => " ++ amino
