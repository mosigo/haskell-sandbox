module Parser.PDB
    ( PDBLine (..)
    , pdbLinesParser
    ) where

import           Data.Text        as T (Text, pack, strip)

import           Control.Monad    (void)
import           Text.Parsec      (anyChar, char, count, eof, many, noneOf,
                                   string, try, (<|>))
import           Text.Parsec.Text (Parser)

data PDBLine = AtomLine  { atomSerial     :: Int     -- Atom  serial number.
                         , atomName       :: Text    -- Atom name.
                         , atomAltLoc     :: Char    -- Alternate location indicator.
                         , atomResName    :: Text    -- Residue name.
                         , atomChainID    :: Char    -- Chain identifier.
                         , atomResSeq     :: Int     -- Residue sequence number.
                         , atomICode      :: Char    -- Code for insertion of residues.
                         , atomX          :: Float   -- Orthogonal coordinates for X in Angstroms.
                         , atomY          :: Float   -- Orthogonal coordinates for Y in Angstroms.
                         , atomZ          :: Float   -- Orthogonal coordinates for Z in Angstroms.
                         , atomOccupancy  :: Float   -- Occupancy.
                         , atomTempFactor :: Float   -- Temperature  factor.
                         , atomElement    :: Text    -- Element symbol, right-justified.
                         , atomCharge     :: Text    -- Charge  on the atom.
                         }
             | OtherLine { otherHeader :: Text
                         , otherText   :: Text
                         }
             deriving (Show)

pdbLinesParser :: Parser [PDBLine]
pdbLinesParser = many (atomLineP <|> otherLineP)

otherLineP :: Parser PDBLine
otherLineP = OtherLine . T.strip . T.pack <$> count 6 anyChar <*> (pack <$> otherLineSymbols)

atomLineP :: Parser PDBLine
atomLineP = AtomLine <$> (try (string "ATOM  ")                                   -- (1 -  6)
                           *> (read <$> count 5 anyChar))                         -- (7 - 11)  atomSerial
                     <*> (T.strip . T.pack <$> (anyChar *> count 4 anyChar))      -- (13 - 16) atomName
                     <*> anyChar                                                  -- (17)      atomAltLoc
                     <*> (T.pack <$> count 3 anyChar)                             -- (18 - 20) atomResName
                     <*> (anyChar *> anyChar)                                     -- (22)      atomChainID
                     <*> (read <$> count 4 anyChar)                               -- (23 - 26) atomResSeq
                     <*> anyChar                                                  -- (27)      atomICode
                     <*> (read <$> (count 3 anyChar *> count 8 anyChar))          -- (31 - 38) atomX
                     <*> (read <$> count 8 anyChar)                               -- (39 - 46) atomY
                     <*> (read <$> count 8 anyChar)                               -- (47 - 54) atomZ
                     <*> (read <$> count 6 anyChar)                               -- (55 - 60) atomOccupancy
                     <*> (read <$> count 6 anyChar)                               -- (61 - 66) atomTempFactor
                     <*> (strip . pack <$> (count 10 anyChar *> count 2 anyChar)) -- (77 - 78) atomElement
                     <*> (strip . pack <$> count 2 anyChar                        -- (79 - 80) atomCharge
                            <* otherLineSymbols)

eol :: Parser ()
eol = void (char '\n') <|> eof

otherLineSymbols :: Parser String
otherLineSymbols = many (noneOf ['\n']) <* eol
