module Parser.PDB
    ( PDBLine (..)
    , pdbLinesParser
    ) where

import Data.Text as T (Text, pack, strip)

import Text.Parsec ( many, count
                   , anyChar, spaces
                   )
import Text.Parsec.Text (Parser)

data PDBLine = PDBLine { lineHeader :: Text
                       , lineText   :: Text
                       } deriving (Show)

pdbLinesParser :: Parser [PDBLine]
pdbLinesParser = many pdbLineP

pdbLineP :: Parser PDBLine
pdbLineP = PDBLine . strip . pack <$> count 6 anyChar <*> (pack <$> count 74 anyChar) <* spaces

-- simpleSpace :: Parser Char
-- simpleSpace = char ' '
--
-- simpleSpaces :: Parser String
-- simpleSpaces = many simpleSpace
