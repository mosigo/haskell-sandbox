module Parser.PDB
    ( PDBLine (..)
    , pdbLinesParser
    ) where

import Data.Text as T (Text, pack, strip)

import Text.Parsec ( many, count
                   , anyChar, spaces
                   , string
                   , (<|>)
                   )
import Text.Parsec.Text (Parser)

data PDBLine = PDBLine { lineHeader :: Text
                       , lineNum    :: Int
                       , lineText   :: Text
                       } deriving (Show)

pdbLinesParser :: Parser [PDBLine]
pdbLinesParser = many pdbLineP

pdbLineP :: Parser PDBLine
pdbLineP = PDBLine . strip . pack <$> count 6 anyChar <*> pdbLineNumP <*> (pack <$> count 70 anyChar) <* spaces

pdbLineNumP :: Parser Int
pdbLineNumP = parseNum <$> (string "    " <|> count 5 anyChar)

parseNum :: String -> Int
parseNum "    " = 1
parseNum s      = read s
