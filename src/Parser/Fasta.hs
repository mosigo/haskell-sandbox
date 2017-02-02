module Parser.Fasta
    ( FastaSequence (..)
    , fastaParser
    ) where

import Text.Parsec (many, char, noneOf, newline, letter)
import Text.Parsec.Text (Parser)
import Data.Text as T (Text, strip, pack, concat)

data FastaSequence = FastaSequence { name :: Text
                                   , seq  :: Text
                                   } deriving (Show)

fastaParser :: Parser [FastaSequence]
fastaParser = many fastaSequence

fastaSequence :: Parser FastaSequence
fastaSequence = FastaSequence <$> nameP <*> sequenceP

nameP :: Parser Text
nameP = strip <$> (T.pack <$> nameP')

nameP' :: Parser String
nameP' =  char '>' *> simpleSpaces *> many (noneOf ['\n']) <* newline

sequenceP :: Parser Text
sequenceP = T.concat <$> many lineP

lineP :: Parser Text
lineP = T.pack <$> many letter <* simpleSpaces <* newline

simpleSpace :: Parser Char
simpleSpace = char ' '

simpleSpaces :: Parser String
simpleSpaces = many simpleSpace
