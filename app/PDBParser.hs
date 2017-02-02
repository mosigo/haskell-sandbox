module Main where

import Parser.PDB as P (PDBLine (..), pdbLinesParser)
import System.Environment (getArgs)
import Text.Parsec (parse)
import qualified Data.Text.IO as TIO (readFile)

main :: IO ()
main = do
  args <- getArgs
  contents <- TIO.readFile $ head args
  result $ parse pdbLinesParser "pdb parser" contents

result :: (Show a) => Either a [PDBLine] -> IO ()
result (Left err)        = putStr "Parse error: " >> print err
result (Right pdbLines) = mconcat $ printLine <$> pdbLines

printLine :: P.PDBLine -> IO ()
printLine pdbLine = do
  putStr "name: "
  print $ P.lineHeader pdbLine
  putStr "seq:  "
  print $ P.lineText pdbLine
  putStr "======\n"
