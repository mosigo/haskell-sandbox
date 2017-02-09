module Main where

import System.Environment (getArgs)
import Text.Parsec (parse)
import qualified Data.Text.IO as TIO (readFile)
import AminoAcidPDB (pdbAminoParser)
import AminoAcid (HydratedAminoAcid)

main :: IO ()
main = do
  args <- getArgs
  contents <- TIO.readFile $ head args
  result $ parse pdbAminoParser "pdb parser" contents

result :: (Show a) => Either a [HydratedAminoAcid] -> IO ()
result (Left err)       = putStr "Parse error: " >> print err
result (Right pdbLines) = mconcat $ printLine <$> pdbLines

printLine :: HydratedAminoAcid -> IO ()
printLine = print
