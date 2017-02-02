module Main where

import Parser.Fasta as F (FastaSequence (..), fastaParser)
import System.Environment (getArgs)
import Text.Parsec (parse)
import qualified Data.Text.IO as TIO (readFile)

main :: IO ()
main = do
  args <- getArgs
  contents <- TIO.readFile $ head args
  result $ parse fastaParser "fasta parser" contents

result :: (Show a) => Either a [FastaSequence] -> IO ()
result (Left err)        = putStr "Parse error: " >> print err
result (Right sequences) = mconcat $ printSequence <$> sequences

printSequence :: F.FastaSequence -> IO ()
printSequence fastaSequence = do
  putStr "name: "
  print $ F.name fastaSequence
  putStr "seq:  "
  print $ F.seq fastaSequence
  putStr "======\n"
