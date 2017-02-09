{-# LANGUAGE OverloadedStrings #-}

module Main where

import Neo4j (loadFromCypher)
import Database.Bolt
import System.Environment (getArgs)
import System.IO (putStrLn)
import qualified Data.Text.IO as TIO (readFile)

main :: IO ()
main = do
  args <- getArgs
  contents <- TIO.readFile $ head args
  putStrLn "Read"
  loadFromCypher conf contents
  putStrLn "Finished"
  where
    conf = def { host="192.168.40.102", user="neo4j", password="lvbm123" }
