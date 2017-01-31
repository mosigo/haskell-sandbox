{-# LANGUAGE OverloadedStrings #-}

module Main where

import Neo4j
import System.IO hiding (hGetContents)
import System.IO.Strict
import Data.Text as T
import Database.Bolt

main :: IO ()
main = do
  handle <- openFile "/tmp/neo4j.init.txt" ReadMode
  contents <- T.pack <$> hGetContents handle
  hClose handle
  loadFromCypher conf contents where
    conf = def { host="192.168.40.166", user="neo4j", password="lvbm123" }
