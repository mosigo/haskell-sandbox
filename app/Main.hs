{-# LANGUAGE OverloadedStrings #-}

module Main where

import Neo4j
import System.IO hiding (hGetContents)
import System.IO.Strict
import Data.Text as T hiding (head)
import Database.Bolt
import System.Environment

main :: IO ()
main = do
  args <- getArgs
  handle <- openFile (head args) ReadMode
  contents <- T.pack <$> hGetContents handle
  hClose handle
  loadFromCypher conf contents where
    conf = def { host="192.168.40.166", user="neo4j", password="lvbm123" }
