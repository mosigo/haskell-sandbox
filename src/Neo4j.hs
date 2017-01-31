{-# LANGUAGE OverloadedStrings #-}

module Neo4j
    ( someFunc
    , loadFromCypher
    ) where

import Database.Bolt
import Data.Text as T

someFunc :: IO ()
someFunc = putStrLn "someFunc"

loadFromCypher :: BoltCfg -> Text -> IO ()
loadFromCypher cfg cypherQuery = do
  pipe <- connect cfg
  run pipe $ query cypherQuery
  close pipe
