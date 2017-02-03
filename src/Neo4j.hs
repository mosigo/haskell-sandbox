{-# LANGUAGE OverloadedStrings #-}

module Neo4j
    ( someFunc
    , loadFromCypher
    ) where

import Database.Bolt
import Data.Text as T hiding (zip)

someFunc :: IO ()
someFunc = putStrLn "someFunc"

loadFromCypher :: BoltCfg -> Text -> IO ()
loadFromCypher cfg cypherQuery = do
  pipe <- connect cfg
  mconcat $ runQuery pipe <$> zip (T.splitOn "@@@@@" cypherQuery) [1..]
  close pipe

runQuery :: Pipe -> (Text, Int) -> IO ()
runQuery pipe (cypherQuery, queryNum)
  | T.length (T.strip cypherQuery) == 0 = putStrLn "Empty string"
  | otherwise                           = do
    _ <- run pipe $ query cypherQuery
    putStrLn $ "query " ++ show queryNum
