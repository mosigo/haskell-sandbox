module Main where

import Test.Hspec
import Test.QuickCheck
import Data.Text (Text (..), pack)

import Text.Parsec (parse)
import AminoAcidPDB (pdbAminoParser)

import Debug.Trace (trace)

main :: IO ()
main = hspec $
  describe "AminoAcidPDB" $ do
    testPDBAminoParser

testPDBAminoParser :: SpecWith ()
testPDBAminoParser =
  it "parse simple pdb" $
    case parse pdbAminoParser "pdb parser" testPDB of
      (Left err)     -> trace (show err) (2 `shouldBe` 1)
      (Right aminos) -> length aminos `shouldBe` 3

testPDB :: Text
testPDB = pack
  "REMARK  99 H H\n\
  \REMARK  98 H EIQLVQSGAEVKKPGSSVKVSCKASGYTFTDYYINWMRQAPGQGLEWIGWIDPGSGNTKYNEKFKGRATLTVDTSTNTAYMELSSLRSEDTAFYFCAREKTTYYYAMDYWGQGTLVTVSS\n\
  \REMARK  97 H EIQLVQSGAEVKKPGSSVKVSCKASGYTFT DYYIN WMRQAPGQGLEWIG WIDPGSGNTKYNEKFKG RATLTVDTSTNTAYMELSSLRSEDTAFYFCAR EKTTYYYAMDY WGQGTLVTVSS\n\
  \REMARK  96 H 0,850574713\n\
  \REMARK  95 H 0,765306122\n\
  \REMARK  99 L L\n\
  \REMARK  98 L DIQMTQSPSTLSASVGDRVTITCRSSKSLLHSNGDTFLYWFQQKPGKAPKLLMYRMSNLASGVPSRFSGSGSGTEFTLTISSLQPDDFATYYCMQHLEYPFTFGQGTKVEVK\n\
  \REMARK  97 L DIQMTQSPSTLSASVGDRVTITC RSSKSLLHSNGDTFLY WFQQKPGKAPKLLMY RMSNLAS GVPSRFSGSGSGTEFTLTISSLQPDDFATYYC MQHLEYPFT FGQGTKVEVK\n\
  \REMARK  96 L 0,925\n\
  \REMARK  95 L 0,784946237\n\
  \ATOM      1  N   ASP L   1      26.302  31.499  -5.832  1.00 65.29           N1+\n\
  \ATOM      2  CA  ASP L   1      26.975  32.715  -5.266  1.00 63.66           C\n\
  \ATOM      3  C   ASP L   1      25.853  33.775  -5.254  1.00 57.26           C\n\
  \ATOM      4  O   ASP L   1      25.285  34.079  -6.306  1.00 59.38           O\n\
  \ATOM      5  CB  ASP L   1      28.210  33.112  -6.118  1.00 64.29           C\n\
  \ATOM      6  CG  ASP L   1      29.149  34.107  -5.423  1.00 68.69           C\n\
  \ATOM      7  OD1 ASP L   1      28.785  35.179  -4.947  1.00 67.27           O\n\
  \ATOM      8  OD2 ASP L   1      30.442  33.675  -5.409  1.00 64.31           O\n\
  \ATOM      9  H1  ASP L   1      27.001  30.806  -6.056  1.00  0.00           H\n\
  \ATOM     10  H2  ASP L   1      25.656  31.125  -5.149  1.00  0.00           H\n\
  \ATOM     11  H3  ASP L   1      25.803  31.767  -6.669  1.00  0.00           H\n\
  \ATOM     12  HA  ASP L   1      27.258  32.418  -4.253  1.00  0.00           H\n\
  \ATOM     13  HB3 ASP L   1      27.914  33.525  -7.084  1.00  0.00           H\n\
  \ATOM     14  HB2 ASP L   1      28.789  32.215  -6.345  1.00  0.00           H\n\
  \ATOM     15  HD2 ASP L   1      31.028  34.300  -4.975  1.00  0.00           H\n\
  \ATOM     16  N   ILE L   2      25.555  34.302  -4.056  1.00 47.12           N\n\
  \ATOM     17  CA  ILE L   2      24.510  35.301  -3.828  1.00 41.15           C\n\
  \ATOM     18  C   ILE L   2      25.172  36.684  -3.727  1.00 38.73           C\n\
  \ATOM     19  O   ILE L   2      26.146  36.833  -2.990  1.00 40.24           O\n\
  \ATOM     20  CB  ILE L   2      23.722  35.009  -2.512  1.00 41.68           C\n\
  \ATOM     21  CG1 ILE L   2      23.089  33.596  -2.580  1.00 42.01           C\n\
  \ATOM     22  CG2 ILE L   2      22.657  36.085  -2.178  1.00 37.08           C\n\
  \ATOM     23  CD1 ILE L   2      22.307  33.143  -1.337  1.00 34.46           C\n\
  \ATOM     24  H   ILE L   2      26.074  34.019  -3.237  1.00  0.00           H\n\
  \ATOM     25  HA  ILE L   2      23.803  35.301  -4.659  1.00  0.00           H\n\
  \ATOM     26  HB  ILE L   2      24.435  35.004  -1.689  1.00  0.00           H\n\
  \ATOM     27 HG13 ILE L   2      23.862  32.847  -2.756  1.00  0.00           H\n\
  \ATOM     28 HG12 ILE L   2      22.443  33.548  -3.452  1.00  0.00           H\n\
  \ATOM     29 HG21 ILE L   2      22.167  35.883  -1.228  1.00  0.00           H\n\
  \ATOM     30 HG22 ILE L   2      23.081  37.083  -2.080  1.00  0.00           H\n\
  \ATOM     31 HG23 ILE L   2      21.881  36.127  -2.943  1.00  0.00           H\n\
  \ATOM     32 HD11 ILE L   2      21.954  32.120  -1.461  1.00  0.00           H\n\
  \ATOM     33 HD12 ILE L   2      22.931  33.156  -0.447  1.00  0.00           H\n\
  \ATOM     34 HD13 ILE L   2      21.428  33.759  -1.156  1.00  0.00           H\n\
  \ATOM     35  N   GLN L   3      24.628  37.666  -4.464  1.00 40.59           N\n\
  \ATOM     36  CA  GLN L   3      25.099  39.050  -4.445  1.00 38.67           C\n\
  \ATOM     37  C   GLN L   3      24.312  39.844  -3.400  1.00 34.84           C\n\
  \ATOM     38  O   GLN L   3      23.084  39.865  -3.453  1.00 34.19           O\n\
  \ATOM     39  CB  GLN L   3      24.986  39.707  -5.841  1.00 43.54           C\n\
  \ATOM     40  CG  GLN L   3      26.056  39.252  -6.856  1.00 56.54           C\n\
  \ATOM     41  CD  GLN L   3      25.866  37.812  -7.337  1.00 67.08           C\n\
  \ATOM     42  OE1 GLN L   3      24.825  37.473  -7.894  1.00 69.48           O\n\
  \ATOM     43  NE2 GLN L   3      26.868  36.957  -7.123  1.00 68.42           N\n\
  \ATOM     44  H   GLN L   3      23.819  37.472  -5.040  1.00  0.00           H\n\
  \ATOM     45  HA  GLN L   3      26.154  39.071  -4.166  1.00  0.00           H\n\
  \ATOM     46  HB3 GLN L   3      25.093  40.787  -5.729  1.00  0.00           H\n\
  \ATOM     47  HB2 GLN L   3      23.990  39.560  -6.258  1.00  0.00           H\n\
  \ATOM     48  HG3 GLN L   3      27.052  39.382  -6.432  1.00  0.00           H\n\
  \ATOM     49  HG2 GLN L   3      26.013  39.897  -7.734  1.00  0.00           H\n\
  \ATOM     50 HE22 GLN L   3      26.772  35.995  -7.413  1.00  0.00           H\n\
  \ATOM     51 HE21 GLN L   3      27.711  37.254  -6.652  1.00  0.00           H\n\
  \ATOM     52  N   MET L   4      25.053  40.490  -2.492  1.00 32.30           N\n\
  \ATOM     53  CA  MET L   4      24.523  41.360  -1.449  1.00 27.06           C\n\
  \ATOM     54  C   MET L   4      24.878  42.802  -1.818  1.00 28.63           C\n\
  \ATOM     55  O   MET L   4      26.049  43.087  -2.070  1.00 30.43           O\n\
  \ATOM     56  CB  MET L   4      25.135  40.961  -0.090  1.00 29.08           C\n\
  \ATOM     57  CG  MET L   4      24.858  39.505   0.319  1.00 27.54           C\n\
  \ATOM     58  SD  MET L   4      23.117  39.068   0.526  1.00 25.60           S\n\
  \ATOM     59  CE  MET L   4      22.830  39.751   2.174  1.00 18.57           C\n\
  \ATOM     60  H   MET L   4      26.060  40.430  -2.536  1.00  0.00           H\n\
  \ATOM     61  HA  MET L   4      23.442  41.260  -1.378  1.00  0.00           H\n\
  \ATOM     62  HB3 MET L   4      24.777  41.635   0.689  1.00  0.00           H\n\
  \ATOM     63  HB2 MET L   4      26.214  41.098  -0.129  1.00  0.00           H\n\
  \ATOM     64  HG3 MET L   4      25.376  39.279   1.250  1.00  0.00           H\n\
  \ATOM     65  HG2 MET L   4      25.264  38.830  -0.430  1.00  0.00           H\n\
  \ATOM     66  HE1 MET L   4      21.762  39.857   2.361  1.00  0.00           H\n\
  \ATOM     67  HE2 MET L   4      23.257  39.105   2.939  1.00  0.00           H\n\
  \ATOM     68  HE3 MET L   4      23.297  40.727   2.262  1.00  0.00           H\n\
  \ATOM     69  N   THR L   5      23.862  43.674  -1.862  1.00 27.42           N\n\
  \ATOM     70  CA  THR L   5      24.001  45.089  -2.197  1.00 31.09           C\n\
  \ATOM     71  C   THR L   5      23.300  45.927  -1.118  1.00 30.29           C\n\
  \ATOM     72  O   THR L   5      22.192  45.589  -0.707  1.00 32.10           O\n\
  \ATOM     73  CB  THR L   5      23.389  45.419  -3.591  1.00 31.96           C\n\
  \ATOM     74  OG1 THR L   5      21.977  45.306  -3.640  1.00 40.13           O\n\
  \ATOM     75  CG2 THR L   5      23.979  44.575  -4.731  1.00 32.00           C\n\
  \ATOM     76  H   THR L   5      22.920  43.356  -1.666  1.00  0.00           H\n\
  \ATOM     77  HA  THR L   5      25.052  45.373  -2.208  1.00  0.00           H\n\
  \ATOM     78  HB  THR L   5      23.618  46.463  -3.812  1.00  0.00           H\n\
  \ATOM     79  HG1 THR L   5      21.602  45.887  -2.972  1.00  0.00           H\n\
  \ATOM     80 HG21 THR L   5      23.597  44.900  -5.699  1.00  0.00           H\n\
  \ATOM     81 HG22 THR L   5      25.066  44.663  -4.758  1.00  0.00           H\n\
  \ATOM     82 HG23 THR L   5      23.735  43.518  -4.620  1.00  0.00           H  "
