# LDA
R code for linear discriminant analysis

---
The set includes R code and data files used in "Quantitative discrimination of flightlessness in fossil Anatidae from skeletal proportions" by Junya Watanabe on The Auk: Ornithological Advances, 134(3): 672-695.

Contents:
README.md: This file.
descriptions: Detailed descriptions of functions included.
examples.txt: Examples showing how to use the functions.
FAN.csv: Measurement data for fossil anatids (species means).
FANI.csv: Measurement data for fossil anatids (individuals).
FANN.csv: Sample sizes for fossil anatids.
functions.R: Includes codes of functions (readable with a text editor).
GR.csv: Group memberships of modern anatid species.
MAN.csv: Measurement data for modern anatids (species means).
MANI.csv: Measurement data for modern anatids (individuals).
MANN.csv: Sample size for modern anatids.

- Most functions were written on R 3.2.0.
- Codes may be slow, although I believe this would not cause practical problems.
- Functions Coef.B and Coef.pc.B depends on the package rrcov.
- A large part of the codes derives from those on Shigenobu Aoki's webpage (in Japanese): http://aoki2.si.gunma-u.ac.jp/R/ (accessed May 2015)
