#' Predict mammal body mass from lower m1 area
#'
#' This is based on the work of S. Legrendre 1986 and is for all mammals.
#'
#' @param area numeric; length x width from lower m1
#' @return mass in grams
#' @export
predmass <- function(area) {
  la <- log(area)
  lnm <- 1.827 * la + 1.810
  # lnm <- 3.757 + la * 1.516
  mass <- exp(lnm)
  mass
}

massM1len <- function(length) {
  # marsupials...
  la <- log(length)
  lnm <- 1.796 + (la * 3.3)
  mass <- exp(lnm)
  mass
}

massm1len <- function(length) {
  # susumu carnivores
  la <- log(length)
  lnm <- 1.681 + (la * 2.97)
  mass <- exp(lnm)
  mass
}

massm1len.lago <- function(length) {
  # susumu lagomorphs
  la <- log(length)
  lnm <- 3.002 + (la * 4.468)
  mass <- exp(lnm)
  mass
}

massM1 <- function(area) {
  # slater from bloch et al 1998
  la <- log(area)
  lnm <- 0.886 + (la * 1.714)
  mass <- exp(lnm)
  mass
}

massm2len <- function(length) {
  # susumu
  la <- log(length)
  lnm <- 2.355 + (la * 3.076)
  mass <- exp(lnm)
  mass
}

massM2area <- function(area) {
  # susumu ungulates
  la <- log(area)
  lnm <- 2.792 + (la * 1.518)
  mass <- exp(lnm)
  mass
}

massML <- function(length) {
  # mandible length
  # slater from foster 2009
  la <- log(length)
  lnm <- la * 2.9677 - 5.6712
  mass <- exp(lnm)
  mass
}

massSL <- function(length) {
  # skull length
  # slater from Luo et al. 2001
  la <- log(length, base = 10)
  lnm <- la * 3.68 - 3.83
  mass <- 10^lnm
  mass
}

mass.ltrl <- function(length) {
  # ungulates
  # susumu
  la <- log(length)
  lnm <- -1.374 + (3.113 * la)
  mass <- exp(lnm)
  mass
}


#' Predict herbivore body mass from lower m1 area
#'
#' This is based on the work of S. Legrendre 1986 and is for herbivores.
#' This is defined as rodents, primates, artiodactyls, perissodactyls.
#'
#' @param area numeric; length x width from lower m1
#' @return mass in grams
#' @export
predherb <- function(area) {
  la <- log(area)
  lnm <- 1.693 * la + 2.694
  mass <- exp(lnm)
  mass
}


#' Predict faunivore body mass from lower m1 area
#'
#' This is based on the work of S. Legrendre 1986 and is for faunivores.
#' This is defined as insectivores, carnivores, marsupials, lipotyphla, bats
#'
#' @param area numeric; length x width from lower m1
#' @return mass in grams
#' @export
predfaun <- function(area) {
  la <- log(area)
  lnm <- 1.739 * la + 1.493
  mass <- exp(lnm)
  mass
}


#' Predict small mammal body mass from lower m1 area
#'
#' This is based on the work of S. Legrendre 1986 and is for small mammals (<500 g) 
#'
#' @param area numeric; length x width from lower m1
#' @return mass in grams
#' @export
predsmal <- function(area) {
  la <- log(area)
  lnm <- 1.621 * la + 1.786
  mass <- exp(lnm)
  mass
}


#' Predict large mammal body mass from lower m1 area
#'
#' This is based on the work of S. Legrendre 1986 and is for small mammals (>500 g) 
#'
#' @param area numeric; length x width from lower m1
#' @return mass in grams
#' @export
predlarg <- function(area) {
  la <- log(area)
  lnm <- 1.538 * la + 3.115
  mass <- exp(lnm)
  mass
}


#' Predict insectivore body mass from lower m1 area
#'
#' This is based on the work of S. Legrendre 1986 and is for insectivores.
#'
#' @param area numeric; length x width from lower m1
#' @return mass in grams
#' @export
predinst <- function(area) {
  la <- log(area)
  lnm <- 1.654 * la + 1.746
  mass <- exp(lnm)
  mass
}


#' Predict bat body mass from lower m1 area
#'
#' This is based on the work of S. Legrendre 1986 and is for bats 
#'
#' @param area numeric; length x width from lower m1
#' @return mass in grams
#' @export
predbats <- function(area) {
  la <- log(area)
  lnm <- 1.289 * la + 1.829
  mass <- exp(lnm)
  mass
}


#' Predict carnivore body mass from lower m1 area
#'
#' This is based on the work of S. Legrendre 1986 and is for carnivores
#'
#' @param area numeric; length x width from lower m1
#' @return mass in grams
#' @export
predcarn <- function(area) {
  la <- log(area)
  lnm <- 1.922 * la + 0.709
  mass <- exp(lnm)
  mass
}


#' Predict primate body mass from lower m1 area
#'
#' This is based on the work of S. Legrendre 1986 and is for primates
#'
#' @param area numeric; length x width from lower m1
#' @return mass in grams
#' @export
predmonk <- function(area) {
  la <- log(area)
  lnm <- 1.490 * la + 3.577
  mass <- exp(lnm)
  mass
}


#' Predict rodent body mass from lower m1 area
#'
#' This is based on the work of S. Legrendre 1986 and is for rodents
#'
#' @param area numeric; length x width from lower m1
#' @return mass in grams
#' @export
predrats <- function(area) {
  la <- log(area)
  lnm <- 1.767 * la + 2.172
  mass <- exp(lnm)
  mass
}


#' Predict ungulate body mass from lower m1 area
#'
#' This is based on the work of S. Legrendre 1986 and is for ungulates
#'
#' @param area numeric; length x width from lower m1
#' @return mass in grams
#' @export
predungl <- function(area) {
  la <- log(area)
  lnm <- 1.564 * la + 3.267
  mass <- exp(lnm)
  mass
}
