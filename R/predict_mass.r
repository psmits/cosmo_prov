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
  mass <- exp(lnm)
  mass
}

massM1len <- function(length) {
  la <- log(length)
  lnm <- 2.475 + (la * 3.004)
  mass <- exp(lnm)
  mass
}

massM1 <- function(area) {
  la <- log(area)
  lnm <- 2.792 + (la * 1.518)
  mass <- exp(lnm)
  mass
}

massm2len <- function(length) {
  la <- log(length)
  lnm <- 2.355 + (la * 3.076)
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
