# 'Predict mammal body mass from lower m1 area
#'
#' This is based on the work of S. Legrendre 1986 and is for all mammals.
#'
#' @param area numeric; length x width from lower m1
#' @return mass in grams
#' @export
general.m1area <- function(area) {
  la <- log(area)
  lnm <- 1.827 * la + 1.810
  mass <- exp(lnm)
  mass
}

general.mandible <- function(length) {
  # mandible length
  # slater from foster 2009
  la <- log(length)
  lnm <- la * 2.9677 - 5.6712
  mass <- exp(lnm)
  mass
}

general.skull <- function(length) {
  # skull length
  # slater from Luo et al. 2001
  la <- log(length, base = 10)
  lnm <- la * 3.68 - 3.83
  mass <- 10^lnm
  mass
}

general.mass <- function(species, measures) {
  gen.mea <- measures[measures$species %in% species, ]

  gen.m1 <- gen.mea[gen.mea$part %in% c('m1', 'Lower m1', 'Lower m1 ') & 
                gen.mea$measure == 'area', ]
  gen.m1$mass <- general.m1area(gen.m1$value)

  gen.mandib <- gen.mea[gen.mea$part %in% c('mandible', 'length of mandible', 
                                            'Mandibular') & 
                        gen.mea$measure == 'length', ]
  gen.mandib$mass <- general.mandible(gen.mandib$value)

  gen.skull <- gen.mea[gen.mea$part %in% c('skull', 'Entire skull length measured 
                                           from I1 to occipital condyles') & 
                       gen.mea$measure == 'length', ]
  gen.skull$mass <- general.skull(gen.skull$value)

  out <- rbind(gen.m1, gen.mandib, gen.skull)
  out <- out[!duplicated(out$species), ]

  out
}

# ungulates
ungulates.m1area <- function(area) {
  la <- log(area)
  lnm <- 3.757 + 1.516 * la
  mass <- exp(lnm)
  mass
}

ungulates.m2length <- function(len) {
  la <- log(len)
  lnm <- 2.366 + 3.076 * la
  mass <- exp(lnm)
  mass
}

ungulates.M2length <- function(len) {
  la <- log(len)
  lnm <- 2.475 + 3.004 * la
  mass <- exp(lnm)
  mass
}

ungulates.M2area <- function(area) {
  la <- log(area)
  lnm <- 2.792 + 1.518 * la
  mass <- exp(lnm)
  mass
}

ungulates.ltrl <- function(length) {
  # ungulates
  # susumu
  la <- log(length)
  lnm <- -1.374 + (3.113 * la)
  mass <- exp(lnm)
  mass
}

ungulate.mass <- function(species, measures) {
  mea <- measures[measures$species %in% species, ]

  m1 <- mea[mea$part %in% c('m1', 'Lower m1', 'Lower m1 ') & 
                mea$measure == 'area', ]
  m1$mass <- ungulates.m1area(m1$value)

  m2 <- mea[mea$part %in% c('m2', 'Lower m2', 'lower M2') & 
                mea$measure == 'length', ]
  m2$mass <- ungulates.m2length(m2$value)

  um2a <- mea[mea$part %in% c('M2') & 
              mea$measure == 'area', ]
  um2a$mass <- ungulates.M2area(um2a$value)

  um2l <- mea[mea$part %in% c('M2') & 
              mea$measure == 'length', ]
  um2l$mass <- ungulates.M2length(um2l$value)

  ltrl <- mea[mea$part %in% c('LTRL') & 
              mea$measure == 'length', ]
  ltrl$mass <- ungulates.ltrl(ltrl$value)

  out <- rbind(m1, m2, um2a, um2l, ltrl)
  out <- out[!duplicated(out$species), ]

  out
}

# carnivores
carnivore.m1length <- function(len) {
  la <- log(len)
  lnm <- 1.681 + 2.97 * la
  mass <- exp(lnm)
  mass
}

carnivore.mass <- function(species, measures) {
  mea <- measures[measures$species %in% species, ]

  m1 <- mea[mea$part %in% c('m1', 'Lower m1', 'Lower m1 ') & 
                mea$measure == 'length', ]
  m1$mass <- carnivore.m1length(m1$value) 

  out <- m1
  out <- out[!duplicated(out$species), ]

  out
}

lagomorph.larl <- function(larl) {
  la <- log(larl)
  lnm <- -2.671 + 3.671 * la
  mass <- exp(lnm)
  mass
}

lagomorph.m1length <- function(len) {
  la <- log(len)
  lnm <- 3.002 + 4.468 * la
  mass <- exp(lnm)
  mass
}

lagomorph.mass <- function(species, measures) {
  mea <- measures[measures$species %in% species, ]

  m1 <- mea[mea$part %in% c('m1', 'Lower m1', 'Lower m1 ') & 
                mea$measure == 'length', ]
  m1$mass <- lagomorph.m1length(m1$value) 

  out <- m1
  out <- out[!duplicated(out$species), ]

  out
}

# insectivores
insectivore.m1area <- function(area) {
  la <- log(area)
  lnm <- 1.726 + 1.628 * la
  mass <- exp(lnm)
  mass
}

insectivore.M1area <- function(area) {
  # slater from bloch et al 1998
  la <- log(area)
  lnm <- 0.886 + (la * 1.714)
  mass <- exp(lnm)
  mass
}

insectivore.mass <- function(species, measures) {
  mea <- measures[measures$species %in% species, ]

  m1 <- mea[mea$part %in% c('m1', 'Lower m1', 'Lower m1 ') & 
                mea$measure == 'area', ]
  m1$mass <- insectivore.m1area(m1$value) 
  um1 <- mea[mea$part %in% c('M1', 'Upper M1', 'Upper M1 ') &
                mea$measure == 'area', ]
  um1$mass <- insectivore.M1area(um1$value) 

  out <- rbind(m1, um1)
  out <- out[!duplicated(out$species), ]

  out
}

# rodents
rodentia.m1area <- function(area) {
  la <- log(area)
  lnm <- 2.172 + 1.767 * la
  mass <- exp(lnm)
  mass
}

rodentia.mass <- function(species, measures) {
  mea <- measures[measures$species %in% species, ]

  m1 <- mea[mea$part %in% c('m1', 'Lower m1', 'Lower m1 ') & 
                mea$measure == 'area', ]
  m1$mass <- rodentia.m1area(m1$value) 

  out <- m1
  out <- out[!duplicated(out$species), ]

  out
}

marsupial.M1length <- function(len) {
  la <- log(len)
  lnm <- 1.83 + (la * 3.284)
  mass <- exp(lnm)
  mass
}

marsupial.M1area <- function(area) {
  la <- log(area)
  lnm <- 1.571 + (la * 1.733)
  mass <- exp(lnm)
  mass
}

marsupial.mass <- function(species, measures) {
  mea <- measures[measures$species %in% species, ]

  m1 <- mea[mea$part %in% c('M1', 'Upper M1', 'Upper M1 ') &
            mea$measure == 'area', ]
  m1$mass <- marsupial.M1area(m1$value) 

  um1 <- mea[mea$part %in% c('M1', 'Upper M1', 'Upper M1 ') &
                mea$measure == 'length', ]
  um1$mass <- marsupial.M1length(um1$value) 

  out <- rbind(m1, um1)
  out <- out[!duplicated(out$species), ]

  out
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
